// Copyright (c) 2015-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "netcdf.h"
#include "polyglot/cf_file.h"

#if POLYMEC_HAVE_DOUBLE_PRECISION
#define NC_REAL NC_DOUBLE
#else
#define NC_REAL NC_FLOAT
#endif

struct cf_file_t 
{
  int file_id;
  int cf_major_version, cf_minor_version, cf_patch_version;
  bool writing;

  // Important identifiers.
  int time_id, lat_id, lon_id, lev_id;
  char lev_name[POLYGLOT_CF_MAX_NAME+1];
};

// Helpers.
static void get_first_attribute(int file_id, 
                                int var_id, 
                                const char* attr,
                                char* value)
{
  int err = nc_get_att_text(file_id, var_id, attr, value);
  if (err == NC_ENOTATT)
    value[0] = '\0';
  else if (err != NC_NOERR)
  {
    polymec_error("cf_file: Error retrieving attribute %s: %s", 
                  attr, nc_strerror(err));
  }
}

static void get_first_global_attribute(int file_id,
                                       const char* attr, 
                                       char* value)
{
  get_first_attribute(file_id, NC_GLOBAL, attr, value);
}

static void put_attribute(int file_id, 
                          int var_id, 
                          const char* attr,
                          const char* value)
{
  int err = nc_put_att_text(file_id, var_id, attr, 1, value);
  if (err != NC_NOERR)
  {
    polymec_error("cf_file: Error setting attribute %s: %s", 
                  attr, nc_strerror(err));
  }
}

// Returns the ID of the variable with the given name, -1 if not found.
static int var_identifier(int file_id, const char* var_name)
{
  int err, id;
  err = nc_inq_varid(file_id, var_name, &id);
  if (err == NC_ENOTVAR)
    return -1;
  else if (err != NC_NOERR)
  {
    polymec_error("cf_file: Error retrieving var %s: %s",
                  var_name, nc_strerror(err));
  }
  else
    return id;
}

void find_vertical_dimension(int file_id, int* lev_id, char* lev_name)
{
  // This name should identify a dimension AND a variable, and the variable should 
  // have a "units" attribute, and a "positive" attribute (OR have a valid 
  // set of pressure units).
  int nvarsp;
  int err = nc_inq(file_id, NULL, &nvarsp, NULL, NULL);
  if (err != NC_NOERR)
    polymec_error("cf_file_open: Error retrieving number of vars: ", nc_strerror(err));
  for (int var_id = 0; var_id < nvarsp; ++var_id)
  {
    // Find the name of this variable.
    char var_name[POLYGLOT_CF_MAX_NAME+1];
    err = nc_inq_varname(file_id, var_id, var_name);
    if (err != NC_NOERR)
      polymec_error("cf_file_open: Error retrieving name of var %d: ", var_id, nc_strerror(err));

    // The vertical coordinate variable should have a single dimension.
    int ndim;
    err = nc_inq_varndims(file_id, var_id, &ndim);
    if (err != NC_NOERR)
      polymec_error("cf_file_open: Error retrieving number of dims for var %d: ", var_id, nc_strerror(err));
    if (ndim != 1) continue;

    // Find its dimension, and verify the dimensions's name is the same.
    int dim_id;
    char dim_name[POLYGLOT_CF_MAX_NAME+1];
    err = nc_inq_vardimid(file_id, var_id, &dim_id);
    if (err != NC_NOERR)
      polymec_error("cf_file_open: Error retrieving dim ID for var %d: ", var_id, nc_strerror(err));
    err = nc_inq_dimname(file_id, dim_id, dim_name);
    if (err != NC_NOERR)
      polymec_error("cf_file_open: Error retrieving dim name for var %d: ", var_id, nc_strerror(err));
    if (strcmp(dim_name, var_name) != 0) continue;

    // If there's an axis attribute set equal to Z, this is it.
    char axis[POLYGLOT_CF_MAX_NAME+1];
    get_first_attribute(file_id, var_id, "axis", axis);
    if (strcmp(axis, "Z") == 0)
    {
      *lev_id = var_id;
      strcpy(lev_name, var_name);
      return;
    }

    // If this has a standard name indicating a vertical coordinate, 
    // we're finished.
    char standard_name[POLYGLOT_CF_MAX_NAME+1];
    get_first_attribute(file_id, var_id, "standard_name", standard_name);
    if ((strcmp(standard_name, "altitude") == 0) ||
        (strcmp(standard_name, "height") == 0))
    {
      *lev_id = var_id;
      strcpy(lev_name, var_name);
      return;
    }

    // Look for a units attribute.
    char units[POLYGLOT_CF_MAX_NAME+1];
    get_first_attribute(file_id, var_id, "units", units);

    // Now look for a positive attribute with a valid value.
    char positive[POLYGLOT_CF_MAX_NAME+1];
    get_first_attribute(file_id, var_id, "positive", positive);
    if ((strlen(positive) > 0) &&
        (string_casecmp(positive, "up") != 0) && 
        (string_casecmp(positive, "down") != 0)) continue;

    // If we didn't find a positive attribute, see whether the units 
    // indicate a pressure.
    if ((strlen(positive) == 0) && 
        ((strcmp(units, "bar") == 0) || 
         (strcmp(units, "millibar") == 0) ||
         (strcmp(units, "decibar") == 0) || 
         (strcmp(units, "atmosphere") == 0) ||
         (strcmp(units, "atm") == 0) ||
         (strcmp(units, "pascal") == 0) ||
         (strcmp(units, "pa") == 0) ||
         (strcmp(units, "hPa") == 0)))
    {
      *lev_id = var_id;
      strcpy(lev_name, var_name);
      return;
    }
  }
}

// Implementation.

cf_file_t* cf_file_new(const char* filename)
{
  int file_id;
  int err = nc_open(filename, NC_WRITE, &file_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_new: Couldn't open file %s: %s", filename, nc_strerror(err));

  // Create our representation.
  cf_file_t* cf = polymec_malloc(sizeof(cf_file_t));
  cf->file_id = file_id;
  cf->cf_major_version = 1;
  cf->cf_minor_version = 6;
  cf->cf_patch_version = 0;
  cf->writing = true;

  cf->time_id = cf->lat_id = cf->lon_id = cf->lev_id = -1;
  strcpy(cf->lev_name, "lev");

  // Write in our conventions.
  char conventions[NC_MAX_NAME+1];
  snprintf(conventions, NC_MAX_NAME, "CF-%d.%d.%d", cf->cf_major_version, 
           cf->cf_minor_version, cf->cf_patch_version);
  err = nc_put_att_text(cf->file_id, NC_GLOBAL, "Conventions", 1, conventions);
  if (err != NC_NOERR)
    polymec_error("cf_file_new: Couldn't write Conventions attribute: %s", nc_strerror(err));

  return cf;
}

cf_file_t* cf_file_open(const char* filename)
{
  int file_id;
  int err = nc_open(filename, NC_NOWRITE, &file_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_open: Couldn't open file %s: %s", filename, nc_strerror(err));

  char conventions[NC_MAX_NAME+1];
  get_first_global_attribute(file_id, "Conventions", conventions);
  if (((conventions[0] != 'c') && (conventions[0] != 'C')) || 
      ((conventions[1] != 'f') && (conventions[1] != 'F')) || 
      (conventions[2] != '-') || (strlen(conventions) < 4))
  {
    nc_close(file_id);
    polymec_error("cf_file_open: File %s is not a CF-compliant NetCDF file.", filename);
  }

  // Create our representation.
  cf_file_t* cf = polymec_malloc(sizeof(cf_file_t));
  cf->file_id = file_id;
  cf->cf_major_version = cf->cf_minor_version = cf->cf_patch_version = 0;
  cf->writing = false;

  // Parse the CF conventions version numbers from the string.
  int num;
  char** versions = string_split(&conventions[3], ".", &num);
  ASSERT(num > 0);
  cf->cf_major_version = atoi(versions[0]);
  if (num > 1)
    cf->cf_minor_version = atoi(versions[1]);
  if (num > 2)
    cf->cf_patch_version = atoi(versions[2]);
  for (int i = 0; i < num; ++i)
    string_free(versions[i]);
  polymec_free(versions);

  // Snoop around and see what's here.
  int dim_id;
  err = nc_inq_dimid(cf->file_id, "lat", &dim_id);
  if (err == NC_NOERR)
    cf->lat_id = dim_id;
  else if (err != NC_EBADDIM)
    polymec_error("cf_file_open: Error retrieving lat dim ID: ", nc_strerror(err));

  err = nc_inq_dimid(cf->file_id, "lon", &dim_id);
  if (err == NC_NOERR)
    cf->lon_id = dim_id;
  else if (err != NC_EBADDIM)
    polymec_error("cf_file_open: Error retrieving lon dim ID: ", nc_strerror(err));

  err = nc_inq_dimid(cf->file_id, "time", &dim_id);
  if (err == NC_NOERR)
    cf->time_id = dim_id;
  else if (err != NC_EBADDIM)
    polymec_error("cf_file_open: Error retrieving time dim ID: ", nc_strerror(err));

  // We have to figure out the vertical dimension / coordinate name and ID. 
  find_vertical_dimension(cf->file_id, &cf->lev_id, cf->lev_name);

  return cf;
}

void cf_file_close(cf_file_t* file)
{
  nc_close(file->file_id);
}

void cf_file_get_version(cf_file_t* file, 
                         int* major_version,
                         int* minor_version,
                         int* patch_version)
{
  *major_version = file->cf_major_version;
  *minor_version = file->cf_minor_version;
  *patch_version = file->cf_patch_version;
}

void cf_file_get_provenance(cf_file_t* file, 
                            char* title,
                            char* institution,
                            char* source,
                            char* history,
                            char* references,
                            char* comment)
{
  cf_file_get_global_attribute(file, "title", title);
  cf_file_get_global_attribute(file, "institution", institution);
  cf_file_get_global_attribute(file, "source", source);
  cf_file_get_global_attribute(file, "history", history);
  cf_file_get_global_attribute(file, "references", references);
  cf_file_get_global_attribute(file, "comment", comment);
}

void cf_file_set_provenance(cf_file_t* file, 
                            const char* title,
                            const char* institution,
                            const char* source,
                            const char* history,
                            const char* references,
                            const char* comment)
{
  ASSERT(file->writing);
  cf_file_set_global_attribute(file, "title", title);
  cf_file_set_global_attribute(file, "institution", institution);
  cf_file_set_global_attribute(file, "source", source);
  cf_file_set_global_attribute(file, "history", history);
  cf_file_set_global_attribute(file, "references", references);
  cf_file_set_global_attribute(file, "comment", comment);
}

void cf_file_set_global_attribute(cf_file_t* file, 
                                  const char* global_attribute_name,
                                  const char* value)
{
  ASSERT(file->writing);
  put_attribute(file->file_id, NC_GLOBAL, global_attribute_name, value);
}

void cf_file_get_global_attribute(cf_file_t* file, 
                                  const char* global_attribute_name,
                                  char* value)
{
  get_first_global_attribute(file->file_id, global_attribute_name, value);
}

void cf_file_define_dimension(cf_file_t* file,
                              const char* dimension_name,
                              int value)
{
  ASSERT(file->writing);
  if (value == -1)
    value = NC_UNLIMITED;
  int id;
  int err = nc_def_dim(file->file_id, dimension_name, value, &id);
  if (err != NC_NOERR)
  {
    polymec_error("cf_file_define_dimension: Could not define dimension %s: %s",
                  dimension_name, nc_strerror(err));
  }
}

int cf_file_dimension(cf_file_t* file, const char* dimension_name)
{
  int id;
  int err = nc_inq_dimid(file->file_id, dimension_name, &id);
  if (err != NC_NOERR)
  {
    polymec_error("cf_file_dimension: Could not retrieve ID for dimension %s: %s",
                  dimension_name, nc_strerror(err));
  }
  size_t dim;
  err = nc_inq_dimlen(file->file_id, id, &dim);
  if (err != NC_NOERR)
  {
    polymec_error("cf_file_dimension: Could not retrieve dimension %s: %s",
                  dimension_name, nc_strerror(err));
  }
  return (int)dim;
}

void cf_file_define_latlon_grid(cf_file_t* file,
                                real_t* latitude_points,
                                int num_latitude_points,
                                const char* latitude_units,
                                real_t* longitude_points,
                                int num_longitude_points,
                                const char* longitude_units,
                                real_t* vertical_points,
                                int num_vertical_points,
                                const char* vertical_units,
                                const char* vertical_orientation)
{
  ASSERT(!cf_file_has_latlon_grid(file));
  ASSERT(num_latitude_points > 0);
  ASSERT(num_longitude_points > 0);
  ASSERT(num_vertical_points > 0);

  // Latitude dimension, data, metadata.
  int dim_id;
  int err = nc_def_dim(file->file_id, "lat", num_latitude_points, &dim_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lat dimension: %s", nc_strerror(err));
  err = nc_def_var(file->file_id, "lat", NC_REAL, 1, &dim_id, &file->lat_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lat variable: %s", nc_strerror(err));
  err = nc_put_var(file->file_id, file->lat_id, latitude_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not set lat data: %s", nc_strerror(err));
  put_attribute(file->file_id, file->lat_id, "long_name", "latitude");
  put_attribute(file->file_id, file->lat_id, "standard_name", "latitude");
  put_attribute(file->file_id, file->lat_id, "units", latitude_units);

  // Longitude metadata.
  err = nc_def_dim(file->file_id, "lon", num_longitude_points, &dim_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lon dimension: %s", nc_strerror(err));
  err = nc_def_var(file->file_id, "lon", NC_REAL, 1, &dim_id, &file->lon_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lon variable: %s", nc_strerror(err));
  err = nc_put_var(file->file_id, file->lon_id, longitude_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not set lon data: %s", nc_strerror(err));
  put_attribute(file->file_id, file->lon_id, "long_name", "longitude");
  put_attribute(file->file_id, file->lon_id, "standard_name", "longitude");
  put_attribute(file->file_id, file->lon_id, "units", longitude_units);

  // Vertical metadata.
  err = nc_def_dim(file->file_id, file->lev_name, num_vertical_points, &dim_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lev dimension: %s", nc_strerror(err));
  err = nc_def_var(file->file_id, file->lev_name, NC_REAL, 1, &dim_id, &file->lev_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lev variable: %s", nc_strerror(err));
  err = nc_put_var(file->file_id, file->lev_id, vertical_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not set lev data: %s", nc_strerror(err));
  put_attribute(file->file_id, file->lev_id, "long_name", "longitude");
  put_attribute(file->file_id, file->lev_id, "standard_name", "longitude");
  put_attribute(file->file_id, file->lev_id, "units", vertical_units);
  put_attribute(file->file_id, file->lev_id, "positive", vertical_orientation);
  put_attribute(file->file_id, file->lev_id, "axis", "Z");
}

bool cf_file_has_latlon_grid(cf_file_t* file)
{
  return ((file->lat_id != -1) && (file->lon_id != -1) && (file->lev_id != -1));
}

void cf_file_get_latlon_grid_metadata(cf_file_t* file,
                                      int* num_latitude_points,
                                      char* latitude_units,
                                      int* num_longitude_points,
                                      char* longitude_units,
                                      int* num_vertical_points,
                                      char* vertical_units,
                                      char* vertical_orientation)
{
  ASSERT(cf_file_has_latlon_grid(file));

  // Latitude.
  *num_latitude_points = cf_file_dimension(file, "lat");
  get_first_attribute(file->file_id, file->lat_id, "units", latitude_units);

  // Longitude.
  *num_latitude_points = cf_file_dimension(file, "lon");
  get_first_attribute(file->file_id, file->lon_id, "units", longitude_units);

  // Vertical.
  *num_latitude_points = cf_file_dimension(file, file->lev_name);
  get_first_attribute(file->file_id, file->lev_id, "units", vertical_units);
  get_first_attribute(file->file_id, file->lev_id, "positive", vertical_orientation);
}

void cf_file_get_latlon_points(cf_file_t* file,
                               real_t* latitude_points,
                               real_t* longitude_points,
                               real_t* vertical_points)
{
  ASSERT(cf_file_has_latlon_grid(file));

  // Latitude.
  int err = nc_get_var(file->file_id, file->lat_id, latitude_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_get_latlon_points: Error retrieving latitudes.");

  // Longitude.
  err = nc_get_var(file->file_id, file->lon_id, longitude_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_get_latlon_points: Error retrieving longitudes.");

  // Vertical.
  err = nc_get_var(file->file_id, file->lev_id, vertical_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_get_latlon_points: Error retrieving vertical coordinates.");
}

void cf_file_define_time(cf_file_t* file,
                         const char* time_units,
                         const char* calendar)
{
  ASSERT(file->time_id == -1);

  // Define the (unlimited) time dimension.
  int t_dim_id;
  int err = nc_def_dim(file->file_id, "time", NC_UNLIMITED, &t_dim_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_time: Could not define time dimension: %s", nc_strerror(err));

  // Now set up the time series.
  err = nc_def_var(file->file_id, "time", NC_REAL, 1, &t_dim_id, &file->time_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_time: Error defining time var: %s", nc_strerror(err));

  // Metadata.
  put_attribute(file->file_id, file->time_id, "long_name", "time");
  put_attribute(file->file_id, file->time_id, "short_name", "time");
  put_attribute(file->file_id, file->time_id, "units", time_units);
  put_attribute(file->file_id, file->time_id, "calendar", calendar);
}

bool cf_file_has_time_series(cf_file_t* file)
{
  return (file->time_id != -1);
}

void cf_file_get_time_metadata(cf_file_t* file,
                               char* time_units,
                               char* calendar)
{
  ASSERT(cf_file_has_time_series(file));
  get_first_attribute(file->file_id, file->time_id, "units", time_units);
  get_first_attribute(file->file_id, file->time_id, "calendar", calendar);
}

int cf_file_append_time(cf_file_t* file, real_t t)
{
  ASSERT(cf_file_has_time_series(file));

  size_t size = (size_t)cf_file_num_times(file);
  int err = nc_put_var1(file->file_id, file->time_id, &size, &t);
  if (err != NC_NOERR)
    polymec_error("cf_file_append_time: Error appending time t = %g: %s", t, nc_strerror(err));

  return (int)size;
}

int cf_file_num_times(cf_file_t* file)
{
  if (!cf_file_has_time_series(file))
    return 0;

  // Find the size of the time series.
  size_t size;
  int err = nc_inq_dimlen(file->file_id, file->time_id, &size);
  if (err != NC_NOERR)
    polymec_error("cf_file_append_time: Error finding length of time series: %s", nc_strerror(err));
  return (int)size;
}

void cf_file_get_times(cf_file_t* file, real_t* times)
{
  int err = nc_get_var(file->file_id, file->time_id, times);
  if (err != NC_NOERR)
    polymec_error("cf_file_get_times: Error retrieving times.");
}

void cf_file_define_latlon_var(cf_file_t* file, 
                               const char* var_name,
                               const char* short_name,
                               const char* long_name,
                               const char* units)
{
  ASSERT(cf_file_has_latlon_grid(file));
  ASSERT(!cf_file_has_latlon_var(file, var_name));

  // Define the variable and its dimensions based on whether we have a time 
  // series.
  int lat_dim, lon_dim, lev_dim, var_id;
  int err = nc_inq_vardimid(file->file_id, file->lat_id, &lat_dim);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_var: Error retrieving lat dim: %s", nc_strerror(err));
  err = nc_inq_vardimid(file->file_id, file->lon_id, &lon_dim);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_var: Error retrieving lon dim: %s", nc_strerror(err));
  err = nc_inq_vardimid(file->file_id, file->lev_id, &lev_dim);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_var: Error retrieving lev dim: %s", nc_strerror(err));
  if (cf_file_has_time_series(file))
  {
    int time_dim;
    err = nc_inq_vardimid(file->file_id, file->time_id, &time_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_define_latlon_var: Error retrieving time dim: %s", nc_strerror(err));
    int dims[4] = {time_dim, lev_dim, lat_dim, lon_dim};
    err = nc_def_var(file->file_id, var_name, NC_REAL, 4, dims, &var_id);
    if (err != NC_NOERR)
      polymec_error("cf_file_define_latlon_var: Error defining var %s: %s", var_name, nc_strerror(err));
  }
  else
  {
    int dims[3] = {lev_dim, lat_dim, lon_dim};
    err = nc_def_var(file->file_id, var_name, NC_REAL, 3, dims, &var_id);
    if (err != NC_NOERR)
      polymec_error("cf_file_define_latlon_var: Error defining var %s: %s", var_name, nc_strerror(err));
  }

  // Metadata.
  put_attribute(file->file_id, var_id, "short_name", short_name);
  put_attribute(file->file_id, var_id, "long_name", long_name);
  put_attribute(file->file_id, var_id, "units", units);
}

void cf_file_get_latlon_var_metadata(cf_file_t* file, 
                                     const char* var_name,
                                     char* short_name,
                                     char* long_name,
                                     char* units)
{
  ASSERT(cf_file_has_latlon_grid(file));
  int var_id = var_identifier(file->file_id, var_name);
  get_first_attribute(file->file_id, var_id, "short_name", short_name);
  get_first_attribute(file->file_id, var_id, "long_name", long_name);
  get_first_attribute(file->file_id, var_id, "units", units);
}

bool cf_file_has_latlon_var(cf_file_t* file,
                            const char* var_name)
{
  if (!cf_file_has_latlon_grid(file))
    return false;
  int var_id = var_identifier(file->file_id, var_name);
  return (var_id != -1);
}

void cf_file_write_latlon_var(cf_file_t* file, 
                              const char* var_name,
                              int time_index, 
                              real_t* var_data)
{
  ASSERT(cf_file_has_latlon_var(file, var_name));

  int var_id = var_identifier(file->file_id, var_name);

  // If we don't have a time series, we just write the whole thing.
  if (!cf_file_has_time_series(file))
  {
    int err = nc_put_var(file->file_id, var_id, var_data);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_var: Error writing data for var %s: %s", var_name, nc_strerror(err));
  }
  // Otherwise, we get fancy and write a hyper slice.
  else
  {
    ASSERT(time_index >= 0);
    ASSERT(time_index < cf_file_num_times(file));

    // Size up the dimensions.
    size_t lat_dim, lon_dim, lev_dim;
    int err = nc_inq_dimlen(file->file_id, file->lat_id, &lat_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_var: Error getting lat dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lon_id, &lon_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_var: Error getting lon dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lev_id, &lev_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_var: Error getting lev dimension %s", nc_strerror(err));

    size_t startp[4] = {time_index, 0, 0, 0};
    size_t countp[4] = {1, lev_dim, lat_dim, lon_dim};
    err = nc_put_vara(file->file_id, var_id, startp, countp, var_data);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_var: Error writing data for var %s: %s", var_name, nc_strerror(err));
  }
}

void cf_file_read_latlon_var(cf_file_t* file, 
                             const char* var_name,
                             int time_index, 
                             real_t* var_data)
{
  ASSERT(cf_file_has_latlon_var(file, var_name));

  int var_id = var_identifier(file->file_id, var_name);

  // If we don't have a time series, we just write the whole thing.
  if (!cf_file_has_time_series(file))
  {
    int err = nc_get_var(file->file_id, var_id, var_data);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_var: Error reading data for var %s: %s", var_name, nc_strerror(err));
  }
  // Otherwise, we get fancy and read a hyper slice.
  else
  {
    ASSERT(time_index >= 0);
    ASSERT(time_index < cf_file_num_times(file));

    // Size up the dimensions.
    size_t lat_dim, lon_dim, lev_dim;
    int err = nc_inq_dimlen(file->file_id, file->lat_id, &lat_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_var: Error getting lat dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lon_id, &lon_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_var: Error getting lon dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lev_id, &lev_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_var: Error getting lev dimension %s", nc_strerror(err));

    size_t startp[4] = {time_index, 0, 0, 0};
    size_t countp[4] = {1, lev_dim, lat_dim, lon_dim};
    err = nc_get_vara(file->file_id, var_id, startp, countp, var_data);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_var: Error writing data for var %s: %s", var_name, nc_strerror(err));
  }
}

void cf_file_define_latlon_surface_var(cf_file_t* file, 
                                       const char* var_name,
                                       const char* short_name,
                                       const char* long_name,
                                       const char* units)
{
  ASSERT(cf_file_has_latlon_grid(file));
  ASSERT(!cf_file_has_latlon_surface_var(file, var_name));

  // Define the variable and its dimensions based on whether we have a time 
  // series.
  int lat_dim, lon_dim, var_id;
  int err = nc_inq_vardimid(file->file_id, file->lat_id, &lat_dim);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_surface_var: Error retrieving lat dim: %s", nc_strerror(err));
  err = nc_inq_vardimid(file->file_id, file->lon_id, &lon_dim);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_surface_var: Error retrieving lon dim: %s", nc_strerror(err));
  if (cf_file_has_time_series(file))
  {
    int time_dim;
    err = nc_inq_vardimid(file->file_id, file->time_id, &time_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_define_latlon_surface_var: Error retrieving time dim: %s", nc_strerror(err));
    int dims[3] = {time_dim, lat_dim, lon_dim};
    err = nc_def_var(file->file_id, var_name, NC_REAL, 3, dims, &var_id);
    if (err != NC_NOERR)
      polymec_error("cf_file_define_latlon_surface_var: Error defining var %s: %s", var_name, nc_strerror(err));
  }
  else
  {
    int dims[2] = {lat_dim, lon_dim};
    err = nc_def_var(file->file_id, var_name, NC_REAL, 2, dims, &var_id);
    if (err != NC_NOERR)
      polymec_error("cf_file_define_latlon_surface_var: Error defining var %s: %s", var_name, nc_strerror(err));
  }

  // Metadata.
  put_attribute(file->file_id, var_id, "short_name", short_name);
  put_attribute(file->file_id, var_id, "long_name", long_name);
  put_attribute(file->file_id, var_id, "units", units);
}

void cf_file_get_latlon_surface_var_metadata(cf_file_t* file, 
                                             const char* var_name,
                                             char* short_name,
                                             char* long_name,
                                             char* units)
{
  cf_file_get_latlon_var_metadata(file, var_name, short_name, long_name, units);
}

bool cf_file_has_latlon_surface_var(cf_file_t* file,
                                    const char* var_name)
{
  return cf_file_has_latlon_var(file, var_name);
}

void cf_file_write_latlon_surface_var(cf_file_t* file, 
                                      const char* var_name,
                                      int time_index, 
                                      real_t* var_data)
{
  ASSERT(cf_file_has_latlon_surface_var(file, var_name));

  int var_id = var_identifier(file->file_id, var_name);

  // If we don't have a time series, we just write the whole thing.
  if (!cf_file_has_time_series(file))
  {
    int err = nc_put_var(file->file_id, var_id, var_data);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_surface_var: Error writing data for var %s: %s", var_name, nc_strerror(err));
  }
  // Otherwise, we get fancy and write a hyper slice.
  else
  {
    ASSERT(time_index >= 0);
    ASSERT(time_index < cf_file_num_times(file));

    // Size up the dimensions.
    size_t lat_dim, lon_dim;
    int err = nc_inq_dimlen(file->file_id, file->lat_id, &lat_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_surface_var: Error getting lat dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lon_id, &lon_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_surface_var: Error getting lon dimension %s", nc_strerror(err));

    size_t startp[4] = {time_index, 0, 0};
    size_t countp[4] = {1, lat_dim, lon_dim};
    err = nc_put_vara(file->file_id, var_id, startp, countp, var_data);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_surface_var: Error writing data for var %s: %s", var_name, nc_strerror(err));
  }
}

void cf_file_read_latlon_surface_var(cf_file_t* file, 
                                     const char* var_name,
                                     int time_index, 
                                     real_t* var_data)
{
  ASSERT(cf_file_has_latlon_surface_var(file, var_name));

  int var_id = var_identifier(file->file_id, var_name);

  // If we don't have a time series, we just write the whole thing.
  if (!cf_file_has_time_series(file))
  {
    int err = nc_get_var(file->file_id, var_id, var_data);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_surface_var: Error reading data for var %s: %s", var_name, nc_strerror(err));
  }
  // Otherwise, we get fancy and read a hyper slice.
  else
  {
    ASSERT(time_index >= 0);
    ASSERT(time_index < cf_file_num_times(file));

    // Size up the dimensions.
    size_t lat_dim, lon_dim;
    int err = nc_inq_dimlen(file->file_id, file->lat_id, &lat_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_surface_var: Error getting lat dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lon_id, &lon_dim);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_surface_var: Error getting lon dimension %s", nc_strerror(err));

    size_t startp[3] = {time_index, 0, 0};
    size_t countp[3] = {1, lat_dim, lon_dim};
    err = nc_get_vara(file->file_id, var_id, startp, countp, var_data);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_surface_var: Error writing data for var %s: %s", var_name, nc_strerror(err));
  }
}

