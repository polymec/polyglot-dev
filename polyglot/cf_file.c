// Copyright (c) 2015-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "netcdf.h"
#include "core/unordered_map.h"
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
  int time_id, time_dim, lat_id, lat_dim, lon_id, lon_dim, lev_id, lev_dim;
  char lev_name[POLYGLOT_CF_MAX_NAME+1];

  // Lat-lon variable metadata/indices.
  int nlat, nlon, nlev;
  string_int_unordered_map_t *ll_vars, *td_ll_vars;
  string_int_unordered_map_t *ll_surface_vars, *td_ll_surface_vars;
};

// Helpers.
static void get_first_attribute(int file_id, 
                                int var_id, 
                                const char* attr,
                                char* value)
{
  size_t len;
  int err = nc_inq_attlen(file_id, var_id, attr, &len);
  if (err == NC_ENOTATT)
  {
    value[0] = '\0';
    return;
  }
  else if (err != NC_NOERR)
  {
    polymec_error("cf_file: Error retrieving attribute %s length: %s", 
                  attr, nc_strerror(err));
  }
  err = nc_get_att_text(file_id, var_id, attr, value);
  if (err == NC_ENOTATT)
    value[0] = '\0';
  else if (err != NC_NOERR)
  {
    polymec_error("cf_file: Error retrieving attribute %s: %s", 
                  attr, nc_strerror(err));
  }
  value[len] = '\0';
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
  int err = nc_put_att_text(file_id, var_id, attr, strlen(value), value);
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

static void find_vertical_coordinate(int file_id, int* lev_id, int* lev_dim, char* lev_name)
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
      *lev_dim = dim_id;
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
        ((string_casecmp(positive, "up") == 0) || 
         (string_casecmp(positive, "down") == 0)))
    {
      *lev_id = var_id;
      strcpy(lev_name, var_name);
      return;
    }

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
  polymec_error("Could not identify vertical coordinate from file metadata.");
}

// Implementation.

cf_file_t* cf_file_new(const char* filename)
{
  int file_id;
  int err = nc_create(filename, NC_CLOBBER | NC_NETCDF4, &file_id);
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
  cf->time_dim = cf->lat_dim = cf->lon_dim = cf->lev_dim = -1;
  strcpy(cf->lev_name, "lev");
  cf->nlat = cf->nlon = cf->nlev = -1;
  cf->ll_vars = string_int_unordered_map_new();
  cf->td_ll_vars = string_int_unordered_map_new();
  cf->ll_surface_vars = string_int_unordered_map_new();
  cf->td_ll_surface_vars = string_int_unordered_map_new();

  // Write in our conventions.
  char conventions[NC_MAX_NAME+1];
  snprintf(conventions, NC_MAX_NAME, "CF-%d.%d.%d", cf->cf_major_version, 
           cf->cf_minor_version, cf->cf_patch_version);
  err = nc_put_att_text(cf->file_id, NC_GLOBAL, "Conventions", strlen(conventions), conventions);
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
  cf->time_id = cf->lat_id = cf->lon_id = cf->lev_id = -1;
  cf->time_dim = cf->lat_dim = cf->lon_dim = cf->lev_dim = -1;
  cf->nlat = cf->nlon = cf->nlev = -1;
  cf->ll_vars = string_int_unordered_map_new();
  cf->td_ll_vars = string_int_unordered_map_new();
  cf->ll_surface_vars = string_int_unordered_map_new();
  cf->td_ll_surface_vars = string_int_unordered_map_new();

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
  {
    cf->lat_dim = dim_id;
    cf->lat_id = var_identifier(cf->file_id, "lat");
  }
  else if (err != NC_EBADDIM)
    polymec_error("cf_file_open: Error retrieving lat dim ID: ", nc_strerror(err));

  err = nc_inq_dimid(cf->file_id, "lon", &dim_id);
  if (err == NC_NOERR)
  {
    cf->lon_dim = dim_id;
    cf->lon_id = var_identifier(cf->file_id, "lon");
  }
  else if (err != NC_EBADDIM)
    polymec_error("cf_file_open: Error retrieving lon dim ID: ", nc_strerror(err));

  err = nc_inq_dimid(cf->file_id, "time", &dim_id);
  if (err == NC_NOERR)
  {
    cf->time_dim = dim_id;
    cf->time_id = var_identifier(cf->file_id, "time");
  }
  else if (err != NC_EBADDIM)
    polymec_error("cf_file_open: Error retrieving time dim ID: ", nc_strerror(err));

  // If we've found a lat/lon grid, feel out the data related to it.
  if ((cf->lat_id != -1) && (cf->lon_id != -1))
  {
    // We have to figure out the vertical dimension / coordinate name and ID. 
    find_vertical_coordinate(cf->file_id, &cf->lev_id, &cf->lev_dim, cf->lev_name);

    // Get the dimensions of the lat/lon grid.
    cf->nlat = cf_file_dimension(cf, "lat");
    cf->nlon = cf_file_dimension(cf, "lon");
    cf->nlev = cf_file_dimension(cf, cf->lev_name);

    // Get all of the lat/lon variables we can find.
    int nvarsp;
    err = nc_inq(file_id, NULL, &nvarsp, NULL, NULL);
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

      // Find the dimension ids.
      int dim_ids[ndim];
      err = nc_inq_vardimid(file_id, var_id, dim_ids);
      if (err != NC_NOERR)
        polymec_error("cf_file_open: Error retrieving dim IDs for var %d: ", var_id, nc_strerror(err));

      // Now we can determine the nature of the variable by inspecting its dimensions.
      if ((ndim == 2) && (dim_ids[0] == cf->lat_dim) && (dim_ids[1] == cf->lon_dim))
        string_int_unordered_map_insert_with_k_dtor(cf->ll_surface_vars, string_dup(var_name), var_id, string_free);
      else if ((ndim == 3) && (dim_ids[0] == cf->time_dim) && (dim_ids[1] == cf->lat_dim) && (dim_ids[2] == cf->lon_dim))
        string_int_unordered_map_insert_with_k_dtor(cf->td_ll_surface_vars, string_dup(var_name), var_id, string_free);
      else if ((ndim == 3) && (dim_ids[0] == cf->lev_dim) && (dim_ids[1] == cf->lat_dim) && (dim_ids[2] == cf->lon_dim))
        string_int_unordered_map_insert_with_k_dtor(cf->ll_vars, string_dup(var_name), var_id, string_free);
      else if ((ndim == 4) && (dim_ids[0] == cf->time_dim) && (dim_ids[1] == cf->lev_dim) && (dim_ids[2] == cf->lat_dim) && (dim_ids[3] == cf->lon_dim))
        string_int_unordered_map_insert_with_k_dtor(cf->td_ll_vars, string_dup(var_name), var_id, string_free);
    }
  }

  return cf;
}

void cf_file_close(cf_file_t* file)
{
  int err = nc_close(file->file_id);
  if (err != NC_NOERR)
    polymec_error("Error closing CF file.", nc_strerror(err));
  string_int_unordered_map_free(file->ll_vars);
  string_int_unordered_map_free(file->td_ll_vars);
  string_int_unordered_map_free(file->ll_surface_vars);
  string_int_unordered_map_free(file->td_ll_surface_vars);
  polymec_free(file);
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
                                int num_latitude_points,
                                const char* latitude_units,
                                int num_longitude_points,
                                const char* longitude_units,
                                int num_vertical_points,
                                const char* vertical_units,
                                const char* vertical_orientation)
{
  ASSERT(!cf_file_has_latlon_grid(file));
  ASSERT(num_latitude_points > 0);
  ASSERT(num_longitude_points > 0);
  ASSERT(num_vertical_points > 0);

  // Validate units and orientation.
  const char* valid_lat_units[] = {"degree_north", "degree_N", "degrees_N",
                                   "degreeN", "degreesN", NULL};
  int index = string_find_in_list(latitude_units, valid_lat_units, true);
  if (index == -1)
    polymec_error("cf_file_define_latlon_grid: Invalid latitude units: %s", latitude_units);

  const char* valid_lon_units[] = {"degree_east", "degree_E", "degrees_E",
                                   "degreeE", "degreesE", NULL};
  index = string_find_in_list(longitude_units, valid_lon_units, true);
  if (index == -1)
    polymec_error("cf_file_define_latlon_grid: Invalid longitude units: %s", longitude_units);

  const char* valid_vert_units[] = {"bar", "millibar", "decibar",
                                    "atmosphere", "atm", "pascal", "Pa", "hPa", 
                                    "meter", "metre", "m", "kilometer", "km", 
                                    "level", "sigma_level", NULL};
  index = string_find_in_list(vertical_units, valid_vert_units, true);
  if (index == -1)
    polymec_error("cf_file_define_latlon_grid: Invalid vertical units: %s", vertical_units);

  const char* valid_vert_orientation[] = {"up", "down", NULL};
  index = string_find_in_list(vertical_orientation, valid_vert_orientation, false);
  if (index == -1)
    polymec_error("cf_file_define_latlon_grid: Invalid vertical orientation: %s", vertical_orientation);

  // Latitude dimension, data, metadata.
  int err = nc_def_dim(file->file_id, "lat", num_latitude_points, &file->lat_dim);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lat dimension: %s", nc_strerror(err));
  err = nc_def_var(file->file_id, "lat", NC_REAL, 1, &file->lat_dim, &file->lat_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lat variable: %s", nc_strerror(err));
  put_attribute(file->file_id, file->lat_id, "long_name", "latitude");
  put_attribute(file->file_id, file->lat_id, "standard_name", "latitude");
  put_attribute(file->file_id, file->lat_id, "units", latitude_units);
  file->nlat = num_latitude_points;

  // Longitude metadata.
  err = nc_def_dim(file->file_id, "lon", num_longitude_points, &file->lon_dim);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lon dimension: %s", nc_strerror(err));
  err = nc_def_var(file->file_id, "lon", NC_REAL, 1, &file->lon_dim, &file->lon_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lon variable: %s", nc_strerror(err));
  put_attribute(file->file_id, file->lon_id, "long_name", "longitude");
  put_attribute(file->file_id, file->lon_id, "standard_name", "longitude");
  put_attribute(file->file_id, file->lon_id, "units", longitude_units);
  file->nlon = num_longitude_points;

  // Vertical metadata.
  err = nc_def_dim(file->file_id, file->lev_name, num_vertical_points, &file->lev_dim);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lev dimension: %s", nc_strerror(err));
  err = nc_def_var(file->file_id, file->lev_name, NC_REAL, 1, &file->lev_dim, &file->lev_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lev variable: %s", nc_strerror(err));
  put_attribute(file->file_id, file->lev_id, "units", vertical_units);
  put_attribute(file->file_id, file->lev_id, "positive", vertical_orientation);
  put_attribute(file->file_id, file->lev_id, "axis", "Z");
  file->nlev = num_vertical_points;
}

void cf_file_write_latlon_grid(cf_file_t* file,
                               real_t* latitude_points,
                               real_t* longitude_points,
                               real_t* vertical_points)
{
  int err = nc_put_var(file->file_id, file->lat_id, latitude_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_write_latlon_grid: Could not set lat data: %s", nc_strerror(err));
  err = nc_put_var(file->file_id, file->lon_id, longitude_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_write_latlon_grid: Could not set lon data: %s", nc_strerror(err));
  err = nc_put_var(file->file_id, file->lev_id, vertical_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_write_latlon_grid: Could not set lev data: %s", nc_strerror(err));

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
  *num_latitude_points = file->nlat;
  get_first_attribute(file->file_id, file->lat_id, "units", latitude_units);

  // Longitude.
  *num_latitude_points = file->nlon;
  get_first_attribute(file->file_id, file->lon_id, "units", longitude_units);

  // Vertical.
  *num_latitude_points = file->nlev;
  get_first_attribute(file->file_id, file->lev_id, "units", vertical_units);
  get_first_attribute(file->file_id, file->lev_id, "positive", vertical_orientation);
}

void cf_file_read_latlon_grid(cf_file_t* file,
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
  int err = nc_def_dim(file->file_id, "time", NC_UNLIMITED, &file->time_dim);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_time: Could not define time dimension: %s", nc_strerror(err));

  // Now set up the time series.
  err = nc_def_var(file->file_id, "time", NC_REAL, 1, &file->time_dim, &file->time_id);
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
                               bool time_dependent,
                               const char* short_name,
                               const char* long_name,
                               const char* units)
{
  ASSERT(cf_file_has_latlon_grid(file));
  ASSERT(!cf_file_has_latlon_var(file, var_name));

  // Define the variable and its dimensions based on whether we have a time 
  // series.
  int var_id;
  if (time_dependent)
  {
    ASSERT(cf_file_has_time_series(file));
    int dims[4] = {file->time_dim, file->lev_dim, file->lat_dim, file->lon_dim};
    int err = nc_def_var(file->file_id, var_name, NC_REAL, 4, dims, &var_id);
    if (err != NC_NOERR)
      polymec_error("cf_file_define_latlon_var: Error defining var %s: %s", var_name, nc_strerror(err));
    string_int_unordered_map_insert_with_k_dtor(file->td_ll_vars, string_dup(var_name), var_id, string_free);
  }
  else
  {
    int dims[3] = {file->lev_dim, file->lat_dim, file->lon_dim};
    int err = nc_def_var(file->file_id, var_name, NC_REAL, 3, dims, &var_id);
    if (err != NC_NOERR)
      polymec_error("cf_file_define_latlon_var: Error defining var %s: %s", var_name, nc_strerror(err));
    string_int_unordered_map_insert_with_k_dtor(file->ll_vars, string_dup(var_name), var_id, string_free);
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
  return (string_int_unordered_map_contains(file->ll_vars, (char*)var_name) ||
          string_int_unordered_map_contains(file->td_ll_vars, (char*)var_name));
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
    size_t lat_len, lon_len, lev_len;
    int err = nc_inq_dimlen(file->file_id, file->lat_dim, &lat_len);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_var: Error getting lat dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lon_dim, &lon_len);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_var: Error getting lon dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lev_dim, &lev_len);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_var: Error getting lev dimension %s", nc_strerror(err));

    size_t startp[4] = {time_index, 0, 0, 0};
    size_t countp[4] = {1, lev_len, lat_len, lon_len};
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

  bool time_dependent = false;
  int var_id;
  int* var_id_p = string_int_unordered_map_get(file->td_ll_vars, (char*)var_name);
  if (var_id_p != NULL)
    var_id = *var_id_p; 
  else
  {
    var_id_p = string_int_unordered_map_get(file->ll_vars, (char*)var_name);
    ASSERT(var_id_p != NULL);
    time_dependent = true;
    var_id = *var_id_p;
  }

  // If we don't have a time series, we just write the whole thing.
  if (!time_dependent)
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
    size_t lat_len, lon_len, lev_len;
    int err = nc_inq_dimlen(file->file_id, file->lat_dim, &lat_len);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_var: Error getting lat dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lon_dim, &lon_len);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_var: Error getting lon dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lev_dim, &lev_len);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_var: Error getting lev dimension %s", nc_strerror(err));

    size_t startp[4] = {time_index, 0, 0, 0};
    size_t countp[4] = {1, lev_len, lat_len, lon_len};
    err = nc_get_vara(file->file_id, var_id, startp, countp, var_data);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_var: Error writing data for var %s: %s", var_name, nc_strerror(err));
  }
}

void cf_file_define_latlon_surface_var(cf_file_t* file, 
                                       const char* var_name,
                                       bool time_dependent,
                                       const char* short_name,
                                       const char* long_name,
                                       const char* units)
{
  ASSERT(cf_file_has_latlon_grid(file));
  ASSERT(!cf_file_has_latlon_surface_var(file, var_name));

  // Define the variable and its dimensions based on whether we have a time 
  // series.
  int var_id;
  if (time_dependent)
  {
    ASSERT(cf_file_has_time_series(file));

    int dims[3] = {file->time_dim, file->lat_dim, file->lon_dim};
    int err = nc_def_var(file->file_id, var_name, NC_REAL, 3, dims, &var_id);
    if (err != NC_NOERR)
      polymec_error("cf_file_define_latlon_surface_var: Error defining var %s: %s", var_name, nc_strerror(err));
    string_int_unordered_map_insert_with_k_dtor(file->td_ll_surface_vars, string_dup(var_name), var_id, string_free);
  }
  else
  {
    int dims[2] = {file->lat_dim, file->lon_dim};
    int err = nc_def_var(file->file_id, var_name, NC_REAL, 2, dims, &var_id);
    if (err != NC_NOERR)
      polymec_error("cf_file_define_latlon_surface_var: Error defining var %s: %s", var_name, nc_strerror(err));
    string_int_unordered_map_insert_with_k_dtor(file->ll_surface_vars, string_dup(var_name), var_id, string_free);
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
  return (string_int_unordered_map_contains(file->ll_surface_vars, (char*)var_name) ||
          string_int_unordered_map_contains(file->td_ll_surface_vars, (char*)var_name));
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
    size_t lat_len, lon_len;
    int err = nc_inq_dimlen(file->file_id, file->lat_dim, &lat_len);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_surface_var: Error getting lat dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lon_dim, &lon_len);
    if (err != NC_NOERR)
      polymec_error("cf_file_write_latlon_surface_var: Error getting lon dimension %s", nc_strerror(err));

    size_t startp[4] = {time_index, 0, 0};
    size_t countp[4] = {1, lat_len, lon_len};
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

  bool time_dependent = false;
  int var_id;
  int* var_id_p = string_int_unordered_map_get(file->td_ll_vars, (char*)var_name);
  if (var_id_p != NULL)
    var_id = *var_id_p; 
  else
  {
    var_id_p = string_int_unordered_map_get(file->ll_vars, (char*)var_name);
    ASSERT(var_id_p != NULL);
    time_dependent = true;
    var_id = *var_id_p;
  }

  // If we don't have a time series, we just write the whole thing.
  if (!time_dependent)
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
    size_t lat_len, lon_len;
    int err = nc_inq_dimlen(file->file_id, file->lat_dim, &lat_len);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_surface_var: Error getting lat dimension %s", nc_strerror(err));
    err = nc_inq_dimlen(file->file_id, file->lon_dim, &lon_len);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_surface_var: Error getting lon dimension %s", nc_strerror(err));

    size_t startp[3] = {time_index, 0, 0};
    size_t countp[3] = {1, lat_len, lon_len};
    err = nc_get_vara(file->file_id, var_id, startp, countp, var_data);
    if (err != NC_NOERR)
      polymec_error("cf_file_read_latlon_surface_var: Error writing data for var %s: %s", var_name, nc_strerror(err));
  }
}

