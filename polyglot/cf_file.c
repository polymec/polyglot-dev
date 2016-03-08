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
};

// Helpers.
static void get_first_attribute(int file_id, 
                                int var_id, 
                                const char* attr,
                                char* value)
{
  size_t attr_len;
  int err = nc_inq_attlen(file_id, var_id, attr, &attr_len);
  if (err != NC_NOERR)
  {
    polymec_error("cf_file: Error retrieving length of attribute %s: %s", 
                  attr, nc_strerror(err));
  }
  if (attr_len > 1)
  {
    char attrs[attr_len][NC_MAX_NAME+1];
    int err = nc_get_att_string(file_id, var_id, attr, (char**)&attrs[0]);
    if (err != NC_NOERR)
    {
      polymec_error("cf_file: Error retrieving attribute %s: %s", 
                    attr, nc_strerror(err));
    }
  }
  else
  {
    int err = nc_get_att_string(file_id, var_id, attr, &value);
    if (err != NC_NOERR)
    {
      polymec_error("cf_file: Error retrieving attribute %s: %s", 
                    attr, nc_strerror(err));
    }
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
  int err = nc_put_att_string(file_id, var_id, attr, 1, &value);
  if (err != NC_NOERR)
  {
    polymec_error("cf_file: Error setting attribute %s: %s", 
                  attr, nc_strerror(err));
  }
}

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

  // Write in our conventions.
  char conventions[NC_MAX_NAME+1];
  snprintf(conventions, NC_MAX_NAME, "CF-%d.%d.%d", cf->cf_major_version, 
           cf->cf_minor_version, cf->cf_patch_version);
  err = nc_put_att_string(cf->file_id, NC_GLOBAL, "Conventions", 1, (const char**)&conventions);
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
  ASSERT(num_latitude_points > 0);
  ASSERT(num_longitude_points > 0);
  ASSERT(num_vertical_points > 0);

  // Latitude dimension, data, metadata.
  int dim_id, var_id;
  int err = nc_def_dim(file->file_id, "lat", num_latitude_points, &dim_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lat dimension: %s", nc_strerror(err));
  err = nc_def_var(file->file_id, "lat", NC_REAL, 1, &dim_id, &var_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lat variable: %s", nc_strerror(err));
  err = nc_put_var(file->file_id, var_id, latitude_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not set lat data: %s", nc_strerror(err));
  put_attribute(file->file_id, var_id, "long_name", "latitude");
  put_attribute(file->file_id, var_id, "standard_name", "latitude");
  put_attribute(file->file_id, var_id, "units", latitude_units);

  // Longitude metadata.
  err = nc_def_dim(file->file_id, "lon", num_longitude_points, &dim_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lon dimension: %s", nc_strerror(err));
  err = nc_def_var(file->file_id, "lon", NC_REAL, 1, &dim_id, &var_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lon variable: %s", nc_strerror(err));
  err = nc_put_var(file->file_id, var_id, longitude_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not set lon data: %s", nc_strerror(err));
  put_attribute(file->file_id, var_id, "long_name", "longitude");
  put_attribute(file->file_id, var_id, "standard_name", "longitude");
  put_attribute(file->file_id, var_id, "units", longitude_units);

  // Vertical metadata.
  err = nc_def_dim(file->file_id, "lev", num_vertical_points, &dim_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lev dimension: %s", nc_strerror(err));
  err = nc_def_var(file->file_id, "lev", NC_REAL, 1, &dim_id, &var_id);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not define lev variable: %s", nc_strerror(err));
  err = nc_put_var(file->file_id, var_id, vertical_points);
  if (err != NC_NOERR)
    polymec_error("cf_file_define_latlon_grid: Could not set lev data: %s", nc_strerror(err));
  put_attribute(file->file_id, var_id, "long_name", "longitude");
  put_attribute(file->file_id, var_id, "standard_name", "longitude");
  put_attribute(file->file_id, var_id, "units", vertical_units);
  put_attribute(file->file_id, var_id, "positive", vertical_orientation);
  put_attribute(file->file_id, var_id, "axis", "Z");
}

bool cf_file_has_latlon_grid(cf_file_t* file)
{
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
}

void cf_file_get_latlon_points(cf_file_t* file,
                               real_t* latitude_points,
                               real_t* longitude_points,
                               real_t* vertical_points)
{
}

void cf_file_define_time(cf_file_t* file,
                         const char* time_units,
                         const char* calendar)
{
}

void cf_file_get_time_metadata(cf_file_t* file,
                               char* time_units,
                               char* calendar)
{
}

int cf_file_set_time(cf_file_t* file, real_t t)
{
}

void cf_file_define_latlon_var(cf_file_t* file, 
                               const char* var_name,
                               const char* short_name,
                               const char* long_name,
                               const char* units)
{
}

void cf_file_get_latlon_var_metadata(cf_file_t* file, 
                                     const char* var_name,
                                     char* short_name,
                                     char* long_name,
                                     char* units)
{
}

void cf_file_write_latlon_var(cf_file_t* file, 
                              const char* var_name,
                              int time_index, 
                              real_t* var_data)
{
}

void cf_file_read_latlon_var(cf_file_t* file, 
                             const char* var_name,
                             int time_index, 
                             real_t* var_data)
{
}

