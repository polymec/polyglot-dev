// Copyright (c) 2015-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYGLOT_CF_FILE_H
#define POLYGLOT_CF_FILE_H

#include "polyglot/polyglot.h"

// The CF file class provides an interface for reading and writing NetCDF 
// files adhering to the Climate/Forecast conventions, the details of which 
// are available at http://cfconventions.org/.

// Maximum length of an identifier name in CF files. Shouldn't exceed name 
// limits in the underlying NetCDF library.
#define POLYGLOT_CF_MAX_NAME 256

// This type provides the interface for CF files.
typedef struct cf_file_t cf_file_t;

// Opens a new CF file for writing simulation data, returning the CF 
// file object. 
cf_file_t* cf_file_new(const char* filename);

// Opens an existing CF file for reading simulation data, returning the CF 
// file object, or NULL if the file is not found or can't be opened, OR if 
// the file is a NetCDF file that doesn't follow the CF conventions.
cf_file_t* cf_file_open(const char* filename);

// Closes and destroys the given CF file handle. If the CF file was opened 
// for writing, this flushes all buffers to disk.
void cf_file_close(cf_file_t* file);

// Retrieves the major, minor, and patch version of the CF conventions 
// that the contents of this file follow. 
void cf_file_get_version(cf_file_t* file, 
                         int* major_version,
                         int* minor_version,
                         int* patch_version);

// Retrieves strings identifying the provenance and contents of this CF file. 
// These identifiers are described at cfconventions.org:
//   title - a succinct description of what is in the dataset.
//   institution - specifies where the original data was produced.
//   source - the method of production of the original data. If it was 
//            model-generated, source should name the model and its version, 
//            as specifically as could be useful. If it is observational, 
//            source should characterize it (e.g., "surface observation" or 
//            "radiosonde").
//   history - provides an audit trail for modifications to the original data. 
//             Well-behaved generic netCDF filters will automatically append 
//             their name and the parameters with which they were invoked to 
//             the global history attribute of an input netCDF file. We 
//             recommend that each line begin with a timestamp indicating the 
//             date and time of day that the program was executed.
//   references - published or web-based references that describe the data or 
//                methods used to produce it.
//   comment - miscellaneous information about the data or methods used to 
//             produce it.
// Each of these strings must be able to hold POLYGLOT_CF_MAX_NAME+1 
// characters. If a given value is not found, its string will be set to the 
// empty string.
void cf_file_get_provenance(cf_file_t* file, 
                            char* title,
                            char* institution,
                            char* source,
                            char* history,
                            char* references,
                            char* comment);

// Sets the provenance information (described above) for the given file, which 
// must have been opened for writing.
void cf_file_set_provenance(cf_file_t* file, 
                            const char* title,
                            const char* institution,
                            const char* source,
                            const char* history,
                            const char* references,
                            const char* comment);

// Sets the value of the given global attribute.
void cf_file_set_global_attribute(cf_file_t* file, 
                                  const char* global_attribute_name,
                                  const char* value);

// Fetches the value of the given global attribute into value, which must be 
// able to hold POLYGLOT_CF_MAX_NAME+1 characters.
void cf_file_get_global_attribute(cf_file_t* file, 
                                  const char* global_attribute_name,
                                  char* value);

// Sets up a generic dimension with the given numeric value, or -1 if 
// the dimension is considered to be "unlimited" (like a time series).
void cf_file_define_dimension(cf_file_t* file,
                              const char* dimension_name,
                              int value);

// Retrieves the value of the given dimension. -1 corresponds to an 
// "unlimited" dimension.
int cf_file_dimension(cf_file_t* file, const char* dimension_name);

// Sets up the metadata for a latitude/longitude grid with the given number 
// of latitudinal/longitudinal/vertical grid points and units of measure.
// Possible units of measure are:
//   latitude: degree_north, degree_N, degrees_N, degreeN, degreesN
//   longitude: degree_east, degree_E, degrees_E, degreeE, degreesE
//   vertical: bar, millibar, decibar, atmosphere (atm), pascal (Pa), hPa,
//             meter (metre, m), kilometer (km), 
//             level, layer, sigma_level
// Additionally, an orientiation (vertical_orientation, "up" or "down") must 
// be specified.
void cf_file_define_latlon_grid(cf_file_t* file,
                                int num_latitude_points,
                                const char* latitude_units,
                                int num_longitude_points,
                                const char* longitude_units,
                                int num_vertical_points,
                                const char* vertical_units,
                                const char* vertical_orientation);

// Returns true if the CF file contains a lat-lon grid, false if not.
bool cf_file_has_latlon_grid(cf_file_t* file);

// Fetches information about the latlon grid in the given file (if present).
// Strings should all be able to store POLYGLOT_CF_MAX_NAME+1 characters.
void cf_file_get_latlon_grid_metadata(cf_file_t* file,
                                      int* num_latitude_points,
                                      char* latitude_units,
                                      int* num_longitude_points,
                                      char* longitude_units,
                                      int* num_vertical_points,
                                      char* vertical_units,
                                      char* vertical_orientation);

// Sets up the time series variable within the CF file, specifying units and 
// a calendar to use. Available units are:
//   day (d), hour (hr, h), minute (min), second (sec, s)
// Plural forms may be used as well, and "since xyz" where xyz is a timestamp
// may be used to specify times relative to an epoch. Available calendars are:
//   gregorian (standard), proleptic_gregorian, noleap (365_day), 
//   all_leap (366_day), 360_day, julian, none
void cf_file_define_time(cf_file_t* file,
                         const char* time_units,
                         const char* calendar);

// Returns true if this file contains a time series, false if not.
bool cf_file_has_time_series(cf_file_t* file);

// Returns the number of time in the file's time series, or 0 if the file has 
// no time series.
int cf_file_num_times(cf_file_t* file);

// Retrieves time information (units and calendar) to strings large enough to 
// hold NC_NAME_MAX+1 characters.
void cf_file_get_time_metadata(cf_file_t* file,
                               char* time_units,
                               char* calendar);

// Defines a (3D) variable that is defined on the points of a lat-lon grid, 
// setting up metadata like short and long names and units. If the variable 
// is time-dependent, its dimensions will be (time, vertical, lat, lon); 
// otherwise they will be (vertical, lat, lon).
void cf_file_define_latlon_var(cf_file_t* file, 
                               const char* var_name,
                               bool time_dependent,
                               const char* short_name,
                               const char* long_name,
                               const char* units);

// Fetches metadata for the given lat-lon variable. All strings must 
// be large enough to hold POLYGLOT_CF_MAX_NAME+1 characters. 
void cf_file_get_latlon_var_metadata(cf_file_t* file, 
                                     const char* var_name,
                                     char* short_name,
                                     char* long_name,
                                     char* units);

// Returns true if this file contains a lat-lon variable with the given name,
// false otherwise.
bool cf_file_has_latlon_var(cf_file_t* file,
                            const char* var_name);

// Defines a 2D surface variable that is defined on the points of a lat-lon 
// grid, setting up metadata like short and long names and units. If the 
// variable is time-dependent, its dimensions will be (time, lat, lon); 
// otherwise they will be (lat, lon).
void cf_file_define_latlon_surface_var(cf_file_t* file, 
                                       const char* var_name,
                                       bool time_dependent,
                                       const char* short_name,
                                       const char* long_name,
                                       const char* units);

// Fetches metadata for the given lat-lon surface variable. All strings must 
// be large enough to hold POLYGLOT_CF_MAX_NAME+1 characters. 
void cf_file_get_latlon_surface_var_metadata(cf_file_t* file, 
                                             const char* var_name,
                                             char* short_name,
                                             char* long_name,
                                             char* units);

// Returns true if this file contains a lat-lon variable with the given name,
// false otherwise.
bool cf_file_has_latlon_surface_var(cf_file_t* file,
                                    const char* var_name);

// Writes the grid's coordinate data to the file.
void cf_file_write_latlon_grid(cf_file_t* file,
                               real_t* latitude_points,
                               real_t* longitude_points,
                               real_t* vertical_points);

// Fetches latitude/longitude/vertical grid points to the given arrays, which 
// must be large enough to fit them.
void cf_file_read_latlon_grid(cf_file_t* file,
                              real_t* latitude_points,
                              real_t* longitude_points,
                              real_t* vertical_points);

// Appends a time to the time series in the grid, returning
// an integer index identifying that time.
int cf_file_append_time(cf_file_t* file, real_t t);

// Writes a variable that is defined on the points of a lat-lon grid, 
// specifying a time index that associates this entry with a given time. This 
// time index is ignored if the variable is not time dependent.
void cf_file_write_latlon_var(cf_file_t* file, 
                              const char* var_name,
                              int time_index, 
                              real_t* var_data);

// Reads a variable that is defined on the points of a lat-lon grid, 
// specifying an index for the time at which the data will be read. This 
// time index is ignored if the file has no time series.
void cf_file_read_latlon_var(cf_file_t* file, 
                             const char* var_name,
                             int time_index, 
                             real_t* var_data);

// Writes a surface variable that is defined on the points of a lat-lon grid, 
// specifying a time index that associates this entry with a given time. This 
// time index is ignored if the variable is not time-dependent.
void cf_file_write_latlon_surface_var(cf_file_t* file, 
                                      const char* var_name,
                                      int time_index, 
                                      real_t* var_data);

// Reads a variable that is defined on the surface of a lat-lon grid, 
// specifying an index for the time at which the data will be read. This 
// time index is ignored if the variable is not time-dependent.
void cf_file_read_latlon_surface_var(cf_file_t* file, 
                                     const char* var_name,
                                     int time_index, 
                                     real_t* var_data);

#endif
