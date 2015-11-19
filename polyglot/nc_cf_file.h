// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYGLOT_NC_CF_FILE_H
#define POLYGLOT_NC_CF_FILE_H

#include "polyglot/polyglot.h"

// The NC CF file class provides an interface for reading and writing NetCDF 
// files adhering to the Climate/Forecast conventions, the details of which 
// are available at http://cfconventions.org/.

// This type provides the interface for CF files.
typedef struct nc_cf_file_t nc_cf_file_t;

// Opens a new CF file for reading simulation data, returning the CF 
// file object. 
nc_cf_file_t* nc_cf_file_new(const char* filename);

// Returns true if the CF file was parsed successfully upon construction 
// and contains valid data, false if not.
bool nc_cf_file_is_valid(nc_cf_file_t* file);

// Returns a pointer to an internal string containing any error message(s) 
// related to the parsing of the CF file. If the CF file was parsed 
// and the object is valid, this returns the empty string.
char* nc_cf_file_error(nc_cf_file_t* file);

// Closes and destroys the given CF file.
void nc_cf_file_close(nc_cf_file_t* file);

#endif