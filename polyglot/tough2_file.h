// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYGLOT_TOUGH2_FILE_H
#define POLYGLOT_TOUGH2_FILE_H

#include "polyglot/polyglot.h"
#include "core/point_cloud.h"

// The TOUGH2 file class provides an interface to data within input files 
// written for TOUGH2. This interface is inherently serial because TOUGH2 
// does not use parallel I/O, so you should create TOUGH2 file objects only 
// on a single process. For information on the various blocks within TOUGH2
// files please refer to the TOUGH2 Users Guide.

// This defines the maximum length of element names for TOUGH2 files.
#define TOUGH2_ELEMNAME_MAX 8

// This defines the maximum length of rock type names for TOUGH2 files.
#define TOUGH2_ROCKNAME_MAX 8

// This defines the maximum length of source names for TOUGH2 files.
#define TOUGH2_SOURCENAME_MAX 8

// This identifies the types of relative permeability functions 
// supported by TOUGH2.
typedef enum
{
  TOUGH2_REL_PERM_LINEAR,
  TOUGH2_REL_PERM_EXPONENTIAL,
  TOUGH2_REL_PERM_COREY,
  TOUGH2_REL_PERM_GRANT,
  TOUGH2_REL_PERM_PERFECT_MOBILITY,
  TOUGH2_REL_PERM_FATT_KLICKOFF,
  TOUGH2_REL_PERM_VAN_GENUCHTEN_MUALEM,
  TOUGH2_REL_PERM_VERMA,
} tough2_file_rel_perm_t;

// This identifies the types of capillary pressure functions 
// supported by TOUGH2.
typedef enum
{
  TOUGH2_CAP_PRES_LINEAR,
  TOUGH2_CAP_PRES_PICKENS,
  TOUGH2_CAP_PRES_TRUST,
  TOUGH2_CAP_PRES_MILLY,
  TOUGH2_CAP_PRES_LEVERETT,
  TOUGH2_CAP_PRES_VAN_GENUCHTEN,
  TOUGH2_CAP_PRES_NONE
} tough2_file_cap_pres_t;

// This is a container type for describing material parameters associated 
// with reservoir domains.
typedef struct
{
  char name[TOUGH2_ROCKNAME_MAX];
  real_t density;
  real_t porosity;
  real_t permeability[3];
  real_t sat_heat_conductivity;
  real_t specific_heat;
  real_t pore_compressibility;
  real_t pore_expansivity;
  real_t unsat_heat_conductivity;
  real_t tortuosity_factor;
  real_t klinkenberg_b;
  real_t xkd3, xkd4;
  tough2_file_rel_perm_t rel_perm_type;
  real_t rel_perm_params[7];
  tough2_file_cap_pres_t cap_pres_t;
  real_t cap_pres_params[7];
} tough2_file_rock_t;

// This is a container type for holding phase-specific diffusion coefficients.
typedef struct
{
  // coeffs[i][j] contains the jth diffusion coefficient for the ith phase.
  real_t coeffs[3][8];
} tough2_file_diffusion_t;

// This is a container type for describing sinks and sources. Objects of this 
// type are garbage-collected.
typedef struct
{
  char elem_name[TOUGH2_ELEMNAME_MAX];
  char name[TOUGH2_SOURCENAME_MAX];
  int ltab;
  char type[4];
  char itab;
  real_t generation_rate;
  real_t specific_enthalpy;
  real_t layer_thickness;
  real_t* generation_times;
  real_t* generation_rates;
  real_t* specific_enthalpies;
} tough2_file_source_t;

// This type provides the interface for TOUGH2 files.
typedef struct tough2_file_t tough2_file_t;

// Opens a new TOUGH2 file for reading simulation data, returning the TOUGH 
// file object. This object must be checked for validity.
tough2_file_t* tough2_file_new(const char* filename);

// Returns true if the TOUGH2 file was parsed successfully upon construction 
// and contains valid data, false if not.
bool tough2_file_is_valid(tough2_file_t* file);

// Returns a pointer to an internal string containing any error message(s) 
// related to the parsing of the TOUGH2 file. If the TOUGH2 file was parsed 
// and the object is valid, this returns the empty string.
char* tough2_file_error(tough2_file_t* file);

// Closes and destroys the given TOUGH2 file.
void tough2_file_close(tough2_file_t* file);

// Returns true if the TOUGH2 file contains a block with the given name, 
// false if not.
bool tough2_file_contains_block(tough2_file_t* file, 
                                const char* block_name);

// Reads a mesh from the given TOUGH2 file, returning a newly-allocated 
// point_cloud object with properties and tags corresponding to the 
// elements and connections within. A TOUGH meshe does not contain enough 
// geometric information to be represented as an arbitrary polyhedral mesh, 
// so we store it as a point cloud with an associated neighbor pairing 
// in a property named "connections". If the file does not contain a mesh or 
// a mesh cannot be generated from it, this function returns NULL.
point_cloud_t* tough2_file_read_mesh(tough2_file_t* file);

// Returns true if the given TOUGH2 file contains a mesh, false if not.
bool tough2_file_contains_mesh(tough2_file_t* file);

// Reads material parameters from the TOUGH2 file, returning them as an 
// array of tough2_file_rock_t objects. The number of entries read from the file
// is stored in num_entries. If the file contains no such data, this function 
// returns NULL, and 0 is stored in num_entries.
tough2_file_rock_t* tough2_file_read_rocks(tough2_file_t* file, int* num_entries);

// Returns true if the given TOUGH2 file contains material parameters, false 
// if not.
bool tough2_file_contains_rocks(tough2_file_t* file);

// Parses information in the file about which equations will be solved, 
// filling in the given variables with this information. If such information is 
// not available, no values will be filled in.
void tough2_file_read_multi(tough2_file_t* file,
                            int* nk, int* neq, int* nph, int* nb, int* nkin);

// Returns true if the given TOUGH2 file contains information about which 
// equations to solve, false if not.
bool tough2_file_contains_multi(tough2_file_t* file);

// Parses information in the file about parameters for computation,
// filling in the given variables with this information. If such information is 
// not available, no values will be filled in. Note that mop should be an
// array that can store at least 24 integers; elst must be a character
// array that can store at least TOUGH2_ELEMNAME_MAX characters; dlt must be 
// an array that can store at least 100 reals; dep must be an array that can 
// store at least 4 reals.
void tough2_file_read_param(tough2_file_t* file,
                            int* noite, int* kdata, int* mcyc, 
                            real_t* msec, int* mcypr, int* mop,
                            real_t* texp, real_t* be,
                            real_t* tstart, real_t* timax,
                            real_t* delten, real_t* deltmx,
                            char* elst, real_t* gf, real_t* redlt,
                            real_t* scale, real_t* dlt, 
                            real_t* re1, real_t* re2, real_t* u,
                            real_t* wup, real_t* wnr, real_t* dfac,
                            real_t* dep);

// Returns true if the given TOUGH file contains computation parameters, 
// false if not.
bool tough2_file_contains_param(tough2_file_t* file);

// Parses domain-specific initial conditions (default values for primary 
// variables). mat must be a character array that can store at least 
// TOUGH2_ELEMNAME_MAX characters; X must be an array that can store at least 
// 4 reals.
void tough2_file_read_indom(tough2_file_t* file,
                            char* mat,
                            real_t* X);

// Returns true if the given TOUGH file contains domain-specific initial 
// conditions, false if not.
bool tough2_file_contains_indom(tough2_file_t* file);

// Parses element-specific initial conditions, returning a newly-allocated
// array containing values for primary variables in component-minor order,
// corresponding to the points in the mesh.
real_t* tough2_file_read_incon(tough2_file_t* file);

// Returns true if the given TOUGH file contains element-specific initial 
// conditions, false if not.
bool tough2_file_contains_incon(tough2_file_t* file);

// Returns a newly-allocated array of names of elements for which it is 
// requested that time-dependent data be written out during the simulation,
// or NULL if no such elements are requested.
char** tough2_file_read_foft(tough2_file_t* file);

// Returns true if the given TOUGH file contains names of elements whose data
// should be written out, false if not.
bool tough2_file_contains_foft(tough2_file_t* file);

// Returns a newly-allocated array of names of connections for which it is 
// requested that time-dependent data be written out during the simulation,
// or NULL if no such connections are requested.
char** tough2_file_read_coft(tough2_file_t* file);

// Returns true if the given TOUGH file contains names of connections whose 
// data should be written out, false if not.
bool tough2_file_contains_coft(tough2_file_t* file);

// Returns a newly-allocated array of names of sources/sinks for which it is 
// requested that time-dependent data be written out during the simulation,
// or NULL if no such sources/sinks are requested.
char** tough2_file_read_goft(tough2_file_t* file);

// Returns true if the given TOUGH file contains names of sources/sinks whose 
// data should be written out, false if not.
bool tough2_file_contains_goft(tough2_file_t* file);

// Reads diffusion coefficients from the TOUGH2 file, returning them as an 
// array of tough2_file_diffusion_t objects. The number of entries read from the 
// file is stored in num_entries. If the file contains no such data, this 
// function returns NULL, and 0 is stored in num_entries.
tough2_file_diffusion_t* tough2_file_read_diffu(tough2_file_t* file, int* num_entries);

// Returns true if the given TOUGH2 file contains diffusion coefficients, 
// false if not.
bool tough2_file_contains_diffu(tough2_file_t* file);

// Reads any miscellaneous floating point parameters into the array selec, which 
// must be able to store at least 64 reals. The number of parameters read into 
// selec is stored in size.
void tough2_file_read_selec(tough2_file_t* file, real_t* selec, int* size);

// Returns true if the given TOUGH2 file contains miscellaneous floating 
// point parameters, false if not.
bool tough2_file_contains_selec(tough2_file_t* file);

// Reads default values for relative permeability and capillary pressure, 
// which is to be assigned to any domains for which there exist no ROCKS
// entries. rel_perm_params and cap_pres_params should each be able to store
// at least 7 reals.
void tough2_file_read_rpcap(tough2_file_t* file, 
                            tough2_file_rel_perm_t* rel_perm_type,
                            real_t* rel_perm_params,
                            tough2_file_cap_pres_t* cap_pres_type,
                            real_t* cap_pres_params);

// Returns true if the given TOUGH2 file contains default values for 
// relative permeability and capillary pressure, false if not.
bool tough2_file_contains_rpcap(tough2_file_t* file);

// Returns a newly-allocated array containing a set of times at which output 
// is requested, or NULL if there are no such requested times in the file.
// The number of times is stored in num_times.
real_t* tough2_file_read_times(tough2_file_t* file, int* num_times);

// Returns true if the given TOUGH2 file contains times at which output is 
// requested, false if not.
bool tough2_file_contains_times(tough2_file_t* file);

// Returns a newly-allocated array containing sources/sinks in the file, 
// or NULL if there are none. The number of sources is stored in num_sources.
tough2_file_source_t* tough2_file_read_gener(tough2_file_t* file, int* num_sources);

// Returns true if the given TOUGH2 file contains sources/sinks,
// false if not.
bool tough2_file_contains_gener(tough2_file_t* file);

#endif
