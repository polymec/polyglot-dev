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

// The TOUGH2 file class provides an interface to mesh and initial condition 
// data within input files written for TOUGH2. This interface is inherently 
// serial because TOUGH2 does not use parallel I/O, so you should create TOUGH2 
// file objects only on a single process. For information on the various blocks 
// within TOUGH2 files please refer to the TOUGH2 Users Guide.

// This type provides the interface for TOUGH2 files.
typedef struct tough2_file_t tough2_file_t;

// Opens a new TOUGH2 file for reading simulation data, returning the TOUGH 
// file object. The type of element name (short or long) must be specified. 
tough2_file_t* tough2_file_new(const char* filename);

// Closes and destroys the given TOUGH2 file.
void tough2_file_close(tough2_file_t* file);

// Reads a mesh from the given TOUGH2 file, returning a newly-allocated 
// point_cloud object with properties and tags corresponding to the 
// elements and connections within. A TOUGH mesh does not contain enough 
// geometric information to be represented as an arbitrary polyhedral mesh, 
// so we store it as a point cloud with an associated neighbor pairing 
// in a property named "connections". 
//
// The following properties can be present in the mesh returned:
//
// volume      - An array containing the volumes of the mesh elements.
// k_mod       - An array containing permeability modifiers (factors) for 
//               each element.
// connections - A neighbor_pairing containing the connections (i, j) between 
//               neighboring elements i and j.
// d1, d2      - Arrays indexed in the same order as the connections, 
//               which assign distances for the connections that identify the 
//               location of a face in non-centroidal geometries. These 
//               properties are only present if they are contained in the 
//               TOUGH2 mesh.
// area        - An array indexed in the same order as the connections,
//               assigning an interfacial area to each connection.
// beta        - An array indexed in the same order as the connections,
//               containing the (signed) cosine of the angle between the 
//               gravity vector g and the line connecting the two elements e1 and 
//               e2. g o beta > 0 implies that e1 is above e2, and g o beta < 0 
//               implies that e2 is above e1.
// sigma       - An array indexed in the same order as the connections,
//               containing the radiant emittance factor for radiative heat 
//               transfer (equal to 1 for a perfect black body and 0 for a 
//               perfect insulator).
// elem_index  - A string->integer map (string_int_unordered_map) that maps 
//               the name of an element in a TOUGH2 mesh to its integer index.
//
// In addition, a tag is created for each reservoir domain identified in 
// the TOUGH2 mesh.
//
// If the file does not contain a mesh or a mesh cannot be generated from it, 
// this function returns NULL.
point_cloud_t* tough2_file_read_mesh(tough2_file_t* file);

// Returns true if the given TOUGH2 file contains a mesh, false if not.
bool tough2_file_contains_mesh(tough2_file_t* file);

// Parses element-specific initial conditions, returning a newly-allocated
// array containing values for primary variables in component-minor order,
// corresponding to the points in the given mesh.
real_t* tough2_file_read_incon(tough2_file_t* file, point_cloud_t* mesh);

// Returns true if the given TOUGH file contains element-specific initial 
// conditions, false if not.
bool tough2_file_contains_incon(tough2_file_t* file);

#endif
