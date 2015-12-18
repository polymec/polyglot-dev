// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyglot/cf_file.h"

struct cf_file_t 
{
  bool valid;
  char error[1024];
};

cf_file_t* cf_file_new(const char* filename)
{
  cf_file_t* cf = polymec_malloc(sizeof(cf_file_t));
  return cf;
}

bool cf_file_is_valid(cf_file_t* file)
{
  return file->valid;
}

char* cf_file_error(cf_file_t* file)
{
  return &file->error[0];
}

void cf_file_close(cf_file_t* file)
{
}

