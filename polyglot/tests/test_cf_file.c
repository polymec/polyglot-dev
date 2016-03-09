// Copyright (c) 2015-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "polyglot/cf_file.h"

void test_cf_file_open(void** state)
{
  cf_file_t* cf = cf_file_open(CMAKE_CURRENT_SOURCE_DIR "/cf_test_data.nc");
  int major, minor, patch;
  cf_file_get_version(cf, &major, &minor, &patch);
  assert_int_equal(1, major);
  assert_int_equal(0, minor);
  assert_int_equal(0, patch);
  assert_true(cf_file_has_latlon_grid(cf));
  assert_true(cf_file_has_latlon_surface_var(cf, "area"));
  assert_true(cf_file_has_latlon_var(cf, "ua"));
  cf_file_close(cf);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  set_log_level(LOG_DEBUG);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_cf_file_open)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
