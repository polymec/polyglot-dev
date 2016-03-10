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
  assert_true(cf_file_has_time_series(cf));
  assert_true(cf_file_has_latlon_surface_var(cf, "area"));
  assert_true(cf_file_has_latlon_surface_var(cf, "pr"));
  assert_true(cf_file_has_latlon_surface_var(cf, "tas"));
  assert_true(cf_file_has_latlon_var(cf, "ua"));
  cf_file_close(cf);
}

void test_cf_file_write(void** state)
{
  cf_file_t* cf = cf_file_new("cf_test_write.nc");
  int major, minor, patch;
  cf_file_get_version(cf, &major, &minor, &patch);
  assert_int_equal(1, major);
  assert_int_equal(6, minor);
  assert_int_equal(0, patch);
  assert_false(cf_file_has_latlon_grid(cf));
  assert_false(cf_file_has_time_series(cf));

  int nlat = 100, nlon = 200, nlev = 10;
  real_t lat[nlat];
  for (int i = 0; i < nlat; ++i)
    lat[i] = -90.0 + 180.0*i/(nlat-1);
  real_t lon[nlon];
  for (int i = 0; i < nlon; ++i)
    lon[i] = 360.0*i/(nlon-1);
  real_t lev[nlev];
  for (int i = 0; i < nlev; ++i)
    lev[i] = 20000.0*i/(nlev-1);
  cf_file_define_latlon_grid(cf, 
                             nlat, "degree_north",
                             nlon, "degree_east",
                             nlev, "meter", "up");
  cf_file_write_latlon_grid(cf, lat, lon, lev);
  assert_true(cf_file_has_latlon_grid(cf));

  cf_file_define_time(cf, "days since 0000-1-1", "noleap");
  assert_true(cf_file_has_time_series(cf));

  cf_file_append_time(cf, 0.0);
  cf_file_append_time(cf, 1.0);
  cf_file_append_time(cf, 2.0);

  cf_file_close(cf);

  // Read the file back in and verify its contents.
  cf = cf_file_open("cf_test_write.nc");
  cf_file_get_version(cf, &major, &minor, &patch);
  assert_int_equal(1, major);
  assert_int_equal(6, minor);
  assert_int_equal(0, patch);
  assert_true(cf_file_has_latlon_grid(cf));
  assert_true(cf_file_has_time_series(cf));
  cf_file_close(cf);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  set_log_level(LOG_DEBUG);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_cf_file_open),
    cmocka_unit_test(test_cf_file_write)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
