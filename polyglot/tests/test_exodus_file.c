// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "polyglot/exodus_file.h"

void test_exodus_file_query(void** state)
{
  size_t real_size;
  float version;
  int num_mpi_processes;
  real_array_t* times = real_array_new();
  bool result = exodus_file_query("testwt-nfaced.exo", &real_size, &version,
                                  &num_mpi_processes, times);
printf("%zd, %f, %d, %zu\n", real_size, version, num_mpi_processes, times->size);
  assert_true(result);
  assert_true(real_size == sizeof(double));
  assert_true(version >= 6.09);
  real_array_free(times);
}

void test_read_exodus_file(void** state)
{
}

void test_write_exodus_file(void** state)
{
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_exodus_file_query),
    unit_test(test_read_exodus_file),
    unit_test(test_write_exodus_file)
  };
  return run_tests(tests);
}
