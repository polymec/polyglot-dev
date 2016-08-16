// Copyright (c) 2015-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This file really only checks that the version of polyglot matches that 
// of polymec.

#include "core/polymec_version.h"
#include "polyglot/polyglot.h"

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-variable"
#include "polyglot/polyglot_version.h"

#if POLYGLOT_MAJOR_VERSION > 0

#if POLYMEC_MAJOR_VERSION != POLYGLOT_MAJOR_VERSION
#error "The installed major version of polymec differs from that of polyglot. Please make sure these versions match."
#endif

#endif

