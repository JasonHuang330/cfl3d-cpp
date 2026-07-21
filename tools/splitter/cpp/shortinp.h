// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#pragma once
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <complex>
#include <cstdint>

namespace shortinp_ns {

void shortinp(int& nb1, int& ne1, int& nb2, int& ne2, int& idm, int& jdm, int& kdm, int& iface);

} // namespace shortinp_ns
