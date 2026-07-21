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

namespace splitter_ns {

void splitter(int& mbloc, int& maxseg, int& npmax, int& mxbli, int& intmax, int& nsub1, int& msegt, int& msegn, int& msplt, int& mtot, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim);

} // namespace splitter_ns
