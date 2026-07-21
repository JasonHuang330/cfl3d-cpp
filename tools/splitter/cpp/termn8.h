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

namespace termn8_ns {

void termn8(int& myid_stop, int& ierrflg, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou);

} // namespace termn8_ns
