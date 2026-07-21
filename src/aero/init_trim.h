// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#pragma once
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <complex>
#include <cstdint>

namespace init_trim_ns {

void init(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, FortranArray3DRef<double> vol, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, int& nummem, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iflagprnt);
void init_trim();

} // namespace init_trim_ns
