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

namespace addx_ns {

void addx(FortranArray4DRef<double> q, FortranArray4DRef<double> qq, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, FortranArray4DRef<double> q1, FortranArray4DRef<double> dq, FortranArray4DRef<double> wq, FortranArray4DRef<double> wqj, int& nbl, FortranArray3DRef<double> blank, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& ll, int& myid);

} // namespace addx_ns
