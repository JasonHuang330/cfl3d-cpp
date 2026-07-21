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

namespace rechk_ns {

void rechk(int& mdim1, int& ndim1, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, int& msub1, int& jmax1, int& kmax1, int& l, FortranArray2DRef<double> x1, FortranArray2DRef<double> y1, FortranArray2DRef<double> z1, int& xie1, int& xie2, int& eta1, int& eta2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl);

} // namespace rechk_ns
