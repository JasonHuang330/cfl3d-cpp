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

namespace cblkk_ns {

void cblkk(int& nbli, int& idimr, int& jdimr, int& kdimr, int& idimt, int& jdimt, int& kdimt, FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva, int& it, int& ir, int& iedge, FortranArray3DRef<double> xr, FortranArray3DRef<double> yr, FortranArray3DRef<double> zr, FortranArray3DRef<double> xt, FortranArray3DRef<double> yt, FortranArray3DRef<double> zt, int& ntime, int& lcnt, FortranArray1DRef<double> geom_miss, int& mxbli);

} // namespace cblkk_ns
