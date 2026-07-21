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

namespace chkrotk_d_ns {

void chkrot(int& nbl, int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& nblp, int& jdimp, int& kdimp, int& idimp, FortranArray3DRef<double> xp, FortranArray3DRef<double> yp, FortranArray3DRef<double> zp, int& nface, FortranArray4DRef<double> bcdata, FortranArray3DRef<double> xyzjp, FortranArray3DRef<double> xyzkp, FortranArray3DRef<double> xyzip, int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend, int& mdim, int& ndim, int& lcnt, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, int& maxbl, FortranArray1DRef<double> period_miss, int& lbcprd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);
void chkrotk_d(int& nbl, int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& nblp, int& jdimp, int& kdimp, int& idimp, int& nface, FortranArray4DRef<double> bcdata, FortranArray3DRef<double> xyzkp, int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend, int& mdim, int& ndim, int& lcnt, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, int& maxbl, FortranArray1DRef<double> period_miss, int& lbcprd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);

} // namespace chkrotk_d_ns
