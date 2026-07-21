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

namespace xe2_ns {

void xe(int& jdim, int& kdim, int& nsub, int& l, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> xmid, FortranArray3DRef<double> ymid, FortranArray3DRef<double> zmid, FortranArray3DRef<double> xmide, FortranArray3DRef<double> ymide, FortranArray3DRef<double> zmide, int& jcell, int& kcell, double& xc, double& yc, double& zc, double& xie, double& eta, int& imiss, int& ifit, int& ic0, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);
void xe2(double& x1, double& x2, double& x3, double& x4, double& x5, double& x6, double& x7, double& x8, double& xc, double& y1, double& y2, double& y3, double& y4, double& y5, double& y6, double& y7, double& y8, double& yc, double& xie, double& eta, int& imiss, int& ifit, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);

} // namespace xe2_ns
