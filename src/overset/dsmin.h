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

namespace dsmin_ns {

void dsmin(int& jdim, int& kdim, int& nsub, FortranArray1DRef<int> jjmax, FortranArray1DRef<int> kkmax, int& lmax, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, double& xc, double& yc, double& zc, int& j1, int& k1, int& l1, FortranArray1DRef<int> lout, int& lflag, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2);

} // namespace dsmin_ns
