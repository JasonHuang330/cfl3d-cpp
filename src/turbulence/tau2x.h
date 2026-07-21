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

namespace tau2x_ns {

void tau(int& mgflag, int& nbl, int& jdim, int& kdim, int& idim, FortranArray2DRef<double> q, FortranArray2DRef<double> res, FortranArray2DRef<double> q1, FortranArray2DRef<double> qr, FortranArray2DRef<int> lw, FortranArray1DRef<double> w, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& nblock, FortranArray1DRef<int> igridg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray1DRef<int> iemg);
void tau2x(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, FortranArray4DRef<double> qr, int& js, int& ks, int& is, int& je, int& ke, int& ie, int& kode);

} // namespace tau2x_ns
