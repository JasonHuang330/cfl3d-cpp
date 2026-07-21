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

namespace aesurf_ns {

void aesurf(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, FortranArray4DRef<double> delti, FortranArray5DRef<double> xmdj, FortranArray5DRef<double> xmdk, FortranArray5DRef<double> xmdi, FortranArray3DRef<double> wkj, FortranArray3DRef<double> wkk, FortranArray3DRef<double> wki, int& maxbl, int& maxseg, int& nmds, int& maxaes, FortranArray2DRef<double> aesrfdat, FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg);

} // namespace aesurf_ns
