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

namespace ae_corr_ns {

void ae_corr(FortranArray3DRef<double> stm, FortranArray3DRef<double> stmi, FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn, FortranArray2DRef<double> gforcn, FortranArray2DRef<double> gforcs, FortranArray2DRef<double> gforcnm, FortranArray2DRef<double> gf0, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, int& maxbl, int& maxseg, FortranArray2DRef<double> aesrfdat, int& nmds, int& maxaes, int& nt, FortranArray1DRef<int> mblk2nd, int& iseqr, FortranArray1DRef<int> levelg, FortranArray1DRef<int> iadvance, int& nblock, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, int& myid, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, FortranArray3DRef<double> perturb, FortranArray4DRef<double> aehist, int& ncycmax, int& maxsegdg, int& myhost, int& mycomm);

} // namespace ae_corr_ns
