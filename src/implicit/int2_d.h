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

namespace int2_d_ns {

void int2(int& iq, int& jdim1, int& kdim1, int& idim1, int& jmax2, int& kmax2, int& nsub1, int& jmax1, int& kmax1, int& l1, int& i1, FortranArray4DRef<double> q1, FortranArray3DRef<double> bc2, int& i2, FortranArray4DRef<double> q2, FortranArray1DRef<double> xie2, FortranArray1DRef<double> eta2, FortranArray2DRef<double> q1g, FortranArray2DRef<double> q2g, FortranArray1DRef<double> temp, FortranArray1DRef<int> nblkpt, int& intmax, int& icheck, int& mtype, FortranArray2DRef<int> iindex, int& ifo, int& ldim, int& npt, int& j21, int& j22, int& k21, int& k22, FortranArray4DRef<double> q1wk, double& dthtx, double& dthty, double& dthtz, int& lim_ptch);
void int2_d(int& iq, int& jdim1, int& kdim1, int& idim1, int& jmax2, int& kmax2, int& nsub1, int& jmax1, int& kmax1, int& l1, int& i1, FortranArray3DRef<double> q1, FortranArray3DRef<double> bc2, int& i2, FortranArray4DRef<double> q2, FortranArray1DRef<double> xie2, FortranArray1DRef<double> eta2, FortranArray2DRef<double> q1g, FortranArray2DRef<double> q2g, FortranArray1DRef<double> temp, FortranArray1DRef<int> nblkpt, int& intmax, int& icheck, int& mtype, FortranArray2DRef<int> iindex, int& ifo, int& ldim, int& npt, int& j21, int& j22, int& k21, int& k22, FortranArray3DRef<double> q1wk, double& dthtx, double& dthty, double& dthtz, int& lim_ptch);

} // namespace int2_d_ns
