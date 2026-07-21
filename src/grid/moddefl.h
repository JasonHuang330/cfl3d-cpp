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

namespace moddefl_ns {

void moddefl(FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn, FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb, double& cmyt, double& cnwt, double& xorig, double& yorig, double& zorig, int& maxaes, int& nmds, int& irbtrim, int& maxbl, int& myid);
void sqrtcumm(int& ittr, FortranArray1DRef<double> dcl, FortranArray1DRef<double> dcm, FortranArray1DRef<double> ddclda, FortranArray1DRef<double> ddcmda, FortranArray1DRef<double> da, FortranArray1DRef<double> dd, double& dlcln, double& dmtrmn, FortranArray2DRef<double> tp, FortranArray2DRef<double> fp, double& a11, double& a12, double& a22, double& r11, double& r22, double& r33, double& r44, double& alf1, double& alf0, double& ddtrm1, double& ddtrm0);

} // namespace moddefl_ns
