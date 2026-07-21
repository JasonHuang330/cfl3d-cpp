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

namespace ae_pred_ns {

void ae_pred(FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> stm, FortranArray3DRef<double> stmi, FortranArray2DRef<double> gforcn, FortranArray2DRef<double> gforcnm, FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn, FortranArray2DRef<double> x0, FortranArray3DRef<double> perturb, double& cmyt, double& cnwt, double& xorig, double& yorig, double& zorig, int& nmds, int& maxaes, int& irbtrim, int& maxbl, int& myid);

} // namespace ae_pred_ns
