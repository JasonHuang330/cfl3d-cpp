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

namespace transmc_ns {

void trans(int& jdim, int& kdim, int& idim, FortranArray2DRef<double> t, FortranArray3DRef<double> tti, FortranArray3DRef<double> ttj, FortranArray3DRef<double> ttk, FortranArray1DRef<double> x, FortranArray1DRef<double> y, FortranArray1DRef<double> z, int& itran, double& rfreqt, double& utran, double& vtran, double& wtran, double& xorg, double& yorg, double& zorg, double& xold, double& yold, double& zold, double& xorg0, double& yorg0, double& zorg0, int& iupdat, double& time2);
void transmc(int& itran, double& rfreqt, double& utran, double& vtran, double& wtran, double& xorg, double& yorg, double& zorg, double& xorg0, double& yorg0, double& zorg0, double& xmc, double& ymc, double& zmc, int& iupdat, double& time2);

} // namespace transmc_ns
