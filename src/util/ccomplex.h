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

namespace ccomplex_ns {

double ccvmgt(double& a, double& b, bool& c);
double ccmax(double& a1, double& a2);
double ccmax4(double& a1, double& a2, double& a3, double& a4);
double ccmaxcr(double& a1, double& a2);
double ccmaxrc(double& a1, double& a2);
double ccmin(double& a1, double& a2);
double ccmin4(double& a1, double& a2, double& a3, double& a4);
double ccmin8(double& a1, double& a2, double& a3, double& a4, double& a5, double& a6, double& a7, double& a8);
double ccmincr(double& a1, double& a2);
double ccminrc(double& a1, double& a2);
double ccabs(double& a);
double ccsign(double& a, double& b);
double ccsignrc(double& a, double& b);
double cctanh(double& a);
double ccatan2(double& csn, double& ccs);
double ccatan(double& z);
double ccacos(double& z);
double ccasin(double& z);
double cctan(double& z);
double cclog10(double& z);
double ccdim(double& a1, double& a2);
double ccerf(double& a1);

} // namespace ccomplex_ns
