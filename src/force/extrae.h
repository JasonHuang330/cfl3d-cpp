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

namespace extrae_ns {

void extra(int& jdim, int& kdim, int& nsub, int& l, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& jcell, int& kcell, int& jcl, int& jcr, double& x5, double& y5, double& z5, int& icase, int& ifit);
void extrae(int& jdim, int& kdim, int& nsub, int& l, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& jcell, int& kcell, int& kcl, int& kcr, double& x7, double& y7, double& z7, int& icase, int& ifit);

} // namespace extrae_ns
