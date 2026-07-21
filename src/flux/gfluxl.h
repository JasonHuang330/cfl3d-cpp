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

namespace gfluxl_ns {

void gfluxl(int& i, int& npl, double& xkap, int& idf, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray3DRef<double> sj, FortranArray3DRef<double> dgp, FortranArray3DRef<double> dgm, FortranArray2DRef<double> t, int& nvtq);

} // namespace gfluxl_ns
