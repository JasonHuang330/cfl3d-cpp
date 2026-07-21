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

namespace loadgr_ns {

void loadgr(FortranArray1DRef<double> w, int& mgwk, int& lx, int& ly, int& lz, int& jindex, FortranArray2DRef<double> x, FortranArray2DRef<double> y, FortranArray2DRef<double> z, int& mdim, int& ndim, int& idimg, int& jdimg, int& kdimg);

} // namespace loadgr_ns
