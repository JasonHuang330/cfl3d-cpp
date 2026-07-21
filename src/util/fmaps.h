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

namespace fmaps_ns {

void fmaps(FortranArray1DRef<double> nx, FortranArray1DRef<double> ny, FortranArray1DRef<double> nz, FortranArray1DRef<double> area, FortranArray1DRef<double> at, FortranArray2DRef<double> f, FortranArray2DRef<double> qr, FortranArray2DRef<double> ql, int& n, int& nvtq);

} // namespace fmaps_ns
