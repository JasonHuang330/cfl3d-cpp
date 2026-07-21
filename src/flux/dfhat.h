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

namespace dfhat_ns {

void dfhat(FortranArray1DRef<double> ax, FortranArray1DRef<double> ay, FortranArray1DRef<double> az, FortranArray1DRef<double> area, FortranArray1DRef<double> at, FortranArray2DRef<double> q, FortranArray3DRef<double> df, int& nn, int& nvtq, int& ipm);

} // namespace dfhat_ns
