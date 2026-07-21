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

namespace amafi_ns {

void amafi(int& kp, int& npl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> ai, FortranArray4DRef<double> bi, FortranArray4DRef<double> ci, FortranArray3DRef<double> dtj, FortranArray2DRef<double> t, int& nvt, FortranArray4DRef<double> dfp, FortranArray4DRef<double> dfm, int& imw);

} // namespace amafi_ns
