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

namespace pre_patch_ns {

void pre_patch(int& nbl, FortranArray2DRef<int> lw, int& icount, int& ninter, FortranArray2DRef<int> iindex, int& intmax, int& nsub1, FortranArray2DRef<int> isav_pat, FortranArray3DRef<int> isav_pat_b, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, FortranArray1DRef<int> iiint1, FortranArray1DRef<int> iiint2, int& maxbl, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, int& ierrflg);

} // namespace pre_patch_ns
