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

namespace u_doubleprime_ns {

void u_doubleprime(int& idim, int& jdim, int& kdim, FortranArray4DRef<double> q, FortranArray4DRef<double> ux, FortranArray3DRef<double> vol, FortranArray4DRef<double> si, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> vx);

} // namespace u_doubleprime_ns
