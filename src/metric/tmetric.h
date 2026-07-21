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

namespace tmetric_ns {

void tmetric(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> sj, FortranArray3DRef<double> sk, FortranArray3DRef<double> si, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray2DRef<double> t, FortranArray2DRef<double> t1, FortranArray2DRef<double> t2, FortranArray3DRef<double> t3, int& nbl);

} // namespace tmetric_ns
