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

namespace resnonin_ns {

void resnonin(int& nbl, int& jdim, int& kdim, int& idim, FortranArray2DRef<double> q, FortranArray1DRef<double> x, FortranArray1DRef<double> y, FortranArray1DRef<double> z, FortranArray2DRef<double> sj, FortranArray2DRef<double> sk, FortranArray2DRef<double> si, FortranArray1DRef<double> vol, FortranArray2DRef<double> res, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim);

} // namespace resnonin_ns
