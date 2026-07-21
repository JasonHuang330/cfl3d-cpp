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

namespace ctime1_ns {

void ctime1(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray3DRef<double> vol, FortranArray3DRef<double> sj, FortranArray3DRef<double> sk, FortranArray3DRef<double> si, FortranArray3DRef<double> dtj, FortranArray2DRef<double> t, double& delt, FortranArray3DRef<double> vist3d, int& itur, double& dtmin, int& iout, int& ntime, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& idef);

} // namespace ctime1_ns
