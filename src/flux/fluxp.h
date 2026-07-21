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

namespace fluxp_ns {

void fluxp(FortranArray1DRef<double> ax, FortranArray1DRef<double> ay, FortranArray1DRef<double> az, FortranArray1DRef<double> are, FortranArray1DRef<double> at, FortranArray2DRef<double> qe, FortranArray2DRef<double> f, int& n, FortranArray2DRef<double> t, int& jkpro, int& nvtq, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim);

} // namespace fluxp_ns
