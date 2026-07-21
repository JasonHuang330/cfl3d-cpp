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

namespace vlutrp_ns {

void vlutr(int& nvmax, int& n, int& nmax, int& il, int& iu, FortranArray4DRef<double> a, FortranArray4DRef<double> b, FortranArray4DRef<double> c, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim);
void vlutrp(int& nvmax, int& n, int& nmax, int& il, int& iu, FortranArray4DRef<double> a, FortranArray4DRef<double> b, FortranArray4DRef<double> c, FortranArray4DRef<double> g, FortranArray4DRef<double> h);

} // namespace vlutrp_ns
