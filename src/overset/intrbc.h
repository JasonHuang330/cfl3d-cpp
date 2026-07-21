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

namespace intrbc_ns {

void intrbc(FortranArray4DRef<double> q, int& jdim, int& kdim, int& idim, int& nbl, int& ldim, int& maxbl, int& iitot, FortranArray1DRef<int> lig, FortranArray1DRef<int> iipntsg, FortranArray1DRef<double> dxintg, FortranArray1DRef<double> dyintg, FortranArray1DRef<double> dzintg, FortranArray1DRef<int> iiig, FortranArray1DRef<int> jjig, FortranArray1DRef<int> kkig, FortranArray3DRef<double> qb, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qq, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& icorr, int& iset);

} // namespace intrbc_ns
