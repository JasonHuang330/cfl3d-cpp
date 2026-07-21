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

namespace avghole_ns {

void avghole(FortranArray4DRef<double> q, FortranArray3DRef<double> blank, int& jdim, int& kdim, int& idim, int& nbl, int& ldim, int& int_updt);
void xupdt(FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, int& jdim, int& kdim, int& idim, int& nbl, int& ldim, FortranArray3DRef<double> bcj0, FortranArray3DRef<double> bck0, FortranArray3DRef<double> bci0, int& maxbl, int& iitot, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg, FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg, FortranArray3DRef<double> qb, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& int_updt);

} // namespace avghole_ns
