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

namespace getibk0_ns {

void getibk(FortranArray3DRef<double> blank, int& jdim, int& kdim, int& idim, int& nbl, FortranArray1DRef<int> intpts, int& nblpts, int& ibpnts, int& iipnts, int& iitot, int& maxbl, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray1DRef<double> dxintg, FortranArray1DRef<double> dyintg, FortranArray1DRef<double> dzintg, FortranArray1DRef<int> iiig, FortranArray1DRef<int> jjig, FortranArray1DRef<int> kkig, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, int& myid, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou);
void getibk0(int& jdim, int& kdim, int& idim, int& nbl, int& itotb, int& itoti, int& maxbl, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& ibufdim, int& nbuf, int& ierrflg, int& myid);

} // namespace getibk0_ns
