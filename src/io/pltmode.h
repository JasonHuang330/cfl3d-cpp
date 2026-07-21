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

namespace pltmode_ns {

void pltmode(int& nm, int& iaes, FortranArray1DRef<int> nblg, int& ngrid, int& maxgr, int& maxbl, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> iaesurf, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, int& maxseg, int& maxsegdg, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray5DRef<double> xmdj, FortranArray5DRef<double> xmdk, FortranArray5DRef<double> xmdi, int& jdim1, int& kdim1, int& idim1, int& maxaes, int& nmds, int& nbl1, int& iunitw, int& iflag);

} // namespace pltmode_ns
