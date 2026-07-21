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

namespace setslave_ns {

void setslave(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& maxbl, int& mxbli, int& maxgr, int& maxseg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& nblock, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray3DRef<int> islavept, int& nslave, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, int& maxsegdg, FortranArray1DRef<int> iwk, int& iwork, int& nmaster, int& ngrid, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray2DRef<int> iskip, FortranArray1DRef<int> nblg, FortranArray1DRef<int> levelg, int& lfgm, FortranArray2DRef<int> nblk, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray2DRef<int> nblelst, int& nnodes, FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax, FortranArray1DRef<int> kskmax, int& nbli);

} // namespace setslave_ns
