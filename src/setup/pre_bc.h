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

namespace pre_bc_ns {

void pre_bc(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray2DRef<int> lw_temp, int& maxbl, int& maxgr, int& maxseg, int& ninter, int& intmax, int& nsub1, FortranArray2DRef<int> iindex, FortranArray2DRef<int> isav_pat, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, FortranArray1DRef<int> iiint1, FortranArray1DRef<int> iiint2, int& mxbli, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, FortranArray2DRef<int> nblk, int& lbcprd, FortranArray2DRef<int> isav_prd, FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk, FortranArray1DRef<int> nblg, int& lbcemb, FortranArray1DRef<int> iemg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> isav_emb, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbcidim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray1DRef<int> iadvance, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray3DRef<int> isav_pat_b, FortranArray1DRef<int> levelg, int& nblock, FortranArray2DRef<int> isav_blk);

} // namespace pre_bc_ns
