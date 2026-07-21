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

namespace lead_ns {

void lead(int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, int& maxbl);
void pointers(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, int& maxl, int& lembed, int& nstart, int& nwork, int& mwork, int& maxbl, int& maxgr, int& maxseg, FortranArray3DRef<int> lwdat, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, int& nblock, int& myhost, int& myid, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray1DRef<int> nblfine, FortranArray1DRef<int> ilamlog, FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> klamlog, FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamhig, FortranArray2DRef<int> idegg, FortranArray2DRef<int> iwfg, FortranArray2DRef<int> idiagg, FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg, FortranArray2DRef<double> rkap0g, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray1DRef<int> memblock, int& icall, int& nmds, int& maxaes, int& mpihost);

} // namespace lead_ns
