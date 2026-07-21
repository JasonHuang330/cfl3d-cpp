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

namespace bc_patch_ns {

void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, double& cl, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, FortranArray3DRef<int> lwdat, int& myid, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, int& nummem);
void bc_patch(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& maxbl, int& maxgr, int& intmax, int& nsub1, int& maxxe, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, int& ninter, FortranArray2DRef<double> windex, FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt, FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz, FortranArray2DRef<int> isav_pat, FortranArray3DRef<int> isav_pat_b, FortranArray1DRef<int> ireq_ar, FortranArray1DRef<int> ireq_snd, FortranArray1DRef<int> index_ar, FortranArray3DRef<int> keep_trac, FortranArray2DRef<int> keep_trac2, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray2DRef<int> istat2, int& istat_size, int& nummem);

} // namespace bc_patch_ns
