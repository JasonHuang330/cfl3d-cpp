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

namespace qinter_ns {

void add2x(FortranArray4DRef<double> q, FortranArray4DRef<double> qc, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, FortranArray4DRef<double> q1, FortranArray4DRef<double> dq, FortranArray4DRef<double> wq, FortranArray3DRef<double> wqj, FortranArray3DRef<double> wqjk, int& js, int& ks, int& is, int& je, int& ke, int& ie, int& ipass, int& nbl, int& nblc, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& ll, int& myid);
void addx(FortranArray4DRef<double> q, FortranArray4DRef<double> qq, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, FortranArray4DRef<double> q1, FortranArray4DRef<double> dq, FortranArray4DRef<double> wq, FortranArray4DRef<double> wqj, int& nbl, FortranArray3DRef<double> blank, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& ll, int& myid);
void qinter(int& iseq, int& lembed, int& maxl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& maxbl, int& maxgr, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, int& ngrid, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0, FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz, FortranArray1DRef<double> time2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae, int& maxsegdg, int& nummem);

} // namespace qinter_ns
