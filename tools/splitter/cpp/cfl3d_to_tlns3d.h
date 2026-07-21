// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#pragma once
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <complex>
#include <cstdint>

namespace cfl3d_to_tlns3d_ns {

void cfl3d_to_tlns3d(int& iver, int& ipatch, int& iptyp, int& ipar, int& nnodes, int& isd, int& mbloc, int& msegn, int& msegt, int& mxbli, FortranArray1DRef<int> il, FortranArray1DRef<int> jl, FortranArray1DRef<int> kl, FortranArray2DRef<double> rkap0g, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<int> idiagg, FortranArray1DRef<int> nblcg, FortranArray2DRef<int> idegg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray2DRef<int> mit, FortranArray1DRef<int> ilamlog, FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamlog, FortranArray1DRef<int> klamhig, FortranArray2DRef<int> iredundant, FortranArray1DRef<int> nseg1, FortranArray1DRef<int> nseg, FortranArray1DRef<int> nsilo, FortranArray1DRef<int> nsihi, FortranArray1DRef<int> nsjlo, FortranArray1DRef<int> nsjhi, FortranArray1DRef<int> nsklo, FortranArray1DRef<int> nskhi, FortranArray3DRef<int> nb1, FortranArray3DRef<int> ne1, FortranArray3DRef<int> nb2, FortranArray3DRef<int> ne2, FortranArray3DRef<int> ibct, FortranArray4DRef<double> bcval, FortranArray1DRef<int> iovrlp, FortranArray3DRef<int> ifoflg, FortranArray3DRef<int> ndat, FortranArray3DRef<double> xmap, FortranArray3DRef<int> imap, FortranArray2DRef<int> ivisb, FortranArray2DRef<double> twotref, FortranArray1DRef<int> itrb1, FortranArray1DRef<int> itrb2, FortranArray1DRef<int> jtrb1, FortranArray1DRef<int> jtrb2, FortranArray1DRef<int> ktrb1, FortranArray1DRef<int> ktrb2, FortranArray1DRef<int> iturbb, FortranArray1DRef<int> ibif1, FortranArray1DRef<int> ibif2, FortranArray1DRef<int> nbi1, FortranArray1DRef<int> nei1, FortranArray1DRef<int> nbj1, FortranArray1DRef<int> nej1, FortranArray1DRef<int> nbk1, FortranArray1DRef<int> nek1, FortranArray1DRef<int> nbi2, FortranArray1DRef<int> nei2, FortranArray1DRef<int> nbj2, FortranArray1DRef<int> nej2, FortranArray1DRef<int> nbk2, FortranArray1DRef<int> nek2, FortranArray1DRef<int> nd11, FortranArray1DRef<int> nd21, FortranArray1DRef<int> nd12, FortranArray1DRef<int> nd22, FortranArray1DRef<int> ifsor, FortranArray1DRef<int> nb1s, FortranArray1DRef<int> ne1s, FortranArray1DRef<int> nb2s, FortranArray1DRef<int> ne2s, int& nrotat, int& ntrans, double& tlref, double& rlref, int& ioflag, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, char (&cflout)[80]);

} // namespace cfl3d_to_tlns3d_ns
