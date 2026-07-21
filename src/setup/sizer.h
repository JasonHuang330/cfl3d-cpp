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

namespace sizer_ns {

void add2x(FortranArray4DRef<double> q, FortranArray4DRef<double> qc, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, FortranArray4DRef<double> q1, FortranArray4DRef<double> dq, FortranArray4DRef<double> wq, FortranArray3DRef<double> wqj, FortranArray3DRef<double> wqjk, int& js, int& ks, int& is, int& je, int& ke, int& ie, int& ipass, int& nbl, int& nblc, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& ll, int& myid);
void addx(FortranArray4DRef<double> q, FortranArray4DRef<double> qq, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, FortranArray4DRef<double> q1, FortranArray4DRef<double> dq, FortranArray4DRef<double> wq, FortranArray4DRef<double> wqj, int& nbl, FortranArray3DRef<double> blank, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& ll, int& myid);
void blomax(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> vor, FortranArray3DRef<double> snj0, FortranArray3DRef<double> snk0, FortranArray3DRef<double> sni0, FortranArray3DRef<double> snjm, FortranArray3DRef<double> snkm, FortranArray3DRef<double> snim, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> eoms, int& iprint, int& inmx, FortranArray1DRef<double> eomui, FortranArray1DRef<double> fbl, FortranArray1DRef<double> rhon, FortranArray1DRef<double> amun, FortranArray1DRef<double> vortn, FortranArray1DRef<double> disn, FortranArray1DRef<double> utot, FortranArray1DRef<double> eomuo, FortranArray1DRef<double> damp, int& nblt, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> blank, int& iover, FortranArray3DRef<double> bci, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl);
void int2(int& iq, int& jdim1, int& kdim1, int& idim1, int& jmax2, int& kmax2, int& nsub1, int& jmax1, int& kmax1, int& l1, int& i1, FortranArray4DRef<double> q1, FortranArray3DRef<double> bc2, int& i2, FortranArray4DRef<double> q2, FortranArray1DRef<double> xie2, FortranArray1DRef<double> eta2, FortranArray2DRef<double> q1g, FortranArray2DRef<double> q2g, FortranArray1DRef<double> temp, FortranArray1DRef<int> nblkpt, int& intmax, int& icheck, int& mtype, FortranArray2DRef<int> iindex, int& ifo, int& ldim, int& npt, int& j21, int& j22, int& k21, int& k22, FortranArray4DRef<double> q1wk, double& dthtx, double& dthty, double& dthtz, int& lim_ptch);
void sizer(int& mwork, int& mworki, int& nplots, int& minnode, int& iitot, int& intmax, int& maxxe, int& mxbli, int& nsub1, int& lbcprd, int& lbcemb, int& lbcrad, int& maxbl, int& maxgr, int& maxseg, int& maxcs, int& ncycmax, int& intmx, int& mxxe, int& mptch, int& msub1, int& nmds, int& maxaes, int& maxsegdg, int& nnodes, int& nslave, int& nmaster, int& myhost, int& myid, int& mycomm, int& nplots0, int& maxnode0, int& mxbli0, int& lbcprd0, int& lbcemb0, int& lbcrad0, int& maxbl0, int& maxgr0, int& maxseg0, int& maxcs0, int& ncycmax0, int& intmax0, int& nsub10, int& intmx0, int& mxxe0, int& mptch0, int& msub10, int& ibufdim0, int& nbuf0, int& mxbcfil0, int& nmds0, int& maxaes0, int& maxsegdg0, int& imode, int& ntr, FortranArray1DRef<char[80]> bcfiles, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou);
void cntfa(int& nbl, int& nwfa, int& maxbl, int& maxseg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> ieg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbcidim, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, int& nblock);
void cntblmx(int& nbl, int& jdim, int& kdim, int& idim, int& iover, int& iprint, int& n14, int& maxbl, int& maxseg, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbcidim, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo);

} // namespace sizer_ns
