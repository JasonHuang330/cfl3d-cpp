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

namespace findmin_new_ns {

void findmin_new(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, FortranArray1DRef<int> iwk, int& iwork, int& nsurf, int& j1, FortranArray1DRef<double> xs, FortranArray1DRef<int> ixs, FortranArray1DRef<int> lsminn, FortranArray1DRef<int> ireq_xs, FortranArray1DRef<int> ireq_bb, int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, FortranArray1DRef<int> nblg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg, FortranArray1DRef<int> mblk2nd);
void getpts(int& imn, int& jmn, int& kmn, int& imx, int& jmx, int& kmx, int& imp1, int& jmp1, int& kmp1, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& nface, int& n1beg, int& n1end, int& n2beg, int& n2end, FortranArray2DRef<double> xs, int& nsurf, int& nbl, FortranArray3DRef<double> blank, FortranArray3DRef<double> bci, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck);
void getptsbb(int& imn, int& jmn, int& kmn, int& imx, int& jmx, int& kmx, int& imp1, int& jmp1, int& kmp1, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& nface, int& n1beg, int& n1end, int& n2beg, int& n2end, FortranArray2DRef<double> xs, FortranArray2DRef<int> ixs, int& nsurf, int& nbl, FortranArray3DRef<double> blank, FortranArray3DRef<double> bci, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck);
void bbdist(int& ng, FortranArray2DRef<double> grid, int& nsurf, FortranArray2DRef<double> surf, int& nbb, FortranArray2DRef<double> bbdef, FortranArray2DRef<int> ipv, FortranArray1DRef<int> vlist, FortranArray1DRef<double> dist, FortranArray1DRef<int> idist, int& ncalc, FortranArray1DRef<double> wk3d5, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);
void bbdst1(int& ng, FortranArray2DRef<double> grid, int& nsurf, FortranArray2DRef<double> surf, int& nbb, FortranArray2DRef<double> bbdef, FortranArray2DRef<int> ipv, FortranArray1DRef<int> vlist, FortranArray1DRef<double> dist, FortranArray1DRef<int> idist, int& ncalc, FortranArray1DRef<double> wrk, FortranArray2DRef<double> tsurf, FortranArray1DRef<double> test);
void calc_dist(int& imn, int& jmn, int& kmn, int& imx, int& jmx, int& kmx, int& imp1, int& jmp1, int& kmp1, int& imp2, int& jmp2, int& kmp2, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> smin, int& nsurf, FortranArray2DRef<double> surf, int& nbb, FortranArray2DRef<double> bbdef, FortranArray2DRef<int> ipv, FortranArray1DRef<int> vlist, FortranArray1DRef<int> ntri, FortranArray2DRef<int> iptri, FortranArray1DRef<double> wk3d5, FortranArray1DRef<int> iwrk, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);
void collect_surf(int& imn, int& jmn, int& kmn, int& imxs, int& jmxs, int& kmxs, int& imp1s, int& jmp1s, int& kmp1s, FortranArray2DRef<double> xs, int& nface, int& n1beg, int& n1end, int& n2beg, int& n2end, int& nsurf, FortranArray2DRef<double> surf, int& isurf, FortranArray1DRef<int> ntri, FortranArray2DRef<int> iptri);
void distcc(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> smingp, FortranArray3DRef<double> smincc);
void distcg(int& jdimf, int& kdimf, int& idimf, FortranArray3DRef<double> sminf, int& jdim, int& kdim, int& idim, FortranArray3DRef<double> sminc);
void initf(int& isize);
void initi(int& isize);
void ifree(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& isize);
void ffree(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& isize);
void makebb(int& minpb0, int& maxlv0, int& nwv0, FortranArray1DRef<double> wx, FortranArray1DRef<double> wy, FortranArray1DRef<double> wz, int& nbb0, int& bblen, FortranArray2DRef<double> bbdef, FortranArray2DRef<int> iv, FortranArray1DRef<int> vlist, FortranArray1DRef<int> iwrk, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);
void calcbb(FortranArray1DRef<double> bbmin, FortranArray1DRef<double> bbmax, int& nvi, FortranArray1DRef<int> vi, FortranArray1DRef<double> wx, FortranArray1DRef<double> wy, FortranArray1DRef<double> wz);
void getvrt(FortranArray1DRef<double> bbmin, FortranArray1DRef<double> bbmax, int& nwv, int& nbbv, FortranArray1DRef<int> bbi, FortranArray1DRef<double> wx, FortranArray1DRef<double> wy, FortranArray1DRef<double> wz);
void shells(int& n, FortranArray1DRef<int> is, FortranArray1DRef<double> s);
void spltbb(int& level0, FortranArray2DRef<double> bbdef, FortranArray2DRef<int> iv, FortranArray1DRef<int> vlist, FortranArray1DRef<int> bbi, FortranArray1DRef<double> wx, FortranArray1DRef<double> wy, FortranArray1DRef<double> wz);
void push(FortranArray2DRef<double> st, int& is, FortranArray1DRef<double> bbmin, FortranArray1DRef<double> bbmax, int& lev);
void pop(FortranArray2DRef<double> st, int& is, FortranArray1DRef<double> bbmin, FortranArray1DRef<double> bbmax, int& lev);
void sort_x(int& nsurf, FortranArray2DRef<double> surf, FortranArray1DRef<int> ntri, FortranArray2DRef<int> iptr, FortranArray1DRef<double> wk3d5, FortranArray1DRef<int> iwrk, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);
void move_real(int& n, FortranArray1DRef<double> x, FortranArray1DRef<int> iperm, FortranArray1DRef<double> temp);
void move_integer(int& n, FortranArray1DRef<int> ix, FortranArray1DRef<int> iperm, FortranArray1DRef<int> itemp);
void triang(double& dist, double& pp1, double& pp2, double& pp3, double& aa1, double& aa2, double& aa3, double& aa4, double& bb1, double& bb2, double& bb3, double& bb4, double& cc1, double& cc2, double& cc3, double& cc4);
void heap_sort(int& n, FortranArray1DRef<double> arrin, FortranArray1DRef<int> indx);
void bbarthcg(int& jdim, int& kdim, int& idim, int& jdimc, int& kdimc, int& idimc, FortranArray3DRef<double> xjb, FortranArray3DRef<double> xkb, FortranArray4DRef<double> xib, FortranArray3DRef<double> blnum, FortranArray3DRef<double> xjbc, FortranArray3DRef<double> xkbc, FortranArray4DRef<double> xibc, FortranArray3DRef<double> blnumc);
void calc_distbb(int& imn, int& jmn, int& kmn, int& imx, int& jmx, int& kmx, int& imp1, int& jmp1, int& kmp1, int& imp2, int& jmp2, int& kmp2, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> smin, int& nsurf, FortranArray2DRef<double> surf, int& nbb, FortranArray2DRef<double> bbdef, FortranArray2DRef<int> ipv, FortranArray1DRef<int> vlist, FortranArray1DRef<int> ntri, FortranArray2DRef<int> iptri, FortranArray2DRef<int> ibbarth, FortranArray3DRef<double> xjb, FortranArray3DRef<double> xkb, FortranArray4DRef<double> xib, FortranArray3DRef<double> blnum, FortranArray1DRef<double> wk3d5, FortranArray1DRef<int> iwrk, int& j1, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);
void collect_surfbb(int& imn, int& jmn, int& kmn, int& imxs, int& jmxs, int& kmxs, int& imp1s, int& jmp1s, int& kmp1s, FortranArray2DRef<double> xs, FortranArray2DRef<int> ixs, int& nface, int& n1beg, int& n1end, int& n2beg, int& n2end, int& nsurf, FortranArray2DRef<double> surf, int& isurf, FortranArray1DRef<int> ntri, FortranArray2DRef<int> iptri, FortranArray2DRef<int> ibbarth, int& nbl);
void sort_xbb(int& nsurf, FortranArray2DRef<double> surf, FortranArray1DRef<int> ntri, FortranArray2DRef<int> iptr, FortranArray2DRef<int> ibbarth, FortranArray1DRef<double> wk3d5, FortranArray1DRef<int> iwrk, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);
int iialloc(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& isize);
int ifalloc(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& isize);
void init(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, FortranArray3DRef<double> vol, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, int& nummem, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iflagprnt);

} // namespace findmin_new_ns
