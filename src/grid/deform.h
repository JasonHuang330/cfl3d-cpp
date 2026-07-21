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

namespace deform_ns {

void deform_el(FortranArray3DRef<int> islavept, int& nslave, int& nmaster, FortranArray1DRef<double> u, FortranArray1DRef<double> xst, FortranArray1DRef<double> yst, FortranArray1DRef<double> zst, int& nt, int& myhost, int& mycomm, int& myid, int& nnodes, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> nblelst, int& maxbl, int& iseqr);
void coincdef(FortranArray1DRef<double> volij, FortranArray1DRef<double> u, FortranArray3DRef<int> islavept, int& nslave, int& nmaster, int& iseqr, FortranArray1DRef<int> n1);
void gaussseidel(int& nc, int& n, FortranArray1DRef<double> b, FortranArray1DRef<double> x, FortranArray1DRef<double> sa, FortranArray1DRef<int> ija, double& tol, int& itmax, double& err, FortranArray1DRef<double> xt, int& myid, int& myhost);
void linbcg(int& nc, int& n, FortranArray1DRef<double> b, FortranArray1DRef<double> x, FortranArray1DRef<double> sa, FortranArray1DRef<int> ija, int& itol, double& tol, int& itmax, int& iter, double& err, double& time, double& eps, int& myid);
void atimes(int& n, int& nc, FortranArray1DRef<double> x, FortranArray1DRef<double> r, FortranArray1DRef<double> sa, FortranArray1DRef<int> ija, int& itrnsp, int& myid);
void dsprstx(FortranArray1DRef<double> sa, FortranArray1DRef<int> ija, FortranArray1DRef<double> x, FortranArray1DRef<double> b, int& n, int& nc, int& myid);
void dprec(FortranArray1DRef<double> sa, FortranArray1DRef<int> ija, FortranArray1DRef<double> b, int& n, int& nc);
void dsprsax(FortranArray1DRef<double> sa, FortranArray1DRef<int> ija, FortranArray1DRef<double> x, FortranArray1DRef<double> b, int& n, int& nc, int& myid);
void asolve(int& n, int& nc, FortranArray1DRef<double> b, FortranArray1DRef<double> x, FortranArray1DRef<double> sa, int& itrnsp);
double snrm(int& n, FortranArray1DRef<double> sx, int& itol);
void invdet(FortranArray2DRef<double> UL, FortranArray2DRef<double> COMIN, int& N, double& DTNRM, double& DETM);
void expdecay(FortranArray1DRef<double> xst, FortranArray1DRef<double> yst, FortranArray1DRef<double> zst, FortranArray1DRef<double> u, FortranArray1DRef<double> ut, FortranArray1DRef<double> r, FortranArray2DRef<double> ripm, double& eps, FortranArray3DRef<int> islavept, int& nslave, int& nmaster, int& ndim, int& nt, int& nnodes, int& myhost, int& myid, int& mycomm, FortranArray1DRef<int> mblk2nd, int& maxbl, int& iseqr);
void hookefe(FortranArray1DRef<double> ei, FortranArray1DRef<double> ej, FortranArray1DRef<double> ek, FortranArray1DRef<double> gij, FortranArray1DRef<double> gjk, FortranArray1DRef<double> gik, FortranArray1DRef<double> xst, FortranArray1DRef<double> yst, FortranArray1DRef<double> zst, FortranArray1DRef<double> volij, FortranArray1DRef<double> volik, double& eps, double& eini, double& gini, double& arg1, double& arg2, FortranArray3DRef<int> islavept, int& nslave, int& nmaster, int& nt, int& nnodes, int& myhost, int& myid, int& mycomm, FortranArray1DRef<int> mblk2nd, int& maxbl, int& iseqr);
void elmetricsfe(FortranArray1DRef<double> xix, FortranArray1DRef<double> xiy, FortranArray1DRef<double> xiz, FortranArray1DRef<double> etax, FortranArray1DRef<double> etay, FortranArray1DRef<double> etaz, FortranArray1DRef<double> zetax, FortranArray1DRef<double> zetay, FortranArray1DRef<double> zetaz, FortranArray1DRef<double> ooj, FortranArray1DRef<double> xs, FortranArray1DRef<double> ys, FortranArray1DRef<double> zs, double& eps, FortranArray3DRef<int> islavept, int& nslave, int& nmaster, int& nnodes, int& myhost, int& myid, int& mycomm, FortranArray1DRef<int> mblk2nd, int& maxbl, int& iseqr);
void elrhs(FortranArray1DRef<double> b, FortranArray1DRef<double> u, FortranArray3DRef<int> islavept, int& nslave, int& nmaster, int& iseqr);
void deform(int& nbl, int& idim, int& jdim, int& kdim, FortranArray1DRef<double> x, FortranArray1DRef<double> y, FortranArray1DRef<double> z, FortranArray1DRef<double> xnm2, FortranArray1DRef<double> ynm2, FortranArray1DRef<double> znm2, FortranArray1DRef<double> xnm1, FortranArray1DRef<double> ynm1, FortranArray1DRef<double> znm1, FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, FortranArray4DRef<double> delti, FortranArray1DRef<double> u, FortranArray4DRef<double> vel, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, int& maxbl, double& dt, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& maxseg, FortranArray1DRef<double> wk, int& nsurf, int& irst, int& iflag, FortranArray3DRef<int> islavept, int& nslave, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg, int& nmaster, int& iseq, FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax, FortranArray1DRef<int> kskmax, int& nt);
void deform_surf(int& nbl, int& idim, int& jdim, int& kdim, FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, FortranArray4DRef<double> delti, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, int& maxbl, int& mseq, double& time, double& dt, int& ita, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& maxseg, FortranArray1DRef<double> wk, FortranArray1DRef<double> u, int& nsurf, int& irst, int& iflag, FortranArray3DRef<int> islavept, int& nslave, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg, int& nmaster, int& iseq);

} // namespace deform_ns
