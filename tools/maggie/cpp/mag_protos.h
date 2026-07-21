// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// mag_protos.h — exact C++ prototypes for every translated maggie subroutine.
// Assembled to match the definitions in the cluster .cpp files.
#pragma once

// ---- driver / top-level (no args) ----
void blckdat();
void setup();
void initia();
void hole();
void outer();
void diagnos();
void output();
void chkxyz();
void chkstn(int* iblank);

// ---- grid / cell-center / io ----
void rgrid(int m,double* x,double* y,double* z,int jd,int kd,int ld);
void rp3d(int m,double* x,double* y,double* z,int jd,int kd,int ld,int& ialph,int nmesh);
void cellcen(double* x,double* y,double* z,int jd,int kd,int ld);
void wcelgrd(int m,double* x,double* y,double* z,int* jimage,int* kimage,int* limage,int jd,int kd,int ld);
void getgrd(int m,double* x,double* y,double* z,int* jimage,int* kimage,int* limage,int jd,int kd,int ld);
void getgrd2(int m,double* x,double* y,double* z,int jd,int kd,int ld);
void branch(int m1,int jd,int kd,int ld,int* jimage,int* kimage,int* limage,double* x,double* y,double* z);

// ---- iblank scratch io ----
void putibl(int m,int* iblank,int jd,int kd,int ld);
void putibl2(int m,int* iblank,int jd,int kd,int ld);
void getibl(int m,int* iblank,int jd,int kd,int ld);
void putint(int nwr,int m,int* ji,int* ki,int* li,int* jbpt,int* kbpt,int* lbpt,double* dxint,double* dyint,double* dzint);
// getint sets the global chkst.nwr (COMMON /chkst/), no nwr arg.
void getint(int m,int* ji,int* ki,int* li,int* jbpt,int* kbpt,int* lbpt,double* dxint,double* dyint,double* dzint);
void wiblnk(int* iblank,int jd,int kd,int ld);
void iblcon(int* iblank,int jd,int kd,int ld);
void iblplt(int m,int* iblank,int jd,int kd,int ld);

// ---- hole cutting ----
void hlocat(int* iblank,double* x,double* y,double* z,int jd,int kd,int ld);
void hlocat2(int m,int jd,int kd,int ld,int jmn,int jmx,int kmn,int kmx,int lmn,int lmx,int* iblank);
void bsurf(int js,int je,int ks,int ke,int ls,int le,int jd,int kd,int ld,double* x,double* y,double* z,int m1);
void frnge(int* iblank,int jd,int kd,int ld,int m1);
void bound(int j,int k,int l,double& x1,double& y1,double& z1,int jd,int kd,int ld,double* x,double* y,double* z,int* ife);

// ---- outer boundary ----
void obibl(int m,int n,int* iblank,int jd,int kd,int ld);
void obtot(int m,double* x,double* y,double* z,int* iblank,int jd,int kd,int ld,int i1);

// ---- search / interpolation core ----
void search2(int i,int m,int m1,int jd1,int kd1,int ld1,double* xm1,double* ym1,double* zm1,
             int* jimage,int* kimage,int* limage,double xp,double yp,double zp,int& iok,int icall);
void topol(int m1,int jdimA,int kdimA,int ldimA,double* x,double* y,double* z,
           int* jimage,int* kimage,int* limage,double xp,double yp,double zp,
           double& xie,double& eta,double& zeta,int& jp,int& kp,int& lp,
           int itmax,int limit,int& iok,int js,int je,int ks,int ke,int ls,int le,int idsmin);
void interp(int jmax1,int kmax1,int lmax1,double* x11,double* y11,double* z11,
            int ji,int ki,int li,double xie,double eta,double zeta,
            double& xintrp,double& yintrp,double& zintrp);
void intpt(int itr,int* jimage,int* kimage,int* limage,int jd,int kd,int ld,int m1,
           double* xm1,double* ym1,double* zm1,int jd1,int kd1,int ld1,int m,int i1);
void extrap(int i,int m,int m1,int jd1,int kd1,int ld1,double* xm1,double* ym1,double* zm1,
            double xp,double yp,double zp,int& iok,int jp,int kp,int lp,
            int* jimage,int* kimage,int* limage);
void dsmin(int jmax,int kmax,int lmax,double* x,double* y,double* z,double xp,double yp,double zp,
           int& jp,int& kp,int& lp,int js,int je,int ks,int ke,int ls,int le,
           int jskip,int kskip,int lskip,double& dmin);
void dsmin2(int jmax,int kmax,int lmax,double* x,double* y,double* z,double xp,double yp,double zp,
            int& jp,int& kp,int& lp,int js,int je,int ks,int ke,int ls,int le,double& dmin);
void reserch(int m,int n,int nserch,int i1,int* iskip,int iorph,int iflg);
void orphan(int m,int n,int nserch,int i1,int iorph,int iflg);

// ---- diagnostics / stencil checks ----
void ibchk(int* ib,int* iblank,int ic,int is,int* ji,int* ki,int* li,int m1,int jd,int kd,int ld,
           double* xb,double* yb,double* zb,int* jinew,int* kinew,int* linew,double* x,double* y,double* z);
void cindex(int* ibc,int& ibpnts,int& iipnts,int* jb,int* kb,int* lb,int mesh);
void reorder(int idimA,int ibpnts,int* intpts,int* jb,int* kb,int* lb,int* iord,int jd,int kd,int ld);
void setptr(int m,int m1,int icase);
void pltpts_(int nmesh,int m,int jmax,int kmax,int lmax,double* xx,double* yy,double* zz,int* ib,int* ibb);
void trace_(int icall,int idum1,int idum2,int idum3,int idum4,double dum1,double dum2,double dum3);

// ---- cell-vertex helpers ----
void xe(int jdimA,int kdimA,int ldimA,double* x,double* y,double* z,int jcell,int kcell,int lcell,
        double& xc,double& yc,double& zc,double& xie,double& eta,double& zeta,int& imiss);
void xe2(double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double& xc,
         double y1,double y2,double y3,double y4,double y5,double y6,double y7,double y8,double& yc,
         double z1,double z2,double z3,double z4,double z5,double z6,double z7,double z8,double& zc,
         double& xie,double& eta,double& zeta,int& imiss,int jinc,int kinc,int linc);

// ---- workstation helpers (wkstn.F) ----
int ismin(int n,double* x,int incx);
int ismax(int n,double* x,int incx);

// ---- timing ----
void cputim(int icall);
