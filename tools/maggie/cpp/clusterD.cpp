// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// clusterD.cpp — maggie "search / interpolation core" (Cluster D)
// Byte-exact F->C++ port. Build with -ffp-contract=off to match gfortran.
//
// Subroutines translated (11):
//   search2, topol, interp, intpt, extrap, dsmin, dsmin2, reserch,
//   orphan, xe, xe2
//
// ============================================================================
// SIGNATURES  (// SIGNATURE: <exact prototype used in this file>)
// ============================================================================
// SIGNATURE: void search2(int i,int m,int m1,int jd1,int kd1,int ld1,double* xm1,double* ym1,double* zm1,int* jimage,int* kimage,int* limage,double xp,double yp,double zp,int& iok,int icall);
// SIGNATURE: void topol(int m1,int jdimA,int kdimA,int ldimA,double* x,double* y,double* z,int* jimage,int* kimage,int* limage,double xp,double yp,double zp,double& xie,double& eta,double& zeta,int& jp,int& kp,int& lp,int itmax,int limit,int& iok,int js,int je,int ks,int ke,int ls,int le,int idsmin);
// SIGNATURE: void interp(int jmax1,int kmax1,int lmax1,double* x11,double* y11,double* z11,int ji,int ki,int li,double xie,double eta,double zeta,double& xintrp,double& yintrp,double& zintrp);
// SIGNATURE: void intpt(int itr,int* jimage,int* kimage,int* limage,int jd,int kd,int ld,int m1,double* xm1,double* ym1,double* zm1,int jd1,int kd1,int ld1,int m,int i1);
// SIGNATURE: void extrap(int i,int m,int m1,int jd1,int kd1,int ld1,double* xm1,double* ym1,double* zm1,double xp,double yp,double zp,int& iok,int jp,int kp,int lp,int* jimage,int* kimage,int* limage);
// SIGNATURE: void dsmin(int jmax,int kmax,int lmax,double* x,double* y,double* z,double xp,double yp,double zp,int& jp,int& kp,int& lp,int js,int je,int ks,int ke,int ls,int le,int jskip,int kskip,int lskip,double& dmin);
// SIGNATURE: void dsmin2(int jmax,int kmax,int lmax,double* x,double* y,double* z,double xp,double yp,double zp,int& jp,int& kp,int& lp,int js,int je,int ks,int ke,int ls,int le,double& dmin);
// SIGNATURE: void reserch(int m,int n,int nserch,int i1,int* iskip,int iorph,int iflg);
// SIGNATURE: void orphan(int m,int n,int nserch,int i1,int iorph,int iflg);
// SIGNATURE: void xe(int jdimA,int kdimA,int ldimA,double* x,double* y,double* z,int jcell,int kcell,int lcell,double& xc,double& yc,double& zc,double& xie,double& eta,double& zeta,int& imiss);
// SIGNATURE: void xe2(double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double& xc,double y1,double y2,double y3,double y4,double y5,double y6,double y7,double y8,double& yc,double z1,double z2,double z3,double z4,double z5,double z6,double z7,double z8,double& zc,double& xie,double& eta,double& zeta,int& imiss,int jinc,int kinc,int linc);
//
// ============================================================================
// SIGNATURE MISMATCHES vs mag_protos.h  (FLAGGED — please reconcile)
// ============================================================================
//  [MISMATCH-1] xe : mag_protos.h declares
//      void xe(...,int lcell,double& xc,double& yc,double& zc);
//    i.e. it STOPS at zc. The Fortran `xe` has FOUR more trailing OUTPUT args:
//      xie, eta, zeta (double&, the computed cell-local coords) and imiss (int&).
//    These are essential outputs (topol/extrap read xie,eta,zeta,imiss back).
//    Matching the proto literally would make xe non-functional and break its
//    only callers (topol, extrap). I therefore implemented the FULL Fortran
//    signature (16 args). mag_protos.h MUST be extended to:
//      void xe(int,int,int,double*,double*,double*,int,int,int,
//              double&,double&,double&,double&,double&,double&,int&);
//
//  [MISMATCH-2] xe2 : mag_protos.h declares
//      void xe2(x1..x8,double& xc, y1..y8,double& yc, z1..z8,double& zc);
//    i.e. it STOPS at zc. The Fortran `xe2` has SEVEN more trailing args:
//      xie,eta,zeta (double& out), imiss (int& out), jinc,kinc,linc (int in).
//    Implemented with the FULL Fortran signature. mag_protos.h MUST be extended.
//
//  [MISMATCH-3] reserch : mag_protos.h declares `int iskip` (scalar) but the
//    Fortran declares `dimension iskip(idim)` — it is an ARRAY. Implemented as
//    `int* iskip` (1-based). Please change mag_protos.h to `int* iskip`.
//
// ----------------------------------------------------------------------------
// NOTE on `trace`: Fortran `call trace(icall,i1,i2,i3,i4,d1,d2,d3)` maps to the
//   proto `trace_(int,int,int,int,int,double,double,double)`. All trace calls
//   in this cluster follow that (int x5, double x3) shape, so trace_ is used
//   directly with the SAME argument order/values as the Fortran.
// ----------------------------------------------------------------------------

#include "mag_common.h"

// ============================================================================
// xe2 — set up trilinear coefficients + Newton iteration for xie,eta,zeta
// ============================================================================
void xe2(double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double& xc,
         double y1,double y2,double y3,double y4,double y5,double y6,double y7,double y8,double& yc,
         double z1,double z2,double z3,double z4,double z5,double z6,double z7,double z8,double& zc,
         double& xie,double& eta,double& zeta,int& imiss,int jinc,int kinc,int linc)
{
    int    idum1 = 0;
    int    idum2 = 0;
    int    idum3 = 0;
    int    idum4 = 0;
    double dum1  = 0.;
    double dum2  = 0.;
    double dum3  = 0.;

    // initial guess
    xie  = .5;
    eta  = .5;
    zeta = .5;

    int limit = 50;
    int iter  = 0;

    double dx2 = x2 - x1;
    double dy2 = y2 - y1;
    double dz2 = z2 - z1;
    double dx3 = x3 - x1;
    double dy3 = y3 - y1;
    double dz3 = z3 - z1;
    double dx4 = x4 - x1;
    double dy4 = y4 - y1;
    double dz4 = z4 - z1;
    double dx5 = x5 - x1;
    double dy5 = y5 - y1;
    double dz5 = z5 - z1;
    double dx6 = x6 - x1;
    double dy6 = y6 - y1;
    double dz6 = z6 - z1;
    double dx7 = x7 - x1;
    double dy7 = y7 - y1;
    double dz7 = z7 - z1;
    double dx8 = x8 - x1;
    double dy8 = y8 - y1;
    double dz8 = z8 - z1;

    // coefficients for tri-linear representation
    double a2 = dx2;
    double a3 = dx4;
    double a4 = dx5;
    double a5 = dx3 - dx2 - dx4;
    double a6 = dx6 - dx2 - dx5;
    double a7 = dx8 - dx4 - dx5;
    double a8 = dx7 - dx3 - dx6 + dx2 - dx8 + dx4 + dx5;

    double b2 = dy2;
    double b3 = dy4;
    double b4 = dy5;
    double b5 = dy3 - dy2 - dy4;
    double b6 = dy6 - dy2 - dy5;
    double b7 = dy8 - dy4 - dy5;
    double b8 = dy7 - dy3 - dy6 + dy2 - dy8 + dy4 + dy5;

    double c2 = dz2;
    double c3 = dz4;
    double c4 = dz5;
    double c5 = dz3 - dz2 - dz4;
    double c6 = dz6 - dz2 - dz5;
    double c7 = dz8 - dz4 - dz5;
    double c8 = dz7 - dz3 - dz6 + dz2 - dz8 + dz4 + dz5;

    // 2d cases
    if (jinc == 0) {
        a2 = 1.;
        a5 = 0.;
        a6 = 0.;
        a8 = 0.;
        b2 = 1.;
        b5 = 0.;
        b6 = 0.;
        b8 = 0.;
        c2 = 1.;
        c5 = 0.;
        c6 = 0.;
        c8 = 0.;
    }
    if (kinc == 0) {
        a3 = 1.;
        a5 = 0.;
        a7 = 0.;
        a8 = 0.;
        b3 = 1.;
        b5 = 0.;
        b7 = 0.;
        b8 = 0.;
        c3 = 1.;
        c5 = 0.;
        c7 = 0.;
        c8 = 0.;
    }
    if (linc == 0) {
        a4 = 1.;
        a6 = 0.;
        a7 = 0.;
        a8 = 0.;
        b4 = 1.;
        b6 = 0.;
        b7 = 0.;
        b8 = 0.;
        c4 = 1.;
        c6 = 0.;
        c7 = 0.;
        c8 = 0.;
    }

    // newton iteration to invert for xie,eta,zeta
    //
    // convergence criterion for Newton iteration...require 4 orders of
    // reduction to starting error=abs(f1)+abs(f2)+abs(f3)

    double f1 = x1 + xie*( a2 + a5*eta + a6*zeta + a8*eta*zeta )
              + eta*( a3 + a7*zeta ) + a4*zeta - xc;
    double f2 = y1 + xie*( b2 + b5*eta + b6*zeta + b8*eta*zeta )
              + eta*( b3 + b7*zeta ) + b4*zeta - yc;
    double f3 = z1 + xie*( c2 + c5*eta + c6*zeta + c8*eta*zeta )
              + eta*( c3 + c7*zeta ) + c4*zeta - zc;

    double error0 = std::abs(f1)+std::abs(f2)+std::abs(f3);
    double epsf   = 1.0e-4*error0;
    double epsf1  = 1.0e-6;
    if (epsf < epsf1) epsf = epsf1;

    double error = 0.;

L2:
    f1 = x1 + xie*( a2 + a5*eta + a6*zeta + a8*eta*zeta )
       + eta*( a3 + a7*zeta ) + a4*zeta - xc;
    f2 = y1 + xie*( b2 + b5*eta + b6*zeta + b8*eta*zeta )
       + eta*( b3 + b7*zeta ) + b4*zeta - yc;
    f3 = z1 + xie*( c2 + c5*eta + c6*zeta + c8*eta*zeta )
       + eta*( c3 + c7*zeta ) + c4*zeta - zc;

    // d(f1)/d(xie)
    double a2b = a2 + a5*eta + zeta*( a6 + a8*eta );
    // d(f1)/d(eta)
    double a3b = a3 + a5*xie + zeta*( a7 + a8*xie );
    // d(f1)/d(zeta)
    double a4b = a4 + a6*xie + eta*( a7 + a8*xie );
    // d(f2)/d(xie)
    double b2b = b2 + b5*eta + zeta*( b6 + b8*eta );
    // d(f2)/d(eta)
    double b3b = b3 + b5*xie + zeta*( b7 + b8*xie );
    // d(f2)/d(zeta)
    double b4b = b4 + b6*xie + eta*( b7 + b8*xie );
    // d(f3)/d(xie)
    double c2b = c2 + c5*eta + zeta*( c6 + c8*eta );
    // d(f3)/d(eta)
    double c3b = c3 + c5*xie + zeta*( c7 + c8*xie );
    // d(f3)/d(zeta)
    double c4b = c4 + c6*xie + eta*( c7 + c8*xie );

    iter = iter + 1;

    double det = 1./( a2b*(b3b*c4b-b4b*c3b)
                    - a3b*(b2b*c4b-b4b*c2b)
                    + a4b*(b2b*c3b-b3b*c2b) );

    xie = xie - det*( f1*(b3b*c4b-b4b*c3b)
                    - f2*(a3b*c4b-a4b*c3b)
                    + f3*(a3b*b4b-a4b*b3b) );

    eta = eta - det*(-f1*(b2b*c4b-b4b*c2b)
                    + f2*(a2b*c4b-a4b*c2b)
                    - f3*(a2b*b4b-a4b*b2b) );

    zeta = zeta - det*(f1*(b2b*c3b-b3b*c2b)
                     - f2*(a2b*c3b-a3b*c2b)
                     + f3*(a2b*b3b-a3b*b2b) );

    // exit newton iteration if xie,eta,zeta become large...the current target
    // cell is a poor guess. set xie,eta,zeta to a fixed large value so as to
    // trigger a minimum distance search upon exiting this routine
    double huge = 1.e10;
    imiss = 0;
    if (std::abs(xie) > huge) {
        xie = huge;
        imiss = 1;
    }
    if (std::abs(eta) > huge) {
        eta = huge;
        imiss = 1;
    }
    if (std::abs(zeta) > huge) {
        zeta = huge;
        imiss = 1;
    }
    if (imiss > 0) return;

    error = std::abs(f1) + std::abs(f2) + std::abs(f3);
    if (error > epsf && iter <= limit) goto L2;

    // check to ensure point is inside cell
    // (label 999)
    imiss = 0;
    if (xie < -tol.epsc || xie > 1.+tol.epsc ||
        eta < -tol.epsc || eta > 1.+tol.epsc ||
        zeta < -tol.epsc || zeta > 1.+tol.epsc) imiss = 1;

    // check for non-converged newton iteration which nevertheless manages
    // to have  0 < xie,eta,zeta < 1. set xie,eta,zeta to a fixed large
    // value so as to trigger a minimum distance search upon exiting this
    // routine
    if (imiss == 0 && error > epsf) {
        trace_(21,idum1,idum2,idum3,idum4,xie,eta,zeta);
        trace_(22,iter,idum2,idum3,idum4,error,epsf,dum3);
        imiss = 1;
        xie   = huge;
        eta   = huge;
        zeta  = huge;
    }
    (void)dum1; (void)dum2;
    return;
}

// ============================================================================
// xe — build the 8 cell corners then call xe2
// ============================================================================
void xe(int jdimA,int kdimA,int ldimA,double* x,double* y,double* z,int jcell,int kcell,int lcell,
        double& xc,double& yc,double& zc,double& xie,double& eta,double& zeta,int& imiss)
{
    View3<double> X(x,jdimA,kdimA,ldimA);
    View3<double> Y(y,jdimA,kdimA,ldimA);
    View3<double> Z(z,jdimA,kdimA,ldimA);

    // index increments for 2d or 3d case
    int jinc = 1;
    int kinc = 1;
    int linc = 1;
    if (jdimA < 2) jinc = 0;
    if (kdimA < 2) kinc = 0;
    if (ldimA < 2) linc = 0;

    double x1 = X(jcell,kcell,lcell);
    double y1 = Y(jcell,kcell,lcell);
    double z1 = Z(jcell,kcell,lcell);
    double x2 = X(jcell+jinc,kcell,lcell);
    double y2 = Y(jcell+jinc,kcell,lcell);
    double z2 = Z(jcell+jinc,kcell,lcell);
    double x4 = X(jcell,kcell+kinc,lcell);
    double y4 = Y(jcell,kcell+kinc,lcell);
    double z4 = Z(jcell,kcell+kinc,lcell);
    double x3 = X(jcell+jinc,kcell+kinc,lcell);
    double y3 = Y(jcell+jinc,kcell+kinc,lcell);
    double z3 = Z(jcell+jinc,kcell+kinc,lcell);
    double x5 = X(jcell,kcell,lcell+linc);
    double y5 = Y(jcell,kcell,lcell+linc);
    double z5 = Z(jcell,kcell,lcell+linc);
    double x6 = X(jcell+jinc,kcell,lcell+linc);
    double y6 = Y(jcell+jinc,kcell,lcell+linc);
    double z6 = Z(jcell+jinc,kcell,lcell+linc);
    double x7 = X(jcell+jinc,kcell+kinc,lcell+linc);
    double y7 = Y(jcell+jinc,kcell+kinc,lcell+linc);
    double z7 = Z(jcell+jinc,kcell+kinc,lcell+linc);
    double x8 = X(jcell,kcell+kinc,lcell+linc);
    double y8 = Y(jcell,kcell+kinc,lcell+linc);
    double z8 = Z(jcell,kcell+kinc,lcell+linc);

    xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,y1,y2,y3,y4,y5,y6,y7,y8,yc,
        z1,z2,z3,z4,z5,z6,z7,z8,zc,xie,eta,zeta,imiss,
        jinc,kinc,linc);
    return;
}

// ============================================================================
// interp — trilinear interpolation to cell-center coords
// ============================================================================
void interp(int jmax1,int kmax1,int lmax1,double* x11,double* y11,double* z11,
            int ji,int ki,int li,double xie,double eta,double zeta,
            double& xintrp,double& yintrp,double& zintrp)
{
    View3<double> X11(x11,jmax1,kmax1,lmax1);
    View3<double> Y11(y11,jmax1,kmax1,lmax1);
    View3<double> Z11(z11,jmax1,kmax1,lmax1);

    // index increments for 2d or 3d case
    int jinc1 = 1;
    int kinc1 = 1;
    int linc1 = 1;
    if (jmax1 < 2) jinc1 = 0;
    if (kmax1 < 2) kinc1 = 0;
    if (lmax1 < 2) linc1 = 0;

    // cell center value of point in mesh m via trilinear
    // interpolation in target cell ji,ki,li of mesh m1
    double x1 = X11(ji,ki,li);
    double y1 = Y11(ji,ki,li);
    double z1 = Z11(ji,ki,li);
    double x2 = X11(ji+jinc1,ki,li);
    double y2 = Y11(ji+jinc1,ki,li);
    double z2 = Z11(ji+jinc1,ki,li);
    double x4 = X11(ji,ki+kinc1,li);
    double y4 = Y11(ji,ki+kinc1,li);
    double z4 = Z11(ji,ki+kinc1,li);
    double x3 = X11(ji+jinc1,ki+kinc1,li);
    double y3 = Y11(ji+jinc1,ki+kinc1,li);
    double z3 = Z11(ji+jinc1,ki+kinc1,li);
    double x5 = X11(ji,ki,li+linc1);
    double y5 = Y11(ji,ki,li+linc1);
    double z5 = Z11(ji,ki,li+linc1);
    double x6 = X11(ji+jinc1,ki,li+linc1);
    double y6 = Y11(ji+jinc1,ki,li+linc1);
    double z6 = Z11(ji+jinc1,ki,li+linc1);
    double x7 = X11(ji+jinc1,ki+kinc1,li+linc1);
    double y7 = Y11(ji+jinc1,ki+kinc1,li+linc1);
    double z7 = Z11(ji+jinc1,ki+kinc1,li+linc1);
    double x8 = X11(ji,ki+kinc1,li+linc1);
    double y8 = Y11(ji,ki+kinc1,li+linc1);
    double z8 = Z11(ji,ki+kinc1,li+linc1);

    double dx2 = x2 - x1;
    double dy2 = y2 - y1;
    double dz2 = z2 - z1;
    double dx3 = x3 - x1;
    double dy3 = y3 - y1;
    double dz3 = z3 - z1;
    double dx4 = x4 - x1;
    double dy4 = y4 - y1;
    double dz4 = z4 - z1;
    double dx5 = x5 - x1;
    double dy5 = y5 - y1;
    double dz5 = z5 - z1;
    double dx6 = x6 - x1;
    double dy6 = y6 - y1;
    double dz6 = z6 - z1;
    double dx7 = x7 - x1;
    double dy7 = y7 - y1;
    double dz7 = z7 - z1;
    double dx8 = x8 - x1;
    double dy8 = y8 - y1;
    double dz8 = z8 - z1;

    // coefficients for tri-linear interpolation
    double a2 = dx2;
    double a3 = dx4;
    double a4 = dx5;
    double a5 = dx3 - dx2 - dx4;
    double a6 = dx6 - dx2 - dx5;
    double a7 = dx8 - dx4 - dx5;
    double a8 = dx7 - dx3 - dx6 + dx2 - dx8 + dx4 + dx5;

    double b2 = dy2;
    double b3 = dy4;
    double b4 = dy5;
    double b5 = dy3 - dy2 - dy4;
    double b6 = dy6 - dy2 - dy5;
    double b7 = dy8 - dy4 - dy5;
    double b8 = dy7 - dy3 - dy6 + dy2 - dy8 + dy4 + dy5;

    double c2 = dz2;
    double c3 = dz4;
    double c4 = dz5;
    double c5 = dz3 - dz2 - dz4;
    double c6 = dz6 - dz2 - dz5;
    double c7 = dz8 - dz4 - dz5;
    double c8 = dz7 - dz3 - dz6 + dz2 - dz8 + dz4 + dz5;

    xintrp = x1 + xie*( a2 + a5*eta + a6*zeta + a8*eta*zeta )
           + eta*( a3 + a7*zeta ) + a4*zeta;
    yintrp = y1 + xie*( b2 + b5*eta + b6*zeta + b8*eta*zeta )
           + eta*( b3 + b7*zeta ) + b4*zeta;
    zintrp = z1 + xie*( c2 + c5*eta + c6*zeta + c8*eta*zeta )
           + eta*( c3 + c7*zeta ) + c4*zeta;
    return;
}

// ============================================================================
// dsmin — closest grid point (coarse skip) to (xp,yp,zp)
// ============================================================================
void dsmin(int jmax,int kmax,int lmax,double* x,double* y,double* z,double xp,double yp,double zp,
           int& jp,int& kp,int& lp,int js,int je,int ks,int ke,int ls,int le,
           int jskip,int kskip,int lskip,double& dmin)
{
    View3<double> X(x,jmax,kmax,lmax);
    View3<double> Y(y,jmax,kmax,lmax);
    View3<double> Z(z,jmax,kmax,lmax);

    dmin = 1.0e+20;

    jp = 0;
    kp = 0;
    lp = 0;
    int num = (je - js)/jskip + 1;
    if (num > iwrdim) {
        std::fprintf(funit(6),
           "failure in dsmin: num exceeds iwrdim num = %5d iwrdim = %5d\n",num,iwrdim);
        std::fprintf(stderr,"dsmin\n"); std::exit(1);
    }

    // do 1234 k=ks,ke,kskip ; do 1234 l=ls,le,lskip
    for (int k=ks; k<=ke; k+=kskip) {
    for (int l=ls; l<=le; l+=lskip) {

        int jj = 0;
        for (int j=js; j<=je; j+=jskip) {
            jj = jj + 1;
            double dx = xp-X(j,k,l);
            double dy = yp-Y(j,k,l);
            double dz = zp-Z(j,k,l);
            workc.rm[jj] = dx*dx + dy*dy + dz*dz;
        }

        int jmin = ismin(num,workc.rm,1);
        double d1 = workc.rm[jmin];
        if (d1 < dmin) {
            jp = (jmin-1)*jskip + js;
            kp = k;
            lp = l;
            dmin = d1;
        }
    }
    }
    return;
}

// ============================================================================
// dsmin2 — closest point on mesh FACES (for orphan extrapolation start)
// ============================================================================
void dsmin2(int jmax,int kmax,int lmax,double* x,double* y,double* z,double xp,double yp,double zp,
            int& jp,int& kp,int& lp,int js,int je,int ks,int ke,int ls,int le,double& dmin)
{
    View3<double> X(x,jmax,kmax,lmax);
    View3<double> Y(y,jmax,kmax,lmax);
    View3<double> Z(z,jmax,kmax,lmax);

    dmin = 1.0e+20;

    // search only faces of mesh, using original grid density
    int jskip = je - js;
    int kskip = ke - ks;
    int lskip = le - ls;
    if (jskip < 1) jskip = 1;
    if (kskip < 1) kskip = 1;
    if (lskip < 1) lskip = 1;

    int num = je - js + 1;
    if (num > iwrdim) {
        std::fprintf(funit(6),
           "failure in dsmin: num exceeds iwrdim num = %5d iwrdim = %5d\n",num,iwrdim);
        std::fprintf(stderr,"dsmin2\n"); std::exit(1);
    }

    // k-faces
    for (int k=ks; k<=ke; k+=kskip) {
    for (int l=ls; l<=le; l+=1) {

        int jj = 0;
        for (int j=js; j<=je; j+=1) {
            jj = jj + 1;
            double dx = xp-X(j,k,l);
            double dy = yp-Y(j,k,l);
            double dz = zp-Z(j,k,l);
            workc.rm[jj] = dx*dx + dy*dy + dz*dz;
        }

        int jmin = ismin(num,workc.rm,1);
        double d1 = workc.rm[jmin];
        if (d1 < dmin) {
            jp = jmin + js - 1;
            kp = k;
            lp = l;
            dmin = d1;
        }
    }
    }

    // l-faces
    for (int l=ls; l<=le; l+=lskip) {
    for (int k=ks; k<=ke; k+=1) {

        int jj = 0;
        for (int j=js; j<=je; j+=1) {
            jj = jj + 1;
            double dx = xp-X(j,k,l);
            double dy = yp-Y(j,k,l);
            double dz = zp-Z(j,k,l);
            workc.rm[jj] = dx*dx + dy*dy + dz*dz;
        }
        int jmin = ismin(num,workc.rm,1);
        double d1 = workc.rm[jmin];
        if (d1 < dmin) {
            jp = jmin + js - 1;
            kp = k;
            lp = l;
            dmin = d1;
        }
    }
    }

    // j-faces
    num = ke - ks + 1;
    if (num > iwrdim) {
        std::fprintf(funit(6),
           "failure in dsmin: num exceeds iwrdim num = %5d iwrdim = %5d\n",num,iwrdim);
        std::fprintf(stderr,"dsmin2\n"); std::exit(1);
    }
    for (int j=js; j<=je; j+=jskip) {
    for (int l=ls; l<=le; l+=1) {

        int kk = 0;
        for (int k=ks; k<=ke; k+=1) {
            kk = kk + 1;
            double dx = xp-X(j,k,l);
            double dy = yp-Y(j,k,l);
            double dz = zp-Z(j,k,l);
            workc.rm[kk] = dx*dx + dy*dy + dz*dz;
        }
        int kmin = ismin(num,workc.rm,1);
        double d1 = workc.rm[kmin];
        if (d1 < dmin) {
            kp = kmin + ks - 1;
            jp = j;
            lp = l;
            dmin = d1;
        }
    }
    }
    return;
}

// ============================================================================
// topol — search mesh m1 for point xp,yp,zp; return xie,eta,zeta,jp,kp,lp
// ============================================================================
void topol(int m1,int jdimA,int kdimA,int ldimA,double* x,double* y,double* z,
           int* jimage,int* kimage,int* limage,double xp,double yp,double zp,
           double& xie,double& eta,double& zeta,int& jp,int& kp,int& lp,
           int itmax,int limit,int& iok,int js,int je,int ks,int ke,int ls,int le,int idsmin)
{
    const int mfroz = 500;

    View3<int> JIMAGE(jimage,jdimA,kdimA,ldimA);
    View3<int> KIMAGE(kimage,jdimA,kdimA,ldimA);
    View3<int> LIMAGE(limage,jdimA,kdimA,ldimA);

    int jfroz[mfroz+1], kfroz[mfroz+1], lfroz[mfroz+1];

    int ifr = 1;
    int ihu = 1;

    int ifroze = 0;
    int ifroz  = 0;
    int ihuge  = 0;
    iok    = 0;

    int    idum1 = 0;
    int    idum2 = 0;
    int    idum3 = 0;
    int    idum4 = 0;
    double dum1  = 0.;
    double dum2  = 0.;
    double dum3  = 0.;

    // find nearest point (with indicies jp,kp,lp) to xp,yp,zp; start with
    // minimum distance search if idsmin > 0
    // search only on coarser version of mesh
    int jskip = 4;
    int kskip = 4;
    int lskip = 4;
    if (je-js < 4) jskip = 2;
    if (ke-ks < 4) kskip = 2;
    if (le-ls < 4) lskip = 2;

    double dmin = 0.;
    if (idsmin > 0) dsmin(jdimA,kdimA,ldimA,x,y,z,
        xp,yp,zp,jp,kp,lp,js,je,ks,ke,ls,le,jskip,kskip,lskip,dmin);

    // keep within bounds 1 < j,k,l < j/k/ldim-1 for interp/extrap stencils
    jp = std::min( jp , jdimA-1 );
    kp = std::min( kp , kdimA-1 );
    lp = std::min( lp , ldimA-1 );
    jp = std::max( 1 , jp );
    kp = std::max( 1 , kp );
    lp = std::max( 1 , lp );

    int imiss = 0;

    for (int intern=1; intern<=itmax; ++intern) {

        jfroz[intern] = jp;
        kfroz[intern] = kp;
        lfroz[intern] = lp;

        trace_(3,intern,idum2,idum3,idum4,dum1,dum2,dum3);

        // find local xie,eta,zeta via Newton iteration in current target cell jp,kp,lp
        trace_(4,jp,kp,lp,m1,dum1,dum2,dum3);

        xe(jdimA,kdimA,ldimA,x,y,z,jp,kp,lp,xp,yp,zp,
           xie,eta,zeta,imiss);

        trace_(5,idum1,idum2,idum3,idum4,xie,eta,zeta);

        // current target cell correct if imiss = 0
        if (imiss == 0) goto L5556;

        // update current guess for target cell based on result of Newton
        // iteration, with max allowable change set by limit
        int jinc, kinc, linc;
        if (xie >= 0) jinc = (int)std::abs(xie);
        if (xie <  0) jinc = (int)std::abs(xie-1);
        if (eta >= 0) kinc = (int)std::abs(eta);
        if (eta <  0) kinc = (int)std::abs(eta-1);
        if (zeta >= 0) linc = (int)std::abs(zeta);
        if (zeta <  0) linc = (int)std::abs(zeta-1);

        jinc = std::min( jinc , limit );
        kinc = std::min( kinc , limit );
        linc = std::min( linc , limit );

        if (xie > 1.0) {
            jp = jp + jinc;
        } else if (xie < 0.) {
            jp = jp - jinc;
        }
        if (eta > 1.0) {
            kp = kp + kinc;
        } else if (eta < 0.) {
            kp = kp - kinc;
        }
        if (zeta > 1.0) {
            lp = lp + linc;
        } else if (zeta < 0.) {
            lp = lp - linc;
        }

        double xieg  = (double)jp;
        double etag  = (double)kp;
        double zetag = (double)lp;

        // keep within bounds of mesh m1
        jp = std::min( jp , jdimA-1 );
        kp = std::min( kp , kdimA-1 );
        lp = std::min( lp , ldimA-1 );
        jp = std::max( 1 , jp );
        kp = std::max( 1 , kp );
        lp = std::max( 1 , lp );

        // account for any branch cuts
        int jpc = jp;
        int kpc = kp;
        int lpc = lp;

        if (xieg < 1. || xieg > jdimA-1) {
            jpc = JIMAGE(jp,kp,lp);
            kpc = KIMAGE(jp,kp,lp);
            lpc = LIMAGE(jp,kp,lp);
        }
        if (etag < 1. || etag > kdimA-1) {
            jpc = JIMAGE(jp,kp,lp);
            kpc = KIMAGE(jp,kp,lp);
            lpc = LIMAGE(jp,kp,lp);
        }
        if (zetag < 1. || zetag > ldimA-1) {
            jpc = JIMAGE(jp,kp,lp);
            kpc = KIMAGE(jp,kp,lp);
            lpc = LIMAGE(jp,kp,lp);
        }

        if (jpc != jp || kpc != kp || lpc != lp) {
            trace_(23,m1,jpc,kpc,lpc,dum1,dum2,dum3);
            trace_(24,idum1,jp,kp,lp,dum1,dum2,dum3);
        }

        jp = jpc;
        kp = kpc;
        lp = lpc;

        // search routine off track if xie,eta or zeta are huge...try to
        // get back on track with minumum distance search
        double huge = 1.e3;
        if (std::abs(xie) >= huge || std::abs(eta) >= huge ||
            std::abs(zeta) >= huge) {
            ihuge = ihuge + 1;
            if (ihuge > ihu) goto L1000;
            if (ifroz > 0) goto L1000;
            trace_(6,idum1,idum2,idum3,idum4,dum1,dum2,dum3);
            dsmin(jdimA,kdimA,ldimA,x,y,z,
                  xp,yp,zp,jp,kp,lp,js,je,ks,ke,ls,le,jskip,kskip,
                  lskip,dmin);

            // keep within bounds 1 < j,k,l < j/k/ldim-1 for interp/extrap stencils
            jp = std::min( jp , jdimA-1 );
            kp = std::min( kp , kdimA-1 );
            lp = std::min( lp , ldimA-1 );
            jp = std::max( 1 , jp );
            kp = std::max( 1 , kp );
            lp = std::max( 1 , lp );

            continue; // go to 5555
        }

        // check for frozen convergence: search routine keeps returning
        // to the same point, without 0 < xie,eta,zeta < 1 at that point. if frozen,
        // attempt to break out of cycle by using minimum distance search
        ifroz = 0;
        for (int ii=1; ii<=intern; ++ii) {
            int intv = intern-ii+1;
            if (jp == jfroz[intv] && kp == kfroz[intv]
                && lfroz[intv] == lp) {
                ifroz  = 1;
                // jpfroz,kpfroz,lpfroz set in Fortran but unused
            }
        }
        if (ifroz == 1) {
            ifroze = ifroze + 1;
            if (ifroze > ifr) goto L1000;
            if (ihuge > 0) goto L1000;
            trace_(7,jp,kp,lp,idum4,dum1,dum2,dum3);
            dsmin(jdimA,kdimA,ldimA,x,y,z,
                  xp,yp,zp,jp,kp,lp,js,je,ks,ke,ls,le,jskip,kskip,
                  lskip,dmin);

            // keep within bounds 1 < j,k,l < j/k/ldim-1 for interp/extrap stencils
            jp = std::min( jp , jdimA-1 );
            kp = std::min( kp , kdimA-1 );
            lp = std::min( lp , ldimA-1 );
            jp = std::max( 1 , jp );
            kp = std::max( 1 , kp );
            lp = std::max( 1 , lp );
        }

    } // 5555 continue (do loop)

    // search routine has been unsuccessful
L1000:
    iok = 0;
    (void)dum1; (void)dum2; (void)dum3;
    return;

    // search routine has been successful
L5556:
    iok = 1;
    return;
}

// ============================================================================
// search2 — top-level search over mesh m1 for the fringe point xp,yp,zp
// ============================================================================
void search2(int i,int m,int m1,int jd1,int kd1,int ld1,double* xm1,double* ym1,double* zm1,
             int* jimage,int* kimage,int* limage,double xp,double yp,double zp,int& iok,int icall)
{
    (void)m;
    int    idum1 = 0;
    int    idum2 = 0;
    int    idum3 = 0;
    int    idum4 = 0;
    double dum1  = 0.;
    double dum2  = 0.;
    double dum3  = 0.;

    int itmax  = 50;
    int limit  = 5;
    int icount = 1;
    (void)icount;

    // search limits for minimum distance search
    int js = grdlim.jsrs[m1];
    int je = grdlim.jsre[m1];
    int ks = grdlim.ksrs[m1];
    int ke = grdlim.ksre[m1];
    int ls = grdlim.lsrs[m1];
    int le = grdlim.lsre[m1];

    // 17085 continue
    trace_(1,i,idum2,idum3,idum4,xp,yp,zp);

    // use last successful solution point in mesh m1 as the initial
    // guess for the cell which surrounds point i if one exists; otherwise,
    // use minimum distance search to get a starting point
    double xiet = 0., etat = 0., zetat = 0.;
    int jp = 0, kp = 0, lp = 0;
    int idsmin;

    if (icall == 1) {
        idsmin = 1;
    } else {
        idsmin = 1; // default if the do-loop finds no match (label 555 fallthrough)
        for (int ii=1; ii<=i-1; ++ii) {
            int iii = i - ii;
            if (intrp2.jn[iii] > 0 &&
                std::abs(wherec.nblkpt[iii])-mdim == m1) {
                xiet  = intrp1.xi[iii];
                etat  = intrp1.yi[iii];
                zetat = intrp1.zi[iii];
                jp    = intrp2.jn[iii];
                kp    = intrp2.kn[iii];
                lp    = intrp2.ln[iii];
                idsmin = 0;
                goto L556;
            }
        }
        // 555 continue then idsmin = 1 (already set above)
        L556: ;
    }

    topol(m1,jd1,kd1,ld1,xm1,ym1,zm1,jimage,kimage,limage,
          xp,yp,zp,xiet,etat,zetat,jp,kp,lp,itmax,
          limit,iok,js,je,ks,ke,ls,le,idsmin);

    if (iok == 1) {
        // search successful
        intrp2.jn[i]  = jp;
        intrp2.kn[i]  = kp;
        intrp2.ln[i]  = lp;
        intrp1.xi[i]  = xiet;
        intrp1.yi[i]  = etat;
        intrp1.zi[i]  = zetat;
    } else {
        trace_(20,idum1,idum2,idum3,idum4,dum1,dum2,dum3);
    }
    return;
}

// ============================================================================
// intpt — find interpolation coeffs for fringe/boundary points of mesh m
// ============================================================================
void intpt(int itr,int* jimage,int* kimage,int* limage,int jd,int kd,int ld,int m1,
           double* xm1,double* ym1,double* zm1,int jd1,int kd1,int ld1,int m,int i1)
{
    (void)jd; (void)kd; (void)ld;
    int    idum1 = 0;
    int    idum2 = 0;
    int    idum3 = 0;
    int    idum4 = 0;
    double dum1  = 0.;
    double dum2  = 0.;
    double dum3  = 0.;
    (void)idum1;

    if (intrp2.itotal > idim) {
        // 601 format('0',' dimension exceeded on interpolation arrays for'
        //   ,1x,'mesh =',i5,2x,'idim =',i5,2x,'itotal = ',i5)
        std::fprintf(funit(6),
            "0 dimension exceeded on interpolation arrays for mesh =%5d  idim =%5d  itotal = %5d\n",
            m, idim, intrp2.itotal);
        std::fprintf(stderr,"intpt\n"); std::exit(1);
    }

    trace_(itr,m1,idum2,idum3,idum4,dum1,dum2,dum3);

    int icall = 0;

    for (int i=i1; i<=intrp2.itotal; ++i) {

        if (wherec.nblkpt[i] == 0) {

            double xp = surf2.xbo[i];
            double yp = surf2.ybo[i];
            double zp = surf2.zbo[i];

            // search over mesh m1 to find the cell surrounding the mesh m
            // point xp,yp,zp and compute the required interpolation data
            icall = icall + 1;
            int iok;
            search2(i,m,m1,jd1,kd1,ld1,xm1,ym1,zm1,jimage,kimage,
                    limage,xp,yp,zp,iok,icall);
            if (iok == 1) wherec.nblkpt[i] = -(m1+mdim);
        }
    }
    return;
}

// ============================================================================
// extrap — search mesh m1 for cell nearest xp,yp,zp; extrapolation data
// ============================================================================
void extrap(int i,int m,int m1,int jd1,int kd1,int ld1,double* xm1,double* ym1,double* zm1,
            double xp,double yp,double zp,int& iok,int jp,int kp,int lp,
            int* jimage,int* kimage,int* limage)
{
    (void)m;
    const int mfroz = 999;

    View3<int> JIMAGE(jimage,jd1,kd1,ld1);
    View3<int> KIMAGE(kimage,jd1,kd1,ld1);
    View3<int> LIMAGE(limage,jd1,kd1,ld1);

    int    jfroz[mfroz+1], kfroz[mfroz+1], lfroz[mfroz+1];
    double xfroz[mfroz+1], yfroz[mfroz+1], zfroz[mfroz+1];

    int    idum1 = 0;
    int    idum2 = 0;
    int    idum3 = 0;
    int    idum4 = 0;
    double dum1  = 0.;
    double dum2  = 0.;
    double dum3  = 0.;

    int    limit  = 1;
    int    itmax  = 200;
    iok    = 0;
    int    ijump  = 0;
    double xie    = .5;
    double eta    = .5;
    double zeta   = .5;

    trace_(11,i,idum2,idum3,idum4,xp,yp,zp);

    int imiss = 0;

    for (int intern=1; intern<=itmax; ++intern) {

        jfroz[intern]  = jp;
        kfroz[intern]  = kp;
        lfroz[intern]  = lp;

        trace_(3,intern,idum2,idum3,idum4,dum1,dum2,dum3);

        // find local xie,eta,zeta via Newton iteration in current target cell
        trace_(4,jp,kp,lp,m1,dum1,dum2,dum3);

        xe(jd1,kd1,ld1,xm1,ym1,zm1,jp,kp,lp,xp,yp,zp,
           xie,eta,zeta,imiss);

        xfroz[intern] = xie;
        yfroz[intern] = eta;
        zfroz[intern] = zeta;

        // current target cell correct if imiss = 0 (interpolation, not extrapolation)
        if (imiss == 0) {
            trace_(5,idum1,idum2,idum3,idum4,xie,eta,zeta);
            goto L5556;
        }

        // current block m1 a poor guess if xie,eta or zeta become large. exit
        // to try another mesh m1
        {
            double huge = 1.e5;
            if (std::abs(xie) >= huge || std::abs(eta) >= huge ||
                std::abs(zeta) >= huge) {
                iok = 0;
                return;
            }
        }

        trace_(5,idum1,idum2,idum3,idum4,xie,eta,zeta);

        // update current guess for target cell based on result of Newton
        // iteration, with max allowable change set by limit
        int jsav = jp;
        int ksav = kp;
        int lsav = lp;
        (void)jsav; (void)ksav; (void)lsav;

        int jinc, kinc, linc;
        if (xie >= 0) jinc = (int)std::abs(xie);
        if (xie <  0) jinc = (int)std::abs(xie-1);
        if (eta >= 0) kinc = (int)std::abs(eta);
        if (eta <  0) kinc = (int)std::abs(eta-1);
        if (zeta >= 0) linc = (int)std::abs(zeta);
        if (zeta <  0) linc = (int)std::abs(zeta-1);

        jinc = std::min( jinc , limit );
        kinc = std::min( kinc , limit );
        linc = std::min( linc , limit );

        if (xie > 1.0) {
            jp = jp + jinc;
        } else if (xie < 0.) {
            jp = jp - jinc;
        }
        if (eta > 1.0) {
            kp = kp + kinc;
        } else if (eta < 0.) {
            kp = kp - kinc;
        }
        if (zeta > 1.0) {
            lp = lp + linc;
        } else if (zeta < 0.) {
            lp = lp - linc;
        }

        double xieg  = (double)jp;
        double etag  = (double)kp;
        double zetag = (double)lp;

        // keep within bounds of mesh m1
        jp = std::min( jp , jd1-1 );
        kp = std::min( kp , kd1-1 );
        lp = std::min( lp , ld1-1 );
        jp = std::max( 1 , jp );
        kp = std::max( 1 , kp );
        lp = std::max( 1 , lp );

        // account for any branch cuts
        int jpc = jp;
        int kpc = kp;
        int lpc = lp;

        if (xieg < 1. || xieg > jd1-1) {
            jpc = JIMAGE(jp,kp,lp);
            kpc = KIMAGE(jp,kp,lp);
            lpc = LIMAGE(jp,kp,lp);
        }
        if (etag < 1. || etag > kd1-1) {
            jpc = JIMAGE(jp,kp,lp);
            kpc = KIMAGE(jp,kp,lp);
            lpc = LIMAGE(jp,kp,lp);
        }
        if (zetag < 1. || zetag > ld1-1) {
            jpc = JIMAGE(jp,kp,lp);
            kpc = KIMAGE(jp,kp,lp);
            lpc = LIMAGE(jp,kp,lp);
        }

        if (jpc != jp || kpc != kp || lpc != lp) {
            ijump = ijump + 1;
            // allow only two jumps across branch cut
            if (ijump <= 2) {
                trace_(23,m1,jpc,kpc,lpc,dum1,dum2,dum3);
                trace_(24,idum1,jp,kp,lp,dum1,dum2,dum3);
                jp = jpc;
                kp = kpc;
                continue; // go to 5555
            }
        }

        // check for frozen convergence: search routine keeps returning to the
        // same point, without 0 < xie,eta,zeta < 1. of those cells searched in
        // the frozen cycle, use the point with the minimum extrapolation
        // coefficient (unless it extrapolates in two or three directions).
        // note: extrapolation measured from 0.0 for neg, 1.0 for pos.
        int ifroz = 0;
        for (int ii=1; ii<=intern; ++ii) {
            int intv = intern-ii+1;
            if (jp == jfroz[intv] && kp == kfroz[intv]
                && lfroz[intv] == lp) ifroz = 1;
        }
        if (ifroz > 0) {
            // frozen convergence...search is cyclic
            double xi0 = xfroz[1];
            double yi0 = yfroz[1];
            double zi0 = zfroz[1];
            int iext = 1;
            int nex0 = 0;

            double ximod0, yimod0, zimod0;
            if (xi0-1.0 > tol.epsc) {
                nex0 = nex0 + 1;
                ximod0 = std::abs(xi0-1.0);
            } else if (xi0 < -tol.epsc) {
                nex0 = nex0 + 1;
                ximod0 = std::abs(xi0);
            } else {
                ximod0 = 0.;
            }
            if (yi0-1.0 > tol.epsc) {
                nex0 = nex0 + 1;
                yimod0 = std::abs(yi0-1.0);
            } else if (yi0 < -tol.epsc) {
                nex0 = nex0 + 1;
                yimod0 = std::abs(yi0);
            } else {
                yimod0 = 0.;
            }
            if (zi0-1.0 > tol.epsc) {
                nex0 = nex0 + 1;
                zimod0 = std::abs(zi0-1.0);
            } else if (zi0 < -tol.epsc) {
                nex0 = nex0 + 1;
                zimod0 = std::abs(zi0);
            } else {
                zimod0 = 0.;
            }
            double exmax0 = std::max(std::max(ximod0,yimod0),zimod0);

            for (int ii=2; ii<=intern; ++ii) {
                double xi1 = xfroz[ii];
                double yi1 = yfroz[ii];
                double zi1 = zfroz[ii];
                int nex1 = 0;
                double ximod1, yimod1, zimod1;
                if (xi1-1.0 > tol.epsc) {
                    nex1 = nex1 + 1;
                    ximod1 = std::abs(xi1-1.0);
                } else if (xi1 < -tol.epsc) {
                    nex1 = nex1 + 1;
                    ximod1 = std::abs(xi1);
                } else {
                    ximod1 = 0.;
                }
                if (yi1-1.0 > tol.epsc) {
                    nex1 = nex1 + 1;
                    yimod1 = std::abs(yi1-1.0);
                } else if (yi1 < -tol.epsc) {
                    nex1 = nex1 + 1;
                    yimod1 = std::abs(yi1);
                } else {
                    yimod1 = 0.;
                }
                if (zi1-1.0 > tol.epsc) {
                    nex1 = nex1 + 1;
                    zimod1 = std::abs(zi1-1.0);
                } else if (zi1 < -tol.epsc) {
                    nex1 = nex1 + 1;
                    zimod1 = std::abs(zi1);
                } else {
                    zimod1 = 0.;
                }
                double exmax1 = std::max(std::max(ximod1,yimod1),zimod1);

                // prefer the stencil with extrapolation in fewest directions
                if (nex1 < nex0) {
                    nex0 = nex1;
                    ximod0 = ximod1;
                    yimod0 = yimod1;
                    zimod0 = zimod1;
                    exmax0 = exmax1;
                    iext = ii;
                // prefer the stencil with smaller extrapolation coefficient,
                // provided it does not extrapolate in more directions
                } else if (exmax1 < exmax0 && nex1 == nex0) {
                    nex0 = nex1;
                    ximod0 = ximod1;
                    yimod0 = yimod1;
                    zimod0 = zimod1;
                    exmax0 = exmax1;
                    iext = ii;
                }
            }

            jp   = jfroz[iext];
            kp   = kfroz[iext];
            lp   = lfroz[iext];
            xie  = xfroz[iext];
            eta  = yfroz[iext];
            zeta = zfroz[iext];
            trace_(17,idum1,idum2,idum3,idum4,dum1,dum2,dum3);
            trace_(4,jp,kp,lp,m1,dum1,dum2,dum3);
            trace_(5,idum1,idum2,idum3,idum4,xie,eta,zeta);
            goto L5556;
        }

    } // 5555 continue

L5556:
    intrp2.jn[i]  = jp;
    intrp2.kn[i]  = kp;
    intrp2.ln[i]  = lp;
    intrp1.xi[i]  = xie;
    intrp1.yi[i]  = eta;
    intrp1.zi[i]  = zeta;
    iok = 1;
    return;
}

// ============================================================================
// reserch — re-check initial-pass stencils for hole/fringe points
// ============================================================================
void reserch(int m,int n,int nserch,int i1,int* iskip,int iorph,int iflg)
{
    // COMMON /book3/ iblank(mlen) -> use book3.iblank (1-based flat)
    int* iblank = book3.iblank;
    double* x = grid1.x;
    double* y = grid1.y;
    double* z = grid1.z;
    int* jimage = image.jimage;
    int* kimage = image.kimage;
    int* limage = image.limage;

    // local automatic arrays ds(idim),nbltmp(idim)  — 1-based
    static double ds[idim+1];
    static int    nbltmp[idim+1];

    int notok = 0;
    for (int i=i1; i<=intrp2.itotal; ++i) {
        iskip[i]  = 9999;
        nbltmp[i] = 0;
    }

    for (int nn=1; nn<=nserch; ++nn) {
        int m1;
        if (iflg == 0) {
            m1 = conec1.LHOLE(m,n,nn);
        } else {
            int nnn = 1;
            m1 = conec1.LOUTR(m,nnn,nn);
        }
        int jd1 = grid2.mjmax[m1];
        int kd1 = grid2.mkmax[m1];
        int ld1 = grid2.mlmax[m1];
        getibl(m1,iblank,jd1,kd1,ld1);
        for (int i=i1; i<=intrp2.itotal; ++i) {
            if (wherec.nblkpt[i] != -(m1+mdim)) continue; // go to 1115

            int j   = intrp2.jn[i];
            int k   = intrp2.kn[i];
            int l   = intrp2.ln[i];
            int jp1 = std::min( j+1,jd1 );
            int kp1 = std::min( k+1,kd1 );
            int lp1 = std::min( l+1,ld1 );
            int ii1 = j   + (k-1)*jd1   + (l-1)*jd1*kd1;
            int ii2 = jp1 + (k-1)*jd1   + (l-1)*jd1*kd1;
            int ii3 = jp1 + (kp1-1)*jd1 + (l-1)*jd1*kd1;
            int ii4 = j   + (kp1-1)*jd1 + (l-1)*jd1*kd1;
            int ii5 = j   + (k-1)*jd1   + (lp1-1)*jd1*kd1;
            int ii6 = jp1 + (k-1)*jd1   + (lp1-1)*jd1*kd1;
            int ii7 = jp1 + (kp1-1)*jd1 + (lp1-1)*jd1*kd1;
            int ii8 = j   + (kp1-1)*jd1 + (lp1-1)*jd1*kd1;

            if (iblank[ii1] <= 0 ||
                iblank[ii2] <= 0 ||
                iblank[ii3] <= 0 ||
                iblank[ii4] <= 0 ||
                iblank[ii5] <= 0 ||
                iblank[ii6] <= 0 ||
                iblank[ii7] <= 0 ||
                iblank[ii8] <= 0) {
                // stencil contains hole/fringe point from mesh m1: do not
                // search mesh m1 in next attempt
                iskip[i] = m1;
                notok = notok + 1;
            }
        }
    }

    if (notok > 0 && nserch > 1) {
        std::fprintf(funit(6),"\n");
        std::fprintf(funit(6),
            "     initial pass through search routine yielded %d illegal stencils\n", notok);
        std::fprintf(funit(6),
            "     ...will attempt to extrapolate/interpolate from another mesh\n");
        for (int i=i1; i<=intrp2.itotal; ++i)
            ds[i] = 1.e30;

        // try again for those points which had illegal stencils
        int nok = 0;
        for (int nn=1; nn<=nserch; ++nn) {
            int m1;
            if (iflg == 0) {
                m1 = conec1.LHOLE(m,n,nn);
            } else {
                int nnn = 1;
                m1 = conec1.LOUTR(m,nnn,nn);
            }
            int jd1 = grid2.mjmax[m1];
            int kd1 = grid2.mkmax[m1];
            int ld1 = grid2.mlmax[m1];
            getibl(m1,iblank,jd1,kd1,ld1);
            getgrd(m1,x,y,z,jimage,kimage,limage,jd1,kd1,ld1);
            for (int i=i1; i<=intrp2.itotal; ++i) {
                if (iskip[i] == 9999) continue; // go to 234
                if (iskip[i] == m1)   continue; // go to 234
                int j  = intrp2.jb[i];
                int k  = intrp2.kb[i];
                int l  = intrp2.lb[i];
                (void)j; (void)k; (void)l;
                double xp = surf2.xbo[i];
                double yp = surf2.ybo[i];
                double zp = surf2.zbo[i];
                int jnsav = intrp2.jn[i];
                int knsav = intrp2.kn[i];
                int lnsav = intrp2.ln[i];
                double xisav = intrp1.xi[i];
                double yisav = intrp1.yi[i];
                double zisav = intrp1.zi[i];

                // minimum distance search
                int js = grdlim.jsrs[m1];
                int je = grdlim.jsre[m1];
                int ks = grdlim.ksrs[m1];
                int ke = grdlim.ksre[m1];
                int ls = grdlim.lsrs[m1];
                int le = grdlim.lsre[m1];
                int jp, kp, lp;
                double dmin1;
                dsmin2(jd1,kd1,ld1,x,y,z,
                       xp,yp,zp,jp,kp,lp,js,je,ks,ke,ls,le,dmin1);

                if (dmin1 < ds[i]) {
                    ds[i] = dmin1;

                    int jpold = jp;
                    int kpold = kp;
                    int lpold = lp;

                    int iok;
                    if (iorph == 0) {
                        // use nearest neighbor (zeroth order interpolation)
                        intrp2.jn[i] = jp;
                        intrp2.kn[i] = kp;
                        intrp2.ln[i] = lp;
                        intrp1.xi[i] = 0.;
                        intrp1.yi[i] = 0.;
                        intrp1.zi[i] = 0.;
                        iok = 1;
                    } else {
                        // keep within bounds 1 < j,k,l < j/k/ld1-1
                        jp = std::min( jp , jd1-1 );
                        kp = std::min( kp , kd1-1 );
                        lp = std::min( lp , ld1-1 );
                        jp = std::max( 1 , jp );
                        kp = std::max( 1 , kp );
                        lp = std::max( 1 , lp );
                        extrap(i,m,m1,jd1,kd1,ld1,x,y,z,xp,yp,zp,iok,
                               jp,kp,lp,jimage,kimage,limage);
                        if (iok == 0) {
                            // linear extrapolation not successful...use zeroth order
                            intrp2.jn[i] = jpold;
                            intrp2.kn[i] = kpold;
                            intrp2.ln[i] = lpold;
                            intrp1.xi[i] = 0.;
                            intrp1.yi[i] = 0.;
                            intrp1.zi[i] = 0.;
                        }
                        // for iorph=2, check extrapolation coefficients; if
                        // outside -0.5..1.5, use nearest neighbor instead
                        if (iorph == 2) {
                            if ((intrp1.xi[i] < -0.5-tol.epsc || intrp1.xi[i] > 1.5+tol.epsc) ||
                                (intrp1.yi[i] < -0.5-tol.epsc || intrp1.yi[i] > 1.5+tol.epsc) ||
                                (intrp1.zi[i] < -0.5-tol.epsc || intrp1.zi[i] > 1.5+tol.epsc)) {
                                if (iok > 0) {
                                    intrp2.jn[i] = jpold;
                                    intrp2.kn[i] = kpold;
                                    intrp2.ln[i] = lpold;
                                    intrp1.xi[i] = 0.;
                                    intrp1.yi[i] = 0.;
                                    intrp1.zi[i] = 0.;
                                }
                            }
                        }
                        iok = 1;
                    }
                    if (iok > 0) {
                        int jj  = intrp2.jn[i];
                        int kk  = intrp2.kn[i];
                        int ll  = intrp2.ln[i];
                        int jp1 = std::min( jj+1,jd1 );
                        int kp1 = std::min( kk+1,kd1 );
                        int lp1 = std::min( ll+1,ld1 );
                        int ii1 = jj  + (kk-1)*jd1  + (ll-1)*jd1*kd1;
                        int ii2 = jp1 + (kk-1)*jd1  + (ll-1)*jd1*kd1;
                        int ii3 = jp1 + (kp1-1)*jd1 + (ll-1)*jd1*kd1;
                        int ii4 = jj  + (kp1-1)*jd1 + (ll-1)*jd1*kd1;
                        int ii5 = jj  + (kk-1)*jd1  + (lp1-1)*jd1*kd1;
                        int ii6 = jp1 + (kk-1)*jd1  + (lp1-1)*jd1*kd1;
                        int ii7 = jp1 + (kp1-1)*jd1 + (lp1-1)*jd1*kd1;
                        int ii8 = jj  + (kp1-1)*jd1 + (lp1-1)*jd1*kd1;

                        // only accept new stencil if no hole/fringe pts.
                        if (iblank[ii1] > 0 &&
                            iblank[ii2] > 0 &&
                            iblank[ii3] > 0 &&
                            iblank[ii4] > 0 &&
                            iblank[ii5] > 0 &&
                            iblank[ii6] > 0 &&
                            iblank[ii7] > 0 &&
                            iblank[ii8] > 0) {
                            if (nbltmp[i] == 0) nok = nok + 1;
                            wherec.nblkpt[i] = -(m1+mdim);
                            nbltmp[i] = 1;
                        } else {
                            intrp2.jn[i] = jnsav;
                            intrp2.kn[i] = knsav;
                            intrp2.ln[i] = lnsav;
                            intrp1.xi[i] = xisav;
                            intrp1.yi[i] = yisav;
                            intrp1.zi[i] = zisav;
                        }
                    }
                }
            }
        }
        std::fprintf(funit(6),
            "        %d of these illegal stencils were successfully replaced\n", nok);
    }
    return;
}

// ============================================================================
// orphan — handle orphan points (nearest neighbor / extrapolation)
// ============================================================================
void orphan(int m,int n,int nserch,int i1,int iorph,int iflg)
{
    double* x = grid1.x;
    double* y = grid1.y;
    double* z = grid1.z;
    int* jimage = image.jimage;
    int* kimage = image.kimage;
    int* limage = image.limage;

    // local automatic arrays ds(idim),nbltmp(idim) — 1-based
    static double ds[idim+1];
    static int    nbltmp[idim+1];

    for (int i=i1; i<=intrp2.itotal; ++i) {
        nbltmp[i] = 0;
        ds[i]     = 1.e30;
    }

    for (int nn=1; nn<=nserch; ++nn) {
        int m1;
        if (iflg == 0) {
            m1 = conec1.LHOLE(m,n,nn);
        } else {
            int nnn = 1;
            m1 = conec1.LOUTR(m,nnn,nn);
        }
        int jd1 = grid2.mjmax[m1];
        int kd1 = grid2.mkmax[m1];
        int ld1 = grid2.mlmax[m1];
        getgrd(m1,x,y,z,jimage,kimage,limage,jd1,kd1,ld1);
        for (int i=i1; i<=intrp2.itotal; ++i) {
            if (wherec.nblkpt[i] == 0) {
                int j  = intrp2.jb[i];
                int k  = intrp2.kb[i];
                int l  = intrp2.lb[i];
                (void)j; (void)k; (void)l;
                double xp = surf2.xbo[i];
                double yp = surf2.ybo[i];
                double zp = surf2.zbo[i];

                // minimum distance search on mesh m1 faces only
                int js = grdlim.jsrs[m1];
                int je = grdlim.jsre[m1];
                int ks = grdlim.ksrs[m1];
                int ke = grdlim.ksre[m1];
                int ls = grdlim.lsrs[m1];
                int le = grdlim.lsre[m1];
                int jp, kp, lp;
                double dmin1;
                dsmin2(jd1,kd1,ld1,x,y,z,
                       xp,yp,zp,jp,kp,lp,js,je,ks,ke,ls,le,dmin1);

                if (dmin1 < ds[i]) {
                    ds[i] = dmin1;

                    int jpold = jp;
                    int kpold = kp;
                    int lpold = lp;

                    int iok;
                    int   isav = 0;
                    int   jnsav = 0, knsav = 0, lnsav = 0, nbsav = 0;
                    double xisav = 0., yisav = 0., zisav = 0.;

                    if (iorph == 0) {
                        // use nearest neighbor (zeroth order interpolation)
                        intrp2.jn[i] = jp;
                        intrp2.kn[i] = kp;
                        intrp2.ln[i] = lp;
                        intrp1.xi[i] = 0.;
                        intrp1.yi[i] = 0.;
                        intrp1.zi[i] = 0.;
                        iok = 1;
                    } else {
                        // extrapolate
                        // keep within bounds 1 < j,k,l < j/k/ld1-1
                        jp = std::min( jp , jd1-1 );
                        kp = std::min( kp , kd1-1 );
                        lp = std::min( lp , ld1-1 );
                        jp = std::max( 1 , jp );
                        kp = std::max( 1 , kp );
                        lp = std::max( 1 , lp );

                        // save extrapolation coefficients from the last block checked
                        isav = 0;
                        if (nbltmp[i] != 0) {
                            isav  = 1;
                            jnsav = intrp2.jn[i];
                            knsav = intrp2.kn[i];
                            lnsav = intrp2.ln[i];
                            xisav = intrp1.xi[i];
                            yisav = intrp1.yi[i];
                            zisav = intrp1.zi[i];
                            nbsav = nbltmp[i];
                        }

                        extrap(i,m,m1,jd1,kd1,ld1,x,y,z,xp,yp,zp,iok,
                               jp,kp,lp,jimage,kimage,limage);
                        if (iok == 0) {
                            // linear extrapolation not successful...use zeroth order
                            intrp2.jn[i] = jpold;
                            intrp2.kn[i] = kpold;
                            intrp2.ln[i] = lpold;
                            intrp1.xi[i] = 0.;
                            intrp1.yi[i] = 0.;
                            intrp1.zi[i] = 0.;
                        }

                        // for iorph=2, check extrapolation coefficients; if
                        // outside -0.5..1.5, use nearest neighbor instead
                        if (iorph == 2) {
                            if ((intrp1.xi[i] < -0.5-tol.epsc || intrp1.xi[i] > 1.5+tol.epsc) ||
                                (intrp1.yi[i] < -0.5-tol.epsc || intrp1.yi[i] > 1.5+tol.epsc) ||
                                (intrp1.zi[i] < -0.5-tol.epsc || intrp1.zi[i] > 1.5+tol.epsc)) {
                                if (iok > 0) {
                                    intrp2.jn[i] = jpold;
                                    intrp2.kn[i] = kpold;
                                    intrp2.ln[i] = lpold;
                                    intrp1.xi[i] = 0.;
                                    intrp1.yi[i] = 0.;
                                    intrp1.zi[i] = 0.;
                                }
                            }
                        }

                        iok = 1;
                    }
                    if (iok > 0) nbltmp[i] = -(m1+mdim);

                    // check current stencil against previous for best (smallest)
                    if (isav > 0) {
                        // saved stencil
                        double ximod0, yimod0, zimod0;
                        if (xisav-1.0 > tol.epsc) {
                            ximod0 = std::abs(xisav-1.0);
                        } else if (xisav < -tol.epsc) {
                            ximod0 = std::abs(xisav);
                        } else {
                            ximod0 = 0.;
                        }
                        if (yisav-1.0 > tol.epsc) {
                            yimod0 = std::abs(yisav-1.0);
                        } else if (yisav < -tol.epsc) {
                            yimod0 = std::abs(yisav);
                        } else {
                            yimod0 = 0.;
                        }
                        if (zisav-1.0 > tol.epsc) {
                            zimod0 = std::abs(zisav-1.0);
                        } else if (zisav < -tol.epsc) {
                            zimod0 = std::abs(zisav);
                        } else {
                            zimod0 = 0.;
                        }
                        double exmax0 = std::max(std::max(ximod0,yimod0),zimod0);

                        // new stencil
                        double ximod1, yimod1, zimod1;
                        if (intrp1.xi[i]-1.0 > tol.epsc) {
                            ximod1 = std::abs(intrp1.xi[i]-1.0);
                        } else if (intrp1.xi[i] < -tol.epsc) {
                            ximod1 = std::abs(intrp1.xi[i]);
                        } else {
                            ximod1 = 0.;
                        }
                        if (intrp1.yi[i]-1.0 > tol.epsc) {
                            yimod1 = std::abs(intrp1.yi[i]-1.0);
                        } else if (intrp1.yi[i] < -tol.epsc) {
                            yimod1 = std::abs(intrp1.yi[i]);
                        } else {
                            yimod1 = 0.;
                        }
                        if (intrp1.zi[i]-1.0 > tol.epsc) {
                            zimod1 = std::abs(intrp1.zi[i]-1.0);
                        } else if (intrp1.zi[i] < -tol.epsc) {
                            zimod1 = std::abs(intrp1.zi[i]);
                        } else {
                            zimod1 = 0.;
                        }
                        double exmax1 = std::max(std::max(ximod1,yimod1),zimod1);

                        if (exmax0 < exmax1) {
                            // old stencil was actually better
                            std::fprintf(funit(99),
                                " i = %d using stencil from mesh  %d not current mesh  %d\n",
                                i, -(mdim+nbsav), m1);
                            std::fprintf(funit(99),
                                "   old xi,yi,zi:   %g %g %g\n", xisav, yisav, zisav);
                            std::fprintf(funit(99),
                                "   new xi,yi,zi:   %g %g %g\n",
                                intrp1.xi[i], intrp1.yi[i], intrp1.zi[i]);
                            intrp1.xi[i] = xisav;
                            intrp1.yi[i] = yisav;
                            intrp1.zi[i] = zisav;
                            intrp2.jn[i] = jnsav;
                            intrp2.kn[i] = knsav;
                            intrp2.ln[i] = lnsav;
                            nbltmp[i] = nbsav;
                        }
                    }
                }
            }
        }
    }

    for (int i=i1; i<=intrp2.itotal; ++i) {
        if (wherec.nblkpt[i] == 0) wherec.nblkpt[i] = nbltmp[i];
    }
    return;
}
