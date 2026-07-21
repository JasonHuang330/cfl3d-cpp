// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// extrae.cpp — faithful C++ translation of extrae.F (CFL3D ronnie).
// Compute extra mid-cell point (x7,y7,z7) at (xie,eta)=(0.,.5) via a
// quadratic least-squares extrapolation along the k (eta) direction.
//
// Sequential real build: real(x)->x.  Preserves exact arithmetic + branch
// structure (Fortran labels 1000/2000/1500 -> C++ gotos).
// Arrays x/y/z (jdim,kdim,nsub): (i,j,k) -> [(i-1)+(j-1)*jdim+(k-1)*jdim*kdim].
#include "ron_common.h"

void extrae(int jdim,int kdim,int nsub,int l,
            double* x,double* y,double* z,int jcell,int kcell,int kcl,int kcr,
            double& x7,double& y7,double& z7,int& icase,int ifit)
{
    (void)nsub; (void)ifit;
    const int d1 = jdim, d2 = kdim;
    #define X(i,j,k) x[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y(i,j,k) y[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z(i,j,k) z[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]

    double x1,y1,z1,x2,y2,z2;
    double dxp,dyp,dzp,xlen;
    double dxq,dyq,dzq,xlenq,check;
    double dxr,dyr,dzr,xlenr;
    double xlens,b1,c1,d1_,b2,c2,d2_,a1,a2,trat,term;

    // check if only bilinear fit reasonable
    icase = 999;
    x1 = X(jcell,kcell,l);
    y1 = Y(jcell,kcell,l);
    z1 = Z(jcell,kcell,l);
    x2 = X(jcell,kcell+1,l);
    y2 = Y(jcell,kcell+1,l);
    z2 = Z(jcell,kcell+1,l);
    x7 = 0.5*( x1 + x2 );
    y7 = 0.5*( y1 + y2 );
    z7 = 0.5*( z1 + z2 );
    if (kcell < kcl || kcell > kcr) goto L1500;
    if (kcl == kcr) goto L1500;
    dxp = x2 - x1;
    dyp = y2 - y1;
    dzp = z2 - z1;
    xlen = std::sqrt( dxp*dxp + dyp*dyp + dzp*dzp );
    if (xlen <= 0.0) goto L1500;
    if (kcell == kcl) goto L1000;
    if (kcell == kcr) goto L2000;

    // interior points
    dxq = x1 - X(jcell,kcell-1,l);
    dyq = y1 - Y(jcell,kcell-1,l);
    dzq = z1 - Z(jcell,kcell-1,l);
    xlenq = std::sqrt( dxq*dxq + dyq*dyq + dzq*dzq );
    if (xlenq <= 0.0) goto L1500;
    check = ( dxq*dxp + dyq*dyp + dzq*dzp )/(xlenq*xlen);
    if (check < 0.707107) goto L2000;
    xlenq = xlenq / xlen;
    if (xlenq < 0.333 || xlenq > 3.0) goto L2000;
    xlens = -1.0;
    b2 = X(jcell,kcell-1,l) - x1 - (x2-x1)*xlens;
    c2 = Y(jcell,kcell-1,l) - y1 - (y2-y1)*xlens;
    d2_ = Z(jcell,kcell-1,l) - z1 - (z2-z1)*xlens;
    a2 = xlens*( 1.0 - xlens );
    x7 = x7 + 0.25*b2/a2;
    y7 = y7 + 0.25*c2/a2;
    z7 = z7 + 0.25*d2_/a2;
    dxr = X(jcell,kcell+2,l) - x2;
    dyr = Y(jcell,kcell+2,l) - y2;
    dzr = Z(jcell,kcell+2,l) - z2;
    xlenr = std::sqrt( dxr*dxr + dyr*dyr + dzr*dzr );
    if (xlenr <= 0.0) goto L1500;
    check = ( dxr*dxp + dyr*dyp + dzr*dzp )/(xlenr*xlen);
    if (check < 0.707107) goto L1500;
    xlenr = xlenr / xlen;
    if (xlenr < 0.333 || xlenr > 3.0) goto L1500;
    xlens = 2.0;
    b1 = X(jcell,kcell+2,l) - x1 - (x2-x1)*xlens;
    c1 = Y(jcell,kcell+2,l) - y1 - (y2-y1)*xlens;
    d1_ = Z(jcell,kcell+2,l) - z1 - (z2-z1)*xlens;
    a1 = xlens*( 1.0 - xlens );
    trat = xlenq*xlenq/(xlenr*xlenr);
    trat = 1.0;
    term = 0.25/( a1*a1 + trat*a2*a2 );
    x7 = 0.5*( x1 + x2 ) + term*( b1*a1 + b2*a2*trat );
    y7 = 0.5*( y1 + y2 ) + term*( c1*a1 + c2*a2*trat );
    z7 = 0.5*( z1 + z2 ) + term*( d1_*a1 + d2_*a2*trat );
    icase = 0;
    goto L1500;

L1000:
    // left edge
    icase = -1;
    dxr = X(jcell,kcell+2,l) - x2;
    dyr = Y(jcell,kcell+2,l) - y2;
    dzr = Z(jcell,kcell+2,l) - z2;
    xlenr = std::sqrt( dxr*dxr + dyr*dyr + dzr*dzr );
    if (xlenr <= 0.0) goto L1500;
    check = ( dxr*dxp + dyr*dyp + dzr*dzp )/(xlenr*xlen);
    if (check < 0.707107) goto L1500;
    xlenr = xlenr / xlen;
    if (xlenr < 0.333 || xlenr > 3.0) goto L1500;
    xlens = 2.0;
    b1 = X(jcell,kcell+2,l) - x1 - (x2-x1)*xlens;
    c1 = Y(jcell,kcell+2,l) - y1 - (y2-y1)*xlens;
    d1_ = Z(jcell,kcell+2,l) - z1 - (z2-z1)*xlens;
    a1 = 1.0/( xlens*(1.0 - xlens) );
    x7 = x7 + 0.25*( b1 * a1 );
    y7 = y7 + 0.25*( c1 * a1 );
    z7 = z7 + 0.25*( d1_ * a1 );
    icase = 1;
    goto L1500;

L2000:
    // right edge
    icase = -2;
    dxq = x1 - X(jcell,kcell-1,l);
    dyq = y1 - Y(jcell,kcell-1,l);
    dzq = z1 - Z(jcell,kcell-1,l);
    xlenq = std::sqrt( dxq*dxq + dyq*dyq + dzq*dzq );
    if (xlenq <= 0.0) goto L1500;
    check = ( dxq*dxp + dyq*dyp + dzq*dzp )/(xlenq*xlen);
    if (check < 0.707107) goto L1500;
    xlenq = xlenq / xlen;
    if (xlenq < 0.333 || xlenq > 3.0) goto L1500;
    xlens = -1.0;
    b2 = X(jcell,kcell-1,l) - x1 - (x2-x1)*xlens;
    c2 = Y(jcell,kcell-1,l) - y1 - (y2-y1)*xlens;
    d2_ = Z(jcell,kcell-1,l) - z1 - (z2-z1)*xlens;
    a2 = 1.0/( xlens*( 1.0 - xlens ) );
    x7 = x7 + 0.25*( b2 * a2 );
    y7 = y7 + 0.25*( c2 * a2 );
    z7 = z7 + 0.25*( d2_ * a2 );
    icase = 2;

L1500:
    return;
    #undef X
    #undef Y
    #undef Z
}
