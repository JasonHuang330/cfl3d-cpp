// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "extra.h"
#include <cmath>

namespace extra_ns {

void extra(int& jdim, int& kdim, int& nsub, int& l,
           FortranArray3DRef<double> x, FortranArray3DRef<double> y,
           FortranArray3DRef<double> z,
           int& jcell, int& kcell, int& jcl, int& jcr,
           double& x5, double& y5, double& z5,
           int& icase, int& ifit)
{
    // Declare ALL locals at top to avoid "jump over initialization" errors with goto
    double x1, y1, z1, x2, y2, z2;
    double dxp, dyp, dzp, xlen;
    double dxq, dyq, dzq, xlenq;
    double dxr, dyr, dzr, xlenr;
    double check, xlens;
    double b1, c1, d1, a1;
    double b2, c2, d2, a2;
    double trat, term;

    //
    // check if only bilinear fit reasonable
    //
    icase = 999;
    x1    = x(jcell,   kcell, l);
    y1    = y(jcell,   kcell, l);
    z1    = z(jcell,   kcell, l);
    x2    = x(jcell+1, kcell, l);
    y2    = y(jcell+1, kcell, l);
    z2    = z(jcell+1, kcell, l);
    x5    = 0.5*( x1 + x2 );
    y5    = 0.5*( y1 + y2 );
    z5    = 0.5*( z1 + z2 );
//  if(ifit.eq.1 .or. ifit.eq.3) return
    if (jcell < jcl || jcell > jcr) goto label1500;
    if (jcl == jcr) goto label1500;
    dxp   = x2 - x1;
    dyp   = y2 - y1;
    dzp   = z2 - z1;
    xlen  = std::sqrt( dxp*dxp + dyp*dyp + dzp*dzp );
    if ((float)xlen <= 0) goto label1500;
    if (jcell == jcl) goto label1000;
    if (jcell == jcr) goto label2000;
    //
    // interior points
    //
    dxq   = x1 - x(jcell-1, kcell, l);
    dyq   = y1 - y(jcell-1, kcell, l);
    dzq   = z1 - z(jcell-1, kcell, l);
    xlenq = std::sqrt( dxq*dxq + dyq*dyq + dzq*dzq );
    if ((float)xlenq <= 0) goto label1500;
    check = ( dxq*dxp + dyq*dyp + dzq*dzp ) / (xlenq*xlen);
    if ((float)check < 0.707107) goto label2000;
    xlenq = xlenq / xlen;
    if ((float)xlenq < 0.333 || (float)xlenq > 3.) goto label2000;
    xlens = -1.;
    b2    = x(jcell-1, kcell, l) - x1 - (x2-x1)*xlens;
    c2    = y(jcell-1, kcell, l) - y1 - (y2-y1)*xlens;
    d2    = z(jcell-1, kcell, l) - z1 - (z2-z1)*xlens;
    a2    = xlens*( 1. - xlens );
    x5    = x5 + .25*b2/a2;
    y5    = y5 + .25*c2/a2;
    z5    = z5 + .25*d2/a2;
    dxr   = x(jcell+2, kcell, l) - x2;
    dyr   = y(jcell+2, kcell, l) - y2;
    dzr   = z(jcell+2, kcell, l) - z2;
    xlenr = std::sqrt( dxr*dxr + dyr*dyr + dzr*dzr );
    if ((float)xlenr <= 0) goto label1500;
    check = ( dxr*dxp + dyr*dyp + dzr*dzp ) / (xlenr*xlen);
    if ((float)check < 0.707107) goto label1500;
    xlenr = xlenr / xlen;
    if ((float)xlenr < 0.333 || (float)xlenr > 3.) goto label1500;
    xlens = 2.;
    b1    = x(jcell+2, kcell, l) - x1 - (x2-x1)*xlens;
    c1    = y(jcell+2, kcell, l) - y1 - (y2-y1)*xlens;
    d1    = z(jcell+2, kcell, l) - z1 - (z2-z1)*xlens;
    a1    = xlens*( 1. - xlens );
    trat  = xlenq*xlenq/(xlenr*xlenr);
    trat  = 1.;
    term  = .25/( a1*a1 + trat*a2*a2 );
    x5    = 0.5*( x1 + x2 ) + term*( b1*a1 + b2*a2*trat );
    y5    = 0.5*( y1 + y2 ) + term*( c1*a1 + c2*a2*trat );
    z5    = 0.5*( z1 + z2 ) + term*( d1*a1 + d2*a2*trat );
    icase = 0;
    goto label1500;
label1000:;
    //
    // left edge
    //
    icase = -1;
    dxr   = x(jcell+2, kcell, l) - x2;
    dyr   = y(jcell+2, kcell, l) - y2;
    dzr   = z(jcell+2, kcell, l) - z2;
    xlenr = std::sqrt( dxr*dxr + dyr*dyr + dzr*dzr );
    if ((float)xlenr <= 0) goto label1500;
    check = ( dxr*dxp + dyr*dyp + dzr*dzp ) / (xlenr*xlen);
    if ((float)check < 0.707107) goto label1500;
    xlenr = xlenr / xlen;
    if ((float)xlenr < 0.333 || (float)xlenr > 3.) goto label1500;
    xlens = 2.;
    b1    = x(jcell+2, kcell, l) - x1 - (x2-x1)*xlens;
    c1    = y(jcell+2, kcell, l) - y1 - (y2-y1)*xlens;
    d1    = z(jcell+2, kcell, l) - z1 - (z2-z1)*xlens;
    a1    = 1./( xlens*(1. - xlens) );
    x5    = x5 + .25*( b1 * a1 );
    y5    = y5 + .25*( c1 * a1 );
    z5    = z5 + .25*( d1 * a1 );
    icase = 1;
    goto label1500;
label2000:;
    //
    // right edge
    //
    icase = -2;
    dxq   = x1 - x(jcell-1, kcell, l);
    dyq   = y1 - y(jcell-1, kcell, l);
    dzq   = z1 - z(jcell-1, kcell, l);
    xlenq = std::sqrt( dxq*dxq + dyq*dyq + dzq*dzq );
    if ((float)xlenq <= 0) goto label1500;
    check = ( dxq*dxp + dyq*dyp + dzq*dzp ) / (xlenq*xlen);
    if ((float)check < 0.707107) goto label1500;
    xlenq = xlenq / xlen;
    if ((float)xlenq < 0.333 || (float)xlenq > 3.) goto label1500;
    xlens = -1.;
    b2    = x(jcell-1, kcell, l) - x1 - (x2-x1)*xlens;
    c2    = y(jcell-1, kcell, l) - y1 - (y2-y1)*xlens;
    d2    = z(jcell-1, kcell, l) - z1 - (z2-z1)*xlens;
    a2    = 1./( xlens*( 1. - xlens ) );
    x5    = x5 + .25*( b2 * a2 );
    y5    = y5 + .25*( c2 * a2 );
    z5    = z5 + .25*( d2 * a2 );
    icase = 2;
label1500:;
    return;
}

} // namespace extra_ns
