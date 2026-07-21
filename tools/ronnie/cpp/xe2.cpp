// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// xe2.cpp — faithful C++ translation of xe2.F (CFL3D ronnie).
// Set up coefficients for a local polynomial fit in the xie/eta directions
// (bilinear / biquadratic / quadratic-linear / linear-quadratic per ifit) and
// invert for xie,eta via Newton iteration using two of the three surface eqns.
//
// Sequential real build: real(x)->x, ccabs->fabs (already inline).
// /tacos/ iretry -> tacos_.iretry
// /tol/  epsc,epsc0,epsreen,epscoll -> tol_.epsc (slot0)  [named epsc here]
// /zero/ iexp    -> zero_.iexp
// bou(ibufdim,nbuf) column-major: (n,iu) -> [(n-1)+(iu-1)*ibufdim]; nou(nbuf).
#include "ron_common.h"

void xe2(double x1,double x2,double x3,double x4,double x5,double x6,
         double x7,double x8,double xc,double y1,double y2,double y3,
         double y4,double y5,double y6,double y7,double y8,double yc,
         double& xie,double& eta,int& imiss,int ifit,
         int* nou,FStr* bou,int nbuf,int ibufdim,int myid)
{
    double dx2 = x2 - x1;
    double dy2 = y2 - y1;
    double dx3 = x3 - x1;
    double dy3 = y3 - y1;
    double dx4 = x4 - x1;
    double dy4 = y4 - y1;
    double dx5 = x5 - x1;
    double dy5 = y5 - y1;
    double dx6 = x6 - x1;
    double dy6 = y6 - y1;
    double dx7 = x7 - x1;
    double dy7 = y7 - y1;
    double dx8 = x8 - x1;
    double dy8 = y8 - y1;

    double a2=0.,a3=0.,a4=0.,a5=0.,a6=0.,a7=0.,a8=0.;
    double b2=0.,b3=0.,b4=0.,b5=0.,b6=0.,b7=0.,b8=0.;
    double df1,df2,df3;

    if (ifit == 1) {
        // coefficients for bi-linear fit
        a2 = dx2;
        a3 = dx4;
        a4 = dx3 - a2 - a3;
        a5 = 0.; a6 = 0.; a7 = 0.; a8 = 0.;
        b2 = dy2;
        b3 = dy4;
        b4 = dy3 - b2 - b3;
        b5 = 0.; b6 = 0.; b7 = 0.; b8 = 0.;
    }

    if (ifit == 2) {
        // coefficients for (degenerate) bi-quadratic fit
        a2  = -dx2 + 4.*dx5;
        a3  = -dx4 + 4.*dx7;
        a5  = 2.*dx2 - 4.*dx5;
        a7  = 2.*dx4 - 4.*dx7;
        df1 = dx3 - a2 - a3 - a5 - a7;
        df2 = dx6 - .5*a2 - a3 - .25*a5 - a7;
        df3 = dx8 - a2 - .5*a3 - a5 - .25*a7;
        a4  = -3.*df1 + 4.*df2 + 4.*df3;
        a6  = 2.*df1 - 4.*df2;
        a8  = 2.*df1 - 4.*df3;
        b2  = -dy2 + 4.*dy5;
        b3  = -dy4 + 4.*dy7;
        b5  = 2.*dy2 - 4.*dy5;
        b7  = 2.*dy4 - 4.*dy7;
        df1 = dy3 - b2 - b3 - b5 - b7;
        df2 = dy6 - .5*b2 - b3 - .25*b5 - b7;
        df3 = dy8 - b2 - .5*b3 - b5 - .25*b7;
        b4  = -3.*df1 + 4.*df2 + 4.*df3;
        b6  = 2.*df1 - 4.*df2;
        b8  = 2.*df1 - 4.*df3;
    }

    if (ifit == 3) {
        // quadratic fit in xie, linear fit in eta
        a3  = dx4;
        b3  = dy4;
        a2  = -dx2 + 4.*dx5;
        b2  = -dy2 + 4.*dy5;
        a5  = 2.*dx2 - 4.*dx5;
        b5  = 2.*dy2 - 4.*dy5;
        df1 = x3 - x2 - a3;
        df2 = x6 - x5 - a3;
        a4  =   -df1 + 4.*df2;
        a6  = 2.*df1 - 4.*df2;
        df1 = y3 - y2 - b3;
        df2 = y6 - y5 - b3;
        b4  =   -df1 + 4.*df2;
        b6  = 2.*df1 - 4.*df2;
        a7  = 0.; a8  = 0.; b7  = 0.; b8  = 0.;
    }

    if (ifit == 4) {
        // linear fit in xie, quadratic fit in eta
        a2  = dx2;
        a3  = -dx4 + 4.*dx7;
        a7  = 2.*dx4 - 4.*dx7;
        df1 = dx3 - a2 - a3 - a7;
        df2 = dx8 - a2 - .5*a3 - .25*a7;
        a4  = -df1 + 4.*df2;
        a8  = 2.*df1 - 4.*df2;
        a5  = 0.; a6  = 0.;
        b2  = dy2;
        b3  = -dy4 + 4.*dy7;
        b7  = 2.*dy4 - 4.*dy7;
        df1 = dy3 - b2 - b3 - b7;
        df2 = dy8 - b2 - .5*b3 - .25*b7;
        b4  = -df1 + 4.*df2;
        b8  = 2.*df1 - 4.*df2;
        b5  = 0.; b6  = 0.;
    }

    // newton iteration to invert for xie and eta
    if (tacos_.iretry == 0) {
        xie = 0.5;
        eta = 0.5;
    } else {
        int jj = f_int(xie);   // jj = xie (real->int truncation toward 0)
        int kk = f_int(eta);   // kk = eta
        xie = xie - jj;
        eta = eta - kk;
    }
    int limit = 50;
    int iter  = 0;

    // convergence criterion for Newton iteration
    double f1 = x1 + a3*eta + eta*( a7*eta + a8*xie*eta )
              + xie*( a2 + a4*eta + a5*xie + a6*xie*eta ) - xc;
    double f2 = y1 + b3*eta + eta*( b7*eta + b8*xie*eta )
              + xie*( b2 + b4*eta + b5*xie + b6*xie*eta ) - yc;
    double error0 = ccabs(f1) + ccabs(f2);
    double epsf   = 1.0e-4*error0;
    // 10.**(-iexp) is machine zero
    double epsf1  = f_max(1.0e-09, std::pow(10.0,(double)(-zero_.iexp+1)));
    if (epsf < epsf1) epsf = epsf1;

    double errm1 = error0;
    double error = error0;
    double errm2;

L2:
    errm2 = errm1;
    errm1 = error;

    f1 = x1 + a3*eta + eta*( a7*eta + a8*xie*eta )
       + xie*( a2 + a4*eta + a5*xie + a6*xie*eta ) - xc;
    f2 = y1 + b3*eta + eta*( b7*eta + b8*xie*eta )
       + xie*( b2 + b4*eta + b5*xie + b6*xie*eta ) - yc;

    // d(f1)/d(xie)
    double a2b = a2 + 2.*xie*( a5 + a6*eta ) + eta*( a4 + a8*eta );
    // d(f1)/d(eta)
    double a3b = a3 + xie*( a4 + a6*xie ) + 2.*eta*( a7 + a8*xie );
    // d(f2)/d(xie)
    double b2b = b2 + 2.*xie*( b5 + b6*eta ) + eta*( b4 + b8*eta );
    // d(f2)/d(eta)
    double b3b = b3 + xie*( b4 + b6*xie ) + 2.*eta*( b7 + b8*xie );

    iter = iter + 1;
    double det = 1./( a2b*b3b - b2b*a3b );
    xie = xie - det*( b3b*f1 - a3b*f2 );
    eta = eta - det*( a2b*f2 - b2b*f1 );

    // exit newton iteration if huge xie or eta encountered
    double huge = 1.e6;
    if (std::fabs(xie) >= huge || std::fabs(eta) >= huge) {
        if (std::fabs(xie) >= huge) xie = huge;
        if (std::fabs(eta) >= huge) eta = huge;
        imiss = 1;
        return;
    }

    error = ccabs(f1) + ccabs(f2);
    if (error > epsf && iter <= limit) goto L2;

    // check to ensure point is inside cell
    imiss = 0;
    if (xie < -tol_.epsc || xie > 1.0+tol_.epsc ||
        eta < -tol_.epsc || eta > 1.0+tol_.epsc) {
        imiss = 1;
    }
    if (iter > limit && imiss == 0) {
        // newton iteration may be hung, but this may be the correct cell.
        double aerrm1 = ccabs(error-errm1);
        double aerrm2 = ccabs(error-errm2);
        if (aerrm1 <= epsf || aerrm2 <= epsf) {
            // iteration is hung
            if (error < epsf*100.) {
                nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
                    "WARNING: Newton iter. hung, but >= 2 orders reduction"
                    " in error was obtained");
                return;
            } else {
                nou[1-1] = i_min(nou[1-1]+1,ibufdim);
                bou[(nou[1-1]-1)+(1-1)*ibufdim].assign(
                    "stopping...Newton iter. hung with < 2 orders reduction"
                    " in error");
                termn8(myid,-1,ibufdim,nbuf,bou,nou);
            }
        } else {
            // not converging...must be wrong cell
            xie = xie + 1.;
            eta = eta + 1.;
            imiss = 1;
        }
    }
    return;
}
