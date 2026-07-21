// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "xe2.h"
#include "xe.h"
#include "termn8.h"
#include "ccomplex.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace xe2_ns {

void xe(int& jdim, int& kdim, int& nsub, int& l,
        FortranArray3DRef<double> x, FortranArray3DRef<double> y,
        FortranArray3DRef<double> z, FortranArray3DRef<double> xmid,
        FortranArray3DRef<double> ymid, FortranArray3DRef<double> zmid,
        FortranArray3DRef<double> xmide, FortranArray3DRef<double> ymide,
        FortranArray3DRef<double> zmide, int& jcell, int& kcell,
        double& xc, double& yc, double& zc, double& xie, double& eta,
        int& imiss, int& ifit, int& ic0,
        FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
        int& nbuf, int& ibufdim, int& myid)
{
    xe_ns::xe(jdim, kdim, nsub, l, x, y, z, xmid, ymid, zmid,
              xmide, ymide, zmide, jcell, kcell, xc, yc, zc,
              xie, eta, imiss, ifit, ic0, nou, bou, nbuf, ibufdim, myid);
}

void xe2(double& x1, double& x2, double& x3, double& x4,
         double& x5, double& x6, double& x7, double& x8, double& xc,
         double& y1, double& y2, double& y3, double& y4,
         double& y5, double& y6, double& y7, double& y8, double& yc,
         double& xie, double& eta, int& imiss, int& ifit,
         FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
         int& nbuf, int& ibufdim, int& myid)
{
    // COMMON blocks
    int32_t& iretry   = cmn_tacos.iretry;
    float&   epsc     = cmn_tol.epsc;
    int32_t& iexp     = cmn_zero.iexp;

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

    double a2, a3, a4, a5, a6, a7, a8;
    double b2, b3, b4, b5, b6, b7, b8;
    double df1, df2, df3;

    if (ifit == 1) {
        // coefficients for bi-linear fit
        a2 = dx2;
        a3 = dx4;
        a4 = dx3 - a2 - a3;
        a5 = 0.;
        a6 = 0.;
        a7 = 0.;
        a8 = 0.;
        b2 = dy2;
        b3 = dy4;
        b4 = dy3 - b2 - b3;
        b5 = 0.;
        b6 = 0.;
        b7 = 0.;
        b8 = 0.;
    }

    if (ifit == 2) {
        // coefficients for (degenerate) bi-quadradratic fit
        // (quadratic in xie and eta, but without the (xie**2)*(eta**2) term)
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

    // coefficients for quadratic fit in xie, linear fit in eta
    if (ifit == 3) {
        a3  = dx4;
        b3  = dy4;
        a2  = -dx2 + 4.*dx5;
        b2  = -dy2 + 4.*dy5;
        a5  = 2.*dx2  - 4.*dx5;
        b5  = 2.*dy2  - 4.*dy5;
        df1 = x3 - x2 - a3;
        df2 = x6 - x5 - a3;
        a4  =   -df1 + 4.*df2;
        a6  = 2.*df1 - 4.*df2;
        df1 = y3 - y2 - b3;
        df2 = y6 - y5 - b3;
        b4  =   -df1 + 4.*df2;
        b6  = 2.*df1 - 4.*df2;
        a7  = 0.;
        a8  = 0.;
        b7  = 0.;
        b8  = 0.;
    }

    // coefficients for linear fit in xie, quadratic fit in eta
    if (ifit == 4) {
        a2  = dx2;
        a3  = -dx4 + 4.*dx7;
        a7  = 2.*dx4 -4.*dx7;
        df1 = dx3 - a2 - a3 - a7;
        df2 = dx8 - a2 -.5*a3 - .25*a7;
        a4  = -df1 + 4.*df2;
        a8  = 2.*df1 - 4.*df2;
        a5  = 0.;
        a6  = 0.;
        b2  = dy2;
        b3  = -dy4 + 4.*dy7;
        b7  = 2.*dy4 -4.*dy7;
        df1 = dy3 - b2 - b3 - b7;
        df2 = dy8 - b2 -.5*b3 - .25*b7;
        b4  = -df1 + 4.*df2;
        b8  = 2.*df1 - 4.*df2;
        b5  = 0.;
        b6  = 0.;
    }

    // newton iteration to invert for xie and eta (using two
    // of the three equations defining the surface patch)
    if (iretry == 0) {
        xie = 0.5;
        eta = 0.5;
    } else {
        int jj = (int)xie;
        int kk = (int)eta;
        xie = xie - jj;
        eta = eta - kk;
    }
    int limit = 50;
    int iter  = 0;
    int itry  = 0;

    // convergence criterion for Newton iteration...require 4 orders of
    // reduction to starting error=abs(f1)+abs(f2)
    double f1 = x1 + a3*eta + eta*( a7*eta + a8*xie*eta )
              + xie*( a2 + a4*eta + a5*xie + a6*xie*eta ) - xc;
    double f2 = y1 + b3*eta + eta*( b7*eta + b8*xie*eta )
              + xie*( b2 + b4*eta + b5*xie + b6*xie*eta ) - yc;
    double error0 = ccomplex_ns::ccabs(f1) + ccomplex_ns::ccabs(f2);
    double epsf   = 1.0e-4 * error0;
    // 10.**(-iexp) is machine zero
    double epsf1  = std::max(1.0e-09, std::pow(10., (double)(-iexp+1)));
    if ((float)epsf < (float)epsf1) epsf = epsf1;

    double errm1 = error0;
    double error = error0;
    double errm2 = error0;

    label2:
    {
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
        double b3b = b3 + xie*( b4 + b6*xie )
                   + 2.*eta*( b7 + b8*xie );

        iter  = iter + 1;
        double det = 1./( a2b*b3b - b2b*a3b );
        xie   = xie - det*( b3b*f1 - a3b*f2 );
        eta   = eta - det*( a2b*f2 - b2b*f1 );

        // exit newton iteration if huge xie or eta encountered
        double huge = 1.e6;
        if (std::abs((float)xie) >= (float)huge || std::abs((float)eta) >= (float)huge) {
            if (std::abs((float)xie) >= (float)huge) xie = huge;
            if (std::abs((float)eta) >= (float)huge) eta = huge;
            imiss = 1;
            return;
        }

        error = ccomplex_ns::ccabs(f1) + ccomplex_ns::ccabs(f2);
        if ((float)error > (float)epsf && iter <= limit) goto label2;

        // check to ensure point is inside cell
        imiss = 0;
        if ((float)xie < -(float)epsc || (float)xie > 1.+(float)epsc ||
            (float)eta < -(float)epsc || (float)eta > 1.+(float)epsc) {
            imiss = 1;
        }
        if (iter > limit && imiss == 0) {
            // newton iteration may be hung, but this may be the correct
            // cell. this behavior is sometimes seen when the code is run in
            // single precision. typically, 4 orders reduction in starting
            // error is required; if iteration is hung and at least 2 orders
            // have been obtained, then just accept it and print a warning. If
            // less than 2 orders reduction in error have been obtained, stop.
            // errm1 = error at iter-1; errm2 = error at iter-2
            double tmp_aerrm1 = error - errm1;
            double tmp_aerrm2 = error - errm2;
            double aerrm1 = ccomplex_ns::ccabs(tmp_aerrm1);
            double aerrm2 = ccomplex_ns::ccabs(tmp_aerrm2);
            if ((float)aerrm1 <= (float)epsf || (float)aerrm2 <= (float)epsf) {
                // iteration is hung
                if ((float)error < (float)epsf * 100.) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120,
                        "WARNING: Newton iter. hung, but"
                        " >= 2 orders reduction in error was obtained");
                    return;
                } else {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "stopping...Newton iter. hung with"
                        " < 2 orders reduction in error");
                    int neg1 = -1;
                    termn8_ns::termn8(myid, neg1, ibufdim, nbuf, bou, nou);
                }
            } else {
                // not converging...must be wrong cell
                xie = xie + 1.;
                eta = eta + 1.;
                imiss = 1;
            }
        }
    }
    return;
}

} // namespace xe2_ns
