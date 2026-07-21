// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "tinvr.h"
#include <cstdio>
#include <cmath>

namespace tinvr_ns {

void tinvr(int& n, FortranArray1DRef<double> r1, FortranArray1DRef<double> r2,
           FortranArray1DRef<double> r3, FortranArray1DRef<double> r4,
           FortranArray1DRef<double> r5, FortranArray1DRef<double> kx,
           FortranArray1DRef<double> ky, FortranArray1DRef<double> kz,
           FortranArray1DRef<double> lx, FortranArray1DRef<double> ly,
           FortranArray1DRef<double> lz, FortranArray1DRef<double> mx,
           FortranArray1DRef<double> my, FortranArray1DRef<double> mz,
           FortranArray1DRef<double> c, FortranArray1DRef<double> ub,
           FortranArray1DRef<double> rho, FortranArray1DRef<double> u,
           FortranArray1DRef<double> v, FortranArray1DRef<double> w,
           int& max, int& itinv, FortranArray1DRef<double> eig2,
           FortranArray1DRef<double> eig3, FortranArray1DRef<double> xm2a)
{
    // Access COMMON block fields
    float& gamma  = cmn_fluid.gamma;
    float& gm1    = cmn_fluid.gm1;
    float& cprec  = cmn_precond.cprec;

    // kx,ky,kz,lx,ly,lz,mx,my,mz are `real` in Fortran; this build is compiled
    // double precision, so they are 8-byte doubles stored in the double t workspace.
    // First loop: compute lx,ly,lz,mx,my,mz from kx,ky,kz
    for (int m = 1; m <= n; m++) {
        double lxi = ky(m) - kz(m);
        double lyi = kz(m) - kx(m);
        double lzi = kx(m) - ky(m);
        double mxi = lxi*lxi + lyi*lyi + lzi*lzi;
        mxi = 1.0e0 / std::sqrt(mxi);

        lx(m) = lxi * mxi;
        ly(m) = lyi * mxi;
        lz(m) = lzi * mxi;
        mx(m) = ky(m)*lz(m) - ly(m)*kz(m);
        my(m) = -kx(m)*lz(m) + lx(m)*kz(m);
        mz(m) = kx(m)*ly(m) - lx(m)*ky(m);
    }

    if ((float)cprec == 0.0f) {

        if (itinv == 0) {
            // T(inverse)*R
            for (int m = 1; m <= n; m++) {
                double vb = u(m)*lx(m) + v(m)*ly(m) + w(m)*lz(m);
                double wb = u(m)*mx(m) + v(m)*my(m) + w(m)*mz(m);

                double q2 = 0.5*(u(m)*u(m) + v(m)*v(m) + w(m)*w(m));

                double t3 = 1.0e0 / c(m);
                double t1 = -q2*r1(m) - r5(m) + u(m)*r2(m) + v(m)*r3(m) + w(m)*r4(m);
                t1 = (double)gm1 * t1 * t3 * t3;
                double t2 = -ub(m)*r1(m) + kx(m)*r2(m) + ky(m)*r3(m) + kz(m)*r4(m);
                t2 = t2 * t3;
                t3 = -vb*r1(m) + lx(m)*r2(m) + ly(m)*r3(m) + lz(m)*r4(m);

                r3(m) = -wb*r1(m) + mx(m)*r2(m) + my(m)*r3(m) + mz(m)*r4(m);
                r1(m) = r1(m) + t1;
                r2(m) = t3;
                r4(m) = 0.5e0*(t2 - t1);
                r5(m) = r4(m) - t2;
            }

        } else {
            // T(inverse)*M*R
            for (int m = 1; m <= n; m++) {
                double t1 = 1.0e0 / c(m);

                r5(m) = r5(m) * t1 * t1;
                r1(m) = r1(m) - r5(m);

                double t2 = rho(m) * r2(m);
                double t3 = rho(m) * r3(m);
                double t4 = rho(m) * r4(m);

                r2(m) =           lx(m)*t2 + ly(m)*t3 + lz(m)*t4;
                r3(m) =           mx(m)*t2 + my(m)*t3 + mz(m)*t4;
                r4(m) = 0.5*(t1*(kx(m)*t2 + ky(m)*t3 + kz(m)*t4) + r5(m));
                r5(m) = -r4(m) + r5(m);
            }
        }

    } else {

        if (itinv == 0) {
            // T(inverse)*M*R (with preconditioning)
            for (int m = 1; m <= n; m++) {
                double t1 = 1.0e0 / c(m);
                double rrho = 1.0e0 / rho(m);
                double xm2 = xm2a(m) * t1;
                double xm2ar = 1.0 / xm2a(m);
                double fplus = (eig2(m) - ub(m)) * xm2ar;
                double fmins = -(eig3(m) - ub(m)) * xm2ar;

                double r11 = r1(m);
                double r21 = r2(m);
                double r31 = r3(m);
                double r41 = r4(m);
                double r51 = r5(m);

                double vmag1 = u(m)*u(m) + v(m)*v(m) + w(m)*w(m);
                double r5t = (double)gm1 * (0.5*vmag1*r11
                           - (u(m)*r21 + v(m)*r31 + w(m)*r41) + r51);

                // multiplication by inverse of precond. matrix
                r1(m) = r11 - (1.0 - xm2)*r5t*t1*t1;
                r2(m) = rrho*(-u(m)*r11 + r21);
                r3(m) = rrho*(-v(m)*r11 + r31);
                r4(m) = rrho*(-w(m)*r11 + r41);
                r5(m) = xm2 * r5t;

                // multiplication by T(inverse)
                r5t   = r5(m) * t1 * t1;
                r1(m) = r1(m) - r5t;

                double t2 = rho(m) * r2(m);
                double t3 = rho(m) * r3(m);
                double t4 = rho(m) * r4(m);

                r2(m) =           lx(m)*t2 + ly(m)*t3 + lz(m)*t4;
                r3(m) =           mx(m)*t2 + my(m)*t3 + mz(m)*t4;
                r4(m) =  0.5*(t1*(kx(m)*t2 + ky(m)*t3 + kz(m)*t4)
                           + r5t*fplus);
                r5(m) = -0.5*(t1*(kx(m)*t2 + ky(m)*t3 + kz(m)*t4)
                           - r5t*fmins);
            }

        } else {
            // T(inverse)*R (with preconditioning)
            for (int m = 1; m <= n; m++) {
                double t1 = 1.0e0 / c(m);
                double xm2ar = 1.0 / xm2a(m);
                double fplus = (eig2(m) - ub(m)) * xm2ar;
                double fmins = -(eig3(m) - ub(m)) * xm2ar;

                double r5t = r5(m) * t1 * t1;
                r1(m) = r1(m) - r5t;

                double t2 = rho(m) * r2(m);
                double t3 = rho(m) * r3(m);
                double t4 = rho(m) * r4(m);

                r2(m) =           lx(m)*t2 + ly(m)*t3 + lz(m)*t4;
                r3(m) =           mx(m)*t2 + my(m)*t3 + mz(m)*t4;
                r4(m) =  0.5*(t1*(kx(m)*t2 + ky(m)*t3 + kz(m)*t4)
                           + r5t*fplus);
                r5(m) = -0.5*(t1*(kx(m)*t2 + ky(m)*t3 + kz(m)*t4)
                           - r5t*fmins);
            }
        }

    }

    return;
}

} // namespace tinvr_ns
