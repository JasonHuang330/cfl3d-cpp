// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "tdq.h"

namespace tdq_ns {

void tdq(int& n, FortranArray1DRef<double> r1, FortranArray1DRef<double> r2,
         FortranArray1DRef<double> r3, FortranArray1DRef<double> r4,
         FortranArray1DRef<double> r5, FortranArray1DRef<double> kx,
         FortranArray1DRef<double> ky, FortranArray1DRef<double> kz,
         FortranArray1DRef<double> lx, FortranArray1DRef<double> ly,
         FortranArray1DRef<double> lz, FortranArray1DRef<double> mx,
         FortranArray1DRef<double> my, FortranArray1DRef<double> mz,
         FortranArray1DRef<double> c, FortranArray1DRef<double> ub,
         FortranArray1DRef<double> rho, FortranArray1DRef<double> u,
         FortranArray1DRef<double> v, FortranArray1DRef<double> w,
         int& max, FortranArray1DRef<double> eig2,
         FortranArray1DRef<double> eig3, FortranArray1DRef<double> xm2a)
{
    // common /precond/ cprec, uref, avn
    float& cprec = cmn_precond.cprec;

    // M(inverse)*T*R
    if (cprec == 0.0f) {
        // cdir$ ivdep
        for (int m = 1; m <= n; m++) {
            double t1 = 1.0 / rho(m);
            double t2 = t1 * r2(m);
            double t3 = t1 * r3(m);
            double t5 = t1 * c(m) * (r4(m) - r5(m));

            r5(m) = r4(m) + r5(m);
            r1(m) = r1(m) + r5(m);
            r5(m) = c(m) * c(m) * r5(m);

            r2(m) = lx(m) * t2 + mx(m) * t3 + kx(m) * t5;
            r3(m) = ly(m) * t2 + my(m) * t3 + ky(m) * t5;
            r4(m) = lz(m) * t2 + mz(m) * t3 + kz(m) * t5;
        }
    } else {
        // cdir$ ivdep
        for (int m = 1; m <= n; m++) {
            // modifications for preconditioning
            double xm2ar = 1.0 / xm2a(m);
            double fplus = (eig2(m) - ub(m)) * xm2ar;
            double fmins = -(eig3(m) - ub(m)) * xm2ar;
            double fsum  = 2.0 / (fplus + fmins);

            double t1 = 1.0 / rho(m);
            double t2 = t1 * r2(m);
            double t3 = t1 * r3(m);
            double t5 = t1 * c(m) * (fmins * r4(m) - fplus * r5(m)) * fsum;

            r5(m) = (r4(m) + r5(m)) * fsum;
            r1(m) = r1(m) + r5(m);
            r5(m) = c(m) * c(m) * r5(m);

            r2(m) = lx(m) * t2 + mx(m) * t3 + kx(m) * t5;
            r3(m) = ly(m) * t2 + my(m) * t3 + ky(m) * t5;
            r4(m) = lz(m) * t2 + mz(m) * t3 + kz(m) * t5;
        }
    }
    return;
}

} // namespace tdq_ns
