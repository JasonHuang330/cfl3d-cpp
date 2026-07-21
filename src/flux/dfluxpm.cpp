// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dfluxpm.h"
#include "ccomplex.h"
#include <cmath>

namespace dfluxpm_ns {

void dfluxpm(FortranArray1DRef<double> s1, FortranArray1DRef<double> s2,
             FortranArray1DRef<double> ax, FortranArray1DRef<double> ay,
             FortranArray1DRef<double> az, FortranArray1DRef<double> area,
             FortranArray1DRef<double> at, FortranArray2DRef<double> q,
             FortranArray3DRef<double> df, int& n, int& jkpro, int& nvtq, int& ipm)
{
    // COMMON /fluid/
    double gamma  = (double)cmn_fluid.gamma;
    double gm1    = (double)cmn_fluid.gm1;
    double gp1    = (double)cmn_fluid.gp1;

    double gam1   = gamma / gm1;
    double gam21i = 1.0 / (gamma * gamma - 1.0);
    double zeroc  = 0.0;

    double sign_val = 0.0;
    if (ipm > 0) sign_val = +1.0;
    if (ipm < 0) sign_val = -1.0;

    for (int l = 1; l <= n; l++) {
        double q2    = q(l,2)*q(l,2) + q(l,3)*q(l,3) + q(l,4)*q(l,4);
        double q22   = 0.5 * q2;
        double epp   = gam1 * q(l,5) + q(l,1) * q22;
        double a     = std::sqrt(gamma * q(l,5) / q(l,1));
        double a2    = a * a;
        double am2   = 2.0 * a;
        double am2i  = 1.0 / am2;
        double ubar  = ax(l)*q(l,2) + ay(l)*q(l,3) + az(l)*q(l,4) + at(l);
        double t1    = at(l) / gamma;

        double xm    = sign_val * ubar / a;

        // s1(l) = ccvmgt(area(l), zeroc, (abs(real(xm)) .lt. 1.0))
        {
            bool cond1 = (std::abs((float)xm) < 1.0f);
            s1(l) = ccomplex_ns::ccvmgt(area(l), zeroc, cond1);
        }
        // s2(l) = ccvmgt(area(l), zeroc, (real(xm) .ge. 1.0))
        {
            bool cond2 = ((float)xm >= 1.0f);
            s2(l) = ccomplex_ns::ccvmgt(area(l), zeroc, cond2);
        }

        double ua    =  ubar + sign_val * a;
        double ub    = -ubar + sign_val * a;

        double u1    = ubar * am2i;
        double u2    = am2 - u1 * ub;
        double u3    = ubar / gp1;

        double ua1   = s1(l) * 0.25 * ua;
        double ua2   = sign_val * ua * am2i;
        double ua3   = 3.0 * ua2 - 1.0;
        double ua4   = 1.0 - ua2;
        double ua5   = s1(l) * q(l,1) * ua2;
        double ua6   = ua1 / a2;
        double ua7   = gamma * ua6;
        double ua8   = 0.5 * ua5;
        double ua8x  = ua8 * ax(l);
        double ua8y  = ua8 * ay(l);
        double ua8z  = ua8 * az(l);
        double ua8a  = ua8 * ua;

        double ub1   = 3.0 * ub / gamma;
        double ub2   = ub1 * ax(l);
        double ub3   = ub2 + 2.0 * q(l,2);
        double ub4   = ub1 * ay(l);
        double ub5   = ub4 + 2.0 * q(l,3);
        double ub6   = ub1 * az(l);
        double ub7   = ub6 + 2.0 * q(l,4);
        double ub8   = 2.0 * a2 / gm1 + ub * (4.0 * u3 - 3.0 * t1) + q2;

        double f1    = s2(l) * ubar;
        double f2    = s2(l) * ax(l);
        double f3    = s2(l) * ay(l);
        double f4    = s2(l) * az(l);
        double f5    = f1 * q(l,1);
        double f6    = f2 * q(l,1);
        double f7    = f3 * q(l,1);
        double f8    = f4 * q(l,1);

        df(l,1,1) = ua1 * ua3                                          + f1;
        df(l,1,2) = ua5 * ax(l)                                        + f6;
        df(l,1,3) = ua5 * ay(l)                                        + f7;
        df(l,1,4) = ua5 * az(l)                                        + f8;
        df(l,1,5) = ua7 * ua4;

        df(l,2,1) = ua1 * (sign_val * ub2 * u1 + q(l,2) * ua3)        + f1 * q(l,2);
        df(l,2,2) = ua8x * ub3 + ua8a                                  + f6 * q(l,2) + f5;
        df(l,2,3) = ua8y * ub3                                         + f7 * q(l,2);
        df(l,2,4) = ua8z * ub3                                         + f8 * q(l,2);
        df(l,2,5) = ua6 * (sign_val * ax(l) * u2 + gamma * q(l,2) * ua4) + f2;

        df(l,3,1) = ua1 * (sign_val * ub4 * u1 + q(l,3) * ua3)        + f1 * q(l,3);
        df(l,3,2) = ua8x * ub5                                         + f6 * q(l,3);
        df(l,3,3) = ua8y * ub5 + ua8a                                  + f7 * q(l,3) + f5;
        df(l,3,4) = ua8z * ub5                                         + f8 * q(l,3);
        df(l,3,5) = ua6 * (sign_val * ay(l) * u2 + gamma * q(l,3) * ua4) + f3;

        df(l,4,1) = ua1 * (sign_val * ub6 * u1 + q(l,4) * ua3)        + f1 * q(l,4);
        df(l,4,2) = ua8x * ub7                                         + f6 * q(l,4);
        df(l,4,3) = ua8y * ub7                                         + f7 * q(l,4);
        df(l,4,4) = ua8z * ub7 + ua8a                                  + f8 * q(l,4) + f5;
        df(l,4,5) = ua6 * (sign_val * az(l) * u2 + gamma * q(l,4) * ua4) + f4;

        df(l,5,1) = ua1 * (sign_val * ub * (3.0 * u1 * (u3 - t1) -
                    a * gam21i) + q22 * ua3)                           + f1 * q22;
        df(l,5,2) = ua8x * ub8 + q(l,2) * ua8a                        + f2 * epp + f5 * q(l,2);
        df(l,5,3) = ua8y * ub8 + q(l,3) * ua8a                        + f3 * epp + f5 * q(l,3);
        df(l,5,4) = ua8z * ub8 + q(l,4) * ua8a                        + f4 * epp + f5 * q(l,4);
        df(l,5,5) = ua7 * (sign_val * (u3 - t1) * u2 +
                    a * gam21i * (am2 + sign_val * ua) + q22 * ua4)    + f1 * gam1 - s2(l) * at(l);
    }
}

} // namespace dfluxpm_ns
