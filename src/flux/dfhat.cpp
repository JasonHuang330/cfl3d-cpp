// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dfhat.h"
#include <cmath>
#include <algorithm>

namespace dfhat_ns {

void dfhat(FortranArray1DRef<double> ax, FortranArray1DRef<double> ay, FortranArray1DRef<double> az, FortranArray1DRef<double> area, FortranArray1DRef<double> at, FortranArray2DRef<double> q, FortranArray3DRef<double> df, int& nn, int& nvtq, int& ipm)
{
    // COMMON block references
    float& gamma  = cmn_fluid.gamma;
    float& gm1    = cmn_fluid.gm1;
    float& cprec  = cmn_precond.cprec;
    float& uref   = cmn_precond.uref;
    float& avn    = cmn_precond.avn;
    float& epsa_l = cmn_entfix.epsa_l;
    float& epsa_r = cmn_entfix.epsa_r;
    int32_t& iexp = cmn_zero.iexp;

    // 10.**(-iexp) is machine zero
    double zero_val = std::pow(10.0, -(double)iexp);
    epsa_l = 2.f * epsa_r;

    double gm1i = 1.0 / (double)gm1;

    double sign_val = 1.0e0;
    if (ipm < 0) sign_val = -1.0e0;

    for (int n = 1; n <= nn; n++) {
        double rhoi  = 1.0 / q(n, 1);
        double c2    = (double)gamma * q(n, 5) * rhoi;
        double c     = std::sqrt(c2);
        double c2i   = 1.0 / c2;
        double q22   = 0.5 * (q(n, 2)*q(n, 2) + q(n, 3)*q(n, 3) + q(n, 4)*q(n, 4));
        double h     = c2 * gm1i + q22;
        double ubar  = q(n, 2)*ax(n) + q(n, 3)*ay(n) + q(n, 4)*az(n) + at(n);
        double ar2   = 0.5 * area(n);

        // preconditioning additions
        double vmag1 = 2.0 * q22;
        double vel2  = std::max(vmag1, (double)avn * (double)uref * (double)uref);
        double vel   = std::sqrt(std::min(c2, vel2));
        vel          = (double)cprec * vel + (1.0 - (double)cprec) * c;
        double xm2   = (vel / c) * (vel / c);
        double xmave = ubar / c;
        double tt1   = 0.5 * (1.0 + xm2);
        double tt2   = 0.5 * std::sqrt(xmave*xmave * (1.0 - xm2)*(1.0 - xm2) + 4.0*xm2);

        double e1u = ubar;
        double e2u = tt1*ubar + tt2*c;
        double e3u = tt1*ubar - tt2*c;

        // limit eigenvalues a la Harten and Gnoffo (NASA TP-2953)
        double cc    = c;
        double uu    = std::abs(q(n, 2));
        double vv    = std::abs(q(n, 3));
        double ww    = std::abs(q(n, 4));
        double epsaa = (double)epsa_l * (cc + uu + vv + ww);
        double epsbb = 0.25 / std::max(epsaa, zero_val);
        double epscc = 2.00 * epsaa;
        double t17   = std::abs(ubar - c);
        double t18   = std::abs(ubar);
        double t19   = std::abs(ubar + c);
        if ((float)t18 < (float)epscc) t18 = t18*t18*epsbb + epsaa;
        if ((float)t17 < (float)epscc) t17 = t17*t17*epsbb + epsaa;
        if ((float)t19 < (float)epscc) t19 = t19*t19*epsbb + epsaa;

        double e1    = ar2 * (e1u + sign_val * t18);
        double e2    = ar2 * (e2u + sign_val * t19);
        double e3    = ar2 * (e3u + sign_val * t17);
        double fplus = (e2u - e1u) / (xm2 * c);
        double fmins = -(e3u - e1u) / (xm2 * c);
        double fsum  = 2.0 / (fplus + fmins) / xm2;

        double ruc1  = q(n, 1) * e1;
        double ruc2  = 0.5 * q(n, 1) * fsum * (e2 - e3) / c;
        double ruc3  = 0.5 * q(n, 1) * (fsum * (fplus*e2 + fmins*e3) - 2.0*e1);
        double ruc4  = 0.5 * q(n, 1) * (fsum * xm2 * (fmins*e2 + fplus*e3) - 2.0*e1);

        double rhor2 = rhoi * ruc2 * xm2 * fplus * fmins;
        double axr3  = ax(n) * ruc4;
        double ayr3  = ay(n) * ruc4;
        double azr3  = az(n) * ruc4;

        df(n, 1, 1) = e1;
        df(n, 1, 2) = ax(n) * ruc2;
        df(n, 1, 3) = ay(n) * ruc2;
        df(n, 1, 4) = az(n) * ruc2;
        df(n, 1, 5) = c2i * rhoi * ruc3;

        df(n, 2, 1) =        q(n, 2)*df(n, 1, 1);
        df(n, 2, 2) = ruc1 + q(n, 2)*df(n, 1, 2) + ax(n)*axr3;
        df(n, 2, 3) =        q(n, 2)*df(n, 1, 3) + ax(n)*ayr3;
        df(n, 2, 4) =        q(n, 2)*df(n, 1, 4) + ax(n)*azr3;
        df(n, 2, 5) =        q(n, 2)*df(n, 1, 5) + ax(n)*rhor2;

        df(n, 3, 1) =        q(n, 3)*df(n, 1, 1);
        df(n, 3, 2) =        q(n, 3)*df(n, 1, 2) + ay(n)*axr3;
        df(n, 3, 3) = ruc1 + q(n, 3)*df(n, 1, 3) + ay(n)*ayr3;
        df(n, 3, 4) =        q(n, 3)*df(n, 1, 4) + ay(n)*azr3;
        df(n, 3, 5) =        q(n, 3)*df(n, 1, 5) + ay(n)*rhor2;

        df(n, 4, 1) =        q(n, 4)*df(n, 1, 1);
        df(n, 4, 2) =        q(n, 4)*df(n, 1, 2) + az(n)*axr3;
        df(n, 4, 3) =        q(n, 4)*df(n, 1, 3) + az(n)*ayr3;
        df(n, 4, 4) = ruc1 + q(n, 4)*df(n, 1, 4) + az(n)*azr3;
        df(n, 4, 5) =        q(n, 4)*df(n, 1, 5) + az(n)*rhor2;

        df(n, 5, 1) =   df(n, 1, 1)*q22;
        df(n, 5, 2) = h*df(n, 1, 2) + (ubar - at(n))*axr3 + q(n, 2)*ruc1;
        df(n, 5, 3) = h*df(n, 1, 3) + (ubar - at(n))*ayr3 + q(n, 3)*ruc1;
        df(n, 5, 4) = h*df(n, 1, 4) + (ubar - at(n))*azr3 + q(n, 4)*ruc1;
        df(n, 5, 5) = h*df(n, 1, 5) + (ubar - at(n))*rhor2 + df(n, 1, 1)*gm1i;
    }
}

} // namespace dfhat_ns
