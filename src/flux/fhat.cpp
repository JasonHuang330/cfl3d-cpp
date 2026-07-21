// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "fhat.h"
#include "ccomplex.h"
#include <cmath>

namespace fhat_ns {

void fhat(FortranArray1DRef<double> ax, FortranArray1DRef<double> ay, FortranArray1DRef<double> az, FortranArray1DRef<double> area, FortranArray1DRef<double> at, FortranArray2DRef<double> f, FortranArray2DRef<double> qr, FortranArray2DRef<double> ql, int& n, int& nvtq)
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

    double zero, x1, c1;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    double t11, t12, t13, t14, t15, t16, t17, t18, t19;
    double f1, f2, f3, f4, f5;
    double cc, uu, vv, ww, epsaa, epsbb, epscc;

    // 10.**(-iexp) is machine zero
    zero   = std::pow(10.0, (double)(-iexp));
    epsa_l = 2.f * epsa_r;

    // delta q across faces t(1-5)
    x1 = (double)gamma / (double)gm1;
    c1 = 1.0e0 / (double)gm1;

    if ((float)cprec == 0.f) {
        // Standard Roe flux (no preconditioning)
        for (int i = 1; i <= n; i++) {
            t1 = qr(i,1) - ql(i,1);
            t2 = qr(i,2) - ql(i,2);
            t3 = qr(i,3) - ql(i,3);
            t4 = qr(i,4) - ql(i,4);

            // pressure and enthalpy
            t16     = 1.0e0 / qr(i,1);
            t5      = qr(i,5);
            qr(i,5) = x1 * qr(i,5) * t16 + 0.5e0 * (qr(i,2)*qr(i,2) + qr(i,3)*qr(i,3)
                                                     + qr(i,4)*qr(i,4));

            t15     = 1.0e0 / ql(i,1);
            t19     = ql(i,5);
            ql(i,5) = x1 * ql(i,5) * t15 + 0.5e0 * (ql(i,2)*ql(i,2) + ql(i,3)*ql(i,3)
                                                     + ql(i,4)*ql(i,4));

            // unsplit contributions  f(r)+f(l)
            t18 = ax(i)*qr(i,2) + ay(i)*qr(i,3) + az(i)*qr(i,4) + at(i);
            t17 = ax(i)*ql(i,2) + ay(i)*ql(i,3) + az(i)*ql(i,4) + at(i);
            t6  = t18 * qr(i,1);
            t7  = t17 * ql(i,1);
            f1  = t6 + t7;
            f2  = t6*qr(i,2) + t7*ql(i,2);
            f3  = t6*qr(i,3) + t7*ql(i,3);
            f4  = t6*qr(i,4) + t7*ql(i,4);
            f5  = t6*qr(i,5) + t7*ql(i,5);
            t8  = t5 + t19;
            f2  = f2 + ax(i)*t8;
            f3  = f3 + ay(i)*t8;
            f4  = f4 + az(i)*t8;
            f5  = f5 - at(i)*t8;

            // roe averaged variables
            t6 = qr(i,1) * t15;
            t7 = std::sqrt(t6);
            t6 = 1.0e0 / (1.0e0 + t7);
            t8 = t7 * t6;

            // average density
            qr(i,1) = ql(i,1) * t7;

            // u,v,w,h average
            t9  = ql(i,2)*t6 + qr(i,2)*t8;
            t10 = ql(i,3)*t6 + qr(i,3)*t8;
            t11 = ql(i,4)*t6 + qr(i,4)*t8;
            t12 = ql(i,5)*t6 + qr(i,5)*t8;

            // extract sound speed
            t6  = (t9*t9 + t10*t10 + t11*t11) * 0.5e0;
            t7  = (double)gm1 * (t12 - t6);
            t8  = std::sqrt(t7);

            t13 = t9*ax(i) + t10*ay(i) + t11*az(i);

            // rhoa*delta(ubar) , delta(p)/c2a
            ql(i,1) = qr(i,1) * (t18 - t17);
            ql(i,2) = (t5 - t19) / t7;

            t18 = t13 + at(i);
            t18 = ccomplex_ns::ccabs(t18);
            t19 = t13 + at(i) + t8;
            t19 = ccomplex_ns::ccabs(t19);
            t17 = t13 + at(i) - t8;
            t17 = ccomplex_ns::ccabs(t17);

            // limit eigenvalues a la Harten and Gnoffo (NASA TP-2953)
            if ((float)epsa_r > 0.f) {
                cc    = t8;
                uu    = ccomplex_ns::ccabs(t9);
                vv    = ccomplex_ns::ccabs(t10);
                ww    = ccomplex_ns::ccabs(t11);
                epsaa = (double)epsa_r * (cc + uu + vv + ww);
                {
                    double tmp_max = ccomplex_ns::ccmax(epsaa, zero);
                    epsbb = 0.25 / tmp_max;
                }
                epscc = 2.00 * epsaa;
                if (t18 < epscc) t18 = t18*t18*epsbb + epsaa;
                if (t17 < epscc) t17 = t17*t17*epsbb + epsaa;
                if (t19 < epscc) t19 = t19*t19*epsbb + epsaa;
            }

            t14 = t18 * (t1 - ql(i,2));
            t15 = 0.5e0 * (ql(i,2) + ql(i,1)/t8);
            t16 = (ql(i,2) - t15) * t17;
            t15 = t15 * t19;

            qr(i,2) = t18 * (t2*qr(i,1) - ax(i)*ql(i,1));
            qr(i,3) = t18 * (t3*qr(i,1) - ay(i)*ql(i,1));
            qr(i,4) = t18 * (t4*qr(i,1) - az(i)*ql(i,1));
            qr(i,5) = t9*qr(i,2) + t10*qr(i,3) + t11*qr(i,4);

            ql(i,1) = t14 + t15 + t16;
            ql(i,2) = t8 * (t15 - t16);

            f1 = f1 - ql(i,1);
            f2 = f2 - ql(i,1)*t9  - ax(i)*ql(i,2) - qr(i,2);
            f3 = f3 - ql(i,1)*t10 - ay(i)*ql(i,2) - qr(i,3);
            f4 = f4 - ql(i,1)*t11 - az(i)*ql(i,2) - qr(i,4);
            f5 = f5 - ql(i,1)*t12 - t13*ql(i,2)   - qr(i,5) + t7*c1*t14;

            // include factor one-half area
            t7     = 0.5e0 * area(i);
            f(i,1) = t7 * f1;
            f(i,2) = t7 * f2;
            f(i,3) = t7 * f3;
            f(i,4) = t7 * f4;
            f(i,5) = t7 * f5;
        }
    } else {
        // Preconditioned Roe flux (Weiss-Smith)
        double t8t, delp;
        double vmag1, vel2, vel, xm2, xmave;
        double tt1, tt2, fplus, fmins, fsum;

        for (int i = 1; i <= n; i++) {
            t1 = qr(i,1) - ql(i,1);
            t2 = qr(i,2) - ql(i,2);
            t3 = qr(i,3) - ql(i,3);
            t4 = qr(i,4) - ql(i,4);

            // pressure and enthalpy
            t16     = 1.0e0 / qr(i,1);
            t5      = qr(i,5);
            qr(i,5) = x1 * qr(i,5) * t16 + 0.5e0 * (qr(i,2)*qr(i,2) + qr(i,3)*qr(i,3)
                                                     + qr(i,4)*qr(i,4));

            t15     = 1.0e0 / ql(i,1);
            t19     = ql(i,5);
            ql(i,5) = x1 * ql(i,5) * t15 + 0.5e0 * (ql(i,2)*ql(i,2) + ql(i,3)*ql(i,3)
                                                     + ql(i,4)*ql(i,4));

            // unsplit contributions  f(r)+f(l)
            t18 = ax(i)*qr(i,2) + ay(i)*qr(i,3) + az(i)*qr(i,4) + at(i);
            t17 = ax(i)*ql(i,2) + ay(i)*ql(i,3) + az(i)*ql(i,4) + at(i);
            t6  = t18 * qr(i,1);
            t7  = t17 * ql(i,1);
            f1  = t6 + t7;
            f2  = t6*qr(i,2) + t7*ql(i,2);
            f3  = t6*qr(i,3) + t7*ql(i,3);
            f4  = t6*qr(i,4) + t7*ql(i,4);
            f5  = t6*qr(i,5) + t7*ql(i,5);
            t8  = t5 + t19;
            f2  = f2 + ax(i)*t8;
            f3  = f3 + ay(i)*t8;
            f4  = f4 + az(i)*t8;
            f5  = f5 - at(i)*t8;
            t8t = t8;
            delp = t5 - t19;

            // roe averaged variables
            t6 = qr(i,1) * t15;
            t7 = std::sqrt(t6);
            t6 = 1.0e0 / (1.0e0 + t7);
            t8 = t7 * t6;

            // average density
            qr(i,1) = ql(i,1) * t7;

            // u,v,w,h average
            t9  = ql(i,2)*t6 + qr(i,2)*t8;
            t10 = ql(i,3)*t6 + qr(i,3)*t8;
            t11 = ql(i,4)*t6 + qr(i,4)*t8;
            t12 = ql(i,5)*t6 + qr(i,5)*t8;

            // extract sound speed
            t6  = (t9*t9 + t10*t10 + t11*t11) * 0.5e0;
            t7  = (double)gm1 * (t12 - t6);
            t8  = std::sqrt(t7);

            t13 = t9*ax(i) + t10*ay(i) + t11*az(i);

            // compute preconditioning parameters
            {
                double tmp_abs_delp = ccomplex_ns::ccabs(delp);
                vmag1 = 2.0*t6 + 2.0*t7*tmp_abs_delp/t8t;
            }
            {
                double uref2 = (double)avn * (double)uref * (double)uref;
                vel2 = ccomplex_ns::ccmax(vmag1, uref2);
            }
            {
                double tmp_min = ccomplex_ns::ccmin(t7, vel2);
                vel = std::sqrt(tmp_min);
            }
            vel  = (double)cprec * vel + (1.0 - (double)cprec) * t8;
            xm2  = (vel/t8) * (vel/t8);
            xmave = t13 / t8;
            tt1  = 0.5 * (1.0 + xm2);
            tt2  = 0.5 * std::sqrt(xmave*xmave*(1.0 - xm2)*(1.0 - xm2) + 4.0*xm2);

            // rhoa*delta(ubar) , delta(p)/c2a
            ql(i,1) = qr(i,1) * (t18 - t17);
            ql(i,2) = (t5 - t19) / t7;

            t18 = t13 + at(i);
            t19 = tt1*t13 + tt2*t8 + at(i);
            t17 = tt1*t13 - tt2*t8 + at(i);
            fplus = (t19 - t18) / (xm2 * t8);
            fmins = -(t17 - t18) / (xm2 * t8);
            fsum  = 2.0 / (fplus + fmins) / xm2;
            t18 = ccomplex_ns::ccabs(t18);
            t19 = ccomplex_ns::ccabs(t19);
            t17 = ccomplex_ns::ccabs(t17);

            // limit eigenvalues a la Harten and Gnoffo (NASA TP-2953)
            if ((float)epsa_r > 0.f) {
                cc    = t8;
                uu    = ccomplex_ns::ccabs(t9);
                vv    = ccomplex_ns::ccabs(t10);
                ww    = ccomplex_ns::ccabs(t11);
                epsaa = (double)epsa_r * (cc + uu + vv + ww);
                {
                    double tmp_max = ccomplex_ns::ccmax(epsaa, zero);
                    epsbb = 0.25 / tmp_max;
                }
                epscc = 2.00 * epsaa;
                if (t18 < epscc) t18 = t18*t18*epsbb + epsaa;
                if (t17 < epscc) t17 = t17*t17*epsbb + epsaa;
                if (t19 < epscc) t19 = t19*t19*epsbb + epsaa;
            }

            t14 = t18 * (t1 - ql(i,2));
            t15 = 0.5e0 * t19 * (fplus*ql(i,2) + ql(i,1)/t8);
            t16 = 0.5e0 * t17 * (fmins*ql(i,2) - ql(i,1)/t8);

            qr(i,2) = t18 * (t2*qr(i,1) - ax(i)*ql(i,1));
            qr(i,3) = t18 * (t3*qr(i,1) - ay(i)*ql(i,1));
            qr(i,4) = t18 * (t4*qr(i,1) - az(i)*ql(i,1));
            qr(i,5) = t9*qr(i,2) + t10*qr(i,3) + t11*qr(i,4);

            ql(i,1) = t14 + fsum*(t15 + t16);
            ql(i,2) = fsum * t8 * xm2 * (fmins*t15 - fplus*t16);

            f1 = f1 - ql(i,1);
            f2 = f2 - ql(i,1)*t9  - ax(i)*ql(i,2) - qr(i,2);
            f3 = f3 - ql(i,1)*t10 - ay(i)*ql(i,2) - qr(i,3);
            f4 = f4 - ql(i,1)*t11 - az(i)*ql(i,2) - qr(i,4);
            f5 = f5 - ql(i,1)*t12 - t13*ql(i,2)   - qr(i,5) + t7*c1*t14;

            // include factor one-half area
            t7     = 0.5e0 * area(i);
            f(i,1) = t7 * f1;
            f(i,2) = t7 * f2;
            f(i,3) = t7 * f3;
            f(i,4) = t7 * f4;
            f(i,5) = t7 * f5;
        }
    }

    return;
}

} // namespace fhat_ns
