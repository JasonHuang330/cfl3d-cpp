// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "fluxp.h"
#include "runtime/fortran_io.h"
#include "ccomplex.h"
#include "q8sdot.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace fluxp_ns {

void fluxp(FortranArray1DRef<double> ax, FortranArray1DRef<double> ay, FortranArray1DRef<double> az, FortranArray1DRef<double> are, FortranArray1DRef<double> at, FortranArray2DRef<double> qe, FortranArray2DRef<double> f, int& n, FortranArray2DRef<double> t, int& jkpro, int& nvtq, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    // COMMON /fluid/
    float& gamma = cmn_fluid.gamma;
    float& gm1   = cmn_fluid.gm1;
    float& gp1   = cmn_fluid.gp1;
    float& gm1g  = cmn_fluid.gm1g;
    float& gp1g  = cmn_fluid.gp1g;
    float& ggm1  = cmn_fluid.ggm1;

    double g1 = 1.e0 / ((double)gamma * (double)gamma - 1.e0);
    double g2 = (double)ggm1 * 0.5e0;
    double c1 = 1.e0 / (double)gamma;
    double c2 = (double)ggm1 * g1;
    double c3 = 2.e0 * g1;
    double c4 = 2.e0 * c1;
    double c5 = 4.e0 * g1;
    double c6 = 2.e0 * g1 * (double)gm1;
    double c7 = 2.e0 - (double)gamma;
    double c8 = .5e0 * (double)gm1;
    double zeroc = 0.;
    double onec  = 1.;

    // Loop 1000
    for (int izz = 1; izz <= n; izz++) {
        t(izz,2) = ax(izz)*qe(izz,2) + ay(izz)*qe(izz,3) + az(izz)*qe(izz,4)
                   + at(izz);
        t(izz,3) = 1.e0 / qe(izz,1);
        t(izz,1) = t(izz,2) * qe(izz,1);

        t(izz,4) = qe(izz,2)*qe(izz,2) + qe(izz,3)*qe(izz,3) +
                   qe(izz,4)*qe(izz,4);
        t(izz,5) = (double)gamma * qe(izz,5) * t(izz,3);
    }

    double dprint = 0.e0;

    // Loop 1001
    for (int izz = 1; izz <= n; izz++) {
        t(izz,6) = 0.e0;
        // ccvmgt(onec, t(izz,6), (real(t(izz,5)).lt.0.e0))
        // returns onec if condition is true, t(izz,6) otherwise
        {
            double tmp_a = onec;
            double tmp_b = t(izz,6);
            bool   tmp_c = ((float)t(izz,5) < 0.e0f);
            t(izz,6) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, tmp_c);
        }
    }

    // dprint = q8sdot(n, t(1,6), n, t(1,6))
    {
        FortranArray1DRef<double> t_col6(&t(1,6), n, 1);
        dprint = q8sdot_ns::q8sdot(n, t_col6, n, t_col6);
    }

    if ((float)dprint > 0.e0f) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " %20.10E%10d%10d",
                      (float)dprint, n, jkpro);
    }

    if ((float)dprint > 0.e0f) {
        // Loop 1002
        for (int izz = 1; izz <= n; izz++) {
            t(izz,6) = t(izz,5);
            bool cond = ((float)t(izz,6) < 0.e0f);

            // qe(izz,1) = ccvmgt(f(izz,1), qe(izz,1), cond)
            {
                double tmp_a = f(izz,1);
                double tmp_b = qe(izz,1);
                qe(izz,1) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond);
            }
            {
                double tmp_a = f(izz,2);
                double tmp_b = qe(izz,2);
                qe(izz,2) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond);
            }
            {
                double tmp_a = f(izz,3);
                double tmp_b = qe(izz,3);
                qe(izz,3) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond);
            }
            {
                double tmp_a = f(izz,4);
                double tmp_b = qe(izz,4);
                qe(izz,4) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond);
            }
            {
                double tmp_a = f(izz,5);
                double tmp_b = qe(izz,5);
                qe(izz,5) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond);
            }
            {
                double tmp_a = ax(izz)*f(izz,2) + ay(izz)*f(izz,3) +
                               az(izz)*f(izz,4);
                double tmp_b = t(izz,2);
                t(izz,2) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond);
            }
            {
                double tmp_a = onec / qe(izz,1);
                double tmp_b = t(izz,3);
                t(izz,3) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond);
            }
            {
                double tmp_a = t(izz,2) * qe(izz,1);
                double tmp_b = t(izz,1);
                t(izz,1) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond);
            }

            {
                double tmp_a = qe(izz,2)*qe(izz,2) + qe(izz,3)*qe(izz,3) +
                               qe(izz,4)*qe(izz,4);
                double tmp_b = t(izz,4);
                t(izz,4) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond);
            }
            {
                double tmp_a = (double)gamma * qe(izz,5) * t(izz,3);
                double tmp_b = t(izz,5);
                t(izz,5) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond);
            }
        }
    }

    // Loop 1003
    for (int izz = 1; izz <= n; izz++) {
        t(izz,6)  = std::sqrt(t(izz,5));
        qe(izz,5) = qe(izz,5) / (double)gm1 + 0.5e0 * t(izz,4) * qe(izz,1);
        t(izz,7)  = t(izz,2) / t(izz,6);

        t(izz,8)  = .5e0 * t(izz,7) + 0.5e0;
        f(izz,1)  = are(izz) * qe(izz,1) * t(izz,6) * t(izz,8) * t(izz,8);

        t(izz,8)  = -c1 * (t(izz,2) - 2.e0 * t(izz,6));
        f(izz,2)  = ax(izz) * t(izz,8) + qe(izz,2);
        f(izz,3)  = ay(izz) * t(izz,8) + qe(izz,3);
        f(izz,4)  = az(izz) * t(izz,8) + qe(izz,4);

        f(izz,5)  = c2 * t(izz,2) * t(izz,8) + c3 * t(izz,5) + .5e0 * t(izz,4)
                    - at(izz) * t(izz,8);

        f(izz,2)  = f(izz,1) * f(izz,2);
        f(izz,3)  = f(izz,1) * f(izz,3);
        f(izz,4)  = f(izz,1) * f(izz,4);
        f(izz,5)  = f(izz,1) * f(izz,5);

        // t(izz,12) = ccabs(t(izz,7))
        {
            double tmp_a = t(izz,7);
            t(izz,12) = ccomplex_ns::ccabs(tmp_a);
        }
        t(izz,11) = 0.e0;
        // t(izz,11) = ccvmgt(onec, t(izz,11), (real(t(izz,12)).ge.+1.e00))
        {
            double tmp_a = onec;
            double tmp_b = t(izz,11);
            bool   tmp_c = ((float)t(izz,12) >= +1.e00f);
            t(izz,11) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, tmp_c);
        }
    }

    // skip = q8sdot(n, t(1,11), n, t(1,11))
    double skip;
    {
        FortranArray1DRef<double> t_col11(&t(1,11), n, 1);
        skip = q8sdot_ns::q8sdot(n, t_col11, n, t_col11);
    }
    if ((float)skip < 0.5e0f) return;

    // Loop 1004
    for (int izz = 1; izz <= n; izz++) {
        bool cond_ge = ((float)t(izz,7) >= +1.e00f);
        bool cond_le = ((float)t(izz,7) <= -1.e00f);

        {
            double tmp_a = ax(izz) * are(izz);
            double tmp_b = t(izz,11);
            t(izz,11) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_ge);
        }
        {
            double tmp_a = ay(izz) * are(izz);
            double tmp_b = t(izz,12);
            t(izz,12) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_ge);
        }
        {
            double tmp_a = az(izz) * are(izz);
            double tmp_b = t(izz,13);
            t(izz,13) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_ge);
        }
        {
            double tmp_a = at(izz) * are(izz);
            double tmp_b = t(izz,14);
            t(izz,14) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_ge);
        }
        {
            double tmp_a = t(izz,2) * are(izz);
            double tmp_b = t(izz,2);
            t(izz,2) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_ge);
        }

        {
            double tmp_a = t(izz,2) * qe(izz,1);
            double tmp_b = f(izz,1);
            f(izz,1) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_ge);
        }
        {
            double tmp_a = qe(izz,1) * t(izz,5) * c1;
            double tmp_b = t(izz,6);
            t(izz,6) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_ge);
        }
        {
            double tmp_a = f(izz,1) * qe(izz,2) + t(izz,6) * t(izz,11);
            double tmp_b = f(izz,2);
            f(izz,2) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_ge);
        }
        {
            double tmp_a = f(izz,1) * qe(izz,3) + t(izz,6) * t(izz,12);
            double tmp_b = f(izz,3);
            f(izz,3) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_ge);
        }
        {
            double tmp_a = f(izz,1) * qe(izz,4) + t(izz,6) * t(izz,13);
            double tmp_b = f(izz,4);
            f(izz,4) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_ge);
        }
        {
            double tmp_a = t(izz,2) * (qe(izz,5) + t(izz,6)) - t(izz,14) * t(izz,6);
            double tmp_b = f(izz,5);
            f(izz,5) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_ge);
        }

        {
            double tmp_a = zeroc;
            double tmp_b = f(izz,1);
            f(izz,1) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_le);
        }
        {
            double tmp_a = zeroc;
            double tmp_b = f(izz,2);
            f(izz,2) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_le);
        }
        {
            double tmp_a = zeroc;
            double tmp_b = f(izz,3);
            f(izz,3) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_le);
        }
        {
            double tmp_a = zeroc;
            double tmp_b = f(izz,4);
            f(izz,4) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_le);
        }
        {
            double tmp_a = zeroc;
            double tmp_b = f(izz,5);
            f(izz,5) = ccomplex_ns::ccvmgt(tmp_a, tmp_b, cond_le);
        }
    }

    return;
}

} // namespace fluxp_ns
