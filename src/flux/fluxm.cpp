// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "fluxm.h"
#include "ccomplex.h"
#include "q8sdot.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace fluxm_ns {

void fluxm(FortranArray1DRef<double> ax, FortranArray1DRef<double> ay, FortranArray1DRef<double> az,
           FortranArray1DRef<double> are, FortranArray1DRef<double> at, FortranArray2DRef<double> qe,
           FortranArray2DRef<double> f, int& n, FortranArray2DRef<double> t, int& jkpro, int& nvtq,
           FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    // COMMON /fluid/
    float& gamma_f = cmn_fluid.gamma;
    float& gm1_f   = cmn_fluid.gm1;
    float& gp1_f   = cmn_fluid.gp1;
    float& gm1g_f  = cmn_fluid.gm1g;
    float& gp1g_f  = cmn_fluid.gp1g;
    float& ggm1_f  = cmn_fluid.ggm1;

    double gamma = (double)gamma_f;
    double gm1   = (double)gm1_f;
    double gp1   = (double)gp1_f;
    double gm1g  = (double)gm1g_f;
    double gp1g  = (double)gp1g_f;
    double ggm1  = (double)ggm1_f;

    double g1 = 1.e0 / (gamma * gamma - 1.e0);
    double g2 = ggm1 * 0.5e0;
    double c1 = 1.e0 / gamma;
    double c2 = ggm1 * g1;
    double c3 = 2.e0 * g1;
    double c4 = 2.e0 * c1;
    double c5 = 4.e0 * g1;
    double c6 = 2.e0 * g1 * gm1;
    double c7 = 2.e0 - gamma;
    double c8 = .5e0 * gm1;
    double zeroc = 0.;
    double onec  = 1.;

    // Loop 1000
    for (int izz = 1; izz <= n; izz++) {
        t(izz, 2) = ax(izz) * qe(izz, 2) + ay(izz) * qe(izz, 3) + az(izz) * qe(izz, 4)
                    + at(izz);
        t(izz, 3) = 1.e0 / qe(izz, 1);
        t(izz, 1) = t(izz, 2) * qe(izz, 1);

        t(izz, 4) = qe(izz, 2) * qe(izz, 2) + qe(izz, 3) * qe(izz, 3) +
                    qe(izz, 4) * qe(izz, 4);
        t(izz, 5) = gamma * qe(izz, 5) * t(izz, 3);
    }

    double dprint = 0.e0;

    // Loop 1001
    for (int izz = 1; izz <= n; izz++) {
        t(izz, 6) = 0.e0;
        bool cond_1001 = ((float)t(izz, 5) < 0.e0);
        t(izz, 6) = ccomplex_ns::ccvmgt(onec, t(izz, 6), cond_1001);
    }

    // q8sdot(n, t(1,6), n, t(1,6))
    {
        FortranArray1DRef<double> t_col6(&t(1, 6), nvtq, 1);
        dprint = q8sdot_ns::q8sdot(n, t_col6, n, t_col6);
    }

    if ((float)dprint > 0.e0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, " %20.10E%10d%10d",
                      (float)dprint, n, jkpro);
    }

    if ((float)dprint > 0.e0) {
        // Loop 1002
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 6) = t(izz, 5);
            bool cond_1002 = ((float)t(izz, 6) < 0.e0);
            qe(izz, 1) = ccomplex_ns::ccvmgt(f(izz, 1), qe(izz, 1), cond_1002);
            qe(izz, 2) = ccomplex_ns::ccvmgt(f(izz, 2), qe(izz, 2), cond_1002);
            qe(izz, 3) = ccomplex_ns::ccvmgt(f(izz, 3), qe(izz, 3), cond_1002);
            qe(izz, 4) = ccomplex_ns::ccvmgt(f(izz, 4), qe(izz, 4), cond_1002);
            qe(izz, 5) = ccomplex_ns::ccvmgt(f(izz, 5), qe(izz, 5), cond_1002);

            double tmp_t2_new = ax(izz) * f(izz, 2) + ay(izz) * f(izz, 3) + az(izz) * f(izz, 4);
            t(izz, 2) = ccomplex_ns::ccvmgt(tmp_t2_new, t(izz, 2), cond_1002);

            double tmp_t3_new = 1.e0 / qe(izz, 1);
            t(izz, 3) = ccomplex_ns::ccvmgt(tmp_t3_new, t(izz, 3), cond_1002);

            double tmp_t1_new = t(izz, 2) * qe(izz, 1);
            t(izz, 1) = ccomplex_ns::ccvmgt(tmp_t1_new, t(izz, 1), cond_1002);

            double tmp_t4_new = qe(izz, 2) * qe(izz, 2) + qe(izz, 3) * qe(izz, 3) +
                                qe(izz, 4) * qe(izz, 4);
            t(izz, 4) = ccomplex_ns::ccvmgt(tmp_t4_new, t(izz, 4), cond_1002);

            double tmp_t5_new = gamma * qe(izz, 5) * t(izz, 3);
            t(izz, 5) = ccomplex_ns::ccvmgt(tmp_t5_new, t(izz, 5), cond_1002);
        }
    }

    // Loop 1003
    for (int izz = 1; izz <= n; izz++) {
        t(izz, 6) = std::sqrt(t(izz, 5));
        qe(izz, 5) = qe(izz, 5) / gm1 + 0.5e0 * t(izz, 4) * qe(izz, 1);
        t(izz, 7) = t(izz, 2) / t(izz, 6);

        t(izz, 8) = .5e0 * t(izz, 7) - 0.5e0;
        f(izz, 1) = -are(izz) * qe(izz, 1) * t(izz, 6) * t(izz, 8) * t(izz, 8);

        t(izz, 8) = -c1 * (t(izz, 2) + 2.e0 * t(izz, 6));
        f(izz, 2) = ax(izz) * t(izz, 8) + qe(izz, 2);
        f(izz, 3) = ay(izz) * t(izz, 8) + qe(izz, 3);
        f(izz, 4) = az(izz) * t(izz, 8) + qe(izz, 4);

        f(izz, 5) = c2 * t(izz, 2) * t(izz, 8) + c3 * t(izz, 5) + .5e0 * t(izz, 4)
                    - at(izz) * t(izz, 8);

        f(izz, 2) = f(izz, 1) * f(izz, 2);
        f(izz, 3) = f(izz, 1) * f(izz, 3);
        f(izz, 4) = f(izz, 1) * f(izz, 4);
        f(izz, 5) = f(izz, 1) * f(izz, 5);

        double tmp_t12 = t(izz, 7);
        t(izz, 12) = ccomplex_ns::ccabs(tmp_t12);
        t(izz, 11) = 0.e0;
        bool cond_1003 = ((float)t(izz, 12) >= +1.e0);
        t(izz, 11) = ccomplex_ns::ccvmgt(onec, t(izz, 11), cond_1003);
    }

    double skip;
    {
        FortranArray1DRef<double> t_col11(&t(1, 11), nvtq, 1);
        skip = q8sdot_ns::q8sdot(n, t_col11, n, t_col11);
    }
    if ((float)skip < 0.5e0) return;

    // Loop 1004
    for (int izz = 1; izz <= n; izz++) {
        bool cond_le = ((float)t(izz, 7) <= -1.e0);
        bool cond_ge = ((float)t(izz, 7) >= +1.e0);

        double tmp_t11_new = ax(izz) * are(izz);
        t(izz, 11) = ccomplex_ns::ccvmgt(tmp_t11_new, t(izz, 11), cond_le);

        double tmp_t12_new = ay(izz) * are(izz);
        t(izz, 12) = ccomplex_ns::ccvmgt(tmp_t12_new, t(izz, 12), cond_le);

        double tmp_t13_new = az(izz) * are(izz);
        t(izz, 13) = ccomplex_ns::ccvmgt(tmp_t13_new, t(izz, 13), cond_le);

        double tmp_t14_new = at(izz) * are(izz);
        t(izz, 14) = ccomplex_ns::ccvmgt(tmp_t14_new, t(izz, 14), cond_le);

        double tmp_t2_new = t(izz, 2) * are(izz);
        t(izz, 2) = ccomplex_ns::ccvmgt(tmp_t2_new, t(izz, 2), cond_le);

        double tmp_f1_new = t(izz, 2) * qe(izz, 1);
        f(izz, 1) = ccomplex_ns::ccvmgt(tmp_f1_new, f(izz, 1), cond_le);

        double tmp_t6_new = qe(izz, 1) * t(izz, 5) * c1;
        t(izz, 6) = ccomplex_ns::ccvmgt(tmp_t6_new, t(izz, 6), cond_le);

        double tmp_f2_new = f(izz, 1) * qe(izz, 2) + t(izz, 6) * t(izz, 11);
        f(izz, 2) = ccomplex_ns::ccvmgt(tmp_f2_new, f(izz, 2), cond_le);

        double tmp_f3_new = f(izz, 1) * qe(izz, 3) + t(izz, 6) * t(izz, 12);
        f(izz, 3) = ccomplex_ns::ccvmgt(tmp_f3_new, f(izz, 3), cond_le);

        double tmp_f4_new = f(izz, 1) * qe(izz, 4) + t(izz, 6) * t(izz, 13);
        f(izz, 4) = ccomplex_ns::ccvmgt(tmp_f4_new, f(izz, 4), cond_le);

        double tmp_f5_new = t(izz, 2) * (qe(izz, 5) + t(izz, 6)) - t(izz, 14) * t(izz, 6);
        f(izz, 5) = ccomplex_ns::ccvmgt(tmp_f5_new, f(izz, 5), cond_le);

        f(izz, 1) = ccomplex_ns::ccvmgt(zeroc, f(izz, 1), cond_ge);
        f(izz, 2) = ccomplex_ns::ccvmgt(zeroc, f(izz, 2), cond_ge);
        f(izz, 3) = ccomplex_ns::ccvmgt(zeroc, f(izz, 3), cond_ge);
        f(izz, 4) = ccomplex_ns::ccvmgt(zeroc, f(izz, 4), cond_ge);
        f(izz, 5) = ccomplex_ns::ccvmgt(zeroc, f(izz, 5), cond_ge);
    }

    return;
}

} // namespace fluxm_ns
