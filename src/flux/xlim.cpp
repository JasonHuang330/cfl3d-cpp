// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "xlim.h"
#include "ccomplex.h"
#include "prolim.h"
#include "prolim2.h"

namespace xlim_ns {

void xlim(double& xkap, int& n, FortranArray1DRef<double> x1, FortranArray1DRef<double> x2, FortranArray1DRef<double> xc, int& iswi, int& npts, int& leq)
{
    // COMMON /chk/ ichk
    int32_t& ichk = cmn_chk.ichk;

    int ibias = 1;
    double eps   = 1.0e-06;
    double phi   = (3. - xkap) / (1. - xkap + eps);
    double eps2  = 0.5e-06;
    double zeroc = 0.;

    if (xkap < -1.0) {
        // first order
        for (int izz = 1; izz <= n; izz++) {
            x2(izz) = 0.;
            x1(izz) = 0.;
        }
        return;
    }

    // bias gradients by average values (leq = 1 or 5)
    if (ibias == 1) {
        if (leq == 1 || leq == 5) {
            for (int izz = 1; izz <= n; izz++) {
                double x2avg   = xc(izz) + x2(izz) * 0.5;
                double x1avg   = xc(izz) - x1(izz) * 0.5;
                x2(izz) = x2(izz) / x2avg;
                x1(izz) = x1(izz) / x1avg;
            }
        }
    }

    if (iswi == -1) {
        // od - second order
        for (int izz = 1; izz <= n; izz++) {
            double x3      = ccomplex_ns::ccabs(x1(izz));
            double x4      = ccomplex_ns::ccabs(x2(izz));
            bool cond      = (x3 < x4);
            double tmp_half = x1(izz) * 0.5e0;
            x2(izz) = ccomplex_ns::ccvmgt(tmp_half, x2(izz), cond);
            x1(izz) = ccomplex_ns::ccvmgt(x2(izz), x1(izz), cond);
            x1(izz) = x2(izz);
        }
    } else if (iswi == 1) {
        // smooth (van albeda) limiter
        double term    = xkap * 4.e0;
        for (int izz = 1; izz <= n; izz++) {
            double x3      = x1(izz) * x1(izz) + x2(izz) * x2(izz) + eps;
            double x4      = (eps2 + x1(izz) * x2(izz)) / x3;
            x4      = x4 * 0.5e0;
            x3      = (x1(izz) + x2(izz)) * x4;
            x4      = (x2(izz) - x1(izz)) * x4 * term * x4;
            x2(izz) = x3 + x4;
            x1(izz) = x3 - x4;
        }
    } else if (iswi == 0) {
        // unlimited (kappa scheme)
        double term    = xkap * 0.25e0;
        for (int izz = 1; izz <= n; izz++) {
            double x3      = .25e0 * (x1(izz) + x2(izz));
            double x4      = term * (x2(izz) - x1(izz));
            x2(izz) = x3 + x4;
            x1(izz) = x3 - x4;
        }
    } else if (iswi == 2) {
        // limited (min-mod scheme)
        double term    = xkap * 0.25e0;
        for (int izz = 1; izz <= n; izz++) {
            double x4      = x1(izz) * x2(izz);
            bool cond4     = (x4 < 0.e0);
            x1(izz) = ccomplex_ns::ccvmgt(zeroc, x1(izz), cond4);
            x2(izz) = ccomplex_ns::ccvmgt(zeroc, x2(izz), cond4);
            x4      = phi * x1(izz);
            x4      = ccomplex_ns::ccabs(x4);
            double x3      = phi * x2(izz);
            x3      = ccomplex_ns::ccabs(x3);
            double x5      = ccomplex_ns::ccabs(x2(izz));
            bool cond_x4x5 = (x4 < x5);
            double tmp_phi1 = x1(izz) * phi;
            x2(izz) = ccomplex_ns::ccvmgt(tmp_phi1, x2(izz), cond_x4x5);
            x5      = ccomplex_ns::ccabs(x1(izz));
            bool cond_x3x5 = (x3 < x5);
            double tmp_phi2 = x2(izz) * phi;
            x1(izz) = ccomplex_ns::ccvmgt(tmp_phi2, x1(izz), cond_x3x5);
            x3      = .25e0 * (x1(izz) + x2(izz));
            x4      = term * (x2(izz) - x1(izz));
            x2(izz) = x3 + x4;
            x1(izz) = x3 - x4;
        }
    } else if (iswi == 3 || iswi == 4) {
        // tuned k=1/3 limiter - Spekreijse - Venkat  AIAA-90-0429
        double delx    = 10. / (double)npts;
        eps2    = delx * delx * delx;
        for (int izz = 1; izz <= n; izz++) {
            double t3      = x1(izz) * x1(izz);
            double t4      = x2(izz) * x2(izz);
            double t5      = x1(izz) * x2(izz);
            double t6      = x1(izz) + x2(izz);
            double term    = 0.5 * (t5 + eps2) / (2. * (t3 + t4) - t5 + 3. * eps2);
            x2(izz) = (x2(izz) + t6) * term;
            x1(izz) = (x1(izz) + t6) * term;
        }
    }

    // bias gradients by average values (leq = 1 or 5)
    if (ibias == 1) {
        if (leq == 1 || leq == 5) {
            for (int izz = 1; izz <= n; izz++) {
                x2(izz) = x2(izz) * xc(izz);
                x1(izz) = x1(izz) * xc(izz);
            }
        }
    }

    // cap density and pressure
    // - ensures that they stay positive
    // - limits their maxima to their stagnation values
    if (ichk == 2) {
        prolim_ns::prolim(n, x1, x2, xc, leq);
    }
    if (ichk == 3) {
        prolim2_ns::prolim2(n, x1, x2, xc, leq);
    }

    return;
}

} // namespace xlim_ns
