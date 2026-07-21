// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rie1d.h"
#include "ccf.h"
#include "ccomplex.h"
#include <cmath>

namespace rie1d_ns {

void rie1d(int& jvdim, FortranArray2DRef<double> t, int& jv, double& cl) {
//
//***********************************************************************
//     Purpose:  Determine far-field boundary data using quasi 1-d
//     characteristic relations.
//     Modified for Weiss-Smith preconditioning by J.R. Edwards, NCSU
//       cprec = 0 ---> original code used
//             > 0 ---> modified code used
//***********************************************************************
//
    // COMMON block aliases
    float& gamma   = cmn_fluid.gamma;
    float& gm1     = cmn_fluid.gm1;
    float& p0      = cmn_ivals.p0;
    float& rho0    = cmn_ivals.rho0;
    float& c0      = cmn_ivals.c0;
    float& u0      = cmn_ivals.u0;
    float& v0      = cmn_ivals.v0;
    float& w0      = cmn_ivals.w0;
    float& xmach   = cmn_info.xmach;
    float& alpha   = cmn_info.alpha;
    int32_t& iipv  = cmn_info.iipv;
    float& xmc     = cmn_fsum.xmc;
    float& zmc     = cmn_fsum.zmc;
    float& cprec   = cmn_precond.cprec;

    if ((float)cprec == 0.f) {
//
// ---- use Riemann invariants to determine conditions
//
        double uf = (double)u0;
        double vf = (double)v0;
        double wf = (double)w0;
        double cf = (double)c0;
        for (int izz = 1; izz <= jv; izz++) {
            t(izz, 21) = (double)u0;
            t(izz, 22) = (double)w0;
            t(izz, 23) = (double)c0;
        }
//   Point vortex hardwire addition:
//   Must be 2-D!  Must be x-z plane!  Uses xmc, zmc for location for vortex
        if (iipv == 1) {
            double xcenter = (double)xmc;
            double zcenter = (double)zmc;
            double cosa = std::cos((double)alpha);
            double sina = std::sin((double)alpha);
            double pi = std::acos(-1.0);
            for (int izz = 1; izz <= jv; izz++) {
                double xa = t(izz, 18) - xcenter;
                double za = t(izz, 19) - zcenter;
                double xmach_d = (double)xmach;
                ccf_ns::ccf(xa, za, cosa, sina, cl, xmach_d,
                            t(izz, 21), t(izz, 22), t(izz, 23), pi);
            }
        }
//
        double x2gm1 = 2.e0 / (double)gm1;
        for (int izz = 1; izz <= jv; izz++) {
            t(izz, 10) = (double)gamma * t(izz, 5) / t(izz, 1);
        }
//
        for (int izz = 1; izz <= jv; izz++) {
            t(izz, 11) = t(izz, 6) * t(izz, 2) + t(izz, 7) * t(izz, 3) + t(izz, 8) * t(izz, 4);
            t(izz, 12) = std::sqrt(t(izz, 10));
            t(izz, 12) = t(izz, 11) + x2gm1 * t(izz, 12);
            t(izz, 13) = t(izz, 6) * t(izz, 21) + t(izz, 7) * vf + t(izz, 8) * t(izz, 22);
            t(izz, 14) = t(izz, 13) - x2gm1 * t(izz, 23);
        }
//
        double x4gm1 = 0.25e0 * (double)gm1;
        for (int izz = 1; izz <= jv; izz++) {
            t(izz, 15) = 0.5e0 * (t(izz, 12) + t(izz, 14));
            t(izz, 16) = x4gm1 * (t(izz, 12) - t(izz, 14));
            t(izz, 17) = t(izz, 15) - t(izz, 13);
//  put unsteady metrics into t(5)
            t(izz, 5)  = t(izz, 20);
            t(izz, 18) = t(izz, 21) + t(izz, 6) * t(izz, 17);
            t(izz, 19) = vf + t(izz, 7) * t(izz, 17);
            t(izz, 20) = t(izz, 22) + t(izz, 8) * t(izz, 17);
        }
//
        double ent0 = (double)gamma * (double)p0 / std::pow((double)rho0, (double)gm1);
        for (int izz = 1; izz <= jv; izz++) {
            t(izz, 9)  = ent0;
            t(izz, 12) = 1.e0 / std::pow(t(izz, 1), (double)gm1);
            t(izz, 13) = t(izz, 15) + t(izz, 5);
            bool cond17 = ((float)t(izz, 13) >= 0.e0);
            double a17 = t(izz, 15) - t(izz, 11);
            t(izz, 17) = ccomplex_ns::ccvmgt(a17, t(izz, 17), cond17);
            bool cond18 = ((float)t(izz, 13) >= 0.e0);
            double a18 = t(izz, 2) + t(izz, 6) * t(izz, 17);
            t(izz, 18) = ccomplex_ns::ccvmgt(a18, t(izz, 18), cond18);
            bool cond19 = ((float)t(izz, 13) >= 0.e0);
            double a19 = t(izz, 3) + t(izz, 7) * t(izz, 17);
            t(izz, 19) = ccomplex_ns::ccvmgt(a19, t(izz, 19), cond19);
            bool cond20 = ((float)t(izz, 13) >= 0.e0);
            double a20 = t(izz, 4) + t(izz, 8) * t(izz, 17);
            t(izz, 20) = ccomplex_ns::ccvmgt(a20, t(izz, 20), cond20);
            bool cond9 = ((float)t(izz, 13) >= 0.e0);
            double a9 = t(izz, 10) * t(izz, 12);
            t(izz, 9)  = ccomplex_ns::ccvmgt(a9, t(izz, 9), cond9);
        }
//
        double xgm1 = 1.e0 / (double)gm1;
        for (int izz = 1; izz <= jv; izz++) {
            t(izz, 16) = t(izz, 16) * t(izz, 16);
            t(izz, 1)  = std::pow(t(izz, 16) / t(izz, 9), xgm1);
            t(izz, 2)  = t(izz, 18);
            t(izz, 3)  = t(izz, 19);
            t(izz, 4)  = t(izz, 20);
            t(izz, 5)  = t(izz, 1) * t(izz, 16) / (double)gamma;
        }
//
    } else {
//
//     use simpler approach (DOES NOT INCLUDE MOVING BOUNDARY
//                           PROVISION AS YET
        for (int izz = 1; izz <= jv; izz++) {
            double ubar = t(izz, 6) * t(izz, 2) + t(izz, 7) * t(izz, 3) + t(izz, 8) * t(izz, 4);
//
//        outflow boundary
//
            if ((float)ubar > 0.0f) {
                t(izz, 5) = (double)p0;
            } else {
                t(izz, 1) = (double)rho0 * t(izz, 5) / (double)p0;
                t(izz, 2) = (double)u0;
                t(izz, 3) = (double)v0;
                t(izz, 4) = (double)w0;
            }
        }
//
    }
//
    return;
}

} // namespace rie1d_ns
