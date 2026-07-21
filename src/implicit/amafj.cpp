// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "amafj.h"
#include <cmath>
#include <cstdint>

namespace amafj_ns {

// Inline helpers for ccmax/ccmin (avoid lvalue-ref issues with ccomplex_ns)
static inline double ccmax(double a, double b) { return (a > b) ? a : b; }
static inline double ccmin(double a, double b) { return (a < b) ? a : b; }

void amafj(int& i, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q,
           FortranArray4DRef<double> aj,
           FortranArray4DRef<double> bj,
           FortranArray4DRef<double> cj,
           FortranArray3DRef<double> dtj,
           FortranArray2DRef<double> t,
           int& nvt,
           FortranArray4DRef<double> dgp,
           FortranArray4DRef<double> dgm)
{
    // COMMON block aliases
    float& gamma   = cmn_fluid.gamma;
    float& gm1     = cmn_fluid.gm1;
    int32_t& ita   = cmn_unst.ita;
    float& cprec   = cmn_precond.cprec;
    float& uref    = cmn_precond.uref;
    float& avn     = cmn_precond.avn;

    // Local scalars
    int kdim1, jdim1, kv, n, n0, jj, kk, k, l, j, izz;
    int ipl, ii, jkv;
    double tfacp1;
    double temp, c2, c, ekin, ho, vmag1, vel2, vel, thet;

    // Assemble matrix equation - interior points
    kdim1 = kdim - 1;
    jdim1 = jdim - 1;
    kv    = npl * kdim1;
    n     = jdim1 * kv;

    if (std::abs(ita) == 1) {
        tfacp1 = 1.0e0;
    } else {
        tfacp1 = 1.5e0;
    }

    // Loop 1232: gather dgp/dgm into t, then compute aj/bj/cj
    for (k = 1; k <= 5; k++) {
        for (l = 1; l <= 5; l++) {
            for (j = 1; j <= jdim; j++) {
                n0 = (j - 1) * kv + 1;

                jj = 1 - jdim;
                for (kk = 1; kk <= kv; kk++) {
                    jj = jj + jdim;
                    t(n0 + kk - 1, 1) = dgp(j + jj - 1, 1, k, l);
                    t(n0 + kk - 1, 2) = dgm(j + jj - 1, 1, k, l);
                }
            }
            // cdir$ ivdep
            for (izz = 1; izz <= n; izz++) {
                bj(izz, 1, k, l) = (t(izz + kv, 1) - t(izz, 2));
                aj(izz, 1, k, l) = -t(izz, 1);
                cj(izz, 1, k, l) =  t(izz + kv, 2);
            }
        }
    }

    // Assemble matrix equation - time terms
    if ((float)cprec == 0.f) {
        for (ipl = 1; ipl <= npl; ipl++) {
            ii  = i + ipl - 1;
            jkv = (ipl - 1) * kdim1;
            for (j = 1; j <= jdim1; j++) {
                n0 = (j - 1) * kv + jkv + 1;

                jj = 1 - jdim;
                // cdir$ ivdep
                for (kk = 1; kk <= kdim1; kk++) {
                    jj = jj + jdim;
                    t(n0 + kk - 1, 1) = q(j + jj - 1, 1, ii, 1);
                    t(n0 + kk - 1, 2) = q(j + jj - 1, 1, ii, 2);
                    t(n0 + kk - 1, 3) = q(j + jj - 1, 1, ii, 3);
                    t(n0 + kk - 1, 4) = q(j + jj - 1, 1, ii, 4);
                    t(n0 + kk - 1, 6) = tfacp1 * dtj(j + jj - 1, 1, ii);
                }
            }
        }
    } else {
        for (ipl = 1; ipl <= npl; ipl++) {
            ii  = i + ipl - 1;
            jkv = (ipl - 1) * kdim1;
            for (j = 1; j <= jdim1; j++) {
                n0 = (j - 1) * kv + jkv + 1;

                jj = 1 - jdim;
                // cdir$ ivdep
                for (kk = 1; kk <= kdim1; kk++) {
                    jj = jj + jdim;
                    t(n0 + kk - 1, 1) = q(j + jj - 1, 1, ii, 1);
                    t(n0 + kk - 1, 2) = q(j + jj - 1, 1, ii, 2);
                    t(n0 + kk - 1, 3) = q(j + jj - 1, 1, ii, 3);
                    t(n0 + kk - 1, 4) = q(j + jj - 1, 1, ii, 4);
                    t(n0 + kk - 1, 5) = q(j + jj - 1, 1, ii, 5);
                    t(n0 + kk - 1, 6) = tfacp1 * dtj(j + jj - 1, 1, ii);
                }
            }
        }
    }

    if ((float)cprec == 0.f) {
        // cdir$ ivdep
        for (izz = 1; izz <= n; izz++) {
            temp              = t(izz, 6) * t(izz, 1);
            bj(izz, 1, 1, 1) = bj(izz, 1, 1, 1) + t(izz, 6);
            bj(izz, 1, 2, 1) = bj(izz, 1, 2, 1) + t(izz, 6) * t(izz, 2);
            bj(izz, 1, 2, 2) = bj(izz, 1, 2, 2) + temp;
            bj(izz, 1, 3, 1) = bj(izz, 1, 3, 1) + t(izz, 6) * t(izz, 3);
            bj(izz, 1, 3, 3) = bj(izz, 1, 3, 3) + temp;
            bj(izz, 1, 4, 1) = bj(izz, 1, 4, 1) + t(izz, 6) * t(izz, 4);
            bj(izz, 1, 4, 4) = bj(izz, 1, 4, 4) + temp;
            bj(izz, 1, 5, 1) = bj(izz, 1, 5, 1)
                              + t(izz, 6) * 0.5 * (t(izz, 2) * t(izz, 2) +
                                                    t(izz, 3) * t(izz, 3) +
                                                    t(izz, 4) * t(izz, 4));
            bj(izz, 1, 5, 2) = bj(izz, 1, 5, 2) + temp * t(izz, 2);
            bj(izz, 1, 5, 3) = bj(izz, 1, 5, 3) + temp * t(izz, 3);
            bj(izz, 1, 5, 4) = bj(izz, 1, 5, 4) + temp * t(izz, 4);
            bj(izz, 1, 5, 5) = bj(izz, 1, 5, 5) + t(izz, 6) / (double)gm1;
        }
    } else {
        // cdir$ ivdep
        for (izz = 1; izz <= n; izz++) {
            c2    = (double)gamma * t(izz, 5) / t(izz, 1);
            c     = std::sqrt(c2);
            ekin  = 0.5 * (t(izz, 2) * t(izz, 2) +
                           t(izz, 3) * t(izz, 3) +
                           t(izz, 4) * t(izz, 4));
            ho    = c2 / (double)gm1 + ekin;
            vmag1 = 2.0 * ekin;
            vel2  = ccmax(vmag1, (double)avn * (double)uref * (double)uref);
            vel   = std::sqrt(ccmin(c2, vel2));
            vel   = (double)cprec * vel + (1.0 - (double)cprec) * c;
            thet  = (1.0 / (vel * vel) - 1.0 / c2);
            temp              = t(izz, 6) * t(izz, 1);
            bj(izz, 1, 1, 1) = bj(izz, 1, 1, 1) + t(izz, 6);
            bj(izz, 1, 1, 5) = bj(izz, 1, 1, 5) + t(izz, 6) * thet;
            bj(izz, 1, 2, 1) = bj(izz, 1, 2, 1) + t(izz, 6) * t(izz, 2);
            bj(izz, 1, 2, 2) = bj(izz, 1, 2, 2) + temp;
            bj(izz, 1, 2, 5) = bj(izz, 1, 2, 5) + t(izz, 6) * thet * t(izz, 2);
            bj(izz, 1, 3, 1) = bj(izz, 1, 3, 1) + t(izz, 6) * t(izz, 3);
            bj(izz, 1, 3, 3) = bj(izz, 1, 3, 3) + temp;
            bj(izz, 1, 3, 5) = bj(izz, 1, 3, 5) + t(izz, 6) * thet * t(izz, 3);
            bj(izz, 1, 4, 1) = bj(izz, 1, 4, 1) + t(izz, 6) * t(izz, 4);
            bj(izz, 1, 4, 4) = bj(izz, 1, 4, 4) + temp;
            bj(izz, 1, 4, 5) = bj(izz, 1, 4, 5) + t(izz, 6) * thet * t(izz, 4);
            bj(izz, 1, 5, 1) = bj(izz, 1, 5, 1) + t(izz, 6) * ekin;
            bj(izz, 1, 5, 2) = bj(izz, 1, 5, 2) + temp * t(izz, 2);
            bj(izz, 1, 5, 3) = bj(izz, 1, 5, 3) + temp * t(izz, 3);
            bj(izz, 1, 5, 4) = bj(izz, 1, 5, 4) + temp * t(izz, 4);
            bj(izz, 1, 5, 5) = bj(izz, 1, 5, 5) + t(izz, 6) * (1.0 / (double)gm1 + thet * ho);
        }
    }

    return;
}

} // namespace amafj_ns
