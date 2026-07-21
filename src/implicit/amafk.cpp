// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "amafk.h"
#include "ccomplex.h"
#include <cmath>

namespace amafk_ns {

void amafk(int& i, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q, FortranArray4DRef<double> ak,
           FortranArray4DRef<double> bk, FortranArray4DRef<double> ck,
           FortranArray3DRef<double> dtj, FortranArray2DRef<double> t,
           int& nvt, FortranArray4DRef<double> dhp, FortranArray4DRef<double> dhm)
{
    // COMMON block aliases
    float& gamma   = cmn_fluid.gamma;
    float& gm1     = cmn_fluid.gm1;
    int32_t& ita   = cmn_unst.ita;
    float& cprec   = cmn_precond.cprec;
    float& uref    = cmn_precond.uref;
    float& avn     = cmn_precond.avn;

    // matrix assembly - interior points
    int kdim1 = kdim - 1;
    int jdim1 = jdim - 1;
    int n     = npl * jdim1 * kdim1;

    double tfacp1;
    if (std::abs(ita) == 1) {
        tfacp1 = 1.e0;
    } else {
        tfacp1 = 1.5e0;
    }

    for (int k = 1; k <= 5; k++) {
        for (int l = 1; l <= 5; l++) {
            // cdir$ ivdep
            for (int izz = 1; izz <= n; izz++) {
                bk(izz, 1, k, l) = (dhp(izz, 2, k, l) - dhm(izz, 1, k, l));
                ak(izz, 1, k, l) = -dhp(izz, 1, k, l);
                ck(izz, 1, k, l) =  dhm(izz, 2, k, l);
            }
        }
    }

    // assemble matrix equation - time terms
    if ((float)cprec == 0.) {
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii  = i + ipl - 1;
            int jkv = (ipl - 1) * jdim1 + 1;
            for (int k = 1; k <= kdim1; k++) {
                int jk1 = jkv + (k - 1) * jdim1 * npl;
                // cdir$ ivdep
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + jk1 - 1, 1) = q(izz, k, ii, 1);
                    t(izz + jk1 - 1, 2) = q(izz, k, ii, 2);
                    t(izz + jk1 - 1, 3) = q(izz, k, ii, 3);
                    t(izz + jk1 - 1, 4) = q(izz, k, ii, 4);
                    t(izz + jk1 - 1, 6) = tfacp1 * dtj(izz, k, ii);
                }
            }
        }
    } else {
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii  = i + ipl - 1;
            int jkv = (ipl - 1) * jdim1 + 1;
            for (int k = 1; k <= kdim1; k++) {
                int jk1 = jkv + (k - 1) * jdim1 * npl;
                // cdir$ ivdep
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + jk1 - 1, 1) = q(izz, k, ii, 1);
                    t(izz + jk1 - 1, 2) = q(izz, k, ii, 2);
                    t(izz + jk1 - 1, 3) = q(izz, k, ii, 3);
                    t(izz + jk1 - 1, 4) = q(izz, k, ii, 4);
                    t(izz + jk1 - 1, 5) = q(izz, k, ii, 5);
                    t(izz + jk1 - 1, 6) = tfacp1 * dtj(izz, k, ii);
                }
            }
        }
    }

    if ((float)cprec == 0.) {
        // cdir$ ivdep
        for (int izz = 1; izz <= n; izz++) {
            double temp          = t(izz, 6) * t(izz, 1);
            bk(izz, 1, 1, 1) = bk(izz, 1, 1, 1) + t(izz, 6);
            bk(izz, 1, 2, 1) = bk(izz, 1, 2, 1) + t(izz, 6) * t(izz, 2);
            bk(izz, 1, 2, 2) = bk(izz, 1, 2, 2) + temp;
            bk(izz, 1, 3, 1) = bk(izz, 1, 3, 1) + t(izz, 6) * t(izz, 3);
            bk(izz, 1, 3, 3) = bk(izz, 1, 3, 3) + temp;
            bk(izz, 1, 4, 1) = bk(izz, 1, 4, 1) + t(izz, 6) * t(izz, 4);
            bk(izz, 1, 4, 4) = bk(izz, 1, 4, 4) + temp;
            bk(izz, 1, 5, 1) = bk(izz, 1, 5, 1)
                              + t(izz, 6) * .5 * (t(izz, 2) * t(izz, 2) +
                                                   t(izz, 3) * t(izz, 3) +
                                                   t(izz, 4) * t(izz, 4));
            bk(izz, 1, 5, 2) = bk(izz, 1, 5, 2) + temp * t(izz, 2);
            bk(izz, 1, 5, 3) = bk(izz, 1, 5, 3) + temp * t(izz, 3);
            bk(izz, 1, 5, 4) = bk(izz, 1, 5, 4) + temp * t(izz, 4);
            bk(izz, 1, 5, 5) = bk(izz, 1, 5, 5) + t(izz, 6) / (double)gm1;
        }
    } else {
        // cdir$ ivdep
        for (int izz = 1; izz <= n; izz++) {
            double c2   = (double)gamma * t(izz, 5) / t(izz, 1);
            double c    = std::sqrt(c2);
            double ekin = 0.5 * (t(izz, 2) * t(izz, 2) +
                                  t(izz, 3) * t(izz, 3) +
                                  t(izz, 4) * t(izz, 4));
            double ho   = c2 / (double)gm1 + ekin;
            double vmag1 = 2.0 * ekin;
            double avn_uref2 = (double)avn * (double)uref * (double)uref;
            double vel2 = ccomplex_ns::ccmax(vmag1, avn_uref2);
            double c2_copy = c2;
            double vel  = std::sqrt(ccomplex_ns::ccmin(c2_copy, vel2));
            vel  = (double)cprec * vel + (1. - (double)cprec) * c;
            double thet = (1.0 / (vel * vel) - 1.0 / c2);
            double temp = t(izz, 6) * t(izz, 1);
            bk(izz, 1, 1, 1) = bk(izz, 1, 1, 1) + t(izz, 6);
            bk(izz, 1, 1, 5) = bk(izz, 1, 1, 5) + t(izz, 6) * thet;
            bk(izz, 1, 2, 1) = bk(izz, 1, 2, 1) + t(izz, 6) * t(izz, 2);
            bk(izz, 1, 2, 2) = bk(izz, 1, 2, 2) + temp;
            bk(izz, 1, 2, 5) = bk(izz, 1, 2, 5) + t(izz, 6) * thet * t(izz, 2);
            bk(izz, 1, 3, 1) = bk(izz, 1, 3, 1) + t(izz, 6) * t(izz, 3);
            bk(izz, 1, 3, 3) = bk(izz, 1, 3, 3) + temp;
            bk(izz, 1, 3, 5) = bk(izz, 1, 3, 5) + t(izz, 6) * thet * t(izz, 3);
            bk(izz, 1, 4, 1) = bk(izz, 1, 4, 1) + t(izz, 6) * t(izz, 4);
            bk(izz, 1, 4, 4) = bk(izz, 1, 4, 4) + temp;
            bk(izz, 1, 4, 5) = bk(izz, 1, 4, 5) + t(izz, 6) * thet * t(izz, 4);
            bk(izz, 1, 5, 1) = bk(izz, 1, 5, 1) + t(izz, 6) * ekin;
            bk(izz, 1, 5, 2) = bk(izz, 1, 5, 2) + temp * t(izz, 2);
            bk(izz, 1, 5, 3) = bk(izz, 1, 5, 3) + temp * t(izz, 3);
            bk(izz, 1, 5, 4) = bk(izz, 1, 5, 4) + temp * t(izz, 4);
            bk(izz, 1, 5, 5) = bk(izz, 1, 5, 5) + t(izz, 6) * (1. / (double)gm1 + thet * ho);
        }
    }

    return;
}

} // namespace amafk_ns
