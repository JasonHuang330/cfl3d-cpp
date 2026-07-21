// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "amafi.h"
#include "q8vrev.h"
#include "ccomplex.h"
#include <cmath>

namespace amafi_ns {

void amafi(int& kp, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q,
           FortranArray4DRef<double> ai,
           FortranArray4DRef<double> bi,
           FortranArray4DRef<double> ci,
           FortranArray3DRef<double> dtj,
           FortranArray2DRef<double> t,
           int& nvt,
           FortranArray4DRef<double> dfp,
           FortranArray4DRef<double> dfm,
           int& imw)
{
    // COMMON block aliases
    double gamma   = (double)cmn_fluid.gamma;
    double gm1     = (double)cmn_fluid.gm1;
    int    ita     = (int)cmn_unst.ita;
    double cprec   = (double)cmn_precond.cprec;
    double uref    = (double)cmn_precond.uref;
    double avn     = (double)cmn_precond.avn;

    // matrix assembly - interior points
    int idim1 = idim - 1;
    int jdim1 = jdim - 1;
    int n     = npl * jdim1 * idim1;

    double tfacp1;
    if (std::abs(ita) == 1) {
        tfacp1 = 1.e0;
    } else {
        tfacp1 = 1.5e0;
    }

    if (imw == 0) {

        for (int m = 1; m <= 5; m++) {
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= n; izz++) {
                bi(izz, 1, m, l) = (dfp(izz, 2, m, l) - dfm(izz, 1, m, l));
                ai(izz, 1, m, l) = -dfp(izz, 1, m, l);
                ci(izz, 1, m, l) =  dfm(izz, 2, m, l);
            }
        }
        }

        // assemble matrix equation - time terms
        if (cprec == 0.) {
            for (int kpl = 1; kpl <= npl; kpl++) {
                int kk  = kp + kpl - 1;
                int jiv = (kpl - 1) * jdim1 + 1;
                for (int i = 1; i <= idim1; i++) {
                    int ji1 = jiv + (i - 1) * jdim1 * npl;
                    for (int izz = 1; izz <= jdim1; izz++) {
                        t(izz + ji1 - 1, 1) = q(izz, kk, i, 1);
                        t(izz + ji1 - 1, 2) = q(izz, kk, i, 2);
                        t(izz + ji1 - 1, 3) = q(izz, kk, i, 3);
                        t(izz + ji1 - 1, 4) = q(izz, kk, i, 4);
                        t(izz + ji1 - 1, 6) = tfacp1 * dtj(izz, kk, i);
                    }
                }
            }
        } else {
            for (int kpl = 1; kpl <= npl; kpl++) {
                int kk  = kp + kpl - 1;
                int jiv = (kpl - 1) * jdim1 + 1;
                for (int i = 1; i <= idim1; i++) {
                    int ji1 = jiv + (i - 1) * jdim1 * npl;
                    for (int izz = 1; izz <= jdim1; izz++) {
                        t(izz + ji1 - 1, 1) = q(izz, kk, i, 1);
                        t(izz + ji1 - 1, 2) = q(izz, kk, i, 2);
                        t(izz + ji1 - 1, 3) = q(izz, kk, i, 3);
                        t(izz + ji1 - 1, 4) = q(izz, kk, i, 4);
                        t(izz + ji1 - 1, 5) = q(izz, kk, i, 5);
                        t(izz + ji1 - 1, 6) = tfacp1 * dtj(izz, kk, i);
                    }
                }
            }
        }

        if (cprec == 0.) {
            for (int izz = 1; izz <= n; izz++) {
                double temp       = t(izz, 6) * t(izz, 1);
                bi(izz, 1, 1, 1) = bi(izz, 1, 1, 1) + t(izz, 6);
                bi(izz, 1, 2, 1) = bi(izz, 1, 2, 1) + t(izz, 6) * t(izz, 2);
                bi(izz, 1, 2, 2) = bi(izz, 1, 2, 2) + temp;
                bi(izz, 1, 3, 1) = bi(izz, 1, 3, 1) + t(izz, 6) * t(izz, 3);
                bi(izz, 1, 3, 3) = bi(izz, 1, 3, 3) + temp;
                bi(izz, 1, 4, 1) = bi(izz, 1, 4, 1) + t(izz, 6) * t(izz, 4);
                bi(izz, 1, 4, 4) = bi(izz, 1, 4, 4) + temp;
                bi(izz, 1, 5, 1) = bi(izz, 1, 5, 1) +
                                   t(izz, 6) * 0.5 * (t(izz, 2) * t(izz, 2) +
                                                       t(izz, 3) * t(izz, 3) +
                                                       t(izz, 4) * t(izz, 4));
                bi(izz, 1, 5, 2) = bi(izz, 1, 5, 2) + temp * t(izz, 2);
                bi(izz, 1, 5, 3) = bi(izz, 1, 5, 3) + temp * t(izz, 3);
                bi(izz, 1, 5, 4) = bi(izz, 1, 5, 4) + temp * t(izz, 4);
                bi(izz, 1, 5, 5) = bi(izz, 1, 5, 5) + t(izz, 6) / gm1;
            }
        } else {
            for (int izz = 1; izz <= n; izz++) {
                double c2    = gamma * t(izz, 5) / t(izz, 1);
                double c     = std::sqrt(c2);
                double ekin  = 0.5 * (t(izz, 2) * t(izz, 2) +
                                      t(izz, 3) * t(izz, 3) +
                                      t(izz, 4) * t(izz, 4));
                double ho    = c2 / gm1 + ekin;
                double vmag1 = 2.0 * ekin;
                double vel2_a = vmag1;
                double vel2_b = avn * uref * uref;
                double vel2  = ccomplex_ns::ccmax(vel2_a, vel2_b);
                double c2_v  = c2;
                double vel2c = ccomplex_ns::ccmin(c2_v, vel2);
                double vel   = std::sqrt(vel2c);
                vel          = cprec * vel + (1. - cprec) * c;
                double thet  = (1.0 / (vel * vel) - 1.0 / c2);
                double temp  = t(izz, 6) * t(izz, 1);
                bi(izz, 1, 1, 1) = bi(izz, 1, 1, 1) + t(izz, 6);
                bi(izz, 1, 1, 5) = bi(izz, 1, 1, 5) + t(izz, 6) * thet;
                bi(izz, 1, 2, 1) = bi(izz, 1, 2, 1) + t(izz, 6) * t(izz, 2);
                bi(izz, 1, 2, 2) = bi(izz, 1, 2, 2) + temp;
                bi(izz, 1, 2, 5) = bi(izz, 1, 2, 5) + t(izz, 6) * thet * t(izz, 2);
                bi(izz, 1, 3, 1) = bi(izz, 1, 3, 1) + t(izz, 6) * t(izz, 3);
                bi(izz, 1, 3, 3) = bi(izz, 1, 3, 3) + temp;
                bi(izz, 1, 3, 5) = bi(izz, 1, 3, 5) + t(izz, 6) * thet * t(izz, 3);
                bi(izz, 1, 4, 1) = bi(izz, 1, 4, 1) + t(izz, 6) * t(izz, 4);
                bi(izz, 1, 4, 4) = bi(izz, 1, 4, 4) + temp;
                bi(izz, 1, 4, 5) = bi(izz, 1, 4, 5) + t(izz, 6) * thet * t(izz, 4);
                bi(izz, 1, 5, 1) = bi(izz, 1, 5, 1) + t(izz, 6) * ekin;
                bi(izz, 1, 5, 2) = bi(izz, 1, 5, 2) + temp * t(izz, 2);
                bi(izz, 1, 5, 3) = bi(izz, 1, 5, 3) + temp * t(izz, 3);
                bi(izz, 1, 5, 4) = bi(izz, 1, 5, 4) + temp * t(izz, 4);
                bi(izz, 1, 5, 5) = bi(izz, 1, 5, 5) + t(izz, 6) * (1. / gm1 + thet * ho);
            }
        }

    } else {

        int jdh = jdim1 / 2;
        int mm  = jdh;
        int idh = idim1 * 2;

        if (cprec == 0.) {
            for (int kpl = 1; kpl <= npl; kpl++) {
                int kk  = kp + kpl - 1;
                int jv1 = (kpl - 1) * jdh + 1;
                int jv2 = (kpl - 1) * jdim1 + 1;
                for (int i = 1; i <= idim1; i++) {
                    int iq = i;
                    for (int k = 1; k <= 5; k++) {
                    for (int l = 1; l <= 5; l++) {
                        for (int izz = 1; izz <= mm; izz++) {
                            bi(izz + jv1 - 1, iq, k, l) = (dfp(izz + jv2 - 1, i + 1, k, l)
                                                           - dfm(izz + jv2 - 1, i,     k, l));
                            ai(izz + jv1 - 1, iq, k, l) = -dfp(izz + jv2 - 1, i,     k, l);
                            ci(izz + jv1 - 1, iq, k, l) =  dfm(izz + jv2 - 1, i + 1, k, l);
                        }
                    }
                    }
                    int jv3 = jv1 + (iq - 1) * jdh * npl;
                    for (int izz = 1; izz <= mm; izz++) {
                        t(izz + jv3 - 1, 1) = q(izz, kk, i, 1);
                        t(izz + jv3 - 1, 2) = q(izz, kk, i, 2);
                        t(izz + jv3 - 1, 3) = q(izz, kk, i, 3);
                        t(izz + jv3 - 1, 4) = q(izz, kk, i, 4);
                        t(izz + jv3 - 1, 6) = dtj(izz, kk, i);
                    }
                }
            }
        } else {
            for (int kpl = 1; kpl <= npl; kpl++) {
                int kk  = kp + kpl - 1;
                int jv1 = (kpl - 1) * jdh + 1;
                int jv2 = (kpl - 1) * jdim1 + 1;
                for (int i = 1; i <= idim1; i++) {
                    int iq = i;
                    for (int k = 1; k <= 5; k++) {
                    for (int l = 1; l <= 5; l++) {
                        for (int izz = 1; izz <= mm; izz++) {
                            bi(izz + jv1 - 1, iq, k, l) = (dfp(izz + jv2 - 1, i + 1, k, l)
                                                           - dfm(izz + jv2 - 1, i,     k, l));
                            ai(izz + jv1 - 1, iq, k, l) = -dfp(izz + jv2 - 1, i,     k, l);
                            ci(izz + jv1 - 1, iq, k, l) =  dfm(izz + jv2 - 1, i + 1, k, l);
                        }
                    }
                    }
                    int jv3 = jv1 + (iq - 1) * jdh * npl;
                    for (int izz = 1; izz <= mm; izz++) {
                        t(izz + jv3 - 1, 1) = q(izz, kk, i, 1);
                        t(izz + jv3 - 1, 2) = q(izz, kk, i, 2);
                        t(izz + jv3 - 1, 3) = q(izz, kk, i, 3);
                        t(izz + jv3 - 1, 4) = q(izz, kk, i, 4);
                        t(izz + jv3 - 1, 5) = q(izz, kk, i, 5);
                        t(izz + jv3 - 1, 6) = dtj(izz, kk, i);
                    }
                }
            }
        }

        for (int kpl = 1; kpl <= npl; kpl++) {
            int kk  = kp + kpl - 1;
            int jv1 = (kpl - 1) * jdh + 1;
            int jv2 = (kpl - 1) * jdim1 + 1 + jdh;
            for (int i = 1; i <= idim1; i++) {
                int iq = idh + 1 - i;
                for (int k = 1; k <= 5; k++) {
                for (int l = 1; l <= 5; l++) {
                    for (int izz = 1; izz <= mm; izz++) {
                        t(izz, 7) = (dfp(izz + jv2 - 1, i + 1, k, l)
                                   - dfm(izz + jv2 - 1, i,     k, l));
                    }
                    FortranArray1DRef<double> t_col7(&t(1, 7), mm, 1);
                    FortranArray1DRef<double> bi_slice(&bi(jv1, iq, k, l), mm, 1);
                    q8vrev_ns::q8vrev(mm, t_col7, mm, bi_slice);

                    FortranArray1DRef<double> dfp_slice(&dfp(jv2, i, k, l), mm, 1);
                    FortranArray1DRef<double> ci_slice(&ci(jv1, iq, k, l), mm, 1);
                    q8vrev_ns::q8vrev(mm, dfp_slice, mm, ci_slice);

                    for (int izz = 1; izz <= mm; izz++) {
                        ci(izz + jv1 - 1, iq, k, l) = -ci(izz + jv1 - 1, iq, k, l);
                    }

                    FortranArray1DRef<double> dfm_slice(&dfm(jv2, i + 1, k, l), mm, 1);
                    FortranArray1DRef<double> ai_slice(&ai(jv1, iq, k, l), mm, 1);
                    q8vrev_ns::q8vrev(mm, dfm_slice, mm, ai_slice);
                }
                }
                int jv3 = jv1 + (iq - 1) * jdh * npl;
                {
                    FortranArray1DRef<double> q_slice1(&q(jdh + 1, kk, i, 1), mm, 1);
                    FortranArray1DRef<double> t_jv3_1(&t(jv3, 1), mm, 1);
                    q8vrev_ns::q8vrev(mm, q_slice1, mm, t_jv3_1);
                }
                {
                    FortranArray1DRef<double> q_slice2(&q(jdh + 1, kk, i, 2), mm, 1);
                    FortranArray1DRef<double> t_jv3_2(&t(jv3, 2), mm, 1);
                    q8vrev_ns::q8vrev(mm, q_slice2, mm, t_jv3_2);
                }
                {
                    FortranArray1DRef<double> q_slice3(&q(jdh + 1, kk, i, 3), mm, 1);
                    FortranArray1DRef<double> t_jv3_3(&t(jv3, 3), mm, 1);
                    q8vrev_ns::q8vrev(mm, q_slice3, mm, t_jv3_3);
                }
                {
                    FortranArray1DRef<double> q_slice4(&q(jdh + 1, kk, i, 4), mm, 1);
                    FortranArray1DRef<double> t_jv3_4(&t(jv3, 4), mm, 1);
                    q8vrev_ns::q8vrev(mm, q_slice4, mm, t_jv3_4);
                }
                {
                    FortranArray1DRef<double> dtj_slice(&dtj(jdh + 1, kk, i), mm, 1);
                    FortranArray1DRef<double> t_jv3_6(&t(jv3, 6), mm, 1);
                    q8vrev_ns::q8vrev(mm, dtj_slice, mm, t_jv3_6);
                }
            }
        }

        // assemble matrix equation - time terms
        if (cprec == 0.) {
            for (int izz = 1; izz <= n; izz++) {
                t(izz, 6)        = tfacp1 * t(izz, 6);
                double temp      = t(izz, 6) * t(izz, 1);
                bi(izz, 1, 1, 1) = bi(izz, 1, 1, 1) + t(izz, 6);
                bi(izz, 1, 2, 1) = bi(izz, 1, 2, 1) + t(izz, 6) * t(izz, 2);
                bi(izz, 1, 2, 2) = bi(izz, 1, 2, 2) + temp;
                bi(izz, 1, 3, 1) = bi(izz, 1, 3, 1) + t(izz, 6) * t(izz, 3);
                bi(izz, 1, 3, 3) = bi(izz, 1, 3, 3) + temp;
                bi(izz, 1, 4, 1) = bi(izz, 1, 4, 1) + t(izz, 6) * t(izz, 4);
                bi(izz, 1, 4, 4) = bi(izz, 1, 4, 4) + temp;
                bi(izz, 1, 5, 1) = bi(izz, 1, 5, 1)
                                 + t(izz, 6) * 0.5 * (t(izz, 2) * t(izz, 2) +
                                                       t(izz, 3) * t(izz, 3) +
                                                       t(izz, 4) * t(izz, 4));
                bi(izz, 1, 5, 2) = bi(izz, 1, 5, 2) + temp * t(izz, 2);
                bi(izz, 1, 5, 3) = bi(izz, 1, 5, 3) + temp * t(izz, 3);
                bi(izz, 1, 5, 4) = bi(izz, 1, 5, 4) + temp * t(izz, 4);
                bi(izz, 1, 5, 5) = bi(izz, 1, 5, 5) + t(izz, 6) / gm1;
            }
        } else {
            for (int izz = 1; izz <= n; izz++) {
                double c2    = gamma * t(izz, 5) / t(izz, 1);
                double c     = std::sqrt(c2);
                double ekin  = 0.5 * (t(izz, 2) * t(izz, 2) +
                                      t(izz, 3) * t(izz, 3) +
                                      t(izz, 4) * t(izz, 4));
                double ho    = c2 / gm1 + ekin;
                double vmag1 = 2.0 * ekin;
                double vel2_a = vmag1;
                double vel2_b = avn * uref * uref;
                double vel2  = ccomplex_ns::ccmax(vel2_a, vel2_b);
                double c2_v  = c2;
                double vel2c = ccomplex_ns::ccmin(c2_v, vel2);
                double vel   = std::sqrt(vel2c);
                vel          = cprec * vel + (1. - cprec) * c;
                double thet  = (1.0 / (vel * vel) - 1.0 / c2);
                t(izz, 6)        = tfacp1 * t(izz, 6);
                double temp      = t(izz, 6) * t(izz, 1);
                bi(izz, 1, 1, 1) = bi(izz, 1, 1, 1) + t(izz, 6);
                bi(izz, 1, 1, 5) = bi(izz, 1, 1, 5) + t(izz, 6) * thet;
                bi(izz, 1, 2, 1) = bi(izz, 1, 2, 1) + t(izz, 6) * t(izz, 2);
                bi(izz, 1, 2, 2) = bi(izz, 1, 2, 2) + temp;
                bi(izz, 1, 2, 5) = bi(izz, 1, 2, 5) + t(izz, 6) * thet * t(izz, 2);
                bi(izz, 1, 3, 1) = bi(izz, 1, 3, 1) + t(izz, 6) * t(izz, 3);
                bi(izz, 1, 3, 3) = bi(izz, 1, 3, 3) + temp;
                bi(izz, 1, 3, 5) = bi(izz, 1, 3, 5) + t(izz, 6) * thet * t(izz, 3);
                bi(izz, 1, 4, 1) = bi(izz, 1, 4, 1) + t(izz, 6) * t(izz, 4);
                bi(izz, 1, 4, 4) = bi(izz, 1, 4, 4) + temp;
                bi(izz, 1, 4, 5) = bi(izz, 1, 4, 5) + t(izz, 6) * thet * t(izz, 4);
                bi(izz, 1, 5, 1) = bi(izz, 1, 5, 1) + t(izz, 6) * ekin;
                bi(izz, 1, 5, 2) = bi(izz, 1, 5, 2) + temp * t(izz, 2);
                bi(izz, 1, 5, 3) = bi(izz, 1, 5, 3) + temp * t(izz, 3);
                bi(izz, 1, 5, 4) = bi(izz, 1, 5, 4) + temp * t(izz, 4);
                bi(izz, 1, 5, 5) = bi(izz, 1, 5, 5) + t(izz, 6) * (1. / gm1 + thet * ho);
            }
        }

    } // end imw != 0

    return;
}

} // namespace amafi_ns
