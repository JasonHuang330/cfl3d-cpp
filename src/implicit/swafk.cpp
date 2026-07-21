// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "swafk.h"
#include "bsub.h"
#include "q8vrev.h"
#include <cmath>

namespace swafk_ns {

// Inline helpers for ccmax/ccmin (ccomplex functions on real values)
static inline double ccmax(double a, double b) { return (a > b) ? a : b; }
static inline double ccmin(double a, double b) { return (a < b) ? a : b; }

void swafk(int& i, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q, FortranArray4DRef<double> ak,
           FortranArray4DRef<double> bk, FortranArray4DRef<double> ck,
           FortranArray3DRef<double> dtj, FortranArray3DRef<double> f,
           int& nvt, FortranArray4DRef<double> res, int& imw)
{
    // COMMON block aliases
    float& gamma   = cmn_fluid.gamma;
    float& gm1     = cmn_fluid.gm1;
    int32_t& ita   = cmn_unst.ita;
    float& cprec   = cmn_precond.cprec;
    float& uref    = cmn_precond.uref;
    float& avn     = cmn_precond.avn;

    // Local scalars
    int kdim1, jdim1, n;
    double gm1i, tfacp1;
    double temp, c2, c, ekin, ho, vmag1, vel2, vel, thet, restmp;
    int ipl, ii, izz, l, k, js, jk, jdh, kq1, kq2, jv;
    int il, iu, id1, id2;

    kdim1 = kdim - 1;
    jdim1 = jdim - 1;
    n     = jdim * kdim1 - 1;
    gm1i  = 1.0 / (double)gm1;

    if (std::abs(ita) == 1) {
        tfacp1 = 1.0e0;
    } else {
        tfacp1 = 1.5e0;
    }

    if (imw == 0) {
        if ((float)cprec == 0.0f) {
            // imw==0, no preconditioning
            for (ipl = 1; ipl <= npl; ipl++) {
                ii = i + ipl - 1;
                // cdir$ ivdep
                for (izz = 1; izz <= n; izz++) {
                    temp              = tfacp1 * dtj(izz,1,ii) * q(izz,1,ii,1);
                    res(izz,1,ii,1)   = tfacp1 * dtj(izz,1,ii) * res(izz,1,ii,1);
                    res(izz,1,ii,5)   = 0.5 * res(izz,1,ii,1)
                                            * (q(izz,1,ii,2) * q(izz,1,ii,2)
                                              +q(izz,1,ii,3) * q(izz,1,ii,3)
                                              +q(izz,1,ii,4) * q(izz,1,ii,4))
                                       + temp * res(izz,1,ii,2) * q(izz,1,ii,2)
                                       + temp * res(izz,1,ii,3) * q(izz,1,ii,3)
                                       + temp * res(izz,1,ii,4) * q(izz,1,ii,4)
                                       + tfacp1 * gm1i * res(izz,1,ii,5) * dtj(izz,1,ii);
                    res(izz,1,ii,2)   = temp * res(izz,1,ii,2)
                                       + q(izz,1,ii,2) * res(izz,1,ii,1);
                    res(izz,1,ii,3)   = temp * res(izz,1,ii,3)
                                       + q(izz,1,ii,3) * res(izz,1,ii,1);
                    res(izz,1,ii,4)   = temp * res(izz,1,ii,4)
                                       + q(izz,1,ii,4) * res(izz,1,ii,1);
                }
                for (l = 1; l <= 5; l++) {
                    js = (ipl - 1) * jdim1 + 1;
                    for (k = 1; k <= kdim1; k++) {
                        jk = (k - 1) * jdim + 1;
                        // cdir$ ivdep
                        for (izz = 1; izz <= jdim1; izz++) {
                            f(izz+js-1, k, l) = res(izz+jk-1, 1, ii, l);
                        }
                    }
                }
            }
        } else {
            // imw==0, with preconditioning
            for (ipl = 1; ipl <= npl; ipl++) {
                ii = i + ipl - 1;
                // cdir$ ivdep
                for (izz = 1; izz <= n; izz++) {
                    // modifications for preconditioning
                    c2     = (double)gamma * q(izz,1,ii,5) / q(izz,1,ii,1);
                    c      = std::sqrt(c2);
                    ekin   = 0.5 * (q(izz,1,ii,2) * q(izz,1,ii,2)
                                  + q(izz,1,ii,3) * q(izz,1,ii,3)
                                  + q(izz,1,ii,4) * q(izz,1,ii,4));
                    ho     = c2 / (double)gm1 + ekin;
                    vmag1  = 2.0 * ekin;
                    vel2   = ccmax(vmag1, (double)avn * (double)uref * (double)uref);
                    vel    = std::sqrt(ccmin(c2, vel2));
                    vel    = (double)cprec * vel + (1.0 - (double)cprec) * c;
                    thet   = 1.0 / (vel * vel) - 1.0 / c2;
                    restmp = tfacp1 * dtj(izz,1,ii) * thet * res(izz,1,ii,5);

                    temp              = tfacp1 * dtj(izz,1,ii) * q(izz,1,ii,1);
                    res(izz,1,ii,1)   = tfacp1 * dtj(izz,1,ii) * res(izz,1,ii,1)
                                       + restmp;
                    res(izz,1,ii,5)   = res(izz,1,ii,1) * ekin
                                       + temp * res(izz,1,ii,2) * q(izz,1,ii,2)
                                       + temp * res(izz,1,ii,3) * q(izz,1,ii,3)
                                       + temp * res(izz,1,ii,4) * q(izz,1,ii,4)
                                       + tfacp1 * gm1i * res(izz,1,ii,5) * dtj(izz,1,ii)
                                       + restmp * ho;
                    res(izz,1,ii,2)   = temp * res(izz,1,ii,2)
                                       + q(izz,1,ii,2) * res(izz,1,ii,1)
                                       + restmp * q(izz,1,ii,2);
                    res(izz,1,ii,3)   = temp * res(izz,1,ii,3)
                                       + q(izz,1,ii,3) * res(izz,1,ii,1)
                                       + restmp * q(izz,1,ii,3);
                    res(izz,1,ii,4)   = temp * res(izz,1,ii,4)
                                       + q(izz,1,ii,4) * res(izz,1,ii,1)
                                       + restmp * q(izz,1,ii,4);
                }
                for (l = 1; l <= 5; l++) {
                    js = (ipl - 1) * jdim1 + 1;
                    for (k = 1; k <= kdim1; k++) {
                        jk = (k - 1) * jdim + 1;
                        // cdir$ ivdep
                        for (izz = 1; izz <= jdim1; izz++) {
                            f(izz+js-1, k, l) = res(izz+jk-1, 1, ii, l);
                        }
                    }
                }
            }
        }
    } else {
        // imw != 0
        if ((float)cprec == 0.0f) {
            // imw!=0, no preconditioning
            jdh = jdim1 / 2;
            for (ipl = 1; ipl <= npl; ipl++) {
                ii = i + ipl - 1;
                // cdir$ ivdep
                for (izz = 1; izz <= n; izz++) {
                    temp              = tfacp1 * dtj(izz,1,ii) * q(izz,1,ii,1);
                    res(izz,1,ii,1)   = tfacp1 * dtj(izz,1,ii) * res(izz,1,ii,1);
                    res(izz,1,ii,5)   = 0.5 * res(izz,1,ii,1)
                                            * (q(izz,1,ii,2) * q(izz,1,ii,2)
                                              +q(izz,1,ii,3) * q(izz,1,ii,3)
                                              +q(izz,1,ii,4) * q(izz,1,ii,4))
                                       + temp * res(izz,1,ii,2) * q(izz,1,ii,2)
                                       + temp * res(izz,1,ii,3) * q(izz,1,ii,3)
                                       + temp * res(izz,1,ii,4) * q(izz,1,ii,4)
                                       + tfacp1 * gm1i * res(izz,1,ii,5) * dtj(izz,1,ii);
                    res(izz,1,ii,2)   = temp * res(izz,1,ii,2)
                                       + q(izz,1,ii,2) * res(izz,1,ii,1);
                    res(izz,1,ii,3)   = temp * res(izz,1,ii,3)
                                       + q(izz,1,ii,3) * res(izz,1,ii,1);
                    res(izz,1,ii,4)   = temp * res(izz,1,ii,4)
                                       + q(izz,1,ii,4) * res(izz,1,ii,1);
                }
                for (l = 1; l <= 5; l++) {
                    for (k = 1; k <= kdim1; k++) {
                        kq1 = kdim1 + k;
                        jv  = (ipl - 1) * jdh + 1;
                        kq2 = kdim - k;
                        jk  = (k - 1) * jdim + 1;
                        // cdir$ ivdep
                        for (izz = 1; izz <= jdh; izz++) {
                            f(izz+jv-1, kq1, l) = res(izz+jk-1, 1, ii, l);
                        }
                        {
                            int jdh_loc = jdh;
                            FortranArray1DRef<double> src(&res(jk+jdh, 1, ii, l), jdh, 1);
                            FortranArray1DRef<double> dst(&f(jv, kq2, l), jdh, 1);
                            q8vrev_ns::q8vrev(jdh_loc, src, jdh_loc, dst);
                        }
                    }
                }
            }
        } else {
            // imw!=0, with preconditioning
            jdh = jdim1 / 2;
            for (ipl = 1; ipl <= npl; ipl++) {
                ii = i + ipl - 1;
                // cdir$ ivdep
                for (izz = 1; izz <= n; izz++) {
                    // modifications for preconditioning
                    c2     = (double)gamma * q(izz,1,ii,5) / q(izz,1,ii,1);
                    c      = std::sqrt(c2);
                    ekin   = 0.5 * (q(izz,1,ii,2) * q(izz,1,ii,2)
                                  + q(izz,1,ii,3) * q(izz,1,ii,3)
                                  + q(izz,1,ii,4) * q(izz,1,ii,4));
                    ho     = c2 / (double)gm1 + ekin;
                    vmag1  = 2.0 * ekin;
                    vel2   = ccmax(vmag1, (double)avn * (double)uref * (double)uref);
                    vel    = std::sqrt(ccmin(c2, vel2));
                    vel    = (double)cprec * vel + (1.0 - (double)cprec) * c;
                    thet   = 1.0 / (vel * vel) - 1.0 / c2;
                    restmp = tfacp1 * dtj(izz,1,ii) * thet * res(izz,1,ii,5);

                    temp              = tfacp1 * dtj(izz,1,ii) * q(izz,1,ii,1);
                    res(izz,1,ii,1)   = tfacp1 * dtj(izz,1,ii) * res(izz,1,ii,1)
                                       + restmp;
                    res(izz,1,ii,5)   = 0.5 * res(izz,1,ii,1)
                                            * (q(izz,1,ii,2) * q(izz,1,ii,2)
                                              +q(izz,1,ii,3) * q(izz,1,ii,3)
                                              +q(izz,1,ii,4) * q(izz,1,ii,4))
                                       + temp * res(izz,1,ii,2) * q(izz,1,ii,2)
                                       + temp * res(izz,1,ii,3) * q(izz,1,ii,3)
                                       + temp * res(izz,1,ii,4) * q(izz,1,ii,4)
                                       + tfacp1 * gm1i * res(izz,1,ii,5) * dtj(izz,1,ii)
                                       + restmp * ho;
                    res(izz,1,ii,2)   = temp * res(izz,1,ii,2)
                                       + q(izz,1,ii,2) * res(izz,1,ii,1)
                                       + restmp * q(izz,1,ii,2);
                    res(izz,1,ii,3)   = temp * res(izz,1,ii,3)
                                       + q(izz,1,ii,3) * res(izz,1,ii,1)
                                       + restmp * q(izz,1,ii,3);
                    res(izz,1,ii,4)   = temp * res(izz,1,ii,4)
                                       + q(izz,1,ii,4) * res(izz,1,ii,1)
                                       + restmp * q(izz,1,ii,4);
                }
                for (l = 1; l <= 5; l++) {
                    for (k = 1; k <= kdim1; k++) {
                        kq1 = kdim1 + k;
                        jv  = (ipl - 1) * jdh + 1;
                        kq2 = kdim - k;
                        jk  = (k - 1) * jdim + 1;
                        // cdir$ ivdep
                        for (izz = 1; izz <= jdh; izz++) {
                            f(izz+jv-1, kq1, l) = res(izz+jk-1, 1, ii, l);
                        }
                        {
                            int jdh_loc = jdh;
                            FortranArray1DRef<double> src(&res(jk+jdh, 1, ii, l), jdh, 1);
                            FortranArray1DRef<double> dst(&f(jv, kq2, l), jdh, 1);
                            q8vrev_ns::q8vrev(jdh_loc, src, jdh_loc, dst);
                        }
                    }
                }
            }
        }
    }

    // solve matrix equation
    il  = 1;
    iu  = kdim1 * (imw + 1);
    n   = npl * jdim1 / (imw + 1);

    id1 = npl * (jdim - 1) / (imw + 1);
    id2 = (kdim - 1) * (imw + 1);
    {
        int i1 = 1;
        bsub_ns::bsub(id1, id2, ak, bk, ck, f, i1, n, il, iu);
    }

    // update delta q
    if (imw == 0) {
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            jv = (ipl - 1) * jdim1 + 1;
            for (k = 1; k <= kdim1; k++) {
                for (l = 1; l <= 5; l++) {
                    // cdir$ ivdep
                    for (izz = 1; izz <= jdim1; izz++) {
                        res(izz, k, ii, l) = f(izz+jv-1, k, l);
                    }
                }
            }
        }
    } else {
        for (ipl = 1; ipl <= npl; ipl++) {
            ii  = i + ipl - 1;
            jv  = (ipl - 1) * jdh + 1;
            for (k = 1; k <= kdim1; k++) {
                kq1 = kdim1 + k;
                kq2 = kdim - k;
                for (l = 1; l <= 5; l++) {
                    // cdir$ ivdep
                    for (izz = 1; izz <= jdh; izz++) {
                        res(izz, k, ii, l) = f(izz+jv-1, kq1, l);
                    }
                    {
                        int jdh_loc = jdh;
                        FortranArray1DRef<double> src(&f(jv, kq2, l), jdh, 1);
                        FortranArray1DRef<double> dst(&res(1+jdh, k, ii, l), jdh, 1);
                        q8vrev_ns::q8vrev(jdh_loc, src, jdh_loc, dst);
                    }
                }
            }
        }
    }

    // label 3000 continue / return
}

} // namespace swafk_ns
