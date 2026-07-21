// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "swafi.h"
#include "q8vrev.h"
#include "bsub.h"
#include "ccomplex.h"
#include <cmath>
#include <cstdlib>

namespace swafi_ns {

void swafi(int& kp, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q, FortranArray4DRef<double> ai,
           FortranArray4DRef<double> bi, FortranArray4DRef<double> ci,
           FortranArray3DRef<double> dtj, FortranArray3DRef<double> f,
           int& nvt, FortranArray4DRef<double> res, int& imw)
{
    // COMMON block aliases
    float& cprec = cmn_precond.cprec;
    float& uref  = cmn_precond.uref;
    float& avn   = cmn_precond.avn;
    float& gamma = cmn_fluid.gamma;
    float& gm1   = cmn_fluid.gm1;
    int32_t& ita = cmn_unst.ita;

    // Local scalars
    int idim1 = idim - 1;
    int jdim1 = jdim - 1;
    double gm1i = 1.0 / (double)gm1;
    double tfacp1;
    if (std::abs(ita) == 1) {
        tfacp1 = 1.0e0;
    } else {
        tfacp1 = 1.5e0;
    }

    if (imw == 0) {
        if ((float)cprec == 0.f) {
            // Loop 1000
            for (int kpl = 1; kpl <= npl; kpl++) {
                int kk = kp + kpl - 1;
                int js = (kpl - 1) * jdim1 + 1;
                for (int i = 1; i <= idim1; i++) {
                    for (int izz = 1; izz <= jdim1; izz++) {
                        double temp = tfacp1 * dtj(izz, kk, i) * q(izz, kk, i, 1);
                        f(izz + js - 1, i, 1) = tfacp1 * dtj(izz, kk, i) * res(izz, kk, i, 1);
                        f(izz + js - 1, i, 5) = 0.5 * f(izz + js - 1, i, 1)
                            * (q(izz, kk, i, 2) * q(izz, kk, i, 2)
                              + q(izz, kk, i, 3) * q(izz, kk, i, 3)
                              + q(izz, kk, i, 4) * q(izz, kk, i, 4))
                            + temp * res(izz, kk, i, 2) * q(izz, kk, i, 2)
                            + temp * res(izz, kk, i, 3) * q(izz, kk, i, 3)
                            + temp * res(izz, kk, i, 4) * q(izz, kk, i, 4)
                            + tfacp1 * gm1i * res(izz, kk, i, 5) * dtj(izz, kk, i);
                        f(izz + js - 1, i, 2) = temp * res(izz, kk, i, 2)
                            + q(izz, kk, i, 2) * f(izz + js - 1, i, 1);
                        f(izz + js - 1, i, 3) = temp * res(izz, kk, i, 3)
                            + q(izz, kk, i, 3) * f(izz + js - 1, i, 1);
                        f(izz + js - 1, i, 4) = temp * res(izz, kk, i, 4)
                            + q(izz, kk, i, 4) * f(izz + js - 1, i, 1);
                    }
                }
            }
        } else {
            // Loop 10001 - preconditioning
            for (int kpl = 1; kpl <= npl; kpl++) {
                int kk = kp + kpl - 1;
                int js = (kpl - 1) * jdim1 + 1;
                for (int i = 1; i <= idim1; i++) {
                    for (int izz = 1; izz <= jdim1; izz++) {
                        double c2 = (double)gamma * q(izz, kk, i, 5) / q(izz, kk, i, 1);
                        double c = std::sqrt(c2);
                        double ekin = 0.5 * (q(izz, kk, i, 2) * q(izz, kk, i, 2)
                                           + q(izz, kk, i, 3) * q(izz, kk, i, 3)
                                           + q(izz, kk, i, 4) * q(izz, kk, i, 4));
                        double ho = c2 / (double)gm1 + ekin;
                        double vmag1 = 2.0 * ekin;
                        double avn_uref2 = (double)avn * (double)uref * (double)uref;
                        double vel2 = ccomplex_ns::ccmax(vmag1, avn_uref2);
                        double vel = std::sqrt(ccomplex_ns::ccmin(c2, vel2));
                        vel = (double)cprec * vel + (1.0 - (double)cprec) * c;
                        double thet = (1.0 / (vel * vel) - 1.0 / c2);
                        double restmp = tfacp1 * dtj(izz, kk, i) * thet * res(izz, kk, i, 5);

                        double temp = tfacp1 * dtj(izz, kk, i) * q(izz, kk, i, 1);
                        f(izz + js - 1, i, 1) = tfacp1 * dtj(izz, kk, i) * res(izz, kk, i, 1)
                            + restmp;
                        f(izz + js - 1, i, 5) = f(izz + js - 1, i, 1) * ekin
                            + temp * res(izz, kk, i, 2) * q(izz, kk, i, 2)
                            + temp * res(izz, kk, i, 3) * q(izz, kk, i, 3)
                            + temp * res(izz, kk, i, 4) * q(izz, kk, i, 4)
                            + tfacp1 * gm1i * res(izz, kk, i, 5) * dtj(izz, kk, i)
                            + restmp * ho;
                        f(izz + js - 1, i, 2) = temp * res(izz, kk, i, 2)
                            + q(izz, kk, i, 2) * f(izz + js - 1, i, 1)
                            + restmp * q(izz, kk, i, 2);
                        f(izz + js - 1, i, 3) = temp * res(izz, kk, i, 3)
                            + q(izz, kk, i, 3) * f(izz + js - 1, i, 1)
                            + restmp * q(izz, kk, i, 3);
                        f(izz + js - 1, i, 4) = temp * res(izz, kk, i, 4)
                            + q(izz, kk, i, 4) * f(izz + js - 1, i, 1)
                            + restmp * q(izz, kk, i, 4);
                    }
                }
            }
        }
    } else {
        // imw != 0
        if ((float)cprec == 0.f) {
            int jdh = jdim1 / 2;
            int idh = idim1 * 2;
            for (int kpl = 1; kpl <= npl; kpl++) {
                int kk = kp + kpl - 1;
                int js = (kpl - 1) * jdh + 1;
                for (int i = 1; i <= idim1; i++) {
                    int iq1 = i;
                    int iq2 = idh + 1 - i;
                    // Loop 1001
                    for (int izz = 1; izz <= jdh; izz++) {
                        double temp = tfacp1
                            * dtj(izz + jdh, kk, i) * q(izz + jdh, kk, i, 1);
                        f(izz + js - 1, iq1, 1) = tfacp1
                            * dtj(izz + jdh, kk, i) * res(izz + jdh, kk, i, 1);
                        f(izz + js - 1, iq1, 5) = 0.5 * f(izz + js - 1, iq1, 1)
                            * (q(izz + jdh, kk, i, 2) * q(izz + jdh, kk, i, 2)
                              + q(izz + jdh, kk, i, 3) * q(izz + jdh, kk, i, 3)
                              + q(izz + jdh, kk, i, 4) * q(izz + jdh, kk, i, 4))
                            + temp * res(izz + jdh, kk, i, 2) * q(izz + jdh, kk, i, 2)
                            + temp * res(izz + jdh, kk, i, 3) * q(izz + jdh, kk, i, 3)
                            + temp * res(izz + jdh, kk, i, 4) * q(izz + jdh, kk, i, 4)
                            + tfacp1 * gm1i * res(izz + jdh, kk, i, 5) * dtj(izz + jdh, kk, i);
                        f(izz + js - 1, iq1, 2) = temp * res(izz + jdh, kk, i, 2)
                            + q(izz + jdh, kk, i, 2) * f(izz + js - 1, iq1, 1);
                        f(izz + js - 1, iq1, 3) = temp * res(izz + jdh, kk, i, 3)
                            + q(izz + jdh, kk, i, 3) * f(izz + js - 1, iq1, 1);
                        f(izz + js - 1, iq1, 4) = temp * res(izz + jdh, kk, i, 4)
                            + q(izz + jdh, kk, i, 4) * f(izz + js - 1, iq1, 1);
                    }
                    // q8vrev calls
                    {
                        FortranArray1DRef<double> fa1(&f(js, iq1, 1), jdh);
                        FortranArray1DRef<double> fb1(&f(js, iq2, 1), jdh);
                        q8vrev_ns::q8vrev(jdh, fa1, jdh, fb1);
                    }
                    {
                        FortranArray1DRef<double> fa2(&f(js, iq1, 2), jdh);
                        FortranArray1DRef<double> fb2(&f(js, iq2, 2), jdh);
                        q8vrev_ns::q8vrev(jdh, fa2, jdh, fb2);
                    }
                    {
                        FortranArray1DRef<double> fa3(&f(js, iq1, 3), jdh);
                        FortranArray1DRef<double> fb3(&f(js, iq2, 3), jdh);
                        q8vrev_ns::q8vrev(jdh, fa3, jdh, fb3);
                    }
                    {
                        FortranArray1DRef<double> fa4(&f(js, iq1, 4), jdh);
                        FortranArray1DRef<double> fb4(&f(js, iq2, 4), jdh);
                        q8vrev_ns::q8vrev(jdh, fa4, jdh, fb4);
                    }
                    {
                        FortranArray1DRef<double> fa5(&f(js, iq1, 5), jdh);
                        FortranArray1DRef<double> fb5(&f(js, iq2, 5), jdh);
                        q8vrev_ns::q8vrev(jdh, fa5, jdh, fb5);
                    }
                    // Loop 1002
                    for (int izz = 1; izz <= jdh; izz++) {
                        double temp = tfacp1 * dtj(izz, kk, i) * q(izz, kk, i, 1);
                        f(izz + js - 1, iq1, 1) = tfacp1 * dtj(izz, kk, i) * res(izz, kk, i, 1);
                        f(izz + js - 1, iq1, 5) = 0.5 * f(izz + js - 1, iq1, 1)
                            * (q(izz, kk, i, 2) * q(izz, kk, i, 2)
                              + q(izz, kk, i, 3) * q(izz, kk, i, 3)
                              + q(izz, kk, i, 4) * q(izz, kk, i, 4))
                            + temp * res(izz, kk, i, 2) * q(izz, kk, i, 2)
                            + temp * res(izz, kk, i, 3) * q(izz, kk, i, 3)
                            + temp * res(izz, kk, i, 4) * q(izz, kk, i, 4)
                            + tfacp1 * gm1i * res(izz, kk, i, 5) * dtj(izz, kk, i);
                        f(izz + js - 1, iq1, 2) = temp * res(izz, kk, i, 2)
                            + q(izz, kk, i, 2) * f(izz + js - 1, iq1, 1);
                        f(izz + js - 1, iq1, 3) = temp * res(izz, kk, i, 3)
                            + q(izz, kk, i, 3) * f(izz + js - 1, iq1, 1);
                        f(izz + js - 1, iq1, 4) = temp * res(izz, kk, i, 4)
                            + q(izz, kk, i, 4) * f(izz + js - 1, iq1, 1);
                    }
                }
            }
        } else {
            // imw != 0, cprec > 0
            int jdh = jdim1 / 2;
            int idh = idim1 * 2;
            for (int kpl = 1; kpl <= npl; kpl++) {
                int kk = kp + kpl - 1;
                int js = (kpl - 1) * jdh + 1;
                for (int i = 1; i <= idim1; i++) {
                    int iq1 = i;
                    int iq2 = idh + 1 - i;
                    // Loop 10011
                    for (int izz = 1; izz <= jdh; izz++) {
                        double c2 = (double)gamma * q(izz + jdh, kk, i, 5) / q(izz + jdh, kk, i, 1);
                        double c = std::sqrt(c2);
                        double ekin = 0.5 * (q(izz + jdh, kk, i, 2) * q(izz + jdh, kk, i, 2)
                                           + q(izz + jdh, kk, i, 3) * q(izz + jdh, kk, i, 3)
                                           + q(izz + jdh, kk, i, 4) * q(izz + jdh, kk, i, 4));
                        double ho = c2 / (double)gm1 + ekin;
                        double vmag1 = 2.0 * ekin;
                        double avn_uref2 = (double)avn * (double)uref * (double)uref;
                        double vel2 = ccomplex_ns::ccmax(vmag1, avn_uref2);
                        double vel = std::sqrt(ccomplex_ns::ccmin(c2, vel2));
                        vel = (double)cprec * vel + (1.0 - (double)cprec) * c;
                        double thet = (1.0 / (vel * vel) - 1.0 / c2);
                        double restmp = tfacp1 * dtj(izz + jdh, kk, i) * thet * res(izz + jdh, kk, i, 5);

                        double temp = tfacp1 * dtj(izz + jdh, kk, i) * q(izz + jdh, kk, i, 1);
                        f(izz + js - 1, iq1, 1) = tfacp1
                            * dtj(izz + jdh, kk, i) * res(izz + jdh, kk, i, 1)
                            + restmp;
                        f(izz + js - 1, iq1, 5) = f(izz + js - 1, iq1, 1) * ekin
                            + temp * res(izz + jdh, kk, i, 2) * q(izz + jdh, kk, i, 2)
                            + temp * res(izz + jdh, kk, i, 3) * q(izz + jdh, kk, i, 3)
                            + temp * res(izz + jdh, kk, i, 4) * q(izz + jdh, kk, i, 4)
                            + tfacp1 * gm1i * res(izz + jdh, kk, i, 5) * dtj(izz + jdh, kk, i)
                            + restmp * ho;
                        f(izz + js - 1, iq1, 2) = temp * res(izz + jdh, kk, i, 2)
                            + q(izz + jdh, kk, i, 2) * f(izz + js - 1, iq1, 1)
                            + restmp * q(izz + jdh, kk, i, 2);
                        f(izz + js - 1, iq1, 3) = temp * res(izz + jdh, kk, i, 3)
                            + q(izz + jdh, kk, i, 3) * f(izz + js - 1, iq1, 1)
                            + restmp * q(izz + jdh, kk, i, 3);
                        f(izz + js - 1, iq1, 4) = temp * res(izz + jdh, kk, i, 4)
                            + q(izz + jdh, kk, i, 4) * f(izz + js - 1, iq1, 1)
                            + restmp * q(izz + jdh, kk, i, 4);
                    }
                    // q8vrev calls
                    {
                        FortranArray1DRef<double> fa1(&f(js, iq1, 1), jdh);
                        FortranArray1DRef<double> fb1(&f(js, iq2, 1), jdh);
                        q8vrev_ns::q8vrev(jdh, fa1, jdh, fb1);
                    }
                    {
                        FortranArray1DRef<double> fa2(&f(js, iq1, 2), jdh);
                        FortranArray1DRef<double> fb2(&f(js, iq2, 2), jdh);
                        q8vrev_ns::q8vrev(jdh, fa2, jdh, fb2);
                    }
                    {
                        FortranArray1DRef<double> fa3(&f(js, iq1, 3), jdh);
                        FortranArray1DRef<double> fb3(&f(js, iq2, 3), jdh);
                        q8vrev_ns::q8vrev(jdh, fa3, jdh, fb3);
                    }
                    {
                        FortranArray1DRef<double> fa4(&f(js, iq1, 4), jdh);
                        FortranArray1DRef<double> fb4(&f(js, iq2, 4), jdh);
                        q8vrev_ns::q8vrev(jdh, fa4, jdh, fb4);
                    }
                    {
                        FortranArray1DRef<double> fa5(&f(js, iq1, 5), jdh);
                        FortranArray1DRef<double> fb5(&f(js, iq2, 5), jdh);
                        q8vrev_ns::q8vrev(jdh, fa5, jdh, fb5);
                    }
                    // Loop 10021
                    for (int izz = 1; izz <= jdh; izz++) {
                        double c2 = (double)gamma * q(izz, kk, i, 5) / q(izz, kk, i, 1);
                        double c = std::sqrt(c2);
                        double ekin = 0.5 * (q(izz, kk, i, 2) * q(izz, kk, i, 2)
                                           + q(izz, kk, i, 3) * q(izz, kk, i, 3)
                                           + q(izz, kk, i, 4) * q(izz, kk, i, 4));
                        double ho = c2 / (double)gm1 + ekin;
                        double vmag1 = 2.0 * ekin;
                        double avn_uref2 = (double)avn * (double)uref * (double)uref;
                        double vel2 = ccomplex_ns::ccmax(vmag1, avn_uref2);
                        double vel = std::sqrt(ccomplex_ns::ccmin(c2, vel2));
                        vel = (double)cprec * vel + (1.0 - (double)cprec) * c;
                        double thet = (1.0 / (vel * vel) - 1.0 / c2);
                        double restmp = tfacp1 * dtj(izz, kk, i) * thet * res(izz, kk, i, 5);

                        double temp = tfacp1 * dtj(izz, kk, i) * q(izz, kk, i, 1);
                        f(izz + js - 1, iq1, 1) = tfacp1 * dtj(izz, kk, i) * res(izz, kk, i, 1)
                            + restmp;
                        f(izz + js - 1, iq1, 5) = f(izz + js - 1, iq1, 1) * ekin
                            + temp * res(izz, kk, i, 2) * q(izz, kk, i, 2)
                            + temp * res(izz, kk, i, 3) * q(izz, kk, i, 3)
                            + temp * res(izz, kk, i, 4) * q(izz, kk, i, 4)
                            + tfacp1 * gm1i * res(izz, kk, i, 5) * dtj(izz, kk, i)
                            + restmp * ho;
                        f(izz + js - 1, iq1, 2) = temp * res(izz, kk, i, 2)
                            + q(izz, kk, i, 2) * f(izz + js - 1, iq1, 1)
                            + restmp * q(izz, kk, i, 2);
                        f(izz + js - 1, iq1, 3) = temp * res(izz, kk, i, 3)
                            + q(izz, kk, i, 3) * f(izz + js - 1, iq1, 1)
                            + restmp * q(izz, kk, i, 3);
                        f(izz + js - 1, iq1, 4) = temp * res(izz, kk, i, 4)
                            + q(izz, kk, i, 4) * f(izz + js - 1, iq1, 1)
                            + restmp * q(izz, kk, i, 4);
                    }
                }
            }
        }
    }

    // Solve matrix equation
    int il = 1;
    int iu = idim1 * (imw + 1);
    int n  = npl * jdim1 / (imw + 1);

    int id1 = npl * (jdim - 1) / (imw + 1);
    int id2 = (idim - 1) * (imw + 1);
    int i1  = 1;
    bsub_ns::bsub(id1, id2, ai, bi, ci, f, i1, n, il, iu);

    // Update delta q
    if (imw == 0) {
        for (int kpl = 1; kpl <= npl; kpl++) {
            int kk = kp + kpl - 1;
            int jv = (kpl - 1) * jdim1 + 1;
            for (int i = 1; i <= idim1; i++) {
                for (int l = 1; l <= 5; l++) {
                    for (int izz = 1; izz <= jdim1; izz++) {
                        res(izz, kk, i, l) = f(izz + jv - 1, i, l);
                    }
                }
            }
        }
    } else {
        int jdh = jdim1 / 2;
        int idh = idim1 * 2;
        for (int kpl = 1; kpl <= npl; kpl++) {
            int kk = kp + kpl - 1;
            int jv = (kpl - 1) * jdh + 1;
            for (int i = 1; i <= idim1; i++) {
                int iq1 = i;
                int iq2 = idh + 1 - i;
                for (int l = 1; l <= 5; l++) {
                    for (int izz = 1; izz <= jdh; izz++) {
                        res(izz, kk, i, l) = f(izz + jv - 1, iq1, l);
                    }
                    {
                        FortranArray1DRef<double> fsrc(&f(jv, iq2, l), jdh);
                        FortranArray1DRef<double> rdst(&res(jdh + 1, kk, i, l), jdh);
                        q8vrev_ns::q8vrev(jdh, fsrc, jdh, rdst);
                    }
                }
            }
        }
    }

    // label 3000 / return
}

} // namespace swafi_ns
