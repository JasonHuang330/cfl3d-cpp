// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "diagk.h"
#include "tinvr.h"
#include "tdq.h"
#include "dlutr.h"
#include "dfbtr.h"
#include "dabckz.h"
#include "ccomplex.h"
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <cmath>
#include <cstdint>

namespace diagk_ns {

void diagk(int& i, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q,
           FortranArray4DRef<double> res,
           FortranArray3DRef<double> dtj,
           FortranArray4DRef<double> sk,
           FortranArray2DRef<double> t,
           FortranArray3DRef<double> vol,
           FortranArray3DRef<double> vist3d,
           FortranArray3DRef<double> blank,
           int& iover)
{
    // COMMON block aliases
    float& gamma   = cmn_fluid.gamma;
    float& epsa_l  = cmn_entfix.epsa_l;
    float& epsa_r  = cmn_entfix.epsa_r;
    float& xmach   = cmn_info.xmach;
    float& reue    = cmn_reyue.reue;
    int32_t& ita   = cmn_unst.ita;
    float& cprec   = cmn_precond.cprec;
    float& uref    = cmn_precond.uref;
    float& avn     = cmn_precond.avn;
    int32_t& iexp  = cmn_zero.iexp;

    // 10.**(-iexp) is machine zero
    double zero   = std::pow(10.0, (double)(-iexp));
    epsa_l        = 2.f * epsa_r;

    // k-implicit j-sweep line inversions af
    int kdim1 = kdim - 1;
    int jdim1 = jdim - 1;
    double xmre = 8.e0 * (double)xmach / (double)reue;

    double tfacp1;
    if (std::abs(ita) == 1) {
        tfacp1 = 1.e0;
    } else {
        tfacp1 = 1.5e0;
    }

    // load rhs (delta q * dtj) into f
    int jv  = npl * jdim1;

    for (int ipl = 1; ipl <= npl; ipl++) {
        int ii  = i + ipl - 1;
        int jv0 = (ipl - 1) * jdim1 + 1;
        for (int k = 1; k <= kdim1; k++) {
            int jk = (k - 1) * jv + jv0;
            for (int l = 1; l <= 5; l++) {
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + jk - 1, 25 + l) = res(izz, k, ii, l);
                    t(izz + jk - 1, l + 5)  = q(izz, k, ii, l);
                }
            }
            for (int izz = 1; izz <= jdim1; izz++) {
                t(izz + jk - 1, 16) = sk(izz, k, ii, 1);
                t(izz + jk - 1, 17) = sk(izz, k, ii, 2);
                t(izz + jk - 1, 18) = sk(izz, k, ii, 3);
                t(izz + jk - 1, 19) = sk(izz, k, ii, 4);
                t(izz + jk - 1, 20) = sk(izz, k, ii, 5);
                t(izz + jk - 1, 21) = tfacp1 * dtj(izz, k, ii);
            }
            if (cmn_reyue.ivisc[2] > 1) {
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + jk - 1, 12) = vol(izz, k, ii);
                    t(izz + jk - 1, 31) = vist3d(izz, k, ii);
                }
            } else if (cmn_reyue.ivisc[2] > 0) {
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + jk - 1, 12) = vol(izz, k, ii);
                }
            }
        }
        int jk = kdim1 * jv + jv0;
        for (int izz = 1; izz <= jdim1; izz++) {
            t(izz + jk - 1, 16) = sk(izz, kdim, ii, 1);
            t(izz + jk - 1, 17) = sk(izz, kdim, ii, 2);
            t(izz + jk - 1, 18) = sk(izz, kdim, ii, 3);
            t(izz + jk - 1, 19) = sk(izz, kdim, ii, 4);
            t(izz + jk - 1, 20) = sk(izz, kdim, ii, 5);
        }
    }

    int n = jv * kdim1;
    for (int l = 1; l <= 5; l++) {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 25 + l) = t(izz, 25 + l) * t(izz, 21);
        }
    }

    // average metric
    for (int izz = 1; izz <= n; izz++) {
        double t1 = t(izz, 16) + t(izz + jv, 16);
        double t2 = t(izz, 17) + t(izz + jv, 17);
        double t3 = t(izz, 18) + t(izz + jv, 18);
        double t4 = t1 * t1 + t2 * t2 + t3 * t3;
        t4           = 1.e0 / std::sqrt(t4);
        t(izz, 1)    = t1 * t4;
        t(izz, 2)    = t2 * t4;
        t(izz, 3)    = t3 * t4;
        t(izz, 13)   = 0.5 * (t(izz, 20) + t(izz + jv, 20));
    }
    for (int izz = 1; izz <= n + jv; izz++) {
        t(izz, 4) = 0.50 * t(izz, 19);
    }

    // recover primitives
    for (int izz = 1; izz <= n; izz++) {
        t(izz, 16) = 1.e0 / t(izz, 6);
    }


    // viscous term
    if (cmn_reyue.ivisc[2] > 0) {
        if (cmn_reyue.ivisc[2] > 1) {
            for (int izz = 1; izz <= n; izz++) {
                t(izz, 32) = t(izz, 16) * (1.e0 + t(izz, 31));
            }
        } else {
            for (int izz = 1; izz <= n; izz++) {
                t(izz, 32) = t(izz, 16);
            }
        }
        for (int izz = 1; izz <= n + jv; izz++) {
            t(izz, 25) = xmre * t(izz, 4) * t(izz, 4);
        }
        for (int izz = 1; izz <= jv; izz++) {
            t(izz, 25) = t(izz, 25) * t(izz, 32) / t(izz, 12);
        }
        int ns = n - jv;
        for (int izz = 1; izz <= ns; izz++) {
            t(izz + jv, 25) = t(izz + jv, 25) * (t(izz, 32)
                            + t(izz + jv, 32)) / (t(izz, 12) + t(izz + jv, 12));
        }
        for (int izz = 1; izz <= jv; izz++) {
            t(izz + n, 25) = t(izz + n, 25) * t(izz + ns, 32) / t(izz + ns, 12);
        }
    } else {
        for (int izz = 1; izz <= n + jv; izz++) {
            t(izz, 25) = 0.e0;
        }
    }

    if ((float)cprec == 0.f) {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 5)  = std::sqrt((double)gamma * t(izz, 10) * t(izz, 16));
            t(izz, 11) = t(izz, 1) * t(izz, 7) + t(izz, 2) * t(izz, 8)
                       + t(izz, 3) * t(izz, 9) + t(izz, 13);
        }
    } else {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 5)  = std::sqrt((double)gamma * t(izz, 10) * t(izz, 16));
            t(izz, 11) = t(izz, 1) * t(izz, 7) + t(izz, 2) * t(izz, 8)
                       + t(izz, 3) * t(izz, 9) + t(izz, 13);

            // calculation of preconditioning quantities
            double vmag1   = t(izz, 7) * t(izz, 7) + t(izz, 8) * t(izz, 8)
                           + t(izz, 9) * t(izz, 9);
            double avn_d   = (double)avn;
            double uref_d  = (double)uref;
            double avn_uref2 = avn_d * uref_d * uref_d;
            double vel2    = ccomplex_ns::ccmax(vmag1, avn_uref2);
            double c5sq    = t(izz, 5) * t(izz, 5);
            double vel     = std::sqrt(ccomplex_ns::ccmin(c5sq, vel2));
            double cprec_d = (double)cprec;
            vel = cprec_d * vel + (1.0 - cprec_d) * t(izz, 5);
            double xm2     = (vel / t(izz, 5)) * (vel / t(izz, 5));
            double xmave   = t(izz, 11) / t(izz, 5);
            double t11     = 0.5 * (1.0 + xm2);
            double t21     = 0.5 * std::sqrt(xmave * xmave * (1.0 - xm2) * (1.0 - xm2) + 4.0 * xm2);
            t(izz, 33)     = t11 * t(izz, 11) + t21 * t(izz, 5);
            t(izz, 34)     = t11 * t(izz, 11) - t21 * t(izz, 5);
            t(izz, 35)     = xm2 * t(izz, 5);
        }
    }

    // t(inverse) r
    int maxf = jv * kdim;
    {
        FortranArray1DRef<double> r1(&t(1, 26), maxf, 1);
        FortranArray1DRef<double> r2(&t(1, 27), maxf, 1);
        FortranArray1DRef<double> r3(&t(1, 28), maxf, 1);
        FortranArray1DRef<double> r4(&t(1, 29), maxf, 1);
        FortranArray1DRef<double> r5(&t(1, 30), maxf, 1);
        FortranArray1DRef<double> kx(&t(1, 1),  maxf, 1);
        FortranArray1DRef<double> ky(&t(1, 2),  maxf, 1);
        FortranArray1DRef<double> kz(&t(1, 3),  maxf, 1);
        FortranArray1DRef<double> lx(&t(1, 14), maxf, 1);
        FortranArray1DRef<double> ly(&t(1, 15), maxf, 1);
        FortranArray1DRef<double> lz(&t(1, 16), maxf, 1);
        FortranArray1DRef<double> mx(&t(1, 17), maxf, 1);
        FortranArray1DRef<double> my(&t(1, 18), maxf, 1);
        FortranArray1DRef<double> mz(&t(1, 19), maxf, 1);
        FortranArray1DRef<double> c_arr(&t(1, 5),   maxf, 1);
        FortranArray1DRef<double> ub_arr(&t(1, 11),  maxf, 1);
        FortranArray1DRef<double> rho_arr(&t(1, 6),  maxf, 1);
        FortranArray1DRef<double> u_arr(&t(1, 7),   maxf, 1);
        FortranArray1DRef<double> v_arr(&t(1, 8),   maxf, 1);
        FortranArray1DRef<double> w_arr(&t(1, 9),   maxf, 1);
        FortranArray1DRef<double> eig2(&t(1, 33),  maxf, 1);
        FortranArray1DRef<double> eig3(&t(1, 34),  maxf, 1);
        FortranArray1DRef<double> xm2a(&t(1, 35),  maxf, 1);
        int itinv = 1;
        tinvr_ns::tinvr(n, r1, r2, r3, r4, r5,
                        kx, ky, kz, lx, ly, lz, mx, my, mz,
                        c_arr, ub_arr, rho_arr, u_arr, v_arr, w_arr,
                        maxf, itinv, eig2, eig3, xm2a);
    }

    // assemble and solve decoupled matrix equations
    int il = 1;
    int iu = kdim1;

    for (int izz = 1; izz <= n; izz++) {
        t(izz, 31) = t(izz, 11);
        t(izz, 32) = ccomplex_ns::ccabs(t(izz, 31));

        // limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
        if ((float)epsa_l > 0.f) {
            double cc    = ccomplex_ns::ccabs(t(izz, 5));
            double uu    = ccomplex_ns::ccabs(t(izz, 7));
            double vv    = ccomplex_ns::ccabs(t(izz, 8));
            double ww    = ccomplex_ns::ccabs(t(izz, 9));
            double epsaa = (double)epsa_l * (cc + uu + vv + ww);
            double epsbb = 0.25 / ccomplex_ns::ccmax(epsaa, zero);
            double epscc = 2.00 * epsaa;
            if ((float)t(izz, 32) < (float)epscc)
                t(izz, 32) = t(izz, 32) * t(izz, 32) * epsbb + epsaa;
        }

        t(izz, 24) = t(izz, 31) + t(izz, 32);
        t(izz, 31) = t(izz, 31) - t(izz, 32);
        t(izz, 23) = t(izz, 21) + t(izz + jv, 4) * t(izz, 24) - t(izz, 4) * t(izz, 31)
                   + t(izz + jv, 25) + t(izz, 25);
    }
    for (int izz = 1; izz <= n - jv; izz++) {
        t(izz + jv, 22) = -t(izz, 24) * t(izz + jv, 4) - t(izz + jv, 25);
        t(izz, 24)      =  t(izz + jv, 31) * t(izz + jv, 4) - t(izz + jv, 25);
    }

    if (iover == 1) {
        FortranArray3DRef<double> a22(&t(1, 22), jdim-1, npl, kdim);
        FortranArray3DRef<double> a23(&t(1, 23), jdim-1, npl, kdim);
        FortranArray3DRef<double> a24(&t(1, 24), jdim-1, npl, kdim);
        dabckz_ns::dabckz(i, npl, jdim, kdim, idim, a22, a23, a24, blank);
    }

    {
        FortranArray2DRef<double> a22(&t(1, 22), jv, kdim);
        FortranArray2DRef<double> a23(&t(1, 23), jv, kdim);
        FortranArray2DRef<double> a24(&t(1, 24), jv, kdim);
        FortranArray2DRef<double> f26(&t(1, 26), jv, kdim);
        FortranArray2DRef<double> f27(&t(1, 27), jv, kdim);
        FortranArray2DRef<double> f28(&t(1, 28), jv, kdim);
        dlutr_ns::dlutr(jv, jv, kdim, il, iu, a22, a23, a24);
        dfbtr_ns::dfbtr(jv, jv, kdim, il, iu, a22, a23, a24, f26);
        dfbtr_ns::dfbtr(jv, jv, kdim, il, iu, a22, a23, a24, f27);
        dfbtr_ns::dfbtr(jv, jv, kdim, il, iu, a22, a23, a24, f28);
    }


    if ((float)cprec == 0.f) {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 31) = t(izz, 11) + t(izz, 5);
            t(izz, 32) = ccomplex_ns::ccabs(t(izz, 31));

            // limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
            if ((float)epsa_l > 0.f) {
                double cc    = ccomplex_ns::ccabs(t(izz, 5));
                double uu    = ccomplex_ns::ccabs(t(izz, 7));
                double vv    = ccomplex_ns::ccabs(t(izz, 8));
                double ww    = ccomplex_ns::ccabs(t(izz, 9));
                double epsaa = (double)epsa_l * (cc + uu + vv + ww);
                double epsbb = 0.25 / ccomplex_ns::ccmax(epsaa, zero);
                double epscc = 2.00 * epsaa;
                if ((float)t(izz, 32) < (float)epscc)
                    t(izz, 32) = t(izz, 32) * t(izz, 32) * epsbb + epsaa;
            }

            t(izz, 24) = t(izz, 31) + t(izz, 32);
            t(izz, 31) = t(izz, 31) - t(izz, 32);
            t(izz, 23) = t(izz, 21) + t(izz + jv, 4) * t(izz, 24) - t(izz, 4) * t(izz, 31)
                       + t(izz + jv, 25) + t(izz, 25);
        }
    } else {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 31) = t(izz, 33);
            t(izz, 32) = ccomplex_ns::ccabs(t(izz, 31));

            // limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
            if ((float)epsa_l > 0.f) {
                double cc    = ccomplex_ns::ccabs(t(izz, 5));
                double uu    = ccomplex_ns::ccabs(t(izz, 7));
                double vv    = ccomplex_ns::ccabs(t(izz, 8));
                double ww    = ccomplex_ns::ccabs(t(izz, 9));
                double epsaa = (double)epsa_l * (cc + uu + vv + ww);
                double epsbb = 0.25 / ccomplex_ns::ccmax(epsaa, zero);
                double epscc = 2.00 * epsaa;
                if ((float)t(izz, 32) < (float)epscc)
                    t(izz, 32) = t(izz, 32) * t(izz, 32) * epsbb + epsaa;
            }

            t(izz, 24) = t(izz, 31) + t(izz, 32);
            t(izz, 31) = t(izz, 31) - t(izz, 32);
            t(izz, 23) = t(izz, 21) + t(izz + jv, 4) * t(izz, 24) - t(izz, 4) * t(izz, 31)
                       + t(izz + jv, 25) + t(izz, 25);
        }
    }
    for (int izz = 1; izz <= n - jv; izz++) {
        t(izz + jv, 22) = -t(izz, 24) * t(izz + jv, 4) - t(izz + jv, 25);
        t(izz, 24)      =  t(izz + jv, 31) * t(izz + jv, 4) - t(izz + jv, 25);
    }

    if (iover == 1) {
        FortranArray3DRef<double> a22(&t(1, 22), jdim-1, npl, kdim);
        FortranArray3DRef<double> a23(&t(1, 23), jdim-1, npl, kdim);
        FortranArray3DRef<double> a24(&t(1, 24), jdim-1, npl, kdim);
        dabckz_ns::dabckz(i, npl, jdim, kdim, idim, a22, a23, a24, blank);
    }

    {
        FortranArray2DRef<double> a22(&t(1, 22), jv, kdim);
        FortranArray2DRef<double> a23(&t(1, 23), jv, kdim);
        FortranArray2DRef<double> a24(&t(1, 24), jv, kdim);
        FortranArray2DRef<double> f29(&t(1, 29), jv, kdim);
        dlutr_ns::dlutr(jv, jv, kdim, il, iu, a22, a23, a24);
        dfbtr_ns::dfbtr(jv, jv, kdim, il, iu, a22, a23, a24, f29);
    }

    if ((float)cprec == 0.f) {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 31) = t(izz, 11) - t(izz, 5);
            t(izz, 32) = ccomplex_ns::ccabs(t(izz, 31));

            // limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
            if ((float)epsa_l > 0.f) {
                double cc    = ccomplex_ns::ccabs(t(izz, 5));
                double uu    = ccomplex_ns::ccabs(t(izz, 7));
                double vv    = ccomplex_ns::ccabs(t(izz, 8));
                double ww    = ccomplex_ns::ccabs(t(izz, 9));
                double epsaa = (double)epsa_l * (cc + uu + vv + ww);
                double epsbb = 0.25 / ccomplex_ns::ccmax(epsaa, zero);
                double epscc = 2.00 * epsaa;
                if ((float)t(izz, 32) < (float)epscc)
                    t(izz, 32) = t(izz, 32) * t(izz, 32) * epsbb + epsaa;
            }

            t(izz, 24) = t(izz, 31) + t(izz, 32);
            t(izz, 31) = t(izz, 31) - t(izz, 32);
            t(izz, 23) = t(izz, 21) + t(izz + jv, 4) * t(izz, 24) - t(izz, 4) * t(izz, 31)
                       + t(izz + jv, 25) + t(izz, 25);
        }
    } else {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 31) = t(izz, 34);
            t(izz, 32) = ccomplex_ns::ccabs(t(izz, 31));

            // limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
            if ((float)epsa_l > 0.f) {
                double cc    = ccomplex_ns::ccabs(t(izz, 5));
                double uu    = ccomplex_ns::ccabs(t(izz, 7));
                double vv    = ccomplex_ns::ccabs(t(izz, 8));
                double ww    = ccomplex_ns::ccabs(t(izz, 9));
                double epsaa = (double)epsa_l * (cc + uu + vv + ww);
                double epsbb = 0.25 / ccomplex_ns::ccmax(epsaa, zero);
                double epscc = 2.00 * epsaa;
                if ((float)t(izz, 32) < (float)epscc)
                    t(izz, 32) = t(izz, 32) * t(izz, 32) * epsbb + epsaa;
            }

            t(izz, 24) = t(izz, 31) + t(izz, 32);
            t(izz, 31) = t(izz, 31) - t(izz, 32);
            t(izz, 23) = t(izz, 21) + t(izz + jv, 4) * t(izz, 24) - t(izz, 4) * t(izz, 31)
                       + t(izz + jv, 25) + t(izz, 25);
        }
    }
    for (int izz = 1; izz <= n - jv; izz++) {
        t(izz + jv, 22) = -t(izz, 24) * t(izz + jv, 4) - t(izz + jv, 25);
        t(izz, 24)      =  t(izz + jv, 31) * t(izz + jv, 4) - t(izz + jv, 25);
    }

    if (iover == 1) {
        FortranArray3DRef<double> a22(&t(1, 22), jdim-1, npl, kdim);
        FortranArray3DRef<double> a23(&t(1, 23), jdim-1, npl, kdim);
        FortranArray3DRef<double> a24(&t(1, 24), jdim-1, npl, kdim);
        dabckz_ns::dabckz(i, npl, jdim, kdim, idim, a22, a23, a24, blank);
    }

    {
        FortranArray2DRef<double> a22(&t(1, 22), jv, kdim);
        FortranArray2DRef<double> a23(&t(1, 23), jv, kdim);
        FortranArray2DRef<double> a24(&t(1, 24), jv, kdim);
        FortranArray2DRef<double> f30(&t(1, 30), jv, kdim);
        dlutr_ns::dlutr(jv, jv, kdim, il, iu, a22, a23, a24);
        dfbtr_ns::dfbtr(jv, jv, kdim, il, iu, a22, a23, a24, f30);
    }

    // t * delta q
    {
        FortranArray1DRef<double> r1(&t(1, 26), maxf, 1);
        FortranArray1DRef<double> r2(&t(1, 27), maxf, 1);
        FortranArray1DRef<double> r3(&t(1, 28), maxf, 1);
        FortranArray1DRef<double> r4(&t(1, 29), maxf, 1);
        FortranArray1DRef<double> r5(&t(1, 30), maxf, 1);
        FortranArray1DRef<double> kx(&t(1, 1),  maxf, 1);
        FortranArray1DRef<double> ky(&t(1, 2),  maxf, 1);
        FortranArray1DRef<double> kz(&t(1, 3),  maxf, 1);
        FortranArray1DRef<double> lx(&t(1, 14), maxf, 1);
        FortranArray1DRef<double> ly(&t(1, 15), maxf, 1);
        FortranArray1DRef<double> lz(&t(1, 16), maxf, 1);
        FortranArray1DRef<double> mx(&t(1, 17), maxf, 1);
        FortranArray1DRef<double> my(&t(1, 18), maxf, 1);
        FortranArray1DRef<double> mz(&t(1, 19), maxf, 1);
        FortranArray1DRef<double> c_arr(&t(1, 5),   maxf, 1);
        FortranArray1DRef<double> ub_arr(&t(1, 11),  maxf, 1);
        FortranArray1DRef<double> rho_arr(&t(1, 6),  maxf, 1);
        FortranArray1DRef<double> u_arr(&t(1, 7),   maxf, 1);
        FortranArray1DRef<double> v_arr(&t(1, 8),   maxf, 1);
        FortranArray1DRef<double> w_arr(&t(1, 9),   maxf, 1);
        FortranArray1DRef<double> eig2(&t(1, 33),  maxf, 1);
        FortranArray1DRef<double> eig3(&t(1, 34),  maxf, 1);
        FortranArray1DRef<double> xm2a(&t(1, 35),  maxf, 1);
        tdq_ns::tdq(n, r1, r2, r3, r4, r5,
                    kx, ky, kz, lx, ly, lz, mx, my, mz,
                    c_arr, ub_arr, rho_arr, u_arr, v_arr, w_arr,
                    maxf, eig2, eig3, xm2a);
    }

    // update delta q
    for (int ipl = 1; ipl <= npl; ipl++) {
        int ii  = i + ipl - 1;
        int jv0 = (ipl - 1) * jdim1 + 1;
        for (int l = 1; l <= 5; l++) {
            for (int k = 1; k <= kdim1; k++) {
                int jk = jv0 + (k - 1) * jv;
                for (int izz = 1; izz <= jdim1; izz++) {
                    res(izz, k, ii, l) = t(izz + jk - 1, 25 + l);
                }
            }
        }
    }

} // end diagk function

} // namespace diagk_ns
