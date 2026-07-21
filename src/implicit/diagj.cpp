// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "diagj.h"
#include <cstdio>
#include "tinvr.h"
#include "tdq.h"
#include "dlutr.h"
#include "dlutrp.h"
#include "dfbtr.h"
#include "dfbtrp.h"
#include "dabcjz.h"
#include "ccomplex.h"
#include <cmath>

namespace diagj_ns {

void diagj(int& i, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q, FortranArray4DRef<double> res,
           FortranArray3DRef<double> dtj, FortranArray4DRef<double> sj,
           FortranArray2DRef<double> t, int& iperd,
           FortranArray3DRef<double> vol, FortranArray3DRef<double> vist3d,
           FortranArray3DRef<double> blank, int& iover)
{
    // COMMON block aliases
    float& gamma    = cmn_fluid.gamma;
    float& epsa_l   = cmn_entfix.epsa_l;
    float& epsa_r   = cmn_entfix.epsa_r;
    float& xmach    = cmn_info.xmach;
    float& reue     = cmn_reyue.reue;
    int32_t* ivisc  = cmn_reyue.ivisc;  // 0-based C array
    int32_t& ita    = cmn_unst.ita;
    float& cprec    = cmn_precond.cprec;
    float& uref     = cmn_precond.uref;
    float& avn      = cmn_precond.avn;
    int32_t& iexp   = cmn_zero.iexp;

    // 10.**(-iexp) is machine zero
    double zero   = std::pow(10.0, (double)(-iexp));
    epsa_l = 2.f * epsa_r;

    // j-implicit k-sweep line inversions af
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    double xmre = 8.e0 * (double)xmach / (double)reue;
    double tfacp1;
    if (std::abs(ita) == 1) {
        tfacp1 = 1.e0;
    } else {
        tfacp1 = 1.5e0;
    }

    // load rhs (-residual) into f
    int kv = npl * kdim;

    for (int j = 1; j <= jdim1; j++) {
        int kj = (j - 1) * kv + 1;
        for (int l = 1; l <= 5; l++) {
            int jj = 1 - jdim;
            for (int ii = 1; ii <= kv; ii++) {
                jj = jj + jdim;
                t(kj + ii - 1, 25 + l) = -res(j + jj - 1, 1, i, l);
            }
            jj = 1 - jdim;
            for (int ii = 1; ii <= kv; ii++) {
                jj = jj + jdim;
                t(kj + ii - 1, l + 5) = q(j + jj - 1, 1, i, l);
            }
        }
        {
            int jj = 1 - jdim;
            for (int ii = 1; ii <= kv; ii++) {
                jj = jj + jdim;
                t(kj + ii - 1, 21) = tfacp1 * dtj(j + jj - 1, 1, i);
            }
        }
        if (ivisc[1] > 0) {
            int jj2 = 1 - jdim;
            for (int ii = 1; ii <= kv; ii++) {
                jj2 = jj2 + jdim;
                t(kj + ii - 1, 12) = vol(j + jj2 - 1, 1, i);
            }
        }
    }

    for (int j = 1; j <= jdim; j++) {
        int kj = (j - 1) * kv + 1;
        for (int l = 1; l <= 5; l++) {
            int jj = 1 - jdim;
            for (int ii = 1; ii <= kv; ii++) {
                jj = jj + jdim;
                t(kj + ii - 1, 15 + l) = sj(j + jj - 1, 1, i, l);
            }
        }
    }

    if (ivisc[1] > 1) {
        int ic = 0;
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii = i + ipl - 1;
            for (int k = 1; k <= kdim; k++) {
                ic = ic + 1;
                if (k != kdim) {
                    for (int j = 1; j <= jdim1; j++) {
                        int kj = (j - 1) * kv;
                        t(kj + ic, 31) = vist3d(j, k, ii);
                    }
                } else {
                    for (int j = 1; j <= jdim1; j++) {
                        int kj = (j - 1) * kv;
                        t(kj + ic, 31) = vist3d(j, kdim1, ii);
                    }
                }
            }
        }
    }

    int n = kv * jdim1;

    // average metric
    for (int izz = 1; izz <= n; izz++) {
        double t1 = t(izz, 16) + t(izz + kv, 16);
        double t2 = t(izz, 17) + t(izz + kv, 17);
        double t3 = t(izz, 18) + t(izz + kv, 18);
        double t4 = t1 * t1 + t2 * t2 + t3 * t3;
        t4 = 1.e0 / std::sqrt(t4);
        t(izz, 1)  = t1 * t4;
        t(izz, 2)  = t2 * t4;
        t(izz, 3)  = t3 * t4;
        t(izz, 13) = 0.5 * (t(izz, 20) + t(izz + kv, 20));
    }
    for (int izz = 1; izz <= n + kv; izz++) {
        t(izz, 4) = 0.50 * t(izz, 19);
    }

    // recover primitives - viscous term
    if (ivisc[1] > 0) {
        if (ivisc[1] > 1) {
            for (int izz = 1; izz <= n; izz++) {
                t(izz, 32) = (1.e0 + t(izz, 31)) / t(izz, 6);
            }
        } else {
            for (int izz = 1; izz <= n; izz++) {
                t(izz, 32) = 1.e0 / t(izz, 6);
            }
        }
        for (int izz = 1; izz <= n + kv; izz++) {
            t(izz, 25) = xmre * t(izz, 4) * t(izz, 4);
        }
        for (int izz = 1; izz <= kv; izz++) {
            t(izz, 25) = t(izz, 25) * t(izz, 32) / t(izz, 12);
        }
        int ns = n - kv;
        for (int izz = 1; izz <= ns; izz++) {
            t(izz + kv, 25) = t(izz + kv, 25) * (t(izz, 32)
                            + t(izz + kv, 32)) / (t(izz, 12) + t(izz + kv, 12));
        }
        for (int izz = 1; izz <= kv; izz++) {
            t(izz + n, 25) = t(izz + n, 25) * t(izz + ns, 32) / t(izz + ns, 12);
        }
    } else {
        for (int izz = 1; izz <= n + kv; izz++) {
            t(izz, 25) = 0.e0;
        }
    }

    if ((float)cprec == 0.f) {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 5)  = std::sqrt((double)gamma * t(izz, 10) / t(izz, 6));
            t(izz, 11) = t(izz, 1) * t(izz, 7) + t(izz, 2) * t(izz, 8)
                       + t(izz, 3) * t(izz, 9) + t(izz, 13);
        }
    } else {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 5)  = std::sqrt((double)gamma * t(izz, 10) / t(izz, 6));
            t(izz, 11) = t(izz, 1) * t(izz, 7) + t(izz, 2) * t(izz, 8)
                       + t(izz, 3) * t(izz, 9) + t(izz, 13);

            // calculation of preconditioning quantities
            double vmag1 = t(izz, 7) * t(izz, 7) + t(izz, 8) * t(izz, 8)
                         + t(izz, 9) * t(izz, 9);
            double tmp1 = vmag1;
            double tmp2 = (double)avn * (double)uref * (double)uref;
            double vel2 = ccomplex_ns::ccmax(tmp1, tmp2);
            double tmp3 = t(izz, 5) * t(izz, 5);
            double tmp4 = vel2;
            double vel = std::sqrt(ccomplex_ns::ccmin(tmp3, tmp4));
            vel = (double)cprec * vel + (1.0 - (double)cprec) * t(izz, 5);
            double xm2 = (vel / t(izz, 5)) * (vel / t(izz, 5));
            double xmave = t(izz, 11) / t(izz, 5);
            double t11 = 0.5 * (1.0 + xm2);
            double t21 = 0.5 * std::sqrt(xmave * xmave * (1.0 - xm2) * (1.0 - xm2) + 4.0 * xm2);
            t(izz, 33) = t11 * t(izz, 11) + t21 * t(izz, 5);
            t(izz, 34) = t11 * t(izz, 11) - t21 * t(izz, 5);
            t(izz, 35) = xm2 * t(izz, 5);
        }
    }

    // t(inverse) r
    int maxf = kv * jdim;
    {
        FortranArray1DRef<double> r1(&t(1, 26), maxf, 1);
        FortranArray1DRef<double> r2(&t(1, 27), maxf, 1);
        FortranArray1DRef<double> r3(&t(1, 28), maxf, 1);
        FortranArray1DRef<double> r4(&t(1, 29), maxf, 1);
        FortranArray1DRef<double> r5(&t(1, 30), maxf, 1);
        FortranArray1DRef<double>    kx(&t(1, 1),  maxf, 1);
        FortranArray1DRef<double>    ky(&t(1, 2),  maxf, 1);
        FortranArray1DRef<double>    kz(&t(1, 3),  maxf, 1);
        FortranArray1DRef<double>    lx(&t(1, 14), maxf, 1);
        FortranArray1DRef<double>    ly(&t(1, 15), maxf, 1);
        FortranArray1DRef<double>    lz(&t(1, 16), maxf, 1);
        FortranArray1DRef<double>    mx(&t(1, 17), maxf, 1);
        FortranArray1DRef<double>    my(&t(1, 18), maxf, 1);
        FortranArray1DRef<double>    mz(&t(1, 19), maxf, 1);
        FortranArray1DRef<double> c_arr(&t(1, 5),  maxf, 1);
        FortranArray1DRef<double> ub(&t(1, 11), maxf, 1);
        FortranArray1DRef<double> rho(&t(1, 6),  maxf, 1);
        FortranArray1DRef<double> u(&t(1, 7),  maxf, 1);
        FortranArray1DRef<double> v(&t(1, 8),  maxf, 1);
        FortranArray1DRef<double> w(&t(1, 9),  maxf, 1);
        FortranArray1DRef<double> eig2(&t(1, 33), maxf, 1);
        FortranArray1DRef<double> eig3(&t(1, 34), maxf, 1);
        FortranArray1DRef<double> xm2a(&t(1, 35), maxf, 1);
        int itinv = 0;
        tinvr_ns::tinvr(n, r1, r2, r3, r4, r5,
                        kx, ky, kz, lx, ly, lz, mx, my, mz,
                        c_arr, ub, rho, u, v, w,
                        maxf, itinv, eig2, eig3, xm2a);
    }

    // assemble and solve decoupled matrix equations
    int il = 1;
    int iu = jdim1;

    double epsi = 0.0;

    // First eigenvalue: u (entropy/shear)
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
            double tmp_epsaa = epsaa;
            double tmp_zero  = zero;
            double epsbb = 0.25 / ccomplex_ns::ccmax(tmp_epsaa, tmp_zero);
            double epscc = 2.00 * epsaa;
            if ((float)t(izz, 32) < (float)epscc)
                t(izz, 32) = t(izz, 32) * t(izz, 32) * epsbb + epsaa;
        }

        t(izz, 24) = t(izz, 31) + t(izz, 32);
        t(izz, 31) = t(izz, 31) - t(izz, 32);
        t(izz, 23) = t(izz, 21) + t(izz + kv, 4) * t(izz, 24) - t(izz, 4) * t(izz, 31)
                   + t(izz + kv, 25) + t(izz, 25);
    }

    if (iperd == 1) {
        for (int izz = 1; izz <= kv; izz++) {
            t(izz, 22)          = -t(izz + n - kv, 24) * t(izz, 4) - t(izz, 25);
            t(izz + n - kv, 24) =  t(izz, 31) * t(izz + n, 4) - t(izz + n, 25);
        }
    }

    for (int izz = 1; izz <= n - kv; izz++) {
        t(izz + kv, 22) = -t(izz, 24) * t(izz + kv, 4) - t(izz + kv, 25);
        t(izz, 24)      =  t(izz + kv, 31) * t(izz + kv, 4) - t(izz + kv, 25);
    }

    if (iover == 1) {
        FortranArray3DRef<double> a22(&t(1, 22), kdim, npl, jdim);
        FortranArray3DRef<double> a23(&t(1, 23), kdim, npl, jdim);
        FortranArray3DRef<double> a24(&t(1, 24), kdim, npl, jdim);
        dabcjz_ns::dabcjz(i, npl, jdim, kdim, idim, a22, a23, a24, blank);
    }

    if (iperd == 0) {
        FortranArray2DRef<double> a22(&t(1, 22), kv, jdim);
        FortranArray2DRef<double> a23(&t(1, 23), kv, jdim);
        FortranArray2DRef<double> a24(&t(1, 24), kv, jdim);
        FortranArray2DRef<double> f26(&t(1, 26), kv, jdim);
        FortranArray2DRef<double> f27(&t(1, 27), kv, jdim);
        FortranArray2DRef<double> f28(&t(1, 28), kv, jdim);
        dlutr_ns::dlutr(kv, kv, jdim, il, iu, a22, a23, a24);
        dfbtr_ns::dfbtr(kv, kv, jdim, il, iu, a22, a23, a24, f26);
        dfbtr_ns::dfbtr(kv, kv, jdim, il, iu, a22, a23, a24, f27);
        dfbtr_ns::dfbtr(kv, kv, jdim, il, iu, a22, a23, a24, f28);
    } else {
        FortranArray2DRef<double> a22(&t(1, 22), kv, jdim);
        FortranArray2DRef<double> a23(&t(1, 23), kv, jdim);
        FortranArray2DRef<double> a24(&t(1, 24), kv, jdim);
        FortranArray2DRef<double> f26(&t(1, 26), kv, jdim);
        FortranArray2DRef<double> f27(&t(1, 27), kv, jdim);
        FortranArray2DRef<double> f28(&t(1, 28), kv, jdim);
        FortranArray2DRef<double> g31(&t(1, 31), kv, jdim);
        FortranArray2DRef<double> h32(&t(1, 32), kv, jdim);
        dlutrp_ns::dlutrp(kv, kv, jdim, il, iu, a22, a23, a24, g31, h32);
        dfbtrp_ns::dfbtrp(kv, kv, jdim, il, iu, a22, a23, a24, f26, g31, h32);
        dfbtrp_ns::dfbtrp(kv, kv, jdim, il, iu, a22, a23, a24, f27, g31, h32);
        dfbtrp_ns::dfbtrp(kv, kv, jdim, il, iu, a22, a23, a24, f28, g31, h32);
    }


    // Second eigenvalue: u+c (or preconditioned eig2)
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
                double tmp_epsaa = epsaa;
                double tmp_zero  = zero;
                double epsbb = 0.25 / ccomplex_ns::ccmax(tmp_epsaa, tmp_zero);
                double epscc = 2.00 * epsaa;
                if ((float)t(izz, 32) < (float)epscc)
                    t(izz, 32) = t(izz, 32) * t(izz, 32) * epsbb + epsaa;
            }

            t(izz, 24) = t(izz, 31) + t(izz, 32);
            t(izz, 31) = t(izz, 31) - t(izz, 32);
            t(izz, 23) = t(izz, 21) + t(izz + kv, 4) * t(izz, 24) - t(izz, 4) * t(izz, 31)
                       + t(izz + kv, 25) + t(izz, 25);
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
                double tmp_epsaa = epsaa;
                double tmp_zero  = zero;
                double epsbb = 0.25 / ccomplex_ns::ccmax(tmp_epsaa, tmp_zero);
                double epscc = 2.00 * epsaa;
                if ((float)t(izz, 32) < (float)epscc)
                    t(izz, 32) = t(izz, 32) * t(izz, 32) * epsbb + epsaa;
            }

            t(izz, 24) = t(izz, 31) + t(izz, 32);
            t(izz, 31) = t(izz, 31) - t(izz, 32);
            t(izz, 23) = t(izz, 21) + t(izz + kv, 4) * t(izz, 24) - t(izz, 4) * t(izz, 31)
                       + t(izz + kv, 25) + t(izz, 25);
        }
    }

    if (iperd == 1) {
        for (int izz = 1; izz <= kv; izz++) {
            t(izz, 22)          = -t(izz + n - kv, 24) * t(izz, 4) - t(izz, 25);
            t(izz + n - kv, 24) =  t(izz, 31) * t(izz + n, 4) - t(izz + n, 25);
        }
    }

    for (int izz = 1; izz <= n - kv; izz++) {
        t(izz + kv, 22) = -t(izz, 24) * t(izz + kv, 4) - t(izz + kv, 25);
        t(izz, 24)      =  t(izz + kv, 31) * t(izz + kv, 4) - t(izz + kv, 25);
    }

    if (iover == 1) {
        FortranArray3DRef<double> a22(&t(1, 22), kdim, npl, jdim);
        FortranArray3DRef<double> a23(&t(1, 23), kdim, npl, jdim);
        FortranArray3DRef<double> a24(&t(1, 24), kdim, npl, jdim);
        dabcjz_ns::dabcjz(i, npl, jdim, kdim, idim, a22, a23, a24, blank);
    }

    if (iperd == 0) {
        FortranArray2DRef<double> a22(&t(1, 22), kv, jdim);
        FortranArray2DRef<double> a23(&t(1, 23), kv, jdim);
        FortranArray2DRef<double> a24(&t(1, 24), kv, jdim);
        FortranArray2DRef<double> f29(&t(1, 29), kv, jdim);
        dlutr_ns::dlutr(kv, kv, jdim, il, iu, a22, a23, a24);
        dfbtr_ns::dfbtr(kv, kv, jdim, il, iu, a22, a23, a24, f29);
    } else {
        FortranArray2DRef<double> a22(&t(1, 22), kv, jdim);
        FortranArray2DRef<double> a23(&t(1, 23), kv, jdim);
        FortranArray2DRef<double> a24(&t(1, 24), kv, jdim);
        FortranArray2DRef<double> f29(&t(1, 29), kv, jdim);
        FortranArray2DRef<double> g31(&t(1, 31), kv, jdim);
        FortranArray2DRef<double> h32(&t(1, 32), kv, jdim);
        dlutrp_ns::dlutrp(kv, kv, jdim, il, iu, a22, a23, a24, g31, h32);
        dfbtrp_ns::dfbtrp(kv, kv, jdim, il, iu, a22, a23, a24, f29, g31, h32);
    }


    // Third eigenvalue: u-c (or preconditioned eig3)
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
                double tmp_epsaa = epsaa;
                double tmp_zero  = zero;
                double epsbb = 0.25 / ccomplex_ns::ccmax(tmp_epsaa, tmp_zero);
                double epscc = 2.00 * epsaa;
                if ((float)t(izz, 32) < (float)epscc)
                    t(izz, 32) = t(izz, 32) * t(izz, 32) * epsbb + epsaa;
            }

            t(izz, 24) = t(izz, 31) + t(izz, 32);
            t(izz, 31) = t(izz, 31) - t(izz, 32);
            t(izz, 23) = t(izz, 21) + t(izz + kv, 4) * t(izz, 24) - t(izz, 4) * t(izz, 31)
                       + t(izz + kv, 25) + t(izz, 25);
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
                double tmp_epsaa = epsaa;
                double tmp_zero  = zero;
                double epsbb = 0.25 / ccomplex_ns::ccmax(tmp_epsaa, tmp_zero);
                double epscc = 2.00 * epsaa;
                if ((float)t(izz, 32) < (float)epscc)
                    t(izz, 32) = t(izz, 32) * t(izz, 32) * epsbb + epsaa;
            }

            t(izz, 24) = t(izz, 31) + t(izz, 32);
            t(izz, 31) = t(izz, 31) - t(izz, 32);
            t(izz, 23) = t(izz, 21) + t(izz + kv, 4) * t(izz, 24) - t(izz, 4) * t(izz, 31)
                       + t(izz + kv, 25) + t(izz, 25);
        }
    }

    if (iperd == 1) {
        for (int izz = 1; izz <= kv; izz++) {
            t(izz, 22)          = -t(izz + n - kv, 24) * t(izz, 4) - t(izz, 25);
            t(izz + n - kv, 24) =  t(izz, 31) * t(izz + n, 4) - t(izz + n, 25);
        }
    }

    for (int izz = 1; izz <= n - kv; izz++) {
        t(izz + kv, 22) = -t(izz, 24) * t(izz + kv, 4) - t(izz + kv, 25);
        t(izz, 24)      =  t(izz + kv, 31) * t(izz + kv, 4) - t(izz + kv, 25);
    }

    if (iover == 1) {
        FortranArray3DRef<double> a22(&t(1, 22), kdim, npl, jdim);
        FortranArray3DRef<double> a23(&t(1, 23), kdim, npl, jdim);
        FortranArray3DRef<double> a24(&t(1, 24), kdim, npl, jdim);
        dabcjz_ns::dabcjz(i, npl, jdim, kdim, idim, a22, a23, a24, blank);
    }

    if (iperd == 0) {
        FortranArray2DRef<double> a22(&t(1, 22), kv, jdim);
        FortranArray2DRef<double> a23(&t(1, 23), kv, jdim);
        FortranArray2DRef<double> a24(&t(1, 24), kv, jdim);
        FortranArray2DRef<double> f30(&t(1, 30), kv, jdim);
        dlutr_ns::dlutr(kv, kv, jdim, il, iu, a22, a23, a24);
        dfbtr_ns::dfbtr(kv, kv, jdim, il, iu, a22, a23, a24, f30);
    } else {
        FortranArray2DRef<double> a22(&t(1, 22), kv, jdim);
        FortranArray2DRef<double> a23(&t(1, 23), kv, jdim);
        FortranArray2DRef<double> a24(&t(1, 24), kv, jdim);
        FortranArray2DRef<double> f30(&t(1, 30), kv, jdim);
        FortranArray2DRef<double> g31(&t(1, 31), kv, jdim);
        FortranArray2DRef<double> h32(&t(1, 32), kv, jdim);
        dlutrp_ns::dlutrp(kv, kv, jdim, il, iu, a22, a23, a24, g31, h32);
        dfbtrp_ns::dfbtrp(kv, kv, jdim, il, iu, a22, a23, a24, f30, g31, h32);
    }

    // t * delta q
    {
        FortranArray1DRef<double> r1(&t(1, 26), maxf, 1);
        FortranArray1DRef<double> r2(&t(1, 27), maxf, 1);
        FortranArray1DRef<double> r3(&t(1, 28), maxf, 1);
        FortranArray1DRef<double> r4(&t(1, 29), maxf, 1);
        FortranArray1DRef<double> r5(&t(1, 30), maxf, 1);
        FortranArray1DRef<double>    kx(&t(1, 1),  maxf, 1);
        FortranArray1DRef<double>    ky(&t(1, 2),  maxf, 1);
        FortranArray1DRef<double>    kz(&t(1, 3),  maxf, 1);
        FortranArray1DRef<double>    lx(&t(1, 14), maxf, 1);
        FortranArray1DRef<double>    ly(&t(1, 15), maxf, 1);
        FortranArray1DRef<double>    lz(&t(1, 16), maxf, 1);
        FortranArray1DRef<double>    mx(&t(1, 17), maxf, 1);
        FortranArray1DRef<double>    my(&t(1, 18), maxf, 1);
        FortranArray1DRef<double>    mz(&t(1, 19), maxf, 1);
        FortranArray1DRef<double> c_arr(&t(1, 5),  maxf, 1);
        FortranArray1DRef<double> ub(&t(1, 11), maxf, 1);
        FortranArray1DRef<double> rho(&t(1, 6),  maxf, 1);
        FortranArray1DRef<double> u(&t(1, 7),  maxf, 1);
        FortranArray1DRef<double> v(&t(1, 8),  maxf, 1);
        FortranArray1DRef<double> w(&t(1, 9),  maxf, 1);
        FortranArray1DRef<double> eig2(&t(1, 33), maxf, 1);
        FortranArray1DRef<double> eig3(&t(1, 34), maxf, 1);
        FortranArray1DRef<double> xm2a(&t(1, 35), maxf, 1);
        tdq_ns::tdq(n, r1, r2, r3, r4, r5,
                    kx, ky, kz, lx, ly, lz, mx, my, mz,
                    c_arr, ub, rho, u, v, w,
                    maxf, eig2, eig3, xm2a);
    }

    // update delta q
    for (int j = 1; j <= jdim1; j++) {
        int kj = (j - 1) * kv + 1;
        for (int l = 1; l <= 5; l++) {
            int jj = 1 - jdim;
            for (int ii = 1; ii <= kv; ii++) {
                jj = jj + jdim;
                res(j + jj - 1, 1, i, l) = t(kj + ii - 1, 25 + l);
            }
        }
    }

    for (int ipl = 1; ipl <= npl; ipl++) {
        int ii = i + ipl - 1;
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= jdim1; izz++) {
                res(izz, kdim, ii, l) = 0.e0;
            }
        }
    }

} // end diagj

} // namespace diagj_ns
