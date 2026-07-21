// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "diagi.h"
#include "tinvr.h"
#include "tdq.h"
#include "dlutr.h"
#include "dfbtr.h"
#include "dabciz.h"
#include "ccomplex.h"
#include <cmath>

namespace diagi_ns {

void diagi(int& k, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q, FortranArray4DRef<double> res,
           FortranArray3DRef<double> dtj, FortranArray4DRef<double> si,
           FortranArray2DRef<double> t, FortranArray3DRef<double> vol,
           FortranArray3DRef<double> vist3d, FortranArray3DRef<double> blank,
           int& iover)
{
    using namespace ccomplex_ns;

    // COMMON block aliases
    float& gamma_f  = cmn_fluid.gamma;
    float& reue_f   = cmn_reyue.reue;
    int32_t* ivisc  = cmn_reyue.ivisc;  // 0-based: ivisc[0]=ivisc(1)
    int32_t& ita    = cmn_unst.ita;
    float& cprec    = cmn_precond.cprec;
    float& uref_f   = cmn_precond.uref;
    float& avn_f    = cmn_precond.avn;
    float& epsa_l   = cmn_entfix.epsa_l;
    float& epsa_r   = cmn_entfix.epsa_r;
    int32_t& iexp   = cmn_zero.iexp;
    float& xmach_f  = cmn_info.xmach;

    // 10.**(-iexp) is machine zero
    double zero_val = std::pow(10.0, (double)(-iexp));
    epsa_l = 2.f * epsa_r;

    // i-implicit j-sweep line inversions af
    int idim1 = idim - 1;
    int jdim1 = jdim - 1;
    double xmre = 8.e0 * (double)xmach_f / (double)reue_f;

    double tfacp1;
    if (std::abs(ita) == 1) {
        tfacp1 = 1.e0;
    } else {
        tfacp1 = 1.5e0;
    }

    // load rhs (delta q * dtj) into f
    int jv = npl * jdim1;

    for (int kpl = 1; kpl <= npl; kpl++) {
        int kk  = k + kpl - 1;
        int jv0 = (kpl - 1) * jdim1 + 1;
        for (int i = 1; i <= idim1; i++) {
            int ji = (i - 1) * jv + jv0;
            for (int l = 1; l <= 5; l++) {
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + ji - 1, 25 + l) = res(izz, kk, i, l);
                    t(izz + ji - 1, l + 5)  = q(izz, kk, i, l);
                }
            }
            for (int izz = 1; izz <= jdim1; izz++) {
                t(izz + ji - 1, 16) = si(izz, kk, i, 1);
                t(izz + ji - 1, 17) = si(izz, kk, i, 2);
                t(izz + ji - 1, 18) = si(izz, kk, i, 3);
                t(izz + ji - 1, 19) = si(izz, kk, i, 4);
                t(izz + ji - 1, 20) = si(izz, kk, i, 5);
                t(izz + ji - 1, 21) = tfacp1 * dtj(izz, kk, i);
            }
            if (ivisc[0] > 1) {
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + ji - 1, 12) = vol(izz, kk, i);
                    t(izz + ji - 1, 31) = vist3d(izz, kk, i);
                }
            } else if (ivisc[0] > 0) {
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + ji - 1, 12) = vol(izz, kk, i);
                }
            }
        }
        // load last i-face metrics
        int ji_last = idim1 * jv + jv0;
        for (int izz = 1; izz <= jdim1; izz++) {
            t(izz + ji_last - 1, 16) = si(izz, kk, idim, 1);
            t(izz + ji_last - 1, 17) = si(izz, kk, idim, 2);
            t(izz + ji_last - 1, 18) = si(izz, kk, idim, 3);
            t(izz + ji_last - 1, 19) = si(izz, kk, idim, 4);
            t(izz + ji_last - 1, 20) = si(izz, kk, idim, 5);
        }
    }

    int n = jv * idim1;

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
        t4 = 1.e0 / std::sqrt(t4);
        t(izz, 1)  = t1 * t4;
        t(izz, 2)  = t2 * t4;
        t(izz, 3)  = t3 * t4;
        t(izz, 13) = 0.5 * (t(izz, 20) + t(izz + jv, 20));
    }
    for (int izz = 1; izz <= n + jv; izz++) {
        t(izz, 4) = 0.50 * t(izz, 19);
    }

    // recover primitives / viscous terms
    if (ivisc[0] > 0) {
        if (ivisc[0] > 1) {
            for (int izz = 1; izz <= n; izz++) {
                t(izz, 32) = (1.e0 + t(izz, 31)) / t(izz, 6);
            }
        } else {
            // Fortran: do 1214 izz=i,n  (i retains last loop value = idim1)
            for (int izz = idim1; izz <= n; izz++) {
                t(izz, 32) = 1.e0 / t(izz, 6);
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
            t(izz + jv, 25) = t(izz + jv, 25) * (t(izz, 32) + t(izz + jv, 32)) /
                              (t(izz, 12) + t(izz + jv, 12));
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
            t(izz, 5)  = std::sqrt((double)gamma_f * t(izz, 10) / t(izz, 6));
            t(izz, 11) = t(izz, 1) * t(izz, 7) + t(izz, 2) * t(izz, 8)
                       + t(izz, 3) * t(izz, 9) + t(izz, 13);
        }
    } else {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 5)  = std::sqrt((double)gamma_f * t(izz, 10) / t(izz, 6));
            t(izz, 11) = t(izz, 1) * t(izz, 7) + t(izz, 2) * t(izz, 8)
                       + t(izz, 3) * t(izz, 9) + t(izz, 13);
            // calculation of preconditioning quantities
            double vmag1 = t(izz, 7) * t(izz, 7) + t(izz, 8) * t(izz, 8)
                         + t(izz, 9) * t(izz, 9);
            double avn_d  = (double)avn_f;
            double uref_d = (double)uref_f;
            double avn_uref2 = avn_d * uref_d * uref_d;
            double vel2 = ccmax(vmag1, avn_uref2);
            double t5sq = t(izz, 5) * t(izz, 5);
            double vel  = std::sqrt(ccmin(t5sq, vel2));
            double cprec_d = (double)cprec;
            vel  = cprec_d * vel + (1.0 - cprec_d) * t(izz, 5);
            double xm2   = (vel / t(izz, 5)) * (vel / t(izz, 5));
            double xmave = t(izz, 11) / t(izz, 5);
            double t11   = 0.5 * (1.0 + xm2);
            double t21   = 0.5 * std::sqrt(xmave * xmave * (1.0 - xm2) * (1.0 - xm2) + 4.0 * xm2);
            t(izz, 33) = t11 * t(izz, 11) + t21 * t(izz, 5);
            t(izz, 34) = t11 * t(izz, 11) - t21 * t(izz, 5);
            t(izz, 35) = xm2 * t(izz, 5);
        }
    }

    // t(inverse) r
    int maxf = jv * idim;

    // Build column refs for tinvr/tdq
    FortranArray1DRef<double> r1(&t(1, 26), maxf, 1);
    FortranArray1DRef<double> r2(&t(1, 27), maxf, 1);
    FortranArray1DRef<double> r3(&t(1, 28), maxf, 1);
    FortranArray1DRef<double> r4(&t(1, 29), maxf, 1);
    FortranArray1DRef<double> r5(&t(1, 30), maxf, 1);
    // Int columns (unit vectors stored as doubles, reinterpret for tinvr/tdq)
    FortranArray1DRef<double> kx_ref(&t(1, 1),  maxf, 1);
    FortranArray1DRef<double> ky_ref(&t(1, 2),  maxf, 1);
    FortranArray1DRef<double> kz_ref(&t(1, 3),  maxf, 1);
    FortranArray1DRef<double> lx_ref(&t(1, 14), maxf, 1);
    FortranArray1DRef<double> ly_ref(&t(1, 15), maxf, 1);
    FortranArray1DRef<double> lz_ref(&t(1, 16), maxf, 1);
    FortranArray1DRef<double> mx_ref(&t(1, 17), maxf, 1);
    FortranArray1DRef<double> my_ref(&t(1, 18), maxf, 1);
    FortranArray1DRef<double> mz_ref(&t(1, 19), maxf, 1);
    FortranArray1DRef<double> c_ref(&t(1, 5),    maxf, 1);
    FortranArray1DRef<double> ub_ref(&t(1, 11),  maxf, 1);
    FortranArray1DRef<double> rho_ref(&t(1, 6),  maxf, 1);
    FortranArray1DRef<double> u_ref(&t(1, 7),    maxf, 1);
    FortranArray1DRef<double> v_ref(&t(1, 8),    maxf, 1);
    FortranArray1DRef<double> w_ref(&t(1, 9),    maxf, 1);
    FortranArray1DRef<double> eig2_ref(&t(1, 33), maxf, 1);
    FortranArray1DRef<double> eig3_ref(&t(1, 34), maxf, 1);
    FortranArray1DRef<double> xm2a_ref(&t(1, 35), maxf, 1);

    int itinv = 1;
    tinvr_ns::tinvr(n, r1, r2, r3, r4, r5,
                    kx_ref, ky_ref, kz_ref,
                    lx_ref, ly_ref, lz_ref,
                    mx_ref, my_ref, mz_ref,
                    c_ref, ub_ref, rho_ref, u_ref, v_ref, w_ref,
                    maxf, itinv,
                    eig2_ref, eig3_ref, xm2a_ref);

    // assemble and solve decoupled matrix equations
    int il = 1;
    int iu = idim1;

    for (int izz = 1; izz <= n; izz++) {
        t(izz, 31) = t(izz, 11);
        t(izz, 32) = ccabs(t(izz, 31));
        // limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
        if ((float)epsa_l > 0.f) {
            double cc    = ccabs(t(izz, 5));
            double uu    = ccabs(t(izz, 7));
            double vv    = ccabs(t(izz, 8));
            double ww    = ccabs(t(izz, 9));
            double epsaa = (double)epsa_l * (cc + uu + vv + ww);
            double epsbb = 0.25 / ccmax(epsaa, zero_val);
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
        FortranArray3DRef<double> a22(&t(1, 22), jdim-1, npl, idim);
        FortranArray3DRef<double> a23(&t(1, 23), jdim-1, npl, idim);
        FortranArray3DRef<double> a24(&t(1, 24), jdim-1, npl, idim);
        dabciz_ns::dabciz(k, npl, jdim, kdim, idim, a22, a23, a24, blank);
    }
    {
        FortranArray2DRef<double> a22(&t(1, 22), jv, idim);
        FortranArray2DRef<double> a23(&t(1, 23), jv, idim);
        FortranArray2DRef<double> a24(&t(1, 24), jv, idim);
        FortranArray2DRef<double> f26(&t(1, 26), jv, idim);
        FortranArray2DRef<double> f27(&t(1, 27), jv, idim);
        FortranArray2DRef<double> f28(&t(1, 28), jv, idim);
        dlutr_ns::dlutr(jv, jv, idim, il, iu, a22, a23, a24);
        dfbtr_ns::dfbtr(jv, jv, idim, il, iu, a22, a23, a24, f26);
        dfbtr_ns::dfbtr(jv, jv, idim, il, iu, a22, a23, a24, f27);
        dfbtr_ns::dfbtr(jv, jv, idim, il, iu, a22, a23, a24, f28);
    }

    if ((float)cprec == 0.f) {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 31) = t(izz, 11) + t(izz, 5);
            t(izz, 32) = ccabs(t(izz, 31));
            // limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
            if ((float)epsa_l > 0.f) {
                double cc    = ccabs(t(izz, 5));
                double uu    = ccabs(t(izz, 7));
                double vv    = ccabs(t(izz, 8));
                double ww    = ccabs(t(izz, 9));
                double epsaa = (double)epsa_l * (cc + uu + vv + ww);
                double epsbb = 0.25 / ccmax(epsaa, zero_val);
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
            t(izz, 32) = ccabs(t(izz, 31));
            // limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
            if ((float)epsa_l > 0.f) {
                double cc    = ccabs(t(izz, 5));
                double uu    = ccabs(t(izz, 7));
                double vv    = ccabs(t(izz, 8));
                double ww    = ccabs(t(izz, 9));
                double epsaa = (double)epsa_l * (cc + uu + vv + ww);
                double epsbb = 0.25 / ccmax(epsaa, zero_val);
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
        FortranArray3DRef<double> a22(&t(1, 22), jdim-1, npl, idim);
        FortranArray3DRef<double> a23(&t(1, 23), jdim-1, npl, idim);
        FortranArray3DRef<double> a24(&t(1, 24), jdim-1, npl, idim);
        dabciz_ns::dabciz(k, npl, jdim, kdim, idim, a22, a23, a24, blank);
    }
    {
        FortranArray2DRef<double> a22(&t(1, 22), jv, idim);
        FortranArray2DRef<double> a23(&t(1, 23), jv, idim);
        FortranArray2DRef<double> a24(&t(1, 24), jv, idim);
        FortranArray2DRef<double> f29(&t(1, 29), jv, idim);
        dlutr_ns::dlutr(jv, jv, idim, il, iu, a22, a23, a24);
        dfbtr_ns::dfbtr(jv, jv, idim, il, iu, a22, a23, a24, f29);
    }

    if ((float)cprec == 0.f) {
        for (int izz = 1; izz <= n; izz++) {
            t(izz, 31) = t(izz, 11) - t(izz, 5);
            t(izz, 32) = ccabs(t(izz, 31));
            // limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
            if ((float)epsa_l > 0.f) {
                double cc    = ccabs(t(izz, 5));
                double uu    = ccabs(t(izz, 7));
                double vv    = ccabs(t(izz, 8));
                double ww    = ccabs(t(izz, 9));
                double epsaa = (double)epsa_l * (cc + uu + vv + ww);
                double epsbb = 0.25 / ccmax(epsaa, zero_val);
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
            t(izz, 32) = ccabs(t(izz, 31));
            // limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
            if ((float)epsa_l > 0.f) {
                double cc    = ccabs(t(izz, 5));
                double uu    = ccabs(t(izz, 7));
                double vv    = ccabs(t(izz, 8));
                double ww    = ccabs(t(izz, 9));
                double epsaa = (double)epsa_l * (cc + uu + vv + ww);
                double epsbb = 0.25 / ccmax(epsaa, zero_val);
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
        FortranArray3DRef<double> a22(&t(1, 22), jdim-1, npl, idim);
        FortranArray3DRef<double> a23(&t(1, 23), jdim-1, npl, idim);
        FortranArray3DRef<double> a24(&t(1, 24), jdim-1, npl, idim);
        dabciz_ns::dabciz(k, npl, jdim, kdim, idim, a22, a23, a24, blank);
    }
    {
        FortranArray2DRef<double> a22(&t(1, 22), jv, idim);
        FortranArray2DRef<double> a23(&t(1, 23), jv, idim);
        FortranArray2DRef<double> a24(&t(1, 24), jv, idim);
        FortranArray2DRef<double> f30(&t(1, 30), jv, idim);
        dlutr_ns::dlutr(jv, jv, idim, il, iu, a22, a23, a24);
        dfbtr_ns::dfbtr(jv, jv, idim, il, iu, a22, a23, a24, f30);
    }

    // t * delta q
    tdq_ns::tdq(n, r1, r2, r3, r4, r5,
                kx_ref, ky_ref, kz_ref,
                lx_ref, ly_ref, lz_ref,
                mx_ref, my_ref, mz_ref,
                c_ref, ub_ref, rho_ref, u_ref, v_ref, w_ref,
                maxf,
                eig2_ref, eig3_ref, xm2a_ref);

    // update delta q
    for (int kpl = 1; kpl <= npl; kpl++) {
        int kk  = k + kpl - 1;
        int jv0 = (kpl - 1) * jdim1 + 1;
        for (int l = 1; l <= 5; l++) {
            for (int i = 1; i <= idim1; i++) {
                int ji = jv0 + (i - 1) * jv;
                for (int izz = 1; izz <= jdim1; izz++) {
                    res(izz, kk, i, l) = t(izz + ji - 1, 25 + l);
                }
            }
        }
    }
} // end diagi

} // namespace diagi_ns
