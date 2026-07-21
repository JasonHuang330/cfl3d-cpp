// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "gfluxv.h"
#include <cstdio>
#include "ccomplex.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace gfluxv_ns {

void gfluxv(int& i, int& npl, int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> res, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray4DRef<double> sj,
            FortranArray3DRef<double> vol, FortranArray2DRef<double> t,
            int& nvtq, FortranArray2DRef<double> wj0,
            FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vmuj,
            FortranArray4DRef<double> vj0, FortranArray3DRef<double> bcj,
            FortranArray4DRef<double> zksav, FortranArray4DRef<double> tj0,
            FortranArray3DRef<double> cmuv, FortranArray3DRef<double> volj0,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& iadv, int& nummem,
            FortranArray4DRef<double> ux)
{
    // COMMON block aliases
    float& gamma    = cmn_fluid.gamma;
    float& gm1      = cmn_fluid.gm1;
    float& pr       = cmn_fluid2.pr;
    float& prt      = cmn_fluid2.prt;
    float& cbar     = cmn_fluid2.cbar;
    float& xmach    = cmn_info.xmach;
    float& reue     = cmn_reyue.reue;
    float& tinf     = cmn_reyue.tinf;
    int*   ivisc    = cmn_reyue.ivisc;   // 0-based: ivisc[1] = ivisc(2)
    int&   isklton  = cmn_sklton.isklton;
    int&   level    = cmn_mgrd.level;
    int&   lglobal  = cmn_mgrd.lglobal;
    float* tur10    = cmn_ivals.tur10;   // 0-based: tur10[0] = tur10(1)
    int&   i_nonlin   = cmn_constit.i_nonlin;
    float& c_nonlin   = cmn_constit.c_nonlin;
    float& snonlin_lim= cmn_constit.snonlin_lim;
    int&   i_tauijs   = cmn_constit.i_tauijs;
    int&   i_qcr2000  = cmn_constit.i_qcr2000;
    int&   i_qcr2013  = cmn_constit.i_qcr2013;
    int&   i_qcr2013v = cmn_constit.i_qcr2013v;
    float& c10      = cmn_easmv.c10;
    float& c11      = cmn_easmv.c11;
    float& c2_e     = cmn_easmv.c2;
    float& c3_e     = cmn_easmv.c3;
    float& c4_e     = cmn_easmv.c4;
    float& c5_e     = cmn_easmv.c5;
    float& sigk1    = cmn_easmv.sigk1;
    float& cmuc1    = cmn_easmv.cmuc1;
    int&   ieasm_type = cmn_easmv.ieasm_type;

    // Local variables
    double coef_eddy, ccr2;
    int kdim1, jdim1;
    int n, l0, kv, ks, nn;
    double xmre, gpr, gm1pr, prtr;
    int idurbinlim;
    int ipl, ii, j, k, izz, l, jkv, jkv2, kc, m, mm, mb;
    int imin, imax, kmin, kmax;
    double t5, t6, t1, ab, bb, wj05, wj06;
    int n1, l1;
    double c2b, c2bp;
    double tr1;
    int iviscc;
    double t25, t24;
    double re, gg, f7;
    double s11, s22, s33, s12, s13, s23;
    double w12, w13, w23;
    double xis, wis;
    double alpa1, alpa2;
    double rho, u, v, w_vel;
    double zksav1, zksav2, bnum, eta, squig, factre;
    double t11, t22, t33, t12, t13, t23;
    double xlnpt;
    double denom;
    double o11, o12, o13, o22, o23, o33;
    double t11p, t22p, t33p, t12p, t13p, t23p;
    double xis2013, s_mod, w_mod;
    double al10, al1, al2, al3, al4;
    double cmuu, eta1, xmut, tau, fmu, zkolmog;
    double tmp;
    double ux1, ux2, ux3, ux4, ux5, ux6, ux7, ux8, ux9;
    double zk_sa;
    double v3d;
    double jp1_d, jm2_d; // not used as double, just int
    int jp1, jm2;
    double c4new, eta1_girimaji;
    double zksav1_tmp, zksav2_tmp; // temps for ccmax calls
    double eta_tmp, squig_tmp;
    double t11_tmp, t22_tmp, t33_tmp;
    double s_mod_tmp, w_mod_tmp;
    double tau_tmp, zkolmog_tmp;
    double denom_tmp;
    double tur10_0_d;



    // coefficient used to turn off eddy-viscosity based
    // turbulent stress calculation for ivmx==70 (full Reynolds stress model)
    coef_eddy = 1.0;
    if (ivisc[1] >= 70) coef_eddy = 0.0;
    if (i_tauijs == 1) coef_eddy = 0.0;
    ccr2 = 2.5;

    kdim1 = kdim - 1;
    jdim1 = jdim - 1;

    // n  : number of cell centers for npl planes
    // l0 : number of cell interfaces for npl planes
    // kv : number of cell centers (and interfaces) on a j=constant plane
    n  = npl * jdim1 * kdim1;
    l0 = npl * kdim1 * jdim;
    kv = npl * kdim1;
    ks = kv + 1;
    nn = n - kv;

    xmre  = (double)xmach / (double)reue;
    gpr   = (double)gamma / (double)pr;
    gm1pr = (double)gm1 * (double)pr;
    prtr  = (double)pr / (double)prt;
    // Durbin TCFD 1991 near-wall limiter (0=off)
    idurbinlim = 0;

    // isklton diagnostic output
    if (isklton > 0 && i == 1 && iadv >= 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "   computing viscous  fluxes, J-direction - central differencing");
        if (i_tauijs == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120, "      using tau_ij method");
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "(note: not done accurately at boundaries for multiblock cases...)");
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "currently one sided; need ux derivs from across boundaries)");
        }
        if (i_qcr2000 == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120, "      using QCR2000");
        }
        if (i_qcr2013 == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120, "      using QCR2013 (limited)");
        }
        if (i_qcr2013v == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120, "      using QCR2013-V");
        }
        if (ivisc[1] == 11 || ivisc[1] == 12 || ivisc[1] == 13 ||
            ivisc[1] == 14 || i_nonlin != 0) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120, "      nonlinear constitutive relation");
        }
        if (i_nonlin != 0) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "      snonlin_lim=%19.8E", (double)snonlin_lim);
        }
    }

    // store selected cell centered data in t array
    // t(1,22): density, t(1,18): u, t(1,19): v, t(1,20): w, t(1,16): c*c/(gm1*pr)
    l1 = jdim * kdim - 1;
    for (ipl = 1; ipl <= npl; ipl++) {
        ii = i + ipl - 1;
        for (j = 1; j <= jdim1; j++) {
            jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
            for (k = 1; k <= kdim1; k++) {
                t(k + jkv, 22) = q(j, k, ii, 1);
                t(k + jkv, 18) = q(j, k, ii, 2);
                t(k + jkv, 19) = q(j, k, ii, 3);
                t(k + jkv, 20) = q(j, k, ii, 4);
                t(k + jkv, 16) = q(j, k, ii, 5);
            }
        }
    }

    // store dependent variables on j=0 and j=jdim boundaries in wj0
    for (m = 0; m <= 10; m += 10) {
        if (m == 0) {
            mm = 1;
            mb = 1;
        } else {
            mm = 3;
            mb = 2;
        }
        imin = i;
        imax = i + npl - 1;
        kmin = 1;
        kmax = kdim - 1;

        for (ii = imin; ii <= imax; ii++) {
            kc = (ii - i) * kdim1 + kmin - 1;
            for (int kk = kmin; kk <= kmax; kk++) {
                kc = kc + 1;
                wj0(kc, m + 12) = qj0(kk, ii, 1, mm);
                wj0(kc, m + 8)  = qj0(kk, ii, 2, mm);
                wj0(kc, m + 9)  = qj0(kk, ii, 3, mm);
                wj0(kc, m + 10) = qj0(kk, ii, 4, mm);
                wj0(kc, m + 6)  = qj0(kk, ii, 5, mm);
                // m+5 flag indicates whether ghost cell or interface values stored in wj0
                wj0(kc, m + 5)  = bcj(kk, ii, mb);
            }
        }
    }

    // t(17)/wj0(m+7): (u*u+v*v+w*w)/2
    // t(1)/wj0(m+1):  1/density
    // t(16)/wj0(m+6): c*c/(gm1*pr)
    for (izz = 1; izz <= n; izz++) {
        t(izz, 17) = 0.5e0 * (t(izz, 18) * t(izz, 18)
                             + t(izz, 19) * t(izz, 19)
                             + t(izz, 20) * t(izz, 20));
        t(izz, 1)  = 1.e0 / t(izz, 22);
        t(izz, 16) = gpr * t(izz, 16) * t(izz, 1) / (double)gm1;
    }
    for (m = 0; m <= 10; m += 10) {
        for (izz = 1; izz <= kv; izz++) {
            wj0(izz, m + 7) = 0.5e0 * (wj0(izz, m + 8)  * wj0(izz, m + 8)
                                      + wj0(izz, m + 9)  * wj0(izz, m + 9)
                                      + wj0(izz, m + 10) * wj0(izz, m + 10));
            wj0(izz, m + 1) = 1.e0 / wj0(izz, m + 12);
            wj0(izz, m + 6) = gpr * wj0(izz, m + 6) * wj0(izz, m + 1) / (double)gm1;
        }
    }



    // t(7): laminar viscosity at cell centers (via sutherland relation)
    c2b  = (double)cbar / (double)tinf;
    c2bp = c2b + 1.e0;

    for (izz = 1; izz <= n; izz++) {
        t5       = gm1pr * t(izz, 16);
        t6       = std::sqrt(t5);
        t(izz, 7) = c2bp * t5 * t6 / (c2b + t5);
    }

    // t(14): laminar viscosity values at cell interfaces
    // interior interfaces
    for (izz = 1; izz <= nn; izz++) {
        t(izz + ks - 1, 14) = (t(izz, 7) + t(izz + ks - 1, 7)) * 0.5e0;
    }

    // j=0 and j=jdim interfaces
    for (izz = 1; izz <= kv; izz++) {
        ab   = 1.0 + wj0(izz, 5);
        bb   = 1.0 - wj0(izz, 5);
        wj05 = gm1pr * 0.5 * (ab * wj0(izz, 6) + bb * t(izz, 16));
        wj06 = std::sqrt(wj05);
        t(izz, 14)     = c2bp * wj05 * wj06 / (c2b + wj05);
        ab   = 1.0 + wj0(izz, 15);
        bb   = 1.0 - wj0(izz, 15);
        wj05 = gm1pr * 0.5 * (ab * wj0(izz, 16) + bb * t(izz + n - kv, 16));
        wj06 = std::sqrt(wj05);
        t(izz + n, 14) = c2bp * wj05 * wj06 / (c2b + wj05);
    }

    // t(15): average jacobian (inverse volume) at cell interface
    n1 = kdim * jdim1 + 1;
    l1 = kdim * jdim - 1;
    for (ipl = 1; ipl <= npl; ipl++) {
        ii = i + ipl - 1;
        for (j = 1; j <= jdim; j++) {
            jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
            if (j == 1) {
                // inverse volume at j=0 interface
                for (k = 1; k <= kdim1; k++) {
                    t(k + jkv, 15) = 2.0 / (volj0(k, ii, 1) + vol(1, k, ii));
                }
            } else if (j == jdim) {
                // inverse volume at j=jdim interface
                for (k = 1; k <= kdim1; k++) {
                    t(k + jkv, 15) = 2.0 / (volj0(k, ii, 3) + vol(jdim1, k, ii));
                }
            } else {
                // inverse volume at interior interfaces
                for (k = 1; k <= kdim1; k++) {
                    t(k + jkv, 15) = 2.0 / (vol(j, k, ii) + vol(j - 1, k, ii));
                }
            }

            // t(25)-t(27): components of grad(eta)
            for (k = 1; k <= kdim1; k++) {
                t1           = sj(j, k, ii, 4) * t(k + jkv, 15);
                t(k + jkv, 25) = sj(j, k, ii, 1) * t1;
                t(k + jkv, 26) = sj(j, k, ii, 2) * t1;
                t(k + jkv, 27) = sj(j, k, ii, 3) * t1;
            }
        }
    }

    // t(6)-t(10): gradients at cell interfaces
    for (l = 1; l <= 5; l++) {
        int k1 = l + 5;
        int k2 = l + 15;
        // interior interfaces
        for (izz = 1; izz <= nn; izz++) {
            t(izz + ks - 1, k1) = t(izz + ks - 1, k2) - t(izz, k2);
        }
        // interfaces at j=0/jdim
        for (izz = 1; izz <= kv; izz++) {
            double z1 = 1.0 + wj0(izz, 5);
            double z2 = 1.0 + wj0(izz, 15);
            t(izz,     k1) = z1 * (t(izz, k2) - wj0(izz, k1));
            t(izz + n, k1) = z2 * (wj0(izz, k2) - t(izz + n - kv, k2));
        }
    }

    // t(24): turbulent viscosity at interfaces (=0 for laminar flow)
    iviscc = ivisc[1];
    if (iviscc > 1) {
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            for (j = 1; j <= jdim; j++) {
                jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
                if (j == 1) {
                    // j=0 interface
                    for (k = 1; k <= kdim1; k++) {
                        t(k + jkv, 24) = bcj(k, ii, 1) * vj0(k, ii, 1, 1) +
                            (1.0 - bcj(k, ii, 1)) * 0.5 * (vj0(k, ii, 1, 1) + vist3d(1, k, ii));
                    }
                } else if (j == jdim) {
                    // j=jdim interface
                    for (k = 1; k <= kdim1; k++) {
                        t(k + jkv, 24) = bcj(k, ii, 2) * vj0(k, ii, 1, 3) +
                            (1.0 - bcj(k, ii, 2)) * 0.5 * (vj0(k, ii, 1, 3) + vist3d(jdim1, k, ii));
                    }
                } else {
                    // interior interfaces
                    for (k = 1; k <= kdim1; k++) {
                        t(k + jkv, 24) = 0.5 * (vist3d(j, k, ii) + vist3d(j - 1, k, ii));
                    }
                }
            }
        }
    } else {
        // laminar
        for (izz = 1; izz <= l0; izz++) {
            t(izz, 24) = 0.0;
        }
    }

    for (izz = 1; izz <= l0; izz++) {
        // ratio of turbulent to laminar viscosity
        t25 = t(izz, 24) / t(izz, 14);
        t24 = 1.e0 + t25 * coef_eddy;
        // t(23): ratio of laminar Prandtl number to total Prandtl number
        t(izz, 23) = (1.e0 + prtr * t25) / t24;
        // t(14): (mach/re)*viscosity/J
        t(izz, 14) = xmre * t24 * t(izz, 14) / t(izz, 15);
        // t(5): t(14)*(grad(zeta))**2
        t(izz, 5)  = t(izz, 14) * (t(izz, 25) * t(izz, 25)
                                  + t(izz, 26) * t(izz, 26)
                                  + t(izz, 27) * t(izz, 27));
    }

    // t(11): t(14)*((d(eta)/dx * d(u)/d(eta)) + ...)/3
    tr1 = 1.e0 / 3.e0;
    for (izz = 1; izz <= l0; izz++) {
        t(izz, 11) = tr1 * t(izz, 14) * (t(izz, 25) * t(izz, 8)
                                        + t(izz, 26) * t(izz, 9)
                                        + t(izz, 27) * t(izz, 10));
    }

    // store off vmuj: viscosity coefficient (laminar + turbulent) on j=0 and j=jdim
    for (ipl = 1; ipl <= npl; ipl++) {
        ii = i + ipl - 1;
        for (k = 1; k <= kdim1; k++) {
            izz = (ipl - 1) * kdim1 + k;
            vmuj(k, ii, 1) = t(izz,     14) * t(izz,     15) / xmre;
            vmuj(k, ii, 2) = t(izz + n, 14) * t(izz + n, 15) / xmre;
        }
    }

    if (iadv < 0) return;



    // t(12): contravariant velocity V at cell interfaces
    // interior interfaces
    for (izz = 1; izz <= nn; izz++) {
        t(izz + ks - 1, 12) = 0.5e0 * (
            t(izz + ks - 1, 25) * (t(izz + ks - 1, 18) + t(izz, 18)) +
            t(izz + ks - 1, 26) * (t(izz + ks - 1, 19) + t(izz, 19)) +
            t(izz + ks - 1, 27) * (t(izz + ks - 1, 20) + t(izz, 20)));
    }
    // interfaces at j=0 and j=jdim
    for (izz = 1; izz <= kv; izz++) {
        ab = 1.0 + wj0(izz, 15);
        bb = 1.0 - wj0(izz, 15);
        t(izz + n, 12) = 0.5 * (
            t(izz + n, 25) * (ab * wj0(izz, 18) + bb * t(izz + n - kv, 18)) +
            t(izz + n, 26) * (ab * wj0(izz, 19) + bb * t(izz + n - kv, 19)) +
            t(izz + n, 27) * (ab * wj0(izz, 20) + bb * t(izz + n - kv, 20)));
        ab = 1.0 + wj0(izz, 5);
        bb = 1.0 - wj0(izz, 5);
        t(izz, 12) = 0.5 * (
            t(izz, 25) * (ab * wj0(izz, 8)  + bb * t(izz, 18)) +
            t(izz, 26) * (ab * wj0(izz, 9)  + bb * t(izz, 19)) +
            t(izz, 27) * (ab * wj0(izz, 10) + bb * t(izz, 20)));
    }

    // calculate fluxes
    for (l = 2; l <= 4; l++) {
        // viscous terms at interfaces (momentum eqs.)
        for (izz = 1; izz <= l0; izz++) {
            t(izz, 16) = t(izz, 5) * t(izz, l + 6) + t(izz, l + 23) * t(izz, 11);
        }
        // (-)viscous flux = Gv(k-1/2)-Gv(k+1/2) (momentum eqs.)
        for (izz = 1; izz <= n; izz++) {
            t(izz, l) = -t(izz + ks - 1, 16) + t(izz, 16);
        }
    }

    // viscous terms at interfaces (energy eq.)
    for (izz = 1; izz <= l0; izz++) {
        t(izz, 16) = t(izz, 5) * (t(izz, 23) * t(izz, 6) + t(izz, 7))
                   + t(izz, 12) * t(izz, 11);
    }
    // (-)viscous flux = Gv(k-1/2)-Gv(k+1/2) (energy eq.)
    for (izz = 1; izz <= n; izz++) {
        t(izz, 5) = -t(izz + ks - 1, 16) + t(izz, 16);
    }

    // calculate residuals
    for (ipl = 1; ipl <= npl; ipl++) {
        ii  = i + ipl - 1;
        for (j = 1; j <= jdim1; j++) {
            jkv = (j - 1) * kdim1 * npl + (ipl - 1) * kdim1;
            for (l = 2; l <= 5; l++) {
                for (k = 1; k <= kdim1; k++) {
                    res(j, k, ii, l) = res(j, k, ii, l) + t(k + jkv, l);
                }
            }
        }
    }



    // Add tau_ij's for general method (i_tauijs==1 branch):
    if (i_tauijs == 1 && ivisc[1] < 70) {
        re = (double)reue / (double)xmach;
        // Add nonlinear terms to RHS
        for (ipl = 1; ipl <= npl; ipl++) {
            ii   = i + ipl - 1;
            jkv2 = (ipl - 1) * kdim1;
            for (j = 1; j <= jdim; j++) {
                jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
                jp1 = std::min(j + 1, jdim - 1);
                jm2 = std::max(j - 2, 1);
                for (k = 1; k <= kdim1; k++) {
                    if (j == 1) {
                        rho = 0.5 * (q(j, k, ii, 1) + qj0(k, ii, 1, 1)) * (1.0 - wj0(k + jkv2, 5)) +
                              qj0(k, ii, 1, 1) * wj0(k + jkv2, 5);
                        u   = 0.5 * (q(j, k, ii, 2) + qj0(k, ii, 2, 1)) * (1.0 - wj0(k + jkv2, 5)) +
                              qj0(k, ii, 2, 1) * wj0(k + jkv2, 5);
                        v   = 0.5 * (q(j, k, ii, 3) + qj0(k, ii, 3, 1)) * (1.0 - wj0(k + jkv2, 5)) +
                              qj0(k, ii, 3, 1) * wj0(k + jkv2, 5);
                        w_vel = 0.5 * (q(j, k, ii, 4) + qj0(k, ii, 4, 1)) * (1.0 - wj0(k + jkv2, 5)) +
                              qj0(k, ii, 4, 1) * wj0(k + jkv2, 5);
                        v3d = 0.5 * (vist3d(j, k, ii) + vj0(k, ii, 1, 1)) * (1.0 - wj0(k + jkv2, 5)) +
                              vj0(k, ii, 1, 1) * wj0(k + jkv2, 5);
                        // use 1st order interpolation at end
                        ux1 = 1.5 * ux(j, k, ii, 1) - 0.5 * ux(jp1, k, ii, 1);
                        ux2 = 1.5 * ux(j, k, ii, 2) - 0.5 * ux(jp1, k, ii, 2);
                        ux3 = 1.5 * ux(j, k, ii, 3) - 0.5 * ux(jp1, k, ii, 3);
                        ux4 = 1.5 * ux(j, k, ii, 4) - 0.5 * ux(jp1, k, ii, 4);
                        ux5 = 1.5 * ux(j, k, ii, 5) - 0.5 * ux(jp1, k, ii, 5);
                        ux6 = 1.5 * ux(j, k, ii, 6) - 0.5 * ux(jp1, k, ii, 6);
                        ux7 = 1.5 * ux(j, k, ii, 7) - 0.5 * ux(jp1, k, ii, 7);
                        ux8 = 1.5 * ux(j, k, ii, 8) - 0.5 * ux(jp1, k, ii, 8);
                        ux9 = 1.5 * ux(j, k, ii, 9) - 0.5 * ux(jp1, k, ii, 9);
                    } else if (j == jdim) {
                        rho = 0.5 * (q(j - 1, k, ii, 1) + qj0(k, ii, 1, 3)) * (1.0 - wj0(k + jkv2, 15)) +
                              qj0(k, ii, 1, 3) * wj0(k + jkv2, 15);
                        u   = 0.5 * (q(j - 1, k, ii, 2) + qj0(k, ii, 2, 3)) * (1.0 - wj0(k + jkv2, 15)) +
                              qj0(k, ii, 2, 3) * wj0(k + jkv2, 15);
                        v   = 0.5 * (q(j - 1, k, ii, 3) + qj0(k, ii, 3, 3)) * (1.0 - wj0(k + jkv2, 15)) +
                              qj0(k, ii, 3, 3) * wj0(k + jkv2, 15);
                        w_vel = 0.5 * (q(j - 1, k, ii, 4) + qj0(k, ii, 4, 3)) * (1.0 - wj0(k + jkv2, 15)) +
                              qj0(k, ii, 4, 3) * wj0(k + jkv2, 15);
                        v3d = 0.5 * (vist3d(j - 1, k, ii) + vj0(k, ii, 1, 3)) * (1.0 - wj0(k + jkv2, 15)) +
                              vj0(k, ii, 1, 3) * wj0(k + jkv2, 15);
                        // use 1st order interpolation at end
                        ux1 = 1.5 * ux(j - 1, k, ii, 1) - 0.5 * ux(jm2, k, ii, 1);
                        ux2 = 1.5 * ux(j - 1, k, ii, 2) - 0.5 * ux(jm2, k, ii, 2);
                        ux3 = 1.5 * ux(j - 1, k, ii, 3) - 0.5 * ux(jm2, k, ii, 3);
                        ux4 = 1.5 * ux(j - 1, k, ii, 4) - 0.5 * ux(jm2, k, ii, 4);
                        ux5 = 1.5 * ux(j - 1, k, ii, 5) - 0.5 * ux(jm2, k, ii, 5);
                        ux6 = 1.5 * ux(j - 1, k, ii, 6) - 0.5 * ux(jm2, k, ii, 6);
                        ux7 = 1.5 * ux(j - 1, k, ii, 7) - 0.5 * ux(jm2, k, ii, 7);
                        ux8 = 1.5 * ux(j - 1, k, ii, 8) - 0.5 * ux(jm2, k, ii, 8);
                        ux9 = 1.5 * ux(j - 1, k, ii, 9) - 0.5 * ux(jm2, k, ii, 9);
                    } else {
                        rho = 0.5 * (q(j, k, ii, 1) + q(j - 1, k, ii, 1));
                        u   = 0.5 * (q(j, k, ii, 2) + q(j - 1, k, ii, 2));
                        v   = 0.5 * (q(j, k, ii, 3) + q(j - 1, k, ii, 3));
                        w_vel = 0.5 * (q(j, k, ii, 4) + q(j - 1, k, ii, 4));
                        v3d = 0.5 * (vist3d(j, k, ii) + vist3d(j - 1, k, ii));
                        ux1 = 0.5 * (ux(j, k, ii, 1) + ux(j - 1, k, ii, 1));
                        ux2 = 0.5 * (ux(j, k, ii, 2) + ux(j - 1, k, ii, 2));
                        ux3 = 0.5 * (ux(j, k, ii, 3) + ux(j - 1, k, ii, 3));
                        ux4 = 0.5 * (ux(j, k, ii, 4) + ux(j - 1, k, ii, 4));
                        ux5 = 0.5 * (ux(j, k, ii, 5) + ux(j - 1, k, ii, 5));
                        ux6 = 0.5 * (ux(j, k, ii, 6) + ux(j - 1, k, ii, 6));
                        ux7 = 0.5 * (ux(j, k, ii, 7) + ux(j - 1, k, ii, 7));
                        ux8 = 0.5 * (ux(j, k, ii, 8) + ux(j - 1, k, ii, 8));
                        ux9 = 0.5 * (ux(j, k, ii, 9) + ux(j - 1, k, ii, 9));
                    }
                    s11 = ux1;
                    s22 = ux5;
                    s33 = ux9;
                    s11 = s11 - ((s11 + s22 + s33) / 3.0);
                    s22 = ux5 - ((ux1 + ux5 + ux9) / 3.0);
                    s33 = ux9 - ((ux1 + ux5 + ux9) / 3.0);
                    // recompute s11star properly
                    {
                        double ss11 = ux1, ss22 = ux5, ss33 = ux9;
                        double trace3 = (ss11 + ss22 + ss33) / 3.0;
                        s11 = ss11 - trace3;
                        s22 = ss22 - trace3;
                        s33 = ss33 - trace3;
                    }
                    s12 = 0.5 * (ux2 + ux4);
                    s13 = 0.5 * (ux3 + ux7);
                    s23 = 0.5 * (ux6 + ux8);
                    w12 = 0.5 * (ux2 - ux4);
                    w13 = 0.5 * (ux3 - ux7);
                    w23 = 0.5 * (ux6 - ux8);
                    // ignoring -2/3*rho*k*delij term
                    t11 = -2.0 * v3d * s11;
                    t22 = -2.0 * v3d * s22;
                    t33 = -2.0 * v3d * s33;
                    t12 = -2.0 * v3d * s12;
                    t13 = -2.0 * v3d * s13;
                    t23 = -2.0 * v3d * s23;
                    if (i_qcr2000 == 1 || i_qcr2013 == 1 || i_qcr2013v == 1) {
                        denom = std::sqrt(ux1*ux1 + ux2*ux2 + ux3*ux3 +
                                          ux4*ux4 + ux5*ux5 + ux6*ux6 +
                                          ux7*ux7 + ux8*ux8 + ux9*ux9) + 1.e-20;
                        o11 = 0.0;
                        o22 = 0.0;
                        o33 = 0.0;
                        o12 = 2.0 * w12 / denom;
                        o13 = 2.0 * w13 / denom;
                        o23 = 2.0 * w23 / denom;
                        t11p = t11 - 0.3 * (o11*t11 + o12*t12 + o13*t13 +
                                             o11*t11 + o12*t12 + o13*t13);
                        t22p = t22 - 0.3 * (-o12*t12 + o22*t22 + o23*t23 +
                                             -o12*t12 + o22*t22 + o23*t23);
                        t33p = t33 - 0.3 * (-o13*t13 - o23*t23 + o33*t33 +
                                             -o13*t13 - o23*t23 + o33*t33);
                        t12p = t12 - 0.3 * (o11*t12 + o12*t22 + o13*t23 +
                                             -o12*t11 + o22*t12 + o23*t13);
                        t13p = t13 - 0.3 * (o11*t13 + o12*t23 + o13*t33 +
                                             -o13*t11 - o23*t12 + o33*t13);
                        t23p = t23 - 0.3 * (-o12*t13 + o22*t23 + o23*t33 +
                                             -o13*t12 - o23*t22 + o33*t23);
                        t11 = t11p;
                        t22 = t22p;
                        t33 = t33p;
                        t12 = t12p;
                        t13 = t13p;
                        t23 = t23p;
                    }
                    if (i_qcr2013 == 1) {
                        xis2013 = s11*s11 + s22*s22 + s33*s33 +
                                  2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        s_mod = std::sqrt(2.0 * xis2013);
                        wis   = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        w_mod = std::sqrt(2.0 * wis);
                        { double a1 = s_mod, a2 = 1.2 * w_mod; s_mod = ccomplex_ns::ccmin(a1, a2); }
                        t11 = t11 + ccr2 * v3d * s_mod;
                        t22 = t22 + ccr2 * v3d * s_mod;
                        t33 = t33 + ccr2 * v3d * s_mod;
                    }
                    if (i_qcr2013v == 1) {
                        wis   = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        w_mod = std::sqrt(2.0 * wis);
                        t11 = t11 + ccr2 * v3d * w_mod;
                        t22 = t22 + ccr2 * v3d * w_mod;
                        t33 = t33 + ccr2 * v3d * w_mod;
                    }
                    t(k + jkv, 28) = xmre * (t(k + jkv, 25)*t11 + t(k + jkv, 26)*t12 +
                                              t(k + jkv, 27)*t13) / t(k + jkv, 15);
                    t(k + jkv, 29) = xmre * (t(k + jkv, 25)*t12 + t(k + jkv, 26)*t22 +
                                              t(k + jkv, 27)*t23) / t(k + jkv, 15);
                    t(k + jkv, 30) = xmre * (t(k + jkv, 25)*t13 + t(k + jkv, 26)*t23 +
                                              t(k + jkv, 27)*t33) / t(k + jkv, 15);
                    t(k + jkv, 31) = xmre / t(k + jkv, 15) * (
                        t(k + jkv, 25) * (u * t11 + v * t12 + w_vel * t13) +
                        t(k + jkv, 26) * (u * t12 + v * t22 + w_vel * t23) +
                        t(k + jkv, 27) * (u * t13 + v * t23 + w_vel * t33));
                } // k loop
            } // j loop
        } // ipl loop
        // Take positive derivative for correct sign
        for (l = 2; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                t(izz, l) = t(izz + ks - 1, l + 26) - t(izz, l + 26);
            }
        }
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            for (j = 1; j <= jdim1; j++) {
                jkv = (j - 1) * kdim1 * npl + (ipl - 1) * kdim1 + 1;
                for (l = 2; l <= 5; l++) {
                    for (k = 1; k <= kdim1; k++) {
                        res(j, k, ii, l) = res(j, k, ii, l) + t(k + jkv - 1, l);
                    }
                }
            }
        }
    } // end i_tauijs==1 branch



    // Add tau_ij's for RSMs (ivisc(2) >= 70):
    if (ivisc[1] >= 70) {
        re = (double)reue / (double)xmach;
        // Add nonlinear terms to RHS
        for (ipl = 1; ipl <= npl; ipl++) {
            ii   = i + ipl - 1;
            jkv2 = (ipl - 1) * kdim1;
            for (j = 1; j <= jdim; j++) {
                jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
                if (j == 1) {
                    for (k = 1; k <= kdim1; k++) {
                        rho = 0.5 * (q(j, k, ii, 1) + qj0(k, ii, 1, 1)) * (1.0 - wj0(k + jkv2, 5)) +
                              qj0(k, ii, 1, 1) * wj0(k + jkv2, 5);
                        u   = 0.5 * (q(j, k, ii, 2) + qj0(k, ii, 2, 1)) * (1.0 - wj0(k + jkv2, 5)) +
                              qj0(k, ii, 2, 1) * wj0(k + jkv2, 5);
                        v   = 0.5 * (q(j, k, ii, 3) + qj0(k, ii, 3, 1)) * (1.0 - wj0(k + jkv2, 5)) +
                              qj0(k, ii, 3, 1) * wj0(k + jkv2, 5);
                        w_vel = 0.5 * (q(j, k, ii, 4) + qj0(k, ii, 4, 1)) * (1.0 - wj0(k + jkv2, 5)) +
                              qj0(k, ii, 4, 1) * wj0(k + jkv2, 5);
                        t11 = -rho * 0.5 * (zksav(j, k, ii, 1) + tj0(k, ii, 1, 1)) * re;
                        t22 = -rho * 0.5 * (zksav(j, k, ii, 2) + tj0(k, ii, 2, 1)) * re;
                        t33 = -rho * 0.5 * (zksav(j, k, ii, 3) + tj0(k, ii, 3, 1)) * re;
                        t12 = -rho * 0.5 * (zksav(j, k, ii, 4) + tj0(k, ii, 4, 1)) * re;
                        t23 = -rho * 0.5 * (zksav(j, k, ii, 5) + tj0(k, ii, 5, 1)) * re;
                        t13 = -rho * 0.5 * (zksav(j, k, ii, 6) + tj0(k, ii, 6, 1)) * re;
                        t(k + jkv, 28) = xmre * (t(k + jkv, 25)*t11 + t(k + jkv, 26)*t12 +
                                                   t(k + jkv, 27)*t13) / t(k + jkv, 15);
                        t(k + jkv, 29) = xmre * (t(k + jkv, 25)*t12 + t(k + jkv, 26)*t22 +
                                                   t(k + jkv, 27)*t23) / t(k + jkv, 15);
                        t(k + jkv, 30) = xmre * (t(k + jkv, 25)*t13 + t(k + jkv, 26)*t23 +
                                                   t(k + jkv, 27)*t33) / t(k + jkv, 15);
                        t(k + jkv, 31) = xmre / t(k + jkv, 15) * (
                            t(k + jkv, 25) * (u * t11 + v * t12 + w_vel * t13) +
                            t(k + jkv, 26) * (u * t12 + v * t22 + w_vel * t23) +
                            t(k + jkv, 27) * (u * t13 + v * t23 + w_vel * t33));
                    }
                } else if (j == jdim) {
                    for (k = 1; k <= kdim1; k++) {
                        rho = 0.5 * (q(j - 1, k, ii, 1) + qj0(k, ii, 1, 3)) * (1.0 - wj0(k + jkv2, 15)) +
                              qj0(k, ii, 1, 3) * wj0(k + jkv2, 15);
                        u   = 0.5 * (q(j - 1, k, ii, 2) + qj0(k, ii, 2, 3)) * (1.0 - wj0(k + jkv2, 15)) +
                              qj0(k, ii, 2, 3) * wj0(k + jkv2, 15);
                        v   = 0.5 * (q(j - 1, k, ii, 3) + qj0(k, ii, 3, 3)) * (1.0 - wj0(k + jkv2, 15)) +
                              qj0(k, ii, 3, 3) * wj0(k + jkv2, 15);
                        w_vel = 0.5 * (q(j - 1, k, ii, 4) + qj0(k, ii, 4, 3)) * (1.0 - wj0(k + jkv2, 15)) +
                              qj0(k, ii, 4, 3) * wj0(k + jkv2, 15);
                        t11 = -rho * 0.5 * (zksav(j - 1, k, ii, 1) + tj0(k, ii, 1, 3)) * re;
                        t22 = -rho * 0.5 * (zksav(j - 1, k, ii, 2) + tj0(k, ii, 2, 3)) * re;
                        t33 = -rho * 0.5 * (zksav(j - 1, k, ii, 3) + tj0(k, ii, 3, 3)) * re;
                        t12 = -rho * 0.5 * (zksav(j - 1, k, ii, 4) + tj0(k, ii, 4, 3)) * re;
                        t23 = -rho * 0.5 * (zksav(j - 1, k, ii, 5) + tj0(k, ii, 5, 3)) * re;
                        t13 = -rho * 0.5 * (zksav(j - 1, k, ii, 6) + tj0(k, ii, 6, 3)) * re;
                        t(k + jkv, 28) = xmre * (t(k + jkv, 25)*t11 + t(k + jkv, 26)*t12 +
                                                   t(k + jkv, 27)*t13) / t(k + jkv, 15);
                        t(k + jkv, 29) = xmre * (t(k + jkv, 25)*t12 + t(k + jkv, 26)*t22 +
                                                   t(k + jkv, 27)*t23) / t(k + jkv, 15);
                        t(k + jkv, 30) = xmre * (t(k + jkv, 25)*t13 + t(k + jkv, 26)*t23 +
                                                   t(k + jkv, 27)*t33) / t(k + jkv, 15);
                        t(k + jkv, 31) = xmre / t(k + jkv, 15) * (
                            t(k + jkv, 25) * (u * t11 + v * t12 + w_vel * t13) +
                            t(k + jkv, 26) * (u * t12 + v * t22 + w_vel * t23) +
                            t(k + jkv, 27) * (u * t13 + v * t23 + w_vel * t33));
                    }
                } else {
                    for (k = 1; k <= kdim1; k++) {
                        rho = 0.5 * (q(j, k, ii, 1) + q(j - 1, k, ii, 1));
                        t11 = -rho * 0.5 * (zksav(j, k, ii, 1) + zksav(j - 1, k, ii, 1)) * re;
                        t22 = -rho * 0.5 * (zksav(j, k, ii, 2) + zksav(j - 1, k, ii, 2)) * re;
                        t33 = -rho * 0.5 * (zksav(j, k, ii, 3) + zksav(j - 1, k, ii, 3)) * re;
                        t12 = -rho * 0.5 * (zksav(j, k, ii, 4) + zksav(j - 1, k, ii, 4)) * re;
                        t23 = -rho * 0.5 * (zksav(j, k, ii, 5) + zksav(j - 1, k, ii, 5)) * re;
                        t13 = -rho * 0.5 * (zksav(j, k, ii, 6) + zksav(j - 1, k, ii, 6)) * re;
                        t(k + jkv, 28) = xmre * (t(k + jkv, 25)*t11 + t(k + jkv, 26)*t12 +
                                                   t(k + jkv, 27)*t13) / t(k + jkv, 15);
                        t(k + jkv, 29) = xmre * (t(k + jkv, 25)*t12 + t(k + jkv, 26)*t22 +
                                                   t(k + jkv, 27)*t23) / t(k + jkv, 15);
                        t(k + jkv, 30) = xmre * (t(k + jkv, 25)*t13 + t(k + jkv, 26)*t23 +
                                                   t(k + jkv, 27)*t33) / t(k + jkv, 15);
                        t(k + jkv, 31) = xmre / t(k + jkv, 15) * (
                            t(k + jkv, 25) * (0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t11
                                             + 0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t12
                                             + 0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t13) +
                            t(k + jkv, 26) * (0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t12
                                             + 0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t22
                                             + 0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t23) +
                            t(k + jkv, 27) * (0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t13
                                             + 0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t23
                                             + 0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t33));
                    }
                }
            } // j loop
        } // ipl loop
        // Take positive derivative for correct sign
        for (l = 2; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                t(izz, l) = t(izz + ks - 1, l + 26) - t(izz, l + 26);
            }
        }
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            for (j = 1; j <= jdim1; j++) {
                jkv = (j - 1) * kdim1 * npl + (ipl - 1) * kdim1 + 1;
                for (l = 2; l <= 5; l++) {
                    for (k = 1; k <= kdim1; k++) {
                        res(j, k, ii, l) = res(j, k, ii, l) + t(k + jkv - 1, l);
                    }
                }
            }
        }
    } // end ivisc(2)>=70 branch



    // EASM v1 branch: ivisc(2)==11 or 12, i_tauijs!=1
    if (i_tauijs != 1 && (ivisc[1] == 11 || ivisc[1] == 12) && level >= lglobal) {
        // f7=factor for k-omega (1) or k-epsilon (0)
        f7 = 0.0;
        if (ivisc[1] == 12) f7 = 1.0;
        re = (double)reue / (double)xmach;
        gg = 1.0 / ((double)c10 + (double)c5_e - 1.0);
        // Add nonlinear terms to RHS
        for (ipl = 1; ipl <= npl; ipl++) {
            ii   = i + ipl - 1;
            jkv2 = (ipl - 1) * kdim1;
            for (j = 1; j <= jdim; j++) {
                jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
                if (j == 1) {
                    for (k = 1; k <= kdim1; k++) {
                        s11 = t(k+jkv,8)*t(k+jkv,25);
                        s22 = t(k+jkv,9)*t(k+jkv,26);
                        s33 = t(k+jkv,10)*t(k+jkv,27);
                        s12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)+t(k+jkv,9)*t(k+jkv,25));
                        s13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,25));
                        s23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,26));
                        w12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)-t(k+jkv,9)*t(k+jkv,25));
                        w13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,25));
                        w23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,26));
                        xis = s11*s11+s22*s22+s33*s33+2.0*s12*s12+2.0*s13*s13+2.0*s23*s23;
                        wis = 2.0*w12*w12+2.0*w13*w13+2.0*w23*w23;
                        alpa1 = (2.0-(double)c4_e)/2.0*gg;
                        alpa2 = (2.0-(double)c3_e)*gg;
                        rho = 0.5*(q(j,k,ii,1)+qj0(k,ii,1,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,1,1)*wj0(k+jkv2,5);
                        u   = 0.5*(q(j,k,ii,2)+qj0(k,ii,2,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,2,1)*wj0(k+jkv2,5);
                        v   = 0.5*(q(j,k,ii,3)+qj0(k,ii,3,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,3,1)*wj0(k+jkv2,5);
                        w_vel=0.5*(q(j,k,ii,4)+qj0(k,ii,4,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,4,1)*wj0(k+jkv2,5);
                        zksav1 = 0.5*(zksav(j,k,ii,1)+tj0(k,ii,1,1));
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        zksav2 = 0.5*(zksav(j,k,ii,2)+tj0(k,ii,2,1));
                        bnum = zksav2*(1.0-f7)+f7;
                        eta   = (2.0-(double)c3_e)*(2.0-(double)c3_e)*(gg*gg/4.0)*xis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        squig = (2.0-(double)c4_e)*(2.0-(double)c4_e)*(gg*gg/4.0)*wis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        { double a1=eta,   a2=10.0; eta   = ccomplex_ns::ccmincr(a1,a2); }
                        { double a1=squig, a2=10.0; squig = ccomplex_ns::ccmincr(a1,a2); }
                        factre = (3.0*(1.0+eta)+0.2e-8*(eta*eta*eta+squig*squig*squig)) /
                                 (3.0*(1.0+eta)+0.2   *(eta*eta*eta+squig*squig*squig));
                        t11 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(-s12*w12-s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t11 = t11*(1.0-bcj(k,ii,1));
                        t22 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s12*w12-s23*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t22 = t22*(1.0-bcj(k,ii,1));
                        t33 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s23*w23+s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t33 = t33*(1.0-bcj(k,ii,1));
                        t12 = -2.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t12 = t12*(1.0-bcj(k,ii,1));
                        t13 = -2.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t13 = t13*(1.0-bcj(k,ii,1));
                        t23 = -2.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t23 = t23*(1.0-bcj(k,ii,1));
                        t(k+jkv,28) = xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+t(k+jkv,27)*t13)/t(k+jkv,15);
                        t(k+jkv,29) = xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+t(k+jkv,27)*t23)/t(k+jkv,15);
                        t(k+jkv,30) = xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+t(k+jkv,27)*t33)/t(k+jkv,15);
                        t(k+jkv,31) = xmre/t(k+jkv,15)*(
                            t(k+jkv,25)*(u*t11+v*t12+w_vel*t13)+
                            t(k+jkv,26)*(u*t12+v*t22+w_vel*t23)+
                            t(k+jkv,27)*(u*t13+v*t23+w_vel*t33))
                            - xmre/t(k+jkv,15)*(t(k+jkv,14)*t(k+jkv,15)/xmre-t(k+jkv,24)+(double)sigk1*t(k+jkv,24))*
                              (t(k+jkv,25)*t(k+jkv,25)+t(k+jkv,26)*t(k+jkv,26)+t(k+jkv,27)*t(k+jkv,27))*
                              (zksav(j,k,ii,2)-tj0(k,ii,2,1));
                    } // k loop j==1



                } else if (j == jdim) {
                    for (k = 1; k <= kdim1; k++) {
                        s11 = t(k+jkv,8)*t(k+jkv,25);
                        s22 = t(k+jkv,9)*t(k+jkv,26);
                        s33 = t(k+jkv,10)*t(k+jkv,27);
                        s12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)+t(k+jkv,9)*t(k+jkv,25));
                        s13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,25));
                        s23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,26));
                        w12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)-t(k+jkv,9)*t(k+jkv,25));
                        w13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,25));
                        w23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,26));
                        xis = s11*s11+s22*s22+s33*s33+2.0*s12*s12+2.0*s13*s13+2.0*s23*s23;
                        wis = 2.0*w12*w12+2.0*w13*w13+2.0*w23*w23;
                        alpa1 = (2.0-(double)c4_e)/2.0*gg;
                        alpa2 = (2.0-(double)c3_e)*gg;
                        rho = 0.5*(q(j-1,k,ii,1)+qj0(k,ii,1,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,1,3)*wj0(k+jkv2,15);
                        u   = 0.5*(q(j-1,k,ii,2)+qj0(k,ii,2,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,2,3)*wj0(k+jkv2,15);
                        v   = 0.5*(q(j-1,k,ii,3)+qj0(k,ii,3,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,3,3)*wj0(k+jkv2,15);
                        w_vel=0.5*(q(j-1,k,ii,4)+qj0(k,ii,4,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,4,3)*wj0(k+jkv2,15);
                        zksav1 = 0.5*(zksav(j-1,k,ii,1)+tj0(k,ii,1,3));
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        zksav2 = 0.5*(zksav(j-1,k,ii,2)+tj0(k,ii,2,3));
                        bnum = zksav2*(1.0-f7)+f7;
                        eta   = (2.0-(double)c3_e)*(2.0-(double)c3_e)*(gg*gg/4.0)*xis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        squig = (2.0-(double)c4_e)*(2.0-(double)c4_e)*(gg*gg/4.0)*wis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        { double a1=eta,   a2=10.0; eta   = ccomplex_ns::ccmincr(a1,a2); }
                        { double a1=squig, a2=10.0; squig = ccomplex_ns::ccmincr(a1,a2); }
                        factre = (3.0*(1.0+eta)+0.2e-8*(eta*eta*eta+squig*squig*squig)) /
                                 (3.0*(1.0+eta)+0.2   *(eta*eta*eta+squig*squig*squig));
                        t11 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(-s12*w12-s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t11 = t11*(1.0-bcj(k,ii,2));
                        t22 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s12*w12-s23*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t22 = t22*(1.0-bcj(k,ii,2));
                        t33 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s23*w23+s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t33 = t33*(1.0-bcj(k,ii,2));
                        t12 = -2.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t12 = t12*(1.0-bcj(k,ii,2));
                        t13 = -2.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t13 = t13*(1.0-bcj(k,ii,2));
                        t23 = -2.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t23 = t23*(1.0-bcj(k,ii,2));
                        t(k+jkv,28) = xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+t(k+jkv,27)*t13)/t(k+jkv,15);
                        t(k+jkv,29) = xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+t(k+jkv,27)*t23)/t(k+jkv,15);
                        t(k+jkv,30) = xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+t(k+jkv,27)*t33)/t(k+jkv,15);
                        t(k+jkv,31) = xmre/t(k+jkv,15)*(
                            t(k+jkv,25)*(u*t11+v*t12+w_vel*t13)+
                            t(k+jkv,26)*(u*t12+v*t22+w_vel*t23)+
                            t(k+jkv,27)*(u*t13+v*t23+w_vel*t33))
                            - xmre/t(k+jkv,15)*(t(k+jkv,14)*t(k+jkv,15)/xmre-t(k+jkv,24)+(double)sigk1*t(k+jkv,24))*
                              (t(k+jkv,25)*t(k+jkv,25)+t(k+jkv,26)*t(k+jkv,26)+t(k+jkv,27)*t(k+jkv,27))*
                              (tj0(k,ii,2,3)-zksav(j-1,k,ii,2));
                    } // k loop j==jdim



                } else {
                    for (k = 1; k <= kdim1; k++) {
                        s11 = t(k+jkv,8)*t(k+jkv,25);
                        s22 = t(k+jkv,9)*t(k+jkv,26);
                        s33 = t(k+jkv,10)*t(k+jkv,27);
                        s12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)+t(k+jkv,9)*t(k+jkv,25));
                        s13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,25));
                        s23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,26));
                        w12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)-t(k+jkv,9)*t(k+jkv,25));
                        w13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,25));
                        w23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,26));
                        xis = s11*s11+s22*s22+s33*s33+2.0*s12*s12+2.0*s13*s13+2.0*s23*s23;
                        wis = 2.0*w12*w12+2.0*w13*w13+2.0*w23*w23;
                        alpa1 = (2.0-(double)c4_e)/2.0*gg;
                        alpa2 = (2.0-(double)c3_e)*gg;
                        rho = 0.5*(q(j,k,ii,1)+q(j-1,k,ii,1));
                        zksav1 = 0.5*(zksav(j,k,ii,1)+zksav(j-1,k,ii,1));
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        zksav2 = 0.5*(zksav(j,k,ii,2)+zksav(j-1,k,ii,2));
                        bnum = zksav2*(1.0-f7)+f7;
                        eta   = (2.0-(double)c3_e)*(2.0-(double)c3_e)*(gg*gg/4.0)*xis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        squig = (2.0-(double)c4_e)*(2.0-(double)c4_e)*(gg*gg/4.0)*wis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        { double a1=eta,   a2=10.0; eta   = ccomplex_ns::ccmincr(a1,a2); }
                        { double a1=squig, a2=10.0; squig = ccomplex_ns::ccmincr(a1,a2); }
                        factre = (3.0*(1.0+eta)+0.2e-8*(eta*eta*eta+squig*squig*squig)) /
                                 (3.0*(1.0+eta)+0.2   *(eta*eta*eta+squig*squig*squig));
                        t11 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(-s12*w12-s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t22 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s12*w12-s23*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t33 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s23*w23+s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t12 = -2.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t13 = -2.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t23 = -2.0*alpa1*t(k+jkv,24)*factre*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*factre*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t(k+jkv,28) = xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+t(k+jkv,27)*t13)/t(k+jkv,15);
                        t(k+jkv,29) = xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+t(k+jkv,27)*t23)/t(k+jkv,15);
                        t(k+jkv,30) = xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+t(k+jkv,27)*t33)/t(k+jkv,15);
                        t(k+jkv,31) = xmre/t(k+jkv,15)*(
                            t(k+jkv,25)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t11
                                        +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t12
                                        +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t13)+
                            t(k+jkv,26)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t12
                                        +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t22
                                        +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t23)+
                            t(k+jkv,27)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t13
                                        +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t23
                                        +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t33))
                            - xmre/t(k+jkv,15)*(t(k+jkv,14)*t(k+jkv,15)/xmre-t(k+jkv,24)+(double)sigk1*t(k+jkv,24))*
                              (t(k+jkv,25)*t(k+jkv,25)+t(k+jkv,26)*t(k+jkv,26)+t(k+jkv,27)*t(k+jkv,27))*
                              (zksav(j,k,ii,2)-zksav(j-1,k,ii,2));
                    } // k loop interior
                }
            } // j loop
        } // ipl loop
    } // end EASM v1 branch



    // EASM v2 branch: ivisc(2)==13 or 14, i_tauijs!=1
    else if (i_tauijs != 1 && (ivisc[1] == 13 || ivisc[1] == 14) && level >= lglobal) {
        f7 = 0.0;
        if (ivisc[1] == 14) f7 = 1.0;
        re = (double)reue / (double)xmach;
        al10 = 0.5*(double)c10 - 1.0;
        al1  = 2.0*(0.5*(double)c11 + 1.0);
        // 2-line improvement for wake
        if (ieasm_type == 0 || ieasm_type == 3 || ieasm_type == 4) {
            al10 = al10 + 1.8864;
            al1  = al1  - 2.0;
        }
        al2 = 0.5*(double)c2_e - 2.0/3.0;
        al3 = 0.5*(double)c3_e - 1.0;
        al4 = 0.5*(double)c4_e - 1.0;
        // Add nonlinear terms to RHS
        for (ipl = 1; ipl <= npl; ipl++) {
            ii   = i + ipl - 1;
            jkv2 = (ipl - 1) * kdim1;
            for (j = 1; j <= jdim; j++) {
                jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
                if (j == 1) {
                    for (k = 1; k <= kdim1; k++) {
                        s11 = t(k+jkv,8)*t(k+jkv,25);
                        s22 = t(k+jkv,9)*t(k+jkv,26);
                        s33 = t(k+jkv,10)*t(k+jkv,27);
                        s12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)+t(k+jkv,9)*t(k+jkv,25));
                        s13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,25));
                        s23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,26));
                        w12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)-t(k+jkv,9)*t(k+jkv,25));
                        w13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,25));
                        w23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,26));
                        xis = s11*s11+s22*s22+s33*s33+2.0*s12*s12+2.0*s13*s13+2.0*s23*s23;
                        // Girimaji JFM 2000 fix to c4
                        if (ieasm_type == 4) {
                            wis = 2.0*w12*w12+2.0*w13*w13+2.0*w23*w23;
                            eta1_girimaji = xis / (xis + wis);
                            if ((float)eta1_girimaji < 0.5f) {
                                c4new = 2.0 - ((2.0-(double)c4_e)*
                                    std::pow(eta1_girimaji/(1.0-eta1_girimaji), 0.75));
                            } else {
                                c4new = (double)c4_e;
                            }
                            al4 = 0.5*c4new - 1.0;
                        }
                        rho = 0.5*(q(j,k,ii,1)+qj0(k,ii,1,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,1,1)*wj0(k+jkv2,5);
                        u   = 0.5*(q(j,k,ii,2)+qj0(k,ii,2,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,2,1)*wj0(k+jkv2,5);
                        v   = 0.5*(q(j,k,ii,3)+qj0(k,ii,3,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,3,1)*wj0(k+jkv2,5);
                        w_vel=0.5*(q(j,k,ii,4)+qj0(k,ii,4,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,4,1)*wj0(k+jkv2,5);
                        zksav1 = 0.5*(zksav(j,k,ii,1)+tj0(k,ii,1,1));
                        zksav2 = 0.5*(zksav(j,k,ii,2)+tj0(k,ii,2,1));
                        cmuu = cmuv(j,k,ii);
                        bnum = zksav2*(1.0-f7)+f7;
                        // Durbin TCFD 1991 near-wall limiter
                        if (idurbinlim != 0 && (ieasm_type == 3 || ieasm_type == 4)) {
                            tau = bnum/zksav1;
                            fmu = t(k+jkv,14)*t(k+jkv,15)/xmre - t(k+jkv,24);
                            zkolmog = 6.0*std::sqrt(fmu*bnum/(rho*zksav2*zksav1));
                            { double a1=tau, a2=zkolmog; tau=ccomplex_ns::ccmax(a1,a2); }
                            zksav1 = bnum/tau;
                        }
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        eta1 = xis*(bnum/zksav1/re)*(bnum/zksav1/re);
                        alpa1 = -al4/(al10-eta1*al1*cmuu);
                        alpa2 = -2.0*al3/(al10-eta1*cmuu*al1);
                        t11 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*bnum/zksav1*(-s12*w12-s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t11 = t11*(1.0-bcj(k,ii,1));
                        t22 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s12*w12-s23*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t22 = t22*(1.0-bcj(k,ii,1));
                        t33 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s23*w23+s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t33 = t33*(1.0-bcj(k,ii,1));
                        t12 = -2.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t12 = t12*(1.0-bcj(k,ii,1));
                        t13 = -2.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t13 = t13*(1.0-bcj(k,ii,1));
                        t23 = -2.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t23 = t23*(1.0-bcj(k,ii,1));
                        t(k+jkv,28) = xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+t(k+jkv,27)*t13)/t(k+jkv,15);
                        t(k+jkv,29) = xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+t(k+jkv,27)*t23)/t(k+jkv,15);
                        t(k+jkv,30) = xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+t(k+jkv,27)*t33)/t(k+jkv,15);
                        xmut = f7*t(k+jkv,24)+(1.0-f7)*(double)cmuc1*rho*zksav2*zksav2/zksav1;
                        t(k+jkv,31) = xmre/t(k+jkv,15)*(
                            t(k+jkv,25)*(u*t11+v*t12+w_vel*t13)+
                            t(k+jkv,26)*(u*t12+v*t22+w_vel*t23)+
                            t(k+jkv,27)*(u*t13+v*t23+w_vel*t33))
                            - xmre/t(k+jkv,15)*(t(k+jkv,14)*t(k+jkv,15)/xmre-t(k+jkv,24)+(double)sigk1*xmut)*
                              (t(k+jkv,25)*t(k+jkv,25)+t(k+jkv,26)*t(k+jkv,26)+t(k+jkv,27)*t(k+jkv,27))*
                              (zksav(j,k,ii,2)-tj0(k,ii,2,1));
                    } // k loop j==1



                } else if (j == jdim) {
                    for (k = 1; k <= kdim1; k++) {
                        s11 = t(k+jkv,8)*t(k+jkv,25);
                        s22 = t(k+jkv,9)*t(k+jkv,26);
                        s33 = t(k+jkv,10)*t(k+jkv,27);
                        s12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)+t(k+jkv,9)*t(k+jkv,25));
                        s13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,25));
                        s23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,26));
                        w12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)-t(k+jkv,9)*t(k+jkv,25));
                        w13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,25));
                        w23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,26));
                        xis = s11*s11+s22*s22+s33*s33+2.0*s12*s12+2.0*s13*s13+2.0*s23*s23;
                        // Girimaji JFM 2000 fix to c4
                        if (ieasm_type == 4) {
                            wis = 2.0*w12*w12+2.0*w13*w13+2.0*w23*w23;
                            eta1_girimaji = xis/(xis+wis);
                            if ((float)eta1_girimaji < 0.5f) {
                                c4new = 2.0-((2.0-(double)c4_e)*
                                    std::pow(eta1_girimaji/(1.0-eta1_girimaji), 0.75));
                            } else {
                                c4new = (double)c4_e;
                            }
                            al4 = 0.5*c4new - 1.0;
                        }
                        rho = 0.5*(q(j-1,k,ii,1)+qj0(k,ii,1,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,1,3)*wj0(k+jkv2,15);
                        u   = 0.5*(q(j-1,k,ii,2)+qj0(k,ii,2,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,2,3)*wj0(k+jkv2,15);
                        v   = 0.5*(q(j-1,k,ii,3)+qj0(k,ii,3,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,3,3)*wj0(k+jkv2,15);
                        w_vel=0.5*(q(j-1,k,ii,4)+qj0(k,ii,4,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,4,3)*wj0(k+jkv2,15);
                        zksav1 = 0.5*(zksav(j-1,k,ii,1)+tj0(k,ii,1,3));
                        zksav2 = 0.5*(zksav(j-1,k,ii,2)+tj0(k,ii,2,3));
                        cmuu = cmuv(j-1,k,ii);
                        bnum = zksav2*(1.0-f7)+f7;
                        if (idurbinlim != 0 && (ieasm_type == 3 || ieasm_type == 4)) {
                            tau = bnum/zksav1;
                            fmu = t(k+jkv,14)*t(k+jkv,15)/xmre - t(k+jkv,24);
                            zkolmog = 6.0*std::sqrt(fmu*bnum/(rho*zksav2*zksav1));
                            { double a1=tau, a2=zkolmog; tau=ccomplex_ns::ccmax(a1,a2); }
                            zksav1 = bnum/tau;
                        }
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        eta1 = xis*(bnum/zksav1/re)*(bnum/zksav1/re);
                        alpa1 = -al4/(al10-eta1*al1*cmuu);
                        alpa2 = -2.0*al3/(al10-eta1*cmuu*al1);
                        t11 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*bnum/zksav1*(-s12*w12-s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t11 = t11*(1.0-bcj(k,ii,2));
                        t22 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s12*w12-s23*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t22 = t22*(1.0-bcj(k,ii,2));
                        t33 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s23*w23+s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t33 = t33*(1.0-bcj(k,ii,2));
                        t12 = -2.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t12 = t12*(1.0-bcj(k,ii,2));
                        t13 = -2.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t13 = t13*(1.0-bcj(k,ii,2));
                        t23 = -2.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t23 = t23*(1.0-bcj(k,ii,2));
                        t(k+jkv,28) = xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+t(k+jkv,27)*t13)/t(k+jkv,15);
                        t(k+jkv,29) = xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+t(k+jkv,27)*t23)/t(k+jkv,15);
                        t(k+jkv,30) = xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+t(k+jkv,27)*t33)/t(k+jkv,15);
                        xmut = f7*t(k+jkv,24)+(1.0-f7)*(double)cmuc1*rho*zksav2*zksav2/zksav1;
                        t(k+jkv,31) = xmre/t(k+jkv,15)*(
                            t(k+jkv,25)*(u*t11+v*t12+w_vel*t13)+
                            t(k+jkv,26)*(u*t12+v*t22+w_vel*t23)+
                            t(k+jkv,27)*(u*t13+v*t23+w_vel*t33))
                            - xmre/t(k+jkv,15)*(t(k+jkv,14)*t(k+jkv,15)/xmre-t(k+jkv,24)+(double)sigk1*xmut)*
                              (t(k+jkv,25)*t(k+jkv,25)+t(k+jkv,26)*t(k+jkv,26)+t(k+jkv,27)*t(k+jkv,27))*
                              (tj0(k,ii,2,3)-zksav(j-1,k,ii,2));
                    } // k loop j==jdim



                } else {
                    for (k = 1; k <= kdim1; k++) {
                        s11 = t(k+jkv,8)*t(k+jkv,25);
                        s22 = t(k+jkv,9)*t(k+jkv,26);
                        s33 = t(k+jkv,10)*t(k+jkv,27);
                        s12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)+t(k+jkv,9)*t(k+jkv,25));
                        s13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,25));
                        s23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,26));
                        w12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)-t(k+jkv,9)*t(k+jkv,25));
                        w13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,25));
                        w23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,26));
                        xis = s11*s11+s22*s22+s33*s33+2.0*s12*s12+2.0*s13*s13+2.0*s23*s23;
                        // Girimaji JFM 2000 fix to c4
                        if (ieasm_type == 4) {
                            wis = 2.0*w12*w12+2.0*w13*w13+2.0*w23*w23;
                            eta1_girimaji = xis/(xis+wis);
                            if ((float)eta1_girimaji < 0.5f) {
                                c4new = 2.0-((2.0-(double)c4_e)*
                                    std::pow(eta1_girimaji/(1.0-eta1_girimaji), 0.75));
                            } else {
                                c4new = (double)c4_e;
                            }
                            al4 = 0.5*c4new - 1.0;
                        }
                        rho = 0.5*(q(j,k,ii,1)+q(j-1,k,ii,1));
                        zksav1 = 0.5*(zksav(j,k,ii,1)+zksav(j-1,k,ii,1));
                        zksav2 = 0.5*(zksav(j,k,ii,2)+zksav(j-1,k,ii,2));
                        cmuu = 0.5*(cmuv(j,k,ii)+cmuv(j-1,k,ii));
                        bnum = zksav2*(1.0-f7)+f7;
                        if (idurbinlim != 0 && (ieasm_type == 3 || ieasm_type == 4)) {
                            tau = bnum/zksav1;
                            fmu = t(k+jkv,14)*t(k+jkv,15)/xmre - t(k+jkv,24);
                            zkolmog = 6.0*std::sqrt(fmu*bnum/(rho*zksav2*zksav1));
                            { double a1=tau, a2=zkolmog; tau=ccomplex_ns::ccmax(a1,a2); }
                            zksav1 = bnum/tau;
                        }
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        eta1 = xis*(bnum/zksav1/re)*(bnum/zksav1/re);
                        alpa1 = -al4/(al10-eta1*al1*cmuu);
                        alpa2 = -2.0*al3/(al10-eta1*cmuu*al1);
                        t11 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*bnum/zksav1*(-s12*w12-s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t22 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s12*w12-s23*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t33 = 2.0*rho*zksav2*re/3.0
                            - 4.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s23*w23+s13*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(k+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t12 = -2.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t13 = -2.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t23 = -2.0*alpa1*t(k+jkv,24)*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                            + 2.0*alpa2*t(k+jkv,24)*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t(k+jkv,28) = xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+t(k+jkv,27)*t13)/t(k+jkv,15);
                        t(k+jkv,29) = xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+t(k+jkv,27)*t23)/t(k+jkv,15);
                        t(k+jkv,30) = xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+t(k+jkv,27)*t33)/t(k+jkv,15);
                        xmut = f7*t(k+jkv,24)+(1.0-f7)*(double)cmuc1*rho*zksav2*zksav2/zksav1;
                        t(k+jkv,31) = xmre/t(k+jkv,15)*(
                            t(k+jkv,25)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t11
                                        +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t12
                                        +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t13)+
                            t(k+jkv,26)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t12
                                        +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t22
                                        +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t23)+
                            t(k+jkv,27)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t13
                                        +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t23
                                        +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t33))
                            - xmre/t(k+jkv,15)*(t(k+jkv,14)*t(k+jkv,15)/xmre-t(k+jkv,24)+(double)sigk1*xmut)*
                              (t(k+jkv,25)*t(k+jkv,25)+t(k+jkv,26)*t(k+jkv,26)+t(k+jkv,27)*t(k+jkv,27))*
                              (zksav(j,k,ii,2)-zksav(j-1,k,ii,2));
                    } // k loop interior
                }
            } // j loop
        } // ipl loop
    } // end EASM v2 branch



    // Generic nonlinear constitutive branch: i_nonlin!=0, i_tauijs!=1
    else if (i_tauijs != 1 && i_nonlin != 0 && level >= lglobal) {
        re = (double)reue / (double)xmach;
        // Add nonlinear terms to RHS
        for (ipl = 1; ipl <= npl; ipl++) {
            ii   = i + ipl - 1;
            jkv2 = (ipl - 1) * kdim1;
            for (j = 1; j <= jdim; j++) {
                jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
                if (j == 1) {
                    for (k = 1; k <= kdim1; k++) {
                        s11 = t(k+jkv,8)*t(k+jkv,25);
                        s22 = t(k+jkv,9)*t(k+jkv,26);
                        s33 = t(k+jkv,10)*t(k+jkv,27);
                        s12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)+t(k+jkv,9)*t(k+jkv,25));
                        s13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,25));
                        s23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,26));
                        w12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)-t(k+jkv,9)*t(k+jkv,25));
                        w13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,25));
                        w23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,26));
                        tmp = (s11+s22+s33)/3.0;
                        s11 = s11 - tmp;
                        s22 = s22 - tmp;
                        s33 = s33 - tmp;
                        // Note: denom term is taken at nearest cell center
                        denom = ux(j,k,ii,1)*ux(j,k,ii,1)+ux(j,k,ii,2)*ux(j,k,ii,2)+ux(j,k,ii,3)*ux(j,k,ii,3)
                              + ux(j,k,ii,4)*ux(j,k,ii,4)+ux(j,k,ii,5)*ux(j,k,ii,5)+ux(j,k,ii,6)*ux(j,k,ii,6)
                              + ux(j,k,ii,7)*ux(j,k,ii,7)+ux(j,k,ii,8)*ux(j,k,ii,8)+ux(j,k,ii,9)*ux(j,k,ii,9);
                        denom = std::sqrt(denom);
                        { double a1=denom, a2=(double)snonlin_lim; denom=ccomplex_ns::ccmax(a1,a2); }
                        rho = 0.5*(q(j,k,ii,1)+qj0(k,ii,1,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,1,1)*wj0(k+jkv2,5);
                        u   = 0.5*(q(j,k,ii,2)+qj0(k,ii,2,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,2,1)*wj0(k+jkv2,5);
                        v   = 0.5*(q(j,k,ii,3)+qj0(k,ii,3,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,3,1)*wj0(k+jkv2,5);
                        w_vel=0.5*(q(j,k,ii,4)+qj0(k,ii,4,1))*(1.0-wj0(k+jkv2,5))+qj0(k,ii,4,1)*wj0(k+jkv2,5);
                        zk_sa = 0.0;
                        if (ivisc[1]==6||ivisc[1]==7||ivisc[1]==10||ivisc[1]==15) {
                            zk_sa = 0.5*(zksav(j,k,ii,2)+tj0(k,ii,2,1));
                        }
                        t11 = 2.0*rho*zk_sa*re/3.0
                            - 8.0*(double)c_nonlin*t(k+jkv,24)/denom*(-s12*w12-s13*w13);
                        t11 = t11*(1.0-bcj(k,ii,1));
                        t22 = 2.0*rho*zk_sa*re/3.0
                            - 8.0*(double)c_nonlin*t(k+jkv,24)/denom*(s12*w12-s23*w23);
                        t22 = t22*(1.0-bcj(k,ii,1));
                        t33 = 2.0*rho*zk_sa*re/3.0
                            - 8.0*(double)c_nonlin*t(k+jkv,24)/denom*(s23*w23+s13*w13);
                        t33 = t33*(1.0-bcj(k,ii,1));
                        t12 = -4.0*(double)c_nonlin*t(k+jkv,24)/denom*(s11*w12-s22*w12-s13*w23-s23*w13);
                        t12 = t12*(1.0-bcj(k,ii,1));
                        t13 = -4.0*(double)c_nonlin*t(k+jkv,24)/denom*(s11*w13+s12*w23-s23*w12-s33*w13);
                        t13 = t13*(1.0-bcj(k,ii,1));
                        t23 = -4.0*(double)c_nonlin*t(k+jkv,24)/denom*(s12*w13+s22*w23+s13*w12-s33*w23);
                        t23 = t23*(1.0-bcj(k,ii,1));
                        t(k+jkv,28) = xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+t(k+jkv,27)*t13)/t(k+jkv,15);
                        t(k+jkv,29) = xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+t(k+jkv,27)*t23)/t(k+jkv,15);
                        t(k+jkv,30) = xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+t(k+jkv,27)*t33)/t(k+jkv,15);
                        t(k+jkv,31) = xmre/t(k+jkv,15)*(
                            t(k+jkv,25)*(u*t11+v*t12+w_vel*t13)+
                            t(k+jkv,26)*(u*t12+v*t22+w_vel*t23)+
                            t(k+jkv,27)*(u*t13+v*t23+w_vel*t33));
                    } // k loop j==1
                } else if (j == jdim) {
                    for (k = 1; k <= kdim1; k++) {
                        s11 = t(k+jkv,8)*t(k+jkv,25);
                        s22 = t(k+jkv,9)*t(k+jkv,26);
                        s33 = t(k+jkv,10)*t(k+jkv,27);
                        s12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)+t(k+jkv,9)*t(k+jkv,25));
                        s13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,25));
                        s23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,26));
                        w12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)-t(k+jkv,9)*t(k+jkv,25));
                        w13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,25));
                        w23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,26));
                        tmp = (s11+s22+s33)/3.0;
                        s11 = s11 - tmp;
                        s22 = s22 - tmp;
                        s33 = s33 - tmp;
                        // Note: denom term is taken at nearest cell center
                        denom = ux(j-1,k,ii,1)*ux(j-1,k,ii,1)+ux(j-1,k,ii,2)*ux(j-1,k,ii,2)
                              + ux(j-1,k,ii,3)*ux(j-1,k,ii,3)+ux(j-1,k,ii,4)*ux(j-1,k,ii,4)
                              + ux(j-1,k,ii,5)*ux(j-1,k,ii,5)+ux(j-1,k,ii,6)*ux(j-1,k,ii,6)
                              + ux(j-1,k,ii,7)*ux(j-1,k,ii,7)+ux(j-1,k,ii,8)*ux(j-1,k,ii,8)
                              + ux(j-1,k,ii,9)*ux(j-1,k,ii,9);
                        denom = std::sqrt(denom);
                        { double a1=denom, a2=(double)snonlin_lim; denom=ccomplex_ns::ccmax(a1,a2); }
                        rho = 0.5*(q(j-1,k,ii,1)+qj0(k,ii,1,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,1,3)*wj0(k+jkv2,15);
                        u   = 0.5*(q(j-1,k,ii,2)+qj0(k,ii,2,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,2,3)*wj0(k+jkv2,15);
                        v   = 0.5*(q(j-1,k,ii,3)+qj0(k,ii,3,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,3,3)*wj0(k+jkv2,15);
                        w_vel=0.5*(q(j-1,k,ii,4)+qj0(k,ii,4,3))*(1.0-wj0(k+jkv2,15))+qj0(k,ii,4,3)*wj0(k+jkv2,15);
                        zk_sa = 0.0;
                        if (ivisc[1]==6||ivisc[1]==7||ivisc[1]==10||ivisc[1]==15) {
                            zk_sa = 0.5*(zksav(j-1,k,ii,2)+tj0(k,ii,2,3));
                        }
                        t11 = 2.0*rho*zk_sa*re/3.0
                            - 8.0*(double)c_nonlin*t(k+jkv,24)/denom*(-s12*w12-s13*w13);
                        t11 = t11*(1.0-bcj(k,ii,2));
                        t22 = 2.0*rho*zk_sa*re/3.0
                            - 8.0*(double)c_nonlin*t(k+jkv,24)/denom*(s12*w12-s23*w23);
                        t22 = t22*(1.0-bcj(k,ii,2));
                        t33 = 2.0*rho*zk_sa*re/3.0
                            - 8.0*(double)c_nonlin*t(k+jkv,24)/denom*(s23*w23+s13*w13);
                        t33 = t33*(1.0-bcj(k,ii,2));
                        t12 = -4.0*(double)c_nonlin*t(k+jkv,24)/denom*(s11*w12-s22*w12-s13*w23-s23*w13);
                        t12 = t12*(1.0-bcj(k,ii,2));
                        t13 = -4.0*(double)c_nonlin*t(k+jkv,24)/denom*(s11*w13+s12*w23-s23*w12-s33*w13);
                        t13 = t13*(1.0-bcj(k,ii,2));
                        t23 = -4.0*(double)c_nonlin*t(k+jkv,24)/denom*(s12*w13+s22*w23+s13*w12-s33*w23);
                        t23 = t23*(1.0-bcj(k,ii,2));
                        t(k+jkv,28) = xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+t(k+jkv,27)*t13)/t(k+jkv,15);
                        t(k+jkv,29) = xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+t(k+jkv,27)*t23)/t(k+jkv,15);
                        t(k+jkv,30) = xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+t(k+jkv,27)*t33)/t(k+jkv,15);
                        t(k+jkv,31) = xmre/t(k+jkv,15)*(
                            t(k+jkv,25)*(u*t11+v*t12+w_vel*t13)+
                            t(k+jkv,26)*(u*t12+v*t22+w_vel*t23)+
                            t(k+jkv,27)*(u*t13+v*t23+w_vel*t33));
                    } // k loop j==jdim



                } else {
                    for (k = 1; k <= kdim1; k++) {
                        s11 = t(k+jkv,8)*t(k+jkv,25);
                        s22 = t(k+jkv,9)*t(k+jkv,26);
                        s33 = t(k+jkv,10)*t(k+jkv,27);
                        s12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)+t(k+jkv,9)*t(k+jkv,25));
                        s13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,25));
                        s23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)+t(k+jkv,10)*t(k+jkv,26));
                        w12 = 0.5*(t(k+jkv,8)*t(k+jkv,26)-t(k+jkv,9)*t(k+jkv,25));
                        w13 = 0.5*(t(k+jkv,8)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,25));
                        w23 = 0.5*(t(k+jkv,9)*t(k+jkv,27)-t(k+jkv,10)*t(k+jkv,26));
                        tmp = (s11+s22+s33)/3.0;
                        s11 = s11 - tmp;
                        s22 = s22 - tmp;
                        s33 = s33 - tmp;
                        // Note: denom term is averaged to cell face
                        ux1 = 0.5*(ux(j-1,k,ii,1)+ux(j,k,ii,1));
                        ux2 = 0.5*(ux(j-1,k,ii,2)+ux(j,k,ii,2));
                        ux3 = 0.5*(ux(j-1,k,ii,3)+ux(j,k,ii,3));
                        ux4 = 0.5*(ux(j-1,k,ii,4)+ux(j,k,ii,4));
                        ux5 = 0.5*(ux(j-1,k,ii,5)+ux(j,k,ii,5));
                        ux6 = 0.5*(ux(j-1,k,ii,6)+ux(j,k,ii,6));
                        ux7 = 0.5*(ux(j-1,k,ii,7)+ux(j,k,ii,7));
                        ux8 = 0.5*(ux(j-1,k,ii,8)+ux(j,k,ii,8));
                        ux9 = 0.5*(ux(j-1,k,ii,9)+ux(j,k,ii,9));
                        denom = ux1*ux1+ux2*ux2+ux3*ux3+ux4*ux4+ux5*ux5+ux6*ux6
                              + ux7*ux7+ux8*ux8+ux9*ux9;
                        denom = std::sqrt(denom);
                        { double a1=denom, a2=(double)snonlin_lim; denom=ccomplex_ns::ccmax(a1,a2); }
                        rho = 0.5*(q(j,k,ii,1)+q(j-1,k,ii,1));
                        zk_sa = 0.0;
                        if (ivisc[1]==6||ivisc[1]==7||ivisc[1]==10||ivisc[1]==15) {
                            zk_sa = 0.5*(zksav(j,k,ii,2)+zksav(j-1,k,ii,2));
                        }
                        t11 = 2.0*rho*zk_sa*re/3.0
                            - 8.0*(double)c_nonlin*t(k+jkv,24)/denom*(-s12*w12-s13*w13);
                        t22 = 2.0*rho*zk_sa*re/3.0
                            - 8.0*(double)c_nonlin*t(k+jkv,24)/denom*(s12*w12-s23*w23);
                        t33 = 2.0*rho*zk_sa*re/3.0
                            - 8.0*(double)c_nonlin*t(k+jkv,24)/denom*(s23*w23+s13*w13);
                        t12 = -4.0*(double)c_nonlin*t(k+jkv,24)/denom*(s11*w12-s22*w12-s13*w23-s23*w13);
                        t13 = -4.0*(double)c_nonlin*t(k+jkv,24)/denom*(s11*w13+s12*w23-s23*w12-s33*w13);
                        t23 = -4.0*(double)c_nonlin*t(k+jkv,24)/denom*(s12*w13+s22*w23+s13*w12-s33*w23);
                        t(k+jkv,28) = xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+t(k+jkv,27)*t13)/t(k+jkv,15);
                        t(k+jkv,29) = xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+t(k+jkv,27)*t23)/t(k+jkv,15);
                        t(k+jkv,30) = xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+t(k+jkv,27)*t33)/t(k+jkv,15);
                        t(k+jkv,31) = xmre/t(k+jkv,15)*(
                            t(k+jkv,25)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t11
                                        +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t12
                                        +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t13)+
                            t(k+jkv,26)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t12
                                        +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t22
                                        +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t23)+
                            t(k+jkv,27)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t13
                                        +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t23
                                        +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t33));
                    } // k loop interior
                }
            } // j loop
        } // ipl loop
    } // end generic nonlinear branch

    // Final residual update for nonlinear branches (ivisc==11,12,13,14 or i_nonlin!=0)
    if (i_tauijs != 1 &&
        (ivisc[1]==11||ivisc[1]==12||ivisc[1]==13||ivisc[1]==14||i_nonlin!=0) &&
        level >= lglobal) {
        for (l = 2; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                // Take positive derivative for correct sign
                t(izz, l) = t(izz + ks - 1, l + 26) - t(izz, l + 26);
            }
        }
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            for (j = 1; j <= jdim1; j++) {
                jkv = (j - 1) * kdim1 * npl + (ipl - 1) * kdim1 + 1;
                for (l = 2; l <= 5; l++) {
                    for (k = 1; k <= kdim1; k++) {
                        res(j, k, ii, l) = res(j, k, ii, l) + t(k + jkv - 1, l);
                    }
                }
            }
        }
    }

    return;
}

} // namespace gfluxv_ns
