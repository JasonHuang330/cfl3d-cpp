// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "hfluxv.h"
#include "ccomplex.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

namespace hfluxv_ns {

void hfluxv(int& i, int& npl, int& jdim, int& kdim, int& idim, int& idf,
            FortranArray4DRef<double> ak, FortranArray4DRef<double> bk,
            FortranArray4DRef<double> ck, FortranArray4DRef<double> res,
            FortranArray4DRef<double> q, FortranArray4DRef<double> qk0,
            FortranArray3DRef<double> sk, FortranArray3DRef<double> vol,
            FortranArray2DRef<double> t, int& nvtq,
            FortranArray2DRef<double> wk0, FortranArray3DRef<double> vist3d,
            FortranArray3DRef<double> vmuk, FortranArray4DRef<double> vk0,
            FortranArray3DRef<double> bck, FortranArray4DRef<double> zksav,
            FortranArray4DRef<double> tk0, FortranArray3DRef<double> cmuv,
            FortranArray3DRef<double> volk0,
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
    int32_t* ivisc  = cmn_reyue.ivisc;  // 0-based: ivisc[2] = ivisc(3)
    int32_t& isklton = cmn_sklton.isklton;
    int32_t& level  = cmn_mgrd.level;
    int32_t& lglobal = cmn_mgrd.lglobal;
    // /easmv/ - rename to avoid conflicts
    float& c10_e    = cmn_easmv.c10;
    float& c11_e    = cmn_easmv.c11;
    float& c2_e     = cmn_easmv.c2;
    float& c3_e     = cmn_easmv.c3;
    float& c4_e     = cmn_easmv.c4;
    float& c5_e     = cmn_easmv.c5;
    float& sigk1_e  = cmn_easmv.sigk1;
    float& cmuc1_e  = cmn_easmv.cmuc1;
    int32_t& ieasm_type_e = cmn_easmv.ieasm_type;
    // /constit/
    int32_t& i_nonlin   = cmn_constit.i_nonlin;
    float& c_nonlin     = cmn_constit.c_nonlin;
    float& snonlin_lim  = cmn_constit.snonlin_lim;
    int32_t& i_tauijs   = cmn_constit.i_tauijs;
    int32_t& i_qcr2000  = cmn_constit.i_qcr2000;
    int32_t& i_qcr2013  = cmn_constit.i_qcr2013;
    int32_t& i_qcr2013v = cmn_constit.i_qcr2013v;
    // /ivals/
    float* tur10 = cmn_ivals.tur10;  // 0-based: tur10[0] = tur10(1)

    // Local variables
    float coef_eddy;
    double ccr2, kdim1_d, jdim1_d;
    int kdim1, jdim1, n, l0, jv, js, nn;
    double xmre, gpr, gm1pr, prtr;
    int idurbinlim;
    double c2b, c2bp;
    int n1, l1;
    int iviscc;
    double tr1;
    int ns, n0;
    double gterm1, gterm2;

    // Scalar temporaries
    double t5, t6, ab, bb, wk05, wk06;
    double t1_tmp, t24_tmp, t25_tmp;
    double z1, z2;
    double rho, u, v, w_vel;
    double v3d;
    double ux1, ux2, ux3, ux4, ux5_v, ux6, ux7, ux8, ux9;
    double s11, s22, s33, s12, s13, s23;
    double w12, w13, w23;
    double t11, t22, t33, t12, t13, t23;
    double s11star, s22star, s33star;
    double denom, o11, o22, o33, o12, o13, o23;
    double t11p, t22p, t33p, t12p, t13p, t23p;
    double xis2013, s_mod, wis, w_mod;
    double xis;
    double re, gg, f7;
    double alpa1, alpa2;
    double zksav1, zksav2, bnum;
    double eta, squig, factre;
    double xlnpt;
    double al10, al1, al2, al3, al4;
    double eta1, cmuu;
    double tau, fmu, zkolmog;
    double c4new, eta1_girimaji;
    double xmut;
    double zk_sa;
    double tmp;
    int ipl, ii, k, j, izz, jkv, jkv2, jk, jc, m, mm, mb;
    int imin, imax, jmin, jmax;
    int kp1, km2;
    int l, k1, k2;


    // coefficient used to turn off eddy-viscosity based
    // turbulent stress calculation for ivmx==70 (full Reynolds stress model)
    coef_eddy = 1.0f;
    if (ivisc[2] >= 70 && idf == 0) coef_eddy = 0.0f;
    if (i_tauijs == 1 && idf == 0) coef_eddy = 0.0f;
    ccr2 = 2.5;

    kdim1 = kdim - 1;
    jdim1 = jdim - 1;

    // n  : number of cell centers for npl planes
    // l0 : number of cell interfaces for npl planes
    // jv : number of cell centers (and interfaces) on a k=constant plane
    n  = npl * jdim1 * kdim1;
    l0 = npl * jdim1 * kdim;
    jv = npl * jdim1;
    js = jv + 1;
    nn = n - jv;

    xmre  = (double)xmach / (double)reue;
    gpr   = (double)gamma / (double)pr;
    gm1pr = (double)gm1 * (double)pr;
    prtr  = (double)pr / (double)prt;
    // Durbin TCFD 1991 near-wall limiter (0=off)
    idurbinlim = 0;

    if (isklton > 0 && i == 1 && iadv >= 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120,
            "   computing viscous  fluxes, K-direction - central differencing");
        if (i_tauijs == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "      using tau_ij method");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "(note: not done accurately at boundaries for multiblock cases...)");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "currently one sided; need ux derivs from across boundaries)");
        }
        if (i_qcr2000 == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "      using QCR2000");
        }
        if (i_qcr2013 == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "      using QCR2013 (limited)");
        }
        if (i_qcr2013v == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "      using QCR2013-V");
        }
        if (ivisc[2] == 11 || ivisc[2] == 12 || ivisc[2] == 13 ||
            ivisc[2] == 14 || i_nonlin != 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "      nonlinear constitutive relation");
        }
        if (i_nonlin != 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "      snonlin_lim=%19.8E", (double)snonlin_lim);
        }
    }

    // store selected cell centered data in t array
    // t(1)  : 1/density
    // t(16) : c*c/(gm1*pr)
    // t(18) : u
    // t(19) : v
    // t(20) : w
    // t(22) : density
    l1 = jdim * kdim - 1;
    for (ipl = 1; ipl <= npl; ipl++) {
        ii = i + ipl - 1;
        for (k = 1; k <= kdim1; k++) {
            int jkv_loc = (k-1)*npl*jdim1 + (ipl-1)*jdim1 + 1;
            for (izz = 1; izz <= jdim1; izz++) {
                t(izz+jkv_loc-1, 22) = q(izz, k, ii, 1);
                t(izz+jkv_loc-1, 18) = q(izz, k, ii, 2);
                t(izz+jkv_loc-1, 19) = q(izz, k, ii, 3);
                t(izz+jkv_loc-1, 20) = q(izz, k, ii, 4);
                t(izz+jkv_loc-1, 16) = q(izz, k, ii, 5);
            }
        }
    }

    // store dependent variables on k=0 and k=kdim boundaries in wk0
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
        jmin = 1;
        jmax = jdim - 1;
        for (ii = imin; ii <= imax; ii++) {
            jc = (ii - i) * jdim1 + jmin - 1;
            for (int jj = jmin; jj <= jmax; jj++) {
                jc = jc + 1;
                wk0(jc, m+12) = qk0(jj, ii, 1, mm);
                wk0(jc, m+8)  = qk0(jj, ii, 2, mm);
                wk0(jc, m+9)  = qk0(jj, ii, 3, mm);
                wk0(jc, m+10) = qk0(jj, ii, 4, mm);
                wk0(jc, m+6)  = qk0(jj, ii, 5, mm);
                // m+5 flag indicates whether ghost cell or interface values stored in wk0
                wk0(jc, m+5)  = bck(jj, ii, mb);
            }
        }
    }

    // t(17)/wk0(m+7)   : (u*u+v*v+w*w)/2
    // t(1)/wk0(m+1)    : 1/density
    // t(16)/wk0(m+6)   : c*c/(gm1*pr)
    for (izz = 1; izz <= n; izz++) {
        t(izz, 17) = 0.5e0*(t(izz,18)*t(izz,18)
                           +t(izz,19)*t(izz,19)
                           +t(izz,20)*t(izz,20));
        t(izz, 1)  = 1.e0 / t(izz, 22);
        t(izz, 16) = gpr * t(izz, 16) * t(izz, 1) / (double)gm1;
    }
    for (m = 0; m <= 10; m += 10) {
        for (izz = 1; izz <= jv; izz++) {
            wk0(izz, m+7) = 0.5e0*(wk0(izz,m+8)*wk0(izz,m+8)
                                  +wk0(izz,m+9)*wk0(izz,m+9)
                                  +wk0(izz,m+10)*wk0(izz,m+10));
            wk0(izz, m+1) = 1.e0 / wk0(izz, m+12);
            wk0(izz, m+6) = gpr * wk0(izz, m+6) * wk0(izz, m+1) / (double)gm1;
        }
    }

    // t(7) : laminar viscosity at cell centers (via sutherland relation)
    c2b  = (double)cbar / (double)tinf;
    c2bp = c2b + 1.e0;

    for (izz = 1; izz <= n; izz++) {
        t5       = gm1pr * t(izz, 16);
        t6       = std::sqrt(t5);
        t(izz,7) = c2bp * t5 * t6 / (c2b + t5);
    }

    // t(14) : laminar viscosity values at cell interfaces
    // interior interfaces
    for (izz = 1; izz <= nn; izz++) {
        t(izz+js-1, 14) = (t(izz,7) + t(izz+js-1,7)) * 0.5e0;
    }
    // k=0 and k=kdim interfaces
    for (izz = 1; izz <= jv; izz++) {
        ab   = 1.0 + wk0(izz, 5);
        bb   = 1.0 - wk0(izz, 5);
        wk05 = gm1pr * 0.5 * (ab*wk0(izz,6) + bb*t(izz,16));
        wk06 = std::sqrt(wk05);
        t(izz, 14)   = c2bp * wk05 * wk06 / (c2b + wk05);
        ab   = 1.0 + wk0(izz, 15);
        bb   = 1.0 - wk0(izz, 15);
        wk05 = gm1pr * 0.5 * (ab*wk0(izz,16) + bb*t(izz+n-jv,16));
        wk06 = std::sqrt(wk05);
        t(izz+n, 14) = c2bp * wk05 * wk06 / (c2b + wk05);
    }


    // t(15) : average jacobian (inverse volume) at cell interface
    n1 = jdim * kdim1 + 1;
    l1 = jdim * kdim - 1;
    for (ipl = 1; ipl <= npl; ipl++) {
        ii  = i + ipl - 1;
        for (k = 1; k <= kdim; k++) {
            jkv = (k-1)*npl*jdim1 + (ipl-1)*jdim1;
            jk  = (k-1)*jdim;
            if (k == 1) {
                // inverse volume at k=0 interface
                for (j = 1; j <= jdim1; j++) {
                    t(j+jkv, 15) = 2.0 / (volk0(j,ii,1) + vol(j,1,ii));
                }
            } else if (k == kdim) {
                // inverse volume at k=kdim interface
                for (j = 1; j <= jdim1; j++) {
                    t(j+jkv, 15) = 2.0 / (volk0(j,ii,3) + vol(j,kdim1,ii));
                }
            } else {
                // inverse volume at interior interfaces
                for (j = 1; j <= jdim1; j++) {
                    t(j+jkv, 15) = 2.0 / (vol(j,k,ii) + vol(j,k-1,ii));
                }
            }

            // t(25) - t(27) : components of grad(zeta)
            // t1    : grad(zeta)/J
            // t(25) : d(zeta)/dx
            // t(26) : d(zeta)/dy
            // t(27) : d(zeta)/dz
            for (j = 1; j <= jdim1; j++) {
                double t1_loc = sk(j+jk, ii, 4) * t(j+jkv, 15);
                t(j+jkv, 25) = sk(j+jk, ii, 1) * t1_loc;
                t(j+jkv, 26) = sk(j+jk, ii, 2) * t1_loc;
                t(j+jkv, 27) = sk(j+jk, ii, 3) * t1_loc;
            }
        }
    }

    // t(6)-t(10) : gradients at cell interfaces
    // t(6)       : d( c*c/(gm1*pr) )/d(zeta)
    // t(7)       : d( (u*u+v*v+w*w)/2 )/d(zeta)
    // t(8)       : d(u)/d(zeta)
    // t(9)       : d(v)/d(zeta)
    // t(10)      : d(w)/d(zeta)
    for (l = 1; l <= 5; l++) {
        k1 = l + 5;
        k2 = l + 15;
        // interior interfaces
        for (izz = 1; izz <= nn; izz++) {
            t(izz+js-1, k1) = t(izz+js-1, k2) - t(izz, k2);
        }
        // interfaces at k=0/kdim
        for (izz = 1; izz <= jv; izz++) {
            z1 = 1.0 + wk0(izz, 5);
            z2 = 1.0 + wk0(izz, 15);
            t(izz,   k1) = z1 * (t(izz, k2) - wk0(izz, k1));
            t(izz+n, k1) = z2 * (wk0(izz, k2) - t(izz+n-jv, k2));
        }
    }

    // t(24) : turbulent viscosity at interfaces (=0 for laminar flow)
    iviscc = ivisc[2];
    if (iviscc > 1) {
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            for (k = 1; k <= kdim; k++) {
                jkv = (k-1)*npl*jdim1 + (ipl-1)*jdim1;
                // k=0 interface
                if (k == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        t(j+jkv, 24) = bck(j,ii,1)*vk0(j,ii,1,1) +
                            (1.0 - bck(j,ii,1))*0.5*(vk0(j,ii,1,1) + vist3d(j,1,ii));
                    }
                // k=kdim interface
                } else if (k == kdim) {
                    for (j = 1; j <= jdim1; j++) {
                        t(j+jkv, 24) = bck(j,ii,2)*vk0(j,ii,1,3) +
                            (1.0 - bck(j,ii,2))*0.5*(vk0(j,ii,1,3) + vist3d(j,kdim1,ii));
                    }
                } else {
                    // interior interfaces
                    for (j = 1; j <= jdim1; j++) {
                        t(j+jkv, 24) = 0.5*(vist3d(j,k,ii) + vist3d(j,k-1,ii));
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
        t25_tmp = t(izz, 24) / t(izz, 14);
        t24_tmp = 1.e0 + t25_tmp * (double)coef_eddy;
        // t(23) : ratio of laminar Prandtl number to total Prandtl number
        t(izz, 23) = (1.e0 + prtr*t25_tmp) / t24_tmp;
        // t(14): (mach/re)*viscosity/J
        t(izz, 14) = xmre * t24_tmp * t(izz, 14) / t(izz, 15);
        // t(5) : t(14)*(grad(zeta))**2
        t(izz, 5) = t(izz, 14) * (t(izz,25)*t(izz,25)
                                  +t(izz,26)*t(izz,26)
                                  +t(izz,27)*t(izz,27));
    }

    // t(11) : t(14)*((d(zeta)/dx * d(u)/d(zeta))
    //               + (d(zeta)/dy * d(v)/d(zeta))
    //               + (d(zeta)/dz * d(w)/d(zeta)))/3
    tr1 = 1.e0 / 3.e0;
    for (izz = 1; izz <= l0; izz++) {
        t(izz, 11) = tr1 * t(izz, 14) * (t(izz,25)*t(izz,8)
                                         +t(izz,26)*t(izz,9)
                                         +t(izz,27)*t(izz,10));
    }

    // store off vmuk....viscosity coefficient (laminar + turbulent) on k=0 and
    // k=kdim boundaries for later use in force integration
    for (ipl = 1; ipl <= npl; ipl++) {
        ii = i + ipl - 1;
        for (j = 1; j <= jdim1; j++) {
            izz = (ipl-1)*jdim1 + j;
            vmuk(j, ii, 1) = t(izz,   14) * t(izz,   15) / xmre;
            vmuk(j, ii, 2) = t(izz+n, 14) * t(izz+n, 15) / xmre;
        }
    }

    if (iadv < 0) return;

    // t(12) : contravariant velocity W at cell interfaces
    // interior interfaces
    for (izz = 1; izz <= nn; izz++) {
        t(izz+js-1, 12) = 0.5e0*(t(izz+js-1,25)*(t(izz+js-1,18)+t(izz,18))+
                                  t(izz+js-1,26)*(t(izz+js-1,19)+t(izz,19))+
                                  t(izz+js-1,27)*(t(izz+js-1,20)+t(izz,20)));
    }
    // interfaces at k=0 and k=kdim
    for (izz = 1; izz <= jv; izz++) {
        ab = 1.0 + wk0(izz, 15);
        bb = 1.0 - wk0(izz, 15);
        t(izz+n, 12) = 0.5*(t(izz+n,25)*(ab*wk0(izz,18)+bb*t(izz+n-jv,18))
                           +t(izz+n,26)*(ab*wk0(izz,19)+bb*t(izz+n-jv,19))
                           +t(izz+n,27)*(ab*wk0(izz,20)+bb*t(izz+n-jv,20)));
        ab = 1.0 + wk0(izz, 5);
        bb = 1.0 - wk0(izz, 5);
        t(izz, 12) = 0.5*(t(izz,25)*(ab*wk0(izz,8) +bb*t(izz,18))
                         +t(izz,26)*(ab*wk0(izz,9) +bb*t(izz,19))
                         +t(izz,27)*(ab*wk0(izz,10)+bb*t(izz,20)));
    }

    if (idf > 0) goto label_3333;


    // calculate fluxes
    for (l = 2; l <= 4; l++) {
        // viscous terms at interfaces (momentum eqs.)
        for (izz = 1; izz <= l0; izz++) {
            t(izz, 16) = t(izz,5)*t(izz,l+6) + t(izz,l+23)*t(izz,11);
        }
        // (-)viscous flux = Hv(k-1/2)-Hv(k+1/2) (momentum eqs.)
        for (izz = 1; izz <= n; izz++) {
            t(izz, l) = -t(izz+js-1, 16) + t(izz, 16);
        }
    }

    // viscous terms at interfaces (energy eq.)
    for (izz = 1; izz <= l0; izz++) {
        t(izz, 16) = t(izz,5)*(t(izz,23)*t(izz,6)+t(izz,7)) +
                     t(izz,12)*t(izz,11);
    }
    // (-)viscous flux = Hv(k-1/2)-Hv(k+1/2) (energy eq.)
    for (izz = 1; izz <= n; izz++) {
        t(izz, 5) = -t(izz+js-1, 16) + t(izz, 16);
    }

    // calculate residuals
    for (ipl = 1; ipl <= npl; ipl++) {
        ii  = i + ipl - 1;
        for (k = 1; k <= kdim1; k++) {
            jkv = (k-1)*jdim1*npl + (ipl-1)*jdim1 + 1;
            for (l = 2; l <= 5; l++) {
                for (izz = 1; izz <= jdim1; izz++) {
                    res(izz, k, ii, l) = res(izz, k, ii, l) + t(izz+jkv-1, l);
                }
            }
        }
    }

    // Add tau_ij's for general method:
    if (i_tauijs == 1 && ivisc[2] < 70) {
        re = (double)reue / (double)xmach;
        for (ipl = 1; ipl <= npl; ipl++) {
            ii   = i + ipl - 1;
            jkv2 = (ipl-1)*jdim1;
            for (k = 1; k <= kdim; k++) {
                jk  = (k-1)*jdim;
                jkv = (k-1)*npl*jdim1 + (ipl-1)*jdim1;
                kp1 = std::min(k+1, kdim-1);
                km2 = std::max(k-2, 1);
                for (j = 1; j <= jdim1; j++) {
                    if (k == 1) {
                        rho   = 0.5*(q(j,k,ii,1)+qk0(j,ii,1,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,1,1)*wk0(j+jkv2,5);
                        u     = 0.5*(q(j,k,ii,2)+qk0(j,ii,2,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,2,1)*wk0(j+jkv2,5);
                        v     = 0.5*(q(j,k,ii,3)+qk0(j,ii,3,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,3,1)*wk0(j+jkv2,5);
                        w_vel = 0.5*(q(j,k,ii,4)+qk0(j,ii,4,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,4,1)*wk0(j+jkv2,5);
                        v3d   = 0.5*(vist3d(j,k,ii)+vk0(j,ii,1,1))*(1.-wk0(j+jkv2,5))+
                                vk0(j,ii,1,1)*wk0(j+jkv2,5);
                        ux1   = 1.5*ux(j,k,ii,1)   - 0.5*ux(j,kp1,ii,1);
                        ux2   = 1.5*ux(j,k,ii,2)   - 0.5*ux(j,kp1,ii,2);
                        ux3   = 1.5*ux(j,k,ii,3)   - 0.5*ux(j,kp1,ii,3);
                        ux4   = 1.5*ux(j,k,ii,4)   - 0.5*ux(j,kp1,ii,4);
                        ux5_v = 1.5*ux(j,k,ii,5)   - 0.5*ux(j,kp1,ii,5);
                        ux6   = 1.5*ux(j,k,ii,6)   - 0.5*ux(j,kp1,ii,6);
                        ux7   = 1.5*ux(j,k,ii,7)   - 0.5*ux(j,kp1,ii,7);
                        ux8   = 1.5*ux(j,k,ii,8)   - 0.5*ux(j,kp1,ii,8);
                        ux9   = 1.5*ux(j,k,ii,9)   - 0.5*ux(j,kp1,ii,9);
                    } else if (k == kdim) {
                        rho   = 0.5*(q(j,k-1,ii,1)+qk0(j,ii,1,3))*
                                (1.-wk0(j+jkv2,15))+
                                qk0(j,ii,1,3)*wk0(j+jkv2,15);
                        u     = 0.5*(q(j,k-1,ii,2)+qk0(j,ii,2,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,2,3)*wk0(j+jkv2,15);
                        v     = 0.5*(q(j,k-1,ii,3)+qk0(j,ii,3,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,3,3)*wk0(j+jkv2,15);
                        w_vel = 0.5*(q(j,k-1,ii,4)+qk0(j,ii,4,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,4,3)*wk0(j+jkv2,15);
                        v3d   = 0.5*(vist3d(j,k-1,ii)+vk0(j,ii,1,3))*(1.-wk0(j+jkv2,15))+
                                vk0(j,ii,1,3)*wk0(j+jkv2,15);
                        ux1   = 1.5*ux(j,k-1,ii,1) - 0.5*ux(j,km2,ii,1);
                        ux2   = 1.5*ux(j,k-1,ii,2) - 0.5*ux(j,km2,ii,2);
                        ux3   = 1.5*ux(j,k-1,ii,3) - 0.5*ux(j,km2,ii,3);
                        ux4   = 1.5*ux(j,k-1,ii,4) - 0.5*ux(j,km2,ii,4);
                        ux5_v = 1.5*ux(j,k-1,ii,5) - 0.5*ux(j,km2,ii,5);
                        ux6   = 1.5*ux(j,k-1,ii,6) - 0.5*ux(j,km2,ii,6);
                        ux7   = 1.5*ux(j,k-1,ii,7) - 0.5*ux(j,km2,ii,7);
                        ux8   = 1.5*ux(j,k-1,ii,8) - 0.5*ux(j,km2,ii,8);
                        ux9   = 1.5*ux(j,k-1,ii,9) - 0.5*ux(j,km2,ii,9);
                    } else {
                        rho   = 0.5*(q(j,k,ii,1)+q(j,k-1,ii,1));
                        u     = 0.5*(q(j,k,ii,2)+q(j,k-1,ii,2));
                        v     = 0.5*(q(j,k,ii,3)+q(j,k-1,ii,3));
                        w_vel = 0.5*(q(j,k,ii,4)+q(j,k-1,ii,4));
                        v3d   = 0.5*(vist3d(j,k,ii)+vist3d(j,k-1,ii));
                        ux1   = 0.5*(ux(j,k,ii,1)+ux(j,k-1,ii,1));
                        ux2   = 0.5*(ux(j,k,ii,2)+ux(j,k-1,ii,2));
                        ux3   = 0.5*(ux(j,k,ii,3)+ux(j,k-1,ii,3));
                        ux4   = 0.5*(ux(j,k,ii,4)+ux(j,k-1,ii,4));
                        ux5_v = 0.5*(ux(j,k,ii,5)+ux(j,k-1,ii,5));
                        ux6   = 0.5*(ux(j,k,ii,6)+ux(j,k-1,ii,6));
                        ux7   = 0.5*(ux(j,k,ii,7)+ux(j,k-1,ii,7));
                        ux8   = 0.5*(ux(j,k,ii,8)+ux(j,k-1,ii,8));
                        ux9   = 0.5*(ux(j,k,ii,9)+ux(j,k-1,ii,9));
                    }
                    s11 = ux1;
                    s22 = ux5_v;
                    s33 = ux9;
                    s11star = s11 - ((s11+s22+s33)/3.0);
                    s22star = s22 - ((s11+s22+s33)/3.0);
                    s33star = s33 - ((s11+s22+s33)/3.0);
                    s12 = 0.5*(ux2+ux4);
                    s13 = 0.5*(ux3+ux7);
                    s23 = 0.5*(ux6+ux8);
                    w12 = 0.5*(ux2-ux4);
                    w13 = 0.5*(ux3-ux7);
                    w23 = 0.5*(ux6-ux8);
                    // ignoring -2/3*rho*k*delij term
                    t11 = -2.0*v3d*s11star;
                    t22 = -2.0*v3d*s22star;
                    t33 = -2.0*v3d*s33star;
                    t12 = -2.0*v3d*s12;
                    t13 = -2.0*v3d*s13;
                    t23 = -2.0*v3d*s23;

                    if (i_qcr2000 == 1 || i_qcr2013 == 1 || i_qcr2013v == 1) {
                        denom = std::sqrt(ux1*ux1 + ux2*ux2 + ux3*ux3 +
                                          ux4*ux4 + ux5_v*ux5_v + ux6*ux6 +
                                          ux7*ux7 + ux8*ux8 + ux9*ux9) + 1.e-20;
                        o11 = 0.0;
                        o22 = 0.0;
                        o33 = 0.0;
                        o12 = 2.0*w12/denom;
                        o13 = 2.0*w13/denom;
                        o23 = 2.0*w23/denom;
                        t11p = t11 - 0.3*( o11*t11+o12*t12+o13*t13+
                                            o11*t11+o12*t12+o13*t13);
                        t22p = t22 - 0.3*(-o12*t12+o22*t22+o23*t23+
                                           -o12*t12+o22*t22+o23*t23);
                        t33p = t33 - 0.3*(-o13*t13-o23*t23+o33*t33+
                                           -o13*t13-o23*t23+o33*t33);
                        t12p = t12 - 0.3*( o11*t12+o12*t22+o13*t23+
                                           -o12*t11+o22*t12+o23*t13);
                        t13p = t13 - 0.3*( o11*t13+o12*t23+o13*t33+
                                           -o13*t11-o23*t12+o33*t13);
                        t23p = t23 - 0.3*(-o12*t13+o22*t23+o23*t33+
                                           -o13*t12-o23*t22+o33*t23);
                        t11 = t11p; t22 = t22p; t33 = t33p;
                        t12 = t12p; t13 = t13p; t23 = t23p;
                    }
                    if (i_qcr2013 == 1) {
                        xis2013 = s11star*s11star + s22star*s22star +
                                  s33star*s33star +
                                  2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        s_mod = std::sqrt(2.0*xis2013);
                        wis   = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        w_mod = std::sqrt(2.0*wis);
                        { double a1=s_mod, a2=1.2*w_mod; s_mod = ccomplex_ns::ccmin(a1,a2); }
                        t11 += ccr2*v3d*s_mod;
                        t22 += ccr2*v3d*s_mod;
                        t33 += ccr2*v3d*s_mod;
                    }
                    if (i_qcr2013v == 1) {
                        wis   = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        w_mod = std::sqrt(2.0*wis);
                        t11 += ccr2*v3d*w_mod;
                        t22 += ccr2*v3d*w_mod;
                        t33 += ccr2*v3d*w_mod;
                    }
                    t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+
                                         t(j+jkv,27)*t13)/t(j+jkv,15);
                    t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+
                                         t(j+jkv,27)*t23)/t(j+jkv,15);
                    t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+
                                         t(j+jkv,27)*t33)/t(j+jkv,15);
                    t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                        t(j+jkv,25)*(u*t11 + v*t12 + w_vel*t13)
                       +t(j+jkv,26)*(u*t12 + v*t22 + w_vel*t23)
                       +t(j+jkv,27)*(u*t13 + v*t23 + w_vel*t33));
                } // j
            } // k
        } // ipl
        for (l = 2; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                t(izz, l) = t(izz+js-1, l+26) - t(izz, l+26);
            }
        }
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            for (k = 1; k <= kdim1; k++) {
                jkv = (k-1)*jdim1*npl + (ipl-1)*jdim1 + 1;
                for (l = 2; l <= 5; l++) {
                    for (izz = 1; izz <= jdim1; izz++) {
                        res(izz, k, ii, l) = res(izz, k, ii, l) + t(izz+jkv-1, l);
                    }
                }
            }
        }
    } // end i_tauijs==1 branch


    // Add tau_ij's for RSMs:
    if (ivisc[2] >= 70) {
        re = (double)reue / (double)xmach;
        for (ipl = 1; ipl <= npl; ipl++) {
            ii   = i + ipl - 1;
            jkv2 = (ipl-1)*jdim1;
            for (k = 1; k <= kdim; k++) {
                jk  = (k-1)*jdim;
                jkv = (k-1)*npl*jdim1 + (ipl-1)*jdim1;
                if (k == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        rho   = 0.5*(q(j,k,ii,1)+qk0(j,ii,1,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,1,1)*wk0(j+jkv2,5);
                        u     = 0.5*(q(j,k,ii,2)+qk0(j,ii,2,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,2,1)*wk0(j+jkv2,5);
                        v     = 0.5*(q(j,k,ii,3)+qk0(j,ii,3,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,3,1)*wk0(j+jkv2,5);
                        w_vel = 0.5*(q(j,k,ii,4)+qk0(j,ii,4,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,4,1)*wk0(j+jkv2,5);
                        t11 = -rho*0.5*(zksav(j,k,ii,1)+tk0(j,ii,1,1))*re;
                        t22 = -rho*0.5*(zksav(j,k,ii,2)+tk0(j,ii,2,1))*re;
                        t33 = -rho*0.5*(zksav(j,k,ii,3)+tk0(j,ii,3,1))*re;
                        t12 = -rho*0.5*(zksav(j,k,ii,4)+tk0(j,ii,4,1))*re;
                        t23 = -rho*0.5*(zksav(j,k,ii,5)+tk0(j,ii,5,1))*re;
                        t13 = -rho*0.5*(zksav(j,k,ii,6)+tk0(j,ii,6,1))*re;
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+
                                             t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+
                                             t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+
                                             t(j+jkv,27)*t33)/t(j+jkv,15);
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(u*t11 + v*t12 + w_vel*t13)
                           +t(j+jkv,26)*(u*t12 + v*t22 + w_vel*t23)
                           +t(j+jkv,27)*(u*t13 + v*t23 + w_vel*t33));
                    }
                } else if (k == kdim) {
                    for (j = 1; j <= jdim1; j++) {
                        rho   = 0.5*(q(j,k-1,ii,1)+qk0(j,ii,1,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,1,3)*wk0(j+jkv2,15);
                        u     = 0.5*(q(j,k-1,ii,2)+qk0(j,ii,2,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,2,3)*wk0(j+jkv2,15);
                        v     = 0.5*(q(j,k-1,ii,3)+qk0(j,ii,3,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,3,3)*wk0(j+jkv2,15);
                        w_vel = 0.5*(q(j,k-1,ii,4)+qk0(j,ii,4,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,4,3)*wk0(j+jkv2,15);
                        t11 = -rho*0.5*(zksav(j,k-1,ii,1)+tk0(j,ii,1,3))*re;
                        t22 = -rho*0.5*(zksav(j,k-1,ii,2)+tk0(j,ii,2,3))*re;
                        t33 = -rho*0.5*(zksav(j,k-1,ii,3)+tk0(j,ii,3,3))*re;
                        t12 = -rho*0.5*(zksav(j,k-1,ii,4)+tk0(j,ii,4,3))*re;
                        t23 = -rho*0.5*(zksav(j,k-1,ii,5)+tk0(j,ii,5,3))*re;
                        t13 = -rho*0.5*(zksav(j,k-1,ii,6)+tk0(j,ii,6,3))*re;
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+
                                             t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+
                                             t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+
                                             t(j+jkv,27)*t33)/t(j+jkv,15);
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(u*t11 + v*t12 + w_vel*t13)
                           +t(j+jkv,26)*(u*t12 + v*t22 + w_vel*t23)
                           +t(j+jkv,27)*(u*t13 + v*t23 + w_vel*t33));
                    }
                } else {
                    for (j = 1; j <= jdim1; j++) {
                        rho = 0.5*(q(j,k,ii,1)+q(j,k-1,ii,1));
                        t11 = -rho*0.5*(zksav(j,k,ii,1)+zksav(j,k-1,ii,1))*re;
                        t22 = -rho*0.5*(zksav(j,k,ii,2)+zksav(j,k-1,ii,2))*re;
                        t33 = -rho*0.5*(zksav(j,k,ii,3)+zksav(j,k-1,ii,3))*re;
                        t12 = -rho*0.5*(zksav(j,k,ii,4)+zksav(j,k-1,ii,4))*re;
                        t23 = -rho*0.5*(zksav(j,k,ii,5)+zksav(j,k-1,ii,5))*re;
                        t13 = -rho*0.5*(zksav(j,k,ii,6)+zksav(j,k-1,ii,6))*re;
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+
                                             t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+
                                             t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+
                                             t(j+jkv,27)*t33)/t(j+jkv,15);
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t11
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t12
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t13)
                           +t(j+jkv,26)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t12
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t22
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t23)
                           +t(j+jkv,27)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t13
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t23
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t33));
                    }
                }
            } // k
        } // ipl
        for (l = 2; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                t(izz, l) = t(izz+js-1, l+26) - t(izz, l+26);
            }
        }
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            for (k = 1; k <= kdim1; k++) {
                jkv = (k-1)*jdim1*npl + (ipl-1)*jdim1 + 1;
                for (l = 2; l <= 5; l++) {
                    for (izz = 1; izz <= jdim1; izz++) {
                        res(izz, k, ii, l) = res(izz, k, ii, l) + t(izz+jkv-1, l);
                    }
                }
            }
        }
    } // end ivisc[2]>=70 branch


    // EASM v1: ivisc(3)==11 or 12
    if (i_tauijs != 1 && (ivisc[2] == 11 || ivisc[2] == 12) && level >= lglobal) {
        // f7=factor for determining whether it's a k-omega or k-epsilon formulation
        // f7=0 for k-epsilon, 1 for k-omega
        f7 = 0.0;
        if (ivisc[2] == 12) f7 = 1.0;
        re = (double)reue / (double)xmach;
        gg = 1.0 / ((double)c10_e + (double)c5_e - 1.0);
        // Add nonlinear terms to RHS
        for (ipl = 1; ipl <= npl; ipl++) {
            ii   = i + ipl - 1;
            jkv2 = (ipl-1)*jdim1;
            for (k = 1; k <= kdim; k++) {
                jkv = (k-1)*npl*jdim1 + (ipl-1)*jdim1;
                if (k == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        // Determine Sij and Wij values (at cell interface)
                        s11 = t(j+jkv,8) *t(j+jkv,25);
                        s22 = t(j+jkv,9) *t(j+jkv,26);
                        s33 = t(j+jkv,10)*t(j+jkv,27);
                        s12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) + t(j+jkv,9)*t(j+jkv,25));
                        s13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,25));
                        s23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,26));
                        w12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) - t(j+jkv,9)*t(j+jkv,25));
                        w13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,25));
                        w23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,26));
                        xis = s11*s11 + s22*s22 + s33*s33 +
                              2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        // "Constants" are function of c3, c4, and gg:
                        alpa1 = (2.0-(double)c4_e)/2.0*gg;
                        alpa2 = (2.0-(double)c3_e)*gg;
                        // Find tauij values:
                        rho   = 0.5*(q(j,k,ii,1)+qk0(j,ii,1,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,1,1)*wk0(j+jkv2,5);
                        u     = 0.5*(q(j,k,ii,2)+qk0(j,ii,2,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,2,1)*wk0(j+jkv2,5);
                        v     = 0.5*(q(j,k,ii,3)+qk0(j,ii,3,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,3,1)*wk0(j+jkv2,5);
                        w_vel = 0.5*(q(j,k,ii,4)+qk0(j,ii,4,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,4,1)*wk0(j+jkv2,5);
                        zksav1 = 0.5*(zksav(j,k,ii,1)+tk0(j,ii,1,1));
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        zksav2 = 0.5*(zksav(j,k,ii,2)+tk0(j,ii,2,1));
                        bnum   = zksav2*(1.0-f7)+f7;
                        eta    = (2.0-(double)c3_e)*(2.0-(double)c3_e)*(gg*gg/4.0)*xis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        squig  = (2.0-(double)c4_e)*(2.0-(double)c4_e)*(gg*gg/4.0)*wis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        { double a1=eta,   a2=10.0; eta   = ccomplex_ns::ccmincr(a1,a2); }
                        { double a1=squig, a2=10.0; squig = ccomplex_ns::ccmincr(a1,a2); }
                        factre = (3.0*(1.0+eta)+.2e-8*(eta*eta*eta+squig*squig*squig)) /
                                 (3.0*(1.0+eta)+   .2*(eta*eta*eta+squig*squig*squig));
                        t11 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(-s12*w12-s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t11 = t11*(1.0-bck(j,ii,1));
                        t22 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s12*w12-s23*w23)/re
                             +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t22 = t22*(1.0-bck(j,ii,1));
                        t33 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s23*w23+s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t33 = t33*(1.0-bck(j,ii,1));
                        t12 = -2.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t12 = t12*(1.0-bck(j,ii,1));
                        t13 = -2.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t13 = t13*(1.0-bck(j,ii,1));
                        t23 = -2.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                              +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t23 = t23*(1.0-bck(j,ii,1));
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+t(j+jkv,27)*t33)/t(j+jkv,15);
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(u*t11 + v*t12 + w_vel*t13)
                           +t(j+jkv,26)*(u*t12 + v*t22 + w_vel*t23)
                           +t(j+jkv,27)*(u*t13 + v*t23 + w_vel*t33))
                           -xmre/t(j+jkv,15)*(t(j+jkv,14)*t(j+jkv,15)/xmre
                            -t(j+jkv,24)+(double)sigk1_e*t(j+jkv,24))*
                            (t(j+jkv,25)*t(j+jkv,25)+t(j+jkv,26)*t(j+jkv,26)+t(j+jkv,27)*t(j+jkv,27))*
                            (zksav(j,k,ii,2)-tk0(j,ii,2,1));
                    } // j
                } else if (k == kdim) {

                    for (j = 1; j <= jdim1; j++) {
                        // Determine Sij and Wij values (at cell interface)
                        s11 = t(j+jkv,8) *t(j+jkv,25);
                        s22 = t(j+jkv,9) *t(j+jkv,26);
                        s33 = t(j+jkv,10)*t(j+jkv,27);
                        s12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) + t(j+jkv,9)*t(j+jkv,25));
                        s13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,25));
                        s23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,26));
                        w12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) - t(j+jkv,9)*t(j+jkv,25));
                        w13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,25));
                        w23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,26));
                        xis = s11*s11 + s22*s22 + s33*s33 +
                              2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        alpa1 = (2.0-(double)c4_e)/2.0*gg;
                        alpa2 = (2.0-(double)c3_e)*gg;
                        rho   = 0.5*(q(j,k-1,ii,1)+qk0(j,ii,1,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,1,3)*wk0(j+jkv2,15);
                        u     = 0.5*(q(j,k-1,ii,2)+qk0(j,ii,2,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,2,3)*wk0(j+jkv2,15);
                        v     = 0.5*(q(j,k-1,ii,3)+qk0(j,ii,3,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,3,3)*wk0(j+jkv2,15);
                        w_vel = 0.5*(q(j,k-1,ii,4)+qk0(j,ii,4,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,4,3)*wk0(j+jkv2,15);
                        zksav1 = 0.5*(zksav(j,k-1,ii,1)+tk0(j,ii,1,3));
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        zksav2 = 0.5*(zksav(j,k-1,ii,2)+tk0(j,ii,2,3));
                        bnum   = zksav2*(1.0-f7)+f7;
                        eta    = (2.0-(double)c3_e)*(2.0-(double)c3_e)*(gg*gg/4.0)*xis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        squig  = (2.0-(double)c4_e)*(2.0-(double)c4_e)*(gg*gg/4.0)*wis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        { double a1=eta,   a2=10.0; eta   = ccomplex_ns::ccmincr(a1,a2); }
                        { double a1=squig, a2=10.0; squig = ccomplex_ns::ccmincr(a1,a2); }
                        factre = (3.0*(1.0+eta)+.2e-8*(eta*eta*eta+squig*squig*squig)) /
                                 (3.0*(1.0+eta)+   .2*(eta*eta*eta+squig*squig*squig));
                        t11 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(-s12*w12-s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t11 = t11*(1.0-bck(j,ii,2));
                        t22 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s12*w12-s23*w23)/re
                             +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t22 = t22*(1.0-bck(j,ii,2));
                        t33 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s23*w23+s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t33 = t33*(1.0-bck(j,ii,2));
                        t12 = -2.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t12 = t12*(1.0-bck(j,ii,2));
                        t13 = -2.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t13 = t13*(1.0-bck(j,ii,2));
                        t23 = -2.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                              +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t23 = t23*(1.0-bck(j,ii,2));
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+t(j+jkv,27)*t33)/t(j+jkv,15);
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(u*t11 + v*t12 + w_vel*t13)
                           +t(j+jkv,26)*(u*t12 + v*t22 + w_vel*t23)
                           +t(j+jkv,27)*(u*t13 + v*t23 + w_vel*t33))
                           -xmre/t(j+jkv,15)*(t(j+jkv,14)*t(j+jkv,15)/xmre
                            -t(j+jkv,24)+(double)sigk1_e*t(j+jkv,24))*
                            (t(j+jkv,25)*t(j+jkv,25)+t(j+jkv,26)*t(j+jkv,26)+t(j+jkv,27)*t(j+jkv,27))*
                            (tk0(j,ii,2,3)-zksav(j,k-1,ii,2));
                    } // j
                } else {

                    // interior interfaces (k != 1 and k != kdim)
                    for (j = 1; j <= jdim1; j++) {
                        s11 = t(j+jkv,8) *t(j+jkv,25);
                        s22 = t(j+jkv,9) *t(j+jkv,26);
                        s33 = t(j+jkv,10)*t(j+jkv,27);
                        s12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) + t(j+jkv,9)*t(j+jkv,25));
                        s13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,25));
                        s23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,26));
                        w12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) - t(j+jkv,9)*t(j+jkv,25));
                        w13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,25));
                        w23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,26));
                        xis = s11*s11 + s22*s22 + s33*s33 +
                              2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        alpa1 = (2.0-(double)c4_e)/2.0*gg;
                        alpa2 = (2.0-(double)c3_e)*gg;
                        rho    = 0.5*(q(j,k,ii,1)+q(j,k-1,ii,1));
                        zksav1 = 0.5*(zksav(j,k,ii,1)+zksav(j,k-1,ii,1));
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        zksav2 = 0.5*(zksav(j,k,ii,2)+zksav(j,k-1,ii,2));
                        bnum   = zksav2*(1.0-f7)+f7;
                        eta    = (2.0-(double)c3_e)*(2.0-(double)c3_e)*(gg*gg/4.0)*xis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        squig  = (2.0-(double)c4_e)*(2.0-(double)c4_e)*(gg*gg/4.0)*wis*(bnum/(zksav1*re))*(bnum/(zksav1*re));
                        { double a1=eta,   a2=10.0; eta   = ccomplex_ns::ccmincr(a1,a2); }
                        { double a1=squig, a2=10.0; squig = ccomplex_ns::ccmincr(a1,a2); }
                        factre = (3.0*(1.0+eta)+.2e-8*(eta*eta*eta+squig*squig*squig)) /
                                 (3.0*(1.0+eta)+   .2*(eta*eta*eta+squig*squig*squig));
                        t11 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(-s12*w12-s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t22 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s12*w12-s23*w23)/re
                             +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t33 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s23*w23+s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t12 = -2.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t13 = -2.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t23 = -2.0*alpa1*t(j+jkv,24)*factre*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                              +2.0*alpa2*t(j+jkv,24)*factre*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+t(j+jkv,27)*t33)/t(j+jkv,15);
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t11
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t12
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t13)
                           +t(j+jkv,26)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t12
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t22
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t23)
                           +t(j+jkv,27)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t13
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t23
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t33))
                           -xmre/t(j+jkv,15)*(t(j+jkv,14)*t(j+jkv,15)/xmre
                            -t(j+jkv,24)+(double)sigk1_e*t(j+jkv,24))*
                            (t(j+jkv,25)*t(j+jkv,25)+t(j+jkv,26)*t(j+jkv,26)+t(j+jkv,27)*t(j+jkv,27))*
                            (zksav(j,k,ii,2)-zksav(j,k-1,ii,2));
                    } // j
                } // k branches
            } // k
        } // ipl
    } // end EASM v1 branch


    // EASM v2 (Girimaji): ivisc(3)==13 or 14
    else if (i_tauijs != 1 && (ivisc[2] == 13 || ivisc[2] == 14) && level >= lglobal) {
        f7 = 0.0;
        if (ivisc[2] == 14) f7 = 1.0;
        re = (double)reue / (double)xmach;
        al10 = 0.5*(double)c10_e - 1.0;
        al1  = 2.0*(0.5*(double)c11_e + 1.0);
        // 2-line improvement for wake
        if (ieasm_type_e == 0 || ieasm_type_e == 3 || ieasm_type_e == 4) {
            al10 = al10 + 1.8864;
            al1  = al1  - 2.0;
        }
        al2 = 0.5*(double)c2_e - 2.0/3.0;
        al3 = 0.5*(double)c3_e - 1.0;
        al4 = 0.5*(double)c4_e - 1.0;
        for (ipl = 1; ipl <= npl; ipl++) {
            ii   = i + ipl - 1;
            jkv2 = (ipl-1)*jdim1;
            for (k = 1; k <= kdim; k++) {
                jkv = (k-1)*npl*jdim1 + (ipl-1)*jdim1;
                if (k == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        s11 = t(j+jkv,8) *t(j+jkv,25);
                        s22 = t(j+jkv,9) *t(j+jkv,26);
                        s33 = t(j+jkv,10)*t(j+jkv,27);
                        s12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) + t(j+jkv,9)*t(j+jkv,25));
                        s13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,25));
                        s23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,26));
                        w12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) - t(j+jkv,9)*t(j+jkv,25));
                        w13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,25));
                        w23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,26));
                        xis = s11*s11 + s22*s22 + s33*s33 +
                              2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        // Girimaji JFM 2000 fix to c4
                        if (ieasm_type_e == 4) {
                            wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                            eta1_girimaji = xis / (xis + wis);
                            if (eta1_girimaji < 0.5) {
                                c4new = 2.0 - ((2.0-(double)c4_e)*
                                        std::pow(eta1_girimaji/(1.0-eta1_girimaji), 0.75));
                            } else {
                                c4new = (double)c4_e;
                            }
                            al4 = 0.5*c4new - 1.0;
                        }
                        rho   = 0.5*(q(j,k,ii,1)+qk0(j,ii,1,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,1,1)*wk0(j+jkv2,5);
                        u     = 0.5*(q(j,k,ii,2)+qk0(j,ii,2,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,2,1)*wk0(j+jkv2,5);
                        v     = 0.5*(q(j,k,ii,3)+qk0(j,ii,3,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,3,1)*wk0(j+jkv2,5);
                        w_vel = 0.5*(q(j,k,ii,4)+qk0(j,ii,4,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,4,1)*wk0(j+jkv2,5);
                        zksav1 = 0.5*(zksav(j,k,ii,1)+tk0(j,ii,1,1));
                        zksav2 = 0.5*(zksav(j,k,ii,2)+tk0(j,ii,2,1));
                        cmuu   = cmuv(j,k,ii);
                        bnum   = zksav2*(1.0-f7)+f7;
                        // Durbin TCFD 1991 near-wall limiter
                        if (idurbinlim != 0 && (ieasm_type_e == 3 || ieasm_type_e == 4)) {
                            tau     = bnum/zksav1;
                            fmu     = t(j+jkv,14)*t(j+jkv,15)/xmre - t(j+jkv,24);
                            zkolmog = 6.0*std::sqrt(fmu*bnum/(rho*zksav2*zksav1));
                            { double a1=tau, a2=zkolmog; tau=ccomplex_ns::ccmax(a1,a2); }
                            zksav1  = bnum/tau;
                        }
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        eta1  = xis*(bnum/zksav1/re)*(bnum/zksav1/re);
                        alpa1 = -al4/(al10-eta1*al1*cmuu);
                        alpa2 = -2.0*al3/(al10-eta1*cmuu*al1);
                        t11 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*bnum/zksav1*(-s12*w12-s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t11 = t11*(1.0-bck(j,ii,1));
                        t22 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s12*w12-s23*w23)/re
                             +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t22 = t22*(1.0-bck(j,ii,1));
                        t33 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s23*w23+s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t33 = t33*(1.0-bck(j,ii,1));
                        t12 = -2.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t12 = t12*(1.0-bck(j,ii,1));
                        t13 = -2.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t13 = t13*(1.0-bck(j,ii,1));
                        t23 = -2.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                              +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t23 = t23*(1.0-bck(j,ii,1));
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+t(j+jkv,27)*t33)/t(j+jkv,15);
                        xmut = f7*t(j+jkv,24)+(1.0-f7)*(double)cmuc1_e*rho*zksav2*zksav2/zksav1;
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(u*t11 + v*t12 + w_vel*t13)
                           +t(j+jkv,26)*(u*t12 + v*t22 + w_vel*t23)
                           +t(j+jkv,27)*(u*t13 + v*t23 + w_vel*t33))
                           -xmre/t(j+jkv,15)*(t(j+jkv,14)*t(j+jkv,15)/xmre
                            -t(j+jkv,24)+(double)sigk1_e*xmut)*
                            (t(j+jkv,25)*t(j+jkv,25)+t(j+jkv,26)*t(j+jkv,26)+t(j+jkv,27)*t(j+jkv,27))*
                            (zksav(j,k,ii,2)-tk0(j,ii,2,1));
                    } // j
                } else if (k == kdim) {

                    for (j = 1; j <= jdim1; j++) {
                        s11 = t(j+jkv,8) *t(j+jkv,25);
                        s22 = t(j+jkv,9) *t(j+jkv,26);
                        s33 = t(j+jkv,10)*t(j+jkv,27);
                        s12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) + t(j+jkv,9)*t(j+jkv,25));
                        s13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,25));
                        s23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,26));
                        w12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) - t(j+jkv,9)*t(j+jkv,25));
                        w13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,25));
                        w23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,26));
                        xis = s11*s11 + s22*s22 + s33*s33 +
                              2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        // Girimaji JFM 2000 fix to c4
                        if (ieasm_type_e == 4) {
                            wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                            eta1_girimaji = xis / (xis + wis);
                            if (eta1_girimaji < 0.5) {
                                c4new = 2.0 - ((2.0-(double)c4_e)*
                                        std::pow(eta1_girimaji/(1.0-eta1_girimaji), 0.75));
                            } else {
                                c4new = (double)c4_e;
                            }
                            al4 = 0.5*c4new - 1.0;
                        }
                        rho   = 0.5*(q(j,k-1,ii,1)+qk0(j,ii,1,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,1,3)*wk0(j+jkv2,15);
                        u     = 0.5*(q(j,k-1,ii,2)+qk0(j,ii,2,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,2,3)*wk0(j+jkv2,15);
                        v     = 0.5*(q(j,k-1,ii,3)+qk0(j,ii,3,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,3,3)*wk0(j+jkv2,15);
                        w_vel = 0.5*(q(j,k-1,ii,4)+qk0(j,ii,4,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,4,3)*wk0(j+jkv2,15);
                        zksav1 = 0.5*(zksav(j,k-1,ii,1)+tk0(j,ii,1,3));
                        zksav2 = 0.5*(zksav(j,k-1,ii,2)+tk0(j,ii,2,3));
                        cmuu   = cmuv(j,k-1,ii);
                        bnum   = zksav2*(1.0-f7)+f7;
                        if (idurbinlim != 0 && (ieasm_type_e == 3 || ieasm_type_e == 4)) {
                            tau     = bnum/zksav1;
                            fmu     = t(j+jkv,14)*t(j+jkv,15)/xmre - t(j+jkv,24);
                            zkolmog = 6.0*std::sqrt(fmu*bnum/(rho*zksav2*zksav1));
                            { double a1=tau, a2=zkolmog; tau=ccomplex_ns::ccmax(a1,a2); }
                            zksav1  = bnum/tau;
                        }
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        eta1  = xis*(bnum/zksav1/re)*(bnum/zksav1/re);
                        alpa1 = -al4/(al10-eta1*al1*cmuu);
                        alpa2 = -2.0*al3/(al10-eta1*cmuu*al1);
                        t11 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*bnum/zksav1*(-s12*w12-s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t11 = t11*(1.0-bck(j,ii,2));
                        t22 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s12*w12-s23*w23)/re
                             +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t22 = t22*(1.0-bck(j,ii,2));
                        t33 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s23*w23+s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t33 = t33*(1.0-bck(j,ii,2));
                        t12 = -2.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t12 = t12*(1.0-bck(j,ii,2));
                        t13 = -2.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t13 = t13*(1.0-bck(j,ii,2));
                        t23 = -2.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                              +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t23 = t23*(1.0-bck(j,ii,2));
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+t(j+jkv,27)*t33)/t(j+jkv,15);
                        xmut = f7*t(j+jkv,24)+(1.0-f7)*(double)cmuc1_e*rho*zksav2*zksav2/zksav1;
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(u*t11 + v*t12 + w_vel*t13)
                           +t(j+jkv,26)*(u*t12 + v*t22 + w_vel*t23)
                           +t(j+jkv,27)*(u*t13 + v*t23 + w_vel*t33))
                           -xmre/t(j+jkv,15)*(t(j+jkv,14)*t(j+jkv,15)/xmre
                            -t(j+jkv,24)+(double)sigk1_e*xmut)*
                            (t(j+jkv,25)*t(j+jkv,25)+t(j+jkv,26)*t(j+jkv,26)+t(j+jkv,27)*t(j+jkv,27))*
                            (tk0(j,ii,2,3)-zksav(j,k-1,ii,2));
                    } // j
                } else {

                    // interior interfaces
                    for (j = 1; j <= jdim1; j++) {
                        s11 = t(j+jkv,8) *t(j+jkv,25);
                        s22 = t(j+jkv,9) *t(j+jkv,26);
                        s33 = t(j+jkv,10)*t(j+jkv,27);
                        s12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) + t(j+jkv,9)*t(j+jkv,25));
                        s13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,25));
                        s23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,26));
                        w12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) - t(j+jkv,9)*t(j+jkv,25));
                        w13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,25));
                        w23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,26));
                        xis = s11*s11 + s22*s22 + s33*s33 +
                              2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        // Girimaji JFM 2000 fix to c4
                        if (ieasm_type_e == 4) {
                            wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                            eta1_girimaji = xis / (xis + wis);
                            if (eta1_girimaji < 0.5) {
                                c4new = 2.0 - ((2.0-(double)c4_e)*
                                        std::pow(eta1_girimaji/(1.0-eta1_girimaji), 0.75));
                            } else {
                                c4new = (double)c4_e;
                            }
                            al4 = 0.5*c4new - 1.0;
                        }
                        rho    = 0.5*(q(j,k,ii,1)+q(j,k-1,ii,1));
                        zksav1 = 0.5*(zksav(j,k,ii,1)+zksav(j,k-1,ii,1));
                        zksav2 = 0.5*(zksav(j,k,ii,2)+zksav(j,k-1,ii,2));
                        cmuu   = 0.5*(cmuv(j,k,ii)+cmuv(j,k-1,ii));
                        bnum   = zksav2*(1.0-f7)+f7;
                        if (idurbinlim != 0 && (ieasm_type_e == 3 || ieasm_type_e == 4)) {
                            tau     = bnum/zksav1;
                            fmu     = t(j+jkv,14)*t(j+jkv,15)/xmre - t(j+jkv,24);
                            zkolmog = 6.0*std::sqrt(fmu*bnum/(rho*zksav2*zksav1));
                            { double a1=tau, a2=zkolmog; tau=ccomplex_ns::ccmax(a1,a2); }
                            zksav1  = bnum/tau;
                        }
                        { double a1=zksav1, a2=(double)tur10[0]; zksav1=ccomplex_ns::ccmax(a1,a2); }
                        eta1  = xis*(bnum/zksav1/re)*(bnum/zksav1/re);
                        alpa1 = -al4/(al10-eta1*al1*cmuu);
                        alpa2 = -2.0*al3/(al10-eta1*cmuu*al1);
                        t11 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*bnum/zksav1*(-s12*w12-s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s11-(s11+s22+s33)/3.0);
                        { double a1=t11, a2=xlnpt; t11=ccomplex_ns::ccmax(a1,a2); }
                        t22 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s12*w12-s23*w23)/re
                             +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s22-(s11+s22+s33)/3.0);
                        { double a1=t22, a2=xlnpt; t22=ccomplex_ns::ccmax(a1,a2); }
                        t33 = 2.0*rho*zksav2*re/3.0
                             -4.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s23*w23+s13*w13)/re
                             +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                        xlnpt = 2.0*t(j+jkv,24)*(s33-(s11+s22+s33)/3.0);
                        { double a1=t33, a2=xlnpt; t33=ccomplex_ns::ccmax(a1,a2); }
                        t12 = -2.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s11*s12+s12*s22+s13*s23)/re;
                        t13 = -2.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                              +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s11*s13+s12*s23+s13*s33)/re;
                        t23 = -2.0*alpa1*t(j+jkv,24)*bnum/zksav1*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                              +2.0*alpa2*t(j+jkv,24)*bnum/zksav1*(s12*s13+s22*s23+s23*s33)/re;
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+t(j+jkv,27)*t33)/t(j+jkv,15);
                        xmut = f7*t(j+jkv,24)+(1.0-f7)*(double)cmuc1_e*rho*zksav2*zksav2/zksav1;
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t11
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t12
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t13)
                           +t(j+jkv,26)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t12
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t22
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t23)
                           +t(j+jkv,27)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t13
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t23
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t33))
                           -xmre/t(j+jkv,15)*(t(j+jkv,14)*t(j+jkv,15)/xmre
                            -t(j+jkv,24)+(double)sigk1_e*xmut)*
                            (t(j+jkv,25)*t(j+jkv,25)+t(j+jkv,26)*t(j+jkv,26)+t(j+jkv,27)*t(j+jkv,27))*
                            (zksav(j,k,ii,2)-zksav(j,k-1,ii,2));
                    } // j
                } // k branches
            } // k
        } // ipl
    } // end EASM v2 branch


    // Generic nonlinear constitutive: i_nonlin != 0
    else if (i_tauijs != 1 && i_nonlin != 0 && level >= lglobal) {
        re = (double)reue / (double)xmach;
        for (ipl = 1; ipl <= npl; ipl++) {
            ii   = i + ipl - 1;
            jkv2 = (ipl-1)*jdim1;
            for (k = 1; k <= kdim; k++) {
                jkv = (k-1)*npl*jdim1 + (ipl-1)*jdim1;
                if (k == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        s11 = t(j+jkv,8) *t(j+jkv,25);
                        s22 = t(j+jkv,9) *t(j+jkv,26);
                        s33 = t(j+jkv,10)*t(j+jkv,27);
                        s12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) + t(j+jkv,9)*t(j+jkv,25));
                        s13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,25));
                        s23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,26));
                        w12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) - t(j+jkv,9)*t(j+jkv,25));
                        w13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,25));
                        w23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,26));
                        tmp = (s11+s22+s33)/3.0;
                        s11 = s11 - tmp;
                        s22 = s22 - tmp;
                        s33 = s33 - tmp;
                        // Note: denom term is taken at nearest cell center
                        denom = ux(j,k,ii,1)*ux(j,k,ii,1) + ux(j,k,ii,2)*ux(j,k,ii,2)
                               +ux(j,k,ii,3)*ux(j,k,ii,3) + ux(j,k,ii,4)*ux(j,k,ii,4)
                               +ux(j,k,ii,5)*ux(j,k,ii,5) + ux(j,k,ii,6)*ux(j,k,ii,6)
                               +ux(j,k,ii,7)*ux(j,k,ii,7) + ux(j,k,ii,8)*ux(j,k,ii,8)
                               +ux(j,k,ii,9)*ux(j,k,ii,9);
                        denom = std::sqrt(denom);
                        { double a1=denom, a2=(double)snonlin_lim; denom=ccomplex_ns::ccmax(a1,a2); }
                        rho   = 0.5*(q(j,k,ii,1)+qk0(j,ii,1,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,1,1)*wk0(j+jkv2,5);
                        u     = 0.5*(q(j,k,ii,2)+qk0(j,ii,2,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,2,1)*wk0(j+jkv2,5);
                        v     = 0.5*(q(j,k,ii,3)+qk0(j,ii,3,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,3,1)*wk0(j+jkv2,5);
                        w_vel = 0.5*(q(j,k,ii,4)+qk0(j,ii,4,1))*(1.-wk0(j+jkv2,5))+
                                qk0(j,ii,4,1)*wk0(j+jkv2,5);
                        zk_sa = 0.0;
                        if (ivisc[2]==6 || ivisc[2]==7 || ivisc[2]==10 || ivisc[2]==15) {
                            zk_sa = 0.5*(zksav(j,k,ii,2)+tk0(j,ii,2,1));
                        }
                        t11 = 2.0*rho*zk_sa*re/3.0
                             -8.0*(double)c_nonlin*t(j+jkv,24)/denom*(-s12*w12-s13*w13);
                        t11 = t11*(1.0-bck(j,ii,1));
                        t22 = 2.0*rho*zk_sa*re/3.0
                             -8.0*(double)c_nonlin*t(j+jkv,24)/denom*(s12*w12-s23*w23);
                        t22 = t22*(1.0-bck(j,ii,1));
                        t33 = 2.0*rho*zk_sa*re/3.0
                             -8.0*(double)c_nonlin*t(j+jkv,24)/denom*(s23*w23+s13*w13);
                        t33 = t33*(1.0-bck(j,ii,1));
                        t12 = -4.0*(double)c_nonlin*t(j+jkv,24)/denom*
                              (s11*w12-s22*w12-s13*w23-s23*w13);
                        t12 = t12*(1.0-bck(j,ii,1));
                        t13 = -4.0*(double)c_nonlin*t(j+jkv,24)/denom*
                              (s11*w13+s12*w23-s23*w12-s33*w13);
                        t13 = t13*(1.0-bck(j,ii,1));
                        t23 = -4.0*(double)c_nonlin*t(j+jkv,24)/denom*
                              (s12*w13+s22*w23+s13*w12-s33*w23);
                        t23 = t23*(1.0-bck(j,ii,1));
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+t(j+jkv,27)*t33)/t(j+jkv,15);
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(u*t11 + v*t12 + w_vel*t13)
                           +t(j+jkv,26)*(u*t12 + v*t22 + w_vel*t23)
                           +t(j+jkv,27)*(u*t13 + v*t23 + w_vel*t33));
                    } // j
                } else if (k == kdim) {
                    for (j = 1; j <= jdim1; j++) {
                        s11 = t(j+jkv,8) *t(j+jkv,25);
                        s22 = t(j+jkv,9) *t(j+jkv,26);
                        s33 = t(j+jkv,10)*t(j+jkv,27);
                        s12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) + t(j+jkv,9)*t(j+jkv,25));
                        s13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,25));
                        s23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,26));
                        w12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) - t(j+jkv,9)*t(j+jkv,25));
                        w13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,25));
                        w23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,26));
                        tmp = (s11+s22+s33)/3.0;
                        s11 = s11 - tmp; s22 = s22 - tmp; s33 = s33 - tmp;
                        // Note: denom term is taken at nearest cell center
                        denom = ux(j,k-1,ii,1)*ux(j,k-1,ii,1) + ux(j,k-1,ii,2)*ux(j,k-1,ii,2)
                               +ux(j,k-1,ii,3)*ux(j,k-1,ii,3) + ux(j,k-1,ii,4)*ux(j,k-1,ii,4)
                               +ux(j,k-1,ii,5)*ux(j,k-1,ii,5) + ux(j,k-1,ii,6)*ux(j,k-1,ii,6)
                               +ux(j,k-1,ii,7)*ux(j,k-1,ii,7) + ux(j,k-1,ii,8)*ux(j,k-1,ii,8)
                               +ux(j,k-1,ii,9)*ux(j,k-1,ii,9);
                        denom = std::sqrt(denom);
                        { double a1=denom, a2=(double)snonlin_lim; denom=ccomplex_ns::ccmax(a1,a2); }
                        rho   = 0.5*(q(j,k-1,ii,1)+qk0(j,ii,1,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,1,3)*wk0(j+jkv2,15);
                        u     = 0.5*(q(j,k-1,ii,2)+qk0(j,ii,2,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,2,3)*wk0(j+jkv2,15);
                        v     = 0.5*(q(j,k-1,ii,3)+qk0(j,ii,3,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,3,3)*wk0(j+jkv2,15);
                        w_vel = 0.5*(q(j,k-1,ii,4)+qk0(j,ii,4,3))*(1.-wk0(j+jkv2,15))+
                                qk0(j,ii,4,3)*wk0(j+jkv2,15);
                        zk_sa = 0.0;
                        if (ivisc[2]==6 || ivisc[2]==7 || ivisc[2]==10 || ivisc[2]==15) {
                            zk_sa = 0.5*(zksav(j,k-1,ii,2)+tk0(j,ii,2,3));
                        }
                        t11 = 2.0*rho*zk_sa*re/3.0
                             -8.0*(double)c_nonlin*t(j+jkv,24)/denom*(-s12*w12-s13*w13);
                        t11 = t11*(1.0-bck(j,ii,2));
                        t22 = 2.0*rho*zk_sa*re/3.0
                             -8.0*(double)c_nonlin*t(j+jkv,24)/denom*(s12*w12-s23*w23);
                        t22 = t22*(1.0-bck(j,ii,2));
                        t33 = 2.0*rho*zk_sa*re/3.0
                             -8.0*(double)c_nonlin*t(j+jkv,24)/denom*(s23*w23+s13*w13);
                        t33 = t33*(1.0-bck(j,ii,2));
                        t12 = -4.0*(double)c_nonlin*t(j+jkv,24)/denom*
                              (s11*w12-s22*w12-s13*w23-s23*w13);
                        t12 = t12*(1.0-bck(j,ii,2));
                        t13 = -4.0*(double)c_nonlin*t(j+jkv,24)/denom*
                              (s11*w13+s12*w23-s23*w12-s33*w13);
                        t13 = t13*(1.0-bck(j,ii,2));
                        t23 = -4.0*(double)c_nonlin*t(j+jkv,24)/denom*
                              (s12*w13+s22*w23+s13*w12-s33*w23);
                        t23 = t23*(1.0-bck(j,ii,2));
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+t(j+jkv,27)*t33)/t(j+jkv,15);
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(u*t11 + v*t12 + w_vel*t13)
                           +t(j+jkv,26)*(u*t12 + v*t22 + w_vel*t23)
                           +t(j+jkv,27)*(u*t13 + v*t23 + w_vel*t33));
                    } // j
                } else {

                    // interior interfaces
                    for (j = 1; j <= jdim1; j++) {
                        s11 = t(j+jkv,8) *t(j+jkv,25);
                        s22 = t(j+jkv,9) *t(j+jkv,26);
                        s33 = t(j+jkv,10)*t(j+jkv,27);
                        s12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) + t(j+jkv,9)*t(j+jkv,25));
                        s13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,25));
                        s23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) + t(j+jkv,10)*t(j+jkv,26));
                        w12 = 0.5*(t(j+jkv,8)*t(j+jkv,26) - t(j+jkv,9)*t(j+jkv,25));
                        w13 = 0.5*(t(j+jkv,8)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,25));
                        w23 = 0.5*(t(j+jkv,9)*t(j+jkv,27) - t(j+jkv,10)*t(j+jkv,26));
                        tmp = (s11+s22+s33)/3.0;
                        s11 = s11 - tmp; s22 = s22 - tmp; s33 = s33 - tmp;
                        // Note: denom term is averaged to cell face
                        ux1   = 0.5*(ux(j,k-1,ii,1)+ux(j,k,ii,1));
                        ux2   = 0.5*(ux(j,k-1,ii,2)+ux(j,k,ii,2));
                        ux3   = 0.5*(ux(j,k-1,ii,3)+ux(j,k,ii,3));
                        ux4   = 0.5*(ux(j,k-1,ii,4)+ux(j,k,ii,4));
                        ux5_v = 0.5*(ux(j,k-1,ii,5)+ux(j,k,ii,5));
                        ux6   = 0.5*(ux(j,k-1,ii,6)+ux(j,k,ii,6));
                        ux7   = 0.5*(ux(j,k-1,ii,7)+ux(j,k,ii,7));
                        ux8   = 0.5*(ux(j,k-1,ii,8)+ux(j,k,ii,8));
                        ux9   = 0.5*(ux(j,k-1,ii,9)+ux(j,k,ii,9));
                        denom = ux1*ux1 + ux2*ux2 + ux3*ux3
                               +ux4*ux4 + ux5_v*ux5_v + ux6*ux6
                               +ux7*ux7 + ux8*ux8 + ux9*ux9;
                        denom = std::sqrt(denom);
                        { double a1=denom, a2=(double)snonlin_lim; denom=ccomplex_ns::ccmax(a1,a2); }
                        rho   = 0.5*(q(j,k,ii,1)+q(j,k-1,ii,1));
                        zk_sa = 0.0;
                        if (ivisc[2]==6 || ivisc[2]==7 || ivisc[2]==10 || ivisc[2]==15) {
                            zk_sa = 0.5*(zksav(j,k,ii,2)+zksav(j,k-1,ii,2));
                        }
                        t11 = 2.0*rho*zk_sa*re/3.0
                             -8.0*(double)c_nonlin*t(j+jkv,24)/denom*(-s12*w12-s13*w13);
                        t22 = 2.0*rho*zk_sa*re/3.0
                             -8.0*(double)c_nonlin*t(j+jkv,24)/denom*(s12*w12-s23*w23);
                        t33 = 2.0*rho*zk_sa*re/3.0
                             -8.0*(double)c_nonlin*t(j+jkv,24)/denom*(s23*w23+s13*w13);
                        t12 = -4.0*(double)c_nonlin*t(j+jkv,24)/denom*
                              (s11*w12-s22*w12-s13*w23-s23*w13);
                        t13 = -4.0*(double)c_nonlin*t(j+jkv,24)/denom*
                              (s11*w13+s12*w23-s23*w12-s33*w13);
                        t23 = -4.0*(double)c_nonlin*t(j+jkv,24)/denom*
                              (s12*w13+s22*w23+s13*w12-s33*w23);
                        t(j+jkv, 28) = xmre*(t(j+jkv,25)*t11+t(j+jkv,26)*t12+t(j+jkv,27)*t13)/t(j+jkv,15);
                        t(j+jkv, 29) = xmre*(t(j+jkv,25)*t12+t(j+jkv,26)*t22+t(j+jkv,27)*t23)/t(j+jkv,15);
                        t(j+jkv, 30) = xmre*(t(j+jkv,25)*t13+t(j+jkv,26)*t23+t(j+jkv,27)*t33)/t(j+jkv,15);
                        t(j+jkv, 31) = xmre/t(j+jkv,15)*(
                            t(j+jkv,25)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t11
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t12
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t13)
                           +t(j+jkv,26)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t12
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t22
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t23)
                           +t(j+jkv,27)*(0.5*(q(j,k,ii,2)+q(j,k-1,ii,2))*t13
                                        +0.5*(q(j,k,ii,3)+q(j,k-1,ii,3))*t23
                                        +0.5*(q(j,k,ii,4)+q(j,k-1,ii,4))*t33));
                    } // j
                } // k branches
            } // k
        } // ipl
    } // end generic nonlinear branch

    // Final residual update for nonlinear branches (11,12,13,14,i_nonlin)
    if (i_tauijs != 1 &&
        (ivisc[2]==11 || ivisc[2]==12 || ivisc[2]==13 || ivisc[2]==14 || i_nonlin!=0) &&
        level >= lglobal) {
        for (l = 2; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                // Take positive derivative, because want negative of all above for correct
                // sign (this term is going into -d/dxj(Hv) term)
                t(izz, l) = t(izz+js-1, l+26) - t(izz, l+26);
            }
        }
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            for (k = 1; k <= kdim1; k++) {
                jkv = (k-1)*jdim1*npl + (ipl-1)*jdim1 + 1;
                for (l = 2; l <= 5; l++) {
                    for (izz = 1; izz <= jdim1; izz++) {
                        res(izz, k, ii, l) = res(izz, k, ii, l) + t(izz+jkv-1, l);
                    }
                }
            }
        }
    }
    return;


    // implicit matrix terms (idf > 0 path)
    label_3333:

    gterm1 = 1.0 / (double)pr;
    gterm2 = (double)gamma / (double)gm1;

    for (izz = 1; izz <= l0; izz++) {
        t(izz, 6) = tr1 * t(izz, 14) * t(izz, 25);
        t(izz, 7) = tr1 * t(izz, 14) * t(izz, 26);
        t(izz, 8) = tr1 * t(izz, 14) * t(izz, 27);
    }

    // momentum equations
    ns = n - jdim1*npl;
    n0 = n;
    for (l = 2; l <= 4; l++) {
        for (izz = 1; izz <= l0; izz++) {
            t(izz, 13) = t(izz, 6) * t(izz, l+23);
            t(izz, 14) = t(izz, 7) * t(izz, l+23);
            t(izz, 15) = t(izz, 8) * t(izz, l+23);
        }
        for (izz = 1; izz <= n0; izz++) {
            bk(izz, 1, l, l) = bk(izz, 1, l, l) + t(izz+js-1, 5) + t(izz, 5);
            ck(izz, 1, l, l) = ck(izz, 1, l, l) - t(izz+js-1, 5);
        }
        for (izz = 1; izz <= ns; izz++) {
            ak(izz, 2, l, l) = ak(izz, 2, l, l) - t(izz+js-1, 5);
        }
        for (k = 2; k <= 4; k++) {
            for (izz = 1; izz <= n0; izz++) {
                bk(izz, 1, l, k) = bk(izz, 1, l, k) + t(izz+js-1, 11+k) + t(izz, 11+k);
                ck(izz, 1, l, k) = ck(izz, 1, l, k) - t(izz+js-1, 11+k);
            }
            for (izz = 1; izz <= ns; izz++) {
                ak(izz, 2, l, k) = ak(izz, 2, l, k) - t(izz+js-1, 11+k);
            }
        }
    }

    // energy equation - (phi)1 term
    for (k = 2; k <= 4; k++) {
        for (izz = 1; izz <= n0; izz++) {
            bk(izz, 1, 5, k) = bk(izz, 1, 5, k) + (t(izz+js-1,5)+t(izz,5))*t(izz,16+k);
            ck(izz, 1, 5, k) = ck(izz, 1, 5, k) - (t(izz+js-1,5))*t(izz+js-1,16+k);
        }
        for (izz = 1; izz <= ns; izz++) {
            ak(izz, 2, 5, k) = ak(izz, 2, 5, k) - (t(izz+js-1,5))*t(izz,16+k);
        }
    }

    // energy equation - heat transfer terms
    for (izz = 1; izz <= l0; izz++) {
        t(izz, 5) = gterm1 * t(izz, 23) * t(izz, 5);
    }

    for (k = 1; k <= 5; k += 4) {
        if (k == 1) {
            for (izz = 1; izz <= n0; izz++) {
                t(izz, 13) = -t(izz, 16) * (double)pr * t(izz, 1);
            }
        } else {
            for (izz = 1; izz <= n0; izz++) {
                t(izz, 13) = gterm2 * t(izz, 1);
            }
        }
        for (izz = n0+1; izz <= l0; izz++) {
            t(izz, 13) = t(izz-js+1, 13);
        }
        for (izz = 1; izz <= n0; izz++) {
            bk(izz, 1, 5, k) = bk(izz, 1, 5, k) + (t(izz+js-1,5)+t(izz,5))*t(izz,13);
            ck(izz, 1, 5, k) = ck(izz, 1, 5, k) - (t(izz+js-1,5))*t(izz+js-1,13);
        }
        for (izz = 1; izz <= ns; izz++) {
            ak(izz, 2, 5, k) = ak(izz, 2, 5, k) - (t(izz+js-1,5))*t(izz,13);
        }
    }

    // energy equation - (phi)2 terms
    for (izz = 1; izz <= l0; izz++) {
        t(izz, 13) = t(izz, 6) * t(izz, 12);
        t(izz, 14) = t(izz, 7) * t(izz, 12);
        t(izz, 15) = t(izz, 8) * t(izz, 12);
    }

    for (k = 2; k <= 4; k++) {
        for (izz = 1; izz <= n0; izz++) {
            bk(izz, 1, 5, k) = bk(izz, 1, 5, k) + t(izz+js-1, 11+k) + t(izz, 11+k);
            ck(izz, 1, 5, k) = ck(izz, 1, 5, k) - t(izz+js-1, 11+k);
        }
        for (izz = 1; izz <= ns; izz++) {
            ak(izz, 2, 5, k) = ak(izz, 2, 5, k) - t(izz+js-1, 11+k);
        }
    }

    for (izz = 1; izz <= l0; izz++) {
        t(izz, 13) = 0.5e0 * t(izz, 25) * t(izz, 11);
        t(izz, 14) = 0.5e0 * t(izz, 26) * t(izz, 11);
        t(izz, 15) = 0.5e0 * t(izz, 27) * t(izz, 11);
    }

    for (k = 2; k <= 4; k++) {
        for (izz = 1; izz <= n0; izz++) {
            bk(izz, 1, 5, k) = bk(izz, 1, 5, k) - t(izz+js-1, 11+k) + t(izz, 11+k);
            ck(izz, 1, 5, k) = ck(izz, 1, 5, k) - t(izz+js-1, 11+k);
        }
        for (izz = 1; izz <= ns; izz++) {
            ak(izz, 2, 5, k) = ak(izz, 2, 5, k) + t(izz+js-1, 11+k);
        }
    }
    return;
}

} // namespace hfluxv_ns
