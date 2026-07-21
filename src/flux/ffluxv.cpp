// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "ffluxv.h"
#include "ccomplex.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

namespace ffluxv_ns {

void ffluxv(int& k, int& npl, int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> res, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qi0, FortranArray4DRef<double> si,
            FortranArray3DRef<double> vol, FortranArray2DRef<double> t,
            int& nvtq, FortranArray2DRef<double> wi0,
            FortranArray3DRef<double> vist3d, FortranArray2DRef<double> vmui,
            FortranArray4DRef<double> vi0, FortranArray3DRef<double> bci,
            FortranArray4DRef<double> zksav, FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> cmuv, FortranArray3DRef<double> voli0,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& iadv, int& nummem,
            FortranArray4DRef<double> ux)
{
    // COMMON block aliases
    float& gamma   = cmn_fluid.gamma;
    float& gm1     = cmn_fluid.gm1;
    float& pr      = cmn_fluid2.pr;
    float& prt     = cmn_fluid2.prt;
    float& cbar    = cmn_fluid2.cbar;
    float& xmach   = cmn_info.xmach;
    float& reue    = cmn_reyue.reue;
    float& tinf    = cmn_reyue.tinf;
    int*   ivisc   = cmn_reyue.ivisc;   // 0-based: ivisc[0]=ivisc(1)
    int&   isklton = cmn_sklton.isklton;
    int&   level   = cmn_mgrd.level;
    int&   lglobal = cmn_mgrd.lglobal;
    float& c10     = cmn_easmv.c10;
    float& c11     = cmn_easmv.c11;
    float& c2_e    = cmn_easmv.c2;
    float& c3_e    = cmn_easmv.c3;
    float& c4_e    = cmn_easmv.c4;
    float& c5_e    = cmn_easmv.c5;
    float& sigk1   = cmn_easmv.sigk1;
    float& cmuc1   = cmn_easmv.cmuc1;
    int&   ieasm_type = cmn_easmv.ieasm_type;
    int&   i_nonlin   = cmn_constit.i_nonlin;
    float& c_nonlin   = cmn_constit.c_nonlin;
    float& snonlin_lim = cmn_constit.snonlin_lim;
    int&   i_tauijs   = cmn_constit.i_tauijs;
    int&   i_qcr2000  = cmn_constit.i_qcr2000;
    int&   i_qcr2013  = cmn_constit.i_qcr2013;
    int&   i_qcr2013v = cmn_constit.i_qcr2013v;
    float* tur10      = cmn_ivals.tur10;  // 0-based: tur10[0]=tur10(1)

    // Local scalars
    float  coef_eddy;
    double ccr2;
    int    idim1, jdim1;
    int    n, l0, jv, js, nn;
    double xmre, gpr, gm1pr, prtr;
    int    idurbinlim;
    int    l1, kpl, kk, i, jiv, j;
    double t5, t6;
    int    mm, mb, jmin, jmax, kmin, kmax, jc;
    int    m;
    int    n1;
    int    ji;
    int    kq, k1, k2;
    int    iviscc;
    double t25, t24;
    double tr1;
    int    izz;
    double ab, bb, wi05, wi06;
    double c2b, c2bp;
    int    jiv2;
    double re, gg, f7;
    double alpa1, alpa2;
    double s11, s22, s33, s12, s13, s23;
    double w12, w13, w23;
    double xis, wis;
    double rho, u_vel, v_vel, w_vel;
    double zksav1, zksav2, bnum;
    double eta, squig, factre;
    double t11, t22_v, t33_v, t12_v, t13_v, t23_v;
    double xlnpt;
    double denom;
    double o11, o22, o33, o12_v, o13_v, o23_v;
    double t11p, t22p, t33p, t12p, t13p, t23p;
    double xis2013, s_mod, w_mod;
    double al10, al1, al2, al3, al4;
    double cmuu, eta1;
    double tau, fmu, zkolmog;
    double xmut;
    double c4new, eta1_girimaji;
    double tmp;
    double ux1, ux2, ux3, ux4, ux5, ux6, ux7, ux8, ux9;
    double zk_sa;
    double v3d;
    int    ip1, im2;
    int    l;


    // coefficient used to turn off eddy-viscosity based
    // turbulent stress calculation for ivmx==70 (full Reynolds stress model)
    coef_eddy = 1.0f;
    if (ivisc[0] >= 70) coef_eddy = 0.0f;
    if (i_tauijs == 1) coef_eddy = 0.0f;
    ccr2 = 2.5;

    idim1 = idim - 1;
    jdim1 = jdim - 1;

    // n  : number of cell centers for npl planes
    // l0 : number of cell interfaces for npl planes
    // jv : number of cell centers (and interfaces) on a i=constant plane
    n  = npl * jdim1 * idim1;
    l0 = npl * jdim1 * idim;
    jv = npl * jdim1;
    js = jv + 1;
    nn = n - jv;

    xmre  = (double)xmach / (double)reue;
    gpr   = (double)gamma / (double)pr;
    gm1pr = (double)gm1 * (double)pr;
    prtr  = (double)pr / (double)prt;
    // Durbin TCFD 1991 near-wall limiter (0=off)
    idurbinlim = 0;

    if (isklton > 0 && k == 1 && iadv >= 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "   computing viscous  fluxes, I-direction - central differencing");
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
        if (ivisc[0] == 11 || ivisc[0] == 12 || ivisc[0] == 13 ||
            ivisc[0] == 14 || i_nonlin != 0) {
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
    // t(1)  : 1/density
    // t(16) : c*c/(gm1*pr)
    // t(18) : u
    // t(19) : v
    // t(20) : w
    // t(22) : density
    l1 = jdim * idim - 1;
    for (kpl = 1; kpl <= npl; kpl++) {
        kk = k + kpl - 1;
        for (i = 1; i <= idim1; i++) {
            jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1 + 1;
            for (j = 1; j <= jdim1; j++) {
                t(j + jiv - 1, 22) = q(j, kk, i, 1);
                t(j + jiv - 1, 18) = q(j, kk, i, 2);
                t(j + jiv - 1, 19) = q(j, kk, i, 3);
                t(j + jiv - 1, 20) = q(j, kk, i, 4);
                t(j + jiv - 1, 16) = q(j, kk, i, 5);
            }
        }
    }

    // store dependent variables on i=0 and i=idim boundaries in wi0
    for (m = 0; m <= 10; m += 10) {
        if (m == 0) {
            mm = 1;
            mb = 1;
        } else {
            mm = 3;
            mb = 2;
        }
        jmin = 1;
        jmax = jdim - 1;
        kmin = k;
        kmax = k + npl - 1;
        for (kk = kmin; kk <= kmax; kk++) {
            jc = (kk - k) * jdim1 + jmin - 1;
            for (int jj = jmin; jj <= jmax; jj++) {
                jc = jc + 1;
                wi0(jc, m + 12) = qi0(jj, kk, 1, mm);
                wi0(jc, m + 8)  = qi0(jj, kk, 2, mm);
                wi0(jc, m + 9)  = qi0(jj, kk, 3, mm);
                wi0(jc, m + 10) = qi0(jj, kk, 4, mm);
                wi0(jc, m + 6)  = qi0(jj, kk, 5, mm);
                // m+5 flag indicates whether ghost cell or interface values stored in wi0
                wi0(jc, m + 5)  = bci(jj, kk, mb);
            }
        }
    }


    // t(17)/wi0(m+7)   : (u*u+v*v+w*w)/2
    // t(1)/wi0(m+1)    : 1/density
    // t(16)/wi0(m+6)   : c*c/(gm1*pr)
    // note: m=0 gives values at i=0, m=10 gives values at i=idim
    for (j = 1; j <= n; j++) {
        t(j, 17) = 0.5e0 * (t(j, 18) * t(j, 18)
                           + t(j, 19) * t(j, 19)
                           + t(j, 20) * t(j, 20));
        t(j, 21) = t(j, 17) * t(j, 22) + t(j, 16) / (double)gm1;
        t(j, 1)  = 1.e0 / t(j, 22);
        t(j, 16) = gpr * t(j, 16) * t(j, 1) / (double)gm1;
    }
    for (m = 0; m <= 10; m += 10) {
        for (j = 1; j <= jv; j++) {
            wi0(j, m + 7) = 0.5e0 * (wi0(j, m + 8)  * wi0(j, m + 8)
                                    + wi0(j, m + 9)  * wi0(j, m + 9)
                                    + wi0(j, m + 10) * wi0(j, m + 10));
            wi0(j, m + 11) = wi0(j, m + 7) * wi0(j, m + 12) + wi0(j, m + 6) / (double)gm1;
            wi0(j, m + 1)  = 1.e0 / wi0(j, m + 12);
            wi0(j, m + 6)  = gpr * wi0(j, m + 6) * wi0(j, m + 1) / (double)gm1;
        }
    }

    // t(7) : laminar viscosity at cell centers (via sutherland relation)
    c2b  = (double)cbar / (double)tinf;
    c2bp = c2b + 1.e0;
    for (j = 1; j <= n; j++) {
        t5     = gm1pr * t(j, 16);
        t6     = std::sqrt(t5);
        t(j, 7) = c2bp * t5 * t6 / (c2b + t5);
    }

    // t(14) : laminar viscosity values at cell interfaces
    // interior interfaces
    for (j = 1; j <= nn; j++) {
        t(j + js - 1, 14) = (t(j, 7) + t(j + js - 1, 7)) * 0.5e0;
    }
    // i=0 and i=idim interfaces
    for (izz = 1; izz <= jv; izz++) {
        ab  = 1.0 + wi0(izz, 5);
        bb  = 1.0 - wi0(izz, 5);
        wi05 = gm1pr * 0.5 * (ab * wi0(izz, 6) + bb * t(izz, 16));
        wi06 = std::sqrt(wi05);
        t(izz, 14)     = c2bp * wi05 * wi06 / (c2b + wi05);
        ab  = 1.0 + wi0(izz, 15);
        bb  = 1.0 - wi0(izz, 15);
        wi05 = gm1pr * 0.5 * (ab * wi0(izz, 16) + bb * t(izz + n - jv, 16));
        wi06 = std::sqrt(wi05);
        t(izz + n, 14) = c2bp * wi05 * wi06 / (c2b + wi05);
    }

    // t(15) : average jacobian (inverse volume) at cell interface
    n1 = jdim * idim1 + 1;
    l1 = jdim * idim - 1;
    for (kpl = 1; kpl <= npl; kpl++) {
        kk  = k + kpl - 1;
        for (i = 1; i <= idim; i++) {
            jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
            ji  = (i - 1) * jdim;
            if (i == 1) {
                // inverse volume at i=0 interface
                for (j = 1; j <= jdim1; j++) {
                    t(j + jiv, 15) = 2.0 / (voli0(j, kk, 1) + vol(j, kk, 1));
                }
            } else if (i == idim) {
                // inverse volume at i=idim interface
                for (j = 1; j <= jdim1; j++) {
                    t(j + jiv, 15) = 2.0 / (voli0(j, kk, 3) + vol(j, kk, idim1));
                }
            } else {
                // inverse volume at interior interfaces
                for (j = 1; j <= jdim1; j++) {
                    t(j + jiv, 15) = 2.0 / (vol(j, kk, i) + vol(j, kk, i - 1));
                }
            }

            // t(25) - t(27) : components of grad(xie)
            // t1    : grad(xie)/J
            // t(25) : d(xie)/dx
            // t(26) : d(xie)/dy
            // t(27) : d(xie)/dz
            for (j = 1; j <= jdim1; j++) {
                double t1_loc  = si(j, kk, i, 4) * t(j + jiv, 15);
                t(j + jiv, 25) = si(j, kk, i, 1) * t1_loc;
                t(j + jiv, 26) = si(j, kk, i, 2) * t1_loc;
                t(j + jiv, 27) = si(j, kk, i, 3) * t1_loc;
            }
        }
    }

    // t(6)-t(10) : gradients at cell interfaces
    // t(6)       : d( c*c/(gm1*pr) )/d(xie)
    // t(7)       : d( (u*u+v*v+w*w)/2 )/d(xie)
    // t(8)       : d(u)/d(xie)
    // t(9)       : d(v)/d(xie)
    // t(10)      : d(w)/d(xie)
    for (kq = 1; kq <= 5; kq++) {
        k1 = kq + 5;
        k2 = kq + 15;
        // interior interfaces
        for (j = 1; j <= nn; j++) {
            t(j + js - 1, k1) = t(j + js - 1, k2) - t(j, k2);
        }
        // interfaces at i=0/idim
        for (j = 1; j <= jv; j++) {
            double z1 = 1.0 + wi0(j, 5);
            double z2 = 1.0 + wi0(j, 15);
            t(j,     k1) = z1 * (t(j, k2) - wi0(j, k1));
            t(j + n, k1) = z2 * (wi0(j, k2) - t(j + n - jv, k2));
        }
    }


    // t(24) : turbulent viscosity at interfaces (=0 for laminar flow)
    iviscc = ivisc[0];
    if (iviscc > 1) {
        for (kpl = 1; kpl <= npl; kpl++) {
            kk = k + kpl - 1;
            for (i = 1; i <= idim; i++) {
                jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
                // i=0 interface
                if (i == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        t(j + jiv, 24) = bci(j, kk, 1) * vi0(j, kk, 1, 1) +
                            (1.0 - bci(j, kk, 1)) * 0.5 * (vi0(j, kk, 1, 1) + vist3d(j, kk, 1));
                    }
                // i=idim interface
                } else if (i == idim) {
                    for (j = 1; j <= jdim1; j++) {
                        t(j + jiv, 24) = bci(j, kk, 2) * vi0(j, kk, 1, 3) +
                            (1.0 - bci(j, kk, 2)) * 0.5 * (vi0(j, kk, 1, 3) + vist3d(j, kk, idim1));
                    }
                } else {
                    // interior interfaces
                    for (j = 1; j <= jdim1; j++) {
                        t(j + jiv, 24) = 0.5 * (vist3d(j, kk, i) + vist3d(j, kk, i - 1));
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

    for (j = 1; j <= l0; j++) {
        // ratio of turbulent to laminar viscosity
        t25 = t(j, 24) / t(j, 14);
        t24 = 1.e0 + t25 * (double)coef_eddy;
        // t(23) : ratio of laminar Prandtl number to total Prandtl number
        t(j, 23) = (1.e0 + prtr * t25) / t24;
        // t(14): (mach/re)*viscosity/J
        t(j, 14) = xmre * t24 * t(j, 14) / t(j, 15);
        // t(5) : t(14)*(grad(zeta))**2
        t(j, 5)  = t(j, 14) * (t(j, 25) * t(j, 25)
                              + t(j, 26) * t(j, 26)
                              + t(j, 27) * t(j, 27));
    }

    // t(11) : t(14)*((d(xie)/dx * d(u)/d(xie))
    //               + (d(xie)/dy * d(v)/d(xie))
    //               + (d(xie)/dz * d(w)/d(xie)))/3
    tr1 = 1.e0 / 3.e0;
    for (j = 1; j <= l0; j++) {
        t(j, 11) = tr1 * t(j, 14) * (t(j, 25) * t(j, 8)
                                    + t(j, 26) * t(j, 9)
                                    + t(j, 27) * t(j, 10));
    }

    // store off vmui....viscosity coefficient (laminar + turbulent) on i=0 and
    // i=idim boundaries for later use in force integration
    for (kpl = 1; kpl <= npl; kpl++) {
        kk = k + kpl - 1;
        for (j = 1; j <= jdim1; j++) {
            izz = (kpl - 1) * jdim1 + j;
            vmui(j, kk)             = t(izz,     14) * t(izz,     15) / xmre;
            vmui(j, kk + (kdim - 1)) = t(izz + n, 14) * t(izz + n, 15) / xmre;
        }
    }

    if (iadv < 0) return;

    // t(12) : contravariant velocity U at cell interfaces
    // interior interfaces
    for (j = 1; j <= nn; j++) {
        t(j + js - 1, 12) = 0.5e0 * (t(j + js - 1, 25) * (t(j + js - 1, 18) + t(j, 18))
                                    + t(j + js - 1, 26) * (t(j + js - 1, 19) + t(j, 19))
                                    + t(j + js - 1, 27) * (t(j + js - 1, 20) + t(j, 20)));
    }
    // interfaces at i=0 and i=idim
    for (izz = 1; izz <= jv; izz++) {
        ab = 1.0 + wi0(izz, 15);
        bb = 1.0 - wi0(izz, 15);
        t(izz + n, 12) = 0.5 * (t(izz + n, 25) * (ab * wi0(izz, 18) + bb * t(izz + n - jv, 18))
                               + t(izz + n, 26) * (ab * wi0(izz, 19) + bb * t(izz + n - jv, 19))
                               + t(izz + n, 27) * (ab * wi0(izz, 20) + bb * t(izz + n - jv, 20)));
        ab = 1.0 + wi0(izz, 5);
        bb = 1.0 - wi0(izz, 5);
        t(izz, 12) = 0.5 * (t(izz, 25) * (ab * wi0(izz, 8)  + bb * t(izz, 18))
                           + t(izz, 26) * (ab * wi0(izz, 9)  + bb * t(izz, 19))
                           + t(izz, 27) * (ab * wi0(izz, 10) + bb * t(izz, 20)));
    }

    // calculate fluxes
    for (l = 2; l <= 4; l++) {
        // viscous terms at interfaces (momentum eqs.)
        for (j = 1; j <= l0; j++) {
            t(j, 16) = t(j, 5) * t(j, l + 6) + t(j, l + 23) * t(j, 11);
        }
        // (-)viscous flux = Fv(k-1/2)-Fv(k+1/2) (momentum eqs.)
        for (j = 1; j <= n; j++) {
            t(j, l) = -t(j + js - 1, 16) + t(j, 16);
        }
    }

    // viscous terms at interfaces (energy eq.)
    for (j = 1; j <= l0; j++) {
        t(j, 16) = t(j, 5) * (t(j, 23) * t(j, 6) + t(j, 7))
                 + t(j, 12) * t(j, 11);
    }
    // (-)viscous flux = Fv(k-1/2)-Fv(k+1/2) (energy eq.)
    for (j = 1; j <= n; j++) {
        t(j, 5) = -t(j + js - 1, 16) + t(j, 16);
    }

    // calculate residuals
    for (kpl = 1; kpl <= npl; kpl++) {
        kk = k + kpl - 1;
        for (i = 1; i <= idim1; i++) {
            jiv = (i - 1) * jdim1 * npl + (kpl - 1) * jdim1 + 1;
            for (l = 2; l <= 5; l++) {
                for (j = 1; j <= jdim1; j++) {
                    res(j, kk, i, l) = res(j, kk, i, l) + t(j + jiv - 1, l);
                }
            }
        }
    }


    // Add tau_ij's for general method:
    if (i_tauijs == 1 && ivisc[0] < 70) {
        re = (double)reue / (double)xmach;
        for (kpl = 1; kpl <= npl; kpl++) {
            kk   = k + kpl - 1;
            jiv2 = (kpl - 1) * jdim1;
            for (i = 1; i <= idim; i++) {
                jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
                ip1 = std::min(i + 1, idim - 1);
                im2 = std::max(i - 2, 1);
                for (j = 1; j <= jdim1; j++) {
                    if (i == 1) {
                        rho   = 0.5 * (q(j, kk, i, 1) + qi0(j, kk, 1, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 1, 1) * wi0(j + jiv2, 5);
                        u_vel = 0.5 * (q(j, kk, i, 2) + qi0(j, kk, 2, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 2, 1) * wi0(j + jiv2, 5);
                        v_vel = 0.5 * (q(j, kk, i, 3) + qi0(j, kk, 3, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 3, 1) * wi0(j + jiv2, 5);
                        w_vel = 0.5 * (q(j, kk, i, 4) + qi0(j, kk, 4, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 4, 1) * wi0(j + jiv2, 5);
                        v3d   = 0.5 * (vist3d(j, kk, i) + vi0(j, kk, 1, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                vi0(j, kk, 1, 1) * wi0(j + jiv2, 5);
                        // use 1st order interpolation at end
                        ux1 = 1.5 * ux(j, kk, i, 1) - 0.5 * ux(j, kk, ip1, 1);
                        ux2 = 1.5 * ux(j, kk, i, 2) - 0.5 * ux(j, kk, ip1, 2);
                        ux3 = 1.5 * ux(j, kk, i, 3) - 0.5 * ux(j, kk, ip1, 3);
                        ux4 = 1.5 * ux(j, kk, i, 4) - 0.5 * ux(j, kk, ip1, 4);
                        ux5 = 1.5 * ux(j, kk, i, 5) - 0.5 * ux(j, kk, ip1, 5);
                        ux6 = 1.5 * ux(j, kk, i, 6) - 0.5 * ux(j, kk, ip1, 6);
                        ux7 = 1.5 * ux(j, kk, i, 7) - 0.5 * ux(j, kk, ip1, 7);
                        ux8 = 1.5 * ux(j, kk, i, 8) - 0.5 * ux(j, kk, ip1, 8);
                        ux9 = 1.5 * ux(j, kk, i, 9) - 0.5 * ux(j, kk, ip1, 9);
                    } else if (i == idim) {
                        rho   = 0.5 * (q(j, kk, i - 1, 1) + qi0(j, kk, 1, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 1, 3) * wi0(j + jiv2, 15);
                        u_vel = 0.5 * (q(j, kk, i - 1, 2) + qi0(j, kk, 2, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 2, 3) * wi0(j + jiv2, 15);
                        v_vel = 0.5 * (q(j, kk, i - 1, 3) + qi0(j, kk, 3, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 3, 3) * wi0(j + jiv2, 15);
                        w_vel = 0.5 * (q(j, kk, i - 1, 4) + qi0(j, kk, 4, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 4, 3) * wi0(j + jiv2, 15);
                        v3d   = 0.5 * (vist3d(j, kk, i - 1) + vi0(j, kk, 1, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                vi0(j, kk, 1, 3) * wi0(j + jiv2, 15);
                        // use 1st order interpolation at end
                        ux1 = 1.5 * ux(j, kk, i - 1, 1) - 0.5 * ux(j, kk, im2, 1);
                        ux2 = 1.5 * ux(j, kk, i - 1, 2) - 0.5 * ux(j, kk, im2, 2);
                        ux3 = 1.5 * ux(j, kk, i - 1, 3) - 0.5 * ux(j, kk, im2, 3);
                        ux4 = 1.5 * ux(j, kk, i - 1, 4) - 0.5 * ux(j, kk, im2, 4);
                        ux5 = 1.5 * ux(j, kk, i - 1, 5) - 0.5 * ux(j, kk, im2, 5);
                        ux6 = 1.5 * ux(j, kk, i - 1, 6) - 0.5 * ux(j, kk, im2, 6);
                        ux7 = 1.5 * ux(j, kk, i - 1, 7) - 0.5 * ux(j, kk, im2, 7);
                        ux8 = 1.5 * ux(j, kk, i - 1, 8) - 0.5 * ux(j, kk, im2, 8);
                        ux9 = 1.5 * ux(j, kk, i - 1, 9) - 0.5 * ux(j, kk, im2, 9);
                    } else {
                        rho   = 0.5 * (q(j, kk, i, 1) + q(j, kk, i - 1, 1));
                        u_vel = 0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2));
                        v_vel = 0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3));
                        w_vel = 0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4));
                        v3d   = 0.5 * (vist3d(j, kk, i) + vist3d(j, kk, i - 1));
                        ux1 = 0.5 * (ux(j, kk, i, 1) + ux(j, kk, i - 1, 1));
                        ux2 = 0.5 * (ux(j, kk, i, 2) + ux(j, kk, i - 1, 2));
                        ux3 = 0.5 * (ux(j, kk, i, 3) + ux(j, kk, i - 1, 3));
                        ux4 = 0.5 * (ux(j, kk, i, 4) + ux(j, kk, i - 1, 4));
                        ux5 = 0.5 * (ux(j, kk, i, 5) + ux(j, kk, i - 1, 5));
                        ux6 = 0.5 * (ux(j, kk, i, 6) + ux(j, kk, i - 1, 6));
                        ux7 = 0.5 * (ux(j, kk, i, 7) + ux(j, kk, i - 1, 7));
                        ux8 = 0.5 * (ux(j, kk, i, 8) + ux(j, kk, i - 1, 8));
                        ux9 = 0.5 * (ux(j, kk, i, 9) + ux(j, kk, i - 1, 9));
                    }
                    s11 = ux1;
                    s22 = ux5;
                    s33 = ux9;
                    double trace3 = (s11 + s22 + s33) / 3.0;
                    double s11star = s11 - trace3;
                    double s22star = s22 - trace3;
                    double s33star = s33 - trace3;
                    s12 = 0.5 * (ux2 + ux4);
                    s13 = 0.5 * (ux3 + ux7);
                    s23 = 0.5 * (ux6 + ux8);
                    w12 = 0.5 * (ux2 - ux4);
                    w13 = 0.5 * (ux3 - ux7);
                    w23 = 0.5 * (ux6 - ux8);
                    // ignoring -2/3*rho*k*delij term
                    t11   = -2.0 * v3d * s11star;
                    t22_v = -2.0 * v3d * s22star;
                    t33_v = -2.0 * v3d * s33star;
                    t12_v = -2.0 * v3d * s12;
                    t13_v = -2.0 * v3d * s13;
                    t23_v = -2.0 * v3d * s23;
                    if (i_qcr2000 == 1 || i_qcr2013 == 1 || i_qcr2013v == 1) {
                        denom = std::sqrt(ux1*ux1 + ux2*ux2 + ux3*ux3 +
                                          ux4*ux4 + ux5*ux5 + ux6*ux6 +
                                          ux7*ux7 + ux8*ux8 + ux9*ux9) + 1.e-20;
                        o11   = 0.0;
                        o22   = 0.0;
                        o33   = 0.0;
                        o12_v = 2.0 * w12 / denom;
                        o13_v = 2.0 * w13 / denom;
                        o23_v = 2.0 * w23 / denom;
                        t11p = t11   - 0.3 * (o11*t11   + o12_v*t12_v + o13_v*t13_v +
                                               o11*t11   + o12_v*t12_v + o13_v*t13_v);
                        t22p = t22_v - 0.3 * (-o12_v*t12_v + o22*t22_v + o23_v*t23_v +
                                               -o12_v*t12_v + o22*t22_v + o23_v*t23_v);
                        t33p = t33_v - 0.3 * (-o13_v*t13_v - o23_v*t23_v + o33*t33_v +
                                               -o13_v*t13_v - o23_v*t23_v + o33*t33_v);
                        t12p = t12_v - 0.3 * (o11*t12_v  + o12_v*t22_v + o13_v*t23_v +
                                               -o12_v*t11  + o22*t12_v  + o23_v*t13_v);
                        t13p = t13_v - 0.3 * (o11*t13_v  + o12_v*t23_v + o13_v*t33_v +
                                               -o13_v*t11  - o23_v*t12_v + o33*t13_v);
                        t23p = t23_v - 0.3 * (-o12_v*t13_v + o22*t23_v + o23_v*t33_v +
                                               -o13_v*t12_v - o23_v*t22_v + o33*t23_v);
                        t11   = t11p;
                        t22_v = t22p;
                        t33_v = t33p;
                        t12_v = t12p;
                        t13_v = t13p;
                        t23_v = t23p;
                    }
                    if (i_qcr2013 == 1) {
                        xis2013 = s11star*s11star + s22star*s22star + s33star*s33star +
                                  2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        s_mod = std::sqrt(2.0 * xis2013);
                        wis   = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        // limiting is typically necessary
                        w_mod = std::sqrt(2.0 * wis);
                        { double a1 = s_mod, a2 = 1.2 * w_mod; s_mod = ccomplex_ns::ccmin(a1, a2); }
                        t11   = t11   + ccr2 * v3d * s_mod;
                        t22_v = t22_v + ccr2 * v3d * s_mod;
                        t33_v = t33_v + ccr2 * v3d * s_mod;
                    }
                    if (i_qcr2013v == 1) {
                        wis   = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        w_mod = std::sqrt(2.0 * wis);
                        t11   = t11   + ccr2 * v3d * w_mod;
                        t22_v = t22_v + ccr2 * v3d * w_mod;
                        t33_v = t33_v + ccr2 * v3d * w_mod;
                    }
                    t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                             t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                    t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                             t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                    t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                             t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                    t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                        t(j + jiv, 25) * (u_vel * t11   + v_vel * t12_v + w_vel * t13_v) +
                        t(j + jiv, 26) * (u_vel * t12_v + v_vel * t22_v + w_vel * t23_v) +
                        t(j + jiv, 27) * (u_vel * t13_v + v_vel * t23_v + w_vel * t33_v));
                } // j loop
            } // i loop
        } // kpl loop
        // Take positive derivative
        for (l = 2; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                t(izz, l) = t(izz + js - 1, l + 26) - t(izz, l + 26);
            }
        }
        for (kpl = 1; kpl <= npl; kpl++) {
            kk = k + kpl - 1;
            for (i = 1; i <= idim1; i++) {
                jiv = (i - 1) * jdim1 * npl + (kpl - 1) * jdim1 + 1;
                for (l = 2; l <= 5; l++) {
                    for (izz = 1; izz <= jdim1; izz++) {
                        res(izz, kk, i, l) = res(izz, kk, i, l) + t(izz + jiv - 1, l);
                    }
                }
            }
        }
    } // end i_tauijs==1 branch


    // Add tau_ij's for RSMs:
    if (ivisc[0] >= 70) {
        re = (double)reue / (double)xmach;
        for (kpl = 1; kpl <= npl; kpl++) {
            kk   = k + kpl - 1;
            jiv2 = (kpl - 1) * jdim1;
            for (i = 1; i <= idim; i++) {
                jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
                if (i == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        rho   = 0.5 * (q(j, kk, i, 1) + qi0(j, kk, 1, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 1, 1) * wi0(j + jiv2, 5);
                        u_vel = 0.5 * (q(j, kk, i, 2) + qi0(j, kk, 2, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 2, 1) * wi0(j + jiv2, 5);
                        v_vel = 0.5 * (q(j, kk, i, 3) + qi0(j, kk, 3, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 3, 1) * wi0(j + jiv2, 5);
                        w_vel = 0.5 * (q(j, kk, i, 4) + qi0(j, kk, 4, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 4, 1) * wi0(j + jiv2, 5);
                        t11   = -rho * 0.5 * (zksav(j, kk, i, 1) + ti0(j, kk, 1, 1)) * re;
                        t22_v = -rho * 0.5 * (zksav(j, kk, i, 2) + ti0(j, kk, 2, 1)) * re;
                        t33_v = -rho * 0.5 * (zksav(j, kk, i, 3) + ti0(j, kk, 3, 1)) * re;
                        t12_v = -rho * 0.5 * (zksav(j, kk, i, 4) + ti0(j, kk, 4, 1)) * re;
                        t23_v = -rho * 0.5 * (zksav(j, kk, i, 5) + ti0(j, kk, 5, 1)) * re;
                        t13_v = -rho * 0.5 * (zksav(j, kk, i, 6) + ti0(j, kk, 6, 1)) * re;
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (u_vel * t11   + v_vel * t12_v + w_vel * t13_v) +
                            t(j + jiv, 26) * (u_vel * t12_v + v_vel * t22_v + w_vel * t23_v) +
                            t(j + jiv, 27) * (u_vel * t13_v + v_vel * t23_v + w_vel * t33_v));
                    }
                } else if (i == idim) {
                    for (j = 1; j <= jdim1; j++) {
                        rho   = 0.5 * (q(j, kk, i - 1, 1) + qi0(j, kk, 1, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 1, 3) * wi0(j + jiv2, 15);
                        u_vel = 0.5 * (q(j, kk, i - 1, 2) + qi0(j, kk, 2, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 2, 3) * wi0(j + jiv2, 15);
                        v_vel = 0.5 * (q(j, kk, i - 1, 3) + qi0(j, kk, 3, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 3, 3) * wi0(j + jiv2, 15);
                        w_vel = 0.5 * (q(j, kk, i - 1, 4) + qi0(j, kk, 4, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 4, 3) * wi0(j + jiv2, 15);
                        t11   = -rho * 0.5 * (zksav(j, kk, i - 1, 1) + ti0(j, kk, 1, 3)) * re;
                        t22_v = -rho * 0.5 * (zksav(j, kk, i - 1, 2) + ti0(j, kk, 2, 3)) * re;
                        t33_v = -rho * 0.5 * (zksav(j, kk, i - 1, 3) + ti0(j, kk, 3, 3)) * re;
                        t12_v = -rho * 0.5 * (zksav(j, kk, i - 1, 4) + ti0(j, kk, 4, 3)) * re;
                        t23_v = -rho * 0.5 * (zksav(j, kk, i - 1, 5) + ti0(j, kk, 5, 3)) * re;
                        t13_v = -rho * 0.5 * (zksav(j, kk, i - 1, 6) + ti0(j, kk, 6, 3)) * re;
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (u_vel * t11   + v_vel * t12_v + w_vel * t13_v) +
                            t(j + jiv, 26) * (u_vel * t12_v + v_vel * t22_v + w_vel * t23_v) +
                            t(j + jiv, 27) * (u_vel * t13_v + v_vel * t23_v + w_vel * t33_v));
                    }
                } else {
                    for (j = 1; j <= jdim1; j++) {
                        rho   = 0.5 * (q(j, kk, i, 1) + q(j, kk, i - 1, 1));
                        t11   = -rho * 0.5 * (zksav(j, kk, i, 1) + zksav(j, kk, i - 1, 1)) * re;
                        t22_v = -rho * 0.5 * (zksav(j, kk, i, 2) + zksav(j, kk, i - 1, 2)) * re;
                        t33_v = -rho * 0.5 * (zksav(j, kk, i, 3) + zksav(j, kk, i - 1, 3)) * re;
                        t12_v = -rho * 0.5 * (zksav(j, kk, i, 4) + zksav(j, kk, i - 1, 4)) * re;
                        t23_v = -rho * 0.5 * (zksav(j, kk, i, 5) + zksav(j, kk, i - 1, 5)) * re;
                        t13_v = -rho * 0.5 * (zksav(j, kk, i, 6) + zksav(j, kk, i - 1, 6)) * re;
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t11 +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t12_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t13_v) +
                            t(j + jiv, 26) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t12_v +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t22_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t23_v) +
                            t(j + jiv, 27) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t13_v +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t23_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t33_v));
                    }
                }
            } // i loop
        } // kpl loop
        // Take positive derivative
        for (l = 2; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                t(izz, l) = t(izz + js - 1, l + 26) - t(izz, l + 26);
            }
        }
        for (kpl = 1; kpl <= npl; kpl++) {
            kk = k + kpl - 1;
            for (i = 1; i <= idim1; i++) {
                jiv = (i - 1) * jdim1 * npl + (kpl - 1) * jdim1 + 1;
                for (l = 2; l <= 5; l++) {
                    for (izz = 1; izz <= jdim1; izz++) {
                        res(izz, kk, i, l) = res(izz, kk, i, l) + t(izz + jiv - 1, l);
                    }
                }
            }
        }
    } // end ivisc[0]>=70 branch


    // EASM v1: ivisc(1)==11 or 12
    if (i_tauijs != 1 && (ivisc[0] == 11 || ivisc[0] == 12) && level >= lglobal) {
        // f7=factor for determining whether it's a k-omega or k-epsilon formulation
        // f7=0 for k-epsilon, 1 for k-omega
        f7 = 0.0;
        if (ivisc[0] == 12) f7 = 1.0;
        re = (double)reue / (double)xmach;
        gg = 1.0 / ((double)c10 + (double)c5_e - 1.0);
        // Add nonlinear terms to RHS
        for (kpl = 1; kpl <= npl; kpl++) {
            kk   = k + kpl - 1;
            jiv2 = (kpl - 1) * jdim1;
            for (i = 1; i <= idim; i++) {
                jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
                if (i == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        // Determine Sij and Wij values (at cell interface)
                        s11 = t(j + jiv, 8)  * t(j + jiv, 25);
                        s22 = t(j + jiv, 9)  * t(j + jiv, 26);
                        s33 = t(j + jiv, 10) * t(j + jiv, 27);
                        s12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) + t(j + jiv, 9)  * t(j + jiv, 25));
                        s13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 25));
                        s23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 26));
                        w12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) - t(j + jiv, 9)  * t(j + jiv, 25));
                        w13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 25));
                        w23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 26));
                        xis = s11*s11 + s22*s22 + s33*s33 + 2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        // "Constants" are function of c3, c4, and gg:
                        alpa1 = (2.0 - (double)c4_e) / 2.0 * gg;
                        alpa2 = (2.0 - (double)c3_e) * gg;
                        // Find tauij values:
                        rho   = 0.5 * (q(j, kk, i, 1) + qi0(j, kk, 1, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 1, 1) * wi0(j + jiv2, 5);
                        u_vel = 0.5 * (q(j, kk, i, 2) + qi0(j, kk, 2, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 2, 1) * wi0(j + jiv2, 5);
                        v_vel = 0.5 * (q(j, kk, i, 3) + qi0(j, kk, 3, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 3, 1) * wi0(j + jiv2, 5);
                        w_vel = 0.5 * (q(j, kk, i, 4) + qi0(j, kk, 4, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 4, 1) * wi0(j + jiv2, 5);
                        zksav1 = 0.5 * (zksav(j, kk, i, 1) + ti0(j, kk, 1, 1));
                        { double a1 = zksav1, a2 = (double)tur10[0]; zksav1 = ccomplex_ns::ccmax(a1, a2); }
                        zksav2 = 0.5 * (zksav(j, kk, i, 2) + ti0(j, kk, 2, 1));
                        bnum   = zksav2 * (1.0 - f7) + f7;
                        eta    = (2.0 - (double)c3_e) * (2.0 - (double)c3_e) * (gg * gg / 4.0) * xis *
                                 (bnum / (zksav1 * re)) * (bnum / (zksav1 * re));
                        squig  = (2.0 - (double)c4_e) * (2.0 - (double)c4_e) * (gg * gg / 4.0) * wis *
                                 (bnum / (zksav1 * re)) * (bnum / (zksav1 * re));
                        { double a1 = eta,   a2 = 10.0; eta   = ccomplex_ns::ccmincr(a1, a2); }
                        { double a1 = squig, a2 = 10.0; squig = ccomplex_ns::ccmincr(a1, a2); }
                        factre = (3.0 * (1.0 + eta) + 0.2e-8 * (eta*eta*eta + squig*squig*squig)) /
                                 (3.0 * (1.0 + eta) + 0.2    * (eta*eta*eta + squig*squig*squig));
                        t11   = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 * (-s12*w12 - s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*s11 + s12*s12 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s11 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t11, a2 = xlnpt; t11 = ccomplex_ns::ccmax(a1, a2); }
                        t11   = t11 * (1.0 - bci(j, kk, 1));
                        t22_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 * (s12*w12 - s23*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s22*s22 + s12*s12 + s23*s23 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s22 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t22_v, a2 = xlnpt; t22_v = ccomplex_ns::ccmax(a1, a2); }
                        t22_v = t22_v * (1.0 - bci(j, kk, 1));
                        t33_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 * (s23*w23 + s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s33*s33 + s23*s23 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s33 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t33_v, a2 = xlnpt; t33_v = ccomplex_ns::ccmax(a1, a2); }
                        t33_v = t33_v * (1.0 - bci(j, kk, 1));
                        t12_v = -2.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*w12 - s22*w12 - s13*w23 - s23*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*s12 + s12*s22 + s13*s23) / re;
                        t12_v = t12_v * (1.0 - bci(j, kk, 1));
                        t13_v = -2.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*w13 + s12*w23 - s23*w12 - s33*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*s13 + s12*s23 + s13*s33) / re;
                        t13_v = t13_v * (1.0 - bci(j, kk, 1));
                        t23_v = -2.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s12*w13 + s22*w23 + s13*w12 - s33*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s12*s13 + s22*s23 + s23*s33) / re;
                        t23_v = t23_v * (1.0 - bci(j, kk, 1));
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (u_vel * t11   + v_vel * t12_v + w_vel * t13_v) +
                            t(j + jiv, 26) * (u_vel * t12_v + v_vel * t22_v + w_vel * t23_v) +
                            t(j + jiv, 27) * (u_vel * t13_v + v_vel * t23_v + w_vel * t33_v))
                          - xmre / t(j + jiv, 15) * (t(j + jiv, 14) * t(j + jiv, 15) / xmre
                            - t(j + jiv, 24) + (double)sigk1 * t(j + jiv, 24)) *
                            (t(j + jiv, 25) * t(j + jiv, 25) + t(j + jiv, 26) * t(j + jiv, 26) +
                             t(j + jiv, 27) * t(j + jiv, 27)) *
                            (zksav(j, kk, i, 2) - ti0(j, kk, 2, 1));
                    } // j loop i==1


                } else if (i == idim) {
                    for (j = 1; j <= jdim1; j++) {
                        // Determine Sij and Wij values (at cell interface)
                        s11 = t(j + jiv, 8)  * t(j + jiv, 25);
                        s22 = t(j + jiv, 9)  * t(j + jiv, 26);
                        s33 = t(j + jiv, 10) * t(j + jiv, 27);
                        s12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) + t(j + jiv, 9)  * t(j + jiv, 25));
                        s13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 25));
                        s23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 26));
                        w12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) - t(j + jiv, 9)  * t(j + jiv, 25));
                        w13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 25));
                        w23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 26));
                        xis = s11*s11 + s22*s22 + s33*s33 + 2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        alpa1 = (2.0 - (double)c4_e) / 2.0 * gg;
                        alpa2 = (2.0 - (double)c3_e) * gg;
                        rho   = 0.5 * (q(j, kk, i - 1, 1) + qi0(j, kk, 1, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 1, 3) * wi0(j + jiv2, 15);
                        u_vel = 0.5 * (q(j, kk, i - 1, 2) + qi0(j, kk, 2, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 2, 3) * wi0(j + jiv2, 15);
                        v_vel = 0.5 * (q(j, kk, i - 1, 3) + qi0(j, kk, 3, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 3, 3) * wi0(j + jiv2, 15);
                        w_vel = 0.5 * (q(j, kk, i - 1, 4) + qi0(j, kk, 4, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 4, 3) * wi0(j + jiv2, 15);
                        zksav1 = 0.5 * (zksav(j, kk, i - 1, 1) + ti0(j, kk, 1, 3));
                        { double a1 = zksav1, a2 = (double)tur10[0]; zksav1 = ccomplex_ns::ccmax(a1, a2); }
                        zksav2 = 0.5 * (zksav(j, kk, i - 1, 2) + ti0(j, kk, 2, 3));
                        bnum   = zksav2 * (1.0 - f7) + f7;
                        eta    = (2.0 - (double)c3_e) * (2.0 - (double)c3_e) * (gg * gg / 4.0) * xis *
                                 (bnum / (zksav1 * re)) * (bnum / (zksav1 * re));
                        squig  = (2.0 - (double)c4_e) * (2.0 - (double)c4_e) * (gg * gg / 4.0) * wis *
                                 (bnum / (zksav1 * re)) * (bnum / (zksav1 * re));
                        { double a1 = eta,   a2 = 10.0; eta   = ccomplex_ns::ccmincr(a1, a2); }
                        { double a1 = squig, a2 = 10.0; squig = ccomplex_ns::ccmincr(a1, a2); }
                        factre = (3.0 * (1.0 + eta) + 0.2e-8 * (eta*eta*eta + squig*squig*squig)) /
                                 (3.0 * (1.0 + eta) + 0.2    * (eta*eta*eta + squig*squig*squig));
                        t11   = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 * (-s12*w12 - s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*s11 + s12*s12 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s11 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t11, a2 = xlnpt; t11 = ccomplex_ns::ccmax(a1, a2); }
                        t11   = t11 * (1.0 - bci(j, kk, 2));
                        t22_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 * (s12*w12 - s23*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s22*s22 + s12*s12 + s23*s23 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s22 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t22_v, a2 = xlnpt; t22_v = ccomplex_ns::ccmax(a1, a2); }
                        t22_v = t22_v * (1.0 - bci(j, kk, 2));
                        t33_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 * (s23*w23 + s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s33*s33 + s23*s23 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s33 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t33_v, a2 = xlnpt; t33_v = ccomplex_ns::ccmax(a1, a2); }
                        t33_v = t33_v * (1.0 - bci(j, kk, 2));
                        t12_v = -2.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*w12 - s22*w12 - s13*w23 - s23*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*s12 + s12*s22 + s13*s23) / re;
                        t12_v = t12_v * (1.0 - bci(j, kk, 2));
                        t13_v = -2.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*w13 + s12*w23 - s23*w12 - s33*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*s13 + s12*s23 + s13*s33) / re;
                        t13_v = t13_v * (1.0 - bci(j, kk, 2));
                        t23_v = -2.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s12*w13 + s22*w23 + s13*w12 - s33*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s12*s13 + s22*s23 + s23*s33) / re;
                        t23_v = t23_v * (1.0 - bci(j, kk, 2));
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (u_vel * t11   + v_vel * t12_v + w_vel * t13_v) +
                            t(j + jiv, 26) * (u_vel * t12_v + v_vel * t22_v + w_vel * t23_v) +
                            t(j + jiv, 27) * (u_vel * t13_v + v_vel * t23_v + w_vel * t33_v))
                          - xmre / t(j + jiv, 15) * (t(j + jiv, 14) * t(j + jiv, 15) / xmre
                            - t(j + jiv, 24) + (double)sigk1 * t(j + jiv, 24)) *
                            (t(j + jiv, 25) * t(j + jiv, 25) + t(j + jiv, 26) * t(j + jiv, 26) +
                             t(j + jiv, 27) * t(j + jiv, 27)) *
                            (ti0(j, kk, 2, 3) - zksav(j, kk, i - 1, 2));
                    } // j loop i==idim


                } else {
                    // interior interfaces
                    for (j = 1; j <= jdim1; j++) {
                        // Determine Sij and Wij values (at cell interfaces):
                        s11 = t(j + jiv, 8)  * t(j + jiv, 25);
                        s22 = t(j + jiv, 9)  * t(j + jiv, 26);
                        s33 = t(j + jiv, 10) * t(j + jiv, 27);
                        s12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) + t(j + jiv, 9)  * t(j + jiv, 25));
                        s13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 25));
                        s23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 26));
                        w12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) - t(j + jiv, 9)  * t(j + jiv, 25));
                        w13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 25));
                        w23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 26));
                        xis = s11*s11 + s22*s22 + s33*s33 + 2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                        alpa1 = (2.0 - (double)c4_e) / 2.0 * gg;
                        alpa2 = (2.0 - (double)c3_e) * gg;
                        rho    = 0.5 * (q(j, kk, i, 1) + q(j, kk, i - 1, 1));
                        zksav1 = 0.5 * (zksav(j, kk, i, 1) + zksav(j, kk, i - 1, 1));
                        { double a1 = zksav1, a2 = (double)tur10[0]; zksav1 = ccomplex_ns::ccmax(a1, a2); }
                        zksav2 = 0.5 * (zksav(j, kk, i, 2) + zksav(j, kk, i - 1, 2));
                        bnum   = zksav2 * (1.0 - f7) + f7;
                        eta    = (2.0 - (double)c3_e) * (2.0 - (double)c3_e) * (gg * gg / 4.0) * xis *
                                 (bnum / (zksav1 * re)) * (bnum / (zksav1 * re));
                        squig  = (2.0 - (double)c4_e) * (2.0 - (double)c4_e) * (gg * gg / 4.0) * wis *
                                 (bnum / (zksav1 * re)) * (bnum / (zksav1 * re));
                        { double a1 = eta,   a2 = 10.0; eta   = ccomplex_ns::ccmincr(a1, a2); }
                        { double a1 = squig, a2 = 10.0; squig = ccomplex_ns::ccmincr(a1, a2); }
                        factre = (3.0 * (1.0 + eta) + 0.2e-8 * (eta*eta*eta + squig*squig*squig)) /
                                 (3.0 * (1.0 + eta) + 0.2    * (eta*eta*eta + squig*squig*squig));
                        t11   = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 * (-s12*w12 - s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*s11 + s12*s12 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s11 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t11, a2 = xlnpt; t11 = ccomplex_ns::ccmax(a1, a2); }
                        t22_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 * (s12*w12 - s23*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s22*s22 + s12*s12 + s23*s23 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s22 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t22_v, a2 = xlnpt; t22_v = ccomplex_ns::ccmax(a1, a2); }
                        t33_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 * (s23*w23 + s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s33*s33 + s23*s23 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s33 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t33_v, a2 = xlnpt; t33_v = ccomplex_ns::ccmax(a1, a2); }
                        t12_v = -2.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*w12 - s22*w12 - s13*w23 - s23*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*s12 + s12*s22 + s13*s23) / re;
                        t13_v = -2.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*w13 + s12*w23 - s23*w12 - s33*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s11*s13 + s12*s23 + s13*s33) / re;
                        t23_v = -2.0 * alpa1 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s12*w13 + s22*w23 + s13*w12 - s33*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * factre * bnum / zksav1 *
                                (s12*s13 + s22*s23 + s23*s33) / re;
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t11 +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t12_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t13_v) +
                            t(j + jiv, 26) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t12_v +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t22_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t23_v) +
                            t(j + jiv, 27) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t13_v +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t23_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t33_v))
                          - xmre / t(j + jiv, 15) * (t(j + jiv, 14) * t(j + jiv, 15) / xmre
                            - t(j + jiv, 24) + (double)sigk1 * t(j + jiv, 24)) *
                            (t(j + jiv, 25) * t(j + jiv, 25) + t(j + jiv, 26) * t(j + jiv, 26) +
                             t(j + jiv, 27) * t(j + jiv, 27)) *
                            (zksav(j, kk, i, 2) - zksav(j, kk, i - 1, 2));
                    } // j loop interior
                } // end i==1/idim/else
            } // i loop
        } // kpl loop
    } // end ivisc[0]==11||12 branch


    // EASM v2: ivisc(1)==13 or 14 (Girimaji fix)
    else if (i_tauijs != 1 && (ivisc[0] == 13 || ivisc[0] == 14) && level >= lglobal) {
        // f7=0 for k-epsilon, 1 for k-omega
        f7 = 0.0;
        if (ivisc[0] == 14) f7 = 1.0;
        re = (double)reue / (double)xmach;
        al10 = 0.5 * (double)c10 - 1.0;
        al1  = 2.0 * (0.5 * (double)c11 + 1.0);
        // 2-line improvement for wake
        if (ieasm_type == 0 || ieasm_type == 3 || ieasm_type == 4) {
            al10 = al10 + 1.8864;
            al1  = al1  - 2.0;
        }
        al2 = 0.5 * (double)c2_e - 2.0 / 3.0;
        al3 = 0.5 * (double)c3_e - 1.0;
        al4 = 0.5 * (double)c4_e - 1.0;
        // Add nonlinear terms to RHS
        for (kpl = 1; kpl <= npl; kpl++) {
            kk   = k + kpl - 1;
            jiv2 = (kpl - 1) * jdim1;
            for (i = 1; i <= idim; i++) {
                jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
                if (i == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        // Determine Sij and Wij values (at cell interface)
                        s11 = t(j + jiv, 8)  * t(j + jiv, 25);
                        s22 = t(j + jiv, 9)  * t(j + jiv, 26);
                        s33 = t(j + jiv, 10) * t(j + jiv, 27);
                        s12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) + t(j + jiv, 9)  * t(j + jiv, 25));
                        s13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 25));
                        s23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 26));
                        w12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) - t(j + jiv, 9)  * t(j + jiv, 25));
                        w13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 25));
                        w23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 26));
                        xis = s11*s11 + s22*s22 + s33*s33 + 2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        // Girimaji JFM 2000 fix to c4
                        if (ieasm_type == 4) {
                            wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                            eta1_girimaji = xis / (xis + wis);
                            if ((float)eta1_girimaji < 0.5f) {
                                c4new = 2.0 - ((2.0 - (double)c4_e) *
                                        std::pow(eta1_girimaji / (1.0 - eta1_girimaji), 0.75));
                            } else {
                                c4new = (double)c4_e;
                            }
                            al4 = 0.5 * c4new - 1.0;
                        }
                        rho   = 0.5 * (q(j, kk, i, 1) + qi0(j, kk, 1, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 1, 1) * wi0(j + jiv2, 5);
                        u_vel = 0.5 * (q(j, kk, i, 2) + qi0(j, kk, 2, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 2, 1) * wi0(j + jiv2, 5);
                        v_vel = 0.5 * (q(j, kk, i, 3) + qi0(j, kk, 3, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 3, 1) * wi0(j + jiv2, 5);
                        w_vel = 0.5 * (q(j, kk, i, 4) + qi0(j, kk, 4, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 4, 1) * wi0(j + jiv2, 5);
                        zksav1 = 0.5 * (zksav(j, kk, i, 1) + ti0(j, kk, 1, 1));
                        zksav2 = 0.5 * (zksav(j, kk, i, 2) + ti0(j, kk, 2, 1));
                        cmuu   = cmuv(j, kk, i);
                        bnum   = zksav2 * (1.0 - f7) + f7;
                        // Durbin TCFD 1991 near-wall limiter
                        if (idurbinlim != 0 && (ieasm_type == 3 || ieasm_type == 4)) {
                            tau     = bnum / zksav1;
                            fmu     = t(j + jiv, 14) * t(j + jiv, 15) / xmre - t(j + jiv, 24);
                            zkolmog = 6.0 * std::sqrt(fmu * bnum / (rho * zksav2 * zksav1));
                            { double a1 = tau, a2 = zkolmog; tau = ccomplex_ns::ccmax(a1, a2); }
                            zksav1  = bnum / tau;
                        }
                        { double a1 = zksav1, a2 = (double)tur10[0]; zksav1 = ccomplex_ns::ccmax(a1, a2); }
                        eta1  = xis * (bnum / zksav1 / re) * (bnum / zksav1 / re);
                        alpa1 = -al4 / (al10 - eta1 * al1 * cmuu);
                        alpa2 = -2.0 * al3 / (al10 - eta1 * cmuu * al1);
                        t11   = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 * (-s12*w12 - s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*s11 + s12*s12 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s11 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t11, a2 = xlnpt; t11 = ccomplex_ns::ccmax(a1, a2); }
                        t11   = t11 * (1.0 - bci(j, kk, 1));
                        t22_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 * (s12*w12 - s23*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s22*s22 + s12*s12 + s23*s23 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s22 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t22_v, a2 = xlnpt; t22_v = ccomplex_ns::ccmax(a1, a2); }
                        t22_v = t22_v * (1.0 - bci(j, kk, 1));
                        t33_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 * (s23*w23 + s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s33*s33 + s23*s23 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s33 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t33_v, a2 = xlnpt; t33_v = ccomplex_ns::ccmax(a1, a2); }
                        t33_v = t33_v * (1.0 - bci(j, kk, 1));
                        t12_v = -2.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*w12 - s22*w12 - s13*w23 - s23*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*s12 + s12*s22 + s13*s23) / re;
                        t12_v = t12_v * (1.0 - bci(j, kk, 1));
                        t13_v = -2.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*w13 + s12*w23 - s23*w12 - s33*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*s13 + s12*s23 + s13*s33) / re;
                        t13_v = t13_v * (1.0 - bci(j, kk, 1));
                        t23_v = -2.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 *
                                (s12*w13 + s22*w23 + s13*w12 - s33*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s12*s13 + s22*s23 + s23*s33) / re;
                        t23_v = t23_v * (1.0 - bci(j, kk, 1));
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        xmut = f7 * t(j + jiv, 24) + (1.0 - f7) * (double)cmuc1 * rho * zksav2 * zksav2 / zksav1;
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (u_vel * t11   + v_vel * t12_v + w_vel * t13_v) +
                            t(j + jiv, 26) * (u_vel * t12_v + v_vel * t22_v + w_vel * t23_v) +
                            t(j + jiv, 27) * (u_vel * t13_v + v_vel * t23_v + w_vel * t33_v))
                          - xmre / t(j + jiv, 15) * (t(j + jiv, 14) * t(j + jiv, 15) / xmre
                            - t(j + jiv, 24) + (double)sigk1 * xmut) *
                            (t(j + jiv, 25) * t(j + jiv, 25) + t(j + jiv, 26) * t(j + jiv, 26) +
                             t(j + jiv, 27) * t(j + jiv, 27)) *
                            (zksav(j, kk, i, 2) - ti0(j, kk, 2, 1));
                    } // j loop i==1


                } else if (i == idim) {
                    for (j = 1; j <= jdim1; j++) {
                        // Determine Sij and Wij values (at cell interface)
                        s11 = t(j + jiv, 8)  * t(j + jiv, 25);
                        s22 = t(j + jiv, 9)  * t(j + jiv, 26);
                        s33 = t(j + jiv, 10) * t(j + jiv, 27);
                        s12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) + t(j + jiv, 9)  * t(j + jiv, 25));
                        s13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 25));
                        s23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 26));
                        w12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) - t(j + jiv, 9)  * t(j + jiv, 25));
                        w13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 25));
                        w23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 26));
                        xis = s11*s11 + s22*s22 + s33*s33 + 2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        // Girimaji JFM 2000 fix to c4
                        if (ieasm_type == 4) {
                            wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                            eta1_girimaji = xis / (xis + wis);
                            if ((float)eta1_girimaji < 0.5f) {
                                c4new = 2.0 - ((2.0 - (double)c4_e) *
                                        std::pow(eta1_girimaji / (1.0 - eta1_girimaji), 0.75));
                            } else {
                                c4new = (double)c4_e;
                            }
                            al4 = 0.5 * c4new - 1.0;
                        }
                        rho   = 0.5 * (q(j, kk, i - 1, 1) + qi0(j, kk, 1, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 1, 3) * wi0(j + jiv2, 15);
                        u_vel = 0.5 * (q(j, kk, i - 1, 2) + qi0(j, kk, 2, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 2, 3) * wi0(j + jiv2, 15);
                        v_vel = 0.5 * (q(j, kk, i - 1, 3) + qi0(j, kk, 3, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 3, 3) * wi0(j + jiv2, 15);
                        w_vel = 0.5 * (q(j, kk, i - 1, 4) + qi0(j, kk, 4, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 4, 3) * wi0(j + jiv2, 15);
                        zksav1 = 0.5 * (zksav(j, kk, i - 1, 1) + ti0(j, kk, 1, 3));
                        zksav2 = 0.5 * (zksav(j, kk, i - 1, 2) + ti0(j, kk, 2, 3));
                        cmuu   = cmuv(j, kk, i - 1);
                        bnum   = zksav2 * (1.0 - f7) + f7;
                        // Durbin TCFD 1991 near-wall limiter
                        if (idurbinlim != 0 && (ieasm_type == 3 || ieasm_type == 4)) {
                            tau     = bnum / zksav1;
                            fmu     = t(j + jiv, 14) * t(j + jiv, 15) / xmre - t(j + jiv, 24);
                            zkolmog = 6.0 * std::sqrt(fmu * bnum / (rho * zksav2 * zksav1));
                            { double a1 = tau, a2 = zkolmog; tau = ccomplex_ns::ccmax(a1, a2); }
                            zksav1  = bnum / tau;
                        }
                        { double a1 = zksav1, a2 = (double)tur10[0]; zksav1 = ccomplex_ns::ccmax(a1, a2); }
                        eta1  = xis * (bnum / zksav1 / re) * (bnum / zksav1 / re);
                        alpa1 = -al4 / (al10 - eta1 * al1 * cmuu);
                        alpa2 = -2.0 * al3 / (al10 - eta1 * cmuu * al1);
                        t11   = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 * (-s12*w12 - s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*s11 + s12*s12 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s11 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t11, a2 = xlnpt; t11 = ccomplex_ns::ccmax(a1, a2); }
                        t11   = t11 * (1.0 - bci(j, kk, 2));
                        t22_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 * (s12*w12 - s23*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s22*s22 + s12*s12 + s23*s23 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s22 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t22_v, a2 = xlnpt; t22_v = ccomplex_ns::ccmax(a1, a2); }
                        t22_v = t22_v * (1.0 - bci(j, kk, 2));
                        t33_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 * (s23*w23 + s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s33*s33 + s23*s23 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s33 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t33_v, a2 = xlnpt; t33_v = ccomplex_ns::ccmax(a1, a2); }
                        t33_v = t33_v * (1.0 - bci(j, kk, 2));
                        t12_v = -2.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*w12 - s22*w12 - s13*w23 - s23*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*s12 + s12*s22 + s13*s23) / re;
                        t12_v = t12_v * (1.0 - bci(j, kk, 2));
                        t13_v = -2.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*w13 + s12*w23 - s23*w12 - s33*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*s13 + s12*s23 + s13*s33) / re;
                        t13_v = t13_v * (1.0 - bci(j, kk, 2));
                        t23_v = -2.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 *
                                (s12*w13 + s22*w23 + s13*w12 - s33*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s12*s13 + s22*s23 + s23*s33) / re;
                        t23_v = t23_v * (1.0 - bci(j, kk, 2));
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        xmut = f7 * t(j + jiv, 24) + (1.0 - f7) * (double)cmuc1 * rho * zksav2 * zksav2 / zksav1;
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (u_vel * t11   + v_vel * t12_v + w_vel * t13_v) +
                            t(j + jiv, 26) * (u_vel * t12_v + v_vel * t22_v + w_vel * t23_v) +
                            t(j + jiv, 27) * (u_vel * t13_v + v_vel * t23_v + w_vel * t33_v))
                          - xmre / t(j + jiv, 15) * (t(j + jiv, 14) * t(j + jiv, 15) / xmre
                            - t(j + jiv, 24) + (double)sigk1 * xmut) *
                            (t(j + jiv, 25) * t(j + jiv, 25) + t(j + jiv, 26) * t(j + jiv, 26) +
                             t(j + jiv, 27) * t(j + jiv, 27)) *
                            (ti0(j, kk, 2, 3) - zksav(j, kk, i - 1, 2));
                    } // j loop i==idim


                } else {
                    // interior interfaces
                    for (j = 1; j <= jdim1; j++) {
                        // Determine Sij and Wij values (at cell interfaces):
                        s11 = t(j + jiv, 8)  * t(j + jiv, 25);
                        s22 = t(j + jiv, 9)  * t(j + jiv, 26);
                        s33 = t(j + jiv, 10) * t(j + jiv, 27);
                        s12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) + t(j + jiv, 9)  * t(j + jiv, 25));
                        s13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 25));
                        s23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 26));
                        w12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) - t(j + jiv, 9)  * t(j + jiv, 25));
                        w13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 25));
                        w23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 26));
                        xis = s11*s11 + s22*s22 + s33*s33 + 2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        // Girimaji JFM 2000 fix to c4
                        if (ieasm_type == 4) {
                            wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                            eta1_girimaji = xis / (xis + wis);
                            if ((float)eta1_girimaji < 0.5f) {
                                c4new = 2.0 - ((2.0 - (double)c4_e) *
                                        std::pow(eta1_girimaji / (1.0 - eta1_girimaji), 0.75));
                            } else {
                                c4new = (double)c4_e;
                            }
                            al4 = 0.5 * c4new - 1.0;
                        }
                        rho    = 0.5 * (q(j, kk, i, 1) + q(j, kk, i - 1, 1));
                        zksav1 = 0.5 * (zksav(j, kk, i, 1) + zksav(j, kk, i - 1, 1));
                        zksav2 = 0.5 * (zksav(j, kk, i, 2) + zksav(j, kk, i - 1, 2));
                        cmuu   = 0.5 * (cmuv(j, kk, i) + cmuv(j, kk, i - 1));
                        bnum   = zksav2 * (1.0 - f7) + f7;
                        // Durbin TCFD 1991 near-wall limiter
                        if (idurbinlim != 0 && (ieasm_type == 3 || ieasm_type == 4)) {
                            tau     = bnum / zksav1;
                            fmu     = t(j + jiv, 14) * t(j + jiv, 15) / xmre - t(j + jiv, 24);
                            zkolmog = 6.0 * std::sqrt(fmu * bnum / (rho * zksav2 * zksav1));
                            { double a1 = tau, a2 = zkolmog; tau = ccomplex_ns::ccmax(a1, a2); }
                            zksav1  = bnum / tau;
                        }
                        { double a1 = zksav1, a2 = (double)tur10[0]; zksav1 = ccomplex_ns::ccmax(a1, a2); }
                        eta1  = xis * (bnum / zksav1 / re) * (bnum / zksav1 / re);
                        alpa1 = -al4 / (al10 - eta1 * al1 * cmuu);
                        alpa2 = -2.0 * al3 / (al10 - eta1 * cmuu * al1);
                        t11   = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 * (-s12*w12 - s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*s11 + s12*s12 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s11 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t11, a2 = xlnpt; t11 = ccomplex_ns::ccmax(a1, a2); }
                        t22_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 * (s12*w12 - s23*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s22*s22 + s12*s12 + s23*s23 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s22 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t22_v, a2 = xlnpt; t22_v = ccomplex_ns::ccmax(a1, a2); }
                        t33_v = 2.0 * rho * zksav2 * re / 3.0
                              - 4.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 * (s23*w23 + s13*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s33*s33 + s23*s23 + s13*s13 - 0.33333 * xis) / re;
                        xlnpt = 2.0 * t(j + jiv, 24) * (s33 - (s11 + s22 + s33) / 3.0);
                        { double a1 = t33_v, a2 = xlnpt; t33_v = ccomplex_ns::ccmax(a1, a2); }
                        t12_v = -2.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*w12 - s22*w12 - s13*w23 - s23*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*s12 + s12*s22 + s13*s23) / re;
                        t13_v = -2.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*w13 + s12*w23 - s23*w12 - s33*w13) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s11*s13 + s12*s23 + s13*s33) / re;
                        t23_v = -2.0 * alpa1 * t(j + jiv, 24) * bnum / zksav1 *
                                (s12*w13 + s22*w23 + s13*w12 - s33*w23) / re
                              + 2.0 * alpa2 * t(j + jiv, 24) * bnum / zksav1 *
                                (s12*s13 + s22*s23 + s23*s33) / re;
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        xmut = f7 * t(j + jiv, 24) + (1.0 - f7) * (double)cmuc1 * rho * zksav2 * zksav2 / zksav1;
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t11 +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t12_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t13_v) +
                            t(j + jiv, 26) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t12_v +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t22_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t23_v) +
                            t(j + jiv, 27) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t13_v +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t23_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t33_v))
                          - xmre / t(j + jiv, 15) * (t(j + jiv, 14) * t(j + jiv, 15) / xmre
                            - t(j + jiv, 24) + (double)sigk1 * xmut) *
                            (t(j + jiv, 25) * t(j + jiv, 25) + t(j + jiv, 26) * t(j + jiv, 26) +
                             t(j + jiv, 27) * t(j + jiv, 27)) *
                            (zksav(j, kk, i, 2) - zksav(j, kk, i - 1, 2));
                    } // j loop interior
                } // end i==1/idim/else
            } // i loop
        } // kpl loop
    } // end ivisc[0]==13||14 branch


    // Generic nonlinear constitutive: i_nonlin != 0
    else if (i_tauijs != 1 && i_nonlin != 0 && level >= lglobal) {
        re = (double)reue / (double)xmach;
        // Add nonlinear terms to RHS
        for (kpl = 1; kpl <= npl; kpl++) {
            kk   = k + kpl - 1;
            jiv2 = (kpl - 1) * jdim1;
            for (i = 1; i <= idim; i++) {
                jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
                if (i == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        // Determine Sij and Wij values (at cell interface)
                        s11 = t(j + jiv, 8)  * t(j + jiv, 25);
                        s22 = t(j + jiv, 9)  * t(j + jiv, 26);
                        s33 = t(j + jiv, 10) * t(j + jiv, 27);
                        s12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) + t(j + jiv, 9)  * t(j + jiv, 25));
                        s13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 25));
                        s23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 26));
                        w12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) - t(j + jiv, 9)  * t(j + jiv, 25));
                        w13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 25));
                        w23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 26));
                        tmp = (s11 + s22 + s33) / 3.0;
                        s11 = s11 - tmp;
                        s22 = s22 - tmp;
                        s33 = s33 - tmp;
                        // Note: denom term is taken at nearest cell center
                        denom = ux(j, kk, i, 1) * ux(j, kk, i, 1) + ux(j, kk, i, 2) * ux(j, kk, i, 2) +
                                ux(j, kk, i, 3) * ux(j, kk, i, 3) + ux(j, kk, i, 4) * ux(j, kk, i, 4) +
                                ux(j, kk, i, 5) * ux(j, kk, i, 5) + ux(j, kk, i, 6) * ux(j, kk, i, 6) +
                                ux(j, kk, i, 7) * ux(j, kk, i, 7) + ux(j, kk, i, 8) * ux(j, kk, i, 8) +
                                ux(j, kk, i, 9) * ux(j, kk, i, 9);
                        denom = std::sqrt(denom);
                        { double a1 = denom, a2 = (double)snonlin_lim; denom = ccomplex_ns::ccmax(a1, a2); }
                        // Find tauij values:
                        rho   = 0.5 * (q(j, kk, i, 1) + qi0(j, kk, 1, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 1, 1) * wi0(j + jiv2, 5);
                        u_vel = 0.5 * (q(j, kk, i, 2) + qi0(j, kk, 2, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 2, 1) * wi0(j + jiv2, 5);
                        v_vel = 0.5 * (q(j, kk, i, 3) + qi0(j, kk, 3, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 3, 1) * wi0(j + jiv2, 5);
                        w_vel = 0.5 * (q(j, kk, i, 4) + qi0(j, kk, 4, 1)) * (1.0 - wi0(j + jiv2, 5)) +
                                qi0(j, kk, 4, 1) * wi0(j + jiv2, 5);
                        zk_sa = 0.0;
                        if (ivisc[0] == 6 || ivisc[0] == 7 || ivisc[0] == 10 || ivisc[0] == 15) {
                            zk_sa = 0.5 * (zksav(j, kk, i, 2) + ti0(j, kk, 2, 1));
                        }
                        t11   = 2.0 * rho * zk_sa * re / 3.0
                              - 8.0 * (double)c_nonlin * t(j + jiv, 24) / denom * (-s12*w12 - s13*w13);
                        t11   = t11 * (1.0 - bci(j, kk, 1));
                        t22_v = 2.0 * rho * zk_sa * re / 3.0
                              - 8.0 * (double)c_nonlin * t(j + jiv, 24) / denom * (s12*w12 - s23*w23);
                        t22_v = t22_v * (1.0 - bci(j, kk, 1));
                        t33_v = 2.0 * rho * zk_sa * re / 3.0
                              - 8.0 * (double)c_nonlin * t(j + jiv, 24) / denom * (s23*w23 + s13*w13);
                        t33_v = t33_v * (1.0 - bci(j, kk, 1));
                        t12_v = -4.0 * (double)c_nonlin * t(j + jiv, 24) / denom *
                                (s11*w12 - s22*w12 - s13*w23 - s23*w13);
                        t12_v = t12_v * (1.0 - bci(j, kk, 1));
                        t13_v = -4.0 * (double)c_nonlin * t(j + jiv, 24) / denom *
                                (s11*w13 + s12*w23 - s23*w12 - s33*w13);
                        t13_v = t13_v * (1.0 - bci(j, kk, 1));
                        t23_v = -4.0 * (double)c_nonlin * t(j + jiv, 24) / denom *
                                (s12*w13 + s22*w23 + s13*w12 - s33*w23);
                        t23_v = t23_v * (1.0 - bci(j, kk, 1));
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (u_vel * t11   + v_vel * t12_v + w_vel * t13_v) +
                            t(j + jiv, 26) * (u_vel * t12_v + v_vel * t22_v + w_vel * t23_v) +
                            t(j + jiv, 27) * (u_vel * t13_v + v_vel * t23_v + w_vel * t33_v));
                    } // j loop i==1


                } else if (i == idim) {
                    for (j = 1; j <= jdim1; j++) {
                        // Determine Sij and Wij values (at cell interface)
                        s11 = t(j + jiv, 8)  * t(j + jiv, 25);
                        s22 = t(j + jiv, 9)  * t(j + jiv, 26);
                        s33 = t(j + jiv, 10) * t(j + jiv, 27);
                        s12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) + t(j + jiv, 9)  * t(j + jiv, 25));
                        s13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 25));
                        s23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 26));
                        w12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) - t(j + jiv, 9)  * t(j + jiv, 25));
                        w13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 25));
                        w23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 26));
                        tmp = (s11 + s22 + s33) / 3.0;
                        s11 = s11 - tmp;
                        s22 = s22 - tmp;
                        s33 = s33 - tmp;
                        // Note: denom term is taken at nearest cell center
                        denom = ux(j, kk, i - 1, 1) * ux(j, kk, i - 1, 1) +
                                ux(j, kk, i - 1, 2) * ux(j, kk, i - 1, 2) +
                                ux(j, kk, i - 1, 3) * ux(j, kk, i - 1, 3) +
                                ux(j, kk, i - 1, 4) * ux(j, kk, i - 1, 4) +
                                ux(j, kk, i - 1, 5) * ux(j, kk, i - 1, 5) +
                                ux(j, kk, i - 1, 6) * ux(j, kk, i - 1, 6) +
                                ux(j, kk, i - 1, 7) * ux(j, kk, i - 1, 7) +
                                ux(j, kk, i - 1, 8) * ux(j, kk, i - 1, 8) +
                                ux(j, kk, i - 1, 9) * ux(j, kk, i - 1, 9);
                        denom = std::sqrt(denom);
                        { double a1 = denom, a2 = (double)snonlin_lim; denom = ccomplex_ns::ccmax(a1, a2); }
                        // Find tauij values:
                        rho   = 0.5 * (q(j, kk, i - 1, 1) + qi0(j, kk, 1, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 1, 3) * wi0(j + jiv2, 15);
                        u_vel = 0.5 * (q(j, kk, i - 1, 2) + qi0(j, kk, 2, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 2, 3) * wi0(j + jiv2, 15);
                        v_vel = 0.5 * (q(j, kk, i - 1, 3) + qi0(j, kk, 3, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 3, 3) * wi0(j + jiv2, 15);
                        w_vel = 0.5 * (q(j, kk, i - 1, 4) + qi0(j, kk, 4, 3)) * (1.0 - wi0(j + jiv2, 15)) +
                                qi0(j, kk, 4, 3) * wi0(j + jiv2, 15);
                        zk_sa = 0.0;
                        if (ivisc[0] == 6 || ivisc[0] == 7 || ivisc[0] == 10 || ivisc[0] == 15) {
                            zk_sa = 0.5 * (zksav(j, kk, i - 1, 2) + ti0(j, kk, 2, 3));
                        }
                        t11   = 2.0 * rho * zk_sa * re / 3.0
                              - 8.0 * (double)c_nonlin * t(j + jiv, 24) / denom * (-s12*w12 - s13*w13);
                        t11   = t11 * (1.0 - bci(j, kk, 2));
                        t22_v = 2.0 * rho * zk_sa * re / 3.0
                              - 8.0 * (double)c_nonlin * t(j + jiv, 24) / denom * (s12*w12 - s23*w23);
                        t22_v = t22_v * (1.0 - bci(j, kk, 2));
                        t33_v = 2.0 * rho * zk_sa * re / 3.0
                              - 8.0 * (double)c_nonlin * t(j + jiv, 24) / denom * (s23*w23 + s13*w13);
                        t33_v = t33_v * (1.0 - bci(j, kk, 2));
                        t12_v = -4.0 * (double)c_nonlin * t(j + jiv, 24) / denom *
                                (s11*w12 - s22*w12 - s13*w23 - s23*w13);
                        t12_v = t12_v * (1.0 - bci(j, kk, 2));
                        t13_v = -4.0 * (double)c_nonlin * t(j + jiv, 24) / denom *
                                (s11*w13 + s12*w23 - s23*w12 - s33*w13);
                        t13_v = t13_v * (1.0 - bci(j, kk, 2));
                        t23_v = -4.0 * (double)c_nonlin * t(j + jiv, 24) / denom *
                                (s12*w13 + s22*w23 + s13*w12 - s33*w23);
                        t23_v = t23_v * (1.0 - bci(j, kk, 2));
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (u_vel * t11   + v_vel * t12_v + w_vel * t13_v) +
                            t(j + jiv, 26) * (u_vel * t12_v + v_vel * t22_v + w_vel * t23_v) +
                            t(j + jiv, 27) * (u_vel * t13_v + v_vel * t23_v + w_vel * t33_v));
                    } // j loop i==idim


                } else {
                    // interior interfaces
                    for (j = 1; j <= jdim1; j++) {
                        // Determine Sij and Wij values (at cell interfaces):
                        s11 = t(j + jiv, 8)  * t(j + jiv, 25);
                        s22 = t(j + jiv, 9)  * t(j + jiv, 26);
                        s33 = t(j + jiv, 10) * t(j + jiv, 27);
                        s12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) + t(j + jiv, 9)  * t(j + jiv, 25));
                        s13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 25));
                        s23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) + t(j + jiv, 10) * t(j + jiv, 26));
                        w12 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 26) - t(j + jiv, 9)  * t(j + jiv, 25));
                        w13 = 0.5 * (t(j + jiv, 8)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 25));
                        w23 = 0.5 * (t(j + jiv, 9)  * t(j + jiv, 27) - t(j + jiv, 10) * t(j + jiv, 26));
                        tmp = (s11 + s22 + s33) / 3.0;
                        s11 = s11 - tmp;
                        s22 = s22 - tmp;
                        s33 = s33 - tmp;
                        // Note: denom term is averaged to cell face
                        ux1 = 0.5 * (ux(j, kk, i - 1, 1) + ux(j, kk, i, 1));
                        ux2 = 0.5 * (ux(j, kk, i - 1, 2) + ux(j, kk, i, 2));
                        ux3 = 0.5 * (ux(j, kk, i - 1, 3) + ux(j, kk, i, 3));
                        ux4 = 0.5 * (ux(j, kk, i - 1, 4) + ux(j, kk, i, 4));
                        ux5 = 0.5 * (ux(j, kk, i - 1, 5) + ux(j, kk, i, 5));
                        ux6 = 0.5 * (ux(j, kk, i - 1, 6) + ux(j, kk, i, 6));
                        ux7 = 0.5 * (ux(j, kk, i - 1, 7) + ux(j, kk, i, 7));
                        ux8 = 0.5 * (ux(j, kk, i - 1, 8) + ux(j, kk, i, 8));
                        ux9 = 0.5 * (ux(j, kk, i - 1, 9) + ux(j, kk, i, 9));
                        denom = ux1*ux1 + ux2*ux2 + ux3*ux3 + ux4*ux4 + ux5*ux5 +
                                ux6*ux6 + ux7*ux7 + ux8*ux8 + ux9*ux9;
                        denom = std::sqrt(denom);
                        { double a1 = denom, a2 = (double)snonlin_lim; denom = ccomplex_ns::ccmax(a1, a2); }
                        // Find tauij values:
                        rho   = 0.5 * (q(j, kk, i, 1) + q(j, kk, i - 1, 1));
                        zk_sa = 0.0;
                        if (ivisc[0] == 6 || ivisc[0] == 7 || ivisc[0] == 10 || ivisc[0] == 15) {
                            zk_sa = 0.5 * (zksav(j, kk, i, 2) + zksav(j, kk, i - 1, 2));
                        }
                        t11   = 2.0 * rho * zk_sa * re / 3.0
                              - 8.0 * (double)c_nonlin * t(j + jiv, 24) / denom * (-s12*w12 - s13*w13);
                        t22_v = 2.0 * rho * zk_sa * re / 3.0
                              - 8.0 * (double)c_nonlin * t(j + jiv, 24) / denom * (s12*w12 - s23*w23);
                        t33_v = 2.0 * rho * zk_sa * re / 3.0
                              - 8.0 * (double)c_nonlin * t(j + jiv, 24) / denom * (s23*w23 + s13*w13);
                        t12_v = -4.0 * (double)c_nonlin * t(j + jiv, 24) / denom *
                                (s11*w12 - s22*w12 - s13*w23 - s23*w13);
                        t13_v = -4.0 * (double)c_nonlin * t(j + jiv, 24) / denom *
                                (s11*w13 + s12*w23 - s23*w12 - s33*w13);
                        t23_v = -4.0 * (double)c_nonlin * t(j + jiv, 24) / denom *
                                (s12*w13 + s22*w23 + s13*w12 - s33*w23);
                        t(j + jiv, 28) = xmre * (t(j + jiv, 25) * t11   + t(j + jiv, 26) * t12_v +
                                                  t(j + jiv, 27) * t13_v) / t(j + jiv, 15);
                        t(j + jiv, 29) = xmre * (t(j + jiv, 25) * t12_v + t(j + jiv, 26) * t22_v +
                                                  t(j + jiv, 27) * t23_v) / t(j + jiv, 15);
                        t(j + jiv, 30) = xmre * (t(j + jiv, 25) * t13_v + t(j + jiv, 26) * t23_v +
                                                  t(j + jiv, 27) * t33_v) / t(j + jiv, 15);
                        t(j + jiv, 31) = xmre / t(j + jiv, 15) * (
                            t(j + jiv, 25) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t11 +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t12_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t13_v) +
                            t(j + jiv, 26) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t12_v +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t22_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t23_v) +
                            t(j + jiv, 27) * (0.5 * (q(j, kk, i, 2) + q(j, kk, i - 1, 2)) * t13_v +
                                              0.5 * (q(j, kk, i, 3) + q(j, kk, i - 1, 3)) * t23_v +
                                              0.5 * (q(j, kk, i, 4) + q(j, kk, i - 1, 4)) * t33_v));
                    } // j loop interior
                } // end i==1/idim/else
            } // i loop
        } // kpl loop
    } // end i_nonlin!=0 branch

    // Final residual update for nonlinear branches (ivisc==11,12,13,14 or i_nonlin!=0)
    if (i_tauijs != 1 &&
        (ivisc[0] == 11 || ivisc[0] == 12 || ivisc[0] == 13 || ivisc[0] == 14 || i_nonlin != 0) &&
        level >= lglobal) {
        for (l = 2; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                // Take positive derivative
                t(izz, l) = t(izz + js - 1, l + 26) - t(izz, l + 26);
            }
        }
        for (kpl = 1; kpl <= npl; kpl++) {
            kk = k + kpl - 1;
            for (i = 1; i <= idim1; i++) {
                jiv = (i - 1) * jdim1 * npl + (kpl - 1) * jdim1 + 1;
                for (l = 2; l <= 5; l++) {
                    for (izz = 1; izz <= jdim1; izz++) {
                        res(izz, kk, i, l) = res(izz, kk, i, l) + t(izz + jiv - 1, l);
                    }
                }
            }
        }
    }

} // end ffluxv

} // namespace ffluxv_ns
