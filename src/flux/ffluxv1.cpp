// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "ffluxv1.h"
#include "ffluxv.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

namespace ffluxv1_ns {

// ffluxv: thin wrapper that delegates to ffluxv_ns::ffluxv
void ffluxv(int& k, int& npl, int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> res, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qi0, FortranArray4DRef<double> si,
            FortranArray3DRef<double> vol, FortranArray2DRef<double> t,
            int& nvtq, FortranArray2DRef<double> wi0,
            FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vmui,
            FortranArray4DRef<double> vi0, FortranArray3DRef<double> bci,
            FortranArray4DRef<double> zksav, FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> cmuv, FortranArray3DRef<double> voli0,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& iadv, int& nummem,
            FortranArray4DRef<double> ux)
{
    // ffluxv declares vmui as a flat 2D view (jdim-1, (kdim-1)*2); reshape the
    // wrapper's 3D (jdim-1,kdim-1,2) storage to match (same contiguous memory).
    FortranArray2DRef<double> vmui2(&vmui(1, 1, 1), jdim - 1, (kdim - 1) * 2);
    ffluxv_ns::ffluxv(k, npl, jdim, kdim, idim, res, q, qi0, si, vol, t,
                      nvtq, wi0, vist3d, vmui2, vi0, bci, zksav, ti0, cmuv,
                      voli0, nou, bou, nbuf, ibufdim, iadv, nummem, ux);
}

// ffluxv1: Calculate right-hand residual contributions in the I-direction
// due to the viscous terms missing in the thin-layer formulation.
void ffluxv1(int& k, int& npl, int& jdim, int& kdim, int& idim,
             FortranArray4DRef<double> res, FortranArray4DRef<double> q,
             FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
             FortranArray4DRef<double> qi0, FortranArray4DRef<double> sj,
             FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
             FortranArray3DRef<double> vol, FortranArray2DRef<double> t,
             int& nvtq, FortranArray2DRef<double> wi0,
             FortranArray3DRef<double> vist3d, FortranArray4DRef<double> vi0,
             FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
             FortranArray3DRef<double> bci, FortranArray3DRef<double> voli0,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& iadv)
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
    int*   ivisc   = cmn_reyue.ivisc;   // 0-based array
    int&   isklton = cmn_sklton.isklton;
    int&   i_tauijs = cmn_constit.i_tauijs;
    int&   level   = cmn_mgrd.level;
    int&   lglobal = cmn_mgrd.lglobal;

    // Local scalar declarations
    float coef_eddy;
    int   idim1, jdim1, kdim1;
    int   n, l0, jv, js, nn;
    double xmre, gpr, gm1pr, prtr, gprgm1;
    double c2b, c2bp;
    int   n1, l1;
    int   iviscc;

    // Loop variables
    int kpl, kk, i, j, l, izz;
    int jiv, ji, jc, mm, mb, m;
    double ab, bb, wi05, wi06, t5, t6, t1, t25, t24;
    double ux_v, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz;


    // coef_eddy setup
    coef_eddy = 1.0f;
    if (ivisc[0] >= 70) coef_eddy = 0.0f;
    if (i_tauijs == 1) coef_eddy = 0.0f;

    // Dimension indices
    idim1 = idim - 1;
    jdim1 = jdim - 1;
    kdim1 = kdim - 1;

    // n  : number of cell centers for npl planes
    // l0 : number of cell interfaces for npl planes
    // jv : number of cell centers (and interfaces) on a i=constant plane
    n  = npl * jdim1 * idim1;
    l0 = npl * jdim1 * idim;
    jv = npl * jdim1;
    js = jv + 1;
    nn = n - jv;

    xmre   = (double)xmach / (double)reue;
    gpr    = (double)gamma / (double)pr;
    gm1pr  = (double)gm1 * (double)pr;
    prtr   = (double)pr / (double)prt;
    gprgm1 = gpr / (double)gm1;

    if (isklton > 0 && k == 1 && iadv >= 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "   computing cross-derivative viscous fluxes, I-dir.");
    }

    // Store selected cell centered data in t array
    // t(j,22): density, t(j,18): u, t(j,19): v, t(j,20): w, t(j,16): c*c/(gm1*pr)
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

    // Store dependent variables on i=0 and i=idim boundaries in wi0
    for (m = 0; m <= 10; m += 10) {
        if (m == 0) {
            mm = 1;
            mb = 1;
        } else {
            mm = 3;
            mb = 2;
        }

        int jmin = 1;
        int jmax = jdim - 1;
        int kmin = k;
        int kmax = k + npl - 1;

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

    // t(17)/wi0(m+7): (u*u+v*v+w*w)/2
    // t(1)/wi0(m+1): 1/density
    // t(16)/wi0(m+6): c*c/(gm1*pr)
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
            wi0(j, m + 7)  = 0.5e0 * (wi0(j, m + 8)  * wi0(j, m + 8)
                                     + wi0(j, m + 9)  * wi0(j, m + 9)
                                     + wi0(j, m + 10) * wi0(j, m + 10));
            wi0(j, m + 11) = wi0(j, m + 7) * wi0(j, m + 12) + wi0(j, m + 6) / (double)gm1;
            wi0(j, m + 1)  = 1.e0 / wi0(j, m + 12);
            wi0(j, m + 6)  = gpr * wi0(j, m + 6) * wi0(j, m + 1) / (double)gm1;
        }
    }


    // t(7): laminar viscosity at cell centers (via sutherland relation)
    c2b  = (double)cbar / (double)tinf;
    c2bp = c2b + 1.e0;

    for (j = 1; j <= n; j++) {
        t5     = gm1pr * t(j, 16);
        t6     = std::sqrt(t5);
        t(j, 7) = c2bp * t5 * t6 / (c2b + t5);
    }

    // t(14): laminar viscosity values at cell interfaces
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
        t(izz, 14)   = c2bp * wi05 * wi06 / (c2b + wi05);
        ab  = 1.0 + wi0(izz, 15);
        bb  = 1.0 - wi0(izz, 15);
        wi05 = gm1pr * 0.5 * (ab * wi0(izz, 16) + bb * t(izz + n - jv, 16));
        wi06 = std::sqrt(wi05);
        t(izz + n, 14) = c2bp * wi05 * wi06 / (c2b + wi05);
    }

    // t(15): average jacobian (inverse volume) at cell interface
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

            // t(25)-t(27): components of grad(xie)
            // t1 = grad(xie)/J
            for (j = 1; j <= jdim1; j++) {
                t1          = si(j, kk, i, 4) * t(j + jiv, 15);
                t(j + jiv, 25) = si(j, kk, i, 1) * t1;
                t(j + jiv, 26) = si(j, kk, i, 2) * t1;
                t(j + jiv, 27) = si(j, kk, i, 3) * t1;
            }
        }
    }

    // t(28)-t(30): Components of grad(eta)
    // grad(eta) = 1/4 [(sj(j,k,i,l)*sj(j,k,i,4) + sj(j,k,i-1,l)*sj(j,k,i-1,4) +
    //             sj(j+1,k,i,l)*sj(j+1,k,i,4) + sj(j+1,k,i-1,l)*sj(j+1,k,i-1,4)] /V_int
    for (kpl = 1; kpl <= npl; kpl++) {
        kk  = k + kpl - 1;
        for (i = 1; i <= idim; i++) {
            jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
            if (i == 1) {
                for (j = 1; j <= jdim1; j++) {
                    t(jiv + j, 28) = 0.25 * t(j + jiv, 15) * (sj(j,   kk, i, 1) * sj(j,   kk, i, 4)
                                                               + sj(j,   kk, i, 1) * sj(j,   kk, i, 4)
                                                               + sj(j+1, kk, i, 1) * sj(j+1, kk, i, 4)
                                                               + sj(j+1, kk, i, 1) * sj(j+1, kk, i, 4));
                    t(jiv + j, 29) = 0.25 * t(j + jiv, 15) * (sj(j,   kk, i, 2) * sj(j,   kk, i, 4)
                                                               + sj(j,   kk, i, 2) * sj(j,   kk, i, 4)
                                                               + sj(j+1, kk, i, 2) * sj(j+1, kk, i, 4)
                                                               + sj(j+1, kk, i, 2) * sj(j+1, kk, i, 4));
                    t(jiv + j, 30) = 0.25 * t(j + jiv, 15) * (sj(j,   kk, i, 3) * sj(j,   kk, i, 4)
                                                               + sj(j,   kk, i, 3) * sj(j,   kk, i, 4)
                                                               + sj(j+1, kk, i, 3) * sj(j+1, kk, i, 4)
                                                               + sj(j+1, kk, i, 3) * sj(j+1, kk, i, 4));
                }
            } else if (i == idim) {
                for (j = 1; j <= jdim1; j++) {
                    t(jiv + j, 28) = 0.25 * t(j + jiv, 15) * (sj(j,   kk, i-1, 1) * sj(j,   kk, i-1, 4)
                                                               + sj(j,   kk, i-1, 1) * sj(j,   kk, i-1, 4)
                                                               + sj(j+1, kk, i-1, 1) * sj(j+1, kk, i-1, 4)
                                                               + sj(j+1, kk, i-1, 1) * sj(j+1, kk, i-1, 4));
                    t(jiv + j, 29) = 0.25 * t(j + jiv, 15) * (sj(j,   kk, i-1, 2) * sj(j,   kk, i-1, 4)
                                                               + sj(j,   kk, i-1, 2) * sj(j,   kk, i-1, 4)
                                                               + sj(j+1, kk, i-1, 2) * sj(j+1, kk, i-1, 4)
                                                               + sj(j+1, kk, i-1, 2) * sj(j+1, kk, i-1, 4));
                    t(jiv + j, 30) = 0.25 * t(j + jiv, 15) * (sj(j,   kk, i-1, 3) * sj(j,   kk, i-1, 4)
                                                               + sj(j,   kk, i-1, 3) * sj(j,   kk, i-1, 4)
                                                               + sj(j+1, kk, i-1, 3) * sj(j+1, kk, i-1, 4)
                                                               + sj(j+1, kk, i-1, 3) * sj(j+1, kk, i-1, 4));
                }
            } else {
                // general case
                for (j = 1; j <= jdim1; j++) {
                    t(jiv + j, 28) = 0.25 * t(j + jiv, 15) * (sj(j,   kk, i,   1) * sj(j,   kk, i,   4)
                                                               + sj(j,   kk, i-1, 1) * sj(j,   kk, i-1, 4)
                                                               + sj(j+1, kk, i,   1) * sj(j+1, kk, i,   4)
                                                               + sj(j+1, kk, i-1, 1) * sj(j+1, kk, i-1, 4));
                    t(jiv + j, 29) = 0.25 * t(j + jiv, 15) * (sj(j,   kk, i,   2) * sj(j,   kk, i,   4)
                                                               + sj(j,   kk, i-1, 2) * sj(j,   kk, i-1, 4)
                                                               + sj(j+1, kk, i,   2) * sj(j+1, kk, i,   4)
                                                               + sj(j+1, kk, i-1, 2) * sj(j+1, kk, i-1, 4));
                    t(jiv + j, 30) = 0.25 * t(j + jiv, 15) * (sj(j,   kk, i,   3) * sj(j,   kk, i,   4)
                                                               + sj(j,   kk, i-1, 3) * sj(j,   kk, i-1, 4)
                                                               + sj(j+1, kk, i,   3) * sj(j+1, kk, i,   4)
                                                               + sj(j+1, kk, i-1, 3) * sj(j+1, kk, i-1, 4));
                }
            }
        }
    }


    // t(31),t(32),t(2): Components of grad(zeta)
    // grad(zeta) = 1/4 [(sk(j,k,i,l)*sk(j,k,i,4) + sk(j,k,i-1,l)*sk(j,k,i-1,4) +
    //              sk(j,k+1,i,l)*sk(j,k+1,i,4) + sk(j,k+1,i-1,l)*sk(j,k+1,i-1,4)] /V_int
    for (kpl = 1; kpl <= npl; kpl++) {
        kk  = k + kpl - 1;
        for (i = 1; i <= idim; i++) {
            jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
            if (i == 1) {
                for (j = 1; j <= jdim1; j++) {
                    t(jiv + j, 31) = 0.25 * t(j + jiv, 15) * (sk(j, kk,   i, 1) * sk(j, kk,   i, 4)
                                                               + sk(j, kk,   i, 1) * sk(j, kk,   i, 4)
                                                               + sk(j, kk+1, i, 1) * sk(j, kk+1, i, 4)
                                                               + sk(j, kk+1, i, 1) * sk(j, kk+1, i, 4));
                    t(jiv + j, 32) = 0.25 * t(j + jiv, 15) * (sk(j, kk,   i, 2) * sk(j, kk,   i, 4)
                                                               + sk(j, kk,   i, 2) * sk(j, kk,   i, 4)
                                                               + sk(j, kk+1, i, 2) * sk(j, kk+1, i, 4)
                                                               + sk(j, kk+1, i, 2) * sk(j, kk+1, i, 4));
                    t(jiv + j, 2)  = 0.25 * t(j + jiv, 15) * (sk(j, kk,   i, 3) * sk(j, kk,   i, 4)
                                                               + sk(j, kk,   i, 3) * sk(j, kk,   i, 4)
                                                               + sk(j, kk+1, i, 3) * sk(j, kk+1, i, 4)
                                                               + sk(j, kk+1, i, 3) * sk(j, kk+1, i, 4));
                }
            } else if (i == idim) {
                for (j = 1; j <= jdim1; j++) {
                    t(jiv + j, 31) = 0.25 * t(j + jiv, 15) * (sk(j, kk,   i-1, 1) * sk(j, kk,   i-1, 4)
                                                               + sk(j, kk,   i-1, 1) * sk(j, kk,   i-1, 4)
                                                               + sk(j, kk+1, i-1, 1) * sk(j, kk+1, i-1, 4)
                                                               + sk(j, kk+1, i-1, 1) * sk(j, kk+1, i-1, 4));
                    t(jiv + j, 32) = 0.25 * t(j + jiv, 15) * (sk(j, kk,   i-1, 2) * sk(j, kk,   i-1, 4)
                                                               + sk(j, kk,   i-1, 2) * sk(j, kk,   i-1, 4)
                                                               + sk(j, kk+1, i-1, 2) * sk(j, kk+1, i-1, 4)
                                                               + sk(j, kk+1, i-1, 2) * sk(j, kk+1, i-1, 4));
                    t(jiv + j, 2)  = 0.25 * t(j + jiv, 15) * (sk(j, kk,   i-1, 3) * sk(j, kk,   i-1, 4)
                                                               + sk(j, kk,   i-1, 3) * sk(j, kk,   i-1, 4)
                                                               + sk(j, kk+1, i-1, 3) * sk(j, kk+1, i-1, 4)
                                                               + sk(j, kk+1, i-1, 3) * sk(j, kk+1, i-1, 4));
                }
            } else {
                // general case
                for (j = 1; j <= jdim1; j++) {
                    t(jiv + j, 31) = 0.25 * t(j + jiv, 15) * (sk(j, kk,   i,   1) * sk(j, kk,   i,   4)
                                                               + sk(j, kk,   i-1, 1) * sk(j, kk,   i-1, 4)
                                                               + sk(j, kk+1, i,   1) * sk(j, kk+1, i,   4)
                                                               + sk(j, kk+1, i-1, 1) * sk(j, kk+1, i-1, 4));
                    t(jiv + j, 32) = 0.25 * t(j + jiv, 15) * (sk(j, kk,   i,   2) * sk(j, kk,   i,   4)
                                                               + sk(j, kk,   i-1, 2) * sk(j, kk,   i-1, 4)
                                                               + sk(j, kk+1, i,   2) * sk(j, kk+1, i,   4)
                                                               + sk(j, kk+1, i-1, 2) * sk(j, kk+1, i-1, 4));
                    t(jiv + j, 2)  = 0.25 * t(j + jiv, 15) * (sk(j, kk,   i,   3) * sk(j, kk,   i,   4)
                                                               + sk(j, kk,   i-1, 3) * sk(j, kk,   i-1, 4)
                                                               + sk(j, kk+1, i,   3) * sk(j, kk+1, i,   4)
                                                               + sk(j, kk+1, i-1, 3) * sk(j, kk+1, i-1, 4));
                }
            }
        }
    }


    // Gradients at cell interfaces (eta direction)
    // t(3): d(c*c/(gm1*pr))/d(eta), t(4): d(u)/d(eta), t(5): d(v)/d(eta), t(6): d(w)/d(eta)
    // Recall bcj=0 => cell data, bcj=1 => face data
    for (kpl = 1; kpl <= npl; kpl++) {
        kk = k + kpl - 1;
        for (i = 1; i <= idim; i++) {
            jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;

            if (i == 1) {
                // replace (i-1) data by boundary data
                for (j = 1; j <= jdim1; j++) {
                    if (j == 1) {
                        t(jiv+j, 3) = 0.25*gprgm1*(
                            q(j+1,kk,i,5)/q(j+1,kk,i,1) -
                            (1.+bcj(kk,i,1))*qj0(kk,i,5,1)/qj0(kk,i,1,1) +
                            bcj(kk,i,1)*q(j,kk,i,5)/q(j,kk,i,1) +
                            2.*(1.+bci(j+1,kk,1))*qi0(j+1,kk,5,1)/qi0(j+1,kk,1,1) -
                            2.*bci(j+1,kk,1)*q(j+1,kk,i,5)/q(j+1,kk,i,1) -
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,5,1)/qi0(j,kk,1,1) +
                            2.*bci(j,kk,1)*q(j,kk,i,5)/q(j,kk,i,1) );
                        t(jiv+j, 4) = 0.25*(
                            q(j+1,kk,i,2) -
                            (1.+bcj(kk,i,1))*qj0(kk,i,2,1) +
                            bcj(kk,i,1)*q(j,kk,i,2) +
                            2.*(1.+bci(j+1,kk,1))*qi0(j+1,kk,2,1) -
                            2.*bci(j+1,kk,1)*q(j+1,kk,i,2) -
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,2,1) +
                            2.*bci(j,kk,1)*q(j,kk,i,2) );
                        t(jiv+j, 5) = 0.25*(
                            q(j+1,kk,i,3) -
                            (1.+bcj(kk,i,1))*qj0(kk,i,3,1) +
                            bcj(kk,i,1)*q(j,kk,i,3) +
                            2.*(1.+bci(j+1,kk,1))*qi0(j+1,kk,3,1) -
                            2.*bci(j+1,kk,1)*q(j+1,kk,i,3) -
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,3,1) +
                            2.*bci(j,kk,1)*q(j,kk,i,3) );
                        t(jiv+j, 6) = 0.25*(
                            q(j+1,kk,i,4) -
                            (1.+bcj(kk,i,1))*qj0(kk,i,4,1) +
                            bcj(kk,i,1)*q(j,kk,i,4) +
                            2.*(1.+bci(j+1,kk,1))*qi0(j+1,kk,4,1) -
                            2.*bci(j+1,kk,1)*q(j+1,kk,i,4) -
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,4,1) +
                            2.*bci(j,kk,1)*q(j,kk,i,4) );
                    } else if (j == jdim1) {
                        t(jiv+j, 3) = 0.25*gprgm1*(
                            (1.+bcj(kk,i,2))*qj0(kk,i,5,3)/qj0(kk,i,1,3) -
                            bcj(kk,i,2)*q(j,kk,i,5)/q(j,kk,i,1) -
                            q(j-1,kk,i,5)/q(j-1,kk,i,1) +
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,5,1)/qi0(j,kk,1,1) -
                            2.*bci(j,kk,1)*q(j,kk,i,5)/q(j,kk,i,1) -
                            2.*(1.+bci(j-1,kk,1))*qi0(j-1,kk,5,1)/qi0(j-1,kk,1,1) +
                            2.*bci(j-1,kk,1)*q(j-1,kk,i,5)/q(j-1,kk,i,1) );
                        t(jiv+j, 4) = 0.25*(
                            (1.+bcj(kk,i,2))*qj0(kk,i,2,3) -
                            bcj(kk,i,2)*q(j,kk,i,2) -
                            q(j-1,kk,i,2) +
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,2,1) -
                            2.*bci(j,kk,1)*q(j,kk,i,2) -
                            2.*(1.+bci(j-1,kk,1))*qi0(j-1,kk,2,1) +
                            2.*bci(j-1,kk,1)*q(j-1,kk,i,2) );
                        t(jiv+j, 5) = 0.25*(
                            (1.+bcj(kk,i,2))*qj0(kk,i,3,3) -
                            bcj(kk,i,2)*q(j,kk,i,3) -
                            q(j-1,kk,i,3) +
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,3,1) -
                            2.*bci(j,kk,1)*q(j,kk,i,3) -
                            2.*(1.+bci(j-1,kk,1))*qi0(j-1,kk,3,1) +
                            2.*bci(j-1,kk,1)*q(j-1,kk,i,3) );
                        t(jiv+j, 6) = 0.25*(
                            (1.+bcj(kk,i,2))*qj0(kk,i,4,3) -
                            bcj(kk,i,2)*q(j,kk,i,4) -
                            q(j-1,kk,i,4) +
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,4,1) -
                            2.*bci(j,kk,1)*q(j,kk,i,4) -
                            2.*(1.+bci(j-1,kk,1))*qi0(j-1,kk,4,1) +
                            2.*bci(j-1,kk,1)*q(j-1,kk,i,4) );
                    } else {
                        // general case for i = 1
                        t(jiv+j, 3) = 0.25*gprgm1*(
                            q(j+1,kk,i,5)/q(j+1,kk,i,1) -
                            q(j-1,kk,i,5)/q(j-1,kk,i,1) +
                            (1.+bci(j+1,kk,1))*qi0(j+1,kk,5,1)/qi0(j+1,kk,1,1) -
                            bci(j+1,kk,1)*q(j+1,kk,i,5)/q(j+1,kk,i,1) -
                            (1.+bci(j-1,kk,1))*qi0(j-1,kk,5,1)/qi0(j-1,kk,1,1) +
                            bci(j-1,kk,1)*q(j-1,kk,i,5)/q(j-1,kk,i,1) );
                        t(jiv+j, 4) = 0.25*(
                            q(j+1,kk,i,2) - q(j-1,kk,i,2) +
                            (1.+bci(j+1,kk,1))*qi0(j+1,kk,2,1) -
                            bci(j+1,kk,1)*q(j+1,kk,i,2) -
                            (1.+bci(j-1,kk,1))*qi0(j-1,kk,2,1) +
                            bci(j-1,kk,1)*q(j-1,kk,i,2) );
                        t(jiv+j, 5) = 0.25*(
                            q(j+1,kk,i,3) - q(j-1,kk,i,3) +
                            (1.+bci(j+1,kk,1))*qi0(j+1,kk,3,1) -
                            bci(j+1,kk,1)*q(j+1,kk,i,3) -
                            (1.+bci(j-1,kk,1))*qi0(j-1,kk,3,1) +
                            bci(j-1,kk,1)*q(j-1,kk,i,3) );
                        t(jiv+j, 6) = 0.25*(
                            q(j+1,kk,i,4) - q(j-1,kk,i,4) +
                            (1.+bci(j+1,kk,1))*qi0(j+1,kk,4,1) -
                            bci(j+1,kk,1)*q(j+1,kk,i,4) -
                            (1.+bci(j-1,kk,1))*qi0(j-1,kk,4,1) +
                            bci(j-1,kk,1)*q(j-1,kk,i,4) );
                    }
                } // j loop, i==1


            } else if (i == idim) {
                // replace i data by boundary data
                for (j = 1; j <= jdim1; j++) {
                    if (j == 1) {
                        t(jiv+j, 3) = 0.25*gprgm1*(
                            2.*(1.+bci(j+1,kk,2))*qi0(j+1,kk,5,3)/qi0(j+1,kk,1,3) -
                            2.*bci(j+1,kk,2)*q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,5,3)/qi0(j,kk,1,3) +
                            2.*bci(j,kk,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) +
                            q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
                            (1.+bcj(kk,i-1,1))*qj0(kk,i-1,5,1)/qj0(kk,i-1,1,1) +
                            bcj(kk,i-1,1)*q(j,kk,i-1,5)/q(j,kk,i-1,1) );
                        t(jiv+j, 4) = 0.25*(
                            2.*(1.+bci(j+1,kk,2))*qi0(j+1,kk,2,3) -
                            2.*bci(j+1,kk,2)*q(j+1,kk,i-1,2) -
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,2,3) +
                            2.*bci(j,kk,2)*q(j,kk,i-1,2) +
                            q(j+1,kk,i-1,2) -
                            (1.+bcj(kk,i-1,1))*qj0(kk,i-1,2,1) +
                            bcj(kk,i-1,1)*q(j,kk,i-1,2) );
                        t(jiv+j, 5) = 0.25*(
                            2.*(1.+bci(j+1,kk,2))*qi0(j+1,kk,3,3) -
                            2.*bci(j+1,kk,2)*q(j+1,kk,i-1,3) -
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,3,3) +
                            2.*bci(j,kk,2)*q(j,kk,i-1,3) +
                            q(j+1,kk,i-1,3) -
                            (1.+bcj(kk,i-1,1))*qj0(kk,i-1,3,1) +
                            bcj(kk,i-1,1)*q(j,kk,i-1,3) );
                        t(jiv+j, 6) = 0.25*(
                            2.*(1.+bci(j+1,kk,2))*qi0(j+1,kk,4,3) -
                            2.*bci(j+1,kk,2)*q(j+1,kk,i-1,4) -
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,4,3) +
                            2.*bci(j,kk,2)*q(j,kk,i-1,4) +
                            q(j+1,kk,i-1,4) -
                            (1.+bcj(kk,i-1,1))*qj0(kk,i-1,4,1) +
                            bcj(kk,i-1,1)*q(j,kk,i-1,4) );
                    } else if (j == jdim1) {
                        t(jiv+j, 3) = 0.25*gprgm1*(
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,5,3)/qi0(j,kk,1,3) -
                            2.*bci(j,kk,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) -
                            2.*(1.+bci(j-1,kk,2))*qi0(j-1,kk,5,3)/qi0(j-1,kk,1,3) +
                            2.*bci(j-1,kk,2)*q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) +
                            (1.+bcj(kk,i-1,2))*qj0(kk,i-1,5,3)/qj0(kk,i-1,1,3) -
                            bcj(kk,i-1,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) -
                            q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) );
                        t(jiv+j, 4) = 0.25*(
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,2,3) -
                            2.*bci(j,kk,2)*q(j,kk,i-1,2) -
                            2.*(1.+bci(j-1,kk,2))*qi0(j-1,kk,2,3) +
                            2.*bci(j-1,kk,2)*q(j-1,kk,i-1,2) +
                            (1.+bcj(kk,i-1,2))*qj0(kk,i-1,2,3) -
                            bcj(kk,i-1,2)*q(j,kk,i-1,2) -
                            q(j-1,kk,i-1,2) );
                        t(jiv+j, 5) = 0.25*(
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,3,3) -
                            2.*bci(j,kk,2)*q(j,kk,i-1,3) -
                            2.*(1.+bci(j-1,kk,2))*qi0(j-1,kk,3,3) +
                            2.*bci(j-1,kk,2)*q(j-1,kk,i-1,3) +
                            (1.+bcj(kk,i-1,2))*qj0(kk,i-1,3,3) -
                            bcj(kk,i-1,2)*q(j,kk,i-1,3) -
                            q(j-1,kk,i-1,3) );
                        t(jiv+j, 6) = 0.25*(
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,4,3) -
                            2.*bci(j,kk,2)*q(j,kk,i-1,4) -
                            2.*(1.+bci(j-1,kk,2))*qi0(j-1,kk,4,3) +
                            2.*bci(j-1,kk,2)*q(j-1,kk,i-1,4) +
                            (1.+bcj(kk,i-1,2))*qj0(kk,i-1,4,3) -
                            bcj(kk,i-1,2)*q(j,kk,i-1,4) -
                            q(j-1,kk,i-1,4) );
                    } else {
                        // general case for i = idim
                        t(jiv+j, 3) = 0.25*gprgm1*(
                            (1.+bci(j+1,kk,2))*qi0(j+1,kk,5,3)/qi0(j+1,kk,1,3) -
                            bci(j+1,kk,2)*q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
                            (1.+bci(j-1,kk,2))*qi0(j-1,kk,5,3)/qi0(j-1,kk,1,3) +
                            bci(j-1,kk,2)*q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) +
                            q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
                            q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) );
                        t(jiv+j, 4) = 0.25*(
                            (1.+bci(j+1,kk,2))*qi0(j+1,kk,2,3) -
                            bci(j+1,kk,2)*q(j+1,kk,i-1,2) -
                            (1.+bci(j-1,kk,2))*qi0(j-1,kk,2,3) +
                            bci(j-1,kk,2)*q(j-1,kk,i-1,2) +
                            q(j+1,kk,i-1,2) - q(j-1,kk,i-1,2) );
                        t(jiv+j, 5) = 0.25*(
                            (1.+bci(j+1,kk,2))*qi0(j+1,kk,3,3) -
                            bci(j+1,kk,2)*q(j+1,kk,i-1,3) -
                            (1.+bci(j-1,kk,2))*qi0(j-1,kk,3,3) +
                            bci(j-1,kk,2)*q(j-1,kk,i-1,3) +
                            q(j+1,kk,i-1,3) - q(j-1,kk,i-1,3) );
                        t(jiv+j, 6) = 0.25*(
                            (1.+bci(j+1,kk,2))*qi0(j+1,kk,4,3) -
                            bci(j+1,kk,2)*q(j+1,kk,i-1,4) -
                            (1.+bci(j-1,kk,2))*qi0(j-1,kk,4,3) +
                            bci(j-1,kk,2)*q(j-1,kk,i-1,4) +
                            q(j+1,kk,i-1,4) - q(j-1,kk,i-1,4) );
                    }
                } // j loop, i==idim


            } else {
                // general case for i between 1 and idim
                for (j = 1; j <= jdim1; j++) {
                    if (j == 1) {
                        t(jiv+j, 3) = 0.25*gprgm1*(
                            q(j+1,kk,i,5)/q(j+1,kk,i,1) -
                            (1.+bcj(kk,i,1))*qj0(kk,i,5,1)/qj0(kk,i,1,1) +
                            bcj(kk,i,1)*q(j,kk,i,5)/q(j,kk,i,1) +
                            q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
                            (1.+bcj(kk,i-1,1))*qj0(kk,i-1,5,1)/qj0(kk,i-1,1,1) +
                            bcj(kk,i-1,1)*q(j,kk,i-1,5)/q(j,kk,i-1,1) );
                        t(jiv+j, 4) = 0.25*(
                            q(j+1,kk,i,2) -
                            (1.+bcj(kk,i,1))*qj0(kk,i,2,1) +
                            bcj(kk,i,1)*q(j,kk,i,2) +
                            q(j+1,kk,i-1,2) -
                            (1.+bcj(kk,i-1,1))*qj0(kk,i-1,2,1) +
                            bcj(kk,i-1,1)*q(j,kk,i-1,2) );
                        t(jiv+j, 5) = 0.25*(
                            q(j+1,kk,i,3) -
                            (1.+bcj(kk,i,1))*qj0(kk,i,3,1) +
                            bcj(kk,i,1)*q(j,kk,i,3) +
                            q(j+1,kk,i-1,3) -
                            (1.+bcj(kk,i-1,1))*qj0(kk,i-1,3,1) +
                            bcj(kk,i-1,1)*q(j,kk,i-1,3) );
                        t(jiv+j, 6) = 0.25*(
                            q(j+1,kk,i,4) -
                            (1.+bcj(kk,i,1))*qj0(kk,i,4,1) +
                            bcj(kk,i,1)*q(j,kk,i,4) +
                            q(j+1,kk,i-1,4) -
                            (1.+bcj(kk,i-1,1))*qj0(kk,i-1,4,1) +
                            bcj(kk,i-1,1)*q(j,kk,i-1,4) );
                    } else if (j == jdim1) {
                        t(jiv+j, 3) = 0.25*gprgm1*(
                            (1.+bcj(kk,i,2))*qj0(kk,i,5,3)/qj0(kk,i,1,3) -
                            bcj(kk,i,2)*q(j,kk,i,5)/q(j,kk,i,1) -
                            q(j-1,kk,i,5)/q(j-1,kk,i,1) +
                            (1.+bcj(kk,i-1,2))*qj0(kk,i-1,5,3)/qj0(kk,i-1,1,3) -
                            bcj(kk,i-1,2)*q(j-1,kk,i,5)/q(j-1,kk,i,1) -
                            q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) );
                        t(jiv+j, 4) = 0.25*(
                            (1.+bcj(kk,i,2))*qj0(kk,i,2,3) -
                            bcj(kk,i,2)*q(j,kk,i,2) -
                            q(j-1,kk,i,2) +
                            (1.+bcj(kk,i-1,2))*qj0(kk,i-1,2,3) -
                            bcj(kk,i-1,2)*q(j-1,kk,i,2) -
                            q(j-1,kk,i-1,2) );
                        t(jiv+j, 5) = 0.25*(
                            (1.+bcj(kk,i,2))*qj0(kk,i,3,3) -
                            bcj(kk,i,2)*q(j,kk,i,3) -
                            q(j-1,kk,i,3) +
                            (1.+bcj(kk,i-1,2))*qj0(kk,i-1,3,3) -
                            bcj(kk,i-1,2)*q(j-1,kk,i,3) -
                            q(j-1,kk,i-1,3) );
                        t(jiv+j, 6) = 0.25*(
                            (1.+bcj(kk,i,2))*qj0(kk,i,4,3) -
                            bcj(kk,i,2)*q(j,kk,i,4) -
                            q(j-1,kk,i,4) +
                            (1.+bcj(kk,i-1,2))*qj0(kk,i-1,4,3) -
                            bcj(kk,i-1,2)*q(j-1,kk,i,4) -
                            q(j-1,kk,i-1,4) );
                    } else {
                        // general case
                        t(jiv+j, 3) = 0.25*gprgm1*(
                            q(j+1,kk,i,5)/q(j+1,kk,i,1) -
                            q(j-1,kk,i,5)/q(j-1,kk,i,1) +
                            q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
                            q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) );
                        t(jiv+j, 4) = 0.25*(
                            q(j+1,kk,i,2)   - q(j-1,kk,i,2) +
                            q(j+1,kk,i-1,2) - q(j-1,kk,i-1,2) );
                        t(jiv+j, 5) = 0.25*(
                            q(j+1,kk,i,3)   - q(j-1,kk,i,3) +
                            q(j+1,kk,i-1,3) - q(j-1,kk,i-1,3) );
                        t(jiv+j, 6) = 0.25*(
                            q(j+1,kk,i,4)   - q(j-1,kk,i,4) +
                            q(j+1,kk,i-1,4) - q(j-1,kk,i-1,4) );
                    }
                } // j loop, general i
            }
        } // i loop
    } // kpl loop (eta gradients)


    // Gradients at cell interfaces (zeta direction)
    // t(7): d(c*c/(gm1*pr))/d(zeta), t(8): d(u)/d(zeta), t(9): d(v)/d(zeta), t(10): d(w)/d(zeta)
    // Recall bck=0 => cell data, bck=1 => face data
    for (kpl = 1; kpl <= npl; kpl++) {
        kk = k + kpl - 1;
        if (kk == 1) {
            // replace (kk-1) data by boundary data
            for (i = 1; i <= idim; i++) {
                jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
                if (i == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        t(jiv+j, 7) = 0.25*gprgm1*(
                            q(j,kk+1,i,5)/q(j,kk+1,i,1) -
                            (1.+bck(j,i,1))*qk0(j,i,5,1)/qk0(j,i,1,1) +
                            bck(j,i,1)*q(j,kk,i,5)/q(j,kk,i,1) +
                            2.*(1.+bci(j,kk+1,1))*qi0(j,kk+1,5,1)/qi0(j,kk+1,1,1) -
                            2.*bci(j,kk+1,1)*q(j,kk+1,i,5)/q(j,kk+1,i,1) -
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,5,1)/qi0(j,kk,1,1) +
                            2.*bci(j,kk,1)*q(j,kk,i,5)/q(j,kk,i,1) );
                        t(jiv+j, 8) = 0.25*(
                            q(j,kk+1,i,2) -
                            (1.+bck(j,i,1))*qk0(j,i,2,1) +
                            bck(j,i,1)*q(j,kk,i,2) +
                            2.*(1.+bci(j,kk+1,1))*qi0(j,kk+1,2,1) -
                            2.*bci(j,kk+1,1)*q(j,kk+1,i,2) -
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,2,1) +
                            2.*bci(j,kk,1)*q(j,kk,i,2) );
                        t(jiv+j, 9) = 0.25*(
                            q(j,kk+1,i,3) -
                            (1.+bck(j,i,1))*qk0(j,i,3,1) +
                            bck(j,i,1)*q(j,kk,i,3) +
                            2.*(1.+bci(j,kk+1,1))*qi0(j,kk+1,3,1) -
                            2.*bci(j,kk+1,1)*q(j,kk+1,i,3) -
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,3,1) +
                            2.*bci(j,kk,1)*q(j,kk,i,3) );
                        t(jiv+j, 10) = 0.25*(
                            q(j,kk+1,i,4) -
                            (1.+bck(j,i,1))*qk0(j,i,4,1) +
                            bck(j,i,1)*q(j,kk,i,4) +
                            2.*(1.+bci(j,kk+1,1))*qi0(j,kk+1,4,1) -
                            2.*bci(j,kk+1,1)*q(j,kk+1,i,4) -
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,4,1) +
                            2.*bci(j,kk,1)*q(j,kk,i,4) );
                    }
                } else if (i == idim) {
                    for (j = 1; j <= jdim1; j++) {
                        t(jiv+j, 7) = 0.25*gprgm1*(
                            2.*(1.+bci(j,kk+1,2))*qi0(j,kk+1,5,3)/qi0(j,kk+1,1,3) -
                            2.*bci(j,kk+1,2)*q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,5,3)/qi0(j,kk,1,3) +
                            2.*bci(j,kk,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) +
                            q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
                            (1.+bck(j,i-1,1))*qk0(j,i-1,5,1)/qk0(j,i-1,1,1) +
                            bck(j,i-1,1)*q(j,kk,i-1,5)/q(j,kk,i-1,1) );
                        t(jiv+j, 8) = 0.25*(
                            2.*(1.+bci(j,kk+1,2))*qi0(j,kk+1,2,3) -
                            2.*bci(j,kk+1,2)*q(j,kk+1,i-1,2) -
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,2,3) +
                            2.*bci(j,kk,2)*q(j,kk,i-1,2) +
                            q(j,kk+1,i-1,2) -
                            (1.+bck(j,i-1,1))*qk0(j,i-1,2,1) +
                            bck(j,i-1,1)*q(j,kk,i-1,2) );
                        t(jiv+j, 9) = 0.25*(
                            2.*(1.+bci(j,kk+1,2))*qi0(j,kk+1,3,3) -
                            2.*bci(j,kk+1,2)*q(j,kk+1,i-1,3) -
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,3,3) +
                            2.*bci(j,kk,2)*q(j,kk,i-1,3) +
                            q(j,kk+1,i-1,3) -
                            (1.+bck(j,i-1,1))*qk0(j,i-1,3,1) +
                            bck(j,i-1,1)*q(j,kk,i-1,3) );
                        t(jiv+j, 10) = 0.25*(
                            2.*(1.+bci(j,kk+1,2))*qi0(j,kk+1,4,3) -
                            2.*bci(j,kk+1,2)*q(j,kk+1,i-1,4) -
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,4,3) +
                            2.*bci(j,kk,2)*q(j,kk,i-1,4) +
                            q(j,kk+1,i-1,4) -
                            (1.+bck(j,i-1,1))*qk0(j,i-1,4,1) +
                            bck(j,i-1,1)*q(j,kk,i-1,4) );
                    }
                } else {
                    // general case for kk = 1
                    for (j = 1; j <= jdim1; j++) {
                        t(jiv+j, 7) = 0.25*gprgm1*(
                            q(j,kk+1,i,5)/q(j,kk+1,i,1) -
                            (1.+bck(j,i,1))*qk0(j,i,5,1)/qk0(j,i,1,1) +
                            bck(j,i,1)*q(j,kk,i,5)/q(j,kk,i,1) +
                            q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
                            (1.+bck(j,i-1,1))*qk0(j,i-1,5,1)/qk0(j,i-1,1,1) +
                            bck(j,i-1,1)*q(j,kk,i-1,5)/q(j,kk,i-1,1) );
                        t(jiv+j, 8) = 0.25*(
                            q(j,kk+1,i,2) -
                            (1.+bck(j,i,1))*qk0(j,i,2,1) +
                            bck(j,i,1)*q(j,kk,i,2) +
                            q(j,kk+1,i-1,2) -
                            (1.+bck(j,i-1,1))*qk0(j,i-1,2,1) +
                            bck(j,i-1,1)*q(j,kk,i-1,2) );
                        t(jiv+j, 9) = 0.25*(
                            q(j,kk+1,i,3) -
                            (1.+bck(j,i,1))*qk0(j,i,3,1) +
                            bck(j,i,1)*q(j,kk,i,3) +
                            q(j,kk+1,i-1,3) -
                            (1.+bck(j,i-1,1))*qk0(j,i-1,3,1) +
                            bck(j,i-1,1)*q(j,kk,i-1,3) );
                        t(jiv+j, 10) = 0.25*(
                            q(j,kk+1,i,4) -
                            (1.+bck(j,i,1))*qk0(j,i,4,1) +
                            bck(j,i,1)*q(j,kk,i,4) +
                            q(j,kk+1,i-1,4) -
                            (1.+bck(j,i-1,1))*qk0(j,i-1,4,1) +
                            bck(j,i-1,1)*q(j,kk,i-1,4) );
                    }
                }
            } // i loop, kk==1


        } else if (kk == kdim1) {
            // replace (kk+1) data by boundary data
            for (i = 1; i <= idim; i++) {
                jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
                if (i == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        t(jiv+j, 7) = 0.25*gprgm1*(
                            (1.+bck(j,i,2))*qk0(j,i,5,3)/qk0(j,i,1,3) -
                            bck(j,i,2)*q(j,kk,i,5)/q(j,kk,i,1) -
                            q(j,kk-1,i,5)/q(j,kk-1,i,1) +
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,5,1)/qi0(j,kk,1,1) -
                            2.*bci(j,kk,1)*q(j,kk,i,5)/q(j,kk,i,1) -
                            2.*(1.+bci(j,kk-1,1))*qi0(j,kk-1,5,1)/qi0(j,kk-1,1,1) +
                            2.*bci(j,kk-1,1)*q(j,kk-1,i,5)/q(j,kk-1,i,1) );
                        t(jiv+j, 8) = 0.25*(
                            (1.+bck(j,i,2))*qk0(j,i,2,3) -
                            bck(j,i,2)*q(j,kk,i,2) -
                            q(j,kk-1,i,2) +
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,2,1) -
                            2.*bci(j,kk,1)*q(j,kk,i,2) -
                            2.*(1.+bci(j,kk-1,1))*qi0(j,kk-1,2,1) +
                            2.*bci(j,kk-1,1)*q(j,kk-1,i,2) );
                        t(jiv+j, 9) = 0.25*(
                            (1.+bck(j,i,2))*qk0(j,i,3,3) -
                            bck(j,i,2)*q(j,kk,i,3) -
                            q(j,kk-1,i,3) +
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,3,1) -
                            2.*bci(j,kk,1)*q(j,kk,i,3) -
                            2.*(1.+bci(j,kk-1,1))*qi0(j,kk-1,3,1) +
                            2.*bci(j,kk-1,1)*q(j,kk-1,i,3) );
                        t(jiv+j, 10) = 0.25*(
                            (1.+bck(j,i,2))*qk0(j,i,4,3) -
                            bck(j,i,2)*q(j,kk,i,4) -
                            q(j,kk-1,i,4) +
                            2.*(1.+bci(j,kk,1))*qi0(j,kk,4,1) -
                            2.*bci(j,kk,1)*q(j,kk,i,4) -
                            2.*(1.+bci(j,kk-1,1))*qi0(j,kk-1,4,1) +
                            2.*bci(j,kk-1,1)*q(j,kk-1,i,4) );
                    }
                } else if (i == idim) {
                    for (j = 1; j <= jdim1; j++) {
                        t(jiv+j, 7) = 0.25*gprgm1*(
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,5,3)/qi0(j,kk,1,3) -
                            2.*bci(j,kk,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) -
                            2.*(1.+bci(j,kk-1,2))*qi0(j,kk-1,5,3)/qi0(j,kk-1,1,3) +
                            2.*bci(j,kk-1,2)*q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) +
                            (1.+bck(j,i-1,2))*qk0(j,i-1,5,3)/qk0(j,i-1,1,3) -
                            bck(j,i-1,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) -
                            q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) );
                        t(jiv+j, 8) = 0.25*(
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,2,3) -
                            2.*bci(j,kk,2)*q(j,kk,i-1,2) -
                            2.*(1.+bci(j,kk-1,2))*qi0(j,kk-1,2,3) +
                            2.*bci(j,kk-1,2)*q(j,kk-1,i-1,2) +
                            (1.+bck(j,i-1,2))*qk0(j,i-1,2,3) -
                            bck(j,i-1,2)*q(j,kk,i-1,2) -
                            q(j,kk-1,i-1,2) );
                        t(jiv+j, 9) = 0.25*(
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,3,3) -
                            2.*bci(j,kk,2)*q(j,kk,i-1,3) -
                            2.*(1.+bci(j,kk-1,2))*qi0(j,kk-1,3,3) +
                            2.*bci(j,kk-1,2)*q(j,kk-1,i-1,3) +
                            (1.+bck(j,i-1,2))*qk0(j,i-1,3,3) -
                            bck(j,i-1,2)*q(j,kk,i-1,3) -
                            q(j,kk-1,i-1,3) );
                        t(jiv+j, 10) = 0.25*(
                            2.*(1.+bci(j,kk,2))*qi0(j,kk,4,3) -
                            2.*bci(j,kk,2)*q(j,kk,i-1,4) -
                            2.*(1.+bci(j,kk-1,2))*qi0(j,kk-1,4,3) +
                            2.*bci(j,kk-1,2)*q(j,kk-1,i-1,4) +
                            (1.+bck(j,i-1,2))*qk0(j,i-1,4,3) -
                            bck(j,i-1,2)*q(j,kk,i-1,4) -
                            q(j,kk-1,i-1,4) );
                    }
                } else {
                    // general case for kk = kdim1
                    for (j = 1; j <= jdim1; j++) {
                        t(jiv+j, 7) = 0.25*gprgm1*(
                            (1.+bck(j,i,2))*qk0(j,i,5,3)/qk0(j,i,1,3) -
                            bck(j,i,2)*q(j,kk,i,5)/q(j,kk,i,1) -
                            q(j,kk-1,i,5)/q(j,kk-1,i,1) +
                            (1.+bck(j,i-1,2))*qk0(j,i-1,5,3)/qk0(j,i-1,1,3) -
                            bck(j,i-1,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) -
                            q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) );
                        t(jiv+j, 8) = 0.25*(
                            (1.+bck(j,i,2))*qk0(j,i,2,3) -
                            bck(j,i,2)*q(j,kk,i,2) -
                            q(j,kk-1,i,2) +
                            (1.+bck(j,i-1,2))*qk0(j,i-1,2,3) -
                            bck(j,i-1,2)*q(j,kk,i-1,2) -
                            q(j,kk-1,i-1,2) );
                        t(jiv+j, 9) = 0.25*(
                            (1.+bck(j,i,2))*qk0(j,i,3,3) -
                            bck(j,i,2)*q(j,kk,i,3) -
                            q(j,kk-1,i,3) +
                            (1.+bck(j,i-1,2))*qk0(j,i-1,3,3) -
                            bck(j,i-1,2)*q(j,kk,i-1,3) -
                            q(j,kk-1,i-1,3) );
                        t(jiv+j, 10) = 0.25*(
                            (1.+bck(j,i,2))*qk0(j,i,4,3) -
                            bck(j,i,2)*q(j,kk,i,4) -
                            q(j,kk-1,i,4) +
                            (1.+bck(j,i-1,2))*qk0(j,i-1,4,3) -
                            bck(j,i-1,2)*q(j,kk,i-1,4) -
                            q(j,kk-1,i-1,4) );
                    }
                }
            } // i loop, kk==kdim1


        } else {
            // general case (kk != 1 and kk != kdim1)
            for (i = 1; i <= idim; i++) {
                jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
                if (i == 1) {
                    // replace (i-1) data by boundary data
                    for (j = 1; j <= jdim1; j++) {
                        t(jiv+j, 7)  = 0.25*gprgm1*(
                            q(j,kk+1,i,5)/q(j,kk+1,i,1) -
                            q(j,kk-1,i,5)/q(j,kk-1,i,1) +
                            (1.+bci(j,kk+1,1))*qi0(j,kk+1,5,1)/qi0(j,kk+1,1,1) -
                            bci(j,kk+1,1)*q(j,kk+1,i,5)/q(j,kk+1,i,1) -
                            (1.+bci(j,kk-1,1))*qi0(j,kk-1,5,1)/qi0(j,kk-1,1,1) +
                            bci(j,kk-1,1)*q(j,kk-1,i,5)/q(j,kk-1,i,1) );
                        t(jiv+j, 8)  = 0.25*(
                            q(j,kk+1,i,2) - q(j,kk-1,i,2) +
                            (1.+bci(j,kk+1,1))*qi0(j,kk+1,2,1) -
                            bci(j,kk+1,1)*q(j,kk+1,i,2) -
                            (1.+bci(j,kk-1,1))*qi0(j,kk-1,2,1) +
                            bci(j,kk-1,1)*q(j,kk-1,i,2) );
                        t(jiv+j, 9)  = 0.25*(
                            q(j,kk+1,i,3) - q(j,kk-1,i,3) +
                            (1.+bci(j,kk+1,1))*qi0(j,kk+1,3,1) -
                            bci(j,kk+1,1)*q(j,kk+1,i,3) -
                            (1.+bci(j,kk-1,1))*qi0(j,kk-1,3,1) +
                            bci(j,kk-1,1)*q(j,kk-1,i,3) );
                        t(jiv+j, 10) = 0.25*(
                            q(j,kk+1,i,4) - q(j,kk-1,i,4) +
                            (1.+bci(j,kk+1,1))*qi0(j,kk+1,4,1) -
                            bci(j,kk+1,1)*q(j,kk+1,i,4) -
                            (1.+bci(j,kk-1,1))*qi0(j,kk-1,4,1) +
                            bci(j,kk-1,1)*q(j,kk-1,i,4) );
                    }
                } else if (i == idim) {
                    // replace i data by boundary data
                    for (j = 1; j <= jdim1; j++) {
                        t(jiv+j, 7)  = 0.25*gprgm1*(
                            (1.+bci(j,kk+1,2))*qi0(j,kk+1,5,3)/qi0(j,kk+1,1,3) -
                            bci(j,kk+1,2)*q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
                            (1.+bci(j,kk-1,2))*qi0(j,kk-1,5,3)/qi0(j,kk-1,1,3) +
                            bci(j,kk-1,2)*q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) +
                            q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
                            q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) );
                        t(jiv+j, 8)  = 0.25*(
                            (1.+bci(j,kk+1,2))*qi0(j,kk+1,2,3) -
                            bci(j,kk+1,2)*q(j,kk+1,i-1,2) -
                            (1.+bci(j,kk-1,2))*qi0(j,kk-1,2,3) +
                            bci(j,kk-1,2)*q(j,kk-1,i-1,2) +
                            q(j,kk+1,i-1,2) - q(j,kk-1,i-1,2) );
                        t(jiv+j, 9)  = 0.25*(
                            (1.+bci(j,kk+1,2))*qi0(j,kk+1,3,3) -
                            bci(j,kk+1,2)*q(j,kk+1,i-1,3) -
                            (1.+bci(j,kk-1,2))*qi0(j,kk-1,3,3) +
                            bci(j,kk-1,2)*q(j,kk-1,i-1,3) +
                            q(j,kk+1,i-1,3) - q(j,kk-1,i-1,3) );
                        t(jiv+j, 10) = 0.25*(
                            (1.+bci(j,kk+1,2))*qi0(j,kk+1,4,3) -
                            bci(j,kk+1,2)*q(j,kk+1,i-1,4) -
                            (1.+bci(j,kk-1,2))*qi0(j,kk-1,4,3) +
                            bci(j,kk-1,2)*q(j,kk-1,i-1,4) +
                            q(j,kk+1,i-1,4) - q(j,kk-1,i-1,4) );
                    }
                } else {
                    // general case
                    for (j = 1; j <= jdim1; j++) {
                        t(jiv+j, 7)  = 0.25*gprgm1*(
                            q(j,kk+1,i,5)/q(j,kk+1,i,1) -
                            q(j,kk-1,i,5)/q(j,kk-1,i,1) +
                            q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
                            q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) );
                        t(jiv+j, 8)  = 0.25*(
                            q(j,kk+1,i,2)   - q(j,kk-1,i,2) +
                            q(j,kk+1,i-1,2) - q(j,kk-1,i-1,2) );
                        t(jiv+j, 9)  = 0.25*(
                            q(j,kk+1,i,3)   - q(j,kk-1,i,3) +
                            q(j,kk+1,i-1,3) - q(j,kk-1,i-1,3) );
                        t(jiv+j, 10) = 0.25*(
                            q(j,kk+1,i,4)   - q(j,kk-1,i,4) +
                            q(j,kk+1,i-1,4) - q(j,kk-1,i-1,4) );
                    }
                }
            } // i loop, general kk
        }
    } // kpl loop (zeta gradients)


    // t(24): turbulent viscosity at interfaces (=0 for laminar flow)
    iviscc = ivisc[0];
    if (iviscc > 1) {
        for (kpl = 1; kpl <= npl; kpl++) {
            kk = k + kpl - 1;
            for (i = 1; i <= idim; i++) {
                jiv = (i - 1) * npl * jdim1 + (kpl - 1) * jdim1;
                if (i == 1) {
                    // i=0 interface
                    for (j = 1; j <= jdim1; j++) {
                        t(j+jiv, 24) = bci(j,kk,1)*vi0(j,kk,1,1) +
                            (1. - bci(j,kk,1)) * 0.5 * (vi0(j,kk,1,1) + vist3d(j,kk,1));
                    }
                } else if (i == idim) {
                    // i=idim interface
                    for (j = 1; j <= jdim1; j++) {
                        t(j+jiv, 24) = bci(j,kk,2)*vi0(j,kk,1,3) +
                            (1. - bci(j,kk,2)) * 0.5 * (vi0(j,kk,1,3) + vist3d(j,kk,idim1));
                    }
                } else {
                    // interior interfaces
                    for (j = 1; j <= jdim1; j++) {
                        t(j+jiv, 24) = 0.5 * (vist3d(j,kk,i) + vist3d(j,kk,i-1));
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
        // t(23): ratio of laminar Prandtl number to total Prandtl number
        t(j, 23) = (1.e0 + prtr * t25) / t24;
        // t(14): (mach/re)*viscosity/J
        t(j, 14) = xmre * t24 * t(j, 14) / t(j, 15);
    }

    if (iadv < 0) return;

    // t(11): u at cell interfaces, t(12): v, t(13): w
    // interior interfaces
    for (j = 1; j <= nn; j++) {
        t(j+js-1, 11) = 0.5e0 * (t(j+js-1, 18) + t(j, 18));
        t(j+js-1, 12) = 0.5e0 * (t(j+js-1, 19) + t(j, 19));
        t(j+js-1, 13) = 0.5e0 * (t(j+js-1, 20) + t(j, 20));
    }
    // interfaces at i=0 and i=idim
    for (izz = 1; izz <= jv; izz++) {
        ab = 1.0 + wi0(izz, 15);
        bb = 1.0 - wi0(izz, 15);
        t(izz+n, 11) = 0.5 * (ab * wi0(izz, 18) + bb * t(izz+n-jv, 18));
        t(izz+n, 12) = 0.5 * (ab * wi0(izz, 19) + bb * t(izz+n-jv, 19));
        t(izz+n, 13) = 0.5 * (ab * wi0(izz, 20) + bb * t(izz+n-jv, 20));
        ab = 1.0 + wi0(izz, 5);
        bb = 1.0 - wi0(izz, 5);
        t(izz, 11) = 0.5 * (ab * wi0(izz, 8)  + bb * t(izz, 18));
        t(izz, 12) = 0.5 * (ab * wi0(izz, 9)  + bb * t(izz, 19));
        t(izz, 13) = 0.5 * (ab * wi0(izz, 10) + bb * t(izz, 20));
    }

    // Calculate fluxes: viscous terms at interfaces
    // form ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz
    // excluding contributions in i/xie direction (already accounted for in thin-layer)
    for (j = 1; j <= l0; j++) {
        ux_v = t(j,4)*t(j,28) + t(j,8)*t(j,31);
        uy   = t(j,4)*t(j,29) + t(j,8)*t(j,32);
        uz   = t(j,4)*t(j,30) + t(j,8)*t(j,2);

        vx   = t(j,5)*t(j,28) + t(j,9)*t(j,31);
        vy   = t(j,5)*t(j,29) + t(j,9)*t(j,32);
        vz   = t(j,5)*t(j,30) + t(j,9)*t(j,2);

        wx   = t(j,6)*t(j,28) + t(j,10)*t(j,31);
        wy   = t(j,6)*t(j,29) + t(j,10)*t(j,32);
        wz   = t(j,6)*t(j,30) + t(j,10)*t(j,2);

        tx   = t(j,3)*t(j,28) + t(j,7)*t(j,31);
        ty   = t(j,3)*t(j,29) + t(j,7)*t(j,32);
        tz   = t(j,3)*t(j,30) + t(j,7)*t(j,2);

        // Note lambda = -2/3 mu, lambda + 2 mu = 4/3 mu
        // u-momentum flux
        t(j, 16) = t(j,14) * (  t(j,25) * (4./3.*ux_v - 2./3.*(vy + wz))
                                + t(j,26) * (uy + vx)
                                + t(j,27) * (uz + wx));

        // v-momentum flux
        t(j, 17) = t(j,14) * (  t(j,25) * (vx + uy)
                                + t(j,26) * (4./3.*vy - 2./3.*(ux_v + wz))
                                + t(j,27) * (vz + wy));

        // w-momentum flux
        t(j, 18) = t(j,14) * (  t(j,25) * (wx + uz)
                                + t(j,26) * (wy + vz)
                                + t(j,27) * (4./3.*wz - 2./3.*(vy + ux_v)));

        // energy flux
        t(j, 19) = t(j,14) * (
            t(j,25) * ((4./3.*ux_v - 2./3.*(vy+wz))*t(j,11) +
                        (uy+vx)*t(j,12) + (uz+wx)*t(j,13)) +
            t(j,26) * ((vx+uy)*t(j,11) + (vz+wy)*t(j,13) +
                        (4./3.*vy - 2./3.*(ux_v+wz))*t(j,12)) +
            t(j,27) * ((wx+uz)*t(j,11) + (wy+vz)*t(j,12) +
                       (4./3.*wz - 2./3.*(ux_v+vy))*t(j,13)) +
            t(j,23) * (t(j,25)*tx + t(j,26)*ty + t(j,27)*tz));
    }

    // (-)viscous flux = Fv(i-1/2) - Fv(i+1/2) (momentum eqs.)
    for (j = 1; j <= n; j++) {
        t(j, 2) = -t(j+js-1, 16) + t(j, 16);
        t(j, 3) = -t(j+js-1, 17) + t(j, 17);
        t(j, 4) = -t(j+js-1, 18) + t(j, 18);
        t(j, 5) = -t(j+js-1, 19) + t(j, 19);
    }

    // Calculate residuals
    for (kpl = 1; kpl <= npl; kpl++) {
        kk  = k + kpl - 1;
        for (i = 1; i <= idim1; i++) {
            jiv = (i - 1) * jdim1 * npl + (kpl - 1) * jdim1 + 1;
            for (l = 2; l <= 5; l++) {
                for (j = 1; j <= jdim1; j++) {
                    res(j, kk, i, l) = res(j, kk, i, l) + t(j + jiv - 1, l);
                }
            }
        }
    }

    return;
}

} // namespace ffluxv1_ns
