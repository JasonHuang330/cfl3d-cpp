// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "gfluxv1.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace gfluxv1_ns {

void gfluxv1(int& i, int& npl, int& jdim, int& kdim, int& idim,
             FortranArray4DRef<double> res, FortranArray4DRef<double> q,
             FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
             FortranArray4DRef<double> qi0, FortranArray4DRef<double> sj,
             FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
             FortranArray3DRef<double> vol, FortranArray2DRef<double> t,
             int& nvtq, FortranArray2DRef<double> wj0,
             FortranArray3DRef<double> vist3d, FortranArray4DRef<double> vj0,
             FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
             FortranArray3DRef<double> bci, FortranArray3DRef<double> volj0,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& iadv)
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
    int32_t* ivisc  = cmn_reyue.ivisc;   // 0-based: ivisc[0]=ivisc(1), ivisc[1]=ivisc(2)
    int32_t& isklton = cmn_sklton.isklton;
    int32_t& i2d    = cmn_twod.i2d;
    int32_t& level  = cmn_mgrd.level;
    int32_t& lglobal = cmn_mgrd.lglobal;
    int32_t& i_tauijs = cmn_constit.i_tauijs;

    // Local variables
    float coef_eddy;
    int kdim1, jdim1, idim1;
    int n, l0, kv, ks, nn;
    double xmre, gpr, gm1pr, prtr, gprgm1;
    int l1, ipl, ii, jkv, j, k, izz, m, mm, mb, imin, imax, kmin, kmax, kc;
    int n1;
    double t1, t5, t6, ab, bb, wj05, wj06;
    int iviscc;
    double t24, t25;
    double ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz;


    // Initialize coef_eddy
    coef_eddy = 1.0f;
    if (ivisc[1] >= 70) coef_eddy = 0.0f;
    if (i_tauijs == 1) coef_eddy = 0.0f;

    kdim1 = kdim - 1;
    jdim1 = jdim - 1;
    idim1 = idim - 1;

    // n  : number of cell centers for npl planes
    // l0 : number of cell interfaces for npl planes
    // kv : number of cell centers (and interfaces) on a j=constant plane
    n  = npl * jdim1 * kdim1;
    l0 = npl * kdim1 * jdim;
    kv = npl * kdim1;
    ks = kv + 1;
    nn = n - kv;

    xmre   = (double)xmach / (double)reue;
    gpr    = (double)gamma / (double)pr;
    gm1pr  = (double)gm1 * (double)pr;
    prtr   = (double)pr / (double)prt;
    gprgm1 = gpr / (double)gm1;

    if (isklton > 0 && i == 1 && iadv >= 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "   computing cross-derivative viscous fluxes, J-dir.");
    }

    // ******************************************************************************
    // store selected cell centered data in t array
    // t(22) : density
    // t(18) : u
    // t(19) : v
    // t(20) : w
    // t(16) : c*c/(gm1*pr)  (stored as q(5) = total energy/density)
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

    // ******************************************************************************
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

    // ******************************************************************************
    // t(17)/wj0(m+7)   : (u*u+v*v+w*w)/2
    // t(1)/wj0(m+1)    : 1/density
    // t(16)/wj0(m+6)   : c*c/(gm1*pr)
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

    // ******************************************************************************
    // t(7) : laminar viscosity at cell centers (via sutherland relation)
    double c2b  = (double)cbar / (double)tinf;
    double c2bp = c2b + 1.e0;

    for (izz = 1; izz <= n; izz++) {
        t5       = gm1pr * t(izz, 16);
        t6       = std::sqrt(t5);
        t(izz, 7) = c2bp * t5 * t6 / (c2b + t5);
    }

    // ******************************************************************************
    // t(14) : laminar viscosity values at cell interfaces
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

    // ******************************************************************************
    // t(15) : average jacobian (inverse volume) at cell interface
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

            // ******************************************************************************
            // t(25) - t(27) : components of grad(eta)
            // t1    : grad(eta)/J
            for (k = 1; k <= kdim1; k++) {
                t1           = sj(j, k, ii, 4) * t(k + jkv, 15);
                t(k + jkv, 25) = sj(j, k, ii, 1) * t1;
                t(k + jkv, 26) = sj(j, k, ii, 2) * t1;
                t(k + jkv, 27) = sj(j, k, ii, 3) * t1;
            }
        }
    }


    // ******************************************************************************
    // t(28) - t(30) : Components of grad(zeta)
    for (ipl = 1; ipl <= npl; ipl++) {
        ii = i + ipl - 1;
        for (j = 1; j <= jdim; j++) {
            jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
            if (j == 1) {
                for (k = 1; k <= kdim1; k++) {
                    t(jkv + k, 28) = 0.25 * t(k + jkv, 15) * (
                        sk(j, k,   ii, 1) * sk(j, k,   ii, 4)
                      + sk(j, k,   ii, 1) * sk(j, k,   ii, 4)
                      + sk(j, k+1, ii, 1) * sk(j, k+1, ii, 4)
                      + sk(j, k+1, ii, 1) * sk(j, k+1, ii, 4));
                    t(jkv + k, 29) = 0.25 * t(k + jkv, 15) * (
                        sk(j, k,   ii, 2) * sk(j, k,   ii, 4)
                      + sk(j, k,   ii, 2) * sk(j, k,   ii, 4)
                      + sk(j, k+1, ii, 2) * sk(j, k+1, ii, 4)
                      + sk(j, k+1, ii, 2) * sk(j, k+1, ii, 4));
                    t(jkv + k, 30) = 0.25 * t(k + jkv, 15) * (
                        sk(j, k,   ii, 3) * sk(j, k,   ii, 4)
                      + sk(j, k,   ii, 3) * sk(j, k,   ii, 4)
                      + sk(j, k+1, ii, 3) * sk(j, k+1, ii, 4)
                      + sk(j, k+1, ii, 3) * sk(j, k+1, ii, 4));
                }
            } else if (j == jdim) {
                for (k = 1; k <= kdim1; k++) {
                    t(jkv + k, 28) = 0.25 * t(k + jkv, 15) * (
                        sk(j-1, k,   ii, 1) * sk(j-1, k,   ii, 4)
                      + sk(j-1, k,   ii, 1) * sk(j-1, k,   ii, 4)
                      + sk(j-1, k+1, ii, 1) * sk(j-1, k+1, ii, 4)
                      + sk(j-1, k+1, ii, 1) * sk(j-1, k+1, ii, 4));
                    t(jkv + k, 29) = 0.25 * t(k + jkv, 15) * (
                        sk(j-1, k,   ii, 2) * sk(j-1, k,   ii, 4)
                      + sk(j-1, k,   ii, 2) * sk(j-1, k,   ii, 4)
                      + sk(j-1, k+1, ii, 2) * sk(j-1, k+1, ii, 4)
                      + sk(j-1, k+1, ii, 2) * sk(j-1, k+1, ii, 4));
                    t(jkv + k, 30) = 0.25 * t(k + jkv, 15) * (
                        sk(j-1, k,   ii, 3) * sk(j-1, k,   ii, 4)
                      + sk(j-1, k,   ii, 3) * sk(j-1, k,   ii, 4)
                      + sk(j-1, k+1, ii, 3) * sk(j-1, k+1, ii, 4)
                      + sk(j-1, k+1, ii, 3) * sk(j-1, k+1, ii, 4));
                }
            } else {
                // general case
                for (k = 1; k <= kdim1; k++) {
                    t(jkv + k, 28) = 0.25 * t(k + jkv, 15) * (
                        sk(j,   k,   ii, 1) * sk(j,   k,   ii, 4)
                      + sk(j-1, k,   ii, 1) * sk(j-1, k,   ii, 4)
                      + sk(j,   k+1, ii, 1) * sk(j,   k+1, ii, 4)
                      + sk(j-1, k+1, ii, 1) * sk(j-1, k+1, ii, 4));
                    t(jkv + k, 29) = 0.25 * t(k + jkv, 15) * (
                        sk(j,   k,   ii, 2) * sk(j,   k,   ii, 4)
                      + sk(j-1, k,   ii, 2) * sk(j-1, k,   ii, 4)
                      + sk(j,   k+1, ii, 2) * sk(j,   k+1, ii, 4)
                      + sk(j-1, k+1, ii, 2) * sk(j-1, k+1, ii, 4));
                    t(jkv + k, 30) = 0.25 * t(k + jkv, 15) * (
                        sk(j,   k,   ii, 3) * sk(j,   k,   ii, 4)
                      + sk(j-1, k,   ii, 3) * sk(j-1, k,   ii, 4)
                      + sk(j,   k+1, ii, 3) * sk(j,   k+1, ii, 4)
                      + sk(j-1, k+1, ii, 3) * sk(j-1, k+1, ii, 4));
                }
            }
        }
    }

    // ******************************************************************************
    // t(31), t(32), t(2) : Components of grad(xie)
    for (ipl = 1; ipl <= npl; ipl++) {
        ii = i + ipl - 1;
        for (j = 1; j <= jdim; j++) {
            jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
            if (j == 1) {
                for (k = 1; k <= kdim1; k++) {
                    t(jkv + k, 31) = 0.25 * t(k + jkv, 15) * (
                        si(j, k, ii,   1) * si(j, k, ii,   4)
                      + si(j, k, ii,   1) * si(j, k, ii,   4)
                      + si(j, k, ii+1, 1) * si(j, k, ii+1, 4)
                      + si(j, k, ii+1, 1) * si(j, k, ii+1, 4));
                    t(jkv + k, 32) = 0.25 * t(k + jkv, 15) * (
                        si(j, k, ii,   2) * si(j, k, ii,   4)
                      + si(j, k, ii,   2) * si(j, k, ii,   4)
                      + si(j, k, ii+1, 2) * si(j, k, ii+1, 4)
                      + si(j, k, ii+1, 2) * si(j, k, ii+1, 4));
                    t(jkv + k, 2)  = 0.25 * t(k + jkv, 15) * (
                        si(j, k, ii,   3) * si(j, k, ii,   4)
                      + si(j, k, ii,   3) * si(j, k, ii,   4)
                      + si(j, k, ii+1, 3) * si(j, k, ii+1, 4)
                      + si(j, k, ii+1, 3) * si(j, k, ii+1, 4));
                }
            } else if (j == jdim) {
                for (k = 1; k <= kdim1; k++) {
                    t(jkv + k, 31) = 0.25 * t(k + jkv, 15) * (
                        si(j-1, k, ii,   1) * si(j-1, k, ii,   4)
                      + si(j-1, k, ii,   1) * si(j-1, k, ii,   4)
                      + si(j-1, k, ii+1, 1) * si(j-1, k, ii+1, 4)
                      + si(j-1, k, ii+1, 1) * si(j-1, k, ii+1, 4));
                    t(jkv + k, 32) = 0.25 * t(k + jkv, 15) * (
                        si(j-1, k, ii,   2) * si(j-1, k, ii,   4)
                      + si(j-1, k, ii,   2) * si(j-1, k, ii,   4)
                      + si(j-1, k, ii+1, 2) * si(j-1, k, ii+1, 4)
                      + si(j-1, k, ii+1, 2) * si(j-1, k, ii+1, 4));
                    t(jkv + k, 2)  = 0.25 * t(k + jkv, 15) * (
                        si(j-1, k, ii,   3) * si(j-1, k, ii,   4)
                      + si(j-1, k, ii,   3) * si(j-1, k, ii,   4)
                      + si(j-1, k, ii+1, 3) * si(j-1, k, ii+1, 4)
                      + si(j-1, k, ii+1, 3) * si(j-1, k, ii+1, 4));
                }
            } else {
                // general case
                for (k = 1; k <= kdim1; k++) {
                    t(jkv + k, 31) = 0.25 * t(k + jkv, 15) * (
                        si(j,   k, ii,   1) * si(j,   k, ii,   4)
                      + si(j-1, k, ii,   1) * si(j-1, k, ii,   4)
                      + si(j,   k, ii+1, 1) * si(j,   k, ii+1, 4)
                      + si(j-1, k, ii+1, 1) * si(j-1, k, ii+1, 4));
                    t(jkv + k, 32) = 0.25 * t(k + jkv, 15) * (
                        si(j,   k, ii,   2) * si(j,   k, ii,   4)
                      + si(j-1, k, ii,   2) * si(j-1, k, ii,   4)
                      + si(j,   k, ii+1, 2) * si(j,   k, ii+1, 4)
                      + si(j-1, k, ii+1, 2) * si(j-1, k, ii+1, 4));
                    t(jkv + k, 2)  = 0.25 * t(k + jkv, 15) * (
                        si(j,   k, ii,   3) * si(j,   k, ii,   4)
                      + si(j-1, k, ii,   3) * si(j-1, k, ii,   4)
                      + si(j,   k, ii+1, 3) * si(j,   k, ii+1, 4)
                      + si(j-1, k, ii+1, 3) * si(j-1, k, ii+1, 4));
                }
            }
        }
    }


    // ******************************************************************************
    // gradients at cell interfaces (needed for full ns terms)
    // t(3) : d( c*c/(gm1*pr) )/d(zeta)
    // t(4) : d(u)/d(zeta)
    // t(5) : d(v)/d(zeta)
    // t(6) : d(w)/d(zeta)
    for (ipl = 1; ipl <= npl; ipl++) {
        ii = i + ipl - 1;
        for (j = 1; j <= jdim; j++) {
            jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;

            if (j == 1) {
                // replace (j-1) data by boundary data
                for (k = 1; k <= kdim1; k++) {
                    if (k == 1) {
                        // replace (k-1) data by boundary data
                        t(jkv + k, 3) = 0.25 * gprgm1 * (
                            q(j, k+1, ii, 5) / q(j, k+1, ii, 1) -
                            (1.0 + bck(j, ii, 1)) * qk0(j, ii, 5, 1) / qk0(j, ii, 1, 1) +
                            bck(j, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) +
                            2.0 * (1.0 + bcj(k+1, ii, 1)) * qj0(k+1, ii, 5, 1) / qj0(k+1, ii, 1, 1) -
                            2.0 * bcj(k+1, ii, 1) * q(j, k+1, ii, 5) / q(j, k+1, ii, 1) -
                            2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 5, 1) / qj0(k, ii, 1, 1) +
                            2.0 * bcj(k, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1));
                        t(jkv + k, 4) = 0.25 * (
                            q(j, k+1, ii, 2) -
                            (1.0 + bck(j, ii, 1)) * qk0(j, ii, 2, 1) +
                            bck(j, ii, 1) * q(j, k, ii, 2) +
                            2.0 * (1.0 + bcj(k+1, ii, 1)) * qj0(k+1, ii, 2, 1) -
                            2.0 * bcj(k+1, ii, 1) * q(j, k+1, ii, 2) -
                            2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 2, 1) +
                            2.0 * bcj(k, ii, 1) * q(j, k, ii, 2));
                        t(jkv + k, 5) = 0.25 * (
                            q(j, k+1, ii, 3) -
                            (1.0 + bck(j, ii, 1)) * qk0(j, ii, 3, 1) +
                            bck(j, ii, 1) * q(j, k, ii, 3) +
                            2.0 * (1.0 + bcj(k+1, ii, 1)) * qj0(k+1, ii, 3, 1) -
                            2.0 * bcj(k+1, ii, 1) * q(j, k+1, ii, 3) -
                            2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 3, 1) +
                            2.0 * bcj(k, ii, 1) * q(j, k, ii, 3));
                        t(jkv + k, 6) = 0.25 * (
                            q(j, k+1, ii, 4) -
                            (1.0 + bck(j, ii, 1)) * qk0(j, ii, 4, 1) +
                            bck(j, ii, 1) * q(j, k, ii, 4) +
                            2.0 * (1.0 + bcj(k+1, ii, 1)) * qj0(k+1, ii, 4, 1) -
                            2.0 * bcj(k+1, ii, 1) * q(j, k+1, ii, 4) -
                            2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 4, 1) +
                            2.0 * bcj(k, ii, 1) * q(j, k, ii, 4));
                    } else if (k == kdim1) {
                        // replace (k+1) data by boundary data
                        t(jkv + k, 3) = 0.25 * gprgm1 * (
                            (1.0 + bck(j, ii, 2)) * qk0(j, ii, 5, 3) / qk0(j, ii, 1, 3) -
                            bck(j, ii, 2) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                            q(j, k-1, ii, 5) / q(j, k-1, ii, 1) +
                            2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 5, 1) / qj0(k, ii, 1, 1) -
                            2.0 * bcj(k, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                            2.0 * (1.0 + bcj(k-1, ii, 1)) * qj0(k-1, ii, 5, 1) / qj0(k-1, ii, 1, 1) +
                            2.0 * bcj(k-1, ii, 1) * q(j, k-1, ii, 5) / q(j, k-1, ii, 1));
                        t(jkv + k, 4) = 0.25 * (
                            (1.0 + bck(j, ii, 2)) * qk0(j, ii, 2, 3) -
                            bck(j, ii, 2) * q(j, k, ii, 2) -
                            q(j, k-1, ii, 2) +
                            2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 2, 1) -
                            2.0 * bcj(k, ii, 1) * q(j, k, ii, 2) -
                            2.0 * (1.0 + bcj(k-1, ii, 1)) * qj0(k-1, ii, 2, 1) +
                            2.0 * bcj(k-1, ii, 1) * q(j, k-1, ii, 2));
                        t(jkv + k, 5) = 0.25 * (
                            (1.0 + bck(j, ii, 2)) * qk0(j, ii, 3, 3) -
                            bck(j, ii, 2) * q(j, k, ii, 3) -
                            q(j, k-1, ii, 3) +
                            2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 3, 1) -
                            2.0 * bcj(k, ii, 1) * q(j, k, ii, 3) -
                            2.0 * (1.0 + bcj(k-1, ii, 1)) * qj0(k-1, ii, 3, 1) +
                            2.0 * bcj(k-1, ii, 1) * q(j, k-1, ii, 3));
                        t(jkv + k, 6) = 0.25 * (
                            (1.0 + bck(j, ii, 2)) * qk0(j, ii, 4, 3) -
                            bck(j, ii, 2) * q(j, k, ii, 4) -
                            q(j, k-1, ii, 4) +
                            2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 4, 1) -
                            2.0 * bcj(k, ii, 1) * q(j, k, ii, 4) -
                            2.0 * (1.0 + bcj(k-1, ii, 1)) * qj0(k-1, ii, 4, 1) +
                            2.0 * bcj(k-1, ii, 1) * q(j, k-1, ii, 4));
                    } else {
                        // general case for j = 1
                        t(jkv + k, 3) = 0.25 * gprgm1 * (
                            q(j, k+1, ii, 5) / q(j, k+1, ii, 1) -
                            q(j, k-1, ii, 5) / q(j, k-1, ii, 1) +
                            (1.0 + bcj(k+1, ii, 1)) * qj0(k+1, ii, 5, 1) / qj0(k+1, ii, 1, 1) -
                            bcj(k+1, ii, 1) * q(j, k+1, ii, 5) / q(j, k+1, ii, 1) -
                            (1.0 + bcj(k-1, ii, 1)) * qj0(k-1, ii, 5, 1) / qj0(k-1, ii, 1, 1) +
                            bcj(k-1, ii, 1) * q(j, k-1, ii, 5) / q(j, k-1, ii, 1));
                        t(jkv + k, 4) = 0.25 * (
                            q(j, k+1, ii, 2) - q(j, k-1, ii, 2) +
                            (1.0 + bcj(k+1, ii, 1)) * qj0(k+1, ii, 2, 1) -
                            bcj(k+1, ii, 1) * q(j, k+1, ii, 2) -
                            (1.0 + bcj(k-1, ii, 1)) * qj0(k-1, ii, 2, 1) +
                            bcj(k-1, ii, 1) * q(j, k-1, ii, 2));
                        t(jkv + k, 5) = 0.25 * (
                            q(j, k+1, ii, 3) - q(j, k-1, ii, 3) +
                            (1.0 + bcj(k+1, ii, 1)) * qj0(k+1, ii, 3, 1) -
                            bcj(k+1, ii, 1) * q(j, k+1, ii, 3) -
                            (1.0 + bcj(k-1, ii, 1)) * qj0(k-1, ii, 3, 1) +
                            bcj(k-1, ii, 1) * q(j, k-1, ii, 3));
                        t(jkv + k, 6) = 0.25 * (
                            q(j, k+1, ii, 4) - q(j, k-1, ii, 4) +
                            (1.0 + bcj(k+1, ii, 1)) * qj0(k+1, ii, 4, 1) -
                            bcj(k+1, ii, 1) * q(j, k+1, ii, 4) -
                            (1.0 + bcj(k-1, ii, 1)) * qj0(k-1, ii, 4, 1) +
                            bcj(k-1, ii, 1) * q(j, k-1, ii, 4));
                    }
                } // end k loop for j==1


            } else if (j == jdim) {
                // replace j data by boundary data
                for (k = 1; k <= kdim1; k++) {
                    if (k == 1) {
                        // replace (k-1) data by boundary data
                        t(jkv + k, 3) = 0.25 * gprgm1 * (
                            2.0 * (1.0 + bcj(k+1, ii, 2)) * qj0(k+1, ii, 5, 3) / qj0(k+1, ii, 1, 3) -
                            2.0 * bcj(k+1, ii, 2) * q(j-1, k+1, ii, 5) / q(j-1, k+1, ii, 1) -
                            2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 5, 3) / qj0(k, ii, 1, 3) +
                            2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1) +
                            q(j-1, k+1, ii, 5) / q(j-1, k+1, ii, 1) -
                            (1.0 + bck(j-1, ii, 1)) * qk0(j-1, ii, 5, 1) / qk0(j-1, ii, 1, 1) +
                            bck(j-1, ii, 1) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1));
                        t(jkv + k, 4) = 0.25 * (
                            2.0 * (1.0 + bcj(k+1, ii, 2)) * qj0(k+1, ii, 2, 3) -
                            2.0 * bcj(k+1, ii, 2) * q(j-1, k+1, ii, 2) -
                            2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 2, 3) +
                            2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 2) +
                            q(j-1, k+1, ii, 2) -
                            (1.0 + bck(j-1, ii, 1)) * qk0(j-1, ii, 2, 1) +
                            bck(j-1, ii, 1) * q(j-1, k, ii, 2));
                        t(jkv + k, 5) = 0.25 * (
                            2.0 * (1.0 + bcj(k+1, ii, 2)) * qj0(k+1, ii, 3, 3) -
                            2.0 * bcj(k+1, ii, 2) * q(j-1, k+1, ii, 3) -
                            2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 3, 3) +
                            2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 3) +
                            q(j-1, k+1, ii, 3) -
                            (1.0 + bck(j-1, ii, 1)) * qk0(j-1, ii, 3, 1) +
                            bck(j-1, ii, 1) * q(j-1, k, ii, 3));
                        t(jkv + k, 6) = 0.25 * (
                            2.0 * (1.0 + bcj(k+1, ii, 2)) * qj0(k+1, ii, 4, 3) -
                            2.0 * bcj(k+1, ii, 2) * q(j-1, k+1, ii, 4) -
                            2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 4, 3) +
                            2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 4) +
                            q(j-1, k+1, ii, 4) -
                            (1.0 + bck(j-1, ii, 1)) * qk0(j-1, ii, 4, 1) +
                            bck(j-1, ii, 1) * q(j-1, k, ii, 4));
                    } else if (k == kdim1) {
                        // replace (k+1) data by boundary data
                        t(jkv + k, 3) = 0.25 * gprgm1 * (
                            2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 5, 3) / qj0(k, ii, 1, 3) -
                            2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1) -
                            2.0 * (1.0 + bcj(k-1, ii, 2)) * qj0(k-1, ii, 5, 3) / qj0(k-1, ii, 1, 3) +
                            2.0 * bcj(k-1, ii, 2) * q(j-1, k-1, ii, 5) / q(j-1, k-1, ii, 1) +
                            (1.0 + bck(j-1, ii, 2)) * qk0(j-1, ii, 5, 3) / qk0(j-1, ii, 1, 3) -
                            bck(j-1, ii, 2) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1) -
                            q(j-1, k-1, ii, 5) / q(j-1, k-1, ii, 1));
                        t(jkv + k, 4) = 0.25 * (
                            2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 2, 3) -
                            2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 2) -
                            2.0 * (1.0 + bcj(k-1, ii, 2)) * qj0(k-1, ii, 2, 3) +
                            2.0 * bcj(k-1, ii, 2) * q(j-1, k-1, ii, 2) +
                            (1.0 + bck(j-1, ii, 2)) * qk0(j-1, ii, 2, 3) -
                            bck(j-1, ii, 2) * q(j-1, k, ii, 2) -
                            q(j-1, k-1, ii, 2));
                        t(jkv + k, 5) = 0.25 * (
                            2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 3, 3) -
                            2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 3) -
                            2.0 * (1.0 + bcj(k-1, ii, 2)) * qj0(k-1, ii, 3, 3) +
                            2.0 * bcj(k-1, ii, 2) * q(j-1, k-1, ii, 3) +
                            (1.0 + bck(j-1, ii, 2)) * qk0(j-1, ii, 3, 3) -
                            bck(j-1, ii, 2) * q(j-1, k, ii, 3) -
                            q(j-1, k-1, ii, 3));
                        t(jkv + k, 6) = 0.25 * (
                            2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 4, 3) -
                            2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 4) -
                            2.0 * (1.0 + bcj(k-1, ii, 2)) * qj0(k-1, ii, 4, 3) +
                            2.0 * bcj(k-1, ii, 2) * q(j-1, k-1, ii, 4) +
                            (1.0 + bck(j-1, ii, 2)) * qk0(j-1, ii, 4, 3) -
                            bck(j-1, ii, 2) * q(j-1, k, ii, 4) -
                            q(j-1, k-1, ii, 4));
                    } else {
                        // general case for j = jdim
                        t(jkv + k, 3) = 0.25 * gprgm1 * (
                            (1.0 + bcj(k+1, ii, 2)) * qj0(k+1, ii, 5, 3) / qj0(k+1, ii, 1, 3) -
                            bcj(k+1, ii, 2) * q(j-1, k+1, ii, 5) / q(j-1, k+1, ii, 1) -
                            (1.0 + bcj(k-1, ii, 2)) * qj0(k-1, ii, 5, 3) / qj0(k-1, ii, 1, 3) +
                            bcj(k-1, ii, 2) * q(j-1, k-1, ii, 5) / q(j-1, k-1, ii, 1) +
                            q(j-1, k+1, ii, 5) / q(j-1, k+1, ii, 1) -
                            q(j-1, k-1, ii, 5) / q(j-1, k-1, ii, 1));
                        t(jkv + k, 4) = 0.25 * (
                            (1.0 + bcj(k+1, ii, 2)) * qj0(k+1, ii, 2, 3) -
                            bcj(k+1, ii, 2) * q(j-1, k+1, ii, 2) -
                            (1.0 + bcj(k-1, ii, 2)) * qj0(k-1, ii, 2, 3) +
                            bcj(k-1, ii, 2) * q(j-1, k-1, ii, 2) +
                            q(j-1, k+1, ii, 2) - q(j-1, k-1, ii, 2));
                        t(jkv + k, 5) = 0.25 * (
                            (1.0 + bcj(k+1, ii, 2)) * qj0(k+1, ii, 3, 3) -
                            bcj(k+1, ii, 2) * q(j-1, k+1, ii, 3) -
                            (1.0 + bcj(k-1, ii, 2)) * qj0(k-1, ii, 3, 3) +
                            bcj(k-1, ii, 2) * q(j-1, k-1, ii, 3) +
                            q(j-1, k+1, ii, 3) - q(j-1, k-1, ii, 3));
                        t(jkv + k, 6) = 0.25 * (
                            (1.0 + bcj(k+1, ii, 2)) * qj0(k+1, ii, 4, 3) -
                            bcj(k+1, ii, 2) * q(j-1, k+1, ii, 4) -
                            (1.0 + bcj(k-1, ii, 2)) * qj0(k-1, ii, 4, 3) +
                            bcj(k-1, ii, 2) * q(j-1, k-1, ii, 4) +
                            q(j-1, k+1, ii, 4) - q(j-1, k-1, ii, 4));
                    }
                } // end k loop for j==jdim


            } else {
                // general case for j between 1 and jdim
                for (k = 1; k <= kdim1; k++) {
                    if (k == 1) {
                        // replace (k-1) data by boundary data
                        t(jkv + k, 3) = 0.25 * gprgm1 * (
                            q(j, k+1, ii, 5) / q(j, k+1, ii, 1) -
                            (1.0 + bck(j, ii, 1)) * qk0(j, ii, 5, 1) / qk0(j, ii, 1, 1) +
                            bck(j, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) +
                            q(j-1, k+1, ii, 5) / q(j-1, k+1, ii, 1) -
                            (1.0 + bck(j-1, ii, 1)) * qk0(j-1, ii, 5, 1) / qk0(j-1, ii, 1, 1) +
                            bck(j-1, ii, 1) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1));
                        t(jkv + k, 4) = 0.25 * (
                            q(j, k+1, ii, 2) -
                            (1.0 + bck(j, ii, 1)) * qk0(j, ii, 2, 1) +
                            bck(j, ii, 1) * q(j, k, ii, 2) +
                            q(j-1, k+1, ii, 2) -
                            (1.0 + bck(j-1, ii, 1)) * qk0(j-1, ii, 2, 1) +
                            bck(j-1, ii, 1) * q(j-1, k, ii, 2));
                        t(jkv + k, 5) = 0.25 * (
                            q(j, k+1, ii, 3) -
                            (1.0 + bck(j, ii, 1)) * qk0(j, ii, 3, 1) +
                            bck(j, ii, 1) * q(j, k, ii, 3) +
                            q(j-1, k+1, ii, 3) -
                            (1.0 + bck(j-1, ii, 1)) * qk0(j-1, ii, 3, 1) +
                            bck(j-1, ii, 1) * q(j-1, k, ii, 3));
                        t(jkv + k, 6) = 0.25 * (
                            q(j, k+1, ii, 4) -
                            (1.0 + bck(j, ii, 1)) * qk0(j, ii, 4, 1) +
                            bck(j, ii, 1) * q(j, k, ii, 4) +
                            q(j-1, k+1, ii, 4) -
                            (1.0 + bck(j-1, ii, 1)) * qk0(j-1, ii, 4, 1) +
                            bck(j-1, ii, 1) * q(j-1, k, ii, 4));
                    } else if (k == kdim1) {
                        // replace (k+1) data by boundary data
                        t(jkv + k, 3) = 0.25 * gprgm1 * (
                            (1.0 + bck(j, ii, 2)) * qk0(j, ii, 5, 3) / qk0(j, ii, 1, 3) -
                            bck(j, ii, 2) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                            q(j, k-1, ii, 5) / q(j, k-1, ii, 1) +
                            (1.0 + bck(j-1, ii, 2)) * qk0(j-1, ii, 5, 3) / qk0(j-1, ii, 1, 3) -
                            bck(j-1, ii, 2) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1) -
                            q(j-1, k-1, ii, 5) / q(j-1, k-1, ii, 1));
                        t(jkv + k, 4) = 0.25 * (
                            (1.0 + bck(j, ii, 2)) * qk0(j, ii, 2, 3) -
                            bck(j, ii, 2) * q(j, k, ii, 2) -
                            q(j, k-1, ii, 2) +
                            (1.0 + bck(j-1, ii, 2)) * qk0(j-1, ii, 2, 3) -
                            bck(j-1, ii, 2) * q(j-1, k, ii, 2) -
                            q(j-1, k-1, ii, 2));
                        t(jkv + k, 5) = 0.25 * (
                            (1.0 + bck(j, ii, 2)) * qk0(j, ii, 3, 3) -
                            bck(j, ii, 2) * q(j, k, ii, 3) -
                            q(j, k-1, ii, 3) +
                            (1.0 + bck(j-1, ii, 2)) * qk0(j-1, ii, 3, 3) -
                            bck(j-1, ii, 2) * q(j-1, k, ii, 3) -
                            q(j-1, k-1, ii, 3));
                        t(jkv + k, 6) = 0.25 * (
                            (1.0 + bck(j, ii, 2)) * qk0(j, ii, 4, 3) -
                            bck(j, ii, 2) * q(j, k, ii, 4) -
                            q(j, k-1, ii, 4) +
                            (1.0 + bck(j-1, ii, 2)) * qk0(j-1, ii, 4, 3) -
                            bck(j-1, ii, 2) * q(j-1, k, ii, 4) -
                            q(j-1, k-1, ii, 4));
                    } else {
                        // general case
                        t(jkv + k, 3) = 0.25 * gprgm1 * (
                            q(j,   k+1, ii, 5) / q(j,   k+1, ii, 1) -
                            q(j,   k-1, ii, 5) / q(j,   k-1, ii, 1) +
                            q(j-1, k+1, ii, 5) / q(j-1, k+1, ii, 1) -
                            q(j-1, k-1, ii, 5) / q(j-1, k-1, ii, 1));
                        t(jkv + k, 4) = 0.25 * (
                            q(j,   k+1, ii, 2) - q(j,   k-1, ii, 2) +
                            q(j-1, k+1, ii, 2) - q(j-1, k-1, ii, 2));
                        t(jkv + k, 5) = 0.25 * (
                            q(j,   k+1, ii, 3) - q(j,   k-1, ii, 3) +
                            q(j-1, k+1, ii, 3) - q(j-1, k-1, ii, 3));
                        t(jkv + k, 6) = 0.25 * (
                            q(j,   k+1, ii, 4) - q(j,   k-1, ii, 4) +
                            q(j-1, k+1, ii, 4) - q(j-1, k-1, ii, 4));
                    }
                } // end k loop general j
            }
        } // end j loop
    } // end ipl loop (do 230)


    // ******************************************************************************
    // gradients at cell interfaces (needed for full ns terms)
    // t(7) : d( c*c/(gm1*pr) )/d(xie)
    // t(8) : d(u)/d(xie)
    // t(9) : d(v)/d(xie)
    // t(10): d(w)/d(xie)
    if (i2d == 0 && idim > 2) {
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            if (ii == 1) {
                // replace (ii-1) data by boundary data
                for (j = 1; j <= jdim; j++) {
                    jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
                    if (j == 1) {
                        // replace (j-1) data by boundary data
                        for (k = 1; k <= kdim1; k++) {
                            t(jkv + k, 7) = 0.25 * gprgm1 * (
                                q(j, k, ii+1, 5) / q(j, k, ii+1, 1) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 5, 1) / qi0(j, k, 1, 1) +
                                bci(j, k, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) +
                                2.0 * (1.0 + bcj(k, ii+1, 1)) * qj0(k, ii+1, 5, 1) / qj0(k, ii+1, 1, 1) -
                                2.0 * bcj(k, ii+1, 1) * q(j, k, ii+1, 5) / q(j, k, ii+1, 1) -
                                2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 5, 1) / qj0(k, ii, 1, 1) +
                                2.0 * bcj(k, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1));
                            t(jkv + k, 8) = 0.25 * (
                                q(j, k, ii+1, 2) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 2, 1) +
                                bci(j, k, 1) * q(j, k, ii, 2) +
                                2.0 * (1.0 + bcj(k, ii+1, 1)) * qj0(k, ii+1, 2, 1) -
                                2.0 * bcj(k, ii+1, 1) * q(j, k, ii+1, 2) -
                                2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 2, 1) +
                                2.0 * bcj(k, ii, 1) * q(j, k, ii, 2));
                            t(jkv + k, 9) = 0.25 * (
                                q(j, k, ii+1, 3) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 3, 1) +
                                bci(j, k, 1) * q(j, k, ii, 3) +
                                2.0 * (1.0 + bcj(k, ii+1, 1)) * qj0(k, ii+1, 3, 1) -
                                2.0 * bcj(k, ii+1, 1) * q(j, k, ii+1, 3) -
                                2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 3, 1) +
                                2.0 * bcj(k, ii, 1) * q(j, k, ii, 3));
                            t(jkv + k, 10) = 0.25 * (
                                q(j, k, ii+1, 4) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 4, 1) +
                                bci(j, k, 1) * q(j, k, ii, 4) +
                                2.0 * (1.0 + bcj(k, ii+1, 1)) * qj0(k, ii+1, 4, 1) -
                                2.0 * bcj(k, ii+1, 1) * q(j, k, ii+1, 4) -
                                2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 4, 1) +
                                2.0 * bcj(k, ii, 1) * q(j, k, ii, 4));
                        }
                    } else if (j == jdim) {
                        // replace j data by boundary data
                        for (k = 1; k <= kdim1; k++) {
                            t(jkv + k, 7) = 0.25 * gprgm1 * (
                                2.0 * (1.0 + bcj(k, ii+1, 2)) * qj0(k, ii+1, 5, 3) / qj0(k, ii+1, 1, 3) -
                                2.0 * bcj(k, ii+1, 2) * q(j-1, k, ii+1, 5) / q(j-1, k, ii+1, 1) -
                                2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 5, 3) / qj0(k, ii, 1, 3) +
                                2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1) +
                                q(j-1, k, ii+1, 5) / q(j-1, k, ii+1, 1) -
                                (1.0 + bci(j-1, k, 1)) * qi0(j-1, k, 5, 1) / qi0(j-1, k, 1, 1) +
                                bci(j-1, k, 1) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1));
                            t(jkv + k, 8) = 0.25 * (
                                2.0 * (1.0 + bcj(k, ii+1, 2)) * qj0(k, ii+1, 2, 3) -
                                2.0 * bcj(k, ii+1, 2) * q(j-1, k, ii+1, 2) -
                                2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 2, 3) +
                                2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 2) +
                                q(j-1, k, ii+1, 2) -
                                (1.0 + bci(j-1, k, 1)) * qi0(j-1, k, 2, 1) +
                                bci(j-1, k, 1) * q(j-1, k, ii, 2));
                            t(jkv + k, 9) = 0.25 * (
                                2.0 * (1.0 + bcj(k, ii+1, 2)) * qj0(k, ii+1, 3, 3) -
                                2.0 * bcj(k, ii+1, 2) * q(j-1, k, ii+1, 3) -
                                2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 3, 3) +
                                2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 3) +
                                q(j-1, k, ii+1, 3) -
                                (1.0 + bci(j-1, k, 1)) * qi0(j-1, k, 3, 1) +
                                bci(j-1, k, 1) * q(j-1, k, ii, 3));
                            t(jkv + k, 10) = 0.25 * (
                                2.0 * (1.0 + bcj(k, ii+1, 2)) * qj0(k, ii+1, 4, 3) -
                                2.0 * bcj(k, ii+1, 2) * q(j-1, k, ii+1, 4) -
                                2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 4, 3) +
                                2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 4) +
                                q(j-1, k, ii+1, 4) -
                                (1.0 + bci(j-1, k, 1)) * qi0(j-1, k, 4, 1) +
                                bci(j-1, k, 1) * q(j-1, k, ii, 4));
                        }
                    } else {
                        // general case for ii = 1
                        for (k = 1; k <= kdim1; k++) {
                            t(jkv + k, 7) = 0.25 * gprgm1 * (
                                q(j, k, ii+1, 5) / q(j, k, ii+1, 1) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 5, 1) / qi0(j, k, 1, 1) +
                                bci(j, k, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) +
                                q(j-1, k, ii+1, 5) / q(j-1, k, ii+1, 1) -
                                (1.0 + bci(j-1, k, 1)) * qi0(j-1, k, 5, 1) / qi0(j-1, k, 1, 1) +
                                bci(j-1, k, 1) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1));
                            t(jkv + k, 8) = 0.25 * (
                                q(j, k, ii+1, 2) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 2, 1) +
                                bci(j, k, 1) * q(j, k, ii, 2) +
                                q(j-1, k, ii+1, 2) -
                                (1.0 + bci(j-1, k, 1)) * qi0(j-1, k, 2, 1) +
                                bci(j-1, k, 1) * q(j-1, k, ii, 2));
                            t(jkv + k, 9) = 0.25 * (
                                q(j, k, ii+1, 3) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 3, 1) +
                                bci(j, k, 1) * q(j, k, ii, 3) +
                                q(j-1, k, ii+1, 3) -
                                (1.0 + bci(j-1, k, 1)) * qi0(j-1, k, 3, 1) +
                                bci(j-1, k, 1) * q(j-1, k, ii, 3));
                            t(jkv + k, 10) = 0.25 * (
                                q(j, k, ii+1, 4) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 4, 1) +
                                bci(j, k, 1) * q(j, k, ii, 4) +
                                q(j-1, k, ii+1, 4) -
                                (1.0 + bci(j-1, k, 1)) * qi0(j-1, k, 4, 1) +
                                bci(j-1, k, 1) * q(j-1, k, ii, 4));
                        }
                    }
                } // end j loop for ii==1


            } else if (ii == idim1) {
                // replace (ii+1) data by boundary data
                for (j = 1; j <= jdim; j++) {
                    jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
                    if (j == 1) {
                        // replace (j-1) data by boundary data
                        for (k = 1; k <= kdim1; k++) {
                            t(jkv + k, 7) = 0.25 * gprgm1 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 5, 3) / qi0(j, k, 1, 3) -
                                bci(j, k, 2) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                                q(j, k, ii-1, 5) / q(j, k, ii-1, 1) +
                                2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 5, 1) / qj0(k, ii, 1, 1) -
                                2.0 * bcj(k, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                                2.0 * (1.0 + bcj(k, ii-1, 1)) * qj0(k, ii-1, 5, 1) / qj0(k, ii-1, 1, 1) +
                                2.0 * bcj(k, ii-1, 1) * q(j, k, ii-1, 5) / q(j, k, ii-1, 1));
                            t(jkv + k, 8) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 2, 3) -
                                bci(j, k, 2) * q(j, k, ii, 2) -
                                q(j, k, ii-1, 2) +
                                2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 2, 1) -
                                2.0 * bcj(k, ii, 1) * q(j, k, ii, 2) -
                                2.0 * (1.0 + bcj(k, ii-1, 1)) * qj0(k, ii-1, 2, 1) +
                                2.0 * bcj(k, ii-1, 1) * q(j, k, ii-1, 2));
                            t(jkv + k, 9) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 3, 3) -
                                bci(j, k, 2) * q(j, k, ii, 3) -
                                q(j, k, ii-1, 3) +
                                2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 3, 1) -
                                2.0 * bcj(k, ii, 1) * q(j, k, ii, 3) -
                                2.0 * (1.0 + bcj(k, ii-1, 1)) * qj0(k, ii-1, 3, 1) +
                                2.0 * bcj(k, ii-1, 1) * q(j, k, ii-1, 3));
                            t(jkv + k, 10) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 4, 3) -
                                bci(j, k, 2) * q(j, k, ii, 4) -
                                q(j, k, ii-1, 4) +
                                2.0 * (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 4, 1) -
                                2.0 * bcj(k, ii, 1) * q(j, k, ii, 4) -
                                2.0 * (1.0 + bcj(k, ii-1, 1)) * qj0(k, ii-1, 4, 1) +
                                2.0 * bcj(k, ii-1, 1) * q(j, k, ii-1, 4));
                        }
                    } else if (j == jdim) {
                        // replace j data by boundary data
                        for (k = 1; k <= kdim1; k++) {
                            t(jkv + k, 7) = 0.25 * gprgm1 * (
                                2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 5, 3) / qj0(k, ii, 1, 3) -
                                2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1) -
                                2.0 * (1.0 + bcj(k, ii-1, 2)) * qj0(k, ii-1, 5, 3) / qj0(k, ii-1, 1, 3) +
                                2.0 * bcj(k, ii-1, 2) * q(j-1, k, ii-1, 5) / q(j-1, k, ii-1, 1) +
                                (1.0 + bci(j-1, k, 2)) * qi0(j-1, k, 5, 3) / qi0(j-1, k, 1, 3) -
                                bci(j-1, k, 2) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1) -
                                q(j-1, k, ii-1, 5) / q(j-1, k, ii-1, 1));
                            t(jkv + k, 8) = 0.25 * (
                                2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 2, 3) -
                                2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 2) -
                                2.0 * (1.0 + bcj(k, ii-1, 2)) * qj0(k, ii-1, 2, 3) +
                                2.0 * bcj(k, ii-1, 2) * q(j-1, k, ii-1, 2) +
                                (1.0 + bci(j-1, k, 2)) * qi0(j-1, k, 2, 3) -
                                bci(j-1, k, 2) * q(j-1, k, ii, 2) -
                                q(j-1, k, ii-1, 4));  // note: Fortran has ii-1,4 here
                            t(jkv + k, 9) = 0.25 * (
                                2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 3, 3) -
                                2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 3) -
                                2.0 * (1.0 + bcj(k, ii-1, 2)) * qj0(k, ii-1, 3, 3) +
                                2.0 * bcj(k, ii-1, 2) * q(j-1, k, ii-1, 3) +
                                (1.0 + bci(j-1, k, 2)) * qi0(j-1, k, 3, 3) -
                                bci(j-1, k, 2) * q(j-1, k, ii, 3) -
                                q(j-1, k, ii-1, 3));
                            t(jkv + k, 10) = 0.25 * (
                                2.0 * (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 4, 3) -
                                2.0 * bcj(k, ii, 2) * q(j-1, k, ii, 4) -
                                2.0 * (1.0 + bcj(k, ii-1, 2)) * qj0(k, ii-1, 4, 3) +
                                2.0 * bcj(k, ii-1, 2) * q(j-1, k, ii-1, 4) +
                                (1.0 + bci(j-1, k, 2)) * qi0(j-1, k, 4, 3) -
                                bci(j-1, k, 2) * q(j-1, k, ii, 4) -
                                q(j-1, k, ii-1, 4));
                        }
                    } else {
                        // general case for ii = idim1
                        for (k = 1; k <= kdim1; k++) {
                            t(jkv + k, 7) = 0.25 * gprgm1 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 5, 3) / qi0(j, k, 1, 3) -
                                bci(j, k, 2) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                                q(j, k, ii-1, 5) / q(j, k, ii-1, 1) +
                                (1.0 + bci(j-1, k, 2)) * qi0(j-1, k, 5, 3) / qi0(j-1, k, 1, 3) -
                                bci(j-1, k, 2) * q(j-1, k, ii, 5) / q(j-1, k, ii, 1) -
                                q(j-1, k, ii-1, 5) / q(j-1, k, ii-1, 1));
                            t(jkv + k, 8) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 2, 3) -
                                bci(j, k, 2) * q(j, k, ii, 2) -
                                q(j, k, ii-1, 2) +
                                (1.0 + bci(j-1, k, 2)) * qi0(j-1, k, 2, 3) -
                                bci(j-1, k, 2) * q(j-1, k, ii, 2) -
                                q(j-1, k, ii-1, 2));
                            t(jkv + k, 9) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 3, 3) -
                                bci(j, k, 2) * q(j, k, ii, 3) -
                                q(j, k, ii-1, 3) +
                                (1.0 + bci(j-1, k, 2)) * qi0(j-1, k, 3, 3) -
                                bci(j-1, k, 2) * q(j-1, k, ii, 3) -
                                q(j-1, k, ii-1, 3));
                            t(jkv + k, 10) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 4, 3) -
                                bci(j, k, 2) * q(j, k, ii, 4) -
                                q(j, k, ii-1, 4) +
                                (1.0 + bci(j-1, k, 2)) * qi0(j-1, k, 4, 3) -
                                bci(j-1, k, 2) * q(j-1, k, ii, 4) -
                                q(j-1, k, ii-1, 4));
                        }
                    }
                } // end j loop for ii==idim1


            } else {
                // general case (ii .ne. 1 and ii .ne. idim1)
                for (j = 1; j <= jdim; j++) {
                    jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
                    if (j == 1) {
                        // replace (j-1) data by boundary data
                        for (k = 1; k <= kdim1; k++) {
                            t(jkv + k, 7) = 0.25 * gprgm1 * (
                                q(j, k, ii+1, 5) / q(j, k, ii+1, 1) -
                                q(j, k, ii-1, 5) / q(j, k, ii-1, 1) +
                                (1.0 + bcj(k, ii+1, 1)) * qj0(k, ii+1, 5, 1) / qj0(k, ii+1, 1, 1) -
                                bcj(k, ii+1, 1) * q(j, k, ii+1, 5) / q(j, k, ii+1, 1) -
                                (1.0 + bcj(k, ii-1, 1)) * qj0(k, ii-1, 5, 1) / qj0(k, ii-1, 1, 1) +
                                bcj(k, ii-1, 1) * q(j, k, ii-1, 5) / q(j, k, ii-1, 1));
                            t(jkv + k, 8) = 0.25 * (
                                q(j, k, ii+1, 2) - q(j, k, ii-1, 2) +
                                (1.0 + bcj(k, ii+1, 1)) * qj0(k, ii+1, 2, 1) -
                                bcj(k, ii+1, 1) * q(j, k, ii+1, 2) -
                                (1.0 + bcj(k, ii-1, 1)) * qj0(k, ii-1, 2, 1) +
                                bcj(k, ii-1, 1) * q(j, k, ii-1, 2));
                            t(jkv + k, 9) = 0.25 * (
                                q(j, k, ii+1, 3) - q(j, k, ii-1, 3) +
                                (1.0 + bcj(k, ii+1, 1)) * qj0(k, ii+1, 3, 1) -
                                bcj(k, ii+1, 1) * q(j, k, ii+1, 3) -
                                (1.0 + bcj(k, ii-1, 1)) * qj0(k, ii-1, 3, 1) +
                                bcj(k, ii-1, 1) * q(j, k, ii-1, 3));
                            t(jkv + k, 10) = 0.25 * (
                                q(j, k, ii+1, 4) - q(j, k, ii-1, 4) +
                                (1.0 + bcj(k, ii+1, 1)) * qj0(k, ii+1, 4, 1) -
                                bcj(k, ii+1, 1) * q(j, k, ii+1, 4) -
                                (1.0 + bcj(k, ii-1, 1)) * qj0(k, ii-1, 4, 1) +
                                bcj(k, ii-1, 1) * q(j, k, ii-1, 4));
                        }
                    } else if (j == jdim) {
                        // replace j data by boundary data
                        for (k = 1; k <= kdim1; k++) {
                            t(jkv + k, 7) = 0.25 * gprgm1 * (
                                (1.0 + bcj(k, ii+1, 2)) * qj0(k, ii+1, 5, 3) / qj0(k, ii+1, 1, 3) -
                                bcj(k, ii+1, 2) * q(j-1, k, ii+1, 5) / q(j-1, k, ii+1, 1) -
                                (1.0 + bcj(k, ii-1, 2)) * qj0(k, ii-1, 5, 3) / qj0(k, ii-1, 1, 3) +
                                bcj(k, ii-1, 2) * q(j-1, k, ii-1, 5) / q(j-1, k, ii-1, 1) +
                                q(j-1, k, ii+1, 5) / q(j-1, k, ii+1, 1) -
                                q(j-1, k, ii-1, 5) / q(j-1, k, ii-1, 1));
                            t(jkv + k, 8) = 0.25 * (
                                (1.0 + bcj(k, ii+1, 2)) * qj0(k, ii+1, 2, 3) -
                                bcj(k, ii+1, 2) * q(j-1, k, ii+1, 2) -
                                (1.0 + bcj(k, ii-1, 2)) * qj0(k, ii-1, 2, 3) +
                                bcj(k, ii-1, 2) * q(j-1, k, ii-1, 2) +
                                q(j-1, k, ii+1, 2) - q(j-1, k, ii-1, 2));
                            t(jkv + k, 9) = 0.25 * (
                                (1.0 + bcj(k, ii+1, 2)) * qj0(k, ii+1, 3, 3) -
                                bcj(k, ii+1, 2) * q(j-1, k, ii+1, 3) -
                                (1.0 + bcj(k, ii-1, 2)) * qj0(k, ii-1, 3, 3) +
                                bcj(k, ii-1, 2) * q(j-1, k, ii-1, 3) +
                                q(j-1, k, ii+1, 3) - q(j-1, k, ii-1, 3));
                            t(jkv + k, 10) = 0.25 * (
                                (1.0 + bcj(k, ii+1, 2)) * qj0(k, ii+1, 4, 3) -
                                bcj(k, ii+1, 2) * q(j-1, k, ii+1, 4) -
                                (1.0 + bcj(k, ii-1, 2)) * qj0(k, ii-1, 4, 3) +
                                bcj(k, ii-1, 2) * q(j-1, k, ii-1, 4) +
                                q(j-1, k, ii+1, 4) - q(j-1, k, ii-1, 4));
                        }
                    } else {
                        // general case
                        for (k = 1; k <= kdim1; k++) {
                            t(jkv + k, 7) = 0.25 * gprgm1 * (
                                q(j,   k, ii+1, 5) / q(j,   k, ii+1, 1) -
                                q(j,   k, ii-1, 5) / q(j,   k, ii-1, 1) +
                                q(j-1, k, ii+1, 5) / q(j-1, k, ii+1, 1) -
                                q(j-1, k, ii-1, 5) / q(j-1, k, ii-1, 1));
                            t(jkv + k, 8) = 0.25 * (
                                q(j,   k, ii+1, 2) - q(j,   k, ii-1, 2) +
                                q(j-1, k, ii+1, 2) - q(j-1, k, ii-1, 2));
                            t(jkv + k, 9) = 0.25 * (
                                q(j,   k, ii+1, 3) - q(j,   k, ii-1, 3) +
                                q(j-1, k, ii+1, 3) - q(j-1, k, ii-1, 3));
                            t(jkv + k, 10) = 0.25 * (
                                q(j,   k, ii+1, 4) - q(j,   k, ii-1, 4) +
                                q(j-1, k, ii+1, 4) - q(j-1, k, ii-1, 4));
                        }
                    }
                } // end j loop general ii
            }
        } // end ipl loop (do 380)
    } // end if (i2d == 0 && idim > 2)


    // ******************************************************************************
    // t(24) : turbulent viscosity at interfaces (=0 for laminar flow)
    iviscc = ivisc[1];  // ivisc(2) -> ivisc[1] (0-based)
    if (iviscc > 1) {
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            for (j = 1; j <= jdim; j++) {
                jkv = (j - 1) * npl * kdim1 + (ipl - 1) * kdim1;
                // j=0 interface
                if (j == 1) {
                    for (k = 1; k <= kdim1; k++) {
                        t(k + jkv, 24) = bcj(k, ii, 1) * vj0(k, ii, 1, 1) +
                            (1.0 - bcj(k, ii, 1)) * 0.5 * (vj0(k, ii, 1, 1) + vist3d(1, k, ii));
                    }
                // j=jdim interface
                } else if (j == jdim) {
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

    // ******************************************************************************
    // compute t(23) and t(14) final
    for (izz = 1; izz <= l0; izz++) {
        // ratio of turbulent to laminar viscosity
        t25 = t(izz, 24) / t(izz, 14);
        t24 = 1.e0 + t25 * (double)coef_eddy;
        // t(23) : ratio of laminar Prandtl number to total Prandtl number
        t(izz, 23) = (1.e0 + prtr * t25) / t24;
        // t(14): (mach/re)*viscosity/J
        t(izz, 14) = xmre * t24 * t(izz, 14) / t(izz, 15);
    }

    if (iadv < 0) return;

    // ******************************************************************************
    // t(11) : u at cell interfaces
    // t(12) : v at cell interfaces
    // t(13) : w at cell interfaces
    // interior interfaces
    for (izz = 1; izz <= nn; izz++) {
        t(izz + ks - 1, 11) = 0.5e0 * (t(izz + ks - 1, 18) + t(izz, 18));
        t(izz + ks - 1, 12) = 0.5e0 * (t(izz + ks - 1, 19) + t(izz, 19));
        t(izz + ks - 1, 13) = 0.5e0 * (t(izz + ks - 1, 20) + t(izz, 20));
    }
    // interfaces at j=0 and j=jdim
    for (izz = 1; izz <= kv; izz++) {
        ab = 1.0 + wj0(izz, 15);
        bb = 1.0 - wj0(izz, 15);
        t(izz + n, 11) = 0.5 * (ab * wj0(izz, 18) + bb * t(izz + n - kv, 18));
        t(izz + n, 12) = 0.5 * (ab * wj0(izz, 19) + bb * t(izz + n - kv, 19));
        t(izz + n, 13) = 0.5 * (ab * wj0(izz, 20) + bb * t(izz + n - kv, 20));
        ab = 1.0 + wj0(izz, 5);
        bb = 1.0 - wj0(izz, 5);
        t(izz, 11) = 0.5 * (ab * wj0(izz, 8)  + bb * t(izz, 18));
        t(izz, 12) = 0.5 * (ab * wj0(izz, 9)  + bb * t(izz, 19));
        t(izz, 13) = 0.5 * (ab * wj0(izz, 10) + bb * t(izz, 20));
    }

    // ******************************************************************************
    // calculate fluxes
    // viscous terms at interfaces
    for (izz = 1; izz <= l0; izz++) {
        // form ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty and tz
        // excluding contributions in j/eta direction
        ux = t(izz, 4) * t(izz, 28) + t(izz, 8) * t(izz, 31);
        uy = t(izz, 4) * t(izz, 29) + t(izz, 8) * t(izz, 32);
        uz = t(izz, 4) * t(izz, 30) + t(izz, 8) * t(izz, 2);

        vx = t(izz, 5) * t(izz, 28) + t(izz, 9) * t(izz, 31);
        vy = t(izz, 5) * t(izz, 29) + t(izz, 9) * t(izz, 32);
        vz = t(izz, 5) * t(izz, 30) + t(izz, 9) * t(izz, 2);

        wx = t(izz, 6) * t(izz, 28) + t(izz, 10) * t(izz, 31);
        wy = t(izz, 6) * t(izz, 29) + t(izz, 10) * t(izz, 32);
        wz = t(izz, 6) * t(izz, 30) + t(izz, 10) * t(izz, 2);

        tx = t(izz, 3) * t(izz, 28) + t(izz, 7) * t(izz, 31);
        ty = t(izz, 3) * t(izz, 29) + t(izz, 7) * t(izz, 32);
        tz = t(izz, 3) * t(izz, 30) + t(izz, 7) * t(izz, 2);

        // Note lambda = -2/3 mu, lambda + 2 mu = 4/3 mu
        // u-momentum flux
        t(izz, 16) = t(izz, 14) * (
            t(izz, 25) * (4.0/3.0 * ux - 2.0/3.0 * (vy + wz))
          + t(izz, 26) * (uy + vx)
          + t(izz, 27) * (uz + wx));

        // v-momentum flux
        t(izz, 17) = t(izz, 14) * (
            t(izz, 25) * (vx + uy)
          + t(izz, 26) * (4.0/3.0 * vy - 2.0/3.0 * (ux + wz))
          + t(izz, 27) * (vz + wy));

        // w-momentum flux
        t(izz, 18) = t(izz, 14) * (
            t(izz, 25) * (wx + uz)
          + t(izz, 26) * (wy + vz)
          + t(izz, 27) * (4.0/3.0 * wz - 2.0/3.0 * (vy + ux)));

        // energy flux
        t(izz, 19) = t(izz, 14) * (
            t(izz, 25) * ((4.0/3.0 * ux - 2.0/3.0 * (vy + wz)) * t(izz, 11) +
                          (uy + vx) * t(izz, 12) + (uz + wx) * t(izz, 13)) +
            t(izz, 26) * ((vx + uy) * t(izz, 11) + (vz + wy) * t(izz, 13) +
                          (4.0/3.0 * vy - 2.0/3.0 * (ux + wz)) * t(izz, 12)) +
            t(izz, 27) * ((wx + uz) * t(izz, 11) + (wy + vz) * t(izz, 12) +
                          (4.0/3.0 * wz - 2.0/3.0 * (ux + vy)) * t(izz, 13)) +
            t(izz, 23) * (t(izz, 25) * tx + t(izz, 26) * ty + t(izz, 27) * tz));
    }

    // (-)viscous flux = Gv(j-1/2) - Gv(j+1/2)
    for (izz = 1; izz <= n; izz++) {
        t(izz, 2) = -t(izz + ks - 1, 16) + t(izz, 16);
        t(izz, 3) = -t(izz + ks - 1, 17) + t(izz, 17);
        t(izz, 4) = -t(izz + ks - 1, 18) + t(izz, 18);
        t(izz, 5) = -t(izz + ks - 1, 19) + t(izz, 19);
    }

    // ******************************************************************************
    // calculate residuals
    for (ipl = 1; ipl <= npl; ipl++) {
        ii  = i + ipl - 1;
        for (j = 1; j <= jdim1; j++) {
            jkv = (j - 1) * kdim1 * npl + (ipl - 1) * kdim1;
            // for 2-D, t(3) should be identically zero
            if (i2d == 1) {
                for (k = 1; k <= kdim1; k++) {
                    t(k + jkv, 3) = 0.0;
                }
            }
            for (int l = 2; l <= 5; l++) {
                for (k = 1; k <= kdim1; k++) {
                    res(j, k, ii, l) = res(j, k, ii, l) + t(k + jkv, l);
                }
            }
        }
    }
    // return
}

} // namespace gfluxv1_ns
