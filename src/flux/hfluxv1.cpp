// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "hfluxv1.h"
#include "hfluxv.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace hfluxv1_ns {

void hfluxv(int& i, int& npl, int& jdim, int& kdim, int& idim, int& idf,
            FortranArray4DRef<double> ak, FortranArray4DRef<double> bk,
            FortranArray4DRef<double> ck, FortranArray4DRef<double> res,
            FortranArray4DRef<double> q, FortranArray4DRef<double> qk0,
            FortranArray3DRef<double> sk, FortranArray3DRef<double> vol,
            FortranArray2DRef<double> t, int& nvtq, FortranArray2DRef<double> wk0,
            FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vmuk,
            FortranArray4DRef<double> vk0, FortranArray3DRef<double> bck,
            FortranArray4DRef<double> zksav, FortranArray4DRef<double> tk0,
            FortranArray3DRef<double> cmuv, FortranArray3DRef<double> volk0,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& iadv, int& nummem,
            FortranArray4DRef<double> ux)
{
    hfluxv_ns::hfluxv(i, npl, jdim, kdim, idim, idf, ak, bk, ck, res, q, qk0,
                      sk, vol, t, nvtq, wk0, vist3d, vmuk, vk0, bck, zksav,
                      tk0, cmuv, volk0, nou, bou, nbuf, ibufdim, iadv, nummem, ux);
}

void hfluxv1(int& i, int& npl, int& jdim, int& kdim, int& idim,
             FortranArray4DRef<double> res, FortranArray4DRef<double> q,
             FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
             FortranArray4DRef<double> qi0, FortranArray4DRef<double> sj,
             FortranArray3DRef<double> sk, FortranArray4DRef<double> si,
             FortranArray3DRef<double> vol, FortranArray2DRef<double> t,
             int& nvtq, FortranArray2DRef<double> wk0,
             FortranArray3DRef<double> vist3d, FortranArray4DRef<double> vk0,
             FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
             FortranArray3DRef<double> bci, FortranArray3DRef<double> volk0,
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
    int&   i2d     = cmn_twod.i2d;
    int&   i_tauijs = cmn_constit.i_tauijs;

    // Local variables
    float  coef_eddy;
    int    kdim1, jdim1, idim1;
    int    n, l0, jv, js, nn;
    double xmre, gpr, gm1pr, prtr, gprgm1;
    double c2b, c2bp;
    int    ipl, ii, k, j, izz, jkv, jk, l;
    int    n1, l1;
    int    iviscc;
    double t5, t6, ab, bb, wk05, wk06;
    double t1;
    double ux_v, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz;
    double t24, t25;
    int    mm, mb, imin, imax, jmin, jmax, jc, m;



    // Lines 75-101: Initialize coef_eddy, dimensions, derived constants
    coef_eddy = 1.0f;
    if (ivisc[2] >= 70) coef_eddy = 0.0f;
    if (i_tauijs == 1) coef_eddy = 0.0f;
    kdim1 = kdim - 1;
    jdim1 = jdim - 1;
    idim1 = idim - 1;

    // n  : number of cell centers for npl planes
    // l0 : number of cell interfaces for npl planes
    // jv : number of cell centers (and interfaces) on a k=constant plane
    n  = npl * jdim1 * kdim1;
    l0 = npl * jdim1 * kdim;
    jv = npl * jdim1;
    js = jv + 1;
    nn = n - jv;

    xmre   = (double)xmach / (double)reue;
    gpr    = (double)gamma / (double)pr;
    gm1pr  = (double)gm1 * (double)pr;
    prtr   = (double)pr / (double)prt;
    gprgm1 = gpr / (double)gm1;

    if (isklton > 0 && i == 1 && iadv >= 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "   computing cross-derivative viscous fluxes, K-dir.");
    }

    // Lines 114-128: Store selected cell centered data in t array
    // t(22)=density, t(18)=u, t(19)=v, t(20)=w, t(16)=c*c/(gm1*pr)
    l1 = jdim * kdim - 1;
    for (ipl = 1; ipl <= npl; ipl++) {
        ii = i + ipl - 1;
        for (k = 1; k <= kdim1; k++) {
            jkv = (k - 1) * npl * jdim1 + (ipl - 1) * jdim1 + 1;
            for (izz = 1; izz <= jdim1; izz++) {
                t(izz + jkv - 1, 22) = q(izz, k, ii, 1);
                t(izz + jkv - 1, 18) = q(izz, k, ii, 2);
                t(izz + jkv - 1, 19) = q(izz, k, ii, 3);
                t(izz + jkv - 1, 20) = q(izz, k, ii, 4);
                t(izz + jkv - 1, 16) = q(izz, k, ii, 5);
            }
        }
    }

    // Lines 134-160: Store dependent variables on k=0 and k=kdim boundaries in wk0
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
        jc = 0;
        for (ii = imin; ii <= imax; ii++) {
            jc = (ii - i) * jdim1 + jmin - 1;
            for (j = jmin; j <= jmax; j++) {
                jc = jc + 1;
                wk0(jc, m + 12) = qk0(j, ii, 1, mm);
                wk0(jc, m + 8)  = qk0(j, ii, 2, mm);
                wk0(jc, m + 9)  = qk0(j, ii, 3, mm);
                wk0(jc, m + 10) = qk0(j, ii, 4, mm);
                wk0(jc, m + 6)  = qk0(j, ii, 5, mm);
                // m+5 flag indicates whether ghost cell or interface values stored in wk0
                wk0(jc, m + 5)  = bck(j, ii, mb);
            }
        }
    }

    // Lines 170-186: Compute t(17), t(1), t(16) and wk0 derived quantities
    for (izz = 1; izz <= n; izz++) {
        t(izz, 17) = 0.5e0 * (t(izz, 18) * t(izz, 18)
                             + t(izz, 19) * t(izz, 19)
                             + t(izz, 20) * t(izz, 20));
        t(izz, 1)  = 1.e0 / t(izz, 22);
        t(izz, 16) = gpr * t(izz, 16) * t(izz, 1) / (double)gm1;
    }
    for (m = 0; m <= 10; m += 10) {
        for (izz = 1; izz <= jv; izz++) {
            wk0(izz, m + 7) = 0.5e0 * (wk0(izz, m + 8) * wk0(izz, m + 8)
                                       + wk0(izz, m + 9) * wk0(izz, m + 9)
                                       + wk0(izz, m + 10) * wk0(izz, m + 10));
            wk0(izz, m + 1) = 1.e0 / wk0(izz, m + 12);
            wk0(izz, m + 6) = gpr * wk0(izz, m + 6) * wk0(izz, m + 1) / (double)gm1;
        }
    }

    // Lines 192-200: t(7) laminar viscosity at cell centers (Sutherland relation)
    c2b  = (double)cbar / (double)tinf;
    c2bp = c2b + 1.e0;
    for (izz = 1; izz <= n; izz++) {
        t5       = gm1pr * t(izz, 16);
        t6       = std::sqrt(t5);
        t(izz, 7) = c2bp * t5 * t6 / (c2b + t5);
    }

    // Lines 207-225: t(14) laminar viscosity at cell interfaces
    // interior interfaces
    for (izz = 1; izz <= nn; izz++) {
        t(izz + js - 1, 14) = (t(izz, 7) + t(izz + js - 1, 7)) * 0.5e0;
    }
    // k=0 and k=kdim interfaces
    for (izz = 1; izz <= jv; izz++) {
        ab   = 1.0 + wk0(izz, 5);
        bb   = 1.0 - wk0(izz, 5);
        wk05 = gm1pr * 0.5 * (ab * wk0(izz, 6) + bb * t(izz, 16));
        wk06 = std::sqrt(wk05);
        t(izz, 14)     = c2bp * wk05 * wk06 / (c2b + wk05);
        ab   = 1.0 + wk0(izz, 15);
        bb   = 1.0 - wk0(izz, 15);
        wk05 = gm1pr * 0.5 * (ab * wk0(izz, 16) + bb * t(izz + n - jv, 16));
        wk06 = std::sqrt(wk05);
        t(izz + n, 14) = c2bp * wk05 * wk06 / (c2b + wk05);
    }



    // Lines 231-269: t(15) average jacobian (inverse volume) at cell interface
    n1 = jdim * kdim1 + 1;
    l1 = jdim * kdim - 1;
    for (ipl = 1; ipl <= npl; ipl++) {
        ii  = i + ipl - 1;
        for (k = 1; k <= kdim; k++) {
            jkv = (k - 1) * npl * jdim1 + (ipl - 1) * jdim1;
            jk  = (k - 1) * jdim;
            if (k == 1) {
                // inverse volume at k=0 interface
                for (j = 1; j <= jdim1; j++) {
                    t(j + jkv, 15) = 2.0 / (volk0(j, ii, 1) + vol(j, 1, ii));
                }
            } else if (k == kdim) {
                // inverse volume at k=kdim interface
                for (j = 1; j <= jdim1; j++) {
                    t(j + jkv, 15) = 2.0 / (volk0(j, ii, 3) + vol(j, kdim1, ii));
                }
            } else {
                // inverse volume at interior interfaces
                for (j = 1; j <= jdim1; j++) {
                    t(j + jkv, 15) = 2.0 / (vol(j, k, ii) + vol(j, k - 1, ii));
                }
            }

            // Lines 263-268: t(25-27) components of grad(zeta)
            for (j = 1; j <= jdim1; j++) {
                t1          = sk(j + jk, ii, 4) * t(j + jkv, 15);
                t(j + jkv, 25) = sk(j + jk, ii, 1) * t1;
                t(j + jkv, 26) = sk(j + jk, ii, 2) * t1;
                t(j + jkv, 27) = sk(j + jk, ii, 3) * t1;
            }
        }
    }

    // Lines 280-333: t(28-30) Components of grad(xie)
    for (ipl = 1; ipl <= npl; ipl++) {
        ii  = i + ipl - 1;
        for (k = 1; k <= kdim; k++) {
            jkv = (k - 1) * npl * jdim1 + (ipl - 1) * jdim1;
            if (k == 1) {
                for (j = 1; j <= jdim1; j++) {
                    t(jkv + j, 28) = 0.25 * t(j + jkv, 15) * (si(j, k, ii, 1) * si(j, k, ii, 4)
                                   + si(j, k, ii, 1) * si(j, k, ii, 4)
                                   + si(j, k, ii + 1, 1) * si(j, k, ii + 1, 4)
                                   + si(j, k, ii + 1, 1) * si(j, k, ii + 1, 4));
                    t(jkv + j, 29) = 0.25 * t(j + jkv, 15) * (si(j, k, ii, 2) * si(j, k, ii, 4)
                                   + si(j, k, ii, 2) * si(j, k, ii, 4)
                                   + si(j, k, ii + 1, 2) * si(j, k, ii + 1, 4)
                                   + si(j, k, ii + 1, 2) * si(j, k, ii + 1, 4));
                    t(jkv + j, 30) = 0.25 * t(j + jkv, 15) * (si(j, k, ii, 3) * si(j, k, ii, 4)
                                   + si(j, k, ii, 3) * si(j, k, ii, 4)
                                   + si(j, k, ii + 1, 3) * si(j, k, ii + 1, 4)
                                   + si(j, k, ii + 1, 3) * si(j, k, ii + 1, 4));
                }
            } else if (k == kdim) {
                for (j = 1; j <= jdim1; j++) {
                    t(jkv + j, 28) = 0.25 * t(j + jkv, 15) * (si(j, k - 1, ii, 1) * si(j, k - 1, ii, 4)
                                   + si(j, k - 1, ii, 1) * si(j, k - 1, ii, 4)
                                   + si(j, k - 1, ii + 1, 1) * si(j, k - 1, ii + 1, 4)
                                   + si(j, k - 1, ii + 1, 1) * si(j, k - 1, ii + 1, 4));
                    t(jkv + j, 29) = 0.25 * t(j + jkv, 15) * (si(j, k - 1, ii, 2) * si(j, k - 1, ii, 4)
                                   + si(j, k - 1, ii, 2) * si(j, k - 1, ii, 4)
                                   + si(j, k - 1, ii + 1, 2) * si(j, k - 1, ii + 1, 4)
                                   + si(j, k - 1, ii + 1, 2) * si(j, k - 1, ii + 1, 4));
                    t(jkv + j, 30) = 0.25 * t(j + jkv, 15) * (si(j, k - 1, ii, 3) * si(j, k - 1, ii, 4)
                                   + si(j, k - 1, ii, 3) * si(j, k - 1, ii, 4)
                                   + si(j, k - 1, ii + 1, 3) * si(j, k - 1, ii + 1, 4)
                                   + si(j, k - 1, ii + 1, 3) * si(j, k - 1, ii + 1, 4));
                }
            } else {
                // general case
                for (j = 1; j <= jdim1; j++) {
                    t(jkv + j, 28) = 0.25 * t(j + jkv, 15) * (si(j, k, ii, 1) * si(j, k, ii, 4)
                                   + si(j, k - 1, ii, 1) * si(j, k - 1, ii, 4)
                                   + si(j, k, ii + 1, 1) * si(j, k, ii + 1, 4)
                                   + si(j, k - 1, ii + 1, 1) * si(j, k - 1, ii + 1, 4));
                    t(jkv + j, 29) = 0.25 * t(j + jkv, 15) * (si(j, k, ii, 2) * si(j, k, ii, 4)
                                   + si(j, k - 1, ii, 2) * si(j, k - 1, ii, 4)
                                   + si(j, k, ii + 1, 2) * si(j, k, ii + 1, 4)
                                   + si(j, k - 1, ii + 1, 2) * si(j, k - 1, ii + 1, 4));
                    t(jkv + j, 30) = 0.25 * t(j + jkv, 15) * (si(j, k, ii, 3) * si(j, k, ii, 4)
                                   + si(j, k - 1, ii, 3) * si(j, k - 1, ii, 4)
                                   + si(j, k, ii + 1, 3) * si(j, k, ii + 1, 4)
                                   + si(j, k - 1, ii + 1, 3) * si(j, k - 1, ii + 1, 4));
                }
            }
        }
    }



    // Lines 344-397: t(31,32,2) Components of grad(eta)
    for (ipl = 1; ipl <= npl; ipl++) {
        ii  = i + ipl - 1;
        for (k = 1; k <= kdim; k++) {
            jkv = (k - 1) * npl * jdim1 + (ipl - 1) * jdim1;
            if (k == 1) {
                for (j = 1; j <= jdim1; j++) {
                    t(jkv + j, 31) = 0.25 * t(j + jkv, 15) * (sj(j, k, ii, 1) * sj(j, k, ii, 4)
                                   + sj(j, k, ii, 1) * sj(j, k, ii, 4)
                                   + sj(j + 1, k, ii, 1) * sj(j + 1, k, ii, 4)
                                   + sj(j + 1, k, ii, 1) * sj(j + 1, k, ii, 4));
                    t(jkv + j, 32) = 0.25 * t(j + jkv, 15) * (sj(j, k, ii, 2) * sj(j, k, ii, 4)
                                   + sj(j, k, ii, 2) * sj(j, k, ii, 4)
                                   + sj(j + 1, k, ii, 2) * sj(j + 1, k, ii, 4)
                                   + sj(j + 1, k, ii, 2) * sj(j + 1, k, ii, 4));
                    t(jkv + j, 2)  = 0.25 * t(j + jkv, 15) * (sj(j, k, ii, 3) * sj(j, k, ii, 4)
                                   + sj(j, k, ii, 3) * sj(j, k, ii, 4)
                                   + sj(j + 1, k, ii, 3) * sj(j + 1, k, ii, 4)
                                   + sj(j + 1, k, ii, 3) * sj(j + 1, k, ii, 4));
                }
            } else if (k == kdim) {
                for (j = 1; j <= jdim1; j++) {
                    t(jkv + j, 31) = 0.25 * t(j + jkv, 15) * (sj(j, k - 1, ii, 1) * sj(j, k - 1, ii, 4)
                                   + sj(j, k - 1, ii, 1) * sj(j, k - 1, ii, 4)
                                   + sj(j + 1, k - 1, ii, 1) * sj(j + 1, k - 1, ii, 4)
                                   + sj(j + 1, k - 1, ii, 1) * sj(j + 1, k - 1, ii, 4));
                    t(jkv + j, 32) = 0.25 * t(j + jkv, 15) * (sj(j, k - 1, ii, 2) * sj(j, k - 1, ii, 4)
                                   + sj(j, k - 1, ii, 2) * sj(j, k - 1, ii, 4)
                                   + sj(j + 1, k - 1, ii, 2) * sj(j + 1, k - 1, ii, 4)
                                   + sj(j + 1, k - 1, ii, 2) * sj(j + 1, k - 1, ii, 4));
                    t(jkv + j, 2)  = 0.25 * t(j + jkv, 15) * (sj(j, k - 1, ii, 3) * sj(j, k - 1, ii, 4)
                                   + sj(j, k - 1, ii, 3) * sj(j, k - 1, ii, 4)
                                   + sj(j + 1, k - 1, ii, 3) * sj(j + 1, k - 1, ii, 4)
                                   + sj(j + 1, k - 1, ii, 3) * sj(j + 1, k - 1, ii, 4));
                }
            } else {
                // general case
                for (j = 1; j <= jdim1; j++) {
                    t(jkv + j, 31) = 0.25 * t(j + jkv, 15) * (sj(j, k, ii, 1) * sj(j, k, ii, 4)
                                   + sj(j, k - 1, ii, 1) * sj(j, k - 1, ii, 4)
                                   + sj(j + 1, k, ii, 1) * sj(j + 1, k, ii, 4)
                                   + sj(j + 1, k - 1, ii, 1) * sj(j + 1, k - 1, ii, 4));
                    t(jkv + j, 32) = 0.25 * t(j + jkv, 15) * (sj(j, k, ii, 2) * sj(j, k, ii, 4)
                                   + sj(j, k - 1, ii, 2) * sj(j, k - 1, ii, 4)
                                   + sj(j + 1, k, ii, 2) * sj(j + 1, k, ii, 4)
                                   + sj(j + 1, k - 1, ii, 2) * sj(j + 1, k - 1, ii, 4));
                    t(jkv + j, 2)  = 0.25 * t(j + jkv, 15) * (sj(j, k, ii, 3) * sj(j, k, ii, 4)
                                   + sj(j, k - 1, ii, 3) * sj(j, k - 1, ii, 4)
                                   + sj(j + 1, k, ii, 3) * sj(j + 1, k, ii, 4)
                                   + sj(j + 1, k - 1, ii, 3) * sj(j + 1, k - 1, ii, 4));
                }
            }
        }
    }



    // Lines 414-768: gradients at cell interfaces (eta direction)
    // t(7): d(c*c/(gm1*pr))/d(zeta), t(8): d(u)/d(eta), t(9): d(v)/d(eta), t(10): d(w)/d(eta)
    for (ipl = 1; ipl <= npl; ipl++) {
        ii = i + ipl - 1;
        for (k = 1; k <= kdim; k++) {
            jkv = (k - 1) * npl * jdim1 + (ipl - 1) * jdim1;

            if (k == 1) {
                // replace (k-1) data by boundary data
                for (j = 1; j <= jdim1; j++) {
                    if (j == 1) {
                        t(jkv + j, 7) = 0.25 * gprgm1 * (
                            q(j + 1, k, ii, 5) / q(j + 1, k, ii, 1) -
                            (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 5, 1) / qj0(k, ii, 1, 1) +
                            bcj(k, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) +
                            2.0 * (1.0 + bck(j + 1, ii, 1)) * qk0(j + 1, ii, 5, 1) / qk0(j + 1, ii, 1, 1) -
                            2.0 * bck(j + 1, ii, 1) * q(j + 1, k, ii, 5) / q(j + 1, k, ii, 1) -
                            2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 5, 1) / qk0(j, ii, 1, 1) +
                            2.0 * bck(j, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1));
                        t(jkv + j, 8) = 0.25 * (
                            q(j + 1, k, ii, 2) -
                            (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 2, 1) +
                            bcj(k, ii, 1) * q(j, k, ii, 2) +
                            2.0 * (1.0 + bck(j + 1, ii, 1)) * qk0(j + 1, ii, 2, 1) -
                            2.0 * bck(j + 1, ii, 1) * q(j + 1, k, ii, 2) -
                            2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 2, 1) +
                            2.0 * bck(j, ii, 1) * q(j, k, ii, 2));
                        t(jkv + j, 9) = 0.25 * (
                            q(j + 1, k, ii, 3) -
                            (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 3, 1) +
                            bcj(k, ii, 1) * q(j, k, ii, 3) +
                            2.0 * (1.0 + bck(j + 1, ii, 1)) * qk0(j + 1, ii, 3, 1) -
                            2.0 * bck(j + 1, ii, 1) * q(j + 1, k, ii, 3) -
                            2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 3, 1) +
                            2.0 * bck(j, ii, 1) * q(j, k, ii, 3));
                        t(jkv + j, 10) = 0.25 * (
                            q(j + 1, k, ii, 4) -
                            (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 4, 1) +
                            bcj(k, ii, 1) * q(j, k, ii, 4) +
                            2.0 * (1.0 + bck(j + 1, ii, 1)) * qk0(j + 1, ii, 4, 1) -
                            2.0 * bck(j + 1, ii, 1) * q(j + 1, k, ii, 4) -
                            2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 4, 1) +
                            2.0 * bck(j, ii, 1) * q(j, k, ii, 4));
                    } else if (j == jdim1) {
                        t(jkv + j, 7) = 0.25 * gprgm1 * (
                            (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 5, 3) / qj0(k, ii, 1, 3) -
                            bcj(k, ii, 2) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                            q(j - 1, k, ii, 5) / q(j - 1, k, ii, 1) +
                            2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 5, 1) / qk0(j, ii, 1, 1) -
                            2.0 * bck(j, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                            2.0 * (1.0 + bck(j - 1, ii, 1)) * qk0(j - 1, ii, 5, 1) / qk0(j - 1, ii, 1, 1) +
                            2.0 * bck(j - 1, ii, 1) * q(j - 1, k, ii, 5) / q(j - 1, k, ii, 1));
                        t(jkv + j, 8) = 0.25 * (
                            (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 2, 3) -
                            bcj(k, ii, 2) * q(j, k, ii, 2) -
                            q(j - 1, k, ii, 2) +
                            2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 2, 1) -
                            2.0 * bck(j, ii, 1) * q(j, k, ii, 2) -
                            2.0 * (1.0 + bck(j - 1, ii, 1)) * qk0(j - 1, ii, 2, 1) +
                            2.0 * bck(j - 1, ii, 1) * q(j - 1, k, ii, 2));
                        t(jkv + j, 9) = 0.25 * (
                            (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 3, 3) -
                            bcj(k, ii, 2) * q(j, k, ii, 3) -
                            q(j - 1, k, ii, 3) +
                            2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 3, 1) -
                            2.0 * bck(j, ii, 1) * q(j, k, ii, 3) -
                            2.0 * (1.0 + bck(j - 1, ii, 1)) * qk0(j - 1, ii, 3, 1) +
                            2.0 * bck(j - 1, ii, 1) * q(j - 1, k, ii, 3));
                        t(jkv + j, 10) = 0.25 * (
                            (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 4, 3) -
                            bcj(k, ii, 2) * q(j, k, ii, 4) -
                            q(j - 1, k, ii, 4) +
                            2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 4, 1) -
                            2.0 * bck(j, ii, 1) * q(j, k, ii, 4) -
                            2.0 * (1.0 + bck(j - 1, ii, 1)) * qk0(j - 1, ii, 4, 1) +
                            2.0 * bck(j - 1, ii, 1) * q(j - 1, k, ii, 4));
                    } else {
                        // general case for k = 1
                        t(jkv + j, 7) = 0.25 * gprgm1 * (
                            q(j + 1, k, ii, 5) / q(j + 1, k, ii, 1) -
                            q(j - 1, k, ii, 5) / q(j - 1, k, ii, 1) +
                            (1.0 + bck(j + 1, ii, 1)) * qk0(j + 1, ii, 5, 1) / qk0(j + 1, ii, 1, 1) -
                            bck(j + 1, ii, 1) * q(j + 1, k, ii, 5) / q(j + 1, k, ii, 1) -
                            (1.0 + bck(j - 1, ii, 1)) * qk0(j - 1, ii, 5, 1) / qk0(j - 1, ii, 1, 1) +
                            bck(j - 1, ii, 1) * q(j - 1, k, ii, 5) / q(j - 1, k, ii, 1));
                        t(jkv + j, 8) = 0.25 * (
                            q(j + 1, k, ii, 2) - q(j - 1, k, ii, 2) +
                            (1.0 + bck(j + 1, ii, 1)) * qk0(j + 1, ii, 2, 1) -
                            bck(j + 1, ii, 1) * q(j + 1, k, ii, 2) -
                            (1.0 + bck(j - 1, ii, 1)) * qk0(j - 1, ii, 2, 1) +
                            bck(j - 1, ii, 1) * q(j - 1, k, ii, 2));
                        t(jkv + j, 9) = 0.25 * (
                            q(j + 1, k, ii, 3) - q(j - 1, k, ii, 3) +
                            (1.0 + bck(j + 1, ii, 1)) * qk0(j + 1, ii, 3, 1) -
                            bck(j + 1, ii, 1) * q(j + 1, k, ii, 3) -
                            (1.0 + bck(j - 1, ii, 1)) * qk0(j - 1, ii, 3, 1) +
                            bck(j - 1, ii, 1) * q(j - 1, k, ii, 3));
                        t(jkv + j, 10) = 0.25 * (
                            q(j + 1, k, ii, 4) - q(j - 1, k, ii, 4) +
                            (1.0 + bck(j + 1, ii, 1)) * qk0(j + 1, ii, 4, 1) -
                            bck(j + 1, ii, 1) * q(j + 1, k, ii, 4) -
                            (1.0 + bck(j - 1, ii, 1)) * qk0(j - 1, ii, 4, 1) +
                            bck(j - 1, ii, 1) * q(j - 1, k, ii, 4));
                    }
                } // end j loop for k==1



            } else if (k == kdim) {
                // replace k data by boundary data
                for (j = 1; j <= jdim1; j++) {
                    if (j == 1) {
                        t(jkv + j, 7) = 0.25 * gprgm1 * (
                            2.0 * (1.0 + bck(j + 1, ii, 2)) * qk0(j + 1, ii, 5, 3) / qk0(j + 1, ii, 1, 3) -
                            2.0 * bck(j + 1, ii, 2) * q(j + 1, k - 1, ii, 5) / q(j + 1, k - 1, ii, 1) -
                            2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 5, 3) / qk0(j, ii, 1, 3) +
                            2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1) +
                            q(j + 1, k - 1, ii, 5) / q(j + 1, k - 1, ii, 1) -
                            (1.0 + bcj(k - 1, ii, 1)) * qj0(k - 1, ii, 5, 1) / qj0(k - 1, ii, 1, 1) +
                            bcj(k - 1, ii, 1) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1));
                        t(jkv + j, 8) = 0.25 * (
                            2.0 * (1.0 + bck(j + 1, ii, 2)) * qk0(j + 1, ii, 2, 3) -
                            2.0 * bck(j + 1, ii, 2) * q(j + 1, k - 1, ii, 2) -
                            2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 2, 3) +
                            2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 2) +
                            q(j + 1, k - 1, ii, 2) -
                            (1.0 + bcj(k - 1, ii, 1)) * qj0(k - 1, ii, 2, 1) +
                            bcj(k - 1, ii, 1) * q(j, k - 1, ii, 2));
                        t(jkv + j, 9) = 0.25 * (
                            2.0 * (1.0 + bck(j + 1, ii, 2)) * qk0(j + 1, ii, 3, 3) -
                            2.0 * bck(j + 1, ii, 2) * q(j + 1, k - 1, ii, 3) -
                            2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 3, 3) +
                            2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 3) +
                            q(j + 1, k - 1, ii, 3) -
                            (1.0 + bcj(k - 1, ii, 1)) * qj0(k - 1, ii, 3, 1) +
                            bcj(k - 1, ii, 1) * q(j, k - 1, ii, 3));
                        t(jkv + j, 10) = 0.25 * (
                            2.0 * (1.0 + bck(j + 1, ii, 2)) * qk0(j + 1, ii, 4, 3) -
                            2.0 * bck(j + 1, ii, 2) * q(j + 1, k - 1, ii, 4) -
                            2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 4, 3) +
                            2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 4) +
                            q(j + 1, k - 1, ii, 4) -
                            (1.0 + bcj(k - 1, ii, 1)) * qj0(k - 1, ii, 4, 1) +
                            bcj(k - 1, ii, 1) * q(j, k - 1, ii, 4));
                    } else if (j == jdim1) {
                        t(jkv + j, 7) = 0.25 * gprgm1 * (
                            2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 5, 3) / qk0(j, ii, 1, 3) -
                            2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1) -
                            2.0 * (1.0 + bck(j - 1, ii, 2)) * qk0(j - 1, ii, 5, 3) / qk0(j - 1, ii, 1, 3) +
                            2.0 * bck(j - 1, ii, 2) * q(j - 1, k - 1, ii, 5) / q(j - 1, k - 1, ii, 1) +
                            (1.0 + bcj(k - 1, ii, 2)) * qj0(k - 1, ii, 5, 3) / qj0(k - 1, ii, 1, 3) -
                            bcj(k - 1, ii, 2) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1) -
                            q(j - 1, k - 1, ii, 5) / q(j - 1, k - 1, ii, 1));
                        t(jkv + j, 8) = 0.25 * (
                            2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 2, 3) -
                            2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 2) -
                            2.0 * (1.0 + bck(j - 1, ii, 2)) * qk0(j - 1, ii, 2, 3) +
                            2.0 * bck(j - 1, ii, 2) * q(j - 1, k - 1, ii, 2) +
                            (1.0 + bcj(k - 1, ii, 2)) * qj0(k - 1, ii, 2, 3) -
                            bcj(k - 1, ii, 2) * q(j, k - 1, ii, 2) -
                            q(j - 1, k - 1, ii, 2));
                        t(jkv + j, 9) = 0.25 * (
                            2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 3, 3) -
                            2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 3) -
                            2.0 * (1.0 + bck(j - 1, ii, 2)) * qk0(j - 1, ii, 3, 3) +
                            2.0 * bck(j - 1, ii, 2) * q(j - 1, k - 1, ii, 3) +
                            (1.0 + bcj(k - 1, ii, 2)) * qj0(k - 1, ii, 3, 3) -
                            bcj(k - 1, ii, 2) * q(j, k - 1, ii, 3) -
                            q(j - 1, k - 1, ii, 3));
                        t(jkv + j, 10) = 0.25 * (
                            2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 4, 3) -
                            2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 4) -
                            2.0 * (1.0 + bck(j - 1, ii, 2)) * qk0(j - 1, ii, 4, 3) +
                            2.0 * bck(j - 1, ii, 2) * q(j - 1, k - 1, ii, 4) +
                            (1.0 + bcj(k - 1, ii, 2)) * qj0(k - 1, ii, 4, 3) -
                            bcj(k - 1, ii, 2) * q(j, k - 1, ii, 4) -
                            q(j - 1, k - 1, ii, 4));
                    } else {
                        // general case for k = kdim
                        t(jkv + j, 7) = 0.25 * gprgm1 * (
                            (1.0 + bck(j + 1, ii, 2)) * qk0(j + 1, ii, 5, 3) / qk0(j + 1, ii, 1, 3) -
                            bck(j + 1, ii, 2) * q(j + 1, k - 1, ii, 5) / q(j + 1, k - 1, ii, 1) -
                            (1.0 + bck(j - 1, ii, 2)) * qk0(j - 1, ii, 5, 3) / qk0(j - 1, ii, 1, 3) +
                            bck(j - 1, ii, 2) * q(j - 1, k - 1, ii, 5) / q(j - 1, k - 1, ii, 1) +
                            q(j + 1, k - 1, ii, 5) / q(j + 1, k - 1, ii, 1) -
                            q(j - 1, k - 1, ii, 5) / q(j - 1, k - 1, ii, 1));
                        t(jkv + j, 8) = 0.25 * (
                            (1.0 + bck(j + 1, ii, 2)) * qk0(j + 1, ii, 2, 3) -
                            bck(j + 1, ii, 2) * q(j + 1, k - 1, ii, 2) -
                            (1.0 + bck(j - 1, ii, 2)) * qk0(j - 1, ii, 2, 3) +
                            bck(j - 1, ii, 2) * q(j - 1, k - 1, ii, 2) +
                            q(j + 1, k - 1, ii, 2) - q(j - 1, k - 1, ii, 2));
                        t(jkv + j, 9) = 0.25 * (
                            (1.0 + bck(j + 1, ii, 2)) * qk0(j + 1, ii, 3, 3) -
                            bck(j + 1, ii, 2) * q(j + 1, k - 1, ii, 3) -
                            (1.0 + bck(j - 1, ii, 2)) * qk0(j - 1, ii, 3, 3) +
                            bck(j - 1, ii, 2) * q(j - 1, k - 1, ii, 3) +
                            q(j + 1, k - 1, ii, 3) - q(j - 1, k - 1, ii, 3));
                        t(jkv + j, 10) = 0.25 * (
                            (1.0 + bck(j + 1, ii, 2)) * qk0(j + 1, ii, 4, 3) -
                            bck(j + 1, ii, 2) * q(j + 1, k - 1, ii, 4) -
                            (1.0 + bck(j - 1, ii, 2)) * qk0(j - 1, ii, 4, 3) +
                            bck(j - 1, ii, 2) * q(j - 1, k - 1, ii, 4) +
                            q(j + 1, k - 1, ii, 4) - q(j - 1, k - 1, ii, 4));
                    }
                } // end j loop for k==kdim



            } else {
                // general case for k between 1 and kdim
                for (j = 1; j <= jdim1; j++) {
                    if (j == 1) {
                        t(jkv + j, 7) = 0.25 * gprgm1 * (
                            q(j + 1, k, ii, 5) / q(j + 1, k, ii, 1) -
                            (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 5, 1) / qj0(k, ii, 1, 1) +
                            bcj(k, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) +
                            q(j + 1, k - 1, ii, 5) / q(j + 1, k - 1, ii, 1) -
                            (1.0 + bcj(k - 1, ii, 1)) * qj0(k - 1, ii, 5, 1) / qj0(k - 1, ii, 1, 1) +
                            bcj(k - 1, ii, 1) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1));
                        t(jkv + j, 8) = 0.25 * (
                            q(j + 1, k, ii, 2) -
                            (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 2, 1) +
                            bcj(k, ii, 1) * q(j, k, ii, 2) +
                            q(j + 1, k - 1, ii, 2) -
                            (1.0 + bcj(k - 1, ii, 1)) * qj0(k - 1, ii, 2, 1) +
                            bcj(k - 1, ii, 1) * q(j, k - 1, ii, 2));
                        t(jkv + j, 9) = 0.25 * (
                            q(j + 1, k, ii, 3) -
                            (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 3, 1) +
                            bcj(k, ii, 1) * q(j, k, ii, 3) +
                            q(j + 1, k - 1, ii, 3) -
                            (1.0 + bcj(k - 1, ii, 1)) * qj0(k - 1, ii, 3, 1) +
                            bcj(k - 1, ii, 1) * q(j, k - 1, ii, 3));
                        t(jkv + j, 10) = 0.25 * (
                            q(j + 1, k, ii, 4) -
                            (1.0 + bcj(k, ii, 1)) * qj0(k, ii, 4, 1) +
                            bcj(k, ii, 1) * q(j, k, ii, 4) +
                            q(j + 1, k - 1, ii, 4) -
                            (1.0 + bcj(k - 1, ii, 1)) * qj0(k - 1, ii, 4, 1) +
                            bcj(k - 1, ii, 1) * q(j, k - 1, ii, 4));
                    } else if (j == jdim1) {
                        t(jkv + j, 7) = 0.25 * gprgm1 * (
                            (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 5, 3) / qj0(k, ii, 1, 3) -
                            bcj(k, ii, 2) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                            q(j - 1, k, ii, 5) / q(j - 1, k, ii, 1) +
                            (1.0 + bcj(k - 1, ii, 2)) * qj0(k - 1, ii, 5, 3) / qj0(k - 1, ii, 1, 3) -
                            bcj(k - 1, ii, 2) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1) -
                            q(j - 1, k - 1, ii, 5) / q(j - 1, k - 1, ii, 1));
                        t(jkv + j, 8) = 0.25 * (
                            (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 2, 3) -
                            bcj(k, ii, 2) * q(j, k, ii, 2) -
                            q(j - 1, k, ii, 2) +
                            (1.0 + bcj(k - 1, ii, 2)) * qj0(k - 1, ii, 2, 3) -
                            bcj(k - 1, ii, 2) * q(j, k - 1, ii, 2) -
                            q(j - 1, k - 1, ii, 2));
                        t(jkv + j, 9) = 0.25 * (
                            (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 3, 3) -
                            bcj(k, ii, 2) * q(j, k, ii, 3) -
                            q(j - 1, k, ii, 3) +
                            (1.0 + bcj(k - 1, ii, 2)) * qj0(k - 1, ii, 3, 3) -
                            bcj(k - 1, ii, 2) * q(j, k - 1, ii, 3) -
                            q(j - 1, k - 1, ii, 3));
                        t(jkv + j, 10) = 0.25 * (
                            (1.0 + bcj(k, ii, 2)) * qj0(k, ii, 4, 3) -
                            bcj(k, ii, 2) * q(j, k, ii, 4) -
                            q(j - 1, k, ii, 4) +
                            (1.0 + bcj(k - 1, ii, 2)) * qj0(k - 1, ii, 4, 3) -
                            bcj(k - 1, ii, 2) * q(j, k - 1, ii, 4) -
                            q(j - 1, k - 1, ii, 4));
                    } else {
                        // general case
                        t(jkv + j, 7) = 0.25 * gprgm1 * (
                            q(j + 1, k, ii, 5) / q(j + 1, k, ii, 1) -
                            q(j - 1, k, ii, 5) / q(j - 1, k, ii, 1) +
                            q(j + 1, k - 1, ii, 5) / q(j + 1, k - 1, ii, 1) -
                            q(j - 1, k - 1, ii, 5) / q(j - 1, k - 1, ii, 1));
                        t(jkv + j, 8) = 0.25 * (
                            q(j + 1, k, ii, 2) - q(j - 1, k, ii, 2) +
                            q(j + 1, k - 1, ii, 2) - q(j - 1, k - 1, ii, 2));
                        t(jkv + j, 9) = 0.25 * (
                            q(j + 1, k, ii, 3) - q(j - 1, k, ii, 3) +
                            q(j + 1, k - 1, ii, 3) - q(j - 1, k - 1, ii, 3));
                        t(jkv + j, 10) = 0.25 * (
                            q(j + 1, k, ii, 4) - q(j - 1, k, ii, 4) +
                            q(j + 1, k - 1, ii, 4) - q(j - 1, k - 1, ii, 4));
                    }
                } // end j loop for general k
            }
        } // end k loop
    } // end ipl loop (loop 230)



    // Lines 783-1163: gradients at cell interfaces (xie direction)
    // t(3): d(c*c/(gm1*pr))/d(xie), t(4): d(u)/d(xie), t(5): d(v)/d(xie), t(6): d(w)/d(xie)
    if (i2d == 0 && idim > 2) {
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            if (ii == 1) {
                // replace (ii-1) data by boundary data
                for (k = 1; k <= kdim; k++) {
                    jkv = (k - 1) * npl * jdim1 + (ipl - 1) * jdim1;
                    if (k == 1) {
                        for (j = 1; j <= jdim1; j++) {
                            t(jkv + j, 3) = 0.25 * gprgm1 * (
                                q(j, k, ii + 1, 5) / q(j, k, ii + 1, 1) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 5, 1) / qi0(j, k, 1, 1) +
                                bci(j, k, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) +
                                2.0 * (1.0 + bck(j, ii + 1, 1)) * qk0(j, ii + 1, 5, 1) / qk0(j, ii + 1, 1, 1) -
                                2.0 * bck(j, ii + 1, 1) * q(j, k, ii + 1, 5) / q(j, k, ii + 1, 1) -
                                2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 5, 1) / qk0(j, ii, 1, 1) +
                                2.0 * bck(j, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1));
                            t(jkv + j, 4) = 0.25 * (
                                q(j, k, ii + 1, 2) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 2, 1) +
                                bci(j, k, 1) * q(j, k, ii, 2) +
                                2.0 * (1.0 + bck(j, ii + 1, 1)) * qk0(j, ii + 1, 2, 1) -
                                2.0 * bck(j, ii + 1, 1) * q(j, k, ii + 1, 2) -
                                2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 2, 1) +
                                2.0 * bck(j, ii, 1) * q(j, k, ii, 2));
                            t(jkv + j, 5) = 0.25 * (
                                q(j, k, ii + 1, 3) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 3, 1) +
                                bci(j, k, 1) * q(j, k, ii, 3) +
                                2.0 * (1.0 + bck(j, ii + 1, 1)) * qk0(j, ii + 1, 3, 1) -
                                2.0 * bck(j, ii + 1, 1) * q(j, k, ii + 1, 3) -
                                2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 3, 1) +
                                2.0 * bck(j, ii, 1) * q(j, k, ii, 3));
                            t(jkv + j, 6) = 0.25 * (
                                q(j, k, ii + 1, 4) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 4, 1) +
                                bci(j, k, 1) * q(j, k, ii, 4) +
                                2.0 * (1.0 + bck(j, ii + 1, 1)) * qk0(j, ii + 1, 4, 1) -
                                2.0 * bck(j, ii + 1, 1) * q(j, k, ii + 1, 4) -
                                2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 4, 1) +
                                2.0 * bck(j, ii, 1) * q(j, k, ii, 4));
                        }
                    } else if (k == kdim) {
                        for (j = 1; j <= jdim1; j++) {
                            t(jkv + j, 3) = 0.25 * gprgm1 * (
                                2.0 * (1.0 + bck(j, ii + 1, 2)) * qk0(j, ii + 1, 5, 3) / qk0(j, ii + 1, 1, 3) -
                                2.0 * bck(j, ii + 1, 2) * q(j, k - 1, ii + 1, 5) / q(j, k - 1, ii + 1, 1) -
                                2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 5, 3) / qk0(j, ii, 1, 3) +
                                2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1) +
                                q(j, k - 1, ii + 1, 5) / q(j, k - 1, ii + 1, 1) -
                                (1.0 + bci(j, k - 1, 1)) * qi0(j, k - 1, 5, 1) / qi0(j, k - 1, 1, 1) +
                                bci(j, k - 1, 1) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1));
                            t(jkv + j, 4) = 0.25 * (
                                2.0 * (1.0 + bck(j, ii + 1, 2)) * qk0(j, ii + 1, 2, 3) -
                                2.0 * bck(j, ii + 1, 2) * q(j, k - 1, ii + 1, 2) -
                                2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 2, 3) +
                                2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 2) +
                                q(j, k - 1, ii + 1, 2) -
                                (1.0 + bci(j, k - 1, 1)) * qi0(j, k - 1, 2, 1) +
                                bci(j, k - 1, 1) * q(j, k - 1, ii, 2));
                            t(jkv + j, 5) = 0.25 * (
                                2.0 * (1.0 + bck(j, ii + 1, 2)) * qk0(j, ii + 1, 3, 3) -
                                2.0 * bck(j, ii + 1, 2) * q(j, k - 1, ii + 1, 3) -
                                2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 3, 3) +
                                2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 3) +
                                q(j, k - 1, ii + 1, 3) -
                                (1.0 + bci(j, k - 1, 1)) * qi0(j, k - 1, 3, 1) +
                                bci(j, k - 1, 1) * q(j, k - 1, ii, 3));
                            t(jkv + j, 6) = 0.25 * (
                                2.0 * (1.0 + bck(j, ii + 1, 2)) * qk0(j, ii + 1, 4, 3) -
                                2.0 * bck(j, ii + 1, 2) * q(j, k - 1, ii + 1, 4) -
                                2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 4, 3) +
                                2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 4) +
                                q(j, k - 1, ii + 1, 4) -
                                (1.0 + bci(j, k - 1, 1)) * qi0(j, k - 1, 4, 1) +
                                bci(j, k - 1, 1) * q(j, k - 1, ii, 4));
                        }
                    } else {
                        // general case for ii = 1
                        for (j = 1; j <= jdim1; j++) {
                            t(jkv + j, 3) = 0.25 * gprgm1 * (
                                q(j, k, ii + 1, 5) / q(j, k, ii + 1, 1) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 5, 1) / qi0(j, k, 1, 1) +
                                bci(j, k, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) +
                                q(j, k - 1, ii + 1, 5) / q(j, k - 1, ii + 1, 1) -
                                (1.0 + bci(j, k - 1, 1)) * qi0(j, k - 1, 5, 1) / qi0(j, k - 1, 1, 1) +
                                bci(j, k - 1, 1) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1));
                            t(jkv + j, 4) = 0.25 * (
                                q(j, k, ii + 1, 2) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 2, 1) +
                                bci(j, k, 1) * q(j, k, ii, 2) +
                                q(j, k - 1, ii + 1, 2) -
                                (1.0 + bci(j, k - 1, 1)) * qi0(j, k - 1, 2, 1) +
                                bci(j, k - 1, 1) * q(j, k - 1, ii, 2));
                            t(jkv + j, 5) = 0.25 * (
                                q(j, k, ii + 1, 3) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 3, 1) +
                                bci(j, k, 1) * q(j, k, ii, 3) +
                                q(j, k - 1, ii + 1, 3) -
                                (1.0 + bci(j, k - 1, 1)) * qi0(j, k - 1, 3, 1) +
                                bci(j, k - 1, 1) * q(j, k - 1, ii, 3));
                            t(jkv + j, 6) = 0.25 * (
                                q(j, k, ii + 1, 4) -
                                (1.0 + bci(j, k, 1)) * qi0(j, k, 4, 1) +
                                bci(j, k, 1) * q(j, k, ii, 4) +
                                q(j, k - 1, ii + 1, 4) -
                                (1.0 + bci(j, k - 1, 1)) * qi0(j, k - 1, 4, 1) +
                                bci(j, k - 1, 1) * q(j, k - 1, ii, 4));
                        }
                    }
                } // end k loop for ii==1



            } else if (ii == idim1) {
                // replace (ii+1) data by boundary data
                for (k = 1; k <= kdim; k++) {
                    jkv = (k - 1) * npl * jdim1 + (ipl - 1) * jdim1;
                    if (k == 1) {
                        for (j = 1; j <= jdim1; j++) {
                            t(jkv + j, 3) = 0.25 * gprgm1 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 5, 3) / qi0(j, k, 1, 3) -
                                bci(j, k, 2) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                                q(j, k, ii - 1, 5) / q(j, k, ii - 1, 1) +
                                2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 5, 1) / qk0(j, ii, 1, 1) -
                                2.0 * bck(j, ii, 1) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                                2.0 * (1.0 + bck(j, ii - 1, 1)) * qk0(j, ii - 1, 5, 1) / qk0(j, ii - 1, 1, 1) +
                                2.0 * bck(j, ii - 1, 1) * q(j, k, ii - 1, 5) / q(j, k, ii - 1, 1));
                            t(jkv + j, 4) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 2, 3) -
                                bci(j, k, 2) * q(j, k, ii, 2) -
                                q(j, k, ii - 1, 2) +
                                2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 2, 1) -
                                2.0 * bck(j, ii, 1) * q(j, k, ii, 2) -
                                2.0 * (1.0 + bck(j, ii - 1, 1)) * qk0(j, ii - 1, 2, 1) +
                                2.0 * bck(j, ii - 1, 1) * q(j, k, ii - 1, 2));
                            t(jkv + j, 5) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 3, 3) -
                                bci(j, k, 2) * q(j, k, ii, 3) -
                                q(j, k, ii - 1, 3) +
                                2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 3, 1) -
                                2.0 * bck(j, ii, 1) * q(j, k, ii, 3) -
                                2.0 * (1.0 + bck(j, ii - 1, 1)) * qk0(j, ii - 1, 3, 1) +
                                2.0 * bck(j, ii - 1, 1) * q(j, k, ii - 1, 3));
                            t(jkv + j, 6) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 4, 3) -
                                bci(j, k, 2) * q(j, k, ii, 4) -
                                q(j, k, ii - 1, 4) +
                                2.0 * (1.0 + bck(j, ii, 1)) * qk0(j, ii, 4, 1) -
                                2.0 * bck(j, ii, 1) * q(j, k, ii, 4) -
                                2.0 * (1.0 + bck(j, ii - 1, 1)) * qk0(j, ii - 1, 4, 1) +
                                2.0 * bck(j, ii - 1, 1) * q(j, k, ii - 1, 4));
                        }
                    } else if (k == kdim) {
                        for (j = 1; j <= jdim1; j++) {
                            t(jkv + j, 3) = 0.25 * gprgm1 * (
                                2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 5, 3) / qk0(j, ii, 1, 3) -
                                2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1) -
                                2.0 * (1.0 + bck(j, ii - 1, 2)) * qk0(j, ii - 1, 5, 3) / qk0(j, ii - 1, 1, 3) +
                                2.0 * bck(j, ii - 1, 2) * q(j, k - 1, ii - 1, 5) / q(j, k - 1, ii - 1, 1) +
                                (1.0 + bci(j, k - 1, 2)) * qi0(j, k - 1, 5, 3) / qi0(j, k - 1, 1, 3) -
                                bci(j, k - 1, 2) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1) -
                                q(j, k - 1, ii - 1, 5) / q(j, k - 1, ii - 1, 1));
                            t(jkv + j, 4) = 0.25 * (
                                2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 2, 3) -
                                2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 2) -
                                2.0 * (1.0 + bck(j, ii - 1, 2)) * qk0(j, ii - 1, 2, 3) +
                                2.0 * bck(j, ii - 1, 2) * q(j, k - 1, ii - 1, 2) +
                                (1.0 + bci(j, k - 1, 2)) * qi0(j, k - 1, 2, 3) -
                                bci(j, k - 1, 2) * q(j, k - 1, ii, 2) -
                                q(j, k - 1, ii - 1, 2));
                            t(jkv + j, 5) = 0.25 * (
                                2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 3, 3) -
                                2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 3) -
                                2.0 * (1.0 + bck(j, ii - 1, 2)) * qk0(j, ii - 1, 3, 3) +
                                2.0 * bck(j, ii - 1, 2) * q(j, k - 1, ii - 1, 3) +
                                (1.0 + bci(j, k - 1, 2)) * qi0(j, k - 1, 3, 3) -
                                bci(j, k - 1, 2) * q(j, k - 1, ii, 3) -
                                q(j, k - 1, ii - 1, 3));
                            t(jkv + j, 6) = 0.25 * (
                                2.0 * (1.0 + bck(j, ii, 2)) * qk0(j, ii, 4, 3) -
                                2.0 * bck(j, ii, 2) * q(j, k - 1, ii, 4) -
                                2.0 * (1.0 + bck(j, ii - 1, 2)) * qk0(j, ii - 1, 4, 3) +
                                2.0 * bck(j, ii - 1, 2) * q(j, k - 1, ii - 1, 4) +
                                (1.0 + bci(j, k - 1, 2)) * qi0(j, k - 1, 4, 3) -
                                bci(j, k - 1, 2) * q(j, k - 1, ii, 4) -
                                q(j, k - 1, ii - 1, 4));
                        }
                    } else {
                        // general case for ii = idim1
                        for (j = 1; j <= jdim1; j++) {
                            t(jkv + j, 3) = 0.25 * gprgm1 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 5, 3) / qi0(j, k, 1, 3) -
                                bci(j, k, 2) * q(j, k, ii, 5) / q(j, k, ii, 1) -
                                q(j, k, ii - 1, 5) / q(j, k, ii - 1, 1) +
                                (1.0 + bci(j, k - 1, 2)) * qi0(j, k - 1, 5, 3) / qi0(j, k - 1, 1, 3) -
                                bci(j, k - 1, 2) * q(j, k - 1, ii, 5) / q(j, k - 1, ii, 1) -
                                q(j, k - 1, ii - 1, 5) / q(j, k - 1, ii - 1, 1));
                            t(jkv + j, 4) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 2, 3) -
                                bci(j, k, 2) * q(j, k, ii, 2) -
                                q(j, k, ii - 1, 2) +
                                (1.0 + bci(j, k - 1, 2)) * qi0(j, k - 1, 2, 3) -
                                bci(j, k - 1, 2) * q(j, k - 1, ii, 2) -
                                q(j, k - 1, ii - 1, 2));
                            t(jkv + j, 5) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 3, 3) -
                                bci(j, k, 2) * q(j, k, ii, 3) -
                                q(j, k, ii - 1, 3) +
                                (1.0 + bci(j, k - 1, 2)) * qi0(j, k - 1, 3, 3) -
                                bci(j, k - 1, 2) * q(j, k - 1, ii, 3) -
                                q(j, k - 1, ii - 1, 3));
                            t(jkv + j, 6) = 0.25 * (
                                (1.0 + bci(j, k, 2)) * qi0(j, k, 4, 3) -
                                bci(j, k, 2) * q(j, k, ii, 4) -
                                q(j, k, ii - 1, 4) +
                                (1.0 + bci(j, k - 1, 2)) * qi0(j, k - 1, 4, 3) -
                                bci(j, k - 1, 2) * q(j, k - 1, ii, 4) -
                                q(j, k - 1, ii - 1, 4));
                        }
                    }
                } // end k loop for ii==idim1



            } else {
                // general case (ii != 1 and ii != idim1)
                for (k = 1; k <= kdim; k++) {
                    jkv = (k - 1) * npl * jdim1 + (ipl - 1) * jdim1;
                    if (k == 1) {
                        for (j = 1; j <= jdim1; j++) {
                            t(jkv + j, 3) = 0.25 * gprgm1 * (
                                q(j, k, ii + 1, 5) / q(j, k, ii + 1, 1) -
                                q(j, k, ii - 1, 5) / q(j, k, ii - 1, 1) +
                                (1.0 + bck(j, ii + 1, 1)) * qk0(j, ii + 1, 5, 1) / qk0(j, ii + 1, 1, 1) -
                                bck(j, ii + 1, 1) * q(j, k, ii + 1, 5) / q(j, k, ii + 1, 1) -
                                (1.0 + bck(j, ii - 1, 1)) * qk0(j, ii - 1, 5, 1) / qk0(j, ii - 1, 1, 1) +
                                bck(j, ii - 1, 1) * q(j, k, ii - 1, 5) / q(j, k, ii - 1, 1));
                            t(jkv + j, 4) = 0.25 * (
                                q(j, k, ii + 1, 2) - q(j, k, ii - 1, 2) +
                                (1.0 + bck(j, ii + 1, 1)) * qk0(j, ii + 1, 2, 1) -
                                bck(j, ii + 1, 1) * q(j, k, ii + 1, 2) -
                                (1.0 + bck(j, ii - 1, 1)) * qk0(j, ii - 1, 2, 1) +
                                bck(j, ii - 1, 1) * q(j, k, ii - 1, 2));
                            t(jkv + j, 5) = 0.25 * (
                                q(j, k, ii + 1, 3) - q(j, k, ii - 1, 3) +
                                (1.0 + bck(j, ii + 1, 1)) * qk0(j, ii + 1, 3, 1) -
                                bck(j, ii + 1, 1) * q(j, k, ii + 1, 3) -
                                (1.0 + bck(j, ii - 1, 1)) * qk0(j, ii - 1, 3, 1) +
                                bck(j, ii - 1, 1) * q(j, k, ii - 1, 3));
                            t(jkv + j, 6) = 0.25 * (
                                q(j, k, ii + 1, 4) - q(j, k, ii - 1, 4) +
                                (1.0 + bck(j, ii + 1, 1)) * qk0(j, ii + 1, 4, 1) -
                                bck(j, ii + 1, 1) * q(j, k, ii + 1, 4) -
                                (1.0 + bck(j, ii - 1, 1)) * qk0(j, ii - 1, 4, 1) +
                                bck(j, ii - 1, 1) * q(j, k, ii - 1, 4));
                        }
                    } else if (k == kdim) {
                        for (j = 1; j <= jdim1; j++) {
                            t(jkv + j, 3) = 0.25 * gprgm1 * (
                                (1.0 + bck(j, ii + 1, 2)) * qk0(j, ii + 1, 5, 3) / qk0(j, ii + 1, 1, 3) -
                                bck(j, ii + 1, 2) * q(j, k - 1, ii + 1, 5) / q(j, k - 1, ii + 1, 1) -
                                (1.0 + bck(j, ii - 1, 2)) * qk0(j, ii - 1, 5, 3) / qk0(j, ii - 1, 1, 3) +
                                bck(j, ii - 1, 2) * q(j, k - 1, ii - 1, 5) / q(j, k - 1, ii - 1, 1) +
                                q(j, k - 1, ii + 1, 5) / q(j, k - 1, ii + 1, 1) -
                                q(j, k - 1, ii - 1, 5) / q(j, k - 1, ii - 1, 1));
                            t(jkv + j, 4) = 0.25 * (
                                (1.0 + bck(j, ii + 1, 2)) * qk0(j, ii + 1, 2, 3) -
                                bck(j, ii + 1, 2) * q(j, k - 1, ii + 1, 2) -
                                (1.0 + bck(j, ii - 1, 2)) * qk0(j, ii - 1, 2, 3) +
                                bck(j, ii - 1, 2) * q(j, k - 1, ii - 1, 2) +
                                q(j, k - 1, ii + 1, 2) - q(j, k - 1, ii - 1, 2));
                            t(jkv + j, 5) = 0.25 * (
                                (1.0 + bck(j, ii + 1, 2)) * qk0(j, ii + 1, 3, 3) -
                                bck(j, ii + 1, 2) * q(j, k - 1, ii + 1, 3) -
                                (1.0 + bck(j, ii - 1, 2)) * qk0(j, ii - 1, 3, 3) +
                                bck(j, ii - 1, 2) * q(j, k - 1, ii - 1, 3) +
                                q(j, k - 1, ii + 1, 3) - q(j, k - 1, ii - 1, 3));
                            t(jkv + j, 6) = 0.25 * (
                                (1.0 + bck(j, ii + 1, 2)) * qk0(j, ii + 1, 4, 3) -
                                bck(j, ii + 1, 2) * q(j, k - 1, ii + 1, 4) -
                                (1.0 + bck(j, ii - 1, 2)) * qk0(j, ii - 1, 4, 3) +
                                bck(j, ii - 1, 2) * q(j, k - 1, ii - 1, 4) +
                                q(j, k - 1, ii + 1, 4) - q(j, k - 1, ii - 1, 4));
                        }
                    } else {
                        // general case
                        for (j = 1; j <= jdim1; j++) {
                            t(jkv + j, 3) = 0.25 * gprgm1 * (
                                q(j, k, ii + 1, 5) / q(j, k, ii + 1, 1) -
                                q(j, k, ii - 1, 5) / q(j, k, ii - 1, 1) +
                                q(j, k - 1, ii + 1, 5) / q(j, k - 1, ii + 1, 1) -
                                q(j, k - 1, ii - 1, 5) / q(j, k - 1, ii - 1, 1));
                            t(jkv + j, 4) = 0.25 * (
                                q(j, k, ii + 1, 2) - q(j, k, ii - 1, 2) +
                                q(j, k - 1, ii + 1, 2) - q(j, k - 1, ii - 1, 2));
                            t(jkv + j, 5) = 0.25 * (
                                q(j, k, ii + 1, 3) - q(j, k, ii - 1, 3) +
                                q(j, k - 1, ii + 1, 3) - q(j, k - 1, ii - 1, 3));
                            t(jkv + j, 6) = 0.25 * (
                                q(j, k, ii + 1, 4) - q(j, k, ii - 1, 4) +
                                q(j, k - 1, ii + 1, 4) - q(j, k - 1, ii - 1, 4));
                        }
                    }
                } // end k loop for general ii
            }
        } // end ipl loop (loop 380)
    } // end if (i2d == 0 && idim > 2)



    // Lines 1169-1199: t(24) turbulent viscosity at interfaces
    iviscc = ivisc[2];
    if (iviscc > 1) {
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            for (k = 1; k <= kdim; k++) {
                jkv = (k - 1) * npl * jdim1 + (ipl - 1) * jdim1;
                // k=0 interface
                if (k == 1) {
                    for (j = 1; j <= jdim1; j++) {
                        t(j + jkv, 24) = bck(j, ii, 1) * vk0(j, ii, 1, 1) +
                            (1.0 - bck(j, ii, 1)) * 0.5 * (vk0(j, ii, 1, 1) + vist3d(j, 1, ii));
                    }
                // k=kdim interface
                } else if (k == kdim) {
                    for (j = 1; j <= jdim1; j++) {
                        t(j + jkv, 24) = bck(j, ii, 2) * vk0(j, ii, 1, 3) +
                            (1.0 - bck(j, ii, 2)) * 0.5 * (vk0(j, ii, 1, 3) + vist3d(j, kdim1, ii));
                    }
                } else {
                    // interior interfaces
                    for (j = 1; j <= jdim1; j++) {
                        t(j + jkv, 24) = 0.5 * (vist3d(j, k, ii) + vist3d(j, k - 1, ii));
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

    // Lines 1201-1218: t(23) ratio of Prandtl numbers, t(14) final viscosity
    for (izz = 1; izz <= l0; izz++) {
        // ratio of turbulent to laminar viscosity
        t25 = t(izz, 24) / t(izz, 14);
        t24 = 1.e0 + t25 * (double)coef_eddy;
        // t(23): ratio of laminar Prandtl number to total Prandtl number
        t(izz, 23) = (1.e0 + prtr * t25) / t24;
        // t(14): (mach/re)*viscosity/J
        t(izz, 14) = xmre * t24 * t(izz, 14) / t(izz, 15);
    }

    // Line 1220: early return if iadv < 0
    if (iadv < 0) return;

    // Lines 1229-1248: t(11-13) velocity at cell interfaces
    // interior interfaces
    for (izz = 1; izz <= nn; izz++) {
        t(izz + js - 1, 11) = 0.5e0 * (t(izz + js - 1, 18) + t(izz, 18));
        t(izz + js - 1, 12) = 0.5e0 * (t(izz + js - 1, 19) + t(izz, 19));
        t(izz + js - 1, 13) = 0.5e0 * (t(izz + js - 1, 20) + t(izz, 20));
    }
    // interfaces at k=0 and k=kdim
    for (izz = 1; izz <= jv; izz++) {
        ab = 1.0 + wk0(izz, 15);
        bb = 1.0 - wk0(izz, 15);
        t(izz + n, 11) = 0.5 * (ab * wk0(izz, 18) + bb * t(izz + n - jv, 18));
        t(izz + n, 12) = 0.5 * (ab * wk0(izz, 19) + bb * t(izz + n - jv, 19));
        t(izz + n, 13) = 0.5 * (ab * wk0(izz, 20) + bb * t(izz + n - jv, 20));
        ab = 1.0 + wk0(izz, 5);
        bb = 1.0 - wk0(izz, 5);
        t(izz, 11) = 0.5 * (ab * wk0(izz, 8)  + bb * t(izz, 18));
        t(izz, 12) = 0.5 * (ab * wk0(izz, 9)  + bb * t(izz, 19));
        t(izz, 13) = 0.5 * (ab * wk0(izz, 10) + bb * t(izz, 20));
    }

    // Lines 1256-1310: Calculate viscous fluxes
    for (izz = 1; izz <= l0; izz++) {
        // form ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz
        // excluding contributions in k/zeta direction (already accounted for in thin-layer)
        ux_v = t(izz, 4) * t(izz, 28) + t(izz, 8) * t(izz, 31);
        uy   = t(izz, 4) * t(izz, 29) + t(izz, 8) * t(izz, 32);
        uz   = t(izz, 4) * t(izz, 30) + t(izz, 8) * t(izz, 2);

        vx   = t(izz, 5) * t(izz, 28) + t(izz, 9) * t(izz, 31);
        vy   = t(izz, 5) * t(izz, 29) + t(izz, 9) * t(izz, 32);
        vz   = t(izz, 5) * t(izz, 30) + t(izz, 9) * t(izz, 2);

        wx   = t(izz, 6) * t(izz, 28) + t(izz, 10) * t(izz, 31);
        wy   = t(izz, 6) * t(izz, 29) + t(izz, 10) * t(izz, 32);
        wz   = t(izz, 6) * t(izz, 30) + t(izz, 10) * t(izz, 2);

        tx   = t(izz, 3) * t(izz, 28) + t(izz, 7) * t(izz, 31);
        ty   = t(izz, 3) * t(izz, 29) + t(izz, 7) * t(izz, 32);
        tz   = t(izz, 3) * t(izz, 30) + t(izz, 7) * t(izz, 2);

        // Note lambda = -2/3 mu, lambda + 2 mu = 4/3 mu
        // u-momentum flux
        t(izz, 16) = t(izz, 14) * (t(izz, 25) * (4.0 / 3.0 * ux_v - 2.0 / 3.0 * (vy + wz))
                                  + t(izz, 26) * (uy + vx)
                                  + t(izz, 27) * (uz + wx));

        // v-momentum flux
        t(izz, 17) = t(izz, 14) * (t(izz, 25) * (vx + uy)
                                  + t(izz, 26) * (4.0 / 3.0 * vy - 2.0 / 3.0 * (ux_v + wz))
                                  + t(izz, 27) * (vz + wy));

        // w-momentum flux
        t(izz, 18) = t(izz, 14) * (t(izz, 25) * (wx + uz)
                                  + t(izz, 26) * (wy + vz)
                                  + t(izz, 27) * (4.0 / 3.0 * wz - 2.0 / 3.0 * (vy + ux_v)));

        // energy flux
        t(izz, 19) = t(izz, 14) * (
            t(izz, 25) * ((4.0 / 3.0 * ux_v - 2.0 / 3.0 * (vy + wz)) * t(izz, 11) +
                          (uy + vx) * t(izz, 12) + (uz + wx) * t(izz, 13)) +
            t(izz, 26) * ((vx + uy) * t(izz, 11) + (vz + wy) * t(izz, 13) +
                          (4.0 / 3.0 * vy - 2.0 / 3.0 * (ux_v + wz)) * t(izz, 12)) +
            t(izz, 27) * ((wx + uz) * t(izz, 11) + (wy + vz) * t(izz, 12) +
                          (4.0 / 3.0 * wz - 2.0 / 3.0 * (ux_v + vy)) * t(izz, 13)) +
            t(izz, 23) * (t(izz, 25) * tx + t(izz, 26) * ty + t(izz, 27) * tz));
    }

    // Lines 1312-1318: (-)viscous flux = Gv(k-1/2) - Gv(k+1/2) (momentum eqs.)
    for (izz = 1; izz <= n; izz++) {
        t(izz, 2) = -t(izz + js - 1, 16) + t(izz, 16);
        t(izz, 3) = -t(izz + js - 1, 17) + t(izz, 17);
        t(izz, 4) = -t(izz + js - 1, 18) + t(izz, 18);
        t(izz, 5) = -t(izz + js - 1, 19) + t(izz, 19);
    }

    // Lines 1324-1339: Calculate residuals
    for (ipl = 1; ipl <= npl; ipl++) {
        ii  = i + ipl - 1;
        for (k = 1; k <= kdim1; k++) {
            jkv = (k - 1) * jdim1 * npl + (ipl - 1) * jdim1 + 1;
            // for 2-D, t(3) should be identically zero
            if (i2d == 1) {
                for (izz = 1; izz <= jdim1; izz++) {
                    t(izz + jkv - 1, 3) = 0.0;
                }
            }
            for (l = 2; l <= 5; l++) {
                for (izz = 1; izz <= jdim1; izz++) {
                    res(izz, k, ii, l) = res(izz, k, ii, l) + t(izz + jkv - 1, l);
                }
            }
        }
    }
} // end hfluxv1

} // namespace hfluxv1_ns
