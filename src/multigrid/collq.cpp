// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "collq.h"
#include "runtime/common_blocks.h"
#include "fill.h"
#include <cstring>
#include <cstdio>
#include <algorithm>

namespace collq_ns {

void collq(FortranArray4DRef<double> q, FortranArray4DRef<double> qq,
           FortranArray3DRef<double> vol, FortranArray3DRef<double> volc,
           int& jdim, int& kdim, int& idim,
           int& jj2, int& kk2, int& ii2,
           FortranArray4DRef<double> res, FortranArray4DRef<double> qr,
           FortranArray4DRef<double> qw,
           FortranArray3DRef<double> vistf, FortranArray3DRef<double> vistc,
           FortranArray4DRef<double> tursavf, FortranArray4DRef<double> tursavc,
           int& nbl, FortranArray1DRef<int> nou,
           FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& nummem)
{
    // COMMON blocks
    int32_t*  ivisc   = cmn_reyue.ivisc;  // 0-based in C++: ivisc[0]=ivisc(1), etc.
    int32_t&  isklton = cmn_sklton.isklton;

    int nbl1 = nbl + 1;

    if (isklton > 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::memset(bou(nou(1), 1), ' ', 120);
        nou(1) = std::min(nou(1) + 1, ibufdim);
        // format 7: 1x,45hrestricting variables and residual from finer,6h block,i4,1x,16hto coarser block,i4
        std::snprintf(bou(nou(1), 1), 120,
            " restricting variables and residual from finer block%4d to coarser block%4d",
            nbl, nbl1);
    }

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;
    int jjl   = jj2 - 1;  (void)jjl;
    int kkl   = kk2 - 1;  (void)kkl;
    int iil   = ii2 - 1;
    int n     = jj2 * kk2;

    // Zero out qr: qr(izz,1,i,l) for i=1..iil, l=1..5, izz=1..n
    for (int i = 1; i <= iil; i++) {
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= n; izz++) {
                qr(izz, 1, i, l) = 0.0;
            }
        }
    }

    // Main restriction loop over 5 variables
    for (int nv_n = 1; nv_n <= 5; nv_n++) {
        int nv = jdim * kdim;

        // qw = q * vol for i=1..idim1
        for (int i = 1; i <= idim1; i++) {
            for (int izz = 1; izz <= nv; izz++) {
                qw(izz, 1, i, nv_n) = q(izz, 1, i, nv_n) * vol(izz, 1, i);
            }
        }

        // qw at idim = q * vol(idim1)
        for (int izz = 1; izz <= nv; izz++) {
            qw(izz, 1, idim, nv_n) = q(izz, 1, idim, nv_n) * vol(izz, 1, idim1);
        }

        if (idim > 2) {
            // 3D restriction
            int kk = 0;
            for (int k = 1; k <= kdim1; k += 2) {
                kk++;
                int jj = 0;
                for (int j = 1; j <= jdim1; j += 2) {
                    jj++;
                    int ii = 0;
                    for (int i = 1; i <= idim1; i += 2) {
                        ii++;
                        qq(jj, kk, ii, nv_n) =
                            (qw(j,   k,   i,   nv_n) + qw(j,   k,   i+1, nv_n)
                           + qw(j+1, k,   i,   nv_n) + qw(j+1, k,   i+1, nv_n)
                           + qw(j,   k+1, i,   nv_n) + qw(j,   k+1, i+1, nv_n)
                           + qw(j+1, k+1, i,   nv_n) + qw(j+1, k+1, i+1, nv_n))
                            / volc(jj, kk, ii);

                        qr(jj, kk, ii, nv_n) =
                            (res(j,   k,   i,   nv_n) + res(j,   k,   i+1, nv_n)
                           + res(j+1, k,   i,   nv_n) + res(j+1, k,   i+1, nv_n)
                           + res(j,   k+1, i,   nv_n) + res(j,   k+1, i+1, nv_n)
                           + res(j+1, k+1, i,   nv_n) + res(j+1, k+1, i+1, nv_n));

                        if ((ivisc[0] > 1 && nv_n == 1) ||
                            (ivisc[1] > 1 && nv_n == 1) ||
                            (ivisc[2] > 1 && nv_n == 1)) {
                            vistc(jj, kk, ii) =
                                0.125e0 * (vistf(j,   k,   i  ) + vistf(j,   k,   i+1)
                                         + vistf(j+1, k,   i  ) + vistf(j+1, k,   i+1)
                                         + vistf(j,   k+1, i  ) + vistf(j,   k+1, i+1)
                                         + vistf(j+1, k+1, i  ) + vistf(j+1, k+1, i+1));
                        }
                    }
                }
            }
        } else {
            // 2D restriction (idim==2, i=1 fixed)
            int ii = 1;
            int i  = 1;
            int kk = 0;
            for (int k = 1; k <= kdim1; k += 2) {
                kk++;
                int jj = 0;
                for (int j = 1; j <= jdim1; j += 2) {
                    jj++;
                    qq(jj, kk, ii, nv_n) =
                        (qw(j,   k,   i, nv_n) + qw(j+1, k,   i, nv_n)
                       + qw(j,   k+1, i, nv_n) + qw(j+1, k+1, i, nv_n))
                        / volc(jj, kk, ii);

                    qr(jj, kk, ii, nv_n) =
                        (res(j,   k,   i, nv_n) + res(j+1, k,   i, nv_n)
                       + res(j,   k+1, i, nv_n) + res(j+1, k+1, i, nv_n));

                    if ((ivisc[0] > 1 && nv_n == 1) ||
                        (ivisc[1] > 1 && nv_n == 1) ||
                        (ivisc[2] > 1 && nv_n == 1)) {
                        vistc(jj, kk, ii) =
                            0.25e0 * (vistf(j,   k,   i) + vistf(j+1, k,   i)
                                    + vistf(j,   k+1, i) + vistf(j+1, k+1, i));
                    }
                }
            }
        }
    } // end do 60 n=1,5

    // call fill(jj2,kk2,ii2,qq,5)
    int ll5 = 5;
    fill_ns::fill(jj2, kk2, ii2, qq, ll5);

    // Turbulence variables for RSMs (ivisc >= 70)
    if (ivisc[0] >= 70 || ivisc[1] >= 70 || ivisc[2] >= 70) {
        for (int n = 1; n <= nummem; n++) {
            if (idim > 2) {
                // 3D restriction of tursav
                int kk = 0;
                for (int k = 1; k <= kdim1; k += 2) {
                    kk++;
                    int jj = 0;
                    for (int j = 1; j <= jdim1; j += 2) {
                        jj++;
                        int ii = 0;
                        for (int i = 1; i <= idim1; i += 2) {
                            ii++;
                            tursavc(jj, kk, ii, n) =
                                0.125e0 * (tursavf(j,   k,   i,   n) + tursavf(j,   k,   i+1, n)
                                         + tursavf(j+1, k,   i,   n) + tursavf(j+1, k,   i+1, n)
                                         + tursavf(j,   k+1, i,   n) + tursavf(j,   k+1, i+1, n)
                                         + tursavf(j+1, k+1, i,   n) + tursavf(j+1, k+1, i+1, n));
                        }
                    }
                }
            } else {
                // 2D restriction of tursav
                int ii = 1;
                int i  = 1;
                int kk = 0;
                for (int k = 1; k <= kdim1; k += 2) {
                    kk++;
                    int jj = 0;
                    for (int j = 1; j <= jdim1; j += 2) {
                        jj++;
                        tursavc(jj, kk, ii, n) =
                            0.25e0 * (tursavf(j,   k,   i, n) + tursavf(j+1, k,   i, n)
                                    + tursavf(j,   k+1, i, n) + tursavf(j+1, k+1, i, n));
                    }
                }
            }
        }
    }

    return;
}

} // namespace collq_ns
