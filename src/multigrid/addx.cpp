// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "addx.h"
#include <cstdio>
#include "termn8.h"
#include "rsmooth.h"
#include "fill.h"
#include "ccomplex.h"
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <cmath>
#include <algorithm>

namespace addx_ns {

void addx(FortranArray4DRef<double> q, FortranArray4DRef<double> qq,
          int& jdim, int& kdim, int& idim,
          int& jj2, int& kk2, int& ii2,
          FortranArray4DRef<double> q1, FortranArray4DRef<double> dq,
          FortranArray4DRef<double> wq, FortranArray4DRef<double> wqj,
          int& nbl, FortranArray3DRef<double> blank,
          FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
          int& nbuf, int& ibufdim, int& ll, int& myid)
{
    // Access COMMON blocks
    int& mode  = cmn_mgrd.mode;
    int& issc  = cmn_mgv.issc;
    // epsssc is float[3] in cmn_mgv; rsmooth expects FortranArray1DRef<double>
    // Build a double version for the call
    FortranArray1D<double> epsssc_d(3);
    epsssc_d(1) = (double)cmn_mgv.epsssc[0];
    epsssc_d(2) = (double)cmn_mgv.epsssc[1];
    epsssc_d(3) = (double)cmn_mgv.epsssc[2];

    // Note: mode=1 must ONLY be used for ll=5 (primitive variables)
    if (mode != 0 && ll != 5) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
                      "mode must = 0 when ll .ne. 5 in addx");
        int ierrflg = -1;
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    // jdim,kdim,idim  finer mesh indices
    // jj2,kk2,ii2   coarser mesh indices
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;
    int jdim2 = jdim - 2;
    int kdim2 = kdim - 2;
    int idim2 = idim - 2;
    int jjl   = jj2 - 1;
    int kkl   = kk2 - 1;
    int iil   = ii2 - 1;

    // Initialize dq to 0
    int nplq = std::min(idim1, 999000 / (jdim * kdim));
    int npl  = nplq;
    for (int i = 1; i <= idim1; i += nplq) {
        npl = nplq;
        if (i + npl - 1 > idim1) npl = idim1 + 1 - i;
        int nnpl = npl * jdim * kdim;
        for (int n = 1; n <= ll; n++) {
            for (int izz = 1; izz <= nnpl; izz++) {
                dq(izz, 1, i, n) = 0.0;
            }
        }
    }

    int nbl1 = nbl - 1;
    (void)nbl1;

    if (mode == 0) {
        // Copy qq -> wq
        nplq = std::min(iil, 999000 / (jj2 * kk2));
        npl  = nplq;
        for (int i = 1; i <= iil; i += nplq) {
            npl = nplq;
            if (i + npl - 1 > iil) npl = iil + 1 - i;
            int nnpl = npl * jj2 * kk2 - jj2 - 1;
            for (int n = 1; n <= ll; n++) {
                for (int izz = 1; izz <= nnpl; izz++) {
                    wq(izz, 1, i, n) = qq(izz, 1, i, n);
                }
            }
        }
    } else {
        // Compute correction qq - q1 -> wq
        nplq = std::min(iil, 999000 / (jj2 * kk2));
        npl  = nplq;
        for (int i = 1; i <= iil; i += nplq) {
            npl = nplq;
            if (i + npl - 1 > iil) npl = iil + 1 - i;
            int nnpl = npl * jj2 * kk2 - jj2 - 1;
            for (int n = 1; n <= ll; n++) {
                for (int izz = 1; izz <= nnpl; izz++) {
                    wq(izz, 1, i, n) = qq(izz, 1, i, n) - q1(izz, 1, i, n);
                }
            }
        }

        // Smooth corrections
        if (issc > 0) {
            FortranArray1DRef<double> epsssc_ref = epsssc_d.ref();
            // rsmooth(epsssc, ii2, jj2, kk2, ii2, wq, wqj, nou, bou, nbuf, ibufdim)
            // signature: rsmooth(eps, idim, jdim, kdim, icall, dq, d, nou, bou, nbuf, ibufdim)
            // rsmooth declares d(jdim,idim+kdim); here jdim=jj2, idim=ii2, kdim=kk2,
            // so d must be dimensioned (jj2, ii2+kk2) -- the leading dim is jj2, NOT
            // jdim*kk2*ii2*ll (a wrong leading dim sends d(j,i>=2) far out of bounds).
            // The underlying wqj (jdim,kk2,ii2,ll) has ample room for jj2*(ii2+kk2).
            FortranArray2DRef<double> wqj_2d(&wqj(1, 1, 1, 1), jj2, ii2 + kk2);
            rsmooth_ns::rsmooth(epsssc_ref, ii2, jj2, kk2, ii2,
                                wq, wqj_2d, nou, bou, nbuf, ibufdim);
        }
    }

    // Interpolate in j
    int jjl1 = jjl - 1;
    int kkl1 = kkl - 1;
    for (int n = 1; n <= ll; n++) {
        for (int i = 1; i <= iil; i++) {
            for (int k = 1; k <= kkl; k++) {
                wqj(1,      k, i, n) = wq(1,   k, i, n);
                wqj(jdim1,  k, i, n) = wq(jjl, k, i, n);
                int jj = 1;
                for (int j = 1; j <= jjl1; j++) {
                    jj = jj + 1;
                    wqj(jj, k, i, n) = 0.75 * wq(j, k, i, n) + 0.25 * wq(j+1, k, i, n);
                    jj = jj + 1;
                    wqj(jj, k, i, n) = 0.25 * wq(j, k, i, n) + 0.75 * wq(j+1, k, i, n);
                }
            }
        }

        // Interpolate in k
        int ii = idim1 - iil;
        for (int i = 1; i <= iil; i++) {
            ii = ii + 1;
            for (int j = 1; j <= jdim1; j++) {
                dq(j, 1, ii, n) = wqj(j, 1, i, n);
            }
            for (int izz = 1; izz <= jdim1; izz++) {
                dq(izz, kdim1, ii, n) = wqj(izz, kkl, i, n);
            }
            int kk = 1;
            for (int k = 1; k <= kkl1; k++) {
                kk = kk + 1;
                for (int izz = 1; izz <= jdim1; izz++) {
                    dq(izz, kk, ii, n) = 0.75 * wqj(izz, k, i, n) + 0.25 * wqj(izz, k+1, i, n);
                }
                kk = kk + 1;
                for (int izz = 1; izz <= jdim1; izz++) {
                    dq(izz, kk, ii, n) = 0.25 * wqj(izz, k, i, n) + 0.75 * wqj(izz, k+1, i, n);
                }
            }
        }

        // Interpolate in i
        if (idim1 > 1) {
            int is = idim1 - iil + 1;
            int np = jdim * kdim1 - 1;

            // i=1
            for (int izz = 1; izz <= np; izz++) {
                dq(izz, 1, 1, n) = dq(izz, 1, is, n);
            }
            int ii2_local = 1;
            for (int i = is; i <= idim2; i++) {
                ii2_local = ii2_local + 1;
                for (int izz = 1; izz <= np; izz++) {
                    dq(izz, 1, ii2_local, n) = 0.75 * dq(izz, 1, i, n) + 0.25 * dq(izz, 1, i+1, n);
                }
                ii2_local = ii2_local + 1;
                for (int izz = 1; izz <= np; izz++) {
                    dq(izz, 1, ii2_local, n) = 0.25 * dq(izz, 1, i, n) + 0.75 * dq(izz, 1, i+1, n);
                }
            }
        }
    } // end do 80 n=1,ll


    if (mode == 0) {
        nplq = std::min(idim1, 999000 / (jdim * kdim));
        npl  = nplq;
        for (int i = 1; i <= idim1; i += nplq) {
            npl = nplq;
            if (i + npl - 1 > idim1) npl = idim1 + 1 - i;
            int nnpl = npl * jdim * kdim - jdim - 1;
            for (int n = 1; n <= ll; n++) {
                for (int izz = 1; izz <= nnpl; izz++) {
                    q(izz, 1, i, n) = dq(izz, 1, i, n);
                }
            }
        }
    } else {
        nplq = std::min(idim1, 999000 / (jdim * kdim));
        npl  = nplq;
        for (int i = 1; i <= idim1; i += nplq) {
            npl = nplq;
            if (i + npl - 1 > idim1) npl = idim1 + 1 - i;
            int nnpl = npl * jdim * kdim - jdim - 1;

            // Update density and pressure to ensure positivity
            // "cut-off" point is determined by alpq
            // minimum value of density is equal to (1/phiq)
            double alpq = -0.2;
            double phiq = 1.0 / 0.5;
            double betq = 1.0 + alpq * phiq;

            for (int izz = 1; izz <= nnpl; izz++) {
                double t1 = dq(izz, 1, i, 1) / q(izz, 1, i, 1);
                double t2 = dq(izz, 1, i, 1) / (betq + std::abs(t1) * phiq);
                bool cond1 = ((float)t1 < (float)alpq);
                dq(izz, 1, i, 1) = cond1 ? t2 : dq(izz, 1, i, 1);

                t1 = dq(izz, 1, i, 5) / q(izz, 1, i, 5);
                t2 = dq(izz, 1, i, 5) / (betq + std::abs(t1) * phiq);
                bool cond5 = ((float)t1 < (float)alpq);
                dq(izz, 1, i, 5) = cond5 ? t2 : dq(izz, 1, i, 5);
            }

            // Update primitive variables
            for (int izz = 1; izz <= nnpl; izz++) {
                q(izz, 1, i, 1) = q(izz, 1, i, 1) + blank(izz, 1, i) * dq(izz, 1, i, 1);
                q(izz, 1, i, 2) = q(izz, 1, i, 2) + blank(izz, 1, i) * dq(izz, 1, i, 2);
                q(izz, 1, i, 3) = q(izz, 1, i, 3) + blank(izz, 1, i) * dq(izz, 1, i, 3);
                q(izz, 1, i, 4) = q(izz, 1, i, 4) + blank(izz, 1, i) * dq(izz, 1, i, 4);
                q(izz, 1, i, 5) = q(izz, 1, i, 5) + blank(izz, 1, i) * dq(izz, 1, i, 5);
            }
        }
    }

    fill_ns::fill(jdim, kdim, idim, q, ll);
    return;
}

} // namespace addx_ns
