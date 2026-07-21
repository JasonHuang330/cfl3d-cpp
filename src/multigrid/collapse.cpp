// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "collapse.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

namespace collapse_ns {

void collapse(int& jdim, int& kdim, int& jmaxo, int& kmaxo,
              FortranArray2DRef<double> xo, FortranArray2DRef<double> yo,
              FortranArray2DRef<double> zo,
              FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
              int& nbuf, int& ibufdim)
{
    // common /tol/ epsc,epsc0,epsreen,eps
    // NOTE: 'eps' is the 4th field = cmn_tol.epscoll
    float& eps    = cmn_tol.epscoll;
    // common /sklt1/ isklt1
    int32_t& isklt1 = cmn_sklt1.isklt1;

    // Local variables
    float fact, epss1, epss2;
    float sx0, sy0, sz0;
    float xoo, yoo, zoo;
    float sx, sy, sz;
    float fact1;
    int k, j, kc, jc, j1, j2, k1, k2, ncount, m, jj, kk;

    fact  = 1.e5f;
    epss1 = fact * eps;

    epss2 = 1.e-01f;

    // ******************************
    // check along lines k=const
    // ******************************

    for (k = 1; k <= kmaxo; k++) {
        kc     = k;
        j1     = 1;
        ncount = 0;

    label_110:
        j2 = j1 + 1 + ncount;
        if (std::abs((float)(xo(j2,k) - xo(j1,k))) < (float)eps &&
            std::abs((float)(yo(j2,k) - yo(j1,k))) < (float)eps &&
            std::abs((float)(zo(j2,k) - zo(j1,k))) < (float)eps) {
            ncount = ncount + 1;
            if (j1 == 1 && j2 == jmaxo) goto label_109;
            goto label_110;
        }
        if (ncount > 0) {
            j2 = j2 - 1;
            if (isklt1 > 0) {
                nou(4) = std::min(nou(4) + 1, ibufdim);
                std::snprintf(bou(nou(4), 4), 120,
                    "      collapsed boundary on k =%3d  between j =%3d  and j =%3d",
                    kc, j1, j2);
            }
            // 1) collapsed segments only over part of the k=const line
            if (j1 > 1) {
                sx0 = (float)(xo(j1,k) - xo(j1-1,k));
                sy0 = (float)(yo(j1,k) - yo(j1-1,k));
                sz0 = (float)(zo(j1,k) - zo(j1-1,k));
            } else {
                sx0 = (float)(xo(j2+1,k) - xo(j2,k));
                sy0 = (float)(yo(j2+1,k) - yo(j2,k));
                sz0 = (float)(zo(j2+1,k) - zo(j2,k));
            }

            xoo = (float)xo(j1,k);
            yoo = (float)yo(j1,k);
            zoo = (float)zo(j1,k);

            if (j1 > 1 && j2 < jmaxo) {
                jc = j1 + ncount / 2;
            } else if (j1 == 1) {
                jc = j1;
            } else {
                jc = j2;
            }

            for (jj = j1; jj <= j2; jj++) {
                fact1     = (float)(jj - jc) * epss1;
                xo(jj,k) = (double)(xoo + fact1 * sx0);
                yo(jj,k) = (double)(yoo + fact1 * sy0);
                zo(jj,k) = (double)(zoo + fact1 * sz0);
            }
        }

        j1     = j1 + ncount + 1;
        ncount = 0;
        if (j1 < jmaxo - 1) goto label_110;

        goto label_111;

        // 2) collapsed segments over the entire k=const line, but k+1
        //    line is not collapsed (or k-1 for k=kmaxo) - polar singularity
    label_109:
        if (isklt1 > 0) {
            nou(4) = std::min(nou(4) + 1, ibufdim);
            std::snprintf(bou(nou(4), 4), 120,
                "      collapsed boundary on k =%3d  between j =%3d  and j =%3d",
                kc, j1, j2);
        }
        m = 1;
        if (k == kmaxo) m = -1;
        for (j = 1; j <= jmaxo; j++) {
            sx      = (float)(xo(j,k+m) - xo(j,k));
            sy      = (float)(yo(j,k+m) - yo(j,k));
            sz      = (float)(zo(j,k+m) - zo(j,k));
            xo(j,k) = xo(j,k) + (double)(sx * epss2);
            yo(j,k) = yo(j,k) + (double)(sy * epss2);
            zo(j,k) = zo(j,k) + (double)(sz * epss2);
        }

    label_111:
        ; // continue outer k loop
    }

    // ******************************
    // check along lines j=const
    // ******************************

    for (j = 1; j <= jmaxo; j++) {
        jc     = j;
        k1     = 1;
        ncount = 0;

    label_1101:
        k2 = k1 + 1 + ncount;
        if (std::abs((float)(xo(j,k2) - xo(j,k1))) < (float)eps &&
            std::abs((float)(yo(j,k2) - yo(j,k1))) < (float)eps &&
            std::abs((float)(zo(j,k2) - zo(j,k1))) < (float)eps) {
            ncount = ncount + 1;
            if (k1 == 1 && k2 == kmaxo) goto label_1009;
            goto label_1101;
        }
        if (ncount > 0) {
            k2 = k2 - 1;
            if (isklt1 > 0) {
                nou(4) = std::min(nou(4) + 1, ibufdim);
                std::snprintf(bou(nou(4), 4), 120,
                    "      collapsed boundary on j =%3d  between k =%3d  and k =%3d",
                    jc, k1, k2);
            }
            // 1) collapsed segments only over part of the j=const line
            if (k1 > 1) {
                sx0 = (float)(xo(j,k1) - xo(j,k1-1));
                sy0 = (float)(yo(j,k1) - yo(j,k1-1));
                sz0 = (float)(zo(j,k1) - zo(j,k1-1));
            } else {
                sx0 = (float)(xo(j,k2+1) - xo(j,k2));
                sy0 = (float)(yo(j,k2+1) - yo(j,k2));
                sz0 = (float)(zo(j,k2+1) - zo(j,k2));
            }

            xoo = (float)xo(j,k1);
            yoo = (float)yo(j,k1);
            zoo = (float)zo(j,k1);

            if (k1 > 1 && k2 < kmaxo) {
                kc = k1 + ncount / 2;
            } else if (k1 == 1) {
                kc = k1;
            } else {
                kc = k2;
            }

            for (kk = k1; kk <= k2; kk++) {
                fact1     = (float)(kk - kc) * epss1;
                xo(j,kk) = (double)(xoo + fact1 * sx0);
                yo(j,kk) = (double)(yoo + fact1 * sy0);
                zo(j,kk) = (double)(zoo + fact1 * sz0);
            }
        }

        k1     = k1 + ncount + 1;
        ncount = 0;
        if (k1 < kmaxo - 1) goto label_1101;

        goto label_1111;

        // 2) collapsed segments over the entire j=const line, but j+1
        //    line is not collapsed (or j-1 for j=jmaxo) - polar singularity
    label_1009:
        if (isklt1 > 0) {
            nou(4) = std::min(nou(4) + 1, ibufdim);
            std::snprintf(bou(nou(4), 4), 120,
                "      collapsed boundary on j =%3d  between k =%3d  and k =%3d",
                jc, k1, k2);
        }
        m = 1;
        if (j == jmaxo) m = -1;
        for (k = 1; k <= kmaxo; k++) {
            sx      = (float)(xo(j+m,k) - xo(j,k));
            sy      = (float)(yo(j+m,k) - yo(j,k));
            sz      = (float)(zo(j+m,k) - zo(j,k));
            xo(j,k) = xo(j,k) + (double)(sx * epss2);
            yo(j,k) = yo(j,k) + (double)(sy * epss2);
            zo(j,k) = zo(j,k) + (double)(sz * epss2);
        }

    label_1111:
        ; // continue outer j loop
    }

    return;
}

} // namespace collapse_ns
