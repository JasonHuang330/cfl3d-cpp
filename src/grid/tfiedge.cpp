// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "tfiedge.h"
#include "termn8.h"
#include "runtime/common_blocks.h"
#include <cmath>
#include <cstring>
#include <algorithm>

namespace tfiedge_ns {

void tfiedge(int& idim, int& jdim, int& kdim,
             FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
             int& i1, int& i2, int& j1, int& j2, int& k1, int& k2,
             FortranArray3DRef<double> arci, FortranArray3DRef<double> arcj, FortranArray3DRef<double> arck,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
             int& myid, int& nbl)
{
    // common /zero/ iexp
    int& iexp = cmn_zero.iexp;

    // local variables
    double psi1, psi2;
    double eta, denom;
    int i, j, k, nvar;

    // tolerance for switch to linear blending function
    // (10.**(-iexp) is machine zero)
    double tol = std::max(1.e-07, std::pow(10.0, (double)(-iexp + 1)));

    // check that one and only one of the pairs i1,i2 j1,j2 k1,k2 varies
    nvar = 0;
    if (i1 != i2) nvar = nvar + 1;
    if (j1 != j2) nvar = nvar + 1;
    if (k1 != k2) nvar = nvar + 1;
    if (nvar == 0 || nvar > 1) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
                      " stopping...input error to subroutine tfiedge");
        int minus1 = -1;
        termn8_ns::termn8(myid, minus1, ibufdim, nbuf, bou, nou);
        return;
    }
    if (i1 != i2) {
        j = j1;
        k = k1;
        denom = arci(j, k, i2) - arci(j, k, i1);
        for (i = i1; i <= i2; i++) {
            if ((float)denom < (float)tol) {
                eta = 0.0;
            } else {
                eta = (arci(j, k, i) - arci(j, k, i1)) / denom;
            }
            psi1 = eta;
            psi2 = 1.0 - eta;
            x(j, k, i) = psi2 * x(j, k, i1) + psi1 * x(j, k, i2);
            y(j, k, i) = psi2 * y(j, k, i1) + psi1 * y(j, k, i2);
            z(j, k, i) = psi2 * z(j, k, i1) + psi1 * z(j, k, i2);
        }
    }

    if (j1 != j2) {
        i = i1;
        k = k1;
        denom = arcj(j2, k, i) - arcj(j1, k, i);
        for (j = j1; j <= j2; j++) {
            if ((float)denom < (float)tol) {
                eta = 0.0;
            } else {
                eta = (arcj(j, k, i) - arcj(j1, k, i)) / denom;
            }
            psi1 = eta;
            psi2 = 1.0 - eta;
            x(j, k, i) = psi2 * x(j1, k, i) + psi1 * x(j2, k, i);
            y(j, k, i) = psi2 * y(j1, k, i) + psi1 * y(j2, k, i);
            z(j, k, i) = psi2 * z(j1, k, i) + psi1 * z(j2, k, i);
        }
    }

    if (k1 != k2) {
        i = i1;
        j = j1;
        denom = arck(j, k2, i) - arck(j, k1, i);
        for (k = k1; k <= k2; k++) {
            if ((float)denom < (float)tol) {
                eta = 0.0;
            } else {
                eta = (arck(j, k, i) - arck(j, k1, i)) / denom;
            }
            psi1 = eta;
            psi2 = 1.0 - eta;
            x(j, k, i) = psi2 * x(j, k1, i) + psi1 * x(j, k2, i);
            y(j, k, i) = psi2 * y(j, k1, i) + psi1 * y(j, k2, i);
            z(j, k, i) = psi2 * z(j, k1, i) + psi1 * z(j, k2, i);
        }
    }
}

} // namespace tfiedge_ns
