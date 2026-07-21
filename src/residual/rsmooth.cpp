// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rsmooth.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace rsmooth_ns {

void rsmooth(FortranArray1DRef<double> eps, int& idim, int& jdim, int& kdim, int& icall,
             FortranArray4DRef<double> dq, FortranArray2DRef<double> d,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim)
{
    // local scalars
    double smoopi, smoopj, smoopk;
    double a, t;
    int i, j, k;

    smoopi = eps(1);
    smoopj = eps(2);
    smoopk = eps(3);

    if (cmn_sklton.isklton == 1) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "smoothing residuals - eps=%12.4E%12.4E%12.4E",
            (float)smoopi, (float)smoopj, (float)smoopk);
    }

    // --------- smoothing in i direction
    if (idim > 2 && std::abs((float)smoopi) > 0.f) {
        for (k = 1; k <= kdim - 1; k++) {
            // forward pass: i=1
            for (j = 1; j <= jdim - 1; j++) {
                a           = smoopi;
                t           = 1.0 / (1.0 + a + a);
                d(j, 1)      = t * a;
                dq(j, k, 1, 1) = t * dq(j, k, 1, 1);
                dq(j, k, 1, 2) = t * dq(j, k, 1, 2);
                dq(j, k, 1, 3) = t * dq(j, k, 1, 3);
                dq(j, k, 1, 4) = t * dq(j, k, 1, 4);
                dq(j, k, 1, 5) = t * dq(j, k, 1, 5);
            }
            // forward pass: i=2..idim-1
            for (i = 2; i <= idim - 1; i++) {
                a = smoopi;
                for (j = 1; j <= jdim - 1; j++) {
                    t           = 1.0 / (1.0 + a + a - a * d(j, i - 1));
                    d(j, i)      = t * a;
                    dq(j, k, i, 1) = t * (dq(j, k, i, 1) + a * dq(j, k, i - 1, 1));
                    dq(j, k, i, 2) = t * (dq(j, k, i, 2) + a * dq(j, k, i - 1, 2));
                    dq(j, k, i, 3) = t * (dq(j, k, i, 3) + a * dq(j, k, i - 1, 3));
                    dq(j, k, i, 4) = t * (dq(j, k, i, 4) + a * dq(j, k, i - 1, 4));
                    dq(j, k, i, 5) = t * (dq(j, k, i, 5) + a * dq(j, k, i - 1, 5));
                }
            }
            // backward pass: i=idim-2..1
            for (i = idim - 2; i >= 1; i--) {
                for (j = 1; j <= jdim - 1; j++) {
                    dq(j, k, i, 1) = dq(j, k, i, 1) + d(j, i) * dq(j, k, i + 1, 1);
                    dq(j, k, i, 2) = dq(j, k, i, 2) + d(j, i) * dq(j, k, i + 1, 2);
                    dq(j, k, i, 3) = dq(j, k, i, 3) + d(j, i) * dq(j, k, i + 1, 3);
                    dq(j, k, i, 4) = dq(j, k, i, 4) + d(j, i) * dq(j, k, i + 1, 4);
                    dq(j, k, i, 5) = dq(j, k, i, 5) + d(j, i) * dq(j, k, i + 1, 5);
                }
            }
        }
    }

    // --------- smoothing in j direction
    if (jdim > 2 && std::abs((float)smoopj) > 0.f) {
        for (i = 1; i <= idim - 1; i++) {
            // forward pass: j=1
            a = smoopj;
            for (k = 1; k <= kdim - 1; k++) {
                t           = 1.0 / (1.0 + a + a);
                d(1, k)      = t * a;
                dq(1, k, i, 1) = t * dq(1, k, i, 1);
                dq(1, k, i, 2) = t * dq(1, k, i, 2);
                dq(1, k, i, 3) = t * dq(1, k, i, 3);
                dq(1, k, i, 4) = t * dq(1, k, i, 4);
                dq(1, k, i, 5) = t * dq(1, k, i, 5);
            }
            // forward pass: j=2..jdim-1
            for (j = 2; j <= jdim - 1; j++) {
                for (k = 1; k <= kdim - 1; k++) {
                    t           = 1.0 / (1.0 + a + a - a * d(j - 1, k));
                    d(j, k)      = t * a;
                    dq(j, k, i, 1) = t * (dq(j, k, i, 1) + a * dq(j - 1, k, i, 1));
                    dq(j, k, i, 2) = t * (dq(j, k, i, 2) + a * dq(j - 1, k, i, 2));
                    dq(j, k, i, 3) = t * (dq(j, k, i, 3) + a * dq(j - 1, k, i, 3));
                    dq(j, k, i, 4) = t * (dq(j, k, i, 4) + a * dq(j - 1, k, i, 4));
                    dq(j, k, i, 5) = t * (dq(j, k, i, 5) + a * dq(j - 1, k, i, 5));
                }
            }
            // backward pass: j=jdim-2..1
            for (j = jdim - 2; j >= 1; j--) {
                for (k = 1; k <= kdim - 1; k++) {
                    dq(j, k, i, 1) = dq(j, k, i, 1) + d(j, k) * dq(j + 1, k, i, 1);
                    dq(j, k, i, 2) = dq(j, k, i, 2) + d(j, k) * dq(j + 1, k, i, 2);
                    dq(j, k, i, 3) = dq(j, k, i, 3) + d(j, k) * dq(j + 1, k, i, 3);
                    dq(j, k, i, 4) = dq(j, k, i, 4) + d(j, k) * dq(j + 1, k, i, 4);
                    dq(j, k, i, 5) = dq(j, k, i, 5) + d(j, k) * dq(j + 1, k, i, 5);
                }
            }
        }
    }

    // --------- smoothing in k direction
    if (kdim > 2 && std::abs((float)smoopk) > 0.f) {
        for (i = 1; i <= idim - 1; i++) {
            // forward pass: k=1
            a = smoopk;
            for (j = 1; j <= jdim - 1; j++) {
                t           = 1.0 / (1.0 + a + a);
                d(j, 1)      = t * a;
                dq(j, 1, i, 1) = t * dq(j, 1, i, 1);
                dq(j, 1, i, 2) = t * dq(j, 1, i, 2);
                dq(j, 1, i, 3) = t * dq(j, 1, i, 3);
                dq(j, 1, i, 4) = t * dq(j, 1, i, 4);
                dq(j, 1, i, 5) = t * dq(j, 1, i, 5);
            }
            // forward pass: k=2..kdim-1
            for (k = 2; k <= kdim - 1; k++) {
                for (j = 1; j <= jdim - 1; j++) {
                    t           = 1.0 / (1.0 + a + a - a * d(j, k - 1));
                    d(j, k)      = t * a;
                    dq(j, k, i, 1) = t * (dq(j, k, i, 1) + a * dq(j, k - 1, i, 1));
                    dq(j, k, i, 2) = t * (dq(j, k, i, 2) + a * dq(j, k - 1, i, 2));
                    dq(j, k, i, 3) = t * (dq(j, k, i, 3) + a * dq(j, k - 1, i, 3));
                    dq(j, k, i, 4) = t * (dq(j, k, i, 4) + a * dq(j, k - 1, i, 4));
                    dq(j, k, i, 5) = t * (dq(j, k, i, 5) + a * dq(j, k - 1, i, 5));
                }
            }
            // backward pass: k=kdim-2..1
            for (k = kdim - 2; k >= 1; k--) {
                for (j = 1; j <= jdim - 1; j++) {
                    dq(j, k, i, 1) = dq(j, k, i, 1) + d(j, k) * dq(j, k + 1, i, 1);
                    dq(j, k, i, 2) = dq(j, k, i, 2) + d(j, k) * dq(j, k + 1, i, 2);
                    dq(j, k, i, 3) = dq(j, k, i, 3) + d(j, k) * dq(j, k + 1, i, 3);
                    dq(j, k, i, 4) = dq(j, k, i, 4) + d(j, k) * dq(j, k + 1, i, 4);
                    dq(j, k, i, 5) = dq(j, k, i, 5) + d(j, k) * dq(j, k + 1, i, 5);
                }
            }
        }
    }
}

} // namespace rsmooth_ns
