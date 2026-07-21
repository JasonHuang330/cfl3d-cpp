// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "tfiface.h"
#include "termn8.h"
#include <cmath>
#include <algorithm>

namespace tfiface_ns {

void tfiface(int& idim, int& jdim, int& kdim,
             FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
             FortranArray3DRef<double> x1, FortranArray3DRef<double> y1, FortranArray3DRef<double> z1,
             FortranArray3DRef<double> x2, FortranArray3DRef<double> y2, FortranArray3DRef<double> z2,
             FortranArray3DRef<double> arci, FortranArray3DRef<double> arcj, FortranArray3DRef<double> arck,
             int& i1, int& i2, int& j1, int& j2, int& k1, int& k2,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& myid)
{
    // local arrays (1-indexed, size 2)
    double phi[3], psi[3], omg[3];  // 1-based: use indices 1 and 2

    // tolerance for switch to linear blending function
    // (10.**(-iexp) is machine zero)
    double tol = std::max(1.e-07, std::pow(10.0, (double)(-cmn_zero.iexp + 1)));

    // -------------------------------------------------------------------------
    // Case: i1 == i2  (face in j-k plane at fixed i)
    // -------------------------------------------------------------------------
    if (i1 == i2) {

        int i = i1;

        // check that j and k ranges span a logically 2d face
        if (j1 == j2 || k1 == k2) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "stopping in tfiface...error in j,k range for face i =%4d", i1);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            return;
        }

        // eta projector  -  1st step
        for (int j = j1; j <= j2; j++) {
            for (int k = k1; k <= k2; k++) {
                double denom = arcj(j2, k, i) - arcj(j1, k, i);
                double eta;
                if ((float)denom < (float)tol) {
                    eta = 0.0;
                } else {
                    eta = (arcj(j, k, i) - arcj(j1, k, i)) / denom;
                }
                psi[1] = eta;
                psi[2] = 1.0 - eta;
                x1(j, k, i) = psi[2] * x(j1, k, i) + psi[1] * x(j2, k, i);
                y1(j, k, i) = psi[2] * y(j1, k, i) + psi[1] * y(j2, k, i);
                z1(j, k, i) = psi[2] * z(j1, k, i) + psi[1] * z(j2, k, i);
            }
        }

        // zeta projector  -  2nd step
        for (int k = k1; k <= k2; k++) {
            for (int j = j1; j <= j2; j++) {
                double denom = arck(j, k2, i) - arck(j, k1, i);
                double zeta;
                if ((float)denom < (float)tol) {
                    zeta = 0.0;
                } else {
                    zeta = (arck(j, k, i) - arck(j, k1, i)) / denom;
                }
                omg[1] = zeta;
                omg[2] = 1.0 - zeta;
                x2(j, k, i) = omg[2] * (x(j, k1, i) - x1(j, k1, i))
                            + omg[1] * (x(j, k2, i) - x1(j, k2, i));
                y2(j, k, i) = omg[2] * (y(j, k1, i) - y1(j, k1, i))
                            + omg[1] * (y(j, k2, i) - y1(j, k2, i));
                z2(j, k, i) = omg[2] * (z(j, k1, i) - z1(j, k1, i))
                            + omg[1] * (z(j, k2, i) - z1(j, k2, i));
            }
        }

        // summation of eta and zeta projectors  -  3rd step
        for (int j = j1; j <= j2; j++) {
            for (int k = k1; k <= k2; k++) {
                x(j, k, i) = x1(j, k, i) + x2(j, k, i);
                y(j, k, i) = y1(j, k, i) + y2(j, k, i);
                z(j, k, i) = z1(j, k, i) + z2(j, k, i);
            }
        }

    }

    // -------------------------------------------------------------------------
    // Case: j1 == j2  (face in i-k plane at fixed j)
    // -------------------------------------------------------------------------
    if (j1 == j2) {

        int j = j1;

        // check that i and k ranges span a logically 2d face
        if (i1 == i2 || k1 == k2) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "stopping in tfiface...error in i,k range for face j =%4d", j1);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            return;
        }

        // xi projector  -  1st step
        for (int i = i1; i <= i2; i++) {
            for (int k = k1; k <= k2; k++) {
                double denom = arci(j, k, i2) - arci(j, k, i1);
                double xi;
                if ((float)denom < (float)tol) {
                    xi = 0.0;
                } else {
                    xi = (arci(j, k, i) - arci(j, k, i1)) / denom;
                }
                phi[1] = xi;
                phi[2] = 1.0 - xi;
                x1(j, k, i) = phi[2] * x(j, k, i1) + phi[1] * x(j, k, i2);
                y1(j, k, i) = phi[2] * y(j, k, i1) + phi[1] * y(j, k, i2);
                z1(j, k, i) = phi[2] * z(j, k, i1) + phi[1] * z(j, k, i2);
            }
        }

        // zeta projector  -  2nd step
        for (int k = k1; k <= k2; k++) {
            for (int i = i1; i <= i2; i++) {
                double denom = arck(j, k2, i) - arck(j, k1, i);
                double zeta;
                if ((float)denom < (float)tol) {
                    zeta = 0.0;
                } else {
                    zeta = (arck(j, k, i) - arck(j, k1, i)) / denom;
                }
                omg[1] = zeta;
                omg[2] = 1.0 - zeta;
                x2(j, k, i) = omg[2] * (x(j, k1, i) - x1(j, k1, i))
                            + omg[1] * (x(j, k2, i) - x1(j, k2, i));
                y2(j, k, i) = omg[2] * (y(j, k1, i) - y1(j, k1, i))
                            + omg[1] * (y(j, k2, i) - y1(j, k2, i));
                z2(j, k, i) = omg[2] * (z(j, k1, i) - z1(j, k1, i))
                            + omg[1] * (z(j, k2, i) - z1(j, k2, i));
            }
        }

        // summation of xi and zeta projectors  -  3rd step
        for (int i = i1; i <= i2; i++) {
            for (int k = k1; k <= k2; k++) {
                x(j, k, i) = x1(j, k, i) + x2(j, k, i);
                y(j, k, i) = y1(j, k, i) + y2(j, k, i);
                z(j, k, i) = z1(j, k, i) + z2(j, k, i);
            }
        }

    }

    // -------------------------------------------------------------------------
    // Case: k1 == k2  (face in i-j plane at fixed k)
    // -------------------------------------------------------------------------
    if (k1 == k2) {

        int k = k1;

        // check that i and j ranges span a logically 2d face
        if (i1 == i2 || j1 == j2) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "stopping in tfiface...error in i,j range for face k =%4d", k1);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            return;
        }

        // xi projector  -  1st step
        for (int i = i1; i <= i2; i++) {
            for (int j = j1; j <= j2; j++) {
                double denom = arci(j, k, i2) - arci(j, k, i1);
                double xi;
                if ((float)denom < (float)tol) {
                    xi = 0.0;
                } else {
                    xi = (arci(j, k, i) - arci(j, k, i1)) / denom;
                }
                phi[1] = xi;
                phi[2] = 1.0 - xi;
                x1(j, k, i) = phi[2] * x(j, k, i1) + phi[1] * x(j, k, i2);
                y1(j, k, i) = phi[2] * y(j, k, i1) + phi[1] * y(j, k, i2);
                z1(j, k, i) = phi[2] * z(j, k, i1) + phi[1] * z(j, k, i2);
            }
        }

        // eta projector  -  2nd step
        for (int j = j1; j <= j2; j++) {
            for (int i = i1; i <= i2; i++) {
                double denom = arcj(j2, k, i) - arcj(j1, k, i);
                double eta;
                if ((float)denom < (float)tol) {
                    eta = 0.0;
                } else {
                    eta = (arcj(j, k, i) - arcj(j1, k, i)) / denom;
                }
                psi[1] = eta;
                psi[2] = 1.0 - eta;
                x2(j, k, i) = psi[2] * (x(j1, k, i) - x1(j1, k, i))
                            + psi[1] * (x(j2, k, i) - x1(j2, k, i));
                y2(j, k, i) = psi[2] * (y(j1, k, i) - y1(j1, k, i))
                            + psi[1] * (y(j2, k, i) - y1(j2, k, i));
                z2(j, k, i) = psi[2] * (z(j1, k, i) - z1(j1, k, i))
                            + psi[1] * (z(j2, k, i) - z1(j2, k, i));
            }
        }

        // summation of xi and eta projectors  -  3rd step
        for (int i = i1; i <= i2; i++) {
            for (int j = j1; j <= j2; j++) {
                x(j, k, i) = x1(j, k, i) + x2(j, k, i);
                y(j, k, i) = y1(j, k, i) + y2(j, k, i);
                z(j, k, i) = z1(j, k, i) + z2(j, k, i);
            }
        }

    }

    return;
}

} // namespace tfiface_ns
