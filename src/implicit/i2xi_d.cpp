// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "i2xi_d.h"

namespace i2xi_d_ns {

void i2x(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
          int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf,
          FortranArray4DRef<double> qkf, FortranArray4DRef<double> qif,
          int& js, int& ks, int& is, int& je, int& ke, int& ie,
          int& nblc, int& ldim, int& nbl,
          FortranArray3DRef<double> bcjf, FortranArray3DRef<double> bckf,
          FortranArray3DRef<double> bcif, int& nface)
{
    (void)qjf; (void)qkf; (void)bcjf; (void)bckf;
    i2xi_d(jc, kc, ic, qc, jf, kf, if_, qif,
           js, ks, is, je, ke, ie, nblc, ldim, nbl, bcif, nface);
}

void i2xi_d(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
            int& jf, int& kf, int& if_, FortranArray4DRef<double> qif,
            int& js, int& ks, int& is, int& je, int& ke, int& ie,
            int& nblc, int& ldim, int& nbl,
            FortranArray3DRef<double> bcif, int& nface)
{
    (void)ic; (void)if_; (void)is; (void)ie; (void)nblc; (void)nbl;

    int kem, jem;
    double f1, f2, f4;
    double q[4]; // 1-based: use indices 1,2,3

    kem = ke - 1;
    jem = je - 1;

    f1 = .75 * .75;
    f2 = .75 * .25;
    f4 = .25 * .25;

    // i = constant plane

    if (nface == 1) {
        // interpolate left boundary
        for (int l = 1; l <= ldim; l++) {
        for (int k = ks; k <= kem; k++) {
        for (int j = js; j <= jem; j++) {
        for (int jl = 1; jl <= 2; jl++) {
            int jj = (j - js) * 2 + jl;
            int j2 = std::max(j - 1 + (jl - 1) * 2, 1);
            j2 = std::min(jc - 1, j2);
            for (int kl = 1; kl <= 2; kl++) {
                int kk = (k - ks) * 2 + kl;
                int k2 = std::max(k - 1 + (kl - 1) * 2, 1);
                k2 = std::min(kc - 1, k2);
                for (int i = 1; i <= 3; i++) {
                    q[i] = f1 * qc(j, k, i, l)
                         + f2 * (qc(j, k2, i, l) + qc(j2, k, i, l))
                         + f4 * qc(j2, k2, i, l);
                }
                qif(jj, kk, l, 1) = .25 * q[1] + .75 * q[2];
                qif(jj, kk, l, 2) = .75 * q[2] + .25 * q[3];
                bcif(jj, kk, 1) = 0.0;
            }
        }
        }
        }
        }
    }

    if (nface == 2) {
        // interpolate right boundary
        for (int l = 1; l <= ldim; l++) {
        for (int k = ks; k <= kem; k++) {
        for (int j = js; j <= jem; j++) {
        for (int jl = 1; jl <= 2; jl++) {
            int jj = (j - js) * 2 + jl;
            int j2 = std::max(j - 1 + (jl - 1) * 2, 1);
            j2 = std::min(jc - 1, j2);
            for (int kl = 1; kl <= 2; kl++) {
                int kk = (k - ks) * 2 + kl;
                int k2 = std::max(k - 1 + (kl - 1) * 2, 1);
                k2 = std::min(kc - 1, k2);
                for (int i = 1; i <= 3; i++) {
                    q[i] = f1 * qc(j, k, i, l)
                         + f2 * (qc(j, k2, i, l) + qc(j2, k, i, l))
                         + f4 * qc(j2, k2, i, l);
                }
                qif(jj, kk, l, 3) = .25 * q[1] + .75 * q[2];
                qif(jj, kk, l, 4) = .75 * q[2] + .25 * q[3];
                bcif(jj, kk, 2) = 0.0;
            }
        }
        }
        }
        }
    }

    // **for safety**
    for (int m = 1; m <= 4; m++) {
    for (int l = 1; l <= ldim; l++) {
        for (int k = 1; k <= kf - 1; k++) {
            qif(jf, k, l, m) = qif(jf - 1, k, l, m);
        }
        for (int j = 1; j <= jf - 1; j++) {
            qif(j, kf, l, m) = qif(j, kf - 1, l, m);
        }
    }
    }

    return;
}

} // namespace i2xi_d_ns
