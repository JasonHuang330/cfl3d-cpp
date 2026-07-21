// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "i2xj_d.h"

namespace i2xj_d_ns {

void i2xj_d(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
             int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf,
             int& js, int& ks, int& is, int& je, int& ke, int& ie,
             int& nblc, int& ldim, int& nbl, FortranArray3DRef<double> bcjf,
             int& nface)
{
    double q[4]; // 1-based: q[1], q[2], q[3]

    int kem = ke - 1;
    int iem = ie - 1;

    double f1 = .75 * .75;
    double f2 = .75 * .25;
    double f4 = .25 * .25;

    //
    //     j = constant planes
    //
    if (nface == 3) {
        //
        //     interpolate left boundary
        //
        for (int l = 1; l <= ldim; l++)
        for (int i = is; i <= iem; i++)
        for (int k = ks; k <= kem; k++)
        for (int kl = 1; kl <= 2; kl++) {
            int kk = (k - ks) * 2 + kl;
            int k2 = std::max(k - 1 + (kl - 1) * 2, 1);
            k2 = std::min(kc - 1, k2);
            for (int il = 1; il <= 2; il++) {
                int ii = (i - is) * 2 + il;
                int i2 = std::max(i - 1 + (il - 1) * 2, 1);
                i2 = std::min(ic - 1, i2);
                for (int j = 1; j <= 3; j++) {
                    q[j] = f1 * qc(j, k,  i,  l)
                         + f2 * (qc(j, k,  i2, l) + qc(j, k2, i,  l))
                         + f4 * qc(j, k2, i2, l);
                }
                qjf(kk, ii, l, 1) = .25 * q[1] + .75 * q[2];
                qjf(kk, ii, l, 2) = .75 * q[2] + .25 * q[3];
                bcjf(kk, ii, 1) = 0.0;
            }
        }
    }

    if (nface == 4) {
        //
        //     interpolate right boundary
        //
        for (int l = 1; l <= ldim; l++)
        for (int i = is; i <= iem; i++)
        for (int k = ks; k <= kem; k++)
        for (int kl = 1; kl <= 2; kl++) {
            int kk = (k - ks) * 2 + kl;
            int k2 = std::max(k - 1 + (kl - 1) * 2, 1);
            k2 = std::min(kc - 1, k2);
            for (int il = 1; il <= 2; il++) {
                int ii = (i - is) * 2 + il;
                int i2 = std::max(i - 1 + (il - 1) * 2, 1);
                i2 = std::min(ic - 1, i2);
                for (int j = 1; j <= 3; j++) {
                    q[j] = f1 * qc(j, k,  i,  l)
                         + f2 * (qc(j, k,  i2, l) + qc(j, k2, i,  l))
                         + f4 * qc(j, k2, i2, l);
                }
                qjf(kk, ii, l, 3) = .25 * q[1] + .75 * q[2];
                qjf(kk, ii, l, 4) = .75 * q[2] + .25 * q[3];
                bcjf(kk, ii, 2) = 0.0;
            }
        }
    }

    //
    //     **for safety**
    //
    for (int m = 1; m <= 4; m++)
    for (int l = 1; l <= ldim; l++)
    for (int i = 1; i <= if_ - 1; i++) {
        qjf(kf, i, l, m) = qjf(kf - 1, i, l, m);
    }
}

void i2x(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
         int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf,
         FortranArray4DRef<double> qkf, FortranArray4DRef<double> qif,
         int& js, int& ks, int& is, int& je, int& ke, int& ie,
         int& nblc, int& ldim, int& nbl,
         FortranArray3DRef<double> bcjf, FortranArray3DRef<double> bckf,
         FortranArray3DRef<double> bcif, int& nface)
{
    (void)qkf;
    (void)qif;
    (void)bckf;
    (void)bcif;
    i2xj_d(jc, kc, ic, qc, jf, kf, if_, qjf,
            js, ks, is, je, ke, ie, nblc, ldim, nbl, bcjf, nface);
}

} // namespace i2xj_d_ns
