// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "i2xsj_d.h"

namespace i2xsj_d_ns {

void i2xsj_d(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
              int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf,
              int& js, int& ks, int& is, int& je, int& ke, int& ie,
              int& nblc, int& ldim, int& nbl,
              FortranArray3DRef<double> bcjf, int& nface)
{
    int kem, iem, l, ii, i, k, kl, kk, k2, j, m;
    double q[4];
    double f1, f2;

    kem = ke - 1;
    iem = ie - 1;

    f1 = .75;
    f2 = .25;

    //
    // j = constant planes
    //

    if (nface == 3) {
        //
        // interpolate left boundary
        //
        for (l = 1; l <= ldim; l++) {
            ii = 0;
            for (i = is; i <= iem; i++) {
                ii = ii + 1;
                for (k = ks; k <= kem; k++) {
                    for (kl = 1; kl <= 2; kl++) {
                        kk = (k - ks) * 2 + kl;
                        k2 = k - 1 + (kl - 1) * 2;
                        if (k2 < 1) k2 = 1;
                        if (k2 > kc - 1) k2 = kc - 1;
                        for (j = 1; j <= 3; j++) {
                            q[j] = f1 * qc(j, k, i, l) + f2 * qc(j, k2, i, l);
                        }
                        qjf(kk, ii, l, 1) = .25 * q[1] + .75 * q[2];
                        qjf(kk, ii, l, 2) = .75 * q[2] + .25 * q[3];
                        bcjf(kk, ii, 1) = 0.0;
                    }
                }
            }
        }
    }

    if (nface == 4) {
        //
        // interpolate right boundary
        //
        for (l = 1; l <= ldim; l++) {
            ii = 0;
            for (i = is; i <= iem; i++) {
                ii = ii + 1;
                for (k = ks; k <= kem; k++) {
                    for (kl = 1; kl <= 2; kl++) {
                        kk = (k - ks) * 2 + kl;
                        k2 = k - 1 + (kl - 1) * 2;
                        if (k2 < 1) k2 = 1;
                        if (k2 > kc - 1) k2 = kc - 1;
                        for (j = 1; j <= 3; j++) {
                            q[j] = f1 * qc(j, k, i, l) + f2 * qc(j, k2, i, l);
                        }
                        qjf(kk, ii, l, 3) = .25 * q[1] + .75 * q[2];
                        qjf(kk, ii, l, 4) = .75 * q[2] + .25 * q[3];
                        bcjf(kk, ii, 2) = 0.0;
                    }
                }
            }
        }
    }

    //
    // **for safety**
    //
    for (m = 1; m <= 4; m++) {
        for (l = 1; l <= ldim; l++) {
            for (i = 1; i <= if_ - 1; i++) {
                qjf(kf, i, l, m) = qjf(kf - 1, i, l, m);
            }
        }
    }

    return;
}

void i2x(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
          int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf,
          FortranArray4DRef<double> qkf, FortranArray4DRef<double> qif,
          int& js, int& ks, int& is, int& je, int& ke, int& ie,
          int& nblc, int& ldim, int& nbl,
          FortranArray3DRef<double> bcjf, FortranArray3DRef<double> bckf,
          FortranArray3DRef<double> bcif, int& nface)
{
    i2xsj_d(jc, kc, ic, qc, jf, kf, if_, qjf,
             js, ks, is, je, ke, ie, nblc, ldim, nbl, bcjf, nface);
}

void i2xs(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
           int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf,
           FortranArray4DRef<double> qkf, FortranArray4DRef<double> qif,
           int& js, int& ks, int& is, int& je, int& ke, int& ie,
           int& nblc, int& ldim, int& nbl,
           FortranArray3DRef<double> bcjf, FortranArray3DRef<double> bckf,
           FortranArray3DRef<double> bcif, int& nface)
{
    i2xsj_d(jc, kc, ic, qc, jf, kf, if_, qjf,
             js, ks, is, je, ke, ie, nblc, ldim, nbl, bcjf, nface);
}

} // namespace i2xsj_d_ns
