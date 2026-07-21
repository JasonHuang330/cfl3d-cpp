// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "i2xsk_d.h"

namespace i2xsk_d_ns {

void i2xsk_d(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
             int& jf, int& kf, int& if_, FortranArray4DRef<double> qkf,
             int& js, int& ks, int& is, int& je, int& ke, int& ie,
             int& nblc, int& ldim, int& nbl,
             FortranArray3DRef<double> bckf, int& nface)
{
    int jem, iem;
    int l, ii, i, j, jl, jj, j2, k, m;
    double f1, f2;
    double q[4]; // 1-based: q[1], q[2], q[3]

    jem = je - 1;
    iem = ie - 1;

    f1 = .75;
    f2 = .25;

    // k = constant planes

    if (nface == 5) {
        // interpolate left boundary
        for (l = 1; l <= ldim; l++) {
            ii = 0;
            for (i = is; i <= iem; i++) {
                ii = ii + 1;
                for (j = js; j <= jem; j++) {
                    for (jl = 1; jl <= 2; jl++) {
                        jj = (j - js) * 2 + jl;
                        j2 = std::max(j - 1 + (jl - 1) * 2, 1);
                        j2 = std::min(jc - 1, j2);
                        for (k = 1; k <= 3; k++) {
                            q[k] = f1 * qc(j, k, i, l) + f2 * qc(j2, k, i, l);
                        }
                        qkf(jj, ii, l, 1) = .25 * q[1] + .75 * q[2];
                        qkf(jj, ii, l, 2) = .75 * q[2] + .25 * q[3];
                        bckf(jj, ii, 1) = 0.0;
                    }
                }
            }
        }
    }

    if (nface == 6) {
        // interpolate right boundary
        for (l = 1; l <= ldim; l++) {
            ii = 0;
            for (i = is; i <= iem; i++) {
                ii = ii + 1;
                for (j = js; j <= jem; j++) {
                    for (jl = 1; jl <= 2; jl++) {
                        jj = (j - js) * 2 + jl;
                        j2 = std::max(j - 1 + (jl - 1) * 2, 1);
                        j2 = std::min(jc - 1, j2);
                        for (k = 1; k <= 3; k++) {
                            q[k] = f1 * qc(j, k, i, l) + f2 * qc(j2, k, i, l);
                        }
                        qkf(jj, ii, l, 3) = .25 * q[1] + .75 * q[2];
                        qkf(jj, ii, l, 4) = .75 * q[2] + .25 * q[3];
                        bckf(jj, ii, 2) = 0.0;
                    }
                }
            }
        }
    }

    // **for safety**
    for (m = 1; m <= 4; m++) {
        for (l = 1; l <= ldim; l++) {
            for (i = 1; i <= if_ - 1; i++) {
                qkf(jf, i, l, m) = qkf(jf - 1, i, l, m);
            }
        }
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
    (void)qjf; (void)qif; (void)bcjf; (void)bcif;
    i2xsk_d(jc, kc, ic, qc, jf, kf, if_, qkf,
            js, ks, is, je, ke, ie, nblc, ldim, nbl, bckf, nface);
}

void i2xs(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
          int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf,
          FortranArray4DRef<double> qkf, FortranArray4DRef<double> qif,
          int& js, int& ks, int& is, int& je, int& ke, int& ie,
          int& nblc, int& ldim, int& nbl,
          FortranArray3DRef<double> bcjf, FortranArray3DRef<double> bckf,
          FortranArray3DRef<double> bcif, int& nface)
{
    (void)qjf; (void)qif; (void)bcjf; (void)bcif;
    i2xsk_d(jc, kc, ic, qc, jf, kf, if_, qkf,
            js, ks, is, je, ke, ie, nblc, ldim, nbl, bckf, nface);
}

} // namespace i2xsk_d_ns
