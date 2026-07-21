// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "i2xk_d.h"
#include <algorithm>

namespace i2xk_d_ns {

void i2x(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
          int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf,
          FortranArray4DRef<double> qkf, FortranArray4DRef<double> qif,
          int& js, int& ks, int& is, int& je, int& ke, int& ie,
          int& nblc, int& ldim, int& nbl,
          FortranArray3DRef<double> bcjf, FortranArray3DRef<double> bckf,
          FortranArray3DRef<double> bcif, int& nface)
{
    (void)qjf; (void)qif; (void)bcjf; (void)bcif;
    i2xk_d(jc, kc, ic, qc, jf, kf, if_, qkf,
            js, ks, is, je, ke, ie, nblc, ldim, nbl, bckf, nface);
}

void i2xk_d(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
             int& jf, int& kf, int& if_, FortranArray4DRef<double> qkf,
             int& js, int& ks, int& is, int& je, int& ke, int& ie,
             int& nblc, int& ldim, int& nbl,
             FortranArray3DRef<double> bckf, int& nface)
{
    (void)kc; (void)kf; (void)ks; (void)ke; (void)nblc; (void)nbl;

    double q[4]; // 1-based: q[1], q[2], q[3]

    int jem = je - 1;
    int iem = ie - 1;

    double f1 = .75 * .75;
    double f2 = .75 * .25;
    double f4 = .25 * .25;

    //
    // k = constant planes
    //

    if (nface == 5) {
        //
        // interpolate left boundary
        //
        for (int l = 1; l <= ldim; l++) {
        for (int i = is; i <= iem; i++) {
        for (int j = js; j <= jem; j++) {
        for (int jl = 1; jl <= 2; jl++) {
            int jj = (j - js) * 2 + jl;
            int j2 = std::max(j - 1 + (jl - 1) * 2, 1);
            j2 = std::min(jc - 1, j2);
            for (int il = 1; il <= 2; il++) {
                int ii = (i - is) * 2 + il;
                int i2 = std::max(i - 1 + (il - 1) * 2, 1);
                i2 = std::min(ic - 1, i2);
                for (int k = 1; k <= 3; k++) {
                    q[k] = f1 * qc(j, k, i, l)
                         + f2 * (qc(j, k, i2, l) + qc(j2, k, i, l))
                         + f4 * qc(j2, k, i2, l);
                }
                qkf(jj, ii, l, 1) = .25 * q[1] + .75 * q[2];
                qkf(jj, ii, l, 2) = .75 * q[2] + .25 * q[3];
                bckf(jj, ii, 1) = 0.0;
            }
        }
        }
        }
        }
    }

    if (nface == 6) {
        //
        // interpolate right boundary
        //
        for (int l = 1; l <= ldim; l++) {
        for (int i = is; i <= iem; i++) {
        for (int j = js; j <= jem; j++) {
        for (int jl = 1; jl <= 2; jl++) {
            int jj = (j - js) * 2 + jl;
            int j2 = std::max(j - 1 + (jl - 1) * 2, 1);
            j2 = std::min(jc - 1, j2);
            for (int il = 1; il <= 2; il++) {
                int ii = (i - is) * 2 + il;
                int i2 = std::max(i - 1 + (il - 1) * 2, 1);
                i2 = std::min(ic - 1, i2);
                for (int k = 1; k <= 3; k++) {
                    q[k] = f1 * qc(j, k, i, l)
                         + f2 * (qc(j, k, i2, l) + qc(j2, k, i, l))
                         + f4 * qc(j2, k, i2, l);
                }
                qkf(jj, ii, l, 3) = .25 * q[1] + .75 * q[2];
                qkf(jj, ii, l, 4) = .75 * q[2] + .25 * q[3];
                bckf(jj, ii, 2) = 0.0;
            }
        }
        }
        }
        }
    }

    //
    // **for safety**
    //
    for (int m = 1; m <= 4; m++) {
    for (int l = 1; l <= ldim; l++) {
        for (int i = 1; i <= if_ - 1; i++) {
            qkf(jf, i, l, m) = qkf(jf - 1, i, l, m);
        }
    }
    }
}

} // namespace i2xk_d_ns
