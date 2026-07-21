// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "i2x.h"
#include <algorithm>

namespace i2x_ns {

void i2x(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
          int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf,
          FortranArray4DRef<double> qkf, FortranArray4DRef<double> qif,
          int& js, int& ks, int& is, int& je, int& ke, int& ie,
          int& nblc, int& ldim, int& nbl,
          FortranArray3DRef<double> bcjf, FortranArray3DRef<double> bckf,
          FortranArray3DRef<double> bcif, int& nface)
{
    // local variables
    int kem, jem, iem;
    int l, i, k, kl, il, jl;
    int kk, k2, ii, i2, jj, j2;
    int j, m;
    double q[4]; // 1-based: q[1], q[2], q[3]

    double f1, f2, f4;

    kem = ke - 1;
    jem = je - 1;
    iem = ie - 1;

    f1 = .75 * .75;
    f2 = .75 * .25;
    f4 = .25 * .25;

    //
    // j = constant planes
    //

    if (nface == 3) {
        //
        // interpolate left boundary
        //
        for (l = 1; l <= ldim; l++)
        for (i = is; i <= iem; i++)
        for (k = ks; k <= kem; k++)
        for (kl = 1; kl <= 2; kl++) {
            kk = (k - ks) * 2 + kl;
            k2 = std::max(k - 1 + (kl - 1) * 2, 1);
            k2 = std::min(kc - 1, k2);
            for (il = 1; il <= 2; il++) {
                ii = (i - is) * 2 + il;
                i2 = std::max(i - 1 + (il - 1) * 2, 1);
                i2 = std::min(ic - 1, i2);
                for (jl = 1; jl <= 3; jl++) {
                    j = js + 1 - jl;
                    q[jl] = f1 * qc(j, k, i, l)
                           + f2 * (qc(j, k, i2, l) + qc(j, k2, i, l))
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
        // interpolate right boundary
        //
        for (l = 1; l <= ldim; l++)
        for (i = is; i <= iem; i++)
        for (k = ks; k <= kem; k++)
        for (kl = 1; kl <= 2; kl++) {
            kk = (k - ks) * 2 + kl;
            k2 = std::max(k - 1 + (kl - 1) * 2, 1);
            k2 = std::min(kc - 1, k2);
            for (il = 1; il <= 2; il++) {
                ii = (i - is) * 2 + il;
                i2 = std::max(i - 1 + (il - 1) * 2, 1);
                i2 = std::min(ic - 1, i2);
                for (jl = 1; jl <= 3; jl++) {
                    j = je + jl - 2;
                    q[jl] = f1 * qc(j, k, i, l)
                           + f2 * (qc(j, k, i2, l) + qc(j, k2, i, l))
                           + f4 * qc(j, k2, i2, l);
                }
                qjf(kk, ii, l, 3) = .25 * q[1] + .75 * q[2];
                qjf(kk, ii, l, 4) = .75 * q[2] + .25 * q[3];
                bcjf(kk, ii, 2) = 0.0;
            }
        }
    }

    //
    // k = constant planes
    //

    if (nface == 5) {
        //
        // interpolate left boundary
        //
        for (l = 1; l <= ldim; l++)
        for (i = is; i <= iem; i++)
        for (j = js; j <= jem; j++)
        for (jl = 1; jl <= 2; jl++) {
            jj = (j - js) * 2 + jl;
            j2 = std::max(j - 1 + (jl - 1) * 2, 1);
            j2 = std::min(jc - 1, j2);
            for (il = 1; il <= 2; il++) {
                ii = (i - is) * 2 + il;
                i2 = std::max(i - 1 + (il - 1) * 2, 1);
                i2 = std::min(ic - 1, i2);
                for (kl = 1; kl <= 3; kl++) {
                    k = ks + 1 - kl;
                    q[kl] = f1 * qc(j, k, i, l)
                           + f2 * (qc(j, k, i2, l) + qc(j2, k, i, l))
                           + f4 * qc(j2, k, i2, l);
                }
                qkf(jj, ii, l, 1) = .25 * q[1] + .75 * q[2];
                qkf(jj, ii, l, 2) = .75 * q[2] + .25 * q[3];
                bckf(jj, ii, 1) = 0.0;
            }
        }
    }

    if (nface == 6) {
        //
        // interpolate right boundary
        //
        for (l = 1; l <= ldim; l++)
        for (i = is; i <= iem; i++)
        for (j = js; j <= jem; j++)
        for (jl = 1; jl <= 2; jl++) {
            jj = (j - js) * 2 + jl;
            j2 = std::max(j - 1 + (jl - 1) * 2, 1);
            j2 = std::min(jc - 1, j2);
            for (il = 1; il <= 2; il++) {
                ii = (i - is) * 2 + il;
                i2 = std::max(i - 1 + (il - 1) * 2, 1);
                i2 = std::min(ic - 1, i2);
                for (kl = 1; kl <= 3; kl++) {
                    k = ke + kl - 2;
                    q[kl] = f1 * qc(j, k, i, l)
                           + f2 * (qc(j, k, i2, l) + qc(j2, k, i, l))
                           + f4 * qc(j2, k, i2, l);
                }
                qkf(jj, ii, l, 3) = .25 * q[1] + .75 * q[2];
                qkf(jj, ii, l, 4) = .75 * q[2] + .25 * q[3];
                bckf(jj, ii, 2) = 0.0;
            }
        }
    }

    //
    // i = constant plane
    //

    if (nface == 1) {
        //
        // interpolate left boundary
        //
        for (l = 1; l <= ldim; l++)
        for (k = ks; k <= kem; k++)
        for (j = js; j <= jem; j++)
        for (jl = 1; jl <= 2; jl++) {
            jj = (j - js) * 2 + jl;
            j2 = std::max(j - 1 + (jl - 1) * 2, 1);
            j2 = std::min(jc - 1, j2);
            for (kl = 1; kl <= 2; kl++) {
                kk = (k - ks) * 2 + kl;
                k2 = std::max(k - 1 + (kl - 1) * 2, 1);
                k2 = std::min(kc - 1, k2);
                for (il = 1; il <= 3; il++) {
                    i = is + 1 - il;
                    q[il] = f1 * qc(j, k, i, l)
                           + f2 * (qc(j, k2, i, l) + qc(j2, k, i, l))
                           + f4 * qc(j2, k2, i, l);
                }
                qif(jj, kk, l, 1) = .25 * q[1] + .75 * q[2];
                qif(jj, kk, l, 2) = .75 * q[2] + .25 * q[3];
                bcif(jj, kk, 1) = 0.0;
            }
        }
    }

    if (nface == 2) {
        //
        // interpolate right boundary
        //
        for (l = 1; l <= ldim; l++)
        for (k = ks; k <= kem; k++)
        for (j = js; j <= jem; j++)
        for (jl = 1; jl <= 2; jl++) {
            jj = (j - js) * 2 + jl;
            j2 = std::max(j - 1 + (jl - 1) * 2, 1);
            j2 = std::min(jc - 1, j2);
            for (kl = 1; kl <= 2; kl++) {
                kk = (k - ks) * 2 + kl;
                k2 = std::max(k - 1 + (kl - 1) * 2, 1);
                k2 = std::min(kc - 1, k2);
                for (il = 1; il <= 3; il++) {
                    i = ie + il - 2;
                    q[il] = f1 * qc(j, k, i, l)
                           + f2 * (qc(j, k2, i, l) + qc(j2, k, i, l))
                           + f4 * qc(j2, k2, i, l);
                }
                qif(jj, kk, l, 3) = .25 * q[1] + .75 * q[2];
                qif(jj, kk, l, 4) = .75 * q[2] + .25 * q[3];
                bcif(jj, kk, 2) = 0.0;
            }
        }
    }

    //
    // **for safety**
    //
    for (m = 1; m <= 4; m++)
    for (l = 1; l <= ldim; l++) {
        for (i = 1; i <= if_ - 1; i++) {
            qjf(kf, i, l, m) = qjf(kf - 1, i, l, m);
            qkf(jf, i, l, m) = qkf(jf - 1, i, l, m);
        }
        for (k = 1; k <= kf - 1; k++) {
            qif(jf, k, l, m) = qif(jf - 1, k, l, m);
        }
        for (j = 1; j <= jf - 1; j++) {
            qif(j, kf, l, m) = qif(j, kf - 1, l, m);
        }
    }
}

} // namespace i2x_ns
