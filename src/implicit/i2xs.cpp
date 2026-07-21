// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "i2xs.h"
#include "i2x.h"
#include <algorithm>

namespace i2xs_ns {

void i2x(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
          int& jf, int& kf, int& if_,
          FortranArray4DRef<double> qjf, FortranArray4DRef<double> qkf,
          FortranArray4DRef<double> qif,
          int& js, int& ks, int& is, int& je, int& ke, int& ie,
          int& nblc, int& ldim, int& nbl,
          FortranArray3DRef<double> bcjf, FortranArray3DRef<double> bckf,
          FortranArray3DRef<double> bcif, int& nface)
{
    i2x_ns::i2x(jc, kc, ic, qc, jf, kf, if_, qjf, qkf, qif,
                 js, ks, is, je, ke, ie, nblc, ldim, nbl,
                 bcjf, bckf, bcif, nface);
}

void i2xs(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc,
           int& jf, int& kf, int& if_,
           FortranArray4DRef<double> qjf, FortranArray4DRef<double> qkf,
           FortranArray4DRef<double> qif,
           int& js, int& ks, int& is, int& je, int& ke, int& ie,
           int& nblc, int& ldim, int& nbl,
           FortranArray3DRef<double> bcjf, FortranArray3DRef<double> bckf,
           FortranArray3DRef<double> bcif, int& nface)
{
    int kem, jem, iem;
    int l, i, k, j, ii, kk, jj, kl, jl, il;
    int k2, j2;
    double f1, f2, f4;
    double q[4]; // 1-based: q[1], q[2], q[3]
    int m;

    kem = ke - 1;
    jem = je - 1;
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
                        k2 = std::max(k - 1 + (kl - 1) * 2, 1);
                        k2 = std::min(kc - 1, k2);
                        for (jl = 1; jl <= 3; jl++) {
                            j = js + 1 - jl;
                            q[jl] = f1 * qc(j, k, i, l) + f2 * qc(j, k2, i, l);
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
                        k2 = std::max(k - 1 + (kl - 1) * 2, 1);
                        k2 = std::min(kc - 1, k2);
                        for (jl = 1; jl <= 3; jl++) {
                            j = je + jl - 2;
                            q[jl] = f1 * qc(j, k, i, l) + f2 * qc(j, k2, i, l);
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
    // k = constant planes
    //

    if (nface == 5) {
        //
        // interpolate left boundary
        //
        for (l = 1; l <= ldim; l++) {
            ii = 0;
            for (i = is; i <= iem; i++) {
                ii = ii + 1;
                for (j = js; j <= jem; j++) {
                    for (jl = 1; jl <= 2; jl++) {
                        jj = (j - js) * 2 + jl;
                        j2 = std::max(j - 1 + (jl - 1) * 2, 1);
                        j2 = std::min(jc - 1, j2);
                        for (kl = 1; kl <= 3; kl++) {
                            k = ks + 1 - kl;
                            q[kl] = f1 * qc(j, k, i, l) + f2 * qc(j2, k, i, l);
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
        //
        // interpolate right boundary
        //
        for (l = 1; l <= ldim; l++) {
            ii = 0;
            for (i = is; i <= iem; i++) {
                ii = ii + 1;
                for (j = js; j <= jem; j++) {
                    for (jl = 1; jl <= 2; jl++) {
                        jj = (j - js) * 2 + jl;
                        j2 = std::max(j - 1 + (jl - 1) * 2, 1);
                        j2 = std::min(jc - 1, j2);
                        for (kl = 1; kl <= 3; kl++) {
                            k = ke + kl - 2;
                            q[kl] = f1 * qc(j, k, i, l) + f2 * qc(j2, k, i, l);
                        }
                        qkf(jj, ii, l, 3) = .25 * q[1] + .75 * q[2];
                        qkf(jj, ii, l, 4) = .75 * q[2] + .25 * q[3];
                        bckf(jj, ii, 2) = 0.0;
                    }
                }
            }
        }
    }

    f1 = .75 * .75;
    f2 = .75 * .25;
    f4 = .25 * .25;

    //
    // i = constant plane
    //

    if (nface == 1) {
        //
        // interpolate left boundary
        //
        for (l = 1; l <= ldim; l++) {
            for (k = ks; k <= kem; k++) {
                for (j = js; j <= jem; j++) {
                    for (jl = 1; jl <= 2; jl++) {
                        jj = (j - js) * 2 + jl;
                        j2 = std::max(j - 1 + (jl - 1) * 2, 1);
                        j2 = std::min(jc - 1, j2);
                        for (kl = 1; kl <= 2; kl++) {
                            kk = (k - ks) * 2 + kl;
                            k2 = std::max(k - 1 + (kl - 1) * 2, 1);
                            k2 = std::min(kc - 1, k2);
                            for (il = 1; il <= 2; il++) {
                                i = is - il;
                                q[il] = f1 * qc(j, k, i, l)
                                      + f2 * (qc(j, k2, i, l) + qc(j2, k, i, l))
                                      + f4 * qc(j2, k2, i, l);
                            }
                            qif(jj, kk, l, 1) = q[1];
                            qif(jj, kk, l, 2) = q[2];
                            bcif(jj, kk, 1) = 0.0;
                        }
                    }
                }
            }
        }
    }

    if (nface == 2) {
        //
        // interpolate right boundary
        //
        for (l = 1; l <= ldim; l++) {
            for (k = ks; k <= kem; k++) {
                for (j = js; j <= jem; j++) {
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
                            qif(jj, kk, l, 3) = q[1];
                            qif(jj, kk, l, 4) = q[2];
                            bcif(jj, kk, 2) = 0.0;
                        }
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
}

} // namespace i2xs_ns
