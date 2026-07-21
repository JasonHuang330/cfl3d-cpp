// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "blocki_d.h"
#include "blocki.h"

namespace blocki_d_ns {

void blocki(FortranArray4DRef<double> q, FortranArray4DRef<double> qi0,
            int& idimr, int& jdimr, int& kdimr, int& jdimt, int& kdimt,
            FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
            int& it, int& ir, int& ldim, FortranArray3DRef<double> bci,
            int& iedge, int& ivolflg)
{
    blocki_ns::blocki(q, qi0, idimr, jdimr, kdimr, jdimt, kdimt,
                      limblk, isva, it, ir, ldim, bci, iedge, ivolflg);
}

void blocki_d(FortranArray3DRef<double> q, FortranArray4DRef<double> qi0,
              int& idimr, int& jdimr, int& kdimr, int& jdimt, int& kdimt,
              FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
              int& it, int& ir, int& nvals, int& ldim,
              FortranArray3DRef<double> bci, int& iedge)
{
    int jst, jet, jinct;
    int kst, ket, kinct;
    int ij, j, k, l;

    jst = limblk(it, 2);
    jet = limblk(it, 5);
    if (jst == jet) {
        jinct = 1;
    } else {
        jinct = (jet - jst) / std::abs(jet - jst);
    }

    kst = limblk(it, 3);
    ket = limblk(it, 6);
    if (kst == ket) {
        kinct = 1;
    } else {
        kinct = (ket - kst) / std::abs(ket - kst);
    }

    // determine the side of the q array to transfer from

    // k = constant side
    if (isva(ir, 1) + isva(ir, 2) == 3) {
        if ((isva(ir, 1) == isva(it, 1)) || (isva(ir, 2) == isva(it, 2))) {
            // i varies with k     and     j varies with j
            ij = 0;
            for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qi0(j, k, l, 1) = q(ij, l, 1);
                        qi0(j, k, l, 2) = q(ij, l, 2);
                        bci(j, k, iedge) = 0.0;
                    }
                }
            }
        } else {
            // j varies with k     and     i varies with j
            ij = 0;
            for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qi0(j, k, l, 1) = q(ij, l, 1);
                        qi0(j, k, l, 2) = q(ij, l, 2);
                        bci(j, k, iedge) = 0.0;
                    }
                }
            }
        }

    // j = constant side
    } else if (isva(ir, 1) + isva(ir, 2) == 4) {
        if ((isva(ir, 1) == isva(it, 1)) || (isva(ir, 2) == isva(it, 2))) {
            // i varies with j    and    k varies with k
            ij = 0;
            for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qi0(j, k, l, 1) = q(ij, l, 1);
                        qi0(j, k, l, 2) = q(ij, l, 2);
                        bci(j, k, iedge) = 0.0;
                    }
                }
            }
        } else {
            // i varies with k    and    k varies with j
            ij = 0;
            for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qi0(j, k, l, 1) = q(ij, l, 1);
                        qi0(j, k, l, 2) = q(ij, l, 2);
                        bci(j, k, iedge) = 0.0;
                    }
                }
            }
        }

    // i = constant side
    } else if (isva(ir, 1) + isva(ir, 2) == 5) {
        if ((isva(ir, 1) == isva(it, 1)) || (isva(ir, 2) == isva(it, 2))) {
            // k varies with k    and    j varies with j
            ij = 0;
            for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qi0(j, k, l, 1) = q(ij, l, 1);
                        qi0(j, k, l, 2) = q(ij, l, 2);
                        bci(j, k, iedge) = 0.0;
                    }
                }
            }
        } else {
            // j varies with k    and    k varies with j
            ij = 0;
            for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qi0(j, k, l, 1) = q(ij, l, 1);
                        qi0(j, k, l, 2) = q(ij, l, 2);
                        bci(j, k, iedge) = 0.0;
                    }
                }
            }
        }
    }

    return;
}

} // namespace blocki_d_ns
