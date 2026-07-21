// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "blockj_d.h"
#include "blockj.h"

namespace blockj_d_ns {

void blockj(FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
            int& idimr, int& jdimr, int& kdimr, int& idimt, int& kdimt,
            FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
            int& it, int& ir, int& ldim,
            FortranArray3DRef<double> bcj, int& iedge, int& ivolflg)
{
    blockj_ns::blockj(q, qj0, idimr, jdimr, kdimr, idimt, kdimt,
                      limblk, isva, it, ir, ldim, bcj, iedge, ivolflg);
}

void blockj_d(FortranArray3DRef<double> q, FortranArray4DRef<double> qj0,
              int& idimr, int& jdimr, int& kdimr, int& idimt, int& kdimt,
              FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
              int& it, int& ir, int& nvals, int& ldim,
              FortranArray3DRef<double> bcj, int& iedge)
{
    int ist, iet, iinct;
    int kst, ket, kinct;
    int ij, i, k, l;

    ist = limblk(it, 1);
    iet = limblk(it, 4);
    if (ist == iet) {
        iinct = 1;
    } else {
        iinct = (iet - ist) / std::abs(iet - ist);
    }

    kst = limblk(it, 3);
    ket = limblk(it, 6);
    if (kst == ket) {
        kinct = 1;
    } else {
        kinct = (ket - kst) / std::abs(ket - kst);
    }

    //
    // determine the side of the q array to transfer from
    //
    // k = constant side
    //
    if (isva(ir, 1) + isva(ir, 2) == 3) {
        if ((isva(ir, 1) == isva(it, 1)) ||
            (isva(ir, 2) == isva(it, 2))) {
            //
            // i varies with i     and     j varies with k
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qj0(k, i, l, 1) = q(ij, l, 1);
                        qj0(k, i, l, 2) = q(ij, l, 2);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        } else {
            //
            // j varies with i     and     i varies with j
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qj0(k, i, l, 1) = q(ij, l, 1);
                        qj0(k, i, l, 2) = q(ij, l, 2);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        }
    //
    // j = constant side
    //
    } else if (isva(ir, 1) + isva(ir, 2) == 4) {
        if ((isva(ir, 1) == isva(it, 1)) ||
            (isva(ir, 2) == isva(it, 2))) {
            //
            // i varies with i    and    k varies with k
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qj0(k, i, l, 1) = q(ij, l, 1);
                        qj0(k, i, l, 2) = q(ij, l, 2);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        } else {
            //
            // k varies with i    and    i varies with k
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qj0(k, i, l, 1) = q(ij, l, 1);
                        qj0(k, i, l, 2) = q(ij, l, 2);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        }
    //
    // i = constant side
    //
    } else if (isva(ir, 1) + isva(ir, 2) == 5) {
        if ((isva(ir, 1) == isva(it, 1)) ||
            (isva(ir, 2) == isva(it, 2))) {
            //
            // k varies with k    and    j varies with i
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qj0(k, i, l, 1) = q(ij, l, 1);
                        qj0(k, i, l, 2) = q(ij, l, 2);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        } else {
            //
            // k varies with i    and    j varies with k
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qj0(k, i, l, 1) = q(ij, l, 1);
                        qj0(k, i, l, 2) = q(ij, l, 2);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        }
    }
}

} // namespace blockj_d_ns
