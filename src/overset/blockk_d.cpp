// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "blockk_d.h"
#include "blockk.h"

namespace blockk_d_ns {

void blockk(FortranArray4DRef<double> q, FortranArray4DRef<double> qk0,
            int& idimr, int& jdimr, int& kdimr, int& idimt, int& jdimt,
            FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
            int& it, int& ir, int& ldim, FortranArray3DRef<double> bck,
            int& iedge, int& ivolflg)
{
    blockk_ns::blockk(q, qk0, idimr, jdimr, kdimr, idimt, jdimt,
                      limblk, isva, it, ir, ldim, bck, iedge, ivolflg);
}

void blockk_d(FortranArray3DRef<double> q, FortranArray4DRef<double> qk0,
              int& idimr, int& jdimr, int& kdimr, int& idimt, int& jdimt,
              FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
              int& it, int& ir, int& nvals, int& ldim,
              FortranArray3DRef<double> bck, int& iedge)
{
    int ist, iet, iinct;
    int jst, jet, jinct;
    int ij, i, j, l;

    ist = limblk(it, 1);
    iet = limblk(it, 4);
    if (ist == iet) {
        iinct = 1;
    } else {
        iinct = (iet - ist) / std::abs(iet - ist);
    }

    jst = limblk(it, 2);
    jet = limblk(it, 5);
    if (jst == jet) {
        jinct = 1;
    } else {
        jinct = (jet - jst) / std::abs(jet - jst);
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
            // i varies with i     and     j varies with j
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qk0(j, i, l, 1) = q(ij, l, 1);
                        qk0(j, i, l, 2) = q(ij, l, 2);
                        bck(j, i, iedge) = 0.0;
                    }
                }
            }
        } else {
            //
            // j varies with i     and     i varies with j
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qk0(j, i, l, 1) = q(ij, l, 1);
                        qk0(j, i, l, 2) = q(ij, l, 2);
                        bck(j, i, iedge) = 0.0;
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
            // i varies with i    and    k varies with j
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qk0(j, i, l, 1) = q(ij, l, 1);
                        qk0(j, i, l, 2) = q(ij, l, 2);
                        bck(j, i, iedge) = 0.0;
                    }
                }
            }
        } else {
            //
            // k varies with i    and    i varies with j
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qk0(j, i, l, 1) = q(ij, l, 1);
                        qk0(j, i, l, 2) = q(ij, l, 2);
                        bck(j, i, iedge) = 0.0;
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
            // k varies with i    and    j varies with j
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qk0(j, i, l, 1) = q(ij, l, 1);
                        qk0(j, i, l, 2) = q(ij, l, 2);
                        bck(j, i, iedge) = 0.0;
                    }
                }
            }
        } else {
            //
            // j varies with i    and    k varies with j
            //
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    for (l = 1; l <= ldim; l++) {
                        qk0(j, i, l, 1) = q(ij, l, 1);
                        qk0(j, i, l, 2) = q(ij, l, 2);
                        bck(j, i, iedge) = 0.0;
                    }
                }
            }
        }
    }

    return;
}

} // namespace blockk_d_ns
