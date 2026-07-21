// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
//  ---------------------------------------------------------------------------
#include <cstdio>
//  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
//  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
//  with point-matched, patched, overset, or embedded connectivities. Both
//  multigrid and mesh sequencing are available in time-accurate or
//  steady-state modes.
//
//  Copyright 2001 United States Government as represented by the Administrator
//  of the National Aeronautics and Space Administration. All Rights Reserved.
//
//  The CFL3D platform is licensed under the Apache License, Version 2.0
//  (the "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0.
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
//  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
//  License for the specific language governing permissions and limitations
//  under the License.
//  ---------------------------------------------------------------------------

#include "blockj.h"

namespace blockj_ns {

//***********************************************************************
//     Purpose:  Transfer information from block (ir) to qj0 array of
//     block (it).
//***********************************************************************
void blockj(FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
            int& idimr, int& jdimr, int& kdimr, int& idimt, int& kdimt,
            FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
            int& it, int& ir, int& ldim,
            FortranArray3DRef<double> bcj, int& iedge, int& ivolflg)
{
    int ist, iet, iinct;
    int kst, ket, kinct;
    int isr, ier_loc, jsr, jer, ksr, ker;
    int iincr, jincr, kincr;
    int icount, jcount, kcount;
    int ilocr, jlocr, klocr;
    int kloc1r, kloc2r;
    int jloc1r, jloc2r;
    int iloc1r, iloc2r;
    int l, i, k;

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

    isr = limblk(ir, 1);
    ier_loc = limblk(ir, 4);
    jsr = limblk(ir, 2);
    jer = limblk(ir, 5);
    ksr = limblk(ir, 3);
    ker = limblk(ir, 6);

    //
    // determine the side of the q array to transfer from
    //

    //
    // k = constant side
    //
    if (isva(ir, 1) + isva(ir, 2) == 3) {
        if (ksr == 1) {
            kloc1r = 1;
            kloc2r = 2;
        } else {
            kloc1r = kdimr - 1;
            kloc2r = kdimr - 2;
        }

        if (kdimr == 2) {
            kloc1r = 1;
            kloc2r = 1;
        }

        if (jer == jsr) {
            jincr = 1;
        } else {
            jincr = (jer - jsr) / std::abs(jer - jsr);
        }

        if (ier_loc == isr) {
            iincr = 1;
        } else {
            iincr = (ier_loc - isr) / std::abs(ier_loc - isr);
        }

        if ((isva(ir, 1) == isva(it, 1)) || (isva(ir, 2) == isva(it, 2))) {
            //
            // i varies with i     and     j varies with k
            //
            for (l = 1; l <= ldim; l++) {
                icount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    icount = icount + 1;
                    jcount = -1;
                    for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                        jcount = jcount + 1;
                        ilocr = isr + iincr * icount;
                        jlocr = jsr + jincr * jcount;
                        qj0(k, i, l, 1) = q(jlocr, kloc1r, ilocr, l);
                        qj0(k, i, l, 2) = q(jlocr, kloc2r, ilocr, l);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        } else {
            for (l = 1; l <= ldim; l++) {
                jcount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    jcount = jcount + 1;
                    icount = -1;
                    for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                        icount = icount + 1;
                        ilocr = isr + iincr * icount;
                        jlocr = jsr + jincr * jcount;
                        qj0(k, i, l, 1) = q(jlocr, kloc1r, ilocr, l);
                        qj0(k, i, l, 2) = q(jlocr, kloc2r, ilocr, l);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        }

    //
    // j = constant side
    //
    } else if (isva(ir, 1) + isva(ir, 2) == 4) {
        if (jsr == 1) {
            jloc1r = 1;
            jloc2r = 2;
        } else {
            jloc1r = jdimr - 1;
            jloc2r = jdimr - 2;
        }

        if (jdimr == 2) {
            jloc1r = 1;
            jloc2r = 1;
        }

        if (ier_loc == isr) {
            iincr = 1;
        } else {
            iincr = (ier_loc - isr) / std::abs(ier_loc - isr);
        }

        if (ker == ksr) {
            kincr = 1;
        } else {
            kincr = (ker - ksr) / std::abs(ker - ksr);
        }

        if ((isva(ir, 1) == isva(it, 1)) || (isva(ir, 2) == isva(it, 2))) {
            //
            // i varies with i    and    k varies with k
            //
            for (l = 1; l <= ldim; l++) {
                icount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    icount = icount + 1;
                    kcount = -1;
                    for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                        kcount = kcount + 1;
                        ilocr = isr + iincr * icount;
                        klocr = ksr + kincr * kcount;
                        qj0(k, i, l, 1) = q(jloc1r, klocr, ilocr, l);
                        qj0(k, i, l, 2) = q(jloc2r, klocr, ilocr, l);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        } else {
            //
            // k varies with i    and    i varies with k
            //
            for (l = 1; l <= ldim; l++) {
                kcount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    kcount = kcount + 1;
                    icount = -1;
                    for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                        icount = icount + 1;
                        ilocr = isr + iincr * icount;
                        klocr = ksr + kincr * kcount;
                        qj0(k, i, l, 1) = q(jloc1r, klocr, ilocr, l);
                        qj0(k, i, l, 2) = q(jloc2r, klocr, ilocr, l);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        }

    //
    // i = constant side
    //
    } else if (isva(ir, 1) + isva(ir, 2) == 5) {
        if (isr == 1) {
            iloc1r = 1;
            iloc2r = 2;
        } else {
            if (ivolflg == 0) {
                iloc1r = idimr - 1;
                iloc2r = idimr - 2;
            } else {
                iloc1r = idimr;
                iloc2r = idimr - 1;
            }
        }

        if (idimr == 2) {
            iloc1r = 1;
            iloc2r = 1;
        }

        if (jer == jsr) {
            jincr = 1;
        } else {
            jincr = (jer - jsr) / std::abs(jer - jsr);
        }

        if (ker == ksr) {
            kincr = 1;
        } else {
            kincr = (ker - ksr) / std::abs(ker - ksr);
        }

        if ((isva(ir, 1) == isva(it, 1)) || (isva(ir, 2) == isva(it, 2))) {
            //
            // k varies with k    and    j varies with i
            //
            for (l = 1; l <= ldim; l++) {
                jcount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    jcount = jcount + 1;
                    kcount = -1;
                    for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                        kcount = kcount + 1;
                        jlocr = jsr + jincr * jcount;
                        klocr = ksr + kincr * kcount;
                        qj0(k, i, l, 1) = q(jlocr, klocr, iloc1r, l);
                        qj0(k, i, l, 2) = q(jlocr, klocr, iloc2r, l);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        } else {
            //
            // k varies with i    and    j varies with k
            //
            for (l = 1; l <= ldim; l++) {
                kcount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    kcount = kcount + 1;
                    jcount = -1;
                    for (k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                        jcount = jcount + 1;
                        jlocr = jsr + jincr * jcount;
                        klocr = ksr + kincr * kcount;
                        qj0(k, i, l, 1) = q(jlocr, klocr, iloc1r, l);
                        qj0(k, i, l, 2) = q(jlocr, klocr, iloc2r, l);
                        bcj(k, i, iedge) = 0.0;
                    }
                }
            }
        }
    }

    return;
}

} // namespace blockj_ns
