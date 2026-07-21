// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "blockk.h"

namespace blockk_ns {

void blockk(FortranArray4DRef<double> q, FortranArray4DRef<double> qk0,
            int& idimr, int& jdimr, int& kdimr, int& idimt, int& jdimt,
            FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
            int& it, int& ir, int& ldim,
            FortranArray3DRef<double> bck, int& iedge, int& ivolflg)
{
//***********************************************************************
//     Purpose:  Transfer information from block (ir) to qk0 array of
//     block (it).
//***********************************************************************

    int ist, iet, iinct;
    int jst, jet, jinct;
    int isr, ier_val, jsr, jer, ksr, ker;
    int iincr, jincr, kincr;
    int kloc1r, kloc2r;
    int jloc1r, jloc2r;
    int iloc1r, iloc2r;
    int icount, jcount, kcount;
    int ilocr, jlocr, klocr;
    int l, i, j;

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

    isr = limblk(ir, 1);
    ier_val = limblk(ir, 4);
    jsr = limblk(ir, 2);
    jer = limblk(ir, 5);
    ksr = limblk(ir, 3);
    ker = limblk(ir, 6);

//
//     determine the side of the q array to transfer from
//
//     k = constant side
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

        if (ier_val == isr) {
            iincr = 1;
        } else {
            iincr = (ier_val - isr) / std::abs(ier_val - isr);
        }

        if ((isva(ir, 1) == isva(it, 1)) || (isva(ir, 2) == isva(it, 2))) {
//
//     i varies with i     and     j varies with j
//
            for (l = 1; l <= ldim; l++) {
                icount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    icount = icount + 1;
                    jcount = -1;
                    for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                        jcount = jcount + 1;
                        ilocr = isr + iincr * icount;
                        jlocr = jsr + jincr * jcount;
                        qk0(j, i, l, 1) = q(jlocr, kloc1r, ilocr, l);
                        qk0(j, i, l, 2) = q(jlocr, kloc2r, ilocr, l);
                        bck(j, i, iedge) = 0.0;
                    }
                }
            }
        } else {
//
//     j varies with i     and     i varies with j
//
            for (l = 1; l <= ldim; l++) {
                jcount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    jcount = jcount + 1;
                    icount = -1;
                    for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                        icount = icount + 1;
                        ilocr = isr + iincr * icount;
                        jlocr = jsr + jincr * jcount;
                        qk0(j, i, l, 1) = q(jlocr, kloc1r, ilocr, l);
                        qk0(j, i, l, 2) = q(jlocr, kloc2r, ilocr, l);
                        bck(j, i, iedge) = 0.0;
                    }
                }
            }
        }

//
//     j = constant side
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

        if (ier_val == isr) {
            iincr = 1;
        } else {
            iincr = (ier_val - isr) / std::abs(ier_val - isr);
        }

        if (ker == ksr) {
            kincr = 1;
        } else {
            kincr = (ker - ksr) / std::abs(ker - ksr);
        }

        if ((isva(ir, 1) == isva(it, 1)) || (isva(ir, 2) == isva(it, 2))) {
//
//     i varies with i    and    k varies with j
//
            for (l = 1; l <= ldim; l++) {
                icount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    icount = icount + 1;
                    kcount = -1;
                    for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                        kcount = kcount + 1;
                        ilocr = isr + iincr * icount;
                        klocr = ksr + kincr * kcount;
                        qk0(j, i, l, 1) = q(jloc1r, klocr, ilocr, l);
                        qk0(j, i, l, 2) = q(jloc2r, klocr, ilocr, l);
                        bck(j, i, iedge) = 0.0;
                    }
                }
            }
        } else {
//
//     k varies with i    and    i varies with j
//
            for (l = 1; l <= ldim; l++) {
                kcount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    kcount = kcount + 1;
                    icount = -1;
                    for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                        icount = icount + 1;
                        ilocr = isr + iincr * icount;
                        klocr = ksr + kincr * kcount;
                        qk0(j, i, l, 1) = q(jloc1r, klocr, ilocr, l);
                        qk0(j, i, l, 2) = q(jloc2r, klocr, ilocr, l);
                        bck(j, i, iedge) = 0.0;
                    }
                }
            }
        }

//
//     i = constant side
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
//     k varies with i    and    j varies with j
//
            for (l = 1; l <= ldim; l++) {
                kcount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    kcount = kcount + 1;
                    jcount = -1;
                    for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                        jcount = jcount + 1;
                        jlocr = jsr + jincr * jcount;
                        klocr = ksr + kincr * kcount;
                        qk0(j, i, l, 1) = q(jlocr, klocr, iloc1r, l);
                        qk0(j, i, l, 2) = q(jlocr, klocr, iloc2r, l);
                        bck(j, i, iedge) = 0.0;
                    }
                }
            }
        } else {
//
//     j varies with i    and    k varies with j
//
            for (l = 1; l <= ldim; l++) {
                jcount = -1;
                for (i = ist; (iinct > 0) ? (i <= iet) : (i >= iet); i += iinct) {
                    jcount = jcount + 1;
                    kcount = -1;
                    for (j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                        kcount = kcount + 1;
                        jlocr = jsr + jincr * jcount;
                        klocr = ksr + kincr * kcount;
                        qk0(j, i, l, 1) = q(jlocr, klocr, iloc1r, l);
                        qk0(j, i, l, 2) = q(jlocr, klocr, iloc2r, l);
                        bck(j, i, iedge) = 0.0;
                    }
                }
            }
        }
    }

    return;
}

} // namespace blockk_ns
