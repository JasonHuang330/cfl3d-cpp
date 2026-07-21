// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include <cstdio>
#include "pre_patch.h"

namespace pre_patch_ns {

void pre_patch(int& nbl, FortranArray2DRef<int> lw, int& icount, int& ninter,
               FortranArray2DRef<int> iindex, int& intmax, int& nsub1,
               FortranArray2DRef<int> isav_pat, FortranArray3DRef<int> isav_pat_b,
               FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1,
               FortranArray1DRef<int> iiint1, FortranArray1DRef<int> iiint2,
               int& maxbl, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg,
               FortranArray1DRef<int> idimg, int& ierrflg)
{
    // local variables
    int icheck, lmax1, nblcc, lst, npt;
    int j21, j22, k21, k22;
    int lqedge, lqedgt, lqedgv, lqedgb;
    int jmax2, kmax2, mint1, mint2;
    int l, mbl, mtype;

    if (std::abs(ninter) > 0) {

        // loop over all patch interfaces
        for (icheck = 1; icheck <= std::abs(ninter); icheck++) {
            lmax1 = iindex(icheck, 1);
            nblcc = iindex(icheck, lmax1 + 2);
            if (nblcc != nbl) continue;
            lst = iindex(icheck, 2*lmax1 + 5);
            npt = iindex(icheck, 2*lmax1 + 4);

            // set range of points requiring interpolation on "to" side of
            // patch interface
            j21 = iindex(icheck, 2*lmax1 + 6);
            j22 = iindex(icheck, 2*lmax1 + 7);
            k21 = iindex(icheck, 2*lmax1 + 8);
            k22 = iindex(icheck, 2*lmax1 + 9);

            // patch surface to be interpolated to is an i=constant surface
            if (iindex(icheck, 2*lmax1 + 3) / 10 == 1) {
                lqedge = lw(4,  nbl);
                lqedgt = lw(25, nbl);
                lqedgv = lw(30, nbl);
                lqedgb = lw(33, nbl);
                jmax2  = jdimg(nbl);
                kmax2  = kdimg(nbl);
                if (iindex(icheck, 2*lmax1 + 3) == 11) {
                    mint1 = 1;
                    mint2 = 2;
                } else {
                    mint1 = 3;
                    mint2 = 4;
                }
            }

            // patch surface to be interpolated to is a j=constant surface
            if (iindex(icheck, 2*lmax1 + 3) / 10 == 2) {
                lqedge = lw(2,  nbl);
                lqedgt = lw(23, nbl);
                lqedgv = lw(28, nbl);
                lqedgb = lw(31, nbl);
                jmax2  = kdimg(nbl);
                kmax2  = idimg(nbl) - 1;
                if (iindex(icheck, 2*lmax1 + 3) == 21) {
                    mint1 = 1;
                    mint2 = 2;
                } else {
                    mint1 = 3;
                    mint2 = 4;
                }
            }

            // patch surface to be interpolated to is a k=constant surface
            if (iindex(icheck, 2*lmax1 + 3) / 10 == 3) {
                lqedge = lw(3,  nbl);
                lqedgt = lw(24, nbl);
                lqedgv = lw(29, nbl);
                lqedgb = lw(32, nbl);
                jmax2  = jdimg(nbl);
                kmax2  = idimg(nbl) - 1;
                if (iindex(icheck, 2*lmax1 + 3) == 31) {
                    mint1 = 1;
                    mint2 = 2;
                } else {
                    mint1 = 3;
                    mint2 = 4;
                }
            }

            // put some relevant values into isav_pat array
            icount = icount + 1;
            isav_pat(icount,  1) = nbl;
            isav_pat(icount,  2) = lmax1;
            isav_pat(icount,  3) = j21;
            isav_pat(icount,  4) = j22;
            isav_pat(icount,  5) = k21;
            isav_pat(icount,  6) = k22;
            isav_pat(icount,  7) = lqedge;
            isav_pat(icount,  8) = lqedgt;
            isav_pat(icount,  9) = lqedgv;
            isav_pat(icount, 10) = lqedgb;
            isav_pat(icount, 11) = jmax2;
            isav_pat(icount, 12) = kmax2;
            isav_pat(icount, 13) = mint1;
            isav_pat(icount, 14) = mint2;
            isav_pat(icount, 15) = lst;
            isav_pat(icount, 16) = npt;
            isav_pat(icount, 17) = icheck;

            for (l = 1; l <= lmax1; l++) {
                mbl   = iindex(icheck, l + 1);
                mtype = iindex(icheck, l + lmax1 + 2);

                // patch surface to be interpolated from is an i=constant surface
                if (mtype / 10 == 1) {
                    jjmax1(l) = jdimg(mbl);
                    kkmax1(l) = kdimg(mbl);
                    if (mtype == 11) {
                        iiint1(l) = 1;
                        iiint2(l) = std::min(2, idimg(mbl) - 1);
                    } else {
                        iiint1(l) = idimg(mbl) - 1;
                        iiint2(l) = std::max(1, idimg(mbl) - 2);
                    }
                }

                // patch surface to be interpolated from is a j=constant surface
                if (mtype / 10 == 2) {
                    jjmax1(l) = kdimg(mbl);
                    kkmax1(l) = idimg(mbl);
                    if (mtype == 21) {
                        iiint1(l) = 1;
                        iiint2(l) = std::min(2, jdimg(mbl) - 1);
                    } else {
                        iiint1(l) = jdimg(mbl) - 1;
                        iiint2(l) = std::max(1, jdimg(mbl) - 2);
                    }
                }

                // patch surface to be interpolated from is a k=constant surface
                if (mtype / 10 == 3) {
                    jjmax1(l) = jdimg(mbl);
                    kkmax1(l) = idimg(mbl);
                    if (mtype == 31) {
                        iiint1(l) = 1;
                        iiint2(l) = std::min(2, kdimg(mbl) - 1);
                    } else {
                        iiint1(l) = kdimg(mbl) - 1;
                        iiint2(l) = std::max(1, kdimg(mbl) - 2);
                    }
                }

                isav_pat_b(icount, l, 1) = mbl;
                isav_pat_b(icount, l, 2) = mtype / 10;
                isav_pat_b(icount, l, 3) = jjmax1(l);
                isav_pat_b(icount, l, 4) = kkmax1(l);
                isav_pat_b(icount, l, 5) = iiint1(l);
                isav_pat_b(icount, l, 6) = iiint2(l);

            } // end do 1705

        } // end do 8001

    } // end if abs(ninter) > 0

    return;
}

} // namespace pre_patch_ns
