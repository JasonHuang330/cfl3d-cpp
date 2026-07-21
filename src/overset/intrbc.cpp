// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "intrbc.h"
#include <cstdio>
#include "augmntq.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstring>
#include <algorithm>

namespace intrbc_ns {

void intrbc(FortranArray4DRef<double> q, int& jdim, int& kdim, int& idim, int& nbl, int& ldim, int& maxbl, int& iitot, FortranArray1DRef<int> lig, FortranArray1DRef<int> iipntsg, FortranArray1DRef<double> dxintg, FortranArray1DRef<double> dyintg, FortranArray1DRef<double> dzintg, FortranArray1DRef<int> iiig, FortranArray1DRef<int> jjig, FortranArray1DRef<int> kkig, FortranArray3DRef<double> qb, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qq, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& icorr, int& iset)
{
    // COMMON blocks
    int& isklton = cmn_sklton.isklton;
    int& i2d     = cmn_twod.i2d;
    int& iexp    = cmn_zero.iexp;

    // 10.**(-iexp) is machine zero
    double qlimit = std::pow(10.0, (double)(-iexp));

    // Note: iset is index for loading into appropriate qb array
    //       (1 for q, 2 for vk0,vj0,vi0, 3 for tj0,tk0,ti0)

    int lsta = lig(nbl);
    int lend = lsta + iipntsg(nbl) - 1;
    if (iipntsg(nbl) == 0) return;

    // qq is "augmented" q array of dimensions jdim+1 x kdim+1 x idim+1
    // containing cell-center data at interior points and cell-face
    // center data at boundary points

    augmntq_ns::augmntq(q, jdim, kdim, idim, nbl, ldim, qj0, qk0, qi0, qq,
                        bcj, bck, bci, nou, bou, nbuf, ibufdim, icorr);

    if (i2d == 0) {
        for (int n = 1; n <= ldim; n++) {
            for (int l = lsta; l <= lend; l++) {
                // set up interpolation coefficients
                double s1 = qq(jjig(l)+1, kkig(l)+1, iiig(l)+1, n);
                double s2 = qq(jjig(l)+1, kkig(l)+1, iiig(l)+2, n);
                double s3 = qq(jjig(l)+2, kkig(l)+1, iiig(l)+2, n);
                double s4 = qq(jjig(l)+2, kkig(l)+1, iiig(l)+1, n);
                double s5 = qq(jjig(l)+1, kkig(l)+2, iiig(l)+1, n);
                double s6 = qq(jjig(l)+1, kkig(l)+2, iiig(l)+2, n);
                double s7 = qq(jjig(l)+2, kkig(l)+2, iiig(l)+2, n);
                double s8 = qq(jjig(l)+2, kkig(l)+2, iiig(l)+1, n);

                double a1 =  s1;
                double a2 = -s1+s2;
                double a3 = -s1+s4;
                double a4 = -s1+s5;
                double a5 =  s1-s2+s3-s4;
                double a6 =  s1-s2-s5+s6;
                double a7 =  s1-s4-s5+s8;
                double a8 = -s1+s2-s3+s4+s5-s6+s7-s8;

                // interpolate and store in qb array
                qb(l, n, iset) = a1 + a2*dxintg(l)
                               + a3*dyintg(l)
                               + a4*dzintg(l)
                               + a5*dxintg(l)*dyintg(l)
                               + a6*dxintg(l)*dzintg(l)
                               + a7*dyintg(l)*dzintg(l)
                               + a8*dxintg(l)*dyintg(l)*dzintg(l);
            }
        }
    } else {
        for (int n = 1; n <= ldim; n++) {
            for (int l = lsta; l <= lend; l++) {
                // set up interpolation coefficients
                double s1 = qq(jjig(l)+1, kkig(l)+1, 2, n);
                double s4 = qq(jjig(l)+2, kkig(l)+1, 2, n);
                double s5 = qq(jjig(l)+1, kkig(l)+2, 2, n);
                double s8 = qq(jjig(l)+2, kkig(l)+2, 2, n);

                double a3 = -s1+s4;
                double a4 = -s1+s5;
                double a7 =  s1-s4-s5+s8;

                // interpolate and store in qb array
                qb(l, n, iset) = s1 + a3*dyintg(l) + a4*dzintg(l) +
                                 a7*dyintg(l)*dzintg(l);
            }
        }
    }

    // limit negative rho,p,vist3d,turb values to a small positive value
    if (iset == 1) {
        for (int n = 1; n <= 5; n += 4) {
            for (int l = lsta; l <= lend; l++) {
                if ((float)qb(l, n, iset) < 0.f) qb(l, n, iset) = qlimit;
            }
        }
    } else {
        for (int n = 1; n <= ldim; n++) {
            for (int l = lsta; l <= lend; l++) {
                if ((float)qb(l, n, iset) < 0.f) {
                    qb(l, n, iset) = qlimit;
                }
            }
        }
    }



    if (isklton == 1 && iset == 1) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        char tmp[120];
        std::snprintf(tmp, sizeof(tmp), "  qb points updated = %d", iipntsg(nbl));
        std::memset(bou(nou(1), 1), ' ', 120);
        int len = (int)std::strlen(tmp);
        std::memcpy(bou(nou(1), 1), tmp, len < 120 ? len : 120);
    }
}

} // namespace intrbc_ns
