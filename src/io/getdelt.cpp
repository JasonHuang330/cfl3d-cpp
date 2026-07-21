// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "getdelt.h"
#include "termn8.h"
#include "runtime/common_blocks.h"
#include <algorithm>
#include <cstdio>

namespace getdelt_ns {

void getdelt(int& maxbl, int& maxsegdg, int& idim, int& jdim, int& kdim,
             FortranArray4DRef<double> delti, FortranArray4DRef<double> deltj,
             FortranArray4DRef<double> deltk,
             FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
             int& nbl, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
             FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
             FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
             int& iseg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim,
             FortranArray3DRef<double> wkj, FortranArray3DRef<double> wkk,
             FortranArray3DRef<double> wki)
{
    // COMMON /mydist2/
    int& myid = cmn_mydist2.myid;

    int is = icsi(nbl, iseg);
    int ie = icsf(nbl, iseg);
    int js = jcsi(nbl, iseg);
    int je = jcsf(nbl, iseg);
    int ks = kcsi(nbl, iseg);
    int ke = kcsf(nbl, iseg);

    if (is == ie) {
        int mm = 1;
        if (is == idim) mm = 2;
        for (int k = ks; k <= ke; k++) {
            for (int j = js; j <= je; j++) {
                delti(j, k, 1, mm) = delti(j, k, 1, mm)
                                   - x(j, k, is) * wki(j, k, mm);
                delti(j, k, 2, mm) = delti(j, k, 2, mm)
                                   - y(j, k, is) * wki(j, k, mm);
                delti(j, k, 3, mm) = delti(j, k, 3, mm)
                                   - z(j, k, is) * wki(j, k, mm);
                wki(j, k, mm) = 0.;
            }
        }
    } else if (js == je) {
        int mm = 1;
        if (js == jdim) mm = 2;
        for (int i = is; i <= ie; i++) {
            for (int k = ks; k <= ke; k++) {
                deltj(k, i, 1, mm) = deltj(k, i, 1, mm)
                                   - x(js, k, i) * wkj(k, i, mm);
                deltj(k, i, 2, mm) = deltj(k, i, 2, mm)
                                   - y(js, k, i) * wkj(k, i, mm);
                deltj(k, i, 3, mm) = deltj(k, i, 3, mm)
                                   - z(js, k, i) * wkj(k, i, mm);
                wkj(k, i, mm) = 0.;
            }
        }
    } else if (ks == ke) {
        int mm = 1;
        if (ks == kdim) mm = 2;
        for (int i = is; i <= ie; i++) {
            for (int j = js; j <= je; j++) {
                deltk(j, i, 1, mm) = deltk(j, i, 1, mm)
                                   - x(j, ks, i) * wkk(j, i, mm);
                deltk(j, i, 2, mm) = deltk(j, i, 2, mm)
                                   - y(j, ks, i) * wkk(j, i, mm);
                deltk(j, i, 3, mm) = deltk(j, i, 3, mm)
                                   - z(j, ks, i) * wkk(j, i, mm);
                wkk(j, i, mm) = 0.;
            }
        }
    } else {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
                      " error in getdelt...one surface dimension must be 1");
        int ierrflg = -1;
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    return;
}

} // namespace getdelt_ns
