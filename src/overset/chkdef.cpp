// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "chkdef.h"
#include "termn8.h"
#include <algorithm>
#include <cstdio>

namespace chkdef_ns {

void chkdef(int& nbl, int& idim, int& jdim, int& kdim,
            FortranArray3DRef<double> bci,
            FortranArray3DRef<double> bcj,
            FortranArray3DRef<double> bck,
            FortranArray2DRef<int> icsi,
            FortranArray2DRef<int> icsf,
            FortranArray2DRef<int> jcsi,
            FortranArray2DRef<int> jcsf,
            FortranArray2DRef<int> kcsi,
            FortranArray2DRef<int> kcsf,
            FortranArray1DRef<int> nsegdfrm,
            FortranArray2DRef<int> idfrmseg,
            int& maxbl, int& maxsegdg,
            FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid)
{
    int iseg, j, k, i, ll;

    for (iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {

        // check i=const surfaces
        if (icsi(nbl,iseg) == icsf(nbl,iseg)) {
            if (icsi(nbl,iseg) == 1) {
                ll = 1;
            } else {
                ll = 2;
            }
            for (j = jcsi(nbl,iseg); j <= jcsf(nbl,iseg)-1; j++) {
                for (k = kcsi(nbl,iseg); k <= kcsf(nbl,iseg)-1; k++) {
                    if (bci(j,k,ll) == 0.) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "stopping...deforming face has interface type bc"
                            " block%4d segment%3d", nbl, iseg);
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "  at j,k%4d%4d", j, k);
                        int m1 = -1;
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }

        // check j=const surfaces
        if (jcsi(nbl,iseg) == jcsf(nbl,iseg)) {
            if (jcsi(nbl,iseg) == 1) {
                ll = 1;
            } else {
                ll = 2;
            }
            for (i = icsi(nbl,iseg); i <= icsf(nbl,iseg)-1; i++) {
                for (k = kcsi(nbl,iseg); k <= kcsf(nbl,iseg)-1; k++) {
                    if (bcj(k,i,ll) == 0.) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "stopping...deforming face has interface type bc"
                            " block%4d segment%3d", nbl, iseg);
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "  at k,i%4d%4d", k, i);
                        int m1 = -1;
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }

        // check k=const surfaces
        if (kcsi(nbl,iseg) == kcsf(nbl,iseg)) {
            if (kcsi(nbl,iseg) == 1) {
                ll = 1;
            } else {
                ll = 2;
            }
            for (j = jcsi(nbl,iseg); j <= jcsf(nbl,iseg)-1; j++) {
                for (i = icsi(nbl,iseg); i <= icsf(nbl,iseg)-1; i++) {
                    if (bck(j,i,ll) == 0.) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "stopping...deforming face has interface type bc"
                            " block%4d segment%3d", nbl, iseg);
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "  at j,i%4d%4d", j, i);
                        int m1 = -1;
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }

    } // end do iseg

    return;
}

} // namespace chkdef_ns
