// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "cntsurf.h"
#include <cmath>

namespace cntsurf_ns {

void cntsurf(int& nsurf, int& maxbl, int& maxgr, int& maxseg, int& ngrid,
             FortranArray1DRef<int> nblg,
             FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0,
             FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim,
             FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim,
             FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo,
             FortranArray4DRef<int> kbcinfo, int& ibctyp)
{
    int igrid, nbl, m, ns, nface, iseg;
    int nbctype, n1beg, n1end, n2beg, n2end;

    nsurf = 0;
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nblg(igrid);

        // i=constant surfaces
        for (m = 1; m <= 2; m++) {
            if (m == 1) {
                ns    = nbci0(nbl);
                nface = 1;
            } else {
                ns    = nbcidim(nbl);
                nface = 2;
            }
            for (iseg = 1; iseg <= ns; iseg++) {
                nbctype = ibcinfo(nbl, iseg, 1, m);
                n1beg   = ibcinfo(nbl, iseg, 2, m);
                n1end   = ibcinfo(nbl, iseg, 3, m);
                n2beg   = ibcinfo(nbl, iseg, 4, m);
                n2end   = ibcinfo(nbl, iseg, 5, m);
                if (std::abs(nbctype) == ibctyp) {
                    nsurf = nsurf + (n2end - n2beg + 1) * (n1end - n1beg + 1);
                }
            }
        }

        // j=constant surfaces
        for (m = 1; m <= 2; m++) {
            if (m == 1) {
                ns    = nbcj0(nbl);
                nface = 3;
            } else {
                ns    = nbcjdim(nbl);
                nface = 4;
            }
            for (iseg = 1; iseg <= ns; iseg++) {
                nbctype = jbcinfo(nbl, iseg, 1, m);
                n1beg   = jbcinfo(nbl, iseg, 4, m);
                n1end   = jbcinfo(nbl, iseg, 5, m);
                n2beg   = jbcinfo(nbl, iseg, 2, m);
                n2end   = jbcinfo(nbl, iseg, 3, m);
                if (std::abs(nbctype) == ibctyp) {
                    nsurf = nsurf + (n2end - n2beg + 1) * (n1end - n1beg + 1);
                }
            }
        }

        // k=constant surfaces
        for (m = 1; m <= 2; m++) {
            if (m == 1) {
                ns    = nbck0(nbl);
                nface = 5;
            } else {
                ns    = nbckdim(nbl);
                nface = 6;
            }
            for (iseg = 1; iseg <= ns; iseg++) {
                nbctype = kbcinfo(nbl, iseg, 1, m);
                n1beg   = kbcinfo(nbl, iseg, 2, m);
                n1end   = kbcinfo(nbl, iseg, 3, m);
                n2beg   = kbcinfo(nbl, iseg, 4, m);
                n2end   = kbcinfo(nbl, iseg, 5, m);
                if (std::abs(nbctype) == ibctyp) {
                    nsurf = nsurf + (n2end - n2beg + 1) * (n1end - n1beg + 1);
                }
            }
        }

    } // end igrid loop
}

} // namespace cntsurf_ns
