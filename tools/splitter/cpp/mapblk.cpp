// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "mapblk.h"

namespace mapblk_ns {

void mapblk(int& nbl, FortranArray1DRef<int> id, FortranArray1DRef<int> imin, FortranArray1DRef<int> imax, FortranArray1DRef<int> jmin, FortranArray1DRef<int> jmax, FortranArray1DRef<int> kmin, FortranArray1DRef<int> kmax, int& mbloc, int& msegt, int& mtot, FortranArray2DRef<int> imap, FortranArray1DRef<int> idbloc, FortranArray2DRef<int> ivisb, FortranArray2DRef<int> itrb, FortranArray1DRef<double> val, FortranArray2DRef<double> xdum, FortranArray2DRef<int> iold, FortranArray1DRef<int> nxtseg, FortranArray1DRef<int> intrfc, FortranArray1DRef<int> ipatch, FortranArray1DRef<int> nsubbl, FortranArray1DRef<int> idobl, FortranArray1DRef<int> nseg, FortranArray1DRef<int> idno, FortranArray2DRef<int> ijk, FortranArray1DRef<int> idseg, FortranArray1DRef<int> idnext)
{
    // common /oldbl/ noldbl
    int& noldbl = cmn_oldbl.noldbl;

    nbl = 0;
    for (int n = 1; n <= noldbl; n++) {
        int ibl = idobl(n);
        for (int n2 = 1; n2 <= nsubbl(n); n2++) {
            nbl = nbl + 1;
            id(nbl) = n;
            imin(nbl) = ijk(1, ibl);
            imax(nbl) = ijk(2, ibl);
            jmin(nbl) = ijk(3, ibl);
            jmax(nbl) = ijk(4, ibl);
            kmin(nbl) = ijk(5, ibl);
            kmax(nbl) = ijk(6, ibl);
            ibl = idnext(ibl);
        }
    }
}

} // namespace mapblk_ns
