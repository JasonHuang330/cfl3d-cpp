// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "movseg.h"
#include "iprvsg.h"

namespace movseg_ns {

void movseg(int& iseg, int& is, int& it, int& mbloc, int& msegt, int& mtot,
            FortranArray2DRef<int> imap, FortranArray1DRef<int> idbloc,
            FortranArray2DRef<int> ivisb, FortranArray2DRef<int> itrb,
            FortranArray1DRef<double> val, FortranArray2DRef<double> xdum,
            FortranArray2DRef<int> iold, FortranArray1DRef<int> nxtseg,
            FortranArray1DRef<int> intrfc, FortranArray1DRef<int> ipatch,
            FortranArray1DRef<int> nsubbl, FortranArray1DRef<int> idobl,
            FortranArray1DRef<int> nseg, FortranArray1DRef<int> idno,
            FortranArray2DRef<int> ijk, FortranArray1DRef<int> idseg,
            FortranArray1DRef<int> idnext)
{
    nseg(is) = nseg(is) - 1;
    nseg(it) = nseg(it) + 1;
    if (iseg == idseg(is)) {
        idseg(is) = nxtseg(iseg);
    } else {
        int itmp = iprvsg_ns::iprvsg(iseg, idseg(is), mbloc, msegt, mtot,
                                     imap, idbloc, ivisb, itrb, val, xdum,
                                     iold, nxtseg, intrfc, ipatch, nsubbl,
                                     idobl, nseg, idno, ijk, idseg, idnext);
        nxtseg(itmp) = nxtseg(iseg);
    }
    nxtseg(iseg) = idseg(it);
    idseg(it) = iseg;
    idbloc(iseg) = it;
    int itmp = intrfc(iseg);
    if (itmp != 0) {
        imap(7, itmp) = it;
    }

    return;
}

} // namespace movseg_ns
