// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "mkseg.h"

namespace mkseg_ns {

void mkseg(int& iseg, int& it, int& ns, int& mbloc, int& msegt, int& mtot,
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
    cmn_segment.nsgtop = cmn_segment.nsgtop + 1;
    ns = cmn_segment.nsgtop;
    idbloc(ns) = it;
    nseg(it) = nseg(it) + 1;
    nxtseg(ns) = idseg(it);
    idseg(it) = ns;
    intrfc(ns) = intrfc(iseg);
    for (int i = 1; i <= msegt; i++) {
        imap(i, ns) = imap(i, iseg);
    }
    for (int i = 1; i <= 4; i++) {
        iold(i, ns) = iold(i, iseg);
    }
    for (int i = 1; i <= msegt; i++) {
        xdum(i, ns) = xdum(i, iseg);
    }
    val(ns) = val(iseg);
}

} // namespace mkseg_ns
