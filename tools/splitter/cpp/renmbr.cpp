// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "renmbr.h"

namespace renmbr_ns {

void renmbr(int& iseg, int& ix1, int& ix2, int& indx, int& mbloc, int& msegt, int& mtot,
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
    imap(ix1, iseg) = imap(ix1, iseg) - (indx - 1);
    imap(ix2, iseg) = imap(ix2, iseg) - (indx - 1);
    if (imap(13, iseg) != 0) {
        imap(ix1 + 11, iseg) = imap(ix1 + 11, iseg) - (indx - 1);
        imap(ix2 + 11, iseg) = imap(ix2 + 11, iseg) - (indx - 1);
    }
    int jseg = intrfc(iseg);
    if (jseg > 0) {
        int jx1 = ix1 + 6;
        int jx2 = ix2 + 6;
        if (imap(jx1, jseg) < imap(jx2, jseg)) {
            imap(jx1, jseg) = imap(ix1, iseg);
            imap(jx2, jseg) = imap(ix2, iseg);
        } else {
            imap(jx1, jseg) = imap(ix2, iseg);
            imap(jx2, jseg) = imap(ix1, iseg);
        }
    }

    return;
}

} // namespace renmbr_ns
