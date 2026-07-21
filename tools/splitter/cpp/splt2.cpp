// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "splt2.h"
#include "mkseg.h"

namespace splt2_ns {

void splt2(int& iseg, int& ix1, int& ix2, int& indx, int& is, int& it, int& mbloc, int& msegt, int& mtot,
           FortranArray2DRef<int> imap, FortranArray1DRef<int> idbloc, FortranArray2DRef<int> ivisb,
           FortranArray2DRef<int> itrb, FortranArray1DRef<double> val, FortranArray2DRef<double> xdum,
           FortranArray2DRef<int> iold, FortranArray1DRef<int> nxtseg, FortranArray1DRef<int> intrfc,
           FortranArray1DRef<int> ipatch, FortranArray1DRef<int> nsubbl, FortranArray1DRef<int> idobl,
           FortranArray1DRef<int> nseg, FortranArray1DRef<int> idno, FortranArray2DRef<int> ijk,
           FortranArray1DRef<int> idseg, FortranArray1DRef<int> idnext)
{
    int ns;

    mkseg_ns::mkseg(iseg, it, ns, mbloc, msegt, mtot,
                    imap, idbloc, ivisb, itrb, val, xdum, iold,
                    nxtseg, intrfc, ipatch, nsubbl, idobl,
                    nseg, idno, ijk, idseg, idnext);

    imap(ix2, iseg) = indx;
    imap(ix1, ns) = indx;
    if (imap(13, iseg) != 0) {
        if (imap(ix1 + 11, iseg) >= indx) {
            imap(13, iseg) = 0;
        } else if (imap(ix2 + 11, iseg) > indx) {
            imap(ix2 + 11, iseg) = indx;
        }
        if (imap(ix2 + 11, ns) <= indx) {
            imap(13, ns) = 0;
        } else if (imap(ix1 + 11, ns) < indx) {
            imap(ix1 + 11, ns) = indx;
        }
    }
    iold(ix2 - 2, iseg) = iold(ix1 - 2, iseg) + (indx - imap(ix1, iseg));
    iold(ix1 - 2, ns) = iold(ix2 - 2, iseg);

    return;
}

} // namespace splt2_ns
