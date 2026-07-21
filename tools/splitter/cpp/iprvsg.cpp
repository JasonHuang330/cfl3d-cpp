// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "iprvsg.h"

namespace iprvsg_ns {

int iprvsg(int& iseg, int& istart, int& mbloc, int& msegt, int& mtot,
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
    int i;
    i = istart;
 label_100:
    if (nxtseg(i) == iseg) goto label_999;
    i = nxtseg(i);
    goto label_100;
 label_999:;
    return i;
}

} // namespace iprvsg_ns
