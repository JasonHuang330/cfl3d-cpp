// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "ipadd.h"

namespace ipadd_ns {

void ipadd(int& ibl, int& iifc, int& mbloc, int& msegt, int& mtot,
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
    int ifc, itop, i, j;

    ifc = std::abs(iifc);
    itop = ipatch(1);
    i = 2;

    // Label 1: loop checking existing patches
label1:
    if (i >= itop) goto label2;
    j = i + 1;
    if (ipatch(i) == ibl && ipatch(j) == ifc) return;
    i = j + 1;
    goto label1;

    // Label 2: add new patch entry
label2:
    ipatch(itop) = ibl;
    itop = itop + 1;
    ipatch(itop) = ifc;
    itop = itop + 1;
    ipatch(1) = itop;
    //      write(6,*) 'adding block ',ibl,', face ',ifc

    return;
}

} // namespace ipadd_ns
