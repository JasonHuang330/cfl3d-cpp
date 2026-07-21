// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "mkintr.h"

namespace mkintr_ns {

void mkintr(int& imn, int& i1, int& i2, int& mbloc, int& msegt, int& mtot,
            FortranArray2DRef<int> imap, FortranArray1DRef<int> idbloc,
            FortranArray2DRef<int> ivisb, FortranArray2DRef<int> itrb,
            FortranArray1DRef<double> val, FortranArray2DRef<double> xdum,
            FortranArray2DRef<int> iold, FortranArray1DRef<int> nxtseg,
            FortranArray1DRef<int> intrfc, FortranArray1DRef<int> ipatch,
            FortranArray1DRef<int> nsubbl, FortranArray1DRef<int> idobl,
            FortranArray1DRef<int> nseg, FortranArray1DRef<int> idno,
            FortranArray2DRef<int> ijk, FortranArray1DRef<int> idseg,
            FortranArray1DRef<int> idnext, FortranArray1DRef<double> xdmold,
            FortranArray1DRef<int> idmold)
{
    // Local variables
    int imx = 0;
    int jmx = 0;
    int iomn = 1000000;
    int iomx = 0;
    int jomn = 1000000;
    int jomx = 0;

    int iseg = idseg(i1);

    // Loop 100: while iseg != 0
    while (iseg != 0) {
        if (imap(2, iseg) == imn) {
            // do 101 i=1,msegt: xdmold(i) = xdum(i,iseg)
            for (int i = 1; i <= msegt; i++) {
                xdmold(i) = xdum(i, iseg);
            }
            // do i=21,msegt: idmold(i) = imap(i,iseg)
            for (int i = 21; i <= msegt; i++) {
                idmold(i) = imap(i, iseg);
            }
            if (imap(4, iseg) > imx) imx = imap(4, iseg);
            if (imap(6, iseg) > jmx) jmx = imap(6, iseg);
            if (iold(1, iseg) < iomn) iomn = iold(1, iseg);
            if (iold(2, iseg) > iomx) iomx = iold(2, iseg);
            if (iold(3, iseg) < jomn) jomn = iold(3, iseg);
            if (iold(4, iseg) > jomx) jomx = iold(4, iseg);
        }
        iseg = nxtseg(iseg);
    }
    // 199 continue

    // Allocate new segment iseg for block i1
    cmn_segment.nsgtop = cmn_segment.nsgtop + 1;
    iseg = cmn_segment.nsgtop;
    nxtseg(iseg) = idseg(i1);
    idseg(i1) = iseg;
    nseg(i1) = nseg(i1) + 1;
    idbloc(iseg) = i1;

    // Allocate new segment jseg for block i2
    cmn_segment.nsgtop = cmn_segment.nsgtop + 1;
    int jseg = cmn_segment.nsgtop;
    nxtseg(jseg) = idseg(i2);
    idseg(i2) = jseg;
    nseg(i2) = nseg(i2) + 1;
    idbloc(jseg) = i2;

    // Link the two interface segments
    intrfc(iseg) = jseg;
    intrfc(jseg) = iseg;

    // Set iold for iseg
    iold(1, iseg) = iomn;
    iold(2, iseg) = iomx;
    iold(3, iseg) = jomn;
    iold(4, iseg) = jomx;

    // Set iold for jseg
    iold(1, jseg) = iomn;
    iold(2, jseg) = iomx;
    iold(3, jseg) = jomn;
    iold(4, jseg) = jomx;

    // Set imap entries 1-12 for iseg
    imap(1,  iseg) = 1;
    imap(2,  iseg) = imn + 1;
    imap(3,  iseg) = 1;
    imap(4,  iseg) = imx;
    imap(5,  iseg) = 1;
    imap(6,  iseg) = jmx;
    imap(7,  iseg) = i2;
    imap(8,  iseg) = imn;
    imap(9,  iseg) = 1;
    imap(10, iseg) = imx;
    imap(11, iseg) = 1;
    imap(12, iseg) = jmx;

    // Set imap entries 1-12 for jseg
    imap(1,  jseg) = 1;
    imap(2,  jseg) = imn;
    imap(3,  jseg) = 1;
    imap(4,  jseg) = imx;
    imap(5,  jseg) = 1;
    imap(6,  jseg) = jmx;
    imap(7,  jseg) = i1;
    imap(8,  jseg) = imn + 1;
    imap(9,  jseg) = 1;
    imap(10, jseg) = imx;
    imap(11, jseg) = 1;
    imap(12, jseg) = jmx;

    // do 200 i=13,20: zero out imap entries 13-20 for both segments
    for (int i = 13; i <= 20; i++) {
        imap(i, iseg) = 0;
        imap(i, jseg) = 0;
    }

    // do 201 i=1,msegt: copy xdmold to xdum for both segments
    for (int i = 1; i <= msegt; i++) {
        xdum(i, iseg) = xdmold(i);
        xdum(i, jseg) = xdmold(i);
    }

    // Set imap entries 21-23 for both segments
    imap(21, iseg) = 0;   // no bc2000 series data on new interface
    imap(21, jseg) = 0;
    imap(22, iseg) = 1;   // don't turn off force on new interface
    imap(22, jseg) = 1;
    imap(23, iseg) = idmold(23);  // preserve iovrlp on new interface
    imap(23, jseg) = idmold(23);

    val(iseg) = 0.0;
    val(jseg) = 0.0;

    return;
}

} // namespace mkintr_ns
