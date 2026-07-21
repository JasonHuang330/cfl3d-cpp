// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "ipfix.h"
#include "mkseg.h"

namespace ipfix_ns {

void ipfix(int& ibl, int& ndir, int& indx, int& inewbl, int& mbloc, int& msegt, int& mtot,
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
    int ix1, ix2, imn, imx, jmn, jmx;

    // Set direction indices based on ndir
    if (ndir == 1) {
        ix1 = 1;
        ix2 = 2;
        imn = 1;
        imx = 2;
        jmn = 3;
        jmx = 4;
    } else if (ndir == 2) {
        ix1 = 3;
        ix2 = 4;
        imn = 3;
        imx = 4;
        jmn = 5;
        jmx = 6;
    } else {
        ix1 = 5;
        ix2 = 6;
        imn = 5;
        imx = 6;
        jmn = 1;
        jmx = 2;
    }

    int itop = ipatch(1);
    int i = 2;

    // Label 1: outer loop
    while (i < itop) {
        int j = i + 1;
        int jbl = ipatch(i);
        int jfc = ipatch(j);

        int is = idseg(jbl);

        // Label 2: inner loop over segment linked list
        while (is != 0) {
            if (imap(1, is) != -1)  goto label_997;
            if (imap(2, is) != jfc) goto label_997;
            if (imap(7, is) != ibl) goto label_997;

            {
                int nft = abs(imap(8, is));
                int jx1, jx2;

                if (imap(8, is) > 0) {
                    if (nft == jmn || nft == jmx) {
                        jx1 = 11;
                        jx2 = 12;
                    } else {
                        jx1 = 9;
                        jx2 = 10;
                    }
                } else {
                    if (nft == jmn || nft == jmx) {
                        jx1 = 9;
                        jx2 = 10;
                    } else {
                        jx1 = 11;
                        jx2 = 12;
                    }
                }

                if (imap(jx1, is) > imap(jx2, is)) {
                    int jtmp = jx1;
                    jx1 = jx2;
                    jx2 = jtmp;
                }

                if (nft == imn) goto label_997;

                if (nft == imx) {
                    imap(7, is) = inewbl;
                } else if (imap(jx1, is) >= indx) {
                    imap(7, is) = inewbl;
                    imap(jx1, is) = (imap(jx1, is) + 1) - indx;
                    imap(jx2, is) = (imap(jx2, is) + 1) - indx;
                } else if (imap(jx2, is) > indx) {
                    int ns_local = 0;
                    mkseg_ns::mkseg(is, jbl, ns_local, mbloc, msegt, mtot,
                                    imap, idbloc, ivisb, itrb, val, xdum, iold,
                                    nxtseg, intrfc, ipatch, nsubbl, idobl,
                                    nseg, idno, ijk, idseg, idnext);
                    imap(jx2, is) = indx;
                    imap(7, ns_local) = inewbl;
                    imap(jx1, ns_local) = 1;
                    imap(jx2, ns_local) = (imap(jx2, ns_local) + 1) - indx;
                }
            }

            label_997:
            is = nxtseg(is);
        }

        // Label 998:
        i = j + 1;
    }

    // Label 999:
    return;
}

} // namespace ipfix_ns
