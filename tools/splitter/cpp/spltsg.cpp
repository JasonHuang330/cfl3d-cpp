// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "spltsg.h"
#include "mkseg.h"
#include "splt2.h"
#include "renmbr.h"

namespace spltsg_ns {

void spltsg(int& iseg, int& ix1, int& ix2, int& indx, int& is, int& it, int& mbloc, int& msegt, int& mtot,
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
    int jseg = intrfc(iseg);
    if (jseg > 0) {
        int jx1, jx2, is1, is2;
        if (imap(8, iseg) > 0) {
            jx1 = ix1;
            jx2 = ix2;
            is1 = ix1 + 6;
            is2 = ix2 + 6;
        } else if (ix1 == 3) {
            jx1 = 5;
            jx2 = 6;
            is1 = 11;
            is2 = 12;
        } else {
            jx1 = 3;
            jx2 = 4;
            is1 = 9;
            is2 = 10;
        }
        int ia = (imap(is2, iseg) - imap(is1, iseg)) /
                 (imap(ix2, ns) - imap(ix1, iseg));
        int ib = imap(is1, iseg) - ia * imap(ix1, iseg);
        int jndx = ia * indx + ib;
        if (iseg == jseg) {
            if (jndx > imap(is1, iseg)) {
                imap(is2, iseg) = jndx;
                imap(is1, ns) = jndx;
                intrfc(ns) = ns;
            } else if (jndx == indx) {
                imap(is2, iseg) = jndx;
                imap(is1, ns) = jndx;
                intrfc(ns) = iseg;
                intrfc(iseg) = ns;
            } else if (jndx > indx) {
                splt2_ns::splt2(ns, ix1, ix2, jndx, it, it, mbloc, msegt, mtot,
                                imap, idbloc, ivisb, itrb, val, xdum, iold,
                                nxtseg, intrfc, ipatch, nsubbl, idobl,
                                nseg, idno, ijk, idseg, idnext);
                imap(ix1, ns) = 1;
                imap(ix2, ns) = jndx - indx + 1;
                imap(is1, ns) = imap(ix2, ns);
                imap(is2, ns) = imap(ix1, ns);
                intrfc(ns) = ns;
                imap(7, ns) = it;
                imap(8, ns) = imap(2, ns);
                ns = cmn_segment.nsgtop;
                intrfc(ns) = iseg;
                intrfc(iseg) = ns;
            } else if (jndx < indx) {
                splt2_ns::splt2(iseg, ix1, ix2, jndx, is, is, mbloc, msegt, mtot,
                                imap, idbloc, ivisb, itrb, val, xdum, iold,
                                nxtseg, intrfc, ipatch, nsubbl, idobl,
                                nseg, idno, ijk, idseg, idnext);
                int nsold = ns;
                ns = cmn_segment.nsgtop;
                imap(ix1, ns) = jndx;
                imap(ix2, ns) = indx;
                imap(is1, ns) = imap(ix2, ns);
                imap(is2, ns) = imap(ix1, ns);
                intrfc(ns) = ns;
                imap(7, ns) = is;
                imap(8, ns) = imap(2, ns);
                ns = nsold;
                intrfc(ns) = iseg;
                intrfc(iseg) = ns;
            }
            int one = 1;
            renmbr_ns::renmbr(iseg, ix1, ix2, one, mbloc, msegt, mtot,
                              imap, idbloc, ivisb, itrb, val, xdum, iold,
                              nxtseg, intrfc, ipatch, nsubbl, idobl,
                              nseg, idno, ijk, idseg, idnext);
        } else {
            imap(is2, iseg) = jndx;
            imap(is1, ns) = jndx;
            int ibl = idbloc(jseg);
            splt2_ns::splt2(jseg, jx1, jx2, jndx, ibl, ibl, mbloc, msegt, mtot,
                            imap, idbloc, ivisb, itrb, val, xdum, iold,
                            nxtseg, intrfc, ipatch, nsubbl, idobl,
                            nseg, idno, ijk, idseg, idnext);
            if (jndx > imap(is1, iseg)) {
                intrfc(ns) = idseg(ibl);
                intrfc(idseg(ibl)) = ns;
                intrfc(jseg) = iseg;
            } else {
                intrfc(ns) = jseg;
                intrfc(idseg(ibl)) = iseg;
                intrfc(jseg) = ns;
                intrfc(iseg) = idseg(ibl);
            }
            int one = 1;
            renmbr_ns::renmbr(iseg, ix1, ix2, one, mbloc, msegt, mtot,
                              imap, idbloc, ivisb, itrb, val, xdum, iold,
                              nxtseg, intrfc, ipatch, nsubbl, idobl,
                              nseg, idno, ijk, idseg, idnext);
        }
    }
    renmbr_ns::renmbr(ns, ix1, ix2, indx, mbloc, msegt, mtot,
                      imap, idbloc, ivisb, itrb, val, xdum, iold,
                      nxtseg, intrfc, ipatch, nsubbl, idobl,
                      nseg, idno, ijk, idseg, idnext);
}

} // namespace spltsg_ns
