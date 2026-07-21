// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "mapsplt.h"
#include "ipadd.h"
#include "movseg.h"
#include "renmbr.h"
#include "spltsg.h"
#include "mkintr.h"
#include "ipfix.h"

namespace mapsplt_ns {

void mapsplt(int& nbl, int& ndir, int& indx, int& mbloc, int& msegt, int& mtot,
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
    int ix1, ix2, imn, imx, jmn, jmx;

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

    int ibl = idobl(nbl);
    int nsub = nsubbl(nbl);

    for (int n = 1; n <= nsub; n++) {
        if ((ijk(ix1, ibl) < indx) && (ijk(ix2, ibl) > indx)) {
            nsubbl(nbl) = nsubbl(nbl) + 1;
            cmn_block.nbltop = cmn_block.nbltop + 1;
            int nbnew = cmn_block.nbltop;
            idnext(nbnew) = idobl(nbl);
            idobl(nbl) = nbnew;
            nseg(nbnew) = 0;
            idseg(nbnew) = 0;

            for (int i = 1; i <= 6; i++) {
                ijk(i, nbnew) = ijk(i, ibl);
            }
            for (int i = 1; i <= msegt; i++) {
                ivisb(i, nbnew) = ivisb(i, ibl);
            }
            for (int i = 1; i <= 7; i++) {
                itrb(i, nbnew) = itrb(i, ibl);
            }

            ijk(ix1, nbnew) = indx;
            ijk(ix2, ibl) = indx;
            int local = (indx + 1) - ijk(ix1, ibl);
            int loc2  = (ijk(ix2, nbnew) + 1) - indx;
            int icnt = 1;
            ipatch(1) = 2;
            int iseg = idseg(ibl);

            // Label 900/999: while loop with goto semantics
            while (true) {
                if (iseg == 0) break; // goto 999

                if (imap(1, iseg) == -1) {
                    ipadd_ns::ipadd(imap(7, iseg), imap(8, iseg),
                                    mbloc, msegt, mtot,
                                    imap, idbloc, ivisb, itrb, val, xdum, iold,
                                    nxtseg, intrfc, ipatch, nsubbl, idobl,
                                    nseg, idno, ijk, idseg, idnext);
                }

                int nft = imap(2, iseg);
                if (nft == imn) {
                    if (imap(13, iseg) != 0) {
                        if (imap(19, iseg) > local) imap(19, iseg) = local;
                        if (imap(20, iseg) > local) imap(20, iseg) = local;
                    }
                    iseg = nxtseg(iseg);
                } else if (nft == imx) {
                    if (imap(13, iseg) != 0) {
                        if (imap(19, iseg) > loc2) imap(19, iseg) = loc2;
                        if (imap(20, iseg) > loc2) imap(20, iseg) = loc2;
                    }
                    int itmp = nxtseg(iseg);
                    movseg_ns::movseg(iseg, ibl, nbnew,
                                      mbloc, msegt, mtot,
                                      imap, idbloc, ivisb, itrb, val, xdum, iold,
                                      nxtseg, intrfc, ipatch, nsubbl, idobl,
                                      nseg, idno, ijk, idseg, idnext);
                    iseg = itmp;
                } else if ((nft == jmn) || (nft == jmx)) {
                    if (imap(6, iseg) <= local) {
                        iseg = nxtseg(iseg);
                    } else if (imap(5, iseg) >= local) {
                        int itmp = nxtseg(iseg);
                        movseg_ns::movseg(iseg, ibl, nbnew,
                                          mbloc, msegt, mtot,
                                          imap, idbloc, ivisb, itrb, val, xdum, iold,
                                          nxtseg, intrfc, ipatch, nsubbl, idobl,
                                          nseg, idno, ijk, idseg, idnext);
                        int lit5 = 5, lit6 = 6;
                        renmbr_ns::renmbr(iseg, lit5, lit6, local,
                                          mbloc, msegt, mtot,
                                          imap, idbloc, ivisb, itrb, val, xdum, iold,
                                          nxtseg, intrfc, ipatch, nsubbl, idobl,
                                          nseg, idno, ijk, idseg, idnext);
                        iseg = itmp;
                    } else {
                        int lit5 = 5, lit6 = 6;
                        spltsg_ns::spltsg(iseg, lit5, lit6, local, ibl, nbnew,
                                          mbloc, msegt, mtot,
                                          imap, idbloc, ivisb, itrb, val, xdum, iold,
                                          nxtseg, intrfc, ipatch, nsubbl, idobl,
                                          nseg, idno, ijk, idseg, idnext);
                        iseg = idseg(ibl);
                    }
                } else {
                    if (imap(4, iseg) <= local) {
                        iseg = nxtseg(iseg);
                    } else if (imap(3, iseg) >= local) {
                        int itmp = nxtseg(iseg);
                        movseg_ns::movseg(iseg, ibl, nbnew,
                                          mbloc, msegt, mtot,
                                          imap, idbloc, ivisb, itrb, val, xdum, iold,
                                          nxtseg, intrfc, ipatch, nsubbl, idobl,
                                          nseg, idno, ijk, idseg, idnext);
                        int lit3 = 3, lit4 = 4;
                        renmbr_ns::renmbr(iseg, lit3, lit4, local,
                                          mbloc, msegt, mtot,
                                          imap, idbloc, ivisb, itrb, val, xdum, iold,
                                          nxtseg, intrfc, ipatch, nsubbl, idobl,
                                          nseg, idno, ijk, idseg, idnext);
                        iseg = itmp;
                    } else {
                        int lit3 = 3, lit4 = 4;
                        spltsg_ns::spltsg(iseg, lit3, lit4, local, ibl, nbnew,
                                          mbloc, msegt, mtot,
                                          imap, idbloc, ivisb, itrb, val, xdum, iold,
                                          nxtseg, intrfc, ipatch, nsubbl, idobl,
                                          nseg, idno, ijk, idseg, idnext);
                        iseg = idseg(ibl);
                    }
                }

                icnt = icnt + 1;
                // goto 900 (continue while loop)
            }
            // label 999: continue

            mkintr_ns::mkintr(imn, ibl, nbnew,
                              mbloc, msegt, mtot,
                              imap, idbloc, ivisb, itrb, val, xdum, iold,
                              nxtseg, intrfc, ipatch, nsubbl, idobl,
                              nseg, idno, ijk, idseg, idnext,
                              xdmold, idmold);
            ipfix_ns::ipfix(ibl, ndir, local, nbnew,
                            mbloc, msegt, mtot,
                            imap, idbloc, ivisb, itrb, val, xdum, iold,
                            nxtseg, intrfc, ipatch, nsubbl, idobl,
                            nseg, idno, ijk, idseg, idnext);
        }
        ibl = idnext(ibl);
    }

    return;
}

} // namespace mapsplt_ns
