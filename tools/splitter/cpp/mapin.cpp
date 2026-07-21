// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "mapin.h"
#include "runtime/fortran_io.h"
#include "termn8.h"
#include <cstdio>
#include <cstring>
#include <algorithm>

namespace mapin_ns {

void mapin(int& imapfil, int& iwrit, int& nbloc, FortranArray1DRef<int> ni, FortranArray1DRef<int> nj, FortranArray1DRef<int> nk, int& mbloc, int& msegt, int& mtot, FortranArray2DRef<int> imap, FortranArray1DRef<int> idbloc, FortranArray2DRef<int> ivisb, FortranArray2DRef<int> itrb, FortranArray1DRef<double> val, FortranArray2DRef<double> xdum, FortranArray2DRef<int> iold, FortranArray1DRef<int> nxtseg, FortranArray1DRef<int> intrfc, FortranArray1DRef<int> ipatch, FortranArray1DRef<int> nsubbl, FortranArray1DRef<int> idobl, FortranArray1DRef<int> nseg, FortranArray1DRef<int> idno, FortranArray2DRef<int> ijk, FortranArray1DRef<int> idseg, FortranArray1DRef<int> idnext, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    // Local variables
    int idum;
    int ibloc, iseg, ns;
    int nbs, nfs, mns1, mxs1, mns2, mxs2;
    int match, nloc, loc;

    // read(imapfil,'(1x)')  -- skip one line
    {
        FILE* fp = fortran_get_unit(imapfil);
        char linebuf[256];
        if (fp) fgets(linebuf, sizeof(linebuf), fp);
    }

    // read(imapfil,*) idum
    fortran_read_list(imapfil, &idum);

    // if(idum.ne.nbloc) then
    if (idum != nbloc) {
        fortran_write_unit(iwrit, "  no. of blocks on map file %3d do not match the no. of blocks on grid file%3d\n", idum, nbloc);
        {
            int myid_stop = 0;
            int ierrflg = -6;
            termn8_ns::termn8(myid_stop, ierrflg, ibufdim, nbuf, bou, nou);
        }
        return;
    }

    // nbltop = nbloc
    cmn_block.nbltop = nbloc;
    // noldbl = nbloc
    cmn_oldbl.noldbl = nbloc;
    // nsgtop = 0
    cmn_segment.nsgtop = 0;

    // do 100 ibloc = 1,nbloc
    for (ibloc = 1; ibloc <= nbloc; ibloc++) {
        nsubbl(ibloc) = 1;
        idobl(ibloc) = ibloc;
        idno(ibloc) = ibloc;
        ijk(1, ibloc) = 1;
        ijk(2, ibloc) = ni(ibloc);
        ijk(3, ibloc) = 1;
        ijk(4, ibloc) = nj(ibloc);
        ijk(5, ibloc) = 1;
        ijk(6, ibloc) = nk(ibloc);
        idseg(ibloc) = cmn_segment.nsgtop + 1;
        idnext(ibloc) = 0;

        // read(imapfil,'(1x)')
        {
            FILE* fp = fortran_get_unit(imapfil);
            char linebuf[256];
            if (fp) fgets(linebuf, sizeof(linebuf), fp);
        }

        // read(imapfil,'(a79)') title1
        {
            FILE* fp = fortran_get_unit(imapfil);
            char linebuf[256];
            if (fp) {
                if (fgets(linebuf, sizeof(linebuf), fp)) {
                    // Remove trailing newline
                    int len = (int)strlen(linebuf);
                    if (len > 0 && linebuf[len-1] == '\n') { linebuf[len-1] = ' '; len--; }
                    // Copy up to 79 chars into title1 (char[80])
                    memset(cmn_titles.title1, ' ', 80);
                    int copy_len = len < 79 ? len : 79;
                    memcpy(cmn_titles.title1, linebuf, copy_len);
                }
            }
        }

        // read(imapfil,*) nseg(ibloc),ivisb(1,ibloc),ivisb(2,ibloc),
        //                 ivisb(3,ibloc),(itrb(n,ibloc),n=1,7)
        // This is a single record with 4 scalars + 7 implied-DO values = 11 values
        fortran_read_list(imapfil,
            &nseg(ibloc), &ivisb(1, ibloc), &ivisb(2, ibloc), &ivisb(3, ibloc),
            &itrb(1, ibloc), &itrb(2, ibloc), &itrb(3, ibloc), &itrb(4, ibloc),
            &itrb(5, ibloc), &itrb(6, ibloc), &itrb(7, ibloc));

        // read(imapfil,*)
        {
            FILE* fp = fortran_get_unit(imapfil);
            char linebuf[256];
            if (fp) fgets(linebuf, sizeof(linebuf), fp);
        }

        // read(imapfil,*) (ivisb(n,ibloc),n=4,15)
        fortran_read_list(imapfil,
            &ivisb(4, ibloc), &ivisb(5, ibloc), &ivisb(6, ibloc), &ivisb(7, ibloc),
            &ivisb(8, ibloc), &ivisb(9, ibloc), &ivisb(10, ibloc), &ivisb(11, ibloc),
            &ivisb(12, ibloc), &ivisb(13, ibloc), &ivisb(14, ibloc), &ivisb(15, ibloc));

        // read(imapfil,*) (ivisb(n,ibloc),n=16,27)
        fortran_read_list(imapfil,
            &ivisb(16, ibloc), &ivisb(17, ibloc), &ivisb(18, ibloc), &ivisb(19, ibloc),
            &ivisb(20, ibloc), &ivisb(21, ibloc), &ivisb(22, ibloc), &ivisb(23, ibloc),
            &ivisb(24, ibloc), &ivisb(25, ibloc), &ivisb(26, ibloc), &ivisb(27, ibloc));

        // read(imapfil,*) (ivisb(n,ibloc),n=28,34)
        fortran_read_list(imapfil,
            &ivisb(28, ibloc), &ivisb(29, ibloc), &ivisb(30, ibloc), &ivisb(31, ibloc),
            &ivisb(32, ibloc), &ivisb(33, ibloc), &ivisb(34, ibloc));

        // read(imapfil,*)
        {
            FILE* fp = fortran_get_unit(imapfil);
            char linebuf[256];
            if (fp) fgets(linebuf, sizeof(linebuf), fp);
        }

        // read(imapfil,'(a79)') title2
        {
            FILE* fp = fortran_get_unit(imapfil);
            char linebuf[256];
            if (fp) {
                if (fgets(linebuf, sizeof(linebuf), fp)) {
                    int len = (int)strlen(linebuf);
                    if (len > 0 && linebuf[len-1] == '\n') { linebuf[len-1] = ' '; len--; }
                    memset(cmn_titles.title2, ' ', 80);
                    int copy_len = len < 79 ? len : 79;
                    memcpy(cmn_titles.title2, linebuf, copy_len);
                }
            }
        }

        // read(imapfil,'(a79)') title3
        {
            FILE* fp = fortran_get_unit(imapfil);
            char linebuf[256];
            if (fp) {
                if (fgets(linebuf, sizeof(linebuf), fp)) {
                    int len = (int)strlen(linebuf);
                    if (len > 0 && linebuf[len-1] == '\n') { linebuf[len-1] = ' '; len--; }
                    memset(cmn_titles.title3, ' ', 80);
                    int copy_len = len < 79 ? len : 79;
                    memcpy(cmn_titles.title3, linebuf, copy_len);
                }
            }
        }

        // do 101 iseg = 1,nseg(ibloc)
        for (iseg = 1; iseg <= nseg(ibloc); iseg++) {
            cmn_segment.nsgtop = cmn_segment.nsgtop + 1;
            int nsgtop = cmn_segment.nsgtop;
            idbloc(nsgtop) = ibloc;
            intrfc(nsgtop) = 0;

            // read(imapfil,'(1x)')
            {
                FILE* fp = fortran_get_unit(imapfil);
                char linebuf[256];
                if (fp) fgets(linebuf, sizeof(linebuf), fp);
            }

            // read(imapfil,*) idum,idum,(imap(n,nsgtop),n= 1,12)
            fortran_read_list(imapfil,
                &idum, &idum,
                &imap(1, nsgtop), &imap(2, nsgtop), &imap(3, nsgtop), &imap(4, nsgtop),
                &imap(5, nsgtop), &imap(6, nsgtop), &imap(7, nsgtop), &imap(8, nsgtop),
                &imap(9, nsgtop), &imap(10, nsgtop), &imap(11, nsgtop), &imap(12, nsgtop));

            // read(imapfil,*) (imap(n,nsgtop),n=13,20)
            fortran_read_list(imapfil,
                &imap(13, nsgtop), &imap(14, nsgtop), &imap(15, nsgtop), &imap(16, nsgtop),
                &imap(17, nsgtop), &imap(18, nsgtop), &imap(19, nsgtop), &imap(20, nsgtop));

            // read(imapfil,*)
            {
                FILE* fp = fortran_get_unit(imapfil);
                char linebuf[256];
                if (fp) fgets(linebuf, sizeof(linebuf), fp);
            }

            // read(imapfil,*)
            {
                FILE* fp = fortran_get_unit(imapfil);
                char linebuf[256];
                if (fp) fgets(linebuf, sizeof(linebuf), fp);
            }

            // read(imapfil,*) (xdum(n,nsgtop),n=1,5)
            fortran_read_list(imapfil,
                &xdum(1, nsgtop), &xdum(2, nsgtop), &xdum(3, nsgtop),
                &xdum(4, nsgtop), &xdum(5, nsgtop));

            // read(imapfil,*) (xdum(n,nsgtop),n=6,10)
            fortran_read_list(imapfil,
                &xdum(6, nsgtop), &xdum(7, nsgtop), &xdum(8, nsgtop),
                &xdum(9, nsgtop), &xdum(10, nsgtop));

            // read(imapfil,*) (xdum(n,nsgtop),n=11,15)
            fortran_read_list(imapfil,
                &xdum(11, nsgtop), &xdum(12, nsgtop), &xdum(13, nsgtop),
                &xdum(14, nsgtop), &xdum(15, nsgtop));

            // read(imapfil,*) (xdum(n,nsgtop),n=16,20)
            fortran_read_list(imapfil,
                &xdum(16, nsgtop), &xdum(17, nsgtop), &xdum(18, nsgtop),
                &xdum(19, nsgtop), &xdum(20, nsgtop));

            // read(imapfil,*) (xdum(n,nsgtop),n=21,25)
            fortran_read_list(imapfil,
                &xdum(21, nsgtop), &xdum(22, nsgtop), &xdum(23, nsgtop),
                &xdum(24, nsgtop), &xdum(25, nsgtop));

            // read(imapfil,*) (xdum(n,nsgtop),n=26,26)
            fortran_read_list(imapfil, &xdum(26, nsgtop));

            // read(imapfil,*) (imap(n,nsgtop),n=21,23)
            fortran_read_list(imapfil,
                &imap(21, nsgtop), &imap(22, nsgtop), &imap(23, nsgtop));

            // read(imapfil,*)
            {
                FILE* fp = fortran_get_unit(imapfil);
                char linebuf[256];
                if (fp) fgets(linebuf, sizeof(linebuf), fp);
            }

            // nxtseg(nsgtop) = nsgtop + 1
            nxtseg(nsgtop) = nsgtop + 1;
            // iold(1,nsgtop) = imap(3,nsgtop)
            iold(1, nsgtop) = imap(3, nsgtop);
            // iold(2,nsgtop) = imap(4,nsgtop)
            iold(2, nsgtop) = imap(4, nsgtop);
            // iold(3,nsgtop) = imap(5,nsgtop)
            iold(3, nsgtop) = imap(5, nsgtop);
            // iold(4,nsgtop) = imap(6,nsgtop)
            iold(4, nsgtop) = imap(6, nsgtop);

        } // end do 101

        // nxtseg(nsgtop) = 0
        nxtseg(cmn_segment.nsgtop) = 0;

    } // end do 100

    // do 200 ibloc=1,nbloc
    for (ibloc = 1; ibloc <= nbloc; ibloc++) {
        loc = idseg(ibloc);
        // do 201 iseg = 1,nseg(ibloc)
        for (iseg = 1; iseg <= nseg(ibloc); iseg++) {
            if ((imap(1, loc) == 0) || (imap(1, loc) == 1)) {
                nbs  = imap(7, loc);
                nfs  = imap(8, loc);
                mns1 = std::min(imap(9, loc),  imap(10, loc));
                mxs1 = std::max(imap(9, loc),  imap(10, loc));
                mns2 = std::min(imap(11, loc), imap(12, loc));
                mxs2 = std::max(imap(11, loc), imap(12, loc));
                nfs  = abs(nfs);
                match = 0;
                nloc = idseg(nbs);
                // do 310 ns=1,nseg(nbs)
                for (ns = 1; ns <= nseg(nbs); ns++) {
                    if ((imap(1, nloc) == 0 || imap(1, nloc) == 1)
                        && imap(2, nloc) == nfs) {
                        if ((imap(3, nloc) == mns1 &&
                             imap(4, nloc) == mxs1) &&
                            (imap(5, nloc) == mns2 &&
                             imap(6, nloc) == mxs2)) {
                            match = nloc;
                        }
                    }
                    nloc = nxtseg(nloc);
                } // end do 310
                if (match == 0) {
                    fortran_write_unit(6, "  segment %3d of block %3d does not match target info\n", iseg, ibloc);
                    fortran_write_unit(6, " %d %d %d %d %d %d\n", nbs, nfs, mns1, mxs1, mns2, mxs2);
                    {
                        int myid_stop = 0;
                        int ierrflg = -6;
                        termn8_ns::termn8(myid_stop, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                }
                intrfc(loc) = match;
            }
            loc = nxtseg(loc);
        } // end do 201
    } // end do 200

    return;
}

} // namespace mapin_ns
