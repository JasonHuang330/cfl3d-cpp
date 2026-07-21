// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "mapout.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>

namespace mapout_ns {

void mapout(int& iwrit, int& nbloc, FortranArray1DRef<int> ni, FortranArray1DRef<int> nj, FortranArray1DRef<int> nk, int& mbloc, int& msegt, int& mtot, FortranArray2DRef<int> imap, FortranArray1DRef<int> idbloc, FortranArray2DRef<int> ivisb, FortranArray2DRef<int> itrb, FortranArray1DRef<double> val, FortranArray2DRef<double> xdum, FortranArray2DRef<int> iold, FortranArray1DRef<int> nxtseg, FortranArray1DRef<int> intrfc, FortranArray1DRef<int> ipatch, FortranArray1DRef<int> nsubbl, FortranArray1DRef<int> idobl, FortranArray1DRef<int> nseg, FortranArray1DRef<int> idno, FortranArray2DRef<int> ijk, FortranArray1DRef<int> idseg, FortranArray1DRef<int> idnext)
{
    // Local variables
    int nbl, n, n2, ibl, ics, iseg, ityp, jbl;
    int imax, jmax, kmax, mx1, mx2, jface;
    int nn, ibloc, loc, l;

    // Loop 100: renumber blocks
    nbl = 0;
    for (n = 1; n <= cmn_oldbl.noldbl; n++) {
        ibl = idobl(n);
        for (n2 = 1; n2 <= nsubbl(n); n2++) {
            nbl = nbl + 1;
            idno(ibl) = nbl;
            ibl = idnext(ibl);
        }
    }

    // Loop 110: update imap entries based on type
    for (n = 1; n <= cmn_oldbl.noldbl; n++) {
        ibl = idobl(n);
        for (n2 = 1; n2 <= nsubbl(n); n2++) {
            ics = idseg(ibl);
            for (iseg = 1; iseg <= nseg(ibl); iseg++) {
                ityp = imap(1, ics);
                if (0 == ityp || 1 == ityp) {
                    // jbl = imap(7,ics);
                    imap(7, intrfc(ics)) = idno(ibl);
                }
                if ((-1) == ityp) {
                    jbl = imap(7, ics);
                    imax = (ijk(2, jbl) + 1) - ijk(1, jbl);
                    jmax = (ijk(4, jbl) + 1) - ijk(3, jbl);
                    kmax = (ijk(6, jbl) + 1) - ijk(5, jbl);
                    imap(7, ics) = idno(jbl);
                    jface = imap(8, ics);
                    if (1 == jface || 2 == jface) {
                        mx1 = jmax;
                        mx2 = kmax;
                    } else if ((-1) == jface || (-2) == jface) {
                        mx1 = kmax;
                        mx2 = jmax;
                    } else if (3 == jface || 4 == jface) {
                        mx1 = kmax;
                        mx2 = imax;
                    } else if ((-3) == jface || (-4) == jface) {
                        mx1 = imax;
                        mx2 = kmax;
                    } else if (5 == jface || 6 == jface) {
                        mx1 = imax;
                        mx2 = jmax;
                    } else if ((-5) == jface || (-6) == jface) {
                        mx1 = jmax;
                        mx2 = imax;
                    }
                    if (imap(9,  ics) > mx1) imap(9,  ics) = mx1;
                    if (imap(10, ics) > mx1) imap(10, ics) = mx1;
                    if (imap(11, ics) > mx2) imap(11, ics) = mx2;
                    if (imap(12, ics) > mx2) imap(12, ics) = mx2;
                }
                ics = nxtseg(ics);
            }
            ibl = idnext(ibl);
        }
    }

    // write(iwrit,1000)  → format 1000: (1x,'nbloc')
    fortran_write_unit(iwrit, " nbloc\n");
    // write(iwrit,'(1x,i5)') nbl
    fortran_write_unit(iwrit, " %5d\n", nbl);

    // write(6,2000)  → format 2000: (/,/,' new block  old block  i0  i1  j0  j1  k0  k1')
    fortran_write_unit(6, "\n\n new block  old block  i0  i1  j0  j1  k0  k1\n");

    // Loop 200
    for (nn = 1; nn <= cmn_oldbl.noldbl; nn++) {
        ibloc = idobl(nn);
        for (n2 = 1; n2 <= nsubbl(nn); n2++) {
            // write(6,2001) idno(ibloc),nn,(ijk(l,ibloc),l=1,6)
            // format 2001: (2x,i4,7x,i4,4x,6i4)
            fortran_write_unit(6, "  %4d       %4d    %4d%4d%4d%4d%4d%4d\n",
                idno(ibloc), nn,
                ijk(1, ibloc), ijk(2, ibloc), ijk(3, ibloc),
                ijk(4, ibloc), ijk(5, ibloc), ijk(6, ibloc));

            // write(iwrit,'(1x)')
            fortran_write_unit(iwrit, " \n");

            // write(iwrit,'(a79)') title1
            // title1 is char[80], print first 79 chars
            fortran_write_unit(iwrit, "%-79.79s\n", cmn_titles.title1);

            // write(iwrit,1002) nseg(ibloc),ivisb(1,ibloc),ivisb(2,ibloc),ivisb(3,ibloc),(itrb(n,ibloc),n=1,7)
            // format 1002: (2i5,i6,11i5)
            fortran_write_unit(iwrit, "%5d%5d%6d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d\n",
                nseg(ibloc), ivisb(1, ibloc), ivisb(2, ibloc), ivisb(3, ibloc),
                itrb(1, ibloc), itrb(2, ibloc), itrb(3, ibloc), itrb(4, ibloc),
                itrb(5, ibloc), itrb(6, ibloc), itrb(7, ibloc));

            // write(iwrit,'("begin non-standard map data")')
            fortran_write_unit(iwrit, "begin non-standard map data\n");

            // write(iwrit,'(12i6)') (ivisb(n,ibloc),n=4,15)
            fortran_write_unit(iwrit, "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
                ivisb(4,  ibloc), ivisb(5,  ibloc), ivisb(6,  ibloc),
                ivisb(7,  ibloc), ivisb(8,  ibloc), ivisb(9,  ibloc),
                ivisb(10, ibloc), ivisb(11, ibloc), ivisb(12, ibloc),
                ivisb(13, ibloc), ivisb(14, ibloc), ivisb(15, ibloc));

            // write(iwrit,'(12i6)') (ivisb(n,ibloc),n=16,27)
            fortran_write_unit(iwrit, "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
                ivisb(16, ibloc), ivisb(17, ibloc), ivisb(18, ibloc),
                ivisb(19, ibloc), ivisb(20, ibloc), ivisb(21, ibloc),
                ivisb(22, ibloc), ivisb(23, ibloc), ivisb(24, ibloc),
                ivisb(25, ibloc), ivisb(26, ibloc), ivisb(27, ibloc));

            // write(iwrit,'(12i6)') (ivisb(n,ibloc),n=28,34)
            // Only 7 values (28..34), but format says 12i6 — Fortran fills remaining with spaces
            // Actually Fortran list-directed: 7 values printed with 12i6 → 7 fields then newline
            fortran_write_unit(iwrit, "%6d%6d%6d%6d%6d%6d%6d\n",
                ivisb(28, ibloc), ivisb(29, ibloc), ivisb(30, ibloc),
                ivisb(31, ibloc), ivisb(32, ibloc), ivisb(33, ibloc),
                ivisb(34, ibloc));

            // write(iwrit,'("end non-standard map data")')
            fortran_write_unit(iwrit, "end non-standard map data\n");

            // write(iwrit,'(a79)') title2
            fortran_write_unit(iwrit, "%-79.79s\n", cmn_titles.title2);

            // write(iwrit,'(a79)') title3
            fortran_write_unit(iwrit, "%-79.79s\n", cmn_titles.title3);

            loc = idseg(ibloc);

            // Loop 202
            for (iseg = 1; iseg <= nseg(ibloc); iseg++) {
                // write(iwrit,'(1x)')
                fortran_write_unit(iwrit, " \n");

                // write(iwrit,1002) idno(ibloc),iseg,(imap(n,loc),n=1,12)
                // format 1002: (2i5,i6,11i5)
                fortran_write_unit(iwrit, "%5d%5d%6d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d\n",
                    idno(ibloc), iseg,
                    imap(1,  loc), imap(2,  loc), imap(3,  loc),
                    imap(4,  loc), imap(5,  loc), imap(6,  loc),
                    imap(7,  loc), imap(8,  loc), imap(9,  loc),
                    imap(10, loc), imap(11, loc), imap(12, loc));

                // write(iwrit,1003) (imap(n,loc),n=13,20),val(loc)
                // format 1003: (15x,8i5,f5.1)
                fortran_write_unit(iwrit, "               %5d%5d%5d%5d%5d%5d%5d%5d%5.1f\n",
                    imap(13, loc), imap(14, loc), imap(15, loc), imap(16, loc),
                    imap(17, loc), imap(18, loc), imap(19, loc), imap(20, loc),
                    val(loc));

                // write(iwrit,'(15x,5f5.1)') 0.0,0.0,0.0,1.0,0.0
                fortran_write_unit(iwrit, "               %5.1f%5.1f%5.1f%5.1f%5.1f\n",
                    0.0, 0.0, 0.0, 1.0, 0.0);

                // write(iwrit,'("begin non-standard map data")')
                fortran_write_unit(iwrit, "begin non-standard map data\n");

                // write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=1,5)
                fortran_write_unit(iwrit, "               %10.4f%10.4f%10.4f%10.4f%10.4f\n",
                    xdum(1, loc), xdum(2, loc), xdum(3, loc), xdum(4, loc), xdum(5, loc));

                // write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=6,10)
                fortran_write_unit(iwrit, "               %10.4f%10.4f%10.4f%10.4f%10.4f\n",
                    xdum(6, loc), xdum(7, loc), xdum(8, loc), xdum(9, loc), xdum(10, loc));

                // write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=11,15)
                fortran_write_unit(iwrit, "               %10.4f%10.4f%10.4f%10.4f%10.4f\n",
                    xdum(11, loc), xdum(12, loc), xdum(13, loc), xdum(14, loc), xdum(15, loc));

                // write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=16,20)
                fortran_write_unit(iwrit, "               %10.4f%10.4f%10.4f%10.4f%10.4f\n",
                    xdum(16, loc), xdum(17, loc), xdum(18, loc), xdum(19, loc), xdum(20, loc));

                // write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=21,25)
                fortran_write_unit(iwrit, "               %10.4f%10.4f%10.4f%10.4f%10.4f\n",
                    xdum(21, loc), xdum(22, loc), xdum(23, loc), xdum(24, loc), xdum(25, loc));

                // write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=26,26)
                // Only 1 value with 5f10.4 format — Fortran prints 1 value then newline
                fortran_write_unit(iwrit, "               %10.4f\n",
                    xdum(26, loc));

                // write(iwrit,'(15x,8i5)') (imap(n,loc),n=21,23)
                // Only 3 values with 8i5 format — Fortran prints 3 values then newline
                fortran_write_unit(iwrit, "               %5d%5d%5d\n",
                    imap(21, loc), imap(22, loc), imap(23, loc));

                // write(iwrit,'("end non-standard map data")')
                fortran_write_unit(iwrit, "end non-standard map data\n");

                loc = nxtseg(loc);
            } // end loop 202

            ibloc = idnext(ibloc);
        } // end loop 201
    } // end loop 200

    // write(iwrit,'(1x)')
    fortran_write_unit(iwrit, " \n");

    return;
}

} // namespace mapout_ns
