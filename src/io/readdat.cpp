// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "readdat.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cstring>
#include <cstdio>

namespace readdat_ns {

void readdat(FortranArray4DRef<double> bcdata, int& mdim, int& ndim, char (&filname)[80], FortranArray1DRef<double> bcdat, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    // initialize data array for check in bc routine
    for (int l = 1; l <= 12; l++)
        for (int ip = 1; ip <= 2; ip++)
            for (int m = 1; m <= mdim; m++)
                for (int n = 1; n <= ndim; n++)
                    bcdata(m, n, ip, l) = -1.e15;

    // Check if filname != 'null'
    // filname is char[80], Fortran blank-padded; compare first 4 chars
    // Fortran: filname.ne.'null' — compare the 80-char string with 'null' (blank-padded)
    char filname_cstr[81];
    std::memcpy(filname_cstr, filname, 80);
    filname_cstr[80] = '\0';
    // Trim trailing spaces to get effective string
    int flen = 80;
    while (flen > 0 && filname_cstr[flen - 1] == ' ') flen--;
    filname_cstr[flen] = '\0';

    bool is_null = (std::strcmp(filname_cstr, "null") == 0);

    if (!is_null) {
        int iunit = 26;
        fortran_open_unit(iunit, filname_cstr, "r");

        fortran_write_unit(11, "\n");
        fortran_write_unit(11, " reading bc data file %-80.80s\n", filname);

        // read(iunit,*) — skip one record
        {
            FILE* fp = fortran_get_unit(iunit);
            char linebuf[256];
            fgets(linebuf, sizeof(linebuf), fp);
        }

        // read(iunit,*) mdum, ndum, np
        int mdum, ndum, np;
        fortran_read_list(iunit, &mdum, &ndum, &np);

        if (mdum != mdim || ndum != ndim) {
            fortran_write_unit(11, " stopping in readdat...data file has dimensions %6d%6d\n", mdum, ndum);
            fortran_write_unit(11, " dimensions set from input file are %6d%6d\n", mdim, ndim);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        if (np < 1 || np > 2) {
            fortran_write_unit(11, " stopping in readdat...data file must have 1 or 2 planes of data - there are %6d\n", np);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        // read(iunit,*) nvar
        int nvar;
        fortran_read_list(iunit, &nvar);

        if (nvar > 12 || nvar < 1) {
            fortran_write_unit(11, " stopping in readdat...no more than 12 variables are allowed - there are %6d\n", nvar);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        // read(iunit,*) ((((bcdata(m,n,ip,l),m=1,mdim),n=1,ndim),ip=1,np),l=1,nvar)
        // Nested implied-DO: innermost is m, then n, then ip, then l
        // Total values = mdim * ndim * np * nvar
        int total = mdim * ndim * np * nvar;
        FortranArray1D<double> tmp(total);
        fortran_read_list_array(iunit, tmp, 1, total);
        // Copy in Fortran column-major order: m varies fastest, then n, then ip, then l
        int idx = 1;
        for (int l = 1; l <= nvar; l++)
            for (int ip = 1; ip <= np; ip++)
                for (int n = 1; n <= ndim; n++)
                    for (int m = 1; m <= mdim; m++) {
                        bcdata(m, n, ip, l) = tmp(idx);
                        idx++;
                    }

        // add 2nd plane of data if not read in (duplicate 1st plane)
        if (np < 2) {
            for (int l = 1; l <= nvar; l++)
                for (int m = 1; m <= mdim; m++)
                    for (int n = 1; n <= ndim; n++)
                        bcdata(m, n, 2, l) = bcdata(m, n, 1, l);
        }

        rewind(fortran_get_unit(iunit));
        fortran_close_unit(iunit);

    } else {
        // filname == 'null': fill bcdata from bcdat
        for (int l = 1; l <= 12; l++)
            for (int ip = 1; ip <= 2; ip++)
                for (int m = 1; m <= mdim; m++)
                    for (int n = 1; n <= ndim; n++)
                        bcdata(m, n, ip, l) = bcdat(l);
    }
}

} // namespace readdat_ns
