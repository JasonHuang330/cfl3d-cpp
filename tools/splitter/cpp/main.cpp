// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// ---------------------------------------------------------------------------
// CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
// Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
// with point-matched, patched, overset, or embedded connectivities. Both
// multigrid and mesh sequencing are available in time-accurate or
// steady-state modes.
//
// Copyright 2001 United States Government as represented by the Administrator
// of the National Aeronautics and Space Administration. All Rights Reserved.
//
// The CFL3D platform is licensed under the Apache License, Version 2.0
// (the "License"); you may not use this file except in compliance with the
// License. You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0.
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
// License for the specific language governing permissions and limitations
// under the License.
// ---------------------------------------------------------------------------
//
// $Id$
//
// Given a grid and a CFL3D input file (and corresponding ronnie
// input file, if applicable) and grid, this program will split the
// grid along user-specified grid lines to create a new grid file
// with more blocks, and will split and create a new CFL3D input
// file (and new ronnie input file) to correspond to the new, split
// grid. For the versions of CFL3D that compute sensitivity derivatives,
// this splitter can also split the "seed matrix" or grid sensitivity
// file.
//
// File formats: 1) grid file is PLOT3D GRID, multiblock style
//                  (formated or unformatted). NOTE: cfl3d style
//                  grid files are currently not supported!
//               2) grid sensitivity file (if applicable)
//                  is PLOT3D FUNCTION, multiblock style
//
// Original coding by E. B. Parlette for TLNS3D (Summer 1995)
// Modifications for CFL3D by R. T. Biedron (Winter 1998)
//
// Files and I/O unit numbers:
//   unit  5...main input file controling the splitting
//   unit 51...copy of unit 5 data (temporary)
//   unit  6...information output file
//   unit 30...input (unsplit) grid file
//   unit 10...input (unsplit) cfl3d input file
//   unit 15...input (unsplit) ronnie input file
//   unit 50...input (unsplit) grid sensitivity file
//   unit  7...output (split) cfl3d input file
//   unit  9...output (split) ronnie input file
//   unit 40...output (split) grid file
//   unit 60...output (split) grid sensitivity file
//   unit 20...input (unplit) tlns3d map file (temporary)
//   unit  8...output (split) tlns3d map file (temporary)
//   (the temporary files above are deleted upon program completion)

#include "runtime/fortran_array.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include "sizer.h"
#include "splitter.h"
#include "termn8.h"
#include <cstdio>
#include <cstring>

int main(int argc, char** argv) {
    // parameters
    const int mbloc0  = 2000;
    const int msplt0  = 200;
    const int ibufdim = 2000;
    const int nbuf    = 4;

    // character*80 inpstring
    char inpstring[81];

    // character*120 bou(ibufdim,nbuf)
    // Fortran column-major: bou(mm,ii) -> bou_storage[ii-1][mm-1]
    // Use plain C storage and wrap in FortranArray2DRef<char[120]>
    static char bou_storage[4][2000][120];

    // dimension nou(nbuf)
    int nou_storage[4];

    // Build FortranArray refs for bou and nou
    // FortranArray2DRef<char[120]>: (data, dim1, dim2, lb1, lb2)
    FortranArray2DRef<char[120]> bou(
        reinterpret_cast<char(*)[120]>(bou_storage),
        ibufdim, nbuf, 1, 1);
    FortranArray1DRef<int> nou(nou_storage, nbuf, 1);

    // common /unit5/ iunit5
    // common /zero/ iexp
    // (accessed via cmn_unit5.iunit5 and cmn_zero.iexp)

    // local variables
    int    ii, mm;
    int    icount;
    float  compare, add, x11;
    int    i, n;
    int    mbloc, maxseg, npmax, mxbli, intmax, nsub1, msegt, msegn, msplt, mtot;

    // Initialize bou and nou arrays
    for (ii = 1; ii <= nbuf; ii++) {
        nou(ii) = 0;
        for (mm = 1; mm <= ibufdim; mm++) {
            // bou(mm,ii) = ' '  (fill with spaces)
            memset(bou(mm, ii), ' ', 120);
        }
    }

    // If a command-line argument is given, open it as unit 5 (the main input file).
    // The Fortran program reads from unit 5 (stdin), but the C++ translation is
    // invoked as: splitter <inputfile.inp>
    if (argc > 1) {
        fortran_open_unit(5, argv[1], "r");
    }

    // Determine machine zero for use in setting tolerances
    // (10.**(-iexp) is machine zero)
    icount  = 0;
    compare = 1.0f;
    for (i = 1; i <= 20; i++) {
        icount = icount + 1;
        add = 1.0f;
        for (n = 1; n <= i; n++) {
            add = add * 0.1f;
        }
        x11 = compare + add;
        if (x11 == compare) {
            cmn_zero.iexp = i - 1;
            goto label_4010;
        }
    }
label_4010:;

    // Copy command-line splitter input file to a temp file
    // that is explicitly opened so that this temp file can be rewound
    cmn_unit5.iunit5 = 51;
    fortran_open_unit(cmn_unit5.iunit5, "tempz3y2x1.inp", "replace");

    // do n=1,9999: read(5,'(a80)',end=999) inpstring; write(iunit5,'(a80)') inpstring
    {
        FILE* fp_in  = fortran_get_unit(5);
        FILE* fp_out = fortran_get_unit(cmn_unit5.iunit5);
        for (n = 1; n <= 9999; n++) {
            // read(5,'(a80)',end=999) inpstring
            // Formatted read: read exactly 80 characters (one record/line)
            if (fgets(inpstring, sizeof(inpstring), fp_in) == nullptr) {
                goto label_999;
            }
            // Strip trailing newline, pad to 80 chars with spaces
            {
                int len = (int)strlen(inpstring);
                // Remove trailing newline if present
                if (len > 0 && inpstring[len - 1] == '\n') {
                    inpstring[len - 1] = ' ';
                    len--;
                }
                if (len > 0 && inpstring[len - 1] == '\r') {
                    inpstring[len - 1] = ' ';
                    len--;
                }
                // Pad to 80 chars with spaces
                for (int k = len; k < 80; k++) {
                    inpstring[k] = ' ';
                }
                inpstring[80] = '\0';
            }
            // write(iunit5,'(a80)') inpstring
            fprintf(fp_out, "%-80.80s\n", inpstring);
        }
    }
label_999:;
    // rewind(iunit5)
    rewind(fortran_get_unit(cmn_unit5.iunit5));

    // open(unit=99, file='splitter.error', form='formatted', status='unknown')
    fortran_open_unit(99, "splitter.error", "replace");

    // call sizer(...)
    sizer_ns::sizer(mbloc, maxseg, npmax, mxbli, intmax, nsub1, msegt,
                    msegn, msplt, mtot, const_cast<int&>(mbloc0),
                    const_cast<int&>(msplt0), nou, bou, const_cast<int&>(nbuf),
                    const_cast<int&>(ibufdim));

    // call splitter(...)
    splitter_ns::splitter(mbloc, maxseg, npmax, mxbli, intmax, nsub1, msegt,
                          msegn, msplt, mtot, nou, bou, const_cast<int&>(nbuf),
                          const_cast<int&>(ibufdim));

    // close(iunit5, status='delete')
    // Note: runtime does not support STATUS='delete'; just close the unit
    fortran_close_unit(cmn_unit5.iunit5);

    // call termn8(0, 0, ibufdim, nbuf, bou, nou)
    {
        int arg_myid_stop = 0;
        int arg_ierrflg   = 0;
        int arg_ibufdim   = ibufdim;
        int arg_nbuf      = nbuf;
        termn8_ns::termn8(arg_myid_stop, arg_ierrflg, arg_ibufdim, arg_nbuf,
                          bou, nou);
    }

    return 0;
}
