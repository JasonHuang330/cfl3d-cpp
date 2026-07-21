// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "parser.h"

namespace parser_ns {

void parser(char (&inpstr)[210], int& npos, int& lc1, int& lc2, int& lcl, int& iflg)
{
    const int linp = 210;

    if (iflg == 1) {
        // normal parsing operation
        lc1 = (lc2 + 1 > 1) ? lc2 + 1 : 1;

        // loop over the sub-strings
        for (int ipos = 1; ipos <= npos; ipos++) {
            if (ipos > 1) lc1 = lc2 + 1;

            // find the first non-blank position
            int lc = lc1;
            for (lc = lc1; lc <= linp; lc++) {
                if (inpstr[lc - 1] != ' ') break;
            }
            // label 1100: continue
            lc1 = lc;

            // find the last non-blank position
            for (lc = lc1; lc <= linp; lc++) {
                if (inpstr[lc - 1] == ' ') break;
            }
            // label 1600: continue
            lc2 = lc - 1;
        }
    } else {
        // find lc2 of the last sub-string (lcl)
        int lc = linp;
        for (lc = linp; lc >= 1; lc--) {
            if (inpstr[lc - 1] != ' ') break;
        }
        // label 2100: continue
        lcl = lc;
    }

    return;
}

} // namespace parser_ns
