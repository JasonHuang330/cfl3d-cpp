// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "parser.h"
#include <algorithm>

namespace parser_ns {

void parser(char (&inpstr)[210], int& npos, int& lc1, int& lc2, int& lcl, int& iflg)
{
    // Purpose: parse 210-character data
    //
    // inpstr = input character string
    // npos   = sub-string position to find after column lc2-on-input
    // lc1    = first column of the sub-string (on output)
    // lc2    = last  column of the sub-string (on output)
    //        = last  column of the previous sub-string (on input)
    // lcl    = last column of the last sub-string
    // iflg   = 1: do normal parsing: start from front and extract
    //         -1: only find lcl

    const int linp = 210;

    // normal parsing operation
    if (iflg == 1) {
        lc1 = std::max(lc2 + 1, 1);

        // loop over the sub-strings
        for (int ipos = 1; ipos <= npos; ipos++) {
            if (ipos > 1) lc1 = lc2 + 1;

            // find the first non-blank position
            int lc = lc1;
            for (lc = lc1; lc <= linp; lc++) {
                if (inpstr[lc - 1] != ' ') goto label1100;
            }
            label1100:
            lc1 = lc;

            // find the last non-blank position
            for (lc = lc1; lc <= linp; lc++) {
                if (inpstr[lc - 1] == ' ') goto label1600;
            }
            label1600:
            lc2 = lc - 1;
        }

    // find lc2 of the last sub-string (lcl)
    } else {
        int lc = linp;
        for (lc = linp; lc >= 1; lc--) {
            if (inpstr[lc - 1] != ' ') goto label2100;
        }
        label2100:
        lcl = lc;
    }

    return;
}

} // namespace parser_ns
