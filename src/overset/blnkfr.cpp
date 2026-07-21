// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "blnkfr.h"

namespace blnkfr_ns {

void blnkfr(int& nbl, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg,
            FortranArray1DRef<int> jjbg, FortranArray2DRef<int> ibpntsg,
            FortranArray1DRef<int> lbg, int& iitot, FortranArray3DRef<double> blank,
            int& jdim, int& kdim, int& idim, int& maxbl, double& blnkval)
{
    // Purpose: Undo/redo blanking at fringe points (but not hole points)
    // blnkval = 1...set blank=1 at fringe points for plotting
    //           0...set blank=0 (normal setting for fringe points)

    int lsta, lend;

    lsta = lbg(nbl);
    lend = lsta - 1;
    if (ibpntsg(nbl, 1) > 0) {
        lend = lsta + ibpntsg(nbl, 1) - 1;
        for (int l = lsta; l <= lend; l++) {
            blank(jjbg(l), kkbg(l), iibg(l)) = blnkval;
        }
    }

    return;
}

} // namespace blnkfr_ns
