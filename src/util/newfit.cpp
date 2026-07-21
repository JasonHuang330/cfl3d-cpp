// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "newfit.h"
#include "runtime/common_blocks.h"

namespace newfit_ns {

void newfit(int& ifits, int& ifit, int& icount)
{
//***********************************************************************
//     Purpose:  Determine new polynomial fit for cells with stubborn
//     convergence.
//***********************************************************************

    // common /hist/ histry(4)
    float* histry = cmn_hist.histry;  // 0-based: histry[n-1] = Fortran histry(n)

    if (icount == 1) {
        histry[0] = (float)ifits;  // histry(1) = ifits
        int iff = ifits;
        for (int n = 2; n <= 4; n++) {
            if (iff == 4) iff = 0;
            iff = iff + 1;
            histry[n - 1] = (float)iff;  // histry(n) = iff
        }
    }

    // ifit = histry(icount+1)
    ifit = (int)histry[icount];  // histry(icount+1) → histry[icount] (0-based)
}

} // namespace newfit_ns
