// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// newfit.cpp — faithful C++ translation of newfit.F (CFL3D ronnie).
// Determine a new polynomial fit for cells with stubborn convergence.
//
// /hist/ histry(4) -> hist_.histry[0..3] (double, per Fortran declaration).
// NOTE: histry is REAL in Fortran; the integer fit codes are stored into it
// and read back with a real->int truncation (ifit = histry(icount+1)).
#include "ron_common.h"

void newfit(int ifits,int& ifit,int icount)
{
    if (icount == 1) {
        hist_.histry[1-1] = (double)ifits;
        int iff = ifits;
        for (int n = 2; n <= 4; ++n) {
            if (iff == 4) iff = 0;
            iff = iff + 1;
            hist_.histry[n-1] = (double)iff;
        }
    }
    // ifit = histry(icount+1): double -> int, Fortran truncates toward zero.
    ifit = f_int(hist_.histry[(icount+1)-1]);
    return;
}
