// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "blkmax.h"

namespace blkmax_ns {

void blkmax(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, double& resmax, int& jm, int& km, int& im)
{
    resmax = res(1,1,1,1);
    jm     = 1;
    km     = 1;
    im     = 1;
    for (int j = 1; j <= jdim; j++) {
    for (int k = 1; k <= kdim; k++) {
    for (int i = 1; i <= idim-1; i++) {
        if ((float)res(j,k,i,1) > (float)resmax) {
            resmax = res(j,k,i,1);
            jm = j;
            km = k;
            im = i;
        }
    }
    }
    }
    return;
}

} // namespace blkmax_ns
