// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "setblk.h"

namespace setblk_ns {

void setblk(FortranArray3DRef<double> blank, int& jdim, int& kdim, int& idim, int& nbl)
{
    // Purpose: Initialize the blank array.
    for (int i = 1; i <= idim; i++) {
        for (int j = 1; j <= jdim; j++) {
            for (int k = 1; k <= kdim; k++) {
                blank(j, k, i) = 1.0;
            }
        }
    }
    return;
}

} // namespace setblk_ns
