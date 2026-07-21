// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "loadgr.h"

namespace loadgr_ns {

void loadgr(FortranArray1DRef<double> w, int& mgwk, int& lx, int& ly, int& lz, int& jindex,
            FortranArray2DRef<double> x, FortranArray2DRef<double> y, FortranArray2DRef<double> z,
            int& mdim, int& ndim, int& idimg, int& jdimg, int& kdimg)
{
    int i, j, k, l;

    //
    // patch surface is an i=constant surface
    //
    if (jindex / 10 == 1) {
        if (jindex == 11) {
            i = 1;
        } else {
            i = idimg;
        }
        for (j = 1; j <= jdimg; j++) {
            for (k = 1; k <= kdimg; k++) {
                l      = (i - 1) * jdimg * kdimg + (k - 1) * jdimg + (j - 1);
                x(j, k) = w(lx + l);
                y(j, k) = w(ly + l);
                z(j, k) = w(lz + l);
            }
        }
    }

    //
    // patch surface is a j=constant surface
    //
    if (jindex / 10 == 2) {
        if (jindex == 21) {
            j = 1;
        } else {
            j = jdimg;
        }
        for (i = 1; i <= idimg; i++) {
            for (k = 1; k <= kdimg; k++) {
                l      = (i - 1) * jdimg * kdimg + (k - 1) * jdimg + (j - 1);
                x(k, i) = w(lx + l);
                y(k, i) = w(ly + l);
                z(k, i) = w(lz + l);
            }
        }
    }

    //
    // patch surface is a k=constant surface
    //
    if (jindex / 10 == 3) {
        if (jindex == 31) {
            k = 1;
        } else {
            k = kdimg;
        }
        for (i = 1; i <= idimg; i++) {
            for (j = 1; j <= jdimg; j++) {
                l      = (i - 1) * jdimg * kdimg + (k - 1) * jdimg + (j - 1);
                x(j, i) = w(lx + l);
                y(j, i) = w(ly + l);
                z(j, i) = w(lz + l);
            }
        }
    }

    return;
}

} // namespace loadgr_ns
