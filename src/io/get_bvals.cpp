// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "get_bvals.h"

namespace get_bvals_ns {

void get_bvals(int& jdim1, int& kdim1, int& idim1, int& jmax1, int& kmax1, int& ldim, FortranArray4DRef<double> q, FortranArray3DRef<double> qbou, int& mtype, int& i1)
{
    if (mtype == 1) {

        for (int j = 1; j <= jmax1; j++) {
            for (int k = 1; k <= kmax1; k++) {
                for (int l = 1; l <= ldim; l++) {
                    qbou(j,k,l) = q(j,k,i1,l);
                }
            }
        }

    } else if (mtype == 2) {

        for (int j = 1; j <= jmax1; j++) {
            for (int k = 1; k <= kmax1; k++) {
                for (int l = 1; l <= ldim; l++) {
                    qbou(j,k,l) = q(i1,j,k,l);
                }
            }
        }

    } else if (mtype == 3) {

        for (int j = 1; j <= jmax1; j++) {
            for (int k = 1; k <= kmax1; k++) {
                for (int l = 1; l <= ldim; l++) {
                    qbou(j,k,l) = q(j,i1,k,l);
                }
            }
        }

    }

    return;
}

} // namespace get_bvals_ns
