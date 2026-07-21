// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "fill.h"

namespace fill_ns {

void fill(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, int& ll)
{
//***********************************************************************
//     Purpose:  Fills the edges of the q array for safety using
//     multi-plane vectorization technique.
//***********************************************************************

    int idim1 = idim - 1;
    int kdim1 = kdim - 1;

    for (int l = 1; l <= ll; l++) {
        for (int i = 1; i <= idim1; i++) {
            // cdir$ ivdep
            for (int k = 1; k <= kdim1; k++) {
                q(jdim, k, i, l) = q(jdim - 1, k, i, l);
            }
            for (int izz = 1; izz <= jdim; izz++) {
                q(izz, kdim, i, l) = q(izz, kdim - 1, i, l);
            }
        }
        int n = jdim * kdim;
        // cdir$ ivdep
        // Fortran flat-index loop: q(izz,1,idim,l) with izz=1..jdim*kdim
        // This iterates over all (j,k) pairs in column-major order
        for (int k = 1; k <= kdim; k++) {
            for (int j = 1; j <= jdim; j++) {
                q(j, k, idim, l) = q(j, k, idim - 1, l);
            }
        }
    }
}

} // namespace fill_ns
