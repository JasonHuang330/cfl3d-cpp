// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "setdqc0.h"

namespace setdqc0_ns {

void setdqc0(int& jdim, int& kdim, int& idim,
             FortranArray4DRef<double> q,
             FortranArray4DRef<double> qc0,
             FortranArray4DRef<double> dqc0)
{
    // COMMON /fluid/
    float& gm1 = cmn_fluid.gm1;

    int idim1 = idim - 1;
    int nt    = jdim * kdim;
    int nplq  = std::min(idim1, 999000 / nt);
    int npl   = nplq;

    for (int i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        int n = nt * npl - jdim - 1;

        // if rho at 1st point is zero, then we know that we do not have qc0
        // data from the restart; no choice but to set dqc0=0 for 1st iteration
        if ((float)qc0(1, 1, 1, 1) == 0.0f) {
            for (int l = 1; l <= 5; l++) {
                for (int izz = 1; izz <= n; izz++) {
                    dqc0(izz, 1, i, l) = 0.0;
                }
            }
        } else {
            for (int izz = 1; izz <= n; izz++) {
                // construct and store conserved variables
                dqc0(izz, 1, i, 1) = q(izz, 1, i, 1) - qc0(izz, 1, i, 1);
                dqc0(izz, 1, i, 2) = q(izz, 1, i, 1) * q(izz, 1, i, 2) - qc0(izz, 1, i, 2);
                dqc0(izz, 1, i, 3) = q(izz, 1, i, 1) * q(izz, 1, i, 3) - qc0(izz, 1, i, 3);
                dqc0(izz, 1, i, 4) = q(izz, 1, i, 1) * q(izz, 1, i, 4) - qc0(izz, 1, i, 4);
                dqc0(izz, 1, i, 5) = q(izz, 1, i, 5) / (double)gm1
                    + 0.5 * q(izz, 1, i, 1)
                    * (q(izz, 1, i, 2) * q(izz, 1, i, 2)
                     + q(izz, 1, i, 3) * q(izz, 1, i, 3)
                     + q(izz, 1, i, 4) * q(izz, 1, i, 4))
                    - qc0(izz, 1, i, 5);
            }
        }
    }
}

} // namespace setdqc0_ns
