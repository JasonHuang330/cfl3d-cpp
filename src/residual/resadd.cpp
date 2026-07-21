// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "resadd.h"
#include "hole.h"
#include "runtime/common_blocks.h"
#include "runtime/fortran_array.h"
#include <cmath>
#include <algorithm>

namespace resadd_ns {

void resadd(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qc0, FortranArray4DRef<double> dqc0,
            FortranArray4DRef<double> res, FortranArray3DRef<double> vol,
            int& iover, FortranArray3DRef<double> blank)
{
    // COMMON block references
    float& dt      = cmn_info.dt;
    int32_t& ita   = cmn_unst.ita;
    int32_t& ncyc  = cmn_mgrd.ncyc;
    float& gm1     = cmn_fluid.gm1;

    // imult: multi-plane vectorization flag
    //        = 0 single plane at a time
    //        > 0 multiple planes at a time
    int imult = 1;

    int idim1 = idim - 1;
    int nt    = jdim * kdim;
    int nplq  = std::min(idim1, 999000 / nt);
    if (imult == 0) nplq = 1;
    int npl   = nplq;

    double tfact;
    if (std::abs(ita) == 1) {
        tfact = 0.e0;
    } else {
        tfact = 0.5e0 / (double)dt;
    }
    double tfacp1 = tfact + 1.e0 / (double)dt;

    // Main loop over i planes
    for (int i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        int n = nt * npl - jdim - 1;

        // ita=2 -- second order temporal differencing
        if (std::abs(ita) == 2) {
            // cdir$ ivdep
            for (int izz = 1; izz <= n; izz++) {
                res(izz, 1, i, 1) = res(izz, 1, i, 1)
                    - tfact * vol(izz, 1, i) * dqc0(izz, 1, i, 1);

                res(izz, 1, i, 2) = res(izz, 1, i, 2)
                    - tfact * vol(izz, 1, i) * dqc0(izz, 1, i, 2);

                res(izz, 1, i, 3) = res(izz, 1, i, 3)
                    - tfact * vol(izz, 1, i) * dqc0(izz, 1, i, 3);

                res(izz, 1, i, 4) = res(izz, 1, i, 4)
                    - tfact * vol(izz, 1, i) * dqc0(izz, 1, i, 4);

                res(izz, 1, i, 5) = res(izz, 1, i, 5)
                    - tfact * vol(izz, 1, i) * dqc0(izz, 1, i, 5);
            }
        }

        // ncyc > 1 -- temporal subiteration
        if (ncyc > 1) {
            // cdir$ ivdep
            for (int izz = 1; izz <= n; izz++) {
                res(izz, 1, i, 1) = res(izz, 1, i, 1)
                    + tfacp1 * vol(izz, 1, i)
                    * (q(izz, 1, i, 1) - qc0(izz, 1, i, 1));

                res(izz, 1, i, 2) = res(izz, 1, i, 2)
                    + tfacp1 * vol(izz, 1, i)
                    * (q(izz, 1, i, 1) * q(izz, 1, i, 2) - qc0(izz, 1, i, 2));

                res(izz, 1, i, 3) = res(izz, 1, i, 3)
                    + tfacp1 * vol(izz, 1, i)
                    * (q(izz, 1, i, 1) * q(izz, 1, i, 3) - qc0(izz, 1, i, 3));

                res(izz, 1, i, 4) = res(izz, 1, i, 4)
                    + tfacp1 * vol(izz, 1, i)
                    * (q(izz, 1, i, 1) * q(izz, 1, i, 4) - qc0(izz, 1, i, 4));

                res(izz, 1, i, 5) = res(izz, 1, i, 5)
                    + tfacp1 * vol(izz, 1, i)
                    * (q(izz, 1, i, 5) / (double)gm1
                       + 0.5 * q(izz, 1, i, 1)
                         * (q(izz, 1, i, 2) * q(izz, 1, i, 2)
                            + q(izz, 1, i, 3) * q(izz, 1, i, 3)
                            + q(izz, 1, i, 4) * q(izz, 1, i, 4))
                       - qc0(izz, 1, i, 5));
            }
        }
    } // end do 50

    // zero out extra layers j=jdim and k=kdim for safety
    npl = nplq;
    for (int i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        for (int l = 1; l <= 5; l++) {
            for (int ipl = 1; ipl <= npl; ipl++) {
                int ii = i + ipl - 1;
                for (int j = 1; j <= jdim; j++) {
                    res(j, kdim, ii, l) = 0.;
                }
                // cdir$ ivdep
                for (int k = 1; k <= kdim - 1; k++) {
                    res(jdim, k, ii, l) = 0.;
                }
            }
        }
    } // end do 500

    // zero out rhs for fringe and hole cells
    if (iover == 1) {
        npl = nplq;
        for (int i = 1; i <= idim1; i += nplq) {
            if (i + npl - 1 > idim1) npl = idim1 - i + 1;
            hole_ns::hole(i, npl, jdim, kdim, idim, res, blank);
        }
    }

    return;
}

} // namespace resadd_ns
