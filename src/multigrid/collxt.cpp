// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "collxt.h"
#include "collx.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <algorithm>

namespace collxt_ns {

void collx(FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> xx, FortranArray3DRef<double> yy, FortranArray3DRef<double> zz, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2)
{
    collx_ns::collx(x, y, z, xx, yy, zz, jdim, kdim, idim, jj2, kk2, ii2);
}

void collxt(FortranArray4DRef<double> xt, FortranArray4DRef<double> xtt, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, int& nbl, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    // common /sklton/ isklton
    int& isklton = cmn_sklton.isklton;

    int nbl1 = nbl + 1;

    if (isklton > 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "     restricting grid speeds from finer block%4d to coarser block%4d",
            nbl, nbl1);
    }

    // restrict xt to coarser mesh
    // jdim,kdim,idim  finer mesh
    // jj2,kk2,ii2     coarser mesh
    for (int ll = 1; ll <= 3; ll++) {
        int ii = 0;
        int iinc = 2;
        if (idim == 2) iinc = 1;
        for (int i = 1; i <= idim; i += iinc) {
            ii = ii + 1;
            int kk = 0;
            for (int k = 1; k <= kdim; k += 2) {
                kk = kk + 1;
                int jj = 0;
                for (int j = 1; j <= jdim; j += 2) {
                    jj = jj + 1;
                    xtt(jj, kk, ii, ll) = xt(j, k, i, ll);
                }
            }
        }
    }
}

} // namespace collxt_ns
