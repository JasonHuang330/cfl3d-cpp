// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "collxtb.h"
#include "collx.h"
#include "collxt.h"

namespace collxtb_ns {

void collx(FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
           FortranArray3DRef<double> xx, FortranArray3DRef<double> yy, FortranArray3DRef<double> zz,
           int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2)
{
    collx_ns::collx(x, y, z, xx, yy, zz, jdim, kdim, idim, jj2, kk2, ii2);
}

void collxt(FortranArray4DRef<double> xt, FortranArray4DRef<double> xtt,
            int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, int& nbl,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    collxt_ns::collxt(xt, xtt, jdim, kdim, idim, jj2, kk2, ii2, nbl, nou, bou, nbuf, ibufdim);
}

void collxtb(FortranArray4DRef<double> xtt, FortranArray4DRef<double> xttt,
             int& mdim, int& ndim, int& mm2, int& nn2, int& nbl)
{
    int ninc, minc;
    int ii, ll, n, m;
    int nn, mm;

    ninc = 2;
    if (ndim == 2) ninc = 1;
    minc = 2;
    if (mdim == 2) minc = 1;

    for (ii = 1; ii <= 2; ii++) {
        for (ll = 1; ll <= 3; ll++) {
            nn = 0;
            for (n = 1; n <= ndim; n += ninc) {
                nn = nn + 1;
                mm = 0;
                for (m = 1; m <= mdim; m += minc) {
                    mm = mm + 1;
                    xttt(mm, nn, ll, ii) = xtt(m, n, ll, ii);
                }
            }
        }
    }
}

} // namespace collxtb_ns
