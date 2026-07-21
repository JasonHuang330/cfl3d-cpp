// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dsmin.h"

namespace dsmin_ns {

void dsmin(int& jdim, int& kdim, int& nsub, FortranArray1DRef<int> jjmax, FortranArray1DRef<int> kkmax, int& lmax, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, double& xc, double& yc, double& zc, int& j1, int& k1, int& l1, FortranArray1DRef<int> lout, int& lflag, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2)
{
    double dmin, d1;
    int ls, le, l, j, k, js, je, ks, ke;

    dmin = 1.0e+20;

    ls = 1;
    le = lmax;
    if (lflag < 0) {
        ls = l1;
        le = l1;
    }

    for (l = ls; l <= le; l++) {

        // skip over blocks already tried
        if (ls != le) {
            if (lout(l) > 0) continue;
        }

        // search only over specified range, and only inside non-expanded limits
        js = xif1(l);
        je = xif2(l);
        ks = etf1(l);
        ke = etf2(l);
        je = std::min(je - 2, jjmax(l) - 2);
        ke = std::min(ke - 2, kkmax(l) - 2);
        js = std::max(js + 1, 2);
        ks = std::max(ks + 1, 2);

        for (k = ks; k <= ke; k++) {
            for (j = js; j <= je; j++) {
                d1 = (xc - x(j, k, l)) * (xc - x(j, k, l))
                   + (yc - y(j, k, l)) * (yc - y(j, k, l))
                   + (zc - z(j, k, l)) * (zc - z(j, k, l));
                if ((float)d1 < (float)dmin) {
                    j1   = j;
                    k1   = k;
                    l1   = l;
                    dmin = d1;
                }
            }
        }
    }

    return;
}

} // namespace dsmin_ns
