// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "cblkk_d.h"
#include "cblkk.h"
#include <cmath>
#include <algorithm>

namespace cblkk_d_ns {

// cblkk: wrapper that delegates to cblkk_ns::cblkk
void cblkk(int& nbli, int& idimr, int& jdimr, int& kdimr, int& idimt, int& jdimt, int& kdimt,
           FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
           int& it, int& ir, int& iedge,
           FortranArray3DRef<double> xr, FortranArray3DRef<double> yr, FortranArray3DRef<double> zr,
           FortranArray3DRef<double> xt, FortranArray3DRef<double> yt, FortranArray3DRef<double> zt,
           int& ntime, int& lcnt, FortranArray1DRef<double> geom_miss, int& mxbli)
{
    cblkk_ns::cblkk(nbli, idimr, jdimr, kdimr, idimt, jdimt, kdimt,
                    limblk, isva, it, ir, iedge,
                    xr, yr, zr, xt, yt, zt,
                    ntime, lcnt, geom_miss, mxbli);
}

// cblkk_d: Check information transferred from block (ir) to qk0 array of block (it).
void cblkk_d(int& nbli, int& idimr, int& jdimr, int& kdimr, int& idimt, int& jdimt, int& kdimt,
             FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
             int& it, int& ir, int& iedge,
             FortranArray2DRef<double> xyzr, int& nvals,
             FortranArray3DRef<double> xt, FortranArray3DRef<double> yt, FortranArray3DRef<double> zt,
             int& ntime, int& lcnt, FortranArray1DRef<double> geom_miss, int& mxbli)
{
    int ist, iet, iinct;
    int jst, jet, jinct;
    int kt;
    int ij, i, j;
    double xr1, yr1, zr1;
    double xt1, yt1, zt1;
    double eps;

    ist = limblk(it, 1);
    iet = limblk(it, 4);
    if (ist == iet) {
        iinct = 1;
    } else {
        iinct = (iet - ist) / std::abs(iet - ist);
    }

    jst = limblk(it, 2);
    jet = limblk(it, 5);
    if (jst == jet) {
        jinct = 1;
    } else {
        jinct = (jet - jst) / std::abs(jet - jst);
    }

    eps = 0.;
    kt  = 1;
    if (iedge == 2) kt = kdimt;

    // determine the side of the q array to transfer from
    // k = constant side
    if (isva(ir, 1) + isva(ir, 2) == 3) {
        if ((isva(ir, 1) == isva(it, 1)) ||
            (isva(ir, 2) == isva(it, 2))) {
            // i varies with i     and     j varies with j
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(j,   kt, i) + xt(j,   kt, i+1)
                                + xt(j+1, kt, i) + xt(j+1, kt, i+1) );
                    yt1 = .25 * ( yt(j,   kt, i) + yt(j,   kt, i+1)
                                + yt(j+1, kt, i) + yt(j+1, kt, i+1) );
                    zt1 = .25 * ( zt(j,   kt, i) + zt(j,   kt, i+1)
                                + zt(j+1, kt, i) + zt(j+1, kt, i+1) );
                    eps = std::max(eps, std::sqrt((xr1-xt1)*(xr1-xt1)
                                               + (yr1-yt1)*(yr1-yt1)
                                               + (zr1-zt1)*(zr1-zt1)));
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }
        } else {
            // j varies with i     and     i varies with j
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(j,   kt, i) + xt(j,   kt, i+1)
                                + xt(j+1, kt, i) + xt(j+1, kt, i+1) );
                    yt1 = .25 * ( yt(j,   kt, i) + yt(j,   kt, i+1)
                                + yt(j+1, kt, i) + yt(j+1, kt, i+1) );
                    zt1 = .25 * ( zt(j,   kt, i) + zt(j,   kt, i+1)
                                + zt(j+1, kt, i) + zt(j+1, kt, i+1) );
                    eps = std::max(eps, std::sqrt((xr1-xt1)*(xr1-xt1)
                                               + (yr1-yt1)*(yr1-yt1)
                                               + (zr1-zt1)*(zr1-zt1)));
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }
        }

    // j = constant side
    } else if (isva(ir, 1) + isva(ir, 2) == 4) {
        if ((isva(ir, 1) == isva(it, 1)) ||
            (isva(ir, 2) == isva(it, 2))) {
            // i varies with i    and    k varies with j
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(j,   kt, i) + xt(j,   kt, i+1)
                                + xt(j+1, kt, i) + xt(j+1, kt, i+1) );
                    yt1 = .25 * ( yt(j,   kt, i) + yt(j,   kt, i+1)
                                + yt(j+1, kt, i) + yt(j+1, kt, i+1) );
                    zt1 = .25 * ( zt(j,   kt, i) + zt(j,   kt, i+1)
                                + zt(j+1, kt, i) + zt(j+1, kt, i+1) );
                    eps = std::max(eps, std::sqrt((xr1-xt1)*(xr1-xt1)
                                               + (yr1-yt1)*(yr1-yt1)
                                               + (zr1-zt1)*(zr1-zt1)));
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }
        } else {
            // k varies with i    and    i varies with j
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(j,   kt, i) + xt(j,   kt, i+1)
                                + xt(j+1, kt, i) + xt(j+1, kt, i+1) );
                    yt1 = .25 * ( yt(j,   kt, i) + yt(j,   kt, i+1)
                                + yt(j+1, kt, i) + yt(j+1, kt, i+1) );
                    zt1 = .25 * ( zt(j,   kt, i) + zt(j,   kt, i+1)
                                + zt(j+1, kt, i) + zt(j+1, kt, i+1) );
                    eps = std::max(eps, std::sqrt((xr1-xt1)*(xr1-xt1)
                                               + (yr1-yt1)*(yr1-yt1)
                                               + (zr1-zt1)*(zr1-zt1)));
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }
        }

    // i = constant side
    } else if (isva(ir, 1) + isva(ir, 2) == 5) {
        if ((isva(ir, 1) == isva(it, 1)) ||
            (isva(ir, 2) == isva(it, 2))) {
            // k varies with i    and    j varies with j
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(j,   kt, i) + xt(j,   kt, i+1)
                                + xt(j+1, kt, i) + xt(j+1, kt, i+1) );
                    yt1 = .25 * ( yt(j,   kt, i) + yt(j,   kt, i+1)
                                + yt(j+1, kt, i) + yt(j+1, kt, i+1) );
                    zt1 = .25 * ( zt(j,   kt, i) + zt(j,   kt, i+1)
                                + zt(j+1, kt, i) + zt(j+1, kt, i+1) );
                    eps = std::max(eps, std::sqrt((xr1-xt1)*(xr1-xt1)
                                               + (yr1-yt1)*(yr1-yt1)
                                               + (zr1-zt1)*(zr1-zt1)));
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }
        } else {
            // j varies with i    and    k varies with j
            ij = 0;
            for (i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (j = jst; jinct > 0 ? j <= jet : j >= jet; j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(j,   kt, i) + xt(j,   kt, i+1)
                                + xt(j+1, kt, i) + xt(j+1, kt, i+1) );
                    yt1 = .25 * ( yt(j,   kt, i) + yt(j,   kt, i+1)
                                + yt(j+1, kt, i) + yt(j+1, kt, i+1) );
                    zt1 = .25 * ( zt(j,   kt, i) + zt(j,   kt, i+1)
                                + zt(j+1, kt, i) + zt(j+1, kt, i+1) );
                    eps = std::max(eps, std::sqrt((xr1-xt1)*(xr1-xt1)
                                               + (yr1-yt1)*(yr1-yt1)
                                               + (zr1-zt1)*(zr1-zt1)));
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }
        }
    }

    return;
}

} // namespace cblkk_d_ns
