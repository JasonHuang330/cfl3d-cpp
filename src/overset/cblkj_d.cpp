// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "cblkj_d.h"
#include "cblkj.h"
#include <cmath>
#include <algorithm>

namespace cblkj_d_ns {

static inline double ccmax(double a, double b) { return std::max(a, b); }

void cblkj(int& nbli, int& idimr, int& jdimr, int& kdimr, int& idimt, int& jdimt, int& kdimt, FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva, int& it, int& ir, int& iedge, FortranArray3DRef<double> xr, FortranArray3DRef<double> yr, FortranArray3DRef<double> zr, FortranArray3DRef<double> xt, FortranArray3DRef<double> yt, FortranArray3DRef<double> zt, int& ntime, int& lcnt, FortranArray1DRef<double> geom_miss, int& mxbli)
{
    cblkj_ns::cblkj(nbli, idimr, jdimr, kdimr, idimt, jdimt, kdimt,
                    limblk, isva, it, ir, iedge,
                    xr, yr, zr, xt, yt, zt,
                    ntime, lcnt, geom_miss, mxbli);
}

void cblkj_d(int& nbli, int& idimt, int& jdimt, int& kdimt, FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva, int& it, int& ir, int& iedge, FortranArray2DRef<double> xyzr, int& nvals, FortranArray3DRef<double> xt, FortranArray3DRef<double> yt, FortranArray3DRef<double> zt, int& ntime, int& lcnt, FortranArray1DRef<double> geom_miss, int& mxbli)
{
    int ist, iet, iinct;
    int kst, ket, kinct;
    int jt, ij;
    double eps;
    double xr1, yr1, zr1;
    double xt1, yt1, zt1;

    ist = limblk(it, 1);
    iet = limblk(it, 4);
    if (ist == iet) {
        iinct = 1;
    } else {
        iinct = (iet - ist) / std::abs(iet - ist);
    }

    kst = limblk(it, 3);
    ket = limblk(it, 6);
    if (kst == ket) {
        kinct = 1;
    } else {
        kinct = (ket - kst) / std::abs(ket - kst);
    }

    eps = 0.;
    jt  = 1;
    if (iedge == 2) jt = jdimt;

    //
    // determine the side of the q array to transfer from
    //

    //
    // k = constant side
    //
    if (isva(ir, 1) + isva(ir, 2) == 3) {
        if ((isva(ir, 1) == isva(it, 1)) ||
            (isva(ir, 2) == isva(it, 2))) {
            //
            // i varies with i     and     j varies with k
            //
            ij = 0;
            for (int i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (int k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(jt, k,   i) + xt(jt, k,   i+1)
                                + xt(jt, k+1, i) + xt(jt, k+1, i+1) );
                    yt1 = .25 * ( yt(jt, k,   i) + yt(jt, k,   i+1)
                                + yt(jt, k+1, i) + yt(jt, k+1, i+1) );
                    zt1 = .25 * ( zt(jt, k,   i) + zt(jt, k,   i+1)
                                + zt(jt, k+1, i) + zt(jt, k+1, i+1) );
                    eps = ccmax(eps, std::sqrt((xr1-xt1)*(xr1-xt1) +
                                              (yr1-yt1)*(yr1-yt1) +
                                              (zr1-zt1)*(zr1-zt1)));
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }

        } else {

            ij = 0;
            for (int i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (int k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(jt, k,   i) + xt(jt, k,   i+1)
                                + xt(jt, k+1, i) + xt(jt, k+1, i+1) );
                    yt1 = .25 * ( yt(jt, k,   i) + yt(jt, k,   i+1)
                                + yt(jt, k+1, i) + yt(jt, k+1, i+1) );
                    zt1 = .25 * ( zt(jt, k,   i) + zt(jt, k,   i+1)
                                + zt(jt, k+1, i) + zt(jt, k+1, i+1) );
                    eps = ccmax(eps, std::sqrt((xr1-xt1)*(xr1-xt1) +
                                              (yr1-yt1)*(yr1-yt1) +
                                              (zr1-zt1)*(zr1-zt1)));
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }

        }

    //
    // j = constant side
    //
    } else if (isva(ir, 1) + isva(ir, 2) == 4) {
        if ((isva(ir, 1) == isva(it, 1)) ||
            (isva(ir, 2) == isva(it, 2))) {
            //
            // i varies with i    and    k varies with k
            //
            ij = 0;
            for (int i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (int k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(jt, k,   i) + xt(jt, k,   i+1)
                                + xt(jt, k+1, i) + xt(jt, k+1, i+1) );
                    yt1 = .25 * ( yt(jt, k,   i) + yt(jt, k,   i+1)
                                + yt(jt, k+1, i) + yt(jt, k+1, i+1) );
                    zt1 = .25 * ( zt(jt, k,   i) + zt(jt, k,   i+1)
                                + zt(jt, k+1, i) + zt(jt, k+1, i+1) );
                    eps = ccmax(eps, std::sqrt((xr1-xt1)*(xr1-xt1) +
                                              (yr1-yt1)*(yr1-yt1) +
                                              (zr1-zt1)*(zr1-zt1)));
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }

        } else {
            //
            // k varies with i    and    i varies with k
            //
            ij = 0;
            for (int i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (int k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(jt, k,   i) + xt(jt, k,   i+1)
                                + xt(jt, k+1, i) + xt(jt, k+1, i+1) );
                    yt1 = .25 * ( yt(jt, k,   i) + yt(jt, k,   i+1)
                                + yt(jt, k+1, i) + yt(jt, k+1, i+1) );
                    zt1 = .25 * ( zt(jt, k,   i) + zt(jt, k,   i+1)
                                + zt(jt, k+1, i) + zt(jt, k+1, i+1) );
                    eps = ccmax(eps, std::sqrt((xr1-xt1)*(xr1-xt1) +
                                              (yr1-yt1)*(yr1-yt1) +
                                              (zr1-zt1)*(zr1-zt1)));
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }

        }

    //
    // i = constant side
    //
    } else if (isva(ir, 1) + isva(ir, 2) == 5) {
        if ((isva(ir, 1) == isva(it, 1)) ||
            (isva(ir, 2) == isva(it, 2))) {
            //
            // k varies with k    and    j varies with i
            //
            ij = 0;
            for (int i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (int k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(jt, k,   i) + xt(jt, k,   i+1)
                                + xt(jt, k+1, i) + xt(jt, k+1, i+1) );
                    yt1 = .25 * ( yt(jt, k,   i) + yt(jt, k,   i+1)
                                + yt(jt, k+1, i) + yt(jt, k+1, i+1) );
                    zt1 = .25 * ( zt(jt, k,   i) + zt(jt, k,   i+1)
                                + zt(jt, k+1, i) + zt(jt, k+1, i+1) );
                    eps = ccmax(eps, std::sqrt((xr1-xt1)*(xr1-xt1) +
                                              (yr1-yt1)*(yr1-yt1) +
                                              (zr1-zt1)*(zr1-zt1)));
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }

        } else {
            //
            // k varies with i    and    j varies with k
            //
            ij = 0;
            for (int i = ist; iinct > 0 ? i <= iet : i >= iet; i += iinct) {
                for (int k = kst; kinct > 0 ? k <= ket : k >= ket; k += kinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * ( xt(jt, k,   i) + xt(jt, k,   i+1)
                                + xt(jt, k+1, i) + xt(jt, k+1, i+1) );
                    yt1 = .25 * ( yt(jt, k,   i) + yt(jt, k,   i+1)
                                + yt(jt, k+1, i) + yt(jt, k+1, i+1) );
                    zt1 = .25 * ( zt(jt, k,   i) + zt(jt, k,   i+1)
                                + zt(jt, k+1, i) + zt(jt, k+1, i+1) );
                    eps = ccmax(eps, std::sqrt((xr1-xt1)*(xr1-xt1) +
                                              (yr1-yt1)*(yr1-yt1) +
                                              (zr1-zt1)*(zr1-zt1)));
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

} // namespace cblkj_d_ns
