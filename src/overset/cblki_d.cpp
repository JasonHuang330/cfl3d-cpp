// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "cblki_d.h"
#include "cblki.h"
#include "ccomplex.h"
#include <cmath>

namespace cblki_d_ns {

void cblki(int& nbli, int& idimr, int& jdimr, int& kdimr, int& idimt, int& jdimt, int& kdimt,
           FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
           int& it, int& ir, int& iedge,
           FortranArray3DRef<double> xr, FortranArray3DRef<double> yr, FortranArray3DRef<double> zr,
           FortranArray3DRef<double> xt, FortranArray3DRef<double> yt, FortranArray3DRef<double> zt,
           int& ntime, int& lcnt, FortranArray1DRef<double> geom_miss, int& mxbli)
{
    cblki_ns::cblki(nbli, idimr, jdimr, kdimr, idimt, jdimt, kdimt,
                    limblk, isva, it, ir, iedge,
                    xr, yr, zr, xt, yt, zt,
                    ntime, lcnt, geom_miss, mxbli);
}

void cblki_d(int& nbli, int& idimr, int& jdimr, int& kdimr, int& idimt, int& jdimt, int& kdimt,
             FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
             int& it, int& ir, int& iedge,
             FortranArray2DRef<double> xyzr, int& nvals,
             FortranArray3DRef<double> xt, FortranArray3DRef<double> yt, FortranArray3DRef<double> zt,
             int& ntime, int& lcnt, FortranArray1DRef<double> geom_miss, int& mxbli)
{
    int jst, jet, jinct;
    int kst, ket, kinct;
    int itn, ij;
    double eps, xr1, yr1, zr1, xt1, yt1, zt1;

    jst = limblk(it, 2);
    jet = limblk(it, 5);
    if (jst == jet) {
        jinct = 1;
    } else {
        jinct = (jet - jst) / std::abs(jet - jst);
    }

    kst = limblk(it, 3);
    ket = limblk(it, 6);
    if (kst == ket) {
        kinct = 1;
    } else {
        kinct = (ket - kst) / std::abs(ket - kst);
    }

    eps = 0.;
    itn = 1;
    if (iedge == 2) itn = idimt;

    // determine the side of the q array to transfer from

    // k = constant side
    if (isva(ir, 1) + isva(ir, 2) == 3) {
        if ((isva(ir, 1) == isva(it, 1)) || (isva(ir, 2) == isva(it, 2))) {
            // i varies with k     and     j varies with j
            ij = 0;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * (xt(j,   k,   itn) + xt(j,   k+1, itn)
                               + xt(j+1, k,   itn) + xt(j+1, k+1, itn));
                    yt1 = .25 * (yt(j,   k,   itn) + yt(j,   k+1, itn)
                               + yt(j+1, k,   itn) + yt(j+1, k+1, itn));
                    zt1 = .25 * (zt(j,   k,   itn) + zt(j,   k+1, itn)
                               + zt(j+1, k,   itn) + zt(j+1, k+1, itn));
                    double val = std::sqrt((xr1-xt1)*(xr1-xt1) + (yr1-yt1)*(yr1-yt1)
                                         + (zr1-zt1)*(zr1-zt1));
                    eps = ccomplex_ns::ccmax(eps, val);
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }
        } else {
            // j varies with k     and     i varies with j
            ij = 0;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * (xt(j,   k,   itn) + xt(j,   k+1, itn)
                               + xt(j+1, k,   itn) + xt(j+1, k+1, itn));
                    yt1 = .25 * (yt(j,   k,   itn) + yt(j,   k+1, itn)
                               + yt(j+1, k,   itn) + yt(j+1, k+1, itn));
                    zt1 = .25 * (zt(j,   k,   itn) + zt(j,   k+1, itn)
                               + zt(j+1, k,   itn) + zt(j+1, k+1, itn));
                    double val = std::sqrt((xr1-xt1)*(xr1-xt1) + (yr1-yt1)*(yr1-yt1)
                                         + (zr1-zt1)*(zr1-zt1));
                    eps = ccomplex_ns::ccmax(eps, val);
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
        if ((isva(ir, 1) == isva(it, 1)) || (isva(ir, 2) == isva(it, 2))) {
            // i varies with j    and    k varies with k
            ij = 0;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * (xt(j,   k,   itn) + xt(j,   k+1, itn)
                               + xt(j+1, k,   itn) + xt(j+1, k+1, itn));
                    yt1 = .25 * (yt(j,   k,   itn) + yt(j,   k+1, itn)
                               + yt(j+1, k,   itn) + yt(j+1, k+1, itn));
                    zt1 = .25 * (zt(j,   k,   itn) + zt(j,   k+1, itn)
                               + zt(j+1, k,   itn) + zt(j+1, k+1, itn));
                    double val = std::sqrt((xr1-xt1)*(xr1-xt1) + (yr1-yt1)*(yr1-yt1)
                                         + (zr1-zt1)*(zr1-zt1));
                    eps = ccomplex_ns::ccmax(eps, val);
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }
        } else {
            // i varies with k    and    k varies with j
            ij = 0;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * (xt(j,   k,   itn) + xt(j,   k+1, itn)
                               + xt(j+1, k,   itn) + xt(j+1, k+1, itn));
                    yt1 = .25 * (yt(j,   k,   itn) + yt(j,   k+1, itn)
                               + yt(j+1, k,   itn) + yt(j+1, k+1, itn));
                    zt1 = .25 * (zt(j,   k,   itn) + zt(j,   k+1, itn)
                               + zt(j+1, k,   itn) + zt(j+1, k+1, itn));
                    double val = std::sqrt((xr1-xt1)*(xr1-xt1) + (yr1-yt1)*(yr1-yt1)
                                         + (zr1-zt1)*(zr1-zt1));
                    eps = ccomplex_ns::ccmax(eps, val);
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
        if ((isva(ir, 1) == isva(it, 1)) || (isva(ir, 2) == isva(it, 2))) {
            // k varies with k    and    j varies with j
            ij = 0;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * (xt(j,   k,   itn) + xt(j,   k+1, itn)
                               + xt(j+1, k,   itn) + xt(j+1, k+1, itn));
                    yt1 = .25 * (yt(j,   k,   itn) + yt(j,   k+1, itn)
                               + yt(j+1, k,   itn) + yt(j+1, k+1, itn));
                    zt1 = .25 * (zt(j,   k,   itn) + zt(j,   k+1, itn)
                               + zt(j+1, k,   itn) + zt(j+1, k+1, itn));
                    double val = std::sqrt((xr1-xt1)*(xr1-xt1) + (yr1-yt1)*(yr1-yt1)
                                         + (zr1-zt1)*(zr1-zt1));
                    eps = ccomplex_ns::ccmax(eps, val);
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }
        } else {
            // j varies with k    and    k varies with j
            ij = 0;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    ij = ij + 1;
                    xr1 = xyzr(ij, 1);
                    yr1 = xyzr(ij, 2);
                    zr1 = xyzr(ij, 3);
                    xt1 = .25 * (xt(j,   k,   itn) + xt(j,   k+1, itn)
                               + xt(j+1, k,   itn) + xt(j+1, k+1, itn));
                    yt1 = .25 * (yt(j,   k,   itn) + yt(j,   k+1, itn)
                               + yt(j+1, k,   itn) + yt(j+1, k+1, itn));
                    zt1 = .25 * (zt(j,   k,   itn) + zt(j,   k+1, itn)
                               + zt(j+1, k,   itn) + zt(j+1, k+1, itn));
                    double val = std::sqrt((xr1-xt1)*(xr1-xt1) + (yr1-yt1)*(yr1-yt1)
                                         + (zr1-zt1)*(zr1-zt1));
                    eps = ccomplex_ns::ccmax(eps, val);
                }
            }
            if (ntime == 1) {
                geom_miss(lcnt) = eps;
            } else {
                geom_miss(lcnt) = 0.;
            }
        }
    }
}

} // namespace cblki_d_ns
