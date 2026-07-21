// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "trans.h"
#include <cmath>

namespace trans_ns {

void trans(int& jdim, int& kdim, int& idim, FortranArray2DRef<double> t,
           FortranArray3DRef<double> tti, FortranArray3DRef<double> ttj,
           FortranArray3DRef<double> ttk, FortranArray1DRef<double> x,
           FortranArray1DRef<double> y, FortranArray1DRef<double> z,
           int& itran, double& rfreqt, double& utran, double& vtran,
           double& wtran, double& xorg, double& yorg, double& zorg,
           double& xold, double& yold, double& zold,
           double& xorg0, double& yorg0, double& zorg0,
           int& iupdat, double& time2)
{
    int idim1, jdim1, kdim1;
    double ft, dfdt, d2fdt2;
    double expt;
    double xnew, ynew, znew;
    double dx, dy, dz;
    double dxdt, dydt, dzdt;
    double d2xdt2, d2ydt2, d2zdt2;
    int n, i, j, k, ii, jj, kk, izz, js;

    idim1 = idim - 1;
    jdim1 = jdim - 1;
    kdim1 = kdim - 1;

    //
    //     ft modulates the displacement
    //     dfdt is the time derivative of ft
    //     d2fdt2 is the second time derivative of ft
    //
    if (itran == 0) {
        return;
    } else if (itran == 1) {
        ft     = time2;
        dfdt   = 1.0;
        d2fdt2 = 0.;
    } else if (itran == 2) {
        ft     = std::sin(rfreqt * time2);
        dfdt   = rfreqt * std::cos(rfreqt * time2);
        d2fdt2 = -(rfreqt) * (rfreqt) * std::sin(rfreqt * time2);
    } else if (itran == 3) {
        expt   = std::exp(-rfreqt * time2);
        ft     = 1. - expt;
        dfdt   = rfreqt * expt;
        d2fdt2 = -(rfreqt) * (rfreqt) * expt;
    } else if (itran == 99) {
        ft     = 0.;
        dfdt   = 1.;
        d2fdt2 = 0.;
    }

    //
    //     xnew = utran*ft + xorg0
    //     ynew = vtran*ft + yorg0
    //     znew = wtran*ft + zorg0
    //
    if (itran != 99) {
        xnew = utran * ft + xorg0;
        ynew = vtran * ft + yorg0;
        znew = wtran * ft + zorg0;
    } else {
        xnew = xorg;
        ynew = yorg;
        znew = zorg;
    }

    dx     = xnew - xold;
    dy     = ynew - yold;
    dz     = znew - zold;
    dxdt   = utran * dfdt;
    dydt   = vtran * dfdt;
    dzdt   = wtran * dfdt;
    d2xdt2 = utran * d2fdt2;
    d2ydt2 = vtran * d2fdt2;
    d2zdt2 = wtran * d2fdt2;

    //
    //     calculate increment to speed of grid points
    //     due to translation and add to current values
    //     t(1)=dx/dt t(2)=dy/dt t(3)=dz/dt
    //
    n = jdim * kdim;
    for (i = 1; i <= idim; i++) {
        js = jdim * kdim * (i - 1) + 1;
        for (izz = 1; izz <= n; izz++) {
            t(izz + js - 1, 1) = t(izz + js - 1, 1) + dxdt;
            t(izz + js - 1, 2) = t(izz + js - 1, 2) + dydt;
            t(izz + js - 1, 3) = t(izz + js - 1, 3) + dzdt;
        }
    }

    //
    //     calculate increment to acceleration of grid
    //     points on the block boundaries due to trans-
    //     lation and add to current values
    //     tti(1)=d2x/dt2 tti(2)=d2y/dt2 tti(3)=d2z/dt2
    //
    //     i0/idim boundaries
    //
    n = jdim * kdim;
    i = 1;
    for (ii = 1; ii <= 2; ii++) {
        for (izz = 1; izz <= n; izz++) {
            // neglect acceleration for now
            // tti(izz,1,ii) = tti(izz,1,ii) + d2xdt2
            // tti(izz,2,ii) = tti(izz,2,ii) + d2ydt2
            // tti(izz,3,ii) = tti(izz,3,ii) + d2zdt2
            tti(izz, 1, ii) = tti(izz, 1, ii) + 0.0;
            tti(izz, 2, ii) = tti(izz, 2, ii) + 0.0;
            tti(izz, 3, ii) = tti(izz, 3, ii) + 0.0;
        }
        i = i + idim - 1;
    }

    //
    //     j0/jdim boundaries
    //
    n = idim * kdim;
    j = 1;
    for (jj = 1; jj <= 2; jj++) {
        for (izz = 1; izz <= n; izz++) {
            // neglect acceleration for now
            // ttj(izz,1,jj) = ttj(izz,1,jj) + d2xdt2
            // ttj(izz,2,jj) = ttj(izz,2,jj) + d2ydt2
            // ttj(izz,3,jj) = ttj(izz,3,jj) + d2zdt2
            ttj(izz, 1, jj) = ttj(izz, 1, jj) + 0.;
            ttj(izz, 2, jj) = ttj(izz, 2, jj) + 0.;
            ttj(izz, 3, jj) = ttj(izz, 3, jj) + 0.;
        }
        j = j + jdim - 1;
    }

    //
    //     k0/kdim boundaries
    //
    n = jdim * idim;
    k = 1;
    for (kk = 1; kk <= 2; kk++) {
        for (izz = 1; izz <= n; izz++) {
            // neglect acceleration for now
            // ttk(izz,1,kk) = ttk(izz,1,kk) + d2xdt2
            // ttk(izz,2,kk) = ttk(izz,2,kk) + d2ydt2
            // ttk(izz,3,kk) = ttk(izz,3,kk) + d2zdt2
            ttk(izz, 1, kk) = ttk(izz, 1, kk) + 0.;
            ttk(izz, 2, kk) = ttk(izz, 2, kk) + 0.;
            ttk(izz, 3, kk) = ttk(izz, 3, kk) + 0.;
        }
        k = k + kdim - 1;
    }

    //
    //     update grid to new position; update rotation
    //     center to new position
    //
    if (iupdat > 0) {
        n = jdim * kdim * idim;
        for (izz = 1; izz <= n; izz++) {
            x(izz) = x(izz) + dx;
            y(izz) = y(izz) + dy;
            z(izz) = z(izz) + dz;
        }

        xorg = xold + dx;
        yorg = yold + dy;
        zorg = zold + dz;
    }

    return;
}

} // namespace trans_ns
