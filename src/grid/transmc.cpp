// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "transmc.h"
#include "trans.h"
#include <cmath>

namespace transmc_ns {

void trans(int& jdim, int& kdim, int& idim, FortranArray2DRef<double> t, FortranArray3DRef<double> tti, FortranArray3DRef<double> ttj, FortranArray3DRef<double> ttk, FortranArray1DRef<double> x, FortranArray1DRef<double> y, FortranArray1DRef<double> z, int& itran, double& rfreqt, double& utran, double& vtran, double& wtran, double& xorg, double& yorg, double& zorg, double& xold, double& yold, double& zold, double& xorg0, double& yorg0, double& zorg0, int& iupdat, double& time2)
{
    trans_ns::trans(jdim, kdim, idim, t, tti, ttj, ttk, x, y, z, itran, rfreqt, utran, vtran, wtran, xorg, yorg, zorg, xold, yold, zold, xorg0, yorg0, zorg0, iupdat, time2);
}

void transmc(int& itran, double& rfreqt, double& utran, double& vtran, double& wtran, double& xorg, double& yorg, double& zorg, double& xorg0, double& yorg0, double& zorg0, double& xmc, double& ymc, double& zmc, int& iupdat, double& time2)
{
    double ft, expt;
    double xold, yold, zold;
    double xnew, ynew, znew;
    double dx, dy, dz;

    //
    //     ft modulates the displacement
    //
    if (itran == 0) {
        return;
    } else if (itran == 1) {
        ft = time2;
    } else if (itran == 2) {
        ft = std::sin(rfreqt * time2);
    } else if (itran == 3) {
        expt = std::exp(-rfreqt * time2);
        ft   = 1. - expt;
    }

    xold = xorg;
    yold = yorg;
    zold = zorg;

    xnew = utran * ft + xorg0;
    ynew = vtran * ft + yorg0;
    znew = wtran * ft + zorg0;

    dx = xnew - xold;
    dy = ynew - yold;
    dz = znew - zold;

    //*************************************************
    //     update moment center to new position;
    //     update rotation point for moment center
    //*************************************************

    if (iupdat > 0) {

        xmc = xmc + dx;
        ymc = ymc + dy;
        zmc = zmc + dz;

        xorg = xorg + dx;
        yorg = yorg + dy;
        zorg = zorg + dz;

    }

    return;
}

} // namespace transmc_ns
