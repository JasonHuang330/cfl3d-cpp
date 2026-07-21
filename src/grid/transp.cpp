// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "transp.h"
#include "trans.h"

namespace transp_ns {

void trans(int& jdim, int& kdim, int& idim, FortranArray2DRef<double> t, FortranArray3DRef<double> tti, FortranArray3DRef<double> ttj, FortranArray3DRef<double> ttk, FortranArray1DRef<double> x, FortranArray1DRef<double> y, FortranArray1DRef<double> z, int& itran, double& rfreqt, double& utran, double& vtran, double& wtran, double& xorg, double& yorg, double& zorg, double& xold, double& yold, double& zold, double& xorg0, double& yorg0, double& zorg0, int& iupdat, double& time2)
{
    trans_ns::trans(jdim, kdim, idim, t, tti, ttj, ttk, x, y, z, itran, rfreqt,
                    utran, vtran, wtran, xorg, yorg, zorg, xold, yold, zold,
                    xorg0, yorg0, zorg0, iupdat, time2);
}

void transp(int& mdim, int& ndim, int& jmax1, int& kmax1, int& msub1, int& l, FortranArray3DRef<double> x1, FortranArray3DRef<double> y1, FortranArray3DRef<double> z1, FortranArray2DRef<double> dx, FortranArray2DRef<double> dy, FortranArray2DRef<double> dz, int& intmx, int& int_)
{
    for (int j = 1; j <= jmax1; j++) {
        for (int k = 1; k <= kmax1; k++) {
            x1(j, k, l) = x1(j, k, l) + dx(int_, l);
            y1(j, k, l) = y1(j, k, l) + dy(int_, l);
            z1(j, k, l) = z1(j, k, l) + dz(int_, l);
        }
    }
}

} // namespace transp_ns
