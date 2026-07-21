// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "project.h"

namespace project_ns {

void project(double& x0, double& y0, double& z0, double& xc, double& yc, double& zc, double& a1, double& a2, double& a3, double& xcp, double& ycp, double& zcp)
{
//***********************************************************************
//     Purpose:  Project point xc,yc,zc into plane containing the point
//     x0,y0,z0 and which has the unit normal with components a1,a2,a3
//***********************************************************************

    double a;
    a = (xc - x0) * a1 + (yc - y0) * a2 + (zc - z0) * a3;
    xcp = xc - a * a1;
    ycp = yc - a * a2;
    zcp = zc - a * a3;
    return;
}

} // namespace project_ns
