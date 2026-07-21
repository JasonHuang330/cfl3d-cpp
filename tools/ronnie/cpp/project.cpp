// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// project.cpp — faithful C++ translation of project.F (CFL3D ronnie).
// Project the point (xc,yc,zc) into the plane containing (x0,y0,z0) with unit
// normal (a1,a2,a3); result returned in (xcp,ycp,zcp).
//
// Sequential real build: no complex.  Preserves exact arithmetic structure.
#include "ron_common.h"

void project(double x0,double y0,double z0,double xc,double yc,double zc,
             double a1,double a2,double a3,
             double& xcp,double& ycp,double& zcp)
{
    double a = (xc-x0)*a1 + (yc-y0)*a2 + (zc-z0)*a3;
    xcp = xc - a*a1;
    ycp = yc - a*a2;
    zcp = zc - a*a3;
    return;
}
