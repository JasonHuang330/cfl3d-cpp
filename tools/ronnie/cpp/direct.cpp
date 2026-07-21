// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// direct.cpp — faithful C++ translation of direct.F (CFL3D ronnie).
// Compute (normalized) directed area components, i.e. components of unit
// normal to a cell face, from the 4 midpoint coordinates (5,6,7,8).
//
// Sequential real build: real(x)->x, ccabs->std::fabs.  bou is the
// character*120 diagnostic buffer (ibufdim,nbuf), nou(nbuf) the line counter.
#include "ron_common.h"

void direct(double x5,double x6,double x7,double x8,
            double y5,double y6,double y7,double y8,
            double z5,double z6,double z7,double z8,
            double& a1,double& a2,double& a3,int& imaxa,
            int* nou,FStr* bou,int nbuf,int ibufdim)
{
    // bou(ibufdim,nbuf): column-major, bou(i,j) -> bou[(i-1)+(j-1)*ibufdim]
    // nou(nbuf) -> nou[j-1]
    double xxie = x8-x7;
    double yxie = y8-y7;
    double zxie = z8-z7;
    double xeta = x6-x5;
    double yeta = y6-y5;
    double zeta = z6-z5;
    a1 = yxie*zeta-zxie*yeta;
    a2 = zxie*xeta-xxie*zeta;
    a3 = xxie*yeta-yxie*xeta;
    double d = std::sqrt(a1*a1+a2*a2+a3*a3);
    if (d <= 0.0) {
        nou[4-1] = i_min(nou[4-1]+1,ibufdim);
        bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
            "  WARNINING: a cell with zero area has been detected in subroutine direct");
        nou[4-1] = i_min(nou[4-1]+1,ibufdim);
        bou[(nou[4-1]-1)+(4-1)*ibufdim].assign("   - severe problem");
        d = 1.0;
    }
    a1 = a1/d;
    a2 = a2/d;
    a3 = a3/d;

    // find coordinate direction with maximum area component
    double amax = ccabs(a1);
    imaxa = 1;
    if (std::fabs(a2) > amax) {
        amax  = ccabs(a2);
        imaxa = 2;
    }
    if (std::fabs(a3) > amax) {
        amax  = ccabs(a3);
        imaxa = 3;
    }
    return;
}
