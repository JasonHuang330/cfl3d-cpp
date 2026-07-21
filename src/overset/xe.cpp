// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "xe.h"
#include "direct.h"
#include "project.h"
#include "xe2.h"

namespace xe_ns {

void xe(int& jdim, int& kdim, int& nsub, int& l,
        FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
        FortranArray3DRef<double> xmid, FortranArray3DRef<double> ymid, FortranArray3DRef<double> zmid,
        FortranArray3DRef<double> xmide, FortranArray3DRef<double> ymide, FortranArray3DRef<double> zmide,
        int& jcell, int& kcell,
        double& xc, double& yc, double& zc,
        double& xie, double& eta,
        int& imiss, int& ifit, int& ic0,
        FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
        int& nbuf, int& ibufdim, int& myid)
{
    int idum1, idum2, idum3, idum4;
    double dum1, dum2, dum3;
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
    double x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8;
    double a1, a2, a3;
    double a[4]; // 1-based: a(1),a(2),a(3)
    int itoss;
    double tol, prod;

    idum1 = 0;
    idum2 = 0;
    idum3 = 0;
    idum4 = 0;
    dum1  = 0.;
    dum2  = 0.;
    dum3  = 0.;

    x1 = x(jcell,   kcell,   l);
    y1 = y(jcell,   kcell,   l);
    z1 = z(jcell,   kcell,   l);
    x2 = x(jcell+1, kcell,   l);
    y2 = y(jcell+1, kcell,   l);
    z2 = z(jcell+1, kcell,   l);
    x4 = x(jcell,   kcell+1, l);
    y4 = y(jcell,   kcell+1, l);
    z4 = z(jcell,   kcell+1, l);
    x3 = x(jcell+1, kcell+1, l);
    y3 = y(jcell+1, kcell+1, l);
    z3 = z(jcell+1, kcell+1, l);
    x5 = xmid(jcell,   kcell,   l);
    y5 = ymid(jcell,   kcell,   l);
    z5 = zmid(jcell,   kcell,   l);
    x6 = xmid(jcell,   kcell+1, l);
    y6 = ymid(jcell,   kcell+1, l);
    z6 = zmid(jcell,   kcell+1, l);
    x7 = xmide(jcell,   kcell,   l);
    y7 = ymide(jcell,   kcell,   l);
    z7 = zmide(jcell,   kcell,   l);
    x8 = xmide(jcell+1, kcell,   l);
    y8 = ymide(jcell+1, kcell,   l);
    z8 = zmide(jcell+1, kcell,   l);

    // compute normalized directed areas/unit normals of "from" cell
    // (corresponding values for "to" cell are stored in common/areas/)
    direct_ns::direct(x5, x6, x7, x8, y5, y6, y7, y8, z5, z6, z7, z8,
                      a1, a2, a3, itoss, nou, bou, nbuf, ibufdim);
    a[1] = a1;
    a[2] = a2;
    a[3] = a3;

    // project current "from" cell node points onto plane defined by "to" cell
    //
    // check inner product of unit normals for compatable orientation of "to"
    // and "from" cells before projection.  If not compatable, exit to try
    // a nearby "from" cell

    tol  = .1;
    prod = (double)cmn_areas.ap[0]*a[1] + (double)cmn_areas.ap[1]*a[2] + (double)cmn_areas.ap[2]*a[3];

    if (std::abs((float)prod) > (float)tol) {
        if (ic0 == 0) {
            // ap(1), ap(2), ap(3) are float in COMMON — cast to double for project
            double ap1 = (double)cmn_areas.ap[0];
            double ap2 = (double)cmn_areas.ap[1];
            double ap3 = (double)cmn_areas.ap[2];
            project_ns::project(xc, yc, zc, x1, y1, z1, ap1, ap2, ap3, x1, y1, z1);
            project_ns::project(xc, yc, zc, x2, y2, z2, ap1, ap2, ap3, x2, y2, z2);
            project_ns::project(xc, yc, zc, x3, y3, z3, ap1, ap2, ap3, x3, y3, z3);
            project_ns::project(xc, yc, zc, x4, y4, z4, ap1, ap2, ap3, x4, y4, z4);
            project_ns::project(xc, yc, zc, x5, y5, z5, ap1, ap2, ap3, x5, y5, z5);
            project_ns::project(xc, yc, zc, x6, y6, z6, ap1, ap2, ap3, x6, y6, z6);
            project_ns::project(xc, yc, zc, x7, y7, z7, ap1, ap2, ap3, x7, y7, z7);
            project_ns::project(xc, yc, zc, x8, y8, z8, ap1, ap2, ap3, x8, y8, z8);
        }
    } else {
        // call trace(8,idum1,idum2,idum3,idum4,dum1,dum2,dum3)
        imiss = 1;
        xie   = xie + 1.;
        eta   = eta + 1.;
        return;
    }

    // use best coordinates for inversion in "to" cell... toss out
    // equation for which directed area ap(i) is maximum.
    itoss = cmn_areas.imaxa;

    if (itoss == 1) {
        // use only y and z equations
        // call trace(9,l,jcell,kcell,idum4,dum1,dum2,dum3)
        xe2_ns::xe2(y1, y2, y3, y4, y5, y6, y7, y8, yc,
                    z1, z2, z3, z4, z5, z6, z7, z8, zc,
                    xie, eta, imiss, ifit, nou, bou, nbuf, ibufdim, myid);
    }

    if (itoss == 2) {
        // use only x and z equations
        // call trace(10,l,jcell,kcell,idum4,dum1,dum2,dum3)
        xe2_ns::xe2(x1, x2, x3, x4, x5, x6, x7, x8, xc,
                    z1, z2, z3, z4, z5, z6, z7, z8, zc,
                    xie, eta, imiss, ifit, nou, bou, nbuf, ibufdim, myid);
    }

    if (itoss == 3) {
        // use only x and y equations
        // call trace(11,l,jcell,kcell,idum4,dum1,dum2,dum3)
        xe2_ns::xe2(x1, x2, x3, x4, x5, x6, x7, x8, xc,
                    y1, y2, y3, y4, y5, y6, y7, y8, yc,
                    xie, eta, imiss, ifit, nou, bou, nbuf, ibufdim, myid);
    }

    return;
}

} // namespace xe_ns
