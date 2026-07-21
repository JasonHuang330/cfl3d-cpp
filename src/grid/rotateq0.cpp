// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rotateq0.h"
#include "rotate.h"
#include "rotateq.h"
#include <cmath>

namespace rotateq0_ns {

void rotate(int& jdim, int& kdim, int& idim, FortranArray2DRef<double> t, FortranArray3DRef<double> tti, FortranArray3DRef<double> ttj, FortranArray3DRef<double> ttk, FortranArray2DRef<double> t1, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& nbl, int& irot, double& rfreqr, double& omegx, double& omegy, double& omegz, double& xorg, double& yorg, double& zorg, double& thetax, double& thetay, double& thetaz, double& thxold, double& thyold, double& thzold, int& iupdat, double& time2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    rotate_ns::rotate(jdim, kdim, idim, t, tti, ttj, ttk, t1, x, y, z, nbl, irot, rfreqr, omegx, omegy, omegz, xorg, yorg, zorg, thetax, thetay, thetaz, thxold, thyold, thzold, iupdat, time2, nou, bou, nbuf, ibufdim);
}

void rotateq(int& jd, int& kd, int& id, FortranArray4DRef<double> q, FortranArray4DRef<double> qrot, int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend, double& dthtx, double& dthty, double& dthtz)
{
    rotateq_ns::rotateq(jd, kd, id, q, qrot, ista, iend, jsta, jend, ksta, kend, dthtx, dthty, dthtz);
}

void rotateq0(int& ld, int& md, FortranArray4DRef<double> q0, FortranArray4DRef<double> q0rot, int& lsta, int& lend, int& msta, int& mend, double& dthtx, double& dthty, double& dthtz)
{
    double ca, sa, temp;

    if (std::abs((float)dthtx) > 0.) {
//
//       rotate q0 about an axis parallel to the x-axis
//
        ca = std::cos(dthtx);
        sa = std::sin(dthtx);
//
        for (int n = 1; n <= 4; n++) {
        for (int l = lsta; l <= lend; l++) {
        for (int m = msta; m <= mend; m++) {
            q0rot(l,m,1,n) =  q0(l,m,1,n);
            q0rot(l,m,2,n) =  q0(l,m,2,n);
            temp           =  q0(l,m,3,n);
            q0rot(l,m,3,n) =  q0(l,m,3,n)*ca - q0(l,m,4,n)*sa;
            q0rot(l,m,4,n) =  temp        *sa + q0(l,m,4,n)*ca;
            q0rot(l,m,5,n) =  q0(l,m,5,n);
        }
        }
        }
//
    } else if (std::abs((float)dthty) > 0.) {
//
//       rotate q0 about an axis parallel to the y-axis
//
        ca = std::cos(dthty);
        sa = std::sin(dthty);
//
        for (int n = 1; n <= 4; n++) {
        for (int l = lsta; l <= lend; l++) {
        for (int m = msta; m <= mend; m++) {
            q0rot(l,m,1,n) =  q0(l,m,1,n);
            temp           =  q0(l,m,2,n);
            q0rot(l,m,2,n) =  q0(l,m,2,n)*ca + q0(l,m,4,n)*sa;
            q0rot(l,m,3,n) =  q0(l,m,3,n);
            q0rot(l,m,4,n) = -temp        *sa + q0(l,m,4,n)*ca;
            q0rot(l,m,5,n) =  q0(l,m,5,n);
        }
        }
        }
//
    } else if (std::abs((float)dthtz) > 0.) {
//
//       rotate q0 about an axis parallel to the z-axis
//
        ca = std::cos(dthtz);
        sa = std::sin(dthtz);
//
        for (int n = 1; n <= 4; n++) {
        for (int l = lsta; l <= lend; l++) {
        for (int m = msta; m <= mend; m++) {
            q0rot(l,m,1,n) =  q0(l,m,1,n);
            temp           =  q0(l,m,2,n);
            q0rot(l,m,2,n) =  q0(l,m,2,n)*ca - q0(l,m,3,n)*sa;
            q0rot(l,m,3,n) =  temp        *sa + q0(l,m,3,n)*ca;
            q0rot(l,m,4,n) =  q0(l,m,4,n);
            q0rot(l,m,5,n) =  q0(l,m,5,n);
        }
        }
        }
//
    } else {
        for (int n = 1; n <= 4; n++) {
        for (int l = lsta; l <= lend; l++) {
        for (int m = msta; m <= mend; m++) {
            q0rot(l,m,1,n) =  q0(l,m,1,n);
            q0rot(l,m,2,n) =  q0(l,m,2,n);
            q0rot(l,m,3,n) =  q0(l,m,3,n);
            q0rot(l,m,4,n) =  q0(l,m,4,n);
            q0rot(l,m,5,n) =  q0(l,m,5,n);
        }
        }
        }
//
    }
//
    return;
}

} // namespace rotateq0_ns
