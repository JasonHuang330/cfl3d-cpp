// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rotateq2_d.h"
#include "rotate.h"
#include "rotateq.h"
#include <cmath>

namespace rotateq2_d_ns {

void rotate(int& jdim, int& kdim, int& idim, FortranArray2DRef<double> t, FortranArray3DRef<double> tti, FortranArray3DRef<double> ttj, FortranArray3DRef<double> ttk, FortranArray2DRef<double> t1, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& nbl, int& irot, double& rfreqr, double& omegx, double& omegy, double& omegz, double& xorg, double& yorg, double& zorg, double& thetax, double& thetay, double& thetaz, double& thxold, double& thyold, double& thzold, int& iupdat, double& time2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    rotate_ns::rotate(jdim, kdim, idim, t, tti, ttj, ttk, t1, x, y, z,
                      nbl, irot, rfreqr, omegx, omegy, omegz,
                      xorg, yorg, zorg, thetax, thetay, thetaz,
                      thxold, thyold, thzold, iupdat, time2,
                      nou, bou, nbuf, ibufdim);
}

void rotateq(int& jd, int& kd, int& id, FortranArray4DRef<double> q, FortranArray4DRef<double> qrot, int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend, double& dthtx, double& dthty, double& dthtz)
{
    rotateq_ns::rotateq(jd, kd, id, q, qrot, ista, iend, jsta, jend, ksta, kend,
                        dthtx, dthty, dthtz);
}

void rotateq2_d(int& jd, int& kd, FortranArray3DRef<double> q, FortranArray3DRef<double> qrot, int& jsta, int& jend, int& ksta, int& kend, double& dthtx, double& dthty, double& dthtz)
{
    double ca, sa, temp;

    if (std::abs((float)dthtx) > 0.f) {
//
//       rotate q about an axis parallel to the x-axis
//
        ca = std::cos(dthtx);
        sa = std::sin(dthtx);
//
        for (int j = jsta; j <= jend; j++) {
        for (int k = ksta; k <= kend; k++) {
            qrot(j,k,1) =  q(j,k,1);
            qrot(j,k,2) =  q(j,k,2);
            temp         =  q(j,k,3);
            qrot(j,k,3) =  q(j,k,3)*ca - q(j,k,4)*sa;
            qrot(j,k,4) =  temp      *sa + q(j,k,4)*ca;
            qrot(j,k,5) =  q(j,k,5);
        }
        }
//
    } else if (std::abs((float)dthty) > 0.f) {
//
//       rotate q about an axis parallel to the y-axis
//
        ca = std::cos(dthty);
        sa = std::sin(dthty);
//
        for (int j = jsta; j <= jend; j++) {
        for (int k = ksta; k <= kend; k++) {
            qrot(j,k,1) =  q(j,k,1);
            temp         =  q(j,k,2);
            qrot(j,k,2) =  q(j,k,2)*ca + q(j,k,4)*sa;
            qrot(j,k,3) =  q(j,k,3);
            qrot(j,k,4) = -temp      *sa + q(j,k,4)*ca;
            qrot(j,k,5) =  q(j,k,5);
        }
        }
//
    } else if (std::abs((float)dthtz) > 0.f) {
//
//       rotate q about an axis parallel to the z-axis
//
        ca = std::cos(dthtz);
        sa = std::sin(dthtz);
//
        for (int j = jsta; j <= jend; j++) {
        for (int k = ksta; k <= kend; k++) {
            qrot(j,k,1) =  q(j,k,1);
            temp         =  q(j,k,2);
            qrot(j,k,2) =  q(j,k,2)*ca - q(j,k,3)*sa;
            qrot(j,k,3) =  temp      *sa + q(j,k,3)*ca;
            qrot(j,k,4) =  q(j,k,4);
            qrot(j,k,5) =  q(j,k,5);
        }
        }
//
    } else {
//
//    no rotation
//
        for (int j = jsta; j <= jend; j++) {
        for (int k = ksta; k <= kend; k++) {
            qrot(j,k,1) =  q(j,k,1);
            qrot(j,k,2) =  q(j,k,2);
            qrot(j,k,3) =  q(j,k,3);
            qrot(j,k,4) =  q(j,k,4);
            qrot(j,k,5) =  q(j,k,5);
        }
        }
//
    }

    return;
}

} // namespace rotateq2_d_ns
