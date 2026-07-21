// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "grdmove.h"
#include <cmath>

namespace grdmove_ns {

void grdmove(int& nbl, int& jdim, int& kdim, int& idim,
             FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
             double& xorig, double& yorig, double& zorig,
             double& xorg, double& yorg, double& zorg,
             double& thetax, double& thetay, double& thetaz)
{
    double sa, ca, xx, yy, zz;

    // rotation about x-axis
    sa = std::sin(thetax);
    ca = std::cos(thetax);
    for (int j = 1; j <= jdim; j++) {
        for (int k = 1; k <= kdim; k++) {
            for (int i = 1; i <= idim; i++) {
                yy = (y(j,k,i)-yorig)*ca-(z(j,k,i)-zorig)*sa+yorig;
                zz = (y(j,k,i)-yorig)*sa+(z(j,k,i)-zorig)*ca+zorig;
                y(j,k,i) = yy;
                z(j,k,i) = zz;
            }
        }
    }

    // rotation about y-axis
    sa = std::sin(thetay);
    ca = std::cos(thetay);
    for (int j = 1; j <= jdim; j++) {
        for (int k = 1; k <= kdim; k++) {
            for (int i = 1; i <= idim; i++) {
                xx =  (x(j,k,i)-xorig)*ca+(z(j,k,i)-zorig)*sa+xorig;
                zz = -(x(j,k,i)-xorig)*sa+(z(j,k,i)-zorig)*ca+zorig;
                x(j,k,i) = xx;
                z(j,k,i) = zz;
            }
        }
    }

    // rotation about z-axis
    sa = std::sin(thetaz);
    ca = std::cos(thetaz);
    for (int j = 1; j <= jdim; j++) {
        for (int k = 1; k <= kdim; k++) {
            for (int i = 1; i <= idim; i++) {
                xx = (x(j,k,i)-xorig)*ca-(y(j,k,i)-yorig)*sa+xorig;
                yy = (x(j,k,i)-xorig)*sa+(y(j,k,i)-yorig)*ca+yorig;
                x(j,k,i) = xx;
                y(j,k,i) = yy;
            }
        }
    }

    // translation
    for (int j = 1; j <= jdim; j++) {
        for (int k = 1; k <= kdim; k++) {
            for (int i = 1; i <= idim; i++) {
                x(j,k,i) = x(j,k,i) + xorg - xorig;
                y(j,k,i) = y(j,k,i) + yorg - yorig;
                z(j,k,i) = z(j,k,i) + zorg - zorig;
            }
        }
    }
}

} // namespace grdmove_ns
