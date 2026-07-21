// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "expand.h"
#include <cmath>

namespace expand_ns {

void expand(int& mdim, int& ndim, int& msub1, int& jmax1, int& kmax1, int& l,
            FortranArray2DRef<double> x1, FortranArray2DRef<double> y1, FortranArray2DRef<double> z1,
            FortranArray3DRef<double> xte, FortranArray3DRef<double> yte, FortranArray3DRef<double> zte,
            double& factjlo, double& factjhi, double& factklo, double& factkhi,
            int& jmax2, int& kmax2,
            FortranArray2DRef<double> x2, FortranArray2DRef<double> y2, FortranArray2DRef<double> z2)
{
    double sk1, sk2, sj1, sj2;
    double factk1, factk2, factj1, factj2;
    double dx, dy, dz, ds;

    sk1 = 0.;
    sk2 = 0.;
    for (int j = 2; j <= jmax1; j++) {
        sk1 = sk1 + std::sqrt((x1(j,1)-x1(j-1,1))*(x1(j,1)-x1(j-1,1))
                             +(y1(j,1)-y1(j-1,1))*(y1(j,1)-y1(j-1,1))
                             +(z1(j,1)-z1(j-1,1))*(z1(j,1)-z1(j-1,1)));
        sk2 = sk2 + std::sqrt((x1(j,kmax1)-x1(j-1,kmax1))*(x1(j,kmax1)-x1(j-1,kmax1))
                             +(y1(j,kmax1)-y1(j-1,kmax1))*(y1(j,kmax1)-y1(j-1,kmax1))
                             +(z1(j,kmax1)-z1(j-1,kmax1))*(z1(j,kmax1)-z1(j-1,kmax1)));
    }
    sj1 = 0.;
    sj2 = 0.;
    for (int k = 2; k <= kmax1; k++) {
        sj1 = sj1 + std::sqrt((x1(1,k)-x1(1,k-1))*(x1(1,k)-x1(1,k-1))
                             +(y1(1,k)-y1(1,k-1))*(y1(1,k)-y1(1,k-1))
                             +(z1(1,k)-z1(1,k-1))*(z1(1,k)-z1(1,k-1)));
        sj2 = sj2 + std::sqrt((x1(jmax1,k)-x1(jmax1,k-1))*(x1(jmax1,k)-x1(jmax1,k-1))
                             +(y1(jmax1,k)-y1(jmax1,k-1))*(y1(jmax1,k)-y1(jmax1,k-1))
                             +(z1(jmax1,k)-z1(jmax1,k-1))*(z1(jmax1,k)-z1(jmax1,k-1)));
    }
    factk1 = factklo * sk1;
    factk2 = factkhi * sk2;
    factj1 = factjlo * sj1;
    factj2 = factjhi * sj2;

    // Loop 210: expand at k=1 boundary
    for (int j = 1; j <= jmax1; j++) {
        dx = x1(j,2) - x1(j,1);
        dy = y1(j,2) - y1(j,1);
        dz = z1(j,2) - z1(j,1);
        ds = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (ds <= 0.) ds = 1.0;
        xte(j+1,1,l) = x1(j,1) - factk1*dx/ds;
        yte(j+1,1,l) = y1(j,1) - factk1*dy/ds;
        zte(j+1,1,l) = z1(j,1) - factk1*dz/ds;
    }

    // Loop 220: expand at k=kmax1+2 boundary
    for (int j = 1; j <= jmax1; j++) {
        dx = x1(j,kmax1) - x1(j,kmax1-1);
        dy = y1(j,kmax1) - y1(j,kmax1-1);
        dz = z1(j,kmax1) - z1(j,kmax1-1);
        ds = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (ds <= 0.) ds = 1.0;
        xte(j+1,kmax1+2,l) = x1(j,kmax1) + factk2*dx/ds;
        yte(j+1,kmax1+2,l) = y1(j,kmax1) + factk2*dy/ds;
        zte(j+1,kmax1+2,l) = z1(j,kmax1) + factk2*dz/ds;
    }

    // Loop 250: copy interior points
    for (int k = 2; k <= kmax1+1; k++) {
        for (int j = 1; j <= jmax1; j++) {
            xte(j+1,k,l) = x1(j,k-1);
            yte(j+1,k,l) = y1(j,k-1);
            zte(j+1,k,l) = z1(j,k-1);
        }
    }

    // Loop 310: expand at j=1 boundary
    for (int k = 1; k <= kmax1+2; k++) {
        dx = xte(3,k,l) - xte(2,k,l);
        dy = yte(3,k,l) - yte(2,k,l);
        dz = zte(3,k,l) - zte(2,k,l);
        ds = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (ds <= 0.) ds = 1.0;
        xte(1,k,l) = xte(2,k,l) - factj1*dx/ds;
        yte(1,k,l) = yte(2,k,l) - factj1*dy/ds;
        zte(1,k,l) = zte(2,k,l) - factj1*dz/ds;
    }

    // Loop 320: expand at j=jmax1+2 boundary
    for (int k = 1; k <= kmax1+2; k++) {
        dx = xte(jmax1+1,k,l) - xte(jmax1,k,l);
        dy = yte(jmax1+1,k,l) - yte(jmax1,k,l);
        dz = zte(jmax1+1,k,l) - zte(jmax1,k,l);
        ds = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (ds <= 0.) ds = 1.0;
        xte(jmax1+2,k,l) = xte(jmax1+1,k,l) + factj2*dx/ds;
        yte(jmax1+2,k,l) = yte(jmax1+1,k,l) + factj2*dy/ds;
        zte(jmax1+2,k,l) = zte(jmax1+1,k,l) + factj2*dz/ds;
    }

    return;
}

} // namespace expand_ns
