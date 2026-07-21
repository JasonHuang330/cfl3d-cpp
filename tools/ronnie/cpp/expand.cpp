// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// expand.cpp — faithful C++ translation of expand.F (CFL3D ronnie).
// Expand "from" grid(s) at boundaries so the "to" grid is completely covered.
//   x1/y1/z1(mdim,ndim)        at (j,k)   -> (j-1)+(k-1)*mdim
//   x2/y2/z2(mdim,ndim)        at (j,k)   -> (j-1)+(k-1)*mdim   (unused here)
//   xte/yte/zte(mdim,ndim,msub1) at (j,k,l)-> (j-1)+(k-1)*mdim+(l-1)*mdim*ndim
#include "ron_common.h"

void expand(int mdim,int ndim,int msub1,int jmax1,int kmax1,int l,
            double* x1,double* y1,double* z1,double* xte,double* yte,double* zte,
            double factjlo,double factjhi,double factklo,double factkhi,
            int jmax2,int kmax2,double* x2,double* y2,double* z2)
{
    (void)msub1; (void)jmax2; (void)kmax2; (void)x2; (void)y2; (void)z2;
    auto I2  = [&](int j,int k){ return (long)(j-1) + (long)(k-1)*mdim; };
    auto I3  = [&](int j,int k){ return (long)(j-1) + (long)(k-1)*mdim
                                        + (long)(l-1)*mdim*ndim; };

    double sk1 = 0.0, sk2 = 0.0;
    for (int j = 2; j <= jmax1; ++j) {
        double dx1 = x1[I2(j,1)]-x1[I2(j-1,1)];
        double dy1 = y1[I2(j,1)]-y1[I2(j-1,1)];
        double dz1 = z1[I2(j,1)]-z1[I2(j-1,1)];
        sk1 += std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
        double dx2 = x1[I2(j,kmax1)]-x1[I2(j-1,kmax1)];
        double dy2 = y1[I2(j,kmax1)]-y1[I2(j-1,kmax1)];
        double dz2 = z1[I2(j,kmax1)]-z1[I2(j-1,kmax1)];
        sk2 += std::sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
    }
    double sj1 = 0.0, sj2 = 0.0;
    for (int k = 2; k <= kmax1; ++k) {
        double dx1 = x1[I2(1,k)]-x1[I2(1,k-1)];
        double dy1 = y1[I2(1,k)]-y1[I2(1,k-1)];
        double dz1 = z1[I2(1,k)]-z1[I2(1,k-1)];
        sj1 += std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
        double dx2 = x1[I2(jmax1,k)]-x1[I2(jmax1,k-1)];
        double dy2 = y1[I2(jmax1,k)]-y1[I2(jmax1,k-1)];
        double dz2 = z1[I2(jmax1,k)]-z1[I2(jmax1,k-1)];
        sj2 += std::sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
    }
    double factk1 = factklo*sk1;
    double factk2 = factkhi*sk2;
    double factj1 = factjlo*sj1;
    double factj2 = factjhi*sj2;

    // do 210: k=1 boundary
    for (int j = 1; j <= jmax1; ++j) {
        double dx = x1[I2(j,2)]-x1[I2(j,1)];
        double dy = y1[I2(j,2)]-y1[I2(j,1)];
        double dz = z1[I2(j,2)]-z1[I2(j,1)];
        double ds = std::sqrt(dx*dx+dy*dy+dz*dz);
        if (ds <= 0.0) ds = 1.0;
        xte[I3(j+1,1)] = x1[I2(j,1)]-factk1*dx/ds;
        yte[I3(j+1,1)] = y1[I2(j,1)]-factk1*dy/ds;
        zte[I3(j+1,1)] = z1[I2(j,1)]-factk1*dz/ds;
    }

    // do 220: k=kmax1 boundary
    for (int j = 1; j <= jmax1; ++j) {
        double dx = x1[I2(j,kmax1)]-x1[I2(j,kmax1-1)];
        double dy = y1[I2(j,kmax1)]-y1[I2(j,kmax1-1)];
        double dz = z1[I2(j,kmax1)]-z1[I2(j,kmax1-1)];
        double ds = std::sqrt(dx*dx+dy*dy+dz*dz);
        if (ds <= 0.0) ds = 1.0;
        xte[I3(j+1,kmax1+2)] = x1[I2(j,kmax1)]+factk2*dx/ds;
        yte[I3(j+1,kmax1+2)] = y1[I2(j,kmax1)]+factk2*dy/ds;
        zte[I3(j+1,kmax1+2)] = z1[I2(j,kmax1)]+factk2*dz/ds;
    }

    // do 250: copy interior
    for (int k = 2; k <= kmax1+1; ++k) {
        for (int j = 1; j <= jmax1; ++j) {
            xte[I3(j+1,k)] = x1[I2(j,k-1)];
            yte[I3(j+1,k)] = y1[I2(j,k-1)];
            zte[I3(j+1,k)] = z1[I2(j,k-1)];
        }
    }

    // do 310: j=1 boundary
    for (int k = 1; k <= kmax1+2; ++k) {
        double dx = xte[I3(3,k)]-xte[I3(2,k)];
        double dy = yte[I3(3,k)]-yte[I3(2,k)];
        double dz = zte[I3(3,k)]-zte[I3(2,k)];
        double ds = std::sqrt(dx*dx+dy*dy+dz*dz);
        if (ds <= 0.0) ds = 1.0;
        xte[I3(1,k)] = xte[I3(2,k)]-factj1*dx/ds;
        yte[I3(1,k)] = yte[I3(2,k)]-factj1*dy/ds;
        zte[I3(1,k)] = zte[I3(2,k)]-factj1*dz/ds;
    }

    // do 320: j=jmax1 boundary
    for (int k = 1; k <= kmax1+2; ++k) {
        double dx = xte[I3(jmax1+1,k)]-xte[I3(jmax1,k)];
        double dy = yte[I3(jmax1+1,k)]-yte[I3(jmax1,k)];
        double dz = zte[I3(jmax1+1,k)]-zte[I3(jmax1,k)];
        double ds = std::sqrt(dx*dx+dy*dy+dz*dz);
        if (ds <= 0.0) ds = 1.0;
        xte[I3(jmax1+2,k)] = xte[I3(jmax1+1,k)]+factj2*dx/ds;
        yte[I3(jmax1+2,k)] = yte[I3(jmax1+1,k)]+factj2*dy/ds;
        zte[I3(jmax1+2,k)] = zte[I3(jmax1+1,k)]+factj2*dz/ds;
    }
}
