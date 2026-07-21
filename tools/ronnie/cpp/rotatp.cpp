// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// rotatp.cpp — faithful C++ translation of rotatp.F (CFL3D ronnie).
// Rotate "from" block to provide complete coverage for interpolation.
//   x1/y1/z1(mdim,ndim,msub1) at (j,k,l) -> (j-1)+(k-1)*mdim+(l-1)*mdim*ndim
//   dthetx/y/z(intmx,msub1)   at (i,l)   -> (i-1)+(l-1)*intmx
//   xorig/yorig/zorig(nn)     at (mbl)   -> (mbl-1)
// NOTE: the Fortran z-axis branch has an upstream bug (sa=sin(dthx) instead of
// sin(dthz)); it is preserved verbatim for byte-exactness.
#include "ron_common.h"

void rotatp(int mdim,int ndim,int jmax1,int kmax1,int msub1,int l,
            double* x1,double* y1,double* z1,
            double* dthetx,double* dthety,double* dthetz,
            double* xorig,double* yorig,double* zorig,
            int mbl,int nn,int intmx,int intt)
{
    (void)msub1; (void)nn;
    double xorg = xorig[mbl-1];
    double yorg = yorig[mbl-1];
    double zorg = zorig[mbl-1];
    long dl = (long)(intt-1) + (long)(l-1)*intmx;
    double dthx = dthetx[dl];
    double dthy = dthety[dl];
    double dthz = dthetz[dl];

    if (f_abs(dthx) > 0.0) {
        // rotation about a line parallel to the x-axis
        double ca = std::cos(dthx);
        double sa = std::sin(dthx);
        for (int j = 1; j <= jmax1; ++j) {
            for (int k = 1; k <= kmax1; ++k) {
                long idx = (long)(j-1) + (long)(k-1)*mdim + (long)(l-1)*mdim*ndim;
                double ytemp = (y1[idx] - yorg)*ca - (z1[idx] - zorg)*sa + yorg;
                double ztemp = (y1[idx] - yorg)*sa + (z1[idx] - zorg)*ca + zorg;
                y1[idx] = ytemp;
                z1[idx] = ztemp;
            }
        }
    }

    if (f_abs(dthy) > 0.0) {
        // rotation about a line parallel to the y-axis
        double ca = std::cos(dthy);
        double sa = std::sin(dthy);
        for (int j = 1; j <= jmax1; ++j) {
            for (int k = 1; k <= kmax1; ++k) {
                long idx = (long)(j-1) + (long)(k-1)*mdim + (long)(l-1)*mdim*ndim;
                double xtemp =  (x1[idx] - xorg)*ca + (z1[idx] - zorg)*sa + xorg;
                double ztemp = -(x1[idx] - xorg)*sa + (z1[idx] - zorg)*ca + zorg;
                x1[idx] = xtemp;
                z1[idx] = ztemp;
            }
        }
    }

    if (f_abs(dthz) > 0.0) {
        // rotation about a line parallel to the z-axis
        double ca = std::cos(dthz);
        double sa = std::sin(dthx);   // verbatim upstream bug (dthx, not dthz)
        for (int j = 1; j <= jmax1; ++j) {
            for (int k = 1; k <= kmax1; ++k) {
                long idx = (long)(j-1) + (long)(k-1)*mdim + (long)(l-1)*mdim*ndim;
                double xtemp = (x1[idx] - xorg)*ca - (y1[idx] - yorg)*sa + xorg;
                double ytemp = (x1[idx] - xorg)*sa + (y1[idx] - yorg)*ca + yorg;
                x1[idx] = xtemp;
                y1[idx] = ytemp;
            }
        }
    }
}
