// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// transp.cpp — faithful C++ translation of transp.F (CFL3D ronnie).
// Translate "from" block to provide complete coverage for interpolation.
//   x1/y1/z1(mdim,ndim,msub1) at (j,k,l) -> (j-1)+(k-1)*mdim+(l-1)*mdim*ndim
//   dx/dy/dz(intmx,msub1)     at (i,l)   -> (i-1)+(l-1)*intmx
// `intt` is Fortran arg `int` (renamed; `int` is a C++ keyword).
#include "ron_common.h"

void transp(int mdim,int ndim,int jmax1,int kmax1,int msub1,int l,
            double* x1,double* y1,double* z1,double* dx,double* dy,double* dz,
            int intmx,int intt)
{
    (void)msub1;
    long dl = (long)(intt-1) + (long)(l-1)*intmx;
    double dxv = dx[dl];
    double dyv = dy[dl];
    double dzv = dz[dl];
    for (int j = 1; j <= jmax1; ++j) {
        for (int k = 1; k <= kmax1; ++k) {
            long idx = (long)(j-1) + (long)(k-1)*mdim + (long)(l-1)*mdim*ndim;
            x1[idx] += dxv;
            y1[idx] += dyv;
            z1[idx] += dzv;
        }
    }
}
