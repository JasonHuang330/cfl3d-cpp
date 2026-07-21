// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// dsmin.cpp — faithful C++ translation of dsmin.F (CFL3D ronnie).
// Find closest point in grid to point (xc,yc,zc).
//   lflag > 0 : search over all "from" blocks not yet searched
//               (lout(l)>0 indicates block l has already been searched)
//   lflag < 0 : search only block l1
//
// Arrays are 1-based column-major (CONTRACT rule 2).
//   x/y/z(jdim,kdim,nsub) at (j,k,l) -> (j-1)+(k-1)*jdim+(l-1)*jdim*kdim
//   jjmax/kkmax/lout/xif1.. (nsub) at (l) -> (l-1)
#include "ron_common.h"

void dsmin(int jdim,int kdim,int nsub,int* jjmax,int* kkmax,int lmax,
           double* x,double* y,double* z,double xc,double yc,double zc,
           int& j1,int& k1,int& l1,int* lout,int lflag,
           int* xif1,int* xif2,int* etf1,int* etf2)
{
    (void)nsub;
    double dmin = 1.0e+20;

    int ls = 1;
    int le = lmax;
    if (lflag < 0) {
        ls = l1;
        le = l1;
    }

    for (int l = ls; l <= le; ++l) {
        // skip over blocks already tried
        if (ls != le) {
            if (lout[l-1] > 0) continue;    // go to 1235
        }

        // search only over specified range, inside non-expanded limits
        int js = xif1[l-1];
        int je = xif2[l-1];
        int ks = etf1[l-1];
        int ke = etf2[l-1];
        je = i_min(je-2, jjmax[l-1]-2);
        ke = i_min(ke-2, kkmax[l-1]-2);
        js = i_max(js+1, 2);
        ks = i_max(ks+1, 2);

        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                long idx = (long)(j-1) + (long)(k-1)*jdim + (long)(l-1)*jdim*kdim;
                double dxv = xc - x[idx];
                double dyv = yc - y[idx];
                double dzv = zc - z[idx];
                double d1 = dxv*dxv + dyv*dyv + dzv*dzv;
                if (d1 < dmin) {
                    j1   = j;
                    k1   = k;
                    l1   = l;
                    dmin = d1;
                }
            }
        }
    }
}
