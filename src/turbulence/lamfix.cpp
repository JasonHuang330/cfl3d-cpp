// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
// Translated from Fortran module 'lamfix' to C++
#include "lamfix.h"
#include <cmath>

namespace lamfix_ns {

static inline double ccabs(double x) { return std::fabs(x); }

void lamfix(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> smin, int& mdim, int& ndim, FortranArray4DRef<double> bcdata, int& ibeg, int& iend, int& jbeg, int& jend, int& kbeg, int& kend, int& nface, int& llev)
{
    int i, j, k, ii, jj, kk;

//  J-dir:
    if (nface == 3) {
        jbeg = 1;
        for (i = ibeg; i <= iend-1; i++) {
            ii = i - ibeg + 1;
            for (k = kbeg; k <= kend-1; k++) {
                kk = k - kbeg + 1;
                if (bcdata(kk,ii,1,3) < 1.) {
                    jend = jdim;
                } else {
                    jend = (int)(bcdata(kk,ii,1,3) + .001);
                    if (llev == 2) jend = (jend+1)/2;
                    if (llev == 3) jend = (jend+3)/4;
                    if (llev == 4) jend = (jend+7)/8;
                    if (llev == 5) jend = (jend+15)/16;
                }
                for (j = jbeg; j <= jend-1; j++) {
                    smin(j,k,i) = -(ccabs(smin(j,k,i)));
                }
            }
        }
    }
    if (nface == 4) {
        jend = jdim;
        for (i = ibeg; i <= iend-1; i++) {
            ii = i - ibeg + 1;
            for (k = kbeg; k <= kend-1; k++) {
                kk = k - kbeg + 1;
                if (bcdata(kk,ii,1,3) < 1.) {
                    jbeg = 1;
                } else {
                    jbeg = (int)(bcdata(kk,ii,1,3) + .001) + 1;
                    if (llev == 2) jbeg = (jbeg+1)/2;
                    if (llev == 3) jbeg = (jbeg+3)/4;
                    if (llev == 4) jbeg = (jbeg+7)/8;
                    if (llev == 5) jbeg = (jbeg+15)/16;
                    jbeg = jdim - jbeg;
                }
                for (j = jbeg; j <= jend-1; j++) {
                    smin(j,k,i) = -(ccabs(smin(j,k,i)));
                }
            }
        }
    }
//  K-dir:
    if (nface == 5) {
        kbeg = 1;
        for (i = ibeg; i <= iend-1; i++) {
            ii = i - ibeg + 1;
            for (j = jbeg; j <= jend-1; j++) {
                jj = j - jbeg + 1;
                if (bcdata(jj,ii,1,3) < 1.) {
                    kend = kdim;
                } else {
                    kend = (int)(bcdata(jj,ii,1,3) + .001);
                    if (llev == 2) kend = (kend+1)/2;
                    if (llev == 3) kend = (kend+3)/4;
                    if (llev == 4) kend = (kend+7)/8;
                    if (llev == 5) kend = (kend+15)/16;
                }
                for (k = kbeg; k <= kend-1; k++) {
                    smin(j,k,i) = -(ccabs(smin(j,k,i)));
                }
            }
        }
    }
    if (nface == 6) {
        kend = kdim;
        for (i = ibeg; i <= iend-1; i++) {
            ii = i - ibeg + 1;
            for (j = jbeg; j <= jend-1; j++) {
                jj = j - jbeg + 1;
                if (bcdata(jj,ii,1,3) < 1.) {
                    kbeg = 1;
                } else {
                    kbeg = (int)(bcdata(jj,ii,1,3) + .001) + 1;
                    if (llev == 2) kbeg = (kbeg+1)/2;
                    if (llev == 3) kbeg = (kbeg+3)/4;
                    if (llev == 4) kbeg = (kbeg+7)/8;
                    if (llev == 5) kbeg = (kbeg+15)/16;
                    kbeg = kdim - kbeg;
                }
                for (k = kbeg; k <= kend-1; k++) {
                    smin(j,k,i) = -(ccabs(smin(j,k,i)));
                }
            }
        }
    }
//  I-dir:
    if (nface == 1) {
        ibeg = 1;
        for (k = kbeg; k <= kend-1; k++) {
            kk = k - kbeg + 1;
            for (j = jbeg; j <= jend-1; j++) {
                jj = j - jbeg + 1;
                if (bcdata(jj,kk,1,3) < 1.) {
                    iend = idim;
                } else {
                    iend = (int)(bcdata(jj,kk,1,3) + .001);
                    if (llev == 2) iend = (iend+1)/2;
                    if (llev == 3) iend = (iend+3)/4;
                    if (llev == 4) iend = (iend+7)/8;
                    if (llev == 5) iend = (iend+15)/16;
                }
                for (i = ibeg; i <= iend-1; i++) {
                    smin(j,k,i) = -(ccabs(smin(j,k,i)));
                }
            }
        }
    }
    if (nface == 2) {
        iend = idim;
        for (k = kbeg; k <= kend-1; k++) {
            kk = k - kbeg + 1;
            for (j = jbeg; j <= jend-1; j++) {
                jj = j - jbeg + 1;
                if (bcdata(jj,kk,1,3) < 1.) {
                    ibeg = 1;
                } else {
                    ibeg = (int)(bcdata(jj,kk,1,3) + .001) + 1;
                    if (llev == 2) ibeg = (ibeg+1)/2;
                    if (llev == 3) ibeg = (ibeg+3)/4;
                    if (llev == 4) ibeg = (ibeg+7)/8;
                    if (llev == 5) ibeg = (ibeg+15)/16;
                    ibeg = idim - ibeg;
                }
                for (i = ibeg; i <= iend-1; i++) {
                    smin(j,k,i) = -(ccabs(smin(j,k,i)));
                }
            }
        }
    }
    return;
}

} // namespace lamfix_ns
