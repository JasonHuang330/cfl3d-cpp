// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "xyzintr.h"

namespace xyzintr_ns {

void xyzintr(int& jdimf, int& kdimf, int& idimf, int& jdimc, int& kdimc, int& idimc,
             FortranArray3DRef<double> xf, FortranArray3DRef<double> yf, FortranArray3DRef<double> zf,
             FortranArray3DRef<double> xc, FortranArray3DRef<double> yc, FortranArray3DRef<double> zc)
{
    // common /twod/ i2d
    int& i2d = cmn_twod.i2d;

    if (i2d == 0) {
        // 3D case

        // Copy coarse grid values to fine grid at odd indices
        for (int jf = 1; jf <= jdimf; jf += 2) {
            int jc = (jf + 1) / 2;
            for (int kf = 1; kf <= kdimf; kf += 2) {
                int kc = (kf + 1) / 2;
                for (int iff = 1; iff <= idimf; iff += 2) {
                    int ic = (iff + 1) / 2;
                    xf(jf, kf, iff) = xc(jc, kc, ic);
                    yf(jf, kf, iff) = yc(jc, kc, ic);
                    zf(jf, kf, iff) = zc(jc, kc, ic);
                }
            }
        }

        // Interpolate in j-direction (even jf, odd kf, odd if)
        for (int jf = 2; jf <= jdimf; jf += 2) {
            for (int kf = 1; kf <= kdimf; kf += 2) {
                for (int iff = 1; iff <= idimf; iff += 2) {
                    xf(jf, kf, iff) = 0.5 * (xf(jf-1, kf, iff) + xf(jf+1, kf, iff));
                    yf(jf, kf, iff) = 0.5 * (yf(jf-1, kf, iff) + yf(jf+1, kf, iff));
                    zf(jf, kf, iff) = 0.5 * (zf(jf-1, kf, iff) + zf(jf+1, kf, iff));
                }
            }
        }

        // Interpolate in k-direction (all jf, even kf, odd if)
        for (int jf = 1; jf <= jdimf; jf++) {
            for (int kf = 2; kf <= kdimf; kf += 2) {
                for (int iff = 1; iff <= idimf; iff += 2) {
                    xf(jf, kf, iff) = 0.5 * (xf(jf, kf-1, iff) + xf(jf, kf+1, iff));
                    yf(jf, kf, iff) = 0.5 * (yf(jf, kf-1, iff) + yf(jf, kf+1, iff));
                    zf(jf, kf, iff) = 0.5 * (zf(jf, kf-1, iff) + zf(jf, kf+1, iff));
                }
            }
        }

        // Interpolate in i-direction (all jf, all kf, even if)
        for (int jf = 1; jf <= jdimf; jf++) {
            for (int kf = 1; kf <= kdimf; kf++) {
                for (int iff = 2; iff <= idimf; iff += 2) {
                    xf(jf, kf, iff) = 0.5 * (xf(jf, kf, iff-1) + xf(jf, kf, iff+1));
                    yf(jf, kf, iff) = 0.5 * (yf(jf, kf, iff-1) + yf(jf, kf, iff+1));
                    zf(jf, kf, iff) = 0.5 * (zf(jf, kf, iff-1) + zf(jf, kf, iff+1));
                }
            }
        }

    } else {
        // 2D case (i2d != 0)

        // Copy coarse grid values to fine grid at odd indices (i=1,2)
        for (int jf = 1; jf <= jdimf; jf += 2) {
            int jc = (jf + 1) / 2;
            for (int kf = 1; kf <= kdimf; kf += 2) {
                int kc = (kf + 1) / 2;
                for (int i = 1; i <= 2; i++) {
                    xf(jf, kf, i) = xc(jc, kc, i);
                    yf(jf, kf, i) = yc(jc, kc, i);
                    zf(jf, kf, i) = zc(jc, kc, i);
                }
            }
        }

        // Interpolate in j-direction (even jf, odd kf, i=1,2)
        for (int jf = 2; jf <= jdimf; jf += 2) {
            for (int kf = 1; kf <= kdimf; kf += 2) {
                for (int i = 1; i <= 2; i++) {
                    xf(jf, kf, i) = 0.5 * (xf(jf-1, kf, i) + xf(jf+1, kf, i));
                    yf(jf, kf, i) = 0.5 * (yf(jf-1, kf, i) + yf(jf+1, kf, i));
                    zf(jf, kf, i) = 0.5 * (zf(jf-1, kf, i) + zf(jf+1, kf, i));
                }
            }
        }

        // Interpolate in k-direction (all jf, even kf, i=1,2)
        for (int jf = 1; jf <= jdimf; jf++) {
            for (int kf = 2; kf <= kdimf; kf += 2) {
                for (int i = 1; i <= 2; i++) {
                    xf(jf, kf, i) = 0.5 * (xf(jf, kf-1, i) + xf(jf, kf+1, i));
                    yf(jf, kf, i) = 0.5 * (yf(jf, kf-1, i) + yf(jf, kf+1, i));
                    zf(jf, kf, i) = 0.5 * (zf(jf, kf-1, i) + zf(jf, kf+1, i));
                }
            }
        }

    }

    return;
}

} // namespace xyzintr_ns
