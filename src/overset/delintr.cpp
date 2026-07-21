// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "delintr.h"

namespace delintr_ns {

void delintr(int& jdimf, int& kdimf, int& idimf, int& jdimc, int& kdimc, int& idimc,
             FortranArray4DRef<double> deltjf, FortranArray4DRef<double> deltkf,
             FortranArray4DRef<double> deltif, FortranArray4DRef<double> deltjc,
             FortranArray4DRef<double> deltkc, FortranArray4DRef<double> deltic)
{
    int& i2d = cmn_twod.i2d;

    if (i2d == 0) {

        // --- deltif interpolation (3D case) ---
        // Copy coarse values to odd fine indices
        for (int jf = 1; jf <= jdimf; jf += 2) {
            int jc = (jf + 1) / 2;
            for (int kf = 1; kf <= kdimf; kf += 2) {
                int kc = (kf + 1) / 2;
                for (int ll = 1; ll <= 2; ll++) {
                    deltif(jf, kf, 1, ll) = deltic(jc, kc, 1, ll);
                    deltif(jf, kf, 2, ll) = deltic(jc, kc, 2, ll);
                    deltif(jf, kf, 3, ll) = deltic(jc, kc, 3, ll);
                }
            }
        }
        // Average even jf indices
        for (int jf = 2; jf <= jdimf; jf += 2) {
            for (int kf = 1; kf <= kdimf; kf += 2) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltif(jf, kf, 1, ll) = 0.5 * (deltif(jf-1, kf, 1, ll)
                                                  + deltif(jf+1, kf, 1, ll));
                    deltif(jf, kf, 2, ll) = 0.5 * (deltif(jf-1, kf, 2, ll)
                                                  + deltif(jf+1, kf, 2, ll));
                    deltif(jf, kf, 3, ll) = 0.5 * (deltif(jf-1, kf, 3, ll)
                                                  + deltif(jf+1, kf, 3, ll));
                }
            }
        }
        // Average even kf indices
        for (int jf = 1; jf <= jdimf; jf++) {
            for (int kf = 2; kf <= kdimf; kf += 2) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltif(jf, kf, 1, ll) = 0.5 * (deltif(jf, kf-1, 1, ll)
                                                  + deltif(jf, kf+1, 1, ll));
                    deltif(jf, kf, 2, ll) = 0.5 * (deltif(jf, kf-1, 2, ll)
                                                  + deltif(jf, kf+1, 2, ll));
                    deltif(jf, kf, 3, ll) = 0.5 * (deltif(jf, kf-1, 3, ll)
                                                  + deltif(jf, kf+1, 3, ll));
                }
            }
        }

        // --- deltjf interpolation (3D case) ---
        // Copy coarse values to odd fine indices
        for (int kf = 1; kf <= kdimf; kf += 2) {
            int kc = (kf + 1) / 2;
            for (int iff = 1; iff <= idimf; iff += 2) {
                int ic = (iff + 1) / 2;
                for (int ll = 1; ll <= 2; ll++) {
                    deltjf(kf, iff, 1, ll) = deltjc(kc, ic, 1, ll);
                    deltjf(kf, iff, 2, ll) = deltjc(kc, ic, 2, ll);
                    deltjf(kf, iff, 3, ll) = deltjc(kc, ic, 3, ll);
                }
            }
        }
        // Average even kf indices
        for (int kf = 2; kf <= kdimf; kf += 2) {
            for (int iff = 1; iff <= idimf; iff += 2) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltjf(kf, iff, 1, ll) = 0.5 * (deltjf(kf-1, iff, 1, ll)
                                                   + deltjf(kf+1, iff, 1, ll));
                    deltjf(kf, iff, 2, ll) = 0.5 * (deltjf(kf-1, iff, 2, ll)
                                                   + deltjf(kf+1, iff, 2, ll));
                    deltjf(kf, iff, 3, ll) = 0.5 * (deltjf(kf-1, iff, 3, ll)
                                                   + deltjf(kf+1, iff, 3, ll));
                }
            }
        }
        // Average even if indices
        for (int kf = 1; kf <= kdimf; kf++) {
            for (int iff = 2; iff <= idimf; iff += 2) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltjf(kf, iff, 1, ll) = 0.5 * (deltjf(kf, iff-1, 1, ll)
                                                   + deltjf(kf, iff+1, 1, ll));
                    deltjf(kf, iff, 2, ll) = 0.5 * (deltjf(kf, iff-1, 2, ll)
                                                   + deltjf(kf, iff+1, 2, ll));
                    deltjf(kf, iff, 3, ll) = 0.5 * (deltjf(kf, iff-1, 3, ll)
                                                   + deltjf(kf, iff+1, 3, ll));
                }
            }
        }

        // --- deltkf interpolation (3D case) ---
        // Copy coarse values to odd fine indices
        for (int jf = 1; jf <= jdimf; jf += 2) {
            int jc = (jf + 1) / 2;
            for (int iff = 1; iff <= idimf; iff += 2) {
                int ic = (iff + 1) / 2;
                for (int ll = 1; ll <= 2; ll++) {
                    deltkf(jf, iff, 1, ll) = deltkc(jc, ic, 1, ll);
                    deltkf(jf, iff, 2, ll) = deltkc(jc, ic, 2, ll);
                    deltkf(jf, iff, 3, ll) = deltkc(jc, ic, 3, ll);
                }
            }
        }
        // Average even jf indices
        for (int jf = 2; jf <= jdimf; jf += 2) {
            for (int iff = 1; iff <= idimf; iff += 2) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltkf(jf, iff, 1, ll) = 0.5 * (deltkf(jf-1, iff, 1, ll)
                                                   + deltkf(jf+1, iff, 1, ll));
                    deltkf(jf, iff, 2, ll) = 0.5 * (deltkf(jf-1, iff, 2, ll)
                                                   + deltkf(jf+1, iff, 2, ll));
                    deltkf(jf, iff, 3, ll) = 0.5 * (deltkf(jf-1, iff, 3, ll)
                                                   + deltkf(jf+1, iff, 3, ll));
                }
            }
        }
        // Average even if indices
        for (int jf = 1; jf <= jdimf; jf++) {
            for (int iff = 2; iff <= idimf; iff += 2) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltkf(jf, iff, 1, ll) = 0.5 * (deltkf(jf, iff-1, 1, ll)
                                                   + deltkf(jf, iff+1, 1, ll));
                    deltkf(jf, iff, 2, ll) = 0.5 * (deltkf(jf, iff-1, 2, ll)
                                                   + deltkf(jf, iff+1, 2, ll));
                    deltkf(jf, iff, 3, ll) = 0.5 * (deltkf(jf, iff-1, 3, ll)
                                                   + deltkf(jf, iff+1, 3, ll));
                }
            }
        }

    } else {

        // --- deltif interpolation (2D case) ---
        // Copy coarse values to odd fine indices
        for (int jf = 1; jf <= jdimf; jf += 2) {
            int jc = (jf + 1) / 2;
            for (int kf = 1; kf <= kdimf; kf += 2) {
                int kc = (kf + 1) / 2;
                for (int ll = 1; ll <= 2; ll++) {
                    deltif(jf, kf, 1, ll) = deltic(jc, kc, 1, ll);
                    deltif(jf, kf, 2, ll) = deltic(jc, kc, 2, ll);
                    deltif(jf, kf, 3, ll) = deltic(jc, kc, 3, ll);
                }
            }
        }
        // Average even jf indices
        for (int jf = 2; jf <= jdimf; jf += 2) {
            for (int kf = 1; kf <= kdimf; kf += 2) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltif(jf, kf, 1, ll) = 0.5 * (deltif(jf-1, kf, 1, ll)
                                                  + deltif(jf+1, kf, 1, ll));
                    deltif(jf, kf, 2, ll) = 0.5 * (deltif(jf-1, kf, 2, ll)
                                                  + deltif(jf+1, kf, 2, ll));
                    deltif(jf, kf, 3, ll) = 0.5 * (deltif(jf-1, kf, 3, ll)
                                                  + deltif(jf+1, kf, 3, ll));
                }
            }
        }
        // Average even kf indices
        for (int jf = 1; jf <= jdimf; jf++) {
            for (int kf = 2; kf <= kdimf; kf += 2) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltif(jf, kf, 1, ll) = 0.5 * (deltif(jf, kf-1, 1, ll)
                                                  + deltif(jf, kf+1, 1, ll));
                    deltif(jf, kf, 2, ll) = 0.5 * (deltif(jf, kf-1, 2, ll)
                                                  + deltif(jf, kf+1, 2, ll));
                    deltif(jf, kf, 3, ll) = 0.5 * (deltif(jf, kf-1, 3, ll)
                                                  + deltif(jf, kf+1, 3, ll));
                }
            }
        }

        // --- deltjf interpolation (2D case) ---
        // Copy coarse values to odd fine kf indices, loop i=1,2
        for (int kf = 1; kf <= kdimf; kf += 2) {
            int kc = (kf + 1) / 2;
            for (int i = 1; i <= 2; i++) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltjf(kf, i, 1, ll) = deltjc(kc, i, 1, ll);
                    deltjf(kf, i, 2, ll) = deltjc(kc, i, 2, ll);
                    deltjf(kf, i, 3, ll) = deltjc(kc, i, 3, ll);
                }
            }
        }
        // Average even kf indices, loop i=1,2
        for (int kf = 2; kf <= kdimf; kf += 2) {
            for (int i = 1; i <= 2; i++) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltjf(kf, i, 1, ll) = 0.5 * (deltjf(kf-1, i, 1, ll)
                                                 + deltjf(kf+1, i, 1, ll));
                    deltjf(kf, i, 2, ll) = 0.5 * (deltjf(kf-1, i, 2, ll)
                                                 + deltjf(kf+1, i, 2, ll));
                    deltjf(kf, i, 3, ll) = 0.5 * (deltjf(kf-1, i, 3, ll)
                                                 + deltjf(kf+1, i, 3, ll));
                }
            }
        }

        // --- deltkf interpolation (2D case) ---
        // Copy coarse values to odd fine jf indices, loop i=1,2
        for (int jf = 1; jf <= jdimf; jf += 2) {
            int jc = (jf + 1) / 2;
            for (int i = 1; i <= 2; i++) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltkf(jf, i, 1, ll) = deltkc(jc, i, 1, ll);
                    deltkf(jf, i, 2, ll) = deltkc(jc, i, 2, ll);
                    deltkf(jf, i, 3, ll) = deltkc(jc, i, 3, ll);
                }
            }
        }
        // Average even jf indices, loop i=1,2
        for (int jf = 2; jf <= jdimf; jf += 2) {
            for (int i = 1; i <= 2; i++) {
                for (int ll = 1; ll <= 2; ll++) {
                    deltkf(jf, i, 1, ll) = 0.5 * (deltkf(jf-1, i, 1, ll)
                                                 + deltkf(jf+1, i, 1, ll));
                    deltkf(jf, i, 2, ll) = 0.5 * (deltkf(jf-1, i, 2, ll)
                                                 + deltkf(jf+1, i, 2, ll));
                    deltkf(jf, i, 3, ll) = 0.5 * (deltkf(jf-1, i, 3, ll)
                                                 + deltkf(jf+1, i, 3, ll));
                }
            }
        }

    } // end if (i2d == 0)

} // delintr

} // namespace delintr_ns
