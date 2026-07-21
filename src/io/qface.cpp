// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "qface.h"

namespace qface_ns {

void qface(int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> dum,
           FortranArray4DRef<double> dumj0,
           FortranArray4DRef<double> dumk0,
           FortranArray4DRef<double> dumi0,
           FortranArray3DRef<double> bcj,
           FortranArray3DRef<double> bck,
           FortranArray3DRef<double> bci,
           FortranArray3DRef<double> blank,
           int& ldim)
{
    // common /twod/ i2d
    int& i2d = cmn_twod.i2d;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    for (int l = 1; l <= ldim; l++) {

        // j - boundaries
        for (int m = 1; m <= 2; m++) {

            int imin = 1;
            int imax = idim - 1;
            int kmin = 1;
            int kmax = kdim - 1;

            if (m == 1) {
                // left boundary
                for (int ii = imin; ii <= imax; ii++) {
                    for (int kk = kmin; kk <= kmax; kk++) {
                        double aa = 1. + bcj(kk, ii, m);
                        double bb = 1. - bcj(kk, ii, m);
                        double cc = 0.;
                        if (blank(1, kk, ii) == 0.) {
                            aa = 0.;
                            bb = 3.;
                            cc = -1.;
                        }
                        dumj0(kk, ii, l, 2) = 0.5 * (aa * dumj0(kk, ii, l, 1)
                                                    + bb * dum(1, kk, ii, l)
                                                    + cc * dum(2, kk, ii, l));
                    }
                }
            } else {
                // right boundary
                for (int ii = imin; ii <= imax; ii++) {
                    for (int kk = kmin; kk <= kmax; kk++) {
                        double aa = 1. + bcj(kk, ii, m);
                        double bb = 1. - bcj(kk, ii, m);
                        double cc = 0.;
                        if (blank(jdim - 1, kk, ii) == 0.) {
                            aa = 0.;
                            bb = 3.;
                            cc = -1.;
                        }
                        dumj0(kk, ii, l, 4) = 0.5 * (aa * dumj0(kk, ii, l, 3)
                                                    + bb * dum(jdim - 1, kk, ii, l)
                                                    + cc * dum(jdim - 2, kk, ii, l));
                    }
                }
            }
        } // end m loop (j-boundaries)

        // k - boundaries
        for (int m = 1; m <= 2; m++) {

            int imin = 1;
            int imax = idim - 1;
            int jmin = 1;
            int jmax = jdim - 1;

            if (m == 1) {
                // left boundary
                for (int ii = imin; ii <= imax; ii++) {
                    for (int jj = jmin; jj <= jmax; jj++) {
                        double aa = 1. + bck(jj, ii, m);
                        double bb = 1. - bck(jj, ii, m);
                        double cc = 0.;
                        if (blank(jj, 1, ii) == 0.) {
                            aa = 0.;
                            bb = 3.;
                            cc = -1.;
                        }
                        dumk0(jj, ii, l, 2) = 0.5 * (aa * dumk0(jj, ii, l, 1)
                                                    + bb * dum(jj, 1, ii, l)
                                                    + cc * dum(jj, 2, ii, l));
                    }
                }
            } else {
                // right boundary
                for (int ii = imin; ii <= imax; ii++) {
                    for (int jj = jmin; jj <= jmax; jj++) {
                        double aa = 1. + bck(jj, ii, m);
                        double bb = 1. - bck(jj, ii, m);
                        double cc = 0.;
                        if (blank(jj, kdim - 1, ii) == 0.) {
                            aa = 0.;
                            bb = 3.;
                            cc = -1.;
                        }
                        dumk0(jj, ii, l, 4) = 0.5 * (aa * dumk0(jj, ii, l, 3)
                                                    + bb * dum(jj, kdim - 1, ii, l)
                                                    + cc * dum(jj, kdim - 2, ii, l));
                    }
                }
            }
        } // end m loop (k-boundaries)

        // i - boundaries
        for (int m = 1; m <= 2; m++) {

            int kmin = 1;
            int kmax = kdim - 1;
            int jmin = 1;
            int jmax = jdim - 1;
            int i2 = 2;
            if (i2d == 1) i2 = 1;
            int id2 = idim - 2;
            if (i2d == 1 || idim == 2) id2 = idim - 1;

            if (m == 1) {
                // left boundary
                for (int kk = kmin; kk <= kmax; kk++) {
                    for (int jj = jmin; jj <= jmax; jj++) {
                        double aa = 1. + bci(jj, kk, m);
                        double bb = 1. - bci(jj, kk, m);
                        double cc = 0.;
                        if (blank(jj, kk, 1) == 0.) {
                            aa = 0.;
                            bb = 3.;
                            cc = -1.;
                        }
                        dumi0(jj, kk, l, 2) = 0.5 * (aa * dumi0(jj, kk, l, 1)
                                                    + bb * dum(jj, kk, 1, l)
                                                    + cc * dum(jj, kk, i2, l));
                    }
                }
            } else {
                // right boundary
                for (int kk = kmin; kk <= kmax; kk++) {
                    for (int jj = jmin; jj <= jmax; jj++) {
                        double aa = 1. + bci(jj, kk, m);
                        double bb = 1. - bci(jj, kk, m);
                        double cc = 0.;
                        if (blank(jj, kk, idim - 1) == 0.) {
                            aa = 0.;
                            bb = 3.;
                            cc = -1.;
                        }
                        dumi0(jj, kk, l, 4) = 0.5 * (aa * dumi0(jj, kk, l, 3)
                                                    + bb * dum(jj, kk, idim - 1, l)
                                                    + cc * dum(jj, kk, id2, l));
                    }
                }
            }
        } // end m loop (i-boundaries)

    } // end l loop

} // qface

} // namespace qface_ns
