// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "delv.h"
#include <cmath>
#include <algorithm>

namespace delv_ns {

void delv(int& jdim, int& kdim, int& idim,
          FortranArray4DRef<double> q,
          FortranArray4DRef<double> sj,
          FortranArray4DRef<double> sk,
          FortranArray4DRef<double> si,
          FortranArray3DRef<double> vol,
          FortranArray4DRef<double> ux,
          FortranArray3DRef<double> wt,
          FortranArray3DRef<double> blank,
          int& iover,
          FortranArray4DRef<double> qj0,
          FortranArray4DRef<double> qk0,
          FortranArray4DRef<double> qi0,
          FortranArray3DRef<double> bcj,
          FortranArray3DRef<double> bck,
          FortranArray3DRef<double> bci,
          int& nbl,
          FortranArray3DRef<double> volj0,
          FortranArray3DRef<double> volk0,
          FortranArray3DRef<double> voli0,
          int& maxbl,
          FortranArray1DRef<double> vormax,
          FortranArray1DRef<int> ivmax,
          FortranArray1DRef<int> jvmax,
          FortranArray1DRef<int> kvmax)
{
    int& i2d = cmn_twod.i2d;

    int ifi = 1;
    int ifj = 1;
    int ifk = 1;
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    // Initialize ux to zero
    for (int l = 1; l <= 9; l++)
        for (int i = 1; i <= idim1; i++)
            for (int k = 1; k <= kdim1; k++)
                for (int j = 1; j <= jdim1; j++)
                    ux(j, k, i, l) = 0.0;

    double onec = 1.0;

    //
    // J-direction contributions
    //
    if (ifj > 0) {
        for (int i = 1; i <= idim1; i++) {
            for (int k = 1; k <= kdim1; k++) {

                // cycle through interfaces (interior)
                for (int j = 2; j <= jdim1; j++) {
                    double term = sj(j, k, i, 4) / (vol(j, k, i) + vol(j-1, k, i));
                    wt(j, 1, 1) = term * ((q(j, k, i, 2) - q(j-1, k, i, 2)) * sj(j, k, i, 1));
                    wt(j, 1, 2) = term * ((q(j, k, i, 2) - q(j-1, k, i, 2)) * sj(j, k, i, 2));
                    wt(j, 1, 3) = term * ((q(j, k, i, 2) - q(j-1, k, i, 2)) * sj(j, k, i, 3));
                    wt(j, 1, 4) = term * ((q(j, k, i, 3) - q(j-1, k, i, 3)) * sj(j, k, i, 1));
                    wt(j, 1, 5) = term * ((q(j, k, i, 3) - q(j-1, k, i, 3)) * sj(j, k, i, 2));
                    wt(j, 1, 6) = term * ((q(j, k, i, 3) - q(j-1, k, i, 3)) * sj(j, k, i, 3));
                    wt(j, 1, 7) = term * ((q(j, k, i, 4) - q(j-1, k, i, 4)) * sj(j, k, i, 1));
                    wt(j, 1, 8) = term * ((q(j, k, i, 4) - q(j-1, k, i, 4)) * sj(j, k, i, 2));
                    wt(j, 1, 9) = term * ((q(j, k, i, 4) - q(j-1, k, i, 4)) * sj(j, k, i, 3));
                }

                // j=1 boundary
                {
                    int j = 1;
                    int jp1 = std::min(jdim - 1, j + 1);
                    double term = sj(j, k, i, 4) / (volj0(k, i, 1) + vol(j, k, i));
                    double factor = bcj(k, i, 1) + 1.0;
                    term = term * factor;
                    wt(j, 1, 1) = term * ((q(j, k, i, 2) - qj0(k, i, 2, 1)) * sj(j, k, i, 1));
                    wt(j, 1, 2) = term * ((q(j, k, i, 2) - qj0(k, i, 2, 1)) * sj(j, k, i, 2));
                    wt(j, 1, 3) = term * ((q(j, k, i, 2) - qj0(k, i, 2, 1)) * sj(j, k, i, 3));
                    wt(j, 1, 4) = term * ((q(j, k, i, 3) - qj0(k, i, 3, 1)) * sj(j, k, i, 1));
                    wt(j, 1, 5) = term * ((q(j, k, i, 3) - qj0(k, i, 3, 1)) * sj(j, k, i, 2));
                    wt(j, 1, 6) = term * ((q(j, k, i, 3) - qj0(k, i, 3, 1)) * sj(j, k, i, 3));
                    wt(j, 1, 7) = term * ((q(j, k, i, 4) - qj0(k, i, 4, 1)) * sj(j, k, i, 1));
                    wt(j, 1, 8) = term * ((q(j, k, i, 4) - qj0(k, i, 4, 1)) * sj(j, k, i, 2));
                    wt(j, 1, 9) = term * ((q(j, k, i, 4) - qj0(k, i, 4, 1)) * sj(j, k, i, 3));
                }

                // j=jdim boundary
                {
                    int j = jdim;
                    int jm2 = std::max(1, j - 2);
                    double term = sj(j, k, i, 4) / (volj0(k, i, 3) + vol(j-1, k, i));
                    double factor = bcj(k, i, 2) + 1.0;
                    term = term * factor;
                    wt(j, 1, 1) = term * ((qj0(k, i, 2, 3) - q(j-1, k, i, 2)) * sj(j, k, i, 1));
                    wt(j, 1, 2) = term * ((qj0(k, i, 2, 3) - q(j-1, k, i, 2)) * sj(j, k, i, 2));
                    wt(j, 1, 3) = term * ((qj0(k, i, 2, 3) - q(j-1, k, i, 2)) * sj(j, k, i, 3));
                    wt(j, 1, 4) = term * ((qj0(k, i, 3, 3) - q(j-1, k, i, 3)) * sj(j, k, i, 1));
                    wt(j, 1, 5) = term * ((qj0(k, i, 3, 3) - q(j-1, k, i, 3)) * sj(j, k, i, 2));
                    wt(j, 1, 6) = term * ((qj0(k, i, 3, 3) - q(j-1, k, i, 3)) * sj(j, k, i, 3));
                    wt(j, 1, 7) = term * ((qj0(k, i, 4, 3) - q(j-1, k, i, 4)) * sj(j, k, i, 1));
                    wt(j, 1, 8) = term * ((qj0(k, i, 4, 3) - q(j-1, k, i, 4)) * sj(j, k, i, 2));
                    wt(j, 1, 9) = term * ((qj0(k, i, 4, 3) - q(j-1, k, i, 4)) * sj(j, k, i, 3));
                }

                // cycle through cell centers
                for (int l = 1; l <= 9; l++)
                    for (int j = 1; j <= jdim1; j++)
                        ux(j, k, i, l) = ux(j, k, i, l) + wt(j, 1, l) + wt(j+1, 1, l);

            } // k
        } // i
    } // ifj

    //
    // K-direction contributions
    //
    if (ifk > 0) {
        for (int i = 1; i <= idim1; i++) {

            // cycle through interfaces (interior)
            for (int k = 2; k <= kdim - 1; k++) {
                for (int j = 1; j <= jdim1; j++) {
                    double term = sk(j, k, i, 4) / (vol(j, k, i) + vol(j, k-1, i));
                    wt(j, k, 1) = term * ((q(j, k, i, 2) - q(j, k-1, i, 2)) * sk(j, k, i, 1));
                    wt(j, k, 2) = term * ((q(j, k, i, 2) - q(j, k-1, i, 2)) * sk(j, k, i, 2));
                    wt(j, k, 3) = term * ((q(j, k, i, 2) - q(j, k-1, i, 2)) * sk(j, k, i, 3));
                    wt(j, k, 4) = term * ((q(j, k, i, 3) - q(j, k-1, i, 3)) * sk(j, k, i, 1));
                    wt(j, k, 5) = term * ((q(j, k, i, 3) - q(j, k-1, i, 3)) * sk(j, k, i, 2));
                    wt(j, k, 6) = term * ((q(j, k, i, 3) - q(j, k-1, i, 3)) * sk(j, k, i, 3));
                    wt(j, k, 7) = term * ((q(j, k, i, 4) - q(j, k-1, i, 4)) * sk(j, k, i, 1));
                    wt(j, k, 8) = term * ((q(j, k, i, 4) - q(j, k-1, i, 4)) * sk(j, k, i, 2));
                    wt(j, k, 9) = term * ((q(j, k, i, 4) - q(j, k-1, i, 4)) * sk(j, k, i, 3));
                }
            }

            // k=1 boundary
            {
                int k = 1;
                int kp1 = std::min(kdim - 1, k + 1);
                for (int j = 1; j <= jdim1; j++) {
                    double term = sk(j, k, i, 4) / (volk0(j, i, 1) + vol(j, k, i));
                    double factor = bck(j, i, 1) + 1.0;
                    term = term * factor;
                    wt(j, k, 1) = term * ((q(j, k, i, 2) - qk0(j, i, 2, 1)) * sk(j, k, i, 1));
                    wt(j, k, 2) = term * ((q(j, k, i, 2) - qk0(j, i, 2, 1)) * sk(j, k, i, 2));
                    wt(j, k, 3) = term * ((q(j, k, i, 2) - qk0(j, i, 2, 1)) * sk(j, k, i, 3));
                    wt(j, k, 4) = term * ((q(j, k, i, 3) - qk0(j, i, 3, 1)) * sk(j, k, i, 1));
                    wt(j, k, 5) = term * ((q(j, k, i, 3) - qk0(j, i, 3, 1)) * sk(j, k, i, 2));
                    wt(j, k, 6) = term * ((q(j, k, i, 3) - qk0(j, i, 3, 1)) * sk(j, k, i, 3));
                    wt(j, k, 7) = term * ((q(j, k, i, 4) - qk0(j, i, 4, 1)) * sk(j, k, i, 1));
                    wt(j, k, 8) = term * ((q(j, k, i, 4) - qk0(j, i, 4, 1)) * sk(j, k, i, 2));
                    wt(j, k, 9) = term * ((q(j, k, i, 4) - qk0(j, i, 4, 1)) * sk(j, k, i, 3));
                }
            }

            // k=kdim boundary
            {
                int k = kdim;
                int km2 = std::max(1, k - 2);
                for (int j = 1; j <= jdim1; j++) {
                    double term = sk(j, k, i, 4) / (volk0(j, i, 3) + vol(j, k-1, i));
                    double factor = bck(j, i, 2) + 1.0;
                    term = term * factor;
                    wt(j, k, 1) = term * ((qk0(j, i, 2, 3) - q(j, k-1, i, 2)) * sk(j, k, i, 1));
                    wt(j, k, 2) = term * ((qk0(j, i, 2, 3) - q(j, k-1, i, 2)) * sk(j, k, i, 2));
                    wt(j, k, 3) = term * ((qk0(j, i, 2, 3) - q(j, k-1, i, 2)) * sk(j, k, i, 3));
                    wt(j, k, 4) = term * ((qk0(j, i, 3, 3) - q(j, k-1, i, 3)) * sk(j, k, i, 1));
                    wt(j, k, 5) = term * ((qk0(j, i, 3, 3) - q(j, k-1, i, 3)) * sk(j, k, i, 2));
                    wt(j, k, 6) = term * ((qk0(j, i, 3, 3) - q(j, k-1, i, 3)) * sk(j, k, i, 3));
                    wt(j, k, 7) = term * ((qk0(j, i, 4, 3) - q(j, k-1, i, 4)) * sk(j, k, i, 1));
                    wt(j, k, 8) = term * ((qk0(j, i, 4, 3) - q(j, k-1, i, 4)) * sk(j, k, i, 2));
                    wt(j, k, 9) = term * ((qk0(j, i, 4, 3) - q(j, k-1, i, 4)) * sk(j, k, i, 3));
                }
            }

            // cycle through cell centers
            for (int l = 1; l <= 9; l++)
                for (int k = 1; k <= kdim1; k++)
                    for (int j = 1; j <= jdim1; j++)
                        ux(j, k, i, l) = ux(j, k, i, l) + wt(j, k, l) + wt(j, k+1, l);

        } // i
    } // ifk

    //
    // I-direction contributions
    //
    if (ifi > 0) {
        // cycle through interfaces
        if (i2d == 0) {
            for (int i = 2; i <= idim1; i++) {
                for (int k = 1; k <= kdim1; k++) {
                    for (int j = 1; j <= jdim1; j++) {
                        double term = si(j, k, i, 4) / (vol(j, k, i) + vol(j, k, i-1));
                        wt(j, k, 1) = term * ((q(j, k, i, 2) - q(j, k, i-1, 2)) * si(j, k, i, 1));
                        wt(j, k, 2) = term * ((q(j, k, i, 2) - q(j, k, i-1, 2)) * si(j, k, i, 2));
                        wt(j, k, 3) = term * ((q(j, k, i, 2) - q(j, k, i-1, 2)) * si(j, k, i, 3));
                        wt(j, k, 4) = term * ((q(j, k, i, 3) - q(j, k, i-1, 3)) * si(j, k, i, 1));
                        wt(j, k, 5) = term * ((q(j, k, i, 3) - q(j, k, i-1, 3)) * si(j, k, i, 2));
                        wt(j, k, 6) = term * ((q(j, k, i, 3) - q(j, k, i-1, 3)) * si(j, k, i, 3));
                        wt(j, k, 7) = term * ((q(j, k, i, 4) - q(j, k, i-1, 4)) * si(j, k, i, 1));
                        wt(j, k, 8) = term * ((q(j, k, i, 4) - q(j, k, i-1, 4)) * si(j, k, i, 2));
                        wt(j, k, 9) = term * ((q(j, k, i, 4) - q(j, k, i-1, 4)) * si(j, k, i, 3));

                        ux(j, k, i-1, 1) = ux(j, k, i-1, 1) + wt(j, k, 1);
                        ux(j, k, i-1, 2) = ux(j, k, i-1, 2) + wt(j, k, 2);
                        ux(j, k, i-1, 3) = ux(j, k, i-1, 3) + wt(j, k, 3);
                        ux(j, k, i-1, 4) = ux(j, k, i-1, 4) + wt(j, k, 4);
                        ux(j, k, i-1, 5) = ux(j, k, i-1, 5) + wt(j, k, 5);
                        ux(j, k, i-1, 6) = ux(j, k, i-1, 6) + wt(j, k, 6);
                        ux(j, k, i-1, 7) = ux(j, k, i-1, 7) + wt(j, k, 7);
                        ux(j, k, i-1, 8) = ux(j, k, i-1, 8) + wt(j, k, 8);
                        ux(j, k, i-1, 9) = ux(j, k, i-1, 9) + wt(j, k, 9);

                        ux(j, k, i, 1) = ux(j, k, i, 1) + wt(j, k, 1);
                        ux(j, k, i, 2) = ux(j, k, i, 2) + wt(j, k, 2);
                        ux(j, k, i, 3) = ux(j, k, i, 3) + wt(j, k, 3);
                        ux(j, k, i, 4) = ux(j, k, i, 4) + wt(j, k, 4);
                        ux(j, k, i, 5) = ux(j, k, i, 5) + wt(j, k, 5);
                        ux(j, k, i, 6) = ux(j, k, i, 6) + wt(j, k, 6);
                        ux(j, k, i, 7) = ux(j, k, i, 7) + wt(j, k, 7);
                        ux(j, k, i, 8) = ux(j, k, i, 8) + wt(j, k, 8);
                        ux(j, k, i, 9) = ux(j, k, i, 9) + wt(j, k, 9);
                    }
                }
            }

            // i=1 boundary (ii=1)
            {
                int ii = 1;
                int iip1 = std::min(idim - 1, ii + 1);
                for (int k = 1; k <= kdim1; k++) {
                    for (int j = 1; j <= jdim1; j++) {
                        double term = si(j, k, ii, 4) / (voli0(j, k, 1) + vol(j, k, ii));
                        double factor = bci(j, k, 1) + 1.0;
                        term = term * factor;
                        wt(j, k, 1) = term * ((q(j, k, ii, 2) - qi0(j, k, 2, 1)) * si(j, k, ii, 1));
                        wt(j, k, 2) = term * ((q(j, k, ii, 2) - qi0(j, k, 2, 1)) * si(j, k, ii, 2));
                        wt(j, k, 3) = term * ((q(j, k, ii, 2) - qi0(j, k, 2, 1)) * si(j, k, ii, 3));
                        wt(j, k, 4) = term * ((q(j, k, ii, 3) - qi0(j, k, 3, 1)) * si(j, k, ii, 1));
                        wt(j, k, 5) = term * ((q(j, k, ii, 3) - qi0(j, k, 3, 1)) * si(j, k, ii, 2));
                        wt(j, k, 6) = term * ((q(j, k, ii, 3) - qi0(j, k, 3, 1)) * si(j, k, ii, 3));
                        wt(j, k, 7) = term * ((q(j, k, ii, 4) - qi0(j, k, 4, 1)) * si(j, k, ii, 1));
                        wt(j, k, 8) = term * ((q(j, k, ii, 4) - qi0(j, k, 4, 1)) * si(j, k, ii, 2));
                        wt(j, k, 9) = term * ((q(j, k, ii, 4) - qi0(j, k, 4, 1)) * si(j, k, ii, 3));
                        ux(j, k, 1, 1) = ux(j, k, 1, 1) + wt(j, k, 1);
                        ux(j, k, 1, 2) = ux(j, k, 1, 2) + wt(j, k, 2);
                        ux(j, k, 1, 3) = ux(j, k, 1, 3) + wt(j, k, 3);
                        ux(j, k, 1, 4) = ux(j, k, 1, 4) + wt(j, k, 4);
                        ux(j, k, 1, 5) = ux(j, k, 1, 5) + wt(j, k, 5);
                        ux(j, k, 1, 6) = ux(j, k, 1, 6) + wt(j, k, 6);
                        ux(j, k, 1, 7) = ux(j, k, 1, 7) + wt(j, k, 7);
                        ux(j, k, 1, 8) = ux(j, k, 1, 8) + wt(j, k, 8);
                        ux(j, k, 1, 9) = ux(j, k, 1, 9) + wt(j, k, 9);
                    }
                }
            }

            // i=idim boundary (ii=idim)
            {
                int ii = idim;
                int iim2 = std::max(1, ii - 2);
                for (int k = 1; k <= kdim1; k++) {
                    for (int j = 1; j <= jdim1; j++) {
                        double term = si(j, k, ii, 4) / (voli0(j, k, 3) + vol(j, k, ii-1));
                        double factor = bci(j, k, 2) + 1.0;
                        term = term * factor;
                        wt(j, k, 1) = term * ((qi0(j, k, 2, 3) - q(j, k, ii-1, 2)) * si(j, k, ii, 1));
                        wt(j, k, 2) = term * ((qi0(j, k, 2, 3) - q(j, k, ii-1, 2)) * si(j, k, ii, 2));
                        wt(j, k, 3) = term * ((qi0(j, k, 2, 3) - q(j, k, ii-1, 2)) * si(j, k, ii, 3));
                        wt(j, k, 4) = term * ((qi0(j, k, 3, 3) - q(j, k, ii-1, 3)) * si(j, k, ii, 1));
                        wt(j, k, 5) = term * ((qi0(j, k, 3, 3) - q(j, k, ii-1, 3)) * si(j, k, ii, 2));
                        wt(j, k, 6) = term * ((qi0(j, k, 3, 3) - q(j, k, ii-1, 3)) * si(j, k, ii, 3));
                        wt(j, k, 7) = term * ((qi0(j, k, 4, 3) - q(j, k, ii-1, 4)) * si(j, k, ii, 1));
                        wt(j, k, 8) = term * ((qi0(j, k, 4, 3) - q(j, k, ii-1, 4)) * si(j, k, ii, 2));
                        wt(j, k, 9) = term * ((qi0(j, k, 4, 3) - q(j, k, ii-1, 4)) * si(j, k, ii, 3));
                        ux(j, k, idim1, 1) = ux(j, k, idim1, 1) + wt(j, k, 1);
                        ux(j, k, idim1, 2) = ux(j, k, idim1, 2) + wt(j, k, 2);
                        ux(j, k, idim1, 3) = ux(j, k, idim1, 3) + wt(j, k, 3);
                        ux(j, k, idim1, 4) = ux(j, k, idim1, 4) + wt(j, k, 4);
                        ux(j, k, idim1, 5) = ux(j, k, idim1, 5) + wt(j, k, 5);
                        ux(j, k, idim1, 6) = ux(j, k, idim1, 6) + wt(j, k, 6);
                        ux(j, k, idim1, 7) = ux(j, k, idim1, 7) + wt(j, k, 7);
                        ux(j, k, idim1, 8) = ux(j, k, idim1, 8) + wt(j, k, 8);
                        ux(j, k, idim1, 9) = ux(j, k, idim1, 9) + wt(j, k, 9);
                    }
                }
            }
        } // i2d == 0
    } // ifi

    //
    // set velocity derivatives for hole cells to one
    //
    if (iover == 1) {
        for (int nk = 1; nk <= 9; nk++) {
            for (int i = 1; i <= idim1; i++) {
                for (int k = 1; k <= kdim1; k++) {
                    for (int j = 1; j <= jdim1; j++) {
                        // ccvmgt(onec, ux(j,k,i,nk), (real(blank(j,k,i)).eq.0.e0))
                        // = if blank==0 then onec else ux
                        ux(j, k, i, nk) = ((float)blank(j, k, i) == 0.0f) ? onec : ux(j, k, i, nk);
                    }
                }
            }
        }
    }
}

} // namespace delv_ns
