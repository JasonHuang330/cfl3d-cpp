// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "arc.h"
#include "topol2.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace arc_ns {

void arc(int& jdim1, int& kdim1, int& msub1, int& msub2,
         FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1,
         int& lmax1,
         FortranArray3DRef<double> x1, FortranArray3DRef<double> y1,
         FortranArray3DRef<double> z1,
         int& limit0,
         FortranArray1DRef<int> jjmax2, FortranArray1DRef<int> kkmax2,
         FortranArray3DRef<double> x2, FortranArray3DRef<double> y2,
         FortranArray3DRef<double> z2,
         FortranArray1DRef<double> xie2, FortranArray1DRef<double> eta2,
         FortranArray1DRef<int> mblkpt,
         int& ifit, int& itmax, int& jcorr, int& kcorr,
         FortranArray3DRef<double> sxie, FortranArray3DRef<double> seta,
         FortranArray3DRef<double> sxie2, FortranArray3DRef<double> seta2,
         FortranArray2DRef<double> xie2s, FortranArray2DRef<double> eta2s,
         int& intmx, int& icheck,
         FortranArray1DRef<int> nblkj, FortranArray1DRef<int> nblkk,
         FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm,
         int& j21, int& j22, int& k21, int& k22, int& npt,
         FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2,
         FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2,
         FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
         int& nbuf, int& ibufdim,
         FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    int& myid = cmn_mydist2.myid;

    int l2, jmax2, kmax2;
    int jmax1, kmax1;
    int l, j, k;
    int km, kup, kup2;
    int ll, ll1, jp, lp, jmax, kmax, kp, lsrch;
    double test1, test2, sc, xie, eta, xiet, etat, phi1, phi2;
    int jm, jup, jup2;

    l2    = 1;
    jmax2 = jjmax2(l2);
    kmax2 = kkmax2(l2);

    // arc lengths for (unexpanded) "from" grid(s)
    for (l = 1; l <= lmax1; l++) {
        jmax1 = jjmax1(l);
        kmax1 = kkmax1(l);
        for (j = 2; j <= jmax1-1; j++) {
            seta(j-1, 1, l) = 0.0;
            for (k = 3; k <= kmax1-1; k++) {
                seta(j-1, k-1, l) = seta(j-1, k-2, l) +
                    std::sqrt((x1(j,k,l)-x1(j,k-1,l))*(x1(j,k,l)-x1(j,k-1,l)) +
                              (y1(j,k,l)-y1(j,k-1,l))*(y1(j,k,l)-y1(j,k-1,l)) +
                              (z1(j,k,l)-z1(j,k-1,l))*(z1(j,k,l)-z1(j,k-1,l)));
            }
        }
        for (k = 2; k <= kmax1-1; k++) {
            sxie(1, k-1, l) = 0.0;
            for (j = 3; j <= jmax1-1; j++) {
                sxie(j-1, k-1, l) = sxie(j-2, k-1, l) +
                    std::sqrt((x1(j,k,l)-x1(j-1,k,l))*(x1(j,k,l)-x1(j-1,k,l)) +
                              (y1(j,k,l)-y1(j-1,k,l))*(y1(j,k,l)-y1(j-1,k,l)) +
                              (z1(j,k,l)-z1(j-1,k,l))*(z1(j,k,l)-z1(j-1,k,l)));
            }
        }
    }

    // arc lengths for "to" grid
    for (j = 1; j <= jmax2; j++) {
        seta2(j, 1, 1) = 0.0;
        for (k = 2; k <= kmax2; k++) {
            seta2(j, k, 1) = seta2(j, k-1, 1) +
                std::sqrt((x2(j,k,1)-x2(j,k-1,1))*(x2(j,k,1)-x2(j,k-1,1)) +
                          (y2(j,k,1)-y2(j,k-1,1))*(y2(j,k,1)-y2(j,k-1,1)) +
                          (z2(j,k,1)-z2(j,k-1,1))*(z2(j,k,1)-z2(j,k-1,1)));
        }
    }
    for (k = 1; k <= kmax2; k++) {
        sxie2(1, k, 1) = 0.0;
        for (j = 2; j <= jmax2; j++) {
            sxie2(j, k, 1) = sxie2(j-1, k, 1) +
                std::sqrt((x2(j,k,1)-x2(j-1,k,1))*(x2(j,k,1)-x2(j-1,k,1)) +
                          (y2(j,k,1)-y2(j-1,k,1))*(y2(j,k,1)-y2(j-1,k,1)) +
                          (z2(j,k,1)-z2(j-1,k,1))*(z2(j,k,1)-z2(j-1,k,1)));
        }
    }

    if (kcorr == 1) {
        // try arc length correction near eta=0
        //
        // Arc length in eta direction now replaces the
        // physical coordinates x,y,z from which computational
        // coordinates away from the boundaries are to be found.
        // The current xie values are assumed correct.

        for (j = j21; j <= j22-1; j++) {
            km = kmm(j);
            if (km < k21) continue;

            test1 = 0.5;
            test2 = 0.75;
            kup  = (int)(test1 * kmm(j));
            kup2 = (int)(test2 * kmm(j));
            if (kup2 <= k21) continue;

            for (k = k21; k <= kup2; k++) {
                // center of "to" grid cell j,k
                sc = 0.25 * (seta2(j,k,1) + seta2(j+1,k,1) +
                             seta2(j+1,k+1,1) + seta2(j,k+1,1));

                ll = (j22-j21)*(k-k21) + (j-j21+1);
                jp   = (int)xie2(ll);
                xie  = xie2(ll) - jp;
                lp   = mblkpt(ll);
                jmax = jjmax1(lp);
                kmax = kkmax1(lp);
                if (k == k21) {
                    kp = 1;
                } else {
                    kp = (int)eta2s(j, k-1) + 1;
                }
                eta   = (double)kp;
                lsrch = 2;
                topol2_ns::topol2(jdim1, kdim1, msub1, jjmax1, kkmax1, lmax1,
                                  xie, eta, seta,
                                  limit0, sc, jp, kp, lp, lsrch, itmax,
                                  xiet, etat,
                                  xif1, xif2, etf1, etf2,
                                  nou, bou, nbuf, ibufdim,
                                  myid, mblk2nd, maxbl);
                eta2s(j, k) = etat;
            }

            // Blend arc length corrected and uncorrected computational
            // coordinates of "to" grid...arc length correction for k .le. kup,
            // uncorrected for k .ge. kup2, and linear interpolation of corrected
            // and uncorrected values in between.

            for (k = k21; k <= kup; k++) {
                ll = (j22-j21)*(k-k21) + (j-j21+1);
                eta2(ll) = eta2s(j, k);
            }
            if (kup2 > kup) {
                for (k = kup; k <= kup2; k++) {
                    phi1 = 1.0 - (float)(k-kup) / (float)(kup2-kup);
                    phi2 = 1.0 - phi1;
                    ll = (j22-j21)*(k-k21) + (j-j21+1);
                    eta2(ll) = eta2(ll)*phi2 + eta2s(j,k)*phi1;
                }
            }

            // check for monotonicity of corrected coordinates
            for (k = k21+1; k <= km; k++) {
                ll  = (j22-j21)*(k-k21)   + (j-j21+1);
                ll1 = (j22-j21)*(k-1-k21) + (j-j21+1);
                if ((float)eta2(ll) <= (float)eta2(ll1)) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120,
                        "  WARNING: corrected eta values are"
                        " not monotonic at j,k = %5d,%5d", j, k);
                }
            }
        }
    }

    if (jcorr == 1) {
        // try arc length correction near xie=0
        //
        // Arc length in xie direction now replaces the
        // physical coordinates x,y,z from which computational
        // coordinates away from the boundaries are to be found.
        // The current eta values are assumed correct.

        for (k = k21; k <= k22-1; k++) {
            jm = jmm(k);
            if (jm < j21) continue;

            test1 = 0.5;
            test2 = 0.75;
            jup  = (int)(test1 * jmm(k));
            jup2 = (int)(test2 * jmm(k));
            if (jup2 <= j21) continue;

            for (j = j21; j <= jup2; j++) {
                // center of "to" grid cell j,k
                sc = 0.25 * (sxie2(j,k,1) + sxie2(j+1,k,1) +
                             sxie2(j+1,k+1,1) + sxie2(j,k+1,1));

                ll = (j22-j21)*(k-k21) + (j-j21+1);
                kp   = (int)eta2(ll);
                eta  = eta2(ll) - kp;
                lp   = mblkpt(ll);
                jmax = jjmax1(lp);
                kmax = kkmax1(lp);
                if (j == j21) {
                    jp = 1;
                } else {
                    jp = (int)xie2s(j-1, k) + 1;
                }
                xie   = (double)jp;
                lsrch = 1;
                topol2_ns::topol2(jdim1, kdim1, msub1, jjmax1, kkmax1, lmax1,
                                  xie, eta, sxie,
                                  limit0, sc, jp, kp, lp, lsrch, itmax,
                                  xiet, etat,
                                  xif1, xif2, etf1, etf2,
                                  nou, bou, nbuf, ibufdim,
                                  myid, mblk2nd, maxbl);
                xie2s(j, k) = xiet;
            }

            // Blend arc length corrected and uncorrected computational
            // coordinates of "to" grid...arc length correction for j .le. jup,
            // uncorrected for j .ge. jup2, and linear interpolation of corrected
            // and uncorrected values in between.

            for (j = j21; j <= jup; j++) {
                ll = (j22-j21)*(k-k21) + (j-j21+1);
                xie2(ll) = xie2s(j, k);
            }
            for (j = jup; j <= jup2; j++) {
                phi1 = 1.0 - (float)(j-jup) / (float)(jup2-jup);
                phi2 = 1.0 - phi1;
                ll = (j22-j21)*(k-k21) + (j-j21+1);
                xie2(ll) = xie2(ll)*phi2 + xie2s(j,k)*phi1;
            }

            // check for monotonicity of corrected coordinates
            for (j = j21+1; j <= jm; j++) {
                ll  = (j22-j21)*(k-k21) + (j-j21+1);
                ll1 = (j22-j21)*(k-k21) + (j-1-j21+1);
                if ((float)xie2(ll) <= (float)xie2(ll1)) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120,
                        "  WARNING: corrected xie values are"
                        " not monotonic at j,k = %5d,%5d", j, k);
                }
            }
        }
    }

    return;
}

} // namespace arc_ns
