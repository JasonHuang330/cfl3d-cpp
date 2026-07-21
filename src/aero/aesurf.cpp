// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "aesurf.h"

namespace aesurf_ns {

void aesurf(int& nbl, int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> deltj,
            FortranArray4DRef<double> deltk,
            FortranArray4DRef<double> delti,
            FortranArray5DRef<double> xmdj,
            FortranArray5DRef<double> xmdk,
            FortranArray5DRef<double> xmdi,
            FortranArray3DRef<double> wkj,
            FortranArray3DRef<double> wkk,
            FortranArray3DRef<double> wki,
            int& maxbl, int& maxseg, int& nmds, int& maxaes,
            FortranArray2DRef<double> aesrfdat,
            FortranArray2DRef<double> xs,
            FortranArray2DRef<double> xxn,
            FortranArray2DRef<int> icsi,
            FortranArray2DRef<int> icsf,
            FortranArray2DRef<int> jcsi,
            FortranArray2DRef<int> jcsf,
            FortranArray2DRef<int> kcsi,
            FortranArray2DRef<int> kcsf,
            FortranArray1DRef<int> nsegdfrm,
            FortranArray2DRef<int> idfrmseg,
            FortranArray2DRef<int> iaesurf,
            int& maxsegdg)
{
    // common /elastic/ ndefrm, naesrf
    int& naesrf = cmn_elastic.naesrf;

    for (int iaes = 1; iaes <= naesrf; iaes++) {

        int    iskyhk  = (int)aesrfdat(1, iaes);
        double grefl   = aesrfdat(2, iaes);
        double uinf    = aesrfdat(3, iaes);
        double qinf    = aesrfdat(4, iaes);
        int    nmodes  = (int)aesrfdat(5, iaes);
        double rgrefl  = 1.0 / grefl;

        for (int is = 1; is <= nsegdfrm(nbl); is++) {

            if (iaesurf(nbl, is) == iaes && idfrmseg(nbl, is) == 99) {

                // convert modal deflections to x,y,z - j=const surfaces
                if (jcsi(nbl, is) == jcsf(nbl, is)) {

                    int j  = jcsi(nbl, is);
                    int ll = 0;
                    int m  = 1;
                    if (jcsi(nbl, is) == jdim) {
                        ll = 3;
                        m  = 2;
                    }
                    int ist = icsi(nbl, is);
                    int ifn = icsf(nbl, is);
                    int kst = kcsi(nbl, is);
                    int kfn = kcsf(nbl, is);

                    for (int i = ist; i <= ifn; i++) {
                        for (int k = kst; k <= kfn; k++) {

                            // save any existing displacements
                            double deltj1 = deltj(k, i, 1, m);
                            double deltj2 = deltj(k, i, 2, m);
                            double deltj3 = deltj(k, i, 3, m);
                            // initialize modal displacements
                            deltj(k, i, 1, m) = 0.0;
                            deltj(k, i, 2, m) = 0.0;
                            deltj(k, i, 3, m) = 0.0;
                            // sum delta contributions from all modes
                            for (int n = 1; n <= nmodes; n++) {
                                deltj(k, i, 1, m) = deltj(k, i, 1, m) +
                                    wkj(k, i, m) * (xs(2*n-1, iaes)
                                    - xxn(2*n-1, iaes)) * xmdj(k, i, ll+1, n, iaes) * rgrefl;
                                deltj(k, i, 2, m) = deltj(k, i, 2, m) +
                                    wkj(k, i, m) * (xs(2*n-1, iaes)
                                    - xxn(2*n-1, iaes)) * xmdj(k, i, ll+2, n, iaes) * rgrefl;
                                deltj(k, i, 3, m) = deltj(k, i, 3, m) +
                                    wkj(k, i, m) * (xs(2*n-1, iaes)
                                    - xxn(2*n-1, iaes)) * xmdj(k, i, ll+3, n, iaes) * rgrefl;
                            }
                            // add modal displacements to existing displacements
                            deltj(k, i, 1, m) = deltj(k, i, 1, m) + deltj1;
                            deltj(k, i, 2, m) = deltj(k, i, 2, m) + deltj2;
                            deltj(k, i, 3, m) = deltj(k, i, 3, m) + deltj3;
                            wkj(k, i, m) = 0.0;
                        }
                    }

                } // end j=const

                // convert modal deflections to x,y,z - k=const surfaces
                if (kcsi(nbl, is) == kcsf(nbl, is)) {

                    int k  = kcsi(nbl, is);
                    int ll = 0;
                    int m  = 1;
                    if (kcsi(nbl, is) == kdim) {
                        ll = 3;
                        m  = 2;
                    }
                    int ist = icsi(nbl, is);
                    int ifn = icsf(nbl, is);
                    int jst = jcsi(nbl, is);
                    int jfn = jcsf(nbl, is);

                    for (int i = ist; i <= ifn; i++) {
                        for (int j = jst; j <= jfn; j++) {

                            // save any existing displacements
                            double deltk1 = deltk(j, i, 1, m);
                            double deltk2 = deltk(j, i, 2, m);
                            double deltk3 = deltk(j, i, 3, m);
                            // initialize modal displacements
                            deltk(j, i, 1, m) = 0.0;
                            deltk(j, i, 2, m) = 0.0;
                            deltk(j, i, 3, m) = 0.0;
                            // sum delta contributions from all modes
                            for (int n = 1; n <= nmodes; n++) {
                                deltk(j, i, 1, m) = deltk(j, i, 1, m) +
                                    wkk(j, i, m) * (xs(2*n-1, iaes)
                                    - xxn(2*n-1, iaes)) * xmdk(j, i, ll+1, n, iaes) * rgrefl;
                                deltk(j, i, 2, m) = deltk(j, i, 2, m) +
                                    wkk(j, i, m) * (xs(2*n-1, iaes)
                                    - xxn(2*n-1, iaes)) * xmdk(j, i, ll+2, n, iaes) * rgrefl;
                                deltk(j, i, 3, m) = deltk(j, i, 3, m) +
                                    wkk(j, i, m) * (xs(2*n-1, iaes)
                                    - xxn(2*n-1, iaes)) * xmdk(j, i, ll+3, n, iaes) * rgrefl;
                            }
                            // add modal displacements to existing displacements
                            deltk(j, i, 1, m) = deltk(j, i, 1, m) + deltk1;
                            deltk(j, i, 2, m) = deltk(j, i, 2, m) + deltk2;
                            deltk(j, i, 3, m) = deltk(j, i, 3, m) + deltk3;
                            wkk(j, i, m) = 0.0;
                        }
                    }

                } // end k=const

                // convert modal deflections to x,y,z - i=const surfaces
                if (icsi(nbl, is) == icsf(nbl, is)) {

                    int i  = icsi(nbl, is);
                    int ll = 0;
                    int m  = 1;
                    if (icsi(nbl, is) == idim) {
                        ll = 3;
                        m  = 2;
                    }
                    int kst = kcsi(nbl, is);
                    int kfn = kcsf(nbl, is);
                    int jst = jcsi(nbl, is);
                    int jfn = jcsf(nbl, is);

                    for (int k = kst; k <= kfn; k++) {
                        for (int j = jst; j <= jfn; j++) {

                            // save any existing displacements
                            double delti1 = delti(j, k, 1, m);
                            double delti2 = delti(j, k, 2, m);
                            double delti3 = delti(j, k, 3, m);
                            // initialize modal displacements
                            delti(j, k, 1, m) = 0.0;
                            delti(j, k, 2, m) = 0.0;
                            delti(j, k, 3, m) = 0.0;
                            // sum delta contributions from all modes
                            for (int n = 1; n <= nmodes; n++) {
                                delti(j, k, 1, m) = delti(j, k, 1, m) +
                                    wki(j, k, m) * (xs(2*n-1, iaes)
                                    - xxn(2*n-1, iaes)) * xmdi(j, k, ll+1, n, iaes) * rgrefl;
                                delti(j, k, 2, m) = delti(j, k, 2, m) +
                                    wki(j, k, m) * (xs(2*n-1, iaes)
                                    - xxn(2*n-1, iaes)) * xmdi(j, k, ll+2, n, iaes) * rgrefl;
                                delti(j, k, 3, m) = delti(j, k, 3, m) +
                                    wki(j, k, m) * (xs(2*n-1, iaes)
                                    - xxn(2*n-1, iaes)) * xmdi(j, k, ll+3, n, iaes) * rgrefl;
                            }
                            // add modal displacements to existing displacements
                            delti(j, k, 1, m) = delti(j, k, 1, m) + delti1;
                            delti(j, k, 2, m) = delti(j, k, 2, m) + delti2;
                            delti(j, k, 3, m) = delti(j, k, 3, m) + delti3;
                            wki(j, k, m) = 0.0;
                        }
                    }

                } // end i=const

            } // end if iaesurf && idfrmseg

        } // end is loop

    } // end iaes loop

    return;
}

} // namespace aesurf_ns
