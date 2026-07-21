// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "int2.h"
#include <cstdio>
#include "rotateq.h"
#include "ccomplex.h"

namespace int2_ns {

void int2(int& iq, int& jdim1, int& kdim1, int& idim1, int& jmax2, int& kmax2,
          int& nsub1, int& jmax1, int& kmax1, int& l1, int& i1,
          FortranArray4DRef<double> q1, FortranArray3DRef<double> bc2,
          int& i2, FortranArray4DRef<double> q2,
          FortranArray1DRef<double> xie2, FortranArray1DRef<double> eta2,
          FortranArray2DRef<double> q1g, FortranArray2DRef<double> q2g,
          FortranArray1DRef<double> temp, FortranArray1DRef<int> nblkpt,
          int& intmax, int& icheck, int& mtype,
          FortranArray2DRef<int> iindex, int& ifo, int& ldim, int& npt,
          int& j21, int& j22, int& k21, int& k22,
          FortranArray4DRef<double> q1wk, double& dthtx, double& dthty,
          double& dthtz, int& lim_ptch)
{
    double eps, factfo, adthtx, adthty, adthtz;
    double q1jk, q1j1k, q1jk1, phi, xiec, etac;
    int i22, jd, kd, id, ista, iend, jsta, jend, ksta, kend;
    int l, k, j, ll, lc, jc, kc, nv, izz;

    eps = 1.0e-06;

    // if non-zero, ensure lim_ptch is 1 and no other value
    if (lim_ptch != 0) lim_ptch = 1;

    // load q into work array q1wk
    nv = jdim1 * kdim1 * idim1;
    for (l = 1; l <= ldim; l++) {
        for (izz = 1; izz <= nv; izz++) {
            q1wk(izz, 1, 1, l) = q1(izz, 1, 1, l);
        }
    }

    // use injection (e.g. piecewise constant interpolation) if
    // solution is first order (ifo=0)
    factfo = 0.;
    if (ifo > 0) factfo = 1.;

    if (i2 == 1 || i2 == 2) i22 = 1;
    if (i2 == 3 || i2 == 4) i22 = 2;

    // interpolation on an i=constant surface
    if (mtype == 1) {

        // rotate "from" q array, if needed - rotated values stored in q1wk
        // (never need to rotate scalar turbulence quantities - iq <= 0)
        if (iq > 0) {
            adthtx = ccomplex_ns::ccabs(dthtx);
            adthty = ccomplex_ns::ccabs(dthty);
            adthtz = ccomplex_ns::ccabs(dthtz);
            if (adthtx > 0. || adthty > 0. || adthtz > 0.) {
                jd   = jdim1;
                kd   = kdim1;
                id   = idim1;
                ista = i1;
                iend = i1;
                jsta = 1;
                jend = jdim1;
                ksta = 1;
                kend = kdim1;
                rotateq_ns::rotateq(jd, kd, id, q1, q1wk, ista, iend,
                                    jsta, jend, ksta, kend,
                                    dthtx, dthty, dthtz);
            }
        }

        for (l = 1; l <= ldim; l++) {

            // determine gradient values in xie
            for (k = 1; k <= kmax1 - 1; k++) {
                for (j = 1; j <= jmax1 - 2; j++) {
                    q1jk    = q1wk(j,   k, i1, l);
                    q1j1k   = q1wk(j+1, k, i1, l);
                    temp(j) = q1j1k - q1jk;
                }
                for (j = 2; j <= jmax1 - 2; j++) {
                    phi      = temp(j) * temp(j-1)
                             / (temp(j)*temp(j) + temp(j-1)*temp(j-1) + eps);
                    phi      = 0.5 + (double)lim_ptch * (phi - 0.5);
                    q1g(j,k) = (temp(j) + temp(j-1)) * phi;
                }
                q1g(1,       k) = temp(1);
                q1g(jmax1-1, k) = temp(jmax1-2);
            }

            // determine gradient values in eta
            for (j = 1; j <= jmax1 - 1; j++) {
                for (k = 1; k <= kmax1 - 2; k++) {
                    q1jk    = q1wk(j, k,   i1, l);
                    q1jk1   = q1wk(j, k+1, i1, l);
                    temp(k) = q1jk1 - q1jk;
                }
                for (k = 2; k <= kmax1 - 2; k++) {
                    phi      = temp(k) * temp(k-1)
                             / (temp(k)*temp(k) + temp(k-1)*temp(k-1) + eps);
                    phi      = 0.5 + (double)lim_ptch * (phi - 0.5);
                    q2g(j,k) = (temp(k) + temp(k-1)) * phi;
                }
                q2g(j, 1)       = temp(1);
                q2g(j, kmax1-1) = temp(kmax1-2);
            }
            // label 4800 (unreferenced, fall through)

            // determine interpolated values
            for (k = k21; k <= k22 - 1; k++) {
                for (j = j21; j <= j22 - 1; j++) {
                    ll = (j22 - j21) * (k - k21) + (j - j21 + 1);
                    lc = nblkpt(ll);
                    if (lc != l1) goto label2001;
                    jc = (int)(xie2(ll));
                    kc = (int)(eta2(ll));

                    // keep within bounds of "from" grid
                    jc = std::max(1, jc);
                    kc = std::max(1, kc);
                    jc = std::min(jc, jmax1 - 1);
                    kc = std::min(kc, kmax1 - 1);
                    {
                        double tmp1 = 1.0;
                        xie2(ll) = ccomplex_ns::ccmaxcr(xie2(ll), tmp1);
                        eta2(ll) = ccomplex_ns::ccmaxcr(eta2(ll), tmp1);
                        double tmp2 = (double)jmax1;
                        xie2(ll) = ccomplex_ns::ccmincr(xie2(ll), tmp2);
                        double tmp3 = (double)kmax1;
                        eta2(ll) = ccomplex_ns::ccmincr(eta2(ll), tmp3);
                    }

                    xiec = (double)jc + 0.5;
                    etac = (double)kc + 0.5;

                    bc2(j, k, i22)  = 0.0;
                    q2(j, k, l, i2) = q1wk(jc, kc, i1, l)
                                    + (q1g(jc, kc) * (xie2(ll) - xiec)
                                    +  q2g(jc, kc) * (eta2(ll) - etac)) * factfo;
                    label2001:;
                }
            }
        } // end do 3800

    // interpolation on a j=constant surface
    } else if (mtype == 2) {

        // rotate "from" q array, if needed - rotated values stored in q1wk
        // (never need to rotate scalar turbulence quantities - iq <= 0)
        if (iq > 0) {
            adthtx = ccomplex_ns::ccabs(dthtx);
            adthty = ccomplex_ns::ccabs(dthty);
            adthtz = ccomplex_ns::ccabs(dthtz);
            if (adthtx > 0. || adthty > 0. || adthtz > 0.) {
                jd   = jdim1;
                kd   = kdim1;
                id   = idim1;
                ista = 1;
                iend = idim1;
                jsta = i1;
                jend = i1;
                ksta = 1;
                kend = kdim1;
                rotateq_ns::rotateq(jd, kd, id, q1, q1wk, ista, iend,
                                    jsta, jend, ksta, kend,
                                    dthtx, dthty, dthtz);
            }
        }

        for (l = 1; l <= ldim; l++) {

            // determine gradient values in xie
            for (k = 1; k <= kmax1 - 1; k++) {
                for (j = 1; j <= jmax1 - 2; j++) {
                    q1jk    = q1wk(i1, j,   k, l);
                    q1j1k   = q1wk(i1, j+1, k, l);
                    temp(j) = q1j1k - q1jk;
                }
                for (j = 2; j <= jmax1 - 2; j++) {
                    phi      = temp(j) * temp(j-1)
                             / (temp(j)*temp(j) + temp(j-1)*temp(j-1) + eps);
                    phi      = 0.5 + (double)lim_ptch * (phi - 0.5);
                    q1g(j,k) = (temp(j) + temp(j-1)) * phi;
                }
                q1g(1,       k) = temp(1);
                q1g(jmax1-1, k) = temp(jmax1-2);
            }

            // determine gradient values in eta
            for (j = 1; j <= jmax1 - 1; j++) {
                for (k = 1; k <= kmax1 - 2; k++) {
                    q1jk    = q1wk(i1, j, k,   l);
                    q1jk1   = q1wk(i1, j, k+1, l);
                    temp(k) = q1jk1 - q1jk;
                }
                for (k = 2; k <= kmax1 - 2; k++) {
                    phi      = temp(k) * temp(k-1)
                             / (temp(k)*temp(k) + temp(k-1)*temp(k-1) + eps);
                    phi      = 0.5 + (double)lim_ptch * (phi - 0.5);
                    q2g(j,k) = (temp(k) + temp(k-1)) * phi;
                }
                q2g(j, 1)       = temp(1);
                q2g(j, kmax1-1) = temp(kmax1-2);
            }
            // label 4810 (unreferenced, fall through)

            // determine interpolated values
            for (k = k21; k <= k22 - 1; k++) {
                for (j = j21; j <= j22 - 1; j++) {
                    ll = (j22 - j21) * (k - k21) + (j - j21 + 1);
                    lc = nblkpt(ll);
                    if (lc != l1) goto label2011;
                    jc = (int)(xie2(ll));
                    kc = (int)(eta2(ll));

                    // keep within bounds of "from" grid
                    jc = std::max(1, jc);
                    kc = std::max(1, kc);
                    jc = std::min(jc, jmax1 - 1);
                    kc = std::min(kc, kmax1 - 1);
                    {
                        double tmp1 = 1.0;
                        xie2(ll) = ccomplex_ns::ccmaxcr(xie2(ll), tmp1);
                        eta2(ll) = ccomplex_ns::ccmaxcr(eta2(ll), tmp1);
                        double tmp2 = (double)jmax1;
                        xie2(ll) = ccomplex_ns::ccmincr(xie2(ll), tmp2);
                        double tmp3 = (double)kmax1;
                        eta2(ll) = ccomplex_ns::ccmincr(eta2(ll), tmp3);
                    }

                    xiec = (double)jc + 0.5;
                    etac = (double)kc + 0.5;

                    bc2(j, k, i22)  = 0.0;
                    q2(j, k, l, i2) = q1wk(i1, jc, kc, l)
                                    + (q1g(jc, kc) * (xie2(ll) - xiec)
                                    +  q2g(jc, kc) * (eta2(ll) - etac)) * factfo;
                    label2011:;
                }
            }
        } // end do 3810

    // interpolation on a k=constant surface
    } else if (mtype == 3) {

        // rotate "from" q array, if needed - rotated values stored in q1wk
        // (never need to rotate scalar turbulence quantities - iq <= 0)
        if (iq > 0) {
            adthtx = ccomplex_ns::ccabs(dthtx);
            adthty = ccomplex_ns::ccabs(dthty);
            adthtz = ccomplex_ns::ccabs(dthtz);
            if (adthtx > 0. || adthty > 0. || adthtz > 0.) {
                jd   = jdim1;
                kd   = kdim1;
                id   = idim1;
                ista = 1;
                iend = idim1;
                jsta = 1;
                jend = jdim1;
                ksta = i1;
                kend = i1;
                rotateq_ns::rotateq(jd, kd, id, q1, q1wk, ista, iend,
                                    jsta, jend, ksta, kend,
                                    dthtx, dthty, dthtz);
            }
        }

        for (l = 1; l <= ldim; l++) {

            // determine gradient values in xie
            for (k = 1; k <= kmax1 - 1; k++) {
                for (j = 1; j <= jmax1 - 2; j++) {
                    q1jk    = q1wk(j,   i1, k, l);
                    q1j1k   = q1wk(j+1, i1, k, l);
                    temp(j) = q1j1k - q1jk;
                }
                for (j = 2; j <= jmax1 - 2; j++) {
                    phi      = temp(j) * temp(j-1)
                             / (temp(j)*temp(j) + temp(j-1)*temp(j-1) + eps);
                    phi      = 0.5 + (double)lim_ptch * (phi - 0.5);
                    q1g(j,k) = (temp(j) + temp(j-1)) * phi;
                }
                q1g(1,       k) = temp(1);
                q1g(jmax1-1, k) = temp(jmax1-2);
            }

            // determine gradient values in eta
            for (j = 1; j <= jmax1 - 1; j++) {
                for (k = 1; k <= kmax1 - 2; k++) {
                    q1jk    = q1wk(j, i1, k,   l);
                    q1jk1   = q1wk(j, i1, k+1, l);
                    temp(k) = q1jk1 - q1jk;
                }
                for (k = 2; k <= kmax1 - 2; k++) {
                    phi      = temp(k) * temp(k-1)
                             / (temp(k)*temp(k) + temp(k-1)*temp(k-1) + eps);
                    phi      = 0.5 + (double)lim_ptch * (phi - 0.5);
                    q2g(j,k) = (temp(k) + temp(k-1)) * phi;
                }
                q2g(j, 1)       = temp(1);
                q2g(j, kmax1-1) = temp(kmax1-2);
            }
            // label 4820 (unreferenced, fall through)

            // determine interpolated values
            for (k = k21; k <= k22 - 1; k++) {
                for (j = j21; j <= j22 - 1; j++) {
                    ll = (j22 - j21) * (k - k21) + (j - j21 + 1);
                    lc = nblkpt(ll);
                    if (lc != l1) goto label2021;
                    jc = (int)(xie2(ll));
                    kc = (int)(eta2(ll));

                    // keep within bounds of "from" grid
                    jc = std::max(1, jc);
                    kc = std::max(1, kc);
                    jc = std::min(jc, jmax1 - 1);
                    kc = std::min(kc, kmax1 - 1);
                    {
                        double tmp1 = 1.0;
                        xie2(ll) = ccomplex_ns::ccmaxcr(xie2(ll), tmp1);
                        eta2(ll) = ccomplex_ns::ccmaxcr(eta2(ll), tmp1);
                        double tmp2 = (double)jmax1;
                        xie2(ll) = ccomplex_ns::ccmincr(xie2(ll), tmp2);
                        double tmp3 = (double)kmax1;
                        eta2(ll) = ccomplex_ns::ccmincr(eta2(ll), tmp3);
                    }

                    xiec = (double)jc + 0.5;
                    etac = (double)kc + 0.5;

                    bc2(j, k, i22)  = 0.0;
                    q2(j, k, l, i2) = q1wk(jc, i1, kc, l)
                                    + (q1g(jc, kc) * (xie2(ll) - xiec)
                                    +  q2g(jc, kc) * (eta2(ll) - etac)) * factfo;
                    label2021:;
                }
            }
        } // end do 3820

    } // end if mtype

    return;
}

} // namespace int2_ns
