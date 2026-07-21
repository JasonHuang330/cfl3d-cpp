// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "metric.h"
#include "termn8.h"
#include "readkey.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace metric_ns {

void metric(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> sj, FortranArray3DRef<double> sk, FortranArray3DRef<double> si, FortranArray2DRef<double> t, FortranArray3DRef<double> t1, int& nbl, int& iflag, int& icnt, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& maxbl, int& maxseg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd)
{
    // local variables
    int isklton = cmn_sklton.isklton;
    int i2d     = cmn_twod.i2d;
    float atol  = cmn_singular.atol;
    int iexp    = cmn_zero.iexp;

    int icntmax = 10;
    double atol2 = std::pow(10.0, (double)(-iexp + 1));

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int n, ns, ne, nseg;
    int i, j, k, m, izz, izz1, izzp1, izzndx;
    int js, je, ks, ke, is, ie, mtyp;
    double asum, temp;
    double delta1, delta2, delta3;


    // *********************************************************************
    // metrics for i=constant surfaces
    // *********************************************************************
    //
    // *** interior faces ***
    //
    n = jdim * kdim - jdim;

    for (i = 2; i <= idim1; i++) {
        // cdir$ ivdep
        for (izz = 1; izz <= n; izz++) {
            // components of vectors connecting opposite corners of cell j,k
            t(izz, 1) = x(izz + 1, 1, i) - x(izz, 2, i);
            t(izz, 2) = y(izz + 1, 1, i) - y(izz, 2, i);
            t(izz, 3) = z(izz + 1, 1, i) - z(izz, 2, i);
            t(izz, 4) = x(izz + 1, 2, i) - x(izz, 1, i);
            t(izz, 5) = y(izz + 1, 2, i) - y(izz, 1, i);
            t(izz, 6) = z(izz + 1, 2, i) - z(izz, 1, i);
            // cross product of vectors
            si(izz, i, 1) =  t(izz, 2) * t(izz, 6) - t(izz, 3) * t(izz, 5);
            si(izz, i, 2) = -t(izz, 1) * t(izz, 6) + t(izz, 3) * t(izz, 4);
            si(izz, i, 3) =  t(izz, 1) * t(izz, 5) - t(izz, 2) * t(izz, 4);
            // magnitude of cross product
            si(izz, i, 4) = si(izz, i, 1) * si(izz, i, 1) + si(izz, i, 2) * si(izz, i, 2) +
                            si(izz, i, 3) * si(izz, i, 3);
        }
        // due to the ordering in izz, the cross product above is incorrect
        // at the fictitious interfaces izz=k*jdim (k=1...kdim-1), and may
        // have a zero value. for safety in the loop below, first set these
        // fictitious interface values to temporary safe values.
        for (k = 1; k <= kdim1; k++) {
            izz = k * jdim;
            si(izz, i, 1) = 0.;
            si(izz, i, 2) = 0.;
            si(izz, i, 3) = 0.;
            si(izz, i, 4) = 4.;
        }
        // cdir$ ivdep
        // store metrics:
        for (izz = 1; izz <= n; izz++) {
            si(izz, i, 4) = 1.e0 / std::sqrt(si(izz, i, 4));
            si(izz, i, 1) = si(izz, i, 1) * si(izz, i, 4);
            si(izz, i, 2) = si(izz, i, 2) * si(izz, i, 4);
            si(izz, i, 3) = si(izz, i, 3) * si(izz, i, 4);
            si(izz, i, 4) = 0.5e0 / si(izz, i, 4);
            si(izz, i, 5) = 0.0;
        }
    } // end do 1040

    // *** i=1/idim faces ***
    for (m = 1; m <= 2; m++) {
        if (m == 1) {
            i = 1;
            nseg = nbci0(nbl);
        } else {
            i = idim;
            nseg = nbcidim(nbl);
        }

        for (int ns_loop = 1; ns_loop <= nseg; ns_loop++) {
            js   = ibcinfo(nbl, ns_loop, 2, m);
            je   = ibcinfo(nbl, ns_loop, 3, m) - 1;
            ks   = ibcinfo(nbl, ns_loop, 4, m);
            ke   = ibcinfo(nbl, ns_loop, 5, m) - 1;
            mtyp = ibcinfo(nbl, ns_loop, 1, m);

            asum = 0.;

            for (j = js; j <= je; j++) {
                for (k = ks; k <= ke; k++) {
                    izz = (k - 1) * jdim + j;
                    // components of vectors connecting opposite corners of cell j,k
                    t(izz, 1) = x(j + 1, k,     i) - x(j, k + 1, i);
                    t(izz, 2) = y(j + 1, k,     i) - y(j, k + 1, i);
                    t(izz, 3) = z(j + 1, k,     i) - z(j, k + 1, i);
                    t(izz, 4) = x(j + 1, k + 1, i) - x(j, k,     i);
                    t(izz, 5) = y(j + 1, k + 1, i) - y(j, k,     i);
                    t(izz, 6) = z(j + 1, k + 1, i) - z(j, k,     i);
                    // cross product of vectors
                    si(izz, i, 1) =  t(izz, 2) * t(izz, 6) - t(izz, 3) * t(izz, 5);
                    si(izz, i, 2) = -t(izz, 1) * t(izz, 6) + t(izz, 3) * t(izz, 4);
                    si(izz, i, 3) =  t(izz, 1) * t(izz, 5) - t(izz, 2) * t(izz, 4);
                    // magnitude of cross product
                    si(izz, i, 4) = si(izz, i, 1) * si(izz, i, 1) + si(izz, i, 2) * si(izz, i, 2) +
                                    si(izz, i, 3) * si(izz, i, 3);
                    asum = asum + std::sqrt(si(izz, i, 4));
                    if (si(izz, i, 4) == 0.) si(izz, i, 4) = atol2;
                }
            }

            asum = 0.5 * asum;
            if ((float)asum < (float)atol) {
                // collapsed metrics
                if (isklton > 0) {
                    if (m == 1) {
                        if (isklton == 1) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   singular metrics: i=1     j=%3d,%3d  k=%3d,%3d",
                                js, je + 1, ks, ke + 1);
                        } else {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "      singular metrics: i=1     j=%3d,%3d  k=%3d,%3d",
                                js, je + 1, ks, ke + 1);
                        }
                    } else {
                        if (isklton == 1) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   singular metrics: i=idim  j=%3d,%3d  k=%3d,%3d",
                                js, je + 1, ks, ke + 1);
                        } else {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "      singular metrics: i=idim  j=%3d,%3d  k=%3d,%3d",
                                js, je + 1, ks, ke + 1);
                        }
                    }
                }
                izzndx = i - (2 * m) + 3;
                for (j = js; j <= je; j++) {
                    for (k = ks; k <= ke; k++) {
                        izz = (k - 1) * jdim + j;
                        // set directions on collapsed face equal to those at closest interior face
                        si(izz, i, 1) = si(izz, izzndx, 1);
                        si(izz, i, 2) = si(izz, izzndx, 2);
                        si(izz, i, 3) = si(izz, izzndx, 3);
                        si(izz, i, 4) = 0.e0;
                        si(izz, i, 5) = si(izz, izzndx, 5);
                    }
                }
            } else {
                // non-singular metrics
                for (j = js; j <= je; j++) {
                    for (k = ks; k <= ke; k++) {
                        izz = (k - 1) * jdim + j;
                        si(izz, i, 4) = 1.e0 / std::sqrt(si(izz, i, 4));
                        si(izz, i, 1) = si(izz, i, 1) * si(izz, i, 4);
                        si(izz, i, 2) = si(izz, i, 2) * si(izz, i, 4);
                        si(izz, i, 3) = si(izz, i, 3) * si(izz, i, 4);
                        si(izz, i, 4) = 0.5e0 / si(izz, i, 4);
                        si(izz, i, 5) = 0.0;
                    }
                }
            }
        } // end do 1010
    } // end do 1000


    // fill in extra values of si for safety
    for (i = 1; i <= idim; i++) {
        // set metrics at jdim
        for (k = 1; k <= kdim - 1; k++) {
            izz  = (k - 1) * jdim + jdim;
            izz1 = izz - 1;
            si(izz, i, 1) = si(izz1, i, 1);
            si(izz, i, 2) = si(izz1, i, 2);
            si(izz, i, 3) = si(izz1, i, 3);
            si(izz, i, 4) = si(izz1, i, 4);
            si(izz, i, 5) = si(izz1, i, 5);
        }
        // set metrics at kdim
        // cdir$ ivdep
        for (j = 1; j <= jdim; j++) {
            izz  = jdim * (kdim - 1) + j;
            izz1 = izz - jdim;
            si(izz, i, 1) = si(izz1, i, 1);
            si(izz, i, 2) = si(izz1, i, 2);
            si(izz, i, 3) = si(izz1, i, 3);
            si(izz, i, 4) = si(izz1, i, 4);
            si(izz, i, 5) = si(izz1, i, 5);
        }
    } // end do 1070

    // Don't do the following checks if iflag=-1:
    if (iflag != -1) {
        // Check for non-planar i-planes when 2-D:
        if (i2d == 1) {
            int izz1_2d = jdim * kdim;
            for (i = 1; i <= idim; i++) {
                for (m = 1; m <= 3; m++) {
                    temp = si(1, i, m);
                    for (izz = 2; izz <= izz1_2d; izz++) {
                        if (std::abs((float)(si(izz, i, m) - temp)) > 1.e-5f) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                " Error. Grid not planar at i,nbl=%5d%5d", i, nbl);
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                " The i-planes need to be planar for 2-D (i2d=1)");
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   ... or else the problem may be crossing lines (negative volumes) in the grid");
                            int m1 = -1;
                            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                        }
                    }
                }
            }
            for (j = 1; j <= jdim; j++) {
                for (k = 1; k <= kdim; k++) {
                    if (std::abs(x(j, k, 2) - x(j, k, 1)) > 1.e-5 ||
                        std::abs(z(j, k, 2) - z(j, k, 1)) > 1.e-5) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120, " ");
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " The 2 i-planes need to be identical (in their planar values) for 2-D (i2d=1)");
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            "   ... or else the problem may be that IALPH is set incorrectly for your grid");
                        int m1 = -1;
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }
        // Check for large changes in metric from one i to the next,
        // indicating problem with the grid
        for (k = 1; k <= kdim - 1; k++) {
            for (j = 1; j <= jdim - 1; j++) {
                izz = (k - 1) * jdim + j;
                for (i = 1; i <= idim - 1; i++) {
                    delta1 = std::abs(si(izz, i + 1, 1) - si(izz, i, 1));
                    delta2 = std::abs(si(izz, i + 1, 2) - si(izz, i, 2));
                    delta3 = std::abs(si(izz, i + 1, 3) - si(izz, i, 3));
                    if ((float)delta1 > 1.99f || (float)delta2 > 1.99f || (float)delta3 > 1.99f) {
                        iflag = 1;
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            "   FATAL si grid normal direction change near j,k,i,i+1=");
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            "     ... suspect bad grid%5d%5d%5d%5d", j, k, i, i + 1);
                    } else if ((float)delta1 > 1.5f || (float)delta2 > 1.5f || (float)delta3 > 1.5f) {
                        icnt = icnt + 1;
                        if (icnt <= icntmax) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   WARNING: Dramatic si grid norm direction change (>120deg)");
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "     near j,k,i,i+1=%5d%5d%5d%5d", j, k, i, i + 1);
                        } else if (icnt == icntmax + 1) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   NOTE: no further WARNINGS will be output...better check this grid!");
                        }
                    }
                }
            }
        }
    } // end if iflag != -1


    // *********************************************************************
    // metrics for j=constant surfaces
    // *********************************************************************
    //
    // *** interior faces ***
    //
    // note: interior metrics for j=constant surfaces are initially
    //       calculated with a reordered izz index (compared to si and
    //       sk metrics) to permit vectorization and yet skip over j=1/jdim,
    //       which may have collapsed metrics. after the metrics are
    //       evaluated, the proper order is reset.
    //
    for (i = 1; i <= idim; i++) {
        for (j = 1; j <= jdim; j++) {
            for (k = 1; k <= kdim; k++) {
                izz = (j - 1) * kdim + k;
                t1(izz, i, 1) = x(j, k, i);
                t1(izz, i, 2) = y(j, k, i);
                t1(izz, i, 3) = z(j, k, i);
            }
        }
    } // end do 2035

    ns = kdim + 1;
    ne = (jdim - 1) * kdim;

    for (i = 1; i <= idim1; i++) {
        // cdir$ ivdep
        for (izz = ns; izz <= ne; izz++) {
            // components of vectors connecting opposite corners of cell i,k
            t(izz, 1) = t1(izz + 1, i,     1) - t1(izz, i + 1, 1);
            t(izz, 2) = t1(izz + 1, i,     2) - t1(izz, i + 1, 2);
            t(izz, 3) = t1(izz + 1, i,     3) - t1(izz, i + 1, 3);
            t(izz, 4) = t1(izz + 1, i + 1, 1) - t1(izz, i,     1);
            t(izz, 5) = t1(izz + 1, i + 1, 2) - t1(izz, i,     2);
            t(izz, 6) = t1(izz + 1, i + 1, 3) - t1(izz, i,     3);
            // cross product of vectors
            sj(izz, i, 1) =  t(izz, 2) * t(izz, 6) - t(izz, 3) * t(izz, 5);
            sj(izz, i, 2) = -t(izz, 1) * t(izz, 6) + t(izz, 3) * t(izz, 4);
            sj(izz, i, 3) =  t(izz, 1) * t(izz, 5) - t(izz, 2) * t(izz, 4);
            // magnitude of cross product
            sj(izz, i, 4) = sj(izz, i, 1) * sj(izz, i, 1) + sj(izz, i, 2) * sj(izz, i, 2) +
                            sj(izz, i, 3) * sj(izz, i, 3);
        }
        // due to the ordering in izz, the cross product above is incorrect
        // at the fictitious interfaces izz=j*kdim (j=1...jdim-1), and may
        // have a zero value. for safety in the loop below, first set these
        // fictitious interface values to temporary safe values.
        for (j = 1; j <= jdim1; j++) {
            izz = kdim * j;
            sj(izz, i, 1) = 0.;
            sj(izz, i, 2) = 0.;
            sj(izz, i, 3) = 0.;
            sj(izz, i, 4) = 4.;
        }
        // store temporary metrics (these are not in correct order):
        for (izz = ns; izz <= ne; izz++) {
            t1(izz, i, 4) = 1.e0 / std::sqrt(sj(izz, i, 4));
            t1(izz, i, 1) = sj(izz, i, 1) * t1(izz, i, 4);
            t1(izz, i, 2) = sj(izz, i, 2) * t1(izz, i, 4);
            t1(izz, i, 3) = sj(izz, i, 3) * t1(izz, i, 4);
            t1(izz, i, 4) = 0.5e0 / t1(izz, i, 4);
            t1(izz, i, 5) = 0.0;
        }
        // store metrics in correct order
        for (j = 2; j <= jdim1; j++) {
            for (k = 1; k <= kdim; k++) {
                izz1 = (j - 1) * kdim + k;
                izz  = (k - 1) * jdim + j;
                sj(izz, i, 1) = t1(izz1, i, 1);
                sj(izz, i, 2) = t1(izz1, i, 2);
                sj(izz, i, 3) = t1(izz1, i, 3);
                sj(izz, i, 4) = t1(izz1, i, 4);
                sj(izz, i, 5) = t1(izz1, i, 5);
            }
        }
    } // end do 2040

    // *** j=1/jdim faces ***
    for (m = 1; m <= 2; m++) {
        if (m == 1) {
            j = 1;
            nseg = nbcj0(nbl);
        } else {
            j = jdim;
            nseg = nbcjdim(nbl);
        }

        for (int ns_loop = 1; ns_loop <= nseg; ns_loop++) {
            is   = jbcinfo(nbl, ns_loop, 2, m);
            ie   = jbcinfo(nbl, ns_loop, 3, m) - 1;
            ks   = jbcinfo(nbl, ns_loop, 4, m);
            ke   = jbcinfo(nbl, ns_loop, 5, m) - 1;
            mtyp = jbcinfo(nbl, ns_loop, 1, m);

            asum = 0.;

            for (i = is; i <= ie; i++) {
                for (k = ks; k <= ke; k++) {
                    izz = (k - 1) * jdim + j;
                    // components of vectors connecting opposite corners of cell i,k
                    t(izz, 1) = x(j, k + 1, i)     - x(j, k,     i + 1);
                    t(izz, 2) = y(j, k + 1, i)     - y(j, k,     i + 1);
                    t(izz, 3) = z(j, k + 1, i)     - z(j, k,     i + 1);
                    t(izz, 4) = x(j, k + 1, i + 1) - x(j, k,     i);
                    t(izz, 5) = y(j, k + 1, i + 1) - y(j, k,     i);
                    t(izz, 6) = z(j, k + 1, i + 1) - z(j, k,     i);
                    // cross product of vectors
                    sj(izz, i, 1) =  t(izz, 2) * t(izz, 6) - t(izz, 3) * t(izz, 5);
                    sj(izz, i, 2) = -t(izz, 1) * t(izz, 6) + t(izz, 3) * t(izz, 4);
                    sj(izz, i, 3) =  t(izz, 1) * t(izz, 5) - t(izz, 2) * t(izz, 4);
                    // magnitude of cross product
                    sj(izz, i, 4) = sj(izz, i, 1) * sj(izz, i, 1) + sj(izz, i, 2) * sj(izz, i, 2) +
                                    sj(izz, i, 3) * sj(izz, i, 3);
                    asum = asum + std::sqrt(sj(izz, i, 4));
                    if (sj(izz, i, 4) == 0.) sj(izz, i, 4) = atol2;
                }
            }

            asum = 0.5 * asum;
            if ((float)asum < (float)atol) {
                // collapsed metrics
                if (isklton > 0) {
                    if (m == 1) {
                        if (isklton == 1) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   singular metrics: j=1     i=%3d,%3d  k=%3d,%3d",
                                is, ie + 1, ks, ke + 1);
                        } else {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "      singular metrics: j=1     i=%3d,%3d  k=%3d,%3d",
                                is, ie + 1, ks, ke + 1);
                        }
                    } else {
                        if (isklton == 1) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   singular metrics: j=jdim  i=%3d,%3d  k=%3d,%3d",
                                is, ie + 1, ks, ke + 1);
                        } else {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "      singular metrics: j=jdim  i=%3d,%3d  k=%3d,%3d",
                                is, ie + 1, ks, ke + 1);
                        }
                    }
                }
                for (i = is; i <= ie; i++) {
                    for (k = ks; k <= ke; k++) {
                        izz    = (k - 1) * jdim + j;
                        izzndx = izz - (2 * m) + 3;
                        // set directions on collapsed face equal to those at closest interior face
                        sj(izz, i, 1) = sj(izzndx, i, 1);
                        sj(izz, i, 2) = sj(izzndx, i, 2);
                        sj(izz, i, 3) = sj(izzndx, i, 3);
                        sj(izz, i, 4) = 0.e0;
                        sj(izz, i, 5) = sj(izzndx, i, 5);
                    }
                }
            } else {
                // non-singular metrics
                for (i = is; i <= ie; i++) {
                    for (k = ks; k <= ke; k++) {
                        izz = (k - 1) * jdim + j;
                        sj(izz, i, 4) = 1.e0 / std::sqrt(sj(izz, i, 4));
                        sj(izz, i, 1) = sj(izz, i, 1) * sj(izz, i, 4);
                        sj(izz, i, 2) = sj(izz, i, 2) * sj(izz, i, 4);
                        sj(izz, i, 3) = sj(izz, i, 3) * sj(izz, i, 4);
                        sj(izz, i, 4) = 0.5e0 / sj(izz, i, 4);
                        sj(izz, i, 5) = 0.0;
                    }
                }
            }
        } // end do 2010
    } // end do 2000


    // fill in extra values of sj for safety
    for (i = 1; i <= idim - 1; i++) {
        // set metrics at kdim
        // cdir$ ivdep
        for (j = 1; j <= jdim; j++) {
            izz  = jdim * (kdim - 1) + j;
            izz1 = izz - jdim;
            sj(izz, i, 1) = sj(izz1, i, 1);
            sj(izz, i, 2) = sj(izz1, i, 2);
            sj(izz, i, 3) = sj(izz1, i, 3);
            sj(izz, i, 4) = sj(izz1, i, 4);
            sj(izz, i, 5) = sj(izz1, i, 5);
        }
    } // end do 2070

    // Don't do the following checks if iflag=-1:
    if (iflag != -1) {
        // Check for large changes in metric from one j to the next,
        // indicating problem with the grid
        for (i = 1; i <= idim - 1; i++) {
            for (k = 1; k <= kdim - 1; k++) {
                for (j = 1; j <= jdim - 1; j++) {
                    izz   = (k - 1) * jdim + j;
                    izzp1 = (k - 1) * jdim + j + 1;
                    delta1 = std::abs(sj(izzp1, i, 1) - sj(izz, i, 1));
                    delta2 = std::abs(sj(izzp1, i, 2) - sj(izz, i, 2));
                    delta3 = std::abs(sj(izzp1, i, 3) - sj(izz, i, 3));
                    if ((float)delta1 > 1.99f || (float)delta2 > 1.99f || (float)delta3 > 1.99f) {
                        iflag = 1;
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            "   FATAL sj grid normal direction change near j,j+1,k,i=");
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            "     ... suspect bad grid%5d%5d%5d%5d", j, j + 1, k, i);
                    } else if ((float)delta1 > 1.5f || (float)delta2 > 1.5f || (float)delta3 > 1.5f) {
                        icnt = icnt + 1;
                        if (icnt <= icntmax) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   WARNING: Dramatic sj grid norm direction change (>120deg)");
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "     near j,j+1,k,i=%5d%5d%5d%5d", j, j + 1, k, i);
                        } else if (icnt == icntmax + 1) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   NOTE: no further WARNINGS will be output...better check this grid!");
                        }
                    }
                }
            }
        }
    } // end if iflag != -1


    // *********************************************************************
    // metrics for k=constant surfaces
    // *********************************************************************
    //
    // *** interior faces ***
    //
    ns = jdim + 1;
    ne = jdim * (kdim - 1);

    for (i = 1; i <= idim1; i++) {
        // cdir$ ivdep
        for (izz = ns; izz <= ne; izz++) {
            // components of vectors connecting opposite corners of cell i,j
            t(izz, 1) = x(izz + 1, 1, i)     - x(izz, 1, i + 1);
            t(izz, 2) = y(izz + 1, 1, i)     - y(izz, 1, i + 1);
            t(izz, 3) = z(izz + 1, 1, i)     - z(izz, 1, i + 1);
            t(izz, 4) = x(izz,     1, i)     - x(izz + 1, 1, i + 1);
            t(izz, 5) = y(izz,     1, i)     - y(izz + 1, 1, i + 1);
            t(izz, 6) = z(izz,     1, i)     - z(izz + 1, 1, i + 1);
            // cross product of vectors
            sk(izz, i, 1) =  t(izz, 2) * t(izz, 6) - t(izz, 3) * t(izz, 5);
            sk(izz, i, 2) = -t(izz, 1) * t(izz, 6) + t(izz, 3) * t(izz, 4);
            sk(izz, i, 3) =  t(izz, 1) * t(izz, 5) - t(izz, 2) * t(izz, 4);
            // magnitude of cross product
            sk(izz, i, 4) = sk(izz, i, 1) * sk(izz, i, 1) + sk(izz, i, 2) * sk(izz, i, 2) +
                            sk(izz, i, 3) * sk(izz, i, 3);
        }
        // due to the ordering in izz, the cross product above is incorrect
        // at the fictitious interfaces izz=j*kdim (j=1...jdim-1), and may
        // have a zero value. for safety in the loop below, first set these
        // fictitious interface values to temporary safe values.
        for (k = 1; k <= kdim; k++) {
            izz = jdim * k;
            sk(izz, i, 1) = 0.e0;
            sk(izz, i, 2) = 0.e0;
            sk(izz, i, 3) = 0.e0;
            sk(izz, i, 4) = 4.e0;
        }
        // store metrics:
        for (izz = ns; izz <= ne; izz++) {
            sk(izz, i, 4) = 1.e0 / std::sqrt(sk(izz, i, 4));
            sk(izz, i, 1) = sk(izz, i, 1) * sk(izz, i, 4);
            sk(izz, i, 2) = sk(izz, i, 2) * sk(izz, i, 4);
            sk(izz, i, 3) = sk(izz, i, 3) * sk(izz, i, 4);
            sk(izz, i, 4) = 0.5e0 / sk(izz, i, 4);
            sk(izz, i, 5) = 0.0;
        }
    } // end do 3040

    // *** k=1/kdim faces ***
    for (m = 1; m <= 2; m++) {
        if (m == 1) {
            k = 1;
            nseg = nbck0(nbl);
        } else {
            k = kdim;
            nseg = nbckdim(nbl);
        }

        for (int ns_loop = 1; ns_loop <= nseg; ns_loop++) {
            is   = kbcinfo(nbl, ns_loop, 2, m);
            ie   = kbcinfo(nbl, ns_loop, 3, m) - 1;
            js   = kbcinfo(nbl, ns_loop, 4, m);
            je   = kbcinfo(nbl, ns_loop, 5, m) - 1;
            mtyp = kbcinfo(nbl, ns_loop, 1, m);

            asum = 0.;

            for (i = is; i <= ie; i++) {
                for (j = js; j <= je; j++) {
                    izz = (k - 1) * jdim + j;
                    // components of vectors connecting opposite corners of cell i,j
                    t(izz, 1) = x(j + 1, k, i)     - x(j, k, i + 1);
                    t(izz, 2) = y(j + 1, k, i)     - y(j, k, i + 1);
                    t(izz, 3) = z(j + 1, k, i)     - z(j, k, i + 1);
                    t(izz, 4) = x(j,     k, i)     - x(j + 1, k, i + 1);
                    t(izz, 5) = y(j,     k, i)     - y(j + 1, k, i + 1);
                    t(izz, 6) = z(j,     k, i)     - z(j + 1, k, i + 1);
                    // cross product of vectors
                    sk(izz, i, 1) =  t(izz, 2) * t(izz, 6) - t(izz, 3) * t(izz, 5);
                    sk(izz, i, 2) = -t(izz, 1) * t(izz, 6) + t(izz, 3) * t(izz, 4);
                    sk(izz, i, 3) =  t(izz, 1) * t(izz, 5) - t(izz, 2) * t(izz, 4);
                    // magnitude of cross product
                    sk(izz, i, 4) = sk(izz, i, 1) * sk(izz, i, 1) + sk(izz, i, 2) * sk(izz, i, 2) +
                                    sk(izz, i, 3) * sk(izz, i, 3);
                    asum = asum + std::sqrt(sk(izz, i, 4));
                    if (sk(izz, i, 4) == 0.) sk(izz, i, 4) = atol2;
                }
            }

            asum = 0.5 * asum;
            if ((float)asum < (float)atol) {
                // collapsed metrics
                if (isklton > 0) {
                    if (m == 1) {
                        if (isklton == 1) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   singular metrics: k=1     i=%3d,%3d  j=%3d,%3d",
                                is, ie + 1, js, je + 1);
                        } else {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "      singular metrics: k=1     i=%3d,%3d  j=%3d,%3d",
                                is, ie + 1, js, je + 1);
                        }
                    } else {
                        if (isklton == 1) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   singular metrics: k=kdim  i=%3d,%3d  j=%3d,%3d",
                                is, ie + 1, js, je + 1);
                        } else {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "      singular metrics: k=kdim  i=%3d,%3d  j=%3d,%3d",
                                is, ie + 1, js, je + 1);
                        }
                    }
                }
                for (i = is; i <= ie; i++) {
                    for (j = js; j <= je; j++) {
                        izz    = (k - 1) * jdim + j;
                        izzndx = izz - (2 * jdim * m) + (3 * jdim);
                        // set directions on collapsed face equal to those at closest interior face
                        sk(izz, i, 1) = sk(izzndx, i, 1);
                        sk(izz, i, 2) = sk(izzndx, i, 2);
                        sk(izz, i, 3) = sk(izzndx, i, 3);
                        sk(izz, i, 4) = 0.e0;
                        sk(izz, i, 5) = sk(izzndx, i, 5);
                    }
                }
            } else {
                // non-singular metrics
                for (i = is; i <= ie; i++) {
                    for (j = js; j <= je; j++) {
                        izz = (k - 1) * jdim + j;
                        sk(izz, i, 4) = 1.e0 / std::sqrt(sk(izz, i, 4));
                        sk(izz, i, 1) = sk(izz, i, 1) * sk(izz, i, 4);
                        sk(izz, i, 2) = sk(izz, i, 2) * sk(izz, i, 4);
                        sk(izz, i, 3) = sk(izz, i, 3) * sk(izz, i, 4);
                        sk(izz, i, 4) = 0.5e0 / sk(izz, i, 4);
                        sk(izz, i, 5) = 0.0;
                    }
                }
            }
        } // end do 3010
    } // end do 3000


    // fill in extra values of sk for safety
    for (i = 1; i <= idim - 1; i++) {
        // set metrics at jdim
        // cdir$ ivdep
        for (k = 1; k <= kdim; k++) {
            izz  = jdim * (k - 1) + jdim;
            izz1 = izz - 1;
            sk(izz, i, 1) = sk(izz1, i, 1);
            sk(izz, i, 2) = sk(izz1, i, 2);
            sk(izz, i, 3) = sk(izz1, i, 3);
            sk(izz, i, 4) = sk(izz1, i, 4);
            sk(izz, i, 5) = sk(izz1, i, 5);
        }
    } // end do 3070

    // Don't do the following checks if iflag=-1:
    if (iflag != -1) {
        // Check for large changes in metric from one k to the next,
        // indicating problem with the grid
        for (i = 1; i <= idim - 1; i++) {
            for (j = 1; j <= jdim - 1; j++) {
                for (k = 1; k <= kdim - 1; k++) {
                    izz   = (k - 1) * jdim + j;
                    izzp1 = k * jdim + j;
                    delta1 = std::abs(sk(izzp1, i, 1) - sk(izz, i, 1));
                    delta2 = std::abs(sk(izzp1, i, 2) - sk(izz, i, 2));
                    delta3 = std::abs(sk(izzp1, i, 3) - sk(izz, i, 3));
                    if ((float)delta1 > 1.99f || (float)delta2 > 1.99f || (float)delta3 > 1.99f) {
                        iflag = 1;
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            "   FATAL sk grid normal direction change near j,k,k+1,i=");
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            "     ... suspect bad grid%5d%5d%5d%5d", j, k, k + 1, i);
                    } else if ((float)delta1 > 1.5f || (float)delta2 > 1.5f || (float)delta3 > 1.5f) {
                        icnt = icnt + 1;
                        if (icnt <= icntmax) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   WARNING: Dramatic sk grid norm direction change (>120deg)");
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "     near j,k,k+1,i=%5d%5d%5d%5d", j, k, k + 1, i);
                        } else if (icnt == icntmax + 1) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "   NOTE: no further WARNINGS will be output...better check this grid!");
                        }
                    }
                }
            }
        }
    } // end if iflag != -1

    return;
} // end metric

void readkey(int& ititr, int& myid, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou, int& iunit11, int& ierrflg)
{
    readkey_ns::readkey(ititr, myid, ibufdim, nbuf, bou, nou, iunit11, ierrflg);
}

} // namespace metric_ns
