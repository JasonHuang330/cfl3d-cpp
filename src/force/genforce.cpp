// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "genforce.h"
#include "force.h"
#include "runtime/common_blocks.h"
#include <cmath>

namespace genforce_ns {

void force(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> sk, FortranArray4DRef<double> sj, FortranArray4DRef<double> si, double& cl, double& cd, double& cz, double& cy, double& cx, double& cmy, double& cmx, double& cmz, double& chd, double& swet, int& i00, FortranArray3DRef<double> ub, FortranArray3DRef<double> vb, FortranArray3DRef<double> wb, FortranArray3DRef<double> vmuk, FortranArray3DRef<double> vmuj, FortranArray3DRef<double> vmui, FortranArray3DRef<double> vol, int& ifo, int& jfo, int& kfo, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray3DRef<double> blank, int& nbl, FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk, FortranArray4DRef<double> xtbi, int& iuns, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& nn, int& maxseg)
{
    force_ns::force(jdim, kdim, idim, x, y, z, sk, sj, si, cl, cd, cz, cy, cx,
                    cmy, cmx, cmz, chd, swet, i00, ub, vb, wb, vmuk, vmuj, vmui,
                    vol, ifo, jfo, kfo, bcj, bck, bci, blank, nbl, xtbj, xtbk,
                    xtbi, iuns, qj0, qk0, qi0, nbci0, nbcj0, nbck0, nbcidim,
                    nbcjdim, nbckdim, ibcinfo, jbcinfo, kbcinfo, nn, maxseg);
}

void genforce(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> sk, FortranArray4DRef<double> sj, FortranArray4DRef<double> si, double& cz, double& cy, double& cx, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray3DRef<double> blank, int& nbl, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, int& maxbl, int& maxseg, int& n2, FortranArray4DRef<double> xmdj, FortranArray4DRef<double> xmdk, FortranArray4DRef<double> xmdi, FortranArray2DRef<double> aesrfdat, int& nmds, int& maxaes, int& maxsegdg, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> idfrmseg, int& iaes, FortranArray2DRef<int> iaesurf)
{
    // COMMON block references
    float& gamma  = cmn_fluid.gamma;
    float& xmach  = cmn_info.xmach;
    float& p0     = cmn_ivals.p0;

    // cpc = 2.0 / (gamma * xmach * xmach)
    double cpc = 2.e0 / ((double)gamma * (double)xmach * (double)xmach);

    cz = 0.e0;
    cy = 0.e0;
    cx = 0.e0;

    // loop over all aeroelastic surface segments (those with idfrmseg=99)
    for (int iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {

        if (iaesurf(nbl, iseg) == iaes &&
            idfrmseg(nbl, iseg) == 99) {

            //******************************************************************
            // generalized forces on k=constant surfaces
            //******************************************************************

            if (kcsi(nbl, iseg) == kcsf(nbl, iseg)) {

                int k  = kcsi(nbl, iseg);
                int kk, kc, l;
                if (k == 1) {
                    kk = 1;
                    kc = 1;
                    l  = 0;
                } else if (k == kdim) {
                    kk = 2;
                    kc = kdim - 1;
                    l  = 3;
                } else {
                    kk = 0; kc = 0; l = 0; // should not happen
                }

                int ist = icsi(nbl, iseg);
                int ifn = icsf(nbl, iseg) - 1;
                int jst = jcsi(nbl, iseg);
                int jfn = jcsf(nbl, iseg) - 1;

                double sgn = 1.0;
                if (kk > 1) sgn = -1.0;

                for (int i = ist; i <= ifn; i++) {
                    double cxl = 0.e0;
                    double cyl = 0.e0;
                    double czl = 0.e0;

                    for (int j = jst; j <= jfn; j++) {

                        double xmdavg = (xmdk(j,   i,   l+1, n2) + xmdk(j+1, i,   l+1, n2)
                                       + xmdk(j,   i+1, l+1, n2) + xmdk(j+1, i+1, l+1, n2)) * .25;
                        double ymdavg = (xmdk(j,   i,   l+2, n2) + xmdk(j+1, i,   l+2, n2)
                                       + xmdk(j,   i+1, l+2, n2) + xmdk(j+1, i+1, l+2, n2)) * .25;
                        double zmdavg = (xmdk(j,   i,   l+3, n2) + xmdk(j+1, i,   l+3, n2)
                                       + xmdk(j,   i+1, l+3, n2) + xmdk(j+1, i+1, l+3, n2)) * .25;

                        double dcp = -(qk0(j, i, 5, kk+kk-1) / (double)p0 - 1.e0) * cpc * sk(j, k, i, 4);
                        double dcx = dcp * sk(j, k, i, 1) * sgn;
                        double dcy = dcp * sk(j, k, i, 2) * sgn;
                        double dcz = dcp * sk(j, k, i, 3) * sgn;

                        // only use contributions from points with interface
                        // (solid surface) boundary conditions, and only those
                        // points not blanked out
                        double fact = bck(j, i, kk) * blank(j, kc, i);
                        dcx = dcx * fact;
                        dcy = dcy * fact;
                        dcz = dcz * fact;

                        cxl = cxl + dcx * xmdavg;
                        cyl = cyl + dcy * ymdavg;
                        czl = czl + dcz * zmdavg;

                    } // j loop (label 15)

                    // integrated values
                    cz = cz + czl;
                    cy = cy + cyl;
                    cx = cx + cxl;

                } // i loop (label 10)

            } // k=constant surfaces

            //******************************************************************
            // generalized forces on j=constant surfaces
            //******************************************************************

            if (jcsi(nbl, iseg) == jcsf(nbl, iseg)) {

                int j  = jcsi(nbl, iseg);
                int jj, jc, l;
                if (j == 1) {
                    jj = 1;
                    jc = 1;
                    l  = 0;
                } else if (j == jdim) {
                    jj = 2;
                    jc = jdim - 1;
                    l  = 3;
                } else {
                    jj = 0; jc = 0; l = 0; // should not happen
                }

                int ist = icsi(nbl, iseg);
                int ifn = icsf(nbl, iseg) - 1;
                int kst = kcsi(nbl, iseg);
                int kfn = kcsf(nbl, iseg) - 1;

                double sgn = 1.0;
                if (jj > 1) sgn = -1.0;

                for (int i = ist; i <= ifn; i++) {
                    double cxl = 0.e0;
                    double cyl = 0.e0;
                    double czl = 0.e0;

                    for (int k = kst; k <= kfn; k++) {

                        double xmdavg = (xmdj(k,   i,   l+1, n2) + xmdj(k+1, i,   l+1, n2)
                                       + xmdj(k,   i+1, l+1, n2) + xmdj(k+1, i+1, l+1, n2)) * .25;
                        double ymdavg = (xmdj(k,   i,   l+2, n2) + xmdj(k+1, i,   l+2, n2)
                                       + xmdj(k,   i+1, l+2, n2) + xmdj(k+1, i+1, l+2, n2)) * .25;
                        double zmdavg = (xmdj(k,   i,   l+3, n2) + xmdj(k+1, i,   l+3, n2)
                                       + xmdj(k,   i+1, l+3, n2) + xmdj(k+1, i+1, l+3, n2)) * .25;

                        double dcp = -(qj0(k, i, 5, jj+jj-1) / (double)p0 - 1.e0) * cpc * sj(j, k, i, 4);
                        double dcx = dcp * sj(j, k, i, 1) * sgn;
                        double dcy = dcp * sj(j, k, i, 2) * sgn;
                        double dcz = dcp * sj(j, k, i, 3) * sgn;

                        // only use contributions from points with interface
                        // (solid surface) boundary conditions, and only those
                        // points not blanked out
                        double fact = bcj(k, i, jj) * blank(jc, k, i);
                        dcx = dcx * fact;
                        dcy = dcy * fact;
                        dcz = dcz * fact;

                        cxl = cxl + dcx * xmdavg;
                        cyl = cyl + dcy * ymdavg;
                        czl = czl + dcz * zmdavg;

                    } // k loop (label 45)

                    // integrated values
                    cz = cz + czl;
                    cy = cy + cyl;
                    cx = cx + cxl;

                } // i loop (label 40)

            } // j=constant surfaces

            //******************************************************************
            // generalized forces on i=constant surfaces
            //******************************************************************

            if (icsi(nbl, iseg) == icsf(nbl, iseg)) {

                int i  = icsi(nbl, iseg);
                int ii, ic, l;
                if (i == 1) {
                    ii = 1;
                    ic = 1;
                    l  = 0;
                } else if (i == idim) {
                    ii = 2;
                    ic = idim - 1;
                    l  = 3;
                } else {
                    ii = 0; ic = 0; l = 0; // should not happen
                }

                int jst = jcsi(nbl, iseg);
                int jfn = jcsf(nbl, iseg) - 1;
                int kst = kcsi(nbl, iseg);
                int kfn = kcsf(nbl, iseg) - 1;

                double sgn = 1.0;
                if (ii > 1) sgn = -1.0;

                for (int j = jst; j <= jfn; j++) {
                    double cxl = 0.e0;
                    double cyl = 0.e0;
                    double czl = 0.e0;

                    for (int k = kst; k <= kfn; k++) {

                        double xmdavg = (xmdi(j,   k,   l+1, n2) + xmdi(j+1, k,   l+1, n2)
                                       + xmdi(j,   k+1, l+1, n2) + xmdi(j+1, k+1, l+1, n2)) * .25;
                        double ymdavg = (xmdi(j,   k,   l+2, n2) + xmdi(j+1, k,   l+2, n2)
                                       + xmdi(j,   k+1, l+2, n2) + xmdi(j+1, k+1, l+2, n2)) * .25;
                        double zmdavg = (xmdi(j,   k,   l+3, n2) + xmdi(j+1, k,   l+3, n2)
                                       + xmdi(j,   k+1, l+3, n2) + xmdi(j+1, k+1, l+3, n2)) * .25;

                        double dcp = -(qi0(j, k, 5, ii+ii-1) / (double)p0 - 1.e0) * cpc * si(j, k, i, 4);
                        double dcx = dcp * si(j, k, i, 1) * sgn;
                        double dcy = dcp * si(j, k, i, 2) * sgn;
                        double dcz = dcp * si(j, k, i, 3) * sgn;

                        // only use contributions from points with interface
                        // (solid surface) boundary conditions, and only those
                        // points not blanked out
                        double fact = bci(j, k, ii) * blank(j, k, ic);
                        dcx = dcx * fact;
                        dcy = dcy * fact;
                        dcz = dcz * fact;

                        cxl = cxl + dcx * xmdavg;
                        cyl = cyl + dcy * ymdavg;
                        czl = czl + dcz * zmdavg;

                    } // k loop (label 75)

                    // integrated values
                    cz = cz + czl;
                    cy = cy + cyl;
                    cx = cx + cxl;

                } // j loop (label 70)

            } // i=constant surfaces

        } // iaesurf/idfrmseg check

    } // iseg loop

    return;
}

} // namespace genforce_ns
