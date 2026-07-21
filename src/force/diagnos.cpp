// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "diagnos.h"
#include "runtime/fortran_io.h"
#include "extra.h"
#include "extrae.h"
#include "direct.h"
#include "termn8.h"
#include "ccomplex.h"
#include <cstdio>
#include <cstring>
#include <cmath>

namespace diagnos_ns {

void diagnos(int& istop, int& iout, int& igap, int& jdim1, int& kdim1, int& msub1, int& msub2,
             FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, int& lmax1,
             FortranArray3DRef<double> x1, FortranArray3DRef<double> y1, FortranArray3DRef<double> z1,
             FortranArray3DRef<double> x1mid, FortranArray3DRef<double> y1mid, FortranArray3DRef<double> z1mid,
             FortranArray3DRef<double> x1mide, FortranArray3DRef<double> y1mide, FortranArray3DRef<double> z1mide,
             FortranArray3DRef<double> x2int, FortranArray3DRef<double> y2int, FortranArray3DRef<double> z2int,
             FortranArray3DRef<double> x2fit, FortranArray2DRef<double> y2fit, FortranArray2DRef<double> z2fit,
             FortranArray1DRef<int> jjmax2, FortranArray1DRef<int> kkmax2,
             FortranArray3DRef<double> x2, FortranArray3DRef<double> y2, FortranArray3DRef<double> z2,
             FortranArray1DRef<double> xie2, FortranArray1DRef<double> eta2,
             FortranArray1DRef<int> mblkpt, int& icheck, int& intmx,
             double& xc, double& yc, double& zc, int& ifit,
             int& j21, int& j22, int& k21, int& k22, int& npt, int& ic0, int& iorph,
             FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2,
             FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2,
             int& itoss0, int& iself,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& myid)
{
    // COMMON block references
    int& ip3dgrd = cmn_igrdtyp.ip3dgrd;
    int& ialph   = cmn_igrdtyp.ialph;
    float& epsc  = cmn_tol.epsc;

    // Local variables
    int jmax2, kmax2;
    int neta, nxie;
    int norph;
    int ierr;
    int nc00, nc0;
    int nct, iflg;
    int nbpt, nnorm;
    double errmx, c0tol;
    double errmax, percnt;
    double epsn, prdsum, dnmx;
    int jerrmx, kerrmx, jdnmx, kdnmx;
    int j, k, ll, ll2, ll3, ll4, l, l1, jc, kc;
    double xiec, etac, xie, eta, xiejk, etajk;
    double axie, aeta, bxie, beta, bj0, bj;
    double x11, x12, x13, x14, x15, x16, x17, x18;
    double y11, y12, y13, y14, y15, y16, y17, y18;
    double z11, z12, z13, z14, z15, z16, z17, z18;
    double dx2, dx3, dx4, dx5, dx6, dx7, dx8;
    double dy2, dy3, dy4, dy5, dy6, dy7, dy8;
    double dz2, dz3, dz4, dz5, dz6, dz7, dz8;
    double a1, a2, a3, a4, a5, a6, a7, a8;
    double b2, b3, b4, b5, b6, b7, b8;
    double c2, c3, c4, c5, c6, c7, c8;
    double df1, df2, df3;
    double x2c, y2c, z2c;
    double x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8;
    double x00, y00, z00;
    double xscal, yscal, zscal, xsc, ysc, zsc;
    double xerr, yerr, zerr, error;
    int itoss, icase;
    double ax, ay, az, bx, by, bz, dnx, dny, dnz, d;
    double xnfit, ynfit, znfit, xnint, ynint, znint, prod, dn;
    int iunit;
    char titlptchgrd[15];
    int len1;
    int js, ks, je, ke;
    int iwhole, maxx;
    int jseg, kseg;
    int jtest, ktest, jdif24, kdif24, jdif31, kdif31, ibrnch;
    int j11, k11;
    int kcall, jcall;
    int neg1 = -1;



    jmax2 = jjmax2(1) - 1;
    kmax2 = kkmax2(1) - 1;

    // diagnostics for case where search/inversion routines have failed
    if (istop == 1) {
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120, "            possible causes:");
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120, "              ...\"from\" block input data incorrectly set");
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120, "              ...\"to\" and \"from\" blocks do not lie on the same surface");
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120, "              ...factj or factk too small");
        goto label_552;
    }

    // diagnostics for completed generalized coordinate interpolation
    nou(4) = std::min(nou(4)+1, ibufdim);
    std::snprintf(bou(nou(4),4), 120, "    beginning diagnostic checks");

    // ****** first check ******
    // check final cell center generalized coordinates
    neta = 0;
    nxie = 0;

    for (k = k21; k <= k22-1; k++) {
    for (j = j21; j <= j22-1; j++) {
        ll = (j22-j21)*(k-k21) + (j-j21+1);
        l1 = mblkpt(ll);
        if (l1 > 0) {
            if ((float)eta2(ll) < 1.f || (float)eta2(ll) > (float)(kkmax1(l1)-2)) {
                nou(2) = std::min(nou(2)+1, ibufdim);
                std::snprintf(bou(nou(2),2), 120,
                    "interp. no.: %4d, \"to\" cell center j,k = %4d%4d found in \"from\" block %4d",
                    icheck, j, k, l1);
                nou(2) = std::min(nou(2)+1, ibufdim);
                std::snprintf(bou(nou(2),2), 120,
                    "   with eta = %14.7e", (float)eta2(ll));
                nou(2) = std::min(nou(2)+1, ibufdim);
                std::snprintf(bou(nou(2),2), 120,
                    "   in this \"from\" block, legal range of eta is 1 to %4d", kkmax1(l1)-2);
                neta = neta + 1;
            }
        }
    }}

    for (j = j21; j <= j22-1; j++) {
    for (k = k21; k <= k22-1; k++) {
        ll = (j22-j21)*(k-k21) + (j-j21+1);
        l1 = mblkpt(ll);
        if (l1 > 0) {
            if ((float)xie2(ll) < 1.f || (float)xie2(ll) > (float)(jjmax1(l1)-2)) {
                nou(2) = std::min(nou(2)+1, ibufdim);
                std::snprintf(bou(nou(2),2), 120,
                    "interp. no.: %4d, \"to\" cell center j,k = %4d%4d found in \"from\" block %4d",
                    icheck, j, k, l1);
                nou(2) = std::min(nou(2)+1, ibufdim);
                std::snprintf(bou(nou(2),2), 120,
                    "   with xie = %14.7e", (float)xie2(ll));
                nou(2) = std::min(nou(2)+1, ibufdim);
                std::snprintf(bou(nou(2),2), 120,
                    "   in this \"from\" block, legal range of xie is 1 to %4d", jjmax1(l1)-2);
                nxie = nxie + 1;
            }
        }
    }}

    nou(4) = std::min(nou(4)+1, ibufdim);
    std::snprintf(bou(nou(4),4), 120, "     %4d points were found outside the legal range of xie", nxie);
    nou(4) = std::min(nou(4)+1, ibufdim);
    std::snprintf(bou(nou(4),4), 120, "     %4d points were found outside the legal range of eta", neta);

    if (neta+nxie > 0) {
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120, "        check fort.9 to make sure values are not far outside legal range");
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120, "            possible causes:");
        if (igap > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "              ...small gaps between \"from\" blocks");
        } else {
            if (nxie > 0) {
                nou(4) = std::min(nou(4)+1, ibufdim);
                std::snprintf(bou(nou(4),4), 120, "              ...mmcxie incorrectly set");
            }
            if (neta > 0) {
                nou(4) = std::min(nou(4)+1, ibufdim);
                std::snprintf(bou(nou(4),4), 120, "              ...mmceta incorrectly set");
            }
        }
        if (iorph > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "              ...these may be the orphan points anticipated by setting iorph > 0");
        }
        iout = 1;
    }



    // ****** second check (grids with orphan points only) ******
    // count number of orphans
    if (iorph > 0) {
        norph = 0;
        for (j = j21; j <= j22-1; j++) {
        for (k = k21; k <= k22-1; k++) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            l1 = mblkpt(ll);
            if (l1 == 0) norph = norph + 1;
        }}
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120, "     %4d points were flagged as orphans", norph);
    }

    // ****** third check (grids that self-connect with a branch cut) ******
    // check for points that get interpolated from themselves (branch cut error)
    if (iself > 0) {
        ierr = 0;
        for (j = j21; j <= j22-1; j++) {
        for (k = k21; k <= k22-1; k++) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            l1 = mblkpt(ll);
            jc = (int)xie2(ll);
            kc = (int)eta2(ll);
            if (jc == j && kc == k) {
                ierr = ierr + 1;
                nou(2) = std::min(nou(2)+1, ibufdim);
                std::snprintf(bou(nou(2),2), 120, "branch cut problem at j,k =%4d%4d", j, k);
            }
        }}
        if (ierr > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "     %4d points are interpolated from themselves", ierr);
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120,
                "          search routine failure due to branch cut - DO NOT run flow solver");
        }
    }

    // ****** fourth check (C-0 continuous grids only) ******
    if (ic0 > 0) {
        c0tol  = (double)epsc;
        nc00   = (j22-j21)*(k22-k21);
        nc0    = nc00;
        errmx  = 0.;
        jerrmx = 1;
        kerrmx = 1;

        for (j = j21; j <= j22-1; j++) {
        for (k = k21; k <= k22-1; k++) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            l1 = mblkpt(ll);
            jc = (int)xie2(ll);
            kc = (int)eta2(ll);
            {
                double xiec_tmp = xie2(ll) - jc - 0.5;
                double etac_tmp = eta2(ll) - kc - 0.5;
                xiec = ccomplex_ns::ccabs(xiec_tmp);
                etac = ccomplex_ns::ccabs(etac_tmp);
            }
            if ((float)xiec > (float)errmx) {
                jerrmx = j;
                kerrmx = k;
                errmx  = xiec;
            }
            if ((float)etac > (float)errmx) {
                jerrmx = j;
                kerrmx = k;
                errmx  = etac;
            }
            if ((float)xiec > (float)c0tol || (float)etac > (float)c0tol) {
                nc0 = nc0 - 1;
            }
        }}
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "     %4d points  flagged as C-0 have |xie-.5| or |eta-.5| >%10.3e",
            nc00-nc0, (float)c0tol);
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "          the max. error %10.3e occurs at j,k = %4d,%4d",
            (float)errmx, jerrmx, kerrmx);
        if (nc00-nc0 != 0) {
            iout = 1;
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "            possible causes:");
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "              ...interface not really C-0, check grid and/or C-0 flag");
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "              ...\"from\" block input data incorrectly set");
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "              ...search routine error...uh-oh");
        } else {
            for (j = j21; j <= j22-1; j++) {
            for (k = k21; k <= k22-1; k++) {
                ll = (j22-j21)*(k-k21) + (j-j21+1);
                jc = (int)xie2(ll);
                kc = (int)eta2(ll);
                xie2(ll) = jc + 0.5;
                eta2(ll) = kc + 0.5;
            }}
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120,
                "          C-0 check acceptable, resetting xie and eta to exact values (0.5,0.5)");
        }
    }



    // ****** fifth check ******
    // check jacobians of grid in generalized coordinates for anomolies (negative areas)
    // not applicable in 2D cases or cases with orphan points!
    jseg = j22 - j21;
    kseg = k22 - k21;
    if (jseg > 1 && kseg > 1 && iorph <= 0) {

        nct  = 0;
        iflg = 0;
        for (l = 1; l <= lmax1; l++) {

            // find j11,k11, the first "to" cell located in "from" block l
            j11 = -1;
            k11 = -1;
            for (j = j21; j <= j22-2; j++) {
            for (k = k21; k <= k22-2; k++) {
                ll  = (j22-j21)*(k-k21)   + (j-j21+1);
                ll2 = (j22-j21)*(k-k21)   + (j+1-j21+1);
                ll3 = (j22-j21)*(k+1-k21) + (j+1-j21+1);
                ll4 = (j22-j21)*(k+1-k21) + (j-j21+1);
                if (mblkpt(ll) == l && mblkpt(ll2) == l &&
                    mblkpt(ll3) == l && mblkpt(ll4) == l) {
                    jtest = jjmax1(mblkpt(ll))/2;
                    ktest = kkmax1(mblkpt(ll))/2;
                    jdif24 = std::abs((int)(xie2(ll2) - xie2(ll4)));
                    kdif24 = std::abs((int)(eta2(ll2) - eta2(ll4)));
                    jdif31 = std::abs((int)(xie2(ll3) - xie2(ll)));
                    kdif31 = std::abs((int)(eta2(ll3) - eta2(ll)));
                    ibrnch = 0;
                    if (jdif24 > jtest || kdif24 > ktest ||
                        jdif31 > jtest || kdif31 > ktest) ibrnch = 1;
                    if (ibrnch == 0) {
                        j11 = j;
                        k11 = k;
                        goto label_1756;
                    }
                }
            }}
            label_1756:;
            if (j11 < 1 || k11 < 1) continue; // go to 1750
            iflg = 1;
            j  = j11;
            k  = k11;
            ll  = (j22-j21)*(k-k21)   + (j-j21+1);
            ll2 = (j22-j21)*(k-k21)   + (j+1-j21+1);
            ll3 = (j22-j21)*(k+1-k21) + (j+1-j21+1);
            ll4 = (j22-j21)*(k+1-k21) + (j-j21+1);
            axie = xie2(ll3) - xie2(ll);
            aeta = eta2(ll3) - eta2(ll);
            bxie = xie2(ll4) - xie2(ll2);
            beta = eta2(ll4) - eta2(ll2);
            bj0  = axie*beta - bxie*aeta;

            for (j = j21; j <= j22-2; j++) {
            for (k = k21; k <= k22-2; k++) {
                ll  = (j22-j21)*(k-k21)   + (j-j21+1);
                ll2 = (j22-j21)*(k-k21)   + (j+1-j21+1);
                ll3 = (j22-j21)*(k+1-k21) + (j+1-j21+1);
                ll4 = (j22-j21)*(k+1-k21) + (j-j21+1);
                if (mblkpt(ll) == l && mblkpt(ll2) == l &&
                    mblkpt(ll3) == l && mblkpt(ll4) == l) {
                    jtest = jjmax1(mblkpt(ll))/2;
                    ktest = kkmax1(mblkpt(ll))/2;
                    jdif24 = std::abs((int)(xie2(ll2) - xie2(ll4)));
                    kdif24 = std::abs((int)(eta2(ll2) - eta2(ll4)));
                    jdif31 = std::abs((int)(xie2(ll3) - xie2(ll)));
                    kdif31 = std::abs((int)(eta2(ll3) - eta2(ll)));
                    ibrnch = 0;
                    if (jdif24 > jtest || kdif24 > ktest ||
                        jdif31 > jtest || kdif31 > ktest) ibrnch = 1;
                    if (ibrnch == 0) {
                        axie = xie2(ll3) - xie2(ll);
                        aeta = eta2(ll3) - eta2(ll);
                        bxie = xie2(ll4) - xie2(ll2);
                        beta = eta2(ll4) - eta2(ll2);
                        bj   = axie*beta - bxie*aeta;
                        if ((float)bj * (float)bj0 <= 0.f) {
                            nct = nct + 1;
                            nou(2) = std::min(nou(2)+1, ibufdim);
                            std::snprintf(bou(nou(2),2), 120,
                                "non-unique point at j,k = %5d,%5d", j, k);
                            nou(2) = std::min(nou(2)+1, ibufdim);
                            std::snprintf(bou(nou(2),2), 120,
                                "xie2(j,k),eta2(j,k)         = %11.5f,%11.5f",
                                (float)xie2(ll), (float)eta2(ll));
                            nou(2) = std::min(nou(2)+1, ibufdim);
                            std::snprintf(bou(nou(2),2), 120,
                                "xie2(j+1,k),eta2(j+1,k)     = %11.5f,%11.5f",
                                (float)xie2(ll2), (float)eta2(ll2));
                            nou(2) = std::min(nou(2)+1, ibufdim);
                            std::snprintf(bou(nou(2),2), 120,
                                "xie2(j+1,k+1),eta2(j+1,k+1) = %11.5f,%11.5f",
                                (float)xie2(ll3), (float)eta2(ll3));
                            nou(2) = std::min(nou(2)+1, ibufdim);
                            std::snprintf(bou(nou(2),2), 120,
                                "xie2(j,k+1),eta2(j,k+1)     = %11.5f,%11.5f",
                                (float)xie2(ll4), (float)eta2(ll4));
                        }
                    }
                }
            }}
        } // end for l

        if (iflg == 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120,
                "  caution...no consecutive \"to\" cells found in any of the \"from\" blocks");
        }
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "     %4d non-unique points were found in the generalized-coordinate ", nct);
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "       mapping between the \"to\" and \"from\" grids");
        if (nct > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "        Caution: this may indicate a serious problem");
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "            possible causes:");
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "              ...mmcxie and/or mmceta incorrectly set");
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120,
                "              ...\"to\" and \"from\" blocks do not lie on the same surface");
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120,
                "              ...branch cut on the \"from\" side may falsely trigger this message");
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120,
                "              ...and the results may in fact be OK - verify with patch_p3d.xx file");
            iout = 1;
        }
    } // end if jseg>1 && kseg>1 && iorph<=0 (fifth check)



    // ****** sixth check ******
    // compute cell centers of "to" grid from interpolation coefficients
    // and compare with direct calculation from "to" grid
    nbpt   = 0;
    percnt = 0.10;
    errmax = 0.;
    jerrmx = 1;
    kerrmx = 1;

    for (j = j21; j <= j22-1; j++) {
    for (k = k21; k <= k22-1; k++) {
        ll = (j22-j21)*(k-k21) + (j-j21+1);
        l1 = mblkpt(ll);

        // skip if orphan point
        if (l1 == 0) continue;

        jc = (int)xie2(ll);
        kc = (int)eta2(ll);

        // add 1 to account for expanded "from" grids
        jc = jc + 1;
        kc = kc + 1;
        xiejk = xie2(ll) + 1;
        etajk = eta2(ll) + 1;
        xiec = (double)jc;
        etac = (double)kc;

        x11 = x1(jc,   kc,   l1);  y11 = y1(jc,   kc,   l1);  z11 = z1(jc,   kc,   l1);
        x12 = x1(jc+1, kc,   l1);  y12 = y1(jc+1, kc,   l1);  z12 = z1(jc+1, kc,   l1);
        x13 = x1(jc+1, kc+1, l1);  y13 = y1(jc+1, kc+1, l1);  z13 = z1(jc+1, kc+1, l1);
        x14 = x1(jc,   kc+1, l1);  y14 = y1(jc,   kc+1, l1);  z14 = z1(jc,   kc+1, l1);
        x15 = x1mid(jc,   kc,   l1); y15 = y1mid(jc,   kc,   l1); z15 = z1mid(jc,   kc,   l1);
        x16 = x1mid(jc,   kc+1, l1); y16 = y1mid(jc,   kc+1, l1); z16 = z1mid(jc,   kc+1, l1);
        x17 = x1mide(jc,   kc,   l1); y17 = y1mide(jc,   kc,   l1); z17 = z1mide(jc,   kc,   l1);
        x18 = x1mide(jc+1, kc,   l1); y18 = y1mide(jc+1, kc,   l1); z18 = z1mide(jc+1, kc,   l1);

        dx2 = x12-x11; dy2 = y12-y11; dz2 = z12-z11;
        dx3 = x13-x11; dy3 = y13-y11; dz3 = z13-z11;
        dx4 = x14-x11; dy4 = y14-y11; dz4 = z14-z11;
        dx5 = x15-x11; dy5 = y15-y11; dz5 = z15-z11;
        dx6 = x16-x11; dy6 = y16-y11; dz6 = z16-z11;
        dx7 = x17-x11; dy7 = y17-y11; dz7 = z17-z11;
        dx8 = x18-x11; dy8 = y18-y11; dz8 = z18-z11;

        if (ifit == 1) {
            // bi-linear fit
            a2=dx2; a3=dx4; a4=dx3-a2-a3; a5=0.; a6=0.; a7=0.; a8=0.;
            b2=dy2; b3=dy4; b4=dy3-b2-b3; b5=0.; b6=0.; b7=0.; b8=0.;
            c2=dz2; c3=dz4; c4=dz3-c2-c3; c5=0.; c6=0.; c7=0.; c8=0.;
        }

        if (ifit == 2) {
            // (degenerate) bi-quadratic fit
            a2=-dx2+4.*dx5; a3=-dx4+4.*dx7; a5=2.*dx2-4.*dx5; a7=2.*dx4-4.*dx7;
            df1=dx3-a2-a3-a5-a7; df2=dx6-.5*a2-a3-.25*a5-a7; df3=dx8-a2-.5*a3-a5-.25*a7;
            a4=-3.*df1+4.*df2+4.*df3; a6=2.*df1-4.*df2; a8=2.*df1-4.*df3;
            b2=-dy2+4.*dy5; b3=-dy4+4.*dy7; b5=2.*dy2-4.*dy5; b7=2.*dy4-4.*dy7;
            df1=dy3-b2-b3-b5-b7; df2=dy6-.5*b2-b3-.25*b5-b7; df3=dy8-b2-.5*b3-b5-.25*b7;
            b4=-3.*df1+4.*df2+4.*df3; b6=2.*df1-4.*df2; b8=2.*df1-4.*df3;
            c2=-dz2+4.*dz5; c3=-dz4+4.*dz7; c5=2.*dz2-4.*dz5; c7=2.*dz4-4.*dz7;
            df1=dz3-c2-c3-c5-c7; df2=dz6-.5*c2-c3-.25*c5-c7; df3=dz8-c2-.5*c3-c5-.25*c7;
            c4=-3.*df1+4.*df2+4.*df3; c6=2.*df1-4.*df2; c8=2.*df1-4.*df3;
        }

        // quadratic fit in xie, linear fit in eta
        if (ifit == 3) {
            a3=dx4; b3=dy4; c3=dz4;
            a2=-dx2+4.*dx5; b2=-dy2+4.*dy5; c2=-dz2+4.*dz5;
            a5=2.*dx2-4.*dx5; b5=2.*dy2-4.*dy5; c5=2.*dz2-4.*dz5;
            df1=x13-x12-a3; df2=x16-x15-a3; a4=-df1+4.*df2; a6=2.*df1-4.*df2;
            df1=y13-y12-b3; df2=y16-y15-b3; b4=-df1+4.*df2; b6=2.*df1-4.*df2;
            df1=z13-z12-c3; df2=z16-z15-c3; c4=-df1+4.*df2; c6=2.*df1-4.*df2;
            a7=0.; a8=0.; b7=0.; b8=0.; c7=0.; c8=0.;
        }

        // linear fit in xie, quadratic fit in eta
        if (ifit == 4) {
            a2=dx2; a3=-dx4+4.*dx7; a7=2.*dx4-4.*dx7;
            df1=dx3-a2-a3-a7; df2=dx8-a2-.5*a3-.25*a7;
            a4=-df1+4.*df2; a8=2.*df1-4.*df2; a5=0.; a6=0.;
            b2=dy2; b3=-dy4+4.*dy7; b7=2.*dy4-4.*dy7;
            df1=dy3-b2-b3-b7; df2=dy8-b2-.5*b3-.25*b7;
            b4=-df1+4.*df2; b8=2.*df1-4.*df2; b5=0.; b6=0.;
            c2=dz2; c3=-dz4+4.*dz7; c7=2.*dz4-4.*dz7;
            df1=dz3-c2-c3-c7; df2=dz8-c2-.5*c3-.25*c7;
            c4=-df1+4.*df2; c8=2.*df1-4.*df2; c5=0.; c6=0.;
        }

        xie = xiejk - xiec;
        eta = etajk - etac;

        x2int(j,k,1) = x11 + a3*eta + eta*(a7*eta + a8*xie*eta)
                     + xie*(a2 + a4*eta + a5*xie + a6*xie*eta);
        y2int(j,k,1) = y11 + b3*eta + eta*(b7*eta + b8*xie*eta)
                     + xie*(b2 + b4*eta + b5*xie + b6*xie*eta);
        z2int(j,k,1) = z11 + c3*eta + eta*(c7*eta + c8*xie*eta)
                     + xie*(c2 + c4*eta + c5*xie + c6*xie*eta);



        // compute center of "to" cell directly from the "to" grid, consistent with ifit
        {
            int jl = 1;
            int jr = jmax2;
            int kl = 1;
            int kr = kmax2;
            kcall = k+1;
            int l_one = 1;
            extra_ns::extra(jdim1, kdim1, msub2, l_one, x2, y2, z2,
                            j, kcall, jl, jr, x6, y6, z6, icase, ifit);
            extra_ns::extra(jdim1, kdim1, msub2, l_one, x2, y2, z2,
                            j, k, jl, jr, x5, y5, z5, icase, ifit);
            extrae_ns::extrae(jdim1, kdim1, msub2, l_one, x2, y2, z2,
                              j, k, kl, kr, x7, y7, z7, icase, ifit);
            jcall = j+1;
            extrae_ns::extrae(jdim1, kdim1, msub2, l_one, x2, y2, z2,
                              jcall, k, kl, kr, x8, y8, z8, icase, ifit);
        }

        // bi-linear
        if (ifit == 1) {
            x2c = 0.25*(x2(j,k,1) + x2(j+1,k,1) + x2(j+1,k+1,1) + x2(j,k+1,1));
            y2c = 0.25*(y2(j,k,1) + y2(j+1,k,1) + y2(j+1,k+1,1) + y2(j,k+1,1));
            z2c = 0.25*(z2(j,k,1) + z2(j+1,k,1) + z2(j+1,k+1,1) + z2(j,k+1,1));
        }
        // bi-quadratic
        if (ifit == 2) {
            x2c = 0.5*(x5+x6+x7+x8) - 0.25*(x2(j,k,1)+x2(j+1,k,1)+x2(j+1,k+1,1)+x2(j,k+1,1));
            y2c = 0.5*(y5+y6+y7+y8) - 0.25*(y2(j,k,1)+y2(j+1,k,1)+y2(j+1,k+1,1)+y2(j,k+1,1));
            z2c = 0.5*(z5+z6+z7+z8) - 0.25*(z2(j,k,1)+z2(j+1,k,1)+z2(j+1,k+1,1)+z2(j,k+1,1));
        }
        // quadratic in xie, linear in eta
        if (ifit == 3) {
            x2c = .5*(x5+x6);
            y2c = .5*(y5+y6);
            z2c = .5*(z5+z6);
        }
        // linear in xie, quadratic in eta
        if (ifit == 4) {
            x2c = .5*(x7+x8);
            y2c = .5*(y7+y8);
            z2c = .5*(z7+z8);
        }

        x2fit(j,k,1) = x2c;
        y2fit(j,k)   = y2c;
        z2fit(j,k)   = z2c;

        // compare interpolated vs direct cell center locations
        x00 = x2(j,k,1);
        y00 = y2(j,k,1);
        z00 = z2(j,k,1);
        xscal = 0.;
        yscal = 0.;
        zscal = 0.;
        for (int jj = 1; jj <= 2; jj++) {
        for (int kk = 1; kk <= 2; kk++) {
            double xsc_tmp = x2(j+jj-1,k+kk-1,1) - x00;
            xsc = ccomplex_ns::ccabs(xsc_tmp);
            if ((float)xsc > (float)xscal) xscal = xsc;
            double ysc_tmp = y2(j+jj-1,k+kk-1,1) - y00;
            ysc = ccomplex_ns::ccabs(ysc_tmp);
            if ((float)ysc > (float)yscal) yscal = ysc;
            double zsc_tmp = z2(j+jj-1,k+kk-1,1) - z00;
            zsc = ccomplex_ns::ccabs(zsc_tmp);
            if ((float)zsc > (float)zscal) zscal = zsc;
        }}
        if ((float)xscal == 0.f) xscal = -1.;
        if ((float)yscal == 0.f) yscal = -1.;
        if ((float)zscal == 0.f) zscal = -1.;

        // check errors only in transverse directions - not projected direction
        xerr = 0.;
        yerr = 0.;
        zerr = 0.;
        if (itoss0 == 0) {
            direct_ns::direct(x15, x16, x17, x18, y15, y16, y17, y18,
                              z15, z16, z17, z18, a1, a2, a3, itoss,
                              nou, bou, nbuf, ibufdim);
        } else {
            itoss = itoss0;
        }
        if (itoss == 1) {
            double tmp1 = y2int(j,k,1) - y2c;
            double tmp2 = z2int(j,k,1) - z2c;
            yerr = ccomplex_ns::ccabs(tmp1) / yscal;
            zerr = ccomplex_ns::ccabs(tmp2) / zscal;
        }
        if (itoss == 2) {
            double tmp1 = x2int(j,k,1) - x2c;
            double tmp2 = z2int(j,k,1) - z2c;
            xerr = ccomplex_ns::ccabs(tmp1) / xscal;
            zerr = ccomplex_ns::ccabs(tmp2) / zscal;
        }
        if (itoss == 3) {
            double tmp1 = x2int(j,k,1) - x2c;
            double tmp2 = y2int(j,k,1) - y2c;
            xerr = ccomplex_ns::ccabs(tmp1) / xscal;
            yerr = ccomplex_ns::ccabs(tmp2) / yscal;
        }
        error = xerr;
        if ((float)yerr > (float)error) error = yerr;
        if ((float)zerr > (float)error) error = zerr;
        if ((float)error > (float)errmax) {
            errmax = error;
            jerrmx = j;
            kerrmx = k;
        }
        if ((float)xerr > (float)percnt ||
            (float)yerr > (float)percnt ||
            (float)zerr > (float)percnt) {
            nbpt = nbpt + 1;
            nou(2) = std::min(nou(2)+1, ibufdim);
            std::snprintf(bou(nou(2),2), 120,
                " interpolation no. %5d for to cell j,k %5d,%5d", icheck, j, k);
            nou(2) = std::min(nou(2)+1, ibufdim);
            std::snprintf(bou(nou(2),2), 120,
                "    xint,yint,zint= %11.5f,%11.5f,%11.5f",
                (float)x2int(j,k,1), (float)y2int(j,k,1), (float)z2int(j,k,1));
            nou(2) = std::min(nou(2)+1, ibufdim);
            std::snprintf(bou(nou(2),2), 120,
                "    x,y,z= %11.5f,%11.5f,%11.5f",
                (float)x2c, (float)y2c, (float)z2c);
            nou(2) = std::min(nou(2)+1, ibufdim);
            std::snprintf(bou(nou(2),2), 120,
                "    xscal,yscal,zscal= %11.5f,%11.5f,%11.5f",
                (float)xscal, (float)yscal, (float)zscal);
        }
    }} // end for j,k (sixth check)



    nou(4) = std::min(nou(4)+1, ibufdim);
    std::snprintf(bou(nou(4),4), 120,
        "     %4d interpolated cell centers differ by more than %3d percent",
        nbpt, (int)((float)percnt*100.f));
    nou(4) = std::min(nou(4)+1, ibufdim);
    std::snprintf(bou(nou(4),4), 120,
        "          from those obtained directly from the grid points");
    nou(4) = std::min(nou(4)+1, ibufdim);
    std::snprintf(bou(nou(4),4), 120,
        "          the max. difference %10.3e percent occurs at j,k = %4d,%4d",
        (float)(errmax*100.), jerrmx, kerrmx);
    if (nbpt > 0) {
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "           note: differences between interpolated and directly obtained cell");
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "           centers often (correctly) arise when boundaries are being rendered");
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "           coincident. however, if differences occur when boundaries are");
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "           not rendered coincident, or if the other diagnostic counts are");
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "           non-zero, a problem in the input or the grid is likely");
    }

    // ****** seventh check ******
    // check normals on (cell center) grids
    // not applicable in 2D cases or for cases with orphans!
    jseg = j22 - j21;
    kseg = k22 - k21;
    if (jseg > 1 && kseg > 1 && iorph <= 0) {

        epsn   = .1;
        nnorm  = 0;
        prdsum = 0.;
        dnmx   = 0.;
        jdnmx  = 1;
        kdnmx  = 1;

        for (j = j21; j <= j22-2; j++) {
        for (k = k21; k <= k22-2; k++) {

            // normals at "to" cell centers calculated directly from "to" grid points
            ax = x2fit(j+1,k+1,1) - x2fit(j,k,1);
            ay = y2fit(j+1,k+1)   - y2fit(j,k);
            az = z2fit(j+1,k+1)   - z2fit(j,k);
            bx = x2fit(j,k+1,1)   - x2fit(j+1,k,1);
            by = y2fit(j,k+1)     - y2fit(j+1,k);
            bz = z2fit(j,k+1)     - z2fit(j+1,k);
            dnx = ay*bz - az*by;
            dny = az*bx - ax*bz;
            dnz = ax*by - ay*bx;
            d = std::sqrt(dnx*dnx + dny*dny + dnz*dnz);
            if ((float)d <= 0.f) d = 1.;
            xnfit = dnx/d;
            ynfit = dny/d;
            znfit = dnz/d;

            // normals at "to" cell centers calculated from "to" cell centers
            ax = x2int(j+1,k,1) - x2int(j,k,1);
            ay = y2int(j+1,k,1) - y2int(j,k,1);
            az = z2int(j+1,k,1) - z2int(j,k,1);
            bx = x2int(j,k+1,1) - x2int(j,k,1);
            by = y2int(j,k+1,1) - y2int(j,k,1);
            bz = z2int(j,k+1,1) - z2int(j,k,1);
            dnx = ay*bz - az*by;
            dny = az*bx - ax*bz;
            dnz = ax*by - ay*bx;
            d = std::sqrt(dnx*dnx + dny*dny + dnz*dnz);
            if ((float)d <= 0.f) d = 1.;
            xnint = dnx/d;
            ynint = dny/d;
            znint = dnz/d;

            // inner product of the two unit normals
            prod   = xnfit*xnint + ynfit*ynint + znfit*znint;
            prdsum = prdsum + prod;
            if ((float)prod < 1.f-(float)epsn || (float)prod > 1.f+(float)epsn) {
                nnorm = nnorm + 1;
            }
            dn = prod - 1.;
            if (std::abs((float)dn) > std::abs((float)dnmx)) {
                dnmx  = dn;
                jdnmx = j;
                kdnmx = k;
            }
        }} // end for j,k (seventh check)

        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "     %4d normals at interpolated cell centers differ by more than %3d percent",
            nnorm, (int)((float)epsn*100.f));
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "        from the normals at cell centers obtained directly from grid points");
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "       the max. difference %10.3e percent occurs at j,k = %4d,%4d",
            (float)(dnmx*100.), jdnmx, kdnmx);
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120,
            "       the average inner product over the interface is %10.3e",
            (float)prdsum / (float)(j22-j21-1) / (float)(k22-k21-1));

        if (nnorm > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "            possible causes:");
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "              ...mmcxie and/or mmceta incorrectly set");
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120,
                "              ...\"to\" and \"from\" blocks do not lie on the same surface");
        }
    } // end if jseg>1 && kseg>1 && iorph<=0 (seventh check)



    // ****** eighth check ******
    // Following output for visual check of the interpolation procedure.
    label_552:;
    if (iout > 0) {

        iunit = 40;

        if (icheck > 99) {
            len1 = 13;
            std::snprintf(titlptchgrd, 15, "patch_p3d.%3d", icheck);
        } else if (icheck > 9) {
            len1 = 12;
            std::snprintf(titlptchgrd, 15, "patch_p3d.%2d", icheck);
        } else {
            len1 = 11;
            std::snprintf(titlptchgrd, 15, "patch_p3d.%1d", icheck);
        }
        for (int i = len1+1; i <= 14; i++) {
            titlptchgrd[i-1] = ' ';
        }
        titlptchgrd[14] = '\0';

        // open(iunit, file=titlptchgrd(1:len1), form='formatted', status='unknown')
        {
            char fname[15];
            std::memcpy(fname, titlptchgrd, len1);
            fname[len1] = '\0';
            fortran_open_unit(iunit, fname, "w");
        }
        // rewind(iunit) - not needed after open for write

        if (iorph <= 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120,
                "   plot3d (/mg/for) file for this patch interface written to%14s",
                titlptchgrd);
        } else {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120,
                "   plot3d (/mg/for/blank) file for this patch interface written to%14s",
                titlptchgrd);
        }

        FILE* fu = fortran_get_unit(iunit);

        if (iorph <= 0) {
            fprintf(fu, " %d\n", lmax1+2);
        } else {
            fprintf(fu, " %d\n", lmax1+3);
        }

        // iwhole = 0: output original "from" grids
        iwhole = 0;
        maxx   = 0;

        if (iwhole == 0) {
            for (l = 1; l <= lmax1; l++) {
                if (xif1(l) == 1) xif1(l) = xif1(l) + 1;
                if (etf1(l) == 1) etf1(l) = etf1(l) + 1;
                if (xif2(l) == jjmax1(l)) xif2(l) = xif2(l) - 1;
                if (etf2(l) == kkmax1(l)) etf2(l) = etf2(l) - 1;
            }
        }

        if (istop != 1) {
            if (iorph <= 0) {
                // write grid dimensions for lmax1 "from" blocks + "to" grid + cell centers
                for (l = 1; l <= lmax1; l++) {
                    fprintf(fu, " %d %d %d",
                        xif2(l)-xif1(l)+1-maxx, etf2(l)-etf1(l)+1-maxx, 1);
                }
                fprintf(fu, " %d %d %d %d %d %d\n",
                    j22-j21+1, k22-k21+1, 1, j22-j21, k22-k21, 1);
            } else {
                for (l = 1; l <= lmax1; l++) {
                    fprintf(fu, " %d %d %d",
                        xif2(l)-xif1(l)+1-maxx, etf2(l)-etf1(l)+1-maxx, 1);
                }
                fprintf(fu, " %d %d %d %d %d %d %d %d %d\n",
                    j22-j21+1, k22-k21+1, 1, j22-j21, k22-k21, 1,
                    j22-j21, k22-k21, 1);
            }
        } else {
            for (l = 1; l <= lmax1; l++) {
                fprintf(fu, " %d %d %d",
                    xif2(l)-xif1(l)+1-maxx, etf2(l)-etf1(l)+1-maxx, 1);
            }
            fprintf(fu, " %d %d %d %d %d %d\n",
                j22-j21+1, k22-k21+1, 1, 1, 1, 1);
        }



        // write "from" grid coordinates
        for (l = 1; l <= lmax1; l++) {
            js = xif1(l);
            ks = etf1(l);
            je = xif2(l);
            ke = etf2(l);
            if (iorph <= 0) {
                if (ialph == 0) {
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)x1(j,k,l));
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)y1(j,k,l));
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)z1(j,k,l));
                } else {
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)x1(j,k,l));
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)z1(j,k,l));
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)(-y1(j,k,l)));
                }
            } else {
                if (ialph == 0) {
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)x1(j,k,l));
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)y1(j,k,l));
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)z1(j,k,l));
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %d", 1);
                } else {
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)x1(j,k,l));
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)z1(j,k,l));
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %g", (float)(-y1(j,k,l)));
                    for (k = ks; k <= ke; k++)
                        for (j = js; j <= je; j++)
                            fprintf(fu, " %d", 1);
                }
            }
            fprintf(fu, "\n");
        } // end for l (from grids)

        // write "to" grid coordinates
        if (iorph <= 0) {
            if (ialph == 0) {
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)x2(j,k,1));
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)y2(j,k,1));
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)z2(j,k,1));
            } else {
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)x2(j,k,1));
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)z2(j,k,1));
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)(-y2(j,k,1)));
            }
        } else {
            if (ialph == 0) {
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)x2(j,k,1));
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)y2(j,k,1));
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)z2(j,k,1));
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %d", 1);
            } else {
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)x2(j,k,1));
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)z2(j,k,1));
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %g", (float)(-y2(j,k,1)));
                for (k = k21; k <= k22; k++)
                    for (j = j21; j <= j22; j++)
                        fprintf(fu, " %d", 1);
            }
        }
        fprintf(fu, "\n");



        // write cell center data (istop==1: single point; else: interpolated centers)
        if (istop == 1) {
            if (ialph == 0) {
                fprintf(fu, " %g %g %g\n", (float)xc, (float)yc, (float)zc);
            } else {
                fprintf(fu, " %g %g %g\n", (float)xc, (float)zc, (float)(-yc));
            }
        } else {
            if (iorph <= 0) {
                if (ialph == 0) {
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)x2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)y2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)z2int(j,k,1));
                } else {
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)x2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)z2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)(-y2int(j,k,1)));
                }
                fprintf(fu, "\n");
            } else {
                // iorph > 0: first write blank array (1=interpolated, 0=orphan)
                for (j = j21; j <= j22-1; j++) {
                for (k = k21; k <= k22-1; k++) {
                    ll = (j22-j21)*(k-k21) + (j-j21+1);
                    l1 = mblkpt(ll);
                    x2fit(j,k,1) = 1.;
                    if (l1 == 0) x2fit(j,k,1) = 0.;
                }}
                if (ialph == 0) {
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)x2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)y2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)z2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %d", (int)x2fit(j,k,1));
                } else {
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)x2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)z2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)(-y2int(j,k,1)));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %d", (int)x2fit(j,k,1));
                }
                fprintf(fu, "\n");

                // second orphan block: set orphan cell centers from "to" grid average
                for (j = j21; j <= j22-1; j++) {
                for (k = k21; k <= k22-1; k++) {
                    ll = (j22-j21)*(k-k21) + (j-j21+1);
                    l1 = mblkpt(ll);
                    x2fit(j,k,1) = 0.;
                    if (l1 == 0) {
                        x2fit(j,k,1) = 1.;
                        x2int(j,k,1) = 0.25*(x2(j,k,1) + x2(j+1,k,1)
                                            + x2(j+1,k+1,1) + x2(j,k+1,1));
                        y2int(j,k,1) = 0.25*(y2(j,k,1) + y2(j+1,k,1)
                                            + y2(j+1,k+1,1) + y2(j,k+1,1));
                        z2int(j,k,1) = 0.25*(z2(j,k,1) + z2(j+1,k,1)
                                            + z2(j+1,k+1,1) + z2(j,k+1,1));
                    }
                }}
                if (ialph == 0) {
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)x2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)y2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)z2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %d", (int)x2fit(j,k,1));
                } else {
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)x2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)z2int(j,k,1));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %g", (float)(-y2int(j,k,1)));
                    for (k = k21; k <= k22-1; k++)
                        for (j = j21; j <= j22-1; j++)
                            fprintf(fu, " %d", (int)x2fit(j,k,1));
                }
                fprintf(fu, "\n");
            } // end iorph > 0
        } // end istop != 1

        fortran_close_unit(iunit);
        if (istop > 0) termn8_ns::termn8(myid, neg1, ibufdim, nbuf, bou, nou);

    } // end if iout > 0

    return;
}

} // namespace diagnos_ns
