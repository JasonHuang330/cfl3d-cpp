// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "invert.h"
#include "extra.h"
#include "extrae.h"
#include "direct.h"
#include "topol.h"
#include "newfit.h"
#include "shear.h"
#include "arc.h"
#include "avgint.h"
#include "diagnos.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace invert_ns {

void invert(int& jdim1, int& kdim1, int& msub1, int& msub2, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, int& lmax1, FortranArray3DRef<double> x1, FortranArray3DRef<double> y1, FortranArray3DRef<double> z1, FortranArray3DRef<double> x1mid, FortranArray3DRef<double> y1mid, FortranArray3DRef<double> z1mid, FortranArray3DRef<double> x1mide, FortranArray3DRef<double> y1mide, FortranArray3DRef<double> z1mide, int& limit0, FortranArray1DRef<int> jjmax2, FortranArray1DRef<int> kkmax2, FortranArray3DRef<double> x2, FortranArray3DRef<double> y2, FortranArray3DRef<double> z2, FortranArray1DRef<double> xie2, FortranArray1DRef<double> eta2, FortranArray1DRef<int> mblkpt, FortranArray1DRef<double> temp, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, int& ifit, int& itmax, FortranArray3DRef<double> sxie, FortranArray3DRef<double> seta, FortranArray3DRef<double> sxie2, FortranArray3DRef<double> seta2, FortranArray2DRef<double> xie2s, FortranArray2DRef<double> eta2s, int& intmx, int& icheck, FortranArray1DRef<int> nblkj, FortranArray1DRef<int> nblkk, FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm, int& mcxie, int& mceta, FortranArray1DRef<int> lout, int& j21, int& j22, int& k21, int& k22, int& npt, int& ic0, int& iorph, int& itoss0, int& ncall, int& ioutpt, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, int& iself, int& ifiner, FortranArray1DRef<double> xie2f, FortranArray1DRef<double> eta2f, FortranArray1DRef<int> mblkptf, int& nptf, int& xi1f, int& xi2f, int& et1f, int& et2f, int& iavg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    // COMMON block aliases
    int32_t& isklt1 = cmn_sklt1.isklt1;
    float* ap       = cmn_areas.ap;   // ap[0]=ap(1), ap[1]=ap(2), ap[2]=ap(3)
    int32_t& imaxa  = cmn_areas.imaxa;
    int32_t& iretry = cmn_tacos.iretry;
    int32_t& itrace = cmn_tracer.itrace;
    // /filenam/ - only grdmov is used here
    char* grdmov    = cmn_filenam.grdmov;

    // Local variables
    int l1, l2, jmax2, kmax2, istop, iout, igap;
    int l, jmax1, kmax1, j, k;
    int jl, jr, kl, kr;
    int icase, iok;
    int nblsav, ll, ifit1;
    int icount, iorfn;
    int ifits;
    int kcall, jcall;
    int ishear, jcorr, kcorr;
    int len1;
    int idum1, idum2, idum3, idum4;
    double dum1, dum2, dum3;
    double x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8;
    double xc, yc, zc;
    double xiet, etat;
    double a1, a2, a3;
    char titlptchgrd[15];


    l2    = 1;
    jmax2 = jjmax2(l2);
    kmax2 = kkmax2(l2);
    istop = 0;
    iout  = 0;
    igap  = 0;
    if (ifit < 0) iout = 1;
    ifit  = std::abs(ifit);

    // itrace < 0, do not write search history for current "to" cell
    // itrace = 0, overwrite history from previous "to" cell with current
    // itrace = 1, retain the search history for ALL cells (may get huge file)
    // trace output found in fort.7
    itrace = -1;

    idum1 = 0;
    idum2 = 0;
    idum3 = 0;
    idum4 = 0;
    dum1  = 0.;
    dum2  = 0.;
    dum3  = 0.;

    if (isklt1 > 0) {
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120, "    beginning inversion for xie, eta of cell centers");
        if (iavg == 0 || ifiner == 0) {
            if (ifit == 1) {
                nou(4) = std::min(nou(4)+1, ibufdim);
                std::snprintf(bou(nou(4),4), 120, "      using bi-linear fit");
            }
            if (ifit == 2) {
                nou(4) = std::min(nou(4)+1, ibufdim);
                std::snprintf(bou(nou(4),4), 120, "      using bi-quadratic fit");
            }
            if (ifit == 3) {
                nou(4) = std::min(nou(4)+1, ibufdim);
                std::snprintf(bou(nou(4),4), 120, "      using quadratic fit in xie, linear fit in eta");
            }
            if (ifit == 4) {
                nou(4) = std::min(nou(4)+1, ibufdim);
                std::snprintf(bou(nou(4),4), 120, "      using linear fit in xie, quadratic fit in eta");
            }
        } else if (ifiner > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "     using averages of finer level: interpolation number%4d", ifiner);
        }
    }

    // call trace(0,icheck,idum2,idum3,idum4,dum1,dum2,dum3)

    // compute cell edge midpoints of "from" grid(s)
    // using quadratic least squares
    for (l = 1; l <= lmax1; l++) {
        jmax1 = jjmax1(l);
        kmax1 = kkmax1(l);
        for (j = 1; j <= jmax1-1; j++) {
            for (k = 1; k <= kmax1; k++) {
                jl = 2;
                jr = jmax1-2;
                extra_ns::extra(jdim1, kdim1, msub1, l, x1, y1, z1,
                                j, k, jl, jr, x5, y5, z5, icase, ifit);
                x1mid(j,k,l) = x5;
                y1mid(j,k,l) = y5;
                z1mid(j,k,l) = z5;
            }
        }
        for (k = 1; k <= kmax1-1; k++) {
            for (j = 1; j <= jmax1; j++) {
                kl = 2;
                kr = kmax1-2;
                extrae_ns::extrae(jdim1, kdim1, msub1, l, x1, y1, z1,
                                  j, k, kl, kr, x7, y7, z7, icase, ifit);
                x1mide(j,k,l) = x7;
                y1mide(j,k,l) = y7;
                z1mide(j,k,l) = z7;
            }
        }
    }

    // don't go through inversion process if coarser level
    // interpolation data are always obtained by averaging
    // fine level data
    if (ifiner > 0 && iavg > 0) goto lbl_552;

    nblsav = 1;
    jl     = 1;
    jr     = jmax2-1;
    kl     = 1;
    kr     = kmax2-1;

    // loop over all "to" cells
    for (j = j21; j <= j22-1; j++) {
    for (k = k22-1; k >= k21; k--) {

        // compute cell edge midpoints of "to" grid
        // using quadratic least squares
        if (k == k22-1) {
            kcall = k+1;
            extra_ns::extra(jdim1, kdim1, msub2, l2, x2, y2, z2,
                            j, kcall, jl, jr, x6, y6, z6, icase, ifit);
        }
        extra_ns::extra(jdim1, kdim1, msub2, l2, x2, y2, z2,
                        j, k, jl, jr, x5, y5, z5, icase, ifit);
        extrae_ns::extrae(jdim1, kdim1, msub2, l2, x2, y2, z2,
                          j, k, kl, kr, x7, y7, z7, icase, ifit);
        jcall = j+1;
        extrae_ns::extrae(jdim1, kdim1, msub2, l2, x2, y2, z2,
                          jcall, k, kl, kr, x8, y8, z8, icase, ifit);

        // compute normalized directed areas/unit normals of "to" cell
        if (itoss0 == 0) {
            direct_ns::direct(x5, x6, x7, x8, y5, y6, y7, y8, z5, z6, z7, z8,
                              a1, a2, a3, imaxa, nou, bou, nbuf, ibufdim);
            ap[0] = (float)a1;
            ap[1] = (float)a2;
            ap[2] = (float)a3;
        }

        ifits  = ifit;
        iretry = 0;
        icount = 1;
        iorfn  = 0;


lbl_17085:
        // compute center of "to" cell, consistent with ifit
        // bi-linear
        if (ifit == 1) {
            xc = 0.25*( x2(j,k,l2) + x2(j+1,k,l2)
                       + x2(j+1,k+1,l2) + x2(j,k+1,l2) );
            yc = 0.25*( y2(j,k,l2) + y2(j+1,k,l2)
                       + y2(j+1,k+1,l2) + y2(j,k+1,l2) );
            zc = 0.25*( z2(j,k,l2) + z2(j+1,k,l2)
                       + z2(j+1,k+1,l2) + z2(j,k+1,l2) );
        }
        // bi-quadratic
        if (ifit == 2) {
            xc = 0.5* ( x5 + x6 + x7 + x8 )
               - 0.25*( x2(j,k,l2)     + x2(j+1,k,l2)
               +        x2(j+1,k+1,l2) + x2(j,k+1,l2) );
            yc = 0.5* ( y5 + y6 + y7 + y8 )
               - 0.25*( y2(j,k,l2)    + y2(j+1,k,l2)
               +        y2(j+1,k+1,l2) + y2(j,k+1,l2) );
            zc = 0.5* ( z5 + z6 + z7 + z8 )
               - 0.25*( z2(j,k,l2)     + z2(j+1,k,l2)
               +        z2(j+1,k+1,l2) + z2(j,k+1,l2) );
        }
        // quadratic in xie, linear in eta
        if (ifit == 3) {
            xc = .5*(x5 + x6);
            yc = .5*(y5 + y6);
            zc = .5*(z5 + z6);
        }
        // linear in xie, quadratic in eta
        if (ifit == 4) {
            xc = .5*(x7 + x8);
            yc = .5*(y7 + y8);
            zc = .5*(z7 + z8);
        }

        // call trace(1,icheck,j,k,ifit,xc,yc,zc)

        // set starting point for search routine
        // first time.........use solution at the previous point, or search
        //                    for minimum distance point if xie=eta=0.
        // subsequent times...use solution at the same point from the
        //                    previous time.
        //
        // xiet,etat are generalized coordinates in expanded "from" grid
        // xie2,eta2 are generalized coordinates in original "from" grid
        if (iretry == 0) {
            if (ncall == 1) {
                if (k < k22-1) {
                    ll = (j22-j21)*(k+1-k21) + (j-j21+1);
                    if (mblkpt(ll) != 0) {
                        xiet = xie2(ll) + 1.;
                        etat = eta2(ll) + 1.;
                        l1   = mblkpt(ll);
                    } else {
                        xiet = 0.;
                        etat = 0.;
                        l1   = nblsav;
                        l1   = std::max(l1, 1);
                    }
                } else if (j == j21) {
                    xiet = 0.;
                    etat = 0.;
                    l1   = 1;
                } else if (j > j21) {
                    ll = (j22-j21)*(k-k21) + (j-1-j21+1);
                    if (mblkpt(ll) != 0) {
                        xiet = xie2(ll) + 1.;
                        etat = eta2(ll) + 1.;
                        l1   = mblkpt(ll);
                    } else {
                        xiet = 0.;
                        etat = 0.;
                        l1   = nblsav;
                        l1   = std::max(l1, 1);
                    }
                }
            } else {
                ll   = (j22-j21)*(k-k21) + (j-j21+1);
                xiet = xie2(ll) + 1.;
                etat = eta2(ll) + 1.;
                l1   = mblkpt(ll);
                l1   = std::max(l1, 1);
            }
        } else {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (mblkpt(ll) != 0) {
                xiet = xie2(ll) + 1.;
                etat = eta2(ll) + 1.;
                l1   = mblkpt(ll);
            } else {
                xiet = 0.;
                etat = 0.;
                l1   = nblsav;
                l1   = std::max(l1, 1);
            }
        }


lbl_17086:
        // search over "from" blocks to find xie,eta associated with "to"
        // cell center xc,yc,zc.
        topol_ns::topol(jdim1, kdim1, msub1, jjmax1, kkmax1, lmax1, l1,
                        x1, y1, z1, x1mid, y1mid, z1mid, x1mide, y1mide, z1mide,
                        limit0, xc, yc, zc, xiet, etat, jimage, kimage, ifit,
                        itmax, igap, iok, lout, ic0, itoss0, j, k, iself,
                        xif1, xif2, etf1, etf2, nou, bou, nbuf, ibufdim, myid);

        // search routine unsuccessful...try an alternative polynomial fit
        if (iretry == 0 && iok != 1) {
            if (icount < 4) {
                newfit_ns::newfit(ifits, ifit, icount);
                icount = icount + 1;
                goto lbl_17085;
            } else {
                if (iorph <= 0) {
                    if (ifiner > 0) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120,
                            "        all attempts to find generalized coordinates");
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120,
                            "        of cell center j,k =%4d,%4d have been unsuccessful...", j, k);
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120,
                            "        will use  averages of finer level: interpolation number%4d", ifiner);
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120,
                            "        for ALL points on this interface");
                        istop = 1;
                        goto lbl_552;
                    } else {
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
                        for (int i = len1; i < 14; i++) {
                            titlptchgrd[i] = ' ';
                        }
                        titlptchgrd[14] = '\0';
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " program terminated in dynamic patching routines - see file %-60.60s",
                            grdmov);
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120,
                            "        stopping...all attempts to find generalized coordinates of");
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120,
                            "        cell center j,k =%4d,%4d have been unsuccessful...check%14.14s",
                            j, k, titlptchgrd);
                        istop = 1;
                        iout  = 1;
                        goto lbl_552;
                    }
                } else {
                    // mark orphan point as interpolated from block "0" - will
                    // not be interpolated in cfl3d - this option used for points
                    // really do lie outside the domain of the from block - for
                    // example, an exposed flap surface.
                    iok = 1;
                    ll = (j22-j21)*(k-k21) + (j-j21+1);
                    nblsav     = mblkpt(ll);
                    xie2(ll)   = 1.;
                    eta2(ll)   = 1.;
                    mblkpt(ll) = 0;
                    x6 = x5;
                    y6 = y5;
                    z6 = z5;
                    // call trace(50,j,k,idum3,idum4,dum1,dum2,dum3)
                    goto lbl_1001;
                }
            }
        }

        // search successful with current polynomial fit
        if (iok == 1) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            xie2(ll)   = xiet - 1.;
            eta2(ll)   = etat - 1.;
            mblkpt(ll) = l1;
            ifit1      = ifit;
        }


        if (iorph > 0) {
            // search routine "succesful", but only because it found what might
            // otherwise be an orphan point in the expanded "from" grid.
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if ((float)xie2(ll) < 1.f || (float)xie2(ll) > (float)(jjmax1(l1)-2) ||
                (float)eta2(ll) < 1.f || (float)eta2(ll) > (float)(kkmax1(l1)-2)) {
                iorfn = iorfn + 1;
                if (iorfn <= 1) {
                    // redo the search for this point, starting with a minimum
                    // distance search to make sure this point should really
                    // be marked as an orphan
                    iok    = 0;
                    icount = 1;
                    iretry = 0;
                    xiet   = 0.;
                    etat   = 0.;
                    goto lbl_17086;
                } else {
                    // mark the point as an orphan
                    iok = 1;
                    nblsav     = mblkpt(ll);
                    xie2(ll)   = 1.;
                    eta2(ll)   = 1.;
                    mblkpt(ll) = 0;
                    x6 = x5;
                    y6 = y5;
                    z6 = z5;
                    // call trace(50,j,k,idum3,idum4,dum1,dum2,dum3)
                    goto lbl_1001;
                }
            }
        }

        // if current polynomial fit is not the one which was input,
        // try again with the input value
        if (ifit != ifits) {
            iretry = 1;
            ifit   = ifits;
            goto lbl_17085;
        }

        // if second try with input value of ifit not successful, use the ifit
        // value for which the search routine was successful.
        if (iretry == 1) {
            if (iok != 1) {
                if (isklt1 > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120,
                        "         iterations using original fit not successful at j,k=  %d %d",
                        j, k);
                    if (ifit == 1) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120,
                            "          used bilinear fit instead");
                    }
                    if (ifit == 2) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120,
                            "          used biquadratic fit instead");
                    }
                    if (ifit == 3) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120,
                            "          used quadratic fit in xie, linear fit in eta instead");
                    }
                    if (ifit == 4) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120,
                            "          used linear fit in xie, quadratic fit in eta instead");
                    }
                }
            }
        }

        x6 = x5;
        y6 = y5;
        z6 = z5;

lbl_1001:;
    } // end k loop
    } // end j loop


    // check and, if necessary, correct values of xie and/or eta
    // near boundaries.
    //
    // ishear > 0: try a shearing correction to make the "to"
    //             and "from" boundaries coincident; if shearing fails, then
    //             try an arc length correction.
    // ishear < 0: use only the arc length correction
    ishear = 1;

    shear_ns::shear(ishear, istop, iout, igap, jdim1, kdim1, msub1,
                    msub2, jjmax1, kkmax1, lmax1, x1, y1, z1, x1mid, y1mid,
                    z1mid, x1mide, y1mide, z1mide, limit0, jjmax2, kkmax2,
                    x2, y2, z2, xie2, eta2, mblkpt, temp, jimage, kimage, ifit,
                    itmax, xc, yc, zc, sxie2, seta2, jcorr, kcorr, intmx, icheck,
                    nblkj, nblkk, jmm, kmm, mcxie, mceta, lout, j21, j22, k21, k22,
                    npt, ic0, iorph, itoss0, xif1, xif2, etf1, etf2, iself, ifiner,
                    nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl);

    if (jcorr != 0 || kcorr != 0)
        arc_ns::arc(jdim1, kdim1, msub1, msub2, jjmax1, kkmax1,
                    lmax1, x1, y1, z1, limit0, jjmax2, kkmax2, x2, y2, z2, xie2,
                    eta2, mblkpt, ifit, itmax, jcorr, kcorr, sxie, seta, sxie2,
                    seta2, xie2s, eta2s, intmx, icheck, nblkj, nblkk, jmm, kmm,
                    j21, j22, k21, k22, npt, xif1, xif2, etf1, etf2,
                    nou, bou, nbuf, ibufdim, mblk2nd, maxbl);

lbl_552:;

    // for coarse level interfaces where the standard search
    // algorithim has failed, use average of finer level data
    if (ifiner > 0 && istop == 1) {
        avgint_ns::avgint(xie2, eta2, mblkpt, npt, xie2f, eta2f, mblkptf,
                          nptf, j21, j22, k21, k22, xi1f, xi2f, et1f, et2f);
        istop = 0;
    }

    // perform diagnostic checks on xie,eta values found for "to" cell
    // centers
    if (ioutpt > 0 || istop > 0) {
        diagnos_ns::diagnos(istop, iout, igap, jdim1, kdim1, msub1,
                            msub2, jjmax1, kkmax1, lmax1, x1, y1, z1, x1mid, y1mid,
                            z1mid, x1mide, y1mide, z1mide, sxie, seta, sxie2, seta2,
                            xie2s, eta2s, jjmax2, kkmax2, x2, y2,
                            z2, xie2, eta2, mblkpt, icheck, intmx, xc, yc, zc, ifit,
                            j21, j22, k21, k22, npt, ic0, iorph, xif1, xif2,
                            etf1, etf2, itoss0, iself,
                            nou, bou, nbuf, ibufdim, myid);
    }
    return;
}

} // namespace invert_ns
