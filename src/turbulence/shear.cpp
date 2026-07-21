// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "shear.h"
#include "extrae.h"
#include "direct.h"
#include "topol.h"
#include "newfit.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace shear_ns {

void shear(int& ishear, int& istop, int& iout, int& igap, int& jdim1, int& kdim1, int& msub1, int& msub2, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, int& lmax1, FortranArray3DRef<double> x1, FortranArray3DRef<double> y1, FortranArray3DRef<double> z1, FortranArray3DRef<double> x1mid, FortranArray3DRef<double> y1mid, FortranArray3DRef<double> z1mid, FortranArray3DRef<double> x1mide, FortranArray3DRef<double> y1mide, FortranArray3DRef<double> z1mide, int& limit0, FortranArray1DRef<int> jjmax2, FortranArray1DRef<int> kkmax2, FortranArray3DRef<double> x2, FortranArray3DRef<double> y2, FortranArray3DRef<double> z2, FortranArray1DRef<double> xie2, FortranArray1DRef<double> eta2, FortranArray1DRef<int> mblkpt, FortranArray1DRef<double> temp, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, int& ifit, int& itmax, double& xc, double& yc, double& zc, FortranArray3DRef<double> sxie2, FortranArray3DRef<double> seta2, int& jcorr, int& kcorr, int& intmx, int& icheck, FortranArray1DRef<int> nblkj, FortranArray1DRef<int> nblkk, FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm, int& mcxie, int& mceta, FortranArray1DRef<int> lout, int& j21, int& j22, int& k21, int& k22, int& npt, int& ic0, int& iorph, int& itoss0, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, int& iself, int& ifiner, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    // COMMON block aliases
    int32_t& isklt1 = cmn_sklt1.isklt1;
    float* ap       = cmn_areas.ap;   // ap[0]=ap(1), ap[1]=ap(2), ap[2]=ap(3)
    int32_t& imaxa  = cmn_areas.imaxa;
    int32_t& iretry = cmn_tacos.iretry;
    char* grdmov    = cmn_filenam.grdmov;

    // Local variables
    int l2, jmax2, kmax2, jl, jr, kl, kr;
    int j, k, ll, jcall, kcall, l1, l11, iok, icount, ncoin, mcoin;
    int js, ks, km, jm, kmaxck, jmaxck;
    int len1, i;
    double x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8;
    double a1, a2, a3;
    double xiet, etat, xiet1, etat1, etac, xiec;
    double fact, fact1, fact2, dis0, dis1, beta, red;
    int ifits, ifit1;
    int icase;
    char titlptchgrd[15];



    l2    = 1;
    jmax2 = jjmax2(l2);
    kmax2 = kkmax2(l2);
    jl = 1;
    jr = jmax2-1;
    kl = 1;
    kr = kmax2-1;

    if (isklt1 > 0) {
        nou(4) = std::min(nou(4)+1, ibufdim);
        std::snprintf(bou(nou(4),4), 120, "    beginning check of boundary values");
    }

    // ***************************************************************************
    //  correct boundary values near eta=0
    // ***************************************************************************

    kcorr = 0;
    ncoin = 0;

    if (mceta == 0) {
        if (isklt1 > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "      eta=0 boundaries not rendered coincident");
        }
        goto label670;
    }
    if (ishear >= 0) {
        if (isklt1 > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "      eta=0 boundaries being rendered coincident via shearing method");
        }
    } else {
        if (isklt1 > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "      eta=0 boundaries being rendered coincident via arc length method");
        }
    }

    // loop over all "to" cell on eta=0 boundary
    for (j = j21; j <= j22-1; j++) {

        // compute edge midpoints of first layer of "to" grid cells
        // along the eta=0 boundary using quadratic least squares
        jcall = j;
        kcall = k21;
        extrae_ns::extra(jdim1, kdim1, msub2, l2, x2, y2, z2,
                         jcall, kcall, jl, jr, x5, y5, z5, icase, ifit);
        if (j == j21) {
            extrae_ns::extrae(jdim1, kdim1, msub2, l2, x2, y2, z2,
                              jcall, kcall, kl, kr, x7, y7, z7, icase, ifit);
        }
        kcall = k21+1;
        extrae_ns::extra(jdim1, kdim1, msub2, l2, x2, y2, z2,
                         jcall, kcall, jl, jr, x6, y6, z6, icase, ifit);
        jcall = j+1;
        kcall = k21;
        extrae_ns::extrae(jdim1, kdim1, msub2, l2, x2, y2, z2,
                          jcall, kcall, kl, kr, x8, y8, z8, icase, ifit);

        // compute normalized directed areas/unit normals of "to" cell
        if (itoss0 == 0) {
            direct_ns::direct(x5, x6, x7, x8, y5, y6, y7, y8, z5, z6, z7, z8,
                              a1, a2, a3, imaxa, nou, bou, nbuf, ibufdim);
            ap[0] = (float)a1;
            ap[1] = (float)a2;
            ap[2] = (float)a3;
        }

        ifits = ifit;
        iretry = 0;
        icount = 1;

label18085:
        // compute midpoint of "to" cell on eta=0 boundary, consistent with ifit
        // bilinear or linear in xie
        if (ifit == 1 || ifit == 4) {
            xc = 0.5*(x2(j,k21,l2)+x2(j+1,k21,l2));
            yc = 0.5*(y2(j,k21,l2)+y2(j+1,k21,l2));
            zc = 0.5*(z2(j,k21,l2)+z2(j+1,k21,l2));
        } else {
            // biquadratic or quadratic in xie
            xc = x5;
            yc = y5;
            zc = z5;
        }

        // set starting point for search routine...use solution at cell center
        // of current cell
        ll = (j-j21+1);
        xiet = xie2(ll) + 1.;
        etat = eta2(ll) + 1.;
        l1   = mblkpt(ll);

        // search over "from" blocks to find xie,eta associated with "to"
        // cell midpoint xc,yc,zc on eta=0
        {
            int jto_val = j;
            int kto_val = 1;
            topol_ns::topol(jdim1, kdim1, msub1, jjmax1, kkmax1, lmax1, l1,
                            x1, y1, z1, x1mid, y1mid, z1mid, x1mide, y1mide, z1mide,
                            limit0, xc, yc, zc, xiet, etat, jimage, kimage,
                            ifit, itmax, igap, iok, lout, ic0, itoss0,
                            jto_val, kto_val, iself, xif1, xif2, etf1, etf2,
                            nou, bou, nbuf, ibufdim, myid);
        }


        // search routine unsuccessful...try alternative polynomial fit
        if (iretry == 0 && iok != 1) {
            if (icount < 4) {
                newfit_ns::newfit(ifits, ifit, icount);
                icount = icount + 1;
                goto label18085;
            } else {
                // don't shear orphan points
                if (iorph > 0) {
                    if (isklt1 > 0) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120, " quitting boundary correction because points have been flagged as orphans");
                    }
                    return;
                }

                if (ifiner > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        all attempts to find generalized coordinates");
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        of cell center j,k =%4d,%4d have been unsuccessful...", j, 1);
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        will use  averages of finer level: interpolation number%4d", ifiner);
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        for ALL points on this interface");
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
                    for (i = len1+1; i <= 14; i++) {
                        titlptchgrd[i-1] = ' ';
                    }
                    titlptchgrd[14] = '\0';
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " program terminated in dynamic patching routines - see file %-60.60s", grdmov);
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "         stopping...all attempts to find generalized coordinates");
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "         of cell midpoint j,k =  %d,1 have been unsuccessful...check%.*s",
                                  j, len1, titlptchgrd);
                }
                istop = 1;
                iout  = 1;
                return;
            }
        }

        // search routine successful with current polynomial fit
        if (iok == 1) {
            etat1 = etat;
            xiet1 = xiet;
            l11   = l1;
            ifit1 = ifit;
        }

        // try again with input polynomial fit
        if (ifit != ifits) {
            iretry = 1;
            ifit = ifits;
            goto label18085;
        }

        // if second try with input value of ifit not successful, use the ifit
        // value for which the search routine was successful
        if (iretry == 1) {
            if (iok != 1) {
                etat = etat1;
                xiet = xiet1;
                l1   = l11;
                if (isklt1 > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "         iterations using original fit  not successful at j,k= %d,%d", j, k21);
                    if (ifit == 1) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120, "          used bilinear fit instead");
                    }
                    if (ifit == 2) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120, "          used biquadratic fit instead");
                    }
                    if (ifit == 3) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120, "          used quadratic fit in xie, linear fit in eta instead");
                    }
                    if (ifit == 4) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120, "          used linear fit in xie, quadratic fit in eta instead");
                    }
                }
            }
        }

        // block locations for k=k21 boundary points are stored in nblkj(j)
        nblkj(j) = l1;
        ifit      = ifits;
        etac      = etat - 2.0;
        fact      = 1.;

        // For some cases, boundaries eta=0 on "to"
        // side and "from" side are not supposed to be coincident.
        mcoin = mceta*(j22-j21)/100;
        if ((float)etac > 1.) {
            ncoin = ncoin+1;
            if (ncoin >= mcoin) {
                if (isklt1 > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "      quitting check of boundary...");
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "      \"to\" and \"from\" eta=0 lines apparently not supposed to be coincident at this boundary");
                }
                goto label670;
            }
        }

        if ((float)etat < 2.0f) {
            // the following correction accounts for the widely varying cell
            // sizes near the boundary of a "from" grid after extension
            js    = (int)xiet;
            dis1  = (x1(js,3,l1)-x1(js,2,l1))*(x1(js,3,l1)-x1(js,2,l1))
                   +(y1(js,3,l1)-y1(js,2,l1))*(y1(js,3,l1)-y1(js,2,l1))
                   +(z1(js,3,l1)-z1(js,2,l1))*(z1(js,3,l1)-z1(js,2,l1));
            dis0  = (x1(js,2,l1)-x1(js,1,l1))*(x1(js,2,l1)-x1(js,1,l1))
                   +(y1(js,2,l1)-y1(js,1,l1))*(y1(js,2,l1)-y1(js,1,l1))
                   +(z1(js,2,l1)-z1(js,1,l1))*(z1(js,2,l1)-z1(js,1,l1));
            fact1 = std::sqrt(dis0/dis1);
            js    = js + 1;
            dis1  = (x1(js,3,l1)-x1(js,2,l1))*(x1(js,3,l1)-x1(js,2,l1))
                   +(y1(js,3,l1)-y1(js,2,l1))*(y1(js,3,l1)-y1(js,2,l1))
                   +(z1(js,3,l1)-z1(js,2,l1))*(z1(js,3,l1)-z1(js,2,l1));
            dis0  = (x1(js,2,l1)-x1(js,1,l1))*(x1(js,2,l1)-x1(js,1,l1))
                   +(y1(js,2,l1)-y1(js,1,l1))*(y1(js,2,l1)-y1(js,1,l1))
                   +(z1(js,2,l1)-z1(js,1,l1))*(z1(js,2,l1)-z1(js,1,l1));
            fact2 = std::sqrt(dis0/dis1);
            fact  = 0.5*(fact1 + fact2);
            etac  = etac*fact;
        }

        temp(j) = etac;

        for (k = k21; k <= k22-1; k++) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (mblkpt(ll) == l1) {
                seta2(j,k,1) = eta2(ll);
                if ((float)eta2(ll) < 1.) seta2(j,k,1) = (eta2(ll)-1.)*fact+1.;
            }
        }

        // find km, the last k-point on the current j-line in the "to" grid that
        // lies in the same "from" block as the k=1 boundary point
        kmm(j) = 0;
        for (k = k21; k <= k22-1; k++) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (mblkpt(ll) != nblkj(j)) goto label2111;
            kmm(j) = k;
        }
label2111:
        x7 = x8;
        y7 = y8;
        z7 = z8;
    } // end do 2000 j=j21,j22-1


    // shearing correction for eta
    if (ishear < 0) {
        kcorr = 1;
        goto label670;
    }

    beta = 25.;
    for (j = j21; j <= j22-1; j++) {
        kmaxck = kkmax1(nblkj(j));
        km     = kmm(j);
        if (km <= 1) continue;
        for (k = k21; k <= km; k++) {
            fact         = 1. - ((float)(km-k)-0.5f)/(float)(km-k21);
            red          = std::exp(-beta*fact*fact);
            seta2(j,k,1) = seta2(j,k,1) - temp(j)*red;

            // check to see if sheared generalized coordinates will remain in the
            // legal range...if not try the arc length correction near the surface
            if ((float)seta2(j,k,1) < 1. ||
                (float)seta2(j,k,1) > (float)(float(kmaxck-2))) {
                if (isklt1 > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        quitting eta shearing correction at j,k=%d,%d", j, k);
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        eta2(j,k),temp(j)= %g,%g",
                                  (float)seta2(j,k,1), (float)temp(j));
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        will try arc length correction near eta boundary ");
                }
                kcorr = 1;
                goto label670;
            }
        }

        for (k = k21; k <= k22-1; k++) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (mblkpt(ll) == nblkj(j)) eta2(ll) = seta2(j,k,1);
        }
    }

label670:

    // ***************************************************************************
    //  correct boundary values near xie=0
    // ***************************************************************************

    jcorr = 0;
    ncoin = 0;

    if (mcxie == 0) {
        if (isklt1 > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "      xie=0 boundaries not rendered coincident");
        }
        goto label770;
    }
    if (ishear >= 0) {
        if (isklt1 > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "      xie=0 boundaries being rendered coincident via shearing method");
        }
    } else {
        if (isklt1 > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120, "      xie=0 boundaries being rendered coincident via arc length method");
        }
    }

    // loop over all "to" cell on xie=0 boundary
    for (k = k21; k <= k22-1; k++) {

        // compute edge midpoints of first layer of "to" grid cells
        // along the xie=0 boundary using quadratic least squares
        kcall = k;
        jcall = j21;
        extrae_ns::extrae(jdim1, kdim1, msub2, l2, x2, y2, z2,
                          jcall, kcall, kl, kr, x7, y7, z7, icase, ifit);
        if (k == k21) {
            extrae_ns::extra(jdim1, kdim1, msub2, l2, x2, y2, z2,
                             jcall, kcall, jl, jr, x5, y5, z5, icase, ifit);
        }
        jcall = j21+1;
        extrae_ns::extrae(jdim1, kdim1, msub2, l2, x2, y2, z2,
                          jcall, kcall, kl, kr, x8, y8, z8, icase, ifit);
        jcall = j21;
        kcall = k+1;
        extrae_ns::extra(jdim1, kdim1, msub2, l2, x2, y2, z2,
                         jcall, kcall, jl, jr, x6, y6, z6, icase, ifit);

        // compute normalized directed areas/unit normals of "to" cell
        if (itoss0 == 0) {
            direct_ns::direct(x5, x6, x7, x8, y5, y6, y7, y8, z5, z6, z7, z8,
                              a1, a2, a3, imaxa, nou, bou, nbuf, ibufdim);
            ap[0] = (float)a1;
            ap[1] = (float)a2;
            ap[2] = (float)a3;
        }

        ifits = ifit;
        iretry = 0;
        icount = 1;

label19085:
        // compute midpoint of "to" cell on xie=0 boundary, consistent with ifit
        // bilinear or linear in eta
        if (ifit == 1 || ifit == 3) {
            xc = 0.5*(x2(j21,k+1,l2)+x2(j21,k,l2));
            yc = 0.5*(y2(j21,k+1,l2)+y2(j21,k,l2));
            zc = 0.5*(z2(j21,k+1,l2)+z2(j21,k,l2));
        } else {
            // biquadratic or quadratic in eta
            xc = x7;
            yc = y7;
            zc = z7;
        }

        // set starting point for search routine...use solution at cell center
        // of current cell
        ll = (j22-j21)*(k-k21) + 1;
        xiet = xie2(ll) + 1.;
        etat = eta2(ll) + 1.;
        l1   = mblkpt(ll);

        // search over "from" blocks to find xie,eta associated with "to"
        // cell midpoint xc,yc,zc on xie=0
        {
            int jto_val = 1;
            int kto_val = k;
            topol_ns::topol(jdim1, kdim1, msub1, jjmax1, kkmax1, lmax1, l1,
                            x1, y1, z1, x1mid, y1mid, z1mid, x1mide, y1mide, z1mide,
                            limit0, xc, yc, zc, xiet, etat, jimage, kimage,
                            ifit, itmax, igap, iok, lout, ic0, itoss0,
                            jto_val, kto_val, iself, xif1, xif2, etf1, etf2,
                            nou, bou, nbuf, ibufdim, myid);
        }


        // search routine unsuccessful...try alternative polynomial fit
        if (iretry == 0 && iok != 1) {
            if (icount < 4) {
                newfit_ns::newfit(ifits, ifit, icount);
                icount = icount + 1;
                goto label19085;
            } else {
                // don't shear orphan points
                if (iorph > 0) {
                    if (isklt1 > 0) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120, " quitting boundary correction because points have been flagged as orphans");
                    }
                    return;
                }

                if (ifiner > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        all attempts to find generalized coordinates");
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        of cell center j,k =%4d,%4d have been unsuccessful...", 1, k);
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        will use  averages of finer level: interpolation number%4d", ifiner);
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        for ALL points on this interface");
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
                    for (i = len1+1; i <= 14; i++) {
                        titlptchgrd[i-1] = ' ';
                    }
                    titlptchgrd[14] = '\0';
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " program terminated in dynamic patching routines - see file %-60.60s", grdmov);
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "         stopping...all attempts to find generalized coordinates");
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "         of cell center j,k =  1,%d,have been unsuccessful...check%.*s",
                                  k, len1, titlptchgrd);
                }
                istop = 1;
                iout  = 1;
                return;
            }
        }

        // search routine successful with current polynomial fit
        if (iok == 1) {
            xiet1 = xiet;
            etat1 = etat;
            l11   = l1;
            ifit1 = ifit;
        }

        // try again with input polynomial fit
        if (ifit != ifits) {
            iretry = 1;
            ifit = ifits;
            goto label19085;
        }

        // if second try with input value of ifit not successful, use the ifit
        // value for which the search routine was successful
        if (iretry == 1) {
            if (iok != 1) {
                xiet = xiet1;
                etat = etat1;
                l1   = l11;
                if (isklt1 > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "         iterations using original fit not successful at j,k= %d,%d", j21, k);
                    if (ifit == 1) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120, "          used bilinear fit instead");
                    }
                    if (ifit == 2) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120, "          used biquadratic fit instead");
                    }
                    if (ifit == 3) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120, "          used quadratic fit in xie, linear fit in eta instead");
                    }
                    if (ifit == 4) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4),4), 120, "          used linear fit in xie, quadratic fit in eta instead");
                    }
                }
            }
        }

        // block locations for j=j21 boundary points stored in nblkk(k)
        nblkk(k) = l1;
        ifit      = ifits;
        xiec      = xiet - 2.0;
        fact      = 1.;


        // For some cases, boundaries xie=0 on "to"
        // side and "from" side are not supposed to be coincident.
        mcoin = mcxie*(k22-k21)/100;
        if ((float)xiec > 1.) {
            ncoin = ncoin+1;
            if (ncoin >= mcoin) {
                if (isklt1 > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "      quitting check of boundary...");
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "      \"to\" and \"from\" xie=0 lines apparently not supposed to be coincident at this boundary");
                }
                goto label770;
            }
        }

        if ((float)xiet < 2.0f) {
            // the following correction accounts for the widely varying cell
            // sizes near the boundary of a "from" grid after extension
            ks    = (int)etat;
            dis1  = (x1(3,ks,l1)-x1(2,ks,l1))*(x1(3,ks,l1)-x1(2,ks,l1))
                   +(y1(3,ks,l1)-y1(2,ks,l1))*(y1(3,ks,l1)-y1(2,ks,l1))
                   +(z1(3,ks,l1)-z1(2,ks,l1))*(z1(3,ks,l1)-z1(2,ks,l1));
            dis0  = (x1(2,ks,l1)-x1(1,ks,l1))*(x1(2,ks,l1)-x1(1,ks,l1))
                   +(y1(2,ks,l1)-y1(1,ks,l1))*(y1(2,ks,l1)-y1(1,ks,l1))
                   +(z1(2,ks,l1)-z1(1,ks,l1))*(z1(2,ks,l1)-z1(1,ks,l1));
            fact1 = std::sqrt(dis0/dis1);
            ks    = ks + 1;
            dis1  = (x1(3,ks,l1)-x1(2,ks,l1))*(x1(3,ks,l1)-x1(2,ks,l1))
                   +(y1(3,ks,l1)-y1(2,ks,l1))*(y1(3,ks,l1)-y1(2,ks,l1))
                   +(z1(3,ks,l1)-z1(2,ks,l1))*(z1(3,ks,l1)-z1(2,ks,l1));
            dis0  = (x1(2,ks,l1)-x1(1,ks,l1))*(x1(2,ks,l1)-x1(1,ks,l1))
                   +(y1(2,ks,l1)-y1(1,ks,l1))*(y1(2,ks,l1)-y1(1,ks,l1))
                   +(z1(2,ks,l1)-z1(1,ks,l1))*(z1(2,ks,l1)-z1(1,ks,l1));
            fact2 = std::sqrt(dis0/dis1);
            fact  = 0.5*(fact1 + fact2);
            xiec  = xiec*fact;
        }

        temp(k) = xiec;

        for (j = j21; j <= j22-1; j++) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (mblkpt(ll) == l1) {
                sxie2(j,k,1) = xie2(ll);
                if ((float)xie2(ll) < (float)j21) sxie2(j,k,1) = (xie2(ll)-1.)*fact+1.;
            }
        }

        // find jm, the last j-point on the current k-line in the "to" grid that
        // lies in the same "from" block as the j=1 boundary point
        jmm(k) = 0;
        for (j = j21; j <= j22-1; j++) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (mblkpt(ll) != nblkk(k)) goto label4111;
            jmm(k) = j;
        }
label4111:
        x5 = x6;
        y5 = y6;
        z5 = z6;
    } // end do 4000 k=k21,k22-1

    // shearing correction for xie
    if (ishear < 0) {
        jcorr = 1;
        goto label770;
    }

    beta = 25.;
    for (k = k21; k <= k22-1; k++) {
        jmaxck = jjmax1(nblkk(k));
        jm     = jmm(k);
        if (jm <= j21) continue;
        for (j = j21; j <= jm; j++) {
            fact          = 1. - ((float)(jm-j)-0.5f)/(float)(jm-j21);
            red           = std::exp(-beta*fact*fact);
            sxie2(j,k,1)  = sxie2(j,k,1) - temp(k)*red;

            // check to see if sheared generalized coordinates will remain in the
            // legal range...if not try the arc length correction near the surface
            if ((float)sxie2(j,k,1) < 1. ||
                (float)sxie2(j,k,1) > (float)(jmaxck-2)) {
                if (isklt1 > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        quitting xie shearing correction at j,k=%d,%d", j, k);
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        xie2(j,k),temp(k)= %g,%g",
                                  (float)sxie2(j,k,1), (float)temp(k));
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    std::snprintf(bou(nou(4),4), 120, "        will try arc length correction near xie boundary");
                }
                jcorr = 1;
                goto label770;
            }
        }
        for (j = j21; j <= j22-1; j++) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (mblkpt(ll) == nblkk(k)) xie2(ll) = sxie2(j,k,1);
        }
    }

label770:
    return;
}

} // namespace shear_ns
