// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "patcher.h"
#include "termn8.h"
#include "loadgr.h"
#include "collapse.h"
#include "rechk.h"
#include "expand.h"
#include "transp.h"
#include "rotatp.h"
#include "avgint.h"
#include "invert.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace patcher_ns {

void patcher(int& nbl, FortranArray2DRef<int> lw, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& ncall, int& ioutpt, int& it_thro, int& maxbl, int& msub1, int& intmx, int& mxxe, int& mptch, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<double> windx, int& nintr, FortranArray2DRef<int> iindx, FortranArray1DRef<int> llimit, FortranArray1DRef<int> iitmax, FortranArray1DRef<int> mmcxie, FortranArray1DRef<int> mmceta, FortranArray1DRef<int> ncheck, FortranArray1DRef<int> iifit, FortranArray1DRef<int> mblkpt, FortranArray1DRef<int> iic0, FortranArray1DRef<int> iiorph, FortranArray1DRef<int> iitoss, FortranArray1DRef<int> ifiner, FortranArray2DRef<double> factjlo, FortranArray2DRef<double> factjhi, FortranArray2DRef<double> factklo, FortranArray2DRef<double> factkhi, FortranArray2DRef<double> dx, FortranArray2DRef<double> dy, FortranArray2DRef<double> dz, FortranArray2DRef<double> dthetx, FortranArray2DRef<double> dthety, FortranArray2DRef<double> dthetz, FortranArray2DRef<int> isav_dpat, FortranArray3DRef<int> isav_dpat_b, FortranArray3DRef<double> xte, FortranArray3DRef<double> yte, FortranArray3DRef<double> zte, FortranArray3DRef<double> xmi, FortranArray3DRef<double> ymi, FortranArray3DRef<double> zmi, FortranArray3DRef<double> xmie, FortranArray3DRef<double> ymie, FortranArray3DRef<double> zmie, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<int> jte, FortranArray1DRef<int> kte, FortranArray3DRef<double> sxie, FortranArray3DRef<double> seta, FortranArray2DRef<double> sxie2, FortranArray2DRef<double> seta2, FortranArray2DRef<double> xie2s, FortranArray2DRef<double> eta2s, FortranArray1DRef<double> temp, FortranArray2DRef<double> x2, FortranArray2DRef<double> y2, FortranArray2DRef<double> z2, FortranArray1DRef<int> nblk1, FortranArray1DRef<int> nblk2, FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm, FortranArray2DRef<double> x1, FortranArray2DRef<double> y1, FortranArray2DRef<double> z1, FortranArray1DRef<int> lout, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, FortranArray1DRef<int> ireq_ar, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    // COMMON block references
    float& epsc     = cmn_tol.epsc;
    float& epsc0    = cmn_tol.epsc0;
    float& epsreen  = cmn_tol.epsreen;
    float& epscoll  = cmn_tol.epscoll;
    int32_t& isklt1 = cmn_sklt1.isklt1;
    int32_t& iexp   = cmn_zero.iexp;
    int32_t& locv   = cmn_save.locv;
    int32_t& ireq   = cmn_save.ireq;
    int32_t& maxdcnt = cmn_is_dpatch.maxdcnt;

    // filenam COMMON aliases
    char* grid      = cmn_filenam.grid;
    char* plt3dg    = cmn_filenam.plt3dg;
    char* plt3dq    = cmn_filenam.plt3dq;
    char* output    = cmn_filenam.output;
    char* residual  = cmn_filenam.residual;
    char* turbres   = cmn_filenam.turbres;
    char* blomx     = cmn_filenam.blomx;
    char* output2   = cmn_filenam.output2;
    char* printout  = cmn_filenam.printout;
    char* pplunge   = cmn_filenam.pplunge;
    char* ovrlap    = cmn_filenam.ovrlap;
    char* patch     = cmn_filenam.patch;
    char* restrt    = cmn_filenam.restrt;
    char* subres    = cmn_filenam.subres;
    char* subtur    = cmn_filenam.subtur;
    char* grdmov    = cmn_filenam.grdmov;
    char* alphahist = cmn_filenam.alphahist;
    char* errfile   = cmn_filenam.errfile;
    char* preout    = cmn_filenam.preout;
    char* aeinp     = cmn_filenam.aeinp;
    char* aeout     = cmn_filenam.aeout;
    char* sdhist    = cmn_filenam.sdhist;
    char* avgg      = cmn_filenam.avgg;
    char* avgq      = cmn_filenam.avgq;

    // Local variables
    int iavg;
    double factj, factk;
    int xi1, xi2, et1, et2;
    int xif1l, xif2l, etf1l, etf2l;
    int xi1f, xi2f, et1f, et2f;
    // jjmax2 and kkmax2 are dimension(1) local arrays
    int jjmax2_arr[2], kkmax2_arr[2];
    FortranArray1DRef<int> jjmax2(jjmax2_arr + 1, 1);  // 1-based, size 1
    FortranArray1DRef<int> kkmax2(kkmax2_arr + 1, 1);  // 1-based, size 1

    int lcnt, nbll, nd_dest, lmax1, icheck, lst, npt, ifit;
    int limit0, itmax, mcxie, mceta, ic0, iorph, itoss0;
    int jmax2, kmax2, lx2, ly2, lz2, idim2g, jdim2g, kdim2g, jindex;
    int iself, mbl, lx1, ly1, lz1, jmax1, kmax1, idim1, jdim1, kdim1;
    int itest1, itest2, jjte, kkte;
    int icheckf, lmax1f, lstf, nptf;
    int l;



    // iavg flag for coarser level data
    // iavg = 0, use finer-level averages only if search on coarser level fails
    //        1, always use finer-level averages for coarser level interpolation data
    iavg = 0;

    // set tolerances, etc.
    // (10.**(-iexp) is machine zero)

    // expansion factors for "from" grids
    factj = 0.01;
    factk = 0.01;

    // convergence tolerance for generalized coordinates
    epsc = (float)std::max(1.e-07, std::pow(10.0, (double)(-iexp + 1)));

    // threshold for collapsed boundaries
    epscoll = (float)std::max(1.0e-10, std::pow(10.0, (double)(-iexp + 1)));

    // threshold for reentrant boundaries
    epsreen = (float)std::max(1.0e-09, std::pow(10.0, (double)(-iexp + 1)));

    // threshold for C-0 continuous interfaces
    epsc0 = (float)std::max(1.0e-07, std::pow(10.0, (double)(-iexp + 1)));

    isklt1 = 0;
    if (ioutpt > 0) isklt1 = 1;
    if (it_thro > 1) goto label_909;

label_909:
    ;

    for (lcnt = 1; lcnt <= maxdcnt; lcnt++) {
        // nbll is the current (to) block
        // mbl is the source (from) block
        // nbl is block coming into this routine
        nbll = isav_dpat(lcnt, 1);
        if (nbll == nbl) {
            nd_dest = mblk2nd(nbl);
            lmax1   = isav_dpat(lcnt, 2);
            icheck  = isav_dpat(lcnt, 17);
            lst     = iindx(icheck, 2*lmax1+5);
            npt     = iindx(icheck, 2*lmax1+4);
            ifit    = iifit(icheck);
            limit0  = llimit(icheck);
            itmax   = iitmax(icheck);
            mcxie   = mmcxie(icheck);
            mceta   = mmceta(icheck);
            ic0     = iic0(icheck);
            iorph   = iiorph(icheck);
            itoss0  = iitoss(icheck);
            xi1     = iindx(icheck, 2*lmax1+6);
            xi2     = iindx(icheck, 2*lmax1+7);
            et1     = iindx(icheck, 2*lmax1+8);
            et2     = iindx(icheck, 2*lmax1+9);
            if (mcxie > 100) mcxie = 200;
            if (mceta > 100) mceta = 200;

            if (isklt1 > 0) {
                nou(4) = std::min(nou(4)+1, ibufdim);
                std::snprintf(bou(nou(4), 4), 120, " ");
                nou(4) = std::min(nou(4)+1, ibufdim);
                std::snprintf(bou(nou(4), 4), 120, " ");
                nou(4) = std::min(nou(4)+1, ibufdim);
                // format 721: ' generalized coordinate interpolation number ',i3
                std::snprintf(bou(nou(4), 4), 120,
                    " generalized coordinate interpolation number %3d", icheck);
            }

            for (l = 1; l <= lmax1; l++) {
                mbl    = iindx(icheck, l+1);
                xif1l  = iindx(icheck, 2*lmax1+9+l);
                xif2l  = iindx(icheck, 3*lmax1+9+l);
                etf1l  = iindx(icheck, 4*lmax1+9+l);
                etf2l  = iindx(icheck, 5*lmax1+9+l);

                // patch surface to be interpolated from is an i=constant surface
                if (iindx(icheck, l+lmax1+2)/10 == 1) {
                    jjmax1(l) = jdimg(mbl);
                    kkmax1(l) = kdimg(mbl);
                    if (isklt1 > 0) {
                        if (iindx(icheck, l+lmax1+2) == 11) {
                            nou(4) = std::min(nou(4)+1, ibufdim);
                            // format 1957: ' interpolation from surface i =',i3,' of block',i3,2x,'(j: ',i3,' to ',i3,', k: ',i3,' to ',i3,')'
                            std::snprintf(bou(nou(4), 4), 120,
                                " interpolation from surface i =%3d of block%3d  (j: %3d to %3d, k: %3d to %3d)",
                                1, mbl, xif1l, xif2l, etf1l, etf2l);
                        } else {
                            nou(4) = std::min(nou(4)+1, ibufdim);
                            std::snprintf(bou(nou(4), 4), 120,
                                " interpolation from surface i =%3d of block%3d  (j: %3d to %3d, k: %3d to %3d)",
                                idimg(mbl), mbl, xif1l, xif2l, etf1l, etf2l);
                        }
                    }
                }

                // patch surface to be interpolated from is a j=constant surface
                if (iindx(icheck, l+lmax1+2)/10 == 2) {
                    jjmax1(l) = kdimg(mbl);
                    kkmax1(l) = idimg(mbl);
                    if (isklt1 > 0) {
                        if (iindx(icheck, l+lmax1+2) == 21) {
                            nou(4) = std::min(nou(4)+1, ibufdim);
                            // format 1958: ' interpolation from surface j =',i3,' of block',i3,2x,'(k: ',i3,' to ',i3,', i: ',i3,' to ',i3,')'
                            std::snprintf(bou(nou(4), 4), 120,
                                " interpolation from surface j =%3d of block%3d  (k: %3d to %3d, i: %3d to %3d)",
                                1, mbl, xif1l, xif2l, etf1l, etf2l);
                        } else {
                            nou(4) = std::min(nou(4)+1, ibufdim);
                            std::snprintf(bou(nou(4), 4), 120,
                                " interpolation from surface j =%3d of block%3d  (k: %3d to %3d, i: %3d to %3d)",
                                jdimg(mbl), mbl, xif1l, xif2l, etf1l, etf2l);
                        }
                    }
                }

                // patch surface to be interpolated from is a k=constant surface
                if (iindx(icheck, l+lmax1+2)/10 == 3) {
                    jjmax1(l) = jdimg(mbl);
                    kkmax1(l) = idimg(mbl);
                    if (isklt1 > 0) {
                        if (iindx(icheck, l+lmax1+2) == 31) {
                            nou(4) = std::min(nou(4)+1, ibufdim);
                            // format 1959: ' interpolation from surface k =',i3,' of block',i3,2x,'(j: ',i3,' to ',i3,', i: ',i3,' to ',i3,')'
                            std::snprintf(bou(nou(4), 4), 120,
                                " interpolation from surface k =%3d of block%3d  (j: %3d to %3d, i: %3d to %3d)",
                                1, mbl, xif1l, xif2l, etf1l, etf2l);
                        } else {
                            nou(4) = std::min(nou(4)+1, ibufdim);
                            std::snprintf(bou(nou(4), 4), 120,
                                " interpolation from surface k =%3d of block%3d  (j: %3d to %3d, i: %3d to %3d)",
                                kdimg(mbl), mbl, xif1l, xif2l, etf1l, etf2l);
                        }
                    }
                }

                if (jjmax1(l) > mptch || kkmax1(l) > mptch) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    // format: ' program terminated in dynamic patching routines - see file ',a60
                    std::snprintf(bou(nou(1), 1), 120,
                        " program terminated in dynamic patching routines - see file %-60.60s", grdmov);
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1), 1), 120,
                        " stopping ... mptch = %d  too small ", mptch);
                    int m1 = -1;
                    termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                }
            } // end do 1605



            // patch surface to be interpolated to is an i=constant surface
            if (iindx(icheck, 2*lmax1+3)/10 == 1) {
                jjmax2(1) = jdimg(nbl);
                kkmax2(1) = kdimg(nbl);
                if (isklt1 > 0) {
                    if (iindx(icheck, 2*lmax1+3) == 11) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        // format 1857: ' ',16x,'to surface i =',i3,' of block',i3,2x,'(j: ',i3,' to ',i3,', k: ',i3,' to ',i3,')'
                        std::snprintf(bou(nou(4), 4), 120,
                            "                  to surface i =%3d of block%3d  (j: %3d to %3d, k: %3d to %3d)",
                            1, nbl, xi1, xi2, et1, et2);
                    } else {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4), 4), 120,
                            "                  to surface i =%3d of block%3d  (j: %3d to %3d, k: %3d to %3d)",
                            idimg(nbl), nbl, xi1, xi2, et1, et2);
                    }
                }
            }

            // patch surface to be interpolated to is a j=constant surface
            if (iindx(icheck, 2*lmax1+3)/10 == 2) {
                jjmax2(1) = kdimg(nbl);
                kkmax2(1) = idimg(nbl);
                if (isklt1 > 0) {
                    if (iindx(icheck, 2*lmax1+3) == 21) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        // format 1858: ' ',16x,'to surface j =',i3,' of block',i3,2x,'(k: ',i3,' to ',i3,', i: ',i3,' to ',i3,')'
                        std::snprintf(bou(nou(4), 4), 120,
                            "                  to surface j =%3d of block%3d  (k: %3d to %3d, i: %3d to %3d)",
                            1, nbl, xi1, xi2, et1, et2);
                    } else {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4), 4), 120,
                            "                  to surface j =%3d of block%3d  (k: %3d to %3d, i: %3d to %3d)",
                            jdimg(nbl), nbl, xi1, xi2, et1, et2);
                    }
                }
            }

            // patch surface to be interpolated to is a k=constant surface
            if (iindx(icheck, 2*lmax1+3)/10 == 3) {
                jjmax2(1) = jdimg(nbl);
                kkmax2(1) = idimg(nbl);
                if (isklt1 > 0) {
                    if (iindx(icheck, 2*lmax1+3) == 31) {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        // format 1859: ' ',16x,'to surface k =',i3,' of block',i3,2x,'(j: ',i3,' to ',i3,', i: ',i3,' to ',i3,')'
                        std::snprintf(bou(nou(4), 4), 120,
                            "                  to surface k =%3d of block%3d  (j: %3d to %3d, i: %3d to %3d)",
                            1, nbl, xi1, xi2, et1, et2);
                    } else {
                        nou(4) = std::min(nou(4)+1, ibufdim);
                        std::snprintf(bou(nou(4), 4), 120,
                            "                  to surface k =%3d of block%3d  (j: %3d to %3d, i: %3d to %3d)",
                            kdimg(nbl), nbl, xi1, xi2, et1, et2);
                    }
                }
            }

            if (jjmax2(1) > mptch || kkmax2(1) > mptch) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " program terminated in dynamic patching routines - see file %-60.60s", grdmov);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " stopping ... mptch = %d  too small ", mptch);
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }

            if (isklt1 > 0) {
                nou(4) = std::min(nou(4)+1, ibufdim);
                std::snprintf(bou(nou(4), 4), 120, " ");
                nou(4) = std::min(nou(4)+1, ibufdim);
                // format 99: ' note: j and k referred to below are indicies local to the patch surface'
                std::snprintf(bou(nou(4), 4), 120,
                    " note: j and k referred to below are indicies local to the patch surface");
            }

            // set-up for "to" side of patch interface
            jmax2  = jjmax2(1);
            kmax2  = kkmax2(1);
            lx2    = lw(10, nbl);
            ly2    = lw(11, nbl);
            lz2    = lw(12, nbl);
            idim2g = idimg(nbl);
            jdim2g = jdimg(nbl);
            kdim2g = kdimg(nbl);
            jindex = iindx(icheck, 2*lmax1+3);

            // load proper grid from 1-d array into 2-d work array
            {
                int mdim_p2 = mptch + 2;
                loadgr_ns::loadgr(w, mgwk, lx2, ly2, lz2, jindex, x2, y2, z2,
                                  mdim_p2, mdim_p2, idim2g, jdim2g, kdim2g);
            }

            // check for collapsed grid lines
            if (isklt1 > 0) {
                nou(4) = std::min(nou(4)+1, ibufdim);
                // format 3958: ' ','   checking for collapsed grid lines on  "to"  side in block ',i3
                std::snprintf(bou(nou(4), 4), 120,
                    "    checking for collapsed grid lines on  \"to\"  side in block %3d", nbl);
            }

            {
                int mdim_p2 = mptch + 2;
                collapse_ns::collapse(mdim_p2, mdim_p2, jmax2, kmax2, x2, y2, z2,
                                      nou, bou, nbuf, ibufdim);
            }



            // set-up for "from" side of patch interface
            iself = 0;
            for (l = 1; l <= lmax1; l++) {
                mbl    = iindx(icheck, l+1);
                lx1    = lw(10, mbl);
                ly1    = lw(11, mbl);
                lz1    = lw(12, mbl);
                jmax1  = jjmax1(l);
                kmax1  = kkmax1(l);
                jindex = iindx(icheck, l+lmax1+2);
                idim1  = idimg(mbl);
                jdim1  = jdimg(mbl);
                kdim1  = kdimg(mbl);
                // iself = 1 if a block face communicates with itself
                itest1 = 100*nbl + iindx(icheck, 2*lmax1+3);
                itest2 = 100*mbl + iindx(icheck, l+lmax1+2);
                if (itest1 == itest2) iself = 1;

                // load proper grid from 1-d array into 2-d work array
                {
                    int mdim_p2 = mptch + 2;
                    loadgr_ns::loadgr(w, mgwk, lx1, ly1, lz1, jindex, x1, y1, z1,
                                      mdim_p2, mdim_p2, idim1, jdim1, kdim1);
                }

                // check for collapsed grid lines
                if (isklt1 > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    // format 3957: ' ','   checking for collapsed grid lines on  "from" side in block ',i3
                    std::snprintf(bou(nou(4), 4), 120,
                        "    checking for collapsed grid lines on \"from\" side in block %3d", mbl);
                }

                {
                    int mdim_p2 = mptch + 2;
                    collapse_ns::collapse(mdim_p2, mdim_p2, jmax1, kmax1, x1, y1, z1,
                                          nou, bou, nbuf, ibufdim);
                }

                // check for branch cuts
                if (isklt1 > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    // format 98: ' ','   checking for branch cuts on "from" side in block ',i3
                    std::snprintf(bou(nou(4), 4), 120,
                        "    checking for branch cuts on \"from\" side in block %3d", mbl);
                }

                xif1l = iindx(icheck, 2*lmax1+9+l);
                xif2l = iindx(icheck, 3*lmax1+9+l);
                etf1l = iindx(icheck, 4*lmax1+9+l);
                etf2l = iindx(icheck, 5*lmax1+9+l);
                {
                    int mdim_p2 = mptch + 2;
                    rechk_ns::rechk(mdim_p2, mdim_p2, jimage, kimage, msub1,
                                    jmax1, kmax1, l, x1, y1, z1,
                                    xif1l, xif2l, etf1l, etf2l,
                                    nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl);
                }

                // expand "from" grid(s) at boundaries to insure that the
                // "to" grid is completely covered
                if (isklt1 > 0) {
                    nou(4) = std::min(nou(4)+1, ibufdim);
                    // format 1755: ' ','   expanding grid boundaries on "from" side in block ',i3
                    std::snprintf(bou(nou(4), 4), 120,
                        "    expanding grid boundaries on \"from\" side in block %3d", mbl);
                }

                // for now, set lo/hi values for factj/factk to single
                // factj/factk values set above...later may add the
                // ability to read individual values on input
                factjlo(icheck, l) = factj;
                factjhi(icheck, l) = factj;
                factklo(icheck, l) = factk;
                factkhi(icheck, l) = factk;

                // should not really need any significant expansion
                // for c-0 grids
                if (ic0 > 0) {
                    factjlo(icheck, l) = 1.e-5;
                    factjhi(icheck, l) = 1.e-5;
                    factklo(icheck, l) = 1.e-5;
                    factkhi(icheck, l) = 1.e-5;
                }

                {
                    int mdim_p2 = mptch + 2;
                    expand_ns::expand(mdim_p2, mdim_p2, msub1, jmax1, kmax1, l,
                                      x1, y1, z1, xte, yte, zte,
                                      factjlo(icheck, l), factjhi(icheck, l),
                                      factklo(icheck, l), factkhi(icheck, l),
                                      jmax2, kmax2, x2, y2, z2);
                }
                jte(l) = jjmax1(l) + 2;
                kte(l) = kkmax1(l) + 2;

                // translate/rotate "from" blocks as needed to provide
                // sufficient coverage on "from" side when only part
                // of the physical domain is modeled
                if (std::abs((float)dx(icheck, l)) > 0. ||
                    std::abs((float)dy(icheck, l)) > 0. ||
                    std::abs((float)dz(icheck, l)) > 0.) {
                    jjte = jte(l);
                    kkte = kte(l);
                    {
                        int mdim_p2 = mptch + 2;
                        transp_ns::transp(mdim_p2, mdim_p2, jjte, kkte, msub1, l,
                                          xte, yte, zte, dx, dy, dz, intmx, icheck);
                    }
                }
                if (std::abs((float)dthetx(icheck, l)) > 0. ||
                    std::abs((float)dthety(icheck, l)) > 0. ||
                    std::abs((float)dthetz(icheck, l)) > 0.) {
                    jjte = jte(l);
                    kkte = kte(l);
                    {
                        int mdim_p2 = mptch + 2;
                        rotatp_ns::rotatp(mdim_p2, mdim_p2, jjte, kkte, msub1, l,
                                          xte, yte, zte, dthetx, dthety, dthetz,
                                          xorig, yorig, zorig, mbl, maxbl, intmx, icheck);
                    }
                }

                // search range on "from" side
                xif1(l) = iindx(icheck, 2*lmax1+9+l);
                if (xif1(l) > 1) {
                    xif1(l) = xif1(l) + 1;
                }
                xif2(l) = iindx(icheck, 3*lmax1+9+l);
                if (xif2(l) == jjmax1(l)) {
                    xif2(l) = xif2(l) + 2;
                } else {
                    xif2(l) = xif2(l) + 1;
                }
                etf1(l) = iindx(icheck, 4*lmax1+9+l);
                if (etf1(l) > 1) {
                    etf1(l) = etf1(l) + 1;
                }
                etf2(l) = iindx(icheck, 5*lmax1+9+l) + 1;
                if (etf2(l) == kkmax1(l)) {
                    etf2(l) = etf2(l) + 2;
                } else {
                    etf2(l) = etf2(l) + 1;
                }
            } // end do 125



            if (ifiner(icheck) != 0) {
                icheckf = ifiner(icheck);
                lmax1f  = iindx(icheckf, 1);
                lstf    = iindx(icheckf, 2*lmax1f+5);
                nptf    = iindx(icheckf, 2*lmax1f+4);
                xi1f    = iindx(icheckf, 2*lmax1f+6);
                xi2f    = iindx(icheckf, 2*lmax1f+7);
                et1f    = iindx(icheckf, 2*lmax1f+8);
                et2f    = iindx(icheckf, 2*lmax1f+9);
            } else {
                lstf    = iindx(icheck, 2*lmax1+5);
                nptf    = iindx(icheck, 2*lmax1+4);
                xi1f    = iindx(icheck, 2*lmax1+6);
                xi2f    = iindx(icheck, 2*lmax1+7);
                et1f    = iindx(icheck, 2*lmax1+8);
                et2f    = iindx(icheck, 2*lmax1+9);
            }

            // obtain interpolation coefficients for coarser levels
            // by averaging finer-level coefficients
            if (ifiner(icheck) != 0 && iavg > 0) {
                // windx(lst,1) → FortranArray1DRef from &windx(lst,1)
                // windx(lst,2) → FortranArray1DRef from &windx(lst,2)
                // mblkpt(lst)  → FortranArray1DRef from &mblkpt(lst)
                // windx(lstf,1) → FortranArray1DRef from &windx(lstf,1)
                // windx(lstf,2) → FortranArray1DRef from &windx(lstf,2)
                // mblkpt(lstf) → FortranArray1DRef from &mblkpt(lstf)
                // Compute sizes: npt and nptf elements from lst/lstf
                int windx_size1 = windx.ubound(1) - lst + 1;
                int windx_size2 = windx.ubound(1) - lstf + 1;
                int mblkpt_size1 = mblkpt.ubound(1) - lst + 1;
                int mblkpt_size2 = mblkpt.ubound(1) - lstf + 1;
                FortranArray1DRef<double> windx_lst1(&windx(lst, 1), windx_size1);
                FortranArray1DRef<double> windx_lst2(&windx(lst, 2), windx_size1);
                FortranArray1DRef<int>    mblkpt_lst(&mblkpt(lst), mblkpt_size1);
                FortranArray1DRef<double> windx_lstf1(&windx(lstf, 1), windx_size2);
                FortranArray1DRef<double> windx_lstf2(&windx(lstf, 2), windx_size2);
                FortranArray1DRef<int>    mblkpt_lstf(&mblkpt(lstf), mblkpt_size2);
                avgint_ns::avgint(windx_lst1, windx_lst2, mblkpt_lst, npt,
                                  windx_lstf1, windx_lstf2, mblkpt_lstf, nptf,
                                  xi1, xi2, et1, et2, xi1f, xi2f, et1f, et2f);
            }


            // determine projection of x2,y2,z2 points onto generalized
            // coordinate system(s) defined by the grid(s) on "from" side
            {
                // Build FortranArray1DRef views for windx columns at lst and lstf
                int windx_n1 = windx.ubound(1) - lst + 1;
                int windx_nf = windx.ubound(1) - lstf + 1;
                int mblkpt_n1 = mblkpt.ubound(1) - lst + 1;
                int mblkpt_nf = mblkpt.ubound(1) - lstf + 1;
                FortranArray1DRef<double> windx_lst1(&windx(lst, 1), windx_n1);
                FortranArray1DRef<double> windx_lst2(&windx(lst, 2), windx_n1);
                FortranArray1DRef<int>    mblkpt_lst(&mblkpt(lst), mblkpt_n1);
                FortranArray1DRef<double> windx_lstf1(&windx(lstf, 1), windx_nf);
                FortranArray1DRef<double> windx_lstf2(&windx(lstf, 2), windx_nf);
                FortranArray1DRef<int>    mblkpt_lstf(&mblkpt(lstf), mblkpt_nf);

                // invert expects sxie2/seta2 as FortranArray3DRef but patcher has them as 2DRef
                // wrap 2D arrays as 3D with last dim = 1
                int sxie2_n1 = sxie2.ubound(1);
                int sxie2_n2 = sxie2.ubound(2);
                FortranArray3DRef<double> sxie2_3d(&sxie2(1, 1), sxie2_n1, sxie2_n2, 1);
                FortranArray3DRef<double> seta2_3d(&seta2(1, 1), sxie2_n1, sxie2_n2, 1);

                // invert expects x2/y2/z2 as FortranArray3DRef but patcher has them as 2DRef
                // wrap 2D arrays as 3D with last dim = 1
                int x2_n1 = x2.ubound(1);
                int x2_n2 = x2.ubound(2);
                FortranArray3DRef<double> x2_3d(&x2(1, 1), x2_n1, x2_n2, 1);
                FortranArray3DRef<double> y2_3d(&y2(1, 1), x2_n1, x2_n2, 1);
                FortranArray3DRef<double> z2_3d(&z2(1, 1), x2_n1, x2_n2, 1);

                // ifiner(icheck) must be passed as int& (not rvalue)
                int ifiner_icheck = ifiner(icheck);

                int mdim_p2 = mptch + 2;
                int msub2_val = 1;
                invert_ns::invert(mdim_p2, mdim_p2, msub1, msub2_val,
                                  jte, kte,
                                  lmax1, xte, yte, zte, xmi, ymi, zmi, xmie, ymie, zmie,
                                  limit0, jjmax2, kkmax2, x2_3d, y2_3d, z2_3d,
                                  windx_lst1, windx_lst2, mblkpt_lst,
                                  temp, jimage, kimage,
                                  ifit, itmax, sxie, seta, sxie2_3d, seta2_3d,
                                  xie2s, eta2s,
                                  intmx, icheck, nblk1, nblk2, jmm, kmm,
                                  mcxie, mceta,
                                  lout, xi1, xi2, et1, et2, npt,
                                  ic0, iorph, itoss0, ncall,
                                  ioutpt, xif1, xif2, etf1, etf2, iself,
                                  ifiner_icheck,
                                  windx_lstf1, windx_lstf2, mblkpt_lstf, nptf,
                                  xi1f, xi2f, et1f, et2f, iavg,
                                  nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl);
            }

        } // end if (nbll == nbl)
    } // end do lcnt

    return;
}

} // namespace patcher_ns
