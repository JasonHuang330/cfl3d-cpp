// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "findmin_new.h"
#include "lead.h"
#include "ccomplex.h"
#include "termn8.h"
#include "init.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <algorithm>

namespace findmin_new_ns {

using namespace ccomplex_ns;

// Local helper: isrcheq(n, arr, inc, val)
// Searches arr(1..n) for first element equal to val, returns 1-based index
static int isrcheq(int n, FortranArray1DRef<double> arr, int /*inc*/, double val) {
    for (int i = 1; i <= n; i++) {
        if (arr(i) == val) return i;
    }
    return 1;
}

// Helper to write into bou(nou(1),1) using snprintf
// bou is FortranArray2DRef<char[120]>, nou is FortranArray1DRef<int>
// Usage: write_bou(bou, nou, ibufdim, fmt, ...)
// We use a macro-like inline approach

void findmin_new(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
                 FortranArray1DRef<double> w, int& mgwk,
                 FortranArray1DRef<double> wk, int& nwork,
                 FortranArray1DRef<int> iwk, int& iwork,
                 int& nsurf, int& j1,
                 FortranArray1DRef<double> xs, FortranArray1DRef<int> ixs,
                 FortranArray1DRef<int> lsminn,
                 FortranArray1DRef<int> ireq_xs, FortranArray1DRef<int> ireq_bb,
                 int& ngrid,
                 FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nbci0,
                 FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
                 FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim,
                 FortranArray1DRef<int> nbckdim,
                 FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
                 FortranArray4DRef<int> ibcinfo,
                 FortranArray1DRef<int> nblg, FortranArray1DRef<int> nou,
                 FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
                 int& maxbl, int& maxgr, int& maxseg,
                 FortranArray1DRef<int> mblk2nd)
{
    // COMMON block references
    int& jdim  = cmn_ginfo.jdim;
    int& kdim  = cmn_ginfo.kdim;
    int& idim  = cmn_ginfo.idim;
    int& jj2   = cmn_ginfo.jj2;
    int& kk2   = cmn_ginfo.kk2;
    int& ii2   = cmn_ginfo.ii2;
    int& nblc  = cmn_ginfo.nblc;
    int& lq    = cmn_ginfo.lq;
    int& lqj0  = cmn_ginfo.lqj0;
    int& lqk0  = cmn_ginfo.lqk0;
    int& lqi0  = cmn_ginfo.lqi0;
    int& lsj   = cmn_ginfo.lsj;
    int& lsk   = cmn_ginfo.lsk;
    int& lsi   = cmn_ginfo.lsi;
    int& lvol  = cmn_ginfo.lvol;
    int& ldtj  = cmn_ginfo.ldtj;
    int& lx    = cmn_ginfo.lx;
    int& ly    = cmn_ginfo.ly;
    int& lz    = cmn_ginfo.lz;
    int& lvis  = cmn_ginfo.lvis;
    int& lsnk0 = cmn_ginfo.lsnk0;
    int& lsni0 = cmn_ginfo.lsni0;
    int& lq1   = cmn_ginfo.lq1;
    int& lqr   = cmn_ginfo.lqr;
    int& lblk  = cmn_ginfo.lblk;
    int& lxib  = cmn_ginfo.lxib;
    int& lsig  = cmn_ginfo.lsig;
    int& lsqtq = cmn_ginfo.lsqtq;
    int& lg    = cmn_ginfo.lg;
    int& ltj0  = cmn_ginfo.ltj0;
    int& ltk0  = cmn_ginfo.ltk0;
    int& lti0  = cmn_ginfo.lti0;
    int& lxkb  = cmn_ginfo.lxkb;
    int& lnbl  = cmn_ginfo.lnbl;
    int& lvj0  = cmn_ginfo.lvj0;
    int& lvk0  = cmn_ginfo.lvk0;
    int& lvi0  = cmn_ginfo.lvi0;
    int& lbcj  = cmn_ginfo.lbcj;
    int& lbck  = cmn_ginfo.lbck;
    int& lbci  = cmn_ginfo.lbci;
    int& lqc0  = cmn_ginfo.lqc0;
    int& ldqc0 = cmn_ginfo.ldqc0;
    int& lxtbi = cmn_ginfo.lxtbi;
    int& lxtbj = cmn_ginfo.lxtbj;
    int& lxtbk = cmn_ginfo.lxtbk;
    int& latbi = cmn_ginfo.latbi;
    int& latbj = cmn_ginfo.latbj;
    int& latbk = cmn_ginfo.latbk;
    int& lbcdj = cmn_ginfo.lbcdj;
    int& lbcdk = cmn_ginfo.lbcdk;
    int& lbcdi = cmn_ginfo.lbcdi;
    int& lxib2 = cmn_ginfo.lxib2;
    int& lux   = cmn_ginfo.lux;
    int& lcmuv = cmn_ginfo.lcmuv;
    int& lvolj0= cmn_ginfo.lvolj0;
    int& lvolk0= cmn_ginfo.lvolk0;
    int& lvoli0= cmn_ginfo.lvoli0;
    int& lxmdj = cmn_ginfo.lxmdj;
    int& lxmdk = cmn_ginfo.lxmdk;
    int& lxmdi = cmn_ginfo.lxmdi;
    int& lvelg = cmn_ginfo.lvelg;
    int& ldeltj= cmn_ginfo.ldeltj;
    int& ldeltk= cmn_ginfo.ldeltk;
    int& ldelti= cmn_ginfo.ldelti;
    int& lxnm2 = cmn_ginfo.lxnm2;
    int& lynm2 = cmn_ginfo.lynm2;
    int& lznm2 = cmn_ginfo.lznm2;
    int& lxnm1 = cmn_ginfo.lxnm1;
    int& lynm1 = cmn_ginfo.lynm1;
    int& lznm1 = cmn_ginfo.lznm1;
    int& lqavg = cmn_ginfo.lqavg;
    int& ivmx  = cmn_maxiv.ivmx;
    int& myid  = cmn_mydist2.myid;

    // Local variables
    int surf, ntri, iptri, isurf, surfi;
    int igrid, nbl, m, iseg, ns, nface;
    int nbctype, n1beg, n1end, n2beg, n2end;
    int npts, nvalxs, nvalbb;
    int minbox, nbb, ntotv;
    int bbdef, ipv, vlist;
    int lsmin, lsminc;
    int ncg;
    int jdimc, kdimc, idimc;
    int ierrflg_m1 = -1;



    // Build xs/ixs as 2D views for passing to getpts/getptsbb/collect_surf etc.
    // xs(4,nsurf), ixs(4,nsurf)
    FortranArray2DRef<double> xs2d(&xs(1), 4, nsurf);
    FortranArray2DRef<int>    ixs2d(&ixs(1), 4, nsurf);

    // memory allocation based on nsurf
    initi(iwork);
    initf(nwork);
    { int tmp = 4*nsurf; surf  = ifalloc(nou, bou, nbuf, ibufdim, myid, tmp); }
    ntri  = iialloc(nou, bou, nbuf, ibufdim, myid, nsurf);
    { int tmp = 8*nsurf; iptri = iialloc(nou, bou, nbuf, ibufdim, myid, tmp); }
    isurf = 0;
    surfi = 0;
    if (ivmx == 4 || ivmx == 25)
        { int tmp = 4*nsurf; surfi = iialloc(nou, bou, nbuf, ibufdim, myid, tmp); }

    // collect all the viscous surface points into an array xsurf(nsurf,3)
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nblg(igrid);
        lead_ns::lead(nbl, lw, lw2, maxbl);

        // i=constant surfaces
        for (m = 1; m <= 2; m++) {
            if (m == 1) { ns = nbci0(nbl);   nface = 1; }
            else         { ns = nbcidim(nbl); nface = 2; }
            for (iseg = 1; iseg <= ns; iseg++) {
                nbctype = ibcinfo(nbl, iseg, 1, m);
                n1beg   = ibcinfo(nbl, iseg, 2, m);
                n1end   = ibcinfo(nbl, iseg, 3, m);
                n2beg   = ibcinfo(nbl, iseg, 4, m);
                n2end   = ibcinfo(nbl, iseg, 5, m);
                npts    = (n1end - n1beg + 1) * (n2end - n2beg + 1);
                nvalxs  = 4 * npts;
                nvalbb  = 4 * npts;
                if (std::abs(nbctype) == 2004 || std::abs(nbctype) == 2014 ||
                    std::abs(nbctype) == 2024 || std::abs(nbctype) == 2034 ||
                    std::abs(nbctype) == 2016) {
                    if (myid == mblk2nd(nbl)) {
                        FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                        FortranArray3DRef<double> wblank(&w(lblk), jdim, kdim, idim);
                        FortranArray3DRef<double> wbci(&w(lbci), jdim, kdim, 2);
                        FortranArray3DRef<double> wbcj(&w(lbcj), kdim, idim-1, 2);
                        FortranArray3DRef<double> wbck(&w(lbck), jdim, idim-1, 2);
                        int imn=1, jmn=1, kmn=1;
                        if (ivmx == 4 || ivmx == 25) {
                            getptsbb(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                     wx, wy, wz, nface, n1beg, n1end, n2beg, n2end,
                                     xs2d, ixs2d, nsurf, nbl, wblank, wbci, wbcj, wbck);
                        } else {
                            getpts(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                   wx, wy, wz, nface, n1beg, n1end, n2beg, n2end,
                                   xs2d, nsurf, nbl, wblank, wbci, wbcj, wbck);
                        }
                    }
                    {
                        FortranArray2DRef<double> wksurf(&wk(surf), nsurf, 4);
                        FortranArray1DRef<int>    iwkntri(&iwk(ntri), nsurf);
                        FortranArray2DRef<int>    iwkiptri(&iwk(iptri), nsurf, 8);
                        int imn=1, jmn=1, kmn=1;
                        if (ivmx == 4 || ivmx == 25) {
                            FortranArray2DRef<int> iwksurfi(&iwk(surfi), nsurf, 4);
                            collect_surfbb(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                           xs2d, ixs2d, nface, n1beg, n1end, n2beg, n2end,
                                           nsurf, wksurf, isurf, iwkntri, iwkiptri, iwksurfi, nbl);
                        } else {
                            collect_surf(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                         xs2d, nface, n1beg, n1end, n2beg, n2end,
                                         nsurf, wksurf, isurf, iwkntri, iwkiptri);
                        }
                    }
                }
            }
        }


        // j=constant surfaces
        for (m = 1; m <= 2; m++) {
            if (m == 1) { ns = nbcj0(nbl);   nface = 3; }
            else         { ns = nbcjdim(nbl); nface = 4; }
            for (iseg = 1; iseg <= ns; iseg++) {
                nbctype = jbcinfo(nbl, iseg, 1, m);
                n1beg   = jbcinfo(nbl, iseg, 4, m);
                n1end   = jbcinfo(nbl, iseg, 5, m);
                n2beg   = jbcinfo(nbl, iseg, 2, m);
                n2end   = jbcinfo(nbl, iseg, 3, m);
                npts    = (n1end - n1beg + 1) * (n2end - n2beg + 1);
                nvalxs  = 4 * npts;
                nvalbb  = 4 * npts;
                if (std::abs(nbctype) == 2004 || std::abs(nbctype) == 2014 ||
                    std::abs(nbctype) == 2024 || std::abs(nbctype) == 2034 ||
                    std::abs(nbctype) == 2016) {
                    if (myid == mblk2nd(nbl)) {
                        FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                        FortranArray3DRef<double> wblank(&w(lblk), jdim, kdim, idim);
                        FortranArray3DRef<double> wbci(&w(lbci), jdim, kdim, 2);
                        FortranArray3DRef<double> wbcj(&w(lbcj), kdim, idim-1, 2);
                        FortranArray3DRef<double> wbck(&w(lbck), jdim, idim-1, 2);
                        int imn=1, jmn=1, kmn=1;
                        if (ivmx == 4 || ivmx == 25) {
                            getptsbb(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                     wx, wy, wz, nface, n1beg, n1end, n2beg, n2end,
                                     xs2d, ixs2d, nsurf, nbl, wblank, wbci, wbcj, wbck);
                        } else {
                            getpts(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                   wx, wy, wz, nface, n1beg, n1end, n2beg, n2end,
                                   xs2d, nsurf, nbl, wblank, wbci, wbcj, wbck);
                        }
                    }
                    {
                        FortranArray2DRef<double> wksurf(&wk(surf), nsurf, 4);
                        FortranArray1DRef<int>    iwkntri(&iwk(ntri), nsurf);
                        FortranArray2DRef<int>    iwkiptri(&iwk(iptri), nsurf, 8);
                        int imn=1, jmn=1, kmn=1;
                        if (ivmx == 4 || ivmx == 25) {
                            FortranArray2DRef<int> iwksurfi(&iwk(surfi), nsurf, 4);
                            collect_surfbb(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                           xs2d, ixs2d, nface, n1beg, n1end, n2beg, n2end,
                                           nsurf, wksurf, isurf, iwkntri, iwkiptri, iwksurfi, nbl);
                        } else {
                            collect_surf(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                         xs2d, nface, n1beg, n1end, n2beg, n2end,
                                         nsurf, wksurf, isurf, iwkntri, iwkiptri);
                        }
                    }
                }
            }
        }

        // k=constant surfaces
        for (m = 1; m <= 2; m++) {
            if (m == 1) { ns = nbck0(nbl);   nface = 5; }
            else         { ns = nbckdim(nbl); nface = 6; }
            for (iseg = 1; iseg <= ns; iseg++) {
                nbctype = kbcinfo(nbl, iseg, 1, m);
                n1beg   = kbcinfo(nbl, iseg, 2, m);
                n1end   = kbcinfo(nbl, iseg, 3, m);
                n2beg   = kbcinfo(nbl, iseg, 4, m);
                n2end   = kbcinfo(nbl, iseg, 5, m);
                npts    = (n1end - n1beg + 1) * (n2end - n2beg + 1);
                nvalxs  = 4 * npts;
                nvalbb  = 4 * npts;
                if (std::abs(nbctype) == 2004 || std::abs(nbctype) == 2014 ||
                    std::abs(nbctype) == 2024 || std::abs(nbctype) == 2034 ||
                    std::abs(nbctype) == 2016) {
                    if (myid == mblk2nd(nbl)) {
                        FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                        FortranArray3DRef<double> wblank(&w(lblk), jdim, kdim, idim);
                        FortranArray3DRef<double> wbci(&w(lbci), jdim, kdim, 2);
                        FortranArray3DRef<double> wbcj(&w(lbcj), kdim, idim-1, 2);
                        FortranArray3DRef<double> wbck(&w(lbck), jdim, idim-1, 2);
                        int imn=1, jmn=1, kmn=1;
                        if (ivmx == 4 || ivmx == 25) {
                            getptsbb(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                     wx, wy, wz, nface, n1beg, n1end, n2beg, n2end,
                                     xs2d, ixs2d, nsurf, nbl, wblank, wbci, wbcj, wbck);
                        } else {
                            getpts(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                   wx, wy, wz, nface, n1beg, n1end, n2beg, n2end,
                                   xs2d, nsurf, nbl, wblank, wbci, wbcj, wbck);
                        }
                    }
                    {
                        FortranArray2DRef<double> wksurf(&wk(surf), nsurf, 4);
                        FortranArray1DRef<int>    iwkntri(&iwk(ntri), nsurf);
                        FortranArray2DRef<int>    iwkiptri(&iwk(iptri), nsurf, 8);
                        int imn=1, jmn=1, kmn=1;
                        if (ivmx == 4 || ivmx == 25) {
                            FortranArray2DRef<int> iwksurfi(&iwk(surfi), nsurf, 4);
                            collect_surfbb(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                           xs2d, ixs2d, nface, n1beg, n1end, n2beg, n2end,
                                           nsurf, wksurf, isurf, iwkntri, iwkiptri, iwksurfi, nbl);
                        } else {
                            collect_surf(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                         xs2d, nface, n1beg, n1end, n2beg, n2end,
                                         nsurf, wksurf, isurf, iwkntri, iwkiptri);
                        }
                    }
                }
            }
        }
    } // end igrid loop


    // determine minimum distance
    // sort with respect to x coordinate
    {
        FortranArray2DRef<double> wksurf(&wk(surf), nsurf, 4);
        FortranArray1DRef<int>    iwkntri(&iwk(ntri), nsurf);
        FortranArray2DRef<int>    iwkiptri(&iwk(iptri), nsurf, 8);
        if (ivmx == 4 || ivmx == 25) {
            FortranArray2DRef<int> iwksurfi(&iwk(surfi), nsurf, 4);
            sort_xbb(nsurf, wksurf, iwkntri, iwkiptri, iwksurfi,
                     wk, iwk, nou, bou, nbuf, ibufdim, myid);
        } else {
            sort_x(nsurf, wksurf, iwkntri, iwkiptri,
                   wk, iwk, nou, bou, nbuf, ibufdim, myid);
        }
    }

    // run Mac Ice routines for making boxes
    minbox = (int)std::sqrt((float)nsurf);
    minbox = std::max(minbox, 50);
    nbb = 3 * nsurf / minbox;

    { int tmp = 6*nbb; bbdef  = ifalloc(nou, bou, nbuf, ibufdim, myid, tmp); }
    { int tmp = 2*nbb; ipv    = iialloc(nou, bou, nbuf, ibufdim, myid, tmp); }
    vlist  = iialloc(nou, bou, nbuf, ibufdim, myid, nsurf);
    {
        FortranArray2DRef<double> wkbbdef(&wk(bbdef), 6, nbb);
        FortranArray2DRef<int>    iwkipv(&iwk(ipv), 2, nbb);
        FortranArray1DRef<int>    iwkvlist(&iwk(vlist), nsurf);
        FortranArray1DRef<double> wksurf_x(&wk(surf), nsurf);
        FortranArray1DRef<double> wksurf_y(&wk(surf + nsurf), nsurf);
        FortranArray1DRef<double> wksurf_z(&wk(surf + 2*nsurf), nsurf);
        int maxlv = 10;
        makebb(minbox, maxlv, nsurf, wksurf_x, wksurf_y, wksurf_z,
               nbb, ntotv, wkbbdef, iwkipv, iwkvlist, iwk,
               nou, bou, nbuf, ibufdim, myid);
    }


    // calculate distance from each point in the field to points on the viscous surfaces
    for (igrid = 1; igrid <= ngrid; igrid++) {
        // fine level
        nbl = nblg(igrid) + j1;
        lead_ns::lead(nbl, lw, lw2, maxbl);
        if (myid == mblk2nd(nbl)) {
            { int tmp = jdim*kdim*idim; lsmin = ifalloc(nou, bou, nbuf, ibufdim, myid, tmp); }
            lsminn(nbl) = lsmin;
            nou(1) = std::min(nou(1)+1, ibufdim);
            // format: 45h computing smin (min dist func) for field eqn, 21h turb model for block, i6
            std::snprintf(bou(nou(1),1), 120,
                "computing smin (min dist func) for field eqn turb model for block%6d", nbl);
            {
                FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                FortranArray3DRef<double> wklsmin(&wk(lsmin), jdim, kdim, idim);
                FortranArray2DRef<double> wksurf(&wk(surf), nsurf, 4);
                FortranArray2DRef<double> wkbbdef(&wk(bbdef), 6, nbb);
                FortranArray2DRef<int>    iwkipv(&iwk(ipv), 2, nbb);
                FortranArray1DRef<int>    iwkvlist(&iwk(vlist), nsurf);
                FortranArray1DRef<int>    iwkntri(&iwk(ntri), nsurf);
                FortranArray2DRef<int>    iwkiptri(&iwk(iptri), nsurf, 8);
                int imn=1, jmn=1, kmn=1;
                if (ivmx == 4 || ivmx == 25) {
                    FortranArray2DRef<int>    iwksurfi(&iwk(surfi), nsurf, 4);
                    FortranArray3DRef<double> wxkb(&w(lxkb), jdim-1, kdim-1, idim-1);
                    FortranArray4DRef<double> wxib(&w(lxib), jdim, kdim, idim, 2);
                    FortranArray3DRef<double> wnbl(&w(lnbl), jdim-1, kdim-1, idim-1);
                    FortranArray3DRef<double> wsni0(&w(lsni0), jdim-1, kdim-1, idim-1);
                    calc_distbb(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                                idim, jdim, kdim,
                                wx, wy, wz, wklsmin,
                                nsurf, wksurf, nbb, wkbbdef, iwkipv,
                                iwkvlist, iwkntri, iwkiptri, iwksurfi,
                                wsni0, wxkb, wxib, wnbl,
                                wk, iwk, j1,
                                nou, bou, nbuf, ibufdim, myid);
                } else {
                    calc_dist(imn, jmn, kmn, idim, jdim, kdim, idim, jdim, kdim,
                              idim, jdim, kdim,
                              wx, wy, wz, wklsmin,
                              nsurf, wksurf, nbb, wkbbdef, iwkipv,
                              iwkvlist, iwkntri, iwkiptri,
                              wk, iwk,
                              nou, bou, nbuf, ibufdim, myid);
                }
            }
        }

        // coarser levels
        ncg = ncgg(igrid) - j1;
        if (ncg > 0) {
            for (m = 1; m <= ncg; m++) {
                nbl = nbl + 1;
                if (myid == mblk2nd(nbl)) {
                    { int tmp = jj2*kk2*ii2; lsminc = ifalloc(nou, bou, nbuf, ibufdim, myid, tmp); }
                    lsminn(nbl) = lsminc;
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "computing smin (min dist func) for field eqn turb model for block%6d", nbl);
                    {
                        FortranArray3DRef<double> wklsmin(&wk(lsmin), jdim, kdim, idim);
                        FortranArray3DRef<double> wklsminc(&wk(lsminc), jj2, kk2, ii2);
                        distcg(jdim, kdim, idim, wklsmin, jj2, kk2, ii2, wklsminc);
                    }
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    lsmin = lsminc;
                }
            }
        }


        if (ncg > 0 && (ivmx == 4 || ivmx == 25)) {
            int nblc_loc = nblg(igrid) + j1;
            for (m = 1; m <= ncg; m++) {
                lead_ns::lead(nblc_loc, lw, lw2, maxbl);
                nbl = nblc_loc;
                nblc_loc = nblc_loc + 1;
                if (myid == mblk2nd(nbl)) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "getting bbarth data on coarser block %d", nblc_loc);
                    jdimc = jdim/2 + 1;
                    kdimc = kdim/2 + 1;
                    idimc = idim/2 + 1;
                    {
                        FortranArray3DRef<double> wxjb_nbl(&w(lw(15,nbl)), jdim-1, kdim-1, idim-1);
                        FortranArray3DRef<double> wxkb_nbl(&w(lw(26,nbl)), jdim-1, kdim-1, idim-1);
                        FortranArray4DRef<double> wxib_nbl(&w(lw(19,nbl)), jdim, kdim, idim, 2);
                        FortranArray3DRef<double> wblnum_nbl(&w(lw(27,nbl)), jdim-1, kdim-1, idim-1);
                        FortranArray3DRef<double> wxjbc(&w(lw(15,nblc_loc)), jdimc-1, kdimc-1, idimc-1);
                        FortranArray3DRef<double> wxkbc(&w(lw(26,nblc_loc)), jdimc-1, kdimc-1, idimc-1);
                        FortranArray4DRef<double> wxibc(&w(lw(19,nblc_loc)), jdimc, kdimc, idimc, 2);
                        FortranArray3DRef<double> wblnumc(&w(lw(27,nblc_loc)), jdimc-1, kdimc-1, idimc-1);
                        bbarthcg(jdim, kdim, idim, jdimc, kdimc, idimc,
                                 wxjb_nbl, wxkb_nbl, wxib_nbl, wblnum_nbl,
                                 wxjbc, wxkbc, wxibc, wblnumc);
                    }
                }
            }
        }
    } // end igrid loop (distance calculation)

    // generate distance to cell-centers and store in permanent array w
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nblg(igrid) + j1;
        if (myid == mblk2nd(nbl)) {
            lead_ns::lead(nbl, lw, lw2, maxbl);
            lsmin = lsminn(nbl);
            {
                FortranArray3DRef<double> wklsmin(&wk(lsmin), jdim, kdim, idim);
                FortranArray3DRef<double> wsnk0(&w(lsnk0), jdim-1, kdim-1, idim-1);
                distcc(jdim, kdim, idim, wklsmin, wsnk0);
            }
            ncg = ncgg(igrid) - j1;
            if (ncg > 0) {
                for (m = 1; m <= ncg; m++) {
                    nbl = nbl + 1;
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    lsmin = lsminn(nbl);
                    FortranArray3DRef<double> wklsmin(&wk(lsmin), jdim, kdim, idim);
                    FortranArray3DRef<double> wsnk0(&w(lsnk0), jdim-1, kdim-1, idim-1);
                    distcc(jdim, kdim, idim, wklsmin, wsnk0);
                }
            }
        }
    }
} // end findmin_new


void getpts(int& imn, int& jmn, int& kmn, int& imx, int& jmx, int& kmx,
            int& imp1, int& jmp1, int& kmp1,
            FortranArray3DRef<double> x, FortranArray3DRef<double> y,
            FortranArray3DRef<double> z,
            int& nface, int& n1beg, int& n1end, int& n2beg, int& n2end,
            FortranArray2DRef<double> xs, int& nsurf, int& nbl,
            FortranArray3DRef<double> blank,
            FortranArray3DRef<double> bci, FortranArray3DRef<double> bcj,
            FortranArray3DRef<double> bck)
{
    int ipts = n1end - n1beg + 1;
    double slarge = 1.e30;
    double hole   = 1.;
    int mm = 1;
    if (nface == 2 || nface == 4 || nface == 6) mm = 2;

    int iw, iw1, jw, jw1, kw, kw1;
    int jj, j1, kk, k1, ii, i1, i;
    double bciuse1, bciuse2, bciuse;
    double blankuse1, blankuse2, blankuse;
    double bcjuse1, bcjuse2, bcjuse;
    double bckuse1, bckuse2, bckuse;
    double factr;

    switch (nface) {
    case 1:
        iw  = imn;
        iw1 = imn;
        goto label201;
    case 2:
        iw  = imx;
        iw1 = imx - 1;
        label201:
        for (kw = n2beg; kw <= n2end; kw++) {
            kk = std::min(kw, n2end-1);
            k1 = std::max(kw-1, n2beg);
            for (jw = n1beg; jw <= n1end; jw++) {
                jj = std::min(jw, n1end-1);
                j1 = std::max(jw-1, n1beg);
                i = 1 + jw - n1beg + ipts*(kw - n2beg);
                bciuse1 = ccmax(bci(jj,kk,mm), bci(j1,kk,mm));
                bciuse2 = ccmax(bci(jj,k1,mm), bci(j1,k1,mm));
                bciuse  = ccmax(bciuse1, bciuse2);
                blankuse1 = ccmax(blank(jj,kk,iw1), blank(j1,kk,iw1));
                blankuse2 = ccmax(blank(jj,k1,iw1), blank(j1,k1,iw1));
                blankuse  = ccmax(blankuse1, blankuse2);
                factr = slarge * ((1. - blankuse)*hole + (1. - bciuse));
                xs(1,i) = x(jw,kw,iw);
                xs(2,i) = y(jw,kw,iw);
                xs(3,i) = z(jw,kw,iw);
                xs(4,i) = factr;
            }
        }
        return;
    case 3:
        jw  = jmn;
        jw1 = jmn;
        goto label401;
    case 4:
        jw  = jmx;
        jw1 = jmx - 1;
        label401:
        for (iw = n2beg; iw <= n2end; iw++) {
            ii = std::min(iw, n2end-1);
            i1 = std::max(iw-1, n2beg);
            for (kw = n1beg; kw <= n1end; kw++) {
                kk = std::min(kw, n1end-1);
                k1 = std::max(kw-1, n1beg);
                i = 1 + kw - n1beg + ipts*(iw - n2beg);
                bcjuse1 = ccmax(bcj(kk,ii,mm), bcj(k1,ii,mm));
                bcjuse2 = ccmax(bcj(kk,i1,mm), bcj(k1,i1,mm));
                bcjuse  = ccmax(bcjuse1, bcjuse2);
                blankuse1 = ccmax(blank(jw1,kk,ii), blank(jw1,k1,ii));
                blankuse2 = ccmax(blank(jw1,kk,i1), blank(jw1,k1,i1));
                blankuse  = ccmax(blankuse1, blankuse2);
                factr = slarge * ((1. - blankuse)*hole + (1. - bcjuse));
                xs(1,i) = x(jw,kw,iw);
                xs(2,i) = y(jw,kw,iw);
                xs(3,i) = z(jw,kw,iw);
                xs(4,i) = factr;
            }
        }
        return;
    case 5:
        kw  = kmn;
        kw1 = kmn;
        goto label601;
    case 6:
        kw  = kmx;
        kw1 = kmx - 1;
        label601:
        for (jw = n2beg; jw <= n2end; jw++) {
            jj = std::min(jw, n2end-1);
            j1 = std::max(jw-1, n2beg);
            for (iw = n1beg; iw <= n1end; iw++) {
                ii = std::min(iw, n1end-1);
                i1 = std::max(iw-1, n1beg);
                i = 1 + iw - n1beg + ipts*(jw - n2beg);
                bckuse1 = ccmax(bck(jj,ii,mm), bck(j1,ii,mm));
                bckuse2 = ccmax(bck(jj,i1,mm), bck(j1,i1,mm));
                bckuse  = ccmax(bckuse1, bckuse2);
                blankuse1 = ccmax(blank(jj,kw1,ii), blank(j1,kw1,ii));
                blankuse2 = ccmax(blank(jj,kw1,i1), blank(j1,kw1,i1));
                blankuse  = ccmax(blankuse1, blankuse2);
                factr = slarge * ((1. - blankuse)*hole + (1. - bckuse));
                xs(1,i) = x(jw,kw,iw);
                xs(2,i) = y(jw,kw,iw);
                xs(3,i) = z(jw,kw,iw);
                xs(4,i) = factr;
            }
        }
        return;
    default:
        break;
    }
}


void getptsbb(int& imn, int& jmn, int& kmn, int& imx, int& jmx, int& kmx,
              int& imp1, int& jmp1, int& kmp1,
              FortranArray3DRef<double> x, FortranArray3DRef<double> y,
              FortranArray3DRef<double> z,
              int& nface, int& n1beg, int& n1end, int& n2beg, int& n2end,
              FortranArray2DRef<double> xs, FortranArray2DRef<int> ixs,
              int& nsurf, int& nbl,
              FortranArray3DRef<double> blank,
              FortranArray3DRef<double> bci, FortranArray3DRef<double> bcj,
              FortranArray3DRef<double> bck)
{
    int ipts = n1end - n1beg + 1;
    int jpts = n2end - n2beg + 1;
    double slarge = 1.e30;
    double hole   = 1.;
    int mm = 1;
    if (nface == 2 || nface == 4 || nface == 6) mm = 2;

    int iw, iw1, jw, jw1, kw, kw1;
    int jj, j1, kk, k1, ii, i1, i;
    double bciuse1, bciuse2, bciuse;
    double blankuse1, blankuse2, blankuse;
    double bcjuse1, bcjuse2, bcjuse;
    double bckuse1, bckuse2, bckuse;
    double factr;

    switch (nface) {
    case 1:
        iw  = imn;
        iw1 = imn;
        goto label201;
    case 2:
        iw  = imx;
        iw1 = imx - 1;
        label201:
        for (kw = n2beg; kw <= n2end; kw++) {
            kk = std::min(kw, n2end-1);
            k1 = std::max(kw-1, n2beg);
            for (jw = n1beg; jw <= n1end; jw++) {
                jj = std::min(jw, n1end-1);
                j1 = std::max(jw-1, n1beg);
                i = 1 + jw - n1beg + ipts*(kw - n2beg);
                bciuse1 = ccmax(bci(jj,kk,mm), bci(j1,kk,mm));
                bciuse2 = ccmax(bci(jj,k1,mm), bci(j1,k1,mm));
                bciuse  = ccmax(bciuse1, bciuse2);
                blankuse1 = ccmax(blank(jj,kk,iw1), blank(j1,kk,iw1));
                blankuse2 = ccmax(blank(jj,k1,iw1), blank(j1,k1,iw1));
                blankuse  = ccmax(blankuse1, blankuse2);
                factr = slarge * ((1. - blankuse)*hole + (1. - bciuse));
                xs(1,i) = x(jw,kw,iw);
                xs(2,i) = y(jw,kw,iw);
                xs(3,i) = z(jw,kw,iw);
                xs(4,i) = factr;
                ixs(1,i) = iw;
                ixs(2,i) = jw;
                ixs(3,i) = kw;
                ixs(4,i) = nbl;
            }
        }
        return;
    case 3:
        jw  = jmn;
        jw1 = jmn;
        goto label401;
    case 4:
        jw  = jmx;
        jw1 = jmx - 1;
        label401:
        for (iw = n2beg; iw <= n2end; iw++) {
            ii = std::min(iw, n2end-1);
            i1 = std::max(iw-1, n2beg);
            for (kw = n1beg; kw <= n1end; kw++) {
                kk = std::min(kw, n1end-1);
                k1 = std::max(kw-1, n1beg);
                i = 1 + kw - n1beg + ipts*(iw - n2beg);
                bcjuse1 = ccmax(bcj(kk,ii,mm), bcj(k1,ii,mm));
                bcjuse2 = ccmax(bcj(kk,i1,mm), bcj(k1,i1,mm));
                bcjuse  = ccmax(bcjuse1, bcjuse2);
                blankuse1 = ccmax(blank(jw1,kk,ii), blank(jw1,k1,ii));
                blankuse2 = ccmax(blank(jw1,kk,i1), blank(jw1,k1,i1));
                blankuse  = ccmax(blankuse1, blankuse2);
                factr = slarge * ((1. - blankuse)*hole + (1. - bcjuse));
                xs(1,i) = x(jw,kw,iw);
                xs(2,i) = y(jw,kw,iw);
                xs(3,i) = z(jw,kw,iw);
                xs(4,i) = factr;
                ixs(1,i) = iw;
                ixs(2,i) = jw;
                ixs(3,i) = kw;
                ixs(4,i) = nbl;
            }
        }
        return;
    case 5:
        kw  = kmn;
        kw1 = kmn;
        goto label601;
    case 6:
        kw  = kmx;
        kw1 = kmx - 1;
        label601:
        for (jw = n2beg; jw <= n2end; jw++) {
            jj = std::min(jw, n2end-1);
            j1 = std::max(jw-1, n2beg);
            for (iw = n1beg; iw <= n1end; iw++) {
                ii = std::min(iw, n1end-1);
                i1 = std::max(iw-1, n1beg);
                i = 1 + iw - n1beg + ipts*(jw - n2beg);
                bckuse1 = ccmax(bck(jj,ii,mm), bck(j1,ii,mm));
                bckuse2 = ccmax(bck(jj,i1,mm), bck(j1,i1,mm));
                bckuse  = ccmax(bckuse1, bckuse2);
                blankuse1 = ccmax(blank(jj,kw1,ii), blank(j1,kw1,ii));
                blankuse2 = ccmax(blank(jj,kw1,i1), blank(j1,kw1,i1));
                blankuse  = ccmax(blankuse1, blankuse2);
                factr = slarge * ((1. - blankuse)*hole + (1. - bckuse));
                xs(1,i) = x(jw,kw,iw);
                xs(2,i) = y(jw,kw,iw);
                xs(3,i) = z(jw,kw,iw);
                xs(4,i) = factr;
                ixs(1,i) = iw;
                ixs(2,i) = jw;
                ixs(3,i) = kw;
                ixs(4,i) = nbl;
            }
        }
        return;
    default:
        break;
    }
}


void bbdist(int& ng, FortranArray2DRef<double> grid, int& nsurf,
            FortranArray2DRef<double> surf, int& nbb,
            FortranArray2DRef<double> bbdef, FortranArray2DRef<int> ipv,
            FortranArray1DRef<int> vlist, FortranArray1DRef<double> dist,
            FortranArray1DRef<int> idist, int& ncalc,
            FortranArray1DRef<double> wk3d5,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid)
{
    int wrk, tsurf, test;
    wrk   = ifalloc(nou, bou, nbuf, ibufdim, myid, nbb);
    { int tmp = 4*nsurf; tsurf = ifalloc(nou, bou, nbuf, ibufdim, myid, tmp); }
    test  = ifalloc(nou, bou, nbuf, ibufdim, myid, nsurf);
    {
        FortranArray1DRef<double> wkwrk(&wk3d5(wrk), nbb);
        FortranArray2DRef<double> wktsurf(&wk3d5(tsurf), nsurf, 4);
        FortranArray1DRef<double> wktest(&wk3d5(test), nsurf);
        bbdst1(ng, grid, nsurf, surf, nbb, bbdef, ipv, vlist,
               dist, idist, ncalc, wkwrk, wktsurf, wktest);
    }
    ffree(nou, bou, nbuf, ibufdim, myid, nsurf);
    { int tmp = 4*nsurf; ffree(nou, bou, nbuf, ibufdim, myid, tmp); }
    ffree(nou, bou, nbuf, ibufdim, myid, nbb);
}

void bbdst1(int& ng, FortranArray2DRef<double> grid, int& nsurf,
            FortranArray2DRef<double> surf, int& nbb,
            FortranArray2DRef<double> bbdef, FortranArray2DRef<int> ipv,
            FortranArray1DRef<int> vlist, FortranArray1DRef<double> dist,
            FortranArray1DRef<int> idist, int& ncalc,
            FortranArray1DRef<double> wrk,
            FortranArray2DRef<double> tsurf, FortranArray1DRef<double> test)
{
    int i, j, k, ii, jj, kk, n, l;
    double x, y, z, smin, bbmin_val, testmin;
    int jp, kmin;

    // copy surf into tsurf reordered by vlist
    for (i = 1; i <= nsurf; i++) {
        kk = vlist(i);
        tsurf(i,1) = surf(kk,1);
        tsurf(i,2) = surf(kk,2);
        tsurf(i,3) = surf(kk,3);
        tsurf(i,4) = surf(kk,4);
    }
    ncalc = 0;
    for (i = 1; i <= ng; i++) {
        x = grid(i,1);
        y = grid(i,2);
        z = grid(i,3);
        smin = 1.0e34;
        jp = 0;
        // build table of distances from grid to bounding boxes
        for (j = 1; j <= nbb; j++) {
            double px = x;
            if ((float)px <= (float)bbdef(1,j)) px = bbdef(1,j);
            if ((float)px >= (float)bbdef(2,j)) px = bbdef(2,j);
            double py = y;
            if ((float)py <= (float)bbdef(3,j)) py = bbdef(3,j);
            if ((float)py >= (float)bbdef(4,j)) py = bbdef(4,j);
            double pz = z;
            if ((float)pz <= (float)bbdef(5,j)) pz = bbdef(5,j);
            if ((float)pz >= (float)bbdef(6,j)) pz = bbdef(6,j);
            wrk(j) = (x-px)*(x-px) + (y-py)*(y-py) + (z-pz)*(z-pz);
        }
        for (j = 1; j <= nbb; j++) {
            // find nearest bounding box that has not been searched
            bbmin_val = wrk(1);
            for (ii = 2; ii <= nbb; ii++) {
                bbmin_val = ccmin(bbmin_val, wrk(ii));
            }
            jj = isrcheq(nbb, wrk, 1, bbmin_val);
            wrk(jj) = 2.0e34;
            // stop searching when nearest bounding box is too far away
            if ((float)bbmin_val > (float)smin) goto label201;
            n = ipv(1,jj);
            l = ipv(2,jj);
            testmin = 1.0e34;
            for (k = 1; k <= n; k++) {
                double xs = tsurf(l,1);
                double ys = tsurf(l,2);
                double zs = tsurf(l,3);
                test(k) = (x-xs)*(x-xs) + (y-ys)*(y-ys) + (z-zs)*(z-zs) + tsurf(l,4);
                testmin = ccmin(testmin, test(k));
                l = l + 1;
            }
            ncalc = ncalc + n;
            if ((float)testmin < (float)smin) {
                FortranArray1DRef<double> test_sub(&test(1), n);
                kmin = isrcheq(n, test_sub, 1, testmin);
                smin = testmin;
                jp = vlist(ipv(2,jj) + kmin - 1);
            }
        }
        label201:
        dist(i) = std::sqrt(smin);
        idist(i) = jp;
    }
}


void calc_dist(int& imn, int& jmn, int& kmn, int& imx, int& jmx, int& kmx,
               int& imp1, int& jmp1, int& kmp1,
               int& imp2, int& jmp2, int& kmp2,
               FortranArray3DRef<double> x, FortranArray3DRef<double> y,
               FortranArray3DRef<double> z, FortranArray3DRef<double> smin,
               int& nsurf, FortranArray2DRef<double> surf, int& nbb,
               FortranArray2DRef<double> bbdef, FortranArray2DRef<int> ipv,
               FortranArray1DRef<int> vlist, FortranArray1DRef<int> ntri,
               FortranArray2DRef<int> iptri,
               FortranArray1DRef<double> wk3d5, FortranArray1DRef<int> iwrk,
               FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
               int& nbuf, int& ibufdim, int& myid)
{
    int ng = (imx-imn+1)*(jmx-jmn+1)*(kmx-kmn+1);
    int grid; { int tmp = 3*ng; grid = ifalloc(nou, bou, nbuf, ibufdim, myid, tmp); }
    int dist  = ifalloc(nou, bou, nbuf, ibufdim, myid, ng);
    int idist = iialloc(nou, bou, nbuf, ibufdim, myid, ng);

    int ig = 0;
    for (int i = imn; i <= imx; i++) {
        for (int k = kmn; k <= kmx; k++) {
            for (int j = jmn; j <= jmx; j++) {
                ig++;
                wk3d5(grid+ig-1) = x(j,k,i);
                wk3d5(grid+ig+ng-1) = y(j,k,i);
                wk3d5(grid+ig+2*ng-1) = z(j,k,i);
            }
        }
    }
    {
        FortranArray2DRef<double> wkgrid(&wk3d5(grid), ng, 3);
        FortranArray1DRef<double> wkdist(&wk3d5(dist), ng);
        FortranArray1DRef<int>    iwkidist(&iwrk(idist), ng);
        int ncalc;
        bbdist(ng, wkgrid, nsurf, surf, nbb, bbdef, ipv, vlist,
               wkdist, iwkidist, ncalc, wk3d5, nou, bou, nbuf, ibufdim, myid);
    }

    // put in calculation to triangles
    for (int i = 1; i <= ng; i++) {
        int isurf_loc = iwrk(idist+i-1);
        double xp = wk3d5(grid+i-1);
        double yp = wk3d5(grid+i+ng-1);
        double zp = wk3d5(grid+i+2*ng-1);
        int numtri = ntri(isurf_loc);
        if (numtri == 4) {
            for (int itri = 1; itri <= numtri; itri++) {
                int it1 = iptri(isurf_loc, 2*itri-1);
                int it2 = iptri(isurf_loc, 2*itri);
                triang(wk3d5(dist+i-1), xp, yp, zp,
                       surf(isurf_loc,1), surf(isurf_loc,2), surf(isurf_loc,3), surf(isurf_loc,4),
                       surf(it1,1), surf(it1,2), surf(it1,3), surf(it1,4),
                       surf(it2,1), surf(it2,2), surf(it2,3), surf(it2,4));
            }
        } else {
            double xps = surf(isurf_loc,1);
            double yps = surf(isurf_loc,2);
            double zps = surf(isurf_loc,3);
            int ixmin;
            for (ixmin = isurf_loc; ixmin >= 1; ixmin--) {
                if ((float)surf(ixmin,1) < (float)xps) goto label10;
            }
            ixmin = 0;
            label10:
            ixmin = ixmin + 1;
            int ixmax;
            for (ixmax = isurf_loc; ixmax <= nsurf; ixmax++) {
                if ((float)surf(ixmax,1) > (float)xps) goto label20;
            }
            ixmax = nsurf + 1;
            label20:
            ixmax = ixmax - 1;
            for (int isurf2 = ixmin; isurf2 <= ixmax; isurf2++) {
                if ((float)surf(isurf2,2) == (float)yps &&
                    (float)surf(isurf2,3) == (float)zps &&
                    ntri(isurf2) != 4) {
                    numtri = ntri(isurf2);
                    for (int itri = 1; itri <= numtri; itri++) {
                        int it1 = iptri(isurf2, 2*itri-1);
                        int it2 = iptri(isurf2, 2*itri);
                        triang(wk3d5(dist+i-1), xp, yp, zp,
                               surf(isurf2,1), surf(isurf2,2), surf(isurf2,3), surf(isurf2,4),
                               surf(it1,1), surf(it1,2), surf(it1,3), surf(it1,4),
                               surf(it2,1), surf(it2,2), surf(it2,3), surf(it2,4));
                    }
                }
            }
        }
    }
    ig = 0;
    for (int i = imn; i <= imx; i++) {
        for (int k = kmn; k <= kmx; k++) {
            for (int j = jmn; j <= jmx; j++) {
                ig++;
                smin(j,k,i) = wk3d5(dist+ig-1);
            }
        }
    }
    ifree(nou, bou, nbuf, ibufdim, myid, ng);
    ffree(nou, bou, nbuf, ibufdim, myid, ng);
    { int tmp = 3*ng; ffree(nou, bou, nbuf, ibufdim, myid, tmp); }
}


void collect_surf(int& imn, int& jmn, int& kmn,
                  int& imxs, int& jmxs, int& kmxs,
                  int& imp1s, int& jmp1s, int& kmp1s,
                  FortranArray2DRef<double> xs,
                  int& nface, int& n1beg, int& n1end, int& n2beg, int& n2end,
                  int& nsurf, FortranArray2DRef<double> surf, int& isurf,
                  FortranArray1DRef<int> ntri, FortranArray2DRef<int> iptri)
{
    int isurf_beg = isurf;
    int ipts = n1end - n1beg + 1;

    // Collect surface points
    for (int k = n2beg; k <= n2end; k++) {
        for (int j = n1beg; j <= n1end; j++) {
            isurf++;
            int i = 1 + (j - n1beg) + ipts*(k - n2beg);
            surf(isurf,1) = xs(1,i);
            surf(isurf,2) = xs(2,i);
            surf(isurf,3) = xs(3,i);
            surf(isurf,4) = xs(4,i);
        }
    }

    // Now add collection of local triangles
    isurf = isurf_beg;
    int n1inc = 1;
    int n2inc = n1end - n1beg + 1;
    for (int j = n2beg; j <= n2end; j++) {
        for (int i = n1beg; i <= n1end; i++) {
            isurf++;
            int itri = 0;
            for (int jj = std::max(j-1, n2beg); jj <= std::min(j+1, n2end); jj++) {
                for (int ii = std::max(i-1, n1beg); ii <= std::min(i+1, n1end); ii++) {
                    if (ii != i && jj != j) {
                        itri++;
                        iptri(isurf, 2*itri-1) = isurf + (ii-i)*n1inc;
                        iptri(isurf, 2*itri)   = isurf + (jj-j)*n2inc;
                    }
                }
            }
            ntri(isurf) = itri;
        }
    }
}

void distcc(int& jdim, int& kdim, int& idim,
            FortranArray3DRef<double> smingp, FortranArray3DRef<double> smincc)
{
    int iw = 0;
    for (int i = 1; i <= idim-1; i++) {
        iw++;
        int kw = 0;
        for (int k = 1; k <= kdim-1; k++) {
            kw++;
            int jw = 0;
            for (int j = 1; j <= jdim-1; j++) {
                jw++;
                smincc(jw,kw,iw) = 0.125*(smingp(j,k,i)     + smingp(j,k,i+1)
                                         + smingp(j+1,k,i)   + smingp(j+1,k,i+1)
                                         + smingp(j,k+1,i)   + smingp(j,k+1,i+1)
                                         + smingp(j+1,k+1,i) + smingp(j+1,k+1,i+1));
            }
        }
    }
}

void distcg(int& jdimf, int& kdimf, int& idimf,
            FortranArray3DRef<double> sminf,
            int& jdim, int& kdim, int& idim,
            FortranArray3DRef<double> sminc)
{
    int iinc = 2;
    if (idimf == 2) iinc = 1;
    int ii = 0;
    for (int i = 1; i <= idimf; i += iinc) {
        ii++;
        int kk = 0;
        for (int k = 1; k <= kdimf; k += 2) {
            kk++;
            int jj = 0;
            for (int j = 1; j <= jdimf; j += 2) {
                jj++;
                sminc(jj,kk,ii) = sminf(j,k,i);
            }
        }
    }
}

void initf(int& isize)
{
    cmn_alloc.ifptr = 1;
    cmn_alloc.ifmax = isize;
}

void initi(int& isize)
{
    cmn_alloc.iptr = 1;
    cmn_alloc.imax = isize;
}

void ifree(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
           int& nbuf, int& ibufdim, int& myid, int& isize)
{
    int& iptr = cmn_alloc.iptr;
    if (isize >= iptr) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " stopping...ifree failed: %d %d", isize, iptr);
        int ierrflg = -1;
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    iptr = iptr - isize;
}

void ffree(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
           int& nbuf, int& ibufdim, int& myid, int& isize)
{
    int& ifptr = cmn_alloc.ifptr;
    if (isize >= ifptr) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " stopping...ffree failed: %d %d", isize, ifptr);
        int ierrflg = -1;
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    ifptr = ifptr - isize;
}



void makebb(int& minpb0, int& maxlv0, int& nwv0,
            FortranArray1DRef<double> wx, FortranArray1DRef<double> wy,
            FortranArray1DRef<double> wz,
            int& nbb0, int& bblen,
            FortranArray2DRef<double> bbdef, FortranArray2DRef<int> iv,
            FortranArray1DRef<int> vlist, FortranArray1DRef<int> iwrk,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid)
{
    int& nbbv   = cmn_bbcom.nbbv;
    int& minpbb = cmn_bbcom.minpbb;
    int& maxlev = cmn_bbcom.maxlev;
    int& nwv    = cmn_bbcom.nwv;
    int& nbb    = cmn_bbcom.nbb;
    int& totbbv = cmn_bbcom.totbbv;
    // bbmin(3), bbmax(3) are float arrays in cmn_bbcom

    minpbb = minpb0;
    maxlev = maxlv0;
    nwv    = nwv0;
    nbb    = 0;
    totbbv = 0;

    // The array bbi is a list of the surface points in wx,wy,wz.
    int bbi = iialloc(nou, bou, nbuf, ibufdim, myid, nwv);
    for (int i = 1; i <= nwv; i++) {
        iwrk(bbi+i-1) = i;
    }

    // Create a bounding box containing all the surface points.
    // calcbb expects FortranArray1DRef<double> bbmin/bbmax
    // but cmn_bbcom has float arrays - use local double arrays
    FortranArray1D<double> bbmin_d(3), bbmax_d(3);
    FortranArray1DRef<int> ibbi(&iwrk(bbi), nwv);
    calcbb(bbmin_d.ref(), bbmax_d.ref(), nwv, ibbi, wx, wy, wz);
    // copy back to cmn_bbcom float arrays
    for (int i = 1; i <= 3; i++) {
        cmn_bbcom.bbmin[i-1] = (float)bbmin_d(i);
        cmn_bbcom.bbmax[i-1] = (float)bbmax_d(i);
    }

    // Now subdivide the bounding box.
    int level = 0;
    spltbb(level, bbdef, iv, vlist, ibbi, wx, wy, wz);
    ifree(nou, bou, nbuf, ibufdim, myid, nwv);
    nbb0  = nbb;
    bblen = totbbv;
}

void calcbb(FortranArray1DRef<double> bbmin, FortranArray1DRef<double> bbmax,
            int& nvi, FortranArray1DRef<int> vi,
            FortranArray1DRef<double> wx, FortranArray1DRef<double> wy,
            FortranArray1DRef<double> wz)
{
    const double bval =  1.0e20;
    const double sval = -1.0e20;

    for (int i = 1; i <= 3; i++) {
        bbmin(i) = bval;
        bbmax(i) = sval;
    }
    for (int j = 1; j <= nvi; j++) {
        int i = vi(j);
        if ((float)wx(i) < (float)bbmin(1)) bbmin(1) = wx(i);
        if ((float)wy(i) < (float)bbmin(2)) bbmin(2) = wy(i);
        if ((float)wz(i) < (float)bbmin(3)) bbmin(3) = wz(i);
        if ((float)wx(i) > (float)bbmax(1)) bbmax(1) = wx(i);
        if ((float)wy(i) > (float)bbmax(2)) bbmax(2) = wy(i);
        if ((float)wz(i) > (float)bbmax(3)) bbmax(3) = wz(i);
    }
}

void getvrt(FortranArray1DRef<double> bbmin, FortranArray1DRef<double> bbmax,
            int& nwv, int& nbbv, FortranArray1DRef<int> bbi,
            FortranArray1DRef<double> wx, FortranArray1DRef<double> wy,
            FortranArray1DRef<double> wz)
{
    nbbv = 0;
    for (int i = 1; i <= nwv; i++) {
        if ((float)wx(i) >= (float)bbmin(1) && (float)wx(i) <= (float)bbmax(1) &&
            (float)wy(i) >= (float)bbmin(2) && (float)wy(i) <= (float)bbmax(2) &&
            (float)wz(i) >= (float)bbmin(3) && (float)wz(i) <= (float)bbmax(3)) {
            nbbv++;
            bbi(nbbv) = i;
        }
    }
}

void shells(int& n, FortranArray1DRef<int> is, FortranArray1DRef<double> s)
{
    int gap = n / 2;
    while ((float)gap > 0) {
        for (int i = gap+1; i <= n; i++) {
            int j = i - gap;
            while (true) {
                if (j <= 0) break;
                if ((float)s(is(j)) <= (float)s(is(j+gap))) break;
                int temp = is(j);
                is(j) = is(j+gap);
                is(j+gap) = temp;
                j = j - gap;
            }
        }
        gap = gap / 2;
    }
}


void spltbb(int& level0, FortranArray2DRef<double> bbdef,
            FortranArray2DRef<int> iv, FortranArray1DRef<int> vlist,
            FortranArray1DRef<int> bbi,
            FortranArray1DRef<double> wx, FortranArray1DRef<double> wy,
            FortranArray1DRef<double> wz)
{
    const int MAXBB = 200;
    int& nbbv   = cmn_bbcom.nbbv;
    int& minpbb = cmn_bbcom.minpbb;
    int& maxlev = cmn_bbcom.maxlev;
    int& nwv    = cmn_bbcom.nwv;
    int& nbb    = cmn_bbcom.nbb;
    int& totbbv = cmn_bbcom.totbbv;

    // stack(7,MAXBB) - local 2D array
    FortranArray2D<double> stack(7, MAXBB);
    FortranArray1D<double> parmin(3), parmax(3);
    FortranArray1D<double> plmin(3), plmax(3);
    FortranArray1D<double> prmin(3), prmax(3);

    // bbmin/bbmax from cmn_bbcom (float) - use local double copies
    FortranArray1D<double> bbmin_d(3), bbmax_d(3);
    for (int i = 1; i <= 3; i++) {
        bbmin_d(i) = cmn_bbcom.bbmin[i-1];
        bbmax_d(i) = cmn_bbcom.bbmax[i-1];
    }

    int isptr = 1;
    push(stack.ref(), isptr, bbmin_d.ref(), bbmax_d.ref(), level0);
    int k = 1;

    while (isptr > 1) {
        int level;
        pop(stack.ref(), isptr, parmin.ref(), parmax.ref(), level);
        getvrt(parmin.ref(), parmax.ref(), nwv, nbbv, bbi, wx, wy, wz);
        if (nbbv >= minpbb && level < maxlev) {
            double bx = parmax(1) - parmin(1);
            double by = parmax(2) - parmin(2);
            double bz = parmax(3) - parmin(3);
            int iright;
            if ((float)bx > (float)by && (float)bx > (float)bz) {
                shells(nbbv, bbi, wx);
                iright = nbbv / 2;
                int i1 = bbi(iright);
                int i2 = bbi(iright+1);
                while (true) {
                    if ((float)wx(i1) != (float)wx(i2) || iright >= nbbv) break;
                    iright++;
                    if (iright >= nbbv) break;
                    i1 = bbi(iright);
                    i2 = bbi(iright+1);
                }
            } else if ((float)by > (float)bx && (float)by > (float)bz) {
                shells(nbbv, bbi, wy);
                iright = nbbv / 2;
                int i1 = bbi(iright);
                int i2 = bbi(iright+1);
                while (true) {
                    if ((float)wy(i1) != (float)wy(i2) || iright >= nbbv) break;
                    iright++;
                    if (iright >= nbbv) break;
                    i1 = bbi(iright);
                    i2 = bbi(iright+1);
                }
            } else {
                shells(nbbv, bbi, wz);
                iright = nbbv / 2;
                int i1 = bbi(iright);
                int i2 = bbi(iright+1);
                while (true) {
                    if ((float)wz(i1) != (float)wz(i2) || iright >= nbbv) break;
                    iright++;
                    if (iright >= nbbv) break;
                    i1 = bbi(iright);
                    i2 = bbi(iright+1);
                }
            }
            calcbb(plmin.ref(), plmax.ref(), iright, bbi, wx, wy, wz);
            int nright = nbbv - iright;
            FortranArray1DRef<int> bbi_right(&bbi(iright+1), nright);
            calcbb(prmin.ref(), prmax.ref(), nright, bbi_right, wx, wy, wz);
            int lp1 = level + 1;
            push(stack.ref(), isptr, prmin.ref(), prmax.ref(), lp1);
            push(stack.ref(), isptr, plmin.ref(), plmax.ref(), lp1);
        } else {
            if (nbbv > 0) {
                nbb++;
                bbdef(1,nbb) = parmin(1);
                bbdef(2,nbb) = parmax(1);
                bbdef(3,nbb) = parmin(2);
                bbdef(4,nbb) = parmax(2);
                bbdef(5,nbb) = parmin(3);
                bbdef(6,nbb) = parmax(3);
                iv(1,nbb) = nbbv;
                iv(2,nbb) = k;
                for (int i = 1; i <= nbbv; i++) {
                    vlist(k) = bbi(i);
                    k++;
                    totbbv++;
                }
            }
        }
    }
}

void push(FortranArray2DRef<double> st, int& is,
          FortranArray1DRef<double> bbmin, FortranArray1DRef<double> bbmax,
          int& lev)
{
    for (int i = 1; i <= 3; i++) {
        int ii = i + 3;
        st(i,is)  = bbmin(i);
        st(ii,is) = bbmax(i);
    }
    st(7,is) = lev + 0.1;
    is = is + 1;
}

void pop(FortranArray2DRef<double> st, int& is,
         FortranArray1DRef<double> bbmin, FortranArray1DRef<double> bbmax,
         int& lev)
{
    is = is - 1;
    for (int i = 1; i <= 3; i++) {
        int ii = i + 3;
        bbmin(i) = st(i,is);
        bbmax(i) = st(ii,is);
    }
    lev = (int)st(7,is);
}


void sort_x(int& nsurf, FortranArray2DRef<double> surf,
            FortranArray1DRef<int> ntri, FortranArray2DRef<int> iptr,
            FortranArray1DRef<double> wk3d5, FortranArray1DRef<int> iwrk,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid)
{
    int iperm = iialloc(nou, bou, nbuf, ibufdim, myid, nsurf);
    int temp  = ifalloc(nou, bou, nbuf, ibufdim, myid, nsurf);
    {
        FortranArray1DRef<int>    iwkiperm(&iwrk(iperm), nsurf);
        FortranArray1DRef<double> wktemp(&wk3d5(temp), nsurf);
        // heap_sort on surf column 1 (x-coordinates)
        FortranArray1DRef<double> surf1(&surf(1,1), nsurf);
        heap_sort(nsurf, surf1, iwkiperm);
        for (int j = 1; j <= 4; j++) {
            FortranArray1DRef<double> surfj(&surf(1,j), nsurf);
            move_real(nsurf, surfj, iwkiperm, wktemp);
        }
    }
    ffree(nou, bou, nbuf, ibufdim, myid, nsurf);
    int itemp = iialloc(nou, bou, nbuf, ibufdim, myid, nsurf);
    {
        FortranArray1DRef<int> iwkiperm(&iwrk(iperm), nsurf);
        FortranArray1DRef<int> iwkitemp(&iwrk(itemp), nsurf);
        move_integer(nsurf, ntri, iwkiperm, iwkitemp);
        for (int j = 1; j <= 8; j++) {
            FortranArray1DRef<int> iptrj(&iptr(1,j), nsurf);
            move_integer(nsurf, iptrj, iwkiperm, iwkitemp);
        }
    }
    ifree(nou, bou, nbuf, ibufdim, myid, nsurf);
    int inv_iperm = iialloc(nou, bou, nbuf, ibufdim, myid, nsurf);
    {
        FortranArray1DRef<int> iwkiperm(&iwrk(iperm), nsurf);
        FortranArray1DRef<int> iwkinv(&iwrk(inv_iperm), nsurf);
        for (int i = 1; i <= nsurf; i++) {
            int ii = iwkiperm(i);
            iwkinv(ii) = i;
        }
        for (int i = 1; i <= nsurf; i++) {
            int numtri = ntri(i);
            for (int j = 1; j <= 2*numtri; j++) {
                int ii = iptr(i,j);
                iptr(i,j) = iwkinv(ii);
            }
        }
    }
    ifree(nou, bou, nbuf, ibufdim, myid, nsurf);
    ifree(nou, bou, nbuf, ibufdim, myid, nsurf);
}

void move_real(int& n, FortranArray1DRef<double> x,
               FortranArray1DRef<int> iperm, FortranArray1DRef<double> temp)
{
    for (int i = 1; i <= n; i++) {
        temp(i) = x(iperm(i));
    }
    for (int i = 1; i <= n; i++) {
        x(i) = temp(i);
    }
}

void move_integer(int& n, FortranArray1DRef<int> ix,
                  FortranArray1DRef<int> iperm, FortranArray1DRef<int> itemp)
{
    for (int i = 1; i <= n; i++) {
        itemp(i) = ix(iperm(i));
    }
    for (int i = 1; i <= n; i++) {
        ix(i) = itemp(i);
    }
}

void triang(double& dist, double& pp1, double& pp2, double& pp3,
            double& aa1, double& aa2, double& aa3, double& aa4,
            double& bb1, double& bb2, double& bb3, double& bb4,
            double& cc1, double& cc2, double& cc3, double& cc4)
{
    if ((float)aa4 != 0. || (float)bb4 != 0. || (float)cc4 != 0.) return;

    double a1 = bb1 - aa1;
    double a2 = bb2 - aa2;
    double a3 = bb3 - aa3;
    double b1 = cc1 - aa1;
    double b2 = cc2 - aa2;
    double b3 = cc3 - aa3;
    double p1 = pp1 - aa1;
    double p2 = pp2 - aa2;
    double p3 = pp3 - aa3;
    double daa = a1*a1 + a2*a2 + a3*a3;
    double dab = a1*b1 + a2*b2 + a3*b3;
    double dbb = b1*b1 + b2*b2 + b3*b3;
    double den = dab*dab - daa*dbb;
    double dap, dbp, s, t, r1, r2, r3, dsq, tmp;
    if ((float)den == 0) goto label100;
    dap = a1*p1 + a2*p2 + a3*p3;
    dbp = b1*p1 + b2*p2 + b3*p3;
    s = (dab*dbp - dbb*dap) / den;
    t = (dab*dap - daa*dbp) / den;
    if ((float)s < 0 || (float)t < 0 || (float)(t+s) > 1) goto label100;
    r1 = p1 - s*a1 - t*b1;
    r2 = p2 - s*a2 - t*b2;
    r3 = p3 - s*a3 - t*b3;
    tmp = std::sqrt(r1*r1 + r2*r2 + r3*r3);
    dist = ccmin(dist, tmp);
    return;
    label100:
    dsq = dist * dist;
    if ((float)daa == 0) goto label200;
    dap = a1*p1 + a2*p2 + a3*p3;
    t = dap / daa;
    if ((float)t < 0 || (float)t > 1) goto label200;
    r1 = p1 - t*a1;
    r2 = p2 - t*a2;
    r3 = p3 - t*a3;
    tmp = r1*r1 + r2*r2 + r3*r3;
    dsq = ccmin(dsq, tmp);
    label200:
    if ((float)dbb == 0) goto label300;
    dbp = b1*p1 + b2*p2 + b3*p3;
    t = dbp / dbb;
    if ((float)t < 0. || (float)t > 1) goto label300;
    r1 = p1 - t*b1;
    r2 = p2 - t*b2;
    r3 = p3 - t*b3;
    tmp = r1*r1 + r2*r2 + r3*r3;
    dsq = ccmin(dsq, tmp);
    label300:
    p1 = pp1 - bb1;
    p2 = pp2 - bb2;
    p3 = pp3 - bb3;
    a1 = cc1 - bb1;
    a2 = cc2 - bb2;
    a3 = cc3 - bb3;
    daa = a1*a1 + a2*a2 + a3*a3;
    if ((float)daa == 0) goto label400;
    dap = a1*p1 + a2*p2 + a3*p3;
    t = dap / daa;
    if ((float)t < 0 || (float)t > 1) goto label400;
    r1 = p1 - t*a1;
    r2 = p2 - t*a2;
    r3 = p3 - t*a3;
    tmp = r1*r1 + r2*r2 + r3*r3;
    dsq = ccmin(dsq, tmp);
    label400:
    if ((float)dsq < (float)dist * (float)dist) dist = std::sqrt(dsq);
}


void heap_sort(int& n, FortranArray1DRef<double> arrin, FortranArray1DRef<int> indx)
{
    for (int j = 1; j <= n; j++) {
        indx(j) = j;
    }
    if (n == 1) return;

    int l = n/2 + 1;
    int ir = n;
    int indxt;
    double q;

    while (true) {
        if (l > 1) {
            l--;
            indxt = indx(l);
            q = arrin(indxt);
        } else {
            indxt = indx(ir);
            q = arrin(indxt);
            indx(ir) = indx(1);
            ir--;
            if (ir == 1) {
                indx(1) = indxt;
                return;
            }
        }
        int i = l;
        int j = l + l;
        while (j <= ir) {
            if (j < ir) {
                if ((float)arrin(indx(j)) < (float)arrin(indx(j+1))) j++;
            }
            if ((float)q < (float)arrin(indx(j))) {
                indx(i) = indx(j);
                i = j;
                j = j + j;
            } else {
                j = ir + 1;
            }
        }
        indx(i) = indxt;
    }
}

void bbarthcg(int& jdim, int& kdim, int& idim,
              int& jdimc, int& kdimc, int& idimc,
              FortranArray3DRef<double> xjb, FortranArray3DRef<double> xkb,
              FortranArray4DRef<double> xib, FortranArray3DRef<double> blnum,
              FortranArray3DRef<double> xjbc, FortranArray3DRef<double> xkbc,
              FortranArray4DRef<double> xibc, FortranArray3DRef<double> blnumc)
{
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int iinc = 2;
    if (idim == 2) iinc = 1;
    int ii = 0;
    for (int i = 1; i <= idim1; i += iinc) {
        ii++;
        int kk = 0;
        for (int k = 1; k <= kdim1; k += 2) {
            kk++;
            int jj = 0;
            for (int j = 1; j <= jdim1; j += 2) {
                jj++;
                blnumc(jj,kk,ii) = blnum(j,k,i) + 1;
                xjbc(jj,kk,ii)   = (int)(xjb(j,k,i)/2 + 1);
                xkbc(jj,kk,ii)   = (int)(xkb(j,k,i)/2 + 1);
                xibc(jj,kk,ii,2) = (int)(xib(j,k,i,2)/2 + 1);
            }
        }
    }
}


void calc_distbb(int& imn, int& jmn, int& kmn, int& imx, int& jmx, int& kmx,
                 int& imp1, int& jmp1, int& kmp1,
                 int& imp2, int& jmp2, int& kmp2,
                 FortranArray3DRef<double> x, FortranArray3DRef<double> y,
                 FortranArray3DRef<double> z, FortranArray3DRef<double> smin,
                 int& nsurf, FortranArray2DRef<double> surf, int& nbb,
                 FortranArray2DRef<double> bbdef, FortranArray2DRef<int> ipv,
                 FortranArray1DRef<int> vlist, FortranArray1DRef<int> ntri,
                 FortranArray2DRef<int> iptri, FortranArray2DRef<int> ibbarth,
                 FortranArray3DRef<double> xjb, FortranArray3DRef<double> xkb,
                 FortranArray4DRef<double> xib, FortranArray3DRef<double> blnum,
                 FortranArray1DRef<double> wk3d5, FortranArray1DRef<int> iwrk,
                 int& j1,
                 FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
                 int& nbuf, int& ibufdim, int& myid)
{
    int ng = (imx-imn+1)*(jmx-jmn+1)*(kmx-kmn+1);
    int grid; { int tmp = 3*ng; grid = ifalloc(nou, bou, nbuf, ibufdim, myid, tmp); }
    int dist  = ifalloc(nou, bou, nbuf, ibufdim, myid, ng);
    int idist = iialloc(nou, bou, nbuf, ibufdim, myid, ng);

    int ig = 0;
    for (int i = imn; i <= imx; i++) {
        for (int k = kmn; k <= kmx; k++) {
            for (int j = jmn; j <= jmx; j++) {
                ig++;
                wk3d5(grid+ig-1)      = x(j,k,i);
                wk3d5(grid+ig+ng-1)   = y(j,k,i);
                wk3d5(grid+ig+2*ng-1) = z(j,k,i);
            }
        }
    }
    {
        FortranArray2DRef<double> wkgrid(&wk3d5(grid), ng, 3);
        FortranArray1DRef<double> wkdist(&wk3d5(dist), ng);
        FortranArray1DRef<int>    iwkidist(&iwrk(idist), ng);
        int ncalc;
        bbdist(ng, wkgrid, nsurf, surf, nbb, bbdef, ipv, vlist,
               wkdist, iwkidist, ncalc, wk3d5, nou, bou, nbuf, ibufdim, myid);
    }

    // set xib, xjb, xkb, blnum for grid points
    int ifactr = 1;
    for (int jj = 0; jj < j1; jj++) ifactr *= 2;
    int imx1 = imx - 1;
    int jmx1 = jmx - 1;
    int kmx1 = kmx - 1;
    ig = 0;
    for (int i = imn; i <= imx1; i++) {
        for (int k = kmn; k <= kmx; k++) {
            for (int j = jmn; j <= jmx; j++) {
                ig++;
                if (j == jmx || k == kmx) continue;
                int isurf_loc = iwrk(idist+ig-1);
                xib(j,k,i,2) = (int)((ibbarth(isurf_loc,1)-1)/ifactr + 1);
                xjb(j,k,i)   = (int)((ibbarth(isurf_loc,2)-1)/ifactr + 1);
                xkb(j,k,i)   = (int)((ibbarth(isurf_loc,3)-1)/ifactr + 1);
                blnum(j,k,i) = ibbarth(isurf_loc,4) + j1;
            }
        }
    }

    // put in calculation to triangles
    for (int i = 1; i <= ng; i++) {
        int isurf_loc = iwrk(idist+i-1);
        double xp = wk3d5(grid+i-1);
        double yp = wk3d5(grid+i+ng-1);
        double zp = wk3d5(grid+i+2*ng-1);
        int numtri = ntri(isurf_loc);
        if (numtri == 4) {
            for (int itri = 1; itri <= numtri; itri++) {
                int it1 = iptri(isurf_loc, 2*itri-1);
                int it2 = iptri(isurf_loc, 2*itri);
                triang(wk3d5(dist+i-1), xp, yp, zp,
                       surf(isurf_loc,1), surf(isurf_loc,2), surf(isurf_loc,3), surf(isurf_loc,4),
                       surf(it1,1), surf(it1,2), surf(it1,3), surf(it1,4),
                       surf(it2,1), surf(it2,2), surf(it2,3), surf(it2,4));
            }
        } else {
            double xps = surf(isurf_loc,1);
            double yps = surf(isurf_loc,2);
            double zps = surf(isurf_loc,3);
            int ixmin;
            for (ixmin = isurf_loc; ixmin >= 1; ixmin--) {
                if ((float)surf(ixmin,1) < (float)xps) goto label10;
            }
            ixmin = 0;
            label10:
            ixmin = ixmin + 1;
            int ixmax;
            for (ixmax = isurf_loc; ixmax <= nsurf; ixmax++) {
                if ((float)surf(ixmax,1) > (float)xps) goto label20;
            }
            ixmax = nsurf + 1;
            label20:
            ixmax = ixmax - 1;
            for (int isurf2 = ixmin; isurf2 <= ixmax; isurf2++) {
                if ((float)surf(isurf2,2) == (float)yps &&
                    (float)surf(isurf2,3) == (float)zps &&
                    ntri(isurf2) != 4) {
                    numtri = ntri(isurf2);
                    for (int itri = 1; itri <= numtri; itri++) {
                        int it1 = iptri(isurf2, 2*itri-1);
                        int it2 = iptri(isurf2, 2*itri);
                        triang(wk3d5(dist+i-1), xp, yp, zp,
                               surf(isurf2,1), surf(isurf2,2), surf(isurf2,3), surf(isurf2,4),
                               surf(it1,1), surf(it1,2), surf(it1,3), surf(it1,4),
                               surf(it2,1), surf(it2,2), surf(it2,3), surf(it2,4));
                    }
                }
            }
        }
    }
    ig = 0;
    for (int i = imn; i <= imx; i++) {
        for (int k = kmn; k <= kmx; k++) {
            for (int j = jmn; j <= jmx; j++) {
                ig++;
                smin(j,k,i) = wk3d5(dist+ig-1);
            }
        }
    }
    ifree(nou, bou, nbuf, ibufdim, myid, ng);
    ffree(nou, bou, nbuf, ibufdim, myid, ng);
    int isize_3ng = 3*ng;
    ffree(nou, bou, nbuf, ibufdim, myid, isize_3ng);
}


void collect_surfbb(int& imn, int& jmn, int& kmn,
                    int& imxs, int& jmxs, int& kmxs,
                    int& imp1s, int& jmp1s, int& kmp1s,
                    FortranArray2DRef<double> xs, FortranArray2DRef<int> ixs,
                    int& nface, int& n1beg, int& n1end, int& n2beg, int& n2end,
                    int& nsurf, FortranArray2DRef<double> surf, int& isurf,
                    FortranArray1DRef<int> ntri, FortranArray2DRef<int> iptri,
                    FortranArray2DRef<int> ibbarth, int& nbl)
{
    int isurf_beg = isurf;
    int ipts = n1end - n1beg + 1;

    // Collect surface points
    for (int k = n2beg; k <= n2end; k++) {
        for (int j = n1beg; j <= n1end; j++) {
            isurf++;
            int i = 1 + (j - n1beg) + ipts*(k - n2beg);
            surf(isurf,1)    = xs(1,i);
            surf(isurf,2)    = xs(2,i);
            surf(isurf,3)    = xs(3,i);
            surf(isurf,4)    = xs(4,i);
            ibbarth(isurf,1) = ixs(1,i);
            ibbarth(isurf,2) = ixs(2,i);
            ibbarth(isurf,3) = ixs(3,i);
            ibbarth(isurf,4) = ixs(4,i);
        }
    }

    // Now add collection of local triangles
    isurf = isurf_beg;
    int n1inc = 1;
    int n2inc = n1end - n1beg + 1;
    for (int j = n2beg; j <= n2end; j++) {
        for (int i = n1beg; i <= n1end; i++) {
            isurf++;
            int itri = 0;
            for (int jj = std::max(j-1, n2beg); jj <= std::min(j+1, n2end); jj++) {
                for (int ii = std::max(i-1, n1beg); ii <= std::min(i+1, n1end); ii++) {
                    if (ii != i && jj != j) {
                        itri++;
                        iptri(isurf, 2*itri-1) = isurf + (ii-i)*n1inc;
                        iptri(isurf, 2*itri)   = isurf + (jj-j)*n2inc;
                    }
                }
            }
            ntri(isurf) = itri;
        }
    }
}

void sort_xbb(int& nsurf, FortranArray2DRef<double> surf,
              FortranArray1DRef<int> ntri, FortranArray2DRef<int> iptr,
              FortranArray2DRef<int> ibbarth,
              FortranArray1DRef<double> wk3d5, FortranArray1DRef<int> iwrk,
              FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
              int& nbuf, int& ibufdim, int& myid)
{
    int iperm = iialloc(nou, bou, nbuf, ibufdim, myid, nsurf);
    int temp  = ifalloc(nou, bou, nbuf, ibufdim, myid, nsurf);
    {
        FortranArray1DRef<int>    iwkiperm(&iwrk(iperm), nsurf);
        FortranArray1DRef<double> wktemp(&wk3d5(temp), nsurf);
        FortranArray1DRef<double> surf1(&surf(1,1), nsurf);
        heap_sort(nsurf, surf1, iwkiperm);
        for (int j = 1; j <= 4; j++) {
            FortranArray1DRef<double> surfj(&surf(1,j), nsurf);
            move_real(nsurf, surfj, iwkiperm, wktemp);
        }
    }
    ffree(nou, bou, nbuf, ibufdim, myid, nsurf);
    int itemp = iialloc(nou, bou, nbuf, ibufdim, myid, nsurf);
    {
        FortranArray1DRef<int> iwkiperm(&iwrk(iperm), nsurf);
        FortranArray1DRef<int> iwkitemp(&iwrk(itemp), nsurf);
        move_integer(nsurf, ntri, iwkiperm, iwkitemp);
        for (int j = 1; j <= 8; j++) {
            FortranArray1DRef<int> iptrj(&iptr(1,j), nsurf);
            move_integer(nsurf, iptrj, iwkiperm, iwkitemp);
        }
        for (int j = 1; j <= 4; j++) {
            FortranArray1DRef<int> ibbarthj(&ibbarth(1,j), nsurf);
            move_integer(nsurf, ibbarthj, iwkiperm, iwkitemp);
        }
    }
    ifree(nou, bou, nbuf, ibufdim, myid, nsurf);
    int inv_iperm = iialloc(nou, bou, nbuf, ibufdim, myid, nsurf);
    {
        FortranArray1DRef<int> iwkiperm(&iwrk(iperm), nsurf);
        FortranArray1DRef<int> iwkinv(&iwrk(inv_iperm), nsurf);
        for (int i = 1; i <= nsurf; i++) {
            int ii = iwkiperm(i);
            iwkinv(ii) = i;
        }
        for (int i = 1; i <= nsurf; i++) {
            int numtri = ntri(i);
            for (int j = 1; j <= 2*numtri; j++) {
                int ii = iptr(i,j);
                iptr(i,j) = iwkinv(ii);
            }
        }
    }
    ifree(nou, bou, nbuf, ibufdim, myid, nsurf);
    ifree(nou, bou, nbuf, ibufdim, myid, nsurf);
}

int iialloc(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid, int& isize)
{
    int& iptr = cmn_alloc.iptr;
    int& imax = cmn_alloc.imax;
    int result = iptr;
    iptr = iptr + isize;
    if (iptr > (imax + 1)) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " stopping...iialloc failed: %d %d",
                      isize, iptr - isize);
        int ierrflg = -1;
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    return result;
}

int ifalloc(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid, int& isize)
{
    int& ifptr = cmn_alloc.ifptr;
    int& ifmax = cmn_alloc.ifmax;
    int result = ifptr;
    ifptr = ifptr + isize;
    if (ifptr > (ifmax + 1)) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " stopping...ifalloc failed: %d %d",
                      isize, ifptr - isize);
        int ierrflg = -1;
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    return result;
}


void init(int& nbl, int& jdim, int& kdim, int& idim,
          FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
          FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
          FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0,
          FortranArray4DRef<double> ti0,
          FortranArray3DRef<double> vol,
          FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0,
          FortranArray3DRef<double> voli0,
          int& nummem,
          FortranArray3DRef<double> x, FortranArray3DRef<double> y,
          FortranArray3DRef<double> z,
          FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
          int& nbuf, int& ibufdim, int& iflagprnt)
{
    init_ns::init(nbl, jdim, kdim, idim, q, qj0, qk0, qi0, tj0, tk0, ti0,
                  vol, volj0, volk0, voli0, nummem, x, y, z,
                  nou, bou, nbuf, ibufdim, iflagprnt);
}

} // namespace findmin_new_ns
