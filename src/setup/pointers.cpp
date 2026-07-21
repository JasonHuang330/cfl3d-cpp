// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "pointers.h"
#include "runtime/fortran_io.h"
#include "module_kwstm.h"
#include "termn8.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

namespace pointers_ns {

// lead: reads back per-block parameters from lw2 and sets COMMON block variables
// This subroutine is the inverse of the lw2 assignments in pointers().
// It extracts ivisc, idiag, iflim, ifds, rkap from lw2 for block nbl.
void lead(int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, int& maxbl)
{
    // /reyue/ ivisc(3)
    cmn_reyue.ivisc[0] = lw2(14, nbl);
    cmn_reyue.ivisc[1] = lw2(15, nbl);
    cmn_reyue.ivisc[2] = lw2(16, nbl);

    // /info/ idiag(3)
    cmn_info.idiag[0] = lw2(23, nbl);
    cmn_info.idiag[1] = lw2(24, nbl);
    cmn_info.idiag[2] = lw2(25, nbl);

    // /info/ iflim(3)
    cmn_info.iflim[0] = lw2(26, nbl);
    cmn_info.iflim[1] = lw2(27, nbl);
    cmn_info.iflim[2] = lw2(28, nbl);

    // /info/ rkap(3) — factor 1e6 must be consistent with pointers()
    cmn_info.rkap[0] = (float)(lw2(32, nbl) / 1.e6);
    cmn_info.rkap[1] = (float)(lw2(33, nbl) / 1.e6);
    cmn_info.rkap[2] = (float)(lw2(34, nbl) / 1.e6);
}

void pointers(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, int& maxl, int& lembed,
              int& nstart, int& nwork, int& mwork, int& maxbl, int& maxgr, int& maxseg,
              FortranArray3DRef<int> lwdat, FortranArray1DRef<int> levelg,
              FortranArray1DRef<int> igridg, FortranArray2DRef<int> iviscg,
              FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
              FortranArray1DRef<int> kdimg, FortranArray1DRef<int> nblcg,
              FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
              FortranArray1DRef<int> idefrm, FortranArray1DRef<int> nbci0,
              FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
              FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim,
              FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo,
              FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
              int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg,
              FortranArray1DRef<int> iemg, int& nblock, int& myhost, int& myid,
              FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou,
              FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
              FortranArray1DRef<int> nblfine, FortranArray1DRef<int> ilamlog,
              FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> klamlog,
              FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamhig,
              FortranArray1DRef<int> klamhig, FortranArray2DRef<int> idegg,
              FortranArray2DRef<int> iwfg, FortranArray2DRef<int> idiagg,
              FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg,
              FortranArray2DRef<double> rkap0g, FortranArray1DRef<int> jsg,
              FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg,
              FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg,
              FortranArray1DRef<int> ieg, FortranArray1DRef<int> memblock,
              int& icall, int& nmds, int& maxaes, int& mpihost)
{
    // COMMON block aliases
    int32_t& ivmx    = cmn_maxiv.ivmx;
    int32_t& level   = cmn_mgrd.level;
    int32_t& lglobal = cmn_mgrd.lglobal;
    int32_t* ivisc   = cmn_reyue.ivisc;
    (void)cmn_elastic.ndefrm; // ndefrm declared in COMMON but not used in this subroutine
    int32_t& naesrf  = cmn_elastic.naesrf;
    int32_t& ita     = cmn_unst.ita;
    float&   dt      = cmn_info.dt;
    int32_t& mgflag  = cmn_info.mgflag;
    int32_t& iteravg = cmn_avgdata.iteravg;
    int32_t& ipertavg= cmn_avgdata.ipertavg;
    int32_t& ikoprod = cmn_konew.ikoprod;
    int32_t& isstdenom = cmn_konew.isstdenom;
    int32_t& isst2003  = cmn_konew.isst2003;
    int32_t& isstrc    = cmn_curvat.isstrc;
    int32_t& isstsf    = cmn_konew.isstsf;
    int32_t& isarc2d   = cmn_curvat.isarc2d;
    int32_t& isarc3d   = cmn_curvat.isarc3d;
    int32_t& isar      = cmn_curvat.isar;
    int32_t& i_wilcox06= cmn_konew.i_wilcox06;
    int32_t& i_wilcox98= cmn_konew.i_wilcox98;
    int32_t& ides      = cmn_des.ides;
    int32_t& i_nonlin  = cmn_constit.i_nonlin;
    int32_t& i_tauijs  = cmn_constit.i_tauijs;
    int32_t& lowmem_ux = cmn_memry.lowmem_ux;

    // local variables
    int ierrflg;
    int igrid, nbl, nblc, iaug, nfiner;
    int jmxx, kmxx, imxx;
    int jmkmim, jm1km1im1, jmkmim1, kmim1, jmim1, jmkm, kmim, jmim;
    int jmkmim_ig, jm1km1im1_ig, jmkmim1_ig, kmim1_ig, jmim1_ig;
    int jmkm_ig, kmim_ig, jmim_ig;
    int jbcmax = 0, kbcmax = 0, ibcmax = 0;
    int jbc, kbc, ibc;
    int nseg;
    int j, k, i, j1, k1, i1;
    int iuns, iem;
    int ns, ns1;
    int ii, jj, nf;
    int jdim, kdim, idim;
    float ftot_w, ftot_wk;



    if (icall == 0) {
        ierrflg = -99;
    } else {
        ierrflg = -1;
    }

    lglobal = 1;
    lembed  = 0;
    for (igrid = 1; igrid <= ngrid; igrid++) {
        lglobal = std::max(lglobal, (int)(ncgg(igrid) + 1));
        lembed  = std::max(lembed,  (int)iemg(igrid));
    }
    maxl = lglobal + lembed;

    if (myid == myhost && mpihost == 1) {

        jmxx = jdimg(nblg(1));
        kmxx = kdimg(nblg(1));
        imxx = idimg(nblg(1));
        jmkmim    = jmxx * kmxx * imxx;
        jm1km1im1 = (jmxx - 1) * (kmxx - 1) * (imxx - 1);
        jmkmim1   = jmxx * kmxx * (imxx - 1);
        kmim1     = kmxx * (imxx - 1);
        jmim1     = jmxx * (imxx - 1);
        jmkm      = jmxx * kmxx;
        kmim      = kmxx * imxx;
        jmim      = jmxx * imxx;

        for (igrid = 2; igrid <= ngrid; igrid++) {
            nbl  = nblg(igrid);
            jmxx = jdimg(nblg(igrid));
            kmxx = kdimg(nblg(igrid));
            imxx = idimg(nblg(igrid));
            jmkmim_ig    = jmxx * kmxx * imxx;
            jm1km1im1_ig = (jmxx - 1) * (kmxx - 1) * (imxx - 1);
            jmkmim1_ig   = jmxx * kmxx * (imxx - 1);
            kmim1_ig     = kmxx * (imxx - 1);
            jmim1_ig     = jmxx * (imxx - 1);
            jmkm_ig      = jmxx * kmxx;
            kmim_ig      = kmxx * imxx;
            jmim_ig      = jmxx * imxx;
            if (jmkmim_ig    > jmkmim)    jmkmim    = jmkmim_ig;
            if (jm1km1im1_ig > jm1km1im1) jm1km1im1 = jm1km1im1_ig;
            if (jmkmim1_ig   > jmkmim1)   jmkmim1   = jmkmim1_ig;
            if (kmim1_ig     > kmim1)     kmim1     = kmim1_ig;
            if (jmim1_ig     > jmim1)     jmim1     = jmim1_ig;
            if (jmkm_ig      > jmkm)      jmkm      = jmkm_ig;
            if (kmim_ig      > kmim)      kmim      = kmim_ig;
            if (jmim_ig      > jmim)      jmim      = jmim_ig;
        }

        // determine max boundary data array dimensions
        jbcmax = 0;
        kbcmax = 0;
        ibcmax = 0;
        for (igrid = 1; igrid <= ngrid; igrid++) {
            nbl = nblg(igrid);

            jbc = 0;
            for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
                if (jbcinfo(nbl, nseg, 7, 1) != 0)
                    jbc = jbc +
                        (jbcinfo(nbl, nseg, 3, 1) - jbcinfo(nbl, nseg, 2, 1)) *
                        (jbcinfo(nbl, nseg, 5, 1) - jbcinfo(nbl, nseg, 4, 1));
            }
            for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
                if (jbcinfo(nbl, nseg, 7, 2) != 0)
                    jbc = jbc +
                        (jbcinfo(nbl, nseg, 3, 2) - jbcinfo(nbl, nseg, 2, 2)) *
                        (jbcinfo(nbl, nseg, 5, 2) - jbcinfo(nbl, nseg, 4, 2));
            }
            if (jbc > jbcmax) jbcmax = jbc;

            kbc = 0;
            for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
                if (kbcinfo(nbl, nseg, 7, 1) != 0)
                    kbc = kbc +
                        (kbcinfo(nbl, nseg, 3, 1) - kbcinfo(nbl, nseg, 2, 1)) *
                        (kbcinfo(nbl, nseg, 5, 1) - kbcinfo(nbl, nseg, 4, 1));
            }
            for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
                if (kbcinfo(nbl, nseg, 7, 2) != 0)
                    kbc = kbc +
                        (kbcinfo(nbl, nseg, 3, 2) - kbcinfo(nbl, nseg, 2, 2)) *
                        (kbcinfo(nbl, nseg, 5, 2) - kbcinfo(nbl, nseg, 4, 2));
            }
            if (kbc > kbcmax) kbcmax = kbc;

            ibc = 0;
            for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
                if (ibcinfo(nbl, nseg, 7, 1) != 0)
                    ibc = ibc +
                        (ibcinfo(nbl, nseg, 3, 1) - ibcinfo(nbl, nseg, 2, 1)) *
                        (ibcinfo(nbl, nseg, 5, 1) - ibcinfo(nbl, nseg, 4, 1));
            }
            for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
                if (ibcinfo(nbl, nseg, 7, 2) != 0)
                    ibc = ibc +
                        (ibcinfo(nbl, nseg, 3, 2) - ibcinfo(nbl, nseg, 2, 2)) *
                        (ibcinfo(nbl, nseg, 5, 2) - ibcinfo(nbl, nseg, 4, 2));
            }
            if (ibc > ibcmax) ibcmax = ibc;

        } // end igrid loop for max bc dims

    } // end if (myid == myhost && mpihost == 1)

    nstart = 1;



    // Main loop over all blocks
    for (nbl = 1; nbl <= nblock; nbl++) {

        ivisc[0] = iviscg(nbl, 1);
        ivisc[1] = iviscg(nbl, 2);
        ivisc[2] = iviscg(nbl, 3);
        igrid    = igridg(nbl);
        level    = levelg(nbl);
        iuns     = std::max({(int)itrans(nbl), (int)irotat(nbl), (int)idefrm(nbl)});
        iem      = iemg(igrid);

        // determine if there are any finer embedded blocks for global block nbl
        nfiner = 0;
        for (int mbl = 1; mbl <= nblock; mbl++) {
            if (levelg(mbl) > lglobal) {
                if (nblcg(mbl) == nbl) {
                    nfiner = nfiner + 1;
                    nblfine(nfiner) = mbl;
                }
            }
        }

        j  = jdimg(nbl);
        k  = kdimg(nbl);
        i  = idimg(nbl);
        j1 = j - 1;
        k1 = k - 1;
        i1 = i - 1;

        if (myid != myhost || (myid == myhost && mpihost == 0)) {
            jmkmim    = j * k * i;
            jm1km1im1 = j1 * k1 * i1;
            jmkmim1   = j * k * i1;
            kmim1     = k * i1;
            jmim1     = j * i1;
            jmkm      = j * k;
            kmim      = k * i;
            jmim      = j * i;
        }

        // augment pointers if host or if block is local to this node,
        // or if this block needs embedded grid data from a non-local block
        nblc = nblcg(nbl);
        iaug = 0;
        if (myid == myhost || myid == mblk2nd(nbl)) iaug = 1;
        if (myid == mblk2nd(nblc)) iaug = 1;
        if (nfiner > 0) {
            for (nf = 1; nf <= nfiner; nf++) {
                if (myid == mblk2nd(nblfine(nf))) iaug = 1;
            }
        }

        if (iaug > 0) {

            // q(jdim,kdim,idim,5) — flow variables at cell centers
            lw(1, nbl) = nstart;
            ns = jmkmim * 5;

            // qj0(kdim,idim-1,5,4) — flow variables at j-boundary
            lw(2, nbl) = lw(1, nbl) + ns;
            ns = kmim1 * 5 * 4;

            // qk0(jdim,idim-1,5,4) — flow variables at k-boundary
            lw(3, nbl) = lw(2, nbl) + ns;
            ns = jmim1 * 5 * 4;

            // qi0(jdim,kdim,5,4) — flow variables at i-boundary
            lw(4, nbl) = lw(3, nbl) + ns;
            ns = jmkm * 5 * 4;

            // sj(jdim,kdim,idim-1,5) — j-face metrics
            lw(5, nbl) = lw(4, nbl) + ns;
            ns = jmkmim1 * 5;

            // sk(jdim,kdim,idim-1,5) — k-face metrics
            lw(6, nbl) = lw(5, nbl) + ns;
            ns = jmkmim1 * 5;

            // si(jdim,kdim,idim,5) — i-face metrics
            lw(7, nbl) = lw(6, nbl) + ns;
            ns = jmkmim * 5;

            // vol(jdim,kdim,idim-1) — cell volume
            lw(8, nbl) = lw(7, nbl) + ns;
            ns = jmkmim1;

            // dtj(jdim,kdim,idim-1) — cell time step
            lw(9, nbl) = lw(8, nbl) + ns;
            ns = jmkmim1;

            // x(jdim,kdim,idim) — grid point x coordinate
            lw(10, nbl) = lw(9, nbl) + ns;
            ns = jmkmim;

            // y(jdim,kdim,idim) — grid point y coordinate
            lw(11, nbl) = lw(10, nbl) + ns;
            ns = jmkmim;

            // z(jdim,kdim,idim) — grid point z coordinate
            lw(12, nbl) = lw(11, nbl) + ns;
            ns = jmkmim;

            // vist3d(jdim,kdim,idim) — eddy viscosity at cell center
            lw(13, nbl) = lw(12, nbl) + ns;
            ns = 0;
            if (ivmx > 1) ns = jmkmim;

            // snk0(jdim-1,kdim-1,idim-1) — distance to nearest solid surface
            lw(14, nbl) = lw(13, nbl) + ns;
            ns = 0;
            if (ivisc[2] > 1 || ivisc[1] >= 4 || ivisc[0] >= 4)
                ns = jm1km1im1;

            // sni0/snj0/xjb(jdim-1,kdim-1,idim-1)
            lw(15, nbl) = lw(14, nbl) + ns;
            ns = 0;
            if (ivisc[0] == 2 || ivisc[1] == 2 || ivisc[2] == 4 ||
                ivisc[1] == 4 || ivisc[0] == 4 || ivisc[2] == 25 ||
                ivisc[1] == 25 || ivisc[0] == 25)
                ns = jm1km1im1;



            // q1(jdim,kdim,idim,5) — cell center flow variables restricted from finer grid
            lw(16, nbl) = lw(15, nbl) + ns;
            ns = jmkmim * 5;
            if (mgflag == 0) ns = 0;
            if (iem == 0) {
                // global grid
                if (mgflag == 1 && level == lglobal) ns = 0;
            } else {
                // embedded grid
                if (mgflag != 2 || level == maxl) ns = 0;
            }

            // qr(jdim,kdim,idim-1,5) — cell center residuals restricted from finer grid
            lw(17, nbl) = lw(16, nbl) + ns;
            ns = jmkmim1 * 5;

            if (iem == 0) {
                // global grid
                if (mgflag == 0)      ns = 0;
                if (level == lglobal) ns = 0;
            }

            // blank(jdim,kdim,idim) — cell center blank values
            lw(18, nbl) = lw(17, nbl) + ns;
            ns = jmkmim;

            // xib/tursav(jdim,kdim,idim,2) — turbulence model parameters
            lw(19, nbl) = lw(18, nbl) + ns;
            ns = 0;
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4)
                ns = jmkmim * 2;
            // extra memory needed for ivisc=30 (3-eqn model)
            if (ivisc[2] == 30 || ivisc[1] == 30 || ivisc[0] == 30)
                ns = jmkmim * 3;
            // extra memory needed for ivisc=40 (4-eqn model)
            if (ivisc[2] == 40 || ivisc[1] == 40 || ivisc[0] == 40)
                ns = jmkmim * 4;
            // extra memory needed for ivisc=72 (stress-omega model, Wilcox 3rd Edition)
            if (ivisc[2] == 72 || ivisc[1] == 72 || ivisc[0] == 72)
                ns = jmkmim * (int)module_kwstm_ns::kws_get_nummem();

            // Johnson-King data no longer stored
            lw(20, nbl) = lw(19, nbl) + ns;
            ns = 0;
            lw(21, nbl) = lw(20, nbl) + ns;
            ns = 0;
            lw(22, nbl) = lw(21, nbl) + ns;
            ns = 0;

            // tj0(kdim,idim-1,2,4) — turbulence model variable at j-boundary
            lw(23, nbl) = lw(22, nbl) + ns;
            ns = 0;
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4)
                ns = kmim1 * 2 * 4;
            // extra memory needed for ivisc=30 (3-eqn model)
            if (ivisc[2] == 30 || ivisc[1] == 30 || ivisc[0] == 30)
                ns = kmim1 * 3 * 4;
            // extra memory needed for ivisc=40 (4-eqn model)
            if (ivisc[2] == 40 || ivisc[1] == 40 || ivisc[0] == 40)
                ns = kmim1 * 4 * 4;
            // extra memory needed for ivisc=72 (stress-omega model)
            if (ivisc[2] == 72 || ivisc[1] == 72 || ivisc[0] == 72)
                ns = kmim1 * (int)module_kwstm_ns::kws_get_nummem() * 4;

            // tk0(jdim,idim-1,2,4) — turbulence model variable at k-boundary
            lw(24, nbl) = lw(23, nbl) + ns;
            ns = 0;
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4)
                ns = jmim1 * 2 * 4;
            // extra memory needed for ivisc=30 (3-eqn model)
            if (ivisc[2] == 30 || ivisc[1] == 30 || ivisc[0] == 30)
                ns = jmim1 * 3 * 4;
            // extra memory needed for ivisc=40 (4-eqn model)
            if (ivisc[2] == 40 || ivisc[1] == 40 || ivisc[0] == 40)
                ns = jmim1 * 4 * 4;
            // extra memory needed for ivisc=72 (4-eqn model)
            if (ivisc[2] == 72 || ivisc[1] == 72 || ivisc[0] == 72)
                ns = jmim1 * (int)module_kwstm_ns::kws_get_nummem() * 4;

            // ti0(jdim,kdim,2,4) — turbulence model variable at i-boundary
            lw(25, nbl) = lw(24, nbl) + ns;
            ns = 0;
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4)
                ns = jmkm * 2 * 4;
            // extra memory needed for ivisc=30 (3-eqn model)
            if (ivisc[2] == 30 || ivisc[1] == 30 || ivisc[0] == 30)
                ns = jmkm * 3 * 4;
            // extra memory needed for ivisc=40 (4-eqn model)
            if (ivisc[2] == 40 || ivisc[1] == 40 || ivisc[0] == 40)
                ns = jmkm * 4 * 4;
            // extra memory needed for ivisc=72 (stress-omega model)
            if (ivisc[2] == 72 || ivisc[1] == 72 || ivisc[0] == 72)
                ns = jmkm * (int)module_kwstm_ns::kws_get_nummem() * 4;

            // xkb(jdim-1,kdim-1,idim-1) — k-index of nearest surface point
            lw(26, nbl) = lw(25, nbl) + ns;
            ns = 0;
            if (ivisc[2] == 2 ||
                ivisc[2] == 4 || ivisc[1] == 4 || ivisc[0] == 4 ||
                ivisc[2] == 25 || ivisc[1] == 25 || ivisc[0] == 25)
                ns = jm1km1im1;

            // blnum(jdim-1,kdim-1,idim-1) — block number containing nearest surface point
            lw(27, nbl) = lw(26, nbl) + ns;
            ns = 0;
            if (ivisc[1] == 2 || ivisc[0] == 2 ||
                ivisc[2] == 4 || ivisc[1] == 4 || ivisc[0] == 4 ||
                ivisc[2] == 25 || ivisc[1] == 25 || ivisc[0] == 25)
                ns = jm1km1im1;

            // vj0(kdim,idim-1,1,4) — eddy viscosity at j-boundary
            lw(28, nbl) = lw(27, nbl) + ns;
            ns = 0;
            if (ivmx > 1) ns = kmim1 * 1 * 4;

            // vk0(jdim,idim-1,1,4) — eddy viscosity at k-boundary
            lw(29, nbl) = lw(28, nbl) + ns;
            ns = 0;
            if (ivmx > 1) ns = jmim1 * 1 * 4;

            // vi0(jdim,kdim,1,4) — eddy viscosity at i-boundary
            lw(30, nbl) = lw(29, nbl) + ns;
            ns = 0;
            if (ivmx > 1) ns = jmkm * 1 * 4;

            // bcj(kdim,idim-1,2) — boundary data type flag at j-boundary
            lw(31, nbl) = lw(30, nbl) + ns;
            ns         = kmim1 * 2;

            // bck(jdim,idim-1,2) — boundary data type flag at k-boundary
            lw(32, nbl) = lw(31, nbl) + ns;
            ns         = jmim1 * 2;

            // bci(jdim,kdim,2) — boundary data type flag at i-boundary
            lw(33, nbl) = lw(32, nbl) + ns;
            ns         = jmkm * 2;

            // qc0(jdim,kdim,idim-1,5) — flow variables at previous time step/subiteration
            lw(34, nbl) = lw(33, nbl) + ns;
            ns = 0;
            if ((float)dt > 0.) ns = jmkmim1 * 5;

            // dqc0(jdim,kdim,idim-1,5) — change in flow variables from previous time step
            lw(35, nbl) = lw(34, nbl) + ns;
            ns = 0;
            if ((float)dt > 0.) ns = jmkmim1 * 5;

            // xtbj(kdim,idim-1,3,2) — j-boundary velocity components on moving grids
            lw(36, nbl) = lw(35, nbl) + ns;
            ns = 0;
            if (iuns > 0) ns = kmim1 * 3 * 2;

            // xtbk(jdim,idim-1,3,2) — k-boundary velocity components on moving grids
            lw(37, nbl) = lw(36, nbl) + ns;
            ns = 0;
            if (iuns > 0) ns = jmim1 * 3 * 2;

            // xtbi(jdim,kdim,3,2) — i-boundary velocity components on moving grids
            lw(38, nbl) = lw(37, nbl) + ns;
            ns = 0;
            if (iuns > 0) ns = jmkm * 3 * 2;

            // atbj(kdim,idim-1,3,2) — j-boundary acceleration components on moving grids
            lw(39, nbl) = lw(38, nbl) + ns;
            ns = 0;
            if (iuns > 0) ns = kmim1 * 3 * 2;

            // atbk(jdim,idim-1,3,2) — k-boundary acceleration components on moving grids
            lw(40, nbl) = lw(39, nbl) + ns;
            ns = 0;
            if (iuns > 0) ns = jmim1 * 3 * 2;

            // atbi(jdim,kdim,3,2) — i-boundary acceleration components on moving grids
            lw(41, nbl) = lw(40, nbl) + ns;
            ns = 0;
            if (iuns > 0) ns = jmkm * 3 * 2;



            // bcdataj(npts,mpts,2,12) — 2000 series bc data for j-boundaries
            lw(42, nbl) = lw(41, nbl) + ns;

            ns  = 0;
            ns1 = 0;
            for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
                lwdat(nbl, nseg, 3) = lw(42, nbl) + ns1;
                if (jbcinfo(nbl, nseg, 7, 1) != 0)
                    ns1 = ns1 +
                        (jbcinfo(nbl, nseg, 3, 1) - jbcinfo(nbl, nseg, 2, 1)) *
                        (jbcinfo(nbl, nseg, 5, 1) - jbcinfo(nbl, nseg, 4, 1)) * 12 * 2;
            }
            for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
                lwdat(nbl, nseg, 4) = lw(42, nbl) + ns1;
                if (jbcinfo(nbl, nseg, 7, 2) != 0)
                    ns1 = ns1 +
                        (jbcinfo(nbl, nseg, 3, 2) - jbcinfo(nbl, nseg, 2, 2)) *
                        (jbcinfo(nbl, nseg, 5, 2) - jbcinfo(nbl, nseg, 4, 2)) * 12 * 2;
            }

            ns = ns + ns1;
            if (myid == myhost && mpihost == 1) {
                ns = jbcmax * 12 * 2;
            }

            // bcdatak(npts,mpts,2,12) — 2000 series bc data for k-boundaries
            lw(43, nbl) = lw(42, nbl) + ns;

            ns  = 0;
            ns1 = 0;
            for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
                lwdat(nbl, nseg, 5) = lw(43, nbl) + ns1;
                if (kbcinfo(nbl, nseg, 7, 1) != 0)
                    ns1 = ns1 +
                        (kbcinfo(nbl, nseg, 3, 1) - kbcinfo(nbl, nseg, 2, 1)) *
                        (kbcinfo(nbl, nseg, 5, 1) - kbcinfo(nbl, nseg, 4, 1)) * 12 * 2;
            }
            for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
                lwdat(nbl, nseg, 6) = lw(43, nbl) + ns1;
                if (kbcinfo(nbl, nseg, 7, 2) != 0)
                    ns1 = ns1 +
                        (kbcinfo(nbl, nseg, 3, 2) - kbcinfo(nbl, nseg, 2, 2)) *
                        (kbcinfo(nbl, nseg, 5, 2) - kbcinfo(nbl, nseg, 4, 2)) * 12 * 2;
            }

            ns = ns + ns1;
            if (myid == myhost && mpihost == 1) {
                ns = kbcmax * 12 * 2;
            }

            // bcdatai(npts,mpts,2,12) — 2000 series bc data for i-boundaries
            lw(44, nbl) = lw(43, nbl) + ns;

            ns  = 0;
            ns1 = 0;
            for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
                lwdat(nbl, nseg, 1) = lw(44, nbl) + ns1;
                if (ibcinfo(nbl, nseg, 7, 1) != 0)
                    ns1 = ns1 +
                        (ibcinfo(nbl, nseg, 3, 1) - ibcinfo(nbl, nseg, 2, 1)) *
                        (ibcinfo(nbl, nseg, 5, 1) - ibcinfo(nbl, nseg, 4, 1)) * 12 * 2;
            }
            for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
                lwdat(nbl, nseg, 2) = lw(44, nbl) + ns1;
                if (ibcinfo(nbl, nseg, 7, 2) != 0)
                    ns1 = ns1 +
                        (ibcinfo(nbl, nseg, 3, 2) - ibcinfo(nbl, nseg, 2, 2)) *
                        (ibcinfo(nbl, nseg, 5, 2) - ibcinfo(nbl, nseg, 4, 2)) * 12 * 2;
            }

            ns = ns + ns1;
            if (myid == myhost && mpihost == 1) {
                ns = ibcmax * 12 * 2;
            }

            // xib2(jdim,kdim,idim,4) — extra storage for subiterations for time-accurate turb. models
            lw(45, nbl) = lw(44, nbl) + ns;
            ns = 0;
            if ((ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) &&
                (float)dt > 0.) ns = jmkmim * 4;
            // extra memory needed for ivisc=30 (3-eqn model)
            if (ivisc[2] == 30 || ivisc[1] == 30 || ivisc[0] == 30)
                ns = jmkmim * 6;
            // extra memory needed for ivisc=40 (4-eqn model)
            if (ivisc[2] == 40 || ivisc[1] == 40 || ivisc[0] == 40)
                ns = jmkmim * 8;
            // extra memory needed for ivisc=72 (stress-omega model)
            if (ivisc[2] == 72 || ivisc[1] == 72 || ivisc[0] == 72)
                ns = jmkmim * 14;

            // qavg(jdim,kdim,idim,5) — iteration-averaged data
            lw(46, nbl) = lw(45, nbl) + ns;
            ns = 0;
            if (iteravg == 1 || iteravg == 2) ns = jmkmim * 5;

            // ux(jdim-1,kdim-1,idim-1,9) — velocity gradients
            lw(47, nbl) = lw(46, nbl) + ns;
            ns = 0;
            if (lowmem_ux == 1) {
                if (ivmx == 8 || ivmx == 9 || ivmx >= 11 ||
                    ivmx == 16 ||
                    ((ivmx == 6 || ivmx == 7 || ivmx == 10) && ikoprod == 1) ||
                    (ivmx == 7 && isstdenom == 1) ||
                    (ivmx == 7 && isst2003 == 1) ||
                    ((ivmx == 6 || ivmx == 7) && isstrc > 0) ||
                    ((ivmx == 6 || ivmx == 7) && isstsf == 1) ||
                    ((ivmx == 5) && isarc2d == 1) ||
                    ((ivmx == 5) && isarc3d == 1) ||
                    ((ivmx == 5) && isar == 1) ||
                    (ivmx == 6 && i_wilcox06 == 1) ||
                    (ivmx == 6 && i_wilcox98 == 1) || ivmx == 72 ||
                    ides >= 2 || i_nonlin != 0 || i_tauijs != 0)
                    ns = jm1km1im1 * 9;
            } else {
                ns = jm1km1im1 * 9;
            }



            // cmuv(jdim-1,kdim-1,idim-1) — EASM, var-g
            lw(48, nbl) = lw(47, nbl) + ns;
            ns = 0;
            if (ivisc[2] == 8 || ivisc[1] == 8 || ivisc[0] == 8 ||
                ivisc[2] == 9 || ivisc[1] == 9 || ivisc[0] == 9 ||
                ivisc[2] == 13 || ivisc[1] == 13 || ivisc[0] == 13 ||
                ivisc[2] == 14 || ivisc[1] == 14 || ivisc[0] == 14)
                ns = jm1km1im1;

            // volj0(kdim,idim-1,1,4) — cell volumes at j-boundary
            lw(49, nbl) = lw(48, nbl) + ns;
            ns = 0;
            if (ivmx >= 1) ns = kmim1 * 1 * 4;

            // volk0(jdim,idim-1,1,4) — cell volumes at k-boundary
            lw(50, nbl) = lw(49, nbl) + ns;
            ns = 0;
            if (ivmx >= 1) ns = jmim1 * 1 * 4;

            // voli0(jdim,kdim,1,4) — cell volumes at i-boundary
            lw(51, nbl) = lw(50, nbl) + ns;
            ns = 0;
            if (ivmx >= 1) ns = jmkm * 1 * 4;

            // aeroelastic j-surface data
            lw(52, nbl) = lw(51, nbl) + ns;
            ns = 0;
            if (naesrf > 0) ns = kmim * 6 * maxaes * nmds;

            // aeroelastic k-surface data
            lw(53, nbl) = lw(52, nbl) + ns;
            ns = 0;
            if (naesrf > 0) ns = jmim * 6 * maxaes * nmds;

            // aeroelastic i-surface data
            lw(54, nbl) = lw(53, nbl) + ns;
            ns = 0;
            if (naesrf > 0) ns = jmkm * 6 * maxaes * nmds;

            // intermediate grid velocity — Note: this data no longer stored
            lw(55, nbl) = lw(54, nbl) + ns;
            ns = 0;

            // x-coordinate at time step n-2
            lw(56, nbl) = lw(55, nbl) + ns;
            ns = 0;
            if (idefrm(nbl) > 0 && std::abs(ita) > 1) ns = jmkmim;

            // y-coordinate at time step n-2
            lw(57, nbl) = lw(56, nbl) + ns;
            ns = 0;
            if (idefrm(nbl) > 0 && std::abs(ita) > 1) ns = jmkmim;

            // z-coordinate at time step n-2
            lw(58, nbl) = lw(57, nbl) + ns;
            ns = 0;
            if (idefrm(nbl) > 0 && std::abs(ita) > 1) ns = jmkmim;

            // displacements of j=const surfaces for deforming grids
            lw(59, nbl) = lw(58, nbl) + ns;
            ns = 0;
            if (idefrm(nbl) > 0) ns = kmim * 3 * 2;

            // displacements of k=const surfaces for deforming grids
            lw(60, nbl) = lw(59, nbl) + ns;
            ns = 0;
            if (idefrm(nbl) > 0) ns = jmim * 3 * 2;

            // displacements of i=const surfaces for deforming grids
            lw(61, nbl) = lw(60, nbl) + ns;
            ns = 0;
            if (idefrm(nbl) > 0) ns = jmkm * 3 * 2;

            // x-coordinate at time step n-1
            lw(62, nbl) = lw(61, nbl) + ns;
            ns = 0;
            if (idefrm(nbl) > 0) ns = jmkmim;

            // y-coordinate at time step n-1
            lw(63, nbl) = lw(62, nbl) + ns;
            ns = 0;
            if (idefrm(nbl) > 0) ns = jmkmim;

            // z-coordinate at time step n-1
            lw(64, nbl) = lw(63, nbl) + ns;
            ns = 0;
            if (idefrm(nbl) > 0) ns = jmkmim;

            // q2avg(jdim,kdim,idim,5) — Average of the square of q
            lw(65, nbl) = lw(64, nbl) + ns;
            ns = 0;
            if (ipertavg == 1 || ipertavg == 2) ns = jmkmim * 5;

            nstart = lw(65, nbl) + ns;

            if (myid == myhost && mpihost == 1) {
                if (nbl < nblock) {
                    nstart = 1;
                }
            }

            // memblock is the memory needed to store the w array for this block
            memblock(nbl) = lw(65, nbl) + ns - lw(1, nbl) + 1;

        } else {

            // if block is not local to this node, set pointers to previous block's starting location
            for (ii = 1; ii <= 65; ii++) {
                lw(ii, nbl) = nstart;
            }
            for (ii = 1; ii <= maxseg; ii++) {
                for (jj = 1; jj <= 6; jj++) {
                    lwdat(nbl, ii, jj) = nstart;
                }
            }

            nstart = lw(65, nbl);

        } // end if (iaug > 0)



        // set up lw2 array
        jdim        = jdimg(nbl);
        kdim        = kdimg(nbl);
        idim        = idimg(nbl);
        lw2(1,  nbl) = jdimg(nbl);
        lw2(2,  nbl) = kdimg(nbl);
        lw2(3,  nbl) = idimg(nbl);
        lw2(4,  nbl) = nblcg(nbl);
        lw2(5,  nbl) = jdim / 2 + 1;
        lw2(6,  nbl) = kdim / 2 + 1;
        lw2(7,  nbl) = idim / 2 + 1;
        lw2(8,  nbl) = ilamlog(nbl);
        lw2(9,  nbl) = ilamhig(nbl);
        lw2(10, nbl) = jlamlog(nbl);
        lw2(11, nbl) = jlamhig(nbl);
        lw2(12, nbl) = klamlog(nbl);
        lw2(13, nbl) = klamhig(nbl);
        lw2(14, nbl) = iviscg(nbl, 1);
        lw2(15, nbl) = iviscg(nbl, 2);
        lw2(16, nbl) = iviscg(nbl, 3);
        lw2(17, nbl) = idegg(nbl, 1);
        lw2(18, nbl) = idegg(nbl, 2);
        lw2(19, nbl) = idegg(nbl, 3);
        lw2(20, nbl) = iwfg(nbl, 1);
        lw2(21, nbl) = iwfg(nbl, 2);
        lw2(22, nbl) = iwfg(nbl, 3);
        lw2(23, nbl) = idiagg(nbl, 1);
        lw2(24, nbl) = idiagg(nbl, 2);
        lw2(25, nbl) = idiagg(nbl, 3);
        lw2(26, nbl) = iflimg(nbl, 1);
        lw2(27, nbl) = iflimg(nbl, 2);
        lw2(28, nbl) = iflimg(nbl, 3);
        lw2(29, nbl) = ifdsg(nbl, 1);
        lw2(30, nbl) = ifdsg(nbl, 2);
        lw2(31, nbl) = ifdsg(nbl, 3);
        // factor multiplying rkap0g must be consistent with subroutine lead!
        lw2(32, nbl) = (int)(rkap0g(nbl, 1) * 1.e6);
        lw2(33, nbl) = (int)(rkap0g(nbl, 2) * 1.e6);
        lw2(34, nbl) = (int)(rkap0g(nbl, 3) * 1.e6);
        lw2(35, nbl) = jsg(nbl);
        lw2(36, nbl) = ksg(nbl);
        lw2(37, nbl) = isg(nbl);
        lw2(38, nbl) = jeg(nbl);
        lw2(39, nbl) = keg(nbl);
        lw2(40, nbl) = ieg(nbl);
        if (iemg(igridg(nbl)) > 0) {
            nblc        = nblcg(nbl);
            lw2(5, nbl) = jdimg(nblc);
            lw2(6, nbl) = kdimg(nblc);
            lw2(7, nbl) = idimg(nblc);
        }

    } // end do 1000 nbl=1,nblock

    ftot_w  = (float)nstart / 1.e+06f;
    nwork   = mwork - nstart;
    ftot_wk = (float)nwork / 1.e+06f;

    if (icall > 0) {
        FILE* f11 = fortran_get_unit(11);
        fprintf(f11, "\n");
        fprintf(f11, " SUMMARY OF PRIMARY MEMORY REQUIREMENTS\n");
        fprintf(f11, "  memory (mw) for w  storage =%12.5E\n", (float)ftot_w);
        fprintf(f11, "  memory (mw) for wk storage =%12.5E\n", (float)ftot_wk);
        if ((float)ftot_wk <= 0.) {
            fprintf(f11, "  stopping...mwork too small\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }

} // end pointers

} // namespace pointers_ns
