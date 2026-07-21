// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "qoutavg.h"
#include "runtime/fortran_io.h"
#include "lead.h"
#include "bc.h"
#include "bc_period.h"
#include "bc_embed.h"
#include "bc_blkint.h"
#include "bc_patch.h"
#include "bc_xmera.h"
#include "qface.h"
#include "blnkfr.h"
#include "plot3d.h"
#include "plot3c.h"
#include "plot3t.h"
#include "plot3davg.h"
#include "termn8.h"
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cmath>

namespace qoutavg_ns {

void qoutavg(int& iseq, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
             FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
             int& nwork, int& nplots, FortranArray1DRef<int> iovrlp,
             FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg,
             FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg,
             FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg,
             FortranArray3DRef<double> qb, FortranArray3DRef<int> lwdat,
             FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0,
             FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim,
             FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim,
             FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
             FortranArray4DRef<int> ibcinfo, FortranArray3DRef<int> bcfilei,
             FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek,
             FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
             FortranArray1DRef<int> idefrm, int& nblock, FortranArray1DRef<int> levelg,
             FortranArray1DRef<int> igridg, FortranArray2DRef<int> iviscg,
             FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg,
             FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblg,
             FortranArray1DRef<double> clw, int& ncycmax, int& nplot3d,
             FortranArray2DRef<int> inpl3d, int& ip3dsurf, int& nprint,
             FortranArray2DRef<int> inpr, FortranArray1DRef<int> iadvance,
             int& mycomm, int& myid, int& myhost, FortranArray1DRef<int> mblk2nd,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg,
             int& iitot, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg,
             FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg,
             FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg,
             int& ninter, FortranArray2DRef<double> windex,
             FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt,
             int& intmax, int& nsub1, int& maxxe, FortranArray2DRef<int> nblk,
             int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva,
             FortranArray1DRef<int> nblon, int& mxbli, FortranArray1DRef<double> thetay,
             FortranArray2DRef<int> ip3ddim, FortranArray1DRef<int> nmap,
             FortranArray1DRef<int> iwk, int& iwork, FortranArray1DRef<double> xorig,
             FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig,
             FortranArray1DRef<double> period_miss, FortranArray1DRef<double> geom_miss,
             double& epsc0, double& epsrot, FortranArray2DRef<int> isav_blk,
             FortranArray2DRef<int> isav_pat, FortranArray3DRef<int> isav_pat_b,
             FortranArray2DRef<int> isav_emb, FortranArray2DRef<int> isav_prd,
             int& lbcprd, int& lbcemb, FortranArray2DRef<double> dthetxx,
             FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz,
             FortranArray1DRef<int> nblcg, int& lfgm,
             FortranArray2DRef<int> istat2_bl, FortranArray2DRef<int> istat2_pa,
             FortranArray2DRef<int> istat2_pe, FortranArray2DRef<int> istat2_em,
             int& istat_size, FortranArray1DRef<double> vormax,
             FortranArray1DRef<int> ivmax, FortranArray1DRef<int> jvmax,
             FortranArray1DRef<int> kvmax, FortranArray1DRef<char[80]> bcfiles,
             int& mxbcfil, int& iprnsurf, int& nt, int& movabs, int& nummem)
{
    // COMMON block aliases
    int& ibin      = cmn_bin.ibin;
    int& iblnk     = cmn_bin.iblnk;
    int& iblnkfr   = cmn_bin.iblnkfr;
    int& ip3dgrad  = cmn_bin.ip3dgrad;

    int& jdim      = cmn_ginfo.jdim;
    int& kdim      = cmn_ginfo.kdim;
    int& idim      = cmn_ginfo.idim;
    int& lq        = cmn_ginfo.lq;
    int& lqj0      = cmn_ginfo.lqj0;
    int& lqk0      = cmn_ginfo.lqk0;
    int& lqi0      = cmn_ginfo.lqi0;
    int& lsj       = cmn_ginfo.lsj;
    int& lsk       = cmn_ginfo.lsk;
    int& lsi       = cmn_ginfo.lsi;
    int& lvol      = cmn_ginfo.lvol;
    int& lx        = cmn_ginfo.lx;
    int& ly        = cmn_ginfo.ly;
    int& lz        = cmn_ginfo.lz;
    int& lvis      = cmn_ginfo.lvis;
    int& lsnk0     = cmn_ginfo.lsnk0;
    int& lblk      = cmn_ginfo.lblk;
    int& lbcj      = cmn_ginfo.lbcj;
    int& lbck      = cmn_ginfo.lbck;
    int& lbci      = cmn_ginfo.lbci;
    int& lvj0      = cmn_ginfo.lvj0;
    int& lvk0      = cmn_ginfo.lvk0;
    int& lvi0      = cmn_ginfo.lvi0;
    int& lux       = cmn_ginfo.lux;
    int& lxib      = cmn_ginfo.lxib;
    int& lcmuv     = cmn_ginfo.lcmuv;
    int& lvolj0    = cmn_ginfo.lvolj0;
    int& lvolk0    = cmn_ginfo.lvolk0;
    int& lvoli0    = cmn_ginfo.lvoli0;
    int& lqavg     = cmn_ginfo.lqavg;
    int& lq2avg    = cmn_ginfo2.lq2avg;

    int& mseq      = cmn_info.mseq;
    int& iipv      = cmn_info.iipv;
    int& ntt       = cmn_info.ntt;
    float* title   = cmn_info.title;
    float& xmach   = cmn_info.xmach;
    float& alpha   = cmn_info.alpha;
    float& beta    = cmn_info.beta;
    int* levelt    = cmn_info.levelt;

    int& ivmx      = cmn_maxiv.ivmx;
    int& level     = cmn_mgrd.level;
    int& lglobal   = cmn_mgrd.lglobal;
    int& i2d       = cmn_twod.i2d;
    int& iunst     = cmn_unst.iunst;
    float& time    = cmn_unst.time;
    int& icall1    = cmn_moov.icall1;
    int& lhdr      = cmn_moov.lhdr;
    float& radtodeg = cmn_conversion.radtodeg;
    int& iteravg   = cmn_avgdata.iteravg;
    int& ifunct    = cmn_plot3dtyp.ifunct;
    int& lowmem_ux = cmn_memry.lowmem_ux;
    float& reue    = cmn_reyue.reue;
    float& tinf    = cmn_reyue.tinf;

    // Local variables
    int level_sav;
    int nttuse;
    double clwuse = 0.0;
    int lres, nsafe, mneed, iwk1, iwk2, iwk3, iwk4, iwk5, iwk6, iwork1;
    int iii, nbl, n, nnn, m;
    int int_updt;
    int ncount, np3d, ifunc, ifuncuse;
    int i1, i2, i3, j1, j2, j3, k1, k2, k3;
    int jdw, kdw, idw;
    int ixwk, ibwk, ixgk, ixvk, ibwk2, ibwk3, nroom, nset;
    int ifuncdim;
    float alphaw, betaw;
    int iflag;



    // if (nplot3d.gt.0 .or. nprint.gt.0 .or. iteravg.gt.0) then
    if (nplot3d > 0 || nprint > 0 || iteravg > 0) {

        // update all bc's for level lglobal (and above, if embedded)
        level_sav = level;
        lglobal = lfgm - (mseq - iseq);

        for (level = lglobal; level <= levelt[iseq-1]; level++) {

            if (iipv > 0) {
                nttuse = std::max(ntt - 1, 1);
                clwuse = clw(nttuse);
            }

            for (nbl = 1; nbl <= nblock; nbl++) {
                // need to call even for blocks not advanced
                if (level == levelg(nbl)) {
                    if (mblk2nd(nbl) == myid) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        int ntime1 = 1;
                        bc_ns::bc(ntime1, nbl, lw, lw2, w, mgwk, wk, nwork,
                                  clwuse, nou, bou, nbuf, ibufdim, maxbl,
                                  maxgr, maxseg, itrans, irotat, idefrm,
                                  igridg, nblg, nbci0, nbcj0, nbck0, nbcidim,
                                  nbcjdim, nbckdim, ibcinfo, jbcinfo,
                                  kbcinfo, bcfilei, bcfilej, bcfilek,
                                  lwdat, myid, idimg, jdimg, kdimg, bcfiles,
                                  mxbcfil, nummem);
                    }
                }
            }

            // update periodic boundary conditions
            lres   = 1;
            nsafe  = nwork - lres + 1;
            mneed  = lbcprd * 5;
            iwk1   = 1;
            iwk2   = iwk1 + mneed;
            iwk3   = iwk2 + mneed;
            iwk4   = iwk3 + mneed;
            iwk5   = iwk4 + mneed * 2;
            iwk6   = iwk5 + mneed;
            iwork1 = iwork - iwk6;
            if (iwork1 < 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "stopping...not enough integer work space for subroutine bc_period");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                int ntime1 = 1;
                // iwk(iwk1..iwk5) are FortranArray1DRef<int> slices
                FortranArray1DRef<int> ireq_ar_ref(&iwk(iwk1), mneed);
                FortranArray1DRef<int> index_ar_ref(&iwk(iwk2), mneed);
                FortranArray1DRef<int> ireq_snd_ref(&iwk(iwk3), mneed);
                FortranArray2DRef<int> keep_trac_ref(&iwk(iwk4), mneed, 2);
                FortranArray1DRef<int> keep_trac2_ref(&iwk(iwk5), mneed);
                FortranArray1DRef<double> wk_lres(&wk(lres), nsafe);
                bc_period_ns::bc_period(ntime1, nbl, lw, lw2, w, mgwk, wk_lres, nsafe,
                                        maxbl, maxgr, maxseg, iadvance,
                                        bcfilei, bcfilej, bcfilek, lwdat,
                                        xorig, yorig, zorig, jdimg, kdimg,
                                        idimg, lbcprd, isav_prd,
                                        period_miss, epsrot,
                                        ireq_ar_ref, index_ar_ref,
                                        ireq_snd_ref, keep_trac_ref,
                                        keep_trac2_ref,
                                        myid, myhost, mycomm, mblk2nd,
                                        nou, bou, nbuf, ibufdim,
                                        istat2_pe, istat_size, bcfiles,
                                        mxbcfil, nummem);
            }

            // update embedded-grid boundary conditions
            lres   = 1;
            nsafe  = nwork - lres + 1;
            mneed  = lbcemb * 3;
            iwk1   = 1;
            iwk2   = iwk1 + mneed;
            iwk3   = iwk2 + mneed;
            iwk4   = iwk3 + mneed;
            iwk5   = iwk4 + mneed * 2;
            iwk6   = iwk5 + mneed;
            iwork1 = iwork - iwk6;
            if (iwork1 < 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "stopping...not enough integer work space for subroutine bc_embed");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                int ntime1 = 1;
                FortranArray1DRef<int> ireq_ar_ref(&iwk(iwk1), mneed);
                FortranArray1DRef<int> ireq_snd_ref(&iwk(iwk2), mneed);
                FortranArray1DRef<int> index_ar_ref(&iwk(iwk3), mneed);
                FortranArray2DRef<int> keep_trac_ref(&iwk(iwk4), mneed, 2);
                FortranArray1DRef<int> keep_trac2_ref(&iwk(iwk5), mneed);
                FortranArray1DRef<double> wk_lres(&wk(lres), nsafe);
                bc_embed_ns::bc_embed(ntime1, nbl, lw, lw2, w, mgwk, wk_lres, nsafe,
                                      maxbl, maxgr, lbcemb, iadvance,
                                      idimg, jdimg, kdimg, isav_emb,
                                      ireq_ar_ref, ireq_snd_ref, index_ar_ref,
                                      keep_trac_ref, keep_trac2_ref,
                                      myid, myhost, mycomm, mblk2nd,
                                      nou, bou, nbuf, ibufdim,
                                      iviscg, istat2_em, istat_size, nummem);
            }



            // update 1-1 block boundary conditions
            lres   = 1;
            nsafe  = nwork - lres + 1;
            mneed  = mxbli * 5;
            iwk1   = 1;
            iwk2   = iwk1 + mneed;
            iwk3   = iwk2 + mneed;
            iwk4   = iwk3 + mneed;
            iwk5   = iwk4 + mneed * 2;
            iwk6   = iwk5 + mneed;
            iwork1 = iwork - iwk6;
            if (iwork1 < 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "stopping...not enough integer work space for subroutine bc_blkint");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                int ntime1 = 1;
                FortranArray1DRef<int> ireq_ar_ref(&iwk(iwk1), mneed);
                FortranArray1DRef<int> index_ar_ref(&iwk(iwk2), mneed);
                FortranArray1DRef<int> ireq_snd_ref(&iwk(iwk3), mneed);
                FortranArray2DRef<int> keep_trac_ref(&iwk(iwk4), mneed, 2);
                FortranArray1DRef<int> keep_trac2_ref(&iwk(iwk5), mneed);
                FortranArray1DRef<double> wk_lres(&wk(lres), nsafe);
                bc_blkint_ns::bc_blkint(ntime1, nbl, lw, lw2, w, mgwk, wk_lres, nsafe,
                                        maxbl, maxgr, mxbli, iadvance,
                                        geom_miss, epsc0, nblk, nbli,
                                        limblk, isva, nblon, jdimg, kdimg,
                                        idimg, mblk2nd, isav_blk,
                                        ireq_ar_ref, index_ar_ref,
                                        ireq_snd_ref, keep_trac_ref,
                                        keep_trac2_ref,
                                        nou, bou, nbuf, ibufdim,
                                        myid, myhost, mycomm,
                                        istat2_bl, istat_size, nummem);
            }

            // update patch-grid boundary conditions
            lres   = 1;
            nsafe  = nwork - lres + 1;
            mneed  = intmax * nsub1 * 3;
            iwk1   = 1;
            iwk2   = iwk1 + mneed;
            iwk3   = iwk2 + mneed;
            iwk4   = iwk3 + mneed;
            iwk5   = iwk4 + mneed * 2;
            iwk6   = iwk5 + mneed * 2;
            iwork1 = iwork - iwk6;
            if (iwork1 < 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "stopping...not enough integer work space for subroutine bc_patch");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                int ntime1 = 1;
                FortranArray1DRef<int> ireq_ar_ref(&iwk(iwk1), mneed);
                FortranArray1DRef<int> ireq_snd_ref(&iwk(iwk2), mneed);
                FortranArray1DRef<int> index_ar_ref(&iwk(iwk3), mneed);
                FortranArray3DRef<int> keep_trac_ref(&iwk(iwk4), mneed, 2, 1);
                FortranArray2DRef<int> keep_trac2_ref(&iwk(iwk5), mneed, 1);
                FortranArray1DRef<double> wk_lres(&wk(lres), nsafe);
                bc_patch_ns::bc_patch(ntime1, nbl, lw, lw2, w, mgwk, wk_lres, nsafe,
                                      maxbl, maxgr, intmax, nsub1, maxxe,
                                      iadvance, jdimg, kdimg, idimg,
                                      ninter, windex, iindex, nblkpt,
                                      dthetxx, dthetyy, dthetzz,
                                      isav_pat, isav_pat_b,
                                      ireq_ar_ref, ireq_snd_ref, index_ar_ref,
                                      keep_trac_ref, keep_trac2_ref,
                                      myid, myhost, mycomm, mblk2nd,
                                      nou, bou, nbuf, ibufdim,
                                      istat2_pa, istat_size, nummem);
            }



            // update chimera boundary conditions
            // don't update interior (fringe) points if cell center data
            // is to be output
            int_updt = 1;
            for (n = 1; n <= nplot3d; n++) {
                if (inpl3d(n,2) > 0) {
                    int_updt = 0;
                }
            }

            for (nbl = 1; nbl <= nblock; nbl++) {
                if (iadvance(nbl) >= 0) {
                    if (level == levelg(nbl)) {
                        if (mblk2nd(nbl) == myid) {
                            lead_ns::lead(nbl, lw, lw2, maxbl);
                            lres  = 1;
                            nsafe = nwork - lres + 1;
                            int ntime1 = 1;
                            FortranArray1DRef<double> wk_lres(&wk(lres), nsafe);
                            bc_xmera_ns::bc_xmera(ntime1, nbl, lw, lw2, w, mgwk,
                                                  wk_lres, nsafe, maxbl, iitot,
                                                  iviscg, iovrlp, lbg, ibpntsg,
                                                  qb, iibg, kkbg, jjbg, ibcg,
                                                  nou, bou, nbuf, ibufdim,
                                                  int_updt, nummem);
                        }
                    }
                }
            }

            // call qface to install face-center values in qi0/qj0/qk0 arrays
            for (nbl = 1; nbl <= nblock; nbl++) {
                if (level == levelg(nbl)) {
                    if (mblk2nd(nbl) == myid) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        int ldim = 5;
                        {
                            FortranArray4DRef<double> q_ref(&w(lq), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> qj0_ref(&w(lqj0), kdim, idim-1, 5, 4);
                            FortranArray4DRef<double> qk0_ref(&w(lqk0), jdim, idim-1, 5, 4);
                            FortranArray4DRef<double> qi0_ref(&w(lqi0), jdim, kdim, 5, 4);
                            FortranArray3DRef<double> bcj_ref(&w(lbcj), kdim, idim-1, 2);
                            FortranArray3DRef<double> bck_ref(&w(lbck), jdim, idim-1, 2);
                            FortranArray3DRef<double> bci_ref(&w(lbci), jdim, kdim, 2);
                            FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                            qface_ns::qface(jdim, kdim, idim, q_ref, qj0_ref, qk0_ref,
                                            qi0_ref, bcj_ref, bck_ref, bci_ref,
                                            blank_ref, ldim);
                        }
                        if (ivmx >= 2) {
                            ldim = 1;
                            FortranArray4DRef<double> vis_ref(&w(lvis), jdim, kdim, idim, 1);
                            FortranArray4DRef<double> vj0_ref(&w(lvj0), kdim, idim-1, 1, 4);
                            FortranArray4DRef<double> vk0_ref(&w(lvk0), jdim, idim-1, 1, 4);
                            FortranArray4DRef<double> vi0_ref(&w(lvi0), jdim, kdim, 1, 4);
                            FortranArray3DRef<double> bcj_ref(&w(lbcj), kdim, idim-1, 2);
                            FortranArray3DRef<double> bck_ref(&w(lbck), jdim, idim-1, 2);
                            FortranArray3DRef<double> bci_ref(&w(lbci), jdim, kdim, 2);
                            FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                            qface_ns::qface(jdim, kdim, idim, vis_ref, vj0_ref, vk0_ref,
                                            vi0_ref, bcj_ref, bck_ref, bci_ref,
                                            blank_ref, ldim);
                        }
                    }
                }
            }

        } // end do level=lglobal,levelt(iseq)

        level = level_sav;

    } // end if (nplot3d.gt.0 .or. nprint.gt.0 .or. iteravg.gt.0)



    // ---- mod for averaging at grid points --------------------------------
    if (iteravg > 0) {
        // prepare to call plot3davg and then call it
        for (n = 1; n <= nblock; n++) {
            nbl = n;

            if (level != levelg(nbl)) goto label999; // skip block if not global
            if (mblk2nd(nbl) == myid || myid == myhost) {

                lead_ns::lead(nbl, lw, lw2, maxbl);

                i1 = 1;
                i2 = idim;
                i3 = 1;
                j1 = 1;
                j2 = jdim;
                j3 = 1;
                k1 = 1;
                k2 = kdim;
                k3 = 1;

                jdw = (j2 - j1) / j3 + 1;
                kdw = (k2 - k1) / k3 + 1;
                idw = (i2 - i1) / i3 + 1;

                // check storage availability
                ixwk = 1;
                nset = 5;
                ibwk = ixwk + jdw * kdw * idw * nset;
                ixgk = ibwk + jdim * kdim * idim * 2;
                ixvk = ixgk + jdw * kdw * idw * 4;
                nroom = nwork - (ixvk + jdw * kdw * idw * 5 - 1);
                if (nroom < 0) {
                    if (myid == myhost) {
                        fortran_write_unit(11, " not enough memory for plot3davg\n");
                        {
                            char buf[120];
                            std::snprintf(buf, 120, " have, need = %12d%12d\n",
                                          nwork, nwork - nroom);
                            fortran_write_unit(11, "%s", buf);
                        }
                        fortran_write_unit(11, " Aborting - not calling plot3davg\n");
                        std::exit(0);
                    }
                } else {
                    int iflag_avg = 1;
                    int iover_avg = iovrlp(nbl);
                    FortranArray4DRef<double> q_ref(&w(lq), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qi0_ref(&w(lqi0), jdim, kdim, 5, 4);
                    FortranArray4DRef<double> qj0_ref(&w(lqj0), kdim, idim-1, 5, 4);
                    FortranArray4DRef<double> qk0_ref(&w(lqk0), jdim, idim-1, 5, 4);
                    FortranArray3DRef<double> x_ref(&w(lx), jdim, kdim, idim);
                    FortranArray3DRef<double> y_ref(&w(ly), jdim, kdim, idim);
                    FortranArray3DRef<double> z_ref(&w(lz), jdim, kdim, idim);
                    FortranArray4DRef<double> xw_ref(&wk(ixwk), jdw, kdw, idw, nset);
                    FortranArray4DRef<double> blank2_ref(&wk(ibwk), jdim, kdim, idim, 2);
                    FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                    FortranArray4DRef<double> xg_ref(&wk(ixgk), jdw, kdw, idw, 4);
                    FortranArray3DRef<double> vis_ref(&w(lvis), jdim, kdim, idim);
                    FortranArray3DRef<double> bcj_ref(&w(lbcj), kdim, idim-1, 2);
                    FortranArray3DRef<double> bck_ref(&w(lbck), jdim, idim-1, 2);
                    FortranArray3DRef<double> bci_ref(&w(lbci), jdim, kdim, 2);
                    FortranArray4DRef<double> vj0_ref(&w(lvj0), kdim, idim-1, 1, 4);
                    FortranArray4DRef<double> vk0_ref(&w(lvk0), jdim, idim-1, 1, 4);
                    FortranArray4DRef<double> vi0_ref(&w(lvi0), jdim, kdim, 1, 4);
                    FortranArray4DRef<double> xv_ref(&wk(ixvk), jdw, kdw, idw, 5);
                    FortranArray4DRef<double> sj_ref(&w(lsj), jdim, kdim, idim-1, 5);
                    FortranArray4DRef<double> sk_ref(&w(lsk), jdim, kdim, idim-1, 5);
                    FortranArray4DRef<double> si_ref(&w(lsi), jdim, kdim, idim, 5);
                    FortranArray3DRef<double> vol_ref(&w(lvol), jdim, kdim, idim-1);
                    FortranArray4DRef<double> qavg_ref(&w(lqavg), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> q2avg_ref(&w(lq2avg), jdim, kdim, idim, 5);
                    int ifunc_avg = 0;
                    plot3davg_ns::plot3davg(jdim, kdim, idim,
                                            i1, i2, i3, j1, j2, j3, k1, k2, k3,
                                            q_ref, qi0_ref, qj0_ref, qk0_ref,
                                            x_ref, y_ref, z_ref,
                                            xw_ref, blank2_ref, blank_ref, xg_ref,
                                            iflag_avg, vis_ref, iover_avg, nbl, nmap,
                                            bcj_ref, bck_ref, bci_ref,
                                            vj0_ref, vk0_ref, vi0_ref,
                                            ifunc_avg, n, jdw, kdw, idw,
                                            nplots, jdimg, kdimg, idimg, nblcg,
                                            jsg, ksg, isg, jeg, keg, ieg,
                                            ninter, iindex, intmax, nsub1,
                                            maxxe, nblk, nbli, limblk, isva,
                                            nblon, mxbli, thetay, maxbl, maxgr,
                                            myid, myhost, mycomm, mblk2nd,
                                            inpl3d, nblock, nblkpt,
                                            xv_ref, sj_ref, sk_ref, si_ref,
                                            vol_ref, nset,
                                            qavg_ref, q2avg_ref, nt, movabs);
                }

            } // end if mblk2nd(nbl).eq.myid or myid.eq.myhost
            label999:;
        } // end do n=1,nblock
    } // end if iteravg.gt.0

    if (movabs == 0) {
        return; // wouldn't have called qout here if not for averaging
    } else if (nt != nt / movabs * movabs) {
        // no need to continue with qout for this time step
        return;
    }
    // End of changes to qout for averaging; remaining portion same as before



    if (nplot3d <= 0) goto label231;

    ncount   = 0;
    np3d     = nplot3d;
    ifunc    = 0;
    ifuncuse = 0;

    // if zone has function file output, all must, and all
    // must have the same function output
    if (std::abs(inpl3d(1,2)) > 2) {
        ifunc    = inpl3d(1,2);
        ifuncuse = 1;
    }
    if (std::abs(inpl3d(1,2)) == 2) {
        ifuncuse = ifunct;
    }

    for (n = 1; n <= nplot3d; n++) {
        if (n == 1) {
            if (myid == myhost) {
                if (ibin == 0) {
                    if (icall1 == 0) {
                        char buf[32];
                        std::snprintf(buf, 32, "%5d%5d%5d", np3d, 0, 0);
                        fortran_write_unit(3, "%5d\n", np3d);
                    }
                    fortran_write_unit(4, "%5d\n", np3d);
                } else {
                    if (icall1 == 0) {
                        FILE* f3 = fortran_get_unit(3);
                        fwrite(&np3d, sizeof(int), 1, f3);
                    }
                    FILE* f4 = fortran_get_unit(4);
                    fwrite(&np3d, sizeof(int), 1, f4);
                }
            }
        }

        nbl = inpl3d(n,1);

        if (nbl > nblock) {
            if (myid == myhost) {
                fortran_write_unit(11, " Block%3d does not exist.  No plot3d output printed.\n", nbl);
            }
            continue; // go to 60
        }

        ncount++;
        i1 = inpl3d(n,3);
        i2 = inpl3d(n,4);
        i3 = inpl3d(n,5);
        j1 = inpl3d(n,6);
        j2 = inpl3d(n,7);
        j3 = inpl3d(n,8);
        k1 = inpl3d(n,9);
        k2 = inpl3d(n,10);
        k3 = inpl3d(n,11);
        if (inpl3d(n,2) > 0) {
            // cell center dimensions
            lead_ns::lead(nbl, lw, lw2, maxbl);
            i2 = std::min(idim-1, i2);
            j2 = std::min(jdim-1, j2);
            k2 = std::min(kdim-1, k2);
            i1 = std::min(idim-1, i1);
            j1 = std::min(jdim-1, j1);
            k1 = std::min(kdim-1, k1);
        }
        ip3ddim(1,ncount) = (i2 - i1) / i3 + 1;
        ip3ddim(2,ncount) = (j2 - j1) / j3 + 1;
        ip3ddim(3,ncount) = (k2 - k1) / k3 + 1;
    } // end do 60

    if (myid == myhost) {
        if (ibin == 0) {
            if (i2d == 0) {
                if (icall1 == 0) {
                    FILE* f3 = fortran_get_unit(3);
                    for (int nn = 1; nn <= ncount; nn++) {
                        fprintf(f3, "%5d%5d%5d", ip3ddim(1,nn), ip3ddim(2,nn), ip3ddim(3,nn));
                    }
                    fprintf(f3, "\n");
                }
                if (ifuncuse == 0) {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        fprintf(f4, "%5d%5d%5d", ip3ddim(1,nn), ip3ddim(2,nn), ip3ddim(3,nn));
                    }
                    fprintf(f4, "\n");
                } else {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        fprintf(f4, "%5d%5d%5d%5d", ip3ddim(1,nn), ip3ddim(2,nn), ip3ddim(3,nn), ifuncuse);
                    }
                    fprintf(f4, "\n");
                }
            } else {
                // i2d != 0
                if (icall1 == 0) {
                    FILE* f3 = fortran_get_unit(3);
                    for (int nn = 1; nn <= ncount; nn++) {
                        fprintf(f3, "%5d%5d%5d", ip3ddim(2,nn), ip3ddim(3,nn), 0);
                    }
                    fprintf(f3, "\n");
                }
                if (ifuncuse == 0) {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        fprintf(f4, "%5d%5d%5d", ip3ddim(2,nn), ip3ddim(3,nn), 0);
                    }
                    fprintf(f4, "\n");
                } else {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        fprintf(f4, "%5d%5d%5d", ip3ddim(2,nn), ip3ddim(3,nn), ifuncuse);
                    }
                    fprintf(f4, "\n");
                }
            }
        } else {
            // binary
            if (i2d == 0) {
                if (icall1 == 0) {
                    FILE* f3 = fortran_get_unit(3);
                    for (int nn = 1; nn <= ncount; nn++) {
                        int v1 = ip3ddim(1,nn), v2 = ip3ddim(2,nn), v3 = ip3ddim(3,nn);
                        fwrite(&v1, sizeof(int), 1, f3);
                        fwrite(&v2, sizeof(int), 1, f3);
                        fwrite(&v3, sizeof(int), 1, f3);
                    }
                }
                if (ifuncuse == 0) {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        int v1 = ip3ddim(1,nn), v2 = ip3ddim(2,nn), v3 = ip3ddim(3,nn);
                        fwrite(&v1, sizeof(int), 1, f4);
                        fwrite(&v2, sizeof(int), 1, f4);
                        fwrite(&v3, sizeof(int), 1, f4);
                    }
                } else {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        int v1 = ip3ddim(1,nn), v2 = ip3ddim(2,nn), v3 = ip3ddim(3,nn);
                        fwrite(&v1, sizeof(int), 1, f4);
                        fwrite(&v2, sizeof(int), 1, f4);
                        fwrite(&v3, sizeof(int), 1, f4);
                        fwrite(&ifuncuse, sizeof(int), 1, f4);
                    }
                }
            } else {
                // i2d != 0, binary
                if (icall1 == 0) {
                    FILE* f3 = fortran_get_unit(3);
                    for (int nn = 1; nn <= ncount; nn++) {
                        int v2 = ip3ddim(2,nn), v3 = ip3ddim(3,nn);
                        fwrite(&v2, sizeof(int), 1, f3);
                        fwrite(&v3, sizeof(int), 1, f3);
                    }
                }
                if (ifuncuse == 0) {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        int v2 = ip3ddim(2,nn), v3 = ip3ddim(3,nn);
                        fwrite(&v2, sizeof(int), 1, f4);
                        fwrite(&v3, sizeof(int), 1, f4);
                    }
                } else {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        int v2 = ip3ddim(2,nn), v3 = ip3ddim(3,nn);
                        fwrite(&v2, sizeof(int), 1, f4);
                        fwrite(&v3, sizeof(int), 1, f4);
                        fwrite(&ifuncuse, sizeof(int), 1, f4);
                    }
                }
            }
        }
    }



    // correspondence between global block number and plot3d block number
    for (n = 1; n <= nblock; n++) {
        nmap(n) = 1;
        for (nnn = 1; nnn <= nplot3d; nnn++) {
            m = inpl3d(nnn,1);
            if (n == m) nmap(n) = nnn;
        }
    }

    if (myid == myhost) {
        if (lhdr > 0) fortran_write_unit(11, "\n");
    }

    for (n = 1; n <= nplot3d; n++) {

        nbl = inpl3d(n,1);

        if (nbl > nblock || nbl <= 0) continue; // go to 70

        if (iblnkfr == 0) {
            // temporarily set blank values at fringe points (not holes)
            // to 1 for plotting purposes
            lead_ns::lead(nbl, lw, lw2, maxbl);
            if (mblk2nd(nbl) == myid && iovrlp(nbl) != 0) {
                double blnkval = 1.0;
                FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                blnkfr_ns::blnkfr(nbl, iibg, kkbg, jjbg, ibpntsg, lbg, iitot,
                                   blank_ref, jdim, kdim, idim, maxbl, blnkval);
            }
        }

        i1 = inpl3d(n,3);
        i2 = inpl3d(n,4);
        i3 = inpl3d(n,5);
        j1 = inpl3d(n,6);
        j2 = inpl3d(n,7);
        j3 = inpl3d(n,8);
        k1 = inpl3d(n,9);
        k2 = inpl3d(n,10);
        k3 = inpl3d(n,11);

        iflag = 1;

        if (inpl3d(n,2) <= 0) {
            // grid point data
            if (mblk2nd(nbl) == myid || myid == myhost) {
                lead_ns::lead(nbl, lw, lw2, maxbl);

                jdw = (j2 - j1) / j3 + 1;
                kdw = (k2 - k1) / k3 + 1;
                idw = (i2 - i1) / i3 + 1;

                ixwk = 1;
                nset = 5;
                ibwk = ixwk + jdw * kdw * idw * nset;
                ixgk = ibwk + jdim * kdim * idim * 2;
                ixvk = ixgk + jdw * kdw * idw * 4;
                nroom = nwork - (ixvk + jdw * kdw * idw * 5 - 1);
                if (nroom < 0) {
                    if (myid == myhost) {
                        fortran_write_unit(11, " not enough memory for plot3d\n");
                        {
                            char buf[120];
                            std::snprintf(buf, 120, " have, need = %12d%12d\n",
                                          nwork, nwork - nroom);
                            fortran_write_unit(11, "%s", buf);
                        }
                        fortran_write_unit(11, " not writing out plot3d files\n");
                    }
                } else {
                    int iover_v = iovrlp(nbl);
                    FortranArray4DRef<double> q_ref(&w(lq), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qi0_ref(&w(lqi0), jdim, kdim, 5, 4);
                    FortranArray4DRef<double> qj0_ref(&w(lqj0), kdim, idim-1, 5, 4);
                    FortranArray4DRef<double> qk0_ref(&w(lqk0), jdim, idim-1, 5, 4);
                    FortranArray3DRef<double> x_ref(&w(lx), jdim, kdim, idim);
                    FortranArray3DRef<double> y_ref(&w(ly), jdim, kdim, idim);
                    FortranArray3DRef<double> z_ref(&w(lz), jdim, kdim, idim);
                    FortranArray4DRef<double> xw_ref(&wk(ixwk), jdw, kdw, idw, nset);
                    FortranArray4DRef<double> blank2_ref(&wk(ibwk), jdim, kdim, idim, 2);
                    FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                    FortranArray4DRef<double> xg_ref(&wk(ixgk), jdw, kdw, idw, 4);
                    FortranArray3DRef<double> vis_ref(&w(lvis), jdim, kdim, idim);
                    FortranArray3DRef<double> bcj_ref(&w(lbcj), kdim, idim-1, 2);
                    FortranArray3DRef<double> bck_ref(&w(lbck), jdim, idim-1, 2);
                    FortranArray3DRef<double> bci_ref(&w(lbci), jdim, kdim, 2);
                    FortranArray4DRef<double> vj0_ref(&w(lvj0), kdim, idim-1, 1, 4);
                    FortranArray4DRef<double> vk0_ref(&w(lvk0), jdim, idim-1, 1, 4);
                    FortranArray4DRef<double> vi0_ref(&w(lvi0), jdim, kdim, 1, 4);
                    FortranArray4DRef<double> xv_ref(&wk(ixvk), jdw, kdw, idw, 5);
                    FortranArray4DRef<double> sj_ref(&w(lsj), jdim, kdim, idim-1, 5);
                    FortranArray4DRef<double> sk_ref(&w(lsk), jdim, kdim, idim-1, 5);
                    FortranArray4DRef<double> si_ref(&w(lsi), jdim, kdim, idim, 5);
                    FortranArray3DRef<double> vol_ref(&w(lvol), jdim, kdim, idim-1);
                    plot3d_ns::plot3d(jdim, kdim, idim, i1, i2, i3, j1, j2, j3,
                                      k1, k2, k3, q_ref, qi0_ref, qj0_ref, qk0_ref,
                                      x_ref, y_ref, z_ref, xw_ref, blank2_ref,
                                      blank_ref, xg_ref, iflag, vis_ref,
                                      iover_v, nbl, nmap, bcj_ref, bck_ref, bci_ref,
                                      vj0_ref, vk0_ref, vi0_ref, ifunc, n,
                                      jdw, kdw, idw, nplots, jdimg, kdimg, idimg,
                                      nblcg, jsg, ksg, isg, jeg, keg, ieg,
                                      ninter, iindex, intmax, nsub1, maxxe,
                                      nblk, nbli, limblk, isva, nblon, mxbli,
                                      thetay, maxbl, maxgr, myid, myhost, mycomm,
                                      mblk2nd, inpl3d, nblock, nblkpt,
                                      xv_ref, sj_ref, sk_ref, si_ref, vol_ref, nset);
                }
            }

        } else if (inpl3d(n,2) == 1 || inpl3d(n,2) > 2) {


            // cell center or face center data
            if (mblk2nd(nbl) == myid || myid == myhost) {
                lead_ns::lead(nbl, lw, lw2, maxbl);

                i2 = std::min(idim-1, i2);
                j2 = std::min(jdim-1, j2);
                k2 = std::min(kdim-1, k2);
                i1 = std::min(idim-1, i1);
                j1 = std::min(jdim-1, j1);
                k1 = std::min(kdim-1, k1);
                jdw = (j2 - j1) / j3 + 1;
                kdw = (k2 - k1) / k3 + 1;
                idw = (i2 - i1) / i3 + 1;

                ixwk = 1;
                ibwk = ixwk + jdw * kdw * idw * 5;
                ixgk = ibwk + jdim * kdim * idim;
                nroom = nwork - (ixgk + jdw * kdw * idw * 4 - 1);
                if (nroom < 0) {
                    if (myid == myhost) {
                        fortran_write_unit(11, " not enough memory for plot3c\n");
                        {
                            char buf[120];
                            std::snprintf(buf, 120, " have, need = %12d%12d\n",
                                          nwork, nwork - nroom);
                            fortran_write_unit(11, "%s", buf);
                        }
                        fortran_write_unit(11, " not writing out plot3d files\n");
                    }
                } else {
                    int iover_v = iovrlp(nbl);
                    FortranArray4DRef<double> q_ref(&w(lq), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qi0_ref(&w(lqi0), jdim, kdim, 5, 4);
                    FortranArray4DRef<double> qj0_ref(&w(lqj0), kdim, idim-1, 5, 4);
                    FortranArray4DRef<double> qk0_ref(&w(lqk0), jdim, idim-1, 5, 4);
                    FortranArray3DRef<double> x_ref(&w(lx), jdim, kdim, idim);
                    FortranArray3DRef<double> y_ref(&w(ly), jdim, kdim, idim);
                    FortranArray3DRef<double> z_ref(&w(lz), jdim, kdim, idim);
                    FortranArray4DRef<double> xw_ref(&wk(ixwk), jdw, kdw, idw, 5);
                    FortranArray3DRef<double> blank2_ref(&wk(ibwk), jdim, kdim, idim);
                    FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                    FortranArray4DRef<double> xg_ref(&wk(ixgk), jdw, kdw, idw, 4);
                    FortranArray3DRef<double> vis_ref(&w(lvis), jdim, kdim, idim);
                    FortranArray4DRef<double> vi0_ref(&w(lvi0), jdim, kdim, 1, 4);
                    FortranArray4DRef<double> vj0_ref(&w(lvj0), kdim, idim-1, 1, 4);
                    FortranArray4DRef<double> vk0_ref(&w(lvk0), jdim, idim-1, 1, 4);
                    FortranArray3DRef<double> smin_ref(&w(lsnk0), jdim, kdim, idim);
                    plot3c_ns::plot3c(jdim, kdim, idim, i1, i2, i3, j1, j2, j3,
                                      k1, k2, k3, q_ref, qi0_ref, qj0_ref, qk0_ref,
                                      x_ref, y_ref, z_ref, xw_ref, blank2_ref,
                                      blank_ref, xg_ref, iflag, vis_ref,
                                      vi0_ref, vj0_ref, vk0_ref,
                                      iover_v, nbl, nmap, smin_ref, ifunc, n,
                                      jdw, kdw, idw, nplots, jdimg, kdimg, idimg,
                                      nblcg, jsg, ksg, isg, jeg, keg, ieg,
                                      ninter, iindex, intmax, nsub1, maxxe,
                                      nblk, nbli, limblk, isva, nblon, mxbli,
                                      thetay, maxbl, maxgr, myid, myhost, mycomm,
                                      mblk2nd, inpl3d, nblock, nblkpt, ip3dsurf);
                }
            }

        } else {
            // cell center turbulence data (plot3d q file format)
            if (mblk2nd(nbl) == myid || myid == myhost) {
                if (ivmx > 1) {
                    lead_ns::lead(nbl, lw, lw2, maxbl);

                    i2 = std::min(idim-1, i2);
                    j2 = std::min(jdim-1, j2);
                    k2 = std::min(kdim-1, k2);
                    i1 = std::min(idim-1, i1);
                    j1 = std::min(jdim-1, j1);
                    k1 = std::min(kdim-1, k1);
                    jdw = (j2 - j1) / j3 + 1;
                    kdw = (k2 - k1) / k3 + 1;
                    idw = (i2 - i1) / i3 + 1;

                    if (ivmx == 8 || ivmx == 9 || ivmx >= 11 || lowmem_ux == 0) {
                        // have permanent storage for ux, starting at lux
                        ixwk     = 1;
                        ifuncdim = std::max(5, ifunct);
                        ibwk     = ixwk + jdw * kdw * idw * ifuncdim;
                        ixgk     = ibwk + jdim * kdim * idim;
                        ibwk3    = ixgk + jdw * kdw * idw * 4;
                        nroom    = nwork - (ibwk3 + jdim * kdim * 9);
                        if (nroom < 0) {
                            if (myid == myhost) {
                                fortran_write_unit(11, " not enough memory for plot3t\n");
                                {
                                    char buf[120];
                                    std::snprintf(buf, 120, " have, need = %12d%12d\n",
                                                  nwork, nwork - nroom);
                                    fortran_write_unit(11, "%s", buf);
                                }
                                fortran_write_unit(11, " not writing out plot3d files\n");
                            }
                        } else {
                            int iover_v = iovrlp(nbl);
                            FortranArray4DRef<double> q_ref(&w(lq), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> x_ref(&w(lx), jdim, kdim, idim);
                            FortranArray3DRef<double> y_ref(&w(ly), jdim, kdim, idim);
                            FortranArray3DRef<double> z_ref(&w(lz), jdim, kdim, idim);
                            FortranArray4DRef<double> xw_ref(&wk(ixwk), jdw, kdw, idw, ifuncdim);
                            FortranArray3DRef<double> blank2_ref(&wk(ibwk), jdim, kdim, idim);
                            FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                            FortranArray4DRef<double> xg_ref(&wk(ixgk), jdw, kdw, idw, 4);
                            FortranArray3DRef<double> vis_ref(&w(lvis), jdim, kdim, idim);
                            FortranArray4DRef<double> sj_ref(&w(lsj), jdim, kdim, idim-1, 5);
                            FortranArray4DRef<double> sk_ref(&w(lsk), jdim, kdim, idim-1, 5);
                            FortranArray4DRef<double> si_ref(&w(lsi), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> smin_ref(&w(lsnk0), jdim, kdim, idim);
                            FortranArray4DRef<double> ux_ref(&w(lux), jdim, kdim, idim, 9);
                            FortranArray4DRef<double> turre_ref(&w(lxib), jdim, kdim, idim, nummem);
                            FortranArray3DRef<double> vol_ref(&w(lvol), jdim, kdim, idim-1);
                            FortranArray4DRef<double> qj0_ref(&w(lqj0), kdim, idim-1, 5, 4);
                            FortranArray4DRef<double> qk0_ref(&w(lqk0), jdim, idim-1, 5, 4);
                            FortranArray4DRef<double> qi0_ref(&w(lqi0), jdim, kdim, 5, 4);
                            FortranArray3DRef<double> bcj_ref(&w(lbcj), kdim, idim-1, 2);
                            FortranArray3DRef<double> bck_ref(&w(lbck), jdim, idim-1, 2);
                            FortranArray3DRef<double> bci_ref(&w(lbci), jdim, kdim, 2);
                            FortranArray3DRef<double> wt_ref(&wk(ibwk3), jdim, kdim, 9);
                            FortranArray3DRef<double> cmuv_ref(&w(lcmuv), jdim, kdim, idim);
                            FortranArray3DRef<double> volj0_ref(&w(lvolj0), kdim, idim-1, 1);
                            FortranArray3DRef<double> volk0_ref(&w(lvolk0), jdim, idim-1, 1);
                            FortranArray3DRef<double> voli0_ref(&w(lvoli0), jdim, kdim, 1);
                            plot3t_ns::plot3t(jdim, kdim, idim, i1, i2, i3, j1, j2, j3,
                                              k1, k2, k3, q_ref, x_ref, y_ref, z_ref,
                                              xw_ref, blank2_ref, blank_ref, xg_ref,
                                              vis_ref, iover_v, nbl, nmap,
                                              sj_ref, sk_ref, si_ref, smin_ref,
                                              ux_ref, turre_ref, vol_ref,
                                              qj0_ref, qk0_ref, qi0_ref,
                                              bcj_ref, bck_ref, bci_ref,
                                              wt_ref, cmuv_ref,
                                              jdw, kdw, idw, nplots,
                                              jdimg, kdimg, idimg, nblcg,
                                              jsg, ksg, isg, jeg, keg, ieg,
                                              ninter, iindex, intmax, nsub1,
                                              maxxe, nblk, nbli, limblk, isva,
                                              nblon, mxbli, thetay, maxbl, maxgr,
                                              myid, myhost, mycomm, mblk2nd,
                                              inpl3d, nblock, nblkpt,
                                              volj0_ref, volk0_ref, voli0_ref,
                                              vormax, ivmax, jvmax, kvmax,
                                              nummem, ifuncdim);
                        }
                    } else {


                        // need to grab ux storage from temporary work array
                        ixwk     = 1;
                        ifuncdim = std::max(5, ifunct);
                        ibwk     = ixwk + jdw * kdw * idw * ifuncdim;
                        ixgk     = ibwk + jdim * kdim * idim;
                        ibwk2    = ixgk + jdw * kdw * idw * 4;
                        ibwk3    = ibwk2 + (jdim-1) * (kdim-1) * (idim-1) * 9;
                        nroom    = nwork - (ibwk3 + jdim * kdim * 9);
                        if (nroom < 0) {
                            if (myid == myhost) {
                                fortran_write_unit(11, " not enough memory for plot3t\n");
                                {
                                    char buf[120];
                                    std::snprintf(buf, 120, " have, need = %12d%12d\n",
                                                  nwork, nwork - nroom);
                                    fortran_write_unit(11, "%s", buf);
                                }
                                fortran_write_unit(11, " not writing out plot3d files\n");
                            }
                        } else {
                            int iover_v = iovrlp(nbl);
                            FortranArray4DRef<double> q_ref(&w(lq), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> x_ref(&w(lx), jdim, kdim, idim);
                            FortranArray3DRef<double> y_ref(&w(ly), jdim, kdim, idim);
                            FortranArray3DRef<double> z_ref(&w(lz), jdim, kdim, idim);
                            FortranArray4DRef<double> xw_ref(&wk(ixwk), jdw, kdw, idw, ifuncdim);
                            FortranArray3DRef<double> blank2_ref(&wk(ibwk), jdim, kdim, idim);
                            FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                            FortranArray4DRef<double> xg_ref(&wk(ixgk), jdw, kdw, idw, 4);
                            FortranArray3DRef<double> vis_ref(&w(lvis), jdim, kdim, idim);
                            FortranArray4DRef<double> sj_ref(&w(lsj), jdim, kdim, idim-1, 5);
                            FortranArray4DRef<double> sk_ref(&w(lsk), jdim, kdim, idim-1, 5);
                            FortranArray4DRef<double> si_ref(&w(lsi), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> smin_ref(&w(lsnk0), jdim, kdim, idim);
                            FortranArray4DRef<double> ux_ref(&wk(ibwk2), (jdim-1), (kdim-1), (idim-1), 9);
                            FortranArray4DRef<double> turre_ref(&w(lxib), jdim, kdim, idim, nummem);
                            FortranArray3DRef<double> vol_ref(&w(lvol), jdim, kdim, idim-1);
                            FortranArray4DRef<double> qj0_ref(&w(lqj0), kdim, idim-1, 5, 4);
                            FortranArray4DRef<double> qk0_ref(&w(lqk0), jdim, idim-1, 5, 4);
                            FortranArray4DRef<double> qi0_ref(&w(lqi0), jdim, kdim, 5, 4);
                            FortranArray3DRef<double> bcj_ref(&w(lbcj), kdim, idim-1, 2);
                            FortranArray3DRef<double> bck_ref(&w(lbck), jdim, idim-1, 2);
                            FortranArray3DRef<double> bci_ref(&w(lbci), jdim, kdim, 2);
                            FortranArray3DRef<double> wt_ref(&wk(ibwk3), jdim, kdim, 9);
                            FortranArray3DRef<double> cmuv_ref(&w(lcmuv), jdim, kdim, idim);
                            FortranArray3DRef<double> volj0_ref(&w(lvolj0), kdim, idim-1, 1);
                            FortranArray3DRef<double> volk0_ref(&w(lvolk0), jdim, idim-1, 1);
                            FortranArray3DRef<double> voli0_ref(&w(lvoli0), jdim, kdim, 1);
                            plot3t_ns::plot3t(jdim, kdim, idim, i1, i2, i3, j1, j2, j3,
                                              k1, k2, k3, q_ref, x_ref, y_ref, z_ref,
                                              xw_ref, blank2_ref, blank_ref, xg_ref,
                                              vis_ref, iover_v, nbl, nmap,
                                              sj_ref, sk_ref, si_ref, smin_ref,
                                              ux_ref, turre_ref, vol_ref,
                                              qj0_ref, qk0_ref, qi0_ref,
                                              bcj_ref, bck_ref, bci_ref,
                                              wt_ref, cmuv_ref,
                                              jdw, kdw, idw, nplots,
                                              jdimg, kdimg, idimg, nblcg,
                                              jsg, ksg, isg, jeg, keg, ieg,
                                              ninter, iindex, intmax, nsub1,
                                              maxxe, nblk, nbli, limblk, isva,
                                              nblon, mxbli, thetay, maxbl, maxgr,
                                              myid, myhost, mycomm, mblk2nd,
                                              inpl3d, nblock, nblkpt,
                                              volj0_ref, volk0_ref, voli0_ref,
                                              vormax, ivmax, jvmax, kvmax,
                                              nummem, ifuncdim);
                        }
                    } // end if ivmx==8 || ...
                } // end if ivmx > 1
            } // end if mblk2nd(nbl)==myid || myid==myhost
        } // end if/else on inpl3d(n,2)

        if (iblnkfr == 0) {
            // reset blank values at fringe points to 0.
            if (mblk2nd(nbl) == myid && iovrlp(nbl) != 0) {
                double blnkval = 0.0;
                FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                blnkfr_ns::blnkfr(nbl, iibg, kkbg, jjbg, ibpntsg, lbg, iitot,
                                   blank_ref, jdim, kdim, idim, maxbl, blnkval);
            }
        }

    } // end do 70 n=1,nplot3d

    if (myid == myhost) {
        if (lhdr > 0) fortran_write_unit(11, "\n");
    }

    // for stationary grid cases, set icall1 flag to prevent output
    // of grid to plot3d file on subsequent calls
    if (iunst == 0) {
        icall1 = 1;
    } else {
        icall1 = 0;
    }

    label231:;



    // ***********************************************************************
    // Print solution data.
    // ***********************************************************************
    if (nprint > 0) {

        if (myid == myhost) {
            if (lhdr > 0 && nplot3d <= 0) fortran_write_unit(11, "\n");
        }

        if (myid == myhost) {
            // write(17,111)(real(title(i)),i=1,20)
            // 111 format(/,20a4)
            {
                FILE* f17 = fortran_get_unit(17);
                fprintf(f17, "\n");
                for (int i = 1; i <= 20; i++) {
                    float tv = title[i-1];
                    char* cp = reinterpret_cast<char*>(&tv);
                    fprintf(f17, "%c%c%c%c", cp[0], cp[1], cp[2], cp[3]);
                }
                fprintf(f17, "\n");
            }
            // write(17,21)
            // 21 format(6x,4hMach,5x,5halpha,6x,4hbeta,6x,4hReUe,3x,7hTinf,dR,6x,4htime)
            fortran_write_unit(17, "      Mach     alpha      beta      ReUe   Tinf,dR      time\n");
            alphaw = (float)cmn_conversion.radtodeg * alpha;
            betaw  = (float)cmn_conversion.radtodeg * beta;
            // write(17,20)real(xmach),real(alphaw),real(betaw),real(reue),real(tinf),real(time)
            // 20 format(3f10.5,e10.3,2f10.5)
            {
                FILE* f17 = fortran_get_unit(17);
                fprintf(f17, "%10.5f%10.5f%10.5f%10.3e%10.5f%10.5f\n",
                        (float)xmach, (float)alphaw, (float)betaw,
                        (float)reue, (float)tinf, (float)time);
            }
        }

        for (n = 1; n <= nprint; n++) {

            nbl = inpr(n,1);

            if (nbl > nblock) {
                if (myid == myhost) {
                    // 77 format(6h Block,i3,36h does not exist.  No output printed.)
                    fortran_write_unit(11, " Block%3d does not exist.  No output printed.\n", nbl);
                }
                continue; // go to 80
            }

            i1 = inpr(n,3);
            i2 = inpr(n,4);
            i3 = inpr(n,5);
            j1 = inpr(n,6);
            j2 = inpr(n,7);
            j3 = inpr(n,8);
            k1 = inpr(n,9);
            k2 = inpr(n,10);
            k3 = inpr(n,11);

            iflag = 2;

            if (inpr(n,2) == 0) {
                // grid point data
                if (mblk2nd(nbl) == myid || myid == myhost) {
                    lead_ns::lead(nbl, lw, lw2, maxbl);

                    if (myid == myhost) {
                        // 348 format(//1x,5hBLOCK,i4,2x,5h(GRID,i4,1h),5x,15hIDIM,JDIM,KDIM=,3I5)
                        FILE* f17 = fortran_get_unit(17);
                        fprintf(f17, "\n\n BLOCK%4d  (GRID%4d)     IDIM,JDIM,KDIM=%5d%5d%5d\n",
                                nbl, igridg(nbl), idim, jdim, kdim);
                        // 349 format(1x,32hNOTE: endpts may not be reliable)
                        fortran_write_unit(17, " NOTE: endpts may not be reliable\n");
                    }

                    jdw = (j2 - j1) / j3 + 1;
                    kdw = (k2 - k1) / k3 + 1;
                    idw = (i2 - i1) / i3 + 1;

                    ixwk = 1;
                    nset = 8;
                    ibwk = ixwk + jdw * kdw * idw * nset;
                    ixgk = ibwk + jdim * kdim * idim * 2;
                    ixvk = ixgk + jdw * kdw * idw * 4;
                    nroom = nwork - (ixvk + jdw * kdw * idw * 5 - 1);
                    if (nroom < 0) {
                        if (myid == myhost) {
                            fortran_write_unit(11, " not enough memory for plot3d\n");
                            {
                                char buf[120];
                                std::snprintf(buf, 120, " have, need = %12d%12d\n",
                                              nwork, nwork - nroom);
                                fortran_write_unit(11, "%s", buf);
                            }
                            fortran_write_unit(11, " not writing out prout file\n");
                        }
                    } else {
                        int iover_v = iovrlp(nbl);
                        FortranArray4DRef<double> q_ref(&w(lq), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qi0_ref(&w(lqi0), jdim, kdim, 5, 4);
                        FortranArray4DRef<double> qj0_ref(&w(lqj0), kdim, idim-1, 5, 4);
                        FortranArray4DRef<double> qk0_ref(&w(lqk0), jdim, idim-1, 5, 4);
                        FortranArray3DRef<double> x_ref(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> y_ref(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> z_ref(&w(lz), jdim, kdim, idim);
                        FortranArray4DRef<double> xw_ref(&wk(ixwk), jdw, kdw, idw, nset);
                        FortranArray4DRef<double> blank2_ref(&wk(ibwk), jdim, kdim, idim, 2);
                        FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                        FortranArray4DRef<double> xg_ref(&wk(ixgk), jdw, kdw, idw, 4);
                        FortranArray3DRef<double> vis_ref(&w(lvis), jdim, kdim, idim);
                        FortranArray3DRef<double> bcj_ref(&w(lbcj), kdim, idim-1, 2);
                        FortranArray3DRef<double> bck_ref(&w(lbck), jdim, idim-1, 2);
                        FortranArray3DRef<double> bci_ref(&w(lbci), jdim, kdim, 2);
                        FortranArray4DRef<double> vj0_ref(&w(lvj0), kdim, idim-1, 1, 4);
                        FortranArray4DRef<double> vk0_ref(&w(lvk0), jdim, idim-1, 1, 4);
                        FortranArray4DRef<double> vi0_ref(&w(lvi0), jdim, kdim, 1, 4);
                        FortranArray4DRef<double> xv_ref(&wk(ixvk), jdw, kdw, idw, 5);
                        FortranArray4DRef<double> sj_ref(&w(lsj), jdim, kdim, idim-1, 5);
                        FortranArray4DRef<double> sk_ref(&w(lsk), jdim, kdim, idim-1, 5);
                        FortranArray4DRef<double> si_ref(&w(lsi), jdim, kdim, idim, 5);
                        FortranArray3DRef<double> vol_ref(&w(lvol), jdim, kdim, idim-1);
                        plot3d_ns::plot3d(jdim, kdim, idim, i1, i2, i3, j1, j2, j3,
                                          k1, k2, k3, q_ref, qi0_ref, qj0_ref, qk0_ref,
                                          x_ref, y_ref, z_ref, xw_ref, blank2_ref,
                                          blank_ref, xg_ref, iflag, vis_ref,
                                          iover_v, nbl, nmap, bcj_ref, bck_ref, bci_ref,
                                          vj0_ref, vk0_ref, vi0_ref, ifunc, n,
                                          jdw, kdw, idw, nplots, jdimg, kdimg, idimg,
                                          nblcg, jsg, ksg, isg, jeg, keg, ieg,
                                          ninter, iindex, intmax, nsub1, maxxe,
                                          nblk, nbli, limblk, isva, nblon, mxbli,
                                          thetay, maxbl, maxgr, myid, myhost, mycomm,
                                          mblk2nd, inpl3d, nblock, nblkpt,
                                          xv_ref, sj_ref, sk_ref, si_ref, vol_ref, nset);
                    }
                }

            } else {


                // cell center data
                if (mblk2nd(nbl) == myid || myid == myhost) {
                    lead_ns::lead(nbl, lw, lw2, maxbl);

                    if (myid == myhost) {
                        FILE* f17 = fortran_get_unit(17);
                        fprintf(f17, "\n\n BLOCK%4d  (GRID%4d)     IDIM,JDIM,KDIM=%5d%5d%5d\n",
                                nbl, igridg(nbl), idim, jdim, kdim);
                    }

                    i2 = std::min(idim-1, i2);
                    j2 = std::min(jdim-1, j2);
                    k2 = std::min(kdim-1, k2);
                    i1 = std::min(idim-1, i1);
                    j1 = std::min(jdim-1, j1);
                    k1 = std::min(kdim-1, k1);
                    jdw = (j2 - j1) / j3 + 1;
                    kdw = (k2 - k1) / k3 + 1;
                    idw = (i2 - i1) / i3 + 1;

                    ixwk = 1;
                    ibwk = ixwk + jdw * kdw * idw * 5;
                    ixgk = ibwk + jdim * kdim * idim;
                    nroom = nwork - (ixgk + jdw * kdw * idw * 4 - 1);
                    if (nroom < 0) {
                        if (myid == myhost) {
                            fortran_write_unit(11, " not enough memory for plot3c\n");
                            {
                                char buf[120];
                                std::snprintf(buf, 120, " have, need = %12d%12d\n",
                                              nwork, nwork - nroom);
                                fortran_write_unit(11, "%s", buf);
                            }
                            fortran_write_unit(11, " not writing out plot3d files\n");
                        }
                    } else {
                        int iover_v = iovrlp(nbl);
                        FortranArray4DRef<double> q_ref(&w(lq), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qi0_ref(&w(lqi0), jdim, kdim, 5, 4);
                        FortranArray4DRef<double> qj0_ref(&w(lqj0), kdim, idim-1, 5, 4);
                        FortranArray4DRef<double> qk0_ref(&w(lqk0), jdim, idim-1, 5, 4);
                        FortranArray3DRef<double> x_ref(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> y_ref(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> z_ref(&w(lz), jdim, kdim, idim);
                        FortranArray4DRef<double> xw_ref(&wk(ixwk), jdw, kdw, idw, 5);
                        FortranArray3DRef<double> blank2_ref(&wk(ibwk), jdim, kdim, idim);
                        FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                        FortranArray4DRef<double> xg_ref(&wk(ixgk), jdw, kdw, idw, 4);
                        FortranArray3DRef<double> vis_ref(&w(lvis), jdim, kdim, idim);
                        FortranArray4DRef<double> vi0_ref(&w(lvi0), jdim, kdim, 1, 4);
                        FortranArray4DRef<double> vj0_ref(&w(lvj0), kdim, idim-1, 1, 4);
                        FortranArray4DRef<double> vk0_ref(&w(lvk0), jdim, idim-1, 1, 4);
                        FortranArray3DRef<double> smin_ref(&w(lsnk0), jdim, kdim, idim);
                        plot3c_ns::plot3c(jdim, kdim, idim, i1, i2, i3, j1, j2, j3,
                                          k1, k2, k3, q_ref, qi0_ref, qj0_ref, qk0_ref,
                                          x_ref, y_ref, z_ref, xw_ref, blank2_ref,
                                          blank_ref, xg_ref, iflag, vis_ref,
                                          vi0_ref, vj0_ref, vk0_ref,
                                          iover_v, nbl, nmap, smin_ref, ifunc, n,
                                          jdw, kdw, idw, nplots, jdimg, kdimg, idimg,
                                          nblcg, jsg, ksg, isg, jeg, keg, ieg,
                                          ninter, iindex, intmax, nsub1, maxxe,
                                          nblk, nbli, limblk, isva, nblon, mxbli,
                                          thetay, maxbl, maxgr, myid, myhost, mycomm,
                                          mblk2nd, inpl3d, nblock, nblkpt, iprnsurf);
                    }
                }
            } // end if inpr(n,2)==0 else

        } // end do 80 n=1,nprint

        if (myid == myhost) {
            if (lhdr > 0) fortran_write_unit(11, "\n");
        }

    } // end if nprint > 0

    return;
}

} // namespace qoutavg_ns
