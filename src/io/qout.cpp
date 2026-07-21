// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "qout.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include "bc.h"
#include "bc_blkint.h"
#include "bc_embed.h"
#include "bc_patch.h"
#include "bc_period.h"
#include "lead.h"
#include "bc_xmera.h"
#include "qface.h"
#include "blnkfr.h"
#include "plot3d.h"
#include "plot3c.h"
#include "plot3t.h"
#include "termn8.h"
#include "module_kwstm.h"
#include "module_profileout.h"
#include "module_contour.h"
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <cmath>

namespace qout_ns {

void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork,
        double& cl, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
        int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg,
        FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
        FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg,
        FortranArray1DRef<int> nblg, FortranArray1DRef<int> nbci0,
        FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
        FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim,
        FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo,
        FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
        FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej,
        FortranArray3DRef<int> bcfilek, FortranArray3DRef<int> lwdat,
        int& myid, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
        FortranArray1DRef<int> kdimg, FortranArray1DRef<char[80]> bcfiles,
        int& mxbcfil, int& nummem)
{
    bc_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, cl, nou, bou, nbuf, ibufdim,
              maxbl, maxgr, maxseg, itrans, irotat, idefrm, igridg, nblg,
              nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
              ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek,
              lwdat, myid, idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
}

void bc_blkint(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
               FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
               int& nwork, int& maxbl, int& maxgr, int& mxbli,
               FortranArray1DRef<int> iadvance, FortranArray1DRef<double> geom_miss,
               double& epsc0, FortranArray2DRef<int> nblk, int& nbli,
               FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva,
               FortranArray1DRef<int> nblon, FortranArray1DRef<int> jdimg,
               FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg,
               FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> isav_blk,
               FortranArray1DRef<int> ireq_ar, FortranArray1DRef<int> index_ar,
               FortranArray1DRef<int> ireq_snd, FortranArray2DRef<int> keep_trac,
               FortranArray1DRef<int> keep_trac2, FortranArray1DRef<int> nou,
               FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
               int& myid, int& myhost, int& mycomm, FortranArray2DRef<int> istat2,
               int& istat_size, int& nummem)
{
    bc_blkint_ns::bc_blkint(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, maxbl, maxgr,
                             mxbli, iadvance, geom_miss, epsc0, nblk, nbli, limblk,
                             isva, nblon, jdimg, kdimg, idimg, mblk2nd, isav_blk,
                             ireq_ar, index_ar, ireq_snd, keep_trac, keep_trac2,
                             nou, bou, nbuf, ibufdim, myid, myhost, mycomm,
                             istat2, istat_size, nummem);
}

void bc_embed(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
              FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
              int& nwork, int& maxbl, int& maxgr, int& lbcemb,
              FortranArray1DRef<int> iadvance, FortranArray1DRef<int> idimg,
              FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg,
              FortranArray2DRef<int> isav_emb, FortranArray1DRef<int> ireq_ar,
              FortranArray1DRef<int> ireq_snd, FortranArray1DRef<int> index_ar,
              FortranArray2DRef<int> keep_trac, FortranArray1DRef<int> keep_trac2,
              int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd,
              FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
              int& nbuf, int& ibufdim, FortranArray2DRef<int> iviscg,
              FortranArray2DRef<int> istat2, int& istat_size, int& nummem)
{
    bc_embed_ns::bc_embed(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, maxbl, maxgr,
                          lbcemb, iadvance, idimg, jdimg, kdimg, isav_emb,
                          ireq_ar, ireq_snd, index_ar, keep_trac, keep_trac2,
                          myid, myhost, mycomm, mblk2nd, nou, bou, nbuf, ibufdim,
                          iviscg, istat2, istat_size, nummem);
}

void bc_patch(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
              FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
              int& nwork, int& maxbl, int& maxgr, int& intmax, int& nsub1, int& maxxe,
              FortranArray1DRef<int> iadvance, FortranArray1DRef<int> jdimg,
              FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg,
              int& ninter, FortranArray2DRef<double> windex, FortranArray2DRef<int> iindex,
              FortranArray1DRef<int> nblkpt, FortranArray2DRef<double> dthetxx,
              FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz,
              FortranArray2DRef<int> isav_pat, FortranArray3DRef<int> isav_pat_b,
              FortranArray1DRef<int> ireq_ar, FortranArray1DRef<int> ireq_snd,
              FortranArray1DRef<int> index_ar, FortranArray3DRef<int> keep_trac,
              FortranArray2DRef<int> keep_trac2, int& myid, int& myhost, int& mycomm,
              FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou,
              FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
              FortranArray2DRef<int> istat2, int& istat_size, int& nummem)
{
    bc_patch_ns::bc_patch(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, maxbl, maxgr,
                          intmax, nsub1, maxxe, iadvance, jdimg, kdimg, idimg,
                          ninter, windex, iindex, nblkpt, dthetxx, dthetyy, dthetzz,
                          isav_pat, isav_pat_b, ireq_ar, ireq_snd, index_ar,
                          keep_trac, keep_trac2, myid, myhost, mycomm, mblk2nd,
                          nou, bou, nbuf, ibufdim, istat2, istat_size, nummem);
}

void bc_period(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
               FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
               int& nwork, int& maxbl, int& maxgr, int& maxseg,
               FortranArray1DRef<int> iadvance, FortranArray3DRef<int> bcfilei,
               FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek,
               FortranArray3DRef<int> lwdat, FortranArray1DRef<double> xorig,
               FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig,
               FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg,
               FortranArray1DRef<int> idimg, int& lbcprd,
               FortranArray2DRef<int> isav_prd, FortranArray1DRef<double> period_miss,
               double& epsrot, FortranArray1DRef<int> ireq_ar,
               FortranArray1DRef<int> index_ar, FortranArray1DRef<int> ireq_snd,
               FortranArray2DRef<int> keep_trac, FortranArray1DRef<int> keep_trac2,
               int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd,
               FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
               int& nbuf, int& ibufdim, FortranArray2DRef<int> istat2,
               int& istat_size, FortranArray1DRef<char[80]> bcfiles,
               int& mxbcfil, int& nummem)
{
    bc_period_ns::bc_period(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, maxbl, maxgr,
                            maxseg, iadvance, bcfilei, bcfilej, bcfilek, lwdat,
                            xorig, yorig, zorig, jdimg, kdimg, idimg, lbcprd,
                            isav_prd, period_miss, epsrot, ireq_ar, index_ar,
                            ireq_snd, keep_trac, keep_trac2, myid, myhost, mycomm,
                            mblk2nd, nou, bou, nbuf, ibufdim, istat2, istat_size,
                            bcfiles, mxbcfil, nummem);
}

void qout(int& iseq, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
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
          FortranArray1DRef<int> idefrm, int& nblock,
          FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg,
          FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg,
          FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg,
          FortranArray1DRef<int> nblg, FortranArray1DRef<double> clw,
          int& ncycmax, int& nplot3d, FortranArray2DRef<int> inpl3d,
          int& ip3dsurf, int& nprint, FortranArray2DRef<int> inpr,
          FortranArray1DRef<int> iadvance, int& mycomm, int& myid, int& myhost,
          FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou,
          FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
          int& maxbl, int& maxgr, int& maxseg, int& iitot,
          FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg,
          FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg,
          FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg,
          int& ninter, FortranArray2DRef<double> windex,
          FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt,
          int& intmax, int& nsub1, int& maxxe,
          FortranArray2DRef<int> nblk, int& nbli,
          FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva,
          FortranArray1DRef<int> nblon, int& mxbli,
          FortranArray1DRef<double> thetay, FortranArray2DRef<int> ip3ddim,
          FortranArray1DRef<int> nmap, FortranArray1DRef<int> iwk, int& iwork,
          FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig,
          FortranArray1DRef<double> zorig, FortranArray1DRef<double> period_miss,
          FortranArray1DRef<double> geom_miss, double& epsc0, double& epsrot,
          FortranArray2DRef<int> isav_blk, FortranArray2DRef<int> isav_pat,
          FortranArray3DRef<int> isav_pat_b, FortranArray2DRef<int> isav_emb,
          FortranArray2DRef<int> isav_prd, int& lbcprd, int& lbcemb,
          FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy,
          FortranArray2DRef<double> dthetzz, FortranArray1DRef<int> nblcg,
          int& lfgm, FortranArray2DRef<int> istat2_bl,
          FortranArray2DRef<int> istat2_pa, FortranArray2DRef<int> istat2_pe,
          FortranArray2DRef<int> istat2_em, int& istat_size,
          FortranArray1DRef<double> vormax, FortranArray1DRef<int> ivmax,
          FortranArray1DRef<int> jvmax, FortranArray1DRef<int> kvmax,
          FortranArray1DRef<char[80]> bcfiles, int& mxbcfil,
          int& iprnsurf, int& nummem)
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
    int& lblk      = cmn_ginfo.lblk;
    int& lsnk0     = cmn_ginfo.lsnk0;
    int& lbcj      = cmn_ginfo.lbcj;
    int& lbck      = cmn_ginfo.lbck;
    int& lbci      = cmn_ginfo.lbci;
    int& lvj0      = cmn_ginfo.lvj0;
    int& lvk0      = cmn_ginfo.lvk0;
    int& lvi0      = cmn_ginfo.lvi0;
    int& lxib      = cmn_ginfo.lxib;
    int& lux       = cmn_ginfo.lux;
    int& lcmuv     = cmn_ginfo.lcmuv;
    int& lvolj0    = cmn_ginfo.lvolj0;
    int& lvolk0    = cmn_ginfo.lvolk0;
    int& lvoli0    = cmn_ginfo.lvoli0;

    float* title   = cmn_info.title;
    float& xmach   = cmn_info.xmach;
    float& alpha_f = cmn_info.alpha;
    float& beta_f  = cmn_info.beta;
    int& iipv      = cmn_info.iipv;
    int& mseq      = cmn_info.mseq;
    int* levelt    = cmn_info.levelt;

    int& ivmx      = cmn_maxiv.ivmx;
    int& level     = cmn_mgrd.level;
    int& lglobal   = cmn_mgrd.lglobal;
    float& reue    = cmn_reyue.reue;
    float& tinf    = cmn_reyue.tinf;
    int& i2d       = cmn_twod.i2d;
    float& time_f  = cmn_unst.time;
    int& iunst     = cmn_unst.iunst;
    int& icall1    = cmn_moov.icall1;
    int& lhdr      = cmn_moov.lhdr;
    int& nframes   = cmn_moov.nframes;
    float& radtodeg = cmn_conversion.radtodeg;
    int& ifunct    = cmn_plot3dtyp.ifunct;
    int& lowmem_ux = cmn_memry.lowmem_ux;

    // Local variables
    int level_sav, nttuse, lres, nsafe, mneed;
    int iwk1, iwk2, iwk3, iwk4, iwk5, iwk6, iwork1;
    int int_updt, ncount, np3d, ifunc, ifuncuse;
    int nbl, n, iii, i1, i2, i3, j1, j2, j3, k1, k2, k3;
    int jdw, kdw, idw, ixwk, ibwk, ixgk, ixvk, ibwk2, ibwk3, nroom, nset;
    int ifuncdim, nnn, m;
    double clwuse = 0.0;
    float alphaw, betaw;
    int ntime_1 = 1;
    int ierrflg_m1 = -1;


    if (nplot3d > 0 || nprint > 0) {
        //
        // update all bc's for level lglobal (and above, if embeded),
        // even if not all blocks at that level get output to plot3d or
        // printout files
        //
        level_sav = level;
        lglobal = lfgm - (mseq - iseq);
        //
        for (level = lglobal; level <= levelt[iseq-1]; level++) {
            //
            if (iipv > 0) {
                nttuse = std::max(cmn_info.ntt - 1, 1);
                clwuse = clw(nttuse);
            }

            for (nbl = 1; nbl <= nblock; nbl++) {
                // need to call even for blocks not advanced
                if (level == levelg(nbl)) {
                    if (mblk2nd(nbl) == myid) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        bc_ns::bc(ntime_1, nbl, lw, lw2, w, mgwk, wk, nwork,
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
            //
            // update periodic boundary conditions
            //
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
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "stopping...not enough integer work space for subroutine bc_period");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                termn8_ns::termn8(myid, ierrflg_m1, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                FortranArray1DRef<int> iwk_1(&iwk(iwk1), mneed, 1);
                FortranArray1DRef<int> iwk_2(&iwk(iwk2), mneed, 1);
                FortranArray1DRef<int> iwk_3(&iwk(iwk3), mneed, 1);
                FortranArray2DRef<int> iwk_4(&iwk(iwk4), lbcprd, 10);
                FortranArray1DRef<int> iwk_5(&iwk(iwk5), mneed, 1);
                bc_period_ns::bc_period(ntime_1, nbl, lw, lw2, w, mgwk, wk, nsafe, maxbl,
                                        maxgr, maxseg, iadvance, bcfilei, bcfilej,
                                        bcfilek, lwdat, xorig, yorig, zorig, jdimg, kdimg,
                                        idimg, lbcprd, isav_prd,
                                        period_miss, epsrot, iwk_1, iwk_2,
                                        iwk_3, iwk_4, iwk_5, myid, myhost,
                                        mycomm, mblk2nd, nou, bou, nbuf, ibufdim,
                                        istat2_pe, istat_size, bcfiles, mxbcfil, nummem);
            }

            //
            // update embeded-grid boundary conditions
            //
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
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "stopping...not enough integer work space for subroutine bc_embed");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                termn8_ns::termn8(myid, ierrflg_m1, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                FortranArray1DRef<int> iwk_1(&iwk(iwk1), mneed, 1);
                FortranArray1DRef<int> iwk_2(&iwk(iwk2), mneed, 1);
                FortranArray1DRef<int> iwk_3(&iwk(iwk3), mneed, 1);
                FortranArray2DRef<int> iwk_4(&iwk(iwk4), lbcemb, 6);
                FortranArray1DRef<int> iwk_5(&iwk(iwk5), mneed, 1);
                bc_embed_ns::bc_embed(ntime_1, nbl, lw, lw2, w, mgwk, wk, nsafe, maxbl,
                                      maxgr, lbcemb, iadvance, idimg, jdimg,
                                      kdimg, isav_emb, iwk_1,
                                      iwk_2, iwk_3, iwk_4, iwk_5,
                                      myid, myhost, mycomm, mblk2nd, nou, bou, nbuf,
                                      ibufdim, iviscg, istat2_em, istat_size, nummem);
            }
            //
            // update 1-1 block boundary conditions
            //
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
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "stopping...not enough integer work space for subroutine bc_blkint");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                termn8_ns::termn8(myid, ierrflg_m1, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                FortranArray1DRef<int> iwk_1(&iwk(iwk1), mneed, 1);
                FortranArray1DRef<int> iwk_2(&iwk(iwk2), mneed, 1);
                FortranArray1DRef<int> iwk_3(&iwk(iwk3), mneed, 1);
                FortranArray2DRef<int> iwk_4(&iwk(iwk4), mxbli, 10);
                FortranArray1DRef<int> iwk_5(&iwk(iwk5), mneed, 1);
                bc_blkint_ns::bc_blkint(ntime_1, nbl, lw, lw2, w, mgwk, wk, nsafe, maxbl,
                                        maxgr, mxbli, iadvance, geom_miss, epsc0, nblk,
                                        nbli, limblk, isva, nblon, jdimg, kdimg, idimg,
                                        mblk2nd, isav_blk, iwk_1,
                                        iwk_2, iwk_3, iwk_4, iwk_5,
                                        nou, bou, nbuf, ibufdim, myid, myhost, mycomm,
                                        istat2_bl, istat_size, nummem);
            }

            //
            // update patch-grid boundary conditions
            //
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
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "stopping...not enough integer work space for subroutine bc_patch");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                termn8_ns::termn8(myid, ierrflg_m1, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                FortranArray1DRef<int> iwk_1(&iwk(iwk1), mneed, 1);
                FortranArray1DRef<int> iwk_2(&iwk(iwk2), mneed, 1);
                FortranArray1DRef<int> iwk_3(&iwk(iwk3), mneed, 1);
                FortranArray3DRef<int> iwk_4(&iwk(iwk4), intmax, nsub1, 6);
                FortranArray2DRef<int> iwk_5(&iwk(iwk5), mneed, 2);
                bc_patch_ns::bc_patch(ntime_1, nbl, lw, lw2, w, mgwk, wk, nsafe, maxbl,
                                      maxgr, intmax, nsub1, maxxe, iadvance, jdimg, kdimg,
                                      idimg, ninter, windex, iindex, nblkpt, dthetxx,
                                      dthetyy, dthetzz, isav_pat, isav_pat_b,
                                      iwk_1, iwk_2, iwk_3,
                                      iwk_4, iwk_5, myid, myhost, mycomm,
                                      mblk2nd, nou, bou, nbuf, ibufdim,
                                      istat2_pa, istat_size, nummem);
            }
            //
            // update chimera boundary conditions
            //
            // don't update interior (fringe) points if cell center data
            // is to be output
            //
            int_updt = 1;
            for (n = 1; n <= nplot3d; n++) {
                if (inpl3d(n, 2) > 0) {
                    int_updt = 0;
                }
            }
            //
            for (nbl = 1; nbl <= nblock; nbl++) {
                if (iadvance(nbl) >= 0) {
                    if (level == levelg(nbl)) {
                        if (mblk2nd(nbl) == myid) {
                            lead_ns::lead(nbl, lw, lw2, maxbl);
                            lres  = 1;
                            nsafe = nwork - lres + 1;
                            bc_xmera_ns::bc_xmera(ntime_1, nbl, lw, lw2, w, mgwk, wk,
                                                  nsafe, maxbl, iitot, iviscg, iovrlp,
                                                  lbg, ibpntsg, qb, iibg, kkbg, jjbg,
                                                  ibcg, nou, bou, nbuf, ibufdim,
                                                  int_updt, nummem);
                        }
                    }
                }
            }
            //
            // call qface to install face-center values in
            // the qi0/qj0/qk0 arrays
            //
            for (nbl = 1; nbl <= nblock; nbl++) {
                if (level == levelg(nbl)) {
                    if (mblk2nd(nbl) == myid) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        int ldim = 5;
                        FortranArray4DRef<double> w_lq(   &w(lq),    jdim, kdim, idim, ldim);
                        FortranArray4DRef<double> w_lqj0( &w(lqj0),  jdim, kdim, 5,    ldim);
                        FortranArray4DRef<double> w_lqk0( &w(lqk0),  jdim, kdim, 5,    ldim);
                        FortranArray4DRef<double> w_lqi0( &w(lqi0),  jdim, kdim, 5,    ldim);
                        FortranArray3DRef<double> w_lbcj( &w(lbcj),  jdim, kdim, idim);
                        FortranArray3DRef<double> w_lbck( &w(lbck),  jdim, kdim, idim);
                        FortranArray3DRef<double> w_lbci( &w(lbci),  jdim, kdim, idim);
                        FortranArray3DRef<double> w_lblk( &w(lblk),  jdim, kdim, idim);
                        qface_ns::qface(jdim, kdim, idim, w_lq, w_lqj0,
                                        w_lqk0, w_lqi0, w_lbcj, w_lbck,
                                        w_lbci, w_lblk, ldim);
                        if (ivmx >= 2) {
                            ldim = 1;
                            FortranArray4DRef<double> w_lvis(&w(lvis), jdim, kdim, idim, ldim);
                            FortranArray4DRef<double> w_lvj0(&w(lvj0), jdim, kdim, 5,    ldim);
                            FortranArray4DRef<double> w_lvk0(&w(lvk0), jdim, kdim, 5,    ldim);
                            FortranArray4DRef<double> w_lvi0(&w(lvi0), jdim, kdim, 5,    ldim);
                            FortranArray3DRef<double> w_lbcj2(&w(lbcj), jdim, kdim, idim);
                            FortranArray3DRef<double> w_lbck2(&w(lbck), jdim, kdim, idim);
                            FortranArray3DRef<double> w_lbci2(&w(lbci), jdim, kdim, idim);
                            FortranArray3DRef<double> w_lblk2(&w(lblk), jdim, kdim, idim);
                            qface_ns::qface(jdim, kdim, idim, w_lvis, w_lvj0,
                                            w_lvk0, w_lvi0, w_lbcj2, w_lbck2,
                                            w_lbci2, w_lblk2, ldim);
                        }
                    }
                }
            }
            //
        } // end do level=lglobal,levelt(iseq)
        //
        level = level_sav;
        //
    } // end if (nplot3d > 0 || nprint > 0)


    if (nplot3d <= 0) goto label_231;
    //
    ncount    = 0;
    np3d      = nplot3d;
    ifunc     = 0;
    ifuncuse  = 0;
    //
    // if zone has function file output, all must, and all
    // must have the same function output
    //
    if (std::abs(inpl3d(1, 2)) > 2) {
        ifunc = inpl3d(1, 2);
        ifuncuse = 1;
    }
    if (std::abs(inpl3d(1, 2)) == 2) {
        ifuncuse = ifunct;
    }
    //
    for (n = 1; n <= nplot3d; n++) {
        if (n == 1) {
            if (myid == myhost) {
                if (ibin == 0) {
                    if (icall1 == 0) {
                        FILE* f3 = fortran_get_unit(3);
                        fprintf(f3, "%5d\n", np3d);
                    }
                    FILE* f4 = fortran_get_unit(4);
                    fprintf(f4, "%5d\n", np3d);
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
        //
        nbl = inpl3d(n, 1);
        //
        if (nbl > nblock) {
            if (myid == myhost) {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, " Block%3d does not exist.  No plot3d output printed.\n", nbl);
            }
            continue; // go to 60
        }
        //
        ncount = ncount + 1;
        i1 = inpl3d(n, 3);
        i2 = inpl3d(n, 4);
        i3 = inpl3d(n, 5);
        j1 = inpl3d(n, 6);
        j2 = inpl3d(n, 7);
        j3 = inpl3d(n, 8);
        k1 = inpl3d(n, 9);
        k2 = inpl3d(n, 10);
        k3 = inpl3d(n, 11);
        if (inpl3d(n, 2) > 0) {
            // cell center dimensions
            lead_ns::lead(nbl, lw, lw2, maxbl);
            i2 = std::min(idim - 1, i2);
            j2 = std::min(jdim - 1, j2);
            k2 = std::min(kdim - 1, k2);
            i1 = std::min(idim - 1, i1);
            j1 = std::min(jdim - 1, j1);
            k1 = std::min(kdim - 1, k1);
        }
        ip3ddim(1, ncount) = (i2 - i1) / i3 + 1;
        ip3ddim(2, ncount) = (j2 - j1) / j3 + 1;
        ip3ddim(3, ncount) = (k2 - k1) / k3 + 1;
    } // end do n=1,nplot3d (loop 60)

    //
    if (myid == myhost) {
        if (ibin == 0) {
            if (i2d == 0) {
                if (icall1 == 0) {
                    FILE* f3 = fortran_get_unit(3);
                    for (int nn = 1; nn <= ncount; nn++)
                        for (int ii = 1; ii <= 3; ii++)
                            fprintf(f3, "%5d", ip3ddim(ii, nn));
                    fprintf(f3, "\n");
                }
                if (ifuncuse == 0) {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++)
                        for (int ii = 1; ii <= 3; ii++)
                            fprintf(f4, "%5d", ip3ddim(ii, nn));
                    fprintf(f4, "\n");
                } else {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        for (int ii = 1; ii <= 3; ii++)
                            fprintf(f4, "%5d", ip3ddim(ii, nn));
                        fprintf(f4, "%5d", ifuncuse);
                    }
                    fprintf(f4, "\n");
                }
            } else {
                // i2d != 0: 2D case, write indices 2 and 3 only
                if (icall1 == 0) {
                    FILE* f3 = fortran_get_unit(3);
                    for (int nn = 1; nn <= ncount; nn++)
                        for (int ii = 2; ii <= 3; ii++)
                            fprintf(f3, "%5d", ip3ddim(ii, nn));
                    fprintf(f3, "\n");
                }
                if (ifuncuse == 0) {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++)
                        for (int ii = 2; ii <= 3; ii++)
                            fprintf(f4, "%5d", ip3ddim(ii, nn));
                    fprintf(f4, "\n");
                } else {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        for (int ii = 2; ii <= 3; ii++)
                            fprintf(f4, "%5d", ip3ddim(ii, nn));
                        fprintf(f4, "%5d", ifuncuse);
                    }
                    fprintf(f4, "\n");
                }
            }
        } else {
            // binary
            if (i2d == 0) {
                if (icall1 == 0) {
                    FILE* f3 = fortran_get_unit(3);
                    for (int nn = 1; nn <= ncount; nn++)
                        for (int ii = 1; ii <= 3; ii++) {
                            int val = ip3ddim(ii, nn);
                            fwrite(&val, sizeof(int), 1, f3);
                        }
                }
                if (ifuncuse == 0) {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++)
                        for (int ii = 1; ii <= 3; ii++) {
                            int val = ip3ddim(ii, nn);
                            fwrite(&val, sizeof(int), 1, f4);
                        }
                } else {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        for (int ii = 1; ii <= 3; ii++) {
                            int val = ip3ddim(ii, nn);
                            fwrite(&val, sizeof(int), 1, f4);
                        }
                        fwrite(&ifuncuse, sizeof(int), 1, f4);
                    }
                }
            } else {
                // i2d != 0: 2D binary
                if (icall1 == 0) {
                    FILE* f3 = fortran_get_unit(3);
                    for (int nn = 1; nn <= ncount; nn++)
                        for (int ii = 2; ii <= 3; ii++) {
                            int val = ip3ddim(ii, nn);
                            fwrite(&val, sizeof(int), 1, f3);
                        }
                }
                if (ifuncuse == 0) {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++)
                        for (int ii = 2; ii <= 3; ii++) {
                            int val = ip3ddim(ii, nn);
                            fwrite(&val, sizeof(int), 1, f4);
                        }
                } else {
                    FILE* f4 = fortran_get_unit(4);
                    for (int nn = 1; nn <= ncount; nn++) {
                        for (int ii = 2; ii <= 3; ii++) {
                            int val = ip3ddim(ii, nn);
                            fwrite(&val, sizeof(int), 1, f4);
                        }
                        fwrite(&ifuncuse, sizeof(int), 1, f4);
                    }
                }
            }
        }
    }

    //
    // correspondence between global block number and plot3d block
    // number stored in nmap(n) for n=1,nblock
    //
    for (n = 1; n <= nblock; n++) {
        nmap(n) = 1;
        for (nnn = 1; nnn <= nplot3d; nnn++) {
            m = inpl3d(nnn, 1);
            if (n == m) nmap(n) = nnn;
        }
    }
    //
    if (myid == myhost) {
        if (lhdr > 0) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, " \n");
        }
    }
    //
    for (n = 1; n <= nplot3d; n++) {
        //
        nbl = inpl3d(n, 1);
        //
        if (nbl > nblock || nbl <= 0) continue; // go to 70
        //
        if (iblnkfr == 0) {
            //
            // temporarily set blank values at fringe points (not holes)
            // to 1 for plotting purposes (helps reduce gaps in plots)
            //
            lead_ns::lead(nbl, lw, lw2, maxbl);
            if (mblk2nd(nbl) == myid && iovrlp(nbl) != 0) {
                double blnkval = 1.0;
                blnkfr_ns::blnkfr(nbl, iibg, kkbg, jjbg, ibpntsg, lbg, iitot,
                                   FortranArray3DRef<double>(&w(lblk), jdim, kdim, idim),
                                   jdim, kdim, idim, maxbl, blnkval);
            }
        }
        //
        i1 = inpl3d(n, 3);
        i2 = inpl3d(n, 4);
        i3 = inpl3d(n, 5);
        j1 = inpl3d(n, 6);
        j2 = inpl3d(n, 7);
        j3 = inpl3d(n, 8);
        k1 = inpl3d(n, 9);
        k2 = inpl3d(n, 10);
        k3 = inpl3d(n, 11);
        //
        int iflag = 1;
        //
        if (inpl3d(n, 2) <= 0) {
            //
            // grid point data
            //
            if (mblk2nd(nbl) == myid || myid == myhost) {
                //
                lead_ns::lead(nbl, lw, lw2, maxbl);
                //
                jdw = (j2 - j1) / j3 + 1;
                kdw = (k2 - k1) / k3 + 1;
                idw = (i2 - i1) / i3 + 1;
                //
                // check storage availability
                //
                ixwk = 1;
                nset = 5;
                ibwk = ixwk + jdw * kdw * idw * nset;
                ixgk = ibwk + jdim * kdim * idim * 2;
                ixvk = ixgk + jdw * kdw * idw * 4;
                nroom = nwork - (ixvk + jdw * kdw * idw * 5 - 1);
                if (nroom < 0) {
                    if (myid == myhost) {
                        FILE* f11 = fortran_get_unit(11);
                        fprintf(f11, " not enough memory for plot3d\n");
                        fprintf(f11, " have, need = %12d%12d\n", nwork, nwork - nroom);
                        fprintf(f11, " Aborting - not calling plot3d\n");
                        std::exit(0);
                    }
                } else {
                    FortranArray4DRef<double> p3d_q  (&w(lq),   jdim, kdim, idim, 5);
                    FortranArray4DRef<double> p3d_qi0 (&w(lqi0), jdim, kdim, 5, 4);
                    FortranArray4DRef<double> p3d_qj0 (&w(lqj0), kdim, idim - 1, 5, 4);
                    FortranArray4DRef<double> p3d_qk0 (&w(lqk0), jdim, idim - 1, 5, 4);
                    FortranArray3DRef<double> p3d_x   (&w(lx),   jdim, kdim, idim);
                    FortranArray3DRef<double> p3d_y   (&w(ly),   jdim, kdim, idim);
                    FortranArray3DRef<double> p3d_z   (&w(lz),   jdim, kdim, idim);
                    FortranArray4DRef<double> p3d_xw  (&wk(ixwk), jdw, kdw, idw, nset);
                    FortranArray4DRef<double> p3d_bk2 (&wk(ibwk), jdim, kdim, idim, 2);
                    FortranArray3DRef<double> p3d_blk (&w(lblk), jdim, kdim, idim);
                    FortranArray4DRef<double> p3d_xg  (&wk(ixgk), jdw, kdw, idw, 4);
                    FortranArray3DRef<double> p3d_vis (&w(lvis), jdim, kdim, idim);
                    FortranArray3DRef<double> p3d_bcj (&w(lbcj), kdim, idim - 1, 2);
                    FortranArray3DRef<double> p3d_bck (&w(lbck), jdim, idim - 1, 2);
                    FortranArray3DRef<double> p3d_bci (&w(lbci), jdim, kdim, 2);
                    FortranArray4DRef<double> p3d_vj0 (&w(lvj0), kdim, idim - 1, 1, 4);
                    FortranArray4DRef<double> p3d_vk0 (&w(lvk0), jdim, idim - 1, 1, 4);
                    FortranArray4DRef<double> p3d_vi0 (&w(lvi0), jdim, kdim, 1, 4);
                    FortranArray4DRef<double> p3d_xv  (&wk(ixvk), jdw, kdw, idw, 5);
                    FortranArray4DRef<double> p3d_sj  (&w(lsj),  jdim, kdim, idim - 1, 5);
                    FortranArray4DRef<double> p3d_sk  (&w(lsk),  jdim, kdim, idim - 1, 5);
                    FortranArray4DRef<double> p3d_si  (&w(lsi),  jdim, kdim, idim, 5);
                    FortranArray3DRef<double> p3d_vol (&w(lvol), jdim, kdim, idim - 1);
                    plot3d_ns::plot3d(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                                     p3d_q, p3d_qi0, p3d_qj0, p3d_qk0, p3d_x,
                                     p3d_y, p3d_z, p3d_xw, p3d_bk2,
                                     p3d_blk, p3d_xg, iflag, p3d_vis, iovrlp(nbl),
                                     nbl, nmap, p3d_bcj, p3d_bck, p3d_bci,
                                     p3d_vj0, p3d_vk0, p3d_vi0, ifunc, n, jdw, kdw, idw,
                                     nplots, jdimg, kdimg, idimg, nblcg, jsg, ksg, isg,
                                     jeg, keg, ieg, ninter, iindex, intmax, nsub1,
                                     maxxe, nblk, nbli, limblk, isva, nblon, mxbli,
                                     thetay, maxbl, maxgr, myid, myhost, mycomm,
                                     mblk2nd, inpl3d, nblock, nblkpt, p3d_xv,
                                     p3d_sj, p3d_sk, p3d_si, p3d_vol, nset);
                }
                //
            }
            //
        } else if (inpl3d(n, 2) == 1 || inpl3d(n, 2) > 2) {

            //
            // cell center or face center data
            //
            if (mblk2nd(nbl) == myid || myid == myhost) {
                //
                lead_ns::lead(nbl, lw, lw2, maxbl);
                //
                i2 = std::min(idim - 1, i2);
                j2 = std::min(jdim - 1, j2);
                k2 = std::min(kdim - 1, k2);
                i1 = std::min(idim - 1, i1);
                j1 = std::min(jdim - 1, j1);
                k1 = std::min(kdim - 1, k1);
                jdw = (j2 - j1) / j3 + 1;
                kdw = (k2 - k1) / k3 + 1;
                idw = (i2 - i1) / i3 + 1;
                //
                // check storage availability
                //
                ixwk = 1;
                ibwk = ixwk + jdw * kdw * idw * 5;
                ixgk = ibwk + jdim * kdim * idim;
                nroom = nwork - (ixgk + jdw * kdw * idw * 4 - 1);
                if (nroom < 0) {
                    if (myid == myhost) {
                        FILE* f11 = fortran_get_unit(11);
                        fprintf(f11, " not enough memory for plot3c\n");
                        fprintf(f11, " have, need = %12d%12d\n", nwork, nwork - nroom);
                        fprintf(f11, " not writing out plot3d files\n");
                    }
                } else {
                    FortranArray4DRef<double> p3c_q  (&w(lq),   jdim, kdim, idim, 5);
                    FortranArray4DRef<double> p3c_qi0 (&w(lqi0), jdim, kdim, 5, 4);
                    FortranArray4DRef<double> p3c_qj0 (&w(lqj0), kdim, idim - 1, 5, 4);
                    FortranArray4DRef<double> p3c_qk0 (&w(lqk0), jdim, idim - 1, 5, 4);
                    FortranArray3DRef<double> p3c_x   (&w(lx),   jdim, kdim, idim);
                    FortranArray3DRef<double> p3c_y   (&w(ly),   jdim, kdim, idim);
                    FortranArray3DRef<double> p3c_z   (&w(lz),   jdim, kdim, idim);
                    FortranArray4DRef<double> p3c_xw  (&wk(ixwk), jdw, kdw, idw, 5);
                    FortranArray3DRef<double> p3c_bk2 (&wk(ibwk), jdim, kdim, idim);
                    FortranArray3DRef<double> p3c_blk (&w(lblk), jdim, kdim, idim);
                    FortranArray4DRef<double> p3c_xg  (&wk(ixgk), jdw, kdw, idw, 4);
                    FortranArray3DRef<double> p3c_vis (&w(lvis), jdim, kdim, idim);
                    FortranArray4DRef<double> p3c_vi0 (&w(lvi0), jdim, kdim, 1, 4);
                    FortranArray4DRef<double> p3c_vj0 (&w(lvj0), kdim, idim - 1, 1, 4);
                    FortranArray4DRef<double> p3c_vk0 (&w(lvk0), jdim, idim - 1, 1, 4);
                    FortranArray3DRef<double> p3c_smin(&w(lsnk0), jdim - 1, kdim - 1, idim - 1);
                    plot3c_ns::plot3c(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                                     p3c_q, p3c_qi0, p3c_qj0, p3c_qk0, p3c_x, p3c_y,
                                     p3c_z, p3c_xw, p3c_bk2, p3c_blk, p3c_xg,
                                     iflag, p3c_vis, p3c_vi0, p3c_vj0, p3c_vk0,
                                     iovrlp(nbl), nbl, nmap, p3c_smin, ifunc, n,
                                     jdw, kdw, idw,
                                     nplots, jdimg, kdimg, idimg, nblcg, jsg, ksg, isg,
                                     jeg, keg, ieg, ninter, iindex, intmax, nsub1,
                                     maxxe, nblk, nbli, limblk, isva, nblon, mxbli,
                                     thetay, maxbl, maxgr, myid, myhost, mycomm,
                                     mblk2nd, inpl3d, nblock, nblkpt, ip3dsurf);
                }
            }
            //
        } else {
            //
            // cell center turbulence data (plot3d q file format)
            //
            if (mblk2nd(nbl) == myid || myid == myhost) {
                //
                if (ivmx > 1) {
                    //
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    //
                    i2 = std::min(idim - 1, i2);
                    j2 = std::min(jdim - 1, j2);
                    k2 = std::min(kdim - 1, k2);
                    i1 = std::min(idim - 1, i1);
                    j1 = std::min(jdim - 1, j1);
                    k1 = std::min(kdim - 1, k1);
                    jdw = (j2 - j1) / j3 + 1;
                    kdw = (k2 - k1) / k3 + 1;
                    idw = (i2 - i1) / i3 + 1;
                    //
                    // check storage availability
                    //
                    if (ivmx == 8 || ivmx == 9 || ivmx >= 11 || lowmem_ux == 0) {
                        // have permanent storage for ux, starting at lux
                        ixwk  = 1;
                        ifuncdim = std::max(5, ifunct);
                        ibwk  = ixwk + jdw * kdw * idw * ifuncdim;
                        ixgk  = ibwk + jdim * kdim * idim;
                        ibwk3 = ixgk + jdw * kdw * idw * 4;
                        nroom = nwork - (ibwk3 + jdim * kdim * 9);
                        if (nroom < 0) {
                            if (myid == myhost) {
                                FILE* f11 = fortran_get_unit(11);
                                fprintf(f11, " not enough memory for plot3t\n");
                                fprintf(f11, " have, need = %12d%12d\n", nwork, nwork - nroom);
                                fprintf(f11, " not writing out plot3d files\n");
                            }
                        } else {
                            FortranArray4DRef<double> p3t_q   (&w(lq),   jdim, kdim, idim, 5);
                            FortranArray3DRef<double> p3t_x   (&w(lx),   jdim, kdim, idim);
                            FortranArray3DRef<double> p3t_y   (&w(ly),   jdim, kdim, idim);
                            FortranArray3DRef<double> p3t_z   (&w(lz),   jdim, kdim, idim);
                            FortranArray4DRef<double> p3t_xw  (&wk(ixwk), jdw, kdw, idw, ifuncdim);
                            FortranArray3DRef<double> p3t_bk2 (&wk(ibwk), jdim, kdim, idim);
                            FortranArray3DRef<double> p3t_blk (&w(lblk), jdim, kdim, idim);
                            FortranArray4DRef<double> p3t_xg  (&wk(ixgk), jdw, kdw, idw, 4);
                            FortranArray3DRef<double> p3t_vis (&w(lvis), jdim, kdim, idim);
                            FortranArray4DRef<double> p3t_sj  (&w(lsj),  jdim, kdim, idim - 1, 5);
                            FortranArray4DRef<double> p3t_sk  (&w(lsk),  jdim, kdim, idim - 1, 5);
                            FortranArray4DRef<double> p3t_si  (&w(lsi),  jdim, kdim, idim, 5);
                            FortranArray3DRef<double> p3t_smin(&w(lsnk0), jdim - 1, kdim - 1, idim - 1);
                            FortranArray4DRef<double> p3t_ux  (&w(lux),  jdim - 1, kdim - 1, idim - 1, 9);
                            FortranArray4DRef<double> p3t_turre(&w(lxib), jdim, kdim, idim, nummem);
                            FortranArray3DRef<double> p3t_vol (&w(lvol), jdim, kdim, idim - 1);
                            FortranArray4DRef<double> p3t_qj0 (&w(lqj0), kdim, idim - 1, 5, 4);
                            FortranArray4DRef<double> p3t_qk0 (&w(lqk0), jdim, idim - 1, 5, 4);
                            FortranArray4DRef<double> p3t_qi0 (&w(lqi0), jdim, kdim, 5, 4);
                            FortranArray3DRef<double> p3t_bcj (&w(lbcj), kdim, idim - 1, 2);
                            FortranArray3DRef<double> p3t_bck (&w(lbck), jdim, idim - 1, 2);
                            FortranArray3DRef<double> p3t_bci (&w(lbci), jdim, kdim, 2);
                            FortranArray3DRef<double> p3t_wt  (&wk(ibwk3), jdim, kdim, 9);
                            FortranArray3DRef<double> p3t_cmuv(&w(lcmuv), jdim - 1, kdim - 1, idim - 1);
                            FortranArray3DRef<double> p3t_volj0(&w(lvolj0), kdim, idim - 1, 4);
                            FortranArray3DRef<double> p3t_volk0(&w(lvolk0), jdim, idim - 1, 4);
                            FortranArray3DRef<double> p3t_voli0(&w(lvoli0), jdim, kdim, 4);
                            plot3t_ns::plot3t(jdim, kdim, idim, i1, i2, i3, j1, j2, j3,
                                             k1, k2, k3, p3t_q, p3t_x, p3t_y, p3t_z,
                                             p3t_xw, p3t_bk2, p3t_blk, p3t_xg,
                                             p3t_vis, iovrlp(nbl), nbl, nmap, p3t_sj,
                                             p3t_sk, p3t_si, p3t_smin, p3t_ux,
                                             p3t_turre, p3t_vol, p3t_qj0, p3t_qk0,
                                             p3t_qi0, p3t_bcj, p3t_bck, p3t_bci,
                                             p3t_wt, p3t_cmuv, jdw, kdw, idw, nplots,
                                             jdimg, kdimg, idimg, nblcg, jsg, ksg, isg,
                                             ieg, jeg, keg, ninter, iindex, intmax, nsub1,
                                             maxxe, nblk, nbli, limblk, isva, nblon,
                                             mxbli, thetay, maxbl, maxgr, myid, myhost,
                                             mycomm, mblk2nd, inpl3d, nblock, nblkpt,
                                             p3t_volj0, p3t_volk0, p3t_voli0, vormax,
                                             ivmax, jvmax, kvmax, nummem, ifuncdim);
                        }
                    } else {

                        // need to grab ux storage from temporary work array,
                        // starting at location ibwk2
                        ixwk  = 1;
                        ifuncdim = std::max(5, ifunct);
                        ibwk  = ixwk + jdw * kdw * idw * ifuncdim;
                        ixgk  = ibwk + jdim * kdim * idim;
                        ibwk2 = ixgk + jdw * kdw * idw * 4;
                        ibwk3 = ibwk2 + (jdim - 1) * (kdim - 1) * (idim - 1) * 9;
                        nroom = nwork - (ibwk3 + jdim * kdim * 9);
                        if (nroom < 0) {
                            if (myid == myhost) {
                                FILE* f11 = fortran_get_unit(11);
                                fprintf(f11, " not enough memory for plot3t\n");
                                fprintf(f11, " have, need = %12d%12d\n", nwork, nwork - nroom);
                                fprintf(f11, " not writing out plot3d files\n");
                            }
                        } else {
                            FortranArray4DRef<double> p3t_q   (&w(lq),   jdim, kdim, idim, 5);
                            FortranArray3DRef<double> p3t_x   (&w(lx),   jdim, kdim, idim);
                            FortranArray3DRef<double> p3t_y   (&w(ly),   jdim, kdim, idim);
                            FortranArray3DRef<double> p3t_z   (&w(lz),   jdim, kdim, idim);
                            FortranArray4DRef<double> p3t_xw  (&wk(ixwk), jdw, kdw, idw, ifuncdim);
                            FortranArray3DRef<double> p3t_bk2 (&wk(ibwk), jdim, kdim, idim);
                            FortranArray3DRef<double> p3t_blk (&w(lblk), jdim, kdim, idim);
                            FortranArray4DRef<double> p3t_xg  (&wk(ixgk), jdw, kdw, idw, 4);
                            FortranArray3DRef<double> p3t_vis (&w(lvis), jdim, kdim, idim);
                            FortranArray4DRef<double> p3t_sj  (&w(lsj),  jdim, kdim, idim - 1, 5);
                            FortranArray4DRef<double> p3t_sk  (&w(lsk),  jdim, kdim, idim - 1, 5);
                            FortranArray4DRef<double> p3t_si  (&w(lsi),  jdim, kdim, idim, 5);
                            FortranArray3DRef<double> p3t_smin(&w(lsnk0), jdim - 1, kdim - 1, idim - 1);
                            FortranArray4DRef<double> p3t_ux  (&wk(ibwk2), jdim - 1, kdim - 1, idim - 1, 9);
                            FortranArray4DRef<double> p3t_turre(&w(lxib), jdim, kdim, idim, nummem);
                            FortranArray3DRef<double> p3t_vol (&w(lvol), jdim, kdim, idim - 1);
                            FortranArray4DRef<double> p3t_qj0 (&w(lqj0), kdim, idim - 1, 5, 4);
                            FortranArray4DRef<double> p3t_qk0 (&w(lqk0), jdim, idim - 1, 5, 4);
                            FortranArray4DRef<double> p3t_qi0 (&w(lqi0), jdim, kdim, 5, 4);
                            FortranArray3DRef<double> p3t_bcj (&w(lbcj), kdim, idim - 1, 2);
                            FortranArray3DRef<double> p3t_bck (&w(lbck), jdim, idim - 1, 2);
                            FortranArray3DRef<double> p3t_bci (&w(lbci), jdim, kdim, 2);
                            FortranArray3DRef<double> p3t_wt  (&wk(ibwk3), jdim, kdim, 9);
                            FortranArray3DRef<double> p3t_cmuv(&w(lcmuv), jdim - 1, kdim - 1, idim - 1);
                            FortranArray3DRef<double> p3t_volj0(&w(lvolj0), kdim, idim - 1, 4);
                            FortranArray3DRef<double> p3t_volk0(&w(lvolk0), jdim, idim - 1, 4);
                            FortranArray3DRef<double> p3t_voli0(&w(lvoli0), jdim, kdim, 4);
                            plot3t_ns::plot3t(jdim, kdim, idim, i1, i2, i3, j1, j2, j3,
                                             k1, k2, k3, p3t_q, p3t_x, p3t_y, p3t_z,
                                             p3t_xw, p3t_bk2, p3t_blk, p3t_xg,
                                             p3t_vis, iovrlp(nbl), nbl, nmap, p3t_sj,
                                             p3t_sk, p3t_si, p3t_smin, p3t_ux,
                                             p3t_turre, p3t_vol, p3t_qj0, p3t_qk0,
                                             p3t_qi0, p3t_bcj, p3t_bck, p3t_bci,
                                             p3t_wt, p3t_cmuv, jdw, kdw, idw, nplots,
                                             jdimg, kdimg, idimg, nblcg, jsg, ksg, isg,
                                             jeg, keg, ieg, ninter, iindex, intmax, nsub1,
                                             maxxe, nblk, nbli, limblk, isva, nblon,
                                             mxbli, thetay, maxbl, maxgr, myid, myhost,
                                             mycomm, mblk2nd, inpl3d, nblock, nblkpt,
                                             p3t_volj0, p3t_volk0, p3t_voli0, vormax,
                                             ivmax, jvmax, kvmax, nummem, ifuncdim);
                        }
                    }
                    //
                } // end if (ivmx > 1)
                //
            } // end if (mblk2nd(nbl)==myid || myid==myhost)
            //
        } // end if/else on inpl3d(n,2)

        // Special turbulence output for certain models
        if (ivmx == 6 || ivmx == 7 || ivmx == 15 || ivmx == 72) {
            lead_ns::lead(nbl, lw, lw2, maxbl);
            {
                FortranArray4DRef<double> q_v  (&w(lq),   jdim,   kdim,   idim,   5);
                FortranArray3DRef<double> x_v  (&w(lx),   jdim,   kdim,   idim);
                FortranArray3DRef<double> y_v  (&w(ly),   jdim,   kdim,   idim);
                FortranArray3DRef<double> z_v  (&w(lz),   jdim,   kdim,   idim);
                FortranArray3DRef<double> vis_v(&w(lvis), jdim,   kdim,   idim);
                FortranArray4DRef<double> ux_v (&w(lux),  jdim-1, kdim-1, idim-1, 9);
                FortranArray4DRef<double> xib_v(&w(lxib), jdim,   kdim,   idim,   nummem);
                FortranArray3DRef<double> snk_v(&w(lsnk0),jdim-1, kdim-1, idim-1);
                if (module_kwstm_ns::need_contplot)
                    module_contour_ns::cont_plot(
                        nframes, jdim, kdim, idim, nummem, nbl,
                        q_v, x_v, y_v, z_v, vis_v, ux_v, xib_v, snk_v);
                if (ivmx == 72) {
                    module_kwstm_ns::kws_plot(
                        nframes, jdim, kdim, idim, nummem, nbl,
                        q_v, x_v, y_v, z_v, vis_v, ux_v, xib_v, snk_v);
                }
                module_profileout_ns::profile_plot(jdim, kdim, idim, nummem, nbl,
                    q_v, x_v, y_v, z_v, vis_v, ux_v, xib_v, snk_v);
                module_profileout_ns::cfcp_plot(jdim, kdim, idim, nummem, nbl,
                    q_v, x_v, y_v, z_v, vis_v, ux_v, xib_v, snk_v);
            }
        }
        //
        if (iblnkfr == 0) {
            //
            // reset blank values at fringe points to 0.
            //
            if (mblk2nd(nbl) == myid && iovrlp(nbl) != 0) {
                double blnkval = 0.0;
                blnkfr_ns::blnkfr(nbl, iibg, kkbg, jjbg, ibpntsg, lbg, iitot,
                                   FortranArray3DRef<double>(&w(lblk), jdim, kdim, idim),
                                   jdim, kdim, idim, maxbl, blnkval);
            }
        }
        //
    } // end do n=1,nplot3d (loop 70)
    //
    if (myid == myhost) {
        if (lhdr > 0) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, " \n");
        }
    }
    //
    // for stationary grid cases, set icall1 flag to prevent output
    // of grid to plot3d file on subsequent calls (which occur only if
    // abs(movie) > 0).  For dynamic grid cases, grid is output every
    // time plot3d routine is is called.
    //
    if (iunst == 0) {
        icall1 = 1;
    } else {
        icall1 = 0;
    }
    //
    label_231:;

    //
    // Print solution data.
    //
    if (nprint > 0) {
        //
        if (myid == myhost) {
            if (lhdr > 0 && nplot3d <= 0) {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, " \n");
            }
        }
        //
        if (myid == myhost) {
            FILE* f17 = fortran_get_unit(17);
            fprintf(f17, "\n");
            for (int ii = 1; ii <= 20; ii++)
                fprintf(f17, "%4.4s", (char*)&title[ii-1]);
            fprintf(f17, "\n");
            fprintf(f17, "      Mach     alpha      beta      ReUe   Tinf,dR      time\n");
            alphaw = (float)radtodeg * alpha_f;
            betaw  = (float)radtodeg * beta_f;
            fprintf(f17, "%10.5f%10.5f%10.5f%10.3e%10.5f%10.5f\n",
                    (double)xmach, (double)alphaw, (double)betaw,
                    (double)reue, (double)tinf, (double)time_f);
        }
        //
        for (n = 1; n <= nprint; n++) {
            //
            nbl = inpr(n, 1);
            //
            if (nbl > nblock) {
                if (myid == myhost) {
                    FILE* f11 = fortran_get_unit(11);
                    fprintf(f11, " Block%3d does not exist.  No output printed.\n", nbl);
                }
                continue; // go to 80
            }
            //
            i1 = inpr(n, 3);
            i2 = inpr(n, 4);
            i3 = inpr(n, 5);
            j1 = inpr(n, 6);
            j2 = inpr(n, 7);
            j3 = inpr(n, 8);
            k1 = inpr(n, 9);
            k2 = inpr(n, 10);
            k3 = inpr(n, 11);
            //
            int iflag = 2;
            //
            if (inpr(n, 2) == 0) {
                //
                // grid point data
                //
                if (mblk2nd(nbl) == myid || myid == myhost) {
                    //
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    //
                    if (myid == myhost) {
                        FILE* f17 = fortran_get_unit(17);
                        fprintf(f17, "\n\n BLOCK%4d  (GRID%4d)     IDIM,JDIM,KDIM=%5d%5d%5d\n",
                                nbl, igridg(nbl), idim, jdim, kdim);
                        fprintf(f17, " NOTE: endpts may not be reliable\n");
                    }
                    //
                    jdw = (j2 - j1) / j3 + 1;
                    kdw = (k2 - k1) / k3 + 1;
                    idw = (i2 - i1) / i3 + 1;
                    //
                    // check storage availability
                    //
                    ixwk = 1;
                    nset = 8;
                    ibwk = ixwk + jdw * kdw * idw * nset;
                    ixgk = ibwk + jdim * kdim * idim * 2;
                    ixvk = ixgk + jdw * kdw * idw * 4;
                    nroom = nwork - (ixvk + jdw * kdw * idw * 5 - 1);
                    if (nroom < 0) {
                        if (myid == myhost) {
                            FILE* f11 = fortran_get_unit(11);
                            fprintf(f11, " not enough memory for plot3d\n");
                            fprintf(f11, " have, need = %12d%12d\n", nwork, nwork - nroom);
                            fprintf(f11, " not writing out prout file\n");
                        }
                    } else {
                        FortranArray4DRef<double> p3d_q  (&w(lq),   jdim, kdim, idim, 5);
                        FortranArray4DRef<double> p3d_qi0 (&w(lqi0), jdim, kdim, 5, 4);
                        FortranArray4DRef<double> p3d_qj0 (&w(lqj0), kdim, idim - 1, 5, 4);
                        FortranArray4DRef<double> p3d_qk0 (&w(lqk0), jdim, idim - 1, 5, 4);
                        FortranArray3DRef<double> p3d_x   (&w(lx),   jdim, kdim, idim);
                        FortranArray3DRef<double> p3d_y   (&w(ly),   jdim, kdim, idim);
                        FortranArray3DRef<double> p3d_z   (&w(lz),   jdim, kdim, idim);
                        FortranArray4DRef<double> p3d_xw  (&wk(ixwk), jdw, kdw, idw, nset);
                        FortranArray4DRef<double> p3d_bk2 (&wk(ibwk), jdim, kdim, idim, 2);
                        FortranArray3DRef<double> p3d_blk (&w(lblk), jdim, kdim, idim);
                        FortranArray4DRef<double> p3d_xg  (&wk(ixgk), jdw, kdw, idw, 4);
                        FortranArray3DRef<double> p3d_vis (&w(lvis), jdim, kdim, idim);
                        FortranArray3DRef<double> p3d_bcj (&w(lbcj), kdim, idim - 1, 2);
                        FortranArray3DRef<double> p3d_bck (&w(lbck), jdim, idim - 1, 2);
                        FortranArray3DRef<double> p3d_bci (&w(lbci), jdim, kdim, 2);
                        FortranArray4DRef<double> p3d_vj0 (&w(lvj0), kdim, idim - 1, 1, 4);
                        FortranArray4DRef<double> p3d_vk0 (&w(lvk0), jdim, idim - 1, 1, 4);
                        FortranArray4DRef<double> p3d_vi0 (&w(lvi0), jdim, kdim, 1, 4);
                        FortranArray4DRef<double> p3d_xv  (&wk(ixvk), jdw, kdw, idw, 5);
                        FortranArray4DRef<double> p3d_sj  (&w(lsj),  jdim, kdim, idim - 1, 5);
                        FortranArray4DRef<double> p3d_sk  (&w(lsk),  jdim, kdim, idim - 1, 5);
                        FortranArray4DRef<double> p3d_si  (&w(lsi),  jdim, kdim, idim, 5);
                        FortranArray3DRef<double> p3d_vol (&w(lvol), jdim, kdim, idim - 1);
                        plot3d_ns::plot3d(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                                         p3d_q, p3d_qi0, p3d_qj0, p3d_qk0, p3d_x,
                                         p3d_y, p3d_z, p3d_xw, p3d_bk2,
                                         p3d_blk, p3d_xg, iflag, p3d_vis, iovrlp(nbl),
                                         nbl, nmap, p3d_bcj, p3d_bck, p3d_bci,
                                         p3d_vj0, p3d_vk0, p3d_vi0, ifunc, n, jdw, kdw, idw,
                                         nplots, jdimg, kdimg, idimg, nblcg, jsg, ksg, isg,
                                         jeg, keg, ieg, ninter, iindex, intmax, nsub1,
                                         maxxe, nblk, nbli, limblk, isva, nblon, mxbli,
                                         thetay, maxbl, maxgr, myid, myhost, mycomm,
                                         mblk2nd, inpl3d, nblock, nblkpt, p3d_xv,
                                         p3d_sj, p3d_sk, p3d_si, p3d_vol, nset);
                    }
                    //
                }
                //
            } else {
                //
                // cell center data
                //
                if (mblk2nd(nbl) == myid || myid == myhost) {
                    //
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    //
                    if (myid == myhost) {
                        FILE* f17 = fortran_get_unit(17);
                        fprintf(f17, "\n\n BLOCK%4d  (GRID%4d)     IDIM,JDIM,KDIM=%5d%5d%5d\n",
                                nbl, igridg(nbl), idim, jdim, kdim);
                    }
                    //
                    i2 = std::min(idim - 1, i2);
                    j2 = std::min(jdim - 1, j2);
                    k2 = std::min(kdim - 1, k2);
                    i1 = std::min(idim - 1, i1);
                    j1 = std::min(jdim - 1, j1);
                    k1 = std::min(kdim - 1, k1);
                    jdw = (j2 - j1) / j3 + 1;
                    kdw = (k2 - k1) / k3 + 1;
                    idw = (i2 - i1) / i3 + 1;
                    //
                    // check storage availability
                    //
                    ixwk = 1;
                    ibwk = ixwk + jdw * kdw * idw * 5;
                    ixgk = ibwk + jdim * kdim * idim;
                    nroom = nwork - (ixgk + jdw * kdw * idw * 4 - 1);
                    if (nroom < 0) {
                        if (myid == myhost) {
                            FILE* f11 = fortran_get_unit(11);
                            fprintf(f11, " not enough memory for plot3c\n");
                            fprintf(f11, " have, need = %12d%12d\n", nwork, nwork - nroom);
                            fprintf(f11, " not writing out plot3d files\n");
                        }
                    } else {
                        FortranArray4DRef<double> p3c_q  (&w(lq),   jdim, kdim, idim, 5);
                        FortranArray4DRef<double> p3c_qi0 (&w(lqi0), jdim, kdim, 5, 4);
                        FortranArray4DRef<double> p3c_qj0 (&w(lqj0), kdim, idim - 1, 5, 4);
                        FortranArray4DRef<double> p3c_qk0 (&w(lqk0), jdim, idim - 1, 5, 4);
                        FortranArray3DRef<double> p3c_x   (&w(lx),   jdim, kdim, idim);
                        FortranArray3DRef<double> p3c_y   (&w(ly),   jdim, kdim, idim);
                        FortranArray3DRef<double> p3c_z   (&w(lz),   jdim, kdim, idim);
                        FortranArray4DRef<double> p3c_xw  (&wk(ixwk), jdw, kdw, idw, 5);
                        FortranArray3DRef<double> p3c_bk2 (&wk(ibwk), jdim, kdim, idim);
                        FortranArray3DRef<double> p3c_blk (&w(lblk), jdim, kdim, idim);
                        FortranArray4DRef<double> p3c_xg  (&wk(ixgk), jdw, kdw, idw, 4);
                        FortranArray3DRef<double> p3c_vis (&w(lvis), jdim, kdim, idim);
                        FortranArray4DRef<double> p3c_vi0 (&w(lvi0), jdim, kdim, 1, 4);
                        FortranArray4DRef<double> p3c_vj0 (&w(lvj0), kdim, idim - 1, 1, 4);
                        FortranArray4DRef<double> p3c_vk0 (&w(lvk0), jdim, idim - 1, 1, 4);
                        FortranArray3DRef<double> p3c_smin(&w(lsnk0), jdim - 1, kdim - 1, idim - 1);
                        plot3c_ns::plot3c(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                                         p3c_q, p3c_qi0, p3c_qj0, p3c_qk0, p3c_x, p3c_y,
                                         p3c_z, p3c_xw, p3c_bk2, p3c_blk, p3c_xg,
                                         iflag, p3c_vis, p3c_vi0, p3c_vj0, p3c_vk0,
                                         iovrlp(nbl), nbl, nmap, p3c_smin, ifunc, n,
                                         jdw, kdw, idw,
                                         nplots, jdimg, kdimg, idimg, nblcg, jsg, ksg, isg,
                                         jeg, keg, ieg, ninter, iindex, intmax, nsub1,
                                         maxxe, nblk, nbli, limblk, isva, nblon, mxbli,
                                         thetay, maxbl, maxgr, myid, myhost, mycomm,
                                         mblk2nd, inpl3d, nblock, nblkpt, iprnsurf);
                    }
                }
                //
            } // end if/else on inpr(n,2)
            //
        } // end do n=1,nprint (loop 80)
        //
        if (myid == myhost) {
            if (lhdr > 0) {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, " \n");
            }
        }
        //
    } // end if (nprint > 0)
    return;
} // end qout

} // namespace qout_ns
