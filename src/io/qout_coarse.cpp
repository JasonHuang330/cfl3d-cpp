// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "qout_coarse.h"
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
#include "plot3d_coarse.h"
#include "plot3c.h"
#include "plot3t.h"
#include "termn8.h"
#include "qout.h"
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <cmath>

namespace qout_coarse_ns {

// Thin wrappers for bc, bc_blkint, bc_embed, bc_patch, bc_period, qout
// (declared in header but implemented here as pass-throughs)

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
          int& mxbcfil, int& iprnsurf, int& nummem)
{
    qout_ns::qout(iseq, lw, lw2, w, mgwk, wk, nwork, nplots, iovrlp, iibg, kkbg,
                  jjbg, ibcg, lbg, ibpntsg, qb, lwdat, nbci0, nbcj0, nbck0,
                  nbcidim, nbcjdim, nbckdim, jbcinfo, kbcinfo, ibcinfo,
                  bcfilei, bcfilej, bcfilek, itrans, irotat, idefrm, nblock,
                  levelg, igridg, iviscg, jdimg, kdimg, idimg, nblg, clw,
                  ncycmax, nplot3d, inpl3d, ip3dsurf, nprint, inpr, iadvance,
                  mycomm, myid, myhost, mblk2nd, nou, bou, nbuf, ibufdim,
                  maxbl, maxgr, maxseg, iitot, jsg, ksg, isg, jeg, keg, ieg,
                  ninter, windex, iindex, nblkpt, intmax, nsub1, maxxe,
                  nblk, nbli, limblk, isva, nblon, mxbli, thetay, ip3ddim,
                  nmap, iwk, iwork, xorig, yorig, zorig, period_miss, geom_miss,
                  epsc0, epsrot, isav_blk, isav_pat, isav_pat_b, isav_emb,
                  isav_prd, lbcprd, lbcemb, dthetxx, dthetyy, dthetzz,
                  nblcg, lfgm, istat2_bl, istat2_pa, istat2_pe, istat2_em,
                  istat_size, vormax, ivmax, jvmax, kvmax, bcfiles, mxbcfil,
                  iprnsurf, nummem);
}

void qout_coarse(int& iseq, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
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
                 FortranArray1DRef<int> iadvance, int& mycomm, int& myid,
                 int& myhost, FortranArray1DRef<int> mblk2nd,
                 FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
                 int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg,
                 int& iitot, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg,
                 FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg,
                 FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg,
                 int& ninter, FortranArray2DRef<double> windex,
                 FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt,
                 int& intmax, int& nsub1, int& maxxe, FortranArray2DRef<int> nblk,
                 int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva,
                 FortranArray1DRef<int> nblon, int& mxbli,
                 FortranArray1DRef<double> thetay, FortranArray2DRef<int> ip3ddim,
                 FortranArray1DRef<int> nmap, FortranArray1DRef<int> iwk,
                 int& iwork, FortranArray1DRef<double> xorig,
                 FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig,
                 FortranArray1DRef<double> period_miss,
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
                 int& iprnsurf, int& nt, int& mov_coarse, int& iinc_coarse,
                 int& jinc_coarse, int& kinc_coarse, int& nummem)
{
    // COMMON block references
    int& ibin      = cmn_bin.ibin;
    int& iblnkfr   = cmn_bin.iblnkfr;
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
    int& lxib      = cmn_ginfo.lxib;
    int& lcmuv     = cmn_ginfo.lcmuv;
    int& lvolj0    = cmn_ginfo.lvolj0;
    int& lvolk0    = cmn_ginfo.lvolk0;
    int& lvoli0    = cmn_ginfo.lvoli0;
    int& lux       = cmn_ginfo.lux;
    int& iipv      = cmn_info.iipv;
    int& ntt       = cmn_info.ntt;
    int& mseq      = cmn_info.mseq;
    int* levelt    = cmn_info.levelt;  // 1-based array, use levelt[iseq-1]
    int& ivmx      = cmn_maxiv.ivmx;
    int& level     = cmn_mgrd.level;
    int& lglobal   = cmn_mgrd.lglobal;
    int& i2d       = cmn_twod.i2d;
    int& iunst     = cmn_unst.iunst;
    int& lhdr      = cmn_moov.lhdr;
    int& icallcrs  = cmn_moovcrs2d.icallcrs;
    float& radtodeg = cmn_conversion.radtodeg;
    float& xmach   = cmn_info.xmach;
    float& alpha   = cmn_info.alpha;
    float& beta    = cmn_info.beta;
    float& reue    = cmn_reyue.reue;
    float& tinf    = cmn_reyue.tinf;
    float& time    = cmn_unst.time;
    float* title   = cmn_info.title;  // 1-based, use title[i-1]

    // Local variables
    int level_sav, nttuse, lres, nsafe, mneed;
    int iwk1, iwk2, iwk3, iwk4, iwk5, iwk6, iwork1;
    int iii, nbl, ncount, np3d, ifunc;
    int i1, i2, i3, j1, j2, j3, k1, k2, k3;
    int jdw, kdw, idw, iflag;
    int ixwk, ibwk, ixgk, ixvk, ibwk2, ibwk3, nroom, nset;
    int n, nnn, m;
    double clwuse;
    int termn8_err = -1;
    int int_updt;



    // if (1.gt.0) then  -- always executes
    {
        // update all bc's for level lglobal (and above, if embedded)
        level_sav = level;
        lglobal = lfgm - (mseq - iseq);

        for (level = lglobal; level <= levelt[iseq-1]; level++) {

            if (iipv > 0) {
                nttuse = std::max(ntt-1, 1);
                clwuse = clw(nttuse);
            }

            for (nbl = 1; nbl <= nblock; nbl++) {
                // need to call even for blocks not advanced
                if (level == levelg(nbl)) {
                    if (mblk2nd(nbl) == myid) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        {
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
            }

            // update periodic boundary conditions
            lres   = 1;
            nsafe  = nwork - lres + 1;
            mneed  = lbcprd * 5;
            iwk1   = 1;
            iwk2   = iwk1 + mneed;
            iwk3   = iwk2 + mneed;
            iwk4   = iwk3 + mneed;
            iwk5   = iwk4 + mneed*2;
            iwk6   = iwk5 + mneed;
            iwork1 = iwork - iwk6;
            if (iwork1 < 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "stopping...not enough integer work space for subroutine bc_period");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                int ntime1 = 1;
                FortranArray1DRef<int> iwk_1(&iwk(iwk1), mneed);
                FortranArray1DRef<int> iwk_2(&iwk(iwk2), mneed);
                FortranArray1DRef<int> iwk_3(&iwk(iwk3), mneed);
                FortranArray2DRef<int> iwk_4(&iwk(iwk4), lbcprd, 10);
                FortranArray1DRef<int> iwk_5(&iwk(iwk5), mneed);
                FortranArray1DRef<double> wk_lres(&wk(lres), nsafe);
                bc_period_ns::bc_period(ntime1, nbl, lw, lw2, w, mgwk, wk_lres, nsafe,
                                        maxbl, maxgr, maxseg, iadvance,
                                        bcfilei, bcfilej, bcfilek, lwdat,
                                        xorig, yorig, zorig, jdimg, kdimg,
                                        idimg, lbcprd, isav_prd,
                                        period_miss, epsrot, iwk_1, iwk_2,
                                        iwk_3, iwk_4, iwk_5, myid, myhost,
                                        mycomm, mblk2nd, nou, bou, nbuf, ibufdim,
                                        istat2_pe, istat_size, bcfiles, mxbcfil, nummem);
            }

            // update embedded-grid boundary conditions
            lres   = 1;
            nsafe  = nwork - lres + 1;
            mneed  = lbcemb * 3;
            iwk1   = 1;
            iwk2   = iwk1 + mneed;
            iwk3   = iwk2 + mneed;
            iwk4   = iwk3 + mneed;
            iwk5   = iwk4 + mneed*2;
            iwk6   = iwk5 + mneed;
            iwork1 = iwork - iwk6;
            if (iwork1 < 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "stopping...not enough integer work space for subroutine bc_embed");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                int ntime1 = 1;
                FortranArray1DRef<int> iwk_1(&iwk(iwk1), mneed);
                FortranArray1DRef<int> iwk_2(&iwk(iwk2), mneed);
                FortranArray1DRef<int> iwk_3(&iwk(iwk3), mneed);
                FortranArray2DRef<int> iwk_4(&iwk(iwk4), lbcemb, 6);
                FortranArray1DRef<int> iwk_5(&iwk(iwk5), mneed);
                FortranArray1DRef<double> wk_lres(&wk(lres), nsafe);
                bc_embed_ns::bc_embed(ntime1, nbl, lw, lw2, w, mgwk, wk_lres, nsafe,
                                      maxbl, maxgr, lbcemb, iadvance, idimg, jdimg,
                                      kdimg, isav_emb, iwk_1,
                                      iwk_2, iwk_3, iwk_4, iwk_5,
                                      myid, myhost, mycomm, mblk2nd, nou, bou, nbuf,
                                      ibufdim, iviscg, istat2_em, istat_size, nummem);
            }



            // update 1-1 block boundary conditions
            lres   = 1;
            nsafe  = nwork - lres + 1;
            mneed  = mxbli * 5;
            iwk1   = 1;
            iwk2   = iwk1 + mneed;
            iwk3   = iwk2 + mneed;
            iwk4   = iwk3 + mneed;
            iwk5   = iwk4 + mneed*2;
            iwk6   = iwk5 + mneed;
            iwork1 = iwork - iwk6;
            if (iwork1 < 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "stopping...not enough integer work space for subroutine bc_blkint");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                int ntime1 = 1;
                FortranArray1DRef<int> iwk_1(&iwk(iwk1), mneed);
                FortranArray1DRef<int> iwk_2(&iwk(iwk2), mneed);
                FortranArray1DRef<int> iwk_3(&iwk(iwk3), mneed);
                FortranArray2DRef<int> iwk_4(&iwk(iwk4), mxbli, 10);
                FortranArray1DRef<int> iwk_5(&iwk(iwk5), mneed);
                FortranArray1DRef<double> wk_lres(&wk(lres), nsafe);
                bc_blkint_ns::bc_blkint(ntime1, nbl, lw, lw2, w, mgwk, wk_lres, nsafe,
                                        maxbl, maxgr, mxbli, iadvance, geom_miss, epsc0,
                                        nblk, nbli, limblk, isva, nblon, jdimg, kdimg,
                                        idimg, mblk2nd, isav_blk, iwk_1,
                                        iwk_2, iwk_3, iwk_4, iwk_5,
                                        nou, bou, nbuf, ibufdim, myid, myhost, mycomm,
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
            iwk5   = iwk4 + mneed*2;
            iwk6   = iwk5 + mneed*2;
            iwork1 = iwork - iwk6;
            if (iwork1 < 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "stopping...not enough integer work space for subroutine bc_patch");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "have, need = %12d%12d", iwork, iwk6);
                termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk6; iii++) {
                iwk(iii) = 0;
            }
            {
                int ntime1 = 1;
                FortranArray1DRef<int> iwk_1(&iwk(iwk1), mneed);
                FortranArray1DRef<int> iwk_2(&iwk(iwk2), mneed);
                FortranArray1DRef<int> iwk_3(&iwk(iwk3), mneed);
                FortranArray3DRef<int> iwk_4(&iwk(iwk4), intmax, nsub1, 6);
                FortranArray2DRef<int> iwk_5(&iwk(iwk5), mneed, 2);
                FortranArray1DRef<double> wk_lres(&wk(lres), nsafe);
                bc_patch_ns::bc_patch(ntime1, nbl, lw, lw2, w, mgwk, wk_lres, nsafe,
                                      maxbl, maxgr, intmax, nsub1, maxxe, iadvance,
                                      jdimg, kdimg, idimg, ninter, windex, iindex,
                                      nblkpt, dthetxx, dthetyy, dthetzz,
                                      isav_pat, isav_pat_b,
                                      iwk_1, iwk_2, iwk_3,
                                      iwk_4, iwk_5, myid, myhost, mycomm,
                                      mblk2nd, nou, bou, nbuf, ibufdim,
                                      istat2_pa, istat_size, nummem);
            }

            // update chimera boundary conditions
            // int_updt = 1 (commented-out loop that could set it to 0 is disabled)
            int_updt = 1;

            for (nbl = 1; nbl <= nblock; nbl++) {
                if (iadvance(nbl) >= 0) {
                    if (level == levelg(nbl)) {
                        if (mblk2nd(nbl) == myid) {
                            lead_ns::lead(nbl, lw, lw2, maxbl);
                            lres  = 1;
                            nsafe = nwork - lres + 1;
                            {
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
            }

            // call qface to install face-center values in qi0/qj0/qk0 arrays
            for (nbl = 1; nbl <= nblock; nbl++) {
                if (level == levelg(nbl)) {
                    if (mblk2nd(nbl) == myid) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        {
                            int ldim = 5;
                            FortranArray4DRef<double> q_arr(&w(lq), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> qj0_arr(&w(lqj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> qk0_arr(&w(lqk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> qi0_arr(&w(lqi0), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> bcj_arr(&w(lbcj), jdim, kdim, 2);
                            FortranArray3DRef<double> bck_arr(&w(lbck), jdim, kdim, 2);
                            FortranArray3DRef<double> bci_arr(&w(lbci), jdim, kdim, 2);
                            FortranArray3DRef<double> lblk_arr(&w(lblk), jdim, kdim, idim);
                            qface_ns::qface(jdim, kdim, idim, q_arr, qj0_arr,
                                            qk0_arr, qi0_arr, bcj_arr, bck_arr,
                                            bci_arr, lblk_arr, ldim);
                        }
                        if (ivmx >= 2) {
                            int ldim = 1;
                            FortranArray4DRef<double> vis_arr(&w(lvis), jdim, kdim, idim, 1);
                            FortranArray4DRef<double> vj0_arr(&w(lvj0), jdim, kdim, idim, 1);
                            FortranArray4DRef<double> vk0_arr(&w(lvk0), jdim, kdim, idim, 1);
                            FortranArray4DRef<double> vi0_arr(&w(lvi0), jdim, kdim, idim, 1);
                            FortranArray3DRef<double> bcj_arr(&w(lbcj), jdim, kdim, 2);
                            FortranArray3DRef<double> bck_arr(&w(lbck), jdim, kdim, 2);
                            FortranArray3DRef<double> bci_arr(&w(lbci), jdim, kdim, 2);
                            FortranArray3DRef<double> lblk_arr(&w(lblk), jdim, kdim, idim);
                            qface_ns::qface(jdim, kdim, idim, vis_arr, vj0_arr,
                                            vk0_arr, vi0_arr, bcj_arr, bck_arr,
                                            bci_arr, lblk_arr, ldim);
                        }
                    }
                }
            }

        } // end do level=lglobal,levelt(iseq)

        level = level_sav;

    } // end if (1.gt.0)



    // !if (nplot3d.le.0) go to 231  -- commented out, so always continue

    ncount = 0;
    np3d   = nblock / lglobal;
    ifunc  = 0;

    // !if (abs(inpl3d(1,2)).gt.2) ifunc = inpl3d(1,2)  -- commented out

    // do 60 n=1,np3d
    for (n = 1; n <= np3d; n++) {
        if (n == 1) {
            if (myid == myhost) {
                if (ibin == 0) {
                    if (icallcrs == 0) {
                        FILE* f101 = fortran_get_unit(101);
                        fprintf(f101, "%5d%5d%5d\n", np3d, 0, 0);
                    }
                    {
                        FILE* f92 = fortran_get_unit(92);
                        fprintf(f92, "%5d%5d%5d\n", np3d, 0, 0);
                    }
                } else {
                    if (icallcrs == 0) {
                        FILE* f101 = fortran_get_unit(101);
                        fwrite(&np3d, sizeof(int), 1, f101);
                    }
                    {
                        FILE* f92 = fortran_get_unit(92);
                        fwrite(&np3d, sizeof(int), 1, f92);
                    }
                }
            }
        }

        // nbl = lglobal*(n-1)+1
        nbl = lglobal * (n-1) + 1;

        if (nbl > nblock) {
            if (myid == myhost) {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, " Block%3d does not exist.  No plot3d output printed.\n", nbl);
            }
            continue; // go to 60
        }

        ncount++;
        lead_ns::lead(nbl, lw, lw2, maxbl);
        // !if incr is too large, allow 1 level of halving it any given block
        i1 = 1;                    // inpl3d(n,3)
        i2 = idim;                 // inpl3d(n,4)
        i3 = iinc_coarse;          // inpl3d(n,5)
        if ((iinc_coarse > 1) && (((idim-1) % iinc_coarse) != 0))
            i3 = iinc_coarse / 2;
        j1 = 1;                    // inpl3d(n,6)
        j2 = jdim;                 // inpl3d(n,7)
        j3 = jinc_coarse;          // inpl3d(n,8)
        if ((jinc_coarse > 1) && (((jdim-1) % jinc_coarse) != 0))
            j3 = jinc_coarse / 2;
        k1 = 1;                    // inpl3d(n,9)
        k2 = kdim;                 // inpl3d(n,10)
        k3 = kinc_coarse;          // inpl3d(n,11)
        if ((kinc_coarse > 1) && (((kdim-1) % kinc_coarse) != 0))
            k3 = kinc_coarse / 2;

        ip3ddim(1,ncount) = (i2-i1)/i3 + 1;
        ip3ddim(2,ncount) = (j2-j1)/j3 + 1;
        ip3ddim(3,ncount) = (k2-k1)/k3 + 1;
    } // end do 60

    if (myid == myhost) {
        if (ibin == 0) {
            if (i2d == 0) {
                if (icallcrs == 0) {
                    FILE* f101 = fortran_get_unit(101);
                    for (n = 1; n <= ncount; n++)
                        fprintf(f101, "%5d%5d%5d", ip3ddim(1,n), ip3ddim(2,n), ip3ddim(3,n));
                    fprintf(f101, "\n");
                }
                {
                    FILE* f92 = fortran_get_unit(92);
                    if (ifunc == 0) {
                        for (n = 1; n <= ncount; n++)
                            fprintf(f92, "%5d%5d%5d", ip3ddim(1,n), ip3ddim(2,n), ip3ddim(3,n));
                        fprintf(f92, "\n");
                    } else {
                        for (n = 1; n <= ncount; n++)
                            fprintf(f92, "%5d%5d%5d%5d", ip3ddim(1,n), ip3ddim(2,n), ip3ddim(3,n), 1);
                        fprintf(f92, "\n");
                    }
                }
            } else {
                // i2d != 0: write indices 2 and 3 only
                if (icallcrs == 0) {
                    FILE* f101 = fortran_get_unit(101);
                    for (n = 1; n <= ncount; n++)
                        fprintf(f101, "%5d%5d%5d", ip3ddim(2,n), ip3ddim(3,n), 0);
                    fprintf(f101, "\n");
                }
                {
                    FILE* f92 = fortran_get_unit(92);
                    if (ifunc == 0) {
                        for (n = 1; n <= ncount; n++)
                            fprintf(f92, "%5d%5d%5d", ip3ddim(2,n), ip3ddim(3,n), 0);
                        fprintf(f92, "\n");
                    } else {
                        for (n = 1; n <= ncount; n++)
                            fprintf(f92, "%5d%5d%5d%5d", ip3ddim(2,n), ip3ddim(3,n), 0, 1);
                        fprintf(f92, "\n");
                    }
                }
            }
        } else {
            // binary
            if (i2d == 0) {
                if (icallcrs == 0) {
                    FILE* f101 = fortran_get_unit(101);
                    for (n = 1; n <= ncount; n++) {
                        int v1 = ip3ddim(1,n), v2 = ip3ddim(2,n), v3 = ip3ddim(3,n);
                        fwrite(&v1, sizeof(int), 1, f101);
                        fwrite(&v2, sizeof(int), 1, f101);
                        fwrite(&v3, sizeof(int), 1, f101);
                    }
                }
                {
                    FILE* f92 = fortran_get_unit(92);
                    if (ifunc == 0) {
                        for (n = 1; n <= ncount; n++) {
                            int v1 = ip3ddim(1,n), v2 = ip3ddim(2,n), v3 = ip3ddim(3,n);
                            fwrite(&v1, sizeof(int), 1, f92);
                            fwrite(&v2, sizeof(int), 1, f92);
                            fwrite(&v3, sizeof(int), 1, f92);
                        }
                    } else {
                        for (n = 1; n <= ncount; n++) {
                            int v1 = ip3ddim(1,n), v2 = ip3ddim(2,n), v3 = ip3ddim(3,n), v4 = 1;
                            fwrite(&v1, sizeof(int), 1, f92);
                            fwrite(&v2, sizeof(int), 1, f92);
                            fwrite(&v3, sizeof(int), 1, f92);
                            fwrite(&v4, sizeof(int), 1, f92);
                        }
                    }
                }
            } else {
                // i2d != 0, binary
                if (icallcrs == 0) {
                    FILE* f101 = fortran_get_unit(101);
                    for (n = 1; n <= ncount; n++) {
                        int v2 = ip3ddim(2,n), v3 = ip3ddim(3,n);
                        fwrite(&v2, sizeof(int), 1, f101);
                        fwrite(&v3, sizeof(int), 1, f101);
                    }
                }
                {
                    FILE* f92 = fortran_get_unit(92);
                    if (ifunc == 0) {
                        for (n = 1; n <= ncount; n++) {
                            int v2 = ip3ddim(2,n), v3 = ip3ddim(3,n);
                            fwrite(&v2, sizeof(int), 1, f92);
                            fwrite(&v3, sizeof(int), 1, f92);
                        }
                    } else {
                        for (n = 1; n <= ncount; n++) {
                            int v2 = ip3ddim(2,n), v3 = ip3ddim(3,n), v4 = 1;
                            fwrite(&v2, sizeof(int), 1, f92);
                            fwrite(&v3, sizeof(int), 1, f92);
                            fwrite(&v4, sizeof(int), 1, f92);
                        }
                    }
                }
            }
        }
    }



    // correspondence between global block number and plot3d block number
    // stored in nmap(n) for n=1,nblock
    for (n = 1; n <= nblock; n++) {
        nmap(n) = 1;
        for (nnn = 1; nnn <= np3d; nnn++) {
            m = lglobal * (nnn-1) + 1;
            if (n == m) nmap(n) = nnn;
        }
    }

    if (myid == myhost) {
        if (lhdr > 0) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, " \n");
        }
    }

    // do 70 n=1,np3d
    for (n = 1; n <= np3d; n++) {

        nbl = lglobal * (n-1) + 1;

        if (nbl > nblock || nbl <= 0) continue; // go to 70

        if (iblnkfr == 0) {
            // temporarily set blank values at fringe points (not holes)
            // to 1 for plotting purposes
            lead_ns::lead(nbl, lw, lw2, maxbl);
            if (mblk2nd(nbl) == myid && iovrlp(nbl) != 0) {
                FortranArray3DRef<double> blank_arr(&w(lblk), jdim, kdim, idim);
                double blnkval1 = 1.0;
                blnkfr_ns::blnkfr(nbl, iibg, kkbg, jjbg, ibpntsg, lbg, iitot,
                                   blank_arr, jdim, kdim, idim, maxbl, blnkval1);
            }
        }

        ncount++;
        lead_ns::lead(nbl, lw, lw2, maxbl);
        i1 = 1;
        i2 = idim;
        i3 = iinc_coarse;
        if ((iinc_coarse > 1) && (((idim-1) % iinc_coarse) != 0))
            i3 = iinc_coarse / 2;
        j1 = 1;
        j2 = jdim;
        j3 = jinc_coarse;
        if ((jinc_coarse > 1) && (((jdim-1) % jinc_coarse) != 0))
            j3 = jinc_coarse / 2;
        k1 = 1;
        k2 = kdim;
        k3 = kinc_coarse;
        if ((kinc_coarse > 1) && (((kdim-1) % kinc_coarse) != 0))
            k3 = kinc_coarse / 2;

        iflag = 1;

        // if (0.le.0) then  -- always executes (grid point data)
        if (0 <= 0) {
            if (mblk2nd(nbl) == myid || myid == myhost) {

                lead_ns::lead(nbl, lw, lw2, maxbl);

                jdw = (j2-j1)/j3 + 1;
                kdw = (k2-k1)/k3 + 1;
                idw = (i2-i1)/i3 + 1;

                // check storage availability
                ixwk = 1;
                nset = 5;
                ibwk = ixwk + jdw*kdw*idw*nset;
                ixgk = ibwk + jdim*kdim*idim*2;
                ixvk = ixgk + jdw*kdw*idw*4;
                nroom = nwork - (ixvk + jdw*kdw*idw*5 - 1);
                if (nroom < 0) {
                    if (myid == myhost) {
                        FILE* f11 = fortran_get_unit(11);
                        fprintf(f11, " not enough memory for plot3d_coarse\n");
                        fprintf(f11, " have, need = %12d%12d\n", nwork, nwork-nroom);
                        fprintf(f11, " Aborting -not calling plot3d_coarse\n");
                        std::exit(0);
                    }
                } else {
                    FortranArray4DRef<double> q_arr(&w(lq), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qi0_arr(&w(lqi0), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qj0_arr(&w(lqj0), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qk0_arr(&w(lqk0), jdim, kdim, idim, 5);
                    FortranArray3DRef<double> x_arr(&w(lx), jdim, kdim, idim);
                    FortranArray3DRef<double> y_arr(&w(ly), jdim, kdim, idim);
                    FortranArray3DRef<double> z_arr(&w(lz), jdim, kdim, idim);
                    FortranArray4DRef<double> xw_arr(&wk(ixwk), jdw, kdw, idw, nset);
                    FortranArray4DRef<double> blank2_arr(&wk(ibwk), jdim, kdim, idim, 2);
                    FortranArray3DRef<double> blank_arr(&w(lblk), jdim, kdim, idim);
                    FortranArray4DRef<double> xg_arr(&wk(ixgk), jdw, kdw, idw, 4);
                    FortranArray3DRef<double> vis_arr(&w(lvis), jdim, kdim, idim);
                    FortranArray3DRef<double> bcj_arr(&w(lbcj), jdim, kdim, 2);
                    FortranArray3DRef<double> bck_arr(&w(lbck), jdim, kdim, 2);
                    FortranArray3DRef<double> bci_arr(&w(lbci), jdim, kdim, 2);
                    FortranArray4DRef<double> vj0_arr(&w(lvj0), jdim, kdim, idim, 1);
                    FortranArray4DRef<double> vk0_arr(&w(lvk0), jdim, kdim, idim, 1);
                    FortranArray4DRef<double> vi0_arr(&w(lvi0), jdim, kdim, idim, 1);
                    FortranArray4DRef<double> xv_arr(&wk(ixvk), jdw, kdw, idw, 5);
                    FortranArray4DRef<double> sj_arr(&w(lsj), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> sk_arr(&w(lsk), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> si_arr(&w(lsi), jdim, kdim, idim, 5);
                    FortranArray3DRef<double> vol_arr(&w(lvol), jdim, kdim, idim);
                    plot3d_coarse_ns::plot3d_coarse(jdim, kdim, idim, i1, i2, i3,
                                                    j1, j2, j3, k1, k2, k3,
                                                    q_arr, qi0_arr, qj0_arr, qk0_arr,
                                                    x_arr, y_arr, z_arr,
                                                    xw_arr, blank2_arr, blank_arr,
                                                    xg_arr, iflag, vis_arr,
                                                    iovrlp(nbl), nbl, nmap,
                                                    bcj_arr, bck_arr, bci_arr,
                                                    vj0_arr, vk0_arr, vi0_arr,
                                                    ifunc, n, jdw, kdw, idw,
                                                    nplots, jdimg, kdimg, idimg,
                                                    nblcg, jsg, ksg, isg,
                                                    jeg, keg, ieg, ninter, iindex,
                                                    intmax, nsub1, maxxe,
                                                    nblk, nbli, limblk, isva, nblon,
                                                    mxbli, thetay, maxbl, maxgr,
                                                    myid, myhost, mycomm, mblk2nd,
                                                    inpl3d, nblock, nblkpt,
                                                    xv_arr, sj_arr, sk_arr, si_arr,
                                                    vol_arr, nset);
                }
            }
        } else if (0 == 1 || 0 > 2) {


            // never executes: cell center or face center data
            if (mblk2nd(nbl) == myid || myid == myhost) {

                lead_ns::lead(nbl, lw, lw2, maxbl);

                i2 = std::min(idim-1, i2);
                j2 = std::min(jdim-1, j2);
                k2 = std::min(kdim-1, k2);
                i1 = std::min(idim-1, i1);
                j1 = std::min(jdim-1, j1);
                k1 = std::min(kdim-1, k1);
                jdw = (j2-j1)/j3 + 1;
                kdw = (k2-k1)/k3 + 1;
                idw = (i2-i1)/i3 + 1;

                ixwk = 1;
                ibwk = ixwk + jdw*kdw*idw*5;
                ixgk = ibwk + jdim*kdim*idim;
                nroom = nwork - (ixgk + jdw*kdw*idw*4 - 1);
                if (nroom < 0) {
                    if (myid == myhost) {
                        FILE* f11 = fortran_get_unit(11);
                        fprintf(f11, " not enough memory for plot3c\n");
                        fprintf(f11, " have, need = %12d%12d\n", nwork, nwork-nroom);
                        fprintf(f11, " not writing out plot3d files\n");
                    }
                } else {
                    FortranArray4DRef<double> q_arr(&w(lq), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qi0_arr(&w(lqi0), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qj0_arr(&w(lqj0), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qk0_arr(&w(lqk0), jdim, kdim, idim, 5);
                    FortranArray3DRef<double> x_arr(&w(lx), jdim, kdim, idim);
                    FortranArray3DRef<double> y_arr(&w(ly), jdim, kdim, idim);
                    FortranArray3DRef<double> z_arr(&w(lz), jdim, kdim, idim);
                    FortranArray4DRef<double> xw_arr(&wk(ixwk), jdw, kdw, idw, 5);
                    FortranArray3DRef<double> blank2_arr(&wk(ibwk), jdim, kdim, idim);
                    FortranArray3DRef<double> blank_arr(&w(lblk), jdim, kdim, idim);
                    FortranArray4DRef<double> xg_arr(&wk(ixgk), jdw, kdw, idw, 4);
                    FortranArray3DRef<double> vis_arr(&w(lvis), jdim, kdim, idim);
                    FortranArray4DRef<double> vi0_arr(&w(lvi0), jdim, kdim, idim, 1);
                    FortranArray4DRef<double> vj0_arr(&w(lvj0), jdim, kdim, idim, 1);
                    FortranArray4DRef<double> vk0_arr(&w(lvk0), jdim, kdim, idim, 1);
                    FortranArray3DRef<double> smin_arr(&w(lsnk0), jdim, kdim, idim);
                    plot3c_ns::plot3c(jdim, kdim, idim, i1, i2, i3, j1, j2, j3,
                                      k1, k2, k3, q_arr, qi0_arr, qj0_arr, qk0_arr,
                                      x_arr, y_arr, z_arr, xw_arr, blank2_arr,
                                      blank_arr, xg_arr, iflag, vis_arr,
                                      vi0_arr, vj0_arr, vk0_arr,
                                      iovrlp(nbl), nbl, nmap, smin_arr,
                                      ifunc, n, jdw, kdw, idw,
                                      nplots, jdimg, kdimg, idimg, nblcg,
                                      jsg, ksg, isg, jeg, keg, ieg,
                                      ninter, iindex, intmax, nsub1, maxxe,
                                      nblk, nbli, limblk, isva, nblon, mxbli,
                                      thetay, maxbl, maxgr, myid, myhost, mycomm,
                                      mblk2nd, inpl3d, nblock, nblkpt, ip3dsurf);
                }
            }
        } else if (0 == 2) {
            // never executes: cell center turbulence data
            if (mblk2nd(nbl) == myid || myid == myhost) {
                if (ivmx > 1) {
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    i2 = std::min(idim-1, i2);
                    j2 = std::min(jdim-1, j2);
                    k2 = std::min(kdim-1, k2);
                    i1 = std::min(idim-1, i1);
                    j1 = std::min(jdim-1, j1);
                    k1 = std::min(kdim-1, k1);
                    jdw = (j2-j1)/j3 + 1;
                    kdw = (k2-k1)/k3 + 1;
                    idw = (i2-i1)/i3 + 1;
                    if (ivmx == 8 || ivmx == 9 || ivmx >= 11) {
                        ixwk  = 1;
                        ibwk  = ixwk + jdw*kdw*idw*5;
                        ixgk  = ibwk + jdim*kdim*idim;
                        ibwk3 = ixgk + jdw*kdw*idw*4;
                        nroom = nwork - (ibwk3 + jdim*kdim*9);
                        if (nroom < 0) {
                            if (myid == myhost) {
                                FILE* f11 = fortran_get_unit(11);
                                fprintf(f11, " not enough memory for plot3t\n");
                                fprintf(f11, " have, need = %12d%12d\n", nwork, nwork-nroom);
                                fprintf(f11, " not writing out plot3d files\n");
                            }
                        } else {
                            FortranArray4DRef<double> q_arr(&w(lq), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> x_arr(&w(lx), jdim, kdim, idim);
                            FortranArray3DRef<double> y_arr(&w(ly), jdim, kdim, idim);
                            FortranArray3DRef<double> z_arr(&w(lz), jdim, kdim, idim);
                            FortranArray4DRef<double> xw_arr(&wk(ixwk), jdw, kdw, idw, 5);
                            FortranArray3DRef<double> blank2_arr(&wk(ibwk), jdim, kdim, idim);
                            FortranArray3DRef<double> blank_arr(&w(lblk), jdim, kdim, idim);
                            FortranArray4DRef<double> xg_arr(&wk(ixgk), jdw, kdw, idw, 4);
                            FortranArray3DRef<double> vis_arr(&w(lvis), jdim, kdim, idim);
                            FortranArray4DRef<double> sj_arr(&w(lsj), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> sk_arr(&w(lsk), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> si_arr(&w(lsi), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> smin_arr(&w(lsnk0), jdim, kdim, idim);
                            FortranArray4DRef<double> ux_arr(&w(lux), jdim, kdim, idim, 9);
                            FortranArray4DRef<double> turre_arr(&w(lxib), jdim, kdim, idim, 1);
                            FortranArray3DRef<double> vol_arr(&w(lvol), jdim, kdim, idim);
                            FortranArray4DRef<double> qj0_arr(&w(lqj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> qk0_arr(&w(lqk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> qi0_arr(&w(lqi0), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> bcj_arr(&w(lbcj), jdim, kdim, 2);
                            FortranArray3DRef<double> bck_arr(&w(lbck), jdim, kdim, 2);
                            FortranArray3DRef<double> bci_arr(&w(lbci), jdim, kdim, 2);
                            FortranArray3DRef<double> wt_arr(&wk(ibwk3), jdim, kdim, 9);
                            FortranArray3DRef<double> cmuv_arr(&w(lcmuv), jdim, kdim, idim);
                            FortranArray3DRef<double> volj0_arr(&w(lvolj0), jdim, kdim, idim);
                            FortranArray3DRef<double> volk0_arr(&w(lvolk0), jdim, kdim, idim);
                            FortranArray3DRef<double> voli0_arr(&w(lvoli0), jdim, kdim, idim);
                            int ifuncdim5 = 5;
                            plot3t_ns::plot3t(jdim, kdim, idim, i1, i2, i3,
                                              j1, j2, j3, k1, k2, k3,
                                              q_arr, x_arr, y_arr, z_arr,
                                              xw_arr, blank2_arr, blank_arr, xg_arr,
                                              vis_arr, iovrlp(nbl), nbl, nmap,
                                              sj_arr, sk_arr, si_arr, smin_arr,
                                              ux_arr, turre_arr, vol_arr,
                                              qj0_arr, qk0_arr, qi0_arr,
                                              bcj_arr, bck_arr, bci_arr,
                                              wt_arr, cmuv_arr,
                                              jdw, kdw, idw, nplots,
                                              jdimg, kdimg, idimg, nblcg,
                                              jsg, ksg, isg, ieg, jeg, keg,
                                              ninter, iindex, intmax, nsub1,
                                              maxxe, nblk, nbli, limblk, isva,
                                              nblon, mxbli, thetay, maxbl, maxgr,
                                              myid, myhost, mycomm, mblk2nd,
                                              inpl3d, nblock, nblkpt,
                                              volj0_arr, volk0_arr, voli0_arr,
                                              vormax, ivmax, jvmax, kvmax,
                                              nummem, ifuncdim5);
                        }
                    } else {


                        // need to grab ux storage from temporary work array
                        ixwk  = 1;
                        ibwk  = ixwk + jdw*kdw*idw*5;
                        ixgk  = ibwk + jdim*kdim*idim;
                        ibwk2 = ixgk + jdw*kdw*idw*4;
                        ibwk3 = ibwk2 + (jdim-1)*(kdim-1)*(idim-1)*9;
                        nroom = nwork - (ibwk3 + jdim*kdim*9);
                        if (nroom < 0) {
                            if (myid == myhost) {
                                FILE* f11 = fortran_get_unit(11);
                                fprintf(f11, " not enough memory for plot3t\n");
                                fprintf(f11, " have, need = %12d%12d\n", nwork, nwork-nroom);
                                fprintf(f11, " not writing out plot3d files\n");
                            }
                        } else {
                            FortranArray4DRef<double> q_arr(&w(lq), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> x_arr(&w(lx), jdim, kdim, idim);
                            FortranArray3DRef<double> y_arr(&w(ly), jdim, kdim, idim);
                            FortranArray3DRef<double> z_arr(&w(lz), jdim, kdim, idim);
                            FortranArray4DRef<double> xw_arr(&wk(ixwk), jdw, kdw, idw, 5);
                            FortranArray3DRef<double> blank2_arr(&wk(ibwk), jdim, kdim, idim);
                            FortranArray3DRef<double> blank_arr(&w(lblk), jdim, kdim, idim);
                            FortranArray4DRef<double> xg_arr(&wk(ixgk), jdw, kdw, idw, 4);
                            FortranArray3DRef<double> vis_arr(&w(lvis), jdim, kdim, idim);
                            FortranArray4DRef<double> sj_arr(&w(lsj), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> sk_arr(&w(lsk), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> si_arr(&w(lsi), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> smin_arr(&w(lsnk0), jdim, kdim, idim);
                            FortranArray4DRef<double> ux_arr(&wk(ibwk2), jdim-1, kdim-1, idim-1, 9);
                            FortranArray4DRef<double> turre_arr(&w(lxib), jdim, kdim, idim, 1);
                            FortranArray3DRef<double> vol_arr(&w(lvol), jdim, kdim, idim);
                            FortranArray4DRef<double> qj0_arr(&w(lqj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> qk0_arr(&w(lqk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> qi0_arr(&w(lqi0), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> bcj_arr(&w(lbcj), jdim, kdim, 2);
                            FortranArray3DRef<double> bck_arr(&w(lbck), jdim, kdim, 2);
                            FortranArray3DRef<double> bci_arr(&w(lbci), jdim, kdim, 2);
                            FortranArray3DRef<double> wt_arr(&wk(ibwk3), jdim, kdim, 9);
                            FortranArray3DRef<double> cmuv_arr(&w(lcmuv), jdim, kdim, idim);
                            FortranArray3DRef<double> volj0_arr(&w(lvolj0), jdim, kdim, idim);
                            FortranArray3DRef<double> volk0_arr(&w(lvolk0), jdim, kdim, idim);
                            FortranArray3DRef<double> voli0_arr(&w(lvoli0), jdim, kdim, idim);
                            int ifuncdim5 = 5;
                            plot3t_ns::plot3t(jdim, kdim, idim, i1, i2, i3,
                                              j1, j2, j3, k1, k2, k3,
                                              q_arr, x_arr, y_arr, z_arr,
                                              xw_arr, blank2_arr, blank_arr, xg_arr,
                                              vis_arr, iovrlp(nbl), nbl, nmap,
                                              sj_arr, sk_arr, si_arr, smin_arr,
                                              ux_arr, turre_arr, vol_arr,
                                              qj0_arr, qk0_arr, qi0_arr,
                                              bcj_arr, bck_arr, bci_arr,
                                              wt_arr, cmuv_arr,
                                              jdw, kdw, idw, nplots,
                                              jdimg, kdimg, idimg, nblcg,
                                              jsg, ksg, isg, ieg, jeg, keg,
                                              ninter, iindex, intmax, nsub1,
                                              maxxe, nblk, nbli, limblk, isva,
                                              nblon, mxbli, thetay, maxbl, maxgr,
                                              myid, myhost, mycomm, mblk2nd,
                                              inpl3d, nblock, nblkpt,
                                              volj0_arr, volk0_arr, voli0_arr,
                                              vormax, ivmax, jvmax, kvmax,
                                              nummem, ifuncdim5);
                        }
                    } // end if ivmx==8||9||>=11 else
                } // end if ivmx > 1
            } // end if mblk2nd
        } // end if/else if/else if

        if (iblnkfr == 0) {
            // reset blank values at fringe points to 0.
            if (mblk2nd(nbl) == myid && iovrlp(nbl) != 0) {
                FortranArray3DRef<double> blank_arr(&w(lblk), jdim, kdim, idim);
                double blnkval0 = 0.0;
                blnkfr_ns::blnkfr(nbl, iibg, kkbg, jjbg, ibpntsg, lbg, iitot,
                                   blank_arr, jdim, kdim, idim, maxbl, blnkval0);
            }
        }

    } // end do 70

    if (myid == myhost) {
        if (lhdr > 0) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, " \n");
        }
    }

    // for stationary grid cases, set icallcrs flag to prevent output
    // of grid to plot3d file on subsequent calls
    if (iunst == 0) {
        icallcrs = 1;
    } else {
        icallcrs = 0;
    }

    // 231 continue (label is a no-op here)



    // if (0.gt.0) then  -- never executes: print solution data section
    if (0 > 0) {

        if (myid == myhost) {
            if (lhdr > 0 && np3d <= 0) {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, " \n");
            }
        }

        if (myid == myhost) {
            FILE* f17 = fortran_get_unit(17);
            // write(17,111)(real(title(i)),i=1,20)
            // 111 format(/,20a4)
            fprintf(f17, "\n");
            for (int i = 1; i <= 20; i++) {
                char buf[5];
                std::memcpy(buf, &title[i-1], 4);
                buf[4] = '\0';
                fprintf(f17, "%4.4s", buf);
            }
            fprintf(f17, "\n");
            // write(17,21)
            // 21 format(6x,4hMach,5x,5halpha,6x,4hbeta,6x,4hReUe,3x,7hTinf,dR,6x,4htime)
            fprintf(f17, "      Mach     alpha      beta      ReUe   Tinf,dR      time\n");
            float alphaw = (float)radtodeg * alpha;
            float betaw  = (float)radtodeg * beta;
            // write(17,20)real(xmach),real(alphaw),real(betaw),real(reue),real(tinf),real(time)
            // 20 format(3f10.5,e10.3,2f10.5)
            fprintf(f17, "%10.5f%10.5f%10.5f%10.3e%10.5f%10.5f\n",
                    (float)xmach, (float)alphaw, (float)betaw,
                    (float)reue, (float)tinf, (float)time);
        }

        // do 80 n=1,nprint
        for (n = 1; n <= nprint; n++) {

            nbl = inpr(n,1);

            if (nbl > nblock) {
                if (myid == myhost) {
                    FILE* f11 = fortran_get_unit(11);
                    // 77 format(6h Block,i3,36h does not exist.  No output printed.)
                    fprintf(f11, " Block%3d does not exist.  No output printed.\n", nbl);
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
                        FILE* f17 = fortran_get_unit(17);
                        // 348 format(//1x,5hBLOCK,i4,2x,5h(GRID,i4,1h),5x,15hIDIM,JDIM,KDIM=,3I5)
                        fprintf(f17, "\n\n BLOCK%4d  (GRID%4d)     IDIM,JDIM,KDIM=%5d%5d%5d\n",
                                nbl, igridg(nbl), idim, jdim, kdim);
                        // 349 format(1x,32hNOTE: endpts may not be reliable)
                        fprintf(f17, " NOTE: endpts may not be reliable\n");
                    }

                    jdw = (j2-j1)/j3 + 1;
                    kdw = (k2-k1)/k3 + 1;
                    idw = (i2-i1)/i3 + 1;

                    ixwk = 1;
                    nset = 8;
                    ibwk = ixwk + jdw*kdw*idw*nset;
                    ixgk = ibwk + jdim*kdim*idim*2;
                    ixvk = ixgk + jdw*kdw*idw*4;
                    nroom = nwork - (ixvk + jdw*kdw*idw*5 - 1);
                    if (nroom < 0) {
                        if (myid == myhost) {
                            FILE* f11 = fortran_get_unit(11);
                            fprintf(f11, " not enough memory for plot3d\n");
                            fprintf(f11, " have, need = %12d%12d\n", nwork, nwork-nroom);
                            fprintf(f11, " not writing out prout file\n");
                        }
                    } else {
                        FortranArray4DRef<double> q_arr(&w(lq), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qi0_arr(&w(lqi0), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qj0_arr(&w(lqj0), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qk0_arr(&w(lqk0), jdim, kdim, idim, 5);
                        FortranArray3DRef<double> x_arr(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> y_arr(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> z_arr(&w(lz), jdim, kdim, idim);
                        FortranArray4DRef<double> xw_arr(&wk(ixwk), jdw, kdw, idw, nset);
                        FortranArray4DRef<double> blank2_arr(&wk(ibwk), jdim, kdim, idim, 2);
                        FortranArray3DRef<double> blank_arr(&w(lblk), jdim, kdim, idim);
                        FortranArray4DRef<double> xg_arr(&wk(ixgk), jdw, kdw, idw, 4);
                        FortranArray3DRef<double> vis_arr(&w(lvis), jdim, kdim, idim);
                        FortranArray3DRef<double> bcj_arr(&w(lbcj), jdim, kdim, 2);
                        FortranArray3DRef<double> bck_arr(&w(lbck), jdim, kdim, 2);
                        FortranArray3DRef<double> bci_arr(&w(lbci), jdim, kdim, 2);
                        FortranArray4DRef<double> vj0_arr(&w(lvj0), jdim, kdim, idim, 1);
                        FortranArray4DRef<double> vk0_arr(&w(lvk0), jdim, kdim, idim, 1);
                        FortranArray4DRef<double> vi0_arr(&w(lvi0), jdim, kdim, idim, 1);
                        FortranArray4DRef<double> xv_arr(&wk(ixvk), jdw, kdw, idw, 5);
                        FortranArray4DRef<double> sj_arr(&w(lsj), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> sk_arr(&w(lsk), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> si_arr(&w(lsi), jdim, kdim, idim, 5);
                        FortranArray3DRef<double> vol_arr(&w(lvol), jdim, kdim, idim);
                        plot3d_coarse_ns::plot3d_coarse(jdim, kdim, idim, i1, i2, i3,
                                                        j1, j2, j3, k1, k2, k3,
                                                        q_arr, qi0_arr, qj0_arr, qk0_arr,
                                                        x_arr, y_arr, z_arr,
                                                        xw_arr, blank2_arr, blank_arr,
                                                        xg_arr, iflag, vis_arr,
                                                        iovrlp(nbl), nbl, nmap,
                                                        bcj_arr, bck_arr, bci_arr,
                                                        vj0_arr, vk0_arr, vi0_arr,
                                                        ifunc, n, jdw, kdw, idw,
                                                        nplots, jdimg, kdimg, idimg,
                                                        nblcg, jsg, ksg, isg,
                                                        jeg, keg, ieg, ninter, iindex,
                                                        intmax, nsub1, maxxe,
                                                        nblk, nbli, limblk, isva, nblon,
                                                        mxbli, thetay, maxbl, maxgr,
                                                        myid, myhost, mycomm, mblk2nd,
                                                        inpl3d, nblock, nblkpt,
                                                        xv_arr, sj_arr, sk_arr, si_arr,
                                                        vol_arr, nset);
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
                    jdw = (j2-j1)/j3 + 1;
                    kdw = (k2-k1)/k3 + 1;
                    idw = (i2-i1)/i3 + 1;

                    ixwk = 1;
                    ibwk = ixwk + jdw*kdw*idw*5;
                    ixgk = ibwk + jdim*kdim*idim;
                    nroom = nwork - (ixgk + jdw*kdw*idw*4 - 1);
                    if (nroom < 0) {
                        if (myid == myhost) {
                            FILE* f11 = fortran_get_unit(11);
                            fprintf(f11, " not enough memory for plot3c\n");
                            fprintf(f11, " have, need = %12d%12d\n", nwork, nwork-nroom);
                            fprintf(f11, " not writing out plot3d files\n");
                        }
                    } else {
                        FortranArray4DRef<double> q_arr(&w(lq), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qi0_arr(&w(lqi0), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qj0_arr(&w(lqj0), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qk0_arr(&w(lqk0), jdim, kdim, idim, 5);
                        FortranArray3DRef<double> x_arr(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> y_arr(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> z_arr(&w(lz), jdim, kdim, idim);
                        FortranArray4DRef<double> xw_arr(&wk(ixwk), jdw, kdw, idw, 5);
                        FortranArray3DRef<double> blank2_arr(&wk(ibwk), jdim, kdim, idim);
                        FortranArray3DRef<double> blank_arr(&w(lblk), jdim, kdim, idim);
                        FortranArray4DRef<double> xg_arr(&wk(ixgk), jdw, kdw, idw, 4);
                        FortranArray3DRef<double> vis_arr(&w(lvis), jdim, kdim, idim);
                        FortranArray4DRef<double> vi0_arr(&w(lvi0), jdim, kdim, idim, 1);
                        FortranArray4DRef<double> vj0_arr(&w(lvj0), jdim, kdim, idim, 1);
                        FortranArray4DRef<double> vk0_arr(&w(lvk0), jdim, kdim, idim, 1);
                        FortranArray3DRef<double> smin_arr(&w(lsnk0), jdim, kdim, idim);
                        plot3c_ns::plot3c(jdim, kdim, idim, i1, i2, i3, j1, j2, j3,
                                          k1, k2, k3, q_arr, qi0_arr, qj0_arr, qk0_arr,
                                          x_arr, y_arr, z_arr, xw_arr, blank2_arr,
                                          blank_arr, xg_arr, iflag, vis_arr,
                                          vi0_arr, vj0_arr, vk0_arr,
                                          iovrlp(nbl), nbl, nmap, smin_arr,
                                          ifunc, n, jdw, kdw, idw,
                                          nplots, jdimg, kdimg, idimg, nblcg,
                                          jsg, ksg, isg, jeg, keg, ieg,
                                          ninter, iindex, intmax, nsub1, maxxe,
                                          nblk, nbli, limblk, isva, nblon, mxbli,
                                          thetay, maxbl, maxgr, myid, myhost, mycomm,
                                          mblk2nd, inpl3d, nblock, nblkpt, iprnsurf);
                    }
                }
            } // end if inpr(n,2)==0

        } // end do 80

        if (myid == myhost) {
            if (lhdr > 0) {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, " \n");
            }
        }

    } // end if (0.gt.0)

    return;
}

} // namespace qout_coarse_ns
