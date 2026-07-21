// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include <cstdio>
#include "mgblk.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include "lead.h"
#include "termn8.h"
#include "updateg.h"
#include "resetg.h"
#include "updatedg.h"
#include "dynptch.h"
#include "pre_bc.h"
#include "qout.h"
#include "qoutavg.h"
#include "qout_coarse.h"
#include "qout_2d.h"
#include "bc.h"
#include "bc_info.h"
#include "bc_xmera.h"
#include "bc_blkint.h"
#include "bc_embed.h"
#include "bc_patch.h"
#include "bc_period.h"
#include "bcchk.h"
#include "chkdef.h"
#include "writ_buf.h"
#include "fa.h"
#include "resid.h"
#include "resadd.h"
#include "tau.h"
#include "resp.h"
#include "rsmooth.h"
#include "update.h"
#include "intrbc.h"
#include "collq.h"
#include "collqc0.h"
#include "coll2q.h"
#include "addx.h"
#include "add2x.h"
#include "newalpha.h"
#include "rotateq.h"
#include "setdqc0.h"
#include "setqc0.h"
#include "wrest.h"
#include "wrestg.h"
#include "ae_pred.h"
#include "ae_corr.h"
#include "rb_pred.h"
#include "rb_corr.h"
#include "mgbl.h"
#include <cstring>
#include <cmath>
#include <algorithm>

namespace mgblk_ns {

void mgblk(int& iseq, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
    FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> timesave,
    FortranArray1DRef<double> wk, int& nwork, FortranArray1DRef<int> iw,
    FortranArray1DRef<int> ireq_qb, FortranArray1DRef<int> iwk, int& iwork,
    int& maxbl, int& maxgr, int& maxseg, int& iitot, int& intmax, int& nsub1,
    int& maxxe, int& ncycmax, FortranArray1DRef<int> iovrlp,
    FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg,
    FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg,
    FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg,
    FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg,
    FortranArray1DRef<double> dxintg, FortranArray1DRef<double> dyintg,
    FortranArray1DRef<double> dzintg, FortranArray1DRef<int> iiig,
    FortranArray1DRef<int> jjig, FortranArray1DRef<int> kkig,
    FortranArray3DRef<double> qb, FortranArray3DRef<int> lwdat,
    FortranArray2DRef<int> nblk, int& nbli, FortranArray3DRef<int> limblk,
    FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, int& nblock,
    FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg,
    FortranArray2DRef<int> iviscg, FortranArray1DRef<int> idimg,
    FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg,
    FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg,
    FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg,
    FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg,
    FortranArray1DRef<int> nblcg, FortranArray2DRef<int> mit,
    FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj,
    FortranArray4DRef<double> bcvalk, FortranArray1DRef<int> nbci0,
    FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0,
    FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0,
    FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo,
    FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
    FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej,
    FortranArray3DRef<int> bcfilek, FortranArray1DRef<double> utrans,
    FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans,
    FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay,
    FortranArray1DRef<double> omegaz, FortranArray1DRef<double> xorig,
    FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig,
    FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx,
    FortranArray1DRef<double> dzmx, FortranArray1DRef<double> dthxmx,
    FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx,
    FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay,
    FortranArray1DRef<double> thetaz, FortranArray1DRef<double> rfreqt,
    FortranArray1DRef<double> rfreqr, FortranArray1DRef<double> xorig0,
    FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0,
    FortranArray1DRef<double> time2, FortranArray1DRef<double> thetaxl,
    FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl,
    FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
    FortranArray1DRef<int> idefrm, FortranArray1DRef<double> rms,
    FortranArray1DRef<double> clw, FortranArray1DRef<double> cdw,
    FortranArray1DRef<double> cdpw, FortranArray1DRef<double> cdvw,
    FortranArray1DRef<double> cxw, FortranArray1DRef<double> cyw,
    FortranArray1DRef<double> czw, FortranArray1DRef<double> cmxw,
    FortranArray1DRef<double> cmyw, FortranArray1DRef<double> cmzw,
    int& n_clcd, FortranArray3DRef<double> clcd, int& nblocks_clcd,
    FortranArray2DRef<int> blocks_clcd, FortranArray1DRef<double> chdw,
    FortranArray1DRef<double> swetw, FortranArray1DRef<double> fmdotw,
    FortranArray1DRef<double> cfttotw, FortranArray1DRef<double> cftmomw,
    FortranArray1DRef<double> cftpw, FortranArray1DRef<double> cftvw,
    FortranArray2DRef<double> rmstr, FortranArray2DRef<int> nneg,
    int& ntr, int& ihstry, FortranArray1DRef<int> iadvance,
    FortranArray1DRef<int> iforce, int& lfgm,
    FortranArray1DRef<double> resmx, FortranArray1DRef<int> imx,
    FortranArray1DRef<int> jmx, FortranArray1DRef<int> kmx,
    FortranArray1DRef<double> vormax, FortranArray1DRef<int> ivmax,
    FortranArray1DRef<int> jvmax, FortranArray1DRef<int> kvmax,
    FortranArray1DRef<double> sx, FortranArray1DRef<double> sy,
    FortranArray1DRef<double> sz, FortranArray1DRef<double> stot,
    FortranArray1DRef<double> pav, FortranArray1DRef<double> ptav,
    FortranArray1DRef<double> tav, FortranArray1DRef<double> ttav,
    FortranArray1DRef<double> xmav, FortranArray1DRef<double> fmdot,
    FortranArray1DRef<double> cfxp, FortranArray1DRef<double> cfyp,
    FortranArray1DRef<double> cfzp, FortranArray1DRef<double> cfdp,
    FortranArray1DRef<double> cflp, FortranArray1DRef<double> cftp,
    FortranArray1DRef<double> cfxv, FortranArray1DRef<double> cfyv,
    FortranArray1DRef<double> cfzv, FortranArray1DRef<double> cfdv,
    FortranArray1DRef<double> cflv, FortranArray1DRef<double> cftv,
    FortranArray1DRef<double> cfxmom, FortranArray1DRef<double> cfymom,
    FortranArray1DRef<double> cfzmom, FortranArray1DRef<double> cfdmom,
    FortranArray1DRef<double> cflmom, FortranArray1DRef<double> cftmom,
    FortranArray1DRef<double> cfxtot, FortranArray1DRef<double> cfytot,
    FortranArray1DRef<double> cfztot, FortranArray1DRef<double> cfdtot,
    FortranArray1DRef<double> cfltot, FortranArray1DRef<double> cfttot,
    FortranArray2DRef<int> icsinfo, int& ncs,
    FortranArray2DRef<double> windex, int& ninter,
    FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt,
    FortranArray2DRef<double> windx, int& nintr,
    FortranArray2DRef<int> iindx, FortranArray1DRef<int> llimit,
    FortranArray1DRef<int> iitmax, FortranArray1DRef<int> mmcxie,
    FortranArray1DRef<int> mmceta, FortranArray1DRef<int> ncheck,
    FortranArray1DRef<int> iifit, FortranArray1DRef<int> mblkpt,
    FortranArray1DRef<int> iic0, FortranArray1DRef<int> iiorph,
    FortranArray1DRef<int> iitoss, FortranArray1DRef<int> ifiner,
    int& msub1, FortranArray2DRef<double> dthetxx,
    FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz,
    FortranArray1DRef<double> swett, FortranArray1DRef<double> clt,
    FortranArray1DRef<double> cdt, FortranArray1DRef<double> cxt,
    FortranArray1DRef<double> cyt, FortranArray1DRef<double> czt,
    FortranArray1DRef<double> cmxt, FortranArray1DRef<double> cmyt,
    FortranArray1DRef<double> cmzt, FortranArray1DRef<double> cdpt,
    FortranArray1DRef<double> cdvt, FortranArray1DRef<int> mblk2nd,
    FortranArray1DRef<double> geom_miss, double& epsc0,
    FortranArray1DRef<double> period_miss, double& epsrot,
    FortranArray2DRef<int> isav_blk, FortranArray2DRef<int> isav_prd,
    int& lbcprd, FortranArray2DRef<int> isav_pat,
    FortranArray3DRef<int> isav_pat_b, FortranArray2DRef<int> isav_emb,
    int& lbcemb, int& mxbli, int& maxcs, int& intmx, int& mxxe, int& mptch,
    FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg,
    FortranArray1DRef<int> iemg, int& ngrid,
    FortranArray2DRef<double> dx, FortranArray2DRef<double> dy,
    FortranArray2DRef<double> dz, FortranArray2DRef<double> dthetx,
    FortranArray2DRef<double> dthety, FortranArray2DRef<double> dthetz,
    FortranArray2DRef<int> isav_dpat, FortranArray3DRef<int> isav_dpat_b,
    FortranArray1DRef<int> lout, FortranArray1DRef<int> ifrom,
    FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2,
    FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2,
    FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1,
    FortranArray1DRef<int> iiint1, FortranArray1DRef<int> iiint2,
    FortranArray1DRef<int> nblk1, FortranArray1DRef<int> nblk2,
    FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage,
    FortranArray1DRef<int> jte, FortranArray1DRef<int> kte,
    FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm,
    FortranArray3DRef<double> xte, FortranArray3DRef<double> yte,
    FortranArray3DRef<double> zte, FortranArray3DRef<double> xmi,
    FortranArray3DRef<double> ymi, FortranArray3DRef<double> zmi,
    FortranArray3DRef<double> xmie, FortranArray3DRef<double> ymie,
    FortranArray3DRef<double> zmie, FortranArray3DRef<double> sxie,
    FortranArray3DRef<double> seta, FortranArray2DRef<double> sxie2,
    FortranArray2DRef<double> seta2, FortranArray2DRef<double> xie2s,
    FortranArray2DRef<double> eta2s, FortranArray1DRef<double> temp,
    FortranArray2DRef<double> x2, FortranArray2DRef<double> y2,
    FortranArray2DRef<double> z2, FortranArray2DRef<double> x1,
    FortranArray2DRef<double> y1, FortranArray2DRef<double> z1,
    int& ip3dsurf, FortranArray2DRef<double> factjlo,
    FortranArray2DRef<double> factjhi, FortranArray2DRef<double> factklo,
    FortranArray2DRef<double> factkhi, int& nplots,
    FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
    int& nbuf, int& ibufdim, FortranArray2DRef<int> istat2_bl,
    FortranArray2DRef<int> istat2_pa, FortranArray2DRef<int> istat2_pe,
    FortranArray2DRef<int> istat2_em, int& istat_size,
    int& nplot3d, int& nprint, FortranArray2DRef<int> inpl3d,
    FortranArray2DRef<int> inpr, FortranArray1DRef<char[80]> bcfiles,
    int& mxbcfil, FortranArray2DRef<double> utrnsae,
    FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae,
    FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae,
    FortranArray2DRef<double> omgzae, FortranArray2DRef<double> xorgae,
    FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae,
    FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae,
    FortranArray2DRef<double> thtzae, FortranArray2DRef<double> rfrqtae,
    FortranArray2DRef<double> rfrqrae, FortranArray2DRef<int> icsi,
    FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi,
    FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi,
    FortranArray2DRef<int> kcsf, FortranArray2DRef<double> freq,
    FortranArray2DRef<double> gmass, FortranArray2DRef<double> damp,
    FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0,
    int& nmds, int& maxaes, FortranArray2DRef<double> aesrfdat,
    FortranArray3DRef<double> perturb, FortranArray3DRef<int> islavept,
    int& nslave, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip,
    FortranArray2DRef<int> kskip, FortranArray3DRef<double> bmat,
    FortranArray3DRef<double> stm, FortranArray3DRef<double> stmi,
    FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn,
    FortranArray2DRef<double> gforcn, FortranArray2DRef<double> gforcnm,
    FortranArray2DRef<double> gforcs, FortranArray1DRef<int> nsegdfrm,
    FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf,
    int& maxsegdg, int& nmaster, FortranArray4DRef<double> aehist,
    FortranArray1DRef<double> timekeep, FortranArray2DRef<double> xorgae0,
    FortranArray2DRef<double> yorgae0, FortranArray2DRef<double> zorgae0,
    FortranArray2DRef<int> icouple, int& iprnsurf,
    FortranArray2DRef<int> nblelst, FortranArray1DRef<int> iskmax,
    FortranArray1DRef<int> jskmax, FortranArray1DRef<int> kskmax,
    FortranArray1DRef<double> ue, int& nummem) {
    // COMMON block aliases
    int32_t& levt_c    = cmn_mgrd.levt;
    int32_t& kode      = cmn_mgrd.kode;
    int32_t& mode      = cmn_mgrd.mode;
    int32_t& ncyc      = cmn_mgrd.ncyc;
    int32_t& mtt       = cmn_mgrd.mtt;
    int32_t& icyc      = cmn_mgrd.icyc;
    int32_t& level     = cmn_mgrd.level;
    int32_t& lglobal   = cmn_mgrd.lglobal;

    int32_t& jdim      = cmn_ginfo.jdim;
    int32_t& kdim      = cmn_ginfo.kdim;
    int32_t& idim      = cmn_ginfo.idim;
    int32_t& jj2       = cmn_ginfo.jj2;
    int32_t& kk2       = cmn_ginfo.kk2;
    int32_t& ii2       = cmn_ginfo.ii2;
    int32_t& nblc      = cmn_ginfo.nblc;
    int32_t& js        = cmn_ginfo.js;
    int32_t& ks        = cmn_ginfo.ks;
    int32_t& is        = cmn_ginfo.is;
    int32_t& je        = cmn_ginfo.je;
    int32_t& ke        = cmn_ginfo.ke;
    int32_t& ie        = cmn_ginfo.ie;
    int32_t& lq        = cmn_ginfo.lq;
    int32_t& lqj0      = cmn_ginfo.lqj0;
    int32_t& lqk0      = cmn_ginfo.lqk0;
    int32_t& lqi0      = cmn_ginfo.lqi0;
    int32_t& lsj       = cmn_ginfo.lsj;
    int32_t& lsk       = cmn_ginfo.lsk;
    int32_t& lsi       = cmn_ginfo.lsi;
    int32_t& lvol      = cmn_ginfo.lvol;
    int32_t& ldtj      = cmn_ginfo.ldtj;
    int32_t& lx        = cmn_ginfo.lx;
    int32_t& ly        = cmn_ginfo.ly;
    int32_t& lz        = cmn_ginfo.lz;
    int32_t& lvis      = cmn_ginfo.lvis;
    int32_t& lsnk0     = cmn_ginfo.lsnk0;
    int32_t& lsni0     = cmn_ginfo.lsni0;
    int32_t& lq1       = cmn_ginfo.lq1;
    int32_t& lqr       = cmn_ginfo.lqr;
    int32_t& lblk      = cmn_ginfo.lblk;
    int32_t& lxib      = cmn_ginfo.lxib;
    int32_t& lsig      = cmn_ginfo.lsig;
    int32_t& lsqtq     = cmn_ginfo.lsqtq;
    int32_t& lg        = cmn_ginfo.lg;
    int32_t& ltj0      = cmn_ginfo.ltj0;
    int32_t& ltk0      = cmn_ginfo.ltk0;
    int32_t& lti0      = cmn_ginfo.lti0;
    int32_t& lxkb      = cmn_ginfo.lxkb;
    int32_t& lnbl      = cmn_ginfo.lnbl;
    int32_t& lvj0      = cmn_ginfo.lvj0;
    int32_t& lvk0      = cmn_ginfo.lvk0;
    int32_t& lvi0      = cmn_ginfo.lvi0;
    int32_t& lbcj      = cmn_ginfo.lbcj;
    int32_t& lbck      = cmn_ginfo.lbck;
    int32_t& lbci      = cmn_ginfo.lbci;
    int32_t& lqc0      = cmn_ginfo.lqc0;
    int32_t& ldqc0     = cmn_ginfo.ldqc0;
    int32_t& lxtbi     = cmn_ginfo.lxtbi;
    int32_t& lxtbj     = cmn_ginfo.lxtbj;
    int32_t& lxtbk     = cmn_ginfo.lxtbk;
    int32_t& latbi     = cmn_ginfo.latbi;
    int32_t& latbj     = cmn_ginfo.latbj;
    int32_t& latbk     = cmn_ginfo.latbk;
    int32_t& lbcdj     = cmn_ginfo.lbcdj;
    int32_t& lbcdk     = cmn_ginfo.lbcdk;
    int32_t& lbcdi     = cmn_ginfo.lbcdi;
    int32_t& lxib2     = cmn_ginfo.lxib2;
    int32_t& lux       = cmn_ginfo.lux;
    int32_t& lcmuv     = cmn_ginfo.lcmuv;
    int32_t& lvolj0    = cmn_ginfo.lvolj0;
    int32_t& lvolk0    = cmn_ginfo.lvolk0;
    int32_t& lvoli0    = cmn_ginfo.lvoli0;
    int32_t& lxmdj     = cmn_ginfo.lxmdj;
    int32_t& lxmdk     = cmn_ginfo.lxmdk;
    int32_t& lxmdi     = cmn_ginfo.lxmdi;
    int32_t& lvelg     = cmn_ginfo.lvelg;
    int32_t& ldeltj    = cmn_ginfo.ldeltj;
    int32_t& ldeltk    = cmn_ginfo.ldeltk;
    int32_t& ldelti    = cmn_ginfo.ldelti;
    int32_t& lxnm2     = cmn_ginfo.lxnm2;
    int32_t& lynm2     = cmn_ginfo.lynm2;
    int32_t& lznm2     = cmn_ginfo.lznm2;
    int32_t& lxnm1     = cmn_ginfo.lxnm1;
    int32_t& lynm1     = cmn_ginfo.lynm1;
    int32_t& lznm1     = cmn_ginfo.lznm1;
    int32_t& lqavg     = cmn_ginfo.lqavg;

    int32_t& lq2avg    = cmn_ginfo2.lq2avg;
    int32_t& iskip_blocks = cmn_ginfo2.iskip_blocks;
    // inc_2d(1..3) and inc_coarse(1..3) are 1-based in Fortran
    // cmn_ginfo2.inc_2d[0] = inc_2d(1), etc.

    float& xnumavg     = cmn_avgdata.xnumavg;
    int32_t& iteravg   = cmn_avgdata.iteravg;
    float& xnumavg2    = cmn_avgdata.xnumavg2;
    int32_t& ipertavg  = cmn_avgdata.ipertavg;
    int32_t& iclcd     = cmn_avgdata.iclcd;
    int32_t& isubit_r  = cmn_avgdata.isubit_r;

    float& xmach       = cmn_info.xmach;
    float& alpha       = cmn_info.alpha;
    float& beta_c      = cmn_info.beta;
    float& dt          = cmn_info.dt;
    float& fmax        = cmn_info.fmax;
    int32_t& nit       = cmn_info.nit;
    int32_t& ntt       = cmn_info.ntt;
    int32_t& nitfo     = cmn_info.nitfo;
    int32_t& iflagts   = cmn_info.iflagts;
    int32_t& nres      = cmn_info.nres;
    int32_t& mgflag    = cmn_info.mgflag;
    int32_t& iconsf    = cmn_info.iconsf;
    int32_t& mseq      = cmn_info.mseq;
    int32_t& ngam      = cmn_info.ngam;
    int32_t& iipv      = cmn_info.iipv;

    int32_t& issc      = cmn_mgv.issc;
    int32_t& issr      = cmn_mgv.issr;

    int32_t& isklton   = cmn_sklton.isklton;

    float& time        = cmn_unst.time;
    float& cfltau      = cmn_unst.cfltau;
    int32_t& ntstep    = cmn_unst.ntstep;
    int32_t& ita       = cmn_unst.ita;
    int32_t& iunst     = cmn_unst.iunst;
    float& cfltau0     = cmn_unst.cfltau0;
    float& cfltaumax   = cmn_unst.cfltaumax;

    int32_t& nwrest    = cmn_wrbl.nwrest;

    int32_t& movie     = cmn_moov.movie;
    int32_t& nframes   = cmn_moov.nframes;
    int32_t& icall1    = cmn_moov.icall1;
    int32_t& lhdr      = cmn_moov.lhdr;
    int32_t& icoarsemovie = cmn_moov.icoarsemovie;
    int32_t& i2dmovie  = cmn_moov.i2dmovie;

    char* restrt       = cmn_filenam.restrt;
    char* avgg         = cmn_filenam.avgg;
    char* avgq         = cmn_filenam.avgq;
    char* avgq2        = cmn_filenam2.avgq2;
    char* avgq2pert    = cmn_filenam2.avgq2pert;
    char* clcds        = cmn_filenam2.clcds;

    int32_t& ialphit   = cmn_alphait.ialphit;
    float& cltarg      = cmn_alphait.cltarg;
    float& rlxalph     = cmn_alphait.rlxalph;
    float& dalim       = cmn_alphait.dalim;
    float& dalpha      = cmn_alphait.dalpha;
    int32_t& icycupdt  = cmn_alphait.icycupdt;

    float& cprec       = cmn_precond.cprec;
    float& uref        = cmn_precond.uref;
    float& avn         = cmn_precond.avn;

    int32_t& nnodes    = cmn_mydist2.nnodes;
    int32_t& myhost    = cmn_mydist2.myhost;
    int32_t& myid      = cmn_mydist2.myid;
    int32_t& mycomm    = cmn_mydist2.mycomm;

    int32_t& ivmx      = cmn_maxiv.ivmx;

    float& time2mc     = cmn_motionmc.time2mc;
    int32_t& itransmc  = cmn_motionmc.itransmc;
    int32_t& irotatmc  = cmn_motionmc.irotatmc;

    int32_t& irigb     = cmn_rigidbody.irigb;
    int32_t& irbtrim   = cmn_rigidbody.irbtrim;

    int32_t& iwarneddy = cmn_turbconv.iwarneddy;
    int32_t& itime2read= cmn_turbconv.itime2read;
    int32_t& itaturb   = cmn_turbconv.itaturb;
    int32_t& nsubturb  = cmn_turbconv.nsubturb;
    int32_t& nfreeze   = cmn_turbconv.nfreeze;

    int32_t& meshdef   = cmn_deformz.meshdef;
    int32_t& ndwrt     = cmn_deformz.ndwrt;
    int32_t& negvol    = cmn_deformz.negvol;

    int32_t& iaccnt_c  = cmn_account.iaccnt;
    int32_t& ioutsub   = cmn_account.ioutsub;

    int32_t& ichk      = cmn_chk.ichk;

    int32_t& icgns     = cmn_cgns.icgns;

    int32_t& iunit5    = cmn_unit5.iunit5;

    int32_t& ndefrm    = cmn_elastic.ndefrm;
    int32_t& naesrf    = cmn_elastic.naesrf;

    float& dt0         = cmn_cfl.dt0;
    float& dtold       = cmn_cfl.dtold;

    // Local variables
    int termn8_err = -1;  // for passing -1 to termn8 by reference
    int writ_buf_11 = 11, writ_buf_14 = 14;  // for passing literals by reference
    bool stop_me = false;
    int levb = 0, nembed = 0, nblstag = 0, nblendg = 0, nblstat = 0, nblendt = 0, nblstat_h = 0;
    int ncall = 0, iupdat = 0, irigb0 = 0, irbtrim0 = 0;
    int iwk1 = 0, iwk2 = 0, iwk3 = 0, iwk4 = 0, iwk5 = 0, iwk6 = 0, iwork1 = 0, mneed = 0;
    int ntt0 = 0, ntime = 0, nit1 = 0, ilc = 0, iaccnt = 0;
    int lres = 0, lwj0 = 0, lwk0 = 0, lwi0 = 0, lvmuk = 0, lvmuj = 0, lvmui = 0, ltot = 0, lsafe = 0;
    int lqc = 0, lvolc = 0, lvtc = 0, lqrc = 0, lxibc = 0, lqc0c = 0, ldqc0c = 0, nroom = 0, mdim = 0;
    int nblend = 0, nblf = 0, nblz = 0, nblcc = 0, igridc = 0, ipass = 0;
    int jc = 0, kc = 0, ic = 0, lq1c = 0;
    int lxtc = 0;
    float deltat[5] = {0.0f};
    int nptsr[20] = {0};
    int nt = 0;
    icyc = 0;
    int kxpand = 0;
    double clwuse = 0.0;
    double dthtx = 0.0, dthty = 0.0, dthtz = 0.0;
    int lqwk = 0;
    int ifluxa = 0, ifamax = 0, nwfa = 0, nifa = 0;
    int ibcflg = 0, istop = 0;
    double rmsb = 0.0, rmst = 0.0;
    int nptsrb = 0;
    FortranArray1D<double> nnegb_arr(nummem), rmstb_arr(nummem), sumn_arr(nummem), negn_arr(nummem);
    FortranArray1D<int> idimt_arr(maxbl), jdimt_arr(maxbl), kdimt_arr(maxbl);
    int igrid = 0, iem = 0, iskipz = 0, idima = 0, jdima = 0, kdima = 0;
    int iflagg = 0, iaes = 0, nmodes = 0;
    int movabs = 0, iflag_coarse = 0, mov_coarse = 0, iinc_coarse = 0, jinc_coarse = 0, kinc_coarse = 0;
    int iflag_2d = 0, mov_2d = 0, iinc_2d = 0, jinc_2d = 0, kinc_2d = 0;
    int icallcrs = 0, icall2d = 0;
    double xorgrb = 0.0, yorgrb = 0.0, zorgrb = 0.0;
    double angx = 0.0, angy = 0.0, angz = 0.0;
    int itatemp = 0;
    int nblocks = nblock;
    int iprerot = 0;
    int nttuse = 0;
    int need = 0;
    int ldim = 0;
    int isf = 0;
    int icall = 0;
    int iuns = 0;
    int nn = 0, iii = 0, lll = 0, l = 0, nbl = 0, n = 0;
    int nblstag_init = 0, nblendg_init = 0, nblstat_init = 0, nblendt_init = 0;
    int mytag_qb[13] = {0};
    int istop_data[10] = {0};
    char basedesired[33] = {0};
    int isubit = 0;
    int nblend_save = 0;
    int nblf_save = 0;
    int nblz_save = 0;
    int nblcc_save = 0;
    int igridc_save = 0;
    int ipass_save = 0;
    int jc_save = 0, kc_save = 0, ic_save = 0;
    int lq1c_save = 0, lqc_save = 0, lvolc_save = 0, lqrc_save = 0, lvtc_save = 0, lxtc_save = 0;
    int lxibc_save = 0, lqc0c_save = 0, ldqc0c_save = 0;
    int nblcc_c = 0;
    int nblf_c = 0;
    int nblz_c = 0;
    int igridc_c = 0;
    int ipass_c = 0;
    int jc_c = 0, kc_c = 0, ic_c = 0;
    int lq1c_c = 0, lqc_c = 0, lvolc_c = 0, lqrc_c = 0, lvtc_c = 0, lxtc_c = 0;
    int lxibc_c = 0, lqc0c_c = 0, ldqc0c_c = 0;
    int nblend_c = 0;
    int nblf_cc = 0;
    int nblz_cc = 0;
    int igridc_cc = 0;
    int ipass_cc = 0;
    int jc_cc = 0, kc_cc = 0, ic_cc = 0;
    int lq1c_cc = 0, lqc_cc = 0, lvolc_cc = 0, lqrc_cc = 0, lvtc_cc = 0, lxtc_cc = 0;
    int lxibc_cc = 0, lqc0c_cc = 0, ldqc0c_cc = 0;
    int nblend_cc = 0;
    // Simplify: use single set of coarse-grid variables (they are reused in loops)
    // The Fortran uses the same variable names for all coarse levels
    // We'll use the same approach
    int nblf2 = 0, nblz2 = 0, nblcc2 = 0, igridc2 = 0, ipass2 = 0;
    int jc2 = 0, kc2 = 0, ic2 = 0;
    int lq1c2 = 0, lqc2 = 0, lvolc2 = 0, lqrc2 = 0, lvtc2 = 0, lxtc2 = 0;
    int lxibc2 = 0, lqc0c2 = 0, ldqc0c2 = 0;
    int nblend2 = 0;
    // For the coarse-grid loop variables, use arrays indexed by level
    // Actually in Fortran they are scalar and reused - we'll do the same
    // Reset all to 0 and use them as scalars in the loops
    // (The Fortran code uses the same scalar names at each level)
    int nblf_l = 0, nblz_l = 0, nblcc_l = 0, igridc_l = 0, ipass_l = 0;
    int jc_l = 0, kc_l = 0, ic_l = 0;
    int lq1c_l = 0, lqc_l = 0, lvolc_l = 0, lqrc_l = 0, lvtc_l = 0, lxtc_l = 0;
    int lxibc_l = 0, lqc0c_l = 0, ldqc0c_l = 0;
    int nblend_l = 0;
    int nblf_u = 0, nblz_u = 0, nblcc_u = 0, igridc_u = 0, ipass_u = 0;
    int jc_u = 0, kc_u = 0, ic_u = 0;
    int lq1c_u = 0, lqc_u = 0, lvolc_u = 0, lqrc_u = 0, lvtc_u = 0, lxtc_u = 0;
    int lxibc_u = 0, lqc0c_u = 0, ldqc0c_u = 0;
    int nblend_u = 0;
    // Use a single set of coarse-grid variables (Fortran reuses same names)
    // These are the "current" coarse grid variables used in the loops
    int nblf_cur = 0, nblz_cur = 0, nblcc_cur = 0, igridc_cur = 0, ipass_cur = 0;
    int jc_cur = 0, kc_cur = 0, ic_cur = 0;
    int lq1c_cur = 0, lqc_cur = 0, lvolc_cur = 0, lqrc_cur = 0, lvtc_cur = 0, lxtc_cur = 0;
    int lxibc_cur = 0, lqc0c_cur = 0, ldqc0c_cur = 0;
    int nblend_cur = 0;
    // Actually, looking at the Fortran more carefully, the coarse-grid variables
    // are just scalars reused in the downward/upward sweeps. Let's just use
    // the same names as in Fortran (nblf, nblz, nblcc, igridc, ipass, jc, kc, ic,
    // lq1c, lqc, lvolc, lqrc, lvtc, lxtc, lxibc, lqc0c, ldqc0c, nblend)
    // which are already declared above. The _cur/_l/_u variants are not needed.
    // We'll use the base names throughout.

    // ===== BEGIN FUNCTION BODY =====
    // Fortran lines 388-400
    levt_c = cmn_info.levelt[iseq-1];  // levelt(iseq)
    levb   = cmn_info.levelb[iseq-1];  // levelb(iseq)
    nitfo  = cmn_info.nitfo1[iseq-1];  // nitfo1(iseq)
    ncyc   = cmn_info.ncyc1[iseq-1];  // ncyc1(iseq)

    stop_me = false;

    if (levt_c - levb + 1 > 5) {
        fprintf(fortran_get_unit(11), "  stopping - too many levels in multigrid - max = 5\n");
         termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
    }

    lglobal = lfgm - (mseq - iseq);
    nembed  = levt_c - lglobal;

    // find first (starting) block at global level, nblstag
    nblstag = 0;
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (mblk2nd(nbl) == myid || myid == myhost) {
            if (lglobal != levelg(nbl)) continue;
            nblstag = nbl;
            break;
        }
    }

    // find last (ending) block at global level, nblendg
    nblendg = 0;
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (mblk2nd(nbl) == myid || myid == myhost) {
            if (lglobal != levelg(nbl)) continue;
            nblendg = nbl;
        }
    }

    // find first (starting) block at top level, nblstat
    nblstat = 0;
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (mblk2nd(nbl) == myid || myid == myhost) {
            if (levt_c != levelg(nbl)) continue;
            nblstat = nbl;
            break;
        }
    }

    nblstat_h = nblstat;

    // find last (ending) block at top level, nblendt
    nblendt = 0;
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (mblk2nd(nbl) == myid || myid == myhost) {
            if (levt_c != levelg(nbl)) continue;
            nblendt = nbl;
        }
    }

    // initialize mismatch check arrays
    for (nn = 1; nn <= 2*mxbli; nn++) geom_miss(nn) = 0.0;
    for (nn = 1; nn <= lbcprd; nn++) period_miss(nn) = 0.0;

    // begin time advancement
    ncall = 1;
    if ((float)dt > 0.0f && movie < 0) {
        isklton = 0;
        lhdr    = 0;
        if (iunst > 0) {
            // calculate grid speeds and metrics for use in boundary conditions
            // but do not change grid position (iupdat=0)
            for (nbl = 1; nbl <= nblock; nbl++) {
                timesave(nbl) = time2(nbl);
                time2(nbl)    = time2(nbl) - (double)dt;
            }
            iupdat   = 0;
            level    = cmn_info.levelt[iseq-1];
            irigb0   = 0;
            irbtrim0 = 0;
            updateg_ns::updateg(lw, lw2, w, mgwk, wk, nwork, iupdat, iseq, maxbl,
                maxgr, maxseg, nbci0, nbcj0, nbck0, nbcidim,
                nbcjdim, nbckdim, ibcinfo, jbcinfo, kbcinfo,
                nblock, levelg, igridg, utrans, vtrans, wtrans,
                omegax, omegay, omegaz, xorig, yorig, zorig,
                thetax, thetay, thetaz, rfreqt, rfreqr, xorig0,
                yorig0, zorig0, time2, thetaxl, thetayl, thetazl,
                itrans, irotat, idefrm, ncgg, iadvance, nou,
                bou, nbuf, ibufdim, myid, myhost, mycomm, mblk2nd,
                irigb0, irbtrim0, nt);
            for (nbl = 1; nbl <= nblock; nbl++) {
                time2(nbl) = timesave(nbl);
            }
            // get interpolation data for dynamic boundaries (for qout)
            iwk1   = 1;
            iwk2   = iwk1 + 65*maxbl;
            iwk3   = iwk2 + 3*intmx;
            iwork1 = iwork - iwk3;
            if (iwork1 < 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine dynptch");
                 termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
            {
                FortranArray1DRef<int> iwk_iwk1_1d = iwk.slice(iwk1, iwk.ubound(1));
                FortranArray1DRef<int> iwk_iwk2_1d = iwk.slice(iwk2, iwk.ubound(1));
                FortranArray2DRef<int> iwk_iwk1(iwk_iwk1_1d.data(), 65, maxbl);
                FortranArray2DRef<int> iwk_iwk2(iwk_iwk2_1d.data(), 65, maxbl);
                dynptch_ns::dynptch(lw, lw2, w, mgwk, wk, nwork, ncall,
                    maxgr, maxbl, msub1, intmx, mxxe, mptch, jdimg,
                    kdimg, idimg, xorig, yorig, zorig, nblock, ngrid,
                    levelg, ncgg, nblg, windex, ninter, iindex, nblkpt,
                    windx, nintr, iindx, mblkpt, llimit, iitmax,
                    mmcxie, mmceta, ncheck, iifit, iic0,
                    iiorph, iitoss, ifiner, factjlo, factjhi,
                    factklo, factkhi, dx, dy, dz, dthetx, dthety,
                    dthetz, dthetxx, dthetyy, dthetzz,
                    isav_dpat, isav_dpat_b, intmax, maxxe, nsub1,
                    iwk_iwk1, iwk_iwk2_1d, lout, ifrom, xif1, xif2, etf1,
                    etf2, jjmax1, kkmax1, iiint1, iiint2, nblk1,
                    nblk2, jimage, kimage, jte, kte, jmm, kmm,
                    xte, yte, zte, xmi, ymi, zmi, xmie, ymie,
                    zmie, sxie, seta, sxie2, seta2, xie2s,
                    eta2s, temp, x2, y2, z2, x1, y1, z1,
                    myid, myhost, mycomm, mblk2nd, nou, bou, nbuf,
                    ibufdim, igridg, iemg);
            }
            ncall = ncall + 1;
        }
        iwk1   = 1;
        iwk2   = iwk1 + 65*maxbl;
        iwk3   = iwk2 + nsub1;
        iwk4   = iwk3 + nsub1;
        iwk5   = iwk4 + nsub1;
        iwk6   = iwk5 + nsub1;
        iwork1 = iwork - iwk6;
        if (iwork1 < 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine pre_bc");
             termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
        }
        for (iii = 1; iii <= iwk6; iii++) iwk(iii) = 0;
        {
            FortranArray2DRef<int> iwk_iwk1_2d(iwk.data()+iwk1-1, 65, maxbl);
            FortranArray1DRef<int> iwk_iwk2 = iwk.slice(iwk2, iwk.ubound(1));
            FortranArray1DRef<int> iwk_iwk3 = iwk.slice(iwk3, iwk.ubound(1));
            FortranArray1DRef<int> iwk_iwk4 = iwk.slice(iwk4, iwk.ubound(1));
            FortranArray1DRef<int> iwk_iwk5 = iwk.slice(iwk5, iwk.ubound(1));
            pre_bc_ns::pre_bc(lw, lw2, iwk_iwk1_2d, maxbl, maxgr, maxseg, ninter, intmax,
                nsub1, iindex, isav_pat, iwk_iwk2,
                iwk_iwk3, iwk_iwk4, iwk_iwk5, mxbli, nbli, limblk,
                isva, nblon, nblk, lbcprd, isav_prd, bcvali, bcvalj,
                bcvalk, nblg, lbcemb, iemg, igridg, isav_emb,
                iviscg, jdimg, kdimg, idimg, nbci0, nbcj0, nbck0,
                nbcjdim, nbckdim, nbcidim, ibcinfo, jbcinfo,
                kbcinfo, iadvance, myid, myhost, mycomm, mblk2nd,
                nou, bou, nbuf, ibufdim, isav_pat_b, levelg, nblock,
                isav_blk);
        }
        iwk1   = 1;
        iwk2   = iwk1 + 3*nplots;
        iwk3   = iwk2 + maxbl;
        iwork1 = iwork - iwk3 + 1;
        if (iwork1 <= 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine qout");
             termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
        }
        for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
        for (iii = 1; iii <= nwork; iii++) wk(iii) = 0.0;
        {
            FortranArray2DRef<int> iwk_ip3ddim(iwk.data()+iwk1-1, 3, nplots);
            FortranArray1DRef<int> iwk_nmap = iwk.slice(iwk2, iwk.ubound(1));
            FortranArray1DRef<int> iwk_iwk3 = iwk.slice(iwk3, iwk.ubound(1));
            qout_ns::qout(iseq, lw, lw2, w, mgwk, wk, nwork,
                nplots, iovrlp, iibg, kkbg, jjbg, ibcg, lbg,
                ibpntsg, qb, lwdat, nbci0, nbcj0, nbck0,
                nbcidim, nbcjdim, nbckdim, jbcinfo, kbcinfo,
                ibcinfo, bcfilei, bcfilej, bcfilek, itrans,
                irotat, idefrm, nblock, levelg, igridg, iviscg, jdimg,
                kdimg, idimg, nblg, clw, ncycmax, nplot3d,
                inpl3d, ip3dsurf, nprint, inpr, iadvance, mycomm,
                myid, myhost, mblk2nd, nou, bou, nbuf, ibufdim, maxbl,
                maxgr, maxseg, iitot, jsg, ksg, isg, jeg, keg, ieg,
                ninter, windex, iindex, nblkpt, intmax, nsub1, maxxe,
                nblk, nbli, limblk, isva, nblon, mxbli, thetay,
                iwk_ip3ddim, iwk_nmap, iwk_iwk3, iwork1, xorig, yorig,
                zorig, period_miss, geom_miss, epsc0, epsrot,
                isav_blk, isav_pat, isav_pat_b, isav_emb, isav_prd,
                lbcprd, lbcemb, dthetxx, dthetyy, dthetzz, nblcg, lfgm,
                istat2_bl, istat2_pa, istat2_pe, istat2_em,
                istat_size, vormax, ivmax, jvmax, kvmax, bcfiles,
                mxbcfil, iprnsurf, nummem);
        }
        for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
        for (iii = 1; iii <= nwork; iii++) wk(iii) = 0.0;
        nframes = nframes + 1;
    }


    // Fortran lines 605-623: output header
    if (myid == myhost) {
        fprintf(fortran_get_unit(11), "\n***** BEGINNING TIME ADVANCEMENT, iseq =%2d *****\n\n", iseq);
        if ((float)dt < 0.0f) {
            fprintf(fortran_get_unit(11), " steady-state computations\n");
        } else {
            if (std::abs(ita) == 1) {
                fprintf(fortran_get_unit(11), " time-accurate computations - 1st order in time, number of time steps: %4d\n", ntstep);
            }
            if (std::abs(ita) == 2) {
                fprintf(fortran_get_unit(11), " time-accurate computations - 2nd order in time, number of time steps: %4d\n", ntstep);
            }
            if (ita < 0) {
                fprintf(fortran_get_unit(11), " pseudo-time term included on LHS for stability with large dt\n");
            }
        }
        if (cprec != 0.0f) {
            fprintf(fortran_get_unit(11), "\n low Mach number preconditioning used, with:\n   cprec, uref, avn =  %7.4f, %7.4f, %7.4f\n", (double)cprec, (double)uref, (double)avn);
        }
    }

    isklton = 1;
    ioutsub = ncyc;

    // loop over time steps: do 6000 nt=1,ntstep
    for (nt = 1; nt <= ntstep; nt++) {
    // label 6000 (loop start)

    // check for user stop file
    if (myid == myhost) {
        stop_me = fortran_file_exists("stop");
        if (stop_me) {
            ntstep = nt;
            for (is = iseq+1; is <= mseq; is++) {
                cmn_info.ncyc1[is-1] = 0;
            }
            fprintf(fortran_get_unit(11), "\n");
            fprintf(fortran_get_unit(11), "Stopping: user-invoked stop file detected at time step %5d\n", ntstep);
            fprintf(fortran_get_unit(11), "\n");
        }
    }

    itatemp = ita;
    // always do 1st order time 1st iteration if restarting and changing time step
    if (nt == 1 && dt != dtold) {
        ita = (ita >= 0) ? 1 : -1;
    }

    if (nt > 1) isklton = 0;

    // start on top level
    level = cmn_info.levelt[iseq-1];

    // aeroelastic predictor step
    if (iunst > 1 && naesrf > 0) {
        ae_pred_ns::ae_pred(aesrfdat, stm, stmi, gforcn, gforcnm, xs, xxn,
            x0, perturb, cmyw(ntt), clw(ntt), xorgrb, yorgrb,
            zorgrb, nmds, maxaes, irbtrim, maxbl, myid);
        // rigid body trim predictor step
        if (irbtrim > 0) irigb = 0;
        if (irigb > 0) {
            rb_pred_ns::rb_pred(nt);
        }
        if (irigb == 1 || irbtrim == 1) {
            for (nbl = 1; nbl <= nblock; nbl++) {
                thetayl(nbl) = thetay(nbl);
                thetay(nbl)  = (double)cmn_trim.alf1;
                omegay(nbl)  = (thetay(nbl) - thetayl(nbl)) / (double)dt;
            }
        }
    }

    // update grid position, speed, and dynamic patch interpolation coefficients
    if (iunst > 0) {
        // store off old grid position for 2nd order (in time) update
        if (meshdef == 1) {
            if (myid == myhost) {
                fprintf(fortran_get_unit(11), " Deforming mesh and calculating volumes and metrics, time step %5d\n", nt);
            }
        }
        for (nbl = 1; nbl <= nblock; nbl++) {
            if (myid == mblk2nd(nbl) && (levelg(nbl) >= lglobal && levelg(nbl) <= cmn_info.levelt[iseq-1])) {
                if (idefrm(nbl) > 0) {
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    mdim = jdim*kdim*idim;
                    for (lll = 1; lll <= mdim; lll++) {
                        w(lxnm1+lll-1) = w(lx+lll-1);
                        w(lynm1+lll-1) = w(ly+lll-1);
                        w(lznm1+lll-1) = w(lz+lll-1);
                    }
                }
            }
        }
        iupdat = 1;
        updateg_ns::updateg(lw, lw2, w, mgwk, wk, nwork, iupdat, iseq, maxbl,
            maxgr, maxseg, nbci0, nbcj0, nbck0, nbcidim,
            nbcjdim, nbckdim, ibcinfo, jbcinfo, kbcinfo,
            nblock, levelg, igridg, utrans, vtrans, wtrans,
            omegax, omegay, omegaz, xorig, yorig, zorig,
            thetax, thetay, thetaz, rfreqt, rfreqr, xorig0,
            yorig0, zorig0, time2, thetaxl, thetayl, thetazl,
            itrans, irotat, idefrm, ncgg, iadvance, nou,
            bou, nbuf, ibufdim, myid, myhost, mycomm, mblk2nd,
            irigb, irbtrim, nt);
        // reset grid position if grid has been translated/rotated past limits
        iwk1   = maxbl;
        iwork1 = iwork - iwk1;
        if (iwork1 < 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine resetg");
             termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
        }
        for (iii = 1; iii <= iwk1; iii++) iwk(iii) = 0;
        resetg_ns::resetg(lw, lw2, w, mgwk, wk, nwork, iupdat, iseq, maxbl,
            maxgr, maxseg, nbci0, nbcj0, nbck0, nbcidim,
            nbcjdim, nbckdim, ibcinfo, jbcinfo, kbcinfo,
            nblock, levelg, igridg, utrans, vtrans, wtrans,
            omegax, omegay, omegaz, xorig, yorig, zorig,
            thetax, thetay, thetaz, rfreqt, rfreqr, xorig0,
            yorig0, zorig0, time2, thetaxl, thetayl, thetazl,
            itrans, irotat, idefrm, ncgg, iadvance,
            dxmx, dymx, dzmx, dthymx, dthzmx, dthxmx,
            iitot, iovrlp, lig, lbg, iipntsg, ibpntsg,
            qb, iibg, kkbg, jjbg, ibcg, nou, bou, nbuf, ibufdim,
            myid, myhost, mycomm, mblk2nd, iwk, nt);
        if (iunst > 1) {
            iwk1   = 2*maxbl;
            iwork1 = iwork - iwk1;
            if (iwork1 < 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine updatedg");
                 termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
            }
            for (iii = 1; iii <= iwk1; iii++) iwk(iii) = 0;
            iupdat = 1;
            {
                FortranArray1DRef<int> iwk_maxbl1 = iwk.slice(maxbl+1, iwk.ubound(1));
                updatedg_ns::updatedg(lw, lw2, w, mgwk, wk, nwork, iupdat, iseq, maxbl,
                    maxgr, maxseg, nbci0, nbcj0, nbck0, nbcidim,
                    nbcjdim, nbckdim, ibcinfo, jbcinfo, kbcinfo,
                    nblock, levelg, igridg, idefrm, ncgg, iadvance, nou,
                    bou, nbuf, ibufdim, myid, myhost, mycomm, mblk2nd,
                    utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae,
                    xorgae, yorgae, zorgae, thtxae, thtyae, thtzae,
                    rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf,
                    kcsi, kcsf, freq, gmass, damp, x0, gf0, nmds, maxaes,
                    aesrfdat, perturb, itrans, irotat, islavept,
                    nslave, iskip, jskip, kskip, xs, xxn, nsegdfrm,
                    idfrmseg, iaesurf, maxsegdg, iwk, nmaster, nt,
                    xorig, yorig, zorig, xorgae0, yorgae0, zorgae0,
                    icouple, iwk_maxbl1, nnodes, nblelst, iskmax,
                    jskmax, kskmax, ue);
            }
        }
    }


    // Fortran lines 799-832: dynptch call if iunst > 0
    if (iunst > 0) {
        iwk1   = 1;
        iwk2   = iwk1 + 65*maxbl;
        iwk3   = iwk2 + 3*intmx;
        iwork1 = iwork - iwk3;
        if (iwork1 < 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine dynptch");
             termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
        }
        for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
        {
            FortranArray2DRef<int> iwk_iwk1_2d(iwk.data()+iwk1-1, 65, maxbl);
            FortranArray1DRef<int> iwk_iwk2 = iwk.slice(iwk2, iwk.ubound(1));
            dynptch_ns::dynptch(lw, lw2, w, mgwk, wk, nwork, ncall,
                maxgr, maxbl, msub1, intmx, mxxe, mptch, jdimg,
                kdimg, idimg, xorig, yorig, zorig, nblock, ngrid,
                levelg, ncgg, nblg, windex, ninter, iindex, nblkpt,
                windx, nintr, iindx, mblkpt, llimit, iitmax,
                mmcxie, mmceta, ncheck, iifit, iic0,
                iiorph, iitoss, ifiner, factjlo, factjhi,
                factklo, factkhi, dx, dy, dz, dthetx, dthety,
                dthetz, dthetxx, dthetyy, dthetzz,
                isav_dpat, isav_dpat_b, intmax, maxxe, nsub1,
                iwk_iwk1_2d, iwk_iwk2, lout, ifrom, xif1, xif2, etf1,
                etf2, jjmax1, kkmax1, iiint1, iiint2, nblk1,
                nblk2, jimage, kimage, jte, kte, jmm, kmm,
                xte, yte, zte, xmi, ymi, zmi, xmie, ymie,
                zmie, sxie, seta, sxie2, seta2, xie2s,
                eta2s, temp, x2, y2, z2, x1, y1, z1,
                myid, myhost, mycomm, mblk2nd, nou, bou, nbuf,
                ibufdim, igridg, iemg);
        }
        ncall = ncall + 1;
    }

    // initial setup for exchange of bc data between blocks/processors
    if (nt == 1) {
        iwk1   = 1;
        iwk2   = iwk1 + 65*maxbl;
        iwk3   = iwk2 + nsub1;
        iwk4   = iwk3 + nsub1;
        iwk5   = iwk4 + nsub1;
        iwk6   = iwk5 + nsub1;
        iwork1 = iwork - iwk6;
        if (iwork1 < 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine pre_bc");
             termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
        }
        for (iii = 1; iii <= iwk6; iii++) iwk(iii) = 0;
        {
            FortranArray2DRef<int> iwk_iwk1(iwk.data()+iwk1-1, 65, maxbl);
            FortranArray1DRef<int> iwk_iwk2 = iwk.slice(iwk2, iwk.ubound(1));
            FortranArray1DRef<int> iwk_iwk3 = iwk.slice(iwk3, iwk.ubound(1));
            FortranArray1DRef<int> iwk_iwk4 = iwk.slice(iwk4, iwk.ubound(1));
            FortranArray1DRef<int> iwk_iwk5 = iwk.slice(iwk5, iwk.ubound(1));
            pre_bc_ns::pre_bc(lw, lw2, iwk_iwk1, maxbl, maxgr, maxseg, ninter, intmax,
                nsub1, iindex, isav_pat, iwk_iwk2,
                iwk_iwk3, iwk_iwk4, iwk_iwk5, mxbli, nbli, limblk,
                isva, nblon, nblk, lbcprd, isav_prd,
                bcvali, bcvalj, bcvalk, nblg, lbcemb,
                iemg, igridg, isav_emb,
                iviscg, jdimg, kdimg, idimg, nbci0, nbcj0, nbck0,
                nbcjdim, nbckdim, nbcidim, ibcinfo, jbcinfo,
                kbcinfo, iadvance, myid, myhost, mycomm, mblk2nd,
                nou, bou, nbuf, ibufdim, isav_pat_b, levelg, nblock,
                isav_blk);
        }
    }

    // begin multigrid cycles: do 7000 icyc=1,ncyc
    for (icyc = 1; icyc <= ncyc; icyc++) {
    // label 7000 (loop start)

    if (icyc > 1) isklton = 0;

    // check for user stop file (steady state only)
    if (myid == myhost) {
        if ((float)dt < 0.0f) {
            stop_me = fortran_file_exists("stop");
            if (stop_me) {
                ncyc = icyc;
                for (is = iseq+1; is <= mseq; is++) {
                    cmn_info.ncyc1[is-1] = 0;
                }
                fprintf(fortran_get_unit(11), "\n");
                fprintf(fortran_get_unit(11), "Stopping: user-invoked stop file detected at cycle %5d\n", icyc);
                fprintf(fortran_get_unit(11), "\n");
            }
        }
    }

    if (isklton == 1 && myid == myhost) {
        fprintf(fortran_get_unit(11), "\n***** BEGINNING MULTIGRID CYCLE *****\n\n");
        fprintf(fortran_get_unit(11), " iseq=%5d\n", iseq);
        fprintf(fortran_get_unit(11), " level top =%3d\n", levt_c);
        fprintf(fortran_get_unit(11), " level bottom =%3d\n", levb);
        if (nembed > 0) fprintf(fortran_get_unit(11), " number of embedded grid levels =%3d\n", nembed);
        {
            int lggl = levt_c - levb - nembed + 1;
            fprintf(fortran_get_unit(11), " number of global grid levels =%3d\n", lggl);
        }
        fprintf(fortran_get_unit(11), " lglobal=%3d\n", lglobal);
    }

    for (int i = 1; i <= level - levb + 1; i++) {
        cmn_info.nsm[i-1] = ngam;
    }

    iaccnt = 1;

    nit = mit(level - levb + 1, iseq);
    if (nit > 20) {
        fprintf(fortran_get_unit(11), " stopping...increase dimension of nptsr\n");
         termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
    }

    kxpand = 1;
    kode   = 1;
    mode   = 0;

    // label 9000: begin/continue downward leg of multigrid cycle
    label_9000:;

    // find last (ending) block on current level, nblend
    nblend = 0;
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (level != levelg(nbl)) continue;
        nblend = nbl;
    }

    // zero out residual and force coefficients
    if (level == levt_c) {
        ntt0 = ntt;
        if (((float)dt < 0.0f && nt == 1) ||
            ((float)dt >= 0.0f && icyc == ioutsub)) {
            ntt = ntt + 1;
            rms(ntt)    = 0.0;
            clw(ntt)    = 0.0;
            cdw(ntt)    = 0.0;
            cyw(ntt)    = 0.0;
            for (int mm = 1; mm <= nummem; mm++) rmstr(ntt, mm) = 0.0;
            for (int mm = 1; mm <= nummem; mm++) nneg(ntt, mm)  = 0;
            chdw(ntt)   = 0.0;
            swetw(ntt)  = 0.0;
            cmxw(ntt)   = 0.0;
            cmyw(ntt)   = 0.0;
            cmzw(ntt)   = 0.0;
            for (int ic2 = 1; ic2 <= n_clcd; ic2++)
                for (int jc2 = 1; jc2 <= 2; jc2++)
                    clcd(jc2, ic2, ntt) = 0.0;
            cxw(ntt)    = 0.0;
            czw(ntt)    = 0.0;
            cdpw(ntt)   = 0.0;
            cdvw(ntt)   = 0.0;
            fmdotw(ntt) = 0.0;
            cfttotw(ntt)= 0.0;
            cftmomw(ntt)= 0.0;
            cftpw(ntt)  = 0.0;
            cftvw(ntt)  = 0.0;
        }
        for (ntime = 1; ntime <= nit; ntime++) {
            nptsr[ntime-1]   = 0;
            nptsr[ntime]     = 0;
        }
        nptsrb = 0;
        for (int mm = 1; mm <= nummem; mm++) nnegb_arr(mm) = 0.0;
        for (int mm = 1; mm <= nummem; mm++) rmstb_arr(mm) = 0.0;
        rmsb = 0.0;
        for (int mm = 1; mm <= nummem; mm++) sumn_arr(mm) = 0.0;
        for (int mm = 1; mm <= nummem; mm++) negn_arr(mm) = 0.0;
    }

    // meshdef = 1 if flow solution is to be bypassed
    if (meshdef == 1) goto label_6998;


    // Fortran lines 1005-1013
    nit1 = nit + std::min(1, mgflag);
    if (level == levb) nit1 = nit;
    if (level > lglobal) nit1 = nit + 1;

    if (isklton > 0 && myid == myhost) {
        fprintf(fortran_get_unit(11), "\n***** BEGINNING RESIDUAL/RESTRICTION LOOP, level, kode, mode =%3d%3d%3d\n", level, kode, mode);
    }
    ilc = 0;

    // do 6500 ntime=1,nit1
    for (ntime = 1; ntime <= nit1; ntime++) {

    // initialize density/pressure bcs for comparison in bcchk
    if (isklton > 0 && ntime == 1) {
        for (nbl = 1; nbl <= nblock; nbl++) {
            if (iadvance(nbl) < 0) continue;
            if (level != levelg(nbl)) continue;
            if (mblk2nd(nbl) == myid) {
                lead_ns::lead(nbl, lw, lw2, maxbl);
                ibcflg = 0;
                istop  = 1;
                {
                    FortranArray4DRef<double> q_bc(w.data()+lq-1, jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qi0_bc(w.data()+lqi0-1, jdim, kdim, 5, 4);
                    FortranArray4DRef<double> qj0_bc(w.data()+lqj0-1, kdim, idim-1, 5, 4);
                    FortranArray4DRef<double> qk0_bc(w.data()+lqk0-1, jdim, idim-1, 5, 4);
                    FortranArray3DRef<double> blank_bc(w.data()+lblk-1, jdim, kdim, idim);
                    bcchk_ns::bcchk(idim, jdim, kdim, q_bc, qi0_bc, qj0_bc, qk0_bc,
                        blank_bc, ibcflg, nbl, nou, bou, nbuf, ibufdim, myid, istop,
                        igridg, maxbl);
                }
            } // end if mblk2nd == myid
        } // end do 6990 nbl
    } // end if isklton > 0 && ntime == 1

    // update physical boundary conditions
    if (iipv > 0) {
        if (level >= lglobal && ntime == nit) {
            nttuse = std::max(ntt-1, 1);
            clwuse = clw(nttuse);
        }
    }
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (iadvance(nbl) < 0) continue;
        if (level != levelg(nbl)) continue;
        if (mblk2nd(nbl) == myid) {
            lead_ns::lead(nbl, lw, lw2, maxbl);
            lres  = 1;
            need  = nwork - lres + 1;
            {
                FortranArray1DRef<double> wk_lres = wk.slice(lres, wk.ubound(1));
                bc_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk_lres, need, clwuse,
                    nou, bou, nbuf, ibufdim, maxbl, maxgr, maxseg, itrans,
                    irotat, idefrm, igridg, nblg, nbci0, nbcj0, nbck0,
                    nbcidim, nbcjdim, nbckdim, ibcinfo, jbcinfo,
                    kbcinfo, bcfilei, bcfilej, bcfilek, lwdat, myid,
                    idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
            }
        }
    }


    // update periodic boundary conditions
    lres   = 1;
    need   = nwork - lres + 1;
    mneed  = lbcprd*5;
    iwk1   = 1;
    iwk2   = iwk1 + mneed;
    iwk3   = iwk2 + mneed;
    iwk4   = iwk3 + mneed;
    iwk5   = iwk4 + mneed*2;
    iwk6   = iwk5 + mneed;
    iwork1 = iwork - iwk6;
    if (iwork1 < 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine bc_period");
         termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
    }
    for (iii = 1; iii <= iwk6; iii++) iwk(iii) = 0;
    {
        FortranArray1DRef<double> wk_lres = wk.slice(lres, wk.ubound(1));
        FortranArray1DRef<int> iwk_iwk1(&iwk(iwk1), mneed, 1);
        FortranArray1DRef<int> iwk_iwk2(&iwk(iwk2), mneed, 1);
        FortranArray1DRef<int> iwk_iwk3(&iwk(iwk3), mneed, 1);
        FortranArray2DRef<int> iwk_iwk4(&iwk(iwk4), lbcprd, 10);
        FortranArray1DRef<int> iwk_iwk5(&iwk(iwk5), mneed, 1);
        mgblk_ns::bc_period(ntime, nbl, lw, lw2, w, mgwk, wk_lres, need, maxbl,
            maxgr, maxseg, iadvance, bcfilei, bcfilej, bcfilek,
            lwdat, xorig, yorig, zorig, jdimg, kdimg, idimg, lbcprd,
            isav_prd, period_miss, epsrot, iwk_iwk1, iwk_iwk2,
            iwk_iwk3, iwk_iwk4, iwk_iwk5, myid, myhost,
            mycomm, mblk2nd, nou, bou, nbuf, ibufdim,
            istat2_pe, istat_size, bcfiles, mxbcfil, nummem);
    }


    // update embedded-grid boundary conditions
    lres   = 1;
    need   = nwork - lres + 1;
    mneed  = lbcemb*3;
    iwk1   = 1;
    iwk2   = iwk1 + mneed;
    iwk3   = iwk2 + mneed;
    iwk4   = iwk3 + mneed;
    iwk5   = iwk4 + mneed*2;
    iwk6   = iwk5 + mneed;
    iwork1 = iwork - iwk6;
    if (iwork1 < 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine bc_embed");
         termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
    }
    for (iii = 1; iii <= iwk6; iii++) iwk(iii) = 0;
    {
        FortranArray1DRef<double> wk_lres = wk.slice(lres, wk.ubound(1));
        FortranArray1DRef<int> iwk_iwk1(&iwk(iwk1), mneed, 1);
        FortranArray1DRef<int> iwk_iwk2(&iwk(iwk2), mneed, 1);
        FortranArray1DRef<int> iwk_iwk3(&iwk(iwk3), mneed, 1);
        FortranArray2DRef<int> iwk_iwk4(&iwk(iwk4), lbcemb, 6);
        FortranArray1DRef<int> iwk_iwk5(&iwk(iwk5), mneed, 1);
        mgblk_ns::bc_embed(ntime, nbl, lw, lw2, w, mgwk, wk_lres, need, maxbl,
            maxgr, lbcemb, iadvance, idimg, jdimg, kdimg,
            isav_emb, iwk_iwk1, iwk_iwk2,
            iwk_iwk3, iwk_iwk4, iwk_iwk5, myid, myhost,
            mycomm, mblk2nd, nou, bou, nbuf, ibufdim, iviscg,
            istat2_em, istat_size, nummem);
    }

    // update 1-1 block boundary conditions
    lres   = 1;
    need   = nwork - lres + 1;
    mneed  = mxbli*5;
    iwk1   = 1;
    iwk2   = iwk1 + mneed;
    iwk3   = iwk2 + mneed;
    iwk4   = iwk3 + mneed;
    iwk5   = iwk4 + mneed*2;
    iwk6   = iwk5 + mneed;
    iwork1 = iwork - iwk6;
    if (iwork1 < 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine bc_blkint");
         termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
    }
    for (iii = 1; iii <= iwk6; iii++) iwk(iii) = 0;
    {
        FortranArray1DRef<double> wk_lres = wk.slice(lres, wk.ubound(1));
        FortranArray1DRef<int> iwk_iwk1(&iwk(iwk1), mneed, 1);
        FortranArray1DRef<int> iwk_iwk2(&iwk(iwk2), mneed, 1);
        FortranArray1DRef<int> iwk_iwk3(&iwk(iwk3), mneed, 1);
        FortranArray2DRef<int> iwk_iwk4(&iwk(iwk4), mxbli, 10);
        FortranArray1DRef<int> iwk_iwk5(&iwk(iwk5), mneed, 1);
        mgblk_ns::bc_blkint(ntime, nbl, lw, lw2, w, mgwk, wk_lres, need, maxbl,
            maxgr, mxbli, iadvance, geom_miss, epsc0, nblk, nbli,
            limblk, isva, nblon, jdimg, kdimg, idimg,
            mblk2nd, isav_blk, iwk_iwk1,
            iwk_iwk2, iwk_iwk3, iwk_iwk4, iwk_iwk5,
            nou, bou, nbuf, ibufdim, myid, myhost, mycomm,
            istat2_bl, istat_size, nummem);
    }

    // update patch-grid boundary conditions
    lres   = 1;
    need   = nwork - lres + 1;
    mneed  = intmax*nsub1*3;
    iwk1   = 1;
    iwk2   = iwk1 + mneed;
    iwk3   = iwk2 + mneed;
    iwk4   = iwk3 + mneed;
    iwk5   = iwk4 + mneed*2;
    iwk6   = iwk5 + mneed*2;
    iwork1 = iwork - iwk6;
    if (iwork1 < 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine bc_patch");
         termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
    }
    for (iii = 1; iii <= iwk6; iii++) iwk(iii) = 0;
    {
        FortranArray1DRef<double> wk_lres = wk.slice(lres, wk.ubound(1));
        FortranArray1DRef<int> iwk_iwk1(&iwk(iwk1), mneed, 1);
        FortranArray1DRef<int> iwk_iwk2(&iwk(iwk2), mneed, 1);
        FortranArray1DRef<int> iwk_iwk3(&iwk(iwk3), mneed, 1);
        FortranArray3DRef<int> iwk_iwk4(&iwk(iwk4), intmax, nsub1, 6);
        FortranArray2DRef<int> iwk_iwk5(&iwk(iwk5), intmax*nsub1*3, 2);
        mgblk_ns::bc_patch(ntime, nbl, lw, lw2, w, mgwk, wk_lres, need, maxbl,
            maxgr, intmax, nsub1, maxxe, iadvance, jdimg, kdimg,
            idimg, ninter, windex, iindex, nblkpt, dthetxx,
            dthetyy, dthetzz, isav_pat, isav_pat_b,
            iwk_iwk1, iwk_iwk2, iwk_iwk3,
            iwk_iwk4, iwk_iwk5, myid, myhost, mycomm,
            mblk2nd, nou, bou, nbuf, ibufdim,
            istat2_pa, istat_size, nummem);
    }

    // update chimera boundary conditions
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (iadvance(nbl) < 0) continue;
        if (level != levelg(nbl)) continue;
        if (mblk2nd(nbl) == myid) {
            lead_ns::lead(nbl, lw, lw2, maxbl);
            lres  = 1;
            need  = nwork - lres + 1;
            {
                FortranArray1DRef<double> wk_lres = wk.slice(lres, wk.ubound(1));
                int int_updt_tmp = 1;
                bc_xmera_ns::bc_xmera(ntime, nbl, lw, lw2, w, mgwk, wk_lres, need, maxbl,
                    iitot, iviscg, iovrlp, lbg, ibpntsg, qb, iibg, kkbg,
                    jjbg, ibcg, nou, bou, nbuf, ibufdim, int_updt_tmp, nummem);
            }
        }
    }


    // output bc info to main output file
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (level != levelg(nbl)) continue;
        if (mblk2nd(nbl) == myid) {
            lead_ns::lead(nbl, lw, lw2, maxbl);
            lres  = 1;
            need  = nwork - lres + 1;
            {
                FortranArray1DRef<double> wk_lres = wk.slice(lres, wk.ubound(1));
                bc_info_ns::bc_info(ntime, nbl, w, mgwk, ibcinfo, jbcinfo,
                    kbcinfo, nbci0, nbcj0, nbck0, nbcidim, nbcjdim,
                    nbckdim, bcfilei, bcfilej, bcfilek, igridg, itrans,
                    irotat, idefrm, idimg, jdimg, kdimg, nblg, iibg,
                    kkbg, jjbg, ibcg, ibpntsg, iipntsg, lig, lbg,
                    isav_blk, isav_prd, iemg, nbli,
                    geom_miss, period_miss, ninter, iindex, nou, bou, nbuf,
                    myid, mblk2nd, ibufdim, maxbl, maxseg, iitot,
                    intmax, nsub1, mxbli, lbcprd, epsc0, epsrot, lwdat,
                    maxgr, iovrlp, bcfiles, mxbcfil);
            }
        }
        if (isklton > 0) {
            writ_buf_ns::writ_buf(nbl, writ_buf_11, nou, bou, nbuf, ibufdim, myhost, myid,
                mycomm, mblk2nd, maxbl);
        }
    }

    // do 6989: bcchk and chkdef
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (iadvance(nbl) < 0) continue;
        if (level != levelg(nbl)) continue;
        if (mblk2nd(nbl) == myid) {
            lead_ns::lead(nbl, lw, lw2, maxbl);
            if (isklton > 0 && ntime == 1) {
                FortranArray4DRef<double> q_v(w.data()+lq-1, jdim, kdim, idim, 5);
                FortranArray4DRef<double> qi0_v(w.data()+lqi0-1, jdim, kdim, 5, 4);
                FortranArray4DRef<double> qj0_v(w.data()+lqj0-1, kdim, idim-1, 5, 4);
                FortranArray4DRef<double> qk0_v(w.data()+lqk0-1, jdim, idim-1, 5, 4);
                FortranArray3DRef<double> blank_v(w.data()+lblk-1, jdim, kdim, idim);
                ibcflg = 1; istop = 1;
                bcchk_ns::bcchk(idim, jdim, kdim, q_v, qi0_v, qj0_v, qk0_v,
                    blank_v, ibcflg, nbl, nou, bou, nbuf, ibufdim, myid, istop,
                    igridg, maxbl);
                ibcflg = 2; istop = 1;
                if (ichk == -1 || ichk > 1) istop = 0;
                bcchk_ns::bcchk(idim, jdim, kdim, q_v, qi0_v, qj0_v, qk0_v,
                    blank_v, ibcflg, nbl, nou, bou, nbuf, ibufdim, myid, istop,
                    igridg, maxbl);
                if (idefrm(nbl) == 1) {
                    FortranArray3DRef<double> bci_v(w.data()+lbci-1, jdim, kdim, 2);
                    FortranArray3DRef<double> bcj_v(w.data()+lbcj-1, kdim, idim-1, 2);
                    FortranArray3DRef<double> bck_v(w.data()+lbck-1, jdim, idim-1, 2);
                    chkdef_ns::chkdef(nbl, idim, jdim, kdim, bci_v, bcj_v, bck_v,
                        icsi, icsf, jcsi, jcsf, kcsi, kcsf, nsegdfrm,
                        idfrmseg, maxbl, maxsegdg, nou, bou, nbuf,
                        ibufdim, myid);
                }
            }
        }
        if (isklton > 0) {
            writ_buf_ns::writ_buf(nbl, writ_buf_11, nou, bou, nbuf, ibufdim, myhost, myid,
                mycomm, mblk2nd, maxbl);
        }
    }

    // cycle through blocks: update time-averaging counters
    if (iteravg == 1 || iteravg == 2) {
        if (level >= lglobal && ntime == nit) {
            if ((float)dt < 0.0f || ((float)dt > 0.0f && icyc == ncyc)) {
                xnumavg  = xnumavg  + 1.0f;
                xnumavg2 = xnumavg2 + 1.0f;
            }
        }
    }

    // do 8000 nbl=1,nblock
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (level != levelg(nbl)) continue;
        lead_ns::lead(nbl, lw, lw2, maxbl);

        if (mblk2nd(nbl) == myid) {
            iprerot = 1;
            // for time accurate calculations, compute and store dqc0 and qc0
            if ((float)dt > 0.0f && iadvance(nbl) >= 0) {
                if (icyc == 1 && level >= lglobal && ntime == 1) {
                    if (std::abs(ita) == 2) {
                        FortranArray4DRef<double> q_v(w.data()+lq-1, jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qc0_v(w.data()+lqc0-1, jdim, kdim, idim-1, 5);
                        FortranArray4DRef<double> dqc0_v(w.data()+ldqc0-1, jdim, kdim, idim-1, 5);
                        setdqc0_ns::setdqc0(jdim, kdim, idim, q_v, qc0_v, dqc0_v);
                    }
                    {
                        FortranArray4DRef<double> q_v(w.data()+lq-1, jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qc0_v(w.data()+lqc0-1, jdim, kdim, idim-1, 5);
                        setqc0_ns::setqc0(jdim, kdim, idim, q_v, qc0_v);
                    }
                    if (irotat(nbl) > 0 && iprerot > 0) {
                        dthtx = thetax(nbl) - thetaxl(nbl);
                        dthty = thetay(nbl) - thetayl(nbl);
                        dthtz = thetaz(nbl) - thetazl(nbl);
                        lqwk  = 1;
                        {
                            FortranArray4DRef<double> q_v(w.data()+lq-1, jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wk_v(wk.data()+lqwk-1, jdim, kdim, idim, 5);
                            int is1=1, ie1=idim, js1=1, je1=jdim, ks1=1, ke1=kdim;
                            rotateq_ns::rotateq(jdim, kdim, idim, q_v, wk_v, is1, ie1,
                                js1, je1, ks1, ke1, dthtx, dthty, dthtz);
                        }
                        for (l = 1; l <= jdim*kdim*idim*5; l++) {
                            w(l+lq-1) = wk(l+lqwk-1);
                        }
                    }
                }
            }

            // accumulate fluxes at embedded grid boundaries
            ifluxa = 0;
            if (level >= lglobal && level != levt_c) ifluxa = std::min(1, iconsf);
            if (ifluxa > 0 && iadvance(nbl) >= 0) {
                ifamax = maxbl*7*3;
                FortranArray1DRef<int> nfajki_ref(cmn_nfablk.nfajki, 3, 1);
                fa_ns::fa(nbl, lw, w, mgwk, wk, nwork, iw, ifamax, nwfa, nifa, nfajki_ref,
                    maxbl, maxseg, jdimg, kdimg, idimg, jsg, ksg, isg, jeg,
                    keg, ieg, jbcinfo, kbcinfo, ibcinfo, nblock, nblcg,
                    nou, bou, nbuf, ibufdim);
            } else {
                nwfa         = 0;
                cmn_nfablk.nfajki[0] = 0;
                cmn_nfablk.nfajki[1] = 0;
                cmn_nfablk.nfajki[2] = 0;
            }

            lres  = nwfa + 1;
            lwj0  = lres  + jdim*kdim*(idim-1)*5;
            lwk0  = lwj0  + kdim*idim*22;
            lwi0  = lwk0  + jdim*idim*22;
            lvmuk = lwi0  + kdim*jdim*22;
            lvmuj = lvmuk + 2*(jdim-1)*(idim-1);
            lvmui = lvmuj + 2*(kdim-1)*(idim-1);
            ltot  = lvmui + 2*(kdim-1)*(jdim-1);
            lsafe = nwork - ltot;
            if (lsafe < 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " insufficient memory in mgblk *stopping*  lsafe=%7d", lsafe);
                 termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
            }

            // compute residual = r(q)
            isf = 0;
            if (level > lglobal && ntime == nit1) isf = std::min(1, iconsf);
            {
                // Create typed views from work array w
                FortranArray4DRef<double> q_v(w.data()+lq-1, jdim, kdim, idim, 5);
                FortranArray3DRef<double> qj0_v(w.data()+lqj0-1, kdim, idim, 5);
                FortranArray3DRef<double> qk0_v(w.data()+lqk0-1, jdim, idim, 5);
                FortranArray3DRef<double> qi0_v(w.data()+lqi0-1, jdim, kdim, 5);
                FortranArray2DRef<double> sj_v(w.data()+lsj-1, jdim, kdim);
                FortranArray2DRef<double> sk_v(w.data()+lsk-1, jdim, kdim);
                FortranArray2DRef<double> si_v(w.data()+lsi-1, jdim, kdim);
                FortranArray1DRef<double> vol_v(w.data()+lvol-1, jdim*kdim);
                FortranArray1DRef<double> dtj_v(w.data()+ldtj-1, jdim*kdim);
                FortranArray1DRef<double> x_v(w.data()+lx-1, jdim*kdim*idim);
                FortranArray1DRef<double> y_v(w.data()+ly-1, jdim*kdim*idim);
                FortranArray1DRef<double> z_v(w.data()+lz-1, jdim*kdim*idim);
                FortranArray3DRef<double> vist3d_v(w.data()+lvis-1, jdim, kdim, idim);
                FortranArray3DRef<double> snk0_v(w.data()+lsnk0-1, jdim, idim, 3);
                FortranArray3DRef<double> sni0_v(w.data()+lsni0-1, jdim, kdim, 3);
                FortranArray2DRef<double> res_v(wk.data()+lres-1, jdim*kdim*(idim-1), 5);
                FortranArray1DRef<double> wj0_v(wk.data()+lwj0-1, kdim*idim*22);
                FortranArray1DRef<double> wk0_v(wk.data()+lwk0-1, jdim*idim*22);
                FortranArray1DRef<double> wi0_v(wk.data()+lwi0-1, kdim*jdim*22);
                FortranArray3DRef<double> vmuk_v(wk.data()+lvmuk-1, jdim-1, idim-1, 2);
                FortranArray3DRef<double> vmuj_v(wk.data()+lvmuj-1, kdim-1, idim-1, 2);
                FortranArray3DRef<double> vmui_v(wk.data()+lvmui-1, jdim-1, kdim-1, 2);
                FortranArray1DRef<double> wk_v(wk.data()+ltot-1, lsafe);
                FortranArray3DRef<double> blank_v(w.data()+lblk-1, jdim, kdim, idim);
                FortranArray4DRef<double> xib_v(w.data()+lxib-1, jdim, kdim, idim, nummem);
                FortranArray2DRef<double> sig_v(w.data()+lsig-1, jdim*kdim, idim-1);
                FortranArray2DRef<double> sqtq_v(w.data()+lsqtq-1, jdim*kdim, idim-1);
                FortranArray2DRef<double> g_v(w.data()+lg-1, jdim*kdim, idim-1);
                // twoeqn: tj0(kdim,idim-1,nummem,4), tk0(jdim,idim-1,nummem,4),
                //         ti0(jdim,kdim,nummem,4). The 2nd/3rd/4th dims were wrong
                //         (idim/5/nummem instead of idim-1/nummem/4), so the k-turb
                //         ghost row tk0(j,i,2,1) read the wrong slot.
                FortranArray4DRef<double> tj0_v(w.data()+ltj0-1, kdim, idim-1, nummem, 4);
                FortranArray4DRef<double> tk0_v(w.data()+ltk0-1, jdim, idim-1, nummem, 4);
                FortranArray4DRef<double> ti0_v(w.data()+lti0-1, jdim, kdim, nummem, 4);
                FortranArray3DRef<double> xkb_v(w.data()+lxkb-1, jdim, kdim, idim);
                FortranArray3DRef<double> blnum_v(w.data()+lnbl-1, jdim, kdim, idim);
                // Fortran resid.f: vj0(kdim,idim-1,1,4), vk0(jdim,idim-1,1,4),
                // vi0(jdim,kdim,1,4).  The old (.,.,5,1) shapes made every
                // vi0(j,k,1,m>1) write land far past vi0's slot, clobbering the
                // adjacent bcj/bck/bci arrays (spurious periodic-face bc flags ->
                // garbage I-direction gradients -> NaN residual).
                FortranArray4DRef<double> vj0_v(w.data()+lvj0-1, kdim, idim-1, 1, 4);
                FortranArray4DRef<double> vk0_v(w.data()+lvk0-1, jdim, idim-1, 1, 4);
                FortranArray4DRef<double> vi0_v(w.data()+lvi0-1, jdim, kdim, 1, 4);
                FortranArray3DRef<double> bcj_v(w.data()+lbcj-1, kdim, idim-1, 2);
                FortranArray3DRef<double> bck_v(w.data()+lbck-1, jdim, idim-1, 2);
                FortranArray3DRef<double> bci_v(w.data()+lbci-1, jdim, kdim, 2);
                FortranArray4DRef<double> ux_v(w.data()+lux-1, jdim-1, kdim-1, idim-1, 9);
                FortranArray4DRef<double> xib2_v(w.data()+lxib2-1, jdim, kdim, idim, nummem);
                FortranArray3DRef<double> cmuv_v(w.data()+lcmuv-1, jdim, kdim, idim);
                FortranArray3DRef<double> volj0_v(w.data()+lvolj0-1, kdim, idim, 1);
                FortranArray3DRef<double> volk0_v(w.data()+lvolk0-1, jdim, idim, 1);
                FortranArray3DRef<double> voli0_v(w.data()+lvoli0-1, jdim, kdim, 1);
                FortranArray4DRef<double> qavg_v(w.data()+lqavg-1, jdim, kdim, idim, 5);
                FortranArray1DRef<double> sumn_ref2 = sumn_arr.ref();
                FortranArray1DRef<int> negn_ref2(reinterpret_cast<int*>(negn_arr.data()), nummem);
                double delt_val = (double)deltat[level-1];
                int iover_val = iovrlp(nbl);
                resid_ns::resid(nbl, ntime, jdim, kdim, idim, q_v, qj0_v, qk0_v,
                    qi0_v, sj_v, sk_v, si_v, vol_v, dtj_v, x_v,
                    y_v, z_v, vist3d_v, sni0_v, snk0_v, sni0_v,
                    res_v, wj0_v, wk0_v, wi0_v, vmuk_v,
                    vmuj_v, vmui_v, wk_v, lsafe,
                    isf, iw, wk, delt_val, blank_v, iover_val,
                    nblendg, nblstat, nblstag, xib_v, sig_v,
                    sqtq_v, g_v, tj0_v, tk0_v, ti0_v,
                    xkb_v, blnum_v, vj0_v, vk0_v, vi0_v, bcj_v,
                    bck_v, bci_v, nt, sumn_ref2, negn_ref2, ux_v,
                    xib2_v, cmuv_v, volj0_v, volk0_v, voli0_v,
                    nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl, maxseg,
                    nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim, ibcinfo,
                    jbcinfo, kbcinfo, vormax, ivmax, jvmax, kvmax, idefrm,
                    iadvance, qavg_v, nummem);
            }
            // add 2nd order time terms and subiteration terms (time-accurate multigrid)
            if ((float)dt > 0.0f && iadvance(nbl) >= 0) {
                FortranArray4DRef<double> q_v(w.data()+lq-1, jdim, kdim, idim, 5);
                FortranArray4DRef<double> qc0_v(w.data()+lqc0-1, jdim, kdim, idim-1, 5);
                FortranArray4DRef<double> dqc0_v(w.data()+ldqc0-1, jdim, kdim, idim-1, 5);
                FortranArray4DRef<double> res_v(wk.data()+lres-1, jdim, kdim, idim-1, 5);
                FortranArray3DRef<double> vol_v(w.data()+lvol-1, jdim, kdim, idim-1);
                FortranArray3DRef<double> blank_v(w.data()+lblk-1, jdim, kdim, idim);
                int iover_val = iovrlp(nbl);
                resadd_ns::resadd(jdim, kdim, idim, q_v, qc0_v, dqc0_v,
                    res_v, vol_v, iover_val, blank_v);
            }
        }

        if (isklton > 0) {
            writ_buf_ns::writ_buf(nbl, writ_buf_11, nou, bou, nbuf, ibufdim, myhost, myid,
                mycomm, mblk2nd, maxbl);
            if (ivmx == 2 || ivmx == 3) {
                writ_buf_ns::writ_buf(nbl, writ_buf_14, nou, bou, nbuf, ibufdim, myhost, myid,
                    mycomm, mblk2nd, maxbl);
            }
        } else {
            // insure ramped time step/cfl number info gets written
            if (iflagts > 0 && std::abs((float)fmax) > 0.0f && (float)dt < 0.0f) {
                if (((icyc == iflagts+1) ||
                    (icyc == ncyc && iflagts+1 > ncyc)) &&
                    level == levt_c && nbl == nblstat_h && nt == 1) {
                    writ_buf_ns::writ_buf(nbl, writ_buf_11, nou, bou, nbuf, ibufdim, myhost, myid,
                        mycomm, mblk2nd, maxbl);
                }
            }
        }


        // add correction from finer levels (tau)
        if (mgflag > 0) {
            if (level != levt_c) {
                if (mblk2nd(nbl) == myid) {
                    if (iadvance(nbl) >= 0) {
                        // tau declares q(jkm_i,5), res(jkm_i1,5), q1(jkm_i,5), qr(jkm_i1,5)
                        // where jkm_i=jdim*kdim*idim, jkm_i1=jdim*kdim*(idim-1)
                        int q_n = jdim*kdim*idim;
                        int r_n = jdim*kdim*(idim-1);
                        FortranArray2DRef<double> q_v2(w.data()+lq-1, q_n, 5);
                        FortranArray2DRef<double> res_v2(wk.data()+lres-1, r_n, 5);
                        FortranArray2DRef<double> q1_v2(w.data()+lq1-1, q_n, 5);
                        FortranArray2DRef<double> qr_v2(w.data()+lqr-1, r_n, 5);
                        tau_ns::tau(mgflag, nbl, jdim, kdim, idim, q_v2, res_v2,
                            q1_v2, qr_v2, lw, w, nou, bou, nbuf, ibufdim, maxbl,
                            maxgr, nblock, igridg, nblcg, jsg, ksg, isg,
                            jeg, keg, ieg, iemg);
                    }
                }
                if (nbl == nblend && kode == 2) kode = 1;
            }
        }

        if (isklton > 0) {
            writ_buf_ns::writ_buf(nbl, writ_buf_11, nou, bou, nbuf, ibufdim, myhost, myid,
                mycomm, mblk2nd, maxbl);
        }
        // compute force coefficients and L2 norm of residual
        if (level >= lglobal && iaccnt > 0 && ntime == nit) {
            if (myid == mblk2nd(nbl) || myid == myhost) {
                // Create typed views for resp
                int q_sz5 = jdim*kdim*idim*5;
                (void)q_sz5;
                FortranArray2DRef<double> q_v2(w.data()+lq-1, jdim*kdim*idim, 5);
                FortranArray2DRef<double> sj_v2(w.data()+lsj-1, jdim*kdim, 1);
                FortranArray2DRef<double> sk_v2(w.data()+lsk-1, jdim*kdim, 1);
                FortranArray2DRef<double> si_v2(w.data()+lsi-1, jdim*kdim, 1);
                FortranArray1DRef<double> vol_v(w.data()+lvol-1, jdim*kdim);
                FortranArray1DRef<double> x_v(w.data()+lx-1, jdim*kdim*idim);
                FortranArray1DRef<double> y_v(w.data()+ly-1, jdim*kdim*idim);
                FortranArray1DRef<double> z_v(w.data()+lz-1, jdim*kdim*idim);
                int res_sz = jdim*kdim*(idim-1)*5;
                FortranArray2DRef<double> res_v2(wk.data()+lres-1, res_sz, 1);
                FortranArray1DRef<double> wk0_v(wk.data()+lwk0-1, jdim*idim*22);
                FortranArray3DRef<double> vmuk_v(wk.data()+lvmuk-1, jdim-1, idim-1, 2);
                FortranArray3DRef<double> vmuj_v(wk.data()+lvmuj-1, kdim-1, idim-1, 2);
                FortranArray3DRef<double> vmui_v(wk.data()+lvmui-1, jdim-1, kdim-1, 2);
                FortranArray1DRef<double> wk_v(wk.data()+ltot-1, lsafe);
                FortranArray3DRef<double> bcj_v(w.data()+lbcj-1, kdim, idim-1, 2);
                FortranArray3DRef<double> bck_v(w.data()+lbck-1, jdim, idim-1, 2);
                FortranArray3DRef<double> bci_v(w.data()+lbci-1, jdim, kdim, 2);
                FortranArray3DRef<double> blank_v(w.data()+lblk-1, jdim, kdim, idim);
                FortranArray4DRef<double> xtbj_v(w.data()+lxtbj-1, kdim, idim-1, 3, 2);
                FortranArray4DRef<double> xtbk_v(w.data()+lxtbk-1, jdim, idim-1, 3, 2);
                FortranArray4DRef<double> xtbi_v(w.data()+lxtbi-1, jdim, kdim, 3, 2);
                FortranArray4DRef<double> qc0_v(w.data()+lqc0-1, jdim, kdim, idim-1, 5);
                FortranArray4DRef<double> dqc0_v(w.data()+ldqc0-1, jdim, kdim, idim-1, 5);
                FortranArray4DRef<double> qj0_v(w.data()+lqj0-1, kdim, idim-1, 5, 4);
                FortranArray4DRef<double> qk0_v(w.data()+lqk0-1, jdim, idim-1, 5, 4);
                FortranArray4DRef<double> qi0_v(w.data()+lqi0-1, jdim, kdim, 5, 4);
                FortranArray1DRef<int> nnegb_ref2(reinterpret_cast<int*>(nnegb_arr.data()), nummem);
                FortranArray1DRef<double> rmstb_ref2 = rmstb_arr.ref();
                FortranArray1DRef<double> sumn_ref2  = sumn_arr.ref();
                FortranArray1DRef<int> negn_ref2(reinterpret_cast<int*>(negn_arr.data()), nummem);
                resp_ns::resp(nbl, ntime, ntt0, jdim, kdim, idim, q_v2, sj_v2,
                    sk_v2, si_v2, vol_v, x_v, y_v, z_v,
                    res_v2, wk0_v, vmuk_v, vmuj_v,
                    vmui_v, wk_v, lsafe, nblstag, nblendg,
                    nblstat, nblendt, nptsr[ntime-1], nptsrb,
                    nnegb_ref2, bcj_v, bck_v, bci_v,
                    blank_v, nt, sumn_ref2, negn_ref2, rmst,
                    xtbj_v, xtbk_v, xtbi_v, qc0_v, dqc0_v,
                    iseq, qj0_v, qk0_v, qi0_v, maxbl, maxgr,
                    maxseg, rms, rmsb, rmstb_ref2, clw, cdw, cdpw, cdvw,
                    cxw, cyw, czw, cmxw, cmyw, cmzw,
                    n_clcd, clcd, nblocks_clcd, blocks_clcd,
                    chdw, swetw, fmdotw,
                    cfttotw, cftmomw, cftpw, cftvw, rmstr,
                    nneg, ihstry, ngrid, nblg, iemg, levelg,
                    iviscg, itrans, irotat, iforce, swett, clt, cdt,
                    cxt, cyt, czt, cmxt, cmyt, cmzt, cdpt, cdvt, nbci0,
                    nbcj0, nbck0, nbcidim, nbcjdim, nbckdim, ibcinfo,
                    jbcinfo, kbcinfo, resmx, imx, jmx, kmx, vormax,
                    ivmax, jvmax, kvmax, sx, sy, sz, stot, pav, ptav,
                    tav, ttav, xmav, fmdot, cfxp, cfyp, cfzp, cfdp,
                    cflp, cftp, cfxv, cfyv, cfzv, cfdv, cflv, cftv,
                    cfxmom, cfymom, cfzmom, cfdmom, cflmom, cftmom,
                    cfxtot, cfytot, cfztot, cfdtot, cfltot, cfttot,
                    ncs, icsinfo, myid, myhost, mycomm, mblk2nd, nou,
                    bou, nbuf, ibufdim, ncycmax, maxcs, thetay, iadvance,
                    idefrm, igridg, nummem);
            }
            if (nbl == nblendg) iaccnt = 0;
        }

        if (iadvance(nbl) < 0) goto label_7999;

        // smooth residual, if desired
        if (mblk2nd(nbl) == myid) {
            if (issr > 0) {
                icall = idim - 1;
                {
                    FortranArray4DRef<double> res_v4(wk.data()+lres-1, jdim, kdim, idim-1, 5);
                    FortranArray2DRef<double> wk_ltot(wk.data()+ltot-1, jdim, idim+kdim);
                    double epsssr_d[3];
                    epsssr_d[0] = (double)cmn_mgv.epsssr[0];
                    epsssr_d[1] = (double)cmn_mgv.epsssr[1];
                    epsssr_d[2] = (double)cmn_mgv.epsssr[2];
                    FortranArray1DRef<double> epsssr_ref2(epsssr_d, 3);
                    rsmooth_ns::rsmooth(epsssr_ref2, idim, jdim, kdim, icall, res_v4, wk_ltot,
                        nou, bou, nbuf, ibufdim);
                }
            }
        }

        // update solution
        if (ntime <= nit) {
            if (mblk2nd(nbl) == myid) {
                // Create typed views for update
                FortranArray4DRef<double> q_v(w.data()+lq-1, jdim, kdim, idim, 5);
                FortranArray2DRef<double> qj0_v2(w.data()+lqj0-1, kdim*idim, 5);
                FortranArray2DRef<double> qk0_v2(w.data()+lqk0-1, jdim*idim, 5);
                FortranArray2DRef<double> qi0_v2(w.data()+lqi0-1, jdim*kdim, 5);
                FortranArray2DRef<double> sj_v2(w.data()+lsj-1, jdim, kdim);
                FortranArray2DRef<double> sk_v2(w.data()+lsk-1, jdim, kdim);
                FortranArray2DRef<double> si_v2(w.data()+lsi-1, jdim, kdim);
                FortranArray1DRef<double> vol_v(w.data()+lvol-1, jdim*kdim);
                FortranArray1DRef<double> dtj_v(w.data()+ldtj-1, jdim*kdim);
                FortranArray3DRef<double> vist3d_v(w.data()+lvis-1, jdim, kdim, idim);
                FortranArray3DRef<double> blank_v(w.data()+lblk-1, jdim, kdim, idim);
                FortranArray1DRef<double> x_v(w.data()+lx-1, jdim*kdim*idim);
                FortranArray1DRef<double> y_v(w.data()+ly-1, jdim*kdim*idim);
                FortranArray1DRef<double> z_v(w.data()+lz-1, jdim*kdim*idim);
                int res_sz = jdim*kdim*(idim-1)*5;
                FortranArray4DRef<double> res_v(wk.data()+lres-1, jdim*kdim, 1, idim-1, 5);
                FortranArray1DRef<double> wk0_v(wk.data()+lwk0-1, jdim*idim*22);
                FortranArray3DRef<double> vmuk_v(wk.data()+lvmuk-1, jdim-1, idim-1, 2);
                FortranArray3DRef<double> vmuj_v(wk.data()+lvmuj-1, kdim-1, idim-1, 2);
                FortranArray3DRef<double> vmui_v(wk.data()+lvmui-1, jdim-1, kdim-1, 2);
                FortranArray1DRef<double> wk_v(wk.data()+ltot-1, lsafe);
                int iover_val = iovrlp(nbl);
                FortranArray4DRef<double> vk0_v(w.data()+lvk0-1, jdim, idim-1, 1, 4);
                FortranArray3DRef<double> bcj_v(w.data()+lbcj-1, kdim, idim-1, 2);
                FortranArray3DRef<double> bck_v(w.data()+lbck-1, jdim, idim-1, 2);
                FortranArray3DRef<double> bci_v(w.data()+lbci-1, jdim, kdim, 2);
                FortranArray3DRef<double> volk0_v(w.data()+lvolk0-1, jdim, idim, 1);
                FortranArray4DRef<double> xib_v(w.data()+lxib-1, jdim, kdim, idim, nummem);
                // update: tk0(jdim,idim-1,nummem,4) — dims were wrong (idim/5/nummem).
                FortranArray4DRef<double> tk0_v(w.data()+ltk0-1, jdim, idim-1, nummem, 4);
                FortranArray3DRef<double> cmuv_v(w.data()+lcmuv-1, jdim, kdim, idim);
                FortranArray4DRef<double> ux_v(w.data()+lux-1, jdim-1, kdim-1, idim-1, 9);
                update_ns::update(jdim, kdim, idim, q_v, qj0_v2, qk0_v2, qi0_v2,
                    sj_v2, sk_v2, si_v2, vol_v, dtj_v, vist3d_v, blank_v,
                    x_v, y_v, z_v, res_v, wk0_v,
                    vmuk_v, vmuj_v, vmui_v, wk_v, lsafe,
                    nbl, iover_val, vk0_v, bcj_v, bck_v, bci_v,
                    nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl,
                    volk0_v, xib_v, tk0_v, cmuv_v,
                    iadvance, nummem, ux_v);

                // update the overlapped values if chimera scheme is used
                if (iovrlp(nbl) == 1) {
                    ldim = 5;
                    if (isklton == 1) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " interpolating for cells overlapped to block%4d", nbl);
                    }
                    need = (idim+1)*(jdim+1)*(kdim+1)*5;
                    if (need > nwork) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "stopping...not enough work space for subroutine intrbc (q)");
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "need, have %11d%11d", need, nwork);
                         termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    // intrbc: q is 4D array (jdim,kdim,idim,ldim)
                    FortranArray4DRef<double> q_intr(w.data()+lq-1, jdim, kdim, idim, ldim);
                    FortranArray4DRef<double> qj0_intr(w.data()+lqj0-1, kdim, idim, ldim, 1);
                    FortranArray4DRef<double> qk0_intr(w.data()+lqk0-1, jdim, idim, ldim, 1);
                    FortranArray4DRef<double> qi0_intr(w.data()+lqi0-1, jdim, kdim, ldim, 1);
                    FortranArray4DRef<double> wk_intr(wk.data(), jdim, kdim, idim, ldim);
                    FortranArray3DRef<double> bcj_v(w.data()+lbcj-1, kdim, idim-1, 2);
                    FortranArray3DRef<double> bck_v(w.data()+lbck-1, jdim, idim-1, 2);
                    FortranArray3DRef<double> bci_v(w.data()+lbci-1, jdim, kdim, 2);
                    int icorr_v = 1, iset_v = 1;
                    intrbc_ns::intrbc(q_intr, jdim, kdim, idim, nbl, ldim, maxbl,
                        iitot, lig, iipntsg, dxintg, dyintg, dzintg,
                        iiig, jjig, kkig, qb, qj0_intr, qk0_intr,
                        qi0_intr, wk_intr, bcj_v, bck_v, bci_v, nou,
                        bou, nbuf, ibufdim, icorr_v, iset_v);
                    // turbulence quantities
                    if (iviscg(nbl,1) >= 2 || iviscg(nbl,2) >= 2 || iviscg(nbl,3) >= 2) {
                        ldim = 1;
                        FortranArray4DRef<double> vis_intr(w.data()+lvis-1, jdim, kdim, idim, ldim);
                        FortranArray4DRef<double> vj0_intr(w.data()+lvj0-1, kdim, idim, ldim, 1);
                        FortranArray4DRef<double> vk0_intr(w.data()+lvk0-1, jdim, idim, ldim, 1);
                        FortranArray4DRef<double> vi0_intr(w.data()+lvi0-1, jdim, kdim, ldim, 1);
                        FortranArray4DRef<double> wk_intr2(wk.data(), jdim, kdim, idim, ldim);
                        FortranArray3DRef<double> bcj_v2(w.data()+lbcj-1, kdim, idim-1, 2);
                        FortranArray3DRef<double> bck_v2(w.data()+lbck-1, jdim, idim-1, 2);
                        FortranArray3DRef<double> bci_v2(w.data()+lbci-1, jdim, kdim, 2);
                        int icorr_v2 = 1, iset_v2 = 2;
                        intrbc_ns::intrbc(vis_intr, jdim, kdim, idim, nbl, ldim, maxbl,
                            iitot, lig, iipntsg, dxintg, dyintg, dzintg,
                            iiig, jjig, kkig, qb, vj0_intr, vk0_intr,
                            vi0_intr, wk_intr2, bcj_v2, bck_v2, bci_v2, nou,
                            bou, nbuf, ibufdim, icorr_v2, iset_v2);
                    }
                    if (iviscg(nbl,1) >= 4 || iviscg(nbl,2) >= 4 || iviscg(nbl,3) >= 4) {
                        ldim = nummem;
                        FortranArray4DRef<double> xib_intr(w.data()+lxib-1, jdim, kdim, idim, ldim);
                        FortranArray4DRef<double> tj0_intr(w.data()+ltj0-1, kdim, idim, ldim, 1);
                        FortranArray4DRef<double> tk0_intr(w.data()+ltk0-1, jdim, idim, ldim, 1);
                        FortranArray4DRef<double> ti0_intr(w.data()+lti0-1, jdim, kdim, ldim, 1);
                        FortranArray4DRef<double> wk_intr3(wk.data(), jdim, kdim, idim, ldim);
                        FortranArray3DRef<double> bcj_v3(w.data()+lbcj-1, kdim, idim-1, 2);
                        FortranArray3DRef<double> bck_v3(w.data()+lbck-1, jdim, idim-1, 2);
                        FortranArray3DRef<double> bci_v3(w.data()+lbci-1, jdim, kdim, 2);
                        int icorr_v3 = 1, iset_v3 = 3;
                        intrbc_ns::intrbc(xib_intr, jdim, kdim, idim, nbl, ldim, maxbl,
                            iitot, lig, iipntsg, dxintg, dyintg, dzintg,
                            iiig, jjig, kkig, qb, tj0_intr, tk0_intr,
                            ti0_intr, wk_intr3, bcj_v3, bck_v3, bci_v3, nou,
                            bou, nbuf, ibufdim, icorr_v3, iset_v3);
                    }
                } // end if iovrlp
            } // end if ntime <= nit && myid
        } // end if iadvance

        // restrict q and r to coarser grid levels
        if (mgflag > 0 && ntime == nit1 && level != levb) {
            if (mblk2nd(nbl) == myid) {
                if (level <= lglobal) {
                    // global grids
                    lqc    = lw(1,  nbl+1);
                    lvolc  = lw(8,  nbl+1);
                    lvtc   = lw(13, nbl+1);
                    lqrc   = lw(17, nbl+1);
                    lxibc  = lw(19, nbl+1);
                    lqc0c  = lw(34, nbl+1);
                    ldqc0c = lw(35, nbl+1);
                    nroom  = nwork - lwj0;
                    mdim   = jdim*kdim*idim*5;
                    if (nroom < mdim) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " insufficient memory in mgblk for collq *stopping*");
                         termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    // collq: typed views
                    {
                        int jj2l = jj2, kk2l = kk2, ii2l = ii2;
                        FortranArray4DRef<double> q_v(w.data()+lq-1, jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qc_v(w.data()+lqc-1, jj2l, kk2l, ii2l, 5);
                        FortranArray3DRef<double> vol_v(w.data()+lvol-1, jdim, kdim, idim);
                        FortranArray3DRef<double> volc_v(w.data()+lvolc-1, jj2l, kk2l, ii2l);
                        FortranArray4DRef<double> res_v(wk.data()+lres-1, jdim, kdim, idim-1, 5);
                        FortranArray4DRef<double> qr_v(w.data()+lqrc-1, jj2l, kk2l, ii2l-1, 5); // collq.f: qr(jj2,kk2,ii2-1,5)
                        FortranArray4DRef<double> qw_v(wk.data()+lwj0-1, jdim, kdim, idim, 5);
                        FortranArray3DRef<double> vistf_v(w.data()+lvis-1, jdim, kdim, idim);
                        FortranArray3DRef<double> vistc_v(w.data()+lvtc-1, jj2l, kk2l, ii2l);
                        FortranArray4DRef<double> tursavf_v(w.data()+lxib-1, jdim, kdim, idim, nummem);
                        FortranArray4DRef<double> tursavc_v(w.data()+lxibc-1, jj2l, kk2l, ii2l, nummem);
                        collq_ns::collq(q_v, qc_v, vol_v, volc_v, jdim, kdim, idim,
                            jj2l, kk2l, ii2l, res_v, qr_v, qw_v,
                            vistf_v, vistc_v, tursavf_v, tursavc_v,
                            nbl, nou, bou, nbuf, ibufdim, nummem);
                    }
                    // restrict qc0 and dqc0 for time-accurate multigrid
                    if ((float)dt > 0.0f) {
                        int jj2l = jj2, kk2l = kk2, ii2l = ii2;
                        FortranArray4DRef<double> qc0_v(w.data()+lqc0-1, jdim, kdim, idim-1, 5);
                        FortranArray4DRef<double> qc0c_v(w.data()+lqc0c-1, jj2l, kk2l, ii2l-1, 5);
                        FortranArray3DRef<double> vol_v(w.data()+lvol-1, jdim, kdim, idim-1);
                        FortranArray3DRef<double> volc_v(w.data()+lvolc-1, jj2l, kk2l, ii2l-1);
                        FortranArray4DRef<double> dqc0_v(w.data()+ldqc0-1, jdim, kdim, idim-1, 5);
                        FortranArray4DRef<double> dqc0c_v(w.data()+ldqc0c-1, jj2l, kk2l, ii2l-1, 5);
                        FortranArray4DRef<double> qw_v(wk.data()+lwj0-1, jdim, kdim, idim, 5);
                        collqc0_ns::collqc0(qc0_v, qc0c_v, vol_v, volc_v,
                            jdim, kdim, idim, jj2l, kk2l, ii2l,
                            dqc0_v, dqc0c_v, qw_v, nbl,
                            nou, bou, nbuf, ibufdim);
                    }
                } else {
                    // embedded grids
                    jc    = jdimg(nblc);
                    kc    = kdimg(nblc);
                    ic    = idimg(nblc);
                    lqc   = lw(1,  nblc);
                    lvolc = lw(8,  nblc);
                    lqrc  = lw(17, nbl);
                    lvtc  = lw(13, nblc);
                    lxtc  = lw(19, nblc);
                    nroom = nwork - lwj0;
                    mdim  = jdim*kdim*idim*5;
                    if (nroom < mdim) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " insufficient memory in mgblk for coll2q *stopping*");
                         termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    {
                        FortranArray4DRef<double> q_v(w.data()+lq-1, jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qc_v(w.data()+lqc-1, jc, kc, ic, 5);
                        FortranArray3DRef<double> vol_v(w.data()+lvol-1, jdim, kdim, idim);
                        FortranArray3DRef<double> volc_v(w.data()+lvolc-1, jc, kc, ic);
                        FortranArray4DRef<double> res_v(wk.data()+lres-1, jdim, kdim, idim-1, 5);
                        FortranArray4DRef<double> qr_v(w.data()+lqrc-1, jdim, kdim, idim-1, 5);
                        FortranArray4DRef<double> qw_v(wk.data()+lwj0-1, jdim, kdim, idim, 5);
                        FortranArray3DRef<double> vistf_v(w.data()+lvis-1, jdim, kdim, idim);
                        FortranArray3DRef<double> vistc_v(w.data()+lvtc-1, jc, kc, ic);
                        FortranArray4DRef<double> tursavf_v(w.data()+lxib-1, jdim, kdim, idim, nummem);
                        FortranArray4DRef<double> tursavc_v(w.data()+lxtc-1, jc, kc, ic, nummem);
                        coll2q_ns::coll2q(q_v, qc_v, vol_v, volc_v, jdim, kdim, idim,
                            jc, kc, ic, res_v, qr_v, qw_v,
                            js, ks, is, je, ke, ie, nbl, nblc,
                            vistf_v, vistc_v, tursavf_v, tursavc_v,
                            nou, bou, nbuf, ibufdim, nummem);
                    }
                }
            }
        }


        label_7999:;
        if (isklton > 0) {
            writ_buf_ns::writ_buf(nbl, writ_buf_11, nou, bou, nbuf, ibufdim, myhost, myid,
                mycomm, mblk2nd, maxbl);
        }
    } // end do 8000 nbl

    } // end do 6500 ntime


    // if using fixed Cl option, update angle of attack
    if (icycupdt > 0) {
        if (ntt % icycupdt == 0) {
            if (ialphit != 0 && level == levt_c) {
                newalpha_ns::newalpha(ncycmax, rms, clw, myid, myhost, mycomm);
                for (nbl = 1; nbl <= nblock; nbl++) {
                    if (mblk2nd(nbl) == myid) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        angx = 0.0;
                        angy = -(double)dalpha;
                        angz = 0.0;
                        {
                            FortranArray4DRef<double> q_rot(w.data()+lq-1, jdim, kdim, idim, 5);
                            int ista=1, iend=idim, jsta=1, jend=jdim, ksta=1, kend=kdim;
                            rotateq_ns::rotateq(jdim, kdim, idim, q_rot, q_rot, ista, iend,
                                jsta, jend, ksta, kend, angx, angy, angz);
                        }
                    }
                }
            }
        }
    }

    if (ilc == 0) {
        cmn_info.nsm[level-levb] = cmn_info.nsm[level-levb] - 1;
        if (cmn_info.nsm[level-levb] < 0) cmn_info.nsm[level-levb] = ngam - 1;
        ilc = 1;
    }

    if (isklton > 0 && myid == myhost) {
        fprintf(fortran_get_unit(11), "\n***** END RESIDUAL/RESTRICTION LOOP *****\n\n");
    }

    kode = 1;

    if (kxpand == -1) {
        // finished with downward leg; begin upward leg
        goto label_7500;
    }

    if (level == levb && levb == levt_c) {
        // single level case; start new cycle
        if (isklton == 1 && myid == myhost) {
            fprintf(fortran_get_unit(11), "\n***** END MULTIGRID CYCLE *****\n\n");
        }
        goto label_6999;
    }

    kode  = 2;
    mode  = 1;
    level = level - 1;
    ntime = 0;
    nit   = mit(level-levb+1, iseq);

    if (level > levb) {
        // not yet at bottom of cycle; continue downward leg
        goto label_9000;
    } else {
        // at bottom of cycle; complete downward leg
        nit    = mit(1, iseq) + mtt;
        kxpand = -1;
        goto label_9000;
    }

    // label 7500: begin/continue upward leg of multigrid cycle
    label_7500:;

    if (mgflag != 0) {
        if (myid == myhost) {
            if (level < lglobal) {
                if (isklton > 0) fprintf(fortran_get_unit(11), "\n***** BEGINNING PROLONGATION/CORRECTION LOOP, level, kode, mode =%3d%3d%3d\n", level, kode, mode);
            } else if (mgflag > 1) {
                if (isklton > 0) fprintf(fortran_get_unit(11), "\n***** BEGINNING PROLONGATION/CORRECTION LOOP, level, kode, mode =%3d%3d%3d\n", level, kode, mode);
            }
        }

        // do 7020 nbl=1,nblock
        for (nbl = 1; nbl <= nblock; nbl++) {
            if (iadvance(nbl) < 0) continue;
            if (level != levelg(nbl)) continue;
            lqc  = lw(1,  nbl);
            lq1c = lw(16, nbl);
            if (level < lglobal) {
                // add corrections to global grids
                nblf = nbl - 1;
                if (mblk2nd(nbl) == myid) {
                    lead_ns::lead(nblf, lw, lw2, maxbl);
                    iwk1  = 1;
                    iwk2  = iwk1 + jdim*kdim*idim*5;
                    iwk3  = iwk2 + jj2*kk2*ii2*5;
                    nroom = nwork - iwk3;
                    mdim  = jdim*kk2*ii2*5;
                    if (nroom < mdim) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " insufficient memory in mgblk for addx *stopping*");
                         termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    if (isklton > 0) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " adding corrections from block%4d to block%4d, igrid=%4d", nbl, nblf, igridg(nblf));
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " fine grid dimensions: jdim,kdim,idim =%5d%5d%5d", jdim, kdim, idim);
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " coarse grid dimensions: jj2,kk2,ii2 =%5d%5d%5d", jj2, kk2, ii2);
                    }
                        FortranArray4DRef<double> q_addx(w.data()+lq-1, jdim, kdim, idim, 5);
                        FortranArray4DRef<double> qq_addx(w.data()+lqc-1, jj2, kk2, ii2, 5);
                        FortranArray4DRef<double> q1_addx(w.data()+lq1c-1, jj2, kk2, ii2, 5);
                        FortranArray4DRef<double> dq_addx(wk.data()+iwk1-1, jdim, kdim, idim, 5);
                        FortranArray4DRef<double> wq_addx(wk.data()+iwk2-1, jj2, kk2, ii2, 5);
                        FortranArray4DRef<double> wqj_addx(wk.data()+iwk3-1, jdim, kk2, ii2, 5);
                        FortranArray3DRef<double> blank_addx(w.data()+lblk-1, jdim, kdim, idim);
                        int ll5 = 5;
                        mgblk_ns::addx(q_addx, qq_addx, jdim, kdim, idim, jj2, kk2, ii2,
                            q1_addx, dq_addx, wq_addx, wqj_addx, nbl, blank_addx,
                            nou, bou, nbuf, ibufdim, ll5, myid);
                    // add corrections to chimera grids
                    if (iovrlp(nblf) == 1) {
                        if (isklton == 1) {
                            nou(1) = std::min(nou(1)+1, ibufdim);
                            std::snprintf(bou(nou(1),1), 120, " interpolating for cells overlapped to block%4d", nblf);
                        }
                        ldim = 5;
                        {
                            FortranArray4DRef<double> q_v(w.data()+lq-1, jdim, kdim, idim, ldim);
                            FortranArray4DRef<double> qj0_v(w.data()+lqj0-1, kdim, idim-1, ldim, 4);
                            FortranArray4DRef<double> qk0_v(w.data()+lqk0-1, jdim, idim-1, ldim, 4);
                            FortranArray4DRef<double> qi0_v(w.data()+lqi0-1, jdim, kdim, ldim, 4);
                            FortranArray4DRef<double> qq_v(wk.data(), jdim+1, kdim+1, idim+1, ldim);
                            FortranArray3DRef<double> bcj_v(w.data()+lbcj-1, kdim, idim-1, 2);
                            FortranArray3DRef<double> bck_v(w.data()+lbck-1, jdim, idim-1, 2);
                            FortranArray3DRef<double> bci_v(w.data()+lbci-1, jdim, kdim, 2);
                            int icorr1 = 1, iset1 = 1;
                            mgblk_ns::intrbc(q_v, jdim, kdim, idim, nblf, ldim, maxbl,
                                iitot, lig, iipntsg, dxintg, dyintg, dzintg,
                                iiig, jjig, kkig, qb, qj0_v, qk0_v,
                                qi0_v, qq_v, bcj_v, bck_v, bci_v, nou,
                                bou, nbuf, ibufdim, icorr1, iset1);
                        }
                        // turbulence quantities
                        if (iviscg(nblf,1) >= 2 || iviscg(nblf,2) >= 2 || iviscg(nblf,3) >= 2) {
                            ldim = 1;
                            FortranArray4DRef<double> q_v(w.data()+lvis-1, jdim, kdim, idim, ldim);
                            FortranArray4DRef<double> qj0_v(w.data()+lvj0-1, kdim, idim-1, ldim, 4);
                            FortranArray4DRef<double> qk0_v(w.data()+lvk0-1, jdim, idim-1, ldim, 4);
                            FortranArray4DRef<double> qi0_v(w.data()+lvi0-1, jdim, kdim, ldim, 4);
                            FortranArray4DRef<double> qq_v(wk.data(), jdim+1, kdim+1, idim+1, ldim);
                            FortranArray3DRef<double> bcj_v(w.data()+lbcj-1, kdim, idim-1, 2);
                            FortranArray3DRef<double> bck_v(w.data()+lbck-1, jdim, idim-1, 2);
                            FortranArray3DRef<double> bci_v(w.data()+lbci-1, jdim, kdim, 2);
                            int icorr1 = 1, iset2 = 2;
                            mgblk_ns::intrbc(q_v, jdim, kdim, idim, nblf, ldim, maxbl,
                                iitot, lig, iipntsg, dxintg, dyintg, dzintg,
                                iiig, jjig, kkig, qb, qj0_v, qk0_v,
                                qi0_v, qq_v, bcj_v, bck_v, bci_v, nou,
                                bou, nbuf, ibufdim, icorr1, iset2);
                        }
                        if (iviscg(nblf,1) >= 4 || iviscg(nblf,2) >= 4 || iviscg(nblf,3) >= 4) {
                            ldim = nummem;
                            FortranArray4DRef<double> q_v(w.data()+lxib-1, jdim, kdim, idim, ldim);
                            FortranArray4DRef<double> qj0_v(w.data()+ltj0-1, kdim, idim-1, ldim, 4);
                            FortranArray4DRef<double> qk0_v(w.data()+ltk0-1, jdim, idim-1, ldim, 4);
                            FortranArray4DRef<double> qi0_v(w.data()+lti0-1, jdim, kdim, ldim, 4);
                            FortranArray4DRef<double> qq_v(wk.data(), jdim+1, kdim+1, idim+1, ldim);
                            FortranArray3DRef<double> bcj_v(w.data()+lbcj-1, kdim, idim-1, 2);
                            FortranArray3DRef<double> bck_v(w.data()+lbck-1, jdim, idim-1, 2);
                            FortranArray3DRef<double> bci_v(w.data()+lbci-1, jdim, kdim, 2);
                            int icorr1 = 1, iset3 = 3;
                            mgblk_ns::intrbc(q_v, jdim, kdim, idim, nblf, ldim, maxbl,
                                iitot, lig, iipntsg, dxintg, dyintg, dzintg,
                                iiig, jjig, kkig, qb, qj0_v, qk0_v,
                                qi0_v, qq_v, bcj_v, bck_v, bci_v, nou,
                                bou, nbuf, ibufdim, icorr1, iset3);
                        }
                    }
                }
            } else {
                // add corrections to embedded grids
                if (mgflag > 1) {
                    ipass = 1;
                    jc    = jdimg(nbl);
                    kc    = kdimg(nbl);
                    ic    = idimg(nbl);
                    for (nblz = 1; nblz <= nblock; nblz++) {
                        if (nblz == nbl) continue;
                        nblcc = nblcg(nblz);
                        if (nblcc == nbl) {
                            igridc = igridg(nblz);
                            if (iemg(igridc) > 0) {
                                lead_ns::lead(nblz, lw, lw2, maxbl);
                                iwk1  = 1;
                                iwk2  = iwk1 + jdim*kdim*idim*5;
                                iwk3  = iwk2 + jc*kc*ic*5;
                                iwk4  = iwk3 + jdim*kc*ic;
                                nroom = nwork - iwk4;
                                mdim  = jdim*kdim*ic;
                                if (nroom < mdim) {
                                    nou(1) = std::min(nou(1)+1, ibufdim);
                                    std::snprintf(bou(nou(1),1), 120, " insufficient memory in mgblk for add2x *stopping*");
                                     termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                                }
                                if (isklton > 0) {
                                    nou(1) = std::min(nou(1)+1, ibufdim);
                                    std::snprintf(bou(nou(1),1), 120, " adding corrections from block%4d to block%4d, igrid=%4d", nbl, nblz, igridg(nblz));
                                    nou(1) = std::min(nou(1)+1, ibufdim);
                                    std::snprintf(bou(nou(1),1), 120, " fine grid dimensions: jdim,kdim,idim =%5d%5d%5d", jdim, kdim, idim);
                                    nou(1) = std::min(nou(1)+1, ibufdim);
                                    std::snprintf(bou(nou(1),1), 120, " coarse grid dimensions: jj2,kk2,ii2 =%5d%5d%5d", jj2, kk2, ii2);
                                }
                                {
                                    FortranArray4DRef<double> q_add2x(w.data()+lq-1, jdim, kdim, idim, 5);
                                    FortranArray4DRef<double> qc_add2x(w.data()+lqc-1, jc, kc, ic, 5);
                                    FortranArray4DRef<double> q1_add2x(w.data()+lq1c-1, jc, kc, ic, 5);
                                    FortranArray4DRef<double> dq_add2x(wk.data()+iwk1-1, jdim, kdim, idim, 5);
                                    FortranArray4DRef<double> wq_add2x(wk.data()+iwk2-1, jc, kc, ic, 5);
                                    FortranArray3DRef<double> wqj_add2x(wk.data()+iwk3-1, jdim, kc, ic);
                                    FortranArray3DRef<double> wqjk_add2x(wk.data()+iwk4-1, jdim, kdim, ic);
                                    int ll5 = 5;
                                    mgblk_ns::add2x(q_add2x, qc_add2x, jdim, kdim, idim,
                                        jc, kc, ic, q1_add2x, dq_add2x, wq_add2x,
                                        wqj_add2x, wqjk_add2x, js, ks, is,
                                        je, ke, ie, ipass, nbl, nblz,
                                        nou, bou, nbuf, ibufdim, ll5, myid);
                                }
                            }
                        }
                        if (isklton > 0) {
                            writ_buf_ns::writ_buf(nblz, writ_buf_11, nou, bou, nbuf, ibufdim, myhost,
                                myid, mycomm, mblk2nd, maxbl);
                        }
                    } // end do 6030 nblz
                }
            }

            if (isklton > 0) {
                writ_buf_ns::writ_buf(nbl, writ_buf_11, nou, bou, nbuf, ibufdim, myhost, myid,
                    mycomm, mblk2nd, maxbl);
            }
        } // end do 7020 nbl

        if (myid == myhost) {
            if (level < lglobal) {
                if (isklton > 0) fprintf(fortran_get_unit(11), "\n***** END PROLONGATION/CORRECTION LOOP *****\n\n");
            } else if (mgflag > 1) {
                if (isklton > 0) fprintf(fortran_get_unit(11), "\n***** END PROLONGATION/CORRECTION LOOP *****\n\n");
            }
        }

        level = level + 1;

        if (level == levt_c) {
            // finished with multigrid cycle
            if (isklton == 1 && myid == myhost) {
                fprintf(fortran_get_unit(11), "\n***** END MULTIGRID CYCLE *****\n\n");
            }
            goto label_6999;
        }

        nit = mtt;
        if (cmn_info.nsm[level-levb] > 0) {
            nit    = mit(level-levb+1, iseq) + mtt;
            kxpand = 1;
        }
        kode  = 1;
        ntime = 0;

        if (nit == 0) {
            // continue with upward leg
            goto label_7500;
        } else {
            // done with upward leg; begin downward leg
            goto label_9000;
        }
    } else {
        // embedded grids without multigrid; set level=levt and start new cycle
        level = levt_c;
    }

    // label 6999: end multigrid cycle
    label_6999:;


    // write restart file if non time-accurate
    if ((float)dt < 0.0f) {
        if (level == levt_c) {
            if (icyc >= ncyc || (nwrest > 0 && icyc/nwrest*nwrest == icyc)) {
                if (myid == myhost) {
                    if (icyc >= ncyc) {
                        iwarneddy = 0;
                    } else {
                        iwarneddy = -1;
                    }
                    if (icgns != 1) {
                        char restrt_c[81]; int rlen = 80;
                        while (rlen > 0 && restrt[rlen-1] == ' ') rlen--;
                        std::memcpy(restrt_c, restrt, rlen); restrt_c[rlen] = '\0';
                        fortran_open_unit(2, restrt_c, "wb");
                        rewind(fortran_get_unit(2));
                        if (ndwrt != 0 && iunst > 1) {
                            if (mseq > 1) {
                                nou(1) = std::min(nou(1)+1, ibufdim);
                                std::snprintf(bou(nou(1),1), 120, "stopping...cannot write dynamic grid file when mseq > 1");
                                 termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                            }
                            fortran_open_unit(98, "dgplot3d.bin", "w");
                            fprintf(fortran_get_unit(98), "%8d\n", ngrid);
                            for (igrid = 1; igrid <= ngrid; igrid++) {
                                nbl = nblg(igrid);
                                iem = iemg(igrid);
                                if (iem == 0) nbl = nbl + (mseq - iseq);
                                lead_ns::lead(nbl, lw, lw2, maxbl);
                                idimt_arr(igrid) = idim;
                                jdimt_arr(igrid) = jdim;
                                kdimt_arr(igrid) = kdim;
                            }
                            for (igrid = 1; igrid <= ngrid; igrid++) {
                                fprintf(fortran_get_unit(98), "%8d%8d%8d", idimt_arr(igrid), jdimt_arr(igrid), kdimt_arr(igrid));
                            }
                            fprintf(fortran_get_unit(98), "\n");
                        }
                        if (iclcd == 1 || iclcd == 2) {
                            char clcds_c[242]; int clen = 241;
                            while (clen > 0 && clcds[clen-1] == ' ') clen--;
                            std::memcpy(clcds_c, clcds, clen); clcds_c[clen] = '\0';
                            fortran_open_unit(102, clcds_c, "wb");
                            rewind(fortran_get_unit(102));
                        }
                    }
                }
                if (iteravg == 1 || iteravg == 2) {
                    if (myid == myhost) {
                        fprintf(fortran_get_unit(11), " writing averaged solution\n");
                        fprintf(fortran_get_unit(11), " averaging counter = %g\n", (double)xnumavg);
                        char avgg_c[81]; int alen = 80;
                        while (alen > 0 && avgg[alen-1] == ' ') alen--;
                        std::memcpy(avgg_c, avgg, alen); avgg_c[alen] = '\0';
                        fprintf(fortran_get_unit(11), " %s\n", avgg_c);
                        char avgq_c[81]; alen = 80;
                        while (alen > 0 && avgq[alen-1] == ' ') alen--;
                        std::memcpy(avgq_c, avgq, alen); avgq_c[alen] = '\0';
                        fprintf(fortran_get_unit(11), " %s\n", avgq_c);
                        rewind(fortran_get_unit(96));
                        rewind(fortran_get_unit(97));
                        fwrite(&ngrid, sizeof(int), 1, fortran_get_unit(96));
                        fwrite(&ngrid, sizeof(int), 1, fortran_get_unit(97));
                        if (ipertavg == 1 || ipertavg == 2) {
                            char avgq2_c[242]; alen = 241;
                            while (alen > 0 && avgq2[alen-1] == ' ') alen--;
                            std::memcpy(avgq2_c, avgq2, alen); avgq2_c[alen] = '\0';
                            fprintf(fortran_get_unit(11), " %s\n", avgq2_c);
                            char avgq2pert_c[242]; alen = 241;
                            while (alen > 0 && avgq2pert[alen-1] == ' ') alen--;
                            std::memcpy(avgq2pert_c, avgq2pert, alen); avgq2pert_c[alen] = '\0';
                            fprintf(fortran_get_unit(11), " %s\n", avgq2pert_c);
                            rewind(fortran_get_unit(95));
                            rewind(fortran_get_unit(98));
                            fwrite(&ngrid, sizeof(int), 1, fortran_get_unit(95));
                            fwrite(&ngrid, sizeof(int), 1, fortran_get_unit(98));
                            for (igrid = 1; igrid <= ngrid; igrid++) {
                                int ib = nblg(igrid) + mseq - iseq;
                                int id = idimg(ib), jd = jdimg(ib), kd = kdimg(ib);
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(95));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(95));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(95));
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(97));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(97));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(97));
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(98));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(98));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(98));
                            }
                        } else {
                            for (igrid = 1; igrid <= ngrid; igrid++) {
                                int ib = nblg(igrid) + mseq - iseq;
                                int id = idimg(ib)-1, jd = jdimg(ib)-1, kd = kdimg(ib)-1;
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(97));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(97));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(97));
                            }
                        }
                    }
                }
                // do 6040 igrid=1,ngrid
                for (igrid = 1; igrid <= ngrid; igrid++) {
                    iskipz = 0;
                    if (igrid == 1) iskipz = 1;
                    nbl = nblg(igrid);
                    iem = iemg(igrid);
                    if (iseq != mseq && iem > 0) continue;
                    if (icgns == 1) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        idima = idim;
                        jdima = jdim;
                        kdima = kdim;
                    }
                    if (iem == 0) nbl = nbl + (mseq - iseq);
                    if (mblk2nd(nbl) == myid || myid == myhost) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        if (icgns == 1) {
                            if (nwork < (idima+1)*(jdima+1)*(kdima+1)) {
                                fprintf(fortran_get_unit(901), " not enough memory for cgns Q write.\n");
                                fprintf(fortran_get_unit(901), " nwork in wk=%6d.  Needed = %6d\n", nwork, (idima+1)*(jdima+1)*(kdima+1));
                                 termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                            }
                        } else {
                            idima = idim;
                            jdima = jdim;
                            kdima = kdim;
                        }
                        {
                            FortranArray4DRef<double> q_wv(w.data()+lq-1, jdim, kdim, idim, 5);
                            FortranArray4DRef<double> qj0_wv(w.data()+lqj0-1, kdim, idim-1, 5, 4);
                            FortranArray4DRef<double> qk0_wv(w.data()+lqk0-1, jdim, idim-1, 5, 4);
                            FortranArray4DRef<double> qi0_wv(w.data()+lqi0-1, jdim, kdim, 5, 4);
                            FortranArray3DRef<double> vist3d_wv(w.data()+lvis-1, jdim, kdim, idim);
                            FortranArray4DRef<double> tursav_wv(w.data()+lxib-1, jdim, kdim, idim, nummem);
                            FortranArray3DRef<double> smin_wv(w.data()+lsnk0-1, jdim-1, kdim-1, idim-1);
                            FortranArray3DRef<double> xjb_wv(w.data()+lsni0-1, jdim-1, kdim-1, idim-1);
                            FortranArray3DRef<double> xkb_wv(w.data()+lxkb-1, jdim-1, kdim-1, idim-1);
                            FortranArray3DRef<double> blnum_wv(w.data()+lnbl-1, jdim-1, kdim-1, idim-1);
                            FortranArray3DRef<double> cmuv_wv(w.data()+lcmuv-1, jdim-1, kdim-1, idim-1);
                            FortranArray3DRef<double> wk_wv(wk.data(), idima+1, jdima+1, kdima+1);
                            FortranArray3DRef<double> bcj_wv(w.data()+lbcj-1, kdim, idim-1, 2);
                            FortranArray3DRef<double> bck_wv(w.data()+lbck-1, jdim, idim-1, 2);
                            FortranArray3DRef<double> bci_wv(w.data()+lbci-1, jdim, kdim, 2);
                            FortranArray4DRef<double> vj0_wv(w.data()+lvj0-1, kdim, idim-1, 1, 4);
                            FortranArray4DRef<double> vk0_wv(w.data()+lvk0-1, jdim, idim-1, 1, 4);
                            FortranArray4DRef<double> vi0_wv(w.data()+lvi0-1, jdim, kdim, 1, 4);
                            FortranArray4DRef<double> tj0_wv(w.data()+ltj0-1, kdim, idim-1, nummem, 4);
                            FortranArray4DRef<double> tk0_wv(w.data()+ltk0-1, jdim, idim-1, nummem, 4);
                            FortranArray4DRef<double> ti0_wv(w.data()+lti0-1, jdim, kdim, nummem, 4);
                            FortranArray3DRef<double> blank_wv(w.data()+lblk-1, jdim, kdim, idim);
                            FortranArray4DRef<double> qavg_wv(w.data()+lqavg-1, jdim, kdim, idim, 5);
                            FortranArray4DRef<double> q2avg_wv(w.data()+lq2avg-1, jdim, kdim, idim, 5);
                            FortranArray3DRef<double> x_wv(w.data()+lx-1, jdim, kdim, idim);
                            FortranArray3DRef<double> y_wv(w.data()+ly-1, jdim, kdim, idim);
                            FortranArray3DRef<double> z_wv(w.data()+lz-1, jdim, kdim, idim);
                            int iover_wv = iovrlp(nbl);
                            wrest_ns::wrest(nbl, jdim, kdim, idim, q_wv, qj0_wv, qk0_wv,
                                qi0_wv, ncycmax, rms, clw, cdw, cdpw, cdvw, cxw, cyw, czw,
                                cmxw, cmyw, cmzw,
                                n_clcd, clcd, nblocks_clcd, blocks_clcd,
                                fmdotw, cftmomw, cftpw, cftvw, cfttotw,
                                rmstr, nneg, iskipz, vist3d_wv, tursav_wv,
                                smin_wv, xjb_wv, xkb_wv, blnum_wv, cmuv_wv,
                                thetay, maxbl, myid, myhost, mycomm, mblk2nd, igrid, wk_wv,
                                idima, jdima, kdima, bcj_wv, bck_wv, bci_wv, vj0_wv,
                                vk0_wv, vi0_wv, tj0_wv, tk0_wv, ti0_wv, blank_wv,
                                iwk, iwork, iover_wv, nou, bou, nbuf, ibufdim,
                                qavg_wv, q2avg_wv, x_wv, y_wv, z_wv, nummem);
                        }
                    }
                } // end do 6040 igrid
                if (myid == myhost) {
                    if (icgns != 1) {
                        fortran_close_unit(2);
                        fortran_close_unit(102);
                    }
                }
            }
        }
    }

    // label 6998: meshdef bypass
    label_6998:;

    if (icyc >= ncyc) goto label_7011;

    label_7011:;
    } // end do 7000 icyc

    // recover original ita
    ita = itatemp;


    // advance time if time-accurate
    if ((float)dt > 0.0f) {
        if (level == levt_c) {
            time     = time + dt;
            time2mc  = time2mc + dt;
            for (nbl = 1; nbl <= nblock; nbl++) {
                time2(nbl) = time2(nbl) + (double)dt;
            }
            timekeep(ntt) = (double)time;
        }
    }

    // aeroelastic corrector step
    if (naesrf > 0) {
        ae_corr_ns::ae_corr(stm, stmi, xs, xxn, gforcn, gforcs, gforcnm,
            gf0, lw, lw2, w, mgwk, maxbl, maxseg,
            aesrfdat, nmds, maxaes, nt, mblk2nd, iseq,
            levelg, iadvance, nblock, icsi, icsf, jcsi, jcsf,
            kcsi, kcsf, myid, nsegdfrm, idfrmseg, iaesurf,
            perturb, aehist, ncycmax, maxsegdg, myhost, mycomm);
    }

    // rigid body mode corrector step
    if (irbtrim != 0) irigb = 0;
    if (irigb == 1) {
        rb_corr_ns::rb_corr(aesrfdat, clw(ntt), cmyw(ntt), maxaes);
    }

    // meshdef = 1 if flow solution is to be bypassed
    if (meshdef == 1) goto label_5990;

    // write restart file if time-accurate
    if ((float)dt > 0.0f) {
        if (level == levt_c) {
            if (nt >= ntstep || (nwrest > 0 && nt/nwrest*nwrest == nt)) {
                if (myid == myhost) {
                    if (nt >= ntstep) {
                        iwarneddy = 0;
                    } else {
                        iwarneddy = -1;
                    }
                    if (icgns != 1) {
                        char restrt_c[81]; int rlen = 80;
                        while (rlen > 0 && restrt[rlen-1] == ' ') rlen--;
                        std::memcpy(restrt_c, restrt, rlen); restrt_c[rlen] = '\0';
                        fortran_open_unit(2, restrt_c, "wb");
                        rewind(fortran_get_unit(2));
                        if (iclcd == 1 || iclcd == 2) {
                            char clcds_c[242]; int clen = 241;
                            while (clen > 0 && clcds[clen-1] == ' ') clen--;
                            std::memcpy(clcds_c, clcds, clen); clcds_c[clen] = '\0';
                            fortran_open_unit(102, clcds_c, "wb");
                            rewind(fortran_get_unit(102));
                        }
                    }
                }
                if (iteravg == 1 || iteravg == 2) {
                    if (myid == myhost) {
                        fprintf(fortran_get_unit(11), " writing averaged solution\n");
                        fprintf(fortran_get_unit(11), " averaging counter = %g\n", (double)xnumavg);
                        fprintf(fortran_get_unit(11), " time step = %d\n", nt);
                        char avgg_c[81]; int alen = 80;
                        while (alen > 0 && avgg[alen-1] == ' ') alen--;
                        std::memcpy(avgg_c, avgg, alen); avgg_c[alen] = '\0';
                        fprintf(fortran_get_unit(11), " %s\n", avgg_c);
                        char avgq_c[81]; alen = 80;
                        while (alen > 0 && avgq[alen-1] == ' ') alen--;
                        std::memcpy(avgq_c, avgq, alen); avgq_c[alen] = '\0';
                        fprintf(fortran_get_unit(11), " %s\n", avgq_c);
                        rewind(fortran_get_unit(96));
                        rewind(fortran_get_unit(97));
                        fwrite(&ngrid, sizeof(int), 1, fortran_get_unit(96));
                        fwrite(&ngrid, sizeof(int), 1, fortran_get_unit(97));
                        if (ipertavg == 1 || ipertavg == 2) {
                            char avgq2_c[242]; alen = 241;
                            while (alen > 0 && avgq2[alen-1] == ' ') alen--;
                            std::memcpy(avgq2_c, avgq2, alen); avgq2_c[alen] = '\0';
                            fprintf(fortran_get_unit(11), " %s\n", avgq2_c);
                            char avgq2pert_c[242]; alen = 241;
                            while (alen > 0 && avgq2pert[alen-1] == ' ') alen--;
                            std::memcpy(avgq2pert_c, avgq2pert, alen); avgq2pert_c[alen] = '\0';
                            fprintf(fortran_get_unit(11), " %s\n", avgq2pert_c);
                            rewind(fortran_get_unit(95));
                            rewind(fortran_get_unit(98));
                            fwrite(&ngrid, sizeof(int), 1, fortran_get_unit(95));
                            fwrite(&ngrid, sizeof(int), 1, fortran_get_unit(98));
                            for (igrid = 1; igrid <= ngrid; igrid++) {
                                int ib = nblg(igrid) + mseq - iseq;
                                int id = idimg(ib), jd = jdimg(ib), kd = kdimg(ib);
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(95));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(95));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(95));
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(97));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(97));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(97));
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(98));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(98));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(98));
                            }
                        } else {
                            for (igrid = 1; igrid <= ngrid; igrid++) {
                                int ib = nblg(igrid) + mseq - iseq;
                                int id = idimg(ib)-1, jd = jdimg(ib)-1, kd = kdimg(ib)-1;
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(96));
                                fwrite(&id, sizeof(int), 1, fortran_get_unit(97));
                                fwrite(&jd, sizeof(int), 1, fortran_get_unit(97));
                                fwrite(&kd, sizeof(int), 1, fortran_get_unit(97));
                            }
                        }
                    }
                }
                // do 6050 igrid=1,ngrid
                for (igrid = 1; igrid <= ngrid; igrid++) {
                    iskipz = 0;
                    if (igrid == 1) iskipz = 1;
                    nbl = nblg(igrid);
                    iem = iemg(igrid);
                    if (iseq != mseq && iem > 0) continue;
                    if (icgns == 1) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        idima = idim;
                        jdima = jdim;
                        kdima = kdim;
                    }
                    if (iem == 0) nbl = nbl + (mseq - iseq);
                    if (mblk2nd(nbl) == myid || myid == myhost) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        if (icgns == 1) {
                            if (nwork < (idima+1)*(jdima+1)*(kdima+1)) {
                                fprintf(fortran_get_unit(901), " not enough memory for cgns Q write.\n");
                                fprintf(fortran_get_unit(901), " nwork in wk=%6d.  Needed = %6d\n", nwork, (idima+1)*(jdima+1)*(kdima+1));
                                 termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                            }
                        } else {
                            idima = idim;
                            jdima = jdim;
                            kdima = kdim;
                        }
                        {
                            FortranArray4DRef<double> q_wv(w.data()+lq-1, jdim, kdim, idim, 5);
                            FortranArray4DRef<double> qj0_wv(w.data()+lqj0-1, kdim, idim-1, 5, 4);
                            FortranArray4DRef<double> qk0_wv(w.data()+lqk0-1, jdim, idim-1, 5, 4);
                            FortranArray4DRef<double> qi0_wv(w.data()+lqi0-1, jdim, kdim, 5, 4);
                            FortranArray3DRef<double> vist3d_wv(w.data()+lvis-1, jdim, kdim, idim);
                            FortranArray4DRef<double> tursav_wv(w.data()+lxib-1, jdim, kdim, idim, nummem);
                            FortranArray3DRef<double> smin_wv(w.data()+lsnk0-1, jdim-1, kdim-1, idim-1);
                            FortranArray3DRef<double> xjb_wv(w.data()+lsni0-1, jdim-1, kdim-1, idim-1);
                            FortranArray3DRef<double> xkb_wv(w.data()+lxkb-1, jdim-1, kdim-1, idim-1);
                            FortranArray3DRef<double> blnum_wv(w.data()+lnbl-1, jdim-1, kdim-1, idim-1);
                            FortranArray3DRef<double> cmuv_wv(w.data()+lcmuv-1, jdim-1, kdim-1, idim-1);
                            FortranArray3DRef<double> wk_wv(wk.data(), idima+1, jdima+1, kdima+1);
                            FortranArray3DRef<double> bcj_wv(w.data()+lbcj-1, kdim, idim-1, 2);
                            FortranArray3DRef<double> bck_wv(w.data()+lbck-1, jdim, idim-1, 2);
                            FortranArray3DRef<double> bci_wv(w.data()+lbci-1, jdim, kdim, 2);
                            FortranArray4DRef<double> vj0_wv(w.data()+lvj0-1, kdim, idim-1, 1, 4);
                            FortranArray4DRef<double> vk0_wv(w.data()+lvk0-1, jdim, idim-1, 1, 4);
                            FortranArray4DRef<double> vi0_wv(w.data()+lvi0-1, jdim, kdim, 1, 4);
                            FortranArray4DRef<double> tj0_wv(w.data()+ltj0-1, kdim, idim-1, nummem, 4);
                            FortranArray4DRef<double> tk0_wv(w.data()+ltk0-1, jdim, idim-1, nummem, 4);
                            FortranArray4DRef<double> ti0_wv(w.data()+lti0-1, jdim, kdim, nummem, 4);
                            FortranArray3DRef<double> blank_wv(w.data()+lblk-1, jdim, kdim, idim);
                            FortranArray4DRef<double> qavg_wv(w.data()+lqavg-1, jdim, kdim, idim, 5);
                            FortranArray4DRef<double> q2avg_wv(w.data()+lq2avg-1, jdim, kdim, idim, 5);
                            FortranArray3DRef<double> x_wv(w.data()+lx-1, jdim, kdim, idim);
                            FortranArray3DRef<double> y_wv(w.data()+ly-1, jdim, kdim, idim);
                            FortranArray3DRef<double> z_wv(w.data()+lz-1, jdim, kdim, idim);
                            int iover_wv = iovrlp(nbl);
                            wrest_ns::wrest(nbl, jdim, kdim, idim, q_wv, qj0_wv, qk0_wv,
                                qi0_wv, ncycmax, rms, clw, cdw, cdpw, cdvw, cxw, cyw, czw,
                                cmxw, cmyw, cmzw,
                                n_clcd, clcd, nblocks_clcd, blocks_clcd,
                                fmdotw, cftmomw, cftpw, cftvw, cfttotw,
                                rmstr, nneg, iskipz, vist3d_wv, tursav_wv,
                                smin_wv, xjb_wv, xkb_wv, blnum_wv, cmuv_wv,
                                thetay, maxbl, myid, myhost, mycomm, mblk2nd, igrid, wk_wv,
                                idima, jdima, kdima, bcj_wv, bck_wv, bci_wv, vj0_wv,
                                vk0_wv, vi0_wv, tj0_wv, tk0_wv, ti0_wv, blank_wv,
                                iwk, iwork, iover_wv, nou, bou, nbuf, ibufdim,
                                qavg_wv, q2avg_wv, x_wv, y_wv, z_wv, nummem);
                        }
                    }
                } // end do 6050 igrid


                // set flag for writing additional data 2nd order (in time) and/or dynamic mesh cases
                iflagg = 0;
                if (std::abs(ita) == 2 && iunst == 0) iflagg = 1;
                if (std::abs(ita) != 2 && iunst > 0)  iflagg = 2;
                if (std::abs(ita) == 2 && iunst > 0)  iflagg = 3;
                if (myid == myhost) {
                    if (icgns != 1) {
                        fwrite(&iflagg, sizeof(int), 1, fortran_get_unit(2));
                    }
                }

                // write 2nd order data
                if (iflagg == 1 || iflagg == 3) {
                    for (igrid = 1; igrid <= ngrid; igrid++) {
                        nbl = nblg(igrid);
                        iem = iemg(igrid);
                        if (iseq != mseq && iem > 0) continue;
                        if (icgns == 1) {
                            lead_ns::lead(nbl, lw, lw2, maxbl);
                            idima = idim; jdima = jdim; kdima = kdim;
                        }
                        if (iem == 0) nbl = nbl + (mseq - iseq);
                        if (mblk2nd(nbl) == myid || myid == myhost) {
                            lead_ns::lead(nbl, lw, lw2, maxbl);
                            if (icgns == 1) {
                                if (nwork < (idima+1)*(jdima+1)*(kdima+1)) {
                                    fprintf(fortran_get_unit(901), " not enough memory for cgns unsteady Q write.\n");
                                    fprintf(fortran_get_unit(901), " nwork in wk=%6d.  Needed = %6d\n", nwork, (idima+1)*(jdima+1)*(kdima+1));
                                     termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                                }
                            } else {
                                idima = idim; jdima = jdim; kdima = kdim;
                            }
                            {
                                FortranArray3DRef<double> x_wv(w.data()+lx-1, jdim, kdim, idim);
                                FortranArray3DRef<double> y_wv(w.data()+ly-1, jdim, kdim, idim);
                                FortranArray3DRef<double> z_wv(w.data()+lz-1, jdim, kdim, idim);
                                FortranArray3DRef<double> xnm2_wv(w.data()+lxnm2-1, jdim, kdim, idim);
                                FortranArray3DRef<double> ynm2_wv(w.data()+lynm2-1, jdim, kdim, idim);
                                FortranArray3DRef<double> znm2_wv(w.data()+lznm2-1, jdim, kdim, idim);
                                FortranArray4DRef<double> deltj_wv(w.data()+ldeltj-1, kdim, idim, 3, 2);
                                FortranArray4DRef<double> deltk_wv(w.data()+ldeltk-1, jdim, idim, 3, 2);
                                FortranArray4DRef<double> delti_wv(w.data()+ldelti-1, jdim, kdim, 3, 2);
                                FortranArray4DRef<double> qc0_wv(w.data()+lqc0-1, jdim, kdim, idim-1, 5);
                                FortranArray4DRef<double> tursav2_wv(w.data()+lxib2-1, jdim, kdim, idim, 2*nummem);
                                int nflagg_wv = 0;
                                int iuns_wv = 0;
                                wrestg_ns::wrestg(nbl, jdim, kdim, idim, x_wv, y_wv, z_wv,
                                    xnm2_wv, ynm2_wv, znm2_wv, deltj_wv, deltk_wv,
                                    delti_wv, qc0_wv, nflagg_wv, iuns_wv, utrans, vtrans, wtrans,
                                    omegax, omegay, omegaz, xorig, yorig, zorig,
                                    dxmx, dymx, dzmx, dthxmx, dthymx,
                                    dthzmx, thetax, thetay, thetaz, rfreqt,
                                    rfreqr, xorig0, yorig0, zorig0, time2,
                                    thetaxl, thetayl, thetazl, itrans, irotat, idefrm,
                                    utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae,
                                    xorgae, yorgae, zorgae, thtxae, thtyae, thtzae,
                                    rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf,
                                    kcsi, kcsf, freq, gmass, damp, x0, gf0, nmds, maxaes,
                                    aesrfdat, perturb, myhost, myid, mycomm, mblk2nd,
                                    maxbl, nsegdfrm, idfrmseg, iaesurf, maxsegdg, wk,
                                    nwork, idima, jdima, kdima, igrid, tursav2_wv, nummem);
                            }
                        }
                    } // end do 6060 igrid
                }

                // write dynamic mesh data
                if (iflagg == 2 || iflagg == 3) {
                    for (igrid = 1; igrid <= ngrid; igrid++) {
                        nbl = nblg(igrid);
                        iem = iemg(igrid);
                        if (iseq != mseq && iem > 0) continue;
                        if (icgns == 1) {
                            lead_ns::lead(nbl, lw, lw2, maxbl);
                            idima = idim; jdima = jdim; kdima = kdim;
                        }
                        if (iem == 0) nbl = nbl + (mseq - iseq);
                        if (mblk2nd(nbl) == myid || myid == myhost) {
                            lead_ns::lead(nbl, lw, lw2, maxbl);
                            if (icgns == 1) {
                                if (nwork < (idima+1)*(jdima+1)*(kdima+1)) {
                                    fprintf(fortran_get_unit(901), " not enough memory for cgns unsteady Q write.\n");
                                    fprintf(fortran_get_unit(901), " nwork in wk=%6d.  Needed = %6d\n", nwork, (idima+1)*(jdima+1)*(kdima+1));
                                     termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                                }
                            } else {
                                idima = idim; jdima = jdim; kdima = kdim;
                            }
                            iuns = std::max(itrans(nbl), std::max(irotat(nbl), idefrm(nbl)));
                            {
                                FortranArray3DRef<double> x_wv(w.data()+lx-1, jdim, kdim, idim);
                                FortranArray3DRef<double> y_wv(w.data()+ly-1, jdim, kdim, idim);
                                FortranArray3DRef<double> z_wv(w.data()+lz-1, jdim, kdim, idim);
                                FortranArray3DRef<double> xnm2_wv(w.data()+lxnm2-1, jdim, kdim, idim);
                                FortranArray3DRef<double> ynm2_wv(w.data()+lynm2-1, jdim, kdim, idim);
                                FortranArray3DRef<double> znm2_wv(w.data()+lznm2-1, jdim, kdim, idim);
                                FortranArray4DRef<double> deltj_wv(w.data()+ldeltj-1, kdim, idim, 3, 2);
                                FortranArray4DRef<double> deltk_wv(w.data()+ldeltk-1, jdim, idim, 3, 2);
                                FortranArray4DRef<double> delti_wv(w.data()+ldelti-1, jdim, kdim, 3, 2);
                                FortranArray4DRef<double> qc0_wv(w.data()+lqc0-1, jdim, kdim, idim-1, 5);
                                FortranArray4DRef<double> tursav2_wv(w.data()+lxib2-1, jdim, kdim, idim, 2*nummem);
                                int nflagg_wv = 1;
                                wrestg_ns::wrestg(nbl, jdim, kdim, idim, x_wv, y_wv, z_wv,
                                    xnm2_wv, ynm2_wv, znm2_wv, deltj_wv,
                                    deltk_wv, delti_wv, qc0_wv, nflagg_wv, iuns, utrans,
                                    vtrans, wtrans, omegax, omegay, omegaz, xorig,
                                    yorig, zorig, dxmx, dymx, dzmx, dthxmx, dthymx,
                                    dthzmx, thetax, thetay, thetaz, rfreqt,
                                    rfreqr, xorig0, yorig0, zorig0, time2,
                                    thetaxl, thetayl, thetazl, itrans, irotat, idefrm,
                                    utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae,
                                    xorgae, yorgae, zorgae, thtxae, thtyae, thtzae,
                                    rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf,
                                    kcsi, kcsf, freq, gmass, damp, x0, gf0, nmds, maxaes,
                                    aesrfdat, perturb, myhost, myid, mycomm, mblk2nd,
                                    maxbl, nsegdfrm, idfrmseg, iaesurf, maxsegdg, wk,
                                    nwork, idima, jdima, kdima, igrid, tursav2_wv, nummem);
                            }
                        }
                    } // end do 6070 igrid
                }

                // write modal data if aeroelastic
                if (myid == myhost) {
                    if (naesrf > 0) {
                        if (icgns != 1) {
                            for (nn = 1; nn <= ntt; nn++) {
                                fwrite(&timekeep(nn), sizeof(double), 1, fortran_get_unit(2));
                            }
                            for (iaes = 1; iaes <= naesrf; iaes++) {
                                nmodes = (int)aesrfdat(5, iaes);
                                for (n = 1; n <= 2*nmodes; n++) fwrite(&xxn(n,iaes), sizeof(double), 1, fortran_get_unit(2));
                                for (n = 1; n <= 2*nmodes; n++) fwrite(&gforcn(n,iaes), sizeof(double), 1, fortran_get_unit(2));
                                for (n = 1; n <= 2*nmodes; n++) fwrite(&gforcnm(n,iaes), sizeof(double), 1, fortran_get_unit(2));
                                for (n = 1; n <= nmodes; n++) {
                                    for (int ll = 1; ll <= 3; ll++) {
                                        for (nn = 1; nn <= ntt; nn++) {
                                            fwrite(&aehist(nn,ll,n,iaes), sizeof(double), 1, fortran_get_unit(2));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (icgns != 1) {
                        fortran_close_unit(2);
                        if (ndwrt != 0) {
                            fortran_close_unit(98);
                        }
                    }
                }
            } // end if nt >= ntstep || nwrest
        } // end if level == levt_c
    } // end if dt > 0


    // output new plot3d file every movie time steps
    movabs = std::abs(movie);
    if (ipertavg == 0) {
        if ((float)dt > 0.0f && movabs > 0) {
            if ((movabs > 0 && (nt == nt/movabs*movabs)) || (iteravg > 0)) {
                if (nt == nt/movabs*movabs) {
                    lhdr   = 0;
                    iwk1   = 1;
                    iwk2   = iwk1 + 3*nplots;
                    iwk3   = iwk2 + maxbl;
                    iwork1 = iwork - iwk3 + 1;
                    if (iwork1 <= 0) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine qout");
                         termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
                    for (iii = 1; iii <= nwork; iii++) wk(iii) = 0.0;
                    {
                        FortranArray2DRef<int> iwk_iwk1(iwk.data()+iwk1-1, 3, nplots);
                        FortranArray1DRef<int> iwk_iwk2 = iwk.slice(iwk2, iwk.ubound(1));
                        FortranArray1DRef<int> iwk_iwk3 = iwk.slice(iwk3, iwk.ubound(1));
                        qout_ns::qout(iseq, lw, lw2, w, mgwk, wk, nwork,
                            nplots, iovrlp, iibg, kkbg, jjbg, ibcg, lbg,
                            ibpntsg, qb, lwdat, nbci0, nbcj0, nbck0,
                            nbcidim, nbcjdim, nbckdim, jbcinfo, kbcinfo,
                            ibcinfo, bcfilei, bcfilej, bcfilek, itrans,
                            irotat, idefrm, nblock, levelg, igridg, iviscg, jdimg,
                            kdimg, idimg, nblg, clw, ncycmax, nplot3d,
                            inpl3d, ip3dsurf, nprint, inpr, iadvance, mycomm,
                            myid, myhost, mblk2nd, nou, bou, nbuf, ibufdim, maxbl,
                            maxgr, maxseg, iitot, jsg, ksg, isg, jeg, keg, ieg,
                            ninter, windex, iindex, nblkpt, intmax, nsub1, maxxe,
                            nblk, nbli, limblk, isva, nblon, mxbli, thetay,
                            iwk_iwk1, iwk_iwk2, iwk_iwk3, iwork1, xorig, yorig,
                            zorig, period_miss, geom_miss, epsc0, epsrot,
                            isav_blk, isav_pat, isav_pat_b, isav_emb, isav_prd,
                            lbcprd, lbcemb, dthetxx, dthetyy, dthetzz, nblcg, lfgm,
                            istat2_bl, istat2_pa, istat2_pe, istat2_em,
                            istat_size, vormax, ivmax, jvmax, kvmax, bcfiles,
                            mxbcfil, iprnsurf, nummem);
                    }
                    for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
                    for (iii = 1; iii <= nwork; iii++) wk(iii) = 0.0;
                    nframes = nframes + 1;
                }
            }
        }
    }
    if (ipertavg == 1 || ipertavg == 2) {
        if (((float)dt > 0.0f && movabs > 0) || (iteravg > 0)) {
            if ((movabs > 0 && (nt == nt/movabs*movabs)) || (iteravg > 0)) {
                lhdr   = 0;
                iwk1   = 1;
                iwk2   = iwk1 + 3*nplots;
                iwk3   = iwk2 + maxbl;
                iwork1 = iwork - iwk3 + 1;
                if (iwork1 <= 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine qoutavg");
                     termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                }
                for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
                for (iii = 1; iii <= nwork; iii++) wk(iii) = 0.0;
                {
                    FortranArray2DRef<int> iwk_iwk1(iwk.data()+iwk1-1, 3, nplots);
                    FortranArray1DRef<int> iwk_iwk2 = iwk.slice(iwk2, iwk.ubound(1));
                    FortranArray1DRef<int> iwk_iwk3 = iwk.slice(iwk3, iwk.ubound(1));
                    qoutavg_ns::qoutavg(iseq, lw, lw2, w, mgwk, wk, nwork,
                        nplots, iovrlp, iibg, kkbg, jjbg, ibcg, lbg,
                        ibpntsg, qb, lwdat, nbci0, nbcj0, nbck0,
                        nbcidim, nbcjdim, nbckdim, jbcinfo, kbcinfo,
                        ibcinfo, bcfilei, bcfilej, bcfilek, itrans,
                        irotat, idefrm, nblock, levelg, igridg, iviscg, jdimg,
                        kdimg, idimg, nblg, clw, ncycmax, nplot3d,
                        inpl3d, ip3dsurf, nprint, inpr, iadvance, mycomm,
                        myid, myhost, mblk2nd, nou, bou, nbuf, ibufdim, maxbl,
                        maxgr, maxseg, iitot, jsg, ksg, isg, jeg, keg, ieg,
                        ninter, windex, iindex, nblkpt, intmax, nsub1, maxxe,
                        nblk, nbli, limblk, isva, nblon, mxbli, thetay,
                        iwk_iwk1, iwk_iwk2, iwk_iwk3, iwork1, xorig, yorig,
                        zorig, period_miss, geom_miss, epsc0, epsrot,
                        isav_blk, isav_pat, isav_pat_b, isav_emb, isav_prd,
                        lbcprd, lbcemb, dthetxx, dthetyy, dthetzz, nblcg, lfgm,
                        istat2_bl, istat2_pa, istat2_pe, istat2_em,
                        istat_size, vormax, ivmax, jvmax, kvmax, bcfiles,
                        mxbcfil, iprnsurf, nt, movabs, nummem);
                }
                for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
                for (iii = 1; iii <= nwork; iii++) wk(iii) = 0.0;
                nframes = nframes + 1;
            }
        }
    }


    // coarse movie output
    if (icoarsemovie != 0) {
        iflag_coarse = icoarsemovie;
        mov_coarse   = icoarsemovie;
        iinc_coarse  = cmn_ginfo2.inc_coarse[0];  // inc_coarse(1)
        jinc_coarse  = cmn_ginfo2.inc_coarse[1];  // inc_coarse(2)
        kinc_coarse  = cmn_ginfo2.inc_coarse[2];  // inc_coarse(3)

        if ((float)dt > 0.0f && iflag_coarse > 0) {
            if (ntt == ntt/mov_coarse*mov_coarse) {
                lhdr   = 0;
                iwk1   = 1;
                iwk2   = iwk1 + 3*nblocks/lglobal;
                iwk3   = iwk2 + maxbl;
                iwork1 = iwork - iwk3 + 1;
                if (iwork1 <= 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine qout_coarse");
                     termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                }
                for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
                for (iii = 1; iii <= nwork; iii++) wk(iii) = 0.0;
                if (nt == 1) { icallcrs = 0; icall2d = 0; }
                {
                    FortranArray2DRef<int> iwk_iwk1(iwk.data()+iwk1-1, 3, nplots);
                    FortranArray1DRef<int> iwk_iwk2 = iwk.slice(iwk2, iwk.ubound(1));
                    FortranArray1DRef<int> iwk_iwk3 = iwk.slice(iwk3, iwk.ubound(1));
                    qout_coarse_ns::qout_coarse(iseq, lw, lw2, w, mgwk, wk, nwork,
                        nplots, iovrlp, iibg, kkbg, jjbg, ibcg, lbg,
                        ibpntsg, qb, lwdat, nbci0, nbcj0, nbck0,
                        nbcidim, nbcjdim, nbckdim, jbcinfo, kbcinfo,
                        ibcinfo, bcfilei, bcfilej, bcfilek, itrans,
                        irotat, idefrm, nblock, levelg, igridg, iviscg, jdimg,
                        kdimg, idimg, nblg, clw, ncycmax, nplot3d,
                        inpl3d, ip3dsurf, nprint, inpr, iadvance, mycomm,
                        myid, myhost, mblk2nd, nou, bou, nbuf, ibufdim, maxbl,
                        maxgr, maxseg, iitot, jsg, ksg, isg, jeg, keg, ieg,
                        ninter, windex, iindex, nblkpt, intmax, nsub1, maxxe,
                        nblk, nbli, limblk, isva, nblon, mxbli, thetay,
                        iwk_iwk1, iwk_iwk2, iwk_iwk3, iwork1, xorig, yorig,
                        zorig, period_miss, geom_miss, epsc0, epsrot,
                        isav_blk, isav_pat, isav_pat_b, isav_emb, isav_prd,
                        lbcprd, lbcemb, dthetxx, dthetyy, dthetzz, nblcg, lfgm,
                        istat2_bl, istat2_pa, istat2_pe, istat2_em,
                        istat_size, vormax, ivmax, jvmax, kvmax, bcfiles,
                        mxbcfil, iprnsurf, nt,
                        mov_coarse, iinc_coarse, jinc_coarse, kinc_coarse,
                        nummem);
                }
                for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
                for (iii = 1; iii <= nwork; iii++) wk(iii) = 0.0;
            }
        }
    }

    // 2D movie output
    if (i2dmovie != 0) {
        iflag_2d = i2dmovie;
        mov_2d   = i2dmovie;
        iinc_2d  = cmn_ginfo2.inc_2d[0];  // inc_2d(1)
        jinc_2d  = cmn_ginfo2.inc_2d[1];  // inc_2d(2)
        kinc_2d  = cmn_ginfo2.inc_2d[2];  // inc_2d(3)

        if ((float)dt > 0.0f && iflag_2d > 0) {
            if (ntt == ntt/mov_2d*mov_2d) {
                lhdr   = 0;
                iwk1   = 1;
                iwk2   = iwk1 + 3*nblocks/lglobal;
                iwk3   = iwk2 + maxbl;
                iwork1 = iwork - iwk3 + 1;
                if (iwork1 <= 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "stopping...not enough integer work space for subroutine qout_2d");
                     termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                }
                for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
                for (iii = 1; iii <= nwork; iii++) wk(iii) = 0.0;
                {
                    FortranArray2DRef<int> iwk_iwk1(iwk.data()+iwk1-1, 3, nplots);
                    FortranArray1DRef<int> iwk_iwk2 = iwk.slice(iwk2, iwk.ubound(1));
                    FortranArray1DRef<int> iwk_iwk3 = iwk.slice(iwk3, iwk.ubound(1));
                    qout_2d_ns::qout_2d(iseq, lw, lw2, w, mgwk, wk, nwork,
                        nplots, iovrlp, iibg, kkbg, jjbg, ibcg, lbg,
                        ibpntsg, qb, lwdat, nbci0, nbcj0, nbck0,
                        nbcidim, nbcjdim, nbckdim, jbcinfo, kbcinfo,
                        ibcinfo, bcfilei, bcfilej, bcfilek, itrans,
                        irotat, idefrm, nblock, levelg, igridg, iviscg, jdimg,
                        kdimg, idimg, nblg, clw, ncycmax, nplot3d,
                        inpl3d, ip3dsurf, nprint, inpr, iadvance, mycomm,
                        myid, myhost, mblk2nd, nou, bou, nbuf, ibufdim, maxbl,
                        maxgr, maxseg, iitot, jsg, ksg, isg, jeg, keg, ieg,
                        ninter, windex, iindex, nblkpt, intmax, nsub1, maxxe,
                        nblk, nbli, limblk, isva, nblon, mxbli, thetay,
                        iwk_iwk1, iwk_iwk2, iwk_iwk3, iwork1, xorig, yorig,
                        zorig, period_miss, geom_miss, epsc0, epsrot,
                        isav_blk, isav_pat, isav_pat_b, isav_emb, isav_prd,
                        lbcprd, lbcemb, dthetxx, dthetyy, dthetzz, nblcg, lfgm,
                        istat2_bl, istat2_pa, istat2_pe, istat2_em,
                        istat_size, vormax, ivmax, jvmax, kvmax, bcfiles,
                        mxbcfil, iprnsurf, nt,
                        mov_2d, iinc_2d, jinc_2d, kinc_2d, nummem);
                }
                for (iii = 1; iii <= iwk3; iii++) iwk(iii) = 0;
                for (iii = 1; iii <= nwork; iii++) wk(iii) = 0.0;
            }
        }
    }

    // label 5990: meshdef bypass target
    label_5990:;

    // check if user has altered ntstep via "stop" file
    if (nt >= ntstep) goto label_6011;

    } // end do 6000 nt (time loop)

    label_6011:;

    if (myid == myhost) {
        fprintf(fortran_get_unit(11), "\n***** ENDING TIME ADVANCEMENT, iseq =%2d *****\n\n", iseq);
    }

} // end mgblk

// ===== Wrapper functions for other mgblk_ns procedures =====
// These delegate to their respective implementation namespaces.

void add2x(FortranArray4DRef<double> q, FortranArray4DRef<double> qc, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, FortranArray4DRef<double> q1, FortranArray4DRef<double> dq, FortranArray4DRef<double> wq, FortranArray3DRef<double> wqj, FortranArray3DRef<double> wqjk, int& js, int& ks, int& is, int& je, int& ke, int& ie, int& ipass, int& nbl, int& nblc, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& ll, int& myid) {
    add2x_ns::add2x(q, qc, jdim, kdim, idim, jj2, kk2, ii2, q1, dq, wq, wqj, wqjk, js, ks, is, je, ke, ie, ipass, nbl, nblc, nou, bou, nbuf, ibufdim, ll, myid);
}

void addx(FortranArray4DRef<double> q, FortranArray4DRef<double> qq, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, FortranArray4DRef<double> q1, FortranArray4DRef<double> dq, FortranArray4DRef<double> wq, FortranArray4DRef<double> wqj, int& nbl, FortranArray3DRef<double> blank, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& ll, int& myid) {
    addx_ns::addx(q, qq, jdim, kdim, idim, jj2, kk2, ii2, q1, dq, wq, wqj, nbl, blank, nou, bou, nbuf, ibufdim, ll, myid);
}

void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, double& cl, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, FortranArray3DRef<int> lwdat, int& myid, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, int& nummem) {
    bc_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, cl, nou, bou, nbuf, ibufdim, maxbl, maxgr, maxseg, itrans, irotat, idefrm, igridg, nblg, nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim, ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek, lwdat, myid, idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
}

void bc_blkint(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& maxbl, int& maxgr, int& mxbli, FortranArray1DRef<int> iadvance, FortranArray1DRef<double> geom_miss, double& epsc0, FortranArray2DRef<int> nblk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> isav_blk, FortranArray1DRef<int> ireq_ar, FortranArray1DRef<int> index_ar, FortranArray1DRef<int> ireq_snd, FortranArray2DRef<int> keep_trac, FortranArray1DRef<int> keep_trac2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& myhost, int& mycomm, FortranArray2DRef<int> istat2, int& istat_size, int& nummem) {
    bc_blkint_ns::bc_blkint(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, maxbl, maxgr, mxbli, iadvance, geom_miss, epsc0, nblk, nbli, limblk, isva, nblon, jdimg, kdimg, idimg, mblk2nd, isav_blk, ireq_ar, index_ar, ireq_snd, keep_trac, keep_trac2, nou, bou, nbuf, ibufdim, myid, myhost, mycomm, istat2, istat_size, nummem);
}

void bc_embed(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& maxbl, int& maxgr, int& lbcemb, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray2DRef<int> isav_emb, FortranArray1DRef<int> ireq_ar, FortranArray1DRef<int> ireq_snd, FortranArray1DRef<int> index_ar, FortranArray2DRef<int> keep_trac, FortranArray1DRef<int> keep_trac2, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray2DRef<int> iviscg, FortranArray2DRef<int> istat2, int& istat_size, int& nummem) {
    bc_embed_ns::bc_embed(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, maxbl, maxgr, lbcemb, iadvance, idimg, jdimg, kdimg, isav_emb, ireq_ar, ireq_snd, index_ar, keep_trac, keep_trac2, myid, myhost, mycomm, mblk2nd, nou, bou, nbuf, ibufdim, iviscg, istat2, istat_size, nummem);
}

void bc_patch(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& maxbl, int& maxgr, int& intmax, int& nsub1, int& maxxe, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, int& ninter, FortranArray2DRef<double> windex, FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt, FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz, FortranArray2DRef<int> isav_pat, FortranArray3DRef<int> isav_pat_b, FortranArray1DRef<int> ireq_ar, FortranArray1DRef<int> ireq_snd, FortranArray1DRef<int> index_ar, FortranArray3DRef<int> keep_trac, FortranArray2DRef<int> keep_trac2, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray2DRef<int> istat2, int& istat_size, int& nummem) {
    bc_patch_ns::bc_patch(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, maxbl, maxgr, intmax, nsub1, maxxe, iadvance, jdimg, kdimg, idimg, ninter, windex, iindex, nblkpt, dthetxx, dthetyy, dthetzz, isav_pat, isav_pat_b, ireq_ar, ireq_snd, index_ar, keep_trac, keep_trac2, myid, myhost, mycomm, mblk2nd, nou, bou, nbuf, ibufdim, istat2, istat_size, nummem);
}

void bc_period(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& maxbl, int& maxgr, int& maxseg, FortranArray1DRef<int> iadvance, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, FortranArray3DRef<int> lwdat, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, int& lbcprd, FortranArray2DRef<int> isav_prd, FortranArray1DRef<double> period_miss, double& epsrot, FortranArray1DRef<int> ireq_ar, FortranArray1DRef<int> index_ar, FortranArray1DRef<int> ireq_snd, FortranArray2DRef<int> keep_trac, FortranArray1DRef<int> keep_trac2, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray2DRef<int> istat2, int& istat_size, FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, int& nummem) {
    bc_period_ns::bc_period(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, maxbl, maxgr, maxseg, iadvance, bcfilei, bcfilej, bcfilek, lwdat, xorig, yorig, zorig, jdimg, kdimg, idimg, lbcprd, isav_prd, period_miss, epsrot, ireq_ar, index_ar, ireq_snd, keep_trac, keep_trac2, myid, myhost, mycomm, mblk2nd, nou, bou, nbuf, ibufdim, istat2, istat_size, bcfiles, mxbcfil, nummem);
}

void coll2q(FortranArray4DRef<double> q, FortranArray4DRef<double> qc, FortranArray3DRef<double> vol, FortranArray3DRef<double> volc, int& jdim, int& kdim, int& idim, int& jc, int& kc, int& ic, FortranArray4DRef<double> res, FortranArray4DRef<double> qr, FortranArray4DRef<double> qw, int& js, int& ks, int& is, int& je, int& ke, int& ie, int& nbl, int& nblc, FortranArray3DRef<double> vistf, FortranArray3DRef<double> vistc, FortranArray4DRef<double> tursavf, FortranArray4DRef<double> tursavc, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& nummem) {
    coll2q_ns::coll2q(q, qc, vol, volc, jdim, kdim, idim, jc, kc, ic, res, qr, qw, js, ks, is, je, ke, ie, nbl, nblc, vistf, vistc, tursavf, tursavc, nou, bou, nbuf, ibufdim, nummem);
}

void collq(FortranArray4DRef<double> q, FortranArray4DRef<double> qq, FortranArray3DRef<double> vol, FortranArray3DRef<double> volc, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, FortranArray4DRef<double> res, FortranArray4DRef<double> qr, FortranArray4DRef<double> qw, FortranArray3DRef<double> vistf, FortranArray3DRef<double> vistc, FortranArray4DRef<double> tursavf, FortranArray4DRef<double> tursavc, int& nbl, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& nummem) {
    collq_ns::collq(q, qq, vol, volc, jdim, kdim, idim, jj2, kk2, ii2, res, qr, qw, vistf, vistc, tursavf, tursavc, nbl, nou, bou, nbuf, ibufdim, nummem);
}

void dynptch(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> work, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& ncall, int& maxgr, int& maxbl, int& msub1, int& intmx, int& mxxe, int& mptch, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, int& nblock, int& ngrid, FortranArray1DRef<int> levelg, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray2DRef<double> windex, int& ninter, FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt, FortranArray2DRef<double> windx, int& nintr, FortranArray2DRef<int> iindx, FortranArray1DRef<int> mblkpt, FortranArray1DRef<int> llimit, FortranArray1DRef<int> iitmax, FortranArray1DRef<int> mmcxie, FortranArray1DRef<int> mmceta, FortranArray1DRef<int> ncheck, FortranArray1DRef<int> iifit, FortranArray1DRef<int> iic0, FortranArray1DRef<int> iiorph, FortranArray1DRef<int> iitoss, FortranArray1DRef<int> ifiner, FortranArray2DRef<double> factjlo, FortranArray2DRef<double> factjhi, FortranArray2DRef<double> factklo, FortranArray2DRef<double> factkhi, FortranArray2DRef<double> dx, FortranArray2DRef<double> dy, FortranArray2DRef<double> dz, FortranArray2DRef<double> dthetx, FortranArray2DRef<double> dthety, FortranArray2DRef<double> dthetz, FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz, FortranArray2DRef<int> isav_dpat, FortranArray3DRef<int> isav_dpat_b, int& intmax, int& maxxe, int& nsub1, FortranArray2DRef<int> lw_temp, FortranArray1DRef<int> ireq_ar, FortranArray1DRef<int> lout, FortranArray1DRef<int> ifrom, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, FortranArray1DRef<int> iiint1, FortranArray1DRef<int> iiint2, FortranArray1DRef<int> nblk1, FortranArray1DRef<int> nblk2, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, FortranArray1DRef<int> jte, FortranArray1DRef<int> kte, FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm, FortranArray3DRef<double> xte, FortranArray3DRef<double> yte, FortranArray3DRef<double> zte, FortranArray3DRef<double> xmi, FortranArray3DRef<double> ymi, FortranArray3DRef<double> zmi, FortranArray3DRef<double> xmie, FortranArray3DRef<double> ymie, FortranArray3DRef<double> zmie, FortranArray3DRef<double> sxie, FortranArray3DRef<double> seta, FortranArray2DRef<double> sxie2, FortranArray2DRef<double> seta2, FortranArray2DRef<double> xie2s, FortranArray2DRef<double> eta2s, FortranArray1DRef<double> temp, FortranArray2DRef<double> x2, FortranArray2DRef<double> y2, FortranArray2DRef<double> z2, FortranArray2DRef<double> x1, FortranArray2DRef<double> y1, FortranArray2DRef<double> z1, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray1DRef<int> igridg, FortranArray1DRef<int> iemg) {
    dynptch_ns::dynptch(lw, lw2, work, mgwk, wk, nwork, ncall, maxgr, maxbl, msub1, intmx, mxxe, mptch, jdimg, kdimg, idimg, xorig, yorig, zorig, nblock, ngrid, levelg, ncgg, nblg, windex, ninter, iindex, nblkpt, windx, nintr, iindx, mblkpt, llimit, iitmax, mmcxie, mmceta, ncheck, iifit, iic0, iiorph, iitoss, ifiner, factjlo, factjhi, factklo, factkhi, dx, dy, dz, dthetx, dthety, dthetz, dthetxx, dthetyy, dthetzz, isav_dpat, isav_dpat_b, intmax, maxxe, nsub1, lw_temp, ireq_ar, lout, ifrom, xif1, xif2, etf1, etf2, jjmax1, kkmax1, iiint1, iiint2, nblk1, nblk2, jimage, kimage, jte, kte, jmm, kmm, xte, yte, zte, xmi, ymi, zmi, xmie, ymie, zmie, sxie, seta, sxie2, seta2, xie2s, eta2s, temp, x2, y2, z2, x1, y1, z1, myid, myhost, mycomm, mblk2nd, nou, bou, nbuf, ibufdim, igridg, iemg);
}

void intrbc(FortranArray4DRef<double> q, int& jdim, int& kdim, int& idim, int& nbl, int& ldim, int& maxbl, int& iitot, FortranArray1DRef<int> lig, FortranArray1DRef<int> iipntsg, FortranArray1DRef<double> dxintg, FortranArray1DRef<double> dyintg, FortranArray1DRef<double> dzintg, FortranArray1DRef<int> iiig, FortranArray1DRef<int> jjig, FortranArray1DRef<int> kkig, FortranArray3DRef<double> qb, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qq, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& icorr, int& iset) {
    intrbc_ns::intrbc(q, jdim, kdim, idim, nbl, ldim, maxbl, iitot, lig, iipntsg, dxintg, dyintg, dzintg, iiig, jjig, kkig, qb, qj0, qk0, qi0, qq, bcj, bck, bci, nou, bou, nbuf, ibufdim, icorr, iset);
}

void mgbl(int& mwork, int& mworki, int& nplots, int& minnode, int& iitot, int& intmax, int& maxxe, int& mxbli, int& nsub1, int& lbcprd, int& lbcemb, int& lbcrad, int& maxbl, int& maxgr, int& maxseg, int& maxcs, int& ncycmax, int& intmx, int& mxxe, int& mptch, int& msub1, int& ibufdim, int& nbuf, int& mxbcfil, int& istat_size, FortranArray1DRef<double> work, FortranArray1DRef<int> iwork, FortranArray3DRef<int> lwdat, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, FortranArray2DRef<int> istat2_bl, FortranArray2DRef<int> istat2_pa, FortranArray2DRef<int> istat2_em, FortranArray2DRef<int> istat2_pe, FortranArray2DRef<int> nblk, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, FortranArray1DRef<double> resmx, FortranArray1DRef<int> imx, FortranArray1DRef<int> jmx, FortranArray1DRef<int> kmx, FortranArray1DRef<double> vormax, FortranArray1DRef<int> ivmax, FortranArray1DRef<int> jvmax, FortranArray1DRef<int> kvmax, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray1DRef<int> iovrlp, FortranArray3DRef<double> qb, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg, FortranArray1DRef<double> dxintg, FortranArray1DRef<double> dyintg, FortranArray1DRef<double> dzintg, FortranArray1DRef<int> iiig, FortranArray1DRef<int> jjig, FortranArray1DRef<int> kkig, FortranArray2DRef<double> rkap0g, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<int> idiagg, FortranArray1DRef<int> nblcg, FortranArray2DRef<int> idegg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray2DRef<int> mit, FortranArray1DRef<int> ilamlog, FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamlog, FortranArray1DRef<int> klamhig, FortranArray2DRef<int> iwfg, FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx, FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz, FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr, FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0, FortranArray1DRef<double> time2, FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg, FortranArray2DRef<int> inpl3d, FortranArray2DRef<int> inpr, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> iforce, FortranArray1DRef<double> rms, FortranArray1DRef<double> clw, FortranArray1DRef<double> cdw, FortranArray1DRef<double> cdpw, FortranArray1DRef<double> cdvw, FortranArray1DRef<double> cxw, FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw, FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw, FortranArray1DRef<double> cmzw, int& n_clcd, FortranArray3DRef<double> clcd, int& nblocks_clcd, FortranArray2DRef<int> blocks_clcd, FortranArray1DRef<double> chdw, FortranArray1DRef<double> swetw, FortranArray1DRef<double> fmdotw, FortranArray1DRef<double> cfttotw, FortranArray1DRef<double> cftmomw, FortranArray1DRef<double> cftpw, FortranArray1DRef<double> cftvw, FortranArray1DRef<double> swett, FortranArray1DRef<double> clt, FortranArray1DRef<double> cdt, FortranArray1DRef<double> cxt, FortranArray1DRef<double> cyt, FortranArray1DRef<double> czt, FortranArray1DRef<double> cmxt, FortranArray1DRef<double> cmyt, FortranArray1DRef<double> cmzt, FortranArray1DRef<double> cdpt, FortranArray1DRef<double> cdvt, FortranArray1DRef<double> sx, FortranArray1DRef<double> sy, FortranArray1DRef<double> sz, FortranArray1DRef<double> stot, FortranArray1DRef<double> pav, FortranArray1DRef<double> ptav, FortranArray1DRef<double> tav, FortranArray1DRef<double> ttav, FortranArray1DRef<double> xmav, FortranArray1DRef<double> fmdot, FortranArray1DRef<double> cfxp, FortranArray1DRef<double> cfyp, FortranArray1DRef<double> cfzp, FortranArray1DRef<double> cfdp, FortranArray1DRef<double> cflp, FortranArray1DRef<double> cftp, FortranArray1DRef<double> cfxv, FortranArray1DRef<double> cfyv, FortranArray1DRef<double> cfzv, FortranArray1DRef<double> cfdv, FortranArray1DRef<double> cflv, FortranArray1DRef<double> cftv, FortranArray1DRef<double> cfxmom, FortranArray1DRef<double> cfymom, FortranArray1DRef<double> cfzmom, FortranArray1DRef<double> cfdmom, FortranArray1DRef<double> cflmom, FortranArray1DRef<double> cftmom, FortranArray1DRef<double> cfxtot, FortranArray1DRef<double> cfytot, FortranArray1DRef<double> cfztot, FortranArray1DRef<double> cfdtot, FortranArray1DRef<double> cfltot, FortranArray1DRef<double> cfttot, FortranArray2DRef<int> icsinfo, FortranArray2DRef<double> windex, FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt, FortranArray2DRef<double> windx, FortranArray2DRef<int> iindx, FortranArray1DRef<int> llimit, FortranArray1DRef<int> iitmax, FortranArray1DRef<int> mmcxie, FortranArray1DRef<int> mmceta, FortranArray1DRef<int> ncheck, FortranArray1DRef<int> iifit, FortranArray1DRef<int> mblkpt, FortranArray1DRef<int> iic0, FortranArray1DRef<int> iiorph, FortranArray1DRef<int> iitoss, FortranArray1DRef<int> ifiner, FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz, FortranArray2DRef<double> dx, FortranArray2DRef<double> dy, FortranArray2DRef<double> dz, FortranArray2DRef<double> dthetx, FortranArray2DRef<double> dthety, FortranArray2DRef<double> dthetz, FortranArray1DRef<int> lout, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, FortranArray1DRef<int> iiint1, FortranArray1DRef<int> iiint2, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, FortranArray1DRef<int> jte, FortranArray1DRef<int> kte, FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm, FortranArray1DRef<int> nblk1, FortranArray1DRef<int> nblk2, FortranArray3DRef<double> xte, FortranArray3DRef<double> yte, FortranArray3DRef<double> zte, FortranArray3DRef<double> xmi, FortranArray3DRef<double> ymi, FortranArray3DRef<double> zmi, FortranArray3DRef<double> xmie, FortranArray3DRef<double> ymie, FortranArray3DRef<double> zmie, FortranArray3DRef<double> sxie, FortranArray3DRef<double> seta, FortranArray2DRef<double> sxie2, FortranArray2DRef<double> seta2, FortranArray2DRef<double> xie2s, FortranArray2DRef<double> eta2s, FortranArray1DRef<double> temp, FortranArray2DRef<double> x2, FortranArray2DRef<double> y2, FortranArray2DRef<double> z2, FortranArray2DRef<double> x1, FortranArray2DRef<double> y1, FortranArray2DRef<double> z1, FortranArray2DRef<double> factjlo, FortranArray2DRef<double> factjhi, FortranArray2DRef<double> factklo, FortranArray2DRef<double> factkhi, FortranArray1DRef<int> ifrom, FortranArray1DRef<double> geom_miss, FortranArray1DRef<double> period_miss, FortranArray2DRef<int> isav_blk, FortranArray2DRef<int> isav_prd, FortranArray2DRef<int> isav_pat, FortranArray3DRef<int> isav_pat_b, FortranArray2DRef<int> isav_dpat, FortranArray3DRef<int> isav_dpat_b, FortranArray2DRef<int> isav_emb, FortranArray1DRef<int> mblk2nd, int& ntr, FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae, FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae, FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass, FortranArray2DRef<double> damp, FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0, int& nmds, int& maxaes, FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb, int& memuse, FortranArray1DRef<char[80]> bcfiles, FortranArray3DRef<int> islavept, int& nslave, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray3DRef<double> bmat, FortranArray3DRef<double> stm, FortranArray3DRef<double> stmi, FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn, FortranArray2DRef<double> gforcn, FortranArray2DRef<double> gforcnm, FortranArray2DRef<double> gforcs, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg, int& nmaster, FortranArray4DRef<double> aehist, FortranArray1DRef<double> timekeep, FortranArray2DRef<double> xorgae0, FortranArray2DRef<double> yorgae0, FortranArray2DRef<double> zorgae0, FortranArray2DRef<int> icouple, FortranArray2DRef<int> nblelst, FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax, FortranArray1DRef<int> kskmax, FortranArray1DRef<double> ue, FortranArray1DRef<int> irdrea) {
    mgbl_ns::mgbl(mwork, mworki, nplots, minnode, iitot, intmax, maxxe, mxbli, nsub1, lbcprd, lbcemb, lbcrad, maxbl, maxgr, maxseg, maxcs, ncycmax, intmx, mxxe, mptch, msub1, ibufdim, nbuf, mxbcfil, istat_size, work, iwork, lwdat, lw, lw2, nou, bou, bcfilei, bcfilej, bcfilek, istat2_bl, istat2_pa, istat2_em, istat2_pe, nblk, limblk, isva, nblon, resmx, imx, jmx, kmx, vormax, ivmax, jvmax, kvmax, lig, lbg, iovrlp, qb, ibpntsg, iipntsg, iibg, kkbg, jjbg, ibcg, dxintg, dyintg, dzintg, iiig, jjig, kkig, rkap0g, levelg, igridg, iflimg, ifdsg, iviscg, jdimg, kdimg, idimg, idiagg, nblcg, idegg, jsg, ksg, isg, jeg, keg, ieg, mit, ilamlog, ilamhig, jlamlog, jlamhig, klamlog, klamhig, iwfg, utrans, vtrans, wtrans, omegax, omegay, omegaz, xorig, yorig, zorig, dxmx, dymx, dzmx, dthxmx, dthymx, dthzmx, thetax, thetay, thetaz, rfreqt, rfreqr, xorig0, yorig0, zorig0, time2, thetaxl, thetayl, thetazl, itrans, irotat, idefrm, bcvali, bcvalj, bcvalk, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim, ibcinfo, jbcinfo, kbcinfo, ncgg, nblg, iemg, inewgg, inpl3d, inpr, iadvance, iforce, rms, clw, cdw, cdpw, cdvw, cxw, cyw, czw, cmxw, cmyw, cmzw, n_clcd, clcd, nblocks_clcd, blocks_clcd, chdw, swetw, fmdotw, cfttotw, cftmomw, cftpw, cftvw, swett, clt, cdt, cxt, cyt, czt, cmxt, cmyt, cmzt, cdpt, cdvt, sx, sy, sz, stot, pav, ptav, tav, ttav, xmav, fmdot, cfxp, cfyp, cfzp, cfdp, cflp, cftp, cfxv, cfyv, cfzv, cfdv, cflv, cftv, cfxmom, cfymom, cfzmom, cfdmom, cflmom, cftmom, cfxtot, cfytot, cfztot, cfdtot, cfltot, cfttot, icsinfo, windex, iindex, nblkpt, windx, iindx, llimit, iitmax, mmcxie, mmceta, ncheck, iifit, mblkpt, iic0, iiorph, iitoss, ifiner, dthetxx, dthetyy, dthetzz, dx, dy, dz, dthetx, dthety, dthetz, lout, xif1, xif2, etf1, etf2, jjmax1, kkmax1, iiint1, iiint2, jimage, kimage, jte, kte, jmm, kmm, nblk1, nblk2, xte, yte, zte, xmi, ymi, zmi, xmie, ymie, zmie, sxie, seta, sxie2, seta2, xie2s, eta2s, temp, x2, y2, z2, x1, y1, z1, factjlo, factjhi, factklo, factkhi, ifrom, geom_miss, period_miss, isav_blk, isav_prd, isav_pat, isav_pat_b, isav_dpat, isav_dpat_b, isav_emb, mblk2nd, ntr, utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae, xorgae, yorgae, zorgae, thtxae, thtyae, thtzae, rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf, kcsi, kcsf, freq, gmass, damp, x0, gf0, nmds, maxaes, aesrfdat, perturb, memuse, bcfiles, islavept, nslave, iskip, jskip, kskip, bmat, stm, stmi, xs, xxn, gforcn, gforcnm, gforcs, nsegdfrm, idfrmseg, iaesurf, maxsegdg, nmaster, aehist, timekeep, xorgae0, yorgae0, zorgae0, icouple, nblelst, iskmax, jskmax, kskmax, ue, irdrea);
}

} // namespace mgblk_ns
