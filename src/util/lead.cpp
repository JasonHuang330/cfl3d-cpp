// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "lead.h"
#include "runtime/common_blocks.h"

namespace lead_ns {

void lead(int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, int& maxbl) {
//
//***********************************************************************
//     Purpose:  Install in common blocks the attributes of a block
//     given a block number.
//***********************************************************************
//
    // Aliases for COMMON /ginfo/
    int32_t& jdim    = cmn_ginfo.jdim;
    int32_t& kdim    = cmn_ginfo.kdim;
    int32_t& idim    = cmn_ginfo.idim;
    int32_t& jj2     = cmn_ginfo.jj2;
    int32_t& kk2     = cmn_ginfo.kk2;
    int32_t& ii2     = cmn_ginfo.ii2;
    int32_t& nblc    = cmn_ginfo.nblc;
    int32_t& js      = cmn_ginfo.js;
    int32_t& ks      = cmn_ginfo.ks;
    int32_t& is      = cmn_ginfo.is;
    int32_t& je      = cmn_ginfo.je;
    int32_t& ke      = cmn_ginfo.ke;
    int32_t& ie      = cmn_ginfo.ie;
    int32_t& lq      = cmn_ginfo.lq;
    int32_t& lqj0    = cmn_ginfo.lqj0;
    int32_t& lqk0    = cmn_ginfo.lqk0;
    int32_t& lqi0    = cmn_ginfo.lqi0;
    int32_t& lsj     = cmn_ginfo.lsj;
    int32_t& lsk     = cmn_ginfo.lsk;
    int32_t& lsi     = cmn_ginfo.lsi;
    int32_t& lvol    = cmn_ginfo.lvol;
    int32_t& ldtj    = cmn_ginfo.ldtj;
    int32_t& lx      = cmn_ginfo.lx;
    int32_t& ly      = cmn_ginfo.ly;
    int32_t& lz      = cmn_ginfo.lz;
    int32_t& lvis    = cmn_ginfo.lvis;
    int32_t& lsnk0   = cmn_ginfo.lsnk0;
    int32_t& lsni0   = cmn_ginfo.lsni0;
    int32_t& lq1     = cmn_ginfo.lq1;
    int32_t& lqr     = cmn_ginfo.lqr;
    int32_t& lblk    = cmn_ginfo.lblk;
    int32_t& lxib    = cmn_ginfo.lxib;
    int32_t& lsig    = cmn_ginfo.lsig;
    int32_t& lsqtq   = cmn_ginfo.lsqtq;
    int32_t& lg      = cmn_ginfo.lg;
    int32_t& ltj0    = cmn_ginfo.ltj0;
    int32_t& ltk0    = cmn_ginfo.ltk0;
    int32_t& lti0    = cmn_ginfo.lti0;
    int32_t& lxkb    = cmn_ginfo.lxkb;
    int32_t& lnbl    = cmn_ginfo.lnbl;
    int32_t& lvj0    = cmn_ginfo.lvj0;
    int32_t& lvk0    = cmn_ginfo.lvk0;
    int32_t& lvi0    = cmn_ginfo.lvi0;
    int32_t& lbcj    = cmn_ginfo.lbcj;
    int32_t& lbck    = cmn_ginfo.lbck;
    int32_t& lbci    = cmn_ginfo.lbci;
    int32_t& lqc0    = cmn_ginfo.lqc0;
    int32_t& ldqc0   = cmn_ginfo.ldqc0;
    int32_t& lxtbi   = cmn_ginfo.lxtbi;
    int32_t& lxtbj   = cmn_ginfo.lxtbj;
    int32_t& lxtbk   = cmn_ginfo.lxtbk;
    int32_t& latbi   = cmn_ginfo.latbi;
    int32_t& latbj   = cmn_ginfo.latbj;
    int32_t& latbk   = cmn_ginfo.latbk;
    int32_t& lbcdj   = cmn_ginfo.lbcdj;
    int32_t& lbcdk   = cmn_ginfo.lbcdk;
    int32_t& lbcdi   = cmn_ginfo.lbcdi;
    int32_t& lxib2   = cmn_ginfo.lxib2;
    int32_t& lux     = cmn_ginfo.lux;
    int32_t& lcmuv   = cmn_ginfo.lcmuv;
    int32_t& lvolj0  = cmn_ginfo.lvolj0;
    int32_t& lvolk0  = cmn_ginfo.lvolk0;
    int32_t& lvoli0  = cmn_ginfo.lvoli0;
    int32_t& lxmdj   = cmn_ginfo.lxmdj;
    int32_t& lxmdk   = cmn_ginfo.lxmdk;
    int32_t& lxmdi   = cmn_ginfo.lxmdi;
    int32_t& lvelg   = cmn_ginfo.lvelg;
    int32_t& ldeltj  = cmn_ginfo.ldeltj;
    int32_t& ldeltk  = cmn_ginfo.ldeltk;
    int32_t& ldelti  = cmn_ginfo.ldelti;
    int32_t& lxnm2   = cmn_ginfo.lxnm2;
    int32_t& lynm2   = cmn_ginfo.lynm2;
    int32_t& lznm2   = cmn_ginfo.lznm2;
    int32_t& lxnm1   = cmn_ginfo.lxnm1;
    int32_t& lynm1   = cmn_ginfo.lynm1;
    int32_t& lznm1   = cmn_ginfo.lznm1;
    int32_t& lqavg   = cmn_ginfo.lqavg;

    // Aliases for COMMON /ginfo2/
    int32_t& lq2avg  = cmn_ginfo2.lq2avg;

    // Aliases for COMMON /lam/
    int32_t& ilamlo  = cmn_lam.ilamlo;
    int32_t& ilamhi  = cmn_lam.ilamhi;
    int32_t& jlamlo  = cmn_lam.jlamlo;
    int32_t& jlamhi  = cmn_lam.jlamhi;
    int32_t& klamlo  = cmn_lam.klamlo;
    int32_t& klamhi  = cmn_lam.klamhi;

    // Aliases for COMMON /reyue/
    int32_t* ivisc   = cmn_reyue.ivisc;  // 0-based

    // Aliases for COMMON /degshf/
    int32_t* ideg    = cmn_degshf.ideg;  // 0-based

    // Aliases for COMMON /wallfun/
    int32_t* iwf     = cmn_wallfun.iwf;  // 0-based

    // Aliases for COMMON /info/
    int32_t* idiag   = cmn_info.idiag;   // 0-based
    int32_t* iflim   = cmn_info.iflim;   // 0-based

    // Aliases for COMMON /fvfds/
    float*   rkap0   = cmn_fvfds.rkap0;  // 0-based
    int32_t* ifds    = cmn_fvfds.ifds;   // 0-based

    // Load from lw2
    jdim      = lw2(1,  nbl);
    kdim      = lw2(2,  nbl);
    idim      = lw2(3,  nbl);
    nblc      = lw2(4,  nbl);
    jj2       = lw2(5,  nbl);
    kk2       = lw2(6,  nbl);
    ii2       = lw2(7,  nbl);
    ilamlo    = lw2(8,  nbl);
    ilamhi    = lw2(9,  nbl);
    jlamlo    = lw2(10, nbl);
    jlamhi    = lw2(11, nbl);
    klamlo    = lw2(12, nbl);
    klamhi    = lw2(13, nbl);
    ivisc[0]  = lw2(14, nbl);
    ivisc[1]  = lw2(15, nbl);
    ivisc[2]  = lw2(16, nbl);
    ideg[0]   = lw2(17, nbl);
    ideg[1]   = lw2(18, nbl);
    ideg[2]   = lw2(19, nbl);
    iwf[0]    = lw2(20, nbl);
    iwf[1]    = lw2(21, nbl);
    iwf[2]    = lw2(22, nbl);
    idiag[0]  = lw2(23, nbl);
    idiag[1]  = lw2(24, nbl);
    idiag[2]  = lw2(25, nbl);
    iflim[0]  = lw2(26, nbl);
    iflim[1]  = lw2(27, nbl);
    iflim[2]  = lw2(28, nbl);
    ifds[0]   = lw2(29, nbl);
    ifds[1]   = lw2(30, nbl);
    ifds[2]   = lw2(31, nbl);
    // factor dividing rkap0g must be consistent with subroutine pointers!
    rkap0[0]  = (float)(lw2(32, nbl) / 1.e6);
    rkap0[1]  = (float)(lw2(33, nbl) / 1.e6);
    rkap0[2]  = (float)(lw2(34, nbl) / 1.e6);
    js        = lw2(35, nbl);
    ks        = lw2(36, nbl);
    is        = lw2(37, nbl);
    je        = lw2(38, nbl);
    ke        = lw2(39, nbl);
    ie        = lw2(40, nbl);

    // Load from lw
    lq     = lw(1,  nbl);
    lqj0   = lw(2,  nbl);
    lqk0   = lw(3,  nbl);
    lqi0   = lw(4,  nbl);
    lsj    = lw(5,  nbl);
    lsk    = lw(6,  nbl);
    lsi    = lw(7,  nbl);
    lvol   = lw(8,  nbl);
    ldtj   = lw(9,  nbl);
    lx     = lw(10, nbl);
    ly     = lw(11, nbl);
    lz     = lw(12, nbl);
    lvis   = lw(13, nbl);
    lsnk0  = lw(14, nbl);
    lsni0  = lw(15, nbl);
    lq1    = lw(16, nbl);
    lqr    = lw(17, nbl);
    lblk   = lw(18, nbl);
    lxib   = lw(19, nbl);
    lsig   = lw(20, nbl);
    lsqtq  = lw(21, nbl);
    lg     = lw(22, nbl);
    ltj0   = lw(23, nbl);
    ltk0   = lw(24, nbl);
    lti0   = lw(25, nbl);
    lxkb   = lw(26, nbl);
    lnbl   = lw(27, nbl);
    lvj0   = lw(28, nbl);
    lvk0   = lw(29, nbl);
    lvi0   = lw(30, nbl);
    lbcj   = lw(31, nbl);
    lbck   = lw(32, nbl);
    lbci   = lw(33, nbl);
    lqc0   = lw(34, nbl);
    ldqc0  = lw(35, nbl);
    lxtbj  = lw(36, nbl);
    lxtbk  = lw(37, nbl);
    lxtbi  = lw(38, nbl);
    latbj  = lw(39, nbl);
    latbk  = lw(40, nbl);
    latbi  = lw(41, nbl);
    lbcdj  = lw(42, nbl);
    lbcdk  = lw(43, nbl);
    lbcdi  = lw(44, nbl);
    lxib2  = lw(45, nbl);
    lqavg  = lw(46, nbl);
    lux    = lw(47, nbl);
    lcmuv  = lw(48, nbl);
    lvolj0 = lw(49, nbl);
    lvolk0 = lw(50, nbl);
    lvoli0 = lw(51, nbl);
    lxmdj  = lw(52, nbl);
    lxmdk  = lw(53, nbl);
    lxmdi  = lw(54, nbl);
    lvelg  = lw(55, nbl);
    lxnm2  = lw(56, nbl);
    lynm2  = lw(57, nbl);
    lznm2  = lw(58, nbl);
    ldeltj = lw(59, nbl);
    ldeltk = lw(60, nbl);
    ldelti = lw(61, nbl);
    lxnm1  = lw(62, nbl);
    lynm1  = lw(63, nbl);
    lznm1  = lw(64, nbl);
    lq2avg = lw(65, nbl);

    return;
}

void pointers(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, int& maxl, int& lembed, int& nstart, int& nwork, int& mwork, int& maxbl, int& maxgr, int& maxseg, FortranArray3DRef<int> lwdat, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, int& nblock, int& myhost, int& myid, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray1DRef<int> nblfine, FortranArray1DRef<int> ilamlog, FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> klamlog, FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamhig, FortranArray2DRef<int> idegg, FortranArray2DRef<int> iwfg, FortranArray2DRef<int> idiagg, FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg, FortranArray2DRef<double> rkap0g, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray1DRef<int> memblock, int& icall, int& nmds, int& maxaes, int& mpihost) {
    // pointers subroutine: not present in this Fortran source file.
    // Body is empty (no Fortran source to translate).
    return;
}

} // namespace lead_ns
