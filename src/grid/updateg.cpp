// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "updateg.h"
#include "lead.h"
#include "trans.h"
#include "rotate.h"
#include "xtbatb.h"
#include "metric.h"
#include "tmetric.h"
#include "collx.h"
#include "collxt.h"
#include "collxtb.h"
#include "termn8.h"
#include "writ_buf.h"
#include "transmc.h"
#include "rotatmc.h"
#include "update.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <algorithm>

namespace updateg_ns {


// Forward wrappers for procedures declared in updateg_ns but implemented elsewhere

void metric(int& jdim, int& kdim, int& idim,
            FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
            FortranArray3DRef<double> sj, FortranArray3DRef<double> sk, FortranArray3DRef<double> si,
            FortranArray2DRef<double> t, FortranArray3DRef<double> t1,
            int& nbl, int& iflag, int& icnt,
            FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
            FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim,
            FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
            int& maxbl, int& maxseg,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
            int& myid, FortranArray1DRef<int> mblk2nd)
{
    metric_ns::metric(jdim, kdim, idim, x, y, z, sj, sk, si, t, t1, nbl, iflag, icnt,
                      nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
                      ibcinfo, jbcinfo, kbcinfo, maxbl, maxseg,
                      nou, bou, nbuf, ibufdim, myid, mblk2nd);
}

void rotate(int& jdim, int& kdim, int& idim,
            FortranArray2DRef<double> t, FortranArray3DRef<double> tti,
            FortranArray3DRef<double> ttj, FortranArray3DRef<double> ttk,
            FortranArray2DRef<double> t1,
            FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
            int& nbl, int& irot, double& rfreqr,
            double& omegx, double& omegy, double& omegz,
            double& xorg, double& yorg, double& zorg,
            double& thetax, double& thetay, double& thetaz,
            double& thxold, double& thyold, double& thzold,
            int& iupdat, double& time2,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    rotate_ns::rotate(jdim, kdim, idim, t, tti, ttj, ttk, t1, x, y, z,
                      nbl, irot, rfreqr, omegx, omegy, omegz,
                      xorg, yorg, zorg, thetax, thetay, thetaz,
                      thxold, thyold, thzold, iupdat, time2,
                      nou, bou, nbuf, ibufdim);
}

void tmetric(int& jdim, int& kdim, int& idim,
             FortranArray3DRef<double> sj, FortranArray3DRef<double> sk, FortranArray3DRef<double> si,
             FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
             FortranArray2DRef<double> t, FortranArray2DRef<double> t1, FortranArray2DRef<double> t2,
             FortranArray3DRef<double> t3, int& nbl)
{
    tmetric_ns::tmetric(jdim, kdim, idim, sj, sk, si, x, y, z, t, t1, t2, t3, nbl);
}

void trans(int& jdim, int& kdim, int& idim,
           FortranArray2DRef<double> t, FortranArray3DRef<double> tti,
           FortranArray3DRef<double> ttj, FortranArray3DRef<double> ttk,
           FortranArray1DRef<double> x, FortranArray1DRef<double> y, FortranArray1DRef<double> z,
           int& itran, double& rfreqt,
           double& utran, double& vtran, double& wtran,
           double& xorg, double& yorg, double& zorg,
           double& xold, double& yold, double& zold,
           double& xorg0, double& yorg0, double& zorg0,
           int& iupdat, double& time2)
{
    trans_ns::trans(jdim, kdim, idim, t, tti, ttj, ttk, x, y, z,
                    itran, rfreqt, utran, vtran, wtran,
                    xorg, yorg, zorg, xold, yold, zold,
                    xorg0, yorg0, zorg0, iupdat, time2);
}

void update(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray2DRef<double> qj0, FortranArray2DRef<double> qk0, FortranArray2DRef<double> qi0,
            FortranArray2DRef<double> sj, FortranArray2DRef<double> sk, FortranArray2DRef<double> si,
            FortranArray1DRef<double> vol, FortranArray1DRef<double> dtj,
            FortranArray3DRef<double> vist3d, FortranArray3DRef<double> blank,
            FortranArray1DRef<double> x, FortranArray1DRef<double> y, FortranArray1DRef<double> z,
            FortranArray4DRef<double> res, FortranArray1DRef<double> wk0,
            FortranArray3DRef<double> vmuk, FortranArray3DRef<double> vmuj, FortranArray3DRef<double> vmui,
            FortranArray1DRef<double> wk, int& nwork, int& nbl, int& iover,
            FortranArray4DRef<double> vk0,
            FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
            int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl,
            FortranArray3DRef<double> volk0, FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tk0, FortranArray3DRef<double> cmuv,
            FortranArray1DRef<int> iadvance, int& nummem, FortranArray4DRef<double> ux)
{
    update_ns::update(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                      vol, dtj, vist3d, blank, x, y, z, res, wk0,
                      vmuk, vmuj, vmui, wk, nwork, nbl, iover, vk0,
                      bcj, bck, bci, nou, bou, nbuf, ibufdim,
                      myid, mblk2nd, maxbl, volk0, tursav, tk0, cmuv,
                      iadvance, nummem, ux);
}



void updateg(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
             FortranArray1DRef<double> w, int& mgwk,
             FortranArray1DRef<double> wk, int& nwork,
             int& iupdat, int& iseqr, int& maxbl, int& maxgr, int& maxseg,
             FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
             FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim,
             FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
             int& nblock, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg,
             FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans,
             FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz,
             FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig,
             FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz,
             FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr,
             FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0,
             FortranArray1DRef<double> time2,
             FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl,
             FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm,
             FortranArray1DRef<int> ncgg, FortranArray1DRef<int> iadvance,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
             int& myid, int& myhost, int& mycomm,
             FortranArray1DRef<int> mblk2nd, int& irigb, int& irbtrim, int& nt)
{
    // COMMON block aliases
    int& jdim    = cmn_ginfo.jdim;
    int& kdim    = cmn_ginfo.kdim;
    int& idim    = cmn_ginfo.idim;
    int& jj2     = cmn_ginfo.jj2;
    int& kk2     = cmn_ginfo.kk2;
    int& ii2     = cmn_ginfo.ii2;
    int& lq      = cmn_ginfo.lq;
    int& lqj0    = cmn_ginfo.lqj0;
    int& lqk0    = cmn_ginfo.lqk0;
    int& lqi0    = cmn_ginfo.lqi0;
    int& lsj     = cmn_ginfo.lsj;
    int& lsk     = cmn_ginfo.lsk;
    int& lsi     = cmn_ginfo.lsi;
    int& lvol    = cmn_ginfo.lvol;
    int& ldtj    = cmn_ginfo.ldtj;
    int& lx      = cmn_ginfo.lx;
    int& ly      = cmn_ginfo.ly;
    int& lz      = cmn_ginfo.lz;
    int& lvis    = cmn_ginfo.lvis;
    int& lsnk0   = cmn_ginfo.lsnk0;
    int& lsni0   = cmn_ginfo.lsni0;
    int& lq1     = cmn_ginfo.lq1;
    int& lqr     = cmn_ginfo.lqr;
    int& lblk    = cmn_ginfo.lblk;
    int& lxib    = cmn_ginfo.lxib;
    int& lsig    = cmn_ginfo.lsig;
    int& lsqtq   = cmn_ginfo.lsqtq;
    int& lg      = cmn_ginfo.lg;
    int& ltj0    = cmn_ginfo.ltj0;
    int& ltk0    = cmn_ginfo.ltk0;
    int& lti0    = cmn_ginfo.lti0;
    int& lxkb    = cmn_ginfo.lxkb;
    int& lnbl    = cmn_ginfo.lnbl;
    int& lvj0    = cmn_ginfo.lvj0;
    int& lvk0    = cmn_ginfo.lvk0;
    int& lvi0    = cmn_ginfo.lvi0;
    int& lbcj    = cmn_ginfo.lbcj;
    int& lbck    = cmn_ginfo.lbck;
    int& lbci    = cmn_ginfo.lbci;
    int& lqc0    = cmn_ginfo.lqc0;
    int& ldqc0   = cmn_ginfo.ldqc0;
    int& lxtbi   = cmn_ginfo.lxtbi;
    int& lxtbj   = cmn_ginfo.lxtbj;
    int& lxtbk   = cmn_ginfo.lxtbk;
    int& latbi   = cmn_ginfo.latbi;
    int& latbj   = cmn_ginfo.latbj;
    int& latbk   = cmn_ginfo.latbk;
    int& lbcdj   = cmn_ginfo.lbcdj;
    int& lbcdk   = cmn_ginfo.lbcdk;
    int& lbcdi   = cmn_ginfo.lbcdi;
    int& lxib2   = cmn_ginfo.lxib2;
    int& lux     = cmn_ginfo.lux;
    int& lcmuv   = cmn_ginfo.lcmuv;
    int& lvolj0  = cmn_ginfo.lvolj0;
    int& lvolk0  = cmn_ginfo.lvolk0;
    int& lvoli0  = cmn_ginfo.lvoli0;
    int& lxmdj   = cmn_ginfo.lxmdj;
    int& lxmdk   = cmn_ginfo.lxmdk;
    int& lxmdi   = cmn_ginfo.lxmdi;
    int& lvelg   = cmn_ginfo.lvelg;
    int& ldeltj  = cmn_ginfo.ldeltj;
    int& ldeltk  = cmn_ginfo.ldeltk;
    int& ldelti  = cmn_ginfo.ldelti;
    int& lxnm2   = cmn_ginfo.lxnm2;
    int& lynm2   = cmn_ginfo.lynm2;
    int& lznm2   = cmn_ginfo.lznm2;
    int& lxnm1   = cmn_ginfo.lxnm1;
    int& lynm1   = cmn_ginfo.lynm1;
    int& lznm1   = cmn_ginfo.lznm1;
    int& lqavg   = cmn_ginfo.lqavg;

    int& mseq    = cmn_info.mseq;
    int& mgflag  = cmn_info.mgflag;
    float& dt    = cmn_info.dt;

    int& lglobal = cmn_mgrd.lglobal;
    int& isklton = cmn_sklton.isklton;

    float& xmc   = cmn_fsum.xmc;
    float& ymc   = cmn_fsum.ymc;
    float& zmc   = cmn_fsum.zmc;

    int32_t& itransmc  = cmn_motionmc.itransmc;
    int32_t& irotatmc  = cmn_motionmc.irotatmc;
    float& rfreqtmc    = cmn_motionmc.rfreqtmc;
    float& rfreqrmc    = cmn_motionmc.rfreqrmc;
    float& utransmc    = cmn_motionmc.utransmc;
    float& vtransmc    = cmn_motionmc.vtransmc;
    float& wtransmc    = cmn_motionmc.wtransmc;
    float& omegaxmc    = cmn_motionmc.omegaxmc;
    float& omegaymc    = cmn_motionmc.omegaymc;
    float& omegazmc    = cmn_motionmc.omegazmc;
    float& xorigmc     = cmn_motionmc.xorigmc;
    float& yorigmc     = cmn_motionmc.yorigmc;
    float& zorigmc     = cmn_motionmc.zorigmc;
    float& xorig0mc    = cmn_motionmc.xorig0mc;
    float& yorig0mc    = cmn_motionmc.yorig0mc;
    float& zorig0mc    = cmn_motionmc.zorig0mc;
    float& thetaxmc    = cmn_motionmc.thetaxmc;
    float& thetaymc    = cmn_motionmc.thetaymc;
    float& thetazmc    = cmn_motionmc.thetazmc;
    float& time2mc     = cmn_motionmc.time2mc;

    float& epstr       = cmn_trim.epstr;
    float& zrg1        = cmn_trim.zrg1;

    // isklton check at entry
    if (isklton == 1) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, " ");
    }

    int icnt = 0;



    for (int nbl = 1; nbl <= nblock; nbl++) {
        if (myid == mblk2nd(nbl) && iadvance(nbl) >= 0 &&
            (levelg(nbl) >= lglobal &&
             levelg(nbl) <= cmn_info.levelt[iseqr - 1])) {

            lead_ns::lead(nbl, lw, lw2, maxbl);
            int iuns1 = std::max(irotat(nbl), itrans(nbl));
            int iuns  = std::max(iuns1, std::max(irigb, irbtrim));

            if (iuns > 0) {
                // temporary storage locations:
                // lvel  = start of grid point velocity array
                // lacci = start of i-boundary point acceleration array
                // laccj = start of j-boundary point acceleration array
                // lacck = start of k-boundary point acceleration array
                // lt1wk = start of work array for subroutines rotate/metric
                // lt2wk = start of work array for subroutine metric
                // lt3wk = start of work array for subroutine metric
                int lvel  = 1;
                int lacci = jdim * kdim * idim * 3 + lvel;
                int laccj = jdim * kdim * 3 * 2 + lacci;
                int lacck = kdim * idim * 3 * 2 + laccj;
                int lt1wk = jdim * idim * 3 * 2 + lacck;
                int lt2wk = jdim * kdim * idim * 3 + lt1wk;
                int lt3wk = jdim * kdim * 6 + lt2wk;

                // zero out velocity/acceleration work arrays
                int mdim = lt1wk - 1;
                for (int izz = 1; izz <= mdim; izz++) {
                    wk(lvel + izz - 1) = 0.0e0;
                }

                // translation/rotation corresponding to rigid-body modes
                if (irigb > 0 && nt > 0) {

                    if (isklton == 1) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                                      " translating rb-block%4d to new position", nbl);
                    }

                    if (nwork < mdim) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                                      " not enough work space for subroutine trans");
                        int m1 = -1;
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    }

                    double xold = xorig(nbl);
                    double yold = yorig(nbl);
                    double zold = zorig(nbl);

                    // This formulation assumes small pitch angles.
                    zorig(nbl)  = (double)zrg1;
                    utrans(nbl) = 0.0;
                    vtrans(nbl) = 0.0;
                    wtrans(nbl) = (zorig(nbl) - zold) / (double)dt;

                    {
                        FortranArray2DRef<double> t_vel(&wk(lvel), jdim * kdim * idim, 3);
                        FortranArray3DRef<double> tti_acci(&wk(lacci), jdim, kdim, idim * 3 * 2);
                        FortranArray3DRef<double> ttj_accj(&wk(laccj), kdim, idim, 3 * 2);
                        FortranArray3DRef<double> ttk_acck(&wk(lacck), jdim, idim, 3 * 2);
                        FortranArray1DRef<double> wx(&w(lx), jdim * kdim * idim);
                        FortranArray1DRef<double> wy(&w(ly), jdim * kdim * idim);
                        FortranArray1DRef<double> wz(&w(lz), jdim * kdim * idim);
                        int itran99 = 99;
                        trans_ns::trans(jdim, kdim, idim,
                                        t_vel, tti_acci, ttj_accj, ttk_acck,
                                        wx, wy, wz,
                                        itran99, rfreqt(nbl),
                                        utrans(nbl), vtrans(nbl), wtrans(nbl),
                                        xorig(nbl), yorig(nbl), zorig(nbl),
                                        xold, yold, zold,
                                        xorig0(nbl), yorig0(nbl), zorig0(nbl),
                                        iupdat, time2(nbl));
                    }

                    // rotation corresponding to rigid-body modes
                    if (isklton == 1) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                                      " rotating rb-block   %4d to new position", nbl);
                    }

                    {
                        int nroom = nwork - lt1wk;
                        mdim = jdim * kdim * idim * 3;
                        if (nroom < mdim) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                          " not enough work space for subroutine rotate");
                            int m1 = -1;
                            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                        }
                    }

                    {
                        FortranArray2DRef<double> t_vel(&wk(lvel), jdim * kdim * idim, 3);
                        FortranArray3DRef<double> tti_acci(&wk(lacci), jdim, kdim, idim * 3 * 2);
                        FortranArray3DRef<double> ttj_accj(&wk(laccj), kdim, idim, 3 * 2);
                        FortranArray3DRef<double> ttk_acck(&wk(lacck), jdim, idim, 3 * 2);
                        FortranArray2DRef<double> t1_lt1wk(&wk(lt1wk), jdim * kdim * idim, 3);
                        FortranArray3DRef<double> wx3(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> wy3(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> wz3(&w(lz), jdim, kdim, idim);
                        int irot99 = 99;
                        rotate_ns::rotate(jdim, kdim, idim,
                                          t_vel, tti_acci, ttj_accj, ttk_acck,
                                          t1_lt1wk, wx3, wy3, wz3,
                                          nbl, irot99, rfreqr(nbl),
                                          omegax(nbl), omegay(nbl), omegaz(nbl),
                                          xorig(nbl), yorig(nbl), zorig(nbl),
                                          thetax(nbl), thetay(nbl), thetaz(nbl),
                                          thetaxl(nbl), thetayl(nbl), thetazl(nbl),
                                          iupdat, time2(nbl),
                                          nou, bou, nbuf, ibufdim);
                    }
                } // end irigb > 0 && nt > 0



                // trim step for rigid-body modes
                if (irbtrim > 0 && nt > 0) {

                    if (isklton == 1) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                                      " rotating rb-block  %4d to trim", nbl);
                    }

                    {
                        int nroom = nwork - lt1wk;
                        int mdim2 = jdim * kdim * idim * 3;
                        if (nroom < mdim2) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                          " not enough work space for subroutine rotate");
                            int m1 = -1;
                            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                        }
                    }

                    omegay(nbl) = (double)epstr * (thetay(nbl) - thetayl(nbl));

                    {
                        FortranArray2DRef<double> t_vel(&wk(lvel), jdim * kdim * idim, 3);
                        FortranArray3DRef<double> tti_acci(&wk(lacci), jdim, kdim, idim * 3 * 2);
                        FortranArray3DRef<double> ttj_accj(&wk(laccj), kdim, idim, 3 * 2);
                        FortranArray3DRef<double> ttk_acck(&wk(lacck), jdim, idim, 3 * 2);
                        FortranArray2DRef<double> t1_lt1wk(&wk(lt1wk), jdim * kdim * idim, 3);
                        FortranArray3DRef<double> wx3(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> wy3(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> wz3(&w(lz), jdim, kdim, idim);
                        int irot99 = 99;
                        rotate_ns::rotate(jdim, kdim, idim,
                                          t_vel, tti_acci, ttj_accj, ttk_acck,
                                          t1_lt1wk, wx3, wy3, wz3,
                                          nbl, irot99, rfreqr(nbl),
                                          omegax(nbl), omegay(nbl), omegaz(nbl),
                                          xorig(nbl), yorig(nbl), zorig(nbl),
                                          thetax(nbl), thetay(nbl), thetaz(nbl),
                                          thetaxl(nbl), thetayl(nbl), thetazl(nbl),
                                          iupdat, time2(nbl),
                                          nou, bou, nbuf, ibufdim);
                    }
                } // end irbtrim > 0 && nt > 0

                // translation
                if (itrans(nbl) > 0 && irbtrim == 0 &&
                    irigb == 0 && nt > 0) {

                    if (isklton == 1) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                                      " translating block%4d to new position", nbl);
                    }

                    if (nwork < mdim) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                                      " not enough work space for subroutine trans");
                        int m1 = -1;
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    }

                    double xold = xorig(nbl);
                    double yold = yorig(nbl);
                    double zold = zorig(nbl);

                    {
                        FortranArray2DRef<double> t_vel(&wk(lvel), jdim * kdim * idim, 3);
                        FortranArray3DRef<double> tti_acci(&wk(lacci), jdim, kdim, idim * 3 * 2);
                        FortranArray3DRef<double> ttj_accj(&wk(laccj), kdim, idim, 3 * 2);
                        FortranArray3DRef<double> ttk_acck(&wk(lacck), jdim, idim, 3 * 2);
                        FortranArray1DRef<double> wx(&w(lx), jdim * kdim * idim);
                        FortranArray1DRef<double> wy(&w(ly), jdim * kdim * idim);
                        FortranArray1DRef<double> wz(&w(lz), jdim * kdim * idim);
                        trans_ns::trans(jdim, kdim, idim,
                                        t_vel, tti_acci, ttj_accj, ttk_acck,
                                        wx, wy, wz,
                                        itrans(nbl), rfreqt(nbl),
                                        utrans(nbl), vtrans(nbl), wtrans(nbl),
                                        xorig(nbl), yorig(nbl), zorig(nbl),
                                        xold, yold, zold,
                                        xorig0(nbl), yorig0(nbl), zorig0(nbl),
                                        iupdat, time2(nbl));
                    }
                } // end itrans > 0

                // rotation
                if (irotat(nbl) > 0 && irbtrim == 0 &&
                    irigb == 0 && nt > 0) {

                    if (isklton == 1) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                                      " rotating block   %4d to new position", nbl);
                    }

                    {
                        int nroom = nwork - lt1wk;
                        int mdim2 = jdim * kdim * idim * 3;
                        if (nroom < mdim2) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                          " not enough work space for subroutine rotate");
                            int m1 = -1;
                            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                        }
                    }

                    thetaxl(nbl) = thetax(nbl);
                    thetayl(nbl) = thetay(nbl);
                    thetazl(nbl) = thetaz(nbl);

                    {
                        FortranArray2DRef<double> t_vel(&wk(lvel), jdim * kdim * idim, 3);
                        FortranArray3DRef<double> tti_acci(&wk(lacci), jdim, kdim, idim * 3 * 2);
                        FortranArray3DRef<double> ttj_accj(&wk(laccj), kdim, idim, 3 * 2);
                        FortranArray3DRef<double> ttk_acck(&wk(lacck), jdim, idim, 3 * 2);
                        FortranArray2DRef<double> t1_lt1wk(&wk(lt1wk), jdim * kdim * idim, 3);
                        FortranArray3DRef<double> wx3(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> wy3(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> wz3(&w(lz), jdim, kdim, idim);
                        rotate_ns::rotate(jdim, kdim, idim,
                                          t_vel, tti_acci, ttj_accj, ttk_acck,
                                          t1_lt1wk, wx3, wy3, wz3,
                                          nbl, irotat(nbl), rfreqr(nbl),
                                          omegax(nbl), omegay(nbl), omegaz(nbl),
                                          xorig(nbl), yorig(nbl), zorig(nbl),
                                          thetax(nbl), thetay(nbl), thetaz(nbl),
                                          thetaxl(nbl), thetayl(nbl), thetazl(nbl),
                                          iupdat, time2(nbl),
                                          nou, bou, nbuf, ibufdim);
                    }
                } // end irotat > 0



                // if the current block will also undergo deformation, defer
                // updating metrics until the changes due to deformation are added.
                if (idefrm(nbl) == 0) {

                    // calculate face-average values of velocity and acceleration
                    // on block boundaries
                    {
                        FortranArray4DRef<double> wxtbj(&w(lxtbj), kdim, idim - 1, 3, 2);
                        FortranArray4DRef<double> wxtbk(&w(lxtbk), jdim, idim - 1, 3, 2);
                        FortranArray4DRef<double> wxtbi(&w(lxtbi), jdim, kdim, 3, 2);
                        FortranArray4DRef<double> watbj(&w(latbj), kdim, idim - 1, 3, 2);
                        FortranArray4DRef<double> watbk(&w(latbk), jdim, idim - 1, 3, 2);
                        FortranArray4DRef<double> watbi(&w(latbi), jdim, kdim, 3, 2);
                        FortranArray2DRef<double> t_vel(&wk(lvel), jdim * kdim * idim, 3);
                        FortranArray3DRef<double> tti_acci(&wk(lacci), jdim, kdim, idim * 3 * 2);
                        FortranArray3DRef<double> ttj_accj(&wk(laccj), kdim, idim, 3 * 2);
                        FortranArray3DRef<double> ttk_acck(&wk(lacck), jdim, idim, 3 * 2);
                        xtbatb_ns::xtbatb(jdim, kdim, idim,
                                          wxtbj, wxtbk, wxtbi,
                                          watbj, watbk, watbi,
                                          t_vel, tti_acci, ttj_accj, ttk_acck);
                    }

                    {
                        int nroom = nwork - lt3wk;
                        int mdim2 = jdim * kdim * idim * 5;
                        if (nroom < mdim2) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                          " not enough work space for metric subroutines");
                            int m1 = -1;
                            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                        }
                    }

                    // calculate spatial metrics for updated grid
                    int iflag = -1;
                    {
                        FortranArray3DRef<double> wx3(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> wy3(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> wz3(&w(lz), jdim, kdim, idim);
                        FortranArray3DRef<double> wsj(&w(lsj), jdim, kdim, idim);
                        FortranArray3DRef<double> wsk(&w(lsk), jdim, kdim, idim);
                        FortranArray3DRef<double> wsi(&w(lsi), jdim, kdim, idim);
                        FortranArray2DRef<double> wt2(&wk(lt2wk), jdim * kdim, idim * 3);
                        FortranArray3DRef<double> wt3(&wk(lt3wk), jdim * kdim, idim, 5);
                        metric_ns::metric(jdim, kdim, idim,
                                          wx3, wy3, wz3,
                                          wsj, wsk, wsi,
                                          wt2, wt3,
                                          nbl, iflag, icnt,
                                          nbci0, nbcj0, nbck0,
                                          nbcidim, nbcjdim, nbckdim,
                                          ibcinfo, jbcinfo, kbcinfo,
                                          maxbl, maxseg,
                                          nou, bou, nbuf, ibufdim,
                                          myid, mblk2nd);
                    }

                    // calculate temporal metrics for updated grid
                    {
                        FortranArray3DRef<double> wsj(&w(lsj), jdim, kdim, idim);
                        FortranArray3DRef<double> wsk(&w(lsk), jdim, kdim, idim);
                        FortranArray3DRef<double> wsi(&w(lsi), jdim, kdim, idim);
                        FortranArray3DRef<double> wx3(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> wy3(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> wz3(&w(lz), jdim, kdim, idim);
                        FortranArray2DRef<double> t_vel(&wk(lvel), jdim * kdim * idim, 3);
                        FortranArray2DRef<double> wt1(&wk(lt1wk), jdim * kdim * idim, 3);
                        FortranArray2DRef<double> wt2(&wk(lt2wk), jdim * kdim, idim * 3);
                        FortranArray3DRef<double> wt3(&wk(lt3wk), jdim * kdim, idim, 5);
                        tmetric_ns::tmetric(jdim, kdim, idim,
                                            wsj, wsk, wsi,
                                            wx3, wy3, wz3,
                                            t_vel, wt1, wt2, wt3,
                                            nbl);
                    }

                    // coarser levels
                    int ncg = ncgg(igridg(nbl)) - (mseq - iseqr);
                    if (ncg > 0 && mgflag > 0) {
                        int nbll = nbl;
                        for (int m = 1; m <= ncg; m++) {
                            nbll = nbll + 1;
                            int nbllm1 = nbll - 1;
                            time2(nbll)  = time2(nbl);
                            xorig(nbll)  = xorig(nbl);
                            yorig(nbll)  = yorig(nbl);
                            zorig(nbll)  = zorig(nbl);
                            thetax(nbll) = thetax(nbl);
                            thetay(nbll) = thetay(nbl);
                            thetaz(nbll) = thetaz(nbl);
                            int lvolc  = lw(8,  nbll);
                            int lxc    = lw(10, nbll);
                            int lyc    = lw(11, nbll);
                            int lzc    = lw(12, nbll);
                            int lxtbjc = lw(36, nbll);
                            int lxtbkc = lw(37, nbll);
                            int lxtbic = lw(38, nbll);
                            int latbjc = lw(39, nbll);
                            int latbkc = lw(40, nbll);
                            int latbic = lw(41, nbll);
                            int lvelc  = lt1wk;

                            if (isklton == 1) {
                                nou(1) = std::min(nou(1) + 1, ibufdim);
                                std::snprintf(bou(nou(1), 1), 120,
                                              "   creating coarser block%4d of dimensions (I/J/K) :%4d%4d%4d",
                                              nbll, ii2, jj2, kk2);
                            }

                            // collocate xyz
                            {
                                FortranArray3DRef<double> wx3(&w(lx), jdim, kdim, idim);
                                FortranArray3DRef<double> wy3(&w(ly), jdim, kdim, idim);
                                FortranArray3DRef<double> wz3(&w(lz), jdim, kdim, idim);
                                FortranArray3DRef<double> wxc(&w(lxc), jj2, kk2, ii2);
                                FortranArray3DRef<double> wyc(&w(lyc), jj2, kk2, ii2);
                                FortranArray3DRef<double> wzc(&w(lzc), jj2, kk2, ii2);
                                collx_ns::collx(wx3, wy3, wz3, wxc, wyc, wzc,
                                                jdim, kdim, idim, jj2, kk2, ii2);
                            }

                            // collocate grid point velocity
                            {
                                FortranArray4DRef<double> wvel(&wk(lvel), jdim, kdim, idim, 3);
                                FortranArray4DRef<double> wvelc(&wk(lvelc), jj2, kk2, ii2, 3);
                                collxt_ns::collxt(wvel, wvelc,
                                                  jdim, kdim, idim, jj2, kk2, ii2,
                                                  nbllm1, nou, bou, nbuf, ibufdim);
                            }
                            int nv = jj2 * kk2 * ii2 * 3;
                            for (int izz = 1; izz <= nv; izz++) {
                                wk(lvel + izz - 1) = wk(lvelc + izz - 1);
                            }



                            // collocate i0/idim boundary velocity/acceleration
                            {
                                FortranArray4DRef<double> wxtbi(&w(lxtbi), jdim, kdim, 3, 2);
                                FortranArray4DRef<double> wxtbic(&w(lxtbic), jj2, kk2, 3, 2);
                                collxtb_ns::collxtb(wxtbi, wxtbic, jdim, kdim, jj2, kk2, nbllm1);
                            }
                            {
                                FortranArray4DRef<double> watbi(&w(latbi), jdim, kdim, 3, 2);
                                FortranArray4DRef<double> watbic(&w(latbic), jj2, kk2, 3, 2);
                                collxtb_ns::collxtb(watbi, watbic, jdim, kdim, jj2, kk2, nbllm1);
                            }

                            // collocate j0/jdim boundary velocity/acceleration
                            {
                                int kdim_v = kdim;
                                int idim_m1 = idim - 1;
                                int kk2_v   = kk2;
                                int ii2_m1  = ii2 - 1;
                                FortranArray4DRef<double> wxtbj(&w(lxtbj), kdim, idim - 1, 3, 2);
                                FortranArray4DRef<double> wxtbjc(&w(lxtbjc), kk2, ii2 - 1, 3, 2);
                                collxtb_ns::collxtb(wxtbj, wxtbjc, kdim_v, idim_m1, kk2_v, ii2_m1, nbllm1);
                            }
                            {
                                int kdim_v = kdim;
                                int idim_m1 = idim - 1;
                                int kk2_v   = kk2;
                                int ii2_m1  = ii2 - 1;
                                FortranArray4DRef<double> watbj(&w(latbj), kdim, idim - 1, 3, 2);
                                FortranArray4DRef<double> watbjc(&w(latbjc), kk2, ii2 - 1, 3, 2);
                                collxtb_ns::collxtb(watbj, watbjc, kdim_v, idim_m1, kk2_v, ii2_m1, nbllm1);
                            }

                            // collocate k0/kdim boundary velocity/acceleration
                            {
                                int jdim_v  = jdim;
                                int idim_m1 = idim - 1;
                                int jj2_v   = jj2;
                                int ii2_m1  = ii2 - 1;
                                FortranArray4DRef<double> wxtbk(&w(lxtbk), jdim, idim - 1, 3, 2);
                                FortranArray4DRef<double> wxtbkc(&w(lxtbkc), jj2, ii2 - 1, 3, 2);
                                collxtb_ns::collxtb(wxtbk, wxtbkc, jdim_v, idim_m1, jj2_v, ii2_m1, nbllm1);
                            }
                            {
                                int jdim_v  = jdim;
                                int idim_m1 = idim - 1;
                                int jj2_v   = jj2;
                                int ii2_m1  = ii2 - 1;
                                FortranArray4DRef<double> watbk(&w(latbk), jdim, idim - 1, 3, 2);
                                FortranArray4DRef<double> watbkc(&w(latbkc), jj2, ii2 - 1, 3, 2);
                                collxtb_ns::collxtb(watbk, watbkc, jdim_v, idim_m1, jj2_v, ii2_m1, nbllm1);
                            }

                            // calculate spatial metrics for updated coarser grid
                            lead_ns::lead(nbll, lw, lw2, maxbl);

                            lvel  = 1;
                            lt1wk = jdim * kdim * idim * 3 + lvel;
                            lt2wk = jdim * kdim * idim * 3 + lt1wk;
                            lt3wk = jdim * kdim * 6 + lt2wk;

                            iflag = -1;
                            {
                                FortranArray3DRef<double> wx3(&w(lx), jdim, kdim, idim);
                                FortranArray3DRef<double> wy3(&w(ly), jdim, kdim, idim);
                                FortranArray3DRef<double> wz3(&w(lz), jdim, kdim, idim);
                                FortranArray3DRef<double> wsj(&w(lsj), jdim, kdim, idim);
                                FortranArray3DRef<double> wsk(&w(lsk), jdim, kdim, idim);
                                FortranArray3DRef<double> wsi(&w(lsi), jdim, kdim, idim);
                                FortranArray2DRef<double> wt2(&wk(lt2wk), jdim * kdim, idim * 3);
                                FortranArray3DRef<double> wt3(&wk(lt3wk), jdim * kdim, idim, 5);
                                metric_ns::metric(jdim, kdim, idim,
                                                  wx3, wy3, wz3,
                                                  wsj, wsk, wsi,
                                                  wt2, wt3,
                                                  nbll, iflag, icnt,
                                                  nbci0, nbcj0, nbck0,
                                                  nbcidim, nbcjdim, nbckdim,
                                                  ibcinfo, jbcinfo, kbcinfo,
                                                  maxbl, maxseg,
                                                  nou, bou, nbuf, ibufdim,
                                                  myid, mblk2nd);
                            }

                            // calculate temporal metrics for updated coarser grid
                            {
                                FortranArray3DRef<double> wsj(&w(lsj), jdim, kdim, idim);
                                FortranArray3DRef<double> wsk(&w(lsk), jdim, kdim, idim);
                                FortranArray3DRef<double> wsi(&w(lsi), jdim, kdim, idim);
                                FortranArray3DRef<double> wx3(&w(lx), jdim, kdim, idim);
                                FortranArray3DRef<double> wy3(&w(ly), jdim, kdim, idim);
                                FortranArray3DRef<double> wz3(&w(lz), jdim, kdim, idim);
                                FortranArray2DRef<double> t_vel(&wk(lvel), jdim * kdim * idim, 3);
                                FortranArray2DRef<double> wt1(&wk(lt1wk), jdim * kdim * idim, 3);
                                FortranArray2DRef<double> wt2(&wk(lt2wk), jdim * kdim, idim * 3);
                                FortranArray3DRef<double> wt3(&wk(lt3wk), jdim * kdim, idim, 5);
                                tmetric_ns::tmetric(jdim, kdim, idim,
                                                    wsj, wsk, wsi,
                                                    wx3, wy3, wz3,
                                                    t_vel, wt1, wt2, wt3,
                                                    nbll);
                            }
                        } // end do m=1,ncg

                        lead_ns::lead(nbl, lw, lw2, maxbl);

                    } // end if ncg > 0 && mgflag > 0
                } // end if idefrm == 0
            } // end if iuns > 0

        } // end if myid == mblk2nd(nbl) ...

        if (isklton == 1) {
            int iunit11 = 11;
            writ_buf_ns::writ_buf(nbl, iunit11, nou, bou, nbuf, ibufdim,
                                  myhost, myid, mycomm, mblk2nd, maxbl);
        }

    } // end do nbl=1,nblock



    // update moment center location
    if (myid == myhost) {
        if (itransmc > 0) {

            if (isklton == 1) {
                std::fprintf(fortran_get_unit(11),
                             " translating moment center to new position\n");
            }

            // Cast float COMMON fields to double for transmc call
            double d_rfreqtmc  = (double)rfreqtmc;
            double d_utransmc  = (double)utransmc;
            double d_vtransmc  = (double)vtransmc;
            double d_wtransmc  = (double)wtransmc;
            double d_xorigmc   = (double)xorigmc;
            double d_yorigmc   = (double)yorigmc;
            double d_zorigmc   = (double)zorigmc;
            double d_xorig0mc  = (double)xorig0mc;
            double d_yorig0mc  = (double)yorig0mc;
            double d_zorig0mc  = (double)zorig0mc;
            double d_xmc       = (double)xmc;
            double d_ymc       = (double)ymc;
            double d_zmc       = (double)zmc;
            double d_time2mc   = (double)time2mc;
            int    i_itransmc  = (int)itransmc;
            transmc_ns::transmc(i_itransmc, d_rfreqtmc,
                                d_utransmc, d_vtransmc, d_wtransmc,
                                d_xorigmc, d_yorigmc, d_zorigmc,
                                d_xorig0mc, d_yorig0mc, d_zorig0mc,
                                d_xmc, d_ymc, d_zmc,
                                iupdat, d_time2mc);
            // Write back modified values
            xorigmc  = (float)d_xorigmc;
            yorigmc  = (float)d_yorigmc;
            zorigmc  = (float)d_zorigmc;
            xmc      = (float)d_xmc;
            ymc      = (float)d_ymc;
            zmc      = (float)d_zmc;
            time2mc  = (float)d_time2mc;
        }

        if (irotatmc > 0) {

            if (isklton == 1) {
                std::fprintf(fortran_get_unit(11),
                             "  rotating moment center to new position\n");
            }

            // Cast float COMMON fields to double for rotatmc call
            double d_rfreqrmc  = (double)rfreqrmc;
            double d_omegaxmc  = (double)omegaxmc;
            double d_omegaymc  = (double)omegaymc;
            double d_omegazmc  = (double)omegazmc;
            double d_xorigmc   = (double)xorigmc;
            double d_yorigmc   = (double)yorigmc;
            double d_zorigmc   = (double)zorigmc;
            double d_thetaxmc  = (double)thetaxmc;
            double d_thetaymc  = (double)thetaymc;
            double d_thetazmc  = (double)thetazmc;
            double d_xmc       = (double)xmc;
            double d_ymc       = (double)ymc;
            double d_zmc       = (double)zmc;
            double d_time2mc   = (double)time2mc;
            int    i_irotatmc  = (int)irotatmc;
            rotatmc_ns::rotatmc(i_irotatmc, d_rfreqrmc,
                                d_omegaxmc, d_omegaymc, d_omegazmc,
                                d_xorigmc, d_yorigmc, d_zorigmc,
                                d_thetaxmc, d_thetaymc, d_thetazmc,
                                d_xmc, d_ymc, d_zmc,
                                iupdat, d_time2mc);
            // Write back modified values
            xorigmc  = (float)d_xorigmc;
            yorigmc  = (float)d_yorigmc;
            zorigmc  = (float)d_zorigmc;
            thetaxmc = (float)d_thetaxmc;
            thetaymc = (float)d_thetaymc;
            thetazmc = (float)d_thetazmc;
            xmc      = (float)d_xmc;
            ymc      = (float)d_ymc;
            zmc      = (float)d_zmc;
            time2mc  = (float)d_time2mc;
        }
    }

    return;
}


} // namespace updateg_ns
