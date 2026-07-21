// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc_blkint.h"
#include "bc.h"
#include "lead.h"
#include "blockk.h"
#include "blockj.h"
#include "blocki.h"
#include "cblkk.h"
#include "cblkj.h"
#include "cblki.h"
#include <cmath>
#include <cstdlib>

namespace bc_blkint_ns {

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
    bc_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, cl, nou, bou,
              nbuf, ibufdim, maxbl, maxgr, maxseg, itrans, irotat, idefrm,
              igridg, nblg, nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
              ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek, lwdat,
              myid, idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
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
    // epsc0....tolerance on 1-1 geometric mismatch
    epsc0 = std::max(1.e-09, std::pow(10.0, (double)(-cmn_zero.iexp + 1)));

    if (ntime > 0 && std::abs(nbli) > 0) {

        // First Case: all data needed to set 1-1 bc lies on the current processor

        int level = cmn_mgrd.level;

        for (int lcnt = cmn_is_blockbc.is_blk[level - 1];
             lcnt <= cmn_is_blockbc.ie_blk[level - 1]; lcnt++) {

            // ic_blk is current (to) block
            // in_blk is neighbor (from) block
            int ic_blk  = isav_blk(lcnt, 4);
            int in_blk  = isav_blk(lcnt, 5);
            int nd_dest = mblk2nd(ic_blk);
            int nd_srce = mblk2nd(in_blk);

            if (nd_srce == myid && nd_dest == myid) {

                if (iadvance(ic_blk) >= 0) {

                    lead_ns::lead(ic_blk, lw, lw2, maxbl);

                    int n      = isav_blk(lcnt, 1);
                    int it     = isav_blk(lcnt, 2);
                    int ir     = isav_blk(lcnt, 3);
                    int jface  = isav_blk(lcnt, 6);
                    int jedge  = isav_blk(lcnt, 7);
                    int jside  = isav_blk(lcnt, 8);
                    int lwt    = isav_blk(lcnt, 9);
                    int iedge  = isav_blk(lcnt, 10);
                    int iss    = isav_blk(lcnt, 11);
                    int ise    = isav_blk(lcnt, 12);
                    int jss    = isav_blk(lcnt, 13);
                    int jse    = isav_blk(lcnt, 14);
                    int kss    = isav_blk(lcnt, 15);
                    int kse    = isav_blk(lcnt, 16);
                    int iti    = isav_blk(lcnt, 17);

                    // set dimensions of blocks involved
                    int idimn = idimg(in_blk);
                    int jdimn = jdimg(in_blk);
                    int kdimn = kdimg(in_blk);
                    int idimc = idimg(ic_blk);
                    int jdimc = jdimg(ic_blk);
                    int kdimc = kdimg(ic_blk);

                    // set pointers for neighboring block
                    int lws  = lw(1,  in_blk);
                    int lwxr = lw(10, in_blk);
                    int lwyr = lw(11, in_blk);
                    int lwzr = lw(12, in_blk);

                    // set pointers for current block (note: lwt has
                    // already been set above, and corresponds to
                    // qi0, qj0, or qk0, depending on the interface)
                    int lwxt = lw(10, ic_blk);
                    int lwyt = lw(11, ic_blk);
                    int lwzt = lw(12, ic_blk);

                    // k = constant interface
                    if (jside == 3) {

                        // check geometric mismatch
                        if (cmn_sklton.isklton == 1) {
                            FortranArray3DRef<double> xr(&w(lwxr), jdimn, kdimn, idimn);
                            FortranArray3DRef<double> yr(&w(lwyr), jdimn, kdimn, idimn);
                            FortranArray3DRef<double> zr(&w(lwzr), jdimn, kdimn, idimn);
                            FortranArray3DRef<double> xt(&w(lwxt), jdimc, kdimc, idimc);
                            FortranArray3DRef<double> yt(&w(lwyt), jdimc, kdimc, idimc);
                            FortranArray3DRef<double> zt(&w(lwzt), jdimc, kdimc, idimc);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            cblkk_ns::cblkk(nbli, idimn, jdimn, kdimn, idimc, jdimc, kdimc,
                                            limblk_n, isva_n, it, ir, iedge,
                                            xr, yr, zr, xt, yt, zt,
                                            iti, lcnt, geom_miss, mxbli);
                        }

                        // interpolate q
                        {
                            int ldim = 5;
                            FortranArray4DRef<double> q_r(&w(lws), jdimn, kdimn, idimn, 5);
                            FortranArray4DRef<double> qk0(&w(lwt), jdimc, idimc - 1, ldim, 4);
                            FortranArray3DRef<double> bck(&w(cmn_ginfo.lbck), jdimc, idimc - 1, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 0;
                            blockk_ns::blockk(q_r, qk0, idimn, jdimn, kdimn,
                                              idimc, jdimc, limblk_n, isva_n,
                                              it, ir, ldim, bck, iedge, ivolflg);
                        }

                        // interpolate turbulent quantities
                        if (cmn_maxiv.ivmx >= 2) {
                            int ldim = 1;
                            int lwst = lw(13, in_blk);
                            int lwtt = lw(29, ic_blk);
                            if (limblk(it, 3, n) != 1) {
                                lwtt = lwtt + jdimc * (idimc - 1) * 1 * 2;
                            }
                            FortranArray4DRef<double> q_r(&w(lwst), jdimn, kdimn, idimn, ldim);
                            FortranArray4DRef<double> qk0(&w(lwtt), jdimc, idimc - 1, ldim, 4);
                            FortranArray3DRef<double> bck(&w(cmn_ginfo.lbck), jdimc, idimc - 1, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 0;
                            blockk_ns::blockk(q_r, qk0, idimn, jdimn, kdimn,
                                              idimc, jdimc, limblk_n, isva_n,
                                              it, ir, ldim, bck, iedge, ivolflg);
                        }
                        if (cmn_maxiv.ivmx >= 4) {
                            int ldim = nummem;
                            int lwst = lw(19, in_blk);
                            int lwtt = lw(24, ic_blk);
                            if (limblk(it, 3, n) != 1) {
                                lwtt = lwtt + jdimc * (idimc - 1) * 2 * nummem;
                            }
                            FortranArray4DRef<double> q_r(&w(lwst), jdimn, kdimn, idimn, ldim);
                            FortranArray4DRef<double> qk0(&w(lwtt), jdimc, idimc - 1, ldim, 4);
                            FortranArray3DRef<double> bck(&w(cmn_ginfo.lbck), jdimc, idimc - 1, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 0;
                            blockk_ns::blockk(q_r, qk0, idimn, jdimn, kdimn,
                                              idimc, jdimc, limblk_n, isva_n,
                                              it, ir, ldim, bck, iedge, ivolflg);
                        }

                        // interpolate cell volumes
                        if (cmn_maxiv.ivmx >= 1 && cmn_is_blockbc.ivolint > 0) {
                            int ldim   = 1;
                            int lvol   = lw(8,  in_blk);
                            int lvolk0 = lw(50, ic_blk);
                            if (limblk(it, 3, n) != 1) {
                                lvolk0 = lvolk0 + jdimc * (idimc - 1) * 1 * 2;
                            }
                            int idimn_m1 = idimn - 1;
                            FortranArray4DRef<double> q_r(&w(lvol), jdimn, kdimn, idimn_m1, ldim);
                            FortranArray4DRef<double> qk0(&w(lvolk0), jdimc, idimc - 1, ldim, 4);
                            FortranArray3DRef<double> bck(&w(cmn_ginfo.lbck), jdimc, idimc - 1, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 1;
                            blockk_ns::blockk(q_r, qk0, idimn_m1, jdimn, kdimn,
                                              idimc, jdimc, limblk_n, isva_n,
                                              it, ir, ldim, bck, iedge, ivolflg);
                        }

                    // j = constant interface
                    } else if (jside == 4) {

                        // check geometric mismatch
                        if (cmn_sklton.isklton == 1) {
                            FortranArray3DRef<double> xr(&w(lwxr), jdimn, kdimn, idimn);
                            FortranArray3DRef<double> yr(&w(lwyr), jdimn, kdimn, idimn);
                            FortranArray3DRef<double> zr(&w(lwzr), jdimn, kdimn, idimn);
                            FortranArray3DRef<double> xt(&w(lwxt), jdimc, kdimc, idimc);
                            FortranArray3DRef<double> yt(&w(lwyt), jdimc, kdimc, idimc);
                            FortranArray3DRef<double> zt(&w(lwzt), jdimc, kdimc, idimc);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            cblkj_ns::cblkj(nbli, idimn, jdimn, kdimn, idimc, jdimc, kdimc,
                                            limblk_n, isva_n, it, ir, iedge,
                                            xr, yr, zr, xt, yt, zt,
                                            iti, lcnt, geom_miss, mxbli);
                        }

                        // interpolate q
                        {
                            int ldim = 5;
                            FortranArray4DRef<double> q_r(&w(lws), jdimn, kdimn, idimn, 5);
                            FortranArray4DRef<double> qj0(&w(lwt), kdimc, idimc - 1, ldim, 4);
                            FortranArray3DRef<double> bcj(&w(cmn_ginfo.lbcj), kdimc, idimc - 1, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 0;
                            blockj_ns::blockj(q_r, qj0, idimn, jdimn, kdimn,
                                              idimc, kdimc, limblk_n, isva_n,
                                              it, ir, ldim, bcj, iedge, ivolflg);
                        }

                        // interpolate turbulent quantities
                        if (cmn_maxiv.ivmx >= 2) {
                            int ldim = 1;
                            int lwst = lw(13, in_blk);
                            int lwtt = lw(28, ic_blk);
                            if (limblk(it, 2, n) != 1) {
                                lwtt = lwtt + kdimc * (idimc - 1) * 1 * 2;
                            }
                            FortranArray4DRef<double> q_r(&w(lwst), jdimn, kdimn, idimn, ldim);
                            FortranArray4DRef<double> qj0(&w(lwtt), kdimc, idimc - 1, ldim, 4);
                            FortranArray3DRef<double> bcj(&w(cmn_ginfo.lbcj), kdimc, idimc - 1, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 0;
                            blockj_ns::blockj(q_r, qj0, idimn, jdimn, kdimn,
                                              idimc, kdimc, limblk_n, isva_n,
                                              it, ir, ldim, bcj, iedge, ivolflg);
                        }
                        if (cmn_maxiv.ivmx >= 4) {
                            int ldim = nummem;
                            int lwst = lw(19, in_blk);
                            int lwtt = lw(23, ic_blk);
                            if (limblk(it, 2, n) != 1) {
                                lwtt = lwtt + kdimc * (idimc - 1) * 2 * nummem;
                            }
                            FortranArray4DRef<double> q_r(&w(lwst), jdimn, kdimn, idimn, ldim);
                            FortranArray4DRef<double> qj0(&w(lwtt), kdimc, idimc - 1, ldim, 4);
                            FortranArray3DRef<double> bcj(&w(cmn_ginfo.lbcj), kdimc, idimc - 1, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 0;
                            blockj_ns::blockj(q_r, qj0, idimn, jdimn, kdimn,
                                              idimc, kdimc, limblk_n, isva_n,
                                              it, ir, ldim, bcj, iedge, ivolflg);
                        }

                        // interpolate cell volumes
                        if (cmn_maxiv.ivmx >= 1 && cmn_is_blockbc.ivolint > 0) {
                            int ldim   = 1;
                            int lvol   = lw(8,  in_blk);
                            int lvolj0 = lw(49, ic_blk);
                            if (limblk(it, 2, n) != 1) {
                                lvolj0 = lvolj0 + kdimc * (idimc - 1) * 1 * 2;
                            }
                            int idimn_m1 = idimn - 1;
                            FortranArray4DRef<double> q_r(&w(lvol), jdimn, kdimn, idimn_m1, ldim);
                            FortranArray4DRef<double> qj0(&w(lvolj0), kdimc, idimc - 1, ldim, 4);
                            FortranArray3DRef<double> bcj(&w(cmn_ginfo.lbcj), kdimc, idimc - 1, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 1;
                            blockj_ns::blockj(q_r, qj0, idimn_m1, jdimn, kdimn,
                                              idimc, kdimc, limblk_n, isva_n,
                                              it, ir, ldim, bcj, iedge, ivolflg);
                        }

                    // i = constant interface
                    } else if (jside == 5) {

                        // check geometric mismatch
                        if (cmn_sklton.isklton == 1) {
                            FortranArray3DRef<double> xr(&w(lwxr), jdimn, kdimn, idimn);
                            FortranArray3DRef<double> yr(&w(lwyr), jdimn, kdimn, idimn);
                            FortranArray3DRef<double> zr(&w(lwzr), jdimn, kdimn, idimn);
                            FortranArray3DRef<double> xt(&w(lwxt), jdimc, kdimc, idimc);
                            FortranArray3DRef<double> yt(&w(lwyt), jdimc, kdimc, idimc);
                            FortranArray3DRef<double> zt(&w(lwzt), jdimc, kdimc, idimc);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            cblki_ns::cblki(nbli, idimn, jdimn, kdimn, idimc, jdimc, kdimc,
                                            limblk_n, isva_n, it, ir, iedge,
                                            xr, yr, zr, xt, yt, zt,
                                            iti, lcnt, geom_miss, mxbli);
                        }

                        // interpolate q
                        {
                            int ldim = 5;
                            FortranArray4DRef<double> q_r(&w(lws), jdimn, kdimn, idimn, 5);
                            FortranArray4DRef<double> qi0(&w(lwt), jdimc, kdimc, ldim, 4);
                            FortranArray3DRef<double> bci(&w(cmn_ginfo.lbci), jdimc, kdimc, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 0;
                            blocki_ns::blocki(q_r, qi0, idimn, jdimn, kdimn,
                                              jdimc, kdimc, limblk_n, isva_n,
                                              it, ir, ldim, bci, iedge, ivolflg);
                        }

                        // interpolate turbulent quantities
                        if (cmn_maxiv.ivmx >= 2) {
                            int ldim = 1;
                            int lwst = lw(13, in_blk);
                            int lwtt = lw(30, ic_blk);
                            if (limblk(it, 1, n) != 1) {
                                lwtt = lwtt + jdimc * kdimc * 1 * 2;
                            }
                            FortranArray4DRef<double> q_r(&w(lwst), jdimn, kdimn, idimn, ldim);
                            FortranArray4DRef<double> qi0(&w(lwtt), jdimc, kdimc, ldim, 4);
                            FortranArray3DRef<double> bci(&w(cmn_ginfo.lbci), jdimc, kdimc, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 0;
                            blocki_ns::blocki(q_r, qi0, idimn, jdimn, kdimn,
                                              jdimc, kdimc, limblk_n, isva_n,
                                              it, ir, ldim, bci, iedge, ivolflg);
                        }
                        if (cmn_maxiv.ivmx >= 4) {
                            int ldim = nummem;
                            int lwst = lw(19, in_blk);
                            int lwtt = lw(25, ic_blk);
                            if (limblk(it, 1, n) != 1) {
                                lwtt = lwtt + jdimc * kdimc * 2 * nummem;
                            }
                            FortranArray4DRef<double> q_r(&w(lwst), jdimn, kdimn, idimn, ldim);
                            FortranArray4DRef<double> qi0(&w(lwtt), jdimc, kdimc, ldim, 4);
                            FortranArray3DRef<double> bci(&w(cmn_ginfo.lbci), jdimc, kdimc, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 0;
                            blocki_ns::blocki(q_r, qi0, idimn, jdimn, kdimn,
                                              jdimc, kdimc, limblk_n, isva_n,
                                              it, ir, ldim, bci, iedge, ivolflg);
                        }

                        // interpolate cell volumes
                        if (cmn_maxiv.ivmx >= 1 && cmn_is_blockbc.ivolint > 0) {
                            int ldim   = 1;
                            int lvol   = lw(8,  in_blk);
                            int lvoli0 = lw(51, ic_blk);
                            if (limblk(it, 1, n) != 1) {
                                lvoli0 = lvoli0 + jdimc * kdimc * 1 * 2;
                            }
                            int idimn_m1 = idimn - 1;
                            FortranArray4DRef<double> q_r(&w(lvol), jdimn, kdimn, idimn_m1, ldim);
                            FortranArray4DRef<double> qi0(&w(lvoli0), jdimc, kdimc, ldim, 4);
                            FortranArray3DRef<double> bci(&w(cmn_ginfo.lbci), jdimc, kdimc, 2);
                            FortranArray2DRef<int> limblk_n(&limblk(1,1,n), 2, 6);
                            FortranArray2DRef<int> isva_n(&isva(1,1,n), 2, 2);
                            int ivolflg = 1;
                            blocki_ns::blocki(q_r, qi0, idimn_m1, jdimn, kdimn,
                                              jdimc, kdimc, limblk_n, isva_n,
                                              it, ir, ldim, bci, iedge, ivolflg);
                        }

                    } // end jside

                } // end iadvance >= 0

            } // end nd_srce == myid && nd_dest == myid

        } // end lcnt loop

        // Second Case: data needed to set 1-1 bc lies on another processor
        // (this case is handled by MPI communication outside this routine)

    } // end ntime > 0 && abs(nbli) > 0

    return;
}

} // namespace bc_blkint_ns
