// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "updatedg.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>

// External procedure headers
#include "lead.h"
#include "writ_buf.h"
#include "trnsurf.h"
#include "rotsurf.h"
#include "aesurf.h"
#include "getsurf.h"
// deform_surf and deform_el are declared in deform.h
#include "deform.h"
#include "xtbatb.h"
#include "metric.h"
#include "tmetric.h"
#include "cellvol.h"
#include "collx.h"
#include "collxt.h"
#include "collv.h"
#include "collxtb.h"
#include "termn8.h"
#include "invert.h"

namespace updatedg_ns {

// ============================================================
// updatedg
// ============================================================
void updatedg(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
              FortranArray1DRef<double> w, int& mgwk,
              FortranArray1DRef<double> wk, int& nwork,
              int& iupdat, int& iseqr, int& maxbl, int& maxgr, int& maxseg,
              FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0,
              FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim,
              FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim,
              FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo,
              FortranArray4DRef<int> kbcinfo,
              int& nblock, FortranArray1DRef<int> levelg,
              FortranArray1DRef<int> igridg, FortranArray1DRef<int> idefrm,
              FortranArray1DRef<int> ncgg, FortranArray1DRef<int> iadvance,
              FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
              int& nbuf, int& ibufdim, int& myid, int& myhost, int& mycomm,
              FortranArray1DRef<int> mblk2nd,
              FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae,
              FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae,
              FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae,
              FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae,
              FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae,
              FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae,
              FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae,
              FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
              FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
              FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
              FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass,
              FortranArray2DRef<double> damp, FortranArray2DRef<double> x0,
              FortranArray2DRef<double> gf0, int& nmds, int& maxaes,
              FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb,
              FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
              FortranArray3DRef<int> islavept, int& nslave,
              FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip,
              FortranArray2DRef<int> kskip,
              FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn,
              FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg,
              FortranArray2DRef<int> iaesurf, int& maxsegdg,
              FortranArray1DRef<int> iwk, int& nmaster, int& nt,
              FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig,
              FortranArray1DRef<double> zorig,
              FortranArray2DRef<double> xorgae0, FortranArray2DRef<double> yorgae0,
              FortranArray2DRef<double> zorgae0,
              FortranArray2DRef<int> icouple, FortranArray1DRef<int> ireq,
              int& nnodes, FortranArray2DRef<int> nblelst,
              FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax,
              FortranArray1DRef<int> kskmax, FortranArray1DRef<double> ue)
{

    // COMMON block aliases
    int& jdim    = cmn_ginfo.jdim;
    int& kdim    = cmn_ginfo.kdim;
    int& idim    = cmn_ginfo.idim;
    int& jj2     = cmn_ginfo.jj2;
    int& kk2     = cmn_ginfo.kk2;
    int& ii2     = cmn_ginfo.ii2;
    int& lglobal = cmn_mgrd.lglobal;
    int& isklton = cmn_sklton.isklton;
    float& time_f  = cmn_unst.time;
    int& iunst   = cmn_unst.iunst;
    int& ita     = cmn_unst.ita;
    float& dt_f  = cmn_info.dt;
    int& mseq    = cmn_info.mseq;
    int& mgflag  = cmn_info.mgflag;
    int& idef_ss = cmn_elastic_ss.idef_ss;
    int* levelt  = cmn_info.levelt; // 0-based C array, Fortran 1-based

    // ginfo pointers
    int& lx    = cmn_ginfo.lx;
    int& ly    = cmn_ginfo.ly;
    int& lz    = cmn_ginfo.lz;
    int& ldeltj= cmn_ginfo.ldeltj;
    int& ldeltk= cmn_ginfo.ldeltk;
    int& ldelti= cmn_ginfo.ldelti;
    int& lxmdj = cmn_ginfo.lxmdj;
    int& lxmdk = cmn_ginfo.lxmdk;
    int& lxmdi = cmn_ginfo.lxmdi;
    int& lxnm1 = cmn_ginfo.lxnm1;
    int& lynm1 = cmn_ginfo.lynm1;
    int& lznm1 = cmn_ginfo.lznm1;
    int& lxnm2 = cmn_ginfo.lxnm2;
    int& lynm2 = cmn_ginfo.lynm2;
    int& lznm2 = cmn_ginfo.lznm2;
    int& lsj   = cmn_ginfo.lsj;
    int& lsk   = cmn_ginfo.lsk;
    int& lsi   = cmn_ginfo.lsi;
    int& lvol  = cmn_ginfo.lvol;
    int& lxtbj = cmn_ginfo.lxtbj;
    int& lxtbk = cmn_ginfo.lxtbk;
    int& lxtbi = cmn_ginfo.lxtbi;
    int& latbj = cmn_ginfo.latbj;
    int& latbk = cmn_ginfo.latbk;
    int& latbi = cmn_ginfo.latbi;

    double time = (double)time_f;
    double dt   = (double)dt_f;

    FortranArray1D<double> xst(nslave);
    FortranArray1D<double> yst(nslave);
    FortranArray1D<double> zst(nslave);

    int icnt=0, nbl, is, lll, nsurf=0, nblmast;
    int lwk1, lwk2, lwk3, lwk4;
    int lvel, lacci, laccj, lacck, lt1wk, lt2wk, lt3wk, lt4wk;
    int mdim, nroom, iflag, iflagv=0, imin=0, imax=0, jmin=0, jmax=0, kmin=0, kmax=0;
    int idim1, jdim1, kdim1;
    int nflag, irst;
    int ncg, nbll, nbllm1, m;
    int lvolc, lxc, lyc, lzc, lxtbjc, lxtbkc, lxtbic, latbjc, latbkc, latbic, lvelc;
    int nv, izz;
    int Isys;
    char cmd1[36], cmd2[36], char1[4];

    if (isklton == 1) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::memset(bou(nou(1),1), ' ', 120);
        int writ_nbl1 = 1, writ_iunit11 = 11;
        writ_buf_ns::writ_buf(writ_nbl1, writ_iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl);
    }

    // update surface positions due to forced motion and/or aeroelastic motion
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (myid == mblk2nd(nbl) &&
            (levelg(nbl) >= lglobal && levelg(nbl) <= levelt[iseqr-1])) {

            lead_ns::lead(nbl, lw, lw2, maxbl);

            // zero out deltj, deltk, delti from previous step
            if (idefrm(nbl) > 0 && idef_ss == 0) {
                for (lll = 1; lll <= kdim*idim*3*2; lll++)
                    w(ldeltj+lll-1) = 0.;
                for (lll = 1; lll <= jdim*idim*3*2; lll++)
                    w(ldeltk+lll-1) = 0.;
                for (lll = 1; lll <= kdim*jdim*3*2; lll++)
                    w(ldelti+lll-1) = 0.;
            }

            if (idefrm(nbl) > 0 && idefrm(nbl) < 999) {
                if (isklton == 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " deforming surface of block  %4d to new shape", nbl);
                }

                lwk1 = 1;
                lwk2 = lwk1 + kdim*idim*2;
                lwk3 = lwk2 + jdim*idim*2;
                lwk4 = lwk3 + jdim*kdim*2;
                if (nwork < lwk4) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " not enough work space for subroutines trnsurf/rotsurf/aesurf");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " Increase keyword memadd allocation by %8d", lwk4 - nwork);
                    { int m1=-1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
                }

                for (lll = lwk1; lll <= lwk4; lll++)
                    wk(lll) = 1.;

                // update segment origin if coupled
                for (is = 1; is <= nsegdfrm(nbl); is++) {
                    if (icouple(nbl,is) != 0) {
                        if (idfrmseg(nbl,is) > 0 && idfrmseg(nbl,is) < 99) {
                            if (icouple(nbl,is) != 0) {
                                nblmast = icouple(nbl,is);
                                xorgae(nbl,is) = xorig(nblmast);
                                yorgae(nbl,is) = yorig(nblmast);
                                zorgae(nbl,is) = zorig(nblmast);
                            }
                        }
                    }
                }

                // surface translation
                {
                    FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                    FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                    FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                    FortranArray4DRef<double> wdeltj(&w(ldeltj), kdim, idim, 3, 2);
                    FortranArray4DRef<double> wdeltk(&w(ldeltk), jdim, idim, 3, 2);
                    FortranArray4DRef<double> wdelti(&w(ldelti), kdim, jdim, 3, 2);
                    FortranArray3DRef<double> wwk1(&wk(lwk1), kdim, idim, 2);
                    FortranArray3DRef<double> wwk2(&wk(lwk2), jdim, idim, 2);
                    FortranArray3DRef<double> wwk3(&wk(lwk3), jdim, kdim, 2);
                    for (is = 1; is <= nsegdfrm(nbl); is++) {
                        if (idfrmseg(nbl,is) == 1) {
                            trnsurf_ns::trnsurf(jdim, kdim, idim,
                                wx, wy, wz,
                                wdeltj, wdeltk, wdelti,
                                nbl, idfrmseg(nbl,is),
                                xorgae(nbl,is), yorgae(nbl,is), zorgae(nbl,is),
                                utrnsae(nbl,is), vtrnsae(nbl,is), wtrnsae(nbl,is),
                                rfrqtae(nbl,is),
                                icsi(nbl,is), icsf(nbl,is), jcsi(nbl,is),
                                jcsf(nbl,is), kcsi(nbl,is), kcsf(nbl,is),
                                time, nou, bou, nbuf, ibufdim, myid,
                                wwk1, wwk2, wwk3,
                                xorgae0(nbl,is), yorgae0(nbl,is), zorgae0(nbl,is));
                        }
                    }
                }

                for (lll = lwk1; lll <= lwk4; lll++)
                    wk(lll) = 1.;

                // surface rotation
                {
                    FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                    FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                    FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                    FortranArray4DRef<double> wdeltj(&w(ldeltj), kdim, idim, 3, 2);
                    FortranArray4DRef<double> wdeltk(&w(ldeltk), jdim, idim, 3, 2);
                    FortranArray4DRef<double> wdelti(&w(ldelti), kdim, jdim, 3, 2);
                    FortranArray3DRef<double> wwk1(&wk(lwk1), kdim, idim, 2);
                    FortranArray3DRef<double> wwk2(&wk(lwk2), jdim, idim, 2);
                    FortranArray3DRef<double> wwk3(&wk(lwk3), jdim, kdim, 2);
                    for (is = 1; is <= nsegdfrm(nbl); is++) {
                        if (idfrmseg(nbl,is) == 2) {
                            rotsurf_ns::rotsurf(jdim, kdim, idim,
                                wx, wy, wz,
                                wdeltj, wdeltk, wdelti,
                                nbl, idfrmseg(nbl,is),
                                xorgae(nbl,is), yorgae(nbl,is), zorgae(nbl,is),
                                omgxae(nbl,is), omgyae(nbl,is), omgzae(nbl,is),
                                thtxae(nbl,is), thtyae(nbl,is), thtzae(nbl,is),
                                rfrqrae(nbl,is),
                                icsi(nbl,is), icsf(nbl,is), jcsi(nbl,is),
                                jcsf(nbl,is), kcsi(nbl,is), kcsf(nbl,is),
                                time, nou, bou, nbuf, ibufdim, myid,
                                wwk1, wwk2, wwk3);
                        }
                    }
                }

                for (lll = lwk1; lll <= lwk4; lll++)
                    wk(lll) = 1.;

                // aeroelastic surface
                {
                    FortranArray4DRef<double> wdeltj(&w(ldeltj), kdim, idim, 3, 2);
                    FortranArray4DRef<double> wdeltk(&w(ldeltk), jdim, idim, 3, 2);
                    FortranArray4DRef<double> wdelti(&w(ldelti), kdim, jdim, 3, 2);
                    FortranArray5DRef<double> wxmdj(&w(lxmdj), jdim, kdim, idim, 1, 1);
                    FortranArray5DRef<double> wxmdk(&w(lxmdk), jdim, kdim, idim, 1, 1);
                    FortranArray5DRef<double> wxmdi(&w(lxmdi), jdim, kdim, idim, 1, 1);
                    FortranArray3DRef<double> wwk1(&wk(lwk1), kdim, idim, 2);
                    FortranArray3DRef<double> wwk2(&wk(lwk2), jdim, idim, 2);
                    FortranArray3DRef<double> wwk3(&wk(lwk3), jdim, kdim, 2);
                    aesurf_ns::aesurf(nbl, jdim, kdim, idim,
                        wdeltj, wdeltk, wdelti,
                        wxmdj, wxmdk, wxmdi,
                        wwk1, wwk2, wwk3,
                        maxbl, maxseg, nmds, maxaes, aesrfdat, xs, xxn,
                        icsi, icsf, jcsi, jcsf, kcsi, kcsf,
                        nsegdfrm, idfrmseg, iaesurf, maxsegdg);
                }

                for (lll = lwk1; lll <= lwk4; lll++)
                    wk(lll) = 1.;
            }
        }

        if (isklton == 1) {
            int iunit_tmp = 11;
            writ_buf_ns::writ_buf(nbl, iunit_tmp, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl);
        }
    } // end do 100

    // update block origin if coupled to forced deformation
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (myid == mblk2nd(nbl) &&
            (levelg(nbl) >= lglobal && levelg(nbl) <= levelt[iseqr-1])) {
            for (is = 1; is <= nsegdfrm(nbl); is++) {
                if (icouple(nbl,is) != 0) {
                    nblmast = icouple(nbl,is);
                    if (idfrmseg(nblmast,is) == 1) {
                        xorig(nbl) = xorgae(nblmast,is);
                        yorig(nbl) = yorgae(nblmast,is);
                        zorig(nbl) = zorgae(nblmast,is);
                    }
                }
            }
        }
    }

    // get list of all points on deforming solid surfaces
    nsurf = 0;
    for (nbl = 1; nbl <= nblock; nbl++) {
        iwk(nbl) = 0;
        if ((levelg(nbl) >= lglobal && levelg(nbl) <= levelt[iseqr-1])) {
            if (idefrm(nbl) > 0 && idefrm(nbl) < 999) {
                lead_ns::lead(nbl, lw, lw2, maxbl);
                if (myid == mblk2nd(nbl)) {
                    FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                    FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                    FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                    FortranArray4DRef<double> wdeltj(&w(ldeltj), kdim, idim, 3, 2);
                    FortranArray4DRef<double> wdeltk(&w(ldeltk), jdim, idim, 3, 2);
                    FortranArray4DRef<double> wdelti(&w(ldelti), kdim, jdim, 3, 2);
                    FortranArray3DRef<double> wxnm1(&w(lxnm1), jdim, kdim, idim);
                    FortranArray3DRef<double> wynm1(&w(lynm1), jdim, kdim, idim);
                    FortranArray3DRef<double> wznm1(&w(lznm1), jdim, kdim, idim);
                    getsurf_ns::getsurf(wx, wy, wz,
                        wdeltj, wdeltk, wdelti,
                        wxnm1, wynm1, wznm1,
                        icsi, icsf, jcsi, jcsf, kcsi, kcsf,
                        wk, nwork, nbl, idim, jdim, kdim, nsurf, iwk(nbl),
                        nsegdfrm, maxbl, idfrmseg, maxsegdg);
                }
                if (myid == myhost) {
                    nsurf = nsurf + iwk(nbl);
                }
            }
        }
    }

    if (iunst > 1 || idef_ss > 0) {
        for (nbl = 1; nbl <= nblock; nbl++) {
            if (myid == mblk2nd(nbl) &&
                (levelg(nbl) >= lglobal && levelg(nbl) <= levelt[iseqr-1])) {

                lead_ns::lead(nbl, lw, lw2, maxbl);

                if (idefrm(nbl) > 0) {
                    lvel  = 9*nsurf + 1;
                    lacci = jdim*kdim*idim*3 + lvel;
                    laccj = jdim*kdim*3*2 + lacci;
                    lacck = kdim*idim*3*2 + laccj;
                    lt1wk = jdim*idim*3*2 + lacck;
                    lt2wk = jdim*kdim*idim*3 + lt1wk;
                    lt3wk = jdim*kdim*6 + lt2wk;
                    lt4wk = jdim*kdim*idim*5 + lt3wk;

                    mdim = jdim*kdim*idim*3;
                    if (nwork < mdim) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " not enough work space for subroutine deform");
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " Increase keyword memadd allocation by %8d", mdim - nwork);
                        { int m1=-1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
                    }

                    if (isklton == 1) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " deforming block  %4d to new shape", nbl);
                    }

                    nflag = 0;
                    irst  = 0;
                    if (idef_ss > 0) irst = 1;
                    FortranArray4DRef<double> wdeltj_ds(&w(ldeltj), kdim, idim, 3, 2);
                    FortranArray4DRef<double> wdeltk_ds(&w(ldeltk), jdim, idim, 3, 2);
                    FortranArray4DRef<double> wdelti_ds(&w(ldelti), kdim, jdim, 3, 2);
                    deform_ns::deform_surf(nbl, idim, jdim, kdim,
                        wdeltj_ds, wdeltk_ds, wdelti_ds,
                        lw, lw2, icsi, icsf, jcsi, jcsf, kcsi, kcsf,
                        maxbl, mseq, time, dt, ita, nou, bou, nbuf, ibufdim,
                        myid, idefrm, nbci0, nbcidim, nbcj0, nbcjdim,
                        nbck0, nbckdim, ibcinfo, jbcinfo, kbcinfo, maxseg,
                        wk, ue, nsurf, irst, nflag, islavept, nslave,
                        nsegdfrm, idfrmseg, iaesurf, maxsegdg, nmaster, iseqr);
                }
            }
        } // end do 150
    }


    if (iunst > 1 || idef_ss > 0) {
        for (int n = 1; n <= nslave; n++) {
            nbl = islavept(n, 9, iseqr);
            int ll = islavept(n, 1, iseqr);
            if (myid == mblk2nd(nbl)) {
                lead_ns::lead(nbl, lw, lw2, maxbl);
                xst(n) = w(lx+ll);
                yst(n) = w(ly+ll);
                zst(n) = w(lz+ll);
            }
        }
        deform_ns::deform_el(islavept, nslave, nmaster, ue,
            xst, yst, zst, nt, myhost, mycomm, myid, nnodes,
            mblk2nd, nblelst, maxbl, iseqr);
    }

    if (iunst > 1 || idef_ss > 0) {
        for (nbl = 1; nbl <= nblock; nbl++) {
            if (myid == mblk2nd(nbl) &&
                (levelg(nbl) >= lglobal && levelg(nbl) <= levelt[iseqr-1])) {

                lead_ns::lead(nbl, lw, lw2, maxbl);

                if (idefrm(nbl) > 0) {
                    lvel  = 9*nsurf + 1;
                    lacci = jdim*kdim*idim*3 + lvel;
                    laccj = jdim*kdim*3*2 + lacci;
                    lacck = kdim*idim*3*2 + laccj;
                    lt1wk = jdim*idim*3*2 + lacck;
                    lt2wk = jdim*kdim*idim*3 + lt1wk;
                    lt3wk = jdim*kdim*6 + lt2wk;
                    lt4wk = jdim*kdim*idim*5 + lt3wk;

                    mdim = jdim*kdim*idim*3;
                    if (nwork < mdim) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " not enough work space for subroutine deform");
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " Increase keyword memadd allocation by %8d", mdim - nwork);
                        { int m1=-1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
                    }

                    if (isklton == 1) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " deforming block  %4d to new shape", nbl);
                    }

                    nflag = 0;
                    irst  = 0;
                    if (idef_ss > 0) irst = 1;
                    {
                        FortranArray1DRef<double> wx(&w(lx), jdim*kdim*idim);
                        FortranArray1DRef<double> wy(&w(ly), jdim*kdim*idim);
                        FortranArray1DRef<double> wz(&w(lz), jdim*kdim*idim);
                        FortranArray1DRef<double> wxnm2(&w(lxnm2), jdim*kdim*idim);
                        FortranArray1DRef<double> wynm2(&w(lynm2), jdim*kdim*idim);
                        FortranArray1DRef<double> wznm2(&w(lznm2), jdim*kdim*idim);
                        FortranArray1DRef<double> wxnm1(&w(lxnm1), jdim*kdim*idim);
                        FortranArray1DRef<double> wynm1(&w(lynm1), jdim*kdim*idim);
                        FortranArray1DRef<double> wznm1(&w(lznm1), jdim*kdim*idim);
                        FortranArray4DRef<double> wdeltj(&w(ldeltj), kdim, idim, 3, 2);
                        FortranArray4DRef<double> wdeltk(&w(ldeltk), jdim, idim, 3, 2);
                        FortranArray4DRef<double> wdelti(&w(ldelti), kdim, jdim, 3, 2);
                        FortranArray4DRef<double> wvel(&wk(lvel), jdim, kdim, idim, 3);
                        deform_ns::deform(nbl, idim, jdim, kdim,
                            wx, wy, wz,
                            wxnm2, wynm2, wznm2,
                            wxnm1, wynm1, wznm1,
                            wdeltj, wdeltk, wdelti,
                            ue, wvel,
                            icsi, icsf, jcsi, jcsf, kcsi, kcsf,
                            maxbl, dt, nou, bou, nbuf, ibufdim, myid,
                            idefrm, nbci0, nbcidim, nbcj0, nbcjdim,
                            nbck0, nbckdim, ibcinfo, jbcinfo, kbcinfo,
                            maxseg, wk, nsurf, irst, nflag,
                            islavept, nslave, iskip, jskip, kskip,
                            nsegdfrm, idfrmseg, iaesurf, maxsegdg,
                            nmaster, iseqr, iskmax, jskmax, kskmax, nt);
                    }

                    if (iunst > 1) {
                        for (lll = lacci; lll <= lt1wk-1; lll++)
                            wk(lll) = 0.;

                        {
                            FortranArray4DRef<double> wxtbj(&w(lxtbj), kdim, idim-1, 1, 1);
                            FortranArray4DRef<double> wxtbk(&w(lxtbk), jdim, idim-1, 1, 1);
                            FortranArray4DRef<double> wxtbi(&w(lxtbi), jdim, kdim, 1, 1);
                            FortranArray4DRef<double> watbj(&w(latbj), kdim, idim-1, 1, 1);
                            FortranArray4DRef<double> watbk(&w(latbk), jdim, idim-1, 1, 1);
                            FortranArray4DRef<double> watbi(&w(latbi), jdim, kdim, 1, 1);
                            FortranArray2DRef<double> wvel(&wk(lvel), jdim*kdim*idim, 3);
                            FortranArray3DRef<double> wacci(&wk(lacci), jdim*kdim, 3, 2);
                            FortranArray3DRef<double> waccj(&wk(laccj), kdim*idim, 3, 2);
                            FortranArray3DRef<double> wacck(&wk(lacck), jdim*idim, 3, 2);
                            xtbatb_ns::xtbatb(jdim, kdim, idim,
                                wxtbj, wxtbk, wxtbi,
                                watbj, watbk, watbi,
                                wvel, wacci, waccj, wacck);
                        }

                        nroom = nwork - lt3wk;
                        mdim  = jdim*kdim*idim*5;
                        if (nroom < mdim) {
                            nou(1) = std::min(nou(1)+1, ibufdim);
                            std::snprintf(bou(nou(1),1), 120,
                                " not enough work space for metric subroutines");
                            nou(1) = std::min(nou(1)+1, ibufdim);
                            std::snprintf(bou(nou(1),1), 120,
                                " Increase keyword memadd allocation by %8d", mdim - nroom);
                            { int m1=-1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
                        }

                        {
                            FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                            FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                            FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                            FortranArray3DRef<double> wsj(&w(lsj), jdim, kdim, idim);
                            FortranArray3DRef<double> wsk(&w(lsk), jdim, kdim, idim);
                            FortranArray3DRef<double> wsi(&w(lsi), jdim, kdim, idim);
                            FortranArray2DRef<double> wt2wk(&wk(lt2wk), jdim * kdim, 6);
                            FortranArray3DRef<double> wt3wk(&wk(lt3wk), jdim * kdim, idim, 5);
                            iflag = -1;
                            metric_ns::metric(jdim, kdim, idim,
                                wx, wy, wz, wsj, wsk, wsi,
                                wt2wk, wt3wk, nbl, iflag, icnt,
                                nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
                                ibcinfo, jbcinfo, kbcinfo, maxbl, maxseg,
                                nou, bou, nbuf, ibufdim, myid, mblk2nd);

                            FortranArray2DRef<double> wvel(&wk(lvel), jdim*kdim*idim, 3);
                            FortranArray2DRef<double> wt1wk(&wk(lt1wk), jdim*kdim*idim, 3);
                            tmetric_ns::tmetric(jdim, kdim, idim,
                                wsj, wsk, wsi,
                                wx, wy, wz,
                                wvel, wt1wk, wt2wk, wt3wk, nbl);
                        }

                        nroom = nwork - lt4wk;
                        mdim  = jdim*kdim*15;
                        if (nroom < mdim) {
                            nou(1) = std::min(nou(1)+1, ibufdim);
                            std::snprintf(bou(nou(1),1), 120,
                                " not enough work space for subroutine cellvol");
                            nou(1) = std::min(nou(1)+1, ibufdim);
                            std::snprintf(bou(nou(1),1), 120,
                                " Increase keyword memadd allocation by %8d", mdim - nroom);
                            { int m1=-1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
                        }

                        {
                            FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                            FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                            FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                            FortranArray3DRef<double> wsj(&w(lsj), jdim, kdim, idim);
                            FortranArray3DRef<double> wsk(&w(lsk), jdim, kdim, idim);
                            FortranArray3DRef<double> wsi(&w(lsi), jdim, kdim, idim);
                            FortranArray3DRef<double> wvol(&w(lvol), jdim, kdim, idim);
                            FortranArray2DRef<double> wt4wk(&wk(lt4wk), jdim, kdim);
                            iflagv = 0;
                            int iflagv1 = 1;
                            cellvol_ns::cellvol(jdim, kdim, idim,
                                wx, wy, wz, wsj, wsk, wsi,
                                wvol, wt4wk, nou, bou, nbuf, ibufdim,
                                myid, mblk2nd, maxbl, nbl,
                                iflagv1, iflagv, imin, imax, jmin, jmax, kmin, kmax);
                        }

                        iflagv = 0;
                        if (iflagv == 1) {
                            idim1 = imax - imin + 1;
                            jdim1 = jmax - jmin + 1;
                            kdim1 = kmax - kmin + 1;
                            {
                                FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                                FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                                FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                                s_gwrite(nbl, idim, jdim, kdim,
                                    imin, imax, jmin, jmax, kmin, kmax,
                                    idim1, jdim1, kdim1, wx, wy, wz);
                            }
                            gridgen_glyf1(nbl);
                            gridgen_glyf2(nbl);
                            std::strcpy(cmd1, "gridgen -b dgplot3d_gridgen0001.glf");
                            std::strcpy(cmd2, "gridgen -b dgplot3d_gridgen0002.glf");
                            if (nbl < 10) {
                                std::snprintf(char1, 4, "%1d", nbl);
                                cmd1[29] = char1[0]; cmd2[29] = char1[0];
                            } else if (nbl < 100) {
                                std::snprintf(char1, 4, "%2d", nbl);
                                cmd1[28] = char1[0]; cmd1[29] = char1[1];
                                cmd2[28] = char1[0]; cmd2[29] = char1[1];
                            } else {
                                std::snprintf(char1, 4, "%3d", nbl);
                                cmd1[27] = char1[0]; cmd1[28] = char1[1]; cmd1[29] = char1[2];
                                cmd2[27] = char1[0]; cmd2[28] = char1[1]; cmd2[29] = char1[2];
                            }
                            Isys = system(cmd1);
                            Isys = system(cmd2);
                            {
                                FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                                FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                                FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                                s_gread(nbl, idim, jdim, kdim,
                                    imin, imax, jmin, jmax, kmin, kmax,
                                    idim1, jdim1, kdim1, wx, wy, wz);
                                FortranArray3DRef<double> wsj(&w(lsj), jdim, kdim, idim);
                                FortranArray3DRef<double> wsk(&w(lsk), jdim, kdim, idim);
                                FortranArray3DRef<double> wsi(&w(lsi), jdim, kdim, idim);
                                FortranArray2DRef<double> wt2wk(&wk(lt2wk), jdim * kdim, 6);
                                FortranArray3DRef<double> wt3wk(&wk(lt3wk), jdim * kdim, idim, 5);
                                FortranArray2DRef<double> wvel(&wk(lvel), jdim*kdim*idim, 3);
                                FortranArray2DRef<double> wt1wk(&wk(lt1wk), jdim*kdim*idim, 3);
                                FortranArray3DRef<double> wvol(&w(lvol), jdim, kdim, idim);
                                FortranArray2DRef<double> wt4wk(&wk(lt4wk), jdim, kdim);
                                iflag = -1;
                                metric_ns::metric(jdim, kdim, idim,
                                    wx, wy, wz, wsj, wsk, wsi,
                                    wt2wk, wt3wk, nbl, iflag, icnt,
                                    nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
                                    ibcinfo, jbcinfo, kbcinfo, maxbl, maxseg,
                                    nou, bou, nbuf, ibufdim, myid, mblk2nd);
                                tmetric_ns::tmetric(jdim, kdim, idim,
                                    wsj, wsk, wsi,
                                    wx, wy, wz,
                                    wvel, wt1wk, wt2wk, wt3wk, nbl);
                                iflagv = 0;
                                int iflagv1 = 1;
                                cellvol_ns::cellvol(jdim, kdim, idim,
                                    wx, wy, wz, wsj, wsk, wsi,
                                    wvol, wt4wk, nou, bou, nbuf, ibufdim,
                                    myid, mblk2nd, maxbl, nbl,
                                    iflagv1, iflagv, imin, imax, jmin, jmax, kmin, kmax);
                            }
                        }

                        // coarser levels
                        ncg = ncgg(igridg(nbl)) - (mseq - iseqr);
                        if (ncg > 0 && mgflag > 0) {
                            nbll = nbl;
                            for (m = 1; m <= ncg; m++) {
                                nbll = nbll + 1;
                                nbllm1 = nbll - 1;
                                for (is = 1; is <= nsegdfrm(nbl); is++) {
                                    xorgae(nbll,is) = xorgae(nbl,is);
                                    yorgae(nbll,is) = yorgae(nbl,is);
                                    zorgae(nbll,is) = zorgae(nbl,is);
                                    thtxae(nbll,is) = thtxae(nbl,is);
                                    thtyae(nbll,is) = thtyae(nbl,is);
                                    thtzae(nbll,is) = thtzae(nbl,is);
                                }
                                lvolc  = lw(8,  nbll);
                                lxc    = lw(10, nbll);
                                lyc    = lw(11, nbll);
                                lzc    = lw(12, nbll);
                                lxtbjc = lw(36, nbll);
                                lxtbkc = lw(37, nbll);
                                lxtbic = lw(38, nbll);
                                latbjc = lw(39, nbll);
                                latbkc = lw(40, nbll);
                                latbic = lw(41, nbll);
                                lvelc  = lt1wk;

                                if (isklton == 1) {
                                    nou(1) = std::min(nou(1)+1, ibufdim);
                                    std::snprintf(bou(nou(1),1), 120,
                                        "   creating coarser block%4d of dimensions (I/J/K) :%4d%4d%4d",
                                        nbll, ii2, jj2, kk2);
                                }

                                {
                                    FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                                    FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                                    FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                                    FortranArray3DRef<double> wxc(&w(lxc), jj2, kk2, ii2);
                                    FortranArray3DRef<double> wyc(&w(lyc), jj2, kk2, ii2);
                                    FortranArray3DRef<double> wzc(&w(lzc), jj2, kk2, ii2);
                                    collx_ns::collx(wx, wy, wz, wxc, wyc, wzc,
                                        jdim, kdim, idim, jj2, kk2, ii2);

                                    FortranArray4DRef<double> wvel(&wk(lvel), jdim, kdim, idim, 3);
                                    FortranArray4DRef<double> wvelc(&wk(lvelc), jj2, kk2, ii2, 3);
                                    collxt_ns::collxt(wvel, wvelc,
                                        jdim, kdim, idim, jj2, kk2, ii2,
                                        nbllm1, nou, bou, nbuf, ibufdim);
                                    nv = jj2*kk2*ii2*3;
                                    for (izz = 1; izz <= nv; izz++)
                                        wk(lvel+izz-1) = wk(lvelc+izz-1);

                                    FortranArray3DRef<double> wvol(&w(lvol), jdim, kdim, idim);
                                    FortranArray3DRef<double> wvolc(&w(lvolc), jj2, kk2, ii2);
                                    collv_ns::collv(wvol, wvolc,
                                        jdim, kdim, idim, jj2, kk2, ii2);

                                    FortranArray4DRef<double> wxtbi(&w(lxtbi), jdim, kdim, 1, 1);
                                    FortranArray4DRef<double> wxtbic(&w(lxtbic), jj2, kk2, 1, 1);
                                    FortranArray4DRef<double> watbi(&w(latbi), jdim, kdim, 1, 1);
                                    FortranArray4DRef<double> watbic(&w(latbic), jj2, kk2, 1, 1);
                                    collxtb_ns::collxtb(wxtbi, wxtbic,
                                        jdim, kdim, jj2, kk2, nbllm1);
                                    collxtb_ns::collxtb(watbi, watbic,
                                        jdim, kdim, jj2, kk2, nbllm1);

                                    int idimm1_tmp = idim-1;
                                    int ii2m1_tmp  = ii2-1;
                                    FortranArray4DRef<double> wxtbj(&w(lxtbj), kdim, idim-1, 1, 1);
                                    FortranArray4DRef<double> wxtbjc(&w(lxtbjc), kk2, ii2-1, 1, 1);
                                    FortranArray4DRef<double> watbj(&w(latbj), kdim, idim-1, 1, 1);
                                    FortranArray4DRef<double> watbjc(&w(latbjc), kk2, ii2-1, 1, 1);
                                    collxtb_ns::collxtb(wxtbj, wxtbjc,
                                        kdim, idimm1_tmp, kk2, ii2m1_tmp, nbllm1);
                                    collxtb_ns::collxtb(watbj, watbjc,
                                        kdim, idimm1_tmp, kk2, ii2m1_tmp, nbllm1);

                                    FortranArray4DRef<double> wxtbk(&w(lxtbk), jdim, idim-1, 1, 1);
                                    FortranArray4DRef<double> wxtbkc(&w(lxtbkc), jj2, ii2-1, 1, 1);
                                    FortranArray4DRef<double> watbk(&w(latbk), jdim, idim-1, 1, 1);
                                    FortranArray4DRef<double> watbkc(&w(latbkc), jj2, ii2-1, 1, 1);
                                    collxtb_ns::collxtb(wxtbk, wxtbkc,
                                        jdim, idimm1_tmp, jj2, ii2m1_tmp, nbllm1);
                                    collxtb_ns::collxtb(watbk, watbkc,
                                        jdim, idimm1_tmp, jj2, ii2m1_tmp, nbllm1);
                                }

                                lead_ns::lead(nbll, lw, lw2, maxbl);

                                lvel  = 1 + 9*nsurf;
                                lt1wk = jdim*kdim*idim*3 + lvel;
                                lt2wk = jdim*kdim*idim*3 + lt1wk;
                                lt3wk = jdim*kdim*6 + lt2wk;

                                {
                                    FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                                    FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                                    FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                                    FortranArray3DRef<double> wsj(&w(lsj), jdim, kdim, idim);
                                    FortranArray3DRef<double> wsk(&w(lsk), jdim, kdim, idim);
                                    FortranArray3DRef<double> wsi(&w(lsi), jdim, kdim, idim);
                                    FortranArray2DRef<double> wt2wk(&wk(lt2wk), jdim * kdim, 6);
                                    FortranArray3DRef<double> wt3wk(&wk(lt3wk), jdim * kdim, idim, 5);
                                    FortranArray2DRef<double> wvel(&wk(lvel), jdim*kdim*idim, 3);
                                    FortranArray2DRef<double> wt1wk(&wk(lt1wk), jdim*kdim*idim, 3);
                                    iflag = -1;
                                    metric_ns::metric(jdim, kdim, idim,
                                        wx, wy, wz, wsj, wsk, wsi,
                                        wt2wk, wt3wk, nbll, iflag, icnt,
                                        nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
                                        ibcinfo, jbcinfo, kbcinfo, maxbl, maxseg,
                                        nou, bou, nbuf, ibufdim, myid, mblk2nd);

                                    tmetric_ns::tmetric(jdim, kdim, idim,
                                        wsj, wsk, wsi,
                                        wx, wy, wz,
                                        wvel, wt1wk, wt2wk, wt3wk, nbll);
                                }
                            } // end do 1820

                            lead_ns::lead(nbl, lw, lw2, maxbl);
                        }
                    } // end if (iunst > 1)
                } // end if (idefrm(nbl) > 0)
            } // end if myid...

            int iunit_tmp2 = 11;
            writ_buf_ns::writ_buf(nbl, iunit_tmp2, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl);
        } // end do 200
    } // end if (iunst > 1 || idef_ss > 0)
}


// ============================================================
// elglobfe
// ============================================================
void elglobfe(FortranArray1DRef<double> sa, FortranArray1DRef<double> xst,
              FortranArray1DRef<double> yst, FortranArray1DRef<double> zst,
              FortranArray1DRef<double> xix, FortranArray1DRef<double> xiy,
              FortranArray1DRef<double> xiz, FortranArray1DRef<double> etax,
              FortranArray1DRef<double> etay, FortranArray1DRef<double> etaz,
              FortranArray1DRef<double> zetax, FortranArray1DRef<double> zetay,
              FortranArray1DRef<double> zetaz, FortranArray1DRef<double> ei,
              FortranArray1DRef<double> ej, FortranArray1DRef<double> ek,
              FortranArray1DRef<double> gij, FortranArray1DRef<double> gjk,
              FortranArray1DRef<double> gik, FortranArray1DRef<double> ooj,
              double& eps, FortranArray2DRef<double> stiffl,
              FortranArray3DRef<int> islavept, FortranArray1DRef<int> ija,
              int& nslave, int& nmaster, int& nnodes,
              int& myhost, int& myid, int& mycomm,
              FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> nblelst,
              int& maxbl, int& iseqr)
{
    int& i2d = cmn_twod.i2d;

    // local arrays
    double tinv[6][6], t[6][6], ul[6][6], c[6][6];
    double cjp1[6][6], cjm1[6][6], cim1[6][6], cip1[6][6];
    double ckm1[6][6], ckp1[6][6];
    // 0-based arrays
    int n1[8], n2[8], n11[8];
    int n11i[20][8], n33i[20][8];
    int iimax[8], n44i[8];

    int isl1, isl2, isl3, isl4, isl5, isl6, isl7;
    int n, nbl, ii, jj, kk, ii3, ii4, j, j1, j2, j3, j2sta, j2end;
    int irow, icol, itst1, itst2, i1, n22i;
    double xi, yi, zi, xj, yj, zj, xk, yk, zk;
    double denomi, denomj;
    double t11, t12, t13, t21, t22, t23, t31, t32, t33;
    double sum, DTNRM, DETM;
    double c0, c1, c2, c3, c4, c5, c6;
    double term1, term2, term3, term4, term5, term6;
    double term7, term8, term9, term10, term11, term12;
    double term13, term14, term15, term16, term17, term18;
    double term19, term20, term21;
    double oj;

    isl1 = 8;
    isl2 = 2;
    isl3 = 3;
    isl4 = 4;
    isl5 = 5;
    if (i2d != 0) {
        isl6 = 4;
        isl7 = 5;
    } else {
        isl6 = 6;
        isl7 = 7;
    }
    n2[0] = 1;
    n2[1] = 4;
    n2[2] = 16;
    n2[3] = 13;
    n2[4] = 7;
    n2[5] = 19;
    n2[6] = 10;
    n2[7] = 22;

    for (n = 1; n <= nslave; n++) {
        n1[0] = n;
        n1[1] = islavept(n, 3, iseqr);
        n1[2] = islavept(n, 5, iseqr);
        n1[3] = islavept(n1[2], 3, iseqr);
        n1[4] = islavept(n, 7, iseqr);
        n1[5] = islavept(n1[2], 7, iseqr);
        n1[6] = islavept(n1[1], 7, iseqr);
        n1[7] = islavept(n1[6], 5, iseqr);
        nbl = islavept(n, 9, iseqr);
        if (n1[1] != n && n1[2] != n && n1[4] != n &&
            islavept(n1[1], 9, iseqr) == nbl &&
            islavept(n1[2], 9, iseqr) == nbl &&
            islavept(n1[4], 9, iseqr) == nbl) {

            xi = .25*(xst(n1[4])-xst(n1[0])+xst(n1[5])-xst(n1[2])
                     +xst(n1[6])-xst(n1[1])+xst(n1[7])-xst(n1[3]));
            zi = .25*(zst(n1[4])-zst(n1[0])+zst(n1[5])-zst(n1[2])
                     +zst(n1[6])-zst(n1[1])+zst(n1[7])-zst(n1[3]));
            yi = .25*(yst(n1[4])-yst(n1[0])+yst(n1[5])-yst(n1[2])
                     +yst(n1[6])-yst(n1[1])+yst(n1[7])-yst(n1[3]));
            xj = .25*(xst(n1[1])-xst(n1[0])+xst(n1[3])-xst(n1[2])
                     +xst(n1[6])-xst(n1[4])+xst(n1[7])-xst(n1[5]));
            yj = .25*(yst(n1[1])-yst(n1[0])+yst(n1[3])-yst(n1[2])
                     +yst(n1[6])-yst(n1[4])+yst(n1[7])-yst(n1[5]));
            zj = .25*(zst(n1[1])-zst(n1[0])+zst(n1[3])-zst(n1[2])
                     +zst(n1[6])-zst(n1[4])+zst(n1[7])-zst(n1[5]));
            denomi = std::sqrt(xi*xi+yi*yi+zi*zi);
            xi = xi/denomi; yi = yi/denomi; zi = zi/denomi;
            denomj = std::sqrt(xj*xj+yj*yj+zj*zj);
            xj = xj/denomj; yj = yj/denomj; zj = zj/denomj;
            xk = yi*zj - zi*yj;
            yk = zi*xj - xi*zj;
            zk = xi*yj - yi*xj;
            xi = yj*zk - zj*yk;
            yi = zj*xk - xj*zk;
            zi = xj*yk - yj*xk;

            t11 = xi; t21 = yi; t31 = zi;
            t12 = xj; t22 = yj; t32 = zj;
            t13 = xk; t23 = yk; t33 = zk;
            t[0][0] = t11*t11; t[0][1] = t12*t12; t[0][2] = t13*t13;
            t[0][3] = 2.*t11*t12; t[0][4] = 2.*t12*t13; t[0][5] = 2.*t11*t13;
            t[1][0] = t21*t21; t[1][1] = t22*t22; t[1][2] = t23*t23;
            t[1][3] = 2.*t21*t22; t[1][4] = 2.*t22*t23; t[1][5] = 2.*t21*t23;
            t[2][0] = t31*t31; t[2][1] = t32*t32; t[2][2] = t33*t33;
            t[2][3] = 2.*t31*t32; t[2][4] = 2.*t32*t33; t[2][5] = 2.*t31*t33;
            t[3][0] = t11*t21; t[3][1] = t12*t22; t[3][2] = t13*t23;
            t[3][3] = t11*t22+t21*t12; t[3][4] = t12*t23+t13*t22; t[3][5] = t11*t23+t13*t21;
            t[4][0] = t31*t21; t[4][1] = t22*t32; t[4][2] = t23*t33;
            t[4][3] = t21*t32+t22*t31; t[4][4] = t22*t33+t32*t23; t[4][5] = t21*t33+t23*t31;
            t[5][0] = t11*t31; t[5][1] = t12*t32; t[5][2] = t13*t33;
            t[5][3] = t11*t32+t12*t31; t[5][4] = t12*t33+t13*t32; t[5][5] = t11*t33+t31*t13;

            // Build FortranArray2D wrappers for INVDET call
            // INVDET(ul, t, 6, DTNRM, DETM)
            // ul and t are local 6x6 C arrays; wrap as FortranArray2D
            {
                FortranArray2D<double> t_fa(6,6), ul_fa(6,6);
                for (int r=1;r<=6;r++) for (int cc=1;cc<=6;cc++) t_fa(r,cc)=t[r-1][cc-1];
                int n6=6;
                deform_ns::invdet(ul_fa, t_fa, n6, DTNRM, DETM);
                for (int r=1;r<=6;r++) for (int cc=1;cc<=6;cc++) ul[r-1][cc-1]=ul_fa(r,cc);
            }

            for (jj = 0; jj < 6; jj++) {
                tinv[0][jj] = ei(n)  * ul[0][jj];
                tinv[1][jj] = ej(n)  * ul[1][jj];
                tinv[2][jj] = ek(n)  * ul[2][jj];
                tinv[3][jj] = gij(n) * ul[3][jj];
                tinv[4][jj] = gjk(n) * ul[4][jj];
                tinv[5][jj] = gik(n) * ul[5][jj];
            }
            for (jj = 0; jj < 6; jj++) {
                for (ii = 0; ii < 6; ii++) {
                    sum = 0.;
                    for (kk = 0; kk < 6; kk++)
                        sum += t[ii][kk] * tinv[kk][jj];
                    c[ii][jj] = sum;
                }
            }

            c0 = 1./2.; c1 = 1./6.; c2 = 1./9.; c3 = 1./12.;
            c4 = 1./18.; c5 = 1./24.; c6 = 1./36.;

// Helper macro: c(i,j) in Fortran 1-based → c[i-1][j-1] in C
#define C(i,j) c[(i)-1][(j)-1]
// stiffl(i,j) in Fortran 1-based → stiffl(i,j) in C++ (FortranArray2D)

            // Rows 1-3
            term1 = xix(n)*xix(n)+etax(n)*etax(n)+zetax(n)*zetax(n);
            term2 = etax(n)*xix(n)+zetax(n)*xix(n)+zetax(n)*etax(n);
            term3 = xiy(n)*xiy(n)+etay(n)*etay(n)+zetay(n)*zetay(n);
            term4 = etay(n)*xiy(n)+zetay(n)*xiy(n)+zetay(n)*etay(n);
            term5 = xiz(n)*xiz(n)+etaz(n)*etaz(n)+zetaz(n)*zetaz(n);
            term6 = etaz(n)*xiz(n)+zetaz(n)*xiz(n)+zetaz(n)*etaz(n);
            term7 = xix(n)*xiy(n)+etay(n)*etax(n)+zetay(n)*zetax(n);
            term8 = xix(n)*etay(n)+xix(n)*zetay(n)+xiy(n)*etax(n)
                   +zetay(n)*etax(n)+xiy(n)*zetax(n)+etay(n)*zetax(n);
            term9 = xiz(n)*xix(n)+etaz(n)*etax(n)+zetaz(n)*zetax(n);
            term10= xix(n)*etaz(n)+xix(n)*zetaz(n)+xiz(n)*etax(n)
                   +zetaz(n)*etax(n)+xiz(n)*zetax(n)+etaz(n)*zetax(n);
            term11= xiz(n)*xiy(n)+etaz(n)*etay(n)+zetaz(n)*zetay(n);
            term12= xiy(n)*etaz(n)+xiy(n)*zetaz(n)+xiz(n)*etay(n)
                   +zetaz(n)*etay(n)+xiz(n)*zetay(n)+etaz(n)*zetay(n);
            term13 = c2*term1+c1*term2;
            term14 = c2*term3+c1*term4;
            term15 = c2*term5+c1*term6;
            term16 = c2*term7+c3*term8;
            term17 = c2*term9+c3*term10;
            term18 = c2*term11+c3*term12;
            stiffl(1,1)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15
                       +(C(1,4)+C(4,1))*term16+(C(1,6)+C(6,1))*term17+(C(4,6)+C(6,4))*term18;
            stiffl(1,2)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15
                       +(C(1,2)+C(4,4))*term16+(C(1,5)+C(6,4))*term17+(C(4,5)+C(6,2))*term18;
            stiffl(1,3)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15
                       +(C(1,5)+C(4,6))*term16+(C(1,3)+C(6,6))*term17+(C(4,3)+C(6,5))*term18;
            stiffl(2,2)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15
                       +(C(2,4)+C(4,2))*term16+(C(2,5)+C(5,2))*term18+(C(4,5)+C(5,4))*term17;
            stiffl(2,3)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15
                       +(C(2,6)+C(4,5))*term16+(C(2,3)+C(5,5))*term18+(C(4,3)+C(5,6))*term17;
            stiffl(3,3)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15
                       +(C(6,5)+C(5,6))*term16+(C(3,5)+C(5,3))*term18+(C(3,6)+C(6,3))*term17;

            term13= -c2*xix(n)*xix(n)+c4*etax(n)*etax(n)+c3*zetax(n)*etax(n)+c4*zetax(n)*zetax(n);
            term14= -c2*xiy(n)*xiy(n)+c4*etay(n)*etay(n)+c3*zetay(n)*etay(n)+c4*zetay(n)*zetay(n);
            term15= -c2*xiz(n)*xiz(n)+c4*etaz(n)*etaz(n)+c3*zetaz(n)*etaz(n)+c4*zetaz(n)*zetaz(n);
            term16= -c2*xiy(n)*xix(n)+c4*etay(n)*etax(n)+c5*zetay(n)*etax(n)+c5*etay(n)*zetax(n)+c4*zetay(n)*zetax(n);
            term17= c3*(etax(n)*xiy(n)+zetax(n)*xiy(n)-xix(n)*etay(n)-xix(n)*zetay(n));
            term18= -c2*xiz(n)*xix(n)+c4*etaz(n)*etax(n)+c5*zetaz(n)*etax(n)+c5*etaz(n)*zetax(n)+c4*zetaz(n)*zetax(n);
            term19= c3*(etax(n)*xiz(n)+zetax(n)*xiz(n)-xix(n)*etaz(n)-xix(n)*zetaz(n));
            term20= -c2*xiz(n)*xiy(n)+c4*etaz(n)*etay(n)+c5*zetaz(n)*etay(n)+c5*etaz(n)*zetay(n)+c4*zetaz(n)*zetay(n);
            term21= c3*(etay(n)*xiz(n)+zetay(n)*xiz(n)-xiy(n)*etaz(n)-xiy(n)*zetaz(n));
            stiffl(1,4)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(4,1)-C(1,4))*term17+(C(1,6)+C(6,1))*term18+(C(6,1)-C(1,6))*term19+(C(4,6)+C(6,4))*term20+(C(6,4)-C(4,6))*term21;
            stiffl(2,4)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(2,1)-C(4,4))*term17+(C(5,1)+C(4,6))*term18+(C(5,1)-C(4,6))*term19+(C(5,4)+C(2,6))*term20+(C(5,4)-C(2,6))*term21;
            stiffl(3,4)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(5,1)-C(6,4))*term17+(C(3,1)+C(6,6))*term18+(C(3,1)-C(6,6))*term19+(C(3,4)+C(5,6))*term20+(C(3,4)-C(5,6))*term21;
            stiffl(3,5)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(5,4)-C(6,2))*term17+(C(3,4)+C(6,5))*term18+(C(3,4)-C(6,5))*term19+(C(3,2)+C(5,5))*term20+(C(3,2)-C(5,5))*term21;
            stiffl(2,5)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(2,4)-C(4,2))*term17+(C(5,4)+C(4,5))*term18+(C(5,4)-C(4,5))*term19+(C(5,2)+C(2,5))*term20+(C(5,2)-C(2,5))*term21;
            stiffl(1,5)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(4,4)-C(1,2))*term17+(C(6,4)+C(1,5))*term18+(C(6,4)-C(1,5))*term19+(C(6,2)+C(4,5))*term20+(C(6,2)-C(4,5))*term21;
            stiffl(3,6)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(5,6)-C(6,5))*term17+(C(3,6)+C(6,3))*term18+(C(3,6)-C(6,3))*term19+(C(3,5)+C(5,3))*term20+(C(3,5)-C(5,3))*term21;
            stiffl(2,6)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(2,6)-C(4,5))*term17+(C(5,6)+C(4,3))*term18+(C(5,6)-C(4,3))*term19+(C(5,5)+C(2,3))*term20+(C(5,5)-C(2,3))*term21;
            stiffl(1,6)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(4,6)-C(1,5))*term17+(C(6,6)+C(1,3))*term18+(C(6,6)-C(1,3))*term19+(C(6,5)+C(4,3))*term20+(C(6,5)-C(4,3))*term21;

            term13=  c4*xix(n)*xix(n)+c4*etax(n)*etax(n)+c3*xix(n)*etax(n)-c2*zetax(n)*zetax(n);
            term14=  c4*xiy(n)*xiy(n)+c4*etay(n)*etay(n)+c3*xiy(n)*etay(n)-c2*zetay(n)*zetay(n);
            term15=  c4*xiz(n)*xiz(n)+c4*etaz(n)*etaz(n)+c3*xiz(n)*etaz(n)-c2*zetaz(n)*zetaz(n);
            term16=  c4*xiy(n)*xix(n)+c5*etay(n)*xix(n)+c5*xiy(n)*etax(n)+c4*etay(n)*etax(n)-c2*zetay(n)*zetax(n);
            term17= c3*(zetax(n)*xiy(n)+zetax(n)*etay(n)-xix(n)*zetay(n)-etax(n)*zetay(n));
            term18=  c4*xiz(n)*xix(n)+c5*etaz(n)*xix(n)+c5*xiz(n)*etax(n)+c4*etaz(n)*etax(n)-c2*zetaz(n)*zetax(n);
            term19= c3*(zetax(n)*xiz(n)+zetax(n)*etaz(n)-xix(n)*zetaz(n)-etax(n)*zetaz(n));
            term20=  c4*xiz(n)*xiy(n)+c5*etaz(n)*xiy(n)+c5*xiz(n)*etay(n)+c4*etaz(n)*etay(n)-c2*zetaz(n)*zetay(n);
            term21= c3*(zetay(n)*xiz(n)+zetay(n)*etaz(n)-xiy(n)*zetaz(n)-etay(n)*zetaz(n));
            stiffl(1,7)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(2,7)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(3,7)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(3,8)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(2,8)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(1,8)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(3,9)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(2,9)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(1,9)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13= -c4*xix(n)*xix(n)+c6*etax(n)*etax(n)-c1*xix(n)*zetax(n)-c4*zetax(n)*zetax(n);
            term14= -c4*xiy(n)*xiy(n)+c6*etay(n)*etay(n)-c1*xiy(n)*zetay(n)-c4*zetay(n)*zetay(n);
            term15= -c4*xiz(n)*xiz(n)+c6*etaz(n)*etaz(n)-c1*xiz(n)*zetaz(n)-c4*zetaz(n)*zetaz(n);
            term16= -c4*xiy(n)*xix(n)-c3*zetay(n)*xix(n)+c6*etay(n)*etax(n)-c3*xiy(n)*zetax(n)-c4*zetay(n)*zetax(n);
            term17= c5*(-etax(n)*xiy(n)+xix(n)*etay(n)+zetax(n)*etay(n)-etax(n)*zetay(n));
            term18= -c4*xiz(n)*xix(n)-c3*zetaz(n)*xix(n)+c6*etaz(n)*etax(n)-c3*xiz(n)*zetax(n)-c4*zetaz(n)*zetax(n);
            term19= c5*(-etax(n)*xiz(n)+xix(n)*etaz(n)+zetax(n)*etaz(n)-etax(n)*zetaz(n));
            term20= -c4*xiz(n)*xiy(n)-c3*zetaz(n)*xiy(n)+c6*etaz(n)*etay(n)-c3*xiz(n)*zetay(n)-c4*zetaz(n)*zetay(n);
            term21= c5*(-etay(n)*xiz(n)+xiy(n)*etaz(n)+zetay(n)*etaz(n)-etay(n)*zetaz(n));
            stiffl(1,10)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(2,10)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(3,10)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(3,11)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(2,11)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(1,11)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(3,12)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(2,12)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(1,12)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13= -c4*xix(n)*xix(n)-c4*etax(n)*etax(n)-c1*xix(n)*etax(n)+c6*zetax(n)*zetax(n);
            term14= -c4*xiy(n)*xiy(n)-c4*etay(n)*etay(n)-c1*xiy(n)*etay(n)+c6*zetay(n)*zetay(n);
            term15= -c4*xiz(n)*xiz(n)-c4*etaz(n)*etaz(n)-c1*xiz(n)*etaz(n)+c6*zetaz(n)*zetaz(n);
            term16= -c4*xiy(n)*xix(n)-c3*etay(n)*xix(n)-c3*xiy(n)*etax(n)-c4*etay(n)*etax(n)+c6*zetay(n)*zetax(n);
            term17= c5*(-zetax(n)*xiy(n)-zetax(n)*etay(n)+zetay(n)*xix(n)+etax(n)*zetay(n));
            term18= -c4*xiz(n)*xix(n)-c3*etaz(n)*xix(n)-c3*xiz(n)*etax(n)-c4*etaz(n)*etax(n)+c6*zetaz(n)*zetax(n);
            term19= c5*(-zetax(n)*xiz(n)-zetax(n)*etaz(n)+zetaz(n)*xix(n)+etax(n)*zetaz(n));
            term20= -c4*xiz(n)*xiy(n)-c3*etaz(n)*xiy(n)-c3*xiz(n)*etay(n)-c4*etaz(n)*etay(n)+c6*zetaz(n)*zetay(n);
            term21= c5*(-zetay(n)*xiz(n)-zetay(n)*etaz(n)+zetaz(n)*xiy(n)+etay(n)*zetaz(n));
            stiffl(1,13)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(2,13)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(3,13)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(3,14)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(2,14)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(1,14)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(3,15)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(2,15)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(1,15)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13=  c4*xix(n)*xix(n)-c2*etax(n)*etax(n)+c3*xix(n)*zetax(n)+c4*zetax(n)*zetax(n);
            term14=  c4*xiy(n)*xiy(n)-c2*etay(n)*etay(n)+c3*xiy(n)*zetay(n)+c4*zetay(n)*zetay(n);
            term15=  c4*xiz(n)*xiz(n)-c2*etaz(n)*etaz(n)+c3*xiz(n)*zetaz(n)+c4*zetaz(n)*zetaz(n);
            term16=  c4*xiy(n)*xix(n)+c5*zetay(n)*xix(n)+c5*xiy(n)*zetax(n)-c2*etay(n)*etax(n)+c4*zetay(n)*zetax(n);
            term17= c3*(-etay(n)*xix(n)+etax(n)*xiy(n)+etax(n)*zetay(n)-zetax(n)*etay(n));
            term18=  c4*xiz(n)*xix(n)+c5*zetaz(n)*xix(n)+c5*xiz(n)*zetax(n)-c2*etaz(n)*etax(n)+c4*zetaz(n)*zetax(n);
            term19= c3*(-etaz(n)*xix(n)+etax(n)*xiz(n)+etax(n)*zetaz(n)-zetax(n)*etaz(n));
            term20=  c4*xiz(n)*xiy(n)+c5*zetaz(n)*xiy(n)+c5*xiz(n)*zetay(n)-c2*etaz(n)*etay(n)+c4*zetaz(n)*zetay(n);
            term21= c3*(-etaz(n)*xiy(n)+etay(n)*xiz(n)+etay(n)*zetaz(n)-zetay(n)*etaz(n));
            stiffl(1,16)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(2,16)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(3,16)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(3,17)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(2,17)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(1,17)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(3,18)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(2,18)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(1,18)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13=  c6*xix(n)*xix(n)-c4*etax(n)*etax(n)-c1*zetax(n)*etax(n)-c4*zetax(n)*zetax(n);
            term14=  c6*xiy(n)*xiy(n)-c4*etay(n)*etay(n)-c1*zetay(n)*etay(n)-c4*zetay(n)*zetay(n);
            term15=  c6*xiz(n)*xiz(n)-c4*etaz(n)*etaz(n)-c1*zetaz(n)*etaz(n)-c4*zetaz(n)*zetaz(n);
            term16=  c6*xiy(n)*xix(n)-c4*etay(n)*etax(n)-c3*zetay(n)*etax(n)-c3*etay(n)*zetax(n)-c4*zetay(n)*zetax(n);
            term17= c5*(etax(n)*xiy(n)+zetax(n)*xiy(n)-xix(n)*etay(n)-xix(n)*zetay(n));
            term18=  c6*xiz(n)*xix(n)-c4*etaz(n)*etax(n)-c3*zetaz(n)*etax(n)-c3*etaz(n)*zetax(n)-c4*zetaz(n)*zetax(n);
            term19= c5*(etax(n)*xiz(n)+zetax(n)*xiz(n)-xix(n)*etaz(n)-xix(n)*zetaz(n));
            term20=  c6*xiz(n)*xiy(n)-c4*etaz(n)*etay(n)-c3*zetaz(n)*etay(n)-c3*etaz(n)*zetay(n)-c4*zetaz(n)*zetay(n);
            term21= c5*(etay(n)*xiz(n)+zetay(n)*xiz(n)-xiy(n)*etaz(n)-xiy(n)*zetaz(n));
            stiffl(1,19)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(2,19)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(3,19)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(3,20)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(2,20)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(1,20)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(3,21)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(2,21)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(1,21)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13= -c6*xix(n)*xix(n)-c3*etax(n)*xix(n)-c3*zetax(n)*xix(n)-c6*etax(n)*etax(n)-c3*zetax(n)*etax(n)-c6*zetax(n)*zetax(n);
            term14= -c6*xiy(n)*xiy(n)-c3*etay(n)*xiy(n)-c3*zetay(n)*xiy(n)-c6*etay(n)*etay(n)-c3*zetay(n)*etay(n)-c6*zetay(n)*zetay(n);
            term15= -c6*xiz(n)*xiz(n)-c3*etaz(n)*xiz(n)-c3*zetaz(n)*xiz(n)-c6*etaz(n)*etaz(n)-c3*zetaz(n)*etaz(n)-c6*zetaz(n)*zetaz(n);
            term16= -c6*xiy(n)*xix(n)-c6*etay(n)*etax(n)-c6*zetay(n)*zetax(n);
            term17= c5*(-etay(n)*xix(n)-zetay(n)*xix(n)-xiy(n)*etax(n)-zetay(n)*etax(n)-xiy(n)*zetax(n)-etay(n)*zetax(n));
            term18= -c6*xiz(n)*xix(n)-c6*etaz(n)*etax(n)-c6*zetaz(n)*zetax(n);
            term19= c5*(-etaz(n)*xix(n)-zetaz(n)*xix(n)-xiz(n)*etax(n)-zetaz(n)*etax(n)-xiz(n)*zetax(n)-etaz(n)*zetax(n));
            term20= -c6*xiz(n)*xiy(n)-c6*etaz(n)*etay(n)-c6*zetaz(n)*zetay(n);
            term21= c5*(-etaz(n)*xiy(n)-zetaz(n)*xiy(n)-xiz(n)*etay(n)-zetaz(n)*etay(n)-xiz(n)*zetay(n)-etaz(n)*zetay(n));
            stiffl(1,22)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)+C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)+C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)+C(6,4))*term21;
            stiffl(2,22)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)+C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)+C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)+C(5,4))*term21;
            stiffl(3,22)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)+C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)+C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)+C(3,4))*term21;
            stiffl(3,23)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)+C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)+C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)+C(3,2))*term21;
            stiffl(2,23)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)+C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)+C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)+C(5,2))*term21;
            stiffl(1,23)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)+C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)+C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)+C(6,2))*term21;
            stiffl(3,24)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)+C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)+C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)+C(3,5))*term21;
            stiffl(2,24)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)+C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)+C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)+C(5,5))*term21;
            stiffl(1,24)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)+C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)+C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)+C(6,5))*term21;
#undef C

#define C(i,j) c[(i)-1][(j)-1]
            // Rows 4-6
            term1 = xix(n)*xix(n)+etax(n)*etax(n)+zetax(n)*zetax(n);
            term2 = -etax(n)*xix(n)-zetax(n)*xix(n)+zetax(n)*etax(n);
            term3 = xiy(n)*xiy(n)+etay(n)*etay(n)+zetay(n)*zetay(n);
            term4 = -etay(n)*xiy(n)-zetay(n)*xiy(n)+zetay(n)*etay(n);
            term5 = xiz(n)*xiz(n)+etaz(n)*etaz(n)+zetaz(n)*zetaz(n);
            term6 = -etaz(n)*xiz(n)-zetaz(n)*xiz(n)+zetaz(n)*etaz(n);
            term7 = xix(n)*xiy(n)+etay(n)*etax(n)+zetay(n)*zetax(n);
            term8 = -xix(n)*etay(n)-xix(n)*zetay(n)-xiy(n)*etax(n)+zetay(n)*etax(n)-xiy(n)*zetax(n)+etay(n)*zetax(n);
            term9 = xiz(n)*xix(n)+etaz(n)*etax(n)+zetaz(n)*zetax(n);
            term10= -xix(n)*etaz(n)-xix(n)*zetaz(n)-xiz(n)*etax(n)+zetaz(n)*etax(n)-xiz(n)*zetax(n)+etaz(n)*zetax(n);
            term11= xiz(n)*xiy(n)+etaz(n)*etay(n)+zetaz(n)*zetay(n);
            term12= -xiz(n)*etay(n)-xiz(n)*zetay(n)-xiy(n)*etaz(n)+zetay(n)*etaz(n)-xiy(n)*zetaz(n)+etay(n)*zetaz(n);
            term13 = c2*term1+c1*term2;
            term14 = c2*term3+c1*term4;
            term15 = c2*term5+c1*term6;
            term16 = c2*term7+c3*term8;
            term17 = c2*term9+c3*term10;
            term18 = c2*term11+c3*term12;
            stiffl(4,4)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,6)+C(6,1))*term17+(C(4,6)+C(6,4))*term18;
            stiffl(4,5)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(1,2)+C(4,4))*term16+(C(1,5)+C(6,4))*term17+(C(4,5)+C(6,2))*term18;
            stiffl(4,6)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(1,5)+C(4,6))*term16+(C(1,3)+C(6,6))*term17+(C(4,3)+C(6,5))*term18;
            stiffl(5,5)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,5)+C(5,4))*term17+(C(2,5)+C(5,2))*term18;
            stiffl(5,6)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,3)+C(5,6))*term17+(C(2,3)+C(5,5))*term18;
            stiffl(6,6)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(6,5)+C(5,6))*term16+(C(3,6)+C(6,3))*term17+(C(3,5)+C(5,3))*term18;

            term13= -c4*xix(n)*xix(n)+c6*etax(n)*etax(n)+c1*xix(n)*zetax(n)-c4*zetax(n)*zetax(n);
            term14= -c4*xiy(n)*xiy(n)+c6*etay(n)*etay(n)+c1*xiy(n)*zetay(n)-c4*zetay(n)*zetay(n);
            term15= -c4*xiz(n)*xiz(n)+c6*etaz(n)*etaz(n)+c1*xiz(n)*zetaz(n)-c4*zetaz(n)*zetaz(n);
            term16= -c4*xiy(n)*xix(n)+c3*zetay(n)*xix(n)+c3*xiy(n)*zetax(n)+c6*etay(n)*etax(n)-c4*zetay(n)*zetax(n);
            term17= c5*(etax(n)*xiy(n)-xix(n)*etay(n)+zetax(n)*etay(n)-etax(n)*zetay(n));
            term18= -c4*xiz(n)*xix(n)+c3*zetaz(n)*xix(n)+c3*xiz(n)*zetax(n)+c6*etaz(n)*etax(n)-c4*zetaz(n)*zetax(n);
            term19= c5*(etax(n)*xiz(n)-xix(n)*etaz(n)+zetax(n)*etaz(n)-etax(n)*zetaz(n));
            term20= -c4*xiz(n)*xiy(n)+c3*zetaz(n)*xiy(n)+c3*xiz(n)*zetay(n)+c6*etaz(n)*etay(n)-c4*zetaz(n)*zetay(n);
            term21= c5*(etay(n)*xiz(n)-xiy(n)*etaz(n)+zetay(n)*etaz(n)-etay(n)*zetaz(n));
            stiffl(4,7)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(5,7)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(6,7)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(6,8)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(5,8)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(4,8)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(6,9)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(5,9)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(4,9)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13=  c4*xix(n)*xix(n)+c4*etax(n)*etax(n)-c3*xix(n)*etax(n)-c2*zetax(n)*zetax(n);
            term14=  c4*xiy(n)*xiy(n)+c4*etay(n)*etay(n)-c3*xiy(n)*etay(n)-c2*zetay(n)*zetay(n);
            term15=  c4*xiz(n)*xiz(n)+c4*etaz(n)*etaz(n)-c3*xiz(n)*etaz(n)-c2*zetaz(n)*zetaz(n);
            term16=  c4*xiy(n)*xix(n)-c5*etax(n)*xiy(n)+c4*etay(n)*etax(n)-c5*xix(n)*etay(n)-c2*zetay(n)*zetax(n);
            term17= c3*(zetay(n)*xix(n)-zetay(n)*etax(n)-zetax(n)*xiy(n)+etay(n)*zetax(n));
            term18=  c4*xiz(n)*xix(n)-c5*etax(n)*xiz(n)+c4*etaz(n)*etax(n)-c5*xix(n)*etaz(n)-c2*zetaz(n)*zetax(n);
            term19= c3*(zetaz(n)*xix(n)-zetaz(n)*etax(n)-zetax(n)*xiz(n)+etaz(n)*zetax(n));
            term20=  c4*xiz(n)*xiy(n)-c5*etay(n)*xiz(n)+c4*etaz(n)*etay(n)-c5*xiy(n)*etaz(n)-c2*zetaz(n)*zetay(n);
            term21= c3*(zetaz(n)*xiy(n)-zetaz(n)*etay(n)-zetay(n)*xiz(n)+etaz(n)*zetay(n));
            stiffl(4,10)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(5,10)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(6,10)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(6,11)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(5,11)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(4,11)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(6,12)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(5,12)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(4,12)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13=  c4*xix(n)*xix(n)-c3*zetax(n)*xix(n)-c2*etax(n)*etax(n)+c4*zetax(n)*zetax(n);
            term14=  c4*xiy(n)*xiy(n)-c3*zetay(n)*xiy(n)-c2*etay(n)*etay(n)+c4*zetay(n)*zetay(n);
            term15=  c4*xiz(n)*xiz(n)-c3*zetaz(n)*xiz(n)-c2*etaz(n)*etaz(n)+c4*zetaz(n)*zetaz(n);
            term16=  c4*xiy(n)*xix(n)-c5*zetay(n)*xix(n)-c2*etay(n)*etax(n)-c5*xiy(n)*zetax(n)+c4*zetay(n)*zetax(n);
            term17= c3*(-etax(n)*xiy(n)+xix(n)*etay(n)-zetax(n)*etay(n)+etax(n)*zetay(n));
            term18=  c4*xiz(n)*xix(n)-c5*zetaz(n)*xix(n)-c2*etaz(n)*etax(n)-c5*xiz(n)*zetax(n)+c4*zetaz(n)*zetax(n);
            term19= c3*(-etax(n)*xiz(n)+xix(n)*etaz(n)-zetax(n)*etaz(n)+etax(n)*zetaz(n));
            term20=  c4*xiz(n)*xiy(n)-c5*zetaz(n)*xiy(n)-c2*etaz(n)*etay(n)-c5*xiz(n)*zetay(n)+c4*zetaz(n)*zetay(n);
            term21= c3*(-etay(n)*xiz(n)+xiy(n)*etaz(n)-zetay(n)*etaz(n)+etay(n)*zetaz(n));
            stiffl(4,13)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(5,13)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(6,13)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(6,14)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(5,14)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(4,14)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(6,15)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(5,15)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(4,15)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13= -c4*xix(n)*xix(n)+c1*etax(n)*xix(n)-c4*etax(n)*etax(n)+c6*zetax(n)*zetax(n);
            term14= -c4*xiy(n)*xiy(n)+c1*etay(n)*xiy(n)-c4*etay(n)*etay(n)+c6*zetay(n)*zetay(n);
            term15= -c4*xiz(n)*xiz(n)+c1*etaz(n)*xiz(n)-c4*etaz(n)*etaz(n)+c6*zetaz(n)*zetaz(n);
            term16= -c4*xiy(n)*xix(n)+c3*etay(n)*xix(n)+c3*xiy(n)*etax(n)-c4*etay(n)*etax(n)+c6*zetay(n)*zetax(n);
            term17= c5*(zetax(n)*xiy(n)-zetax(n)*etay(n)-xix(n)*zetay(n)+etax(n)*zetay(n));
            term18= -c4*xiz(n)*xix(n)+c3*etaz(n)*xix(n)+c3*xiz(n)*etax(n)-c4*etaz(n)*etax(n)+c6*zetaz(n)*zetax(n);
            term19= c5*(zetax(n)*xiz(n)-zetax(n)*etaz(n)-xix(n)*zetaz(n)+etax(n)*zetaz(n));
            term20= -c4*xiz(n)*xiy(n)+c3*etaz(n)*xiy(n)+c3*xiz(n)*etay(n)-c4*etaz(n)*etay(n)+c6*zetaz(n)*zetay(n);
            term21= c5*(zetay(n)*xiz(n)-zetay(n)*etaz(n)-xiy(n)*zetaz(n)+etay(n)*zetaz(n));
            stiffl(4,16)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(5,16)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(6,16)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(6,17)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(5,17)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(4,17)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(6,18)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(5,18)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(4,18)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13= -c6*(xix(n)*xix(n)+etax(n)*etax(n)+zetax(n)*zetax(n))+c3*(etax(n)*xix(n)+zetax(n)*xix(n)-zetax(n)*etax(n));
            term14= -c6*(xiy(n)*xiy(n)+etay(n)*etay(n)+zetay(n)*zetay(n))+c3*(etay(n)*xiy(n)+zetay(n)*xiy(n)-zetay(n)*etay(n));
            term15= -c6*(xiz(n)*xiz(n)+etaz(n)*etaz(n)+zetaz(n)*zetaz(n))+c3*(etaz(n)*xiz(n)+zetaz(n)*xiz(n)-zetaz(n)*etaz(n));
            term16= -c6*xiy(n)*xix(n)-c4*etay(n)*etax(n)-c6*zetay(n)*zetax(n)+c5*etay(n)*xix(n)+c5*zetay(n)*xix(n);
            term17= c5*(etax(n)*xiy(n)-etax(n)*zetay(n)+xiy(n)*zetax(n)-etay(n)*zetax(n));
            term18= -c6*xiz(n)*xix(n)-c4*etaz(n)*etax(n)-c6*zetaz(n)*zetax(n)+c5*etaz(n)*xix(n)+c5*zetaz(n)*xix(n);
            term19= c5*(etax(n)*xiz(n)-etax(n)*zetaz(n)+xiz(n)*zetax(n)-etaz(n)*zetax(n));
            term20= -c6*xiz(n)*xiy(n)-c4*etaz(n)*etay(n)-c6*zetaz(n)*zetay(n)+c5*etay(n)*xiz(n)+c5*zetay(n)*xiz(n);
            term21= c5*(etaz(n)*xiy(n)-etaz(n)*zetay(n)+xiy(n)*zetaz(n)-etay(n)*zetaz(n));
            stiffl(4,19)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)+C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)+C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)+C(6,4))*term21;
            stiffl(5,19)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)+C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)+C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)+C(5,4))*term21;
            stiffl(6,19)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)+C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)+C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)+C(3,4))*term21;
            stiffl(6,20)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)+C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)+C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)+C(3,2))*term21;
            stiffl(5,20)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)+C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)+C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)+C(5,2))*term21;
            stiffl(4,20)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)+C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)+C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)+C(6,2))*term21;
            stiffl(6,21)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)+C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)+C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)+C(3,5))*term21;
            stiffl(5,21)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)+C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)+C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)+C(5,5))*term21;
            stiffl(4,21)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)+C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)+C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)+C(6,5))*term21;

            term13=  c6*xix(n)*xix(n)-c4*etax(n)*etax(n)-c1*zetax(n)*etax(n)-c4*zetax(n)*zetax(n);
            term14=  c6*xiy(n)*xiy(n)-c4*etay(n)*etay(n)-c1*zetay(n)*etay(n)-c4*zetay(n)*zetay(n);
            term15=  c6*xiz(n)*xiz(n)-c4*etaz(n)*etaz(n)-c1*zetaz(n)*etaz(n)-c4*zetaz(n)*zetaz(n);
            term16=  c6*xiy(n)*xix(n)-c4*etay(n)*etax(n)-c3*zetay(n)*etax(n)-c3*etay(n)*zetax(n)-c4*zetay(n)*zetax(n);
            term17= c5*(-etax(n)*xiy(n)-zetax(n)*xiy(n)+xix(n)*etay(n)+zetay(n)*xix(n));
            term18=  c6*xiz(n)*xix(n)-c4*etaz(n)*etax(n)-c3*zetaz(n)*etax(n)-c3*etaz(n)*zetax(n)-c4*zetaz(n)*zetax(n);
            term19= c5*(-etax(n)*xiz(n)-zetax(n)*xiz(n)+xix(n)*etaz(n)+zetaz(n)*xix(n));
            term20=  c6*xiz(n)*xiy(n)-c4*etaz(n)*etay(n)-c3*zetaz(n)*etay(n)-c3*etaz(n)*zetay(n)-c4*zetaz(n)*zetay(n);
            term21= c5*(-etay(n)*xiz(n)-zetay(n)*xiz(n)+xiy(n)*etaz(n)+zetaz(n)*xiy(n));
            stiffl(4,22)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(5,22)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(6,22)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(6,23)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(5,23)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(4,23)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(6,24)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(5,24)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(4,24)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;
#undef C

#define C(i,j) c[(i)-1][(j)-1]
            // Rows 7-9
            term13=  c2*xix(n)*xix(n)+c2*etax(n)*etax(n)+c2*zetax(n)*zetax(n)+c1*etax(n)*xix(n)-c1*zetax(n)*xix(n)-c1*zetax(n)*etax(n);
            term14=  c2*xiy(n)*xiy(n)+c2*etay(n)*etay(n)+c2*zetay(n)*zetay(n)+c1*etay(n)*xiy(n)-c1*zetay(n)*xiy(n)-c1*zetay(n)*etay(n);
            term15=  c2*xiz(n)*xiz(n)+c2*etaz(n)*etaz(n)+c2*zetaz(n)*zetaz(n)+c1*etaz(n)*xiz(n)-c1*zetaz(n)*xiz(n)-c1*zetaz(n)*etaz(n);
            term16=  c2*xiy(n)*xix(n)+c2*etay(n)*etax(n)+c2*zetay(n)*zetax(n)+c3*etay(n)*xix(n)-c3*zetay(n)*xix(n)+c3*xiy(n)*etax(n)-c3*zetay(n)*etax(n)-c3*xiy(n)*zetax(n)-c3*etay(n)*zetax(n);
            term18=  c2*xiz(n)*xix(n)+c2*etaz(n)*etax(n)+c2*zetaz(n)*zetax(n)+c3*etaz(n)*xix(n)-c3*zetaz(n)*xix(n)+c3*xiz(n)*etax(n)-c3*zetaz(n)*etax(n)-c3*xiz(n)*zetax(n)-c3*etaz(n)*zetax(n);
            term20=  c2*xiz(n)*xiy(n)+c2*etaz(n)*etay(n)+c2*zetaz(n)*zetay(n)+c3*etaz(n)*xiy(n)-c3*zetaz(n)*xiy(n)+c3*xiz(n)*etay(n)-c3*zetaz(n)*etay(n)-c3*xiz(n)*zetay(n)-c3*etaz(n)*zetay(n);
            stiffl(7,7)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,6)+C(6,1))*term18+(C(4,6)+C(6,4))*term20;
            stiffl(8,8)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(5,4)+C(4,5))*term18+(C(5,2)+C(2,5))*term20;
            stiffl(7,8)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(6,4)+C(1,5))*term18+(C(6,2)+C(4,5))*term20;
            stiffl(9,9)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(3,6)+C(6,3))*term18+(C(3,5)+C(5,3))*term20;
            stiffl(8,9)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(5,6)+C(4,3))*term18+(C(5,5)+C(2,3))*term20;
            stiffl(7,9)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(6,6)+C(1,3))*term18+(C(6,5)+C(4,3))*term20;

            term13= -c2*xix(n)*xix(n)+c4*etax(n)*etax(n)-c3*zetax(n)*etax(n)+c4*zetax(n)*zetax(n);
            term14= -c2*xiy(n)*xiy(n)+c4*etay(n)*etay(n)-c3*zetay(n)*etay(n)+c4*zetay(n)*zetay(n);
            term15= -c2*xiz(n)*xiz(n)+c4*etaz(n)*etaz(n)-c3*zetaz(n)*etaz(n)+c4*zetaz(n)*zetaz(n);
            term16= -c2*xiy(n)*xix(n)+c4*etay(n)*etax(n)-c5*zetay(n)*etax(n)-c5*etay(n)*zetax(n)+c4*zetay(n)*zetax(n);
            term17= -c3*etax(n)*xiy(n)+c3*zetax(n)*xiy(n)+c3*xix(n)*etay(n)-c3*xix(n)*zetay(n);
            term18= -c2*xiz(n)*xix(n)+c4*etaz(n)*etax(n)-c5*zetaz(n)*etax(n)-c5*etaz(n)*zetax(n)+c4*zetaz(n)*zetax(n);
            term19= -c3*etax(n)*xiz(n)+c3*zetax(n)*xiz(n)+c3*xix(n)*etaz(n)-c3*xix(n)*zetaz(n);
            term20= -c2*xiz(n)*xiy(n)+c4*etaz(n)*etay(n)-c5*zetaz(n)*etay(n)-c5*etaz(n)*zetay(n)+c4*zetaz(n)*zetay(n);
            term21= -c3*etay(n)*xiz(n)+c3*zetay(n)*xiz(n)+c3*xiy(n)*etaz(n)-c3*xiy(n)*zetaz(n);
            stiffl(7,10)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(8,10)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(9,10)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(9,11)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(8,11)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(7,11)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(9,12)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(8,12)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(7,12)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;
#undef C

#define C(i,j) c[(i)-1][(j)-1]
            // Rows 7-9 continued (13-24)
            term13= -c6*xix(n)*xix(n)-c6*etax(n)*etax(n)-c6*zetax(n)*zetax(n)-c3*etax(n)*xix(n)+c3*zetax(n)*xix(n)+c3*zetax(n)*etax(n);
            term14= -c6*xiy(n)*xiy(n)-c6*etay(n)*etay(n)-c6*zetay(n)*zetay(n)-c3*etay(n)*xiy(n)+c3*zetay(n)*xiy(n)+c3*zetay(n)*etay(n);
            term15= -c6*xiz(n)*xiz(n)-c6*etaz(n)*etaz(n)-c6*zetaz(n)*zetaz(n)-c3*etaz(n)*xiz(n)+c3*zetaz(n)*xiz(n)+c3*zetaz(n)*etaz(n);
            term16= -c6*xiy(n)*xix(n)-c6*etay(n)*etax(n)-c6*zetay(n)*zetax(n)-c5*etay(n)*xix(n)+c5*zetay(n)*xix(n)-c5*xiy(n)*etax(n)+c5*zetay(n)*etax(n)+c5*xiy(n)*zetax(n)+c5*etay(n)*zetax(n);
            term18= -c6*xiz(n)*xix(n)-c6*etaz(n)*etax(n)-c6*zetaz(n)*zetax(n)-c5*etaz(n)*xix(n)+c5*zetaz(n)*xix(n)-c5*xiz(n)*etax(n)+c5*zetaz(n)*etax(n)+c5*xiz(n)*zetax(n)+c5*etaz(n)*zetax(n);
            term20= -c6*xiz(n)*xiy(n)-c6*etaz(n)*etay(n)-c6*zetaz(n)*zetay(n)-c5*etaz(n)*xiy(n)+c5*zetaz(n)*xiy(n)-c5*xiz(n)*etay(n)+c5*zetaz(n)*etay(n)+c5*xiz(n)*zetay(n)+c5*etaz(n)*zetay(n);
            stiffl(7,13)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,6)+C(6,1))*term18+(C(4,6)+C(6,4))*term20;
            stiffl(8,13)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(5,1)+C(4,6))*term18+(C(5,4)+C(2,6))*term20;
            stiffl(9,13)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(3,1)+C(6,6))*term18+(C(3,4)+C(5,6))*term20;
            stiffl(9,14)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(3,4)+C(6,5))*term18+(C(3,2)+C(5,5))*term20;
            stiffl(8,14)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(5,4)+C(4,5))*term18+(C(5,2)+C(2,5))*term20;
            stiffl(7,14)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(6,4)+C(1,5))*term18+(C(6,2)+C(4,5))*term20;
            stiffl(9,15)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(3,6)+C(6,3))*term18+(C(3,5)+C(5,3))*term20;
            stiffl(8,15)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(5,6)+C(4,3))*term18+(C(5,5)+C(2,3))*term20;
            stiffl(7,15)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(6,6)+C(1,3))*term18+(C(6,5)+C(4,3))*term20;

            term13=  c6*xix(n)*xix(n)-c4*etax(n)*etax(n)+c1*zetax(n)*etax(n)-c4*zetax(n)*zetax(n);
            term14=  c6*xiy(n)*xiy(n)-c4*etay(n)*etay(n)+c1*zetay(n)*etay(n)-c4*zetay(n)*zetay(n);
            term15=  c6*xiz(n)*xiz(n)-c4*etaz(n)*etaz(n)+c1*zetaz(n)*etaz(n)-c4*zetaz(n)*zetaz(n);
            term16=  c6*xiy(n)*xix(n)-c4*etay(n)*etax(n)+c3*zetay(n)*etax(n)+c3*etay(n)*zetax(n)-c4*zetay(n)*zetax(n);
            term17= -c5*etay(n)*xix(n)+c5*zetay(n)*xix(n)+c5*xiy(n)*etax(n)-c5*xiy(n)*zetax(n);
            term18=  c6*xiz(n)*xix(n)-c4*etaz(n)*etax(n)+c3*zetaz(n)*etax(n)+c3*etaz(n)*zetax(n)-c4*zetaz(n)*zetax(n);
            term19= -c5*etaz(n)*xix(n)+c5*zetaz(n)*xix(n)+c5*xiz(n)*etax(n)-c5*xiz(n)*zetax(n);
            term20=  c6*xiz(n)*xiy(n)-c4*etaz(n)*etay(n)+c3*zetaz(n)*etay(n)+c3*etaz(n)*zetay(n)-c4*zetaz(n)*zetay(n);
            term21= -c5*etaz(n)*xiy(n)+c5*zetaz(n)*xiy(n)+c5*xiz(n)*etay(n)-c5*xiz(n)*zetay(n);
            stiffl(7,16)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(8,16)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(9,16)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(9,17)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(8,17)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(7,17)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(9,18)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(8,18)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(7,18)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13=  c4*xix(n)*xix(n)-c3*zetax(n)*xix(n)-c2*etax(n)*etax(n)+c4*zetax(n)*zetax(n);
            term14=  c4*xiy(n)*xiy(n)-c3*zetay(n)*xiy(n)-c2*etay(n)*etay(n)+c4*zetay(n)*zetay(n);
            term15=  c4*xiz(n)*xiz(n)-c3*zetaz(n)*xiz(n)-c2*etaz(n)*etaz(n)+c4*zetaz(n)*zetaz(n);
            term16=  c4*xiy(n)*xix(n)-c5*zetay(n)*xix(n)-c2*etay(n)*etax(n)-c5*xiy(n)*zetax(n)+c4*zetay(n)*zetax(n);
            term17= c3*(xiy(n)*etax(n)-xix(n)*etay(n)+zetax(n)*etay(n)-etax(n)*zetay(n));
            term18=  c4*xiz(n)*xix(n)-c5*zetaz(n)*xix(n)-c2*etaz(n)*etax(n)-c5*xiz(n)*zetax(n)+c4*zetaz(n)*zetax(n);
            term19= c3*(xiz(n)*etax(n)-xix(n)*etaz(n)+zetax(n)*etaz(n)-etax(n)*zetaz(n));
            term20=  c4*xiz(n)*xiy(n)-c5*zetaz(n)*xiy(n)-c2*etaz(n)*etay(n)-c5*xiz(n)*zetay(n)+c4*zetaz(n)*zetay(n);
            term21= c3*(xiz(n)*etay(n)-xiy(n)*etaz(n)+zetay(n)*etaz(n)-etay(n)*zetaz(n));
            stiffl(7,19)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(8,19)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(9,19)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(9,20)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(8,20)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(7,20)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(9,21)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(8,21)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(7,21)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13= -c4*xix(n)*xix(n)-c4*etax(n)*etax(n)-c1*etax(n)*xix(n)+c6*zetax(n)*zetax(n);
            term14= -c4*xiy(n)*xiy(n)-c4*etay(n)*etay(n)-c1*etay(n)*xiy(n)+c6*zetay(n)*zetay(n);
            term15= -c4*xiz(n)*xiz(n)-c4*etaz(n)*etaz(n)-c1*etaz(n)*xiz(n)+c6*zetaz(n)*zetaz(n);
            term16= -c4*xiy(n)*xix(n)-c3*etay(n)*xix(n)-c3*xiy(n)*etax(n)-c4*etay(n)*etax(n)+c6*zetay(n)*zetax(n);
            term17= c5*(-xix(n)*zetay(n)-etax(n)*zetay(n)+zetax(n)*xiy(n)+zetax(n)*etay(n));
            term18= -c4*xiz(n)*xix(n)-c3*etaz(n)*xix(n)-c3*xiz(n)*etax(n)-c4*etaz(n)*etax(n)+c6*zetaz(n)*zetax(n);
            term19= c5*(-xix(n)*zetaz(n)-etax(n)*zetaz(n)+zetax(n)*xiz(n)+zetax(n)*etaz(n));
            term20= -c4*xiz(n)*xiy(n)-c3*etaz(n)*xiy(n)-c3*xiz(n)*etay(n)-c4*etaz(n)*etay(n)+c6*zetaz(n)*zetay(n);
            term21= c5*(-xiy(n)*zetaz(n)-etay(n)*zetaz(n)+zetay(n)*xiz(n)+zetay(n)*etaz(n));
            stiffl(7,22)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(8,22)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(9,22)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(9,23)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(8,23)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(7,23)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(9,24)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(8,24)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(7,24)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;
#undef C

#define C(i,j) c[(i)-1][(j)-1]
            // Rows 10-12
            term13=  c2*xix(n)*xix(n)+c2*etax(n)*etax(n)+c2*zetax(n)*zetax(n)-c1*etax(n)*xix(n)+c1*zetax(n)*xix(n)-c1*zetax(n)*etax(n);
            term14=  c2*xiy(n)*xiy(n)+c2*etay(n)*etay(n)+c2*zetay(n)*zetay(n)-c1*etay(n)*xiy(n)+c1*zetay(n)*xiy(n)-c1*zetay(n)*etay(n);
            term15=  c2*xiz(n)*xiz(n)+c2*etaz(n)*etaz(n)+c2*zetaz(n)*zetaz(n)-c1*etaz(n)*xiz(n)+c1*zetaz(n)*xiz(n)-c1*zetaz(n)*etaz(n);
            term16=  c2*xiy(n)*xix(n)+c2*etay(n)*etax(n)+c2*zetay(n)*zetax(n)-c3*etay(n)*xix(n)+c3*zetay(n)*xix(n)-c3*xiy(n)*etax(n)-c3*zetay(n)*etax(n)+c3*xiy(n)*zetax(n)-c3*etay(n)*zetax(n);
            term17 = 0.;
            term18=  c2*xiz(n)*xix(n)+c2*etaz(n)*etax(n)+c2*zetaz(n)*zetax(n)-c3*etaz(n)*xix(n)+c3*zetaz(n)*xix(n)-c3*xiz(n)*etax(n)-c3*zetaz(n)*etax(n)+c3*xiz(n)*zetax(n)-c3*etaz(n)*zetax(n);
            term19 = 0.;
            term20=  c2*xiz(n)*xiy(n)+c2*etaz(n)*etay(n)+c2*zetaz(n)*zetay(n)-c3*etaz(n)*xiy(n)+c3*zetaz(n)*xiy(n)-c3*xiz(n)*etay(n)-c3*zetaz(n)*etay(n)+c3*xiz(n)*zetay(n)-c3*etaz(n)*zetay(n);
            term21 = 0.;
            stiffl(10,10)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,6)+C(6,1))*term18+(C(4,6)+C(6,4))*term20;
            stiffl(11,11)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(5,4)+C(4,5))*term18+(C(5,2)+C(2,5))*term20;
            stiffl(10,11)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(6,4)+C(1,5))*term18+(C(6,2)+C(4,5))*term20;
            stiffl(12,12)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(3,6)+C(6,3))*term18+(C(3,5)+C(5,3))*term20;
            stiffl(11,12)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(5,6)+C(4,3))*term18+(C(5,5)+C(2,3))*term20;
            stiffl(10,12)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(6,6)+C(1,3))*term18+(C(6,5)+C(4,3))*term20;

            term13=  c6*xix(n)*xix(n)-c4*etax(n)*etax(n)-c4*zetax(n)*zetax(n)+c1*zetax(n)*etax(n);
            term14=  c6*xiy(n)*xiy(n)-c4*etay(n)*etay(n)-c4*zetay(n)*zetay(n)+c1*zetay(n)*etay(n);
            term15=  c6*xiz(n)*xiz(n)-c4*etaz(n)*etaz(n)-c4*zetaz(n)*zetaz(n)+c1*zetaz(n)*etaz(n);
            term16=  c6*xiy(n)*xix(n)-c4*etay(n)*etax(n)+c3*zetay(n)*etax(n)+c3*etay(n)*zetax(n)-c4*zetay(n)*zetax(n);
            term17= -c5*etax(n)*xiy(n)+c5*zetax(n)*xiy(n)+c5*xix(n)*etay(n)-c5*xix(n)*zetay(n);
            term18=  c6*xiz(n)*xix(n)-c4*etaz(n)*etax(n)+c3*zetaz(n)*etax(n)+c3*etaz(n)*zetax(n)-c4*zetaz(n)*zetax(n);
            term19= -c5*etax(n)*xiz(n)+c5*zetax(n)*xiz(n)+c5*xix(n)*etaz(n)-c5*xix(n)*zetaz(n);
            term20=  c6*xiz(n)*xiy(n)-c4*etaz(n)*etay(n)+c3*zetaz(n)*etay(n)+c3*etaz(n)*zetay(n)-c4*zetaz(n)*zetay(n);
            term21= -c5*etay(n)*xiz(n)+c5*zetay(n)*xiz(n)+c5*xiy(n)*etaz(n)-c5*xiy(n)*zetaz(n);
            stiffl(10,13)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(11,13)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(4,4)+C(2,1))*term16+(C(4,4)-C(2,1))*term17+(C(4,6)+C(5,1))*term18+(C(4,6)-C(5,1))*term19+(C(2,6)+C(5,4))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(12,13)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(6,4)+C(5,1))*term16+(C(6,4)-C(5,1))*term17+(C(6,6)+C(3,1))*term18+(C(6,6)-C(3,1))*term19+(C(5,6)+C(3,4))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(12,14)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(6,2)+C(5,4))*term16+(C(6,2)-C(5,4))*term17+(C(6,5)+C(3,4))*term18+(C(6,5)-C(3,4))*term19+(C(5,5)+C(3,2))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(11,14)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(4,2)+C(2,4))*term16+(C(4,2)-C(2,4))*term17+(C(4,5)+C(5,4))*term18+(C(4,5)-C(5,4))*term19+(C(2,5)+C(5,2))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(10,14)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(1,2)+C(4,4))*term16+(C(1,2)-C(4,4))*term17+(C(1,5)+C(6,4))*term18+(C(1,5)-C(6,4))*term19+(C(4,5)+C(6,2))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(12,15)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(6,5)+C(5,6))*term16+(C(6,5)-C(5,6))*term17+(C(6,3)+C(3,6))*term18+(C(6,3)-C(3,6))*term19+(C(5,3)+C(3,5))*term20+(C(5,3)-C(3,5))*term20;
            stiffl(11,15)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(4,5)+C(2,6))*term16+(C(4,5)-C(2,6))*term17+(C(4,3)+C(5,6))*term18+(C(4,3)-C(5,6))*term19+(C(2,3)+C(5,5))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(10,15)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(1,5)+C(4,6))*term16+(C(1,5)-C(4,6))*term17+(C(1,3)+C(6,6))*term18+(C(1,3)-C(6,6))*term19+(C(4,3)+C(6,5))*term20+(C(4,3)-C(6,5))*term21;

            term13= -c6*xix(n)*xix(n)+c3*etax(n)*xix(n)-c3*zetax(n)*xix(n)-c6*etax(n)*etax(n)+c3*zetax(n)*etax(n)-c6*zetax(n)*zetax(n);
            term14= -c6*xiy(n)*xiy(n)+c3*etay(n)*xiy(n)-c3*zetay(n)*xiy(n)-c6*etay(n)*etay(n)+c3*zetay(n)*etay(n)-c6*zetay(n)*zetay(n);
            term15= -c6*xiz(n)*xiz(n)+c3*etaz(n)*xiz(n)-c3*zetaz(n)*xiz(n)-c6*etaz(n)*etaz(n)+c3*zetaz(n)*etaz(n)-c6*zetaz(n)*zetaz(n);
            term16= -c6*xiy(n)*xix(n)+c5*etay(n)*xix(n)-c5*zetay(n)*xix(n)+c5*xiy(n)*etax(n)-c6*etay(n)*etax(n)+c5*zetay(n)*etax(n)-c5*xiy(n)*zetax(n)+c5*etay(n)*zetax(n)-c6*zetay(n)*zetax(n);
            term17 = 0.;
            term18= -c6*xiz(n)*xix(n)+c5*etaz(n)*xix(n)-c5*zetaz(n)*xix(n)+c5*xiz(n)*etax(n)-c6*etaz(n)*etax(n)+c5*zetaz(n)*etax(n)-c5*xiz(n)*zetax(n)+c5*etaz(n)*zetax(n)-c6*zetaz(n)*zetax(n);
            term19 = 0.;
            term20= -c6*xiz(n)*xiy(n)+c5*etaz(n)*xiy(n)-c5*zetaz(n)*xiy(n)+c5*xiz(n)*etay(n)-c6*etaz(n)*etay(n)+c5*zetaz(n)*etay(n)-c5*xiz(n)*zetay(n)+c5*etaz(n)*zetay(n)-c6*zetaz(n)*zetay(n);
            term21 = 0.;
            stiffl(10,16)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,6)+C(6,1))*term18+(C(4,6)+C(6,4))*term20;
            stiffl(11,16)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(5,1)+C(4,6))*term18+(C(5,4)+C(2,6))*term20;
            stiffl(12,16)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(3,1)+C(6,6))*term18+(C(3,4)+C(5,6))*term20;
            stiffl(12,17)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(3,4)+C(6,5))*term18+(C(3,2)+C(5,5))*term20;
            stiffl(11,17)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(5,4)+C(4,5))*term18+(C(5,2)+C(2,5))*term20;
            stiffl(10,17)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(6,4)+C(1,5))*term18+(C(6,2)+C(4,5))*term20;
            stiffl(12,18)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(3,6)+C(6,3))*term18+(C(3,5)+C(5,3))*term20;
            stiffl(11,18)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(5,6)+C(4,3))*term18+(C(5,5)+C(2,3))*term20;
            stiffl(10,18)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(6,6)+C(1,3))*term18+(C(6,5)+C(4,3))*term20;

            term13= -c4*xix(n)*xix(n)+c1*etax(n)*xix(n)-c4*etax(n)*etax(n)+c6*zetax(n)*zetax(n);
            term14= -c4*xiy(n)*xiy(n)+c1*etay(n)*xiy(n)-c4*etay(n)*etay(n)+c6*zetay(n)*zetay(n);
            term15= -c4*xiz(n)*xiz(n)+c1*etaz(n)*xiz(n)-c4*etaz(n)*etaz(n)+c6*zetaz(n)*zetaz(n);
            term16= -c4*xiy(n)*xix(n)+c3*etay(n)*xix(n)+c3*xiy(n)*etax(n)-c4*etay(n)*etax(n)+c6*zetay(n)*zetax(n);
            term17= c5*(-zetax(n)*xiy(n)+zetax(n)*etay(n)+xix(n)*zetay(n)-etax(n)*zetay(n));
            term18= -c4*xiz(n)*xix(n)+c3*etaz(n)*xix(n)+c3*xiz(n)*etax(n)-c4*etaz(n)*etax(n)+c6*zetaz(n)*zetax(n);
            term19= c5*(-zetax(n)*xiz(n)+zetax(n)*etaz(n)+xix(n)*zetaz(n)-etax(n)*zetaz(n));
            term20= -c4*xiz(n)*xiy(n)+c3*etaz(n)*xiy(n)+c3*xiz(n)*etay(n)-c4*etaz(n)*etay(n)+c6*zetaz(n)*zetay(n);
            term21= c5*(-zetay(n)*xiz(n)+zetay(n)*etaz(n)+xiy(n)*zetaz(n)-etay(n)*zetaz(n));
            stiffl(10,19)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(11,19)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(12,19)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(12,20)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(11,20)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(10,20)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(12,21)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(11,21)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(10,21)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13=  c4*xix(n)*xix(n)+c3*zetax(n)*xix(n)-c2*etax(n)*etax(n)+c4*zetax(n)*zetax(n);
            term14=  c4*xiy(n)*xiy(n)+c3*zetay(n)*xiy(n)-c2*etay(n)*etay(n)+c4*zetay(n)*zetay(n);
            term15=  c4*xiz(n)*xiz(n)+c3*zetaz(n)*xiz(n)-c2*etaz(n)*etaz(n)+c4*zetaz(n)*zetaz(n);
            term16=  c4*xiy(n)*xix(n)+c5*zetay(n)*xix(n)-c2*etay(n)*etax(n)+c5*xiy(n)*zetax(n)+c4*zetay(n)*zetax(n);
            term17= c3*(-etax(n)*xiy(n)+xix(n)*etay(n)+zetax(n)*etay(n)-etax(n)*zetay(n));
            term18=  c4*xiz(n)*xix(n)+c5*zetaz(n)*xix(n)-c2*etaz(n)*etax(n)+c5*xiz(n)*zetax(n)+c4*zetaz(n)*zetax(n);
            term19= c3*(-etax(n)*xiz(n)+xix(n)*etaz(n)+zetax(n)*etaz(n)-etax(n)*zetaz(n));
            term20=  c4*xiz(n)*xiy(n)+c5*zetaz(n)*xiy(n)-c2*etaz(n)*etay(n)+c5*xiz(n)*zetay(n)+c4*zetaz(n)*zetay(n);
            term21= c3*(-etay(n)*xiz(n)+xiy(n)*etaz(n)+zetay(n)*etaz(n)-etay(n)*zetaz(n));
            stiffl(10,22)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(11,22)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(12,22)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(12,23)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(11,23)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(10,23)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(12,24)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(11,24)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(10,24)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;
#undef C

#define C(i,j) c[(i)-1][(j)-1]
            // Rows 13-15
            term13=  c2*xix(n)*xix(n)+c2*etax(n)*etax(n)+c2*zetax(n)*zetax(n)+c1*etax(n)*xix(n)-c1*zetax(n)*xix(n)-c1*zetax(n)*etax(n);
            term14=  c2*xiy(n)*xiy(n)+c2*etay(n)*etay(n)+c2*zetay(n)*zetay(n)+c1*etay(n)*xiy(n)-c1*zetay(n)*xiy(n)-c1*zetay(n)*etay(n);
            term15=  c2*xiz(n)*xiz(n)+c2*etaz(n)*etaz(n)+c2*zetaz(n)*zetaz(n)+c1*etaz(n)*xiz(n)-c1*zetaz(n)*xiz(n)-c1*zetaz(n)*etaz(n);
            term16=  c2*xiy(n)*xix(n)+c3*etay(n)*xix(n)-c3*zetay(n)*xix(n)+c3*xiy(n)*etax(n)+c2*etay(n)*etax(n)-c3*zetay(n)*etax(n)-c3*xiy(n)*zetax(n)-c3*etay(n)*zetax(n)+c2*zetay(n)*zetax(n);
            term17 = 0.;
            term18=  c2*xiz(n)*xix(n)+c3*etaz(n)*xix(n)-c3*zetaz(n)*xix(n)+c3*xiz(n)*etax(n)+c2*etaz(n)*etax(n)-c3*zetaz(n)*etax(n)-c3*xiz(n)*zetax(n)-c3*etaz(n)*zetax(n)+c2*zetaz(n)*zetax(n);
            term19 = 0.;
            term20=  c2*xiz(n)*xiy(n)+c3*etaz(n)*xiy(n)-c3*zetaz(n)*xiy(n)+c3*xiz(n)*etay(n)+c2*etaz(n)*etay(n)-c3*zetaz(n)*etay(n)-c3*xiz(n)*zetay(n)-c3*etaz(n)*zetay(n)+c2*zetaz(n)*zetay(n);
            term21 = 0.;
            stiffl(13,13)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,6)+C(6,1))*term18+(C(4,6)+C(6,4))*term20;
            stiffl(14,14)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(4,2)+C(2,4))*term16+(C(4,5)+C(5,4))*term18+(C(2,5)+C(5,2))*term20;
            stiffl(13,14)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(1,2)+C(4,4))*term16+(C(1,5)+C(6,4))*term18+(C(4,5)+C(6,2))*term20;
            stiffl(15,15)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(6,5)+C(5,6))*term16+(C(6,3)+C(3,6))*term18+(C(5,3)+C(3,5))*term20;
            stiffl(14,15)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(4,5)+C(2,6))*term16+(C(4,3)+C(5,6))*term18+(C(2,3)+C(5,5))*term20;
            stiffl(13,15)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(1,5)+C(4,6))*term16+(C(1,3)+C(6,6))*term18+(C(4,3)+C(6,5))*term20;

            term13= -c2*xix(n)*xix(n)+c4*etax(n)*etax(n)-c3*zetax(n)*etax(n)+c4*zetax(n)*zetax(n);
            term14= -c2*xiy(n)*xiy(n)+c4*etay(n)*etay(n)-c3*zetay(n)*etay(n)+c4*zetay(n)*zetay(n);
            term15= -c2*xiz(n)*xiz(n)+c4*etaz(n)*etaz(n)-c3*zetaz(n)*etaz(n)+c4*zetaz(n)*zetaz(n);
            term16= -c2*xiy(n)*xix(n)+c4*etay(n)*etax(n)-c5*zetay(n)*etax(n)-c5*etay(n)*zetax(n)+c4*zetay(n)*zetax(n);
            term17= -c3*etax(n)*xiy(n)+c3*zetax(n)*xiy(n)+c3*xix(n)*etay(n)-c3*xix(n)*zetay(n);
            term18= -c2*xiz(n)*xix(n)+c4*etaz(n)*etax(n)-c5*zetaz(n)*etax(n)-c5*etaz(n)*zetax(n)+c4*zetaz(n)*zetax(n);
            term19= -c3*etax(n)*xiz(n)+c3*zetax(n)*xiz(n)+c3*xix(n)*etaz(n)-c3*xix(n)*zetaz(n);
            term20= -c2*xiz(n)*xiy(n)+c4*etaz(n)*etay(n)-c5*zetaz(n)*etay(n)-c5*etaz(n)*zetay(n)+c4*zetaz(n)*zetay(n);
            term21= -c3*etay(n)*xiz(n)+c3*zetay(n)*xiz(n)+c3*xiy(n)*etaz(n)-c3*xiy(n)*zetaz(n);
            stiffl(13,16)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(14,16)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(15,16)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(15,17)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(14,17)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(13,17)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(15,18)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(14,18)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(13,18)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;
#undef C

#define C(i,j) c[(i)-1][(j)-1]
            // Rows 13-15 continued (19-24)
            term13= -c4*xix(n)*xix(n)+c1*zetax(n)*xix(n)+c6*etax(n)*etax(n)-c4*zetax(n)*zetax(n);
            term14= -c4*xiy(n)*xiy(n)+c1*zetay(n)*xiy(n)+c6*etay(n)*etay(n)-c4*zetay(n)*zetay(n);
            term15= -c4*xiz(n)*xiz(n)+c1*zetaz(n)*xiz(n)+c6*etaz(n)*etaz(n)-c4*zetaz(n)*zetaz(n);
            term16= -c4*xiy(n)*xix(n)+c3*zetay(n)*xix(n)+c6*etay(n)*etax(n)+c3*xiy(n)*zetax(n)-c4*zetay(n)*zetax(n);
            term17= c5*(-etax(n)*xiy(n)+xix(n)*etay(n)-zetax(n)*etay(n)+etax(n)*zetay(n));
            term18= -c4*xiz(n)*xix(n)+c3*zetaz(n)*xix(n)+c6*etaz(n)*etax(n)+c3*xiz(n)*zetax(n)-c4*zetaz(n)*zetax(n);
            term19= c5*(-etax(n)*xiz(n)+xix(n)*etaz(n)-zetax(n)*etaz(n)+etax(n)*zetaz(n));
            term20= -c4*xiz(n)*xiy(n)+c3*zetaz(n)*xiy(n)+c6*etaz(n)*etay(n)+c3*xiz(n)*zetay(n)-c4*zetaz(n)*zetay(n);
            term21= c5*(-etay(n)*xiz(n)+xiy(n)*etaz(n)-zetay(n)*etaz(n)+etay(n)*zetaz(n));
            stiffl(13,19)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(14,19)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(15,19)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(15,20)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(14,20)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(13,20)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(15,21)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(14,21)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(13,21)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13=  c4*xix(n)*xix(n)+c3*etax(n)*xix(n)+c4*etax(n)*etax(n)-c2*zetax(n)*zetax(n);
            term14=  c4*xiy(n)*xiy(n)+c3*etay(n)*xiy(n)+c4*etay(n)*etay(n)-c2*zetay(n)*zetay(n);
            term15=  c4*xiz(n)*xiz(n)+c3*etaz(n)*xiz(n)+c4*etaz(n)*etaz(n)-c2*zetaz(n)*zetaz(n);
            term16=  c4*xiy(n)*xix(n)+c5*etay(n)*xix(n)+c5*xiy(n)*etax(n)+c4*etay(n)*etax(n)-c2*zetay(n)*zetax(n);
            term17= c3*(-zetax(n)*xiy(n)-zetax(n)*etay(n)+xix(n)*zetay(n)+etax(n)*zetay(n));
            term18=  c4*xiz(n)*xix(n)+c5*etaz(n)*xix(n)+c5*xiz(n)*etax(n)+c4*etaz(n)*etax(n)-c2*zetaz(n)*zetax(n);
            term19= c3*(-zetax(n)*xiz(n)-zetax(n)*etaz(n)+xix(n)*zetaz(n)+etax(n)*zetaz(n));
            term20=  c4*xiz(n)*xiy(n)+c5*etaz(n)*xiy(n)+c5*xiz(n)*etay(n)+c4*etaz(n)*etay(n)-c2*zetaz(n)*zetay(n);
            term21= c3*(-zetay(n)*xiz(n)-zetay(n)*etaz(n)+xiy(n)*zetaz(n)+etay(n)*zetaz(n));
            stiffl(13,22)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(14,22)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(15,22)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(15,23)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(14,23)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(13,23)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(15,24)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(14,24)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(13,24)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            // Rows 16-18
            term13=  c2*xix(n)*xix(n)-c1*etax(n)*xix(n)+c1*zetax(n)*xix(n)+c2*etax(n)*etax(n)-c1*zetax(n)*etax(n)+c2*zetax(n)*zetax(n);
            term14=  c2*xiy(n)*xiy(n)-c1*etay(n)*xiy(n)+c1*zetay(n)*xiy(n)+c2*etay(n)*etay(n)-c1*zetay(n)*etay(n)+c2*zetay(n)*zetay(n);
            term15=  c2*xiz(n)*xiz(n)-c1*etaz(n)*xiz(n)+c1*zetaz(n)*xiz(n)+c2*etaz(n)*etaz(n)-c1*zetaz(n)*etaz(n)+c2*zetaz(n)*zetaz(n);
            term16=  c2*xiy(n)*xix(n)-c3*etay(n)*xix(n)+c3*zetay(n)*xix(n)-c3*xiy(n)*etax(n)+c2*etay(n)*etax(n)-c3*zetay(n)*etax(n)+c3*xiy(n)*zetax(n)-c3*etay(n)*zetax(n)+c2*zetay(n)*zetax(n);
            term17 = 0.;
            term18=  c2*xiz(n)*xix(n)-c3*etaz(n)*xix(n)+c3*zetaz(n)*xix(n)-c3*xiz(n)*etax(n)+c2*etaz(n)*etax(n)-c3*zetaz(n)*etax(n)+c3*xiz(n)*zetax(n)-c3*etaz(n)*zetax(n)+c2*zetaz(n)*zetax(n);
            term19 = 0.;
            term20=  c2*xiy(n)*xiz(n)-c3*etay(n)*xiz(n)+c3*zetay(n)*xiz(n)-c3*xiy(n)*etaz(n)+c2*etay(n)*etaz(n)-c3*zetay(n)*etaz(n)+c3*xiy(n)*zetaz(n)-c3*etay(n)*zetaz(n)+c2*zetay(n)*zetaz(n);
            term21 = 0.;
            stiffl(16,16)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,6)+C(6,1))*term18+(C(4,6)+C(6,4))*term20;
            stiffl(17,17)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(5,4)+C(4,5))*term18+(C(5,2)+C(2,5))*term20;
            stiffl(16,17)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(6,4)+C(1,5))*term18+(C(6,2)+C(4,5))*term20;
            stiffl(18,18)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(3,6)+C(6,3))*term18+(C(3,5)+C(5,3))*term20;
            stiffl(17,18)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(5,6)+C(4,3))*term18+(C(5,5)+C(2,3))*term20;
            stiffl(16,18)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(6,6)+C(1,3))*term18+(C(6,5)+C(4,3))*term20;

            term13=  c4*xix(n)*xix(n)-c3*etax(n)*xix(n)+c4*etax(n)*etax(n)-c2*zetax(n)*zetax(n);
            term14=  c4*xiy(n)*xiy(n)-c3*etay(n)*xiy(n)+c4*etay(n)*etay(n)-c2*zetay(n)*zetay(n);
            term15=  c4*xiz(n)*xiz(n)-c3*etaz(n)*xiz(n)+c4*etaz(n)*etaz(n)-c2*zetaz(n)*zetaz(n);
            term16=  c4*xiy(n)*xix(n)-c5*etay(n)*xix(n)-c5*xiy(n)*etax(n)+c4*etay(n)*etax(n)-c2*zetay(n)*zetax(n);
            term17= c3*(zetax(n)*xiy(n)-zetax(n)*etay(n)-xix(n)*zetay(n)+etax(n)*zetay(n));
            term18=  c4*xiz(n)*xix(n)-c5*etaz(n)*xix(n)-c5*xiz(n)*etax(n)+c4*etaz(n)*etax(n)-c2*zetaz(n)*zetax(n);
            term19= c3*(zetax(n)*xiz(n)-zetax(n)*etaz(n)-xix(n)*zetaz(n)+etax(n)*zetaz(n));
            term20=  c4*xiz(n)*xiy(n)-c5*etaz(n)*xiy(n)-c5*xiz(n)*etay(n)+c4*etaz(n)*etay(n)-c2*zetaz(n)*zetay(n);
            term21= c3*(zetay(n)*xiz(n)-zetay(n)*etaz(n)-xiy(n)*zetaz(n)+etay(n)*zetaz(n));
            stiffl(16,19)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(17,19)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(18,19)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(18,20)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(17,20)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(16,20)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(18,21)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(17,21)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(16,21)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            term13= -c4*xix(n)*xix(n)-c1*zetax(n)*xix(n)+c6*etax(n)*etax(n)-c4*zetax(n)*zetax(n);
            term14= -c4*xiy(n)*xiy(n)-c1*zetay(n)*xiy(n)+c6*etay(n)*etay(n)-c4*zetay(n)*zetay(n);
            term15= -c4*xiz(n)*xiz(n)-c1*zetaz(n)*xiz(n)+c6*etaz(n)*etaz(n)-c4*zetaz(n)*zetaz(n);
            term16= -c4*xiy(n)*xix(n)-c3*zetay(n)*xix(n)+c6*etay(n)*etax(n)-c3*xiy(n)*zetax(n)-c4*zetay(n)*zetax(n);
            term17= c5*(etax(n)*xiy(n)-xix(n)*etay(n)-zetax(n)*etay(n)+etax(n)*zetay(n));
            term18= -c4*xiz(n)*xix(n)-c3*zetaz(n)*xix(n)+c6*etaz(n)*etax(n)-c3*xiz(n)*zetax(n)-c4*zetaz(n)*zetax(n);
            term19= c5*(etax(n)*xiz(n)-xix(n)*etaz(n)-zetax(n)*etaz(n)+etax(n)*zetaz(n));
            term20= -c4*xiz(n)*xiy(n)-c3*zetaz(n)*xiy(n)+c6*etaz(n)*etay(n)-c3*xiz(n)*zetay(n)-c4*zetaz(n)*zetay(n);
            term21= c5*(etay(n)*xiz(n)-xiy(n)*etaz(n)-zetay(n)*etaz(n)+etay(n)*zetaz(n));
            stiffl(16,22)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(17,22)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(18,22)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(18,23)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(17,23)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(16,23)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(18,24)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(17,24)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(16,24)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;
#undef C

#define C(i,j) c[(i)-1][(j)-1]
            // Rows 19-21
            term13=  c2*xix(n)*xix(n)-c1*etax(n)*xix(n)-c1*zetax(n)*xix(n)+c2*etax(n)*etax(n)+c1*zetax(n)*etax(n)+c2*zetax(n)*zetax(n);
            term14=  c2*xiy(n)*xiy(n)-c1*etay(n)*xiy(n)-c1*zetay(n)*xiy(n)+c2*etay(n)*etay(n)+c1*zetay(n)*etay(n)+c2*zetay(n)*zetay(n);
            term15=  c2*xiz(n)*xiz(n)-c1*etaz(n)*xiz(n)-c1*zetaz(n)*xiz(n)+c2*etaz(n)*etaz(n)+c1*zetaz(n)*etaz(n)+c2*zetaz(n)*zetaz(n);
            term16=  c2*xiy(n)*xix(n)-c3*etay(n)*xix(n)-c3*zetay(n)*xix(n)-c3*xiy(n)*etax(n)+c2*etay(n)*etax(n)+c3*zetay(n)*etax(n)-c3*xiy(n)*zetax(n)+c3*etay(n)*zetax(n)+c2*zetay(n)*zetax(n);
            term17 = 0.;
            term18=  c2*xiz(n)*xix(n)-c3*etaz(n)*xix(n)-c3*zetaz(n)*xix(n)-c3*xiz(n)*etax(n)+c2*etaz(n)*etax(n)+c3*zetaz(n)*etax(n)-c3*xiz(n)*zetax(n)+c3*etaz(n)*zetax(n)+c2*zetaz(n)*zetax(n);
            term19 = 0.;
            term20=  c2*xiz(n)*xiy(n)-c3*etaz(n)*xiy(n)-c3*zetaz(n)*xiy(n)-c3*xiz(n)*etay(n)+c2*etaz(n)*etay(n)+c3*zetaz(n)*etay(n)-c3*xiz(n)*zetay(n)+c3*etaz(n)*zetay(n)+c2*zetaz(n)*zetay(n);
            term21 = 0.;
            stiffl(19,19)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,6)+C(6,1))*term18+(C(4,6)+C(6,4))*term20;
            stiffl(20,20)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(5,4)+C(4,5))*term18+(C(5,2)+C(2,5))*term20;
            stiffl(19,20)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(6,4)+C(1,5))*term18+(C(6,2)+C(4,5))*term20;
            stiffl(21,21)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(3,6)+C(6,3))*term18+(C(3,5)+C(5,3))*term20;
            stiffl(20,21)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(5,6)+C(4,3))*term18+(C(5,5)+C(2,3))*term20;
            stiffl(19,21)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(6,6)+C(1,3))*term18+(C(6,5)+C(4,3))*term20;

            term13= -c2*xix(n)*xix(n)+c4*etax(n)*etax(n)+c3*zetax(n)*etax(n)+c4*zetax(n)*zetax(n);
            term14= -c2*xiy(n)*xiy(n)+c4*etay(n)*etay(n)+c3*zetay(n)*etay(n)+c4*zetay(n)*zetay(n);
            term15= -c2*xiz(n)*xiz(n)+c4*etaz(n)*etaz(n)+c3*zetaz(n)*etaz(n)+c4*zetaz(n)*zetaz(n);
            term16= -c2*xiy(n)*xix(n)+c4*etay(n)*etax(n)+c5*zetay(n)*etax(n)+c5*etay(n)*zetax(n)+c4*zetay(n)*zetax(n);
            term17= c3*(etax(n)*xiy(n)+zetax(n)*xiy(n)-xix(n)*etay(n)-xix(n)*zetay(n));
            term18= -c2*xiz(n)*xix(n)+c4*etaz(n)*etax(n)+c5*zetaz(n)*etax(n)+c5*etaz(n)*zetax(n)+c4*zetaz(n)*zetax(n);
            term19= c3*(etax(n)*xiz(n)+zetax(n)*xiz(n)-xix(n)*etaz(n)-xix(n)*zetaz(n));
            term20= -c2*xiz(n)*xiy(n)+c4*etaz(n)*etay(n)+c5*zetaz(n)*etay(n)+c5*etaz(n)*zetay(n)+c4*zetaz(n)*zetay(n);
            term21= c3*(etay(n)*xiz(n)+zetay(n)*xiz(n)-xiy(n)*etaz(n)-xiy(n)*zetaz(n));
            stiffl(19,22)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,4)-C(4,1))*term17+(C(1,6)+C(6,1))*term18+(C(1,6)-C(6,1))*term19+(C(4,6)+C(6,4))*term20+(C(4,6)-C(6,4))*term21;
            stiffl(20,22)=C(4,1)*term13+C(2,4)*term14+C(5,6)*term15+(C(2,1)+C(4,4))*term16+(C(4,4)-C(2,1))*term17+(C(5,1)+C(4,6))*term18+(C(4,6)-C(5,1))*term19+(C(5,4)+C(2,6))*term20+(C(2,6)-C(5,4))*term21;
            stiffl(21,22)=C(6,1)*term13+C(5,4)*term14+C(3,6)*term15+(C(5,1)+C(6,4))*term16+(C(6,4)-C(5,1))*term17+(C(3,1)+C(6,6))*term18+(C(6,6)-C(3,1))*term19+(C(3,4)+C(5,6))*term20+(C(5,6)-C(3,4))*term21;
            stiffl(21,23)=C(6,4)*term13+C(5,2)*term14+C(3,5)*term15+(C(5,4)+C(6,2))*term16+(C(6,2)-C(5,4))*term17+(C(3,4)+C(6,5))*term18+(C(6,5)-C(3,4))*term19+(C(3,2)+C(5,5))*term20+(C(5,5)-C(3,2))*term21;
            stiffl(20,23)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(4,2)-C(2,4))*term17+(C(5,4)+C(4,5))*term18+(C(4,5)-C(5,4))*term19+(C(5,2)+C(2,5))*term20+(C(2,5)-C(5,2))*term21;
            stiffl(19,23)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(1,2)-C(4,4))*term17+(C(6,4)+C(1,5))*term18+(C(1,5)-C(6,4))*term19+(C(6,2)+C(4,5))*term20+(C(4,5)-C(6,2))*term21;
            stiffl(21,24)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(6,5)-C(5,6))*term17+(C(3,6)+C(6,3))*term18+(C(6,3)-C(3,6))*term19+(C(3,5)+C(5,3))*term20+(C(5,3)-C(3,5))*term21;
            stiffl(20,24)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(4,5)-C(2,6))*term17+(C(5,6)+C(4,3))*term18+(C(4,3)-C(5,6))*term19+(C(5,5)+C(2,3))*term20+(C(2,3)-C(5,5))*term21;
            stiffl(19,24)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(1,5)-C(4,6))*term17+(C(6,6)+C(1,3))*term18+(C(1,3)-C(6,6))*term19+(C(6,5)+C(4,3))*term20+(C(4,3)-C(6,5))*term21;

            // Rows 22-24
            term13=  c2*xix(n)*xix(n)+c1*etax(n)*xix(n)+c1*zetax(n)*xix(n)+c2*etax(n)*etax(n)+c1*zetax(n)*etax(n)+c2*zetax(n)*zetax(n);
            term14=  c2*xiy(n)*xiy(n)+c1*etay(n)*xiy(n)+c1*zetay(n)*xiy(n)+c2*etay(n)*etay(n)+c1*zetay(n)*etay(n)+c2*zetay(n)*zetay(n);
            term15=  c2*xiz(n)*xiz(n)+c1*etaz(n)*xiz(n)+c1*zetaz(n)*xiz(n)+c2*etaz(n)*etaz(n)+c1*zetaz(n)*etaz(n)+c2*zetaz(n)*zetaz(n);
            term16=  c2*xiy(n)*xix(n)+c3*etay(n)*xix(n)+c3*zetay(n)*xix(n)+c3*xiy(n)*etax(n)+c2*etay(n)*etax(n)+c3*zetay(n)*etax(n)+c3*xiy(n)*zetax(n)+c3*etay(n)*zetax(n)+c2*zetay(n)*zetax(n);
            term17 = 0.;
            term18=  c2*xiz(n)*xix(n)+c3*etaz(n)*xix(n)+c3*zetaz(n)*xix(n)+c3*xiz(n)*etax(n)+c2*etaz(n)*etax(n)+c3*zetaz(n)*etax(n)+c3*xiz(n)*zetax(n)+c3*etaz(n)*zetax(n)+c2*zetaz(n)*zetax(n);
            term19 = 0.;
            term20=  c2*xiz(n)*xiy(n)+c3*etaz(n)*xiy(n)+c3*zetaz(n)*xiy(n)+c3*xiz(n)*etay(n)+c2*etaz(n)*etay(n)+c3*zetaz(n)*etay(n)+c3*xiz(n)*zetay(n)+c3*etaz(n)*zetay(n)+c2*zetaz(n)*zetay(n);
            term21 = 0.;
            stiffl(22,22)=C(1,1)*term13+C(4,4)*term14+C(6,6)*term15+(C(1,4)+C(4,1))*term16+(C(1,6)+C(6,1))*term18+(C(4,6)+C(6,4))*term20;
            stiffl(23,23)=C(4,4)*term13+C(2,2)*term14+C(5,5)*term15+(C(2,4)+C(4,2))*term16+(C(5,4)+C(4,5))*term18+(C(5,2)+C(2,5))*term20;
            stiffl(22,23)=C(1,4)*term13+C(4,2)*term14+C(6,5)*term15+(C(4,4)+C(1,2))*term16+(C(6,4)+C(1,5))*term18+(C(6,2)+C(4,5))*term20;
            stiffl(24,24)=C(6,6)*term13+C(5,5)*term14+C(3,3)*term15+(C(5,6)+C(6,5))*term16+(C(3,6)+C(6,3))*term18+(C(3,5)+C(5,3))*term20;
            stiffl(23,24)=C(4,6)*term13+C(2,5)*term14+C(5,3)*term15+(C(2,6)+C(4,5))*term16+(C(5,6)+C(4,3))*term18+(C(5,5)+C(2,3))*term20;
            stiffl(22,24)=C(1,6)*term13+C(4,5)*term14+C(6,3)*term15+(C(4,6)+C(1,5))*term16+(C(6,6)+C(1,3))*term18+(C(6,5)+C(4,3))*term20;
#undef C

            // symmetrize
            for (ii = 1; ii <= 24; ii++)
                for (jj = ii+1; jj <= 24; jj++)
                    stiffl(jj,ii) = stiffl(ii,jj);

            // assembly
            for (ii = 0; ii < 8; ii++) {
                n11i[0][ii] = n1[ii];
                iimax[ii]   = 1;
                n44i[ii]    = 1;
                if (islavept(n1[ii], 8, iseqr) == 0) n44i[ii] = 0;
                if (islavept(n1[ii], 11, iseqr) > 1) {
                    iimax[ii] = islavept(n1[ii], 11, iseqr);
                    for (int ii4 = 2; ii4 <= iimax[ii]; ii4++) {
                        n11i[ii4-1][ii] = islavept(n1[ii], 12+ii4-2, iseqr);
                        if (islavept(n11i[ii4-1][ii], 8, iseqr) == 0)
                            n44i[ii] = 0;
                    }
                }
                for (int ii4 = 1; ii4 <= iimax[ii]; ii4++)
                    n33i[ii4-1][ii] = 3*(n11i[ii4-1][ii]-1);
            }
            oj = 1./ooj(n);

            for (ii = 0; ii < 8; ii++) {
                if (n44i[ii] != 0) {
                    for (int ii3 = 1; ii3 <= iimax[ii]; ii3++) {
                        for (j = 1; j <= 3; j++) {
                            irow  = n2[ii]+j-1;
                            itst1 = n33i[ii3-1][ii]+j;
                            itst2 = ija(itst1);
                            sa(itst1) = sa(itst1) + stiffl(irow,irow)*oj;
                            i1 = 0;
                            for (j1 = 1; j1 <= 3; j1++) {
                                if (j1 != j) {
                                    sa(itst2+i1)  = sa(itst2+i1) + stiffl(irow, n2[ii]+j1-1)*oj;
                                    ija(itst2+i1) = n33i[ii3-1][ii]+j1;
                                    i1++;
                                }
                            }
                            j2sta = ija(itst1)+2;
                            j2end = ija(itst1+1)-1;
                            for (kk = 0; kk < 8; kk++) {
                                if (kk != ii) {
                                    for (j2 = j2sta; j2 <= j2end; j2 += 3) {
                                        if (ija(j2) == 0) {
                                            for (j1 = 1; j1 <= 3; j1++) {
                                                icol = n2[kk]+j1-1;
                                                j3   = j2+j1-1;
                                                sa(j3)  = sa(j3) + stiffl(irow,icol)*oj;
                                                ija(j3) = n33i[0][kk]+j1;
                                            }
                                            goto label2200;
                                        } else {
                                            for (int ii4 = 1; ii4 <= iimax[kk]; ii4++) {
                                                itst1 = n33i[ii4-1][kk]+1;
                                                if (ija(j2) == itst1) {
                                                    for (j1 = 1; j1 <= 3; j1++) {
                                                        icol = n2[kk]+j1-1;
                                                        j3   = j2+j1-1;
                                                        sa(j3) = sa(j3) + stiffl(irow,icol)*oj;
                                                    }
                                                    goto label2200;
                                                }
                                            }
                                        }
                                    }
                                    label2200:;
                                }
                            }
                        }
                    }
                } else {
                    for (j = 1; j <= 3; j++)
                        sa(3*(n1[ii]-1)+j) = 1.0;
                }
            }
        } // end if n1[1]!=n ...
    } // end do 200 (n loop)
}


// ============================================================
// s_gwrite
// ============================================================
void s_gwrite(int& nbl, int& idim, int& jdim, int& kdim,
              int& ista, int& iend, int& jsta, int& jend,
              int& ksta, int& kend, int& idim1, int& jdim1, int& kdim1,
              FortranArray3DRef<double> x, FortranArray3DRef<double> y,
              FortranArray3DRef<double> z)
{
    int iunit = 5000 + nbl;
    char fname[32];
    std::snprintf(fname, 32, "dgplot3d%04d.xyz", nbl);
    fortran_open_unit(iunit, fname, "w");
    FILE* fp = fortran_get_unit(iunit);
    fprintf(fp, " %d\n", 1);
    fprintf(fp, " %d %d %d\n", idim1, jdim1, kdim1);
    for (int k = ksta; k <= kend; k++)
        for (int j = jsta; j <= jend; j++)
            for (int i = ista; i <= iend; i++)
                fprintf(fp, " %20.12E\n", x(j,k,i));
    for (int k = ksta; k <= kend; k++)
        for (int j = jsta; j <= jend; j++)
            for (int i = ista; i <= iend; i++)
                fprintf(fp, " %20.12E\n", y(j,k,i));
    for (int k = ksta; k <= kend; k++)
        for (int j = jsta; j <= jend; j++)
            for (int i = ista; i <= iend; i++)
                fprintf(fp, " %20.12E\n", z(j,k,i));
    fortran_close_unit(iunit);
}

// ============================================================
// s_gread
// ============================================================
void s_gread(int& nbl, int& idim, int& jdim, int& kdim,
             int& ista, int& iend, int& jsta, int& jend,
             int& ksta, int& kend, int& idim1, int& jdim1, int& kdim1,
             FortranArray3DRef<double> x, FortranArray3DRef<double> y,
             FortranArray3DRef<double> z)
{
    int iunit = 5000 + nbl;
    char fname[32];
    std::snprintf(fname, 32, "dgplot3d%04d.xyz", nbl);
    fortran_open_unit(iunit, fname, "r");
    int nblks, id1, jd1, kd1;
    fortran_read_list(iunit, &nblks);
    fortran_read_list(iunit, &id1, &jd1, &kd1);
    for (int k = ksta; k <= kend; k++)
        for (int j = jsta; j <= jend; j++)
            for (int i = ista; i <= iend; i++)
                fortran_read_list(iunit, &x(j,k,i));
    for (int k = ksta; k <= kend; k++)
        for (int j = jsta; j <= jend; j++)
            for (int i = ista; i <= iend; i++)
                fortran_read_list(iunit, &y(j,k,i));
    for (int k = ksta; k <= kend; k++)
        for (int j = jsta; j <= jend; j++)
            for (int i = ista; i <= iend; i++)
                fortran_read_list(iunit, &z(j,k,i));
    fortran_close_unit(iunit);
}


// ============================================================
// gridgen_glyf1
// ============================================================
void gridgen_glyf1(int& nbl)
{
    char bgcontrl[9], fgcontrl[9];
    char filnm[6];
    char char1[4];
    char filnm2[25];
    double xmgrlx, xmgprlng;
    int iter, iflnum;

    fortran_open_unit(9700, "gridgen.options1", "r");
    // read(9700,12030) bgcontrl,fgcontrl  -- format(2a8)
    {
        FILE* fp = fortran_get_unit(9700);
        char line[80];
        if (fgets(line, sizeof(line), fp)) {
            std::memset(bgcontrl, ' ', 8); bgcontrl[8] = '\0';
            std::memset(fgcontrl, ' ', 8); fgcontrl[8] = '\0';
            int len = (int)strlen(line);
            if (len > 0 && line[len-1] == '\n') line[--len] = '\0';
            for (int i = 0; i < 8 && i < len; i++) bgcontrl[i] = line[i];
            for (int i = 0; i < 8 && 8+i < len; i++) fgcontrl[i] = line[8+i];
        }
    }
    fortran_read_list(9700, &xmgrlx, &xmgprlng, &iter);
    fortran_close_unit(9700);

    std::strcpy(filnm2, "dgplot3d_gridgen0001.glf");
    if (nbl < 10) {
        std::snprintf(char1, 4, "%1d", nbl);
        filnm2[18] = char1[0];
    } else if (nbl < 100) {
        std::snprintf(char1, 4, "%2d", nbl);
        filnm2[17] = char1[0]; filnm2[18] = char1[1];
    } else {
        std::snprintf(char1, 4, "%3d", nbl);
        filnm2[16] = char1[0]; filnm2[17] = char1[1]; filnm2[18] = char1[2];
    }
    fortran_open_unit(9900, filnm2, "w");
    std::strcpy(filnm, "fort.");
    iflnum = 5000 + nbl;
    FILE* fp = fortran_get_unit(9900);
    // format 11021
    fprintf(fp,
        "# Gridgen Journal File V1 (Gridgen 15.09 REL 3)\n"
        "# Created Sun Aug 13 14:47:02 2006\n"
        "#\n"
        "package require PWI_Glyph 1.6.9\n"
        "\n"
        "# Delete any existing grids and database entities.  Reset AS/W, defaults, and\n"
        "# tolerances.\n"
        "gg::memClear\n"
        "gg::aswDeleteBC -glob \"*\"\n"
        "gg::aswDeleteVC -glob \"*\"\n"
        "gg::aswSet GENERIC -dim 3\n"
        "gg::defReset\n"
        "gg::tolReset\n"
        "# Check for user input after each function, but delay screen updates until\n"
        "# script is finished.\n"
        "gg::updatePolicy INPUT_ONLY\n"
        "\n"
        "set _ggTemp_(1) [gg::blkImport \"%s%4d\" \\ \n"
        "   -style PLOT3D -form ASCII -precision DOUBLE]\n",
        filnm, iflnum);
    // format 11022: set _BL(1) [lindex $_ggTemp_(1) 0]
    fprintf(fp, "set _BL(%d) [lindex $_ggTemp_(%d) %d]\n", 1, 1, 0);
    // format 11025
    fprintf(fp, "unset _ggTemp_(1)\n");
    // format 11026
    fprintf(fp, "set _ggTemp_(2) [list $_BL(%d)]\n", 1);
    // format 11039
    fprintf(fp,
        "gg::blkEllSolverBegin $_ggTemp_(2)\n"
        "# gg::blkTFISolverRun $_ggTemp_(2)\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -bg_control %s\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -fg_control %s\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -spacing_calc INTERPOLATE\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -angle_calc INTERPOLATE\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -spacing_blend 2\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -angle_blend 2\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -mg_relax %3.2f\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -mg_prolongation %3.2f\n"
        "  gg::blkEllSolverStep -iterations %5d -nodisplay\n"
        "gg::blkEllSolverEnd\n"
        "unset _ggTemp_(2)\n"
        "gg::blkExport ALL \"%s%4d\" \\ \n"
        "   -style PLOT3D -form ASCII -precision DOUBLE\n"
        "gg::terminate\n",
        bgcontrl, fgcontrl, xmgrlx, xmgprlng, iter, filnm, iflnum);
    fortran_close_unit(9900);
}

// ============================================================
// gridgen_glyf2
// ============================================================
void gridgen_glyf2(int& nbl)
{
    char bgcontrl[9], fgcontrl[9];
    char filnm[6];
    char char1[4];
    char filnm2[25];
    double xmgrlx, xmgprlng;
    int iter, iflnum;

    fortran_open_unit(9700, "gridgen.options2", "r");
    {
        FILE* fp = fortran_get_unit(9700);
        char line[80];
        if (fgets(line, sizeof(line), fp)) {
            std::memset(bgcontrl, ' ', 8); bgcontrl[8] = '\0';
            std::memset(fgcontrl, ' ', 8); fgcontrl[8] = '\0';
            int len = (int)strlen(line);
            if (len > 0 && line[len-1] == '\n') line[--len] = '\0';
            for (int i = 0; i < 8 && i < len; i++) bgcontrl[i] = line[i];
            for (int i = 0; i < 8 && 8+i < len; i++) fgcontrl[i] = line[8+i];
        }
    }
    fortran_read_list(9700, &xmgrlx, &xmgprlng, &iter);
    fortran_close_unit(9700);

    std::strcpy(filnm2, "dgplot3d_gridgen0002.glf");
    if (nbl < 10) {
        std::snprintf(char1, 4, "%1d", nbl);
        filnm2[18] = char1[0];
    } else if (nbl < 100) {
        std::snprintf(char1, 4, "%2d", nbl);
        filnm2[17] = char1[0]; filnm2[18] = char1[1];
    } else {
        std::snprintf(char1, 4, "%3d", nbl);
        filnm2[16] = char1[0]; filnm2[17] = char1[1]; filnm2[18] = char1[2];
    }
    fortran_open_unit(9900, filnm2, "w");
    std::strcpy(filnm, "fort.");
    iflnum = 5000 + nbl;
    FILE* fp = fortran_get_unit(9900);
    fprintf(fp,
        "# Gridgen Journal File V1 (Gridgen 15.09 REL 3)\n"
        "# Created Sun Aug 13 14:47:02 2006\n"
        "#\n"
        "package require PWI_Glyph 1.6.9\n"
        "\n"
        "# Delete any existing grids and database entities.  Reset AS/W, defaults, and\n"
        "# tolerances.\n"
        "gg::memClear\n"
        "gg::aswDeleteBC -glob \"*\"\n"
        "gg::aswDeleteVC -glob \"*\"\n"
        "gg::aswSet GENERIC -dim 3\n"
        "gg::defReset\n"
        "gg::tolReset\n"
        "# Check for user input after each function, but delay screen updates until\n"
        "# script is finished.\n"
        "gg::updatePolicy INPUT_ONLY\n"
        "\n"
        "set _ggTemp_(1) [gg::blkImport \"%s%4d\" \\ \n"
        "   -style PLOT3D -form ASCII -precision DOUBLE]\n",
        filnm, iflnum);
    fprintf(fp, "set _BL(%d) [lindex $_ggTemp_(%d) %d]\n", 1, 1, 0);
    fprintf(fp, "unset _ggTemp_(1)\n");
    fprintf(fp, "set _ggTemp_(2) [list $_BL(%d)]\n", 1);
    fprintf(fp,
        "gg::blkEllSolverBegin $_ggTemp_(2)\n"
        "# gg::blkTFISolverRun $_ggTemp_(2)\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -bg_control %s\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -fg_control %s\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -spacing_calc INTERPOLATE\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -angle_calc INTERPOLATE\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -spacing_blend 2\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -angle_blend 2\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -mg_relax %3.2f\n"
        "  gg::blkEllSolverAtt $_ggTemp_(2) -mg_prolongation %3.2f\n"
        "  gg::blkEllSolverStep -iterations %5d -nodisplay\n"
        "gg::blkEllSolverEnd\n"
        "unset _ggTemp_(2)\n"
        "gg::blkExport ALL \"%s%4d\" \\ \n"
        "   -style PLOT3D -form ASCII -precision DOUBLE\n"
        "gg::terminate\n",
        bgcontrl, fgcontrl, xmgrlx, xmgprlng, iter, filnm, iflnum);
    fortran_close_unit(9900);
}

} // namespace updatedg_ns
