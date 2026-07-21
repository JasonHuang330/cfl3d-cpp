// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "qinter.h"
#include "add2x.h"
#include "addx.h"
#include "lead.h"
#include "termn8.h"
#include "writ_buf.h"
#include "grdmove.h"
#include "xyzintr.h"
#include "delintr.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <algorithm>

namespace qinter_ns {

// add2x and addx are implemented in their own modules; delegate to them.
void add2x(FortranArray4DRef<double> q, FortranArray4DRef<double> qc, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, FortranArray4DRef<double> q1, FortranArray4DRef<double> dq, FortranArray4DRef<double> wq, FortranArray3DRef<double> wqj, FortranArray3DRef<double> wqjk, int& js, int& ks, int& is, int& je, int& ke, int& ie, int& ipass, int& nbl, int& nblc, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& ll, int& myid)
{
    add2x_ns::add2x(q, qc, jdim, kdim, idim, jj2, kk2, ii2, q1, dq, wq, wqj, wqjk, js, ks, is, je, ke, ie, ipass, nbl, nblc, nou, bou, nbuf, ibufdim, ll, myid);
}

void addx(FortranArray4DRef<double> q, FortranArray4DRef<double> qq, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, FortranArray4DRef<double> q1, FortranArray4DRef<double> dq, FortranArray4DRef<double> wq, FortranArray4DRef<double> wqj, int& nbl, FortranArray3DRef<double> blank, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& ll, int& myid)
{
    addx_ns::addx(q, qq, jdim, kdim, idim, jj2, kk2, ii2, q1, dq, wq, wqj, nbl, blank, nou, bou, nbuf, ibufdim, ll, myid);
}

void qinter(int& iseq, int& lembed, int& maxl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& maxbl, int& maxgr, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, int& ngrid, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0, FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz, FortranArray1DRef<double> time2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae, int& maxsegdg, int& nummem)
{
    // COMMON block aliases
    int& jdim   = cmn_ginfo.jdim;
    int& kdim   = cmn_ginfo.kdim;
    int& idim   = cmn_ginfo.idim;
    int& jj2    = cmn_ginfo.jj2;
    int& kk2    = cmn_ginfo.kk2;
    int& ii2    = cmn_ginfo.ii2;
    int& nblc   = cmn_ginfo.nblc;
    int& js     = cmn_ginfo.js;
    int& ks     = cmn_ginfo.ks;
    int& is     = cmn_ginfo.is;
    int& je     = cmn_ginfo.je;
    int& ke     = cmn_ginfo.ke;
    int& ie     = cmn_ginfo.ie;
    int& lq     = cmn_ginfo.lq;
    int& lx     = cmn_ginfo.lx;
    int& ly     = cmn_ginfo.ly;
    int& lz     = cmn_ginfo.lz;
    int& lvis   = cmn_ginfo.lvis;
    int& lblk   = cmn_ginfo.lblk;
    int& lxib   = cmn_ginfo.lxib;
    int& ldeltj = cmn_ginfo.ldeltj;
    int& ldeltk = cmn_ginfo.ldeltk;
    int& ldelti = cmn_ginfo.ldelti;
    int& lxnm2  = cmn_ginfo.lxnm2;
    int& lynm2  = cmn_ginfo.lynm2;
    int& lznm2  = cmn_ginfo.lznm2;
    int& mseq   = cmn_info.mseq;
    int& mode   = cmn_mgrd.mode;
    int& ivmx   = cmn_maxiv.ivmx;

    // Local variables
    int iflg, igrid, nbl, iem, nblf, nblz;
    int iwk1, iwk2, iwk3, iwk4;
    int lqc, lq1c, lturc, lvisc;
    int nroom, mdim;
    int iuns, iseg;
    int inewg;
    int nblout;
    int ipass;
    int lxc, lyc, lzc;
    int lxnm2c, lynm2c, lznm2c;
    int ldeltjc, ldeltkc, ldeltic;
    int iunit11 = 11;
    int ineg1   = -1;

    mode = 0;
    iflg = 0;


    //
    // interpolate solution to finer mesh - global meshes
    //
    if (iseq <= mseq && iseq != 1) {

        if (myid == myhost) {
            fortran_write_unit(11, "***** BEGINNING SEQUENCING TO FINER LEVEL *****\n\n");
        }
        iflg = 1;

        for (igrid = 1; igrid <= ngrid; igrid++) {

            nbl = nblg(igrid);
            iem = iemg(igrid);
            if (iem > 0) continue;  // go to 1000

            // finer mesh nblf, coarser mesh nblz
            nblf = nblg(igrid) + (mseq - iseq);
            nblz = nblf + 1;

            if (mblk2nd(nblf) == myid) {

                lead_ns::lead(nblf, lw, lw2, maxbl);
                iwk1  = 1;
                iwk2  = iwk1 + jdim * kdim * idim * 7;
                iwk3  = iwk2 + jj2 * kk2 * ii2 * 7;
                lqc   = lw(1,  nblz);
                lq1c  = lw(16, nblz);
                lturc = lw(19, nblz);
                lvisc = lw(13, nblz);
                nroom = nwork - iwk3;
                mdim  = jdim * kk2 * ii2 * 7;
                if (nroom < mdim) {
                    nou(1) = std::min(nou(1) + 1, ibufdim);
                    std::snprintf(bou(nou(1), 1), 120,
                        " not enough work space for subroutine addx");
                    termn8_ns::termn8(myid, ineg1, ibufdim, nbuf, bou, nou);
                }

                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " interpolating solution on coarser block%4d  to   finer block%4d (grid%4d)",
                    nblz, nblf, igrid);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "   jdim,kdim,idim (finer grid)=%5d%5d%5d",
                    jdim, kdim, idim);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "   jj2,kk2,ii2  (coarser grid)=%5d%5d%5d",
                    jj2, kk2, ii2);

                // call addx(w(lq),w(lqc),...,5,myid)
                {
                    int ll5 = 5;
                    FortranArray4DRef<double> q_ref(&w(lq),      jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qc_ref(&w(lqc),    jj2,  kk2,  ii2,  5);
                    FortranArray4DRef<double> q1c_ref(&w(lq1c),  jj2,  kk2,  ii2,  5);
                    FortranArray4DRef<double> dq_ref(&wk(iwk1),  jdim, kdim, idim, 5);
                    FortranArray4DRef<double> wq_ref(&wk(iwk2),  jj2,  kk2,  ii2,  5);
                    FortranArray4DRef<double> wqj_ref(&wk(iwk3), jdim, kk2, ii2, 5);
                    FortranArray3DRef<double> blank_ref(&w(lblk), jdim, kdim, idim);
                    addx_ns::addx(q_ref, qc_ref, jdim, kdim, idim, jj2, kk2, ii2,
                                  q1c_ref, dq_ref, wq_ref, wqj_ref, nblz, blank_ref,
                                  nou, bou, nbuf, ibufdim, ll5, myid);
                }
                if (ivmx >= 4) {
                    nou(1) = std::min(nou(1) + 1, ibufdim);
                    std::snprintf(bou(nou(1), 1), 120,
                        "   interpolating turb quantities from coarser to finer block");
                    // call addx(w(lxib),w(lturc),...,nummem,myid)
                    {
                        int llnm = nummem;
                        FortranArray4DRef<double> xib_ref(&w(lxib),   jdim, kdim, idim, nummem);
                        FortranArray4DRef<double> turc_ref(&w(lturc),  jj2,  kk2,  ii2,  nummem);
                        FortranArray4DRef<double> q1c_ref(&w(lq1c),   jj2,  kk2,  ii2,  nummem);
                        FortranArray4DRef<double> dq_ref(&wk(iwk1),   jdim, kdim, idim, nummem);
                        FortranArray4DRef<double> wq_ref(&wk(iwk2),   jj2,  kk2,  ii2,  nummem);
                        FortranArray4DRef<double> wqj_ref(&wk(iwk3),  jdim, kk2, ii2, nummem);
                        FortranArray3DRef<double> blank_ref(&w(lblk),  jdim, kdim, idim);
                        addx_ns::addx(xib_ref, turc_ref, jdim, kdim, idim, jj2, kk2, ii2,
                                      q1c_ref, dq_ref, wq_ref, wqj_ref, nblz, blank_ref,
                                      nou, bou, nbuf, ibufdim, llnm, myid);
                    }
                    // call addx(w(lvis),w(lvisc),...,1,myid)
                    {
                        int ll1 = 1;
                        FortranArray4DRef<double> vis_ref(&w(lvis),    jdim, kdim, idim, 1);
                        FortranArray4DRef<double> visc_ref(&w(lvisc),  jj2,  kk2,  ii2,  1);
                        FortranArray4DRef<double> q1c_ref(&w(lq1c),   jj2,  kk2,  ii2,  1);
                        FortranArray4DRef<double> dq_ref(&wk(iwk1),   jdim, kdim, idim, 1);
                        FortranArray4DRef<double> wq_ref(&wk(iwk2),   jj2,  kk2,  ii2,  1);
                        FortranArray4DRef<double> wqj_ref(&wk(iwk3),  jdim, kk2, ii2, 1);
                        FortranArray3DRef<double> blank_ref(&w(lblk),  jdim, kdim, idim);
                        addx_ns::addx(vis_ref, visc_ref, jdim, kdim, idim, jj2, kk2, ii2,
                                      q1c_ref, dq_ref, wq_ref, wqj_ref, nblz, blank_ref,
                                      nou, bou, nbuf, ibufdim, ll1, myid);
                    }
                }

                iuns = std::max({itrans(nblf), irotat(nblf), idefrm(nblf)});
                if (iuns > 0) {
                    if (idefrm(nblf) > 0) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            "   interpolating deforming mesh data from coarser to finer block");
                        for (iseg = 1; iseg <= nsegdfrm(nblf); iseg++) {
                            xorgae(nblf, iseg) = xorgae(nblz, iseg);
                            yorgae(nblf, iseg) = yorgae(nblz, iseg);
                            zorgae(nblf, iseg) = zorgae(nblz, iseg);
                            thtxae(nblf, iseg) = thtxae(nblz, iseg);
                            thtyae(nblf, iseg) = thtyae(nblz, iseg);
                            thtzae(nblf, iseg) = thtzae(nblz, iseg);
                        }
                        lxc     = lw(10, nblz);
                        lyc     = lw(11, nblz);
                        lzc     = lw(12, nblz);
                        lxnm2c  = lw(56, nblz);
                        lynm2c  = lw(57, nblz);
                        lznm2c  = lw(58, nblz);
                        ldeltjc = lw(59, nblz);
                        ldeltkc = lw(60, nblz);
                        ldeltic = lw(61, nblz);

                        // call xyzintr for x,y,z
                        {
                            FortranArray3DRef<double> xf_ref(&w(lx),   jdim, kdim, idim);
                            FortranArray3DRef<double> yf_ref(&w(ly),   jdim, kdim, idim);
                            FortranArray3DRef<double> zf_ref(&w(lz),   jdim, kdim, idim);
                            FortranArray3DRef<double> xc_ref(&w(lxc),  jj2,  kk2,  ii2);
                            FortranArray3DRef<double> yc_ref(&w(lyc),  jj2,  kk2,  ii2);
                            FortranArray3DRef<double> zc_ref(&w(lzc),  jj2,  kk2,  ii2);
                            xyzintr_ns::xyzintr(jdim, kdim, idim, jj2, kk2, ii2,
                                                xf_ref, yf_ref, zf_ref,
                                                xc_ref, yc_ref, zc_ref);
                        }
                        // call xyzintr for xnm2,ynm2,znm2
                        {
                            FortranArray3DRef<double> xnm2f_ref(&w(lxnm2),  jdim, kdim, idim);
                            FortranArray3DRef<double> ynm2f_ref(&w(lynm2),  jdim, kdim, idim);
                            FortranArray3DRef<double> znm2f_ref(&w(lznm2),  jdim, kdim, idim);
                            FortranArray3DRef<double> xnm2c_ref(&w(lxnm2c), jj2,  kk2,  ii2);
                            FortranArray3DRef<double> ynm2c_ref(&w(lynm2c), jj2,  kk2,  ii2);
                            FortranArray3DRef<double> znm2c_ref(&w(lznm2c), jj2,  kk2,  ii2);
                            xyzintr_ns::xyzintr(jdim, kdim, idim, jj2, kk2, ii2,
                                                xnm2f_ref, ynm2f_ref, znm2f_ref,
                                                xnm2c_ref, ynm2c_ref, znm2c_ref);
                        }
                        // call delintr for deltj,deltk,delti
                        {
                            int idim1 = idim - 1;
                            int kdim1 = kdim - 1;
                            int jdim1 = jdim - 1;
                            int ii21  = ii2 - 1;
                            int kk21  = kk2 - 1;
                            int jj21  = jj2 - 1;
                            FortranArray4DRef<double> deltjf_ref(&w(ldeltj),  jdim,  kdim,  idim1, 5);
                            FortranArray4DRef<double> deltkf_ref(&w(ldeltk),  jdim,  kdim1, idim,  5);
                            FortranArray4DRef<double> deltif_ref(&w(ldelti),  jdim1, kdim,  idim,  5);
                            FortranArray4DRef<double> deltjc_ref(&w(ldeltjc), jj2,   kk2,   ii21,  5);
                            FortranArray4DRef<double> deltkc_ref(&w(ldeltkc), jj2,   kk21,  ii2,   5);
                            FortranArray4DRef<double> deltic_ref(&w(ldeltic), jj21,  kk2,   ii2,   5);
                            delintr_ns::delintr(jdim, kdim, idim, jj2, kk2, ii2,
                                                deltjf_ref, deltkf_ref, deltif_ref,
                                                deltjc_ref, deltkc_ref, deltic_ref);
                        }
                    } // end if idefrm(nblf) > 0

                    // update position of finer mesh to current position of
                    // coarser mesh to account for rigid grid motion
                    if (itrans(nblf) > 0 || irotat(nblf) > 0) {
                        xorig(nblf)  = xorig(nblz);
                        yorig(nblf)  = yorig(nblz);
                        zorig(nblf)  = zorig(nblz);
                        thetax(nblf) = thetax(nblz);
                        thetay(nblf) = thetay(nblz);
                        thetaz(nblf) = thetaz(nblz);
                        time2(nblf)  = time2(nblz);
                        FortranArray3DRef<double> xf_ref(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> yf_ref(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> zf_ref(&w(lz), jdim, kdim, idim);
                        grdmove_ns::grdmove(nblf, jdim, kdim, idim,
                                            xf_ref, yf_ref, zf_ref,
                                            xorig0(nblf), yorig0(nblf), zorig0(nblf),
                                            xorig(nblf),  yorig(nblf),  zorig(nblf),
                                            thetax(nblf), thetay(nblf), thetaz(nblf));
                    }
                } // end if iuns > 0

            } // end if mblk2nd(nblf) == myid

            nblout = nblg(igrid);
            writ_buf_ns::writ_buf(nblout, iunit11, nou, bou, nbuf, ibufdim,
                                  myhost, myid, mycomm, mblk2nd, maxbl);

        } // end do 1000 (igrid loop for global meshes)

        if (myid == myhost) {
            if (lembed == 0 || iseq != mseq) {
                fortran_write_unit(11, "\n***** ENDING SEQUENCING TO FINER LEVEL *****\n");
                return;
            }
        }

    } // end if iseq <= mseq && iseq != 1


    //
    // embedded grids - either all new  (iseq>1)
    //                  or some new     (iseq=1   inewg=1)
    //
    // cycle through levels  coarsest to finest
    //
    if (lembed > 0) {

        for (int level = mseq + 1; level <= maxl; level++) {
            for (igrid = 1; igrid <= ngrid; igrid++) {
                nbl   = nblg(igrid);
                iem   = iemg(igrid);
                inewg = inewgg(igrid);

                if (iem == 0) continue;                  // go to 2100
                if (iseq == 1 && inewg == 0) continue;  // go to 2100
                if (levelg(nbl) != level) continue;      // go to 2100

                // embedded mesh nbl; coarser mesh nblc
                lead_ns::lead(nbl, lw, lw2, maxbl);
                jj2 = jdimg(nblc);
                kk2 = kdimg(nblc);
                ii2 = idimg(nblc);
                lqc   = lw(1,  nblc);
                lq1c  = lw(16, nblc);
                lturc = lw(19, nblc);
                lvisc = lw(13, nblc);
                iwk1 = 1;
                iwk2 = iwk1 + jdim * kdim * idim * 7;
                iwk3 = iwk2 + jj2 * kk2 * ii2 * 7;
                iwk4 = iwk3 + jdim * kk2 * ii2;
                nroom = nwork - iwk4;
                mdim  = jdim * kdim * ii2;
                if (nroom < mdim) {
                    nou(1) = std::min(nou(1) + 1, ibufdim);
                    std::snprintf(bou(nou(1), 1), 120,
                        " not enough work space for subroutine add2x");
                    termn8_ns::termn8(myid, ineg1, ibufdim, nbuf, bou, nou);
                }

                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " interpolating solution on coarser block%4d  to embeded block%4d (grid%4d)",
                    nblc, nbl, igridg(nbl));
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "   jdim,kdim,idim (finer grid)=%5d%5d%5d",
                    jdim, kdim, idim);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "   jj2,kk2,ii2  (coarser grid)=%5d%5d%5d",
                    jj2, kk2, ii2);

                ipass = 1;
                // call add2x(w(lq),w(lqc),...,5,myid)
                {
                    int ll5 = 5;
                    FortranArray4DRef<double> q_ref(&w(lq),       jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qc_ref(&w(lqc),     jj2,  kk2,  ii2,  5);
                    FortranArray4DRef<double> q1c_ref(&w(lq1c),   jj2,  kk2,  ii2,  5);
                    FortranArray4DRef<double> dq_ref(&wk(iwk1),   jdim, kdim, idim, 5);
                    FortranArray4DRef<double> wq_ref(&wk(iwk2),   jj2,  kk2,  ii2,  5);
                    FortranArray3DRef<double> wqj_ref(&wk(iwk3),  jdim, kk2,  ii2);
                    FortranArray3DRef<double> wqjk_ref(&wk(iwk4), jdim, kdim, ii2);
                    add2x_ns::add2x(q_ref, qc_ref, jdim, kdim, idim, jj2, kk2, ii2,
                                    q1c_ref, dq_ref, wq_ref, wqj_ref, wqjk_ref,
                                    js, ks, is, je, ke, ie, ipass, nbl, nblc,
                                    nou, bou, nbuf, ibufdim, ll5, myid);
                }
                if (ivmx >= 4) {
                    nou(1) = std::min(nou(1) + 1, ibufdim);
                    std::snprintf(bou(nou(1), 1), 120,
                        "   interpolating turb quantities from coarser to embeded block");
                    // call add2x(w(lxib),w(lturc),...,nummem,myid)
                    {
                        int llnm = nummem;
                        FortranArray4DRef<double> xib_ref(&w(lxib),    jdim, kdim, idim, nummem);
                        FortranArray4DRef<double> turc_ref(&w(lturc),   jj2,  kk2,  ii2,  nummem);
                        FortranArray4DRef<double> q1c_ref(&w(lq1c),    jj2,  kk2,  ii2,  nummem);
                        FortranArray4DRef<double> dq_ref(&wk(iwk1),    jdim, kdim, idim, nummem);
                        FortranArray4DRef<double> wq_ref(&wk(iwk2),    jj2,  kk2,  ii2,  nummem);
                        FortranArray3DRef<double> wqj_ref(&wk(iwk3),   jdim, kk2,  ii2);
                        FortranArray3DRef<double> wqjk_ref(&wk(iwk4),  jdim, kdim, ii2);
                        add2x_ns::add2x(xib_ref, turc_ref, jdim, kdim, idim, jj2, kk2, ii2,
                                        q1c_ref, dq_ref, wq_ref, wqj_ref, wqjk_ref,
                                        js, ks, is, je, ke, ie, ipass, nbl, nblc,
                                        nou, bou, nbuf, ibufdim, llnm, myid);
                    }
                    // call add2x(w(lvis),w(lvisc),...,1,myid)
                    {
                        int ll1 = 1;
                        FortranArray4DRef<double> vis_ref(&w(lvis),    jdim, kdim, idim, 1);
                        FortranArray4DRef<double> visc_ref(&w(lvisc),  jj2,  kk2,  ii2,  1);
                        FortranArray4DRef<double> q1c_ref(&w(lq1c),   jj2,  kk2,  ii2,  1);
                        FortranArray4DRef<double> dq_ref(&wk(iwk1),   jdim, kdim, idim, 1);
                        FortranArray4DRef<double> wq_ref(&wk(iwk2),   jj2,  kk2,  ii2,  1);
                        FortranArray3DRef<double> wqj_ref(&wk(iwk3),  jdim, kk2,  ii2);
                        FortranArray3DRef<double> wqjk_ref(&wk(iwk4), jdim, kdim, ii2);
                        add2x_ns::add2x(vis_ref, visc_ref, jdim, kdim, idim, jj2, kk2, ii2,
                                        q1c_ref, dq_ref, wq_ref, wqj_ref, wqjk_ref,
                                        js, ks, is, je, ke, ie, ipass, nbl, nblc,
                                        nou, bou, nbuf, ibufdim, ll1, myid);
                    }
                }

                iuns = std::max({itrans(nbl), irotat(nbl), idefrm(nbl)});
                if (iuns > 0) {
                    // update position of embedded mesh to current position of
                    // coarser (global) mesh
                    xorig(nbl)  = xorig(nblc);
                    yorig(nbl)  = yorig(nblc);
                    zorig(nbl)  = zorig(nblc);
                    thetax(nbl) = thetax(nblc);
                    thetay(nbl) = thetay(nblc);
                    thetaz(nbl) = thetaz(nblc);
                    time2(nbl)  = time2(nblc);

                    if (idefrm(nbl) == 0) {
                        FortranArray3DRef<double> xf_ref(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> yf_ref(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> zf_ref(&w(lz), jdim, kdim, idim);
                        grdmove_ns::grdmove(nbl, jdim, kdim, idim,
                                            xf_ref, yf_ref, zf_ref,
                                            xorig0(nbl), yorig0(nbl), zorig0(nbl),
                                            xorig(nbl),  yorig(nbl),  zorig(nbl),
                                            thetax(nbl), thetay(nbl), thetaz(nbl));
                    }
                }

                nblout = nblg(igrid);
                writ_buf_ns::writ_buf(nblout, iunit11, nou, bou, nbuf, ibufdim,
                                      myhost, myid, mycomm, mblk2nd, maxbl);

            } // end do 2100 (igrid loop for embedded grids)
        } // end do 2000 (level loop)

        if (myid == myhost) {
            if (iflg > 0) {
                fortran_write_unit(11, "\n***** ENDING SEQUENCING TO FINER LEVEL *****\n");
            }
        }

    } // end if lembed > 0

    return;

} // end qinter

} // namespace qinter_ns
