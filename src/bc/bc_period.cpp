// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc_period.h"
#include "lead.h"
#include "termn8.h"
#include "chkrot.h"
#include "bc2005.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstring>
#include <algorithm>

namespace bc_period_ns {

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
    return;
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
               int& myid, int& myhost, int& mycomm,
               FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou,
               FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
               FortranArray2DRef<int> istat2, int& istat_size,
               FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, int& nummem)
{
    // COMMON block references
    int& jdim    = cmn_ginfo.jdim;
    int& kdim    = cmn_ginfo.kdim;
    int& idim    = cmn_ginfo.idim;
    int& lq      = cmn_ginfo.lq;
    int& lqj0    = cmn_ginfo.lqj0;
    int& lqk0    = cmn_ginfo.lqk0;
    int& lqi0    = cmn_ginfo.lqi0;
    int& lx      = cmn_ginfo.lx;
    int& ly      = cmn_ginfo.ly;
    int& lz      = cmn_ginfo.lz;
    int& lvis    = cmn_ginfo.lvis;
    int& lxib    = cmn_ginfo.lxib;
    int& ltj0    = cmn_ginfo.ltj0;
    int& ltk0    = cmn_ginfo.ltk0;
    int& lti0    = cmn_ginfo.lti0;
    int& lvj0    = cmn_ginfo.lvj0;
    int& lvk0    = cmn_ginfo.lvk0;
    int& lvi0    = cmn_ginfo.lvi0;
    int& level   = cmn_mgrd.level;
    int& isklton = cmn_sklton.isklton;
    int& iexp    = cmn_zero.iexp;

    // is_perbc COMMON (1-based Fortran arrays stored as C 0-based)
    int* is_prd  = cmn_is_perbc.is_prd;
    int* ie_prd  = cmn_is_perbc.ie_prd;
    int& nbcprd  = cmn_is_perbc.nbcprd;

    // epsrot = max(1.e-09, 10.**(-iexp+1))
    epsrot = std::max(1.e-09, std::pow(10.0, (double)(-iexp + 1)));

    if (ntime > 0 && nbcprd > 0) {

        int iwk_indx = 1;

        int lcnt_start = is_prd[level - 1];
        int lcnt_end   = ie_prd[level - 1];

        for (int lcnt = lcnt_start; lcnt <= lcnt_end; lcnt++) {
            int nbll    = isav_prd(lcnt, 1);
            int nblp    = isav_prd(lcnt, 12);
            int nd_recv = mblk2nd(nbll);
            int nd_srce = mblk2nd(nblp);

            if (nd_srce == myid && nd_recv == myid) {

                if (iadvance(nbll) >= 0) {

                    lead_ns::lead(nbll, lw, lw2, maxbl);

                    int n     = lcnt;
                    int nface = isav_prd(lcnt, 2);
                    int ista  = isav_prd(lcnt, 3);
                    int iend  = isav_prd(lcnt, 4);
                    int jsta  = isav_prd(lcnt, 5);
                    int jend  = isav_prd(lcnt, 6);
                    int ksta  = isav_prd(lcnt, 7);
                    int kend  = isav_prd(lcnt, 8);
                    int mdim  = isav_prd(lcnt, 9);
                    int ndim  = isav_prd(lcnt, 10);
                    int nseg  = isav_prd(lcnt, 11);
                    int ldata = lwdat(nbll, nseg, nface);

                    char filname[80];
                    if (nface == 1) {
                        std::memcpy(filname, bcfiles(bcfilei(nbll, nseg, 1)), 80);
                    } else if (nface == 2) {
                        std::memcpy(filname, bcfiles(bcfilei(nbll, nseg, 2)), 80);
                    } else if (nface == 3) {
                        std::memcpy(filname, bcfiles(bcfilej(nbll, nseg, 1)), 80);
                    } else if (nface == 4) {
                        std::memcpy(filname, bcfiles(bcfilej(nbll, nseg, 2)), 80);
                    } else if (nface == 5) {
                        std::memcpy(filname, bcfiles(bcfilek(nbll, nseg, 1)), 80);
                    } else if (nface == 6) {
                        std::memcpy(filname, bcfiles(bcfilek(nbll, nseg, 2)), 80);
                    }

                    int idimp = idimg(nblp);
                    int jdimp = jdimg(nblp);
                    int kdimp = kdimg(nblp);

                    int maxdims = 0;
                    if (nface == 1 || nface == 2) maxdims = jdimp * kdimp;
                    if (nface == 3 || nface == 4) maxdims = kdimp * idimp;
                    if (nface == 5 || nface == 6) maxdims = jdimp * idimp;


                    // k = constant interface (nface 5 or 6)
                    if (nface == 5 || nface == 6) {
                        if (isklton == 1) {
                            int kcheck = iwk_indx + maxdims * 9 - 1;
                            if (kcheck > nwork) {
                                nou(1) = std::min(nou(1) + 1, ibufdim);
                                std::snprintf(bou(nou(1), 1), 120,
                                    " stop in bc_period...insufficient wk storage for chkrot");
                                int m1 = -1;
                                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                            }
                            FortranArray3DRef<double> wx(&w(lx),   jdim, kdim, idim);
                            FortranArray3DRef<double> wy(&w(ly),   jdim, kdim, idim);
                            FortranArray3DRef<double> wz(&w(lz),   jdim, kdim, idim);
                            FortranArray3DRef<double> xp(&w(lw(10, nblp)), jdimp, kdimp, idimp);
                            FortranArray3DRef<double> yp(&w(lw(11, nblp)), jdimp, kdimp, idimp);
                            FortranArray3DRef<double> zp(&w(lw(12, nblp)), jdimp, kdimp, idimp);
                            FortranArray4DRef<double> bcdata(&w(ldata), mdim, ndim, 2, 12);
                            FortranArray3DRef<double> xyzjp(&wk(iwk_indx),               maxdims, 3, 1);
                            FortranArray3DRef<double> xyzkp(&wk(iwk_indx + maxdims * 3),  maxdims, 3, 1);
                            FortranArray3DRef<double> xyzip(&wk(iwk_indx + maxdims * 6),  maxdims, 3, 1);
                            chkrot_ns::chkrot(nbll, jdim, kdim, idim,
                                              wx, wy, wz,
                                              nblp, jdimp, kdimp, idimp,
                                              xp, yp, zp,
                                              nface, bcdata,
                                              xyzjp, xyzkp, xyzip,
                                              ista, iend, jsta, jend, ksta, kend,
                                              mdim, ndim, lcnt,
                                              xorig, yorig, zorig,
                                              maxbl, period_miss, lbcprd,
                                              nou, bou, nbuf, ibufdim, myid);
                        }
                        if (isklton == 1) {
                            int kcheck = iwk_indx + maxdims * 36 - 1;
                            if (kcheck > nwork) {
                                nou(1) = std::min(nou(1) + 1, ibufdim);
                                std::snprintf(bou(nou(1), 1), 120,
                                    " stop in bc_period...insufficient wk storage for bctype 2005");
                                int m1 = -1;
                                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                            }
                        }
                        int iskltsav = isklton;
                        isklton = 0;
                        {
                            FortranArray4DRef<double> wq   (&w(lq),   jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wqj0 (&w(lqj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wqk0 (&w(lqk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wqi0 (&w(lqi0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wxib (&w(lxib), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wtj0 (&w(ltj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wtk0 (&w(ltk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wti0 (&w(lti0), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> wvis (&w(lvis), jdim, kdim, idim);
                            FortranArray4DRef<double> wvj0 (&w(lvj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wvk0 (&w(lvk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wvi0 (&w(lvi0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> bcdata(&w(ldata), mdim, ndim, 2, 12);
                            FortranArray4DRef<double> qp   (&w(lw(1,  nblp)), jdimp, kdimp, idimp, 5);
                            FortranArray3DRef<double> vp   (&w(lw(13, nblp)), jdimp, kdimp, idimp);
                            FortranArray4DRef<double> tp   (&w(lw(19, nblp)), jdimp, kdimp, idimp, 2);
                            FortranArray4DRef<double> qrotj(&wk(iwk_indx),                maxdims, 10, 1, 1);
                            FortranArray4DRef<double> qrotk(&wk(iwk_indx + maxdims * 10), maxdims, 10, 1, 1);
                            FortranArray4DRef<double> qroti(&wk(iwk_indx + maxdims * 20), maxdims, 10, 1, 1);
                            bc2005_ns::bc2005(jdim, kdim, idim,
                                              wq, wqj0, wqk0, wqi0,
                                              ista, iend, jsta, jend, ksta, kend,
                                              nface, wxib, wtj0, wtk0, wti0,
                                              wvis, wvj0, wvk0, wvi0,
                                              mdim, ndim, bcdata, filname,
                                              qp, vp, tp,
                                              jdimp, kdimp, idimp,
                                              qrotj, qrotk, qroti,
                                              nbll, nblp,
                                              nou, bou, nbuf, ibufdim,
                                              myid, mblk2nd, maxbl, nummem);
                        }
                        isklton = iskltsav;
                    } // end nface 5/6

                    // j = constant interface (nface 3 or 4)
                    if (nface == 3 || nface == 4) {
                        if (isklton == 1) {
                            int kcheck = iwk_indx + maxdims * 9 - 1;
                            if (kcheck > nwork) {
                                nou(1) = std::min(nou(1) + 1, ibufdim);
                                std::snprintf(bou(nou(1), 1), 120,
                                    " stop in bc_period...insufficient wk storage for chkrot");
                                int m1 = -1;
                                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                            }
                            FortranArray3DRef<double> wx(&w(lx),   jdim, kdim, idim);
                            FortranArray3DRef<double> wy(&w(ly),   jdim, kdim, idim);
                            FortranArray3DRef<double> wz(&w(lz),   jdim, kdim, idim);
                            FortranArray3DRef<double> xp(&w(lw(10, nblp)), jdimp, kdimp, idimp);
                            FortranArray3DRef<double> yp(&w(lw(11, nblp)), jdimp, kdimp, idimp);
                            FortranArray3DRef<double> zp(&w(lw(12, nblp)), jdimp, kdimp, idimp);
                            FortranArray4DRef<double> bcdata(&w(ldata), mdim, ndim, 2, 12);
                            FortranArray3DRef<double> xyzjp(&wk(iwk_indx),               maxdims, 3, 1);
                            FortranArray3DRef<double> xyzkp(&wk(iwk_indx + maxdims * 3),  maxdims, 3, 1);
                            FortranArray3DRef<double> xyzip(&wk(iwk_indx + maxdims * 6),  maxdims, 3, 1);
                            chkrot_ns::chkrot(nbll, jdim, kdim, idim,
                                              wx, wy, wz,
                                              nblp, jdimp, kdimp, idimp,
                                              xp, yp, zp,
                                              nface, bcdata,
                                              xyzjp, xyzkp, xyzip,
                                              ista, iend, jsta, jend, ksta, kend,
                                              mdim, ndim, lcnt,
                                              xorig, yorig, zorig,
                                              maxbl, period_miss, lbcprd,
                                              nou, bou, nbuf, ibufdim, myid);
                        }
                        if (isklton == 1) {
                            int kcheck = iwk_indx + maxdims * 36 - 1;
                            if (kcheck > nwork) {
                                nou(1) = std::min(nou(1) + 1, ibufdim);
                                std::snprintf(bou(nou(1), 1), 120,
                                    " stop in bc_period...insufficient wk storage for bctype 2005");
                                int m1 = -1;
                                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                            }
                        }
                        int iskltsav = isklton;
                        isklton = 0;
                        {
                            FortranArray4DRef<double> wq   (&w(lq),   jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wqj0 (&w(lqj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wqk0 (&w(lqk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wqi0 (&w(lqi0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wxib (&w(lxib), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wtj0 (&w(ltj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wtk0 (&w(ltk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wti0 (&w(lti0), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> wvis (&w(lvis), jdim, kdim, idim);
                            FortranArray4DRef<double> wvj0 (&w(lvj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wvk0 (&w(lvk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wvi0 (&w(lvi0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> bcdata(&w(ldata), mdim, ndim, 2, 12);
                            FortranArray4DRef<double> qp   (&w(lw(1,  nblp)), jdimp, kdimp, idimp, 5);
                            FortranArray3DRef<double> vp   (&w(lw(13, nblp)), jdimp, kdimp, idimp);
                            FortranArray4DRef<double> tp   (&w(lw(19, nblp)), jdimp, kdimp, idimp, 2);
                            FortranArray4DRef<double> qrotj(&wk(iwk_indx),                maxdims, 10, 1, 1);
                            FortranArray4DRef<double> qrotk(&wk(iwk_indx + maxdims * 10), maxdims, 10, 1, 1);
                            FortranArray4DRef<double> qroti(&wk(iwk_indx + maxdims * 20), maxdims, 10, 1, 1);
                            bc2005_ns::bc2005(jdim, kdim, idim,
                                              wq, wqj0, wqk0, wqi0,
                                              ista, iend, jsta, jend, ksta, kend,
                                              nface, wxib, wtj0, wtk0, wti0,
                                              wvis, wvj0, wvk0, wvi0,
                                              mdim, ndim, bcdata, filname,
                                              qp, vp, tp,
                                              jdimp, kdimp, idimp,
                                              qrotj, qrotk, qroti,
                                              nbll, nblp,
                                              nou, bou, nbuf, ibufdim,
                                              myid, mblk2nd, maxbl, nummem);
                        }
                        isklton = iskltsav;
                    } // end nface 3/4

                    // i = constant interface (nface 1 or 2)
                    if (nface == 1 || nface == 2) {
                        if (isklton == 1) {
                            int kcheck = iwk_indx + maxdims * 9 - 1;
                            if (kcheck > nwork) {
                                nou(1) = std::min(nou(1) + 1, ibufdim);
                                std::snprintf(bou(nou(1), 1), 120,
                                    " stop in bc_period...insufficient wk storage for chkrot");
                                int m1 = -1;
                                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                            }
                            FortranArray3DRef<double> wx(&w(lx),   jdim, kdim, idim);
                            FortranArray3DRef<double> wy(&w(ly),   jdim, kdim, idim);
                            FortranArray3DRef<double> wz(&w(lz),   jdim, kdim, idim);
                            FortranArray3DRef<double> xp(&w(lw(10, nblp)), jdimp, kdimp, idimp);
                            FortranArray3DRef<double> yp(&w(lw(11, nblp)), jdimp, kdimp, idimp);
                            FortranArray3DRef<double> zp(&w(lw(12, nblp)), jdimp, kdimp, idimp);
                            FortranArray4DRef<double> bcdata(&w(ldata), mdim, ndim, 2, 12);
                            FortranArray3DRef<double> xyzjp(&wk(iwk_indx),               maxdims, 3, 1);
                            FortranArray3DRef<double> xyzkp(&wk(iwk_indx + maxdims * 3),  maxdims, 3, 1);
                            FortranArray3DRef<double> xyzip(&wk(iwk_indx + maxdims * 6),  maxdims, 3, 1);
                            chkrot_ns::chkrot(nbll, jdim, kdim, idim,
                                              wx, wy, wz,
                                              nblp, jdimp, kdimp, idimp,
                                              xp, yp, zp,
                                              nface, bcdata,
                                              xyzjp, xyzkp, xyzip,
                                              ista, iend, jsta, jend, ksta, kend,
                                              mdim, ndim, lcnt,
                                              xorig, yorig, zorig,
                                              maxbl, period_miss, lbcprd,
                                              nou, bou, nbuf, ibufdim, myid);
                        }
                        if (isklton == 1) {
                            int kcheck = iwk_indx + maxdims * 36 - 1;
                            if (kcheck > nwork) {
                                nou(1) = std::min(nou(1) + 1, ibufdim);
                                std::snprintf(bou(nou(1), 1), 120,
                                    " stop in bc_period...insufficient wk storage for bctype 2005");
                                int m1 = -1;
                                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                            }
                        }
                        int iskltsav = isklton;
                        isklton = 0;
                        {
                            FortranArray4DRef<double> wq   (&w(lq),   jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wqj0 (&w(lqj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wqk0 (&w(lqk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wqi0 (&w(lqi0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wxib (&w(lxib), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wtj0 (&w(ltj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wtk0 (&w(ltk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wti0 (&w(lti0), jdim, kdim, idim, 5);
                            FortranArray3DRef<double> wvis (&w(lvis), jdim, kdim, idim);
                            FortranArray4DRef<double> wvj0 (&w(lvj0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wvk0 (&w(lvk0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> wvi0 (&w(lvi0), jdim, kdim, idim, 5);
                            FortranArray4DRef<double> bcdata(&w(ldata), mdim, ndim, 2, 12);
                            FortranArray4DRef<double> qp   (&w(lw(1,  nblp)), jdimp, kdimp, idimp, 5);
                            FortranArray3DRef<double> vp   (&w(lw(13, nblp)), jdimp, kdimp, idimp);
                            FortranArray4DRef<double> tp   (&w(lw(19, nblp)), jdimp, kdimp, idimp, 2);
                            FortranArray4DRef<double> qrotj(&wk(iwk_indx),                maxdims, 10, 1, 1);
                            FortranArray4DRef<double> qrotk(&wk(iwk_indx + maxdims * 10), maxdims, 10, 1, 1);
                            FortranArray4DRef<double> qroti(&wk(iwk_indx + maxdims * 20), maxdims, 10, 1, 1);
                            bc2005_ns::bc2005(jdim, kdim, idim,
                                              wq, wqj0, wqk0, wqi0,
                                              ista, iend, jsta, jend, ksta, kend,
                                              nface, wxib, wtj0, wtk0, wti0,
                                              wvis, wvj0, wvk0, wvi0,
                                              mdim, ndim, bcdata, filname,
                                              qp, vp, tp,
                                              jdimp, kdimp, idimp,
                                              qrotj, qrotk, qroti,
                                              nbll, nblp,
                                              nou, bou, nbuf, ibufdim,
                                              myid, mblk2nd, maxbl, nummem);
                        }
                        isklton = iskltsav;
                    } // end nface 1/2

                } // end if iadvance(nbll) >= 0

            } // end if nd_srce == myid && nd_recv == myid

        } // end for lcnt

        // Second Case: data needed to set periodic bc lies on another
        //              processor
        // (no code in Fortran source for this case)

    } // end if ntime > 0 && nbcprd > 0

    return;
}

} // namespace bc_period_ns
