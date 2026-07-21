// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "resetg.h"
#include "runtime/fortran_io.h"
#include "lead.h"
#include "grdmove.h"
#include "rotateq.h"
#include "rotateq0.h"
#include "rotateqb.h"
#include "updateg.h"
#include "writ_buf.h"
#include "ccomplex.h"

#include <cmath>
#include <cstdio>
#include <cstring>

namespace resetg_ns {

void resetg(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& iupdat, int& iseqr, int& maxbl, int& maxgr, int& maxseg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& nblock, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz, FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr, FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0, FortranArray1DRef<double> time2, FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> iadvance, FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx, FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> dthxmx, int& iitot, FortranArray1DRef<int> iovrlp, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray1DRef<int> iipntsg, FortranArray2DRef<int> ibpntsg, FortranArray3DRef<double> qb, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> iresetb, int& nt)
{
    // COMMON block aliases
    int& jdim    = cmn_ginfo.jdim;

    int& kdim    = cmn_ginfo.kdim;
    int& idim    = cmn_ginfo.idim;
    int& lq      = cmn_ginfo.lq;
    int& lqj0    = cmn_ginfo.lqj0;
    int& lqk0    = cmn_ginfo.lqk0;
    int& lqi0    = cmn_ginfo.lqi0;
    int& lqc0    = cmn_ginfo.lqc0;
    int& lx      = cmn_ginfo.lx;
    int& ly      = cmn_ginfo.ly;
    int& lz      = cmn_ginfo.lz;
    int& lglobal = cmn_mgrd.lglobal;
    double dt    = (double)cmn_info.dt;
    int& ita     = cmn_unst.ita;
    float& radtodeg = cmn_conversion.radtodeg;
    // /fsum/ aliases
    float& xmc   = cmn_fsum.xmc;
    float& ymc   = cmn_fsum.ymc;
    float& zmc   = cmn_fsum.zmc;

    // /motionmc/ aliases
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
    float& dxmxmc      = cmn_motionmc.dxmxmc;
    float& dymxmc      = cmn_motionmc.dymxmc;
    float& dzmxmc      = cmn_motionmc.dzmxmc;
    float& dthxmxmc    = cmn_motionmc.dthxmxmc;
    float& dthymxmc    = cmn_motionmc.dthymxmc;
    float& dthzmxmc    = cmn_motionmc.dthzmxmc;
    int32_t& itransmc  = cmn_motionmc.itransmc;
    int32_t& irotatmc  = cmn_motionmc.irotatmc;
    float& time2mc     = cmn_motionmc.time2mc;
    // Local variables
    double adxmx, adymx, adzmx, adthxmx, adthymx, adthzmx;
    double axorig, ayorig, azorig, athetx, athety, athetz;
    double t2max, t2maxmc;
    double xorg, yorg, zorg;
    double thtx_loc, thty_loc, thtz_loc;
    double thx, thy, thz, thxmx, thymx, thzmx;
    double dthtx, dthty, dthtz;
    double thtxmc, thtymc, thtzmc;
    double dthtxmc, dthtymc, dthtzmc;
    double xml, yml, zml, ca, sa;
    int ireset, iresetmc;
    int iuns;
    int lqwk;
    int irigb0, irbtrim0;



    // do 5 nbl = 1,nblock
    for (int nbl = 1; nbl <= nblock; nbl++) {
        iresetb(nbl) = 0;
        ireset = 0;
        if (iadvance(nbl) < 0) continue;  // go to 5
        iuns = std::max(irotat(nbl), itrans(nbl));
        if (levelg(nbl) >= lglobal && iuns > 0) {

            // check for need to reset
            lead_ns::lead(nbl, lw, lw2, maxbl);

            {
                double tmp;
                tmp = dxmx(nbl);   adxmx   = ccomplex_ns::ccabs(tmp);
                tmp = dymx(nbl);   adymx   = ccomplex_ns::ccabs(tmp);
                tmp = dzmx(nbl);   adzmx   = ccomplex_ns::ccabs(tmp);
                tmp = dthxmx(nbl); adthxmx = ccomplex_ns::ccabs(tmp);
                tmp = dthymx(nbl); adthymx = ccomplex_ns::ccabs(tmp);
                tmp = dthzmx(nbl); adthzmx = ccomplex_ns::ccabs(tmp);
                tmp = xorig(nbl) - xorig0(nbl); axorig = ccomplex_ns::ccabs(tmp);
                tmp = yorig(nbl) - yorig0(nbl); ayorig = ccomplex_ns::ccabs(tmp);
                tmp = zorig(nbl) - zorig0(nbl); azorig = ccomplex_ns::ccabs(tmp);
                tmp = thetax(nbl); athetx = ccomplex_ns::ccabs(tmp);
                tmp = thetay(nbl); athety = ccomplex_ns::ccabs(tmp);
                tmp = thetaz(nbl); athetz = ccomplex_ns::ccabs(tmp);
            }

            if (itrans(nbl) == 1) {
                if ((float)adxmx > 0.f && (float)axorig > (float)adxmx) ireset = 1;
                if ((float)adymx > 0.f && (float)ayorig > (float)adymx) ireset = 1;
                if ((float)adzmx > 0.f && (float)azorig > (float)adzmx) ireset = 1;
            } else {
                if ((float)adxmx > 0.f && (float)axorig > (float)adxmx) ireset = -1;
                if ((float)adymx > 0.f && (float)ayorig > (float)adymx) ireset = -1;
                if ((float)adzmx > 0.f && (float)azorig > (float)adzmx) ireset = -1;
            }
            if (irotat(nbl) == 1) {
                if ((float)adthxmx > 0.f && (float)athetx > (float)adthxmx) ireset = 1;
                if ((float)adthymx > 0.f && (float)athety > (float)adthymx) ireset = 1;
                if ((float)adthzmx > 0.f && (float)athetz > (float)adthzmx) ireset = 1;
            } else {
                if ((float)adthxmx > 0.f && (float)athetx > (float)adthxmx) ireset = -1;
                if ((float)adthymx > 0.f && (float)athety > (float)adthymx) ireset = -1;
                if ((float)adthzmx > 0.f && (float)athetz > (float)adthzmx) ireset = -1;
            }

            if (ireset >= 0) iresetb(nbl) = ireset;
            if (ireset < 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, " ");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                // format 99
                std::snprintf(bou(nou(1), 1), 120,
                    " WARNING: block%4d will not be reset - resetting allowed only if itrans/irotat = 1",
                    nbl);
            }

            // reset block nbl
            if (ireset > 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, " ");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                // format 100
                std::snprintf(bou(nou(1), 1), 120,
                    "resetting position of block%4d (grid%4d)",
                    nbl, igridg(nbl));

                // reset translation
                if (itrans(nbl) > 0) {
                    nou(1) = std::min(nou(1) + 1, ibufdim);
                    // format 101
                    std::snprintf(bou(nou(1), 1), 120,
                        "   current x,y,z displacements:        %8.3f  %8.3f  %8.3f",
                        (float)xorig(nbl), (float)yorig(nbl), (float)zorig(nbl));
                    nou(1) = std::min(nou(1) + 1, ibufdim);
                    // format 102
                    std::snprintf(bou(nou(1), 1), 120,
                        "   max. allowable x,y,z displacements: %8.3f  %8.3f  %8.3f",
                        (float)dxmx(nbl), (float)dymx(nbl), (float)dzmx(nbl));
                    xorg = xorig(nbl);
                    yorg = yorig(nbl);
                    zorg = zorig(nbl);
                    if (std::abs((float)utrans(nbl)) > 0.f) {
                        double tmp_u  = utrans(nbl);
                        double tmp_dx = dxmx(nbl);
                        t2max       = ccomplex_ns::ccabs(tmp_dx) / ccomplex_ns::ccabs(tmp_u);
                        time2(nbl)  = time2(nbl) - t2max - dt;
                        xorig(nbl)  = xorig0(nbl) + time2(nbl) * utrans(nbl);
                    }
                    if (std::abs((float)vtrans(nbl)) > 0.f) {
                        double tmp_v  = vtrans(nbl);
                        double tmp_dy = dymx(nbl);
                        t2max       = ccomplex_ns::ccabs(tmp_dy) / ccomplex_ns::ccabs(tmp_v);
                        time2(nbl)  = time2(nbl) - t2max - dt;
                        yorig(nbl)  = yorig0(nbl) + time2(nbl) * vtrans(nbl);
                    }
                    if (std::abs((float)wtrans(nbl)) > 0.f) {
                        double tmp_w  = wtrans(nbl);
                        double tmp_dz = dzmx(nbl);
                        t2max       = ccomplex_ns::ccabs(tmp_dz) / ccomplex_ns::ccabs(tmp_w);
                        time2(nbl)  = time2(nbl) - t2max - dt;
                        zorig(nbl)  = zorig0(nbl) + time2(nbl) * wtrans(nbl);
                    }
                    {
                        FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                        grdmove_ns::grdmove(nbl, jdim, kdim, idim,
                                            wx, wy, wz,
                                            xorg, yorg, zorg,
                                            xorig(nbl), yorig(nbl), zorig(nbl),
                                            thetax(nbl), thetay(nbl), thetaz(nbl));
                    }
                }


                // reset rotation
                if (irotat(nbl) > 0) {
                    thx   = thetax(nbl) * (double)radtodeg;
                    thy   = thetay(nbl) * (double)radtodeg;
                    thz   = thetaz(nbl) * (double)radtodeg;
                    thxmx = dthxmx(nbl) * (double)radtodeg;
                    thymx = dthymx(nbl) * (double)radtodeg;
                    thzmx = dthzmx(nbl) * (double)radtodeg;
                    nou(1) = std::min(nou(1) + 1, ibufdim);
                    // format 103
                    std::snprintf(bou(nou(1), 1), 120,
                        "   current x,y,z rotational displacements:        %8.3f  %8.3f  %8.3f",
                        (float)thx, (float)thy, (float)thz);
                    nou(1) = std::min(nou(1) + 1, ibufdim);
                    // format 104
                    std::snprintf(bou(nou(1), 1), 120,
                        "   max. allowable x,y,z rotational displacements: %8.3f  %8.3f  %8.3f",
                        (float)thxmx, (float)thymx, (float)thzmx);
                    xorg = xorig(nbl);
                    yorg = yorig(nbl);
                    zorg = zorig(nbl);
                    thtx_loc = thetax(nbl);
                    thty_loc = thetay(nbl);
                    thtz_loc = thetaz(nbl);
                    if (std::abs((float)omegax(nbl)) > 0.f) {
                        double tmp_ox   = omegax(nbl);
                        double tmp_dthx = dthxmx(nbl);
                        t2max        = ccomplex_ns::ccabs(tmp_dthx) / ccomplex_ns::ccabs(tmp_ox);
                        time2(nbl)   = time2(nbl) - t2max - dt;
                        thetax(nbl)  = time2(nbl) * omegax(nbl);
                    }
                    if (std::abs((float)omegay(nbl)) > 0.f) {
                        double tmp_oy   = omegay(nbl);
                        double tmp_dthy = dthymx(nbl);
                        t2max        = ccomplex_ns::ccabs(tmp_dthy) / ccomplex_ns::ccabs(tmp_oy);
                        time2(nbl)   = time2(nbl) - t2max - dt;
                        thetay(nbl)  = time2(nbl) * omegay(nbl);
                    }
                    if (std::abs((float)omegaz(nbl)) > 0.f) {
                        double tmp_oz   = omegaz(nbl);
                        double tmp_dthz = dthzmx(nbl);
                        t2max        = ccomplex_ns::ccabs(tmp_dthz) / ccomplex_ns::ccabs(tmp_oz);
                        time2(nbl)   = time2(nbl) - t2max - dt;
                        thetaz(nbl)  = time2(nbl) * omegaz(nbl);
                    }
                    dthtx = thetax(nbl) - thtx_loc;
                    dthty = thetay(nbl) - thty_loc;
                    dthtz = thetaz(nbl) - thtz_loc;
                    {
                        FortranArray3DRef<double> wx(&w(lx), jdim, kdim, idim);
                        FortranArray3DRef<double> wy(&w(ly), jdim, kdim, idim);
                        FortranArray3DRef<double> wz(&w(lz), jdim, kdim, idim);
                        grdmove_ns::grdmove(nbl, jdim, kdim, idim,
                                            wx, wy, wz,
                                            xorg, yorg, zorg,
                                            xorig(nbl), yorig(nbl), zorig(nbl),
                                            dthtx, dthty, dthtz);
                    }
                    // We rotated the x,y,zs thru an extra delta t, since we
                    // will be updating them again by delta t at the end of this
                    // routine (in the call to updateg).  For the Qs, however,
                    // we don't want to rotate them by extra amount:
                    dthtx = dthtx + dt * omegax(nbl);
                    dthty = dthty + dt * omegay(nbl);
                    dthtz = dthtz + dt * omegaz(nbl);

                    // interior values
                    lqwk = 1;
                    {
                        FortranArray4DRef<double> wq(&w(lq), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> wqrot(&wk(lqwk), jdim, kdim, idim, 5);
                        int ista = 1, iend = idim, jsta = 1, jend = jdim, ksta = 1, kend = kdim;
                        rotateq_ns::rotateq(jdim, kdim, idim, wq, wqrot,
                                            ista, iend, jsta, jend, ksta, kend,
                                            dthtx, dthty, dthtz);
                    }
                    for (int l = 1; l <= jdim * kdim * idim * 5; l++) {
                        w(l + lq - 1) = wk(l + lqwk - 1);
                    }
                    // qj0 values
                    {
                        int idim_m1 = idim - 1;
                        FortranArray4DRef<double> wqj0(&w(lqj0), kdim, idim_m1, 5, 4);
                        FortranArray4DRef<double> wqj0rot(&wk(lqwk), kdim, idim_m1, 5, 4);
                        int lsta = 1, lend = kdim, msta = 1, mend = idim_m1;
                        rotateq0_ns::rotateq0(kdim, idim_m1, wqj0, wqj0rot,
                                              lsta, lend, msta, mend,
                                              dthtx, dthty, dthtz);
                    }
                    for (int l = 1; l <= kdim * (idim - 1) * 5 * 4; l++) {
                        w(l + lqj0 - 1) = wk(l + lqwk - 1);
                    }
                    // qk0 values
                    {
                        int idim_m1 = idim - 1;
                        FortranArray4DRef<double> wqk0(&w(lqk0), jdim, idim_m1, 5, 4);
                        FortranArray4DRef<double> wqk0rot(&wk(lqwk), jdim, idim_m1, 5, 4);
                        int lsta = 1, lend = jdim, msta = 1, mend = idim_m1;
                        rotateq0_ns::rotateq0(jdim, idim_m1, wqk0, wqk0rot,
                                              lsta, lend, msta, mend,
                                              dthtx, dthty, dthtz);
                    }
                    for (int l = 1; l <= jdim * (idim - 1) * 5 * 4; l++) {
                        w(l + lqk0 - 1) = wk(l + lqwk - 1);
                    }
                    // qi0 values
                    {
                        FortranArray4DRef<double> wqi0(&w(lqi0), jdim, kdim, 5, 4);
                        FortranArray4DRef<double> wqi0rot(&wk(lqwk), jdim, kdim, 5, 4);
                        int lsta = 1, lend = jdim, msta = 1, mend = kdim;
                        rotateq0_ns::rotateq0(jdim, kdim, wqi0, wqi0rot,
                                              lsta, lend, msta, mend,
                                              dthtx, dthty, dthtz);
                    }
                    for (int l = 1; l <= jdim * kdim * 5 * 4; l++) {
                        w(l + lqi0 - 1) = wk(l + lqwk - 1);
                    }
                    // qc0 values for second-order time advancement
                    if (std::abs(ita) == 2) {
                        int idim_m1 = idim - 1;
                        FortranArray4DRef<double> wqc0(&w(lqc0), jdim, kdim, idim_m1, 5);
                        FortranArray4DRef<double> wqc0rot(&wk(lqwk), jdim, kdim, idim_m1, 5);
                        int ista = 1, iend = idim_m1, jsta = 1, jend = jdim, ksta = 1, kend = kdim;
                        rotateq_ns::rotateq(jdim, kdim, idim_m1, wqc0, wqc0rot,
                                            ista, iend, jsta, jend, ksta, kend,
                                            dthtx, dthty, dthtz);
                        for (int l = 1; l <= jdim * kdim * (idim - 1) * 5; l++) {
                            w(l + lqc0 - 1) = wk(l + lqwk - 1);
                        }
                    }
                    // qb values for chimera scheme
                    if (iovrlp(nbl) == 1) {
                        rotateqb_ns::rotateqb(nbl, dthtx, dthty, dthtz,
                                              maxbl, iitot, ibcg,
                                              lig, lbg, ibpntsg, iipntsg, qb);
                    }
                } // end if (irotat(nbl) > 0)

                nou(1) = std::min(nou(1) + 1, ibufdim);
                // format 105: ' '
                std::snprintf(bou(nou(1), 1), 120, " ");

            } // end if (ireset > 0)

            int iunit11 = 11;
            writ_buf_ns::writ_buf(nbl, iunit11, nou, bou, nbuf, ibufdim,
                                  myhost, myid, mycomm, mblk2nd, maxbl);

        } // end if (levelg(nbl) >= lglobal && iuns > 0)
    } // end do 5 nbl = 1,nblock


    // check to see if moment center needs to be reset
    iresetmc = 0;

    if (itransmc == 1 || irotatmc == 1) {
        {
            double tmp;
            tmp = (double)dxmxmc;   adxmx   = ccomplex_ns::ccabs(tmp);
            tmp = (double)dymxmc;   adymx   = ccomplex_ns::ccabs(tmp);
            tmp = (double)dzmxmc;   adzmx   = ccomplex_ns::ccabs(tmp);
            tmp = (double)dthxmxmc; adthxmx = ccomplex_ns::ccabs(tmp);
            tmp = (double)dthymxmc; adthymx = ccomplex_ns::ccabs(tmp);
            tmp = (double)dthzmxmc; adthzmx = ccomplex_ns::ccabs(tmp);
            tmp = (double)xorigmc;  axorig  = ccomplex_ns::ccabs(tmp);
            tmp = (double)yorigmc;  ayorig  = ccomplex_ns::ccabs(tmp);
            tmp = (double)zorigmc;  azorig  = ccomplex_ns::ccabs(tmp);
            tmp = (double)thetaxmc; athetx  = ccomplex_ns::ccabs(tmp);
            tmp = (double)thetaymc; athety  = ccomplex_ns::ccabs(tmp);
            tmp = (double)thetazmc; athetz  = ccomplex_ns::ccabs(tmp);
        }
        if ((float)adxmx > 0.f && (float)axorig > (float)adxmx)   iresetmc = 1;
        if ((float)adymx > 0.f && (float)ayorig > (float)adymx)   iresetmc = 1;
        if ((float)adzmx > 0.f && (float)azorig > (float)adzmx)   iresetmc = 1;
        if ((float)adthxmx > 0.f && (float)athetx > (float)adthxmx) iresetmc = 1;
        if ((float)adthymx > 0.f && (float)athety > (float)adthymx) iresetmc = 1;
        if ((float)adthzmx > 0.f && (float)athetz > (float)adthzmx) iresetmc = 1;

        if (iresetmc > 0) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, "\n");
            // format 106: "resetting moment center"
            fprintf(f11, "resetting moment center\n");

            double xorg_mc = 0.0, yorg_mc = 0.0, zorg_mc = 0.0;

            if (itransmc > 0) {
                // format 101
                fprintf(f11,
                    "   current x,y,z displacements:        %8.3f  %8.3f  %8.3f\n",
                    (float)xorigmc, (float)yorigmc, (float)zorigmc);
                // format 102
                fprintf(f11,
                    "   max. allowable x,y,z displacements: %8.3f  %8.3f  %8.3f\n",
                    (float)dxmxmc, (float)dymxmc, (float)dzmxmc);
                if (std::abs((float)utransmc) > 0.f) {
                    double tmp_dx = (double)dxmxmc;
                    double tmp_u  = (double)utransmc;
                    t2maxmc  = ccomplex_ns::ccabs(tmp_dx) / ccomplex_ns::ccabs(tmp_u);
                    time2mc  = (float)((double)time2mc - t2maxmc - dt);
                    xorigmc  = (float)((double)xorig0mc + (double)time2mc * (double)utransmc);
                    xmc      = (float)((double)cmn_motionmc.xmc0 + (double)time2mc * (double)utransmc);
                }
                if (std::abs((float)vtransmc) > 0.f) {
                    double tmp_dy = (double)dymxmc;
                    double tmp_v  = (double)vtransmc;
                    t2maxmc  = ccomplex_ns::ccabs(tmp_dy) / ccomplex_ns::ccabs(tmp_v);
                    time2mc  = (float)((double)time2mc - t2maxmc - dt);
                    yorigmc  = (float)((double)yorig0mc + (double)time2mc * (double)vtransmc);
                    ymc      = (float)((double)cmn_motionmc.ymc0 + (double)time2mc * (double)vtransmc);
                }
                if (std::abs((float)wtransmc) > 0.f) {
                    double tmp_dz = (double)dzmxmc;
                    double tmp_w  = (double)wtransmc;
                    t2maxmc  = ccomplex_ns::ccabs(tmp_dz) / ccomplex_ns::ccabs(tmp_w);
                    time2mc  = (float)((double)time2mc - t2maxmc - dt);
                    zorigmc  = (float)((double)zorig0mc + (double)time2mc * (double)wtransmc);
                    zmc      = (float)((double)cmn_motionmc.zmc0 + (double)time2mc * (double)wtransmc);
                }
                xorg_mc = (double)xorigmc;
                yorg_mc = (double)yorigmc;
                zorg_mc = (double)zorigmc;
            }

            if (irotatmc > 0) {
                thx   = (double)thetaxmc * (double)radtodeg;
                thy   = (double)thetaymc * (double)radtodeg;
                thz   = (double)thetazmc * (double)radtodeg;
                thxmx = (double)dthxmxmc * (double)radtodeg;
                thymx = (double)dthymxmc * (double)radtodeg;
                thzmx = (double)dthzmxmc * (double)radtodeg;
                thtxmc = (double)thetaxmc;
                thtymc = (double)thetaymc;
                thtzmc = (double)thetazmc;
                // format 103
                fprintf(f11,
                    "   current x,y,z rotational displacements:        %8.3f  %8.3f  %8.3f\n",
                    (float)thx, (float)thy, (float)thz);
                // format 104
                fprintf(f11,
                    "   max. allowable x,y,z rotational displacements: %8.3f  %8.3f  %8.3f\n",
                    (float)thxmx, (float)thymx, (float)thzmx);

                if (std::abs((float)omegaxmc) > 0.f) {
                    double tmp_dthx = (double)dthxmxmc;
                    double tmp_ox   = (double)omegaxmc;
                    t2maxmc  = ccomplex_ns::ccabs(tmp_dthx) / ccomplex_ns::ccabs(tmp_ox);
                    time2mc  = (float)((double)time2mc - t2maxmc - dt);
                    thetaxmc = (float)((double)time2mc * (double)omegaxmc);
                    dthtxmc  = (double)thetaxmc - thtxmc;
                    yml = (double)ymc;
                    zml = (double)zmc;
                    ca  = std::cos(dthtxmc);
                    sa  = std::sin(dthtxmc);
                    ymc = (float)((yml - (double)yorigmc) * ca - (zml - (double)zorigmc) * sa + yorg_mc);
                    zmc = (float)((yml - (double)yorigmc) * sa + (zml - (double)zorigmc) * ca + zorg_mc);
                }
                if (std::abs((float)omegaymc) > 0.f) {
                    double tmp_dthy = (double)dthymxmc;
                    double tmp_oy   = (double)omegaymc;
                    t2maxmc  = ccomplex_ns::ccabs(tmp_dthy) / ccomplex_ns::ccabs(tmp_oy);
                    time2mc  = (float)((double)time2mc - t2maxmc - dt);
                    thetaymc = (float)((double)time2mc * (double)omegaymc);
                    dthtymc  = (double)thetaymc - thtymc;
                    xml = (double)xmc;
                    zml = (double)zmc;
                    ca  = std::cos(dthtymc);
                    sa  = std::sin(dthtymc);
                    xmc = (float)( (xml - (double)xorigmc) * ca + (zml - (double)zorigmc) * sa + xorg_mc);
                    zmc = (float)(-(xml - (double)xorigmc) * sa + (zml - (double)zorigmc) * ca + zorg_mc);
                }
                if (std::abs((float)omegazmc) > 0.f) {
                    double tmp_dthz = (double)dthzmxmc;
                    double tmp_oz   = (double)omegazmc;
                    t2maxmc  = ccomplex_ns::ccabs(tmp_dthz) / ccomplex_ns::ccabs(tmp_oz);
                    time2mc  = (float)((double)time2mc - t2maxmc - dt);
                    thetazmc = (float)((double)time2mc * (double)omegazmc);
                    dthtzmc  = (double)thetazmc - thtzmc;
                    xml = (double)xmc;
                    yml = (double)ymc;
                    ca  = std::cos(dthtzmc);
                    sa  = std::sin(dthtzmc);
                    xmc = (float)((xml - (double)xorigmc) * ca - (yml - (double)yorigmc) * sa + xorg_mc);
                    ymc = (float)((xml - (double)xorigmc) * sa + (yml - (double)yorigmc) * ca + yorg_mc);
                }
            }

            // format 105: ' '
            fprintf(f11, " \n");

        } // end if (iresetmc > 0)

    } // end if (itransmc == 1 || irotatmc == 1)

    // update time, grid and metrics
    ireset = 0;
    for (int nbl = 1; nbl <= nblock; nbl++) {
        if (iresetb(nbl) > 0) {
            time2(nbl) = time2(nbl) + dt;
            ireset = 1;
        }
    }
    if (iresetmc > 0) {
        time2mc = (float)((double)time2mc + dt);
        ireset = 1;
    }

    if (ireset > 0) {
        irigb0   = 0;
        irbtrim0 = 0;
        updateg_ns::updateg(lw, lw2, w, mgwk, wk, nwork, iupdat, iseqr, maxbl,
                            maxgr, maxseg, nbci0, nbcj0, nbck0, nbcidim,
                            nbcjdim, nbckdim, ibcinfo, jbcinfo, kbcinfo,
                            nblock, levelg, igridg, utrans, vtrans, wtrans,
                            omegax, omegay, omegaz, xorig, yorig, zorig,
                            thetax, thetay, thetaz, rfreqt, rfreqr, xorig0,
                            yorig0, zorig0, time2, thetaxl, thetayl, thetazl,
                            itrans, irotat, idefrm, ncgg, iadvance,
                            nou, bou, nbuf, ibufdim, myid, myhost, mycomm,
                            mblk2nd, irigb0, irbtrim0, nt);
    }

} // end resetg

} // namespace resetg_ns
