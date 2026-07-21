// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "force.h"
#include <cstdio>
#include "rp3d.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdint>

namespace force_ns {

void force(int& jdim, int& kdim, int& idim,
           FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
           FortranArray4DRef<double> sk, FortranArray4DRef<double> sj, FortranArray4DRef<double> si,
           double& cl, double& cd, double& cz, double& cy, double& cx,
           double& cmy, double& cmx, double& cmz, double& chd, double& swet,
           int& i00,
           FortranArray3DRef<double> ub, FortranArray3DRef<double> vb, FortranArray3DRef<double> wb,
           FortranArray3DRef<double> vmuk, FortranArray3DRef<double> vmuj, FortranArray3DRef<double> vmui,
           FortranArray3DRef<double> vol,
           int& ifo, int& jfo, int& kfo,
           FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci,
           FortranArray3DRef<double> blank,
           int& nbl,
           FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk, FortranArray4DRef<double> xtbi,
           int& iuns,
           FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
           FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
           FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim,
           FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
           int& nn, int& maxseg)
{
    // COMMON block aliases
    float& cdv    = cmn_drag.cdv;
    float& cdp    = cmn_drag.cdp;
    float& gamma  = cmn_fluid.gamma;
    float& xmc    = cmn_fsum.xmc;
    float& ymc    = cmn_fsum.ymc;
    float& zmc    = cmn_fsum.zmc;
    float& xmach  = cmn_info.xmach;
    float& alpha  = cmn_info.alpha;
    float& beta   = cmn_info.beta;
    int32_t& i2d  = cmn_twod.i2d;
    float& p0     = cmn_ivals.p0;
    float& reue   = cmn_reyue.reue;
    int32_t* ivisc = cmn_reyue.ivisc;  // 0-based C array, ivisc[0]=ivisc(1), etc.
    int32_t& ialph = cmn_igrdtyp.ialph;

    // Local variables
    double cpc, cnst;
    double cosa, sina, cosb, sinb;
    double xa, ya, za;
    double dcp, dcx, dcy, dcz;
    double dcxp, dcyp, dczp;
    double dcdpl, dcdvl;
    double urel, vrel, wrel;
    double tau, vnorm;
    double fact, dsw;
    double cxl, cyl, czl, cmyl, cmxl, cmzl, chdl, swetl, cdvl, cdpl;
    double cds, cls, cmys, cmxs, cmzs, swets;
    double sgn;
    double temp;
    double xas, yas, zas;
    int nseg, k, kc, j, jc, i, ic;
    int ist, ifn, jst, jfn, kst, kfn;
    int kk, jj, ii, ns;


    cpc   = 2.e0 / ((double)gamma * (double)xmach * (double)xmach);
    cnst  = 4.0 / ((double)reue * (double)xmach);

    cosa  = std::cos((double)alpha);
    sina  = std::sin((double)alpha);
    cosb  = std::cos((double)beta);
    sinb  = std::sin((double)beta);

    cl    = 0.e0;
    cd    = 0.e0;
    cz    = 0.e0;
    cy    = 0.e0;
    cx    = 0.e0;
    cmy   = 0.e0;
    cmx   = 0.e0;
    cmz   = 0.e0;
    chd   = 0.e0;
    swet  = 0.e0;
    cdp   = 0.e0;
    cdv   = 0.e0;


    //******************************************************************
    //     forces on k=constant surfaces
    //******************************************************************
    for (kk = 1; kk <= 2; kk++) {

        if (kk == 1) {
            nseg = nbck0(nbl);
            k    = 1;
            kc   = 1;
        } else {
            nseg = nbckdim(nbl);
            k    = kdim;
            kc   = kdim - 1;
        }

        for (ns = 1; ns <= nseg; ns++) {

            if (kbcinfo(nbl,ns,6,kk) > 0) {

                ist = kbcinfo(nbl,ns,2,kk);
                ifn = kbcinfo(nbl,ns,3,kk) - 1;
                jst = kbcinfo(nbl,ns,4,kk);
                jfn = kbcinfo(nbl,ns,5,kk) - 1;

                sgn = 1.0;
                if (kk > 1) sgn = -1.0;

                for (i = ist; i <= ifn; i++) {
                    cxl   = 0.e0;
                    cyl   = 0.e0;
                    czl   = 0.e0;
                    cmyl  = 0.e0;
                    cmxl  = 0.e0;
                    cmzl  = 0.e0;
                    chdl  = 0.e0;
                    xas   = 0.e0;
                    yas   = 0.e0;
                    zas   = 0.e0;
                    swetl = 0.e0;
                    cdvl  = 0.e0;
                    cdpl  = 0.e0;

                    for (j = jst; j <= jfn; j++) {

                        xa = .25e0 * (x(j,k,i) + x(j+1,k,i) + x(j,k,i+1) + x(j+1,k,i+1));
                        ya = .25e0 * (y(j,k,i) + y(j+1,k,i) + y(j,k,i+1) + y(j+1,k,i+1));
                        za = .25e0 * (z(j,k,i) + z(j+1,k,i) + z(j,k,i+1) + z(j+1,k,i+1));

                        dcp = -(qk0(j,i,5,kk+kk-1) / (double)p0 - 1.e0) * cpc * sk(j,k,i,4);
                        dcx = dcp * sk(j,k,i,1) * sgn;
                        dcy = dcp * sk(j,k,i,2) * sgn;
                        dcz = dcp * sk(j,k,i,3) * sgn;
                        // pressure drag
                        dcdpl = dcx*cosa*cosb - dcy*sinb + dcz*sina*cosb;

                        dcdvl = 0.;
                        if (ivisc[2] > 0 &&
                            (std::abs(kbcinfo(nbl,ns,1,kk)) == 2004 ||
                             std::abs(kbcinfo(nbl,ns,1,kk)) == 2014 ||
                             std::abs(kbcinfo(nbl,ns,1,kk)) == 2024 ||
                             std::abs(kbcinfo(nbl,ns,1,kk)) == 2034 ||
                             std::abs(kbcinfo(nbl,ns,1,kk)) == 2016)) {
                            dcxp  = dcx;
                            dcyp  = dcy;
                            dczp  = dcz;
                            urel  = ub(j,kc,i);
                            vrel  = vb(j,kc,i);
                            wrel  = wb(j,kc,i);
                            if (iuns > 0 && kbcinfo(nbl,ns,1,kk) > 0) {
                                urel = ub(j,kc,i) - xtbk(j,i,1,kk);
                                vrel = vb(j,kc,i) - xtbk(j,i,2,kk);
                                wrel = wb(j,kc,i) - xtbk(j,i,3,kk);
                            }
                            tau   = vmuk(j,i,kk) * cnst / vol(j,kc,i) * sk(j,k,i,4) * sk(j,k,i,4);
                            vnorm = (urel*sk(j,k,i,1) + vrel*sk(j,k,i,2)
                                     + wrel*sk(j,k,i,3)) * sgn;
                            dcx   = tau * (urel - vnorm*sk(j,k,i,1)*sgn);
                            dcz   = tau * (wrel - vnorm*sk(j,k,i,3)*sgn);
                            dcy   = tau * (vrel - vnorm*sk(j,k,i,2)*sgn);
                            // viscous drag
                            dcdvl = dcx*cosa*cosb - dcy*sinb + dcz*sina*cosb;
                            dcx = dcxp + dcx;
                            dcy = dcyp + dcy;
                            dcz = dczp + dcz;
                        }

                        // only use contributions from points with interface (solid surface)
                        // boundary conditions, and only from those points not blanked out
                        fact = bck(j,i,kk) * blank(j,kc,i);
                        dcx   = dcx   * fact;
                        dcy   = dcy   * fact;
                        dcz   = dcz   * fact;
                        dcdpl = dcdpl * fact;
                        dcdvl = dcdvl * fact;
                        dsw   = sk(j,k,i,4) * fact;

                        swetl = swetl + dsw;
                        cxl   = cxl + dcx;
                        cyl   = cyl + dcy;
                        czl   = czl + dcz;
                        cdpl  = cdpl + dcdpl;
                        cdvl  = cdvl + dcdvl;
                        cmyl  = cmyl - dcz*(xa - (double)xmc) + dcx*(za - (double)zmc);
                        cmxl  = cmxl + dcz*(ya - (double)ymc) - dcy*(za - (double)zmc);
                        cmzl  = cmzl + dcy*(xa - (double)xmc) - dcx*(ya - (double)ymc);

                    } // j loop

                    // integrated values
                    cds   = cxl*cosa*cosb - cyl*sinb + czl*sina*cosb;
                    cls   = -cxl*sina + czl*cosa;
                    cmys  = cmyl;
                    cmxs  = cmxl;
                    cmzs  = cmzl;
                    swets = swetl;
                    cl    = cl + cls;
                    cd    = cd + cds;
                    cz    = cz + czl;
                    cy    = cy + cyl;
                    cx    = cx + cxl;
                    cmy   = cmy + cmys;
                    cmx   = cmx + cmxs;
                    cmz   = cmz + cmzs;
                    chd   = 1.0;
                    swet  = swet + swets;
                    cdp   = (float)((double)cdp + cdpl);
                    cdv   = (float)((double)cdv + cdvl);

                } // i loop
            } // if kbcinfo > 0
        } // ns loop
    } // kk loop


    //******************************************************************
    //     forces on j=constant surfaces
    //******************************************************************

    for (jj = 1; jj <= 2; jj++) {

        if (jj == 1) {
            nseg = nbcj0(nbl);
            j    = 1;
            jc   = 1;
        } else {
            nseg = nbcjdim(nbl);
            j    = jdim;
            jc   = jdim - 1;
        }

        for (ns = 1; ns <= nseg; ns++) {

            if (jbcinfo(nbl,ns,6,jj) > 0) {

                ist = jbcinfo(nbl,ns,2,jj);
                ifn = jbcinfo(nbl,ns,3,jj) - 1;
                kst = jbcinfo(nbl,ns,4,jj);
                kfn = jbcinfo(nbl,ns,5,jj) - 1;

                sgn = 1.0;
                if (jj > 1) sgn = -1.0;

                for (i = ist; i <= ifn; i++) {
                    cxl   = 0.e0;
                    cyl   = 0.e0;
                    czl   = 0.e0;
                    cmyl  = 0.e0;
                    cmxl  = 0.e0;
                    cmzl  = 0.e0;
                    chdl  = 0.e0;
                    swetl = 0.e0;
                    cdvl  = 0.e0;
                    cdpl  = 0.e0;

                    for (k = kst; k <= kfn; k++) {

                        xa = .25e0 * (x(j,k,i) + x(j,k+1,i) + x(j,k,i+1) + x(j,k+1,i+1));
                        ya = .25e0 * (y(j,k,i) + y(j,k+1,i) + y(j,k,i+1) + y(j,k+1,i+1));
                        za = .25e0 * (z(j,k,i) + z(j,k+1,i) + z(j,k,i+1) + z(j,k+1,i+1));

                        dcp = -(qj0(k,i,5,jj+jj-1) / (double)p0 - 1.e0) * cpc * sj(j,k,i,4);
                        dcx = dcp * sj(j,k,i,1) * sgn;
                        dcy = dcp * sj(j,k,i,2) * sgn;
                        dcz = dcp * sj(j,k,i,3) * sgn;
                        // pressure drag
                        dcdpl = dcx*cosa*cosb - dcy*sinb + dcz*sina*cosb;

                        dcdvl = 0.;
                        if (ivisc[1] > 0 &&
                            (std::abs(jbcinfo(nbl,ns,1,jj)) == 2004 ||
                             std::abs(jbcinfo(nbl,ns,1,jj)) == 2014 ||
                             std::abs(jbcinfo(nbl,ns,1,jj)) == 2024 ||
                             std::abs(jbcinfo(nbl,ns,1,jj)) == 2034 ||
                             std::abs(jbcinfo(nbl,ns,1,jj)) == 2016)) {
                            dcxp  = dcx;
                            dcyp  = dcy;
                            dczp  = dcz;
                            urel  = ub(jc,k,i);
                            vrel  = vb(jc,k,i);
                            wrel  = wb(jc,k,i);
                            if (iuns > 0 && jbcinfo(nbl,ns,1,jj) > 0) {
                                urel = ub(jc,k,i) - xtbj(k,i,1,jj);
                                vrel = vb(jc,k,i) - xtbj(k,i,2,jj);
                                wrel = wb(jc,k,i) - xtbj(k,i,3,jj);
                            }
                            tau   = vmuj(k,i,jj) * cnst / vol(jc,k,i) * sj(j,k,i,4) * sj(j,k,i,4);
                            vnorm = (urel*sj(j,k,i,1) + vrel*sj(j,k,i,2)
                                     + wrel*sj(j,k,i,3)) * sgn;
                            dcx   = tau * (urel - vnorm*sj(j,k,i,1)*sgn);
                            dcz   = tau * (wrel - vnorm*sj(j,k,i,3)*sgn);
                            dcy   = tau * (vrel - vnorm*sj(j,k,i,2)*sgn);
                            // viscous drag
                            dcdvl = dcx*cosa*cosb - dcy*sinb + dcz*sina*cosb;
                            dcx = dcxp + dcx;
                            dcy = dcyp + dcy;
                            dcz = dczp + dcz;
                        }

                        // only use contributions from points with interface (solid surface)
                        // boundary conditions, and only from those points not blanked out
                        fact = bcj(k,i,jj) * blank(jc,k,i);
                        dcx   = dcx   * fact;
                        dcy   = dcy   * fact;
                        dcz   = dcz   * fact;
                        dcdpl = dcdpl * fact;
                        dcdvl = dcdvl * fact;
                        dsw   = sj(j,k,i,4) * fact;

                        swetl = swetl + dsw;
                        cxl   = cxl + dcx;
                        cyl   = cyl + dcy;
                        czl   = czl + dcz;
                        cdpl  = cdpl + dcdpl;
                        cdvl  = cdvl + dcdvl;
                        cmyl  = cmyl - dcz*(xa - (double)xmc) + dcx*(za - (double)zmc);
                        cmxl  = cmxl + dcz*(ya - (double)ymc) - dcy*(za - (double)zmc);
                        cmzl  = cmzl + dcy*(xa - (double)xmc) - dcx*(ya - (double)ymc);

                    } // k loop

                    // integrated values
                    cds   = cxl*cosa*cosb - cyl*sinb + czl*sina*cosb;
                    cls   = -cxl*sina + czl*cosa;
                    cmys  = cmyl;
                    cmxs  = cmxl;
                    cmzs  = cmzl;
                    swets = swetl;
                    cl    = cl + cls;
                    cd    = cd + cds;
                    cz    = cz + czl;
                    cy    = cy + cyl;
                    cx    = cx + cxl;
                    cmy   = cmy + cmys;
                    cmx   = cmx + cmxs;
                    cmz   = cmz + cmzs;
                    chd   = 1.0;
                    swet  = swet + swets;
                    cdp   = (float)((double)cdp + cdpl);
                    cdv   = (float)((double)cdv + cdvl);

                } // i loop
            } // if jbcinfo > 0
        } // ns loop
    } // jj loop


    //******************************************************************
    //     forces on i=constant surfaces
    //******************************************************************

    for (ii = 1; ii <= 2; ii++) {

        if (ii == 1) {
            nseg = nbci0(nbl);
            i    = 1;
            ic   = 1;
        } else {
            nseg = nbcidim(nbl);
            i    = idim;
            ic   = idim - 1;
        }

        for (ns = 1; ns <= nseg; ns++) {

            if (ibcinfo(nbl,ns,6,ii) > 0) {

                jst = ibcinfo(nbl,ns,2,ii);
                jfn = ibcinfo(nbl,ns,3,ii) - 1;
                kst = ibcinfo(nbl,ns,4,ii);
                kfn = ibcinfo(nbl,ns,5,ii) - 1;

                sgn = 1.0;
                if (ii > 1) sgn = -1.0;

                for (j = jst; j <= jfn; j++) {
                    cxl   = 0.e0;
                    cyl   = 0.e0;
                    czl   = 0.e0;
                    cmyl  = 0.e0;
                    cmxl  = 0.e0;
                    cmzl  = 0.e0;
                    chdl  = 0.e0;
                    swetl = 0.e0;
                    cdvl  = 0.e0;
                    cdpl  = 0.e0;

                    for (k = kst; k <= kfn; k++) {

                        xa = .25e0 * (x(j,k,i) + x(j,k+1,i) + x(j+1,k,i) + x(j+1,k+1,i));
                        ya = .25e0 * (y(j,k,i) + y(j,k+1,i) + y(j+1,k,i) + y(j+1,k+1,i));
                        za = .25e0 * (z(j,k,i) + z(j,k+1,i) + z(j+1,k,i) + z(j+1,k+1,i));

                        dcp = -(qi0(j,k,5,ii+ii-1) / (double)p0 - 1.e0) * cpc * si(j,k,i,4);
                        dcx = dcp * si(j,k,i,1) * sgn;
                        dcy = dcp * si(j,k,i,2) * sgn;
                        dcz = dcp * si(j,k,i,3) * sgn;
                        // pressure drag
                        dcdpl = dcx*cosa*cosb - dcy*sinb + dcz*sina*cosb;

                        dcdvl = 0.;
                        if (ivisc[0] > 0 &&
                            (std::abs(ibcinfo(nbl,ns,1,ii)) == 2004 ||
                             std::abs(ibcinfo(nbl,ns,1,ii)) == 2014 ||
                             std::abs(ibcinfo(nbl,ns,1,ii)) == 2024 ||
                             std::abs(ibcinfo(nbl,ns,1,ii)) == 2034 ||
                             std::abs(ibcinfo(nbl,ns,1,ii)) == 2016)) {
                            dcxp  = dcx;
                            dcyp  = dcy;
                            dczp  = dcz;
                            urel  = ub(j,k,ic);
                            vrel  = vb(j,k,ic);
                            wrel  = wb(j,k,ic);
                            if (iuns > 0 && ibcinfo(nbl,ns,1,ii) > 0) {
                                urel = ub(j,k,ic) - xtbi(j,k,1,ii);
                                vrel = vb(j,k,ic) - xtbi(j,k,2,ii);
                                wrel = wb(j,k,ic) - xtbi(j,k,3,ii);
                            }
                            tau   = vmui(j,k,ii) * cnst / vol(j,k,ic) * si(j,k,i,4) * si(j,k,i,4);
                            vnorm = (urel*si(j,k,i,1) + vrel*si(j,k,i,2)
                                     + wrel*si(j,k,i,3)) * sgn;
                            dcx   = tau * (urel - vnorm*si(j,k,i,1)*sgn);
                            dcz   = tau * (wrel - vnorm*si(j,k,i,3)*sgn);
                            dcy   = tau * (vrel - vnorm*si(j,k,i,2)*sgn);
                            // viscous drag
                            dcdvl = dcx*cosa*cosb - dcy*sinb + dcz*sina*cosb;
                            dcx = dcxp + dcx;
                            dcy = dcyp + dcy;
                            dcz = dczp + dcz;
                        }

                        // only use contributions from points with interface (solid surface)
                        // boundary conditions, and only from those points not blanked out
                        fact = bci(j,k,ii) * blank(j,k,ic);
                        dcx   = dcx   * fact;
                        dcy   = dcy   * fact;
                        dcz   = dcz   * fact;
                        dcdpl = dcdpl * fact;
                        dcdvl = dcdvl * fact;
                        dsw   = si(j,k,i,4) * fact;

                        swetl = swetl + dsw;
                        cxl   = cxl + dcx;
                        cyl   = cyl + dcy;
                        czl   = czl + dcz;
                        cdpl  = cdpl + dcdpl;
                        cdvl  = cdvl + dcdvl;
                        cmyl  = cmyl - dcz*(xa - (double)xmc) + dcx*(za - (double)zmc);
                        cmxl  = cmxl + dcz*(ya - (double)ymc) - dcy*(za - (double)zmc);
                        cmzl  = cmzl + dcy*(xa - (double)xmc) - dcx*(ya - (double)ymc);

                    } // k loop

                    // integrated values
                    cds   = cxl*cosa*cosb - cyl*sinb + czl*sina*cosb;
                    cls   = -cxl*sina + czl*cosa;
                    cmys  = cmyl;
                    cmxs  = cmxl;
                    cmzs  = cmzl;
                    swets = swetl;
                    cl    = cl + cls;
                    cd    = cd + cds;
                    cz    = cz + czl;
                    cy    = cy + cyl;
                    cx    = cx + cxl;
                    cmy   = cmy + cmys;
                    cmx   = cmx + cmxs;
                    cmz   = cmz + cmzs;
                    chd   = 1.0;
                    swet  = swet + swets;
                    cdp   = (float)((double)cdp + cdpl);
                    cdv   = (float)((double)cdv + cdvl);

                } // j loop
            } // if ibcinfo > 0
        } // ns loop
    } // ii loop

    // force cmx and cmy=0 if 2-D; then user doesn't have to worry about
    // inputting "correct" cell center location for the non-2-D-direction
    if (i2d == 1) {
        cmx = 0.;
        cmz = 0.;
    }

    // note: the forces above are computed in the standard cfl3d coordinate
    // system where z is "up". if ialph = 1, then must map these into a system
    // with y "up" so that when output, they are consistent with the input
    // grid. If x,y,z is the standard cfl3d system and x',y',z' is the
    // input grid system (with x'=x), then for ialph .ne. 0 the transform
    // z = y', y = -z' maps the input grid onto the standard cfl3d system
    // (see also subroutine rp3d). Conversely z' = -y, y' = z' maps the
    // cfl3d system back into the input system (this is what is used below,
    // where the ' has been dropped from the left hand side)
    if (ialph > 0) {
        temp = cmz;
        cmz  = -cmy;
        cmy  = temp;
        temp = cz;
        cz   = -cy;
        cy   = temp;
    }

    return;
}

void rp3d(FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
          int& jdim, int& kdim, int& idim, int& igrid, int& ialph, int& igeom_img, int& irr)
{
    rp3d_ns::rp3d(x, y, z, jdim, kdim, idim, igrid, ialph, igeom_img, irr);
}

} // namespace force_ns
