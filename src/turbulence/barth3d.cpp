// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "barth3d.h"
#include "triv.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstring>


namespace barth3d_ns {

static inline double ccabs_bb(double x) { return std::fabs(x); }
static inline double ccmincr_bb(double a, double b) { return std::min(a, b); }
static inline double ccmaxcr_bb(double a, double b) { return std::max(a, b); }
static inline double ccsignrc_bb(double a, double b) { return (b >= 0.0) ? std::fabs(a) : -std::fabs(a); }

void barth3d(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, FortranArray3DRef<double> dtj, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vor, FortranArray3DRef<double> smin, FortranArray3DRef<double> xjb, FortranArray4DRef<double> tursav, FortranArray3DRef<double> xkb, FortranArray3DRef<double> turre, FortranArray3DRef<double> damp1, FortranArray3DRef<double> damp2, FortranArray3DRef<double> timestp, FortranArray3DRef<double> fnu, FortranArray2DRef<double> bx, FortranArray2DRef<double> cx, FortranArray2DRef<double> dx, FortranArray2DRef<double> fx, FortranArray2DRef<double> workx, FortranArray2DRef<double> by, FortranArray2DRef<double> cy, FortranArray2DRef<double> dy, FortranArray2DRef<double> fy, FortranArray2DRef<double> worky, FortranArray2DRef<double> bz, FortranArray2DRef<double> cz, FortranArray2DRef<double> dz, FortranArray2DRef<double> fz, FortranArray2DRef<double> workz, int& ntime, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, int& nbl, FortranArray3DRef<double> blnum, FortranArray3DRef<double> blank, int& iover, double& sumn1, double& sumn2, int& negn1, int& negn2, FortranArray4DRef<double> tursav2, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iex, int& iex2, int& nummem)
{
    // COMMON block aliases
    int& isklton         = cmn_sklton.isklton;
    int& icyc            = cmn_mgrd.icyc;
    int& nfreeze         = cmn_turbconv.nfreeze;
    int& ncyc            = cmn_mgrd.ncyc;
    int& ita             = cmn_unst.ita;
    float& dt            = cmn_info.dt;
    int& itaturb         = cmn_turbconv.itaturb;
    float* cflturb       = cmn_turbconv.cflturb;  // 0-based C array
    float& reue          = cmn_reyue.reue;
    float& xmach         = cmn_info.xmach;
    float& cbar          = cmn_fluid2.cbar;
    float& tinf          = cmn_reyue.tinf;
    float& gamma_f       = cmn_fluid.gamma;
    int& iturbord        = cmn_turbconv.iturbord;
    int& iaxi2planeturb  = cmn_axisym.iaxi2planeturb;
    int& istrongturbdis  = cmn_axisym.istrongturbdis;
    int& i2d             = cmn_twod.i2d;
    int& ilamlo          = cmn_lam.ilamlo;
    int& ilamhi          = cmn_lam.ilamhi;
    int& jlamlo          = cmn_lam.jlamlo;
    int& jlamhi          = cmn_lam.jlamhi;
    int& klamlo          = cmn_lam.klamlo;
    int& klamhi          = cmn_lam.klamhi;
    int& i_lam_forcezero = cmn_lam.i_lam_forcezero;
    float* tur10         = cmn_ivals.tur10;  // 0-based C array

    // Local variables
    int nsubit, nnit, iwrite, jd2, nss;
    int i, j, k, il, iu, jl, ju, kl, ku;
    int nblb, jbb, kbb, ibb;
    int not_;
    int negn;
    double phi, factor, re, c2b, c2bp;
    double akarman, cmu, c1e, c2e, ratt, sige, sigr, cr1, cr2, cmax, b1;
    double aplus1, aplus2;
    double tt, ypls, ssw, wnu, ra, exp1, exp2, ddy, dd, sdd;
    double xp, yp, zp, xm, ym, zm, xa, ya, za;
    double ttpo, ttmo, ttpn, ttmn, ttp, ttm;
    double cnud, cdp, cdm, trep, trem, cap, cam;
    double byy, cyy, dyy, bxx, cxx, dxx, bzz, czz, dzz;
    double volku, volkl, volju, voljl, voliu, volil;
    double xc, yc, zc, tc, uu, sgnu, app, apm;
    double cutoff, fact, sumn;


    // Lines 89-126: isklton messages and freeze check
    if (isklton > 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   computing turbulent viscosity using Baldwin-Barth, block =%4d", nbl);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     Freestream tur10 = %19.8e", (double)tur10[0]);
        if (iturbord == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     1st order advection on RHS");
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     2nd order advection on RHS");
        }
        if (iaxi2planeturb == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     BB model ignoring i-dir");
        }
        if (istrongturbdis == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     strong conserv - diss terms");
        }
    }
    if (icyc <= nfreeze) {
        if (isklton > 0) {
            nss = std::min(ncyc, nfreeze);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " ");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " turbulence model is frozen for %5d iterations or subits", nss);
        }
        sumn1 = 0.;
        sumn2 = 0.;
        negn1 = 0;
        negn2 = 0;
        return;
    }
    phi = 0.;
    if ((float)dt > 0.) {
        if (std::abs(ita) == 2) {
            phi = 0.5;
        } else {
            phi = 0.;
        }
        // revert to old way (always 1st order for turb model) if itaturb=0
        if (itaturb == 0) {
            phi = 0.;
            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   turb model is 1st order in time");
            }
        } else {
            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   turb model is same order in time as mean flow eqns");
            }
        }
    }

    // Lines 156-194: timestep setup
    nsubit = 1;
    factor = 10.;
    if ((float)cflturb[0] != 0.) {
        factor = (double)cflturb[0];
    }
    if ((float)dt < 0.) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            timestp(j,k,i) = factor*vol(j,k,i)/dtj(j,k,i);
            timestp(j,k,i) = ccmincr_bb(timestp(j,k,i), 100.);
        }
    } else {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            timestp(j,k,i) = (double)dt;
        }
    }

    // Lines 197-227: constants and fnu computation
    akarman = .41;
    cmu     = .09;
    c1e     = 1.20;
    c2e     = 2.00;
    ratt    = c1e/c2e;
    sige    = akarman*akarman/((c2e-c1e)*std::sqrt(cmu));
    sigr    = sige;
    cr1     = 1.0;
    cr2     = .01;
    cmax    = 1.0;
    b1      = 1.0;
    aplus1  = 26.;
    aplus2  = 10.;
    jd2     = (jdim-1)/2;
    re      = (double)reue/(double)xmach;
    c2b     = (double)cbar/(double)tinf;
    c2bp    = c2b + 1.0;

    iwrite = 0;
    nnit   = nsubit;

    // Get laminar viscosity (divided by rho) at cell centers
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++)
    for (j = 1; j <= jdim-1; j++) {
        tt = (double)gamma_f * q(j,k,i,5) / q(j,k,i,1);
        fnu(j,k,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt)/q(j,k,i,1);
    }

    // Lines 230-266: damping factors
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++)
    for (j = 1; j <= jdim-1; j++) {
        nblb = (int)(blnum(j,k,i) + 0.1);
        if (nblb != nbl) {
            ypls = 1000.;
        } else {
            jbb = (int)(xjb(j,k,i) + 0.1);
            kbb = (int)(xkb(j,k,i) + 0.1);
            ibb = (int)(tursav(j,k,i,2) + 0.1);
            if (jbb == jdim) jbb = jbb - 1;
            if (kbb == kdim) kbb = kbb - 1;
            if (ibb == idim) ibb = ibb - 1;
            ibb = std::max(ibb, 1);
            ssw = vor(jbb,kbb,ibb);
            wnu = fnu(jbb,kbb,ibb);
            ra  = std::sqrt(re*(ssw/wnu + 1.e-10*(double)xmach*(double)xmach));
            ypls = ra * ccabs_bb(smin(j,k,i));
            ypls = ccmaxcr_bb(ypls, .0001);
        }
        exp1 = std::exp(-(ypls/aplus1));
        exp2 = std::exp(-(ypls/aplus2));
        damp1(j,k,i) = (1.-exp1)*(1.-exp2);
        ddy = exp1/aplus1*(1.-exp2) + exp2/aplus2*(1.-exp1);
        dd  = damp1(j,k,i);
        sdd = std::sqrt(dd);
        damp2(j,k,i) = ratt + (1.-ratt)*(1./(akarman*ypls)+dd)*(sdd+ypls*ddy/sdd);
    }



    // Lines 272-310: tursav2 setup and turre initialization
    if ((float)dt > 0. && icyc == 1) {
        if (std::abs(ita) == 2) {
            if ((float)tursav2(1,1,1,1) == 0.) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    tursav2(j,k,i,3) = 0.;
                }
            } else {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    tursav2(j,k,i,3) = tursav(j,k,i,1) - tursav2(j,k,i,1);
                }
            }
        }
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            tursav2(j,k,i,1) = tursav(j,k,i,1);
        }
    }
    // Get TURRE values
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++)
    for (j = 1; j <= jdim-1; j++) {
        turre(j,k,i) = tursav(j,k,i,1);
    }

    // Main iteration loop: do 500 not=1,nnit
    for (not_ = 1; not_ <= nnit; not_++) {

        // (1) k=0 boundary
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            turre(j,0,i) = tk0(j,i,1,1);
        }
        // (2) k=kdim boundary
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            turre(j,kdim,i) = tk0(j,i,1,3);
        }
        // (3) j=0 boundary
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            turre(0,k,i) = tj0(k,i,1,1);
        }
        // (4) j=jdim boundary
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            turre(jdim,k,i) = tj0(k,i,1,3);
        }
        if (i2d != 1 && iaxi2planeturb != 1) {
            // (5) i=0 boundary
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,k,0) = ti0(j,k,1,1);
            }
            // (6) i=idim boundary
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,k,idim) = ti0(j,k,1,3);
            }
        }
        if (iturbord != 1) {
            // (1) k=0 boundary 2nd ghost
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,-1,i) = tk0(j,i,1,2);
            }
            // (2) k=kdim boundary 2nd ghost
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,kdim+1,i) = tk0(j,i,1,4);
            }
            // (3) j=0 boundary 2nd ghost
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                turre(-1,k,i) = tj0(k,i,1,2);
            }
            // (4) j=jdim boundary 2nd ghost
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                turre(jdim+1,k,i) = tj0(k,i,1,4);
            }
            if (i2d != 1 && iaxi2planeturb != 1) {
                // (5) i=0 boundary 2nd ghost
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    turre(j,k,-1) = ti0(j,k,1,2);
                }
                // (6) i=idim boundary 2nd ghost
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    turre(j,k,idim+1) = ti0(j,k,1,4);
                }
            }
        }



        // F_eta_eta viscous terms - Interior points
        for (k = 2; k <= kdim-2; k++) {
            kl = k-1;
            ku = k+1;
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                volku = vol(j,ku,i);
                xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                volkl = vol(j,kl,i);
                xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                cdp  = ttp*cnud;
                cdm  = ttm*cnud;
                trep = .5*(turre(j,k+1,i)+turre(j,k,i));
                trem = .5*(turre(j,k-1,i)+turre(j,k,i));
                cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/(sige*re);
                cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/(sige*re);
                byy  = ccmincr_bb(-cdm+cam, 0.);
                cyy  = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                dyy  = ccmincr_bb(-cdp+cap, 0.);
                vist3d(j,k,i) = -byy*turre(j,k-1,i) - cyy*turre(j,k,i) - dyy*turre(j,k+1,i);
            }
        }

        // K0 boundary points
        {
            k  = 1;
            kl = 1;
            ku = std::min(2, kdim-1);
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                volku = vol(j,ku,i);
                xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                volkl = volk0(j,i,1);
                xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                cdp  = ttp*cnud;
                cdm  = ttm*cnud;
                trep = .5*(turre(j,k+1,i)+turre(j,k,i));
                trem = .5*(turre(j,k-1,i)+turre(j,k,i));
                cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/(sige*re);
                cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/(sige*re);
                byy  = ccmincr_bb(-cdm+cam, 0.);
                cyy  = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                dyy  = ccmincr_bb(-cdp+cap, 0.);
                vist3d(j,k,i) = -byy*turre(j,k-1,i) - cyy*turre(j,k,i) - dyy*turre(j,k+1,i);
            }
        }

        // KDIM boundary points
        {
            k  = kdim-1;
            kl = kdim-2;
            ku = kdim-1;
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                volku = volk0(j,i,3);
                xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                volkl = vol(j,kl,i);
                xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                cdp  = ttp*cnud;
                cdm  = ttm*cnud;
                trep = .5*(turre(j,k+1,i)+turre(j,k,i));
                trem = .5*(turre(j,k-1,i)+turre(j,k,i));
                cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/(sige*re);
                cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/(sige*re);
                byy  = ccmincr_bb(-cdm+cam, 0.);
                cyy  = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                dyy  = ccmincr_bb(-cdp+cap, 0.);
                vist3d(j,k,i) = -byy*turre(j,k-1,i) - cyy*turre(j,k,i) - dyy*turre(j,k+1,i);
            }
        }



        // Advective terms in eta
        if (iturbord == 1) {
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                sgnu = ccsignrc_bb(1., uu);
                app  = 0.5*(1.+sgnu);
                apm  = 0.5*(1.-sgnu);
                vist3d(j,k,i) = vist3d(j,k,i) - uu*(app*(turre(j,k,i)-turre(j,k-1,i))
                                                    + apm*(turre(j,k+1,i)-turre(j,k,i)));
            }
        } else {
            // 2nd order upwind; LHS remains 1st order everywhere
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                sgnu = ccsignrc_bb(1., uu);
                app  = 0.5*(1.+sgnu);
                apm  = 0.5*(1.-sgnu);
                vist3d(j,k,i) = vist3d(j,k,i)
                    - 0.5*uu*app*turre(j,k-2,i)
                    + 2.*uu*app*turre(j,k-1,i)
                    - 1.5*uu*app*turre(j,k,i)
                    + 1.5*uu*apm*turre(j,k,i)
                    - 2.*uu*apm*turre(j,k+1,i)
                    + 0.5*uu*apm*turre(j,k+2,i);
            }
        }

        // F_xi_xi viscous terms - Interior points
        for (j = 2; j <= jdim-2; j++) {
            jl = j-1;
            ju = j+1;
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                volju = vol(ju,k,i);
                xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                voljl = vol(jl,k,i);
                xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                cdp  = ttp*cnud;
                cdm  = ttm*cnud;
                trep = .5*(turre(j+1,k,i)+turre(j,k,i));
                trem = .5*(turre(j-1,k,i)+turre(j,k,i));
                cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/(sige*re);
                cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/(sige*re);
                bxx  = ccmincr_bb(-cdm+cam, 0.);
                cxx  = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                dxx  = ccmincr_bb(-cdp+cap, 0.);
                vist3d(j,k,i) = vist3d(j,k,i) - bxx*turre(j-1,k,i) - cxx*turre(j,k,i) - dxx*turre(j+1,k,i);
            }
        }



        // J0 boundary points
        {
            j  = 1;
            jl = 1;
            ju = std::min(2, jdim-1);
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                volju = vol(ju,k,i);
                xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                voljl = volj0(k,i,1);
                xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                cdp  = ttp*cnud;
                cdm  = ttm*cnud;
                trep = .5*(turre(j+1,k,i)+turre(j,k,i));
                trem = .5*(turre(j-1,k,i)+turre(j,k,i));
                cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/(sige*re);
                cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/(sige*re);
                bxx  = ccmincr_bb(-cdm+cam, 0.);
                cxx  = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                dxx  = ccmincr_bb(-cdp+cap, 0.);
                vist3d(j,k,i) = vist3d(j,k,i) - bxx*turre(j-1,k,i) - cxx*turre(j,k,i) - dxx*turre(j+1,k,i);
            }
        }

        // JDIM boundary points
        {
            j  = jdim-1;
            jl = jdim-2;
            ju = jdim-1;
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                volju = volj0(k,i,3);
                xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                voljl = vol(jl,k,i);
                xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                cdp  = ttp*cnud;
                cdm  = ttm*cnud;
                trep = .5*(turre(j+1,k,i)+turre(j,k,i));
                trem = .5*(turre(j-1,k,i)+turre(j,k,i));
                cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/(sige*re);
                cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/(sige*re);
                bxx  = ccmincr_bb(-cdm+cam, 0.);
                cxx  = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                dxx  = ccmincr_bb(-cdp+cap, 0.);
                vist3d(j,k,i) = vist3d(j,k,i) - bxx*turre(j-1,k,i) - cxx*turre(j,k,i) - dxx*turre(j+1,k,i);
            }
        }

        // Advective terms in xi
        if (iturbord == 1) {
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                sgnu = ccsignrc_bb(1., uu);
                app  = 0.5*(1.+sgnu);
                apm  = 0.5*(1.-sgnu);
                vist3d(j,k,i) = vist3d(j,k,i) - uu*(app*(turre(j,k,i)-turre(j-1,k,i))
                                                    + apm*(turre(j+1,k,i)-turre(j,k,i)));
            }
        } else {
            // 2nd order upwind
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                sgnu = ccsignrc_bb(1., uu);
                app  = 0.5*(1.+sgnu);
                apm  = 0.5*(1.-sgnu);
                vist3d(j,k,i) = vist3d(j,k,i)
                    - 0.5*uu*app*turre(j-2,k,i)
                    + 2.*uu*app*turre(j-1,k,i)
                    - 1.5*uu*app*turre(j,k,i)
                    + 1.5*uu*apm*turre(j,k,i)
                    - 2.*uu*apm*turre(j+1,k,i)
                    + 0.5*uu*apm*turre(j+2,k,i);
            }
        }



        // F_zeta_zeta viscous terms
        if (i2d != 1 && iaxi2planeturb != 1) {
            // Interior points
            for (i = 2; i <= idim-2; i++) {
                il = i-1;
                iu = i+1;
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    voliu = vol(j,k,iu);
                    xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    volil = vol(j,k,il);
                    xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                    cdp  = ttp*cnud;
                    cdm  = ttm*cnud;
                    trep = .5*(turre(j,k,i+1)+turre(j,k,i));
                    trem = .5*(turre(j,k,i-1)+turre(j,k,i));
                    cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/(sige*re);
                    cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/(sige*re);
                    bzz  = ccmincr_bb(-cdm+cam, 0.);
                    czz  = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                    dzz  = ccmincr_bb(-cdp+cap, 0.);
                    vist3d(j,k,i) = vist3d(j,k,i) - bzz*turre(j,k,i-1) - czz*turre(j,k,i) - dzz*turre(j,k,i+1);
                }
            }

            // I0 boundary points
            {
                i  = 1;
                il = 1;
                iu = std::min(2, idim-1);
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    voliu = vol(j,k,iu);
                    xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    volil = voli0(j,k,1);
                    xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                    cdp  = ttp*cnud;
                    cdm  = ttm*cnud;
                    trep = .5*(turre(j,k,i+1)+turre(j,k,i));
                    trem = .5*(turre(j,k,i-1)+turre(j,k,i));
                    cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/(sige*re);
                    cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/(sige*re);
                    bzz  = ccmincr_bb(-cdm+cam, 0.);
                    czz  = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                    dzz  = ccmincr_bb(-cdp+cap, 0.);
                    vist3d(j,k,i) = vist3d(j,k,i) - bzz*turre(j,k,i-1) - czz*turre(j,k,i) - dzz*turre(j,k,i+1);
                }
            }

            // IDIM boundary points
            {
                i  = idim-1;
                il = idim-2;
                iu = idim-1;
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    voliu = voli0(j,k,3);
                    xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    volil = vol(j,k,il);
                    xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                    cdp  = ttp*cnud;
                    cdm  = ttm*cnud;
                    trep = .5*(turre(j,k,i+1)+turre(j,k,i));
                    trem = .5*(turre(j,k,i-1)+turre(j,k,i));
                    cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/(sige*re);
                    cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/(sige*re);
                    bzz  = ccmincr_bb(-cdm+cam, 0.);
                    czz  = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                    dzz  = ccmincr_bb(-cdp+cap, 0.);
                    vist3d(j,k,i) = vist3d(j,k,i) - bzz*turre(j,k,i-1) - czz*turre(j,k,i) - dzz*turre(j,k,i+1);
                }
            }

            // Advective terms in zeta
            if (iturbord == 1) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                    sgnu = ccsignrc_bb(1., uu);
                    app  = 0.5*(1.+sgnu);
                    apm  = 0.5*(1.-sgnu);
                    vist3d(j,k,i) = vist3d(j,k,i) - uu*(app*(turre(j,k,i)-turre(j,k,i-1))
                                                        + apm*(turre(j,k,i+1)-turre(j,k,i)));
                }
            } else {
                // 2nd order upwind
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                    sgnu = ccsignrc_bb(1., uu);
                    app  = 0.5*(1.+sgnu);
                    apm  = 0.5*(1.-sgnu);
                    vist3d(j,k,i) = vist3d(j,k,i)
                        - 0.5*uu*app*turre(j,k,i-2)
                        + 2.*uu*app*turre(j,k,i-1)
                        - 1.5*uu*app*turre(j,k,i)
                        + 1.5*uu*apm*turre(j,k,i)
                        - 2.*uu*apm*turre(j,k,i+1)
                        + 0.5*uu*apm*turre(j,k,i+2);
                }
            }
        } // end if (i2d != 1 && iaxi2planeturb != 1) for F_zeta_zeta



        // Add source term to RHS
        if (isklton > 0) {
            if (ilamlo == 0 || jlamlo == 0 || klamlo == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   block %4d in B-B turb model has no laminar regions", nbl);
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   block %4d in B-B turb model - laminar region is:", nbl);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   i=%5d to%5d, j=%5d to%5d, k=%5d to%5d",
                    ilamlo, ilamhi, jlamlo, jlamhi, klamlo, klamhi);
                if (i_lam_forcezero == 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "    ...forcing vist3d=0");
                }
            }
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   NOTE:  This particular model <<transitions>> on its own, but there is");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   no guarantee that it will transition at all.  Check vist3d levels if unsure.");
        }
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++)
        for (k = 1; k <= kdim-1; k++) {
            if ((i >= ilamlo && i < ilamhi &&
                 j >= jlamlo && j < jlamhi &&
                 k >= klamlo && k < klamhi) ||
                (float)smin(j,k,i) < 0.) {
                cutoff = 0.;
            } else {
                cutoff = 1.;
            }
            vist3d(j,k,i) = vist3d(j,k,i) + (c2e*damp2(j,k,i)-c1e)*
                std::sqrt(cmu*damp1(j,k,i))*vor(j,k,i)*cutoff*turre(j,k,i);
        }

        // Implicit F_eta_eta viscous terms. Do over all i's
        for (i = 1; i <= idim-1; i++) {
            // Interior points
            for (k = 2; k <= kdim-2; k++) {
                kl = k-1;
                ku = k+1;
                for (j = 1; j <= jdim-1; j++) {
                    volku = vol(j,ku,i);
                    xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    volkl = vol(j,kl,i);
                    xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                    ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                    cdp  = ttp*cnud;
                    cdm  = ttm*cnud;
                    trep = .5*(turre(j,k+1,i)+turre(j,k,i));
                    trem = .5*(turre(j,k-1,i)+turre(j,k,i));
                    cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/(sige*re);
                    cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/(sige*re);
                    by(j,k) = ccmincr_bb(-cdm+cam, 0.);
                    cy(j,k) = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                    dy(j,k) = ccmincr_bb(-cdp+cap, 0.);
                }
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                    uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                    sgnu = ccsignrc_bb(1., uu);
                    app  = 0.5*(1.+sgnu);
                    apm  = 0.5*(1.-sgnu);
                    by(j,k) = by(j,k) - uu*app;
                    cy(j,k) = cy(j,k) + uu*(app-apm);
                    dy(j,k) = dy(j,k) + uu*apm;
                }
                for (j = 1; j <= jdim-1; j++) {
                    fact    = timestp(j,k,i);
                    by(j,k) = by(j,k)*fact;
                    cy(j,k) = cy(j,k)*fact + 1.0*(1.+phi);
                    dy(j,k) = dy(j,k)*fact;
                    fy(j,k) = vist3d(j,k,i)*fact;
                }
                if ((float)dt > 0.) {
                    for (j = 1; j <= jdim-1; j++) {
                        fy(j,k) = fy(j,k) + (1.+phi)*(tursav2(j,k,i,1)-turre(j,k,i))
                                           + phi*tursav2(j,k,i,3);
                    }
                }
            } // end interior k loop



            // K0 boundary points (implicit eta)
            {
                k  = 1;
                kl = 1;
                ku = std::min(2, kdim-1);
                for (j = 1; j <= jdim-1; j++) {
                    volku = vol(j,ku,i);
                    xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    volkl = volk0(j,i,1);
                    xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                    ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                    cdp  = ttp*cnud;
                    cdm  = ttm*cnud;
                    trep = .5*(turre(j,k+1,i)+turre(j,k,i));
                    trem = .5*(turre(j,k-1,i)+turre(j,k,i));
                    cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/(sige*re);
                    cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/(sige*re);
                    by(j,k) = ccmincr_bb(-cdm+cam, 0.);
                    cy(j,k) = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                    dy(j,k) = ccmincr_bb(-cdp+cap, 0.);
                }
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                    uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                    sgnu = ccsignrc_bb(1., uu);
                    app  = 0.5*(1.+sgnu);
                    apm  = 0.5*(1.-sgnu);
                    by(j,k) = by(j,k) - uu*app;
                    cy(j,k) = cy(j,k) + uu*(app-apm);
                    dy(j,k) = dy(j,k) + uu*apm;
                }
                for (j = 1; j <= jdim-1; j++) {
                    fact    = timestp(j,k,i);
                    by(j,k) = by(j,k)*fact;
                    cy(j,k) = cy(j,k)*fact + 1.0*(1.+phi);
                    dy(j,k) = dy(j,k)*fact;
                    fy(j,k) = vist3d(j,k,i)*fact;
                }
                if ((float)dt > 0.) {
                    for (j = 1; j <= jdim-1; j++) {
                        fy(j,k) = fy(j,k) + (1.+phi)*(tursav2(j,k,i,1)-turre(j,k,i))
                                           + phi*tursav2(j,k,i,3);
                    }
                }
            }

            // KDIM boundary points (implicit eta)
            {
                k  = kdim-1;
                kl = kdim-2;
                ku = kdim-1;
                for (j = 1; j <= jdim-1; j++) {
                    volku = vol(j,i,3);  // Fortran: vol(j,i,3) - faithful translation
                    xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    volkl = vol(j,kl,i);
                    xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                    ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                    cdp  = ttp*cnud;
                    cdm  = ttm*cnud;
                    trep = .5*(turre(j,k+1,i)+turre(j,k,i));
                    trem = .5*(turre(j,k-1,i)+turre(j,k,i));
                    cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/(sige*re);
                    cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/(sige*re);
                    by(j,k) = ccmincr_bb(-cdm+cam, 0.);
                    cy(j,k) = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                    dy(j,k) = ccmincr_bb(-cdp+cap, 0.);
                }
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                    uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                    sgnu = ccsignrc_bb(1., uu);
                    app  = 0.5*(1.+sgnu);
                    apm  = 0.5*(1.-sgnu);
                    by(j,k) = by(j,k) - uu*app;
                    cy(j,k) = cy(j,k) + uu*(app-apm);
                    dy(j,k) = dy(j,k) + uu*apm;
                }
                for (j = 1; j <= jdim-1; j++) {
                    fact    = timestp(j,k,i);
                    by(j,k) = by(j,k)*fact;
                    cy(j,k) = cy(j,k)*fact + 1.0*(1.+phi);
                    dy(j,k) = dy(j,k)*fact;
                    fy(j,k) = vist3d(j,k,i)*fact;
                }
                if ((float)dt > 0.) {
                    for (j = 1; j <= jdim-1; j++) {
                        fy(j,k) = fy(j,k) + (1.+phi)*(tursav2(j,k,i,1)-turre(j,k,i))
                                           + phi*tursav2(j,k,i,3);
                    }
                }
            }

            if (iover == 1) {
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    fy(j,k) = fy(j,k)*blank(j,k,i);
                    by(j,k) = by(j,k)*blank(j,k,i);
                    dy(j,k) = dy(j,k)*blank(j,k,i);
                    cy(j,k) = cy(j,k)*blank(j,k,i) + (1.-blank(j,k,i));
                }
            }
            { int n1=jdim-1, n2=kdim-1, i1=1, i2=jdim-1, j1=1, j2=kdim-1;
              triv_ns::triv(n1, n2, i1, i2, j1, j2, worky, by, cy, dy, fy); }
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                vist3d(j,k,i) = fy(j,k);
            }
        } // end do i (implicit eta)



        // Implicit F_xi_xi viscous terms. Do over all i's
        for (i = 1; i <= idim-1; i++) {
            // Interior points
            for (j = 2; j <= jdim-2; j++) {
                jl = j-1;
                ju = j+1;
                for (k = 1; k <= kdim-1; k++) {
                    volju = vol(ju,k,i);
                    xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    voljl = vol(jl,k,i);
                    xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                    ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                    cdp  = ttp*cnud;
                    cdm  = ttm*cnud;
                    trep = .5*(turre(j+1,k,i)+turre(j,k,i));
                    trem = .5*(turre(j-1,k,i)+turre(j,k,i));
                    cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/(sige*re);
                    cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/(sige*re);
                    bx(k,j) = ccmincr_bb(-cdm+cam, 0.);
                    cx(k,j) = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                    dx(k,j) = ccmincr_bb(-cdp+cap, 0.);
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                    uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                    sgnu = ccsignrc_bb(1., uu);
                    app  = 0.5*(1.+sgnu);
                    apm  = 0.5*(1.-sgnu);
                    bx(k,j) = bx(k,j) - uu*app;
                    cx(k,j) = cx(k,j) + uu*(app-apm);
                    dx(k,j) = dx(k,j) + uu*apm;
                }
                for (k = 1; k <= kdim-1; k++) {
                    fact    = timestp(j,k,i);
                    bx(k,j) = bx(k,j)*fact;
                    cx(k,j) = cx(k,j)*fact + 1.0*(1.+phi);
                    dx(k,j) = dx(k,j)*fact;
                    fx(k,j) = vist3d(j,k,i)*(1.+phi);
                }
            } // end interior j loop

            // J0 boundary points (implicit xi)
            {
                j  = 1;
                jl = 1;
                ju = std::min(2, jdim-1);
                for (k = 1; k <= kdim-1; k++) {
                    volju = vol(ju,k,i);
                    xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    voljl = volj0(k,i,1);
                    xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                    ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                    cdp  = ttp*cnud;
                    cdm  = ttm*cnud;
                    trep = .5*(turre(j+1,k,i)+turre(j,k,i));
                    trem = .5*(turre(j-1,k,i)+turre(j,k,i));
                    cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/(sige*re);
                    cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/(sige*re);
                    bx(k,j) = ccmincr_bb(-cdm+cam, 0.);
                    cx(k,j) = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                    dx(k,j) = ccmincr_bb(-cdp+cap, 0.);
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                    uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                    sgnu = ccsignrc_bb(1., uu);
                    app  = 0.5*(1.+sgnu);
                    apm  = 0.5*(1.-sgnu);
                    bx(k,j) = bx(k,j) - uu*app;
                    cx(k,j) = cx(k,j) + uu*(app-apm);
                    dx(k,j) = dx(k,j) + uu*apm;
                }
                for (k = 1; k <= kdim-1; k++) {
                    fact    = timestp(j,k,i);
                    bx(k,j) = bx(k,j)*fact;
                    cx(k,j) = cx(k,j)*fact + 1.0*(1.+phi);
                    dx(k,j) = dx(k,j)*fact;
                    fx(k,j) = vist3d(j,k,i)*(1.+phi);
                }
            }



            // JDIM boundary points (implicit xi)
            {
                j  = jdim-1;
                jl = jdim-2;
                ju = jdim-1;
                for (k = 1; k <= kdim-1; k++) {
                    volju = volj0(k,i,3);
                    xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    voljl = vol(jl,k,i);
                    xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                    ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                    cdp  = ttp*cnud;
                    cdm  = ttm*cnud;
                    trep = .5*(turre(j+1,k,i)+turre(j,k,i));
                    trem = .5*(turre(j-1,k,i)+turre(j,k,i));
                    cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/(sige*re);
                    cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/(sige*re);
                    bx(k,j) = ccmincr_bb(-cdm+cam, 0.);
                    cx(k,j) = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                    dx(k,j) = ccmincr_bb(-cdp+cap, 0.);
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                    uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                    sgnu = ccsignrc_bb(1., uu);
                    app  = 0.5*(1.+sgnu);
                    apm  = 0.5*(1.-sgnu);
                    bx(k,j) = bx(k,j) - uu*app;
                    cx(k,j) = cx(k,j) + uu*(app-apm);
                    dx(k,j) = dx(k,j) + uu*apm;
                }
                for (k = 1; k <= kdim-1; k++) {
                    fact    = timestp(j,k,i);
                    bx(k,j) = bx(k,j)*fact;
                    cx(k,j) = cx(k,j)*fact + 1.0*(1.+phi);
                    dx(k,j) = dx(k,j)*fact;
                    fx(k,j) = vist3d(j,k,i)*(1.+phi);
                }
            }

            if (iover == 1) {
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    fx(k,j) = fx(k,j)*blank(j,k,i);
                    bx(k,j) = bx(k,j)*blank(j,k,i);
                    dx(k,j) = dx(k,j)*blank(j,k,i);
                    cx(k,j) = cx(k,j)*blank(j,k,i) + (1.-blank(j,k,i));
                }
            }
            { int n1=kdim-1, n2=jdim-1, i1=1, i2=kdim-1, j1=1, j2=jdim-1;
              triv_ns::triv(n1, n2, i1, i2, j1, j2, workx, bx, cx, dx, fx); }
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                vist3d(j,k,i) = fx(k,j);
            }
        } // end do i (implicit xi)



        // Implicit F_zeta_zeta viscous terms. Do over all j's
        if (i2d != 1 && iaxi2planeturb != 1) {
            for (j = 1; j <= jdim-1; j++) {
                // Interior points
                for (i = 2; i <= idim-2; i++) {
                    il = i-1;
                    iu = i+1;
                    for (k = 1; k <= kdim-1; k++) {
                        voliu = vol(j,k,iu);
                        xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        volil = vol(j,k,il);
                        xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        ttpo = xp*xa+yp*ya+zp*za;
                        ttmo = xm*xa+ym*ya+zm*za;
                        ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                        ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                        ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                        ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                        cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                        cdp  = ttp*cnud;
                        cdm  = ttm*cnud;
                        trep = .5*(turre(j,k,i+1)+turre(j,k,i));
                        trem = .5*(turre(j,k,i-1)+turre(j,k,i));
                        cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/(sige*re);
                        cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/(sige*re);
                        bz(k,i) = ccmincr_bb(-cdm+cam, 0.);
                        cz(k,i) = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                        dz(k,i) = ccmincr_bb(-cdp+cap, 0.);
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                        uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                        sgnu = ccsignrc_bb(1., uu);
                        app  = 0.5*(1.+sgnu);
                        apm  = 0.5*(1.-sgnu);
                        bz(k,i) = bz(k,i) - uu*app;
                        cz(k,i) = cz(k,i) + uu*(app-apm);
                        dz(k,i) = dz(k,i) + uu*apm;
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        fact    = timestp(j,k,i);
                        bz(k,i) = bz(k,i)*fact;
                        cz(k,i) = cz(k,i)*fact + 1.0*(1.+phi);
                        dz(k,i) = dz(k,i)*fact;
                        fz(k,i) = vist3d(j,k,i)*(1.+phi);
                    }
                } // end interior i loop

                // I0 boundary points (implicit zeta)
                {
                    i  = 1;
                    il = idim-2;
                    iu = idim-1;
                    for (k = 1; k <= kdim-1; k++) {
                        voliu = vol(j,k,iu);
                        xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        volil = voli0(j,k,1);
                        xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        ttpo = xp*xa+yp*ya+zp*za;
                        ttmo = xm*xa+ym*ya+zm*za;
                        ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                        ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                        ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                        ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                        cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                        cdp  = ttp*cnud;
                        cdm  = ttm*cnud;
                        trep = .5*(turre(j,k,i+1)+turre(j,k,i));
                        trem = .5*(turre(j,k,i-1)+turre(j,k,i));
                        cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/(sige*re);
                        cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/(sige*re);
                        bz(k,i) = ccmincr_bb(-cdm+cam, 0.);
                        cz(k,i) = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                        dz(k,i) = ccmincr_bb(-cdp+cap, 0.);
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                        uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                        sgnu = ccsignrc_bb(1., uu);
                        app  = 0.5*(1.+sgnu);
                        apm  = 0.5*(1.-sgnu);
                        bz(k,i) = bz(k,i) - uu*app;
                        cz(k,i) = cz(k,i) + uu*(app-apm);
                        dz(k,i) = dz(k,i) + uu*apm;
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        fact    = timestp(j,k,i);
                        bz(k,i) = bz(k,i)*fact;
                        cz(k,i) = cz(k,i)*fact + 1.0*(1.+phi);
                        dz(k,i) = dz(k,i)*fact;
                        fz(k,i) = vist3d(j,k,i)*(1.+phi);
                    }
                }



                // IDIM boundary points (implicit zeta)
                {
                    i  = idim-1;
                    il = idim-2;
                    iu = idim-1;
                    for (k = 1; k <= kdim-1; k++) {
                        voliu = voli0(j,k,3);
                        xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        volil = vol(j,k,il);
                        xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        ttpo = xp*xa+yp*ya+zp*za;
                        ttmo = xm*xa+ym*ya+zm*za;
                        ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                        ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                        ttp  = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                        ttm  = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                        cnud = (fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*(2./sige))/re;
                        cdp  = ttp*cnud;
                        cdm  = ttm*cnud;
                        trep = .5*(turre(j,k,i+1)+turre(j,k,i));
                        trem = .5*(turre(j,k,i-1)+turre(j,k,i));
                        cap  = cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/(sige*re);
                        cam  = cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/(sige*re);
                        bz(k,i) = ccmincr_bb(-cdm+cam, 0.);
                        cz(k,i) = ccmaxcr_bb(cdp-cap, 0.) + ccmaxcr_bb(cdm-cam, 0.);
                        dz(k,i) = ccmincr_bb(-cdp+cap, 0.);
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                        uu   = b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc;
                        sgnu = ccsignrc_bb(1., uu);
                        app  = 0.5*(1.+sgnu);
                        apm  = 0.5*(1.-sgnu);
                        bz(k,i) = bz(k,i) - uu*app;
                        cz(k,i) = cz(k,i) + uu*(app-apm);
                        dz(k,i) = dz(k,i) + uu*apm;
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        fact    = timestp(j,k,i);
                        bz(k,i) = bz(k,i)*fact;
                        cz(k,i) = cz(k,i)*fact + 1.0*(1.+phi);
                        dz(k,i) = dz(k,i)*fact;
                        fz(k,i) = vist3d(j,k,i)*(1.+phi);
                    }
                }

                if (iover == 1) {
                    for (i = 1; i <= idim-1; i++)
                    for (k = 1; k <= kdim-1; k++) {
                        fz(k,i) = fz(k,i)*blank(j,k,i);
                        bz(k,i) = bz(k,i)*blank(j,k,i);
                        dz(k,i) = dz(k,i)*blank(j,k,i);
                        cz(k,i) = cz(k,i)*blank(j,k,i) + (1.-blank(j,k,i));
                    }
                }
                { int n1=kdim-1, n2=idim-1, i1=1, i2=kdim-1, j1=1, j2=idim-1;
                  triv_ns::triv(n1, n2, i1, i2, j1, j2, workz, bz, cz, dz, fz); }
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++) {
                    vist3d(j,k,i) = fz(k,i);
                }
            } // end do j (implicit zeta)
        } // end if (i2d != 1 && iaxi2planeturb != 1) for implicit zeta

        // Update TURRE
        sumn = 0.;
        negn = 0;
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            if ((float)(turre(j,k,i)+vist3d(j,k,i)) < 1.e-12f) {
                negn++;
                turre(j,k,i) = 1.e-12;
                vist3d(j,k,i) = 0.;
            } else {
                turre(j,k,i) = turre(j,k,i) + vist3d(j,k,i);
            }
            sumn = sumn + vist3d(j,k,i)*vist3d(j,k,i);
        }

        sumn = std::sqrt(sumn) / (double)((kdim-1)*(jdim-1)*(idim-1));

        if (iwrite == 1) {
            // write(15,...) icyc,not_,log10(sumn),negn,nbl  -- commented out in Fortran
        }

    } // end do 500 not=1,nnit

    sumn1 = sumn;
    sumn2 = 1.;
    negn1 = negn;
    negn2 = 0;

    // Update VIST3D
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++)
    for (j = 1; j <= jdim-1; j++) {
        tursav(j,k,i,1) = turre(j,k,i);
        vist3d(j,k,i)   = cmu*damp1(j,k,i)*turre(j,k,i)*q(j,k,i,1);
    }

    if (i_lam_forcezero == 1) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            if ((i >= ilamlo && i < ilamhi &&
                 j >= jlamlo && j < jlamhi &&
                 k >= klamlo && k < klamhi) ||
                (float)smin(j,k,i) < 0.) {
                vist3d(j,k,i) = 0.;
            }
        }
    }

} // end barth3d

} // namespace barth3d_ns
