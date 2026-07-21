// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "threeeqn.h"
#include "ccomplex.h"
#include "triv.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstring>

namespace threeeqn_ns {

void threeeqn(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, FortranArray3DRef<double> dtj, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vor, FortranArray3DRef<double> smin, FortranArray4DRef<double> zksav, FortranArray4DRef<double> turre, FortranArray3DRef<double> damp1, FortranArray3DRef<double> blend, FortranArray3DRef<double> timestp, FortranArray3DRef<double> fnu, FortranArray2DRef<double> bx, FortranArray2DRef<double> bx2, FortranArray2DRef<double> cx, FortranArray2DRef<double> cx2, FortranArray2DRef<double> dx, FortranArray2DRef<double> dx2, FortranArray2DRef<double> fx, FortranArray2DRef<double> fx2, FortranArray2DRef<double> workx, FortranArray2DRef<double> by, FortranArray2DRef<double> by2, FortranArray2DRef<double> cy, FortranArray2DRef<double> cy2, FortranArray2DRef<double> dy, FortranArray2DRef<double> dy2, FortranArray2DRef<double> fy, FortranArray2DRef<double> fy2, FortranArray2DRef<double> worky, FortranArray2DRef<double> bz, FortranArray2DRef<double> bz2, FortranArray2DRef<double> cz, FortranArray2DRef<double> cz2, FortranArray2DRef<double> dz, FortranArray2DRef<double> dz2, FortranArray2DRef<double> fz, FortranArray2DRef<double> fz2, FortranArray2DRef<double> workz, int& ntime, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, int& nbl, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, FortranArray3DRef<double> blank, int& iover, double& sumn1, double& sumn2, double& sumn3, int& negn1, int& negn2, int& negn3, FortranArray4DRef<double> ux, FortranArray4DRef<double> rhside, FortranArray4DRef<double> zksav2, FortranArray3DRef<double> v3dtmp, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& maxbl, int& maxseg, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iex, int& iex2, int& iex3, FortranArray2DRef<double> bx3, FortranArray2DRef<double> cx3, FortranArray2DRef<double> dx3, FortranArray2DRef<double> fx3, FortranArray2DRef<double> by3, FortranArray2DRef<double> cy3, FortranArray2DRef<double> dy3, FortranArray2DRef<double> fy3, FortranArray2DRef<double> bz3, FortranArray2DRef<double> cz3, FortranArray2DRef<double> dz3, FortranArray2DRef<double> fz3, int& nummem)
{
    // COMMON block aliases
    float& gamma    = cmn_fluid.gamma;
    float& gm1      = cmn_fluid.gm1;
    float& pr       = cmn_fluid2.pr;
    float& prt      = cmn_fluid2.prt;
    float& cbar     = cmn_fluid2.cbar;
    float& dt       = cmn_info.dt;
    int32_t& icyc   = cmn_mgrd.icyc;
    int32_t& ncyc   = cmn_mgrd.ncyc;
    int32_t& level  = cmn_mgrd.level;
    int32_t& lglobal= cmn_mgrd.lglobal;
    int32_t* ncyc1  = cmn_info.ncyc1;  // 0-based
    float& reue     = cmn_reyue.reue;
    float& tinf     = cmn_reyue.tinf;
    int32_t* ivisc  = cmn_reyue.ivisc; // 0-based
    int32_t& isklton= cmn_sklton.isklton;
    int32_t& i2d    = cmn_twod.i2d;
    int32_t& iexp   = cmn_zero.iexp;
    int32_t* iwf    = cmn_wallfun.iwf;  // 0-based
    float* cflturb  = cmn_turbconv.cflturb; // 0-based
    float& edvislim = cmn_turbconv.edvislim;
    int32_t& iturbprod = cmn_turbconv.iturbprod;
    int32_t& nsubturb  = cmn_turbconv.nsubturb;
    int32_t& nfreeze   = cmn_turbconv.nfreeze;
    int32_t& itaturb   = cmn_turbconv.itaturb;
    float& tur1cut     = cmn_turbconv.tur1cut;
    float& tur2cut     = cmn_turbconv.tur2cut;
    int32_t& iturbord  = cmn_turbconv.iturbord;
    float& tur1cutlev  = cmn_turbconv.tur1cutlev;
    float& tur2cutlev  = cmn_turbconv.tur2cutlev;
    float& time        = cmn_unst.time;
    int32_t& ita       = cmn_unst.ita;
    int32_t& ivmx      = cmn_maxiv.ivmx;
    float* tur10       = cmn_ivals.tur10; // 0-based
    float& rho0        = cmn_ivals.rho0;
    int32_t& ilamlo    = cmn_lam.ilamlo;
    int32_t& ilamhi    = cmn_lam.ilamhi;
    int32_t& jlamlo    = cmn_lam.jlamlo;
    int32_t& jlamhi    = cmn_lam.jlamhi;
    int32_t& klamlo    = cmn_lam.klamlo;
    int32_t& klamhi    = cmn_lam.klamhi;
    int32_t& i_lam_forcezero = cmn_lam.i_lam_forcezero;
    int32_t& ikoprod   = cmn_konew.ikoprod;
    int32_t& isstdenom = cmn_konew.isstdenom;
    float& pklimterm   = cmn_konew.pklimterm;
    int32_t& keepambient = cmn_konew.keepambient;
    float& re_thetat0  = cmn_konew.re_thetat0;
    float& prod2d3dtrace = cmn_konew.prod2d3dtrace;
    int32_t& iaxi2planeturb = cmn_axisym.iaxi2planeturb;
    int32_t& istrongturbdis = cmn_axisym.istrongturbdis;
    float& xmach       = cmn_info.xmach;

    // Local variables
    int itrans_on;
    double tu_percent, re_thetat, re_thetac, flength, ce2;
    double factor1, factor2, factor3;
    int nsubit;
    double vk, a1, cmuc1, beta1, sigo1, sigk1, alp1;
    double sigg, sigkmu, cmuc2, beta2, sigo2, sigk2, alp2;
    int jd2;
    double re, c2b, c2bp;
    int iwrite;
    double phi, xminn;
    int nss;
    int i, j, k, il, iu, jl, ju, kl, ku;
    double xa, ya, za, xp, yp, zp, xm, ym, zm;
    double ttpo, ttmo, ttpn, ttmn, ttp, ttm;
    double anutp, anutm, fnup, fnum, cdp, cdm;
    double byy, cyy, dyy, bxx, cxx, dxx, bzz, czz, dzz;
    double dfacep, dfacem, sigkp, sigkm, sigop, sigom;
    double volku, volkl, volju, voljl, voliu, volil;
    double xc, yc, zc, tc, uu, sgnu, app, apm;
    double cyadd, cy2add, betax, cmuc, alp;
    double s11, s22, s33, s12, s13, s23, tracepart, s11t, s22t, s33t, xis;
    double pk, dk, dkfactor;
    double re_nu, re_t, fonset1, fonset2, fonset3, fonset, fturb;
    double ftrnlength, trnlength;
    double re_thetat_tilde, uuu, theta_bl, delta_bl, delta, re_omega;
    double fwake, f_thetat1, f_thetat2, f_thetat, f_reattach;
    double gamma_sep1, gamma_sep, gamma_eff, gamma_lim;
    double denom1, denom2, arg2, f2, denom, xxx;
    double sumno, sumnk, sumni;
    int negno, negnk, negni;
    double arg1, arga, argb, arg, temp, ca, tt;
    double cutoff;
    // iwriteaux section
    int iwriteaux, nnumb, jset;
    double qset, utau, ypl, uplus, zkplus, eplus, uvplus;
    int not_; // loop variable (not is reserved in C++)



    // itrans_on=1
    itrans_on = 1;
    // Set transition-related parameters
    tu_percent = 100.*std::sqrt(2.*tur10[1]/(3.*(double)xmach*(double)xmach));
    if ((float)re_thetat0 < 0.) {
        if (ivmx == 30) {
            re_thetat = 400.*std::pow(tu_percent, -5./8.);
        } else if (ivmx == 31) {
            re_thetat = 803.73*std::pow(tu_percent+0.6067, -1.027);
        } else {
            re_thetat = 0.0; // uninitialized in Fortran for other ivmx
        }
    } else {
        re_thetat = (double)re_thetat0;
    }
    re_thetac = 0.0; flength = 0.0; ce2 = 0.0;
    if (ivmx == 30) {
        re_thetac = 0.85*re_thetat;
    } else if (ivmx == 31) {
        re_thetac = re_thetat;
        { double tmp = tu_percent; flength = 163.*std::log((float)tmp)+3.625; }
        ce2 = 50.;
    }

    if (isklton > 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     Computing turbulent viscosity using 3-eqns, block=%5d", nbl);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     WARNING: ivisc=30 still under development... use at your own risk!");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     Freestream tur10,tur20,tur30 = %19.8e%19.8e%19.8e",
            (double)tur10[0], (double)tur10[1], (double)tur10[2]);
        if (iturbord == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     1st order advection on RHS");
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     2nd order advection on RHS");
        }
    }
    if (isklton > 0) {
        if (ivmx == 30) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-omega SST (Menter) + transition model (Langtry unpublished version)");
        } else if (ivmx == 31) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-omega SST (Menter) + transition model (JTurbomach 128:413 2006 and IJHFF 29:48 2008)");
        }
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     Tu (percent)=%10.4f", tu_percent);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     mut_inf/muref=%10.4f, Re_thetat=%10.4f",
            (double)rho0*(double)tur10[1]/(double)tur10[0], re_thetat);
        if (ikoprod == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     strain-based production term");
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     approx (vort) production term");
        }
        if (isstdenom == 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     vort in denom of mut term");
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     S in denom of mut term");
        }
        if (std::abs((double)prod2d3dtrace-0.5) < 0.01) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     Sij used in 2SijSij prod term forced to be traceless in 2-D sense");
        } else if (std::abs((double)prod2d3dtrace-0.33333333) < 0.01) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     Sij used in 2SijSij prod term forced to be traceless in 3-D sense");
        }
        if (keepambient == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     3-eqn ambient turbulence levels not allowed to decay");
        }
        if (iaxi2planeturb == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     3-eqn model ignoring i-dir");
        }
        if (istrongturbdis == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     strong conserv - diss terms");
        }
    }

    // Note: (10.**(-iexp) is machine zero)
    xminn = std::pow(10., (double)(-iexp+1));

    // Set number of subiterations
    nsubit = nsubturb;

    // Set CFL factors
    factor1 = 10.;
    factor2 = 10.;
    factor3 = 100.;

    if ((float)cflturb[0] != 0.) factor1 = (double)cflturb[0];
    if ((float)cflturb[1] != 0.) factor2 = (double)cflturb[1];
    if ((float)cflturb[2] != 0.) factor3 = (double)cflturb[2];
    // factor2 and 3 are set relative to factor1
    factor2 = factor2/factor1;
    factor3 = factor3/factor1;

    // Timestep for turb model
    if ((float)dt < 0) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            { double a1t=factor1*vol(j,k,i)/dtj(j,k,i), a2t=100.; timestp(j,k,i)=ccomplex_ns::ccmincr(a1t,a2t); }
        }
    } else {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            timestp(j,k,i) = (double)dt;
            factor2 = 1.;
            factor3 = 1.;
        }
    }

    // Set up constants
    vk  = .41;
    a1  = .31;
    cmuc1 = 0.09;
    beta1 = 0.0; sigo1 = 0.0; sigk1 = 0.0; alp1 = 0.0;
    if (ivmx == 30 || ivmx == 31) beta1 = 0.075;
    if (ivmx == 30 || ivmx == 31) sigo1 = 0.5;
    if (ivmx == 30 || ivmx == 31) sigk1 = 0.85;
    if (ivmx == 30 || ivmx == 31) alp1  = beta1/cmuc1 - sigo1*vk*vk/std::sqrt(cmuc1);
    sigg   = 1.0;
    sigkmu = 1.0;
    cmuc2 = 0.09;
    beta2 = 0.0828;
    sigo2 = 0.856;
    sigk2 = 1.00;
    alp2  = beta2/cmuc2 - sigo2*vk*vk/std::sqrt(cmuc2);

    jd2  = (jdim-1)/2;
    re   = (double)reue/(double)xmach;
    c2b  = (double)cbar/(double)tinf;
    c2bp = c2b + 1.0;

    iwrite = 0;

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
        sumn3 = 0.;
        negn1 = 0;
        negn2 = 0;
        negn3 = 0;
        return;
    }
    phi = 0.;
    if ((float)dt > 0.) {
        if (std::abs(ita) == 2) {
            phi = 0.5;
        } else {
            phi = 0.;
        }
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



    // Get laminar viscosity at cell centers
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++)
    for (j = 1; j <= jdim-1; j++) {
        tt = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
        fnu(j,k,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
    }
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++) {
        tt = (double)gamma*qj0(k,i,5,1)/qj0(k,i,1,1);
        fnu(0,k,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
        tt = (double)gamma*qj0(k,i,5,3)/qj0(k,i,1,3);
        fnu(jdim,k,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
    }
    for (i = 1; i <= idim-1; i++)
    for (j = 1; j <= jdim-1; j++) {
        tt = (double)gamma*qk0(j,i,5,1)/qk0(j,i,1,1);
        fnu(j,0,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
        tt = (double)gamma*qk0(j,i,5,3)/qk0(j,i,1,3);
        fnu(j,kdim,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
    }
    if (i2d != 1 && iaxi2planeturb != 1) {
        for (j = 1; j <= jdim-1; j++)
        for (k = 1; k <= kdim-1; k++) {
            tt = (double)gamma*qi0(j,k,5,1)/qi0(j,k,1,1);
            fnu(j,k,0) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
            tt = (double)gamma*qi0(j,k,5,3)/qi0(j,k,1,3);
            fnu(j,k,idim) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
        }
    }

    // Load appropriate turb viscosity at cell centers
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++)
    for (j = 1; j <= jdim-1; j++) {
        v3dtmp(j,k,i) = vist3d(j,k,i);
    }
    // Load appropriate vist3d value into ghost cells
    for (i = 1; i <= idim-1; i++)
    for (j = 1; j <= jdim-1; j++) {
        v3dtmp(j,0,i) = bck(j,i,1)*(iwf[2]*v3dtmp(j,1,i) +
            (1-iwf[2])*2.*vk0(j,i,1,1)-v3dtmp(j,1,i))+
            (1.-bck(j,i,1))*vk0(j,i,1,1);
        v3dtmp(j,kdim,i) = bck(j,i,2)*(iwf[2]*v3dtmp(j,kdim-1,i) +
            (1-iwf[2])*2.*vk0(j,i,1,3)-v3dtmp(j,kdim-1,i))+
            (1.-bck(j,i,2))*vk0(j,i,1,3);
    }
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++) {
        v3dtmp(0,k,i) = bcj(k,i,1)*(iwf[1]*v3dtmp(1,k,i) +
            (1-iwf[1])*2.*vj0(k,i,1,1)-v3dtmp(1,k,i))+
            (1.-bcj(k,i,1))*vj0(k,i,1,1);
        v3dtmp(jdim,k,i) = bcj(k,i,2)*(iwf[1]*v3dtmp(jdim-1,k,i) +
            (1-iwf[1])*2.*vj0(k,i,1,3)-v3dtmp(jdim-1,k,i))+
            (1.-bcj(k,i,2))*vj0(k,i,1,3);
    }
    if (i2d != 1 && iaxi2planeturb != 1) {
        for (j = 1; j <= jdim-1; j++)
        for (k = 1; k <= kdim-1; k++) {
            v3dtmp(j,k,0) = bci(j,k,1)*(iwf[0]*v3dtmp(j,k,1) +
                (1-iwf[0])*2.*vi0(j,k,1,1)-v3dtmp(j,k,1))+
                (1.-bci(j,k,1))*vi0(j,k,1,1);
            v3dtmp(j,k,idim) = bci(j,k,2)*(iwf[0]*v3dtmp(j,k,idim-1) +
                (1-iwf[0])*2.*vi0(j,k,1,3)-v3dtmp(j,k,idim-1))+
                (1.-bci(j,k,2))*vi0(j,k,1,3);
        }
    }

    // If 1st global subiteration for time-accurate computation, save zksav
    if ((float)dt > 0. && icyc == 1) {
        if (std::abs(ita) == 2) {
            if ((float)zksav2(1,1,1,1) == 0.) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    zksav2(j,k,i,4) = 0.;
                    zksav2(j,k,i,5) = 0.;
                    zksav2(j,k,i,6) = 0.;
                }
            } else {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    zksav2(j,k,i,4) = zksav(j,k,i,1)-zksav2(j,k,i,1);
                    zksav2(j,k,i,5) = zksav(j,k,i,2)-zksav2(j,k,i,2);
                    zksav2(j,k,i,6) = zksav(j,k,i,3)-zksav2(j,k,i,3);
                }
            }
        }
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            zksav2(j,k,i,1) = zksav(j,k,i,1);
            zksav2(j,k,i,2) = zksav(j,k,i,2);
            zksav2(j,k,i,3) = zksav(j,k,i,3);
        }
    }

    // Get TURRE values
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++)
    for (j = 1; j <= jdim-1; j++) {
        turre(j,k,i,1) = zksav(j,k,i,1);
        turre(j,k,i,2) = zksav(j,k,i,2);
        turre(j,k,i,3) = zksav(j,k,i,3);
    }



    // Iterate to solve the equations (do 500 not=1,nsubit)
    for (not_ = 1; not_ <= nsubit; not_++) {

        // Set up boundary conditions (ghost cells)
        // (1) k=0 boundary:
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            turre(j,0,i,1) = tk0(j,i,1,1);
            turre(j,0,i,2) = tk0(j,i,2,1);
            turre(j,0,i,3) = tk0(j,i,3,1);
        }
        // (2) k=kdim boundary:
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            turre(j,kdim,i,1) = tk0(j,i,1,3);
            turre(j,kdim,i,2) = tk0(j,i,2,3);
            turre(j,kdim,i,3) = tk0(j,i,3,3);
        }
        // (3) j=0 boundary:
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            turre(0,k,i,1) = tj0(k,i,1,1);
            turre(0,k,i,2) = tj0(k,i,2,1);
            turre(0,k,i,3) = tj0(k,i,3,1);
        }
        // (4) j=jdim boundary:
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            turre(jdim,k,i,1) = tj0(k,i,1,3);
            turre(jdim,k,i,2) = tj0(k,i,2,3);
            turre(jdim,k,i,3) = tj0(k,i,3,3);
        }
        if (i2d != 1 && iaxi2planeturb != 1) {
            // (5) i=0 boundary:
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,k,0,1) = ti0(j,k,1,1);
                turre(j,k,0,2) = ti0(j,k,2,1);
                turre(j,k,0,3) = ti0(j,k,3,1);
            }
            // (6) i=idim boundary:
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,k,idim,1) = ti0(j,k,1,3);
                turre(j,k,idim,2) = ti0(j,k,2,3);
                turre(j,k,idim,3) = ti0(j,k,3,3);
            }
        }
        if (iturbord != 1) {
            // (1) k=0 boundary (2nd ghost):
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,-1,i,1) = tk0(j,i,1,2);
                turre(j,-1,i,2) = tk0(j,i,2,2);
                turre(j,-1,i,3) = tk0(j,i,3,2);
            }
            // (2) k=kdim boundary (2nd ghost):
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,kdim+1,i,1) = tk0(j,i,1,4);
                turre(j,kdim+1,i,2) = tk0(j,i,2,4);
                turre(j,kdim+1,i,3) = tk0(j,i,3,4);
            }
            // (3) j=0 boundary (2nd ghost):
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                turre(-1,k,i,1) = tj0(k,i,1,2);
                turre(-1,k,i,2) = tj0(k,i,2,2);
                turre(-1,k,i,3) = tj0(k,i,3,2);
            }
            // (4) j=jdim boundary (2nd ghost):
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                turre(jdim+1,k,i,1) = tj0(k,i,1,4);
                turre(jdim+1,k,i,2) = tj0(k,i,2,4);
                turre(jdim+1,k,i,3) = tj0(k,i,3,4);
            }
            if (i2d != 1 && iaxi2planeturb != 1) {
                // (5) i=0 boundary (2nd ghost):
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    turre(j,k,-1,1) = ti0(j,k,1,2);
                    turre(j,k,-1,2) = ti0(j,k,2,2);
                    turre(j,k,-1,3) = ti0(j,k,3,2);
                }
                // (6) i=idim boundary (2nd ghost):
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    turre(j,k,idim+1,1) = ti0(j,k,1,4);
                    turre(j,k,idim+1,2) = ti0(j,k,2,4);
                    turre(j,k,idim+1,3) = ti0(j,k,3,4);
                }
            }
        }



        // Get damp1 = CD = cross derivative term for SST
        if (ivmx == 30 || ivmx == 31) {
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4) +
                    sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4) +
                    sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4) +
                    sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tt = xa*xa+ya*ya+za*za;
                ca = 2.*sigo2*tt/(turre(j,k,i,1)*re);
                damp1(j,k,i) = 0.25*ca*(turre(j,k+1,i,1)-turre(j,k-1,i,1))*
                                        (turre(j,k+1,i,2)-turre(j,k-1,i,2));
            }
            for (j = 1; j <= jdim-1; j++)
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4) +
                    sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4) +
                    sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4) +
                    sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tt = xa*xa+ya*ya+za*za;
                ca = 2.*sigo2*tt/(turre(j,k,i,1)*re);
                damp1(j,k,i) = damp1(j,k,i)+
                    0.25*ca*(turre(j+1,k,i,1)-turre(j-1,k,i,1))*
                             (turre(j+1,k,i,2)-turre(j-1,k,i,2));
            }
            if (i2d != 1 && iaxi2planeturb != 1) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4) +
                        si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4) +
                        si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4) +
                        si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tt = xa*xa+ya*ya+za*za;
                    ca = 2.*sigo2*tt/(turre(j,k,i,1)*re);
                    damp1(j,k,i) = damp1(j,k,i)+
                        0.25*ca*(turre(j,k,i+1,1)-turre(j,k,i-1,1))*
                                 (turre(j,k,i+1,2)-turre(j,k,i-1,2));
                }
            }
        }

        // Get blend = F1 factor
        if (ivmx == 30 || ivmx == 31) {
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                arg1 = std::sqrt(turre(j,k,i,2));
                { double a1t=smin(j,k,i); arg1 = arg1/(.09*re*turre(j,k,i,1)*ccomplex_ns::ccabs(a1t)); }
                arg2 = 500.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*re*re*
                    smin(j,k,i)*turre(j,k,i,1));
                { double a1t=arg1, a2t=arg2; arga=ccomplex_ns::ccmax(a1t,a2t); }
                { double a1t=damp1(j,k,i)*re, a2t=1.e-20; temp=ccomplex_ns::ccmaxcr(a1t,a2t); }
                argb = 4.*sigo2*turre(j,k,i,2)/(temp*smin(j,k,i)*smin(j,k,i));
                { double a1t=arga, a2t=argb; arg=ccomplex_ns::ccmin(a1t,a2t); }
                { double a1t=arg*arg*arg*arg; blend(j,k,i)=ccomplex_ns::cctanh(a1t); }
                if (itrans_on == 1) {
                    { double a1t=smin(j,k,i); re_nu = q(j,k,i,1)*ccomplex_ns::ccabs(a1t)*std::sqrt(turre(j,k,i,2))/fnu(j,k,i)*re; }
                    double fff3 = std::exp(-std::pow(re_nu/120., 8));
                    { double a1t=blend(j,k,i), a2t=fff3; blend(j,k,i)=ccomplex_ns::ccmax(a1t,a2t); }
                }
            }
        }





        // F_eta_eta viscous terms - Interior points
        for (k = 2; k <= kdim-2; k++) {
            kl = k-1; ku = k+1;
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
                dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                fnum = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,1) = -byy*turre(j,k-1,i,1)-cyy*turre(j,k,i,1)-dyy*turre(j,k+1,i,1);
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,2) = -byy*turre(j,k-1,i,2)-cyy*turre(j,k,i,2)-dyy*turre(j,k+1,i,2);
                cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,3) = -byy*turre(j,k-1,i,3)-cyy*turre(j,k,i,3)-dyy*turre(j,k+1,i,3);
            }
        }

        // K0 boundary points
        k = 1; kl = 1; ku = std::min(2, kdim-1);
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
            dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
            sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
            sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
            sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
            sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
            ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
            ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
            anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
            anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
            fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
            fnum = .5*(fnu(j,k-1,i)+fnu(j,k,i));
            cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,1) = -byy*turre(j,k-1,i,1)-cyy*turre(j,k,i,1)-dyy*turre(j,k+1,i,1);
            cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,2) = -byy*turre(j,k-1,i,2)-cyy*turre(j,k,i,2)-dyy*turre(j,k+1,i,2);
            cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,3) = -byy*turre(j,k-1,i,3)-cyy*turre(j,k,i,3)-dyy*turre(j,k+1,i,3);
        }

        // KDIM boundary points
        k = kdim-1; kl = kdim-2; ku = kdim-1;
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
            dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
            sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
            sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
            sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
            sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
            ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
            ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
            anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
            anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
            fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
            fnum = .5*(fnu(j,k-1,i)+fnu(j,k,i));
            cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,1) = -byy*turre(j,k-1,i,1)-cyy*turre(j,k,i,1)-dyy*turre(j,k+1,i,1);
            cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,2) = -byy*turre(j,k-1,i,2)-cyy*turre(j,k,i,2)-dyy*turre(j,k+1,i,2);
            cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,3) = -byy*turre(j,k-1,i,3)-cyy*turre(j,k,i,3)-dyy*turre(j,k+1,i,3);
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
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                rhside(j,k,i,1) = rhside(j,k,i,1)-uu*(app*(turre(j,k,i,1)-turre(j,k-1,i,1))+apm*(turre(j,k+1,i,1)-turre(j,k,i,1)));
                rhside(j,k,i,2) = rhside(j,k,i,2)-uu*(app*(turre(j,k,i,2)-turre(j,k-1,i,2))+apm*(turre(j,k+1,i,2)-turre(j,k,i,2)));
                rhside(j,k,i,3) = rhside(j,k,i,3)-uu*(app*(turre(j,k,i,3)-turre(j,k-1,i,3))+apm*(turre(j,k+1,i,3)-turre(j,k,i,3)));
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
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                rhside(j,k,i,1) = rhside(j,k,i,1)-0.5*uu*app*turre(j,k-2,i,1)
                    +2.*uu*app*turre(j,k-1,i,1)-1.5*uu*app*turre(j,k,i,1)
                    +1.5*uu*apm*turre(j,k,i,1)-2.*uu*apm*turre(j,k+1,i,1)
                    +0.5*uu*apm*turre(j,k+2,i,1);
                rhside(j,k,i,2) = rhside(j,k,i,2)-0.5*uu*app*turre(j,k-2,i,2)
                    +2.*uu*app*turre(j,k-1,i,2)-1.5*uu*app*turre(j,k,i,2)
                    +1.5*uu*apm*turre(j,k,i,2)-2.*uu*apm*turre(j,k+1,i,2)
                    +0.5*uu*apm*turre(j,k+2,i,2);
                rhside(j,k,i,3) = rhside(j,k,i,3)-0.5*uu*app*turre(j,k-2,i,3)
                    +2.*uu*app*turre(j,k-1,i,3)-1.5*uu*app*turre(j,k,i,3)
                    +1.5*uu*apm*turre(j,k,i,3)-2.*uu*apm*turre(j,k+1,i,3)
                    +0.5*uu*apm*turre(j,k+2,i,3);
            }
        }

        // F_xi_xi viscous terms - Interior points
        for (j = 2; j <= jdim-2; j++) {
            jl = j-1; ju = j+1;
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
                dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                fnum = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,1) = rhside(j,k,i,1)-bxx*turre(j-1,k,i,1)-cxx*turre(j,k,i,1)-dxx*turre(j+1,k,i,1);
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,2) = rhside(j,k,i,2)-bxx*turre(j-1,k,i,2)-cxx*turre(j,k,i,2)-dxx*turre(j+1,k,i,2);
                cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,3) = rhside(j,k,i,3)-bxx*turre(j-1,k,i,3)-cxx*turre(j,k,i,3)-dxx*turre(j+1,k,i,3);
            }
        }



        // J0 boundary points
        j = 1; jl = 1; ju = std::min(2, jdim-1);
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
            dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
            sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
            sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
            sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
            sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
            ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
            ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
            anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
            anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
            fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
            fnum = .5*(fnu(j-1,k,i)+fnu(j,k,i));
            cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,1) = rhside(j,k,i,1)-bxx*turre(j-1,k,i,1)-cxx*turre(j,k,i,1)-dxx*turre(j+1,k,i,1);
            cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,2) = rhside(j,k,i,2)-bxx*turre(j-1,k,i,2)-cxx*turre(j,k,i,2)-dxx*turre(j+1,k,i,2);
            cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,3) = rhside(j,k,i,3)-bxx*turre(j-1,k,i,3)-cxx*turre(j,k,i,3)-dxx*turre(j+1,k,i,3);
        }

        // JDIM boundary points
        j = jdim-1; jl = jdim-2; ju = jdim-1;
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
            dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
            sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
            sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
            sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
            sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
            ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
            ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
            anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
            anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
            fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
            fnum = .5*(fnu(j-1,k,i)+fnu(j,k,i));
            cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,1) = rhside(j,k,i,1)-bxx*turre(j-1,k,i,1)-cxx*turre(j,k,i,1)-dxx*turre(j+1,k,i,1);
            cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,2) = rhside(j,k,i,2)-bxx*turre(j-1,k,i,2)-cxx*turre(j,k,i,2)-dxx*turre(j+1,k,i,2);
            cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,3) = rhside(j,k,i,3)-bxx*turre(j-1,k,i,3)-cxx*turre(j,k,i,3)-dxx*turre(j+1,k,i,3);
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
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                rhside(j,k,i,1) = rhside(j,k,i,1)-uu*(app*(turre(j,k,i,1)-turre(j-1,k,i,1))+apm*(turre(j+1,k,i,1)-turre(j,k,i,1)));
                rhside(j,k,i,2) = rhside(j,k,i,2)-uu*(app*(turre(j,k,i,2)-turre(j-1,k,i,2))+apm*(turre(j+1,k,i,2)-turre(j,k,i,2)));
                rhside(j,k,i,3) = rhside(j,k,i,3)-uu*(app*(turre(j,k,i,3)-turre(j-1,k,i,3))+apm*(turre(j+1,k,i,3)-turre(j,k,i,3)));
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
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                rhside(j,k,i,1) = rhside(j,k,i,1)-0.5*uu*app*turre(j-2,k,i,1)
                    +2.*uu*app*turre(j-1,k,i,1)-1.5*uu*app*turre(j,k,i,1)
                    +1.5*uu*apm*turre(j,k,i,1)-2.*uu*apm*turre(j+1,k,i,1)
                    +0.5*uu*apm*turre(j+2,k,i,1);
                rhside(j,k,i,2) = rhside(j,k,i,2)-0.5*uu*app*turre(j-2,k,i,2)
                    +2.*uu*app*turre(j-1,k,i,2)-1.5*uu*app*turre(j,k,i,2)
                    +1.5*uu*apm*turre(j,k,i,2)-2.*uu*apm*turre(j+1,k,i,2)
                    +0.5*uu*apm*turre(j+2,k,i,2);
                rhside(j,k,i,3) = rhside(j,k,i,3)-0.5*uu*app*turre(j-2,k,i,3)
                    +2.*uu*app*turre(j-1,k,i,3)-1.5*uu*app*turre(j,k,i,3)
                    +1.5*uu*apm*turre(j,k,i,3)-2.*uu*apm*turre(j+1,k,i,3)
                    +0.5*uu*apm*turre(j+2,k,i,3);
            }
        }



        // F_zeta_zeta viscous terms
        if (i2d != 1 && iaxi2planeturb != 1) {
            // Interior points
            for (i = 2; i <= idim-2; i++) {
                il = i-1; iu = i+1;
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                    dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                    sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                    ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                    anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                    fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                    fnum = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                    cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,1) = rhside(j,k,i,1)-bzz*turre(j,k,i-1,1)-czz*turre(j,k,i,1)-dzz*turre(j,k,i+1,1);
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,2) = rhside(j,k,i,2)-bzz*turre(j,k,i-1,2)-czz*turre(j,k,i,2)-dzz*turre(j,k,i+1,2);
                    cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,3) = rhside(j,k,i,3)-bzz*turre(j,k,i-1,3)-czz*turre(j,k,i,3)-dzz*turre(j,k,i+1,3);
                }
            }

            // I0 boundary points
            i = 1; il = 1; iu = std::min(2, idim-1);
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                fnum = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,1) = rhside(j,k,i,1)-bzz*turre(j,k,i-1,1)-czz*turre(j,k,i,1)-dzz*turre(j,k,i+1,1);
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,2) = rhside(j,k,i,2)-bzz*turre(j,k,i-1,2)-czz*turre(j,k,i,2)-dzz*turre(j,k,i+1,2);
                cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,3) = rhside(j,k,i,3)-bzz*turre(j,k,i-1,3)-czz*turre(j,k,i,3)-dzz*turre(j,k,i+1,3);
            }

            // IDIM boundary points
            i = idim-1; il = idim-2; iu = idim-1;
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                fnum = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,1) = rhside(j,k,i,1)-bzz*turre(j,k,i-1,1)-czz*turre(j,k,i,1)-dzz*turre(j,k,i+1,1);
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,2) = rhside(j,k,i,2)-bzz*turre(j,k,i-1,2)-czz*turre(j,k,i,2)-dzz*turre(j,k,i+1,2);
                cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,3) = rhside(j,k,i,3)-bzz*turre(j,k,i-1,3)-czz*turre(j,k,i,3)-dzz*turre(j,k,i+1,3);
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
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                    app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                    rhside(j,k,i,1) = rhside(j,k,i,1)-uu*(app*(turre(j,k,i,1)-turre(j,k,i-1,1))+apm*(turre(j,k,i+1,1)-turre(j,k,i,1)));
                    rhside(j,k,i,2) = rhside(j,k,i,2)-uu*(app*(turre(j,k,i,2)-turre(j,k,i-1,2))+apm*(turre(j,k,i+1,2)-turre(j,k,i,2)));
                    rhside(j,k,i,3) = rhside(j,k,i,3)-uu*(app*(turre(j,k,i,3)-turre(j,k,i-1,3))+apm*(turre(j,k,i+1,3)-turre(j,k,i,3)));
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
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                    app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                    rhside(j,k,i,1) = rhside(j,k,i,1)-0.5*uu*app*turre(j,k,i-2,1)
                        +2.*uu*app*turre(j,k,i-1,1)-1.5*uu*app*turre(j,k,i,1)
                        +1.5*uu*apm*turre(j,k,i,1)-2.*uu*apm*turre(j,k,i+1,1)
                        +0.5*uu*apm*turre(j,k,i+2,1);
                    rhside(j,k,i,2) = rhside(j,k,i,2)-0.5*uu*app*turre(j,k,i-2,2)
                        +2.*uu*app*turre(j,k,i-1,2)-1.5*uu*app*turre(j,k,i,2)
                        +1.5*uu*apm*turre(j,k,i,2)-2.*uu*apm*turre(j,k,i+1,2)
                        +0.5*uu*apm*turre(j,k,i+2,2);
                    rhside(j,k,i,3) = rhside(j,k,i,3)-0.5*uu*app*turre(j,k,i-2,3)
                        +2.*uu*app*turre(j,k,i-1,3)-1.5*uu*app*turre(j,k,i,3)
                        +1.5*uu*apm*turre(j,k,i,3)-2.*uu*apm*turre(j,k,i+1,3)
                        +0.5*uu*apm*turre(j,k,i+2,3);
                }
            }
        } // end if (i2d != 1 && iaxi2planeturb != 1) for F_zeta_zeta



        // ADD SOURCE TERMS TO RHS
        if (isklton > 0) {
            if (ilamlo == 0 || jlamlo == 0 || klamlo == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " Block #%5d in 3-eqn turb model has no laminar regions", nbl);
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " Block #%5d in 3-eqn turb model - laminar region is:", nbl);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " i=%5d to%5d, j=%5d to%5d, k=%5d to%5d",
                    ilamlo, ilamhi, jlamlo, jlamhi, klamlo, klamhi);
                if (i_lam_forcezero == 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "    ...forcing vist3d=0");
                }
            }
        }

        // Source terms for ivmx==30
        if (ivmx == 30) {
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
                betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                cmuc  = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                alp   = blend(j,k,i)*alp1 +(1.-blend(j,k,i))*alp2;
                s11 = ux(j,k,i,1); s22 = ux(j,k,i,5); s33 = ux(j,k,i,9);
                s12 = 0.5*(ux(j,k,i,2)+ux(j,k,i,4));
                s13 = 0.5*(ux(j,k,i,3)+ux(j,k,i,7));
                s23 = 0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                tracepart = (s11+s22+s33)*(double)prod2d3dtrace;
                s11t = s11-tracepart; s22t = s22-tracepart; s33t = s33-tracepart;
                xis = s11t*s11t+s22t*s22t+s33t*s33t+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                if (ikoprod == 1) {
                    rhside(j,k,i,1) = rhside(j,k,i,1)+cutoff*alp/re*2.*xis
                        -re*betax*turre(j,k,i,1)*turre(j,k,i,1)
                        +(1.-blend(j,k,i))*damp1(j,k,i);
                    pk = vist3d(j,k,i)/(q(j,k,i,1)*re)*2.*xis;
                } else {
                    rhside(j,k,i,1) = rhside(j,k,i,1)+cutoff*alp/re*
                        vor(j,k,i)*vor(j,k,i)-re*betax*turre(j,k,i,1)*turre(j,k,i,1)
                        +(1.-blend(j,k,i))*damp1(j,k,i);
                    pk = vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)*vor(j,k,i);
                }
                dk = re*cmuc*turre(j,k,i,1)*turre(j,k,i,2);
                { double a1t=pk, a2t=(double)pklimterm*dk; pk=ccomplex_ns::ccmin(a1t,a2t); }
                { double a1t=turre(j,k,i,3), a2t=0.1; dkfactor=ccomplex_ns::ccmaxcr(a1t,a2t); }
                { double a1t=dkfactor, a2t=1.0; rhside(j,k,i,2) = rhside(j,k,i,2)+cutoff*pk*turre(j,k,i,3)-dk*ccomplex_ns::ccmincr(a1t,a2t); }
                rhside(j,k,i,1) = rhside(j,k,i,1)+(double)keepambient*re*betax*(double)tur10[0]*(double)tur10[0];
                rhside(j,k,i,2) = rhside(j,k,i,2)+(double)keepambient*re*cmuc*(double)tur10[0]*(double)tur10[1];
                { double a1t=smin(j,k,i); re_nu = q(j,k,i,1)*smin(j,k,i)*smin(j,k,i)*vor(j,k,i)*re/fnu(j,k,i); }
                re_t = q(j,k,i,1)*turre(j,k,i,2)/(fnu(j,k,i)*turre(j,k,i,1));
                fonset1 = re_nu/(2.193*re_thetac);
                { double a1t=fonset1, a2t=fonset1*fonset1*fonset1*fonset1; fonset2=ccomplex_ns::ccmax(a1t,a2t); }
                { double a1t=fonset2, a2t=2.0; fonset2=ccomplex_ns::ccmincr(a1t,a2t); }
                { double a1t=re_nu/(2.193*0.35*300.); ftrnlength=ccomplex_ns::cctanh(a1t); }
                trnlength = (0.014*re_thetac+1.25)*ftrnlength+2.*(1.-ftrnlength);
                { double a1t=(1.-re_t/trnlength), a2t=0.; fonset3=ccomplex_ns::ccmaxcr(a1t,a2t); }
                { double a1t=fonset2-fonset3, a2t=0.; fonset=ccomplex_ns::ccmaxcr(a1t,a2t); }
                fturb = std::exp(-std::pow(re_t/2., 2));
                { double a1t=turre(j,k,i,3), a2t=1.0; double tmp=ccomplex_ns::ccmincr(a1t,a2t);
                  double a1t2=turre(j,k,i,3)/0.35; double tanh_val=ccomplex_ns::cctanh(a1t2);
                  rhside(j,k,i,3) = rhside(j,k,i,3)+
                    5.*std::sqrt(2.*xis)*fonset*(1.0-tmp)-
                    0.1*vor(j,k,i)*fturb*tanh_val; }
            }
        } else if (ivmx == 31) {
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
                betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                cmuc  = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                alp   = blend(j,k,i)*alp1 +(1.-blend(j,k,i))*alp2;
                re_thetat_tilde = re_thetat;
                { double a1t=smin(j,k,i); re_nu = q(j,k,i,1)*smin(j,k,i)*smin(j,k,i)*vor(j,k,i)*re/fnu(j,k,i); }
                re_t = q(j,k,i,1)*turre(j,k,i,2)/(fnu(j,k,i)*turre(j,k,i,1));
                uuu = std::sqrt(q(j,k,i,2)*q(j,k,i,2)+q(j,k,i,3)*q(j,k,i,3)+q(j,k,i,4)*q(j,k,i,4));
                { double a1t=uuu, a2t=1.e-20; uuu=ccomplex_ns::ccmaxcr(a1t,a2t); }
                theta_bl = re_thetat_tilde*fnu(j,k,i)/(q(j,k,i,1)*uuu*re);
                delta_bl = 15.*theta_bl/2.;
                { double a1t=smin(j,k,i); delta = 50.*vor(j,k,i)*ccomplex_ns::ccabs(a1t)*delta_bl/uuu; }
                re_omega = q(j,k,i,1)*turre(j,k,i,1)*smin(j,k,i)*smin(j,k,i)/fnu(j,k,i)*re*re;
                fwake = std::exp(-std::pow(re_omega/1.e5, 2));
                { double a1t=smin(j,k,i); f_thetat1 = fwake*std::exp(-std::pow(ccomplex_ns::ccabs(a1t)/delta, 4)); }
                { double a1t=turre(j,k,i,3), a2t=1.0; gamma_lim=ccomplex_ns::ccmincr(a1t,a2t); }
                f_thetat2 = 1.-std::pow((ce2*gamma_lim-1.)/(ce2-1.), 2);
                { double a1t=f_thetat1, a2t=f_thetat2; f_thetat=ccomplex_ns::ccmax(a1t,a2t); }
                { double a1t=f_thetat, a2t=1.0; f_thetat=ccomplex_ns::ccmincr(a1t,a2t); }
                f_reattach = std::exp(-std::pow(re_t/15., 4));
                fonset1 = re_nu/(2.193*re_thetac);
                { double a1t=fonset1-1., a2t=0.; gamma_sep1=ccomplex_ns::ccmaxcr(a1t,a2t); }
                { double a1t=8.0*gamma_sep1*f_reattach, a2t=5.0; gamma_sep=ccomplex_ns::ccmincr(a1t,a2t); }
                gamma_sep = gamma_sep*f_thetat;
                { double a1t=turre(j,k,i,3), a2t=gamma_sep; gamma_eff=ccomplex_ns::ccmax(a1t,a2t); }
                s11 = ux(j,k,i,1); s22 = ux(j,k,i,5); s33 = ux(j,k,i,9);
                s12 = 0.5*(ux(j,k,i,2)+ux(j,k,i,4));
                s13 = 0.5*(ux(j,k,i,3)+ux(j,k,i,7));
                s23 = 0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                tracepart = (s11+s22+s33)*(double)prod2d3dtrace;
                s11t = s11-tracepart; s22t = s22-tracepart; s33t = s33-tracepart;
                xis = s11t*s11t+s22t*s22t+s33t*s33t+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                if (ikoprod == 1) {
                    rhside(j,k,i,1) = rhside(j,k,i,1)+cutoff*alp/re*2.*xis
                        -re*betax*turre(j,k,i,1)*turre(j,k,i,1)
                        +(1.-blend(j,k,i))*damp1(j,k,i);
                    pk = vist3d(j,k,i)/(q(j,k,i,1)*re)*2.*xis;
                } else {
                    rhside(j,k,i,1) = rhside(j,k,i,1)+cutoff*alp/re*
                        vor(j,k,i)*vor(j,k,i)-re*betax*turre(j,k,i,1)*turre(j,k,i,1)
                        +(1.-blend(j,k,i))*damp1(j,k,i);
                    pk = vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)*vor(j,k,i);
                }
                dk = re*cmuc*turre(j,k,i,1)*turre(j,k,i,2);
                { double a1t=pk, a2t=(double)pklimterm*dk; pk=ccomplex_ns::ccmin(a1t,a2t); }
                { double a1t=gamma_eff, a2t=0.1; dkfactor=ccomplex_ns::ccmaxcr(a1t,a2t); }
                { double a1t=dkfactor, a2t=1.0; rhside(j,k,i,2) = rhside(j,k,i,2)+cutoff*pk*gamma_eff-dk*ccomplex_ns::ccmincr(a1t,a2t); }
                rhside(j,k,i,1) = rhside(j,k,i,1)+(double)keepambient*re*betax*(double)tur10[0]*(double)tur10[0];
                rhside(j,k,i,2) = rhside(j,k,i,2)+(double)keepambient*re*cmuc*(double)tur10[0]*(double)tur10[1];
                { double a1t=fonset1, a2t=fonset1*fonset1*fonset1*fonset1; fonset2=ccomplex_ns::ccmax(a1t,a2t); }
                { double a1t=fonset2, a2t=2.0; fonset2=ccomplex_ns::ccmincr(a1t,a2t); }
                { double a1t=(1.-std::pow(re_t/2.5,3)), a2t=0.; fonset3=ccomplex_ns::ccmaxcr(a1t,a2t); }
                { double a1t=fonset2-fonset3, a2t=0.; fonset=ccomplex_ns::ccmaxcr(a1t,a2t); }
                fturb = std::exp(-std::pow(re_t/4., 4));
                rhside(j,k,i,3) = rhside(j,k,i,3)+
                    (flength*2.0*std::sqrt(2.*xis)*std::pow(gamma_lim*fonset, 0.5))*
                    (1.0-gamma_lim)-
                    (0.06*vor(j,k,i)*gamma_lim*fturb*(ce2*gamma_lim-1.0));
            }
        }




        // Implicit F_eta_eta viscous terms. Do over all i's
        for (i = 1; i <= idim-1; i++) {
            // Build by/cy/dy/fy for all k values
            // Interior points
            for (k = 2; k <= kdim-2; k++) {
                kl = k-1; ku = k+1;
                for (j = 1; j <= jdim-1; j++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
                    dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
                    sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                    ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                    anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                    fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                    fnum = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                    cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    by(j,k) = -cdm; cy(j,k) = cdp+cdm; dy(j,k) = -cdp;
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    by2(j,k) = -cdm; cy2(j,k) = cdp+cdm; dy2(j,k) = -cdp;
                    cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                    by3(j,k) = -cdm; cy3(j,k) = cdp+cdm; dy3(j,k) = -cdp;
                }
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                    app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                    by(j,k) = by(j,k)-uu*app; by2(j,k) = by2(j,k)-uu*app; by3(j,k) = by3(j,k)-uu*app;
                    cy(j,k) = cy(j,k)+uu*(app-apm); cy2(j,k) = cy2(j,k)+uu*(app-apm); cy3(j,k) = cy3(j,k)+uu*(app-apm);
                    dy(j,k) = dy(j,k)+uu*apm; dy2(j,k) = dy2(j,k)+uu*apm; dy3(j,k) = dy3(j,k)+uu*apm;
                }
                for (j = 1; j <= jdim-1; j++) {
                    double cyadd = 0., cy2add = 0.;
                    if (ivmx == 30 || ivmx == 31) {
                        cmuc  = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                        betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                        { double a1t=damp1(j,k,i); cyadd = 2.*re*betax*turre(j,k,i,1)+ccomplex_ns::ccabs(a1t)/turre(j,k,i,1); }
                        cy2add = re*cmuc*turre(j,k,i,1);
                    }
                    cy(j,k) = cy(j,k)+cyadd; cy2(j,k) = cy2(j,k)+cy2add;
                }
                for (j = 1; j <= jdim-1; j++) {
                    by(j,k) = by(j,k)*timestp(j,k,i);
                    by2(j,k) = by2(j,k)*timestp(j,k,i)*factor2;
                    by3(j,k) = by3(j,k)*timestp(j,k,i)*factor3;
                    cy(j,k) = cy(j,k)*timestp(j,k,i)+1.0*(1.+phi);
                    cy2(j,k) = cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                    cy3(j,k) = cy3(j,k)*timestp(j,k,i)*factor3+1.0*(1.+phi);
                    dy(j,k) = dy(j,k)*timestp(j,k,i);
                    dy2(j,k) = dy2(j,k)*timestp(j,k,i)*factor2;
                    dy3(j,k) = dy3(j,k)*timestp(j,k,i)*factor3;
                    fy(j,k) = rhside(j,k,i,1)*timestp(j,k,i);
                    fy2(j,k) = rhside(j,k,i,2)*timestp(j,k,i)*factor2;
                    fy3(j,k) = rhside(j,k,i,3)*timestp(j,k,i)*factor3;
                }
                if ((float)dt > 0.) {
                    for (j = 1; j <= jdim-1; j++) {
                        fy(j,k) = fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-turre(j,k,i,1))+phi*zksav2(j,k,i,4);
                        fy2(j,k) = fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-turre(j,k,i,2))+phi*zksav2(j,k,i,5);
                        fy3(j,k) = fy3(j,k)+(1.+phi)*(zksav2(j,k,i,3)-turre(j,k,i,3))+phi*zksav2(j,k,i,6);
                    }
                }
            } // end interior k loop

            // K0 boundary points
            k = 1; kl = 1; ku = std::min(2, kdim-1);
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
                dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                fnum = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                by(j,k) = -cdm; cy(j,k) = cdp+cdm; dy(j,k) = -cdp;
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                by2(j,k) = -cdm; cy2(j,k) = cdp+cdm; dy2(j,k) = -cdp;
                cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                by3(j,k) = -cdm; cy3(j,k) = cdp+cdm; dy3(j,k) = -cdp;
            }
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                by(j,k) = by(j,k)-uu*app; by2(j,k) = by2(j,k)-uu*app; by3(j,k) = by3(j,k)-uu*app;
                cy(j,k) = cy(j,k)+uu*(app-apm); cy2(j,k) = cy2(j,k)+uu*(app-apm); cy3(j,k) = cy3(j,k)+uu*(app-apm);
                dy(j,k) = dy(j,k)+uu*apm; dy2(j,k) = dy2(j,k)+uu*apm; dy3(j,k) = dy3(j,k)+uu*apm;
            }
            for (j = 1; j <= jdim-1; j++) {
                double cyadd = 0., cy2add = 0.;
                if (ivmx == 30 || ivmx == 31) {
                    cmuc  = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                    betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                    { double a1t=damp1(j,k,i); cyadd = 2.*re*betax*turre(j,k,i,1)+ccomplex_ns::ccabs(a1t)/turre(j,k,i,1); }
                    cy2add = re*cmuc*turre(j,k,i,1);
                }
                cy(j,k) = cy(j,k)+cyadd; cy2(j,k) = cy2(j,k)+cy2add;
            }
            for (j = 1; j <= jdim-1; j++) {
                by(j,k) = by(j,k)*timestp(j,k,i);
                by2(j,k) = by2(j,k)*timestp(j,k,i)*factor2;
                by3(j,k) = by3(j,k)*timestp(j,k,i)*factor3;
                cy(j,k) = cy(j,k)*timestp(j,k,i)+1.0*(1.+phi);
                cy2(j,k) = cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                cy3(j,k) = cy3(j,k)*timestp(j,k,i)*factor3+1.0*(1.+phi);
                dy(j,k) = dy(j,k)*timestp(j,k,i);
                dy2(j,k) = dy2(j,k)*timestp(j,k,i)*factor2;
                dy3(j,k) = dy3(j,k)*timestp(j,k,i)*factor3;
                fy(j,k) = rhside(j,k,i,1)*timestp(j,k,i);
                fy2(j,k) = rhside(j,k,i,2)*timestp(j,k,i)*factor2;
                fy3(j,k) = rhside(j,k,i,3)*timestp(j,k,i)*factor3;
            }
            if ((float)dt > 0.) {
                for (j = 1; j <= jdim-1; j++) {
                    fy(j,k) = fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-turre(j,k,i,1))+phi*zksav2(j,k,i,4);
                    fy2(j,k) = fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-turre(j,k,i,2))+phi*zksav2(j,k,i,5);
                    fy3(j,k) = fy3(j,k)+(1.+phi)*(zksav2(j,k,i,3)-turre(j,k,i,3))+phi*zksav2(j,k,i,6);
                }
            }

            // KDIM boundary points
            k = kdim-1; kl = kdim-2; ku = kdim-1;
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
                dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                fnum = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                by(j,k) = -cdm; cy(j,k) = cdp+cdm; dy(j,k) = -cdp;
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                by2(j,k) = -cdm; cy2(j,k) = cdp+cdm; dy2(j,k) = -cdp;
                cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                by3(j,k) = -cdm; cy3(j,k) = cdp+cdm; dy3(j,k) = -cdp;
            }
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                by(j,k) = by(j,k)-uu*app; by2(j,k) = by2(j,k)-uu*app; by3(j,k) = by3(j,k)-uu*app;
                cy(j,k) = cy(j,k)+uu*(app-apm); cy2(j,k) = cy2(j,k)+uu*(app-apm); cy3(j,k) = cy3(j,k)+uu*(app-apm);
                dy(j,k) = dy(j,k)+uu*apm; dy2(j,k) = dy2(j,k)+uu*apm; dy3(j,k) = dy3(j,k)+uu*apm;
            }
            for (j = 1; j <= jdim-1; j++) {
                double cyadd = 0., cy2add = 0.;
                if (ivmx == 30 || ivmx == 31) {
                    cmuc  = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                    betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                    { double a1t=damp1(j,k,i); cyadd = 2.*re*betax*turre(j,k,i,1)+ccomplex_ns::ccabs(a1t)/turre(j,k,i,1); }
                    cy2add = re*cmuc*turre(j,k,i,1);
                }
                cy(j,k) = cy(j,k)+cyadd; cy2(j,k) = cy2(j,k)+cy2add;
            }
            for (j = 1; j <= jdim-1; j++) {
                by(j,k) = by(j,k)*timestp(j,k,i);
                by2(j,k) = by2(j,k)*timestp(j,k,i)*factor2;
                by3(j,k) = by3(j,k)*timestp(j,k,i)*factor3;
                cy(j,k) = cy(j,k)*timestp(j,k,i)+1.0*(1.+phi);
                cy2(j,k) = cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                cy3(j,k) = cy3(j,k)*timestp(j,k,i)*factor3+1.0*(1.+phi);
                dy(j,k) = dy(j,k)*timestp(j,k,i);
                dy2(j,k) = dy2(j,k)*timestp(j,k,i)*factor2;
                dy3(j,k) = dy3(j,k)*timestp(j,k,i)*factor3;
                fy(j,k) = rhside(j,k,i,1)*timestp(j,k,i);
                fy2(j,k) = rhside(j,k,i,2)*timestp(j,k,i)*factor2;
                fy3(j,k) = rhside(j,k,i,3)*timestp(j,k,i)*factor3;
            }
            if ((float)dt > 0.) {
                for (j = 1; j <= jdim-1; j++) {
                    fy(j,k) = fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-turre(j,k,i,1))+phi*zksav2(j,k,i,4);
                    fy2(j,k) = fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-turre(j,k,i,2))+phi*zksav2(j,k,i,5);
                    fy3(j,k) = fy3(j,k)+(1.+phi)*(zksav2(j,k,i,3)-turre(j,k,i,3))+phi*zksav2(j,k,i,6);
                }
            }

            // Apply blank (overset)
            if (iover == 1) {
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    fy(j,k) = fy(j,k)*blank(j,k,i);
                    by(j,k) = by(j,k)*blank(j,k,i);
                    dy(j,k) = dy(j,k)*blank(j,k,i);
                    cy(j,k) = cy(j,k)*blank(j,k,i)+(1.-blank(j,k,i));
                    fy2(j,k) = fy2(j,k)*blank(j,k,i);
                    by2(j,k) = by2(j,k)*blank(j,k,i);
                    dy2(j,k) = dy2(j,k)*blank(j,k,i);
                    cy2(j,k) = cy2(j,k)*blank(j,k,i)+(1.-blank(j,k,i));
                    fy3(j,k) = fy3(j,k)*blank(j,k,i);
                    by3(j,k) = by3(j,k)*blank(j,k,i);
                    dy3(j,k) = dy3(j,k)*blank(j,k,i);
                    cy3(j,k) = cy3(j,k)*blank(j,k,i)+(1.-blank(j,k,i));
                }
            }

            // Solve tridiagonal systems in k direction
            { int jd1=jdim-1, kd1=kdim-1, one1=1, one2=1;
              triv_ns::triv(jd1, kd1, one1, jd1, one2, kd1, worky, by, cy, dy, fy);
              triv_ns::triv(jd1, kd1, one1, jd1, one2, kd1, worky, by2, cy2, dy2, fy2);
              triv_ns::triv(jd1, kd1, one1, jd1, one2, kd1, worky, by3, cy3, dy3, fy3); }

            // Store result back in rhside
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                rhside(j,k,i,1) = fy(j,k);
                rhside(j,k,i,2) = fy2(j,k);
                rhside(j,k,i,3) = fy3(j,k);
            }
        } // end i loop for k sweep

        // Implicit F_xi_xi viscous terms. Do over all i's
        for (i = 1; i <= idim-1; i++) {
            // Interior points
            for (j = 2; j <= jdim-2; j++) {
                jl = j-1; ju = j+1;
                for (k = 1; k <= kdim-1; k++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
                    dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
                    sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                    ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                    anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                    fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                    fnum = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                    cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    bx(k,j) = -cdm; cx(k,j) = cdp+cdm; dx(k,j) = -cdp;
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    bx2(k,j) = -cdm; cx2(k,j) = cdp+cdm; dx2(k,j) = -cdp;
                    cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                    bx3(k,j) = -cdm; cx3(k,j) = cdp+cdm; dx3(k,j) = -cdp;
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                    app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                    bx(k,j) = bx(k,j)-uu*app; bx2(k,j) = bx2(k,j)-uu*app; bx3(k,j) = bx3(k,j)-uu*app;
                    cx(k,j) = cx(k,j)+uu*(app-apm); cx2(k,j) = cx2(k,j)+uu*(app-apm); cx3(k,j) = cx3(k,j)+uu*(app-apm);
                    dx(k,j) = dx(k,j)+uu*apm; dx2(k,j) = dx2(k,j)+uu*apm; dx3(k,j) = dx3(k,j)+uu*apm;
                }
                for (k = 1; k <= kdim-1; k++) {
                    double cxadd = 0., cx2add = 0.;
                    if (ivmx == 30 || ivmx == 31) {
                        cmuc  = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                        betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                        { double a1t=damp1(j,k,i); cxadd = 2.*re*betax*turre(j,k,i,1)+ccomplex_ns::ccabs(a1t)/turre(j,k,i,1); }
                        cx2add = re*cmuc*turre(j,k,i,1);
                    }
                    cx(k,j) = cx(k,j)+cxadd; cx2(k,j) = cx2(k,j)+cx2add;
                }
                for (k = 1; k <= kdim-1; k++) {
                    bx(k,j) = bx(k,j)*timestp(j,k,i);
                    bx2(k,j) = bx2(k,j)*timestp(j,k,i)*factor2;
                    bx3(k,j) = bx3(k,j)*timestp(j,k,i)*factor3;
                    cx(k,j) = cx(k,j)*timestp(j,k,i)+1.0*(1.+phi);
                    cx2(k,j) = cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                    cx3(k,j) = cx3(k,j)*timestp(j,k,i)*factor3+1.0*(1.+phi);
                    dx(k,j) = dx(k,j)*timestp(j,k,i);
                    dx2(k,j) = dx2(k,j)*timestp(j,k,i)*factor2;
                    dx3(k,j) = dx3(k,j)*timestp(j,k,i)*factor3;
                    fx(k,j) = rhside(j,k,i,1)*(1.+phi);
                    fx2(k,j) = rhside(j,k,i,2)*(1.+phi);
                    fx3(k,j) = rhside(j,k,i,3)*(1.+phi);
                }
            } // end interior j loop

            // J0 boundary points
            j = 1; jl = 1; ju = std::min(2, jdim-1);
            for (k = 1; k <= kdim-1; k++) {
                dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
                dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                fnum = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                bx(k,j) = -cdm; cx(k,j) = cdp+cdm; dx(k,j) = -cdp;
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                bx2(k,j) = -cdm; cx2(k,j) = cdp+cdm; dx2(k,j) = -cdp;
                cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                bx3(k,j) = -cdm; cx3(k,j) = cdp+cdm; dx3(k,j) = -cdp;
            }
            for (k = 1; k <= kdim-1; k++) {
                xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                bx(k,j) = bx(k,j)-uu*app; bx2(k,j) = bx2(k,j)-uu*app; bx3(k,j) = bx3(k,j)-uu*app;
                cx(k,j) = cx(k,j)+uu*(app-apm); cx2(k,j) = cx2(k,j)+uu*(app-apm); cx3(k,j) = cx3(k,j)+uu*(app-apm);
                dx(k,j) = dx(k,j)+uu*apm; dx2(k,j) = dx2(k,j)+uu*apm; dx3(k,j) = dx3(k,j)+uu*apm;
            }
            for (k = 1; k <= kdim-1; k++) {
                double cxadd = 0., cx2add = 0.;
                if (ivmx == 30 || ivmx == 31) {
                    cmuc  = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                    betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                    { double a1t=damp1(j,k,i); cxadd = 2.*re*betax*turre(j,k,i,1)+ccomplex_ns::ccabs(a1t)/turre(j,k,i,1); }
                    cx2add = re*cmuc*turre(j,k,i,1);
                }
                cx(k,j) = cx(k,j)+cxadd; cx2(k,j) = cx2(k,j)+cx2add;
            }
            for (k = 1; k <= kdim-1; k++) {
                bx(k,j) = bx(k,j)*timestp(j,k,i);
                bx2(k,j) = bx2(k,j)*timestp(j,k,i)*factor2;
                bx3(k,j) = bx3(k,j)*timestp(j,k,i)*factor3;
                cx(k,j) = cx(k,j)*timestp(j,k,i)+1.0*(1.+phi);
                cx2(k,j) = cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                cx3(k,j) = cx3(k,j)*timestp(j,k,i)*factor3+1.0*(1.+phi);
                dx(k,j) = dx(k,j)*timestp(j,k,i);
                dx2(k,j) = dx2(k,j)*timestp(j,k,i)*factor2;
                dx3(k,j) = dx3(k,j)*timestp(j,k,i)*factor3;
                fx(k,j) = rhside(j,k,i,1)*(1.+phi);
                fx2(k,j) = rhside(j,k,i,2)*(1.+phi);
                fx3(k,j) = rhside(j,k,i,3)*(1.+phi);
            }

            // JDIM boundary points
            j = jdim-1; jl = jdim-2; ju = jdim-1;
            for (k = 1; k <= kdim-1; k++) {
                dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
                dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                fnum = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                bx(k,j) = -cdm; cx(k,j) = cdp+cdm; dx(k,j) = -cdp;
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                bx2(k,j) = -cdm; cx2(k,j) = cdp+cdm; dx2(k,j) = -cdp;
                cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                bx3(k,j) = -cdm; cx3(k,j) = cdp+cdm; dx3(k,j) = -cdp;
            }
            for (k = 1; k <= kdim-1; k++) {
                xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                bx(k,j) = bx(k,j)-uu*app; bx2(k,j) = bx2(k,j)-uu*app; bx3(k,j) = bx3(k,j)-uu*app;
                cx(k,j) = cx(k,j)+uu*(app-apm); cx2(k,j) = cx2(k,j)+uu*(app-apm); cx3(k,j) = cx3(k,j)+uu*(app-apm);
                dx(k,j) = dx(k,j)+uu*apm; dx2(k,j) = dx2(k,j)+uu*apm; dx3(k,j) = dx3(k,j)+uu*apm;
            }
            for (k = 1; k <= kdim-1; k++) {
                double cxadd = 0., cx2add = 0.;
                if (ivmx == 30 || ivmx == 31) {
                    cmuc  = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                    betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                    { double a1t=damp1(j,k,i); cxadd = 2.*re*betax*turre(j,k,i,1)+ccomplex_ns::ccabs(a1t)/turre(j,k,i,1); }
                    cx2add = re*cmuc*turre(j,k,i,1);
                }
                cx(k,j) = cx(k,j)+cxadd; cx2(k,j) = cx2(k,j)+cx2add;
            }
            for (k = 1; k <= kdim-1; k++) {
                bx(k,j) = bx(k,j)*timestp(j,k,i);
                bx2(k,j) = bx2(k,j)*timestp(j,k,i)*factor2;
                bx3(k,j) = bx3(k,j)*timestp(j,k,i)*factor3;
                cx(k,j) = cx(k,j)*timestp(j,k,i)+1.0*(1.+phi);
                cx2(k,j) = cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                cx3(k,j) = cx3(k,j)*timestp(j,k,i)*factor3+1.0*(1.+phi);
                dx(k,j) = dx(k,j)*timestp(j,k,i);
                dx2(k,j) = dx2(k,j)*timestp(j,k,i)*factor2;
                dx3(k,j) = dx3(k,j)*timestp(j,k,i)*factor3;
                fx(k,j) = rhside(j,k,i,1)*(1.+phi);
                fx2(k,j) = rhside(j,k,i,2)*(1.+phi);
                fx3(k,j) = rhside(j,k,i,3)*(1.+phi);
            }

            // Apply blank (overset)
            if (iover == 1) {
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    fx(k,j) = fx(k,j)*blank(j,k,i);
                    bx(k,j) = bx(k,j)*blank(j,k,i);
                    dx(k,j) = dx(k,j)*blank(j,k,i);
                    cx(k,j) = cx(k,j)*blank(j,k,i)+(1.-blank(j,k,i));
                    fx2(k,j) = fx2(k,j)*blank(j,k,i);
                    bx2(k,j) = bx2(k,j)*blank(j,k,i);
                    dx2(k,j) = dx2(k,j)*blank(j,k,i);
                    cx2(k,j) = cx2(k,j)*blank(j,k,i)+(1.-blank(j,k,i));
                    fx3(k,j) = fx3(k,j)*blank(j,k,i);
                    bx3(k,j) = bx3(k,j)*blank(j,k,i);
                    dx3(k,j) = dx3(k,j)*blank(j,k,i);
                    cx3(k,j) = cx3(k,j)*blank(j,k,i)+(1.-blank(j,k,i));
                }
            }

            // Solve tridiagonal systems in j direction
            { int kd1=kdim-1, jd1=jdim-1, one1=1, one2=1;
              triv_ns::triv(kd1, jd1, one1, kd1, one2, jd1, workx, bx, cx, dx, fx);
              triv_ns::triv(kd1, jd1, one1, kd1, one2, jd1, workx, bx2, cx2, dx2, fx2);
              triv_ns::triv(kd1, jd1, one1, kd1, one2, jd1, workx, bx3, cx3, dx3, fx3); }

            // Store result back in rhside
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                rhside(j,k,i,1) = fx(k,j);
                rhside(j,k,i,2) = fx2(k,j);
                rhside(j,k,i,3) = fx3(k,j);
            }
        } // end i loop for j sweep

        // Implicit F_zeta_zeta viscous terms. Do over all j's
        if (i2d != 1 && iaxi2planeturb != 1) {
            for (j = 1; j <= jdim-1; j++) {
                // Interior points
                for (i = 2; i <= idim-2; i++) {
                    il = i-1; iu = i+1;
                    for (k = 1; k <= kdim-1; k++) {
                        dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                        dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                        sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                        sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                        sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                        sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                        ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                        ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                        anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                        anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                        fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                        fnum = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                        cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                        bz(k,i) = -cdm; cz(k,i) = cdp+cdm; dz(k,i) = -cdp;
                        cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                        bz2(k,i) = -cdm; cz2(k,i) = cdp+cdm; dz2(k,i) = -cdp;
                        cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                        bz3(k,i) = -cdm; cz3(k,i) = cdp+cdm; dz3(k,i) = -cdp;
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                        uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                        { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                        app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                        bz(k,i) = bz(k,i)-uu*app; bz2(k,i) = bz2(k,i)-uu*app; bz3(k,i) = bz3(k,i)-uu*app;
                        cz(k,i) = cz(k,i)+uu*(app-apm); cz2(k,i) = cz2(k,i)+uu*(app-apm); cz3(k,i) = cz3(k,i)+uu*(app-apm);
                        dz(k,i) = dz(k,i)+uu*apm; dz2(k,i) = dz2(k,i)+uu*apm; dz3(k,i) = dz3(k,i)+uu*apm;
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        double czadd = 0., cz2add = 0.;
                        if (ivmx == 30 || ivmx == 31) {
                            cmuc  = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                            betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                            { double a1t=damp1(j,k,i); czadd = 2.*re*betax*turre(j,k,i,1)+ccomplex_ns::ccabs(a1t)/turre(j,k,i,1); }
                            cz2add = re*cmuc*turre(j,k,i,1);
                        }
                        cz(k,i) = cz(k,i)+czadd; cz2(k,i) = cz2(k,i)+cz2add;
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        bz(k,i) = bz(k,i)*timestp(j,k,i);
                        bz2(k,i) = bz2(k,i)*timestp(j,k,i)*factor2;
                        bz3(k,i) = bz3(k,i)*timestp(j,k,i)*factor3;
                        cz(k,i) = cz(k,i)*timestp(j,k,i)+1.0*(1.+phi);
                        cz2(k,i) = cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                        cz3(k,i) = cz3(k,i)*timestp(j,k,i)*factor3+1.0*(1.+phi);
                        dz(k,i) = dz(k,i)*timestp(j,k,i);
                        dz2(k,i) = dz2(k,i)*timestp(j,k,i)*factor2;
                        dz3(k,i) = dz3(k,i)*timestp(j,k,i)*factor3;
                        fz(k,i) = rhside(j,k,i,1)*(1.+phi);
                        fz2(k,i) = rhside(j,k,i,2)*(1.+phi);
                        fz3(k,i) = rhside(j,k,i,3)*(1.+phi);
                    }
                } // end interior i loop

                // I0 boundary points
                i = 1; il = 1; iu = std::min(2, idim-1);
                for (k = 1; k <= kdim-1; k++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                    dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                    sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                    ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                    anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                    fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                    fnum = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                    cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    bz(k,i) = -cdm; cz(k,i) = cdp+cdm; dz(k,i) = -cdp;
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    bz2(k,i) = -cdm; cz2(k,i) = cdp+cdm; dz2(k,i) = -cdp;
                    cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                    bz3(k,i) = -cdm; cz3(k,i) = cdp+cdm; dz3(k,i) = -cdp;
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                    app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                    bz(k,i) = bz(k,i)-uu*app; bz2(k,i) = bz2(k,i)-uu*app; bz3(k,i) = bz3(k,i)-uu*app;
                    cz(k,i) = cz(k,i)+uu*(app-apm); cz2(k,i) = cz2(k,i)+uu*(app-apm); cz3(k,i) = cz3(k,i)+uu*(app-apm);
                    dz(k,i) = dz(k,i)+uu*apm; dz2(k,i) = dz2(k,i)+uu*apm; dz3(k,i) = dz3(k,i)+uu*apm;
                }
                for (k = 1; k <= kdim-1; k++) {
                    double czadd = 0., cz2add = 0.;
                    if (ivmx == 30 || ivmx == 31) {
                        cmuc  = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                        betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                        { double a1t=damp1(j,k,i); czadd = 2.*re*betax*turre(j,k,i,1)+ccomplex_ns::ccabs(a1t)/turre(j,k,i,1); }
                        cz2add = re*cmuc*turre(j,k,i,1);
                    }
                    cz(k,i) = cz(k,i)+czadd; cz2(k,i) = cz2(k,i)+cz2add;
                }
                for (k = 1; k <= kdim-1; k++) {
                    bz(k,i) = bz(k,i)*timestp(j,k,i);
                    bz2(k,i) = bz2(k,i)*timestp(j,k,i)*factor2;
                    bz3(k,i) = bz3(k,i)*timestp(j,k,i)*factor3;
                    cz(k,i) = cz(k,i)*timestp(j,k,i)+1.0*(1.+phi);
                    cz2(k,i) = cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                    cz3(k,i) = cz3(k,i)*timestp(j,k,i)*factor3+1.0*(1.+phi);
                    dz(k,i) = dz(k,i)*timestp(j,k,i);
                    dz2(k,i) = dz2(k,i)*timestp(j,k,i)*factor2;
                    dz3(k,i) = dz3(k,i)*timestp(j,k,i)*factor3;
                    fz(k,i) = rhside(j,k,i,1)*(1.+phi);
                    fz2(k,i) = rhside(j,k,i,2)*(1.+phi);
                    fz3(k,i) = rhside(j,k,i,3)*(1.+phi);
                }

                // IDIM boundary points
                i = idim-1; il = idim-2; iu = idim-1;
                for (k = 1; k <= kdim-1; k++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                    dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                    sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
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
                    ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                    anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                    fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                    fnum = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                    cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    bz(k,i) = -cdm; cz(k,i) = cdp+cdm; dz(k,i) = -cdp;
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    bz2(k,i) = -cdm; cz2(k,i) = cdp+cdm; dz2(k,i) = -cdp;
                    cdp = (sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                    bz3(k,i) = -cdm; cz3(k,i) = cdp+cdm; dz3(k,i) = -cdp;
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double a1t=1., a2t=uu; sgnu=ccomplex_ns::ccsignrc(a1t,a2t); }
                    app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                    bz(k,i) = bz(k,i)-uu*app; bz2(k,i) = bz2(k,i)-uu*app; bz3(k,i) = bz3(k,i)-uu*app;
                    cz(k,i) = cz(k,i)+uu*(app-apm); cz2(k,i) = cz2(k,i)+uu*(app-apm); cz3(k,i) = cz3(k,i)+uu*(app-apm);
                    dz(k,i) = dz(k,i)+uu*apm; dz2(k,i) = dz2(k,i)+uu*apm; dz3(k,i) = dz3(k,i)+uu*apm;
                }
                for (k = 1; k <= kdim-1; k++) {
                    double czadd = 0., cz2add = 0.;
                    if (ivmx == 30 || ivmx == 31) {
                        cmuc  = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                        betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                        { double a1t=damp1(j,k,i); czadd = 2.*re*betax*turre(j,k,i,1)+ccomplex_ns::ccabs(a1t)/turre(j,k,i,1); }
                        cz2add = re*cmuc*turre(j,k,i,1);
                    }
                    cz(k,i) = cz(k,i)+czadd; cz2(k,i) = cz2(k,i)+cz2add;
                }
                for (k = 1; k <= kdim-1; k++) {
                    bz(k,i) = bz(k,i)*timestp(j,k,i);
                    bz2(k,i) = bz2(k,i)*timestp(j,k,i)*factor2;
                    bz3(k,i) = bz3(k,i)*timestp(j,k,i)*factor3;
                    cz(k,i) = cz(k,i)*timestp(j,k,i)+1.0*(1.+phi);
                    cz2(k,i) = cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                    cz3(k,i) = cz3(k,i)*timestp(j,k,i)*factor3+1.0*(1.+phi);
                    dz(k,i) = dz(k,i)*timestp(j,k,i);
                    dz2(k,i) = dz2(k,i)*timestp(j,k,i)*factor2;
                    dz3(k,i) = dz3(k,i)*timestp(j,k,i)*factor3;
                    fz(k,i) = rhside(j,k,i,1)*(1.+phi);
                    fz2(k,i) = rhside(j,k,i,2)*(1.+phi);
                    fz3(k,i) = rhside(j,k,i,3)*(1.+phi);
                }

                // Apply blank (overset)
                if (iover == 1) {
                    for (i = 1; i <= idim-1; i++)
                    for (k = 1; k <= kdim-1; k++) {
                        fz(k,i) = fz(k,i)*blank(j,k,i);
                        bz(k,i) = bz(k,i)*blank(j,k,i);
                        dz(k,i) = dz(k,i)*blank(j,k,i);
                        cz(k,i) = cz(k,i)*blank(j,k,i)+(1.-blank(j,k,i));
                        fz2(k,i) = fz2(k,i)*blank(j,k,i);
                        bz2(k,i) = bz2(k,i)*blank(j,k,i);
                        dz2(k,i) = dz2(k,i)*blank(j,k,i);
                        cz2(k,i) = cz2(k,i)*blank(j,k,i)+(1.-blank(j,k,i));
                        fz3(k,i) = fz3(k,i)*blank(j,k,i);
                        bz3(k,i) = bz3(k,i)*blank(j,k,i);
                        dz3(k,i) = dz3(k,i)*blank(j,k,i);
                        cz3(k,i) = cz3(k,i)*blank(j,k,i)+(1.-blank(j,k,i));
                    }
                }

                // Solve tridiagonal systems in i direction
                { int kd1=kdim-1, id1=idim-1, one1=1, one2=1;
                  triv_ns::triv(kd1, id1, one1, kd1, one2, id1, workz, bz, cz, dz, fz);
                  triv_ns::triv(kd1, id1, one1, kd1, one2, id1, workz, bz2, cz2, dz2, fz2);
                  triv_ns::triv(kd1, id1, one1, kd1, one2, id1, workz, bz3, cz3, dz3, fz3); }

                // Store result back in rhside
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++) {
                    rhside(j,k,i,1) = fz(k,i);
                    rhside(j,k,i,2) = fz2(k,i);
                    rhside(j,k,i,3) = fz3(k,i);
                }
            } // end j loop for i sweep
        } // end if (i2d != 1 && iaxi2planeturb != 1) for implicit zeta

        // Update TURRE
        sumno = 0.; sumnk = 0.; sumni = 0.;
        negno = 0; negnk = 0; negni = 0;
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            sumno = sumno + rhside(j,k,i,1)*rhside(j,k,i,1);
            if ((float)(turre(j,k,i,1)+rhside(j,k,i,1)) <= (float)tur1cutlev) {
                negno++;
                if ((float)tur1cut > 0.) turre(j,k,i,1) = (double)tur1cut;
            } else {
                turre(j,k,i,1) = turre(j,k,i,1)+rhside(j,k,i,1);
            }
            sumnk = sumnk + rhside(j,k,i,2)*rhside(j,k,i,2);
            if ((float)(turre(j,k,i,2)+rhside(j,k,i,2)) <= (float)tur2cutlev) {
                negnk++;
                if ((float)tur2cut > 0.) turre(j,k,i,2) = (double)tur2cut;
            } else {
                turre(j,k,i,2) = turre(j,k,i,2)+rhside(j,k,i,2);
            }
            sumni = sumni + rhside(j,k,i,3)*rhside(j,k,i,3);
            if (itrans_on == 1) {
                if ((float)(turre(j,k,i,3)+rhside(j,k,i,3)) <= 0.) {
                    negni++;
                    turre(j,k,i,3) = 0.;
                } else if ((float)(turre(j,k,i,3)+rhside(j,k,i,3)) >= 100.) {
                    turre(j,k,i,3) = 100.;
                } else {
                    turre(j,k,i,3) = turre(j,k,i,3)+rhside(j,k,i,3);
                }
            }
        }
        sumno = std::sqrt(sumno)/(float)((kdim-1)*(jdim-1)*(idim-1));
        sumnk = std::sqrt(sumnk)/(float)((kdim-1)*(jdim-1)*(idim-1));
        sumni = std::sqrt(sumni)/(float)((kdim-1)*(jdim-1)*(idim-1));

    } // end do 500 not=1,nsubit

    sumn1 = sumno;
    sumn2 = sumnk;
    sumn3 = sumni;
    negn1 = negno;
    negn2 = negnk;
    negn3 = negni;

    // Update VIST3D and save omega and k values
    if (ivmx == 30 || ivmx == 31) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            { double a1t=smin(j,k,i); denom1 = 2./re*std::sqrt(turre(j,k,i,2))/(.09*turre(j,k,i,1)*ccomplex_ns::ccabs(a1t)); }
            denom2 = 500.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*re*re*
                smin(j,k,i)*turre(j,k,i,1));
            { double a1t=denom1, a2t=denom2; arg2=ccomplex_ns::ccmax(a1t,a2t); }
            { double a1t=arg2*arg2; f2=ccomplex_ns::cctanh(a1t); }
            if (isstdenom == 1) {
                s11 = ux(j,k,i,1); s22 = ux(j,k,i,5); s33 = ux(j,k,i,9);
                s12 = 0.5*(ux(j,k,i,2)+ux(j,k,i,4));
                s13 = 0.5*(ux(j,k,i,3)+ux(j,k,i,7));
                s23 = 0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                xis = s11*s11+s22*s22+s33*s33+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                xxx = std::sqrt(2.*xis);
                { double a1t=a1*turre(j,k,i,1), a2t=xxx*f2/re; denom=ccomplex_ns::ccmax(a1t,a2t); }
            } else {
                { double a1t=a1*turre(j,k,i,1), a2t=vor(j,k,i)*f2/re; denom=ccomplex_ns::ccmax(a1t,a2t); }
            }
            vist3d(j,k,i) = a1*q(j,k,i,1)*turre(j,k,i,2)/denom;
            { double a1t=vist3d(j,k,i), a2t=(double)edvislim; vist3d(j,k,i)=ccomplex_ns::ccmin(a1t,a2t); }
            zksav(j,k,i,1) = turre(j,k,i,1);
            zksav(j,k,i,2) = turre(j,k,i,2);
            zksav(j,k,i,3) = turre(j,k,i,3);
        }
    }

    // i_lam_forcezero: force vist3d to zero in laminar region
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


    // iwriteaux section (iwriteaux is hardwired to 0, but translate faithfully)
    iwriteaux = 0;
    if (iwriteaux == 1) {
        if (icyc == ncyc1[0] || icyc == ncyc1[1] || icyc == ncyc1[2] ||
            icyc == ncyc1[3] || icyc == ncyc1[4]) {
            if (ntime == 1) {
                jset = (int)(.779*(float)jdim);
                fortran_write_unit(92, "%5d jset,x=%5d%12.5e  u+,k+,e+,-uv+,logy+\n",
                    kdim-1, jset, (float)(0.5*(x(jset,1,1)+x(jset+1,1,1))));
                { double a1t=smin(jset,1,1); qset = std::sqrt(
                    (q(jset,1,1,2)-qk0(jset,1,2,1))*(q(jset,1,1,2)-qk0(jset,1,2,1))+
                    (q(jset,1,1,3)-qk0(jset,1,3,1))*(q(jset,1,1,3)-qk0(jset,1,3,1))+
                    (q(jset,1,1,4)-qk0(jset,1,4,1))*(q(jset,1,1,4)-qk0(jset,1,4,1)));
                  utau = std::sqrt((fnu(jset,1,1)+vk0(jset,1,1,1))*qset/
                    (ccomplex_ns::ccabs(a1t)*q(jset,1,1,1)*re)); }
                for (k = 1; k <= kdim-1; k++) {
                    { double a1t=smin(jset,k,1); ypl = re*q(jset,1,1,1)*utau*ccomplex_ns::ccabs(a1t)/fnu(jset,1,1); }
                    { double a1t=ypl; ypl = ccomplex_ns::cclog10(a1t); }
                    uplus = std::sqrt(q(jset,k,1,2)*q(jset,k,1,2)+q(jset,k,1,3)*q(jset,k,1,3)+q(jset,k,1,4)*q(jset,k,1,4))/utau;
                    zkplus = turre(jset,k,1,2)/(utau*utau);
                    { double a1t=smin(jset,1,1); eplus = turre(jset,k,1,1)*re/(utau*utau*qset/ccomplex_ns::ccabs(a1t)); }
                    uvplus = vist3d(jset,k,1)*vor(jset,k,1)/(re*utau*utau*q(jset,k,1,1));
                    fortran_write_unit(92, "%15.5e%15.5e%15.5e%15.5e%15.5e\n",
                        (float)uplus, (float)zkplus, (float)eplus, (float)uvplus, (float)ypl);
                }
                nnumb = 2;
                fortran_write_unit(91, "%5d\n", nnumb);
                fortran_write_unit(91, "   uv\n");
                fortran_write_unit(91, "   y\n");
                int jsets[] = {68, 93, 107, 113, 123, 129, 134};
                for (int jj = 0; jj < 7; jj++) {
                    jset = jsets[jj];
                    fortran_write_unit(91, "%5d jset=%5d  -uv/uinf**2,y  x=%15.5e\n",
                        kdim-1, jset, (float)(0.5*(x(jset,1,1)+x(jset+1,1,1))));
                    for (k = 1; k <= kdim-1; k++) {
                        { double a1t=smin(jset,k,1);
                          fortran_write_unit(91, "%15.5e%15.5e\n",
                            (float)(vist3d(jset,k,1)*vor(jset,k,1)/q(jset,k,1,1)/((double)reue*(double)xmach)),
                            (float)ccomplex_ns::ccabs(a1t)); }
                    }
                }
            }
        }
    }

} // end threeeqn

} // namespace threeeqn_ns
