// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "foureqn.h"
#include "ccomplex.h"
#include "triv.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace ccomplex_ns;
using namespace triv_ns;

namespace foureqn_ns {

void foureqn(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, FortranArray3DRef<double> dtj, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vor, FortranArray3DRef<double> smin, FortranArray4DRef<double> zksav, FortranArray4DRef<double> turre, FortranArray3DRef<double> damp1, FortranArray3DRef<double> blend, FortranArray3DRef<double> timestp, FortranArray3DRef<double> wrks, FortranArray3DRef<double> fnu, FortranArray2DRef<double> bx, FortranArray2DRef<double> bx2, FortranArray2DRef<double> cx, FortranArray2DRef<double> cx2, FortranArray2DRef<double> dx, FortranArray2DRef<double> dx2, FortranArray2DRef<double> fx, FortranArray2DRef<double> fx2, FortranArray2DRef<double> workx, FortranArray2DRef<double> by, FortranArray2DRef<double> by2, FortranArray2DRef<double> cy, FortranArray2DRef<double> cy2, FortranArray2DRef<double> dy, FortranArray2DRef<double> dy2, FortranArray2DRef<double> fy, FortranArray2DRef<double> fy2, FortranArray2DRef<double> worky, FortranArray2DRef<double> bz, FortranArray2DRef<double> bz2, FortranArray2DRef<double> cz, FortranArray2DRef<double> cz2, FortranArray2DRef<double> dz, FortranArray2DRef<double> dz2, FortranArray2DRef<double> fz, FortranArray2DRef<double> fz2, FortranArray2DRef<double> workz, int& ntime, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, int& nbl, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, FortranArray3DRef<double> blank, int& iover, double& sumn1, double& sumn2, double& sumn3, double& sumn4, int& negn1, int& negn2, int& negn3, int& negn4, FortranArray4DRef<double> ux, FortranArray4DRef<double> rhside, FortranArray4DRef<double> zksav2, FortranArray3DRef<double> v3dtmp, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& maxbl, int& maxseg, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iex, int& iex2, int& iex3, FortranArray2DRef<double> bx3, FortranArray2DRef<double> cx3, FortranArray2DRef<double> dx3, FortranArray2DRef<double> fx3, FortranArray2DRef<double> by3, FortranArray2DRef<double> cy3, FortranArray2DRef<double> dy3, FortranArray2DRef<double> fy3, FortranArray2DRef<double> bz3, FortranArray2DRef<double> cz3, FortranArray2DRef<double> dz3, FortranArray2DRef<double> fz3, FortranArray2DRef<double> bx4, FortranArray2DRef<double> cx4, FortranArray2DRef<double> dx4, FortranArray2DRef<double> fx4, FortranArray2DRef<double> by4, FortranArray2DRef<double> cy4, FortranArray2DRef<double> dy4, FortranArray2DRef<double> fy4, FortranArray2DRef<double> bz4, FortranArray2DRef<double> cz4, FortranArray2DRef<double> dz4, FortranArray2DRef<double> fz4, FortranArray3DRef<double> xlscale, FortranArray3DRef<double> fdsav, int& nummem)
{
    // COMMON block aliases
    float& cdes        = cmn_des.cdes;
    int32_t& ides      = cmn_des.ides;
    float& cddes       = cmn_des.cddes;
    float& xmach       = cmn_info.xmach;
    float& dt          = cmn_info.dt;
    int32_t& icyc      = cmn_mgrd.icyc;
    int32_t& ncyc      = cmn_mgrd.ncyc;
    int32_t& level     = cmn_mgrd.level;
    int32_t& lglobal   = cmn_mgrd.lglobal;
    float& gamma       = cmn_fluid.gamma;
    float& cbar        = cmn_fluid2.cbar;
    float& reue        = cmn_reyue.reue;
    float& tinf        = cmn_reyue.tinf;
    int32_t& isklton   = cmn_sklton.isklton;
    int32_t& i2d       = cmn_twod.i2d;
    int32_t& iexp      = cmn_zero.iexp;
    int32_t& ivmx      = cmn_maxiv.ivmx;
    int32_t& ita       = cmn_unst.ita;
    int32_t& itaturb   = cmn_turbconv.itaturb;
    int32_t& nsubturb  = cmn_turbconv.nsubturb;
    int32_t& nfreeze   = cmn_turbconv.nfreeze;
    int32_t& iturbord  = cmn_turbconv.iturbord;
    float& tur1cut     = cmn_turbconv.tur1cut;
    float& tur2cut     = cmn_turbconv.tur2cut;
    float& tur1cutlev  = cmn_turbconv.tur1cutlev;
    float& tur2cutlev  = cmn_turbconv.tur2cutlev;
    float& edvislim    = cmn_turbconv.edvislim;
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
    float& prod2d3dtrace = cmn_konew.prod2d3dtrace;
    int32_t& isstrc    = cmn_curvat.isstrc;
    float& sstrc_crc   = cmn_curvat.sstrc_crc;
    int32_t& iaxi2planeturb = cmn_axisym.iaxi2planeturb;
    int32_t& istrongturbdis = cmn_axisym.istrongturbdis;
    float& rho0        = cmn_ivals.rho0;
    // tur10 array (0-based C array, 1-based Fortran)
    // cmn_ivals.tur10[0] = tur10(1), etc.
    // ncyc1 array
    // cmn_info.ncyc1[0] = ncyc1(1), etc.

    // Local variables
    int itrans_on;
    double tu_percent, ca1, ca2, ce1, ce2, sigma_f, s1, cthetat, sigma_thetat, ce2inv;
    double xminn;
    int nsubit;
    double factor1, factor2, factor3, factor4;
    double vk, a1;
    double cmuc1, beta1, sigo1, sigk1, alp1;
    double sigkmu, sigg;
    double cmuc2, beta2, sigo2, sigk2, alp2;
    int jd2;
    double re, c2b, c2bp;
    int iwrite;
    int not_; // loop variable (not is reserved in some contexts)
    int i, j, k, kl, ku, jl, ju, il, iu;
    double dfacep, dfacem, sigkp, sigkm, sigop, sigom;
    double volku, volkl, volju, voljl, voliu, volil;
    double xp, yp, zp, xm, ym, zm, xa, ya, za;
    double ttpo, ttmo, ttpn, ttmn, ttp, ttm;
    double anutp, anutm, fnup, fnum;
    double cdp, cdm;
    double byy, cyy, dyy, bxx, cxx, dxx, bzz, czz, dzz;
    double xc, yc, zc, tc, uu, sgnu, app, apm;
    double cyadd, cy2add, cy3add, cy4add;
    double phi;
    double tt;
    double sumno, sumnk, sumni, sumnr;
    int negno, negnk, negni, negnr;
    // source term variables
    double cutoff, betax, cmuc, alp;
    double s11, s22, s33, s12, s13, s23, tracepart, s11t, s22t, s33t, xis;
    double dist2, dist, fnuinv, re_v, re_t, re_omega;
    double f_turb, f_sublayer, f_length, re_thetac;
    double f_onset1, f_onset2, f_onset3, f_onset;
    double f_reattach;
    double uuu, vortpluseps, dist_delta;
    double f_thetat1, f_thetat2, f_thetat;
    double gamma_sep, gamma_eff;
    double du_dx, du_dy, du_dz, du_ds;
    double tu_percent_ltd, rey, reth_part, thetat;
    double xlam, dlam, eff_lambda, eff_deriv, resr, dresr, dthetat;
    double re_thetat;
    double f4, sij, ri;
    double pk, dk, dkfactor;
    double p_gamma, d_gamma;
    double p_re, d_re;
    double deltaj, deltak, deltai, delta, ell, rd, fd, term;
    double denom1, denom2, arg2, f2, denom, xxx;
    double re_y, fff3;
    double arga, argb, arg, temp;
    int nmo;
    // iwriteaux section
    int iwriteaux;
    int jset, nnumb;
    double qset, utau, ypl, uplus, zkplus, eplus, uvplus;
    // nss
    int nss;
    // damp1 CD term
    double ca_cd;
    // arg1 for blend
    double arg1;


    // Parameter: itrans_on=1
    itrans_on = 1;
    // Set transition-related parameters
    tu_percent = 100.0*std::sqrt(2.0*(double)cmn_ivals.tur10[1]/(3.0*(double)xmach*(double)xmach));
    ca1 = 2.0;
    ca2 = 0.06;
    ce1 = 1.0;
    ce2 = 50.0;
    sigma_f = 1.0;
    s1 = 2.0;
    cthetat = 0.03;
    sigma_thetat = 2.0;
    ce2inv = 1.0/ce2;

    if (isklton > 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     Computing turbulent viscosity using 4-eqns, block=%5d", nbl);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     nummem=%5d", nummem);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     WARNING: ivisc=40 still under development... use at your own risk!");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     WARNING: Results may be sensitive to I.C.s and how case is run!!!!");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     Freestream tur10-tur40 = %19.8e%19.8e%19.8e%19.8e",
            (double)cmn_ivals.tur10[0], (double)cmn_ivals.tur10[1],
            (double)cmn_ivals.tur10[2], (double)cmn_ivals.tur10[3]);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     cflturb(1)-(4) = %12.5f%12.5f%12.5f%12.5f",
            (double)cmn_turbconv.cflturb[0], (double)cmn_turbconv.cflturb[1],
            (double)cmn_turbconv.cflturb[2], (double)cmn_turbconv.cflturb[3]);
        if (iturbord == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     1st order advection on RHS");
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     2nd order advection on RHS");
        }
    }
    if (isklton > 0) {
        if (ivmx == 40) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-omega SST (Menter) + transition model (AIAA J 47(12) 2009, p2894-2906)");
        }
        if (itrans_on == 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     transition part turned OFF");
        }
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     Tu (percent)=%10.4f", tu_percent);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     mut_inf/muref=%10.4f",
            (double)rho0*(double)cmn_ivals.tur10[1]/(double)cmn_ivals.tur10[0]);
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
        if (ides == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   using model in conjunction with DES, cdes=%7.3f", (double)cdes);
        } else if (ides == 2) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   using model in conjunction with DDES, cdes=%7.3f", (double)cdes);
        } else if (ides == 3) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   using model in conjunction with MDDES, cdes=%7.3f, cddes=%7.3f", (double)cdes, (double)cddes);
        }
        if (std::abs((double)prod2d3dtrace - 0.5) < 0.01) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     Sij used in 2SijSij prod term forced to be traceless in 2-D sense");
        } else if (std::abs((double)prod2d3dtrace - 0.33333333) < 0.01) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     Sij used in 2SijSij prod term forced to be traceless in 3-D sense");
        }
        if (isstrc == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     SSTRC-type curvature correction employed (AIAA 98-2554), sstrc_crc=%5.2f", (double)sstrc_crc);
        }
        if (keepambient == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     4-eqn ambient turbulence levels not allowed to decay");
        }
        if (iaxi2planeturb == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     4-eqn model ignoring i-dir");
        }
        if (istrongturbdis == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     strong conserv - diss terms");
        }
    }


    // Note: (10.**(-iexp) is machine zero)
    xminn = std::pow(10.0, (double)(-iexp+1));

    // Set number of subiterations
    nsubit = nsubturb;

    // Set CFL factors
    factor1 = 10.0;
    factor2 = 10.0;
    factor3 = 0.1;
    factor4 = 0.1;

    if ((float)cmn_turbconv.cflturb[0] != 0.0f) factor1 = (double)cmn_turbconv.cflturb[0];
    if ((float)cmn_turbconv.cflturb[1] != 0.0f) factor2 = (double)cmn_turbconv.cflturb[1];
    if ((float)cmn_turbconv.cflturb[2] != 0.0f) factor3 = (double)cmn_turbconv.cflturb[2];
    if ((float)cmn_turbconv.cflturb[3] != 0.0f) factor4 = (double)cmn_turbconv.cflturb[3];
    factor2 = factor2/factor1;
    factor3 = factor3/factor1;
    factor4 = factor4/factor1;

    // Timestep for turb model
    if ((float)dt < 0.0f) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            double a1t = factor1*vol(j,k,i)/dtj(j,k,i), a2t = 100.0;
            timestp(j,k,i) = ccmincr(a1t, a2t);
        }
    } else {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            timestp(j,k,i) = (double)dt;
            factor2 = 1.0;
            factor3 = 1.0;
            factor4 = 1.0;
        }
    }

    // Set up constants
    vk = 0.41;
    a1 = 0.31;
    cmuc1 = 0.09;
    beta1 = 0.0; sigo1 = 0.0; sigk1 = 0.0; alp1 = 0.0;
    if (ivmx == 40) beta1 = 0.075;
    if (ivmx == 40) sigo1 = 0.5;
    if (ivmx == 40) sigk1 = 0.85;
    if (ivmx == 40) alp1 = beta1/cmuc1 - sigo1*vk*vk/std::sqrt(cmuc1);
    sigkmu = 1.0;
    sigg = 1.0/sigma_f;

    cmuc2 = 0.09;
    beta2 = 0.0828;
    sigo2 = 0.856;
    sigk2 = 1.00;
    alp2 = beta2/cmuc2 - sigo2*vk*vk/std::sqrt(cmuc2);

    jd2 = (jdim-1)/2;
    re = (double)reue/(double)xmach;
    c2b = (double)cbar/(double)tinf;
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
        sumn1 = 0.0; sumn2 = 0.0; sumn3 = 0.0; sumn4 = 0.0;
        negn1 = 0; negn2 = 0; negn3 = 0; negn4 = 0;
        return;
    }
    phi = 0.0;
    if ((float)dt > 0.0f) {
        if (std::abs(ita) == 2) {
            phi = 0.5;
        } else {
            phi = 0.0;
        }
        if (itaturb == 0) {
            phi = 0.0;
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
        v3dtmp(j,0,i) = bck(j,i,1)*((double)cmn_wallfun.iwf[2]*v3dtmp(j,1,i) +
            (1-(double)cmn_wallfun.iwf[2])*2.0*vk0(j,i,1,1)-v3dtmp(j,1,i)) +
            (1.0-bck(j,i,1))*vk0(j,i,1,1);
        v3dtmp(j,kdim,i) = bck(j,i,2)*((double)cmn_wallfun.iwf[2]*v3dtmp(j,kdim-1,i) +
            (1-(double)cmn_wallfun.iwf[2])*2.0*vk0(j,i,1,3)-v3dtmp(j,kdim-1,i)) +
            (1.0-bck(j,i,2))*vk0(j,i,1,3);
    }
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++) {
        v3dtmp(0,k,i) = bcj(k,i,1)*((double)cmn_wallfun.iwf[1]*v3dtmp(1,k,i) +
            (1-(double)cmn_wallfun.iwf[1])*2.0*vj0(k,i,1,1)-v3dtmp(1,k,i)) +
            (1.0-bcj(k,i,1))*vj0(k,i,1,1);
        v3dtmp(jdim,k,i) = bcj(k,i,2)*((double)cmn_wallfun.iwf[1]*v3dtmp(jdim-1,k,i) +
            (1-(double)cmn_wallfun.iwf[1])*2.0*vj0(k,i,1,3)-v3dtmp(jdim-1,k,i)) +
            (1.0-bcj(k,i,2))*vj0(k,i,1,3);
    }
    if (i2d != 1 && iaxi2planeturb != 1) {
        for (j = 1; j <= jdim-1; j++)
        for (k = 1; k <= kdim-1; k++) {
            v3dtmp(j,k,0) = bci(j,k,1)*((double)cmn_wallfun.iwf[0]*v3dtmp(j,k,1) +
                (1-(double)cmn_wallfun.iwf[0])*2.0*vi0(j,k,1,1)-v3dtmp(j,k,1)) +
                (1.0-bci(j,k,1))*vi0(j,k,1,1);
            v3dtmp(j,k,idim) = bci(j,k,2)*((double)cmn_wallfun.iwf[0]*v3dtmp(j,k,idim-1) +
                (1-(double)cmn_wallfun.iwf[0])*2.0*vi0(j,k,1,3)-v3dtmp(j,k,idim-1)) +
                (1.0-bci(j,k,2))*vi0(j,k,1,3);
        }
    }


    // If this is 1st global subiteration for time-accurate computation, save zksav
    if ((float)dt > 0.0f && icyc == 1) {
        if (std::abs(ita) == 2) {
            if ((float)zksav2(1,1,1,1) == 0.0f) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    zksav2(j,k,i,5) = 0.0;
                    zksav2(j,k,i,6) = 0.0;
                    zksav2(j,k,i,7) = 0.0;
                    zksav2(j,k,i,8) = 0.0;
                }
            } else {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    zksav2(j,k,i,5) = zksav(j,k,i,1) - zksav2(j,k,i,1);
                    zksav2(j,k,i,6) = zksav(j,k,i,2) - zksav2(j,k,i,2);
                    zksav2(j,k,i,7) = zksav(j,k,i,3) - zksav2(j,k,i,3);
                    zksav2(j,k,i,8) = zksav(j,k,i,4) - zksav2(j,k,i,4);
                }
            }
        }
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            zksav2(j,k,i,1) = zksav(j,k,i,1);
            zksav2(j,k,i,2) = zksav(j,k,i,2);
            zksav2(j,k,i,3) = zksav(j,k,i,3);
            zksav2(j,k,i,4) = zksav(j,k,i,4);
        }
    }

    // Get TURRE values
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++)
    for (j = 1; j <= jdim-1; j++) {
        turre(j,k,i,1) = zksav(j,k,i,1);
        turre(j,k,i,2) = zksav(j,k,i,2);
        turre(j,k,i,3) = zksav(j,k,i,3);
        turre(j,k,i,4) = zksav(j,k,i,4);
    }

    // Iterate to solve the equations
    for (not_ = 1; not_ <= nsubit; not_++) {

        // set up boundary conditions (they are in ghost cells everywhere)
        // (1) k=0 boundary:
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            turre(j,0,i,1) = tk0(j,i,1,1);
            turre(j,0,i,2) = tk0(j,i,2,1);
            turre(j,0,i,3) = tk0(j,i,3,1);
            turre(j,0,i,4) = tk0(j,i,4,1);
        }
        // (2) k=kdim boundary:
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            turre(j,kdim,i,1) = tk0(j,i,1,3);
            turre(j,kdim,i,2) = tk0(j,i,2,3);
            turre(j,kdim,i,3) = tk0(j,i,3,3);
            turre(j,kdim,i,4) = tk0(j,i,4,3);
        }
        // (3) j=0 boundary:
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            turre(0,k,i,1) = tj0(k,i,1,1);
            turre(0,k,i,2) = tj0(k,i,2,1);
            turre(0,k,i,3) = tj0(k,i,3,1);
            turre(0,k,i,4) = tj0(k,i,4,1);
        }
        // (4) j=jdim boundary:
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            turre(jdim,k,i,1) = tj0(k,i,1,3);
            turre(jdim,k,i,2) = tj0(k,i,2,3);
            turre(jdim,k,i,3) = tj0(k,i,3,3);
            turre(jdim,k,i,4) = tj0(k,i,4,3);
        }
        if (i2d != 1 && iaxi2planeturb != 1) {
            // (5) i=0 boundary:
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,k,0,1) = ti0(j,k,1,1);
                turre(j,k,0,2) = ti0(j,k,2,1);
                turre(j,k,0,3) = ti0(j,k,3,1);
                turre(j,k,0,4) = ti0(j,k,4,1);
            }
            // (6) i=idim boundary:
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,k,idim,1) = ti0(j,k,1,3);
                turre(j,k,idim,2) = ti0(j,k,2,3);
                turre(j,k,idim,3) = ti0(j,k,3,3);
                turre(j,k,idim,4) = ti0(j,k,4,3);
            }
        }


        if (iturbord != 1) {
            // (1) k=0 boundary (2nd ghost):
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,-1,i,1) = tk0(j,i,1,2);
                turre(j,-1,i,2) = tk0(j,i,2,2);
                turre(j,-1,i,3) = tk0(j,i,3,2);
                turre(j,-1,i,4) = tk0(j,i,4,2);
            }
            // (2) k=kdim boundary (2nd ghost):
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,kdim+1,i,1) = tk0(j,i,1,4);
                turre(j,kdim+1,i,2) = tk0(j,i,2,4);
                turre(j,kdim+1,i,3) = tk0(j,i,3,4);
                turre(j,kdim+1,i,4) = tk0(j,i,4,4);
            }
            // (3) j=0 boundary (2nd ghost):
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                turre(-1,k,i,1) = tj0(k,i,1,2);
                turre(-1,k,i,2) = tj0(k,i,2,2);
                turre(-1,k,i,3) = tj0(k,i,3,2);
                turre(-1,k,i,4) = tj0(k,i,4,2);
            }
            // (4) j=jdim boundary (2nd ghost):
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                turre(jdim+1,k,i,1) = tj0(k,i,1,4);
                turre(jdim+1,k,i,2) = tj0(k,i,2,4);
                turre(jdim+1,k,i,3) = tj0(k,i,3,4);
                turre(jdim+1,k,i,4) = tj0(k,i,4,4);
            }
            if (i2d != 1 && iaxi2planeturb != 1) {
                // (5) i=0 boundary (2nd ghost):
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    turre(j,k,-1,1) = ti0(j,k,1,2);
                    turre(j,k,-1,2) = ti0(j,k,2,2);
                    turre(j,k,-1,3) = ti0(j,k,3,2);
                    turre(j,k,-1,4) = ti0(j,k,4,2);
                }
                // (6) i=idim boundary (2nd ghost):
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    turre(j,k,idim+1,1) = ti0(j,k,1,4);
                    turre(j,k,idim+1,2) = ti0(j,k,2,4);
                    turre(j,k,idim+1,3) = ti0(j,k,3,4);
                    turre(j,k,idim+1,4) = ti0(j,k,4,4);
                }
            }
        }

        // Get damp1 = CD = cross derivative term for SST:
        if (ivmx == 40) {
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                double xa2 = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4) + sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                double ya2 = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4) + sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                double za2 = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4) + sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                double tt2 = xa2*xa2 + ya2*ya2 + za2*za2;
                double ca_cd2 = 2.0*sigo2*tt2/(turre(j,k,i,1)*re);
                damp1(j,k,i) = 0.25*ca_cd2*(turre(j,k+1,i,1)-turre(j,k-1,i,1))*
                                             (turre(j,k+1,i,2)-turre(j,k-1,i,2));
            }
            for (j = 1; j <= jdim-1; j++)
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                double xa2 = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4) + sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                double ya2 = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4) + sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                double za2 = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4) + sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                double tt2 = xa2*xa2 + ya2*ya2 + za2*za2;
                double ca_cd2 = 2.0*sigo2*tt2/(turre(j,k,i,1)*re);
                damp1(j,k,i) = damp1(j,k,i) +
                    0.25*ca_cd2*(turre(j+1,k,i,1)-turre(j-1,k,i,1))*
                                (turre(j+1,k,i,2)-turre(j-1,k,i,2));
            }
            if (i2d != 1 && iaxi2planeturb != 1) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    double xa2 = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4) + si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    double ya2 = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4) + si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    double za2 = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4) + si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    double tt2 = xa2*xa2 + ya2*ya2 + za2*za2;
                    double ca_cd2 = 2.0*sigo2*tt2/(turre(j,k,i,1)*re);
                    damp1(j,k,i) = damp1(j,k,i) +
                        0.25*ca_cd2*(turre(j,k,i+1,1)-turre(j,k,i-1,1))*
                                    (turre(j,k,i+1,2)-turre(j,k,i-1,2));
                }
            }
        }



        // get blend = F1 factor
        if (ivmx == 40) {
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                double smin_jki = smin(j,k,i);
                double smin_abs_v = smin_jki; double smin_abs = ccabs(smin_abs_v);
                double a1t = std::sqrt(turre(j,k,i,2))/(.09*re*turre(j,k,i,1)*smin_abs);
                double a2t = 500.0*fnu(j,k,i)/(q(j,k,i,1)*smin_jki*re*re*
                    smin_jki*turre(j,k,i,1));
                double arga2; { double av1=a1t, av2=a2t; arga2 = ccmax(av1,av2); }
                double temp2; { double av=damp1(j,k,i)*re, av2=1.e-20; temp2 = ccmaxcr(av,av2); }
                double argb2 = 4.0*sigo2*turre(j,k,i,2)/(temp2*smin_jki*smin_jki);
                double arg3; { double av1=arga2, av2=argb2; arg3 = ccmin(av1,av2); }
                double arg4 = arg3*arg3*arg3*arg3;
                { double av=arg4; blend(j,k,i) = cctanh(av); }
                if (itrans_on == 1) {
                    double re_y2 = q(j,k,i,1)*smin_abs*std::sqrt(turre(j,k,i,2))/fnu(j,k,i)*re;
                    double fff3_2 = std::exp(-std::pow(re_y2/120.0, 8.0));
                    double bv1 = blend(j,k,i), bv2 = fff3_2;
                    blend(j,k,i) = ccmax(bv1, bv2);
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
                sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
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
                double ttpo = xp*xa+yp*ya+zp*za;
                double ttmo = xm*xa+ym*ya+zm*za;
                double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                fnup = 0.5*(fnu(j,k+1,i)+fnu(j,k,i));
                fnum = 0.5*(fnu(j,k-1,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,1) = -byy*turre(j,k-1,i,1)-cyy*turre(j,k,i,1)-dyy*turre(j,k+1,i,1);
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,2) = -byy*turre(j,k-1,i,2)-cyy*turre(j,k,i,2)-dyy*turre(j,k+1,i,2);
                cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,3) = -byy*turre(j,k-1,i,3)-cyy*turre(j,k,i,3)-dyy*turre(j,k+1,i,3);
                cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,4) = -byy*turre(j,k-1,i,4)-cyy*turre(j,k,i,4)-dyy*turre(j,k+1,i,4);
            }
        }


        // K0 boundary points
        k = 1; kl = 1; ku = std::min(2, kdim-1);
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
            dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
            sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
            sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
            sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
            sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
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
            double ttpo = xp*xa+yp*ya+zp*za;
            double ttmo = xm*xa+ym*ya+zm*za;
            double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
            double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
            ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
            ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
            anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
            anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
            fnup = 0.5*(fnu(j,k+1,i)+fnu(j,k,i));
            fnum = 0.5*(fnu(j,k-1,i)+fnu(j,k,i));
            cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,1) = -byy*turre(j,k-1,i,1)-cyy*turre(j,k,i,1)-dyy*turre(j,k+1,i,1);
            cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,2) = -byy*turre(j,k-1,i,2)-cyy*turre(j,k,i,2)-dyy*turre(j,k+1,i,2);
            cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,3) = -byy*turre(j,k-1,i,3)-cyy*turre(j,k,i,3)-dyy*turre(j,k+1,i,3);
            cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
            cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,4) = -byy*turre(j,k-1,i,4)-cyy*turre(j,k,i,4)-dyy*turre(j,k+1,i,4);
        }


        // KDIM boundary points
        k = kdim-1; kl = kdim-2; ku = kdim-1;
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
            dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
            sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
            sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
            sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
            sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
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
            double ttpo = xp*xa+yp*ya+zp*za;
            double ttmo = xm*xa+ym*ya+zm*za;
            double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
            double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
            ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
            ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
            anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
            anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
            fnup = 0.5*(fnu(j,k+1,i)+fnu(j,k,i));
            fnum = 0.5*(fnu(j,k-1,i)+fnu(j,k,i));
            cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,1) = -byy*turre(j,k-1,i,1)-cyy*turre(j,k,i,1)-dyy*turre(j,k+1,i,1);
            cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,2) = -byy*turre(j,k-1,i,2)-cyy*turre(j,k,i,2)-dyy*turre(j,k+1,i,2);
            cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,3) = -byy*turre(j,k-1,i,3)-cyy*turre(j,k,i,3)-dyy*turre(j,k+1,i,3);
            cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
            cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
            byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
            rhside(j,k,i,4) = -byy*turre(j,k-1,i,4)-cyy*turre(j,k,i,4)-dyy*turre(j,k+1,i,4);
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
                { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                rhside(j,k,i,1) = rhside(j,k,i,1)-uu*(app*(turre(j,k,i,1)-turre(j,k-1,i,1))+apm*(turre(j,k+1,i,1)-turre(j,k,i,1)));
                rhside(j,k,i,2) = rhside(j,k,i,2)-uu*(app*(turre(j,k,i,2)-turre(j,k-1,i,2))+apm*(turre(j,k+1,i,2)-turre(j,k,i,2)));
                rhside(j,k,i,3) = rhside(j,k,i,3)-uu*(app*(turre(j,k,i,3)-turre(j,k-1,i,3))+apm*(turre(j,k+1,i,3)-turre(j,k,i,3)));
                rhside(j,k,i,4) = rhside(j,k,i,4)-uu*(app*(turre(j,k,i,4)-turre(j,k-1,i,4))+apm*(turre(j,k+1,i,4)-turre(j,k,i,4)));
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
                { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                rhside(j,k,i,1) = rhside(j,k,i,1)
                    -0.5*uu*app*turre(j,k-2,i,1)+2.0*uu*app*turre(j,k-1,i,1)
                    -1.5*uu*app*turre(j,k,i,1)+1.5*uu*apm*turre(j,k,i,1)
                    -2.0*uu*apm*turre(j,k+1,i,1)+0.5*uu*apm*turre(j,k+2,i,1);
                rhside(j,k,i,2) = rhside(j,k,i,2)
                    -0.5*uu*app*turre(j,k-2,i,2)+2.0*uu*app*turre(j,k-1,i,2)
                    -1.5*uu*app*turre(j,k,i,2)+1.5*uu*apm*turre(j,k,i,2)
                    -2.0*uu*apm*turre(j,k+1,i,2)+0.5*uu*apm*turre(j,k+2,i,2);
                rhside(j,k,i,3) = rhside(j,k,i,3)
                    -0.5*uu*app*turre(j,k-2,i,3)+2.0*uu*app*turre(j,k-1,i,3)
                    -1.5*uu*app*turre(j,k,i,3)+1.5*uu*apm*turre(j,k,i,3)
                    -2.0*uu*apm*turre(j,k+1,i,3)+0.5*uu*apm*turre(j,k+2,i,3);
                rhside(j,k,i,4) = rhside(j,k,i,4)
                    -0.5*uu*app*turre(j,k-2,i,4)+2.0*uu*app*turre(j,k-1,i,4)
                    -1.5*uu*app*turre(j,k,i,4)+1.5*uu*apm*turre(j,k,i,4)
                    -2.0*uu*apm*turre(j,k+1,i,4)+0.5*uu*apm*turre(j,k+2,i,4);
            }
        }


        // F_xi_xi viscous terms - Interior points
        for (j = 2; j <= jdim-2; j++) {
            jl = j-1; ju = j+1;
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
                dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
                sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
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
                double ttpo = xp*xa+yp*ya+zp*za;
                double ttmo = xm*xa+ym*ya+zm*za;
                double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                fnup = 0.5*(fnu(j+1,k,i)+fnu(j,k,i));
                fnum = 0.5*(fnu(j-1,k,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,1) = rhside(j,k,i,1)-bxx*turre(j-1,k,i,1)-cxx*turre(j,k,i,1)-dxx*turre(j+1,k,i,1);
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,2) = rhside(j,k,i,2)-bxx*turre(j-1,k,i,2)-cxx*turre(j,k,i,2)-dxx*turre(j+1,k,i,2);
                cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,3) = rhside(j,k,i,3)-bxx*turre(j-1,k,i,3)-cxx*turre(j,k,i,3)-dxx*turre(j+1,k,i,3);
                cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,4) = rhside(j,k,i,4)-bxx*turre(j-1,k,i,4)-cxx*turre(j,k,i,4)-dxx*turre(j+1,k,i,4);
            }
        }


        // J0 boundary points
        j = 1; jl = 1; ju = std::min(2, jdim-1);
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
            dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
            sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
            sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
            sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
            sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
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
            double ttpo = xp*xa+yp*ya+zp*za;
            double ttmo = xm*xa+ym*ya+zm*za;
            double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
            double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
            ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
            ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
            anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
            anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
            fnup = 0.5*(fnu(j+1,k,i)+fnu(j,k,i));
            fnum = 0.5*(fnu(j-1,k,i)+fnu(j,k,i));
            cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,1) = rhside(j,k,i,1)-bxx*turre(j-1,k,i,1)-cxx*turre(j,k,i,1)-dxx*turre(j+1,k,i,1);
            cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,2) = rhside(j,k,i,2)-bxx*turre(j-1,k,i,2)-cxx*turre(j,k,i,2)-dxx*turre(j+1,k,i,2);
            cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,3) = rhside(j,k,i,3)-bxx*turre(j-1,k,i,3)-cxx*turre(j,k,i,3)-dxx*turre(j+1,k,i,3);
            cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
            cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,4) = rhside(j,k,i,4)-bxx*turre(j-1,k,i,4)-cxx*turre(j,k,i,4)-dxx*turre(j+1,k,i,4);
        }


        // JDIM boundary points
        j = jdim-1; jl = jdim-2; ju = jdim-1;
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
            dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
            sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
            sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
            sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
            sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
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
            double ttpo = xp*xa+yp*ya+zp*za;
            double ttmo = xm*xa+ym*ya+zm*za;
            double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
            double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
            ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
            ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
            anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
            anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
            fnup = 0.5*(fnu(j+1,k,i)+fnu(j,k,i));
            fnum = 0.5*(fnu(j-1,k,i)+fnu(j,k,i));
            cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,1) = rhside(j,k,i,1)-bxx*turre(j-1,k,i,1)-cxx*turre(j,k,i,1)-dxx*turre(j+1,k,i,1);
            cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,2) = rhside(j,k,i,2)-bxx*turre(j-1,k,i,2)-cxx*turre(j,k,i,2)-dxx*turre(j+1,k,i,2);
            cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
            cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,3) = rhside(j,k,i,3)-bxx*turre(j-1,k,i,3)-cxx*turre(j,k,i,3)-dxx*turre(j+1,k,i,3);
            cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
            cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
            bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
            rhside(j,k,i,4) = rhside(j,k,i,4)-bxx*turre(j-1,k,i,4)-cxx*turre(j,k,i,4)-dxx*turre(j+1,k,i,4);
        }


        // Advective terms in xi
        if (iturbord == 1) {
            for (j = 1; j <= jdim-1; j++)
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                rhside(j,k,i,1) = rhside(j,k,i,1)-uu*(app*(turre(j,k,i,1)-turre(j-1,k,i,1))+apm*(turre(j+1,k,i,1)-turre(j,k,i,1)));
                rhside(j,k,i,2) = rhside(j,k,i,2)-uu*(app*(turre(j,k,i,2)-turre(j-1,k,i,2))+apm*(turre(j+1,k,i,2)-turre(j,k,i,2)));
                rhside(j,k,i,3) = rhside(j,k,i,3)-uu*(app*(turre(j,k,i,3)-turre(j-1,k,i,3))+apm*(turre(j+1,k,i,3)-turre(j,k,i,3)));
                rhside(j,k,i,4) = rhside(j,k,i,4)-uu*(app*(turre(j,k,i,4)-turre(j-1,k,i,4))+apm*(turre(j+1,k,i,4)-turre(j,k,i,4)));
            }
        } else {
            for (j = 1; j <= jdim-1; j++)
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                rhside(j,k,i,1) = rhside(j,k,i,1)
                    -0.5*uu*app*turre(j-2,k,i,1)+2.0*uu*app*turre(j-1,k,i,1)
                    -1.5*uu*app*turre(j,k,i,1)+1.5*uu*apm*turre(j,k,i,1)
                    -2.0*uu*apm*turre(j+1,k,i,1)+0.5*uu*apm*turre(j+2,k,i,1);
                rhside(j,k,i,2) = rhside(j,k,i,2)
                    -0.5*uu*app*turre(j-2,k,i,2)+2.0*uu*app*turre(j-1,k,i,2)
                    -1.5*uu*app*turre(j,k,i,2)+1.5*uu*apm*turre(j,k,i,2)
                    -2.0*uu*apm*turre(j+1,k,i,2)+0.5*uu*apm*turre(j+2,k,i,2);
                rhside(j,k,i,3) = rhside(j,k,i,3)
                    -0.5*uu*app*turre(j-2,k,i,3)+2.0*uu*app*turre(j-1,k,i,3)
                    -1.5*uu*app*turre(j,k,i,3)+1.5*uu*apm*turre(j,k,i,3)
                    -2.0*uu*apm*turre(j+1,k,i,3)+0.5*uu*apm*turre(j+2,k,i,3);
                rhside(j,k,i,4) = rhside(j,k,i,4)
                    -0.5*uu*app*turre(j-2,k,i,4)+2.0*uu*app*turre(j-1,k,i,4)
                    -1.5*uu*app*turre(j,k,i,4)+1.5*uu*apm*turre(j,k,i,4)
                    -2.0*uu*apm*turre(j+1,k,i,4)+0.5*uu*apm*turre(j+2,k,i,4);
            }
        }


        // F_zeta_zeta viscous terms - Interior points
        if (i2d != 1 && iaxi2planeturb != 1) {
            for (i = 2; i <= idim-2; i++) {
                il = i-1; iu = i+1;
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                    dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                    sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
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
                    double ttpo = xp*xa+yp*ya+zp*za;
                    double ttmo = xm*xa+ym*ya+zm*za;
                    double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                    anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                    fnup = 0.5*(fnu(j,k,i+1)+fnu(j,k,i));
                    fnum = 0.5*(fnu(j,k,i-1)+fnu(j,k,i));
                    cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,1) = rhside(j,k,i,1)-bzz*turre(j,k,i-1,1)-czz*turre(j,k,i,1)-dzz*turre(j,k,i+1,1);
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,2) = rhside(j,k,i,2)-bzz*turre(j,k,i-1,2)-czz*turre(j,k,i,2)-dzz*turre(j,k,i+1,2);
                    cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,3) = rhside(j,k,i,3)-bzz*turre(j,k,i-1,3)-czz*turre(j,k,i,3)-dzz*turre(j,k,i+1,3);
                    cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                    cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,4) = rhside(j,k,i,4)-bzz*turre(j,k,i-1,4)-czz*turre(j,k,i,4)-dzz*turre(j,k,i+1,4);
                }
            }
        }


        // I0 boundary points
        if (i2d != 1 && iaxi2planeturb != 1) {
            i = 1; il = 1; iu = std::min(2, idim-1);
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
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
                double ttpo = xp*xa+yp*ya+zp*za;
                double ttmo = xm*xa+ym*ya+zm*za;
                double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                fnup = 0.5*(fnu(j,k,i+1)+fnu(j,k,i));
                fnum = 0.5*(fnu(j,k,i-1)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,1) = rhside(j,k,i,1)-bzz*turre(j,k,i-1,1)-czz*turre(j,k,i,1)-dzz*turre(j,k,i+1,1);
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,2) = rhside(j,k,i,2)-bzz*turre(j,k,i-1,2)-czz*turre(j,k,i,2)-dzz*turre(j,k,i+1,2);
                cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,3) = rhside(j,k,i,3)-bzz*turre(j,k,i-1,3)-czz*turre(j,k,i,3)-dzz*turre(j,k,i+1,3);
                cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,4) = rhside(j,k,i,4)-bzz*turre(j,k,i-1,4)-czz*turre(j,k,i,4)-dzz*turre(j,k,i+1,4);
            }
        }


        // IDIM boundary points
        if (i2d != 1 && iaxi2planeturb != 1) {
            i = idim-1; il = idim-2; iu = idim-1;
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
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
                double ttpo = xp*xa+yp*ya+zp*za;
                double ttmo = xm*xa+ym*ya+zm*za;
                double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                fnup = 0.5*(fnu(j,k,i+1)+fnu(j,k,i));
                fnum = 0.5*(fnu(j,k,i-1)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,1) = rhside(j,k,i,1)-bzz*turre(j,k,i-1,1)-czz*turre(j,k,i,1)-dzz*turre(j,k,i+1,1);
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,2) = rhside(j,k,i,2)-bzz*turre(j,k,i-1,2)-czz*turre(j,k,i,2)-dzz*turre(j,k,i+1,2);
                cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,3) = rhside(j,k,i,3)-bzz*turre(j,k,i-1,3)-czz*turre(j,k,i,3)-dzz*turre(j,k,i+1,3);
                cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                rhside(j,k,i,4) = rhside(j,k,i,4)-bzz*turre(j,k,i-1,4)-czz*turre(j,k,i,4)-dzz*turre(j,k,i+1,4);
            }
        }


        // Advective terms in zeta
        if (i2d != 1 && iaxi2planeturb != 1) {
            if (iturbord == 1) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                    app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                    rhside(j,k,i,1) = rhside(j,k,i,1)-uu*(app*(turre(j,k,i,1)-turre(j,k,i-1,1))+apm*(turre(j,k,i+1,1)-turre(j,k,i,1)));
                    rhside(j,k,i,2) = rhside(j,k,i,2)-uu*(app*(turre(j,k,i,2)-turre(j,k,i-1,2))+apm*(turre(j,k,i+1,2)-turre(j,k,i,2)));
                    rhside(j,k,i,3) = rhside(j,k,i,3)-uu*(app*(turre(j,k,i,3)-turre(j,k,i-1,3))+apm*(turre(j,k,i+1,3)-turre(j,k,i,3)));
                    rhside(j,k,i,4) = rhside(j,k,i,4)-uu*(app*(turre(j,k,i,4)-turre(j,k,i-1,4))+apm*(turre(j,k,i+1,4)-turre(j,k,i,4)));
                }
            } else {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                    app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                    rhside(j,k,i,1) = rhside(j,k,i,1)
                        -0.5*uu*app*turre(j,k,i-2,1)+2.0*uu*app*turre(j,k,i-1,1)
                        -1.5*uu*app*turre(j,k,i,1)+1.5*uu*apm*turre(j,k,i,1)
                        -2.0*uu*apm*turre(j,k,i+1,1)+0.5*uu*apm*turre(j,k,i+2,1);
                    rhside(j,k,i,2) = rhside(j,k,i,2)
                        -0.5*uu*app*turre(j,k,i-2,2)+2.0*uu*app*turre(j,k,i-1,2)
                        -1.5*uu*app*turre(j,k,i,2)+1.5*uu*apm*turre(j,k,i,2)
                        -2.0*uu*apm*turre(j,k,i+1,2)+0.5*uu*apm*turre(j,k,i+2,2);
                    rhside(j,k,i,3) = rhside(j,k,i,3)
                        -0.5*uu*app*turre(j,k,i-2,3)+2.0*uu*app*turre(j,k,i-1,3)
                        -1.5*uu*app*turre(j,k,i,3)+1.5*uu*apm*turre(j,k,i,3)
                        -2.0*uu*apm*turre(j,k,i+1,3)+0.5*uu*apm*turre(j,k,i+2,3);
                    rhside(j,k,i,4) = rhside(j,k,i,4)
                        -0.5*uu*app*turre(j,k,i-2,4)+2.0*uu*app*turre(j,k,i-1,4)
                        -1.5*uu*app*turre(j,k,i,4)+1.5*uu*apm*turre(j,k,i,4)
                        -2.0*uu*apm*turre(j,k,i+1,4)+0.5*uu*apm*turre(j,k,i+2,4);
                }
            }
        }


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
        // DES implementation
        if (ides == 1 && ivmx == 40) {
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                double deltaj2 = 2.0*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4));
                double deltak2 = 2.0*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4));
                double deltai2 = 2.0*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4));
                double delta2; { double av1=deltaj2, av2=deltak2; delta2 = ccmax(av1,av2); }
                if (i2d != 1 && iaxi2planeturb != 1) {
                    double av1=delta2, av2=deltai2; delta2 = ccmax(av1,av2);
                }
                double ell2 = std::sqrt(turre(j,k,i,2))/(cmuc1*turre(j,k,i,1)*re);
                { double av1=ell2, av2=(double)cdes*delta2; xlscale(j,k,i) = ccmin(av1,av2); }
            }
        } else if (ides >= 2 && ivmx == 40) {
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                double deltaj2 = 2.0*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4));
                double deltak2 = 2.0*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4));
                double deltai2 = 2.0*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4));
                double delta2; { double av1=deltaj2, av2=deltak2; delta2 = ccmax(av1,av2); }
                if (i2d != 1 && iaxi2planeturb != 1) {
                    double av1=delta2, av2=deltai2; delta2 = ccmax(av1,av2);
                }
                double ell2 = std::sqrt(turre(j,k,i,2))/(cmuc1*turre(j,k,i,1)*re);
                { double av=smin(j,k,i); dist = ccabs(av); }
                double velterm = ux(j,k,i,1)*ux(j,k,i,1)+ux(j,k,i,2)*ux(j,k,i,2)+ux(j,k,i,3)*ux(j,k,i,3)+
                    ux(j,k,i,4)*ux(j,k,i,4)+ux(j,k,i,5)*ux(j,k,i,5)+ux(j,k,i,6)*ux(j,k,i,6)+
                    ux(j,k,i,7)*ux(j,k,i,7)+ux(j,k,i,8)*ux(j,k,i,8)+ux(j,k,i,9)*ux(j,k,i,9);
                rd = (vist3d(j,k,i)+fnu(j,k,i))/(q(j,k,i,1)*
                    std::sqrt(velterm)*vk*vk*dist*dist*re);
                { double av=8.0*rd; double av3=av*av*av; fd = 1.0-cctanh(av3); }
                { double av=ell2-(double)cdes*delta2; double zero=0.0; term = ccmaxrc(zero,av); }
                xlscale(j,k,i) = ell2 - fd*term;
                if (ides == 3) {
                    fdsav(j,k,i) = fd;
                }
            }
        }


        // Source terms for ivmx==40
        if (ivmx == 40) {
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                double cutoff2;
                if ((i >= ilamlo && i < ilamhi &&
                     j >= jlamlo && j < jlamhi &&
                     k >= klamlo && k < klamhi) ||
                    (float)smin(j,k,i) < 0.0f) {
                    cutoff2 = 0.0;
                } else if (ides == 3 && (float)fdsav(j,k,i) > (float)cddes) {
                    cutoff2 = (1.0 - fdsav(j,k,i))/(1.0-(double)cddes);
                } else {
                    cutoff2 = 1.0;
                }
                double betax2 = blend(j,k,i)*beta1+(1.0-blend(j,k,i))*beta2;
                double cmuc2v = blend(j,k,i)*cmuc1+(1.0-blend(j,k,i))*cmuc2;
                double alp2v  = blend(j,k,i)*alp1 +(1.0-blend(j,k,i))*alp2;
                // Sij values
                s11 = ux(j,k,i,1);
                s22 = ux(j,k,i,5);
                s33 = ux(j,k,i,9);
                s12 = 0.5*(ux(j,k,i,2)+ux(j,k,i,4));
                s13 = 0.5*(ux(j,k,i,3)+ux(j,k,i,7));
                s23 = 0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                tracepart = (s11+s22+s33)*(double)prod2d3dtrace;
                s11t = s11-tracepart; s22t = s22-tracepart; s33t = s33-tracepart;
                xis = s11t*s11t+s22t*s22t+s33t*s33t+2.0*s12*s12+2.0*s13*s13+2.0*s23*s23;
                dist2 = smin(j,k,i)*smin(j,k,i);
                dist = std::sqrt(dist2);
                fnuinv = 1.0/fnu(j,k,i);
                re_v = q(j,k,i,1)*dist2*std::sqrt(2.0*xis)*re*fnuinv;
                re_t = q(j,k,i,1)*turre(j,k,i,2)*fnuinv/turre(j,k,i,1);
                re_omega = q(j,k,i,1)*dist2*turre(j,k,i,1)*re*re*fnuinv;
                f_turb = std::exp(-std::pow(0.25*re_t, 4.0));
                f_sublayer = std::exp(-std::pow(0.005*re_omega, 2.0));
                if ((float)turre(j,k,i,4) < 400.0f) {
                    f_length = 39.8189-(119.27e-4)*turre(j,k,i,4)-(132.567e-6)*turre(j,k,i,4)*turre(j,k,i,4);
                } else if ((float)turre(j,k,i,4) >= 400.0f && (float)turre(j,k,i,4) < 596.0f) {
                    f_length = 263.404-(123.939e-2)*turre(j,k,i,4)+(194.548e-5)*turre(j,k,i,4)*turre(j,k,i,4)-(101.695e-8)*turre(j,k,i,4)*turre(j,k,i,4)*turre(j,k,i,4);
                } else if ((float)turre(j,k,i,4) >= 596.0f && (float)turre(j,k,i,4) < 1200.0f) {
                    f_length = 0.5-((turre(j,k,i,4)-596.0)*3.0e-4);
                } else {
                    f_length = 0.3188;
                }
                f_length = f_length*(1.0-f_sublayer)+40.0*f_sublayer;
                if ((float)turre(j,k,i,4) <= 1870.0f) {
                    re_thetac = turre(j,k,i,4)-(396.035e-2-(120.656e-4)*turre(j,k,i,4)+(868.23e-6)*turre(j,k,i,4)*turre(j,k,i,4)-(696.506e-9)*turre(j,k,i,4)*turre(j,k,i,4)*turre(j,k,i,4)+(174.105e-12)*turre(j,k,i,4)*turre(j,k,i,4)*turre(j,k,i,4)*turre(j,k,i,4));
                } else {
                    re_thetac = turre(j,k,i,4)-(593.11+(turre(j,k,i,4)-1870.0)*0.482);
                }
                f_onset1 = re_v/(2.193*re_thetac);
                { double av1=f_onset1, av2=f_onset1*f_onset1*f_onset1*f_onset1; f_onset2 = ccmax(av1,av2); }
                { double av=f_onset2; double two=2.0; f_onset2 = ccmincr(av, two); }
                { double av=1.0-std::pow(0.4*re_t,3.0); double zero=0.0; f_onset3 = ccmaxcr(av, zero); }
                { double av=f_onset2-f_onset3; double zero=0.0; f_onset = ccmaxcr(av, zero); }
                f_reattach = std::exp(-std::pow(0.05*re_t, 4.0));
                uuu = std::sqrt(q(j,k,i,2)*q(j,k,i,2)+q(j,k,i,3)*q(j,k,i,3)+q(j,k,i,4)*q(j,k,i,4));
                { double av=uuu; double eps=1.0e-20; uuu = ccmaxcr(av, eps); }
                vortpluseps = vor(j,k,i)+1.0e-20;
                dist_delta = q(j,k,i,1)*re*uuu*uuu/(375.0*vortpluseps*turre(j,k,i,4)*fnu(j,k,i));
                f_thetat1 = std::exp(-std::pow(1.0e-5*re_omega, 2.0))*std::exp(-std::pow(dist_delta, 4.0));
                f_thetat2 = 1.0-std::pow((turre(j,k,i,3)-ce2inv)/(1.0-ce2inv), 2.0);
                { double av1=f_thetat1, av2=f_thetat2; f_thetat = ccmax(av1,av2); }
                { double av=f_thetat; double one=1.0; f_thetat = ccmincr(av, one); }
                { double av=re_v/(3.235*re_thetac)-1.0; double zero=0.0; gamma_sep = ccmaxcr(av, zero); }
                { double av=gamma_sep*f_reattach*s1; double two=2.0; gamma_sep = ccmincr(av, two); }
                gamma_sep = gamma_sep*f_thetat;
                { double av1=turre(j,k,i,3), av2=gamma_sep; gamma_eff = ccmax(av1,av2); }
                du_dx = (q(j,k,i,2)*ux(j,k,i,1)+q(j,k,i,3)*ux(j,k,i,4)+q(j,k,i,4)*ux(j,k,i,7))/uuu;
                du_dy = (q(j,k,i,2)*ux(j,k,i,2)+q(j,k,i,3)*ux(j,k,i,5)+q(j,k,i,4)*ux(j,k,i,8))/uuu;
                du_dz = (q(j,k,i,2)*ux(j,k,i,3)+q(j,k,i,3)*ux(j,k,i,6)+q(j,k,i,4)*ux(j,k,i,9))/uuu;
                du_ds = q(j,k,i,2)/uuu*du_dx+q(j,k,i,3)/uuu*du_dy+q(j,k,i,4)/uuu*du_dz;
                tu_percent_ltd = 100.0*std::sqrt(2.0*turre(j,k,i,2)/3.0)/uuu;
                { double av=tu_percent_ltd; double mintu=0.027; tu_percent_ltd = ccmaxcr(av, mintu); }
                rey = q(j,k,i,1)*uuu*re/fnu(j,k,i);
                if ((float)tu_percent_ltd > 1.3f) {
                    reth_part = 331.5*std::pow(tu_percent_ltd-0.5658, -0.671);
                } else {
                    reth_part = 1173.51-589.428*tu_percent_ltd+0.2196/(tu_percent_ltd*tu_percent_ltd);
                }
                thetat = reth_part/rey;
                for (nmo = 1; nmo <= 10; nmo++) {
                    xlam = q(j,k,i,1)*thetat*thetat*fnuinv*du_ds*re;
                    dlam = 2.0*q(j,k,i,1)*thetat*fnuinv*du_ds*re;
                    if ((float)xlam <= -0.1f) { xlam = -0.1; dlam = 0.0; }
                    else if ((float)xlam >= 0.1f) { xlam = 0.1; dlam = 0.0; }
                    if ((float)xlam <= 0.0f) {
                        eff_lambda = 1.0-((-12.986*xlam-123.66*xlam*xlam-405.689*xlam*xlam*xlam)*std::exp(-std::pow(tu_percent_ltd/1.5, 1.5)));
                        eff_deriv  = (12.986+2.0*123.66*xlam+3.0*405.689*xlam*xlam)*dlam*std::exp(-std::pow(tu_percent_ltd/1.5, 1.5));
                    } else {
                        eff_lambda = 1.0+(0.275*(1.0-std::exp(-35.0*xlam))*std::exp(-2.0*tu_percent_ltd));
                        eff_deriv  = 35.0*0.275*std::exp(-35.0*xlam)*dlam*std::exp(-2.0*tu_percent_ltd);
                    }
                    resr = reth_part*eff_lambda-rey*thetat;
                    dresr = reth_part*eff_deriv-rey;
                    dthetat = -resr/dresr;
                    { double av1=dthetat, av2=-0.9*thetat; dthetat = ccmax(av1,av2); }
                    { double av1=dthetat, av2=0.9*thetat; dthetat = ccmin(av1,av2); }
                    thetat = thetat+dthetat;
                }
                re_thetat = rey*thetat;
                { double av=re_thetat; double twenty=20.0; re_thetat = ccmaxcr(av, twenty); }
                f4 = 1.0;
                if (isstrc == 1) {
                    sij = std::sqrt(2.0*xis)+1.0e-20;
                    ri = (vor(j,k,i)/sij)*(vor(j,k,i)/sij-1.0);
                    f4 = 1.0/(1.0+(double)sstrc_crc*ri);
                }
                if (ikoprod == 1) {
                    rhside(j,k,i,1) = rhside(j,k,i,1)+cutoff2*alp2v/re*2.0*xis-f4*re*betax2*turre(j,k,i,1)*turre(j,k,i,1)+(1.0-blend(j,k,i))*damp1(j,k,i);
                    pk = vist3d(j,k,i)/(q(j,k,i,1)*re)*2.0*xis;
                } else {
                    rhside(j,k,i,1) = rhside(j,k,i,1)+cutoff2*alp2v/re*vor(j,k,i)*vor(j,k,i)-f4*re*betax2*turre(j,k,i,1)*turre(j,k,i,1)+(1.0-blend(j,k,i))*damp1(j,k,i);
                    pk = vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)*vor(j,k,i);
                }
                if (ides != 0) {
                    dk = std::pow(turre(j,k,i,2), 1.5)/xlscale(j,k,i);
                } else {
                    dk = re*cmuc2v*turre(j,k,i,1)*turre(j,k,i,2);
                }
                if (itrans_on == 1) {
                    { double av=gamma_eff; double pt1=0.1; dkfactor = ccmaxcr(av, pt1); }
                    { double av=dkfactor; double one=1.0; dkfactor = ccmincr(av, one); }
                    dk = dk*dkfactor;
                    pk = pk*gamma_eff;
                }
                { double av1=pk, av2=(double)pklimterm*dk; pk = ccmin(av1,av2); }
                rhside(j,k,i,2) = rhside(j,k,i,2)+cutoff2*pk-dk;
                rhside(j,k,i,1) = rhside(j,k,i,1)+(double)keepambient*re*betax2*(double)cmn_ivals.tur10[0]*(double)cmn_ivals.tur10[0];
                rhside(j,k,i,2) = rhside(j,k,i,2)+(double)keepambient*re*cmuc2v*(double)cmn_ivals.tur10[0]*(double)cmn_ivals.tur10[1];
                p_gamma = f_length*ca1*std::sqrt(2.0*xis)*std::sqrt(turre(j,k,i,3)*f_onset)*(1.0-(ce1*turre(j,k,i,3)));
                d_gamma = ca2*vor(j,k,i)*turre(j,k,i,3)*f_turb*(ce2*turre(j,k,i,3)-1.0);
                rhside(j,k,i,3) = rhside(j,k,i,3)+p_gamma-d_gamma;
                wrks(j,k,i) = cthetat*q(j,k,i,1)*uuu*uuu*0.002*fnuinv*(1.0-f_thetat)*re;
                p_re = wrks(j,k,i)*re_thetat;
                d_re = wrks(j,k,i)*turre(j,k,i,4);
                rhside(j,k,i,4) = rhside(j,k,i,4)+p_re-d_re;
            }
        }



        // Implicit F_eta_eta viscous terms. Do over all i's
        for (i = 1; i <= idim-1; i++) {
            for (k = 2; k <= kdim-2; k++) {
                kl = k-1; ku = k+1;
                for (j = 1; j <= jdim-1; j++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
                    dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
                    sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
                    volku = vol(j,ku,i); volkl = vol(j,kl,i);
                    xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    double ttpo2 = xp*xa+yp*ya+zp*za;
                    double ttmo2 = xm*xa+ym*ya+zm*za;
                    double ttpn2 = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                    double ttmn2 = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                    ttp = ttpo2*(1-istrongturbdis)+ttpn2*istrongturbdis;
                    ttm = ttmo2*(1-istrongturbdis)+ttmn2*istrongturbdis;
                    anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                    anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                    fnup = 0.5*(fnu(j,k+1,i)+fnu(j,k,i));
                    fnum = 0.5*(fnu(j,k-1,i)+fnu(j,k,i));
                    cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    by(j,k) = -cdm; cy(j,k) = cdp+cdm; dy(j,k) = -cdp;
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    by2(j,k) = -cdm; cy2(j,k) = cdp+cdm; dy2(j,k) = -cdp;
                    cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                    by3(j,k) = -cdm; cy3(j,k) = cdp+cdm; dy3(j,k) = -cdp;
                    cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                    cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                    by4(j,k) = -cdm; cy4(j,k) = cdp+cdm; dy4(j,k) = -cdp;
                }
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                    app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                    by(j,k) -= uu*app; by2(j,k) -= uu*app; by3(j,k) -= uu*app; by4(j,k) -= uu*app;
                    cy(j,k) += uu*(app-apm); cy2(j,k) += uu*(app-apm); cy3(j,k) += uu*(app-apm); cy4(j,k) += uu*(app-apm);
                    dy(j,k) += uu*apm; dy2(j,k) += uu*apm; dy3(j,k) += uu*apm; dy4(j,k) += uu*apm;
                }
                for (j = 1; j <= jdim-1; j++) {
                    double cyadd2=0.0, cy2add2=0.0, cy3add2=0.0, cy4add2=0.0;
                    if (ivmx == 40) {
                        double cmuc_l = blend(j,k,i)*cmuc1+(1.0-blend(j,k,i))*cmuc2;
                        double betax_l = blend(j,k,i)*beta1+(1.0-blend(j,k,i))*beta2;
                        double damp1_abs; { double av=damp1(j,k,i); damp1_abs = ccabs(av); }
                        cyadd2 = 2.0*re*betax_l*turre(j,k,i,1)+damp1_abs/turre(j,k,i,1);
                        cy2add2 = re*cmuc_l*turre(j,k,i,1);
                        double vor_abs; { double av=2.0*ce2*turre(j,k,i,3)-1.0; vor_abs = ccabs(av); }
                        cy3add2 = ca2*vor(j,k,i)*vor_abs;
                        cy4add2 = wrks(j,k,i);
                    }
                    cy(j,k) += cyadd2; cy2(j,k) += cy2add2; cy3(j,k) += cy3add2; cy4(j,k) += cy4add2;
                }
                for (j = 1; j <= jdim-1; j++) {
                    by(j,k) *= timestp(j,k,i); by2(j,k) *= timestp(j,k,i)*factor2;
                    by3(j,k) *= timestp(j,k,i)*factor3; by4(j,k) *= timestp(j,k,i)*factor4;
                    cy(j,k) = cy(j,k)*timestp(j,k,i)+1.0*(1.0+phi);
                    cy2(j,k) = cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.0+phi);
                    cy3(j,k) = cy3(j,k)*timestp(j,k,i)*factor3+1.0*(1.0+phi);
                    cy4(j,k) = cy4(j,k)*timestp(j,k,i)*factor4+1.0*(1.0+phi);
                    dy(j,k) *= timestp(j,k,i); dy2(j,k) *= timestp(j,k,i)*factor2;
                    dy3(j,k) *= timestp(j,k,i)*factor3; dy4(j,k) *= timestp(j,k,i)*factor4;
                    fy(j,k) = rhside(j,k,i,1)*timestp(j,k,i);
                    fy2(j,k) = rhside(j,k,i,2)*timestp(j,k,i)*factor2;
                    fy3(j,k) = rhside(j,k,i,3)*timestp(j,k,i)*factor3;
                    fy4(j,k) = rhside(j,k,i,4)*timestp(j,k,i)*factor4;
                }
                if ((float)dt > 0.0f) {
                    for (j = 1; j <= jdim-1; j++) {
                        fy(j,k) += (1.0+phi)*(zksav2(j,k,i,1)-turre(j,k,i,1))+phi*zksav2(j,k,i,5);
                        fy2(j,k) += (1.0+phi)*(zksav2(j,k,i,2)-turre(j,k,i,2))+phi*zksav2(j,k,i,6);
                        fy3(j,k) += (1.0+phi)*(zksav2(j,k,i,3)-turre(j,k,i,3))+phi*zksav2(j,k,i,7);
                        fy4(j,k) += (1.0+phi)*(zksav2(j,k,i,4)-turre(j,k,i,4))+phi*zksav2(j,k,i,8);
                    }
                }
            }


            // K0 boundary points
            k = 1; kl = 1; ku = std::min(2, kdim-1);
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
                dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
                sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
                volku = vol(j,ku,i); volkl = volk0(j,i,1);
                xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                double ttpo2 = xp*xa+yp*ya+zp*za;
                double ttmo2 = xm*xa+ym*ya+zm*za;
                double ttpn2 = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                double ttmn2 = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                ttp = ttpo2*(1-istrongturbdis)+ttpn2*istrongturbdis;
                ttm = ttmo2*(1-istrongturbdis)+ttmn2*istrongturbdis;
                anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                fnup = 0.5*(fnu(j,k+1,i)+fnu(j,k,i));
                fnum = 0.5*(fnu(j,k-1,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                by(j,k) = -cdm; cy(j,k) = cdp+cdm; dy(j,k) = -cdp;
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                by2(j,k) = -cdm; cy2(j,k) = cdp+cdm; dy2(j,k) = -cdp;
                cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                by3(j,k) = -cdm; cy3(j,k) = cdp+cdm; dy3(j,k) = -cdp;
                cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                by4(j,k) = -cdm; cy4(j,k) = cdp+cdm; dy4(j,k) = -cdp;
            }
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                by(j,k) -= uu*app; by2(j,k) -= uu*app; by3(j,k) -= uu*app; by4(j,k) -= uu*app;
                cy(j,k) += uu*(app-apm); cy2(j,k) += uu*(app-apm); cy3(j,k) += uu*(app-apm); cy4(j,k) += uu*(app-apm);
                dy(j,k) += uu*apm; dy2(j,k) += uu*apm; dy3(j,k) += uu*apm; dy4(j,k) += uu*apm;
            }
            for (j = 1; j <= jdim-1; j++) {
                double cyadd2=0.0, cy2add2=0.0, cy3add2=0.0, cy4add2=0.0;
                if (ivmx == 40) {
                    double cmuc_l = blend(j,k,i)*cmuc1+(1.0-blend(j,k,i))*cmuc2;
                    double betax_l = blend(j,k,i)*beta1+(1.0-blend(j,k,i))*beta2;
                    double damp1_abs; { double av=damp1(j,k,i); damp1_abs = ccabs(av); }
                    cyadd2 = 2.0*re*betax_l*turre(j,k,i,1)+damp1_abs/turre(j,k,i,1);
                    cy2add2 = re*cmuc_l*turre(j,k,i,1);
                    double vor_abs; { double av=2.0*ce2*turre(j,k,i,3)-1.0; vor_abs = ccabs(av); }
                    cy3add2 = ca2*vor(j,k,i)*vor_abs;
                    cy4add2 = wrks(j,k,i);
                }
                cy(j,k) += cyadd2; cy2(j,k) += cy2add2; cy3(j,k) += cy3add2; cy4(j,k) += cy4add2;
            }
            for (j = 1; j <= jdim-1; j++) {
                by(j,k) *= timestp(j,k,i); by2(j,k) *= timestp(j,k,i)*factor2;
                by3(j,k) *= timestp(j,k,i)*factor3; by4(j,k) *= timestp(j,k,i)*factor4;
                cy(j,k) = cy(j,k)*timestp(j,k,i)+1.0*(1.0+phi);
                cy2(j,k) = cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.0+phi);
                cy3(j,k) = cy3(j,k)*timestp(j,k,i)*factor3+1.0*(1.0+phi);
                cy4(j,k) = cy4(j,k)*timestp(j,k,i)*factor4+1.0*(1.0+phi);
                dy(j,k) *= timestp(j,k,i); dy2(j,k) *= timestp(j,k,i)*factor2;
                dy3(j,k) *= timestp(j,k,i)*factor3; dy4(j,k) *= timestp(j,k,i)*factor4;
                fy(j,k) = rhside(j,k,i,1)*timestp(j,k,i);
                fy2(j,k) = rhside(j,k,i,2)*timestp(j,k,i)*factor2;
                fy3(j,k) = rhside(j,k,i,3)*timestp(j,k,i)*factor3;
                fy4(j,k) = rhside(j,k,i,4)*timestp(j,k,i)*factor4;
            }
            if ((float)dt > 0.0f) {
                for (j = 1; j <= jdim-1; j++) {
                    fy(j,k) += (1.0+phi)*(zksav2(j,k,i,1)-turre(j,k,i,1))+phi*zksav2(j,k,i,5);
                    fy2(j,k) += (1.0+phi)*(zksav2(j,k,i,2)-turre(j,k,i,2))+phi*zksav2(j,k,i,6);
                    fy3(j,k) += (1.0+phi)*(zksav2(j,k,i,3)-turre(j,k,i,3))+phi*zksav2(j,k,i,7);
                    fy4(j,k) += (1.0+phi)*(zksav2(j,k,i,4)-turre(j,k,i,4))+phi*zksav2(j,k,i,8);
                }
            }


            // KDIM boundary points
            k = kdim-1; kl = kdim-2; ku = kdim-1;
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
                dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
                sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
                volku = volk0(j,i,3); volkl = vol(j,kl,i);
                xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                double ttpo2 = xp*xa+yp*ya+zp*za;
                double ttmo2 = xm*xa+ym*ya+zm*za;
                double ttpn2 = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                double ttmn2 = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                ttp = ttpo2*(1-istrongturbdis)+ttpn2*istrongturbdis;
                ttm = ttmo2*(1-istrongturbdis)+ttmn2*istrongturbdis;
                anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                fnup = 0.5*(fnu(j,k+1,i)+fnu(j,k,i));
                fnum = 0.5*(fnu(j,k-1,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                by(j,k) = -cdm; cy(j,k) = cdp+cdm; dy(j,k) = -cdp;
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                by2(j,k) = -cdm; cy2(j,k) = cdp+cdm; dy2(j,k) = -cdp;
                cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                by3(j,k) = -cdm; cy3(j,k) = cdp+cdm; dy3(j,k) = -cdp;
                cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                by4(j,k) = -cdm; cy4(j,k) = cdp+cdm; dy4(j,k) = -cdp;
            }
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                by(j,k) -= uu*app; by2(j,k) -= uu*app; by3(j,k) -= uu*app; by4(j,k) -= uu*app;
                cy(j,k) += uu*(app-apm); cy2(j,k) += uu*(app-apm); cy3(j,k) += uu*(app-apm); cy4(j,k) += uu*(app-apm);
                dy(j,k) += uu*apm; dy2(j,k) += uu*apm; dy3(j,k) += uu*apm; dy4(j,k) += uu*apm;
            }
            for (j = 1; j <= jdim-1; j++) {
                double cyadd2=0.0, cy2add2=0.0, cy3add2=0.0, cy4add2=0.0;
                if (ivmx == 40) {
                    double cmuc_l = blend(j,k,i)*cmuc1+(1.0-blend(j,k,i))*cmuc2;
                    double betax_l = blend(j,k,i)*beta1+(1.0-blend(j,k,i))*beta2;
                    double damp1_abs; { double av=damp1(j,k,i); damp1_abs = ccabs(av); }
                    cyadd2 = 2.0*re*betax_l*turre(j,k,i,1)+damp1_abs/turre(j,k,i,1);
                    cy2add2 = re*cmuc_l*turre(j,k,i,1);
                    double vor_abs; { double av=2.0*ce2*turre(j,k,i,3)-1.0; vor_abs = ccabs(av); }
                    cy3add2 = ca2*vor(j,k,i)*vor_abs;
                    cy4add2 = wrks(j,k,i);
                }
                cy(j,k) += cyadd2; cy2(j,k) += cy2add2; cy3(j,k) += cy3add2; cy4(j,k) += cy4add2;
            }
            for (j = 1; j <= jdim-1; j++) {
                by(j,k) *= timestp(j,k,i); by2(j,k) *= timestp(j,k,i)*factor2;
                by3(j,k) *= timestp(j,k,i)*factor3; by4(j,k) *= timestp(j,k,i)*factor4;
                cy(j,k) = cy(j,k)*timestp(j,k,i)+1.0*(1.0+phi);
                cy2(j,k) = cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.0+phi);
                cy3(j,k) = cy3(j,k)*timestp(j,k,i)*factor3+1.0*(1.0+phi);
                cy4(j,k) = cy4(j,k)*timestp(j,k,i)*factor4+1.0*(1.0+phi);
                dy(j,k) *= timestp(j,k,i); dy2(j,k) *= timestp(j,k,i)*factor2;
                dy3(j,k) *= timestp(j,k,i)*factor3; dy4(j,k) *= timestp(j,k,i)*factor4;
                fy(j,k) = rhside(j,k,i,1)*timestp(j,k,i);
                fy2(j,k) = rhside(j,k,i,2)*timestp(j,k,i)*factor2;
                fy3(j,k) = rhside(j,k,i,3)*timestp(j,k,i)*factor3;
                fy4(j,k) = rhside(j,k,i,4)*timestp(j,k,i)*factor4;
            }
            if ((float)dt > 0.0f) {
                for (j = 1; j <= jdim-1; j++) {
                    fy(j,k) += (1.0+phi)*(zksav2(j,k,i,1)-turre(j,k,i,1))+phi*zksav2(j,k,i,5);
                    fy2(j,k) += (1.0+phi)*(zksav2(j,k,i,2)-turre(j,k,i,2))+phi*zksav2(j,k,i,6);
                    fy3(j,k) += (1.0+phi)*(zksav2(j,k,i,3)-turre(j,k,i,3))+phi*zksav2(j,k,i,7);
                    fy4(j,k) += (1.0+phi)*(zksav2(j,k,i,4)-turre(j,k,i,4))+phi*zksav2(j,k,i,8);
                }
            }
            // Overset blanking
            if (iover == 1) {
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    fy(j,k) *= blank(j,k,i); by(j,k) *= blank(j,k,i); dy(j,k) *= blank(j,k,i);
                    cy(j,k) = cy(j,k)*blank(j,k,i)+(1.0-blank(j,k,i));
                    fy2(j,k) *= blank(j,k,i); by2(j,k) *= blank(j,k,i); dy2(j,k) *= blank(j,k,i);
                    cy2(j,k) = cy2(j,k)*blank(j,k,i)+(1.0-blank(j,k,i));
                    fy3(j,k) *= blank(j,k,i); by3(j,k) *= blank(j,k,i); dy3(j,k) *= blank(j,k,i);
                    cy3(j,k) = cy3(j,k)*blank(j,k,i)+(1.0-blank(j,k,i));
                    fy4(j,k) *= blank(j,k,i); by4(j,k) *= blank(j,k,i); dy4(j,k) *= blank(j,k,i);
                    cy4(j,k) = cy4(j,k)*blank(j,k,i)+(1.0-blank(j,k,i));
                }
            }
            {
                int jl2=1, ju2=jdim-1, kl2=1, ku2=kdim-1;
                int jdim1=jdim-1, kdim1=kdim-1;
                triv(jdim1, kdim1, jl2, ju2, kl2, ku2, worky, by, cy, dy, fy);
                triv(jdim1, kdim1, jl2, ju2, kl2, ku2, worky, by2, cy2, dy2, fy2);
                triv(jdim1, kdim1, jl2, ju2, kl2, ku2, worky, by3, cy3, dy3, fy3);
                triv(jdim1, kdim1, jl2, ju2, kl2, ku2, worky, by4, cy4, dy4, fy4);
            }
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                rhside(j,k,i,1) = fy(j,k);
                rhside(j,k,i,2) = fy2(j,k);
                rhside(j,k,i,3) = fy3(j,k);
                rhside(j,k,i,4) = fy4(j,k);
            }
        } // end i loop for implicit eta


        // Implicit F_xi_xi viscous terms. Do over all i's
        for (i = 1; i <= idim-1; i++) {
            // Interior points
            for (j = 2; j <= jdim-2; j++) {
                jl = j-1; ju = j+1;
                for (k = 1; k <= kdim-1; k++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
                    dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
                    sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
                    volju = vol(ju,k,i); voljl = vol(jl,k,i);
                    xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    double ttpo2 = xp*xa+yp*ya+zp*za;
                    double ttmo2 = xm*xa+ym*ya+zm*za;
                    double ttpn2 = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                    double ttmn2 = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                    ttp = ttpo2*(1-istrongturbdis)+ttpn2*istrongturbdis;
                    ttm = ttmo2*(1-istrongturbdis)+ttmn2*istrongturbdis;
                    anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                    anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                    fnup = 0.5*(fnu(j+1,k,i)+fnu(j,k,i));
                    fnum = 0.5*(fnu(j-1,k,i)+fnu(j,k,i));
                    cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    bx(k,j) = -cdm; cx(k,j) = cdp+cdm; dx(k,j) = -cdp;
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    bx2(k,j) = -cdm; cx2(k,j) = cdp+cdm; dx2(k,j) = -cdp;
                    cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                    bx3(k,j) = -cdm; cx3(k,j) = cdp+cdm; dx3(k,j) = -cdp;
                    cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                    cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                    bx4(k,j) = -cdm; cx4(k,j) = cdp+cdm; dx4(k,j) = -cdp;
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                    app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                    bx(k,j) -= uu*app; bx2(k,j) -= uu*app; bx3(k,j) -= uu*app; bx4(k,j) -= uu*app;
                    cx(k,j) += uu*(app-apm); cx2(k,j) += uu*(app-apm); cx3(k,j) += uu*(app-apm); cx4(k,j) += uu*(app-apm);
                    dx(k,j) += uu*apm; dx2(k,j) += uu*apm; dx3(k,j) += uu*apm; dx4(k,j) += uu*apm;
                }
                for (k = 1; k <= kdim-1; k++) {
                    bx(k,j) *= timestp(j,k,i); bx2(k,j) *= timestp(j,k,i)*factor2;
                    bx3(k,j) *= timestp(j,k,i)*factor3; bx4(k,j) *= timestp(j,k,i)*factor4;
                    cx(k,j) = cx(k,j)*timestp(j,k,i)+1.0*(1.0+phi);
                    cx2(k,j) = cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.0+phi);
                    cx3(k,j) = cx3(k,j)*timestp(j,k,i)*factor3+1.0*(1.0+phi);
                    cx4(k,j) = cx4(k,j)*timestp(j,k,i)*factor4+1.0*(1.0+phi);
                    dx(k,j) *= timestp(j,k,i); dx2(k,j) *= timestp(j,k,i)*factor2;
                    dx3(k,j) *= timestp(j,k,i)*factor3; dx4(k,j) *= timestp(j,k,i)*factor4;
                    fx(k,j) = rhside(j,k,i,1)*(1.0+phi);
                    fx2(k,j) = rhside(j,k,i,2)*(1.0+phi);
                    fx3(k,j) = rhside(j,k,i,3)*(1.0+phi);
                    fx4(k,j) = rhside(j,k,i,4)*(1.0+phi);
                }
            }


            // J0 boundary points
            j = 1; jl = 1; ju = std::min(2, jdim-1);
            for (k = 1; k <= kdim-1; k++) {
                dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
                dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
                sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
                volju = vol(ju,k,i); voljl = volj0(k,i,1);
                xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                double ttpo2 = xp*xa+yp*ya+zp*za;
                double ttmo2 = xm*xa+ym*ya+zm*za;
                double ttpn2 = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                double ttmn2 = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                ttp = ttpo2*(1-istrongturbdis)+ttpn2*istrongturbdis;
                ttm = ttmo2*(1-istrongturbdis)+ttmn2*istrongturbdis;
                anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                fnup = 0.5*(fnu(j+1,k,i)+fnu(j,k,i));
                fnum = 0.5*(fnu(j-1,k,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                bx(k,j) = -cdm; cx(k,j) = cdp+cdm; dx(k,j) = -cdp;
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                bx2(k,j) = -cdm; cx2(k,j) = cdp+cdm; dx2(k,j) = -cdp;
                cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                bx3(k,j) = -cdm; cx3(k,j) = cdp+cdm; dx3(k,j) = -cdp;
                cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                bx4(k,j) = -cdm; cx4(k,j) = cdp+cdm; dx4(k,j) = -cdp;
            }
            for (k = 1; k <= kdim-1; k++) {
                xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                bx(k,j) -= uu*app; bx2(k,j) -= uu*app; bx3(k,j) -= uu*app; bx4(k,j) -= uu*app;
                cx(k,j) += uu*(app-apm); cx2(k,j) += uu*(app-apm); cx3(k,j) += uu*(app-apm); cx4(k,j) += uu*(app-apm);
                dx(k,j) += uu*apm; dx2(k,j) += uu*apm; dx3(k,j) += uu*apm; dx4(k,j) += uu*apm;
            }
            for (k = 1; k <= kdim-1; k++) {
                bx(k,j) *= timestp(j,k,i); bx2(k,j) *= timestp(j,k,i)*factor2;
                bx3(k,j) *= timestp(j,k,i)*factor3; bx4(k,j) *= timestp(j,k,i)*factor4;
                cx(k,j) = cx(k,j)*timestp(j,k,i)+1.0*(1.0+phi);
                cx2(k,j) = cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.0+phi);
                cx3(k,j) = cx3(k,j)*timestp(j,k,i)*factor3+1.0*(1.0+phi);
                cx4(k,j) = cx4(k,j)*timestp(j,k,i)*factor4+1.0*(1.0+phi);
                dx(k,j) *= timestp(j,k,i); dx2(k,j) *= timestp(j,k,i)*factor2;
                dx3(k,j) *= timestp(j,k,i)*factor3; dx4(k,j) *= timestp(j,k,i)*factor4;
                fx(k,j) = rhside(j,k,i,1)*(1.0+phi);
                fx2(k,j) = rhside(j,k,i,2)*(1.0+phi);
                fx3(k,j) = rhside(j,k,i,3)*(1.0+phi);
                fx4(k,j) = rhside(j,k,i,4)*(1.0+phi);
            }


            // JDIM boundary points
            j = jdim-1; jl = jdim-2; ju = jdim-1;
            for (k = 1; k <= kdim-1; k++) {
                dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
                dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
                sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
                volju = volj0(k,i,3); voljl = vol(jl,k,i);
                xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                double ttpo2 = xp*xa+yp*ya+zp*za;
                double ttmo2 = xm*xa+ym*ya+zm*za;
                double ttpn2 = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                double ttmn2 = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                ttp = ttpo2*(1-istrongturbdis)+ttpn2*istrongturbdis;
                ttm = ttmo2*(1-istrongturbdis)+ttmn2*istrongturbdis;
                anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                fnup = 0.5*(fnu(j+1,k,i)+fnu(j,k,i));
                fnum = 0.5*(fnu(j-1,k,i)+fnu(j,k,i));
                cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                bx(k,j) = -cdm; cx(k,j) = cdp+cdm; dx(k,j) = -cdp;
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                bx2(k,j) = -cdm; cx2(k,j) = cdp+cdm; dx2(k,j) = -cdp;
                cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                bx3(k,j) = -cdm; cx3(k,j) = cdp+cdm; dx3(k,j) = -cdp;
                cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                bx4(k,j) = -cdm; cx4(k,j) = cdp+cdm; dx4(k,j) = -cdp;
            }
            for (k = 1; k <= kdim-1; k++) {
                xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                bx(k,j) -= uu*app; bx2(k,j) -= uu*app; bx3(k,j) -= uu*app; bx4(k,j) -= uu*app;
                cx(k,j) += uu*(app-apm); cx2(k,j) += uu*(app-apm); cx3(k,j) += uu*(app-apm); cx4(k,j) += uu*(app-apm);
                dx(k,j) += uu*apm; dx2(k,j) += uu*apm; dx3(k,j) += uu*apm; dx4(k,j) += uu*apm;
            }
            for (k = 1; k <= kdim-1; k++) {
                bx(k,j) *= timestp(j,k,i); bx2(k,j) *= timestp(j,k,i)*factor2;
                bx3(k,j) *= timestp(j,k,i)*factor3; bx4(k,j) *= timestp(j,k,i)*factor4;
                cx(k,j) = cx(k,j)*timestp(j,k,i)+1.0*(1.0+phi);
                cx2(k,j) = cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.0+phi);
                cx3(k,j) = cx3(k,j)*timestp(j,k,i)*factor3+1.0*(1.0+phi);
                cx4(k,j) = cx4(k,j)*timestp(j,k,i)*factor4+1.0*(1.0+phi);
                dx(k,j) *= timestp(j,k,i); dx2(k,j) *= timestp(j,k,i)*factor2;
                dx3(k,j) *= timestp(j,k,i)*factor3; dx4(k,j) *= timestp(j,k,i)*factor4;
                fx(k,j) = rhside(j,k,i,1)*(1.0+phi);
                fx2(k,j) = rhside(j,k,i,2)*(1.0+phi);
                fx3(k,j) = rhside(j,k,i,3)*(1.0+phi);
                fx4(k,j) = rhside(j,k,i,4)*(1.0+phi);
            }
            // Overset blanking
            if (iover == 1) {
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    fx(k,j) *= blank(j,k,i); bx(k,j) *= blank(j,k,i); dx(k,j) *= blank(j,k,i);
                    cx(k,j) = cx(k,j)*blank(j,k,i)+(1.0-blank(j,k,i));
                    fx2(k,j) *= blank(j,k,i); bx2(k,j) *= blank(j,k,i); dx2(k,j) *= blank(j,k,i);
                    cx2(k,j) = cx2(k,j)*blank(j,k,i)+(1.0-blank(j,k,i));
                    fx3(k,j) *= blank(j,k,i); bx3(k,j) *= blank(j,k,i); dx3(k,j) *= blank(j,k,i);
                    cx3(k,j) = cx3(k,j)*blank(j,k,i)+(1.0-blank(j,k,i));
                    fx4(k,j) *= blank(j,k,i); bx4(k,j) *= blank(j,k,i); dx4(k,j) *= blank(j,k,i);
                    cx4(k,j) = cx4(k,j)*blank(j,k,i)+(1.0-blank(j,k,i));
                }
            }
            {
                int kl2=1, ku2=kdim-1, jl2=1, ju2=jdim-1;
                int kdim1=kdim-1, jdim1=jdim-1;
                triv(kdim1, jdim1, kl2, ku2, jl2, ju2, workx, bx, cx, dx, fx);
                triv(kdim1, jdim1, kl2, ku2, jl2, ju2, workx, bx2, cx2, dx2, fx2);
                triv(kdim1, jdim1, kl2, ku2, jl2, ju2, workx, bx3, cx3, dx3, fx3);
                triv(kdim1, jdim1, kl2, ku2, jl2, ju2, workx, bx4, cx4, dx4, fx4);
            }
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                rhside(j,k,i,1) = fx(k,j);
                rhside(j,k,i,2) = fx2(k,j);
                rhside(j,k,i,3) = fx3(k,j);
                rhside(j,k,i,4) = fx4(k,j);
            }
        } // end i loop for implicit xi


        // Implicit F_zeta_zeta viscous terms. Do over all j's
        if (i2d != 1 && iaxi2planeturb != 1) {
            for (j = 1; j <= jdim-1; j++) {
                // Interior points
                for (i = 2; i <= idim-2; i++) {
                    il = i-1; iu = i+1;
                    for (k = 1; k <= kdim-1; k++) {
                        dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                        dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                        sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                        sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                        sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                        sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
                        voliu = vol(j,k,iu); volil = vol(j,k,il);
                        xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        double ttpo2 = xp*xa+yp*ya+zp*za;
                        double ttmo2 = xm*xa+ym*ya+zm*za;
                        double ttpn2 = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                        double ttmn2 = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                        ttp = ttpo2*(1-istrongturbdis)+ttpn2*istrongturbdis;
                        ttm = ttmo2*(1-istrongturbdis)+ttmn2*istrongturbdis;
                        anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                        anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                        fnup = 0.5*(fnu(j,k,i+1)+fnu(j,k,i));
                        fnum = 0.5*(fnu(j,k,i-1)+fnu(j,k,i));
                        cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                        bz(k,i) = -cdm; cz(k,i) = cdp+cdm; dz(k,i) = -cdp;
                        cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                        bz2(k,i) = -cdm; cz2(k,i) = cdp+cdm; dz2(k,i) = -cdp;
                        cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                        bz3(k,i) = -cdm; cz3(k,i) = cdp+cdm; dz3(k,i) = -cdp;
                        cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                        cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                        bz4(k,i) = -cdm; cz4(k,i) = cdp+cdm; dz4(k,i) = -cdp;
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                        uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                        { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                        app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                        bz(k,i) -= uu*app; bz2(k,i) -= uu*app; bz3(k,i) -= uu*app; bz4(k,i) -= uu*app;
                        cz(k,i) += uu*(app-apm); cz2(k,i) += uu*(app-apm); cz3(k,i) += uu*(app-apm); cz4(k,i) += uu*(app-apm);
                        dz(k,i) += uu*apm; dz2(k,i) += uu*apm; dz3(k,i) += uu*apm; dz4(k,i) += uu*apm;
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        bz(k,i) *= timestp(j,k,i); bz2(k,i) *= timestp(j,k,i)*factor2;
                        bz3(k,i) *= timestp(j,k,i)*factor3; bz4(k,i) *= timestp(j,k,i)*factor4;
                        cz(k,i) = cz(k,i)*timestp(j,k,i)+1.0*(1.0+phi);
                        cz2(k,i) = cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.0+phi);
                        cz3(k,i) = cz3(k,i)*timestp(j,k,i)*factor3+1.0*(1.0+phi);
                        cz4(k,i) = cz4(k,i)*timestp(j,k,i)*factor4+1.0*(1.0+phi);
                        dz(k,i) *= timestp(j,k,i); dz2(k,i) *= timestp(j,k,i)*factor2;
                        dz3(k,i) *= timestp(j,k,i)*factor3; dz4(k,i) *= timestp(j,k,i)*factor4;
                        fz(k,i) = rhside(j,k,i,1)*(1.0+phi);
                        fz2(k,i) = rhside(j,k,i,2)*(1.0+phi);
                        fz3(k,i) = rhside(j,k,i,3)*(1.0+phi);
                        fz4(k,i) = rhside(j,k,i,4)*(1.0+phi);
                    }
                }


                // I0 boundary points
                i = 1; il = 1; iu = std::min(2, idim-1);
                for (k = 1; k <= kdim-1; k++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                    dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                    sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
                    voliu = vol(j,k,iu); volil = voli0(j,k,1);
                    xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    double ttpo2 = xp*xa+yp*ya+zp*za;
                    double ttmo2 = xm*xa+ym*ya+zm*za;
                    double ttpn2 = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    double ttmn2 = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    ttp = ttpo2*(1-istrongturbdis)+ttpn2*istrongturbdis;
                    ttm = ttmo2*(1-istrongturbdis)+ttmn2*istrongturbdis;
                    anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                    anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                    fnup = 0.5*(fnu(j,k,i+1)+fnu(j,k,i));
                    fnum = 0.5*(fnu(j,k,i-1)+fnu(j,k,i));
                    cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    bz(k,i) = -cdm; cz(k,i) = cdp+cdm; dz(k,i) = -cdp;
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    bz2(k,i) = -cdm; cz2(k,i) = cdp+cdm; dz2(k,i) = -cdp;
                    cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                    bz3(k,i) = -cdm; cz3(k,i) = cdp+cdm; dz3(k,i) = -cdp;
                    cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                    cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                    bz4(k,i) = -cdm; cz4(k,i) = cdp+cdm; dz4(k,i) = -cdp;
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                    app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                    bz(k,i) -= uu*app; bz2(k,i) -= uu*app; bz3(k,i) -= uu*app; bz4(k,i) -= uu*app;
                    cz(k,i) += uu*(app-apm); cz2(k,i) += uu*(app-apm); cz3(k,i) += uu*(app-apm); cz4(k,i) += uu*(app-apm);
                    dz(k,i) += uu*apm; dz2(k,i) += uu*apm; dz3(k,i) += uu*apm; dz4(k,i) += uu*apm;
                }
                for (k = 1; k <= kdim-1; k++) {
                    bz(k,i) *= timestp(j,k,i); bz2(k,i) *= timestp(j,k,i)*factor2;
                    bz3(k,i) *= timestp(j,k,i)*factor3; bz4(k,i) *= timestp(j,k,i)*factor4;
                    cz(k,i) = cz(k,i)*timestp(j,k,i)+1.0*(1.0+phi);
                    cz2(k,i) = cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.0+phi);
                    cz3(k,i) = cz3(k,i)*timestp(j,k,i)*factor3+1.0*(1.0+phi);
                    cz4(k,i) = cz4(k,i)*timestp(j,k,i)*factor4+1.0*(1.0+phi);
                    dz(k,i) *= timestp(j,k,i); dz2(k,i) *= timestp(j,k,i)*factor2;
                    dz3(k,i) *= timestp(j,k,i)*factor3; dz4(k,i) *= timestp(j,k,i)*factor4;
                    fz(k,i) = rhside(j,k,i,1)*(1.0+phi);
                    fz2(k,i) = rhside(j,k,i,2)*(1.0+phi);
                    fz3(k,i) = rhside(j,k,i,3)*(1.0+phi);
                    fz4(k,i) = rhside(j,k,i,4)*(1.0+phi);
                }


                // IDIM boundary points
                i = idim-1; il = idim-2; iu = idim-1;
                for (k = 1; k <= kdim-1; k++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                    dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                    sigkp = dfacep*sigk1+(1.0-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.0-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.0-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.0-dfacem)*sigo2;
                    voliu = voli0(j,k,3); volil = vol(j,k,il);
                    xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    double ttpo2 = xp*xa+yp*ya+zp*za;
                    double ttmo2 = xm*xa+ym*ya+zm*za;
                    double ttpn2 = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    double ttmn2 = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    ttp = ttpo2*(1-istrongturbdis)+ttpn2*istrongturbdis;
                    ttm = ttmo2*(1-istrongturbdis)+ttmn2*istrongturbdis;
                    anutp = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                    anutm = 0.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                    fnup = 0.5*(fnu(j,k,i+1)+fnu(j,k,i));
                    fnum = 0.5*(fnu(j,k,i-1)+fnu(j,k,i));
                    cdp = (fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    bz(k,i) = -cdm; cz(k,i) = cdp+cdm; dz(k,i) = -cdp;
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    bz2(k,i) = -cdm; cz2(k,i) = cdp+cdm; dz2(k,i) = -cdp;
                    cdp = (fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm = (fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re);
                    bz3(k,i) = -cdm; cz3(k,i) = cdp+cdm; dz3(k,i) = -cdp;
                    cdp = (sigma_thetat*(fnup+anutp))*ttp/(q(j,k,i,1)*re);
                    cdm = (sigma_thetat*(fnum+anutm))*ttm/(q(j,k,i,1)*re);
                    bz4(k,i) = -cdm; cz4(k,i) = cdp+cdm; dz4(k,i) = -cdp;
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    { double av1=1.0, av2=uu; sgnu = ccsignrc(av1,av2); }
                    app = 0.5*(1.0+sgnu); apm = 0.5*(1.0-sgnu);
                    bz(k,i) -= uu*app; bz2(k,i) -= uu*app; bz3(k,i) -= uu*app; bz4(k,i) -= uu*app;
                    cz(k,i) += uu*(app-apm); cz2(k,i) += uu*(app-apm); cz3(k,i) += uu*(app-apm); cz4(k,i) += uu*(app-apm);
                    dz(k,i) += uu*apm; dz2(k,i) += uu*apm; dz3(k,i) += uu*apm; dz4(k,i) += uu*apm;
                }
                for (k = 1; k <= kdim-1; k++) {
                    bz(k,i) *= timestp(j,k,i); bz2(k,i) *= timestp(j,k,i)*factor2;
                    bz3(k,i) *= timestp(j,k,i)*factor3; bz4(k,i) *= timestp(j,k,i)*factor4;
                    cz(k,i) = cz(k,i)*timestp(j,k,i)+1.0*(1.0+phi);
                    cz2(k,i) = cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.0+phi);
                    cz3(k,i) = cz3(k,i)*timestp(j,k,i)*factor3+1.0*(1.0+phi);
                    cz4(k,i) = cz4(k,i)*timestp(j,k,i)*factor4+1.0*(1.0+phi);
                    dz(k,i) *= timestp(j,k,i); dz2(k,i) *= timestp(j,k,i)*factor2;
                    dz3(k,i) *= timestp(j,k,i)*factor3; dz4(k,i) *= timestp(j,k,i)*factor4;
                    fz(k,i) = rhside(j,k,i,1)*(1.0+phi);
                    fz2(k,i) = rhside(j,k,i,2)*(1.0+phi);
                    fz3(k,i) = rhside(j,k,i,3)*(1.0+phi);
                    fz4(k,i) = rhside(j,k,i,4)*(1.0+phi);
                }
                // Overset blanking
                if (iover == 1) {
                    for (k = 1; k <= kdim-1; k++)
                    for (i = 1; i <= idim-1; i++) {
                        fz(k,i) *= blank(j,k,i); bz(k,i) *= blank(j,k,i); dz(k,i) *= blank(j,k,i);
                        cz(k,i) = cz(k,i)*blank(j,k,i)+(1.0-blank(j,k,i));
                        fz2(k,i) *= blank(j,k,i); bz2(k,i) *= blank(j,k,i); dz2(k,i) *= blank(j,k,i);
                        cz2(k,i) = cz2(k,i)*blank(j,k,i)+(1.0-blank(j,k,i));
                        fz3(k,i) *= blank(j,k,i); bz3(k,i) *= blank(j,k,i); dz3(k,i) *= blank(j,k,i);
                        cz3(k,i) = cz3(k,i)*blank(j,k,i)+(1.0-blank(j,k,i));
                        fz4(k,i) *= blank(j,k,i); bz4(k,i) *= blank(j,k,i); dz4(k,i) *= blank(j,k,i);
                        cz4(k,i) = cz4(k,i)*blank(j,k,i)+(1.0-blank(j,k,i));
                    }
                }
                {
                    int kl2=1, ku2=kdim-1, il2=1, iu2=idim-1;
                    int kdim1=kdim-1, idim1=idim-1;
                    triv(kdim1, idim1, kl2, ku2, il2, iu2, workz, bz, cz, dz, fz);
                    triv(kdim1, idim1, kl2, ku2, il2, iu2, workz, bz2, cz2, dz2, fz2);
                    triv(kdim1, idim1, kl2, ku2, il2, iu2, workz, bz3, cz3, dz3, fz3);
                    triv(kdim1, idim1, kl2, ku2, il2, iu2, workz, bz4, cz4, dz4, fz4);
                }
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++) {
                    rhside(j,k,i,1) = fz(k,i);
                    rhside(j,k,i,2) = fz2(k,i);
                    rhside(j,k,i,3) = fz3(k,i);
                    rhside(j,k,i,4) = fz4(k,i);
                }
            } // end j loop for implicit zeta
        } // end if i2d


        // Update TURRE
        sumno = 0.0; sumnk = 0.0; sumni = 0.0; sumnr = 0.0;
        negno = 0; negnk = 0; negni = 0; negnr = 0;
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            sumno += rhside(j,k,i,1)*rhside(j,k,i,1);
            if ((float)(turre(j,k,i,1)+rhside(j,k,i,1)) <= (float)tur1cutlev) {
                negno++;
                if ((float)tur1cut > 0.0f) turre(j,k,i,1) = (double)tur1cut;
            } else {
                turre(j,k,i,1) = turre(j,k,i,1)+rhside(j,k,i,1);
            }
            sumnk += rhside(j,k,i,2)*rhside(j,k,i,2);
            if ((float)(turre(j,k,i,2)+rhside(j,k,i,2)) <= (float)tur2cutlev) {
                negnk++;
                if ((float)tur2cut > 0.0f) turre(j,k,i,2) = (double)tur2cut;
            } else {
                turre(j,k,i,2) = turre(j,k,i,2)+rhside(j,k,i,2);
            }
            sumni += rhside(j,k,i,3)*rhside(j,k,i,3);
            if (itrans_on == 1) {
                if ((float)(turre(j,k,i,3)+rhside(j,k,i,3)) <= 0.0f) {
                    negni++;
                    turre(j,k,i,3) = 0.0;
                } else {
                    turre(j,k,i,3) = turre(j,k,i,3)+rhside(j,k,i,3);
                }
            }
            sumnr += rhside(j,k,i,4)*rhside(j,k,i,4);
            if (itrans_on == 1) {
                if ((float)(turre(j,k,i,4)+rhside(j,k,i,4)) <= 20.0f) {
                    negnr++;
                    turre(j,k,i,4) = 20.0;
                } else {
                    turre(j,k,i,4) = turre(j,k,i,4)+rhside(j,k,i,4);
                }
            }
        }
        sumno = std::sqrt(sumno)/(double)((kdim-1)*(jdim-1)*(idim-1));
        sumnk = std::sqrt(sumnk)/(double)((kdim-1)*(jdim-1)*(idim-1));
        sumni = std::sqrt(sumni)/(double)((kdim-1)*(jdim-1)*(idim-1));
        sumnr = std::sqrt(sumnr)/(double)((kdim-1)*(jdim-1)*(idim-1));

    } // end not_ loop (nsubit)

    sumn1 = sumno; sumn2 = sumnk; sumn3 = sumni; sumn4 = sumnr;
    negn1 = negno; negn2 = negnk; negn3 = negni; negn4 = negnr;


    // Update VIST3D and save omega and k values
    if (ivmx == 40) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            double smin_abs_v = smin(j,k,i); double smin_abs = ccabs(smin_abs_v);
            denom1 = 2.0/re*std::sqrt(turre(j,k,i,2))/(.09*turre(j,k,i,1)*smin_abs);
            denom2 = 500.0*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*re*re*smin(j,k,i)*turre(j,k,i,1));
            { double av1=denom1, av2=denom2; arg2 = ccmax(av1,av2); }
            { double av=arg2*arg2; f2 = cctanh(av); }
            if (isstdenom == 1) {
                s11 = ux(j,k,i,1); s22 = ux(j,k,i,5); s33 = ux(j,k,i,9);
                s12 = 0.5*(ux(j,k,i,2)+ux(j,k,i,4));
                s13 = 0.5*(ux(j,k,i,3)+ux(j,k,i,7));
                s23 = 0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                xis = s11*s11+s22*s22+s33*s33+2.0*s12*s12+2.0*s13*s13+2.0*s23*s23;
                xxx = std::sqrt(2.0*xis);
                { double av1=a1*turre(j,k,i,1), av2=xxx*f2/re; denom = ccmax(av1,av2); }
            } else {
                { double av1=a1*turre(j,k,i,1), av2=vor(j,k,i)*f2/re; denom = ccmax(av1,av2); }
            }
            vist3d(j,k,i) = a1*q(j,k,i,1)*turre(j,k,i,2)/denom;
            { double av1=vist3d(j,k,i), av2=(double)edvislim; vist3d(j,k,i) = ccmin(av1,av2); }
            zksav(j,k,i,1) = turre(j,k,i,1);
            zksav(j,k,i,2) = turre(j,k,i,2);
            zksav(j,k,i,3) = turre(j,k,i,3);
            zksav(j,k,i,4) = turre(j,k,i,4);
        }
    }

    if (i_lam_forcezero == 1) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            if ((i >= ilamlo && i < ilamhi &&
                 j >= jlamlo && j < jlamhi &&
                 k >= klamlo && k < klamhi) ||
                (float)smin(j,k,i) < 0.0f) {
                vist3d(j,k,i) = 0.0;
            }
        }
    }

    // Hardwired output section (iwriteaux=0, so this is always skipped)
    iwriteaux = 0;
    if (iwriteaux == 1) {
        if (icyc == cmn_info.ncyc1[0] || icyc == cmn_info.ncyc1[1] ||
            icyc == cmn_info.ncyc1[2] || icyc == cmn_info.ncyc1[3] ||
            icyc == cmn_info.ncyc1[4]) {
            if (ntime == 1) {
                jset = (int)(.779*(float)(jdim));
                FILE* f92 = fortran_get_unit(92);
                fprintf(f92, "%5d jset,x=%5d%12.5e  u+,k+,e+,-uv+,logy+\n",
                    kdim-1, jset, (float)(0.5*(x(jset,1,1)+x(jset+1,1,1))));
                double qset2 = std::sqrt((q(jset,1,1,2)-qk0(jset,1,2,1))*(q(jset,1,1,2)-qk0(jset,1,2,1))+
                    (q(jset,1,1,3)-qk0(jset,1,3,1))*(q(jset,1,1,3)-qk0(jset,1,3,1))+
                    (q(jset,1,1,4)-qk0(jset,1,4,1))*(q(jset,1,1,4)-qk0(jset,1,4,1)));
                double smin_abs_v2 = smin(jset,1,1); double smin_abs2 = ccabs(smin_abs_v2);
                utau = std::sqrt((fnu(jset,1,1)+vk0(jset,1,1,1))*qset2/(smin_abs2*q(jset,1,1,1)*re));
                for (k = 1; k <= kdim-1; k++) {
                    double smin_abs_k = smin(jset,k,1); double smin_abs_kv = ccabs(smin_abs_k);
                    ypl = re*q(jset,1,1,1)*utau*smin_abs_kv/fnu(jset,1,1);
                    { double av=ypl; ypl = cclog10(av); }
                    uplus = std::sqrt(q(jset,k,1,2)*q(jset,k,1,2)+q(jset,k,1,3)*q(jset,k,1,3)+q(jset,k,1,4)*q(jset,k,1,4))/utau;
                    zkplus = turre(jset,k,1,2)/(utau*utau);
                    eplus = turre(jset,k,1,1)*re/(utau*utau*qset2/smin_abs2);
                    uvplus = vist3d(jset,k,1)*vor(jset,k,1)/(re*utau*utau*q(jset,k,1,1));
                    fprintf(f92, "%15.5e%15.5e%15.5e%15.5e%15.5e\n",
                        (float)uplus, (float)zkplus, (float)eplus, (float)uvplus, (float)ypl);
                }
                FILE* f91 = fortran_get_unit(91);
                nnumb = 2;
                fprintf(f91, "%5d\n", nnumb);
                fprintf(f91, "   uv\n");
                fprintf(f91, "   y\n");
                int jsets[] = {68, 93, 107, 113, 123, 129, 134};
                for (int jsi = 0; jsi < 7; jsi++) {
                    jset = jsets[jsi];
                    fprintf(f91, "%5d jset=%5d  -uv/uinf**2,y  x=%15.5e\n",
                        kdim-1, jset, (float)(0.5*(x(jset,1,1)+x(jset+1,1,1))));
                    for (k = 1; k <= kdim-1; k++) {
                        double smin_abs_kv2 = smin(jset,k,1); double smin_abs_k2 = ccabs(smin_abs_kv2);
                        fprintf(f91, "%15.5e%15.5e\n",
                            (float)(vist3d(jset,k,1)*vor(jset,k,1)/q(jset,k,1,1)/((double)reue*(double)xmach)),
                            (float)smin_abs_k2);
                    }
                }
            }
        }
    }
}

} // namespace foureqn_ns
