// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "spalart.h"
#include <cstdio>
#include "runtime/fortran_io.h"
#include "ccomplex.h"
#include "triv.h"
#include "sijrate2d.h"
#include "sijrate3d.h"
#include "mms.h"
#include "init.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

// Import ccomplex_ns functions that don't need rvalue wrappers
using ccomplex_ns::ccvmgt;
using ccomplex_ns::ccmax4;
using ccomplex_ns::ccmin;
using ccomplex_ns::ccmin4;
using ccomplex_ns::ccmin8;

using ccomplex_ns::ccabs;
using ccomplex_ns::ccsign;
using ccomplex_ns::ccsignrc;
using ccomplex_ns::cctanh;
using ccomplex_ns::ccatan2;
using ccomplex_ns::ccacos;
using ccomplex_ns::ccasin;
using ccomplex_ns::cctan;
using ccomplex_ns::cclog10;
using ccomplex_ns::ccdim;
using ccomplex_ns::ccerf;

// Local by-value wrappers to allow rvalue (literal) arguments for these functions
static inline double ccmaxcr(double a1, double a2) { return std::max(a1, a2); }
static inline double ccmax(double a1, double a2) { return std::max(a1, a2); }
static inline double ccmaxrc(double a1, double a2) { return std::max(a1, a2); }
static inline double ccmincr(double a1, double a2) { return std::min(a1, a2); }
static inline double ccminrc(double a1, double a2) { return std::min(a1, a2); }
static inline double ccatan(double z) { return std::atan(z); }

namespace spalart_ns {

void init(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, FortranArray3DRef<double> vol, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, int& nummem, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iflagprnt)
{
    init_ns::init(nbl, jdim, kdim, idim, q, qj0, qk0, qi0, tj0, tk0, ti0, vol, volj0, volk0, voli0, nummem, x, y, z, nou, bou, nbuf, ibufdim, iflagprnt);
}

void spalart(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, FortranArray3DRef<double> dtj, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vor, FortranArray3DRef<double> smin, FortranArray4DRef<double> tursav, FortranArray3DRef<double> turre, FortranArray3DRef<double> damp1, FortranArray3DRef<double> timestp, FortranArray3DRef<double> fnu, FortranArray2DRef<double> bx, FortranArray2DRef<double> cx, FortranArray2DRef<double> dx, FortranArray2DRef<double> fx, FortranArray2DRef<double> workx, FortranArray2DRef<double> by, FortranArray2DRef<double> cy, FortranArray2DRef<double> dy, FortranArray2DRef<double> fy, FortranArray2DRef<double> worky, FortranArray2DRef<double> bz, FortranArray2DRef<double> cz, FortranArray2DRef<double> dz, FortranArray2DRef<double> fz, FortranArray2DRef<double> workz, int& ntime, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, int& nbl, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray3DRef<double> blank, int& iover, double& sumn1, double& sumn2, int& negn1, int& negn2, FortranArray4DRef<double> tursav2, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iex, int& iex2, int& iex3, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray4DRef<double> ux, FortranArray4DRef<double> vx, int& nummem, FortranArray3DRef<double> fd, int& isarcnum)
{
    // COMMON block aliases
    float& cdes       = cmn_des.cdes;
    int32_t& ides     = cmn_des.ides;
    float& cddes      = cmn_des.cddes;
    float& dt         = cmn_info.dt;
    int32_t& icyc     = cmn_mgrd.icyc;
    int32_t& ncyc     = cmn_mgrd.ncyc;
    float& gamma      = cmn_fluid.gamma;
    float& cbar       = cmn_fluid2.cbar;
    float& tinf       = cmn_reyue.tinf;
    float& reue       = cmn_reyue.reue;
    float& xmach      = cmn_info.xmach;
    int32_t& isklton  = cmn_sklton.isklton;
    int32_t& i2d      = cmn_twod.i2d;
    int32_t& iexp     = cmn_zero.iexp;
    int32_t& nsubturb = cmn_turbconv.nsubturb;
    int32_t& nfreeze  = cmn_turbconv.nfreeze;
    int32_t& itaturb  = cmn_turbconv.itaturb;
    int32_t& iturbord = cmn_turbconv.iturbord;
    float& tur1cut    = cmn_turbconv.tur1cut;
    float& tur2cut    = cmn_turbconv.tur2cut;
    int32_t& ita      = cmn_unst.ita;
    int32_t& isarc2d  = cmn_curvat.isarc2d;
    float& sarccr3    = cmn_curvat.sarccr3;
    int32_t& isar     = cmn_curvat.isar;
    float& crot       = cmn_curvat.crot;
    int32_t& isarc3d  = cmn_curvat.isarc3d;
    int32_t& iexact_trunc = cmn_mms.iexact_trunc;
    int32_t& iexact_disc  = cmn_mms.iexact_disc;
    int32_t& iexact_ring  = cmn_mms.iexact_ring;
    int32_t& iaxi2planeturb = cmn_axisym.iaxi2planeturb;
    int32_t& istrongturbdis = cmn_axisym.istrongturbdis;
    int32_t& i_saneg    = cmn_sa_options.i_saneg;
    int32_t& i_sanoft2  = cmn_sa_options.i_sanoft2;
    float& sa_cw2       = cmn_sa_options.sa_cw2;
    float& sa_cw3       = cmn_sa_options.sa_cw3;
    float& sa_cv1       = cmn_sa_options.sa_cv1;
    float& sa_ct3       = cmn_sa_options.sa_ct3;
    float& sa_ct4       = cmn_sa_options.sa_ct4;
    float& sa_cb1       = cmn_sa_options.sa_cb1;
    float& sa_cb2       = cmn_sa_options.sa_cb2;
    float& sa_sigma     = cmn_sa_options.sa_sigma;
    float& sa_karman    = cmn_sa_options.sa_karman;
    int32_t& ilamlo     = cmn_lam.ilamlo;
    int32_t& ilamhi     = cmn_lam.ilamhi;
    int32_t& jlamlo     = cmn_lam.jlamlo;
    int32_t& jlamhi     = cmn_lam.jlamhi;
    int32_t& klamlo     = cmn_lam.klamlo;
    int32_t& klamhi     = cmn_lam.klamhi;
    int32_t& i_lam_forcezero = cmn_lam.i_lam_forcezero;
    int32_t* ncyc1      = cmn_info.ncyc1; // 1-based array[5]

    // Local variables
    double phi, xminn, factor, re, c2b, c2bp;
    double akarman, cb1, sigma, cb2, cw1, cw2, cw3, cv1, ct3, ct4;
    double tt, tt2;
    int nsubit, nnit, jd2, iwrite, nss;
    double sumn, negn_d;
    int negn;
    double sumn_loc;


    // isklton > 0: print info
    if (isklton > 0) {
        if (i_saneg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   Using SA-neg form");
        }
        if (i_sanoft2 == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   Using SA-noft2 form");
        }
        if (iexact_trunc != 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            if (iexact_trunc == 1) {
                std::snprintf(bou(nou(1),1), 120, "   Using exact MMS source in Spalart routine to evaluate truncation error, MS1");
            } else if (iexact_trunc == 2) {
                std::snprintf(bou(nou(1),1), 120, "   Using exact MMS source in Spalart routine to evaluate truncation error, MS2");
            } else if (iexact_trunc == 4) {
                std::snprintf(bou(nou(1),1), 120, "   Using exact MMS source in Spalart routine to evaluate truncation error, MS4");
            }
        } else {
            if (ides == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   computing turbulent viscosity using Spalart DES (cdes=%7.3f), block =%4d", (double)cdes, nbl);
            } else if (ides == 2) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   computing turbulent viscosity using Spalart DDES (cdes=%7.3f), block =%4d", (double)cdes, nbl);
            } else if (ides == 3) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   computing turbulent viscosity using Spalart MDDES (cdes=%7.3f,cddes=%7.3f), block =%4d", (double)cdes, (double)cddes, nbl);
            } else if (ides >= 4 && ides <= 6) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   computing turbulent viscosity using Spalart MuDDES (cdes=%7.3f,cddes=%7.3f), block =%4d", (double)cdes, (double)cddes, nbl);
            } else if (ides >= 7) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   computing turbulent viscosity using S-A fd on nu (cdes=%7.3f,cddes=%7.3f), block =%4d", (double)cdes, (double)cddes, nbl);
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   computing turbulent viscosity using Spalart, block =%4d", nbl);
            }
            if (isarc2d == 1 || isarc3d == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "      Rotation-curvature (SARC) terms are ON; sarccr3=%7.3f", (double)sarccr3);
            }
            if (isarc2d == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "      WARNING: curv terms are active in 2-D sense only!!!");
            }
            if (isar == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "      Dacles-Mariani rotation (SAR) correction is ON, crot=%7.3f", (double)crot);
            }
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     Freestream tur10 = %19.8e", (double)(float)cmn_ivals.tur10[0]);
            if (iexact_disc != 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                if (iexact_disc == 1) {
                    std::snprintf(bou(nou(1),1), 120, "     using exact MMS source in Spalart routine to evaluate discretization error, MS1");
                } else if (iexact_disc == 2) {
                    std::snprintf(bou(nou(1),1), 120, "     using exact MMS source in Spalart routine to evaluate discretization error, MS2");
                } else if (iexact_disc == 4) {
                    std::snprintf(bou(nou(1),1), 120, "     using exact MMS source in Spalart routine to evaluate discretization error, MS4");
                }
            }
            if (iaxi2planeturb == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     SA model ignoring i-dir");
            }
            if (istrongturbdis == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     strong conserv - diss terms");
            }
        }
        if (iturbord == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     1st order advection on RHS");
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     2nd order advection on RHS");
        }
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   SA constants:");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     cw2=%15.5f", (double)(float)sa_cw2);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     cw3=%15.5f", (double)(float)sa_cw3);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     cv1=%15.5f", (double)(float)sa_cv1);
        if (i_sanoft2 == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     ct3=0");
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     ct3=%15.5f", (double)(float)sa_ct3);
        }
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     ct4=%15.5f", (double)(float)sa_ct4);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     cb1=%15.5f", (double)(float)sa_cb1);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     cb2=%15.5f", (double)(float)sa_cb2);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     sigma=%15.5f", (double)(float)sa_sigma);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     akarman=%15.5f", (double)(float)sa_karman);
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
    if ((float)dt > 0.0f) {
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

    // Note: (10.**(-iexp) is machine zero)
    xminn = std::pow(10., (double)(-iexp+1));

    // Set number of subiterations
    nsubit = nsubturb;

    // Set factor that multiplies the N-S CFL number
    factor = 10.;
    // Overwrite with keyword value if nonzero
    if ((float)cmn_turbconv.cflturb[0] != 0.0f) {
        factor = (double)cmn_turbconv.cflturb[0];
    }

    // Timestep for turb model
    if ((float)dt < 0.0f) {
        for (int i = 1; i <= idim-1; i++)
        for (int k = 1; k <= kdim-1; k++)
        for (int j = 1; j <= jdim-1; j++) {
            timestp(j,k,i) = factor*vol(j,k,i)/dtj(j,k,i);
            double tmp = timestp(j,k,i); double c100 = 100.;
            timestp(j,k,i) = ccmincr(tmp, c100);
        }
    } else {
        // turbulence model advanced with physical time only
        for (int i = 1; i <= idim-1; i++)
        for (int k = 1; k <= kdim-1; k++)
        for (int j = 1; j <= jdim-1; j++) {
            timestp(j,k,i) = (double)dt;
        }
    }

    // Set up constants
    akarman = (double)sa_karman;
    cb1     = (double)sa_cb1;
    sigma   = (double)sa_sigma;
    cb2     = (double)sa_cb2;
    cw1     = cb1/akarman/akarman + (1.+cb2)/sigma;
    cw2     = (double)sa_cw2;
    cw3     = (double)sa_cw3;
    cv1     = (double)sa_cv1;
    if (i_sanoft2 == 1) {
        ct3 = 0.0;
    } else {
        ct3 = (double)sa_ct3;
    }
    ct4 = (double)sa_ct4;

    // Set up some other needed parameters
    jd2  = (jdim-1)/2;
    re   = (double)reue / (double)xmach;
    c2b  = (double)cbar / (double)tinf;
    c2bp = c2b + 1.0;
    nnit = nsubit;
    iwrite = 0;

    // Get laminar viscosity (divided by rho) at cell centers
    for (int i = 1; i <= idim-1; i++)
    for (int k = 1; k <= kdim-1; k++)
    for (int j = 1; j <= jdim-1; j++) {
        tt = (double)gamma * q(j,k,i,5) / q(j,k,i,1);
        fnu(j,k,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt)/q(j,k,i,1);
    }
    for (int i = 1; i <= idim-1; i++)
    for (int k = 1; k <= kdim-1; k++) {
        tt = (double)gamma * qj0(k,i,5,1) / qj0(k,i,1,1);
        fnu(0,k,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt)/qj0(k,i,1,1);
        tt = (double)gamma * qj0(k,i,5,3) / qj0(k,i,1,3);
        fnu(jdim,k,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt)/qj0(k,i,1,3);
    }
    for (int i = 1; i <= idim-1; i++)
    for (int j = 1; j <= jdim-1; j++) {
        tt = (double)gamma * qk0(j,i,5,1) / qk0(j,i,1,1);
        fnu(j,0,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt)/qk0(j,i,1,1);
        tt = (double)gamma * qk0(j,i,5,3) / qk0(j,i,1,3);
        fnu(j,kdim,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt)/qk0(j,i,1,3);
    }
    if (i2d != 1 && iaxi2planeturb != 1) {
        for (int j = 1; j <= jdim-1; j++)
        for (int k = 1; k <= kdim-1; k++) {
            tt = (double)gamma * qi0(j,k,5,1) / qi0(j,k,1,1);
            fnu(j,k,0) = c2bp*tt*std::sqrt(tt)/(c2b+tt)/qi0(j,k,1,1);
            tt = (double)gamma * qi0(j,k,5,3) / qi0(j,k,1,3);
            fnu(j,k,idim) = c2bp*tt*std::sqrt(tt)/(c2b+tt)/qi0(j,k,1,3);
        }
    }

    // If this is 1st global subiteration for time-accurate computation,
    // save tursav (at time step n)
    if ((float)dt > 0.0f && icyc == 1) {
        if (std::abs(ita) == 2) {
            if ((float)tursav2(1,1,1,1) == 0.0f) {
                for (int i = 1; i <= idim-1; i++)
                for (int k = 1; k <= kdim-1; k++)
                for (int j = 1; j <= jdim-1; j++) {
                    tursav2(j,k,i,3) = 0.;
                }
            } else {
                for (int i = 1; i <= idim-1; i++)
                for (int k = 1; k <= kdim-1; k++)
                for (int j = 1; j <= jdim-1; j++) {
                    tursav2(j,k,i,3) = tursav(j,k,i,1) - tursav2(j,k,i,1);
                }
            }
        }
        for (int i = 1; i <= idim-1; i++)
        for (int k = 1; k <= kdim-1; k++)
        for (int j = 1; j <= jdim-1; j++) {
            tursav2(j,k,i,1) = tursav(j,k,i,1);
        }
    }

    // Get TURRE values
    for (int i = 1; i <= idim-1; i++)
    for (int k = 1; k <= kdim-1; k++)
    for (int j = 1; j <= jdim-1; j++) {
        turre(j,k,i) = tursav(j,k,i,1);
    }

    // Now iterate to solve the equation
    sumn = 0.;
    negn = 0;
    for (int not_ = 1; not_ <= nnit; not_++) {

        // set up boundary conditions (they are in ghost cells everywhere)
        // (1) k=0 boundary:
        for (int i = 1; i <= idim-1; i++)
        for (int j = 1; j <= jdim-1; j++) {
            turre(j,0,i) = tk0(j,i,1,1);
        }
        // (2) k=kdim boundary:
        for (int i = 1; i <= idim-1; i++)
        for (int j = 1; j <= jdim-1; j++) {
            turre(j,kdim,i) = tk0(j,i,1,3);
        }
        // (3) j=0 boundary:
        for (int i = 1; i <= idim-1; i++)
        for (int k = 1; k <= kdim-1; k++) {
            turre(0,k,i) = tj0(k,i,1,1);
        }
        // (4) j=jdim boundary:
        for (int i = 1; i <= idim-1; i++)
        for (int k = 1; k <= kdim-1; k++) {
            turre(jdim,k,i) = tj0(k,i,1,3);
        }
        if (i2d != 1 && iaxi2planeturb != 1) {
            // (5) i=0 boundary:
            for (int k = 1; k <= kdim-1; k++)
            for (int j = 1; j <= jdim-1; j++) {
                turre(j,k,0) = ti0(j,k,1,1);
            }
            // (6) i=idim boundary:
            for (int k = 1; k <= kdim-1; k++)
            for (int j = 1; j <= jdim-1; j++) {
                turre(j,k,idim) = ti0(j,k,1,3);
            }
        }
        if (iturbord != 1) {
            // (1) k=0 boundary:
            for (int i = 1; i <= idim-1; i++)
            for (int j = 1; j <= jdim-1; j++) {
                turre(j,-1,i) = tk0(j,i,1,2);
            }
            // (2) k=kdim boundary:
            for (int i = 1; i <= idim-1; i++)
            for (int j = 1; j <= jdim-1; j++) {
                turre(j,kdim+1,i) = tk0(j,i,1,4);
            }
            // (3) j=0 boundary:
            for (int i = 1; i <= idim-1; i++)
            for (int k = 1; k <= kdim-1; k++) {
                turre(-1,k,i) = tj0(k,i,1,2);
            }
            // (4) j=jdim boundary:
            for (int i = 1; i <= idim-1; i++)
            for (int k = 1; k <= kdim-1; k++) {
                turre(jdim+1,k,i) = tj0(k,i,1,4);
            }
            if (i2d != 1 && iaxi2planeturb != 1) {
                // (5) i=0 boundary:
                for (int k = 1; k <= kdim-1; k++)
                for (int j = 1; j <= jdim-1; j++) {
                    turre(j,k,-1) = ti0(j,k,1,2);
                }
                // (6) i=idim boundary:
                for (int k = 1; k <= kdim-1; k++)
                for (int j = 1; j <= jdim-1; j++) {
                    turre(j,k,idim+1) = ti0(j,k,1,4);
                }
            }
        }

        // F_eta_eta viscous terms - Interior points
        for (int k = 2; k <= kdim-2; k++) {
            int kl = k-1;
            int ku = k+1;
            for (int i = 1; i <= idim-1; i++)
            for (int j = 1; j <= jdim-1; j++) {
                double volku = vol(j,ku,i);
                double xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                double yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                double zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                double volkl = vol(j,kl,i);
                double xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                double ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                double zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                double xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                double ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                double za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                double ttpo = xp*xa+yp*ya+zp*za;
                double ttmo = xm*xa+ym*ya+zm*za;
                double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                double cnud = -cb2*turre(j,k,i)/(sigma*re);
                double cap = ttp*cnud;
                double cam = ttm*cnud;
                double anutp = .5*(turre(j,k+1,i)+turre(j,k,i));
                double anutm = .5*(turre(j,k-1,i)+turre(j,k,i));
                double fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                double fnum_v = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                double diff_facp, diff_facm;
                if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                    double chi_p = anutp/fnup;
                    double chi_m = anutm/fnum_v;
                    diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                    diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                } else {
                    diff_facp = 1.0;
                    diff_facm = 1.0;
                }
                double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                double byy = -ccmaxcr(tmp1, 0.);
                double cyy = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                double dyy = -ccmaxcr(tmp2, 0.);
                vist3d(j,k,i) = -byy*turre(j,k-1,i) - cyy*turre(j,k,i) - dyy*turre(j,k+1,i);
            }
        }

        // K0 boundary points
        {
            int k = 1;
            int ku = std::min(2, kdim-1);
            for (int i = 1; i <= idim-1; i++)
            for (int j = 1; j <= jdim-1; j++) {
                double volku = vol(j,ku,i);
                double xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                double yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                double zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                double volkl = volk0(j,i,1);
                double xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                double ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                double zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                double xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                double ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                double za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                double ttpo = xp*xa+yp*ya+zp*za;
                double ttmo = xm*xa+ym*ya+zm*za;
                double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                double cnud = -cb2*turre(j,k,i)/(sigma*re);
                double cap = ttp*cnud;
                double cam = ttm*cnud;
                double anutp = .5*(turre(j,k+1,i)+turre(j,k,i));
                double anutm = .5*(turre(j,k-1,i)+turre(j,k,i));
                double fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                double fnum_v = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                double diff_facp, diff_facm;
                if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                    double chi_p = anutp/fnup;
                    double chi_m = anutm/fnum_v;
                    diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                    diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                } else {
                    diff_facp = 1.0;
                    diff_facm = 1.0;
                }
                double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                double byy = -ccmaxcr(tmp1, 0.);
                double cyy = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                double dyy = -ccmaxcr(tmp2, 0.);
                vist3d(j,k,i) = -byy*turre(j,k-1,i) - cyy*turre(j,k,i) - dyy*turre(j,k+1,i);
            }
        }

        // KDIM boundary points
        {
            int k = kdim-1;
            int kl = kdim-2;
            for (int i = 1; i <= idim-1; i++)
            for (int j = 1; j <= jdim-1; j++) {
                double volku = volk0(j,i,3);
                double xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                double yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                double zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                double volkl = vol(j,kl,i);
                double xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                double ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                double zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                double xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                double ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                double za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                double ttpo = xp*xa+yp*ya+zp*za;
                double ttmo = xm*xa+ym*ya+zm*za;
                double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                double cnud = -cb2*turre(j,k,i)/(sigma*re);
                double cap = ttp*cnud;
                double cam = ttm*cnud;
                double anutp = .5*(turre(j,k+1,i)+turre(j,k,i));
                double anutm = .5*(turre(j,k-1,i)+turre(j,k,i));
                double fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                double fnum_v = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                double diff_facp, diff_facm;
                if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                    double chi_p = anutp/fnup;
                    double chi_m = anutm/fnum_v;
                    diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                    diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                } else {
                    diff_facp = 1.0;
                    diff_facm = 1.0;
                }
                double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                double byy = -ccmaxcr(tmp1, 0.);
                double cyy = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                double dyy = -ccmaxcr(tmp2, 0.);
                vist3d(j,k,i) = -byy*turre(j,k-1,i) - cyy*turre(j,k,i) - dyy*turre(j,k+1,i);
            }
        }

        // Advective terms in eta
        if (iturbord == 1) {
            for (int k = 1; k <= kdim-1; k++)
            for (int i = 1; i <= idim-1; i++)
            for (int j = 1; j <= jdim-1; j++) {
                double xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                double yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                double zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                double tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                double one = 1.; double sgnu = ccsignrc(one, uu);
                double app = 0.5*(1.+sgnu);
                double apm = 0.5*(1.-sgnu);
                vist3d(j,k,i) = vist3d(j,k,i) - uu*(app*(turre(j,k,i)-turre(j,k-1,i)) + apm*(turre(j,k+1,i)-turre(j,k,i)));
            }
        } else {
            // 2nd order upwind; LHS remains 1st order everywhere
            for (int k = 1; k <= kdim-1; k++)
            for (int i = 1; i <= idim-1; i++)
            for (int j = 1; j <= jdim-1; j++) {
                double xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                double yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                double zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                double tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                double one = 1.; double sgnu = ccsignrc(one, uu);
                double app = 0.5*(1.+sgnu);
                double apm = 0.5*(1.-sgnu);
                vist3d(j,k,i) = vist3d(j,k,i)
                    -0.5*uu*app*turre(j,k-2,i)
                    +2.*uu*app*turre(j,k-1,i)
                    -1.5*uu*app*turre(j,k,i)
                    +1.5*uu*apm*turre(j,k,i)
                    -2.*uu*apm*turre(j,k+1,i)
                    +0.5*uu*apm*turre(j,k+2,i);
            }
        }

        // F_xi_xi viscous terms - interior points
        for (int j = 2; j <= jdim-2; j++) {
            int jl = j-1;
            int ju = j+1;
            for (int i = 1; i <= idim-1; i++)
            for (int k = 1; k <= kdim-1; k++) {
                double volju = vol(ju,k,i);
                double xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                double yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                double zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                double voljl = vol(jl,k,i);
                double xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                double ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                double zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                double xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                double ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                double za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                double ttpo = xp*xa+yp*ya+zp*za;
                double ttmo = xm*xa+ym*ya+zm*za;
                double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                double cnud = -cb2*turre(j,k,i)/(sigma*re);
                double cap = ttp*cnud;
                double cam = ttm*cnud;
                double anutp = .5*(turre(j+1,k,i)+turre(j,k,i));
                double anutm = .5*(turre(j-1,k,i)+turre(j,k,i));
                double fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                double fnum_v = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                double diff_facp, diff_facm;
                if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                    double chi_p = anutp/fnup;
                    double chi_m = anutm/fnum_v;
                    diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                    diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                } else {
                    diff_facp = 1.0;
                    diff_facm = 1.0;
                }
                double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                double bxx = -ccmaxcr(tmp1, 0.);
                double cxx = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                double dxx = -ccmaxcr(tmp2, 0.);
                vist3d(j,k,i) = vist3d(j,k,i) - bxx*turre(j-1,k,i) - cxx*turre(j,k,i) - dxx*turre(j+1,k,i);
            }
        }

        // J0 boundary points
        {
            int j = 1;
            int ju = std::min(2, jdim-1);
            for (int i = 1; i <= idim-1; i++)
            for (int k = 1; k <= kdim-1; k++) {
                double volju = vol(ju,k,i);
                double xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                double yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                double zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                double voljl = volj0(k,i,1);
                double xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                double ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                double zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                double xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                double ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                double za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                double ttpo = xp*xa+yp*ya+zp*za;
                double ttmo = xm*xa+ym*ya+zm*za;
                double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                double cnud = -cb2*turre(j,k,i)/(sigma*re);
                double cap = ttp*cnud;
                double cam = ttm*cnud;
                double anutp = .5*(turre(j+1,k,i)+turre(j,k,i));
                double anutm = .5*(turre(j-1,k,i)+turre(j,k,i));
                double fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                double fnum_v = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                double diff_facp, diff_facm;
                if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                    double chi_p = anutp/fnup;
                    double chi_m = anutm/fnum_v;
                    diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                    diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                } else {
                    diff_facp = 1.0;
                    diff_facm = 1.0;
                }
                double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                double bxx = -ccmaxcr(tmp1, 0.);
                double cxx = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                double dxx = -ccmaxcr(tmp2, 0.);
                vist3d(j,k,i) = vist3d(j,k,i) - bxx*turre(j-1,k,i) - cxx*turre(j,k,i) - dxx*turre(j+1,k,i);
            }
        }

        // JDIM boundary points
        {
            int j = jdim-1;
            int jl = jdim-2;
            for (int i = 1; i <= idim-1; i++)
            for (int k = 1; k <= kdim-1; k++) {
                double volju = volj0(k,i,3);
                double xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                double yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                double zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                double voljl = vol(jl,k,i);
                double xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                double ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                double zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                double xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                double ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                double za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                double ttpo = xp*xa+yp*ya+zp*za;
                double ttmo = xm*xa+ym*ya+zm*za;
                double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                double cnud = -cb2*turre(j,k,i)/(sigma*re);
                double cap = ttp*cnud;
                double cam = ttm*cnud;
                double anutp = .5*(turre(j+1,k,i)+turre(j,k,i));
                double anutm = .5*(turre(j-1,k,i)+turre(j,k,i));
                double fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                double fnum_v = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                double diff_facp, diff_facm;
                if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                    double chi_p = anutp/fnup;
                    double chi_m = anutm/fnum_v;
                    diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                    diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                } else {
                    diff_facp = 1.0;
                    diff_facm = 1.0;
                }
                double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                double bxx = -ccmaxcr(tmp1, 0.);
                double cxx = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                double dxx = -ccmaxcr(tmp2, 0.);
                vist3d(j,k,i) = vist3d(j,k,i) - bxx*turre(j-1,k,i) - cxx*turre(j,k,i) - dxx*turre(j+1,k,i);
            }
        }

        // Advective terms in xi
        if (iturbord == 1) {
            for (int i = 1; i <= idim-1; i++)
            for (int k = 1; k <= kdim-1; k++)
            for (int j = 1; j <= jdim-1; j++) {
                double xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                double yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                double zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                double tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                double one = 1.; double sgnu = ccsignrc(one, uu);
                double app = 0.5*(1.+sgnu);
                double apm = 0.5*(1.-sgnu);
                vist3d(j,k,i) = vist3d(j,k,i) - uu*(app*(turre(j,k,i)-turre(j-1,k,i)) + apm*(turre(j+1,k,i)-turre(j,k,i)));
            }
        } else {
            // 2nd order upwind
            for (int i = 1; i <= idim-1; i++)
            for (int k = 1; k <= kdim-1; k++)
            for (int j = 1; j <= jdim-1; j++) {
                double xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                double yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                double zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                double tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                double one = 1.; double sgnu = ccsignrc(one, uu);
                double app = 0.5*(1.+sgnu);
                double apm = 0.5*(1.-sgnu);
                vist3d(j,k,i) = vist3d(j,k,i)
                    -0.5*uu*app*turre(j-2,k,i)
                    +2.*uu*app*turre(j-1,k,i)
                    -1.5*uu*app*turre(j,k,i)
                    +1.5*uu*apm*turre(j,k,i)
                    -2.*uu*apm*turre(j+1,k,i)
                    +0.5*uu*apm*turre(j+2,k,i);
            }
        }

        // F_zeta_zeta viscous terms
        if (i2d != 1 && iaxi2planeturb != 1) {
            // Interior points
            for (int i = 2; i <= idim-2; i++) {
                int il = i-1;
                int iu = i+1;
                for (int k = 1; k <= kdim-1; k++)
                for (int j = 1; j <= jdim-1; j++) {
                    double voliu = vol(j,k,iu);
                    double xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    double yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    double zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    double volil = vol(j,k,il);
                    double xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    double ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    double zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    double xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    double ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    double za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    double ttpo = xp*xa+yp*ya+zp*za;
                    double ttmo = xm*xa+ym*ya+zm*za;
                    double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    double cnud = -cb2*turre(j,k,i)/(sigma*re);
                    double cap = ttp*cnud;
                    double cam = ttm*cnud;
                    double anutp = .5*(turre(j,k,i+1)+turre(j,k,i));
                    double anutm = .5*(turre(j,k,i-1)+turre(j,k,i));
                    double fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                    double fnum_v = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                    double diff_facp, diff_facm;
                    if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                        double chi_p = anutp/fnup;
                        double chi_m = anutm/fnum_v;
                        diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                        diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                    } else {
                        diff_facp = 1.0;
                        diff_facm = 1.0;
                    }
                    double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                    double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                    double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                    double bzz = -ccmaxcr(tmp1, 0.);
                    double czz = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                    double dzz = -ccmaxcr(tmp2, 0.);
                    vist3d(j,k,i) = vist3d(j,k,i) - bzz*turre(j,k,i-1) - czz*turre(j,k,i) - dzz*turre(j,k,i+1);
                }
            }

            // I0 boundary points
            {
                int i = 1;
                int iu = std::min(2, idim-1);
                for (int k = 1; k <= kdim-1; k++)
                for (int j = 1; j <= jdim-1; j++) {
                    double voliu = vol(j,k,iu);
                    double xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    double yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    double zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    double volil = voli0(j,k,1);
                    double xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    double ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    double zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    double xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    double ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    double za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    double ttpo = xp*xa+yp*ya+zp*za;
                    double ttmo = xm*xa+ym*ya+zm*za;
                    double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    double cnud = -cb2*turre(j,k,i)/(sigma*re);
                    double cap = ttp*cnud;
                    double cam = ttm*cnud;
                    double anutp = .5*(turre(j,k,i+1)+turre(j,k,i));
                    double anutm = .5*(turre(j,k,i-1)+turre(j,k,i));
                    double fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                    double fnum_v = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                    double diff_facp, diff_facm;
                    if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                        double chi_p = anutp/fnup;
                        double chi_m = anutm/fnum_v;
                        diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                        diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                    } else {
                        diff_facp = 1.0;
                        diff_facm = 1.0;
                    }
                    double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                    double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                    double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                    double bzz = -ccmaxcr(tmp1, 0.);
                    double czz = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                    double dzz = -ccmaxcr(tmp2, 0.);
                    vist3d(j,k,i) = vist3d(j,k,i) - bzz*turre(j,k,i-1) - czz*turre(j,k,i) - dzz*turre(j,k,i+1);
                }
            }

            // IDIM boundary points
            {
                int i = idim-1;
                int il = idim-2;
                for (int k = 1; k <= kdim-1; k++)
                for (int j = 1; j <= jdim-1; j++) {
                    double voliu = voli0(j,k,3);
                    double xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    double yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    double zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    double volil = vol(j,k,il);
                    double xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    double ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    double zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    double xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    double ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    double za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    double ttpo = xp*xa+yp*ya+zp*za;
                    double ttmo = xm*xa+ym*ya+zm*za;
                    double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    double cnud = -cb2*turre(j,k,i)/(sigma*re);
                    double cap = ttp*cnud;
                    double cam = ttm*cnud;
                    double anutp = .5*(turre(j,k,i+1)+turre(j,k,i));
                    double anutm = .5*(turre(j,k,i-1)+turre(j,k,i));
                    double fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                    double fnum_v = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                    double diff_facp, diff_facm;
                    if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                        double chi_p = anutp/fnup;
                        double chi_m = anutm/fnum_v;
                        diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                        diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                    } else {
                        diff_facp = 1.0;
                        diff_facm = 1.0;
                    }
                    double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                    double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                    double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                    double bzz = -ccmaxcr(tmp1, 0.);
                    double czz = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                    double dzz = -ccmaxcr(tmp2, 0.);
                    vist3d(j,k,i) = vist3d(j,k,i) - bzz*turre(j,k,i-1) - czz*turre(j,k,i) - dzz*turre(j,k,i+1);
                }
            }

            // Advective terms in zeta
            if (iturbord == 1) {
                for (int i = 1; i <= idim-1; i++)
                for (int k = 1; k <= kdim-1; k++)
                for (int j = 1; j <= jdim-1; j++) {
                    double xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    double yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    double zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    double tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    double one = 1.; double sgnu = ccsignrc(one, uu);
                    double app = 0.5*(1.+sgnu);
                    double apm = 0.5*(1.-sgnu);
                    vist3d(j,k,i) = vist3d(j,k,i) - uu*(app*(turre(j,k,i)-turre(j,k,i-1)) + apm*(turre(j,k,i+1)-turre(j,k,i)));
                }
            } else {
                // 2nd order upwind
                for (int i = 1; i <= idim-1; i++)
                for (int k = 1; k <= kdim-1; k++)
                for (int j = 1; j <= jdim-1; j++) {
                    double xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    double yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    double zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    double tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    double one = 1.; double sgnu = ccsignrc(one, uu);
                    double app = 0.5*(1.+sgnu);
                    double apm = 0.5*(1.-sgnu);
                    vist3d(j,k,i) = vist3d(j,k,i)
                        -0.5*uu*app*turre(j,k,i-2)
                        +2.*uu*app*turre(j,k,i-1)
                        -1.5*uu*app*turre(j,k,i)
                        +1.5*uu*apm*turre(j,k,i)
                        -2.*uu*apm*turre(j,k,i+1)
                        +0.5*uu*apm*turre(j,k,i+2);
                }
            }
        } // end if i2d != 1 && iaxi2planeturb != 1 (zeta terms)

        // Add source term to RHS
        if (isklton > 0) {
            if (ilamlo == 0 || jlamlo == 0 || klamlo == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   block %4d in S-A turb model has no laminar regions", nbl);
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   block %4d in S-A turb model - laminar region is:", nbl);
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

        // can use damp1 to temporarily store distance variable (different for DES)
        if (ides == 1) {
            // DES
            for (int i = 1; i <= idim-1; i++)
            for (int j = 1; j <= jdim-1; j++)
            for (int k = 1; k <= kdim-1; k++) {
                double deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4));
                double deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4));
                double deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4));
                double delta_v = deltaj; double tmp2 = deltak;
                delta_v = ccmax(delta_v, tmp2);
                if (i2d != 1 && iaxi2planeturb != 1) {
                    delta_v = ccmax(delta_v, deltai);
                }
                double tmp_abs = ccabs(smin(j,k,i));
                double tmp_cdes = (double)cdes * delta_v;
                damp1(j,k,i) = ccmin(tmp_abs, tmp_cdes);
            }
        } else if (ides >= 2) {
            // DDES
            for (int i = 1; i <= idim-1; i++)
            for (int j = 1; j <= jdim-1; j++)
            for (int k = 1; k <= kdim-1; k++) {
                double deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4));
                double deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4));
                double deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4));
                double delta_v = deltaj; double tmp2 = deltak;
                delta_v = ccmax(delta_v, tmp2);
                if (i2d != 1 && iaxi2planeturb != 1) {
                    delta_v = ccmax(delta_v, deltai);
                }
                double dist_v = ccabs(smin(j,k,i));
                double velterm = ux(j,k,i,1)*ux(j,k,i,1) + ux(j,k,i,2)*ux(j,k,i,2) + ux(j,k,i,3)*ux(j,k,i,3)
                               + ux(j,k,i,4)*ux(j,k,i,4) + ux(j,k,i,5)*ux(j,k,i,5) + ux(j,k,i,6)*ux(j,k,i,6)
                               + ux(j,k,i,7)*ux(j,k,i,7) + ux(j,k,i,8)*ux(j,k,i,8) + ux(j,k,i,9)*ux(j,k,i,9);
                double rd = turre(j,k,i)/(std::sqrt(velterm)*akarman*akarman*dist_v*dist_v*re);
                double arg8rd = (8.0*rd)*(8.0*rd)*(8.0*rd);
                fd(j,k,i) = 1.0 - cctanh(arg8rd);
                if (ides != 7) {
                    double tmp_zero = 0.0;
                    double term_v = ccmaxrc(tmp_zero, dist_v - (double)cdes*delta_v);
                    damp1(j,k,i) = dist_v - fd(j,k,i)*term_v;
                } else {
                    damp1(j,k,i) = ccabs(smin(j,k,i));
                }
            }
        } else {
            for (int i = 1; i <= idim-1; i++)
            for (int j = 1; j <= jdim-1; j++)
            for (int k = 1; k <= kdim-1; k++) {
                damp1(j,k,i) = ccabs(smin(j,k,i));
            }
        }

        // For SARC model, compute rtilde as fn of Dalpha/Dt
        if (isarc2d == 1) {
            sijrate2d_ns::sijrate2d(idim, jdim, kdim, q, qj0, qk0, bcj, bck, vol, sj, sk, vx);
            // compute curvature term, part of Dalpha/Dt: store in vx(1)
            for (int i = 1; i <= idim-1; i++)
            for (int j = 1; j <= jdim-1; j++)
            for (int k = 1; k <= kdim-1; k++) {
                double s11 = ux(j,k,i,1) - (ux(j,k,i,1)+ux(j,k,i,9))/2.;
                double s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7));
                vx(j,k,i,1) = s11*vx(j,k,i,4) - s13*vx(j,k,i,3);
            }
        }
        if (isarc3d == 1) {
            sijrate3d_ns::sijrate3d(idim, jdim, kdim, q, ux, vol, si, sj, sk, vx);
        }

        // Main source term loop
        for (int i = 1; i <= idim-1; i++)
        for (int j = 1; j <= jdim-1; j++)
        for (int k = 1; k <= kdim-1; k++) {
            double cutoff;
            if ((i >= ilamlo && i < ilamhi &&
                 j >= jlamlo && j < jlamhi &&
                 k >= klamlo && k < klamhi) ||
                (float)smin(j,k,i) < 0.0f) {
                cutoff = 0.;
            } else if (ides == 3 && (float)fd(j,k,i) > (float)cddes) {
                cutoff = (1.0 - fd(j,k,i))/(1.0-(double)cddes);
            } else if (ides == 4 && (float)fd(j,k,i) > (float)cddes) {
                cutoff = (1.0 - fd(j,k,i))/(1.0-(double)cddes);
            } else if (ides == 5 && (float)fd(j,k,i) > 0.999f) {
                cutoff = (1.0 - fd(j,k,i))/(1.0-0.999);
            } else {
                cutoff = 1.;
            }
            double ss = vor(j,k,i);
            if (isar == 1) {
                double s11 = ux(j,k,i,1);
                double s22 = ux(j,k,i,5);
                double s33 = ux(j,k,i,9);
                double s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4));
                double s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7));
                double s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8));
                double xis = s11*s11 + s22*s22 + s33*s33 + 2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23;
                double xisabs = std::sqrt(2.*xis);
                double tmp_zero = 0.; double tmp_diff = xisabs - vor(j,k,i);
                ss = vor(j,k,i) + (double)crot * ccminrc(tmp_zero, tmp_diff);
            }
            double chi = turre(j,k,i)/fnu(j,k,i);
            double fv1, fv2, sst, rr, gg, fw, ft2;
            double sbar;
            if (i_saneg != 1) {
                fv1 = chi*chi*chi/(chi*chi*chi+cv1*cv1*cv1);
                fv2 = 1.-(chi/(1.+chi*fv1));
                sst = ss + turre(j,k,i)*fv2/(re*(akarman*damp1(j,k,i))*(akarman*damp1(j,k,i)));
                sst = ccmax(sst, xminn);
                rr = turre(j,k,i)/(sst*re*(akarman*damp1(j,k,i))*(akarman*damp1(j,k,i)));
                rr = ccmincr(rr, 10.);
                gg = rr + cw2*(std::pow(rr,6.)-rr);
                gg = ccmax(gg, xminn);
                fw = gg*std::pow((1.+std::pow(cw3,6.))/(std::pow(gg,6.)+std::pow(cw3,6.)), 1./6.);
                ft2 = ct3*std::exp(-ct4*chi*chi);
            } else if (i_saneg == 1 && (float)turre(j,k,i) >= 0.0f) {
                fv1 = chi*chi*chi/(chi*chi*chi+cv1*cv1*cv1);
                fv2 = 1.-(chi/(1.+chi*fv1));
                sbar = turre(j,k,i)*fv2/(re*(akarman*damp1(j,k,i))*(akarman*damp1(j,k,i)));
                if ((float)sbar >= -0.7f*(float)ss) {
                    sst = ss + sbar;
                } else {
                    sst = ss + (ss*(0.7*0.7*ss+0.9*sbar)/((0.9-2.*0.7)*ss-sbar));
                }
                if (sst == 0.) {
                    rr = 10.;
                } else {
                    rr = turre(j,k,i)/(sst*re*(akarman*damp1(j,k,i))*(akarman*damp1(j,k,i)));
                    rr = ccmincr(rr, 10.);
                }
                gg = rr + cw2*(std::pow(rr,6.)-rr);
                gg = ccmax(gg, xminn);
                fw = gg*std::pow((1.+std::pow(cw3,6.))/(std::pow(gg,6.)+std::pow(cw3,6.)), 1./6.);
                ft2 = ct3*std::exp(-ct4*chi*chi);
            } else {
                // i_saneg==1 and turre<0: set dummy values (used below only in else branch)
                fv1 = 0.; fv2 = 0.; sst = 0.; rr = 0.; gg = xminn; fw = 0.; ft2 = 0.; sbar = 0.;
            }

            double term1, term2;
            if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                term1 = cb1*(1.-ct3)*ss;
                term2 = cw1;
            } else {
                double tmp_zero = 0.; double tmp_ss = ss;
                term1 = cb1*(1.-ft2)*ccmaxrc(tmp_zero, tmp_ss);
                term2 = cb1*((1.-ft2)*fv2+ft2)/akarman/akarman - cw1*fw;
                // For SARC
                if (isarc2d == 1 || isarc3d == 1) {
                    double s11 = ux(j,k,i,1);
                    double s22 = ux(j,k,i,5);
                    double s33 = ux(j,k,i,9);
                    double s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4));
                    double s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7));
                    double s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8));
                    double w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4));
                    double w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7));
                    double w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8));
                    double xis = s11*s11 + s22*s22 + s33*s33 + 2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23;
                    double ss_tmp = ss; ss_tmp = ccmax(ss_tmp, xminn);
                    ss = ss_tmp;
                    double xisabs = std::sqrt(2.*xis);
                    double rstar = xisabs/ss;
                    xisabs = ccmax(xisabs, xminn);
                    double rtilde;
                    if (isarc2d == 1) {
                        rtilde = -4.*vx(j,k,i,1)*w13 / std::pow(0.5*(ss*ss+xisabs*xisabs), 2.);
                    } else {
                        rtilde = 2./std::pow(0.5*(ss*ss+xisabs*xisabs), 2.) *
                            ( -w12*vx(j,k,i,2)*(s11-s22)
                              -w13*vx(j,k,i,3)*(s11-s33)
                              -w23*vx(j,k,i,5)*(s22-s33)
                              +s12*(-w12*(vx(j,k,i,4)-vx(j,k,i,1)) - w13*vx(j,k,i,5) - w23*vx(j,k,i,3))
                              +s13*(-w13*(vx(j,k,i,6)-vx(j,k,i,1)) - w12*vx(j,k,i,5) + w23*vx(j,k,i,2))
                              +s23*(-w23*(vx(j,k,i,6)-vx(j,k,i,4)) + w12*vx(j,k,i,3) + w13*vx(j,k,i,2)) );
                    }
                    double fr1 = 4.*rstar/(1.+rstar)*(1.-(double)sarccr3*ccatan(12.*rtilde)) - 1.;
                    term1 = cb1*(fr1-ft2)*ss;
                    term2 = cb1*((fr1-ft2)*fv2+ft2)/akarman/akarman - cw1*fw;
                }
            }
            double dist2i = 1./(re*damp1(j,k,i)*damp1(j,k,i)+1.e-20);
            double tt_src = cutoff*term1*turre(j,k,i) + term2*turre(j,k,i)*turre(j,k,i)*dist2i;
            // Store quantity to be added to certain implicit LHS terms:
            damp1(j,k,i) = 2.*term2*turre(j,k,i)*dist2i;
            double dfv1 = (fv1-fv1*fv1)*3./turre(j,k,i);
            double dfv2 = (fv2-1.)/turre(j,k,i) + ((1.-fv2)*(1.-fv2))*(fv1/turre(j,k,i)+dfv1);
            double dft2 = -(2.*ct4*turre(j,k,i)/(fnu(j,k,i)*fnu(j,k,i)))*ft2;
            double drr = rr/turre(j,k,i) - rr*rr*(fv2/turre(j,k,i)+dfv2);
            double dgg = (1.-cw2+6.*cw2*std::pow(rr,5.))*drr;
            gg = ccmax(gg, xminn*10.);
            // fix for single prec. via TLNS3D via FUN3D via Tim Barth
            double dfw = std::pow((1.+std::pow(cw3,6.))/(std::pow(gg,6.)+std::pow(cw3,6.)), 1./6.)
                       - (std::pow(1.+std::pow(cw3,6.), 1./6.) / std::pow(std::pow(gg,6.)+std::pow(cw3,6.), 7./6.)) * std::pow(gg,6.);
            dfw = dfw*dgg;
            if (!(i_saneg == 1 && (float)turre(j,k,i) < 0.0f)) {
                damp1(j,k,i) = damp1(j,k,i) +
                    dist2i*(turre(j,k,i)*turre(j,k,i))*(cb1/(akarman*akarman)*
                    (dfv2-ft2*dfv2-fv2*dft2+dft2) - cw1*dfw);
            }
            // Add to RHS:
            vist3d(j,k,i) = vist3d(j,k,i) + tt_src;
        } // end main source term loop

        if (iexact_trunc != 0 || iexact_disc != 0) {
            // note vol already included in resid and sign is reverse of mean flow
            mms_ns::exact_turb_force(jdim, kdim, idim, x, y, z, vol, vist3d, smin, iexact_trunc, iexact_disc);
            if (icyc == ncyc1[0] || icyc == ncyc1[1] || icyc == ncyc1[2] ||
                icyc == ncyc1[3] || icyc == ncyc1[4]) {
                if (ntime == 1) {
                    mms_ns::exact_turb_res(jdim, kdim, idim, vol, vist3d, iexact_trunc, iexact_disc);
                }
            }
        }

        // Implicit F_eta_eta viscous terms. Do over all i's
        for (int i = 1; i <= idim-1; i++) {
            // Interior points
            for (int k = 2; k <= kdim-2; k++) {
                int kl = k-1;
                int ku = k+1;
                for (int j = 1; j <= jdim-1; j++) {
                    double volku = vol(j,ku,i);
                    double xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    double yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    double zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    double volkl = vol(j,kl,i);
                    double xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    double ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    double zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    double xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    double ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    double za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    double ttpo = xp*xa+yp*ya+zp*za;
                    double ttmo = xm*xa+ym*ya+zm*za;
                    double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                    double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                    double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    double cnud = -cb2*turre(j,k,i)/(sigma*re);
                    double cap = ttp*cnud;
                    double cam = ttm*cnud;
                    double anutp = .5*(turre(j,k+1,i)+turre(j,k,i));
                    double anutm = .5*(turre(j,k-1,i)+turre(j,k,i));
                    double fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                    double fnum_v = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                    double diff_facp, diff_facm;
                    if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                        double chi_p = anutp/fnup;
                        double chi_m = anutm/fnum_v;
                        diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                        diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                    } else {
                        diff_facp = 1.0;
                        diff_facm = 1.0;
                    }
                    double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                    double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                    double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                    by(j,k) = -ccmaxcr(tmp1, 0.);
                    cy(j,k) = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                    dy(j,k) = -ccmaxcr(tmp2, 0.);
                }
                for (int j = 1; j <= jdim-1; j++) {
                    double xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    double yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    double zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    double tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                    double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    double one = 1.; double sgnu = ccsignrc(one, uu);
                    double app = 0.5*(1.+sgnu);
                    double apm = 0.5*(1.-sgnu);
                    by(j,k) = by(j,k) - uu*app;
                    if ((float)damp1(j,k,i) < 0.0f) {
                        cy(j,k) = cy(j,k) + uu*(app-apm) - damp1(j,k,i);
                    } else {
                        cy(j,k) = cy(j,k) + uu*(app-apm);
                    }
                    dy(j,k) = dy(j,k) + uu*apm;
                }
                for (int j = 1; j <= jdim-1; j++) {
                    double fact = timestp(j,k,i);
                    by(j,k) = by(j,k)*fact;
                    cy(j,k) = cy(j,k)*fact + 1.0*(1.+phi);
                    dy(j,k) = dy(j,k)*fact;
                    fy(j,k) = vist3d(j,k,i)*fact;
                }
                if ((float)dt > 0.0f) {
                    for (int j = 1; j <= jdim-1; j++) {
                        fy(j,k) = fy(j,k) + (1.+phi)*(tursav2(j,k,i,1)-turre(j,k,i)) + phi*tursav2(j,k,i,3);
                    }
                }
            } // end interior k loop

            // K0 boundary points
            {
                int k = 1;
                int ku = std::min(2, kdim-1);
                for (int j = 1; j <= jdim-1; j++) {
                    double volku = vol(j,ku,i);
                    double xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    double yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    double zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    double volkl = volk0(j,i,1);
                    double xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    double ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    double zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    double xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    double ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    double za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    double ttpo = xp*xa+yp*ya+zp*za;
                    double ttmo = xm*xa+ym*ya+zm*za;
                    double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                    double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                    double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    double cnud = -cb2*turre(j,k,i)/(sigma*re);
                    double cap = ttp*cnud;
                    double cam = ttm*cnud;
                    double anutp = .5*(turre(j,k+1,i)+turre(j,k,i));
                    double anutm = .5*(turre(j,k-1,i)+turre(j,k,i));
                    double fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                    double fnum_v = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                    double diff_facp, diff_facm;
                    if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                        double chi_p = anutp/fnup;
                        double chi_m = anutm/fnum_v;
                        diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                        diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                    } else {
                        diff_facp = 1.0;
                        diff_facm = 1.0;
                    }
                    double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                    double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                    double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                    by(j,k) = -ccmaxcr(tmp1, 0.);
                    cy(j,k) = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                    dy(j,k) = -ccmaxcr(tmp2, 0.);
                }
                for (int j = 1; j <= jdim-1; j++) {
                    double xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    double yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    double zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    double tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                    double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    double one = 1.; double sgnu = ccsignrc(one, uu);
                    double app = 0.5*(1.+sgnu);
                    double apm = 0.5*(1.-sgnu);
                    by(j,k) = by(j,k) - uu*app;
                    if ((float)damp1(j,k,i) < 0.0f) {
                        cy(j,k) = cy(j,k) + uu*(app-apm) - damp1(j,k,i);
                    } else {
                        cy(j,k) = cy(j,k) + uu*(app-apm);
                    }
                    dy(j,k) = dy(j,k) + uu*apm;
                }
                for (int j = 1; j <= jdim-1; j++) {
                    double fact = timestp(j,k,i);
                    by(j,k) = by(j,k)*fact;
                    cy(j,k) = cy(j,k)*fact + 1.0*(1.+phi);
                    dy(j,k) = dy(j,k)*fact;
                    fy(j,k) = vist3d(j,k,i)*fact;
                }
                if ((float)dt > 0.0f) {
                    for (int j = 1; j <= jdim-1; j++) {
                        fy(j,k) = fy(j,k) + (1.+phi)*(tursav2(j,k,i,1)-turre(j,k,i)) + phi*tursav2(j,k,i,3);
                    }
                }
            } // end K0 boundary

            // KDIM boundary points
            {
                int k = kdim-1;
                int kl = kdim-2;
                for (int j = 1; j <= jdim-1; j++) {
                    double volku = volk0(j,i,3);
                    double xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    double yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    double zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    double volkl = vol(j,kl,i);
                    double xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    double ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    double zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    double xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    double ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    double za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    double ttpo = xp*xa+yp*ya+zp*za;
                    double ttmo = xm*xa+ym*ya+zm*za;
                    double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                    double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                    double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    double cnud = -cb2*turre(j,k,i)/(sigma*re);
                    double cap = ttp*cnud;
                    double cam = ttm*cnud;
                    double anutp = .5*(turre(j,k+1,i)+turre(j,k,i));
                    double anutm = .5*(turre(j,k-1,i)+turre(j,k,i));
                    double fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                    double fnum_v = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                    double diff_facp, diff_facm;
                    if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                        double chi_p = anutp/fnup;
                        double chi_m = anutm/fnum_v;
                        diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                        diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                    } else {
                        diff_facp = 1.0;
                        diff_facm = 1.0;
                    }
                    double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                    double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                    double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                    by(j,k) = -ccmaxcr(tmp1, 0.);
                    cy(j,k) = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                    dy(j,k) = -ccmaxcr(tmp2, 0.);
                }
                for (int j = 1; j <= jdim-1; j++) {
                    double xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    double yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    double zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    double tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                    double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    double one = 1.; double sgnu = ccsignrc(one, uu);
                    double app = 0.5*(1.+sgnu);
                    double apm = 0.5*(1.-sgnu);
                    by(j,k) = by(j,k) - uu*app;
                    if ((float)damp1(j,k,i) < 0.0f) {
                        cy(j,k) = cy(j,k) + uu*(app-apm) - damp1(j,k,i);
                    } else {
                        cy(j,k) = cy(j,k) + uu*(app-apm);
                    }
                    dy(j,k) = dy(j,k) + uu*apm;
                }
                for (int j = 1; j <= jdim-1; j++) {
                    double fact = timestp(j,k,i);
                    by(j,k) = by(j,k)*fact;
                    cy(j,k) = cy(j,k)*fact + 1.0*(1.+phi);
                    dy(j,k) = dy(j,k)*fact;
                    fy(j,k) = vist3d(j,k,i)*fact;
                }
                if ((float)dt > 0.0f) {
                    for (int j = 1; j <= jdim-1; j++) {
                        fy(j,k) = fy(j,k) + (1.+phi)*(tursav2(j,k,i,1)-turre(j,k,i)) + phi*tursav2(j,k,i,3);
                    }
                }
            } // end KDIM boundary
            if (iover == 1) {
                for (int k = 1; k <= kdim-1; k++)
                for (int j = 1; j <= jdim-1; j++) {
                    fy(j,k) = fy(j,k)*blank(j,k,i);
                    by(j,k) = by(j,k)*blank(j,k,i);
                    dy(j,k) = dy(j,k)*blank(j,k,i);
                    cy(j,k) = cy(j,k)*blank(j,k,i) + (1.-blank(j,k,i));
                }
            }
            {
                int jl2=1, ju2=jdim-1, kl2=1, ku2=kdim-1;
                int jdim1=jdim-1, kdim1=kdim-1;
                triv_ns::triv(jdim1, kdim1, jl2, ju2, kl2, ku2, worky, by, cy, dy, fy);
            }
            for (int k = 1; k <= kdim-1; k++)
            for (int j = 1; j <= jdim-1; j++) {
                vist3d(j,k,i) = fy(j,k);
            }
        } // end implicit eta loop over i

        // Implicit F_xi_xi viscous terms. Do over all i's
        for (int i = 1; i <= idim-1; i++) {
            // Interior points
            for (int j = 2; j <= jdim-2; j++) {
                int jl = j-1;
                int ju = j+1;
                for (int k = 1; k <= kdim-1; k++) {
                    double volju = vol(ju,k,i);
                    double xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    double yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    double zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    double voljl = vol(jl,k,i);
                    double xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    double ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    double zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    double xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    double ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    double za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    double ttpo = xp*xa+yp*ya+zp*za;
                    double ttmo = xm*xa+ym*ya+zm*za;
                    double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                    double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                    double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    double cnud = -cb2*turre(j,k,i)/(sigma*re);
                    double cap = ttp*cnud;
                    double cam = ttm*cnud;
                    double anutp = .5*(turre(j+1,k,i)+turre(j,k,i));
                    double anutm = .5*(turre(j-1,k,i)+turre(j,k,i));
                    double fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                    double fnum_v = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                    double diff_facp, diff_facm;
                    if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                        double chi_p = anutp/fnup;
                        double chi_m = anutm/fnum_v;
                        diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                        diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                    } else {
                        diff_facp = 1.0;
                        diff_facm = 1.0;
                    }
                    double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                    double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                    double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                    bx(k,j) = -ccmaxcr(tmp1, 0.);
                    cx(k,j) = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                    dx(k,j) = -ccmaxcr(tmp2, 0.);
                }
                for (int k = 1; k <= kdim-1; k++) {
                    double xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    double yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    double zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    double tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                    double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    double one = 1.; double sgnu = ccsignrc(one, uu);
                    double app = 0.5*(1.+sgnu);
                    double apm = 0.5*(1.-sgnu);
                    bx(k,j) = bx(k,j) - uu*app;
                    cx(k,j) = cx(k,j) + uu*(app-apm);
                    dx(k,j) = dx(k,j) + uu*apm;
                }
                for (int k = 1; k <= kdim-1; k++) {
                    double fact = timestp(j,k,i);
                    bx(k,j) = bx(k,j)*fact;
                    cx(k,j) = cx(k,j)*fact + 1.0*(1.+phi);
                    dx(k,j) = dx(k,j)*fact;
                    fx(k,j) = vist3d(j,k,i)*(1.+phi);
                }
            } // end interior j loop

            // J0 boundary points
            {
                int j = 1;
                int ju = std::min(2, jdim-1);
                for (int k = 1; k <= kdim-1; k++) {
                    double volju = vol(ju,k,i);
                    double xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    double yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    double zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    double voljl = volj0(k,i,1);
                    double xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    double ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    double zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    double xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    double ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    double za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    double ttpo = xp*xa+yp*ya+zp*za;
                    double ttmo = xm*xa+ym*ya+zm*za;
                    double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                    double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                    double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    double cnud = -cb2*turre(j,k,i)/(sigma*re);
                    double cap = ttp*cnud;
                    double cam = ttm*cnud;
                    double anutp = .5*(turre(j+1,k,i)+turre(j,k,i));
                    double anutm = .5*(turre(j-1,k,i)+turre(j,k,i));
                    double fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                    double fnum_v = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                    double diff_facp, diff_facm;
                    if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                        double chi_p = anutp/fnup;
                        double chi_m = anutm/fnum_v;
                        diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                        diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                    } else {
                        diff_facp = 1.0;
                        diff_facm = 1.0;
                    }
                    double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                    double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                    double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                    bx(k,j) = -ccmaxcr(tmp1, 0.);
                    cx(k,j) = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                    dx(k,j) = -ccmaxcr(tmp2, 0.);
                }
                for (int k = 1; k <= kdim-1; k++) {
                    double xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    double yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    double zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    double tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                    double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    double one = 1.; double sgnu = ccsignrc(one, uu);
                    double app = 0.5*(1.+sgnu);
                    double apm = 0.5*(1.-sgnu);
                    bx(k,j) = bx(k,j) - uu*app;
                    cx(k,j) = cx(k,j) + uu*(app-apm);
                    dx(k,j) = dx(k,j) + uu*apm;
                }
                for (int k = 1; k <= kdim-1; k++) {
                    double fact = timestp(j,k,i);
                    bx(k,j) = bx(k,j)*fact;
                    cx(k,j) = cx(k,j)*fact + 1.0*(1.+phi);
                    dx(k,j) = dx(k,j)*fact;
                    fx(k,j) = vist3d(j,k,i)*(1.+phi);
                }
            } // end J0 boundary

            // JDIM boundary points
            {
                int j = jdim-1;
                int jl = jdim-2;
                for (int k = 1; k <= kdim-1; k++) {
                    double volju = volj0(k,i,3);
                    double xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    double yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    double zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    double voljl = vol(jl,k,i);
                    double xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    double ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    double zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    double xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    double ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    double za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    double ttpo = xp*xa+yp*ya+zp*za;
                    double ttmo = xm*xa+ym*ya+zm*za;
                    double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                    double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                    double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    double cnud = -cb2*turre(j,k,i)/(sigma*re);
                    double cap = ttp*cnud;
                    double cam = ttm*cnud;
                    double anutp = .5*(turre(j+1,k,i)+turre(j,k,i));
                    double anutm = .5*(turre(j-1,k,i)+turre(j,k,i));
                    double fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                    double fnum_v = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                    double diff_facp, diff_facm;
                    if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                        double chi_p = anutp/fnup;
                        double chi_m = anutm/fnum_v;
                        diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                        diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                    } else {
                        diff_facp = 1.0;
                        diff_facm = 1.0;
                    }
                    double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                    double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                    double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                    bx(k,j) = -ccmaxcr(tmp1, 0.);
                    cx(k,j) = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                    dx(k,j) = -ccmaxcr(tmp2, 0.);
                }
                for (int k = 1; k <= kdim-1; k++) {
                    double xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    double yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    double zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    double tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                    double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    double one = 1.; double sgnu = ccsignrc(one, uu);
                    double app = 0.5*(1.+sgnu);
                    double apm = 0.5*(1.-sgnu);
                    bx(k,j) = bx(k,j) - uu*app;
                    cx(k,j) = cx(k,j) + uu*(app-apm);
                    dx(k,j) = dx(k,j) + uu*apm;
                }
                for (int k = 1; k <= kdim-1; k++) {
                    double fact = timestp(j,k,i);
                    bx(k,j) = bx(k,j)*fact;
                    cx(k,j) = cx(k,j)*fact + 1.0*(1.+phi);
                    dx(k,j) = dx(k,j)*fact;
                    fx(k,j) = vist3d(j,k,i)*(1.+phi);
                }
            } // end JDIM boundary
            if (iover == 1) {
                for (int k = 1; k <= kdim-1; k++)
                for (int j = 1; j <= jdim-1; j++) {
                    fx(k,j) = fx(k,j)*blank(j,k,i);
                    bx(k,j) = bx(k,j)*blank(j,k,i);
                    dx(k,j) = dx(k,j)*blank(j,k,i);
                    cx(k,j) = cx(k,j)*blank(j,k,i) + (1.-blank(j,k,i));
                }
            }
            {
                int kl2=1, ku2=kdim-1, jl2=1, ju2=jdim-1;
                int kdim1=kdim-1, jdim1=jdim-1;
                triv_ns::triv(kdim1, jdim1, kl2, ku2, jl2, ju2, workx, bx, cx, dx, fx);
            }
            for (int j = 1; j <= jdim-1; j++)
            for (int k = 1; k <= kdim-1; k++) {
                vist3d(j,k,i) = fx(k,j);
            }
        } // end implicit xi loop over i

        // Implicit F_zeta_zeta viscous terms
        if (i2d != 1 && iaxi2planeturb != 1) {
            for (int j = 1; j <= jdim-1; j++) {
                // Interior points
                for (int i = 2; i <= idim-2; i++) {
                    int il = i-1;
                    int iu = i+1;
                    for (int k = 1; k <= kdim-1; k++) {
                        double voliu = vol(j,k,iu);
                        double xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        double yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        double zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        double volil = vol(j,k,il);
                        double xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        double ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        double zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        double xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        double ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        double za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        double ttpo = xp*xa+yp*ya+zp*za;
                        double ttmo = xm*xa+ym*ya+zm*za;
                        double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                        double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                        double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                        double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                        double cnud = -cb2*turre(j,k,i)/(sigma*re);
                        double cap = ttp*cnud;
                        double cam = ttm*cnud;
                        double anutp = .5*(turre(j,k,i+1)+turre(j,k,i));
                        double anutm = .5*(turre(j,k,i-1)+turre(j,k,i));
                        double fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                        double fnum_v = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                        double diff_facp, diff_facm;
                        if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                            double chi_p = anutp/fnup;
                            double chi_m = anutm/fnum_v;
                            diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                            diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                        } else {
                            diff_facp = 1.0;
                            diff_facm = 1.0;
                        }
                        double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                        double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                        double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                        bz(k,i) = -ccmaxcr(tmp1, 0.);
                        cz(k,i) = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                        dz(k,i) = -ccmaxcr(tmp2, 0.);
                    }
                    for (int k = 1; k <= kdim-1; k++) {
                        double xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        double yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        double zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        double tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                        double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                        double one = 1.; double sgnu = ccsignrc(one, uu);
                        double app = 0.5*(1.+sgnu);
                        double apm = 0.5*(1.-sgnu);
                        bz(k,i) = bz(k,i) - uu*app;
                        cz(k,i) = cz(k,i) + uu*(app-apm);
                        dz(k,i) = dz(k,i) + uu*apm;
                    }
                    for (int k = 1; k <= kdim-1; k++) {
                        double fact = timestp(j,k,i);
                        bz(k,i) = bz(k,i)*fact;
                        cz(k,i) = cz(k,i)*fact + 1.0*(1.+phi);
                        dz(k,i) = dz(k,i)*fact;
                        fz(k,i) = vist3d(j,k,i)*(1.+phi);
                    }
                } // end interior i loop

                // I0 boundary points
                {
                    int i = 1;
                    int iu = std::min(2, idim-1);
                    for (int k = 1; k <= kdim-1; k++) {
                        double voliu = vol(j,k,iu);
                        double xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        double yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        double zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        double volil = voli0(j,k,1);
                        double xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        double ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        double zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        double xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        double ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        double za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        double ttpo = xp*xa+yp*ya+zp*za;
                        double ttmo = xm*xa+ym*ya+zm*za;
                        double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                        double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                        double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                        double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                        double cnud = -cb2*turre(j,k,i)/(sigma*re);
                        double cap = ttp*cnud;
                        double cam = ttm*cnud;
                        double anutp = .5*(turre(j,k,i+1)+turre(j,k,i));
                        double anutm = .5*(turre(j,k,i-1)+turre(j,k,i));
                        double fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                        double fnum_v = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                        double diff_facp, diff_facm;
                        if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                            double chi_p = anutp/fnup;
                            double chi_m = anutm/fnum_v;
                            diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                            diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                        } else {
                            diff_facp = 1.0;
                            diff_facm = 1.0;
                        }
                        double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                        double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                        double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                        bz(k,i) = -ccmaxcr(tmp1, 0.);
                        cz(k,i) = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                        dz(k,i) = -ccmaxcr(tmp2, 0.);
                    }
                    for (int k = 1; k <= kdim-1; k++) {
                        double xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        double yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        double zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        double tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                        double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                        double one = 1.; double sgnu = ccsignrc(one, uu);
                        double app = 0.5*(1.+sgnu);
                        double apm = 0.5*(1.-sgnu);
                        bz(k,i) = bz(k,i) - uu*app;
                        cz(k,i) = cz(k,i) + uu*(app-apm);
                        dz(k,i) = dz(k,i) + uu*apm;
                    }
                    for (int k = 1; k <= kdim-1; k++) {
                        double fact = timestp(j,k,i);
                        bz(k,i) = bz(k,i)*fact;
                        cz(k,i) = cz(k,i)*fact + 1.0*(1.+phi);
                        dz(k,i) = dz(k,i)*fact;
                        fz(k,i) = vist3d(j,k,i)*(1.+phi);
                    }
                } // end I0 boundary

                // IDIM boundary points
                {
                    int i = idim-1;
                    int il = idim-2;
                    for (int k = 1; k <= kdim-1; k++) {
                        double voliu = voli0(j,k,3);
                        double xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        double yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        double zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        double volil = vol(j,k,il);
                        double xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        double ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        double zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        double xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        double ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        double za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        double ttpo = xp*xa+yp*ya+zp*za;
                        double ttmo = xm*xa+ym*ya+zm*za;
                        double ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                        double ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                        double ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                        double ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                        double cnud = -cb2*turre(j,k,i)/(sigma*re);
                        double cap = ttp*cnud;
                        double cam = ttm*cnud;
                        double anutp = .5*(turre(j,k,i+1)+turre(j,k,i));
                        double anutm = .5*(turre(j,k,i-1)+turre(j,k,i));
                        double fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                        double fnum_v = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                        double diff_facp, diff_facm;
                        if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
                            double chi_p = anutp/fnup;
                            double chi_m = anutm/fnum_v;
                            diff_facp = (16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p);
                            diff_facm = (16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m);
                        } else {
                            diff_facp = 1.0;
                            diff_facm = 1.0;
                        }
                        double cdp = (fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re);
                        double cdm = (fnum_v+(diff_facm+cb2)*anutm)*ttm/(sigma*re);
                        double tmp1 = cdm+cam; double tmp2 = cdp+cap;
                        bz(k,i) = -ccmaxcr(tmp1, 0.);
                        cz(k,i) = ccmaxcr(tmp2, 0.) + ccmaxcr(tmp1, 0.);
                        dz(k,i) = -ccmaxcr(tmp2, 0.);
                    }
                    for (int k = 1; k <= kdim-1; k++) {
                        double xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        double yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        double zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        double tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                        double uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                        double one = 1.; double sgnu = ccsignrc(one, uu);
                        double app = 0.5*(1.+sgnu);
                        double apm = 0.5*(1.-sgnu);
                        bz(k,i) = bz(k,i) - uu*app;
                        cz(k,i) = cz(k,i) + uu*(app-apm);
                        dz(k,i) = dz(k,i) + uu*apm;
                    }
                    for (int k = 1; k <= kdim-1; k++) {
                        double fact = timestp(j,k,i);
                        bz(k,i) = bz(k,i)*fact;
                        cz(k,i) = cz(k,i)*fact + 1.0*(1.+phi);
                        dz(k,i) = dz(k,i)*fact;
                        fz(k,i) = vist3d(j,k,i)*(1.+phi);
                    }
                } // end IDIM boundary
                if (iover == 1) {
                    for (int i = 1; i <= idim-1; i++)
                    for (int k = 1; k <= kdim-1; k++) {
                        fz(k,i) = fz(k,i)*blank(j,k,i);
                        bz(k,i) = bz(k,i)*blank(j,k,i);
                        dz(k,i) = dz(k,i)*blank(j,k,i);
                        cz(k,i) = cz(k,i)*blank(j,k,i) + (1.-blank(j,k,i));
                    }
                }
                {
                    int kl2=1, ku2=kdim-1, il2=1, iu2=idim-1;
                    int kdim1=kdim-1, idim1=idim-1;
                    triv_ns::triv(kdim1, idim1, kl2, ku2, il2, iu2, workz, bz, cz, dz, fz);
                }
                for (int i = 1; i <= idim-1; i++)
                for (int k = 1; k <= kdim-1; k++) {
                    vist3d(j,k,i) = fz(k,i);
                }
            } // end implicit zeta loop over j
        } // end if i2d != 1 && iaxi2planeturb != 1 (implicit zeta)

        // Update turre (Fortran spalart.F ~2466: when NOT the SA-neg variant,
        // floor the updated value to 1.e-12 — NOT max(turre,0). Flooring to 0
        // lets turre reach ~1e-20 and go NaN the next cycle on coarse mesh-seq
        // levels; the 1.e-12 floor keeps the SA field alive so it can regrow.)
        // Fortran spalart.F ~2462-2489: reset sumn/negn, update turre, and
        // accumulate sumn = sum(vist3d(increment)^2) over the field. In the
        // sub-1e-12 branch vist3d (the SA increment) is zeroed BEFORE the sum,
        // so those cells contribute nothing to the turbulence residual.
        sumn = 0.;
        negn = 0;
        for (int i = 1; i <= idim-1; i++)
        for (int k = 1; k <= kdim-1; k++)
        for (int j = 1; j <= jdim-1; j++) {
            double newval = turre(j,k,i) + vist3d(j,k,i);
            if (i_saneg != 1 && (float)newval < 1.0e-12f) {
                negn = negn + 1;
                turre(j,k,i) = 1.0e-12;
                vist3d(j,k,i) = 0.;
            } else if (i_saneg == 1 && (float)newval < 0.0f) {
                negn = negn + 1;
                turre(j,k,i) = newval;
            } else {
                turre(j,k,i) = newval;
            }
            sumn = sumn + vist3d(j,k,i)*vist3d(j,k,i);
        }
        sumn = std::sqrt(sumn)/(double)((kdim-1)*(jdim-1)*(idim-1));

        // NOTE: Fortran spalart.F applies NO tur1cut/tur2cut clamp on turre
        // (those are the 2-eqn models' default eps/omega/k seed values, and
        // for the SA model default to 1.e-20). A spurious translated
        // `turre=min(turre,tur1cut=1e-20)` collapsed the whole SA field to
        // ~1e-20 on coarse mesh-seq levels → vist3d≈0 → divergence → NaN.
        // Removed to match Fortran.

        // Laminar region: force turre to zero
        if (ilamlo != 0 && jlamlo != 0 && klamlo != 0) {
            if (i_lam_forcezero == 1) {
                for (int i = ilamlo; i <= ilamhi-1; i++)
                for (int k = klamlo; k <= klamhi-1; k++)
                for (int j = jlamlo; j <= jlamhi-1; j++) {
                    turre(j,k,i) = 0.;
                }
            }
        }

    } // end not_ loop (subiterations)

    // Compute vist3d from turre
    for (int i = 1; i <= idim-1; i++)
    for (int k = 1; k <= kdim-1; k++)
    for (int j = 1; j <= jdim-1; j++) {
        double chi = turre(j,k,i)/fnu(j,k,i);
        double fv1;
        if (i_saneg == 1 && (float)turre(j,k,i) < 0.0f) {
            fv1 = 0.;
        } else {
            fv1 = chi*chi*chi/(chi*chi*chi+cv1*cv1*cv1);
        }
        // Fortran spalart.F 2513: vist3d = fv1*turre*q(j,k,i,1)  (mu_t = rho*fv1*nu~).
        // The earlier /re form made the eddy viscosity ~re(=reue/xmach) times too
        // small — freestream 3e-10 instead of 0.009 — i.e. effectively laminar.
        vist3d(j,k,i) = fv1*turre(j,k,i)*q(j,k,i,1);
        double tmp_zero = 0.;
        vist3d(j,k,i) = ccmax(vist3d(j,k,i), tmp_zero);
    }

    // Limit vist3d
    if ((float)cmn_turbconv.edvislim > 0.0f) {
        for (int i = 1; i <= idim-1; i++)
        for (int k = 1; k <= kdim-1; k++)
        for (int j = 1; j <= jdim-1; j++) {
            double tmp_lim = (double)cmn_turbconv.edvislim;
            vist3d(j,k,i) = ccmincr(vist3d(j,k,i), tmp_lim);
        }
    }

    // Laminar region: force vist3d to zero
    if (ilamlo != 0 && jlamlo != 0 && klamlo != 0) {
        if (i_lam_forcezero == 1) {
            for (int i = ilamlo; i <= ilamhi-1; i++)
            for (int k = klamlo; k <= klamhi-1; k++)
            for (int j = jlamlo; j <= jlamhi-1; j++) {
                vist3d(j,k,i) = 0.;
            }
        }
    }

    // Save turre to tursav
    for (int i = 1; i <= idim-1; i++)
    for (int k = 1; k <= kdim-1; k++)
    for (int j = 1; j <= jdim-1; j++) {
        tursav(j,k,i,1) = turre(j,k,i);
    }
    // Fortran spalart.F 2494-2497: sumn1 = sumn (RMS of the SA increment,
    // accumulated in the update loop), sumn2 = 1, negn1 = negn (clipped-cell
    // count), negn2 = 0. The turbulence residual (turres) is |sumn1|.
    sumn1 = sumn;
    sumn2 = 1.;
    negn1 = negn;
    negn2 = 0;

    // Compute exact MMS solution if needed
    if (iexact_trunc != 0 || iexact_disc != 0) {
        // tursav(j,k,i,1) is the SA variable
        // exact_turb_q expects 3D array for q
        // In Fortran: call exact_turb_q(jdim,kdim,idim,x,y,z,tursav,smin,vist3d,iexact_trunc,iexact_disc)
        // tursav is 4D but passed as 3D - use a 3D view of tursav(:,:,:,1)
        // We create a temporary 3D array for this
        FortranArray3D<double> tursav_3d(jdim, kdim, idim);
        for (int i = 1; i <= idim-1; i++)
        for (int k = 1; k <= kdim-1; k++)
        for (int j = 1; j <= jdim-1; j++) {
            tursav_3d(j,k,i) = tursav(j,k,i,1);
        }
        FortranArray3DRef<double> tursav_3d_ref = tursav_3d.ref();
        mms_ns::exact_turb_q(jdim, kdim, idim, x, y, z, tursav_3d_ref, smin, vist3d, iexact_trunc, iexact_disc);
        if (iexact_ring != 0) {
            mms_ns::exact_turb_q_ring(jdim, kdim, idim, x, y, z, tursav_3d_ref, smin, vist3d, iexact_trunc, iexact_disc);
        }
    }

    // zero_resid_ring: vist3d passed as 4D
    if (iexact_ring != 0) {
        // vist3d is 3D(jdim,kdim,idim), zero_resid_ring expects 4D
        // In Fortran: call zero_resid_ring(jdim,kdim,idim,vist3d,jdim,kdim,idim,1,2,iexact_trunc,iexact_disc)
        // We create a 4D view of vist3d data
        FortranArray4D<double> vist3d_4d(jdim, kdim, idim, 1);
        for (int i = 1; i <= idim-1; i++)
        for (int k = 1; k <= kdim-1; k++)
        for (int j = 1; j <= jdim-1; j++) {
            vist3d_4d(j,k,i,1) = vist3d(j,k,i);
        }
        FortranArray4DRef<double> vist3d_4d_ref = vist3d_4d.ref();
        int numeq = 1, numrows = 2;
        mms_ns::zero_resid_ring(jdim, kdim, idim, vist3d_4d_ref, jdim, kdim, idim, numeq, numrows, iexact_trunc, iexact_disc);
        for (int i = 1; i <= idim-1; i++)
        for (int k = 1; k <= kdim-1; k++)
        for (int j = 1; j <= jdim-1; j++) {
            vist3d(j,k,i) = vist3d_4d(j,k,i,1);
        }
    }

    // Compute vor (vorticity magnitude) for output
    // (already computed by caller, but update vist3d for output)
    // Print warning if negative values
    if (negn1 > 0) {
        if (cmn_turbconv.iwarneddy == 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   WARNING: %8d negative SA values in block %4d, sum=%15.5e", negn1, nbl, sumn1);
        }
    }

    // (Fortran sets sumn2=1, negn2=0 above; no negative-vist3d accumulation.)

} // end spalart()

} // namespace spalart_ns
