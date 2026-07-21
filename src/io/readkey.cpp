// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "readkey.h"
#include "termn8.h"
#include "parser.h"
#include "trnsfr_vals.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <string>

namespace readkey_ns {

void readkey(int& ititr, int& myid, int& ibufdim, int& nbuf,
             FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou,
             int& iunit11, int& ierrflg)
{
    // local variables
    char inpstr[210];
    float realval[21]; // 1-based: use indices 1..20
    bool isthere;
    float dclda, dcldd, dcmda, dcmdd;
    float turbintensity_inf_percent, eddy_visc_inf;
    int idummy, idummy1, idummy2, idummy3;
    float errornum, errorden, error_val, zkinf;
    int n;
    int npos, lc1, lc2, lcl, lclw;

    // COMMON block aliases
    float& gamma         = cmn_fluid.gamma;
    float& pr            = cmn_fluid2.pr;
    float& prt           = cmn_fluid2.prt;
    float& cbar          = cmn_fluid2.cbar;
    float& atol          = cmn_singular.atol;
    int&   iexp          = cmn_zero.iexp;
    float& xmach_img     = cmn_complx.xmach_img;
    float& alpha_img     = cmn_complx.alpha_img;
    float& beta_img      = cmn_complx.beta_img;
    float& reue_img      = cmn_complx.reue_img;
    float& tinf_img      = cmn_complx.tinf_img;
    float& geom_img      = cmn_complx.geom_img;
    float& surf_img      = cmn_complx.surf_img;
    float& xrotrate_img  = cmn_complx.xrotrate_img;
    float& yrotrate_img  = cmn_complx.yrotrate_img;
    float& zrotrate_img  = cmn_complx.zrotrate_img;
    int&   icgns         = cmn_cgns.icgns;
    float& cprec         = cmn_precond.cprec;
    float& uref          = cmn_precond.uref;
    float& avn           = cmn_precond.avn;
    int&   ialphit       = cmn_alphait.ialphit;
    float& cltarg        = cmn_alphait.cltarg;
    float& rlxalph       = cmn_alphait.rlxalph;
    float& dalim         = cmn_alphait.dalim;
    int&   icycupdt      = cmn_alphait.icycupdt;
    int&   nsubturb      = cmn_turbconv.nsubturb;
    float* cflturb       = cmn_turbconv.cflturb; // 0-based C array, use [0]..[6]
    float& edvislim      = cmn_turbconv.edvislim;
    int&   iturbprod     = cmn_turbconv.iturbprod;
    int&   nfreeze       = cmn_turbconv.nfreeze;
    int&   itime2read    = cmn_turbconv.itime2read;
    int&   itaturb       = cmn_turbconv.itaturb;
    float& tur1cut       = cmn_turbconv.tur1cut;
    float& tur2cut       = cmn_turbconv.tur2cut;
    int&   iturbord      = cmn_turbconv.iturbord;
    float& tur1cutlev    = cmn_turbconv.tur1cutlev;
    float& tur2cutlev    = cmn_turbconv.tur2cutlev;
    float& epsa_r        = cmn_entfix.epsa_r;
    float& epsa_l        = cmn_entfix.epsa_l;
    int&   nkey          = cmn_key.nkey;
    int&   ivolint       = cmn_is_blockbc.ivolint;
    int&   idef_ss       = cmn_elastic_ss.idef_ss;
    int&   ibin          = cmn_bin.ibin;
    int&   iblnk         = cmn_bin.iblnk;
    int&   iblnkfr       = cmn_bin.iblnkfr;
    int&   ip3dgrad      = cmn_bin.ip3dgrad;
    int&   memadd        = cmn_memory.memadd;
    int&   memaddi       = cmn_memory.memaddi;
    int&   negvol        = cmn_deformz.negvol;
    int&   meshdef       = cmn_deformz.meshdef;
    int&   ndgrd         = cmn_deformz.ndgrd;
    int&   ndwrt         = cmn_deformz.ndwrt;
    int&   irghost       = cmn_ghost.irghost;
    int&   iwghost       = cmn_ghost.iwghost;
    float& xcentrot      = cmn_noninertial.xcentrot;
    float& ycentrot      = cmn_noninertial.ycentrot;
    float& zcentrot      = cmn_noninertial.zcentrot;
    float& xrotrate      = cmn_noninertial.xrotrate;
    float& yrotrate      = cmn_noninertial.yrotrate;
    float& zrotrate      = cmn_noninertial.zrotrate;
    int&   noninflag     = cmn_noninertial.noninflag;
    float* tur10         = cmn_ivals.tur10; // 0-based C array, use [0]..[6]
    float& roll_angle    = cmn_gridtrans.roll_angle;
    int&   ikoprod       = cmn_konew.ikoprod;
    int&   isstdenom     = cmn_konew.isstdenom;
    float& pklimterm     = cmn_konew.pklimterm;
    int&   iaxi2plane    = cmn_axisym.iaxi2plane;
    int&   iturbord_ref  = cmn_turbconv.iturbord; // already aliased above
    int&   ifullns       = cmn_fullns.ifullns;
    int&   ibeta8kzeta   = cmn_konew.ibeta8kzeta;
    int&   isarc2d       = cmn_curvat.isarc2d;
    int&   isarc3d       = cmn_curvat.isarc3d;
    float& sarccr3       = cmn_curvat.sarccr3;
    int&   ieasmcc2d     = cmn_curvat.ieasmcc2d;
    int&   ipertavg      = cmn_avgdata.ipertavg;
    int&   iteravg       = cmn_avgdata.iteravg;
    int&   icoarsemovie  = cmn_moov.icoarsemovie;
    int&   i2dmovie      = cmn_moov.i2dmovie;
    int&   iclcd         = cmn_avgdata.iclcd;
    int&   iskip_blocks  = cmn_ginfo2.iskip_blocks;
    float& cfltauMax     = cmn_unst.cfltaumax;
    float& cfltau0       = cmn_unst.cfltau0;
    int&   irbtrim       = cmn_rigidbody.irbtrim;
    int&   irigb         = cmn_rigidbody.irigb;
    float& greflrb       = cmn_rbstmt2.greflrb;
    float& tmass         = cmn_rbstmt2.tmass;
    float& yinert        = cmn_rbstmt2.yinert;
    float& gaccel        = cmn_rbstmt2.gaccel;
    float& relax         = cmn_trim.relax;
    int&   itrminc       = cmn_trim.itrminc;
    int&   i_bsl         = cmn_konew.i_bsl;
    int&   keepambient   = cmn_konew.keepambient;
    float& re_thetat0    = cmn_konew.re_thetat0;
    float& cs_smagorinsky= cmn_lesinfo.cs_smagorinsky;
    float& xdir_only_source = cmn_sourceterm.xdir_only_source;
    float& randomize     = cmn_random_input.randomize;
    int&   iexact_trunc  = cmn_mms.iexact_trunc;
    int&   iexact_disc   = cmn_mms.iexact_disc;
    int&   iexact_ring   = cmn_mms.iexact_ring;
    int&   i_wilcox06    = cmn_konew.i_wilcox06;
    int&   i_wilcox06_chiw = cmn_konew.i_wilcox06_chiw;
    int&   i_turbprod_kterm = cmn_konew.i_turbprod_kterm;
    int&   i_catris_kw   = cmn_konew.i_catris_kw;
    int&   ismincforce   = cmn_sminn.ismincforce;
    float& prod2d3dtrace = cmn_konew.prod2d3dtrace;
    int&   i_compress_correct = cmn_konew.i_compress_correct;
    int&   les_model     = cmn_lesinfo.les_model;
    int&   les_wallscale = cmn_lesinfo.les_wallscale;
    float& cs_wale       = cmn_lesinfo.cs_wale;
    float& cs_vreman     = cmn_lesinfo.cs_vreman;
    int&   isstrc        = cmn_curvat.isstrc;
    float& sstrc_crc     = cmn_curvat.sstrc_crc;
    int&   isstsf        = cmn_konew.isstsf;
    float& scal_ic       = cmn_initfac.scal_ic;
    int&   ifunct        = cmn_plot3dtyp.ifunct;
    int&   lowmem_ux     = cmn_memry.lowmem_ux;
    int&   isar          = cmn_curvat.isar;
    float& crot          = cmn_curvat.crot;
    int&   i_nonlin      = cmn_constit.i_nonlin;
    float& c_nonlin      = cmn_constit.c_nonlin;
    float& snonlin_lim   = cmn_constit.snonlin_lim;
    int&   isubit_r      = cmn_avgdata.isubit_r;
    int&   i_wilcox98    = cmn_konew.i_wilcox98;
    int&   i_wilcox98_chiw = cmn_konew.i_wilcox98_chiw;
    float& cmulim        = cmn_easmlim.cmulim;
    int&   iaxi2planeturb= cmn_axisym.iaxi2planeturb;
    int&   istrongturbdis= cmn_axisym.istrongturbdis;
    int&   ieasm_type    = cmn_easmv.ieasm_type;
    int&   ipatch1st     = cmn_is_patch.ipatch1st;
    int&   isst2003      = cmn_konew.isst2003;
    int&   issglrrw2012  = cmn_reystressmodel.issglrrw2012;
    int&   ifort50write  = cmn_writestuff.ifort50write;
    int&   j_ifort50write= cmn_writestuff.j_ifort50write;
    int&   i_ifort50write= cmn_writestuff.i_ifort50write;
    int&   i_sas_rsm     = cmn_reystressmodel.i_sas_rsm;
    int&   iforcev0      = cmn_axisym.iforcev0;
    int&   i_saneg       = cmn_sa_options.i_saneg;
    int&   i_sanoft2     = cmn_sa_options.i_sanoft2;
    int&   i_lam_forcezero = cmn_lam.i_lam_forcezero;
    int&   i_specialtop_kmax1001 = cmn_specialtop_kmax1001.i_specialtop_kmax1001;
    float& a_specialtop_kmax1001 = cmn_specialtop_kmax1001.a_specialtop_kmax1001;
    float& xc_specialtop_kmax1001 = cmn_specialtop_kmax1001.xc_specialtop_kmax1001;
    float& sig_specialtop_kmax1001 = cmn_specialtop_kmax1001.sig_specialtop_kmax1001;
    float& vtp_specialtop_kmax1001 = cmn_specialtop_kmax1001.vtp_specialtop_kmax1001;
    float& wc_specialtop_kmax1001 = cmn_specialtop_kmax1001.wc_specialtop_kmax1001;
    float& fac_specialtop_kmax1001 = cmn_specialtop_kmax1001.fac_specialtop_kmax1001;
    float& cc_specialtop_kmax1001 = cmn_specialtop_kmax1001.cc_specialtop_kmax1001;
    float& xerf_specialtop_kmax1001 = cmn_specialtop_kmax1001.xerf_specialtop_kmax1001;
    float& sigerf_specialtop_kmax1001 = cmn_specialtop_kmax1001.sigerf_specialtop_kmax1001;
    float& sa_cw2        = cmn_sa_options.sa_cw2;
    float& sa_cw3        = cmn_sa_options.sa_cw3;
    float& sa_cv1        = cmn_sa_options.sa_cv1;
    float& sa_ct3        = cmn_sa_options.sa_ct3;
    float& sa_ct4        = cmn_sa_options.sa_ct4;
    float& sa_cb1        = cmn_sa_options.sa_cb1;
    float& sa_cb2        = cmn_sa_options.sa_cb2;
    float& sa_sigma      = cmn_sa_options.sa_sigma;
    float& sa_karman     = cmn_sa_options.sa_karman;
    int&   iupdatemean   = cmn_iupdate.iupdatemean;
    float& uub_2034      = cmn_wallvel_2034.uub_2034;
    float& vvb_2034      = cmn_wallvel_2034.vvb_2034;
    float& wwb_2034      = cmn_wallvel_2034.wwb_2034;
    float& vvb_xaxisrot2034 = cmn_wallvel_2034.vvb_xaxisrot2034;
    int&   i_tauijs      = cmn_constit.i_tauijs;
    int&   i_qcr2000     = cmn_constit.i_qcr2000;
    int&   i_qcr2013     = cmn_constit.i_qcr2013;
    int&   i_qcr2013v    = cmn_constit.i_qcr2013v;
    int&   i_yapterm     = cmn_reystressmodel.i_yapterm;
    float& xmach         = cmn_info.xmach;
    int&   ivmx          = cmn_maxiv.ivmx;
    int&   iunit5        = cmn_unit5.iunit5;
    char*  output_dir    = cmn_filenam2.output_dir;
    // tp is column-major 4x4 stored as float[16] in cmn_trim.tp
    // tp(i,j) in Fortran = cmn_trim.tp[(j-1)*4 + (i-1)] in C
    float* tp_arr        = cmn_trim.tp;



    // ititr = 0
    ititr = 0;
    nkey  = 0;

    // ratio of specific heats
    gamma  = 1.4f;
    nkey   = nkey + 1;
    // Prandtl number
    pr     = 0.72f;
    nkey   = nkey + 1;
    // Turbulent Prandtl number
    prt    = 0.9f;
    nkey   = nkey + 1;
    // reference temperature for Sutherland's Law
    cbar   = 198.6f;
    nkey   = nkey + 1;
    // tolerance for collapsed metrics
    atol   = std::max(1.e-07f, (float)std::pow(10.0, (double)(-iexp+1)));
    nkey   = nkey + 1;
    // complex perturbation to the Mach number
    xmach_img = 0.f;
    nkey      = nkey + 1;
    // complex perturbation to angle of attack, alpha
    alpha_img = 0.f;
    nkey      = nkey + 1;
    // complex perturbation to yaw angle, beta
    beta_img = 0.f;
    nkey     = nkey + 1;
    // complex perturbation to the unit Reynolds number, reue
    reue_img = 0.f;
    nkey     = nkey + 1;
    // complex perturbation to the freestream temperature, tinf
    tinf_img = 0.f;
    nkey     = nkey + 1;
    // size of complex perturbation to grid
    geom_img  = 0.f;
    nkey      = nkey + 1;
    // CGNS file flag
    icgns     = 0;
    nkey      = nkey + 1;
    // relative amount of preconditioning
    cprec     = 0.f;
    nkey      = nkey + 1;
    // limiting velocity for preconditioning
    uref      = xmach;
    nkey      = nkey + 1;
    // multiplicative factor for uref**2
    avn      = 1.0f;
    nkey     = nkey + 1;
    // specified Cl
    cltarg   = 99999.0f;
    ialphit  = 0;
    nkey     = nkey + 1;
    // angle-of-attack relaxation factor for specified Cl
    rlxalph  = 1.f;
    nkey     = nkey + 1;
    // no. of subiterations in turbulence model per outer iteration
    nsubturb = 1;
    nkey     = nkey + 1;
    // factor governing cfl number for turbulence model
    cflturb[0]  = 0.f;
    cflturb[1]  = 0.f;
    cflturb[2]  = 0.f;
    cflturb[3]  = 0.f;
    cflturb[4]  = 0.f;
    cflturb[5]  = 0.f;
    cflturb[6]  = 0.f;
    nkey     = nkey + 7;
    // limit on ratio of minimum to maximum |eigenvalue|
    epsa_r   = 0.f;
    epsa_l   = 2.f*epsa_r;
    nkey     = nkey + 1;
    // number of cycles over which to freeze turbulence model
    nfreeze  = 0;
    nkey     = nkey + 1;
    // flag for using/not using the exact volume terms on 1-1 block boundaries
    ivolint  = 1;
    nkey     = nkey + 1;
    // flag for mesh deformation in steady-state mode
    idef_ss  = 0;
    nkey     = nkey + 1;
    // flag for writing unformated/formatted plot3d files
    ibin     = 1;
    nkey     = nkey + 1;
    // flag for enabling/disabling the writing of the iblank array
    iblnk    = 1;
    nkey     = nkey + 1;
    // flag for enabling/disabling the blanking of fringe points
    iblnkfr  = 1;
    nkey     = nkey + 1;
    // flag to switch from solution output to derivative output
    ip3dgrad = 0;
    nkey     = nkey + 1;
    // size of complex perturbation to surface grid
    surf_img  = 0.f;
    nkey      = nkey + 1;
    // additional user-specified memory
    memadd    = 0;
    nkey      = nkey + 1;
    memaddi   = 0;
    nkey      = nkey + 1;
    // flag to turn off stops when negative volumes/bad metrics
    negvol    = 0;
    nkey      = nkey + 1;
    // flag to turn on/off time marching flow solver
    meshdef   = 0;
    nkey      = nkey + 1;
    // eddy viscosity limiter for two eqn. turbulence models
    edvislim  = 1.e10f;
    nkey      = nkey + 1;
    // flag to set whether approximate or full production term is used
    iturbprod = 0;
    nkey      = nkey + 1;
    // flag to read ghost-cell data from restart file
    irghost   = 1;
    nkey      = nkey + 1;
    // flag to write ghost-cell data to restart file
    iwghost   = 1;
    nkey      = nkey + 1;
    // limit of alpha change (deg.) per update
    dalim     = 0.2f;
    nkey      = nkey + 1;
    // no. of cycles between alpha updates
    icycupdt  = 1;
    nkey      = nkey + 1;
    // non-inertial reference frame flag
    noninflag = 0;
    nkey      = nkey + 1;
    // rotation center x-coordinate for non-inertial reference frame
    xcentrot = 0.f;
    nkey     = nkey + 1;
    // rotation center y-coordinate for non-inertial reference frame
    ycentrot = 0.f;
    nkey     = nkey + 1;
    // rotation center z-coordinate for non-inertial reference frame
    zcentrot = 0.f;
    nkey     = nkey + 1;
    // rotation rate in x-direction for non-inertial reference frame
    xrotrate = 0.f;
    nkey     = nkey + 1;
    // rotation rate in y-direction for non-inertial reference frame
    yrotrate = 0.f;
    nkey     = nkey + 1;
    // rotation rate in z-direction for non-inertial reference frame
    zrotrate = 0.f;
    nkey     = nkey + 1;
    // complex perturbation to the rotation rate in x-direction
    xrotrate_img = 0.f;
    nkey      = nkey + 1;
    // complex perturbation to the rotation rate in y-direction
    yrotrate_img = 0.f;
    nkey      = nkey + 1;
    // complex perturbation to the rotation rate in z-direction
    zrotrate_img = 0.f;
    nkey      = nkey + 1;
    // flag to read (1) or skip reading (0) 2nd order-time turbulence terms
    itime2read = 1;
    nkey     = nkey + 1;
    // flag to control time-accuracy of turbulence model
    itaturb = 1;
    nkey     = nkey + 1;

    // ititr = 0
    ititr = 0;
    nkey  = 0;
    // ratio of specific heats
    gamma  = 1.4f;
    nkey   = nkey + 1;
    // Prandtl number
    pr     = 0.72f;
    nkey   = nkey + 1;
    // Turbulent Prandtl number
    prt    = 0.9f;
    nkey   = nkey + 1;
    // reference temperature for Sutherland's Law
    cbar   = 198.6f;
    nkey   = nkey + 1;
    // tolerance for collapsed metrics
    atol   = std::max(1.e-07f, (float)std::pow(10.0, (double)(-iexp+1)));
    nkey   = nkey + 1;
    // complex perturbation to the Mach number
    xmach_img = 0.f;
    nkey      = nkey + 1;
    // complex perturbation to angle of attack, alpha
    alpha_img = 0.f;
    nkey      = nkey + 1;
    // complex perturbation to yaw angle, beta
    beta_img = 0.f;
    nkey     = nkey + 1;
    // complex perturbation to the unit Reynolds number, reue
    reue_img = 0.f;
    nkey     = nkey + 1;
    // complex perturbation to the freestream temperature, tinf
    tinf_img = 0.f;
    nkey     = nkey + 1;
    // size of complex perturbation to grid
    geom_img  = 0.f;
    nkey      = nkey + 1;
    // CGNS file flag
    icgns     = 0;
    nkey      = nkey + 1;
    // relative amount of preconditioning
    cprec     = 0.f;
    nkey      = nkey + 1;
    // limiting velocity for preconditioning
    uref      = xmach;
    nkey      = nkey + 1;
    // multiplicative factor for uref**2
    avn      = 1.0f;
    nkey     = nkey + 1;
    // specified Cl
    cltarg   = 99999.0f;
    ialphit  = 0;
    nkey     = nkey + 1;
    // angle-of-attack relaxation factor for specified Cl
    rlxalph  = 1.f;
    nkey     = nkey + 1;
    // no. of subiterations in turbulence model per outer iteration
    nsubturb = 1;
    nkey     = nkey + 1;
    // factor governing cfl number for turbulence model
    cflturb[0]  = 0.f;
    cflturb[1]  = 0.f;
    cflturb[2]  = 0.f;
    cflturb[3]  = 0.f;
    cflturb[4]  = 0.f;
    cflturb[5]  = 0.f;
    cflturb[6]  = 0.f;
    nkey     = nkey + 7;
    // limit on ratio of minimum to maximum |eigenvalue|
    epsa_r   = 0.f;
    epsa_l   = 2.f*epsa_r;
    nkey     = nkey + 1;
    // number of cycles over which to freeze turbulence model
    nfreeze  = 0;
    nkey     = nkey + 1;
    // flag for using/not using the exact volume terms on 1-1 block boundaries
    ivolint  = 1;
    nkey     = nkey + 1;
    // flag for mesh deformation in steady-state mode
    idef_ss  = 0;
    nkey     = nkey + 1;
    // flag for writing unformated/formatted plot3d files
    ibin     = 1;
    nkey     = nkey + 1;
    // flag for enabling/disabling the writing of the iblank array
    iblnk    = 1;
    nkey     = nkey + 1;
    // flag for enabling/disabling the blanking of fringe points
    iblnkfr  = 1;
    nkey     = nkey + 1;
    // flag to switch from solution output to derivative output
    ip3dgrad = 0;
    nkey     = nkey + 1;
    // size of complex perturbation to surface grid
    surf_img  = 0.f;
    nkey      = nkey + 1;
    // additional user-specified memory
    memadd    = 0;
    nkey      = nkey + 1;
    memaddi   = 0;
    nkey      = nkey + 1;
    // flag to turn off stops when negative volumes/bad metrics
    negvol    = 0;
    nkey      = nkey + 1;
    // flag to turn on/off time marching flow solver
    meshdef   = 0;
    nkey      = nkey + 1;
    // eddy viscosity limiter for two eqn. turbulence models
    edvislim  = 1.e10f;
    nkey      = nkey + 1;
    // flag to set whether approximate or full production term is used
    iturbprod = 0;
    nkey      = nkey + 1;
    // flag to read ghost-cell data from restart file
    irghost   = 1;
    nkey      = nkey + 1;
    // flag to write ghost-cell data to restart file
    iwghost   = 1;
    nkey      = nkey + 1;
    // limit of alpha change (deg.) per update
    dalim     = 0.2f;
    nkey      = nkey + 1;
    // no. of cycles between alpha updates
    icycupdt  = 1;
    nkey      = nkey + 1;
    // non-inertial reference frame flag
    noninflag = 0;
    nkey      = nkey + 1;
    // rotation center x-coordinate for non-inertial reference frame
    xcentrot = 0.f;
    nkey     = nkey + 1;
    // rotation center y-coordinate for non-inertial reference frame
    ycentrot = 0.f;
    nkey     = nkey + 1;
    // rotation center z-coordinate for non-inertial reference frame
    zcentrot = 0.f;
    nkey     = nkey + 1;

    // rotation rate in x-direction for non-inertial reference frame
    xrotrate = 0.f;
    nkey     = nkey + 1;
    // rotation rate in y-direction for non-inertial reference frame
    yrotrate = 0.f;
    nkey     = nkey + 1;
    // rotation rate in z-direction for non-inertial reference frame
    zrotrate = 0.f;
    nkey     = nkey + 1;
    // complex perturbation to the rotation rate in x-direction
    xrotrate_img = 0.f;
    nkey      = nkey + 1;
    // complex perturbation to the rotation rate in y-direction
    yrotrate_img = 0.f;
    nkey      = nkey + 1;
    // complex perturbation to the rotation rate in z-direction
    zrotrate_img = 0.f;
    nkey      = nkey + 1;
    // flag to read (1) or skip reading (0) 2nd order-time turbulence terms
    itime2read = 1;
    nkey     = nkey + 1;
    // flag to control time-accuracy of turbulence model
    itaturb = 1;
    nkey     = nkey + 1;
    // flag to perform DES with turbulence model
    cmn_des.ides = 0;
    nkey     = nkey + 1;
    // constant associated with DES or DDES
    cmn_des.cdes = 0.65f;
    nkey     = nkey + 1;
    // additional constant associated with Modified DDES
    cmn_des.cddes = 0.975f;
    nkey     = nkey + 1;
    // flag to store iteration-averaged q values in a PLOT3D file
    iteravg = 0;
    nkey     = nkey + 1;
    // turbulent quantity freestream levels
    tur10[0] = -1.f;
    tur10[1] = -1.f;
    tur10[2] = -1.f;
    tur10[3] = -1.f;
    tur10[4] = -1.f;
    tur10[5] = -1.f;
    tur10[6] = -1.f;
    nkey     = nkey + 7;
    // value that epsilon or omega is set to for 2-eqn models
    if (ivmx == 15) {
        tur1cut = -1.f;
    } else {
        tur1cut = 1.e-20f;
    }
    nkey     = nkey + 1;
    // roll angle
    roll_angle = 0.f;
    nkey       = nkey + 1;
    // k-omega/sst/k-epsilon parameter controlling production term
    ikoprod = 0;
    nkey       = nkey + 1;
    // sst parameter controlling denominator of mut term
    isstdenom = 0;
    nkey       = nkey + 1;
    // limiter term on Pk in the two-eqn models
    pklimterm = 20.f;
    nkey       = nkey + 1;
    // parameter for use with particular axi cases
    iaxi2plane = 0;
    nkey       = nkey + 1;
    // parameter controls whether turbulence model advection terms are 1st or 2nd order
    iturbord = 1;
    nkey       = nkey + 1;
    // value that k is set to for 2-eqn models
    tur2cut = 1.e-20f;
    nkey     = nkey + 1;
    // lower limit on epsilon or omega for 2-eqn models
    tur1cutlev = 0.f;
    nkey     = nkey + 1;
    // lower limit on k for 2-eqn models
    tur2cutlev = 0.f;
    nkey     = nkey + 1;
    // full Navier-Stokes
    ifullns = 0;
    nkey     = nkey + 1;
    // beta8 parameter in k-enstrophy model
    ibeta8kzeta = 0;
    nkey     = nkey + 1;
    // SARC - Spalart-Allmaras curvature correction in 2D sense only
    isarc2d = 0;
    nkey     = nkey + 1;
    // SARC - Spalart-Allmaras curvature correction in 3D
    isarc3d = 0;
    nkey     = nkey + 1;
    // parameter cr3 in SARC model
    sarccr3 = 1.0f;
    nkey     = nkey + 1;
    // EASMCC - EASM curvature correction in 2D sense only
    ieasmcc2d = 0;
    nkey     = nkey + 1;
    // ipertavg - Average r,u,v,w,p and second moments
    ipertavg = 0;
    nkey     = nkey + 1;
    // icoarsemovie - write movie files of coarse planar and volume data
    icoarsemovie = 0;
    nkey     = nkey + 1;
    // i2dmovie - write 2D movie files of coarse planar and volume data
    i2dmovie = 0;
    nkey     = nkey + 1;
    // iclcd - write cl and cd history for separate blocks using clcd.inp
    iclcd    = 0;
    nkey     = nkey + 1;
    // iskip_blocks - skip factor for 2d coarse movie
    iskip_blocks    = 1;
    nkey     = nkey + 1;
    // cfltauMax - Maximum value for the cfltau during subiterations
    cfltauMax= -1.0f;
    nkey     = nkey + 1;
    // cfltau0 - Exponent for cfltau growth with subiteration
    cfltau0= 1.0f;
    nkey     = nkey + 1;
    // parameter to perform trim (1) or not to perform trim (0)
    irbtrim = 0;
    nkey     = nkey + 1;
    // parameter to perform rigid body sim (1) or not to perform (0)
    irigb    = 0;
    nkey     = nkey + 1;
    // Conversion factor from rigid body dynamics length scale to CFD grid length scale
    greflrb  = 1.f;
    nkey     = nkey + 1;
    // Total aircraft mass for trim and rigid body motion
    tmass = 1.f;
    nkey     = nkey + 1;
    // Total aircraft y-mass moment of inertia for rigid body motion
    yinert = 1.f;
    nkey     = nkey + 1;
    // Acceleration of gravity for rigid/flex trim and motion
    gaccel = 1.f;
    nkey     = nkey + 1;
    // Relaxation parameter for the trim algorithm
    relax = 0.5f;
    nkey     = nkey + 1;
    // Iteration increment to update trim equations
    itrminc = 5;
    nkey     = nkey + 1;
    // dclda to be used for trim
    dclda    = 6.f;
    nkey     = nkey + 1;
    // dcldd to be used for trim
    dcldd    = 1.4f;
    nkey     = nkey + 1;
    // dcmda to be used for trim
    dcmda    = -0.2f;
    nkey     = nkey + 1;
    // dcmdd to be used for trim
    dcmdd    = -0.88f;
    nkey     = nkey + 1;
    // Parameter controling dynamic grid input
    ndgrd    = 0;
    nkey     = nkey + 1;
    // Parameter controling dynamic grid output
    ndwrt    = 0;
    nkey     = nkey + 1;
    // i_bsl used to turn on Menter's BSL model
    i_bsl     = 0;
    nkey     = nkey + 1;
    // keepambient used to maintain ambient turb levels for 2-eq models
    keepambient = 0;
    nkey     = nkey + 1;
    // re_thetat0 used for transition model
    re_thetat0 = -1.0f;
    nkey     = nkey + 1;
    // turbintensity_inf_percent is freestream turb intensity
    turbintensity_inf_percent = -1.0f;
    nkey     = nkey + 1;
    // eddy_visc_inf is freestream mu_t/mu_inf
    eddy_visc_inf = -1.0f;
    nkey     = nkey + 1;
    // cs_smagorinsky is the (non-dynamic) Smagorinsky constant for LES
    cs_smagorinsky = 0.0f;
    nkey     = nkey + 1;
    // xdir_only_source is imposed source term
    xdir_only_source = 0.0f;
    nkey     = nkey + 1;
    // randomize adds random perturbation to restart
    randomize = 0.0f;
    nkey     = nkey + 1;
    // iexact_trunc is used to check truncation error against exact solution
    iexact_trunc = 0;
    nkey     = nkey + 1;
    // iexact_disc is used to check discretization error against exact solution
    iexact_disc = 0;
    nkey     = nkey + 1;
    // iexact_ring is used to overwrite the exact solution on the outer 2 ring layers
    iexact_ring = 0;
    nkey     = nkey + 1;
    // i_wilcox06 is used to change Wilcox88 model to Wilcox06
    i_wilcox06 = 0;
    nkey     = nkey + 1;
    // i_wilcox06_chiw is used to turn on/off Wilcox06 k-omega vortex stretching parameter
    i_wilcox06_chiw = 1;
    nkey     = nkey + 1;
    // i_turbprod_kterm is used to determine whether 2/3*rho*k term gets subtracted
    i_turbprod_kterm = 0;
    nkey     = nkey + 1;
    // i_catris_kw set to 1 alters the k-omega turb diffusion terms
    i_catris_kw = 0;
    nkey     = nkey + 1;
    // ismincforce overrides normal smin/initvist restart usage
    ismincforce = -1;
    nkey     = nkey + 1;
    // prod2d3dtrace forces Sij used in 2SijSij to be traceless
    prod2d3dtrace = 0.f;
    nkey     = nkey + 1;
    // i_compress_correct adds dilatation-dissipation type compressibility correction
    i_compress_correct = 0;
    nkey     = nkey + 1;
    // les_model determines LES subgrid model to use with ivisc=25
    les_model = 0;
    nkey     = nkey + 1;
    // les_wallscale turns on van Driest type wall scaling of Delta in LES model
    les_wallscale = 0;
    nkey     = nkey + 1;
    // cs_wale is the (non-dynamic) WALE constant for LES
    cs_wale = 0.0f;
    nkey     = nkey + 1;
    // cs_vreman is the (non-dynamic) Vreman constant for LES
    cs_vreman = 0.0f;
    nkey     = nkey + 1;
    // isstrc =1 turns on curvature correction for ivmx=6 or 7
    isstrc = 0;
    nkey     = nkey + 1;
    // sstrc_crc is the constant for isstrc=1
    sstrc_crc = 1.4f;
    nkey     = nkey + 1;
    // isstsf turns on separation fix for ivmx=6 or 7
    isstsf = 0;
    nkey     = nkey + 1;
    // scal_ic is a scaling factor for 2-eqn model BL-type approx ICs
    scal_ic = 5.e6f;
    nkey     = nkey + 1;
    // ifunct will output a PLOT3D function file
    ifunct = 0;
    nkey     = nkey + 1;
    // lowmem_ux=1 allows reversion to lower memory usage of ux array
    lowmem_ux=0;
    nkey     = nkey + 1;
    // SAR - Spalart-Allmaras rotation correction
    isar = 0;
    nkey     = nkey + 1;
    // constant for SAR - Spalart-Allmaras rotation correction
    crot = 2.0f;
    nkey     = nkey + 1;
    // nonlinear constitutive relation for use with linear models
    i_nonlin = 0;
    nkey     = nkey + 1;
    // nonlinear constitutive relation constant
    c_nonlin = 0.3f;
    nkey     = nkey + 1;
    // limiter for nonlinear constitutive relation constant
    snonlin_lim = 1.e-10f;
    nkey     = nkey + 1;
    // isubit_r - write out the subiteration residual for all variables
    isubit_r    = 0;
    nkey     = nkey + 1;
    // i_wilcox98 is used to change Wilcox88 model to Wilcox98
    i_wilcox98 = 0;
    nkey     = nkey + 1;
    // i_wilcox98_chiw is used to turn on/off Wilcox98 k-omega vortex stretching parameter
    i_wilcox98_chiw = 1;
    nkey     = nkey + 1;
    // cmulim limits the abs min computed value of cmu for EASM
    cmulim = .0005f;
    nkey     = nkey + 1;
    // iaxi2planeturb used to force solve of field turb model eqns in j-k plane only
    iaxi2planeturb=0;
    nkey     = nkey + 1;
    // istrongturbdis for strong conservation of diss terms in field turb eqns
    istrongturbdis=0;
    nkey     = nkey + 1;
    // ieasm_type used to change pressure-strain models
    ieasm_type=0;
    nkey     = nkey + 1;
    // ipatch1st=1 forces patching interpolation to be 1st order
    ipatch1st= 0;
    nkey     = nkey + 1;
    // isst2003=1 when ivisc=7 forces SST-2003 model
    isst2003=0;
    nkey     = nkey + 1;
    // issglrrw2012=1 when using ivisc=72 and want SSG/LRR-RSM-w2012
    issglrrw2012=0;
    nkey     = nkey + 1;

    // optional printout to fort.50
    ifort50write=0;
    nkey     = nkey + 1;
    j_ifort50write=1;
    nkey     = nkey + 1;
    i_ifort50write=1;
    nkey     = nkey + 1;
    // i_sas_rsm adds/subtracts SAS-like term to RSM omega eqn
    i_sas_rsm=0;
    nkey     = nkey + 1;
    // iforcev0=1 forces v=0 on update
    iforcev0=0;
    nkey     = nkey + 1;
    // i_saneg=1 uses SA-neg version
    i_saneg=0;
    nkey     = nkey + 1;
    // i_sanoft2=1 uses SA-noft2 version
    i_sanoft2=0;
    nkey     = nkey + 1;
    // i_lam_forcezero=1 forces laminar regions to have zero vist3d
    i_lam_forcezero=0;
    nkey     = nkey + 1;
    // i_specialtop_kmax1001
    i_specialtop_kmax1001=0;
    a_specialtop_kmax1001=0.f;
    xc_specialtop_kmax1001=0.f;
    sig_specialtop_kmax1001=1.f;
    vtp_specialtop_kmax1001=0.f;
    wc_specialtop_kmax1001=0.f;
    fac_specialtop_kmax1001=1.0f;
    cc_specialtop_kmax1001=0.5f;
    xerf_specialtop_kmax1001=100000.f;
    sigerf_specialtop_kmax1001=1.e-9f;
    nkey     = nkey + 10;
    // SA constants
    sa_cw2=0.3f;
    sa_cw3=2.0f;
    sa_cv1=7.1f;
    sa_ct3=1.2f;
    sa_ct4=0.5f;
    sa_cb1=0.1355f;
    sa_cb2=0.622f;
    sa_sigma=2.f/3.f;
    sa_karman=0.41f;
    nkey=nkey + 9;
    // iupdatemean = 0 forces no mean flow update
    iupdatemean=1;
    nkey     = nkey + 1;
    // uub_2034,vvb_2034,wwb_2034 are wall velocities used with bc2034
    uub_2034=0.f;
    vvb_2034=0.f;
    wwb_2034=0.f;
    vvb_xaxisrot2034=0.f;
    nkey     = nkey + 4;
    // i_tauijs adds turbulence via tauij terms
    i_tauijs=0;
    nkey     = nkey + 1;
    // i_qcr2000 turns on QCR2000
    i_qcr2000=0;
    nkey     = nkey + 1;
    // i_qcr2013 turns on QCR2013 (limited)
    i_qcr2013=0;
    nkey     = nkey + 1;
    // i_qcr2013v turns on QCR2013-V
    i_qcr2013v=0;
    nkey     = nkey + 1;
    // i_yapterm=1 adds Yap-like term to get rid of backbending in RSM
    i_yapterm=0;
    nkey     = nkey + 1;


    // check for keyword-driven inputs
    // read(iunit5,1593) inpstr  -- read 210-char line
    {
        FILE* f5 = fortran_get_unit(iunit5);
        std::memset(inpstr, ' ', 210);
        if (fgets(inpstr, 211, f5) == nullptr) {
            std::memset(inpstr, ' ', 210);
        } else {
            int len = (int)std::strlen(inpstr);
            if (len > 0 && inpstr[len-1] == '\n') { inpstr[len-1] = ' '; len--; }
            // pad to 210 with spaces
            for (int ii = len; ii < 210; ii++) inpstr[ii] = ' ';
        }
    }
    if (inpstr[0] == '>') {
        goto label_1000;
    } else if (inpstr[0] == '<') {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120,
            "   ERROR: keyword input must start with a '>'-line in the input file");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    } else {
        ititr = 1;
        // read(inpstr,10)(realval(i),i=1,20)  -- format 20a4
        // Each realval(i) is a 4-char chunk interpreted as float
        // In Fortran, format 20a4 reads 20 groups of 4 chars into real array
        // This is a binary/character transfer: copy 4 bytes per element
        for (int i = 1; i <= 20; i++) {
            std::memcpy(&realval[i], &inpstr[(i-1)*4], 4);
        }
        for (int i = 1; i <= 20; i++) {
            cmn_info.title[i-1] = realval[i];
        }
        return;
    }

label_1000:
    if (iunit11 > 0) {
        FILE* f11 = fortran_get_unit(iunit11);
        fprintf(f11, ">---------------------begin keyword-driven input section --------------------->\n");
    }

label_1100:
    {
        FILE* f5 = fortran_get_unit(iunit5);
        std::memset(inpstr, ' ', 210);
        bool got_eof = false;
        if (fgets(inpstr, 211, f5) == nullptr) {
            got_eof = true;
        } else {
            int len = (int)std::strlen(inpstr);
            if (len > 0 && inpstr[len-1] == '\n') { inpstr[len-1] = ' '; len--; }
            for (int ii = len; ii < 210; ii++) inpstr[ii] = ' ';
        }
        if (got_eof) goto label_1900;
    }
    if (inpstr[0] == '<') {
        if (iunit11 > 0) {
            FILE* f11 = fortran_get_unit(iunit11);
            fprintf(f11, " <----------------------end keyword-driven input section ----------------------<\n");
        }
        // do any order-dependent things here before exiting

        // Check for the existence of the averaging files
        if (iteravg == 2) {
            isthere = false;
            char filename[242];
            std::memset(filename, 0, sizeof(filename));
            // output_dir(1:len_trim(output_dir)) // "/" // "cfl3d_avgq.p3d" or "cfl3d_avg_ruvwp.p3d"
            {
                int ltrim = 240;
                while (ltrim > 0 && output_dir[ltrim-1] == ' ') ltrim--;
                if (ipertavg == 2) {
                    std::snprintf(filename, 242, "%.*s/cfl3d_avg_ruvwp.p3d", ltrim, output_dir);
                } else {
                    std::snprintf(filename, 242, "%.*s/cfl3d_avgq.p3d", ltrim, output_dir);
                }
            }
            // inquire(file=filename, exist=isthere)
            {
                FILE* ftmp = fopen(filename, "r");
                if (ftmp != nullptr) { isthere = true; fclose(ftmp); }
            }
            if (isthere) {
                // open(150, file=filename, form='unformatted', status='old')
                fortran_open_unit(150, filename, "rb");
                // read(150,end=1011,err=1011) idummy
                {
                    FILE* f150 = fortran_get_unit(150);
                    if (fread(&idummy, sizeof(int), 1, f150) != 1) {
                        fortran_close_unit(150);
                        goto label_1012;
                    }
                }
                fortran_close_unit(150);
                goto label_1013;
            } else {
                goto label_1012;
            }
        } else {
            goto label_1013;
        }
label_1012:
        iteravg = 1;
        if (ipertavg == 2) ipertavg = 1;
label_1013:;

        // Check for the existence of the clcd.bin file
        if (iclcd == 2) {
            isthere = false;
            char filename2[242];
            std::memset(filename2, 0, sizeof(filename2));
            {
                int ltrim = 240;
                while (ltrim > 0 && output_dir[ltrim-1] == ' ') ltrim--;
                std::snprintf(filename2, 242, "%.*s/clcd.bin", ltrim, output_dir);
            }
            {
                FILE* ftmp = fopen(filename2, "r");
                if (ftmp != nullptr) { isthere = true; fclose(ftmp); }
            }
            if (isthere) {
                fortran_open_unit(150, filename2, "rb");
                {
                    FILE* f150 = fortran_get_unit(150);
                    if (fread(&idummy1, sizeof(int), 1, f150) != 1 ||
                        fread(&idummy2, sizeof(int), 1, f150) != 1 ||
                        fread(&idummy3, sizeof(int), 1, f150) != 1) {
                        fortran_close_unit(150);
                        goto label_1015;
                    }
                }
                fortran_close_unit(150);
                goto label_1016;
            } else {
                goto label_1015;
            }
        } else {
            goto label_1016;
        }
label_1015:
        iclcd = 1;
label_1016:;


        // the variables turbintensity_inf_percent and eddy_visc_inf
        // do not go to the rest of the code; instead, if set, they
        // modify tur10(2) and tur10(1), respectively
        if (turbintensity_inf_percent >= 0.f) {
            if (ivmx == 72) {
                if (iunit11 > 0) {
                    FILE* f11 = fortran_get_unit(iunit11);
                    fprintf(f11, " turbintensity_inf_percent always overrides tur10,20,30 for ivmx=72\n");
                }
                zkinf = 1.5f * xmach*xmach * (turbintensity_inf_percent/100.f) * (turbintensity_inf_percent/100.f);
                tur10[0] = -2.f/3.f * zkinf;
                tur10[1] = -2.f/3.f * zkinf;
                tur10[2] = -2.f/3.f * zkinf;
            } else {
                if (iunit11 > 0) {
                    FILE* f11 = fortran_get_unit(iunit11);
                    fprintf(f11, " turbintensity_inf_percent always overrides tur20\n");
                }
                tur10[1] = 1.5f * xmach*xmach * (turbintensity_inf_percent/100.f) * (turbintensity_inf_percent/100.f);
            }
        }
        if (eddy_visc_inf >= 0.f) {
            if (ivmx == 72) {
                if (iunit11 > 0) {
                    FILE* f11 = fortran_get_unit(iunit11);
                    fprintf(f11, " eddy_visc_inf always overrides tur70 for ivmx=72\n");
                }
            } else {
                if (iunit11 > 0) {
                    FILE* f11 = fortran_get_unit(iunit11);
                    fprintf(f11, " eddy_visc_inf always overrides tur10\n");
                }
            }
            if (ivmx == 4) {
                tur10[0] = eddy_visc_inf / 0.09f;
            } else if (ivmx == 6 || ivmx == 7 || ivmx == 30 || ivmx == 40) {
                tur10[0] = tur10[1] / eddy_visc_inf;
            } else if (ivmx == 8 || ivmx == 14) {
                tur10[0] = 0.0895f * tur10[1] / eddy_visc_inf;
            } else if (ivmx == 12) {
                tur10[0] = 0.081f * tur10[1] / eddy_visc_inf;
            } else if (ivmx == 10 || ivmx == 15) {
                tur10[0] = 0.09f * tur10[1]*tur10[1] / eddy_visc_inf;
            } else if (ivmx == 9 || ivmx == 13) {
                tur10[0] = 0.0885f * tur10[1]*tur10[1] / eddy_visc_inf;
            } else if (ivmx == 11) {
                tur10[0] = 0.081f * tur10[1]*tur10[1] / eddy_visc_inf;
            } else if (ivmx == 5) {
                // for SA model, need to iterate to find nuwiggle:
                tur10[0] = eddy_visc_inf;
                for (n = 1; n <= 100; n++) {
                    errornum = tur10[0]*tur10[0]*tur10[0]*tur10[0]
                             - eddy_visc_inf*tur10[0]*tur10[0]*tur10[0]
                             - 357.911f*eddy_visc_inf;
                    errorden = 4.f*tur10[0]*tur10[0]*tur10[0]
                             - 3.f*eddy_visc_inf*tur10[0]*tur10[0];
                    error_val = errornum / errorden;
                    if (std::abs(error_val) < 1.e-6f) goto label_105;
                    tur10[0] = tur10[0] - error_val;
                }
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "stopping...Newton iteration using eddy_visc_inf not converged in readkey");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
label_105:;
            } else if (ivmx == 72) {
                tur10[6] = -3.f/2.f * tur10[0] / eddy_visc_inf;
            }
        }
        return;
    }

    // if line starts with '/' or ' ', skip it
    if (inpstr[0] == '/' || inpstr[0] == ' ') goto label_1100;

    // echo the keyword character strings
    npos = 1;
    lc2  = 0;
    {
        int neg1 = -1;
        parser_ns::parser(inpstr, npos, lc1, lc2, lcl, neg1);
    }
    lclw = std::min(lcl, 80);
    if (iunit11 > 0) {
        FILE* f11 = fortran_get_unit(iunit11);
        // write(iunit11,1594) inpstr(1:lclw)  -- format(a)
        fprintf(f11, "%.*s\n", lclw, inpstr);
    }

    npos = 1;
    lc2  = 0;
    lcl  = 210;
    {
        int pos1 = 1;
        parser_ns::parser(inpstr, npos, lc1, lc2, lcl, pos1);
    }


    // keyword comparison: inpstr(lc1:lc2) is the keyword token
    {
        int klen = lc2 - lc1 + 1;
        std::string kw(inpstr + lc1 - 1, klen);
        int rest_start = lc2; // 0-based index after keyword (lc2+1 in Fortran = lc2 in C 0-based)

        if (kw == "gamma") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            gamma = realval[1];
        } else if (kw == "pr") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            pr = realval[1];
        } else if (kw == "prt") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            prt = realval[1];
        } else if (kw == "cbar") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cbar = realval[1];
        } else if (kw == "atol") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            atol = realval[1];
        } else if (kw == "xmach_img") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            xmach_img = realval[1];
        } else if (kw == "alpha_img") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            alpha_img = realval[1];
        } else if (kw == "beta_img") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            beta_img = realval[1];
        } else if (kw == "reue_img") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            reue_img = realval[1];
        } else if (kw == "tinf_img") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tinf_img = realval[1];
        } else if (kw == "geom_img") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            geom_img = realval[1];
        } else if (kw == "icgns") {
            // CGNS grid I/O is enabled (linked against libcgns); accept icgns=1.
            sscanf(inpstr + rest_start, "%d", &icgns);
        } else if (kw == "cprec") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cprec = realval[1];
        } else if (kw == "uref") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            uref = realval[1];
        } else if (kw == "avn") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            avn = realval[1];
        } else if (kw == "cltarg") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cltarg = realval[1];
            if (cltarg != 99999.0f) { ialphit = 1; }
        } else if (kw == "rlxalph") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            rlxalph = realval[1];
        } else if (kw == "nsubturb") {
            sscanf(inpstr + rest_start, "%d", &nsubturb);
        } else if (kw == "cflturb") {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "stopping...keyword cflturb no longer allowed... use cflturb1, 2, etc instead");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        } else if (kw == "cflturb1") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cflturb[0] = std::abs(realval[1]);
        } else if (kw == "cflturb2") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cflturb[1] = std::abs(realval[1]);
        } else if (kw == "cflturb3") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cflturb[2] = std::abs(realval[1]);
        } else if (kw == "cflturb4") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cflturb[3] = std::abs(realval[1]);
        } else if (kw == "cflturb5") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cflturb[4] = std::abs(realval[1]);
        } else if (kw == "cflturb6") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cflturb[5] = std::abs(realval[1]);
        } else if (kw == "cflturb7") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cflturb[6] = std::abs(realval[1]);
        } else if (kw == "epsa_r") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            epsa_r = realval[1];
            epsa_l = 2.f * epsa_r;
        } else if (kw == "nfreeze") {
            sscanf(inpstr + rest_start, "%d", &nfreeze);
        } else if (kw == "ivolint") {
            sscanf(inpstr + rest_start, "%d", &ivolint);
        } else if (kw == "idef_ss") {
            sscanf(inpstr + rest_start, "%d", &idef_ss);
        } else if (kw == "ibin") {
            sscanf(inpstr + rest_start, "%d", &ibin);
        } else if (kw == "iblnk") {
            sscanf(inpstr + rest_start, "%d", &iblnk);
        } else if (kw == "iblnkfr") {
            sscanf(inpstr + rest_start, "%d", &iblnkfr);
        } else if (kw == "ip3dgrad") {
            sscanf(inpstr + rest_start, "%d", &ip3dgrad);
        } else if (kw == "surf_img") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            surf_img = realval[1];
            if (idef_ss == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "stopping...must have idef_ss = 1 if surf_img .ne. 0.");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        } else if (kw == "memadd") {
            sscanf(inpstr + rest_start, "%d", &memadd);
        } else if (kw == "memaddi") {
            sscanf(inpstr + rest_start, "%d", &memaddi);
        } else if (kw == "negvol") {
            sscanf(inpstr + rest_start, "%d", &negvol);
        } else if (kw == "meshdef") {
            sscanf(inpstr + rest_start, "%d", &meshdef);
        } else if (kw == "edvislim") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            edvislim = realval[1];
        } else if (kw == "iturbprod") {
            sscanf(inpstr + rest_start, "%d", &iturbprod);
        } else if (kw == "irghost") {
            sscanf(inpstr + rest_start, "%d", &irghost);
        } else if (kw == "iwghost") {
            sscanf(inpstr + rest_start, "%d", &iwghost);
        } else if (kw == "dalim") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            dalim = realval[1];
        } else if (kw == "icycupdt") {
            sscanf(inpstr + rest_start, "%d", &icycupdt);
        } else if (kw == "noninflag") {
            sscanf(inpstr + rest_start, "%d", &noninflag);
        } else if (kw == "xcentrot") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            xcentrot = realval[1];
        } else if (kw == "ycentrot") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            ycentrot = realval[1];
        } else if (kw == "zcentrot") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            zcentrot = realval[1];
        } else if (kw == "xrotrate") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            xrotrate = realval[1];
        } else if (kw == "yrotrate") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            yrotrate = realval[1];
        } else if (kw == "zrotrate") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            zrotrate = realval[1];
        } else if (kw == "xrotrate_img") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            xrotrate_img = realval[1];
        } else if (kw == "yrotrate_img") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            yrotrate_img = realval[1];
        } else if (kw == "zrotrate_img") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            zrotrate_img = realval[1];
        } else if (kw == "itime2read") {
            sscanf(inpstr + rest_start, "%d", &itime2read);
        } else if (kw == "itaturb") {
            sscanf(inpstr + rest_start, "%d", &itaturb);
        } else if (kw == "ides") {
            sscanf(inpstr + rest_start, "%d", &cmn_des.ides);
        } else if (kw == "cdes") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cmn_des.cdes = realval[1];
        } else if (kw == "cddes") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cmn_des.cddes = realval[1];
        } else if (kw == "iteravg") {
            sscanf(inpstr + rest_start, "%d", &iteravg);
        } else if (kw == "tur10") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tur10[0] = realval[1];
        } else if (kw == "tur20") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tur10[1] = realval[1];
        } else if (kw == "tur30") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tur10[2] = realval[1];
        } else if (kw == "tur40") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tur10[3] = realval[1];
        } else if (kw == "tur50") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tur10[4] = realval[1];
        } else if (kw == "tur60") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tur10[5] = realval[1];
        } else if (kw == "tur70") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tur10[6] = realval[1];
        } else if (kw == "tur1cut") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tur1cut = realval[1];
        } else if (kw == "roll_angle") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            roll_angle = realval[1];
        } else if (kw == "ikoprod") {
            sscanf(inpstr + rest_start, "%d", &ikoprod);
        } else if (kw == "isstdenom") {
            sscanf(inpstr + rest_start, "%d", &isstdenom);
        } else if (kw == "pklimterm") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            pklimterm = realval[1];
        } else if (kw == "iaxi2plane") {
            sscanf(inpstr + rest_start, "%d", &iaxi2plane);
        } else if (kw == "iturbord") {
            sscanf(inpstr + rest_start, "%d", &iturbord);
        } else if (kw == "tur2cut") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tur2cut = realval[1];
        } else if (kw == "tur1cutlev") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tur1cutlev = realval[1];
        } else if (kw == "tur2cutlev") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tur2cutlev = realval[1];
        } else if (kw == "ifullns") {
            sscanf(inpstr + rest_start, "%d", &ifullns);
        } else if (kw == "ibeta8kzeta") {
            sscanf(inpstr + rest_start, "%d", &ibeta8kzeta);
        } else if (kw == "isarc2d") {
            sscanf(inpstr + rest_start, "%d", &isarc2d);
        } else if (kw == "isarc3d") {
            sscanf(inpstr + rest_start, "%d", &isarc3d);
        } else if (kw == "sarccr3") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sarccr3 = realval[1];
        } else if (kw == "ieasmcc2d") {
            sscanf(inpstr + rest_start, "%d", &ieasmcc2d);
        } else if (kw == "ipertavg") {
            sscanf(inpstr + rest_start, "%d", &ipertavg);
            if (ipertavg != 0) iteravg = ipertavg;
        } else if (kw == "icoarsemovie") {
            sscanf(inpstr + rest_start, "%d", &icoarsemovie);
        } else if (kw == "i2dmovie") {
            sscanf(inpstr + rest_start, "%d", &i2dmovie);
        } else if (kw == "iclcd") {
            sscanf(inpstr + rest_start, "%d", &iclcd);
        } else if (kw == "iskip_blocks") {
            sscanf(inpstr + rest_start, "%d", &iskip_blocks);
        } else if (kw == "cfltauMax") {
            sscanf(inpstr + rest_start, "%f", &cfltauMax);
        } else if (kw == "cfltau0") {
            sscanf(inpstr + rest_start, "%f", &cfltau0);
        } else if (kw == "irbtrim") {
            sscanf(inpstr + rest_start, "%d", &irbtrim);
        } else if (kw == "irigb") {
            sscanf(inpstr + rest_start, "%d", &irigb);
        } else if (kw == "greflrb") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            greflrb = realval[1];
        } else if (kw == "tmass") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            tmass = realval[1];
        } else if (kw == "yinert") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            yinert = realval[1];
        } else if (kw == "gaccel") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            gaccel = realval[1];
        } else if (kw == "relax") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            relax = realval[1];
        } else if (kw == "itrminc") {
            sscanf(inpstr + rest_start, "%d", &itrminc);
        } else if (kw == "dclda") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            dclda = realval[1];
        } else if (kw == "dcldd") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            dcldd = realval[1];
        } else if (kw == "dcmda") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            dcmda = realval[1];
        } else if (kw == "dcmdd") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            dcmdd = realval[1];
        } else if (kw == "ndgrd") {
            sscanf(inpstr + rest_start, "%d", &ndgrd);
        } else if (kw == "ndwrt") {
            sscanf(inpstr + rest_start, "%d", &ndwrt);
        } else if (kw == "i_bsl") {
            sscanf(inpstr + rest_start, "%d", &i_bsl);
        } else if (kw == "keepambient") {
            sscanf(inpstr + rest_start, "%d", &keepambient);
        } else if (kw == "re_thetat0") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            re_thetat0 = realval[1];
        } else if (kw == "turbintensity_inf_percent") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            turbintensity_inf_percent = realval[1];
        } else if (kw == "eddy_visc_inf") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            eddy_visc_inf = realval[1];
        } else if (kw == "cs_smagorinsky") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cs_smagorinsky = realval[1];
        } else if (kw == "xdir_only_source") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            xdir_only_source = realval[1];
        } else if (kw == "randomize") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            randomize = std::abs(realval[1]);
        } else if (kw == "iexact_trunc") {
            sscanf(inpstr + rest_start, "%d", &iexact_trunc);
        } else if (kw == "iexact_disc") {
            sscanf(inpstr + rest_start, "%d", &iexact_disc);
        } else if (kw == "iexact_ring") {
            sscanf(inpstr + rest_start, "%d", &iexact_ring);
        } else if (kw == "i_wilcox06") {
            sscanf(inpstr + rest_start, "%d", &i_wilcox06);
        } else if (kw == "i_wilcox06_chiw") {
            sscanf(inpstr + rest_start, "%d", &i_wilcox06_chiw);
        } else if (kw == "i_turbprod_kterm") {
            sscanf(inpstr + rest_start, "%d", &i_turbprod_kterm);
        } else if (kw == "i_catris_kw") {
            sscanf(inpstr + rest_start, "%d", &i_catris_kw);
        } else if (kw == "ismincforce") {
            sscanf(inpstr + rest_start, "%d", &ismincforce);
        } else if (kw == "prod2d3dtrace") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            prod2d3dtrace = std::abs(realval[1]);
        } else if (kw == "i_compress_correct") {
            sscanf(inpstr + rest_start, "%d", &i_compress_correct);
        } else if (kw == "les_model") {
            sscanf(inpstr + rest_start, "%d", &les_model);
        } else if (kw == "les_wallscale") {
            sscanf(inpstr + rest_start, "%d", &les_wallscale);
        } else if (kw == "cs_wale") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cs_wale = std::abs(realval[1]);
        } else if (kw == "cs_vreman") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cs_vreman = std::abs(realval[1]);
        } else if (kw == "isstrc") {
            sscanf(inpstr + rest_start, "%d", &isstrc);
        } else if (kw == "sstrc_crc") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sstrc_crc = std::abs(realval[1]);
        } else if (kw == "isstsf") {
            sscanf(inpstr + rest_start, "%d", &isstsf);
        } else if (kw == "scal_ic") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            scal_ic = std::abs(realval[1]);
        } else if (kw == "ifunct") {
            sscanf(inpstr + rest_start, "%d", &ifunct);
        } else if (kw == "lowmem_ux") {
            sscanf(inpstr + rest_start, "%d", &lowmem_ux);
        } else if (kw == "isar") {
            sscanf(inpstr + rest_start, "%d", &isar);
        } else if (kw == "crot") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            crot = std::abs(realval[1]);
        } else if (kw == "i_nonlin") {
            sscanf(inpstr + rest_start, "%d", &i_nonlin);
        } else if (kw == "c_nonlin") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            c_nonlin = std::abs(realval[1]);
        } else if (kw == "snonlin_lim") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            snonlin_lim = std::abs(realval[1]);
        } else if (kw == "isubit_r") {
            sscanf(inpstr + rest_start, "%d", &isubit_r);
        } else if (kw == "i_wilcox98") {
            sscanf(inpstr + rest_start, "%d", &i_wilcox98);
        } else if (kw == "i_wilcox98_chiw") {
            sscanf(inpstr + rest_start, "%d", &i_wilcox98_chiw);
        } else if (kw == "cmulim") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cmulim = std::abs(realval[1]);
        } else if (kw == "iaxi2planeturb") {
            sscanf(inpstr + rest_start, "%d", &iaxi2planeturb);
        } else if (kw == "istrongturbdis") {
            sscanf(inpstr + rest_start, "%d", &istrongturbdis);
        } else if (kw == "ieasm_type") {
            sscanf(inpstr + rest_start, "%d", &ieasm_type);
        } else if (kw == "ipatch1st") {
            sscanf(inpstr + rest_start, "%d", &ipatch1st);
        } else if (kw == "isst2003") {
            sscanf(inpstr + rest_start, "%d", &isst2003);
        } else if (kw == "issglrrw2012") {
            sscanf(inpstr + rest_start, "%d", &issglrrw2012);
        } else if (kw == "ifort50write") {
            sscanf(inpstr + rest_start, "%d", &ifort50write);
        } else if (kw == "j_ifort50write") {
            sscanf(inpstr + rest_start, "%d", &j_ifort50write);
        } else if (kw == "i_ifort50write") {
            sscanf(inpstr + rest_start, "%d", &i_ifort50write);
        } else if (kw == "i_sas_rsm") {
            sscanf(inpstr + rest_start, "%d", &i_sas_rsm);
        } else if (kw == "iforcev0") {
            sscanf(inpstr + rest_start, "%d", &iforcev0);
        } else if (kw == "i_saneg") {
            sscanf(inpstr + rest_start, "%d", &i_saneg);
        } else if (kw == "i_sanoft2") {
            sscanf(inpstr + rest_start, "%d", &i_sanoft2);
        } else if (kw == "i_lam_forcezero") {
            sscanf(inpstr + rest_start, "%d", &i_lam_forcezero);
        } else if (kw == "i_specialtop_kmax1001") {
            sscanf(inpstr + rest_start, "%d", &i_specialtop_kmax1001);
        } else if (kw == "a_specialtop_kmax1001") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            a_specialtop_kmax1001 = realval[1];
        } else if (kw == "xc_specialtop_kmax1001") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            xc_specialtop_kmax1001 = realval[1];
        } else if (kw == "sig_specialtop_kmax1001") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sig_specialtop_kmax1001 = realval[1];
        } else if (kw == "vtp_specialtop_kmax1001") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            vtp_specialtop_kmax1001 = realval[1];
        } else if (kw == "wc_specialtop_kmax1001") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            wc_specialtop_kmax1001 = realval[1];
        } else if (kw == "fac_specialtop_kmax1001") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            fac_specialtop_kmax1001 = realval[1];
        } else if (kw == "cc_specialtop_kmax1001") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            cc_specialtop_kmax1001 = realval[1];
        } else if (kw == "xerf_specialtop_kmax1001") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            xerf_specialtop_kmax1001 = realval[1];
        } else if (kw == "sigerf_specialtop_kmax1001") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sigerf_specialtop_kmax1001 = realval[1];
        } else if (kw == "sa_cw2") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sa_cw2 = realval[1];
        } else if (kw == "sa_cw3") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sa_cw3 = realval[1];
        } else if (kw == "sa_cv1") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sa_cv1 = realval[1];
        } else if (kw == "sa_ct3") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sa_ct3 = realval[1];
        } else if (kw == "sa_ct4") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sa_ct4 = realval[1];
        } else if (kw == "sa_cb1") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sa_cb1 = realval[1];
        } else if (kw == "sa_cb2") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sa_cb2 = realval[1];
        } else if (kw == "sa_sigma") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sa_sigma = realval[1];
        } else if (kw == "sa_karman") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            sa_karman = realval[1];
        } else if (kw == "iupdatemean") {
            sscanf(inpstr + rest_start, "%d", &iupdatemean);
        } else if (kw == "uub_2034") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            uub_2034 = realval[1];
        } else if (kw == "vvb_2034") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            vvb_2034 = realval[1];
        } else if (kw == "wwb_2034") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            wwb_2034 = realval[1];
        } else if (kw == "vvb_xaxisrot2034") {
            sscanf(inpstr + rest_start, "%f", &realval[1]);
            vvb_xaxisrot2034 = realval[1];
        } else if (kw == "i_tauijs") {
            sscanf(inpstr + rest_start, "%d", &i_tauijs);
        } else if (kw == "i_qcr2000") {
            sscanf(inpstr + rest_start, "%d", &i_qcr2000);
        } else if (kw == "i_qcr2013") {
            sscanf(inpstr + rest_start, "%d", &i_qcr2013);
        } else if (kw == "i_qcr2013v") {
            sscanf(inpstr + rest_start, "%d", &i_qcr2013v);
        } else if (kw == "i_yapterm") {
            sscanf(inpstr + rest_start, "%d", &i_yapterm);
        } else {
            // unknown keyword
            if (iunit11 > 0) {
                FILE* f11 = fortran_get_unit(iunit11);
                fprintf(f11, "*** STOPPING: The keyword above is not supported. ***\n");
            }
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    } // end keyword block

    // Set trim data
    tp_arr[0] = dclda;
    tp_arr[4] = dcldd;
    tp_arr[1] = dcmda;
    tp_arr[5] = dcmdd;

    goto label_1100;

label_1900:
    // Make sure we are running ddes if cddes is set
    if (std::abs(cmn_des.cddes - 0.975f) > 0.001f) cmn_des.ides = 3;
    if ((float)cmn_des.cddes > 0.999f) cmn_des.cddes = 0.999f;
    if (cmn_des.ides == 3) {
        printf(" Running ides with cddes =  %g\n", (double)cmn_des.cddes);
    }

    nou(1) = std::min(nou(1)+1, ibufdim);
    std::snprintf(bou(nou(1),1), 120,
        "   ERROR: There was no '<'-line in the input file");
    nou(1) = std::min(nou(1)+1, ibufdim);
    std::snprintf(bou(nou(1),1), 120,
        "          to exit the keyword-driven input section.");
    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);

    return;
}

// trnsfr_vals: declared in readkey_ns (readkey.h), forwards to trnsfr_vals_ns
void trnsfr_vals(FortranArray1DRef<double> work, int& mwork, FortranArray1DRef<int> iwork, int& mworki, int& maxgr, int& maxbl, int& maxseg, int& nbli, int& nblock, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray1DRef<int> iovrlp, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblcg, FortranArray2DRef<int> idegg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> ieg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> klamlog, FortranArray1DRef<int> ilamlog, FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamhig, FortranArray2DRef<int> mit, FortranArray2DRef<int> iwfg, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> iforce, int& lfgm, int& ihstry, int& ninter, int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg, FortranArray2DRef<int> idiagg, int& nplot3d, int& nprint, int& ncs, int& ip3dsurf, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx, FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz, FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr, FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0, FortranArray1DRef<double> time2, FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae, FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae, FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass, FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0, FortranArray2DRef<double> damp, FortranArray3DRef<double> perturb, FortranArray2DRef<double> aesrfdat, FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray2DRef<double> rkap0g, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, int& mxbli, FortranArray2DRef<int> nblk, FortranArray3DRef<int> limblk, FortranArray1DRef<int> nblon, FortranArray3DRef<int> isva, FortranArray1DRef<char[80]> bcfiles, FortranArray2DRef<int> icsinfo, int& mxbcfil, int& maxcs, int& nmds, int& maxaes, FortranArray2DRef<int> inpl3d, FortranArray2DRef<int> inpr, int& nplots, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg, FortranArray2DRef<double> xorgae0, FortranArray2DRef<double> yorgae0, FortranArray2DRef<double> zorgae0, FortranArray2DRef<int> icouple, int& iprnsurf)
{
    trnsfr_vals_ns::trnsfr_vals(work, mwork, iwork, mworki, maxgr, maxbl, maxseg, nbli, nblock, levelg, igridg, iovrlp, nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim, iflimg, ifdsg, iviscg, jdimg, kdimg, idimg, nblcg, idegg, isg, jsg, ksg, ieg, jeg, keg, jlamlog, klamlog, ilamlog, ilamhig, jlamhig, klamhig, mit, iwfg, iadvance, iforce, lfgm, ihstry, ninter, ngrid, ncgg, nblg, iemg, inewgg, idiagg, nplot3d, nprint, ncs, ip3dsurf, mblk2nd, utrans, vtrans, wtrans, omegax, omegay, omegaz, xorig, yorig, zorig, dxmx, dymx, dzmx, dthxmx, dthymx, dthzmx, thetax, thetay, thetaz, rfreqt, rfreqr, xorig0, yorig0, zorig0, time2, thetaxl, thetayl, thetazl, itrans, irotat, idefrm, utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae, xorgae, yorgae, zorgae, thtxae, thtyae, thtzae, rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf, kcsi, kcsf, freq, gmass, x0, gf0, damp, perturb, aesrfdat, bcvali, bcvalj, bcvalk, ibcinfo, jbcinfo, kbcinfo, rkap0g, bcfilei, bcfilej, bcfilek, mxbli, nblk, limblk, nblon, isva, bcfiles, icsinfo, mxbcfil, maxcs, nmds, maxaes, inpl3d, inpr, nplots, nou, bou, nbuf, ibufdim, iskip, jskip, kskip, nsegdfrm, idfrmseg, iaesurf, maxsegdg, xorgae0, yorgae0, zorgae0, icouple, iprnsurf);
}

} // namespace readkey_ns
