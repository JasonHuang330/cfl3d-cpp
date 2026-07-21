// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "readkey.h"
#include "parser.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>

namespace readkey_ns {

void readkey(int& ititr, int& myid, int& ibufdim, int& nbuf,
             FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou,
             int& iunit11, int& ierrflg)
{
    // Local variables
    char inpstr[210];
    float realval[20];
    bool isthere;
    float turbintensity_inf_percent;
    float eddy_visc_inf;
    float dclda, dcldd, dcmda, dcmdd;
    int idummy, idummy1, idummy2, idummy3;
    int npos, lc1, lc2, lcl, lclw;
    int i, n;
    float errornum, errorden, error_val;
    float zkinf;

    // Helper lambdas for keyword parsing (defined early to avoid goto-over-init issues)
    auto read_float_kw = [&](float& val) {
        int start = lc2 - 1; if (start < 0) start = 0; if (start >= 210) { val = 0.f; return; }
        char tmp[211]; int len = 210 - start; memcpy(tmp, inpstr + start, len); tmp[len] = '\0';
        sscanf(tmp, "%f", &val);
    };
    auto read_int_kw = [&](int& val) {
        int start = lc2 - 1; if (start < 0) start = 0; if (start >= 210) { val = 0; return; }
        char tmp[211]; int len = 210 - start; memcpy(tmp, inpstr + start, len); tmp[len] = '\0';
        sscanf(tmp, "%d", &val);
    };
    auto kw_eq = [&](const char* kw) -> bool {
        int klen = (int)strlen(kw); int flen = lc2 - lc1 + 1;
        if (flen != klen) return false;
        return (strncmp(inpstr + lc1 - 1, kw, klen) == 0);
    };

    // Common block references (shorthand aliases)
    auto& iunit5         = cmn_unit5.iunit5;
    auto& time           = cmn_unst.time;
    auto& cfltau         = cmn_unst.cfltau;
    auto& ntstep         = cmn_unst.ntstep;
    auto& ita            = cmn_unst.ita;
    auto& iunst          = cmn_unst.iunst;
    auto& cfltau0        = cmn_unst.cfltau0;
    auto& cfltauMax      = cmn_unst.cfltaumax;
    auto& cdes           = cmn_des.cdes;
    auto& ides           = cmn_des.ides;
    auto& cddes          = cmn_des.cddes;
    auto& xnumavg        = cmn_avgdata.xnumavg;
    auto& iteravg        = cmn_avgdata.iteravg;
    auto& xnumavg2       = cmn_avgdata.xnumavg2;
    auto& ipertavg       = cmn_avgdata.ipertavg;
    auto& iclcd          = cmn_avgdata.iclcd;
    auto& isubit_r       = cmn_avgdata.isubit_r;
    auto& movie          = cmn_moov.movie;
    auto& nframes        = cmn_moov.nframes;
    auto& icall1         = cmn_moov.icall1;
    auto& lhdr           = cmn_moov.lhdr;
    auto& icoarsemovie   = cmn_moov.icoarsemovie;
    auto& i2dmovie       = cmn_moov.i2dmovie;
    auto& lq2avg         = cmn_ginfo2.lq2avg;
    auto& iskip_blocks   = cmn_ginfo2.iskip_blocks;
    // inc_2d and inc_coarse accessed via cmn_ginfo2 directly
    auto& gamma          = cmn_fluid.gamma;
    auto& gm1            = cmn_fluid.gm1;
    auto& gp1            = cmn_fluid.gp1;
    auto& gm1g           = cmn_fluid.gm1g;
    auto& gp1g           = cmn_fluid.gp1g;
    auto& ggm1           = cmn_fluid.ggm1;
    auto& pr             = cmn_fluid2.pr;
    auto& prt            = cmn_fluid2.prt;
    auto& cbar           = cmn_fluid2.cbar;
    // info common
    auto& xmach          = cmn_info.xmach;
    auto& alpha          = cmn_info.alpha;
    auto& beta           = cmn_info.beta;
    auto& dt             = cmn_info.dt;
    auto& fmax           = cmn_info.fmax;
    auto& nit            = cmn_info.nit;
    auto& ntt            = cmn_info.ntt;
    auto& nitfo          = cmn_info.nitfo;
    auto& iflagts        = cmn_info.iflagts;
    auto& nres           = cmn_info.nres;
    auto& mgflag         = cmn_info.mgflag;
    auto& iconsf         = cmn_info.iconsf;
    auto& mseq           = cmn_info.mseq;
    auto& ngam           = cmn_info.ngam;
    auto& iipv           = cmn_info.iipv;
    auto& iexp           = cmn_zero.iexp;
    auto& atol           = cmn_singular.atol;
    auto& xmach_img      = cmn_complx.xmach_img;
    auto& alpha_img      = cmn_complx.alpha_img;
    auto& beta_img       = cmn_complx.beta_img;
    auto& reue_img       = cmn_complx.reue_img;
    auto& tinf_img       = cmn_complx.tinf_img;
    auto& geom_img       = cmn_complx.geom_img;
    auto& surf_img       = cmn_complx.surf_img;
    auto& xrotrate_img   = cmn_complx.xrotrate_img;
    auto& yrotrate_img   = cmn_complx.yrotrate_img;
    auto& zrotrate_img   = cmn_complx.zrotrate_img;
    auto& icgns          = cmn_cgns.icgns;
    auto& iccg           = cmn_cgns.iccg;
    auto& ibase          = cmn_cgns.ibase;
    auto& nzones         = cmn_cgns.nzones;
    auto& nsoluse        = cmn_cgns.nsoluse;
    auto& irind          = cmn_cgns.irind;
    auto& jrind          = cmn_cgns.jrind;
    auto& krind          = cmn_cgns.krind;
    auto& cprec          = cmn_precond.cprec;
    auto& uref           = cmn_precond.uref;
    auto& avn            = cmn_precond.avn;
    // alphait: ialphit, cltarg, rlxalph(=resupdt), dalim, dalpha(=dcldal), icycupdt
    auto& cltarg         = cmn_alphait.cltarg;
    auto& rlxalph        = cmn_alphait.resupdt;
    auto& dalim          = cmn_alphait.dalim;
    auto& ialphit        = cmn_alphait.ialphit;
    auto& icycupdt       = cmn_alphait.icycupdt;
    // turbconv
    auto& edvislim       = cmn_turbconv.edvislim;
    auto& iturbprod      = cmn_turbconv.iturbprod;
    auto& nsubturb       = cmn_turbconv.nsubturb;
    auto& nfreeze        = cmn_turbconv.nfreeze;
    auto& iwarneddy      = cmn_turbconv.iwarneddy;
    auto& itime2read     = cmn_turbconv.itime2read;
    auto& itaturb        = cmn_turbconv.itaturb;
    auto& tur1cut        = cmn_turbconv.tur1cut;
    auto& tur2cut        = cmn_turbconv.tur2cut;
    auto& iturbord       = cmn_turbconv.iturbord;
    auto& tur1cutlev     = cmn_turbconv.tur1cutlev;
    auto& tur2cutlev     = cmn_turbconv.tur2cutlev;
    auto& epsa_l         = cmn_entfix.epsa_l;
    auto& epsa_r         = cmn_entfix.epsa_r;
    auto& nkey           = cmn_key.nkey;
    auto& ivolint        = cmn_is_blockbc.ivolint;
    auto& idef_ss        = cmn_elastic_ss.idef_ss;
    auto& ibin           = cmn_bin.ibin;
    auto& iblnk          = cmn_bin.iblnk;
    auto& iblnkfr        = cmn_bin.iblnkfr;
    auto& ip3dgrad       = cmn_bin.ip3dgrad;
    auto& memadd         = cmn_memory.memadd;
    auto& memaddi        = cmn_memory.memaddi;
    auto& beta1          = cmn_deformz.beta1;
    auto& beta2          = cmn_deformz.beta2;
    auto& alpha1         = cmn_deformz.alpha1;
    auto& alpha2         = cmn_deformz.alpha2;
    auto& isktyp         = cmn_deformz.isktyp;
    auto& negvol         = cmn_deformz.negvol;
    auto& meshdef        = cmn_deformz.meshdef;
    auto& nsprgit        = cmn_deformz.nsprgit;
    auto& ndgrd          = cmn_deformz.ndgrd;
    auto& ndwrt          = cmn_deformz.ndwrt;
    auto& irghost        = cmn_ghost.irghost;
    auto& iwghost        = cmn_ghost.iwghost;
    auto& xcentrot       = cmn_noninertial.xcentrot;
    auto& ycentrot       = cmn_noninertial.ycentrot;
    auto& zcentrot       = cmn_noninertial.zcentrot;
    auto& xrotrate       = cmn_noninertial.xrotrate;
    auto& yrotrate       = cmn_noninertial.yrotrate;
    auto& zrotrate       = cmn_noninertial.zrotrate;
    auto& noninflag      = cmn_noninertial.noninflag;
    auto& roll_angle     = cmn_gridtrans.roll_angle;
    auto& ikoprod        = cmn_konew.ikoprod;
    auto& isstdenom      = cmn_konew.isstdenom;
    auto& pklimterm      = cmn_konew.pklimterm;
    auto& ibeta8kzeta    = cmn_konew.ibeta8kzeta;
    auto& i_bsl          = cmn_konew.i_bsl;
    auto& keepambient    = cmn_konew.keepambient;
    auto& re_thetat0     = cmn_konew.re_thetat0;
    auto& i_wilcox06     = cmn_konew.i_wilcox06;
    auto& i_wilcox06_chiw= cmn_konew.i_wilcox06_chiw;
    auto& i_turbprod_kterm = cmn_konew.i_turbprod_kterm;
    auto& i_catris_kw    = cmn_konew.i_catris_kw;
    auto& prod2d3dtrace  = cmn_konew.prod2d3dtrace;
    auto& i_compress_correct = cmn_konew.i_compress_correct;
    auto& isstsf         = cmn_konew.isstsf;
    auto& i_wilcox98     = cmn_konew.i_wilcox98;
    auto& i_wilcox98_chiw= cmn_konew.i_wilcox98_chiw;
    auto& isst2003       = cmn_konew.isst2003;
    auto& iaxi2plane     = cmn_axisym.iaxi2plane;
    auto& iaxi2planeturb = cmn_axisym.iaxi2planeturb;
    auto& istrongturbdis = cmn_axisym.istrongturbdis;
    auto& iforcev0       = cmn_axisym.iforcev0;
    auto& ifullns        = cmn_fullns.ifullns;
    auto& ivmx           = cmn_maxiv.ivmx;
    auto& isarc2d        = cmn_curvat.isarc2d;
    auto& sarccr3        = cmn_curvat.sarccr3;
    auto& ieasmcc2d      = cmn_curvat.ieasmcc2d;
    auto& isstrc         = cmn_curvat.isstrc;
    auto& sstrc_crc      = cmn_curvat.sstrc_crc;
    auto& isar           = cmn_curvat.isar;
    auto& crot           = cmn_curvat.crot;
    auto& isarc3d        = cmn_curvat.isarc3d;
    auto& irigb          = cmn_rigidbody.irigb;
    auto& irbtrim        = cmn_rigidbody.irbtrim;
    auto& tmass          = cmn_rbstmt2.tmass;
    auto& yinert         = cmn_rbstmt2.yinert;
    auto& uinfrb         = cmn_rbstmt2.uinfrb;
    auto& qinfrb         = cmn_rbstmt2.qinfrb;
    auto& greflrb        = cmn_rbstmt2.greflrb;
    auto& gaccel         = cmn_rbstmt2.gaccel;
    auto& crefrb         = cmn_rbstmt2.crefrb;
    auto& xtmref         = cmn_rbstmt2.xtmref;
    auto& areat          = cmn_rbstmt2.areat;
    // filenam2
    char* avgq2          = cmn_filenam2.avgq2;
    char* avgq2pert      = cmn_filenam2.avgq2pert;
    char* clcds          = cmn_filenam2.clcds;
    char* clcdp          = cmn_filenam2.clcdp;
    char* output_dir     = cmn_filenam2.output_dir;
    // trim
    auto& dmtrmn         = cmn_trim.dmtrmn;
    auto& dmtrmnm        = cmn_trim.dmtrmnm;
    auto& dlcln          = cmn_trim.dlcln;
    auto& dlclnm         = cmn_trim.dlclnm;
    auto& trtol          = cmn_trim.trtol;
    auto& cmy            = cmn_trim.cmy;
    auto& cnw            = cmn_trim.cnw;
    auto& alf0           = cmn_trim.alf0;
    auto& alf1           = cmn_trim.alf1;
    auto& dzdt           = cmn_trim.dzdt;
    auto& thtd0          = cmn_trim.thtd0;
    auto& thtd1          = cmn_trim.thtd1;
    auto& zrg0           = cmn_trim.zrg0;
    auto& zrg1           = cmn_trim.zrg1;
    auto& dtrmsmx        = cmn_trim.dtrmsmx;
    auto& dtrmsmn        = cmn_trim.dtrmsmn;
    auto& dalfmx         = cmn_trim.dalfmx;
    auto& ddtmx          = cmn_trim.ddtmx;
    auto& ddtrm0         = cmn_trim.ddtrm0;
    auto& ddtrm1         = cmn_trim.ddtrm1;
    auto& itrmt          = cmn_trim.itrmt;
    auto& itrminc        = cmn_trim.itrminc;
    // tp(4,4) is cmn_trim.tp[16], column-major: tp(i,j) -> tp[(j-1)*4+(i-1)]
    auto& zlfct          = cmn_trim.zlfct;
    auto& epstr          = cmn_trim.epstr;
    auto& relax          = cmn_trim.relax;
    auto& ittrst         = cmn_trim.ittrst;
    // lesinfo
    auto& les_model      = cmn_lesinfo.les_model;
    auto& les_wallscale  = cmn_lesinfo.les_wallscale;
    auto& cs_smagorinsky = cmn_lesinfo.cs_smagorinsky;
    auto& cs_wale        = cmn_lesinfo.cs_wale;
    auto& cs_vreman      = cmn_lesinfo.cs_vreman;
    auto& xdir_only_source = cmn_sourceterm.xdir_only_source;
    auto& randomize      = cmn_random_input.randomize;
    auto& iexact_trunc   = cmn_mms.iexact_trunc;
    auto& iexact_disc    = cmn_mms.iexact_disc;
    auto& iexact_ring    = cmn_mms.iexact_ring;
    auto& isminc         = cmn_sminn.isminc;
    auto& ismincforce    = cmn_sminn.ismincforce;
    auto& scal_ic        = cmn_initfac.scal_ic;
    auto& ifunct         = cmn_plot3dtyp.ifunct;
    auto& lowmem_ux      = cmn_memry.lowmem_ux;
    auto& i_nonlin       = cmn_constit.i_nonlin;
    auto& c_nonlin       = cmn_constit.c_nonlin;
    auto& snonlin_lim    = cmn_constit.snonlin_lim;
    auto& i_tauijs       = cmn_constit.i_tauijs;
    auto& i_qcr2000      = cmn_constit.i_qcr2000;
    auto& i_qcr2013      = cmn_constit.i_qcr2013;
    auto& i_qcr2013v     = cmn_constit.i_qcr2013v;
    auto& cmulim         = cmn_easmlim.cmulim;
    auto& c10            = cmn_easmv.c10;
    auto& c11            = cmn_easmv.c11;
    auto& c2             = cmn_easmv.c2;
    auto& c3             = cmn_easmv.c3;
    auto& c4             = cmn_easmv.c4;
    auto& c5             = cmn_easmv.c5;
    auto& sigk1          = cmn_easmv.sigk1;
    auto& cmuc1          = cmn_easmv.cmuc1;
    auto& ieasm_type     = cmn_easmv.ieasm_type;
    auto& ipatch1st      = cmn_is_patch.ipatch1st;
    auto& issglrrw2012   = cmn_reystressmodel.issglrrw2012;
    auto& i_sas_rsm      = cmn_reystressmodel.i_sas_rsm;
    auto& i_yapterm      = cmn_reystressmodel.i_yapterm;
    auto& ifort50write   = cmn_writestuff.ifort50write;
    auto& j_ifort50write = cmn_writestuff.j_ifort50write;
    auto& i_ifort50write = cmn_writestuff.i_ifort50write;
    auto& i_saneg        = cmn_sa_options.i_saneg;
    auto& i_sanoft2      = cmn_sa_options.i_sanoft2;
    auto& sa_cw2         = cmn_sa_options.sa_cw2;
    auto& sa_cw3         = cmn_sa_options.sa_cw3;
    auto& sa_cv1         = cmn_sa_options.sa_cv1;
    auto& sa_ct3         = cmn_sa_options.sa_ct3;
    auto& sa_ct4         = cmn_sa_options.sa_ct4;
    auto& sa_cb1         = cmn_sa_options.sa_cb1;
    auto& sa_cb2         = cmn_sa_options.sa_cb2;
    auto& sa_sigma       = cmn_sa_options.sa_sigma;
    auto& sa_karman      = cmn_sa_options.sa_karman;
    auto& ilamlo         = cmn_lam.ilamlo;
    auto& ilamhi         = cmn_lam.ilamhi;
    auto& jlamlo         = cmn_lam.jlamlo;
    auto& jlamhi         = cmn_lam.jlamhi;
    auto& klamlo         = cmn_lam.klamlo;
    auto& klamhi         = cmn_lam.klamhi;
    auto& i_lam_forcezero= cmn_lam.i_lam_forcezero;
    auto& i_specialtop_kmax1001   = cmn_specialtop_kmax1001.i_specialtop_kmax1001;
    auto& a_specialtop_kmax1001   = cmn_specialtop_kmax1001.a_specialtop_kmax1001;
    auto& xc_specialtop_kmax1001  = cmn_specialtop_kmax1001.xc_specialtop_kmax1001;
    auto& sig_specialtop_kmax1001 = cmn_specialtop_kmax1001.sig_specialtop_kmax1001;
    auto& vtp_specialtop_kmax1001 = cmn_specialtop_kmax1001.vtp_specialtop_kmax1001;
    auto& wc_specialtop_kmax1001  = cmn_specialtop_kmax1001.wc_specialtop_kmax1001;
    auto& fac_specialtop_kmax1001 = cmn_specialtop_kmax1001.fac_specialtop_kmax1001;
    auto& cc_specialtop_kmax1001  = cmn_specialtop_kmax1001.cc_specialtop_kmax1001;
    auto& xerf_specialtop_kmax1001  = cmn_specialtop_kmax1001.xerf_specialtop_kmax1001;
    auto& sigerf_specialtop_kmax1001= cmn_specialtop_kmax1001.sigerf_specialtop_kmax1001;
    auto& iupdatemean    = cmn_iupdate.iupdatemean;
    auto& uub_2034       = cmn_wallvel_2034.uub_2034;
    auto& vvb_2034       = cmn_wallvel_2034.vvb_2034;
    auto& wwb_2034       = cmn_wallvel_2034.wwb_2034;
    auto& vvb_xaxisrot2034 = cmn_wallvel_2034.vvb_xaxisrot2034;



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
    cmn_turbconv.cflturb[0]  = 0.f;
    cmn_turbconv.cflturb[1]  = 0.f;
    cmn_turbconv.cflturb[2]  = 0.f;
    cmn_turbconv.cflturb[3]  = 0.f;
    cmn_turbconv.cflturb[4]  = 0.f;
    cmn_turbconv.cflturb[5]  = 0.f;
    cmn_turbconv.cflturb[6]  = 0.f;
    nkey     = nkey + 7;
    // entropy fix parameters
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
    // flag for writing unformatted/formatted plot3d files
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
    // additional user-specified memory to allocate
    memadd    = 0;
    nkey      = nkey + 1;
    memaddi   = 0;
    nkey      = nkey + 1;
    // flag to turn off stops when negative volumes/bad metrics are encountered
    negvol    = 0;
    nkey      = nkey + 1;
    // flag to turn on/off time marching flow solver
    meshdef   = 0;
    nkey      = nkey + 1;
    // eddy viscosity limiter for two eqn. turbulence models
    edvislim  = 1.e10f;
    nkey      = nkey + 1;
    // flag to set whether approximate or full production term is used in EASM models
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
    ides = 0;
    nkey     = nkey + 1;
    // constant associated with DES or DDES
    cdes = 0.65f;
    nkey     = nkey + 1;
    // additional constant associated with Modified DDES
    cddes = 0.975f;
    nkey     = nkey + 1;
    // flag to store iteration-averaged q values
    iteravg = 0;
    nkey     = nkey + 1;


    // ititr = 0
    ititr = 0;
    nkey  = 0;
    // ratio of specific heats
    gamma  = 1.4f; nkey = nkey + 1;
    // Prandtl number
    pr     = 0.72f; nkey = nkey + 1;
    // Turbulent Prandtl number
    prt    = 0.9f; nkey = nkey + 1;
    // reference temperature for Sutherland's Law
    cbar   = 198.6f; nkey = nkey + 1;
    // tolerance for collapsed metrics
    atol   = std::max(1.e-07f, (float)std::pow(10.0, (double)(-iexp+1))); nkey = nkey + 1;
    // complex perturbations
    xmach_img = 0.f; nkey = nkey + 1;
    alpha_img = 0.f; nkey = nkey + 1;
    beta_img  = 0.f; nkey = nkey + 1;
    reue_img  = 0.f; nkey = nkey + 1;
    tinf_img  = 0.f; nkey = nkey + 1;
    geom_img  = 0.f; nkey = nkey + 1;
    // CGNS file flag
    icgns     = 0; nkey = nkey + 1;
    // preconditioning
    cprec     = 0.f; nkey = nkey + 1;
    uref      = xmach; nkey = nkey + 1;
    avn       = 1.0f; nkey = nkey + 1;
    // specified Cl
    cltarg    = 99999.0f; ialphit = 0; nkey = nkey + 1;
    // angle-of-attack relaxation factor
    rlxalph   = 1.f; nkey = nkey + 1;
    // no. of subiterations in turbulence model
    nsubturb  = 1; nkey = nkey + 1;
    // cfl number for turbulence model
    cmn_turbconv.cflturb[0] = 0.f;
    cmn_turbconv.cflturb[1] = 0.f;
    cmn_turbconv.cflturb[2] = 0.f;
    cmn_turbconv.cflturb[3] = 0.f;
    cmn_turbconv.cflturb[4] = 0.f;
    cmn_turbconv.cflturb[5] = 0.f;
    cmn_turbconv.cflturb[6] = 0.f;
    nkey = nkey + 7;
    // entropy fix
    epsa_r = 0.f; epsa_l = 2.f*epsa_r; nkey = nkey + 1;
    // freeze turbulence model
    nfreeze  = 0; nkey = nkey + 1;
    // exact volume terms on 1-1 block boundaries
    ivolint  = 1; nkey = nkey + 1;
    // mesh deformation in steady-state mode
    idef_ss  = 0; nkey = nkey + 1;
    // unformatted/formatted plot3d files
    ibin     = 1; nkey = nkey + 1;
    // iblank array
    iblnk    = 1; nkey = nkey + 1;
    // blanking of fringe points
    iblnkfr  = 1; nkey = nkey + 1;
    // solution output to derivative output
    ip3dgrad = 0; nkey = nkey + 1;
    // complex perturbation to surface grid
    surf_img  = 0.f; nkey = nkey + 1;
    // additional user-specified memory
    memadd    = 0; nkey = nkey + 1;
    memaddi   = 0; nkey = nkey + 1;
    // negative volumes
    negvol    = 0; nkey = nkey + 1;
    // time marching flow solver
    meshdef   = 0; nkey = nkey + 1;
    // eddy viscosity limiter
    edvislim  = 1.e10f; nkey = nkey + 1;
    // approximate or full production term in EASM models
    iturbprod = 0; nkey = nkey + 1;
    // ghost-cell data from restart file
    irghost   = 1; nkey = nkey + 1;
    // ghost-cell data to restart file
    iwghost   = 1; nkey = nkey + 1;
    // limit of alpha change
    dalim     = 0.2f; nkey = nkey + 1;
    // cycles between alpha updates
    icycupdt  = 1; nkey = nkey + 1;
    // non-inertial reference frame flag
    noninflag = 0; nkey = nkey + 1;
    // rotation center
    xcentrot = 0.f; nkey = nkey + 1;
    ycentrot = 0.f; nkey = nkey + 1;
    zcentrot = 0.f; nkey = nkey + 1;
    // rotation rate
    xrotrate = 0.f; nkey = nkey + 1;
    yrotrate = 0.f; nkey = nkey + 1;
    zrotrate = 0.f; nkey = nkey + 1;
    // complex perturbation to rotation rate
    xrotrate_img = 0.f; nkey = nkey + 1;
    yrotrate_img = 0.f; nkey = nkey + 1;
    zrotrate_img = 0.f; nkey = nkey + 1;
    // 2nd order-time turbulence terms
    itime2read = 1; nkey = nkey + 1;
    // time-accuracy of turbulence model
    itaturb = 1; nkey = nkey + 1;
    // DES with turbulence model
    ides = 0; nkey = nkey + 1;
    cdes = 0.65f; nkey = nkey + 1;
    cddes = 0.975f; nkey = nkey + 1;
    // iteration-averaged q values
    iteravg = 0; nkey = nkey + 1;


    // turbulent quantity freestream levels
    cmn_ivals.tur10[0] = -1.f;
    cmn_ivals.tur10[1] = -1.f;
    cmn_ivals.tur10[2] = -1.f;
    cmn_ivals.tur10[3] = -1.f;
    cmn_ivals.tur10[4] = -1.f;
    cmn_ivals.tur10[5] = -1.f;
    cmn_ivals.tur10[6] = -1.f;
    nkey = nkey + 7;
    // value that epsilon or omega is set to for 2-eqn models
    if (ivmx == 15) { tur1cut = -1.f; } else { tur1cut = 1.e-20f; }
    nkey = nkey + 1;
    // roll angle
    roll_angle = 0.f; nkey = nkey + 1;
    // k-omega/sst/k-epsilon parameter controlling production term
    ikoprod = 0; nkey = nkey + 1;
    // sst parameter controlling denominator of mut term
    isstdenom = 0; nkey = nkey + 1;
    // limiter term on Pk in the two-eqn models
    pklimterm = 20.f; nkey = nkey + 1;
    // parameter for use with particular axi cases
    iaxi2plane = 0; nkey = nkey + 1;
    // turbulence model advection terms order
    iturbord = 1; nkey = nkey + 1;
    // value that k is set to for 2-eqn models
    tur2cut = 1.e-20f; nkey = nkey + 1;
    // lower limit on epsilon or omega for 2-eqn models
    tur1cutlev = 0.f; nkey = nkey + 1;
    // lower limit on k for 2-eqn models
    tur2cutlev = 0.f; nkey = nkey + 1;
    // full Navier-Stokes
    ifullns = 0; nkey = nkey + 1;
    // beta8 parameter in k-enstrophy model
    ibeta8kzeta = 0; nkey = nkey + 1;
    // SARC curvature correction in 2D
    isarc2d = 0; nkey = nkey + 1;
    // SARC curvature correction in 3D
    isarc3d = 0; nkey = nkey + 1;
    // parameter cr3 in SARC model
    sarccr3 = 1.0f; nkey = nkey + 1;
    // EASMCC curvature correction in 2D
    ieasmcc2d = 0; nkey = nkey + 1;
    // ipertavg
    ipertavg = 0; nkey = nkey + 1;
    // icoarsemovie
    icoarsemovie = 0; nkey = nkey + 1;
    // i2dmovie
    i2dmovie = 0; nkey = nkey + 1;
    // iclcd
    iclcd    = 0; nkey = nkey + 1;
    // iskip_blocks
    iskip_blocks = 1; nkey = nkey + 1;
    // cfltauMax
    cfltauMax = -1.0f; nkey = nkey + 1;
    // cfltau0
    cfltau0 = 1.0f; nkey = nkey + 1;
    // irbtrim
    irbtrim = 0; nkey = nkey + 1;
    // irigb
    irigb    = 0; nkey = nkey + 1;
    // greflrb
    greflrb  = 1.f; nkey = nkey + 1;
    // tmass
    tmass = 1.f; nkey = nkey + 1;
    // yinert
    yinert = 1.f; nkey = nkey + 1;
    // gaccel
    gaccel = 1.f; nkey = nkey + 1;
    // relax
    relax = 0.5f; nkey = nkey + 1;
    // itrminc
    itrminc = 5; nkey = nkey + 1;
    // dclda, dcldd, dcmda, dcmdd
    dclda    = 6.f; nkey = nkey + 1;
    dcldd    = 1.4f; nkey = nkey + 1;
    dcmda    = -0.2f; nkey = nkey + 1;
    dcmdd    = -0.88f; nkey = nkey + 1;
    // ndgrd, ndwrt
    ndgrd    = 0; nkey = nkey + 1;
    ndwrt    = 0; nkey = nkey + 1;
    // i_bsl
    i_bsl     = 0; nkey = nkey + 1;
    // keepambient
    keepambient = 0; nkey = nkey + 1;
    // re_thetat0
    re_thetat0 = -1.0f; nkey = nkey + 1;
    // turbintensity_inf_percent
    turbintensity_inf_percent = -1.0f; nkey = nkey + 1;
    // eddy_visc_inf
    eddy_visc_inf = -1.0f; nkey = nkey + 1;
    // cs_smagorinsky
    cs_smagorinsky = 0.0f; nkey = nkey + 1;
    // xdir_only_source
    xdir_only_source = 0.0f; nkey = nkey + 1;
    // randomize
    randomize = 0.0f; nkey = nkey + 1;
    // iexact_trunc, iexact_disc, iexact_ring
    iexact_trunc = 0; nkey = nkey + 1;
    iexact_disc  = 0; nkey = nkey + 1;
    iexact_ring  = 0; nkey = nkey + 1;
    // i_wilcox06, i_wilcox06_chiw
    i_wilcox06 = 0; nkey = nkey + 1;
    i_wilcox06_chiw = 1; nkey = nkey + 1;
    // i_turbprod_kterm
    i_turbprod_kterm = 0; nkey = nkey + 1;
    // i_catris_kw
    i_catris_kw = 0; nkey = nkey + 1;
    // ismincforce
    ismincforce = -1; nkey = nkey + 1;
    // prod2d3dtrace
    prod2d3dtrace = 0.f; nkey = nkey + 1;
    // i_compress_correct
    i_compress_correct = 0; nkey = nkey + 1;
    // les_model, les_wallscale
    les_model = 0; nkey = nkey + 1;
    les_wallscale = 0; nkey = nkey + 1;
    // cs_wale, cs_vreman
    cs_wale = 0.0f; nkey = nkey + 1;
    cs_vreman = 0.0f; nkey = nkey + 1;
    // isstrc, sstrc_crc
    isstrc = 0; nkey = nkey + 1;
    sstrc_crc = 1.4f; nkey = nkey + 1;
    // isstsf
    isstsf = 0; nkey = nkey + 1;
    // scal_ic
    scal_ic = 5.e6f; nkey = nkey + 1;
    // ifunct
    ifunct = 0; nkey = nkey + 1;
    // lowmem_ux
    lowmem_ux = 0; nkey = nkey + 1;
    // isar, crot
    isar = 0; nkey = nkey + 1;
    crot = 2.0f; nkey = nkey + 1;
    // i_nonlin, c_nonlin, snonlin_lim
    i_nonlin = 0; nkey = nkey + 1;
    c_nonlin = 0.3f; nkey = nkey + 1;
    snonlin_lim = 1.e-10f; nkey = nkey + 1;
    // isubit_r
    isubit_r = 0; nkey = nkey + 1;
    // i_wilcox98, i_wilcox98_chiw
    i_wilcox98 = 0; nkey = nkey + 1;
    i_wilcox98_chiw = 1; nkey = nkey + 1;
    // cmulim
    cmulim = .0005f; nkey = nkey + 1;
    // iaxi2planeturb, istrongturbdis
    iaxi2planeturb = 0; nkey = nkey + 1;
    istrongturbdis = 0; nkey = nkey + 1;
    // ieasm_type
    ieasm_type = 0; nkey = nkey + 1;
    // ipatch1st
    ipatch1st = 0; nkey = nkey + 1;
    // isst2003
    isst2003 = 0; nkey = nkey + 1;
    // issglrrw2012
    issglrrw2012 = 0; nkey = nkey + 1;
    // ifort50write, j_ifort50write, i_ifort50write
    ifort50write = 0; nkey = nkey + 1;
    j_ifort50write = 1; nkey = nkey + 1;
    i_ifort50write = 1; nkey = nkey + 1;
    // i_sas_rsm
    i_sas_rsm = 0; nkey = nkey + 1;
    // iforcev0
    iforcev0 = 0; nkey = nkey + 1;
    // i_saneg, i_sanoft2
    i_saneg = 0; nkey = nkey + 1;
    i_sanoft2 = 0; nkey = nkey + 1;
    // i_lam_forcezero
    i_lam_forcezero = 0; nkey = nkey + 1;
    // i_specialtop_kmax1001
    i_specialtop_kmax1001 = 0;
    a_specialtop_kmax1001 = 0.f;
    xc_specialtop_kmax1001 = 0.f;
    sig_specialtop_kmax1001 = 1.f;
    vtp_specialtop_kmax1001 = 0.f;
    wc_specialtop_kmax1001 = 0.f;
    fac_specialtop_kmax1001 = 1.0f;
    cc_specialtop_kmax1001 = 0.5f;
    xerf_specialtop_kmax1001 = 100000.f;
    sigerf_specialtop_kmax1001 = 1.e-9f;
    nkey = nkey + 10;
    // SA constants
    sa_cw2 = 0.3f; sa_cw3 = 2.0f; sa_cv1 = 7.1f;
    sa_ct3 = 1.2f; sa_ct4 = 0.5f; sa_cb1 = 0.1355f;
    sa_cb2 = 0.622f; sa_sigma = 2.f/3.f; sa_karman = 0.41f;
    nkey = nkey + 9;
    // iupdatemean
    iupdatemean = 1; nkey = nkey + 1;
    // wall velocities for bc2034
    uub_2034 = 0.f; vvb_2034 = 0.f; wwb_2034 = 0.f; vvb_xaxisrot2034 = 0.f;
    nkey = nkey + 4;
    // i_tauijs, i_qcr2000, i_qcr2013, i_qcr2013v, i_yapterm
    i_tauijs = 0; nkey = nkey + 1;
    i_qcr2000 = 0; nkey = nkey + 1;
    i_qcr2013 = 0; nkey = nkey + 1;
    i_qcr2013v = 0; nkey = nkey + 1;
    i_yapterm = 0; nkey = nkey + 1;


    // check for keyword-driven inputs
    // read(iunit5,1593) inpstr  -- format a210: read full line
    {
        FILE* fp5 = fortran_get_unit(iunit5);
        char linebuf[256];
        memset(inpstr, ' ', 210);
        if (fgets(linebuf, sizeof(linebuf), fp5)) {
            int len = (int)strlen(linebuf);
            if (len > 0 && linebuf[len-1] == '\n') { linebuf[--len] = '\0'; }
            int copylen = (len < 210) ? len : 210;
            memcpy(inpstr, linebuf, copylen);
        }
    }
    if (inpstr[0] == '>') {
        goto label_1000;
    } else if (inpstr[0] == '<') {
        nou(1) = std::min(nou(1)+1, ibufdim);
        snprintf(bou(nou(1),1), 120,
            "   ERROR: keyword input must start with a '>'-line in the input file");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    } else {
        ititr = 1;
        // read(inpstr,10)(realval(i),i=1,20)  -- format 20a4
        // Copy raw 4-byte chunks from inpstr into title array as float bits
        for (i = 0; i < 20; i++) {
            float fval = 0.f;
            memcpy(&fval, inpstr + i*4, 4);
            cmn_info.title[i] = fval;
        }
        return;
    }

label_1000:
    if (iunit11 > 0) {
        fortran_write_unit(iunit11,
            ">--------------------- begin keyword-driven input section --------------------->");
    }
label_1100:
    {
        // read(iunit5,1593,end=1900) inpstr
        FILE* fp5 = fortran_get_unit(iunit5);
        char linebuf[256];
        memset(inpstr, ' ', 210);
        if (!fgets(linebuf, sizeof(linebuf), fp5)) {
            goto label_1900;
        }
        int len = (int)strlen(linebuf);
        if (len > 0 && linebuf[len-1] == '\n') { linebuf[--len] = '\0'; }
        int copylen = (len < 210) ? len : 210;
        memcpy(inpstr, linebuf, copylen);
    }
    if (inpstr[0] == '<') {
        if (iunit11 > 0) {
            fortran_write_unit(iunit11,
                "<---------------------- end keyword-driven input section ----------------------<");
        }
        goto label_close_section;
    }
    // not '<': check for comment or blank line
    if (inpstr[0] == '/' || inpstr[0] == ' ') goto label_1100;




    // echo the keyword character strings
    npos = 1; lc2 = 0;
    { int iflg_m1 = -1; parser_ns::parser(inpstr, npos, lc1, lc2, lcl, iflg_m1); }
    lclw = std::min(lcl, 80);
    if (iunit11 > 0) {
        char tmpbuf[211]; memset(tmpbuf, 0, sizeof(tmpbuf));
        int copylen = (lclw < 210) ? lclw : 210;
        memcpy(tmpbuf, inpstr, copylen);
        int tlen = copylen;
        while (tlen > 0 && tmpbuf[tlen-1] == ' ') tlen--;
        tmpbuf[tlen] = '\0';
        fortran_write_unit(iunit11, "%s", tmpbuf);
    }
    npos = 1; lc2 = 0; lcl = 210;
    { int iflg_1 = 1; parser_ns::parser(inpstr, npos, lc1, lc2, lcl, iflg_1); }



label_kw_dispatch:
    if (kw_eq("gamma")) { lc2++; read_float_kw(realval[0]); gamma = realval[0]; goto label_kw_done; }
    if (kw_eq("pr")) { lc2++; read_float_kw(realval[0]); pr = realval[0]; goto label_kw_done; }
    if (kw_eq("prt")) { lc2++; read_float_kw(realval[0]); prt = realval[0]; goto label_kw_done; }
    if (kw_eq("cbar")) { lc2++; read_float_kw(realval[0]); cbar = realval[0]; goto label_kw_done; }
    if (kw_eq("atol")) { lc2++; read_float_kw(realval[0]); atol = realval[0]; goto label_kw_done; }
    if (kw_eq("xmach_img")) { lc2++; read_float_kw(realval[0]); xmach_img = realval[0]; goto label_kw_done; }
    if (kw_eq("alpha_img")) { lc2++; read_float_kw(realval[0]); alpha_img = realval[0]; goto label_kw_done; }
    if (kw_eq("beta_img")) { lc2++; read_float_kw(realval[0]); beta_img = realval[0]; goto label_kw_done; }
    if (kw_eq("reue_img")) { lc2++; read_float_kw(realval[0]); reue_img = realval[0]; goto label_kw_done; }
    if (kw_eq("tinf_img")) { lc2++; read_float_kw(realval[0]); tinf_img = realval[0]; goto label_kw_done; }
    if (kw_eq("geom_img")) { lc2++; read_float_kw(realval[0]); geom_img = realval[0]; goto label_kw_done; }
    if (kw_eq("icgns")) {
        lc2++; read_int_kw(icgns);
        if (icgns == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            snprintf(bou(nou(1),1), 120, "stopping...must install with -cgnsdir=... in order to use _CGNS_ (icgns=1)");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        goto label_kw_done;
    }
    if (kw_eq("cprec")) { lc2++; read_float_kw(realval[0]); cprec = realval[0]; goto label_kw_done; }
    if (kw_eq("uref")) { lc2++; read_float_kw(realval[0]); uref = realval[0]; goto label_kw_done; }
    if (kw_eq("avn")) { lc2++; read_float_kw(realval[0]); avn = realval[0]; goto label_kw_done; }
    if (kw_eq("cltarg")) {
        lc2++; read_float_kw(realval[0]); cltarg = realval[0];
        if (cltarg != 99999.0f) ialphit = 1;
        goto label_kw_done;
    }
    if (kw_eq("rlxalph")) { lc2++; read_float_kw(realval[0]); rlxalph = realval[0]; goto label_kw_done; }
    if (kw_eq("nsubturb")) { lc2++; read_int_kw(nsubturb); goto label_kw_done; }
    if (kw_eq("cflturb")) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        snprintf(bou(nou(1),1), 120, "stopping...keyword cflturb no longer allowed... use cflturb1, 2, etc instead");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        goto label_kw_done;
    }
    if (kw_eq("cflturb1")) { lc2++; read_float_kw(realval[0]); cmn_turbconv.cflturb[0]=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("cflturb2")) { lc2++; read_float_kw(realval[0]); cmn_turbconv.cflturb[1]=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("cflturb3")) { lc2++; read_float_kw(realval[0]); cmn_turbconv.cflturb[2]=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("cflturb4")) { lc2++; read_float_kw(realval[0]); cmn_turbconv.cflturb[3]=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("cflturb5")) { lc2++; read_float_kw(realval[0]); cmn_turbconv.cflturb[4]=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("cflturb6")) { lc2++; read_float_kw(realval[0]); cmn_turbconv.cflturb[5]=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("cflturb7")) { lc2++; read_float_kw(realval[0]); cmn_turbconv.cflturb[6]=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("epsa_r")) { lc2++; read_float_kw(realval[0]); epsa_r=realval[0]; epsa_l=2.f*epsa_r; goto label_kw_done; }
    if (kw_eq("nfreeze")) { lc2++; read_int_kw(nfreeze); goto label_kw_done; }
    if (kw_eq("ivolint")) { lc2++; read_int_kw(ivolint); goto label_kw_done; }
    if (kw_eq("idef_ss")) { lc2++; read_int_kw(idef_ss); goto label_kw_done; }
    if (kw_eq("ibin")) { lc2++; read_int_kw(ibin); goto label_kw_done; }
    if (kw_eq("iblnk")) { lc2++; read_int_kw(iblnk); goto label_kw_done; }
    if (kw_eq("iblnkfr")) { lc2++; read_int_kw(iblnkfr); goto label_kw_done; }
    if (kw_eq("ip3dgrad")) { lc2++; read_int_kw(ip3dgrad); goto label_kw_done; }
    if (kw_eq("surf_img")) {
        lc2++; read_float_kw(realval[0]); surf_img = realval[0];
        if (idef_ss == 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            snprintf(bou(nou(1),1), 120, "stopping...must have idef_ss = 1 if surf_img .ne. 0.");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        goto label_kw_done;
    }
    if (kw_eq("memadd")) { lc2++; read_int_kw(memadd); goto label_kw_done; }
    if (kw_eq("memaddi")) { lc2++; read_int_kw(memaddi); goto label_kw_done; }
    if (kw_eq("negvol")) { lc2++; read_int_kw(negvol); goto label_kw_done; }
    if (kw_eq("meshdef")) { lc2++; read_int_kw(meshdef); goto label_kw_done; }
    if (kw_eq("edvislim")) { lc2++; read_float_kw(realval[0]); edvislim=realval[0]; goto label_kw_done; }
    if (kw_eq("iturbprod")) { lc2++; read_int_kw(iturbprod); goto label_kw_done; }
    if (kw_eq("irghost")) { lc2++; read_int_kw(irghost); goto label_kw_done; }
    if (kw_eq("iwghost")) { lc2++; read_int_kw(iwghost); goto label_kw_done; }
    if (kw_eq("dalim")) { lc2++; read_float_kw(realval[0]); dalim=realval[0]; goto label_kw_done; }
    if (kw_eq("icycupdt")) { lc2++; read_int_kw(icycupdt); goto label_kw_done; }
    if (kw_eq("noninflag")) { lc2++; read_int_kw(noninflag); goto label_kw_done; }
    if (kw_eq("xcentrot")) { lc2++; read_float_kw(realval[0]); xcentrot=realval[0]; goto label_kw_done; }
    if (kw_eq("ycentrot")) { lc2++; read_float_kw(realval[0]); ycentrot=realval[0]; goto label_kw_done; }
    if (kw_eq("zcentrot")) { lc2++; read_float_kw(realval[0]); zcentrot=realval[0]; goto label_kw_done; }
    if (kw_eq("xrotrate")) { lc2++; read_float_kw(realval[0]); xrotrate=realval[0]; goto label_kw_done; }
    if (kw_eq("yrotrate")) { lc2++; read_float_kw(realval[0]); yrotrate=realval[0]; goto label_kw_done; }
    if (kw_eq("zrotrate")) { lc2++; read_float_kw(realval[0]); zrotrate=realval[0]; goto label_kw_done; }
    if (kw_eq("xrotrate_img")) { lc2++; read_float_kw(realval[0]); xrotrate_img=realval[0]; goto label_kw_done; }
    if (kw_eq("yrotrate_img")) { lc2++; read_float_kw(realval[0]); yrotrate_img=realval[0]; goto label_kw_done; }
    if (kw_eq("zrotrate_img")) { lc2++; read_float_kw(realval[0]); zrotrate_img=realval[0]; goto label_kw_done; }
    if (kw_eq("itime2read")) { lc2++; read_int_kw(itime2read); goto label_kw_done; }
    if (kw_eq("itaturb")) { lc2++; read_int_kw(itaturb); goto label_kw_done; }
    if (kw_eq("ides")) { lc2++; read_int_kw(ides); goto label_kw_done; }
    if (kw_eq("cdes")) { lc2++; read_float_kw(realval[0]); cdes=realval[0]; goto label_kw_done; }
    if (kw_eq("cddes")) { lc2++; read_float_kw(realval[0]); cddes=realval[0]; goto label_kw_done; }
    if (kw_eq("iteravg")) { lc2++; read_int_kw(iteravg); goto label_kw_done; }
    if (kw_eq("tur10")) { lc2++; read_float_kw(realval[0]); cmn_ivals.tur10[0]=realval[0]; goto label_kw_done; }
    if (kw_eq("tur20")) { lc2++; read_float_kw(realval[0]); cmn_ivals.tur10[1]=realval[0]; goto label_kw_done; }
    if (kw_eq("tur30")) { lc2++; read_float_kw(realval[0]); cmn_ivals.tur10[2]=realval[0]; goto label_kw_done; }
    if (kw_eq("tur40")) { lc2++; read_float_kw(realval[0]); cmn_ivals.tur10[3]=realval[0]; goto label_kw_done; }
    if (kw_eq("tur50")) { lc2++; read_float_kw(realval[0]); cmn_ivals.tur10[4]=realval[0]; goto label_kw_done; }
    if (kw_eq("tur60")) { lc2++; read_float_kw(realval[0]); cmn_ivals.tur10[5]=realval[0]; goto label_kw_done; }
    if (kw_eq("tur70")) { lc2++; read_float_kw(realval[0]); cmn_ivals.tur10[6]=realval[0]; goto label_kw_done; }
    if (kw_eq("tur1cut")) { lc2++; read_float_kw(realval[0]); tur1cut=realval[0]; goto label_kw_done; }
    if (kw_eq("roll_angle")) { lc2++; read_float_kw(realval[0]); roll_angle=realval[0]; goto label_kw_done; }
    if (kw_eq("ikoprod")) { lc2++; read_int_kw(ikoprod); goto label_kw_done; }
    if (kw_eq("isstdenom")) { lc2++; read_int_kw(isstdenom); goto label_kw_done; }
    if (kw_eq("pklimterm")) { lc2++; read_float_kw(realval[0]); pklimterm=realval[0]; goto label_kw_done; }
    if (kw_eq("iaxi2plane")) { lc2++; read_int_kw(iaxi2plane); goto label_kw_done; }
    if (kw_eq("iturbord")) { lc2++; read_int_kw(iturbord); goto label_kw_done; }
    if (kw_eq("tur2cut")) { lc2++; read_float_kw(realval[0]); tur2cut=realval[0]; goto label_kw_done; }
    if (kw_eq("tur1cutlev")) { lc2++; read_float_kw(realval[0]); tur1cutlev=realval[0]; goto label_kw_done; }
    if (kw_eq("tur2cutlev")) { lc2++; read_float_kw(realval[0]); tur2cutlev=realval[0]; goto label_kw_done; }
    if (kw_eq("ifullns")) { lc2++; read_int_kw(ifullns); goto label_kw_done; }
    if (kw_eq("ibeta8kzeta")) { lc2++; read_int_kw(ibeta8kzeta); goto label_kw_done; }
    if (kw_eq("isarc2d")) { lc2++; read_int_kw(isarc2d); goto label_kw_done; }
    if (kw_eq("isarc3d")) { lc2++; read_int_kw(isarc3d); goto label_kw_done; }
    if (kw_eq("sarccr3")) { lc2++; read_float_kw(realval[0]); sarccr3=realval[0]; goto label_kw_done; }
    if (kw_eq("ieasmcc2d")) { lc2++; read_int_kw(ieasmcc2d); goto label_kw_done; }
    if (kw_eq("ipertavg")) { lc2++; read_int_kw(ipertavg); if (ipertavg != 0) iteravg=ipertavg; goto label_kw_done; }
    if (kw_eq("icoarsemovie")) { lc2++; read_int_kw(icoarsemovie); goto label_kw_done; }
    if (kw_eq("i2dmovie")) { lc2++; read_int_kw(i2dmovie); goto label_kw_done; }
    if (kw_eq("iclcd")) { lc2++; read_int_kw(iclcd); goto label_kw_done; }
    if (kw_eq("iskip_blocks")) { lc2++; read_int_kw(iskip_blocks); goto label_kw_done; }
    if (kw_eq("cfltauMax")) { lc2++; read_float_kw(realval[0]); cfltauMax=realval[0]; goto label_kw_done; }
    if (kw_eq("cfltau0")) { lc2++; read_float_kw(realval[0]); cfltau0=realval[0]; goto label_kw_done; }
    if (kw_eq("irbtrim")) { lc2++; read_int_kw(irbtrim); goto label_kw_done; }
    if (kw_eq("irigb")) { lc2++; read_int_kw(irigb); goto label_kw_done; }
    if (kw_eq("greflrb")) { lc2++; read_float_kw(realval[0]); greflrb=realval[0]; goto label_kw_done; }
    if (kw_eq("tmass")) { lc2++; read_float_kw(realval[0]); tmass=realval[0]; goto label_kw_done; }
    if (kw_eq("yinert")) { lc2++; read_float_kw(realval[0]); yinert=realval[0]; goto label_kw_done; }
    if (kw_eq("gaccel")) { lc2++; read_float_kw(realval[0]); gaccel=realval[0]; goto label_kw_done; }
    if (kw_eq("relax")) { lc2++; read_float_kw(realval[0]); relax=realval[0]; goto label_kw_done; }
    if (kw_eq("itrminc")) { lc2++; read_int_kw(itrminc); goto label_kw_done; }
    if (kw_eq("dclda")) { lc2++; read_float_kw(realval[0]); dclda=realval[0]; goto label_kw_done; }
    if (kw_eq("dcldd")) { lc2++; read_float_kw(realval[0]); dcldd=realval[0]; goto label_kw_done; }
    if (kw_eq("dcmda")) { lc2++; read_float_kw(realval[0]); dcmda=realval[0]; goto label_kw_done; }
    if (kw_eq("dcmdd")) { lc2++; read_float_kw(realval[0]); dcmdd=realval[0]; goto label_kw_done; }
    if (kw_eq("ndgrd")) { lc2++; read_int_kw(ndgrd); goto label_kw_done; }
    if (kw_eq("ndwrt")) { lc2++; read_int_kw(ndwrt); goto label_kw_done; }
    if (kw_eq("i_bsl")) { lc2++; read_int_kw(i_bsl); goto label_kw_done; }
    if (kw_eq("keepambient")) { lc2++; read_int_kw(keepambient); goto label_kw_done; }
    if (kw_eq("re_thetat0")) { lc2++; read_float_kw(realval[0]); re_thetat0=realval[0]; goto label_kw_done; }
    if (kw_eq("turbintensity_inf_percent")) { lc2++; read_float_kw(realval[0]); turbintensity_inf_percent=realval[0]; goto label_kw_done; }
    if (kw_eq("eddy_visc_inf")) { lc2++; read_float_kw(realval[0]); eddy_visc_inf=realval[0]; goto label_kw_done; }
    if (kw_eq("cs_smagorinsky")) { lc2++; read_float_kw(realval[0]); cs_smagorinsky=realval[0]; goto label_kw_done; }
    if (kw_eq("xdir_only_source")) { lc2++; read_float_kw(realval[0]); xdir_only_source=realval[0]; goto label_kw_done; }
    if (kw_eq("randomize")) { lc2++; read_float_kw(realval[0]); randomize=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("iexact_trunc")) { lc2++; read_int_kw(iexact_trunc); goto label_kw_done; }
    if (kw_eq("iexact_disc")) { lc2++; read_int_kw(iexact_disc); goto label_kw_done; }
    if (kw_eq("iexact_ring")) { lc2++; read_int_kw(iexact_ring); goto label_kw_done; }
    if (kw_eq("i_wilcox06")) { lc2++; read_int_kw(i_wilcox06); goto label_kw_done; }
    if (kw_eq("i_wilcox06_chiw")) { lc2++; read_int_kw(i_wilcox06_chiw); goto label_kw_done; }
    if (kw_eq("i_turbprod_kterm")) { lc2++; read_int_kw(i_turbprod_kterm); goto label_kw_done; }
    if (kw_eq("i_catris_kw")) { lc2++; read_int_kw(i_catris_kw); goto label_kw_done; }
    if (kw_eq("ismincforce")) { lc2++; read_int_kw(ismincforce); goto label_kw_done; }
    if (kw_eq("prod2d3dtrace")) { lc2++; read_float_kw(realval[0]); prod2d3dtrace=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("i_compress_correct")) { lc2++; read_int_kw(i_compress_correct); goto label_kw_done; }
    if (kw_eq("les_model")) { lc2++; read_int_kw(les_model); goto label_kw_done; }
    if (kw_eq("les_wallscale")) { lc2++; read_int_kw(les_wallscale); goto label_kw_done; }
    if (kw_eq("cs_wale")) { lc2++; read_float_kw(realval[0]); cs_wale=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("cs_vreman")) { lc2++; read_float_kw(realval[0]); cs_vreman=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("isstrc")) { lc2++; read_int_kw(isstrc); goto label_kw_done; }
    if (kw_eq("sstrc_crc")) { lc2++; read_float_kw(realval[0]); sstrc_crc=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("isstsf")) { lc2++; read_int_kw(isstsf); goto label_kw_done; }
    if (kw_eq("scal_ic")) { lc2++; read_float_kw(realval[0]); scal_ic=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("ifunct")) { lc2++; read_int_kw(ifunct); goto label_kw_done; }
    if (kw_eq("lowmem_ux")) { lc2++; read_int_kw(lowmem_ux); goto label_kw_done; }
    if (kw_eq("isar")) { lc2++; read_int_kw(isar); goto label_kw_done; }
    if (kw_eq("crot")) { lc2++; read_float_kw(realval[0]); crot=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("i_nonlin")) { lc2++; read_int_kw(i_nonlin); goto label_kw_done; }
    if (kw_eq("c_nonlin")) { lc2++; read_float_kw(realval[0]); c_nonlin=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("snonlin_lim")) { lc2++; read_float_kw(realval[0]); snonlin_lim=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("isubit_r")) { lc2++; read_int_kw(isubit_r); goto label_kw_done; }
    if (kw_eq("i_wilcox98")) { lc2++; read_int_kw(i_wilcox98); goto label_kw_done; }
    if (kw_eq("i_wilcox98_chiw")) { lc2++; read_int_kw(i_wilcox98_chiw); goto label_kw_done; }
    if (kw_eq("cmulim")) { lc2++; read_float_kw(realval[0]); cmulim=std::abs(realval[0]); goto label_kw_done; }
    if (kw_eq("iaxi2planeturb")) { lc2++; read_int_kw(iaxi2planeturb); goto label_kw_done; }
    if (kw_eq("istrongturbdis")) { lc2++; read_int_kw(istrongturbdis); goto label_kw_done; }
    if (kw_eq("ieasm_type")) { lc2++; read_int_kw(ieasm_type); goto label_kw_done; }
    if (kw_eq("ipatch1st")) { lc2++; read_int_kw(ipatch1st); goto label_kw_done; }
    if (kw_eq("isst2003")) { lc2++; read_int_kw(isst2003); goto label_kw_done; }
    if (kw_eq("issglrrw2012")) { lc2++; read_int_kw(issglrrw2012); goto label_kw_done; }
    if (kw_eq("ifort50write")) { lc2++; read_int_kw(ifort50write); goto label_kw_done; }
    if (kw_eq("j_ifort50write")) { lc2++; read_int_kw(j_ifort50write); goto label_kw_done; }
    if (kw_eq("i_ifort50write")) { lc2++; read_int_kw(i_ifort50write); goto label_kw_done; }
    if (kw_eq("i_sas_rsm")) { lc2++; read_int_kw(i_sas_rsm); goto label_kw_done; }
    if (kw_eq("iforcev0")) { lc2++; read_int_kw(iforcev0); goto label_kw_done; }
    if (kw_eq("i_saneg")) { lc2++; read_int_kw(i_saneg); goto label_kw_done; }
    if (kw_eq("i_sanoft2")) { lc2++; read_int_kw(i_sanoft2); goto label_kw_done; }
    if (kw_eq("i_lam_forcezero")) { lc2++; read_int_kw(i_lam_forcezero); goto label_kw_done; }
    if (kw_eq("i_specialtop_kmax1001")) { lc2++; read_int_kw(i_specialtop_kmax1001); goto label_kw_done; }
    if (kw_eq("a_specialtop_kmax1001")) { lc2++; read_float_kw(realval[0]); a_specialtop_kmax1001=realval[0]; goto label_kw_done; }
    if (kw_eq("xc_specialtop_kmax1001")) { lc2++; read_float_kw(realval[0]); xc_specialtop_kmax1001=realval[0]; goto label_kw_done; }
    if (kw_eq("sig_specialtop_kmax1001")) { lc2++; read_float_kw(realval[0]); sig_specialtop_kmax1001=realval[0]; goto label_kw_done; }
    if (kw_eq("vtp_specialtop_kmax1001")) { lc2++; read_float_kw(realval[0]); vtp_specialtop_kmax1001=realval[0]; goto label_kw_done; }
    if (kw_eq("wc_specialtop_kmax1001")) { lc2++; read_float_kw(realval[0]); wc_specialtop_kmax1001=realval[0]; goto label_kw_done; }
    if (kw_eq("fac_specialtop_kmax1001")) { lc2++; read_float_kw(realval[0]); fac_specialtop_kmax1001=realval[0]; goto label_kw_done; }
    if (kw_eq("cc_specialtop_kmax1001")) { lc2++; read_float_kw(realval[0]); cc_specialtop_kmax1001=realval[0]; goto label_kw_done; }
    if (kw_eq("xerf_specialtop_kmax1001")) { lc2++; read_float_kw(realval[0]); xerf_specialtop_kmax1001=realval[0]; goto label_kw_done; }
    if (kw_eq("sigerf_specialtop_kmax1001")) { lc2++; read_float_kw(realval[0]); sigerf_specialtop_kmax1001=realval[0]; goto label_kw_done; }
    if (kw_eq("sa_cw2")) { lc2++; read_float_kw(realval[0]); sa_cw2=realval[0]; goto label_kw_done; }
    if (kw_eq("sa_cw3")) { lc2++; read_float_kw(realval[0]); sa_cw3=realval[0]; goto label_kw_done; }
    if (kw_eq("sa_cv1")) { lc2++; read_float_kw(realval[0]); sa_cv1=realval[0]; goto label_kw_done; }
    if (kw_eq("sa_ct3")) { lc2++; read_float_kw(realval[0]); sa_ct3=realval[0]; goto label_kw_done; }
    if (kw_eq("sa_ct4")) { lc2++; read_float_kw(realval[0]); sa_ct4=realval[0]; goto label_kw_done; }
    if (kw_eq("sa_cb1")) { lc2++; read_float_kw(realval[0]); sa_cb1=realval[0]; goto label_kw_done; }
    if (kw_eq("sa_cb2")) { lc2++; read_float_kw(realval[0]); sa_cb2=realval[0]; goto label_kw_done; }
    if (kw_eq("sa_sigma")) { lc2++; read_float_kw(realval[0]); sa_sigma=realval[0]; goto label_kw_done; }
    if (kw_eq("sa_karman")) { lc2++; read_float_kw(realval[0]); sa_karman=realval[0]; goto label_kw_done; }
    if (kw_eq("iupdatemean")) { lc2++; read_int_kw(iupdatemean); goto label_kw_done; }
    if (kw_eq("uub_2034")) { lc2++; read_float_kw(realval[0]); uub_2034=realval[0]; goto label_kw_done; }
    if (kw_eq("vvb_2034")) { lc2++; read_float_kw(realval[0]); vvb_2034=realval[0]; goto label_kw_done; }
    if (kw_eq("wwb_2034")) { lc2++; read_float_kw(realval[0]); wwb_2034=realval[0]; goto label_kw_done; }
    if (kw_eq("vvb_xaxisrot2034")) { lc2++; read_float_kw(realval[0]); vvb_xaxisrot2034=realval[0]; goto label_kw_done; }
    if (kw_eq("i_tauijs")) { lc2++; read_int_kw(i_tauijs); goto label_kw_done; }
    if (kw_eq("i_qcr2000")) { lc2++; read_int_kw(i_qcr2000); goto label_kw_done; }
    if (kw_eq("i_qcr2013")) { lc2++; read_int_kw(i_qcr2013); goto label_kw_done; }
    if (kw_eq("i_qcr2013v")) { lc2++; read_int_kw(i_qcr2013v); goto label_kw_done; }
    if (kw_eq("i_yapterm")) { lc2++; read_int_kw(i_yapterm); goto label_kw_done; }
    // else: unknown keyword
    if (iunit11 > 0) {
        fortran_write_unit(iunit11, "*** STOPPING: The keyword above is not supported. ***");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    goto label_kw_done;

label_kw_done:
    goto label_1100;

label_close_section:
    // do any order-dependent things here before exiting
    // Set trim matrix tp from dclda, dcldd, dcmda, dcmdd
    // tp(i,j) column-major: flat index = (j-1)*4+(i-1)
    cmn_trim.tp[0] = dclda;  // tp(1,1)
    cmn_trim.tp[4] = dcldd;  // tp(1,2)
    cmn_trim.tp[1] = dcmda;  // tp(2,1)
    cmn_trim.tp[5] = dcmdd;  // tp(2,2)


    // Check for the existence of the averaging files
    if (iteravg == 2) {
        isthere = false;
        char filename[512];
        {
            int tlen = (int)strlen(output_dir);
            while (tlen > 0 && output_dir[tlen-1] == ' ') tlen--;
            if (ipertavg == 2) {
                snprintf(filename, sizeof(filename), "%.*s/cfl3d_avg_ruvwp.p3d", tlen, output_dir);
            } else {
                snprintf(filename, sizeof(filename), "%.*s/cfl3d_avgq.p3d", tlen, output_dir);
            }
        }
        {
            FILE* ftmp = fopen(filename, "r");
            isthere = (ftmp != nullptr);
            if (ftmp) fclose(ftmp);
        }
        if (isthere) {
            fortran_open_unit(150, filename, "rb");
            {
                FILE* f150 = fortran_get_unit(150);
                if (!f150 || fread(&idummy, sizeof(int), 1, f150) != 1) {
                    fortran_close_unit(150);
                    goto label_cs_1011;
                }
            }
            fortran_close_unit(150);
            goto label_cs_1013;
        } else {
            goto label_cs_1012;
        }
    } else {
        goto label_cs_1013;
    }

label_cs_1011:
    fortran_close_unit(150);
label_cs_1012:
    iteravg = 1;
    if (ipertavg == 2) ipertavg = 1;
label_cs_1013:
    ;

    // Check for the existence of the clcd.bin file
    if (iclcd == 2) {
        isthere = false;
        char filename2[512];
        {
            int tlen = (int)strlen(output_dir);
            while (tlen > 0 && output_dir[tlen-1] == ' ') tlen--;
            snprintf(filename2, sizeof(filename2), "%.*s/clcd.bin", tlen, output_dir);
        }
        {
            FILE* ftmp = fopen(filename2, "r");
            isthere = (ftmp != nullptr);
            if (ftmp) fclose(ftmp);
        }
        if (isthere) {
            fortran_open_unit(150, filename2, "rb");
            {
                FILE* f150 = fortran_get_unit(150);
                if (!f150 || fread(&idummy1, sizeof(int), 1, f150) != 1 ||
                    fread(&idummy2, sizeof(int), 1, f150) != 1 ||
                    fread(&idummy3, sizeof(int), 1, f150) != 1) {
                    fortran_close_unit(150);
                    goto label_cs_1014;
                }
            }
            fortran_close_unit(150);
            goto label_cs_1016;
        } else {
            goto label_cs_1015;
        }
    } else {
        goto label_cs_1016;
    }

label_cs_1014:
    fortran_close_unit(150);
label_cs_1015:
    iclcd = 1;
label_cs_1016:
    ;

    // the variables turbintensity_inf_percent and eddy_visc_inf
    // do not go to the rest of the code; instead, if set, they
    // modify tur10(2) and tur10(1), respectively
    if (turbintensity_inf_percent >= 0.f) {
        if (ivmx == 72) {
            if (iunit11 > 0) fortran_write_unit(iunit11,
                " turbintensity_inf_percent always overrides tur10,20,30 for ivmx=72");
            zkinf = 1.5f * xmach*xmach * (turbintensity_inf_percent/100.f) * (turbintensity_inf_percent/100.f);
            cmn_ivals.tur10[0] = -2.f/3.f * zkinf;
            cmn_ivals.tur10[1] = -2.f/3.f * zkinf;
            cmn_ivals.tur10[2] = -2.f/3.f * zkinf;
        } else {
            if (iunit11 > 0) fortran_write_unit(iunit11,
                " turbintensity_inf_percent always overrides tur20");
            cmn_ivals.tur10[1] = 1.5f * xmach*xmach * (turbintensity_inf_percent/100.f) * (turbintensity_inf_percent/100.f);
        }
    }
    if (eddy_visc_inf >= 0.f) {
        if (ivmx == 72) {
            if (iunit11 > 0) fortran_write_unit(iunit11,
                " eddy_visc_inf always overrides tur70 for ivmx=72");
        } else {
            if (iunit11 > 0) fortran_write_unit(iunit11,
                " eddy_visc_inf always overrides tur10");
        }
        if (ivmx == 4) {
            cmn_ivals.tur10[0] = eddy_visc_inf / 0.09f;
        } else if (ivmx == 6 || ivmx == 7 || ivmx == 30 || ivmx == 40) {
            cmn_ivals.tur10[0] = cmn_ivals.tur10[1] / eddy_visc_inf;
        } else if (ivmx == 8 || ivmx == 14) {
            cmn_ivals.tur10[0] = 0.0895f * cmn_ivals.tur10[1] / eddy_visc_inf;
        } else if (ivmx == 12) {
            cmn_ivals.tur10[0] = 0.081f * cmn_ivals.tur10[1] / eddy_visc_inf;
        } else if (ivmx == 10 || ivmx == 15) {
            cmn_ivals.tur10[0] = 0.09f * cmn_ivals.tur10[1]*cmn_ivals.tur10[1] / eddy_visc_inf;
        } else if (ivmx == 9 || ivmx == 13) {
            cmn_ivals.tur10[0] = 0.0885f * cmn_ivals.tur10[1]*cmn_ivals.tur10[1] / eddy_visc_inf;
        } else if (ivmx == 11) {
            cmn_ivals.tur10[0] = 0.081f * cmn_ivals.tur10[1]*cmn_ivals.tur10[1] / eddy_visc_inf;
        } else if (ivmx == 5) {
            // for SA model, need to iterate to find nuwiggle
            cmn_ivals.tur10[0] = eddy_visc_inf;
            for (n = 1; n <= 100; n++) {
                float t1 = cmn_ivals.tur10[0];
                errornum = t1*t1*t1*t1 - eddy_visc_inf*t1*t1*t1 - 357.911f*eddy_visc_inf;
                errorden = 4.f*t1*t1*t1 - 3.f*eddy_visc_inf*t1*t1;
                error_val = errornum / errorden;
                if (std::abs(error_val) < 1.e-6f) goto label_105;
                cmn_ivals.tur10[0] = cmn_ivals.tur10[0] - error_val;
            }
            nou(1) = std::min(nou(1)+1, ibufdim);
            snprintf(bou(nou(1),1), 120,
                "stopping...Newton iteration using eddy_visc_inf not converged in readkey");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
label_105:
            ;
        } else if (ivmx == 72) {
            cmn_ivals.tur10[6] = -3.f/2.f * cmn_ivals.tur10[0] / eddy_visc_inf;
        }
    }
    return;



label_1900:
    // EOF encountered before '<' terminator
    nou(1) = std::min(nou(1)+1, ibufdim);
    snprintf(bou(nou(1),1), 120,
        "   ERROR: keyword input section must end with a '<'-line in the input file");
    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    return;
}

} // namespace readkey_ns
