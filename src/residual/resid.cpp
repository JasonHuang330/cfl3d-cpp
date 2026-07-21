// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "resid.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include "module_kwstm.h"
#include "barth3d.h"
#include "blomax.h"
#include "ctime1.h"
#include "delv.h"
#include "dthole.h"
#include "ffluxr.h"
#include "ffluxv.h"
#include "ffluxv1.h"
#include "foureqn.h"
#include "gfluxr.h"
#include "gfluxv.h"
#include "gfluxv1.h"
#include "hfluxr.h"
#include "hfluxv.h"
#include "hfluxv1.h"
#include "wmag.h"
#include "spalart.h"
#include "twoeqn.h"
#include "threeeqn.h"
#include "lesdiag.h"
#include "termn8.h"
#include "l2norm.h"
#include "hole.h"
#include "resnonin.h"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstring>

// Forward declarations aligned to the authoritative definitions in the linked
// objects (mms.cpp / ffluxv.cpp), matching namespace and argument ranks exactly.
namespace mms_ns {
    void exact_flow_force(int& jdim, int& kdim, int& idim,
                          FortranArray3DRef<double> x,
                          FortranArray3DRef<double> y,
                          FortranArray3DRef<double> z,
                          FortranArray3DRef<double> vol,
                          FortranArray4DRef<double> res,
                          FortranArray3DRef<double> smin,
                          int& iexact_trunc, int& iexact_disc);
    void exact_norm(int& jdim, int& kdim, int& idim,
                    FortranArray3DRef<double> x,
                    FortranArray3DRef<double> y,
                    FortranArray3DRef<double> z,
                    FortranArray4DRef<double> q,
                    FortranArray3DRef<double> tursav,
                    FortranArray3DRef<double> vol,
                    FortranArray4DRef<double> res,
                    FortranArray3DRef<double> smin,
                    int& iexact_trunc, int& iexact_disc);
    void zero_resid_ring(int& jdim, int& kdim, int& idim,
                         FortranArray4DRef<double> res,
                         int& jj, int& kk, int& ii,
                         int& numeq, int& numrows,
                         int& iexact_trunc, int& iexact_disc);
}
// ffluxv.h wrongly declares vmui as a 3D ref; the real definition in ffluxv.cpp
// takes vmui as a 2D ref. Declare the authoritative 2D-ref overload here so the
// call below resolves to (and links against) the actual definition.
namespace ffluxv_ns {
    void ffluxv(int& k, int& npl, int& jdim, int& kdim, int& idim,
                FortranArray4DRef<double> res, FortranArray4DRef<double> q,
                FortranArray4DRef<double> qi0, FortranArray4DRef<double> si,
                FortranArray3DRef<double> vol, FortranArray2DRef<double> t,
                int& nvtq, FortranArray2DRef<double> wi0,
                FortranArray3DRef<double> vist3d, FortranArray2DRef<double> vmui,
                FortranArray4DRef<double> vi0, FortranArray3DRef<double> bci,
                FortranArray4DRef<double> zksav, FortranArray4DRef<double> ti0,
                FortranArray3DRef<double> cmuv, FortranArray3DRef<double> voli0,
                FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
                int& nbuf, int& ibufdim, int& iadv, int& nummem,
                FortranArray4DRef<double> ux);
}

namespace resid_ns {

void resid(int& nbl, int& ntime, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q,
           FortranArray3DRef<double> qj0,
           FortranArray3DRef<double> qk0,
           FortranArray3DRef<double> qi0,
           FortranArray2DRef<double> sj,
           FortranArray2DRef<double> sk,
           FortranArray2DRef<double> si,
           FortranArray1DRef<double> vol,
           FortranArray1DRef<double> dtj,
           FortranArray1DRef<double> x,
           FortranArray1DRef<double> y,
           FortranArray1DRef<double> z,
           FortranArray3DRef<double> vist3d,
           FortranArray3DRef<double> snj0,
           FortranArray3DRef<double> snk0,
           FortranArray3DRef<double> sni0,
           FortranArray2DRef<double> res,
           FortranArray1DRef<double> wj0,
           FortranArray1DRef<double> wk0,
           FortranArray1DRef<double> wi0,
           FortranArray3DRef<double> vmuk,
           FortranArray3DRef<double> vmuj,
           FortranArray3DRef<double> vmui,
           FortranArray1DRef<double> wk,
           int& nwork,
           int& isf,
           FortranArray1DRef<int> iwfa,
           FortranArray1DRef<double> wfa,
           double& delt,
           FortranArray3DRef<double> blank,
           int& iover,
           int& nblendg,
           int& nblstat,
           int& nblstag,
           FortranArray4DRef<double> xib,
           FortranArray2DRef<double> sig,
           FortranArray2DRef<double> sqtq,
           FortranArray2DRef<double> g,
           FortranArray4DRef<double> tj0,
           FortranArray4DRef<double> tk0,
           FortranArray4DRef<double> ti0,
           FortranArray3DRef<double> xkb,
           FortranArray3DRef<double> blnum,
           FortranArray4DRef<double> vj0,
           FortranArray4DRef<double> vk0,
           FortranArray4DRef<double> vi0,
           FortranArray3DRef<double> bcj,
           FortranArray3DRef<double> bck,
           FortranArray3DRef<double> bci,
           int& nt,
           FortranArray1DRef<double> sumn,
           FortranArray1DRef<int> negn,
           FortranArray4DRef<double> ux,
           FortranArray4DRef<double> xib2,
           FortranArray3DRef<double> cmuv,
           FortranArray3DRef<double> volj0,
           FortranArray3DRef<double> volk0,
           FortranArray3DRef<double> voli0,
           FortranArray1DRef<int> nou,
           FortranArray2DRef<char[120]> bou,
           int& nbuf,
           int& ibufdim,
           int& myid,
           FortranArray1DRef<int> mblk2nd,
           int& maxbl,
           int& maxseg,
           FortranArray1DRef<int> nbci0,
           FortranArray1DRef<int> nbcj0,
           FortranArray1DRef<int> nbck0,
           FortranArray1DRef<int> nbcidim,
           FortranArray1DRef<int> nbcjdim,
           FortranArray1DRef<int> nbckdim,
           FortranArray4DRef<int> ibcinfo,
           FortranArray4DRef<int> jbcinfo,
           FortranArray4DRef<int> kbcinfo,
           FortranArray1DRef<double> vormax,
           FortranArray1DRef<int> ivmax,
           FortranArray1DRef<int> jvmax,
           FortranArray1DRef<int> kvmax,
           FortranArray1DRef<int> idefrm,
           FortranArray1DRef<int> iadvance,
           FortranArray4DRef<double> qavg,
           int& nummem)
{
    // COMMON block aliases
    float& xnumavg     = cmn_avgdata.xnumavg;
    int32_t& iteravg   = cmn_avgdata.iteravg;
    int32_t& ipertavg  = cmn_avgdata.ipertavg;
    int32_t& ioutsub   = cmn_account.ioutsub;
    float& gm1         = cmn_fluid.gm1;
    float& dt0         = cmn_cfl.dt0;
    float* rkap0       = cmn_fvfds.rkap0;  // 0-based C array, Fortran 1-based
    float* rkap        = cmn_info.rkap;    // 0-based C array
    float& dt          = cmn_info.dt;
    float& fmax        = cmn_info.fmax;
    int32_t& nit       = cmn_info.nit;
    int32_t& nitfo     = cmn_info.nitfo;
    int32_t& iflagts   = cmn_info.iflagts;
    int32_t& nres      = cmn_info.nres;
    int32_t* ncyc1     = cmn_info.ncyc1;   // 0-based
    int32_t& levt      = cmn_mgrd.levt;
    int32_t& kode      = cmn_mgrd.kode;
    int32_t& ncyc      = cmn_mgrd.ncyc;
    int32_t& icyc      = cmn_mgrd.icyc;
    int32_t& level     = cmn_mgrd.level;
    int32_t& lglobal   = cmn_mgrd.lglobal;
    int32_t* nfajki    = cmn_nfablk.nfajki; // 0-based
    // ivisc_i is the int array from /reyue/
    int32_t* ivisc_i   = cmn_reyue.ivisc;  // 0-based
    int32_t& isklton   = cmn_sklton.isklton;
    int32_t& i2d       = cmn_twod.i2d;
    float& time        = cmn_unst.time;
    int32_t& ntstep    = cmn_unst.ntstep;
    int32_t& ita       = cmn_unst.ita;
    int32_t& iunst     = cmn_unst.iunst;
    int32_t& nwrest    = cmn_wrbl.nwrest;
    int32_t& iaccnt    = cmn_account.iaccnt;
    float& xcentrot    = cmn_noninertial.xcentrot;
    float& ycentrot    = cmn_noninertial.ycentrot;
    float& zcentrot    = cmn_noninertial.zcentrot;
    float& xrotrate    = cmn_noninertial.xrotrate;
    float& yrotrate    = cmn_noninertial.yrotrate;
    float& zrotrate    = cmn_noninertial.zrotrate;
    int32_t& noninflag = cmn_noninertial.noninflag;
    int32_t& ikoprod   = cmn_konew.ikoprod;
    int32_t& isstdenom = cmn_konew.isstdenom;
    int32_t& isst2003  = cmn_konew.isst2003;
    int32_t& isstsf    = cmn_konew.isstsf;
    int32_t& i_wilcox06= cmn_konew.i_wilcox06;
    int32_t& i_wilcox98= cmn_konew.i_wilcox98;
    int32_t& ivmx      = cmn_maxiv.ivmx;
    int32_t& iturbprod = cmn_turbconv.iturbprod;
    int32_t& nsubturb  = cmn_turbconv.nsubturb;
    int32_t& nfreeze   = cmn_turbconv.nfreeze;
    int32_t& iturbord  = cmn_turbconv.iturbord;
    int32_t& ifullns   = cmn_fullns.ifullns;
    int32_t& isarc2d   = cmn_curvat.isarc2d;
    int32_t& isarc3d   = cmn_curvat.isarc3d;
    int32_t& isar      = cmn_curvat.isar;
    int32_t& isstrc    = cmn_curvat.isstrc;
    int32_t& ieasmcc2d = cmn_curvat.ieasmcc2d;
    float& xdir_only_source = cmn_sourceterm.xdir_only_source;
    int32_t& iexact_trunc = cmn_mms.iexact_trunc;
    int32_t& iexact_disc  = cmn_mms.iexact_disc;
    int32_t& iexact_ring  = cmn_mms.iexact_ring;
    int32_t& ides      = cmn_des.ides;
    int32_t& lowmem_ux = cmn_memry.lowmem_ux;
    int32_t& i_nonlin  = cmn_constit.i_nonlin;
    int32_t& i_tauijs  = cmn_constit.i_tauijs;
    int32_t& issglrrw2012 = cmn_reystressmodel.issglrrw2012;
    int32_t& i_sas_rsm    = cmn_reystressmodel.i_sas_rsm;
    int32_t& i_yapterm    = cmn_reystressmodel.i_yapterm;

    // Local variables
    int jdim1, kdim1, idim1;
    int maxmem, nroom;
    double dttol, t2, fact;
    int iout, iterm;
    double dtmin;
    int imult;
    int nv, nvtq, n, nplq, npl, niwfac;
    int i, k, l, izz, ii, jj, jkv, kdww2;
    int nnpl, ist;
    int ipw, nblt, inmx;
    int iwk1, iwk2, iwk3, iwk4, iwk5, iwk6, iwk7, iwk8, iwk9, iwk10;
    int iwk11, iwk12, iwk13, iwk14, iwk15, iwk16, iwk17, iwk18, iwk19, iwk20;
    int iwk21, iwk22, iwk23, iwk24, iwk25, iwk26, iwk27, iwk28, iwk29, iwk30;
    int iwk31, iwk32, iwk33, iwk34, iwk35, iwk36, iwk37, iwk38, iwk39, iwk40;
    int iwk41, iwk42, iwk43, iwk44, iwk45, iwk46, iwk47, iwk48, iwk49, iwk50;
    int iwk51, iwk52, iwk53, iwk54, iwk55, iwk56, iwk57, iwk58, iwk59, iwk60;
    int iwk61, iwk62, iwk63, iwk64, iwk65, iwk66, iwk67, iwk68, iwk69, iwk70;
    int iwk71, iwk72, iwk73, iwk74, iwk75, iwk76, iwk77, iwk78, iwk79, iwk80;
    int iwk81, iwk82;
    int iex, iex2, iex3;
    int iccnum, isarcnum;
    int icall;
    double resd;
    int jbctyp[2], kbctyp[2], ibctyp[2];
    int nss;
    double fsrce;
    int myhost = 0; // myhost not passed, use 0

    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;

    if (jdim==17) { static int c=0; c++; int qb=0,vb=0; int qj=-1,qk=-1,qi=-1,qc=-1;
      for(int cc=1;cc<=5;cc++)for(int ii=1;ii<=idim;ii++)for(int kk=1;kk<=kdim;kk++)for(int jj=1;jj<=jdim;jj++){ if(std::isnan((double)q(jj,kk,ii,cc))){qb++; if(qj<0){qj=jj;qk=kk;qi=ii;qc=cc;}} }
      for(int ii=1;ii<=idim-1;ii++)for(int kk=1;kk<=kdim-1;kk++)for(int jj=1;jj<=jdim-1;jj++){ if(std::isnan((double)vist3d(jj,kk,ii)))vb++; }
    }

    //
    // check storage
    //
    maxmem = std::max(jdim*kdim*35, jdim*idim*35);
    if (nwork < maxmem) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " Stopping.  Insufficient memory in resid.  You must increase mwork.");
        {
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }

    if (icyc >= nitfo+1 && level >= lglobal) {
        //
        // finest global grids or embedded grids
        //
        rkap[0] = rkap0[0];
        rkap[1] = rkap0[1];
        rkap[2] = rkap0[2];
    } else if (level >= lglobal) {
        rkap[0] = -3.f;
        rkap[1] = -3.f;
        rkap[2] = -3.f;
    } else {
        //
        // coarser grids multigrid
        //
        rkap[0] = rkap0[0];
        rkap[1] = rkap0[1];
        rkap[2] = rkap0[2];
    }

    //
    // time step  kode.ne.-1
    //
    if (kode == -1 || ntime > 1) goto label130;

    if (isklton == 1 && level == levt) {
        if ((float)dt > 0.0f) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " ");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " initial time step, icyc=%12.3E%5d", (float)dt, icyc);
        }
        if ((float)dt < 0.0f) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " ");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " initial cfl number, icyc=%12.3E%5d", std::abs((float)dt), icyc);
        }
    }

    //
    // local time step: ramp initial CFL number to a value fmax times larger.
    // The time step distribution is recalculated for each cycle.
    //
    dttol = 1.e-6;
    if (iflagts > 0 && (float)dt < 0.f && icyc > 1) {
        if (nt == 1 && std::abs((float)dt) <
            (float)fmax * std::abs((float)dt0) - (float)dttol) {
            // ramp CFL on first block on top level
            if (level == levt && nbl == nblstat) {
                t2   = 1.0 / (double)iflagts;
                fact = std::pow((double)fmax, t2);
                dt   = (float)((double)dt * fact);
            }
        }
    }

    //
    // global time step: ramp initial time step to a value fmax times larger.
    // The CFL number distribution is recalculated for each time step.
    //
    dttol = 1.e-6;
    if (iflagts > 0 && (float)dt > 0.f && nt > 1) {
        if (icyc == 1 && std::abs((float)dt) <
            (float)fmax * std::abs((float)dt0) - (float)dttol) {
            // ramp time step on first block on top level
            if (level == levt && nbl == nblstat) {
                t2   = 1.0 / (double)iflagts;
                fact = std::pow((double)fmax, t2);
                dt   = (float)((double)dt * fact);
            }
        }
    }

    //
    // calculate time step/CFL based on CFL/timestep
    iout = 0;
    if (isklton == 1) iout = 1;
    iterm = std::max({ivisc_i[0], ivisc_i[1], ivisc_i[2]});

    //
    // don't call ctime1 during time-accurate subiterations if
    // only physical time terms are used...actually really only
    // need to call ctime1 once in that case
    //
    if ((float)dt > 0.f && ita > 0 && icyc > 1) goto label130;

    {
        // Build 3D refs for ctime1
        FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1); // ctime1 reads vol(izz,1,i): 3rd index i needs stride jdim*kdim (was (jkm,idim1,1) → i-stride jkm*idim1, OOB garbage for idim1>1; invisible in 2D where idim1==1)
        FortranArray3DRef<double> sj3d(&sj(1,1), jdim*kdim, idim1, 5);
        FortranArray3DRef<double> sk3d(&sk(1,1), jdim*kdim, idim1, 5);
        FortranArray3DRef<double> si3d(&si(1,1), jdim*kdim, idim, 5);
        FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim1); // ctime1.F: dtj(jdim,kdim,idim-1) — 2nd index is k
        // ctime1 expects: FortranArray2DRef<double> t (wk as 2D)
        FortranArray2DRef<double> wk_t(&wk(1), jdim*kdim, 17);
        int idef_nbl = idefrm(nbl);
        ctime1_ns::ctime1(nbl, jdim, kdim, idim, q, vol3d, sj3d, sk3d, si3d,
                          dtj3d, wk_t, delt, vist3d, iterm, dtmin, iout, ntime,
                          nou, bou, nbuf, ibufdim, idef_nbl);
    }

    //
    // zero out time step at hole points for chimera scheme
    //
    if (iover == 1) {
        FortranArray3DRef<double> dtj3d(&dtj(1), jdim*kdim, idim1, 1);
        FortranArray3DRef<double> vol3d(&vol(1), jdim*kdim, idim1, 1);
        dthole_ns::dthole(jdim, kdim, idim, dtj3d, vol3d, blank, dtmin,
                          nou, bou, nbuf, ibufdim);
    }

    if (iflagts > 0 && std::abs((float)fmax) > 0.f &&
        (float)dt < 0.f && myid <= 1) {
        if (((icyc == iflagts+1) || (icyc == ncyc && iflagts+1 > ncyc))
            && level == levt && nbl == nblstat && nt == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " ");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " final cfl number, icyc=%12.3E%5d", std::abs((float)dt), icyc);
            if (icyc != ncyc) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " ");
            }
        }
    }

    if (iflagts > 0 && std::abs((float)fmax) > 0.f &&
        (float)dt > 0.f && myid <= 1) {
        if (((nt == iflagts+1) || (nt == ntstep && iflagts+1 > ntstep))
            && level == levt && nbl == nblstat && icyc == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " ");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " final dt, time step=%12.3E%5d", (float)dt, nt);
            if (nt != ntstep) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " ");
            }
        }
    }

label130:;


    //
    // turbulent viscosity - finest grids only
    //
    if (iadvance(nbl) >= 0) {

    if (ivisc_i[2] > 1 || ivisc_i[1] > 1 || ivisc_i[0] > 1) {
        icall = 0;
        if (icyc == (icyc/1)*1 && level >= lglobal && ntime <= nit)
            icall = 1;
        if (lowmem_ux == 1) {
            if (ivmx == 8 || ivmx == 9 || ivmx >= 11 ||
                ivmx == 16 || ivmx == 72 ||
                ((ivmx == 6 || ivmx == 7 || ivmx == 10) && ikoprod == 1) ||
                (ivmx == 7 && isstdenom == 1) ||
                (ivmx == 7 && isst2003 == 1) ||
                ((ivmx == 6 || ivmx == 7) && isstrc > 0) ||
                ((ivmx == 6 || ivmx == 7) && isstsf == 1) ||
                ((ivmx == 5) && isarc2d == 1) ||
                ((ivmx == 5) && isarc3d == 1) ||
                ((ivmx == 5) && isar == 1) ||
                (ivmx == 6 && i_wilcox06 == 1) ||
                (ivmx == 6 && i_wilcox98 == 1) ||
                ides >= 2 || i_nonlin != 0 || i_tauijs != 0) {
                // Build 4D refs for delv
                FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
                FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                FortranArray3DRef<double> wt3d(&wk(1), jdim, kdim, idim1);
                FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                delv_ns::delv(jdim, kdim, idim, q, sj4d, sk4d, si4d, vol3d,
                              ux, wt3d, blank, iover,
                              qj04d, qk04d, qi04d, bcj, bck, bci, nbl,
                              volj0, volk0, voli0, maxbl, vormax, ivmax, jvmax, kvmax);
            }
        } else {
            FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
            FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
            FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
            FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
            FortranArray3DRef<double> wt3d(&wk(1), jdim, kdim, idim1);
            FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
            FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
            FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
            delv_ns::delv(jdim, kdim, idim, q, sj4d, sk4d, si4d, vol3d,
                          ux, wt3d, blank, iover,
                          qj04d, qk04d, qi04d, bcj, bck, bci, nbl,
                          volj0, volk0, voli0, maxbl, vormax, ivmax, jvmax, kvmax);
        }

        if (icall > 0) {
            //
            // vorticity magnitude
            //
            if (isklton == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " ");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " computing vorticity magnitude");
            }

            ipw = 0;
            if ((float)dt < 0.0f) {
                if (icyc/nwrest*nwrest == icyc || icyc == ncyc) {
                    ipw = 1;
                }
            } else {
                if (nt/nwrest*nwrest == nt || nt == ntstep) {
                    if (icyc == ioutsub) ipw = 1;
                }
            }
            {
                FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
                FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                // Fortran: call wmag(...,res,res(1,2),wk,...)
                //   vor = res              dim (jdim-1,kdim-1,idim-1)
                //   w   = res(1,2)         dim (jdim-1,kdim-1,idim-1,3)  (3 components)
                //   wt  = wk (scratch)     dim (jdim,kdim,3)
                FortranArray3DRef<double> vor3d(&res(1,1), jdim, kdim, idim1);
                FortranArray4DRef<double> w4d(&res(1,2), jdim, kdim, idim1, 3);
                FortranArray3DRef<double> wt3d(&wk(1), jdim, kdim, 3);
                FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                wmag_ns::wmag(jdim, kdim, idim, q, sj4d, sk4d, si4d, vol3d,
                              vor3d, w4d, wt3d, ipw, blank, iover,
                              qj04d, qk04d, qi04d, bcj, bck, bci, nbl,
                              volj0, volk0, voli0, vormax, ivmax, jvmax, kvmax, maxbl);
            }

            //
            // wall/wake turbulence model
            //
            if (isklton == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " ");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " evaluating turbulence model");
            }


            if (ivmx == 2) {
                // Baldwin-Lomax
                inmx = std::max({jdim, kdim, idim});
                iwk1 = 1;
                iwk2 = iwk1 + inmx;
                iwk3 = iwk2 + inmx;
                iwk4 = iwk3 + inmx;
                iwk5 = iwk4 + inmx;
                iwk6 = iwk5 + inmx;
                iwk7 = iwk6 + inmx;
                iwk8 = iwk7 + inmx;
                iwk9 = iwk8 + inmx;
                nroom = nwork - (iwk9 + inmx);
                if (nroom < 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " not enough memory for B-L turb model.");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " nroom=%12d", nroom);
                    termn8_ns::termn8(myid, const_cast<int&>(*(int*)"\xff\xff\xff\xff"), ibufdim, nbuf, bou, nou);
                    { int m1=-1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
                }
                nblt = nbl;
                ipw = 0;
                if ((float)dt < 0.f) {
                    if (icyc == ncyc) ipw = 1;
                    if (icyc == ncyc && icyc <= 3) ipw = 2;
                } else {
                    if (icyc == ioutsub) {
                        if (nt == ntstep) ipw = 1;
                        if (nt == ntstep && ntstep <= 3) ipw = 2;
                    }
                }
                {
                    FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                    FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                    FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                    FortranArray3DRef<double> vor3d(&res(1,1), jdim, kdim, idim1);
                    FortranArray3DRef<double> eoms3d(&res(1,2), jdim, kdim, idim1);
                    FortranArray3DRef<double> x3d(&x(1), jdim, kdim, idim);
                    FortranArray3DRef<double> y3d(&y(1), jdim, kdim, idim);
                    FortranArray3DRef<double> z3d(&z(1), jdim, kdim, idim);
                    FortranArray1DRef<double> eomui_ref(&wk(iwk1), inmx);
                    FortranArray1DRef<double> fbl_ref(&wk(iwk2), inmx);
                    FortranArray1DRef<double> rhon_ref(&wk(iwk3), inmx);
                    FortranArray1DRef<double> amun_ref(&wk(iwk4), inmx);
                    FortranArray1DRef<double> vortn_ref(&wk(iwk5), inmx);
                    FortranArray1DRef<double> disn_ref(&wk(iwk6), inmx);
                    FortranArray1DRef<double> utot_ref(&wk(iwk7), inmx);
                    FortranArray1DRef<double> eomuo_ref(&wk(iwk8), inmx);
                    FortranArray1DRef<double> damp_ref(&wk(iwk9), inmx);
                    blomax_ns::blomax(jdim, kdim, idim, q, qi04d, qj04d, qk04d,
                                      vor3d, snj0, snk0, sni0,
                                      blnum, xkb, blnum,
                                      vist3d, eoms3d, ipw, inmx,
                                      eomui_ref, fbl_ref, rhon_ref, amun_ref,
                                      vortn_ref, disn_ref, utot_ref, eomuo_ref, damp_ref,
                                      nblt, x3d, y3d, z3d, blank, iover,
                                      bci, bcj, bck,
                                      nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl);
                }
            } else if (ivmx == 4) {
                // Baldwin-Barth
                iwk1 = 1;
                if (iturbord == 1) {
                    iwk5 = iwk1 + (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                    iex = 0;
                    iex2 = -i2d;
                } else {
                    iwk5 = iwk1 + (jdim+3)*(kdim+3)*(idim+3-(4*i2d));
                    iex = 1;
                    iex2 = 1 - (2*i2d);
                }
                inmx = (kdim-1)*(jdim-1);
                iwk6  = iwk5  + inmx;
                iwk7  = iwk6  + inmx;
                iwk8  = iwk7  + inmx;
                iwk9  = iwk8  + inmx;
                iwk10 = iwk9  + inmx;
                inmx = (jdim-1)*(kdim-1);
                iwk11 = iwk10 + inmx;
                iwk12 = iwk11 + inmx;
                iwk13 = iwk12 + inmx;
                iwk14 = iwk13 + inmx;
                iwk15 = iwk14 + inmx;
                inmx = (kdim-1)*(idim-1);
                iwk16 = iwk15 + inmx;
                iwk17 = iwk16 + inmx;
                iwk18 = iwk17 + inmx;
                iwk19 = iwk18 + inmx;
                nroom = nwork - (iwk19 + inmx);
                if (nroom < 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " not enough memory for B-B turb model.");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " nroom=%12d", nroom);
                    { int m1=-1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
                }
                {
                    FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                    FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> x3d(&x(1), jdim, kdim, idim);
                    FortranArray3DRef<double> y3d(&y(1), jdim, kdim, idim);
                    FortranArray3DRef<double> z3d(&z(1), jdim, kdim, idim);
                    // turre = wk(iwk1) as 3D
                    int turre_sz = (iturbord==1) ? (jdim+1)*(kdim+1)*(idim+1-(2*i2d))
                                                 : (jdim+3)*(kdim+3)*(idim+3-(4*i2d));
                    FortranArray3DRef<double> turre3d(&wk(iwk1), turre_sz, 1, 1);
                    FortranArray3DRef<double> vor3d(&res(1,1), jdim, kdim, idim1);
                    FortranArray3DRef<double> damp1_3d(&res(1,2), jdim, kdim, idim1);
                    FortranArray3DRef<double> damp2_3d(&res(1,3), jdim, kdim, idim1);
                    FortranArray3DRef<double> timestp3d(&res(1,4), jdim, kdim, idim1);
                    FortranArray3DRef<double> fnu3d(&res(1,5), jdim, kdim, idim1);
                    FortranArray3DRef<double> wk5_3d(&wk(iwk5), (kdim-1)*(jdim-1), 1, 1);
                    // 2D work arrays for barth3d
                    int bx_sz = (kdim-1)*(jdim-1);
                    int by_sz = (jdim-1)*(kdim-1);
                    int bz_sz = (kdim-1)*(idim-1);
                    FortranArray2DRef<double> bx_ref(&wk(iwk6),  kdim-1, jdim-1);
                    FortranArray2DRef<double> cx_ref(&wk(iwk7),  kdim-1, jdim-1);
                    FortranArray2DRef<double> dx_ref(&wk(iwk8),  kdim-1, jdim-1);
                    FortranArray2DRef<double> fx_ref(&wk(iwk9),  kdim-1, jdim-1);
                    FortranArray2DRef<double> wx_ref(&wk(iwk10), kdim-1, jdim-1);
                    FortranArray2DRef<double> by_ref(&wk(iwk11), jdim-1, kdim-1);
                    FortranArray2DRef<double> cy_ref(&wk(iwk12), jdim-1, kdim-1);
                    FortranArray2DRef<double> dy_ref(&wk(iwk13), jdim-1, kdim-1);
                    FortranArray2DRef<double> fy_ref(&wk(iwk14), jdim-1, kdim-1);
                    FortranArray2DRef<double> wy_ref(&wk(iwk15), jdim-1, kdim-1);
                    FortranArray2DRef<double> bz_ref(&wk(iwk16), kdim-1, idim-1);
                    FortranArray2DRef<double> cz_ref(&wk(iwk17), kdim-1, idim-1);
                    FortranArray2DRef<double> dz_ref(&wk(iwk18), kdim-1, idim-1);
                    FortranArray2DRef<double> fz_ref(&wk(iwk19), kdim-1, idim-1);
                    // workz uses wk(iwk5) as scratch (same memory as bx area, reused)
                    FortranArray2DRef<double> wz_ref(&wk(iwk5), kdim-1, idim-1);
                    double sumn1_d = (double)sumn(1);
                    double sumn2_d = (double)sumn(2);
                    barth3d_ns::barth3d(jdim, kdim, idim, q, sj4d, sk4d, si4d,
                                        vol3d, dtj3d, x3d, y3d, z3d, vist3d,
                                        vor3d, snk0, snj0, xib, xkb, turre3d,
                                        damp1_3d, damp2_3d, timestp3d, fnu3d,
                                        bx_ref, cx_ref, dx_ref, fx_ref, wx_ref,
                                        by_ref, cy_ref, dy_ref, fy_ref, wy_ref,
                                        bz_ref, cz_ref, dz_ref, fz_ref, wz_ref,
                                        ntime, tj0, tk0, ti0, nbl, blnum, blank, iover,
                                        sumn1_d, sumn2_d, negn(1), negn(2),
                                        xib2, volj0, volk0, voli0,
                                        nou, bou, nbuf, ibufdim, iex, iex2, nummem);
                    sumn(1) = sumn1_d;
                    sumn(2) = sumn2_d;
                }


            } else if (ivmx == 5) {
                // Spalart
                iwk1 = 1;
                if (iturbord == 1) {
                    iwk3 = iwk1 + (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                    iex = 0;
                    iex2 = -i2d;
                } else {
                    iwk3 = iwk1 + (jdim+3)*(kdim+3)*(idim+3-(4*i2d));
                    iex = 1;
                    iex2 = 1 - (2*i2d);
                }
                iwk4 = iwk3 + (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                iex3 = -i2d;
                inmx = (kdim-1)*(jdim-1);
                iwk5  = iwk4  + inmx;
                iwk6  = iwk5  + inmx;
                iwk7  = iwk6  + inmx;
                iwk8  = iwk7  + inmx;
                iwk9  = iwk8  + inmx;
                inmx = (jdim-1)*(kdim-1);
                iwk10 = iwk9  + inmx;
                iwk11 = iwk10 + inmx;
                iwk12 = iwk11 + inmx;
                iwk13 = iwk12 + inmx;
                iwk14 = iwk13 + inmx;
                inmx = (kdim-1)*(idim-1);
                iwk15 = iwk14 + inmx;
                iwk16 = iwk15 + inmx;
                iwk17 = iwk16 + inmx;
                iwk18 = iwk17 + inmx;
                iwk19 = iwk18 + inmx;
                if (isarc2d == 1) {
                    isarcnum = 4;
                    inmx = (jdim+1)*(kdim+1)*(idim-1)*isarcnum;
                } else if (isarc3d == 1) {
                    isarcnum = 6;
                    inmx = (jdim+1)*(kdim+1)*(idim-1)*isarcnum;
                } else {
                    isarcnum = 1;
                    inmx = 0;
                }
                nroom = nwork - (iwk19 + inmx);
                if (nroom < 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " not enough memory for Spalart turb model.");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " nroom=%12d", nroom);
                    { int m1=-1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
                }
                {
                    FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                    FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> x3d(&x(1), jdim, kdim, idim);
                    FortranArray3DRef<double> y3d(&y(1), jdim, kdim, idim);
                    FortranArray3DRef<double> z3d(&z(1), jdim, kdim, idim);
                    FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                    FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                    FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                    int tursav_sz = (iturbord==1) ? (jdim+1)*(kdim+1)*(idim+1-(2*i2d))
                                                  : (jdim+3)*(kdim+3)*(idim+3-(4*i2d));
                    // Fortran spalart arg storage (call in resid.F):
                    //   tursav=xib, turre=wk(iwk1), damp1=res(1,2), timestp=res(1,3),
                    //   fnu=wk(iwk3), bx..workz=wk(iwk4)..wk(iwk18), vx=wk(iwk19),
                    //   fd=res(1,4).  (The earlier C++ shifted every wk view by one
                    //   slot and put tursav at iwk1/fd at iwk19 — never exercised
                    //   because grdflat5 uses the two-equation model, not SA.)
                    // turre(0-iex:jdim+iex, 0-iex:kdim+iex, 0-iex2:idim+iex2), base wk(iwk1)
                    FortranArray3DRef<double> turre3d(&wk(iwk1),
                        jdim + 2*iex + 1, kdim + 2*iex + 1, idim + 2*iex2 + 1,
                        -iex, -iex, -iex2);
                    FortranArray3DRef<double> vor3d(&res(1,1), jdim, kdim, idim1);
                    FortranArray3DRef<double> damp1_3d(&res(1,2), jdim, kdim, idim1);
                    FortranArray3DRef<double> timestp3d(&res(1,3), jdim, kdim, idim1);
                    // fnu(0:jdim, 0:kdim, 0-iex3:idim+iex3), base wk(iwk3)
                    FortranArray3DRef<double> fnu3d(&wk(iwk3),
                        jdim + 1, kdim + 1, idim + 2*iex3 + 1,
                        0, 0, -iex3);
                    FortranArray2DRef<double> bx_ref(&wk(iwk4),  kdim-1, jdim-1);
                    FortranArray2DRef<double> cx_ref(&wk(iwk5),  kdim-1, jdim-1);
                    FortranArray2DRef<double> dx_ref(&wk(iwk6),  kdim-1, jdim-1);
                    FortranArray2DRef<double> fx_ref(&wk(iwk7),  kdim-1, jdim-1);
                    FortranArray2DRef<double> wx_ref(&wk(iwk8),  kdim-1, jdim-1);
                    FortranArray2DRef<double> by_ref(&wk(iwk9),  jdim-1, kdim-1);
                    FortranArray2DRef<double> cy_ref(&wk(iwk10), jdim-1, kdim-1);
                    FortranArray2DRef<double> dy_ref(&wk(iwk11), jdim-1, kdim-1);
                    FortranArray2DRef<double> fy_ref(&wk(iwk12), jdim-1, kdim-1);
                    FortranArray2DRef<double> wy_ref(&wk(iwk13), jdim-1, kdim-1);
                    FortranArray2DRef<double> bz_ref(&wk(iwk14), kdim-1, idim-1);
                    FortranArray2DRef<double> cz_ref(&wk(iwk15), kdim-1, idim-1);
                    FortranArray2DRef<double> dz_ref(&wk(iwk16), kdim-1, idim-1);
                    FortranArray2DRef<double> fz_ref(&wk(iwk17), kdim-1, idim-1);
                    FortranArray2DRef<double> wz_ref(&wk(iwk18), kdim-1, idim-1);
                    // fd = res(1,4)
                    FortranArray3DRef<double> fd3d(&res(1,4), jdim, kdim, idim1);
                    // vx(0:jdim,0:kdim,idim-1,isarcnum), base wk(iwk19)
                    FortranArray4DRef<double> vx4d(&wk(iwk19),
                        jdim + 1, kdim + 1, idim - 1, (isarcnum > 0 ? isarcnum : 1),
                        0, 0, 1, 1);
                    double sumn1_d = (double)sumn(1);
                    double sumn2_d = (double)sumn(2);
                    // snk0 arrives dimensioned (jdim,idim,3) from mgblk; the wall-distance
                    // (smin) array underneath is (jdim-1,kdim-1,idim-1). Re-base like the
                    // other turb models (blomax/twoeqn smin3d) so spalart indexes correctly.
                    FortranArray3DRef<double> smin3d_sa(&snk0(1,1,1), jdim-1, kdim-1, idim-1);
                    spalart_ns::spalart(jdim, kdim, idim, q, sj4d, sk4d, si4d,
                                        vol3d, dtj3d, x3d, y3d, z3d, vist3d,
                                        vor3d, smin3d_sa, xib, turre3d,
                                        damp1_3d, timestp3d, fnu3d,
                                        bx_ref, cx_ref, dx_ref, fx_ref, wx_ref,
                                        by_ref, cy_ref, dy_ref, fy_ref, wy_ref,
                                        bz_ref, cz_ref, dz_ref, fz_ref, wz_ref,
                                        ntime, tj0, tk0, ti0, nbl,
                                        qj04d, qk04d, qi04d, blank, iover,
                                        sumn1_d, sumn2_d, negn(1), negn(2),
                                        xib2, volj0, volk0, voli0,
                                        nou, bou, nbuf, ibufdim, iex, iex2, iex3,
                                        bcj, bck, bci, ux, vx4d, nummem, fd3d, isarcnum);
                    sumn(1) = sumn1_d;
                    sumn(2) = sumn2_d;
                }


            } else if (ivmx >= 6 && ivmx <= 16) {
                // two-equation turb model
                iwk1 = 1;
                if (iturbord == 1) {
                    iwk4 = iwk1 + (jdim+1)*(kdim+1)*(idim+1-(2*i2d))*2;
                    iex = 0;
                    iex2 = -i2d;
                } else {
                    iwk4 = iwk1 + (jdim+3)*(kdim+3)*(idim+3-(4*i2d))*2;
                    iex = 1;
                    iex2 = 1 - (2*i2d);
                }
                iwk5 = iwk4 + (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                iex3 = -i2d;
                inmx = (kdim-1)*(jdim-1);
                iwk6  = iwk5  + inmx; iwk7  = iwk6  + inmx; iwk8  = iwk7  + inmx;
                iwk9  = iwk8  + inmx; iwk10 = iwk9  + inmx; iwk11 = iwk10 + inmx;
                iwk12 = iwk11 + inmx; iwk13 = iwk12 + inmx; iwk14 = iwk13 + inmx;
                inmx = (jdim-1)*(kdim-1);
                iwk15 = iwk14 + inmx; iwk16 = iwk15 + inmx; iwk17 = iwk16 + inmx;
                iwk18 = iwk17 + inmx; iwk19 = iwk18 + inmx; iwk20 = iwk19 + inmx;
                iwk21 = iwk20 + inmx; iwk22 = iwk21 + inmx; iwk23 = iwk22 + inmx;
                inmx = (kdim-1)*(idim-1);
                iwk24 = iwk23 + inmx; iwk25 = iwk24 + inmx; iwk26 = iwk25 + inmx;
                iwk27 = iwk26 + inmx; iwk28 = iwk27 + inmx; iwk29 = iwk28 + inmx;
                iwk30 = iwk29 + inmx; iwk31 = iwk30 + inmx; iwk32 = iwk31 + inmx;
                inmx = (jdim-1)*(kdim-1)*(idim-1)*2;
                iwk33 = iwk32 + inmx;
                inmx = (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                iwk34 = iwk33 + inmx;
                if (ivmx == 15) { inmx = jdim; } else { inmx = 0; }
                iwk35 = iwk34 + inmx; iwk36 = iwk35 + inmx; iwk37 = iwk36 + inmx;
                iwk38 = iwk37 + inmx; iwk39 = iwk38 + inmx; iwk40 = iwk39 + inmx;
                iwk41 = iwk40 + inmx; iwk42 = iwk41 + inmx; iwk43 = iwk42 + inmx;
                iwk44 = iwk43 + inmx; iwk45 = iwk44 + inmx; iwk46 = iwk45 + inmx;
                iwk47 = iwk46 + inmx; iwk48 = iwk47 + inmx; iwk49 = iwk48 + inmx;
                iwk50 = iwk49 + inmx; iwk51 = iwk50 + inmx; iwk52 = iwk51 + inmx;
                iwk53 = iwk52 + inmx; iwk54 = iwk53 + inmx; iwk55 = iwk54 + inmx;
                iwk56 = iwk55 + inmx; iwk57 = iwk56 + inmx; iwk58 = iwk57 + inmx;
                iwk59 = iwk58 + inmx; iwk60 = iwk59 + inmx; iwk61 = iwk60 + inmx;
                iwk62 = iwk61 + inmx; iwk63 = iwk62 + inmx; iwk64 = iwk63 + inmx;
                iwk65 = iwk64 + inmx; iwk66 = iwk65 + inmx; iwk67 = iwk66 + inmx;
                iwk68 = iwk67 + inmx; iwk69 = iwk68 + inmx; iwk70 = iwk69 + inmx;
                iwk71 = iwk70 + inmx; iwk72 = iwk71 + inmx; iwk73 = iwk72 + inmx;
                iwk74 = iwk73 + inmx; iwk75 = iwk74 + inmx; iwk76 = iwk75 + inmx;
                if (ivmx == 15) { inmx = (jdim+1)*(kdim+1)*(idim+1); } else { inmx = 0; }
                iwk77 = iwk76 + inmx; iwk78 = iwk77 + inmx; iwk79 = iwk78 + inmx;
                if (ivmx == 15) { inmx = (jdim+1)*(kdim+1)*(idim+1)*2; } else { inmx = 0; }
                iwk80 = iwk79 + inmx;
                if (ieasmcc2d == 1) {
                    iccnum = 4; inmx = (jdim+1)*(kdim+1)*(idim-1)*4;
                } else if (isstrc == 2) {
                    iccnum = 6; inmx = (jdim+1)*(kdim+1)*(idim-1)*6;
                } else if (ivmx == 16) {
                    iccnum = 4; inmx = (jdim+1)*(kdim+1)*(idim-1)*4;
                } else {
                    iccnum = 1; inmx = 0;
                }
                iwk81 = iwk80 + inmx;
                if (ides != 0) { inmx = (jdim-1)*(kdim-1)*(idim-1); } else { inmx = 0; }
                iwk82 = iwk81 + inmx;
                if (ides == 3) { inmx = (jdim-1)*(kdim-1)*(idim-1); } else { inmx = 0; }
                nroom = nwork - (iwk82 + inmx);
                if (nroom < 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " not enough memory for 2-eqn turb model.");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " nroom=%12d", nroom);
                    { int m1=-1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
                }


                {
                    FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                    FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> x3d(&x(1), jdim, kdim, idim);
                    FortranArray3DRef<double> y3d(&y(1), jdim, kdim, idim);
                    FortranArray3DRef<double> z3d(&z(1), jdim, kdim, idim);
                    FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                    FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                    FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                    // turre(0-iex:jdim+iex, 0-iex:kdim+iex, 0-iex2:idim+iex2, 2)
                    // base wk(iwk1) = turre(-iex,-iex,-iex2,1)
                    FortranArray4DRef<double> turre4d(&wk(iwk1),
                        jdim + 2*iex + 1, kdim + 2*iex + 1, idim + 2*iex2 + 1, 2,
                        -iex, -iex, -iex2, 1);
                    // fnu(0:jdim, 0:kdim, 0-iex3:idim+iex3), base wk(iwk4)
                    FortranArray3DRef<double> fnu3d(&wk(iwk4),
                        jdim + 1, kdim + 1, idim + 2*iex3 + 1,
                        0, 0, -iex3);
                    FortranArray3DRef<double> vor3d(&res(1,1), jdim, kdim, idim1);
                    FortranArray3DRef<double> smin3d(&snk0(1,1,1), jdim-1, kdim-1, idim-1);
                    FortranArray3DRef<double> damp1_3d(&res(1,2), jdim, kdim, idim1);
                    FortranArray3DRef<double> blend3d(&res(1,3), jdim, kdim, idim1);
                    FortranArray3DRef<double> timestp3d(&res(1,4), jdim, kdim, idim1);
                    FortranArray3DRef<double> wrks3d(&res(1,5), jdim, kdim, idim1);
                    // Fortran twoeqn arg mapping: bx=iwk5, bx2=iwk6, ... workx=iwk13,
                    //   by=iwk14 ... worky=iwk22, bz=iwk23 ... fz2=iwk30, workz=iwk31,
                    //   rhside=iwk32, v3dtmp=iwk33.  (The coefficient views were all
                    //   shifted +1, which made workz alias rhside at iwk32.)
                    FortranArray2DRef<double> bx_ref(&wk(iwk5),  kdim-1, jdim-1);
                    FortranArray2DRef<double> bx2_ref(&wk(iwk6), kdim-1, jdim-1);
                    FortranArray2DRef<double> cx_ref(&wk(iwk7),  kdim-1, jdim-1);
                    FortranArray2DRef<double> cx2_ref(&wk(iwk8), kdim-1, jdim-1);
                    FortranArray2DRef<double> dx_ref(&wk(iwk9),  kdim-1, jdim-1);
                    FortranArray2DRef<double> dx2_ref(&wk(iwk10),kdim-1, jdim-1);
                    FortranArray2DRef<double> fx_ref(&wk(iwk11), kdim-1, jdim-1);
                    FortranArray2DRef<double> fx2_ref(&wk(iwk12),kdim-1, jdim-1);
                    FortranArray2DRef<double> wx_ref(&wk(iwk13), kdim-1, jdim-1);
                    FortranArray2DRef<double> by_ref(&wk(iwk14), jdim-1, kdim-1);
                    FortranArray2DRef<double> by2_ref(&wk(iwk15),jdim-1, kdim-1);
                    FortranArray2DRef<double> cy_ref(&wk(iwk16), jdim-1, kdim-1);
                    FortranArray2DRef<double> cy2_ref(&wk(iwk17),jdim-1, kdim-1);
                    FortranArray2DRef<double> dy_ref(&wk(iwk18), jdim-1, kdim-1);
                    FortranArray2DRef<double> dy2_ref(&wk(iwk19),jdim-1, kdim-1);
                    FortranArray2DRef<double> fy_ref(&wk(iwk20), jdim-1, kdim-1);
                    FortranArray2DRef<double> fy2_ref(&wk(iwk21),jdim-1, kdim-1);
                    FortranArray2DRef<double> wy_ref(&wk(iwk22), jdim-1, kdim-1);
                    FortranArray2DRef<double> bz_ref(&wk(iwk23), kdim-1, idim-1);
                    FortranArray2DRef<double> bz2_ref(&wk(iwk24),kdim-1, idim-1);
                    FortranArray2DRef<double> cz_ref(&wk(iwk25), kdim-1, idim-1);
                    FortranArray2DRef<double> cz2_ref(&wk(iwk26),kdim-1, idim-1);
                    FortranArray2DRef<double> dz_ref(&wk(iwk27), kdim-1, idim-1);
                    FortranArray2DRef<double> dz2_ref(&wk(iwk28),kdim-1, idim-1);
                    FortranArray2DRef<double> fz_ref(&wk(iwk29), kdim-1, idim-1);
                    FortranArray2DRef<double> fz2_ref(&wk(iwk30),kdim-1, idim-1);
                    FortranArray2DRef<double> wz_ref(&wk(iwk31), kdim-1, idim-1);
                    int rhside_sz = (jdim-1)*(kdim-1)*(idim-1)*2;
                    // rhside(jdim-1,kdim-1,idim-1,2) — must use true 4D dims, not a
                    // flat (rhside_sz,1,1,1) view (that aliases k-slot of cell N onto
                    // omega-slot of cell N+1, corrupting the omega equation RHS).
                    FortranArray4DRef<double> rhside4d(&wk(iwk32), jdim-1, kdim-1, idim-1, 2);
                    (void)rhside_sz;
                    int v3d_sz = (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                    // v3dtmp(0:jdim, 0:kdim, 0-iex3:idim+iex3) — needs the true
                    // ghost-offset 3D shape, not a flat (v3d_sz,1,1) view (which
                    // mislocates every k>1 / i>1 cell and has no ghost cells,
                    // corrupting the turbulence diffusion coeffs at all boundaries).
                    FortranArray3DRef<double> v3dtmp3d(&wk(iwk33),
                        jdim + 1, kdim + 1, idim + 2*iex3 + 1,
                        0, 0, -iex3);
                    (void)v3d_sz;
                    int grd_sz_x = (kdim-1)*(jdim-1);
                    int grd_sz_y = (jdim-1)*(kdim-1);
                    int grd_sz_z = (kdim-1)*(idim-1);
                    FortranArray1DRef<double> dkdx_r(&wk(iwk35), std::max(1,grd_sz_x));
                    FortranArray1DRef<double> dkdy_r(&wk(iwk36), std::max(1,grd_sz_x));
                    FortranArray1DRef<double> dkdz_r(&wk(iwk37), std::max(1,grd_sz_x));
                    FortranArray1DRef<double> dzdx_r(&wk(iwk38), std::max(1,grd_sz_x));
                    FortranArray1DRef<double> dzdy_r(&wk(iwk39), std::max(1,grd_sz_y));
                    FortranArray1DRef<double> dzdz_r(&wk(iwk40), std::max(1,grd_sz_y));
                    FortranArray1DRef<double> dw1dx_r(&wk(iwk41),std::max(1,grd_sz_y));
                    FortranArray1DRef<double> dw1dy_r(&wk(iwk42),std::max(1,grd_sz_y));
                    FortranArray1DRef<double> dw1dz_r(&wk(iwk43),std::max(1,grd_sz_y));
                    FortranArray1DRef<double> dw2dx_r(&wk(iwk44),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dw2dy_r(&wk(iwk45),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dw2dz_r(&wk(iwk46),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dw3dx_r(&wk(iwk47),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dw3dy_r(&wk(iwk48),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dw3dz_r(&wk(iwk49),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> drdx_r(&wk(iwk50), std::max(1,grd_sz_z));
                    FortranArray1DRef<double> drdy_r(&wk(iwk51), std::max(1,grd_sz_z));
                    FortranArray1DRef<double> drdz_r(&wk(iwk52), std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dpdx_r(&wk(iwk53), std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dpdy_r(&wk(iwk54), std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dpdz_r(&wk(iwk55), std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dkdj_r(&wk(iwk56), std::max(1,grd_sz_x));
                    FortranArray1DRef<double> dkdk_r(&wk(iwk57), std::max(1,grd_sz_x));
                    FortranArray1DRef<double> dkdi_r(&wk(iwk58), std::max(1,grd_sz_x));
                    FortranArray1DRef<double> dzdj_r(&wk(iwk59), std::max(1,grd_sz_x));
                    FortranArray1DRef<double> dzdk_r(&wk(iwk60), std::max(1,grd_sz_y));
                    FortranArray1DRef<double> dzdi_r(&wk(iwk61), std::max(1,grd_sz_y));
                    FortranArray1DRef<double> dw1dj_r(&wk(iwk62),std::max(1,grd_sz_y));
                    FortranArray1DRef<double> dw1dk_r(&wk(iwk63),std::max(1,grd_sz_y));
                    FortranArray1DRef<double> dw1di_r(&wk(iwk64),std::max(1,grd_sz_y));
                    FortranArray1DRef<double> dw2dj_r(&wk(iwk65),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dw2dk_r(&wk(iwk66),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dw2di_r(&wk(iwk67),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dw3dj_r(&wk(iwk68),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dw3dk_r(&wk(iwk69),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dw3di_r(&wk(iwk70),std::max(1,grd_sz_z));
                    FortranArray1DRef<double> drdj_r(&wk(iwk71), std::max(1,grd_sz_z));
                    FortranArray1DRef<double> drdk_r(&wk(iwk72), std::max(1,grd_sz_z));
                    FortranArray1DRef<double> drdi_r(&wk(iwk73), std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dpdj_r(&wk(iwk74), std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dpdk_r(&wk(iwk75), std::max(1,grd_sz_z));
                    FortranArray1DRef<double> dpdi_r(&wk(iwk76), std::max(1,grd_sz_z));
                    int w1_sz = (ivmx==15) ? (jdim+1)*(kdim+1)*(idim+1) : 1;
                    FortranArray3DRef<double> w1_3d(&wk(iwk77), w1_sz, 1, 1);
                    FortranArray3DRef<double> w2_3d(&wk(iwk78), w1_sz, 1, 1);
                    FortranArray3DRef<double> w3_3d(&wk(iwk79), w1_sz, 1, 1);
                    int srce_sz = (ivmx==15) ? (jdim+1)*(kdim+1)*(idim+1)*2 : 1;
                    FortranArray4DRef<double> srce4d(&wk(iwk80), srce_sz, 1, 1, 1);
                    int xlscale_sz = (iccnum > 1) ? (jdim+1)*(kdim+1)*(idim-1)*iccnum : 1;
                    FortranArray3DRef<double> xlscale3d(&wk(iwk81), xlscale_sz, 1, 1);
                    int fdsav_sz = (ides != 0) ? (jdim-1)*(kdim-1)*(idim-1) : 1;
                    FortranArray3DRef<double> fdsav3d(&wk(iwk82), fdsav_sz, 1, 1);
                    double sumn1_d = (double)sumn(1);
                    double sumn2_d = (double)sumn(2);
                    twoeqn_ns::twoeqn(jdim, kdim, idim, q, sj4d, sk4d, si4d,
                                      vol3d, dtj3d, x3d, y3d, z3d, vist3d,
                                      vor3d, smin3d, xib, turre4d,
                                      damp1_3d, blend3d, timestp3d, fnu3d,
                                      bx_ref, bx2_ref, cx_ref, cx2_ref,
                                      dx_ref, dx2_ref, fx_ref, fx2_ref, wx_ref,
                                      by_ref, by2_ref, cy_ref, cy2_ref,
                                      dy_ref, dy2_ref, fy_ref, fy2_ref, wy_ref,
                                      bz_ref, bz2_ref, cz_ref, cz2_ref,
                                      dz_ref, dz2_ref, fz_ref, fz2_ref, wz_ref,
                                      ntime, tj0, tk0, ti0, nbl,
                                      qj04d, qk04d, qi04d, vj0, vk0, vi0,
                                      blank, iover, sumn1_d, sumn2_d, negn(1), negn(2),
                                      ux, rhside4d, xib2, v3dtmp3d, cmuv,
                                      bcj, bck, bci,
                                      nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                                      ibcinfo, jbcinfo, kbcinfo,
                                      maxbl, maxseg, volj0, volk0, voli0,
                                      nou, bou, nbuf, ibufdim, iex, iex2, iex3,
                                      dkdx_r, dkdy_r, dkdz_r, dzdx_r, dzdy_r, dzdz_r,
                                      dw1dx_r, dw1dy_r, dw1dz_r, dw2dx_r, dw2dy_r, dw2dz_r,
                                      dw3dx_r, dw3dy_r, dw3dz_r, drdx_r, drdy_r, drdz_r,
                                      dpdx_r, dpdy_r, dpdz_r,
                                      dkdj_r, dkdk_r, dkdi_r, dzdj_r, dzdk_r, dzdi_r,
                                      dw1dj_r, dw1dk_r, dw1di_r, dw2dj_r, dw2dk_r, dw2di_r,
                                      dw3dj_r, dw3dk_r, dw3di_r, drdj_r, drdk_r, drdi_r,
                                      dpdj_r, dpdk_r, dpdi_r,
                                      w1_3d, w2_3d, w3_3d, srce4d, srce4d,
                                      xlscale3d, fdsav3d, nummem, iccnum);
                    sumn(1) = sumn1_d;
                    sumn(2) = sumn2_d;
                }
            } else if (ivmx == 30) {
                // three-equation turb model
                iwk1 = 1;
                if (iturbord == 1) {
                    iwk4 = iwk1 + (jdim+1)*(kdim+1)*(idim+1-(2*i2d))*3;
                    iex = 0; iex2 = -i2d;
                } else {
                    iwk4 = iwk1 + (jdim+3)*(kdim+3)*(idim+3-(4*i2d))*3;
                    iex = 1; iex2 = 1-(2*i2d);
                }
                iwk5 = iwk4 + (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                iex3 = -i2d;
                inmx = (kdim-1)*(jdim-1);
                iwk6  = iwk5  + inmx; iwk7  = iwk6  + inmx; iwk8  = iwk7  + inmx;
                iwk9  = iwk8  + inmx; iwk10 = iwk9  + inmx; iwk11 = iwk10 + inmx;
                iwk12 = iwk11 + inmx; iwk13 = iwk12 + inmx; iwk14 = iwk13 + inmx;
                inmx = (jdim-1)*(kdim-1);
                iwk15 = iwk14 + inmx; iwk16 = iwk15 + inmx; iwk17 = iwk16 + inmx;
                iwk18 = iwk17 + inmx; iwk19 = iwk18 + inmx; iwk20 = iwk19 + inmx;
                iwk21 = iwk20 + inmx; iwk22 = iwk21 + inmx; iwk23 = iwk22 + inmx;
                inmx = (kdim-1)*(idim-1);
                iwk24 = iwk23 + inmx; iwk25 = iwk24 + inmx; iwk26 = iwk25 + inmx;
                iwk27 = iwk26 + inmx; iwk28 = iwk27 + inmx; iwk29 = iwk28 + inmx;
                iwk30 = iwk29 + inmx; iwk31 = iwk30 + inmx; iwk32 = iwk31 + inmx;
                inmx = (jdim-1)*(kdim-1)*(idim-1)*3;
                iwk33 = iwk32 + inmx;
                inmx = (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                iwk34 = iwk33 + inmx;
                inmx = (kdim-1)*(jdim-1);
                iwk35 = iwk34 + inmx; iwk36 = iwk35 + inmx; iwk37 = iwk36 + inmx; iwk38 = iwk37 + inmx;
                inmx = (jdim-1)*(kdim-1);
                iwk39 = iwk38 + inmx; iwk40 = iwk39 + inmx; iwk41 = iwk40 + inmx; iwk42 = iwk41 + inmx;
                inmx = (kdim-1)*(idim-1);
                iwk43 = iwk42 + inmx; iwk44 = iwk43 + inmx; iwk45 = iwk44 + inmx;
                nroom = nwork - (iwk45 + inmx);
                if (nroom < 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " not enough memory for 3-eqn turb model.");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " nroom=%12d", nroom);
                    { int m1=-1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
                }
                {
                    FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                    FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> x3d(&x(1), jdim, kdim, idim);
                    FortranArray3DRef<double> y3d(&y(1), jdim, kdim, idim);
                    FortranArray3DRef<double> z3d(&z(1), jdim, kdim, idim);
                    FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                    FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                    FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                    int zksav_sz = (iturbord==1) ? (jdim+1)*(kdim+1)*(idim+1-(2*i2d))*3
                                                 : (jdim+3)*(kdim+3)*(idim+3-(4*i2d))*3;
                    FortranArray4DRef<double> zksav4d(&wk(iwk1), zksav_sz, 1, 1, 1);
                    int fnu_sz = (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                    FortranArray3DRef<double> fnu3d(&wk(iwk5), fnu_sz, 1, 1);
                    FortranArray3DRef<double> vor3d(&res(1,1), jdim, kdim, idim1);
                    FortranArray3DRef<double> smin3d(&snk0(1,1,1), jdim-1, kdim-1, idim-1);
                    FortranArray3DRef<double> damp1_3d(&res(1,2), jdim, kdim, idim1);
                    FortranArray3DRef<double> blend3d(&res(1,3), jdim, kdim, idim1);
                    FortranArray3DRef<double> timestp3d(&res(1,4), jdim, kdim, idim1);
                    FortranArray3DRef<double> wrks3d(&res(1,5), jdim, kdim, idim1);
                    FortranArray2DRef<double> bx_ref(&wk(iwk6),  kdim-1, jdim-1);
                    FortranArray2DRef<double> bx2_ref(&wk(iwk7), kdim-1, jdim-1);
                    FortranArray2DRef<double> cx_ref(&wk(iwk8),  kdim-1, jdim-1);
                    FortranArray2DRef<double> cx2_ref(&wk(iwk9), kdim-1, jdim-1);
                    FortranArray2DRef<double> dx_ref(&wk(iwk10), kdim-1, jdim-1);
                    FortranArray2DRef<double> dx2_ref(&wk(iwk11),kdim-1, jdim-1);
                    FortranArray2DRef<double> fx_ref(&wk(iwk12), kdim-1, jdim-1);
                    FortranArray2DRef<double> fx2_ref(&wk(iwk13),kdim-1, jdim-1);
                    FortranArray2DRef<double> wx_ref(&wk(iwk14), kdim-1, jdim-1);
                    FortranArray2DRef<double> by_ref(&wk(iwk15), jdim-1, kdim-1);
                    FortranArray2DRef<double> by2_ref(&wk(iwk16),jdim-1, kdim-1);
                    FortranArray2DRef<double> cy_ref(&wk(iwk17), jdim-1, kdim-1);
                    FortranArray2DRef<double> cy2_ref(&wk(iwk18),jdim-1, kdim-1);
                    FortranArray2DRef<double> dy_ref(&wk(iwk19), jdim-1, kdim-1);
                    FortranArray2DRef<double> dy2_ref(&wk(iwk20),jdim-1, kdim-1);
                    FortranArray2DRef<double> fy_ref(&wk(iwk21), jdim-1, kdim-1);
                    FortranArray2DRef<double> fy2_ref(&wk(iwk22),jdim-1, kdim-1);
                    FortranArray2DRef<double> wy_ref(&wk(iwk23), jdim-1, kdim-1);
                    FortranArray2DRef<double> bz_ref(&wk(iwk24), kdim-1, idim-1);
                    FortranArray2DRef<double> bz2_ref(&wk(iwk25),kdim-1, idim-1);
                    FortranArray2DRef<double> cz_ref(&wk(iwk26), kdim-1, idim-1);
                    FortranArray2DRef<double> cz2_ref(&wk(iwk27),kdim-1, idim-1);
                    FortranArray2DRef<double> dz_ref(&wk(iwk28), kdim-1, idim-1);
                    FortranArray2DRef<double> dz2_ref(&wk(iwk29),kdim-1, idim-1);
                    FortranArray2DRef<double> fz_ref(&wk(iwk30), kdim-1, idim-1);
                    FortranArray2DRef<double> fz2_ref(&wk(iwk31),kdim-1, idim-1);
                    FortranArray2DRef<double> wz_ref(&wk(iwk32), kdim-1, idim-1);
                    int rhside_sz = (jdim-1)*(kdim-1)*(idim-1)*3;
                    // rhside(jdim-1,kdim-1,idim-1,3) — true 4D dims (see 2-eqn note above).
                    FortranArray4DRef<double> rhside4d(&wk(iwk32), jdim-1, kdim-1, idim-1, 3);
                    (void)rhside_sz;
                    int v3d_sz = (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                    // v3dtmp(0:jdim, 0:kdim, 0-iex3:idim+iex3) — needs the true
                    // ghost-offset 3D shape, not a flat (v3d_sz,1,1) view (which
                    // mislocates every k>1 / i>1 cell and has no ghost cells,
                    // corrupting the turbulence diffusion coeffs at all boundaries).
                    FortranArray3DRef<double> v3dtmp3d(&wk(iwk33),
                        jdim + 1, kdim + 1, idim + 2*iex3 + 1,
                        0, 0, -iex3);
                    (void)v3d_sz;
                    FortranArray2DRef<double> bx3_ref(&wk(iwk35), kdim-1, jdim-1);
                    FortranArray2DRef<double> cx3_ref(&wk(iwk36), kdim-1, jdim-1);
                    FortranArray2DRef<double> dx3_ref(&wk(iwk37), kdim-1, jdim-1);
                    FortranArray2DRef<double> fx3_ref(&wk(iwk38), kdim-1, jdim-1);
                    FortranArray2DRef<double> by3_ref(&wk(iwk39), jdim-1, kdim-1);
                    FortranArray2DRef<double> cy3_ref(&wk(iwk40), jdim-1, kdim-1);
                    FortranArray2DRef<double> dy3_ref(&wk(iwk41), jdim-1, kdim-1);
                    FortranArray2DRef<double> fy3_ref(&wk(iwk42), jdim-1, kdim-1);
                    FortranArray2DRef<double> bz3_ref(&wk(iwk43), kdim-1, idim-1);
                    FortranArray2DRef<double> cz3_ref(&wk(iwk44), kdim-1, idim-1);
                    FortranArray2DRef<double> dz3_ref(&wk(iwk45), kdim-1, idim-1);
                    FortranArray2DRef<double> fz3_ref(&wk(iwk45), kdim-1, idim-1); // same as dz3
                    double sumn1_d = (double)sumn(1);
                    double sumn2_d = (double)sumn(2);
                    double sumn3_d = (double)sumn(3);
                    threeeqn_ns::threeeqn(jdim, kdim, idim, q, sj4d, sk4d, si4d,
                                          vol3d, dtj3d, x3d, y3d, z3d, vist3d,
                                          vor3d, smin3d, zksav4d, zksav4d,
                                          damp1_3d, blend3d, timestp3d, fnu3d,
                                          bx_ref, bx2_ref, cx_ref, cx2_ref,
                                          dx_ref, dx2_ref, fx_ref, fx2_ref, wx_ref,
                                          by_ref, by2_ref, cy_ref, cy2_ref,
                                          dy_ref, dy2_ref, fy_ref, fy2_ref, wy_ref,
                                          bz_ref, bz2_ref, cz_ref, cz2_ref,
                                          dz_ref, dz2_ref, fz_ref, fz2_ref, wz_ref,
                                          ntime, tj0, tk0, ti0, nbl,
                                          qj04d, qk04d, qi04d, vj0, vk0, vi0,
                                          blank, iover, sumn1_d, sumn2_d, sumn3_d,
                                          negn(1), negn(2), negn(3),
                                          ux, rhside4d, xib2, v3dtmp3d,
                                          bcj, bck, bci,
                                          nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                                          ibcinfo, jbcinfo, kbcinfo,
                                          maxbl, maxseg, volj0, volk0, voli0,
                                          nou, bou, nbuf, ibufdim, iex, iex2, iex3,
                                          bx3_ref, cx3_ref, dx3_ref, fx3_ref,
                                          by3_ref, cy3_ref, dy3_ref, fy3_ref,
                                          bz3_ref, cz3_ref, dz3_ref, fz3_ref, nummem);
                    sumn(1) = sumn1_d;
                    sumn(2) = sumn2_d;
                    sumn(3) = sumn3_d;
                }


            } else if (ivmx == 40) {
                // four-equation turb model
                iwk1 = 1;
                if (iturbord == 1) {
                    iwk4 = iwk1 + (jdim+1)*(kdim+1)*(idim+1-(2*i2d))*4;
                    iex = 0; iex2 = -i2d;
                } else {
                    iwk4 = iwk1 + (jdim+3)*(kdim+3)*(idim+3-(4*i2d))*4;
                    iex = 1; iex2 = 1-(2*i2d);
                }
                iwk5 = iwk4 + (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                iex3 = -i2d;
                inmx = (kdim-1)*(jdim-1);
                iwk6  = iwk5  + inmx; iwk7  = iwk6  + inmx; iwk8  = iwk7  + inmx;
                iwk9  = iwk8  + inmx; iwk10 = iwk9  + inmx; iwk11 = iwk10 + inmx;
                iwk12 = iwk11 + inmx; iwk13 = iwk12 + inmx; iwk14 = iwk13 + inmx;
                inmx = (jdim-1)*(kdim-1);
                iwk15 = iwk14 + inmx; iwk16 = iwk15 + inmx; iwk17 = iwk16 + inmx;
                iwk18 = iwk17 + inmx; iwk19 = iwk18 + inmx; iwk20 = iwk19 + inmx;
                iwk21 = iwk20 + inmx; iwk22 = iwk21 + inmx; iwk23 = iwk22 + inmx;
                inmx = (kdim-1)*(idim-1);
                iwk24 = iwk23 + inmx; iwk25 = iwk24 + inmx; iwk26 = iwk25 + inmx;
                iwk27 = iwk26 + inmx; iwk28 = iwk27 + inmx; iwk29 = iwk28 + inmx;
                iwk30 = iwk29 + inmx; iwk31 = iwk30 + inmx; iwk32 = iwk31 + inmx;
                inmx = (jdim-1)*(kdim-1)*(idim-1)*4;
                iwk33 = iwk32 + inmx;
                inmx = (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                iwk34 = iwk33 + inmx;
                inmx = (kdim-1)*(jdim-1);
                iwk35 = iwk34 + inmx; iwk36 = iwk35 + inmx; iwk37 = iwk36 + inmx; iwk38 = iwk37 + inmx;
                inmx = (jdim-1)*(kdim-1);
                iwk39 = iwk38 + inmx; iwk40 = iwk39 + inmx; iwk41 = iwk40 + inmx; iwk42 = iwk41 + inmx;
                inmx = (kdim-1)*(idim-1);
                iwk43 = iwk42 + inmx; iwk44 = iwk43 + inmx; iwk45 = iwk44 + inmx; iwk46 = iwk45 + inmx;
                inmx = (kdim-1)*(jdim-1);
                iwk47 = iwk46 + inmx; iwk48 = iwk47 + inmx; iwk49 = iwk48 + inmx; iwk50 = iwk49 + inmx;
                inmx = (jdim-1)*(kdim-1);
                iwk51 = iwk50 + inmx; iwk52 = iwk51 + inmx; iwk53 = iwk52 + inmx; iwk54 = iwk53 + inmx;
                inmx = (kdim-1)*(idim-1);
                iwk55 = iwk54 + inmx; iwk56 = iwk55 + inmx; iwk57 = iwk56 + inmx; iwk58 = iwk57 + inmx;
                if (ides != 0) { inmx = (jdim-1)*(kdim-1)*(idim-1); } else { inmx = 0; }
                iwk59 = iwk58 + inmx;
                if (ides == 3) { inmx = (jdim-1)*(kdim-1)*(idim-1); } else { inmx = 0; }
                nroom = nwork - (iwk59 + inmx);
                if (nroom < 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " not enough memory for 4-eqn turb model.");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " nroom=%12d", nroom);
                    { int m1=-1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
                }
                {
                    FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                    FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> x3d(&x(1), jdim, kdim, idim);
                    FortranArray3DRef<double> y3d(&y(1), jdim, kdim, idim);
                    FortranArray3DRef<double> z3d(&z(1), jdim, kdim, idim);
                    FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                    FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                    FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                    int zksav_sz = (iturbord==1) ? (jdim+1)*(kdim+1)*(idim+1-(2*i2d))*4
                                                 : (jdim+3)*(kdim+3)*(idim+3-(4*i2d))*4;
                    FortranArray4DRef<double> zksav4d(&wk(iwk1), zksav_sz, 1, 1, 1);
                    int fnu_sz = (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                    FortranArray3DRef<double> fnu3d(&wk(iwk5), fnu_sz, 1, 1);
                    FortranArray3DRef<double> vor3d(&res(1,1), jdim, kdim, idim1);
                    FortranArray3DRef<double> smin3d(&snk0(1,1,1), jdim-1, kdim-1, idim-1);
                    FortranArray3DRef<double> damp1_3d(&res(1,2), jdim, kdim, idim1);
                    FortranArray3DRef<double> blend3d(&res(1,3), jdim, kdim, idim1);
                    FortranArray3DRef<double> timestp3d(&res(1,4), jdim, kdim, idim1);
                    FortranArray3DRef<double> wrks3d(&res(1,5), jdim, kdim, idim1);
                    FortranArray2DRef<double> bx_ref(&wk(iwk6),  kdim-1, jdim-1);
                    FortranArray2DRef<double> bx2_ref(&wk(iwk7), kdim-1, jdim-1);
                    FortranArray2DRef<double> cx_ref(&wk(iwk8),  kdim-1, jdim-1);
                    FortranArray2DRef<double> cx2_ref(&wk(iwk9), kdim-1, jdim-1);
                    FortranArray2DRef<double> dx_ref(&wk(iwk10), kdim-1, jdim-1);
                    FortranArray2DRef<double> dx2_ref(&wk(iwk11),kdim-1, jdim-1);
                    FortranArray2DRef<double> fx_ref(&wk(iwk12), kdim-1, jdim-1);
                    FortranArray2DRef<double> fx2_ref(&wk(iwk13),kdim-1, jdim-1);
                    FortranArray2DRef<double> wx_ref(&wk(iwk14), kdim-1, jdim-1);
                    FortranArray2DRef<double> by_ref(&wk(iwk15), jdim-1, kdim-1);
                    FortranArray2DRef<double> by2_ref(&wk(iwk16),jdim-1, kdim-1);
                    FortranArray2DRef<double> cy_ref(&wk(iwk17), jdim-1, kdim-1);
                    FortranArray2DRef<double> cy2_ref(&wk(iwk18),jdim-1, kdim-1);
                    FortranArray2DRef<double> dy_ref(&wk(iwk19), jdim-1, kdim-1);
                    FortranArray2DRef<double> dy2_ref(&wk(iwk20),jdim-1, kdim-1);
                    FortranArray2DRef<double> fy_ref(&wk(iwk21), jdim-1, kdim-1);
                    FortranArray2DRef<double> fy2_ref(&wk(iwk22),jdim-1, kdim-1);
                    FortranArray2DRef<double> wy_ref(&wk(iwk23), jdim-1, kdim-1);
                    FortranArray2DRef<double> bz_ref(&wk(iwk24), kdim-1, idim-1);
                    FortranArray2DRef<double> bz2_ref(&wk(iwk25),kdim-1, idim-1);
                    FortranArray2DRef<double> cz_ref(&wk(iwk26), kdim-1, idim-1);
                    FortranArray2DRef<double> cz2_ref(&wk(iwk27),kdim-1, idim-1);
                    FortranArray2DRef<double> dz_ref(&wk(iwk28), kdim-1, idim-1);
                    FortranArray2DRef<double> dz2_ref(&wk(iwk29),kdim-1, idim-1);
                    FortranArray2DRef<double> fz_ref(&wk(iwk30), kdim-1, idim-1);
                    FortranArray2DRef<double> fz2_ref(&wk(iwk31),kdim-1, idim-1);
                    FortranArray2DRef<double> wz_ref(&wk(iwk32), kdim-1, idim-1);
                    int rhside_sz = (jdim-1)*(kdim-1)*(idim-1)*4;
                    // rhside(jdim-1,kdim-1,idim-1,4) — true 4D dims (see 2-eqn note above).
                    FortranArray4DRef<double> rhside4d(&wk(iwk32), jdim-1, kdim-1, idim-1, 4);
                    (void)rhside_sz;
                    int v3d_sz = (jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                    // v3dtmp(0:jdim, 0:kdim, 0-iex3:idim+iex3) — needs the true
                    // ghost-offset 3D shape, not a flat (v3d_sz,1,1) view (which
                    // mislocates every k>1 / i>1 cell and has no ghost cells,
                    // corrupting the turbulence diffusion coeffs at all boundaries).
                    FortranArray3DRef<double> v3dtmp3d(&wk(iwk33),
                        jdim + 1, kdim + 1, idim + 2*iex3 + 1,
                        0, 0, -iex3);
                    (void)v3d_sz;
                    FortranArray2DRef<double> bx3_ref(&wk(iwk35), kdim-1, jdim-1);
                    FortranArray2DRef<double> cx3_ref(&wk(iwk36), kdim-1, jdim-1);
                    FortranArray2DRef<double> dx3_ref(&wk(iwk37), kdim-1, jdim-1);
                    FortranArray2DRef<double> fx3_ref(&wk(iwk38), kdim-1, jdim-1);
                    FortranArray2DRef<double> by3_ref(&wk(iwk39), jdim-1, kdim-1);
                    FortranArray2DRef<double> cy3_ref(&wk(iwk40), jdim-1, kdim-1);
                    FortranArray2DRef<double> dy3_ref(&wk(iwk41), jdim-1, kdim-1);
                    FortranArray2DRef<double> fy3_ref(&wk(iwk42), jdim-1, kdim-1);
                    FortranArray2DRef<double> bz3_ref(&wk(iwk43), kdim-1, idim-1);
                    FortranArray2DRef<double> cz3_ref(&wk(iwk44), kdim-1, idim-1);
                    FortranArray2DRef<double> dz3_ref(&wk(iwk45), kdim-1, idim-1);
                    FortranArray2DRef<double> fz3_ref(&wk(iwk46), kdim-1, idim-1);
                    FortranArray2DRef<double> bx4_ref(&wk(iwk47), kdim-1, jdim-1);
                    FortranArray2DRef<double> cx4_ref(&wk(iwk48), kdim-1, jdim-1);
                    FortranArray2DRef<double> dx4_ref(&wk(iwk49), kdim-1, jdim-1);
                    FortranArray2DRef<double> fx4_ref(&wk(iwk50), kdim-1, jdim-1);
                    FortranArray2DRef<double> by4_ref(&wk(iwk51), jdim-1, kdim-1);
                    FortranArray2DRef<double> cy4_ref(&wk(iwk52), jdim-1, kdim-1);
                    FortranArray2DRef<double> dy4_ref(&wk(iwk53), jdim-1, kdim-1);
                    FortranArray2DRef<double> fy4_ref(&wk(iwk54), jdim-1, kdim-1);
                    FortranArray2DRef<double> bz4_ref(&wk(iwk55), kdim-1, idim-1);
                    FortranArray2DRef<double> cz4_ref(&wk(iwk56), kdim-1, idim-1);
                    FortranArray2DRef<double> dz4_ref(&wk(iwk57), kdim-1, idim-1);
                    FortranArray2DRef<double> fz4_ref(&wk(iwk58), kdim-1, idim-1);
                    int fdsav_sz = (ides != 0) ? (jdim-1)*(kdim-1)*(idim-1) : 1;
                    FortranArray3DRef<double> xlscale3d(&wk(iwk59), fdsav_sz, 1, 1);
                    FortranArray3DRef<double> fdsav3d(&wk(iwk59), fdsav_sz, 1, 1);
                    double sumn1_d = (double)sumn(1);
                    double sumn2_d = (double)sumn(2);
                    double sumn3_d = (double)sumn(3);
                    double sumn4_d = (double)sumn(4);
                    foureqn_ns::foureqn(jdim, kdim, idim, q, sj4d, sk4d, si4d,
                                        vol3d, dtj3d, x3d, y3d, z3d, vist3d,
                                        vor3d, smin3d, zksav4d, zksav4d,
                                        damp1_3d, blend3d, timestp3d, wrks3d, fnu3d,
                                        bx_ref, bx2_ref, cx_ref, cx2_ref,
                                        dx_ref, dx2_ref, fx_ref, fx2_ref, wx_ref,
                                        by_ref, by2_ref, cy_ref, cy2_ref,
                                        dy_ref, dy2_ref, fy_ref, fy2_ref, wy_ref,
                                        bz_ref, bz2_ref, cz_ref, cz2_ref,
                                        dz_ref, dz2_ref, fz_ref, fz2_ref, wz_ref,
                                        ntime, tj0, tk0, ti0, nbl,
                                        qj04d, qk04d, qi04d, vj0, vk0, vi0,
                                        blank, iover,
                                        sumn1_d, sumn2_d, sumn3_d, sumn4_d,
                                        negn(1), negn(2), negn(3), negn(4),
                                        ux, rhside4d, xib2, v3dtmp3d,
                                        bcj, bck, bci,
                                        nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                                        ibcinfo, jbcinfo, kbcinfo,
                                        maxbl, maxseg, volj0, volk0, voli0,
                                        nou, bou, nbuf, ibufdim, iex, iex2, iex3,
                                        bx3_ref, cx3_ref, dx3_ref, fx3_ref,
                                        by3_ref, cy3_ref, dy3_ref, fy3_ref,
                                        bz3_ref, cz3_ref, dz3_ref, fz3_ref,
                                        bx4_ref, cx4_ref, dx4_ref, fx4_ref,
                                        by4_ref, cy4_ref, dy4_ref, fy4_ref,
                                        bz4_ref, cz4_ref, dz4_ref, fz4_ref,
                                        xlscale3d, fdsav3d, nummem);
                    sumn(1) = sumn1_d; sumn(2) = sumn2_d;
                    sumn(3) = sumn3_d; sumn(4) = sumn4_d;
                }


            } else if (ivmx == 72) {
                // full stress model, k-omega based stress model
                if (isklton == 1) {
                    if (issglrrw2012 == 1) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " Computing turbulence using Eisfeld SSG/LRR-RSM-w2012 7-eqn RSM, block=%5d", nbl);
                    } else if (issglrrw2012 == 2) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " Computing turbulence using Eisfeld SSG/LRR-RSM-w2012 7-eqn RSM with F1=1, block=%5d", nbl);
                    } else if (issglrrw2012 == 3) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " Computing turbulence using Eisfeld SSG/LRR-RSM-w2012 7-eqn RSM with simple diffusion model, block=%5d", nbl);
                    } else if (issglrrw2012 == 4) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " Computing turbulence using Eisfeld SSG/LRR-RSM-w2012 7-eqn RSM with F1=1 and simple diffusion model, block=%5d", nbl);
                    } else if (issglrrw2012 == 5) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " Computing turbulence using Eisfeld SSG/LRR-RSM-w2012 7-eqn RSM with Wilcox simple diffusion (non-blended), block=%5d", nbl);
                    } else if (issglrrw2012 == 6) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " Computing turbulence using Eisfeld SSG/LRR-RSM-w2012-g 7-eqn RSM, block=%5d", nbl);
                    } else {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " Computing turbulence using Wilcox Stress-Omega 7-eqn RSM, block=%5d", nbl);
                    }
                    if (i_sas_rsm == 1) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "   ...adding SAS-like term (AIAA-2014-0586)");
                    } else if (i_sas_rsm == -1) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "   ...subtracting SAS-like term (AIAA-2014-0586)");
                    }
                    if (i_yapterm == 1) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "   ...adding Yap term");
                    }
                }
                if (icyc <= nfreeze) {
                    if (isklton > 0) {
                        nss = std::min(ncyc, nfreeze);
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " ");
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " RSM turbulence model is frozen for %5d iterations or subits", nss);
                    }
                    for (int n = 1; n <= nummem; n++) {
                        sumn(n) = 0.;
                        negn(n) = 0;
                    }
                } else {
                    FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                    FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> vj03d(volj0);
                    FortranArray3DRef<double> vk03d(volk0);
                    FortranArray3DRef<double> vi03d(voli0);
                    FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                    FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                    FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                    FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim1);
                    FortranArray3DRef<double> x3d(&x(1), jdim, kdim, idim);
                    FortranArray3DRef<double> y3d(&y(1), jdim, kdim, idim);
                    FortranArray3DRef<double> z3d(&z(1), jdim, kdim, idim);
                    FortranArray1DRef<int> negn_ref(&negn(1), nummem);
                    module_kwstm_ns::kws_main(jdim, kdim, idim, myid, myhost, nummem,
                                              tj0, tk0, ti0, xib,
                                              ux,
                                              sj4d, sk4d, si4d, vol3d,
                                              vj03d, vk03d, vi03d,
                                              q, qj04d, qk04d, qi04d,
                                              dtj3d, vist3d,
                                              level, icyc, sumn, negn_ref,
                                              xib2, snk0, x3d, y3d, z3d, nbl,
                                              issglrrw2012, i_sas_rsm, i_yapterm);
                }
            } else if (ivmx == 25) {
                // LES-type model and diagnostics
                FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
                FortranArray3DRef<double> vor3d(&res(1,1), jdim, kdim, idim1);
                FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                FortranArray3DRef<double> x3d(&x(1), jdim, kdim, idim);
                FortranArray3DRef<double> y3d(&y(1), jdim, kdim, idim);
                FortranArray3DRef<double> z3d(&z(1), jdim, kdim, idim);
                // xib passed as xjb (FortranArray3DRef) - reshape 4D to 3D
                FortranArray3DRef<double> xib3d(&xib(1,1,1,1), jdim*kdim*idim, nummem, 1);
                // xkb passed as tursav (FortranArray4DRef) - reshape 3D to 4D
                FortranArray4DRef<double> xkb4d(&xkb(1,1,1), jdim-1, kdim-1, idim-1, 1);
                lesdiag_ns::lesdiag(myid, jdim, kdim, idim, q, ux, vist3d,
                                    vol3d, si4d, sj4d, sk4d,
                                    vor3d, snk0, xib3d, xkb4d, xkb, blnum,
                                    nou, bou, nbuf, ibufdim, nbl, nummem,
                                    x3d, y3d, z3d);
            }
            // end if (icall > 0)
        } // end if (icall > 0)
    } // end if (ivisc > 1)
    } // end if (iadvance >= 0) for turbulence section


    //
    // residuals  res = r(q)
    //
    if (isklton > 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " computing residuals for block%6d", nbl);
        if (ifullns != 0 || ivmx >= 70) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Full N-S terms are ON - meanflow only");
        }
    }

    // end if (iadvance >= 0) for residuals message
    // (the message is inside the iadvance block in Fortran)

    iwk4 = 1;
    iwk5 = 1;
    iwk6 = 1;

    // zero residuals
    n    = jdim*kdim;
    nplq = std::min(idim1, 999000/n);
    npl  = nplq;
    for (i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        nnpl = n * npl;
        ist  = (i-1)*n + 1;
        for (l = 1; l <= 5; l++) {
            for (izz = 1; izz <= nnpl; izz++) {
                res(izz+ist-1, l) = 0.0;
            }
        }
    }

    //
    // imult: multi-plane vectorization flag
    //        = 0 single plane at a time
    //        > 0 multiple planes at a time
    //
    imult = 1;


    //
    // residuals   J direction
    //
    if (iadvance(nbl) >= 0) {
        nv = 35;
        if (idefrm(nbl) > 0) nv = 41;
        nvtq  = std::min(999000, nwork/nv);
        n     = jdim * kdim;
        nplq  = std::min(idim1, nvtq/n);
        if (imult == 0) nplq = 1;
        npl   = nplq;
        niwfac = 1;
        for (i = 1; i <= idim1; i += nplq) {
            if (i + npl - 1 > idim1) npl = idim1 - i + 1;
            nvtq  = npl * jdim * kdim;
            jbctyp[0] = jbcinfo(nbl, 1, 1, 1);
            jbctyp[1] = jbcinfo(nbl, 1, 1, 2);
            {
                FortranArray4DRef<double> res4d(&res(1,1), jdim*kdim, 1, idim1, 5);
                FortranArray3DRef<double> sj3d(&sj(1,1), jdim*kdim, idim1, 5);
                FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                FortranArray2DRef<double> wk2d(&wk(1), nvtq, nv);
                FortranArray1DRef<int> jbctyp_ref(jbctyp, 2);
                int idef_nbl = idefrm(nbl);
                double rkap2 = (double)rkap[1];
                gfluxr_ns::gfluxr(i, npl, rkap2, jdim, kdim, idim,
                                  res4d, q, qj04d, sj3d, wk2d, nvtq, nv,
                                  nfajki[0], wfa, FortranArray1DRef<int>(&iwfa(niwfac), 1),
                                  jbctyp_ref, isf, nbl, bcj,
                                  nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl, idef_nbl);
            }
        }
    }

    //
    // residuals   J direction - viscous terms
    //
    if (ivisc_i[1] > 0) {
        nvtq = std::min(999000, nwork/32);
        n    = jdim * kdim;
        nplq = std::min(idim1, nvtq/n);
        if (imult == 0) nplq = 1;
        npl = nplq;
        for (i = 1; i <= idim1; i += nplq) {
            if (i + npl - 1 > idim1) npl = idim1 - i + 1;
            nvtq = npl * jdim * kdim;
            {
                FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
                FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                FortranArray2DRef<double> wk2d(&wk(1), nvtq, 32);
                FortranArray2DRef<double> wj02d(&wj0(1), npl*(kdim-1), 22); // gfluxv.F declares wj0(npl*(kdim-1),22); jdim*kdim was wrong (over-sized for npl==1, but UNDER-sized -> column overflow when npl>1, e.g. ONERA_M6)
                int iadv_nbl = iadvance(nbl);
                gfluxv_ns::gfluxv(i, npl, jdim, kdim, idim,
                                  res4d, q, qj04d, sj4d, vol3d, wk2d, nvtq,
                                  wj02d, vist3d, vmuj, vj0, bcj, xib, tj0, cmuv,
                                  volj0, nou, bou, nbuf, ibufdim, iadv_nbl, nummem, ux);
            }
            //
            // full viscous terms
            //
            if (ifullns != 0) {
                FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
                FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
                FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                FortranArray2DRef<double> wk2d(&wk(1), nvtq, 32);
                FortranArray2DRef<double> wj02d(&wj0(1), npl*(kdim-1), 22); // gfluxv.F declares wj0(npl*(kdim-1),22); jdim*kdim was wrong (over-sized for npl==1, but UNDER-sized -> column overflow when npl>1, e.g. ONERA_M6)
                int iadv_nbl = iadvance(nbl);
                gfluxv1_ns::gfluxv1(i, npl, jdim, kdim, idim,
                                    res4d, q, qj04d, qk04d, qi04d,
                                    sj4d, sk4d, si4d, vol3d, wk2d, nvtq,
                                    wj02d, vist3d, vj0, bcj, bck, bci, volj0,
                                    nou, bou, nbuf, ibufdim, iadv_nbl);
            }
        }
    }

    if (icyc == 1 && iadvance(nbl) >= 0) {
        FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
        FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
        int ntime0 = 0; int irdq1 = 1;
        l2norm_ns::l2norm(nbl, ntime0, resd, irdq1, jdim, kdim, idim, res4d, vol3d);
    }


    //
    // residuals   K direction
    //
    if (iadvance(nbl) >= 0) {
        nv = 35;
        if (idefrm(nbl) > 0) nv = 41;
        nvtq = std::min(999000, nwork/nv);
        n    = jdim * kdim;
        nplq = std::min(idim1, nvtq/n);
        if (imult == 0) nplq = 1;
        npl  = nplq;
        if (nfajki[1] > 0) {
            niwfac = nfajki[0]*7 + 1;
        }
        for (i = 1; i <= idim1; i += nplq) {
            if (i + npl - 1 > idim1) npl = idim1 - i + 1;
            nvtq = npl * jdim * kdim;
            kbctyp[0] = kbcinfo(nbl, 1, 1, 1);
            kbctyp[1] = kbcinfo(nbl, 1, 1, 2);
            {
                FortranArray4DRef<double> res4d(&res(1,1), jdim*kdim, 1, idim1, 5);
                FortranArray3DRef<double> sk3d(&sk(1,1), jdim*kdim, idim1, 5);
                FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                FortranArray2DRef<double> wk2d(&wk(1), nvtq, nv);
                FortranArray1DRef<int> kbctyp_ref(kbctyp, 2);
                int idef_nbl = idefrm(nbl);
                double rkap3 = (double)rkap[2];
                hfluxr_ns::hfluxr(i, npl, rkap3, jdim, kdim, idim,
                                  res4d, q, qk04d, sk3d, wk2d, nvtq, nv,
                                  nfajki[1], wfa, FortranArray1DRef<int>(&iwfa(niwfac), 1),
                                  kbctyp_ref, isf, nbl, bck,
                                  nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl, idef_nbl);
            }
        }
    }

    //
    // residuals   K direction - viscous terms
    //
    if (ivisc_i[2] > 0) {
        nvtq = std::min(999000, nwork/32);
        n    = jdim * kdim;
        nplq = std::min(idim1, nvtq/n);
        if (imult == 0) nplq = 1;
        npl = nplq;
        for (i = 1; i <= idim1; i += nplq) {
            if (i + npl - 1 > idim1) npl = idim1 - i + 1;
            nvtq = npl * jdim * kdim;
            {
                FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
                FortranArray3DRef<double> sk3d(&sk(1,1), jdim*kdim, idim1, 5);
                FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                FortranArray2DRef<double> wk2d(&wk(1), nvtq, 32);
                FortranArray2DRef<double> wk02d(&wk0(1), npl*(jdim-1), 22); // hfluxv.F: wk0(npl*(jdim-1),22)
                int idf0 = 0;
                int iadv_nbl = iadvance(nbl);
                // ak, bk, ck are wk(iwk4), wk(iwk5), wk(iwk6) = wk(1), wk(1), wk(1)
                FortranArray4DRef<double> ak4d(&wk(iwk4), 1, 1, 1, 1);
                FortranArray4DRef<double> bk4d(&wk(iwk5), 1, 1, 1, 1);
                FortranArray4DRef<double> ck4d(&wk(iwk6), 1, 1, 1, 1);
                hfluxv_ns::hfluxv(i, npl, jdim, kdim, idim, idf0,
                                  ak4d, bk4d, ck4d,
                                  res4d, q, qk04d, sk3d, vol3d, wk2d, nvtq,
                                  wk02d, vist3d, vmuk, vk0, bck, xib, tk0, cmuv,
                                  volk0, nou, bou, nbuf, ibufdim, iadv_nbl, nummem, ux);
            }
            // full viscous terms (cross-derivatives)
            if (ifullns != 0) {
                FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
                FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                FortranArray3DRef<double> sk3d(&sk(1,1), jdim*kdim, idim1, 5);
                FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                FortranArray2DRef<double> wk2d(&wk(1), nvtq, 32);
                FortranArray2DRef<double> wk02d(&wk0(1), npl*(jdim-1), 22); // hfluxv.F: wk0(npl*(jdim-1),22)
                int iadv_nbl = iadvance(nbl);
                hfluxv1_ns::hfluxv1(i, npl, jdim, kdim, idim,
                                    res4d, q, qj04d, qk04d, qi04d,
                                    sj4d, sk3d, si4d, vol3d, wk2d, nvtq,
                                    wk02d, vist3d, vk0, bcj, bck, bci, volk0,
                                    nou, bou, nbuf, ibufdim, iadv_nbl);
            }
        }
    }

    if (icyc == 1 && iadvance(nbl) >= 0) {
        FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
        FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
        int ntime0 = 0; int irdq1 = 1;
        l2norm_ns::l2norm(nbl, ntime0, resd, irdq1, jdim, kdim, idim, res4d, vol3d);
    }


    //
    // residuals in  I direction
    //
    if (i2d == 0) {

        if (iadvance(nbl) >= 0) {
            n  = jdim * idim;
            nv = 35;
            if (idefrm(nbl) > 0) nv = 41;
            nvtq = std::min(999000, nwork/nv);
            nplq = std::min(kdim, nvtq/n);
            if (imult == 0) nplq = 1;
            npl  = nplq;
            if (nfajki[2] > 0) {
                niwfac = (nfajki[0] + nfajki[1])*7 + 1;
            }
            for (k = 1; k <= kdim1; k += nplq) {
                if (k + npl - 1 > kdim1) npl = kdim1 - k + 1;
                nvtq = nplq * jdim * idim;
                ibctyp[0] = ibcinfo(nbl, 1, 1, 1);
                ibctyp[1] = ibcinfo(nbl, 1, 1, 2);
                {
                    FortranArray4DRef<double> res4d(&res(1,1), jdim*kdim, 1, idim1, 5);
                    FortranArray3DRef<double> si3d(&si(1,1), jdim*kdim, idim, 5);
                    FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                    FortranArray2DRef<double> wk2d(&wk(1), nvtq, nv);
                    FortranArray1DRef<int> ibctyp_ref(ibctyp, 2);
                    int idef_nbl = idefrm(nbl);
                    double rkap1 = (double)rkap[0];
                    ffluxr_ns::ffluxr(k, npl, rkap1, jdim, kdim, idim,
                                      res4d, q, qi04d, si3d, wk2d, nvtq, nv,
                                      nfajki[2], wfa, FortranArray1DRef<int>(&iwfa(niwfac), 1),
                                      ibctyp_ref, isf, nbl, bci,
                                      nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl, idef_nbl);
                }
            }
        }

        //
        // residuals   I direction - viscous terms
        //
        if (ivisc_i[0] > 0) {
            nvtq = std::min(999000, nwork/32);
            n    = jdim * idim;
            nplq = std::min(kdim1, nvtq/n);
            if (imult == 0) nplq = 1;
            npl  = nplq;
            for (k = 1; k <= kdim1; k += nplq) {
                if (k + npl - 1 > kdim1) npl = kdim1 - k + 1;
                nvtq = npl * jdim * idim;
                {
                    FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                    FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                    FortranArray2DRef<double> wk2d(&wk(1), nvtq, 32);
                    FortranArray2DRef<double> wi02d(&wi0(1), npl*(jdim-1), 22); // ffluxv.F: wi0(npl*(jdim-1),22)
                    // vmui declared 3D at resid level (jdim-1,kdim-1,2); ffluxv's
                    // real definition takes it as 2D (jdim-1, (kdim-1)*2).
                    FortranArray2DRef<double> vmui2d(&vmui(1,1,1), jdim-1, (kdim-1)*2);
                    int iadv_nbl = iadvance(nbl);
                    ffluxv_ns::ffluxv(k, npl, jdim, kdim, idim,
                                      res4d, q, qi04d, si4d, vol3d, wk2d, nvtq,
                                      wi02d, vist3d, vmui2d, vi0, bci, xib, ti0, cmuv,
                                      voli0, nou, bou, nbuf, ibufdim, iadv_nbl, nummem, ux);
                }
                //
                // full viscous terms
                //
                if (ifullns != 0) {
                    FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> sj4d(&sj(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> sk4d(&sk(1,1), jdim, kdim, idim1, 5);
                    FortranArray4DRef<double> si4d(&si(1,1), jdim, kdim, idim, 5);
                    FortranArray4DRef<double> qj04d(&qj0(1,1,1), kdim, idim1, 5, 4);
                    FortranArray4DRef<double> qk04d(&qk0(1,1,1), jdim, idim1, 5, 4);
                    FortranArray4DRef<double> qi04d(&qi0(1,1,1), jdim, kdim, 5, 4);
                    FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                    FortranArray2DRef<double> wk2d(&wk(1), nvtq, 32);
                    FortranArray2DRef<double> wi02d(&wi0(1), npl*(jdim-1), 22); // ffluxv.F: wi0(npl*(jdim-1),22)
                    int iadv_nbl = iadvance(nbl);
                    ffluxv1_ns::ffluxv1(k, npl, jdim, kdim, idim,
                                        res4d, q, qj04d, qk04d, qi04d,
                                        sj4d, sk4d, si4d, vol3d, wk2d, nvtq,
                                        wi02d, vist3d, vi0, bcj, bck, bci, voli0,
                                        nou, bou, nbuf, ibufdim, iadv_nbl);
                }
            }
        }

        if (icyc == 1 && iadvance(nbl) >= 0) {
            FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
            FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
            int ntime0 = 0; int irdq1 = 1;
            l2norm_ns::l2norm(nbl, ntime0, resd, irdq1, jdim, kdim, idim, res4d, vol3d);
        }

    } // end if (i2d == 0)


    //
    // call resnonin to add rotating noninertial source term to res
    //
    if (noninflag > 0) {
        // resnonin expects 2D q (flat), 1D x/y/z/vol, 2D sj/sk/si/res
        FortranArray2DRef<double> q2d(&q(1,1,1,1), jdim*kdim*idim, 5);
        resnonin_ns::resnonin(nbl, jdim, kdim, idim, q2d, x, y, z, sj, sk, si, vol, res,
                              nou, bou, nbuf, ibufdim);
    }

    //
    // add source term in x-dir for problems with inflow/outflow bdys
    // that are periodic with each other
    //
    if (std::abs((float)xdir_only_source) > 1.e-12f) {
        if (isklton == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Adding x-dir source term to eqns");
        }
        izz = 0;
        for (i = 1; i <= idim-1; i++) {
            for (k = 1; k <= kdim; k++) {
                for (int j = 1; j <= jdim; j++) {
                    izz++;
                    fsrce = -(double)xdir_only_source * vol(izz);
                    res(izz, 2) = res(izz, 2) + fsrce;
                    res(izz, 5) = res(izz, 5) + fsrce * q(j, k, i, 2);
                }
            }
        }
    }

    if (iadvance(nbl) >= 0) {
        //
        // zero out residuals in edges of arrays (multi-plane vectorization)
        //
        for (i = 1; i <= idim1; i++) {
            jkv = (i-1) * jdim * kdim;
            for (l = 1; l <= 5; l++) {
                jj    = 1 - jdim;
                kdww2 = kdim - 2;
                for (ii = 1; ii <= kdww2; ii++) {
                    jj = jj + jdim;
                    res(jkv + jdim + jj - 1, l) = 0.;
                }
                for (izz = 1; izz <= jdim+1; izz++) {
                    res(izz + jkv + jdim*kdim1 - 1, l) = 0.;
                }
            }
        }

        if (iover == 1) {
            nvtq = std::min(999000, nwork/35);
            n    = jdim * kdim;
            nplq = std::min(idim1, nvtq/n);
            if (imult == 0) nplq = 1;
            npl  = nplq;

            for (i = 1; i <= idim1; i += nplq) {
                if (i + npl - 1 > idim1) npl = idim1 - i + 1;
                if (isklton > 0 && i == 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " blanking the rhs for the cells in the hole  block%5d", nbl);
                }
                {
                    FortranArray4DRef<double> res4d(&res(1,1), jdim*kdim, 1, idim1, 5);
                    hole_ns::hole(i, npl, jdim, kdim, idim, res4d, blank);
                }
            }
        }

    } // end if (iadvance >= 0)

    //
    // exact solution
    //
    if (iexact_trunc != 0 || iexact_disc != 0) {
        if (isklton > 0 && iexact_trunc != 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            if (iexact_trunc == 1) {
                std::snprintf(bou(nou(1),1), 120, "   Using exact MMS source in resid routine to evaluate truncation error, MS1");
            } else if (iexact_trunc == 2) {
                std::snprintf(bou(nou(1),1), 120, "   Using exact MMS source in resid routine to evaluate truncation error, MS2");
            } else if (iexact_trunc == 4) {
                std::snprintf(bou(nou(1),1), 120, "   Using exact MMS source in resid routine to evaluate truncation error, MS4");
            }
        } else if (isklton > 0 && iexact_disc != 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            if (iexact_disc == 1) {
                std::snprintf(bou(nou(1),1), 120, "   Using exact MMS source in resid routine to evaluate discretization error, MS1");
            } else if (iexact_disc == 2) {
                std::snprintf(bou(nou(1),1), 120, "   Using exact MMS source in resid routine to evaluate discretization error, MS2");
            } else if (iexact_disc == 4) {
                std::snprintf(bou(nou(1),1), 120, "   Using exact MMS source in resid routine to evaluate discretization error, MS4");
            }
        }
        {
            // dims from mms.f exact_flow_force:
            //   x,y,z(jdim,kdim,idim); vol(jdim,kdim,idim-1);
            //   res(jdim,kdim,idim-1,5); smin(jdim-1,kdim-1,idim-1)
            FortranArray3DRef<double> x3d(&x(1), jdim, kdim, idim);
            FortranArray3DRef<double> y3d(&y(1), jdim, kdim, idim);
            FortranArray3DRef<double> z3d(&z(1), jdim, kdim, idim);
            FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
            FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
            FortranArray3DRef<double> smin3d(&snk0(1,1,1), jdim-1, kdim-1, idim-1);
            mms_ns::exact_flow_force(jdim, kdim, idim, x3d, y3d, z3d, vol3d, res4d,
                                     smin3d, iexact_trunc, iexact_disc);
        }
        if (iexact_ring == 1) {
            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   Setting res=0 in outer ring");
            }
            {
                // zero_resid_ring(jdim,kdim,idim,res,jdim,kdim,idim-1,5,2,...)
                // res dimensioned res(jj,kk,ii,numeq) = (jdim,kdim,idim-1,5)
                int jj_ = jdim, kk_ = kdim, ii_ = idim-1;
                int numeq = 5, numrows = 2;
                FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
                mms_ns::zero_resid_ring(jdim, kdim, idim, res4d,
                                        jj_, kk_, ii_,
                                        numeq, numrows,
                                        iexact_trunc, iexact_disc);
            }
        }

        if (icyc == ncyc1[0] || icyc == ncyc1[1] ||
            icyc == ncyc1[2] || icyc == ncyc1[3] ||
            icyc == ncyc1[4]) {
            if (ntime == 1 && level >= lglobal) {
                {
                    // dims from mms.f exact_norm:
                    //   x,y,z,tursav(jdim,kdim,idim); q(jdim,kdim,idim,5);
                    //   vol(jdim,kdim,idim-1); res(jdim,kdim,idim-1,5);
                    //   smin(jdim-1,kdim-1,idim-1)
                    FortranArray3DRef<double> x3d(&x(1), jdim, kdim, idim);
                    FortranArray3DRef<double> y3d(&y(1), jdim, kdim, idim);
                    FortranArray3DRef<double> z3d(&z(1), jdim, kdim, idim);
                    FortranArray3DRef<double> tursav3d(&xib(1,1,1,1), jdim, kdim, idim);
                    FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim1);
                    FortranArray4DRef<double> res4d(&res(1,1), jdim, kdim, idim1, 5);
                    FortranArray3DRef<double> smin3d(&snk0(1,1,1), jdim-1, kdim-1, idim-1);
                    mms_ns::exact_norm(jdim, kdim, idim, x3d, y3d, z3d, q, tursav3d,
                                       vol3d, res4d, smin3d, iexact_trunc, iexact_disc);
                }
            }
        }
    }

    if (ipertavg == 1 || ipertavg == 2) {
        return;
    }

    if (iteravg == 1 || iteravg == 2) {
        if (level >= lglobal && ntime == nit) {
            if ((float)dt < 0.f || ((float)dt > 0.f && icyc == ncyc)) {
                // get iteration-averaged Q values and increment xnumavg
                for (i = 1; i <= idim-1; i++) {
                    for (int j = 1; j <= jdim-1; j++) {
                        for (k = 1; k <= kdim-1; k++) {
                            qavg(j,k,i,1) = (qavg(j,k,i,1)*(xnumavg-1.f) +
                                             q(j,k,i,1)) / xnumavg;
                            qavg(j,k,i,2) = (qavg(j,k,i,2)*(xnumavg-1.f) +
                                             q(j,k,i,1)*q(j,k,i,2)) / xnumavg;
                            qavg(j,k,i,3) = (qavg(j,k,i,3)*(xnumavg-1.f) +
                                             q(j,k,i,1)*q(j,k,i,3)) / xnumavg;
                            qavg(j,k,i,4) = (qavg(j,k,i,4)*(xnumavg-1.f) +
                                             q(j,k,i,1)*q(j,k,i,4)) / xnumavg;
                            qavg(j,k,i,5) = (qavg(j,k,i,5)*(xnumavg-1.f) +
                                             (q(j,k,i,5)/gm1 + 0.5*(q(j,k,i,2)*q(j,k,i,2) +
                                              q(j,k,i,3)*q(j,k,i,3) + q(j,k,i,4)*q(j,k,i,4))
                                              *q(j,k,i,1))) / xnumavg;
                        }
                    }
                }
            }
        }
    }
    return;
} // end resid

} // namespace resid_ns
