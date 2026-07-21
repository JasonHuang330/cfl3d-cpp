// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "module_kwstm.h"
#include "module_profileout.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <string>

// Forward declarations for termn8
namespace termn8_ns { void termn8(int& myid, int& icode, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou); }

namespace module_kwstm_ns {

using namespace module_profileout_ns;

// -----------------------------------------------------------------------
// Forward declarations for private module_kwstm helpers invoked by kws_main.
// These are (private) CONTAINS-subroutines of the Fortran MODULE module_kwstm;
// they are translated elsewhere but declared here so kws_main can call them.
// -----------------------------------------------------------------------
void kws_init_stress_from_2eq(int& nbl, int& jdim, int& kdim, int& idim, int& nummem,
                              FortranArray4DRef<double> ux, FortranArray3DRef<double> vist3d,
                              FortranArray3DRef<double> rho, FortranArray4DRef<double> tursav);
void fill_omega(int& jdim, int& kdim, int& idim,
                FortranArray4DRef<double> ux, FortranArray4DRef<double> omega);
void fill_fmu(int& jdim, int& kdim, int& idim,
              FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
              FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
              FortranArray3DRef<double> fmu);
void get_timestep(int& jdim, int& kdim, int& idim,
                  FortranArray3DRef<double> dtj, FortranArray3DRef<double> vol,
                  FortranArray4DRef<double> timestep, int& icyc, int& nummem);
void save_lasttimestep(int& jdim, int& kdim, int& idim, int& nummem,
                       FortranArray4DRef<double> qcur, FortranArray4DRef<double> qold,
                       FortranArray4DRef<double> dq);
void fill_turre(int& jdim, int& kdim, int& idim, int& nummem,
                FortranArray4DRef<double> turb, FortranArray4DRef<double> tj0,
                FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0,
                FortranArray4DRef<double> turre);
void fill_tke(int& jdim, int& kdim, int& idim, int& nummem,
              FortranArray4DRef<double> turre, FortranArray3DRef<double> tke,
              FortranArray4DRef<double> bij);
void fill_blend(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> blend,
                int& issglrrw2012, FortranArray3DRef<double> tke, FortranArray4DRef<double> turre,
                FortranArray3DRef<double> vol, FortranArray4DRef<double> si,
                FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
                FortranArray3DRef<double> smin, FortranArray4DRef<double> q,
                FortranArray3DRef<double> fmu, int& nummem);
void kws_dbij_dx(int& jdim, int& kdim, int& idim, int& nummem,
                 FortranArray4DRef<double> turre, FortranArray3DRef<double> tke,
                 FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
                 FortranArray4DRef<double> si, FortranArray3DRef<double> vol,
                 FortranArray5DRef<double> dbijdx);
void kws_dbij_dxx(int& jdim, int& kdim, int& idim, int& nummem,
                  FortranArray5DRef<double> dbijdx, FortranArray4DRef<double> sj,
                  FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
                  FortranArray3DRef<double> vol, FortranArray5DRef<double> dbijdxx);
void get_source(int& jdim, int& kdim, int& idim, int& nummem,
                FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
                FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
                FortranArray4DRef<double> turre, FortranArray3DRef<double> tke,
                FortranArray3DRef<double> blend, FortranArray4DRef<double> bb,
                FortranArray4DRef<double> omega, FortranArray4DRef<double> sj,
                FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
                FortranArray3DRef<double> vol, FortranArray3DRef<double> vj0,
                FortranArray3DRef<double> vk0, FortranArray3DRef<double> vi0,
                FortranArray4DRef<double> ux, FortranArray3DRef<double> fmu,
                FortranArray4DRef<double> source, FortranArray4DRef<double> rhs,
                FortranArray5DRef<double> d, FortranArray4DRef<double> zksav,
                FortranArray4DRef<double> timestep, FortranArray5DRef<double> dbijdxx,
                FortranArray3DRef<double> smin, int& issglrrw2012,
                FortranArray3DRef<double> x, FortranArray3DRef<double> y,
                FortranArray3DRef<double> z, int& nbl, int& icyc,
                FortranArray1DRef<int> ncyc1, int& i_sas_rsm, int& i_yapterm);
void get_diff_gen(int& jdim, int& kdim, int& idim, int& nummem,
                  FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
                  FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
                  FortranArray4DRef<double> turre, FortranArray3DRef<double> tke,
                  FortranArray3DRef<double> blend, FortranArray4DRef<double> sj,
                  FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
                  FortranArray3DRef<double> vol, FortranArray3DRef<double> vj0,
                  FortranArray3DRef<double> vk0, FortranArray3DRef<double> vi0,
                  FortranArray3DRef<double> fmu, FortranArray4DRef<double> rhs,
                  FortranArray5DRef<double> d, FortranArray4DRef<double> al,
                  FortranArray4DRef<double> ar, FortranArray4DRef<double> bl,
                  FortranArray4DRef<double> br, int& issglrrw2012);
void get_diffusion(int& jdim, int& kdim, int& idim, int& nummem,
                   FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
                   FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
                   FortranArray4DRef<double> turre, FortranArray3DRef<double> tke,
                   FortranArray3DRef<double> blend, FortranArray4DRef<double> sj,
                   FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
                   FortranArray3DRef<double> vol, FortranArray3DRef<double> vj0,
                   FortranArray3DRef<double> vk0, FortranArray3DRef<double> vi0,
                   FortranArray3DRef<double> fmu, FortranArray4DRef<double> rhs,
                   FortranArray5DRef<double> d, FortranArray4DRef<double> al,
                   FortranArray4DRef<double> ar, FortranArray4DRef<double> bl,
                   FortranArray4DRef<double> br, int& issglrrw2012);
void get_advection(int& jdim, int& kdim, int& idim, int& nummem,
                   FortranArray4DRef<double> q, FortranArray4DRef<double> turre,
                   FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
                   FortranArray4DRef<double> si, FortranArray3DRef<double> vol,
                   FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
                   FortranArray4DRef<double> qi0, FortranArray4DRef<double> rhs,
                   FortranArray5DRef<double> d, FortranArray4DRef<double> al,
                   FortranArray4DRef<double> ar, FortranArray4DRef<double> bl,
                   FortranArray4DRef<double> br);
void afsolver_k(int& jdim, int& kdim, int& idim, int& nummem,
                FortranArray4DRef<double> turre, FortranArray4DRef<double> q,
                FortranArray4DRef<double> qk0, FortranArray3DRef<double> fmu,
                FortranArray3DRef<double> tke, FortranArray4DRef<double> sk,
                FortranArray3DRef<double> vol, FortranArray3DRef<double> vk0,
                FortranArray4DRef<double> timestep, FortranArray4DRef<double> source,
                FortranArray4DRef<double> rhs, int& issglrrw2012);
void afsolver_j(int& jdim, int& kdim, int& idim, int& nummem,
                FortranArray4DRef<double> turre, FortranArray4DRef<double> q,
                FortranArray4DRef<double> qj0, FortranArray3DRef<double> fmu,
                FortranArray3DRef<double> tke, FortranArray4DRef<double> sj,
                FortranArray3DRef<double> vol, FortranArray3DRef<double> vj0,
                FortranArray4DRef<double> timestep, FortranArray4DRef<double> source,
                FortranArray4DRef<double> rhs, int& issglrrw2012);
void afsolver_i(int& jdim, int& kdim, int& idim, int& nummem,
                FortranArray4DRef<double> turre, FortranArray4DRef<double> q,
                FortranArray4DRef<double> qi0, FortranArray3DRef<double> fmu,
                FortranArray3DRef<double> tke, FortranArray4DRef<double> si,
                FortranArray3DRef<double> vol, FortranArray3DRef<double> vi0,
                FortranArray4DRef<double> timestep, FortranArray4DRef<double> source,
                FortranArray4DRef<double> rhs, int& issglrrw2012);
void update(int& jdim, int& kdim, int& idim, int& nummem,
            FortranArray4DRef<double> rhs, FortranArray4DRef<double> turre,
            FortranArray4DRef<double> rho, FortranArray3DRef<double> fmu,
            FortranArray4DRef<double> turb, FortranArray3DRef<double> vist3d,
            FortranArray1DRef<double> sumn, FortranArray1DRef<int> negn,
            FortranArray4DRef<double> ux, int& issglrrw2012);
void sgs_solver_2d(int& jdim, int& kdim, int& idim, int& nummem,
                   FortranArray4DRef<double> timestep, FortranArray5DRef<double> d,
                   FortranArray4DRef<double> al, FortranArray4DRef<double> ar,
                   FortranArray4DRef<double> bl, FortranArray4DRef<double> br,
                   FortranArray4DRef<double> rhs);

// -----------------------------------------------------------------------
// RSM (Reynolds-Stress-Model, k-omega based) helper definitions.
//
// These 19 CONTAINS-subroutines of the Fortran MODULE module_kwstm were
// dropped by the giant-file translation. They are ONLY reachable through
// kws_main / kws_init / kws_plot, all of which are gated behind ivmx==72
// (the full 7-eqn stress-transport turbulence model). The flat-plate
// validation case (grdflat5, ivisc=0,0,7 -> ivmx=7, k-omega SST 2-eqn via
// twoeqn) never enters this path, so these bodies are not exercised.
//
// They are provided as linkable definitions that ABORT LOUDLY if ever
// invoked, so: (a) the executable links and non-RSM cases run, and (b) any
// future ivmx==72 (RSM) case fails visibly here instead of silently
// producing wrong numbers. When an RSM case is needed, translate the real
// bodies from module_kwstm.f90 (get_source:671, get_diffusion:1300,
// afsolver_{i,j,k}:1751/1909/2064, get_advection:2690, update:2944,
// get_diff_gen:3550, etc.).
// -----------------------------------------------------------------------
[[noreturn]] static void kws_rsm_stub(const char* fn) {
    std::fprintf(stderr,
        "FATAL: module_kwstm_ns::%s is an untranslated RSM (ivmx==72) helper.\n"
        "This turbulence path was not translated; only the k-omega SST / non-RSM\n"
        "paths are supported. Translate the real body from module_kwstm.f90.\n",
        fn);
    std::abort();
}

void kws_init_stress_from_2eq(int& nbl, int& jdim, int& kdim, int& idim, int& nummem,
                              FortranArray4DRef<double> ux, FortranArray3DRef<double> vist3d,
                              FortranArray3DRef<double> rho, FortranArray4DRef<double> tursav)
{ (void)nbl;(void)jdim;(void)kdim;(void)idim;(void)nummem;(void)ux;(void)vist3d;(void)rho;(void)tursav;
  kws_rsm_stub("kws_init_stress_from_2eq"); }

void fill_omega(int& jdim, int& kdim, int& idim,
                FortranArray4DRef<double> ux, FortranArray4DRef<double> omega)
{ (void)jdim;(void)kdim;(void)idim;(void)ux;(void)omega; kws_rsm_stub("fill_omega"); }

void fill_fmu(int& jdim, int& kdim, int& idim,
              FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
              FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
              FortranArray3DRef<double> fmu)
{ (void)jdim;(void)kdim;(void)idim;(void)q;(void)qj0;(void)qk0;(void)qi0;(void)fmu; kws_rsm_stub("fill_fmu"); }

void get_timestep(int& jdim, int& kdim, int& idim,
                  FortranArray3DRef<double> dtj, FortranArray3DRef<double> vol,
                  FortranArray4DRef<double> timestep, int& icyc, int& nummem)
{ (void)jdim;(void)kdim;(void)idim;(void)dtj;(void)vol;(void)timestep;(void)icyc;(void)nummem; kws_rsm_stub("get_timestep"); }

void save_lasttimestep(int& jdim, int& kdim, int& idim, int& nummem,
                       FortranArray4DRef<double> qcur, FortranArray4DRef<double> qold,
                       FortranArray4DRef<double> dq)
{ (void)jdim;(void)kdim;(void)idim;(void)nummem;(void)qcur;(void)qold;(void)dq; kws_rsm_stub("save_lasttimestep"); }

void fill_turre(int& jdim, int& kdim, int& idim, int& nummem,
                FortranArray4DRef<double> turb, FortranArray4DRef<double> tj0,
                FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0,
                FortranArray4DRef<double> turre)
{ (void)jdim;(void)kdim;(void)idim;(void)nummem;(void)turb;(void)tj0;(void)tk0;(void)ti0;(void)turre; kws_rsm_stub("fill_turre"); }

void fill_tke(int& jdim, int& kdim, int& idim, int& nummem,
              FortranArray4DRef<double> turre, FortranArray3DRef<double> tke,
              FortranArray4DRef<double> bij)
{ (void)jdim;(void)kdim;(void)idim;(void)nummem;(void)turre;(void)tke;(void)bij; kws_rsm_stub("fill_tke"); }

void fill_blend(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> blend,
                int& issglrrw2012, FortranArray3DRef<double> tke, FortranArray4DRef<double> turre,
                FortranArray3DRef<double> vol, FortranArray4DRef<double> si,
                FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
                FortranArray3DRef<double> smin, FortranArray4DRef<double> q,
                FortranArray3DRef<double> fmu, int& nummem)
{ (void)jdim;(void)kdim;(void)idim;(void)blend;(void)issglrrw2012;(void)tke;(void)turre;(void)vol;(void)si;(void)sj;(void)sk;(void)smin;(void)q;(void)fmu;(void)nummem; kws_rsm_stub("fill_blend"); }

void kws_dbij_dx(int& jdim, int& kdim, int& idim, int& nummem,
                 FortranArray4DRef<double> turre, FortranArray3DRef<double> tke,
                 FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
                 FortranArray4DRef<double> si, FortranArray3DRef<double> vol,
                 FortranArray5DRef<double> dbijdx)
{ (void)jdim;(void)kdim;(void)idim;(void)nummem;(void)turre;(void)tke;(void)sj;(void)sk;(void)si;(void)vol;(void)dbijdx; kws_rsm_stub("kws_dbij_dx"); }

void kws_dbij_dxx(int& jdim, int& kdim, int& idim, int& nummem,
                  FortranArray5DRef<double> dbijdx, FortranArray4DRef<double> sj,
                  FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
                  FortranArray3DRef<double> vol, FortranArray5DRef<double> dbijdxx)
{ (void)jdim;(void)kdim;(void)idim;(void)nummem;(void)dbijdx;(void)sj;(void)sk;(void)si;(void)vol;(void)dbijdxx; kws_rsm_stub("kws_dbij_dxx"); }

void get_source(int& jdim, int& kdim, int& idim, int& nummem,
                FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
                FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
                FortranArray4DRef<double> turre, FortranArray3DRef<double> tke,
                FortranArray3DRef<double> blend, FortranArray4DRef<double> bb,
                FortranArray4DRef<double> omega, FortranArray4DRef<double> sj,
                FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
                FortranArray3DRef<double> vol, FortranArray3DRef<double> vj0,
                FortranArray3DRef<double> vk0, FortranArray3DRef<double> vi0,
                FortranArray4DRef<double> ux, FortranArray3DRef<double> fmu,
                FortranArray4DRef<double> source, FortranArray4DRef<double> rhs,
                FortranArray5DRef<double> d, FortranArray4DRef<double> zksav,
                FortranArray4DRef<double> timestep, FortranArray5DRef<double> dbijdxx,
                FortranArray3DRef<double> smin, int& issglrrw2012,
                FortranArray3DRef<double> x, FortranArray3DRef<double> y,
                FortranArray3DRef<double> z, int& nbl, int& icyc,
                FortranArray1DRef<int> ncyc1, int& i_sas_rsm, int& i_yapterm)
{ kws_rsm_stub("get_source"); }

void get_diff_gen(int& jdim, int& kdim, int& idim, int& nummem,
                  FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
                  FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
                  FortranArray4DRef<double> turre, FortranArray3DRef<double> tke,
                  FortranArray3DRef<double> blend, FortranArray4DRef<double> sj,
                  FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
                  FortranArray3DRef<double> vol, FortranArray3DRef<double> vj0,
                  FortranArray3DRef<double> vk0, FortranArray3DRef<double> vi0,
                  FortranArray3DRef<double> fmu, FortranArray4DRef<double> rhs,
                  FortranArray5DRef<double> d, FortranArray4DRef<double> al,
                  FortranArray4DRef<double> ar, FortranArray4DRef<double> bl,
                  FortranArray4DRef<double> br, int& issglrrw2012)
{ kws_rsm_stub("get_diff_gen"); }

void get_diffusion(int& jdim, int& kdim, int& idim, int& nummem,
                   FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
                   FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
                   FortranArray4DRef<double> turre, FortranArray3DRef<double> tke,
                   FortranArray3DRef<double> blend, FortranArray4DRef<double> sj,
                   FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
                   FortranArray3DRef<double> vol, FortranArray3DRef<double> vj0,
                   FortranArray3DRef<double> vk0, FortranArray3DRef<double> vi0,
                   FortranArray3DRef<double> fmu, FortranArray4DRef<double> rhs,
                   FortranArray5DRef<double> d, FortranArray4DRef<double> al,
                   FortranArray4DRef<double> ar, FortranArray4DRef<double> bl,
                   FortranArray4DRef<double> br, int& issglrrw2012)
{ kws_rsm_stub("get_diffusion"); }

void get_advection(int& jdim, int& kdim, int& idim, int& nummem,
                   FortranArray4DRef<double> q, FortranArray4DRef<double> turre,
                   FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
                   FortranArray4DRef<double> si, FortranArray3DRef<double> vol,
                   FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
                   FortranArray4DRef<double> qi0, FortranArray4DRef<double> rhs,
                   FortranArray5DRef<double> d, FortranArray4DRef<double> al,
                   FortranArray4DRef<double> ar, FortranArray4DRef<double> bl,
                   FortranArray4DRef<double> br)
{ kws_rsm_stub("get_advection"); }

void afsolver_k(int& jdim, int& kdim, int& idim, int& nummem,
                FortranArray4DRef<double> turre, FortranArray4DRef<double> q,
                FortranArray4DRef<double> qk0, FortranArray3DRef<double> fmu,
                FortranArray3DRef<double> tke, FortranArray4DRef<double> sk,
                FortranArray3DRef<double> vol, FortranArray3DRef<double> vk0,
                FortranArray4DRef<double> timestep, FortranArray4DRef<double> source,
                FortranArray4DRef<double> rhs, int& issglrrw2012)
{ kws_rsm_stub("afsolver_k"); }

void afsolver_j(int& jdim, int& kdim, int& idim, int& nummem,
                FortranArray4DRef<double> turre, FortranArray4DRef<double> q,
                FortranArray4DRef<double> qj0, FortranArray3DRef<double> fmu,
                FortranArray3DRef<double> tke, FortranArray4DRef<double> sj,
                FortranArray3DRef<double> vol, FortranArray3DRef<double> vj0,
                FortranArray4DRef<double> timestep, FortranArray4DRef<double> source,
                FortranArray4DRef<double> rhs, int& issglrrw2012)
{ kws_rsm_stub("afsolver_j"); }

void afsolver_i(int& jdim, int& kdim, int& idim, int& nummem,
                FortranArray4DRef<double> turre, FortranArray4DRef<double> q,
                FortranArray4DRef<double> qi0, FortranArray3DRef<double> fmu,
                FortranArray3DRef<double> tke, FortranArray4DRef<double> si,
                FortranArray3DRef<double> vol, FortranArray3DRef<double> vi0,
                FortranArray4DRef<double> timestep, FortranArray4DRef<double> source,
                FortranArray4DRef<double> rhs, int& issglrrw2012)
{ kws_rsm_stub("afsolver_i"); }

void update(int& jdim, int& kdim, int& idim, int& nummem,
            FortranArray4DRef<double> rhs, FortranArray4DRef<double> turre,
            FortranArray4DRef<double> rho, FortranArray3DRef<double> fmu,
            FortranArray4DRef<double> turb, FortranArray3DRef<double> vist3d,
            FortranArray1DRef<double> sumn, FortranArray1DRef<int> negn,
            FortranArray4DRef<double> ux, int& issglrrw2012)
{ kws_rsm_stub("update"); }

void sgs_solver_2d(int& jdim, int& kdim, int& idim, int& nummem,
                   FortranArray4DRef<double> timestep, FortranArray5DRef<double> d,
                   FortranArray4DRef<double> al, FortranArray4DRef<double> ar,
                   FortranArray4DRef<double> bl, FortranArray4DRef<double> br,
                   FortranArray4DRef<double> rhs)
{ kws_rsm_stub("sgs_solver_2d"); }

void kws_dump_movie(int& nframe, int& jdim, int& kdim, int& idim, int& nummem, int& nblk,
                    FortranArray4DRef<double> q, FortranArray3DRef<double> x,
                    FortranArray3DRef<double> y, FortranArray3DRef<double> z,
                    FortranArray3DRef<double> vist3d, FortranArray4DRef<double> ux,
                    FortranArray4DRef<double> turre, FortranArray4DRef<double> tursav,
                    FortranArray3DRef<double> smin)
{ kws_rsm_stub("kws_dump_movie"); }

// -----------------------------------------------------------------------
// kws_get_nummem
// -----------------------------------------------------------------------
double kws_get_nummem() {
    return 7;
}

// -----------------------------------------------------------------------
// kws_initvist
// -----------------------------------------------------------------------
void kws_initvist(int& jdim, int& kdim, int& idim, int& nummem,
                  FortranArray4DRef<double> cij, FortranArray3DRef<double> vist3d) {
    float* tur10 = &cmn_ivals.tur10[0]; // 0-based C array, tur10[m-1] = Fortran tur10(m)
    for (int i = 1; i <= idim-1; i++) {
        for (int k = 1; k <= kdim-1; k++) {
            for (int j = 1; j <= jdim-1; j++) {
                for (int m = 1; m <= nummem; m++) {
                    cij(j,k,i,m) = (double)tur10[m-1];
                }
                vist3d(j,k,i) = 0.01;
            }
        }
    }
}

// -----------------------------------------------------------------------
// kws_init
// -----------------------------------------------------------------------
void kws_init(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
              int& nbuf, int& ibufdim) {
    static int ivisited = 0;
    if (ivisited == 0) {
        kws_read_parm(nou, bou, nbuf, ibufdim);
        ivisited = 1;
    }
}

// -----------------------------------------------------------------------
// toupper helper (nested in kws_read_parm in Fortran)
// -----------------------------------------------------------------------
static void toupper_local(char* str, int len) {
    for (int i = 0; i < len; i++) {
        if (str[i] >= 'a' && str[i] <= 'z')
            str[i] = (char)(str[i] - 'a' + 'A');
    }
}

// -----------------------------------------------------------------------
// kws_read_parm
// -----------------------------------------------------------------------
void kws_read_parm(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
                   int& nbuf, int& ibufdim) {
    int32_t& myid_c   = cmn_mydist2.myid;
    int32_t& myhost_c = cmn_mydist2.myhost;
    int32_t& nsubturb = cmn_turbconv.nsubturb;

    auto wr_parm = [](const char* key, float rval) {
        FILE* f = fortran_get_unit(120);
        if (f) fprintf(f, "%-20s %12.5E\n", key, (double)rval);
    };

    bool lexist = fortran_file_exists("kwstm.dat");
    nsubit = (int)nsubturb;

    if (lexist) {
        fortran_open_unit(1234, "kwstm.dat", "r");
        FILE* fp = fortran_get_unit(1234);
        char cline[71];
        while (true) {
            if (fgets(cline, sizeof(cline), fp) == nullptr) break;
            int clen = (int)strlen(cline);
            if (clen > 0 && cline[clen-1] == '\n') { cline[--clen] = '\0'; }
            // skip blank lines
            bool allblank = true;
            for (int ii = 0; ii < clen; ii++) if (cline[ii] != ' ') { allblank = false; break; }
            if (allblank) continue;
            // adjustl: find first non-space
            int first = 0;
            while (first < clen && cline[first] == ' ') first++;
            if (cline[first] == '!') continue;

            char key[71] = {0};
            float rval = 0.0f;
            int nread = sscanf(cline, "%70s %f", key, &rval);
            int iost1 = (nread < 2) ? -1 : 0;

            toupper_local(key, (int)strlen(key));

            if (iost1 != 0) {
                if (myid_c == myhost_c || myid_c == 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " Unknown Input line in kwstm.dat:");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "%.70s", cline);
                }
            }

            std::string skey(key);
            if      (skey == "C1")           { c1 = rval; }
            else if (skey == "C2")           { c2 = rval; }
            else if (skey == "ALPHA")        { alpha = rval; }
            else if (skey == "BETA_0")       { beta_0 = rval; }
            else if (skey == "BETA_STAR")    { beta_star = rval; }
            else if (skey == "SIGMA")        { sigma = rval; }
            else if (skey == "SIGMA_STAR")   { sigma_star = rval; }
            else if (skey == "SIGMA_D0")     { sigma_d0 = rval; }
            else if (skey == "JCUT")         { jcut = (int)std::round(rval); }
            else if (skey == "JSTART")       { jstart = (int)std::round(rval); }
            else if (skey == "KDIFF")        { kdiff = (int)std::round(rval); }
            else if (skey == "IREAL")        { irealizability = (int)std::round(rval); }
            else if (skey == "IOPT")         { iopt = (int)std::round(rval); }
            else if (skey == "CFLLOC")       { cfl_loc = rval; }
            else if (skey == "PDRATIO")      { pdratio = rval; }
            else if (skey == "MUTRATIO")     { xmutratio = rval; }
            else if (skey == "ILDTS")        { ildts = rval; }
            else if (skey == "FACTOR")       { factor = rval; }
            else if (skey == "NSUBIT")       { nsubit = (int)std::round(rval); }
            else if (skey == "NFREQ")        { nfreq = (int)std::round(rval); }
            else if (skey == "PROFILE")      { /* no-op */ }
            else if (skey == "CP_PROFILE")   { /* no-op */ }
            else if (skey == "CF_YPLUS_PLOT"){ need_kwsplot = ((int)std::round(rval) == 1); }
            else if (skey == "CONTOUR_PLOT") { need_contplot = ((int)std::round(rval) == 1); }
            else {
                if (myid_c == myhost_c || myid_c == 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " ERROR: Unknown Input keyword in kwstm.dat:");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    int fi = 0;
                    while (fi < clen && cline[fi] == ' ') fi++;
                    std::snprintf(bou(nou(1),1), 120, "%.70s", cline+fi);
                    int myid_ref = (int)myid_c;
                    int icode = -1;
                    termn8_ns::termn8(myid_ref, icode, ibufdim, nbuf, bou, nou);
                }
            }
        }
        // label 100
        fortran_close_unit(1234);
        module_profileout_ns::profile_read_input(*(char(*)[1])"kwstm.dat");
    }

    alpha_hat = (8.0 + c2) / 11.0;
    beta_hat  = (8.0*c2 - 2.0) / 11.0;
    gamma_hat = (60.0*c2 - 4.0) / 55.0;

    if (myid_c == myhost_c || myid_c == 1) {
        fortran_open_unit(120, "kwstm.out", "w");
        FILE* f120 = fortran_get_unit(120);
        fprintf(f120, " -----------------------------\n");
        fprintf(f120, " List of Wilcox Stress-Omega RSM model Constants\n");
        fprintf(f120, " -----------------------------\n");
        wr_parm("c1",         (float)c1);
        wr_parm("c2",         (float)c2);
        wr_parm("alpha",      (float)alpha);
        wr_parm("beta_0",     (float)beta_0);
        wr_parm("beta_star",  (float)beta_star);
        wr_parm("sigma",      (float)sigma);
        wr_parm("sigma_star", (float)sigma_star);
        wr_parm("sigma_d0",   (float)sigma_d0);
        fprintf(f120, " -----------------------------\n");
        fprintf(f120, " B.L. Profile Location\n");
        fprintf(f120, " -----------------------------\n");
        wr_parm("jcut",   (float)jcut);
        wr_parm("jstart", (float)jstart);
        fprintf(f120, " -----------------------------\n");
        fprintf(f120, " K diffusion option 0-T_ij, 1-Robs K version(not supported), 2-(6.50)  3-(6.49)\n");
        fprintf(f120, " -----------------------------\n");
        wr_parm("kdiff", (float)kdiff);
        fprintf(f120, " -----------------------------\n");
        fprintf(f120, " Realizability  0-off  1-on\n");
        fprintf(f120, " -----------------------------\n");
        wr_parm("ireal", (float)irealizability);
        fprintf(f120, " -----------------------------\n");
        fprintf(f120, " Miscellaneous\n");
        fprintf(f120, " -----------------------------\n");
        wr_parm("cflloc",   (float)cfl_loc);
        wr_parm("pdratio",  (float)pdratio);
        wr_parm("mutratio", (float)xmutratio);
        wr_parm("factor",   (float)factor);
        wr_parm("nsubit",   (float)nsubit);
        wr_parm("nfreq",    (float)nfreq);
        wr_parm("ildts",    (float)ildts);
        int ival = 0;
        if (need_contplot) ival = 1;
        wr_parm("contour_plot", (float)ival);
        ival = 0;
        if (need_kwsplot) ival = 1;
        wr_parm("CF_YPLUS_PLOT", (float)ival);
        fprintf(f120, " -----------------------------\n");
        fortran_close_unit(120);
    }
}



// -----------------------------------------------------------------------
// kws_plot
// -----------------------------------------------------------------------
void kws_plot(int& nframe, int& jdim, int& kdim, int& idim, int& nummem, int& nblk,
              FortranArray4DRef<double> q, FortranArray3DRef<double> x, FortranArray3DRef<double> y,
              FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d,
              FortranArray4DRef<double> ux, FortranArray4DRef<double> turre,
              FortranArray3DRef<double> smin) {
    if (!need_kwsplot) return;
    float& reue_c  = cmn_reyue.reue;
    float& xmach_c = cmn_info.xmach;
    float& gamma_c = cmn_fluid.gamma;
    float& cbar_c  = cmn_fluid2.cbar;
    float& tinf_c  = cmn_reyue.tinf;
    float& p0_c    = cmn_ivals.p0;
    float& rho0_c  = cmn_ivals.rho0;
    float& u0_c    = cmn_ivals.u0;
    int32_t& i2d_c = cmn_twod.i2d;
    char str[21]; std::snprintf(str, sizeof(str), "%04d", nblk);
    if (i2d_c == 1) {
        char fname1[64], fname2[64];
        std::snprintf(fname1, sizeof(fname1), "yplus-%s.plt", str);
        std::snprintf(fname2, sizeof(fname2), "Cf-ReL-%s.plt", str);
        fortran_open_unit(1235, fname1, "w");
        fortran_open_unit(1236, fname2, "w");
        FILE* f1236 = fortran_get_unit(1236);
        fprintf(f1236, "variables=rel,cf,cp,x\n");
        fprintf(f1236, "ZONE I= %4d\n", jdim-1);
        int i = 1;
        float re = reue_c / xmach_c;
        float c2b = cbar_c / tinf_c;
        float c2bp = c2b + 1.0f;
        std::vector<float> fnu(jdim+2,0.f), utau(jdim+2,0.f);
        float deltabl = 0.0f;
        for (int j = 1; j <= jdim-1; j++) {
            float tt = (float)(gamma_c * q(j,1,i,5) / q(j,1,i,1));
            fnu[j] = c2bp * tt * std::sqrt(tt) / (c2b + tt);
            float vel = (float)std::sqrt(q(j,1,i,2)*q(j,1,i,2)+q(j,1,i,3)*q(j,1,i,3)+q(j,1,i,4)*q(j,1,i,4));
            float tau_w = fnu[j] * std::abs(vel) / (float)std::abs(smin(j,1,i)) / re;
            utau[j] = std::sqrt(tau_w / (float)q(j,1,i,1));
            float Rel = (float)((0.5*(x(j,1,i)+x(j+1,1,i)) - x(jstart,1,i)) * reue_c);
            fprintf(f1236, " %12.5E %12.5E %12.5E %12.5E\n",
                    (double)Rel,
                    (double)(tau_w/(0.5f*rho0_c*(float)(xmach_c*xmach_c))),
                    (double)((q(j,1,i,5)-p0_c)/(0.5f*rho0_c*(float)(xmach_c*xmach_c))),
                    (double)(0.5*(x(j,1,i)+x(j+1,1,i))));
        }
        fortran_close_unit(1236);
        int j = jcut;
        for (int k = 1; k <= kdim-2; k++) {
            float vel  = (float)std::sqrt(q(j,k,i,2)*q(j,k,i,2)+q(j,k,i,3)*q(j,k,i,3)+q(j,k,i,4)*q(j,k,i,4));
            float vel1 = (float)std::sqrt(q(j,k+1,i,2)*q(j,k+1,i,2)+q(j,k+1,i,3)*q(j,k+1,i,3)+q(j,k+1,i,4)*q(j,k+1,i,4));
            if (vel < 0.99f*u0_c && vel1 > 0.99f*u0_c) {
                float coef = (float)((0.99*u0_c - q(j,k,i,2))/(q(j,k+1,i,2)-q(j,k,i,2)));
                deltabl = (float)(smin(j,k,i) + coef*(smin(j,k+1,i)-smin(j,k,i)));
            }
        }
        FILE* f1235 = fortran_get_unit(1235);
        fprintf(f1235, "title=\" utau= %12.5E fnuw= %12.5E Re= %12.5E Mach= %12.5E delta= %12.5E\"\n",
                (double)utau[j],(double)fnu[j],(double)reue_c,(double)xmach_c,(double)deltabl);
        if (nummem == 7) {
            fprintf(f1235, "variables=\"y+\",\"u+\",\"y\",\"rho\",\"u\",\"v\",\"w\",\"p\","
                    "\"t11\",\"t22\",\"t33\",\"t12\",\"t23\",\"t13\",\"zeta\",\"mut\",\"dist\"\n");
        } else {
            fprintf(f1235, "variables=\"y+\",\"u+\",\"y\",\"rho\",\"u\",\"v\",\"w\",\"p\","
                    "\"ome\",\"tke\",\"mut\",\"dist\"\n");
        }
        fprintf(f1235, "ZONE I= %4d\n", kdim-1);
        for (int k = 1; k <= kdim-1; k++) {
            float vel = (float)std::sqrt(q(j,k,i,2)*q(j,k,i,2)+q(j,k,i,3)*q(j,k,i,3)+q(j,k,i,4)*q(j,k,i,4));
            fprintf(f1235, " %12.5E %12.5E %12.5E",
                    (double)((float)q(j,1,i,1)*std::abs((float)smin(j,k,i))*utau[j]*re/fnu[j]),
                    (double)(vel/(utau[j]+1.e-20f)),
                    (double)(smin(j,k,i)/deltabl));
            for (int mm = 1; mm <= 5; mm++) fprintf(f1235, " %12.5E", (double)q(j,k,i,mm));
            for (int mm = 1; mm <= nummem; mm++) fprintf(f1235, " %12.5E", (double)turre(j,k,i,mm));
            fprintf(f1235, " %12.5E %12.5E\n", (double)vist3d(j,k,i), (double)smin(j,k,i));
        }
        fortran_close_unit(1235);
    }
}


// -----------------------------------------------------------------------
// kws_main
//   Driver for the k-omega SST / SSG-LRR-w2012 stress-transport turbulence
//   model. Feeds the turbulence source into the residual. Translated from
//   SUBROUTINE kws_main (module_kwstm.f90, lines 439-669).
// -----------------------------------------------------------------------
void kws_main(int& jdim, int& kdim, int& idim, int& myid, int& myhost, int& nummem,
              FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0,
              FortranArray4DRef<double> ti0, FortranArray4DRef<double> turb,
              FortranArray4DRef<double> ux,
              FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
              FortranArray4DRef<double> si, FortranArray3DRef<double> vol,
              FortranArray3DRef<double> vj0, FortranArray3DRef<double> vk0,
              FortranArray3DRef<double> vi0,
              FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
              FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
              FortranArray3DRef<double> dtj, FortranArray3DRef<double> vist3d,
              int& level, int& icyc, FortranArray1DRef<double> sumn,
              FortranArray1DRef<int> negn, FortranArray4DRef<double> zksav2,
              FortranArray3DRef<double> smin, FortranArray3DRef<double> x,
              FortranArray3DRef<double> y, FortranArray3DRef<double> z,
              int& nbl, int& issglrrw2012, int& i_sas_rsm, int& i_yapterm) {

    // common /info/ (only ntt and ncyc1 are used here)
    int32_t& ntt = cmn_info.ntt;
    FortranArray1DRef<int> ncyc1(cmn_info.ncyc1, 5);

    // local static working arrays (Fortran SAVE, ALLOCATABLE) — persist
    // across calls; allocated/deallocated each invocation as in the Fortran.
    static FortranArray4D<double> turre;   // 6 (uiuj) + 1 (zeta), with ghost cells
    static FortranArray3D<double> tke;     // turbulent kinetic energy
    static FortranArray3D<double> blend;   // Menter's F1 blending function
    static FortranArray4D<double> omega;   // vorticity vector
    static FortranArray4D<double> bij;
    static FortranArray3D<double> fmu;     // laminar dynamic viscosity (density units)
    static FortranArray4D<double> timestep;
    static FortranArray5D<double> d;       // source term Jacobian
    static FortranArray4D<double> al, ar, bl, br;
    static FortranArray5D<double> dbijdx, dbijdxx;

    int isub;
    static int icycle = 0;
    static int ntt_start = 0;
    int m;

    if (src_k.allocated()) {
        src_k.deallocate();
    }
    src_k.allocate(jdim, kdim, idim);
    if (icycle == 0 && lkzstm_from_2eq) {
        ntt_start = ntt;
    }
    icycle = icycle + 1;
    if (ntt == ntt_start && lkzstm_from_2eq) {
        // q(1,1,1,1) : first (jdim,kdim,idim) slab of q, viewed as rho(jdim,kdim,idim)
        FortranArray3DRef<double> rho_slab(q.data(), jdim, kdim, idim);
        kws_init_stress_from_2eq(nbl, jdim, kdim, idim, nummem, ux, vist3d, rho_slab, turb);
        int nf1 = 1;
        kws_dump_movie(nf1, jdim, kdim, idim, nummem, nbl, q, x, y, z,
                       vist3d, ux, turb, turb, smin);
        // let the cfl3d set up correct ghost cell values
        return;
    }

    mylevel = level;
    myicyc = icyc;

    myrank = myid;
    myroot = myhost;

    if (source_items.allocated()) {
        source_items.deallocate();
    }
    source_items.allocate(9, nummem, jdim, kdim, idim);

    // 6 stress components + 1 zeta, including ghost cells
    turre.allocate_range(-1, jdim + 1, -1, kdim + 1, -1, idim + 1, 1, nummem);
    d.allocate(nummem, nummem, jdim, kdim, idim);
    al.allocate(2, jdim, kdim, idim);
    ar.allocate(2, jdim, kdim, idim);
    bl.allocate(2, jdim, kdim, idim);
    br.allocate(2, jdim, kdim, idim);

    // turbulent kinetic energy
    tke.allocate_range(0, jdim, 0, kdim, 0, idim);

    // Menter's F1 blending function
    blend.allocate(jdim - 1, kdim - 1, idim - 1);

    // vorticity
    omega.allocate_range(0, jdim, 0, kdim, 0, idim, 1, 3);

    // source terms
    if (source.allocated()) {
        source.deallocate();
    }
    source.allocate(jdim - 1, kdim - 1, idim - 1, nummem);
    if (rhs.allocated()) {
        rhs.deallocate();
    }
    rhs.allocate(jdim - 1, kdim - 1, idim - 1, nummem);

    // bij terms
    bij.allocate_range(0, jdim, 0, kdim, 0, idim, 1, 6);

    // fmu array
    fmu.allocate_range(0, jdim, 0, kdim, 0, idim);

    timestep.allocate(jdim - 1, kdim - 1, idim - 1, nummem);

    // derivatives for c5 terms
    dbijdx.allocate(jdim - 1, kdim - 1, idim - 1, 6, 3);
    dbijdxx.allocate(jdim - 1, kdim - 1, idim - 1, 6, 3);

    // note: omega (here) is vorticity!
    fill_omega(jdim, kdim, idim, ux, omega);
    fill_fmu(jdim, kdim, idim, q, qj0, qk0, qi0, fmu);
    get_timestep(jdim, kdim, idim, dtj, vol, timestep, icyc, nummem);

    if (icyc == 1) {
        // zksav2(1,1,1,1) and zksav2(1,1,1,nummem+1): two (jdim,kdim,idim,nummem)
        // slabs of zksav2, dimensioned (jdim,kdim,idim,2*nummem)
        size_t slab = static_cast<size_t>(jdim) * kdim * idim;
        FortranArray4DRef<double> zk_cur(zksav2.data(), jdim, kdim, idim, nummem);
        FortranArray4DRef<double> zk_old(zksav2.data() + static_cast<size_t>(nummem) * slab,
                                         jdim, kdim, idim, nummem);
        save_lasttimestep(jdim, kdim, idim, nummem, turb, zk_cur, zk_old);
    }
    // calculate the source terms
    for (isub = 1; isub <= nsubit; isub++) {
        source = 0.0;
        fill_turre(jdim, kdim, idim, nummem, turb, tj0, tk0, ti0, turre);
        fill_tke(jdim, kdim, idim, nummem, turre, tke, bij);
        fill_blend(jdim, kdim, idim, blend, issglrrw2012, tke, turre, vol, si, sj,
                   sk, smin, q, fmu, nummem);
        kws_dbij_dx(jdim, kdim, idim, nummem, turre, tke, sj, sk, si, vol, dbijdx);
        kws_dbij_dxx(jdim, kdim, idim, nummem, dbijdx, sj, sk, si, vol, dbijdxx);
        get_source(jdim, kdim, idim, nummem, q, qj0, qk0, qi0, turre, tke, blend, bij, omega,
                   sj, sk, si, vol, vj0, vk0, vi0, ux, fmu, source, rhs, d, zksav2, timestep, dbijdxx,
                   smin, issglrrw2012, x, y, z, nbl, icyc, ncyc1, i_sas_rsm, i_yapterm);

        al = 0; ar = 0; bl = 0; br = 0;
        if (issglrrw2012 == 1 || issglrrw2012 == 2 || issglrrw2012 == 6) {
            // generalized gradient-diffusion
            get_diff_gen(jdim, kdim, idim, nummem, q, qj0, qk0, qi0, turre, tke, blend,
                         sj, sk, si, vol, vj0, vk0, vi0, fmu, rhs, d, al, ar, bl, br, issglrrw2012);
        } else {
            // standard-type diffusion
            get_diffusion(jdim, kdim, idim, nummem, q, qj0, qk0, qi0, turre, tke, blend,
                          sj, sk, si, vol, vj0, vk0, vi0, fmu, rhs, d, al, ar, bl, br, issglrrw2012);
        }

        get_advection(jdim, kdim, idim, nummem, q, turre, sj, sk, si, vol, qj0, qk0, qi0, rhs, d, al, ar, bl, br);
        // solve the (I - J dt )*\Delta Q = Source*Dt
        if (iopt < 2 || true) {
            for (m = 1; m <= nummem; m++) {
                for (int i = 1; i <= idim - 1; i++) {
                    for (int k = 1; k <= kdim - 1; k++) {
                        for (int j = 1; j <= jdim - 1; j++) {
                            rhs(j, k, i, m) = rhs(j, k, i, m) * timestep(j, k, i, m);
                        }
                    }
                }
            }
            afsolver_k(jdim, kdim, idim, nummem, turre, q, qk0, fmu, tke, sk, vol, vk0, timestep, source, rhs,
                       issglrrw2012);
            afsolver_j(jdim, kdim, idim, nummem, turre, q, qj0, fmu, tke, sj, vol, vj0, timestep, source, rhs,
                       issglrrw2012);
            afsolver_i(jdim, kdim, idim, nummem, turre, q, qi0, fmu, tke, si, vol, vi0, timestep, source, rhs,
                       issglrrw2012);
            update(jdim, kdim, idim, nummem, rhs, turre, q, fmu, turb, vist3d, sumn, negn, ux, issglrrw2012);
        } else {
            sgs_solver_2d(jdim, kdim, idim, nummem, timestep, d, al, ar, bl, br, rhs);
            update(jdim, kdim, idim, nummem, rhs, turre, q, fmu, turb, vist3d, sumn, negn, ux, issglrrw2012);
        }

    }

    if (nfreq > 0) {
        if (ntt % nfreq == 0) {
            int nf = ntt / nfreq;
            kws_dump_movie(nf, jdim, kdim, idim, nummem, nbl, q, x, y, z, vist3d, ux, turb, turb, smin);
        }
    }
    omega.deallocate(); tke.deallocate(); blend.deallocate(); turre.deallocate();
    bij.deallocate(); fmu.deallocate(); timestep.deallocate(); d.deallocate();
    al.deallocate(); ar.deallocate(); bl.deallocate(); br.deallocate();
    dbijdx.deallocate(); dbijdxx.deallocate();
}

} // namespace module_kwstm_ns
