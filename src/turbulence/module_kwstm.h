// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#pragma once
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <complex>
#include <cstdint>

namespace module_kwstm_ns {

// module variables (shared global state)
inline double c1 = 9./5.;
inline double c2 = 10./19.;
inline double alpha = 0.52;
inline double beta_0 = 0.0708;
inline double beta_star = 0.09;
inline double sigma = 0.5;
inline double sigma_star = 0.6;
inline double sigma_d0 = 0.125;
inline double alpha_o = 0.5556;
inline double beta_o = 0.075;
inline double sigma_w_o = 0.5;
inline double sigma_d_o = 0.;
inline double c1_o = 1.8;
inline double c1star_o = 0.;
inline double c2_o = 0.;
inline double c3_o = 0.8;
inline double c3star_o = 0.;
inline double c4_o = 0.97091;
inline double c5_o = 0.57818;
inline double d_o = 0.0675;
inline double d_sd_o = 0.5;
inline double alpha_e = 0.44;
inline double beta_e = 0.0828;
inline double sigma_w_e = 0.856;
inline double sigma_d_e = 1.712;
inline double c1_e = 1.7;
inline double c1star_e = 0.9;
inline double c2_e = 1.05;
inline double c3_e = 0.8;
inline double c3star_e = 0.65;
inline double c4_e = 0.625;
inline double c5_e = 0.2;
inline double d_e = 0.22;
inline double d_sd_e = 1.62963;
inline double alpha_hat;
inline double beta_hat;
inline double gamma_hat;
inline int nsubit = 1;
inline double factor = 1.0;
inline int irealizability = 1;
inline double pdratio = 10.0;
inline int nfreq = 100000;
inline bool need_contplot = false;
inline bool need_kwsplot = false;
inline int ux_signal = 0;
inline int myrank;
inline int myroot;
inline int mylevel;
inline int myicyc;
inline int level_o = -1;
inline int icyc_o = 1;
inline int debug_j = -11;
inline int jcut = 64;
inline int jstart = 17;
inline int kdiff = 0;
inline int iopt = 0;
inline double cfl_psd = 1.0;
inline double cfl_loc = 10.0;
inline double xmutratio = 0.0;
inline double ildts = 0;
inline bool lkzstm_from_2eq = false;

// module arrays (owning storage)
inline FortranArray3D<double> src_k;
inline FortranArray5D<double> source_items;
inline FortranArray4D<double> source;
inline FortranArray4D<double> rhs;

double kws_get_nummem();
void kws_initvist(int& jdim, int& kdim, int& idim, int& nummem, FortranArray4DRef<double> cij, FortranArray3DRef<double> vist3d);
void kws_init(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim);
void kws_read_parm(FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim);
void kws_main(int& jdim, int& kdim, int& idim, int& myid, int& myhost, int& nummem, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, FortranArray4DRef<double> turb, FortranArray4DRef<double> ux, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, FortranArray3DRef<double> vj0, FortranArray3DRef<double> vk0, FortranArray3DRef<double> vi0, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray3DRef<double> dtj, FortranArray3DRef<double> vist3d, int& level, int& icyc, FortranArray1DRef<double> sumn, FortranArray1DRef<int> negn, FortranArray4DRef<double> zksav2, FortranArray3DRef<double> smin, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& nbl, int& issglrrw2012, int& i_sas_rsm, int& i_yapterm);
void kws_plot(int& nframe, int& jdim, int& kdim, int& idim, int& nummem, int& nblk, FortranArray4DRef<double> q, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> ux, FortranArray4DRef<double> turre, FortranArray3DRef<double> smin);
void kws_dump_movie(int& nframe, int& jdim, int& kdim, int& idim, int& nummem, int& nblk, FortranArray4DRef<double> q, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> ux, FortranArray4DRef<double> turre, FortranArray4DRef<double> tursav, FortranArray3DRef<double> smin);
void kws_init_from_2eq(bool& lfrom2eq, int& nbl, int& jdim, int& kdim, int& idim, int& nummem, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> tursav);
void kws_save_2eq(int& nbl, int& jdim, int& kdim, int& idim, int& nummem, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> tursav);

} // namespace module_kwstm_ns
