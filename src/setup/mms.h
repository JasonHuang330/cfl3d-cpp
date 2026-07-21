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

namespace mms_ns {

void exact_norm(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> q, FortranArray3DRef<double> tursav, FortranArray3DRef<double> vol, FortranArray4DRef<double> res, FortranArray3DRef<double> smin, int& iexact_trunc, int& iexact_disc);
void exact_turb_res(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> vol, FortranArray3DRef<double> res, int& iexact_trunc, int& iexact_disc);
void exact_flow_force(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vol, FortranArray4DRef<double> res, FortranArray3DRef<double> smin, int& iexact_trunc, int& iexact_disc);
void exact_turb_force(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vol, FortranArray3DRef<double> res, FortranArray3DRef<double> smin, int& iexact_trunc, int& iexact_disc);
void exact_flow_q(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> q, int& iexact_trunc, int& iexact_disc);
void zero_resid_ring(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, int& jj, int& kk, int& ii, int& numeq, int& numrows, int& iexact_trunc, int& iexact_disc);
void exact_flow_q_ring(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> q, int& iexact_trunc, int& iexact_disc);
void exact_turb_q(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> q, FortranArray3DRef<double> smin, FortranArray3DRef<double> vist3d, int& iexact_trunc, int& iexact_disc);
void exact_turb_q_ring(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> q, FortranArray3DRef<double> smin, FortranArray3DRef<double> vist3d, int& iexact_trunc, int& iexact_disc);
void exact_flow_bc(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, int& nface, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray3DRef<double> smin, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, int& nummem, int& iexact_trunc, int& iexact_disc);
void analytic_compressible(double& x, double& y, double& z, int& neq, FortranArray1DRef<double> q, int& i_convert_q, int& i_forcing, int& i_gradient, double& distf, double& xmut, int& iexact_trunc, int& iexact_disc);
void exact_zero(double& q, double& qx, double& qy, double& qz, double& qxx, double& qyy, double& qzz, double& qxy, double& qxz, double& qyz);
void total_temperature_2d(double& t, double& tx, double& ty, double& txx, double& tyy, double& txy, double& u, double& ux, double& uy, double& uxx, double& uyy, double& uxy, double& v, double& vx, double& vy, double& vxx, double& vyy, double& vxy);
double viscosity_law(double& t);
double dviscosity_law(double& t);
void exact_lisbon_ms2_forcing(double& beta, double& x, double& y, FortranArray1DRef<double> f);
void exact_lisbon_ms2_u(double& scale, double& x, double& y, double& q, double& qx, double& qy, double& qxx, double& qyy, double& qxy);
void exact_lisbon_ms2_v(double& scale, double& x, double& y, double& q, double& qx, double& qy, double& qxx, double& qyy, double& qxy);
void exact_lisbon_ms2_p(double& scale, double& x, double& y, double& q, double& qx, double& qy);
void exact_lisbon_ms1_vt(double& scale, double& x, double& y, double& q, double& qx, double& qy, double& qxx, double& qyy);
void exact_lisbon_ms2_vt(double& scale, double& x, double& y, double& q, double& qx, double& qy, double& qxx, double& qyy);
void exact_lisbon_ms4_vt(double& scale, double& x, double& y, double& q, double& qx, double& qy, double& qxx, double& qyy);
double ums(double& x, double& y);
double dudxms(double& x, double& y);
double dudyms(double& x, double& y);
double dudx2ms(double& x, double& y);
double dudy2ms(double& x, double& y);
double dudxyms(double& x, double& y);
double dudx2yms(double& x, double& y);
double dudy3ms(double& x, double& y);
double vms(double& x, double& y);
double dvdxms(double& x, double& y);
double dvdyms(double& x, double& y);
double dvdx2ms(double& x, double& y);
double dvdy2ms(double& x, double& y);
double dvdxyms(double& x, double& y);
double dvdxy2ms(double& x, double& y);
double dvdx3ms(double& x, double& y);
double pms(double& x, double& y);
double dpdxms(double& x, double& y);
double dpdyms(double& x, double& y);
double strainms(double& x, double& y);
double vortms(double& x, double& y);
double dvodxms(double& x, double& y);
double dvodyms(double& x, double& y);
double dvodx2ms(double& x, double& y);
double dvody2ms(double& x, double& y);
double smxsams(double& x, double& y);
double smysams(double& x, double& y);
double smxmtms(double& x, double& y);
double smymtms(double& x, double& y);
double em1ms(double& x, double& y);
double em2ms(double& x, double& y);
double em4ms(double& x, double& y);
double dem1dxms(double& x, double& y);
double dem2dxms(double& x, double& y);
double dem4dxms(double& x, double& y);
double dem1dyms(double& x, double& y);
double dem2dyms(double& x, double& y);
double dem4dyms(double& x, double& y);
double dem1dx2ms(double& x, double& y);
double dem2dx2ms(double& x, double& y);
double dem4dx2ms(double& x, double& y);
double dem1dy2ms(double& x, double& y);
double dem2dy2ms(double& x, double& y);
double dem4dy2ms(double& x, double& y);
double ssams(double& x, double& y);
double eddysams(double& x, double& y);
double desadxms(double& x, double& y);
double desadyms(double& x, double& y);
double fv1sams(double& em);
double dfv1sams(double& em);
double smtms(double& x, double& y);
double eddymtms(double& x, double& y);
double demtdxms(double& x, double& y);
double demtdyms(double& x, double& y);
double d2mtms(double& em);
double dd2mtms(double& em);
double user_derfc(double& x);
double user_derf(double& x);
void exact_polyfg(double& scale, FortranArray1DRef<double> polyfc, double& x, double& y, double& z, double& q, double& qx, double& qy, double& qz, double& qxx, double& qyy, double& qzz, double& qxy, double& qxz, double& qyz);

} // namespace mms_ns
