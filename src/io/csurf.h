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

namespace csurf_ns {

void csurf(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> sk, FortranArray4DRef<double> sj, FortranArray4DRef<double> si, FortranArray4DRef<double> q, int& ics, FortranArray3DRef<double> ub, FortranArray3DRef<double> vb, FortranArray3DRef<double> wb, FortranArray3DRef<double> vmuk, FortranArray3DRef<double> vmuj, FortranArray3DRef<double> vmui, FortranArray3DRef<double> vol, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray3DRef<double> blank, FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk, FortranArray4DRef<double> xtbi, int& iuns, int& ncs, FortranArray2DRef<int> icsinfo, double& sx, double& sy, double& sz, double& stot, double& pav, double& ptav, double& tav, double& ttav, double& xmav, double& fmdot, double& cfxp, double& cfyp, double& cfzp, double& cfdp, double& cflp, double& cftp, double& cfxv, double& cfyv, double& cfzv, double& cfdv, double& cflv, double& cftv, double& cfxmom, double& cfymom, double& cfzmom, double& cfdmom, double& cflmom, double& cftmom, double& cfxtot, double& cfytot, double& cfztot, double& cfdtot, double& cfltot, double& cfttot, int& maxcs, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0);

} // namespace csurf_ns
