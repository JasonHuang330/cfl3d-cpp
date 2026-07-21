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

namespace csout_ns {

void csout(int& iseq, int& maxbl, int& maxcs, FortranArray1DRef<int> igridg, FortranArray1DRef<int> levelg, int& ncs, FortranArray1DRef<double> sx, FortranArray1DRef<double> sy, FortranArray1DRef<double> sz, FortranArray1DRef<double> stot, FortranArray1DRef<double> pav, FortranArray1DRef<double> ptav, FortranArray1DRef<double> tav, FortranArray1DRef<double> ttav, FortranArray1DRef<double> xmav, FortranArray1DRef<double> fmdot, FortranArray1DRef<double> cfxp, FortranArray1DRef<double> cfyp, FortranArray1DRef<double> cfzp, FortranArray1DRef<double> cfdp, FortranArray1DRef<double> cflp, FortranArray1DRef<double> cftp, FortranArray1DRef<double> cfxv, FortranArray1DRef<double> cfyv, FortranArray1DRef<double> cfzv, FortranArray1DRef<double> cfdv, FortranArray1DRef<double> cflv, FortranArray1DRef<double> cftv, FortranArray1DRef<double> cfxmom, FortranArray1DRef<double> cfymom, FortranArray1DRef<double> cfzmom, FortranArray1DRef<double> cfdmom, FortranArray1DRef<double> cflmom, FortranArray1DRef<double> cftmom, FortranArray1DRef<double> cfxtot, FortranArray1DRef<double> cfytot, FortranArray1DRef<double> cfztot, FortranArray1DRef<double> cfdtot, FortranArray1DRef<double> cfltot, FortranArray1DRef<double> cfttot, FortranArray2DRef<int> icsinfo);

} // namespace csout_ns
