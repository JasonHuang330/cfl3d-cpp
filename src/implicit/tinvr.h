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

namespace tinvr_ns {

void tinvr(int& n, FortranArray1DRef<double> r1, FortranArray1DRef<double> r2, FortranArray1DRef<double> r3, FortranArray1DRef<double> r4, FortranArray1DRef<double> r5, FortranArray1DRef<double> kx, FortranArray1DRef<double> ky, FortranArray1DRef<double> kz, FortranArray1DRef<double> lx, FortranArray1DRef<double> ly, FortranArray1DRef<double> lz, FortranArray1DRef<double> mx, FortranArray1DRef<double> my, FortranArray1DRef<double> mz, FortranArray1DRef<double> c, FortranArray1DRef<double> ub, FortranArray1DRef<double> rho, FortranArray1DRef<double> u, FortranArray1DRef<double> v, FortranArray1DRef<double> w, int& max, int& itinv, FortranArray1DRef<double> eig2, FortranArray1DRef<double> eig3, FortranArray1DRef<double> xm2a);

} // namespace tinvr_ns
