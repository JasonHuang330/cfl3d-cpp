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

namespace collmod_ns {

void collmod(FortranArray5DRef<double> xmdj, FortranArray5DRef<double> xmdk, FortranArray5DRef<double> xmdi, FortranArray5DRef<double> xmdjc, FortranArray5DRef<double> xmdkc, FortranArray5DRef<double> xmdic, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, int& nm, int& nmds, int& iaes, int& maxaes);

} // namespace collmod_ns
