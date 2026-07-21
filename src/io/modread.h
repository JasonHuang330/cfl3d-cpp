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

namespace modread_ns {

void modread(int& idim, int& jdim, int& kdim, int& n, int& nbl, int& iunit, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, int& maxbl, int& maxseg, int& nmds, FortranArray5DRef<double> xmdj, FortranArray5DRef<double> xmdk, FortranArray5DRef<double> xmdi, int& iaes, int& maxaes);

} // namespace modread_ns
