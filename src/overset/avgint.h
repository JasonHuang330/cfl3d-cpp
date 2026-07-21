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

namespace avgint_ns {

void avgint(FortranArray1DRef<double> xiec, FortranArray1DRef<double> etac, FortranArray1DRef<int> nblkc, int& nptc, FortranArray1DRef<double> xief, FortranArray1DRef<double> etaf, FortranArray1DRef<int> nblkf, int& nptf, int& j1c, int& j2c, int& k1c, int& k2c, int& j1f, int& j2f, int& k1f, int& k2f);

} // namespace avgint_ns
