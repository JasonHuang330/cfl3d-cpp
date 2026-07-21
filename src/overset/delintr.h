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

namespace delintr_ns {

void delintr(int& jdimf, int& kdimf, int& idimf, int& jdimc, int& kdimc, int& idimc, FortranArray4DRef<double> deltjf, FortranArray4DRef<double> deltkf, FortranArray4DRef<double> deltif, FortranArray4DRef<double> deltjc, FortranArray4DRef<double> deltkc, FortranArray4DRef<double> deltic);

} // namespace delintr_ns
