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

namespace blnkfr_ns {

void blnkfr(int& nbl, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> lbg, int& iitot, FortranArray3DRef<double> blank, int& jdim, int& kdim, int& idim, int& maxbl, double& blnkval);

} // namespace blnkfr_ns
