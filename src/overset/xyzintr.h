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

namespace xyzintr_ns {

void xyzintr(int& jdimf, int& kdimf, int& idimf, int& jdimc, int& kdimc, int& idimc, FortranArray3DRef<double> xf, FortranArray3DRef<double> yf, FortranArray3DRef<double> zf, FortranArray3DRef<double> xc, FortranArray3DRef<double> yc, FortranArray3DRef<double> zc);

} // namespace xyzintr_ns
