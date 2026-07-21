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

namespace i2xs_ns {

void i2x(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc, int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf, FortranArray4DRef<double> qkf, FortranArray4DRef<double> qif, int& js, int& ks, int& is, int& je, int& ke, int& ie, int& nblc, int& ldim, int& nbl, FortranArray3DRef<double> bcjf, FortranArray3DRef<double> bckf, FortranArray3DRef<double> bcif, int& nface);
void i2xs(int& jc, int& kc, int& ic, FortranArray4DRef<double> qc, int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf, FortranArray4DRef<double> qkf, FortranArray4DRef<double> qif, int& js, int& ks, int& is, int& je, int& ke, int& ie, int& nblc, int& ldim, int& nbl, FortranArray3DRef<double> bcjf, FortranArray3DRef<double> bckf, FortranArray3DRef<double> bcif, int& nface);

} // namespace i2xs_ns
