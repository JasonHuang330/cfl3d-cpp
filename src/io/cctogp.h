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

namespace cctogp_ns {

void cctogp(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3, int& j1, int& j2, int& j3, int& k1, int& k2, int& k3, FortranArray4DRef<double> dum, FortranArray4DRef<double> dumi0, FortranArray4DRef<double> dumj0, FortranArray4DRef<double> dumk0, int& jdw, int& kdw, int& idw, FortranArray4DRef<double> dumgp, int& ldim);

} // namespace cctogp_ns
