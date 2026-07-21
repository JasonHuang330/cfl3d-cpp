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

namespace rpatch_ns {

void rpatch(int& maxbl, int& maxxe, int& intmax, int& nsub1, FortranArray2DRef<double> windex, int& ninter, FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt, FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);

} // namespace rpatch_ns
