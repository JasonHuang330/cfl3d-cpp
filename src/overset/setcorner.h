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

namespace setcorner_ns {

void setcorner(int& j, int& k, int& i, FortranArray3DRef<double> xnm1, FortranArray3DRef<double> ynm1, FortranArray3DRef<double> znm1, FortranArray3DRef<double> dx, FortranArray3DRef<double> dy, FortranArray3DRef<double> dz, int& jdim, int& kdim, int& idim, FortranArray1DRef<double> wk, int& nsurf, int& iflag, int& ivert, FortranArray3DRef<int> islavept, int& nslave, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& ibufdim, int& nbuf, int& myid, FortranArray1DRef<int> ibl, int& nmaster, int& iseq);

} // namespace setcorner_ns
