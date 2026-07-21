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

namespace tfivol_ns {

void tfivol(int& idim, int& jdim, int& kdim, int& iskp, int& jskp, int& kskp, FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax, FortranArray1DRef<int> kskmax, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, int& isktyp, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> x1, FortranArray3DRef<double> y1, FortranArray3DRef<double> z1, FortranArray3DRef<double> x2, FortranArray3DRef<double> y2, FortranArray3DRef<double> z2, FortranArray3DRef<double> x3, FortranArray3DRef<double> y3, FortranArray3DRef<double> z3, FortranArray3DRef<double> x4, FortranArray3DRef<double> y4, FortranArray3DRef<double> z4, FortranArray3DRef<double> x5, FortranArray3DRef<double> y5, FortranArray3DRef<double> z5, FortranArray3DRef<double> x6, FortranArray3DRef<double> y6, FortranArray3DRef<double> z6, FortranArray3DRef<double> x7, FortranArray3DRef<double> y7, FortranArray3DRef<double> z7, FortranArray3DRef<double> arci, FortranArray3DRef<double> arcj, FortranArray3DRef<double> arck, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& maxbl, int& nbl);

} // namespace tfivol_ns
