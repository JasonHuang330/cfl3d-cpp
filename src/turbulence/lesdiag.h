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

namespace lesdiag_ns {

void lesdiag(int& myid, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> ux, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vol, FortranArray4DRef<double> si, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray3DRef<double> vor, FortranArray3DRef<double> smin, FortranArray3DRef<double> xjb, FortranArray4DRef<double> tursav, FortranArray3DRef<double> xkb, FortranArray3DRef<double> blnum, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& nbl, int& nummem, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z);

} // namespace lesdiag_ns
