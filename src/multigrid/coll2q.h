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

namespace coll2q_ns {

void coll2q(FortranArray4DRef<double> q, FortranArray4DRef<double> qc, FortranArray3DRef<double> vol, FortranArray3DRef<double> volc, int& jdim, int& kdim, int& idim, int& jc, int& kc, int& ic, FortranArray4DRef<double> res, FortranArray4DRef<double> qr, FortranArray4DRef<double> qw, int& js, int& ks, int& is, int& je, int& ke, int& ie, int& nbl, int& nblc, FortranArray3DRef<double> vistf, FortranArray3DRef<double> vistc, FortranArray4DRef<double> tursavf, FortranArray4DRef<double> tursavc, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& nummem);

} // namespace coll2q_ns
