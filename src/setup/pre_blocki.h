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

namespace pre_blocki_ns {

void pre_blocki(int& idimr, int& jdimr, int& kdimr, FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva, int& it, int& ir, FortranArray4DRef<double> q, FortranArray3DRef<double> qlocal, int& maxdims, int& ldim, int& ij, int& ivolflg, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, FortranArray2DRef<char[120]> bou, int& ibufdim, int& nbuf, FortranArray1DRef<int> nou);

} // namespace pre_blocki_ns
