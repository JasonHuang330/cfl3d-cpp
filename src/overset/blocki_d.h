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

namespace blocki_d_ns {

void blocki(FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, int& idimr, int& jdimr, int& kdimr, int& jdimt, int& kdimt, FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva, int& it, int& ir, int& ldim, FortranArray3DRef<double> bci, int& iedge, int& ivolflg);
void blocki_d(FortranArray3DRef<double> q, FortranArray4DRef<double> qi0, int& idimr, int& jdimr, int& kdimr, int& jdimt, int& kdimt, FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva, int& it, int& ir, int& nvals, int& ldim, FortranArray3DRef<double> bci, int& iedge);

} // namespace blocki_d_ns
