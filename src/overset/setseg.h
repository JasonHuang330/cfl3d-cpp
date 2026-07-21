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

namespace setseg_ns {

void setseg(int& maxgr, int& maxbl, int& maxseg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, int& ndefrm0, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, int& maxsegdg, int& ngrid);

} // namespace setseg_ns
