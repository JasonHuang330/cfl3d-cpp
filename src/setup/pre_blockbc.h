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

namespace pre_blockbc_ns {

void pre_blockbc(int& nbl, FortranArray2DRef<int> lw, int& icount, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray2DRef<int> isav_blk, FortranArray2DRef<int> nblk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, int& mxbli, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& maxbl, int& ierrflg);

} // namespace pre_blockbc_ns
