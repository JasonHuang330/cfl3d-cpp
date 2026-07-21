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

namespace blomax_ns {

void blomax(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> vor, FortranArray3DRef<double> snj0, FortranArray3DRef<double> snk0, FortranArray3DRef<double> sni0, FortranArray3DRef<double> snjm, FortranArray3DRef<double> snkm, FortranArray3DRef<double> snim, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> eoms, int& iprint, int& inmx, FortranArray1DRef<double> eomui, FortranArray1DRef<double> fbl, FortranArray1DRef<double> rhon, FortranArray1DRef<double> amun, FortranArray1DRef<double> vortn, FortranArray1DRef<double> disn, FortranArray1DRef<double> utot, FortranArray1DRef<double> eomuo, FortranArray1DRef<double> damp, int& nblt, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> blank, int& iover, FortranArray3DRef<double> bci, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl);

} // namespace blomax_ns
