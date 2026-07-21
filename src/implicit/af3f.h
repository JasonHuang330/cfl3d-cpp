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

namespace af3f_ns {

void af3f(int& nbl, int& jdim, int& kdim, int& idim, FortranArray2DRef<double> q, FortranArray1DRef<double> vol, FortranArray3DRef<double> qj0, FortranArray3DRef<double> qk0, FortranArray3DRef<double> qi0, FortranArray1DRef<double> dtj, FortranArray2DRef<double> sj, FortranArray2DRef<double> sk, FortranArray2DRef<double> si, FortranArray2DRef<double> res, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> blank, FortranArray3DRef<double> vmuk, double& resid, FortranArray1DRef<double> wk, int& nwork, FortranArray1DRef<double> wk0, int& nwk0, int& iover, FortranArray4DRef<double> vk0, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, FortranArray3DRef<double> volk0, FortranArray4DRef<double> tursav, FortranArray4DRef<double> tk0, FortranArray3DRef<double> cmuv, FortranArray1DRef<int> iadvance, int& nummem, FortranArray4DRef<double> ux);

} // namespace af3f_ns
