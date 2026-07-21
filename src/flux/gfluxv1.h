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

namespace gfluxv1_ns {

void gfluxv(int& i, int& npl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> sj, FortranArray3DRef<double> vol, FortranArray2DRef<double> t, int& nvtq, FortranArray2DRef<double> wj0, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vmuj, FortranArray4DRef<double> vj0, FortranArray3DRef<double> bcj, FortranArray4DRef<double> zksav, FortranArray4DRef<double> tj0, FortranArray3DRef<double> cmuv, FortranArray3DRef<double> volj0, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iadv, int& nummem, FortranArray4DRef<double> ux);
void gfluxv1(int& i, int& npl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, FortranArray2DRef<double> t, int& nvtq, FortranArray2DRef<double> wj0, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> vj0, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray3DRef<double> volj0, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iadv);

} // namespace gfluxv1_ns
