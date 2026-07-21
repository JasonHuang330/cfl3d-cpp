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

namespace l2norm2_ns {

void l2norm(int& nbl, int& ntime, double& rmsl, int& irdq, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, FortranArray3DRef<double> vol);
void l2norm2(int& nbl, int& ntime, double& rmsl, int& irdq, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, FortranArray3DRef<double> vol, FortranArray4DRef<double> qc0, FortranArray4DRef<double> dqc0, FortranArray4DRef<double> q, FortranArray3DRef<double> blank);
void resadd(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qc0, FortranArray4DRef<double> dqc0, FortranArray4DRef<double> res, FortranArray3DRef<double> vol, int& iover, FortranArray3DRef<double> blank);

} // namespace l2norm2_ns
