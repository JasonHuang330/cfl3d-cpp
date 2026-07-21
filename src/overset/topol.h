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

namespace topol_ns {

void topol(int& jdim, int& kdim, int& nsub, FortranArray1DRef<int> jjmax, FortranArray1DRef<int> kkmax, int& lmax, int& l, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> xmid, FortranArray3DRef<double> ymid, FortranArray3DRef<double> zmid, FortranArray3DRef<double> xmide, FortranArray3DRef<double> ymide, FortranArray3DRef<double> zmide, int& limit, double& xc, double& yc, double& zc, double& xie, double& eta, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, int& ifit, int& itmax, int& igap, int& iok, FortranArray1DRef<int> lout, int& ic0, int& itoss0, int& jto, int& kto, int& iself, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);

} // namespace topol_ns
