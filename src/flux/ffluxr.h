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

namespace ffluxr_ns {

void ffluxr(int& k, int& npl, double& xkap, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray3DRef<double> si, FortranArray2DRef<double> t, int& nvtq, int& nv, int& nfa, FortranArray1DRef<double> wfa, FortranArray1DRef<int> iwfa, FortranArray1DRef<int> ibctyp, int& isf, int& nbl, FortranArray3DRef<double> bci, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, int& idef);

} // namespace ffluxr_ns
