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

namespace spalart_ns {

void init(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, FortranArray3DRef<double> vol, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, int& nummem, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iflagprnt);
void spalart(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, FortranArray3DRef<double> dtj, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vor, FortranArray3DRef<double> smin, FortranArray4DRef<double> tursav, FortranArray3DRef<double> turre, FortranArray3DRef<double> damp1, FortranArray3DRef<double> timestp, FortranArray3DRef<double> fnu, FortranArray2DRef<double> bx, FortranArray2DRef<double> cx, FortranArray2DRef<double> dx, FortranArray2DRef<double> fx, FortranArray2DRef<double> workx, FortranArray2DRef<double> by, FortranArray2DRef<double> cy, FortranArray2DRef<double> dy, FortranArray2DRef<double> fy, FortranArray2DRef<double> worky, FortranArray2DRef<double> bz, FortranArray2DRef<double> cz, FortranArray2DRef<double> dz, FortranArray2DRef<double> fz, FortranArray2DRef<double> workz, int& ntime, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, int& nbl, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray3DRef<double> blank, int& iover, double& sumn1, double& sumn2, int& negn1, int& negn2, FortranArray4DRef<double> tursav2, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iex, int& iex2, int& iex3, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray4DRef<double> ux, FortranArray4DRef<double> vx, int& nummem, FortranArray3DRef<double> fd, int& isarcnum);

} // namespace spalart_ns
