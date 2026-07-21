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

namespace threeeqn_ns {

void threeeqn(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, FortranArray3DRef<double> dtj, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vor, FortranArray3DRef<double> smin, FortranArray4DRef<double> zksav, FortranArray4DRef<double> turre, FortranArray3DRef<double> damp1, FortranArray3DRef<double> blend, FortranArray3DRef<double> timestp, FortranArray3DRef<double> fnu, FortranArray2DRef<double> bx, FortranArray2DRef<double> bx2, FortranArray2DRef<double> cx, FortranArray2DRef<double> cx2, FortranArray2DRef<double> dx, FortranArray2DRef<double> dx2, FortranArray2DRef<double> fx, FortranArray2DRef<double> fx2, FortranArray2DRef<double> workx, FortranArray2DRef<double> by, FortranArray2DRef<double> by2, FortranArray2DRef<double> cy, FortranArray2DRef<double> cy2, FortranArray2DRef<double> dy, FortranArray2DRef<double> dy2, FortranArray2DRef<double> fy, FortranArray2DRef<double> fy2, FortranArray2DRef<double> worky, FortranArray2DRef<double> bz, FortranArray2DRef<double> bz2, FortranArray2DRef<double> cz, FortranArray2DRef<double> cz2, FortranArray2DRef<double> dz, FortranArray2DRef<double> dz2, FortranArray2DRef<double> fz, FortranArray2DRef<double> fz2, FortranArray2DRef<double> workz, int& ntime, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, int& nbl, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, FortranArray3DRef<double> blank, int& iover, double& sumn1, double& sumn2, double& sumn3, int& negn1, int& negn2, int& negn3, FortranArray4DRef<double> ux, FortranArray4DRef<double> rhside, FortranArray4DRef<double> zksav2, FortranArray3DRef<double> v3dtmp, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& maxbl, int& maxseg, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iex, int& iex2, int& iex3, FortranArray2DRef<double> bx3, FortranArray2DRef<double> cx3, FortranArray2DRef<double> dx3, FortranArray2DRef<double> fx3, FortranArray2DRef<double> by3, FortranArray2DRef<double> cy3, FortranArray2DRef<double> dy3, FortranArray2DRef<double> fy3, FortranArray2DRef<double> bz3, FortranArray2DRef<double> cz3, FortranArray2DRef<double> dz3, FortranArray2DRef<double> fz3, int& nummem);

} // namespace threeeqn_ns
