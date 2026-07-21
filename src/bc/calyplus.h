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

namespace calyplus_ns {

void calyplus(int& jdim, int& kdim, int& idim, int& nbl, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, int& iover, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> yplusj, FortranArray3DRef<double> yplusk, FortranArray3DRef<double> yplusi, FortranArray3DRef<double> blankj, FortranArray3DRef<double> blankk, FortranArray3DRef<double> blanki, FortranArray3DRef<double> dnj, FortranArray3DRef<double> dnk, FortranArray3DRef<double> dni, FortranArray3DRef<double> vistj, FortranArray3DRef<double> vistk, FortranArray3DRef<double> visti, int& ifunc, int& iunit, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, int& maxbl, int& maxseg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray3DRef<double> vol);

} // namespace calyplus_ns
