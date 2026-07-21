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

namespace rrest_ns {

void rrest(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, int& ncycmax, int& ntr, FortranArray1DRef<double> rms, FortranArray1DRef<double> clw, FortranArray1DRef<double> cdw, FortranArray1DRef<double> cdpw, FortranArray1DRef<double> cdvw, FortranArray1DRef<double> cxw, FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw, FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw, FortranArray1DRef<double> cmzw, int& n_clcd, FortranArray3DRef<double> clcd, int& nblocks_clcd, FortranArray2DRef<int> blocks_clcd, FortranArray1DRef<double> fmdotw, FortranArray1DRef<double> cftmomw, FortranArray1DRef<double> cftpw, FortranArray1DRef<double> cftvw, FortranArray1DRef<double> cfttotw, FortranArray2DRef<double> rmstr, FortranArray2DRef<int> nneg, int& iskip, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> tursav, FortranArray3DRef<double> smin, FortranArray3DRef<double> xjb, FortranArray3DRef<double> xkb, FortranArray3DRef<double> blnum, FortranArray3DRef<double> cmuv, int& maxbl, FortranArray1DRef<int> mblk2nd, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& igrid, FortranArray3DRef<double> wk, int& idima, int& jdima, int& kdima, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, FortranArray4DRef<double> qavg, FortranArray4DRef<double> q2avg, int& nummem);

} // namespace rrest_ns
