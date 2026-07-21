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

namespace invert_ns {

void invert(int& jdim1, int& kdim1, int& msub1, int& msub2, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, int& lmax1, FortranArray3DRef<double> x1, FortranArray3DRef<double> y1, FortranArray3DRef<double> z1, FortranArray3DRef<double> x1mid, FortranArray3DRef<double> y1mid, FortranArray3DRef<double> z1mid, FortranArray3DRef<double> x1mide, FortranArray3DRef<double> y1mide, FortranArray3DRef<double> z1mide, int& limit0, FortranArray1DRef<int> jjmax2, FortranArray1DRef<int> kkmax2, FortranArray3DRef<double> x2, FortranArray3DRef<double> y2, FortranArray3DRef<double> z2, FortranArray1DRef<double> xie2, FortranArray1DRef<double> eta2, FortranArray1DRef<int> mblkpt, FortranArray1DRef<double> temp, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, int& ifit, int& itmax, FortranArray3DRef<double> sxie, FortranArray3DRef<double> seta, FortranArray3DRef<double> sxie2, FortranArray3DRef<double> seta2, FortranArray2DRef<double> xie2s, FortranArray2DRef<double> eta2s, int& intmx, int& icheck, FortranArray1DRef<int> nblkj, FortranArray1DRef<int> nblkk, FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm, int& mcxie, int& mceta, FortranArray1DRef<int> lout, int& j21, int& j22, int& k21, int& k22, int& npt, int& ic0, int& iorph, int& itoss0, int& ncall, int& ioutpt, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, int& iself, int& ifiner, FortranArray1DRef<double> xie2f, FortranArray1DRef<double> eta2f, FortranArray1DRef<int> mblkptf, int& nptf, int& xi1f, int& xi2f, int& et1f, int& et2f, int& iavg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl);

} // namespace invert_ns
