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

namespace genforce_ns {

void force(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> sk, FortranArray4DRef<double> sj, FortranArray4DRef<double> si, double& cl, double& cd, double& cz, double& cy, double& cx, double& cmy, double& cmx, double& cmz, double& chd, double& swet, int& i00, FortranArray3DRef<double> ub, FortranArray3DRef<double> vb, FortranArray3DRef<double> wb, FortranArray3DRef<double> vmuk, FortranArray3DRef<double> vmuj, FortranArray3DRef<double> vmui, FortranArray3DRef<double> vol, int& ifo, int& jfo, int& kfo, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray3DRef<double> blank, int& nbl, FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk, FortranArray4DRef<double> xtbi, int& iuns, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& nn, int& maxseg);
void genforce(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> sk, FortranArray4DRef<double> sj, FortranArray4DRef<double> si, double& cz, double& cy, double& cx, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray3DRef<double> blank, int& nbl, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, int& maxbl, int& maxseg, int& n2, FortranArray4DRef<double> xmdj, FortranArray4DRef<double> xmdk, FortranArray4DRef<double> xmdi, FortranArray2DRef<double> aesrfdat, int& nmds, int& maxaes, int& maxsegdg, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> idfrmseg, int& iaes, FortranArray2DRef<int> iaesurf);

} // namespace genforce_ns
