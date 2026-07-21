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

namespace fa2xi_ns {

void fa(int& nbl, FortranArray2DRef<int> lw, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, FortranArray1DRef<int> iw, int& ifamax, int& nwfa, int& nifa, FortranArray1DRef<int> nfajki, int& maxbl, int& maxseg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, int& nblock, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim);
void fa2xi(int& jf, int& kf, int& if_, FortranArray4DRef<double> qif, int& js, int& ks, int& is, int& je, int& ke, int& ie, FortranArray1DRef<int> ibctyp, FortranArray1DRef<double> w, FortranArray1DRef<int> iw, int& nw, int& niw, int& neta, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim);

} // namespace fa2xi_ns
