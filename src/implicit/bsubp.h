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

namespace bsubp_ns {

void bsub(int& id1, int& id2, FortranArray4DRef<double> a, FortranArray4DRef<double> b, FortranArray4DRef<double> c, FortranArray3DRef<double> f, int& i1, int& i2, int& il, int& iu);
void bsubp(int& id1, int& id2, FortranArray4DRef<double> a, FortranArray4DRef<double> b, FortranArray4DRef<double> c, FortranArray3DRef<double> f, int& i1, int& i2, int& il, int& iu, FortranArray4DRef<double> g, FortranArray4DRef<double> h);

} // namespace bsubp_ns
