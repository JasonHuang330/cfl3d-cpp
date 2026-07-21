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

namespace histout_ns {

void histout(int& ihstry, FortranArray1DRef<double> rms, FortranArray1DRef<double> clw, FortranArray1DRef<double> cdw, FortranArray1DRef<double> cdpw, FortranArray1DRef<double> cdvw, FortranArray1DRef<double> cxw, FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw, FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw, FortranArray1DRef<double> cmzw, int& n_clcd, FortranArray3DRef<double> clcd, int& nblocks_clcd, FortranArray2DRef<int> blocks_clcd, FortranArray1DRef<double> chdw, FortranArray1DRef<double> swetw, FortranArray1DRef<double> fmdotw, FortranArray1DRef<double> cfttotw, FortranArray1DRef<double> cftmomw, FortranArray1DRef<double> cftpw, FortranArray1DRef<double> cftvw, FortranArray2DRef<double> rmstr, FortranArray2DRef<int> nneg, int& ncycmax, FortranArray4DRef<double> aehist, FortranArray2DRef<double> aesrfdat, int& nmds, int& maxaes, FortranArray1DRef<double> timekeep, int& nummem);

} // namespace histout_ns
