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

namespace chkrap_ns {

void chkrap(int& nbl, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& nface, int& idim, int& jdim, int& kdim, int& nn, int& maxxseg, int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend, int& iwrap, int& jwrap, int& kwrap, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid);

} // namespace chkrap_ns
