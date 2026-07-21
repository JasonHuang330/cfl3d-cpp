// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#pragma once
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <complex>
#include <cstdint>

namespace readkey_ns {

void readkey(int& ititr, int& myid, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou, int& iunit11, int& ierrflg);

} // namespace readkey_ns
