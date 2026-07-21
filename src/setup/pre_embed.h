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

namespace pre_embed_ns {

void pre_embed(int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, int& icount, int& maxbl, int& maxseg, int& lbcemb, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbcidim, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, FortranArray1DRef<int> iemg, FortranArray1DRef<int> igridg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<int> isav_emb, FortranArray1DRef<int> is_emb, FortranArray1DRef<int> ie_emb, int& nbcemb, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& maxgr, int& ierrflg);

} // namespace pre_embed_ns
