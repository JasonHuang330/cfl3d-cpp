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

namespace forceout_ns {

void force(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> sk, FortranArray4DRef<double> sj, FortranArray4DRef<double> si, double& cl, double& cd, double& cz, double& cy, double& cx, double& cmy, double& cmx, double& cmz, double& chd, double& swet, int& i00, FortranArray3DRef<double> ub, FortranArray3DRef<double> vb, FortranArray3DRef<double> wb, FortranArray3DRef<double> vmuk, FortranArray3DRef<double> vmuj, FortranArray3DRef<double> vmui, FortranArray3DRef<double> vol, int& ifo, int& jfo, int& kfo, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray3DRef<double> blank, int& nbl, FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk, FortranArray4DRef<double> xtbi, int& iuns, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& nn, int& maxseg);
void forceout(int& iseq, int& maxbl, int& maxgr, int& maxseg, int& nblock, FortranArray1DRef<int> iforce, FortranArray1DRef<int> igridg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> levelg, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray1DRef<double> swett, FortranArray1DRef<double> clt, FortranArray1DRef<double> cdt, FortranArray1DRef<double> cxt, FortranArray1DRef<double> cyt, FortranArray1DRef<double> czt, FortranArray1DRef<double> cmxt, FortranArray1DRef<double> cmyt, FortranArray1DRef<double> cmzt, FortranArray1DRef<double> cdpt, FortranArray1DRef<double> cdvt, FortranArray1DRef<double> swetw, FortranArray1DRef<double> clw, FortranArray1DRef<double> cdw, FortranArray1DRef<double> cxw, FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw, FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw, FortranArray1DRef<double> cmzw, FortranArray1DRef<double> cdpw, FortranArray1DRef<double> cdvw, int& ncycmax, int& myhost, int& myid, int& mycomm, FortranArray1DRef<int> mblk2nd);

} // namespace forceout_ns
