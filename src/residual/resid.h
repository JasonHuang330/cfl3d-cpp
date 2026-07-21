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

namespace resid_ns {

void resid(int& nbl, int& ntime, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q,
           FortranArray3DRef<double> qj0,
           FortranArray3DRef<double> qk0,
           FortranArray3DRef<double> qi0,
           FortranArray2DRef<double> sj,
           FortranArray2DRef<double> sk,
           FortranArray2DRef<double> si,
           FortranArray1DRef<double> vol,
           FortranArray1DRef<double> dtj,
           FortranArray1DRef<double> x,
           FortranArray1DRef<double> y,
           FortranArray1DRef<double> z,
           FortranArray3DRef<double> vist3d,
           FortranArray3DRef<double> snj0,
           FortranArray3DRef<double> snk0,
           FortranArray3DRef<double> sni0,
           FortranArray2DRef<double> res,
           FortranArray1DRef<double> wj0,
           FortranArray1DRef<double> wk0,
           FortranArray1DRef<double> wi0,
           FortranArray3DRef<double> vmuk,
           FortranArray3DRef<double> vmuj,
           FortranArray3DRef<double> vmui,
           FortranArray1DRef<double> wk,
           int& nwork,
           int& isf,
           FortranArray1DRef<int> iwfa,
           FortranArray1DRef<double> wfa,
           double& delt,
           FortranArray3DRef<double> blank,
           int& iover,
           int& nblendg,
           int& nblstat,
           int& nblstag,
           FortranArray4DRef<double> xib,
           FortranArray2DRef<double> sig,
           FortranArray2DRef<double> sqtq,
           FortranArray2DRef<double> g,
           FortranArray4DRef<double> tj0,
           FortranArray4DRef<double> tk0,
           FortranArray4DRef<double> ti0,
           FortranArray3DRef<double> xkb,
           FortranArray3DRef<double> blnum,
           FortranArray4DRef<double> vj0,
           FortranArray4DRef<double> vk0,
           FortranArray4DRef<double> vi0,
           FortranArray3DRef<double> bcj,
           FortranArray3DRef<double> bck,
           FortranArray3DRef<double> bci,
           int& nt,
           FortranArray1DRef<double> sumn,
           FortranArray1DRef<int> negn,
           FortranArray4DRef<double> ux,
           FortranArray4DRef<double> xib2,
           FortranArray3DRef<double> cmuv,
           FortranArray3DRef<double> volj0,
           FortranArray3DRef<double> volk0,
           FortranArray3DRef<double> voli0,
           FortranArray1DRef<int> nou,
           FortranArray2DRef<char[120]> bou,
           int& nbuf,
           int& ibufdim,
           int& myid,
           FortranArray1DRef<int> mblk2nd,
           int& maxbl,
           int& maxseg,
           FortranArray1DRef<int> nbci0,
           FortranArray1DRef<int> nbcj0,
           FortranArray1DRef<int> nbck0,
           FortranArray1DRef<int> nbcidim,
           FortranArray1DRef<int> nbcjdim,
           FortranArray1DRef<int> nbckdim,
           FortranArray4DRef<int> ibcinfo,
           FortranArray4DRef<int> jbcinfo,
           FortranArray4DRef<int> kbcinfo,
           FortranArray1DRef<double> vormax,
           FortranArray1DRef<int> ivmax,
           FortranArray1DRef<int> jvmax,
           FortranArray1DRef<int> kvmax,
           FortranArray1DRef<int> idefrm,
           FortranArray1DRef<int> iadvance,
           FortranArray4DRef<double> qavg,
           int& nummem);

} // namespace resid_ns
