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

namespace patcher_ns {

void patcher(int& nbl, FortranArray2DRef<int> lw, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& ncall, int& ioutpt, int& it_thro, int& maxbl, int& msub1, int& intmx, int& mxxe, int& mptch, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<double> windx, int& nintr, FortranArray2DRef<int> iindx, FortranArray1DRef<int> llimit, FortranArray1DRef<int> iitmax, FortranArray1DRef<int> mmcxie, FortranArray1DRef<int> mmceta, FortranArray1DRef<int> ncheck, FortranArray1DRef<int> iifit, FortranArray1DRef<int> mblkpt, FortranArray1DRef<int> iic0, FortranArray1DRef<int> iiorph, FortranArray1DRef<int> iitoss, FortranArray1DRef<int> ifiner, FortranArray2DRef<double> factjlo, FortranArray2DRef<double> factjhi, FortranArray2DRef<double> factklo, FortranArray2DRef<double> factkhi, FortranArray2DRef<double> dx, FortranArray2DRef<double> dy, FortranArray2DRef<double> dz, FortranArray2DRef<double> dthetx, FortranArray2DRef<double> dthety, FortranArray2DRef<double> dthetz, FortranArray2DRef<int> isav_dpat, FortranArray3DRef<int> isav_dpat_b, FortranArray3DRef<double> xte, FortranArray3DRef<double> yte, FortranArray3DRef<double> zte, FortranArray3DRef<double> xmi, FortranArray3DRef<double> ymi, FortranArray3DRef<double> zmi, FortranArray3DRef<double> xmie, FortranArray3DRef<double> ymie, FortranArray3DRef<double> zmie, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<int> jte, FortranArray1DRef<int> kte, FortranArray3DRef<double> sxie, FortranArray3DRef<double> seta, FortranArray2DRef<double> sxie2, FortranArray2DRef<double> seta2, FortranArray2DRef<double> xie2s, FortranArray2DRef<double> eta2s, FortranArray1DRef<double> temp, FortranArray2DRef<double> x2, FortranArray2DRef<double> y2, FortranArray2DRef<double> z2, FortranArray1DRef<int> nblk1, FortranArray1DRef<int> nblk2, FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm, FortranArray2DRef<double> x1, FortranArray2DRef<double> y1, FortranArray2DRef<double> z1, FortranArray1DRef<int> lout, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, FortranArray1DRef<int> ireq_ar, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim);

} // namespace patcher_ns
