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

namespace rotateqb_ns {

void rotate(int& jdim, int& kdim, int& idim, FortranArray2DRef<double> t, FortranArray3DRef<double> tti, FortranArray3DRef<double> ttj, FortranArray3DRef<double> ttk, FortranArray2DRef<double> t1, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& nbl, int& irot, double& rfreqr, double& omegx, double& omegy, double& omegz, double& xorg, double& yorg, double& zorg, double& thetax, double& thetay, double& thetaz, double& thxold, double& thyold, double& thzold, int& iupdat, double& time2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim);
void rotateq(int& jd, int& kd, int& id, FortranArray4DRef<double> q, FortranArray4DRef<double> qrot, int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend, double& dthtx, double& dthty, double& dthtz);
void rotateqb(int& nbl, double& dthtx, double& dthty, double& dthtz, int& maxbl, int& iitot, FortranArray1DRef<int> ibcg, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, FortranArray3DRef<double> qb);

} // namespace rotateqb_ns
