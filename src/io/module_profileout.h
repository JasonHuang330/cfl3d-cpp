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

namespace module_profileout_ns {

struct output_parm {
    char name[40];
    int iblk;
    int kend;
    int kstr;
    int jend;
    int jstr;
    int iend;
    int istr;
    int wall_dir;
    int wall_end;
    int wall_str;
};

// module arrays (owning storage)
inline FortranArray1D<double> oparm;
inline FortranArray1D<double> cfparm;

void profile_read_input(char (&filename)[1]);
void profile_plot(int& jdim, int& kdim, int& idim, int& nummem, int& nblk, FortranArray4DRef<double> q, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> ux, FortranArray4DRef<double> turre, FortranArray3DRef<double> smin);
void cfcp_plot(int& jdim, int& kdim, int& idim, int& nummem, int& nblk, FortranArray4DRef<double> q, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> ux, FortranArray4DRef<double> turre, FortranArray3DRef<double> smin);

} // namespace module_profileout_ns
