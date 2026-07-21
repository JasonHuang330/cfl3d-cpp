// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "q8smax.h"

namespace q8smax_ns {

double q8smax(int& n, FortranArray1DRef<double> a)
{
    double result = a(1);
    for (int i = 2; i <= n; i++) {
        if ((float)result < (float)a(i)) result = a(i);
    }
    return result;
}

} // namespace q8smax_ns
