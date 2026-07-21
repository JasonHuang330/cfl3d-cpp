// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "q8vrev.h"

namespace q8vrev_ns {

void q8vrev(int& n, FortranArray1DRef<double> a, int& ndum, FortranArray1DRef<double> b) {
    int np = n + 1;
    for (int i = 1; i <= n; i++) {
        b(i) = a(np - i);
    }
}

} // namespace q8vrev_ns
