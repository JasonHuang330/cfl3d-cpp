// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
//  ---------------------------------------------------------------------------
//  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
//  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
//  with point-matched, patched, overset, or embedded connectivities. Both
//  multigrid and mesh sequencing are available in time-accurate or
//  steady-state modes.
//
//  Copyright 2001 United States Government as represented by the Administrator
//  of the National Aeronautics and Space Administration. All Rights Reserved.
//
//  The CFL3D platform is licensed under the Apache License, Version 2.0
//  (the "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0.
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
//  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
//  License for the specific language governing permissions and limitations
//  under the License.
//  ---------------------------------------------------------------------------

#include "wkstn.h"

namespace wkstn_ns {

// Purpose: To find the first occurance in the array x that is equal
// to target; a replacement for the standard Cray function
int isrcheq(int& n, FortranArray1DRef<double> x, int& incx, double& target)
{
    int isrcheq_result = 0;
    int j = 1;
    if (n <= 0) return isrcheq_result;
    if (incx < 0) j = 1 - (n - 1) * incx;
    int i = 1;
    for (i = 1; i <= n; i++) {
        if ((float)x(j) == (float)target) goto label200;
        j = j + incx;
    }
label200:
    isrcheq_result = i;
    return isrcheq_result;
}

} // namespace wkstn_ns
