// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// ---------------------------------------------------------------------------
// CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
// Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
// with point-matched, patched, overset, or embedded connectivities. Both
// multigrid and mesh sequencing are available in time-accurate or
// steady-state modes.
//
// Copyright 2001 United States Government as represented by the Administrator
// of the National Aeronautics and Space Administration. All Rights Reserved.
//
// The CFL3D platform is licensed under the Apache License, Version 2.0
// (the "License"); you may not use this file except in compliance with the
// License. You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0.
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
// License for the specific language governing permissions and limitations
// under the License.
// ---------------------------------------------------------------------------

#include "my_flush.h"
#include "runtime/fortran_io.h"
#include <cstdio>

namespace my_flush_ns {

// Purpose: Generic routine to flush an output buffer
void my_flush(int& iunit) {
    // Flush the FILE* associated with the given Fortran unit number.
    // This covers all platform-specific flush variants from the original
    // Fortran (#ifdef IBM, SGI, SUN, HP, ALPHA, PG, LAHEY, GENERIC, etc.)
    FILE* fp = fortran_get_unit(iunit);
    if (fp != nullptr) {
        fflush(fp);
    }
}

} // namespace my_flush_ns
