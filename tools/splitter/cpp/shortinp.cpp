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

#include "shortinp.h"

namespace shortinp_ns {

// converts cfl3d "shortcut" input to full range of index values
void shortinp(int& nb1, int& ne1, int& nb2, int& ne2, int& idm, int& jdm, int& kdm, int& iface)
{
    if (nb1 == 0) nb1 = 1;
    if (nb2 == 0) nb2 = 1;
    if (ne1 == 0) {
        if (iface == 1 || iface == 2) ne1 = jdm;
        if (iface == 3 || iface == 4) ne1 = kdm;
        if (iface == 5 || iface == 6) ne1 = idm;
    }
    if (ne2 == 0) {
        if (iface == 1 || iface == 2) ne2 = kdm;
        if (iface == 3 || iface == 4) ne2 = idm;
        if (iface == 5 || iface == 6) ne2 = jdm;
    }

    return;
}

} // namespace shortinp_ns
