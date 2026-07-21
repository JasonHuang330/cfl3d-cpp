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

#include "mreal.h"

namespace mreal_ns {

// Purpose: Unload grid data from transfer array xyz to x,y,z work arrays
void mreal(FortranArray3DRef<double> xyz, int& mdim1, int& ndim1, int& mdim2, int& ndim2,
           FortranArray2DRef<double> x, FortranArray2DRef<double> y, FortranArray2DRef<double> z)
{
    // initialize the individual x,y,z arrays, since they are typically
    // larger that the corresponding sections in xyz (i.e. mdim2*ndim2
    // vs. mdim1*ndim1)
    for (int m = 1; m <= mdim2; m++) {
        for (int n = 1; n <= ndim2; n++) {
            x(m, n) = 0.;
            y(m, n) = 0.;
            z(m, n) = 0.;
        }
    }

    for (int m = 1; m <= mdim1; m++) {
        for (int n = 1; n <= ndim1; n++) {
            x(m, n) = xyz(m, n, 1);
            y(m, n) = xyz(m, n, 2);
            z(m, n) = xyz(m, n, 3);
        }
    }

    return;
}

} // namespace mreal_ns
