// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "ld_datk.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cstring>
#include <algorithm>

namespace ld_datk_ns {

void ld_datk(FortranArray4DRef<double> data, int& jdim, int& kdim, int& idim,
             FortranArray4DRef<double> datak, int& ldim,
             FortranArray1DRef<int> mp, int& np,
             int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& myid)
{
    if (np > 4) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "stopping...increase dimension of mp array in ld_datk and its calling routine");
        int ierrflg = -1;
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    // load data values into datak array:
    for (int i = ista; i <= iend; i++) {
        for (int j = jsta; j <= jend; j++) {
            for (int m = 1; m <= np; m++) {
                for (int l = 1; l <= ldim; l++) {
                    datak(j, m, i, l) = data(j, mp(m), i, l);
                }
            }
        }
    }
}

} // namespace ld_datk_ns
