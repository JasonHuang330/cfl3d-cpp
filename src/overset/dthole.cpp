// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dthole.h"
#include "ccomplex.h"

namespace dthole_ns {

void dthole(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> dtj, FortranArray3DRef<double> vol, FortranArray3DRef<double> blank, double& dtmin, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
//
//     nou(1) = min(nou(1)+1,ibufdim)
//     if (isklton.gt.0) then
//        nou(1) = min(nou(1)+1,ibufdim)
//        write(bou(nou(1),1),*)' in dthole, dtmin= ',real(dtmin)
//     end if
//
    int kn = jdim * kdim * (idim - 1);
//cdir$ ivdep
    for (int n = 1; n <= kn; n++) {
        bool cond = ((float)blank(n, 1, 1) == 0.0f);
        double a = vol(n, 1, 1) / dtmin;
        double b = dtj(n, 1, 1);
        dtj(n, 1, 1) = ccomplex_ns::ccvmgt(a, b, cond);
    }
    return;
}

} // namespace dthole_ns
