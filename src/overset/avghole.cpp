// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "avghole.h"
#include "xupdt.h"
#include <algorithm>

namespace avghole_ns {

void avghole(FortranArray4DRef<double> q, FortranArray3DRef<double> blank, int& jdim, int& kdim, int& idim, int& nbl, int& ldim, int& int_updt)
{
    for (int l = 1; l <= ldim; l++) {
        for (int i = 1; i <= idim - 1; i++) {
            for (int k = 1; k <= kdim - 1; k++) {
                for (int j = 1; j <= jdim - 1; j++) {
                    if (blank(j, k, i) == 0.) {
                        int jp = std::min(j + 1, jdim - 1);
                        int jm = std::max(j - 1, 1);
                        int kp = std::min(k + 1, kdim - 1);
                        int km = std::max(k - 1, 1);
                        int ip = std::min(i + 1, idim - 1);
                        int im = std::max(i - 1, 1);
                        q(j, k, i, l) = (q(jp,  k, ip, l) + q(jm,  k, ip, l)
                                      +  q( j, kp, ip, l) + q( j, km, ip, l)
                                      +  q(jp,  k, im, l) + q(jm,  k, im, l)
                                      +  q( j, kp, im, l) + q( j, km, im, l)) * 0.125;
                    }
                }
            }
        }
    }
}

void xupdt(FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, int& jdim, int& kdim, int& idim, int& nbl, int& ldim, FortranArray3DRef<double> bcj0, FortranArray3DRef<double> bck0, FortranArray3DRef<double> bci0, int& maxbl, int& iitot, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg, FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg, FortranArray3DRef<double> qb, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& int_updt)
{
    xupdt_ns::xupdt(q, qj0, qk0, qi0, jdim, kdim, idim, nbl, ldim, bcj0, bck0, bci0, maxbl, iitot, iibg, kkbg, jjbg, ibcg, lbg, ibpntsg, qb, nou, bou, nbuf, ibufdim, int_updt);
}

} // namespace avghole_ns
