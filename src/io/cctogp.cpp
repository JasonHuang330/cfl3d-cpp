// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "cctogp.h"

namespace cctogp_ns {

void cctogp(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3,
            int& j1, int& j2, int& j3, int& k1, int& k2, int& k3,
            FortranArray4DRef<double> dum,
            FortranArray4DRef<double> dumi0,
            FortranArray4DRef<double> dumj0,
            FortranArray4DRef<double> dumk0,
            int& jdw, int& kdw, int& idw,
            FortranArray4DRef<double> dumgp,
            int& ldim)
{
    int jdim1, kdim1, idim1;
    int iw, kw, jw;
    int i, k, j, l;
    int id, id1, jd, kd, m;

    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;

    // determine values at grid points

    iw = 0;
    for (i = i1; (i3 > 0) ? (i <= i2) : (i >= i2); i += i3) {
        iw = iw + 1;
        kw = 0;
        for (k = k1; (k3 > 0) ? (k <= k2) : (k >= k2); k += k3) {
            kw = kw + 1;
            jw = 0;
            for (j = j1; (j3 > 0) ? (j <= j2) : (j >= j2); j += j3) {
                jw = jw + 1;
                for (l = 1; l <= ldim; l++) {
                    id  = i;
                    id1 = id - 1;
                    if (id1 <= 0)    id1 = 1;
                    if (id > idim1)  id  = idim1;
                    jd = j;
                    kd = k;

                    if (k == 1 || k == kdim) {
                        // k=1/k=kdim faces
                        if (k == kdim) kd = kdim1;
                        if (j == 1 || j == jdim) {
                            // edge points and corner points
                            if (j == jdim) jd = jdim1;
                            dumgp(jw, kw, iw, l) = 0.5 * (dum(jd, kd, id, l)
                                                         + dum(jd, kd, id1, l));
                        } else {
                            // interior points on k=1/k=kdim faces
                            m = 2;
                            if (k == kdim) m = 4;
                            dumgp(jw, kw, iw, l) = .25 * (dumk0(jd,   id,  l, m)
                                                         + dumk0(jd-1, id,  l, m)
                                                         + dumk0(jd,   id1, l, m)
                                                         + dumk0(jd-1, id1, l, m));
                        }
                    } else if (j == 1 || j == jdim) {
                        // interior points on j=1/j=jdim faces
                        if (j == jdim) jd = jdim1;
                        m = 2;
                        if (j == jdim) m = 4;
                        dumgp(jw, kw, iw, l) = .25 * (dumj0(kd,   id,  l, m)
                                                     + dumj0(kd-1, id,  l, m)
                                                     + dumj0(kd,   id1, l, m)
                                                     + dumj0(kd-1, id1, l, m));
                    } else {
                        // interior points
                        if (i == 1 || i == idim) continue;
                        dumgp(jw, kw, iw, l) = .125 * (dum(j,   k,   i,   l)
                                                      + dum(j,   k,   i-1, l)
                                                      + dum(j-1, k,   i,   l)
                                                      + dum(j-1, k,   i-1, l)
                                                      + dum(j,   k-1, i,   l)
                                                      + dum(j,   k-1, i-1, l)
                                                      + dum(j-1, k-1, i,   l)
                                                      + dum(j-1, k-1, i-1, l));
                    }
                } // l
            } // j
        } // k
    } // i

    // interior points on i=1/i=idim faces
    iw = 0;
    for (i = i1; (i3 > 0) ? (i <= i2) : (i >= i2); i += i3) {
        iw = iw + 1;
        if (i != 1 && i != idim) continue;
        m  = 2;
        if (i == idim) m = 4;
        kw = 0;
        for (k = k1; (k3 > 0) ? (k <= k2) : (k >= k2); k += k3) {
            kw = kw + 1;
            jw = 0;
            for (j = j1; (j3 > 0) ? (j <= j2) : (j >= j2); j += j3) {
                jw = jw + 1;
                if (j == 1 || j == jdim || k == 1 || k == kdim) continue;
                // interior points
                for (l = 1; l <= ldim; l++) {
                    dumgp(jw, kw, iw, l) = .25 * (dumi0(j,   k,   l, m)
                                                 + dumi0(j,   k-1, l, m)
                                                 + dumi0(j-1, k,   l, m)
                                                 + dumi0(j-1, k-1, l, m));
                } // l
            } // j
        } // k
    } // i

    return;
}

} // namespace cctogp_ns
