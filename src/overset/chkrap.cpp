// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "chkrap.h"
#include "termn8.h"
#include <algorithm>
#include <cstdio>

namespace chkrap_ns {

void chkrap(int& nbl, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& nface, int& idim, int& jdim, int& kdim, int& nn, int& maxxseg, int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend, int& iwrap, int& jwrap, int& kwrap, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid)
{
    int iokrap, nseg;
    int ista1, iend1, jsta1, jend1, ksta1, kend1;
    int jbc, ibc, kbc;

    iokrap = 0;
    iwrap = 0;
    jwrap = 0;
    kwrap = 0;

    // ***********************************************************************
    // k=1 or k=kdim is a singular axis
    // ***********************************************************************

    if (nface == 5 || nface == 6) {

        // check for wrap around in j direction, adjacent to k-axis,
        // within same i-range as the singular axis segment

        jwrap = 1;
        for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
            ista1 = jbcinfo(nbl, nseg, 2, 1);
            iend1 = jbcinfo(nbl, nseg, 3, 1);
            ksta1 = jbcinfo(nbl, nseg, 4, 1);
            kend1 = jbcinfo(nbl, nseg, 5, 1);
            if ((nface == 5 && ksta1 == 1) ||
                (nface == 6 && kend1 == kdim)) {
                if (ista1 < iend && iend1 > ista) {
                    // this segment on the j=1 face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    jbc = jbcinfo(nbl, nseg, 1, 1);
                    if (jbc != 0) jwrap = 0;
                }
            }
        }

        // check for wrap around in j direction, adjacent to k-axis,
        // within same i-range as the singular axis segment

        for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
            ista1 = jbcinfo(nbl, nseg, 2, 2);
            iend1 = jbcinfo(nbl, nseg, 3, 2);
            ksta1 = jbcinfo(nbl, nseg, 4, 2);
            kend1 = jbcinfo(nbl, nseg, 5, 2);
            if ((nface == 5 && ksta1 == 1) ||
                (nface == 6 && kend1 == kdim)) {
                if (ista1 < iend && iend1 > ista) {
                    // this segment on the j=jdim face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    jbc = jbcinfo(nbl, nseg, 1, 2);
                    if (jbc != 0) jwrap = 0;
                }
            }
        }

        // check for wrap around in i direction, adjacent to k-axis,
        // within same j-range as the singular axis segment

        iwrap = 1;
        for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
            jsta1 = ibcinfo(nbl, nseg, 2, 1);
            jend1 = ibcinfo(nbl, nseg, 3, 1);
            ksta1 = ibcinfo(nbl, nseg, 4, 1);
            kend1 = ibcinfo(nbl, nseg, 5, 1);
            if ((nface == 5 && ksta1 == 1) ||
                (nface == 6 && kend1 == kdim)) {
                if (jsta1 < jend && jend1 > jsta) {
                    // this segment on the i=1 face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    ibc = ibcinfo(nbl, nseg, 1, 1);
                    if (ibc != 0) iwrap = 0;
                }
            }
        }

        // check for wrap around in i direction, adjacent to k-axis,
        // within same j-range as the singular axis segment

        for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
            jsta1 = ibcinfo(nbl, nseg, 2, 2);
            jend1 = ibcinfo(nbl, nseg, 3, 2);
            ksta1 = ibcinfo(nbl, nseg, 4, 2);
            kend1 = ibcinfo(nbl, nseg, 5, 2);
            if ((nface == 5 && ksta1 == 1) ||
                (nface == 6 && kend1 == kdim)) {
                if (jsta1 < jend && jend1 > jsta) {
                    // this segment on the i=idim face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    ibc = ibcinfo(nbl, nseg, 1, 2);
                    if (ibc != 0) iwrap = 0;
                }
            }
        }

        // MUST have proper wrap around set in either i or j
        // directions in order to use bctype 1012 on k=1 or k=kdim

        if (jwrap > 0 || iwrap > 0) {
            iokrap = 1;
        } else {
            if (nface == 5) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, " error in applying bctype 1012 on k=1 boundary:");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, "    must have wrap around in i or j directions");
            } else {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, " error in applying bctype 1012 on k=kdim boundary:");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, "    must have wrap around in i or j directions");
            }
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }

    }

    // ***********************************************************************
    // j=1 or j=jdim is a singular axis
    // ***********************************************************************

    if (nface == 3 || nface == 4) {

        // check for wrap around in k direction, adjacent to j-axis,
        // within same i-range as the singular axis segment

        kwrap = 1;
        for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
            ista1 = kbcinfo(nbl, nseg, 2, 1);
            iend1 = kbcinfo(nbl, nseg, 3, 1);
            jsta1 = kbcinfo(nbl, nseg, 4, 1);
            jend1 = kbcinfo(nbl, nseg, 5, 1);
            if ((nface == 3 && jsta1 == 1) ||
                (nface == 4 && jend1 == jdim)) {
                if (ista1 < iend && iend1 > ista) {
                    // this segment on the k=1 face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    kbc = kbcinfo(nbl, nseg, 1, 1);
                    if (kbc != 0) kwrap = 0;
                }
            }
        }

        // check for wrap around in k direction, adjacent to j-axis,
        // within same i-range as the singular axis segment

        for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
            ista1 = kbcinfo(nbl, nseg, 2, 2);
            iend1 = kbcinfo(nbl, nseg, 3, 2);
            jsta1 = kbcinfo(nbl, nseg, 4, 2);
            jend1 = kbcinfo(nbl, nseg, 5, 2);
            if ((nface == 3 && jsta1 == 1) ||
                (nface == 4 && jend1 == jdim)) {
                if (ista1 < iend && iend1 > ista) {
                    // this segment on the k=1 face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    kbc = kbcinfo(nbl, nseg, 1, 2);
                    if (kbc != 0) kwrap = 0;
                }
            }
        }

        // check for wrap around in i direction, adjacent to j-axis,
        // within same k-range as the singular axis segment

        iwrap = 1;
        for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
            jsta1 = ibcinfo(nbl, nseg, 2, 1);
            jend1 = ibcinfo(nbl, nseg, 3, 1);
            ksta1 = ibcinfo(nbl, nseg, 4, 1);
            kend1 = ibcinfo(nbl, nseg, 5, 1);
            if ((nface == 3 && jsta1 == 1) ||
                (nface == 4 && jend1 == jdim)) {
                if (ksta1 < kend && kend1 > ksta) {
                    // this segment on the i=1 face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    ibc = ibcinfo(nbl, nseg, 1, 1);
                    if (ibc != 0) iwrap = 0;
                }
            }
        }

        // check for wrap around in i direction, adjacent to j-axis,
        // within same k-range as the singular axis segment

        for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
            jsta1 = ibcinfo(nbl, nseg, 2, 2);
            jend1 = ibcinfo(nbl, nseg, 3, 2);
            ksta1 = ibcinfo(nbl, nseg, 4, 2);
            kend1 = ibcinfo(nbl, nseg, 5, 2);
            if ((nface == 3 && jsta1 == 1) ||
                (nface == 4 && jend1 == jdim)) {
                if (ksta1 < kend && kend1 > ksta) {
                    // this segment on the i=1 face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    ibc = ibcinfo(nbl, nseg, 1, 2);
                    if (ibc != 0) iwrap = 0;
                }
            }
        }

        // MUST have proper wrap around set in either i or k
        // directions in order to use bctype 1012 on j=1 or j=jdim

        if (kwrap > 0 || iwrap > 0) {
            iokrap = 1;
        } else {
            if (nface == 3) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, " error in applying bctype 1012 on j=1 boundary:");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, "    must have wrap around in i or k directions");
            } else {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, " error in applying bctype 1012 on j=jdim boundary:");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, "    must have wrap around in i or k directions");
            }
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }

    }

    // ***********************************************************************
    // i=1 or i=idim is a singular axis
    // ***********************************************************************

    if (nface == 1 || nface == 2) {

        // check for wrap around in j direction, adjacent to i-axis,
        // within same k-range as the singular axis segment

        jwrap = 1;
        for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
            ista1 = jbcinfo(nbl, nseg, 2, 1);
            iend1 = jbcinfo(nbl, nseg, 3, 1);
            ksta1 = jbcinfo(nbl, nseg, 4, 1);
            kend1 = jbcinfo(nbl, nseg, 5, 1);
            if ((nface == 1 && ista1 == 1) ||
                (nface == 2 && iend1 == idim)) {
                if (ksta1 < kend && kend1 > ksta) {
                    // this segment on the j=1 face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    jbc = jbcinfo(nbl, nseg, 1, 1);
                    if (jbc != 0) jwrap = 0;
                }
            }
        }

        // check for wrap around in j direction, adjacent to i-axis,
        // within same k-range as the singular axis segment

        for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
            ista1 = jbcinfo(nbl, nseg, 2, 2);
            iend1 = jbcinfo(nbl, nseg, 3, 2);
            ksta1 = jbcinfo(nbl, nseg, 4, 2);
            kend1 = jbcinfo(nbl, nseg, 5, 2);
            if ((nface == 1 && ista1 == 1) ||
                (nface == 2 && iend1 == idim)) {
                if (ksta1 < kend && kend1 > ksta) {
                    // this segment on the j=jdim face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    jbc = jbcinfo(nbl, nseg, 1, 2);
                    if (jbc != 0) jwrap = 0;
                }
            }
        }

        // check for wrap around in k direction, adjacent to i-axis,
        // within same j-range as the singular axis segment

        kwrap = 1;
        for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
            ista1 = kbcinfo(nbl, nseg, 2, 1);
            iend1 = kbcinfo(nbl, nseg, 3, 1);
            jsta1 = kbcinfo(nbl, nseg, 4, 1);
            jend1 = kbcinfo(nbl, nseg, 5, 1);
            if ((nface == 1 && ista1 == 1) ||
                (nface == 2 && iend1 == idim)) {
                if (jsta1 < jend && jend1 > jsta) {
                    // this segment on the k=1 face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    kbc = kbcinfo(nbl, nseg, 1, 1);
                    if (kbc != 0) kwrap = 0;
                }
            }
        }

        // check for wrap around in k direction, adjacent to i-axis,
        // within same j-range as the singular axis segment

        for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
            ista1 = kbcinfo(nbl, nseg, 2, 2);
            iend1 = kbcinfo(nbl, nseg, 3, 2);
            jsta1 = kbcinfo(nbl, nseg, 4, 2);
            jend1 = kbcinfo(nbl, nseg, 5, 2);
            if ((nface == 1 && ista1 == 1) ||
                (nface == 2 && iend1 == idim)) {
                if (jsta1 < jend && jend1 > jsta) {
                    // this segment on the k=kdim face must have a wrap around in
                    // order to apply bctype 1012 - turn OFF wraparound flag if not
                    kbc = kbcinfo(nbl, nseg, 1, 2);
                    if (kbc != 0) kwrap = 0;
                }
            }
        }

        // MUST have proper wrap around set in either k or j
        // directions in order to use bctype 1012 on i=1 or i=idim

        if (jwrap > 0 || kwrap > 0) {
            iokrap = 1;
        } else {
            if (nface == 1) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, " error in applying bctype 1012 on i=1 boundary:");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, "    must have wrap around in k or j directions");
            } else {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, " error in applying bctype 1012 on i=idim boundary:");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, "    must have wrap around in k or j directions");
            }
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }

    }

    return;
}

} // namespace chkrap_ns
