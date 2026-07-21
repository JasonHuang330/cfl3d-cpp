// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "getdhdr.h"
#include <cstring>

namespace getdhdr_ns {

void getdhdr(FortranArray1DRef<char[10]> datahdr, int& ibctyp, int& ndata)
{
    int m;

    // do 1 m=1,10
    // 1 datahdr(m) = '          '
    for (m = 1; m <= 10; m++) {
        std::memcpy(datahdr(m), "          ", 10);
    }

    if (ibctyp == 2002) {
        std::memcpy(datahdr(1), "    p/pinf", 10);
    }

    if (ibctyp == 2003) {
        std::memcpy(datahdr(1), "      mach", 10);
        std::memcpy(datahdr(2), "   pt/pinf", 10);
        std::memcpy(datahdr(3), "   tt/tinf", 10);
        std::memcpy(datahdr(4), "     alpha", 10);
        std::memcpy(datahdr(5), "      beta", 10);
        if (std::abs(ndata) > 5) {
            std::memcpy(datahdr(6), "     turb1", 10);
        }
        if (std::abs(ndata) > 6) {
            std::memcpy(datahdr(7), "     turb2", 10);
        }
    }

    if (ibctyp == 2009 || ibctyp == 2010) {
        std::memcpy(datahdr(1), "   pt/pinf", 10);
        std::memcpy(datahdr(2), "   tt/tinf", 10);
        std::memcpy(datahdr(3), "     alpha", 10);
        std::memcpy(datahdr(4), "      beta", 10);
        if (std::abs(ndata) > 4) {
            std::memcpy(datahdr(5), "     turb1", 10);
        }
        if (std::abs(ndata) > 5) {
            std::memcpy(datahdr(6), "     turb2", 10);
        }
    }

    if (ibctyp == 2019) {
        std::memcpy(datahdr(1), "  pt/ptinf", 10);
        std::memcpy(datahdr(2), "  tt/ttinf", 10);
        if (std::abs(ndata) > 2) {
            std::memcpy(datahdr(3), "     turb1", 10);
        }
        if (std::abs(ndata) > 3) {
            std::memcpy(datahdr(4), "     turb2", 10);
        }
    }

    if (std::abs(ibctyp) == 2004 || std::abs(ibctyp) == 2034) {
        std::memcpy(datahdr(1), "   tw/tinf", 10);
        std::memcpy(datahdr(2), "        cq", 10);
    }

    if (std::abs(ibctyp) == 2014) {
        std::memcpy(datahdr(1), "   tw/tinf", 10);
        std::memcpy(datahdr(2), "        cq", 10);
        std::memcpy(datahdr(3), "     index", 10);
    }

    if (std::abs(ibctyp) == 2024) {
        std::memcpy(datahdr(1), "   tw/tinf", 10);
        std::memcpy(datahdr(2), "        cq", 10);
        std::memcpy(datahdr(3), " intermtbc", 10);
    }

    if (ibctyp == 2005) {
        std::memcpy(datahdr(1), "      nblp", 10);
        std::memcpy(datahdr(2), "     dthtx", 10);
        std::memcpy(datahdr(3), "     dthty", 10);
        std::memcpy(datahdr(4), "     dthtz", 10);
    }

    if (ibctyp == 2006) {
        std::memcpy(datahdr(1), "      nblc", 10);
        std::memcpy(datahdr(2), "    p/pinf", 10);
        std::memcpy(datahdr(3), " integ dir", 10);
        std::memcpy(datahdr(4), " axial dir", 10);
    }

    if (ibctyp == 2007) {
        std::memcpy(datahdr(1), "       rho", 10);
        std::memcpy(datahdr(2), "         u", 10);
        std::memcpy(datahdr(3), "         v", 10);
        std::memcpy(datahdr(4), "         w", 10);
        std::memcpy(datahdr(5), "         p", 10);
        if (std::abs(ndata) > 5) {
            std::memcpy(datahdr(6), "     turb1", 10);
        }
        if (std::abs(ndata) > 6) {
            std::memcpy(datahdr(7), "     turb2", 10);
        }
    }

    if (ibctyp == 2102) {
        std::memcpy(datahdr(1), "    p/pinf", 10);
        std::memcpy(datahdr(2), "   dp/pinf", 10);
        std::memcpy(datahdr(3), "    rfreqp", 10);
        std::memcpy(datahdr(4), "      lref", 10);
    }

    if (ibctyp == 2103) {
        std::memcpy(datahdr(1), "    p/pinf", 10);
        std::memcpy(datahdr(2), "   dp/pinf", 10);
        std::memcpy(datahdr(3), "    rfreqp", 10);
        std::memcpy(datahdr(4), "      lref", 10);
        std::memcpy(datahdr(5), "    phioff", 10);
    }

    if (ibctyp == 2008 || ibctyp == 2038) {
        std::memcpy(datahdr(1), "       rho", 10);
        std::memcpy(datahdr(2), "         u", 10);
        std::memcpy(datahdr(3), "         v", 10);
        std::memcpy(datahdr(4), "         w", 10);
        if (std::abs(ndata) > 4) {
            std::memcpy(datahdr(5), "     turb1", 10);
        }
        if (std::abs(ndata) > 5) {
            std::memcpy(datahdr(6), "     turb2", 10);
        }
    }

    if (std::abs(ibctyp) == 2016) {
        std::memcpy(datahdr(1), "   tw/tinf", 10);
        std::memcpy(datahdr(2), "        cq", 10);
        std::memcpy(datahdr(3), "       cqu", 10);
        std::memcpy(datahdr(4), "     sjetx", 10);
        std::memcpy(datahdr(5), "     sjety", 10);
        std::memcpy(datahdr(6), "     sjetz", 10);
        std::memcpy(datahdr(7), "     rfreq", 10);
    }

    if (std::abs(ibctyp) == 2026) {
        std::memcpy(datahdr(1), "      vmag", 10);
        std::memcpy(datahdr(2), "     rfreq", 10);
        std::memcpy(datahdr(3), "  sideangj", 10);
        std::memcpy(datahdr(4), "       sxa", 10);
        std::memcpy(datahdr(5), "       sya", 10);
        std::memcpy(datahdr(6), "       sza", 10);
        std::memcpy(datahdr(7), "       sxb", 10);
        std::memcpy(datahdr(8), "       syb", 10);
        std::memcpy(datahdr(9), "       szb", 10);
        if (std::abs(ndata) > 9) {
            std::memcpy(datahdr(10), "     turb1", 10);
        }
        if (std::abs(ndata) > 10) {
            std::memcpy(datahdr(11), "     turb2", 10);
        }
    }

    if (ibctyp == 2018) {
        std::memcpy(datahdr(1), "         T", 10);
        std::memcpy(datahdr(2), "      rhou", 10);
        std::memcpy(datahdr(3), "      rhov", 10);
        std::memcpy(datahdr(4), "      rhow", 10);
        if (std::abs(ndata) > 4) {
            std::memcpy(datahdr(5), "     turb1", 10);
        }
        if (std::abs(ndata) > 5) {
            std::memcpy(datahdr(6), "     turb2", 10);
        }
    }

    if (ibctyp == 2028) {
        std::memcpy(datahdr(1), "      freq", 10);
        std::memcpy(datahdr(2), "   rhoumax", 10);
        std::memcpy(datahdr(3), "   rhovmax", 10);
        std::memcpy(datahdr(4), "   rhowmax", 10);
        if (std::abs(ndata) > 4) {
            std::memcpy(datahdr(5), "     turb1", 10);
        }
        if (std::abs(ndata) > 5) {
            std::memcpy(datahdr(6), "     turb2", 10);
        }
    }

    return;
}

} // namespace getdhdr_ns
