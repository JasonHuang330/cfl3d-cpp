// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "getdhdr.h"
#include <cstring>

namespace getdhdr_ns {

// Helper: assign a 10-char Fortran string literal to a char[10] element
// (no null terminator — fixed-length character semantics)
static inline void assign10(char (&dst)[10], const char* src) {
    std::memcpy(dst, src, 10);
}

void getdhdr(FortranArray1DRef<char[10]> datahdr, int& ibctyp, int& iver, int& ndata)
{
    // Initialize all 10 entries to 10 spaces
    for (int m = 1; m <= 10; m++) {
        assign10(datahdr(m), "          ");
    }

    if (ibctyp == 2002) {
        assign10(datahdr(1), "    p/pinf");
    }

    if (ibctyp == 2003) {
        assign10(datahdr(1), "      mach");
        assign10(datahdr(2), "   pt/pinf");
        assign10(datahdr(3), "   tt/tinf");
        assign10(datahdr(4), "     alpha");
        assign10(datahdr(5), "      beta");
        if (iver >= 5) {
            if (std::abs(ndata) > 5) {
                assign10(datahdr(6), "     turb1");
            }
            if (std::abs(ndata) > 6) {
                assign10(datahdr(7), "     turb2");
            }
        }
    }

    if (iver >= 5) {
        if (ibctyp == 2009 || ibctyp == 2010) {
            assign10(datahdr(1), "   pt/pinf");
            assign10(datahdr(2), "   tt/tinf");
            assign10(datahdr(3), "     alpha");
            assign10(datahdr(4), "      beta");
            if (std::abs(ndata) > 4) {
                assign10(datahdr(5), "     turb1");
            }
            if (std::abs(ndata) > 5) {
                assign10(datahdr(6), "     turb2");
            }
        }
    }

    if (std::abs(ibctyp) == 2004 || std::abs(ibctyp) == 2014) {
        if (iver == 4) {
            assign10(datahdr(1), "        cq");
        } else {
            assign10(datahdr(1), "   tw/tinf");
            assign10(datahdr(2), "        cq");
        }
    }

    if (std::abs(ibctyp) == 2016) {
        assign10(datahdr(1), "   tw/tinf");
        assign10(datahdr(2), "        cq");
        assign10(datahdr(3), "       cqu");
        assign10(datahdr(4), "     sjetx");
        assign10(datahdr(5), "     sjety");
        assign10(datahdr(6), "     sjetz");
        assign10(datahdr(7), "     rfreq");
    }

    if (ibctyp == 2005) {
        assign10(datahdr(1), "      nblp");
        assign10(datahdr(2), "     dthtx");
        assign10(datahdr(3), "     dthty");
        assign10(datahdr(4), "     dthtz");
    }

    if (ibctyp == 2006) {
        assign10(datahdr(1), "      nblc");
        assign10(datahdr(2), "    p/pinf");
        assign10(datahdr(3), " integ dir");
        assign10(datahdr(4), " axial dir");
    }

    if (iver >= 5) {
        if (ibctyp == 2007) {
            assign10(datahdr(1), "       rho");
            assign10(datahdr(2), "         u");
            assign10(datahdr(3), "         v");
            assign10(datahdr(4), "         w");
            assign10(datahdr(5), "         p");
            if (std::abs(ndata) > 5) {
                assign10(datahdr(6), "     turb1");
            }
            if (std::abs(ndata) > 6) {
                assign10(datahdr(7), "     turb2");
            }
        }
    }

    if (ibctyp == 2102) {
        assign10(datahdr(1), "    p/pinf");
        assign10(datahdr(2), "   dp/pinf");
        assign10(datahdr(3), "    rfreqp");
        assign10(datahdr(4), "      lref");
    }

    if (iver >= 5) {
        if (ibctyp == 2008) {
            assign10(datahdr(1), "       rho");
            assign10(datahdr(2), "         u");
            assign10(datahdr(3), "         v");
            assign10(datahdr(4), "         w");
            if (std::abs(ndata) > 4) {
                assign10(datahdr(5), "     turb1");
            }
            if (std::abs(ndata) > 5) {
                assign10(datahdr(6), "     turb2");
            }
        }

        if (ibctyp == 2018) {
            assign10(datahdr(1), "         T");
            assign10(datahdr(2), "      rhou");
            assign10(datahdr(3), "      rhov");
            assign10(datahdr(4), "      rhow");
            if (std::abs(ndata) > 4) {
                assign10(datahdr(5), "     turb1");
            }
            if (std::abs(ndata) > 5) {
                assign10(datahdr(6), "     turb2");
            }
        }

        if (ibctyp == 2028) {
            assign10(datahdr(1), "      freq");
            assign10(datahdr(2), "   rhoumax");
            assign10(datahdr(3), "   rhovmax");
            assign10(datahdr(4), "   rhowmax");
            if (std::abs(ndata) > 4) {
                assign10(datahdr(5), "     turb1");
            }
            if (std::abs(ndata) > 5) {
                assign10(datahdr(6), "     turb2");
            }
        }
    }

    return;
}

} // namespace getdhdr_ns
