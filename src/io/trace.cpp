// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "trace.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cmath>

namespace trace_ns {

void trace(int& icall, int& idum1, int& idum2, int& idum3, int& idum4, double& dum1, double& dum2, double& dum3)
{
    int& itrace = cmn_tracer.itrace;

    int iunit = 7;

    if (itrace < 0) return;

    if (icall == 0) {
        fortran_write_unit(iunit,
            " %2s*** Search History During Generalized"
            " Coordinate Interpolation Number %3d ***\n\n",
            "", idum1);
    }

    if (icall == 1 || icall == 2) {
        if (itrace == 0) {
            FILE* f = fortran_get_unit(iunit);
            rewind(f);
            fortran_write_unit(iunit,
                " %2s*** Search History During Generalized"
                " Coordinate Interpolation Number %3d ***\n\n",
                "", idum1);
        }
    }

    if (icall == 1) {
        fortran_write_unit(iunit,
            "   iterations for \"to\" grid:  cell j,k = %3d,%3d"
            "  with center xc,yc,zc = %11.4E,%11.4E,%11.4E\n",
            idum2, idum3, (float)dum1, (float)dum2, (float)dum3);
        fortran_write_unit(iunit,
            "    using ifit = %3d\n",
            idum4);
    }

    if (icall == 2) {
        fortran_write_unit(iunit,
            "   iterations for \"to\" grid:  cell j,k=%3d,%3d"
            "  with edge mid-point xc,yc,zc=%11.4E,%11.4E,%11.4E\n",
            idum2, idum3, (float)dum1, (float)dum2, (float)dum3);
        fortran_write_unit(iunit,
            "    using ifit = %3d\n",
            idum4);
    }

    if (icall == 3) {
        fortran_write_unit(iunit,
            "     intern=%3d\n",
            idum1);
    }

    if (icall == 4) {
        fortran_write_unit(iunit,
            "         xie,eta=%10.3E,%10.3E (local values)\n",
            (float)dum1, (float)dum2);
    }

    if (icall == 5) {
        fortran_write_unit(iunit,
            "   could not find \"to\" cell center in \"from\" block%3d\n"
            "   will try another block\n",
            idum1);
    }

    if (icall == 6) {
        fortran_write_unit(iunit,
            "     checking to see if this point really lies in block %3d\n",
            idum1);
    }

    if (icall == 7) {
        fortran_write_unit(iunit,
            "    point actually lies in block %3d\n",
            idum1);
        fortran_write_unit(iunit,
            "    future searches to be done in this block\n");
    }

    if (icall == 8) {
        fortran_write_unit(iunit,
            "    incompatable orientation of \"from\" cell and \"to\" cell\n");
        fortran_write_unit(iunit,
            "      ...will try in a nearby cell\n");
    }

    if (icall == 9) {
        fortran_write_unit(iunit,
            "       using y-z to invert in \"from\" grid%3d  cell j,k=%3d,%3d\n",
            idum1, idum2, idum3);
    }

    if (icall == 10) {
        fortran_write_unit(iunit,
            "       using x-z to invert in \"from\" grid%3d  cell j,k=%3d,%3d\n",
            idum1, idum2, idum3);
    }

    if (icall == 11) {
        fortran_write_unit(iunit,
            "       using x-y to invert in \"from\" grid%3d  cell j,k=%3d,%3d\n",
            idum1, idum2, idum3);
    }

    if (icall == 41) {
        fortran_write_unit(iunit,
            "     search routine off track...using minimum"
            " distance search to get back on track\n");
    }

    if (icall == 42) {
        fortran_write_unit(iunit,
            "     frozen convergence...attempting to break"
            " cycle via minimum distance search\n");
    }

    if (icall == 50) {
        fortran_write_unit(iunit,
            " marking point j,k=%3d%3d as an orphan \n",
            idum1, idum2);
    }

    if (icall == 99) {
        if (idum1 != idum3 || idum2 != idum4) {
            fortran_write_unit(iunit,
                " branch cut modification: before, jp,kp=%3d%3d"
                " after, jpc,kpc=%3d%3d\n",
                idum1, idum2, idum3, idum4);
        }
    }

    return;
}

} // namespace trace_ns
