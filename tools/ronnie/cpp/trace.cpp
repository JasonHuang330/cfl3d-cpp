// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// trace.cpp — faithful C++ translation of trace.F (CFL3D ronnie).
// Writes the search-routine history for the current "to" cell to unit 7.
//
// /tracer/ itrace -> tracer_.itrace.  In the point-matched test itrace<0 so
// this returns immediately, but the full formatted output is reproduced.
//
// Formatted output emulation:
//   i3   -> f_iw(v,3)          (right-justified width-3 integer)
//   e11.4-> tc::fmtE(v,11,4)   (Fortran 0.dddd E+xx mantissa form)
//   e10.3-> tc::fmtE(v,10,3)
// The leading blank of every FORMAT (carriage-control ' ') is preserved.
#include "ron_common.h"

void trace(int icall,int idum1,int idum2,int idum3,int idum4,
           double dum1,double dum2,double dum3)
{
    int iunit = 7;
    FILE* f = FUNIT(iunit);

    if (tracer_.itrace < 0) return;

    if (icall == 0) {
        // 100 format(' ',2x,'*** Search History ... Number ',i3,' ***',/)
        std::fprintf(f, "   *** Search History During Generalized"
                        " Coordinate Interpolation Number %s ***\n\n",
                     f_iw(idum1,3).c_str());
    }

    if (icall == 1 || icall == 2) {
        if (tracer_.itrace == 0) {
            std::rewind(f);
            std::fprintf(f, "   *** Search History During Generalized"
                            " Coordinate Interpolation Number %s ***\n\n",
                         f_iw(idum1,3).c_str());
        }
    }

    if (icall == 1) {
        // 101 format(' ',2x,'iterations for "to" grid:  cell j,k = ',i3,
        //   ',',i3,'  with center xc,yc,zc = ',e11.4,',',e11.4,',',e11.4)
        std::fprintf(f, "   iterations for \"to\" grid:  cell j,k = %s,%s"
                        "  with center xc,yc,zc = %s,%s,%s\n",
                     f_iw(idum2,3).c_str(), f_iw(idum3,3).c_str(),
                     tc::fmtE(dum1,11,4).c_str(), tc::fmtE(dum2,11,4).c_str(),
                     tc::fmtE(dum3,11,4).c_str());
        // 111 format(' ',3x,'using ifit = ',i3)
        std::fprintf(f, "    using ifit = %s\n", f_iw(idum4,3).c_str());
    }

    if (icall == 2) {
        // 102 format(' ',2x,'iterations for "to" grid:  cell j,k=',i3,
        //   ',',i3,'  with edge mid-point xc,yc,zc=',e11.4,',',e11.4,',',e11.4)
        std::fprintf(f, "   iterations for \"to\" grid:  cell j,k=%s,%s"
                        "  with edge mid-point xc,yc,zc=%s,%s,%s\n",
                     f_iw(idum2,3).c_str(), f_iw(idum3,3).c_str(),
                     tc::fmtE(dum1,11,4).c_str(), tc::fmtE(dum2,11,4).c_str(),
                     tc::fmtE(dum3,11,4).c_str());
        std::fprintf(f, "    using ifit = %s\n", f_iw(idum4,3).c_str());
    }

    if (icall == 3) {
        // 103 format(' ',4x,'intern=',i3)
        std::fprintf(f, "     intern=%s\n", f_iw(idum1,3).c_str());
    }

    if (icall == 4) {
        // 104 format(' ',8x,'xie,eta=',e10.3,',',e10.3,' (local values)')
        std::fprintf(f, "         xie,eta=%s,%s (local values)\n",
                     tc::fmtE(dum1,10,3).c_str(), tc::fmtE(dum2,10,3).c_str());
    }

    if (icall == 5) {
        // 105 format(' ','  could not find "to" cell center in "from" block',i3,
        //   /,'   will try another block')
        std::fprintf(f, "   could not find \"to\" cell center in \"from\" block%s\n"
                        "   will try another block\n", f_iw(idum1,3).c_str());
    }

    if (icall == 6) {
        // 106 format(' ',4x,'checking to see if this point really lies in block ',i3)
        std::fprintf(f, "     checking to see if this point really lies"
                        " in block %s\n", f_iw(idum1,3).c_str());
    }

    if (icall == 7) {
        // 107 format(' ','    point actually lies in block ',i3)
        std::fprintf(f, "     point actually lies in block %s\n",
                     f_iw(idum1,3).c_str());
        // 1071 format(' ','    future searches to be done in this block')
        std::fprintf(f, "     future searches to be done in this block\n");
    }

    if (icall == 8) {
        // 108 format(' ','    incompatable orientation of "from" cell and "to" cell')
        std::fprintf(f, "     incompatable orientation of \"from\" cell"
                        " and \"to\" cell\n");
        // 1081 format(' ','      ...will try in a nearby cell')
        std::fprintf(f, "       ...will try in a nearby cell\n");
    }

    if (icall == 9) {
        // 109 format(' ',6x,'using y-z to invert in "from" grid',i3,'  cell j,k=',i3,',',i3)
        std::fprintf(f, "       using y-z to invert in \"from\" grid%s"
                        "  cell j,k=%s,%s\n", f_iw(idum1,3).c_str(),
                     f_iw(idum2,3).c_str(), f_iw(idum3,3).c_str());
    }

    if (icall == 10) {
        // 1010 format(' ',6x,'using x-z to invert in "from" grid',i3,'  cell j,k=',i3,',',i3)
        std::fprintf(f, "       using x-z to invert in \"from\" grid%s"
                        "  cell j,k=%s,%s\n", f_iw(idum1,3).c_str(),
                     f_iw(idum2,3).c_str(), f_iw(idum3,3).c_str());
    }

    if (icall == 11) {
        // 1011 format(' ',6x,'using x-y to invert in "from" grid',i3,'  cell j,k=',i3,',',i3)
        std::fprintf(f, "       using x-y to invert in \"from\" grid%s"
                        "  cell j,k=%s,%s\n", f_iw(idum1,3).c_str(),
                     f_iw(idum2,3).c_str(), f_iw(idum3,3).c_str());
    }

    if (icall == 41) {
        // 1041 format(' ',4x,'search routine off track...using minimum distance search to get back on track')
        std::fprintf(f, "     search routine off track...using minimum"
                        " distance search to get back on track\n");
    }

    if (icall == 42) {
        // 1042 format(' ',4x,'frozen convergence...attempting to break cycle via minimum distance search')
        std::fprintf(f, "     frozen convergence...attempting to break"
                        " cycle via minimum distance search\n");
    }

    if (icall == 50) {
        // 1050 format(' ','marking point j,k=',i3,i3,' as an orphan ')
        std::fprintf(f, " marking point j,k=%s%s as an orphan \n",
                     f_iw(idum1,3).c_str(), f_iw(idum2,3).c_str());
    }

    if (icall == 99) {
        if (idum1 != idum3 || idum2 != idum4) {
            // 1099 format(' ','branch cut modification: before, jp,kp=',i3,i3,
            //   ' after, jpc,kpc=',i3,i3)
            std::fprintf(f, " branch cut modification: before, jp,kp=%s%s"
                            " after, jpc,kpc=%s%s\n",
                         f_iw(idum1,3).c_str(), f_iw(idum2,3).c_str(),
                         f_iw(idum3,3).c_str(), f_iw(idum4,3).c_str());
        }
    }
    return;
}
