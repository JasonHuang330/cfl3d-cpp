// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "usrint.h"
#include "my_flush.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstdlib>

namespace usrint_ns {

void usrint() {
    // character*80 errfile
    char errfile[81] = "cfl3d.error";

    // open(unit=99,file=errfile,form='formatted',status='unknown')
    fortran_open_unit(99, errfile, "w");

    // ierrflg = -999
    int ierrflg = -999;

    // write(99,99) ierrflg
    // 99 format(' error code:',/,i4)
    {
        FILE* f99 = fortran_get_unit(99);
        fprintf(f99, " error code:\n%4d\n", ierrflg);
    }

    // write(99,1)
    // 1 format(/,' abnormal termination due to receipt of',
    //  .         ' system signal',/,
    //  .         ' (kill, floating pt. exception,',
    //  .         ' segmentation fault, etc.)',/)
    {
        FILE* f99 = fortran_get_unit(99);
        fprintf(f99, "\n abnormal termination due to receipt of system signal\n"
                     " (kill, floating pt. exception, segmentation fault, etc.)\n\n");
    }

    // call my_flush(3), my_flush(4), ..., my_flush(99)
    {
        int u;
        u = 3;  my_flush_ns::my_flush(u);
        u = 4;  my_flush_ns::my_flush(u);
        u = 11; my_flush_ns::my_flush(u);
        u = 12; my_flush_ns::my_flush(u);
        u = 13; my_flush_ns::my_flush(u);
        u = 14; my_flush_ns::my_flush(u);
        u = 15; my_flush_ns::my_flush(u);
        u = 17; my_flush_ns::my_flush(u);
        u = 20; my_flush_ns::my_flush(u);
        u = 23; my_flush_ns::my_flush(u);
        u = 24; my_flush_ns::my_flush(u);
        u = 25; my_flush_ns::my_flush(u);
        u = 99; my_flush_ns::my_flush(u);
    }

    // stop
    std::exit(0);
}

} // namespace usrint_ns
