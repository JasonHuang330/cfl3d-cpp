// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "umalloc.h"
#include "my_flush.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstdlib>

namespace umalloc_ns {

void umalloc(int& n_words, int& intflag, char (&text)[1], int& memuse, int& status)
{
    // parameter(len_i=4,len_r4=4)
    const int len_i  = 4;
    const int len_r4 = 4;

    int len_r = 4;

    int len = 0;
    if (intflag == 1) {
        len = len_i;
    } else if (intflag == 0) {
        len = len_r;
    } else if (intflag == -1) {
        len = len_r4;
    }

    int n_bytes = n_words * len;

    // check status of allocate call
    if (status != 0) {
        fprintf(fortran_get_unit(6),
                "stopping...failed trying to allocate memory for array %-20.20s\n",
                text);

        if (n_bytes < 1e10) {
            // write(6,12) n_words, n_bytes
            // 12 format('requested ', i10, ' words, (',i11, ' bytes)')
            fprintf(fortran_get_unit(6),
                    "requested %10d words, (%11d bytes)\n",
                    n_words, n_bytes);
        } else if (n_bytes < 1e15) {
            // write(6,13) n_words, n_bytes
            // 13 format('requested ', i15, ' words, (',i16, ' bytes)')
            fprintf(fortran_get_unit(6),
                    "requested %15d words, (%16d bytes)\n",
                    n_words, n_bytes);
        }

        // write(6,14) memuse
        // 14 format('prior to this call, total allocation was', i16, ' bytes')
        fprintf(fortran_get_unit(6),
                "prior to this call, total allocation was%16d bytes\n",
                memuse);

        int iunit6 = 6;
        my_flush_ns::my_flush(iunit6);

        std::exit(1);

    } else {
        // write(6,'(''allocating:'',a20,i15,'' bytes'')') text,n_bytes
        // call my_flush(6)
        // (commented out in original)
    }

    memuse = memuse + n_bytes;
}

} // namespace umalloc_ns
