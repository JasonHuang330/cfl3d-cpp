// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "umalloc_r.h"
#include "umalloc.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include <cstdio>
#include <cstdlib>

namespace umalloc_r_ns {

// umalloc: forward to umalloc_ns::umalloc
void umalloc(int& n_words, int& intflag, char (&text)[1], int& memuse, int& status)
{
    umalloc_ns::umalloc(n_words, intflag, text, memuse, status);
}

void umalloc_r(int& n_words, int& intflag, char (&text)[1], int& memuse, int& status)
{
    // parameter(len_i=4,len_r4=4)
    const int len_i  = 4;
    const int len_r4 = 4;

    // common /mydist2/ nnodes,myhost,myid,mycomm
    int& nnodes = cmn_mydist2.nnodes;
    int& myhost = cmn_mydist2.myhost;
    int& myid   = cmn_mydist2.myid;
    int& mycomm = cmn_mydist2.mycomm;

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
        FILE* f6 = fortran_get_unit(6);
        fprintf(f6, "stopping...failed trying to allocate memory for array %-20.20s\n", text);
        if (n_bytes < (int)1e10) {
            // format 12: 'requested ', i10, ' words, (',i11, ' bytes)'
            fprintf(f6, "requested %10d words, (%11d bytes)\n", n_words, n_bytes);
        } else if (n_bytes < (int)1e15) {
            // format 13: 'requested ', i15, ' words, (',i16, ' bytes)'
            fprintf(f6, "requested %15d words, (%16d bytes)\n", n_words, n_bytes);
        }
        // format 14: 'prior to this call, total allocation was', i16, ' bytes'
        fprintf(f6, "prior to this call, total allocation was%16d bytes\n", memuse);
        // call my_flush(6)  -- commented out in Fortran
        std::exit(1);
    } else {
        // write(6,...) commented out in Fortran
        // call my_flush(6)  -- commented out in Fortran
    }

    memuse = memuse + n_bytes;
}

} // namespace umalloc_r_ns
