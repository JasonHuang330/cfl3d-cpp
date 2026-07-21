// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "umalloc.h"
#include "my_flush.h"
#include <cstdio>
#include <cstdlib>

namespace umalloc_ns {

void umalloc(int& n_words, int& intflag, char (&text)[1], int& memuse, int& status)
{
    // PARAMETER constants
    const int len_i  = 4;
    const int len_r4 = 4;

    // DBLE_PRECSN conditional
#ifdef DBLE_PRECSN
    int len_r = 8;
#else
    int len_r = 4;
#endif

    // Set len based on intflag
    int len = 0;
    if (intflag == 1) {
        len = len_i;
    } else if (intflag == 0) {
        len = len_r;
    } else if (intflag == -1) {
        len = len_r4;
    }

    int n_bytes = n_words * len;

    // Check status of allocate call
    if (status != 0) {
        fprintf(stdout, "stopping...failed trying to allocate memory for array %-20s\n", text);
        if (n_bytes < 1e10) {
            // Format 12: 'requested ', i10, ' words, (',i11, ' bytes)'
            fprintf(stdout, "requested %10d words, (%11d bytes)\n", n_words, n_bytes);
        } else if (n_bytes < 1e15) {
            // Format 13: 'requested ', i15, ' words, (',i16, ' bytes)'
            fprintf(stdout, "requested %15d words, (%16d bytes)\n", n_words, n_bytes);
        }
        // Format 14: 'prior to this call, total allocation was', i16, ' bytes'
        fprintf(stdout, "prior to this call, total allocation was%16d bytes\n", memuse);
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
