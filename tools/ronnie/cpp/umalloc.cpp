// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// umalloc.cpp — track memory allocated (sequential, single precision build).
// intflag: 1 -> integer (4 bytes), 0 -> real (single-precision build: 4 bytes),
//         -1 -> real*4 (4 bytes).  status is always 0 here (allocate succeeds).
#include "ron_common.h"

void umalloc(long n_words,int intflag,const char* text,long& memuse,int status){
    const int len_i=4, len_r4=4;
    const int len_r=4;               // not DBLE_PRECSN, not CRAY, not CMPLX
    int len=0;
    if      (intflag== 1) len=len_i;
    else if (intflag== 0) len=len_r;
    else if (intflag==-1) len=len_r4;
    long n_bytes = n_words * (long)len;

    if (status!=0){
        std::fprintf(stdout,
          "stopping...failed trying to allocate memory for array %-20s\n",text);
        if (n_bytes < (long)1e10)
            std::fprintf(stdout,"requested %10ld words, (%11ld bytes)\n",
                         n_words,n_bytes);
        else if (n_bytes < (long)1e15)
            std::fprintf(stdout,"requested %15ld words, (%16ld bytes)\n",
                         n_words,n_bytes);
        std::fprintf(stdout,
          "prior to this call, total allocation was%16ld bytes\n",memuse);
        std::fflush(stdout);
        std::exit(1);
    }
    memuse += n_bytes;
}
