// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// outbuf.cpp — write one buffered character*120 line, trailing blanks trimmed.
#include "ron_common.h"

void outbuf(const FStr& str,int iunit){
    // find last non-blank char in the 120-wide field (Fortran loop i=120..1)
    std::string s=str.s; s.resize(120,' ');
    int i=120;
    for(; i>=1; --i){ if(s[i-1]!=' ') break; }
    // write chars 1..i via FORMAT(120a1) then flush. If all blank (i==0),
    // Fortran writes an empty record (just a newline).
    FILE* f=FUNIT(iunit);
    if(f){
        for(int k=1;k<=i;++k) std::fputc(s[k-1],f);
        std::fputc('\n',f);
    }
    my_flush(iunit);
}
