// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// usrint.cpp — signal handler: abnormal termination on system signal.
#include "ron_common.h"

void usrint(int /*sig*/){
    // errfile = 'ronnie.error'; open unit 99
    FILE* f99=std::fopen("ronnie.error","w");
    if(f99){
        int ierrflg=-999;
        // format 99: ' error code:',/,i4
        std::fprintf(f99," error code:\n%s\n",f_iw(ierrflg,4).c_str());
        // format 1
        std::fprintf(f99,
          "\n abnormal termination due to receipt of system signal\n"
          " (kill, floating pt. exception, segmentation fault, etc.)\n\n");
        std::fflush(f99);
    }
    my_flush(66);
    my_flush(25);
    my_flush(99);
    if(FUNIT(66)){ std::fclose(FUNIT(66)); FUNIT(66)=nullptr; }
    if(FUNIT(25)){ std::fclose(FUNIT(25)); FUNIT(25)=nullptr; }
    if(f99) std::fclose(f99);
    if(FUNIT(99)){ FUNIT(99)=nullptr; }
    std::exit(1);
}
