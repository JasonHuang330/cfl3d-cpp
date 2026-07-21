// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// my_flush.cpp — flush an output buffer (sequential build: GENERIC -> flush).
#include "ron_common.h"

void my_flush(int iunit){
    FILE* f=FUNIT(iunit);
    if (f) std::fflush(f);
}
