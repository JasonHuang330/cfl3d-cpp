// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "outbuf.h"
#include "my_flush.h"
#include "runtime/fortran_io.h"
#include <cstdio>

namespace outbuf_ns {

void outbuf(FortranArray1DRef<char[1]> str, int& iunit) {
    int i = 0;
    for (int k = 120; k >= 1; k--) {
        if (str(k)[0] != ' ') {
            i = k;
            break;
        }
    }
    // write(iunit,1) (str(k),k=1,i)
    // format(120a1) — write i characters to iunit
    FILE* fp = fortran_get_unit(iunit);
    for (int k = 1; k <= i; k++) {
        fputc(str(k)[0], fp);
    }
    fputc('\n', fp);
    my_flush_ns::my_flush(iunit);
}

} // namespace outbuf_ns
