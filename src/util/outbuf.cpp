// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "outbuf.h"
#include "my_flush.h"
#include "runtime/fortran_io.h"
#include <string>

namespace outbuf_ns {

void outbuf(FortranArray1DRef<char[1]> str, int& iunit)
{
    // Find the last non-blank character (loop from 120 down to 1)
    int i = 120;
    for (; i >= 1; --i) {
        if (str(i)[0] != ' ') {
            break;
        }
    }

    // Write characters 1..i to iunit using format (120a1)
    std::string output_str;
    for (int k = 1; k <= i; ++k) {
        output_str += str(k)[0];
    }
    fortran_write_unit(iunit, "%s", output_str.c_str());

    // Call my_flush(iunit)
    my_flush_ns::my_flush(iunit);
}

} // namespace outbuf_ns
