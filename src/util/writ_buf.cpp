// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "writ_buf.h"
#include "outbuf.h"
#include "runtime/fortran_io.h"
#include <cstring>

namespace writ_buf_ns {

void writ_buf(int& nbl, int& iunit, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myhost, int& myid, int& mycomm, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    int iou = 1;
    if (iunit == 11) iou = 1;
    if (iunit ==  9) iou = 2;
    if (iunit == 14) iou = 3;
    if (iunit == 25) iou = 4;

    if (nou(iou) > 0) {
        for (int kou = 1; kou <= nou(iou); kou++) {
            // bou(kou,iou) is char[120]; outbuf expects FortranArray1DRef<char[1]>
            // Treat the 120-char element as a 1D array of 120 char[1] elements
            FortranArray1DRef<char[1]> str_ref(reinterpret_cast<char(*)[1]>(&bou(kou, iou)), 1, 120);
            outbuf_ns::outbuf(str_ref, iunit);
        }

        if (nou(iou) >= ibufdim) {
            fortran_write_unit(iunit, "WARNING: internal buffer length exceeded -  make parameter ibufdim > %6d\n", nou(iou));
            fortran_write_unit(iunit, "continuing, but you will be missing some output data \n");
        }

        for (int kou = 1; kou <= nou(iou); kou++) {
            std::memset(&bou(kou, iou), ' ', 120);
        }
        nou(iou) = 0;
    }
}

} // namespace writ_buf_ns
