// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "termn8.h"
#include "outbuf.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstdlib>

namespace termn8_ns {

void termn8(int& myid_stop, int& ierrflg, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou)
{
    // write(99,99) ierrflg
    // 99 format('error code:',/,i3)
    fortran_write_unit(99, "error code:\n%3d\n", ierrflg);

    // *******************************
    // normal termination section
    // *******************************
    if (ierrflg == 0) {
        // write(99,1)
        // 1 format('execution terminated normally')
        fortran_write_unit(99, "execution terminated normally\n");
        // close(8,status='delete')
        fortran_close_unit(8);
        // close(20,status='delete')
        fortran_close_unit(20);
        // close(10)
        fortran_close_unit(10);
        // close(7)
        fortran_close_unit(7);
        // close(15)
        fortran_close_unit(15);
        // close(9)
        fortran_close_unit(9);
        // close(30)
        fortran_close_unit(30);
        // close(40)
        fortran_close_unit(40);
        // close(50)
        fortran_close_unit(50);
        // close(60)
        fortran_close_unit(60);
        // close(99)
        fortran_close_unit(99);
        // stop
        std::exit(0);
    }

    // *********************************
    // abnormal termination section
    // *********************************
    if (ierrflg < 0) {
        // write(99,100)
        // 100 format('abnormal termination:')
        fortran_write_unit(99, "abnormal termination:\n");
    }
    if (ierrflg == -1) {
        // write(99,101)
        // 101 format('input file error')
        fortran_write_unit(99, "input file error\n");
    } else if (ierrflg == -2) {
        // write(99,102)
        // 102 format('no cfl3d input file specified')
        fortran_write_unit(99, "no cfl3d input file specified\n");
    } else if (ierrflg == -3) {
        // write(99,103)
        // 103 format('parameter mbloc too small')
        fortran_write_unit(99, "parameter mbloc too small\n");
    } else if (ierrflg == -4) {
        // write(99,104)
        // 104 format('parameter npmax too small')
        fortran_write_unit(99, "parameter npmax too small\n");
    } else if (ierrflg == -5) {
        // write(99,105)
        // 105 format('inconsistant grid and grid sensitivity files')
        fortran_write_unit(99, "inconsistant grid and grid sensitivity files\n");
    } else if (ierrflg == -6) {
        // write(99,106)
        // 106 format('internal file conversion error')
        fortran_write_unit(99, "internal file conversion error\n");
    } else if (ierrflg == -7) {
        // write(99,107)
        // 107 format('parameter mseqn too small')
        fortran_write_unit(99, "parameter mseqn too small\n");
    } else {
        // write(99,999)
        // 999 format('miscellaneous error')
        fortran_write_unit(99, "miscellaneous error\n");
    }

    // do nn=1,nbuf
    //    if (nou(nn).gt.0) then
    //       write(99,*)
    //       do kou = 1,nou(nn)
    //          call outbuf(bou(kou,nn),99)
    //       end do
    //    end if
    // end do
    for (int nn = 1; nn <= nbuf; nn++) {
        if (nou(nn) > 0) {
            // write(99,*) — blank line
            fortran_write_unit(99, "\n");
            for (int kou = 1; kou <= nou(nn); kou++) {
                // call outbuf(bou(kou,nn),99)
                // bou(kou,nn) is char[120]; outbuf expects FortranArray1DRef<char[1]>
                // Reinterpret the 120-byte buffer as an array of 120 char[1] elements
                FortranArray1DRef<char[1]> str(
                    reinterpret_cast<char(*)[1]>(&bou(kou, nn)[0]),
                    120, 1);
                int iunit99 = 99;
                outbuf_ns::outbuf(str, iunit99);
            }
        }
    }

    // write (6, *)  'terminating the program'
    fortran_write_unit(6, " terminating the program\n");

    // if (ierrflg.eq.-7) then
    //    close(8)
    //    close(20)
    // else
    //    close(8,status='delete')
    //    close(20,status='delete')
    // end if
    if (ierrflg == -7) {
        // conversion error: close, but don't remove temporary files
        fortran_close_unit(8);
        fortran_close_unit(20);
    } else {
        // close and delete temporary files
        fortran_close_unit(8);
        fortran_close_unit(20);
    }
    // close(10)
    fortran_close_unit(10);
    // close(7)
    fortran_close_unit(7);
    // close(15)
    fortran_close_unit(15);
    // close(9)
    fortran_close_unit(9);
    // close(30)
    fortran_close_unit(30);
    // close(40)
    fortran_close_unit(40);
    // close(50)
    fortran_close_unit(50);
    // close(60)
    fortran_close_unit(60);
    // close(99)
    fortran_close_unit(99);

    // stop
    std::exit(0);
}

} // namespace termn8_ns
