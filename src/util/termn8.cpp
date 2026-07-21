// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "termn8.h"
#include "outbuf.h"
#include "my_flush.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstdlib>

namespace termn8_ns {

void termn8(int& myid_stop, int& ierrflg, int& ibufdim, int& nbuf,
            FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou)
{
    // common /mydist2/ nnodes,myhost,myid,mycomm
    int& myid = cmn_mydist2.myid;

    // write(99,99) ierrflg
    // 99 format(' error code:',/,i4)
    fortran_write_unit(99, " error code:\n%4d\n", ierrflg);

    if (ierrflg == 0) {
        // write(99,1)
        // 1 format(/,' execution terminated normally')
        fortran_write_unit(99, "\n execution terminated normally\n");
    } else {
        if (ierrflg != -999) {
            if (ierrflg == -99) {
                if (nou(1) > 0) {
                    // write(99,22)
                    // 22 format(/,' abnormal termination during array',
                    //  .        ' sizing (precfl3d)',/,
                    //  .        ' (error message follows)',/)
                    fortran_write_unit(99,
                        "\n abnormal termination during array sizing (precfl3d)\n"
                        " (error message follows)\n\n");
                } else {
                    // write(99,33)
                    // 33 format(/,' abnormal termination during array',
                    //  .        ' sizing (precfl3d)',/,
                    //  .        ' (see precfl3d.out for error message)',/)
                    fortran_write_unit(99,
                        "\n abnormal termination during array sizing (precfl3d)\n"
                        " (see precfl3d.out for error message)\n\n");
                }
            } else {
                int iflag = 0;
                for (int nn = 1; nn <= nbuf; nn++) {
                    if (nou(nn) > 0) iflag = 1;
                }
                if (iflag > 0) {
                    // write(99,2)
                    // 2 format(/,' abnormal termination due to cfl3d',
                    //  .       ' error check',/,
                    //  .       ' (error message follows)',/)
                    fortran_write_unit(99,
                        "\n abnormal termination due to cfl3d error check\n"
                        " (error message follows)\n\n");
                } else {
                    // write(99,3)
                    // 3 format(/,' abnormal termination due to cfl3d',
                    //  .       ' error check',/,
                    //  .       ' (see main output file for error message)',/,
                    //  .       ' (also check dynamic patch output file if',
                    //  .       ' applicable)',/)
                    fortran_write_unit(99,
                        "\n abnormal termination due to cfl3d error check\n"
                        " (see main output file for error message)\n"
                        " (also check dynamic patch output file if applicable)\n\n");
                }
            }
        } else {
            // write(99,4)
            // 4 format(/,' abnormal termination due to receipt of',
            //  .       ' system signal',/,
            //  .       ' (kill, floating pt. exception,',
            //  .       ' segmentation fault, etc.)',/)
            fortran_write_unit(99,
                "\n abnormal termination due to receipt of system signal\n"
                " (kill, floating pt. exception, segmentation fault, etc.)\n\n");
        }

        // dump all internal buffers
        for (int nn = 1; nn <= nbuf; nn++) {
            if (nou(nn) > 0) {
                // write(99,*)
                fortran_write_unit(99, "\n");
                if (nn == 1) {
                    // write(99,*)'dump of unit 11 (main output) buffer:'
                    fortran_write_unit(99, " dump of unit 11 (main output) buffer:\n");
                } else if (nn == 2) {
                    // write(99,*)'dump of unit 09 (fort.9, dyn. patch) buffer:'
                    fortran_write_unit(99, " dump of unit 09 (fort.9, dyn. patch) buffer:\n");
                } else if (nn == 3) {
                    // write(99,*)'dump of unit 14 (baldwin-lomax) buffer:'
                    fortran_write_unit(99, " dump of unit 14 (baldwin-lomax) buffer:\n");
                } else if (nn == 4) {
                    // write(99,*)'dump of unit 25 (cfl3d.dynamic_patch) buffer:'
                    fortran_write_unit(99, " dump of unit 25 (cfl3d.dynamic_patch) buffer:\n");
                } else {
                    // write(99,*)'need a message for buffer',nn
                    fortran_write_unit(99, " need a message for buffer%d\n", nn);
                }
                // write(99,*)
                fortran_write_unit(99, "\n");
                for (int kou = 1; kou <= nou(nn); kou++) {
                    // call outbuf(bou(kou,nn),99)
                    // bou(kou,nn) is char[120]; outbuf expects FortranArray1DRef<char[1]>
                    FortranArray1DRef<char[1]> bou_row(
                        reinterpret_cast<char(*)[1]>(&bou(kou, nn)),
                        120);
                    int unit99 = 99;
                    outbuf_ns::outbuf(bou_row, unit99);
                }
            }
        }
    }

    // call my_flush for each unit
    { int u =  3; my_flush_ns::my_flush(u); }
    { int u =  4; my_flush_ns::my_flush(u); }
    { int u = 11; my_flush_ns::my_flush(u); }
    { int u = 12; my_flush_ns::my_flush(u); }
    { int u = 13; my_flush_ns::my_flush(u); }
    { int u = 14; my_flush_ns::my_flush(u); }
    { int u = 15; my_flush_ns::my_flush(u); }
    { int u = 17; my_flush_ns::my_flush(u); }
    { int u = 20; my_flush_ns::my_flush(u); }
    { int u = 23; my_flush_ns::my_flush(u); }
    { int u = 24; my_flush_ns::my_flush(u); }
    { int u = 25; my_flush_ns::my_flush(u); }
    { int u = 66; my_flush_ns::my_flush(u); }
    { int u = 96; my_flush_ns::my_flush(u); }
    { int u = 97; my_flush_ns::my_flush(u); }
    { int u = 98; my_flush_ns::my_flush(u); }
    { int u = 99; my_flush_ns::my_flush(u); }

    if (myid == myid_stop) {
        fortran_close_unit(1);
        fortran_close_unit(3);
        fortran_close_unit(4);
        fortran_close_unit(11);
        fortran_close_unit(12);
        fortran_close_unit(13);
        fortran_close_unit(15);
        fortran_close_unit(17);
        fortran_close_unit(20);
        fortran_close_unit(23);
        fortran_close_unit(24);
        fortran_close_unit(25);
        fortran_close_unit(66);
        fortran_close_unit(96);
        fortran_close_unit(97);
        fortran_close_unit(98);
        fortran_close_unit(99);
    }

    if (ierrflg == 0) {
        // stop
        std::exit(0);
    } else if (ierrflg == -99) {
        // write (6,*) 'program termination due to a precfl3d error check'
        // write (6,*) 'see file precfl3d.error'
        // stop
        fortran_write_unit(6, " program termination due to a precfl3d error check\n");
        fortran_write_unit(6, " see file precfl3d.error\n");
        std::exit(1);
    } else if (ierrflg > -99) {
        // write (6,*) 'program termination due to a cfl3d error check'
        // write (6,*) 'see file cfl3d.error'
        // stop
        fortran_write_unit(6, " program termination due to a cfl3d error check\n");
        fortran_write_unit(6, " see file cfl3d.error\n");
        std::exit(1);
    } else if (ierrflg == -999) {
        // write (6,*) 'program termination due to receipt of a system signal'
        // write (6,*) 'see files precfl3d.error and cfl3d.error'
        // stop
        fortran_write_unit(6, " program termination due to receipt of a system signal\n");
        fortran_write_unit(6, " see files precfl3d.error and cfl3d.error\n");
        std::exit(1);
    }

    // stop (final, catches any remaining case)
    std::exit(0);
}

} // namespace termn8_ns
