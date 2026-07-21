// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "cputim.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include <cstdio>
#include <cmath>
#include <ctime>
#include <sys/times.h>
#include <unistd.h>
#include <algorithm>

namespace cputim_ns {

// Fortran COMMON /time1/ tim(3,3), tm(3) — accessed via cmn_time1
// tim(i,j) Fortran 1-based column-major → cmn_time1.tim[(j-1)*3 + (i-1)]
// tm(n)    Fortran 1-based              → cmn_time1.tm[n-1]
#define TIM(i,j) (cmn_time1.tim[((j)-1)*3 + ((i)-1)])
#define TM(n)    (cmn_time1.tm[(n)-1])

//
// cputim(icall, nodes, string, myhost, myid, mycomm, iunit)
//
// Purpose: Determine timings for cpu usage.
//   icall = 0  for initial call
//         < 0  for final call; total run time is output
//         + 1  for a call preceding an individual subroutine (nothing output)
//         + 2  for a call following the individual subroutine (time since last call output)
//   nodes  = number of nodes, if distributed
//   string = character string up to 50 characters identifying item being timed
//
// The authoritative C++ header declares this as void cputim() with no parameters.
// The full parameterized implementation is provided below as an internal helper.
//
static void cputim_impl(int icall, int /*nodes*/, const char* string,
                        int /*myhost*/, int /*myid*/, int /*mycomm*/, int iunit)
{
    float walltime, totaltime;
    int ia[3];  // ia(1)=hour, ia(2)=minute, ia(3)=second

    // initialize
    if (icall == 0) {
        for (int j = 1; j <= 3; j++) {
            for (int i = 1; i <= 3; i++) {
                TIM(i,j) = 0.0f;
            }
        }
    }

    // timing for workstations
    // etime(tm): returns elapsed time, fills tm(1)=user time, tm(2)=system time
    {
        struct tms tms_buf;
        times(&tms_buf);
        long clk_tck = sysconf(_SC_CLK_TCK);
        TM(1) = (float)((double)tms_buf.tms_utime / (double)clk_tck);
        TM(2) = (float)((double)tms_buf.tms_stime / (double)clk_tck);
        totaltime = TM(1) + TM(2);
    }

    // itime(ia): fills ia(1)=hour (0-23), ia(2)=minute (0-59), ia(3)=second (0-59)
    {
        time_t now = time(nullptr);
        struct tm* lt = localtime(&now);
        ia[0] = lt->tm_hour;   // ia(1) = hour
        ia[1] = lt->tm_min;    // ia(2) = minute
        ia[2] = lt->tm_sec;    // ia(3) = second
    }

    // tm(3) = ia(3) + 60*(ia(2) + 60*ia(1))
    TM(3) = (float)(ia[2] + 60 * (ia[1] + 60 * ia[0]));

    // first time here
    if ((float)TIM(1,3) == 0.0f) {
        for (int n = 1; n <= 3; n++) {
            TIM(n,3) = TM(n);
        }
    }

    // check if wall clock passed 24:00 (midnight) since the last call to this routine
    if ((float)TM(3) < (float)TIM(3,3)) {
        TM(3) = TM(3) + (float)(24 * 3600);
    }

    // get totals
    for (int n = 1; n <= 3; n++) {
        TIM(n,2) = TM(n) - TIM(n,3);
        TIM(n,1) = TIM(n,1) + TIM(n,2);
        TIM(n,3) = TM(n);
    }

    // output the collective timings for this run
    if (icall < 0) {
        //
        // best attempt at correct wall time for > 24 hour runs:
        // the problem seems to be that tim(3,2), supposedly the wall
        // time, gets reset after 24 hours. have to be creative in the case
        // the total time - system+user - is not *quite* 24 hours and
        // the true wall time is just over 24 hours, so that ndays=0
        // but tim(3,2) is small. in that case use the total time as
        // the wall time.
        //
        walltime  = TIM(3,1);
        totaltime = TIM(1,1) + TIM(2,1);
        int ndays = (int)((float)totaltime) / 86400;
        if (ndays > 0) {
            walltime = walltime + (float)(ndays * 86400);
        } else {
            walltime = std::max(walltime, totaltime);
        }
        FILE* fp = fortran_get_unit(iunit);
        fprintf(fp, "\n%-50.50s\n\n   node      user    system     total    wall clock\n",
                string);
        fprintf(fp, " %6d%10.2f%10.2f%10.2f%12.2f\n",
                0,
                (float)TIM(1,1),
                (float)TIM(2,1),
                (float)TIM(1,1) + (float)TIM(2,1),
                (float)walltime);
        float totu = TIM(1,1);
        float tots = TIM(2,1);
        float tott = TIM(1,1) + TIM(2,1);
        fprintf(fp, " ------------------------------------\n");
        fprintf(fp, " total:%10.2f%10.2f%10.2f\n",
                (float)totu, (float)tots, (float)tott);
        int nhrs      = (int)((float)walltime / 3600.0f);
        float nsecs_f = walltime - (float)(nhrs * 3600);
        int nmins     = (int)(nsecs_f / 60.0f);
        int nsecs     = (int)(nsecs_f - (float)(nmins * 60));
        fprintf(fp, "\n total run (wall) time = %4d hours %4d minutes %4d seconds\n",
                nhrs, nmins, nsecs);
    }

    // print time since last call to this routine
    if (icall == 2) {
        //
        // best attempt at correct wall time for > 24 hour runs
        //
        walltime  = TIM(3,2);
        totaltime = TIM(1,1) + TIM(2,1);
        int ndays = (int)((float)totaltime) / 86400;
        if (ndays > 0) {
            walltime = walltime + (float)(ndays * 86400);
        } else {
            walltime = std::max(walltime, totaltime);
        }
        FILE* fp = fortran_get_unit(iunit);
        fprintf(fp, "\n%-50.50s\n\n   node      user    system     total    wall clock\n",
                string);
        fprintf(fp, " %6d%10.2f%10.2f%10.2f%12.2f\n",
                0,
                (float)TIM(1,2),
                (float)TIM(2,2),
                (float)TIM(1,2) + (float)TIM(2,2),
                (float)walltime);
        float totu = TIM(1,2);
        float tots = TIM(2,2);
        float tott = TIM(1,2) + TIM(2,2);
        fprintf(fp, " ------------------------------------\n");
        fprintf(fp, " total:%10.2f%10.2f%10.2f\n",
                (float)totu, (float)tots, (float)tott);
    }

    return;
}

// Declared entry point matching the authoritative header: void cputim()
// Parameters (icall, nodes, string, myhost, myid, mycomm, iunit) are not
// available in the no-arg form. The full implementation is in cputim_impl().
void cputim()
{
    // No-argument entry point as declared in the authoritative header.
    // The Fortran subroutine cputim(icall,nodes,string,myhost,myid,mycomm,iunit)
    // requires parameters that are not available here. The full parameterized
    // implementation is in cputim_impl() above for use within this translation unit.
}

// Parameterized entry point for callers that pass arguments (e.g. mgbl.cpp)
void cputim(int& icall, int& nodes, char* string, int& myhost, int& myid, int& mycomm, int& iunit)
{
    cputim_impl(icall, nodes, string, myhost, myid, mycomm, iunit);
}

#undef TIM
#undef TM

} // namespace cputim_ns
