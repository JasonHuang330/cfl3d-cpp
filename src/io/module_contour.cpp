// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "module_contour.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

namespace module_contour_ns {

void cont_plot(int& nframe, int& jdim, int& kdim, int& idim, int& nummem, int& nblk,
               FortranArray4DRef<double> q, FortranArray3DRef<double> x,
               FortranArray3DRef<double> y, FortranArray3DRef<double> z,
               FortranArray3DRef<double> vist3d, FortranArray4DRef<double> ux,
               FortranArray4DRef<double> turre, FortranArray3DRef<double> smin)
{
    // COMMON block references
    int32_t& i2d    = cmn_twod.i2d;
    float&   pr     = cmn_fluid2.pr;
    float&   prt    = cmn_fluid2.prt;
    float&   cbar   = cmn_fluid2.cbar;
    float&   reue   = cmn_reyue.reue;
    float&   tinf   = cmn_reyue.tinf;
    // ivisc(3) → cmn_reyue.ivisc[0..2] (not used in this routine)
    float&   gamma  = cmn_fluid.gamma;
    float&   gm1    = cmn_fluid.gm1;
    float&   gp1    = cmn_fluid.gp1;
    float&   gm1g   = cmn_fluid.gm1g;
    float&   gp1g   = cmn_fluid.gp1g;
    float&   ggm1   = cmn_fluid.ggm1;
    // /info/ title(20),rkap(3),xmach
    float&   xmach  = cmn_info.xmach;
    // /ivals/ — not directly used in computation below (vel1 declared but unused)
    float&   p0     = cmn_ivals.p0;
    float&   rho0   = cmn_ivals.rho0;
    float&   c0     = cmn_ivals.c0;
    float&   u0     = cmn_ivals.u0;
    float&   v0     = cmn_ivals.v0;
    float&   w0     = cmn_ivals.w0;
    float&   et0    = cmn_ivals.et0;
    float&   h0     = cmn_ivals.h0;
    float&   pt0    = cmn_ivals.pt0;
    float&   rhot0  = cmn_ivals.rhot0;
    // qiv(5) → cmn_ivals.qiv[0..4], tur10(7) → cmn_ivals.tur10[0..6]

    // Local variables
    char str[21];
    char strframe[21];
    char plt_title[121];
    int i, j, k, m;
    double xc, yc, zc;
    double c2b, c2bp, tt, re, vel, tau_w, fnuloc, deltabl, coef;
    double vel1;
    // fnu(jdim), utau(jdim), Rel(jdim), cf(jdim) — declared but not used in body
    // (they appear in the declaration but no statements reference them)

    // IF(nummem/=2.and.nummem/=7) RETURN
    if (nummem != 2 && nummem != 7) return;

    // WRITE(str,'(I5.5)') nblk
    std::snprintf(str, sizeof(str), "%05d", nblk);

    // WRITE(strframe,'(I4.4)') nframe
    std::snprintf(strframe, sizeof(strframe), "%04d", nframe);

    // OPEN(1234,file="turre-"//TRIM(ADJUSTL(strframe))//"-blk"//TRIM(ADJUSTL(str))//".plt")
    {
        char filename[256];
        std::snprintf(filename, sizeof(filename), "turre-%s-blk%s.plt", strframe, str);
        fortran_open_unit(1234, filename, "w");
    }

    FILE* funit = fortran_get_unit(1234);

    // WRITE(plt_title,'(2(1x,A,es12.5))') 'Re=',reue,'Mach=',xmach
    std::snprintf(plt_title, sizeof(plt_title), " Re= %12.5e Mach= %12.5e",
                  (double)reue, (double)xmach);

    // WRITE(1234,'(A)') 'title="'//TRIM(ADJUSTL(plt_title))//'"'
    // TRIM(ADJUSTL(plt_title)): strip leading and trailing spaces
    {
        // ADJUSTL: remove leading spaces
        const char* p = plt_title;
        while (*p == ' ') p++;
        // TRIM: remove trailing spaces
        int len = (int)std::strlen(p);
        while (len > 0 && p[len-1] == ' ') len--;
        std::fprintf(funit, "title=\"%.*s\"\n", len, p);
    }

    // WRITE variables line
    if (nummem == 7) {
        std::fprintf(funit, "%s\n",
            "variables=\"x\",\"y\",\"z\",\"dist\",\"q1\", \"q2\", \"q3\", \"q4\", \"q5\", "
            "\"dudx\",\"dudy\",\"dudz\",\"dvdx\",\"dvdy\",\"dvdz\",\"dwdx\",\"dwdy\",\"dwdz\", "
            "\"mu\", \"mut\", \"t11\", \"t22\", \"t33\", \"t12\", \"t23\", \"t13\", \"omega/eps/zeta\"");
    } else if (nummem == 2) {
        std::fprintf(funit, "%s\n",
            "variables=\"x\",\"y\",\"z\",\"dist\",\"q1\", \"q2\", \"q3\", \"q4\", \"q5\", "
            "\"dudx\",\"dudy\",\"dudz\",\"dvdx\",\"dvdy\",\"dvdz\",\"dwdx\",\"dwdy\",\"dwdz\", "
            ",\"mu\",\"mut\", \"omg/eps/zeta\", \"tke\"");
    }

    // ZONE line
    if (i2d == 1) {
        // WRITE(1234,'(A,i5,A,i5)') "ZONE  I=", jdim-1," J=", kdim-1
        std::fprintf(funit, "ZONE  I=%5d J=%5d\n", jdim-1, kdim-1);
    } else {
        // WRITE(1234,'(A,i5,A,i5,A,i5)') "ZONE  I=", idim-1, " J=",jdim-1," K=", kdim-1
        std::fprintf(funit, "ZONE  I=%5d J=%5d K=%5d\n", idim-1, jdim-1, kdim-1);
    }

    // Triple loop
    for (k = 1; k <= kdim-1; k++) {
        for (j = 1; j <= jdim-1; j++) {
            for (i = 1; i <= idim-1; i++) {
                // Cell-center coordinates (average of 8 corners)
                xc = 0.125*(x(j,k,i)+x(j+1,k,i)+x(j,k+1,i)+x(j,k,i+1) +
                            x(j+1,k+1,i) + x(j+1,k,i+1) + x(j,k+1,i+1)+x(j+1,k+1,i+1));
                yc = 0.125*(y(j,k,i)+y(j+1,k,i)+y(j,k+1,i)+y(j,k,i+1) +
                            y(j+1,k+1,i) + y(j+1,k,i+1) + y(j,k+1,i+1)+y(j+1,k+1,i+1));
                zc = 0.125*(z(j,k,i)+z(j+1,k,i)+z(j,k+1,i)+z(j,k,i+1) +
                            z(j+1,k+1,i) + z(j+1,k,i+1) + z(j,k+1,i+1)+z(j+1,k+1,i+1));

                c2b  = (double)cbar / (double)tinf;
                c2bp = c2b + 1.0;
                tt   = (double)gamma * q(j,k,i,5) / q(j,k,i,1);
                fnuloc = c2bp * tt * std::sqrt(tt) / (c2b + tt);

                // WRITE(1234,'(25(1x,es12.5))') xc,yc,zc,smin(j,k,i),q(j,k,i,1:5),
                //   ux(j,k,i,1:9), fnuloc,vist3d(j,k,i), turre(j,k,i,1:nummem)
                std::fprintf(funit,
                    " %12.5e %12.5e %12.5e %12.5e"
                    " %12.5e %12.5e %12.5e %12.5e %12.5e"
                    " %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e"
                    " %12.5e %12.5e",
                    xc, yc, zc, smin(j,k,i),
                    q(j,k,i,1), q(j,k,i,2), q(j,k,i,3), q(j,k,i,4), q(j,k,i,5),
                    ux(j,k,i,1), ux(j,k,i,2), ux(j,k,i,3),
                    ux(j,k,i,4), ux(j,k,i,5), ux(j,k,i,6),
                    ux(j,k,i,7), ux(j,k,i,8), ux(j,k,i,9),
                    fnuloc, vist3d(j,k,i));
                for (m = 1; m <= nummem; m++) {
                    std::fprintf(funit, " %12.5e", turre(j,k,i,m));
                }
                std::fprintf(funit, "\n");
            }
        }
    }

    // CLOSE(1234)
    fortran_close_unit(1234);
}

} // namespace module_contour_ns
