// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "csout.h"
#include "runtime/fortran_io.h"

namespace csout_ns {

void csout(int& iseq, int& maxbl, int& maxcs, FortranArray1DRef<int> igridg, FortranArray1DRef<int> levelg, int& ncs, FortranArray1DRef<double> sx, FortranArray1DRef<double> sy, FortranArray1DRef<double> sz, FortranArray1DRef<double> stot, FortranArray1DRef<double> pav, FortranArray1DRef<double> ptav, FortranArray1DRef<double> tav, FortranArray1DRef<double> ttav, FortranArray1DRef<double> xmav, FortranArray1DRef<double> fmdot, FortranArray1DRef<double> cfxp, FortranArray1DRef<double> cfyp, FortranArray1DRef<double> cfzp, FortranArray1DRef<double> cfdp, FortranArray1DRef<double> cflp, FortranArray1DRef<double> cftp, FortranArray1DRef<double> cfxv, FortranArray1DRef<double> cfyv, FortranArray1DRef<double> cfzv, FortranArray1DRef<double> cfdv, FortranArray1DRef<double> cflv, FortranArray1DRef<double> cftv, FortranArray1DRef<double> cfxmom, FortranArray1DRef<double> cfymom, FortranArray1DRef<double> cfzmom, FortranArray1DRef<double> cfdmom, FortranArray1DRef<double> cflmom, FortranArray1DRef<double> cftmom, FortranArray1DRef<double> cfxtot, FortranArray1DRef<double> cfytot, FortranArray1DRef<double> cfztot, FortranArray1DRef<double> cfdtot, FortranArray1DRef<double> cfltot, FortranArray1DRef<double> cfttot, FortranArray2DRef<int> icsinfo)
{
    // common /mgrd/
    int& levt    = cmn_mgrd.levt;
    int& lglobal = cmn_mgrd.lglobal;

    // print results for each surface on the current global/embedded level
    int ics1    = 0;
    int iprttot = 0;

    for (int ics = 1; ics <= ncs; ics++) {
        int nbl1 = icsinfo(ics, 1);
        if (levelg(nbl1) >= lglobal && levelg(nbl1) <= levt) {
            ics1 = ics1 + 1;
            if (icsinfo(ics, 9) != 0) iprttot = 1;

            // format 1000: /,/,' Control Surface',i3,/,' ===================',/
            fortran_write_unit(17, "\n\n Control Surface%3d\n ===================\n\n", ics1);

            int igrid = igridg(icsinfo(ics, 1));
            // format 1100: ' Grid',i3,' (Block',i3,')  i =',i3,',',i3,'  j =',i3,',',i3,'  k =',i3,',',i3,'  iwall = ',i3,'  Normal = ',i3,/
            fortran_write_unit(17,
                " Grid%3d (Block%3d)  i =%3d,%3d  j =%3d,%3d  k =%3d,%3d  iwall = %3d  Normal = %3d\n\n",
                igrid,
                icsinfo(ics, 1), icsinfo(ics, 2), icsinfo(ics, 3),
                icsinfo(ics, 4), icsinfo(ics, 5), icsinfo(ics, 6),
                icsinfo(ics, 7), icsinfo(ics, 8), icsinfo(ics, 9));

            // format 1113: ' x-area = ',e12.4,' y-area = ',e12.4,' z-area = ',e12.4,' total-area = ',e12.4
            fortran_write_unit(17,
                " x-area = %12.4e y-area = %12.4e z-area = %12.4e total-area = %12.4e\n",
                (float)sx(ics), (float)sy(ics), (float)sz(ics), (float)stot(ics));

            if (icsinfo(ics, 8) == 0) {
                // format 1114: /,' Mass averaged properties'
                fortran_write_unit(17, "\n Mass averaged properties\n");

                // format 1115: ' P/Pinf      = ',e12.4,'         Pt/Pinf = ',e12.4,/
                //              ' T/Tinf      = ',e12.4,'         Tt/Tinf = ',e12.4,/
                //              ' Mach number = ',e12.4
                fortran_write_unit(17,
                    " P/Pinf      = %12.4e         Pt/Pinf = %12.4e\n"
                    " T/Tinf      = %12.4e         Tt/Tinf = %12.4e\n"
                    " Mach number = %12.4e\n",
                    (float)pav(ics), (float)ptav(ics),
                    (float)tav(ics), (float)ttav(ics),
                    (float)xmav(ics));

                // format 1116: /,' Mass flow / (rhoinf*vinf*(L_R)**2) = ',e14.5
                fortran_write_unit(17,
                    "\n Mass flow / (rhoinf*vinf*(L_R)**2) = %14.5e\n",
                    (float)fmdot(ics));
            }

            // format 1200: /,22x,'x-force',7x,'y-force',7x,'z-force',3x,'resultant-force',1x,'lift-force',4x,'drag-force'
            fortran_write_unit(17,
                "\n                      x-force       y-force       z-force   resultant-force lift-force    drag-force\n");

            // format 1300: ' Pressure force  ',6e14.5
            fortran_write_unit(17,
                " Pressure force  %14.5e%14.5e%14.5e%14.5e%14.5e%14.5e\n",
                (float)cfxp(ics), (float)cfyp(ics), (float)cfzp(ics),
                (float)cftp(ics), (float)cflp(ics), (float)cfdp(ics));

            if (icsinfo(ics, 8) == 1) {
                // format 1400: ' Viscous force   ',6e14.5
                fortran_write_unit(17,
                    " Viscous force   %14.5e%14.5e%14.5e%14.5e%14.5e%14.5e\n",
                    (float)cfxv(ics), (float)cfyv(ics), (float)cfzv(ics),
                    (float)cftv(ics), (float)cflv(ics), (float)cfdv(ics));
            }

            if (icsinfo(ics, 8) == 0) {
                // format 1500: ' Thrust force    ',6e14.5
                fortran_write_unit(17,
                    " Thrust force    %14.5e%14.5e%14.5e%14.5e%14.5e%14.5e\n",
                    (float)cfxmom(ics), (float)cfymom(ics), (float)cfzmom(ics),
                    (float)cftmom(ics), (float)cflmom(ics), (float)cfdmom(ics));
            }

            // format 1600: ' Total force     ',6e14.5,/
            fortran_write_unit(17,
                " Total force     %14.5e%14.5e%14.5e%14.5e%14.5e%14.5e\n\n",
                (float)cfxtot(ics), (float)cfytot(ics), (float)cfztot(ics),
                (float)cfttot(ics), (float)cfltot(ics), (float)cfdtot(ics));
        }
    }

    // print out totals
    if (iprttot > 0) {
        // format 1700: /,/,' Totals for all sufaces (global level) with normal .ne. 0',/,
        //              ' =========================================','==============='/
        fortran_write_unit(17,
            "\n\n Totals for all sufaces (global level) with normal .ne. 0\n"
            " ========================================================\n\n");

        // format 1113
        fortran_write_unit(17,
            " x-area = %12.4e y-area = %12.4e z-area = %12.4e total-area = %12.4e\n",
            (float)sx(ncs + 1), (float)sy(ncs + 1), (float)sz(ncs + 1), (float)stot(ncs + 1));

        // format 1116
        fortran_write_unit(17,
            "\n Mass flow / (rhoinf*vinf*(L_R)**2) = %14.5e\n",
            (float)fmdot(ncs + 1));

        // format 1200
        fortran_write_unit(17,
            "\n                      x-force       y-force       z-force   resultant-force lift-force    drag-force\n");

        // format 1300
        fortran_write_unit(17,
            " Pressure force  %14.5e%14.5e%14.5e%14.5e%14.5e%14.5e\n",
            (float)cfxp(ncs + 1), (float)cfyp(ncs + 1), (float)cfzp(ncs + 1),
            (float)cftp(ncs + 1), (float)cflp(ncs + 1), (float)cfdp(ncs + 1));

        // format 1400
        fortran_write_unit(17,
            " Viscous force   %14.5e%14.5e%14.5e%14.5e%14.5e%14.5e\n",
            (float)cfxv(ncs + 1), (float)cfyv(ncs + 1), (float)cfzv(ncs + 1),
            (float)cftv(ncs + 1), (float)cflv(ncs + 1), (float)cfdv(ncs + 1));

        // format 1500
        fortran_write_unit(17,
            " Thrust force    %14.5e%14.5e%14.5e%14.5e%14.5e%14.5e\n",
            (float)cfxmom(ncs + 1), (float)cfymom(ncs + 1), (float)cfzmom(ncs + 1),
            (float)cftmom(ncs + 1), (float)cflmom(ncs + 1), (float)cfdmom(ncs + 1));

        // format 1600
        fortran_write_unit(17,
            " Total force     %14.5e%14.5e%14.5e%14.5e%14.5e%14.5e\n\n",
            (float)cfxtot(ncs + 1), (float)cfytot(ncs + 1), (float)cfztot(ncs + 1),
            (float)cfttot(ncs + 1), (float)cfltot(ncs + 1), (float)cfdtot(ncs + 1));
    }
}

} // namespace csout_ns
