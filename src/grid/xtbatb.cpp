// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "xtbatb.h"

namespace xtbatb_ns {

void xtbatb(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk,
            FortranArray4DRef<double> xtbi, FortranArray4DRef<double> atbj,
            FortranArray4DRef<double> atbk, FortranArray4DRef<double> atbi,
            FortranArray2DRef<double> t,    FortranArray3DRef<double> tti,
            FortranArray3DRef<double> ttj,  FortranArray3DRef<double> ttk)
{
    int kdim1, idim1, jdkd, jdid, idkd;
    int i, j, k, nn;
    int izz, izzk1, izzi1, izzi1k1, izz0, izz1, izzj1, izzi1j1, izzj1k1;

    kdim1 = kdim - 1;
    idim1 = idim - 1;
    jdkd  = jdim * kdim;
    jdid  = jdim * idim;
    idkd  = idim * kdim;

    for (i = 1; i <= idim1; i++) {

        // **********************************************
        // j=1/jdim boundary velocities
        // use face-averaged values
        // **********************************************

        izz1 = (i - 1) * jdkd;
        for (nn = 1; nn <= 2; nn++) {
            j = 1;
            if (nn == 2) j = jdim;
            for (k = 1; k <= kdim1; k++) {
                izz     = izz1 + (k - 1) * jdim + j;
                izzk1   = izz  + jdim;
                izzi1   = izz  + jdkd;
                izzi1k1 = izz  + jdkd + jdim;
                xtbj(k, i, 1, nn) = 0.25 * ( t(izz, 1)   + t(izzk1, 1)  +
                                              t(izzi1, 1) + t(izzi1k1, 1) );
                xtbj(k, i, 2, nn) = 0.25 * ( t(izz, 2)   + t(izzk1, 2)  +
                                              t(izzi1, 2) + t(izzi1k1, 2) );
                xtbj(k, i, 3, nn) = 0.25 * ( t(izz, 3)   + t(izzk1, 3)  +
                                              t(izzi1, 3) + t(izzi1k1, 3) );
            }

            // zero out extra layer of boundary faces for saftey
            xtbj(kdim, i, 1, nn) = 0.;
            xtbj(kdim, i, 2, nn) = 0.;
            xtbj(kdim, i, 3, nn) = 0.;
        }

        // **********************************************
        // j=1/jdim boundary accelerations
        // use face-averaged values
        // **********************************************

        izz0 = (i - 1) * kdim;
        for (nn = 1; nn <= 2; nn++) {
            for (k = 1; k <= kdim - 1; k++) {
                izz     = izz0 + k;
                izzk1   = izz  + 1;
                izzi1   = izz  + kdim;
                izzi1k1 = izz  + kdim + 1;
                atbj(k, i, 1, nn) = 0.25 * ( ttj(izz, 1, nn)   + ttj(izzk1, 1, nn)  +
                                              ttj(izzi1, 1, nn) + ttj(izzi1k1, 1, nn) );
                atbj(k, i, 2, nn) = 0.25 * ( ttj(izz, 2, nn)   + ttj(izzk1, 2, nn)  +
                                              ttj(izzi1, 2, nn) + ttj(izzi1k1, 2, nn) );
                atbj(k, i, 3, nn) = 0.25 * ( ttj(izz, 3, nn)   + ttj(izzk1, 3, nn)  +
                                              ttj(izzi1, 3, nn) + ttj(izzi1k1, 3, nn) );
            }

            // zero out extra layer of boundary faces for saftey
            atbj(kdim, i, 1, nn) = 0.;
            atbj(kdim, i, 2, nn) = 0.;
            atbj(kdim, i, 3, nn) = 0.;
        }

        // **********************************************
        // k=1/kdim boundary velocities
        // use face-averaged values
        // **********************************************

        izz0 = (i - 1) * jdkd;
        for (nn = 1; nn <= 2; nn++) {
            k = 1;
            if (nn == 2) k = kdim;
            for (j = 1; j <= jdim - 1; j++) {
                izz     = izz0 + (k - 1) * jdim + j;
                izzi1   = izz  + jdkd;
                izzj1   = izz  + 1;
                izzi1j1 = izz  + jdkd + 1;
                xtbk(j, i, 1, nn) = 0.25 * ( t(izz, 1)   + t(izzj1, 1)  +
                                              t(izzi1, 1) + t(izzi1j1, 1) );
                xtbk(j, i, 2, nn) = 0.25 * ( t(izz, 2)   + t(izzj1, 2)  +
                                              t(izzi1, 2) + t(izzi1j1, 2) );
                xtbk(j, i, 3, nn) = 0.25 * ( t(izz, 3)   + t(izzj1, 3)  +
                                              t(izzi1, 3) + t(izzi1j1, 3) );
            }

            // zero out extra layer of boundary faces for saftey
            xtbk(jdim, i, 1, nn) = 0.;
            xtbk(jdim, i, 2, nn) = 0.;
            xtbk(jdim, i, 3, nn) = 0.;
        }

        // **********************************************
        // k=1/kdim boundary accelerations
        // use face-averaged values
        // **********************************************

        for (nn = 1; nn <= 2; nn++) {
            for (j = 1; j <= jdim - 1; j++) {
                izz     = (i - 1) * jdim + j;
                izzj1   = izz + 1;
                izzi1   = izz + jdim;
                izzi1j1 = izz + jdim + 1;
                atbk(j, i, 1, nn) = 0.25 * ( ttk(izz, 1, nn)   + ttk(izzj1, 1, nn)  +
                                              ttk(izzi1, 1, nn) + ttk(izzi1j1, 1, nn) );
                atbk(j, i, 2, nn) = 0.25 * ( ttk(izz, 2, nn)   + ttk(izzj1, 2, nn)  +
                                              ttk(izzi1, 2, nn) + ttk(izzi1j1, 2, nn) );
                atbk(j, i, 3, nn) = 0.25 * ( ttk(izz, 3, nn)   + ttk(izzj1, 3, nn)  +
                                              ttk(izzi1, 3, nn) + ttk(izzi1j1, 3, nn) );
            }

            // zero out extra layer of boundary faces for saftey
            atbk(jdim, i, 1, nn) = 0.;
            atbk(jdim, i, 2, nn) = 0.;
            atbk(jdim, i, 3, nn) = 0.;
        }

    } // end loop i=1..idim1

    // **********************************************
    // i=1/idim boundary velocities
    // use face-averaged values
    // **********************************************

    for (nn = 1; nn <= 2; nn++) {
        i = 1;
        if (nn == 2) i = idim;
        izz0 = (i - 1) * jdkd;
        for (j = 1; j <= jdim - 1; j++) {
            for (k = 1; k <= kdim - 1; k++) {
                izz     = izz0 + (k - 1) * jdim + j;
                izzj1   = izz + 1;
                izzk1   = izz + jdim;
                izzj1k1 = izz + 1 + jdim;
                xtbi(j, k, 1, nn) = 0.25 * ( t(izz, 1)   + t(izzj1, 1)  +
                                              t(izzk1, 1) + t(izzj1k1, 1) );
                xtbi(j, k, 2, nn) = 0.25 * ( t(izz, 2)   + t(izzj1, 2)  +
                                              t(izzk1, 2) + t(izzj1k1, 2) );
                xtbi(j, k, 3, nn) = 0.25 * ( t(izz, 3)   + t(izzj1, 3)  +
                                              t(izzk1, 3) + t(izzj1k1, 3) );
            }
        }

        // zero out extra layer of boundary faces for saftey
        for (j = 1; j <= jdim; j++) {
            xtbi(j, kdim, 1, nn) = 0.;
            xtbi(j, kdim, 2, nn) = 0.;
            xtbi(j, kdim, 3, nn) = 0.;
        }

        // zero out extra layer of boundary faces for saftey
        for (k = 1; k <= kdim; k++) {
            xtbi(jdim, k, 1, nn) = 0.;
            xtbi(jdim, k, 2, nn) = 0.;
            xtbi(jdim, k, 3, nn) = 0.;
        }
    }

    // **********************************************
    // i=1/idim boundary accelerations
    // use face-averaged values
    // **********************************************

    for (nn = 1; nn <= 2; nn++) {
        for (j = 1; j <= jdim - 1; j++) {
            for (k = 1; k <= kdim - 1; k++) {
                izz     = (k - 1) * jdim + j;
                izzj1   = izz + 1;
                izzk1   = izz + jdim;
                izzj1k1 = izz + jdim + 1;
                atbi(j, k, 1, nn) = 0.25 * ( tti(izz, 1, nn)   + tti(izzj1, 1, nn)  +
                                              tti(izzk1, 1, nn) + tti(izzj1k1, 1, nn) );
                atbi(j, k, 2, nn) = 0.25 * ( tti(izz, 2, nn)   + tti(izzj1, 2, nn)  +
                                              tti(izzk1, 2, nn) + tti(izzj1k1, 2, nn) );
                atbi(j, k, 3, nn) = 0.25 * ( tti(izz, 3, nn)   + tti(izzj1, 3, nn)  +
                                              tti(izzk1, 3, nn) + tti(izzj1k1, 3, nn) );
            }
        }

        // zero out extra layer of boundary faces for saftey
        for (k = 1; k <= kdim; k++) {
            atbi(jdim, k, 1, nn) = 0.;
            atbi(jdim, k, 2, nn) = 0.;
            atbi(jdim, k, 3, nn) = 0.;
        }

        // zero out extra layer of boundary faces for saftey
        for (j = 1; j <= jdim; j++) {
            atbi(j, kdim, 1, nn) = 0.;
            atbi(j, kdim, 2, nn) = 0.;
            atbi(j, kdim, 3, nn) = 0.;
        }
    }

    return;
}

} // namespace xtbatb_ns
