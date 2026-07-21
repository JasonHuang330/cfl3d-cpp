// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "wmag.h"
#include <cstdio>
#include <cmath>

namespace wmag_ns {

void wmag(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q,
          FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
          FortranArray4DRef<double> si, FortranArray3DRef<double> vol,
          FortranArray3DRef<double> vor, FortranArray4DRef<double> w,
          FortranArray3DRef<double> wt, int& ip, FortranArray3DRef<double> blank,
          int& iover, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
          FortranArray4DRef<double> qi0, FortranArray3DRef<double> bcj,
          FortranArray3DRef<double> bck, FortranArray3DRef<double> bci,
          int& nbl, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0,
          FortranArray3DRef<double> voli0, FortranArray1DRef<double> vormax,
          FortranArray1DRef<int> ivmax, FortranArray1DRef<int> jvmax,
          FortranArray1DRef<int> kvmax, int& maxbl)
{
    int jdim1, kdim1, idim1;
    int i, j, k, l;
    int jp1, jm2, kp1, km2, ii, iip1, iim2;
    double term, factor, onec, vmax;

    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;

    onec = 1.;

    //
    // J-direction contributions
    //
    for (i = 1; i <= idim1; i++) {
    for (k = 1; k <= kdim1; k++) {

        // cycle through interfaces
        for (j = 2; j <= jdim1; j++) {
            term      = sj(j,k,i,4) / (vol(j,k,i) + vol(j-1,k,i));
            wt(j,1,1) = term * ( (q(j,k,i,4) - q(j-1,k,i,4)) * sj(j,k,i,2)
                                -(q(j,k,i,3) - q(j-1,k,i,3)) * sj(j,k,i,3) );
            wt(j,1,2) = term * ( (q(j,k,i,2) - q(j-1,k,i,2)) * sj(j,k,i,3)
                                -(q(j,k,i,4) - q(j-1,k,i,4)) * sj(j,k,i,1) );
            wt(j,1,3) = term * ( (q(j,k,i,3) - q(j-1,k,i,3)) * sj(j,k,i,1)
                                -(q(j,k,i,2) - q(j-1,k,i,2)) * sj(j,k,i,2) );
        }

        // contribution for j = 1
        j         = 1;
        jp1       = std::min(jdim-1, j+1);
        term      = sj(j,k,i,4) / (volj0(k,i,1) + vol(j,k,i));
        factor    = bcj(k,i,1) + 1.0;
        term      = term * factor;
        wt(j,1,1) = term * ( (q(j,k,i,4) - qj0(k,i,4,1)) * sj(j,k,i,2)
                             -(q(j,k,i,3) - qj0(k,i,3,1)) * sj(j,k,i,3) );
        wt(j,1,2) = term * ( (q(j,k,i,2) - qj0(k,i,2,1)) * sj(j,k,i,3)
                             -(q(j,k,i,4) - qj0(k,i,4,1)) * sj(j,k,i,1) );
        wt(j,1,3) = term * ( (q(j,k,i,3) - qj0(k,i,3,1)) * sj(j,k,i,1)
                             -(q(j,k,i,2) - qj0(k,i,2,1)) * sj(j,k,i,2) );

        // contribution for j = jdim
        j         = jdim;
        jm2       = std::max(1, j-2);
        term      = sj(j,k,i,4) / (volj0(k,i,3) + vol(j-1,k,i));
        factor    = bcj(k,i,2) + 1.0;
        term      = term * factor;
        wt(j,1,1) = term * ( (qj0(k,i,4,3) - q(j-1,k,i,4)) * sj(j,k,i,2)
                             -(qj0(k,i,3,3) - q(j-1,k,i,3)) * sj(j,k,i,3) );
        wt(j,1,2) = term * ( (qj0(k,i,2,3) - q(j-1,k,i,2)) * sj(j,k,i,3)
                             -(qj0(k,i,4,3) - q(j-1,k,i,4)) * sj(j,k,i,1) );
        wt(j,1,3) = term * ( (qj0(k,i,3,3) - q(j-1,k,i,3)) * sj(j,k,i,1)
                             -(qj0(k,i,2,3) - q(j-1,k,i,2)) * sj(j,k,i,2) );

        // cycle through cell centers
        for (l = 1; l <= 3; l++) {
        for (j = 1; j <= jdim1; j++) {
            w(j,k,i,l) = wt(j,1,l) + wt(j+1,1,l);
        }}

    }} // end i, k loops for J-direction

    //
    // K-direction contributions
    //
    for (i = 1; i <= idim1; i++) {

        // cycle through interfaces
        for (k = 2; k <= kdim-1; k++) {
        for (j = 1; j <= jdim1; j++) {
            term      = sk(j,k,i,4) / (vol(j,k,i) + vol(j,k-1,i));
            wt(j,k,1) = term * ( (q(j,k,i,4) - q(j,k-1,i,4)) * sk(j,k,i,2)
                                -(q(j,k,i,3) - q(j,k-1,i,3)) * sk(j,k,i,3) );
            wt(j,k,2) = term * ( (q(j,k,i,2) - q(j,k-1,i,2)) * sk(j,k,i,3)
                                -(q(j,k,i,4) - q(j,k-1,i,4)) * sk(j,k,i,1) );
            wt(j,k,3) = term * ( (q(j,k,i,3) - q(j,k-1,i,3)) * sk(j,k,i,1)
                                -(q(j,k,i,2) - q(j,k-1,i,2)) * sk(j,k,i,2) );
        }}

        // contribution for k = 1
        k         = 1;
        kp1       = std::min(kdim-1, k+1);
        for (j = 1; j <= jdim1; j++) {
            term      = sk(j,k,i,4) / (volk0(j,i,1) + vol(j,k,i));
            factor    = bck(j,i,1) + 1.0;
            term      = term * factor;
            wt(j,k,1) = term * ( (q(j,k,i,4) - qk0(j,i,4,1)) * sk(j,k,i,2)
                                 -(q(j,k,i,3) - qk0(j,i,3,1)) * sk(j,k,i,3) );
            wt(j,k,2) = term * ( (q(j,k,i,2) - qk0(j,i,2,1)) * sk(j,k,i,3)
                                 -(q(j,k,i,4) - qk0(j,i,4,1)) * sk(j,k,i,1) );
            wt(j,k,3) = term * ( (q(j,k,i,3) - qk0(j,i,3,1)) * sk(j,k,i,1)
                                 -(q(j,k,i,2) - qk0(j,i,2,1)) * sk(j,k,i,2) );
        }

        // contribution for k = kdim
        k         = kdim;
        km2       = std::max(1, k-2);
        for (j = 1; j <= jdim1; j++) {
            term      = sk(j,k,i,4) / (volk0(j,i,3) + vol(j,k-1,i));
            factor    = bck(j,i,2) + 1.0;
            term      = term * factor;
            wt(j,k,1) = term * ( (qk0(j,i,4,3) - q(j,k-1,i,4)) * sk(j,k,i,2)
                                 -(qk0(j,i,3,3) - q(j,k-1,i,3)) * sk(j,k,i,3) );
            wt(j,k,2) = term * ( (qk0(j,i,2,3) - q(j,k-1,i,2)) * sk(j,k,i,3)
                                 -(qk0(j,i,4,3) - q(j,k-1,i,4)) * sk(j,k,i,1) );
            wt(j,k,3) = term * ( (qk0(j,i,3,3) - q(j,k-1,i,3)) * sk(j,k,i,1)
                                 -(qk0(j,i,2,3) - q(j,k-1,i,2)) * sk(j,k,i,2) );
        }

        // cycle through cell centers
        for (l = 1; l <= 3; l++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            w(j,k,i,l) = w(j,k,i,l) + wt(j,k,l) + wt(j,k+1,l);
        }}}

    } // end i loop for K-direction

    //
    // I-direction contributions
    //
    if (idim > 2) {

        for (i = 2; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            term         = si(j,k,i,4) / (vol(j,k,i) + vol(j,k,i-1));
            wt(j,k,1)    = term * ( (q(j,k,i,4) - q(j,k,i-1,4)) * si(j,k,i,2)
                                   -(q(j,k,i,3) - q(j,k,i-1,3)) * si(j,k,i,3) );
            wt(j,k,2)    = term * ( (q(j,k,i,2) - q(j,k,i-1,2)) * si(j,k,i,3)
                                   -(q(j,k,i,4) - q(j,k,i-1,4)) * si(j,k,i,1) );
            wt(j,k,3)    = term * ( (q(j,k,i,3) - q(j,k,i-1,3)) * si(j,k,i,1)
                                   -(q(j,k,i,2) - q(j,k,i-1,2)) * si(j,k,i,2) );
            w(j,k,i-1,1) = w(j,k,i-1,1) + wt(j,k,1);
            w(j,k,i-1,2) = w(j,k,i-1,2) + wt(j,k,2);
            w(j,k,i-1,3) = w(j,k,i-1,3) + wt(j,k,3);
            w(j,k,i,1)   = w(j,k,i,1)   + wt(j,k,1);
            w(j,k,i,2)   = w(j,k,i,2)   + wt(j,k,2);
            w(j,k,i,3)   = w(j,k,i,3)   + wt(j,k,3);
        }}}

        // additional contribution for i = 1
        ii   = 1;
        iip1 = std::min(idim-1, ii+1);
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            term       = si(j,k,ii,4) / (voli0(j,k,1) + vol(j,k,ii));
            factor     = bci(j,k,1) + 1.0;
            term       = term * factor;
            wt(j,k,1)  = term * ( (q(j,k,ii,4) - qi0(j,k,4,1)) * si(j,k,ii,2)
                                  -(q(j,k,ii,3) - qi0(j,k,3,1)) * si(j,k,ii,3) );
            wt(j,k,2)  = term * ( (q(j,k,ii,2) - qi0(j,k,2,1)) * si(j,k,ii,3)
                                  -(q(j,k,ii,4) - qi0(j,k,4,1)) * si(j,k,ii,1) );
            wt(j,k,3)  = term * ( (q(j,k,ii,3) - qi0(j,k,3,1)) * si(j,k,ii,1)
                                  -(q(j,k,ii,2) - qi0(j,k,2,1)) * si(j,k,ii,2) );
            w(j,k,1,1) = w(j,k,1,1) + wt(j,k,1);
            w(j,k,1,2) = w(j,k,1,2) + wt(j,k,2);
            w(j,k,1,3) = w(j,k,1,3) + wt(j,k,3);
        }}

        // additional contribution for i = idim1
        ii   = idim;
        iim2 = std::max(1, ii-2);
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            term      = si(j,k,ii,4) / (voli0(j,k,3) + vol(j,k,ii-1));
            factor    = bci(j,k,2) + 1.0;
            term      = term * factor;
            wt(j,k,1) = term * ( (qi0(j,k,4,3) - q(j,k,ii-1,4)) * si(j,k,ii,2)
                                 -(qi0(j,k,3,3) - q(j,k,ii-1,3)) * si(j,k,ii,3) );
            wt(j,k,2) = term * ( (qi0(j,k,2,3) - q(j,k,ii-1,2)) * si(j,k,ii,3)
                                 -(qi0(j,k,4,3) - q(j,k,ii-1,4)) * si(j,k,ii,1) );
            wt(j,k,3) = term * ( (qi0(j,k,3,3) - q(j,k,ii-1,3)) * si(j,k,ii,1)
                                 -(qi0(j,k,2,3) - q(j,k,ii-1,2)) * si(j,k,ii,2) );
            w(j,k,idim1,1) = w(j,k,idim1,1) + wt(j,k,1);
            w(j,k,idim1,2) = w(j,k,idim1,2) + wt(j,k,2);
            w(j,k,idim1,3) = w(j,k,idim1,3) + wt(j,k,3);
        }}

    } // end if (idim > 2)
    //
    // vorticity magnitude
    //
    for (i = 1; i <= idim1; i++) {
    for (k = 1; k <= kdim1; k++) {
    for (j = 1; j <= jdim1; j++) {
        vor(j,k,i) = std::sqrt( w(j,k,i,1)*w(j,k,i,1)
                               + w(j,k,i,2)*w(j,k,i,2)
                               + w(j,k,i,3)*w(j,k,i,3) );
    }}}

    //
    // set vorticity magnitude for hole cells to one
    //
    if (iover == 1) {
        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            // ccvmgt(onec, vor(j,k,i), (real(blank(j,k,i)).eq.0.e0))
            // returns onec if condition true, vor(j,k,i) otherwise
            vor(j,k,i) = (blank(j,k,i) == 0.0) ? onec : vor(j,k,i);
        }}}
    }

    if (ip > 0) {
        vmax = 0.;
        jvmax(nbl) = 1;
        kvmax(nbl) = 1;
        ivmax(nbl) = 1;
        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            if (vor(j,k,i) > vmax) {
                jvmax(nbl)  = j;
                kvmax(nbl)  = k;
                ivmax(nbl)  = i;
                vormax(nbl) = vor(j,k,i);
                vmax        = vor(j,k,i);
            }
        }}}
    }

    return;
}

} // namespace wmag_ns
