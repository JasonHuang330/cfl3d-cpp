// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc1011.h"
#include "runtime/common_blocks.h"
#include "runtime/fortran_array.h"
#include <algorithm>
#include <cstring>
#include <cstdio>

namespace bc1011_ns {

void bc1011(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0, FortranArray4DRef<double> sj,
            FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci, FortranArray4DRef<double> xtbj,
            FortranArray4DRef<double> xtbk, FortranArray4DRef<double> xtbi,
            FortranArray4DRef<double> atbj, FortranArray4DRef<double> atbk,
            FortranArray4DRef<double> atbi, int& ista, int& iend,
            int& jsta, int& jend, int& ksta, int& kend, int& nface,
            FortranArray4DRef<double> tursav, FortranArray4DRef<double> tj0,
            FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> vist3d, FortranArray4DRef<double> vj0,
            FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0,
            int& isym, int& jsym, int& ksym, int& iuns,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& nummem)
{
    // COMMON block references
    int& level   = cmn_mgrd.level;
    int& lglobal = cmn_mgrd.lglobal;
    int* ivisc   = cmn_reyue.ivisc;  // 0-based: ivisc[0]=ivisc(1), [1]=ivisc(2), [2]=ivisc(3)

    // Local variables
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    double vcont1, vcont2, wcont1, wcont2, ucont1, ucont2;
    int i2, l, i, j, k;



    //******************************************************************************
    //     j=1 boundary            singular axis - half plane           bctype 1011
    //******************************************************************************

    if (nface == 3) {

        // symmetry in k
        if (ksym > 0) {
            for (i = ista; i <= iend1; i++) {
            for (k = ksta; k <= kend1; k++) {

                vcont1 = q(1,kend-k,i,2)*sk(1,kend,i,1) +
                         q(1,kend-k,i,3)*sk(1,kend,i,2) +
                         q(1,kend-k,i,4)*sk(1,kend,i,3) + sk(1,kend,i,5);
                vcont2 = q(2,kend-k,i,2)*sk(2,kend,i,1) +
                         q(2,kend-k,i,3)*sk(2,kend,i,2) +
                         q(2,kend-k,i,4)*sk(2,kend,i,3) + sk(2,kend,i,5);

                qj0(k,i,1,1) = q(1,kend-k,i,1);
                qj0(k,i,2,1) = q(1,kend-k,i,2) - 2.*vcont1*sk(1,kend,i,1);
                qj0(k,i,3,1) = q(1,kend-k,i,3) - 2.*vcont1*sk(1,kend,i,2);
                qj0(k,i,4,1) = q(1,kend-k,i,4) - 2.*vcont1*sk(1,kend,i,3);
                qj0(k,i,5,1) = q(1,kend-k,i,5);

                qj0(k,i,1,2) = q(2,kend-k,i,1);
                qj0(k,i,2,2) = q(2,kend-k,i,2) - 2.*vcont2*sk(2,kend,i,1);
                qj0(k,i,3,2) = q(2,kend-k,i,3) - 2.*vcont2*sk(2,kend,i,2);
                qj0(k,i,4,2) = q(2,kend-k,i,4) - 2.*vcont2*sk(2,kend,i,3);
                qj0(k,i,5,2) = q(2,kend-k,i,5);

                bcj(k,i,1) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (i = ista; i <= iend1; i++) {
                for (k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,1) = vist3d(1,kend-k,i);
                    vj0(k,i,1,2) = vist3d(2,kend-k,i);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (i = ista; i <= iend1; i++) {
                    for (k = ksta; k <= kend1; k++) {
                        tj0(k,i,l,1) = tursav(1,kend-k,i,l);
                        tj0(k,i,l,2) = tursav(2,kend-k,i,l);
                    }}}
                }
            }

        } // end ksym > 0

        // symmetry in i
        if (isym > 0) {
            for (i = ista; i <= iend1; i++) {
            for (k = ksta; k <= kend1; k++) {

                vcont1 = q(1,k,iend-i,2)*si(1,k,iend,1) +
                         q(1,k,iend-i,3)*si(1,k,iend,2) +
                         q(1,k,iend-i,4)*si(1,k,iend,3) + si(1,k,iend,5);
                vcont2 = q(2,k,iend-i,2)*si(2,k,iend,1) +
                         q(2,k,iend-i,3)*si(2,k,iend,2) +
                         q(2,k,iend-i,4)*si(2,k,iend,3) + si(2,k,iend,5);

                qj0(k,i,1,1) = q(1,k,iend-i,1);
                qj0(k,i,2,1) = q(1,k,iend-i,2) - 2.*vcont1*si(1,k,iend,1);
                qj0(k,i,3,1) = q(1,k,iend-i,3) - 2.*vcont1*si(1,k,iend,2);
                qj0(k,i,4,1) = q(1,k,iend-i,4) - 2.*vcont1*si(1,k,iend,3);
                qj0(k,i,5,1) = q(1,k,iend-i,5);

                qj0(k,i,1,2) = q(2,k,iend-i,1);
                qj0(k,i,2,2) = q(2,k,iend-i,2) - 2.*vcont2*si(2,k,iend,1);
                qj0(k,i,3,2) = q(2,k,iend-i,3) - 2.*vcont2*si(2,k,iend,2);
                qj0(k,i,4,2) = q(2,k,iend-i,4) - 2.*vcont2*si(2,k,iend,3);
                qj0(k,i,5,2) = q(2,k,iend-i,5);

                bcj(k,i,1) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (i = ista; i <= iend1; i++) {
                for (k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,1) = vist3d(1,k,iend-i);
                    vj0(k,i,1,2) = vist3d(2,k,iend-i);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (i = ista; i <= iend1; i++) {
                    for (k = ksta; k <= kend1; k++) {
                        tj0(k,i,l,1) = tursav(1,k,iend-i,l);
                        tj0(k,i,l,2) = tursav(2,k,iend-i,l);
                    }}}
                }
            }

        } // end isym > 0

    } // end nface == 3



    //******************************************************************************
    //     j=jdim boundary         singular axis - half plane           bctype 1011
    //******************************************************************************

    if (nface == 4) {

        // symmetry in k
        if (ksym > 0) {
            for (i = ista; i <= iend1; i++) {
            for (k = ksta; k <= kend1; k++) {

                vcont1 = q(jdim-1,kend-k,i,2)*sk(jdim-1,kend,i,1) +
                         q(jdim-1,kend-k,i,3)*sk(jdim-1,kend,i,2) +
                         q(jdim-1,kend-k,i,4)*sk(jdim-1,kend,i,3) +
                                              sk(jdim-1,kend,i,5);
                vcont2 = q(jdim-2,kend-k,i,2)*sk(jdim-2,kend,i,1) +
                         q(jdim-2,kend-k,i,3)*sk(jdim-2,kend,i,2) +
                         q(jdim-2,kend-k,i,4)*sk(jdim-2,kend,i,3) +
                                              sk(jdim-2,kend,i,5);

                qj0(k,i,1,3) = q(jdim-1,kend-k,i,1);
                qj0(k,i,2,3) = q(jdim-1,kend-k,i,2) - 2.*vcont1*sk(jdim-1,kend,i,1);
                qj0(k,i,3,3) = q(jdim-1,kend-k,i,3) - 2.*vcont1*sk(jdim-1,kend,i,2);
                qj0(k,i,4,3) = q(jdim-1,kend-k,i,4) - 2.*vcont1*sk(jdim-1,kend,i,3);
                qj0(k,i,5,3) = q(jdim-1,kend-k,i,5);

                qj0(k,i,1,4) = q(jdim-2,kend-k,i,1);
                qj0(k,i,2,4) = q(jdim-2,kend-k,i,2) - 2.*vcont2*sk(jdim-2,kend,i,1);
                qj0(k,i,3,4) = q(jdim-2,kend-k,i,3) - 2.*vcont2*sk(jdim-2,kend,i,2);
                qj0(k,i,4,4) = q(jdim-2,kend-k,i,4) - 2.*vcont2*sk(jdim-2,kend,i,3);
                qj0(k,i,5,4) = q(jdim-2,kend-k,i,5);

                bcj(k,i,2) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (i = ista; i <= iend1; i++) {
                for (k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,3) = vist3d(jdim-1,kend-k,i);
                    vj0(k,i,1,4) = vist3d(jdim-2,kend-k,i);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (i = ista; i <= iend1; i++) {
                    for (k = ksta; k <= kend1; k++) {
                        tj0(k,i,l,3) = tursav(jdim-1,kend-k,i,l);
                        tj0(k,i,l,4) = tursav(jdim-2,kend-k,i,l);
                    }}}
                }
            }

        } // end ksym > 0

        // symmetry in i
        if (isym > 0) {
            for (i = ista; i <= iend1; i++) {
            for (k = ksta; k <= kend1; k++) {

                vcont1 = q(jdim-1,k,iend-i,2)*si(jdim-1,k,iend,1) +
                         q(jdim-1,k,iend-i,3)*si(jdim-1,k,iend,2) +
                         q(jdim-1,k,iend-i,4)*si(jdim-1,k,iend,3) +
                                              si(jdim-1,k,iend,5);
                vcont2 = q(jdim-2,k,iend-i,2)*si(jdim-2,k,iend,1) +
                         q(jdim-2,k,iend-i,3)*si(jdim-2,k,iend,2) +
                         q(jdim-2,k,iend-i,4)*si(jdim-2,k,iend,3) +
                                              si(jdim-2,k,iend,5);

                qj0(k,i,1,3) = q(jdim-1,k,iend-i,1);
                qj0(k,i,2,3) = q(jdim-1,k,iend-i,2) - 2.*vcont1*si(jdim-1,k,iend,1);
                qj0(k,i,3,3) = q(jdim-1,k,iend-i,3) - 2.*vcont1*si(jdim-1,k,iend,2);
                qj0(k,i,4,3) = q(jdim-1,k,iend-i,4) - 2.*vcont1*si(jdim-1,k,iend,3);
                qj0(k,i,5,3) = q(jdim-1,k,iend-i,5);

                qj0(k,i,1,4) = q(jdim-2,k,iend-i,1);
                qj0(k,i,2,4) = q(jdim-2,k,iend-i,2) - 2.*vcont2*si(jdim-2,k,iend,1);
                qj0(k,i,3,4) = q(jdim-2,k,iend-i,3) - 2.*vcont2*si(jdim-2,k,iend,2);
                qj0(k,i,4,4) = q(jdim-2,k,iend-i,4) - 2.*vcont2*si(jdim-2,k,iend,3);
                qj0(k,i,5,4) = q(jdim-2,k,iend-i,5);

                bcj(k,i,2) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (i = ista; i <= iend1; i++) {
                for (k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,3) = vist3d(jdim-1,k,iend-i);
                    vj0(k,i,1,4) = vist3d(jdim-2,k,iend-i);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (i = ista; i <= iend1; i++) {
                    for (k = ksta; k <= kend1; k++) {
                        tj0(k,i,l,3) = tursav(jdim-1,k,iend-i,l);
                        tj0(k,i,l,4) = tursav(jdim-2,k,iend-i,l);
                    }}}
                }
            }

        } // end isym > 0

    } // end nface == 4



    //******************************************************************************
    //     k=1 boundary            singular axis - half plane           bctype 1011
    //******************************************************************************

    if (nface == 5) {

        // symmetry in j
        if (jsym > 0) {
            for (i = ista; i <= iend1; i++) {
            for (j = jsta; j <= jend1; j++) {

                wcont1 = q(jend-j,1,i,2)*sj(jend,1,i,1) +
                         q(jend-j,1,i,3)*sj(jend,1,i,2) +
                         q(jend-j,1,i,4)*sj(jend,1,i,3) + sj(jend,1,i,5);
                wcont2 = q(jend-j,2,i,2)*sj(jend,2,i,1) +
                         q(jend-j,2,i,3)*sj(jend,2,i,2) +
                         q(jend-j,2,i,4)*sj(jend,2,i,3) + sj(jend,2,i,5);

                qk0(j,i,1,1) = q(jend-j,1,i,1);
                qk0(j,i,2,1) = q(jend-j,1,i,2) - 2.*wcont1*sj(jend,1,i,1);
                qk0(j,i,3,1) = q(jend-j,1,i,3) - 2.*wcont1*sj(jend,1,i,2);
                qk0(j,i,4,1) = q(jend-j,1,i,4) - 2.*wcont1*sj(jend,1,i,3);
                qk0(j,i,5,1) = q(jend-j,1,i,5);

                qk0(j,i,1,2) = q(jend-j,2,i,1);
                qk0(j,i,2,2) = q(jend-j,2,i,2) - 2.*wcont2*sj(jend,2,i,1);
                qk0(j,i,3,2) = q(jend-j,2,i,3) - 2.*wcont2*sj(jend,2,i,2);
                qk0(j,i,4,2) = q(jend-j,2,i,4) - 2.*wcont2*sj(jend,2,i,3);
                qk0(j,i,5,2) = q(jend-j,2,i,5);

                bck(j,i,1) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (i = ista; i <= iend1; i++) {
                for (j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,1) = vist3d(jend-j,1,i);
                    vk0(j,i,1,2) = vist3d(jend-j,2,i);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (i = ista; i <= iend1; i++) {
                    for (j = jsta; j <= jend1; j++) {
                        tk0(j,i,l,1) = tursav(jend-j,1,i,l);
                        tk0(j,i,l,2) = tursav(jend-j,2,i,l);
                    }}}
                }
            }

        } // end jsym > 0

        // symmetry in i
        if (isym > 0) {
            for (i = ista; i <= iend1; i++) {
            for (j = jsta; j <= jend1; j++) {

                wcont1 = q(j,1,iend-i,2)*si(j,1,iend,1) +
                         q(j,1,iend-i,3)*si(j,1,iend,2) +
                         q(j,1,iend-i,4)*si(j,1,iend,3) + si(j,1,iend,5);
                wcont2 = q(j,2,iend-i,2)*si(j,2,iend,1) +
                         q(j,2,iend-i,3)*si(j,2,iend,2) +
                         q(j,2,iend-i,4)*si(j,2,iend,3) + si(j,2,iend,5);

                qk0(j,i,1,1) = q(j,1,iend-i,1);
                qk0(j,i,2,1) = q(j,1,iend-i,2) - 2.*wcont1*si(j,1,iend,1);
                qk0(j,i,3,1) = q(j,1,iend-i,3) - 2.*wcont1*si(j,1,iend,2);
                qk0(j,i,4,1) = q(j,1,iend-i,4) - 2.*wcont1*si(j,1,iend,3);
                qk0(j,i,5,1) = q(j,1,iend-i,5);

                qk0(j,i,1,2) = q(j,2,iend-i,1);
                qk0(j,i,2,2) = q(j,2,iend-i,2) - 2.*wcont2*si(j,2,iend,1);
                qk0(j,i,3,2) = q(j,2,iend-i,3) - 2.*wcont2*si(j,2,iend,2);
                qk0(j,i,4,2) = q(j,2,iend-i,4) - 2.*wcont2*si(j,2,iend,3);
                qk0(j,i,5,2) = q(j,2,iend-i,5);

                bck(j,i,1) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (i = ista; i <= iend1; i++) {
                for (j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,1) = vist3d(j,1,iend-i);
                    vk0(j,i,1,2) = vist3d(j,2,iend-i);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (i = ista; i <= iend1; i++) {
                    for (j = jsta; j <= jend1; j++) {
                        tk0(j,i,l,1) = tursav(j,1,iend-i,l);
                        tk0(j,i,l,2) = tursav(j,2,iend-i,l);
                    }}}
                }
            }

        } // end isym > 0

    } // end nface == 5



    //******************************************************************************
    //     k=kdim boundary         singular axis - half plane           bctype 1011
    //******************************************************************************

    if (nface == 6) {

        // symmetry in j
        if (jsym > 0) {
            for (i = ista; i <= iend1; i++) {
            for (j = jsta; j <= jend1; j++) {

                wcont1 = q(jend-j,kdim-1,i,2)*sj(jend,kdim-1,i,1) +
                         q(jend-j,kdim-1,i,3)*sj(jend,kdim-1,i,2) +
                         q(jend-j,kdim-1,i,4)*sj(jend,kdim-1,i,3) +
                                              sj(jend,kdim-1,i,5);
                wcont2 = q(jend-j,kdim-2,i,2)*sj(jend,kdim-2,i,1) +
                         q(jend-j,kdim-2,i,3)*sj(jend,kdim-2,i,2) +
                         q(jend-j,kdim-2,i,4)*sj(jend,kdim-2,i,3) +
                                              sj(jend,kdim-2,i,5);

                qk0(j,i,1,3) = q(jend-j,kdim-1,i,1);
                qk0(j,i,2,3) = q(jend-j,kdim-1,i,2) - 2.*wcont1*sj(jend,kdim-1,i,1);
                qk0(j,i,3,3) = q(jend-j,kdim-1,i,3) - 2.*wcont1*sj(jend,kdim-1,i,2);
                qk0(j,i,4,3) = q(jend-j,kdim-1,i,4) - 2.*wcont1*sj(jend,kdim-1,i,3);
                qk0(j,i,5,3) = q(jend-j,kdim-1,i,5);

                qk0(j,i,1,4) = q(jend-j,kdim-2,i,1);
                qk0(j,i,2,4) = q(jend-j,kdim-2,i,2) - 2.*wcont2*sj(jend,kdim-2,i,1);
                qk0(j,i,3,4) = q(jend-j,kdim-2,i,3) - 2.*wcont2*sj(jend,kdim-2,i,2);
                qk0(j,i,4,4) = q(jend-j,kdim-2,i,4) - 2.*wcont2*sj(jend,kdim-2,i,3);
                qk0(j,i,5,4) = q(jend-j,kdim-2,i,5);

                bck(j,i,2) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (i = ista; i <= iend1; i++) {
                for (j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,3) = vist3d(jend-j,kdim-1,i);
                    vk0(j,i,1,4) = vist3d(jend-j,kdim-2,i);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (i = ista; i <= iend1; i++) {
                    for (j = jsta; j <= jend1; j++) {
                        tk0(j,i,l,3) = tursav(jend-j,kdim-1,i,l);
                        tk0(j,i,l,4) = tursav(jend-j,kdim-2,i,l);
                    }}}
                }
            }

        } // end jsym > 0

        // symmetry in i
        if (isym > 0) {
            for (i = ista; i <= iend1; i++) {
            for (j = jsta; j <= jend1; j++) {

                wcont1 = q(j,kdim-1,iend-i,2)*si(j,kdim-1,iend,1) +
                         q(j,kdim-1,iend-i,3)*si(j,kdim-1,iend,2) +
                         q(j,kdim-1,iend-i,4)*si(j,kdim-1,iend,3) +
                                              si(j,kdim-1,iend,5);
                wcont2 = q(j,kdim-2,iend-i,2)*si(j,kdim-2,iend,1) +
                         q(j,kdim-2,iend-i,3)*si(j,kdim-2,iend,2) +
                         q(j,kdim-2,iend-i,4)*si(j,kdim-2,iend,3) +
                                              si(j,kdim-2,iend,5);

                qk0(j,i,1,3) = q(j,kdim-1,iend-i,1);
                qk0(j,i,2,3) = q(j,kdim-1,iend-i,2) - 2.*wcont1*si(j,kdim-1,iend,1);
                qk0(j,i,3,3) = q(j,kdim-1,iend-i,3) - 2.*wcont1*si(j,kdim-1,iend,2);
                qk0(j,i,4,3) = q(j,kdim-1,iend-i,4) - 2.*wcont1*si(j,kdim-1,iend,3);
                qk0(j,i,5,3) = q(j,kdim-1,iend-i,5);

                qk0(j,i,1,4) = q(j,kdim-2,iend-i,1);
                qk0(j,i,2,4) = q(j,kdim-2,iend-i,2) - 2.*wcont2*si(j,kdim-2,iend,1);
                qk0(j,i,3,4) = q(j,kdim-2,iend-i,3) - 2.*wcont2*si(j,kdim-2,iend,2);
                qk0(j,i,4,4) = q(j,kdim-2,iend-i,4) - 2.*wcont2*si(j,kdim-2,iend,3);
                qk0(j,i,5,4) = q(j,kdim-2,iend-i,5);

                bck(j,i,2) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (i = ista; i <= iend1; i++) {
                for (j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,3) = vist3d(j,kdim-1,iend-i);
                    vk0(j,i,1,4) = vist3d(j,kdim-2,iend-i);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (i = ista; i <= iend1; i++) {
                    for (j = jsta; j <= jend1; j++) {
                        tk0(j,i,l,3) = tursav(j,kdim-1,iend-i,l);
                        tk0(j,i,l,4) = tursav(j,kdim-2,iend-i,l);
                    }}}
                }
            }

        } // end isym > 0

    } // end nface == 6



    //******************************************************************************
    //     i=1 boundary            singular axis - half plane           bctype 1011
    //******************************************************************************

    if (nface == 1) {

        i2 = std::min(2, idim1);

        // symmetry in j
        if (jsym > 0) {
            for (k = ksta; k <= kend1; k++) {
            for (j = jsta; j <= jend1; j++) {

                ucont1 = q(jend-j,k,1,2)*sj(jend,k,1,1) +
                         q(jend-j,k,1,3)*sj(jend,k,1,2) +
                         q(jend-j,k,1,4)*sj(jend,k,1,3) + sj(jend,k,1,5);
                ucont2 = q(jend-j,k,i2,2)*sj(jend,k,i2,1) +
                         q(jend-j,k,i2,3)*sj(jend,k,i2,2) +
                         q(jend-j,k,i2,4)*sj(jend,k,i2,3) + sj(jend,k,i2,5);

                qi0(j,k,1,1) = q(jend-j,k,1,1);
                qi0(j,k,2,1) = q(jend-j,k,1,2) - 2.*ucont1*sj(jend,k,1,1);
                qi0(j,k,3,1) = q(jend-j,k,1,3) - 2.*ucont1*sj(jend,k,1,2);
                qi0(j,k,4,1) = q(jend-j,k,1,4) - 2.*ucont1*sj(jend,k,1,3);
                qi0(j,k,5,1) = q(jend-j,k,1,5);

                qi0(j,k,1,2) = q(jend-j,k,i2,1);
                qi0(j,k,2,2) = q(jend-j,k,i2,2) - 2.*ucont2*sj(jend,k,i2,1);
                qi0(j,k,3,2) = q(jend-j,k,i2,3) - 2.*ucont2*sj(jend,k,i2,2);
                qi0(j,k,4,2) = q(jend-j,k,i2,4) - 2.*ucont2*sj(jend,k,i2,3);
                qi0(j,k,5,2) = q(jend-j,k,i2,5);

                bci(j,k,1) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (k = ksta; k <= kend1; k++) {
                for (j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,1) = vist3d(jend-j,k,1);
                    vi0(j,k,1,2) = vist3d(jend-j,k,i2);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (k = ksta; k <= kend1; k++) {
                    for (j = jsta; j <= jend1; j++) {
                        ti0(j,k,l,1) = tursav(jend-j,k,1,l);
                        ti0(j,k,l,2) = tursav(jend-j,k,i2,l);
                    }}}
                }
            }

        } // end jsym > 0

        // symmetry in k
        if (ksym > 0) {
            for (k = ksta; k <= kend1; k++) {
            for (j = jsta; j <= jend1; j++) {

                ucont1 = q(j,kend-k,1,2)*sk(j,kend,1,1) +
                         q(j,kend-k,1,3)*sk(j,kend,1,2) +
                         q(j,kend-k,1,4)*sk(j,kend,1,3) + sk(j,kend,1,5);
                ucont2 = q(j,kend-k,i2,2)*sk(j,kend,i2,1) +
                         q(j,kend-k,i2,3)*sk(j,kend,i2,2) +
                         q(j,kend-k,i2,4)*sk(j,kend,i2,3) + sk(j,kend,i2,5);

                qi0(j,k,1,1) = q(j,kend-k,1,1);
                qi0(j,k,2,1) = q(j,kend-k,1,2) - 2.*ucont1*sk(j,kend,1,1);
                qi0(j,k,3,1) = q(j,kend-k,1,3) - 2.*ucont1*sk(j,kend,1,2);
                qi0(j,k,4,1) = q(j,kend-k,1,4) - 2.*ucont1*sk(j,kend,1,3);
                qi0(j,k,5,1) = q(j,kend-k,1,5);

                qi0(j,k,1,2) = q(j,kend-k,i2,1);
                qi0(j,k,2,2) = q(j,kend-k,i2,2) - 2.*ucont2*sk(j,kend,i2,1);
                qi0(j,k,3,2) = q(j,kend-k,i2,3) - 2.*ucont2*sk(j,kend,i2,2);
                qi0(j,k,4,2) = q(j,kend-k,i2,4) - 2.*ucont2*sk(j,kend,i2,3);
                qi0(j,k,5,2) = q(j,kend-k,i2,5);

                bci(j,k,1) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (k = ksta; k <= kend1; k++) {
                for (j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,1) = vist3d(j,kend-k,1);
                    vi0(j,k,1,2) = vist3d(j,kend-k,i2);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (k = ksta; k <= kend1; k++) {
                    for (j = jsta; j <= jend1; j++) {
                        ti0(j,k,l,1) = tursav(j,kend-k,1,l);
                        ti0(j,k,l,2) = tursav(j,kend-k,i2,l);
                    }}}
                }
            }

        } // end ksym > 0

    } // end nface == 1



    //******************************************************************************
    //     i=idim boundary         singular axis - half plane           bctype 1011
    //******************************************************************************

    if (nface == 2) {

        i2 = std::max(1, idim - 2);

        // symmetry in j
        if (jsym > 0) {
            for (k = ksta; k <= kend1; k++) {
            for (j = jsta; j <= jend1; j++) {

                ucont1 = q(jend-j,k,idim-1,2)*sj(jend,k,idim-1,1) +
                         q(jend-j,k,idim-1,3)*sj(jend,k,idim-1,2) +
                         q(jend-j,k,idim-1,4)*sj(jend,k,idim-1,3) +
                                              sj(jend,k,idim-1,5);
                ucont2 = q(jend-j,k,i2,2)*sj(jend,k,i2,1) +
                         q(jend-j,k,i2,3)*sj(jend,k,i2,2) +
                         q(jend-j,k,i2,4)*sj(jend,k,i2,3) +
                                         sj(jend,k,i2,5);

                qi0(j,k,1,3) = q(jend-j,k,idim-1,1);
                qi0(j,k,2,3) = q(jend-j,k,idim-1,2) - 2.*ucont1*sj(jend,k,idim-1,1);
                qi0(j,k,3,3) = q(jend-j,k,idim-1,3) - 2.*ucont1*sj(jend,k,idim-1,2);
                qi0(j,k,4,3) = q(jend-j,k,idim-1,4) - 2.*ucont1*sj(jend,k,idim-1,3);
                qi0(j,k,5,3) = q(jend-j,k,idim-1,5);

                qi0(j,k,1,4) = q(jend-j,k,i2,1);
                qi0(j,k,2,4) = q(jend-j,k,i2,2) - 2.*ucont2*sj(jend,k,i2,1);
                qi0(j,k,3,4) = q(jend-j,k,i2,3) - 2.*ucont2*sj(jend,k,i2,2);
                qi0(j,k,4,4) = q(jend-j,k,i2,4) - 2.*ucont2*sj(jend,k,i2,3);
                qi0(j,k,5,4) = q(jend-j,k,i2,5);

                bci(j,k,2) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (k = ksta; k <= kend1; k++) {
                for (j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,3) = vist3d(jend-j,k,idim-1);
                    vi0(j,k,1,4) = vist3d(jend-j,k,i2);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (k = ksta; k <= kend1; k++) {
                    for (j = jsta; j <= jend1; j++) {
                        ti0(j,k,l,3) = tursav(jend-j,k,idim-1,l);
                        ti0(j,k,l,4) = tursav(jend-j,k,i2,l);
                    }}}
                }
            }

        } // end jsym > 0

        // symmetry in k
        if (ksym > 0) {
            for (k = ksta; k <= kend1; k++) {
            for (j = jsta; j <= jend1; j++) {

                ucont1 = q(j,kend-k,idim-1,2)*sk(j,kend,idim-1,1) +
                         q(j,kend-k,idim-1,3)*sk(j,kend,idim-1,2) +
                         q(j,kend-k,idim-1,4)*sk(j,kend,idim-1,3) +
                                              sk(j,kend,idim-1,5);
                ucont2 = q(j,kend-k,i2,2)*sk(j,kend,i2,1) +
                         q(j,kend-k,i2,3)*sk(j,kend,i2,2) +
                         q(j,kend-k,i2,4)*sk(j,kend,i2,3) +
                                         sk(j,kend,i2,5);

                qi0(j,k,1,3) = q(j,kend-k,idim-1,1);
                qi0(j,k,2,3) = q(j,kend-k,idim-1,2) - 2.*ucont1*sk(j,kend,idim-1,1);
                qi0(j,k,3,3) = q(j,kend-k,idim-1,3) - 2.*ucont1*sk(j,kend,idim-1,2);
                qi0(j,k,4,3) = q(j,kend-k,idim-1,4) - 2.*ucont1*sk(j,kend,idim-1,3);
                qi0(j,k,5,3) = q(j,kend-k,idim-1,5);

                qi0(j,k,1,4) = q(j,kend-k,i2,1);
                qi0(j,k,2,4) = q(j,kend-k,i2,2) - 2.*ucont2*sk(j,kend,i2,1);
                qi0(j,k,3,4) = q(j,kend-k,i2,3) - 2.*ucont2*sk(j,kend,i2,2);
                qi0(j,k,4,4) = q(j,kend-k,i2,4) - 2.*ucont2*sk(j,kend,i2,3);
                qi0(j,k,5,4) = q(j,kend-k,i2,5);

                bci(j,k,2) = 0.0;

            }}

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (k = ksta; k <= kend1; k++) {
                for (j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,3) = vist3d(j,kend-k,idim-1);
                    vi0(j,k,1,4) = vist3d(j,kend-k,i2);
                }}
            }
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (l = 1; l <= nummem; l++) {
                    for (k = ksta; k <= kend1; k++) {
                    for (j = jsta; j <= jend1; j++) {
                        ti0(j,k,l,3) = tursav(j,kend-k,idim-1,l);
                        ti0(j,k,l,4) = tursav(j,kend-k,i2,l);
                    }}}
                }
            }

        } // end ksym > 0

    } // end nface == 2

    return;
} // end bc1011



void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
        int& nwork, double& cl, FortranArray1DRef<int> nou,
        FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
        int& maxbl, int& maxgr, int& maxseg,
        FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
        FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg,
        FortranArray1DRef<int> nblg, FortranArray1DRef<int> nbci0,
        FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
        FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim,
        FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo,
        FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
        FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej,
        FortranArray3DRef<int> bcfilek, FortranArray3DRef<int> lwdat,
        int& myid, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
        FortranArray1DRef<int> kdimg, FortranArray1DRef<char[80]> bcfiles,
        int& mxbcfil, int& nummem)
{
    int iuns = 0;

    for (int ibl = 1; ibl <= nbl; ibl++) {
        int jdim = jdimg(ibl);
        int kdim = kdimg(ibl);
        int idim = idimg(ibl);

        // Build array views from w using lw offsets
        FortranArray4DRef<double> q    (&w(lw(1,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray4DRef<double> qj0  (&w(lw(2,ibl)),  kdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qk0  (&w(lw(3,ibl)),  jdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qi0  (&w(lw(4,ibl)),  jdim,   kdim,   5,      4);
        FortranArray4DRef<double> sj   (&w(lw(5,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> sk   (&w(lw(6,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> si   (&w(lw(7,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray3DRef<double> bcj  (&w(lw(8,ibl)),  kdim,   idim-1, 2);
        FortranArray3DRef<double> bck  (&w(lw(9,ibl)),  jdim,   idim-1, 2);
        FortranArray3DRef<double> bci  (&w(lw(10,ibl)), jdim,   kdim,   2);
        FortranArray4DRef<double> xtbj (&w(lw(11,ibl)), kdim,   idim-1, 3,      2);
        FortranArray4DRef<double> xtbk (&w(lw(12,ibl)), jdim,   idim-1, 3,      2);
        FortranArray4DRef<double> xtbi (&w(lw(13,ibl)), jdim,   kdim,   3,      2);
        FortranArray4DRef<double> atbj (&w(lw(14,ibl)), kdim,   idim-1, 3,      2);
        FortranArray4DRef<double> atbk (&w(lw(15,ibl)), jdim,   idim-1, 3,      2);
        FortranArray4DRef<double> atbi (&w(lw(16,ibl)), jdim,   kdim,   3,      2);

        // Build array views from wk using lw2 offsets
        FortranArray4DRef<double> tursav(&wk(lw2(1,ibl)), jdim,   kdim,   idim,   nummem);
        FortranArray4DRef<double> tj0   (&wk(lw2(2,ibl)), kdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> tk0   (&wk(lw2(3,ibl)), jdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> ti0   (&wk(lw2(4,ibl)), jdim,   kdim,   nummem, 4);
        FortranArray3DRef<double> vist3d(&wk(lw2(5,ibl)), jdim,   kdim,   idim);
        FortranArray4DRef<double> vj0   (&wk(lw2(6,ibl)), kdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vk0   (&wk(lw2(7,ibl)), jdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vi0   (&wk(lw2(8,ibl)), jdim,   kdim,   1,      4);

        // Loop over j-face boundary segments
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl,iseg,1,1) == 1011) {
                int ista  = jbcinfo(ibl,iseg,1,2);
                int iend  = jbcinfo(ibl,iseg,1,3);
                int jsta  = jbcinfo(ibl,iseg,1,4);
                int jend  = jbcinfo(ibl,iseg,1,5);
                int ksta  = jbcinfo(ibl,iseg,1,6);
                int kend  = jbcinfo(ibl,iseg,1,7);
                int nface = jbcinfo(ibl,iseg,1,8);
                int isym  = jbcinfo(ibl,iseg,1,9);
                int jsym  = jbcinfo(ibl,iseg,1,10);
                int ksym  = jbcinfo(ibl,iseg,1,11);
                bc1011(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       isym, jsym, ksym, iuns, nou, bou, nbuf, ibufdim, nummem);
            }
        }

        // Loop over k-face boundary segments
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) == 1011) {
                int ista  = kbcinfo(ibl,iseg,1,2);
                int iend  = kbcinfo(ibl,iseg,1,3);
                int jsta  = kbcinfo(ibl,iseg,1,4);
                int jend  = kbcinfo(ibl,iseg,1,5);
                int ksta  = kbcinfo(ibl,iseg,1,6);
                int kend  = kbcinfo(ibl,iseg,1,7);
                int nface = kbcinfo(ibl,iseg,1,8);
                int isym  = kbcinfo(ibl,iseg,1,9);
                int jsym  = kbcinfo(ibl,iseg,1,10);
                int ksym  = kbcinfo(ibl,iseg,1,11);
                bc1011(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       isym, jsym, ksym, iuns, nou, bou, nbuf, ibufdim, nummem);
            }
        }

        // Loop over i-face boundary segments
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) == 1011) {
                int ista  = ibcinfo(ibl,iseg,1,2);
                int iend  = ibcinfo(ibl,iseg,1,3);
                int jsta  = ibcinfo(ibl,iseg,1,4);
                int jend  = ibcinfo(ibl,iseg,1,5);
                int ksta  = ibcinfo(ibl,iseg,1,6);
                int kend  = ibcinfo(ibl,iseg,1,7);
                int nface = ibcinfo(ibl,iseg,1,8);
                int isym  = ibcinfo(ibl,iseg,1,9);
                int jsym  = ibcinfo(ibl,iseg,1,10);
                int ksym  = ibcinfo(ibl,iseg,1,11);
                bc1011(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       isym, jsym, ksym, iuns, nou, bou, nbuf, ibufdim, nummem);
            }
        }

    } // end ibl loop

    return;
} // end bc

} // namespace bc1011_ns
