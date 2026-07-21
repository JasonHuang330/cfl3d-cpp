// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rsurf.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <algorithm>
#include <cstdio>
#include <cstring>

namespace rsurf_ns {

void rsurf(int& maxbl, int& maxsegdg, int& idim, int& jdim, int& kdim,
           FortranArray4DRef<double> delti, FortranArray4DRef<double> deltj,
           FortranArray4DRef<double> deltk, int& nbl,
           FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
           FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
           FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
           int& iseg, FortranArray1DRef<int> nou,
           FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
           int& iunit)
{
    // COMMON blocks
    int& ialph  = cmn_igrdtyp.ialph;
    int& myid   = cmn_mydist2.myid;

    // local variables
    int isurf_img = 0;
    int is, ie, js, je, ks, ke, mm;
    double temp;

    is = icsi(nbl, iseg);
    ie = icsf(nbl, iseg);
    js = jcsi(nbl, iseg);
    je = jcsf(nbl, iseg);
    ks = kcsi(nbl, iseg);
    ke = kcsf(nbl, iseg);

    if (is == ie) {
        mm = 1;
        if (is == idim) mm = 2;
        // read(iunit,*)((delti(j,k,1,mm),j=js,je),k=ks,ke),
        //              ((delti(j,k,2,mm),j=js,je),k=ks,ke),
        //              ((delti(j,k,3,mm),j=js,je),k=ks,ke)
        // Single list-directed READ with nested implied-DO: outer k, inner j
        {
            int nj = je - js + 1;
            int nk = ke - ks + 1;
            int ntot = nj * nk;
            FortranArray1D<double> tmp(ntot);
            // component 1
            fortran_read_list_array(iunit, tmp, 1, ntot);
            {
                int idx = 1;
                for (int k = ks; k <= ke; k++)
                    for (int j = js; j <= je; j++)
                        delti(j, k, 1, mm) = tmp(idx++);
            }
            // component 2
            fortran_read_list_array(iunit, tmp, 1, ntot);
            {
                int idx = 1;
                for (int k = ks; k <= ke; k++)
                    for (int j = js; j <= je; j++)
                        delti(j, k, 2, mm) = tmp(idx++);
            }
            // component 3
            fortran_read_list_array(iunit, tmp, 1, ntot);
            {
                int idx = 1;
                for (int k = ks; k <= ke; k++)
                    for (int j = js; j <= je; j++)
                        delti(j, k, 3, mm) = tmp(idx++);
            }
        }
        if (ialph != 0) {
            for (int k = ks; k <= ke; k++) {
                for (int j = js; j <= je; j++) {
                    temp               = delti(j, k, 3, mm);
                    delti(j, k, 3, mm) = delti(j, k, 2, mm);
                    delti(j, k, 2, mm) = -temp;
                }
            }
        }
    } else if (js == je) {
        mm = 1;
        if (js == jdim) mm = 2;
        // read(iunit,*)((deltj(k,i,1,mm),i=is,ie),k=ks,ke),
        //              ((deltj(k,i,2,mm),i=is,ie),k=ks,ke),
        //              ((deltj(k,i,3,mm),i=is,ie),k=ks,ke)
        // Single list-directed READ with nested implied-DO: outer k, inner i
        {
            int ni = ie - is + 1;
            int nk = ke - ks + 1;
            int ntot = ni * nk;
            FortranArray1D<double> tmp(ntot);
            // component 1
            fortran_read_list_array(iunit, tmp, 1, ntot);
            {
                int idx = 1;
                for (int k = ks; k <= ke; k++)
                    for (int i = is; i <= ie; i++)
                        deltj(k, i, 1, mm) = tmp(idx++);
            }
            // component 2
            fortran_read_list_array(iunit, tmp, 1, ntot);
            {
                int idx = 1;
                for (int k = ks; k <= ke; k++)
                    for (int i = is; i <= ie; i++)
                        deltj(k, i, 2, mm) = tmp(idx++);
            }
            // component 3
            fortran_read_list_array(iunit, tmp, 1, ntot);
            {
                int idx = 1;
                for (int k = ks; k <= ke; k++)
                    for (int i = is; i <= ie; i++)
                        deltj(k, i, 3, mm) = tmp(idx++);
            }
        }
        if (ialph != 0) {
            for (int i = is; i <= ie; i++) {
                for (int k = ks; k <= ke; k++) {
                    temp               = deltj(k, i, 3, mm);
                    deltj(k, i, 3, mm) = deltj(k, i, 2, mm);
                    deltj(k, i, 2, mm) = -temp;
                }
            }
        }
    } else if (ks == ke) {
        mm = 1;
        if (ks == kdim) mm = 2;
        // read(iunit,*)((deltk(j,i,1,mm),i=is,ie),j=js,je),
        //              ((deltk(j,i,2,mm),i=is,ie),j=js,je),
        //              ((deltk(j,i,3,mm),i=is,ie),j=js,je)
        // Single list-directed READ with nested implied-DO: outer j, inner i
        {
            int ni = ie - is + 1;
            int nj = je - js + 1;
            int ntot = ni * nj;
            FortranArray1D<double> tmp(ntot);
            // component 1
            fortran_read_list_array(iunit, tmp, 1, ntot);
            {
                int idx = 1;
                for (int j = js; j <= je; j++)
                    for (int i = is; i <= ie; i++)
                        deltk(j, i, 1, mm) = tmp(idx++);
            }
            // component 2
            fortran_read_list_array(iunit, tmp, 1, ntot);
            {
                int idx = 1;
                for (int j = js; j <= je; j++)
                    for (int i = is; i <= ie; i++)
                        deltk(j, i, 2, mm) = tmp(idx++);
            }
            // component 3
            fortran_read_list_array(iunit, tmp, 1, ntot);
            {
                int idx = 1;
                for (int j = js; j <= je; j++)
                    for (int i = is; i <= ie; i++)
                        deltk(j, i, 3, mm) = tmp(idx++);
            }
        }
        if (ialph != 0) {
            for (int i = is; i <= ie; i++) {
                for (int j = js; j <= je; j++) {
                    temp               = deltk(j, i, 3, mm);
                    deltk(j, i, 3, mm) = deltk(j, i, 2, mm);
                    deltk(j, i, 2, mm) = -temp;
                }
            }
        }
    } else {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            " error in rsurf...one surface dimension must be 1");
        int m1 = -1;
        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
    }

    return;
}

} // namespace rsurf_ns
