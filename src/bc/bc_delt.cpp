// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc_delt.h"
#include "bc.h"

namespace bc_delt_ns {

void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, double& cl, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, FortranArray3DRef<int> lwdat, int& myid, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, int& nummem)
{
    bc_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, cl, nou, bou, nbuf, ibufdim, maxbl, maxgr, maxseg, itrans, irotat, idefrm, igridg, nblg, nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim, ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek, lwdat, myid, idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
}

void bc_delt(int& nbl, FortranArray3DRef<double> dx, FortranArray3DRef<double> dy, FortranArray3DRef<double> dz, FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, FortranArray4DRef<double> delti, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, int& jdim, int& kdim, int& idim, int& maxbl, int& maxsegdg, FortranArray1DRef<int> nsegdfrm)
{
    int i, j, k;
    int ii, jj, kk;

    for (int iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {

        if (icsi(nbl, iseg) == icsf(nbl, iseg)) {
            i = icsi(nbl, iseg);
            if (icsi(nbl, iseg) == 1) {
                ii = 1;
            } else {
                ii = 2;
            }
            for (j = jcsi(nbl, iseg); j <= jcsf(nbl, iseg); j++) {
                for (k = kcsi(nbl, iseg); k <= kcsf(nbl, iseg); k++) {
                    dx(j, k, i) = delti(j, k, 1, ii);
                    dy(j, k, i) = delti(j, k, 2, ii);
                    dz(j, k, i) = delti(j, k, 3, ii);
                }
            }
        }

        if (jcsi(nbl, iseg) == jcsf(nbl, iseg)) {
            j = jcsi(nbl, iseg);
            if (jcsi(nbl, iseg) == 1) {
                jj = 1;
            } else {
                jj = 2;
            }
            for (k = kcsi(nbl, iseg); k <= kcsf(nbl, iseg); k++) {
                for (i = icsi(nbl, iseg); i <= icsf(nbl, iseg); i++) {
                    dx(j, k, i) = deltj(k, i, 1, jj);
                    dy(j, k, i) = deltj(k, i, 2, jj);
                    dz(j, k, i) = deltj(k, i, 3, jj);
                }
            }
        }

        if (kcsi(nbl, iseg) == kcsf(nbl, iseg)) {
            k = kcsi(nbl, iseg);
            if (kcsi(nbl, iseg) == 1) {
                kk = 1;
            } else {
                kk = 2;
            }
            for (j = jcsi(nbl, iseg); j <= jcsf(nbl, iseg); j++) {
                for (i = icsi(nbl, iseg); i <= icsf(nbl, iseg); i++) {
                    dx(j, k, i) = deltk(j, i, 1, kk);
                    dy(j, k, i) = deltk(j, i, 2, kk);
                    dz(j, k, i) = deltk(j, i, 3, kk);
                }
            }
        }

    }
}

} // namespace bc_delt_ns
