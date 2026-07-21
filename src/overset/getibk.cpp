// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "getibk.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstdint>

namespace getibk_ns {

void getibk(FortranArray3DRef<double> blank, int& jdim, int& kdim, int& idim, int& nbl, FortranArray1DRef<int> intpts, int& nblpts, int& ibpnts, int& iipnts, int& iitot, int& maxbl, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray1DRef<double> dxintg, FortranArray1DRef<double> dyintg, FortranArray1DRef<double> dzintg, FortranArray1DRef<int> iiig, FortranArray1DRef<int> jjig, FortranArray1DRef<int> kkig, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, int& myid, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou)
{
    int idim1 = idim - 1;
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;

    // initialize all iblank values to 1
    for (int i = 1; i <= idim; i++)
        for (int j = 1; j <= jdim; j++)
            for (int k = 1; k <= kdim; k++)
                blank(j, k, i) = 1.0;

    // read(21) jchk,kchk,lchk
    FILE* fp21 = fortran_get_unit(21);
    int jchk, kchk, lchk;
    {
        int32_t rec_start, rec_end;
        fread(&rec_start, sizeof(int32_t), 1, fp21);
        fread(&jchk, sizeof(int), 1, fp21);
        fread(&kchk, sizeof(int), 1, fp21);
        fread(&lchk, sizeof(int), 1, fp21);
        fread(&rec_end, sizeof(int32_t), 1, fp21);
    }

    if (jchk != jdim1 || kchk != kdim1 || lchk != idim1) {
        fortran_write_unit(11, "  mismatch in indices.....stopping in getibk\n");
        fortran_write_unit(11, "  jdim,kdim,idim= %d %d %d\n", jdim, kdim, idim);
        fortran_write_unit(11, "  jmax,kmax,lmax= %d %d %d\n", jchk, kchk, lchk);
        int ierrflg = -1;
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        return;
    }

    // read(21) ibpnts,intpts,iipnts,idum,idum
    {
        int idum1, idum2;
        int32_t rec_start, rec_end;
        fread(&rec_start, sizeof(int32_t), 1, fp21);
        fread(&ibpnts, sizeof(int), 1, fp21);
        fread(&intpts(1), sizeof(int), 1, fp21);
        fread(&intpts(2), sizeof(int), 1, fp21);
        fread(&intpts(3), sizeof(int), 1, fp21);
        fread(&intpts(4), sizeof(int), 1, fp21);
        fread(&iipnts, sizeof(int), 1, fp21);
        fread(&idum1, sizeof(int), 1, fp21);
        fread(&idum2, sizeof(int), 1, fp21);
        fread(&rec_end, sizeof(int32_t), 1, fp21);
    }
    ibpntsg(nbl, 1) = intpts(1);
    ibpntsg(nbl, 2) = intpts(2);
    ibpntsg(nbl, 3) = intpts(3);
    ibpntsg(nbl, 4) = intpts(4);
    iipntsg(nbl)    = iipnts;

    int lsta = lig(nbl);
    int lend = lsta + iipnts - 1;
    if (lend > iitot) {
        fortran_write_unit(11, " stopping in getibk - insufficient memory allocation\n");
        fortran_write_unit(11, " lend,iitot = %20d%20d\n", lend, iitot);
        int minus1 = -1;
        termn8_ns::termn8(myid, minus1, ibufdim, nbuf, bou, nou);
        return;
    }

    // read(21)(jjig(l),kkig(l),iiig(l),dyintg(l),dzintg(l),dxintg(l),l=lsta,lend)
    {
        int32_t rec_start, rec_end;
        fread(&rec_start, sizeof(int32_t), 1, fp21);
        for (int l = lsta; l <= lend; l++) {
            fread(&jjig(l),   sizeof(int),    1, fp21);
            fread(&kkig(l),   sizeof(int),    1, fp21);
            fread(&iiig(l),   sizeof(int),    1, fp21);
            fread(&dyintg(l), sizeof(double), 1, fp21);
            fread(&dzintg(l), sizeof(double), 1, fp21);
            fread(&dxintg(l), sizeof(double), 1, fp21);
        }
        fread(&rec_end, sizeof(int32_t), 1, fp21);
    }

    lsta = lbg(nbl);
    lend = lsta + ibpnts - 1;
    if (lend > iitot) {
        fortran_write_unit(11, " stopping in getibk - insufficient memory allocation\n");
        fortran_write_unit(11, " lend,iitot = %20d%20d\n", lend, iitot);
        int minus1 = -1;
        termn8_ns::termn8(myid, minus1, ibufdim, nbuf, bou, nou);
        return;
    }

    // read(21)(jjbg(l),kkbg(l),iibg(l),ibcg(l),l=lsta,lend)
    {
        int32_t rec_start, rec_end;
        fread(&rec_start, sizeof(int32_t), 1, fp21);
        for (int l = lsta; l <= lend; l++) {
            fread(&jjbg(l), sizeof(int), 1, fp21);
            fread(&kkbg(l), sizeof(int), 1, fp21);
            fread(&iibg(l), sizeof(int), 1, fp21);
            fread(&ibcg(l), sizeof(int), 1, fp21);
        }
        fread(&rec_end, sizeof(int32_t), 1, fp21);
    }

    // read(21)(((blank(j,k,i),j=1,jdim1),k=1,kdim1),i=1,idim1)
    {
        int32_t rec_start, rec_end;
        fread(&rec_start, sizeof(int32_t), 1, fp21);
        for (int i = 1; i <= idim1; i++)
            for (int k = 1; k <= kdim1; k++)
                for (int j = 1; j <= jdim1; j++) {
                    double val;
                    fread(&val, sizeof(double), 1, fp21);
                    blank(j, k, i) = val;
                }
        fread(&rec_end, sizeof(int32_t), 1, fp21);
    }

    // check iblank values (commented-out writes preserved as no-ops)
    for (int i = 1; i <= idim1; i++) {
        int isum = 0;
        for (int k = 1; k <= kdim1; k++)
            for (int j = 1; j <= jdim1; j++)
                if ((int)blank(j, k, i) != 1) isum = isum + 1;
        if (isum == 0) {
            // write(15,*) ' for i=',i,' all iblank values are 1'
        } else {
            // write(15,*) ' iblank array...i=',i
            for (int k = 1; k <= kdim1; k++) {
                // write(15,8833) (int(blank(j,k,i)),j=1,jdim1)
            }
        }
    }
    // 8833 format(1x,71i1)

    nblpts = 0;
    for (int i = 1; i <= idim1; i++)
        for (int k = 1; k <= kdim1; k++)
            for (int j = 1; j <= jdim1; j++)
                if (blank(j, k, i) != 1.0)
                    nblpts = nblpts + 1;

    if (nbl < maxbl) {
        lig(nbl + 1) = lig(nbl) + iipntsg(nbl);
        lbg(nbl + 1) = lbg(nbl) + ibpnts;
    }
}

} // namespace getibk_ns
