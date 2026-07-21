// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "getibk0.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstdint>
#include <algorithm>

namespace getibk0_ns {



void getibk0(int& jdim, int& kdim, int& idim, int& nbl, int& itotb, int& itoti, int& maxbl, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& ibufdim, int& nbuf, int& ierrflg, int& myid)
{
    // Local variables
    int idim1, jdim1, kdim1;
    int lsta, lend;
    int itotmx;
    int jchk, kchk, lchk;
    int ibpnts, iipnts, idum;
    int intpts[4];  // local dimension intpts(4)
    // Scalars read in implied-DO loops (discarded)
    int jji, kki, iii;
    float dyint, dzint, dxint;
    int jjb, kkb, iib, ibc;
    float blank_val;

    // COMMON /params/
    int& liitot = cmn_params.liitot;

    idim1 = idim - 1;
    jdim1 = jdim - 1;
    kdim1 = kdim - 1;

    // write(66,101) nbl
    // 101 format(1x,37hreading overlap information for block,1x,i3)
    FILE* f66 = fortran_get_unit(66);
    fprintf(f66, " reading overlap information for block  %3d\n", nbl);

    // read(21) jchk,kchk,lchk
    FILE* f21 = fortran_get_unit(21);
    {
        int32_t rec_start, rec_end;
        fread(&rec_start, sizeof(int32_t), 1, f21);
        fread(&jchk, sizeof(int32_t), 1, f21);
        fread(&kchk, sizeof(int32_t), 1, f21);
        fread(&lchk, sizeof(int32_t), 1, f21);
        fread(&rec_end, sizeof(int32_t), 1, f21);
    }

    if (jchk != jdim1 || kchk != kdim1 || lchk != idim1) {
        fprintf(f66, "  mismatch in indices.....stopping in getibk\n");
        fprintf(f66, "  jdim,kdim,idim= %d %d %d\n", jdim, kdim, idim);
        fprintf(f66, "  jmax,kmax,lmax= %d %d %d\n", jchk, kchk, lchk);
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    // read(21) ibpnts,intpts,iipnts,idum,idum
    {
        int32_t rec_start, rec_end;
        fread(&rec_start, sizeof(int32_t), 1, f21);
        fread(&ibpnts, sizeof(int32_t), 1, f21);
        fread(&intpts[0], sizeof(int32_t), 1, f21);
        fread(&intpts[1], sizeof(int32_t), 1, f21);
        fread(&intpts[2], sizeof(int32_t), 1, f21);
        fread(&intpts[3], sizeof(int32_t), 1, f21);
        fread(&iipnts, sizeof(int32_t), 1, f21);
        fread(&idum, sizeof(int32_t), 1, f21);
        fread(&idum, sizeof(int32_t), 1, f21);
        fread(&rec_end, sizeof(int32_t), 1, f21);
    }

    // write(66,*) '  ibpnts,intpts,iipnts = ',ibpnts,intpts,iipnts
    fprintf(f66, "  ibpnts,intpts,iipnts =  %d %d %d %d %d %d\n",
            ibpnts, intpts[0], intpts[1], intpts[2], intpts[3], iipnts);

    ibpntsg(nbl, 1) = intpts[0];
    ibpntsg(nbl, 2) = intpts[1];
    ibpntsg(nbl, 3) = intpts[2];
    ibpntsg(nbl, 4) = intpts[3];
    iipntsg(nbl)    = iipnts;

    lsta = lig(nbl);
    lend = lsta + iipnts - 1;

    // read(21)(jji,kki,iii,dyint,dzint,dxint, l=lsta,lend)
    // Implied-DO: reads iipnts records worth of data; variables are discarded
    {
        int32_t rec_start, rec_end;

        fread(&rec_start, sizeof(int32_t), 1, f21);
        // sizing pass discards these values; skip the whole record by its
        // byte-count marker so we're agnostic to real*4 vs real*8. The coeffs
        // are real*8 (maggie built -fdefault-real-8); reading them as float
        // consumed half the record and desynced the next block's header.
        (void)lsta; (void)lend;
        fseek(f21, (long)rec_start, SEEK_CUR);
        fread(&rec_end, sizeof(int32_t), 1, f21);
    }

    lsta = lbg(nbl);
    lend = lsta + ibpnts - 1;

    // read(21)(jjb,kkb,iib,ibc,l=lsta,lend)
    {
        int32_t rec_start, rec_end;
        fread(&rec_start, sizeof(int32_t), 1, f21);
        for (int l = lsta; l <= lend; l++) {
            fread(&jjb, sizeof(int32_t), 1, f21);
            fread(&kkb, sizeof(int32_t), 1, f21);
            fread(&iib, sizeof(int32_t), 1, f21);
            fread(&ibc, sizeof(int32_t), 1, f21);
        }
        fread(&rec_end, sizeof(int32_t), 1, f21);
    }

    // read(21)(((blank,j=1,jdim1),k=1,kdim1),i=1,idim1)
    // Reads jdim1*kdim1*idim1 real values; blank is discarded
    {
        int32_t rec_start, rec_end;
        fread(&rec_start, sizeof(int32_t), 1, f21);
        // blank array is discarded in the sizing pass; skip by byte count
        // (real*8, precision-agnostic).
        (void)blank_val;
        fseek(f21, (long)rec_start, SEEK_CUR);
        fread(&rec_end, sizeof(int32_t), 1, f21);
    }

    if (nbl < maxbl) {
        lig(nbl + 1) = lig(nbl) + iipntsg(nbl);
        lbg(nbl + 1) = lbg(nbl) + ibpnts;
    }

    itoti = itoti + iipntsg(nbl);
    itotb = itotb + ibpnts;

    itotmx = std::max(itoti, itotb);
    liitot = std::max(itotmx, liitot);
}

} // namespace getibk0_ns
