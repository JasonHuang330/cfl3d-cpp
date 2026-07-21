// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc9999.h"
#include "mms.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <algorithm>
#include <cstdio>
#include <cstring>

namespace bc9999_ns {

void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork,
        double& cl, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
        int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg,
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

        FortranArray4DRef<double> q(&w(lw(1,ibl)), jdim, kdim, idim, 5);
        FortranArray4DRef<double> qj0(&w(lw(2,ibl)), kdim, idim-1, 5, 4);
        FortranArray4DRef<double> qk0(&w(lw(3,ibl)), jdim, idim-1, 5, 4);
        FortranArray4DRef<double> qi0(&w(lw(4,ibl)), jdim, kdim, 5, 4);
        FortranArray4DRef<double> sj(&w(lw(5,ibl)), jdim, kdim, idim-1, 5);
        FortranArray4DRef<double> sk(&w(lw(6,ibl)), jdim, kdim, idim-1, 5);
        FortranArray4DRef<double> si(&w(lw(7,ibl)), jdim, kdim, idim, 5);
        FortranArray3DRef<double> bcj(&w(lw(8,ibl)), kdim, idim-1, 2);
        FortranArray3DRef<double> bck(&w(lw(9,ibl)), jdim, idim-1, 2);
        FortranArray3DRef<double> bci(&w(lw(10,ibl)), jdim, kdim, 2);
        FortranArray4DRef<double> xtbj(&w(lw(11,ibl)), kdim, idim-1, 3, 2);
        FortranArray4DRef<double> xtbk(&w(lw(12,ibl)), jdim, idim-1, 3, 2);
        FortranArray4DRef<double> xtbi(&w(lw(13,ibl)), jdim, kdim, 3, 2);
        FortranArray4DRef<double> atbj(&w(lw(14,ibl)), kdim, idim-1, 3, 2);
        FortranArray4DRef<double> atbk(&w(lw(15,ibl)), jdim, idim-1, 3, 2);
        FortranArray4DRef<double> atbi(&w(lw(16,ibl)), jdim, kdim, 3, 2);

        FortranArray4DRef<double> tursav(&w(lw2(1,ibl)), jdim, kdim, idim, nummem);
        FortranArray4DRef<double> tj0(&w(lw2(2,ibl)), kdim, idim-1, nummem, 4);
        FortranArray4DRef<double> tk0(&w(lw2(3,ibl)), jdim, idim-1, nummem, 4);
        FortranArray4DRef<double> ti0(&w(lw2(4,ibl)), jdim, kdim, nummem, 4);
        FortranArray3DRef<double> vist3d(&w(lw2(5,ibl)), jdim, kdim, idim);
        FortranArray4DRef<double> vj0(&w(lw2(6,ibl)), kdim, idim-1, 1, 4);
        FortranArray4DRef<double> vk0(&w(lw2(7,ibl)), jdim, idim-1, 1, 4);
        FortranArray4DRef<double> vi0(&w(lw2(8,ibl)), jdim, kdim, 1, 4);

        FortranArray3DRef<double> x(&w(lw(17,ibl)), jdim, kdim, idim);
        FortranArray3DRef<double> y(&w(lw(18,ibl)), jdim, kdim, idim);
        FortranArray3DRef<double> z(&w(lw(19,ibl)), jdim, kdim, idim);
        FortranArray3DRef<double> smin(&w(lw(20,ibl)), jdim-1, kdim-1, idim-1);

        // Loop over j-face BC segments
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl,iseg,1,1) == 9999) {
                int ista  = jbcinfo(ibl,iseg,1,2);
                int iend  = jbcinfo(ibl,iseg,1,3);
                int jsta  = jbcinfo(ibl,iseg,1,4);
                int jend  = jbcinfo(ibl,iseg,1,5);
                int ksta  = jbcinfo(ibl,iseg,1,6);
                int kend  = jbcinfo(ibl,iseg,1,7);
                int nface = jbcinfo(ibl,iseg,1,8);
                bc9999(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem,
                       x, y, z, smin);
            }
        }

        // Loop over k-face BC segments
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) == 9999) {
                int ista  = kbcinfo(ibl,iseg,1,2);
                int iend  = kbcinfo(ibl,iseg,1,3);
                int jsta  = kbcinfo(ibl,iseg,1,4);
                int jend  = kbcinfo(ibl,iseg,1,5);
                int ksta  = kbcinfo(ibl,iseg,1,6);
                int kend  = kbcinfo(ibl,iseg,1,7);
                int nface = kbcinfo(ibl,iseg,1,8);
                bc9999(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem,
                       x, y, z, smin);
            }
        }

        // Loop over i-face BC segments
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) == 9999) {
                int ista  = ibcinfo(ibl,iseg,1,2);
                int iend  = ibcinfo(ibl,iseg,1,3);
                int jsta  = ibcinfo(ibl,iseg,1,4);
                int jend  = ibcinfo(ibl,iseg,1,5);
                int ksta  = ibcinfo(ibl,iseg,1,6);
                int kend  = ibcinfo(ibl,iseg,1,7);
                int nface = ibcinfo(ibl,iseg,1,8);
                bc9999(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem,
                       x, y, z, smin);
            }
        }
    }
}

void bc9999(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
            FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si, FortranArray3DRef<double> bcj,
            FortranArray3DRef<double> bck, FortranArray3DRef<double> bci,
            FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk,
            FortranArray4DRef<double> xtbi, FortranArray4DRef<double> atbj,
            FortranArray4DRef<double> atbk, FortranArray4DRef<double> atbi,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            int& nface, FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0, FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0, int& iuns,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& nummem,
            FortranArray3DRef<double> x, FortranArray3DRef<double> y,
            FortranArray3DRef<double> z, FortranArray3DRef<double> smin)
{
    // COMMON /sklton/ isklton
    int& isklton = cmn_sklton.isklton;
    // COMMON /maxiv/ ivmx
    int& ivmx = cmn_maxiv.ivmx;
    // COMMON /mms/ iexact_trunc, iexact_disc, iexact_ring
    int& iexact_trunc = cmn_mms.iexact_trunc;
    int& iexact_disc  = cmn_mms.iexact_disc;
    int& iexact_ring  = cmn_mms.iexact_ring;

    // myid from COMMON /mpicomm/ → cmn_mydist2
    int& myid = cmn_mydist2.myid;

    if (iexact_trunc != 0 || iexact_disc != 0) {
        mms_ns::exact_flow_bc(jdim, kdim, idim, x, y, z, qi0, qj0, qk0,
                              nface, bcj, bck, bci, smin,
                              tj0, tk0, ti0, vj0, vk0, vi0,
                              nummem, iexact_trunc, iexact_disc);
    } else {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120,
            " Error... bc9999 only usable with"
            " iexact_trunc .ne. 0 or iexact_disc .ne. 0");
        int m1 = -1;
        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
    }
}

} // namespace bc9999_ns
