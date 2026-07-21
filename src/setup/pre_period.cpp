// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "pre_period.h"
#include "lead.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cstring>
#include <algorithm>

namespace pre_period_ns {

void pre_period(int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, int& icount, int& maxbl, int& maxseg, int& lbcprd, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbcidim, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, FortranArray1DRef<int> igridg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<int> isav_prd, FortranArray1DRef<int> is_prd, FortranArray1DRef<int> ie_prd, int& nbcprd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk, int& myid, FortranArray1DRef<int> nblg, int& maxgr, int& ierrflg)
{
    // Local variables from COMMON /ginfo/
    int& jdim  = cmn_ginfo.jdim;
    int& kdim  = cmn_ginfo.kdim;
    int& idim  = cmn_ginfo.idim;

    // Local scalars
    int nseg, ista, iend, jsta, jend, ksta, kend, mdim, ndim, nblp, nface;

    lead_ns::lead(nbl, lw, lw2, maxbl);

    // check for 2005 type bcs

    // Loop 802: i=1 face (ista=1, iend=1)
    for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
        ista = 1;
        iend = 1;
        jsta = ibcinfo(nbl, nseg, 2, 1);
        jend = ibcinfo(nbl, nseg, 3, 1);
        ksta = ibcinfo(nbl, nseg, 4, 1);
        kend = ibcinfo(nbl, nseg, 5, 1);
        mdim = jend - jsta;
        ndim = kend - ksta;
        if (ibcinfo(nbl, nseg, 1, 1) == 2005) {
            nblp = nblg((int)bcvali(nbl, nseg, 1, 1))
                 + (nbl - nblg(igridg(nbl)));
            nface  = 1;
            icount = icount + 1;
            isav_prd(icount, 1)  = nbl;
            isav_prd(icount, 2)  = nface;
            isav_prd(icount, 3)  = ista;
            isav_prd(icount, 4)  = iend;
            isav_prd(icount, 5)  = jsta;
            isav_prd(icount, 6)  = jend;
            isav_prd(icount, 7)  = ksta;
            isav_prd(icount, 8)  = kend;
            isav_prd(icount, 9)  = mdim;
            isav_prd(icount, 10) = ndim;
            isav_prd(icount, 11) = nseg;
            isav_prd(icount, 12) = nblp;
        }
    }

    // Loop 803: i=idim face (ista=idim, iend=idim)
    for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
        ista = idim;
        iend = idim;
        jsta = ibcinfo(nbl, nseg, 2, 2);
        jend = ibcinfo(nbl, nseg, 3, 2);
        ksta = ibcinfo(nbl, nseg, 4, 2);
        kend = ibcinfo(nbl, nseg, 5, 2);
        mdim = jend - jsta;
        ndim = kend - ksta;
        if (ibcinfo(nbl, nseg, 1, 2) == 2005) {
            nblp = nblg((int)bcvali(nbl, nseg, 1, 2))
                 + (nbl - nblg(igridg(nbl)));
            nface  = 2;
            icount = icount + 1;
            isav_prd(icount, 1)  = nbl;
            isav_prd(icount, 2)  = nface;
            isav_prd(icount, 3)  = ista;
            isav_prd(icount, 4)  = iend;
            isav_prd(icount, 5)  = jsta;
            isav_prd(icount, 6)  = jend;
            isav_prd(icount, 7)  = ksta;
            isav_prd(icount, 8)  = kend;
            isav_prd(icount, 9)  = mdim;
            isav_prd(icount, 10) = ndim;
            isav_prd(icount, 11) = nseg;
            isav_prd(icount, 12) = nblp;
        }
    }

    // Loop 804: j=1 face (jsta=1, jend=1)
    for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
        ista = jbcinfo(nbl, nseg, 2, 1);
        iend = jbcinfo(nbl, nseg, 3, 1);
        jsta = 1;
        jend = 1;
        ksta = jbcinfo(nbl, nseg, 4, 1);
        kend = jbcinfo(nbl, nseg, 5, 1);
        mdim = kend - ksta;
        ndim = iend - ista;
        if (jbcinfo(nbl, nseg, 1, 1) == 2005) {
            nblp = nblg((int)bcvalj(nbl, nseg, 1, 1))
                 + (nbl - nblg(igridg(nbl)));
            nface  = 3;
            icount = icount + 1;
            isav_prd(icount, 1)  = nbl;
            isav_prd(icount, 2)  = nface;
            isav_prd(icount, 3)  = ista;
            isav_prd(icount, 4)  = iend;
            isav_prd(icount, 5)  = jsta;
            isav_prd(icount, 6)  = jend;
            isav_prd(icount, 7)  = ksta;
            isav_prd(icount, 8)  = kend;
            isav_prd(icount, 9)  = mdim;
            isav_prd(icount, 10) = ndim;
            isav_prd(icount, 11) = nseg;
            isav_prd(icount, 12) = nblp;
        }
    }

    // Loop 805: j=jdim face (jsta=jdim, jend=jdim)
    for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
        ista = jbcinfo(nbl, nseg, 2, 2);
        iend = jbcinfo(nbl, nseg, 3, 2);
        jsta = jdim;
        jend = jdim;
        ksta = jbcinfo(nbl, nseg, 4, 2);
        kend = jbcinfo(nbl, nseg, 5, 2);
        mdim = kend - ksta;
        ndim = iend - ista;
        if (jbcinfo(nbl, nseg, 1, 2) == 2005) {
            nblp = nblg((int)bcvalj(nbl, nseg, 1, 2))
                 + (nbl - nblg(igridg(nbl)));
            nface  = 4;
            icount = icount + 1;
            isav_prd(icount, 1)  = nbl;
            isav_prd(icount, 2)  = nface;
            isav_prd(icount, 3)  = ista;
            isav_prd(icount, 4)  = iend;
            isav_prd(icount, 5)  = jsta;
            isav_prd(icount, 6)  = jend;
            isav_prd(icount, 7)  = ksta;
            isav_prd(icount, 8)  = kend;
            isav_prd(icount, 9)  = mdim;
            isav_prd(icount, 10) = ndim;
            isav_prd(icount, 11) = nseg;
            isav_prd(icount, 12) = nblp;
        }
    }

    // Loop 806: k=1 face (ksta=1, kend=1)
    for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
        ista = kbcinfo(nbl, nseg, 2, 1);
        iend = kbcinfo(nbl, nseg, 3, 1);
        jsta = kbcinfo(nbl, nseg, 4, 1);
        jend = kbcinfo(nbl, nseg, 5, 1);
        ksta = 1;
        kend = 1;
        mdim = jend - jsta;
        ndim = iend - ista;
        if (kbcinfo(nbl, nseg, 1, 1) == 2005) {
            nblp = nblg((int)bcvalk(nbl, nseg, 1, 1))
                 + (nbl - nblg(igridg(nbl)));
            nface  = 5;
            icount = icount + 1;
            isav_prd(icount, 1)  = nbl;
            isav_prd(icount, 2)  = nface;
            isav_prd(icount, 3)  = ista;
            isav_prd(icount, 4)  = iend;
            isav_prd(icount, 5)  = jsta;
            isav_prd(icount, 6)  = jend;
            isav_prd(icount, 7)  = ksta;
            isav_prd(icount, 8)  = kend;
            isav_prd(icount, 9)  = mdim;
            isav_prd(icount, 10) = ndim;
            isav_prd(icount, 11) = nseg;
            isav_prd(icount, 12) = nblp;
        }
    }

    // Loop 807: k=kdim face (ksta=kdim, kend=kdim)
    for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
        ista = kbcinfo(nbl, nseg, 2, 2);
        iend = kbcinfo(nbl, nseg, 3, 2);
        jsta = kbcinfo(nbl, nseg, 4, 2);
        jend = kbcinfo(nbl, nseg, 5, 2);
        ksta = kdim;
        kend = kdim;
        mdim = jend - jsta;
        ndim = iend - ista;
        if (kbcinfo(nbl, nseg, 1, 2) == 2005) {
            nblp = nblg((int)bcvalk(nbl, nseg, 1, 2))
                 + (nbl - nblg(igridg(nbl)));
            nface  = 6;
            icount = icount + 1;
            isav_prd(icount, 1)  = nbl;
            isav_prd(icount, 2)  = nface;
            isav_prd(icount, 3)  = ista;
            isav_prd(icount, 4)  = iend;
            isav_prd(icount, 5)  = jsta;
            isav_prd(icount, 6)  = jend;
            isav_prd(icount, 7)  = ksta;
            isav_prd(icount, 8)  = kend;
            isav_prd(icount, 9)  = mdim;
            isav_prd(icount, 10) = ndim;
            isav_prd(icount, 11) = nseg;
            isav_prd(icount, 12) = nblp;
        }
    }

    nbcprd = icount;

    if (nbcprd > lbcprd) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, "  Stopping in pre_period:");
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, "  nbcprd too small");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    return;
}

} // namespace pre_period_ns
