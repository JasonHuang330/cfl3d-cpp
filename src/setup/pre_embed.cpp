// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "pre_embed.h"
#include "lead.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include <cstdio>
#include <algorithm>

namespace pre_embed_ns {

void pre_embed(int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, int& icount, int& maxbl, int& maxseg, int& lbcemb, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbcidim, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, FortranArray1DRef<int> iemg, FortranArray1DRef<int> igridg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<int> isav_emb, FortranArray1DRef<int> is_emb, FortranArray1DRef<int> ie_emb, int& nbcemb, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& maxgr, int& ierrflg)
{
    // local aliases for ginfo common block
    int& jdim  = cmn_ginfo.jdim;
    int& kdim  = cmn_ginfo.kdim;
    int& idim  = cmn_ginfo.idim;
    int& nblc  = cmn_ginfo.nblc;
    int& js    = cmn_ginfo.js;
    int& ks    = cmn_ginfo.ks;
    int& is    = cmn_ginfo.is;
    int& je    = cmn_ginfo.je;
    int& ke    = cmn_ginfo.ke;
    int& ie    = cmn_ginfo.ie;

    lead_ns::lead(nbl, lw, lw2, maxbl);

    // check to see if this block is embedded
    int igrid = igridg(nbl);
    // igridg(nbl) may be unset (0) during the sizing pass for the coarsest
    // block of the last grid; an unset grid index means this block's grid
    // isn't determined yet, so it cannot be an embedded block -> skip.
    // (Guards against an out-of-bounds iemg(0) read; the Fortran relies on
    // iemg(0) happening to read a benign 0 in its memory layout.)
    if (igrid < 1) return;
    int iem   = iemg(igrid);
    if (iem <= 0) return;

    // nblc is the coarser block in which the current block is embedded
    // (note: nblc,js,je,...ke are set by the call to lead)
    int jc = jdimg(nblc);
    int kc = kdimg(nblc);
    int ic = idimg(nblc);

    int nsi = (idim - 1) / (ie - is);

    if (jdim != (je - js) * 2 + 1 || je > jc || js < 1) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, " embedding error:  jdim,je,js= %d %d %d", jdim, je, js);
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    } else if (kdim != (ke - ks) * 2 + 1 || ke > kc || ks < 1) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, " embedding error:  kdim,ke,ks= %d %d %d", kdim, ke, ks);
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    } else if (ie > ic || is < 1) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, " embedding error:  idim,ie,is= %d %d %d", idim, ie, is);
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    if (nsi < 1 || nsi > 2) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, " embedding error:  idim,ie,is,nsi= %d %d %d %d", idim, ie, is, nsi);
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    // determine all the embedded type bcs for this block

    // loop 802
    for (int nseg = 1; nseg <= nbci0(nbl); nseg++) {
        if (ibcinfo(nbl, nseg, 1, 1) == 21) {
            int nface  = 1;
            icount = icount + 1;
            isav_emb(icount, 1)  = nbl;
            isav_emb(icount, 2)  = nface;
            isav_emb(icount, 3)  = is;
            isav_emb(icount, 4)  = ie;
            isav_emb(icount, 5)  = js;
            isav_emb(icount, 6)  = je;
            isav_emb(icount, 7)  = ks;
            isav_emb(icount, 8)  = ke;
            isav_emb(icount, 9)  = nblc;
            isav_emb(icount, 10) = nsi;
        }
    }

    // loop 803
    for (int nseg = 1; nseg <= nbcidim(nbl); nseg++) {
        if (ibcinfo(nbl, nseg, 1, 2) == 21) {
            int nface  = 2;
            icount = icount + 1;
            isav_emb(icount, 1)  = nbl;
            isav_emb(icount, 2)  = nface;
            isav_emb(icount, 3)  = is;
            isav_emb(icount, 4)  = ie;
            isav_emb(icount, 5)  = js;
            isav_emb(icount, 6)  = je;
            isav_emb(icount, 7)  = ks;
            isav_emb(icount, 8)  = ke;
            isav_emb(icount, 9)  = nblc;
            isav_emb(icount, 10) = nsi;
        }
    }

    // loop 804
    for (int nseg = 1; nseg <= nbcj0(nbl); nseg++) {
        if (jbcinfo(nbl, nseg, 1, 1) == 21) {
            int nface  = 3;
            icount = icount + 1;
            isav_emb(icount, 1)  = nbl;
            isav_emb(icount, 2)  = nface;
            isav_emb(icount, 3)  = is;
            isav_emb(icount, 4)  = ie;
            isav_emb(icount, 5)  = js;
            isav_emb(icount, 6)  = je;
            isav_emb(icount, 7)  = ks;
            isav_emb(icount, 8)  = ke;
            isav_emb(icount, 9)  = nblc;
            isav_emb(icount, 10) = nsi;
        }
    }

    // loop 805
    for (int nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
        if (jbcinfo(nbl, nseg, 1, 2) == 21) {
            int nface  = 4;
            icount = icount + 1;
            isav_emb(icount, 1)  = nbl;
            isav_emb(icount, 2)  = nface;
            isav_emb(icount, 3)  = is;
            isav_emb(icount, 4)  = ie;
            isav_emb(icount, 5)  = js;
            isav_emb(icount, 6)  = je;
            isav_emb(icount, 7)  = ks;
            isav_emb(icount, 8)  = ke;
            isav_emb(icount, 9)  = nblc;
            isav_emb(icount, 10) = nsi;
        }
    }

    // loop 806
    for (int nseg = 1; nseg <= nbck0(nbl); nseg++) {
        if (kbcinfo(nbl, nseg, 1, 1) == 21) {
            int nface  = 5;
            icount = icount + 1;
            isav_emb(icount, 1)  = nbl;
            isav_emb(icount, 2)  = nface;
            isav_emb(icount, 3)  = is;
            isav_emb(icount, 4)  = ie;
            isav_emb(icount, 5)  = js;
            isav_emb(icount, 6)  = je;
            isav_emb(icount, 7)  = ks;
            isav_emb(icount, 8)  = ke;
            isav_emb(icount, 9)  = nblc;
            isav_emb(icount, 10) = nsi;
        }
    }

    // loop 807
    for (int nseg = 1; nseg <= nbckdim(nbl); nseg++) {
        if (kbcinfo(nbl, nseg, 1, 2) == 21) {
            int nface  = 6;
            icount = icount + 1;
            isav_emb(icount, 1)  = nbl;
            isav_emb(icount, 2)  = nface;
            isav_emb(icount, 3)  = is;
            isav_emb(icount, 4)  = ie;
            isav_emb(icount, 5)  = js;
            isav_emb(icount, 6)  = je;
            isav_emb(icount, 7)  = ks;
            isav_emb(icount, 8)  = ke;
            isav_emb(icount, 9)  = nblc;
            isav_emb(icount, 10) = nsi;
        }
    }

    nbcemb = icount;

    if (nbcemb > lbcemb) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, "  Stopping in pre_embed:");
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, "  nbcemb too small");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    return;
}

} // namespace pre_embed_ns
