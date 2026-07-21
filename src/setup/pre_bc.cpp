// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "pre_bc.h"
#include <cstdio>
#include "pre_patch.h"
#include "pre_blockbc.h"
#include "pre_period.h"
#include "pre_embed.h"

namespace pre_bc_ns {

void pre_bc(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray2DRef<int> lw_temp, int& maxbl, int& maxgr, int& maxseg, int& ninter, int& intmax, int& nsub1, FortranArray2DRef<int> iindex, FortranArray2DRef<int> isav_pat, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, FortranArray1DRef<int> iiint1, FortranArray1DRef<int> iiint2, int& mxbli, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, FortranArray2DRef<int> nblk, int& lbcprd, FortranArray2DRef<int> isav_prd, FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk, FortranArray1DRef<int> nblg, int& lbcemb, FortranArray1DRef<int> iemg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> isav_emb, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbcidim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray1DRef<int> iadvance, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray3DRef<int> isav_pat_b, FortranArray1DRef<int> levelg, int& nblock, FortranArray2DRef<int> isav_blk)
{
    // Local variables
    int ierrflg = -1;
    int icount_pat, icount_blk, icount_prd, icount_emb;
    int icount_pat1, icount_blk1, icount_prd1, icount_emb1;
    int levl, nnn, lk, lm, nbl;

    // Aliases for COMMON block fields
    int& levt    = cmn_mgrd.levt;

    // initialize arrays required for patch boundaries
    for (nnn = 1; nnn <= intmax; nnn++) {
        for (lk = 1; lk <= 17; lk++) {
            isav_pat(nnn, lk) = 0;
        }
        for (lk = 1; lk <= nsub1; lk++) {
            for (lm = 1; lm <= 6; lm++) {
                isav_pat_b(nnn, lk, lm) = 0;
            }
        }
    }

    // initialize arrays required for 1:1 block boundaries
    for (nnn = 1; nnn <= 2*mxbli; nnn++) {
        for (lk = 1; lk <= 17; lk++) {
            isav_blk(nnn, lk) = 0;
        }
    }

    // initialize arrays required for periodic block boundaries
    for (nnn = 1; nnn <= lbcprd; nnn++) {
        for (lk = 1; lk <= 12; lk++) {
            isav_prd(nnn, lk) = 0;
        }
    }

    // initialize arrays required for embedded block boundaries
    for (nnn = 1; nnn <= lbcemb; nnn++) {
        for (lk = 1; lk <= 12; lk++) {
            isav_emb(nnn, lk) = 0;
        }
    }

    icount_pat = 0;
    icount_blk = 0;
    icount_prd = 0;
    icount_emb = 0;

    // Build FortranArray1DRef views over the COMMON block arrays for is_prd/ie_prd/is_emb/ie_emb
    // These are passed to pre_period and pre_embed
    FortranArray1DRef<int> is_prd_ref(cmn_is_perbc.is_prd, 5);
    FortranArray1DRef<int> ie_prd_ref(cmn_is_perbc.ie_prd, 5);
    FortranArray1DRef<int> is_emb_ref(cmn_is_embedbc.is_emb, 5);
    FortranArray1DRef<int> ie_emb_ref(cmn_is_embedbc.ie_emb, 5);

    for (levl = 1; levl <= levt; levl++) {
        cmn_is_blockbc.is_blk[levl-1] = icount_blk + 1;
        cmn_is_patch.is_pat[levl-1]   = icount_pat + 1;
        cmn_is_perbc.is_prd[levl-1]   = icount_prd + 1;
        cmn_is_embedbc.is_emb[levl-1] = icount_emb + 1;
        cmn_is_blockbc.ie_blk[levl-1] = cmn_is_blockbc.is_blk[levl-1] - 1;
        cmn_is_patch.ie_pat[levl-1]   = cmn_is_patch.is_pat[levl-1]   - 1;
        cmn_is_perbc.ie_prd[levl-1]   = cmn_is_perbc.is_prd[levl-1]   - 1;
        cmn_is_embedbc.ie_emb[levl-1] = cmn_is_embedbc.is_emb[levl-1] - 1;

        for (nbl = 1; nbl <= nblock; nbl++) {
            if (levl != levelg(nbl)) continue;  // go to 6909

            icount_pat1 = icount_pat;
            pre_patch_ns::pre_patch(nbl, lw, icount_pat, ninter,
                                    iindex, intmax, nsub1, isav_pat,
                                    isav_pat_b, jjmax1, kkmax1,
                                    iiint1, iiint2, maxbl, jdimg, kdimg, idimg,
                                    ierrflg);

            icount_blk1 = icount_blk;
            pre_blockbc_ns::pre_blockbc(nbl, lw, icount_blk, idimg,
                                        jdimg, kdimg, isav_blk, nblk, nbli,
                                        limblk, isva, nblon, mxbli, nou, bou, nbuf,
                                        ibufdim, myid, maxbl, ierrflg);

            icount_prd1 = icount_prd;
            pre_period_ns::pre_period(nbl, lw, lw2, icount_prd,
                                      maxbl, maxseg, lbcprd,
                                      nbcj0, nbck0, nbci0, nbcjdim, nbckdim,
                                      nbcidim, jbcinfo, kbcinfo, ibcinfo,
                                      igridg, jdimg, kdimg, idimg, isav_prd,
                                      is_prd_ref, ie_prd_ref, cmn_is_perbc.nbcprd,
                                      nou, bou, nbuf, ibufdim,
                                      bcvali, bcvalj, bcvalk, myid, nblg, maxgr,
                                      ierrflg);

            icount_emb1 = icount_emb;
            pre_embed_ns::pre_embed(nbl, lw, lw2, icount_emb,
                                    maxbl, maxseg, lbcemb,
                                    nbcj0, nbck0, nbci0, nbcjdim, nbckdim,
                                    nbcidim, jbcinfo, kbcinfo, ibcinfo, iemg,
                                    igridg, jdimg, kdimg, idimg, isav_emb,
                                    is_emb_ref, ie_emb_ref, cmn_is_embedbc.nbcemb,
                                    nou, bou, nbuf, ibufdim,
                                    myid, maxgr, ierrflg);

            if (icount_blk > icount_blk1) {
                cmn_is_blockbc.ie_blk[levl-1] = icount_blk;
            }
            if (icount_pat > icount_pat1) {
                cmn_is_patch.ie_pat[levl-1] = icount_pat;
            }
            if (icount_prd > icount_prd1) {
                cmn_is_perbc.ie_prd[levl-1] = icount_prd;
            }
            if (icount_emb > icount_emb1) {
                cmn_is_embedbc.ie_emb[levl-1] = icount_emb;
            }
            // 6909 continue
        }
    }

    return;
}

} // namespace pre_bc_ns
