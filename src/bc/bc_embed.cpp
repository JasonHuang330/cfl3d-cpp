// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc_embed.h"
#include "bc.h"
#include "lead.h"
#include "i2x.h"
#include "i2xs.h"
#include "runtime/common_blocks.h"

namespace bc_embed_ns {

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
    bc_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, cl, nou, bou, nbuf, ibufdim,
              maxbl, maxgr, maxseg, itrans, irotat, idefrm, igridg, nblg,
              nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
              ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek, lwdat,
              myid, idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
}

void bc_embed(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
              FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
              int& nwork, int& maxbl, int& maxgr, int& lbcemb,
              FortranArray1DRef<int> iadvance, FortranArray1DRef<int> idimg,
              FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg,
              FortranArray2DRef<int> isav_emb, FortranArray1DRef<int> ireq_ar,
              FortranArray1DRef<int> ireq_snd, FortranArray1DRef<int> index_ar,
              FortranArray2DRef<int> keep_trac, FortranArray1DRef<int> keep_trac2,
              int& myid, int& myhost, int& mycomm,
              FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou,
              FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
              FortranArray2DRef<int> iviscg, FortranArray2DRef<int> istat2,
              int& istat_size, int& nummem)
{
    // COMMON block aliases
    int& level  = cmn_mgrd.level;
    int& nbcemb = cmn_is_embedbc.nbcemb;

    // Local aliases for ginfo (set by lead())
    int& jdim  = cmn_ginfo.jdim;
    int& kdim  = cmn_ginfo.kdim;
    int& idim  = cmn_ginfo.idim;
    int& lqj0  = cmn_ginfo.lqj0;
    int& lqk0  = cmn_ginfo.lqk0;
    int& lqi0  = cmn_ginfo.lqi0;
    int& lvj0  = cmn_ginfo.lvj0;
    int& lvk0  = cmn_ginfo.lvk0;
    int& lvi0  = cmn_ginfo.lvi0;
    int& lbcj  = cmn_ginfo.lbcj;
    int& lbck  = cmn_ginfo.lbck;
    int& lbci  = cmn_ginfo.lbci;

    // Check loop bounds
    if (cmn_is_embedbc.ie_emb[level - 1] < cmn_is_embedbc.is_emb[level - 1]) return;

    if (ntime > 0 && nbcemb > 0) {

        // First Case: all data needed to set embedded bc lies on the
        //             current processor

        for (int lcnt = cmn_is_embedbc.is_emb[level - 1];
                 lcnt <= cmn_is_embedbc.ie_emb[level - 1]; lcnt++) {

            // nblf is finer (embedded) block
            // nblc is coarser block in which block nblf is embedded
            int nblf    = isav_emb(lcnt, 1);
            int nblc    = isav_emb(lcnt, 9);
            int nd_recv = mblk2nd(nblf);
            int nd_srce = mblk2nd(nblc);

            if (nd_srce == myid && nd_recv == myid) {

                if (iadvance(nblf) >= 0) {

                    lead_ns::lead(nblf, lw, lw2, maxbl);

                    int n       = lcnt;
                    int nface   = isav_emb(lcnt, 2);
                    int is      = isav_emb(lcnt, 3);
                    int ie      = isav_emb(lcnt, 4);
                    int js      = isav_emb(lcnt, 5);
                    int je      = isav_emb(lcnt, 6);
                    int ks      = isav_emb(lcnt, 7);
                    int ke      = isav_emb(lcnt, 8);
                    nblc        = isav_emb(lcnt, 9);
                    int nsi     = isav_emb(lcnt, 10);
                    nd_srce     = mblk2nd(nblc);
                    int idimc   = idimg(nblc);
                    int jdimc   = jdimg(nblc);
                    int kdimc   = kdimg(nblc);

                    int maxdims = 0;
                    if (nface == 1 || nface == 2) maxdims = jdimc * kdimc;
                    if (nface == 3 || nface == 4) maxdims = kdimc * idimc;
                    if (nface == 5 || nface == 6) maxdims = jdimc * idimc;

                    int lqc  = lw(1,  nblc);
                    int lqct = lw(19, nblc);
                    int lqcv = lw(13, nblc);

                    if (nsi == 2) {

                        // full coarsening in i-direction

                        // interpolate q
                        int ldim = 5;
                        {
                            FortranArray4DRef<double> qc(&w(lqc),  jdimc, kdimc, idimc, 5);
                            FortranArray4DRef<double> qjf(&w(lqj0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qkf(&w(lqk0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qif(&w(lqi0), 1, 1, 1, 1);
                            FortranArray3DRef<double> bcjf(&w(lbcj), 1, 1, 1);
                            FortranArray3DRef<double> bckf(&w(lbck), 1, 1, 1);
                            FortranArray3DRef<double> bcif(&w(lbci), 1, 1, 1);
                            i2x_ns::i2x(jdimc, kdimc, idimc, qc,
                                        jdim, kdim, idim,
                                        qjf, qkf, qif,
                                        js, ks, is, je, ke, ie,
                                        nblc, ldim, nblf,
                                        bcjf, bckf, bcif, nface);
                        }

                        // interpolate vist3d
                        if (iviscg(nblc, 1) >= 2 || iviscg(nblc, 2) >= 2
                                                  || iviscg(nblc, 3) >= 2) {
                            ldim = 1;
                            FortranArray4DRef<double> qcv(&w(lqcv), jdimc, kdimc, idimc, 1);
                            FortranArray4DRef<double> qjf(&w(lvj0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qkf(&w(lvk0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qif(&w(lvi0), 1, 1, 1, 1);
                            FortranArray3DRef<double> bcjf(&w(lbcj), 1, 1, 1);
                            FortranArray3DRef<double> bckf(&w(lbck), 1, 1, 1);
                            FortranArray3DRef<double> bcif(&w(lbci), 1, 1, 1);
                            i2x_ns::i2x(jdimc, kdimc, idimc, qcv,
                                        jdim, kdim, idim,
                                        qjf, qkf, qif,
                                        js, ks, is, je, ke, ie,
                                        nblc, ldim, nblf,
                                        bcjf, bckf, bcif, nface);
                        }

                        // interpolate turb. data
                        if (iviscg(nblc, 1) >= 4 || iviscg(nblc, 2) >= 4
                                                  || iviscg(nblc, 3) >= 4) {
                            ldim = nummem;
                            FortranArray4DRef<double> qct(&w(lqct), jdimc, kdimc, idimc, nummem);
                            FortranArray4DRef<double> qjf(&w(lqj0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qkf(&w(lqk0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qif(&w(lqi0), 1, 1, 1, 1);
                            FortranArray3DRef<double> bcjf(&w(lbcj), 1, 1, 1);
                            FortranArray3DRef<double> bckf(&w(lbck), 1, 1, 1);
                            FortranArray3DRef<double> bcif(&w(lbci), 1, 1, 1);
                            i2x_ns::i2x(jdimc, kdimc, idimc, qct,
                                        jdim, kdim, idim,
                                        qjf, qkf, qif,
                                        js, ks, is, je, ke, ie,
                                        nblc, ldim, nblf,
                                        bcjf, bckf, bcif, nface);
                        }

                    } else if (nsi == 1) {

                        // semi coarsening in i-direction

                        // interpolate q
                        int ldim = 5;
                        {
                            FortranArray4DRef<double> qc(&w(lqc),  jdimc, kdimc, idimc, 5);
                            FortranArray4DRef<double> qjf(&w(lqj0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qkf(&w(lqk0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qif(&w(lqi0), 1, 1, 1, 1);
                            FortranArray3DRef<double> bcjf(&w(lbcj), 1, 1, 1);
                            FortranArray3DRef<double> bckf(&w(lbck), 1, 1, 1);
                            FortranArray3DRef<double> bcif(&w(lbci), 1, 1, 1);
                            i2xs_ns::i2xs(jdimc, kdimc, idimc, qc,
                                          jdim, kdim, idim,
                                          qjf, qkf, qif,
                                          js, ks, is, je, ke, ie,
                                          nblc, ldim, nblf,
                                          bcjf, bckf, bcif, nface);
                        }

                        // interpolate vist3d
                        if (iviscg(nblc, 1) >= 2 || iviscg(nblc, 2) >= 2
                                                  || iviscg(nblc, 3) >= 2) {
                            ldim = 1;
                            FortranArray4DRef<double> qcv(&w(lqcv), jdimc, kdimc, idimc, 1);
                            FortranArray4DRef<double> qjf(&w(lvj0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qkf(&w(lvk0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qif(&w(lvi0), 1, 1, 1, 1);
                            FortranArray3DRef<double> bcjf(&w(lbcj), 1, 1, 1);
                            FortranArray3DRef<double> bckf(&w(lbck), 1, 1, 1);
                            FortranArray3DRef<double> bcif(&w(lbci), 1, 1, 1);
                            i2xs_ns::i2xs(jdimc, kdimc, idimc, qcv,
                                          jdim, kdim, idim,
                                          qjf, qkf, qif,
                                          js, ks, is, je, ke, ie,
                                          nblc, ldim, nblf,
                                          bcjf, bckf, bcif, nface);
                        }

                        // interpolate turb. data
                        if (iviscg(nblc, 1) >= 4 || iviscg(nblc, 2) >= 4
                                                  || iviscg(nblc, 3) >= 4) {
                            ldim = nummem;
                            FortranArray4DRef<double> qct(&w(lqct), jdimc, kdimc, idimc, nummem);
                            FortranArray4DRef<double> qjf(&w(lqj0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qkf(&w(lqk0), 1, 1, 1, 1);
                            FortranArray4DRef<double> qif(&w(lqi0), 1, 1, 1, 1);
                            FortranArray3DRef<double> bcjf(&w(lbcj), 1, 1, 1);
                            FortranArray3DRef<double> bckf(&w(lbck), 1, 1, 1);
                            FortranArray3DRef<double> bcif(&w(lbci), 1, 1, 1);
                            i2xs_ns::i2xs(jdimc, kdimc, idimc, qct,
                                          jdim, kdim, idim,
                                          qjf, qkf, qif,
                                          js, ks, is, je, ke, ie,
                                          nblc, ldim, nblf,
                                          bcjf, bckf, bcif, nface);
                        }
                    }

                } // end if (iadvance(nblf) >= 0)

            } // end if (nd_srce == myid && nd_recv == myid)

        } // end do lcnt

        // Second Case: data needed to set periodic bc lies on another
        //              processor

    } // end if (ntime > 0 && nbcemb > 0)

    return;
}

} // namespace bc_embed_ns
