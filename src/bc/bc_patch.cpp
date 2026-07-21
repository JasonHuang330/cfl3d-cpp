// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc_patch.h"
#include <cstdio>
#include "lead.h"
#include "int2.h"
#include "termn8.h"
#include <cmath>
#include <algorithm>

namespace bc_patch_ns {

// bc: standard CFL3D bc dispatcher for the bc_patch module.
// This module does not contain a bc dispatcher in the Fortran source;
// the empty body is required by the header interface.
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
    return;
}

void bc_patch(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
              FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork,
              int& maxbl, int& maxgr, int& intmax, int& nsub1, int& maxxe,
              FortranArray1DRef<int> iadvance, FortranArray1DRef<int> jdimg,
              FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg,
              int& ninter, FortranArray2DRef<double> windex, FortranArray2DRef<int> iindex,
              FortranArray1DRef<int> nblkpt, FortranArray2DRef<double> dthetxx,
              FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz,
              FortranArray2DRef<int> isav_pat, FortranArray3DRef<int> isav_pat_b,
              FortranArray1DRef<int> ireq_ar, FortranArray1DRef<int> ireq_snd,
              FortranArray1DRef<int> index_ar, FortranArray3DRef<int> keep_trac,
              FortranArray2DRef<int> keep_trac2, int& myid, int& myhost, int& mycomm,
              FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou,
              FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
              FortranArray2DRef<int> istat2, int& istat_size, int& nummem)
{
    // COMMON block references
    int32_t& icyc      = cmn_mgrd.icyc;
    int32_t& nitfo     = cmn_info.nitfo;
    int32_t& level     = cmn_mgrd.level;
    int32_t& ivmx      = cmn_maxiv.ivmx;
    int32_t* is_pat    = cmn_is_patch.is_pat;   // 1-based: is_pat[0..4]
    int32_t* ie_pat    = cmn_is_patch.ie_pat;   // 1-based: ie_pat[0..4]
    int32_t& ipatch1st = cmn_is_patch.ipatch1st;

    // lim_ptch = 0 for patching without limiting; = 1 with limiting
    int lim_ptch = 0;

    if (ntime > 0 && std::abs(ninter) > 0) {

        // ***********************************************************************
        //    First Case: some or all data needed to set patch bc lies on
        //                the current processor
        // ***********************************************************************

        // first order interpolation for q if icyc .le. nitfo
        // or if ipatch1st is specified as 1
        int ifo = 0;
        if (icyc > nitfo && ipatch1st != 1) ifo = 1;

        // Loop over patch entries for this level
        // Fortran: is_pat(level) to ie_pat(level) — 1-based array, level is 1-based index
        for (int lcnt = is_pat[level - 1]; lcnt <= ie_pat[level - 1]; lcnt++) {
            // nb_to is current (to) block
            // nb_from is neighbor (from) block
            int nb_to   = isav_pat(lcnt, 1);
            int nfb     = isav_pat(lcnt, 2);
            int nd_recv = mblk2nd(nb_to);

            if (iadvance(nb_to) >= 0) {
                for (int l = 1; l <= nfb; l++) {
                    int nb_from = isav_pat_b(lcnt, l, 1);
                    int nd_srce = mblk2nd(nb_from);

                    if (nd_srce == myid && nd_recv == myid) {

                        lead_ns::lead(nb_to, lw, lw2, maxbl);

                        int j21    = isav_pat(lcnt, 3);
                        int j22    = isav_pat(lcnt, 4);
                        int k21    = isav_pat(lcnt, 5);
                        int k22    = isav_pat(lcnt, 6);
                        int lqedge = isav_pat(lcnt, 7);
                        int lqedgt = isav_pat(lcnt, 8);
                        int lqedgv = isav_pat(lcnt, 9);
                        int lqedgb = isav_pat(lcnt, 10);
                        int jmax2  = isav_pat(lcnt, 11);
                        int kmax2  = isav_pat(lcnt, 12);
                        int mint1  = isav_pat(lcnt, 13);
                        int mint2  = isav_pat(lcnt, 14);
                        int lst    = isav_pat(lcnt, 15);
                        int npt    = isav_pat(lcnt, 16);
                        int icheck = isav_pat(lcnt, 17);
                        int mtype  = isav_pat_b(lcnt, l, 2);
                        int jmax1  = isav_pat_b(lcnt, l, 3);
                        int kmax1  = isav_pat_b(lcnt, l, 4);
                        int idim1  = idimg(nb_from);
                        int jdim1  = jdimg(nb_from);
                        int kdim1  = kdimg(nb_from);
                        int iint1  = isav_pat_b(lcnt, l, 5);
                        int iint2  = isav_pat_b(lcnt, l, 6);
                        int lqintl = lw(1, nb_from);
                        int lqintt = lw(19, nb_from);
                        int lqintv = lw(13, nb_from);
                        double dthtx = dthetxx(icheck, l);
                        double dthty = dthetyy(icheck, l);
                        double dthtz = dthetzz(icheck, l);

                        int lwk1   = 1;
                        int lwk2   = lwk1 + jmax1 * kmax1;
                        int lwk3   = lwk2 + jmax1 * kmax1;
                        int lwk4   = lwk3 + jmax1 * kmax1;
                        int klimit = lwk4 + jdim1 * kdim1 * idim1 * 5 - 1;

                        if (klimit > nwork) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                " stopping in bc_patch... insufficient wk storage");
                            int termn8_err = -1;
                            termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                        }

                        // interpolate q data
                        int ldim   = 5;
                        int iq     = 1;
                        int jkldim = ldim * jmax1 * kmax1;

                        // Build array refs from w/wk sub-arrays for int2 calls
                        // q1: FortranArray4DRef<double> from w(lqintl)
                        // bc2: FortranArray3DRef<double> from w(lqedgb)
                        // q2: FortranArray4DRef<double> from w(lqedge)
                        // xie2: FortranArray1DRef<double> from windex(lst,1)
                        // eta2: FortranArray1DRef<double> from windex(lst,2)
                        // q1g: FortranArray2DRef<double> from wk(lwk1)
                        // q2g: FortranArray2DRef<double> from wk(lwk2)
                        // temp: FortranArray1DRef<double> from wk(lwk3)
                        // nblkpt_sub: FortranArray1DRef<int> from nblkpt(lst)
                        // q1wk: FortranArray4DRef<double> from wk(lwk4)

                        // interpolate to first layer of ghost cells (q)
                        {
                            FortranArray4DRef<double> q1_ref(&w(lqintl), jdim1, kdim1, idim1, 5);
                            FortranArray3DRef<double> bc2_ref(&w(lqedgb), jmax2, kmax2, 2);
                            FortranArray4DRef<double> q2_ref(&w(lqedge), jmax2, kmax2, ldim, 4);
                            FortranArray1DRef<double> xie2_ref(&windex(lst, 1), npt);
                            FortranArray1DRef<double> eta2_ref(&windex(lst, 2), npt);
                            FortranArray2DRef<double> q1g_ref(&wk(lwk1), jmax1, kmax1);
                            FortranArray2DRef<double> q2g_ref(&wk(lwk2), jmax1, kmax1);
                            FortranArray1DRef<double> temp_ref(&wk(lwk3), jmax1 * kmax1);
                            FortranArray1DRef<int>    nblkpt_ref(&nblkpt(lst), npt);
                            FortranArray4DRef<double> q1wk_ref(&wk(lwk4), jdim1, kdim1, idim1, 5);

                            int2_ns::int2(iq, jdim1, kdim1, idim1, jmax2, kmax2,
                                          nsub1, jmax1, kmax1, l, iint1,
                                          q1_ref, bc2_ref, mint1, q2_ref,
                                          xie2_ref, eta2_ref,
                                          q1g_ref, q2g_ref, temp_ref,
                                          nblkpt_ref, intmax, icheck, mtype,
                                          iindex, ifo, ldim, npt,
                                          j21, j22, k21, k22,
                                          q1wk_ref, dthtx, dthty, dthtz, lim_ptch);
                        }

                        // interpolate to second layer of ghost cells (q)
                        {
                            FortranArray4DRef<double> q1_ref(&w(lqintl), jdim1, kdim1, idim1, 5);
                            FortranArray3DRef<double> bc2_ref(&w(lqedgb), jmax2, kmax2, 2);
                            FortranArray4DRef<double> q2_ref(&w(lqedge), jmax2, kmax2, ldim, 4);
                            FortranArray1DRef<double> xie2_ref(&windex(lst, 1), npt);
                            FortranArray1DRef<double> eta2_ref(&windex(lst, 2), npt);
                            FortranArray2DRef<double> q1g_ref(&wk(lwk1), jmax1, kmax1);
                            FortranArray2DRef<double> q2g_ref(&wk(lwk2), jmax1, kmax1);
                            FortranArray1DRef<double> temp_ref(&wk(lwk3), jmax1 * kmax1);
                            FortranArray1DRef<int>    nblkpt_ref(&nblkpt(lst), npt);
                            FortranArray4DRef<double> q1wk_ref(&wk(lwk4), jdim1, kdim1, idim1, 5);

                            int2_ns::int2(iq, jdim1, kdim1, idim1, jmax2, kmax2,
                                          nsub1, jmax1, kmax1, l, iint2,
                                          q1_ref, bc2_ref, mint2, q2_ref,
                                          xie2_ref, eta2_ref,
                                          q1g_ref, q2g_ref, temp_ref,
                                          nblkpt_ref, intmax, icheck, mtype,
                                          iindex, ifo, ldim, npt,
                                          j21, j22, k21, k22,
                                          q1wk_ref, dthtx, dthty, dthtz, lim_ptch);
                        }

                        // interpolate vist3d data
                        if (ivmx >= 2) {
                            ldim   = 1;
                            iq     = 0;
                            jkldim = ldim * jmax1 * kmax1;

                            // interpolate to first layer of ghost cells (vist3d)
                            {
                                FortranArray4DRef<double> q1_ref(&w(lqintv), jdim1, kdim1, idim1, 1);
                                FortranArray3DRef<double> bc2_ref(&w(lqedgb), jmax2, kmax2, 2);
                                FortranArray4DRef<double> q2_ref(&w(lqedgv), jmax2, kmax2, ldim, 4);
                                FortranArray1DRef<double> xie2_ref(&windex(lst, 1), npt);
                                FortranArray1DRef<double> eta2_ref(&windex(lst, 2), npt);
                                FortranArray2DRef<double> q1g_ref(&wk(lwk1), jmax1, kmax1);
                                FortranArray2DRef<double> q2g_ref(&wk(lwk2), jmax1, kmax1);
                                FortranArray1DRef<double> temp_ref(&wk(lwk3), jmax1 * kmax1);
                                FortranArray1DRef<int>    nblkpt_ref(&nblkpt(lst), npt);
                                FortranArray4DRef<double> q1wk_ref(&wk(lwk4), jdim1, kdim1, idim1, 5);

                                int2_ns::int2(iq, jdim1, kdim1, idim1, jmax2, kmax2,
                                              nsub1, jmax1, kmax1, l, iint1,
                                              q1_ref, bc2_ref, mint1, q2_ref,
                                              xie2_ref, eta2_ref,
                                              q1g_ref, q2g_ref, temp_ref,
                                              nblkpt_ref, intmax, icheck, mtype,
                                              iindex, ifo, ldim, npt,
                                              j21, j22, k21, k22,
                                              q1wk_ref, dthtx, dthty, dthtz, lim_ptch);
                            }

                            // interpolate to second layer of ghost cells (vist3d)
                            {
                                FortranArray4DRef<double> q1_ref(&w(lqintv), jdim1, kdim1, idim1, 1);
                                FortranArray3DRef<double> bc2_ref(&w(lqedgb), jmax2, kmax2, 2);
                                FortranArray4DRef<double> q2_ref(&w(lqedgv), jmax2, kmax2, ldim, 4);
                                FortranArray1DRef<double> xie2_ref(&windex(lst, 1), npt);
                                FortranArray1DRef<double> eta2_ref(&windex(lst, 2), npt);
                                FortranArray2DRef<double> q1g_ref(&wk(lwk1), jmax1, kmax1);
                                FortranArray2DRef<double> q2g_ref(&wk(lwk2), jmax1, kmax1);
                                FortranArray1DRef<double> temp_ref(&wk(lwk3), jmax1 * kmax1);
                                FortranArray1DRef<int>    nblkpt_ref(&nblkpt(lst), npt);
                                FortranArray4DRef<double> q1wk_ref(&wk(lwk4), jdim1, kdim1, idim1, 5);

                                int2_ns::int2(iq, jdim1, kdim1, idim1, jmax2, kmax2,
                                              nsub1, jmax1, kmax1, l, iint2,
                                              q1_ref, bc2_ref, mint2, q2_ref,
                                              xie2_ref, eta2_ref,
                                              q1g_ref, q2g_ref, temp_ref,
                                              nblkpt_ref, intmax, icheck, mtype,
                                              iindex, ifo, ldim, npt,
                                              j21, j22, k21, k22,
                                              q1wk_ref, dthtx, dthty, dthtz, lim_ptch);
                            }
                        } // end if (ivmx >= 2)

                        // interpolate turb. data
                        if (ivmx >= 4) {
                            ldim   = nummem;
                            iq     = 0;
                            jkldim = ldim * jmax1 * kmax1;

                            // interpolate to first layer of ghost cells (turb)
                            {
                                FortranArray4DRef<double> q1_ref(&w(lqintt), jdim1, kdim1, idim1, nummem);
                                FortranArray3DRef<double> bc2_ref(&w(lqedgb), jmax2, kmax2, 2);
                                FortranArray4DRef<double> q2_ref(&w(lqedgt), jmax2, kmax2, ldim, 4);
                                FortranArray1DRef<double> xie2_ref(&windex(lst, 1), npt);
                                FortranArray1DRef<double> eta2_ref(&windex(lst, 2), npt);
                                FortranArray2DRef<double> q1g_ref(&wk(lwk1), jmax1, kmax1);
                                FortranArray2DRef<double> q2g_ref(&wk(lwk2), jmax1, kmax1);
                                FortranArray1DRef<double> temp_ref(&wk(lwk3), jmax1 * kmax1);
                                FortranArray1DRef<int>    nblkpt_ref(&nblkpt(lst), npt);
                                FortranArray4DRef<double> q1wk_ref(&wk(lwk4), jdim1, kdim1, idim1, 5);

                                int2_ns::int2(iq, jdim1, kdim1, idim1, jmax2, kmax2,
                                              nsub1, jmax1, kmax1, l, iint1,
                                              q1_ref, bc2_ref, mint1, q2_ref,
                                              xie2_ref, eta2_ref,
                                              q1g_ref, q2g_ref, temp_ref,
                                              nblkpt_ref, intmax, icheck, mtype,
                                              iindex, ifo, ldim, npt,
                                              j21, j22, k21, k22,
                                              q1wk_ref, dthtx, dthty, dthtz, lim_ptch);
                            }

                            // interpolate to second layer of ghost cells (turb)
                            {
                                FortranArray4DRef<double> q1_ref(&w(lqintt), jdim1, kdim1, idim1, nummem);
                                FortranArray3DRef<double> bc2_ref(&w(lqedgb), jmax2, kmax2, 2);
                                FortranArray4DRef<double> q2_ref(&w(lqedgt), jmax2, kmax2, ldim, 4);
                                FortranArray1DRef<double> xie2_ref(&windex(lst, 1), npt);
                                FortranArray1DRef<double> eta2_ref(&windex(lst, 2), npt);
                                FortranArray2DRef<double> q1g_ref(&wk(lwk1), jmax1, kmax1);
                                FortranArray2DRef<double> q2g_ref(&wk(lwk2), jmax1, kmax1);
                                FortranArray1DRef<double> temp_ref(&wk(lwk3), jmax1 * kmax1);
                                FortranArray1DRef<int>    nblkpt_ref(&nblkpt(lst), npt);
                                FortranArray4DRef<double> q1wk_ref(&wk(lwk4), jdim1, kdim1, idim1, 5);

                                int2_ns::int2(iq, jdim1, kdim1, idim1, jmax2, kmax2,
                                              nsub1, jmax1, kmax1, l, iint2,
                                              q1_ref, bc2_ref, mint2, q2_ref,
                                              xie2_ref, eta2_ref,
                                              q1g_ref, q2g_ref, temp_ref,
                                              nblkpt_ref, intmax, icheck, mtype,
                                              iindex, ifo, ldim, npt,
                                              j21, j22, k21, k22,
                                              q1wk_ref, dthtx, dthty, dthtz, lim_ptch);
                            }
                        } // end if (ivmx >= 4)

                    } // end if (nd_srce == myid && nd_recv == myid)

                } // end do l = 1, nfb

            } // end if (iadvance(nb_to) >= 0)

        } // end do lcnt

        // ***********************************************************************
        //    Second Case: some or all data needed to set patch bc lies on
        //                 another processor
        // ***********************************************************************

    } // end if (ntime > 0 && abs(ninter) > 0)

    return;
}

} // namespace bc_patch_ns
