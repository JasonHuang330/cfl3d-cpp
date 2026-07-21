// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "ae_corr.h"
#include "lead.h"
#include "genforce.h"
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"

namespace ae_corr_ns {

void ae_corr(FortranArray3DRef<double> stm, FortranArray3DRef<double> stmi,
             FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn,
             FortranArray2DRef<double> gforcn, FortranArray2DRef<double> gforcs,
             FortranArray2DRef<double> gforcnm, FortranArray2DRef<double> gf0,
             FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
             FortranArray1DRef<double> w, int& mgwk, int& maxbl, int& maxseg,
             FortranArray2DRef<double> aesrfdat, int& nmds, int& maxaes,
             int& nt, FortranArray1DRef<int> mblk2nd, int& iseqr,
             FortranArray1DRef<int> levelg, FortranArray1DRef<int> iadvance,
             int& nblock, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
             FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
             FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
             int& myid, FortranArray1DRef<int> nsegdfrm,
             FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf,
             FortranArray3DRef<double> perturb, FortranArray4DRef<double> aehist,
             int& ncycmax, int& maxsegdg, int& myhost, int& mycomm)
{
    // COMMON block aliases
    int& jdim    = cmn_ginfo.jdim;
    int& kdim    = cmn_ginfo.kdim;
    int& idim    = cmn_ginfo.idim;
    int& lsk     = cmn_ginfo.lsk;
    int& lsj     = cmn_ginfo.lsj;
    int& lsi     = cmn_ginfo.lsi;
    int& lbcj    = cmn_ginfo.lbcj;
    int& lbck    = cmn_ginfo.lbck;
    int& lbci    = cmn_ginfo.lbci;
    int& lblk    = cmn_ginfo.lblk;
    int& lqj0    = cmn_ginfo.lqj0;
    int& lqk0    = cmn_ginfo.lqk0;
    int& lqi0    = cmn_ginfo.lqi0;
    int& lxmdj   = cmn_ginfo.lxmdj;
    int& lxmdk   = cmn_ginfo.lxmdk;
    int& lxmdi   = cmn_ginfo.lxmdi;
    int& lglobal = cmn_mgrd.lglobal;
    int& naesrf  = cmn_elastic.naesrf;
    int& ntt     = cmn_info.ntt;
    // levelt is int32_t[5] in cmn_info (1-based, so levelt[iseqr-1])

    //
    // compute aeroelastic generalized forces
    //
    for (int iaes = 1; iaes <= naesrf; iaes++) {
        int    iskyhk = (int)aesrfdat(1, iaes);
        double grefl  = aesrfdat(2, iaes);
        double uinf   = aesrfdat(3, iaes);
        double qinf   = aesrfdat(4, iaes);
        int    nmodes = (int)aesrfdat(5, iaes);

        for (int n = 2; n <= 2 * nmodes; n += 2) {
            int n2 = n / 2;
            double cx = 0.0;
            double cy = 0.0;
            double cz = 0.0;

            for (int nbl = 1; nbl <= nblock; nbl++) {
                double cxb = 0.0;
                double cyb = 0.0;
                double czb = 0.0;

                if (iadvance(nbl) >= 0 &&
                    levelg(nbl)   >= lglobal &&
                    levelg(nbl)   <= cmn_info.levelt[iseqr - 1]) {

                    if (myid == mblk2nd(nbl)) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);

                        // Wrap work array slices for genforce
                        FortranArray4DRef<double> wsk  (&w(lsk),  jdim, kdim, idim, 5);
                        FortranArray4DRef<double> wsj  (&w(lsj),  jdim, kdim, idim, 5);
                        FortranArray4DRef<double> wsi  (&w(lsi),  jdim, kdim, idim, 5);
                        FortranArray3DRef<double> wbcj (&w(lbcj), jdim, kdim, idim);
                        FortranArray3DRef<double> wbck (&w(lbck), jdim, kdim, idim);
                        FortranArray3DRef<double> wbci (&w(lbci), jdim, kdim, idim);
                        FortranArray3DRef<double> wblk (&w(lblk), jdim, kdim, idim);
                        FortranArray4DRef<double> wqj0 (&w(lqj0), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> wqk0 (&w(lqk0), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> wqi0 (&w(lqi0), jdim, kdim, idim, 5);
                        FortranArray4DRef<double> wxmdj(&w(lxmdj), jdim, kdim, idim, 3);
                        FortranArray4DRef<double> wxmdk(&w(lxmdk), jdim, kdim, idim, 3);
                        FortranArray4DRef<double> wxmdi(&w(lxmdi), jdim, kdim, idim, 3);

                        genforce_ns::genforce(
                            jdim, kdim, idim,
                            wsk, wsj, wsi,
                            czb, cyb, cxb,
                            wbcj, wbck, wbci, wblk,
                            nbl,
                            wqj0, wqk0, wqi0,
                            maxbl, maxseg, n2,
                            wxmdj, wxmdk, wxmdi,
                            aesrfdat, nmds, maxaes,
                            maxsegdg, nsegdfrm,
                            jcsi, jcsf, kcsi, kcsf,
                            icsi, icsf,
                            idfrmseg, iaes, iaesurf);
                    }

                    if (myid == myhost) {
                        cx += cxb;
                        cy += cyb;
                        cz += czb;
                    }
                }
            } // nbl

            if (iskyhk == 1 && nt == 1) {
                gf0(n, iaes) = -qinf * grefl * grefl * (cx + cy + cz) + gf0(n, iaes);
            }
            gforcs(n - 1, iaes) = 0.0;
            gforcs(n,     iaes) = qinf * grefl * grefl * (cx + cy + cz) - gf0(n, iaes);
        } // n
    } // iaes

    //
    // modal displacement and velocity correction
    //
    for (int iaes = 1; iaes <= naesrf; iaes++) {
        int nmodes = (int)aesrfdat(5, iaes);
        for (int n = 1; n <= 2 * nmodes; n++) {
            // don't update if the modal time variation is specified
            double moddfl = perturb((n + 1) / 2, iaes, 1);
            if (moddfl == 0) {
                xs(n, iaes) = 0.0;
                for (int j = 1; j <= 2 * nmodes; j++) {
                    xs(n, iaes) = xs(n, iaes)
                        + stm(n, j, iaes) * xxn(j, iaes)
                        + 0.5 * stmi(n, j, iaes) * (gforcs(j, iaes) + gforcn(j, iaes));
                }
            }
        }
    }

    //
    // update variables (set current values to old for next time step)
    //
    for (int iaes = 1; iaes <= naesrf; iaes++) {
        int nmodes = (int)aesrfdat(5, iaes);
        for (int n = 1; n <= 2 * nmodes; n++) {
            gforcnm(n, iaes) = gforcn(n, iaes);
            gforcn(n, iaes)  = gforcs(n, iaes);
            xxn(n, iaes)     = xs(n, iaes);
        }
    }

    //
    // store off the variables into the time-history array
    //
    if (myid == myhost) {
        for (int iaes = 1; iaes <= naesrf; iaes++) {
            int nmodes = (int)aesrfdat(5, iaes);
            for (int n = 1; n <= nmodes; n++) {
                aehist(ntt, 1, n, iaes) = xs(2 * n - 1, iaes);
                aehist(ntt, 2, n, iaes) = xs(2 * n,     iaes);
                aehist(ntt, 3, n, iaes) = gforcs(2 * n, iaes);
            }
        }
    }
}

} // namespace ae_corr_ns
