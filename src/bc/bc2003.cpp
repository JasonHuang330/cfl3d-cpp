// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2003.h"
#include "rie1de.h"
#include <cmath>
#include <cstring>
#include <algorithm>

namespace bc2003_ns {

// -----------------------------------------------------------------------
// bc() dispatcher — standard CFL3D pattern for bctype 2003
// -----------------------------------------------------------------------
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

        // Build array views from w using lw offsets
        FortranArray4DRef<double> q    (&w(lw(1,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray4DRef<double> qj0  (&w(lw(2,ibl)),  kdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qk0  (&w(lw(3,ibl)),  jdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qi0  (&w(lw(4,ibl)),  jdim,   kdim,   5,      4);
        FortranArray4DRef<double> sj   (&w(lw(5,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> sk   (&w(lw(6,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> si   (&w(lw(7,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray3DRef<double> bcj  (&w(lw(8,ibl)),  kdim,   idim-1, 2);
        FortranArray3DRef<double> bck  (&w(lw(9,ibl)),  jdim,   idim-1, 2);
        FortranArray3DRef<double> bci  (&w(lw(10,ibl)), jdim,   kdim,   2);

        // Build array views from wk using lw2 offsets
        FortranArray4DRef<double> tursav(&wk(lw2(1,ibl)), jdim,   kdim,   idim,   nummem);
        FortranArray4DRef<double> tj0   (&wk(lw2(2,ibl)), kdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> tk0   (&wk(lw2(3,ibl)), jdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> ti0   (&wk(lw2(4,ibl)), jdim,   kdim,   nummem, 4);
        FortranArray3DRef<double> vist3d(&wk(lw2(5,ibl)), jdim,   kdim,   idim);
        FortranArray4DRef<double> vj0   (&wk(lw2(6,ibl)), kdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vk0   (&wk(lw2(7,ibl)), jdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vi0   (&wk(lw2(8,ibl)), jdim,   kdim,   1,      4);

        // jvdim for t array: max of (jdim-1)*(kdim-1) and (jdim-1)*(idim-1) and (kdim-1)*(idim-1)
        int jvdim_val = std::max({(jdim-1)*(kdim-1), (jdim-1)*(idim-1), (kdim-1)*(idim-1)});
        FortranArray2D<double> t_arr(jvdim_val, 25);
        FortranArray2DRef<double> t = t_arr.ref();

        // --- j-face segments ---
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl,iseg,1,1) != 2003) continue;
            int ista  = jbcinfo(ibl,iseg,1,2);
            int iend  = jbcinfo(ibl,iseg,1,3);
            int jsta  = jbcinfo(ibl,iseg,1,4);
            int jend  = jbcinfo(ibl,iseg,1,5);
            int ksta  = jbcinfo(ibl,iseg,1,6);
            int kend  = jbcinfo(ibl,iseg,1,7);
            int nface = jbcinfo(ibl,iseg,1,8);
            int mdim  = lwdat(ibl,iseg,3);
            int ndim  = lwdat(ibl,iseg,4);
            int ldata = lw2(9,ibl);
            FortranArray4DRef<double> bcdata(&wk(ldata), mdim, ndim, 2, 12);
            char filname[80];
            int ifile = bcfilej(ibl,iseg,1);
            std::memcpy(filname, bcfiles(ifile), 80);
            bc2003(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si, bcj, bck, bci,
                   ista, iend, jsta, jend, ksta, kend, t, jvdim_val, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   mdim, ndim, bcdata, filname, iuns,
                   nou, bou, nbuf, ibufdim, myid, nummem);
        }

        // --- k-face segments ---
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) != 2003) continue;
            int ista  = kbcinfo(ibl,iseg,1,2);
            int iend  = kbcinfo(ibl,iseg,1,3);
            int jsta  = kbcinfo(ibl,iseg,1,4);
            int jend  = kbcinfo(ibl,iseg,1,5);
            int ksta  = kbcinfo(ibl,iseg,1,6);
            int kend  = kbcinfo(ibl,iseg,1,7);
            int nface = kbcinfo(ibl,iseg,1,8);
            int mdim  = lwdat(ibl,iseg,3);
            int ndim  = lwdat(ibl,iseg,4);
            int ldata = lw2(9,ibl);
            FortranArray4DRef<double> bcdata(&wk(ldata), mdim, ndim, 2, 12);
            char filname[80];
            int ifile = bcfilek(ibl,iseg,1);
            std::memcpy(filname, bcfiles(ifile), 80);
            bc2003(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si, bcj, bck, bci,
                   ista, iend, jsta, jend, ksta, kend, t, jvdim_val, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   mdim, ndim, bcdata, filname, iuns,
                   nou, bou, nbuf, ibufdim, myid, nummem);
        }

        // --- i-face segments ---
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) != 2003) continue;
            int ista  = ibcinfo(ibl,iseg,1,2);
            int iend  = ibcinfo(ibl,iseg,1,3);
            int jsta  = ibcinfo(ibl,iseg,1,4);
            int jend  = ibcinfo(ibl,iseg,1,5);
            int ksta  = ibcinfo(ibl,iseg,1,6);
            int kend  = ibcinfo(ibl,iseg,1,7);
            int nface = ibcinfo(ibl,iseg,1,8);
            int mdim  = lwdat(ibl,iseg,3);
            int ndim  = lwdat(ibl,iseg,4);
            int ldata = lw2(9,ibl);
            FortranArray4DRef<double> bcdata(&wk(ldata), mdim, ndim, 2, 12);
            char filname[80];
            int ifile = bcfilei(ibl,iseg,1);
            std::memcpy(filname, bcfiles(ifile), 80);
            bc2003(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si, bcj, bck, bci,
                   ista, iend, jsta, jend, ksta, kend, t, jvdim_val, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   mdim, ndim, bcdata, filname, iuns,
                   nou, bou, nbuf, ibufdim, myid, nummem);
        }
    } // end ibl loop
}

// -----------------------------------------------------------------------
// bc2003() — characteristic inlet BC for engine inlet
// -----------------------------------------------------------------------
void bc2003(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            FortranArray2DRef<double> t, int& jvdim, int& nface,
            FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0,
            int& mdim, int& ndim, FortranArray4DRef<double> bcdata,
            char (&filname)[80], int& iuns,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid, int& nummem)
{
    // COMMON block aliases
    float& gamma    = cmn_fluid.gamma;
    float& gm1      = cmn_fluid.gm1;
    int&   level    = cmn_mgrd.level;
    int&   lglobal  = cmn_mgrd.lglobal;
    float* tur10    = cmn_ivals.tur10;   // 0-based: tur10[l-1]
    int*   ivisc    = cmn_reyue.ivisc;   // 0-based: ivisc[0..2]
    float& radtodeg = cmn_conversion.radtodeg;

    // Local dimension variables
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    // This bc makes use of only one plane of data
    int ip = 1;


    // ******************************************************************
    // j=1 boundary   inflow for engine   bctype 2003
    // ******************************************************************
    if (nface == 3) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,6) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            int js = (i - ista) * (kend - ksta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    t(js+k-ksta, l) = q(1, k, i, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    t(js+k-ksta, 5+l) = -sj(1, k, i, l);
                }
            }
            for (int k = ksta; k <= kend1; k++) {
                int kk = k - ksta + 1;
                double xme   = bcdata(kk, ii, ip, 1);
                double pte   = bcdata(kk, ii, ip, 2);
                double tte   = bcdata(kk, ii, ip, 3);
                double alpe  = bcdata(kk, ii, ip, 4);
                double betae = bcdata(kk, ii, ip, 5);

                alpe  = alpe  / (double)radtodeg;
                betae = betae / (double)radtodeg;
                double te = tte / (1.0 + 0.5 * (double)gm1 * xme * xme);
                double ae = std::sqrt(te);
                double ue = ae * xme * std::cos(alpe) * std::cos(betae);
                double ve = -ae * xme * std::sin(betae);
                double we = ae * xme * std::sin(alpe) * std::cos(betae);

                t(js+k-ksta, 21) = ue;
                t(js+k-ksta, 22) = ve;
                t(js+k-ksta, 23) = we;
                t(js+k-ksta, 24) = ae;
                t(js+k-ksta, 25) = pte;
                t(js+k-ksta, 20) = -sj(1, k, i, 5);
            }
        }

        int jv = (kend - ksta) * (iend - ista);
        rie1de_ns::rie1de(jvdim, t, jv);

        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista) * (kend - ksta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    qj0(k, i, l, 1) = t(k-ksta+js, l);
                    qj0(k, i, l, 2) = qj0(k, i, l, 1);
                    bcj(k, i, 1)    = 0.0;
                }
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k, i, 1, 1) = vist3d(1, k, i);
                    vj0(k, i, 1, 2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            int kk = k - ksta + 1;
                            double ubar = -(qj0(k,i,2,1)*sj(1,k,i,1) +
                                            qj0(k,i,3,1)*sj(1,k,i,2) +
                                            qj0(k,i,4,1)*sj(1,k,i,3));
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(kk, ii, ip, 5+l);
                            if ((float)ubar < 0.f) {
                                tj0(k, i, l, 1) = t1;
                                tj0(k, i, l, 2) = t1;
                            } else {
                                tj0(k, i, l, 1) = tursav(1, k, i, l);
                                tj0(k, i, l, 2) = tj0(k, i, l, 1);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==3


    // ******************************************************************
    // j=jdim boundary   inflow for engine   bctype 2003
    // ******************************************************************
    if (nface == 4) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,6) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            int js = (i - ista) * (kend - ksta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    t(js+k-ksta, l) = q(jdim1, k, i, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    t(js+k-ksta, 5+l) = sj(jdim, k, i, l);
                }
            }
            for (int k = ksta; k <= kend1; k++) {
                int kk = k - ksta + 1;
                double xme   = bcdata(kk, ii, ip, 1);
                double pte   = bcdata(kk, ii, ip, 2);
                double tte   = bcdata(kk, ii, ip, 3);
                double alpe  = bcdata(kk, ii, ip, 4);
                double betae = bcdata(kk, ii, ip, 5);

                alpe  = alpe  / (double)radtodeg;
                betae = betae / (double)radtodeg;
                double te = tte / (1.0 + 0.5 * (double)gm1 * xme * xme);
                double ae = std::sqrt(te);
                double ue = ae * xme * std::cos(alpe) * std::cos(betae);
                double ve = -ae * xme * std::sin(betae);
                double we = ae * xme * std::sin(alpe) * std::cos(betae);

                t(js+k-ksta, 21) = ue;
                t(js+k-ksta, 22) = ve;
                t(js+k-ksta, 23) = we;
                t(js+k-ksta, 24) = ae;
                t(js+k-ksta, 25) = pte;
                t(js+k-ksta, 20) = sj(jdim, k, i, 5);
            }
        }

        int jv = (kend - ksta) * (iend - ista);
        rie1de_ns::rie1de(jvdim, t, jv);

        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista) * (kend - ksta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    qj0(k, i, l, 3) = t(k-ksta+js, l);
                    qj0(k, i, l, 4) = qj0(k, i, l, 3);
                    bcj(k, i, 2)    = 0.0;
                }
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k, i, 1, 3) = vist3d(jdim-1, k, i);
                    vj0(k, i, 1, 4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            int kk = k - ksta + 1;
                            double ubar = qj0(k,i,2,3)*sj(jdim,k,i,1) +
                                          qj0(k,i,3,3)*sj(jdim,k,i,2) +
                                          qj0(k,i,4,3)*sj(jdim,k,i,3);
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(kk, ii, ip, 5+l);
                            if ((float)ubar < 0.f) {
                                tj0(k, i, l, 3) = t1;
                                tj0(k, i, l, 4) = t1;
                            } else {
                                tj0(k, i, l, 3) = tursav(jdim-1, k, i, l);
                                tj0(k, i, l, 4) = tj0(k, i, l, 3);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==4


    // ******************************************************************
    // k=1 boundary   inflow for engine   bctype 2003
    // ******************************************************************
    if (nface == 5) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,6) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            int js = (i - ista) * (jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    t(js+j-jsta, l) = q(j, 1, i, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    t(js+j-jsta, 5+l) = -sk(j, 1, i, l);
                }
            }
            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;
                double xme   = bcdata(jj, ii, ip, 1);
                double pte   = bcdata(jj, ii, ip, 2);
                double tte   = bcdata(jj, ii, ip, 3);
                double alpe  = bcdata(jj, ii, ip, 4);
                double betae = bcdata(jj, ii, ip, 5);

                alpe  = alpe  / (double)radtodeg;
                betae = betae / (double)radtodeg;
                double te = tte / (1.0 + 0.5 * (double)gm1 * xme * xme);
                double ae = std::sqrt(te);
                double ue = ae * xme * std::cos(alpe) * std::cos(betae);
                double ve = -ae * xme * std::sin(betae);
                double we = ae * xme * std::sin(alpe) * std::cos(betae);

                t(js+j-jsta, 21) = ue;
                t(js+j-jsta, 22) = ve;
                t(js+j-jsta, 23) = we;
                t(js+j-jsta, 24) = ae;
                t(js+j-jsta, 25) = pte;
                t(js+j-jsta, 20) = -sk(j, 1, i, 5);
            }
        }

        int jv = (jend - jsta) * (iend - ista);
        rie1de_ns::rie1de(jvdim, t, jv);

        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista) * (jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    qk0(j, i, l, 1) = t(j-jsta+js, l);
                    qk0(j, i, l, 2) = qk0(j, i, l, 1);
                    bck(j, i, 1)    = 0.0;
                }
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j, i, 1, 1) = vist3d(j, 1, i);
                    vk0(j, i, 1, 2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            double ubar = -(qk0(j,i,2,1)*sk(j,1,i,1) +
                                            qk0(j,i,3,1)*sk(j,1,i,2) +
                                            qk0(j,i,4,1)*sk(j,1,i,3));
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(jj, ii, ip, 5+l);
                            if ((float)ubar < 0.f) {
                                tk0(j, i, l, 1) = t1;
                                tk0(j, i, l, 2) = t1;
                            } else {
                                tk0(j, i, l, 1) = tursav(j, 1, i, l);
                                tk0(j, i, l, 2) = tk0(j, i, l, 1);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==5


    // ******************************************************************
    // k=kdim boundary   inflow for engine   bctype 2003
    // ******************************************************************
    if (nface == 6) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,6) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            int js = (i - ista) * (jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    t(js+j-jsta, l) = q(j, kdim1, i, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    t(js+j-jsta, 5+l) = sk(j, kdim, i, l);
                }
            }
            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;
                double xme   = bcdata(jj, ii, ip, 1);
                double pte   = bcdata(jj, ii, ip, 2);
                double tte   = bcdata(jj, ii, ip, 3);
                double alpe  = bcdata(jj, ii, ip, 4);
                double betae = bcdata(jj, ii, ip, 5);

                alpe  = alpe  / (double)radtodeg;
                betae = betae / (double)radtodeg;
                double te = tte / (1.0 + 0.5 * (double)gm1 * xme * xme);
                double ae = std::sqrt(te);
                double ue = ae * xme * std::cos(alpe) * std::cos(betae);
                double ve = -ae * xme * std::sin(betae);
                double we = ae * xme * std::sin(alpe) * std::cos(betae);

                t(js+j-jsta, 21) = ue;
                t(js+j-jsta, 22) = ve;
                t(js+j-jsta, 23) = we;
                t(js+j-jsta, 24) = ae;
                t(js+j-jsta, 25) = pte;
                t(js+j-jsta, 20) = sk(j, kdim, i, 5);
            }
        }

        int jv = (jend - jsta) * (iend - ista);
        rie1de_ns::rie1de(jvdim, t, jv);

        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista) * (jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    qk0(j, i, l, 3) = t(j-jsta+js, l);
                    qk0(j, i, l, 4) = qk0(j, i, l, 3);
                    bck(j, i, 2)    = 0.0;
                }
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j, i, 1, 3) = vist3d(j, kdim-1, i);
                    vk0(j, i, 1, 4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            double ubar = qk0(j,i,2,3)*sk(j,kdim,i,1) +
                                          qk0(j,i,3,3)*sk(j,kdim,i,2) +
                                          qk0(j,i,4,3)*sk(j,kdim,i,3);
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(jj, ii, ip, 5+l);
                            if ((float)ubar < 0.f) {
                                tk0(j, i, l, 3) = t1;
                                tk0(j, i, l, 4) = t1;
                            } else {
                                tk0(j, i, l, 3) = tursav(j, kdim-1, i, l);
                                tk0(j, i, l, 4) = tk0(j, i, l, 3);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==6


    // ******************************************************************
    // i=1 boundary   inflow for engine   bctype 2003
    // ******************************************************************
    if (nface == 1) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,6) > -1.e10f) itrflg1 = 1;

        for (int k = ksta; k <= kend1; k++) {
            int kk = k - ksta + 1;
            int js = (k - ksta) * (jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    t(js+j-jsta, l) = q(j, k, 1, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    t(js+j-jsta, 5+l) = -si(j, k, 1, l);
                }
            }
            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;
                double xme   = bcdata(jj, kk, ip, 1);
                double pte   = bcdata(jj, kk, ip, 2);
                double tte   = bcdata(jj, kk, ip, 3);
                double alpe  = bcdata(jj, kk, ip, 4);
                double betae = bcdata(jj, kk, ip, 5);

                alpe  = alpe  / (double)radtodeg;
                betae = betae / (double)radtodeg;
                double te = tte / (1.0 + 0.5 * (double)gm1 * xme * xme);
                double ae = std::sqrt(te);
                double ue = ae * xme * std::cos(alpe) * std::cos(betae);
                double ve = -ae * xme * std::sin(betae);
                double we = ae * xme * std::sin(alpe) * std::cos(betae);

                t(js+j-jsta, 21) = ue;
                t(js+j-jsta, 22) = ve;
                t(js+j-jsta, 23) = we;
                t(js+j-jsta, 24) = ae;
                t(js+j-jsta, 25) = pte;
                t(js+j-jsta, 20) = -si(j, k, 1, 5);
            }
        }

        int jv = (jend - jsta) * (kend - ksta);
        rie1de_ns::rie1de(jvdim, t, jv);

        for (int k = ksta; k <= kend1; k++) {
            int js = (k - ksta) * (jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    qi0(j, k, l, 1) = t(j-jsta+js, l);
                    qi0(j, k, l, 2) = qi0(j, k, l, 1);
                    bci(j, k, 1)    = 0.0;
                }
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j, k, 1, 1) = vist3d(j, k, 1);
                    vi0(j, k, 1, 2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            double ubar = -(qi0(j,k,2,1)*si(j,k,1,1) +
                                            qi0(j,k,3,1)*si(j,k,1,2) +
                                            qi0(j,k,4,1)*si(j,k,1,3));
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(jj, kk, ip, 5+l);
                            if ((float)ubar < 0.f) {
                                ti0(j, k, l, 1) = t1;
                                ti0(j, k, l, 2) = t1;
                            } else {
                                ti0(j, k, l, 1) = tursav(j, k, 1, l);
                                ti0(j, k, l, 2) = ti0(j, k, l, 1);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==1


    // ******************************************************************
    // i=idim boundary   inflow for engine   bctype 2003
    // ******************************************************************
    if (nface == 2) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,6) > -1.e10f) itrflg1 = 1;

        for (int k = ksta; k <= kend1; k++) {
            int kk = k - ksta + 1;
            int js = (k - ksta) * (jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    t(js+j-jsta, l) = q(j, k, idim1, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    t(js+j-jsta, 5+l) = si(j, k, idim, l);
                }
            }
            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;
                double xme   = bcdata(jj, kk, ip, 1);
                double pte   = bcdata(jj, kk, ip, 2);
                double tte   = bcdata(jj, kk, ip, 3);
                double alpe  = bcdata(jj, kk, ip, 4);
                double betae = bcdata(jj, kk, ip, 5);

                alpe  = alpe  / (double)radtodeg;
                betae = betae / (double)radtodeg;
                double te = tte / (1.0 + 0.5 * (double)gm1 * xme * xme);
                double ae = std::sqrt(te);
                double ue = ae * xme * std::cos(alpe) * std::cos(betae);
                double ve = -ae * xme * std::sin(betae);
                double we = ae * xme * std::sin(alpe) * std::cos(betae);

                t(js+j-jsta, 21) = ue;
                t(js+j-jsta, 22) = ve;
                t(js+j-jsta, 23) = we;
                t(js+j-jsta, 24) = ae;
                t(js+j-jsta, 25) = pte;
                t(js+j-jsta, 20) = si(j, k, idim, 5);
            }
        }

        int jv = (jend - jsta) * (kend - ksta);
        rie1de_ns::rie1de(jvdim, t, jv);

        for (int k = ksta; k <= kend1; k++) {
            int js = (k - ksta) * (jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    qi0(j, k, l, 3) = t(j-jsta+js, l);
                    qi0(j, k, l, 4) = qi0(j, k, l, 3);
                    bci(j, k, 2)    = 0.0;
                }
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j, k, 1, 3) = vist3d(j, k, idim-1);
                    vi0(j, k, 1, 4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            double ubar = qi0(j,k,2,3)*si(j,k,idim,1) +
                                          qi0(j,k,3,3)*si(j,k,idim,2) +
                                          qi0(j,k,4,3)*si(j,k,idim,3);
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(jj, kk, ip, 5+l);
                            if ((float)ubar < 0.f) {
                                ti0(j, k, l, 3) = t1;
                                ti0(j, k, l, 4) = t1;
                            } else {
                                ti0(j, k, l, 3) = tursav(j, k, idim-1, l);
                                ti0(j, k, l, 4) = ti0(j, k, l, 3);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==2

} // end bc2003()

} // namespace bc2003_ns
