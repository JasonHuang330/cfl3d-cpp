// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2009.h"
#include <cmath>
#include <cstring>
#include <algorithm>

namespace bc2009_ns {

// -----------------------------------------------------------------------
// bc dispatcher: loops over all blocks, finds segments with bctype==2009,
// extracts per-block arrays from global work arrays, and calls bc2009.
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

        // Build turbulence array views from wk using lw2 offsets
        FortranArray4DRef<double> tursav(&wk(lw2(1,ibl)), jdim,   kdim,   idim,   nummem);
        FortranArray4DRef<double> tj0   (&wk(lw2(2,ibl)), kdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> tk0   (&wk(lw2(3,ibl)), jdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> ti0   (&wk(lw2(4,ibl)), jdim,   kdim,   nummem, 4);
        FortranArray3DRef<double> vist3d(&wk(lw2(5,ibl)), jdim,   kdim,   idim);
        FortranArray4DRef<double> vj0   (&wk(lw2(6,ibl)), kdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vk0   (&wk(lw2(7,ibl)), jdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vi0   (&wk(lw2(8,ibl)), jdim,   kdim,   1,      4);

        // ---- j-face segments ----
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl,iseg,1,1) != 2009) continue;
            int ista  = jbcinfo(ibl,iseg,1,2);
            int iend  = jbcinfo(ibl,iseg,1,3);
            int jsta  = jbcinfo(ibl,iseg,1,4);
            int jend  = jbcinfo(ibl,iseg,1,5);
            int ksta  = jbcinfo(ibl,iseg,1,6);
            int kend  = jbcinfo(ibl,iseg,1,7);
            int nface = jbcinfo(ibl,iseg,1,8);
            int mdim  = nbcjdim(ibl);
            int ndim  = jbcinfo(ibl,iseg,1,3) - jbcinfo(ibl,iseg,1,2) + 1;
            // bcdata stored in wk at offset lwdat(ibl,iseg,2) for j-faces
            int ldata = lwdat(ibl,iseg,2);
            int jfil  = bcfilej(ibl,iseg,1);
            // mdim/ndim for bcdata: use nbcjdim for mdim
            // The bcdata array is (mdim,ndim,2,12)
            // For j-faces: mdim = kend-ksta+1 (or nbcjdim), ndim = iend-ista+1
            int bmdim = kend - ksta + 1;
            int bndim = iend - ista + 1;
            FortranArray4DRef<double> bcdata(&wk(ldata), bmdim, bndim, 2, 12);
            // filname from bcfiles
            char filname[80];
            std::memcpy(filname, bcfiles(jfil), 80);
            bc2009(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                   bcj, bck, bci, ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   bmdim, bndim, bcdata, filname, iuns,
                   nou, bou, nbuf, ibufdim, myid, nummem);
        }

        // ---- k-face segments ----
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) != 2009) continue;
            int ista  = kbcinfo(ibl,iseg,1,2);
            int iend  = kbcinfo(ibl,iseg,1,3);
            int jsta  = kbcinfo(ibl,iseg,1,4);
            int jend  = kbcinfo(ibl,iseg,1,5);
            int ksta  = kbcinfo(ibl,iseg,1,6);
            int kend  = kbcinfo(ibl,iseg,1,7);
            int nface = kbcinfo(ibl,iseg,1,8);
            int kfil  = bcfilek(ibl,iseg,1);
            int bmdim = jend - jsta + 1;
            int bndim = iend - ista + 1;
            int ldata = lwdat(ibl,iseg,3);
            FortranArray4DRef<double> bcdata(&wk(ldata), bmdim, bndim, 2, 12);
            char filname[80];
            std::memcpy(filname, bcfiles(kfil), 80);
            bc2009(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                   bcj, bck, bci, ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   bmdim, bndim, bcdata, filname, iuns,
                   nou, bou, nbuf, ibufdim, myid, nummem);
        }

        // ---- i-face segments ----
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) != 2009) continue;
            int ista  = ibcinfo(ibl,iseg,1,2);
            int iend  = ibcinfo(ibl,iseg,1,3);
            int jsta  = ibcinfo(ibl,iseg,1,4);
            int jend  = ibcinfo(ibl,iseg,1,5);
            int ksta  = ibcinfo(ibl,iseg,1,6);
            int kend  = ibcinfo(ibl,iseg,1,7);
            int nface = ibcinfo(ibl,iseg,1,8);
            int ifil  = bcfilei(ibl,iseg,1);
            int bmdim = jend - jsta + 1;
            int bndim = kend - ksta + 1;
            int ldata = lwdat(ibl,iseg,1);
            FortranArray4DRef<double> bcdata(&wk(ldata), bmdim, bndim, 2, 12);
            char filname[80];
            std::memcpy(filname, bcfiles(ifil), 80);
            bc2009(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                   bcj, bck, bci, ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   bmdim, bndim, bcdata, filname, iuns,
                   nou, bou, nbuf, ibufdim, myid, nummem);
        }
    }
}

// -----------------------------------------------------------------------
// bc2009: Set inflow boundary conditions (nozzle total BCs)
//         given total pressure ratio, total temperature, and flow angle.
// -----------------------------------------------------------------------
void bc2009(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj,  FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            int& nface,
            FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0,
            int& mdim, int& ndim,
            FortranArray4DRef<double> bcdata,
            char (&filname)[80],
            int& iuns,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid, int& nummem)
{

    // COMMON block aliases
    float& gamma    = cmn_fluid.gamma;
    float& gm1      = cmn_fluid.gm1;
    float& gm1g     = cmn_fluid.gm1g;
    int&   level    = cmn_mgrd.level;
    int&   lglobal  = cmn_mgrd.lglobal;
    float& p0       = cmn_ivals.p0;
    float* tur10    = cmn_ivals.tur10;   // 1-based: tur10[l-1]
    int*   ivisc    = cmn_reyue.ivisc;   // 1-based: ivisc[0..2]
    float& radtodeg = cmn_conversion.radtodeg;

    // Dimension helpers
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    // This bc makes use of only one plane of data
    int ip = 1;

    // Local working variables
    double pte, tte, alpe, betae;
    double pe, plimit, rme2, te, vtote, rhoe;
    double ue, ve, we;
    double ubar, t1;
    int    itrflg1;
    int    ii, kk, jj, js;
    int    l;

    // ******************************************************************
    // j=1 boundary             nozzle total BCs                type 2009
    // ******************************************************************
    if (nface == 3) {

        // check to see if turbulence data is input (itrflg1 = 1) or
        // if freestream values are to be used (itrflg1 = 0)
        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;
            js = (i - ista) * (kend - ksta) + 1;

            for (int k = ksta; k <= kend1; k++) {
                kk = k - ksta + 1;

                pte   = bcdata(kk,ii,ip,1);
                tte   = bcdata(kk,ii,ip,2);
                alpe  = bcdata(kk,ii,ip,3);
                betae = bcdata(kk,ii,ip,4);

                pte   = pte * (double)p0;
                alpe  = alpe / (double)radtodeg;
                betae = betae / (double)radtodeg;

                pe = q(1,k,i,5);
                plimit = 0.999 * pte;
                if ((float)pe >= (float)plimit) pe = plimit;
                rme2  = 2.0 / (double)gm1 * (std::pow(pte/pe, (double)gm1g) - 1.0);
                te    = tte / (1.0 + (double)gm1 * 0.5 * rme2);
                vtote = std::sqrt(rme2 * te);
                rhoe  = pe * (double)gamma / te;

                ue = vtote * std::cos(alpe) * std::cos(betae);
                ve = -vtote * std::sin(betae);
                we = vtote * std::sin(alpe) * std::cos(betae);

                qj0(k,i,1,1) = rhoe;
                qj0(k,i,1,2) = qj0(k,i,1,1);
                qj0(k,i,2,1) = ue;
                qj0(k,i,2,2) = qj0(k,i,2,1);
                qj0(k,i,3,1) = ve;
                qj0(k,i,3,2) = qj0(k,i,3,1);
                qj0(k,i,4,1) = we;
                qj0(k,i,4,2) = qj0(k,i,4,1);
                qj0(k,i,5,1) = pe;
                qj0(k,i,5,2) = qj0(k,i,5,1);

                bcj(k,i,1) = 0.0;
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,1) = vist3d(1,k,i);
                    vj0(k,i,1,2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            kk = k - ksta + 1;
                            ubar = -(qj0(k,i,2,1)*sj(1,k,i,1) +
                                     qj0(k,i,3,1)*sj(1,k,i,2) +
                                     qj0(k,i,4,1)*sj(1,k,i,3));
                            t1 = (double)(1 - itrflg1) * (double)tur10[l-1]
                               + (double)itrflg1 * bcdata(kk,ii,ip,4+l);
                            if ((float)ubar < 0.f) {
                                tj0(k,i,l,1) = t1;
                                tj0(k,i,l,2) = t1;
                            } else {
                                tj0(k,i,l,1) = tursav(1,k,i,l);
                                tj0(k,i,l,2) = tj0(k,i,l,1);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==3


    // ******************************************************************
    // j=jdim boundary          nozzle total BCs                type 2009
    // ******************************************************************
    if (nface == 4) {

        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;
            js = (i - ista) * (kend - ksta) + 1;

            for (int k = ksta; k <= kend1; k++) {
                kk = k - ksta + 1;

                pte   = bcdata(kk,ii,ip,1);
                tte   = bcdata(kk,ii,ip,2);
                alpe  = bcdata(kk,ii,ip,3);
                betae = bcdata(kk,ii,ip,4);

                pte   = pte * (double)p0;
                alpe  = alpe / (double)radtodeg;
                betae = betae / (double)radtodeg;

                pe = q(jdim1,k,i,5);
                plimit = 0.999 * pte;
                if ((float)pe >= (float)plimit) pe = plimit;
                rme2  = 2.0 / (double)gm1 * (std::pow(pte/pe, (double)gm1g) - 1.0);
                te    = tte / (1.0 + (double)gm1 * 0.5 * rme2);
                vtote = std::sqrt(rme2 * te);
                rhoe  = pe * (double)gamma / te;

                ue = vtote * std::cos(alpe) * std::cos(betae);
                ve = -vtote * std::sin(betae);
                we = vtote * std::sin(alpe) * std::cos(betae);

                qj0(k,i,1,3) = rhoe;
                qj0(k,i,1,4) = qj0(k,i,1,3);
                qj0(k,i,2,3) = ue;
                qj0(k,i,2,4) = qj0(k,i,2,3);
                qj0(k,i,3,3) = ve;
                qj0(k,i,3,4) = qj0(k,i,3,3);
                qj0(k,i,4,3) = we;
                qj0(k,i,4,4) = qj0(k,i,4,3);
                qj0(k,i,5,3) = pe;
                qj0(k,i,5,4) = qj0(k,i,5,3);

                bcj(k,i,2) = 0.0;
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,3) = vist3d(jdim-1,k,i);
                    vj0(k,i,1,4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            kk = k - ksta + 1;
                            ubar = qj0(k,i,2,3)*sj(jdim,k,i,1) +
                                   qj0(k,i,3,3)*sj(jdim,k,i,2) +
                                   qj0(k,i,4,3)*sj(jdim,k,i,3);
                            t1 = (double)(1 - itrflg1) * (double)tur10[l-1]
                               + (double)itrflg1 * bcdata(kk,ii,ip,4+l);
                            if ((float)ubar < 0.f) {
                                tj0(k,i,l,3) = t1;
                                tj0(k,i,l,4) = t1;
                            } else {
                                tj0(k,i,l,3) = tursav(jdim-1,k,i,l);
                                tj0(k,i,l,4) = tj0(k,i,l,3);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==4


    // ******************************************************************
    // k=1 boundary             nozzle total BCs                type 2009
    // ******************************************************************
    if (nface == 5) {

        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;
            js = (i - ista) * (jend - jsta) + 1;

            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;

                pte   = bcdata(jj,ii,ip,1);
                tte   = bcdata(jj,ii,ip,2);
                alpe  = bcdata(jj,ii,ip,3);
                betae = bcdata(jj,ii,ip,4);

                pte   = pte * (double)p0;
                alpe  = alpe / (double)radtodeg;
                betae = betae / (double)radtodeg;

                pe = q(j,1,i,5);
                plimit = 0.999 * pte;
                if ((float)pe >= (float)plimit) pe = plimit;
                rme2  = 2.0 / (double)gm1 * (std::pow(pte/pe, (double)gm1g) - 1.0);
                te    = tte / (1.0 + (double)gm1 * 0.5 * rme2);
                vtote = std::sqrt(rme2 * te);
                rhoe  = pe * (double)gamma / te;

                ue = vtote * std::cos(alpe) * std::cos(betae);
                ve = -vtote * std::sin(betae);
                we = vtote * std::sin(alpe) * std::cos(betae);

                qk0(j,i,1,1) = rhoe;
                qk0(j,i,1,2) = qk0(j,i,1,1);
                qk0(j,i,2,1) = ue;
                qk0(j,i,2,2) = qk0(j,i,2,1);
                qk0(j,i,3,1) = ve;
                qk0(j,i,3,2) = qk0(j,i,3,1);
                qk0(j,i,4,1) = we;
                qk0(j,i,4,2) = qk0(j,i,4,1);
                qk0(j,i,5,1) = pe;
                qk0(j,i,5,2) = qk0(j,i,5,1);

                bck(j,i,1) = 0.0;
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,1) = vist3d(j,1,i);
                    vk0(j,i,1,2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            jj = j - jsta + 1;
                            ubar = -(qk0(j,i,2,1)*sk(j,1,i,1) +
                                     qk0(j,i,3,1)*sk(j,1,i,2) +
                                     qk0(j,i,4,1)*sk(j,1,i,3));
                            t1 = (double)(1 - itrflg1) * (double)tur10[l-1]
                               + (double)itrflg1 * bcdata(jj,ii,ip,4+l);
                            if ((float)ubar < 0.f) {
                                tk0(j,i,l,1) = t1;
                                tk0(j,i,l,2) = t1;
                            } else {
                                tk0(j,i,l,1) = tursav(j,1,i,l);
                                tk0(j,i,l,2) = tk0(j,i,l,1);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==5


    // ******************************************************************
    // k=kdim boundary          nozzle total BCs                type 2009
    // ******************************************************************
    if (nface == 6) {

        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;
            js = (i - ista) * (jend - jsta) + 1;

            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;

                pte   = bcdata(jj,ii,ip,1);
                tte   = bcdata(jj,ii,ip,2);
                alpe  = bcdata(jj,ii,ip,3);
                betae = bcdata(jj,ii,ip,4);

                pte   = pte * (double)p0;
                alpe  = alpe / (double)radtodeg;
                betae = betae / (double)radtodeg;

                pe = q(j,kdim1,i,5);
                plimit = 0.999 * pte;
                if ((float)pe >= (float)plimit) pe = plimit;
                rme2  = 2.0 / (double)gm1 * (std::pow(pte/pe, (double)gm1g) - 1.0);
                te    = tte / (1.0 + (double)gm1 * 0.5 * rme2);
                vtote = std::sqrt(rme2 * te);
                rhoe  = pe * (double)gamma / te;

                ue = vtote * std::cos(alpe) * std::cos(betae);
                ve = -vtote * std::sin(betae);
                we = vtote * std::sin(alpe) * std::cos(betae);

                qk0(j,i,1,3) = rhoe;
                qk0(j,i,1,4) = qk0(j,i,1,3);
                qk0(j,i,2,3) = ue;
                qk0(j,i,2,4) = qk0(j,i,2,3);
                qk0(j,i,3,3) = ve;
                qk0(j,i,3,4) = qk0(j,i,3,3);
                qk0(j,i,4,3) = we;
                qk0(j,i,4,4) = qk0(j,i,4,3);
                qk0(j,i,5,3) = pe;
                qk0(j,i,5,4) = qk0(j,i,5,3);

                bck(j,i,2) = 0.0;
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,3) = vist3d(j,kdim-1,i);
                    vk0(j,i,1,4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            jj = j - jsta + 1;
                            ubar = qk0(j,i,2,3)*sk(j,kdim,i,1) +
                                   qk0(j,i,3,3)*sk(j,kdim,i,2) +
                                   qk0(j,i,4,3)*sk(j,kdim,i,3);
                            t1 = (double)(1 - itrflg1) * (double)tur10[l-1]
                               + (double)itrflg1 * bcdata(jj,ii,ip,4+l);
                            if ((float)ubar < 0.f) {
                                tk0(j,i,l,3) = t1;
                                tk0(j,i,l,4) = t1;
                            } else {
                                tk0(j,i,l,3) = tursav(j,kdim-1,i,l);
                                tk0(j,i,l,4) = tk0(j,i,l,3);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==6


    // ******************************************************************
    // i=1 boundary             nozzle total BCs                type 2009
    // ******************************************************************
    if (nface == 1) {

        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        for (int k = ksta; k <= kend1; k++) {
            kk = k - ksta + 1;
            js = (k - ksta) * (jend - jsta) + 1;

            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;

                pte   = bcdata(jj,kk,ip,1);
                tte   = bcdata(jj,kk,ip,2);
                alpe  = bcdata(jj,kk,ip,3);
                betae = bcdata(jj,kk,ip,4);

                pte   = pte * (double)p0;
                alpe  = alpe / (double)radtodeg;
                betae = betae / (double)radtodeg;

                pe = q(j,k,1,5);
                plimit = 0.999 * pte;
                if ((float)pe >= (float)plimit) pe = plimit;
                rme2  = 2.0 / (double)gm1 * (std::pow(pte/pe, (double)gm1g) - 1.0);
                te    = tte / (1.0 + (double)gm1 * 0.5 * rme2);
                vtote = std::sqrt(rme2 * te);
                rhoe  = pe * (double)gamma / te;

                ue = vtote * std::cos(alpe) * std::cos(betae);
                ve = -vtote * std::sin(betae);
                we = vtote * std::sin(alpe) * std::cos(betae);

                qi0(j,k,1,1) = rhoe;
                qi0(j,k,1,2) = qi0(j,k,1,1);
                qi0(j,k,2,1) = ue;
                qi0(j,k,2,2) = qi0(j,k,2,1);
                qi0(j,k,3,1) = ve;
                qi0(j,k,3,2) = qi0(j,k,3,1);
                qi0(j,k,4,1) = we;
                qi0(j,k,4,2) = qi0(j,k,4,1);
                qi0(j,k,5,1) = pe;
                qi0(j,k,5,2) = qi0(j,k,5,1);

                bci(j,k,1) = 0.0;
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,1) = vist3d(j,k,1);
                    vi0(j,k,1,2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            jj = j - jsta + 1;
                            ubar = -(qi0(j,k,2,1)*si(j,k,1,1) +
                                     qi0(j,k,3,1)*si(j,k,1,2) +
                                     qi0(j,k,4,1)*si(j,k,1,3));
                            t1 = (double)(1 - itrflg1) * (double)tur10[l-1]
                               + (double)itrflg1 * bcdata(jj,kk,ip,4+l);
                            if ((float)ubar < 0.f) {
                                ti0(j,k,l,1) = t1;
                                ti0(j,k,l,2) = t1;
                            } else {
                                ti0(j,k,l,1) = tursav(j,k,1,l);
                                ti0(j,k,l,2) = ti0(j,k,l,1);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==1


    // ******************************************************************
    // i=idim boundary          nozzle total BCs                type 2009
    // ******************************************************************
    if (nface == 2) {

        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        for (int k = ksta; k <= kend1; k++) {
            kk = k - ksta + 1;
            js = (k - ksta) * (jend - jsta) + 1;

            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;

                pte   = bcdata(jj,kk,ip,1);
                tte   = bcdata(jj,kk,ip,2);
                alpe  = bcdata(jj,kk,ip,3);
                betae = bcdata(jj,kk,ip,4);

                pte   = pte * (double)p0;
                alpe  = alpe / (double)radtodeg;
                betae = betae / (double)radtodeg;

                pe = q(j,k,idim1,5);
                plimit = 0.999 * pte;
                if ((float)pe >= (float)plimit) pe = plimit;
                rme2  = 2.0 / (double)gm1 * (std::pow(pte/pe, (double)gm1g) - 1.0);
                te    = tte / (1.0 + (double)gm1 * 0.5 * rme2);
                vtote = std::sqrt(rme2 * te);
                rhoe  = pe * (double)gamma / te;

                ue = vtote * std::cos(alpe) * std::cos(betae);
                ve = -vtote * std::sin(betae);
                we = vtote * std::sin(alpe) * std::cos(betae);

                qi0(j,k,1,3) = rhoe;
                qi0(j,k,1,4) = qi0(j,k,1,3);
                qi0(j,k,2,3) = ue;
                qi0(j,k,2,4) = qi0(j,k,2,3);
                qi0(j,k,3,3) = ve;
                qi0(j,k,3,4) = qi0(j,k,3,3);
                qi0(j,k,4,3) = we;
                qi0(j,k,4,4) = qi0(j,k,4,3);
                qi0(j,k,5,3) = pe;
                qi0(j,k,5,4) = qi0(j,k,5,3);

                bci(j,k,2) = 0.0;
            }
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,3) = vist3d(j,k,idim-1);
                    vi0(j,k,1,4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            jj = j - jsta + 1;
                            ubar = qi0(j,k,2,3)*si(j,k,idim,1) +
                                   qi0(j,k,3,3)*si(j,k,idim,2) +
                                   qi0(j,k,4,3)*si(j,k,idim,3);
                            t1 = (double)(1 - itrflg1) * (double)tur10[l-1]
                               + (double)itrflg1 * bcdata(jj,kk,ip,4+l);
                            if ((float)ubar < 0.f) {
                                ti0(j,k,l,3) = t1;
                                ti0(j,k,l,4) = t1;
                            } else {
                                ti0(j,k,l,3) = tursav(j,k,idim-1,l);
                                ti0(j,k,l,4) = ti0(j,k,l,3);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==2

} // bc2009

} // namespace bc2009_ns
