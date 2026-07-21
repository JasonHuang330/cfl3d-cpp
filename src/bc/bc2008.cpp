// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2008.h"
#include "umalloc.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>

namespace bc2008_ns {

// -----------------------------------------------------------------------
// bc() dispatcher — loops over all blocks and calls bc2008 for each
// segment whose bctype is 2008, 2018, 2028, or 2038
// -----------------------------------------------------------------------
void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
        int& nwork, double& cl,
        FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
        int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg,
        FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
        FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg,
        FortranArray1DRef<int> nblg,
        FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0,
        FortranArray1DRef<int> nbck0,
        FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim,
        FortranArray1DRef<int> nbckdim,
        FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo,
        FortranArray4DRef<int> kbcinfo,
        FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej,
        FortranArray3DRef<int> bcfilek,
        FortranArray3DRef<int> lwdat,
        int& myid,
        FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
        FortranArray1DRef<int> kdimg,
        FortranArray1DRef<char[80]> bcfiles,
        int& mxbcfil, int& nummem)
{
    int iuns = 0;

    for (int ibl = 1; ibl <= nbl; ibl++) {
        int jdim = jdimg(ibl);
        int kdim = kdimg(ibl);
        int idim = idimg(ibl);
        int jdim1 = jdim - 1;
        int kdim1 = kdim - 1;
        int idim1 = idim - 1;

        // Build array views from w[] using lw offsets
        FortranArray4DRef<double> q   (&w(lw(1,ibl)), jdim,  kdim,  idim,  5);
        FortranArray4DRef<double> qj0 (&w(lw(2,ibl)), kdim,  idim1, 5,     4);
        FortranArray4DRef<double> qk0 (&w(lw(3,ibl)), jdim,  idim1, 5,     4);
        FortranArray4DRef<double> qi0 (&w(lw(4,ibl)), jdim,  kdim,  5,     4);
        FortranArray4DRef<double> sj  (&w(lw(5,ibl)), jdim,  kdim,  idim1, 5);
        FortranArray4DRef<double> sk  (&w(lw(6,ibl)), jdim,  kdim,  idim1, 5);
        FortranArray4DRef<double> si  (&w(lw(7,ibl)), jdim,  kdim,  idim,  5);
        FortranArray3DRef<double> bcj (&w(lw(8,ibl)), kdim,  idim1, 2);
        FortranArray3DRef<double> bck (&w(lw(9,ibl)), jdim,  idim1, 2);
        FortranArray3DRef<double> bci (&w(lw(10,ibl)),jdim,  kdim,  2);
        FortranArray4DRef<double> xtbj(&w(lw(11,ibl)),kdim,  idim1, 3,     2);
        FortranArray4DRef<double> xtbk(&w(lw(12,ibl)),jdim,  idim1, 3,     2);
        FortranArray4DRef<double> xtbi(&w(lw(13,ibl)),jdim,  kdim,  3,     2);
        FortranArray4DRef<double> atbj(&w(lw(14,ibl)),kdim,  idim1, 3,     2);
        FortranArray4DRef<double> atbk(&w(lw(15,ibl)),jdim,  idim1, 3,     2);
        FortranArray4DRef<double> atbi(&w(lw(16,ibl)),jdim,  kdim,  3,     2);

        // Turbulence / viscosity arrays from lw2
        FortranArray4DRef<double> tursav(&w(lw2(1,ibl)), jdim,  kdim,  idim,  nummem);
        FortranArray4DRef<double> tj0   (&w(lw2(2,ibl)), kdim,  idim1, nummem, 4);
        FortranArray4DRef<double> tk0   (&w(lw2(3,ibl)), jdim,  idim1, nummem, 4);
        FortranArray4DRef<double> ti0   (&w(lw2(4,ibl)), jdim,  kdim,  nummem, 4);
        FortranArray3DRef<double> vist3d(&w(lw2(5,ibl)), jdim,  kdim,  idim);
        FortranArray4DRef<double> vj0   (&w(lw2(6,ibl)), kdim,  idim1, 1,      4);
        FortranArray4DRef<double> vk0   (&w(lw2(7,ibl)), jdim,  idim1, 1,      4);
        FortranArray4DRef<double> vi0   (&w(lw2(8,ibl)), jdim,  kdim,  1,      4);

        // ---- j-face segments ----
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            int bctype = jbcinfo(ibl,iseg,1,1);
            if (bctype != 2008 && bctype != 2018 && bctype != 2028 && bctype != 2038)
                continue;
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
            FortranArray4DRef<double> bcdata(&w(ldata), mdim, ndim, 2, 12);
            int ifilno = bcfilej(ibl,iseg,1);
            char filname[80];
            std::memcpy(filname, bcfiles(ifilno), 80);
            int iflag = bctype;
            bc2008(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                   bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                   ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   mdim, ndim, bcdata, filname, iuns,
                   nou, bou, nbuf, ibufdim, myid, iflag, nummem);
        }

        // ---- k-face segments ----
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            int bctype = kbcinfo(ibl,iseg,1,1);
            if (bctype != 2008 && bctype != 2018 && bctype != 2028 && bctype != 2038)
                continue;
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
            FortranArray4DRef<double> bcdata(&w(ldata), mdim, ndim, 2, 12);
            int ifilno = bcfilek(ibl,iseg,1);
            char filname[80];
            std::memcpy(filname, bcfiles(ifilno), 80);
            int iflag = bctype;
            bc2008(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                   bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                   ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   mdim, ndim, bcdata, filname, iuns,
                   nou, bou, nbuf, ibufdim, myid, iflag, nummem);
        }

        // ---- i-face segments ----
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            int bctype = ibcinfo(ibl,iseg,1,1);
            if (bctype != 2008 && bctype != 2018 && bctype != 2028 && bctype != 2038)
                continue;
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
            FortranArray4DRef<double> bcdata(&w(ldata), mdim, ndim, 2, 12);
            int ifilno = bcfilei(ibl,iseg,1);
            char filname[80];
            std::memcpy(filname, bcfiles(ifilno), 80);
            int iflag = bctype;
            bc2008(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                   bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                   ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   mdim, ndim, bcdata, filname, iuns,
                   nou, bou, nbuf, ibufdim, myid, iflag, nummem);
        }
    } // end ibl loop
}

// -----------------------------------------------------------------------
// bc2008() — specific boundary condition routine
// Sets rho/rho_inf, u/a_inf, v/a_inf, w/a_inf, and up to two turbulence
// quantities; pressure is extrapolated from the interior
// -----------------------------------------------------------------------
void bc2008(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj,  FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci,
            FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk,
            FortranArray4DRef<double> xtbi,
            FortranArray4DRef<double> atbj, FortranArray4DRef<double> atbk,
            FortranArray4DRef<double> atbi,
            int& ista, int& iend, int& jsta, int& jend,
            int& ksta, int& kend, int& nface,
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
            int& nbuf, int& ibufdim,
            int& myid, int& iflag, int& nummem)
{
    // COMMON block aliases
    float& gamma   = cmn_fluid.gamma;
    float& xmach   = cmn_info.xmach;
    float  time    = cmn_unst.time;
    int&   level   = cmn_mgrd.level;
    int&   lglobal = cmn_mgrd.lglobal;
    int*   ivisc   = cmn_reyue.ivisc;   // 0-based
    float* tur10   = cmn_ivals.tur10;   // 0-based
    int&   ialph   = cmn_igrdtyp.ialph;

    // Local variables
    double twopi = 8.0 * std::atan(1.0);

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    int stats_val = 0;
    int jsize = 0;

    // Allocate harvest array if needed (iflag == 2038)
    std::vector<double> harvest_vec;
    if (iflag == 2038) {
        if (nface == 1 || nface == 2)
            jsize = (jend1 - jsta + 1) * (kend1 - ksta + 1);
        if (nface == 3 || nface == 4)
            jsize = (iend1 - ista + 1) * (kend1 - ksta + 1);
        if (nface == 5 || nface == 6)
            jsize = (iend1 - ista + 1) * (jend1 - jsta + 1);
        harvest_vec.assign(jsize, 0.0);
        // call umalloc(jsize,0,'harvest',memuse,stats)
        {
            int memuse = 0;
            char text1[1] = {0};
            umalloc_ns::umalloc(jsize, stats_val, text1, memuse, stats_val);
        }
        // call random_seed; call random_number(harvest)
        // Use C++ random number generation (uniform [0,1))
        for (int idx = 0; idx < jsize; idx++) {
            harvest_vec[idx] = (double)std::rand() / ((double)RAND_MAX + 1.0);
        }
    }
    // Provide 1-based access to harvest via lambda
    auto harvest = [&](int idx) -> double& { return harvest_vec[idx - 1]; };



    // =========================================================================
    //  j=1 boundary   (nface == 3)
    // =========================================================================
    if (nface == 3) {
        int ipp     = 1;
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,5) > -1.e10f) itrflg1 = 1;

        if (ialph == 0) {
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int i = ista; i <= iend1; i++) {
                    int ii = i - ista + 1;
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qj0(k,i,1,ip) = bcdata(kk,ii,ip,1);
                            qj0(k,i,2,ip) = bcdata(kk,ii,ip,2);
                            qj0(k,i,3,ip) = bcdata(kk,ii,ip,3);
                            qj0(k,i,4,ip) = bcdata(kk,ii,ip,4);
                            qj0(k,i,5,ip) = q(1,k,i,5);
                        } else if (iflag == 2018) {
                            qj0(k,i,5,ip) = q(1,k,i,5);
                            qj0(k,i,1,ip) = (double)gamma * qj0(k,i,5,ip) / bcdata(kk,ii,ip,1);
                            qj0(k,i,2,ip) = bcdata(kk,ii,ip,2) / qj0(k,i,1,ip);
                            qj0(k,i,3,ip) = bcdata(kk,ii,ip,3) / qj0(k,i,1,ip);
                            qj0(k,i,4,ip) = bcdata(kk,ii,ip,4) / qj0(k,i,1,ip);
                        } else if (iflag == 2028) {
                            qj0(k,i,1,ip) = q(1,k,i,1);
                            qj0(k,i,2,ip) = bcdata(kk,ii,ip,2) / qj0(k,i,1,ip) *
                                            std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,3,ip) = bcdata(kk,ii,ip,3) / qj0(k,i,1,ip) *
                                            std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,4,ip) = bcdata(kk,ii,ip,4) / qj0(k,i,1,ip) *
                                            std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,5,ip) = q(1,k,i,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(kk,ii,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(kk,ii,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(kk,ii,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(kk,ii,ip,2) + 8.5) * rnd_nbr;
                            }
                            qj0(k,i,1,ip) = bcdata(kk,ii,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,2,ip) = bcdata(kk,ii,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,3,ip) = bcdata(kk,ii,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,4,ip) = bcdata(kk,ii,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,5,ip) = q(1,k,i,5);
                        }
                        bcj(k,i,1) = 0.0;
                    }
                }
            }
        } else {
            // ialph != 0
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int i = ista; i <= iend1; i++) {
                    int ii = i - ista + 1;
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qj0(k,i,1,ip) =  bcdata(kk,ii,ip,1);
                            qj0(k,i,2,ip) =  bcdata(kk,ii,ip,2);
                            qj0(k,i,3,ip) = -bcdata(kk,ii,ip,4);
                            qj0(k,i,4,ip) =  bcdata(kk,ii,ip,3);
                            qj0(k,i,5,ip) =  q(1,k,i,5);
                        } else if (iflag == 2018) {
                            qj0(k,i,5,ip) =  q(1,k,i,5);
                            qj0(k,i,1,ip) =  (double)gamma * qj0(k,i,5,ip) / bcdata(kk,ii,ip,1);
                            qj0(k,i,2,ip) =  bcdata(kk,ii,ip,2) / qj0(k,i,1,ip);
                            qj0(k,i,3,ip) = -bcdata(kk,ii,ip,4) / qj0(k,i,1,ip);
                            qj0(k,i,4,ip) =  bcdata(kk,ii,ip,3) / qj0(k,i,1,ip);
                        } else if (iflag == 2028) {
                            qj0(k,i,1,ip) =  q(1,k,i,1);
                            qj0(k,i,2,ip) =  bcdata(kk,ii,ip,2) / qj0(k,i,1,ip) *
                                             std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,3,ip) = -bcdata(kk,ii,ip,4) / qj0(k,i,1,ip) *
                                             std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,4,ip) =  bcdata(kk,ii,ip,3) / qj0(k,i,1,ip) *
                                             std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,5,ip) =  q(1,k,i,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(kk,ii,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(kk,ii,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(kk,ii,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(kk,ii,ip,2) + 8.5) * rnd_nbr;
                            }
                            qj0(k,i,1,ip) = bcdata(kk,ii,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,2,ip) = bcdata(kk,ii,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,3,ip) =-bcdata(kk,ii,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,4,ip) = bcdata(kk,ii,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,5,ip) = q(1,k,i,5);
                        }
                        bcj(k,i,1) = 0.0;
                    }
                }
            }
        } // end ialph

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,1) = vist3d(1,k,i);
                    vj0(k,i,1,2) = vist3d(1,k,i);
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
                            int ip;
                            ip = 1;
                            double t11 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(kk,ii,ip,4+l);
                            ip = 2;
                            double t12 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(kk,ii,ip,4+l);
                            tj0(k,i,l,1) = t11;
                            tj0(k,i,l,2) = t12;
                        }
                    }
                }
            }
        }
    } // end nface==3



    // =========================================================================
    //  j=jdim boundary   (nface == 4)
    // =========================================================================
    if (nface == 4) {
        int ipp     = 1;
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,5) > -1.e10f) itrflg1 = 1;

        if (ialph == 0) {
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int i = ista; i <= iend1; i++) {
                    int ii = i - ista + 1;
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qj0(k,i,1,ip+2) = bcdata(kk,ii,ip,1);
                            qj0(k,i,2,ip+2) = bcdata(kk,ii,ip,2);
                            qj0(k,i,3,ip+2) = bcdata(kk,ii,ip,3);
                            qj0(k,i,4,ip+2) = bcdata(kk,ii,ip,4);
                            qj0(k,i,5,ip+2) = q(jdim1,k,i,5);
                        } else if (iflag == 2018) {
                            qj0(k,i,5,ip+2) = q(jdim1,k,i,5);
                            qj0(k,i,1,ip+2) = (double)gamma * qj0(k,i,5,ip+2) / bcdata(kk,ii,ip,1);
                            qj0(k,i,2,ip+2) = bcdata(kk,ii,ip,2) / qj0(k,i,1,ip+2);
                            qj0(k,i,3,ip+2) = bcdata(kk,ii,ip,3) / qj0(k,i,1,ip+2);
                            qj0(k,i,4,ip+2) = bcdata(kk,ii,ip,4) / qj0(k,i,1,ip+2);
                        } else if (iflag == 2028) {
                            qj0(k,i,1,ip+2) = q(jdim1,k,i,1);
                            qj0(k,i,2,ip+2) = bcdata(kk,ii,ip,2) / qj0(k,i,1,ip+2) *
                                              std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,3,ip+2) = bcdata(kk,ii,ip,3) / qj0(k,i,1,ip+2) *
                                              std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,4,ip+2) = bcdata(kk,ii,ip,4) / qj0(k,i,1,ip+2) *
                                              std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,5,ip+2) = q(jdim1,k,i,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(kk,ii,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(kk,ii,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(kk,ii,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(kk,ii,ip,2) + 8.5) * rnd_nbr;
                            }
                            qj0(k,i,1,ip+2) = bcdata(kk,ii,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,2,ip+2) = bcdata(kk,ii,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,3,ip+2) = bcdata(kk,ii,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,4,ip+2) = bcdata(kk,ii,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,5,ip+2) = q(jdim1,k,i,5);
                        }
                        bcj(k,i,2) = 0.0;
                    }
                }
            }
        } else {
            // ialph != 0
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int i = ista; i <= iend1; i++) {
                    int ii = i - ista + 1;
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qj0(k,i,1,ip+2) =  bcdata(kk,ii,ip,1);
                            qj0(k,i,2,ip+2) =  bcdata(kk,ii,ip,2);
                            qj0(k,i,3,ip+2) = -bcdata(kk,ii,ip,4);
                            qj0(k,i,4,ip+2) =  bcdata(kk,ii,ip,3);
                            qj0(k,i,5,ip+2) =  q(jdim1,k,i,5);
                        } else if (iflag == 2018) {
                            qj0(k,i,5,ip+2) =  q(jdim1,k,i,5);
                            qj0(k,i,1,ip+2) =  (double)gamma * qj0(k,i,5,ip+2) / bcdata(kk,ii,ip,1);
                            qj0(k,i,2,ip+2) =  bcdata(kk,ii,ip,2) / qj0(k,i,1,ip+2);
                            qj0(k,i,3,ip+2) = -bcdata(kk,ii,ip,4) / qj0(k,i,1,ip+2);
                            qj0(k,i,4,ip+2) =  bcdata(kk,ii,ip,3) / qj0(k,i,1,ip+2);
                        } else if (iflag == 2028) {
                            qj0(k,i,1,ip+2) =  q(jdim1,k,i,1);
                            qj0(k,i,2,ip+2) =  bcdata(kk,ii,ip,2) / qj0(k,i,1,ip+2) *
                                               std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,3,ip+2) = -bcdata(kk,ii,ip,4) / qj0(k,i,1,ip+2) *
                                               std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,4,ip+2) =  bcdata(kk,ii,ip,3) / qj0(k,i,1,ip+2) *
                                               std::cos(twopi * bcdata(kk,ii,ip,1) * (double)time);
                            qj0(k,i,5,ip+2) =  q(jdim1,k,i,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(kk,ii,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(kk,ii,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(kk,ii,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(kk,ii,ip,2) + 8.5) * rnd_nbr;
                            }
                            qj0(k,i,1,ip+2) = bcdata(kk,ii,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,2,ip+2) = bcdata(kk,ii,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,3,ip+2) =-bcdata(kk,ii,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,4,ip+2) = bcdata(kk,ii,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qj0(k,i,5,ip+2) = q(jdim1,k,i,5);
                        }
                        bcj(k,i,2) = 0.0;
                    }
                }
            }
        } // end ialph

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,3) = vist3d(jdim1,k,i);
                    vj0(k,i,1,4) = vist3d(jdim1,k,i);
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
                            int ip;
                            ip = 1;
                            double t13 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(kk,ii,ip,4+l);
                            ip = 2;
                            double t14 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(kk,ii,ip,4+l);
                            tj0(k,i,l,3) = t13;
                            tj0(k,i,l,4) = t14;
                        }
                    }
                }
            }
        }
    } // end nface==4



    // =========================================================================
    //  k=1 boundary   (nface == 5)
    // =========================================================================
    if (nface == 5) {
        int ipp     = 1;
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,5) > -1.e10f) itrflg1 = 1;

        if (ialph == 0) {
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int i = ista; i <= iend1; i++) {
                    int ii = i - ista + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qk0(j,i,1,ip) = bcdata(jj,ii,ip,1);
                            qk0(j,i,2,ip) = bcdata(jj,ii,ip,2);
                            qk0(j,i,3,ip) = bcdata(jj,ii,ip,3);
                            qk0(j,i,4,ip) = bcdata(jj,ii,ip,4);
                            qk0(j,i,5,ip) = q(j,1,i,5);
                        } else if (iflag == 2018) {
                            qk0(j,i,5,ip) = q(j,1,i,5);
                            qk0(j,i,1,ip) = (double)gamma * qk0(j,i,5,ip) / bcdata(jj,ii,ip,1);
                            qk0(j,i,2,ip) = bcdata(jj,ii,ip,2) / qk0(j,i,1,ip);
                            qk0(j,i,3,ip) = bcdata(jj,ii,ip,3) / qk0(j,i,1,ip);
                            qk0(j,i,4,ip) = bcdata(jj,ii,ip,4) / qk0(j,i,1,ip);
                        } else if (iflag == 2028) {
                            qk0(j,i,1,ip) = q(j,1,i,1);
                            qk0(j,i,2,ip) = bcdata(jj,ii,ip,2) / qk0(j,i,1,ip) *
                                            std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,3,ip) = bcdata(jj,ii,ip,3) / qk0(j,i,1,ip) *
                                            std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,4,ip) = bcdata(jj,ii,ip,4) / qk0(j,i,1,ip) *
                                            std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,5,ip) = q(j,1,i,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(jj,ii,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(jj,ii,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(jj,ii,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(jj,ii,ip,2) + 8.5) * rnd_nbr;
                            }
                            qk0(j,i,1,ip) = bcdata(jj,ii,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,2,ip) = bcdata(jj,ii,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,3,ip) = bcdata(jj,ii,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,4,ip) = bcdata(jj,ii,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,5,ip) = q(j,1,i,5);
                        }
                        bck(j,i,1) = 0.0;
                    }
                }
            }
        } else {
            // ialph != 0
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int i = ista; i <= iend1; i++) {
                    int ii = i - ista + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qk0(j,i,1,ip) =  bcdata(jj,ii,ip,1);
                            qk0(j,i,2,ip) =  bcdata(jj,ii,ip,2);
                            qk0(j,i,3,ip) = -bcdata(jj,ii,ip,4);
                            qk0(j,i,4,ip) =  bcdata(jj,ii,ip,3);
                            qk0(j,i,5,ip) =  q(j,1,i,5);
                        } else if (iflag == 2018) {
                            qk0(j,i,5,ip) =  q(j,1,i,5);
                            qk0(j,i,1,ip) =  (double)gamma * qk0(j,i,5,ip) / bcdata(jj,ii,ip,1);
                            qk0(j,i,2,ip) =  bcdata(jj,ii,ip,2) / qk0(j,i,1,ip);
                            qk0(j,i,3,ip) = -bcdata(jj,ii,ip,4) / qk0(j,i,1,ip);
                            qk0(j,i,4,ip) =  bcdata(jj,ii,ip,3) / qk0(j,i,1,ip);
                        } else if (iflag == 2028) {
                            qk0(j,i,1,ip) =  q(j,1,i,1);
                            qk0(j,i,2,ip) =  bcdata(jj,ii,ip,2) / qk0(j,i,1,ip) *
                                             std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,3,ip) = -bcdata(jj,ii,ip,4) / qk0(j,i,1,ip) *
                                             std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,4,ip) =  bcdata(jj,ii,ip,3) / qk0(j,i,1,ip) *
                                             std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,5,ip) =  q(j,1,i,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(jj,ii,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(jj,ii,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(jj,ii,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(jj,ii,ip,2) + 8.5) * rnd_nbr;
                            }
                            qk0(j,i,1,ip) = bcdata(jj,ii,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,2,ip) = bcdata(jj,ii,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,3,ip) =-bcdata(jj,ii,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,4,ip) = bcdata(jj,ii,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,5,ip) = q(j,1,i,5);
                        }
                        bck(j,i,1) = 0.0;
                    }
                }
            }
        } // end ialph

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,1) = vist3d(j,1,i);
                    vk0(j,i,1,2) = vist3d(j,1,i);
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
                            int ip;
                            ip = 1;
                            double t11 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(jj,ii,ip,4+l);
                            ip = 2;
                            double t12 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(jj,ii,ip,4+l);
                            tk0(j,i,l,1) = t11;
                            tk0(j,i,l,2) = t12;
                        }
                    }
                }
            }
        }
    } // end nface==5



    // =========================================================================
    //  k=kdim boundary   (nface == 6)
    // =========================================================================
    if (nface == 6) {
        int ipp     = 1;
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,5) > -1.e10f) itrflg1 = 1;

        if (ialph == 0) {
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int i = ista; i <= iend1; i++) {
                    int ii = i - ista + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qk0(j,i,1,ip+2) = bcdata(jj,ii,ip,1);
                            qk0(j,i,2,ip+2) = bcdata(jj,ii,ip,2);
                            qk0(j,i,3,ip+2) = bcdata(jj,ii,ip,3);
                            qk0(j,i,4,ip+2) = bcdata(jj,ii,ip,4);
                            qk0(j,i,5,ip+2) = q(j,kdim1,i,5);
                        } else if (iflag == 2018) {
                            qk0(j,i,5,ip+2) = q(j,kdim1,i,5);
                            qk0(j,i,1,ip+2) = (double)gamma * qk0(j,i,5,ip+2) / bcdata(jj,ii,ip,1);
                            qk0(j,i,2,ip+2) = bcdata(jj,ii,ip,2) / qk0(j,i,1,ip+2);
                            qk0(j,i,3,ip+2) = bcdata(jj,ii,ip,3) / qk0(j,i,1,ip+2);
                            qk0(j,i,4,ip+2) = bcdata(jj,ii,ip,4) / qk0(j,i,1,ip+2);
                        } else if (iflag == 2028) {
                            qk0(j,i,1,ip+2) = q(j,kdim1,i,1);
                            qk0(j,i,2,ip+2) = bcdata(jj,ii,ip,2) / qk0(j,i,1,ip+2) *
                                              std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,3,ip+2) = bcdata(jj,ii,ip,3) / qk0(j,i,1,ip+2) *
                                              std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,4,ip+2) = bcdata(jj,ii,ip,4) / qk0(j,i,1,ip+2) *
                                              std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,5,ip+2) = q(j,kdim1,i,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(jj,ii,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(jj,ii,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(jj,ii,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(jj,ii,ip,2) + 8.5) * rnd_nbr;
                            }
                            qk0(j,i,1,ip+2) = bcdata(jj,ii,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,2,ip+2) = bcdata(jj,ii,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,3,ip+2) = bcdata(jj,ii,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,4,ip+2) = bcdata(jj,ii,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,5,ip+2) = q(j,kdim1,i,5);
                        }
                        bck(j,i,2) = 0.0;
                    }
                }
            }
        } else {
            // ialph != 0
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int i = ista; i <= iend1; i++) {
                    int ii = i - ista + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qk0(j,i,1,ip+2) =  bcdata(jj,ii,ip,1);
                            qk0(j,i,2,ip+2) =  bcdata(jj,ii,ip,2);
                            qk0(j,i,3,ip+2) = -bcdata(jj,ii,ip,4);
                            qk0(j,i,4,ip+2) =  bcdata(jj,ii,ip,3);
                            qk0(j,i,5,ip+2) =  q(j,kdim1,i,5);
                        } else if (iflag == 2018) {
                            qk0(j,i,5,ip+2) =  q(j,kdim1,i,5);
                            qk0(j,i,1,ip+2) =  (double)gamma * qk0(j,i,5,ip+2) / bcdata(jj,ii,ip,1);
                            qk0(j,i,2,ip+2) =  bcdata(jj,ii,ip,2) / qk0(j,i,1,ip+2);
                            qk0(j,i,3,ip+2) = -bcdata(jj,ii,ip,4) / qk0(j,i,1,ip+2);
                            qk0(j,i,4,ip+2) =  bcdata(jj,ii,ip,3) / qk0(j,i,1,ip+2);
                        } else if (iflag == 2028) {
                            qk0(j,i,1,ip+2) =  q(j,kdim1,i,1);
                            qk0(j,i,2,ip+2) =  bcdata(jj,ii,ip,2) / qk0(j,i,1,ip+2) *
                                               std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,3,ip+2) = -bcdata(jj,ii,ip,4) / qk0(j,i,1,ip+2) *
                                               std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,4,ip+2) =  bcdata(jj,ii,ip,3) / qk0(j,i,1,ip+2) *
                                               std::cos(twopi * bcdata(jj,ii,ip,1) * (double)time);
                            qk0(j,i,5,ip+2) =  q(j,kdim1,i,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(jj,ii,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(jj,ii,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(jj,ii,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(jj,ii,ip,2) + 8.5) * rnd_nbr;
                            }
                            qk0(j,i,1,ip+2) = bcdata(jj,ii,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,2,ip+2) = bcdata(jj,ii,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,3,ip+2) =-bcdata(jj,ii,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,4,ip+2) = bcdata(jj,ii,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qk0(j,i,5,ip+2) = q(j,kdim1,i,5);
                        }
                        bck(j,i,2) = 0.0;
                    }
                }
            }
        } // end ialph

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,3) = vist3d(j,kdim1,i);
                    vk0(j,i,1,4) = vist3d(j,kdim1,i);
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
                            int ip;
                            ip = 1;
                            double t13 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(jj,ii,ip,4+l);
                            ip = 2;
                            double t14 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(jj,ii,ip,4+l);
                            tk0(j,i,l,3) = t13;
                            tk0(j,i,l,4) = t14;
                        }
                    }
                }
            }
        }
    } // end nface==6



    // =========================================================================
    //  i=1 boundary   (nface == 1)
    // =========================================================================
    if (nface == 1) {
        int ipp     = 1;
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,5) > -1.e10f) itrflg1 = 1;

        if (ialph == 0) {
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int k = ksta; k <= kend1; k++) {
                    int kk = k - ksta + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qi0(j,k,1,ip) = bcdata(jj,kk,ip,1);
                            qi0(j,k,2,ip) = bcdata(jj,kk,ip,2);
                            qi0(j,k,3,ip) = bcdata(jj,kk,ip,3);
                            qi0(j,k,4,ip) = bcdata(jj,kk,ip,4);
                            qi0(j,k,5,ip) = q(j,k,1,5);
                        } else if (iflag == 2018) {
                            qi0(j,k,5,ip) = q(j,k,1,5);
                            qi0(j,k,1,ip) = (double)gamma * qi0(j,k,5,ip) / bcdata(jj,kk,ip,1);
                            qi0(j,k,2,ip) = bcdata(jj,kk,ip,2) / qi0(j,k,1,ip);
                            qi0(j,k,3,ip) = bcdata(jj,kk,ip,3) / qi0(j,k,1,ip);
                            qi0(j,k,4,ip) = bcdata(jj,kk,ip,4) / qi0(j,k,1,ip);
                        } else if (iflag == 2028) {
                            qi0(j,k,1,ip) = q(j,k,1,1);
                            qi0(j,k,2,ip) = bcdata(jj,kk,ip,2) / qi0(j,k,1,ip) *
                                            std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,3,ip) = bcdata(jj,kk,ip,3) / qi0(j,k,1,ip) *
                                            std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,4,ip) = bcdata(jj,kk,ip,4) / qi0(j,k,1,ip) *
                                            std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,5,ip) = q(j,k,1,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(jj,kk,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(jj,kk,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(jj,kk,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(jj,kk,ip,2) + 8.5) * rnd_nbr;
                            }
                            qi0(j,k,1,ip) = bcdata(jj,kk,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,2,ip) = bcdata(jj,kk,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,3,ip) = bcdata(jj,kk,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,4,ip) = bcdata(jj,kk,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,5,ip) = q(j,k,1,5);
                        }
                        bci(j,k,1) = 0.0;
                    }
                }
            }
        } else {
            // ialph != 0
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int k = ksta; k <= kend1; k++) {
                    int kk = k - ksta + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qi0(j,k,1,ip) =  bcdata(jj,kk,ip,1);
                            qi0(j,k,2,ip) =  bcdata(jj,kk,ip,2);
                            qi0(j,k,3,ip) = -bcdata(jj,kk,ip,4);
                            qi0(j,k,4,ip) =  bcdata(jj,kk,ip,3);
                            qi0(j,k,5,ip) =  q(j,k,1,5);
                        } else if (iflag == 2018) {
                            qi0(j,k,5,ip) =  q(j,k,1,5);
                            qi0(j,k,1,ip) =  (double)gamma * qi0(j,k,5,ip) / bcdata(jj,kk,ip,1);
                            qi0(j,k,2,ip) =  bcdata(jj,kk,ip,2) / qi0(j,k,1,ip);
                            qi0(j,k,3,ip) = -bcdata(jj,kk,ip,4) / qi0(j,k,1,ip);
                            qi0(j,k,4,ip) =  bcdata(jj,kk,ip,3) / qi0(j,k,1,ip);
                        } else if (iflag == 2028) {
                            qi0(j,k,1,ip) =  q(j,k,1,1);
                            qi0(j,k,2,ip) =  bcdata(jj,kk,ip,2) / qi0(j,k,1,ip) *
                                             std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,3,ip) = -bcdata(jj,kk,ip,4) / qi0(j,k,1,ip) *
                                             std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,4,ip) =  bcdata(jj,kk,ip,3) / qi0(j,k,1,ip) *
                                             std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,5,ip) =  q(j,k,1,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(jj,kk,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(jj,kk,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(jj,kk,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(jj,kk,ip,2) + 8.5) * rnd_nbr;
                            }
                            qi0(j,k,1,ip) = bcdata(jj,kk,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,2,ip) = bcdata(jj,kk,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,3,ip) =-bcdata(jj,kk,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,4,ip) = bcdata(jj,kk,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,5,ip) = q(j,k,1,5);
                        }
                        bci(j,k,1) = 0.0;
                    }
                }
            }
        } // end ialph

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,1) = vist3d(j,k,1);
                    vi0(j,k,1,2) = vist3d(j,k,1);
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
                            int ip;
                            ip = 1;
                            double t11 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(jj,kk,ip,4+l);
                            ip = 2;
                            double t12 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(jj,kk,ip,4+l);
                            ti0(j,k,l,1) = t11;
                            ti0(j,k,l,2) = t12;
                        }
                    }
                }
            }
        }
    } // end nface==1



    // =========================================================================
    //  i=idim boundary   (nface == 2)
    // =========================================================================
    if (nface == 2) {
        int ipp     = 1;
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,5) > -1.e10f) itrflg1 = 1;

        if (ialph == 0) {
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int k = ksta; k <= kend1; k++) {
                    int kk = k - ksta + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qi0(j,k,1,ip+2) = bcdata(jj,kk,ip,1);
                            qi0(j,k,2,ip+2) = bcdata(jj,kk,ip,2);
                            qi0(j,k,3,ip+2) = bcdata(jj,kk,ip,3);
                            qi0(j,k,4,ip+2) = bcdata(jj,kk,ip,4);
                            qi0(j,k,5,ip+2) = q(j,k,idim1,5);
                        } else if (iflag == 2018) {
                            qi0(j,k,5,ip+2) = q(j,k,idim1,5);
                            qi0(j,k,1,ip+2) = (double)gamma * qi0(j,k,5,ip+2) / bcdata(jj,kk,ip,1);
                            qi0(j,k,2,ip+2) = bcdata(jj,kk,ip,2) / qi0(j,k,1,ip+2);
                            qi0(j,k,3,ip+2) = bcdata(jj,kk,ip,3) / qi0(j,k,1,ip+2);
                            qi0(j,k,4,ip+2) = bcdata(jj,kk,ip,4) / qi0(j,k,1,ip+2);
                        } else if (iflag == 2028) {
                            qi0(j,k,1,ip+2) = q(j,k,idim1,1);
                            qi0(j,k,2,ip+2) = bcdata(jj,kk,ip,2) / qi0(j,k,1,ip+2) *
                                              std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,3,ip+2) = bcdata(jj,kk,ip,3) / qi0(j,k,1,ip+2) *
                                              std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,4,ip+2) = bcdata(jj,kk,ip,4) / qi0(j,k,1,ip+2) *
                                              std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,5,ip+2) = q(j,k,idim1,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(jj,kk,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(jj,kk,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(jj,kk,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(jj,kk,ip,2) + 8.5) * rnd_nbr;
                            }
                            qi0(j,k,1,ip+2) = bcdata(jj,kk,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,2,ip+2) = bcdata(jj,kk,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,3,ip+2) = bcdata(jj,kk,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,4,ip+2) = bcdata(jj,kk,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,5,ip+2) = q(j,k,idim1,5);
                        }
                        bci(j,k,2) = 0.0;
                    }
                }
            }
        } else {
            // ialph != 0
            for (int ip = 1; ip <= 2; ip++) {
                int irnd_count = 0;
                for (int k = ksta; k <= kend1; k++) {
                    int kk = k - ksta + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        irnd_count++;
                        if (iflag == 2008) {
                            qi0(j,k,1,ip+2) =  bcdata(jj,kk,ip,1);
                            qi0(j,k,2,ip+2) =  bcdata(jj,kk,ip,2);
                            qi0(j,k,3,ip+2) = -bcdata(jj,kk,ip,4);
                            qi0(j,k,4,ip+2) =  bcdata(jj,kk,ip,3);
                            qi0(j,k,5,ip+2) =  q(j,k,idim1,5);
                        } else if (iflag == 2018) {
                            qi0(j,k,5,ip+2) =  q(j,k,idim1,5);
                            qi0(j,k,1,ip+2) =  (double)gamma * qi0(j,k,5,ip+2) / bcdata(jj,kk,ip,1);
                            qi0(j,k,2,ip+2) =  bcdata(jj,kk,ip,2) / qi0(j,k,1,ip+2);
                            qi0(j,k,3,ip+2) = -bcdata(jj,kk,ip,4) / qi0(j,k,1,ip+2);
                            qi0(j,k,4,ip+2) =  bcdata(jj,kk,ip,3) / qi0(j,k,1,ip+2);
                        } else if (iflag == 2028) {
                            qi0(j,k,1,ip+2) =  q(j,k,idim1,1);
                            qi0(j,k,2,ip+2) =  bcdata(jj,kk,ip,2) / qi0(j,k,1,ip+2) *
                                               std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,3,ip+2) = -bcdata(jj,kk,ip,4) / qi0(j,k,1,ip+2) *
                                               std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,4,ip+2) =  bcdata(jj,kk,ip,3) / qi0(j,k,1,ip+2) *
                                               std::cos(twopi * bcdata(jj,kk,ip,1) * (double)time);
                            qi0(j,k,5,ip+2) =  q(j,k,idim1,5);
                        } else if (iflag == 2038) {
                            double rnd_nbr = harvest(irnd_count) - 0.5;
                            if ((float)bcdata(jj,kk,ip,2) > 0.85f * (float)xmach) rnd_nbr = 0.0;
                            if ((float)bcdata(jj,kk,ip,2) > 0.75f * (float)xmach &&
                                (float)bcdata(jj,kk,ip,2) <= 0.85f * (float)xmach) {
                                rnd_nbr = (-10.0 / (double)xmach * bcdata(jj,kk,ip,2) + 8.5) * rnd_nbr;
                            }
                            qi0(j,k,1,ip+2) = bcdata(jj,kk,ip,1) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,2,ip+2) = bcdata(jj,kk,ip,2) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,3,ip+2) =-bcdata(jj,kk,ip,4) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,4,ip+2) = bcdata(jj,kk,ip,3) * (1.0 + rnd_nbr * 0.10);
                            qi0(j,k,5,ip+2) = q(j,k,idim1,5);
                        }
                        bci(j,k,2) = 0.0;
                    }
                }
            }
        } // end ialph

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,3) = vist3d(j,k,idim1);
                    vi0(j,k,1,4) = vist3d(j,k,idim1);
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
                            int ip;
                            ip = 1;
                            double t13 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(jj,kk,ip,4+l);
                            ip = 2;
                            double t14 = (1 - itrflg1) * (double)tur10[l-1] + itrflg1 * bcdata(jj,kk,ip,4+l);
                            ti0(j,k,l,3) = t13;
                            ti0(j,k,l,4) = t14;
                        }
                    }
                }
            }
        }
    } // end nface==2

    // harvest is automatically freed when harvest_vec goes out of scope
    // (corresponds to Fortran: if (iflag .eq. 2038) deallocate(harvest))

} // end bc2008()

} // namespace bc2008_ns
