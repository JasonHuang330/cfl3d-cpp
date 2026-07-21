// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc1008.h"
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <cmath>
#include <cstdint>

namespace bc1008_ns {


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
        FortranArray4DRef<double> xtbj (&w(lw(11,ibl)), kdim,   idim-1, 3,      2);
        FortranArray4DRef<double> xtbk (&w(lw(12,ibl)), jdim,   idim-1, 3,      2);
        FortranArray4DRef<double> xtbi (&w(lw(13,ibl)), jdim,   kdim,   3,      2);
        FortranArray4DRef<double> atbj (&w(lw(14,ibl)), kdim,   idim-1, 3,      2);
        FortranArray4DRef<double> atbk (&w(lw(15,ibl)), jdim,   idim-1, 3,      2);
        FortranArray4DRef<double> atbi (&w(lw(16,ibl)), jdim,   kdim,   3,      2);

        // Build array views from wk using lw2 offsets
        FortranArray4DRef<double> tursav(&wk(lw2(1,ibl)), jdim,   kdim,   idim,   nummem);
        FortranArray4DRef<double> tj0   (&wk(lw2(2,ibl)), kdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> tk0   (&wk(lw2(3,ibl)), jdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> ti0   (&wk(lw2(4,ibl)), jdim,   kdim,   nummem, 4);
        FortranArray3DRef<double> vist3d(&wk(lw2(5,ibl)), jdim,   kdim,   idim);
        FortranArray4DRef<double> vj0   (&wk(lw2(6,ibl)), kdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vk0   (&wk(lw2(7,ibl)), jdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vi0   (&wk(lw2(8,ibl)), jdim,   kdim,   1,      4);

        // Loop over j-face boundary segments
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl,iseg,1,1) == 1008) {
                int ista  = jbcinfo(ibl,iseg,1,2);
                int iend  = jbcinfo(ibl,iseg,1,3);
                int jsta  = jbcinfo(ibl,iseg,1,4);
                int jend  = jbcinfo(ibl,iseg,1,5);
                int ksta  = jbcinfo(ibl,iseg,1,6);
                int kend  = jbcinfo(ibl,iseg,1,7);
                int nface = jbcinfo(ibl,iseg,1,8);
                bc1008(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem);
            }
        }

        // Loop over k-face boundary segments
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) == 1008) {
                int ista  = kbcinfo(ibl,iseg,1,2);
                int iend  = kbcinfo(ibl,iseg,1,3);
                int jsta  = kbcinfo(ibl,iseg,1,4);
                int jend  = kbcinfo(ibl,iseg,1,5);
                int ksta  = kbcinfo(ibl,iseg,1,6);
                int kend  = kbcinfo(ibl,iseg,1,7);
                int nface = kbcinfo(ibl,iseg,1,8);
                bc1008(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem);
            }
        }

        // Loop over i-face boundary segments
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) == 1008) {
                int ista  = ibcinfo(ibl,iseg,1,2);
                int iend  = ibcinfo(ibl,iseg,1,3);
                int jsta  = ibcinfo(ibl,iseg,1,4);
                int jend  = ibcinfo(ibl,iseg,1,5);
                int ksta  = ibcinfo(ibl,iseg,1,6);
                int kend  = ibcinfo(ibl,iseg,1,7);
                int nface = ibcinfo(ibl,iseg,1,8);
                bc1008(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem);
            }
        }
    }
}


void bc1008(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0,
            FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj,
            FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj,
            FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci,
            FortranArray4DRef<double> xtbj,
            FortranArray4DRef<double> xtbk,
            FortranArray4DRef<double> xtbi,
            FortranArray4DRef<double> atbj,
            FortranArray4DRef<double> atbk,
            FortranArray4DRef<double> atbi,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            int& nface,
            FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0,
            FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0,
            FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0,
            int& iuns,
            FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& nummem)
{
    // COMMON block aliases
    float& gamma  = cmn_fluid.gamma;
    float& gm1g   = cmn_fluid.gm1g;
    float& p0     = cmn_ivals.p0;
    float& rho0   = cmn_ivals.rho0;
    float& h0     = cmn_ivals.h0;
    float* tur10  = cmn_ivals.tur10;   // 0-based: tur10[l-1]
    int32_t& level   = cmn_mgrd.level;
    int32_t& lglobal = cmn_mgrd.lglobal;
    int32_t* ivisc   = cmn_reyue.ivisc; // 0-based: ivisc[0..2]

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    double ent0, gami;

    // **************************************************************************
    //      j=1 boundary            tunnel inflow                        bctype 1008
    // **************************************************************************

    if (nface == 3) {

        ent0 = (double)p0 / std::pow((double)rho0, (double)gamma);
        gami = 1.0 / (double)gamma;

        for (int i = ista; i <= iend1; i++) {
        for (int k = ksta; k <= kend1; k++) {

            qj0(k,i,5,1) = q(1,k,i,5);
            qj0(k,i,1,1) = std::pow(qj0(k,i,5,1)/ent0, gami);
            qj0(k,i,2,1) = std::sqrt(2.0*((double)h0 - qj0(k,i,5,1)/qj0(k,i,1,1)/(double)gm1g));
            qj0(k,i,3,1) = 0.0;
            qj0(k,i,4,1) = 0.0;

            qj0(k,i,1,2) = qj0(k,i,1,1);
            qj0(k,i,2,2) = qj0(k,i,2,1);
            qj0(k,i,3,2) = qj0(k,i,3,1);
            qj0(k,i,4,2) = qj0(k,i,4,1);
            qj0(k,i,5,2) = qj0(k,i,5,1);

            bcj(k,i,1) = 0.0;

        }}

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int k = ksta; k <= kend1; k++) {
                vj0(k,i,1,1) = vist3d(1,k,i);
                vj0(k,i,1,2) = 0.0;
            }}
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    tj0(k,i,l,1) = (double)tur10[l-1];
                    tj0(k,i,l,2) = (double)tur10[l-1];
                }}}
            }
        }

    }

    // **************************************************************************
    //      j=jdim boundary         tunnel inflow                        bctype 1008
    // **************************************************************************

    if (nface == 4) {

        ent0 = (double)p0 / std::pow((double)rho0, (double)gamma);
        gami = 1.0 / (double)gamma;

        for (int i = ista; i <= iend1; i++) {
        for (int k = ksta; k <= kend1; k++) {

            qj0(k,i,5,3) = q(jdim1,k,i,5);
            qj0(k,i,1,3) = std::pow(qj0(k,i,5,3)/ent0, gami);
            qj0(k,i,2,3) = std::sqrt(2.0*((double)h0 - qj0(k,i,5,3)/qj0(k,i,1,3)/(double)gm1g));
            qj0(k,i,3,3) = 0.0;
            qj0(k,i,4,3) = 0.0;

            qj0(k,i,1,4) = qj0(k,i,1,3);
            qj0(k,i,2,4) = qj0(k,i,2,3);
            qj0(k,i,3,4) = qj0(k,i,3,3);
            qj0(k,i,4,4) = qj0(k,i,4,3);
            qj0(k,i,5,4) = qj0(k,i,5,3);

            bcj(k,i,2) = 0.0;

        }}

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int k = ksta; k <= kend1; k++) {
                vj0(k,i,1,3) = vist3d(jdim-1,k,i);
                vj0(k,i,1,4) = 0.0;
            }}
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    tj0(k,i,l,3) = (double)tur10[l-1];
                    tj0(k,i,l,4) = (double)tur10[l-1];
                }}}
            }
        }

    }


    // **************************************************************************
    //      k=1 boundary            tunnel inflow                        bctype 1008
    // **************************************************************************

    if (nface == 5) {

        ent0 = (double)p0 / std::pow((double)rho0, (double)gamma);
        gami = 1.0 / (double)gamma;

        for (int i = ista; i <= iend1; i++) {
        for (int j = jsta; j <= jend1; j++) {

            qk0(j,i,5,1) = q(j,1,i,5);
            qk0(j,i,1,1) = std::pow(qk0(j,i,5,1)/ent0, gami);
            qk0(j,i,2,1) = std::sqrt(2.0*((double)h0 - qk0(j,i,5,1)/qk0(j,i,1,1)/(double)gm1g));
            qk0(j,i,3,1) = 0.0;
            qk0(j,i,4,1) = 0.0;

            qk0(j,i,1,2) = qk0(j,i,1,1);
            qk0(j,i,2,2) = qk0(j,i,2,1);
            qk0(j,i,3,2) = qk0(j,i,3,1);
            qk0(j,i,4,2) = qk0(j,i,4,1);
            qk0(j,i,5,2) = qk0(j,i,5,1);

            bck(j,i,1) = 0.0;

        }}

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int j = jsta; j <= jend1; j++) {
                vk0(j,i,1,1) = vist3d(j,1,i);
                vk0(j,i,1,2) = 0.0;
            }}
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    tk0(j,i,l,1) = (double)tur10[l-1];
                    tk0(j,i,l,2) = (double)tur10[l-1];
                }}}
            }
        }

    }

    // **************************************************************************
    //      k=kdim boundary         tunnel inflow                        bctype 1008
    // **************************************************************************

    if (nface == 6) {

        ent0 = (double)p0 / std::pow((double)rho0, (double)gamma);
        gami = 1.0 / (double)gamma;

        for (int i = ista; i <= iend1; i++) {
        for (int j = jsta; j <= jend1; j++) {

            qk0(j,i,5,3) = q(j,kdim1,i,5);
            qk0(j,i,1,3) = std::pow(qk0(j,i,5,3)/ent0, gami);
            qk0(j,i,2,3) = std::sqrt(2.0*((double)h0 - qk0(j,i,5,3)/qk0(j,i,1,3)/(double)gm1g));
            qk0(j,i,3,3) = 0.0;
            qk0(j,i,4,3) = 0.0;

            qk0(j,i,1,4) = qk0(j,i,1,3);
            qk0(j,i,2,4) = qk0(j,i,2,3);
            qk0(j,i,3,4) = qk0(j,i,3,3);
            qk0(j,i,4,4) = qk0(j,i,4,3);
            qk0(j,i,5,4) = qk0(j,i,5,3);

            bck(j,i,2) = 0.0;

        }}

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int j = jsta; j <= jend1; j++) {
                vk0(j,i,1,3) = vist3d(j,kdim-1,i);
                vk0(j,i,1,4) = 0.0;
            }}
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    tk0(j,i,l,3) = (double)tur10[l-1];
                    tk0(j,i,l,4) = (double)tur10[l-1];
                }}}
            }
        }

    }

    // **************************************************************************
    //      i=1 boundary            tunnel inflow                        bctype 1008
    // **************************************************************************

    if (nface == 1) {

        ent0 = (double)p0 / std::pow((double)rho0, (double)gamma);
        gami = 1.0 / (double)gamma;

        for (int j = jsta; j <= jend1; j++) {
        for (int k = ksta; k <= kend1; k++) {

            qi0(j,k,5,1) = q(j,k,1,5);
            qi0(j,k,1,1) = std::pow(qi0(j,k,5,1)/ent0, gami);
            qi0(j,k,2,1) = std::sqrt(2.0*((double)h0 - qi0(j,k,5,1)/qi0(j,k,1,1)/(double)gm1g));
            qi0(j,k,3,1) = 0.0;
            qi0(j,k,4,1) = 0.0;

            qi0(j,k,1,2) = qi0(j,k,1,1);
            qi0(j,k,2,2) = qi0(j,k,2,1);
            qi0(j,k,3,2) = qi0(j,k,3,1);
            qi0(j,k,4,2) = qi0(j,k,4,1);
            qi0(j,k,5,2) = qi0(j,k,5,1);

            bci(j,k,1) = 0.0;

        }}

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
            for (int j = jsta; j <= jend1; j++) {
                vi0(j,k,1,1) = vist3d(j,k,1);
                vi0(j,k,1,2) = 0.0;
            }}
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    ti0(j,k,l,1) = (double)tur10[l-1];
                    ti0(j,k,l,2) = (double)tur10[l-1];
                }}}
            }
        }

    }

    // **************************************************************************
    //      i=idim boundary          tunnel inflow                       bctype 1008
    // **************************************************************************

    if (nface == 2) {

        ent0 = (double)p0 / std::pow((double)rho0, (double)gamma);
        gami = 1.0 / (double)gamma;

        for (int j = jsta; j <= jend1; j++) {
        for (int k = ksta; k <= kend1; k++) {

            qi0(j,k,5,3) = q(j,k,idim1,5);
            qi0(j,k,1,3) = std::pow(qi0(j,k,5,3)/ent0, gami);
            qi0(j,k,2,3) = std::sqrt(2.0*((double)h0 - qi0(j,k,5,3)/qi0(j,k,1,3)/(double)gm1g));
            qi0(j,k,3,3) = 0.0;
            qi0(j,k,4,3) = 0.0;

            qi0(j,k,1,4) = qi0(j,k,1,1);
            qi0(j,k,2,4) = qi0(j,k,2,1);
            qi0(j,k,3,4) = qi0(j,k,3,1);
            qi0(j,k,4,4) = qi0(j,k,4,1);
            qi0(j,k,5,4) = qi0(j,k,5,1);

            bci(j,k,2) = 0.0;

        }}

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
            for (int j = jsta; j <= jend1; j++) {
                vi0(j,k,1,3) = vist3d(j,k,idim-1);
                vi0(j,k,1,4) = 0.0;
            }}
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    ti0(j,k,l,3) = (double)tur10[l-1];
                    ti0(j,k,l,4) = (double)tur10[l-1];
                }}}
            }
        }

    }

    return;
}

} // namespace bc1008_ns
