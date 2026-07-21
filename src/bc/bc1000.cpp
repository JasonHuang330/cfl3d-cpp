// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc1000.h"

namespace bc1000_ns {

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

        // Extract array views from work array w using lw offsets
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

        // Extract array views from work array wk using lw2 offsets
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
            if (jbcinfo(ibl,iseg,1,1) == 1000) {
                int ista  = jbcinfo(ibl,iseg,1,2);
                int iend  = jbcinfo(ibl,iseg,1,3);
                int jsta  = jbcinfo(ibl,iseg,1,4);
                int jend  = jbcinfo(ibl,iseg,1,5);
                int ksta  = jbcinfo(ibl,iseg,1,6);
                int kend  = jbcinfo(ibl,iseg,1,7);
                int nface = jbcinfo(ibl,iseg,1,8);
                bc1000(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem);
            }
        }

        // Loop over k-face boundary segments
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) == 1000) {
                int ista  = kbcinfo(ibl,iseg,1,2);
                int iend  = kbcinfo(ibl,iseg,1,3);
                int jsta  = kbcinfo(ibl,iseg,1,4);
                int jend  = kbcinfo(ibl,iseg,1,5);
                int ksta  = kbcinfo(ibl,iseg,1,6);
                int kend  = kbcinfo(ibl,iseg,1,7);
                int nface = kbcinfo(ibl,iseg,1,8);
                bc1000(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem);
            }
        }

        // Loop over i-face boundary segments
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) == 1000) {
                int ista  = ibcinfo(ibl,iseg,1,2);
                int iend  = ibcinfo(ibl,iseg,1,3);
                int jsta  = ibcinfo(ibl,iseg,1,4);
                int jend  = ibcinfo(ibl,iseg,1,5);
                int ksta  = ibcinfo(ibl,iseg,1,6);
                int kend  = ibcinfo(ibl,iseg,1,7);
                int nface = ibcinfo(ibl,iseg,1,8);
                bc1000(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem);
            }
        }
    }
}

void bc1000(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
            FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si, FortranArray3DRef<double> bcj,
            FortranArray3DRef<double> bck, FortranArray3DRef<double> bci,
            FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk,
            FortranArray4DRef<double> xtbi, FortranArray4DRef<double> atbj,
            FortranArray4DRef<double> atbk, FortranArray4DRef<double> atbi,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            int& nface, FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0, FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0, int& iuns,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& nummem)
{
    // COMMON block aliases
    float* qiv   = cmn_ivals.qiv;    // qiv(5), 0-based
    float* tur10 = cmn_ivals.tur10;  // tur10(7), 0-based
    int32_t& level   = cmn_mgrd.level;
    int32_t& lglobal = cmn_mgrd.lglobal;
    int32_t* ivisc   = cmn_reyue.ivisc; // ivisc(3), 0-based

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    //******************************************************************************
    //      j=1 boundary            freestream                           bctype 1000
    //******************************************************************************

    if (nface == 3) {
        for (int l = 1; l <= 5; l++) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    qj0(k,i,l,1) = qiv[l-1];
                    qj0(k,i,l,2) = qj0(k,i,l,1);
                    bcj(k,i,1)   = 0.0;
                }
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
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        for (int k = ksta; k <= kend1; k++) {
                            tj0(k,i,l,1) = tur10[l-1];
                            tj0(k,i,l,2) = tur10[l-1];
                        }
                    }
                }
            }
        }
    }

    //******************************************************************************
    //      j=jdim boundary             freestream                       bctype 1000
    //******************************************************************************

    if (nface == 4) {
        for (int l = 1; l <= 5; l++) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    qj0(k,i,l,3) = qiv[l-1];
                    qj0(k,i,l,4) = qj0(k,i,l,3);
                    bcj(k,i,2)   = 0.0;
                }
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
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        for (int k = ksta; k <= kend1; k++) {
                            tj0(k,i,l,3) = tur10[l-1];
                            tj0(k,i,l,4) = tur10[l-1];
                        }
                    }
                }
            }
        }
    }

    //******************************************************************************
    //      k=1 boundary            freestream                           bctype 1000
    //******************************************************************************

    if (nface == 5) {
        for (int l = 1; l <= 5; l++) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    qk0(j,i,l,1) = qiv[l-1];
                    qk0(j,i,l,2) = qk0(j,i,l,1);
                    bck(j,i,1)   = 0.0;
                }
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
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        for (int j = jsta; j <= jend1; j++) {
                            tk0(j,i,l,1) = tur10[l-1];
                            tk0(j,i,l,2) = tur10[l-1];
                        }
                    }
                }
            }
        }
    }

    //******************************************************************************
    //      k=kdim boundary             freestream                       bctype 1000
    //******************************************************************************

    if (nface == 6) {
        for (int l = 1; l <= 5; l++) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    qk0(j,i,l,3) = qiv[l-1];
                    qk0(j,i,l,4) = qk0(j,i,l,3);
                    bck(j,i,2)   = 0.0;
                }
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
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        for (int j = jsta; j <= jend1; j++) {
                            tk0(j,i,l,3) = tur10[l-1];
                            tk0(j,i,l,4) = tur10[l-1];
                        }
                    }
                }
            }
        }
    }

    //******************************************************************************
    //      i=1 boundary            freestream                           bctype 1000
    //******************************************************************************

    if (nface == 1) {
        for (int l = 1; l <= 5; l++) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    qi0(j,k,l,1) = qiv[l-1];
                    qi0(j,k,l,2) = qi0(j,k,l,1);
                    bci(j,k,1)   = 0.0;
                }
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
                for (int l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        for (int j = jsta; j <= jend1; j++) {
                            ti0(j,k,l,1) = tur10[l-1];
                            ti0(j,k,l,2) = tur10[l-1];
                        }
                    }
                }
            }
        }
    }

    //******************************************************************************
    //      i=idim boundary             freestream                       bctype 1000
    //******************************************************************************

    if (nface == 2) {
        for (int l = 1; l <= 5; l++) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    qi0(j,k,l,3) = qiv[l-1];
                    qi0(j,k,l,4) = qi0(j,k,l,3);
                    bci(j,k,2)   = 0.0;
                }
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
                for (int l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        for (int j = jsta; j <= jend1; j++) {
                            ti0(j,k,l,3) = tur10[l-1];
                            ti0(j,k,l,4) = tur10[l-1];
                        }
                    }
                }
            }
        }
    }
}

} // namespace bc1000_ns
