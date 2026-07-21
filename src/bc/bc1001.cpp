// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc1001.h"
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <limits>

namespace bc1001_ns {

// -----------------------------------------------------------------------
// bc() dispatcher
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

        FortranArray4DRef<double> tursav(&w(lw2(1,ibl)), jdim,   kdim,   idim,   nummem);
        FortranArray4DRef<double> tj0   (&w(lw2(2,ibl)), kdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> tk0   (&w(lw2(3,ibl)), jdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> ti0   (&w(lw2(4,ibl)), jdim,   kdim,   nummem, 4);
        FortranArray3DRef<double> vist3d(&w(lw2(5,ibl)), jdim,   kdim,   idim);
        FortranArray4DRef<double> vj0   (&w(lw2(6,ibl)), kdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vk0   (&w(lw2(7,ibl)), jdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vi0   (&w(lw2(8,ibl)), jdim,   kdim,   1,      4);

        // x and z grid coordinates (lw(17) and lw(18) if available, else use zeros)
        // In CFL3D standard layout, x is at lw(17,ibl) and z at lw(18,ibl)
        FortranArray3DRef<double> x_arr (&w(lw(17,ibl)), jdim, kdim, idim);
        FortranArray3DRef<double> z_arr (&w(lw(18,ibl)), jdim, kdim, idim);

        // j-direction boundary conditions
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl,iseg,1,1) == 1001) {
                int ista  = jbcinfo(ibl,iseg,1,2);
                int iend  = jbcinfo(ibl,iseg,1,3);
                int jsta  = jbcinfo(ibl,iseg,1,4);
                int jend  = jbcinfo(ibl,iseg,1,5);
                int ksta  = jbcinfo(ibl,iseg,1,6);
                int kend  = jbcinfo(ibl,iseg,1,7);
                int nface = jbcinfo(ibl,iseg,1,8);
                bc1001(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem, x_arr, z_arr);
            }
        }

        // k-direction boundary conditions
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) == 1001) {
                int ista  = kbcinfo(ibl,iseg,1,2);
                int iend  = kbcinfo(ibl,iseg,1,3);
                int jsta  = kbcinfo(ibl,iseg,1,4);
                int jend  = kbcinfo(ibl,iseg,1,5);
                int ksta  = kbcinfo(ibl,iseg,1,6);
                int kend  = kbcinfo(ibl,iseg,1,7);
                int nface = kbcinfo(ibl,iseg,1,8);
                bc1001(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem, x_arr, z_arr);
            }
        }

        // i-direction boundary conditions
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) == 1001) {
                int ista  = ibcinfo(ibl,iseg,1,2);
                int iend  = ibcinfo(ibl,iseg,1,3);
                int jsta  = ibcinfo(ibl,iseg,1,4);
                int jend  = ibcinfo(ibl,iseg,1,5);
                int ksta  = ibcinfo(ibl,iseg,1,6);
                int kend  = ibcinfo(ibl,iseg,1,7);
                int nface = ibcinfo(ibl,iseg,1,8);
                bc1001(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, nou, bou, nbuf, ibufdim, nummem, x_arr, z_arr);
            }
        }
    }
}

// -----------------------------------------------------------------------
// Internal helper: get_n1n2n3
// Builds a local orthogonal coordinate system (rotation matrix rn)
// given the normal vector components (x1,x2,x3).
// rn is a 3x3 float array stored column-major: rn[col][row] (Fortran order)
// -----------------------------------------------------------------------
static void get_n1n2n3(float x1, float x2, float x3, float rn[3][3])
{
    // rn(:,1) = (x1,x2,x3) normalized  [column 1 = rows 1..3]
    rn[0][0] = x1; rn[1][0] = x2; rn[2][0] = x3;
    float rm = std::sqrt(rn[0][0]*rn[0][0] + rn[1][0]*rn[1][0] + rn[2][0]*rn[2][0]);
    rn[0][0] /= rm; rn[1][0] /= rm; rn[2][0] /= rm;

    // find index of smallest absolute value in column 1
    float rmin = std::numeric_limits<float>::max();
    int i_min = 0;
    for (int i = 0; i < 3; i++) {
        if (rmin > std::abs(rn[i][0])) {
            rmin = std::abs(rn[i][0]);
            i_min = i;
        }
    }

    // rn(:,3) = 0; rn(i_min,3) = 1
    rn[0][2] = 0.f; rn[1][2] = 0.f; rn[2][2] = 0.f;
    rn[i_min][2] = 1.f;

    // cross(rn(:,1), rn(:,3), rn(:,2))  → rn(:,2) = rn(:,1) x rn(:,3)
    {
        float v1[3] = {rn[0][0], rn[1][0], rn[2][0]};
        float v2[3] = {rn[0][2], rn[1][2], rn[2][2]};
        rn[0][1] = v1[1]*v2[2] - v1[2]*v2[1];
        rn[1][1] = v1[2]*v2[0] - v1[0]*v2[2];
        rn[2][1] = v1[0]*v2[1] - v1[1]*v2[0];
    }
    rm = std::sqrt(rn[0][1]*rn[0][1] + rn[1][1]*rn[1][1] + rn[2][1]*rn[2][1]);
    rn[0][1] /= rm; rn[1][1] /= rm; rn[2][1] /= rm;

    // cross(rn(:,1), rn(:,2), rn(:,3))  → rn(:,3) = rn(:,1) x rn(:,2)
    {
        float v1[3] = {rn[0][0], rn[1][0], rn[2][0]};
        float v2[3] = {rn[0][1], rn[1][1], rn[2][1]};
        rn[0][2] = v1[1]*v2[2] - v1[2]*v2[1];
        rn[1][2] = v1[2]*v2[0] - v1[0]*v2[2];
        rn[2][2] = v1[0]*v2[1] - v1[1]*v2[0];
    }
}

// -----------------------------------------------------------------------
// bc1001: symmetry plane boundary conditions
// -----------------------------------------------------------------------
void bc1001(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci,
            FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk,
            FortranArray4DRef<double> xtbi,
            FortranArray4DRef<double> atbj, FortranArray4DRef<double> atbk,
            FortranArray4DRef<double> atbi,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            int& nface,
            FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0,
            int& iuns,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& nummem,
            FortranArray3DRef<double> x, FortranArray3DRef<double> z)
{
    // COMMON block aliases
    int32_t& level   = cmn_mgrd.level;
    int32_t& lglobal = cmn_mgrd.lglobal;
    int32_t* ivisc   = cmn_reyue.ivisc;   // 0-based: ivisc[0..2]
    float&   xmach   = cmn_info.xmach;

    // specialtop_kmax1001 aliases
    int32_t& i_specialtop = cmn_specialtop_kmax1001.i_specialtop_kmax1001;
    float&   a_st         = cmn_specialtop_kmax1001.a_specialtop_kmax1001;
    float&   xc_st        = cmn_specialtop_kmax1001.xc_specialtop_kmax1001;
    float&   sig_st       = cmn_specialtop_kmax1001.sig_specialtop_kmax1001;
    float&   vtp_st       = cmn_specialtop_kmax1001.vtp_specialtop_kmax1001;
    float&   wc_st        = cmn_specialtop_kmax1001.wc_specialtop_kmax1001;
    float&   fac_st       = cmn_specialtop_kmax1001.fac_specialtop_kmax1001;
    float&   cc_st        = cmn_specialtop_kmax1001.cc_specialtop_kmax1001;
    float&   xerf_st      = cmn_specialtop_kmax1001.xerf_specialtop_kmax1001;
    float&   sigerf_st    = cmn_specialtop_kmax1001.sigerf_specialtop_kmax1001;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    // local rotation matrices (float, Fortran real)
    // rn[col][row] in C (column-major like Fortran)
    float rn[3][3], rnt[3][3];
    float tloc[6], tghost[6];



    // ******************************************************************
    // j=1 boundary   symmetry plane   bctype 1001
    // ******************************************************************
    if (nface == 3) {

        for (int i = ista; i <= iend1; i++) {
        for (int k = ksta; k <= kend1; k++) {

            double vcont1 = q(1,k,i,2)*sj(1,k,i,1) +
                            q(1,k,i,3)*sj(1,k,i,2) +
                            q(1,k,i,4)*sj(1,k,i,3) + sj(1,k,i,5);
            double vcont2 = q(2,k,i,2)*sj(1,k,i,1) +
                            q(2,k,i,3)*sj(1,k,i,2) +
                            q(2,k,i,4)*sj(1,k,i,3) + sj(1,k,i,5);

            qj0(k,i,1,1) = q(1,k,i,1);
            qj0(k,i,2,1) = q(1,k,i,2) - 2.*vcont1*sj(1,k,i,1);
            qj0(k,i,3,1) = q(1,k,i,3) - 2.*vcont1*sj(1,k,i,2);
            qj0(k,i,4,1) = q(1,k,i,4) - 2.*vcont1*sj(1,k,i,3);
            qj0(k,i,5,1) = q(1,k,i,5);

            qj0(k,i,1,2) = q(2,k,i,1);
            qj0(k,i,2,2) = q(2,k,i,2) - 2.*vcont2*sj(1,k,i,1);
            qj0(k,i,3,2) = q(2,k,i,3) - 2.*vcont2*sj(1,k,i,2);
            qj0(k,i,4,2) = q(2,k,i,4) - 2.*vcont2*sj(1,k,i,3);
            qj0(k,i,5,2) = q(2,k,i,5);

            bcj(k,i,1) = 0.0;

        } // k
        } // i

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int k = ksta; k <= kend1; k++) {
                vj0(k,i,1,1) = vist3d(1,k,i);
                vj0(k,i,1,2) = vist3d(2,k,i);
            }
            }
        }

        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    // get local orthogonal coordinate system of sym-plane
                    get_n1n2n3((float)sj(1,k,i,1), (float)sj(1,k,i,2),
                               (float)sj(1,k,i,3), rn);
                    // rnt = transpose of rn: rnt(:,col) = rn(col,:)
                    // Fortran: rnt(:,1)=rn(1,:); rnt(:,2)=rn(2,:); rnt(:,3)=rn(3,:)
                    // In our storage rn[col][row], rn(i,j) = rn[j-1][i-1]
                    // rnt(:,1) = rn(1,:) means rnt(row,1) = rn(1,row) for row=1..3
                    // i.e. rnt[0][row-1] = rn[row-1][0]
                    for (int ii = 0; ii < 3; ii++) {
                        rnt[0][ii] = rn[ii][0];
                        rnt[1][ii] = rn[ii][1];
                        rnt[2][ii] = rn[ii][2];
                    }

                    // tensor_rotate(tursav(1,k,i,1:6), rn, rnt, tloc)
                    {
                        float t6[6];
                        for (int ll = 1; ll <= 6; ll++) t6[ll-1] = (float)tursav(1,k,i,ll);
                        FortranArray1DRef<double> t6_ref(reinterpret_cast<double*>(t6), 6);
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> tloc_d(6), tghost_d(6);
                        FortranArray2DRef<double> rn_ref  = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref = rnt_d.ref();
                        FortranArray1DRef<double> tloc_ref  = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        // build t6 as double ref
                        FortranArray1D<double> t6_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = (double)t6[ll-1];
                        FortranArray1DRef<double> t6_d_ref = t6_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) tj0(k,i,ll,1) = tloc_d(ll);
                    }
                    tj0(k,i,7,1) = tursav(1,k,i,7);

                    // tensor_rotate(tursav(2,k,i,1:6), rn, rnt, tloc)
                    {
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> t6_d(6), tloc_d(6), tghost_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = tursav(2,k,i,ll);
                        FortranArray1DRef<double> t6_d_ref   = t6_d.ref();
                        FortranArray1DRef<double> tloc_ref   = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        FortranArray2DRef<double> rn_ref     = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref    = rnt_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) tj0(k,i,ll,2) = tloc_d(ll);
                    }
                    tj0(k,i,7,2) = tursav(2,k,i,7);

                } // k
                } // i

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    tj0(k,i,l,1) = tursav(1,k,i,l);
                    tj0(k,i,l,2) = tursav(2,k,i,l);
                }
                }
                }
            }
        } // level >= lglobal

    } // nface == 3



    // ******************************************************************
    // j=jdim boundary   symmetry plane   bctype 1001
    // ******************************************************************
    if (nface == 4) {

        for (int i = ista; i <= iend1; i++) {
        for (int k = ksta; k <= kend1; k++) {

            double vcont1 = q(jdim-1,k,i,2)*sj(jdim,k,i,1) +
                            q(jdim-1,k,i,3)*sj(jdim,k,i,2) +
                            q(jdim-1,k,i,4)*sj(jdim,k,i,3) + sj(jdim,k,i,5);
            double vcont2 = q(jdim-2,k,i,2)*sj(jdim,k,i,1) +
                            q(jdim-2,k,i,3)*sj(jdim,k,i,2) +
                            q(jdim-2,k,i,4)*sj(jdim,k,i,3) + sj(jdim,k,i,5);

            qj0(k,i,1,3) = q(jdim-1,k,i,1);
            qj0(k,i,2,3) = q(jdim-1,k,i,2) - 2.*vcont1*sj(jdim,k,i,1);
            qj0(k,i,3,3) = q(jdim-1,k,i,3) - 2.*vcont1*sj(jdim,k,i,2);
            qj0(k,i,4,3) = q(jdim-1,k,i,4) - 2.*vcont1*sj(jdim,k,i,3);
            qj0(k,i,5,3) = q(jdim-1,k,i,5);

            qj0(k,i,1,4) = q(jdim-2,k,i,1);
            qj0(k,i,2,4) = q(jdim-2,k,i,2) - 2.*vcont2*sj(jdim,k,i,1);
            qj0(k,i,3,4) = q(jdim-2,k,i,3) - 2.*vcont2*sj(jdim,k,i,2);
            qj0(k,i,4,4) = q(jdim-2,k,i,4) - 2.*vcont2*sj(jdim,k,i,3);
            qj0(k,i,5,4) = q(jdim-2,k,i,5);

            bcj(k,i,2) = 0.0;

        } // k
        } // i

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int k = ksta; k <= kend1; k++) {
                vj0(k,i,1,3) = vist3d(jdim-1,k,i);
                vj0(k,i,1,4) = vist3d(jdim-2,k,i);
            }
            }
        }

        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    get_n1n2n3((float)sj(jdim,k,i,1), (float)sj(jdim,k,i,2),
                               (float)sj(jdim,k,i,3), rn);
                    for (int ii = 0; ii < 3; ii++) {
                        rnt[0][ii] = rn[ii][0];
                        rnt[1][ii] = rn[ii][1];
                        rnt[2][ii] = rn[ii][2];
                    }

                    {
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> t6_d(6), tloc_d(6), tghost_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = tursav(jdim-1,k,i,ll);
                        FortranArray1DRef<double> t6_d_ref   = t6_d.ref();
                        FortranArray1DRef<double> tloc_ref   = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        FortranArray2DRef<double> rn_ref     = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref    = rnt_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) tj0(k,i,ll,3) = tloc_d(ll);
                    }
                    tj0(k,i,7,3) = tursav(jdim-1,k,i,7);

                    {
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> t6_d(6), tloc_d(6), tghost_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = tursav(jdim-2,k,i,ll);
                        FortranArray1DRef<double> t6_d_ref   = t6_d.ref();
                        FortranArray1DRef<double> tloc_ref   = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        FortranArray2DRef<double> rn_ref     = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref    = rnt_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) tj0(k,i,ll,4) = tloc_d(ll);
                    }
                    tj0(k,i,7,4) = tursav(jdim-2,k,i,7);

                } // k
                } // i

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    tj0(k,i,l,3) = tursav(jdim-1,k,i,l);
                    tj0(k,i,l,4) = tursav(jdim-2,k,i,l);
                }
                }
                }
            }
        } // level >= lglobal

    } // nface == 4



    // ******************************************************************
    // k=1 boundary   symmetry plane   bctype 1001
    // ******************************************************************
    if (nface == 5) {

        for (int i = ista; i <= iend1; i++) {
        for (int j = jsta; j <= jend1; j++) {

            double wcont1 = q(j,1,i,2)*sk(j,1,i,1) +
                            q(j,1,i,3)*sk(j,1,i,2) +
                            q(j,1,i,4)*sk(j,1,i,3) + sk(j,1,i,5);
            double wcont2 = q(j,2,i,2)*sk(j,1,i,1) +
                            q(j,2,i,3)*sk(j,1,i,2) +
                            q(j,2,i,4)*sk(j,1,i,3) + sk(j,1,i,5);

            qk0(j,i,1,1) = q(j,1,i,1);
            qk0(j,i,2,1) = q(j,1,i,2) - 2.*wcont1*sk(j,1,i,1);
            qk0(j,i,3,1) = q(j,1,i,3) - 2.*wcont1*sk(j,1,i,2);
            qk0(j,i,4,1) = q(j,1,i,4) - 2.*wcont1*sk(j,1,i,3);
            qk0(j,i,5,1) = q(j,1,i,5);

            qk0(j,i,1,2) = q(j,2,i,1);
            qk0(j,i,2,2) = q(j,2,i,2) - 2.*wcont2*sk(j,1,i,1);
            qk0(j,i,3,2) = q(j,2,i,3) - 2.*wcont2*sk(j,1,i,2);
            qk0(j,i,4,2) = q(j,2,i,4) - 2.*wcont2*sk(j,1,i,3);
            qk0(j,i,5,2) = q(j,2,i,5);

            bck(j,i,1) = 0.0;

        } // j
        } // i

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int j = jsta; j <= jend1; j++) {
                vk0(j,i,1,1) = vist3d(j,1,i);
                vk0(j,i,1,2) = vist3d(j,2,i);
            }
            }
        }

        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    get_n1n2n3((float)sk(j,1,i,1), (float)sk(j,1,i,2),
                               (float)sk(j,1,i,3), rn);
                    for (int ii = 0; ii < 3; ii++) {
                        rnt[0][ii] = rn[ii][0];
                        rnt[1][ii] = rn[ii][1];
                        rnt[2][ii] = rn[ii][2];
                    }

                    {
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> t6_d(6), tloc_d(6), tghost_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = tursav(j,1,i,ll);
                        FortranArray1DRef<double> t6_d_ref   = t6_d.ref();
                        FortranArray1DRef<double> tloc_ref   = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        FortranArray2DRef<double> rn_ref     = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref    = rnt_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) tk0(j,i,ll,1) = tloc_d(ll);
                    }
                    tk0(j,i,7,1) = tursav(j,1,i,7);

                    {
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> t6_d(6), tloc_d(6), tghost_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = tursav(j,2,i,ll);
                        FortranArray1DRef<double> t6_d_ref   = t6_d.ref();
                        FortranArray1DRef<double> tloc_ref   = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        FortranArray2DRef<double> rn_ref     = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref    = rnt_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) tk0(j,i,ll,2) = tloc_d(ll);
                    }
                    tk0(j,i,7,2) = tursav(j,2,i,7);

                } // j
                } // i

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    tk0(j,i,l,1) = tursav(j,1,i,l);
                    tk0(j,i,l,2) = tursav(j,2,i,l);
                }
                }
                }
            }
        } // level >= lglobal

    } // nface == 5



    // ******************************************************************
    // k=kdim boundary   symmetry plane   bctype 1001
    // ******************************************************************
    if (nface == 6) {

        double a     = (double)a_st;
        double xc    = (double)xc_st;
        double sig   = (double)sig_st;
        double vtp   = (double)vtp_st;
        double wc    = (double)wc_st;
        double fac   = (double)fac_st;
        double cc    = (double)cc_st;
        double xerf  = (double)xerf_st;
        double sigerf= (double)sigerf_st;

        for (int i = ista; i <= iend1; i++) {
        for (int j = jsta; j <= jend1; j++) {

            double xxx = 0.0, b = 0.0, dz = 0.0;
            if (i_specialtop == 1) {
                xxx = 0.5*(x(j,kdim,i) + x(j+1,kdim,i));
                b   = (xc - xxx) / sig;
                dz  = z(j,kdim,i) - z(j,kdim-1,i);
            }

            double wcont1 = q(j,kdim-1,i,2)*sk(j,kdim,i,1) +
                            q(j,kdim-1,i,3)*sk(j,kdim,i,2) +
                            q(j,kdim-1,i,4)*sk(j,kdim,i,3) + sk(j,kdim,i,5);
            double wcont2 = q(j,kdim-2,i,2)*sk(j,kdim,i,1) +
                            q(j,kdim-2,i,3)*sk(j,kdim,i,2) +
                            q(j,kdim-2,i,4)*sk(j,kdim,i,3) + sk(j,kdim,i,5);

            qk0(j,i,1,3) = q(j,kdim-1,i,1);
            if (i_specialtop == 1) {
                double vtop = (fac*sig*a*b + (1.-fac)*a)*std::exp(cc-(b*b))*(double)xmach +
                              vtp*(0.5*(1.-std::erf((xxx-xerf)/sigerf)))*(double)xmach;
                double dudy = fac*a*(2.0*(b*b)-1.)*std::exp(cc-(b*b))*(double)xmach +
                              (1.-fac)*a*2.0/sig*b*std::exp(cc-(b*b))*(double)xmach;
                qk0(j,i,2,3) = q(j,kdim-1,i,2) + dz*dudy;
                if (wc == 0.) {
                    qk0(j,i,3,3) = q(j,kdim-1,i,3) - 2.*wcont1*sk(j,kdim,i,2);
                } else {
                    qk0(j,i,3,3) = 2.*wc*(double)xmach - q(j,kdim-1,i,3);
                }
                qk0(j,i,4,3) = 2.*vtop - q(j,kdim-1,i,4);
                qk0(j,i,5,3) = q(j,kdim-1,i,5);
            } else {
                qk0(j,i,2,3) = q(j,kdim-1,i,2) - 2.*wcont1*sk(j,kdim,i,1);
                qk0(j,i,3,3) = q(j,kdim-1,i,3) - 2.*wcont1*sk(j,kdim,i,2);
                qk0(j,i,4,3) = q(j,kdim-1,i,4) - 2.*wcont1*sk(j,kdim,i,3);
                qk0(j,i,5,3) = q(j,kdim-1,i,5);
            }

            qk0(j,i,1,4) = q(j,kdim-2,i,1);
            if (i_specialtop == 1) {
                qk0(j,i,2,4) = 2.*qk0(j,i,2,3) - q(j,kdim-1,i,2);
                if (wc == 0.) {
                    qk0(j,i,3,4) = q(j,kdim-2,i,3) - 2.*wcont2*sk(j,kdim,i,2);
                } else {
                    qk0(j,i,3,4) = 2.*qk0(j,i,3,3) - q(j,kdim-1,i,3);
                }
                qk0(j,i,4,4) = 2.*qk0(j,i,4,3) - q(j,kdim-1,i,4);
            } else {
                qk0(j,i,2,4) = q(j,kdim-2,i,2) - 2.*wcont2*sk(j,kdim,i,1);
                qk0(j,i,3,4) = q(j,kdim-2,i,3) - 2.*wcont2*sk(j,kdim,i,2);
                qk0(j,i,4,4) = q(j,kdim-2,i,4) - 2.*wcont2*sk(j,kdim,i,3);
            }
            qk0(j,i,5,4) = q(j,kdim-2,i,5);

            bck(j,i,2) = 0.0;

        } // j
        } // i

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int j = jsta; j <= jend1; j++) {
                vk0(j,i,1,3) = vist3d(j,kdim-1,i);
                vk0(j,i,1,4) = vist3d(j,kdim-2,i);
            }
            }
        }

        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    get_n1n2n3((float)sk(j,kdim,i,1), (float)sk(j,kdim,i,2),
                               (float)sk(j,kdim,i,3), rn);
                    for (int ii = 0; ii < 3; ii++) {
                        rnt[0][ii] = rn[ii][0];
                        rnt[1][ii] = rn[ii][1];
                        rnt[2][ii] = rn[ii][2];
                    }

                    {
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> t6_d(6), tloc_d(6), tghost_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = tursav(j,kdim-1,i,ll);
                        FortranArray1DRef<double> t6_d_ref   = t6_d.ref();
                        FortranArray1DRef<double> tloc_ref   = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        FortranArray2DRef<double> rn_ref     = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref    = rnt_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) tk0(j,i,ll,3) = tloc_d(ll);
                    }
                    tk0(j,i,7,3) = tursav(j,kdim-1,i,7);

                    {
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> t6_d(6), tloc_d(6), tghost_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = tursav(j,kdim-2,i,ll);
                        FortranArray1DRef<double> t6_d_ref   = t6_d.ref();
                        FortranArray1DRef<double> tloc_ref   = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        FortranArray2DRef<double> rn_ref     = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref    = rnt_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) tk0(j,i,ll,4) = tloc_d(ll);
                    }
                    tk0(j,i,7,4) = tursav(j,kdim-2,i,7);

                } // j
                } // i

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    tk0(j,i,l,3) = tursav(j,kdim-1,i,l);
                    tk0(j,i,l,4) = tursav(j,kdim-2,i,l);
                }
                }
                }
            }
        } // level >= lglobal

    } // nface == 6



    // ******************************************************************
    // i=1 boundary   symmetry plane   bctype 1001
    // ******************************************************************
    if (nface == 1) {

        int i2 = std::min(2, idim1);

        for (int k = ksta; k <= kend1; k++) {
        for (int j = jsta; j <= jend1; j++) {

            double ucont1 = q(j,k,1,2)*si(j,k,1,1) +
                            q(j,k,1,3)*si(j,k,1,2) +
                            q(j,k,1,4)*si(j,k,1,3) + si(j,k,1,5);
            double ucont2 = q(j,k,i2,2)*si(j,k,1,1) +
                            q(j,k,i2,3)*si(j,k,1,2) +
                            q(j,k,i2,4)*si(j,k,1,3) + si(j,k,1,5);

            qi0(j,k,1,1) = q(j,k,1,1);
            qi0(j,k,2,1) = q(j,k,1,2) - 2.*ucont1*si(j,k,1,1);
            qi0(j,k,3,1) = q(j,k,1,3) - 2.*ucont1*si(j,k,1,2);
            qi0(j,k,4,1) = q(j,k,1,4) - 2.*ucont1*si(j,k,1,3);
            qi0(j,k,5,1) = q(j,k,1,5);

            qi0(j,k,1,2) = q(j,k,i2,1);
            qi0(j,k,2,2) = q(j,k,i2,2) - 2.*ucont2*si(j,k,1,1);
            qi0(j,k,3,2) = q(j,k,i2,3) - 2.*ucont2*si(j,k,1,2);
            qi0(j,k,4,2) = q(j,k,i2,4) - 2.*ucont2*si(j,k,1,3);
            qi0(j,k,5,2) = q(j,k,i2,5);

            bci(j,k,1) = 0.0;

        } // j
        } // k

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
            for (int j = jsta; j <= jend1; j++) {
                vi0(j,k,1,1) = vist3d(j,k,1);
                vi0(j,k,1,2) = vist3d(j,k,i2);
            }
            }
        }

        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    get_n1n2n3((float)si(j,k,1,1), (float)si(j,k,1,2),
                               (float)si(j,k,1,3), rn);
                    for (int ii = 0; ii < 3; ii++) {
                        rnt[0][ii] = rn[ii][0];
                        rnt[1][ii] = rn[ii][1];
                        rnt[2][ii] = rn[ii][2];
                    }

                    {
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> t6_d(6), tloc_d(6), tghost_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = tursav(j,k,1,ll);
                        FortranArray1DRef<double> t6_d_ref   = t6_d.ref();
                        FortranArray1DRef<double> tloc_ref   = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        FortranArray2DRef<double> rn_ref     = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref    = rnt_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) ti0(j,k,ll,1) = tloc_d(ll);
                    }
                    ti0(j,k,7,1) = tursav(j,k,1,7);

                    {
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> t6_d(6), tloc_d(6), tghost_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = tursav(j,k,i2,ll);
                        FortranArray1DRef<double> t6_d_ref   = t6_d.ref();
                        FortranArray1DRef<double> tloc_ref   = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        FortranArray2DRef<double> rn_ref     = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref    = rnt_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) ti0(j,k,ll,2) = tloc_d(ll);
                    }
                    ti0(j,k,7,2) = tursav(j,k,i2,7);

                } // j
                } // k

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    ti0(j,k,l,1) = tursav(j,k,1,l);
                    ti0(j,k,l,2) = tursav(j,k,i2,l);
                }
                }
                }
            }
        } // level >= lglobal

    } // nface == 1



    // ******************************************************************
    // i=idim boundary   symmetry plane   bctype 1001
    // ******************************************************************
    if (nface == 2) {

        int i2 = std::max(1, idim-2);

        for (int k = ksta; k <= kend1; k++) {
        for (int j = jsta; j <= jend1; j++) {

            double ucont1 = q(j,k,idim-1,2)*si(j,k,idim,1) +
                            q(j,k,idim-1,3)*si(j,k,idim,2) +
                            q(j,k,idim-1,4)*si(j,k,idim,3) + si(j,k,idim,5);
            double ucont2 = q(j,k,i2,2)*si(j,k,idim,1) +
                            q(j,k,i2,3)*si(j,k,idim,2) +
                            q(j,k,i2,4)*si(j,k,idim,3) + si(j,k,idim,5);

            qi0(j,k,1,3) = q(j,k,idim-1,1);
            qi0(j,k,2,3) = q(j,k,idim-1,2) - 2.*ucont1*si(j,k,idim,1);
            qi0(j,k,3,3) = q(j,k,idim-1,3) - 2.*ucont1*si(j,k,idim,2);
            qi0(j,k,4,3) = q(j,k,idim-1,4) - 2.*ucont1*si(j,k,idim,3);
            qi0(j,k,5,3) = q(j,k,idim-1,5);

            qi0(j,k,1,4) = q(j,k,i2,1);
            qi0(j,k,2,4) = q(j,k,i2,2) - 2.*ucont2*si(j,k,idim,1);
            qi0(j,k,3,4) = q(j,k,i2,3) - 2.*ucont2*si(j,k,idim,2);
            qi0(j,k,4,4) = q(j,k,i2,4) - 2.*ucont2*si(j,k,idim,3);
            qi0(j,k,5,4) = q(j,k,i2,5);

            bci(j,k,2) = 0.0;

        } // j
        } // k

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
            for (int j = jsta; j <= jend1; j++) {
                vi0(j,k,1,3) = vist3d(j,k,idim-1);
                vi0(j,k,1,4) = vist3d(j,k,i2);
            }
            }
        }

        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    get_n1n2n3((float)si(j,k,idim-1,1), (float)si(j,k,idim-1,2),
                               (float)si(j,k,idim-1,3), rn);
                    for (int ii = 0; ii < 3; ii++) {
                        rnt[0][ii] = rn[ii][0];
                        rnt[1][ii] = rn[ii][1];
                        rnt[2][ii] = rn[ii][2];
                    }

                    {
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> t6_d(6), tloc_d(6), tghost_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = tursav(j,k,idim-1,ll);
                        FortranArray1DRef<double> t6_d_ref   = t6_d.ref();
                        FortranArray1DRef<double> tloc_ref   = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        FortranArray2DRef<double> rn_ref     = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref    = rnt_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) ti0(j,k,ll,3) = tloc_d(ll);
                    }
                    ti0(j,k,7,3) = tursav(j,k,idim-1,7);

                    {
                        FortranArray2D<double> rn_d(3,3), rnt_d(3,3);
                        for (int c=0;c<3;c++) for (int r=0;r<3;r++) {
                            rn_d(r+1,c+1)  = (double)rn[c][r];
                            rnt_d(r+1,c+1) = (double)rnt[c][r];
                        }
                        FortranArray1D<double> t6_d(6), tloc_d(6), tghost_d(6);
                        for (int ll=1;ll<=6;ll++) t6_d(ll) = tursav(j,k,i2,ll);
                        FortranArray1DRef<double> t6_d_ref   = t6_d.ref();
                        FortranArray1DRef<double> tloc_ref   = tloc_d.ref();
                        FortranArray1DRef<double> tghost_ref = tghost_d.ref();
                        FortranArray2DRef<double> rn_ref     = rn_d.ref();
                        FortranArray2DRef<double> rnt_ref    = rnt_d.ref();
                        tensor_rotate(t6_d_ref, rn_ref, rnt_ref, tloc_ref);
                        tensor_sym(tloc_ref, tghost_ref);
                        tensor_rotate(tghost_ref, rnt_ref, rn_ref, tloc_ref);
                        for (int ll=1;ll<=6;ll++) ti0(j,k,ll,4) = tloc_d(ll);
                    }
                    ti0(j,k,7,4) = tursav(j,k,i2,7);

                } // j
                } // k

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    ti0(j,k,l,3) = tursav(j,k,idim-1,l);
                    ti0(j,k,l,4) = tursav(j,k,i2,l);
                }
                }
                }
            }
        } // level >= lglobal

    } // nface == 2

} // bc1001



// -----------------------------------------------------------------------
// cross: cross product v3 = v1 x v2
// -----------------------------------------------------------------------
void cross(FortranArray1DRef<double> v1, FortranArray1DRef<double> v2,
           FortranArray1DRef<double> v3)
{
    v3(1) = v1(2)*v2(3) - v1(3)*v2(2);
    v3(2) = v1(3)*v2(1) - v1(1)*v2(3);
    v3(3) = v1(1)*v2(2) - v1(2)*v2(1);
}

// -----------------------------------------------------------------------
// tensor_rotate: rotate stress tensor t(6) using rotation matrices rn, rnt
// Result stored in tloc(6)
// Fortran: s(:,j) = sum_m t33(:,m)*rn(m,j)  (matrix multiply T*RN)
//          v(:,j) = sum_m rnt(:,m)*s(m,j)    (matrix multiply RNT*S)
// -----------------------------------------------------------------------
void tensor_rotate(FortranArray1DRef<double> t, FortranArray2DRef<double> rn,
                   FortranArray2DRef<double> rnt, FortranArray1DRef<double> tloc)
{
    // Build symmetric 3x3 tensor from Voigt notation
    // Fortran column-major: t33(row,col) → C: t33[col-1][row-1]
    double t33[3][3];
    t33[0][0] = t(1);  // t33(1,1)
    t33[1][1] = t(2);  // t33(2,2)
    t33[2][2] = t(3);  // t33(3,3)
    t33[1][0] = t(4);  // t33(1,2) — col=2, row=1
    t33[2][1] = t(5);  // t33(2,3) — col=3, row=2
    t33[2][0] = t(6);  // t33(1,3) — col=3, row=1
    t33[0][1] = t33[1][0];  // t33(2,1)
    t33[1][2] = t33[2][1];  // t33(3,2)
    t33[0][2] = t33[2][0];  // t33(3,1)

    // s(:,col) = t33(:,1)*rn(1,col) + t33(:,2)*rn(2,col) + t33(:,3)*rn(3,col)
    // s[col-1][row-1] = sum_m t33[m][row-1] * rn(m+1,col)
    double s[3][3];
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            s[col][row] = 0.0;
            for (int m = 0; m < 3; m++) {
                s[col][row] += t33[m][row] * rn(m+1, col+1);
            }
        }
    }

    // v(:,col) = rnt(:,1)*s(1,col) + rnt(:,2)*s(2,col) + rnt(:,3)*s(3,col)
    // v[col-1][row-1] = sum_m rnt(row+1,m+1) * s[col][m]
    double v[3][3];
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            v[col][row] = 0.0;
            for (int m = 0; m < 3; m++) {
                v[col][row] += rnt(row+1, m+1) * s[col][m];
            }
        }
    }

    tloc(1) = v[0][0];  // v(1,1)
    tloc(2) = v[1][1];  // v(2,2)
    tloc(3) = v[2][2];  // v(3,3)
    tloc(4) = v[1][0];  // v(1,2)
    tloc(5) = v[2][1];  // v(2,3)
    tloc(6) = v[2][0];  // v(1,3)
}

// -----------------------------------------------------------------------
// tensor_sym: apply symmetry condition to stress tensor (Voigt notation)
// -----------------------------------------------------------------------
void tensor_sym(FortranArray1DRef<double> tloc, FortranArray1DRef<double> tghost)
{
    tghost(1) =  tloc(1);
    tghost(2) =  tloc(2);
    tghost(3) =  tloc(3);
    tghost(4) = -tloc(4);
    tghost(5) =  tloc(5);
    tghost(6) = -tloc(6);
}

} // namespace bc1001_ns
