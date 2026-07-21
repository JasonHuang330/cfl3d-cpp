// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2005.h"
#include "rotateq.h"
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <cstring>
#include <algorithm>

// Forward-declare module_stm_2005 functions (header is empty but functions exist)
namespace module_stm_2005_ns {
    void stm2k5_get_rotmat(float& thx, float& thy, float& thz,
                           FortranArray2DRef<float> rn,
                           FortranArray2DRef<float> rnt);
    void stm2k5_bc(FortranArray1DRef<float> vin,
                   FortranArray2DRef<float> rn,
                   FortranArray2DRef<float> rnt,
                   FortranArray1DRef<float> vout);
}

namespace bc2005_ns {

// ============================================================
// bc() — dispatcher: loops over all blocks, finds bctype==2005
// ============================================================
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
    for (int ibl = 1; ibl <= nbl; ibl++) {
        int jdim = jdimg(ibl);
        int kdim = kdimg(ibl);
        int idim = idimg(ibl);

        // Build array views from w using lw offsets
        FortranArray4DRef<double> q   (&w(lw(1,ibl)), jdim,   kdim,   idim,   5);
        FortranArray4DRef<double> qj0 (&w(lw(2,ibl)), kdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qk0 (&w(lw(3,ibl)), jdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qi0 (&w(lw(4,ibl)), jdim,   kdim,   5,      4);

        // Build array views from wk using lw2 offsets
        FortranArray4DRef<double> tursav(&wk(lw2(1,ibl)), jdim,   kdim,   idim,   nummem);
        FortranArray4DRef<double> tj0   (&wk(lw2(2,ibl)), kdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> tk0   (&wk(lw2(3,ibl)), jdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> ti0   (&wk(lw2(4,ibl)), jdim,   kdim,   nummem, 4);
        FortranArray3DRef<double> vist3d(&wk(lw2(5,ibl)), jdim,   kdim,   idim);
        FortranArray4DRef<double> vj0   (&wk(lw2(6,ibl)), kdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vk0   (&wk(lw2(7,ibl)), jdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vi0   (&wk(lw2(8,ibl)), jdim,   kdim,   1,      4);

        // bcdata and filname from lwdat/bcfiles
        // mdim=nbcjdim(ibl), ndim=nbcjdim(ibl) — use maxseg for safety
        int mdim_j = nbcjdim(ibl);
        int ndim_j = nbcjdim(ibl);
        int mdim_k = nbckdim(ibl);
        int ndim_k = nbckdim(ibl);
        int mdim_i = nbcidim(ibl);
        int ndim_i = nbcidim(ibl);

        // ---- j-face segments ----
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl,iseg,1,1) != 2005) continue;
            int ista  = jbcinfo(ibl,iseg,1,2);
            int iend  = jbcinfo(ibl,iseg,1,3);
            int jsta  = jbcinfo(ibl,iseg,1,4);
            int jend  = jbcinfo(ibl,iseg,1,5);
            int ksta  = jbcinfo(ibl,iseg,1,6);
            int kend  = jbcinfo(ibl,iseg,1,7);
            int nface = jbcinfo(ibl,iseg,1,8);
            int nblp  = jbcinfo(ibl,iseg,1,9);

            int jdimp = jdimg(nblp);
            int kdimp = kdimg(nblp);
            int idimp = idimg(nblp);

            FortranArray4DRef<double> qp  (&w(lw(1,nblp)),    jdimp, kdimp, idimp, 5);
            FortranArray3DRef<double> vp  (&wk(lw2(5,nblp)),  jdimp, kdimp, idimp);
            FortranArray4DRef<double> tp  (&wk(lw2(1,nblp)),  jdimp, kdimp, idimp, nummem);

            // Work arrays for rotation
            FortranArray4D<double> qrotj(2, kdimp, idimp, 5);
            FortranArray4D<double> qrotk(jdimp, 2, idimp, 5);
            FortranArray4D<double> qroti(jdimp, kdimp, 2, 5);

            // bcdata
            int ibcf = bcfilej(ibl,iseg,1);
            FortranArray4DRef<double> bcdata(&w(lwdat(ibl,iseg,1)), mdim_j, ndim_j, 2, 12);
            char* filname_ptr = bcfiles(ibcf);

            bc2005(jdim, kdim, idim, q, qj0, qk0, qi0,
                   ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   mdim_j, ndim_j, bcdata, *reinterpret_cast<char(*)[80]>(filname_ptr),
                   qp, vp, tp, jdimp, kdimp, idimp,
                   qrotj.ref(), qrotk.ref(), qroti.ref(),
                   ibl, nblp, nou, bou, nbuf, ibufdim, myid,
                   nblg, maxbl, nummem);
        }

        // ---- k-face segments ----
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) != 2005) continue;
            int ista  = kbcinfo(ibl,iseg,1,2);
            int iend  = kbcinfo(ibl,iseg,1,3);
            int jsta  = kbcinfo(ibl,iseg,1,4);
            int jend  = kbcinfo(ibl,iseg,1,5);
            int ksta  = kbcinfo(ibl,iseg,1,6);
            int kend  = kbcinfo(ibl,iseg,1,7);
            int nface = kbcinfo(ibl,iseg,1,8);
            int nblp  = kbcinfo(ibl,iseg,1,9);

            int jdimp = jdimg(nblp);
            int kdimp = kdimg(nblp);
            int idimp = idimg(nblp);

            FortranArray4DRef<double> qp  (&w(lw(1,nblp)),    jdimp, kdimp, idimp, 5);
            FortranArray3DRef<double> vp  (&wk(lw2(5,nblp)),  jdimp, kdimp, idimp);
            FortranArray4DRef<double> tp  (&wk(lw2(1,nblp)),  jdimp, kdimp, idimp, nummem);

            FortranArray4D<double> qrotj(2, kdimp, idimp, 5);
            FortranArray4D<double> qrotk(jdimp, 2, idimp, 5);
            FortranArray4D<double> qroti(jdimp, kdimp, 2, 5);

            int ibcf = bcfilek(ibl,iseg,1);
            FortranArray4DRef<double> bcdata(&w(lwdat(ibl,iseg,2)), mdim_k, ndim_k, 2, 12);
            char* filname_ptr = bcfiles(ibcf);

            bc2005(jdim, kdim, idim, q, qj0, qk0, qi0,
                   ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   mdim_k, ndim_k, bcdata, *reinterpret_cast<char(*)[80]>(filname_ptr),
                   qp, vp, tp, jdimp, kdimp, idimp,
                   qrotj.ref(), qrotk.ref(), qroti.ref(),
                   ibl, nblp, nou, bou, nbuf, ibufdim, myid,
                   nblg, maxbl, nummem);
        }

        // ---- i-face segments ----
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) != 2005) continue;
            int ista  = ibcinfo(ibl,iseg,1,2);
            int iend  = ibcinfo(ibl,iseg,1,3);
            int jsta  = ibcinfo(ibl,iseg,1,4);
            int jend  = ibcinfo(ibl,iseg,1,5);
            int ksta  = ibcinfo(ibl,iseg,1,6);
            int kend  = ibcinfo(ibl,iseg,1,7);
            int nface = ibcinfo(ibl,iseg,1,8);
            int nblp  = ibcinfo(ibl,iseg,1,9);

            int jdimp = jdimg(nblp);
            int kdimp = kdimg(nblp);
            int idimp = idimg(nblp);

            FortranArray4DRef<double> qp  (&w(lw(1,nblp)),    jdimp, kdimp, idimp, 5);
            FortranArray3DRef<double> vp  (&wk(lw2(5,nblp)),  jdimp, kdimp, idimp);
            FortranArray4DRef<double> tp  (&wk(lw2(1,nblp)),  jdimp, kdimp, idimp, nummem);

            FortranArray4D<double> qrotj(2, kdimp, idimp, 5);
            FortranArray4D<double> qrotk(jdimp, 2, idimp, 5);
            FortranArray4D<double> qroti(jdimp, kdimp, 2, 5);

            int ibcf = bcfilei(ibl,iseg,1);
            FortranArray4DRef<double> bcdata(&w(lwdat(ibl,iseg,3)), mdim_i, ndim_i, 2, 12);
            char* filname_ptr = bcfiles(ibcf);

            bc2005(jdim, kdim, idim, q, qj0, qk0, qi0,
                   ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                   mdim_i, ndim_i, bcdata, *reinterpret_cast<char(*)[80]>(filname_ptr),
                   qp, vp, tp, jdimp, kdimp, idimp,
                   qrotj.ref(), qrotk.ref(), qroti.ref(),
                   ibl, nblp, nou, bou, nbuf, ibufdim, myid,
                   nblg, maxbl, nummem);
        }
    }
}

// ============================================================
// bc2005() — periodic BC with rotation
// ============================================================
void bc2005(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0,
            FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0,
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
            int& mdim, int& ndim,
            FortranArray4DRef<double> bcdata,
            char (&filname)[80],
            FortranArray4DRef<double> qp,
            FortranArray3DRef<double> vp,
            FortranArray4DRef<double> tp,
            int& jdimp, int& kdimp, int& idimp,
            FortranArray4DRef<double> qrotj,
            FortranArray4DRef<double> qrotk,
            FortranArray4DRef<double> qroti,
            int& nbl, int& nblp,
            FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim,
            int& myid,
            FortranArray1DRef<int> mblk2nd,
            int& maxbl, int& nummem)
{
    // COMMON block aliases
    int&   level   = cmn_mgrd.level;
    int&   lglobal = cmn_mgrd.lglobal;
    int*   ivisc   = cmn_reyue.ivisc;   // 0-based: ivisc[0],ivisc[1],ivisc[2]
    float& radtodeg = cmn_conversion.radtodeg;
    int&   ialph   = cmn_igrdtyp.ialph;

    // Local variables
    int jdim1, kdim1, idim1;
    int jend1, kend1, iend1;
    int ip;
    double dthtx, dthty, dthtz;
    FortranArray2D<float> rn(3,3), rnt(3,3);
    int jmx, kmx, imx;
    int jdm2, kdm2, idm2;
    int jdp2, kdp2, idp2;
    int l, i, j, k;

    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;

    jend1 = jend - 1;
    kend1 = kend - 1;
    iend1 = iend - 1;

    // The caller hands qrotj/qrotk/qroti in as flat work-array slices; this
    // routine's Fortran contract re-dimensions them per its own
    // `dimension qrotj(2,kdimp,idimp,5),qrotk(jdimp,2,idimp,5),
    //            qroti(jdimp,kdimp,2,5)`.
    // Re-view here so qrot*(...) index with the correct strides; otherwise the
    // rotated periodic-face values are scattered/overflowed and the ghost qi0
    // ends up zero (negative-density abort).
    qrotj = FortranArray4DRef<double>(qrotj.data(), 2,     kdimp, idimp, 5);
    qrotk = FortranArray4DRef<double>(qrotk.data(), jdimp, 2,     idimp, 5);
    qroti = FortranArray4DRef<double>(qroti.data(), jdimp, kdimp, 2,     5);

    // Likewise re-view the turbulent-viscosity boundary arrays to their Fortran
    // shapes: vj0(kdim,idim-1,1,4), vk0(jdim,idim-1,1,4), vi0(jdim,kdim,1,4).
    // bc_period hands these in with an over-large (jdim,kdim,idim,5) shape; a
    // write through that shape spills vi0 past its slot into the adjacent
    // bcj/bck/bci flag arrays, corrupting the I-direction periodic-face flags.
    vj0 = FortranArray4DRef<double>(vj0.data(), kdim, idim-1, 1, 4);
    vk0 = FortranArray4DRef<double>(vk0.data(), jdim, idim-1, 1, 4);
    vi0 = FortranArray4DRef<double>(vi0.data(), jdim, kdim, 1, 4);

    // Same for the mean-flow and turbulence ghost arrays — bc_period passes
    // qi0/qj0/qk0 as (jdim,kdim,idim,5) and tj0/tk0/ti0 with mismatched shapes,
    // but bc2005's Fortran contract is qi0(jdim,kdim,5,4), qj0(kdim,idim-1,5,4),
    // qk0(jdim,idim-1,5,4), tj0(kdim,idim-1,nummem,4), tk0(jdim,idim-1,nummem,4),
    // ti0(jdim,kdim,nummem,4). Without this, bc2005's qi0(...,3/4) writes (the
    // i=idim periodic ghost) land at the wrong offsets and never reach the
    // solver's reads — leaving a stale/garbage i=idim ghost, a spurious
    // I-direction inviscid flux, and a non-zero freestream residual that
    // diverges the run.
    qi0 = FortranArray4DRef<double>(qi0.data(), jdim, kdim, 5, 4);
    qj0 = FortranArray4DRef<double>(qj0.data(), kdim, idim-1, 5, 4);
    qk0 = FortranArray4DRef<double>(qk0.data(), jdim, idim-1, 5, 4);
    tj0 = FortranArray4DRef<double>(tj0.data(), kdim, idim-1, nummem, 4);
    tk0 = FortranArray4DRef<double>(tk0.data(), jdim, idim-1, nummem, 4);
    ti0 = FortranArray4DRef<double>(ti0.data(), jdim, kdim, nummem, 4);

    // this bc makes use of only one plane of data
    ip    = 1;
    dthtx = (double)bcdata(1,1,ip,2) / (double)radtodeg;
    if (ialph == 0) {
        dthty = (double)bcdata(1,1,ip,3) / (double)radtodeg;
        dthtz = (double)bcdata(1,1,ip,4) / (double)radtodeg;
    } else {
        dthty = -(double)bcdata(1,1,ip,4) / (double)radtodeg;
        dthtz =  (double)bcdata(1,1,ip,3) / (double)radtodeg;
    }



    // ******************************************************************************
    //      j=1 boundary             periodic boundary                   bctype 2005
    // ******************************************************************************

    if (nface == 3) {

        // Load qp values into qrotj array:
        for (l = 1; l <= 5; l++)
        for (i = ista; i <= iend1; i++)
        for (k = ksta; k <= kend1; k++) {
            qrotj(1,k,i,l) = qp(jdimp-1,k,i,l);
            jmx = std::max(jdimp-2, 1);
            qrotj(2,k,i,l) = qp(jmx,k,i,l);
        }

        // Rotate qrotj values
        {
            int jd2 = 2, kd_p = kdimp, id_p = idimp;
            if (jdimp == 2) {
                int is1=ista, ie1=iend1, js1=1, je1=1, ks1=ksta, ke1=kend1;
                rotateq_ns::rotateq(jd2, kd_p, id_p, qrotj, qrotj,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
                double dthtx2 = 2.0*dthtx, dthty2 = 2.0*dthty, dthtz2 = 2.0*dthtz;
                int js2=2, je2=2;
                rotateq_ns::rotateq(jd2, kd_p, id_p, qrotj, qrotj,
                                    is1, ie1, js2, je2, ks1, ke1,
                                    dthtx2, dthty2, dthtz2);
            } else {
                int is1=ista, ie1=iend1, js1=1, je1=2, ks1=ksta, ke1=kend1;
                rotateq_ns::rotateq(jd2, kd_p, id_p, qrotj, qrotj,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
            }
        }

        // Apply periodic BCs
        for (l = 1; l <= 5; l++)
        for (i = ista; i <= iend1; i++)
        for (k = ksta; k <= kend1; k++) {
            qj0(k,i,l,1) = qrotj(1,k,i,l);
            qj0(k,i,l,2) = qrotj(2,k,i,l);
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                vj0(k,i,1,1) = vp(jdimp-1,k,i);
                vj0(k,i,1,2) = 0.0;
            }
        }

        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                float fthx = (float)dthtx, fthy = (float)dthty, fthz = (float)dthtz;
                FortranArray2DRef<float> rn_ref = rn.ref(), rnt_ref = rnt.ref();
                module_stm_2005_ns::stm2k5_get_rotmat(fthx, fthy, fthz, rn_ref, rnt_ref);

                // rotate the stress tensors
                for (i = ista; i <= iend1; i++) {
                    for (k = ksta; k <= kend1; k++) {
                        {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(jdimp-1,k,i,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tj0(k,i,ll,1) = (double)tmp_out[ll-1];
                        }
                        if (jdimp == 2) {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tj0(k,i,ll,1);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tj0(k,i,ll,2) = (double)tmp_out[ll-1];
                        } else {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(jdimp-2,k,i,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tj0(k,i,ll,2) = (double)tmp_out[ll-1];
                        }
                    }
                }

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                jdm2 = std::max(jdimp-2, 1);
                for (l = 1; l <= nummem; l++)
                for (i = ista; i <= iend1; i++)
                for (k = ksta; k <= kend1; k++) {
                    tj0(k,i,l,1) = tp(jdimp-1,k,i,l);
                    tj0(k,i,l,2) = tp(jdm2,k,i,l);
                }
            }
        }

    } // end nface==3



    // ******************************************************************************
    //      j=jdim boundary          periodic boundary                   bctype 2005
    // ******************************************************************************

    if (nface == 4) {

        // Load qp values into qrotj array:
        for (l = 1; l <= 5; l++)
        for (i = ista; i <= iend1; i++)
        for (k = ksta; k <= kend1; k++) {
            qrotj(1,k,i,l) = qp(1,k,i,l);
            jmx = std::min(jdimp-1, 2);
            qrotj(2,k,i,l) = qp(jmx,k,i,l);
        }

        // Rotate qrotj values
        {
            int jd2 = 2, kd_p = kdimp, id_p = idimp;
            if (jdimp == 2) {
                int is1=ista, ie1=iend1, js1=1, je1=1, ks1=ksta, ke1=kend1;
                rotateq_ns::rotateq(jd2, kd_p, id_p, qrotj, qrotj,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
                double dthtx2 = 2.0*dthtx, dthty2 = 2.0*dthty, dthtz2 = 2.0*dthtz;
                int js2=2, je2=2;
                rotateq_ns::rotateq(jd2, kd_p, id_p, qrotj, qrotj,
                                    is1, ie1, js2, je2, ks1, ke1,
                                    dthtx2, dthty2, dthtz2);
            } else {
                int is1=ista, ie1=iend1, js1=1, je1=2, ks1=ksta, ke1=kend1;
                rotateq_ns::rotateq(jd2, kd_p, id_p, qrotj, qrotj,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
            }
        }

        // Apply periodic BCs
        for (l = 1; l <= 5; l++)
        for (i = ista; i <= iend1; i++)
        for (k = ksta; k <= kend1; k++) {
            qj0(k,i,l,3) = qrotj(1,k,i,l);
            qj0(k,i,l,4) = qrotj(2,k,i,l);
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                vj0(k,i,1,3) = vp(1,k,i);
                vj0(k,i,1,4) = 0.0;
            }
        }

        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                // obtain the rotation matrix and its transpose
                float fthx = (float)dthtx, fthy = (float)dthty, fthz = (float)dthtz;
                FortranArray2DRef<float> rn_ref = rn.ref(), rnt_ref = rnt.ref();
                module_stm_2005_ns::stm2k5_get_rotmat(fthx, fthy, fthz, rn_ref, rnt_ref);

                // rotate the stress tensors
                for (i = ista; i <= iend1; i++) {
                    for (k = ksta; k <= kend1; k++) {
                        {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(1,k,i,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tj0(k,i,ll,3) = (double)tmp_out[ll-1];
                        }
                        if (jdimp == 2) {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tj0(k,i,ll,3);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tj0(k,i,ll,4) = (double)tmp_out[ll-1];
                        } else {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(2,k,i,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tj0(k,i,ll,4) = (double)tmp_out[ll-1];
                        }
                    }
                }

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                jdp2 = std::min(2, jdimp-1);
                for (l = 1; l <= nummem; l++)
                for (i = ista; i <= iend1; i++)
                for (k = ksta; k <= kend1; k++) {
                    tj0(k,i,l,3) = tp(1,k,i,l);
                    tj0(k,i,l,4) = tp(jdp2,k,i,l);
                }
            }
        }

    } // end nface==4



    // ******************************************************************************
    //      k=1 boundary             periodic boundary                   bctype 2005
    // ******************************************************************************

    if (nface == 5) {

        // Load qp values into qrotk array:
        for (l = 1; l <= 5; l++)
        for (i = ista; i <= iend1; i++)
        for (j = jsta; j <= jend1; j++) {
            qrotk(j,1,i,l) = qp(j,kdimp-1,i,l);
            kmx = std::max(kdimp-2, 1);
            qrotk(j,2,i,l) = qp(j,kmx,i,l);
        }

        // Rotate qrotk values
        {
            int jd_p = jdimp, kd2 = 2, id_p = idimp;
            if (kdimp == 2) {
                int is1=ista, ie1=iend1, js1=jsta, je1=jend1, ks1=1, ke1=1;
                rotateq_ns::rotateq(jd_p, kd2, id_p, qrotk, qrotk,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
                double dthtx2 = 2.0*dthtx, dthty2 = 2.0*dthty, dthtz2 = 2.0*dthtz;
                int ks2=2, ke2=2;
                rotateq_ns::rotateq(jd_p, kd2, id_p, qrotk, qrotk,
                                    is1, ie1, js1, je1, ks2, ke2,
                                    dthtx2, dthty2, dthtz2);
            } else {
                int is1=ista, ie1=iend1, js1=jsta, je1=jend1, ks1=1, ke1=2;
                rotateq_ns::rotateq(jd_p, kd2, id_p, qrotk, qrotk,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
            }
        }

        // Apply periodic BCs
        for (l = 1; l <= 5; l++)
        for (i = ista; i <= iend1; i++)
        for (j = jsta; j <= jend1; j++) {
            qk0(j,i,l,1) = qrotk(j,1,i,l);
            qk0(j,i,l,2) = qrotk(j,2,i,l);
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                vk0(j,i,1,1) = vp(j,kdimp-1,i);
                vk0(j,i,1,2) = 0.0;
            }
        }

        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                // obtain the rotation matrix and its transpose
                float fthx = (float)dthtx, fthy = (float)dthty, fthz = (float)dthtz;
                FortranArray2DRef<float> rn_ref = rn.ref(), rnt_ref = rnt.ref();
                module_stm_2005_ns::stm2k5_get_rotmat(fthx, fthy, fthz, rn_ref, rnt_ref);

                // rotate the stress tensors
                for (i = ista; i <= iend1; i++) {
                    for (j = jsta; j <= jend1; j++) {
                        {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(j,kdimp-1,i,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tk0(j,i,ll,1) = (double)tmp_out[ll-1];
                        }
                        if (kdimp == 2) {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tk0(j,i,ll,1);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tk0(j,i,ll,2) = (double)tmp_out[ll-1];
                        } else {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(j,kdimp-2,i,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tk0(j,i,ll,2) = (double)tmp_out[ll-1];
                        }
                    }
                }

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                kdm2 = std::max(kdimp-2, 1);
                for (l = 1; l <= nummem; l++)
                for (i = ista; i <= iend1; i++)
                for (j = jsta; j <= jend1; j++) {
                    tk0(j,i,l,1) = tp(j,kdimp-1,i,l);
                    tk0(j,i,l,2) = tp(j,kdm2,i,l);
                }
            }
        }

    } // end nface==5



    // ******************************************************************************
    //      k=kdim boundary          periodic boundary                   bctype 2005
    // ******************************************************************************

    if (nface == 6) {

        // Load qp values into qrotk array:
        for (l = 1; l <= 5; l++)
        for (i = ista; i <= iend1; i++)
        for (j = jsta; j <= jend1; j++) {
            qrotk(j,1,i,l) = qp(j,1,i,l);
            kmx = std::min(kdimp-1, 2);
            qrotk(j,2,i,l) = qp(j,kmx,i,l);
        }

        // Rotate qrotk values
        {
            int jd_p = jdimp, kd2 = 2, id_p = idimp;
            if (kdimp == 2) {
                int is1=ista, ie1=iend1, js1=jsta, je1=jend1, ks1=1, ke1=1;
                rotateq_ns::rotateq(jd_p, kd2, id_p, qrotk, qrotk,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
                double dthtx2 = 2.0*dthtx, dthty2 = 2.0*dthty, dthtz2 = 2.0*dthtz;
                int ks2=2, ke2=2;
                rotateq_ns::rotateq(jd_p, kd2, id_p, qrotk, qrotk,
                                    is1, ie1, js1, je1, ks2, ke2,
                                    dthtx2, dthty2, dthtz2);
            } else {
                int is1=ista, ie1=iend1, js1=jsta, je1=jend1, ks1=1, ke1=2;
                rotateq_ns::rotateq(jd_p, kd2, id_p, qrotk, qrotk,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
            }
        }

        // Apply periodic BCs
        for (l = 1; l <= 5; l++)
        for (i = ista; i <= iend1; i++)
        for (j = jsta; j <= jend1; j++) {
            qk0(j,i,l,3) = qrotk(j,1,i,l);
            qk0(j,i,l,4) = qrotk(j,2,i,l);
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                vk0(j,i,1,3) = vp(j,1,i);
                vk0(j,i,1,4) = 0.0;
            }
        }

        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                // obtain the rotation matrix and its transpose
                float fthx = (float)dthtx, fthy = (float)dthty, fthz = (float)dthtz;
                FortranArray2DRef<float> rn_ref = rn.ref(), rnt_ref = rnt.ref();
                module_stm_2005_ns::stm2k5_get_rotmat(fthx, fthy, fthz, rn_ref, rnt_ref);

                // rotate the stress tensors
                for (i = ista; i <= iend1; i++) {
                    for (j = jsta; j <= jend1; j++) {
                        {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(j,1,i,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tk0(j,i,ll,3) = (double)tmp_out[ll-1];
                        }
                        if (kdimp == 2) {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tk0(j,i,ll,3);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tk0(j,i,ll,4) = (double)tmp_out[ll-1];
                        } else {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(j,2,i,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) tk0(j,i,ll,4) = (double)tmp_out[ll-1];
                        }
                    }
                }

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                kdp2 = std::min(2, kdimp-1);
                for (l = 1; l <= nummem; l++)
                for (i = ista; i <= iend1; i++)
                for (j = jsta; j <= jend1; j++) {
                    tk0(j,i,l,3) = tp(j,1,i,l);
                    tk0(j,i,l,4) = tp(j,kdp2,i,l);
                }
            }
        }

    } // end nface==6



    // ******************************************************************************
    //      i=1 boundary             periodic boundary                   bctype 2005
    // ******************************************************************************

    if (nface == 1) {

        // Load qp values into qroti array:
        for (l = 1; l <= 5; l++)
        for (k = ksta; k <= kend1; k++)
        for (j = jsta; j <= jend1; j++) {
            qroti(j,k,1,l) = qp(j,k,idimp-1,l);
            imx = std::max(idimp-2, 1);
            qroti(j,k,2,l) = qp(j,k,imx,l);
        }

        // Rotate qroti values
        {
            int jd_p = jdimp, kd_p = kdimp, id2 = 2;
            if (idimp == 2) {
                int is1=1, ie1=1, js1=jsta, je1=jend1, ks1=ksta, ke1=kend1;
                rotateq_ns::rotateq(jd_p, kd_p, id2, qroti, qroti,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
                double dthtx2 = 2.0*dthtx, dthty2 = 2.0*dthty, dthtz2 = 2.0*dthtz;
                int is2=2, ie2=2;
                rotateq_ns::rotateq(jd_p, kd_p, id2, qroti, qroti,
                                    is2, ie2, js1, je1, ks1, ke1,
                                    dthtx2, dthty2, dthtz2);
            } else {
                int is1=1, ie1=2, js1=jsta, je1=jend1, ks1=ksta, ke1=kend1;
                rotateq_ns::rotateq(jd_p, kd_p, id2, qroti, qroti,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
            }
        }

        // Apply periodic BCs
        for (l = 1; l <= 5; l++)
        for (k = ksta; k <= kend1; k++)
        for (j = jsta; j <= jend1; j++) {
            qi0(j,k,l,1) = qroti(j,k,1,l);
            qi0(j,k,l,2) = qroti(j,k,2,l);
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                vi0(j,k,1,1) = vp(j,k,idimp-1);
                vi0(j,k,1,2) = 0.0;
            }
        }

        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                // obtain the rotation matrix and its transpose
                float fthx = (float)dthtx, fthy = (float)dthty, fthz = (float)dthtz;
                FortranArray2DRef<float> rn_ref = rn.ref(), rnt_ref = rnt.ref();
                module_stm_2005_ns::stm2k5_get_rotmat(fthx, fthy, fthz, rn_ref, rnt_ref);

                // rotate the stress tensors
                for (k = ksta; k <= kend1; k++) {
                    for (j = jsta; j <= jend1; j++) {
                        {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(j,k,idimp-1,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) ti0(j,k,ll,1) = (double)tmp_out[ll-1];
                        }
                        if (idimp == 2) {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)ti0(j,k,ll,1);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) ti0(j,k,ll,2) = (double)tmp_out[ll-1];
                        } else {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(j,k,idimp-2,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) ti0(j,k,ll,2) = (double)tmp_out[ll-1];
                        }
                    }
                }

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                idm2 = std::max(idimp-2, 1);
                for (l = 1; l <= nummem; l++)
                for (k = ksta; k <= kend1; k++)
                for (j = jsta; j <= jend1; j++) {
                    ti0(j,k,l,1) = tp(j,k,idimp-1,l);
                    ti0(j,k,l,2) = tp(j,k,idm2,l);
                }
            }
        }

    } // end nface==1



    // ******************************************************************************
    //      i=idim boundary          periodic boundary                   bctype 2005
    // ******************************************************************************

    if (nface == 2) {

        // Load qp values into qroti array:
        for (l = 1; l <= 5; l++)
        for (k = ksta; k <= kend1; k++)
        for (j = jsta; j <= jend1; j++) {
            qroti(j,k,1,l) = qp(j,k,1,l);
            imx = std::min(idimp-1, 2);
            qroti(j,k,2,l) = qp(j,k,imx,l);
        }

        // Rotate qroti values
        {
            int jd_p = jdimp, kd_p = kdimp, id2 = 2;
            if (idimp == 2) {
                int is1=1, ie1=1, js1=jsta, je1=jend1, ks1=ksta, ke1=kend1;
                rotateq_ns::rotateq(jd_p, kd_p, id2, qroti, qroti,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
                double dthtx2 = 2.0*dthtx, dthty2 = 2.0*dthty, dthtz2 = 2.0*dthtz;
                int is2=2, ie2=2;
                rotateq_ns::rotateq(jd_p, kd_p, id2, qroti, qroti,
                                    is2, ie2, js1, je1, ks1, ke1,
                                    dthtx2, dthty2, dthtz2);
            } else {
                int is1=1, ie1=2, js1=jsta, je1=jend1, ks1=ksta, ke1=kend1;
                rotateq_ns::rotateq(jd_p, kd_p, id2, qroti, qroti,
                                    is1, ie1, js1, je1, ks1, ke1,
                                    dthtx, dthty, dthtz);
            }
        }

        // Apply periodic BCs
        for (l = 1; l <= 5; l++)
        for (k = ksta; k <= kend1; k++)
        for (j = jsta; j <= jend1; j++) {
            qi0(j,k,l,3) = qroti(j,k,1,l);
            qi0(j,k,l,4) = qroti(j,k,2,l);
        }

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                vi0(j,k,1,3) = vp(j,k,1);
                vi0(j,k,1,4) = 0.0;
            }
        }

        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                // obtain the rotation matrix and its transpose
                float fthx = (float)dthtx, fthy = (float)dthty, fthz = (float)dthtz;
                FortranArray2DRef<float> rn_ref = rn.ref(), rnt_ref = rnt.ref();
                module_stm_2005_ns::stm2k5_get_rotmat(fthx, fthy, fthz, rn_ref, rnt_ref);

                // rotate the stress tensors
                for (k = ksta; k <= kend1; k++) {
                    for (j = jsta; j <= jend1; j++) {
                        {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(j,k,1,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) ti0(j,k,ll,3) = (double)tmp_out[ll-1];
                        }
                        if (idimp == 2) {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)ti0(j,k,ll,3);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) ti0(j,k,ll,4) = (double)tmp_out[ll-1];
                        } else {
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)tp(j,k,2,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) ti0(j,k,ll,4) = (double)tmp_out[ll-1];
                        }
                    }
                }

            } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                idp2 = std::min(2, idimp-1);
                for (l = 1; l <= nummem; l++)
                for (k = ksta; k <= kend1; k++)
                for (j = jsta; j <= jend1; j++) {
                    ti0(j,k,l,3) = tp(j,k,1,l);
                    ti0(j,k,l,4) = tp(j,k,idp2,l);
                }
            }
        }

    } // end nface==2

} // end bc2005()

} // namespace bc2005_ns
