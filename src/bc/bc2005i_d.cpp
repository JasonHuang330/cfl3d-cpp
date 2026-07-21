// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2005i_d.h"
#include "bc2005.h"
#include "rotateq.h"
#include "module_stm_2005.h"
#include "runtime/common_blocks.h"
#include <cstring>
#include <algorithm>

// Forward declarations for module_stm_2005 functions (header is empty but functions exist in module_stm_2005.cpp)
namespace module_stm_2005_ns {
    void stm2k5_get_rotmat(float& thx, float& thy, float& thz,
                           FortranArray2DRef<float> rn,
                           FortranArray2DRef<float> rnt);
    void stm2k5_bc(FortranArray1DRef<float> vin,
                   FortranArray2DRef<float> rn,
                   FortranArray2DRef<float> rnt,
                   FortranArray1DRef<float> vout);
}

namespace bc2005i_d_ns {

// bc() delegates to bc2005_ns::bc() — same pattern as bc2005j_d and bc2005k_d
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
    bc2005_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, cl, nou, bou, nbuf, ibufdim,
                  maxbl, maxgr, maxseg, itrans, irotat, idefrm, igridg, nblg,
                  nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
                  ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek,
                  lwdat, myid, idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
}

// bc2005() delegates to bc2005_ns::bc2005()
void bc2005(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0, int& ista, int& iend, int& jsta, int& jend,
            int& ksta, int& kend, int& nface, FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0, FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0, int& mdim, int& ndim,
            FortranArray4DRef<double> bcdata, char (&filname)[80],
            FortranArray4DRef<double> qp, FortranArray3DRef<double> vp,
            FortranArray4DRef<double> tp, int& jdimp, int& kdimp, int& idimp,
            FortranArray4DRef<double> qrotj, FortranArray4DRef<double> qrotk,
            FortranArray4DRef<double> qroti, int& nbl, int& nblp,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd,
            int& maxbl, int& nummem)
{
    bc2005_ns::bc2005(jdim, kdim, idim, q, qj0, qk0, qi0, ista, iend, jsta, jend,
                      ksta, kend, nface, tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                      mdim, ndim, bcdata, filname, qp, vp, tp, jdimp, kdimp, idimp,
                      qrotj, qrotk, qroti, nbl, nblp, nou, bou, nbuf, ibufdim,
                      myid, mblk2nd, maxbl, nummem);
}

// bc2005i_d: Set periodic boundary conditions for i-faces (nface=1 or nface=2)
void bc2005i_d(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> datai0,
               int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
               int& nface, int& mdim, int& ndim, FortranArray4DRef<double> bcdata,
               char (&filname)[80], int& jdimp, int& kdimp, int& idimp,
               FortranArray4DRef<double> datapi, int& nbl, int& nblp, int& ldim,
               FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
               int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd,
               int& maxbl, int& iflag)
{
    // COMMON block references
    int&   level   = cmn_mgrd.level;
    int&   lglobal = cmn_mgrd.lglobal;
    int*   ivisc   = cmn_reyue.ivisc;   // 0-based: ivisc[0]=ivisc(1), ivisc[1]=ivisc(2), ivisc[2]=ivisc(3)
    float& radtodeg = cmn_conversion.radtodeg;
    int&   ialph   = cmn_igrdtyp.ialph;

    // Local variables
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    // This bc makes use of only one plane of data
    int ip = 1;
    double dthtx = (double)bcdata(1,1,ip,2) / (double)radtodeg;
    double dthty, dthtz;
    if (ialph == 0) {
        dthty = (double)bcdata(1,1,ip,3) / (double)radtodeg;
        dthtz = (double)bcdata(1,1,ip,4) / (double)radtodeg;
    } else {
        dthty = -(double)bcdata(1,1,ip,4) / (double)radtodeg;
        dthtz =  (double)bcdata(1,1,ip,3) / (double)radtodeg;
    }

    // Rotation matrix (float, as Fortran real) — 1-based 3x3
    FortranArray2D<float> rn(1,3,1,3), rnt(1,3,1,3);

    //******************************************************************************
    //  i=1 boundary             periodic boundary                   bctype 2005
    //******************************************************************************

    if (nface == 1) {

        if (iflag == 0) {
            // rotate datapi values before setting datai0 values if data
            // corresponds to q
            if (idimp == 2) {
                int ip1 = 1, ip2 = 1;
                rotateq_ns::rotateq(jdimp, kdimp, idimp, datapi, datapi,
                                    ip1, ip2, jsta, jend1, ksta, kend1,
                                    dthtx, dthty, dthtz);
                ip1 = 2; ip2 = 2;
                double dthtx2 = 2.0 * dthtx;
                double dthty2 = 2.0 * dthty;
                double dthtz2 = 2.0 * dthtz;
                rotateq_ns::rotateq(jdimp, kdimp, idimp, datapi, datapi,
                                    ip1, ip2, jsta, jend1, ksta, kend1,
                                    dthtx2, dthty2, dthtz2);
            } else {
                int ip1 = 1, ip2 = 2;
                rotateq_ns::rotateq(jdimp, kdimp, idimp, datapi, datapi,
                                    ip1, ip2, jsta, jend1, ksta, kend1,
                                    dthtx, dthty, dthtz);
            }

            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    for (int l = 1; l <= 5; l++) {
                        datai0(j,k,l,1) = datapi(j,k,1,l);
                        datai0(j,k,l,2) = datapi(j,k,2,l);
                    }
                }
            }

        } // end iflag==0

        if (iflag == 1) {

            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (int k = ksta; k <= kend1; k++) {
                    for (int j = jsta; j <= jend1; j++) {
                        datai0(j,k,1,1) = datapi(j,k,1,1);
                        datai0(j,k,1,2) = 0.0;
                    }
                }
            }

        } // end iflag==1

        if (iflag == 2) {

            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                    float fthx = (float)dthtx, fthy = (float)dthty, fthz = (float)dthtz;
                    FortranArray2DRef<float> rn_ref = rn.ref(), rnt_ref = rnt.ref();
                    module_stm_2005_ns::stm2k5_get_rotmat(fthx, fthy, fthz, rn_ref, rnt_ref);
                    for (int k = ksta; k <= kend1; k++) {
                        for (int j = jsta; j <= jend1; j++) {
                            // datapi(j,k,1,1:7) → non-contiguous, use temp
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)datapi(j,k,1,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) datai0(j,k,ll,1) = (double)tmp_out[ll-1];

                            if (idimp == 2) {
                                // datai0(j,k,1:7,1) → non-contiguous, use temp
                                float tmp_in2[7], tmp_out2[7];
                                for (int ll = 1; ll <= 7; ll++) tmp_in2[ll-1] = (float)datai0(j,k,ll,1);
                                FortranArray1DRef<float> in2_ref(tmp_in2, 7);
                                FortranArray1DRef<float> out2_ref(tmp_out2, 7);
                                module_stm_2005_ns::stm2k5_bc(in2_ref, rn_ref, rnt_ref, out2_ref);
                                for (int ll = 1; ll <= 7; ll++) datai0(j,k,ll,2) = (double)tmp_out2[ll-1];
                            } else {
                                // datapi(j,k,2,1:7) → non-contiguous, use temp
                                float tmp_in3[7], tmp_out3[7];
                                for (int ll = 1; ll <= 7; ll++) tmp_in3[ll-1] = (float)datapi(j,k,2,ll);
                                FortranArray1DRef<float> in3_ref(tmp_in3, 7);
                                FortranArray1DRef<float> out3_ref(tmp_out3, 7);
                                module_stm_2005_ns::stm2k5_bc(in3_ref, rn_ref, rnt_ref, out3_ref);
                                for (int ll = 1; ll <= 7; ll++) datai0(j,k,ll,2) = (double)tmp_out3[ll-1];
                            }
                        }
                    }

                } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {

                    for (int k = ksta; k <= kend1; k++) {
                        for (int j = jsta; j <= jend1; j++) {
                            for (int l = 1; l <= ldim; l++) {
                                datai0(j,k,l,1) = datapi(j,k,1,l);
                                datai0(j,k,l,2) = datapi(j,k,2,l);
                            }
                        }
                    }

                }
            }

        } // end iflag==2

    } // end nface==1

    //******************************************************************************
    //  i=idim boundary          periodic boundary                   bctype 2005
    //******************************************************************************

    if (nface == 2) {

        if (iflag == 0) {
            // rotate datapi values before setting datai0 values if data
            // corresponds to q
            if (idimp == 2) {
                int ip1 = 1, ip2 = 1;
                rotateq_ns::rotateq(jdimp, kdimp, idimp, datapi, datapi,
                                    ip1, ip2, jsta, jend1, ksta, kend1,
                                    dthtx, dthty, dthtz);
                ip1 = 2; ip2 = 2;
                double dthtx2 = 2.0 * dthtx;
                double dthty2 = 2.0 * dthty;
                double dthtz2 = 2.0 * dthtz;
                rotateq_ns::rotateq(jdimp, kdimp, idimp, datapi, datapi,
                                    ip1, ip2, jsta, jend1, ksta, kend1,
                                    dthtx2, dthty2, dthtz2);
            } else {
                int ip1 = 1, ip2 = 2;
                rotateq_ns::rotateq(jdimp, kdimp, idimp, datapi, datapi,
                                    ip1, ip2, jsta, jend1, ksta, kend1,
                                    dthtx, dthty, dthtz);
            }

            for (int l = 1; l <= 5; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    for (int j = jsta; j <= jend1; j++) {
                        datai0(j,k,l,3) = datapi(j,k,1,l);
                        datai0(j,k,l,4) = datapi(j,k,2,l);
                    }
                }
            }

        } // end iflag==0

        if (iflag == 1) {

            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int k = ksta; k <= kend1; k++) {
                    for (int j = jsta; j <= jend1; j++) {
                        datai0(j,k,1,3) = datapi(j,k,1,1);
                        datai0(j,k,1,4) = 0.0;
                    }
                }
            }

        } // end iflag==1

        if (iflag == 2) {

            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {

                    float fthx = (float)dthtx, fthy = (float)dthty, fthz = (float)dthtz;
                    FortranArray2DRef<float> rn_ref = rn.ref(), rnt_ref = rnt.ref();
                    module_stm_2005_ns::stm2k5_get_rotmat(fthx, fthy, fthz, rn_ref, rnt_ref);
                    for (int k = ksta; k <= kend1; k++) {
                        for (int j = jsta; j <= jend1; j++) {
                            // datapi(j,k,1,1:7) → non-contiguous, use temp
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)datapi(j,k,1,ll);
                            FortranArray1DRef<float> in_ref(tmp_in, 7);
                            FortranArray1DRef<float> out_ref(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(in_ref, rn_ref, rnt_ref, out_ref);
                            for (int ll = 1; ll <= 7; ll++) datai0(j,k,ll,3) = (double)tmp_out[ll-1];

                            if (idimp == 2) {
                                // datai0(j,k,1:7,3) → non-contiguous, use temp
                                float tmp_in2[7], tmp_out2[7];
                                for (int ll = 1; ll <= 7; ll++) tmp_in2[ll-1] = (float)datai0(j,k,ll,3);
                                FortranArray1DRef<float> in2_ref(tmp_in2, 7);
                                FortranArray1DRef<float> out2_ref(tmp_out2, 7);
                                module_stm_2005_ns::stm2k5_bc(in2_ref, rn_ref, rnt_ref, out2_ref);
                                for (int ll = 1; ll <= 7; ll++) datai0(j,k,ll,4) = (double)tmp_out2[ll-1];
                            } else {
                                // datapi(j,k,2,1:7) → non-contiguous, use temp
                                float tmp_in3[7], tmp_out3[7];
                                for (int ll = 1; ll <= 7; ll++) tmp_in3[ll-1] = (float)datapi(j,k,2,ll);
                                FortranArray1DRef<float> in3_ref(tmp_in3, 7);
                                FortranArray1DRef<float> out3_ref(tmp_out3, 7);
                                module_stm_2005_ns::stm2k5_bc(in3_ref, rn_ref, rnt_ref, out3_ref);
                                for (int ll = 1; ll <= 7; ll++) datai0(j,k,ll,4) = (double)tmp_out3[ll-1];
                            }
                        }
                    }

                } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {

                    for (int k = ksta; k <= kend1; k++) {
                        for (int j = jsta; j <= jend1; j++) {
                            for (int l = 1; l <= ldim; l++) {
                                datai0(j,k,l,3) = datapi(j,k,1,l);
                                datai0(j,k,l,4) = datapi(j,k,2,l);
                            }
                        }
                    }

                }
            }

        } // end iflag==2

    } // end nface==2

    return;
}

} // namespace bc2005i_d_ns
