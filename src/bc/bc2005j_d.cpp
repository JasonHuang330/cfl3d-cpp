// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2005j_d.h"
#include "bc2005.h"
#include "rotateq.h"
#include "module_stm_2005.h"
#include <cstring>
#include <algorithm>

// Forward declarations for stm2k5 functions (module_stm_2005_ns header is empty)
namespace module_stm_2005_ns {
    void stm2k5_get_rotmat(float& thx, float& thy, float& thz,
                           FortranArray2DRef<float> rn, FortranArray2DRef<float> rnt);
    void stm2k5_bc(FortranArray1DRef<float> vin, FortranArray2DRef<float> rn,
                   FortranArray2DRef<float> rnt, FortranArray1DRef<float> vout);
}

namespace bc2005j_d_ns {

// bc() delegates to bc2005_ns::bc()
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

// bc2005j_d: periodic BC for j-faces (nface==3 or nface==4)
void bc2005j_d(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> dataj0,
               int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
               int& nface, int& mdim, int& ndim, FortranArray4DRef<double> bcdata,
               char (&filname)[80], int& jdimp, int& kdimp, int& idimp,
               FortranArray4DRef<double> datapj, int& nbl, int& nblp, int& ldim,
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

    // Local rotation matrix (Fortran real = float), wrapped as FortranArray2DRef
    float rn_data[3][3], rnt_data[3][3];
    FortranArray2DRef<float> rn(reinterpret_cast<float*>(rn_data), 3, 3);
    FortranArray2DRef<float> rnt(reinterpret_cast<float*>(rnt_data), 3, 3);

    // Derived dimensions
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    // Extract rotation angles from bcdata
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

    // **************************************************************************
    //  j=1 boundary             periodic boundary                   bctype 2005
    // **************************************************************************

    if (nface == 3) {

        if (iflag == 0) {
            // rotate datapj values before setting dataj0 values if data
            // corresponds to q
            int jd2 = 2;
            if (jdimp == 2) {
                int jsta1 = 1, jend1_r = 1;
                rotateq_ns::rotateq(jd2, kdimp, idimp, datapj, datapj,
                                    ista, iend1, jsta1, jend1_r, ksta, kend1,
                                    dthtx, dthty, dthtz);
                int jsta2 = 2, jend2 = 2;
                double dthtx2 = 2.0*dthtx, dthty2 = 2.0*dthty, dthtz2 = 2.0*dthtz;
                rotateq_ns::rotateq(jd2, kdimp, idimp, datapj, datapj,
                                    ista, iend1, jsta2, jend2, ksta, kend1,
                                    dthtx2, dthty2, dthtz2);
            } else {
                int jsta1 = 1, jend2 = 2;
                rotateq_ns::rotateq(jd2, kdimp, idimp, datapj, datapj,
                                    ista, iend1, jsta1, jend2, ksta, kend1,
                                    dthtx, dthty, dthtz);
            }

            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    for (int l = 1; l <= 5; l++) {
                        dataj0(k,i,l,1) = datapj(1,k,i,l);
                        dataj0(k,i,l,2) = datapj(2,k,i,l);
                    }
                }
            }
        } // iflag==0

        if (iflag == 1) {
            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (int i = ista; i <= iend1; i++) {
                    for (int k = ksta; k <= kend1; k++) {
                        dataj0(k,i,1,1) = datapj(1,k,i,1);
                        dataj0(k,i,1,2) = 0.0;
                    }
                }
            }
        } // iflag==1

        if (iflag == 2) {
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {
                    float fthx = (float)dthtx, fthy = (float)dthty, fthz = (float)dthtz;
                    module_stm_2005_ns::stm2k5_get_rotmat(fthx, fthy, fthz, rn, rnt);
                    for (int i = ista; i <= iend1; i++) {
                        for (int k = ksta; k <= kend1; k++) {
                            // datapj(1,k,i,1:7) — varying 4th index, non-contiguous
                            // dataj0(k,i,1:7,1) — varying 3rd index, non-contiguous
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)datapj(1,k,i,ll);
                            FortranArray1DRef<float> vin(tmp_in, 7);
                            FortranArray1DRef<float> vout(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(vin, rn, rnt, vout);
                            for (int ll = 1; ll <= 7; ll++) dataj0(k,i,ll,1) = tmp_out[ll-1];

                            if (jdimp == 2) {
                                // dataj0(k,i,1:7,1) → dataj0(k,i,1:7,2)
                                float tmp_in2[7], tmp_out2[7];
                                for (int ll = 1; ll <= 7; ll++) tmp_in2[ll-1] = (float)dataj0(k,i,ll,1);
                                FortranArray1DRef<float> vin2(tmp_in2, 7);
                                FortranArray1DRef<float> vout2(tmp_out2, 7);
                                module_stm_2005_ns::stm2k5_bc(vin2, rn, rnt, vout2);
                                for (int ll = 1; ll <= 7; ll++) dataj0(k,i,ll,2) = tmp_out2[ll-1];
                            } else {
                                // datapj(2,k,i,1:7) → dataj0(k,i,1:7,2)
                                float tmp_in3[7], tmp_out3[7];
                                for (int ll = 1; ll <= 7; ll++) tmp_in3[ll-1] = (float)datapj(2,k,i,ll);
                                FortranArray1DRef<float> vin3(tmp_in3, 7);
                                FortranArray1DRef<float> vout3(tmp_out3, 7);
                                module_stm_2005_ns::stm2k5_bc(vin3, rn, rnt, vout3);
                                for (int ll = 1; ll <= 7; ll++) dataj0(k,i,ll,2) = tmp_out3[ll-1];
                            }
                        }
                    }
                } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (int i = ista; i <= iend1; i++) {
                        for (int k = ksta; k <= kend1; k++) {
                            for (int l = 1; l <= ldim; l++) {
                                dataj0(k,i,l,1) = datapj(1,k,i,l);
                                dataj0(k,i,l,2) = datapj(2,k,i,l);
                            }
                        }
                    }
                }
            }
        } // iflag==2

    } // nface==3

    // **************************************************************************
    //  j=jdim boundary          periodic boundary                   bctype 2005
    // **************************************************************************

    if (nface == 4) {

        if (iflag == 0) {
            // rotate datapk values before setting dataj0 values if data
            // corresponds to q
            int jd2 = 2;
            if (jdimp == 2) {
                int jsta1 = 1, jend1_r = 1;
                rotateq_ns::rotateq(jd2, kdimp, idimp, datapj, datapj,
                                    ista, iend1, jsta1, jend1_r, ksta, kend1,
                                    dthtx, dthty, dthtz);
                int jsta2 = 2, jend2 = 2;
                double dthtx2 = 2.0*dthtx, dthty2 = 2.0*dthty, dthtz2 = 2.0*dthtz;
                rotateq_ns::rotateq(jd2, kdimp, idimp, datapj, datapj,
                                    ista, iend1, jsta2, jend2, ksta, kend1,
                                    dthtx2, dthty2, dthtz2);
            } else {
                int jsta1 = 1, jend2 = 2;
                rotateq_ns::rotateq(jd2, kdimp, idimp, datapj, datapj,
                                    ista, iend1, jsta1, jend2, ksta, kend1,
                                    dthtx, dthty, dthtz);
            }

            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    for (int l = 1; l <= 5; l++) {
                        dataj0(k,i,l,3) = datapj(1,k,i,l);
                        dataj0(k,i,l,4) = datapj(2,k,i,l);
                    }
                }
            }
        } // iflag==0

        if (iflag == 1) {
            if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
                for (int i = ista; i <= iend1; i++) {
                    for (int k = ksta; k <= kend1; k++) {
                        dataj0(k,i,1,3) = datapj(1,k,i,1);
                        dataj0(k,i,1,4) = 0.0;
                    }
                }
            }
        } // iflag==1

        if (iflag == 2) {
            // only need to do advanced model turbulence B.C.s on finest grid
            if (level >= lglobal) {
                if (ivisc[2] >= 70 || ivisc[1] >= 70 || ivisc[0] >= 70) {
                    float fthx = (float)dthtx, fthy = (float)dthty, fthz = (float)dthtz;
                    module_stm_2005_ns::stm2k5_get_rotmat(fthx, fthy, fthz, rn, rnt);
                    for (int i = ista; i <= iend1; i++) {
                        for (int k = ksta; k <= kend1; k++) {
                            // datapj(1,k,i,1:7) → dataj0(k,i,1:7,3)
                            float tmp_in[7], tmp_out[7];
                            for (int ll = 1; ll <= 7; ll++) tmp_in[ll-1] = (float)datapj(1,k,i,ll);
                            FortranArray1DRef<float> vin(tmp_in, 7);
                            FortranArray1DRef<float> vout(tmp_out, 7);
                            module_stm_2005_ns::stm2k5_bc(vin, rn, rnt, vout);
                            for (int ll = 1; ll <= 7; ll++) dataj0(k,i,ll,3) = tmp_out[ll-1];

                            if (jdimp == 2) {
                                // dataj0(k,i,1:7,3) → dataj0(k,i,1:7,4)
                                float tmp_in2[7], tmp_out2[7];
                                for (int ll = 1; ll <= 7; ll++) tmp_in2[ll-1] = (float)dataj0(k,i,ll,3);
                                FortranArray1DRef<float> vin2(tmp_in2, 7);
                                FortranArray1DRef<float> vout2(tmp_out2, 7);
                                module_stm_2005_ns::stm2k5_bc(vin2, rn, rnt, vout2);
                                for (int ll = 1; ll <= 7; ll++) dataj0(k,i,ll,4) = tmp_out2[ll-1];
                            } else {
                                // datapj(2,k,i,1:7) → dataj0(k,i,1:7,4)
                                float tmp_in3[7], tmp_out3[7];
                                for (int ll = 1; ll <= 7; ll++) tmp_in3[ll-1] = (float)datapj(2,k,i,ll);
                                FortranArray1DRef<float> vin3(tmp_in3, 7);
                                FortranArray1DRef<float> vout3(tmp_out3, 7);
                                module_stm_2005_ns::stm2k5_bc(vin3, rn, rnt, vout3);
                                for (int ll = 1; ll <= 7; ll++) dataj0(k,i,ll,4) = tmp_out3[ll-1];
                            }
                        }
                    }
                } else if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                    for (int i = ista; i <= iend1; i++) {
                        for (int k = ksta; k <= kend1; k++) {
                            for (int l = 1; l <= ldim; l++) {
                                dataj0(k,i,l,3) = datapj(1,k,i,l);
                                dataj0(k,i,l,4) = datapj(2,k,i,l);
                            }
                        }
                    }
                }
            }
        } // iflag==2

    } // nface==4

    return;
}

} // namespace bc2005j_d_ns
