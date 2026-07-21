// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rrestg.h"
#include "runtime/fortran_io.h"
#include "termn8.h"
#include <cmath>
#include <cstdio>
#include <cstring>

namespace rrestg_ns {

// Helper: read Fortran unformatted sequential record begin/end markers
static void read_rec_begin(FILE* fp) {
    int32_t marker;
    fread(&marker, sizeof(int32_t), 1, fp);
}
static void read_rec_end(FILE* fp) {
    int32_t marker;
    fread(&marker, sizeof(int32_t), 1, fp);
}

void rrestg(int& nbl, int& igrid, int& jdim, int& kdim, int& idim,
            FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
            FortranArray3DRef<double> xnm2, FortranArray3DRef<double> ynm2, FortranArray3DRef<double> znm2,
            FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, FortranArray4DRef<double> delti,
            FortranArray4DRef<double> qc0,
            int& nflagg, int& iuns,
            FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans,
            FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz,
            FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig,
            FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx,
            FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx,
            FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz,
            FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr,
            FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0,
            FortranArray1DRef<double> time2,
            FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl,
            FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm,
            FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae,
            FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae,
            FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae,
            FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae,
            FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae,
            FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
            FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
            FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
            FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass,
            FortranArray2DRef<double> damp, FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0,
            int& nmds, int& maxaes,
            FortranArray2DRef<double> aesrfdat,
            FortranArray3DRef<double> perturb,
            int& myhost, int& myid, int& mycomm,
            FortranArray1DRef<int> mblk2nd, int& maxbl,
            int& ibufdim, int& nbuf,
            FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou,
            FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg,
            FortranArray2DRef<int> iaesurf, int& maxsegdg,
            FortranArray1DRef<double> wk, int& nwork,
            int& idima, int& jdima, int& kdima,
            FortranArray4DRef<double> tursav2, int& nummem)
{
    // COMMON block aliases
    float& dtold    = cmn_cfl.dtold;
    float& dt       = cmn_info.dt;
    float* qiv      = cmn_ivals.qiv;   // 1-based: qiv[0..4]
    float* tur10    = cmn_ivals.tur10; // 1-based: tur10[0..6]
    int32_t& ivmx   = cmn_maxiv.ivmx;
    float& xmc      = cmn_fsum.xmc;
    float& ymc      = cmn_fsum.ymc;
    float& zmc      = cmn_fsum.zmc;
    float& omegaxmc  = cmn_motionmc.omegaxmc;
    float& omegaymc  = cmn_motionmc.omegaymc;
    float& omegazmc  = cmn_motionmc.omegazmc;
    float& xorigmc   = cmn_motionmc.xorigmc;
    float& yorigmc   = cmn_motionmc.yorigmc;
    float& zorigmc   = cmn_motionmc.zorigmc;
    float& thetaxmc  = cmn_motionmc.thetaxmc;
    float& thetaymc  = cmn_motionmc.thetaymc;
    float& thetazmc  = cmn_motionmc.thetazmc;
    float& time2mc   = cmn_motionmc.time2mc;
    int32_t& ita     = cmn_unst.ita;
    int32_t& iunst   = cmn_unst.iunst;
    int32_t& icgns   = cmn_cgns.icgns;
    int32_t& itime2read = cmn_turbconv.itime2read;
    int32_t& ndgrd   = cmn_deformz.ndgrd;
    int32_t& irotatmc = cmn_motionmc.irotatmc;

    // Local variables
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;
    int jdum, kdum, idum;
    int nummem_read;
    double dum;
    double xdum, ydum, zdum;
    int itrans1, irotat1;
    double rfreqt1, rfreqr1;
    double xorig1, yorig1, zorig1;
    double xorig01, yorig01, zorig01;
    double utrans1, vtrans1, wtrans1;
    double dxmx1, dymx1, dzmx1;
    int itransmc1, irotatmc1;
    double rfreqtmc1, rfreqrmc1;
    double xorigmc1, yorigmc1, zorigmc1;
    double xorig0mc1, yorig0mc1, zorig0mc1;
    double utransmc1, vtransmc1, wtransmc1;
    double xmc1, ymc1, zmc1;
    double dxmxmc1, dymxmc1, dzmxmc1;
    double thetax1, thetay1, thetaz1;
    double omegax1, omegay1, omegaz1;
    double dthxmx1, dthymx1, dthzmx1;
    double thetaxmc1, thetaymc1, thetazmc1;
    double omegaxmc1, omegaymc1, omegazmc1;
    double dthxmxmc1, dthymxmc1, dthzmxmc1;
    double time21, time2mc1;
    double time211, time2mc11;
    double dt1;
    int idfrmsdum;


    if (nflagg == 0) {
        // *******************************
        //       read 2nd order time data
        // *******************************

        if (myid == myhost) {
            if (icgns != 1) {
                FILE* fp2  = fortran_get_unit(2);
                FILE* fp11 = fortran_get_unit(11);

                // write(11,50) nbl
                // 50 format(33h reading 2nd order time data from,23h restart file for block,i4)
                fprintf(fp11, " reading 2nd order time data from restart file for block%4d\n", nbl);

                // read(2) jdum,kdum,idum
                read_rec_begin(fp2);
                fread(&jdum, sizeof(int32_t), 1, fp2);
                fread(&kdum, sizeof(int32_t), 1, fp2);
                fread(&idum, sizeof(int32_t), 1, fp2);
                read_rec_end(fp2);

                if (jdum != jdim && kdum != kdim && idum != idim) {
                    // write(11,55) nbl
                    // 55 format(43h error in 2nd order time section of restart,15h file for block,i3)
                    fprintf(fp11, " error in 2nd order time section of restart file for block%3d\n", nbl);
                    // write(11,65) jdum,kdum,idum,jdim,kdim,idim
                    // 65 format(29h restart file has j x k x i :,3i5,/,29h should be        j x k x i :,3i5)
                    fprintf(fp11, " restart file has j x k x i :%5d%5d%5d\n", jdum, kdum, idum);
                    fprintf(fp11, " should be        j x k x i :%5d%5d%5d\n", jdim, kdim, idim);
                    int ierrflg = -1;
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }

                // read(2) ((((qc0(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,5)
                {
                    read_rec_begin(fp2);
                    for (int l = 1; l <= 5; l++)
                        for (int i = 1; i <= idim1; i++)
                            for (int k = 1; k <= kdim1; k++)
                                for (int j = 1; j <= jdim1; j++) {
                                    double v;
                                    fread(&v, sizeof(double), 1, fp2);
                                    qc0(j, k, i, l) = v;
                                }
                    read_rec_end(fp2);
                }

                // itime2read=1=default=YES, read 2nd order time-accurate turb info
                // itime2read=0=NO, do not read
                if (itime2read != 0) {
                    // read(2) dtold
                    {
                        read_rec_begin(fp2);
                        float dtold_val;
                        fread(&dtold_val, sizeof(float), 1, fp2);
                        dtold = dtold_val;
                        read_rec_end(fp2);
                    }

                    if (ivmx >= 4) {
                        // read(2) ((((tursav2(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,2)
                        read_rec_begin(fp2);
                        for (int l = 1; l <= 2; l++)
                            for (int i = 1; i <= idim1; i++)
                                for (int k = 1; k <= kdim1; k++)
                                    for (int j = 1; j <= jdim1; j++) {
                                        double v;
                                        fread(&v, sizeof(double), 1, fp2);
                                        tursav2(j, k, i, l) = v;
                                    }
                        read_rec_end(fp2);
                    } else {
                        // read(2) ((((dum,j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,2)
                        read_rec_begin(fp2);
                        for (int l = 1; l <= 2; l++)
                            for (int i = 1; i <= idim1; i++)
                                for (int k = 1; k <= kdim1; k++)
                                    for (int j = 1; j <= jdim1; j++) {
                                        fread(&dum, sizeof(double), 1, fp2);
                                    }
                        read_rec_end(fp2);
                    }

                    // NOTE: for ivmx .ge. 30, currently you CANNOT restart to/from a
                    // different ivmx; if this is ever allowed, then the following may
                    // read/write incorrectly!!!
                    if (ivmx >= 30) {
                        // read(2) nummem_read
                        {
                            read_rec_begin(fp2);
                            fread(&nummem_read, sizeof(int32_t), 1, fp2);
                            read_rec_end(fp2);
                        }
                        if (nummem_read != nummem) {
                            fprintf(fp11, " ERROR: nummem_read .ne. nummem in rrestg\n");
                            fprintf(fp11, "   cannot change ivisc with time-accurate\n");
                            fprintf(fp11, "   run for ivisc.ge.30\n");
                            int ierrflg = -1;
                            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                        }
                        // read(2) ((((tursav2(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=3,nummem_read)
                        read_rec_begin(fp2);
                        for (int l = 3; l <= nummem_read; l++)
                            for (int i = 1; i <= idim1; i++)
                                for (int k = 1; k <= kdim1; k++)
                                    for (int j = 1; j <= jdim1; j++) {
                                        double v;
                                        fread(&v, sizeof(double), 1, fp2);
                                        tursav2(j, k, i, l) = v;
                                    }
                        read_rec_end(fp2);
                    }
                } // end if (itime2read != 0)
            } else {
                // icgns == 1: empty branch
            }
        } // end if (myid == myhost)

        if (dtold != dt) {
            FILE* fp11 = fortran_get_unit(11);
            fprintf(fp11, " time step has changed: old=%12.6f, new=%12.6f\n",
                    (double)dtold, (double)dt);
            fprintf(fp11, " ...1st step will be forced to be 1st order"
                          " in time (to maintain order property)\n");
        }

        // fill in edge values of qc0 array for safety before passing
        // the data to the appropriate node
        for (int l = 1; l <= 5; l++) {
            for (int i = 1; i <= idim1; i++) {
                int k = kdim;
                for (int j = 1; j <= jdim1; j++) {
                    qc0(j, k, i, l) = (double)qiv[l-1];
                }
                int j = jdim;
                for (int k2 = 1; k2 <= kdim; k2++) {
                    qc0(j, k2, i, l) = (double)qiv[l-1];
                }
            }
        }

        // fill in edge values of tursav2 array for safety before passing
        // the data to the appropriate node
        if (ivmx >= 4) {
            for (int i = 1; i <= idim; i++) {
                int k = kdim;
                for (int j = 1; j <= jdim1; j++) {
                    for (int l = 1; l <= nummem; l++) {
                        tursav2(j, k, i, l)         = (double)tur10[l-1];
                        tursav2(j, k, i, nummem + l) = 0.0;
                    }
                }
                int j = jdim;
                for (int k2 = 1; k2 <= kdim; k2++) {
                    for (int l = 1; l <= nummem; l++) {
                        tursav2(j, k2, i, l)         = (double)tur10[l-1];
                        tursav2(j, k2, i, nummem + l) = 0.0;
                    }
                }
            }
        }

        return;

    } else if (nflagg == 1) {

        // *******************************
        //        read dynamic mesh data
        // *******************************

        if (myid == myhost) {
            if (icgns != 1) {
                // read(2) iuns
                FILE* fp2 = fortran_get_unit(2);
                read_rec_begin(fp2);
                fread(&iuns, sizeof(int32_t), 1, fp2);
                read_rec_end(fp2);
            } else {
                // empty branch
            }
        }

        if (iuns != 0) {

            if (myid == myhost) {
                if (icgns != 1) {
                    FILE* fp2  = fortran_get_unit(2);
                    FILE* fp11 = fortran_get_unit(11);

                    // write(11,105) nbl
                    // 105 format(31h reading dynamic mesh data from,23h restart file for block,i4)
                    fprintf(fp11, " reading dynamic mesh data from restart file for block%4d\n", nbl);

                    // read(2) jdum,kdum,idum
                    read_rec_begin(fp2);
                    fread(&jdum, sizeof(int32_t), 1, fp2);
                    fread(&kdum, sizeof(int32_t), 1, fp2);
                    fread(&idum, sizeof(int32_t), 1, fp2);
                    read_rec_end(fp2);

                    if (jdum != jdim && kdum != kdim && idum != idim) {
                        // write(11,100) nbl
                        // 100 format(41h error in dynamic mesh section of restart,15h file for block,i3)
                        fprintf(fp11, " error in dynamic mesh section of restart file for block%3d\n", nbl);
                        // write(11,101) jdum,kdum,idum,jdim,kdim,idim
                        // 101 format(29h restart file has j x k x i :,3i5,/,29h should be        j x k x i :,3i5)
                        fprintf(fp11, " restart file has j x k x i :%5d%5d%5d\n", jdum, kdum, idum);
                        fprintf(fp11, " should be        j x k x i :%5d%5d%5d\n", jdim, kdim, idim);
                        int ierrflg = -1;
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                } // end if (icgns != 1)
            } // end if (myid == myhost)

            if (myid == myhost) {
                if (icgns != 1) {
                    FILE* fp2  = fortran_get_unit(2);
                    FILE* fp11 = fortran_get_unit(11);

                    // read(2) itrans1,rfreqt1,...,dt1
                    read_rec_begin(fp2);
                    fread(&itrans1,    sizeof(int32_t), 1, fp2);
                    fread(&rfreqt1,    sizeof(double),  1, fp2);
                    fread(&xorig1,     sizeof(double),  1, fp2);
                    fread(&yorig1,     sizeof(double),  1, fp2);
                    fread(&zorig1,     sizeof(double),  1, fp2);
                    fread(&xorig01,    sizeof(double),  1, fp2);
                    fread(&yorig01,    sizeof(double),  1, fp2);
                    fread(&zorig01,    sizeof(double),  1, fp2);
                    fread(&utrans1,    sizeof(double),  1, fp2);
                    fread(&vtrans1,    sizeof(double),  1, fp2);
                    fread(&wtrans1,    sizeof(double),  1, fp2);
                    fread(&dxmx1,      sizeof(double),  1, fp2);
                    fread(&dymx1,      sizeof(double),  1, fp2);
                    fread(&dzmx1,      sizeof(double),  1, fp2);
                    fread(&itransmc1,  sizeof(int32_t), 1, fp2);
                    fread(&rfreqtmc1,  sizeof(double),  1, fp2);
                    fread(&xorigmc1,   sizeof(double),  1, fp2);
                    fread(&yorigmc1,   sizeof(double),  1, fp2);
                    fread(&zorigmc1,   sizeof(double),  1, fp2);
                    fread(&xorig0mc1,  sizeof(double),  1, fp2);
                    fread(&yorig0mc1,  sizeof(double),  1, fp2);
                    fread(&zorig0mc1,  sizeof(double),  1, fp2);
                    fread(&utransmc1,  sizeof(double),  1, fp2);
                    fread(&vtransmc1,  sizeof(double),  1, fp2);
                    fread(&wtransmc1,  sizeof(double),  1, fp2);
                    fread(&xmc1,       sizeof(double),  1, fp2);
                    fread(&ymc1,       sizeof(double),  1, fp2);
                    fread(&zmc1,       sizeof(double),  1, fp2);
                    fread(&dxmxmc1,    sizeof(double),  1, fp2);
                    fread(&dymxmc1,    sizeof(double),  1, fp2);
                    fread(&dzmxmc1,    sizeof(double),  1, fp2);
                    fread(&irotat1,    sizeof(int32_t), 1, fp2);
                    fread(&rfreqr1,    sizeof(double),  1, fp2);
                    fread(&thetax1,    sizeof(double),  1, fp2);
                    fread(&thetay1,    sizeof(double),  1, fp2);
                    fread(&thetaz1,    sizeof(double),  1, fp2);
                    fread(&omegax1,    sizeof(double),  1, fp2);
                    fread(&omegay1,    sizeof(double),  1, fp2);
                    fread(&omegaz1,    sizeof(double),  1, fp2);
                    fread(&dthxmx1,    sizeof(double),  1, fp2);
                    fread(&dthymx1,    sizeof(double),  1, fp2);
                    fread(&dthzmx1,    sizeof(double),  1, fp2);
                    fread(&irotatmc1,  sizeof(int32_t), 1, fp2);
                    fread(&rfreqrmc1,  sizeof(double),  1, fp2);
                    fread(&thetaxmc1,  sizeof(double),  1, fp2);
                    fread(&thetaymc1,  sizeof(double),  1, fp2);
                    fread(&thetazmc1,  sizeof(double),  1, fp2);
                    fread(&omegaxmc1,  sizeof(double),  1, fp2);
                    fread(&omegaymc1,  sizeof(double),  1, fp2);
                    fread(&omegazmc1,  sizeof(double),  1, fp2);
                    fread(&dthxmxmc1,  sizeof(double),  1, fp2);
                    fread(&dthymxmc1,  sizeof(double),  1, fp2);
                    fread(&dthzmxmc1,  sizeof(double),  1, fp2);
                    fread(&time21,     sizeof(double),  1, fp2);
                    fread(&time2mc1,   sizeof(double),  1, fp2);
                    fread(&dt1,        sizeof(double),  1, fp2);
                    read_rec_end(fp2);

                } else {
                    // empty branch (icgns == 1)
                }

                // fix-up for change in moving grid parameters
                if (irotat(nbl) == 1) {
                    if (std::abs((float)omegax1) > 0.0f &&
                        (float)omegax(nbl) != (float)omegax1) {
                        if ((float)omegay1 == 0.0f && (float)omegaz1 == 0.0f) {
                            time21 = time21 - dt1;
                            time211 = omegax1 / omegax(nbl) * time21;
                            FILE* fp11 = fortran_get_unit(11);
                            // write(11,201) nbl,igrid
                            // 201 format(/,3x,16h WARNING: block ,i3,7h (grid ,i3,1h),21h rotational speed was,22h changed upon restart:)
                            fprintf(fp11, "\n   WARNING: block %3d (grid %3d) rotational speed was changed upon restart:\n", nbl, igrid);
                            // write(11,202) real(omegax(nbl)),real(omegax1)
                            // 202 format(3x,14h new omegax = ,f10.4,3x,14h old omegax = ,f10.4)
                            fprintf(fp11, "   new omegax = %10.4f   old omegax = %10.4f\n", (float)omegax(nbl), (float)omegax1);
                            // write(11,203) real(time211),real(time21)
                            // 203 format(3x,31h to compensate, adjusting block,21h time counter, time2:,/,3x,14h new time2  = ,f10.4,3x,14h old time2  = ,f10.4,/)
                            fprintf(fp11, "   to compensate, adjusting block time counter, time2:\n   new time2  = %10.4f   old time2  = %10.4f\n\n", (float)time211, (float)time21);
                            time21 = time211 + dt1;
                        } else {
                            FILE* fp11 = fortran_get_unit(11);
                            // write(11,204)
                            // 204 format(/,3x,26h stopping...rotation speed,34h change in more than one direction,/,3x,29h not allowed between restarts)
                            fprintf(fp11, "\n   stopping...rotation speed change in more than one direction\n   not allowed between restarts\n");
                            int ierrflg = -1;
                            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                        }
                    }
                    if (std::abs((float)omegay1) > 0.0f &&
                        (float)omegay(nbl) != (float)omegay1) {
                        if ((float)omegax1 == 0.0f && (float)omegaz1 == 0.0f) {
                            time21 = time21 - dt1;
                            time211 = omegay1 / omegay(nbl) * time21;
                            FILE* fp11 = fortran_get_unit(11);
                            fprintf(fp11, "\n   WARNING: block %3d (grid %3d) rotational speed was changed upon restart:\n", nbl, igrid);
                            // write(11,205) real(omegay(nbl)),real(omegay1)
                            // 205 format(3x,14h new omegay = ,f10.4,3x,14h old omegay = ,f10.4)
                            fprintf(fp11, "   new omegay = %10.4f   old omegay = %10.4f\n", (float)omegay(nbl), (float)omegay1);
                            fprintf(fp11, "   to compensate, adjusting block time counter, time2:\n   new time2  = %10.4f   old time2  = %10.4f\n\n", (float)time211, (float)time21);
                            time21 = time211 + dt1;
                        } else {
                            FILE* fp11 = fortran_get_unit(11);
                            fprintf(fp11, "\n   stopping...rotation speed change in more than one direction\n   not allowed between restarts\n");
                            int ierrflg = -1;
                            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                        }
                    }
                    if (std::abs((float)omegaz1) > 0.0f &&
                        (float)omegaz(nbl) != (float)omegaz1) {
                        if ((float)omegax1 == 0.0f && (float)omegay1 == 0.0f) {
                            time21 = time21 - dt1;
                            time211 = omegaz1 / omegaz(nbl) * time21;
                            FILE* fp11 = fortran_get_unit(11);
                            fprintf(fp11, "\n   WARNING: block %3d (grid %3d) rotational speed was changed upon restart:\n", nbl, igrid);
                            // write(11,206) real(omegaz(nbl)),real(omegaz1)
                            // 206 format(3x,14h new omegaz = ,f10.4,3x,14h old omegaz = ,f10.4)
                            fprintf(fp11, "   new omegaz = %10.4f   old omegaz = %10.4f\n", (float)omegaz(nbl), (float)omegaz1);
                            fprintf(fp11, "   to compensate, adjusting block time counter, time2:\n   new time2  = %10.4f   old time2  = %10.4f\n\n", (float)time211, (float)time21);
                            time21 = time211 + dt1;
                        } else {
                            FILE* fp11 = fortran_get_unit(11);
                            fprintf(fp11, "\n   stopping...rotation speed change in more than one direction\n   not allowed between restarts\n");
                            int ierrflg = -1;
                            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                        }
                    }
                } // end if (irotat(nbl) == 1)


                if (irotatmc == 1) {
                    if (std::abs((float)omegaxmc1) > 0.0f &&
                        (float)omegaxmc != (float)omegaxmc1) {
                        if ((float)omegaymc1 == 0.0f && (float)omegazmc1 == 0.0f) {
                            time2mc1 = time2mc1 - dt1;
                            time2mc11 = omegaxmc1 / omegaxmc * time2mc1;
                            FILE* fp11 = fortran_get_unit(11);
                            // write(11,301)
                            // 301 format(/,3x,23h WARNING: moment center,34h rotational speed was changed upon,8h restart:)
                            fprintf(fp11, "\n   WARNING: moment center rotational speed was changed upon restart:\n");
                            // write(11,302) real(omegaxmc),real(omegaxmc1)
                            // 302 format(3x,16h new omegaxmc = ,f10.4,3x,16h old omegaxmc = ,f10.4)
                            fprintf(fp11, "   new omegaxmc = %10.4f   old omegaxmc = %10.4f\n", (float)omegaxmc, (float)omegaxmc1);
                            // write(11,303) real(time2mc11),real(time2mc1)
                            // 303 format(3x,32h to compensate, adjusting moment,30h center time counter, time2mc:,/,3x,16h new time2mc  = ,f10.4,3x,16h old time2mc  = ,f10.4,/)
                            fprintf(fp11, "   to compensate, adjusting moment center time counter, time2mc:\n   new time2mc  = %10.4f   old time2mc  = %10.4f\n\n", (float)time2mc11, (float)time2mc1);
                            time2mc1 = time2mc11 + dt;
                        } else {
                            FILE* fp11 = fortran_get_unit(11);
                            // write(11,304)
                            // 304 format(/,3x,25h stopping...moment center,25h rotation speed change in,28h more than one direction not,/,3x,25h allowed between restarts)
                            fprintf(fp11, "\n   stopping...moment center rotation speed change in more than one direction not\n   allowed between restarts\n");
                            int ierrflg = -1;
                            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                        }
                    }
                    if (std::abs((float)omegaymc1) > 0.0f &&
                        (float)omegaymc != (float)omegaymc1) {
                        if ((float)omegaxmc1 == 0.0f && (float)omegazmc1 == 0.0f) {
                            time2mc1 = time2mc1 - dt1;
                            time2mc11 = omegaymc1 / omegaymc * time2mc1;
                            FILE* fp11 = fortran_get_unit(11);
                            fprintf(fp11, "\n   WARNING: moment center rotational speed was changed upon restart:\n");
                            // write(11,305) real(omegaymc),real(omegaymc1)
                            // 305 format(3x,16h new omegaymc = ,f10.4,3x,16h old omegaymc = ,f10.4)
                            fprintf(fp11, "   new omegaymc = %10.4f   old omegaymc = %10.4f\n", (float)omegaymc, (float)omegaymc1);
                            fprintf(fp11, "   to compensate, adjusting moment center time counter, time2mc:\n   new time2mc  = %10.4f   old time2mc  = %10.4f\n\n", (float)time2mc11, (float)time2mc1);
                            time2mc1 = time2mc11 + dt1;
                        } else {
                            FILE* fp11 = fortran_get_unit(11);
                            fprintf(fp11, "\n   stopping...moment center rotation speed change in more than one direction not\n   allowed between restarts\n");
                            int ierrflg = -1;
                            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                        }
                    }
                    if (std::abs((float)omegazmc1) > 0.0f &&
                        (float)omegazmc != (float)omegazmc1) {
                        if ((float)omegaxmc1 == 0.0f && (float)omegaymc1 == 0.0f) {
                            time2mc1 = time2mc1 - dt1;
                            time2mc11 = omegazmc1 / omegazmc * time2mc1;
                            FILE* fp11 = fortran_get_unit(11);
                            fprintf(fp11, "\n   WARNING: moment center rotational speed was changed upon restart:\n");
                            // write(11,306) real(omegazmc),real(omegazmc1)
                            // 306 format(3x,16h new omegazmc = ,f10.4,3x,16h old omegazmc = ,f10.4)
                            fprintf(fp11, "   new omegazmc = %10.4f   old omegazmc = %10.4f\n", (float)omegazmc, (float)omegazmc1);
                            fprintf(fp11, "   to compensate, adjusting moment center time counter, time2mc:\n   new time2mc  = %10.4f   old time2mc  = %10.4f\n\n", (float)time2mc11, (float)time2mc1);
                            time2mc1 = time2mc11 + dt1;
                        } else {
                            FILE* fp11 = fortran_get_unit(11);
                            fprintf(fp11, "\n   stopping...moment center rotation speed change in more than one direction not\n   allowed between restarts\n");
                            int ierrflg = -1;
                            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                        }
                    }
                } // end if (irotatmc == 1)

                // set old position/time values to current values
                time2(nbl)   = time21;
                xorig(nbl)   = xorig1;
                yorig(nbl)   = yorig1;
                zorig(nbl)   = zorig1;
                thetax(nbl)  = thetax1;
                thetay(nbl)  = thetay1;
                thetaz(nbl)  = thetaz1;
                time2mc      = (float)time2mc1;
                xorigmc      = (float)xorigmc1;
                yorigmc      = (float)yorigmc1;
                zorigmc      = (float)zorigmc1;
                xmc          = (float)xmc1;
                ymc          = (float)ymc1;
                zmc          = (float)zmc1;
                thetaxmc     = (float)thetaxmc1;
                thetaymc     = (float)thetaymc1;
                thetazmc     = (float)thetazmc1;

            } // end if (myid == myhost) [second block]


            if (idefrm(nbl) > 0) {
                if (myid == myhost) {
                    if (icgns != 1) {
                        FILE* fp2 = fortran_get_unit(2);

                        // read(2) idefrm(nbl),nsegdfrm(nbl)
                        {
                            read_rec_begin(fp2);
                            int32_t idefrm_val, nsegdfrm_val;
                            fread(&idefrm_val,   sizeof(int32_t), 1, fp2);
                            fread(&nsegdfrm_val, sizeof(int32_t), 1, fp2);
                            read_rec_end(fp2);
                            idefrm(nbl)   = idefrm_val;
                            nsegdfrm(nbl) = nsegdfrm_val;
                        }

                        for (int is = 1; is <= nsegdfrm(nbl); is++) {
                            // read(2) idfrmsdum,utrnsae(nbl,is),...,kcsf(nbl,is)
                            read_rec_begin(fp2);
                            {
                                double utrnsae_v, vtrnsae_v, wtrnsae_v;
                                double omgxae_v, omgyae_v, omgzae_v;
                                double xorgae_v, yorgae_v, zorgae_v;
                                double thtxae_v, thtyae_v, thtzae_v;
                                double rfrqtae_v, rfrqrae_v;
                                int32_t icsi_v, icsf_v, jcsi_v, jcsf_v, kcsi_v, kcsf_v;

                                fread(&idfrmsdum,  sizeof(int32_t), 1, fp2);
                                fread(&utrnsae_v,  sizeof(double),  1, fp2);
                                fread(&vtrnsae_v,  sizeof(double),  1, fp2);
                                fread(&wtrnsae_v,  sizeof(double),  1, fp2);
                                fread(&omgxae_v,   sizeof(double),  1, fp2);
                                fread(&omgyae_v,   sizeof(double),  1, fp2);
                                fread(&omgzae_v,   sizeof(double),  1, fp2);
                                fread(&xorgae_v,   sizeof(double),  1, fp2);
                                fread(&yorgae_v,   sizeof(double),  1, fp2);
                                fread(&zorgae_v,   sizeof(double),  1, fp2);
                                fread(&thtxae_v,   sizeof(double),  1, fp2);
                                fread(&thtyae_v,   sizeof(double),  1, fp2);
                                fread(&thtzae_v,   sizeof(double),  1, fp2);
                                fread(&rfrqtae_v,  sizeof(double),  1, fp2);
                                fread(&rfrqrae_v,  sizeof(double),  1, fp2);
                                fread(&icsi_v,     sizeof(int32_t), 1, fp2);
                                fread(&icsf_v,     sizeof(int32_t), 1, fp2);
                                fread(&jcsi_v,     sizeof(int32_t), 1, fp2);
                                fread(&jcsf_v,     sizeof(int32_t), 1, fp2);
                                fread(&kcsi_v,     sizeof(int32_t), 1, fp2);
                                fread(&kcsf_v,     sizeof(int32_t), 1, fp2);

                                utrnsae(nbl, is) = utrnsae_v;
                                vtrnsae(nbl, is) = vtrnsae_v;
                                wtrnsae(nbl, is) = wtrnsae_v;
                                omgxae(nbl, is)  = omgxae_v;
                                omgyae(nbl, is)  = omgyae_v;
                                omgzae(nbl, is)  = omgzae_v;
                                xorgae(nbl, is)  = xorgae_v;
                                yorgae(nbl, is)  = yorgae_v;
                                zorgae(nbl, is)  = zorgae_v;
                                thtxae(nbl, is)  = thtxae_v;
                                thtyae(nbl, is)  = thtyae_v;
                                thtzae(nbl, is)  = thtzae_v;
                                rfrqtae(nbl, is) = rfrqtae_v;
                                rfrqrae(nbl, is) = rfrqrae_v;
                                icsi(nbl, is)    = icsi_v;
                                icsf(nbl, is)    = icsf_v;
                                jcsi(nbl, is)    = jcsi_v;
                                jcsf(nbl, is)    = jcsf_v;
                                kcsi(nbl, is)    = kcsi_v;
                                kcsf(nbl, is)    = kcsf_v;
                            }
                            read_rec_end(fp2);
                        } // end do is=1,nsegdfrm(nbl)

                        if (std::abs(ita) > 1) {
                            // read(2) (((xnm2(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
                            //         (((ynm2(j,k,i),...),znm2(j,k,i)...)
                            read_rec_begin(fp2);
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v; fread(&v, sizeof(double), 1, fp2);
                                        xnm2(j, k, i) = v;
                                    }
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v; fread(&v, sizeof(double), 1, fp2);
                                        ynm2(j, k, i) = v;
                                    }
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v; fread(&v, sizeof(double), 1, fp2);
                                        znm2(j, k, i) = v;
                                    }
                            read_rec_end(fp2);
                        } else {
                            // read(2) (((xdum,j=1,jdim),k=1,kdim),i=1,idim),
                            //         (((ydum,...),zdum,...)
                            read_rec_begin(fp2);
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++)
                                        fread(&xdum, sizeof(double), 1, fp2);
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++)
                                        fread(&ydum, sizeof(double), 1, fp2);
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++)
                                        fread(&zdum, sizeof(double), 1, fp2);
                            read_rec_end(fp2);
                        }

                        // read(2) (((x(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
                        //         (((y(j,k,i),...),z(j,k,i)...)
                        read_rec_begin(fp2);
                        for (int i = 1; i <= idim; i++)
                            for (int k = 1; k <= kdim; k++)
                                for (int j = 1; j <= jdim; j++) {
                                    double v; fread(&v, sizeof(double), 1, fp2);
                                    x(j, k, i) = v;
                                }
                        for (int i = 1; i <= idim; i++)
                            for (int k = 1; k <= kdim; k++)
                                for (int j = 1; j <= jdim; j++) {
                                    double v; fread(&v, sizeof(double), 1, fp2);
                                    y(j, k, i) = v;
                                }
                        for (int i = 1; i <= idim; i++)
                            for (int k = 1; k <= kdim; k++)
                                for (int j = 1; j <= jdim; j++) {
                                    double v; fread(&v, sizeof(double), 1, fp2);
                                    z(j, k, i) = v;
                                }
                        read_rec_end(fp2);

                        if (ndgrd != 0 && iunst > 1) {
                            // read(98,*)(((x(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
                            //           (((y(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
                            //           (((z(j,k,i),i=1,idim),j=1,jdim),k=1,kdim)
                            // List-directed read: i varies fastest, then j, then k
                            int ntot = jdim * kdim * idim;
                            FortranArray1D<double> tmp_x(ntot), tmp_y(ntot), tmp_z(ntot);
                            fortran_read_list_array(98, tmp_x, 1, ntot);
                            fortran_read_list_array(98, tmp_y, 1, ntot);
                            fortran_read_list_array(98, tmp_z, 1, ntot);
                            int idx = 1;
                            for (int k = 1; k <= kdim; k++)
                                for (int j = 1; j <= jdim; j++)
                                    for (int i = 1; i <= idim; i++) {
                                        x(j, k, i) = tmp_x(idx);
                                        idx++;
                                    }
                            idx = 1;
                            for (int k = 1; k <= kdim; k++)
                                for (int j = 1; j <= jdim; j++)
                                    for (int i = 1; i <= idim; i++) {
                                        y(j, k, i) = tmp_y(idx);
                                        idx++;
                                    }
                            idx = 1;
                            for (int k = 1; k <= kdim; k++)
                                for (int j = 1; j <= jdim; j++)
                                    for (int i = 1; i <= idim; i++) {
                                        z(j, k, i) = tmp_z(idx);
                                        idx++;
                                    }
                        }

                    } else {
                        // empty branch (icgns == 1)
                    }
                } // end if (myid == myhost)
            } // end if (idefrm(nbl) > 0)

        } // end if (iuns != 0)

    } // end else if (nflagg == 1)

    return;
} // end rrestg


} // namespace rrestg_ns
