// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "wrestg.h"
#include "wrest.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cmath>
#include <cstdint>

namespace wrestg_ns {

// Helper: write Fortran unformatted record begin/end markers
static void write_rec_begin(FILE* fp, int32_t nbytes) {
    fwrite(&nbytes, sizeof(int32_t), 1, fp);
}
static void write_rec_end(FILE* fp, int32_t nbytes) {
    fwrite(&nbytes, sizeof(int32_t), 1, fp);
}

// wrestg_ns::wrest is a forwarding wrapper to wrest_ns::wrest
void wrest(int& nbl, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
           FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
           int& ncycmax, FortranArray1DRef<double> rms,
           FortranArray1DRef<double> clw, FortranArray1DRef<double> cdw,
           FortranArray1DRef<double> cdpw, FortranArray1DRef<double> cdvw,
           FortranArray1DRef<double> cxw, FortranArray1DRef<double> cyw,
           FortranArray1DRef<double> czw, FortranArray1DRef<double> cmxw,
           FortranArray1DRef<double> cmyw, FortranArray1DRef<double> cmzw,
           int& n_clcd, FortranArray3DRef<double> clcd,
           int& nblocks_clcd, FortranArray2DRef<int> blocks_clcd,
           FortranArray1DRef<double> fmdotw, FortranArray1DRef<double> cftmomw,
           FortranArray1DRef<double> cftpw, FortranArray1DRef<double> cftvw,
           FortranArray1DRef<double> cfttotw, FortranArray2DRef<double> rmstr,
           FortranArray2DRef<int> nneg, int& iskip,
           FortranArray3DRef<double> vist3d, FortranArray4DRef<double> tursav,
           FortranArray3DRef<double> smin, FortranArray3DRef<double> xjb,
           FortranArray3DRef<double> xkb, FortranArray3DRef<double> blnum,
           FortranArray3DRef<double> cmuv, FortranArray1DRef<double> thetay,
           int& maxbl, int& myid, int& myhost, int& mycomm,
           FortranArray1DRef<int> mblk2nd, int& igrid,
           FortranArray3DRef<double> wk, int& idima, int& jdima, int& kdima,
           FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
           FortranArray3DRef<double> bci, FortranArray4DRef<double> vj0,
           FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0,
           FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0,
           FortranArray4DRef<double> ti0, FortranArray3DRef<double> blank,
           FortranArray1DRef<int> iwk, int& iwork, int& iover,
           FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
           int& nbuf, int& ibufdim, FortranArray4DRef<double> qavg,
           FortranArray4DRef<double> q2avg, FortranArray3DRef<double> x,
           FortranArray3DRef<double> y, FortranArray3DRef<double> z,
           int& nummem)
{
    wrest_ns::wrest(nbl, jdim, kdim, idim, q, qj0, qk0, qi0, ncycmax,
                    rms, clw, cdw, cdpw, cdvw, cxw, cyw, czw, cmxw, cmyw, cmzw,
                    n_clcd, clcd, nblocks_clcd, blocks_clcd,
                    fmdotw, cftmomw, cftpw, cftvw, cfttotw,
                    rmstr, nneg, iskip, vist3d, tursav, smin, xjb, xkb,
                    blnum, cmuv, thetay, maxbl, myid, myhost, mycomm,
                    mblk2nd, igrid, wk, idima, jdima, kdima,
                    bcj, bck, bci, vj0, vk0, vi0, tj0, tk0, ti0,
                    blank, iwk, iwork, iover, nou, bou, nbuf, ibufdim,
                    qavg, q2avg, x, y, z, nummem);
}

void wrestg(int& nbl, int& jdim, int& kdim, int& idim,
            FortranArray3DRef<double> x, FortranArray3DRef<double> y,
            FortranArray3DRef<double> z, FortranArray3DRef<double> xnm2,
            FortranArray3DRef<double> ynm2, FortranArray3DRef<double> znm2,
            FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk,
            FortranArray4DRef<double> delti, FortranArray4DRef<double> qc0,
            int& nflagg, int& iuns,
            FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans,
            FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax,
            FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz,
            FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig,
            FortranArray1DRef<double> zorig, FortranArray1DRef<double> dxmx,
            FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx,
            FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx,
            FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> thetax,
            FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz,
            FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr,
            FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0,
            FortranArray1DRef<double> zorig0, FortranArray1DRef<double> time2,
            FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl,
            FortranArray1DRef<double> thetazl, FortranArray1DRef<int> itrans,
            FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm,
            FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae,
            FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae,
            FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae,
            FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae,
            FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae,
            FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae,
            FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae,
            FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
            FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
            FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
            FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass,
            FortranArray2DRef<double> damp, FortranArray2DRef<double> x0,
            FortranArray2DRef<double> gf0, int& nmds, int& maxaes,
            FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb,
            int& myhost, int& myid, int& mycomm,
            FortranArray1DRef<int> mblk2nd, int& maxbl,
            FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg,
            FortranArray2DRef<int> iaesurf, int& maxsegdg,
            FortranArray1DRef<double> wk, int& nwork, int& idima, int& jdima,
            int& kdima, int& igrid, FortranArray4DRef<double> tursav2,
            int& nummem)
{
    // COMMON block aliases
    int32_t& ivmx      = cmn_maxiv.ivmx;
    float&   dt        = cmn_info.dt;
    int32_t& ita       = cmn_unst.ita;
    int32_t& icgns     = cmn_cgns.icgns;
    int32_t& itransmc  = cmn_motionmc.itransmc;
    int32_t& irotatmc  = cmn_motionmc.irotatmc;
    float&   rfreqtmc  = cmn_motionmc.rfreqtmc;
    float&   rfreqrmc  = cmn_motionmc.rfreqrmc;
    float&   xorigmc   = cmn_motionmc.xorigmc;
    float&   yorigmc   = cmn_motionmc.yorigmc;
    float&   zorigmc   = cmn_motionmc.zorigmc;
    float&   xorig0mc  = cmn_motionmc.xorig0mc;
    float&   yorig0mc  = cmn_motionmc.yorig0mc;
    float&   zorig0mc  = cmn_motionmc.zorig0mc;
    float&   utransmc  = cmn_motionmc.utransmc;
    float&   vtransmc  = cmn_motionmc.vtransmc;
    float&   wtransmc  = cmn_motionmc.wtransmc;
    float&   thetaxmc  = cmn_motionmc.thetaxmc;
    float&   thetaymc  = cmn_motionmc.thetaymc;
    float&   thetazmc  = cmn_motionmc.thetazmc;
    float&   omegaxmc  = cmn_motionmc.omegaxmc;
    float&   omegaymc  = cmn_motionmc.omegaymc;
    float&   omegazmc  = cmn_motionmc.omegazmc;
    float&   dxmxmc    = cmn_motionmc.dxmxmc;
    float&   dymxmc    = cmn_motionmc.dymxmc;
    float&   dzmxmc    = cmn_motionmc.dzmxmc;
    float&   dthxmxmc  = cmn_motionmc.dthxmxmc;
    float&   dthymxmc  = cmn_motionmc.dthymxmc;
    float&   dthzmxmc  = cmn_motionmc.dthzmxmc;
    float&   time2mc   = cmn_motionmc.time2mc;
    float&   xmc       = cmn_fsum.xmc;
    float&   ymc       = cmn_fsum.ymc;
    float&   zmc       = cmn_fsum.zmc;
    int32_t& ndwrt     = cmn_deformz.ndwrt;

    // Local variables
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    FILE* fp2  = fortran_get_unit(2);
    FILE* fp11 = fortran_get_unit(11);
    FILE* fp98 = fortran_get_unit(98);

    if (nflagg == 0) {
        //
        // ************************************
        //        write out 2nd order time data
        // ************************************
        //
        if (myid == myhost) {
            if (icgns != 1) {
                // write(11,101) nbl
                // 101 format(45h writing 2nd order time data to restart file,,
                //            6h block,i4)
                fprintf(fp11, "writing 2nd order time data to restart file, block%4d\n", nbl);

                // write(2) jdim,kdim,idim
                {
                    int32_t rec = 3 * (int32_t)sizeof(int32_t);
                    write_rec_begin(fp2, rec);
                    int32_t v;
                    v = (int32_t)jdim; fwrite(&v, sizeof(int32_t), 1, fp2);
                    v = (int32_t)kdim; fwrite(&v, sizeof(int32_t), 1, fp2);
                    v = (int32_t)idim; fwrite(&v, sizeof(int32_t), 1, fp2);
                    write_rec_end(fp2, rec);
                }

                // write(2) ((((qc0(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,5)
                {
                    int32_t rec = (int32_t)(5 * jdim1 * kdim1 * idim1) * (int32_t)sizeof(double);
                    write_rec_begin(fp2, rec);
                    for (int l = 1; l <= 5; l++)
                        for (int i = 1; i <= idim1; i++)
                            for (int k = 1; k <= kdim1; k++)
                                for (int j = 1; j <= jdim1; j++) {
                                    double v = qc0(j, k, i, l);
                                    fwrite(&v, sizeof(double), 1, fp2);
                                }
                    write_rec_end(fp2, rec);
                }

                // write(2) dt
                {
                    int32_t rec = (int32_t)sizeof(float);
                    write_rec_begin(fp2, rec);
                    fwrite(&dt, sizeof(float), 1, fp2);
                    write_rec_end(fp2, rec);
                }

                // only need to store old values of turb quantities (j,k,i,1) and
                // (j,k,i,2), NOT the deltaQ's = (j,k,i,3) and (j,k,i,4):
                if (ivmx >= 4) {
                    // write(2) ((((tursav2(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,2)
                    int32_t rec = (int32_t)(2 * jdim1 * kdim1 * idim1) * (int32_t)sizeof(double);
                    write_rec_begin(fp2, rec);
                    for (int l = 1; l <= 2; l++)
                        for (int i = 1; i <= idim1; i++)
                            for (int k = 1; k <= kdim1; k++)
                                for (int j = 1; j <= jdim1; j++) {
                                    double v = tursav2(j, k, i, l);
                                    fwrite(&v, sizeof(double), 1, fp2);
                                }
                    write_rec_end(fp2, rec);
                } else {
                    // write(2) ((((zero,j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,2)
                    // zero=0. is a Fortran real (float)
                    float zero = 0.0f;
                    int32_t rec = (int32_t)(2 * jdim1 * kdim1 * idim1) * (int32_t)sizeof(float);
                    write_rec_begin(fp2, rec);
                    for (int l = 1; l <= 2; l++)
                        for (int i = 1; i <= idim1; i++)
                            for (int k = 1; k <= kdim1; k++)
                                for (int j = 1; j <= jdim1; j++)
                                    fwrite(&zero, sizeof(float), 1, fp2);
                    write_rec_end(fp2, rec);
                }

                // NOTE: for ivmx .ge. 30, currently you CANNOT restart to/from a
                // different ivmx; if this is ever allowed, then the following may
                // read/write incorrectly!!!
                if (ivmx >= 30) {
                    // write(2) nummem
                    {
                        int32_t rec = (int32_t)sizeof(int32_t);
                        write_rec_begin(fp2, rec);
                        int32_t v = (int32_t)nummem;
                        fwrite(&v, sizeof(int32_t), 1, fp2);
                        write_rec_end(fp2, rec);
                    }
                    // write(2) ((((tursav2(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=3,nummem)
                    {
                        int count = nummem - 2;
                        int32_t rec = (int32_t)(count * jdim1 * kdim1 * idim1) * (int32_t)sizeof(double);
                        write_rec_begin(fp2, rec);
                        for (int l = 3; l <= nummem; l++)
                            for (int i = 1; i <= idim1; i++)
                                for (int k = 1; k <= kdim1; k++)
                                    for (int j = 1; j <= jdim1; j++) {
                                        double v = tursav2(j, k, i, l);
                                        fwrite(&v, sizeof(double), 1, fp2);
                                    }
                        write_rec_end(fp2, rec);
                    }
                }
            } // icgns != 1
        } // myid == myhost

    } else if (nflagg == 1) {
        //
        // ************************************
        //        write out dynamic mesh data
        // ************************************
        //
        if (myid == myhost) {
            if (icgns != 1) {
                // write(2) iuns
                int32_t rec = (int32_t)sizeof(int32_t);
                write_rec_begin(fp2, rec);
                int32_t v = (int32_t)iuns;
                fwrite(&v, sizeof(int32_t), 1, fp2);
                write_rec_end(fp2, rec);
            }
        }

        if (iuns != 0) {

            if (myid == myhost) {
                if (icgns != 1) {
                    // write(11,102) nbl
                    // 102 format(43h writing dynamic mesh data to restart file,,
                    //            8h   block,i4)
                    fprintf(fp11, "writing dynamic mesh data to restart file,    block%4d\n", nbl);

                    // write(2) jdim,kdim,idim
                    int32_t rec = 3 * (int32_t)sizeof(int32_t);
                    write_rec_begin(fp2, rec);
                    int32_t vj = (int32_t)jdim, vk = (int32_t)kdim, vi = (int32_t)idim;
                    fwrite(&vj, sizeof(int32_t), 1, fp2);
                    fwrite(&vk, sizeof(int32_t), 1, fp2);
                    fwrite(&vi, sizeof(int32_t), 1, fp2);
                    write_rec_end(fp2, rec);
                }
            }

            if (myid == myhost) {
                if (icgns != 1) {
                    // write(2) itrans(nbl),rfreqt(nbl),...,dt
                    // Sizes: itrans[int32], rfreqt..dzmx[14 doubles], itransmc[int32],
                    //        rfreqtmc..dzmxmc[15 floats], irotat[int32],
                    //        rfreqr..dthzmx[11 doubles], irotatmc[int32],
                    //        rfreqrmc..dthzmxmc[9 floats], time2[double], time2mc[float], dt[float]
                    int32_t rec =
                        (int32_t)sizeof(int32_t)          // itrans(nbl)
                        + 14 * (int32_t)sizeof(double)    // rfreqt..dzmx
                        + (int32_t)sizeof(int32_t)        // itransmc
                        + 15 * (int32_t)sizeof(float)     // rfreqtmc..dzmxmc
                        + (int32_t)sizeof(int32_t)        // irotat(nbl)
                        + 11 * (int32_t)sizeof(double)    // rfreqr..dthzmx
                        + (int32_t)sizeof(int32_t)        // irotatmc
                        + 9 * (int32_t)sizeof(float)      // rfreqrmc..dthzmxmc
                        + (int32_t)sizeof(double)         // time2(nbl)
                        + (int32_t)sizeof(float)          // time2mc
                        + (int32_t)sizeof(float);         // dt
                    write_rec_begin(fp2, rec);
                    {
                        int32_t vi2; double vd; float vf;
                        vi2 = (int32_t)itrans(nbl);         fwrite(&vi2, sizeof(int32_t), 1, fp2);
                        vd = rfreqt(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = xorig(nbl);                    fwrite(&vd, sizeof(double), 1, fp2);
                        vd = yorig(nbl);                    fwrite(&vd, sizeof(double), 1, fp2);
                        vd = zorig(nbl);                    fwrite(&vd, sizeof(double), 1, fp2);
                        vd = xorig0(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = yorig0(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = zorig0(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = utrans(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = vtrans(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = wtrans(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = dxmx(nbl);                     fwrite(&vd, sizeof(double), 1, fp2);
                        vd = dymx(nbl);                     fwrite(&vd, sizeof(double), 1, fp2);
                        vd = dzmx(nbl);                     fwrite(&vd, sizeof(double), 1, fp2);
                        vi2 = (int32_t)itransmc;            fwrite(&vi2, sizeof(int32_t), 1, fp2);
                        vf = rfreqtmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = xorigmc;                       fwrite(&vf, sizeof(float), 1, fp2);
                        vf = yorigmc;                       fwrite(&vf, sizeof(float), 1, fp2);
                        vf = zorigmc;                       fwrite(&vf, sizeof(float), 1, fp2);
                        vf = xorig0mc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = yorig0mc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = zorig0mc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = utransmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = vtransmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = wtransmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = xmc;                           fwrite(&vf, sizeof(float), 1, fp2);
                        vf = ymc;                           fwrite(&vf, sizeof(float), 1, fp2);
                        vf = zmc;                           fwrite(&vf, sizeof(float), 1, fp2);
                        vf = dxmxmc;                        fwrite(&vf, sizeof(float), 1, fp2);
                        vf = dymxmc;                        fwrite(&vf, sizeof(float), 1, fp2);
                        vf = dzmxmc;                        fwrite(&vf, sizeof(float), 1, fp2);
                        vi2 = (int32_t)irotat(nbl);         fwrite(&vi2, sizeof(int32_t), 1, fp2);
                        vd = rfreqr(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = thetax(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = thetay(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = thetaz(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = omegax(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = omegay(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = omegaz(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = dthxmx(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = dthymx(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vd = dthzmx(nbl);                   fwrite(&vd, sizeof(double), 1, fp2);
                        vi2 = (int32_t)irotatmc;            fwrite(&vi2, sizeof(int32_t), 1, fp2);
                        vf = rfreqrmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = thetaxmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = thetaymc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = thetazmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = omegaxmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = omegaymc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = omegazmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = dthxmxmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = dthymxmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vf = dthzmxmc;                      fwrite(&vf, sizeof(float), 1, fp2);
                        vd = time2(nbl);                    fwrite(&vd, sizeof(double), 1, fp2);
                        vf = time2mc;                       fwrite(&vf, sizeof(float), 1, fp2);
                        vf = dt;                            fwrite(&vf, sizeof(float), 1, fp2);
                    }
                    write_rec_end(fp2, rec);
                } // icgns != 1
            } // myid == myhost

            if (idefrm(nbl) > 0) {
                //
                // the restart file is written out at time step n; we need
                // to reconstruct the grid at time step n, and for second
                // order time accurate solutions, we also need the grid
                // at time n-1 (first order does not require the grid at
                // n-1). However, for the first order case the grid
                // at n is written out twice, so that the restart file
                // always contains the same amount of data.
                //
                if (myid == myhost) {
                    if (icgns != 1) {
                        // write(2) idefrm(nbl),nsegdfrm(nbl)
                        {
                            int32_t rec = 2 * (int32_t)sizeof(int32_t);
                            write_rec_begin(fp2, rec);
                            int32_t v1 = (int32_t)idefrm(nbl);
                            int32_t v2 = (int32_t)nsegdfrm(nbl);
                            fwrite(&v1, sizeof(int32_t), 1, fp2);
                            fwrite(&v2, sizeof(int32_t), 1, fp2);
                            write_rec_end(fp2, rec);
                        }

                        // do is=1,nsegdfrm(nbl)
                        //    write(2) idfrmseg(nbl,is),...,kcsf(nbl,is)
                        for (int is = 1; is <= nsegdfrm(nbl); is++) {
                            // Fields: idfrmseg[int32], utrnsae..rfrqrae[15 doubles],
                            //         icsi,icsf,jcsi,jcsf,kcsi,kcsf[6 int32]
                            int32_t rec =
                                (int32_t)sizeof(int32_t)       // idfrmseg
                                + 15 * (int32_t)sizeof(double) // utrnsae..rfrqrae
                                + 6 * (int32_t)sizeof(int32_t);// icsi..kcsf
                            write_rec_begin(fp2, rec);
                            {
                                int32_t vi2; double vd;
                                vi2 = (int32_t)idfrmseg(nbl, is); fwrite(&vi2, sizeof(int32_t), 1, fp2);
                                vd = utrnsae(nbl, is);            fwrite(&vd, sizeof(double), 1, fp2);
                                vd = vtrnsae(nbl, is);            fwrite(&vd, sizeof(double), 1, fp2);
                                vd = wtrnsae(nbl, is);            fwrite(&vd, sizeof(double), 1, fp2);
                                vd = omgxae(nbl, is);             fwrite(&vd, sizeof(double), 1, fp2);
                                vd = omgyae(nbl, is);             fwrite(&vd, sizeof(double), 1, fp2);
                                vd = omgzae(nbl, is);             fwrite(&vd, sizeof(double), 1, fp2);
                                vd = xorgae(nbl, is);             fwrite(&vd, sizeof(double), 1, fp2);
                                vd = yorgae(nbl, is);             fwrite(&vd, sizeof(double), 1, fp2);
                                vd = zorgae(nbl, is);             fwrite(&vd, sizeof(double), 1, fp2);
                                vd = thtxae(nbl, is);             fwrite(&vd, sizeof(double), 1, fp2);
                                vd = thtyae(nbl, is);             fwrite(&vd, sizeof(double), 1, fp2);
                                vd = thtzae(nbl, is);             fwrite(&vd, sizeof(double), 1, fp2);
                                vd = rfrqtae(nbl, is);            fwrite(&vd, sizeof(double), 1, fp2);
                                vd = rfrqrae(nbl, is);            fwrite(&vd, sizeof(double), 1, fp2);
                                vi2 = (int32_t)icsi(nbl, is);     fwrite(&vi2, sizeof(int32_t), 1, fp2);
                                vi2 = (int32_t)icsf(nbl, is);     fwrite(&vi2, sizeof(int32_t), 1, fp2);
                                vi2 = (int32_t)jcsi(nbl, is);     fwrite(&vi2, sizeof(int32_t), 1, fp2);
                                vi2 = (int32_t)jcsf(nbl, is);     fwrite(&vi2, sizeof(int32_t), 1, fp2);
                                vi2 = (int32_t)kcsi(nbl, is);     fwrite(&vi2, sizeof(int32_t), 1, fp2);
                                vi2 = (int32_t)kcsf(nbl, is);     fwrite(&vi2, sizeof(int32_t), 1, fp2);
                            }
                            write_rec_end(fp2, rec);
                        }

                        // xnm2, etc. only needed (and stored) for 2nd order, so
                        // write out x,y,z twice in 1st order case in order to
                        // simplify restart logic
                        if (std::abs(ita) > 1) {
                            // write(2) (((xnm2(j,k,i),...),i=1,idim),
                            //           (((ynm2(j,k,i),...),i=1,idim),
                            //           (((znm2(j,k,i),...),i=1,idim)
                            int32_t rec = (int32_t)(3 * jdim * kdim * idim) * (int32_t)sizeof(double);
                            write_rec_begin(fp2, rec);
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v = xnm2(j, k, i);
                                        fwrite(&v, sizeof(double), 1, fp2);
                                    }
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v = ynm2(j, k, i);
                                        fwrite(&v, sizeof(double), 1, fp2);
                                    }
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v = znm2(j, k, i);
                                        fwrite(&v, sizeof(double), 1, fp2);
                                    }
                            write_rec_end(fp2, rec);
                        } else {
                            // write(2) (((x(j,k,i),...),i=1,idim),
                            //           (((y(j,k,i),...),i=1,idim),
                            //           (((z(j,k,i),...),i=1,idim)
                            int32_t rec = (int32_t)(3 * jdim * kdim * idim) * (int32_t)sizeof(double);
                            write_rec_begin(fp2, rec);
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v = x(j, k, i);
                                        fwrite(&v, sizeof(double), 1, fp2);
                                    }
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v = y(j, k, i);
                                        fwrite(&v, sizeof(double), 1, fp2);
                                    }
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v = z(j, k, i);
                                        fwrite(&v, sizeof(double), 1, fp2);
                                    }
                            write_rec_end(fp2, rec);
                        }

                        // write(2) (((x(j,k,i),...),i=1,idim),
                        //           (((y(j,k,i),...),i=1,idim),
                        //           (((z(j,k,i),...),i=1,idim)
                        {
                            int32_t rec = (int32_t)(3 * jdim * kdim * idim) * (int32_t)sizeof(double);
                            write_rec_begin(fp2, rec);
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v = x(j, k, i);
                                        fwrite(&v, sizeof(double), 1, fp2);
                                    }
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v = y(j, k, i);
                                        fwrite(&v, sizeof(double), 1, fp2);
                                    }
                            for (int i = 1; i <= idim; i++)
                                for (int k = 1; k <= kdim; k++)
                                    for (int j = 1; j <= jdim; j++) {
                                        double v = z(j, k, i);
                                        fwrite(&v, sizeof(double), 1, fp2);
                                    }
                            write_rec_end(fp2, rec);
                        }

                        if (ndwrt != 0) {
                            // write(98,'(5e21.12)')
                            //   (((x(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
                            //   (((y(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
                            //   (((z(j,k,i),i=1,idim),j=1,jdim),k=1,kdim)
                            // Note: loop order is k outer, j middle, i inner
                            int count = 0;
                            for (int k = 1; k <= kdim; k++)
                                for (int j = 1; j <= jdim; j++)
                                    for (int i = 1; i <= idim; i++) {
                                        fprintf(fp98, "%21.12e", x(j, k, i));
                                        count++;
                                        if (count % 5 == 0) fprintf(fp98, "\n");
                                    }
                            for (int k = 1; k <= kdim; k++)
                                for (int j = 1; j <= jdim; j++)
                                    for (int i = 1; i <= idim; i++) {
                                        fprintf(fp98, "%21.12e", y(j, k, i));
                                        count++;
                                        if (count % 5 == 0) fprintf(fp98, "\n");
                                    }
                            for (int k = 1; k <= kdim; k++)
                                for (int j = 1; j <= jdim; j++)
                                    for (int i = 1; i <= idim; i++) {
                                        fprintf(fp98, "%21.12e", z(j, k, i));
                                        count++;
                                        if (count % 5 == 0) fprintf(fp98, "\n");
                                    }
                            if (count % 5 != 0) fprintf(fp98, "\n");
                        }

                    } // icgns != 1
                } // myid == myhost
            } // idefrm(nbl) > 0
        } // iuns != 0
    } // nflagg == 1

    return;
} // end wrestg

} // namespace wrestg_ns
