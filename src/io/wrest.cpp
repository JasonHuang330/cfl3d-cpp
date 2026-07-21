// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "wrest.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <algorithm>

namespace wrest_ns {

static void write_rec_mark(FILE* fp, int32_t nbytes) {
    fwrite(&nbytes, sizeof(int32_t), 1, fp);
}

void wrest(int& nbl, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
           FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
           int& ncycmax, FortranArray1DRef<double> rms, FortranArray1DRef<double> clw,
           FortranArray1DRef<double> cdw, FortranArray1DRef<double> cdpw,
           FortranArray1DRef<double> cdvw, FortranArray1DRef<double> cxw,
           FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw,
           FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw,
           FortranArray1DRef<double> cmzw, int& n_clcd,
           FortranArray3DRef<double> clcd, int& nblocks_clcd,
           FortranArray2DRef<int> blocks_clcd, FortranArray1DRef<double> fmdotw,
           FortranArray1DRef<double> cftmomw, FortranArray1DRef<double> cftpw,
           FortranArray1DRef<double> cftvw, FortranArray1DRef<double> cfttotw,
           FortranArray2DRef<double> rmstr, FortranArray2DRef<int> nneg,
           int& iskip, FortranArray3DRef<double> vist3d,
           FortranArray4DRef<double> tursav, FortranArray3DRef<double> smin,
           FortranArray3DRef<double> xjb, FortranArray3DRef<double> xkb,
           FortranArray3DRef<double> blnum, FortranArray3DRef<double> cmuv,
           FortranArray1DRef<double> thetay, int& maxbl,
           int& myid, int& myhost, int& mycomm,
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
    // COMMON block aliases
    float& radtodeg  = cmn_conversion.radtodeg;
    float* title     = cmn_info.title;
    float& xmach     = cmn_info.xmach;
    float& alpha     = cmn_info.alpha;
    float& beta      = cmn_info.beta;
    int&   ntt       = cmn_info.ntt;
    float& time      = cmn_unst.time;
    float& reue      = cmn_reyue.reue;
    int*   ivisc     = cmn_reyue.ivisc;
    int&   iwghost   = cmn_ghost.iwghost;
    float& edvislim  = cmn_turbconv.edvislim;
    int&   iwarneddy = cmn_turbconv.iwarneddy;
    int&   icgns     = cmn_cgns.icgns;
    int&   ialph     = cmn_igrdtyp.ialph;
    float& xnumavg   = cmn_avgdata.xnumavg;
    float& xnumavg2  = cmn_avgdata.xnumavg2;
    int&   iteravg   = cmn_avgdata.iteravg;
    int&   ipertavg  = cmn_avgdata.ipertavg;
    int&   iclcd     = cmn_avgdata.iclcd;

    // Local variables
    float alphw_f = (float)((double)radtodeg * ((double)alpha + thetay(nbl)));
    int idim1 = idim - 1;
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int jki   = jdim * kdim * idim;
    int jkim  = jdim1 * kdim1 * idim1;
    int jki2  = 2 * jki;
    int jki5  = 5 * jki;
    int jkim5 = jkim * 5;
    int jk2   = jdim * kdim * 2;
    int ki2   = kdim * idim1 * 2;
    int ji2   = jdim * idim1 * 2;
    int jk20  = jdim * kdim * 20;
    int ki20  = kdim * idim1 * 20;
    int ji20  = jdim * idim1 * 20;
    int jk4   = jdim * kdim * 4;
    int ki4   = kdim * idim1 * 4;
    int ji4   = jdim * idim1 * 4;
    int jk8   = jdim * kdim * 8;
    int ki8   = kdim * idim1 * 8;
    int ji8   = jdim * idim1 * 8;

    if (ivisc[0] == 30 || ivisc[1] == 30 || ivisc[2] == 30) {
        jki2 = 3 * jki;
        jk8  = jdim * kdim * 12;
        ki8  = kdim * idim1 * 12;
        ji8  = jdim * idim1 * 12;
    }
    if (ivisc[0] == 40 || ivisc[1] == 40 || ivisc[2] == 40) {
        jki2 = 4 * jki;
        jk8  = jdim * kdim * 16;
        ki8  = kdim * idim1 * 16;
        ji8  = jdim * idim1 * 16;
    }
    if (ivisc[0] == 72 || ivisc[1] == 72 || ivisc[2] == 72) {
        jki2 = 7 * jki;
        jk8  = jdim * kdim * 28;
        ki8  = kdim * idim1 * 28;
        ji8  = jdim * idim1 * 28;
    }
    (void)jki2; (void)jki5; (void)jkim5; (void)jk2; (void)ki2; (void)ji2;
    (void)jk20; (void)ki20; (void)ji20; (void)jk4; (void)ki4; (void)ji4;
    (void)jk8; (void)ki8; (void)ji8; (void)jkim;

    // write(11,...) log messages and first write(2) header record
    if (myid == myhost) {
        FILE* f11 = fortran_get_unit(11);
        if (iskip > 0) {
            fprintf(f11, "\n writing restart file for block%6d\n", nbl);
        } else {
            fprintf(f11, " writing restart file for block%6d\n", nbl);
        }

        if (iskip > 0) {
            for (int n = 1; n <= ntt; n++) {
                for (int l = 1; l <= nummem; l++) {
                    if ((float)rmstr(n, l) <= 0.f) rmstr(n, l) = 1.;
                }
            }
        }

        if (icgns != 1) {
            int jt = jdim;
            int kt = kdim;
            int it = idim;
            // write(2) title,xmach,jt,kt,it,alphw,reue,ntt,time
            FILE* fp2 = fortran_get_unit(2);
            int32_t rec_bytes = (int32_t)(
                20 * sizeof(float) + sizeof(float) +
                3 * sizeof(int32_t) + sizeof(float) +
                sizeof(float) + sizeof(int32_t) + sizeof(float));
            write_rec_mark(fp2, rec_bytes);
            fwrite(title, sizeof(float), 20, fp2);
            fwrite(&xmach,   sizeof(float),   1, fp2);
            fwrite(&jt,      sizeof(int32_t), 1, fp2);
            fwrite(&kt,      sizeof(int32_t), 1, fp2);
            fwrite(&it,      sizeof(int32_t), 1, fp2);
            fwrite(&alphw_f, sizeof(float),   1, fp2);
            fwrite(&reue,    sizeof(float),   1, fp2);
            fwrite(&ntt,     sizeof(int32_t), 1, fp2);
            fwrite(&time,    sizeof(float),   1, fp2);
            write_rec_mark(fp2, rec_bytes);

            if (iskip > 0) {
                // write(2) (rms(n),n=1,ntt),(clw(n),n=1,ntt),...16 arrays
                FILE* fp2b = fortran_get_unit(2);
                int32_t rb = (int32_t)(16 * ntt * sizeof(double));
                write_rec_mark(fp2b, rb);
                for (int n=1;n<=ntt;n++){double v=rms(n);     fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=clw(n);     fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cdw(n);     fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cdpw(n);    fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cdvw(n);    fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cxw(n);     fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cyw(n);     fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=czw(n);     fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cmxw(n);    fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cmyw(n);    fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cmzw(n);    fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=fmdotw(n);  fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cftmomw(n); fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cftpw(n);   fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cftvw(n);   fwrite(&v,sizeof(double),1,fp2b);}
                for (int n=1;n<=ntt;n++){double v=cfttotw(n); fwrite(&v,sizeof(double),1,fp2b);}
                write_rec_mark(fp2b, rb);
            }
        } // icgns != 1
    } // myid == myhost (first block)


    // Second myid==myhost block: write q, ghost q, ivisc, rmstr/nneg
    if (myid == myhost) {
        if (icgns != 1) {
            // write(2) ((((q(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,5)
            {
                FILE* fp2 = fortran_get_unit(2);
                int32_t rb = (int32_t)(5 * jdim1 * kdim1 * idim1 * sizeof(double));
                write_rec_mark(fp2, rb);
                for (int l=1;l<=5;l++)
                    for (int i=1;i<=idim1;i++)
                        for (int k=1;k<=kdim1;k++)
                            for (int j=1;j<=jdim1;j++){
                                double v=q(j,k,i,l); fwrite(&v,sizeof(double),1,fp2);}
                write_rec_mark(fp2, rb);
            }
            if (iwghost != 0) {
                // write(2) qi0,qj0,qk0 ghost data
                FILE* fp2 = fortran_get_unit(2);
                int32_t rb = (int32_t)(
                    (4*5*jdim1*kdim1 + 4*5*kdim1*idim1 + 4*5*jdim1*idim1)*sizeof(double));
                write_rec_mark(fp2, rb);
                // qi0(j,k,l,m): j=1,jdim1; k=1,kdim1; l=1,5; m=1,4
                for (int m=1;m<=4;m++) for (int l=1;l<=5;l++)
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++){
                        double v=qi0(j,k,l,m); fwrite(&v,sizeof(double),1,fp2);}
                // qj0(k,i,l,m): k=1,kdim1; i=1,idim1; l=1,5; m=1,4
                for (int m=1;m<=4;m++) for (int l=1;l<=5;l++)
                    for (int i=1;i<=idim1;i++) for (int k=1;k<=kdim1;k++){
                        double v=qj0(k,i,l,m); fwrite(&v,sizeof(double),1,fp2);}
                // qk0(j,i,l,m): j=1,jdim1; i=1,idim1; l=1,5; m=1,4
                for (int m=1;m<=4;m++) for (int l=1;l<=5;l++)
                    for (int i=1;i<=idim1;i++) for (int j=1;j<=jdim1;j++){
                        double v=qk0(j,i,l,m); fwrite(&v,sizeof(double),1,fp2);}
                write_rec_mark(fp2, rb);
            }

            // write(2) ivisc(1),ivisc(2),ivisc(3)
            {
                FILE* fp2 = fortran_get_unit(2);
                int32_t rb = (int32_t)(3 * sizeof(int32_t));
                write_rec_mark(fp2, rb);
                fwrite(&ivisc[0], sizeof(int32_t), 1, fp2);
                fwrite(&ivisc[1], sizeof(int32_t), 1, fp2);
                fwrite(&ivisc[2], sizeof(int32_t), 1, fp2);
                write_rec_mark(fp2, rb);
            }

            if (iskip > 0) {
                // write(2) (rmstr(n,1),n=1,ntt),(rmstr(n,2),n=1,ntt),
                //          (nneg(n,1),n=1,ntt),(nneg(n,2),n=1,ntt)
                FILE* fp2 = fortran_get_unit(2);
                int32_t rb = (int32_t)(2*ntt*sizeof(double) + 2*ntt*sizeof(int32_t));
                write_rec_mark(fp2, rb);
                for (int n=1;n<=ntt;n++){double v=rmstr(n,1); fwrite(&v,sizeof(double),1,fp2);}
                for (int n=1;n<=ntt;n++){double v=rmstr(n,2); fwrite(&v,sizeof(double),1,fp2);}
                for (int n=1;n<=ntt;n++){int32_t v=nneg(n,1); fwrite(&v,sizeof(int32_t),1,fp2);}
                for (int n=1;n<=ntt;n++){int32_t v=nneg(n,2); fwrite(&v,sizeof(int32_t),1,fp2);}
                write_rec_mark(fp2, rb);

                if (ivisc[0]>=30 || ivisc[1]>=30 || ivisc[2]>=30) {
                    // write(2) nummem
                    {
                        FILE* fp2b = fortran_get_unit(2);
                        int32_t rb2 = (int32_t)sizeof(int32_t);
                        write_rec_mark(fp2b, rb2);
                        fwrite(&nummem, sizeof(int32_t), 1, fp2b);
                        write_rec_mark(fp2b, rb2);
                    }
                    for (int l=3; l<=nummem; l++) {
                        // write(2) (rmstr(n,l),n=1,ntt),(nneg(n,l),n=1,ntt)
                        FILE* fp2b = fortran_get_unit(2);
                        int32_t rb2 = (int32_t)(ntt*sizeof(double) + ntt*sizeof(int32_t));
                        write_rec_mark(fp2b, rb2);
                        for (int n=1;n<=ntt;n++){double v=rmstr(n,l); fwrite(&v,sizeof(double),1,fp2b);}
                        for (int n=1;n<=ntt;n++){int32_t v=nneg(n,l); fwrite(&v,sizeof(int32_t),1,fp2b);}
                        write_rec_mark(fp2b, rb2);
                    }
                }
            }
        } // icgns != 1
    } // myid == myhost (second block)


    // ivisc >= 2: write vist3d and ghost viscosity data
    if (ivisc[0]>=2 || ivisc[1]>=2 || ivisc[2]>=2) {
        if (myid == myhost) {
            // check eddy viscosity limit
            if (iwarneddy != -1) {
                for (int j=1;j<=jdim1;j++)
                    for (int k=1;k<=kdim1;k++)
                        for (int i=1;i<=idim1;i++)
                            if ((float)vist3d(j,k,i) == (float)edvislim)
                                iwarneddy = 1;
            }
            {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, "   writing vist3d data to restart file, block%6d\n", nbl);
            }
            if (icgns != 1) {
                // write(2) (((vist3d(j,k,i),j=1,jdim1),k=1,kdim1),i=1,idim1)
                {
                    FILE* fp2 = fortran_get_unit(2);
                    int32_t rb = (int32_t)(jdim1*kdim1*idim1*sizeof(double));
                    write_rec_mark(fp2, rb);
                    for (int i=1;i<=idim1;i++) for (int k=1;k<=kdim1;k++)
                        for (int j=1;j<=jdim1;j++){
                            double v=vist3d(j,k,i); fwrite(&v,sizeof(double),1,fp2);}
                    write_rec_mark(fp2, rb);
                }
                if (iwghost != 0) {
                    // write(2) vi0,vj0,vk0 ghost viscosity
                    FILE* fp2 = fortran_get_unit(2);
                    int32_t rb = (int32_t)(
                        (4*1*jdim*kdim + 4*1*kdim*idim1 + 4*1*jdim*idim1)*sizeof(double));
                    write_rec_mark(fp2, rb);
                    // vi0(j,k,l,m): j=1,jdim; k=1,kdim; l=1,1; m=1,4
                    for (int m=1;m<=4;m++) for (int l=1;l<=1;l++)
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++){
                            double v=vi0(j,k,l,m); fwrite(&v,sizeof(double),1,fp2);}
                    // vj0(k,i,l,m): k=1,kdim; i=1,idim1; l=1,1; m=1,4
                    for (int m=1;m<=4;m++) for (int l=1;l<=1;l++)
                        for (int i=1;i<=idim1;i++) for (int k=1;k<=kdim;k++){
                            double v=vj0(k,i,l,m); fwrite(&v,sizeof(double),1,fp2);}
                    // vk0(j,i,l,m): j=1,jdim; i=1,idim1; l=1,1; m=1,4
                    for (int m=1;m<=4;m++) for (int l=1;l<=1;l++)
                        for (int i=1;i<=idim1;i++) for (int j=1;j<=jdim;j++){
                            double v=vk0(j,i,l,m); fwrite(&v,sizeof(double),1,fp2);}
                    write_rec_mark(fp2, rb);
                }
            }
        } // myid == myhost
    } // ivisc >= 2


    // ivisc >= 4: write turbulence quantities and smin
    if (ivisc[0]>=4 || ivisc[1]>=4 || ivisc[2]>=4) {
        if (myid == myhost) {
            {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, "   writing field eqn turb quantities to restart file, block%6d\n", nbl);
            }
            if (icgns != 1) {
                // write(2) ((((tursav(j,k,i,m),j=1,jdim1),k=1,kdim1),i=1,idim1),m=1,2)
                {
                    FILE* fp2 = fortran_get_unit(2);
                    int32_t rb = (int32_t)(2*jdim1*kdim1*idim1*sizeof(double));
                    write_rec_mark(fp2, rb);
                    for (int m=1;m<=2;m++) for (int i=1;i<=idim1;i++)
                        for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++){
                            double v=tursav(j,k,i,m); fwrite(&v,sizeof(double),1,fp2);}
                    write_rec_mark(fp2, rb);
                }
                if (ivisc[0]>=30 || ivisc[1]>=30 || ivisc[2]>=30) {
                    // write(2) nummem
                    {
                        FILE* fp2 = fortran_get_unit(2);
                        int32_t rb = (int32_t)sizeof(int32_t);
                        write_rec_mark(fp2, rb);
                        fwrite(&nummem, sizeof(int32_t), 1, fp2);
                        write_rec_mark(fp2, rb);
                    }
                    // write(2) ((((tursav(j,k,i,m),...),m=3,nummem)
                    {
                        FILE* fp2 = fortran_get_unit(2);
                        int cnt = (nummem >= 3) ? (nummem - 3 + 1) : 0;
                        int32_t rb = (int32_t)(cnt*jdim1*kdim1*idim1*sizeof(double));
                        write_rec_mark(fp2, rb);
                        for (int m=3;m<=nummem;m++) for (int i=1;i<=idim1;i++)
                            for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++){
                                double v=tursav(j,k,i,m); fwrite(&v,sizeof(double),1,fp2);}
                        write_rec_mark(fp2, rb);
                    }
                }
                if (iwghost != 0) {
                    // write(2) ti0,tj0,tk0 (l=1,2)
                    {
                        FILE* fp2 = fortran_get_unit(2);
                        int32_t rb = (int32_t)(
                            (4*2*jdim*kdim + 4*2*kdim*idim1 + 4*2*jdim*idim1)*sizeof(double));
                        write_rec_mark(fp2, rb);
                        // ti0(j,k,l,m): j=1,jdim; k=1,kdim; l=1,2; m=1,4
                        for (int m=1;m<=4;m++) for (int l=1;l<=2;l++)
                            for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++){
                                double v=ti0(j,k,l,m); fwrite(&v,sizeof(double),1,fp2);}
                        // tj0(k,i,l,m): k=1,kdim; i=1,idim1; l=1,2; m=1,4
                        for (int m=1;m<=4;m++) for (int l=1;l<=2;l++)
                            for (int i=1;i<=idim1;i++) for (int k=1;k<=kdim;k++){
                                double v=tj0(k,i,l,m); fwrite(&v,sizeof(double),1,fp2);}
                        // tk0(j,i,l,m): j=1,jdim; i=1,idim1; l=1,2; m=1,4
                        for (int m=1;m<=4;m++) for (int l=1;l<=2;l++)
                            for (int i=1;i<=idim1;i++) for (int j=1;j<=jdim;j++){
                                double v=tk0(j,i,l,m); fwrite(&v,sizeof(double),1,fp2);}
                        write_rec_mark(fp2, rb);
                    }
                    if (ivisc[0]>=30 || ivisc[1]>=30 || ivisc[2]>=30) {
                        // write(2) ti0,tj0,tk0 (l=3,nummem)
                        FILE* fp2 = fortran_get_unit(2);
                        int lcount = (nummem >= 3) ? (nummem - 3 + 1) : 0;
                        int32_t rb = (int32_t)(
                            (4*lcount*jdim*kdim + 4*lcount*kdim*idim1 + 4*lcount*jdim*idim1)*sizeof(double));
                        write_rec_mark(fp2, rb);
                        // ti0(j,k,l,m): j=1,jdim; k=1,kdim; l=3,nummem; m=1,4
                        for (int m=1;m<=4;m++) for (int l=3;l<=nummem;l++)
                            for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++){
                                double v=ti0(j,k,l,m); fwrite(&v,sizeof(double),1,fp2);}
                        // tj0(k,i,l,m): k=1,kdim; i=1,idim1; l=3,nummem; m=1,4
                        for (int m=1;m<=4;m++) for (int l=3;l<=nummem;l++)
                            for (int i=1;i<=idim1;i++) for (int k=1;k<=kdim;k++){
                                double v=tj0(k,i,l,m); fwrite(&v,sizeof(double),1,fp2);}
                        // tk0(j,i,l,m): j=1,jdim; i=1,idim1; l=3,nummem; m=1,4
                        for (int m=1;m<=4;m++) for (int l=3;l<=nummem;l++)
                            for (int i=1;i<=idim1;i++) for (int j=1;j<=jdim;j++){
                                double v=tk0(j,i,l,m); fwrite(&v,sizeof(double),1,fp2);}
                        write_rec_mark(fp2, rb);
                    }
                } // iwghost != 0
            } // icgns != 1
        } // myid == myhost

        // minimum distance function
        if (myid == myhost) {
            if (icgns != 1) {
                // write(2) (((ccabs(smin(j,k,i)),j=1,jdim1),k=1,kdim1),i=1,idim1)
                FILE* fp2 = fortran_get_unit(2);
                int32_t rb = (int32_t)(jdim1*kdim1*idim1*sizeof(double));
                write_rec_mark(fp2, rb);
                for (int i=1;i<=idim1;i++) for (int k=1;k<=kdim1;k++)
                    for (int j=1;j<=jdim1;j++){
                        double v=std::abs(smin(j,k,i)); fwrite(&v,sizeof(double),1,fp2);}
                write_rec_mark(fp2, rb);
            }
        }
    } // ivisc >= 4


    // ivisc == 4 or 25: additional smin-related data (xjb, xkb, blnum)
    if (ivisc[0]==4  || ivisc[1]==4  || ivisc[2]==4  ||
        ivisc[0]==25 || ivisc[1]==25 || ivisc[2]==25) {
        if (myid == myhost) {
            if (icgns != 1) {
                // write(2) (((xjb(j,k,i),j=1,jdim-1),k=1,kdim-1),i=1,idim-1)
                FILE* fp2 = fortran_get_unit(2);
                int32_t rb = (int32_t)(jdim1*kdim1*idim1*sizeof(double));
                write_rec_mark(fp2, rb);
                for (int i=1;i<=idim-1;i++) for (int k=1;k<=kdim-1;k++)
                    for (int j=1;j<=jdim-1;j++){
                        double v=xjb(j,k,i); fwrite(&v,sizeof(double),1,fp2);}
                write_rec_mark(fp2, rb);
            }
        }
        if (myid == myhost) {
            if (icgns != 1) {
                // write(2) (((xkb(j,k,i),j=1,jdim-1),k=1,kdim-1),i=1,idim-1)
                FILE* fp2 = fortran_get_unit(2);
                int32_t rb = (int32_t)(jdim1*kdim1*idim1*sizeof(double));
                write_rec_mark(fp2, rb);
                for (int i=1;i<=idim-1;i++) for (int k=1;k<=kdim-1;k++)
                    for (int j=1;j<=jdim-1;j++){
                        double v=xkb(j,k,i); fwrite(&v,sizeof(double),1,fp2);}
                write_rec_mark(fp2, rb);
            }
        }
        if (myid == myhost) {
            if (icgns != 1) {
                // write(2) (((blnum(j,k,i),j=1,jdim-1),k=1,kdim-1),i=1,idim-1)
                FILE* fp2 = fortran_get_unit(2);
                int32_t rb = (int32_t)(jdim1*kdim1*idim1*sizeof(double));
                write_rec_mark(fp2, rb);
                for (int i=1;i<=idim-1;i++) for (int k=1;k<=kdim-1;k++)
                    for (int j=1;j<=jdim-1;j++){
                        double v=blnum(j,k,i); fwrite(&v,sizeof(double),1,fp2);}
                write_rec_mark(fp2, rb);
            }
        }
    }

    // ivisc == 8,9,13,14: write cmuv
    if (ivisc[0]==8  || ivisc[1]==8  || ivisc[2]==8  ||
        ivisc[0]==9  || ivisc[1]==9  || ivisc[2]==9  ||
        ivisc[0]==13 || ivisc[1]==13 || ivisc[2]==13 ||
        ivisc[0]==14 || ivisc[1]==14 || ivisc[2]==14) {
        if (myid == myhost) {
            if (icgns != 1) {
                // write(2) (((cmuv(j,k,i),j=1,jdim-1),k=1,kdim-1),i=1,idim-1)
                FILE* fp2 = fortran_get_unit(2);
                int32_t rb = (int32_t)(jdim1*kdim1*idim1*sizeof(double));
                write_rec_mark(fp2, rb);
                for (int i=1;i<=idim-1;i++) for (int k=1;k<=kdim-1;k++)
                    for (int j=1;j<=jdim-1;j++){
                        double v=cmuv(j,k,i); fwrite(&v,sizeof(double),1,fp2);}
                write_rec_mark(fp2, rb);
            }
        }
    }

    // icgns == 1 and myid == myhost: empty block in Fortran (no-op)
    // if (icgns == 1 && myid == myhost) { /* empty */ }


    // iteravg == 1 or 2: write averaged data
    if (iteravg == 1 || iteravg == 2) {
        if (ipertavg == 0) {
            // write cell-centered plot3d files
            if (myid == myhost) {
                if (ialph == 0) {
                    // write(96) cell-centered x,y,z + blank(int)
                    FILE* fp96 = fortran_get_unit(96);
                    int32_t rb = (int32_t)(3*jdim1*kdim1*idim1*sizeof(double)
                                          + jdim1*kdim1*idim1*sizeof(int32_t));
                    write_rec_mark(fp96, rb);
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){
                            double v=0.125*(x(j,k,i)+x(j+1,k,i)+x(j,k+1,i)+x(j,k,i+1)
                                +x(j+1,k+1,i)+x(j+1,k,i+1)+x(j,k+1,i+1)+x(j+1,k+1,i+1));
                            fwrite(&v,sizeof(double),1,fp96);}
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){
                            double v=0.125*(y(j,k,i)+y(j+1,k,i)+y(j,k+1,i)+y(j,k,i+1)
                                +y(j+1,k+1,i)+y(j+1,k,i+1)+y(j,k+1,i+1)+y(j+1,k+1,i+1));
                            fwrite(&v,sizeof(double),1,fp96);}
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){
                            double v=0.125*(z(j,k,i)+z(j+1,k,i)+z(j,k+1,i)+z(j,k,i+1)
                                +z(j+1,k+1,i)+z(j+1,k,i+1)+z(j,k+1,i+1)+z(j+1,k+1,i+1));
                            fwrite(&v,sizeof(double),1,fp96);}
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){
                            int32_t v=(int32_t)blank(j,k,i); fwrite(&v,sizeof(int32_t),1,fp96);}
                    write_rec_mark(fp96, rb);
                } else {
                    // ialph != 0: write x, z, -y
                    FILE* fp96 = fortran_get_unit(96);
                    int32_t rb = (int32_t)(3*jdim1*kdim1*idim1*sizeof(double)
                                          + jdim1*kdim1*idim1*sizeof(int32_t));
                    write_rec_mark(fp96, rb);
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){
                            double v=0.125*(x(j,k,i)+x(j+1,k,i)+x(j,k+1,i)+x(j,k,i+1)
                                +x(j+1,k+1,i)+x(j+1,k,i+1)+x(j,k+1,i+1)+x(j+1,k+1,i+1));
                            fwrite(&v,sizeof(double),1,fp96);}
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){
                            double v=0.125*(z(j,k,i)+z(j+1,k,i)+z(j,k+1,i)+z(j,k,i+1)
                                +z(j+1,k+1,i)+z(j+1,k,i+1)+z(j,k+1,i+1)+z(j+1,k+1,i+1));
                            fwrite(&v,sizeof(double),1,fp96);}
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){
                            double v=-0.125*(y(j,k,i)+y(j+1,k,i)+y(j,k+1,i)+y(j,k,i+1)
                                +y(j+1,k+1,i)+y(j+1,k,i+1)+y(j,k+1,i+1)+y(j+1,k+1,i+1));
                            fwrite(&v,sizeof(double),1,fp96);}
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){
                            int32_t v=(int32_t)blank(j,k,i); fwrite(&v,sizeof(int32_t),1,fp96);}
                    write_rec_mark(fp96, rb);
                }

                // write(97) xmach,alpha,reue,xnumavg
                {
                    FILE* fp97 = fortran_get_unit(97);
                    int32_t rb = (int32_t)(4*sizeof(float));
                    write_rec_mark(fp97, rb);
                    fwrite(&xmach,   sizeof(float),1,fp97);
                    fwrite(&alpha,   sizeof(float),1,fp97);
                    fwrite(&reue,    sizeof(float),1,fp97);
                    fwrite(&xnumavg, sizeof(float),1,fp97);
                    write_rec_mark(fp97, rb);
                }

                if (ialph == 0) {
                    // write(97) ((((qavg(j,k,i,m),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),m=1,5)
                    FILE* fp97 = fortran_get_unit(97);
                    int32_t rb = (int32_t)(5*jdim1*kdim1*idim1*sizeof(double));
                    write_rec_mark(fp97, rb);
                    for (int m=1;m<=5;m++) for (int k=1;k<=kdim1;k++)
                        for (int j=1;j<=jdim1;j++) for (int i=1;i<=idim1;i++){
                            double v=qavg(j,k,i,m); fwrite(&v,sizeof(double),1,fp97);}
                    write_rec_mark(fp97, rb);
                } else {
                    // ialph != 0: components 1,2,4,-3,5
                    FILE* fp97 = fortran_get_unit(97);
                    int32_t rb = (int32_t)(5*jdim1*kdim1*idim1*sizeof(double));
                    write_rec_mark(fp97, rb);
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){double v=qavg(j,k,i,1); fwrite(&v,sizeof(double),1,fp97);}
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){double v=qavg(j,k,i,2); fwrite(&v,sizeof(double),1,fp97);}
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){double v=qavg(j,k,i,4); fwrite(&v,sizeof(double),1,fp97);}
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){double v=-qavg(j,k,i,3); fwrite(&v,sizeof(double),1,fp97);}
                    for (int k=1;k<=kdim1;k++) for (int j=1;j<=jdim1;j++)
                        for (int i=1;i<=idim1;i++){double v=qavg(j,k,i,5); fwrite(&v,sizeof(double),1,fp97);}
                    write_rec_mark(fp97, rb);
                }
            } // myid == myhost (ipertavg==0)
        } else {
            // ipertavg != 0
            if (myid == myhost) {
                if (ialph == 0) {
                    // write(96) x,y,z node coordinates
                    FILE* fp96 = fortran_get_unit(96);
                    int32_t rb = (int32_t)(3*jdim*kdim*idim*sizeof(double));
                    write_rec_mark(fp96, rb);
                    for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                        for (int i=1;i<=idim;i++){double v=x(j,k,i); fwrite(&v,sizeof(double),1,fp96);}
                    for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                        for (int i=1;i<=idim;i++){double v=y(j,k,i); fwrite(&v,sizeof(double),1,fp96);}
                    for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                        for (int i=1;i<=idim;i++){double v=z(j,k,i); fwrite(&v,sizeof(double),1,fp96);}
                    write_rec_mark(fp96, rb);
                } else {
                    // ialph != 0: write x, z, -y
                    FILE* fp96 = fortran_get_unit(96);
                    int32_t rb = (int32_t)(3*jdim*kdim*idim*sizeof(double));
                    write_rec_mark(fp96, rb);
                    for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                        for (int i=1;i<=idim;i++){double v=x(j,k,i); fwrite(&v,sizeof(double),1,fp96);}
                    for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                        for (int i=1;i<=idim;i++){double v=z(j,k,i); fwrite(&v,sizeof(double),1,fp96);}
                    for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                        for (int i=1;i<=idim;i++){double v=-y(j,k,i); fwrite(&v,sizeof(double),1,fp96);}
                    write_rec_mark(fp96, rb);
                }

                // write(97) xmach,alpha,reue,xnumavg
                {
                    FILE* fp97 = fortran_get_unit(97);
                    int32_t rb = (int32_t)(4*sizeof(float));
                    write_rec_mark(fp97, rb);
                    fwrite(&xmach,   sizeof(float),1,fp97);
                    fwrite(&alpha,   sizeof(float),1,fp97);
                    fwrite(&reue,    sizeof(float),1,fp97);
                    fwrite(&xnumavg, sizeof(float),1,fp97);
                    write_rec_mark(fp97, rb);
                }
                // write(98) xmach,alpha,reue,xnumavg2
                {
                    FILE* fp98 = fortran_get_unit(98);
                    int32_t rb = (int32_t)(4*sizeof(float));
                    write_rec_mark(fp98, rb);
                    fwrite(&xmach,    sizeof(float),1,fp98);
                    fwrite(&alpha,    sizeof(float),1,fp98);
                    fwrite(&reue,     sizeof(float),1,fp98);
                    fwrite(&xnumavg2, sizeof(float),1,fp98);
                    write_rec_mark(fp98, rb);
                }
                // write(95) xmach,alpha,reue,xnumavg2
                {
                    FILE* fp95 = fortran_get_unit(95);
                    int32_t rb = (int32_t)(4*sizeof(float));
                    write_rec_mark(fp95, rb);
                    fwrite(&xmach,    sizeof(float),1,fp95);
                    fwrite(&alpha,    sizeof(float),1,fp95);
                    fwrite(&reue,     sizeof(float),1,fp95);
                    fwrite(&xnumavg2, sizeof(float),1,fp95);
                    write_rec_mark(fp95, rb);
                }


                if (ialph == 0) {
                    // write(97) ((((qavg(j,k,i,m),i=1,idim),j=1,jdim),k=1,kdim),m=1,5)
                    {
                        FILE* fp97 = fortran_get_unit(97);
                        int32_t rb = (int32_t)(5*jdim*kdim*idim*sizeof(double));
                        write_rec_mark(fp97, rb);
                        for (int m=1;m<=5;m++) for (int k=1;k<=kdim;k++)
                            for (int j=1;j<=jdim;j++) for (int i=1;i<=idim;i++){
                                double v=qavg(j,k,i,m); fwrite(&v,sizeof(double),1,fp97);}
                        write_rec_mark(fp97, rb);
                    }
                    // write(98) ((((q2avg(j,k,i,m),i=1,idim),j=1,jdim),k=1,kdim),m=1,5)
                    {
                        FILE* fp98 = fortran_get_unit(98);
                        int32_t rb = (int32_t)(5*jdim*kdim*idim*sizeof(double));
                        write_rec_mark(fp98, rb);
                        for (int m=1;m<=5;m++) for (int k=1;k<=kdim;k++)
                            for (int j=1;j<=jdim;j++) for (int i=1;i<=idim;i++){
                                double v=q2avg(j,k,i,m); fwrite(&v,sizeof(double),1,fp98);}
                        write_rec_mark(fp98, rb);
                    }
                    // write(95) ((((q2avg(j,k,i,m)-qavg(j,k,i,m)**2,...),m=1,5)
                    {
                        FILE* fp95 = fortran_get_unit(95);
                        int32_t rb = (int32_t)(5*jdim*kdim*idim*sizeof(double));
                        write_rec_mark(fp95, rb);
                        for (int m=1;m<=5;m++) for (int k=1;k<=kdim;k++)
                            for (int j=1;j<=jdim;j++) for (int i=1;i<=idim;i++){
                                double v=q2avg(j,k,i,m)-qavg(j,k,i,m)*qavg(j,k,i,m);
                                fwrite(&v,sizeof(double),1,fp95);}
                        write_rec_mark(fp95, rb);
                    }
                } else {
                    // ialph != 0
                    // write(97) qavg components 1,2,4,-3,5
                    {
                        FILE* fp97 = fortran_get_unit(97);
                        int32_t rb = (int32_t)(5*jdim*kdim*idim*sizeof(double));
                        write_rec_mark(fp97, rb);
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){double v=qavg(j,k,i,1); fwrite(&v,sizeof(double),1,fp97);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){double v=qavg(j,k,i,2); fwrite(&v,sizeof(double),1,fp97);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){double v=qavg(j,k,i,4); fwrite(&v,sizeof(double),1,fp97);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){double v=-qavg(j,k,i,3); fwrite(&v,sizeof(double),1,fp97);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){double v=qavg(j,k,i,5); fwrite(&v,sizeof(double),1,fp97);}
                        write_rec_mark(fp97, rb);
                    }
                    // write(98) q2avg components 1,2,4,3,5
                    {
                        FILE* fp98 = fortran_get_unit(98);
                        int32_t rb = (int32_t)(5*jdim*kdim*idim*sizeof(double));
                        write_rec_mark(fp98, rb);
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){double v=q2avg(j,k,i,1); fwrite(&v,sizeof(double),1,fp98);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){double v=q2avg(j,k,i,2); fwrite(&v,sizeof(double),1,fp98);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){double v=q2avg(j,k,i,4); fwrite(&v,sizeof(double),1,fp98);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){double v=q2avg(j,k,i,3); fwrite(&v,sizeof(double),1,fp98);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){double v=q2avg(j,k,i,5); fwrite(&v,sizeof(double),1,fp98);}
                        write_rec_mark(fp98, rb);
                    }
                    // write(95) q2avg-qavg**2 components 1,2,4,3,5
                    {
                        FILE* fp95 = fortran_get_unit(95);
                        int32_t rb = (int32_t)(5*jdim*kdim*idim*sizeof(double));
                        write_rec_mark(fp95, rb);
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){
                                double v=q2avg(j,k,i,1)-qavg(j,k,i,1)*qavg(j,k,i,1);
                                fwrite(&v,sizeof(double),1,fp95);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){
                                double v=q2avg(j,k,i,2)-qavg(j,k,i,2)*qavg(j,k,i,2);
                                fwrite(&v,sizeof(double),1,fp95);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){
                                double v=q2avg(j,k,i,4)-qavg(j,k,i,4)*qavg(j,k,i,4);
                                fwrite(&v,sizeof(double),1,fp95);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){
                                double v=q2avg(j,k,i,3)-qavg(j,k,i,3)*qavg(j,k,i,3);
                                fwrite(&v,sizeof(double),1,fp95);}
                        for (int k=1;k<=kdim;k++) for (int j=1;j<=jdim;j++)
                            for (int i=1;i<=idim;i++){
                                double v=q2avg(j,k,i,5)-qavg(j,k,i,5)*qavg(j,k,i,5);
                                fwrite(&v,sizeof(double),1,fp95);}
                        write_rec_mark(fp95, rb);
                    }
                } // end ialph != 0
            } // end myid == myhost (ipertavg != 0)
        } // end ipertavg != 0
    } // end iteravg == 1 || 2

    // clcd information
    if (myid == myhost && iskip > 0) {
        if (iclcd == 1 || iclcd == 2) {
            // write(102) n_clcd,nblocks_clcd,ntt
            {
                FILE* fp102 = fortran_get_unit(102);
                int32_t rb = (int32_t)(3*sizeof(int32_t));
                write_rec_mark(fp102, rb);
                fwrite(&n_clcd,       sizeof(int32_t),1,fp102);
                fwrite(&nblocks_clcd, sizeof(int32_t),1,fp102);
                fwrite(&ntt,          sizeof(int32_t),1,fp102);
                write_rec_mark(fp102, rb);
            }
            // write(102) clcd(1:2,1:n_clcd,1:ntt)
            {
                FILE* fp102 = fortran_get_unit(102);
                int32_t rb = (int32_t)(2*n_clcd*ntt*sizeof(double));
                write_rec_mark(fp102, rb);
                for (int n=1;n<=ntt;n++) for (int j=1;j<=n_clcd;j++)
                    for (int l=1;l<=2;l++){
                        double v=clcd(l,j,n); fwrite(&v,sizeof(double),1,fp102);}
                write_rec_mark(fp102, rb);
            }
        }
    }

    return;
} // end wrest


} // namespace wrest_ns
