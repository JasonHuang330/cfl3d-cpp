// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rrest.h"
#include "termn8.h"
#include "init_mast.h"
#include "umalloc.h"
#include "ccomplex.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>

namespace rrest_ns {

// Helper: read one Fortran unformatted record: skip leading marker, read data, skip trailing marker
static inline void read_unformatted_raw(FILE* f, void* buf, size_t nbytes) {
    int32_t reclen;
    fread(&reclen, 4, 1, f);
    fread(buf, 1, nbytes, f);
    fread(&reclen, 4, 1, f);
}

// Helper: skip an entire Fortran unformatted record (read length, skip data, read trailing length)
static inline void skip_unformatted_record(FILE* f) {
    int32_t reclen;
    fread(&reclen, 4, 1, f);
    fseek(f, reclen, SEEK_CUR);
    fread(&reclen, 4, 1, f);
}

void rrest(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, int& ncycmax, int& ntr, FortranArray1DRef<double> rms, FortranArray1DRef<double> clw, FortranArray1DRef<double> cdw, FortranArray1DRef<double> cdpw, FortranArray1DRef<double> cdvw, FortranArray1DRef<double> cxw, FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw, FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw, FortranArray1DRef<double> cmzw, int& n_clcd, FortranArray3DRef<double> clcd, int& nblocks_clcd, FortranArray2DRef<int> blocks_clcd, FortranArray1DRef<double> fmdotw, FortranArray1DRef<double> cftmomw, FortranArray1DRef<double> cftpw, FortranArray1DRef<double> cftvw, FortranArray1DRef<double> cfttotw, FortranArray2DRef<double> rmstr, FortranArray2DRef<int> nneg, int& iskip, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> tursav, FortranArray3DRef<double> smin, FortranArray3DRef<double> xjb, FortranArray3DRef<double> xkb, FortranArray3DRef<double> blnum, FortranArray3DRef<double> cmuv, int& maxbl, FortranArray1DRef<int> mblk2nd, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& igrid, FortranArray3DRef<double> wk, int& idima, int& jdima, int& kdima, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, FortranArray4DRef<double> qavg, FortranArray4DRef<double> q2avg, int& nummem)
{
    // COMMON block aliases
    int32_t& icgns       = cmn_cgns.icgns;

    float& xnumavg       = cmn_avgdata.xnumavg;
    int32_t& iteravg     = cmn_avgdata.iteravg;
    float& xnumavg2      = cmn_avgdata.xnumavg2;
    int32_t& ipertavg    = cmn_avgdata.ipertavg;
    int32_t& iclcd       = cmn_avgdata.iclcd;

    float& xmach         = cmn_info.xmach;
    float& alpha         = cmn_info.alpha;
    float& dt            = cmn_info.dt;
    int32_t& ntt         = cmn_info.ntt;

    int32_t& ivmx        = cmn_maxiv.ivmx;

    float& time          = cmn_unst.time;

    int32_t& irest2      = cmn_wrestq.irest2;

    float& reue          = cmn_reyue.reue;
    int32_t* ivisc       = cmn_reyue.ivisc;  // ivisc[0..2] = ivisc(1..3)

    int32_t& isminc      = cmn_sminn.isminc;
    int32_t& ismincforce = cmn_sminn.ismincforce;

    int32_t& ialphit     = cmn_alphait.ialphit;
    float& radtodeg      = cmn_conversion.radtodeg;

    int32_t& irghost     = cmn_ghost.irghost;
    int32_t& ialph       = cmn_igrdtyp.ialph;

    float& randomize     = cmn_random_input.randomize;

    // Local variables
    int idim1, jdim1, kdim1;
    float titlw[20];
    float xmachw = 0.0f, alphw = 0.0f, reuew = 0.0f;
    int jt = 0, kt = 0, it = 0;
    float dum = 0.0f;
    int idum = 0;
    int iv1 = 0, iv2 = 0, iv3 = 0;
    int ivmxold;
    int iread;
    int irinddata[6];
    int n;
    int irnd_count;
    float rnd_nbr;
    int jki_size;
    int stats_var = 0;
    int nummem_read = 0;
    double temp;
    int i_ntt = 0;
    int i_clcd_r = 0;
    int iblocks_clcd_r = 0;
    int memuse = 0;
    std::vector<float> harvest;

    FILE* f2 = fortran_get_unit(2);

    idim1 = idim - 1;
    jdim1 = jdim - 1;
    kdim1 = kdim - 1;



    // Lines 102-114: read header record from unit 2 (if not CGNS)
    if (icgns == 1) {
        // CGNS path: no binary read here
    } else {
        // read(2) titlw, xmachw, jt, kt, it, alphw, reuew, ntr, time
        {
            int32_t reclen;
            fread(&reclen, 4, 1, f2);
            fread(titlw, sizeof(float), 20, f2);
            fread(&xmachw, sizeof(float), 1, f2);
            fread(&jt, sizeof(int), 1, f2);
            fread(&kt, sizeof(int), 1, f2);
            fread(&it, sizeof(int), 1, f2);
            fread(&alphw, sizeof(float), 1, f2);
            fread(&reuew, sizeof(float), 1, f2);
            fread(&ntr, sizeof(int), 1, f2);
            fread(&time, sizeof(float), 1, f2);
            fread(&reclen, 4, 1, f2);
        }
        if (jt != jdim || kt != kdim || it != idim) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, " stopping.inconsistent restart grid indices\n");
            fprintf(f11, " block indices  idim,jdim,kdim = %d %d %d\n", idim, jdim, kdim);
            fprintf(f11, " restart indices  idim,jdim,kdim = %d %d %d\n", it, jt, kt);
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }

    // Line 120: ntt = ntr
    ntt = ntr;

    // Lines 124-128: alpha iteration
    if (ialphit > 0) {
        alpha = alphw / radtodeg;
        // reset freestream for new alpha on host
        init_mast_ns::init_mast();
    }

    // Lines 130-172: convergence data read
    if (icgns != 1) {
        if (ntr <= ncycmax) {
            // Convergence data (residual, force coefficients, mass flow, etc.)
            if (iskip > 0) {
                // read(2) (rms(n),n=1,ntr),(clw(n),n=1,ntr),...
                {
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); rms(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); clw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cdw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cdpw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cdvw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cxw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cyw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); czw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cmxw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cmyw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cmzw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); fmdotw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cftmomw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cftpw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cftvw(nn) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); cfttotw(nn) = (double)v; }
                    fread(&reclen, 4, 1, f2);
                }
            }
        } else {
            if (irest2 == 1) {
                if (iskip > 0) {
                    // read(2) (dum,n=1,ntr)*16 - skip 16*ntr floats in one record
                    {
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        for (int nn = 0; nn < 16 * ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                        fread(&reclen, 4, 1, f2);
                    }
                }
            } else {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, "\n stopping...");
                fprintf(f11, "previous number of iterations computed > dimension ncycmax\n");
                fprintf(f11, " ntr,ncycmax =  %d %d\n", ntr, ncycmax);
                fprintf(f11, " increase value of ncycmax to at LEAST  %d\n", ntr + ncycmax);
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
        }
    }



    // Lines 176-184: read primitive variables q
    if (icgns != 1) {
        // read(2) ((((q(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,5)
        // q is stored column-major: j varies fastest, then k, then i, then l
        {
            int32_t reclen;
            fread(&reclen, 4, 1, f2);
            for (int l = 1; l <= 5; l++) {
                for (int i = 1; i <= idim1; i++) {
                    for (int k = 1; k <= kdim1; k++) {
                        for (int j = 1; j <= jdim1; j++) {
                            float v;
                            fread(&v, sizeof(float), 1, f2);
                            q(j, k, i, l) = (double)v;
                        }
                    }
                }
            }
            fread(&reclen, 4, 1, f2);
        }
        if (irghost != 0) {
            // read(2) ((((qi0(j,k,l,m),j=1,jdim1),k=1,kdim1),l=1,5),m=1,4),
            //         ((((qj0(k,i,l,m),k=1,kdim1),i=1,idim1),l=1,5),m=1,4),
            //         ((((qk0(j,i,l,m),j=1,jdim1),i=1,idim1),l=1,5),m=1,4)
            {
                int32_t reclen;
                fread(&reclen, 4, 1, f2);
                for (int m = 1; m <= 4; m++) {
                    for (int l = 1; l <= 5; l++) {
                        for (int k = 1; k <= kdim1; k++) {
                            for (int j = 1; j <= jdim1; j++) {
                                float v;
                                fread(&v, sizeof(float), 1, f2);
                                qi0(j, k, l, m) = (double)v;
                            }
                        }
                    }
                }
                for (int m = 1; m <= 4; m++) {
                    for (int l = 1; l <= 5; l++) {
                        for (int i = 1; i <= idim1; i++) {
                            for (int k = 1; k <= kdim1; k++) {
                                float v;
                                fread(&v, sizeof(float), 1, f2);
                                qj0(k, i, l, m) = (double)v;
                            }
                        }
                    }
                }
                for (int m = 1; m <= 4; m++) {
                    for (int l = 1; l <= 5; l++) {
                        for (int i = 1; i <= idim1; i++) {
                            for (int j = 1; j <= jdim1; j++) {
                                float v;
                                fread(&v, sizeof(float), 1, f2);
                                qk0(j, i, l, m) = (double)v;
                            }
                        }
                    }
                }
                fread(&reclen, 4, 1, f2);
            }
        }
    }

    // Lines 186-236: Mach number fix-up
    n = jdim * kdim;
    if (xmachw != xmach) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::memset(bou(nou(1), 1), ' ', 120);
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "WARNING: fix-up taken for Mach number inconsistancy in restart file");
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "  scaling u,v,w by M(new)/M(old); cp held fixed: M(new), M(old)= %12.4e%12.4e",
            (float)xmach, (float)xmachw);
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::memset(bou(nou(1), 1), ' ', 120);
        for (int i = 1; i <= idim; i++) {
            // u/cref, v/cref, w/cref
            for (int izz = 1; izz <= n; izz++) {
                q(izz, 1, i, 2) = q(izz, 1, i, 2) * xmach / xmachw;
                q(izz, 1, i, 3) = q(izz, 1, i, 3) * xmach / xmachw;
                q(izz, 1, i, 4) = q(izz, 1, i, 4) * xmach / xmachw;
                // cp=constant
                q(izz, 1, i, 5) = q(izz, 1, i, 5) * (double)cmn_fluid.gamma;
                q(izz, 1, i, 5) = (q(izz, 1, i, 5) - 1.0) * (double)xmach *
                                   (double)xmach / (double)xmachw / (double)xmachw;
                q(izz, 1, i, 5) = (q(izz, 1, i, 5) + 1.0) / (double)cmn_fluid.gamma;
            }
            // check for negative density, pressure
            for (int j = 1; j <= jdim; j++) {
                for (int k = 1; k <= kdim; k++) {
                    if ((float)q(j, k, i, 1) <= 0.0f || (float)q(j, k, i, 5) <= 0.0f) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120, " neg. density,pressure");
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120, " %g %g",
                            (float)q(j, k, i, 1), (float)q(j, k, i, 5));
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120, " %d %d %d", i, j, k);
                        int m1 = -1;
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }
    }

    // Lines 237-249: ReUe inconsistency warning
    if (reuew != reue) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::memset(bou(nou(1), 1), ' ', 120);
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "WARNING:  ReUe inconsistency in restart file data Solution will proceed with new value!");
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "  ReUe(new), ReUe(old)= %12.4e%12.4e", (float)reue, (float)reuew);
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::memset(bou(nou(1), 1), ' ', 120);
    }



    // Lines 251-274: randomize
    if ((float)randomize > 1.e-12f) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, " adding random component to restart values");
        jki_size = (jdim - 1) * (kdim - 1) * (idim - 1);
        harvest.resize(jki_size);
        {
            char text1[1] = {' '};
            umalloc_ns::umalloc(jki_size, stats_var, text1, memuse, stats_var);
        }
        // call random_seed; call random_number(harvest)
        // Use C++ random number generation
        {
            unsigned int seed = 12345;
            srand(seed);
            for (int idx = 0; idx < jki_size; idx++) {
                harvest[idx] = (float)rand() / (float)RAND_MAX;
            }
        }
        irnd_count = 0;
        for (int i = 1; i <= idim - 1; i++) {
            for (int k = 1; k <= kdim - 1; k++) {
                for (int j = 1; j <= jdim - 1; j++) {
                    irnd_count = irnd_count + 1;
                    rnd_nbr = harvest[irnd_count - 1] - 0.5f;
                    q(j, k, i, 1) = q(j, k, i, 1) * (1.0 + (double)rnd_nbr * (double)randomize);
                    q(j, k, i, 2) = q(j, k, i, 2) * (1.0 + (double)rnd_nbr * (double)randomize);
                    q(j, k, i, 3) = q(j, k, i, 3) * (1.0 + (double)rnd_nbr * (double)randomize);
                    q(j, k, i, 4) = q(j, k, i, 4) * (1.0 + (double)rnd_nbr * (double)randomize);
                    q(j, k, i, 5) = q(j, k, i, 5) * (1.0 + (double)rnd_nbr * (double)randomize);
                }
            }
        }
    }

    // Lines 279-282: read turbulence viscosity flags
    if (icgns == 1) {
        // CGNS path: no binary read
    } else {
        // read(2) iv1, iv2, iv3
        {
            int32_t reclen;
            fread(&reclen, 4, 1, f2);
            fread(&iv1, sizeof(int), 1, f2);
            fread(&iv2, sizeof(int), 1, f2);
            fread(&iv3, sizeof(int), 1, f2);
            fread(&reclen, 4, 1, f2);
        }
    }

    // Lines 283-292: ivisc inconsistency warning
    if (iv1 != ivisc[0] || iv2 != ivisc[1] || iv3 != ivisc[2]) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::memset(bou(nou(1), 1), ' ', 120);
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            " WARNING: restarting block%6d with DIFFERENT ivisc: old=%3d%3d%3d, new=%3d%3d%3d",
            nbl, iv1, iv2, iv3, ivisc[0], ivisc[1], ivisc[2]);
    }

    // Lines 293-300: determine iread
    ivmxold = std::max(iv1, iv2);
    ivmxold = std::max(ivmxold, iv3);
    iread = 1;
    if (ivmxold != ivmx) iread = 0;
    if ((ivmxold >= 6 && ivmxold <= 15) && (ivmx >= 6 && ivmx <= 15)) iread = 1;
    if ((ivmxold == 4 && ivmx == 25) || (ivmxold == 25 && ivmx == 4)) iread = 1;

    // Lines 304-325: time-accurate ivisc=30/40 checks
    if ((float)dt > 0.0f) {
        if (ivmxold == 30 && ivmx != 30) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, " old ivisc=30; new must also be... stopping\n");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (ivmxold != 30 && ivmx == 30) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, " new ivisc=30; old must also be... stopping\n");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }
    if ((float)dt > 0.0f) {
        if (ivmxold == 40 && ivmx != 40) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, " old ivisc=40; new must also be... stopping\n");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (ivmxold != 40 && ivmx == 40) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, " new ivisc=40; old must also be... stopping\n");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }

    // Lines 327-347: special checks for ivisc=72
    if (iread == 0) {
        if (ivmx == 72) {
            if (ivmxold != 6 && ivmxold != 7 && ivmxold != 8 &&
                ivmxold != 12 && ivmxold != 14) {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, "cannot restart ivisc=72 from turb model %d\n", ivmxold);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " ERROR: cannot restart ivisc=72 from turb model%6d", ivmxold);
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
        }
        if (ivmxold == 72 && ivmx != 72) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, "cannot restart from turb model 72\n");
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                " ERROR: cannot restart from turb model 72");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }



    // Lines 351-375: turbulence model convergence data
    if (iskip > 0) {
        if (irest2 == 1) {
            if (icgns != 1) {
                // read(2) (dum,n=1,ntr),(dum,n=1,ntr),(idum,n=1,ntr),(idum,n=1,ntr)
                {
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    for (int nn = 0; nn < ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                    for (int nn = 0; nn < ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                    for (int nn = 0; nn < ntr; nn++) { int v; fread(&v, sizeof(int), 1, f2); }
                    for (int nn = 0; nn < ntr; nn++) { int v; fread(&v, sizeof(int), 1, f2); }
                    fread(&reclen, 4, 1, f2);
                }
                if (iv1 >= 30 || iv2 >= 30 || iv3 >= 30) {
                    // read(2) nummem_read
                    {
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        fread(&nummem_read, sizeof(int), 1, f2);
                        fread(&reclen, 4, 1, f2);
                    }
                    for (int l = 3; l <= nummem_read; l++) {
                        // read(2) (dum,n=1,ntt),(idum,n=1,ntt)
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        for (int nn = 0; nn < ntt; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                        for (int nn = 0; nn < ntt; nn++) { int v; fread(&v, sizeof(int), 1, f2); }
                        fread(&reclen, 4, 1, f2);
                    }
                }
            }
        } else {
            if (icgns != 1) {
                // read(2) (rmstr(n,1),n=1,ntr),(rmstr(n,2),n=1,ntr),
                //         (nneg(n,1),n=1,ntr),(nneg(n,2),n=1,ntr)
                {
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); rmstr(nn, 1) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { float v; fread(&v, sizeof(float), 1, f2); rmstr(nn, 2) = (double)v; }
                    for (int nn = 1; nn <= ntr; nn++) { int v; fread(&v, sizeof(int), 1, f2); nneg(nn, 1) = v; }
                    for (int nn = 1; nn <= ntr; nn++) { int v; fread(&v, sizeof(int), 1, f2); nneg(nn, 2) = v; }
                    fread(&reclen, 4, 1, f2);
                }
                if (iv1 >= 30 || iv2 >= 30 || iv3 >= 30) {
                    // read(2) nummem_read
                    {
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        fread(&nummem_read, sizeof(int), 1, f2);
                        fread(&reclen, 4, 1, f2);
                    }
                    for (int l = 3; l <= nummem_read; l++) {
                        // read(2) (rmstr(n,l),n=1,ntt),(nneg(n,l),n=1,ntt)
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        for (int nn = 1; nn <= ntt; nn++) { float v; fread(&v, sizeof(float), 1, f2); rmstr(nn, l) = (double)v; }
                        for (int nn = 1; nn <= ntt; nn++) { int v; fread(&v, sizeof(int), 1, f2); nneg(nn, l) = v; }
                        fread(&reclen, 4, 1, f2);
                    }
                }
            }
        }
    }

    // Lines 377-378: CGNS + myid==myhost (empty block)
    if (icgns == 1 && myid == myhost) {
        // nothing
    }



    // Lines 379-674: iread == 0 branch (different turb model)
    if (iread == 0) {
        // Lines 392-440: skip turb data from file
        if (iv1 >= 2 || iv2 >= 2 || iv3 >= 2) {
            if (icgns != 1) {
                // read(2) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
                {
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    int cnt = jdim1 * kdim1 * idim1;
                    for (int nn = 0; nn < cnt; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                    fread(&reclen, 4, 1, f2);
                }
                if (irghost != 0) {
                    // read(2) ((((dum,j=1,jdim),k=1,kdim),l=1,1),m=1,4),
                    //         ((((dum,k=1,kdim),i=1,idim1),l=1,1),m=1,4),
                    //         ((((dum,j=1,jdim),i=1,idim1),l=1,1),m=1,4)
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    int cnt1 = jdim * kdim * 1 * 4;
                    int cnt2 = kdim * idim1 * 1 * 4;
                    int cnt3 = jdim * idim1 * 1 * 4;
                    for (int nn = 0; nn < cnt1; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                    for (int nn = 0; nn < cnt2; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                    for (int nn = 0; nn < cnt3; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                    fread(&reclen, 4, 1, f2);
                }
            }
        }

        if (iv1 >= 4 || iv2 >= 4 || iv3 >= 4) {
            if (icgns != 1) {
                // read(2) ((((dum,j=1,jdim1),k=1,kdim1),i=1,idim1),m=1,2)
                {
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    int cnt = jdim1 * kdim1 * idim1 * 2;
                    for (int nn = 0; nn < cnt; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                    fread(&reclen, 4, 1, f2);
                }
                if (irghost != 0) {
                    // read(2) ((((dum,j=1,jdim),k=1,kdim),l=1,2),m=1,4),
                    //         ((((dum,k=1,kdim),i=1,idim1),l=1,2),m=1,4),
                    //         ((((dum,j=1,jdim),i=1,idim1),l=1,2),m=1,4)
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    int cnt1 = jdim * kdim * 2 * 4;
                    int cnt2 = kdim * idim1 * 2 * 4;
                    int cnt3 = jdim * idim1 * 2 * 4;
                    for (int nn = 0; nn < cnt1; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                    for (int nn = 0; nn < cnt2; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                    for (int nn = 0; nn < cnt3; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                    fread(&reclen, 4, 1, f2);
                }
                // read(2) (((smin(j,k,i),j=1,jdim1),k=1,kdim1),i=1,idim1)
                {
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    for (int i = 1; i <= idim1; i++) {
                        for (int k = 1; k <= kdim1; k++) {
                            for (int j = 1; j <= jdim1; j++) {
                                float v;
                                fread(&v, sizeof(float), 1, f2);
                                smin(j, k, i) = (double)v;
                            }
                        }
                    }
                    fread(&reclen, 4, 1, f2);
                }
            }
            if (iv1 == 4 || iv2 == 4 || iv3 == 4 ||
                iv1 == 25 || iv2 == 25 || iv3 == 25) {
                if (icgns != 1) {
                    // read(2) (((dum,...))) x3
                    for (int rr = 0; rr < 3; rr++) {
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        int cnt = jdim1 * kdim1 * idim1;
                        for (int nn = 0; nn < cnt; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                        fread(&reclen, 4, 1, f2);
                    }
                }
            }
            if (iv1 == 8 || iv2 == 8 || iv3 == 8 ||
                iv1 == 9 || iv2 == 9 || iv3 == 9 ||
                iv1 == 13 || iv2 == 13 || iv3 == 13 ||
                iv1 == 14 || iv2 == 14 || iv3 == 14) {
                if (icgns != 1) {
                    // read(2) (((dum,...)))
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    int cnt = jdim1 * kdim1 * idim1;
                    for (int nn = 0; nn < cnt; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                    fread(&reclen, 4, 1, f2);
                }
            }
            // Lines 431-443: isminc logic
            if ((ivisc[0] >= 5 || ivisc[1] >= 5 || ivisc[2] >= 5) &&
                isminc != 2 && (ivmxold != 4 && ivmxold != 25)) isminc = 3;
            if ((ivisc[0] == 25 || ivisc[1] == 25 || ivisc[2] == 25) &&
                (ivmxold != 4 && ivmxold != 25)) isminc = 1;
            if ((ivisc[0] == 25 || ivisc[1] == 25 || ivisc[2] == 25) &&
                ivmxold == 4) isminc = 3;
            if ((ivisc[0] == 4 || ivisc[1] == 4 || ivisc[2] == 4) &&
                ivmxold == 25) isminc = 3;
        }
        if (isminc == 2) {
            isminc = 1;
        }
        if (ismincforce != -1) isminc = ismincforce;
    } else {
        // Lines 446-674: iread != 0 branch (same turb model)
        if (isminc != 2) {
            isminc = 0;
        }
        if (ismincforce != -1) isminc = ismincforce;

        if (iv1 >= 2 || iv2 >= 2 || iv3 >= 2) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                " reading vist3d data from  restart file, block%5d", nbl);
            if (icgns != 1) {
                // read(2) (((vist3d(j,k,i),j=1,jdim1),k=1,kdim1),i=1,idim1)
                {
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    for (int i = 1; i <= idim1; i++) {
                        for (int k = 1; k <= kdim1; k++) {
                            for (int j = 1; j <= jdim1; j++) {
                                float v;
                                fread(&v, sizeof(float), 1, f2);
                                vist3d(j, k, i) = (double)v;
                            }
                        }
                    }
                    fread(&reclen, 4, 1, f2);
                }
                if (irghost != 0) {
                    // read(2) ((((vi0(j,k,l,m),j=1,jdim),k=1,kdim),l=1,1),m=1,4),
                    //         ((((vj0(k,i,l,m),k=1,kdim),i=1,idim1),l=1,1),m=1,4),
                    //         ((((vk0(j,i,l,m),j=1,jdim),i=1,idim1),l=1,1),m=1,4)
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    for (int m = 1; m <= 4; m++) {
                        for (int l = 1; l <= 1; l++) {
                            for (int k = 1; k <= kdim; k++) {
                                for (int j = 1; j <= jdim; j++) {
                                    float v;
                                    fread(&v, sizeof(float), 1, f2);
                                    vi0(j, k, l, m) = (double)v;
                                }
                            }
                        }
                    }
                    for (int m = 1; m <= 4; m++) {
                        for (int l = 1; l <= 1; l++) {
                            for (int i = 1; i <= idim1; i++) {
                                for (int k = 1; k <= kdim; k++) {
                                    float v;
                                    fread(&v, sizeof(float), 1, f2);
                                    vj0(k, i, l, m) = (double)v;
                                }
                            }
                        }
                    }
                    for (int m = 1; m <= 4; m++) {
                        for (int l = 1; l <= 1; l++) {
                            for (int i = 1; i <= idim1; i++) {
                                for (int j = 1; j <= jdim; j++) {
                                    float v;
                                    fread(&v, sizeof(float), 1, f2);
                                    vk0(j, i, l, m) = (double)v;
                                }
                            }
                        }
                    }
                    fread(&reclen, 4, 1, f2);
                }
            }
        }



        if (iv1 >= 4 || iv2 >= 4 || iv3 >= 4) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                " reading field eqn turb  quantities from restart file, block%6d", nbl);
            if (icgns != 1) {
                // read(2) ((((tursav(j,k,i,m),j=1,jdim1),k=1,kdim1),i=1,idim1),m=1,2)
                {
                    int32_t reclen;
                    fread(&reclen, 4, 1, f2);
                    for (int m = 1; m <= 2; m++) {
                        for (int i = 1; i <= idim1; i++) {
                            for (int k = 1; k <= kdim1; k++) {
                                for (int j = 1; j <= jdim1; j++) {
                                    float v;
                                    fread(&v, sizeof(float), 1, f2);
                                    tursav(j, k, i, m) = (double)v;
                                }
                            }
                        }
                    }
                    fread(&reclen, 4, 1, f2);
                }
                if (iv1 >= 30 || iv2 >= 30 || iv3 >= 30) {
                    // read(2) nummem_read
                    {
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        fread(&nummem_read, sizeof(int), 1, f2);
                        fread(&reclen, 4, 1, f2);
                    }
                    // read(2) ((((tursav(j,k,i,m),j=1,jdim1),k=1,kdim1),i=1,idim1),m=3,nummem_read)
                    {
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        for (int m = 3; m <= nummem_read; m++) {
                            for (int i = 1; i <= idim1; i++) {
                                for (int k = 1; k <= kdim1; k++) {
                                    for (int j = 1; j <= jdim1; j++) {
                                        float v;
                                        fread(&v, sizeof(float), 1, f2);
                                        tursav(j, k, i, m) = (double)v;
                                    }
                                }
                            }
                        }
                        fread(&reclen, 4, 1, f2);
                    }
                }
                if (irghost != 0) {
                    // read(2) ((((ti0(j,k,l,m),j=1,jdim),k=1,kdim),l=1,2),m=1,4),
                    //         ((((tj0(k,i,l,m),k=1,kdim),i=1,idim1),l=1,2),m=1,4),
                    //         ((((tk0(j,i,l,m),j=1,jdim),i=1,idim1),l=1,2),m=1,4)
                    {
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        for (int m = 1; m <= 4; m++) {
                            for (int l = 1; l <= 2; l++) {
                                for (int k = 1; k <= kdim; k++) {
                                    for (int j = 1; j <= jdim; j++) {
                                        float v;
                                        fread(&v, sizeof(float), 1, f2);
                                        ti0(j, k, l, m) = (double)v;
                                    }
                                }
                            }
                        }
                        for (int m = 1; m <= 4; m++) {
                            for (int l = 1; l <= 2; l++) {
                                for (int i = 1; i <= idim1; i++) {
                                    for (int k = 1; k <= kdim; k++) {
                                        float v;
                                        fread(&v, sizeof(float), 1, f2);
                                        tj0(k, i, l, m) = (double)v;
                                    }
                                }
                            }
                        }
                        for (int m = 1; m <= 4; m++) {
                            for (int l = 1; l <= 2; l++) {
                                for (int i = 1; i <= idim1; i++) {
                                    for (int j = 1; j <= jdim; j++) {
                                        float v;
                                        fread(&v, sizeof(float), 1, f2);
                                        tk0(j, i, l, m) = (double)v;
                                    }
                                }
                            }
                        }
                        fread(&reclen, 4, 1, f2);
                    }
                    if (iv1 >= 30 || iv2 >= 30 || iv3 >= 30) {
                        // read(2) ((((ti0(j,k,l,m),j=1,jdim),k=1,kdim),l=3,nummem_read),m=1,4),
                        //         ((((tj0(k,i,l,m),k=1,kdim),i=1,idim1),l=3,nummem_read),m=1,4),
                        //         ((((tk0(j,i,l,m),j=1,jdim),i=1,idim1),l=3,nummem_read),m=1,4)
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        for (int m = 1; m <= 4; m++) {
                            for (int l = 3; l <= nummem_read; l++) {
                                for (int k = 1; k <= kdim; k++) {
                                    for (int j = 1; j <= jdim; j++) {
                                        float v;
                                        fread(&v, sizeof(float), 1, f2);
                                        ti0(j, k, l, m) = (double)v;
                                    }
                                }
                            }
                        }
                        for (int m = 1; m <= 4; m++) {
                            for (int l = 3; l <= nummem_read; l++) {
                                for (int i = 1; i <= idim1; i++) {
                                    for (int k = 1; k <= kdim; k++) {
                                        float v;
                                        fread(&v, sizeof(float), 1, f2);
                                        tj0(k, i, l, m) = (double)v;
                                    }
                                }
                            }
                        }
                        for (int m = 1; m <= 4; m++) {
                            for (int l = 3; l <= nummem_read; l++) {
                                for (int i = 1; i <= idim1; i++) {
                                    for (int j = 1; j <= jdim; j++) {
                                        float v;
                                        fread(&v, sizeof(float), 1, f2);
                                        tk0(j, i, l, m) = (double)v;
                                    }
                                }
                            }
                        }
                        fread(&reclen, 4, 1, f2);
                    }
                }
            }



            // Lines 499-615: turbulence model conversion factors
            // iv old=6/7, new=8/12/14: multiply tursav(:,:,:,1) by 0.09
            if ((iv1 == 6 || iv2 == 6 || iv3 == 6 || iv1 == 7 || iv2 == 7 || iv3 == 7) &&
                (ivisc[0] == 8 || ivisc[1] == 8 || ivisc[2] == 8 ||
                 ivisc[0] == 12 || ivisc[1] == 12 || ivisc[2] == 12 ||
                 ivisc[0] == 14 || ivisc[1] == 14 || ivisc[2] == 14)) {
                for (int i = 1; i <= idim1; i++) {
                    for (int k = 1; k <= kdim1; k++) {
                        for (int j = 1; j <= jdim1; j++) {
                            tursav(j, k, i, 1) = tursav(j, k, i, 1) * 0.09;
                        }
                    }
                }
            }
            // iv old=8/12/14, new=6/7: divide tursav(:,:,:,1) by 0.09
            if ((ivisc[0] == 6 || ivisc[1] == 6 || ivisc[2] == 6 ||
                 ivisc[0] == 7 || ivisc[1] == 7 || ivisc[2] == 7) &&
                (iv1 == 8 || iv2 == 8 || iv3 == 8 ||
                 iv1 == 12 || iv2 == 12 || iv3 == 12 ||
                 iv1 == 14 || iv2 == 14 || iv3 == 14)) {
                for (int i = 1; i <= idim1; i++) {
                    for (int k = 1; k <= kdim1; k++) {
                        for (int j = 1; j <= jdim1; j++) {
                            tursav(j, k, i, 1) = tursav(j, k, i, 1) / 0.09;
                        }
                    }
                }
            }
            // iv old=9/10/11/13/15, new=8/12/14: tursav1 = tursav1/tursav2, limit to >=1
            if ((ivisc[0] == 8 || ivisc[1] == 8 || ivisc[2] == 8 ||
                 ivisc[0] == 12 || ivisc[1] == 12 || ivisc[2] == 12 ||
                 ivisc[0] == 14 || ivisc[1] == 14 || ivisc[2] == 14) &&
                (iv1 == 9 || iv2 == 9 || iv3 == 9 ||
                 iv1 == 10 || iv2 == 10 || iv3 == 10 ||
                 iv1 == 11 || iv2 == 11 || iv3 == 11 ||
                 iv1 == 13 || iv2 == 13 || iv3 == 13 ||
                 iv1 == 15 || iv2 == 15 || iv3 == 15)) {
                for (int i = 1; i <= idim1; i++) {
                    for (int k = 1; k <= kdim1; k++) {
                        for (int j = 1; j <= jdim1; j++) {
                            tursav(j, k, i, 1) = tursav(j, k, i, 1) / tursav(j, k, i, 2);
                            // Limit new omega to handle extremes
                            double a1 = tursav(j, k, i, 1);
                            double a2 = 1.0;
                            tursav(j, k, i, 1) = ccomplex_ns::ccmincr(a1, a2);
                        }
                    }
                }
            }
            // iv old=9/10/11/13/15, new=6/7: tursav1 = tursav1/(0.09*tursav2), limit
            if ((ivisc[0] == 6 || ivisc[1] == 6 || ivisc[2] == 6 ||
                 ivisc[0] == 7 || ivisc[1] == 7 || ivisc[2] == 7) &&
                (iv1 == 9 || iv2 == 9 || iv3 == 9 ||
                 iv1 == 10 || iv2 == 10 || iv3 == 10 ||
                 iv1 == 11 || iv2 == 11 || iv3 == 11 ||
                 iv1 == 13 || iv2 == 13 || iv3 == 13 ||
                 iv1 == 15 || iv2 == 15 || iv3 == 15)) {
                for (int i = 1; i <= idim1; i++) {
                    for (int k = 1; k <= kdim1; k++) {
                        for (int j = 1; j <= jdim1; j++) {
                            tursav(j, k, i, 1) = tursav(j, k, i, 1) / (0.09 * tursav(j, k, i, 2));
                            // Limit new omega to handle extremes
                            double a1 = tursav(j, k, i, 1);
                            double a2 = 1.0;
                            tursav(j, k, i, 1) = ccomplex_ns::ccmincr(a1, a2);
                        }
                    }
                }
            }
            // iv old=8/12/14, new=9/10/11/13/15: tursav1 = tursav1*tursav2, limit
            if ((ivisc[0] == 9 || ivisc[1] == 9 || ivisc[2] == 9 ||
                 ivisc[0] == 10 || ivisc[1] == 10 || ivisc[2] == 10 ||
                 ivisc[0] == 11 || ivisc[1] == 11 || ivisc[2] == 11 ||
                 ivisc[0] == 13 || ivisc[1] == 13 || ivisc[2] == 13 ||
                 ivisc[0] == 15 || ivisc[1] == 15 || ivisc[2] == 15) &&
                (iv1 == 8 || iv2 == 8 || iv3 == 8 ||
                 iv1 == 12 || iv2 == 12 || iv3 == 12 ||
                 iv1 == 14 || iv2 == 14 || iv3 == 14)) {
                for (int i = 1; i <= idim1; i++) {
                    for (int k = 1; k <= kdim1; k++) {
                        for (int j = 1; j <= jdim1; j++) {
                            tursav(j, k, i, 1) = tursav(j, k, i, 1) * tursav(j, k, i, 2);
                            // Limit new epsilon to handle extremes
                            double a1 = tursav(j, k, i, 1);
                            double a2 = 1.0;
                            tursav(j, k, i, 1) = ccomplex_ns::ccmincr(a1, a2);
                        }
                    }
                }
            }
            // iv old=6/7, new=9/10/11/13/15: tursav1 = tursav1*0.09*tursav2, limit
            if ((ivisc[0] == 9 || ivisc[1] == 9 || ivisc[2] == 9 ||
                 ivisc[0] == 10 || ivisc[1] == 10 || ivisc[2] == 10 ||
                 ivisc[0] == 11 || ivisc[1] == 11 || ivisc[2] == 11 ||
                 ivisc[0] == 13 || ivisc[1] == 13 || ivisc[2] == 13 ||
                 ivisc[0] == 15 || ivisc[1] == 15 || ivisc[2] == 15) &&
                (iv1 == 6 || iv2 == 6 || iv3 == 6 ||
                 iv1 == 7 || iv2 == 7 || iv3 == 7)) {
                for (int i = 1; i <= idim1; i++) {
                    for (int k = 1; k <= kdim1; k++) {
                        for (int j = 1; j <= jdim1; j++) {
                            tursav(j, k, i, 1) = tursav(j, k, i, 1) * 0.09 * tursav(j, k, i, 2);
                            // Limit new epsilon to handle extremes
                            double a1 = tursav(j, k, i, 1);
                            double a2 = 1.0;
                            tursav(j, k, i, 1) = ccomplex_ns::ccmincr(a1, a2);
                        }
                    }
                }
            }



            // Lines 616-628: ivisc=72, old=6/7: tursav7=tursav1, tursav1=-tursav2*2/3
            if ((ivisc[0] == 72 || ivisc[1] == 72 || ivisc[2] == 72) &&
                (iv1 == 6 || iv2 == 6 || iv3 == 6 ||
                 iv1 == 7 || iv2 == 7 || iv3 == 7)) {
                for (int i = 1; i <= idim1; i++) {
                    for (int k = 1; k <= kdim1; k++) {
                        for (int j = 1; j <= jdim1; j++) {
                            tursav(j, k, i, 7) = tursav(j, k, i, 1);
                            tursav(j, k, i, 1) = -tursav(j, k, i, 2) * 2.0 / 3.0;
                        }
                    }
                }
            }
            // Lines 629-642: ivisc=72, old=8/12/14: tursav7=tursav1/0.09, tursav1=-tursav2*2/3
            if ((ivisc[0] == 72 || ivisc[1] == 72 || ivisc[2] == 72) &&
                (iv1 == 8 || iv2 == 8 || iv3 == 8 ||
                 iv1 == 12 || iv2 == 12 || iv3 == 12 ||
                 iv1 == 14 || iv2 == 14 || iv3 == 14)) {
                for (int i = 1; i <= idim1; i++) {
                    for (int k = 1; k <= kdim1; k++) {
                        for (int j = 1; j <= jdim1; j++) {
                            tursav(j, k, i, 7) = tursav(j, k, i, 1) / 0.09;
                            tursav(j, k, i, 1) = -tursav(j, k, i, 2) * 2.0 / 3.0;
                        }
                    }
                }
            }

            // Lines 643-646: read smin
            if (icgns != 1) {
                // read(2) (((smin(j,k,i),j=1,jdim1),k=1,kdim1),i=1,idim1)
                int32_t reclen;
                fread(&reclen, 4, 1, f2);
                for (int i = 1; i <= idim1; i++) {
                    for (int k = 1; k <= kdim1; k++) {
                        for (int j = 1; j <= jdim1; j++) {
                            float v;
                            fread(&v, sizeof(float), 1, f2);
                            smin(j, k, i) = (double)v;
                        }
                    }
                }
                fread(&reclen, 4, 1, f2);
            }

            // Lines 647-657: read xjb, xkb, blnum for iv=4 or iv=25
            if (iv1 == 4 || iv2 == 4 || iv3 == 4 ||
                iv1 == 25 || iv2 == 25 || iv3 == 25) {
                if (icgns != 1) {
                    // read(2) (((xjb(j,k,i),...)))
                    {
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        for (int i = 1; i <= idim1; i++) {
                            for (int k = 1; k <= kdim1; k++) {
                                for (int j = 1; j <= jdim1; j++) {
                                    float v;
                                    fread(&v, sizeof(float), 1, f2);
                                    xjb(j, k, i) = (double)v;
                                }
                            }
                        }
                        fread(&reclen, 4, 1, f2);
                    }
                    // read(2) (((xkb(j,k,i),...)))
                    {
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        for (int i = 1; i <= idim1; i++) {
                            for (int k = 1; k <= kdim1; k++) {
                                for (int j = 1; j <= jdim1; j++) {
                                    float v;
                                    fread(&v, sizeof(float), 1, f2);
                                    xkb(j, k, i) = (double)v;
                                }
                            }
                        }
                        fread(&reclen, 4, 1, f2);
                    }
                    // read(2) (((blnum(j,k,i),...)))
                    {
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        for (int i = 1; i <= idim1; i++) {
                            for (int k = 1; k <= kdim1; k++) {
                                for (int j = 1; j <= jdim1; j++) {
                                    float v;
                                    fread(&v, sizeof(float), 1, f2);
                                    blnum(j, k, i) = (double)v;
                                }
                            }
                        }
                        fread(&reclen, 4, 1, f2);
                    }
                }
            }

            // Lines 658-672: read cmuv for iv=8/9/13/14
            if (iv1 == 8 || iv2 == 8 || iv3 == 8 ||
                iv1 == 9 || iv2 == 9 || iv3 == 9 ||
                iv1 == 13 || iv2 == 13 || iv3 == 13 ||
                iv1 == 14 || iv2 == 14 || iv3 == 14) {
                if (icgns != 1) {
                    if (ivmx == 8 || ivmx == 9 || ivmx == 13 || ivmx == 14) {
                        // read(2) (((cmuv(j,k,i),...)))
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        for (int i = 1; i <= idim1; i++) {
                            for (int k = 1; k <= kdim1; k++) {
                                for (int j = 1; j <= jdim1; j++) {
                                    float v;
                                    fread(&v, sizeof(float), 1, f2);
                                    cmuv(j, k, i) = (double)v;
                                }
                            }
                        }
                        fread(&reclen, 4, 1, f2);
                    } else {
                        // read(2) (((dum,...))) - skip
                        int32_t reclen;
                        fread(&reclen, 4, 1, f2);
                        int cnt = jdim1 * kdim1 * idim1;
                        for (int nn = 0; nn < cnt; nn++) { float v; fread(&v, sizeof(float), 1, f2); }
                        fread(&reclen, 4, 1, f2);
                    }
                }
            }
        } // end if (iv1>=4 || iv2>=4 || iv3>=4)
    } // end else (iread != 0)



    // Lines 679-688: irest2 handling
    if (irest2 == 1) {
        ntt = 0;
        ntr = 0;
    } else {
        // ntt = ntr (was defined earlier for message passing)
    }

    // Lines 690-740: iteravg == 2 block
    if (iteravg == 2) {
        if (ipertavg == 0) {
            if (myid == myhost) {
                FILE* f97 = fortran_get_unit(97);
                // read(97,end=1010,err=1010) dum,dum,dum,xnumavg
                {
                    int32_t reclen;
                    if (fread(&reclen, 4, 1, f97) != 1) goto label_1010;
                    float d1, d2, d3, xna;
                    if (fread(&d1, sizeof(float), 1, f97) != 1) goto label_1010;
                    if (fread(&d2, sizeof(float), 1, f97) != 1) goto label_1010;
                    if (fread(&d3, sizeof(float), 1, f97) != 1) goto label_1010;
                    if (fread(&xna, sizeof(float), 1, f97) != 1) goto label_1010;
                    if (fread(&reclen, 4, 1, f97) != 1) goto label_1010;
                    xnumavg = xna;
                }
                // read(97,end=1010,err=1010) ((((qavg(j,k,i,m),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),m=1,5)
                {
                    int32_t reclen;
                    if (fread(&reclen, 4, 1, f97) != 1) goto label_1010;
                    for (int m = 1; m <= 5; m++) {
                        for (int k = 1; k <= kdim - 1; k++) {
                            for (int j = 1; j <= jdim - 1; j++) {
                                for (int i = 1; i <= idim - 1; i++) {
                                    float v;
                                    if (fread(&v, sizeof(float), 1, f97) != 1) goto label_1010;
                                    qavg(j, k, i, m) = (double)v;
                                }
                            }
                        }
                    }
                    if (fread(&reclen, 4, 1, f97) != 1) goto label_1010;
                }
                // need to switch y and z if ialph != 0
                if (ialph != 0) {
                    for (int i = 1; i <= idim - 1; i++) {
                        for (int j = 1; j <= jdim - 1; j++) {
                            for (int k = 1; k <= kdim - 1; k++) {
                                temp = qavg(j, k, i, 4);
                                qavg(j, k, i, 4) = qavg(j, k, i, 3);
                                qavg(j, k, i, 3) = -temp;
                            }
                        }
                    }
                }
            }
        } else {
            if (myid == myhost) {
                FILE* f97 = fortran_get_unit(97);
                // read(97,end=1010,err=1010) dum,dum,dum,xnumavg
                {
                    int32_t reclen;
                    if (fread(&reclen, 4, 1, f97) != 1) goto label_1010;
                    float d1, d2, d3, xna;
                    if (fread(&d1, sizeof(float), 1, f97) != 1) goto label_1010;
                    if (fread(&d2, sizeof(float), 1, f97) != 1) goto label_1010;
                    if (fread(&d3, sizeof(float), 1, f97) != 1) goto label_1010;
                    if (fread(&xna, sizeof(float), 1, f97) != 1) goto label_1010;
                    if (fread(&reclen, 4, 1, f97) != 1) goto label_1010;
                    xnumavg = xna;
                }
                // read(97,end=1010,err=1010) ((((qavg(j,k,i,m),i=1,idim),j=1,jdim),k=1,kdim),m=1,5)
                {
                    int32_t reclen;
                    if (fread(&reclen, 4, 1, f97) != 1) goto label_1010;
                    for (int m = 1; m <= 5; m++) {
                        for (int k = 1; k <= kdim; k++) {
                            for (int j = 1; j <= jdim; j++) {
                                for (int i = 1; i <= idim; i++) {
                                    float v;
                                    if (fread(&v, sizeof(float), 1, f97) != 1) goto label_1010;
                                    qavg(j, k, i, m) = (double)v;
                                }
                            }
                        }
                    }
                    if (fread(&reclen, 4, 1, f97) != 1) goto label_1010;
                }

                FILE* f98 = fortran_get_unit(98);
                // read(98,end=1010,err=1010) dum,dum,dum,xnumavg2
                {
                    int32_t reclen;
                    if (fread(&reclen, 4, 1, f98) != 1) goto label_1010;
                    float d1, d2, d3, xna2;
                    if (fread(&d1, sizeof(float), 1, f98) != 1) goto label_1010;
                    if (fread(&d2, sizeof(float), 1, f98) != 1) goto label_1010;
                    if (fread(&d3, sizeof(float), 1, f98) != 1) goto label_1010;
                    if (fread(&xna2, sizeof(float), 1, f98) != 1) goto label_1010;
                    if (fread(&reclen, 4, 1, f98) != 1) goto label_1010;
                    xnumavg2 = xna2;
                }
                // read(98,end=1010,err=1010) ((((q2avg(j,k,i,m),i=1,idim),j=1,jdim),k=1,kdim),m=1,5)
                {
                    int32_t reclen;
                    if (fread(&reclen, 4, 1, f98) != 1) goto label_1010;
                    for (int m = 1; m <= 5; m++) {
                        for (int k = 1; k <= kdim; k++) {
                            for (int j = 1; j <= jdim; j++) {
                                for (int i = 1; i <= idim; i++) {
                                    float v;
                                    if (fread(&v, sizeof(float), 1, f98) != 1) goto label_1010;
                                    q2avg(j, k, i, m) = (double)v;
                                }
                            }
                        }
                    }
                    if (fread(&reclen, 4, 1, f98) != 1) goto label_1010;
                }

                // need to switch y and z if ialph != 0
                if (ialph != 0) {
                    for (int i = 1; i <= idim; i++) {
                        for (int j = 1; j <= jdim; j++) {
                            for (int k = 1; k <= kdim; k++) {
                                temp = qavg(j, k, i, 4);
                                qavg(j, k, i, 4) = qavg(j, k, i, 3);
                                qavg(j, k, i, 3) = -temp;

                                temp = q2avg(j, k, i, 4);
                                q2avg(j, k, i, 4) = q2avg(j, k, i, 3);
                                q2avg(j, k, i, 3) = temp;
                            }
                        }
                    }
                }
            }
        } // end ipertavg
    } // end iteravg



    // Lines 745-764: clcd information
    if (iskip > 0) {
        if (myid == myhost) {
            if (iclcd == 1 || iclcd == 2) {
                // clcd(1,:,:) = 1.e21
                for (int nn3 = clcd.lbound(3); nn3 <= clcd.ubound(3); nn3++) {
                    for (int nn2 = clcd.lbound(2); nn2 <= clcd.ubound(2); nn2++) {
                        clcd(1, nn2, nn3) = 1.e21;
                    }
                }
                // clcd(2,:,:) = 0.e0
                for (int nn3 = clcd.lbound(3); nn3 <= clcd.ubound(3); nn3++) {
                    for (int nn2 = clcd.lbound(2); nn2 <= clcd.ubound(2); nn2++) {
                        clcd(2, nn2, nn3) = 0.0;
                    }
                }
                i_ntt = 0;
            }
            if (iclcd == 2) {
                FILE* f102 = fortran_get_unit(102);
                // read(102,end=1012,err=1012) i_clcd, iblocks_clcd, i_ntt
                {
                    int32_t reclen;
                    if (fread(&reclen, 4, 1, f102) != 1) goto label_1012;
                    if (fread(&i_clcd_r, sizeof(int), 1, f102) != 1) goto label_1012;
                    if (fread(&iblocks_clcd_r, sizeof(int), 1, f102) != 1) goto label_1012;
                    if (fread(&i_ntt, sizeof(int), 1, f102) != 1) goto label_1012;
                    if (fread(&reclen, 4, 1, f102) != 1) goto label_1012;
                }
                if (n_clcd != i_clcd_r) {
                    FILE* f_out = stdout;
                    fprintf(f_out, "Inconsistent number of ClCd calcs in restart file  %d %d\n",
                        n_clcd, i_clcd_r);
                    int m1 = -1;
                    termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                }
                // read(102,end=1012,err=1012) clcd(1:2,1:i_clcd,1:i_ntt)
                {
                    int32_t reclen;
                    if (fread(&reclen, 4, 1, f102) != 1) goto label_1012;
                    for (int nn3 = 1; nn3 <= i_ntt; nn3++) {
                        for (int nn2 = 1; nn2 <= i_clcd_r; nn2++) {
                            for (int nn1 = 1; nn1 <= 2; nn1++) {
                                float v;
                                if (fread(&v, sizeof(float), 1, f102) != 1) goto label_1012;
                                clcd(nn1, nn2, nn3) = (double)v;
                            }
                        }
                    }
                    if (fread(&reclen, 4, 1, f102) != 1) goto label_1012;
                }
            }
        }
    }

    return;

label_1010:
    {
        FILE* f11 = fortran_get_unit(11);
        fprintf(f11, "\n stopping... flag iteravg=2, but cannot read running-average Q file, nbl=%6d\n", nbl);
        fprintf(f11, " ... it either does not exist or it is an incorrect file\n");
        int m1 = -1;
        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
    }
    return;

label_1012:
    {
        FILE* f11 = fortran_get_unit(11);
        fprintf(f11, "\n stopping... flag iclcd=2, but cannot read clcd.bin file, nbl=%6d\n", nbl);
        fprintf(f11, " ... it either does not exist or it is an incorrect file\n");
        int m1 = -1;
        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
    }
    return;
}

} // namespace rrest_ns
