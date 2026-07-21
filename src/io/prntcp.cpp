// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "prntcp.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>

namespace prntcp_ns {

void prntcp(int& jdim, int& kdim, int& idim,
            FortranArray1DRef<double> cpi, FortranArray1DRef<double> cpj, FortranArray1DRef<double> cpk,
            FortranArray4DRef<double> q, int& nbl, int& maxbl, int& maxseg,
            FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
            FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
            FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim,
            FortranArray1DRef<double> thetay, FortranArray1DRef<int> mblk2nd,
            int& myid, int& myhost, int& mycomm, int& irite)
{

    // COMMON block aliases
    float& gamma     = cmn_fluid.gamma;
    float& xmach     = cmn_info.xmach;
    float& alpha     = cmn_info.alpha;
    float& p0        = cmn_ivals.p0;
    float& time      = cmn_unst.time;
    float& radtodeg  = cmn_conversion.radtodeg;

    // Local variables
    double alphd = (double)radtodeg * ((double)alpha + thetay(nbl));
    double cpc   = 2.0e0 / ((double)gamma * (double)xmach * (double)xmach);

    int idim1 = idim - 1;
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;

    // format 3: (10h     block,10h      face,10h   segment,/,3i10)
    auto fmt3 = [](int b, int f, int s) {
        FILE* fp = fortran_get_unit(20);
        fprintf(fp, "     block      face   segment\n");
        fprintf(fp, "%10d%10d%10d\n", b, f, s);
    };
    // format 4: (10h      time,10h     alpha,10h      jlen,10h      klen)
    auto fmt4 = []() {
        FILE* fp = fortran_get_unit(20);
        fprintf(fp, "      time     alpha      jlen      klen\n");
    };
    // format 5: (10h      time,10h     alpha,10h      klen,10h      ilen)
    auto fmt5 = []() {
        FILE* fp = fortran_get_unit(20);
        fprintf(fp, "      time     alpha      klen      ilen\n");
    };
    // format 6: (10h      time,10h     alpha,10h      jlen,10h      ilen)
    auto fmt6 = []() {
        FILE* fp = fortran_get_unit(20);
        fprintf(fp, "      time     alpha      jlen      ilen\n");
    };
    // format 1: (2f10.5,2i10)
    auto fmt1 = [](float t, float a, int d1, int d2) {
        FILE* fp = fortran_get_unit(20);
        fprintf(fp, "%10.5f%10.5f%10d%10d\n", (double)t, (double)a, d1, d2);
    };

    //
    // pressures on i0 segments
    //
    if (nbci0(nbl) > 0) {
        for (int ns = 1; ns <= nbci0(nbl); ns++) {
            if (std::abs(ibcinfo(nbl,ns,1,1)) == 1005 ||
                std::abs(ibcinfo(nbl,ns,1,1)) == 1006 ||
                std::abs(ibcinfo(nbl,ns,1,1)) == 2004 ||
                std::abs(ibcinfo(nbl,ns,1,1)) == 2014 ||
                std::abs(ibcinfo(nbl,ns,1,1)) == 2024 ||
                std::abs(ibcinfo(nbl,ns,1,1)) == 2034 ||
                std::abs(ibcinfo(nbl,ns,1,1)) == 2016) {
                if (myid == mblk2nd(nbl)) {
                    int jst  = ibcinfo(nbl,ns,2,1);
                    int jfn  = ibcinfo(nbl,ns,3,1) - 1;
                    int kst  = ibcinfo(nbl,ns,4,1);
                    int kfn  = ibcinfo(nbl,ns,5,1) - 1;
                    int jj   = 0;
                    int jlen = jfn - jst + 1;
                    int klen = kfn - kst + 1;
                    for (int j = jst; j <= jfn; j++) {
                        jj = jj + 1;
                        int kk = 0;
                        for (int k = kst; k <= kfn; k++) {
                            kk = kk + 1;
                            int izz = (kk - 1) * jlen + jj;
                            cpi(izz) = (q(j,k,1,5) / (double)p0 - 1.0e0) * cpc;
                        }
                    }
                }
                if (myid == myhost) {
                    int jst  = ibcinfo(nbl,ns,2,1);
                    int jfn  = ibcinfo(nbl,ns,3,1) - 1;
                    int kst  = ibcinfo(nbl,ns,4,1);
                    int kfn  = ibcinfo(nbl,ns,5,1) - 1;
                    int jlen = jfn - jst + 1;
                    int klen = kfn - kst + 1;
                    fmt3(nbl, 1, ns);
                    fmt4();
                    fmt1((float)time, (float)alphd, jlen, klen);
                    for (int j = 1; j <= jlen; j++) {
                        FILE* fp = fortran_get_unit(20);
                        for (int k = 1; k <= klen; k++) {
                            fprintf(fp, "%10.6f", (double)(float)cpi((k-1)*jlen + j));
                            if (k % 10 == 0 || k == klen) fprintf(fp, "\n");
                        }
                    }
                }
            }
        }
    }

    //
    // pressures on idim segments
    //
    if (nbcidim(nbl) > 0) {
        for (int ns = 1; ns <= nbcidim(nbl); ns++) {
            if (std::abs(ibcinfo(nbl,ns,1,2)) == 1005 ||
                std::abs(ibcinfo(nbl,ns,1,2)) == 1006 ||
                std::abs(ibcinfo(nbl,ns,1,2)) == 2004 ||
                std::abs(ibcinfo(nbl,ns,1,2)) == 2014 ||
                std::abs(ibcinfo(nbl,ns,1,2)) == 2024 ||
                std::abs(ibcinfo(nbl,ns,1,2)) == 2034 ||
                std::abs(ibcinfo(nbl,ns,1,2)) == 2016) {
                if (myid == mblk2nd(nbl)) {
                    int jst  = ibcinfo(nbl,ns,2,2);
                    int jfn  = ibcinfo(nbl,ns,3,2) - 1;
                    int kst  = ibcinfo(nbl,ns,4,2);
                    int kfn  = ibcinfo(nbl,ns,5,2) - 1;
                    int jj   = 0;
                    int jlen = jfn - jst + 1;
                    int klen = kfn - kst + 1;
                    for (int j = jst; j <= jfn; j++) {
                        jj = jj + 1;
                        int kk = 0;
                        for (int k = kst; k <= kfn; k++) {
                            kk = kk + 1;
                            int izz = (kk - 1) * jlen + jj;
                            cpi(izz) = (q(j,k,idim1,5) / (double)p0 - 1.0e0) * cpc;
                        }
                    }
                }
                if (myid == myhost) {
                    int jst  = ibcinfo(nbl,ns,2,2);
                    int jfn  = ibcinfo(nbl,ns,3,2) - 1;
                    int kst  = ibcinfo(nbl,ns,4,2);
                    int kfn  = ibcinfo(nbl,ns,5,2) - 1;
                    int jlen = jfn - jst + 1;
                    int klen = kfn - kst + 1;
                    fmt3(nbl, 2, ns);
                    fmt4();
                    fmt1((float)time, (float)alphd, jlen, klen);
                    for (int j = 1; j <= jlen; j++) {
                        FILE* fp = fortran_get_unit(20);
                        for (int k = 1; k <= klen; k++) {
                            fprintf(fp, "%10.6f", (double)(float)cpi((k-1)*jlen + j));
                            if (k % 10 == 0 || k == klen) fprintf(fp, "\n");
                        }
                    }
                }
            }
        }
    }

    //
    // pressures on j0 segments
    //
    if (nbcj0(nbl) > 0) {
        for (int ns = 1; ns <= nbcj0(nbl); ns++) {
            if (std::abs(jbcinfo(nbl,ns,1,1)) == 1005 ||
                std::abs(jbcinfo(nbl,ns,1,1)) == 1006 ||
                std::abs(jbcinfo(nbl,ns,1,1)) == 2004 ||
                std::abs(jbcinfo(nbl,ns,1,1)) == 2014 ||
                std::abs(jbcinfo(nbl,ns,1,1)) == 2024 ||
                std::abs(jbcinfo(nbl,ns,1,1)) == 2034 ||
                std::abs(jbcinfo(nbl,ns,1,1)) == 2016) {
                if (myid == mblk2nd(nbl)) {
                    int ist  = jbcinfo(nbl,ns,2,1);
                    int ifn  = jbcinfo(nbl,ns,3,1) - 1;
                    int kst  = jbcinfo(nbl,ns,4,1);
                    int kfn  = jbcinfo(nbl,ns,5,1) - 1;
                    int ii   = 0;
                    int ilen = ifn - ist + 1;
                    int klen = kfn - kst + 1;
                    for (int i = ist; i <= ifn; i++) {
                        ii = ii + 1;
                        int kk = 0;
                        for (int k = kst; k <= kfn; k++) {
                            kk = kk + 1;
                            int izz = (ii - 1) * klen + kk;
                            cpj(izz) = (q(1,k,i,5) / (double)p0 - 1.0e0) * cpc;
                        }
                    }
                }
                if (myid == myhost) {
                    int ist  = jbcinfo(nbl,ns,2,1);
                    int ifn  = jbcinfo(nbl,ns,3,1) - 1;
                    int kst  = jbcinfo(nbl,ns,4,1);
                    int kfn  = jbcinfo(nbl,ns,5,1) - 1;
                    int ilen = ifn - ist + 1;
                    int klen = kfn - kst + 1;
                    fmt3(nbl, 3, ns);
                    fmt5();
                    fmt1((float)time, (float)alphd, klen, ilen);
                    for (int k = 1; k <= klen; k++) {
                        FILE* fp = fortran_get_unit(20);
                        for (int i = 1; i <= ilen; i++) {
                            fprintf(fp, "%10.6f", (double)(float)cpj((i-1)*klen + k));
                            if (i % 10 == 0 || i == ilen) fprintf(fp, "\n");
                        }
                    }
                }
            }
        }
    }

    //
    // pressures on jdim segments
    //
    if (nbcjdim(nbl) > 0) {
        for (int ns = 1; ns <= nbcjdim(nbl); ns++) {
            if (std::abs(jbcinfo(nbl,ns,1,2)) == 1005 ||
                std::abs(jbcinfo(nbl,ns,1,2)) == 1006 ||
                std::abs(jbcinfo(nbl,ns,1,2)) == 2004 ||
                std::abs(jbcinfo(nbl,ns,1,2)) == 2014 ||
                std::abs(jbcinfo(nbl,ns,1,2)) == 2024 ||
                std::abs(jbcinfo(nbl,ns,1,2)) == 2034 ||
                std::abs(jbcinfo(nbl,ns,1,2)) == 2016) {
                if (myid == mblk2nd(nbl)) {
                    int ist  = jbcinfo(nbl,ns,2,2);
                    int ifn  = jbcinfo(nbl,ns,3,2) - 1;
                    int kst  = jbcinfo(nbl,ns,4,2);
                    int kfn  = jbcinfo(nbl,ns,5,2) - 1;
                    int ii   = 0;
                    int ilen = ifn - ist + 1;
                    int klen = kfn - kst + 1;
                    for (int i = ist; i <= ifn; i++) {
                        ii = ii + 1;
                        int kk = 0;
                        for (int k = kst; k <= kfn; k++) {
                            kk = kk + 1;
                            int izz = (ii - 1) * klen + kk;
                            cpj(izz) = (q(jdim1,k,i,5) / (double)p0 - 1.0e0) * cpc;
                        }
                    }
                }
                if (myid == myhost) {
                    int ist  = jbcinfo(nbl,ns,2,2);
                    int ifn  = jbcinfo(nbl,ns,3,2) - 1;
                    int kst  = jbcinfo(nbl,ns,4,2);
                    int kfn  = jbcinfo(nbl,ns,5,2) - 1;
                    int ilen = ifn - ist + 1;
                    int klen = kfn - kst + 1;
                    fmt3(nbl, 4, ns);
                    fmt5();
                    fmt1((float)time, (float)alphd, klen, ilen);
                    for (int k = 1; k <= klen; k++) {
                        FILE* fp = fortran_get_unit(20);
                        for (int i = 1; i <= ilen; i++) {
                            fprintf(fp, "%10.6f", (double)(float)cpj((i-1)*klen + k));
                            if (i % 10 == 0 || i == ilen) fprintf(fp, "\n");
                        }
                    }
                }
            }
        }
    }

    //
    // pressures on k0 segments
    //
    if (nbck0(nbl) > 0) {
        for (int ns = 1; ns <= nbck0(nbl); ns++) {
            if (std::abs(kbcinfo(nbl,ns,1,1)) == 1005 ||
                std::abs(kbcinfo(nbl,ns,1,1)) == 1006 ||
                std::abs(kbcinfo(nbl,ns,1,1)) == 2004 ||
                std::abs(kbcinfo(nbl,ns,1,1)) == 2014 ||
                std::abs(kbcinfo(nbl,ns,1,1)) == 2024 ||
                std::abs(kbcinfo(nbl,ns,1,1)) == 2034 ||
                std::abs(kbcinfo(nbl,ns,1,1)) == 2016) {
                if (myid == mblk2nd(nbl)) {
                    int ist  = kbcinfo(nbl,ns,2,1);
                    int ifn  = kbcinfo(nbl,ns,3,1) - 1;
                    int jst  = kbcinfo(nbl,ns,4,1);
                    int jfn  = kbcinfo(nbl,ns,5,1) - 1;
                    int ii   = 0;
                    int ilen = ifn - ist + 1;
                    int jlen = jfn - jst + 1;
                    for (int i = ist; i <= ifn; i++) {
                        ii = ii + 1;
                        int jj = 0;
                        for (int j = jst; j <= jfn; j++) {
                            jj = jj + 1;
                            int izz = (ii - 1) * jlen + jj;
                            cpk(izz) = (q(j,1,i,5) / (double)p0 - 1.0e0) * cpc;
                        }
                    }
                }
                if (myid == myhost) {
                    int ist  = kbcinfo(nbl,ns,2,1);
                    int ifn  = kbcinfo(nbl,ns,3,1) - 1;
                    int jst  = kbcinfo(nbl,ns,4,1);
                    int jfn  = kbcinfo(nbl,ns,5,1) - 1;
                    int ilen = ifn - ist + 1;
                    int jlen = jfn - jst + 1;
                    fmt3(nbl, 5, ns);
                    fmt6();
                    fmt1((float)time, (float)alphd, jlen, ilen);
                    for (int i = 1; i <= ilen; i++) {
                        FILE* fp = fortran_get_unit(20);
                        for (int j = 1; j <= jlen; j++) {
                            fprintf(fp, "%10.6f", (double)(float)cpk((i-1)*jlen + j));
                            if (j % 10 == 0 || j == jlen) fprintf(fp, "\n");
                        }
                    }
                }
            }
        }
    }

    //
    // pressures on kdim segments
    //
    if (nbckdim(nbl) > 0) {
        for (int ns = 1; ns <= nbckdim(nbl); ns++) {
            if (std::abs(kbcinfo(nbl,ns,1,2)) == 1005 ||
                std::abs(kbcinfo(nbl,ns,1,2)) == 1006 ||
                std::abs(kbcinfo(nbl,ns,1,2)) == 2004 ||
                std::abs(kbcinfo(nbl,ns,1,2)) == 2014 ||
                std::abs(kbcinfo(nbl,ns,1,2)) == 2024 ||
                std::abs(kbcinfo(nbl,ns,1,2)) == 2034 ||
                std::abs(kbcinfo(nbl,ns,1,2)) == 2016) {
                if (myid == mblk2nd(nbl)) {
                    int ist  = kbcinfo(nbl,ns,2,2);
                    int ifn  = kbcinfo(nbl,ns,3,2) - 1;
                    int jst  = kbcinfo(nbl,ns,4,2);
                    int jfn  = kbcinfo(nbl,ns,5,2) - 1;
                    int ii   = 0;
                    int ilen = ifn - ist + 1;
                    int jlen = jfn - jst + 1;
                    for (int i = ist; i <= ifn; i++) {
                        ii = ii + 1;
                        int jj = 0;
                        for (int j = jst; j <= jfn; j++) {
                            jj = jj + 1;
                            int izz = (ii - 1) * jlen + jj;
                            cpk(izz) = (q(j,kdim1,i,5) / (double)p0 - 1.0e0) * cpc;
                        }
                    }
                }
                if (myid == myhost) {
                    int ist  = kbcinfo(nbl,ns,2,2);
                    int ifn  = kbcinfo(nbl,ns,3,2) - 1;
                    int jst  = kbcinfo(nbl,ns,4,2);
                    int jfn  = kbcinfo(nbl,ns,5,2) - 1;
                    int ilen = ifn - ist + 1;
                    int jlen = jfn - jst + 1;
                    fmt3(nbl, 6, ns);
                    fmt6();
                    fmt1((float)time, (float)alphd, jlen, ilen);
                    for (int i = 1; i <= ilen; i++) {
                        FILE* fp = fortran_get_unit(20);
                        for (int j = 1; j <= jlen; j++) {
                            fprintf(fp, "%10.6f", (double)(float)cpk((i-1)*jlen + j));
                            if (j % 10 == 0 || j == jlen) fprintf(fp, "\n");
                        }
                    }
                }
            }
        }
    }
} // end prntcp

} // namespace prntcp_ns
