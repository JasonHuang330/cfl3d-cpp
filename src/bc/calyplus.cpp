// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "calyplus.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>

namespace calyplus_ns {

void calyplus(int& jdim, int& kdim, int& idim, int& nbl, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, int& iover, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> yplusj, FortranArray3DRef<double> yplusk, FortranArray3DRef<double> yplusi, FortranArray3DRef<double> blankj, FortranArray3DRef<double> blankk, FortranArray3DRef<double> blanki, FortranArray3DRef<double> dnj, FortranArray3DRef<double> dnk, FortranArray3DRef<double> dni, FortranArray3DRef<double> vistj, FortranArray3DRef<double> vistk, FortranArray3DRef<double> visti, int& ifunc, int& iunit, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, int& maxbl, int& maxseg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray3DRef<double> vol)
{
    // COMMON block aliases
    float& gamma   = cmn_fluid.gamma;
    float& cbar    = cmn_fluid2.cbar;
    float& xmach   = cmn_info.xmach;
    float& reue    = cmn_reyue.reue;
    float& tinf    = cmn_reyue.tinf;
    int32_t* ivisc = cmn_reyue.ivisc;  // 1-based: ivisc[0]=ivisc(1), etc.
    int32_t* iwf   = cmn_wallfun.iwf;  // 1-based: iwf[0]=iwf(1), etc.
    int32_t& i2d   = cmn_twod.i2d;

    // /ypinfo/ aliases
    float& ypsumb   = cmn_ypinfo.ypsumb;
    float& ypsumsqb = cmn_ypinfo.ypsumsqb;
    float& ypmaxb   = cmn_ypinfo.ypmaxb;
    float& ypminb   = cmn_ypinfo.ypminb;
    float& dnmaxb   = cmn_ypinfo.dnmaxb;
    float& dnminb   = cmn_ypinfo.dnminb;
    float& ypchk    = cmn_ypinfo.ypchk;
    int32_t& nptsb    = cmn_ypinfo.nptsb;
    int32_t& jypmaxb  = cmn_ypinfo.jypmaxb;
    int32_t& kypmaxb  = cmn_ypinfo.kypmaxb;
    int32_t& iypmaxb  = cmn_ypinfo.iypmaxb;
    int32_t& jypminb  = cmn_ypinfo.jypminb;
    int32_t& kypminb  = cmn_ypinfo.kypminb;
    int32_t& iypminb  = cmn_ypinfo.iypminb;
    int32_t& jdnmaxb  = cmn_ypinfo.jdnmaxb;
    int32_t& kdnmaxb  = cmn_ypinfo.kdnmaxb;
    int32_t& idnmaxb  = cmn_ypinfo.idnmaxb;
    int32_t& jdnminb  = cmn_ypinfo.jdnminb;
    int32_t& kdnminb  = cmn_ypinfo.kdnminb;
    int32_t& idnminb  = cmn_ypinfo.idnminb;
    int32_t& nypchkb  = cmn_ypinfo.nypchkb;

    // Local variables
    int jdim1, kdim1, idim1;
    int m, mm, nbc, kd, kdx, kd1;
    int iflag, npts, nypchk;
    int i, j, k, nseg;
    int i1, i2, j1, j2, k1, k2;
    int id, id1, jd, jd1;
    int iypmax, jypmax, kypmax, iypmin, jypmin, kypmin;
    int idnmax, jdnmax, kdnmax, idnmin, jdnmin, kdnmin;
    int kypmaxin;
    int jj, kk, ii;
    double ypsum, ypsumsq, ypmax, ypmin, ypavg, ypsdev, diff;
    double dnmax, dnmin;
    double q1k1, q2k1, q3k1, q4k1, q5k1, t1k1;
    double q1k2, q2k2, q3k2, q4k2, q5k2, t1k2;
    double q1j1, q2j1, q3j1, q4j1, q5j1, t1j1;
    double q1j2, q2j2, q3j2, q4j2, q5j2, t1j2;
    double q1i1, q2i1, q3i1, q4i1, q5i1, t1i1;
    double q1i2, q2i2, q3i2, q4i2, q5i2, t1i2;
    double vist, avgmut, emuka, pres1, cfy, dn;
    double dx, dy, dz, qa1, qa2;
    int id_loop, jd_loop;  // loop variables that shadow outer scope



    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;

    // Initialize /ypinfo/ block
    ypsumb   = 0.f;
    ypsumsqb = 0.f;
    ypmaxb   = 0.f;
    iypmaxb  = 0;
    jypmaxb  = 0;
    kypmaxb  = 0;
    ypminb   = 1.e9f;
    iypminb  = 0;
    jypminb  = 0;
    kypminb  = 0;
    dnmaxb   = 0.f;
    idnmaxb  = 0;
    jdnmaxb  = 0;
    kdnmaxb  = 0;
    dnminb   = 1.e9f;
    idnminb  = 0;
    jdnminb  = 0;
    kdnminb  = 0;
    nypchkb  = 0;
    nptsb    = 0;

    if (ivisc[2] > 1) {
        // k=1 and/or k=kdim surfaces
        for (m = 1; m <= 2; m++) {

            iflag   = 0;
            ypsum   = 0.;
            ypsumsq = 0.;
            ypmax   = 0.;
            iypmax  = 0;
            jypmax  = 0;
            ypmin   = 1.e9;
            iypmin  = 0;
            jypmin  = 0;
            dnmax   = 0.;
            idnmax  = 0;
            jdnmax  = 0;
            dnmin   = 1.e9;
            idnmin  = 0;
            jdnmin  = 0;
            nypchk  = 0;
            npts    = 0;

            if (m == 1) {
                mm  = 2;
                nbc = nbck0(nbl);
                kd  = 1;
                kdx = 1;
                kd1 = kd + 1;
            } else {
                mm  = 4;
                nbc = nbckdim(nbl);
                kd  = kdim1;
                kdx = kdim;
                kd1 = kd - 1;
            }

            for (i = 1; i <= idim; i++) {
                for (j = 1; j <= jdim; j++) {
                    blankk(j, i, m) = 0.;
                    yplusk(j, i, m) = 0.;
                    dnk(j, i, m)    = 0.;
                    vistk(j, i, m)  = 0.;
                }
            }

            for (nseg = 1; nseg <= nbc; nseg++) {

                if (std::abs(kbcinfo(nbl, nseg, 1, m)) == 2004 ||
                    std::abs(kbcinfo(nbl, nseg, 1, m)) == 2014 ||
                    std::abs(kbcinfo(nbl, nseg, 1, m)) == 2024 ||
                    std::abs(kbcinfo(nbl, nseg, 1, m)) == 2034 ||
                    std::abs(kbcinfo(nbl, nseg, 1, m)) == 2016) {

                    iflag = 1;

                    i1 = kbcinfo(nbl, nseg, 2, m) + 1;
                    i2 = kbcinfo(nbl, nseg, 3, m) - 1;
                    j1 = kbcinfo(nbl, nseg, 4, m) + 1;
                    j2 = kbcinfo(nbl, nseg, 5, m) - 1;
                    if (i1 > i2) i1 = i2;
                    if (j1 > j2) j1 = j2;

                    for (i = i1; i <= i2; i++) {
                        id  = i;
                        id1 = id - 1;
                        if (id1 < 1)    id1 = 1;
                        if (id > idim1) id  = idim1;
                        for (j = j1; j <= j2; j++) {
                            jd  = j;
                            jd1 = jd - 1;
                            if (jd1 < 1)    jd1 = 1;
                            if (jd > jdim1) jd  = jdim1;

                            blankk(j, i, m) = 1.;

                            // wall value - average in J and I
                            q1k1 = .25*(qk0(jd,id,1,mm)  + qk0(jd1,id,1,mm)
                                      + qk0(jd,id1,1,mm) + qk0(jd1,id1,1,mm));
                            q2k1 = .25*(qk0(jd,id,2,mm)  + qk0(jd1,id,2,mm)
                                      + qk0(jd,id1,2,mm) + qk0(jd1,id1,2,mm));
                            q3k1 = .25*(qk0(jd,id,3,mm)  + qk0(jd1,id,3,mm)
                                      + qk0(jd,id1,3,mm) + qk0(jd1,id1,3,mm));
                            q4k1 = .25*(qk0(jd,id,4,mm)  + qk0(jd1,id,4,mm)
                                      + qk0(jd,id1,4,mm) + qk0(jd1,id1,4,mm));
                            q5k1 = .25*(qk0(jd,id,5,mm)  + qk0(jd1,id,5,mm)
                                      + qk0(jd,id1,5,mm) + qk0(jd1,id1,5,mm)) * (double)gamma;
                            t1k1 = q5k1 / q1k1;



                            // first cell center location - average in J and I
                            q1k2 = .25*(q(jd,kd,id,1)  + q(jd1,kd,id,1)
                                      + q(jd,kd,id1,1) + q(jd1,kd,id1,1));
                            q2k2 = .25*(q(jd,kd,id,2)  + q(jd1,kd,id,2)
                                      + q(jd,kd,id1,2) + q(jd1,kd,id1,2));
                            q3k2 = .25*(q(jd,kd,id,3)  + q(jd1,kd,id,3)
                                      + q(jd,kd,id1,3) + q(jd1,kd,id1,3));
                            q4k2 = .25*(q(jd,kd,id,4)  + q(jd1,kd,id,4)
                                      + q(jd,kd,id1,4) + q(jd1,kd,id1,4));
                            q5k2 = .25*(q(jd,kd,id,5)  + q(jd1,kd,id,5)
                                      + q(jd,kd,id1,5) + q(jd1,kd,id1,5)) * (double)gamma;
                            vist = .25*(vist3d(jd,kd,id)  + vist3d(jd1,kd,id)
                                      + vist3d(jd,kd,id1) + vist3d(jd1,kd,id1));
                            t1k2 = q5k2 / q1k2;

                            if (sk(jd ,kdx,id ,4) == 0. ||
                                sk(jd1,kdx,id ,4) == 0. ||
                                sk(jd ,kdx,id1,4) == 0. ||
                                sk(jd1,kdx,id1,4) == 0.) {
                                dx = x(j, kd+1, i) - x(j, kd, i);
                                dy = y(j, kd+1, i) - y(j, kd, i);
                                dz = z(j, kd+1, i) - z(j, kd, i);
                                dn = std::sqrt(dx*dx + dy*dy + dz*dz);
                            } else {
                                dn = (vol(jd ,kd,id )/sk(jd ,kdx,id ,4) +
                                      vol(jd1,kd,id )/sk(jd1,kdx,id ,4) +
                                      vol(jd ,kd,id1)/sk(jd ,kdx,id1,4) +
                                      vol(jd1,kd,id1)/sk(jd1,kdx,id1,4)) / 4.;
                            }
                            if ((float)dn > (float)dnmax) {
                                dnmax  = dn;
                                idnmax = i;
                                jdnmax = j;
                            }
                            if ((float)dn < (float)dnmin) {
                                dnmin  = dn;
                                idnmin = i;
                                jdnmin = j;
                            }
                            dnk(j, i, m)   = dn;
                            vistk(j, i, m) = vist;
                            if ((float)dn == 0.) {
                                dn = 1.;
                            }
                            qa1 = std::sqrt(q2k1*q2k1 + q3k1*q3k1 + q4k1*q4k1);
                            qa2 = std::sqrt(q2k2*q2k2 + q3k2*q3k2 + q4k2*q4k2);
                            // Get turb viscosity at wall (0 unless wall fn used)
                            if (iwf[2] == 1) {
                                avgmut = .25*(vk0(jd,id,1,mm)  + vk0(jd1,id,1,mm)
                                            + vk0(jd,id1,1,mm) + vk0(jd1,id1,1,mm));
                            } else {
                                avgmut = 0.;
                            }
                            emuka = std::pow(t1k1, 1.5) * ((1.0 + (double)cbar/(double)tinf) / (t1k1 + (double)cbar/(double)tinf));
                            pres1 = q5k1;
                            cfy   = 2.*(qa2 - qa1) / dn;
                            cfy   = 2.*(emuka + avgmut)*cfy / ((double)reue * (double)xmach);

                            yplusk(j, i, m) = 0.;
                            if (std::abs((float)cfy) > 0.)
                                yplusk(j, i, m) =
                                    dn*(double)reue*q1k1*std::sqrt(std::abs((float)cfy*0.5f/(float)q1k1)) / emuka;
                            ypsum = ypsum + yplusk(j, i, m);
                            if ((float)yplusk(j,i,m) > (float)ypmax) {
                                ypmax  = yplusk(j, i, m);
                                iypmax = i;
                                jypmax = j;
                            }
                            if ((float)yplusk(j,i,m) > 0. &&
                                (float)yplusk(j,i,m) < (float)ypmin) {
                                ypmin  = yplusk(j, i, m);
                                iypmin = i;
                                jypmin = j;
                            }
                            if ((float)yplusk(j,i,m) > (float)ypchk) {
                                nypchk = nypchk + 1;
                            }

                        } // j loop
                    } // i loop

                    npts = npts + (i2 - i1 + 1)*(j2 - j1 + 1);

                    // 2d output on segment by segment basis
                    if (ifunc > 0 && i2d == 1) {
                        FILE* fu = fortran_get_unit(iunit);
                        fprintf(fu, " zone t = \"  \", i = %d, j = %d, k = %d\n",
                                j2-j1+1, 1, 1);
                        i = 1;
                        for (j = j1; j <= j2; j++) {
                            fprintf(fu, "  %11.5E  %11.5E  %11.5E  %11.5E  %11.5E\n",
                                    (float)x(j,kd,i), (float)z(j,kd,i),
                                    (float)yplusk(j,i,m),
                                    (float)dnk(j,i,m), (float)vistk(j,i,m));
                        }
                    }

                } // if bc type matches
            } // nseg loop



            if (iflag > 0) {

                ypavg = ypsum / npts;
                for (i = 1; i <= idim; i++) {
                    for (j = 1; j <= jdim; j++) {
                        diff    = yplusk(j, i, m) - ypavg;
                        ypsumsq = ypsumsq + diff*diff;
                    }
                }
                ypsdev = std::sqrt(ypsumsq / (npts - 1));

                if (m == 1) {
                    fortran_write_unit(11, "\n  K=1    SURFACE:\n");
                } else {
                    fortran_write_unit(11, "\n  K=KDIM SURFACE:\n");
                }
                fortran_write_unit(11, "     Y+ MAX   JLOC   ILOC       Y+ MIN   JLOC   ILOC\n");
                {
                    FILE* f11 = fortran_get_unit(11);
                    fprintf(f11, " %10.3E    %3d    %3d   %10.3E    %3d    %3d\n",
                            (float)ypmax, jypmax, iypmax,
                            (float)ypmin, jypmin, iypmin);
                }
                fortran_write_unit(11, "     DN MAX   JLOC   ILOC       DN MIN   JLOC   ILOC\n");
                {
                    FILE* f11 = fortran_get_unit(11);
                    fprintf(f11, " %10.3E    %3d    %3d   %10.3E    %3d    %3d\n",
                            (float)dnmax, jdnmax, idnmax,
                            (float)dnmin, jdnmin, idnmin);
                    fprintf(f11, "     Y+ AVG    Y+ STD DEV      NY+ >%2d   NPTS\n",
                            (int)(float)ypchk);
                    fprintf(f11, " %10.3E    %10.3E       %6d %6d\n",
                            (float)ypavg, (float)ypsdev, nypchk, npts);
                }

                if (ifunc > 0 && i2d == 0) {
                    kk = 1;
                    if (m == 2) kk = kdim;
                    {
                        FILE* fu = fortran_get_unit(iunit);
                        for (j = 1; j <= jdim; j++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu, " %E", (float)x(j,kk,i));
                        for (j = 1; j <= jdim; j++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu, " %E", (float)y(j,kk,i));
                        for (j = 1; j <= jdim; j++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu, " %E", (float)z(j,kk,i));
                        for (j = 1; j <= jdim; j++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu, " %d", (int)(float)blankk(j,i,m));
                        fprintf(fu, "\n");
                    }
                    {
                        FILE* fu1 = fortran_get_unit(iunit+1);
                        for (j = 1; j <= jdim; j++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu1, " %E", (float)yplusk(j,i,m));
                        fprintf(fu1, "\n");
                    }
                    {
                        FILE* fu2 = fortran_get_unit(iunit+2);
                        for (j = 1; j <= jdim; j++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu2, " %E", (float)dnk(j,i,m));
                        fprintf(fu2, "\n");
                    }
                    {
                        FILE* fu3 = fortran_get_unit(iunit+3);
                        for (j = 1; j <= jdim; j++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu3, " %E", (float)vistk(j,i,m));
                        fprintf(fu3, "\n");
                    }
                }

                // block totals
                ypsumb   = ypsumb   + (float)ypsum;
                ypsumsqb = ypsumsqb + (float)ypsumsq;
                nptsb    = nptsb    + npts;
                nypchkb  = nypchkb  + nypchk;
                if ((float)ypmax > (float)ypmaxb) {
                    ypmaxb   = (float)ypmax;
                    iypmaxb  = iypmax;
                    jypmaxb  = jypmax;
                    if (m == 1) {
                        kypmaxb = 1;
                    } else {
                        kypmaxb = kdim;
                    }
                }
                if ((float)ypmin < (float)ypminb) {
                    ypminb   = (float)ypmin;
                    iypminb  = iypmin;
                    jypminb  = jypmin;
                    if (m == 1) {
                        kypminb = 1;
                    } else {
                        kypminb = kdim;
                    }
                }
                if ((float)dnmax > (float)dnmaxb) {
                    dnmaxb   = (float)dnmax;
                    idnmaxb  = idnmax;
                    jdnmaxb  = jdnmax;
                    if (m == 1) {
                        kdnmaxb = 1;
                    } else {
                        kdnmaxb = kdim;
                    }
                }
                if ((float)dnmin < (float)dnminb) {
                    dnminb   = (float)dnmin;
                    idnminb  = idnmin;
                    jdnminb  = jdnmin;
                    if (m == 1) {
                        kdnminb = 1;
                    } else {
                        kdnminb = kdim;
                    }
                }
            } // iflag > 0

        } // m loop
    } // ivisc[2] > 1



    if (ivisc[1] > 1) {
        // j=1 and/or j=jdim surfaces
        int jd_surf, jdx_surf, jd1_surf;
        for (m = 1; m <= 2; m++) {

            iflag   = 0;
            ypsum   = 0.;
            ypsumsq = 0.;
            ypmax   = 0.;
            iypmax  = 0;
            kypmax  = 0;
            ypmin   = 1.e9;
            iypmin  = 0;
            kypmaxin = 0;
            dnmax   = 0.;
            idnmax  = 0;
            kdnmax  = 0;
            dnmin   = 1.e9;
            idnmin  = 0;
            kdnmin  = 0;
            nypchk  = 0;
            npts    = 0;

            if (m == 1) {
                mm  = 2;
                nbc = nbcj0(nbl);
                jd_surf  = 1;
                jdx_surf = 1;
                jd1_surf = jd_surf + 1;
            } else {
                mm  = 4;
                nbc = nbcjdim(nbl);
                jd_surf  = jdim1;
                jdx_surf = jdim;
                jd1_surf = jd_surf - 1;
            }

            for (i = 1; i <= idim; i++) {
                for (k = 1; k <= kdim; k++) {
                    blankj(k, i, m) = 0.;
                    yplusj(k, i, m) = 0.;
                    dnj(k, i, m)    = 0.;
                    vistj(k, i, m)  = 0.;
                }
            }

            for (nseg = 1; nseg <= nbc; nseg++) {

                if (std::abs(jbcinfo(nbl, nseg, 1, m)) == 2004 ||
                    std::abs(jbcinfo(nbl, nseg, 1, m)) == 2014 ||
                    std::abs(jbcinfo(nbl, nseg, 1, m)) == 2024 ||
                    std::abs(jbcinfo(nbl, nseg, 1, m)) == 2034 ||
                    std::abs(jbcinfo(nbl, nseg, 1, m)) == 2016) {

                    iflag = 1;

                    i1 = jbcinfo(nbl, nseg, 2, m) + 1;
                    i2 = jbcinfo(nbl, nseg, 3, m) - 1;
                    k1 = jbcinfo(nbl, nseg, 4, m) + 1;
                    k2 = jbcinfo(nbl, nseg, 5, m) - 1;
                    if (i1 > i2) i1 = i2;
                    if (k1 > k2) k1 = k2;

                    for (i = i1; i <= i2; i++) {
                        id  = i;
                        id1 = id - 1;
                        if (id1 < 1)    id1 = 1;
                        if (id > idim1) id  = idim1;
                        for (k = k1; k <= k2; k++) {
                            kd  = k;
                            kd1 = kd - 1;
                            if (kd1 < 1)    kd1 = 1;
                            if (kd > kdim1) kd  = kdim1;

                            blankj(k, i, m) = 1.;

                            // wall value - average in K and I
                            q1j1 = .25*(qj0(kd,id,1,mm)  + qj0(kd1,id,1,mm)
                                      + qj0(kd,id1,1,mm) + qj0(kd1,id1,1,mm));
                            q2j1 = .25*(qj0(kd,id,2,mm)  + qj0(kd1,id,2,mm)
                                      + qj0(kd,id1,2,mm) + qj0(kd1,id1,2,mm));
                            q3j1 = .25*(qj0(kd,id,3,mm)  + qj0(kd1,id,3,mm)
                                      + qj0(kd,id1,3,mm) + qj0(kd1,id1,3,mm));
                            q4j1 = .25*(qj0(kd,id,4,mm)  + qj0(kd1,id,4,mm)
                                      + qj0(kd,id1,4,mm) + qj0(kd1,id1,4,mm));
                            q5j1 = .25*(qj0(kd,id,5,mm)  + qj0(kd1,id,5,mm)
                                      + qj0(kd,id1,5,mm) + qj0(kd1,id1,5,mm)) * (double)gamma;
                            t1j1 = q5j1 / q1j1;

                            // first cell center location - average in K and I
                            q1j2 = .25*(q(jd_surf,kd,id,1)  + q(jd_surf,kd1,id,1)
                                      + q(jd_surf,kd,id1,1) + q(jd_surf,kd1,id1,1));
                            q2j2 = .25*(q(jd_surf,kd,id,2)  + q(jd_surf,kd1,id,2)
                                      + q(jd_surf,kd,id1,2) + q(jd_surf,kd1,id1,2));
                            q3j2 = .25*(q(jd_surf,kd,id,3)  + q(jd_surf,kd1,id,3)
                                      + q(jd_surf,kd,id1,3) + q(jd_surf,kd1,id1,3));
                            q4j2 = .25*(q(jd_surf,kd,id,4)  + q(jd_surf,kd1,id,4)
                                      + q(jd_surf,kd,id1,4) + q(jd_surf,kd1,id1,4));
                            q5j2 = .25*(q(jd_surf,kd,id,5)  + q(jd_surf,kd1,id,5)
                                      + q(jd_surf,kd,id1,5) + q(jd_surf,kd1,id1,5)) * (double)gamma;
                            vist = .25*(vist3d(jd_surf,kd,id)  + vist3d(jd_surf,kd1,id)
                                      + vist3d(jd_surf,kd,id1) + vist3d(jd_surf,kd1,id1));
                            t1j2 = q5j2 / q1j2;



                            if (sj(jdx_surf,kd ,id ,4) == 0. ||
                                sj(jdx_surf,kd1,id ,4) == 0. ||
                                sj(jdx_surf,kd ,id1,4) == 0. ||
                                sj(jdx_surf,kd1,id1,4) == 0.) {
                                dx = x(jd_surf+1, k, i) - x(jd_surf, k, i);
                                dy = y(jd_surf+1, k, i) - y(jd_surf, k, i);
                                dz = z(jd_surf+1, k, i) - z(jd_surf, k, i);
                                dn = std::sqrt(dx*dx + dy*dy + dz*dz);
                            } else {
                                dn = (vol(jd_surf,kd ,id )/sj(jdx_surf,kd ,id ,4) +
                                      vol(jd_surf,kd1,id )/sj(jdx_surf,kd1,id ,4) +
                                      vol(jd_surf,kd ,id1)/sj(jdx_surf,kd ,id1,4) +
                                      vol(jd_surf,kd1,id1)/sj(jdx_surf,kd1,id1,4)) / 4.;
                            }
                            if ((float)dn > (float)dnmax) {
                                dnmax  = dn;
                                idnmax = i;
                                kdnmax = k;
                            }
                            if ((float)dn < (float)dnmin) {
                                dnmin  = dn;
                                idnmin = i;
                                kdnmin = k;
                            }
                            dnj(k, i, m)   = dn;
                            vistj(k, i, m) = vist;
                            if ((float)dn == 0.) {
                                dn = 1.;
                            }
                            qa1 = std::sqrt(q2j1*q2j1 + q3j1*q3j1 + q4j1*q4j1);
                            qa2 = std::sqrt(q2j2*q2j2 + q3j2*q3j2 + q4j2*q4j2);
                            // Get turb viscosity at wall (0 unless wall fn used)
                            if (iwf[1] == 1) {
                                avgmut = .25*(vj0(kd,id,1,mm)  + vj0(kd1,id,1,mm)
                                            + vj0(kd,id1,1,mm) + vj0(kd1,id1,1,mm));
                            } else {
                                avgmut = 0.;
                            }
                            emuka = std::pow(t1j1, 1.5) * ((1.0 + (double)cbar/(double)tinf) / (t1j1 + (double)cbar/(double)tinf));
                            pres1 = q5j1;
                            cfy   = 2.*(qa2 - qa1) / dn;
                            cfy   = 2.*(emuka + avgmut)*cfy / ((double)reue * (double)xmach);

                            yplusj(k, i, m) = 0.;
                            if (std::abs((float)cfy) > 0.)
                                yplusj(k, i, m) =
                                    dn*(double)reue*q1j1*std::sqrt(std::abs((float)cfy*0.5f/(float)q1j1)) / emuka;
                            ypsum = ypsum + yplusj(k, i, m);
                            if ((float)yplusj(k,i,m) > (float)ypmax) {
                                ypmax  = yplusj(k, i, m);
                                iypmax = i;
                                kypmax = k;
                            }
                            if ((float)yplusj(k,i,m) > 0. &&
                                (float)yplusj(k,i,m) < (float)ypmin) {
                                ypmin  = yplusj(k, i, m);
                                iypmin = i;
                                kypmin = k;
                            }
                            if ((float)yplusj(k,i,m) > (float)ypchk) {
                                nypchk = nypchk + 1;
                            }

                        } // k loop
                    } // i loop

                    npts = npts + (i2 - i1 + 1)*(k2 - k1 + 1);

                    // 2d output on segment by segment basis
                    if (ifunc > 0 && i2d == 1) {
                        FILE* fu = fortran_get_unit(iunit);
                        fprintf(fu, " zone t = \"  \", i = %d, j = %d, k = %d\n",
                                k2-k1+1, 1, 1);
                        i = 1;
                        fprintf(fu, " npts = %d\n", k2-k1+1);
                        fprintf(fu, "       X            Y            Y+           DN      TURB. VISC.\n");
                        for (k = k1; k <= k2; k++) {
                            fprintf(fu, "  %11.5E  %11.5E  %11.5E  %11.5E  %11.5E\n",
                                    (float)x(jd_surf,k,i), (float)z(jd_surf,k,i),
                                    (float)yplusj(k,i,m),
                                    (float)dnj(k,i,m), (float)vistj(k,i,m));
                        }
                    }

                } // if bc type matches
            } // nseg loop



            if (iflag > 0) {

                ypavg = ypsum / npts;
                for (i = i1; i <= i2; i++) {
                    for (k = k1; k <= k2; k++) {
                        diff    = yplusj(k, i, m) - ypavg;
                        ypsumsq = ypsumsq + diff*diff;
                    }
                }
                ypsdev = std::sqrt(ypsumsq / (npts - 1));

                if (m == 1) {
                    fortran_write_unit(11, "\n  J=1    SURFACE:\n");
                } else {
                    fortran_write_unit(11, "\n  J=JDIM SURFACE:\n");
                }
                fortran_write_unit(11, "     Y+ MAX   KLOC   ILOC       Y+ MIN   KLOC   ILOC\n");
                {
                    FILE* f11 = fortran_get_unit(11);
                    fprintf(f11, " %10.3E    %3d    %3d   %10.3E    %3d    %3d\n",
                            (float)ypmax, kypmax, iypmax,
                            (float)ypmin, kypmin, iypmin);
                }
                fortran_write_unit(11, "     DN MAX   KLOC   ILOC       DN MIN   KLOC   ILOC\n");
                {
                    FILE* f11 = fortran_get_unit(11);
                    fprintf(f11, " %10.3E    %3d    %3d   %10.3E    %3d    %3d\n",
                            (float)dnmax, kdnmax, idnmax,
                            (float)dnmin, kdnmin, idnmin);
                    fprintf(f11, "     Y+ AVG    Y+ STD DEV      NY+ >%2d   NPTS\n",
                            (int)(float)ypchk);
                    fprintf(f11, " %10.3E    %10.3E       %6d %6d\n",
                            (float)ypavg, (float)ypsdev, nypchk, npts);
                }

                if (ifunc > 0 && i2d == 0) {
                    jj = 1;
                    if (m == 2) jj = jdim;
                    {
                        FILE* fu = fortran_get_unit(iunit);
                        for (k = 1; k <= kdim; k++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu, " %E", (float)x(jj,k,i));
                        for (k = 1; k <= kdim; k++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu, " %E", (float)y(jj,k,i));
                        for (k = 1; k <= kdim; k++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu, " %E", (float)z(jj,k,i));
                        for (k = 1; k <= kdim; k++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu, " %d", (int)(float)blankj(k,i,m));
                        fprintf(fu, "\n");
                    }
                    {
                        FILE* fu1 = fortran_get_unit(iunit+1);
                        for (k = 1; k <= kdim; k++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu1, " %E", (float)yplusj(k,i,m));
                        fprintf(fu1, "\n");
                    }
                    {
                        FILE* fu2 = fortran_get_unit(iunit+2);
                        for (k = 1; k <= kdim; k++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu2, " %E", (float)dnj(k,i,m));
                        fprintf(fu2, "\n");
                    }
                    {
                        FILE* fu3 = fortran_get_unit(iunit+3);
                        for (k = 1; k <= kdim; k++)
                            for (i = 1; i <= idim; i++)
                                fprintf(fu3, " %E", (float)vistj(k,i,m));
                        fprintf(fu3, "\n");
                    }
                }

                // block totals
                ypsumb   = ypsumb   + (float)ypsum;
                ypsumsqb = ypsumsqb + (float)ypsumsq;
                nptsb    = nptsb    + npts;
                nypchkb  = nypchkb  + nypchk;
                if ((float)ypmax > (float)ypmaxb) {
                    ypmaxb   = (float)ypmax;
                    iypmaxb  = iypmax;
                    kypmaxb  = kypmax;
                    if (m == 1) {
                        jypmaxb = 1;
                    } else {
                        jypmaxb = jdim;
                    }
                }
                if ((float)ypmin < (float)ypminb) {
                    ypminb   = (float)ypmin;
                    iypminb  = iypmin;
                    kypminb  = kypmin;
                    if (m == 1) {
                        jypminb = 1;
                    } else {
                        jypminb = jdim;
                    }
                }
                if ((float)dnmax > (float)dnmaxb) {
                    dnmaxb   = (float)dnmax;
                    idnmaxb  = idnmax;
                    kdnmaxb  = kdnmax;
                    if (m == 1) {
                        jdnmaxb = 1;
                    } else {
                        jdnmaxb = jdim;
                    }
                }
                if ((float)dnmin < (float)dnminb) {
                    dnminb   = (float)dnmin;
                    idnminb  = idnmin;
                    kdnminb  = kdnmin;
                    if (m == 1) {
                        jdnminb = 1;
                    } else {
                        jdnminb = jdim;
                    }
                }
            } // iflag > 0

        } // m loop
    } // ivisc[1] > 1



    if (idim == 2) goto label_999;

    if (ivisc[0] > 1) {
        // i=1 and/or i=idim surfaces
        int id_surf, idx_surf, id1_surf;
        for (m = 1; m <= 2; m++) {

            iflag   = 0;
            ypsum   = 0.;
            ypsumsq = 0.;
            ypmax   = 0.;
            jypmax  = 0;
            kypmax  = 0;
            ypmin   = 1.e9;
            jypmin  = 0;
            kypmin  = 0;
            dnmax   = 0.;
            jdnmax  = 0;
            kdnmax  = 0;
            dnmin   = 1.e9;
            jdnmin  = 0;
            kdnmin  = 0;
            nypchk  = 0;
            npts    = 0;

            if (m == 1) {
                mm  = 2;
                nbc = nbci0(nbl);
                id_surf  = 1;
                idx_surf = 1;
                id1_surf = id_surf + 1;
            } else {
                mm  = 4;
                nbc = nbcidim(nbl);
                id_surf  = idim1;
                idx_surf = idim;
                id1_surf = id_surf - 1;
            }

            for (k = 1; k <= kdim; k++) {
                for (j = 1; j <= jdim; j++) {
                    blanki(j, k, m) = 0.;
                    yplusi(j, k, m) = 0.;
                    dni(j, k, m)    = 0.;
                    visti(j, k, m)  = 0.;
                }
            }

            for (nseg = 1; nseg <= nbc; nseg++) {

                if (std::abs(ibcinfo(nbl, nseg, 1, m)) == 2004 ||
                    std::abs(ibcinfo(nbl, nseg, 1, m)) == 2014 ||
                    std::abs(ibcinfo(nbl, nseg, 1, m)) == 2024 ||
                    std::abs(ibcinfo(nbl, nseg, 1, m)) == 2034 ||
                    std::abs(ibcinfo(nbl, nseg, 1, m)) == 2016) {

                    iflag = 1;

                    j1 = ibcinfo(nbl, nseg, 2, m) + 1;
                    j2 = ibcinfo(nbl, nseg, 3, m) - 1;
                    k1 = ibcinfo(nbl, nseg, 4, m) + 1;
                    k2 = ibcinfo(nbl, nseg, 5, m) - 1;
                    if (j1 > j2) j1 = j2;
                    if (k1 > k2) k1 = k2;

                    for (j = j1; j <= j2; j++) {
                        jd  = j;
                        jd1 = jd - 1;
                        if (jd1 < 1)    jd1 = 1;
                        if (jd > jdim1) jd  = jdim1;
                        for (k = k1; k <= k2; k++) {
                            kd  = k;
                            kd1 = kd - 1;
                            if (kd1 < 1)    kd1 = 1;
                            if (kd > kdim1) kd  = kdim1;

                            blanki(j, k, m) = 1.;

                            // wall value - average in J and K
                            q1i1 = .25*(qi0(jd,kd,1,mm)  + qi0(jd1,kd,1,mm)
                                      + qi0(jd,kd1,1,mm) + qi0(jd1,kd1,1,mm));
                            q2i1 = .25*(qi0(jd,kd,2,mm)  + qi0(jd1,kd,2,mm)
                                      + qi0(jd,kd1,2,mm) + qi0(jd1,kd1,2,mm));
                            q3i1 = .25*(qi0(jd,kd,3,mm)  + qi0(jd1,kd,3,mm)
                                      + qi0(jd,kd1,3,mm) + qi0(jd1,kd1,3,mm));
                            q4i1 = .25*(qi0(jd,kd,4,mm)  + qi0(jd1,kd,4,mm)
                                      + qi0(jd,kd1,4,mm) + qi0(jd1,kd1,4,mm));
                            q5i1 = .25*(qi0(jd,kd,5,mm)  + qi0(jd1,kd,5,mm)
                                      + qi0(jd,kd1,5,mm) + qi0(jd1,kd1,5,mm)) * (double)gamma;
                            t1i1 = q5i1 / q1i1;

                            // first cell center location - average in J and K
                            q1i2 = .25*(q(jd,kd,id_surf,1)  + q(jd,kd1,id_surf,1)
                                      + q(jd1,kd,id_surf,1) + q(jd1,kd1,id_surf,1));
                            q2i2 = .25*(q(jd,kd,id_surf,2)  + q(jd,kd1,id_surf,2)
                                      + q(jd1,kd,id_surf,2) + q(jd1,kd1,id_surf,2));
                            q3i2 = .25*(q(jd,kd,id_surf,3)  + q(jd,kd1,id_surf,3)
                                      + q(jd1,kd,id_surf,3) + q(jd1,kd1,id_surf,3));
                            q4i2 = .25*(q(jd,kd,id_surf,4)  + q(jd,kd1,id_surf,4)
                                      + q(jd1,kd,id_surf,4) + q(jd1,kd1,id_surf,4));
                            q5i2 = .25*(q(jd,kd,id_surf,5)  + q(jd,kd1,id_surf,5)
                                      + q(jd1,kd,id_surf,5) + q(jd1,kd1,id_surf,5)) * (double)gamma;
                            vist = .25*(vist3d(jd,kd,id_surf)  + vist3d(jd,kd1,id_surf)
                                      + vist3d(jd1,kd,id_surf) + vist3d(jd1,kd1,id_surf));
                            t1i2 = q5i2 / q1i2;



                            if (si(jd ,kd ,idx_surf,4) == 0. ||
                                si(jd ,kd1,idx_surf,4) == 0. ||
                                si(jd1,kd ,idx_surf,4) == 0. ||
                                si(jd1,kd1,idx_surf,4) == 0.) {
                                dx = x(j, k, id_surf+1) - x(j, k, id_surf);
                                dy = y(j, k, id_surf+1) - y(j, k, id_surf);
                                dz = z(j, k, id_surf+1) - z(j, k, id_surf);
                                dn = std::sqrt(dx*dx + dy*dy + dz*dz);
                            } else {
                                dn = (vol(jd ,kd ,id_surf)/si(jd ,kd ,idx_surf,4) +
                                      vol(jd ,kd1,id_surf)/si(jd ,kd1,idx_surf,4) +
                                      vol(jd1,kd ,id_surf)/si(jd1,kd ,idx_surf,4) +
                                      vol(jd1,kd1,id_surf)/si(jd1,kd1,idx_surf,4)) / 4.;
                            }
                            if ((float)dn > (float)dnmax) {
                                dnmax  = dn;
                                jdnmax = j;
                                kdnmax = k;
                            }
                            if ((float)dn < (float)dnmin) {
                                dnmin  = dn;
                                jdnmin = j;
                                kdnmin = k;
                            }
                            dni(j, k, m)   = dn;
                            visti(j, k, m) = vist;
                            if ((float)dn == 0.) {
                                dn = 1.;
                            }
                            qa1 = std::sqrt(q2i1*q2i1 + q3i1*q3i1 + q4i1*q4i1);
                            qa2 = std::sqrt(q2i2*q2i2 + q3i2*q3i2 + q4i2*q4i2);
                            // Get turb viscosity at wall (0 unless wall fn used)
                            if (iwf[0] == 1) {
                                avgmut = .25*(vi0(jd,kd,1,mm)  + vi0(jd1,kd,1,mm)
                                            + vi0(jd,kd1,1,mm) + vi0(jd1,kd1,1,mm));
                            } else {
                                avgmut = 0.;
                            }
                            emuka = std::pow(t1i1, 1.5) * ((1.0 + (double)cbar/(double)tinf) / (t1i1 + (double)cbar/(double)tinf));
                            pres1 = q5i1;
                            cfy   = 2.*(qa2 - qa1) / dn;
                            cfy   = 2.*(emuka + avgmut)*cfy / ((double)reue * (double)xmach);

                            yplusi(j, k, m) = 0.;
                            if (std::abs((float)cfy) > 0.)
                                yplusi(j, k, m) =
                                    dn*(double)reue*q1i1*std::sqrt(std::abs((float)cfy*0.5f/(float)q1i1)) / emuka;
                            ypsum = ypsum + yplusi(j, k, m);
                            if ((float)yplusi(j,k,m) > (float)ypmax) {
                                ypmax  = yplusi(j, k, m);
                                jypmax = j;
                                kypmax = k;
                            }
                            if ((float)yplusi(j,k,m) > 0. &&
                                (float)yplusi(j,k,m) < (float)ypmin) {
                                ypmin  = yplusi(j, k, m);
                                jypmin = j;
                                kypmin = k;
                            }
                            if ((float)yplusi(j,k,m) > (float)ypchk) {
                                nypchk = nypchk + 1;
                            }

                        } // k loop
                    } // j loop

                    npts = npts + (k2 - k1 + 1)*(j2 - j1 + 1);

                } // if bc type matches
            } // nseg loop

            if (iflag > 0) {

                ypavg = ypsum / npts;
                for (j = j1; j <= j2; j++) {
                    for (k = k1; k <= k2; k++) {
                        diff    = yplusi(j, k, m) - ypavg;
                        ypsumsq = ypsumsq + diff*diff;
                    }
                }
                ypsdev = std::sqrt(ypsumsq / (npts - 1));

                if (m == 1) {
                    fortran_write_unit(11, "\n  I=1    SURFACE:\n");
                } else {
                    fortran_write_unit(11, "\n  I=IDIM SURFACE:\n");
                }
                fortran_write_unit(11, "     Y+ MAX   JLOC   KLOC       Y+ MIN   JLOC   KLOC\n");
                {
                    FILE* f11 = fortran_get_unit(11);
                    fprintf(f11, " %10.3E    %3d    %3d   %10.3E    %3d    %3d\n",
                            (float)ypmax, jypmax, kypmax,
                            (float)ypmin, jypmin, kypmin);
                }
                fortran_write_unit(11, "     DN MAX   JLOC   KLOC       DN MIN   JLOC   KLOC\n");
                {
                    FILE* f11 = fortran_get_unit(11);
                    fprintf(f11, " %10.3E    %3d    %3d   %10.3E    %3d    %3d\n",
                            (float)dnmax, jdnmax, kdnmax,
                            (float)dnmin, jdnmin, kdnmin);
                    fprintf(f11, "     Y+ AVG    Y+ STD DEV      NY+ >%2d   NPTS\n",
                            (int)(float)ypchk);
                    fprintf(f11, " %10.3E    %10.3E       %6d %6d\n",
                            (float)ypavg, (float)ypsdev, nypchk, npts);
                }



                if (ifunc > 0) {
                    ii = 1;
                    if (m == 2) ii = idim;
                    {
                        FILE* fu = fortran_get_unit(iunit);
                        for (k = 1; k <= kdim; k++)
                            for (j = 1; j <= jdim; j++)
                                fprintf(fu, " %E", (float)x(j,k,ii));
                        for (k = 1; k <= kdim; k++)
                            for (j = 1; j <= jdim; j++)
                                fprintf(fu, " %E", (float)y(j,k,ii));
                        for (k = 1; k <= kdim; k++)
                            for (j = 1; j <= jdim; j++)
                                fprintf(fu, " %E", (float)z(j,k,ii));
                        for (k = 1; k <= kdim; k++)
                            for (j = 1; j <= jdim; j++)
                                fprintf(fu, " %d", (int)(float)blanki(j,k,m));
                        fprintf(fu, "\n");
                    }
                    {
                        FILE* fu1 = fortran_get_unit(iunit+1);
                        for (k = 1; k <= kdim; k++)
                            for (j = 1; j <= jdim; j++)
                                fprintf(fu1, " %E", (float)yplusi(j,k,m));
                        fprintf(fu1, "\n");
                    }
                    {
                        FILE* fu2 = fortran_get_unit(iunit+2);
                        for (k = 1; k <= kdim; k++)
                            for (j = 1; j <= jdim; j++)
                                fprintf(fu2, " %E", (float)dni(j,k,m));
                        fprintf(fu2, "\n");
                    }
                    {
                        FILE* fu3 = fortran_get_unit(iunit+3);
                        for (k = 1; k <= kdim; k++)
                            for (j = 1; j <= jdim; j++)
                                fprintf(fu3, " %E", (float)visti(j,k,m));
                        fprintf(fu3, "\n");
                    }
                }

                // block totals
                ypsumb   = ypsumb   + (float)ypsum;
                ypsumsqb = ypsumsqb + (float)ypsumsq;
                nptsb    = nptsb    + npts;
                nypchkb  = nypchkb  + nypchk;
                if ((float)ypmax > (float)ypmaxb) {
                    ypmaxb   = (float)ypmax;
                    kypmaxb  = kypmax;
                    jypmaxb  = jypmax;
                    if (m == 1) {
                        iypmaxb = 1;
                    } else {
                        iypmaxb = idim;
                    }
                }
                if ((float)ypmin < (float)ypminb) {
                    ypminb   = (float)ypmin;
                    kypminb  = kypmin;
                    jypminb  = jypmin;
                    if (m == 1) {
                        iypminb = 1;
                    } else {
                        iypminb = idim;
                    }
                }
                if ((float)dnmax > (float)dnmaxb) {
                    dnmaxb   = (float)dnmax;
                    kdnmaxb  = kdnmax;
                    jdnmaxb  = jdnmax;
                    if (m == 1) {
                        idnmaxb = 1;
                    } else {
                        idnmaxb = idim;
                    }
                }
                if ((float)dnmin < (float)dnminb) {
                    dnminb   = (float)dnmin;
                    kdnminb  = kdnmin;
                    jdnminb  = jdnmin;
                    if (m == 1) {
                        idnminb = 1;
                    } else {
                        idnminb = idim;
                    }
                }
            } // iflag > 0

        } // m loop
    } // ivisc[0] > 1

label_999:
    return;
}

} // namespace calyplus_ns
