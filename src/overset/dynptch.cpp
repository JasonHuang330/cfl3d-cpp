// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dynptch.h"
#include <cstdio>
#include "runtime/fortran_io.h"
#include "global2.h"
#include "lead.h"
#include "patcher.h"
#include "writ_buf.h"
#include "pre_patch.h"
#include "termn8.h"
#include "my_flush.h"
#include <cstdio>
#include <algorithm>
#include <cstring>

namespace dynptch_ns {

void dynptch(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> work, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& ncall, int& maxgr, int& maxbl, int& msub1, int& intmx, int& mxxe, int& mptch, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, int& nblock, int& ngrid, FortranArray1DRef<int> levelg, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray2DRef<double> windex, int& ninter, FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt, FortranArray2DRef<double> windx, int& nintr, FortranArray2DRef<int> iindx, FortranArray1DRef<int> mblkpt, FortranArray1DRef<int> llimit, FortranArray1DRef<int> iitmax, FortranArray1DRef<int> mmcxie, FortranArray1DRef<int> mmceta, FortranArray1DRef<int> ncheck, FortranArray1DRef<int> iifit, FortranArray1DRef<int> iic0, FortranArray1DRef<int> iiorph, FortranArray1DRef<int> iitoss, FortranArray1DRef<int> ifiner, FortranArray2DRef<double> factjlo, FortranArray2DRef<double> factjhi, FortranArray2DRef<double> factklo, FortranArray2DRef<double> factkhi, FortranArray2DRef<double> dx, FortranArray2DRef<double> dy, FortranArray2DRef<double> dz, FortranArray2DRef<double> dthetx, FortranArray2DRef<double> dthety, FortranArray2DRef<double> dthetz, FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz, FortranArray2DRef<int> isav_dpat, FortranArray3DRef<int> isav_dpat_b, int& intmax, int& maxxe, int& nsub1, FortranArray2DRef<int> lw_temp, FortranArray1DRef<int> ireq_ar, FortranArray1DRef<int> lout, FortranArray1DRef<int> ifrom, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, FortranArray1DRef<int> iiint1, FortranArray1DRef<int> iiint2, FortranArray1DRef<int> nblk1, FortranArray1DRef<int> nblk2, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, FortranArray1DRef<int> jte, FortranArray1DRef<int> kte, FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm, FortranArray3DRef<double> xte, FortranArray3DRef<double> yte, FortranArray3DRef<double> zte, FortranArray3DRef<double> xmi, FortranArray3DRef<double> ymi, FortranArray3DRef<double> zmi, FortranArray3DRef<double> xmie, FortranArray3DRef<double> ymie, FortranArray3DRef<double> zmie, FortranArray3DRef<double> sxie, FortranArray3DRef<double> seta, FortranArray2DRef<double> sxie2, FortranArray2DRef<double> seta2, FortranArray2DRef<double> xie2s, FortranArray2DRef<double> eta2s, FortranArray1DRef<double> temp, FortranArray2DRef<double> x2, FortranArray2DRef<double> y2, FortranArray2DRef<double> z2, FortranArray2DRef<double> x1, FortranArray2DRef<double> y1, FortranArray2DRef<double> z1, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray1DRef<int> igridg, FortranArray1DRef<int> iemg)
{
    // COMMON block references
    int& ninter0  = cmn_savnint.ninter0;
    int& isklton  = cmn_sklton.isklton;
    char* grdmov  = cmn_filenam.grdmov;
    int& maxdcnt  = cmn_is_dpatch.maxdcnt;
    int& itrace   = cmn_tracer.itrace;

    // Local variables
    int icallgl, iunit, ioflag, imode;
    int ioutpt, it_thro;
    int ierrflg, icount_dpat;
    int len0, nfb0, nfb, nfb1, n, ll, mm, nbl, nnn, lst, len, l;

    // read in the dynamic patch data at the bottom of the input file
    if (myid == myhost) {
        if (ncall == 1) {
            icallgl = 1;
            iunit   = 25;
            ioflag  = 2;
            imode   = 1;
            global2_ns::global2(maxbl, maxgr, msub1, nintr, intmx, ngrid,
                                idimg, jdimg, kdimg, levelg, ncgg, nblg,
                                iindx, llimit, iitmax, mmcxie, mmceta,
                                ncheck, iifit, iic0, iiorph, iitoss, ifiner,
                                dx, dy, dz, dthetx, dthety, dthetz,
                                myid, mptch, mxxe, icallgl, iunit,
                                nou, bou, ibufdim, nbuf, ifrom, xif1, etf1,
                                xif2, etf2, igridg, iemg, nblock, ioflag, imode);
        }
    }

    // set tracing flag for debugging problem patch interfaces; note:
    // 1) itrace should be set < 0 unless debugging is required
    // 2) trace output is found in fort.7
    // 3) currently trace output is NOT supported in parallel mode
    //
    //   itrace < 0, do not write search history for current "to" cell
    //   itrace = 0, overwrite history from previous "to" cell with current
    //   itrace = 1, retain the search history for ALL cells (may get huge file)
    //
    itrace = -1;

    if (nintr == 0) return;

    // ioutpt controls output for patch diagnostics...output patch
    // diagnostics only on first call to patcher
    ioutpt = 0;
    if (ncall == 1) ioutpt = 1;

    // save number of regular (not dynamic) patch interpolations
    if (ncall == 1) {
        ninter0 = ninter;
    }

    if (myid == myhost) {
        if (isklton == 1) {
            fortran_write_unit(11, "\n calculating dynamic grid interpolation coefficients\n");
        }
    }

    // zero out dynamic patch interpolation coefficients (first call to
    // patcher only - for subsequent calls, the previous values are used
    // as starting points for the search routine
    if (ncall == 1) {
        for (ll = 1; ll <= 2; ll++) {
            for (mm = 1; mm <= mxxe; mm++) {
                windx(mm, ll) = 0.;
            }
        }
    }

    // set up pre_patch data for dynamic grids
    ierrflg = -1;

    if (ncall == 1) {
        icount_dpat = 0;
        for (nbl = 1; nbl <= nblock; nbl++) {
            pre_patch_ns::pre_patch(nbl, lw, icount_dpat, nintr,
                                    iindx, intmx, msub1, isav_dpat,
                                    isav_dpat_b, jjmax1, kkmax1,
                                    iiint1, iiint2, maxbl, jdimg, kdimg, idimg,
                                    ierrflg);
        }

        if (myid == myhost) {
            maxdcnt = icount_dpat;
        }
    }

    // cycle through blocks, check for and establishing connection
    // information for those blocks with patching
    if (myid == myhost) {
        if (ioutpt > 0) {
            fortran_write_unit(25, "\n\n BEGINNING GENERALIZED-COORDINATE INTERPOLATION\n");
        }
    }

    it_thro = 0;
    for (nbl = 1; nbl <= nblock; nbl++) {
        lead_ns::lead(nbl, lw, lw2, maxbl);
        it_thro = it_thro + 1;
        patcher_ns::patcher(nbl, lw, work, mgwk, wk, nwork, ncall, ioutpt,
                            it_thro, maxbl, msub1, intmx, mxxe, mptch,
                            jdimg, kdimg, idimg, windx, nintr, iindx,
                            llimit, iitmax, mmcxie, mmceta, ncheck, iifit,
                            mblkpt, iic0, iiorph, iitoss, ifiner,
                            factjlo, factjhi, factklo, factkhi,
                            dx, dy, dz, dthetx, dthety, dthetz,
                            isav_dpat, isav_dpat_b,
                            xte, yte, zte, xmi, ymi, zmi, xmie, ymie, zmie,
                            jjmax1, kkmax1, jimage, kimage,
                            xorig, yorig, zorig,
                            jte, kte, sxie, seta, sxie2, seta2, xie2s, eta2s,
                            temp, x2, y2, z2, nblk1, nblk2, jmm, kmm,
                            x1, y1, z1, lout, xif1, xif2, etf1, etf2,
                            ireq_ar, myid, myhost, mycomm, mblk2nd,
                            nou, bou, nbuf, ibufdim);
        if (ioutpt > 0) {
            int iunit25 = 25;
            int iunit9  = 9;
            writ_buf_ns::writ_buf(nbl, iunit25, nou, bou, nbuf, ibufdim,
                                  myhost, myid, mycomm, mblk2nd, maxbl);
            writ_buf_ns::writ_buf(nbl, iunit9, nou, bou, nbuf, ibufdim,
                                  myhost, myid, mycomm, mblk2nd, maxbl);
        }
    }

    // append dynamic patch data to regular patch data
    // (see top of subroutine patcher for definition of iindex array)
    len0 = 0;
    if (ninter0 > 0) {
        nfb0 = iindex(ninter0, 1);
        len0 = iindex(ninter0, 2*nfb0+4)
             + iindex(ninter0, 2*nfb0+5) - 1;
    }
    for (n = 1; n <= nintr; n++) {
        iindex(n+ninter0, 1) = iindx(n, 1);
        nfb = iindx(n, 1);
        if (myid == myhost) {
            for (ll = 1; ll <= nfb; ll++) {
                iindex(n+ninter0, 1+ll) = iindx(n, 1+ll);
            }
            iindex(n+ninter0, nfb+2) = iindx(n, nfb+2);
            for (ll = 1; ll <= nfb; ll++) {
                iindex(n+ninter0, nfb+2+ll) = iindx(n, nfb+2+ll);
            }
            iindex(n+ninter0, 2*nfb+3) = iindx(n, 2*nfb+3);
            iindex(n+ninter0, 2*nfb+4) = iindx(n, 2*nfb+4);
            if (n == 1) {
                iindex(n+ninter0, 2*nfb+5) = len0 + 1;
            } else {
                nfb1 = iindex(n+ninter0-1, 1);
                iindex(n+ninter0, 2*nfb+5) = iindex(n+ninter0-1, 2*nfb1+4)
                                           + iindex(n+ninter0-1, 2*nfb1+5);
            }
            iindex(n+ninter0, 2*nfb+6) = iindx(n, 2*nfb+6);
            iindex(n+ninter0, 2*nfb+7) = iindx(n, 2*nfb+7);
            iindex(n+ninter0, 2*nfb+8) = iindx(n, 2*nfb+8);
            iindex(n+ninter0, 2*nfb+9) = iindx(n, 2*nfb+9);
        }

        lst = iindx(n, 2*nfb+5);
        len = lst + iindx(n, 2*nfb+4) - 1;

        for (nnn = lst; nnn <= len; nnn++) {
            nblkpt(len0+nnn) = mblkpt(nnn);
        }
        for (ll = 1; ll <= 2; ll++) {
            for (nnn = lst; nnn <= len; nnn++) {
                windex(len0+nnn, ll) = windx(nnn, ll);
            }
        }

        if (myid == myhost) {
            if (n == nintr && (len0+len-1) > maxxe) {
                // trim grdmov to null-terminated string for output
                char grdmov_str[61];
                std::memcpy(grdmov_str, grdmov, 60);
                grdmov_str[60] = '\0';
                // trim trailing spaces
                for (int ii = 59; ii >= 0; ii--) {
                    if (grdmov_str[ii] == ' ') grdmov_str[ii] = '\0';
                    else break;
                }
                fortran_write_unit(11, " program terminated in dynamic patching routines - see file %-60s\n", grdmov_str);
                nou(4) = std::min(nou(4)+1, ibufdim);
                std::snprintf(bou(nou(4), 4), 120, " ");
                nou(4) = std::min(nou(4)+1, ibufdim);
                std::snprintf(bou(nou(4), 4), 120,
                    " stopping...parameter maxxe is too small...must be at least%6d",
                    len0+len-1);
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
            for (ll = 1; ll <= nfb; ll++) {
                dthetxx(n+ninter0, ll) = dthetx(n, ll);
                dthetyy(n+ninter0, ll) = dthety(n, ll);
                dthetzz(n+ninter0, ll) = dthetz(n, ll);
                iindex(n+ninter0, 2*nfb+ll+9) = iindx(n, 2*nfb+ll+9);
                iindex(n+ninter0, 3*nfb+ll+9) = iindx(n, 3*nfb+ll+9);
                iindex(n+ninter0, 4*nfb+ll+9) = iindx(n, 4*nfb+ll+9);
                iindex(n+ninter0, 5*nfb+ll+9) = iindx(n, 5*nfb+ll+9);
            }
        }
    }

    if (myid == myhost) {

        ninter = ninter0 + nintr;

        if (ioutpt > 0) {

            // check iindex array
            fortran_write_unit(25, "\n");
            fortran_write_unit(25, "\n");
            fortran_write_unit(25, " ***** IINDEX ARRAY *****\n");
            fortran_write_unit(25, "\n");
            fortran_write_unit(25, " ninter: %d\n", ninter);
            fortran_write_unit(25, "\n");
            for (n = 1; n <= ninter; n++) {
                fortran_write_unit(25, "\n");
                fortran_write_unit(25, " interp. no.: %d\n", n);
                fortran_write_unit(25, "\n");
                fortran_write_unit(25, " nfb: %d\n", iindex(n, 1));
                nfb = iindex(n, 1);
                {
                    // write ' from block(s): ' followed by iindex(n,l+1) for l=1,nfb
                    FILE* f25 = fortran_get_unit(25);
                    fprintf(f25, " from block(s): ");
                    for (l = 1; l <= nfb; l++) {
                        fprintf(f25, " %d", iindex(n, l+1));
                    }
                    fprintf(f25, "\n");
                }
                fortran_write_unit(25, " to: %d\n", iindex(n, nfb+2));
                {
                    // write ' topology (from): ' followed by iindex(n,l+nfb+2) for l=1,nfb
                    FILE* f25 = fortran_get_unit(25);
                    fprintf(f25, " topology (from): ");
                    for (l = 1; l <= nfb; l++) {
                        fprintf(f25, " %d", iindex(n, l+nfb+2));
                    }
                    fprintf(f25, "\n");
                }
                fortran_write_unit(25, " topology (  to): %d\n", iindex(n, 2*nfb+3));
                fortran_write_unit(25, " number of points: %d\n", iindex(n, 2*nfb+4));
                fortran_write_unit(25, " starting index: %d\n", iindex(n, 2*nfb+5));
                fortran_write_unit(25, " xie range: %d %d\n",
                                   iindex(n, 2*nfb+6), iindex(n, 2*nfb+7));
                fortran_write_unit(25, " eta range: %d %d\n",
                                   iindex(n, 2*nfb+8), iindex(n, 2*nfb+9));
                for (l = 1; l <= nfb; l++) {
                    fortran_write_unit(25, " xie search range in from block %d : %d %d\n",
                                       l, iindex(n, 2*nfb+9+l), iindex(n, 3*nfb+9+l));
                    fortran_write_unit(25, " eta search range, from block %d : %d %d\n",
                                       l, iindex(n, 4*nfb+9+l), iindex(n, 5*nfb+9+l));
                }
            }
            fortran_write_unit(25, "\n");
            nfb = iindx(nintr, 1);
            lst = iindx(nintr, 2*nfb+5);
            len = lst + iindx(nintr, 2*nfb+4) - 1;
            fortran_write_unit(25, "  minimum dimension for parameter MAXXE in CFL3D:%6d\n", len);
            int iunit25 = 25;
            my_flush_ns::my_flush(iunit25);
        }

    }

    return;
}

void patcher(int& nbl, FortranArray2DRef<int> lw, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& ncall, int& ioutpt, int& it_thro, int& maxbl, int& msub1, int& intmx, int& mxxe, int& mptch, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<double> windx, int& nintr, FortranArray2DRef<int> iindx, FortranArray1DRef<int> llimit, FortranArray1DRef<int> iitmax, FortranArray1DRef<int> mmcxie, FortranArray1DRef<int> mmceta, FortranArray1DRef<int> ncheck, FortranArray1DRef<int> iifit, FortranArray1DRef<int> mblkpt, FortranArray1DRef<int> iic0, FortranArray1DRef<int> iiorph, FortranArray1DRef<int> iitoss, FortranArray1DRef<int> ifiner, FortranArray2DRef<double> factjlo, FortranArray2DRef<double> factjhi, FortranArray2DRef<double> factklo, FortranArray2DRef<double> factkhi, FortranArray2DRef<double> dx, FortranArray2DRef<double> dy, FortranArray2DRef<double> dz, FortranArray2DRef<double> dthetx, FortranArray2DRef<double> dthety, FortranArray2DRef<double> dthetz, FortranArray2DRef<int> isav_dpat, FortranArray3DRef<int> isav_dpat_b, FortranArray3DRef<double> xte, FortranArray3DRef<double> yte, FortranArray3DRef<double> zte, FortranArray3DRef<double> xmi, FortranArray3DRef<double> ymi, FortranArray3DRef<double> zmi, FortranArray3DRef<double> xmie, FortranArray3DRef<double> ymie, FortranArray3DRef<double> zmie, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<int> jte, FortranArray1DRef<int> kte, FortranArray3DRef<double> sxie, FortranArray3DRef<double> seta, FortranArray2DRef<double> sxie2, FortranArray2DRef<double> seta2, FortranArray2DRef<double> xie2s, FortranArray2DRef<double> eta2s, FortranArray1DRef<double> temp, FortranArray2DRef<double> x2, FortranArray2DRef<double> y2, FortranArray2DRef<double> z2, FortranArray1DRef<int> nblk1, FortranArray1DRef<int> nblk2, FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm, FortranArray2DRef<double> x1, FortranArray2DRef<double> y1, FortranArray2DRef<double> z1, FortranArray1DRef<int> lout, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, FortranArray1DRef<int> ireq_ar, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    patcher_ns::patcher(nbl, lw, w, mgwk, wk, nwork, ncall, ioutpt, it_thro,
                        maxbl, msub1, intmx, mxxe, mptch,
                        jdimg, kdimg, idimg, windx, nintr, iindx,
                        llimit, iitmax, mmcxie, mmceta, ncheck, iifit,
                        mblkpt, iic0, iiorph, iitoss, ifiner,
                        factjlo, factjhi, factklo, factkhi,
                        dx, dy, dz, dthetx, dthety, dthetz,
                        isav_dpat, isav_dpat_b,
                        xte, yte, zte, xmi, ymi, zmi, xmie, ymie, zmie,
                        jjmax1, kkmax1, jimage, kimage,
                        xorig, yorig, zorig,
                        jte, kte, sxie, seta, sxie2, seta2, xie2s, eta2s,
                        temp, x2, y2, z2, nblk1, nblk2, jmm, kmm,
                        x1, y1, z1, lout, xif1, xif2, etf1, etf2,
                        ireq_ar, myid, myhost, mycomm, mblk2nd,
                        nou, bou, nbuf, ibufdim);
}

} // namespace dynptch_ns
