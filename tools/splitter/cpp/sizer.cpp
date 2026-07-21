// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "sizer.h"
#include "runtime/fortran_io.h"
#include "umalloc.h"
#include "readkey.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <vector>

namespace sizer_ns {

// Helper: read a60 format from unit into char[81] buffer
static void read_a60(int unit, char* buf) {
    FILE* fp = fortran_get_unit(unit);
    char linebuf[1024];
    memset(buf, ' ', 80); buf[80] = '\0';
    if (!fp) return;
    if (fgets(linebuf, sizeof(linebuf), fp)) {
        int len = (int)strlen(linebuf);
        if (len > 0 && linebuf[len-1] == '\n') { linebuf[len-1] = '\0'; len--; }
        int copy = len < 60 ? len : 60;
        memcpy(buf, linebuf, copy);
        buf[copy] = '\0';
    }
}

// Helper: skip a line (blank read)
static void skip_line(int unit) {
    FILE* fp = fortran_get_unit(unit);
    if (!fp) return; // Fortran: READ on unopened unit is a no-op (or error)
    char linebuf[1024];
    fgets(linebuf, sizeof(linebuf), fp);
}

// Helper: Fortran-style string equality (trim trailing spaces)
static bool str_eq_f(const char* s, const char* ref) {
    int slen = (int)strlen(s);
    while (slen > 0 && s[slen-1] == ' ') slen--;
    int rlen = (int)strlen(ref);
    if (slen != rlen) return false;
    return strncmp(s, ref, rlen) == 0;
}

// Helper: trim trailing spaces and null-terminate for use as C string
static void trim_to_cstr(const char* src, char* dst, int maxlen) {
    strncpy(dst, src, maxlen);
    dst[maxlen] = '\0';
    int l = (int)strlen(dst);
    while (l > 0 && dst[l-1] == ' ') l--;
    dst[l] = '\0';
}

// Helper: read one block dimension line from unit (separate READ per block)
static void read_block_dims(int unit, int& id, int& jd, int& kd) {
    fortran_read_list(unit, &id, &jd, &kd);
}

// Helper: read one split-input block number from unit (separate READ per split)
static void read_split_blk(int unit, int& isbloc) {
    fortran_read_list(unit, &isbloc);
}

// Helper: read one split-input direction from unit (separate READ per split)
static void read_split_dir(int unit, int& ndir) {
    fortran_read_list(unit, &ndir);
}

void sizer(int& mbloc, int& maxseg, int& npmax, int& mxbli, int& intmax, int& nsub1, int& msegt, int& msegn, int& msplt, int& mtot, int& mbloc0, int& msplt0, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    char cflinp[81], cflout[81], roninp[81], ronout[81], sdgridin[81];
    char sdgridout[81], dovrlap[81], dpatch[81], dresid[81], gridin[81], gridout[81];
    int stats, memuse;
    int icflin, ironin, ibin, isdin;
    int icflout, ironout, ibout, isdout;
    int nbl, nout, nn, nb, no;
    int nskip, n, ibloc;
    int isblocdum;
    int iunit5sav, ititr, myid;
    int ii, jj, kk;
    int isp1, jsp1, ksp1, mspltnn;

    FortranArray1D<int> idimg, isblocin, isplits, jdimg, jsplits, kdimg, ksplits, ndirin, ndirin0;

    // common /unit5/ iunit5
    int& iunit5 = cmn_unit5.iunit5;

    memuse = 0;

    // allocate( idimg(mbloc0), stat=stats )
    idimg.allocate(mbloc0); stats = 0;
    { int one = 1; char txt[6] = "idimg"; umalloc_ns::umalloc(mbloc0, one, reinterpret_cast<char(&)[1]>(txt), memuse, stats); }

    // allocate( isblocin(msplt0*mbloc0), stat=stats )
    isblocin.allocate(msplt0 * mbloc0); stats = 0;
    { int sz = msplt0 * mbloc0; int one = 1; char txt[9] = "isblocin"; umalloc_ns::umalloc(sz, one, reinterpret_cast<char(&)[1]>(txt), memuse, stats); }

    // allocate( isplits(mbloc0), stat=stats )
    isplits.allocate(mbloc0); stats = 0;
    { int one = 1; char txt[8] = "isplits"; umalloc_ns::umalloc(mbloc0, one, reinterpret_cast<char(&)[1]>(txt), memuse, stats); }

    // allocate( jdimg(mbloc0), stat=stats )
    jdimg.allocate(mbloc0); stats = 0;
    { int one = 1; char txt[6] = "jdimg"; umalloc_ns::umalloc(mbloc0, one, reinterpret_cast<char(&)[1]>(txt), memuse, stats); }

    // allocate( jsplits(mbloc0), stat=stats )
    jsplits.allocate(mbloc0); stats = 0;
    { int one = 1; char txt[8] = "jsplits"; umalloc_ns::umalloc(mbloc0, one, reinterpret_cast<char(&)[1]>(txt), memuse, stats); }

    // allocate( kdimg(mbloc0), stat=stats )
    kdimg.allocate(mbloc0); stats = 0;
    { int one = 1; char txt[6] = "kdimg"; umalloc_ns::umalloc(mbloc0, one, reinterpret_cast<char(&)[1]>(txt), memuse, stats); }

    // allocate( ksplits(mbloc0), stat=stats )
    ksplits.allocate(mbloc0); stats = 0;
    { int one = 1; char txt[8] = "ksplits"; umalloc_ns::umalloc(mbloc0, one, reinterpret_cast<char(&)[1]>(txt), memuse, stats); }

    // allocate( ndirin(msplt0*mbloc0), stat=stats )
    ndirin.allocate(msplt0 * mbloc0); stats = 0;
    { int sz = msplt0 * mbloc0; int one = 1; char txt[7] = "ndirin"; umalloc_ns::umalloc(sz, one, reinterpret_cast<char(&)[1]>(txt), memuse, stats); }

    // allocate( ndirin0(msplt0*mbloc0), stat=stats )
    ndirin0.allocate(msplt0 * mbloc0); stats = 0;
    { int sz = msplt0 * mbloc0; int one = 1; char txt[8] = "ndirin0"; umalloc_ns::umalloc(sz, one, reinterpret_cast<char(&)[1]>(txt), memuse, stats); }

    skip_line(iunit5);                  // read(iunit5,*)
    read_a60(iunit5, cflinp);           // read(iunit5,'(a60)') cflinp
    read_a60(iunit5, roninp);           // read(iunit5,'(a60)') roninp
    read_a60(iunit5, gridin);           // read(iunit5,'(a60)') gridin
    read_a60(iunit5, sdgridin);         // read(iunit5,'(a60)') sdgridin
    skip_line(iunit5);                  // read(iunit5,*)
    fortran_read_list(iunit5, &icflin, &ironin, &ibin, &isdin); // read(iunit5,*) icflin,ironin,ibin,isdin
    icflin = abs(icflin);
    skip_line(iunit5);                  // read(iunit5,*)
    read_a60(iunit5, cflout);           // read(iunit5,'(a60)') cflout
    read_a60(iunit5, ronout);           // read(iunit5,'(a60)') ronout
    read_a60(iunit5, gridout);          // read(iunit5,'(a60)') gridout
    read_a60(iunit5, sdgridout);        // read(iunit5,'(a60)') sdgridout
    skip_line(iunit5);                  // read(iunit5,*)
    fortran_read_list(iunit5, &icflout, &ironout, &ibout, &isdout); // read(iunit5,*) icflout,ironout,ibout,isdout
    skip_line(iunit5);                  // read(iunit5,*)

    //
    // read the either the cfl3d input file or the (plot3d) grid
    // file to get the number of blocks and grid dimensions
    //
    if (str_eq_f(cflinp, "null")) {
        if (!str_eq_f(gridin, "null")) {
            char gridin_cstr[81];
            trim_to_cstr(gridin, gridin_cstr, 80);
            if (ibin == 1) {
                // the following line is needed to use ieee binary file
                // call asnfile (gridin, '-F f77 -N ieee' , IER)
                fortran_open_unit(30, gridin_cstr, "rb");
            } else {
                fortran_open_unit(30, gridin_cstr, "r");
            }
            rewind(fortran_get_unit(30)); // rewind 30
        }
        if (ibin == 1) {
            // read(30) nbl  -- unformatted binary with Fortran record markers
            {
                FILE* fp = fortran_get_unit(30);
                int32_t reclen;
                fread(&reclen, sizeof(int32_t), 1, fp);
                fread(&nbl, sizeof(int), 1, fp);
                fread(&reclen, sizeof(int32_t), 1, fp);
            }
            // read(30) (idimg(n),jdimg(n),kdimg(n),n=1,nbl)
            {
                FILE* fp = fortran_get_unit(30);
                int32_t reclen;
                fread(&reclen, sizeof(int32_t), 1, fp);
                for (int n2 = 1; n2 <= nbl; n2++) {
                    fread(&idimg(n2), sizeof(int), 1, fp);
                    fread(&jdimg(n2), sizeof(int), 1, fp);
                    fread(&kdimg(n2), sizeof(int), 1, fp);
                }
                fread(&reclen, sizeof(int32_t), 1, fp);
            }
        } else {
            // read(30,*) nbl
            fortran_read_list(30, &nbl);
            // read(30,*) (idimg(n),jdimg(n),kdimg(n),n=1,nbl)
            // Single implied-DO read: all 3*nbl values from one record
            // Use a temporary interleaved buffer then unpack
            {
                std::vector<int> tmp(3 * nbl);
                fortran_read_list_array(30, tmp.data(), 3 * nbl);
                for (int n2 = 1; n2 <= nbl; n2++) {
                    idimg(n2) = tmp[3*(n2-1)+0];
                    jdimg(n2) = tmp[3*(n2-1)+1];
                    kdimg(n2) = tmp[3*(n2-1)+2];
                }
            }
        }
    } else {
        // open(unit=10,file=cflinp,form='formatted',status='old')
        {
            char cflinp_cstr[81];
            trim_to_cstr(cflinp, cflinp_cstr, 80);
            fortran_open_unit(10, cflinp_cstr, "r");
        }
        nskip = 14;
        n = 1; while (n <= nskip) { skip_line(10); n++; } // do n=1,nskip; read(10,*); end do

        // check for keyword input
        iunit5sav = iunit5;
        iunit5 = 10;
        {
            int iunit11 = 99;
            int ierrflg_val = -1;
            readkey_ns::readkey(ititr, myid, ibufdim, nbuf, bou, nou, iunit11, ierrflg_val);
        }
        iunit5 = iunit5sav;

        if (ititr == 0) { nskip = 8; } else { nskip = 7; }
        n = 1; while (n <= nskip) { skip_line(10); n++; } // do n=1,nskip; read(10,*); end do

        fortran_read_list(10, &nbl); // read(10,*) nbl
        nbl = abs(nbl);
        skip_line(10); // read(10,*)
        ibloc = 1; while (ibloc <= nbl) { skip_line(10); ibloc++; } // do ibloc=1,nbl; read(10,*); end do
        skip_line(10); // read(10,*)
        // do ibloc=1,nbl; read(10,*) idimg(ibloc),jdimg(ibloc),kdimg(ibloc); end do
        // Each iteration is a separate READ statement (one record per block)
        ibloc = 1;
        while (ibloc <= nbl) {
            read_block_dims(10, idimg(ibloc), jdimg(ibloc), kdimg(ibloc));
            ibloc++;
        }
    }

    fortran_read_list(iunit5, &nout); // read(iunit5,*) nout

    // shortcut to split all blocks the same as the single block input
    if (nout < 0) {
        nout = abs(nout);
        // Each iteration: 3 separate READ statements (one record each)
        no = 1;
        while (no <= nout) {
            read_split_blk(iunit5, isblocdum);   // read(iunit5,*) isblocdum
            read_split_dir(iunit5, ndirin0(no)); // read(iunit5,*) ndirin0(no)
            skip_line(iunit5);                    // read(iunit5,*)
            no++;
        }
        nn = 0;
        for (nb = 1; nb <= nbl; nb++) {
            for (no = 1; no <= nout; no++) {
                nn = nn + 1;
                isblocin(nn) = nb;
                ndirin(nn)   = ndirin0(no);
            }
        }
        nout = nn;
    } else {
        // Each iteration: 3 separate READ statements (one record each)
        nn = 1;
        while (nn <= nout) {
            read_split_blk(iunit5, isblocin(nn)); // read(iunit5,*) isblocin(nn)
            read_split_dir(iunit5, ndirin(nn));   // read(iunit5,*) ndirin(nn)
            skip_line(iunit5);                     // read(iunit5,*)
            nn++;
        }
    }

    //
    // count the number of splits in each direction of each block
    // and determine the maximum number of splits per block
    //
    msplt = 0;
    for (nn = 1; nn <= nbl; nn++) {
        isplits(nn) = 0;
        jsplits(nn) = 0;
        ksplits(nn) = 0;
        for (no = 1; no <= nout; no++) {
            if (isblocin(no) == nn) {
                if (ndirin(no) == 1) {
                    isplits(nn) = isplits(nn) + 1;
                } else if (ndirin(no) == 2) {
                    jsplits(nn) = jsplits(nn) + 1;
                } else if (ndirin(no) == 3) {
                    ksplits(nn) = ksplits(nn) + 1;
                }
            }
        }
        mspltnn = isplits(nn) + jsplits(nn) + ksplits(nn);
        if (mspltnn > msplt) msplt = mspltnn;
    }

    // use a min. of msplt=1 to avoid allocating 0 bytes
    msplt = std::max(1, msplt);

    //
    // determine the number of grid points after splitting
    //
    npmax = 0;
    for (nn = 1; nn <= nbl; nn++) {
        npmax = npmax + jdimg(nn) * kdimg(nn) * idimg(nn);
        if (isplits(nn) > 0) {
            for (ii = 1; ii <= isplits(nn); ii++) {
                npmax = npmax + jdimg(nn) * kdimg(nn);
            }
            idimg(nn) = idimg(nn) + isplits(nn);
        }
        if (jsplits(nn) > 0) {
            for (jj = 1; jj <= jsplits(nn); jj++) {
                npmax = npmax + idimg(nn) * kdimg(nn);
            }
            jdimg(nn) = jdimg(nn) + jsplits(nn);
        }
        if (ksplits(nn) > 0) {
            for (kk = 1; kk <= ksplits(nn); kk++) {
                npmax = npmax + jdimg(nn) * idimg(nn);
            }
            kdimg(nn) = kdimg(nn) + ksplits(nn);
        }
    }

    //
    // determine number of blocks after splitting
    //
    mbloc = 0;
    for (nn = 1; nn <= nbl; nn++) {
        isp1 = isplits(nn) + 1;
        jsp1 = jsplits(nn) + 1;
        ksp1 = ksplits(nn) + 1;
        mbloc = mbloc + isp1 * jsp1 * ksp1;
    }

    // the following parameters should be quite sufficient for any case
    mxbli  = 2500;
    intmax = 1500;
    msegt  = 75;
    maxseg = 200;

    // set the rest of the parameters
    nsub1 = mbloc;
    msegn = maxseg * 6;
    mtot  = mbloc * msegn;

    { FILE* fp = fortran_get_unit(iunit5); if (fp) rewind(fp); }
    { FILE* fp = fortran_get_unit(10);     if (fp) rewind(fp); }
    { FILE* fp = fortran_get_unit(30);     if (fp) rewind(fp); }

    // free up memory
    isplits.deallocate();
    jsplits.deallocate();
    ksplits.deallocate();
    idimg.deallocate();
    jdimg.deallocate();
    kdimg.deallocate();
    ndirin0.deallocate();
    ndirin.deallocate();
    isblocin.deallocate();

    return;
}

} // namespace sizer_ns
