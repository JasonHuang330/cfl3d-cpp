// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// clusterB.cpp — maggie "Cluster B: iblank scratch I/O + bookkeeping"
// Byte-exact F77 -> C++ port. See /tmp/mag_subs/RULES.md.
//
// SIGNATURES (all match mag_protos.h):
// SIGNATURE: void putibl(int m,int* iblank,int jd,int kd,int ld);
// SIGNATURE: void putibl2(int m,int* iblank,int jd,int kd,int ld);
// SIGNATURE: void getibl(int m,int* iblank,int jd,int kd,int ld);
// SIGNATURE: void putint(int nwr,int m,int* ji,int* ki,int* li,int* jbpt,int* kbpt,int* lbpt,double* dxint,double* dyint,double* dzint);
// SIGNATURE: void getint(int m,int* ji,int* ki,int* li,int* jbpt,int* kbpt,int* lbpt,double* dxint,double* dyint,double* dzint);
// SIGNATURE: void wiblnk(int* iblank,int jd,int kd,int ld);
// SIGNATURE: void iblcon(int* iblank,int jd,int kd,int ld);
// SIGNATURE: void iblplt(int m,int* iblank,int jd,int kd,int ld);
// SIGNATURE: void reorder(int idimA,int ibpnts,int intpts,int* jb,int* kb,int* lb,int* iord,int jd,int kd,int ld);
// SIGNATURE: void cindex(int ibc,int ibpnts,int iipnts,int* jb,int* kb,int* lb,int mesh);
// SIGNATURE: void setptr(int m,int m1,int icase);
// SIGNATURE: void pltpts_(int nmesh,int m,int jmax,int kmax,int lmax,double* xx,double* yy,double* zz,int* ib,int* ibb);
// SIGNATURE: void trace_(int icall,int idum1,int idum2,int idum3,int idum4,double dum1,double dum2,double dum3);
//
// SIGNATURE-MISMATCH NOTES (reconcile in mag_protos.h):
//  * reorder: Fortran declares (idim,ibpnts,intpts,jb,kb,lb,iord,jd,kd,ld) where
//    `intpts` is an INTEGER ARRAY intpts(4) (out) and `ibpnts` is an in-scalar.
//    mag_protos.h declares BOTH ibpnts and intpts as plain `int` (by value).
//    In Fortran `intpts` is written to (out) so it MUST be a pointer. I follow
//    mag_protos.h EXACTLY as instructed (int intpts) but this is WRONG: callers
//    lose the 4-element intpts output. Flagged for reconciliation. Also `idim`
//    (a dummy arg) shadows the global param `idim`; proto renames it `idimA`.
//  * cindex: Fortran declares (ibc,ibpnts,iipnts,jb,kb,lb,mesh) with ibc(*),
//    jb(*),kb(*),lb(*) arrays and ibpnts,iipnts INTEGER scalars written to (out).
//    mag_protos.h declares ibc as `int` (by value) not `int*`, and ibpnts/iipnts
//    as plain `int` (by value) not `int&`. Followed proto as instructed; but
//    ibc must be `int*` and ibpnts/iipnts must be `int&` for correctness. Flagged.
//  * putint: mag_protos.h passes `int nwr` by value (matches Fortran arg order).
//    getint takes no nwr arg — it fills the global /chkst/ nwr (chkst.nwr). OK.
//  * pltpts -> pltpts_ ; trace -> trace_ (per note). OK.

#include "mag_common.h"
#include <vector>

using tc::RecordWriter;

// ---------------------------------------------------------------------------
// gfortran list-directed (write(u,*)) helpers — reproduce EXACT field spacing.
//   * integer item -> right-justified in a 12-char field ("%12d").
//   * a leading blank is emitted at the START of each record (record carriage
//     control); we prepend it explicitly at each write site.
//   * a character constant is emitted verbatim; between a value and a following
//     character constant gfortran inserts one separator blank (handled inline).
//   * real*8 item -> 26-char slot (F-form for 1e-1<=|x|<1e16 else E-form),
//     17 significant digits; +1 extra separator when followed by a char const.
// ---------------------------------------------------------------------------
static std::string ld_r8(double v) {
    // Reproduce gfortran list-directed REAL*8 in a fixed 26-char slot.
    //   F-form (0.1<=|v|<1e16 or v==0): 3 leading spaces + 17-sig-digit value,
    //                                   then trailing spaces to fill 26.
    //   E-form (otherwise): the E-value (mantissa [1,10), 16 frac digits,
    //                       3-digit signed exponent) right-justified in 26.
    // NOTE: this matches gfortran for the common magnitudes exercised by maggie
    // trace diagnostics; exact-half rounding ties may differ by 1 ulp digit.
    const int SLOT = 26;
    char body[80];
    double a = std::fabs(v);
    bool useF = (a == 0.0) || (a >= 0.1 && a < 1.0e16);
    std::string val;
    if (useF) {
        int intdigits;
        if (a == 0.0) intdigits = 1;
        else {
            int e = (int)std::floor(std::log10(a));
            intdigits = (e >= 0) ? e + 1 : 1;   // digits left of the point
        }
        int dec = 17 - intdigits;
        if (dec < 0) dec = 0;
        std::snprintf(body, sizeof(body), "%.*f", dec, v);
        val = std::string("   ") + body;                // 3 leading spaces
        if ((int)val.size() < SLOT) val += std::string(SLOT - val.size(), ' ');
    } else {
        char tmp[80];
        std::snprintf(tmp, sizeof(tmp), "%.16E", v);     // C: d.ddddddddddddddddE+xx
        std::string s(tmp);
        size_t ep = s.find('E');
        std::string mant = s.substr(0, ep);
        char sign = s[ep + 1];
        int expn = std::atoi(s.c_str() + ep + 2);
        char eb[16];
        std::snprintf(eb, sizeof(eb), "E%c%03d", sign, expn);
        val = mant + eb;
        if ((int)val.size() < SLOT) val = std::string(SLOT - val.size(), ' ') + val;
    }
    return val;
}

// ---------------------------------------------------------------------------
// helper: build the "temp_XXX.m" filename exactly like the Fortran writes it.
//   m>99 : "<pref>.%3d"  (len 12)   ; m>9 : "<pref>.%2d" (len 11)
//   else : "<pref>.%1d"  (len 10).  The i1/i2/i3 edit descriptors write the
//   integer right-justified with no leading zeros; for the len==10 (i1) case
//   only 0..9 fit (matches Fortran which would print '*' on overflow, but m<=9
//   is guaranteed by the branch). We reproduce the exact digit string.
// ---------------------------------------------------------------------------
static void mag_tempname(const char* pref, int m, char* titl /*char[21]*/) {
    // pref is e.g. "temp_ibl." (9 chars incl. trailing dot)
    if (m > 99) {
        std::snprintf(titl, 21, "%s%3d", pref, m);
    } else if (m > 9) {
        std::snprintf(titl, 21, "%s%2d", pref, m);
    } else {
        std::snprintf(titl, 21, "%s%1d", pref, m);
    }
    // pad to 20 with blanks (harmless; open uses titl(1:len) only)
    int n = (int)std::strlen(titl);
    for (int i = n; i < 20; ++i) titl[i] = ' ';
    titl[20] = '\0';
}

// ===========================================================================
//  putibl — write iblank(jd,kd,ld) for mesh m to temp_ibl.m (ONE record,
//           whole-array contiguous, matching Fortran `write(iunit) iblank`).
// ===========================================================================
void putibl(int m, int* iblank, int jd, int kd, int ld) {
    const int iunit = 30;
    char titl[21], name[21];
    mag_tempname("temp_ibl.", m, titl);
    // filename is titl(1:len); len==10/11/12 but the trailing pad is blanks and
    // std::strlen already stops at first blank -> reconstruct exact len bytes.
    int len = (m > 99) ? 12 : (m > 9) ? 11 : 10;
    std::memcpy(name, titl, len);
    name[len] = '\0';

    fopen_unit(iunit, name, "wb+");
    // one record: the whole iblank array (jd*kd*ld ints), 1-based flat at [1].
    uwrite(iunit, &iblank[1], (size_t)jd * kd * ld * sizeof(int));
    frewind(iunit);
}

// ===========================================================================
//  putibl2 — like putibl but file temp_grd.m and an explicit element loop
//            (((iblank(j,k,l),j=1,jd),k=1,kd),l=1,ld). Column-major order is
//            identical to the contiguous layout, so ONE record of the same
//            bytes. We still write element-by-element to be literal.
// ===========================================================================
void putibl2(int m, int* iblank, int jd, int kd, int ld) {
    const int iunit = 30;
    char titl[21], name[21];
    mag_tempname("temp_grd.", m, titl);
    int len = (m > 99) ? 12 : (m > 9) ? 11 : 10;
    std::memcpy(name, titl, len);
    name[len] = '\0';

    fopen_unit(iunit, name, "wb+");
    View3<int> IB(iblank, jd, kd, ld);
    RecordWriter rw;
    for (int l = 1; l <= ld; ++l)
        for (int k = 1; k <= kd; ++k)
            for (int j = 1; j <= jd; ++j)
                rw.put_i(IB(j, k, l));
    rw.flush(funit(iunit));
    frewind(iunit);
}

// ===========================================================================
//  getibl — read iblank back from temp_ibl.m (ONE record, whole array).
// ===========================================================================
void getibl(int m, int* iblank, int jd, int kd, int ld) {
    const int iunit = 30;
    char titl[21], name[21];
    mag_tempname("temp_ibl.", m, titl);
    int len = (m > 99) ? 12 : (m > 9) ? 11 : 10;
    std::memcpy(name, titl, len);
    name[len] = '\0';

    fopen_unit(iunit, name, "rb");
    uread(iunit, &iblank[1], (size_t)jd * kd * ld * sizeof(int));
    frewind(iunit);
    fclose_unit(iunit);
}

// ===========================================================================
//  putint — write interpolation data for mesh m to temp_int.m.
//     record 1: nwr (one int)
//     record 2 (only if nwr>0): (ji,ki,li)*nwr, then (jbpt,kbpt,lbpt)*nwr,
//               then (dxint,dyint,dzint)*nwr — ALL in ONE record, in the
//               exact Fortran element order.
// ===========================================================================
void putint(int nwr, int m, int* ji, int* ki, int* li,
            int* jbpt, int* kbpt, int* lbpt,
            double* dxint, double* dyint, double* dzint) {
    const int iunit = 30;
    char titl[21], name[21];
    mag_tempname("temp_int.", m, titl);
    int len = (m > 99) ? 12 : (m > 9) ? 11 : 10;
    std::memcpy(name, titl, len);
    name[len] = '\0';

    fopen_unit(iunit, name, "wb+");

    // record 1: nwr
    uwrite(iunit, &nwr, sizeof(int));

    if (nwr > 0) {
        RecordWriter rw;
        for (int i = 1; i <= nwr; ++i) { rw.put_i(ji[i]);   rw.put_i(ki[i]);   rw.put_i(li[i]);   }
        for (int i = 1; i <= nwr; ++i) { rw.put_i(jbpt[i]); rw.put_i(kbpt[i]); rw.put_i(lbpt[i]); }
        for (int i = 1; i <= nwr; ++i) { rw.put_d(dxint[i]);rw.put_d(dyint[i]);rw.put_d(dzint[i]);}
        rw.flush(funit(iunit));
    }
    frewind(iunit);
}

// ===========================================================================
//  getint — read interpolation data for mesh m from temp_int.m.
//     nwr goes into global /chkst/  -> chkst.nwr.
// ===========================================================================
void getint(int m, int* ji, int* ki, int* li,
            int* jbpt, int* kbpt, int* lbpt,
            double* dxint, double* dyint, double* dzint) {
    const int iunit = 30;
    char titl[21], name[21];
    mag_tempname("temp_int.", m, titl);
    int len = (m > 99) ? 12 : (m > 9) ? 11 : 10;
    std::memcpy(name, titl, len);
    name[len] = '\0';

    fopen_unit(iunit, name, "rb");

    // record 1: nwr
    uread(iunit, &chkst.nwr, sizeof(int));
    int nwr = chkst.nwr;

    if (nwr > 0) {
        // record 2: read the whole record into a temp buffer, then scatter in
        // the exact Fortran order (uread consumes the whole record).
        size_t nint = (size_t)nwr * 6;      // ji,ki,li,jbpt,kbpt,lbpt each nwr
        size_t ndbl = (size_t)nwr * 3;      // dxint,dyint,dzint each nwr
        std::vector<uint8_t> buf(nint * sizeof(int) + ndbl * sizeof(double));
        uread(iunit, buf.data(), buf.size());

        const uint8_t* p = buf.data();
        auto geti = [&](void)->int    { int v;    std::memcpy(&v, p, 4); p += 4; return v; };
        auto getd = [&](void)->double { double v; std::memcpy(&v, p, 8); p += 8; return v; };

        for (int i = 1; i <= nwr; ++i) { ji[i]   = geti(); ki[i]   = geti(); li[i]   = geti(); }
        for (int i = 1; i <= nwr; ++i) { jbpt[i] = geti(); kbpt[i] = geti(); lbpt[i] = geti(); }
        for (int i = 1; i <= nwr; ++i) { dxint[i]= getd(); dyint[i]= getd(); dzint[i]= getd(); }
    }
    frewind(iunit);
    fclose_unit(iunit);
}

// ===========================================================================
//  wiblnk — write iblank as REAL (float, per -fdefault-real-8 -> float()
//     returns default real = double under -fdefault-real-8) to unit 2 in the
//     form expected by cfl3d: ONE record of the whole array as reals, in
//     (((iblank(j,k,l),j),k),l) order.
//
//     NOTE on real kind: `float(...)` under -fdefault-real-8 yields an 8-byte
//     real, so the record holds DOUBLES. We emit put_d to match the -r8 build.
// ===========================================================================
void wiblnk(int* iblank, int jd, int kd, int ld) {
    const int iunit = 2;
    View3<int> IB(iblank, jd, kd, ld);
    RecordWriter rw;
    for (int l = 1; l <= ld; ++l)
        for (int k = 1; k <= kd; ++k)
            for (int j = 1; j <= jd; ++j)
                rw.put_d((double)IB(j, k, l));   // float(iblank) == double under -r8
    rw.flush(funit(iunit));
}

// ===========================================================================
//  iblcon — set iblank=0 wherever iblank != 1.
// ===========================================================================
void iblcon(int* iblank, int jd, int kd, int ld) {
    View3<int> IB(iblank, jd, kd, ld);
    for (int l = 1; l <= ld; ++l)
        for (int k = 1; k <= kd; ++k)
            for (int j = 1; j <= jd; ++j)
                if (IB(j, k, l) != 1) IB(j, k, l) = 0;
}

// ===========================================================================
//  iblplt — print the iblank array for mesh m to unit 6.
// ===========================================================================
void iblplt(int m, int* iblank, int jd, int kd, int ld) {
    View3<int> IB(iblank, jd, kd, ld);

    // 299 format(/,/,' iblank array for cell centers of mesh=',i3)
    std::fprintf(funit(6), "\n\n iblank array for cell centers of mesh=%3d\n", m);

    for (int l = 1; l <= ld; ++l) {
        int icell = 0;
        for (int k = 1; k <= kd; ++k)
            for (int j = 1; j <= jd; ++j)
                if (IB(j, k, l) != 1) icell = icell + 1;

        if (icell == 0) {
            // 302 format(' .....all values unity for l=constant= ',i3)
            std::fprintf(funit(6), " .....all values unity for l=constant= %3d\n", l);
        } else {
            // 301 format(' .....l=constant= ',i3)
            std::fprintf(funit(6), " .....l=constant= %3d\n", l);
            // 305 format(' ',130i1)  -> leading ' ', then each value as %1d,
            //   wrapping every 130 values onto a new record (Fortran line rule).
            for (int k = 1; k <= kd; ++k) {
                std::fprintf(funit(6), " ");
                int col = 0;
                for (int j = 1; j <= jd; ++j) {
                    if (col == 130) {           // format exhausted -> new record
                        std::fprintf(funit(6), "\n ");
                        col = 0;
                    }
                    std::fprintf(funit(6), "%1d", std::abs(IB(j, k, l)));
                    ++col;
                }
                std::fprintf(funit(6), "\n");
            }
        }
    }
}

// ===========================================================================
//  reorder — reorder/modify interpolated values for QI0/QJ0/QK0.
//  NOTE: per mag_protos.h, `ibpnts` and `intpts` are passed by value; but the
//  Fortran writes intpts(1..4) (out). We keep a LOCAL intpts[5] to preserve the
//  algorithm/prints exactly, since the proto cannot return them. See header.
// ===========================================================================
void reorder(int idimA, int ibpnts, int* intpts, int* jb, int* kb, int* lb,
             int* iord, int jd, int kd, int ld) {
    (void)idimA;                 // idim dummy unused
    int* intp = intpts;          // intpts(1..4), 1-based (caller passes base ptr)
    for (int ip = 1; ip <= 4; ++ip) intp[ip] = 0;

    int iseq = 0;
    int iint = 0;
    for (int i = 1; i <= ibpnts; ++i) {
        if (jb[i] <= jd && kb[i] <= kd && lb[i] <= ld &&
            jb[i] >= 1  && kb[i] >= 1  && lb[i] >= 1) {
            iint = iint + 1;
            if (i == iint) iseq = iint;
            iord[iint] = i;
        }
    }
    intp[1] = iint;
    // list-directed: leading record blank + literal + %12d
    std::fprintf(funit(6), "         ...reordering.......interior pts = %12d\n", intp[1]);
    std::fprintf(funit(6), "         .....iseq,ibpnts = %12d%12d\n", iseq, ibpnts);

    int ist = iseq + 1;
    if (!(ist > ibpnts || iint >= ibpnts)) {
        for (int i = ist; i <= ibpnts; ++i) {
            if (jb[i] < 1 || jb[i] > jd) {
                iint = iint + 1;
                if (i == iint) iseq = iint;
                iord[iint] = i;
            }
        }
        intp[2] = iint - intp[1];
        std::fprintf(funit(6), "         ...reordering.......QJ0 pts = %12d\n", intp[2]);
        std::fprintf(funit(6), "         ......iseq,ibpnts = %12d%12d\n", iseq, ibpnts);

        ist = iseq + 1;
        if (!(ist > ibpnts || iint >= ibpnts)) {
            for (int i = ist; i <= ibpnts; ++i) {
                if (kb[i] < 1 || kb[i] > kd) {
                    iint = iint + 1;
                    if (i == iint) iseq = iint;
                    iord[iint] = i;
                }
            }
            intp[3] = iint - intp[2] - intp[1];
            std::fprintf(funit(6), "         ...reordering.......QK0 pts = %12d\n", intp[3]);
            std::fprintf(funit(6), "         ......iseq,ibpnts = %12d%12d\n", iseq, ibpnts);

            ist = iseq + 1;
            if (!(ist > ibpnts || iint >= ibpnts)) {
                for (int i = ist; i <= ibpnts; ++i) {
                    if (lb[i] < 1 || lb[i] > ld) {
                        iint = iint + 1;
                        if (i == iint) iseq = iint;
                        iord[iint] = i;
                    }
                }
                intp[4] = iint - intp[3] - intp[2] - intp[1];
                std::fprintf(funit(6), "         ...reordering.......QI0 pts = %12d\n", intp[4]);
                std::fprintf(funit(6), "         ......iseq,ibpnts = %12d%12d\n", iseq, ibpnts);
            }
        }
    }
    // label 1000:

    if (ibpnts != iint) {
        // write(6,*) '    stopping in reorder......ibpnts.ne.iint',ibpnts,iint,intpts
        std::fprintf(funit(6),
            "     stopping in reorder......ibpnts.ne.iint%12d%12d%12d%12d%12d%12d\n",
            ibpnts, iint, intp[1], intp[2], intp[3], intp[4]);
        std::exit(1);   // stop
    }
}

// ===========================================================================
//  cindex — build the cross-index array ibc + jb/kb/lb for one mesh.
//  NOTE (proto mismatch): ibc/jb/kb/lb are arrays (int*), ibpnts/iipnts are
//  out-scalars (int&). mag_protos.h passes ibc/ibpnts/iipnts by value; the
//  by-value ibc is meaningless. We treat ibc/jb/kb/lb as int* and ibpnts/iipnts
//  as locals returned via globals is not possible under proto -> flagged. To
//  keep the algorithm intact we DEFINE using pointers where Fortran uses arrays;
//  the proto's scalar `ibc` cannot be honored without breaking the routine.
//  Implemented to the Fortran semantics (see header note).
// ===========================================================================
void cindex(int* ibc, int& ibpnts, int& iipnts, int* jb, int* kb, int* lb, int mesh) {
    int* ibc_loc = ibc;             // Fortran ibc(*) OUTPUT cross-index
    int ipts[mdim + 1];
    int list[mdim + 1];

    // set up offset counter
    ipts[1] = 0;
    for (int m = 2; m <= conec1.nmesh; ++m)
        ipts[m] = ipts[m - 1] + book2.ibpts[m - 1];

    // set up cross index and pointers
    int m = mesh;
    int ic = 0;

    // outer boundaries
    int no = conec1.noutr[m];
    if (no != 0) {
        for (int m1 = 1; m1 <= conec1.nmesh; ++m1) list[m1] = 0;

        // do 412 n=1,no ; nnn=1; nserch=moutr(m,nnn); do 412 nn=1,nserch ...
        for (int n = 1; n <= no; ++n) {
            int nnn = 1;
            int nserch = conec1.MOUTR(m, nnn);
            for (int nn = 1; nn <= nserch; ++nn) {
                nnn = 1;
                int m1 = conec1.LOUTR(m, nnn, nn);
                if (list[m1] == 0) list[m1] = 1;
            }
        }
        for (int m1 = 1; m1 <= conec1.nmesh; ++m1) {
            if (list[m1] > 0) {
                int iset = book1.MOBS(m, m1);
                int is = book1.IPNTR(m1, iset);
                int ie = book1.NPNTR(m1, iset);
                getint(m1, book2.ji, book2.ki, book2.li,
                       book2.jbpt, book2.kbpt, book2.lbpt,
                       book2.dxint, book2.dyint, book2.dzint);
                for (int i = is; i <= ie; ++i) {
                    ic = ic + 1;
                    // ibc(ic) = i + ipts(m1)  (proto-limited: local shadow)
                    ibc_loc[ic] = i + ipts[m1];
                    jb[ic] = book2.jbpt[i];
                    kb[ic] = book2.kbpt[i];
                    lb[ic] = book2.lbpt[i];
                }
            }
        }
    }

    // hole boundary (fringe point)
    int nh = conec1.nhole[m];
    if (nh != 0) {
        for (int m1 = 1; m1 <= conec1.nmesh; ++m1) list[m1] = 0;

        for (int n = 1; n <= nh; ++n) {
            int nserch = conec1.MHOLE(m, n);
            for (int nn = 1; nn <= nserch; ++nn) {
                int m1 = conec1.LHOLE(m, n, nn);
                if (list[m1] == 0) list[m1] = 1;
            }
        }
        for (int m1 = 1; m1 <= conec1.nmesh; ++m1) {
            if (list[m1] > 0) {
                int iset = book1.MHBS(m, m1);
                int is = book1.IPNTR(m1, iset);
                int ie = book1.NPNTR(m1, iset);
                getint(m1, book2.ji, book2.ki, book2.li,
                       book2.jbpt, book2.kbpt, book2.lbpt,
                       book2.dxint, book2.dyint, book2.dzint);
                for (int i = is; i <= ie; ++i) {
                    ic = ic + 1;
                    ibc_loc[ic] = i + ipts[m1];
                    jb[ic] = book2.jbpt[i];
                    kb[ic] = book2.kbpt[i];
                    lb[ic] = book2.lbpt[i];
                }
            }
        }
    }

    // set pointers (out): ibpnts=ic, iipnts=ibpts(m)
    ibpnts = ic;
    iipnts = book2.ibpts[m];

    // check dimensions
    if (ibpnts > idim) {
        // 601 format(/,/,' ',10x,'failure in cindex the number of boundary',
        //  1x,'points in mesh ',i3,1x,'is ibpnts = ',i5,
        //  1x,'and the dimension is idim = ',i5)
        std::fprintf(funit(6),
            "\n\n           failure in cindex the number of boundary points in mesh %3d is ibpnts = %5d and the dimension is idim = %5d\n",
            mesh, ibpnts, idim);
        std::fprintf(stderr, "cindex\n");
        std::exit(1);   // stop 'cindex'
    }
}

// ===========================================================================
//  setptr — set interpolation pointers and load lists.
// ===========================================================================
void setptr(int m, int m1, int icase) {
    std::fprintf(funit(6), "\n");   // write(6,*) empty -> blank line
    std::fprintf(funit(6), "     updating pointers for mesh %12d\n", m);

    // update pointers
    getint(m1, book2.ji, book2.ki, book2.li,
           book2.jbpt, book2.kbpt, book2.lbpt,
           book2.dxint, book2.dyint, book2.dzint);

    book1.nsets[m1] = book1.nsets[m1] + 1;
    int n = book1.nsets[m1];
    if (n > mhldim * mdim) {
        // 601 format('0',10x,'failure in setptr. the number of sets, n = ',
        //  i5,1x,'exceeds the max dimension = ',i5,1x,'for mesh ',i3)
        std::fprintf(funit(6),
            "0          failure in setptr. the number of sets, n = %5d exceeds the max dimension = %5d for mesh %3d\n",
            n, mhldim * mdim, m1);
        std::fprintf(stderr, "setptr \n");
        std::exit(1);
    }

    book1.IPNTR(m1, n) = book2.ibpts[m1] + 1;
    book1.NPNTR(m1, n) = book1.IPNTR(m1, n) + intrp2.itotal - 1;
    book2.ibpts[m1]    = book2.ibpts[m1] + intrp2.itotal;

    std::fprintf(funit(6), "       nsets(m1),m1= %12d%12d\n", n, m1);
    std::fprintf(funit(6), "       ibpts,ipntr(m1,n),npntr(m1,n)= %12d%12d%12d\n",
                 book2.ibpts[m1], book1.IPNTR(m1, n), book1.NPNTR(m1, n));

    if (book2.ibpts[m1] > idim) {
        // 602 format('0',10x,'failure in setup. the number of points,',
        //  ' ibpts =',i5,1x,'exceeds the max dimension = ',i5,1x,'for mesh ',i3)
        std::fprintf(funit(6),
            "0          failure in setup. the number of points, ibpts =%5d exceeds the max dimension = %5d for mesh %3d\n",
            book2.ibpts[m1], idim, m1);
        std::fprintf(stderr, "setptr\n");
        std::exit(1);
    }

    if (icase == 1) {
        // hole boundaries
        book1.MHBS(m, m1) = n;
        std::fprintf(funit(6), "       m,m1,mhbs(m,m1)= %12d%12d%12d\n", m, m1, book1.MHBS(m, m1));
    } else if (icase == 2) {
        // outer boundaries
        book1.MOBS(m, m1) = n;
        std::fprintf(funit(6), "       m,m1,mobs(m,m1)= %12d%12d%12d\n", m, m1, book1.MOBS(m, m1));
    } else {
        std::fprintf(stderr, "setptr\n");
        std::exit(1);
    }

    // update lists
    int ib = book1.IPNTR(m1, n) - 1;
    for (int i = 1; i <= intrp2.itotal; ++i) {
        ib = ib + 1;
        book2.jbpt[ib]  = intrp2.jb[i];
        book2.kbpt[ib]  = intrp2.kb[i];
        book2.lbpt[ib]  = intrp2.lb[i];
        book2.ji[ib]    = intrp2.jn[i];
        book2.ki[ib]    = intrp2.kn[i];
        book2.li[ib]    = intrp2.ln[i];
        book2.dxint[ib] = intrp1.xi[i];
        book2.dyint[ib] = intrp1.yi[i];
        book2.dzint[ib] = intrp1.zi[i];
    }
    putint(ib, m1, book2.ji, book2.ki, book2.li,
           book2.jbpt, book2.kbpt, book2.lbpt,
           book2.dxint, book2.dyint, book2.dzint);
}

// ===========================================================================
//  pltpts_ — output PLOT3D file (unit 9) for visual check of point classes.
//   real*4 x,y,z etc. stay float -> written as 4-byte floats to unit 9.
// ===========================================================================
void pltpts_(int nmesh, int m, int jmax, int kmax, int lmax,
             double* xx, double* yy, double* zz, int* ib, int* ibb) {
    // real*4 local coordinate copies (jdim,kdim,ldim); huge -> heap-allocate.
    // 1-based flat, size (jdim*kdim*ldim)+1.
    const long nxyz = (long)jdim * kdim * ldim;
    std::vector<float> xbuf(nxyz + 1), ybuf(nxyz + 1), zbuf(nxyz + 1);
    float* x = xbuf.data();
    float* y = ybuf.data();
    float* z = zbuf.data();
    // real*4 xbnd(idim) etc.
    std::vector<float> xbnd(idim + 1), ybnd(idim + 1), zbnd(idim + 1);
    std::vector<float> xorph(idim + 1), yorph(idim + 1), zorph(idim + 1);
    std::vector<float> xill(idim + 1), yill(idim + 1), zill(idim + 1);

    View3<double> XX(xx, jmax, kmax, lmax), YY(yy, jmax, kmax, lmax), ZZ(zz, jmax, kmax, lmax);
    View3<int>    IB(ib, jmax, kmax, lmax), IBB(ibb, jmax, kmax, lmax);
    // local x/y/z are dimensioned (jdim,kdim,ldim)
    View3<float>  X(x, jdim, kdim, ldim), Y(y, jdim, kdim, ldim), Z(z, jdim, kdim, ldim);

    for (int j = 1; j <= jmax; ++j)
        for (int k = 1; k <= kmax; ++k)
            for (int l = 1; l <= lmax; ++l) {
                X(j, k, l) = (float)XX(j, k, l);
                Y(j, k, l) = (float)YY(j, k, l);
                Z(j, k, l) = (float)ZZ(j, k, l);
            }

    // (iflag1/2/3 are set but unused in the Fortran)

    if (m == 1) {
        int nmesh3 = 3 * nmesh;
        int i = 0;
        for (int imesh = 1; imesh <= nmesh3; imesh += 3) {
            i = i + 1;
            diagno.jp[imesh] = grid2.mjmax[i];
            diagno.kp[imesh] = grid2.mkmax[i];
            diagno.lp[imesh] = grid2.mlmax[i];
        }
        for (int imesh = 1; imesh <= nmesh3; imesh += 3) {
            diagno.jp[imesh + 1] = diagno.jp[imesh];
            diagno.jp[imesh + 2] = diagno.jp[imesh];
            diagno.kp[imesh + 1] = diagno.kp[imesh];
            diagno.kp[imesh + 2] = diagno.kp[imesh];
            diagno.lp[imesh + 1] = diagno.lp[imesh];
            diagno.lp[imesh + 2] = diagno.lp[imesh];
        }

        if (pltpt.nbnd > 0) {
            nmesh3 = nmesh3 + 1;
            diagno.jp[nmesh3] = pltpt.nbnd;
            diagno.kp[nmesh3] = 1;
            diagno.lp[nmesh3] = 1;
        }
        if (pltpt.norph > 0) {
            nmesh3 = nmesh3 + 1;
            diagno.jp[nmesh3] = pltpt.norph;
            diagno.kp[nmesh3] = 1;
            diagno.lp[nmesh3] = 1;
        }
        if (pltpt.nill > 0) {
            nmesh3 = nmesh3 + 1;
            diagno.jp[nmesh3] = pltpt.nill;
            diagno.kp[nmesh3] = 1;
            diagno.lp[nmesh3] = 1;
        }

        // write(9) nmesh3
        uwrite(9, &nmesh3, sizeof(int));
        // write(9) (jp(i),kp(i),lp(i),i=1,nmesh3)   -> ONE record
        {
            RecordWriter rw;
            for (int ii = 1; ii <= nmesh3; ++ii) {
                rw.put_i(diagno.jp[ii]);
                rw.put_i(diagno.kp[ii]);
                rw.put_i(diagno.lp[ii]);
            }
            rw.flush(funit(9));
        }
    }

    // (ibl0/ibl1/ibln are counters, set but unused in output)

    // ----- output FIELD points -----
    for (int l = 1; l <= lmax; ++l)
        for (int k = 1; k <= kmax; ++k)
            for (int j = 1; j <= jmax; ++j) {
                int iblk = IBB(j, k, l);
                if (iblk == 1)       IB(j, k, l) = 1;
                else if (iblk == 0)  IB(j, k, l) = 0;
                else                 IB(j, k, l) = 0;
            }

    // helper to emit a PLOT3D solution record: three coord planes + ib, ONE record.
    // For ialph==0: x,y,z ; for ialph!=0: x,z,-y. Values are real*4, ib is int.
    auto write_block = [&](void) {
        RecordWriter rw;
        if (igrdtyp.ialph == 0) {
            for (int l = 1; l <= lmax; ++l) for (int k = 1; k <= kmax; ++k) for (int j = 1; j <= jmax; ++j) rw.put_f(X(j, k, l));
            for (int l = 1; l <= lmax; ++l) for (int k = 1; k <= kmax; ++k) for (int j = 1; j <= jmax; ++j) rw.put_f(Y(j, k, l));
            for (int l = 1; l <= lmax; ++l) for (int k = 1; k <= kmax; ++k) for (int j = 1; j <= jmax; ++j) rw.put_f(Z(j, k, l));
            for (int l = 1; l <= lmax; ++l) for (int k = 1; k <= kmax; ++k) for (int j = 1; j <= jmax; ++j) rw.put_i(IB(j, k, l));
        } else {
            for (int l = 1; l <= lmax; ++l) for (int k = 1; k <= kmax; ++k) for (int j = 1; j <= jmax; ++j) rw.put_f(X(j, k, l));
            for (int l = 1; l <= lmax; ++l) for (int k = 1; k <= kmax; ++k) for (int j = 1; j <= jmax; ++j) rw.put_f(Z(j, k, l));
            for (int l = 1; l <= lmax; ++l) for (int k = 1; k <= kmax; ++k) for (int j = 1; j <= jmax; ++j) rw.put_f(-Y(j, k, l));
            for (int l = 1; l <= lmax; ++l) for (int k = 1; k <= kmax; ++k) for (int j = 1; j <= jmax; ++j) rw.put_i(IB(j, k, l));
        }
        rw.flush(funit(9));
    };

    write_block();

    // ----- output HOLE points -----
    for (int l = 1; l <= lmax; ++l)
        for (int k = 1; k <= kmax; ++k)
            for (int j = 1; j <= jmax; ++j) {
                int iblk = IBB(j, k, l);
                if (iblk == 1)       IB(j, k, l) = 0;
                else if (iblk == 0)  IB(j, k, l) = 1;
                else                 IB(j, k, l) = 0;
            }
    write_block();

    // ----- output FRINGE points -----
    for (int l = 1; l <= lmax; ++l)
        for (int k = 1; k <= kmax; ++k)
            for (int j = 1; j <= jmax; ++j) {
                int iblk = IBB(j, k, l);
                if (iblk == 0 || iblk == 1) IB(j, k, l) = 0;
                else                        IB(j, k, l) = 1;
            }
    write_block();

    if (m == nmesh) {
        View3<int> IB1(ib, jmax, kmax, lmax); // ib(i,1,1) addressing == ib[i]

        // boundary points
        if (pltpt.nbnd > 0) {
            for (int i = 1; i <= pltpt.nbnd; ++i) {
                xbnd[i] = (float)pltpt.xbnd[i];
                ybnd[i] = (float)pltpt.ybnd[i];
                zbnd[i] = (float)pltpt.zbnd[i];
                IB1(i, 1, 1) = 1;
            }
            RecordWriter rw;
            if (igrdtyp.ialph == 0) {
                for (int i = 1; i <= pltpt.nbnd; ++i) rw.put_f(xbnd[i]);
                for (int i = 1; i <= pltpt.nbnd; ++i) rw.put_f(ybnd[i]);
                for (int i = 1; i <= pltpt.nbnd; ++i) rw.put_f(zbnd[i]);
                for (int i = 1; i <= pltpt.nbnd; ++i) rw.put_i(IB1(i, 1, 1));
            } else {
                for (int i = 1; i <= pltpt.nbnd; ++i) rw.put_f(xbnd[i]);
                for (int i = 1; i <= pltpt.nbnd; ++i) rw.put_f(zbnd[i]);
                for (int i = 1; i <= pltpt.nbnd; ++i) rw.put_f(-ybnd[i]);
                for (int i = 1; i <= pltpt.nbnd; ++i) rw.put_i(IB1(i, 1, 1));
            }
            rw.flush(funit(9));
        }

        // orphan (extrapolated) points
        if (pltpt.norph > 0) {
            for (int i = 1; i <= pltpt.norph; ++i) {
                xorph[i] = (float)pltpt.xorph[i];
                yorph[i] = (float)pltpt.yorph[i];
                zorph[i] = (float)pltpt.zorph[i];
                IB1(i, 1, 1) = 1;
            }
            RecordWriter rw;
            if (igrdtyp.ialph == 0) {
                for (int i = 1; i <= pltpt.norph; ++i) rw.put_f(xorph[i]);
                for (int i = 1; i <= pltpt.norph; ++i) rw.put_f(yorph[i]);
                for (int i = 1; i <= pltpt.norph; ++i) rw.put_f(zorph[i]);
                for (int i = 1; i <= pltpt.norph; ++i) rw.put_i(IB1(i, 1, 1));
            } else {
                for (int i = 1; i <= pltpt.norph; ++i) rw.put_f(xorph[i]);
                for (int i = 1; i <= pltpt.norph; ++i) rw.put_f(zorph[i]);
                for (int i = 1; i <= pltpt.norph; ++i) rw.put_f(-yorph[i]);
                for (int i = 1; i <= pltpt.norph; ++i) rw.put_i(IB1(i, 1, 1));
            }
            rw.flush(funit(9));
        }

        // illegal-stencil points
        if (pltpt.nill > 0) {
            for (int i = 1; i <= pltpt.nill; ++i) {
                xill[i] = (float)pltpt.xill[i];
                yill[i] = (float)pltpt.yill[i];
                zill[i] = (float)pltpt.zill[i];
                IB1(i, 1, 1) = 1;
            }
            RecordWriter rw;
            if (igrdtyp.ialph == 0) {
                for (int i = 1; i <= pltpt.nill; ++i) rw.put_f(xill[i]);
                for (int i = 1; i <= pltpt.nill; ++i) rw.put_f(yill[i]);
                for (int i = 1; i <= pltpt.nill; ++i) rw.put_f(zill[i]);
                for (int i = 1; i <= pltpt.nill; ++i) rw.put_i(IB1(i, 1, 1));
            } else {
                for (int i = 1; i <= pltpt.nill; ++i) rw.put_f(xill[i]);
                for (int i = 1; i <= pltpt.nill; ++i) rw.put_f(zill[i]);
                for (int i = 1; i <= pltpt.nill; ++i) rw.put_f(-yill[i]);
                for (int i = 1; i <= pltpt.nill; ++i) rw.put_i(IB1(i, 1, 1));
            }
            rw.flush(funit(9));
        }
    }
}

// ===========================================================================
//  trace_ — write search-routine history for the current point to unit 7.
//   Formats replicated byte-for-byte. e11.4 via tc::fmtE(v,11,4).
//   List-directed write(7,*) spacing matched (leading blank + one blank
//   between items; integers/reals in Fortran default list-directed widths).
// ===========================================================================
void trace_(int icall, int idum1, int idum2, int idum3, int idum4,
            double dum1, double dum2, double dum3) {
    const int iunit = 7;

    if (trace1.itrace < 0) return;

    if (icall == -1) {
        // 99 format(' ',2x,'*** Search History For Fringe Points ***',
        //          /,11x,'*** Searching Mesh ',i2,' ***')
        // ' ',2x = 3 leading spaces; continuation line 11x = 11 spaces.
        std::fprintf(funit(iunit),
            "   *** Search History For Fringe Points ***\n"
            "           *** Searching Mesh %2d ***\n", idum1);
    }

    if (icall == 0) {
        // three blank list-directed writes then format 100
        std::fprintf(funit(iunit), "\n");
        std::fprintf(funit(iunit), "\n");
        std::fprintf(funit(iunit), "\n");
        // 100 format(' ',2x,'*** Search History For Boundary Points ***',
        //           /,12x,'*** Searching Mesh ',i2,' ***')
        std::fprintf(funit(iunit),
            "   *** Search History For Boundary Points ***\n"
            "            *** Searching Mesh %2d ***\n", idum1);
    }

    if (icall == 1) {
        if (trace1.itrace == 0) {
            frewind(iunit);
            // write(iunit,100) with EMPTY i/o list: gfortran writes the format up
            // to the i2 data descriptor, then stops (no i2 output, no trailing
            // ' ***'). Result: first line + '            *** Searching Mesh '
            // (trailing space from the literal preceding i2), no ' ***'.
            std::fprintf(funit(iunit),
                "   *** Search History For Boundary Points ***\n"
                "            *** Searching Mesh \n");
        }
    }

    if (icall == 1) {
        // 101 format(' ',2x,'iterations for interpolation to point  i =',i6,
        //   /,2x,'  with coordinates xp,yp,zp=',e11.4,',',e11.4,',',e11.4)
        std::fprintf(funit(iunit),
            "   iterations for interpolation to point  i =%6d\n", idum1);
        std::fprintf(funit(iunit),
            "    with coordinates xp,yp,zp=%s,%s,%s\n",
            tc::fmtE(dum1, 11, 4).c_str(),
            tc::fmtE(dum2, 11, 4).c_str(),
            tc::fmtE(dum3, 11, 4).c_str());
    }

    if (icall == 11) {
        // 1101 format(' ',2x,'iterations for extrapolation to point  i =',i6,
        //   /,2x,'  with coordinates xp,yp,zp=',e11.4,',',e11.4,',',e11.4)
        std::fprintf(funit(iunit),
            "   iterations for extrapolation to point  i =%6d\n", idum1);
        std::fprintf(funit(iunit),
            "    with coordinates xp,yp,zp=%s,%s,%s\n",
            tc::fmtE(dum1, 11, 4).c_str(),
            tc::fmtE(dum2, 11, 4).c_str(),
            tc::fmtE(dum3, 11, 4).c_str());
    }

    if (icall == 3) {
        // 103 format(' ',4x,'intern=',i2)
        std::fprintf(funit(iunit), "     intern=%2d\n", idum1);
    }

    if (icall == 4) {
        // 104 format(' ',4x,'searching in cell j,k,i = ',i3,i3,i3,' of block ',i3)
        std::fprintf(funit(iunit),
            "     searching in cell j,k,i = %3d%3d%3d of block %3d\n",
            idum1, idum2, idum3, idum4);
    }

    if (icall == 5) {
        // 105 format(' ',8x,'xie,eta,zeta=',e10.3,',',e10.3,',',e10.3,' (local values)')
        std::fprintf(funit(iunit),
            "         xie,eta,zeta=%s,%s,%s (local values)\n",
            tc::fmtE(dum1, 10, 3).c_str(),
            tc::fmtE(dum2, 10, 3).c_str(),
            tc::fmtE(dum3, 10, 3).c_str());
    }

    if (icall == 6) {
        std::fprintf(funit(iunit),
            " search off track...using dsmin to get back on track\n");
    }

    if (icall == 7) {
        std::fprintf(funit(iunit),
            " frozen convergence...using dsmin to attempt to break cycle\n");
    }

    if (icall == 17) {
        std::fprintf(funit(iunit),
            " frozen convergence...will extrapolate from best cell\n");
    }

    if (icall == 20) {
        // 120 format(' ',4x,'giving up...will try in another mesh or',/,
        //           6x,'attempt extrapolation if that fails')
        std::fprintf(funit(iunit),
            "     giving up...will try in another mesh or\n"
            "      attempt extrapolation if that fails\n");
    }

    if (icall == 21) {
        // write(7,*) '      miscue in xe2: newton iteration not',
        //  ' converged but 0 < xie,eta,zeta < 1'
        std::fprintf(funit(iunit),
            "       miscue in xe2: newton iteration not converged but 0 < xie,eta,zeta < 1\n");
        // write(7,*) '      xie,eta,zeta = ',dum1,dum2,dum3  (list-directed reals)
        std::fprintf(funit(iunit),
            "       xie,eta,zeta = %s%s%s\n",
            ld_r8(dum1).c_str(), ld_r8(dum2).c_str(), ld_r8(dum3).c_str());
    }

    if (icall == 22) {
        // write(7,*) '      iteration = ',idum1,' error = ',dum1,
        //  ' error tolerance = ',dum2
        // list-directed: after the first real (followed by a char constant),
        // gfortran inserts one extra separator blank -> two spaces before 'error'.
        std::fprintf(funit(iunit),
            "       iteration = %12d  error = %s  error tolerance = %s\n",
            idum1, ld_r8(dum1).c_str(), ld_r8(dum2).c_str());
        std::fprintf(funit(iunit),
            "       will set xie,eta,zeta to large value totrigger min. distance search\n");
    }

    if (icall == 23) {
        // write(7,*)'   mesh ',idum1,' branch cut: jpc,kpc,lpc = ',idum2,idum3,idum4
        std::fprintf(funit(iunit),
            "    mesh %12d  branch cut: jpc,kpc,lpc = %12d%12d%12d\n",
            idum1, idum2, idum3, idum4);
    }

    if (icall == 24) {
        // write(7,*)'                        jp,kp,lp = ',idum2,idum3,idum4
        std::fprintf(funit(iunit),
            "                         jp,kp,lp = %12d%12d%12d\n",
            idum2, idum3, idum4);
    }
}
