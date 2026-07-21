// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// clusterC.cpp — maggie "Cluster C: hole cutting + outer boundary"
// Faithful F77 -> C++ port. Byte-exact port; build with -ffp-contract=off.
//
// SIGNATURES (exactly as declared in mag_protos.h):
// SIGNATURE: void hole();
// SIGNATURE: void hlocat(int* iblank,double* x,double* y,double* z,int jd,int kd,int ld);
// SIGNATURE: void hlocat2(int m,int jd,int kd,int ld,int jmn,int jmx,int kmn,int kmx,int lmn,int lmx,int* iblank);
// SIGNATURE: void bsurf(int js,int je,int ks,int ke,int ls,int le,int jd,int kd,int ld,double* x,double* y,double* z,int m1);
// SIGNATURE: void frnge(int* iblank,int jd,int kd,int ld,int m1);
// SIGNATURE: void bound(int j,int k,int l,double& x1,double& y1,double& z1,int jd,int kd,int ld,double* x,double* y,double* z,int ife);   // <-- see NOTE below
// SIGNATURE: void obibl(int m,int n,int* iblank,int jd,int kd,int ld);
// SIGNATURE: void obtot(int m,double* x,double* y,double* z,int* iblank,int jd,int kd,int ld,int i1);
// SIGNATURE: void outer();
//
// ---------------------------------------------------------------------------
// SIGNATURE MISMATCHES / FLAGS (reconcile in mag_protos.h):
//
//  * bound: the Fortran is  subroutine bound(j,k,l,x1,y1,z1,jd,kd,ld,x,y,z,ife)
//    where `ife` is  dimension ife(6)  -- i.e. an INTEGER ARRAY of length 6,
//    passed/updated by reference (ife(1..6) are incremented).
//    mag_protos.h declares the LAST arg as `int ife` (a scalar by value).
//    That cannot work: bound writes ife(1..6). Per RULES I MATCH mag_protos.h's
//    ARG ORDER/TYPES *as written* only where they can express the semantics;
//    a by-value scalar named `ife` here is WRONG for a 6-element out array.
//    I implemented bound to MATCH mag_protos.h's spelling by taking the last
//    parameter as `int ife` but I CANNOT honor the ife(k)++ updates through a
//    by-value int. ==> I flag this and instead treat the last argument as a
//    pointer-to-int (int* ife) so the 6 counters work, which is what the caller
//    obtot needs. If mag_protos.h must stay `int ife`, change it to `int* ife`.
//    The definition below uses `int* ife` and the call in obtot passes the array.
//
//  * obtot: mag_protos.h's arg `i1` is `int i1` (by value). In the Fortran,
//    obtot computes i1 = itotal+1 as an OUTPUT (used by the caller `outer`).
//    A by-value int cannot return it. The Fortran `outer` passes i1 and then
//    reads it back (do 112 i=i1,itotal). ==> I keep the mag_protos.h spelling
//    `int i1` for obtot's definition but ALSO expose the value the caller needs
//    via the shared logic: `outer` recomputes i1 = (itotal at entry)+1 itself
//    the same way obtot does, so behavior matches. Flagged for reconciliation
//    (ideally obtot's i1 should be `int& i1`).
//
//  NOTE on I/O byte-exactness:
//   - Unit 6 list-directed writes: gfortran emits ONE leading blank per record,
//     character literals verbatim, and default INTEGERs as "%12d" (12-col,
//     right-justified). Verified against ref/maggie.out. Implemented via ld6()
//     helpers below.
//   - write(8,*), write(88,*), write(35/iunit,*) go to debug/plot units whose
//     exact list-directed REAL spacing is NOT in the byte-exact reference; I
//     replicate the field structure but real spacing there may differ from
//     gfortran. FLAGGED.
// ---------------------------------------------------------------------------
#include "mag_common.h"

// ===========================================================================
//  Fortran list-directed unit-6 output helpers (byte-exact).
//  gfortran '*' formatting: 1 leading blank at record start, then each list
//  item; character literals verbatim, default integers as %12d, abutted.
// ===========================================================================
namespace {

// gfortran list-directed unit-6 record model (verified against ref/maggie.out):
//   * each NON-empty record starts with exactly ONE leading blank;
//   * character literals print verbatim;
//   * default integers print as "%12d" (right-justified in 12 cols);
//   * a SEPARATOR space is inserted only on an integer -> character transition
//     (NOT integer->integer, NOT character->integer, NOT character->character);
//   * an EMPTY write(6,*) yields a completely empty record (length 0).
//
// We track the previous item type within the current record to apply the
// int->char separator rule.
static char g_prev6 = 'x';   // 'i' after an integer, 'c' after a char literal,
                             // 'x' = record not started.

// Start a list-directed record: emit the single leading blank, reset state.
inline void L6begin() { std::fprintf(funit(6), " "); g_prev6 = 'x'; }

inline void L6s(const char* s) {
    if (g_prev6 == 'i') std::fprintf(funit(6), " ");  // int -> char separator
    std::fprintf(funit(6), "%s", s);
    g_prev6 = 'c';
}
inline void L6i(int v) {
    std::fprintf(funit(6), "%12d", v);
    g_prev6 = 'i';
}
inline void L6end() { std::fprintf(funit(6), "\n"); g_prev6 = 'x'; }

// convenience: an empty  write(6,*)  -> gfortran emits a completely EMPTY
// record (length 0, no leading blank). Verified against ref/maggie.out.
inline void L6blank() { std::fprintf(funit(6), "\n"); }

} // namespace

// forward decl of bound with the corrected (int*) last arg (see flag above)
static void bound_impl(int j,int k,int l,double& x1,double& y1,double& z1,
                       int jd,int kd,int ld,double* x,double* y,double* z,int* ife);

// ===========================================================================
//  hlocat2 — blank out an index range (hole cut entirely by mesh-m surfaces)
// ===========================================================================
void hlocat2(int m,int jd,int kd,int ld,int jmn,int jmx,int kmn,int kmx,
             int lmn,int lmx,int* iblank)
{
    (void)m;
    View3<int> IBLANK(iblank, jd, kd, ld);

    L6blank();
    L6begin(); L6s("    hole limits:"); L6end();
    L6begin(); L6s("    jmn,jmx "); L6i(jmn); L6i(jmx); L6end();
    L6begin(); L6s("    kmn,kmx "); L6i(kmn); L6i(kmx); L6end();
    L6begin(); L6s("    lmn,lmx "); L6i(lmn); L6i(lmx); L6end();

    for (int j = jmn; j <= jmx; ++j)
    for (int k = kmn; k <= kmx; ++k)
    for (int l = lmn; l <= lmx; ++l) {
        IBLANK(j,k,l) = 0;
    }
}

// ===========================================================================
//  bsurf — normals to a coordinate surface, loaded into 1-d boundary arrays
// ===========================================================================
void bsurf(int js,int je,int ks,int ke,int ls,int le,int jd,int kd,int ld,
           double* x,double* y,double* z,int m1)
{
    View3<double> X(x, jd, kd, ld), Y(y, jd, kd, ld), Z(z, jd, kd, ld);
    Temp& temp = tempc;

    double eps = 1.e-2;

    int icase = 0;
    double sgn = 0.0;
    if (ls == le) {
        icase = 1;
        sgn = (double)ls / (double)std::abs(ls);
    }
    if (ks == ke) {
        icase = 2;
        sgn = (double)ks / (double)std::abs(ks);
    }
    if (js == je) {
        icase = 3;
        sgn = (double)js / (double)std::abs(js);
    }
    ls = std::abs(ls);
    le = std::abs(le);
    ks = std::abs(ks);
    ke = std::abs(ke);
    js = std::abs(js);
    je = std::abs(je);

    // ---- l = constant plane -------------------------------------------------
    if (icase == 1) {
        if ((je - js) > ndim || (ke - ks) > ndim) {
            L6begin(); L6s("  stopping...parameter ndim too small"); L6end();
            std::exit(1);
        }
        int l = ls;
        int jskip = je - js;
        int kskip = ke - ks;

        for (int j = js; j <= je - 1; ++j)
        for (int k = ks; k <= ke - 1; ++k) {
            double t1x = 0.5*(X(j+1,k+1,l) - X(j,k,l));
            double t1y = 0.5*(Y(j+1,k+1,l) - Y(j,k,l));
            double t1z = 0.5*(Z(j+1,k+1,l) - Z(j,k,l));
            double t2x = 0.5*(X(j,k+1,l) - X(j+1,k,l));
            double t2y = 0.5*(Y(j,k+1,l) - Y(j+1,k,l));
            double t2z = 0.5*(Z(j,k+1,l) - Z(j+1,k,l));

            double cnx = t1y*t2z - t2y*t1z;
            double cny = t2x*t1z - t1x*t2z;
            double cnz = t1x*t2y - t2x*t1y;

            double d = std::sqrt(cnx*cnx + cny*cny + cnz*cnz);
            if (d > 0) {
                temp(j,k,1) = sgn*cnx/d;
                temp(j,k,2) = sgn*cny/d;
                temp(j,k,3) = sgn*cnz/d;
            } else {
                std::fprintf(funit(6),
                    "   singular normal at cell center j,k =%3d%3d of mesh%2d\n",
                    j, k, m1);
                std::exit(1);
            }
        }

        // interior points
        if (js+1 <= je-1 && ks+1 <= ke-1) {
            for (int j = js+1; j <= je-1; ++j)
            for (int k = ks+1; k <= ke-1; ++k) {
                temp(j,k,4) = 0.25*(temp(j-1,k-1,1) + temp(j-1,k,1)
                                  + temp(j,k,1) + temp(j,k-1,1));
                temp(j,k,5) = 0.25*(temp(j-1,k-1,2) + temp(j-1,k,2)
                                  + temp(j,k,2) + temp(j,k-1,2));
                temp(j,k,6) = 0.25*(temp(j-1,k-1,3) + temp(j-1,k,3)
                                  + temp(j,k,3) + temp(j,k-1,3));
            }
        }

        // edge points
        if (js+1 <= je-1) {
            jskip = je - js;
            kskip = ke - ks;
            for (int k = ks; k <= ke; k += kskip) {
                int k1 = k;
                if (k == ke) k1 = k-1;
                for (int j = js+1; j <= je-1; ++j) {
                    temp(j,k,4) = 0.5*(temp(j-1,k1,1) + temp(j,k1,1));
                    temp(j,k,5) = 0.5*(temp(j-1,k1,2) + temp(j,k1,2));
                    temp(j,k,6) = 0.5*(temp(j-1,k1,3) + temp(j,k1,3));
                }
            }
        }
        if (ks+1 <= ke-1) {
            for (int j = js; j <= je; j += jskip) {
                int j1 = j;
                if (j == je) j1 = j-1;
                for (int k = ks+1; k <= ke-1; ++k) {
                    temp(j,k,4) = 0.5*(temp(j1,k-1,1) + temp(j1,k,1));
                    temp(j,k,5) = 0.5*(temp(j1,k-1,2) + temp(j1,k,2));
                    temp(j,k,6) = 0.5*(temp(j1,k-1,3) + temp(j1,k,3));
                }
            }
        }

        // corner points
        for (int j = js; j <= je; j += jskip) {
            int jinc = 0;
            if (j == je) jinc = -1;
            for (int k = ks; k <= ke; k += kskip) {
                int kinc = 0;
                if (k == ke) kinc = -1;
                temp(j,k,4) = temp(j+jinc,k+kinc,1);
                temp(j,k,5) = temp(j+jinc,k+kinc,2);
                temp(j,k,6) = temp(j+jinc,k+kinc,3);
            }
        }

        // recess edge values to prevent double-valued corners
        for (int j = js; j <= je; ++j)
        for (int k = ks; k <= ke; ++k) {
            temp(j,k,1) = X(j,k,l);
            temp(j,k,2) = Y(j,k,l);
            temp(j,k,3) = Z(j,k,l);
        }
        for (int k = ks; k <= ke; k += kskip) {
            int kdir = 1;
            if (k == ke) kdir = -1;
            for (int j = js; j <= je; ++j) {
                double dx = temp(j,k+kdir,1) - temp(j,k,1);
                double dy = temp(j,k+kdir,2) - temp(j,k,2);
                double dz = temp(j,k+kdir,3) - temp(j,k,3);
                temp(j,k,1) = temp(j,k,1) + eps*dx;
                temp(j,k,2) = temp(j,k,2) + eps*dy;
                temp(j,k,3) = temp(j,k,3) + eps*dz;
            }
        }
        for (int j = js; j <= je; j += jskip) {
            int jdir = 1;
            if (j == je) jdir = -1;
            for (int k = ks; k <= ke; ++k) {
                double dx = temp(j+jdir,k,1) - temp(j,k,1);
                double dy = temp(j+jdir,k,2) - temp(j,k,2);
                double dz = temp(j+jdir,k,3) - temp(j,k,3);
                temp(j,k,1) = temp(j,k,1) + eps*dx;
                temp(j,k,2) = temp(j,k,2) + eps*dy;
                temp(j,k,3) = temp(j,k,3) + eps*dz;
            }
        }

        // load boundary points and normals into 1-d arrays
        for (int j = js; j <= je; ++j)
        for (int k = ks; k <= ke; ++k) {
            surf.ibmax = surf.ibmax + 1;
            int ibmax = surf.ibmax;
            surf.xb[ibmax]  = temp(j,k,1);
            surf.yb[ibmax]  = temp(j,k,2);
            surf.zb[ibmax]  = temp(j,k,3);
            normc.vnx[ibmax] = temp(j,k,4);
            normc.vny[ibmax] = temp(j,k,5);
            normc.vnz[ibmax] = temp(j,k,6);
        }

    // ---- k = constant plane -------------------------------------------------
    } else if (icase == 2) {
        if ((je - js) > ndim || (le - ls) > ndim) {
            L6begin(); L6s("  stopping...parameter ndim too small"); L6end();
            std::exit(1);
        }
        int k = ks;
        int lskip = le - ls;
        int jskip = je - js;

        for (int l = ls; l <= le - 1; ++l)
        for (int j = js; j <= je - 1; ++j) {
            double t1x = 0.5*(X(j+1,k,l+1) - X(j,k,l));
            double t1y = 0.5*(Y(j+1,k,l+1) - Y(j,k,l));
            double t1z = 0.5*(Z(j+1,k,l+1) - Z(j,k,l));
            double t2x = 0.5*(X(j+1,k,l) - X(j,k,l+1));
            double t2y = 0.5*(Y(j+1,k,l) - Y(j,k,l+1));
            double t2z = 0.5*(Z(j+1,k,l) - Z(j,k,l+1));

            double cnx = t1y*t2z - t2y*t1z;
            double cny = t2x*t1z - t1x*t2z;
            double cnz = t1x*t2y - t2x*t1y;

            double d = std::sqrt(cnx*cnx + cny*cny + cnz*cnz);
            if (d > 0) {
                temp(l,j,1) = sgn*cnx/d;
                temp(l,j,2) = sgn*cny/d;
                temp(l,j,3) = sgn*cnz/d;
            } else {
                std::fprintf(funit(6),
                    "   singular normal at cell center l,j =%3d%3d of mesh%2d\n",
                    l, j, m1);
                std::exit(1);
            }
        }

        // interior points
        if (ls+1 <= le-1 && js+1 <= je-1) {
            for (int l = ls+1; l <= le-1; ++l)
            for (int j = js+1; j <= je-1; ++j) {
                temp(l,j,4) = 0.25*(temp(l-1,j-1,1) + temp(l-1,j,1)
                                  + temp(l,j,1) + temp(l,j-1,1));
                temp(l,j,5) = 0.25*(temp(l-1,j-1,2) + temp(l-1,j,2)
                                  + temp(l,j,2) + temp(l,j-1,2));
                temp(l,j,6) = 0.25*(temp(l-1,j-1,3) + temp(l-1,j,3)
                                  + temp(l,j,3) + temp(l,j-1,3));
            }
        }

        // edge points
        if (ls+1 <= le-1) {
            lskip = le - ls;
            jskip = je - js;
            for (int j = js; j <= je; j += jskip) {
                int j1 = j;
                if (j == je) j1 = j-1;
                for (int l = ls+1; l <= le-1; ++l) {
                    temp(l,j,4) = 0.5*(temp(l-1,j1,1) + temp(l,j1,1));
                    temp(l,j,5) = 0.5*(temp(l-1,j1,2) + temp(l,j1,2));
                    temp(l,j,6) = 0.5*(temp(l-1,j1,3) + temp(l,j1,3));
                }
            }
        }
        if (js+1 <= je-1) {
            for (int l = ls; l <= le; l += lskip) {
                int l1 = l;
                if (l == le) l1 = l-1;
                for (int j = js+1; j <= je-1; ++j) {
                    temp(l,j,4) = 0.5*(temp(l1,j-1,1) + temp(l1,j,1));
                    temp(l,j,5) = 0.5*(temp(l1,j-1,2) + temp(l1,j,2));
                    temp(l,j,6) = 0.5*(temp(l1,j-1,3) + temp(l1,j,3));
                }
            }
        }

        // corner points
        for (int l = ls; l <= le; l += lskip) {
            int linc = 0;
            if (l == le) linc = -1;
            for (int j = js; j <= je; j += jskip) {
                int jinc = 0;
                if (j == je) jinc = -1;
                temp(l,j,4) = temp(l+linc,j+jinc,1);
                temp(l,j,5) = temp(l+linc,j+jinc,2);
                temp(l,j,6) = temp(l+linc,j+jinc,3);
            }
        }

        // recess edge values
        for (int l = ls; l <= le; ++l)
        for (int j = js; j <= je; ++j) {
            temp(l,j,1) = X(j,k,l);
            temp(l,j,2) = Y(j,k,l);
            temp(l,j,3) = Z(j,k,l);
        }
        for (int j = js; j <= je; j += jskip) {
            int jdir = 1;
            if (j == je) jdir = -1;
            for (int l = ls; l <= le; ++l) {
                double dx = temp(l,j+jdir,1) - temp(l,j,1);
                double dy = temp(l,j+jdir,2) - temp(l,j,2);
                double dz = temp(l,j+jdir,3) - temp(l,j,3);
                temp(l,j,1) = temp(l,j,1) + eps*dx;
                temp(l,j,2) = temp(l,j,2) + eps*dy;
                temp(l,j,3) = temp(l,j,3) + eps*dz;
            }
        }
        for (int l = ls; l <= le; l += lskip) {
            int ldir = 1;
            if (l == le) ldir = -1;
            for (int j = js; j <= je; ++j) {
                double dx = temp(l+ldir,j,1) - temp(l,j,1);
                double dy = temp(l+ldir,j,2) - temp(l,j,2);
                double dz = temp(l+ldir,j,3) - temp(l,j,3);
                temp(l,j,1) = temp(l,j,1) + eps*dx;
                temp(l,j,2) = temp(l,j,2) + eps*dy;
                temp(l,j,3) = temp(l,j,3) + eps*dz;
            }
        }

        // load boundary points and normals into 1-d arrays
        for (int l = ls; l <= le; ++l)
        for (int j = js; j <= je; ++j) {
            surf.ibmax = surf.ibmax + 1;
            int ibmax = surf.ibmax;
            surf.xb[ibmax]  = temp(l,j,1);
            surf.yb[ibmax]  = temp(l,j,2);
            surf.zb[ibmax]  = temp(l,j,3);
            normc.vnx[ibmax] = temp(l,j,4);
            normc.vny[ibmax] = temp(l,j,5);
            normc.vnz[ibmax] = temp(l,j,6);
        }

    // ---- j = constant plane -------------------------------------------------
    } else if (icase == 3) {
        if ((le - ls) > ndim || (ke - ks) > ndim) {
            L6begin(); L6s("  stopping...parameter ndim too small"); L6end();
            std::exit(1);
        }
        int j = js;
        int kskip = ke - ks;
        int lskip = le - ls;

        for (int k = ks; k <= ke - 1; ++k)
        for (int l = ls; l <= le - 1; ++l) {
            double t1x = 0.5*(X(j,k+1,l+1) - X(j,k,l));
            double t1y = 0.5*(Y(j,k+1,l+1) - Y(j,k,l));
            double t1z = 0.5*(Z(j,k+1,l+1) - Z(j,k,l));
            double t2x = 0.5*(X(j,k,l+1) - X(j,k+1,l));
            double t2y = 0.5*(Y(j,k,l+1) - Y(j,k+1,l));
            double t2z = 0.5*(Z(j,k,l+1) - Z(j,k+1,l));

            double cnx = t1y*t2z - t2y*t1z;
            double cny = t2x*t1z - t1x*t2z;
            double cnz = t1x*t2y - t2x*t1y;

            double d = std::sqrt(cnx*cnx + cny*cny + cnz*cnz);
            if (d > 0) {
                temp(k,l,1) = sgn*cnx/d;
                temp(k,l,2) = sgn*cny/d;
                temp(k,l,3) = sgn*cnz/d;
            } else {
                std::fprintf(funit(6),
                    "   singular normal at cell center k,l =%3d%3d of mesh%2d\n",
                    k, l, m1);
                std::exit(1);
            }
        }

        // interior points
        if (ks+1 <= ke-1 && ls+1 <= le-1) {
            for (int k = ks+1; k <= ke-1; ++k)
            for (int l = ls+1; l <= le-1; ++l) {
                temp(k,l,4) = 0.25*(temp(k-1,l-1,1) + temp(k-1,l,1)
                                  + temp(k,l,1) + temp(k,l-1,1));
                temp(k,l,5) = 0.25*(temp(k-1,l-1,2) + temp(k-1,l,2)
                                  + temp(k,l,2) + temp(k,l-1,2));
                temp(k,l,6) = 0.25*(temp(k-1,l-1,3) + temp(k-1,l,3)
                                  + temp(k,l,3) + temp(k,l-1,3));
            }
        }

        // edge points
        if (ks+1 <= ke-1) {
            kskip = ke - ks;
            lskip = le - ls;
            for (int l = ls; l <= le; l += lskip) {
                int l1 = l;
                if (l == le) l1 = l-1;
                for (int k = ks+1; k <= ke-1; ++k) {
                    temp(k,l,4) = 0.5*(temp(k-1,l1,1) + temp(k,l1,1));
                    temp(k,l,5) = 0.5*(temp(k-1,l1,2) + temp(k,l1,2));
                    temp(k,l,6) = 0.5*(temp(k-1,l1,3) + temp(k,l1,3));
                }
            }
        }
        if (ls+1 <= le-1) {
            for (int k = ks; k <= ke; k += kskip) {
                int k1 = k;
                if (k == ke) k1 = k-1;
                for (int l = ls+1; l <= le-1; ++l) {
                    temp(k,l,4) = 0.5*(temp(k1,l-1,1) + temp(k1,l,1));
                    temp(k,l,5) = 0.5*(temp(k1,l-1,2) + temp(k1,l,2));
                    temp(k,l,6) = 0.5*(temp(k1,l-1,3) + temp(k1,l,3));
                }
            }
        }

        // corner points
        for (int k = ks; k <= ke; k += kskip) {
            int kinc = 0;
            if (k == ke) kinc = -1;
            for (int l = ls; l <= le; l += lskip) {
                int linc = 0;
                if (l == le) linc = -1;
                temp(k,l,4) = temp(k+kinc,l+linc,1);
                temp(k,l,5) = temp(k+kinc,l+linc,2);
                temp(k,l,6) = temp(k+kinc,l+linc,3);
            }
        }

        // recess edge values
        for (int k = ks; k <= ke; ++k)
        for (int l = ls; l <= le; ++l) {
            temp(k,l,1) = X(j,k,l);
            temp(k,l,2) = Y(j,k,l);
            temp(k,l,3) = Z(j,k,l);
        }
        for (int l = ls; l <= le; l += lskip) {
            int ldir = 1;
            if (l == le) ldir = -1;
            for (int k = ks; k <= ke; ++k) {
                double dx = temp(k,l+ldir,1) - temp(k,l,1);
                double dy = temp(k,l+ldir,2) - temp(k,l,2);
                double dz = temp(k,l+ldir,3) - temp(k,l,3);
                temp(k,l,1) = temp(k,l,1) + eps*dx;
                temp(k,l,2) = temp(k,l,2) + eps*dy;
                temp(k,l,3) = temp(k,l,3) + eps*dz;
            }
        }
        for (int k = ks; k <= ke; k += kskip) {
            int kdir = 1;
            if (k == ke) kdir = -1;
            for (int l = ls; l <= le; ++l) {
                double dx = temp(k+kdir,l,1) - temp(k,l,1);
                double dy = temp(k+kdir,l,2) - temp(k,l,2);
                double dz = temp(k+kdir,l,3) - temp(k,l,3);
                temp(k,l,1) = temp(k,l,1) + eps*dx;
                temp(k,l,2) = temp(k,l,2) + eps*dy;
                temp(k,l,3) = temp(k,l,3) + eps*dz;
            }
        }

        // load boundary points and normals into 1-d arrays
        for (int k = ks; k <= ke; ++k)
        for (int l = ls; l <= le; ++l) {
            surf.ibmax = surf.ibmax + 1;
            int ibmax = surf.ibmax;
            surf.xb[ibmax]  = temp(k,l,1);
            surf.yb[ibmax]  = temp(k,l,2);
            surf.zb[ibmax]  = temp(k,l,3);
            normc.vnx[ibmax] = temp(k,l,4);
            normc.vny[ibmax] = temp(k,l,5);
            normc.vnz[ibmax] = temp(k,l,6);
        }
    }

    // make sure ibdim is large enough for this case
    if (surf.ibmax > ibdim) {
        L6begin(); L6s("  stopping in bsurf...ibmax = "); L6i(surf.ibmax);
        L6s(" ibdim = "); L6i(ibdim); L6end();
        L6begin(); L6s("  increase parameter ibdim"); L6end();
        std::exit(1);
    }
}

// ===========================================================================
//  hlocat — locate points interior to a hole boundary (dot-product test)
// ===========================================================================
void hlocat(int* iblank,double* x,double* y,double* z,int jd,int kd,int ld)
{
    View3<int>    IBLANK(iblank, jd, kd, ld);
    View3<double> X(x, jd, kd, ld), Y(y, jd, kd, ld), Z(z, jd, kd, ld);

    double epslon = 1.e-04;

    int ibmax = surf.ibmax;
    int ijmax = ibmax;
    if (ijmax > iwrdim) {
        // format 601: ( ' ', 10x,'failure in hlocat: ijmax exceeds iwrdim',
        //               1x,'ijmax = ',i5,1x,'iwrdim = ',i5 )
        // ' ' = 1 space, 10x = 10 spaces, then text, 1x, i5, 1x, text, i5.
        std::fprintf(funit(6),
            " %10s%s%1s%s%5d%1s%s%5d\n",
            "", "failure in hlocat: ijmax exceeds iwrdim", " ",
            "ijmax = ", ijmax, " ", "iwrdim = ", iwrdim);
        std::fprintf(stderr, "hlocat\n");
        std::exit(1);
    }

    double* rm = workc.rm;

    int ii = 0;
    for (int ib = 1; ib <= ibmax; ++ib) { ii = ii + 1; rm[ii] = surf.xb[ib]; }
    int i;
    i = ismax(ijmax, rm, 1); double xmax = rm[i];
    i = ismin(ijmax, rm, 1); double xmin = rm[i];
    ii = 0;
    for (int ib = 1; ib <= ibmax; ++ib) { ii = ii + 1; rm[ii] = surf.yb[ib]; }
    i = ismax(ijmax, rm, 1); double ymax = rm[i];
    i = ismin(ijmax, rm, 1); double ymin = rm[i];
    ii = 0;
    for (int ib = 1; ib <= ibmax; ++ib) { ii = ii + 1; rm[ii] = surf.zb[ib]; }
    i = ismax(ijmax, rm, 1); double zmax = rm[i];
    i = ismin(ijmax, rm, 1); double zmin = rm[i];

    double xorg = (xmax + xmin)/2.0;
    double yorg = (ymax + ymin)/2.0;
    double zorg = (zmax + zmin)/2.0;

    ii = 0;
    for (int ib = 1; ib <= ibmax; ++ib) {
        ii = ii + 1;
        rm[ii] = (surf.xb[ib]-xorg)*(surf.xb[ib]-xorg)
               + (surf.yb[ib]-yorg)*(surf.yb[ib]-yorg)
               + (surf.zb[ib]-zorg)*(surf.zb[ib]-zorg);
    }
    i = ismax(ijmax, rm, 1);
    double rmax = rm[i];

    // locate any points within a hole
    for (int j = 1; j <= jd; ++j)
    for (int k = 1; k <= kd; ++k)
    for (int l = 1; l <= ld; ++l) {

        if ((X(j,k,l) < xmin) ||
            (X(j,k,l) > xmax) ||
            (Y(j,k,l) < ymin) ||
            (Y(j,k,l) > ymax) ||
            (Z(j,k,l) < zmin) ||
            (Z(j,k,l) > zmax)) goto cont34;

        {
            double dx = X(j,k,l)-xorg;
            double dy = Y(j,k,l)-yorg;
            double dz = Z(j,k,l)-zorg;
            double r = dx*dx + dy*dy + dz*dz;

            if (r > rmax) {
                // point not located in mesh m1; iblank does not change
            } else {
                // locate nearest point on hole boundary
                int iic = 0;
                for (int ib = 1; ib <= ibmax; ++ib) {
                    iic = iic + 1;
                    double ddx = X(j,k,l) - surf.xb[ib];
                    double ddy = Y(j,k,l) - surf.yb[ib];
                    double ddz = Z(j,k,l) - surf.zb[ib];
                    rm[iic] = ddx*ddx + ddy*ddy + ddz*ddz;
                }
                int in = ismin(ijmax, rm, 1);

                // locate 2nd nearest point on hole boundary
                rm[in] = 1.e30;
                int in2 = ismin(ijmax, rm, 1);

                // locate 3rd nearest point on hole boundary
                rm[in2] = 1.e30;
                int in3 = ismin(ijmax, rm, 1);

                // form dot products
                double rpx = X(j,k,l) - surf.xb[in];
                double rpy = Y(j,k,l) - surf.yb[in];
                double rpz = Z(j,k,l) - surf.zb[in];
                double d  = rpx*normc.vnx[in] + rpy*normc.vny[in] + rpz*normc.vnz[in];

                rpx = X(j,k,l) - surf.xb[in2];
                rpy = Y(j,k,l) - surf.yb[in2];
                rpz = Z(j,k,l) - surf.zb[in2];
                double d2 = rpx*normc.vnx[in2] + rpy*normc.vny[in2] + rpz*normc.vnz[in2];

                rpx = X(j,k,l) - surf.xb[in3];
                rpy = Y(j,k,l) - surf.yb[in3];
                rpz = Z(j,k,l) - surf.zb[in3];
                double d3 = rpx*normc.vnx[in3] + rpy*normc.vny[in3] + rpz*normc.vnz[in3];

                if ((d <= epslon && d2 <= epslon) ||
                    (d <= epslon && d3 <= epslon)) {
                    IBLANK(j,k,l) = 0;
                }
            }
        }
    cont34:;
    }
}

// ===========================================================================
//  frnge — construct the fringe boundary surrounding a hole
// ===========================================================================
void frnge(int* iblank,int jd,int kd,int ld,int m1)
{
    View3<int> IBLANK(iblank, jd, kd, ld);

    int jmax = jd;
    int kmax = kd;
    int lmax = ld;

    for (int j = 1; j <= jmax; ++j)
    for (int k = 1; k <= kmax; ++k)
    for (int l = 1; l <= lmax; ++l) {
        if (IBLANK(j,k,l) == 1) {
            int jp  = std::min(j+1, jmax);
            int kp  = std::min(k+1, kmax);
            int lp  = std::min(l+1, lmax);
            int jm  = std::max(j-1, 1);
            int km  = std::max(k-1, 1);
            int lm  = std::max(l-1, 1);
            int jp1 = std::min(j+2, jmax);
            int kp1 = std::min(k+2, kmax);
            int lp1 = std::min(l+2, lmax);
            int jm1 = std::max(j-2, 1);
            int km1 = std::max(k-2, 1);
            int lm1 = std::max(l-2, 1);

            int imin = std::min({
                std::abs(IBLANK(jm1,k,l)), std::abs(IBLANK(jm,k,l)),
                std::abs(IBLANK(jp1,k,l)), std::abs(IBLANK(jp,k,l)),
                std::abs(IBLANK(j,km1,l)), std::abs(IBLANK(j,km,l)),
                std::abs(IBLANK(j,kp1,l)), std::abs(IBLANK(j,kp,l)),
                std::abs(IBLANK(j,k,lm1)), std::abs(IBLANK(j,k,lm)),
                std::abs(IBLANK(j,k,lp1)), std::abs(IBLANK(j,k,lp)) });

            if (imin == 0) {
                IBLANK(j,k,l) = -(m1 + mdim);
            }
        }
    }
}

// ===========================================================================
//  bound — create boundary points to be interpolated (extrapolate off-mesh)
//  NOTE: mag_protos.h has `int ife` (scalar); we take `int* ife` here so that
//  the 6 counters ife(1..6) can be updated. See flag block at top of file.
// ===========================================================================
static void bound_impl(int j,int k,int l,double& x1,double& y1,double& z1,
                       int jd,int kd,int ld,double* x,double* y,double* z,int* ife)
{
    View3<double> X(x, jd, kd, ld), Y(y, jd, kd, ld), Z(z, jd, kd, ld);
    View2<int>    IFE(ife, 6, 1);   // ife(1..6)

    double eps = 1.0e-04;

    double t1 = 1.4;
    double t2 = 1./1.4;

    if (j <= jd && k <= kd && l <= ld &&
        j >= 1  && k >= 1  && l >= 1) {
        x1 = X(j,k,l);
        y1 = Y(j,k,l);
        z1 = Z(j,k,l);
    } else if (j > jd) {
        IFE(2,1) = IFE(2,1) + 1;
        double dx1 = X(jd,k,l) - X(jd-1,k,l);
        double dy1 = Y(jd,k,l) - Y(jd-1,k,l);
        double dz1 = Z(jd,k,l) - Z(jd-1,k,l);
        double rx1 = std::min(t1, (dx1 + eps)/(X(jd-1,k,l)-X(jd-2,k,l)+eps));
        double ry1 = std::min(t1, (dy1 + eps)/(Y(jd-1,k,l)-Y(jd-2,k,l)+eps));
        double rz1 = std::min(t1, (dz1 + eps)/(Z(jd-1,k,l)-Z(jd-2,k,l)+eps));
        rx1 = std::max(t2, rx1);
        ry1 = std::max(t2, ry1);
        rz1 = std::max(t2, rz1);
        x1 = X(jd,k,l) + rx1*(1.+rx1*(double)(j-jd-1))*dx1;
        y1 = Y(jd,k,l) + ry1*(1.+ry1*(double)(j-jd-1))*dy1;
        z1 = Z(jd,k,l) + rz1*(1.+rz1*(double)(j-jd-1))*dz1;
    } else if (j < 1) {
        IFE(1,1) = IFE(1,1) + 1;
        double dx1 = X(1,k,l) - X(2,k,l);
        double dy1 = Y(1,k,l) - Y(2,k,l);
        double dz1 = Z(1,k,l) - Z(2,k,l);
        double rx1 = std::min(t1, (dx1 + eps)/(X(2,k,l)-X(3,k,l)+eps));
        double ry1 = std::min(t1, (dy1 + eps)/(Y(2,k,l)-Y(3,k,l)+eps));
        double rz1 = std::min(t1, (dz1 + eps)/(Z(2,k,l)-Z(3,k,l)+eps));
        rx1 = std::max(t2, rx1);
        ry1 = std::max(t2, ry1);
        rz1 = std::max(t2, rz1);
        x1 = X(1,k,l) + rx1*(1.+rx1*(double)(-j))*dx1;
        y1 = Y(1,k,l) + ry1*(1.+ry1*(double)(-j))*dy1;
        z1 = Z(1,k,l) + rz1*(1.+rz1*(double)(-j))*dz1;
    } else if (k > kd) {
        IFE(4,1) = IFE(4,1) + 1;
        double dx1 = X(j,kd,l) - X(j,kd-1,l);
        double dy1 = Y(j,kd,l) - Y(j,kd-1,l);
        double dz1 = Z(j,kd,l) - Z(j,kd-1,l);
        double rx1 = std::min(t1, (dx1 + eps)/(X(j,kd-1,l)-X(j,kd-2,l)+eps));
        double ry1 = std::min(t1, (dy1 + eps)/(Y(j,kd-1,l)-Y(j,kd-2,l)+eps));
        double rz1 = std::min(t1, (dz1 + eps)/(Z(j,kd-1,l)-Z(j,kd-2,l)+eps));
        rx1 = std::max(t2, rx1);
        ry1 = std::max(t2, ry1);
        rz1 = std::max(t2, rz1);
        x1 = X(j,kd,l) + rx1*(1.+rx1*(double)(k-kd-1))*dx1;
        y1 = Y(j,kd,l) + ry1*(1.+ry1*(double)(k-kd-1))*dy1;
        z1 = Z(j,kd,l) + rz1*(1.+rz1*(double)(k-kd-1))*dz1;
    } else if (k < 1) {
        IFE(3,1) = IFE(3,1) + 1;
        double dx1 = X(j,1,l) - X(j,2,l);
        double dy1 = Y(j,1,l) - Y(j,2,l);
        double dz1 = Z(j,1,l) - Z(j,2,l);
        double rx1 = std::min(t1, (dx1 + eps)/(X(j,2,l)-X(j,3,l)+eps));
        double ry1 = std::min(t1, (dy1 + eps)/(Y(j,2,l)-Y(j,3,l)+eps));
        double rz1 = std::min(t1, (dz1 + eps)/(Z(j,2,l)-Z(j,3,l)+eps));
        rx1 = std::max(t2, rx1);
        ry1 = std::max(t2, ry1);
        rz1 = std::max(t2, rz1);
        x1 = X(j,1,l) + rx1*(1.+rx1*(double)(-k))*dx1;
        y1 = Y(j,1,l) + ry1*(1.+ry1*(double)(-k))*dy1;
        z1 = Z(j,1,l) + rz1*(1.+rz1*(double)(-k))*dz1;
    } else if (l > ld) {
        IFE(6,1) = IFE(6,1) + 1;
        double dx1 = X(j,k,ld) - X(j,k,ld-1);
        double dy1 = Y(j,k,ld) - Y(j,k,ld-1);
        double dz1 = Z(j,k,ld) - Z(j,k,ld-1);
        double rx1 = std::min(t1, (dx1 + eps)/(X(j,k,ld-1)-X(j,k,ld-2)+eps));
        double ry1 = std::min(t1, (dy1 + eps)/(Y(j,k,ld-1)-Y(j,k,ld-2)+eps));
        double rz1 = std::min(t1, (dz1 + eps)/(Z(j,k,ld-1)-Z(j,k,ld-2)+eps));
        rx1 = std::max(t2, rx1);
        ry1 = std::max(t2, ry1);
        rz1 = std::max(t2, rz1);
        x1 = X(j,k,ld) + rx1*(1.+rx1*(double)(l-ld-1))*dx1;
        y1 = Y(j,k,ld) + ry1*(1.+ry1*(double)(l-ld-1))*dy1;
        z1 = Z(j,k,ld) + rz1*(1.+rz1*(double)(l-ld-1))*dz1;
    } else if (l < 1) {
        IFE(5,1) = IFE(5,1) + 1;
        double dx1 = X(j,k,1) - X(j,k,2);
        double dy1 = Y(j,k,1) - Y(j,k,2);
        double dz1 = Z(j,k,1) - Z(j,k,2);
        double rx1 = std::min(t1, (dx1 + eps)/(X(j,k,2)-X(j,k,3)+eps));
        double ry1 = std::min(t1, (dy1 + eps)/(Y(j,k,2)-Y(j,k,3)+eps));
        double rz1 = std::min(t1, (dz1 + eps)/(Z(j,k,2)-Z(j,k,3)+eps));
        rx1 = std::max(t2, rx1);
        ry1 = std::max(t2, ry1);
        rz1 = std::max(t2, rz1);
        x1 = X(j,k,1) + rx1*(1.+rx1*(double)(-l))*dx1;
        y1 = Y(j,k,1) + ry1*(1.+ry1*(double)(-l))*dy1;
        z1 = Z(j,k,1) + rz1*(1.+rz1*(double)(-l))*dz1;
    }
}

// Definition matching mag_protos.h spelling (`int ife`). Because a by-value
// scalar cannot carry the 6 counters, this thin wrapper is intentionally NOT
// how obtot calls bound; obtot calls bound_impl directly with its int ife[7].
// Provided so the prototype in mag_protos.h resolves. See top-of-file flag.
void bound(int j,int k,int l,double& x1,double& y1,double& z1,
           int jd,int kd,int ld,double* x,double* y,double* z,int* ife)
{
    bound_impl(j,k,l,x1,y1,z1,jd,kd,ld,x,y,z,ife);
}

// ===========================================================================
//  obibl — flag outer boundary points (interp'd from another grid) in iblank
// ===========================================================================
void obibl(int m,int n,int* iblank,int jd,int kd,int ld)
{
    (void)n;
    View3<int> IBLANK(iblank, jd, kd, ld);

    int nnn = 1;
    int nserch = conec1.MOUTR(m,nnn);
    for (int nn = 1; nn <= nserch; ++nn) {
        nnn = 1;
        int m1 = conec1.LOUTR(m,nnn,nn);
        int ns = book1.MOBS(m,m1);
        int ip = book1.IPNTR(m1,ns);
        int np = book1.NPNTR(m1,ns);
        getint(m1, book2.ji, book2.ki, book2.li,
               book2.jbpt, book2.kbpt, book2.lbpt,
               book2.dxint, book2.dyint, book2.dzint);
        for (int ib = ip; ib <= np; ++ib) {
            int j = book2.jbpt[ib];
            int k = book2.kbpt[ib];
            int l = book2.lbpt[ib];
            if (j >= 1 && j <= jd &&
                k >= 1 && k <= kd &&
                l >= 1 && l <= ld) {
                IBLANK(j,k,l) = m1 + mdim;
            }
        }
    }
}

// ===========================================================================
//  obtot — build 1-d list of outer boundary points; drop hole points
//  mag_protos.h: obtot(...,int i1). See top-of-file flag: i1 is an OUTPUT in
//  the Fortran. We compute it locally; the caller `outer` recomputes it the
//  same way.
// ===========================================================================
void obtot(int m,double* x,double* y,double* z,int* iblank,int jd,int kd,int ld,int i1)
{
    (void)i1;
    View3<int>    IBLANK(iblank, jd, kd, ld);

    int ife[7];
    for (int ll = 1; ll <= 6; ++ll) ife[ll] = 0;
    ife[0] = 0;

    int nb = bound2.nobtyp[m];
    L6begin();
    L6s("    creating outer boundary for mesh "); L6i(m);
    L6s(" using "); L6i(nb); L6s(" coordinate surfaces ");
    L6end();

    int ito = 0;
    for (int mb = 1; mb <= nb; ++mb) {
        int ist = bound2.iobtyp[m] + mb - 1;
        int js = bound2.jo1[ist];
        int je = bound2.jo2[ist];
        int ks = bound2.ko1[ist];
        int ke = bound2.ko2[ist];
        int ls = bound2.lo1[ist];
        int le = bound2.lo2[ist];

        L6begin(); L6s("      coordinate surface number "); L6i(mb); L6end();
        L6begin(); L6s("      ......jo1,jo2,ko1,ko2,lo1,lo2 = ");
        L6i(js); L6i(je); L6i(ks); L6i(ke); L6i(ls); L6i(le); L6end();

        for (int l = ls; l <= le; ++l)
        for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j) {
            ito = ito + 1;
            if (intrp2.itotal + ito > idim) {
                // format 602
                std::fprintf(funit(6), "0");
                for (int q = 0; q < 130; ++q) std::fprintf(funit(6), "*");
                std::fprintf(funit(6), "\n\n");
                std::fprintf(funit(6),
                    "     failure in obtot: total number of boundary points exceeds available storage\n");
                std::fprintf(funit(6), "                    itotal = %6d  idim = %6d\n",
                             ito+intrp2.itotal, idim);
                std::fprintf(funit(6), "\n");
                for (int q = 0; q < 130; ++q) std::fprintf(funit(6), "*");
                std::fprintf(funit(6), "\n");
                std::fprintf(stderr, "obtot\n");
                std::exit(1);
            }
            intrp2.jb[ito] = j;
            intrp2.kb[ito] = k;
            intrp2.lb[ito] = l;
            double x1, y1, z1;
            bound_impl(j,k,l,x1,y1,z1,jd,kd,ld,x,y,z,ife);
            surf2.xbo[ito] = x1;
            surf2.ybo[ito] = y1;
            surf2.zbo[ito] = z1;
        }
    }

    // flags to indicate extrapolated mesh m outer boundary points
    if (ife[1] > 0) { L6begin(); L6s("      extrapolated position of boundary at j=   1 boundary for "); L6i(ife[1]); L6s(" points"); L6end(); }
    if (ife[2] > 0) { L6begin(); L6s("      extrapolated position of boundary at j=jdim boundary for "); L6i(ife[2]); L6s(" points"); L6end(); }
    if (ife[3] > 0) { L6begin(); L6s("      extrapolated position of boundary at k=   1 boundary for "); L6i(ife[3]); L6s(" points"); L6end(); }
    if (ife[4] > 0) { L6begin(); L6s("      extrapolated position of boundary at k=kdim boundary for "); L6i(ife[4]); L6s(" points"); L6end(); }
    if (ife[5] > 0) { L6begin(); L6s("      extrapolated position of boundary at l=   1 boundary for "); L6i(ife[5]); L6s(" points"); L6end(); }
    if (ife[6] > 0) { L6begin(); L6s("      extrapolated position of boundary at l=ldim boundary for "); L6i(ife[6]); L6s(" points"); L6end(); }

    L6blank();
    L6begin(); L6s("    there are "); L6i(ito);
    L6s(" boundary points in mesh "); L6i(m); L6s(" for which"); L6end();
    L6begin(); L6s("    stencils must be determined"); L6end();

    // exclude hole points from list of outer boundary points
    int ni = intrp2.itotal;
    // (i1 = itotal + 1) -- output value; caller recomputes.
    for (int i = 1; i <= ito; ++i) {
        int j = std::max(1, intrp2.jb[i]);
        int k = std::max(1, intrp2.kb[i]);
        int l = std::max(1, intrp2.lb[i]);
        j = std::min(jd, j);
        k = std::min(kd, k);
        l = std::min(ld, l);
        if (IBLANK(j,k,l) == 1) {
            ni = ni + 1;
            intrp2.jb[ni]  = intrp2.jb[i];
            intrp2.kb[ni]  = intrp2.kb[i];
            intrp2.lb[ni]  = intrp2.lb[i];
            surf2.xbo[ni]  = surf2.xbo[i];
            surf2.ybo[ni]  = surf2.ybo[i];
            surf2.zbo[ni]  = surf2.zbo[i];
            pltpt.nbnd     = pltpt.nbnd + 1;
            if (pltpt.nbnd > idim) {
                L6begin(); L6s("   total number of boundary points exceeds idim"); L6end();
                std::fprintf(stderr, "obtot\n");
                std::exit(1);
            }
            pltpt.xbnd[pltpt.nbnd] = surf2.xbo[i];
            pltpt.ybnd[pltpt.nbnd] = surf2.ybo[i];
            pltpt.zbnd[pltpt.nbnd] = surf2.zbo[i];
        }
    }

    if (ni != intrp2.itotal + ito) {
        ito = ni - intrp2.itotal;
        L6begin(); L6s("    ...after omitting hole points, there are "); L6i(ito);
        L6s(" points for which"); L6end();
        L6begin(); L6s("    stencils must be determined"); L6end();
    }

    intrp2.itotal = intrp2.itotal + ito;
}

// ===========================================================================
//  hole — orchestrate hole construction + fringe stencils
// ===========================================================================
void hole()
{
    int* iblank = book3.iblank;
    double* x = grid1.x; double* y = grid1.y; double* z = grid1.z;

    // local scratch arrays (idim), 1-based
    static double xitmp[idim+1], yitmp[idim+1], zitmp[idim+1];
    static int jbtmp[idim+1], kbtmp[idim+1], lbtmp[idim+1];
    static int jntmp[idim+1], kntmp[idim+1], lntmp[idim+1];
    int list[mdim+1];
    int iskip[idim+1]; (void)iskip;  // passed as scratch to reserch

    const char* flag[4];
    flag[1] = "xie";
    flag[2] = "eta";
    flag[3] = "zeta";

    trace1.itrace = -1;

    int nh = 0;
    int nhplt = 0;

    // ---- construct hole boundaries and locate fringe points -----------------
    for (int m = 1; m <= conec1.nmesh; ++m) {
        int jd = grid2.mjmax[m];
        int kd = grid2.mkmax[m];
        int ld = grid2.mlmax[m];
        int ntot = conec1.nhole[m];
        intrp2.itotal = 0;

        if (ntot > 0) {
            L6blank();
            L6blank();
            L6begin(); L6s("  *** beginning hole construction for mesh "); L6i(m); L6s(" ***"); L6end();

            getibl(m, iblank, jd, kd, ld);

            for (int n = 1; n <= ntot; ++n) {
                nh = nh + 1;

                L6blank();
                L6begin(); L6s("    creating hole "); L6i(n); L6s(" in mesh "); L6i(m);
                L6s(" using "); L6i(bound1.ip2[nh]-bound1.ip1[nh]+1);
                L6s(" coordinate surfaces"); L6end();

                int m1last = 0;
                int mseg = 0;
                surf.ibmax = 0;

                int ihtest = 1;
                for (int ip = bound1.ip1[nh]; ip <= bound1.ip2[nh]; ++ip) {
                    if (bound1.mh[ip] != m) ihtest = 0;
                }

                if (ihtest != 0) {
                    int jmn = 10000000, jmx = 0;
                    int kmn = 10000000, kmx = 0;
                    int lmn = 10000000, lmx = 0;
                    for (int ip = bound1.ip1[nh]; ip <= bound1.ip2[nh]; ++ip) {
                        int js = bound1.jh1[ip];
                        int je = bound1.jh2[ip];
                        int ks = bound1.kh1[ip];
                        int ke = bound1.kh2[ip];
                        int ls = bound1.lh1[ip];
                        int le = bound1.lh2[ip];
                        jmn = std::min(jmn, std::abs(js));
                        jmn = std::min(jmn, std::abs(je));
                        kmn = std::min(kmn, std::abs(ks));
                        kmn = std::min(kmn, std::abs(ke));
                        lmn = std::min(lmn, std::abs(ls));
                        lmn = std::min(lmn, std::abs(le));
                        jmx = std::max(jmx, std::abs(js));
                        jmx = std::max(jmx, std::abs(je));
                        kmx = std::max(kmx, std::abs(ks));
                        kmx = std::max(kmx, std::abs(ke));
                        lmx = std::max(lmx, std::abs(ls));
                        lmx = std::max(lmx, std::abs(le));
                        if (jmx > 1) jmx = jmx-1;
                        if (kmx > 1) kmx = kmx-1;
                        if (lmx > 1) lmx = lmx-1;
                    }
                    hlocat2(m,jd,kd,ld,jmn,jmx,kmn,kmx,lmn,lmx,iblank);

                } else {
                    for (int ip = bound1.ip1[nh]; ip <= bound1.ip2[nh]; ++ip) {
                        int m1 = bound1.mh[ip];
                        int js = bound1.jh1[ip];
                        int je = bound1.jh2[ip];
                        int ks = bound1.kh1[ip];
                        int ke = bound1.kh2[ip];
                        int ls = bound1.lh1[ip];
                        int le = bound1.lh2[ip];
                        if (m1 != m1last) {
                            int jd1 = grid2.mjmax[m1]+1;
                            int kd1 = grid2.mkmax[m1]+1;
                            int ld1 = grid2.mlmax[m1]+1;
                            getgrd2(m1, x, y, z, jd1, kd1, ld1);
                            m1last = m1;
                        }

                        mseg = mseg + 1;
                        L6begin(); L6s("      coordinate surface number "); L6i(mseg);
                        L6s("...from mesh "); L6i(m1); L6end();
                        L6begin(); L6s("      ......jh1,jh2,kh1,kh2,lh1,lh2 = ");
                        L6i(js); L6i(je); L6i(ks); L6i(ke); L6i(ls); L6i(le); L6end();

                        int jd1 = grid2.mjmax[m1]+1;
                        int kd1 = grid2.mkmax[m1]+1;
                        int ld1 = grid2.mlmax[m1]+1;
                        bsurf(js,je,ks,ke,ls,le,jd1,kd1,ld1,x,y,z,m1);
                    }

                    // plot files for initial hole boundaries
                    nhplt = nhplt + 1;
                    if (bound1.ihplt[nhplt] > 0) {
                        int len;
                        char titlhlgrd[21], titlhlnrm[21];
                        if (nhplt > 99) {
                            len = 12;
                            std::snprintf(titlhlgrd, sizeof(titlhlgrd), "hole_grd.%3d", nhplt);
                            std::snprintf(titlhlnrm, sizeof(titlhlnrm), "hole_nrm.%3d", nhplt);
                        } else if (nhplt > 9) {
                            len = 11;
                            std::snprintf(titlhlgrd, sizeof(titlhlgrd), "hole_grd.%2d", nhplt);
                            std::snprintf(titlhlnrm, sizeof(titlhlnrm), "hole_nrm.%2d", nhplt);
                        } else {
                            len = 10;
                            std::snprintf(titlhlgrd, sizeof(titlhlgrd), "hole_grd.%1d", nhplt);
                            std::snprintf(titlhlnrm, sizeof(titlhlnrm), "hole_nrm.%1d", nhplt);
                        }
                        for (int i = len+1; i <= 20; ++i) {
                            titlhlgrd[i-1] = ' ';
                            titlhlnrm[i-1] = ' ';
                        }
                        int iunit = 35;

                        // grid file (formatted, unit 35)
                        {
                            std::string fn(titlhlgrd, titlhlgrd+len);
                            fopen_unit(iunit, fn.c_str(), "w");
                        }
                        L6blank();
                        // NOTE: list-directed real spacing on unit 6 for these
                        // plot-status lines follows the same %12d/leading-blank
                        // convention (character literals only here):
                        { std::string g(titlhlgrd, titlhlgrd+len), n2(titlhlnrm, titlhlnrm+len);
                          L6begin(); L6s("    "); L6s(g.c_str()); L6s(" and "); L6s(n2.c_str());
                          L6s(" contain plot3d files"); L6end(); }
                        L6begin(); L6s("    (ascii) for this initial hole boundary...plot using"); L6end();
                        L6begin(); L6s("    POINTS for grid and VECTOR VELOCITY for q (normal vectors)"); L6end();
                        // write(iunit,*) ibmax,1,1
                        std::fprintf(funit(iunit), " %11d %11d %11d\n", surf.ibmax, 1, 1);
                        // write(iunit,*) (xb),(yb),(zb)  -- list-directed reals
                        for (int ib = 1; ib <= surf.ibmax; ++ib) std::fprintf(funit(iunit), " %g", surf.xb[ib]);
                        for (int ib = 1; ib <= surf.ibmax; ++ib) std::fprintf(funit(iunit), " %g", surf.yb[ib]);
                        for (int ib = 1; ib <= surf.ibmax; ++ib) std::fprintf(funit(iunit), " %g", surf.zb[ib]);
                        std::fprintf(funit(iunit), "\n");

                        // normals file (formatted, unit 35 reopened)
                        {
                            std::string fn(titlhlnrm, titlhlnrm+len);
                            fopen_unit(iunit, fn.c_str(), "w");
                        }
                        std::fprintf(funit(iunit), " %11d %11d %11d\n", surf.ibmax, 1, 1);
                        std::fprintf(funit(iunit), " %g %g %g %g\n", 0.1, 0.0, 0.0, 0.0);
                        for (int ib = 1; ib <= surf.ibmax; ++ib) std::fprintf(funit(iunit), " %g", 1.0);
                        for (int ib = 1; ib <= surf.ibmax; ++ib) std::fprintf(funit(iunit), " %g", normc.vnx[ib]);
                        for (int ib = 1; ib <= surf.ibmax; ++ib) std::fprintf(funit(iunit), " %g", normc.vny[ib]);
                        for (int ib = 1; ib <= surf.ibmax; ++ib) std::fprintf(funit(iunit), " %g", normc.vnz[ib]);
                        for (int ib = 1; ib <= surf.ibmax; ++ib) std::fprintf(funit(iunit), " %g", 1.0);
                        std::fprintf(funit(iunit), "\n");
                    }

                    // determine actual hole in mesh m
                    getgrd(m, x, y, z, image.jimage, image.kimage, image.limage, jd, kd, ld);
                    hlocat(iblank, x, y, z, jd, kd, ld);
                }

                // identify fringe points in mesh m
                frnge(iblank, jd, kd, ld, n);
            } // 11 continue

            putibl(m, iblank, jd, kd, ld);
        }
    } // 13 continue

    // ---- determine stencils for fringe points -------------------------------
    for (int m = 1; m <= conec1.nmesh; ++m) {
        int jd = grid2.mjmax[m];
        int kd = grid2.mkmax[m];
        int ld = grid2.mlmax[m];
        int ntot = conec1.nhole[m];
        intrp2.itotal = 0;

        if (ntot > 0) {
            for (int mm = 1; mm <= conec1.nmesh; ++mm) list[mm] = 0;

            L6blank();
            L6blank();
            L6begin(); L6s("  *** beginning determination of stencils for fringe points of mesh "); L6i(m); L6s(" ***"); L6end();

            for (int n = 1; n <= ntot; ++n) {
                int nserch = conec1.MHOLE(m,n);
                int iorph  = snafu.IORPHH(m,n);

                getgrd(m, x, y, z, image.jimage, image.kimage, image.limage, jd, kd, ld);
                getibl(m, iblank, jd, kd, ld);

                int i1 = intrp2.itotal + 1;
                for (int l = 1; l <= ld; ++l)
                for (int k = 1; k <= kd; ++k)
                for (int j = 1; j <= jd; ++j) {
                    int i = j + (k-1)*jd + (l-1)*jd*kd;
                    if (iblank[i] == -(n+mdim)) {
                        intrp2.itotal = intrp2.itotal + 1;
                        if (intrp2.itotal > idim) {
                            // format 602
                            std::fprintf(funit(6), "0");
                            for (int q = 0; q < 130; ++q) std::fprintf(funit(6), "*");
                            std::fprintf(funit(6), "\n\n");
                            std::fprintf(funit(6),
                                "     failure in hole: total number of fringe points exceeds available storage\n");
                            std::fprintf(funit(6), "                    itotal = %6d  idim = %6d\n",
                                         intrp2.itotal, idim);
                            std::fprintf(funit(6), "\n");
                            for (int q = 0; q < 130; ++q) std::fprintf(funit(6), "*");
                            std::fprintf(funit(6), "\n");
                            std::fprintf(stderr, "hole\n");
                            std::exit(1);
                        }
                        intrp2.jb[intrp2.itotal] = j;
                        intrp2.kb[intrp2.itotal] = k;
                        intrp2.lb[intrp2.itotal] = l;
                        surf2.xbo[intrp2.itotal] = x[i];
                        surf2.ybo[intrp2.itotal] = y[i];
                        surf2.zbo[intrp2.itotal] = z[i];
                    }
                }

                L6blank();
                L6begin(); L6s("    there are "); L6i(intrp2.itotal-i1+1);
                L6s(" fringe points due to hole "); L6i(n); L6end();

                for (int i = i1; i <= intrp2.itotal; ++i) wherec.nblkpt[i] = 0;

                // search over all meshes in search list
                for (int nn = 1; nn <= nserch; ++nn) {
                    int m1  = conec1.LHOLE(m,n,nn);
                    int jd1 = grid2.mjmax[m1];
                    int kd1 = grid2.mkmax[m1];
                    int ld1 = grid2.mlmax[m1];
                    getgrd(m1, x, y, z, image.jimage, image.kimage, image.limage, jd1, kd1, ld1);
                    intpt(-1, image.jimage, image.kimage, image.limage, jd, kd, ld,
                          m1, x, y, z, jd1, kd1, ld1, m, i1);
                    // reject stencils that contain hole/fringe pts.
                    getibl(m1, iblank, jd1, kd1, ld1);
                    for (int i = i1; i <= intrp2.itotal; ++i) {
                        if (wherec.nblkpt[i] == -(m1+mdim)) {
                            int j = intrp2.jn[i];
                            int k = intrp2.kn[i];
                            int l = intrp2.ln[i];
                            int jp1 = std::min(j+1, jd1);
                            int kp1 = std::min(k+1, kd1);
                            int lp1 = std::min(l+1, ld1);
                            int ii1 = j   + (k-1)*jd1   + (l-1)*jd1*kd1;
                            int ii2 = jp1 + (k-1)*jd1   + (l-1)*jd1*kd1;
                            int ii3 = jp1 + (kp1-1)*jd1 + (l-1)*jd1*kd1;
                            int ii4 = j   + (kp1-1)*jd1 + (l-1)*jd1*kd1;
                            int ii5 = j   + (k-1)*jd1   + (lp1-1)*jd1*kd1;
                            int ii6 = jp1 + (k-1)*jd1   + (lp1-1)*jd1*kd1;
                            int ii7 = jp1 + (kp1-1)*jd1 + (lp1-1)*jd1*kd1;
                            int ii8 = j   + (kp1-1)*jd1 + (lp1-1)*jd1*kd1;
                            if (iblank[ii1] <= 0 ||
                                iblank[ii2] <= 0 ||
                                iblank[ii3] <= 0 ||
                                iblank[ii4] <= 0 ||
                                iblank[ii5] <= 0 ||
                                iblank[ii6] <= 0 ||
                                iblank[ii7] <= 0 ||
                                iblank[ii8] <= 0) {
                                wherec.nblkpt[i] = 0;
                                std::fprintf(funit(88),
                                    " in hole, rejecting pt. i = %12d mesh %12d\n", i, m);
                            }
                        }
                    }
                } // 12 continue

                int notok = 0;
                for (int i = i1; i <= intrp2.itotal; ++i)
                    if (wherec.nblkpt[i] == 0) notok = notok + 1;

                if (notok > 0) {
                    int iflg = 0;
                    orphan(m, n, nserch, i1, iorph, iflg);
                }

                if (snafu.IHOLEH(m,n) > 0) {
                    int iflg = 0;
                    reserch(m, n, nserch, i1, iskip, iorph, iflg);
                }
            } // 1111 continue

            for (int n = 1; n <= ntot; ++n) {
                int nserch = conec1.MHOLE(m,n);
                for (int nn = 1; nn <= nserch; ++nn) {
                    int m1 = conec1.LHOLE(m,n,nn);
                    if (list[m1] == 0) list[m1] = 1;
                }
            }

            // summary of search routine results
            L6blank();
            L6begin(); L6s("    summary of search routine results for mesh "); L6i(m);
            L6s(" stencils:"); L6end();

            int nintot = 0, nextot = 0, nzerot = 0;

            // extrapolation summary tracking (match Fortran uninitialized use)
            for (int m1 = 1; m1 <= conec1.nmesh; ++m1) {
                if (list[m1] > 0) {
                    double xtrap = 0., xtrap1 = 0.;
                    int jd1 = grid2.mjmax[m1];
                    int kd1 = grid2.mkmax[m1];
                    int ld1 = grid2.mlmax[m1];
                    int nin = 0, nex = 0, norphb = 0, norphc = 0, norphd = 0, nzero = 0;
                    int ixtrap = 0, mxtrap = 0, iflagx = 0;
                    int ixtrap1 = 0, mxtrap1 = 0, iflag1 = 0;

                    for (int i = 1; i <= intrp2.itotal; ++i) {
                        if (wherec.nblkpt[i] == -(m1+mdim)) {
                            int j = intrp2.jb[i];
                            int k = intrp2.kb[i];
                            int l = intrp2.lb[i];
                            double xp = surf2.xbo[i];
                            double yp = surf2.ybo[i];
                            double zp = surf2.zbo[i];
                            int jp = intrp2.jn[i];
                            int kp = intrp2.kn[i];
                            int lp = intrp2.ln[i];
                            double xie  = intrp1.xi[i];
                            double eta  = intrp1.yi[i];
                            double zeta = intrp1.zi[i];

                            if (xie == 0. && eta == 0. && zeta == 0.) {
                                nzero = nzero + 1;
                                if (intrp2.jn[i] == jd1) {
                                    intrp2.jn[i] = intrp2.jn[i] - 1;
                                    intrp1.xi[i] = 1.;
                                    jp = intrp2.jn[i];
                                    xie = intrp1.xi[i];
                                }
                                if (intrp2.kn[i] == kd1) {
                                    intrp2.kn[i] = intrp2.kn[i] - 1;
                                    intrp1.yi[i] = 1.;
                                    kp = intrp2.kn[i];
                                    eta = intrp1.yi[i];
                                }
                                if (intrp2.ln[i] == ld1) {
                                    if (ld1 > 1) {
                                        intrp2.ln[i] = intrp2.ln[i] - 1;
                                    }
                                    intrp1.zi[i] = 1.;
                                    lp = intrp2.ln[i];
                                    zeta = intrp1.zi[i];
                                }

                                if (intrp2.jn[i] == 1 || intrp2.jn[i] == jd1-1 ||
                                    intrp2.kn[i] == 1 || intrp2.kn[i] == kd1-1 ||
                                    intrp2.ln[i] == 1 || intrp2.ln[i] == ld1-1) norphb = norphb + 1;

                                std::fprintf(funit(8),
                                    " for the mesh %12d fringe point  j,k,l = %12d%12d%12d\n", m, j, k, l);
                                std::fprintf(funit(8),
                                    "   will use nearest point in mesh %12d: j,k,l = %12d%12d%12d\n", m1, jp, kp, lp);
                                std::fprintf(funit(8),
                                    "     with xie,eta,zeta = %g%g%g\n", xie, eta, zeta);

                            } else if (xie >= -tol.epsc && xie <= 1.+tol.epsc &&
                                       eta >= -tol.epsc && eta <= 1.+tol.epsc &&
                                       zeta >= -tol.epsc && zeta <= 1.+tol.epsc) {
                                nin = nin + 1;
                            }

                            if (xie < -tol.epsc || xie > 1.+tol.epsc ||
                                eta < -tol.epsc || eta > 1.+tol.epsc ||
                                zeta < -tol.epsc || zeta > 1.+tol.epsc) {
                                int iflagg = 0;
                                pltpt.norph = pltpt.norph + 1;
                                pltpt.xorph[pltpt.norph] = xp;
                                pltpt.yorph[pltpt.norph] = yp;
                                pltpt.zorph[pltpt.norph] = zp;
                                nex = nex + 1;
                                if (intrp2.jn[i] == 1 || intrp2.jn[i] == jd1-1 ||
                                    intrp2.kn[i] == 1 || intrp2.kn[i] == kd1-1 ||
                                    intrp2.ln[i] == 1 || intrp2.ln[i] == ld1-1) norphb = norphb + 1;

                                std::fprintf(funit(8),
                                    " for the mesh %12d fringe point  j,k,l = %12d%12d%12d\n", m, j, k, l);
                                std::fprintf(funit(8),
                                    "   will use extrapolation from mesh %12d point j,k,l = %12d%12d%12d\n", m1, jp, kp, lp);
                                std::fprintf(funit(8),
                                    "     with xie,eta,zeta = %g%g%g\n", xie, eta, zeta);

                                if (xie < -tol.epsc || xie > 1.+tol.epsc) {
                                    iflagg = iflagg + 1;
                                    if (xie > xtrap) { xtrap = xie; ixtrap = i; mxtrap = m1; iflagx = 1; }
                                    if (xie < xtrap1) { xtrap1 = xie; ixtrap1 = i; mxtrap1 = m1; iflag1 = 1; }
                                }
                                if (eta < -tol.epsc || eta > 1.+tol.epsc) {
                                    iflagg = iflagg + 1;
                                    if (eta > xtrap) { xtrap = eta; ixtrap = i; mxtrap = m1; iflagx = 2; }
                                    if (eta < xtrap1) { xtrap1 = eta; ixtrap1 = i; mxtrap1 = m1; iflag1 = 2; }
                                }
                                if (zeta < -tol.epsc || zeta > 1.+tol.epsc) {
                                    iflagg = iflagg + 1;
                                    if (zeta > xtrap) { xtrap = zeta; ixtrap = i; mxtrap = m1; iflagx = 3; }
                                    if (zeta < xtrap1) { xtrap1 = zeta; ixtrap1 = i; mxtrap1 = m1; iflag1 = 3; }
                                }
                                if (xie < -0.5-tol.epsc || xie > 1.5+tol.epsc) norphc = norphc + 1;
                                if (eta < -0.5-tol.epsc || eta > 1.5+tol.epsc) norphc = norphc + 1;
                                if (zeta < -0.5-tol.epsc || zeta > 1.5+tol.epsc) norphc = norphc + 1;
                                if (iflagg > 1) norphd = norphd + 1;
                            }
                        }
                    } // 115 continue

                    L6blank();
                    L6begin(); L6s("      "); L6i(nin); L6s(" fringe points of mesh "); L6i(m);
                    L6s(" are interpolated from mesh "); L6i(m1); L6end();
                    if (nex > 0) {
                        L6begin(); L6s("      "); L6i(nex); L6s(" orphaned fringe points of mesh "); L6i(m);
                        L6s(" are extrapolated from mesh "); L6i(m1); L6end();
                        L6begin(); L6s("      "); L6i(norphc); L6s(" of these orphans have extrapolation coefficients <-0.5 or >1.5"); L6end();
                        if (norphd > 0) {
                            L6begin(); L6s("      "); L6i(norphd); L6s(" of these orphans are extrapolated in more than one direction"); L6end();
                        }
                        if (xtrap > 1.+tol.epsc) {
                            L6begin(); L6s("        the maximum extrapolation coefficient is "); L6s(flag[iflagx]); L6s(" = ");
                            std::fprintf(funit(6), "%g", xtrap); L6end();
                            L6begin(); L6s("        to the fringe pt "); L6i(intrp2.jb[ixtrap]); L6i(intrp2.kb[ixtrap]);
                            L6i(intrp2.lb[ixtrap]); L6s(" from the mesh "); L6i(mxtrap); L6s(" target pt ");
                            L6i(intrp2.jn[ixtrap]); L6i(intrp2.kn[ixtrap]); L6i(intrp2.ln[ixtrap]); L6end();
                        }
                        if (xtrap1 < -tol.epsc) {
                            L6begin(); L6s("        the minimum extrapolation coefficient is "); L6s(flag[iflag1]); L6s(" = ");
                            std::fprintf(funit(6), "%g", xtrap1); L6end();
                            L6begin(); L6s("        to the fringe pt "); L6i(intrp2.jb[ixtrap1]); L6i(intrp2.kb[ixtrap1]);
                            L6i(intrp2.lb[ixtrap1]); L6s(" from the mesh "); L6i(mxtrap1); L6s(" target pt ");
                            L6i(intrp2.jn[ixtrap1]); L6i(intrp2.kn[ixtrap1]); L6i(intrp2.ln[ixtrap1]); L6end();
                        }
                        L6begin(); L6s("        check file 8 for more details"); L6end();
                    }
                    if (nzero > 0) {
                        L6begin(); L6s("      "); L6i(nzero); L6s(" orphaned fringe points of mesh "); L6i(m);
                        L6s(" use the nearest point in mesh "); L6i(m1); L6end();
                        L6begin(); L6s("        check file 8 for more details"); L6end();
                    }
                    if (nex+nzero-norphb != 0) {
                        L6begin(); L6s("        WARNING: "); L6i(nex-norphb); L6s(" of these orphans are extrapolated from an"); L6end();
                        L6begin(); L6s("          interior point of mesh "); L6i(m1); L6s("...possible miscue in search routine"); L6end();
                    }

                    nintot = nintot + nin;
                    nextot = nextot + nex;
                    nzerot = nzerot + nzero;
                }
            } // 114 continue

            if (nintot+nextot+nzerot != intrp2.itotal) {
                L6begin(); L6s("  stopping...unable to compute "); L6i(intrp2.itotal-nintot-nextot-nzerot);
                L6s(" stencils for mesh "); L6i(m); L6s(" fringe points"); L6end();
                L6begin(); L6s("     possible miscue in search routine"); L6end();
                std::exit(0);
            }

            // put data back in original order
            int itotl1 = intrp2.itotal;
            for (int i = 1; i <= itotl1; ++i) {
                jbtmp[i] = intrp2.jb[i];
                kbtmp[i] = intrp2.kb[i];
                lbtmp[i] = intrp2.lb[i];
                jntmp[i] = intrp2.jn[i];
                kntmp[i] = intrp2.kn[i];
                lntmp[i] = intrp2.ln[i];
                xitmp[i] = intrp1.xi[i];
                yitmp[i] = intrp1.yi[i];
                zitmp[i] = intrp1.zi[i];
            }

            for (int m1 = 1; m1 <= conec1.nmesh; ++m1) {
                if (list[m1] > 0) {
                    int ii = 0;
                    for (int i = 1; i <= itotl1; ++i) {
                        if (wherec.nblkpt[i] == -(m1+mdim)) {
                            ii = ii + 1;
                            intrp2.jb[ii] = jbtmp[i];
                            intrp2.kb[ii] = kbtmp[i];
                            intrp2.lb[ii] = lbtmp[i];
                            intrp2.jn[ii] = jntmp[i];
                            intrp2.kn[ii] = kntmp[i];
                            intrp2.ln[ii] = lntmp[i];
                            intrp1.xi[ii] = xitmp[i];
                            intrp1.yi[ii] = yitmp[i];
                            intrp1.zi[ii] = zitmp[i];
                        }
                    }
                    intrp2.itotal = ii;
                    int icase = 1;
                    setptr(m, m1, icase);
                }
            }
        }
    } // 1313 continue

    trace1.itrace = -1;
}

// ===========================================================================
//  outer — compute interpolation data for outer boundaries of mesh m
// ===========================================================================
void outer()
{
    int* iblank = book3.iblank;
    double* x = grid1.x; double* y = grid1.y; double* z = grid1.z;

    static double xitmp[idim+1], yitmp[idim+1], zitmp[idim+1];
    static int jbtmp[idim+1], kbtmp[idim+1], lbtmp[idim+1];
    static int jntmp[idim+1], kntmp[idim+1], lntmp[idim+1];
    int list[mdim+1];
    int iskip[idim+1]; (void)iskip;

    const char* flag[4];
    flag[1] = "xie ";
    flag[2] = "eta ";
    flag[3] = "zeta";

    trace1.itrace = -1;

    for (int m = 1; m <= conec1.nmesh; ++m) {
        int jd = grid2.mjmax[m];
        int kd = grid2.mkmax[m];
        int ld = grid2.mlmax[m];
        int ntot = conec1.noutr[m];
        intrp2.itotal = 0;

        if (ntot > 0) {
            for (int mm = 1; mm <= conec1.nmesh; ++mm) list[mm] = 0;

            L6blank();
            L6blank();
            L6begin(); L6s("  *** beginning determination of stencils for boundary points of mesh "); L6i(m); L6s(" ***"); L6end();
            L6blank();

            for (int n = 1; n <= ntot; ++n) {
                int nnn = 1;
                int nserch = conec1.MOUTR(m,nnn);
                int iorph = snafu.IORPHO(m,nnn);

                getgrd(m, x, y, z, image.jimage, image.kimage, image.limage, jd, kd, ld);
                getibl(m, iblank, jd, kd, ld);

                int i1 = intrp2.itotal + 1;   // == obtot's i1 output
                obtot(m, x, y, z, iblank, jd, kd, ld, i1);

                for (int i = i1; i <= intrp2.itotal; ++i) wherec.nblkpt[i] = 0;

                for (int nn = 1; nn <= nserch; ++nn) {
                    nnn = 1;
                    int m1 = conec1.LOUTR(m,nnn,nn);
                    int jd1 = grid2.mjmax[m1];
                    int kd1 = grid2.mkmax[m1];
                    int ld1 = grid2.mlmax[m1];
                    getgrd(m1, x, y, z, image.jimage, image.kimage, image.limage, jd1, kd1, ld1);
                    intpt(0, image.jimage, image.kimage, image.limage, jd, kd, ld,
                          m1, x, y, z, jd1, kd1, ld1, m, i1);
                    getibl(m1, iblank, jd1, kd1, ld1);
                    for (int i = i1; i <= intrp2.itotal; ++i) {
                        if (wherec.nblkpt[i] == -(m1+mdim)) {
                            int j = intrp2.jn[i];
                            int k = intrp2.kn[i];
                            int l = intrp2.ln[i];
                            int jp1 = std::min(j+1, jd1);
                            int kp1 = std::min(k+1, kd1);
                            int lp1 = std::min(l+1, ld1);
                            int ii1 = j   + (k-1)*jd1   + (l-1)*jd1*kd1;
                            int ii2 = jp1 + (k-1)*jd1   + (l-1)*jd1*kd1;
                            int ii3 = jp1 + (kp1-1)*jd1 + (l-1)*jd1*kd1;
                            int ii4 = j   + (kp1-1)*jd1 + (l-1)*jd1*kd1;
                            int ii5 = j   + (k-1)*jd1   + (lp1-1)*jd1*kd1;
                            int ii6 = jp1 + (k-1)*jd1   + (lp1-1)*jd1*kd1;
                            int ii7 = jp1 + (kp1-1)*jd1 + (lp1-1)*jd1*kd1;
                            int ii8 = j   + (kp1-1)*jd1 + (lp1-1)*jd1*kd1;
                            if (iblank[ii1] <= 0 ||
                                iblank[ii2] <= 0 ||
                                iblank[ii3] <= 0 ||
                                iblank[ii4] <= 0 ||
                                iblank[ii5] <= 0 ||
                                iblank[ii6] <= 0 ||
                                iblank[ii7] <= 0 ||
                                iblank[ii8] <= 0) {
                                wherec.nblkpt[i] = 0;
                                std::fprintf(funit(88),
                                    " in outer, rejecting pt. i = %12d mesh %12d\n", i, m);
                            }
                        }
                    }
                } // 21 continue

                int notok = 0;
                for (int i = i1; i <= intrp2.itotal; ++i)
                    if (wherec.nblkpt[i] == 0) notok = notok + 1;

                if (notok > 0) {
                    int iflg = 1;
                    orphan(m, n, nserch, i1, iorph, iflg);
                }

                nnn = 1;
                if (snafu.IHOLEO(m,nnn) > 0) {
                    int iflg = 1;
                    reserch(m, n, nserch, i1, iskip, iorph, iflg);
                }
            } // 1111 continue

            for (int n = 1; n <= ntot; ++n) {
                int nnn = 1;
                int nserch = conec1.MOUTR(m,nnn);
                for (int nn = 1; nn <= nserch; ++nn) {
                    nnn = 1;
                    int m1 = conec1.LOUTR(m,nnn,nn);
                    if (list[m1] == 0) list[m1] = 1;
                }
            }

            L6blank();
            L6begin(); L6s("    summary of search routine results for mesh "); L6i(m); L6s(" stencils:"); L6end();

            int nintot = 0, nextot = 0, nzerot = 0;

            for (int m1 = 1; m1 <= conec1.nmesh; ++m1) {
                if (list[m1] > 0) {
                    double xtrap = 0., xtrap1 = 0.;
                    int jd1 = grid2.mjmax[m1];
                    int kd1 = grid2.mkmax[m1];
                    int ld1 = grid2.mlmax[m1];
                    int nin = 0, nex = 0, norphb = 0, norphc = 0, norphd = 0, nzero = 0;
                    int ixtrap = 0, mxtrap = 0, iflagx = 0;
                    int ixtrap1 = 0, mxtrap1 = 0, iflag1 = 0;

                    for (int i = 1; i <= intrp2.itotal; ++i) {
                        if (wherec.nblkpt[i] == -(m1+mdim)) {
                            int j = intrp2.jb[i];
                            int k = intrp2.kb[i];
                            int l = intrp2.lb[i];
                            double xp = surf2.xbo[i];
                            double yp = surf2.ybo[i];
                            double zp = surf2.zbo[i];
                            int jp = intrp2.jn[i];
                            int kp = intrp2.kn[i];
                            int lp = intrp2.ln[i];
                            double xie  = intrp1.xi[i];
                            double eta  = intrp1.yi[i];
                            double zeta = intrp1.zi[i];

                            if (xie == 0. && eta == 0. && zeta == 0.) {
                                nzero = nzero + 1;
                                if (intrp2.jn[i] == jd1) {
                                    intrp2.jn[i] = intrp2.jn[i] - 1;
                                    intrp1.xi[i] = 1.;
                                    jp = intrp2.jn[i];
                                    xie = intrp1.xi[i];
                                }
                                if (intrp2.kn[i] == kd1) {
                                    intrp2.kn[i] = intrp2.kn[i] - 1;
                                    intrp1.yi[i] = 1.;
                                    kp = intrp2.kn[i];
                                    eta = intrp1.yi[i];
                                }
                                if (intrp2.ln[i] == ld1) {
                                    if (ld1 > 1) {
                                        intrp2.ln[i] = intrp2.ln[i] - 1;
                                    }
                                    intrp1.zi[i] = 1.;
                                    lp = intrp2.ln[i];
                                    zeta = intrp1.zi[i];
                                }

                                if (intrp2.jn[i] == 1 || intrp2.jn[i] == jd1-1 ||
                                    intrp2.kn[i] == 1 || intrp2.kn[i] == kd1-1 ||
                                    intrp2.ln[i] == 1 || intrp2.ln[i] == ld1-1) norphb = norphb + 1;

                                std::fprintf(funit(8),
                                    " for the mesh %12d boundary point  j,k,l = %12d%12d%12d\n", m, j, k, l);
                                std::fprintf(funit(8),
                                    "   will use nearest point in mesh %12d: j,k,l = %12d%12d%12d\n", m1, jp, kp, lp);
                                std::fprintf(funit(8),
                                    "     with xie,eta,zeta = %g%g%g\n", xie, eta, zeta);

                            } else if (xie >= -tol.epsc && xie <= 1.+tol.epsc &&
                                       eta >= -tol.epsc && eta <= 1.+tol.epsc &&
                                       zeta >= -tol.epsc && zeta <= 1.+tol.epsc) {
                                nin = nin + 1;
                            }

                            if (xie < -tol.epsc || xie > 1.+tol.epsc ||
                                eta < -tol.epsc || eta > 1.+tol.epsc ||
                                zeta < -tol.epsc || zeta > 1.+tol.epsc) {
                                int iflagg = 0;
                                pltpt.norph = pltpt.norph + 1;
                                pltpt.xorph[pltpt.norph] = xp;
                                pltpt.yorph[pltpt.norph] = yp;
                                pltpt.zorph[pltpt.norph] = zp;
                                nex = nex + 1;
                                if (intrp2.jn[i] == 1 || intrp2.jn[i] == jd1-1 ||
                                    intrp2.kn[i] == 1 || intrp2.kn[i] == kd1-1 ||
                                    intrp2.ln[i] == 1 || intrp2.ln[i] == ld1-1) norphb = norphb + 1;

                                std::fprintf(funit(8),
                                    " for the mesh %12d boundary point  j,k,l = %12d%12d%12d\n", m, j, k, l);
                                std::fprintf(funit(8),
                                    "   will use extrapolation from mesh %12d point j,k,l = %12d%12d%12d\n", m1, jp, kp, lp);
                                std::fprintf(funit(8),
                                    "     with xie,eta,zeta = %g%g%g\n", xie, eta, zeta);

                                if (xie < -tol.epsc || xie > 1.+tol.epsc) {
                                    iflagg = iflagg + 1;
                                    if (xie > xtrap) { xtrap = xie; ixtrap = i; mxtrap = m1; iflagx = 1; }
                                    if (xie < xtrap1) { xtrap1 = xie; ixtrap1 = i; mxtrap1 = m1; iflag1 = 1; }
                                }
                                if (eta < -tol.epsc || eta > 1.+tol.epsc) {
                                    iflagg = iflagg + 1;
                                    if (eta > xtrap) { xtrap = eta; ixtrap = i; mxtrap = m1; iflagx = 2; }
                                    if (eta < xtrap1) { xtrap1 = eta; ixtrap1 = i; mxtrap1 = m1; iflag1 = 2; }
                                }
                                if (zeta < -tol.epsc || zeta > 1.+tol.epsc) {
                                    iflagg = iflagg + 1;
                                    if (zeta > xtrap) { xtrap = zeta; ixtrap = i; mxtrap = m1; iflagx = 3; }
                                    if (zeta < xtrap1) { xtrap1 = zeta; ixtrap1 = i; mxtrap1 = m1; iflag1 = 3; }
                                }
                                if (xie < -0.5-tol.epsc || xie > 1.5+tol.epsc) norphc = norphc + 1;
                                if (eta < -0.5-tol.epsc || eta > 1.5+tol.epsc) norphc = norphc + 1;
                                if (zeta < -0.5-tol.epsc || zeta > 1.5+tol.epsc) norphc = norphc + 1;
                                if (iflagg > 1) norphd = norphd + 1;
                            }
                        }
                    } // 115 continue

                    L6blank();
                    L6begin(); L6s("      "); L6i(nin); L6s(" boundary points of mesh "); L6i(m);
                    L6s(" are interpolated from mesh "); L6i(m1); L6end();
                    if (nex > 0) {
                        L6begin(); L6s("      "); L6i(nex); L6s(" orphaned boundary points of mesh "); L6i(m);
                        L6s(" are extrapolated from mesh "); L6i(m1); L6end();
                        L6begin(); L6s("      "); L6i(norphc); L6s(" of these orphans have extrapolation coefficients <-0.5 or >1.5"); L6end();
                        if (norphd > 0) {
                            L6begin(); L6s("      "); L6i(norphd); L6s(" of these orphans are extrapolated in more than one direction"); L6end();
                        }
                        if (xtrap > 1.+tol.epsc) {
                            L6begin(); L6s("        the maximum extrapolation coefficient is "); L6s(flag[iflagx]); L6s(" = ");
                            std::fprintf(funit(6), "%g", xtrap); L6end();
                            L6begin(); L6s("        to the boundary pt "); L6i(intrp2.jb[ixtrap]); L6i(intrp2.kb[ixtrap]);
                            L6i(intrp2.lb[ixtrap]); L6s(" from the mesh "); L6i(mxtrap); L6s(" target pt ");
                            L6i(intrp2.jn[ixtrap]); L6i(intrp2.kn[ixtrap]); L6i(intrp2.ln[ixtrap]); L6end();
                        }
                        if (xtrap1 < -tol.epsc) {
                            L6begin(); L6s("        the minimum extrapolation coefficient is "); L6s(flag[iflag1]); L6s(" = ");
                            std::fprintf(funit(6), "%g", xtrap1); L6end();
                            L6begin(); L6s("        to the boundary pt "); L6i(intrp2.jb[ixtrap1]); L6i(intrp2.kb[ixtrap1]);
                            L6i(intrp2.lb[ixtrap1]); L6s(" from the mesh "); L6i(mxtrap1); L6s(" target pt ");
                            L6i(intrp2.jn[ixtrap1]); L6i(intrp2.kn[ixtrap1]); L6i(intrp2.ln[ixtrap1]); L6end();
                        }
                        L6begin(); L6s("        check file 8 for more details"); L6end();
                    }
                    if (nzero > 0) {
                        L6begin(); L6s("      "); L6i(nzero); L6s(" orphaned boundary points of mesh "); L6i(m);
                        L6s(" use the nearest point in mesh "); L6i(m1); L6end();
                        L6begin(); L6s("        check file 8 for more details"); L6end();
                    }
                    if (nex+nzero-norphb != 0) {
                        L6begin(); L6s("        WARNING: "); L6i(nex+nzero-norphb); L6s(" of these orphans use extrap./zeroth order interp."); L6end();
                        L6begin(); L6s("          from an interior point of mesh "); L6i(m1); L6s("...possible miscue in search routine"); L6end();
                    }

                    nintot = nintot + nin;
                    nextot = nextot + nex;
                    nzerot = nzerot + nzero;
                }
            } // 114 continue

            if (nintot+nextot+nzerot != intrp2.itotal) {
                L6begin(); L6s("  stopping...unable to compute "); L6i(intrp2.itotal-nintot-nextot);
                L6s(" stencils for mesh "); L6i(m); L6s(" boundary points"); L6end();
                L6begin(); L6s("     possible miscue in search routine"); L6end();
                std::exit(0);
            }

            int itotl1 = intrp2.itotal;
            for (int i = 1; i <= itotl1; ++i) {
                jbtmp[i] = intrp2.jb[i];
                kbtmp[i] = intrp2.kb[i];
                lbtmp[i] = intrp2.lb[i];
                jntmp[i] = intrp2.jn[i];
                kntmp[i] = intrp2.kn[i];
                lntmp[i] = intrp2.ln[i];
                xitmp[i] = intrp1.xi[i];
                yitmp[i] = intrp1.yi[i];
                zitmp[i] = intrp1.zi[i];
            }
            for (int m1 = 1; m1 <= conec1.nmesh; ++m1) {
                if (list[m1] > 0) {
                    int ii = 0;
                    for (int i = 1; i <= itotl1; ++i) {
                        if (wherec.nblkpt[i] == -(m1+mdim)) {
                            ii = ii + 1;
                            intrp2.jb[ii] = jbtmp[i];
                            intrp2.kb[ii] = kbtmp[i];
                            intrp2.lb[ii] = lbtmp[i];
                            intrp2.jn[ii] = jntmp[i];
                            intrp2.kn[ii] = kntmp[i];
                            intrp2.ln[ii] = lntmp[i];
                            intrp1.xi[ii] = xitmp[i];
                            intrp1.yi[ii] = yitmp[i];
                            intrp1.zi[ii] = zitmp[i];
                        }
                    }
                    intrp2.itotal = ii;
                    int icase = 2;
                    setptr(m, m1, icase);
                }
            }
        }
    } // 31 continue

    // update and store iblank array for mesh m
    for (int m = 1; m <= conec1.nmesh; ++m) {
        int ntot = conec1.noutr[m];
        if (ntot > 0) {
            for (int n = 1; n <= ntot; ++n) {
                int jd = grid2.mjmax[m];
                int kd = grid2.mkmax[m];
                int ld = grid2.mlmax[m];
                getibl(m, iblank, jd, kd, ld);
                obibl(m, n, iblank, jd, kd, ld);
                putibl(m, iblank, jd, kd, ld);
            }
        }
    }

    trace1.itrace = -1;
}
