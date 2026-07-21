// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// mag_common.h — COMMON-block globals + Fortran-array helpers for the maggie
// C++ port. Mirrors CFL3D/header/mag1.h parameters and every COMMON block in
// maggie.F. Column-major, 1-based indexing throughout (to match Fortran).
#pragma once
#include "../../common/tools_common.h"
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <string>
#include <algorithm>

// ---- mag1.h parameters -----------------------------------------------------
static const int mdim   = 14;
static const int jdim   = 300;
static const int kdim   = 300;
static const int ldim   = 2;
static const int mhldim = 10;
static const int ibdim  = 2000;
static const int idim   = 8000;
static const int ndim   = 300;
static const int ipmax  = 28;
static const int iwrdim = idim + ibdim;
static const int mlen   = jdim * kdim * ldim;

// ---- 1-based fixed array wrapper (Fortran DIMENSION a(n)) -------------------
// F1[i] with i in 1..N. Underlying storage is N+1 with [0] unused.
template <class T, int N>
struct F1 {
    T d[N + 1];
    T& operator()(int i) { return d[i]; }
    const T& operator()(int i) const { return d[i]; }
    T* raw() { return d; }               // d[0..N]
    T* base1() { return d; }             // pass &base1()[0] -> element(1) at [1]
    void zero() { std::memset(d, 0, sizeof(d)); }
};

// ---- column-major N-D views over a flat 1-based pointer --------------------
// Fortran a(i1,i2) with dims (n1,n2): offset = i1 + n1*(i2-1), 1-based flat.
// We store a pointer p where p[1] is the first element (1-based flat index).
template <class T>
struct View2 {
    T* p; int n1, n2;
    View2(T* p_, int n1_, int n2_) : p(p_), n1(n1_), n2(n2_) {}
    T& operator()(int i, int j) { return p[i + n1 * (j - 1)]; }
};
template <class T>
struct View3 {
    T* p; int n1, n2, n3;
    View3(T* p_, int n1_, int n2_, int n3_) : p(p_), n1(n1_), n2(n2_), n3(n3_) {}
    T& operator()(int i, int j, int k) { return p[i + n1 * (j - 1) + (long)n1 * n2 * (k - 1)]; }
};

// index helper: 1-based flat index of a(i,j,k) with dims (n1,n2,n3)
inline long fidx3(int i, int j, int k, int n1, int n2) {
    return (long)i + (long)n1 * (j - 1) + (long)n1 * n2 * (k - 1);
}
inline long fidx2(int i, int j, int n1) { return (long)i + (long)n1 * (j - 1); }

// ===========================================================================
//  COMMON blocks — one global struct instance each. Field names & shapes
//  match maggie.F exactly. All arrays are 1-based (size dim+1, [0] unused);
//  multi-dim arrays are stored column-major in a flat buffer and accessed
//  via a per-COMMON accessor, OR as explicit F1 with manual index where the
//  original code only uses linear COMMON storage.
// ===========================================================================

// /book1/ ipntr(mdim,mhldim*mdim), npntr(mdim,mhldim*mdim),
//         mhbs(mdim,mdim), mobs(mdim,mdim), nsets(mdim)
struct Book1 {
    int ipntr[mdim*(mhldim*mdim) + 1];
    int npntr[mdim*(mhldim*mdim) + 1];
    int mhbs [mdim*mdim + 1];
    int mobs [mdim*mdim + 1];
    int nsets[mdim + 1];
    int& IPNTR(int i,int j){ return ipntr[fidx2(i,j,mdim)]; }
    int& NPNTR(int i,int j){ return npntr[fidx2(i,j,mdim)]; }
    int& MHBS (int i,int j){ return mhbs [fidx2(i,j,mdim)]; }
    int& MOBS (int i,int j){ return mobs [fidx2(i,j,mdim)]; }
};
extern Book1 book1;

// /book2/ dxint(idim),dyint(idim),dzint(idim),ibpts(mdim),
//         jbpt(idim),kbpt(idim),lbpt(idim),ji(idim),ki(idim),li(idim)
struct Book2 {
    double dxint[idim+1], dyint[idim+1], dzint[idim+1];
    int    ibpts[mdim+1];
    int    jbpt[idim+1], kbpt[idim+1], lbpt[idim+1];
    int    ji[idim+1], ki[idim+1], li[idim+1];
};
extern Book2 book2;

// /book3/ iblank(mlen)
struct Book3 { int iblank[mlen+1]; };
extern Book3 book3;

// /bound1/ ihbtyp(mdim,mdim),jh1(ipmax),jh2(ipmax),kh1,kh2,lh1,lh2,
//          ip1(ipmax),ip2(ipmax),mh(ipmax),ihplt(ipmax)
struct Bound1 {
    int ihbtyp[mdim*mdim + 1];
    int jh1[ipmax+1], jh2[ipmax+1], kh1[ipmax+1], kh2[ipmax+1],
        lh1[ipmax+1], lh2[ipmax+1], ip1[ipmax+1], ip2[ipmax+1],
        mh[ipmax+1], ihplt[ipmax+1];
    int& IHBTYP(int i,int j){ return ihbtyp[fidx2(i,j,mdim)]; }
};
extern Bound1 bound1;

// /bound2/ iobtyp(mdim),nobtyp(mdim),jo1(ipmax),jo2,ko1,ko2,lo1,lo2
struct Bound2 {
    int iobtyp[mdim+1], nobtyp[mdim+1];
    int jo1[ipmax+1], jo2[ipmax+1], ko1[ipmax+1], ko2[ipmax+1],
        lo1[ipmax+1], lo2[ipmax+1];
};
extern Bound2 bound2;

// /grdlim/ jsrs(mdim),jsre,ksrs,ksre,lsrs,lsre
struct Grdlim {
    int jsrs[mdim+1], jsre[mdim+1], ksrs[mdim+1], ksre[mdim+1],
        lsrs[mdim+1], lsre[mdim+1];
};
extern Grdlim grdlim;

// /grid1/ x(mlen),y(mlen),z(mlen)
struct Grid1 { double x[mlen+1], y[mlen+1], z[mlen+1]; };
extern Grid1 grid1;

// /grid2/ mjmax(mdim),mkmax(mdim),mlmax(mdim)
struct Grid2 { int mjmax[mdim+1], mkmax[mdim+1], mlmax[mdim+1]; };
extern Grid2 grid2;

// /conec1/ nmesh,nhole(mdim),noutr(mdim),mhole(mdim,mhldim),moutr(mdim,mhldim),
//          lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
struct Conec1 {
    int nmesh;
    int nhole[mdim+1], noutr[mdim+1];
    int mhole[mdim*mhldim + 1], moutr[mdim*mhldim + 1];
    int lhole[mdim*mhldim*mdim + 1], loutr[mdim*mhldim*mdim + 1];
    int& MHOLE(int i,int j){ return mhole[fidx2(i,j,mdim)]; }
    int& MOUTR(int i,int j){ return moutr[fidx2(i,j,mdim)]; }
    int& LHOLE(int i,int j,int k){ return lhole[fidx3(i,j,k,mdim,mhldim)]; }
    int& LOUTR(int i,int j,int k){ return loutr[fidx3(i,j,k,mdim,mhldim)]; }
};
extern Conec1 conec1;

// /image/ jimage(mlen),kimage(mlen),limage(mlen)
struct Image { int jimage[mlen+1], kimage[mlen+1], limage[mlen+1]; };
extern Image image;

// /igrdtyp/ ip3dgrd,ialph
struct Igrdtyp { int ip3dgrd, ialph; };
extern Igrdtyp igrdtyp;

// /iseq/ iseq(mdim),icount(mdim)
struct Iseq { int iseq[mdim+1], icount[mdim+1]; };
extern Iseq iseqc;

// /diagno/ jp(mdim*3+3),kp(mdim*3+3),lp(mdim*3+3),iplt3d
struct Diagno { int jp[mdim*3+3+1], kp[mdim*3+3+1], lp[mdim*3+3+1], iplt3d; };
extern Diagno diagno;

// /snafu/ iholeh(mdim,mhldim),iholeo(mdim,mhldim),iorphh(mdim,mhldim),iorpho(mdim,mhldim)
struct Snafu {
    int iholeh[mdim*mhldim+1], iholeo[mdim*mhldim+1],
        iorphh[mdim*mhldim+1], iorpho[mdim*mhldim+1];
    int& IHOLEH(int i,int j){ return iholeh[fidx2(i,j,mdim)]; }
    int& IHOLEO(int i,int j){ return iholeo[fidx2(i,j,mdim)]; }
    int& IORPHH(int i,int j){ return iorphh[fidx2(i,j,mdim)]; }
    int& IORPHO(int i,int j){ return iorpho[fidx2(i,j,mdim)]; }
};
extern Snafu snafu;

// /pltpt/ xbnd(idim)..zill(idim),nbnd,norph,nill
struct Pltpt {
    double xbnd[idim+1], ybnd[idim+1], zbnd[idim+1];
    double xorph[idim+1], yorph[idim+1], zorph[idim+1];
    double xill[idim+1], yill[idim+1], zill[idim+1];
    int nbnd, norph, nill;
};
extern Pltpt pltpt;

// /files/ grid,outpt,ovrlp,plt3d  (character*80)
struct Files { char grid[81], outpt[81], ovrlp[81], plt3d[81]; };
extern Files files;

// /tol/ epsc
struct Tol { double epsc; };
extern Tol tol;

// /intrp1/ xi(idim),yi(idim),zi(idim)
struct Intrp1 { double xi[idim+1], yi[idim+1], zi[idim+1]; };
extern Intrp1 intrp1;

// /intrp2/ jb(idim),kb(idim),lb(idim),jn(idim),kn(idim),ln(idim),itotal
struct Intrp2 {
    int jb[idim+1], kb[idim+1], lb[idim+1],
        jn[idim+1], kn[idim+1], ln[idim+1], itotal;
};
extern Intrp2 intrp2;

// /where/ nblkpt(idim)
struct Where { int nblkpt[idim+1]; };
extern Where wherec;

// /surf/ xb(ibdim),yb(ibdim),zb(ibdim),ibmax
struct Surf { double xb[ibdim+1], yb[ibdim+1], zb[ibdim+1]; int ibmax; };
extern Surf surf;

// /surf2/ xbo(idim),ybo(idim),zbo(idim)
struct Surf2 { double xbo[idim+1], ybo[idim+1], zbo[idim+1]; };
extern Surf2 surf2;

// /norm/ vnx(ibdim),vny(ibdim),vnz(ibdim)
struct Norm { double vnx[ibdim+1], vny[ibdim+1], vnz[ibdim+1]; };
extern Norm normc;

// /trace1/ itrace
struct Trace1 { int itrace; };
extern Trace1 trace1;

// /chkst/ nwr
struct Chkst { int nwr; };
extern Chkst chkst;

// /chkpt/ xcc0(idim)..zscale(idim)
struct Chkpt {
    double xcc0[idim+1], ycc0[idim+1], zcc0[idim+1];
    double xccint[idim+1], yccint[idim+1], zccint[idim+1];
    double xscale[idim+1], yscale[idim+1], zscale[idim+1];
};
extern Chkpt chkpt;

// /dimen/ ltest(mdim),jtest(mdim),ktest(mdim)
struct Dimen { int ltest[mdim+1], jtest[mdim+1], ktest[mdim+1]; };
extern Dimen dimenc;

// /work/ rm(iwrdim)   /iwork/ ibck(iwrdim)
struct Work  { double rm[iwrdim+1]; };
struct Iwork { int ibck[iwrdim+1]; };
extern Work workc;
extern Iwork iworkc;

// /temp/ temp(ndim,ndim,6)
struct Temp {
    double* p; // heap: too big for stack (300*300*6 doubles ~ 4.3MB)
    Temp();
    ~Temp();
    double& operator()(int i,int j,int k){ return p[fidx3(i,j,k,ndim,ndim)]; }
};
extern Temp tempc;

// /temp1/ x1(jdim,kdim,ldim),y1,z1  — heap-backed (huge)
struct Temp1 {
    double *x1p,*y1p,*z1p;
    Temp1(); ~Temp1();
    double& X1(int i,int j,int k){ return x1p[fidx3(i,j,k,jdim,kdim)]; }
    double& Y1(int i,int j,int k){ return y1p[fidx3(i,j,k,jdim,kdim)]; }
    double& Z1(int i,int j,int k){ return z1p[fidx3(i,j,k,jdim,kdim)]; }
};
extern Temp1 temp1c;

// ---- Fortran runtime unit table (unit N -> FILE*) --------------------------
// maggie opens units 2,6,9,10,30,35 and uses 30 as scratch for temp_*.n files.
FILE* funit(int u);                 // returns current FILE* for a unit
void  fopen_unit(int u, const char* name, const char* mode);
void  frewind(int u);
void  fclose_unit(int u);

// unformatted record helpers bound to a unit
void uwrite(int u, const void* p, size_t n);
void uread (int u, void* p, size_t n);

// ---- subroutine prototypes ----
// The exact C++ signatures for every translated maggie subroutine live in
// mag_protos.h (assembled from the actual translations so signatures cannot
// drift from their definitions). Every implementation .cpp includes it.
#include "mag_protos.h"
