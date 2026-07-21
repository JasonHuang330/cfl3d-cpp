// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// ron_common.h — shared globals + helpers for the C++ translation of CFL3D
// `ronnie` (patched-grid preprocessor). Sequential build (no MPI).
//
// Conventions (MUST be followed by every translated .cpp):
//  * Fortran default REAL is 8 bytes here (compiled with -fdefault-real-8),
//    so every Fortran `real`/`double precision` -> C++ `double`.
//  * Fortran arrays are column-major and 1-based. We keep 1-based indexing by
//    over-allocating and using helper accessors, OR by explicit (idx-1)
//    arithmetic that mirrors the Fortran storage formula. Multi-dim Fortran
//    array A(d1,d2,...) at (i,j,..) maps to flat offset
//        (i-1) + (j-1)*d1 + (k-1)*d1*d2 + ...
//  * COMMON blocks are single global struct instances (definitions live in
//    ron_common.cpp; declared extern here).
//
#pragma once
#include "../../common/tools_common.h"
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

// ---- Fortran intrinsics as inline helpers (double) -------------------------
inline double f_abs(double x){ return std::fabs(x); }
inline int    f_iabs(int x){ return std::abs(x); }
inline double f_sign(double a,double b){ return (b>=0.0)? std::fabs(a): -std::fabs(a); }
inline int    i_sign(int a,int b){ return (b>=0)? std::abs(a): -std::abs(a); }
inline double f_dim(double a,double b){ double d=a-b; return d>0.0?d:0.0; }
inline double f_min(double a,double b){ return a<b?a:b; }
inline double f_max(double a,double b){ return a>b?a:b; }
inline int    i_min(int a,int b){ return a<b?a:b; }
inline int    i_max(int a,int b){ return a>b?a:b; }
inline double f_float(int i){ return (double)i; }
inline int    f_int(double x){ return (int)x; }       // Fortran int() truncates toward 0
inline int    f_nint(double x){ return (int)std::llround(x); }
inline double ccabs(double x){ return std::fabs(x); } // ccomplex real-build stub

// ---- Fortran character string field (fixed-width, space-padded) ------------
// character*n -> std::string kept exactly n chars where it matters for output.
struct FStr {
    std::string s;
    int len;
    explicit FStr(int n=80): s(std::string((size_t)n,' ')), len(n) {}
    void assign(const std::string& v){ s.assign((size_t)len,' ');
        for(int i=0;i<len && i<(int)v.size();++i) s[i]=v[i]; }
    const char* c() const { return s.c_str(); }
};

// ---- COMMON blocks (declared here, defined in ron_common.cpp) ---------------
// /ginfo/ — block geometry + 1-d storage indices (only jdim,kdim,idim,jj2,kk2,ii2,lx,ly,lz used)
struct Ginfo {
    int jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie;
    int lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis;
    int lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg;
    int ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci;
    int lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk;
    int lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0;
    int lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti;
    int lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg;
};
extern Ginfo ginfo_;

struct Info { // /info/
    double title[20]; double rkap[3]; double xmach,alpha,beta,dt,fmax;
    int nit,ntt,idiag[3],nitfo,iflagts,iflim[3],nres,levelb[5],mgflag;
    int iconsf,mseq,ncyc1[5],levelt[5],nitfo1[5],ngam,nsm[5],iipv;
};
extern Info info_;

struct Params { // /params/
    int lmaxgr,lmaxbl,lmxseg,lmaxcs,lnplts,lmxbli,lmaxxe,lnsub1,lintmx,lmxxe;
    int liitot,isum,lncycm,isum_n,lminnode,isumi,isumi_n,lmptch,lmsub1,lintmax;
    int libufdim,lnbuf,llbcprd,llbcemb,llbcrad,lnmds,lmaxaes,lnslave,lmxsegdg,lnmaster;
};
extern Params params_;

struct Filenam { // /filenam/  (24 names, character*80)
    FStr grid,plt3dg,plt3dq,output,residual,turbres,blomx,output2,printout,
         pplunge,ovrlap,patch,restrt,subres,subtur,grdmov,alphahist,errfile,
         preout,aeinp,aeout,sdhist,avgg,avgq;
    Filenam(): grid(80),plt3dg(80),plt3dq(80),output(80),residual(80),turbres(80),
        blomx(80),output2(80),printout(80),pplunge(80),ovrlap(80),patch(80),
        restrt(80),subres(80),subtur(80),grdmov(80),alphahist(80),errfile(80),
        preout(80),aeinp(80),aeout(80),sdhist(80),avgg(80),avgq(80) {}
};
extern Filenam filenam_;

struct Mydist2 { int nnodes,myhost,myid,mycomm; }; extern Mydist2 mydist2_;
struct Igrdtyp { int ip3dgrd,ialph; };             extern Igrdtyp igrdtyp_;
struct Conversion { double radtodeg; };            extern Conversion conversion_;
struct Unit5 { int iunit5; };                      extern Unit5 unit5_;
struct Tracer { int itrace; };                     extern Tracer tracer_;
struct Zero { int iexp; };                         extern Zero zero_;
struct Sklt1 { int isklt1; };                      extern Sklt1 sklt1_;
struct IsDpatch { int maxdcnt; };                  extern IsDpatch is_dpatch_;
struct Save { int locv,ireq; };                    extern Save save_;
struct Tacos { int iretry; };                      extern Tacos tacos_;
struct Areas { double ap[3]; int imaxa; };         extern Areas areas_;
struct Hist { double histry[4]; };                 extern Hist hist_;
struct Degshf { int ideg[3]; };                    extern Degshf degshf_;
struct Fvfds { double rkap0[3]; int ifds[3]; };    extern Fvfds fvfds_;
struct Ginfo2 { int lq2avg,iskip_blocks,inc_2d[3],inc_coarse[3]; }; extern Ginfo2 ginfo2_;
struct Lam { int ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi,i_lam_forcezero; }; extern Lam lam_;
struct Reyue { double reue,tinf; int ivisc[3]; };  extern Reyue reyue_;
struct Wallfun { int iwf[3]; };                    extern Wallfun wallfun_;
// /tol/ has POSITIONALLY-aliased names across files:
//   collapse/rechk: epsc,epsc0,epsreen,eps    patcher: epsc,epsc0,epsreen,epscoll
//   invert(some):   epsc,epsc0,eps,epscoll  -> we expose all 4 slots by position.
struct Tol { double epsc,epsc0,epsreen,eps; }; extern Tol tol_;
// NOTE: where a file names the 3rd/4th slots differently, map by POSITION:
//   slot0=epsc, slot1=epsc0, slot2=(epsreen or eps), slot3=(eps or epscoll).
struct Time1 { double tim,tm; };                   extern Time1 time1_;
struct Time2 { double tim_nd; };                   extern Time2 time2_;
struct Rjbdbgi { int lunfio0; };                   extern Rjbdbgi rjbdbgi_;

// ---- global Fortran unit -> FILE* map (units used: 1,2,5,7,9,10,11,25,66,99)
extern FILE* g_units[128];
inline FILE*& FUNIT(int u){ return g_units[u]; }

// ---- forward declarations of all translated subroutines --------------------
// (signatures mirror the Fortran argument lists; arrays passed as double*/int*)
void umalloc(long numbytes,int itype,const char* varname,long& memuse,int stats);
void cputim(int istart,int nnodes,FStr& string,int myhost,int myid,int mycomm,int iunit);
void termn8(int myid_stop,int ierrflg,int ibufdim,int nbuf,FStr* bou,int* nou);
void my_flush(int iunit);
void outbuf(const FStr& line,int iunit);
void writ_buf(int nbl,int iunit,int* nou,FStr* bou,int nbuf,int ibufdim,
              int myhost,int myid,int mycomm,int* mblk2nd,int maxbl);
void lead(int nbl,int* lw,int* lw2,int maxbl);
void collx(double* x,double* y,double* z,double* xc,double* yc,double* zc,
           int jdim,int kdim,int idim,int jj2,int kk2,int ii2);
void rp3d(double* x,double* y,double* z,int jdim,int kdim,int idim,int igrid,
          int ialph,int idum,int& irr);
void loadgr(double* w,int mgwk,int lx,int ly,int lz,int jindex,
            double* x,double* y,double* z,int mdim,int ndim,
            int idimg,int jdimg,int kdimg);

// ---- diagnostic buffer (character*120 bou(ibufdim,nbuf)) helper -------------
// The Fortran pattern
//     nou(iu) = min(nou(iu)+1,ibufdim)
//     write(bou(nou(iu),iu),FMT) ...
// appends one 120-char line to buffer column `iu` (1-based). bou is passed as a
// raw FStr* base pointer; element (n,iu) (1-based) lives at (n-1)+(iu-1)*ibufdim.
// `bou_put` reproduces one such formatted line. Pass the already-formatted text
// (built with snprintf reproducing the FORTRAN FORMAT); it is space-padded to
// width 120 like a character*120 record. Returns nothing.
inline void bou_put(int* nou,FStr* bou,int ibufdim,int iu,const std::string& txt){
    nou[iu-1] = i_min(nou[iu-1]+1, ibufdim);
    int n = nou[iu-1];                       // 1-based line index
    FStr& dst = bou[(n-1) + (iu-1)*ibufdim];
    dst = FStr(120);
    dst.assign(txt);
}

// ---- Group B (patcher orchestration + grid-geometry prep) ------------------
void dsmin(int jdim,int kdim,int nsub,int* jjmax,int* kkmax,int lmax,
           double* x,double* y,double* z,double xc,double yc,double zc,
           int& j1,int& k1,int& l1,int* lout,int lflag,
           int* xif1,int* xif2,int* etf1,int* etf2);
void avgint(double* xiec,double* etac,int* nblkc,int nptc,
            double* xief,double* etaf,int* nblkf,int nptf,
            int j1c,int j2c,int k1c,int k2c,int j1f,int j2f,int k1f,int k2f);
void transp(int mdim,int ndim,int jmax1,int kmax1,int msub1,int l,
            double* x1,double* y1,double* z1,double* dx,double* dy,double* dz,
            int intmx,int intt);
void rotatp(int mdim,int ndim,int jmax1,int kmax1,int msub1,int l,
            double* x1,double* y1,double* z1,
            double* dthetx,double* dthety,double* dthetz,
            double* xorig,double* yorig,double* zorig,
            int mbl,int nn,int intmx,int intt);
void expand(int mdim,int ndim,int msub1,int jmax1,int kmax1,int l,
            double* x1,double* y1,double* z1,double* xte,double* yte,double* zte,
            double factjlo,double factjhi,double factklo,double factkhi,
            int jmax2,int kmax2,double* x2,double* y2,double* z2);
void collapse(int jdim,int kdim,int jmaxo,int kmaxo,
              double* xo,double* yo,double* zo,
              int* nou,FStr* bou,int nbuf,int ibufdim);
void rechk(int mdim1,int ndim1,int* jimage,int* kimage,int msub1,
           int jmax1,int kmax1,int l,double* x1,double* y1,double* z1,
           int xie1,int xie2,int eta1,int eta2,
           int* nou,FStr* bou,int nbuf,int ibufdim,int myid,int* mblk2nd,int maxbl);
void arc(int jdim1,int kdim1,int msub1,int msub2,
         int* jjmax1,int* kkmax1,int lmax1,double* x1,double* y1,double* z1,
         int limit0,int* jjmax2,int* kkmax2,
         double* x2,double* y2,double* z2,double* xie2,double* eta2,int* mblkpt,
         int ifit,int itmax,int jcorr,int kcorr,
         double* sxie,double* seta,double* sxie2,double* seta2,
         double* xie2s,double* eta2s,
         int intmx,int icheck,int* nblkj,int* nblkk,int* jmm,int* kmm,
         int j21,int j22,int k21,int k22,
         int npt,int* xif1,int* xif2,int* etf1,int* etf2,
         int* nou,FStr* bou,int nbuf,int ibufdim,int* mblk2nd,int maxbl);
void topol2(int jdim1,int kdim1,int msub1,int* jjmax1,int* kkmax1,int lmax1,
            double xie,double eta,double* s,int limit0,double sc,
            int jp,int kp,int lp,int lsrch,int itmax,double& xiet,double& etat,
            int* xif1,int* xif2,int* etf1,int* etf2,
            int* nou,FStr* bou,int nbuf,int ibufdim,int myid,int* mblk2nd,int maxbl);
void invert(int mdim,int ndim,int msub1,int msub2,int* jte,int* kte,
            int lmax1,double* xte,double* yte,double* zte,
            double* xmi,double* ymi,double* zmi,double* xmie,double* ymie,double* zmie,
            int limit0,int* jjmax2,int* kkmax2,double* x2,double* y2,double* z2,
            double* xie2,double* eta2,int* mblkpt,double* temp,int* jimage,int* kimage,
            int ifit,int itmax,double* sxie,double* seta,double* sxie2,double* seta2,
            double* xie2s,double* eta2s,int intmx,int icheck,int* nblk1,int* nblk2,
            int* jmm,int* kmm,int mcxie,int mceta,int* lout,
            int xi1,int xi2,int et1,int et2,int npt,int ic0,int iorph,int itoss0,
            int ncall,int ioutpt,int* xif1,int* xif2,int* etf1,int* etf2,
            int iself,int ifiner,double* xie2f,double* eta2f,int* mblkptf,int nptf,
            int xi1f,int xi2f,int et1f,int et2f,int iavg,
            int* nou,FStr* bou,int nbuf,int ibufdim,int myid,int* mblk2nd,int maxbl);
void patcher(int nbl,int* lw,double* w,int mgwk,double* wk,int nwork,int ncall,
             int ioutpt,int it_thro,int maxbl,int msub1,int intmx,int mxxe,int mptch,
             int* jdimg,int* kdimg,int* idimg,double* windx,int nintr,int* iindx,
             int* llimit,int* iitmax,int* mmcxie,int* mmceta,int* ncheck,int* iifit,
             int* mblkpt,int* iic0,int* iiorph,int* iitoss,int* ifiner,
             double* factjlo,double* factjhi,double* factklo,double* factkhi,
             double* dx,double* dy,double* dz,double* dthetx,double* dthety,double* dthetz,
             int* isav_dpat,int* isav_dpat_b,
             double* xte,double* yte,double* zte,double* xmi,double* ymi,double* zmi,
             double* xmie,double* ymie,double* zmie,int* jjmax1,int* kkmax1,
             int* jimage,int* kimage,double* xorig,double* yorig,double* zorig,
             int* jte,int* kte,double* sxie,double* seta,double* sxie2,double* seta2,
             double* xie2s,double* eta2s,double* temp,double* x2,double* y2,double* z2,
             int* nblk1,int* nblk2,int* jmm,int* kmm,double* x1,double* y1,double* z1,
             int* lout,int* xif1,int* xif2,int* etf1,int* etf2,int* ireq_ar,
             int myid,int myhost,int mycomm,int* mblk2nd,int* nou,FStr* bou,
             int nbuf,int ibufdim);

void setup(int* lw,int* lw2,double* w,int mgwk,int* itest,int* jtest,int* ktest,
           int maxbl,int maxgr,int ngrid,int* ncgg,int* iemg,int* nblg);

void usrint(int sig);

void global2(int maxbl,int maxgr,int msub1,int& nintr,int intmx,int ngrid,
             int* idimg,int* jdimg,int* kdimg,int* levelg,int* ncgg,int* nblg,
             int* iindx,int* llimit,int* iitmax,int* mmcxie,int* mmceta,
             int* ncheck,int* iifit,int* iic0,int* iiorph,int* iitoss,int* ifiner,
             double* dx,double* dy,double* dz,double* dthetx,double* dthety,
             double* dthetz,int myid,int mptch,int mxxe,int icall,int iunit25,
             int* nou,FStr* bou,int ibufdim,int nbuf,int* ifrom,int* xif1,int* etf1,
             int* xif2,int* etf2,int* igridg,int* iemg,int nblock,int ioflag,int imode);

void pre_patch(int nbl,int* lw,int& icount,int ninter,int* iindex,int intmax,
               int nsub1,int* isav_pat,int* isav_pat_b,int* jjmax1,int* kkmax1,
               int* iiint1,int* iiint2,int maxbl,int* jdimg,int* kdimg,int* idimg,
               int ierrflg);

void sizer(int& maxgr,int& maxbl,int& maxxe,int& nsub1,int& intmax,int& mpatch,
           int& nwork,int nbuf0,int ibufdim0,int maxgr0,int maxbl0,int maxxe0,
           int nsub10,int intmax0,int mpatch0,int imode);

void ronnie(int maxbl,int maxgr,int nwork,int intmax,int nsub1,int maxxe,
            int mpatch,int nbuf,int ibufdim);
// (patcher is declared above with the Group B math routines.)

// ---- Fortran formatted-write helpers (defined in ron_format.cpp) -----------
// Emulate specific FORMAT statements byte-exactly. All write to FUNIT(iunit).
std::string f_iw(long v,int w);              // Fortran Iw right-justified integer
std::string f_ild(long v);                   // Fortran list-directed integer (write(u,*))
std::string f_glist(double v);               // Fortran list-directed real (write(u,*))
std::string f_fw(double v,int w,int d);      // Fortran Fw.d fixed-point
void ron_banner(int iunit,int isron,double memuse); // preronnie/ronnie banner+mem
void ron_iofiles(int iunit,const FStr& grid,const FStr& output,const FStr& patch);
