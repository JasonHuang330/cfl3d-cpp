// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// splitter.cpp - translated from splitter.F
#include "splitter.h"
#include "runtime/fortran_io.h"
#include "cfl3d_to_tlns3d.h"
#include "tlns3d_to_cfl3d.h"
#include "termn8.h"
#include "umalloc.h"
#include "mapin.h"
#include "mapout.h"
#include "mapblk.h"
#include "mapsplt.h"
#include "ingrd.h"
#include "outgrd.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

using namespace cfl3d_to_tlns3d_ns;
using namespace tlns3d_to_cfl3d_ns;
using namespace termn8_ns;
using namespace umalloc_ns;
using namespace mapin_ns;
using namespace mapout_ns;
using namespace mapblk_ns;
using namespace mapsplt_ns;
using namespace ingrd_ns;
using namespace outgrd_ns;

namespace splitter_ns {

void splitter(int& mbloc, int& maxseg, int& npmax, int& mxbli, int& intmax, int& nsub1, int& msegt, int& msegn, int& msplt, int& mtot, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    // Common block references
    int& iunit5  = cmn_unit5.iunit5;
    int& nbloc   = cmn_grdinfo.nbloc;
    int& ncgmax  = cmn_chk1.ncgmax;
    int& ip3dgrd = cmn_info3.ip3dgrd;
    int& ndv     = cmn_info3.ndv;
    char (&gridin)[80]    = cmn_cflfiles.gridin;
    char (&gridout)[80]   = cmn_cflfiles.gridout;
    char (&sdgridin)[80]  = cmn_cflfiles.sdgridin;
    char (&sdgridout)[80] = cmn_cflfiles.sdgridout;
    char (&dovrlap)[80]   = cmn_cflfiles.dovrlap;
    char (&dpatch)[80]    = cmn_cflfiles.dpatch;
    char (&dresid)[80]    = cmn_cflfiles.dresid;

    // Local character variables
    char header[80];
    char imapin[80];
    char imapout[80];
    char cflinp[80];
    char cflout[80];
    char roninp[80];
    char ronout[80];

    // Local scalar integers
    int memuse, stats = 0;
    int icflin, ironin, ibin, isdin;
    int icflout, ironout, ibout, isdout;
    int iparin, iparout;
    int nbl, iflag, ntot, nout;
    int ipatch0, isd;
    int nblk, nblsd;
    int nnodes;
    int nrotat, ntrans;
    double tlref, rlref;
    int ioflag;
    int iflg, ncgmax1;

    // Local scalar reals
    double xdum0, ydum0, zdum0;

    // Allocatable arrays - integer
    FortranArray4D<int>    bcval_i;   // not used as int - actually double below
    FortranArray4D<double> bcval;
    FortranArray4D<double> bcvali;
    FortranArray4D<double> bcvalj;
    FortranArray4D<double> bcvalk;
    FortranArray2D<double> cpchamber;
    FortranArray2D<double> cporous;
    FortranArray2D<double> cq;
    FortranArray2D<double> cradiation;
    FortranArray2D<int>    etf1;
    FortranArray2D<int>    etf2;
    FortranArray1D<int>    ett1;
    FortranArray1D<int>    ett2;
    FortranArray1D<int>    iadvance;
    FortranArray4D<int>    ibcinfo;
    FortranArray3D<int>    ibct;
    FortranArray1D<int>    ibeg;
    FortranArray1D<int>    ibif1;
    FortranArray1D<int>    ibif2;
    FortranArray1D<int>    iconcat;
    FortranArray1D<int>    idbloc;
    FortranArray2D<int>    idegg;
    FortranArray2D<int>    idiagg;
    FortranArray1D<int>    idimg;
    FortranArray1D<int>    idmold;
    FortranArray1D<int>    idnext;
    FortranArray1D<int>    idno;
    FortranArray1D<int>    idobl;
    FortranArray1D<int>    idseg;
    FortranArray1D<int>    ieg;
    FortranArray1D<int>    iemg;
    FortranArray2D<int>    ifdsg;
    FortranArray2D<int>    iflimg;
    FortranArray3D<int>    ifoflg;
    FortranArray1D<int>    iforce;
    FortranArray2D<int>    ifrom;
    FortranArray1D<int>    ifsor;
    FortranArray1D<int>    igridg;
    FortranArray1D<int>    iic0;
    FortranArray1D<int>    iifit;
    FortranArray1D<int>    iiorph;
    FortranArray1D<int>    iitmax;
    FortranArray2D<int>    ijk;
    FortranArray1D<int>    il;
    FortranArray1D<int>    ilamhig;
    FortranArray1D<int>    ilamlog;
    FortranArray1D<int>    iln;
    FortranArray1D<int>    ilo;
    FortranArray1D<int>    ilosd;
    FortranArray3D<int>    imap;
    FortranArray1D<int>    imax;
    FortranArray1D<int>    imin;
    FortranArray1D<int>    inewgg;
    FortranArray1D<int>    intrfc;
    FortranArray2D<int>    iold;
    FortranArray1D<int>    iovrlp;
    FortranArray1D<int>    ipatch;
    FortranArray2D<int>    iredundant;
    FortranArray1D<int>    isblocin;
    FortranArray1D<int>    isg;
    FortranArray1D<int>    isin;
    FortranArray1D<int>    isin0;
    FortranArray1D<int>    iskp;
    FortranArray3D<int>    isva;
    FortranArray1D<int>    ito;
    FortranArray2D<int>    itrb;
    FortranArray1D<int>    itrb1;
    FortranArray1D<int>    itrb2;
    FortranArray1D<int>    iturbb;
    FortranArray2D<int>    ivisb;
    FortranArray2D<int>    iviscg;
    FortranArray4D<int>    jbcinfo;
    FortranArray1D<int>    jdimg;
    FortranArray1D<int>    jeg;
    FortranArray1D<int>    jl;
    FortranArray1D<int>    jlamhig;
    FortranArray1D<int>    jlamlog;
    FortranArray1D<int>    jln;
    FortranArray1D<int>    jlo;
    FortranArray1D<int>    jlosd;
    FortranArray1D<int>    jmax;
    FortranArray1D<int>    jmin;
    FortranArray1D<int>    jsg;
    FortranArray1D<int>    jskp;
    FortranArray1D<int>    jtrb1;
    FortranArray1D<int>    jtrb2;
    FortranArray4D<int>    kbcinfo;
    FortranArray1D<int>    kdimg;
    FortranArray1D<int>    keg;
    FortranArray1D<int>    kl;
    FortranArray1D<int>    klamhig;
    FortranArray1D<int>    klamlog;
    FortranArray1D<int>    kln;
    FortranArray1D<int>    klo;
    FortranArray1D<int>    klosd;
    FortranArray1D<int>    kmax;
    FortranArray1D<int>    kmin;
    FortranArray1D<int>    ksg;
    FortranArray1D<int>    kskp;
    FortranArray1D<int>    ktrb1;
    FortranArray1D<int>    ktrb2;
    FortranArray1D<int>    levelg;
    FortranArray3D<int>    limblk;
    FortranArray1D<int>    llimit;
    FortranArray2D<int>    mit;
    FortranArray1D<int>    mmceta;
    FortranArray1D<int>    mmcxie;
    FortranArray3D<int>    nb1;
    FortranArray1D<int>    nb1s;
    FortranArray3D<int>    nb2;
    FortranArray1D<int>    nb2s;
    FortranArray1D<int>    nbci0;
    FortranArray1D<int>    nbcidim;
    FortranArray1D<int>    nbcj0;
    FortranArray1D<int>    nbcjdim;
    FortranArray1D<int>    nbck0;
    FortranArray1D<int>    nbckdim;
    FortranArray1D<int>    nbi1;
    FortranArray1D<int>    nbi2;
    FortranArray1D<int>    nbj1;
    FortranArray1D<int>    nbj2;
    FortranArray1D<int>    nbk1;
    FortranArray1D<int>    nbk2;
    FortranArray1D<int>    nblcg;
    FortranArray1D<int>    nblg;
    FortranArray2D<int>    nblkk;
    FortranArray1D<int>    nblon;
    FortranArray1D<int>    nbs;
    FortranArray1D<int>    ncgg;
    FortranArray1D<int>    nd11;
    FortranArray1D<int>    nd12;
    FortranArray1D<int>    nd21;
    FortranArray1D<int>    nd22;
    FortranArray3D<int>    ndat;
    FortranArray1D<int>    ndirin;
    FortranArray1D<int>    ndirin0;
    FortranArray3D<int>    ne1;
    FortranArray1D<int>    ne1s;
    FortranArray3D<int>    ne2;
    FortranArray1D<int>    ne2s;
    FortranArray1D<int>    nei1;
    FortranArray1D<int>    nei2;
    FortranArray1D<int>    nej1;
    FortranArray1D<int>    nej2;
    FortranArray1D<int>    nek1;
    FortranArray1D<int>    nek2;
    FortranArray1D<int>    nfb;
    FortranArray1D<int>    ni;
    FortranArray1D<int>    nj;
    FortranArray1D<int>    nk;
    FortranArray1D<int>    npts;
    FortranArray1D<int>    nseg;
    FortranArray1D<int>    nseg1;
    FortranArray1D<int>    nsihi;
    FortranArray1D<int>    nsilo;
    FortranArray1D<int>    nsjhi;
    FortranArray1D<int>    nsjlo;
    FortranArray1D<int>    nskhi;
    FortranArray1D<int>    nsklo;
    FortranArray1D<int>    nsubbl;
    FortranArray1D<int>    nxtseg;
    FortranArray2D<double> rkap0g;
    FortranArray2D<double> twotref;
    FortranArray1D<double> val;
    FortranArray1D<double> x;
    FortranArray1D<double> xdmold;
    FortranArray2D<double> xdum;
    FortranArray2D<int>    xif1;
    FortranArray2D<int>    xif2;
    FortranArray1D<int>    xit1;
    FortranArray1D<int>    xit2;
    FortranArray3D<double> xmap;
    FortranArray1D<double> y;
    FortranArray1D<double> z;



    // allocate memory
    memuse = 0;

    // Helper lambda for umalloc calls
    auto do_umalloc = [&](int n, int iflag, const char* txt) {
        char tbuf[1];
        tbuf[0] = txt[0];
        stats = 0;  // allocate() succeeded, so stats=0 (matches Fortran allocate stat=stats)
        umalloc(n, iflag, tbuf, memuse, stats);
    };

    // Helper: check if a char[80] string equals "null" (Fortran string comparison)
    auto is_null_str = [](const char* s) -> bool {
        return std::strncmp(s, "null", 4) == 0 && (s[4] == '\0' || s[4] == ' ');
    };

    bcval.allocate(mbloc, 6, msegn, 7);
    do_umalloc(mbloc*6*msegn*7, 0, "bcval");
    bcvali.allocate(mbloc, maxseg, 7, 2);
    do_umalloc(mbloc*maxseg*7*2, 0, "bcvali");
    bcvalj.allocate(mbloc, maxseg, 7, 2);
    do_umalloc(mbloc*maxseg*7*2, 0, "bcvalj");
    bcvalk.allocate(mbloc, maxseg, 7, 2);
    do_umalloc(mbloc*maxseg*7*2, 0, "bcvalk");
    cpchamber.allocate(msegn, mbloc);
    do_umalloc(msegn*mbloc, 0, "cpchamber");
    cporous.allocate(msegn, mbloc);
    do_umalloc(msegn*mbloc, 0, "cporous");
    cq.allocate(msegn, mbloc);
    do_umalloc(msegn*mbloc, 0, "cq");
    cradiation.allocate(msegn, mbloc);
    do_umalloc(msegn*mbloc, 0, "cradiation");
    etf1.allocate(intmax, nsub1);
    do_umalloc(intmax*nsub1, 1, "etf1");
    etf2.allocate(intmax, nsub1);
    do_umalloc(intmax*nsub1, 1, "etf2");
    ett1.allocate(intmax);
    do_umalloc(intmax, 1, "ett1");
    ett2.allocate(intmax);
    do_umalloc(intmax, 1, "ett2");
    iadvance.allocate(mbloc);
    do_umalloc(mbloc, 1, "iadvance");
    ibcinfo.allocate(mbloc, maxseg, 7, 2);
    do_umalloc(mbloc*maxseg*7*2, 1, "ibcinfo");
    ibct.allocate(mbloc, 6, msegn);
    do_umalloc(mbloc*6*msegn, 1, "ibct");
    ibeg.allocate(mbloc);
    do_umalloc(mbloc, 1, "ibeg");
    ibif1.allocate(mxbli);
    do_umalloc(mxbli, 1, "ibif1");
    ibif2.allocate(mxbli);
    do_umalloc(mxbli, 1, "ibif2");
    iconcat.allocate(intmax);
    do_umalloc(intmax, 1, "iconcat");
    idbloc.allocate(mtot);
    do_umalloc(mtot, 1, "idbloc");
    idegg.allocate(mbloc, 3);
    do_umalloc(mbloc*3, 1, "idegg");
    idiagg.allocate(mbloc, 3);
    do_umalloc(mbloc*3, 1, "idiagg");
    idimg.allocate(mbloc);
    do_umalloc(mbloc, 1, "idimg");
    idmold.allocate(msegt);
    do_umalloc(msegt, 1, "idmold");
    idnext.allocate(mbloc);
    do_umalloc(mbloc, 1, "idnext");
    idno.allocate(mbloc);
    do_umalloc(mbloc, 1, "dno");
    idobl.allocate(mbloc);
    do_umalloc(mbloc, 1, "idobl");
    idseg.allocate(mbloc);
    do_umalloc(mbloc, 1, "idseg");
    ieg.allocate(mbloc);
    do_umalloc(mbloc, 1, "ieg");
    iemg.allocate(mbloc);
    do_umalloc(mbloc, 1, "iemg");
    ifdsg.allocate(mbloc, 3);
    do_umalloc(mbloc*3, 1, "ifdsg");
    iflimg.allocate(mbloc, 3);
    do_umalloc(mbloc*3, 1, "iflimg");
    ifoflg.allocate(mbloc, 6, msegn);
    do_umalloc(mbloc*6*msegn, 1, "ifoflg");
    iforce.allocate(mbloc);
    do_umalloc(mbloc, 1, "iforce");
    ifrom.allocate(intmax, nsub1);
    do_umalloc(intmax*nsub1, 1, "ifrom");
    ifsor.allocate(mbloc);
    do_umalloc(mbloc, 1, "ifsor");
    igridg.allocate(mbloc);
    do_umalloc(mbloc, 1, "igridg");
    iic0.allocate(intmax);
    do_umalloc(intmax, 1, "iic0");
    iifit.allocate(intmax);
    do_umalloc(intmax, 1, "iifit");
    iiorph.allocate(intmax);
    do_umalloc(intmax, 1, "iiorph");
    iitmax.allocate(intmax);
    do_umalloc(intmax, 1, "iitmax");
    ijk.allocate(6, mbloc);
    do_umalloc(6*mbloc, 1, "ijk");
    il.allocate(mbloc);
    do_umalloc(mbloc, 1, "il");
    ilamhig.allocate(mbloc);
    do_umalloc(mbloc, 1, "ilamhig");
    ilamlog.allocate(mbloc);
    do_umalloc(mbloc, 1, "ilamlog");
    iln.allocate(mbloc);
    do_umalloc(mbloc, 1, "iln");
    ilo.allocate(mbloc);
    do_umalloc(mbloc, 1, "ilo");
    ilosd.allocate(mbloc);
    do_umalloc(mbloc, 1, "ilosd");
    imap.allocate(msegt, msegn, mbloc);
    do_umalloc(msegt*msegn*mbloc, 1, "imap");
    imax.allocate(mbloc);
    do_umalloc(mbloc, 1, "imax");
    imin.allocate(mbloc);
    do_umalloc(mbloc, 1, "imin");
    inewgg.allocate(mbloc);
    do_umalloc(mbloc, 1, "inewgg");
    intrfc.allocate(mtot);
    do_umalloc(mtot, 1, "intrfc");
    iold.allocate(4, mtot);
    do_umalloc(4*mtot, 1, "iold");
    iovrlp.allocate(mbloc);
    do_umalloc(mbloc, 1, "iovrlp");
    ipatch.allocate(mtot);
    do_umalloc(mtot, 1, "ipatch");
    iredundant.allocate(mbloc, msegn);
    do_umalloc(mbloc*msegn, 1, "iredundant");
    isblocin.allocate(msplt*mbloc);
    do_umalloc(msplt*mbloc, 1, "isblocin");
    isg.allocate(mbloc);
    do_umalloc(mbloc, 1, "isg");
    isin.allocate(msplt*mbloc);
    do_umalloc(msplt*mbloc, 1, "isin");
    isin0.allocate(msplt);
    do_umalloc(msplt, 1, "isin0");
    iskp.allocate(mbloc);
    do_umalloc(mbloc, 1, "iskp");
    isva.allocate(2, 2, mxbli);
    do_umalloc(2*2*mxbli, 1, "isva");
    ito.allocate(intmax);
    do_umalloc(intmax, 1, "ito");
    itrb.allocate(7, mbloc);
    do_umalloc(7*mbloc, 1, "itrb");
    itrb1.allocate(mbloc);
    do_umalloc(mbloc, 1, "itrb1");
    itrb2.allocate(mbloc);
    do_umalloc(mbloc, 1, "itrb2");
    iturbb.allocate(mbloc);
    do_umalloc(mbloc, 1, "iturbb");
    ivisb.allocate(msegt, mbloc);
    do_umalloc(msegt*mbloc, 1, "ivisb");
    iviscg.allocate(mbloc, 3);
    do_umalloc(mbloc*3, 1, "iviscg");
    jbcinfo.allocate(mbloc, maxseg, 7, 2);
    do_umalloc(mbloc*maxseg*7*2, 1, "jbcinfo");
    jdimg.allocate(mbloc);
    do_umalloc(mbloc, 1, "jdimg");
    jeg.allocate(mbloc);
    do_umalloc(mbloc, 1, "jeg");
    jl.allocate(mbloc);
    do_umalloc(mbloc, 1, "jl");
    jlamhig.allocate(mbloc);
    do_umalloc(mbloc, 1, "jlamhig");
    jlamlog.allocate(mbloc);
    do_umalloc(mbloc, 1, "jlamlog");
    jln.allocate(mbloc);
    do_umalloc(mbloc, 1, "jln");
    jlo.allocate(mbloc);
    do_umalloc(mbloc, 1, "jlo");
    jlosd.allocate(mbloc);
    do_umalloc(mbloc, 1, "jlosd");
    jmax.allocate(mbloc);
    do_umalloc(mbloc, 1, "jmax");
    jmin.allocate(mbloc);
    do_umalloc(mbloc, 1, "jmin");
    jsg.allocate(mbloc);
    do_umalloc(mbloc, 1, "jsg");
    jskp.allocate(mbloc);
    do_umalloc(mbloc, 1, "jskp");
    jtrb1.allocate(mbloc);
    do_umalloc(mbloc, 1, "jtrb1");
    jtrb2.allocate(mbloc);
    do_umalloc(mbloc, 1, "jtrb2");
    kbcinfo.allocate(mbloc, maxseg, 7, 2);
    do_umalloc(mbloc*maxseg*7*2, 1, "kbcinfo");
    kdimg.allocate(mbloc);
    do_umalloc(mbloc, 1, "kdimg");
    keg.allocate(mbloc);
    do_umalloc(mbloc, 1, "keg");
    kl.allocate(mbloc);
    do_umalloc(mbloc, 1, "kl");
    klamhig.allocate(mbloc);
    do_umalloc(mbloc, 1, "klamhig");
    klamlog.allocate(mbloc);
    do_umalloc(mbloc, 1, "klamlog");
    kln.allocate(mbloc);
    do_umalloc(mbloc, 1, "kln");
    klo.allocate(mbloc);
    do_umalloc(mbloc, 1, "klo");
    klosd.allocate(mbloc);
    do_umalloc(mbloc, 1, "klosd");
    kmax.allocate(mbloc);
    do_umalloc(mbloc, 1, "kmax");
    kmin.allocate(mbloc);
    do_umalloc(mbloc, 1, "kmin");
    ksg.allocate(mbloc);
    do_umalloc(mbloc, 1, "ksg");
    kskp.allocate(mbloc);
    do_umalloc(mbloc, 1, "kskp");
    ktrb1.allocate(mbloc);
    do_umalloc(mbloc, 1, "ktrb1");
    ktrb2.allocate(mbloc);
    do_umalloc(mbloc, 1, "ktrb2");
    levelg.allocate(mbloc);
    do_umalloc(mbloc, 1, "levelg");
    limblk.allocate(2, 6, mxbli);
    do_umalloc(2*6*mxbli, 1, "limblk");
    llimit.allocate(intmax);
    do_umalloc(intmax, 1, "llimit");
    mit.allocate(5, mbloc);
    do_umalloc(5*mbloc, 1, "mit");
    mmceta.allocate(intmax);
    do_umalloc(intmax, 1, "mmceta");
    mmcxie.allocate(intmax);
    do_umalloc(intmax, 1, "mmcxie");

    nb1.allocate(mbloc, 6, msegn);
    do_umalloc(mbloc*6*msegn, 1, "nb1");
    nb1s.allocate(mbloc);
    do_umalloc(mbloc, 1, "nb1s");
    nb2.allocate(mbloc, 6, msegn);
    do_umalloc(mbloc*6*msegn, 1, "nb2");
    nb2s.allocate(mbloc);
    do_umalloc(mbloc, 1, "nb2s");
    nbci0.allocate(mbloc);
    do_umalloc(mbloc, 1, "nbci0");
    nbcidim.allocate(mbloc);
    do_umalloc(mbloc, 1, "nbcidim");
    nbcj0.allocate(mbloc);
    do_umalloc(mbloc, 1, "nbcj0");
    nbcjdim.allocate(mbloc);
    do_umalloc(mbloc, 1, "nbcjdim");
    nbck0.allocate(mbloc);
    do_umalloc(mbloc, 1, "nbck0");
    nbckdim.allocate(mbloc);
    do_umalloc(mbloc, 1, "nbckdim");
    nbi1.allocate(mxbli);
    do_umalloc(mxbli, 1, "nbi1");
    nbi2.allocate(mxbli);
    do_umalloc(mxbli, 1, "nbi2");
    nbj1.allocate(mxbli);
    do_umalloc(mxbli, 1, "nbj1");
    nbj2.allocate(mxbli);
    do_umalloc(mxbli, 1, "nbj2");
    nbk1.allocate(mxbli);
    do_umalloc(mxbli, 1, "nbk1");
    nbk2.allocate(mxbli);
    do_umalloc(mxbli, 1, "nbk2");
    nblcg.allocate(mbloc);
    do_umalloc(mbloc, 1, "nblcg");
    nblg.allocate(mbloc);
    do_umalloc(mbloc, 1, "nblg");
    nblkk.allocate(2, mxbli);
    do_umalloc(2*mxbli, 1, "nblkk");
    nblon.allocate(mxbli);
    do_umalloc(mxbli, 1, "nblon");
    nbs.allocate(mbloc);
    do_umalloc(mbloc, 1, "nbs");
    ncgg.allocate(mbloc);
    do_umalloc(mbloc, 1, "ncgg");
    nd11.allocate(mxbli);
    do_umalloc(mxbli, 1, "nd11");
    nd12.allocate(mxbli);
    do_umalloc(mxbli, 1, "nd12");
    nd21.allocate(mxbli);
    do_umalloc(mxbli, 1, "nd21");
    nd22.allocate(mxbli);
    do_umalloc(mxbli, 1, "nd22");
    ndat.allocate(mbloc, 6, msegn);
    do_umalloc(mbloc*6*msegn, 1, "ndat");
    ndirin.allocate(msplt*mbloc);
    do_umalloc(msplt*mbloc, 1, "ndirin");
    ndirin0.allocate(msplt);
    do_umalloc(msplt, 1, "ndirin0");
    ne1.allocate(mbloc, 6, msegn);
    do_umalloc(mbloc*6*msegn, 1, "ne1");
    ne1s.allocate(mbloc);
    do_umalloc(mbloc, 1, "ne1s");
    ne2.allocate(mbloc, 6, msegn);
    do_umalloc(mbloc*6*msegn, 1, "ne2");
    ne2s.allocate(mbloc);
    do_umalloc(mbloc, 1, "ne2s");
    nei1.allocate(mxbli);
    do_umalloc(mxbli, 1, "nei1");
    nei2.allocate(mxbli);
    do_umalloc(mxbli, 1, "nei2");
    nej1.allocate(mxbli);
    do_umalloc(mxbli, 1, "nej1");
    nej2.allocate(mxbli);
    do_umalloc(mxbli, 1, "nej2");
    nek1.allocate(mxbli);
    do_umalloc(mxbli, 1, "nek1");
    nek2.allocate(mxbli);
    do_umalloc(mxbli, 1, "nek2");
    nfb.allocate(intmax);
    do_umalloc(intmax, 1, "nfb");
    ni.allocate(mbloc);
    do_umalloc(mbloc, 1, "ni");
    nj.allocate(mbloc);
    do_umalloc(mbloc, 1, "nj");
    nk.allocate(mbloc);
    do_umalloc(mbloc, 1, "nk");
    npts.allocate(mbloc);
    do_umalloc(mbloc, 1, "npts");
    nseg.allocate(mbloc);
    do_umalloc(mbloc, 1, "nseg");
    nseg1.allocate(mbloc);
    do_umalloc(mbloc, 1, "nseg1");
    nsihi.allocate(mbloc);
    do_umalloc(mbloc, 1, "nsihi");
    nsilo.allocate(mbloc);
    do_umalloc(mbloc, 1, "nsilo");
    nsjhi.allocate(mbloc);
    do_umalloc(mbloc, 1, "nsjhi");
    nsjlo.allocate(mbloc);
    do_umalloc(mbloc, 1, "nsjlo");
    nskhi.allocate(mbloc);
    do_umalloc(mbloc, 1, "nskhi");
    nsklo.allocate(mbloc);
    do_umalloc(mbloc, 1, "nsklo");
    nsubbl.allocate(mbloc);
    do_umalloc(mbloc, 1, "nsubbl");
    nxtseg.allocate(mtot);
    do_umalloc(mtot, 1, "nxtseg");
    rkap0g.allocate(mbloc, 3);
    do_umalloc(mbloc*3, 0, "rkap0g");
    twotref.allocate(msegn, mbloc);
    do_umalloc(msegn*mbloc, 0, "twotref");
    val.allocate(mtot);
    do_umalloc(mtot, 0, "val");
    x.allocate(npmax);
    do_umalloc(npmax, 0, "x");
    xdmold.allocate(msegt);
    do_umalloc(msegt, 0, "xdmold");
    xdum.allocate(msegt, mtot);
    do_umalloc(msegt*mtot, 0, "xdum");
    xif1.allocate(intmax, nsub1);
    do_umalloc(intmax*nsub1, 1, "xif1");
    xif2.allocate(intmax, nsub1);
    do_umalloc(intmax*nsub1, 1, "xif2");
    xit1.allocate(intmax);
    do_umalloc(intmax, 1, "xit1");
    xit2.allocate(intmax);
    do_umalloc(intmax, 1, "xit2");
    xmap.allocate(msegt, msegn, mbloc);
    do_umalloc(msegt*msegn*mbloc, 0, "xmap");
    y.allocate(npmax);
    do_umalloc(npmax, 0, "y");
    z.allocate(npmax);
    do_umalloc(npmax, 0, "z");

    // set names for the temporary tlns3d map files
    std::strncpy(imapin,  "tlnsmap.in",  79); imapin[79]  = '\0';
    std::strncpy(imapout, "tlnsmap.out", 79); imapout[79] = '\0';

    // output banner
    printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf(" * *                                                               * *\n");
    printf(" * *        SPLITTER - CFL3D BLOCK AND INPUT FILE SPLITTER         * *\n");
    printf(" * *                                                               * *\n");
    printf(" * *   VERSION 6.7 :  Computational Fluids Lab, Mail Stop 128,     * *\n");
    printf(" * *                  NASA Langley Research Center, Hampton, VA    * *\n");
    printf(" * *                  Release Date:  February  1, 2017.            * *\n");
    printf(" * *                                                               * *\n");
    printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");



#ifdef DBLE_PRECSN
    printf("\n memory allocation: %12.6f Mbytes, double precision\n", (float)memuse/1.e6f);
#else
    printf("\n memory allocation: %12.6f Mbytes, single precision\n", (float)memuse/1.e6f);
#endif

    // read preliminary info (file names, etc.)
    {
        char line[256];
        FILE* fu5 = fortran_get_unit(iunit5);
        // read(iunit5,*)
        fgets(line, sizeof(line), fu5);
        // read(iunit5,'(a60)') cflinp
        fgets(line, sizeof(line), fu5);
        std::strncpy(cflinp, line, 60); cflinp[60] = '\0';
        // strip trailing newline/spaces
        for (int _i = 59; _i >= 0 && (cflinp[_i] == '\n' || cflinp[_i] == '\r' || cflinp[_i] == ' '); _i--) cflinp[_i] = '\0';
        // read(iunit5,'(a60)') roninp
        fgets(line, sizeof(line), fu5);
        std::strncpy(roninp, line, 60); roninp[60] = '\0';
        for (int _i = 59; _i >= 0 && (roninp[_i] == '\n' || roninp[_i] == '\r' || roninp[_i] == ' '); _i--) roninp[_i] = '\0';
        // read(iunit5,'(a60)') gridin
        fgets(line, sizeof(line), fu5);
        std::strncpy(gridin, line, 60); gridin[60] = '\0';
        for (int _i = 59; _i >= 0 && (gridin[_i] == '\n' || gridin[_i] == '\r' || gridin[_i] == ' '); _i--) gridin[_i] = '\0';
        // read(iunit5,'(a60)') sdgridin
        fgets(line, sizeof(line), fu5);
        std::strncpy(sdgridin, line, 60); sdgridin[60] = '\0';
        for (int _i = 59; _i >= 0 && (sdgridin[_i] == '\n' || sdgridin[_i] == '\r' || sdgridin[_i] == ' '); _i--) sdgridin[_i] = '\0';
        // read(iunit5,*)
        fgets(line, sizeof(line), fu5);
        // read(iunit5,*) icflin,ironin,ibin,isdin
        fortran_read_list(iunit5, &icflin, &ironin, &ibin, &isdin);
        if (icflin == -5) icflin = 5;
        // read(iunit5,*)
        fgets(line, sizeof(line), fu5);
        // read(iunit5,'(a60)') cflout
        fgets(line, sizeof(line), fu5);
        std::strncpy(cflout, line, 60); cflout[60] = '\0';
        for (int _i = 59; _i >= 0 && (cflout[_i] == '\n' || cflout[_i] == '\r' || cflout[_i] == ' '); _i--) cflout[_i] = '\0';
        // read(iunit5,'(a60)') ronout
        fgets(line, sizeof(line), fu5);
        std::strncpy(ronout, line, 60); ronout[60] = '\0';
        for (int _i = 59; _i >= 0 && (ronout[_i] == '\n' || ronout[_i] == '\r' || ronout[_i] == ' '); _i--) ronout[_i] = '\0';
        // read(iunit5,'(a60)') gridout
        fgets(line, sizeof(line), fu5);
        std::strncpy(gridout, line, 60); gridout[60] = '\0';
        for (int _i = 59; _i >= 0 && (gridout[_i] == '\n' || gridout[_i] == '\r' || gridout[_i] == ' '); _i--) gridout[_i] = '\0';
        // read(iunit5,'(a60)') sdgridout
        fgets(line, sizeof(line), fu5);
        std::strncpy(sdgridout, line, 60); sdgridout[60] = '\0';
        for (int _i = 59; _i >= 0 && (sdgridout[_i] == '\n' || sdgridout[_i] == '\r' || sdgridout[_i] == ' '); _i--) sdgridout[_i] = '\0';
        // read(iunit5,*)
        fgets(line, sizeof(line), fu5);
        // read(iunit5,*) icflout,ironout,ibout,isdout
        fortran_read_list(iunit5, &icflout, &ironout, &ibout, &isdout);
        if (icflout == -5) icflout = 5;
        // read(iunit5,*)
        fgets(line, sizeof(line), fu5);
    }

    printf("\n");
    printf("input (unsplit) files\n");
    printf("  %-60.60s\n", cflinp);
    printf("  %-60.60s\n", roninp);
    printf("  %-60.60s\n", gridin);
    printf("  %-60.60s\n", sdgridin);
    printf("icflver   ironver   igrdfmt    isdfmt\n");
    printf("%7d%10d%10d%10d\n", icflin, ironin, ibin, isdin);
    printf("output (split) files\n");

    // if any one of the input files is 'null', then corresponding
    // output file must be 'null' too! also, can't split ronnie
    // input files without also splitting cfl3d input file
    if (is_null_str(cflinp)) {
        printf("STOPPING: must provide a cfl3d input file!\n");
        { int _t0=0, _tm2=-2; termn8(_t0, _tm2, ibufdim, nbuf, bou, nou); }
    }
    if (is_null_str(roninp))
        std::strncpy(ronout, "null", 79);
    if (is_null_str(gridin))
        std::strncpy(gridout, "null", 79);
    if (is_null_str(sdgridin)) {
        std::strncpy(sdgridout, "null", 79);
        std::strncpy(dovrlap,   "null", 79);
        std::strncpy(dpatch,    "null", 79);
        std::strncpy(dresid,    "null", 79);
    }
    printf("  %-60.60s\n", cflout);
    printf("  %-60.60s\n", ronout);
    printf("  %-60.60s\n", gridout);
    printf("  %-60.60s\n", sdgridout);
    printf("icflver   ironver   igrdfmt    isdfmt\n");
    printf("%7d%10d%10d%10d\n", icflout, ironout, ibout, isdout);
    printf("\n");

    // parallel input version 4 of cfl3d denoted with icflin < 0
    // and or icflout < 0 on input; there is no distinction between
    // parallel and sequential input files for version 6/5
    iparin  = 0;
    iparout = 0;
    if (icflin < 0) {
        iparin = 1;
        icflin = std::abs(icflin);
    }
    if (icflout < 0) {
        iparout = 1;
        icflout = std::abs(icflout);
    }
    if (std::abs(icflin) > 4) {
        icflin = std::abs(icflin);
    }
    if (std::abs(icflout) > 4) {
        icflout = std::abs(icflout);
    }



    // open files
    fortran_open_unit(10, cflinp, "r");
    rewind(10);
    if (!is_null_str(cflout)) {
        fortran_open_unit(7, cflout, "w");
        rewind(7);
    }
    ipatch0 = 0;
    if (!is_null_str(roninp)) {
        ipatch0 = 1;
        fortran_open_unit(15, roninp, "r");
        rewind(15);
    }
    if (!is_null_str(ronout)) {
        fortran_open_unit(9, ronout, "w");
        rewind(9);
    }
    fortran_open_unit(20, imapin, "w");
    rewind(20);
    if (!is_null_str(gridin)) {
        if (ibin == 1) {
            // the following line is needed to use ieee binary file
            // call asnfile (gridin, '-F f77 -N ieee' , IER)
            fortran_open_unit(30, gridin, "rb");
        } else {
            fortran_open_unit(30, gridin, "r");
        }
        rewind(30);
    }
    if (!is_null_str(gridout)) {
        if (ibout == 1) {
            // the following line is needed to use ieee binary file
            // call asnfile (gridout, '-F f77 -N ieee' , IER)
            fortran_open_unit(40, gridout, "wb");
        } else {
            fortran_open_unit(40, gridout, "w");
        }
        rewind(40);
    }
    fortran_open_unit(8, imapout, "w");
    rewind(8);
    isd = 1;
    if (!is_null_str(sdgridin)) {
        if (isdin == 1) {
            fortran_open_unit(50, sdgridin, "rb");
        } else {
            fortran_open_unit(50, sdgridin, "r");
        }
        rewind(50);
    }
    if (!is_null_str(sdgridout)) {
        if (isdout == 1) {
            fortran_open_unit(60, sdgridout, "wb");
        } else {
            fortran_open_unit(60, sdgridout, "w");
        }
        rewind(60);
    }

    // convert cfl3d input file (and ronnie input file, if
    // applicable) to a tlns3d map file
    printf("converting unsplit cfl3d input file to tlns3d map file\n");
    printf("\n");

    // Build 2D ref for imap (sequence association: 3D -> 2D)
    // imap(msegt,msegn,mbloc) passed as imap(msegt,msegn) to cfl3d_to_tlns3d
    // cfl3d_to_tlns3d takes FortranArray3DRef<int> imap - OK, pass directly

    cfl3d_to_tlns3d(icflin, ipatch0, ironin, iparin, nnodes, isd,
                    mbloc, msegn, msegt, mxbli, il, jl, kl,
                    rkap0g, levelg, igridg, iflimg, ifdsg,
                    iviscg, jdimg, kdimg, idimg, idiagg,
                    nblcg, idegg, jsg, ksg, jeg, keg,
                    ieg, mit, ilamlog, ilamhig, jlamlog,
                    jlamhig, klamlog, klamhig,
                    iredundant, nseg1, nseg, nsilo, nsihi,
                    nsjlo, nsjhi, nsklo, nskhi, nb1, ne1,
                    nb2, ne2, ibct, bcval, iovrlp, ifoflg,
                    ndat, xmap, imap, ivisb, twotref, itrb1,
                    itrb2, jtrb1, jtrb2, ktrb1, ktrb2, iturbb,
                    ibif1, ibif2, nbi1, nei1, nbj1, nej1, nbk1,
                    nek1, nbi2, nei2, nbj2, nej2, nbk2, nek2,
                    nd11, nd21, nd12, nd22, ifsor, nb1s, ne1s,
                    nb2s, ne2s, nrotat, ntrans, tlref, rlref, ioflag,
                    nou, bou, nbuf, ibufdim, cflout);
    rewind(20);

    nbl = nbloc;

    // check dimensions
    printf("checking dimensions...\n");
    printf("\n");

    if (nbl > mbloc) {
        printf("Number of blocks (%3d) is larger than the dimensioned value (%3d)\n", nbl, mbloc);
        { int _t0=0, _tm3=-3; termn8(_t0, _tm3, ibufdim, nbuf, bou, nou); }
    }

    iflag = 0;
    ibeg(1) = 1;
    npts(1) = il(1) * jl(1) * kl(1);
    ntot = ibeg(1) + npts(1) - 1;
    for (int n = 2; n <= nbl; n++) {
        ibeg(n) = ibeg(n-1) + npts(n-1);
        npts(n) = il(n) * jl(n) * kl(n);
        ntot = ibeg(n) + npts(n) - 1;
    }
    if (ntot > npmax) {
        printf("number of points %8d is larger than the dimensioned value %8d\n", ntot, npmax);
        iflag = 1;
    }

    if (iflag == 1) { int _t0=0, _tm4=-4; termn8(_t0, _tm4, ibufdim, nbuf, bou, nou); }



    // read in the input grid file; if input grid file is 'null', get
    // needed block dimensions from unsplit cfl3d input file
    if (!is_null_str(gridin)) {

        if (ip3dgrd > 0) {
            if (ibin == 1) {
                int nbl1;
                FILE* f30 = fortran_get_unit(30);
                fread(&nbl1, sizeof(int), 1, f30);
                if (nbl1 != nbl) {
                    { int _t0=0, _tm1=-1; termn8(_t0, _tm1, ibufdim, nbuf, bou, nou); }
                }
                for (int n = 1; n <= nbl; n++) {
                    fread(&ilo(n), sizeof(int), 1, f30);
                    fread(&jlo(n), sizeof(int), 1, f30);
                    fread(&klo(n), sizeof(int), 1, f30);
                }
            } else {
                int nbl1;
                fortran_read_list(30, &nbl1);
                if (nbl1 != nbl) {
                    { int _t0=0, _tm1=-1; termn8(_t0, _tm1, ibufdim, nbuf, bou, nou); }
                }
                // read(30,*) (ilo(n),jlo(n),klo(n),n=1,nbl) - ONE implied-DO read
                // reads 3*nbl interleaved values from one (possibly multi-line) record
                {
                    FortranArray1D<int> tmp_dims(3*nbl);
                    fortran_read_list_array(30, tmp_dims, 1, 3*nbl);
                    for (int n = 1; n <= nbl; n++) {
                        ilo(n) = tmp_dims(3*(n-1)+1);
                        jlo(n) = tmp_dims(3*(n-1)+2);
                        klo(n) = tmp_dims(3*(n-1)+3);
                    }
                }
            }
            for (int n = 1; n <= nbl; n++) {
                if (ilo(n) != il(n) || jlo(n) != jl(n) || klo(n) != kl(n)) {
                    { int _t0=0, _tm1=-1; termn8(_t0, _tm1, ibufdim, nbuf, bou, nou); }
                }
            }
        }

        printf("reading grid...\n");
        printf("grid: %-60.60s\n", gridin);
        printf("\n");

        if (ibin == 1) {
            for (int n = 1; n <= nbl; n++) {
                if (ip3dgrd == 0) {
                    FILE* f30 = fortran_get_unit(30);
                    fread(&jlo(n), sizeof(int), 1, f30);
                    fread(&klo(n), sizeof(int), 1, f30);
                    fread(&ilo(n), sizeof(int), 1, f30);
                    if (ilo(n) != il(n) || jlo(n) != jl(n) || klo(n) != kl(n)) {
                        { int _t0=0, _tm1=-1; termn8(_t0, _tm1, ibufdim, nbuf, bou, nou); }
                    }
                }
                // call ingrd(30,ibin,x(ibeg(n)),y(ibeg(n)),z(ibeg(n)),ilo(n),jlo(n),klo(n),ip3dgrd)
                // x,y,z are 1D but ingrd expects 3D - use sequence association via FortranArray3DRef
                {
                    int ibgn = ibeg(n);
                    FortranArray3DRef<double> xref(x.data() + ibgn - 1, ilo(n), jlo(n), klo(n));
                    FortranArray3DRef<double> yref(y.data() + ibgn - 1, ilo(n), jlo(n), klo(n));
                    FortranArray3DRef<double> zref(z.data() + ibgn - 1, ilo(n), jlo(n), klo(n));
                    { int _nf30=30; ingrd(_nf30, ibin, xref, yref, zref, ilo(n), jlo(n), klo(n), ip3dgrd); }
                }
            }
        } else {
            for (int n = 1; n <= nbl; n++) {
                if (ip3dgrd == 0) {
                    fortran_read_list(30, &jlo(n), &klo(n), &ilo(n));
                    if (ilo(n) != il(n) || jlo(n) != jl(n) || klo(n) != kl(n)) {
                        { int _t0=0, _tm1=-1; termn8(_t0, _tm1, ibufdim, nbuf, bou, nou); }
                    }
                }
                {
                    int ibgn = ibeg(n);
                    FortranArray3DRef<double> xref(x.data() + ibgn - 1, ilo(n), jlo(n), klo(n));
                    FortranArray3DRef<double> yref(y.data() + ibgn - 1, ilo(n), jlo(n), klo(n));
                    FortranArray3DRef<double> zref(z.data() + ibgn - 1, ilo(n), jlo(n), klo(n));
                    { int _nf30=30; ingrd(_nf30, ibin, xref, yref, zref, ilo(n), jlo(n), klo(n), ip3dgrd); }
                }
            }
        }
        ntot = ibeg(nbl) + npts(nbl) - 1;

    } else {

        for (int n = 1; n <= nbl; n++) {
            ilo(n) = il(n);
            jlo(n) = jl(n);
            klo(n) = kl(n);
        }

    }

    for (int n = 1; n <= nbl; n++) {
        imin(n) = 0;
        imax(n) = 0;
        iskp(n) = 0;
        jmin(n) = 0;
        jmax(n) = 0;
        jskp(n) = 0;
        kmin(n) = 0;
        kmax(n) = 0;
        kskp(n) = 0;
    }

    for (int n = 1; n <= nbl; n++) {
        printf("block #%4d: il=%4d, jl=%4d, kl=%4d\n", n, ilo(n), jlo(n), klo(n));
        iln(n) = ilo(n);
        jln(n) = jlo(n);
        kln(n) = klo(n);
        nbs(n) = n;
    }



    // read in the tlns3d map file
    // mapin takes FortranArray2DRef<int> imap - but imap is 3D(msegt,msegn,mbloc)
    // Fortran sequence association: pass as 2D(msegt*msegn, mbloc)
    {
        // mapin declares imap(msegt,mtot) and indexes imap(1..23,nsgtop): the
        // leading dim it sees is msegt, NOT msegt*msegn. The caller's 3D
        // imap(msegt,msegn,mbloc) sequence-associates as imap(msegt, msegn*mbloc).
        // Passing leading dim msegt*msegn sent imap(n,nsgtop>=2) far out of
        // bounds (into the y grid array). Leading dim must be msegt.
        FortranArray2DRef<int> imap_2d(imap.data(), msegt, msegn*mbloc);
        int _u20=20, _u6=6;
        mapin(_u20, _u6, nbl, iln, jln, kln, mbloc, msegt, mtot,
              imap_2d, idbloc, ivisb, itrb, val, xdum, iold,
              nxtseg, intrfc, ipatch, nsubbl, idobl, nseg, idno,
              ijk, idseg, idnext, nou, bou, nbuf, ibufdim);
    }

    // read splitting directives and do the splitting
    fortran_read_list(iunit5, &nout);
    // shortcut to split all blocks the same as the single block input
    if (nout < 0) {
        nout = std::abs(nout);
        int isblocdum;
        for (int no = 1; no <= nout; no++) {
            fortran_read_list(iunit5, &isblocdum);
            fortran_read_list(iunit5, &ndirin0(no));
            fortran_read_list(iunit5, &isin0(no));
        }
        int nn = 0;
        for (int nb = 1; nb <= nbl; nb++) {
            for (int no = 1; no <= nout; no++) {
                nn = nn + 1;
                isblocin(nn) = nb;
                ndirin(nn)   = ndirin0(no);
                isin(nn)     = isin0(no);
            }
        }
        nout = nn;
    } else {
        for (int nn = 1; nn <= nout; nn++) {
            fortran_read_list(iunit5, &isblocin(nn));
            fortran_read_list(iunit5, &ndirin(nn));
            fortran_read_list(iunit5, &isin(nn));
        }
    }
    printf("\n");
    printf("number of splits = %3d\n", nout);
    int ntotal = nbl + nout;
    if (ntotal > mbloc) {
        printf("Number of output blocks (%3d) is larger than the dimensioned value (%3d)\n", nbl, mbloc);
        { int _t0=0, _tm3=-3; termn8(_t0, _tm3, ibufdim, nbuf, bou, nou); }
    }
    printf("  split  block  coord  index\n");
    for (int no = 1; no <= nout; no++) {

        int isbloc = isblocin(no);
        int ndir   = ndirin(no);
        int is     = isin(no);
        if (is == 0) {
            if (ndir == 1) {
                is = (iln(isbloc)+1) / 2;
            } else if (ndir == 2) {
                is = (jln(isbloc)+1) / 2;
            } else {
                is = (kln(isbloc)+1) / 2;
            }
        }
        if (ndir == 1) {
            printf("%7d%7d      I%7d\n", no, isbloc, is);
        } else if (ndir == 2) {
            printf("%7d%7d      J%7d\n", no, isbloc, is);
        } else {
            printf("%7d%7d      K%7d\n", no, isbloc, is);
        }

        {
            FortranArray2DRef<int> imap_2d(imap.data(), msegt*msegn, mbloc);
            mapsplt(isbloc, ndir, is, mbloc, msegt, mtot,
                    imap_2d, idbloc, ivisb, itrb, val, xdum, iold,
                    nxtseg, intrfc, ipatch, nsubbl, idobl,
                    nseg, idno, ijk, idseg, idnext, xdmold, idmold);
        }

    }

    {
        FortranArray2DRef<int> imap_2d(imap.data(), msegt*msegn, mbloc);
        int _u8=8;
        mapout(_u8, nbloc, ni, nj, nk, mbloc, msegt, mtot,
               imap_2d, idbloc, ivisb, itrb, val, xdum, iold,
               nxtseg, intrfc, ipatch, nsubbl, idobl, nseg, idno,
               ijk, idseg, idnext);
    }

    {
        FortranArray2DRef<int> imap_2d(imap.data(), msegt*msegn, mbloc);
        mapblk(nblk, nbs, imin, imax, jmin, jmax, kmin, kmax, mbloc,
               msegt, mtot, imap_2d, idbloc, ivisb, itrb, val, xdum,
               iold, nxtseg, intrfc, ipatch, nsubbl, idobl,
               nseg, idno, ijk, idseg, idnext);
    }
    for (int n = 1; n <= nblk; n++) {
        iskp(n) = 1;
        jskp(n) = 1;
        kskp(n) = 1;
        iln(n) = (imax(n)-imin(n))/iskp(n) + 1;
        jln(n) = (jmax(n)-jmin(n))/jskp(n) + 1;
        kln(n) = (kmax(n)-kmin(n))/kskp(n) + 1;
    }



    // check multigridability of split-grid basic dimensions
    ncgmax = 0;
    iflg = 1;
    for (int nnn = 1; nnn <= 5; nnn++) {
        int nn = nnn - 1;
        for (int n = 1; n <= nblk; n++) {
            int idtest, jdtest, kdtest;
            if (iln(n) > 2) {
                idtest = (iln(n)-1) / (1 << nn) + 1;
            } else {
                idtest = 1;
            }
            jdtest = (jln(n)-1) / (1 << nn) + 1;
            kdtest = (kln(n)-1) / (1 << nn) + 1;
            if ((idtest/2*2 == idtest) ||
                (jdtest/2*2 == jdtest) ||
                (kdtest/2*2 == kdtest)) {
                iflg = 0;
            }
        }
        if (iflg != 0) ncgmax = ncgmax + 1;
    }

    // check that coarsened dimensions don't drop below 3 in
    // j or k directions (or i if not a 2d case) - 2 cells
    // in any direction other than i (for 2d) causes problems
    // with 2nd order scheme since 2 ghost cells cannot be
    // filled properly at block interfaces
    iflg = 1;
    ncgmax1 = 0;
    if (ncgmax > 0) {
        for (int nnn = 1; nnn <= ncgmax; nnn++) {
            int nn = nnn;
            for (int n = 1; n <= nblk; n++) {
                int idtest = (iln(n)-1) / (1 << nn) + 1;
                int jdtest = (jln(n)-1) / (1 << nn) + 1;
                int kdtest = (kln(n)-1) / (1 << nn) + 1;
                if (iln(n) > 2) {
                    if (idtest <= 2) iflg = 0;
                }
                if (jdtest <= 2 || kdtest <= 2) iflg = 0;
            }
            if (iflg != 0) ncgmax1 = ncgmax1 + 1;
        }
    }

    ncgmax = std::min(ncgmax, ncgmax1);

    printf("\n");
    printf("split-grid basic dimensions are multigridable to ncg = %2d\n", ncgmax);
    if (iflg == 0) {
        printf("NOTE: ncg was limited to prevent coarsened dimensions less than 3\n");
    }

    // print out new grid
    if (!is_null_str(gridout)) {
        if (ibout == 0) {

            if (ip3dgrd > 0) {
                fortran_write_unit(40, "%5d\n", nblk);
                for (int n = 1; n <= nblk; n++) {
                    fortran_write_unit(40, "%5d%5d%5d", iln(n), jln(n), kln(n));
                }
                fortran_write_unit(40, "\n");
            }
            for (int n = 1; n <= nblk; n++) {
                if (ip3dgrd == 0) {
                    fortran_write_unit(40, "%5d%5d%5d\n", jln(n), kln(n), iln(n));
                }
                int isk = iskp(n);
                int jsk = jskp(n);
                int ksk = kskp(n);
                int nn  = nbs(n);
                {
                    int ibgn = ibeg(nn);
                    FortranArray3DRef<double> xref(x.data() + ibgn - 1, ilo(nn), jlo(nn), klo(nn));
                    FortranArray3DRef<double> yref(y.data() + ibgn - 1, ilo(nn), jlo(nn), klo(nn));
                    FortranArray3DRef<double> zref(z.data() + ibgn - 1, ilo(nn), jlo(nn), klo(nn));
                    int _u40=40;
                    outgrd(_u40, ibout, xref, yref, zref, ilo(nn), jlo(nn), klo(nn),
                           imin(n), imax(n), isk,
                           jmin(n), jmax(n), jsk,
                           kmin(n), kmax(n), ksk, ip3dgrd);
                }
            }

        } else {

            if (ip3dgrd > 0) {
                // Fortran: write(40) nblk        (record 1)
                //          write(40) (iln,jln,kln,n=1,nblk)  (record 2)
                // Each unformatted sequential record carries leading+trailing
                // 4-byte length markers — emit them so the grid is a valid
                // Fortran-unformatted file (byte-exact with the reference).
                FILE* f40 = fortran_get_unit(40);
                int rl1 = (int)sizeof(int);
                fwrite(&rl1, sizeof(int), 1, f40);
                fwrite(&nblk, sizeof(int), 1, f40);
                fwrite(&rl1, sizeof(int), 1, f40);
                int rl2 = 3 * nblk * (int)sizeof(int);
                fwrite(&rl2, sizeof(int), 1, f40);
                for (int n = 1; n <= nblk; n++) {
                    fwrite(&iln(n), sizeof(int), 1, f40);
                    fwrite(&jln(n), sizeof(int), 1, f40);
                    fwrite(&kln(n), sizeof(int), 1, f40);
                }
                fwrite(&rl2, sizeof(int), 1, f40);
            }
            for (int n = 1; n <= nblk; n++) {
                if (ip3dgrd == 0) {
                    // Fortran: write(40) jln(n),kln(n),iln(n)  — one record.
                    FILE* f40 = fortran_get_unit(40);
                    int rl = 3 * (int)sizeof(int);
                    fwrite(&rl, sizeof(int), 1, f40);
                    fwrite(&jln(n), sizeof(int), 1, f40);
                    fwrite(&kln(n), sizeof(int), 1, f40);
                    fwrite(&iln(n), sizeof(int), 1, f40);
                    fwrite(&rl, sizeof(int), 1, f40);
                }
                int isk = iskp(n);
                int jsk = jskp(n);
                int ksk = kskp(n);
                int nn  = nbs(n);
                {
                    int ibgn = ibeg(nn);
                    FortranArray3DRef<double> xref(x.data() + ibgn - 1, ilo(nn), jlo(nn), klo(nn));
                    FortranArray3DRef<double> yref(y.data() + ibgn - 1, ilo(nn), jlo(nn), klo(nn));
                    FortranArray3DRef<double> zref(z.data() + ibgn - 1, ilo(nn), jlo(nn), klo(nn));
                    int _u40=40;
                    outgrd(_u40, ibout, xref, yref, zref, ilo(nn), jlo(nn), klo(nn),
                           imin(n), imax(n), isk,
                           jmin(n), jmax(n), jsk,
                           kmin(n), kmax(n), ksk, ip3dgrd);
                }
            }

        }
    }

    {
        int np = 0;
        for (int n = 1; n <= nbl; n++) {
            np = np + ilo(n)*jlo(n)*klo(n);
        }
        printf("\n");
        printf("Input  points: %8d\n", np);
        np = 0;
        for (int n = 1; n <= nblk; n++) {
            np = np + iln(n)*jln(n)*kln(n);
        }
        printf("Ouput  points: %8d\n", np);
        printf("\n");
    }



    // create split cfl3d input file (and split ronnie
    // input file if applicable)
    if (!is_null_str(cflout)) {
        printf("converting split tlns3d map file to cfl3d input file\n");
        // set number of compute nodes = number of (split) blocks
        // for parallel computations
        nnodes = nblk;
        tlns3d_to_cfl3d(ibout, icflout, ipatch0, ironout, nblk,
                        iln, jln, kln, iparout, nnodes, isd, ifrom, xif1, xif2, etf1,
                        etf2, ito, xit1, xit2, ett1, ett2, nfb, iredundant, iconcat,
                        xmap, msegn, mbloc, iovrlp, nseg, imap, ivisb, twotref,
                        itrb1, itrb2, jtrb1, jtrb2, ktrb1, ktrb2, iturbb,
                        cq, cporous, cpchamber, cradiation,
                        llimit, iitmax, mmcxie, mmceta, iifit, iic0, iiorph,
                        bcvali, bcvalj, bcvalk, nbci0, nbcidim, nbcj0, nbcjdim,
                        nbck0, nbckdim, ibcinfo, jbcinfo, kbcinfo, nblkk, limblk,
                        isva, nblon, rkap0g, levelg, igridg, iflimg, ifdsg, iviscg,
                        jdimg, kdimg, idimg, idiagg, nblcg, idegg, jsg, ksg, isg,
                        jeg, keg, ieg, mit, ilamlog, ilamhig, jlamlog, jlamhig,
                        klamlog, klamhig, ncgg, nblg, iemg, inewgg, iadvance, iforce,
                        intmax, nsub1, msegt, maxseg, mxbli, nrotat, ntrans, tlref,
                        rlref, ioflag, nou, bou, nbuf, ibufdim);
    }

    // read and split SD grid, if needed
    if (is_null_str(sdgridin)) isd = 0;
    if (is_null_str(sdgridout)) isd = 0;

    if (isd > 0) {

        printf("\n");
        printf("reading SD grid...\n");

        if (isdin == 0) {

            FILE* f50 = fortran_get_unit(50);
            // read(50,'(a80)') header
            fgets(header, 81, f50);
            // read(50,*) nblsd,ndv
            fortran_read_list(50, &nblsd, &ndv);
            // read(50,*) (ilosd(n),jlosd(n),klosd(n),n=1,nblsd) - ONE implied-DO read
            {
                FortranArray1D<int> tmp_dims(3*nblsd);
                fortran_read_list_array(50, tmp_dims, 1, 3*nblsd);
                for (int n = 1; n <= nblsd; n++) {
                    ilosd(n) = tmp_dims(3*(n-1)+1);
                    jlosd(n) = tmp_dims(3*(n-1)+2);
                    klosd(n) = tmp_dims(3*(n-1)+3);
                }
            }

        } else {

            FILE* f50 = fortran_get_unit(50);
            // read(50) header
            fread(header, 1, 80, f50);
            // read(50) nblsd,ndv
            fread(&nblsd, sizeof(int), 1, f50);
            fread(&ndv, sizeof(int), 1, f50);
            // read(50) (ilosd(n),jlosd(n),klosd(n),n=1,nblsd)
            for (int n = 1; n <= nblsd; n++) {
                fread(&ilosd(n), sizeof(int), 1, f50);
                fread(&jlosd(n), sizeof(int), 1, f50);
                fread(&klosd(n), sizeof(int), 1, f50);
            }

        }

        printf("SD grid: %-60.60s\n", sdgridin);
        printf("\n");
        printf("number of design variables = %4d\n", ndv);
        for (int n = 1; n <= nblsd; n++) {
            printf("block #%4d: il=%4d, jl=%4d, kl=%4d\n", n, ilosd(n), jlosd(n), klosd(n));
        }
        printf("\n");

        // check consistancy with grid file
        if (nblsd != nbl) {
            fprintf(fortran_get_unit(6), " stopping: number of blocks in grid file    = %d\n", nbl);
            fprintf(fortran_get_unit(6), " stopping: number of blocks in grid SD file = %d\n", nblsd);
            { int _t0=0, _tm5=-5; termn8(_t0, _tm5, ibufdim, nbuf, bou, nou); }
        }
        for (int n = 1; n <= nbl; n++) {
            if (ilo(n) != ilosd(n) || jlo(n) != jlosd(n) || klo(n) != klosd(n)) {
                fprintf(fortran_get_unit(6), " stopping: mismatch in grid dimensions and SD dimensions in block %d\n", nbl);
                fprintf(fortran_get_unit(6), " grid i,j,k = %d %d %d\n", ilo(n), jlo(n), klo(n));
                fprintf(fortran_get_unit(6), " SD   i,j,k = %d %d %d\n", ilosd(n), jlosd(n), klosd(n));
                { int _t0=0, _tm5=-5; termn8(_t0, _tm5, ibufdim, nbuf, bou, nou); }
            }
        }

        if (!is_null_str(sdgridout)) {
            if (isdout == 0) {
                fortran_write_unit(60, "%-80.80s\n", header);
                fortran_write_unit(60, " %d %d\n", nblk, ndv);
                for (int n = 1; n <= nblk; n++) {
                    fortran_write_unit(60, " %d %d %d", iln(n), jln(n), kln(n));
                }
                fortran_write_unit(60, "\n");
            } else {
                FILE* f60 = fortran_get_unit(60);
                fwrite(header, 1, 80, f60);
                fwrite(&nblk, sizeof(int), 1, f60);
                fwrite(&ndv, sizeof(int), 1, f60);
                for (int n = 1; n <= nblk; n++) {
                    fwrite(&iln(n), sizeof(int), 1, f60);
                    fwrite(&jln(n), sizeof(int), 1, f60);
                    fwrite(&kln(n), sizeof(int), 1, f60);
                }
            }
        }



        for (int n = 1; n <= nblk; n++) {
            // write(*,'(1x,8i5)') n,imin(n),imax(n),jmin(n),jmax(n),
            //                     kmin(n),kmax(n),nbs(n)
            int isk = iskp(n);
            int jsk = jskp(n);
            int ksk = kskp(n);
            int nn  = nbs(n);

            // read grid sensitivity derivatives for (unsplit) block nn
            rewind(50);
            if (isdin == 0) {
                FILE* f50 = fortran_get_unit(50);
                char lbuf[256];
                fgets(lbuf, sizeof(lbuf), f50); // read(50,*)
                fgets(lbuf, sizeof(lbuf), f50); // read(50,*)
                fgets(lbuf, sizeof(lbuf), f50); // read(50,*)
            } else {
                FILE* f50 = fortran_get_unit(50);
                char lbuf[256];
                fread(lbuf, 1, 80, f50); // read(50)
                int tmp; fread(&tmp, sizeof(int), 1, f50); // read(50)
                fread(&tmp, sizeof(int), 1, f50);
                fread(lbuf, 1, 1, f50); // read(50)
            }
            for (int mm = 1; mm <= nbl; mm++) {
                if (mm != nn) {
                    int ilos = ilosd(mm);
                    int jlos = jlosd(mm);
                    int klos = klosd(mm);
                    for (int idv = 1; idv <= ndv; idv++) {
                        if (isdin == 0) {
                            // read(50,*) (((xdum0,...),...),...), (((ydum0,...),...),...), (((zdum0,...),...),...) 
                            // Each is ONE read statement reading ilos*jlos*klos values - use temp array
                            int ntmp = ilos*jlos*klos;
                            FortranArray1D<double> tmp_arr(ntmp);
                            fortran_read_list_array(50, tmp_arr, 1, ntmp);
                            fortran_read_list_array(50, tmp_arr, 1, ntmp);
                            fortran_read_list_array(50, tmp_arr, 1, ntmp);
                        } else {
                            FILE* f50 = fortran_get_unit(50);
                            for (int _cnt = 0; _cnt < ilos*jlos*klos; _cnt++) {
                                fread(&xdum0, sizeof(double), 1, f50);
                            }
                            for (int _cnt = 0; _cnt < ilos*jlos*klos; _cnt++) {
                                fread(&ydum0, sizeof(double), 1, f50);
                            }
                            for (int _cnt = 0; _cnt < ilos*jlos*klos; _cnt++) {
                                fread(&zdum0, sizeof(double), 1, f50);
                            }
                        }
                    }
                } else {
                    for (int idv = 1; idv <= ndv; idv++) {
                        if (isdin == 0) {
                            // read(50,*) (x(i),i=ibeg(nn),ibeg(nn)+npts(nn)-1), ...
                            // ONE read statement reads all values - use fortran_read_list_array
                            int ibgn = ibeg(nn);
                            int iend = ibeg(nn) + npts(nn) - 1;
                            fortran_read_list_array(50, x, ibgn, iend);
                            fortran_read_list_array(50, y, ibgn, iend);
                            fortran_read_list_array(50, z, ibgn, iend);
                        } else {
                            FILE* f50 = fortran_get_unit(50);
                            int ibgn = ibeg(nn);
                            int iend = ibeg(nn) + npts(nn) - 1;
                            for (int _i = ibgn; _i <= iend; _i++)
                                fread(&x(_i), sizeof(double), 1, f50);
                            for (int _i = ibgn; _i <= iend; _i++)
                                fread(&y(_i), sizeof(double), 1, f50);
                            for (int _i = ibgn; _i <= iend; _i++)
                                fread(&z(_i), sizeof(double), 1, f50);
                        }
                        {
                            int ibgn = ibeg(nn);
                            FortranArray3DRef<double> xref(x.data() + ibgn - 1, ilo(nn), jlo(nn), klo(nn));
                            FortranArray3DRef<double> yref(y.data() + ibgn - 1, ilo(nn), jlo(nn), klo(nn));
                            FortranArray3DRef<double> zref(z.data() + ibgn - 1, ilo(nn), jlo(nn), klo(nn));
                            int _u60=60, _ip3d=1;
                            outgrd(_u60, isdout, xref, yref, zref, ilo(nn), jlo(nn), klo(nn),
                                   imin(n), imax(n), isk,
                                   jmin(n), jmax(n), jsk,
                                   kmin(n), kmax(n), ksk, _ip3d);
                        }
                    }
                }
            }

        } // end do 740 n=1,nblk

    } // end if (isd > 0)



    // free up memory
    ifrom.deallocate();
    xif1.deallocate();
    xif2.deallocate();
    etf1.deallocate();
    etf2.deallocate();
    ito.deallocate();
    xit1.deallocate();
    xit2.deallocate();
    ett1.deallocate();
    ett2.deallocate();
    nfb.deallocate();
    iadvance.deallocate();
    iforce.deallocate();
    nblkk.deallocate();
    limblk.deallocate();
    isva.deallocate();
    nblon.deallocate();
    ncgg.deallocate();
    nblg.deallocate();
    iemg.deallocate();
    inewgg.deallocate();
    bcvali.deallocate();
    bcvalj.deallocate();
    bcvalk.deallocate();
    nbci0.deallocate();
    nbcidim.deallocate();
    nbcj0.deallocate();
    nbcjdim.deallocate();
    nbck0.deallocate();
    nbckdim.deallocate();
    ibcinfo.deallocate();
    jbcinfo.deallocate();
    kbcinfo.deallocate();
    llimit.deallocate();
    iitmax.deallocate();
    mmcxie.deallocate();
    mmceta.deallocate();
    iifit.deallocate();
    iic0.deallocate();
    iiorph.deallocate();
    iredundant.deallocate();
    iconcat.deallocate();
    cq.deallocate();
    cporous.deallocate();
    cpchamber.deallocate();
    cradiation.deallocate();
    isblocin.deallocate();
    ndirin.deallocate();
    isin.deallocate();
    ndirin0.deallocate();
    isin0.deallocate();
    nbs.deallocate();
    imin.deallocate();
    jmin.deallocate();
    kmin.deallocate();
    imax.deallocate();
    jmax.deallocate();
    kmax.deallocate();
    iskp.deallocate();
    jskp.deallocate();
    kskp.deallocate();
    ilo.deallocate();
    jlo.deallocate();
    klo.deallocate();
    iln.deallocate();
    jln.deallocate();
    kln.deallocate();
    ibeg.deallocate();
    npts.deallocate();
    x.deallocate();
    y.deallocate();
    z.deallocate();
    ilosd.deallocate();
    jlosd.deallocate();
    klosd.deallocate();
    il.deallocate();
    jl.deallocate();
    kl.deallocate();
    rkap0g.deallocate();
    levelg.deallocate();
    igridg.deallocate();
    iflimg.deallocate();
    ifdsg.deallocate();
    iviscg.deallocate();
    jdimg.deallocate();
    kdimg.deallocate();
    idimg.deallocate();
    idiagg.deallocate();
    nblcg.deallocate();
    idegg.deallocate();
    jsg.deallocate();
    ksg.deallocate();
    isg.deallocate();
    jeg.deallocate();
    keg.deallocate();
    ieg.deallocate();
    mit.deallocate();
    ilamlog.deallocate();
    ilamhig.deallocate();
    jlamlog.deallocate();
    jlamhig.deallocate();
    klamlog.deallocate();
    klamhig.deallocate();
    imap.deallocate();
    idbloc.deallocate();
    ivisb.deallocate();
    itrb.deallocate();
    val.deallocate();
    xdum.deallocate();
    iold.deallocate();
    nxtseg.deallocate();
    intrfc.deallocate();
    ipatch.deallocate();
    ni.deallocate();
    nj.deallocate();
    nk.deallocate();
    xdmold.deallocate();
    idmold.deallocate();
    nseg1.deallocate();
    nsilo.deallocate();
    nsihi.deallocate();
    nsjlo.deallocate();
    nsjhi.deallocate();
    nsklo.deallocate();
    nskhi.deallocate();
    nb1.deallocate();
    ne1.deallocate();
    nb2.deallocate();
    ne2.deallocate();
    ibct.deallocate();
    bcval.deallocate();
    iovrlp.deallocate();
    ifoflg.deallocate();
    ndat.deallocate();
    xmap.deallocate();
    twotref.deallocate();
    itrb1.deallocate();
    itrb2.deallocate();
    jtrb1.deallocate();
    jtrb2.deallocate();
    ktrb1.deallocate();
    ktrb2.deallocate();
    iturbb.deallocate();
    ibif1.deallocate();
    ibif2.deallocate();
    nbi1.deallocate();
    nei1.deallocate();
    nbj1.deallocate();
    nej1.deallocate();
    nbk1.deallocate();
    nek1.deallocate();
    nbi2.deallocate();
    nei2.deallocate();
    nbj2.deallocate();
    nej2.deallocate();
    nbk2.deallocate();
    nek2.deallocate();
    nd11.deallocate();
    nd21.deallocate();
    nd12.deallocate();
    nd22.deallocate();
    ifsor.deallocate();
    nb1s.deallocate();
    ne1s.deallocate();
    nb2s.deallocate();
    ne2s.deallocate();
    nsubbl.deallocate();
    idobl.deallocate();
    nseg.deallocate();
    idno.deallocate();
    ijk.deallocate();
    idseg.deallocate();
    idnext.deallocate();

    printf("SPLITTING/CONVERSION COMPLETED\n");

    return;
}

} // namespace splitter_ns
