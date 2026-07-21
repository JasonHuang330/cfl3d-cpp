// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "tlns3d_to_cfl3d.h"
#include "runtime/fortran_io.h"
#include "termn8.h"
#include "readkey.h"
#include "getdhdr.h"
#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>

namespace tlns3d_to_cfl3d_ns {

void tlns3d_to_cfl3d(int& iform, int& iver, int& ipatch, int& ironver, int& ngrid, FortranArray1DRef<int> iln, FortranArray1DRef<int> jln, FortranArray1DRef<int> kln, int& ipar, int& nnodes, int& isd, FortranArray2DRef<int> ifrom, FortranArray2DRef<int> xif1, FortranArray2DRef<int> xif2, FortranArray2DRef<int> etf1, FortranArray2DRef<int> etf2, FortranArray1DRef<int> ito, FortranArray1DRef<int> xit1, FortranArray1DRef<int> xit2, FortranArray1DRef<int> ett1, FortranArray1DRef<int> ett2, FortranArray1DRef<int> nfb, FortranArray2DRef<int> iredundant, FortranArray1DRef<int> iconcat, FortranArray3DRef<double> xmap, int& msegn, int& mbloc, FortranArray1DRef<int> iovrlp, FortranArray1DRef<int> nseg, FortranArray3DRef<int> imap, FortranArray2DRef<int> ivisb, FortranArray2DRef<double> twotref, FortranArray1DRef<int> itrb1, FortranArray1DRef<int> itrb2, FortranArray1DRef<int> jtrb1, FortranArray1DRef<int> jtrb2, FortranArray1DRef<int> ktrb1, FortranArray1DRef<int> ktrb2, FortranArray1DRef<int> iturbb, FortranArray2DRef<double> cq, FortranArray2DRef<double> cporous, FortranArray2DRef<double> cpchamber, FortranArray2DRef<double> cradiation, FortranArray1DRef<int> llimit, FortranArray1DRef<int> iitmax, FortranArray1DRef<int> mmcxie, FortranArray1DRef<int> mmceta, FortranArray1DRef<int> iifit, FortranArray1DRef<int> iic0, FortranArray1DRef<int> iiorph, FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray2DRef<int> nblk, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, FortranArray2DRef<double> rkap0g, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<int> idiagg, FortranArray1DRef<int> nblcg, FortranArray2DRef<int> idegg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray2DRef<int> mit, FortranArray1DRef<int> ilamlog, FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamlog, FortranArray1DRef<int> klamhig, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> iforce, int& intmax, int& nsub1, int& msegt, int& maxseg, int& mxbli, int& nrotat, int& ntrans, double& tlref, double& rlref, int& ioflag, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    // Local variables
    char datahdr_storage[10][10] = {};
    FortranArray1DRef<char[10]> datahdr(datahdr_storage, 10, 1);

    // Local scalars
    int iold, nbloc, ibloc, iseg_var, n, m;
    int ibloc1, isegg;
    float rkapi, rkapj, rkapk;
    int irkapi, irkapj, irkapk;
    int ivv;
    int iflag1, nbli, ninter;
    int ifo1, ifo2, jfo1, jfo2, kfo1, kfo2;
    int ifo, jfo, kfo;
    int ibctt, ibctyp, ift, i1tmn, i1tmx, i2tmn, i2tmx;
    int ndata, ns, mm_var;
    int ifs, i1smn, i1smx, i2smn, i2smx;
    int nbt, nbs;
    int ibt, ibs, mseg_loc;
    int mft, m1tmn, m1tmx, m2tmn, m2tmx, mfs, m1smn, m1smx, m2smn, m2smx;
    int mbt, mbs;
    int i1tmn_1, i1tmx_1, i2tmn_1, i2tmx_1;
    int i1smn_1, i1smx_1, i2smn_1, i2smx_1;
    int m1tmn_1, m1tmx_1, m2tmn_1, m2tmx_1;
    int m1smn_1, m1smx_1, m2smn_1, m2smx_1;
    int iflag, itemp;
    int irange1, irange2, jrange1, jrange2, krange1, krange2;
    int llrange11, llrange12, llrange21, llrange22;
    int nbtlast, iftlast, i1tmnlast, i1tmxlast, i2tmnlast, i2tmxlast;
    int ichk_var, ito1, ito2, ifrom1, ifrom2, nfbb;
    int idebug, ninter0, ninter1;
    double xit1n, xit2n, ett1n, ett2n;
    int iton, itom;
    double xit1m, xit2m, ett1m, ett2m;
    int nfbm, nf_var;
    int ncgmax2, iflg, nnnn, nnn, igrid_var;
    int jtest1, jtest2, ktest1, ktest2, itest1, itest2;
    int icgflg, icgflg2, igrd;
    double twall;
    int ngridout, nskip, iunit5sav, ititr, myid_var;
    int i2dd, inewg, igridc, is_var, js_var, ks_var, ie_var, je_var, ke_var;
    int movie, ncs;
    int itrans, irotat;
    double rfreqt, xmag, ymag, zmag, dxmax, dymax, dzmax;
    double rfreqr, thxmag, thymag, thzmag, xorig, yorig, zorig;
    double thxmax, thymax, thzmax;
    int ibcfile;
    int in_var;
    int iseg;

    // Aliases for COMMON block members
    int& iunit5 = cmn_unit5.iunit5;
    char* gridin    = cmn_cflfiles.gridin;
    char* plt3dg    = cmn_cflfiles.plt3dg;
    char* plt3dq    = cmn_cflfiles.plt3dq;
    char* output    = cmn_cflfiles.output;
    char* resid     = cmn_cflfiles.resid;
    char* turbres   = cmn_cflfiles.turbres;
    char* blomx     = cmn_cflfiles.blomx;
    char* output2   = cmn_cflfiles.output2;
    char* printout  = cmn_cflfiles.printout;
    char* pplunge   = cmn_cflfiles.pplunge;
    char* ovrlap    = cmn_cflfiles.ovrlap;
    char* patch     = cmn_cflfiles.patch;
    char* restrt    = cmn_cflfiles.restrt;
    char* gridout   = cmn_cflfiles.gridout;
    char* sdgridin  = cmn_cflfiles.sdgridin;
    char* sdgridout = cmn_cflfiles.sdgridout;
    char* dovrlap   = cmn_cflfiles.dovrlap;
    char* dpatch    = cmn_cflfiles.dpatch;
    char* dresid    = cmn_cflfiles.dresid;
    char* rout      = cmn_ronfiles.rout;
    float* titleron = cmn_ron.titleron;
    int& ncgmax     = cmn_chk1.ncgmax;
    float* title    = cmn_info.title;
    float* rkap     = cmn_info.rkap;
    float& xmach    = cmn_info.xmach;
    float& alpha    = cmn_info.alpha;
    float& beta     = cmn_info.beta;
    float& dt       = cmn_info.dt;
    float& fmax     = cmn_info.fmax;
    int& nit        = cmn_info.nit;
    int& ntt        = cmn_info.ntt;
    int* idiag      = cmn_info.idiag;
    int& nitfo      = cmn_info.nitfo;
    int& iflagts    = cmn_info.iflagts;
    int* iflim      = cmn_info.iflim;
    int& nres       = cmn_info.nres;
    int* levelb     = cmn_info.levelb;
    int& mgflag     = cmn_info.mgflag;
    int& iconsf     = cmn_info.iconsf;
    int& mseq       = cmn_info.mseq;
    int* ncyc1      = cmn_info.ncyc1;
    int* levelt     = cmn_info.levelt;
    int* nitfo1     = cmn_info.nitfo1;
    int& ngam       = cmn_info.ngam;
    int* nsm        = cmn_info.nsm;
    int& iipv       = cmn_info.iipv;
    float& tinf     = cmn_info2.tinf;
    float& reue     = cmn_info2.reue;
    float& c2spe    = cmn_info2.c2spe;
    float& bref     = cmn_info2.bref;
    float& cref     = cmn_info2.cref;
    float& sref     = cmn_info2.sref;
    float& xmc      = cmn_info2.xmc;
    float& ymc      = cmn_info2.ymc;
    float& zmc      = cmn_info2.zmc;
    float& rfreq    = cmn_info2.rfreq;
    float& alphau   = cmn_info2.alphau;
    float& cloc     = cmn_info2.cloc;
    float& cfltau   = cmn_info2.cfltau;
    int& ncg        = cmn_info3.ncg;
    int& isnd       = cmn_info3.isnd;
    int& ialph      = cmn_info3.ialph;
    int& irest      = cmn_info3.irest;
    int& iunst      = cmn_info3.iunst;
    int& ntstep     = cmn_info3.ntstep;
    int& ita        = cmn_info3.ita;
    int& ihstry     = cmn_info3.ihstry;
    int& nplot3d    = cmn_info3.nplot3d;
    int& nprint     = cmn_info3.nprint;
    int& nwrest     = cmn_info3.nwrest;
    int& ichk       = cmn_info3.ichk;
    int& i2d        = cmn_info3.i2d;
    int* mglev      = cmn_info3.mglev;
    int* nem        = cmn_info3.nem;
    int* mitl       = cmn_info3.mitl;
    int& mtt        = cmn_info3.mtt;
    int& ndv        = cmn_info3.ndv;
    int& isdform    = cmn_info3.isdform;
    int& ip3dgrd    = cmn_info3.ip3dgrd;
    int& iplt3dtyp  = cmn_info3.iplt3dtyp;
    float* epsssc   = cmn_mgv.epsssc;
    float* epsssr   = cmn_mgv.epsssr;
    int& issc       = cmn_mgv.issc;
    int& issr       = cmn_mgv.issr;
    float& cltarg   = cmn_alphait.cltarg;
    float& resupdt  = cmn_alphait.resupdt;
    int& ialphit    = cmn_alphait.ialphit;
    float& cprec    = cmn_precon1.cprec;
    float& uref     = cmn_precon1.uref;
    float& avn      = cmn_precon1.avn;
    int& iprecon    = cmn_precon1.iprecon;

    // set idimg,jdimg,kdim to the new (split) values
    for (n = 1; n <= ngrid; n++) {
        idimg(n) = iln(n);
        jdimg(n) = jln(n);
        kdimg(n) = kln(n);
    }
    nbloc = ngrid;

    nbloc = std::abs(nbloc);
    if (nbloc > mbloc) {
        printf("stopping: parameter mbloc must be at least %d\n", nbloc);
        { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
    }

    // always use map file with third line (iold=1)
    iold = 1;

    // read in tlns3d map file
    rewind(fortran_get_unit(8));
    fortran_read_list(8);
    fortran_read_list(8, &nbloc);
    // DO 600 ibloc=1,nbloc / DO 600 iseg=1,nseg(ibloc)
    for (ibloc = 1; ibloc <= nbloc; ibloc++) {
        fortran_read_list(8);
        fortran_read_list(8);
        {
            int _nseg_ibloc, _ivisb4, _ivisb5, _ivisb6, _ivisb7, _ivisb8, _ivisb9, _ivisb10;
            int _iturbb, _itrb1, _itrb2, _jtrb1, _jtrb2, _ktrb1, _ktrb2;
            fortran_read_list(8, &_nseg_ibloc,
                &_ivisb4, &_ivisb5, &_ivisb6,
                &_iturbb, &_itrb1, &_itrb2, &_jtrb1, &_jtrb2, &_ktrb1, &_ktrb2);
            nseg(ibloc)      = _nseg_ibloc;
            ivisb(1,ibloc)   = _ivisb4;
            ivisb(2,ibloc)   = _ivisb5;
            ivisb(3,ibloc)   = _ivisb6;
            iturbb(ibloc)    = _iturbb;
            itrb1(ibloc)     = _itrb1;
            itrb2(ibloc)     = _itrb2;
            jtrb1(ibloc)     = _jtrb1;
            jtrb2(ibloc)     = _jtrb2;
            ktrb1(ibloc)     = _ktrb1;
            ktrb2(ibloc)     = _ktrb2;
        }
        fortran_read_list(8);
        {
            int _v4,_v5,_v6,_v7,_v8,_v9,_v10,_v11,_v12,_v13,_v14,_v15;
            fortran_read_list(8, &_v4,&_v5,&_v6,&_v7,&_v8,&_v9,&_v10,&_v11,&_v12,&_v13,&_v14,&_v15);
            ivisb(4,ibloc)=_v4; ivisb(5,ibloc)=_v5; ivisb(6,ibloc)=_v6;
            ivisb(7,ibloc)=_v7; ivisb(8,ibloc)=_v8; ivisb(9,ibloc)=_v9;
            ivisb(10,ibloc)=_v10; ivisb(11,ibloc)=_v11; ivisb(12,ibloc)=_v12;
            ivisb(13,ibloc)=_v13; ivisb(14,ibloc)=_v14; ivisb(15,ibloc)=_v15;
        }
        {
            int _v16,_v17,_v18,_v19,_v20,_v21,_v22,_v23,_v24,_v25,_v26,_v27;
            fortran_read_list(8, &_v16,&_v17,&_v18,&_v19,&_v20,&_v21,&_v22,&_v23,&_v24,&_v25,&_v26,&_v27);
            ivisb(16,ibloc)=_v16; ivisb(17,ibloc)=_v17; ivisb(18,ibloc)=_v18;
            ivisb(19,ibloc)=_v19; ivisb(20,ibloc)=_v20; ivisb(21,ibloc)=_v21;
            ivisb(22,ibloc)=_v22; ivisb(23,ibloc)=_v23; ivisb(24,ibloc)=_v24;
            ivisb(25,ibloc)=_v25; ivisb(26,ibloc)=_v26; ivisb(27,ibloc)=_v27;
        }
        {
            int _v28,_v29,_v30,_v31,_v32,_v33,_v34;
            fortran_read_list(8, &_v28,&_v29,&_v30,&_v31,&_v32,&_v33,&_v34);
            ivisb(28,ibloc)=_v28; ivisb(29,ibloc)=_v29; ivisb(30,ibloc)=_v30;
            ivisb(31,ibloc)=_v31; ivisb(32,ibloc)=_v32; ivisb(33,ibloc)=_v33;
            ivisb(34,ibloc)=_v34;
        }
        fortran_read_list(8);
        fortran_read_list(8);
        fortran_read_list(8);

        for (iseg = 1; iseg <= nseg(ibloc); iseg++) {
            if (iold == 0) {
                fortran_read_list(8);
                {
                    int _im1,_im2,_im3,_im4,_im5,_im6,_im7,_im8,_im9,_im10,_im11,_im12;
                    fortran_read_list(8, &ibloc1, &isegg,
                        &_im1,&_im2,&_im3,&_im4,&_im5,&_im6,
                        &_im7,&_im8,&_im9,&_im10,&_im11,&_im12);
                    imap(1,iseg,ibloc)=_im1; imap(2,iseg,ibloc)=_im2;
                    imap(3,iseg,ibloc)=_im3; imap(4,iseg,ibloc)=_im4;
                    imap(5,iseg,ibloc)=_im5; imap(6,iseg,ibloc)=_im6;
                    imap(7,iseg,ibloc)=_im7; imap(8,iseg,ibloc)=_im8;
                    imap(9,iseg,ibloc)=_im9; imap(10,iseg,ibloc)=_im10;
                    imap(11,iseg,ibloc)=_im11; imap(12,iseg,ibloc)=_im12;
                }
                {
                    int _im13,_im14,_im15,_im16,_im17,_im18,_im19,_im20;
                    double _twotref;
                    fortran_read_list(8, &_im13,&_im14,&_im15,&_im16,
                        &_im17,&_im18,&_im19,&_im20, &_twotref);
                    imap(13,iseg,ibloc)=_im13; imap(14,iseg,ibloc)=_im14;
                    imap(15,iseg,ibloc)=_im15; imap(16,iseg,ibloc)=_im16;
                    imap(17,iseg,ibloc)=_im17; imap(18,iseg,ibloc)=_im18;
                    imap(19,iseg,ibloc)=_im19; imap(20,iseg,ibloc)=_im20;
                    twotref(iseg,ibloc) = _twotref;
                }
            } else {
                fortran_read_list(8);
                {
                    int _im1,_im2,_im3,_im4,_im5,_im6,_im7,_im8,_im9,_im10,_im11,_im12;
                    fortran_read_list(8, &ibloc1, &isegg,
                        &_im1,&_im2,&_im3,&_im4,&_im5,&_im6,
                        &_im7,&_im8,&_im9,&_im10,&_im11,&_im12);
                    imap(1,iseg,ibloc)=_im1; imap(2,iseg,ibloc)=_im2;
                    imap(3,iseg,ibloc)=_im3; imap(4,iseg,ibloc)=_im4;
                    imap(5,iseg,ibloc)=_im5; imap(6,iseg,ibloc)=_im6;
                    imap(7,iseg,ibloc)=_im7; imap(8,iseg,ibloc)=_im8;
                    imap(9,iseg,ibloc)=_im9; imap(10,iseg,ibloc)=_im10;
                    imap(11,iseg,ibloc)=_im11; imap(12,iseg,ibloc)=_im12;
                }
                {
                    int _im13,_im14,_im15,_im16,_im17,_im18,_im19,_im20;
                    fortran_read_list(8, &_im13,&_im14,&_im15,&_im16,
                        &_im17,&_im18,&_im19,&_im20);
                    imap(13,iseg,ibloc)=_im13; imap(14,iseg,ibloc)=_im14;
                    imap(15,iseg,ibloc)=_im15; imap(16,iseg,ibloc)=_im16;
                    imap(17,iseg,ibloc)=_im17; imap(18,iseg,ibloc)=_im18;
                    imap(19,iseg,ibloc)=_im19; imap(20,iseg,ibloc)=_im20;
                }
            }

            if (iold > 0) {
                double _twotref, _cq, _cporous, _cpchamber, _cradiation;
                fortran_read_list(8, &_twotref, &_cq, &_cporous, &_cpchamber, &_cradiation);
                twotref(iseg,ibloc)   = _twotref;
                cq(iseg,ibloc)        = _cq;
                cporous(iseg,ibloc)   = _cporous;
                cpchamber(iseg,ibloc) = _cpchamber;
                cradiation(iseg,ibloc)= _cradiation;
            }

            // additional data for cfl3d
            fortran_read_list(8);
            {
                double _x1,_x2,_x3,_x4,_x5;
                fortran_read_list(8, &_x1,&_x2,&_x3,&_x4,&_x5);
                xmap(1,iseg,ibloc)=_x1; xmap(2,iseg,ibloc)=_x2;
                xmap(3,iseg,ibloc)=_x3; xmap(4,iseg,ibloc)=_x4;
                xmap(5,iseg,ibloc)=_x5;
            }
            {
                double _x6,_x7,_x8,_x9,_x10;
                fortran_read_list(8, &_x6,&_x7,&_x8,&_x9,&_x10);
                xmap(6,iseg,ibloc)=_x6; xmap(7,iseg,ibloc)=_x7;
                xmap(8,iseg,ibloc)=_x8; xmap(9,iseg,ibloc)=_x9;
                xmap(10,iseg,ibloc)=_x10;
            }
            {
                double _x11,_x12,_x13,_x14,_x15;
                fortran_read_list(8, &_x11,&_x12,&_x13,&_x14,&_x15);
                xmap(11,iseg,ibloc)=_x11; xmap(12,iseg,ibloc)=_x12;
                xmap(13,iseg,ibloc)=_x13; xmap(14,iseg,ibloc)=_x14;
                xmap(15,iseg,ibloc)=_x15;
            }
            {
                double _x16,_x17,_x18,_x19,_x20;
                fortran_read_list(8, &_x16,&_x17,&_x18,&_x19,&_x20);
                xmap(16,iseg,ibloc)=_x16; xmap(17,iseg,ibloc)=_x17;
                xmap(18,iseg,ibloc)=_x18; xmap(19,iseg,ibloc)=_x19;
                xmap(20,iseg,ibloc)=_x20;
            }
            {
                double _x21,_x22,_x23,_x24,_x25;
                fortran_read_list(8, &_x21,&_x22,&_x23,&_x24,&_x25);
                xmap(21,iseg,ibloc)=_x21; xmap(22,iseg,ibloc)=_x22;
                xmap(23,iseg,ibloc)=_x23; xmap(24,iseg,ibloc)=_x24;
                xmap(25,iseg,ibloc)=_x25;
            }
            {
                double _x26;
                fortran_read_list(8, &_x26);
                xmap(26,iseg,ibloc)=_x26;
            }
            {
                int _im21,_im22,_im23;
                fortran_read_list(8, &_im21,&_im22,&_im23);
                imap(21,iseg,ibloc)=_im21;
                imap(22,iseg,ibloc)=_im22;
                imap(23,iseg,ibloc)=_im23;
            }
            fortran_read_list(8);
        } // end iseg loop
    } // end ibloc loop (600 continue)



    // translate tlns3d mapfile data to cfl3d data
    for (n = 1; n <= ngrid; n++) {
        iviscg(n,1) = ivisb(1,n);
        iviscg(n,2) = ivisb(2,n);
        iviscg(n,3) = ivisb(3,n);

        // ivv = 0 if there are no viscous blocks
        ivv = 1;
        if (iviscg(n,1) == 0 && iviscg(n,2) == 0 && iviscg(n,3) == 0)
            ivv = 0;
    }

    iflag1 = 0;
    nbli   = 0;
    ninter = 0;

    for (n = 1; n <= ngrid; n++) {
        for (iseg = 1; iseg <= nseg(n); iseg++) {
            iredundant(n,iseg) = 0;
        }
    }

    // DO 1000 n=1,ngrid
    for (n = 1; n <= ngrid; n++) {
        // some items carried in the ivisb array are independent
        // of whether the original block gets split or not
        ncgg(n)     = ivisb( 4,n);
        iemg(n)     = ivisb( 5,n);
        iadvance(n) = ivisb( 6,n);
        iforce(n)   = ivisb( 7,n);
        ilamlog(n)  = ivisb( 8,n);
        ilamhig(n)  = ivisb( 9,n);
        jlamlog(n)  = ivisb(10,n);
        jlamhig(n)  = ivisb(11,n);
        klamlog(n)  = ivisb(12,n);
        klamhig(n)  = ivisb(13,n);
        inewgg(n)   = ivisb(14,n);
        igridg(n)   = ivisb(15,n);
        isg(n)      = ivisb(16,n);
        jsg(n)      = ivisb(17,n);
        ksg(n)      = ivisb(18,n);
        ieg(n)      = ivisb(19,n);
        jeg(n)      = ivisb(20,n);
        keg(n)      = ivisb(21,n);
        idiagg(n,1) = ivisb(22,n);
        idiagg(n,2) = ivisb(23,n);
        idiagg(n,3) = ivisb(24,n);
        iflimg(n,1) = ivisb(25,n);
        iflimg(n,2) = ivisb(26,n);
        iflimg(n,3) = ivisb(27,n);
        ifdsg(n,1)  = ivisb(28,n);
        ifdsg(n,2)  = ivisb(29,n);
        ifdsg(n,3)  = ivisb(30,n);
        irkapi      = ivisb(31,n);
        irkapj      = ivisb(32,n);
        irkapk      = ivisb(33,n);
        if (irkapi == -1) rkapi = -1.0f;
        if (irkapi ==  0) rkapi = 0.0f;
        if (irkapi ==  1) rkapi = 1.0f;
        if (irkapi ==  3) rkapi = 0.3333f;
        if (irkapj == -1) rkapj = -1.0f;
        if (irkapj ==  0) rkapj = 0.0f;
        if (irkapj ==  1) rkapj = 1.0f;
        if (irkapj ==  3) rkapj = 0.3333f;
        if (irkapk == -1) rkapk = -1.0f;
        if (irkapk ==  0) rkapk = 0.0f;
        if (irkapk ==  1) rkapk = 1.0f;
        if (irkapk ==  3) rkapk = 0.3333f;
        rkap0g(n,1) = rkapi;
        rkap0g(n,2) = rkapj;
        rkap0g(n,3) = rkapk;
        iovrlp(n)   = ivisb(34,n);

        // questionable data now overwritten
        iforce(n)  = 0;
        iemg(n)    = 0;
        isg(n)     = 0;
        jsg(n)     = 0;
        ksg(n)     = 0;
        ieg(n)     = 0;
        jeg(n)     = 0;
        keg(n)     = 0;
        igridg(n)  = 0;
        ilamlog(n) = 0;
        ilamhig(n) = 0;
        jlamlog(n) = 0;
        jlamhig(n) = 0;
        klamlog(n) = 0;
        klamhig(n) = 0;

        // initialize number of bc segments and force flags
        nbci0(n)   = 0;
        nbcidim(n) = 0;
        nbcj0(n)   = 0;
        nbcjdim(n) = 0;
        nbck0(n)   = 0;
        nbckdim(n) = 0;
        ifo1       = 0;
        ifo2       = 0;
        jfo1       = 0;
        jfo2       = 0;
        kfo1       = 0;
        kfo2       = 0;

        // ivv for this block
        ivv = 1;
        if (iviscg(n,1) == 0 && iviscg(n,2) == 0 && iviscg(n,3) == 0)
            ivv = 0;



        // DO 1010 iseg=1,nseg(n)
        for (iseg = 1; iseg <= nseg(n); iseg++) {
            ibctt   = imap(1,iseg,n);
            if (ibctt == 0 || ibctt == 1 || ibctt == -1) {
                ibctyp = 0;
            } else {
                ibctyp = ibctt;
            }

            ift   = imap(2,iseg,n);
            i1tmn = imap(3,iseg,n);
            i1tmx = imap(4,iseg,n);
            i2tmn = imap(5,iseg,n);
            i2tmx = imap(6,iseg,n);

            // i0 boundary data
            if (ift == 1) {
                // don't include duplicate bc data
                if (nbci0(n) > 0 && ibctyp == 0) {
                    m = 1;
                    if ((ibctyp == ibcinfo(n,nbci0(n),1,m)) &&
                        (i1tmn  == ibcinfo(n,nbci0(n),2,m)) &&
                        (i1tmx  == ibcinfo(n,nbci0(n),3,m)) &&
                        (i2tmn  == ibcinfo(n,nbci0(n),4,m)) &&
                        (i2tmx  == ibcinfo(n,nbci0(n),5,m))) {
                        goto label_6001;
                    }
                }
                nbci0(n) = nbci0(n) + 1;
                ns = nbci0(n);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
                m  = 1;
                ibcinfo(n,ns,1,m) = ibctyp;
                ibcinfo(n,ns,2,m) = i1tmn;
                ibcinfo(n,ns,3,m) = i1tmx;
                ibcinfo(n,ns,4,m) = i2tmn;
                ibcinfo(n,ns,5,m) = i2tmx;
                ibcinfo(n,ns,6,m) = imap(22,iseg,n);
                if (ivv == 0) {
                    if (ibctyp == 1005 || ibctyp == 1006) {
                        ifo1 = 1;
                    }
                } else {
                    if (ibctyp == 1004 || std::abs(ibctyp) == 2004 ||
                        std::abs(ibctyp) == 2014) {
                        ifo1 = 1;
                    }
                }
                ndata = imap(21,iseg,n);
                ibcinfo(n,ns,7,m) = ndata;
                if (std::abs(ndata) > 0) {
                    for (mm_var = 1; mm_var <= std::abs(ndata); mm_var++) {
                        bcvali(n,ns,mm_var,m) = xmap(mm_var,iseg,n);
                    }
                }
            }
            label_6001:;

            // idim boundary data
            if (ift == 2) {
                if (nbcidim(n) > 0 && ibctyp == 0) {
                    m = 2;
                    if ((ibctyp == ibcinfo(n,nbcidim(n),1,m)) &&
                        (i1tmn  == ibcinfo(n,nbcidim(n),2,m)) &&
                        (i1tmx  == ibcinfo(n,nbcidim(n),3,m)) &&
                        (i2tmn  == ibcinfo(n,nbcidim(n),4,m)) &&
                        (i2tmx  == ibcinfo(n,nbcidim(n),5,m))) {
                        goto label_6002;
                    }
                }
                nbcidim(n) = nbcidim(n) + 1;
                ns = nbcidim(n);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
                m  = 2;
                ibcinfo(n,ns,1,m) = ibctyp;
                ibcinfo(n,ns,2,m) = i1tmn;
                ibcinfo(n,ns,3,m) = i1tmx;
                ibcinfo(n,ns,4,m) = i2tmn;
                ibcinfo(n,ns,5,m) = i2tmx;
                ibcinfo(n,ns,6,m) = imap(22,iseg,n);
                if (ivv == 0) {
                    if (ibctyp == 1005 || ibctyp == 1006) {
                        ifo2 = 1;
                    }
                } else {
                    if (ibctyp == 1004 || std::abs(ibctyp) == 2004 ||
                        std::abs(ibctyp) == 2014) {
                        ifo2 = 1;
                    }
                }
                ndata = imap(21,iseg,n);
                ibcinfo(n,ns,7,m) = ndata;
                if (std::abs(ndata) > 0) {
                    for (mm_var = 1; mm_var <= std::abs(ndata); mm_var++) {
                        bcvali(n,ns,mm_var,m) = xmap(mm_var,iseg,n);
                    }
                }
            }
            label_6002:;

            // j0 boundary data
            if (ift == 3) {
                if (nbcj0(n) > 0 && ibctyp == 0) {
                    m = 1;
                    if ((ibctyp == jbcinfo(n,nbcj0(n),1,m)) &&
                        (i1tmn  == jbcinfo(n,nbcj0(n),2,m)) &&
                        (i1tmx  == jbcinfo(n,nbcj0(n),3,m)) &&
                        (i2tmn  == jbcinfo(n,nbcj0(n),4,m)) &&
                        (i2tmx  == jbcinfo(n,nbcj0(n),5,m))) {
                        goto label_6003;
                    }
                }
                nbcj0(n) = nbcj0(n) + 1;
                ns = nbcj0(n);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
                m  = 1;
                jbcinfo(n,ns,1,m) = ibctyp;
                jbcinfo(n,ns,2,m) = i2tmn;
                jbcinfo(n,ns,3,m) = i2tmx;
                jbcinfo(n,ns,4,m) = i1tmn;
                jbcinfo(n,ns,5,m) = i1tmx;
                jbcinfo(n,ns,6,m) = imap(22,iseg,n);
                if (ivv == 0) {
                    if (ibctyp == 1005 || ibctyp == 1006) {
                        jfo1 = 1;
                    }
                } else {
                    if (ibctyp == 1004 || std::abs(ibctyp) == 2004 ||
                        std::abs(ibctyp) == 2014) {
                        jfo1 = 1;
                    }
                }
                ndata = imap(21,iseg,n);
                jbcinfo(n,ns,7,m) = ndata;
                if (std::abs(ndata) > 0) {
                    for (mm_var = 1; mm_var <= std::abs(ndata); mm_var++) {
                        bcvalj(n,ns,mm_var,m) = xmap(mm_var,iseg,n);
                    }
                }
            }
            label_6003:;

            // jdim boundary data
            if (ift == 4) {
                if (nbcjdim(n) > 0 && ibctyp == 0) {
                    m = 2;
                    if ((ibctyp == jbcinfo(n,nbcjdim(n),1,m)) &&
                        (i1tmn  == jbcinfo(n,nbcjdim(n),2,m)) &&
                        (i1tmx  == jbcinfo(n,nbcjdim(n),3,m)) &&
                        (i2tmn  == jbcinfo(n,nbcjdim(n),4,m)) &&
                        (i2tmx  == jbcinfo(n,nbcjdim(n),5,m))) {
                        goto label_6004;
                    }
                }
                nbcjdim(n) = nbcjdim(n) + 1;
                ns = nbcjdim(n);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
                m  = 2;
                jbcinfo(n,ns,1,m) = ibctyp;
                jbcinfo(n,ns,2,m) = i2tmn;
                jbcinfo(n,ns,3,m) = i2tmx;
                jbcinfo(n,ns,4,m) = i1tmn;
                jbcinfo(n,ns,5,m) = i1tmx;
                jbcinfo(n,ns,6,m) = imap(22,iseg,n);
                if (ivv == 0) {
                    if (ibctyp == 1005 || ibctyp == 1006) {
                        jfo2 = 1;
                    }
                } else {
                    if (ibctyp == 1004 || std::abs(ibctyp) == 2004 ||
                        std::abs(ibctyp) == 2014) {
                        jfo2 = 1;
                    }
                }
                ndata = imap(21,iseg,n);
                jbcinfo(n,ns,7,m) = ndata;
                if (std::abs(ndata) > 0) {
                    for (mm_var = 1; mm_var <= std::abs(ndata); mm_var++) {
                        bcvalj(n,ns,mm_var,m) = xmap(mm_var,iseg,n);
                    }
                }
            }
            label_6004:;



            // k0 boundary data
            if (ift == 5) {
                if (nbck0(n) > 0 && ibctyp == 0) {
                    m = 1;
                    if ((ibctyp == kbcinfo(n,nbck0(n),1,m)) &&
                        (i1tmn  == kbcinfo(n,nbck0(n),2,m)) &&
                        (i1tmx  == kbcinfo(n,nbck0(n),3,m)) &&
                        (i2tmn  == kbcinfo(n,nbck0(n),4,m)) &&
                        (i2tmx  == kbcinfo(n,nbck0(n),5,m))) {
                        goto label_6005;
                    }
                }
                nbck0(n) = nbck0(n) + 1;
                ns = nbck0(n);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
                m  = 1;
                kbcinfo(n,ns,1,m) = ibctyp;
                kbcinfo(n,ns,2,m) = i1tmn;
                kbcinfo(n,ns,3,m) = i1tmx;
                kbcinfo(n,ns,4,m) = i2tmn;
                kbcinfo(n,ns,5,m) = i2tmx;
                kbcinfo(n,ns,6,m) = imap(22,iseg,n);
                if (ivv == 0) {
                    if (ibctyp == 1005 || ibctyp == 1006) {
                        kfo1 = 1;
                    }
                } else {
                    if (ibctyp == 1004 || std::abs(ibctyp) == 2004 ||
                        std::abs(ibctyp) == 2014) {
                        kfo1 = 1;
                    }
                }
                ndata = imap(21,iseg,n);
                kbcinfo(n,ns,7,m) = ndata;
                if (std::abs(ndata) > 0) {
                    for (mm_var = 1; mm_var <= std::abs(ndata); mm_var++) {
                        bcvalk(n,ns,mm_var,m) = xmap(mm_var,iseg,n);
                    }
                }
            }
            label_6005:;

            // kdim boundary data
            if (ift == 6) {
                if (nbckdim(n) > 0 && ibctyp == 0) {
                    m = 2;
                    if ((ibctyp == kbcinfo(n,nbckdim(n),1,m)) &&
                        (i1tmn  == kbcinfo(n,nbckdim(n),2,m)) &&
                        (i1tmx  == kbcinfo(n,nbckdim(n),3,m)) &&
                        (i2tmn  == kbcinfo(n,nbckdim(n),4,m)) &&
                        (i2tmx  == kbcinfo(n,nbckdim(n),5,m))) {
                        goto label_6006;
                    }
                }
                nbckdim(n) = nbckdim(n) + 1;
                ns = nbckdim(n);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
                m  = 2;
                kbcinfo(n,ns,1,m) = ibctyp;
                kbcinfo(n,ns,2,m) = i1tmn;
                kbcinfo(n,ns,3,m) = i1tmx;
                kbcinfo(n,ns,4,m) = i2tmn;
                kbcinfo(n,ns,5,m) = i2tmx;
                kbcinfo(n,ns,6,m) = imap(22,iseg,n);
                if (ivv == 0) {
                    if (ibctyp == 1005 || ibctyp == 1006) {
                        kfo2 = 1;
                    }
                } else {
                    if (ibctyp == 1004 || std::abs(ibctyp) == 2004 ||
                        std::abs(ibctyp) == 2014) {
                        kfo2 = 1;
                    }
                }
                ndata = imap(21,iseg,n);
                kbcinfo(n,ns,7,m) = ndata;
                if (std::abs(ndata) > 0) {
                    for (mm_var = 1; mm_var <= std::abs(ndata); mm_var++) {
                        bcvalk(n,ns,mm_var,m) = xmap(mm_var,iseg,n);
                    }
                }
            }
            label_6006:;

        } // end DO 1010 iseg loop

        // set iforce for this block based on solid surface bc information
        ifo = 0;
        jfo = 0;
        kfo = 0;
        if (ifo1 == 1 && ifo2 == 0) ifo = 1;
        if (ifo1 == 0 && ifo2 == 1) ifo = 2;
        if (ifo1 == 1 && ifo2 == 1) ifo = 3;
        if (jfo1 == 1 && jfo2 == 0) jfo = 1;
        if (jfo1 == 0 && jfo2 == 1) jfo = 2;
        if (jfo1 == 1 && jfo2 == 1) jfo = 3;
        if (kfo1 == 1 && kfo2 == 0) kfo = 1;
        if (kfo1 == 0 && kfo2 == 1) kfo = 2;
        if (kfo1 == 1 && kfo2 == 1) kfo = 3;
        iforce(n) = ifo*100 + jfo*10 + kfo;



        // extract 1-1 data
        // on the first pass, check for matching partners of
        // 1-1 segments and mark one as redundant
        for (iseg = 1; iseg <= nseg(n); iseg++) {
            // t = target ("to") s = source ("from")
            if (iredundant(n,iseg) == 0) {
            if (imap(1,iseg,n) == 0 || imap(1,iseg,n) == 1) {
                ift   = imap(2,iseg,n);
                i1tmn = imap(3,iseg,n);
                i1tmx = imap(4,iseg,n);
                i2tmn = imap(5,iseg,n);
                i2tmx = imap(6,iseg,n);
                ifs   = imap(8,iseg,n);
                i1smn = imap(9,iseg,n);
                i1smx = imap(10,iseg,n);
                i2smn = imap(11,iseg,n);
                i2smx = imap(12,iseg,n);
                ibt = n;
                ibs = imap(7,iseg,n);
                for (mseg_loc = 1; mseg_loc <= nseg(ibs); mseg_loc++) {
                    if (ibs == ibt && mseg_loc == iseg) goto label_1616;
                    if (imap(1,mseg_loc,ibs) == 0 || imap(1,mseg_loc,ibs) == 1
                        && iredundant(ibs,mseg_loc) == 0) {
                        mft   = imap(2,mseg_loc,ibs);
                        m1tmn = imap(3,mseg_loc,ibs);
                        m1tmx = imap(4,mseg_loc,ibs);
                        m2tmn = imap(5,mseg_loc,ibs);
                        m2tmx = imap(6,mseg_loc,ibs);
                        mfs   = imap(8,mseg_loc,ibs);
                        m1smn = imap(9,mseg_loc,ibs);
                        m1smx = imap(10,mseg_loc,ibs);
                        m2smn = imap(11,mseg_loc,ibs);
                        m2smx = imap(12,mseg_loc,ibs);
                        mbt = ibs;
                        mbs = imap(7,mseg_loc,ibs);
                        if (mbt == ibs && mbs == ibt &&
                            mft == std::abs(ifs) && ift == std::abs(mfs)) {
                            i1tmn_1 = i1tmn; i1tmx_1 = i1tmx;
                            if (i1tmn > i1tmx) { i1tmx_1 = i1tmn; i1tmn_1 = i1tmx; }
                            i2tmn_1 = i2tmn; i2tmx_1 = i2tmx;
                            if (i2tmn > i2tmx) { i2tmx_1 = i2tmn; i2tmn_1 = i2tmx; }
                            i1smn_1 = i1smn; i1smx_1 = i1smx;
                            if (i1smn > i1smx) { i1smx_1 = i1smn; i1smn_1 = i1smx; }
                            i2smn_1 = i2smn; i2smx_1 = i2smx;
                            if (i2smn > i2smx) { i2smx_1 = i2smn; i2smn_1 = i2smx; }
                            m1tmn_1 = m1tmn; m1tmx_1 = m1tmx;
                            if (m1tmn > m1tmx) { m1tmx_1 = m1tmn; m1tmn_1 = m1tmx; }
                            m2tmn_1 = m2tmn; m2tmx_1 = m2tmx;
                            if (m2tmn > m2tmx) { m2tmx_1 = m2tmn; m2tmn_1 = m2tmx; }
                            m1smn_1 = m1smn; m1smx_1 = m1smx;
                            if (m1smn > m1smx) { m1smx_1 = m1smn; m1smn_1 = m1smx; }
                            m2smn_1 = m2smn; m2smx_1 = m2smx;
                            if (m2smn > m2smx) { m2smx_1 = m2smn; m2smn_1 = m2smx; }
                            iflag = 0;
                            iflag = iflag + (i1tmn_1 - m1smn_1);
                            iflag = iflag + (i1tmx_1 - m1smx_1);
                            iflag = iflag + (i2tmn_1 - m2smn_1);
                            iflag = iflag + (i2tmx_1 - m2smx_1);
                            if (iflag == 0) {
                                iredundant(ibs,mseg_loc) = 1;
                            }
                        }
                    }
                    label_1616:;
                } // end mseg_loc loop
            }
            } // end iredundant check
        } // end iseg loop (first pass)



        // DO 1020 iseg=1,nseg(n) - set cfl3d 1-1 interface data
        for (iseg = 1; iseg <= nseg(n); iseg++) {
            if (imap(1,iseg,n) == 0 || imap(1,iseg,n) == 1) {
            if (iredundant(n,iseg) == 0) {
                nbli = nbli + 1;
                nblon(nbli) = 0;
                ift   = imap(2,iseg,n);
                i1tmn = imap(3,iseg,n);
                i1tmx = imap(4,iseg,n);
                i2tmn = imap(5,iseg,n);
                i2tmx = imap(6,iseg,n);
                ifs   = imap(8,iseg,n);
                i1smn = imap(9,iseg,n);
                i1smx = imap(10,iseg,n);
                i2smn = imap(11,iseg,n);
                i2smx = imap(12,iseg,n);
                nbt = n;
                nbs = imap(7,iseg,n);
                nblk(1,nbli) = nbt;
                nblk(2,nbli) = nbs;
                if (ift == 1) {
                    limblk(1,1,nbli)=1; limblk(1,4,nbli)=1;
                    limblk(1,2,nbli)=i1tmn; limblk(1,5,nbli)=i1tmx;
                    limblk(1,3,nbli)=i2tmn; limblk(1,6,nbli)=i2tmx;
                    isva(1,1,nbli)=2; isva(1,2,nbli)=3;
                }
                if (ift == 2) {
                    limblk(1,1,nbli)=idimg(nbt); limblk(1,4,nbli)=idimg(nbt);
                    limblk(1,2,nbli)=i1tmn; limblk(1,5,nbli)=i1tmx;
                    limblk(1,3,nbli)=i2tmn; limblk(1,6,nbli)=i2tmx;
                    isva(1,1,nbli)=2; isva(1,2,nbli)=3;
                }
                if (ift == 3) {
                    limblk(1,1,nbli)=i2tmn; limblk(1,4,nbli)=i2tmx;
                    limblk(1,2,nbli)=1; limblk(1,5,nbli)=1;
                    limblk(1,3,nbli)=i1tmn; limblk(1,6,nbli)=i1tmx;
                    isva(1,1,nbli)=1; isva(1,2,nbli)=3;
                }
                if (ift == 4) {
                    limblk(1,1,nbli)=i2tmn; limblk(1,4,nbli)=i2tmx;
                    limblk(1,2,nbli)=jdimg(nbt); limblk(1,5,nbli)=jdimg(nbt);
                    limblk(1,3,nbli)=i1tmn; limblk(1,6,nbli)=i1tmx;
                    isva(1,1,nbli)=1; isva(1,2,nbli)=3;
                }
                if (ift == 5) {
                    limblk(1,1,nbli)=i1tmn; limblk(1,4,nbli)=i1tmx;
                    limblk(1,2,nbli)=i2tmn; limblk(1,5,nbli)=i2tmx;
                    limblk(1,3,nbli)=1; limblk(1,6,nbli)=1;
                    isva(1,1,nbli)=1; isva(1,2,nbli)=2;
                }
                if (ift == 6) {
                    limblk(1,1,nbli)=i1tmn; limblk(1,4,nbli)=i1tmx;
                    limblk(1,2,nbli)=i2tmn; limblk(1,5,nbli)=i2tmx;
                    limblk(1,3,nbli)=kdimg(nbt); limblk(1,6,nbli)=kdimg(nbt);
                    isva(1,1,nbli)=1; isva(1,2,nbli)=2;
                }
                if (std::abs(ifs) == 1) {
                    limblk(2,1,nbli)=1; limblk(2,4,nbli)=1;
                    limblk(2,2,nbli)=i1smn; limblk(2,5,nbli)=i1smx;
                    limblk(2,3,nbli)=i2smn; limblk(2,6,nbli)=i2smx;
                    isva(2,1,nbli)=2; isva(2,2,nbli)=3;
                }
                if (std::abs(ifs) == 2) {
                    limblk(2,1,nbli)=idimg(nbs); limblk(2,4,nbli)=idimg(nbs);
                    limblk(2,2,nbli)=i1smn; limblk(2,5,nbli)=i1smx;
                    limblk(2,3,nbli)=i2smn; limblk(2,6,nbli)=i2smx;
                    isva(2,1,nbli)=2; isva(2,2,nbli)=3;
                }
                if (std::abs(ifs) == 3) {
                    limblk(2,1,nbli)=i2smn; limblk(2,4,nbli)=i2smx;
                    limblk(2,2,nbli)=1; limblk(2,5,nbli)=1;
                    limblk(2,3,nbli)=i1smn; limblk(2,6,nbli)=i1smx;
                    isva(2,1,nbli)=1; isva(2,2,nbli)=3;
                }
                if (std::abs(ifs) == 4) {
                    limblk(2,1,nbli)=i2smn; limblk(2,4,nbli)=i2smx;
                    limblk(2,2,nbli)=jdimg(nbs); limblk(2,5,nbli)=jdimg(nbs);
                    limblk(2,3,nbli)=i1smn; limblk(2,6,nbli)=i1smx;
                    isva(2,1,nbli)=1; isva(2,2,nbli)=3;
                }
                if (std::abs(ifs) == 5) {
                    limblk(2,1,nbli)=i1smn; limblk(2,4,nbli)=i1smx;
                    limblk(2,2,nbli)=i2smn; limblk(2,5,nbli)=i2smx;
                    limblk(2,3,nbli)=1; limblk(2,6,nbli)=1;
                    isva(2,1,nbli)=1; isva(2,2,nbli)=2;
                }
                if (std::abs(ifs) == 6) {
                    limblk(2,1,nbli)=i1smn; limblk(2,4,nbli)=i1smx;
                    limblk(2,2,nbli)=i2smn; limblk(2,5,nbli)=i2smx;
                    limblk(2,3,nbli)=kdimg(nbs); limblk(2,6,nbli)=kdimg(nbs);
                    isva(2,1,nbli)=1; isva(2,2,nbli)=2;
                }
                if (ifs < 0) {
                    itemp = isva(2,1,nbli);
                    isva(2,1,nbli) = isva(2,2,nbli);
                    isva(2,2,nbli) = itemp;
                }
                irange1 = std::abs(limblk(1,4,nbli) - limblk(1,1,nbli));
                irange2 = std::abs(limblk(2,4,nbli) - limblk(2,1,nbli));
                jrange1 = std::abs(limblk(1,5,nbli) - limblk(1,2,nbli));
                jrange2 = std::abs(limblk(2,5,nbli) - limblk(2,2,nbli));
                krange1 = std::abs(limblk(1,6,nbli) - limblk(1,3,nbli));
                krange2 = std::abs(limblk(2,6,nbli) - limblk(2,3,nbli));
                if (isva(1,1,nbli)==1)      llrange11=irange1;
                else if (isva(1,1,nbli)==2) llrange11=jrange1;
                else if (isva(1,1,nbli)==3) llrange11=krange1;
                else                         llrange11=0;
                if (isva(1,2,nbli)==1)      llrange12=irange1;
                else if (isva(1,2,nbli)==2) llrange12=jrange1;
                else if (isva(1,2,nbli)==3) llrange12=krange1;
                else                         llrange12=0;
                if (isva(2,1,nbli)==1)      llrange21=irange2;
                else if (isva(2,1,nbli)==2) llrange21=jrange2;
                else if (isva(2,1,nbli)==3) llrange21=krange2;
                else                         llrange21=0;
                if (isva(2,2,nbli)==1)      llrange22=irange2;
                else if (isva(2,2,nbli)==2) llrange22=jrange2;
                else if (isva(2,2,nbli)==3) llrange22=krange2;
                else                         llrange22=0;
                if (llrange11 != llrange21 || llrange12 != llrange22) {
                    itemp = isva(2,1,nbli);
                    isva(2,1,nbli) = isva(2,2,nbli);
                    isva(2,2,nbli) = itemp;
                }
            }
            } // end iredundant/imap check
        } // end DO 1020 iseg loop



        // extract patched interface data
        // DO 1030 iseg=1,nseg(n)
        nbtlast   = 0;
        iftlast   = 0;
        i1tmnlast = 0;
        i1tmxlast = 0;
        i2tmnlast = 0;
        i2tmxlast = 0;
        for (iseg = 1; iseg <= nseg(n); iseg++) {
            if (imap(1,iseg,n) == -1) {
                ift   = imap(2,iseg,n);
                i1tmn = imap(3,iseg,n);
                i1tmx = imap(4,iseg,n);
                i2tmn = imap(5,iseg,n);
                i2tmx = imap(6,iseg,n);
                ifs   = std::abs(imap(8,iseg,n));
                i1smn = imap(9,iseg,n);
                i1smx = imap(10,iseg,n);
                i2smn = imap(11,iseg,n);
                i2smx = imap(12,iseg,n);
                // reverse order on source zone indices to go from low to high
                if (i1smn > i1smx) { itemp=i1smn; i1smn=i1smx; i1smx=itemp; }
                if (i2smn > i2smx) { itemp=i2smn; i2smn=i2smx; i2smx=itemp; }
                nbt = n;
                nbs = imap(7,iseg,n);
                // check if same target face as last segment
                ichk_var = 0;
                ichk_var += std::abs(nbt   - nbtlast);
                ichk_var += std::abs(ift   - iftlast);
                ichk_var += std::abs(i1tmn - i1tmnlast);
                ichk_var += std::abs(i1tmx - i1tmxlast);
                ichk_var += std::abs(i2tmn - i2tmnlast);
                ichk_var += std::abs(i2tmx - i2tmxlast);
                // target ("to") data
                if (ichk_var != 0) {
                    ninter = ninter + 1;
                    nfb(ninter) = 0;
                    ito1 = nbt;
                    if (ift == 1 || ift == 2) {
                        ito2 = 11;
                        if (ift == 2) ito2 = 12;
                        xit1(ninter) = i1tmn;
                        xit2(ninter) = i1tmx;
                        ett1(ninter) = i2tmn;
                        ett2(ninter) = i2tmx;
                        ito(ninter)  = 100*ito1 + ito2;
                    }
                    if (ift == 3 || ift == 4) {
                        ito2 = 21;
                        if (ift == 4) ito2 = 22;
                        xit1(ninter) = i1tmn;
                        xit2(ninter) = i1tmx;
                        ett1(ninter) = i2tmn;
                        ett2(ninter) = i2tmx;
                        ito(ninter)  = 100*ito1 + ito2;
                    }
                    if (ift == 5 || ift == 6) {
                        ito2 = 31;
                        if (ift == 6) ito2 = 32;
                        xit1(ninter) = i2tmn;
                        xit2(ninter) = i2tmx;
                        ett1(ninter) = i1tmn;
                        ett2(ninter) = i1tmx;
                        ito(ninter)  = 100*ito1 + ito2;
                    }
                }
                // source ("from") data
                nfb(ninter) = nfb(ninter) + 1;
                nfbb = nfb(ninter);
                ifrom1 = nbs;
                if (ifs == 1 || ifs == 2) {
                    ifrom2 = 11;
                    if (ifs == 2) ifrom2 = 12;
                    if (i1smx > jdimg(nbs)) {
                        printf("correcting error in tlns3d map file:\n");
                        printf("i1smx>jdimg(nbs): i1smx,jdimg(nbs) = %d %d\n", i1smx, jdimg(nbs));
                        i1smx = jdimg(nbs);
                    }
                    if (i2smx > kdimg(nbs)) {
                        printf("correcting error in tlns3d map file:\n");
                        printf("i2smx>kdimg(nbs): i2smx,kdimg(nbs) = %d %d\n", i2smx, kdimg(nbs));
                        i2smx = kdimg(nbs);
                    }
                    xif1(ninter,nfbb) = i1smn;
                    xif2(ninter,nfbb) = i1smx;
                    etf1(ninter,nfbb) = i2smn;
                    etf2(ninter,nfbb) = i2smx;
                    ifrom(ninter,nfbb) = 100*ifrom1 + ifrom2;
                }
                if (ifs == 3 || ifs == 4) {
                    ifrom2 = 21;
                    if (i1smx > kdimg(nbs)) {
                        printf("correcting error in tlns3d map file:\n");
                        printf("i1smx>kdimg(nbs): i1smx,kdimg(nbs) = %d %d\n", i1smx, kdimg(nbs));
                        i1smx = kdimg(nbs);
                    }
                    if (i2smx > idimg(nbs)) {
                        printf("correcting error in tlns3d map file:\n");
                        printf("i2smx>idimg(nbs): i2smx,idimg(nbs) = %d %d\n", i2smx, idimg(nbs));
                        i2smx = idimg(nbs);
                    }
                    if (ifs == 4) ifrom2 = 22;
                    xif1(ninter,nfbb) = i1smn;
                    xif2(ninter,nfbb) = i1smx;
                    etf1(ninter,nfbb) = i2smn;
                    etf2(ninter,nfbb) = i2smx;
                    ifrom(ninter,nfbb) = 100*ifrom1 + ifrom2;
                }
                if (ifs == 5 || ifs == 6) {
                    ifrom2 = 31;
                    if (i2smx > jdimg(nbs)) {
                        printf("correcting error in tlns3d map file:\n");
                        printf("i2smx>jdimg(nbs): i2smx,jdimg(nbs) = %d %d\n", i2smx, jdimg(nbs));
                        i2smx = jdimg(nbs);
                    }
                    if (i1smx > idimg(nbs)) {
                        printf("correcting error in tlns3d map file:\n");
                        printf("i1smx>idimg(nbs): i1smx,idimg(nbs) = %d %d\n", i1smx, idimg(nbs));
                        i1smx = idimg(nbs);
                    }
                    if (ifs == 6) ifrom2 = 32;
                    xif1(ninter,nfbb) = i2smn;
                    xif2(ninter,nfbb) = i2smx;
                    etf1(ninter,nfbb) = i1smn;
                    etf2(ninter,nfbb) = i1smx;
                    ifrom(ninter,nfbb) = 100*ifrom1 + ifrom2;
                }
                nbtlast   = nbt;
                iftlast   = ift;
                i1tmnlast = i1tmn;
                i1tmxlast = i1tmx;
                i2tmnlast = i2tmn;
                i2tmxlast = i2tmx;
            }
        } // end DO 1030 iseg loop

    } // end DO 1000 n loop



    // further concatenate any applicable patches
    idebug = 0;

    if (idebug > 0) {
        fortran_write_unit(21, "%s\n", "first pass through patch interfaces");
        for (int int_var = 1; int_var <= ninter; int_var++) {
            fortran_write_unit(21, "%6d%6d%6d%6d%6d%6d%6d\n",
                int_var, ito(int_var), xit1(int_var), xit2(int_var),
                ett1(int_var), ett2(int_var), nfb(int_var));
            for (int l = 1; l <= nfb(int_var); l++) {
                fortran_write_unit(21, "      %6d%6d%6d%6d%6d\n",
                    ifrom(int_var,l), xif1(int_var,l), xif2(int_var,l),
                    etf1(int_var,l), etf2(int_var,l));
            }
        }
    }

    ninter0 = ninter;

    for (n = 1; n <= ninter0; n++) {
        iconcat(n) = 0;
    }

    ninter = 0;
    for (n = 1; n <= ninter0; n++) {
        if (iconcat(n) == 0) {
            ninter = ninter + 1;
            xit1n = xit1(n);
            xit2n = xit2(n);
            ett1n = ett1(n);
            ett2n = ett2(n);
            iton  = ito(n);
            for (m = 1; m <= ninter0; m++) {
                if (n != m) {
                    xit1m = xit1(m);
                    xit2m = xit2(m);
                    ett1m = ett1(m);
                    ett2m = ett2(m);
                    itom  = ito(m);
                    ichk_var = 0;
                    ichk_var += std::abs((int)(xit1m-xit1n));
                    ichk_var += std::abs((int)(xit2m-xit2n));
                    ichk_var += std::abs((int)(ett1m-ett1n));
                    ichk_var += std::abs((int)(ett2m-ett2n));
                    ichk_var += std::abs(itom-iton);
                    if (ichk_var == 0) {
                        // patches m and n can be concatenated
                        nfbm = nfb(m);
                        for (nf_var = 1; nf_var <= nfbm; nf_var++) {
                            nfb(n) = nfb(n) + 1;
                            nfbb   = nfb(n);
                            xif1(n,nfbb)  = xif1(m,nf_var);
                            xif2(n,nfbb)  = xif2(m,nf_var);
                            etf1(n,nfbb)  = etf1(m,nf_var);
                            etf2(n,nfbb)  = etf2(m,nf_var);
                            ifrom(n,nfbb) = ifrom(m,nf_var);
                        }
                        iconcat(m) = 1;
                        if (idebug > 0) fortran_write_unit(21, " patch %d will concatenate with patch %d\n", n, m);
                    }
                }
            }
        }
    }
    if (idebug > 0) {
        fortran_write_unit(21, "old ninter = %d\n", ninter0);
        fortran_write_unit(21, "new ninter = %d\n", ninter);
    }

    // set default patching parameters iifit,llimit,iitmax,iico,iiorph
    if (ninter > 0) {
        for (n = 1; n <= ninter; n++) {
            iifit(n)  = 1;
            llimit(n) = 1;
            iitmax(n) = 100;
            iic0(n)   = 0;
            iiorph(n) = 0;
            // for now, set mmcxie and mmceta = 50
            mmcxie(n) = 50;
            mmceta(n) = 50;
        }
    }



    // check multigridability of boundary-condition segments
    ncgmax2 = 0;
    iflg = 1;
    for (nnnn = 1; nnnn <= ncgmax+1; nnnn++) {
        nnn = nnnn - 1;
        for (igrid_var = 1; igrid_var <= ngrid; igrid_var++) {
            for (m = 1; m <= 2; m++) {
                if (m == 1) {
                    ns = nbci0(igrid_var);
                    if (ns > maxseg) {
                        printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                        { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                    }
                } else {
                    ns = nbcidim(igrid_var);
                    if (ns > maxseg) {
                        printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                        { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                    }
                }
                for (iseg = 1; iseg <= ns; iseg++) {
                    int pow2nnn = 1; for (int _p=0;_p<nnn;_p++) pow2nnn*=2;
                    jtest1 = (ibcinfo(igrid_var,iseg,2,m)-1)/pow2nnn + 1;
                    jtest2 = (ibcinfo(igrid_var,iseg,3,m)-1)/pow2nnn + 1;
                    ktest1 = (ibcinfo(igrid_var,iseg,4,m)-1)/pow2nnn + 1;
                    ktest2 = (ibcinfo(igrid_var,iseg,5,m)-1)/pow2nnn + 1;
                    if ((jtest1/2*2 == jtest1) || (ktest1/2*2 == ktest1) ||
                        (jtest2/2*2 == jtest2) || (ktest2/2*2 == ktest2)) {
                        iflg = 0;
                    }
                }
            }
            for (m = 1; m <= 2; m++) {
                if (m == 1) {
                    ns = nbcj0(igrid_var);
                    if (ns > maxseg) {
                        printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                        { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                    }
                } else {
                    ns = nbcjdim(igrid_var);
                    if (ns > maxseg) {
                        printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                        { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                    }
                }
                for (iseg = 1; iseg <= ns; iseg++) {
                    int pow2nnn = 1; for (int _p=0;_p<nnn;_p++) pow2nnn*=2;
                    itest1 = (jbcinfo(igrid_var,iseg,2,m)-1)/pow2nnn + 1;
                    if (idimg(igrid_var) == 2) {
                        itest2 = 1;
                    } else {
                        itest2 = (jbcinfo(igrid_var,iseg,3,m)-1)/pow2nnn + 1;
                    }
                    ktest1 = (jbcinfo(igrid_var,iseg,4,m)-1)/pow2nnn + 1;
                    ktest2 = (jbcinfo(igrid_var,iseg,5,m)-1)/pow2nnn + 1;
                    if ((itest1/2*2 == itest1) || (ktest1/2*2 == ktest1) ||
                        (itest2/2*2 == itest2) || (ktest2/2*2 == ktest2)) {
                        iflg = 0;
                    }
                }
            }
            for (m = 1; m <= 2; m++) {
                if (m == 1) {
                    ns = nbck0(igrid_var);
                    if (ns > maxseg) {
                        printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                        { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                    }
                } else {
                    ns = nbckdim(igrid_var);
                    if (ns > maxseg) {
                        printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                        { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                    }
                }
                for (iseg = 1; iseg <= ns; iseg++) {
                    int pow2nnn = 1; for (int _p=0;_p<nnn;_p++) pow2nnn*=2;
                    itest1 = (kbcinfo(igrid_var,iseg,2,m)-1)/pow2nnn + 1;
                    if (idimg(igrid_var) == 2) {
                        itest2 = 1;
                    } else {
                        itest2 = (kbcinfo(igrid_var,iseg,3,m)-1)/pow2nnn + 1;
                    }
                    ktest1 = (kbcinfo(igrid_var,iseg,4,m)-1)/pow2nnn + 1;
                    ktest2 = (kbcinfo(igrid_var,iseg,5,m)-1)/pow2nnn + 1;
                    if ((itest1/2*2 == itest1) || (ktest1/2*2 == ktest1) ||
                        (itest2/2*2 == itest2) || (ktest2/2*2 == ktest2)) {
                        iflg = 0;
                    }
                }
            }
        }
        if (iflg != 0) ncgmax2 = ncgmax2 + 1;
    }

    printf(" \n");
    printf("bc segment dimensions are multigridable to ncg = %2d\n", ncgmax2);
    printf(" \n");

    ncgmax = std::min(ncgmax, ncgmax2);

    printf("from block segment dimensions and overall block dimensions, ncgmax = %2d\n", ncgmax);
    icgflg = 0;
    for (igrd = 1; igrd <= ngrid; igrd++) {
        if (ncgmax < ncgg(igrd)) icgflg = 1;
    }
    if (icgflg > 0) {
        printf("this value will be used for ncg in the split input file\n");
        printf(" \n");
        printf("   NOTE: this is a smaller ncg than used in the unsplit input file\n");
        printf("   you may want to consider an alternative splitting to maintain\n");
        printf("   the same level of multigridability\n");
        for (igrd = 1; igrd <= ngrid; igrd++) {
            ncgg(igrd) = ncgmax;
        }
    }
    icgflg2 = 0;
    for (igrd = 1; igrd <= ngrid; igrd++) {
        if (ncgmax > ncgg(igrd)) icgflg2 = 1;
    }
    if (icgflg2 > 0) {
        printf(" \n");
        printf("   NOTE: this is a larger ncg than used in the unsplit input file\n");
        printf("   the unsplit ncg value will be maintained but you might want to\n");
        printf("   increase the value of ncg and alter the  mulitgrid parameters\n");
        printf("   to reflect this level of multigridability\n");
    }
    printf(" \n");

    ncg = ncgmax;



    // begin generation of cfl3d input file

    // set wall temp data for bc2004 in version 5 based on isnd and c2spe
    if (isnd == 0) {
        twall = 0.0;
    } else {
        if (c2spe > 0.0f) {
            twall = c2spe;
        } else {
            twall = -1.0;
        }
    }

    // write(7,18) format(19hinput/output files:)
    fortran_write_unit(7, "input/output files:\n");

    // grid file name already read in is used - use generic names for other files
    fortran_write_unit(7, "%-60.60s\n", gridout);
    fortran_write_unit(7, "%-60.60s\n", plt3dg);
    fortran_write_unit(7, "%-60.60s\n", plt3dq);
    fortran_write_unit(7, "%-60.60s\n", output);
    fortran_write_unit(7, "%-60.60s\n", resid);
    fortran_write_unit(7, "%-60.60s\n", turbres);
    fortran_write_unit(7, "%-60.60s\n", blomx);
    fortran_write_unit(7, "%-60.60s\n", output2);
    fortran_write_unit(7, "%-60.60s\n", printout);
    fortran_write_unit(7, "%-60.60s\n", pplunge);
    fortran_write_unit(7, "%-60.60s\n", ovrlap);
    fortran_write_unit(7, "%-60.60s\n", patch);
    fortran_write_unit(7, "%-60.60s\n", restrt);

    // check for keyword input
    nskip = 14;
    for (n = 1; n <= nskip; n++) {
        fortran_read_list(10);
    }

    iunit5sav = iunit5;
    iunit5    = 10;
    { int _ititr=0, _myid=0, _iunit11=7, _ierrflg=-1;
      readkey_ns::readkey(_ititr, _myid, ibufdim, nbuf, bou, nou, _iunit11, _ierrflg); }
    iunit5    = iunit5sav;

    // write(7,11)(title(i),i=1,20)  format(1h ,20a4)
    // title is float[20] used as 4-char chunks
    fortran_write_unit(7, " %-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s"
                          "%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s\n",
        (char*)&title[0],(char*)&title[1],(char*)&title[2],(char*)&title[3],
        (char*)&title[4],(char*)&title[5],(char*)&title[6],(char*)&title[7],
        (char*)&title[8],(char*)&title[9],(char*)&title[10],(char*)&title[11],
        (char*)&title[12],(char*)&title[13],(char*)&title[14],(char*)&title[15],
        (char*)&title[16],(char*)&title[17],(char*)&title[18],(char*)&title[19]);

    // handle Cl-specified option
    if (ialphit > 0) {
        if (iver == 4) {
            xmach = -std::abs(xmach);
        } else {
            sref  = -std::abs(sref);
        }
    }
    if (iprecon > 0) {
        if (iver > 4) {
            xmach = -std::abs(xmach);
        }
    }

    if (iver == 4) {
        if (ipar == 0) {
            // format(6x,4hMach,5x,5halpha,6x,4hbeta,6x,4hReUe,3x,7hTinf,dR,6x,4hisnd,5x,5hc2spe)
            fortran_write_unit(7, "      Mach     alpha      beta      ReUe   Tinf,dR      isnd    c2spe\n");
            // format(3f10.5,e10.3,f10.5,i10,f10.5,i10)
            fortran_write_unit(7, "%10.5f%10.5f%10.5f%10.3e%10.5f%10d%10.5f\n",
                (double)xmach, (double)alpha, (double)beta, (double)reue, (double)tinf, isnd, (double)c2spe);
        } else {
            // format(6x,4hMach,5x,5halpha,6x,4hbeta,6x,4hReUe,3x,7hTinf,dR,6x,4hisnd,5x,5hc2spe,5x,5hnodes)
            fortran_write_unit(7, "      Mach     alpha      beta      ReUe   Tinf,dR      isnd    c2spe    nodes\n");
            // format(3f10.5,e10.3,f10.5,i10,f10.5,i10)
            fortran_write_unit(7, "%10.5f%10.5f%10.5f%10.3e%10.5f%10d%10.5f%10d\n",
                (double)xmach, (double)alpha, (double)beta, (double)reue, (double)tinf, isnd, (double)c2spe, nnodes);
        }
    } else {
        // format(6x,4hMach,5x,5halpha,6x,4hbeta,6x,4hReUe,3x,7hTinf,dR,5x,5hialph,4x,6hihstry)
        fortran_write_unit(7, "      Mach     alpha      beta      ReUe   Tinf,dR     ialph    ihstry\n");
        // format(3f10.5,e10.3,f10.5,i10,i10,i10)
        fortran_write_unit(7, "%10.5f%10.5f%10.5f%10.3e%10.5f%10d%10d\n",
            (double)xmach, (double)alpha, (double)beta, (double)reue, (double)tinf, ialph, ihstry);
    }

    xmach = std::abs(xmach);

    if (iver == 4) {
        if (ialphit == 0) {
            // format(6x,4hsref,6x,4hcref,6x,4hbref,7x,3hxmc,7x,3hymc,7x,3hzmc)
            fortran_write_unit(7, "      sref      cref      bref       xmc       ymc       zmc\n");
            // format(f10.3,7f10.4)
            fortran_write_unit(7, "%10.3f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
                (double)sref, (double)cref, (double)bref, (double)xmc, (double)ymc, (double)zmc);
        } else {
            // format(6x,4hsref,6x,4hcref,6x,4hbref,7x,3hxmc,7x,3hymc,7x,3hzmc,4x,6hcltarg,3x,7hresupdt)
            fortran_write_unit(7, "      sref      cref      bref       xmc       ymc       zmc    cltarg   resupdt\n");
            // format(f10.3,7f10.4)
            fortran_write_unit(7, "%10.3f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
                (double)sref, (double)cref, (double)bref, (double)xmc, (double)ymc, (double)zmc,
                (double)cltarg, (double)resupdt);
        }
    } else {
        fortran_write_unit(7, "      sref      cref      bref       xmc       ymc       zmc\n");
        fortran_write_unit(7, "%10.3f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
            (double)sref, (double)cref, (double)bref, (double)xmc, (double)ymc, (double)zmc);
    }

    if (iver == 4) {
        // format(8x,2hdt,5x,5hirest,3x,7hiflagts,6x,4hfmax,5x,5hiunst,5x,5hrfreq,4x,6halphau,6x,4hcloc)
        fortran_write_unit(7, "        dt     irest   iflagts      fmax     iunst     rfreq    alphau      cloc\n");
        // format(f10.5,2i10,f10.5,i10,3f10.5)
        fortran_write_unit(7, "%10.5f%10d%10d%10.5f%10d%10.5f%10.5f%10.5f\n",
            (double)dt, irest, iflagts, (double)fmax, iunst, (double)rfreq, (double)alphau, (double)cloc);
    } else {
        // format(8x,2hdt,5x,5hirest,3x,7hiflagts,6x,4hfmax,5x,5hiunst,3x,7hcfl_tau)
        fortran_write_unit(7, "        dt     irest   iflagts      fmax     iunst   cfl_tau\n");
        // format(f10.5,2i10,f10.5,i10,3f10.5)
        fortran_write_unit(7, "%10.5f%10d%10d%10.5f%10d%10.5f\n",
            (double)dt, irest, iflagts, (double)fmax, iunst, (double)cfltau);
    }



    ngridout = ngrid;
    if (ip3dgrd > 0) ngridout = -ngrid;

    // in version 4.1, assume no plot3d files are output
    if (iver == 4 && ipar == 0) {
        nplot3d = 0;
        nprint  = 0;
    } else {
        if (std::abs(nplot3d) > 0) nplot3d = -1;
        if (std::abs(nprint)  > 0) nprint  = 0;
    }

    nwrest = 9999;
    ichk   = 0;

    // format(5x,5hngrid,3x,7hnplot3d,4x,6hnprint,4x,6hnwrest,6x,4hichk,7x,3hi2d,3x,7h ntstep,4x,6h   ita)
    fortran_write_unit(7, "     ngrid   nplot3d    nprint    nwrest      ichk       i2d  ntstep      ita\n");
    // format(10i10)
    fortran_write_unit(7, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
        ngridout, nplot3d, nprint, nwrest, ichk, i2d, ntstep, ita);

    // format(7x,3hncg,7x,3hiem,2x,8hiadvance,4x,6hiforce,2x,8hivisc(i),2x,8hivisc(j),2x,8hivisc(k))
    fortran_write_unit(7, "       ncg       iem  iadvance    iforce  ivisc(i)  ivisc(j)  ivisc(k)\n");
    for (igrd = 1; igrd <= ngrid; igrd++) {
        // format(10i10)
        fortran_write_unit(7, "%10d%10d%10d%10d%10d%10d%10d\n",
            ncgg(igrd), iemg(igrd), iadvance(igrd), iforce(igrd),
            iviscg(igrd,1), iviscg(igrd,2), iviscg(igrd,3));
    }

    // format(6x,4hidim,6x,4hjdim,6x,4hkdim)
    fortran_write_unit(7, "      idim      jdim      kdim\n");
    i2dd = 0;
    for (igrd = 1; igrd <= ngrid; igrd++) {
        if (i2d > 0 && idimg(igrd) != 2) i2dd = 1;
        // format(3i10)
        fortran_write_unit(7, "%10d%10d%10d\n", idimg(igrd), jdimg(igrd), kdimg(igrd));
    }

    // format(4x,6hilamlo,4x,6hilamhi,4x,6hjlamlo,4x,6hjlamhi,4x,6hklamlo,4x,6hklamhi)
    fortran_write_unit(7, "    ilamlo    ilamhi    jlamlo    jlamhi    klamlo    klamhi\n");
    for (igrd = 1; igrd <= ngrid; igrd++) {
        // format(10i10)
        fortran_write_unit(7, "%10d%10d%10d%10d%10d%10d\n",
            ilamlog(igrd), ilamhig(igrd), jlamlog(igrd),
            jlamhig(igrd), klamlog(igrd), klamhig(igrd));
    }

    // format(5x,5hinewg,4x,6higridc,8x,2his,8x,2hjs,8x,2hks,8x,2hie,8x,2hje,8x,2hke)
    fortran_write_unit(7, "     inewg    igridc        is        js        ks        ie        je        ke\n");
    for (igrd = 1; igrd <= ngrid; igrd++) {
        inewg  = inewgg(igrd);
        igridc = igridg(igrd);
        is_var = isg(igrd);
        js_var = jsg(igrd);
        ks_var = ksg(igrd);
        ie_var = ieg(igrd);
        je_var = jeg(igrd);
        ke_var = keg(igrd);
        // format(10i10)
        fortran_write_unit(7, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
            inewg, igridc, is_var, js_var, ks_var, ie_var, je_var, ke_var);
    }

    // format(2x,8hidiag(i),2x,8hidiag(j),2x,8hidiag(k),2x,8hiflim(i),2x,8hiflim(j),2x,8hiflim(k))
    fortran_write_unit(7, "  idiag(i)  idiag(j)  idiag(k)  iflim(i)  iflim(j)  iflim(k)\n");
    for (igrd = 1; igrd <= ngrid; igrd++) {
        // format(10i10)
        fortran_write_unit(7, "%10d%10d%10d%10d%10d%10d\n",
            idiagg(igrd,1), idiagg(igrd,2), idiagg(igrd,3),
            iflimg(igrd,1), iflimg(igrd,2), iflimg(igrd,3));
    }

    // format(3x,7hifds(i),3x,7hifds(j),3x,7hifds(k),2x,8hrkap0(i),2x,8hrkap0(j),2x,8hrkap0(k))
    fortran_write_unit(7, "   ifds(i)   ifds(j)   ifds(k)  rkap0(i)  rkap0(j)  rkap0(k)\n");
    for (igrd = 1; igrd <= ngrid; igrd++) {
        // format(3i10,3f10.4)
        fortran_write_unit(7, "%10d%10d%10d%10.4f%10.4f%10.4f\n",
            ifdsg(igrd,1), ifdsg(igrd,2), ifdsg(igrd,3),
            rkap0g(igrd,1), rkap0g(igrd,2), rkap0g(igrd,3));
    }

    // format(6x,4hgrid,5x,5hnbci0,3x,7hnbcidim,5x,5hnbcj0,3x,7hnbcjdim,5x,5hnbck0,3x,7hnbckdim,4x,6hiovrlp)
    fortran_write_unit(7, "      grid     nbci0   nbcidim     nbcj0   nbcjdim     nbck0   nbckdim    iovrlp\n");
    for (igrd = 1; igrd <= ngrid; igrd++) {
        // format(10i10)
        fortran_write_unit(7, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
            igrd, nbci0(igrd), nbcidim(igrd), nbcj0(igrd),
            nbcjdim(igrd), nbck0(igrd), nbckdim(igrd), iovrlp(igrd));
    }



    // I0/IDIM Boundary Data
    for (m = 1; m <= 2; m++) {
        if (m == 1) {
            // format(3hi0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hjsta,6x,4hjend,6x,4hksta,6x,4hkend,5x,5hndata)
            fortran_write_unit(7, "i0:   grid   segment    bctype      jsta      jend      ksta      kend     ndata\n");
        } else {
            // format(5hidim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hjsta,6x,4hjend,6x,4hksta,6x,4hkend,5x,5hndata)
            fortran_write_unit(7, "idim: grid   segment    bctype      jsta      jend      ksta      kend     ndata\n");
        }
        for (igrd = 1; igrd <= ngrid; igrd++) {
            if (m == 1) {
                ns = nbci0(igrd);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
            } else {
                ns = nbcidim(igrd);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
            }
            for (iseg = 1; iseg <= ns; iseg++) {
                mseg_loc = iseg;
                if (ibcinfo(igrd,iseg,6,m) == 0) {
                    mseg_loc = -iseg;
                }
                if (iver != 4 && ibcinfo(igrd,iseg,1,m) == 1004) {
                    ibcinfo(igrd,iseg,1,m) = 2004;
                    ibcinfo(igrd,iseg,7,m) = 2;
                    bcvali(igrd,iseg,1,m)  = twall;
                    bcvali(igrd,iseg,2,m)  = 0.0;
                }
                // format(10i10)
                fortran_write_unit(7, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    igrd, mseg_loc, ibcinfo(igrd,iseg,1,m),
                    ibcinfo(igrd,iseg,2,m), ibcinfo(igrd,iseg,3,m),
                    ibcinfo(igrd,iseg,4,m), ibcinfo(igrd,iseg,5,m),
                    ibcinfo(igrd,iseg,7,m));
                if (ibcinfo(igrd,iseg,7,m) > 0) {
                    ndata  = ibcinfo(igrd,iseg,7,m);
                    ibctyp = ibcinfo(igrd,iseg,1,m);
                    {
                        getdhdr_ns::getdhdr(datahdr, ibctyp, iver, ndata);
                    }
                    // format(10a10)
                    { char buf[101]=""; for(int _mm=1;_mm<=ndata;_mm++) { char tmp[11]; snprintf(tmp,11,"%10.10s",datahdr(_mm)); strncat(buf,tmp,10); } fortran_write_unit(7,"%s\n",buf); }
                    // format(10f10.5)
                    { char buf[201]=""; for(int _mm=1;_mm<=ndata;_mm++) { char tmp[11]; snprintf(tmp,11,"%10.5f",(double)bcvali(igrd,iseg,_mm,m)); strncat(buf,tmp,10); } fortran_write_unit(7,"%s\n",buf); }
                }
                if (ibcinfo(igrd,iseg,7,m) < 0) {
                    // format(17h     bc data file)
                    fortran_write_unit(7, "     bc data file\n");
                    ibcfile = (int)bcvali(igrd,iseg,1,m);
                    fortran_write_unit(7, "bc_file%2d\n", ibcfile);
                }
            }
        }
    }

    // J0/JDIM Boundary Data
    for (m = 1; m <= 2; m++) {
        if (m == 1) {
            // format(3hj0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,6x,4hiend,6x,4hksta,6x,4hkend,5x,5hndata)
            fortran_write_unit(7, "j0:   grid   segment    bctype      ista      iend      ksta      kend     ndata\n");
        } else {
            // format(5hjdim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,6x,4hiend,6x,4hksta,6x,4hkend,5x,5hndata)
            fortran_write_unit(7, "jdim: grid   segment    bctype      ista      iend      ksta      kend     ndata\n");
        }
        for (igrd = 1; igrd <= ngrid; igrd++) {
            if (m == 1) {
                ns = nbcj0(igrd);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
            } else {
                ns = nbcjdim(igrd);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
            }
            for (iseg = 1; iseg <= ns; iseg++) {
                mseg_loc = iseg;
                if (jbcinfo(igrd,iseg,6,m) == 0) {
                    mseg_loc = -iseg;
                }
                if (iver != 4 && jbcinfo(igrd,iseg,1,m) == 1004) {
                    jbcinfo(igrd,iseg,1,m) = 2004;
                    jbcinfo(igrd,iseg,7,m) = 2;
                    bcvalj(igrd,iseg,1,m)  = twall;
                    bcvalj(igrd,iseg,2,m)  = 0.0;
                }
                // format(10i10)
                fortran_write_unit(7, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    igrd, mseg_loc, jbcinfo(igrd,iseg,1,m),
                    jbcinfo(igrd,iseg,2,m), jbcinfo(igrd,iseg,3,m),
                    jbcinfo(igrd,iseg,4,m), jbcinfo(igrd,iseg,5,m),
                    jbcinfo(igrd,iseg,7,m));
                if (jbcinfo(igrd,iseg,7,m) > 0) {
                    ndata  = jbcinfo(igrd,iseg,7,m);
                    ibctyp = jbcinfo(igrd,iseg,1,m);
                    {
                        getdhdr_ns::getdhdr(datahdr, ibctyp, iver, ndata);
                    }
                    { char buf[101]=""; for(int _mm=1;_mm<=ndata;_mm++) { char tmp[11]; snprintf(tmp,11,"%10.10s",datahdr(_mm)); strncat(buf,tmp,10); } fortran_write_unit(7,"%s\n",buf); }
                    { char buf[201]=""; for(int _mm=1;_mm<=ndata;_mm++) { char tmp[11]; snprintf(tmp,11,"%10.5f",(double)bcvalj(igrd,iseg,_mm,m)); strncat(buf,tmp,10); } fortran_write_unit(7,"%s\n",buf); }
                }
                if (jbcinfo(igrd,iseg,7,m) < 0) {
                    fortran_write_unit(7, "     bc data file\n");
                    ibcfile = (int)bcvalj(igrd,iseg,1,m);
                    fortran_write_unit(7, "bc_file%2d\n", ibcfile);
                }
            }
        }
    }



    // K0/KDIM Boundary Data
    for (m = 1; m <= 2; m++) {
        if (m == 1) {
            // format(3hk0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,6x,4hiend,6x,4hjsta,6x,4hjend,5x,5hndata)
            fortran_write_unit(7, "k0:   grid   segment    bctype      ista      iend      jsta      jend     ndata\n");
        } else {
            // format(5hkdim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,6x,4hiend,6x,4hjsta,6x,4hjend,5x,5hndata)
            fortran_write_unit(7, "kdim: grid   segment    bctype      ista      iend      jsta      jend     ndata\n");
        }
        for (igrd = 1; igrd <= ngrid; igrd++) {
            if (m == 1) {
                ns = nbck0(igrd);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
            } else {
                ns = nbckdim(igrd);
                if (ns > maxseg) {
                    printf("stopping...must increase parameter maxseg to at least %d\n", ns);
                    { int _z=0, _m3=-3; termn8_ns::termn8(_z,_m3,ibufdim,nbuf,bou,nou); }
                }
            }
            for (iseg = 1; iseg <= ns; iseg++) {
                mseg_loc = iseg;
                if (kbcinfo(igrd,iseg,6,m) == 0) {
                    mseg_loc = -iseg;
                }
                if (iver != 4 && kbcinfo(igrd,iseg,1,m) == 1004) {
                    kbcinfo(igrd,iseg,1,m) = 2004;
                    kbcinfo(igrd,iseg,7,m) = 2;
                    bcvalk(igrd,iseg,1,m)  = twall;
                    bcvalk(igrd,iseg,2,m)  = 0.0;
                }
                // format(10i10)
                fortran_write_unit(7, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    igrd, mseg_loc, kbcinfo(igrd,iseg,1,m),
                    kbcinfo(igrd,iseg,2,m), kbcinfo(igrd,iseg,3,m),
                    kbcinfo(igrd,iseg,4,m), kbcinfo(igrd,iseg,5,m),
                    kbcinfo(igrd,iseg,7,m));
                if (kbcinfo(igrd,iseg,7,m) > 0) {
                    ndata  = kbcinfo(igrd,iseg,7,m);
                    ibctyp = kbcinfo(igrd,iseg,1,m);
                    {
                        getdhdr_ns::getdhdr(datahdr, ibctyp, iver, ndata);
                    }
                    { char buf[101]=""; for(int _mm=1;_mm<=ndata;_mm++) { char tmp[11]; snprintf(tmp,11,"%10.10s",datahdr(_mm)); strncat(buf,tmp,10); } fortran_write_unit(7,"%s\n",buf); }
                    { char buf[201]=""; for(int _mm=1;_mm<=ndata;_mm++) { char tmp[11]; snprintf(tmp,11,"%10.5f",(double)bcvalk(igrd,iseg,_mm,m)); strncat(buf,tmp,10); } fortran_write_unit(7,"%s\n",buf); }
                }
                if (kbcinfo(igrd,iseg,7,m) < 0) {
                    fortran_write_unit(7, "     bc data file\n");
                    ibcfile = (int)bcvalk(igrd,iseg,1,m);
                    fortran_write_unit(7, "bc_file%2d\n", ibcfile);
                }
            }
        }
    }

    // check to see if ncg has been lowered due to the splitting process
    if (icgflg != 0) {
        if (ncgmax == 0) {
            mgflag = 0;
            mseq   = 1;
        }
        if (mseq == 1) {
            if (mgflag != 0) {
                // set mglevg to upper limit consistent with new ncg
                mglev[0] = ncgmax + 1;
                for (int ii = 1; ii <= mglev[0]; ii++) {
                    mitl[(1-1)*5+(ii-1)] = 1;
                }
            } else {
                mglev[0] = 1;
                mitl[(1-1)*5+(1-1)] = 1;
            }
        } else {
            // set mseq to upper limit consistent with new ncg
            if (mseq > ncgmax + 1) mseq = ncgmax + 1;
            if (mgflag != 0) {
                // set mglevg to upper limit consistent with new ncg
                for (int mm_v = 1; mm_v <= mseq; mm_v++) {
                    mglev[mm_v-1] = mm_v;
                    for (int ii = 1; ii <= mglev[mm_v-1]; ii++) {
                        mitl[(mm_v-1)*5+(ii-1)] = 1;
                    }
                }
            } else {
                for (int mm_v = 1; mm_v <= mseq; mm_v++) {
                    mglev[mm_v-1] = 1;
                    mitl[(mm_v-1)*5+(1-1)] = 1;
                }
            }
        }
    }

    if (iver == 4) {
        // format(6x,4hmseq,4x,6hmgflag,4x,6hiconsf,7x,3hmtt,6x,4hngam)
        fortran_write_unit(7, "      mseq    mgflag    iconsf       mtt      ngam\n");
        // format(10i10)
        fortran_write_unit(7, "%10d%10d%10d%10d%10d\n", mseq, mgflag, iconsf, mtt, ngam);
    } else {
        if (iprecon == 0) {
            fortran_write_unit(7, "      mseq    mgflag    iconsf       mtt      ngam\n");
            fortran_write_unit(7, "%10d%10d%10d%10d%10d\n", mseq, mgflag, iconsf, mtt, ngam);
        } else {
            // format(6x,4hmseq,4x,6hmgflag,4x,6hiconsf,7x,3hmtt,6x,4hngam,5x,5hcprec,6x,4huref,7x,3havn)
            fortran_write_unit(7, "      mseq    mgflag    iconsf       mtt      ngam     cprec      uref       avn\n");
            // format(5i10,3f10.4)
            fortran_write_unit(7, "%10d%10d%10d%10d%10d%10.4f%10.4f%10.4f\n",
                mseq, mgflag, iconsf, mtt, ngam, (double)cprec, (double)uref, (double)avn);
        }
    }

    // format(6x,4hissc,1x,9hepsssc(1),1x,9hepsssc(2),1x,9hepsssc(3),6x,4hissr,1x,9hepsssr(1),1x,9hepsssr(2),1x,9hepsssr(3))
    fortran_write_unit(7, "      issc epsssc(1) epsssc(2) epsssc(3)      issr epsssr(1) epsssr(2) epsssr(3)\n");
    // format(2(i10,3f10.4))
    fortran_write_unit(7, "%10d%10.4f%10.4f%10.4f%10d%10.4f%10.4f%10.4f\n",
        issc, (double)epsssc[0], (double)epsssc[1], (double)epsssc[2],
        issr, (double)epsssr[0], (double)epsssr[1], (double)epsssr[2]);

    // format(6x,4hncyc,4x,6hmglevg,5x,5hnemgl,5x,5hnitfo)
    fortran_write_unit(7, "      ncyc    mglevg     nemgl     nitfo\n");
    for (m = 1; m <= mseq; m++) {
        // format(10i10)
        fortran_write_unit(7, "%10d%10d%10d%10d\n",
            ncyc1[m-1], mglev[m-1], nem[m-1], nitfo1[m-1]);
    }

    // format(6x,4hmit1,6x,4hmit2,6x,4hmit3,6x,4hmit4,6x,9hmit5  ...)
    fortran_write_unit(7, "      mit1      mit2      mit3      mit4      mit5  ...\n");
    for (m = 1; m <= mseq; m++) {
        // write(7,36)(mitl(i,m),i=1,(mglev(m)+nem(m)))
        int nlev = mglev[m-1] + nem[m-1];
        char buf[201] = "";
        for (int ii = 1; ii <= nlev; ii++) {
            char tmp[11];
            snprintf(tmp, 11, "%10d", mitl[(m-1)*5+(ii-1)]);
            strncat(buf, tmp, 10);
        }
        fortran_write_unit(7, "%s\n", buf);
    }



    // 1-1 blocking data - version 4.x
    if (iver <= 4) {
        // format(19h 1-1 blocking data:)
        fortran_write_unit(7, " 1-1 blocking data:\n");
        // format(6x,4hnbli)
        fortran_write_unit(7, "      nbli\n");
        // format(10i10)
        fortran_write_unit(7, "%10d\n", nbli);
        // format(5x,5hnblon,3x,7hnblk(1),3x,7hnblk(2))
        fortran_write_unit(7, "     nblon   nblk(1)   nblk(2)\n");
        if (nbli > 0) {
            for (n = 1; n <= nbli; n++) {
                // format(10i10)
                fortran_write_unit(7, "%10d%10d%10d\n", nblon(n), nblk(1,n), nblk(2,n));
            }
        }
        // format(6x,4hista,6x,4hjsta,6x,4hksta,6x,4hiend,6x,4hjend,6x,4hkend,1x,9hisva(1,1),1x,9hisva(1,2))
        fortran_write_unit(7, "      ista      jsta      ksta      iend      jend      kend isva(1,1) isva(1,2)\n");
        if (nbli > 0) {
            for (n = 1; n <= nbli; n++) {
                // format(10i10) - (limblk(1,l,n),l=1,6),(isva(1,ind,n),ind=1,2)
                fortran_write_unit(7, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    limblk(1,1,n), limblk(1,2,n), limblk(1,3,n),
                    limblk(1,4,n), limblk(1,5,n), limblk(1,6,n),
                    isva(1,1,n), isva(1,2,n));
            }
        }
        // format(6x,4hista,6x,4hjsta,6x,4hksta,6x,4hiend,6x,4hjend,6x,4hkend,1x,9hisva(2,1),1x,9hisva(2,2))
        fortran_write_unit(7, "      ista      jsta      ksta      iend      jend      kend isva(2,1) isva(2,2)\n");
        if (nbli > 0) {
            for (n = 1; n <= nbli; n++) {
                // format(10i10) - (limblk(2,l,n),l=1,6),(isva(2,ind,n),ind=1,2)
                fortran_write_unit(7, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    limblk(2,1,n), limblk(2,2,n), limblk(2,3,n),
                    limblk(2,4,n), limblk(2,5,n), limblk(2,6,n),
                    isva(2,1,n), isva(2,2,n));
            }
        }
    } else {
        // 1-1 blocking data - version 5.x
        fortran_write_unit(7, " 1-1 blocking data:\n");
        fortran_write_unit(7, "      nbli\n");
        // format(10i8)
        fortran_write_unit(7, "%8d\n", nbli);
        // format(2x,6hnumber,4x,4hgrid,4x,4hista,4x,4hjsta,4x,4hksta,4x,4hiend,4x,4hjend,4x,4hkend,3x,5hisva1,3x,5hisva2)
        fortran_write_unit(7, "  number    grid    ista    jsta    ksta    iend    jend    kend   isva1   isva2\n");
        if (nbli > 0) {
            int mm_v = 0;
            for (n = 1; n <= nbli; n++) {
                mm_v = mm_v + 1;
                // format(10i8)
                fortran_write_unit(7, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
                    mm_v, nblk(1,n),
                    limblk(1,1,n), limblk(1,2,n), limblk(1,3,n),
                    limblk(1,4,n), limblk(1,5,n), limblk(1,6,n),
                    isva(1,1,n), isva(1,2,n));
            }
        }
        fortran_write_unit(7, "  number    grid    ista    jsta    ksta    iend    jend    kend   isva1   isva2\n");
        if (nbli > 0) {
            int mm_v = 0;
            for (n = 1; n <= nbli; n++) {
                mm_v = mm_v + 1;
                // format(10i8)
                fortran_write_unit(7, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
                    mm_v, nblk(2,n),
                    limblk(2,1,n), limblk(2,2,n), limblk(2,3,n),
                    limblk(2,4,n), limblk(2,5,n), limblk(2,6,n),
                    isva(2,1,n), isva(2,2,n));
            }
        }
    }

    // patched grid flag
    ninter1 = 0;
    if (ninter > 0) ninter1 = -1;
    // format(22h patch interface data:)
    fortran_write_unit(7, " patch interface data:\n");
    // format(4x,6hninter)
    fortran_write_unit(7, "    ninter\n");
    // format(10i10)
    fortran_write_unit(7, "%10d\n", ninter1);

    // plot3d output data
    fortran_write_unit(7, " plot3d output:\n");
    // format(2x,4hgrid,1x,5hiptyp,2x,4hista,2x,4hiend,2x,4hiinc,2x,4hjsta,2x,4hjend,2x,4hjinc,2x,4hksta,2x,4hkend,2x,4hkinc)
    fortran_write_unit(7, "  grid iptyp  ista  iend  iinc  jsta  jend  jinc  ksta  kend  kinc\n");
    if (iver == 4 && ipar > 0 && std::abs(nplot3d) > 0) {
        // format(11i6)
        fortran_write_unit(7, "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
            1, iplt3dtyp, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    }
    if (iver >= 5) {
        // format(11i6)
        fortran_write_unit(7, "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
            1, iplt3dtyp, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    // no movie output - just output a header line in the input file
    movie = 0;
    // format(1x,5hmovie)
    fortran_write_unit(7, " movie\n");
    if ((iver == 4 && nplot3d != 0) || iver >= 5) {
        // format(i6)
        fortran_write_unit(7, "%6d\n", movie);
    }

    // print out data
    fortran_write_unit(7, " print out:\n");
    fortran_write_unit(7, "  grid iptyp  ista  iend  iinc  jsta  jend  jinc  ksta  kend  kinc\n");
    if (iver == 4 && ipar > 0 && std::abs(nprint) > 0) {
        // format(11i6)
        fortran_write_unit(7, "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    // no control surface data - just output a header line in the input file
    ncs = 0;
    if (iver > 4) {
        fortran_write_unit(7, " control surfaces:\n");
        fortran_write_unit(7, "   ncs\n");
        // format(i6)
        fortran_write_unit(7, "%6d\n", ncs);
        // format(2x,4hgrid,2x,4hista,2x,4hiend,2x,4hjsta,2x,4hjend,2x,4hksta,2x,4hkend,1x,5hiwall,1x,5hinorm)
        fortran_write_unit(7, "  grid  ista  iend  jsta  jend  ksta  kend iwall inorm\n");
    }



    // dynamic mesh data
    // count number of translating/rotating zones in the split grid
    ntrans = 0;
    nrotat = 0;
    for (igrd = 1; igrd <= ngrid; igrd++) {
        itrans = (int)xmap( 8,1,igrd);
        irotat = (int)xmap(16,1,igrd);
        if (itrans > 0) {
            ntrans = ntrans + 1;
        }
        if (irotat > 0) {
            nrotat = nrotat + 1;
        }
    }

    fortran_write_unit(7, " moving grid data - translation\n");
    fortran_write_unit(7, " ntrans\n");
    // format(i7)
    fortran_write_unit(7, "%7d\n", ntrans);
    fortran_write_unit(7, "   lref\n");
    if (ntrans > 0) {
        // format(f9.3)
        fortran_write_unit(7, "%9.3f\n", tlref);
    }
    fortran_write_unit(7, "   grid   itrans    rfreq     xmag     ymag     zmag\n");
    if (ntrans > 0) {
        for (igrd = 1; igrd <= ngrid; igrd++) {
            itrans = (int)xmap( 8,1,igrd);
            rfreqt = xmap( 9,1,igrd);
            xmag   = xmap(10,1,igrd);
            ymag   = xmap(11,1,igrd);
            zmag   = xmap(12,1,igrd);
            if (itrans > 0) {
                // format(i7,i9,4f9.3)
                fortran_write_unit(7, "%7d%9d%9.3f%9.3f%9.3f%9.3f\n",
                    igrd, itrans, rfreqt, xmag, ymag, zmag);
            }
        }
    }
    fortran_write_unit(7, "   grid    dxmax    dymax    dzmax\n");
    if (ntrans > 0) {
        for (igrd = 1; igrd <= ngrid; igrd++) {
            itrans = (int)xmap( 8,1,igrd);
            dxmax  = xmap(13,1,igrd);
            dymax  = xmap(14,1,igrd);
            dzmax  = xmap(15,1,igrd);
            if (itrans > 0) {
                // format(i7,4f9.3)
                fortran_write_unit(7, "%7d%9.3f%9.3f%9.3f\n", igrd, dxmax, dymax, dzmax);
            }
        }
    }

    fortran_write_unit(7, " moving grid data - rotation\n");
    fortran_write_unit(7, " nrotat\n");
    // format(i7)
    fortran_write_unit(7, "%7d\n", nrotat);
    fortran_write_unit(7, "      lref\n");
    if (nrotat > 0) {
        // format(f10.4)
        fortran_write_unit(7, "%10.4f\n", rlref);
    }
    fortran_write_unit(7, "   grid   irotat    rfreq   thxmag   thymag   thzmag   xorig   yorig   zorig\n");
    if (nrotat > 0) {
        for (igrd = 1; igrd <= ngrid; igrd++) {
            irotat = (int)xmap(16,1,igrd);
            rfreqr = xmap(17,1,igrd);
            thxmag = xmap(18,1,igrd);
            thymag = xmap(19,1,igrd);
            thzmag = xmap(20,1,igrd);
            xorig  = xmap(22,1,igrd);
            yorig  = xmap(23,1,igrd);
            zorig  = xmap(23,1,igrd);
            if (irotat > 0) {
                // format(i7,i9,7f9.3)
                fortran_write_unit(7, "%7d%9d%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f\n",
                    igrd, irotat, rfreqr, thxmag, thymag, thzmag, xorig, yorig, zorig);
            }
        }
    }
    fortran_write_unit(7, "   grid   thxmax   thymax   thzmax\n");
    if (nrotat > 0) {
        for (igrd = 1; igrd <= ngrid; igrd++) {
            irotat = (int)xmap(16,1,igrd);
            thxmax = xmap(24,1,igrd);
            thymax = xmap(25,1,igrd);
            thzmax = xmap(26,1,igrd);
            if (irotat > 0) {
                // format(i7,3f9.3)
                fortran_write_unit(7, "%7d%9.3f%9.3f%9.3f\n", igrd, thxmax, thymax, thzmax);
            }
        }
    }

    // ignore dynamic mesh data for now!

    // sensitivity data (version 4 only)
    if (iver == 4 && isd > 0) {
        fortran_write_unit(7, " sensitivity data:\n");
        fortran_write_unit(7, "        ndv   isdform\n");
        // format(2i10)
        fortran_write_unit(7, "%10d%10d\n", ndv, isdform);
        fortran_write_unit(7, " sensitivity i/o files\n");
        fortran_write_unit(7, "%-60.60s\n", sdgridout);
        fortran_write_unit(7, "%-60.60s\n", dovrlap);
        fortran_write_unit(7, "%-60.60s\n", dpatch);
        fortran_write_unit(7, "%-60.60s\n", dresid);
    }



    // begin generation of ronnie patched-grid input file
    if (ninter > 0) {
        // write(9,18) format(19hinput/output files:)
        fortran_write_unit(9, "input/output files:\n");
        fortran_write_unit(9, "%-60.60s\n", gridout);
        fortran_write_unit(9, "%-60.60s\n", rout);
        fortran_write_unit(9, "%-60.60s\n", patch);

        if (ironver == 1) {
            int itrace = -1;
            fortran_write_unit(9, "ioflag  itrace\n");
            // format(i6,2x,i6)
            fortran_write_unit(9, "%6d  %6d\n", ioflag, itrace);
        } else {
            ioflag = 0;
        }
        // write(9,11)(titleron(i),i=1,20)  format(1h ,20a4)
        fortran_write_unit(9, " %-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s"
                              "%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s%-4.4s\n",
            (char*)&titleron[0],(char*)&titleron[1],(char*)&titleron[2],(char*)&titleron[3],
            (char*)&titleron[4],(char*)&titleron[5],(char*)&titleron[6],(char*)&titleron[7],
            (char*)&titleron[8],(char*)&titleron[9],(char*)&titleron[10],(char*)&titleron[11],
            (char*)&titleron[12],(char*)&titleron[13],(char*)&titleron[14],(char*)&titleron[15],
            (char*)&titleron[16],(char*)&titleron[17],(char*)&titleron[18],(char*)&titleron[19]);

        // format(1x,5hngrid)
        fortran_write_unit(9, " ngrid\n");
        // format(13i6)
        fortran_write_unit(9, "%6d\n", ngridout);

        // format(3x,3hncg,3x,3hiem,2x,4hidim,2x,4hjdim,2x,4hkdim)
        fortran_write_unit(9, "   ncg   iem  idim  jdim  kdim\n");
        for (n = 1; n <= ngrid; n++) {
            // format(13i6)
            fortran_write_unit(9, "%6d%6d%6d%6d%6d\n",
                ncg, iemg(n), idimg(n), jdimg(n), kdimg(n));
        }

        // format(1x,6hninter)
        fortran_write_unit(9, " ninter\n");
        // format(13i6)
        fortran_write_unit(9, "%6d\n", ninter);
        // format(3x,3hint,1x,5hiifit,1x,5hlimit,1x,5hitmax,1x,5hmcxie,1x,5hmceta,3x,3hc-0,1x,5hiorph)
        fortran_write_unit(9, "   int iifit limit itmax mcxie mceta  c-0 iorph\n");
        for (int int_var = 1; int_var <= ninter; int_var++) {
            // format(13i6)
            fortran_write_unit(9, "%6d%6d%6d%6d%6d%6d%6d%6d\n",
                int_var, iifit(int_var), llimit(int_var), iitmax(int_var),
                mmcxie(int_var), mmceta(int_var), iic0(int_var), iiorph(int_var));
        }

        if (ioflag == 0) {
            // format(3x,3hint,4x,2hto,2x,4hxie1,2x,4hxie2,2x,4heta1,2x,4heta2,3x,3hnfb,4(2x,4hfrom),3x,3h...)
            fortran_write_unit(9, "   int    to  xie1  xie2  eta1  eta2  nfb    from    from    from    from ...\n");
            in_var = 0;
            for (int int_var = 1; int_var <= ninter0; int_var++) {
                if (iconcat(int_var) == 0) {
                    in_var = in_var + 1;
                    // format(13i6): in,ito(int),xit1(int),xit2(int),ett1(int),ett2(int),nfb(int),(ifrom(int,l),l=1,nfb(int))
                    char buf[500] = "";
                    char tmp[7];
                    snprintf(tmp,7,"%6d",in_var);          strncat(buf,tmp,6);
                    snprintf(tmp,7,"%6d",ito(int_var));    strncat(buf,tmp,6);
                    snprintf(tmp,7,"%6d",xit1(int_var));   strncat(buf,tmp,6);
                    snprintf(tmp,7,"%6d",xit2(int_var));   strncat(buf,tmp,6);
                    snprintf(tmp,7,"%6d",ett1(int_var));   strncat(buf,tmp,6);
                    snprintf(tmp,7,"%6d",ett2(int_var));   strncat(buf,tmp,6);
                    snprintf(tmp,7,"%6d",nfb(int_var));    strncat(buf,tmp,6);
                    for (int l = 1; l <= nfb(int_var); l++) {
                        snprintf(tmp,7,"%6d",ifrom(int_var,l)); strncat(buf,tmp,6);
                    }
                    fortran_write_unit(9, "%s\n", buf);
                }
            }
        } else {
            // format(3x,3hint,4x,2hto,2x,4hxie1,2x,4hxie2,2x,4heta1,2x,4heta2,3x,3hnfb,/,8x,4hfrom,2x,4hxie1,2x,4hxie2,2x,4heta1,2x,4heta2)
            fortran_write_unit(9, "   int    to  xie1  xie2  eta1  eta2  nfb\n");
            fortran_write_unit(9, "        from  xie1  xie2  eta1  eta2\n");
            in_var = 0;
            for (int int_var = 1; int_var <= ninter0; int_var++) {
                if (iconcat(int_var) == 0) {
                    in_var = in_var + 1;
                    // format(13i6)
                    fortran_write_unit(9, "%6d%6d%6d%6d%6d%6d%6d\n",
                        in_var, ito(int_var), xit1(int_var), xit2(int_var),
                        ett1(int_var), ett2(int_var), nfb(int_var));
                    for (int l = 1; l <= nfb(int_var); l++) {
                        // format(6x,13i6)
                        fortran_write_unit(9, "      %6d%6d%6d%6d%6d\n",
                            ifrom(int_var,l), xif1(int_var,l), xif2(int_var,l),
                            etf1(int_var,l), etf2(int_var,l));
                    }
                }
            }
        }
    } // end if (ninter > 0)

    return;
} // end tlns3d_to_cfl3d

} // namespace tlns3d_to_cfl3d_ns
