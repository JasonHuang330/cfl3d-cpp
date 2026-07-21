// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "global2.h"
#include "global.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace global2_ns {

// global() in global2_ns forwards to global_ns::global
void global(int& myid, int& maxbl, int& maxgr, int& maxseg, int& maxcs, int& nplots, int& mxbli, FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, FortranArray2DRef<int> nblk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, FortranArray2DRef<double> rkap0g, int& nblock, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<int> idiagg, FortranArray1DRef<int> nblcg, FortranArray2DRef<int> idegg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray2DRef<int> mit, FortranArray1DRef<int> ilamlog, FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamlog, FortranArray1DRef<int> klamhig, FortranArray2DRef<int> iwfg, FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx, FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz, FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr, FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg, FortranArray1DRef<int> iovrlp, int& ninter, int& nplot3d, FortranArray2DRef<int> inpl3d, int& ip3dsurf, int& nprint, FortranArray2DRef<int> inpr, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> iforce, int& lfgm, int& ncs, FortranArray2DRef<int> icsinfo, int& ihstry, int& ncycmax, FortranArray1DRef<int> iv, FortranArray1DRef<double> time2, FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl, int& intmax, int& nsub1, FortranArray2DRef<int> iindex, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, int& icall, int& iunit11, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& ibufdim, int& nbuf, FortranArray1DRef<int> mglevg, FortranArray1DRef<int> nemgl, FortranArray2DRef<int> ipl3dtmp, int& ntr, FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae, FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae, FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass, FortranArray2DRef<double> damp, FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0, int& nmds, int& maxaes, FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg, FortranArray2DRef<double> xorgae0, FortranArray2DRef<double> yorgae0, FortranArray2DRef<double> zorgae0, FortranArray2DRef<int> icouple, int& iprnsurf)
{
    global_ns::global(myid, maxbl, maxgr, maxseg, maxcs, nplots, mxbli, bcvali, bcvalj, bcvalk, nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim, ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek, nblk, nbli, limblk, isva, nblon, rkap0g, nblock, levelg, igridg, iflimg, ifdsg, iviscg, jdimg, kdimg, idimg, idiagg, nblcg, idegg, jsg, ksg, isg, jeg, keg, ieg, mit, ilamlog, ilamhig, jlamlog, jlamhig, klamlog, klamhig, iwfg, utrans, vtrans, wtrans, omegax, omegay, omegaz, xorig, yorig, zorig, dxmx, dymx, dzmx, dthxmx, dthymx, dthzmx, thetax, thetay, thetaz, rfreqt, rfreqr, xorig0, yorig0, zorig0, itrans, irotat, idefrm, ngrid, ncgg, nblg, iemg, inewgg, iovrlp, ninter, nplot3d, inpl3d, ip3dsurf, nprint, inpr, iadvance, iforce, lfgm, ncs, icsinfo, ihstry, ncycmax, iv, time2, thetaxl, thetayl, thetazl, intmax, nsub1, iindex, lig, lbg, ibpntsg, iipntsg, icall, iunit11, nou, bou, ibufdim, nbuf, mglevg, nemgl, ipl3dtmp, ntr, bcfiles, mxbcfil, utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae, xorgae, yorgae, zorgae, thtxae, thtyae, thtzae, rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf, kcsi, kcsf, freq, gmass, damp, x0, gf0, nmds, maxaes, aesrfdat, perturb, iskip, jskip, kskip, nsegdfrm, idfrmseg, iaesurf, maxsegdg, xorgae0, yorgae0, zorgae0, icouple, iprnsurf);
}

void global2(int& maxbl, int& maxgr, int& msub1, int& nintr, int& intmx, int& ngrid,
             FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg,
             FortranArray1DRef<int> levelg, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg,
             FortranArray2DRef<int> iindx, FortranArray1DRef<int> llimit, FortranArray1DRef<int> iitmax,
             FortranArray1DRef<int> mmcxie, FortranArray1DRef<int> mmceta, FortranArray1DRef<int> ncheck,
             FortranArray1DRef<int> iifit, FortranArray1DRef<int> iic0, FortranArray1DRef<int> iiorph,
             FortranArray1DRef<int> iitoss, FortranArray1DRef<int> ifiner,
             FortranArray2DRef<double> dx, FortranArray2DRef<double> dy, FortranArray2DRef<double> dz,
             FortranArray2DRef<double> dthetx, FortranArray2DRef<double> dthety, FortranArray2DRef<double> dthetz,
             int& myid, int& mptch, int& mxxe, int& icall, int& iunit25,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& ibufdim, int& nbuf,
             FortranArray1DRef<int> ifrom, FortranArray1DRef<int> xif1, FortranArray1DRef<int> etf1,
             FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf2,
             FortranArray1DRef<int> igridg, FortranArray1DRef<int> iemg,
             int& nblock, int& ioflag, int& imode)
{
    // COMMON block aliases
    int32_t& lmaxgr   = cmn_params.lmaxgr;
    int32_t& lmaxbl   = cmn_params.lmaxbl;
    int32_t& lmxseg   = cmn_params.lmxseg;
    int32_t& lmaxcs   = cmn_params.lmaxcs;
    int32_t& lnplts   = cmn_params.lnplts;
    int32_t& lmxbli   = cmn_params.lmxbli;
    int32_t& lmaxxe   = cmn_params.lmaxxe;
    int32_t& lnsub1   = cmn_params.lnsub1;
    int32_t& lintmx   = cmn_params.lintmx;
    int32_t& lmxxe    = cmn_params.lmxxe;
    int32_t& liitot   = cmn_params.liitot;
    int32_t& isum     = cmn_params.isum;
    int32_t& lncycm   = cmn_params.lncycm;
    int32_t& isum_n   = cmn_params.isum_n;
    int32_t& lminnode = cmn_params.lminnode;
    int32_t& isumi    = cmn_params.isumi;
    int32_t& isumi_n  = cmn_params.isumi_n;
    int32_t& lmptch   = cmn_params.lmptch;
    int32_t& lmsub1   = cmn_params.lmsub1;
    int32_t& lintmax  = cmn_params.lintmax;
    int32_t& libufdim = cmn_params.libufdim;
    int32_t& lnbuf    = cmn_params.lnbuf;
    int32_t& llbcprd  = cmn_params.llbcprd;
    int32_t& llbcemb  = cmn_params.llbcemb;
    int32_t& llbcrad  = cmn_params.llbcrad;
    int32_t& lnmds    = cmn_params.lnmds;
    int32_t& lmaxaes  = cmn_params.lmaxaes;
    int32_t& lnslave  = cmn_params.lnslave;
    int32_t& lmxsegdg = cmn_params.lmxsegdg;
    int32_t& lnmaster = cmn_params.lnmaster;

    // /info/ aliases
    float* title   = cmn_info.title;   // 1-based: title[0..19]
    float* rkap    = cmn_info.rkap;
    float& xmach   = cmn_info.xmach;
    float& alpha   = cmn_info.alpha;
    float& beta    = cmn_info.beta;
    float& dt      = cmn_info.dt;
    float& fmax    = cmn_info.fmax;
    int32_t& nit   = cmn_info.nit;
    int32_t& ntt   = cmn_info.ntt;
    int32_t* idiag = cmn_info.idiag;
    int32_t& nitfo = cmn_info.nitfo;
    int32_t& iflagts = cmn_info.iflagts;
    int32_t* iflim = cmn_info.iflim;
    int32_t& nres  = cmn_info.nres;
    int32_t* levelb = cmn_info.levelb;
    int32_t& mgflag = cmn_info.mgflag;
    int32_t& iconsf = cmn_info.iconsf;
    int32_t& mseq  = cmn_info.mseq;
    int32_t* ncyc1 = cmn_info.ncyc1;
    int32_t* levelt = cmn_info.levelt;
    int32_t* nitfo1 = cmn_info.nitfo1;
    int32_t& ngam  = cmn_info.ngam;
    int32_t* nsm   = cmn_info.nsm;
    int32_t& iipv  = cmn_info.iipv;

    // /conversion/
    float& radtodeg = cmn_conversion.radtodeg;

    // /unit5/
    int32_t& iunit5 = cmn_unit5.iunit5;

    // /filenam/ - grdmov
    char* grdmov = cmn_filenam.grdmov;

    // /igrdtyp/
    int32_t& ip3dgrd = cmn_igrdtyp.ip3dgrd;
    int32_t& ialph   = cmn_igrdtyp.ialph;

    // Local variables
    int ierrflg;
    int ncgmax = 0, ntest, nfbmax;
    int ntemp, ifirst;
    int igrid, ncg, nbl, mbl;
    int n, l, m, i;
    int int_val, ito, itonew, nfb, nfb1;
    int xi1, xi2, et1, et2;
    int xif11, xif22, etf11, etf22;
    int id1, id2, jd1, jd2, kd1, kd2;
    int itop, mpa, mxxchk;
    int lprev, nbl_local;
    int lfem, lfgm_local, lcgm, igptot, igpts;
    int lf, levelc;
    int istop, ichk;
    int nfbl, npts, nst, ntot;
    float realval[21]; // 1-based: use indices 1..20

    FILE* fu25 = fortran_get_unit(iunit25);
    FILE* fu11 = fortran_get_unit(11);
    FILE* fu8  = fortran_get_unit(8);



    // Lines 108-122: init ierrflg, read nintr
    if (icall == 0) {
        ierrflg = -99;
    } else {
        ierrflg = -1;
    }

    nintr = 0;

    if (imode == 1) {
        fortran_read_list(iunit5);
    }
    fortran_read_list(iunit5);
    fortran_read_list(iunit5, &nintr);

    if (nintr == 0) return;

    // Lines 124-153
    lnbuf = lnbuf + 2;

    if (ioflag == 2) {
        fprintf(fu25, "\n %s\n", "dynamic patch input data");
        fprintf(fu25, " %s\n", "nintr");
    } else {
        fprintf(fu25, "\n%s\n", "patch input data");
        fprintf(fu25, " %s\n", "ninter");
    }

    // do 680 igrid=1,ngrid
    for (igrid = 1; igrid <= ngrid; igrid++) {
        ncg = ncgg(igrid);
        if (igrid == 1) {
            ncgmax = ncg;
        } else {
            ncgmax = std::max(ncgmax, ncg);
        }
    }
    ntest = (ncgmax + 1) * nintr;
    if (ntest > intmx) {
        fprintf(fu25, " too many blocks to interpolate; intmx should be at least %4d\n", ntest);
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    // write(iunit25,36)nintr  -- format 36: 13i6
    fprintf(fu25, "%6d\n", nintr);

    // Lines 213-262: read iifit/llimit/iitmax etc.
    fortran_read_list(iunit5);
    if (ioflag == 2) {
        fprintf(fu25, "   int iifit    limit    itmax    mcxie    mceta       c-0    iorph    itoss\n");
    } else {
        fprintf(fu25, "   int iifit    limit    itmax    mcxie    mceta       c-0    iorph\n");
    }

    if (nintr > 0) {
        // do 7030 n=1,nintr: iifit(n) = -99
        for (n = 1; n <= nintr; n++) {
            iifit(n) = -99;
        }
        // do 7031 n=1,nintr
        for (n = 1; n <= nintr; n++) {
            int iif;
            // Fortran: read(iunit5,*) int_val,iif,llimit(int_val),iitmax(int_val),...
            // int_val is read FIRST and then indexes the arrays in the SAME record.
            // A single variadic call would compute &arr(int_val) with the STALE
            // int_val; read into temps, then store using the freshly-read int_val.
            int t_llim=0, t_itmax=0, t_mcxie=0, t_mceta=0, t_ic0=0, t_iorph=0, t_itoss=0;
            if (ioflag == 2) {
                fortran_read_list(iunit5, &int_val, &iif, &t_llim, &t_itmax,
                                  &t_mcxie, &t_mceta, &t_ic0, &t_iorph, &t_itoss);
            } else {
                fortran_read_list(iunit5, &int_val, &iif, &t_llim, &t_itmax,
                                  &t_mcxie, &t_mceta, &t_ic0, &t_iorph);
                t_itoss = 0;
            }
            llimit(int_val) = t_llim;   iitmax(int_val) = t_itmax;
            mmcxie(int_val) = t_mcxie;  mmceta(int_val) = t_mceta;
            iic0(int_val)   = t_ic0;    iiorph(int_val) = t_iorph;
            iitoss(int_val) = t_itoss;
            if (iifit(int_val) != -99) {
                // write(11,...) grdmov
                char grdmov_buf[61];
                std::memcpy(grdmov_buf, grdmov, 60);
                grdmov_buf[60] = '\0';
                fprintf(fu11, " program terminated in dynamic patching routines - see file %60s\n", grdmov_buf);
                fprintf(fu25, " stopping...attempting to set data for int %d more than once\n", int_val);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            iifit(int_val) = iif;
            if (iifit(int_val) == 0) iifit(int_val) = 1;
            if (iic0(int_val) > 0) {
                if (iifit(int_val) >= 0) {
                    iifit(int_val) = 1;
                } else {
                    iifit(int_val) = -1;
                }
                mmcxie(int_val) = 0;
                mmceta(int_val) = 0;
            }
            if (ioflag == 2) {
                // format 38: i6,i6,10i9
                fprintf(fu25, "%6d%6d%9d%9d%9d%9d%9d%9d%9d\n",
                        int_val, iifit(int_val), llimit(int_val), iitmax(int_val),
                        mmcxie(int_val), mmceta(int_val), iic0(int_val),
                        iiorph(int_val), iitoss(int_val));
            } else {
                fprintf(fu25, "%6d%6d%9d%9d%9d%9d%9d%9d\n",
                        int_val, iifit(int_val), llimit(int_val), iitmax(int_val),
                        mmcxie(int_val), mmceta(int_val), iic0(int_val),
                        iiorph(int_val));
            }
        }
    }



    // Lines 264-683: read "to"/"from" block data, set iindx
    nfbmax = 0;
    // do 1702 n=1,nintr: iindx(n,1) = -99
    for (n = 1; n <= nintr; n++) {
        iindx(n, 1) = -99;
    }

    // format strings (defined as labels in Fortran):
    // 8146: int,to,xie1,xie2,eta1,eta2,nfb,4*(from),...
    // 8147: int,to,xie1,xie2,eta1,eta2,nfb
    // 8148: dx,dy,dz,dthetx,dthety,dthetz
    // 8149: from,xie1,xie2,eta1,eta2

    if (ioflag == 0) {
        fortran_read_list(iunit5);
        fprintf(fu25, "   int    to  xie1  xie2  eta1  eta2   nfb  from  from  from  from   ...\n");
    } else if (ioflag == 1) {
        fortran_read_list(iunit5);
        fortran_read_list(iunit5);
        fprintf(fu25, "   int    to     xie1     xie2     eta1     eta2      nfb\n");
        fprintf(fu25, "        from     xie1     xie2     eta1     eta2\n");
    }

    // do 1703 n=1,nintr
    for (n = 1; n <= nintr; n++) {
        if (ioflag == 2) {
            fortran_read_list(iunit5);
            fortran_read_list(iunit5, &int_val, &ito, &xi1, &xi2, &et1, &et2, &nfb);
        } else if (ioflag == 1) {
            fortran_read_list(iunit5, &int_val, &ito, &xi1, &xi2, &et1, &et2, &nfb);
        } else if (ioflag == 0) {
            // read(iunit5,*) int,ito,xi1,xi2,et1,et2,nfb,(ifrom(l),l=1,nfb)
            // This is ONE Fortran read statement. First read the 7 scalars,
            // then use fortran_read_list_array for the implied-DO ifrom values.
            fortran_read_list(iunit5, &int_val, &ito, &xi1, &xi2, &et1, &et2, &nfb);
            // ifrom(1..nfb) are on the same record (or continuation lines)
            fortran_read_list_array(iunit5, &ifrom(1), nfb);
            for (l = 1; l <= nfb; l++) {
                xif1(l)            = 0;
                xif2(l)            = 0;
                etf1(l)            = 0;
                etf2(l)            = 0;
                dx(int_val, l)     = 0.0;
                dy(int_val, l)     = 0.0;
                dz(int_val, l)     = 0.0;
                dthetx(int_val, l) = 0.0;
                dthety(int_val, l) = 0.0;
                dthetz(int_val, l) = 0.0;
            }
        }


        // make sure parameter msub1 is big enough
        if (nfb > msub1) {
            fprintf(fu25, " msub1 is too small; should be least: %2d\n", nfb);
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        nfbmax = std::max(nfbmax, nfb);

        if (ioflag == 1) {
            // Each iteration reads one line: ifrom(l),xif1(l),xif2(l),etf1(l),etf2(l)
            // This is a genuine do-loop with one read per l (not an implied-DO)
            l = 1;
            while (l <= nfb) {
                fortran_read_list(iunit5, &ifrom(l), &xif1(l), &xif2(l), &etf1(l), &etf2(l));
                dx(int_val, l)     = 0.0;
                dy(int_val, l)     = 0.0;
                dz(int_val, l)     = 0.0;
                dthetx(int_val, l) = 0.0;
                dthety(int_val, l) = 0.0;
                dthetz(int_val, l) = 0.0;
                l++;
            }
        } else if (ioflag == 2) {
            for (l = 1; l <= nfb; l++) {
                fortran_read_list(iunit5);
                fortran_read_list(iunit5, &ifrom(l), &xif1(l), &xif2(l), &etf1(l), &etf2(l));
                fortran_read_list(iunit5);
                // read(iunit5,*) (realval(i),i=1,6)
                fortran_read_list(iunit5, &realval[1], &realval[2], &realval[3],
                                  &realval[4], &realval[5], &realval[6]);
                dx(int_val, l)     = realval[1];
                dthetx(int_val, l) = realval[4];
                if (ialph > 0) {
                    dy(int_val, l)     = -realval[3];
                    dthety(int_val, l) = -realval[6];
                    dz(int_val, l)     = realval[2];
                    dthetz(int_val, l) = realval[5];
                } else {
                    dy(int_val, l)     = realval[2];
                    dthety(int_val, l) = realval[5];
                    dz(int_val, l)     = realval[3];
                    dthetz(int_val, l) = realval[6];
                }
            }
        }

        if (iindx(int_val, 1) != -99) {
            char grdmov_buf[61];
            std::memcpy(grdmov_buf, grdmov, 60);
            grdmov_buf[60] = '\0';
            fprintf(fu11, " program terminated in dynamic patching routines - see file %60s\n", grdmov_buf);
            fprintf(fu25, "\n stopping...attempting to set data for int %d more than once\n", int_val);
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        iindx(int_val, 1) = nfb;
        // do 1701 l=1,nfb
        for (l = 1; l <= nfb; l++) {
            iindx(int_val, l + 1)       = ifrom(l) / 100;
            iindx(int_val, l + nfb + 2) = ifrom(l) - iindx(int_val, l + 1) * 100;
        }

        iindx(int_val, nfb + 2)   = ito / 100;
        nbl                        = nblg(iindx(int_val, nfb + 2));
        iindx(int_val, 2 * nfb + 3) = ito - iindx(int_val, nfb + 2) * 100;



        // set full ranges if zeros are input
        if (iindx(int_val, 2 * nfb + 3) / 10 == 1) {
            if (xi1 == 0 && xi2 == 0) {
                xi1 = 1;
                xi2 = jdimg(nbl);
            }
            if (et1 == 0 && et2 == 0) {
                et1 = 1;
                et2 = kdimg(nbl);
            }
        }
        if (iindx(int_val, 2 * nfb + 3) / 10 == 2) {
            if (xi1 == 0 && xi2 == 0) {
                xi1 = 1;
                xi2 = kdimg(nbl);
            }
            if (et1 == 0 && et2 == 0) {
                et1 = 1;
                et2 = idimg(nbl);
            }
        }
        if (iindx(int_val, 2 * nfb + 3) / 10 == 3) {
            if (xi1 == 0 && xi2 == 0) {
                xi1 = 1;
                xi2 = jdimg(nbl);
            }
            if (et1 == 0 && et2 == 0) {
                et1 = 1;
                et2 = idimg(nbl);
            }
        }

        // set full ranges if zeros are input for search range in "from" block(s)
        // do 1705 l=1,nfb
        for (l = 1; l <= nfb; l++) {
            mbl = nblg(iindx(int_val, l + 1));
            if (iindx(int_val, l + nfb + 2) / 10 == 1) {
                if (xif1(l) == 0 && xif2(l) == 0) {
                    xif1(l) = 1;
                    xif2(l) = jdimg(mbl);
                }
                if (etf1(l) == 0 && etf2(l) == 0) {
                    etf1(l) = 1;
                    etf2(l) = kdimg(mbl);
                }
            }
            if (iindx(int_val, l + nfb + 2) / 10 == 2) {
                if (xif1(l) == 0 && xif2(l) == 0) {
                    xif1(l) = 1;
                    xif2(l) = kdimg(mbl);
                }
                if (etf1(l) == 0 && etf2(l) == 0) {
                    etf1(l) = 1;
                    etf2(l) = idimg(mbl);
                }
            }
            if (iindx(int_val, l + nfb + 2) / 10 == 3) {
                if (xif1(l) == 0 && xif2(l) == 0) {
                    xif1(l) = 1;
                    xif2(l) = jdimg(mbl);
                }
                if (etf1(l) == 0 && etf2(l) == 0) {
                    etf1(l) = 1;
                    etf2(l) = idimg(mbl);
                }
            }
        }

        if (ioflag == 0 || ioflag == 1) {
            // write(iunit25,38) int,ito,xi1,xi2,et1,et2,nfb
            fprintf(fu25, "%6d%6d%9d%9d%9d%9d%9d\n",
                    int_val, ito, xi1, xi2, et1, et2, nfb);
            for (l = 1; l <= nfb; l++) {
                // write(iunit25,40) ifrom(l),xif1(l),xif2(l),etf1(l),etf2(l)
                // format 40: 6x,i6,10i9
                fprintf(fu25, "      %6d%9d%9d%9d%9d\n",
                        ifrom(l), xif1(l), xif2(l), etf1(l), etf2(l));
            }
        } else {
            // ioflag==2
            fprintf(fu25, "   int    to     xie1     xie2     eta1     eta2      nfb\n");
            fprintf(fu25, "%6d%6d%9d%9d%9d%9d%9d\n",
                    int_val, ito, xi1, xi2, et1, et2, nfb);
            for (l = 1; l <= nfb; l++) {
                fprintf(fu25, "        from     xie1     xie2     eta1     eta2\n");
                fprintf(fu25, "      %6d%9d%9d%9d%9d\n",
                        ifrom(l), xif1(l), xif2(l), etf1(l), etf2(l));
                fprintf(fu25, "        dx       dy       dz   dthetx   dthety   dthetz\n");
                // write(iunit25,39) real(dx(int,l)),...
                // format 39: 3x,6f9.4
                fprintf(fu25, "   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
                        (float)dx(int_val, l), (float)dy(int_val, l),
                        (float)dz(int_val, l), (float)dthetx(int_val, l),
                        (float)dthety(int_val, l), (float)dthetz(int_val, l));
                dthetx(int_val, l) = dthetx(int_val, l) / (double)radtodeg;
                dthety(int_val, l) = dthety(int_val, l) / (double)radtodeg;
                dthetz(int_val, l) = dthetz(int_val, l) / (double)radtodeg;
            }
        }



        // check input ranges of xie and eta
        if (iindx(int_val, 2 * nfb + 3) / 10 == 1) {
            if (xi1 < 1 || xi1 > jdimg(nbl)) {
                fprintf(fu25, " incorrect xie range input on this i=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (xi2 < 1 || xi2 > jdimg(nbl)) {
                fprintf(fu25, " incorrect xie range input on this i=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (xi2 <= xi1) {
                fprintf(fu25, " incorrect xie range input on this i=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (et1 < 1 || et1 > kdimg(nbl)) {
                fprintf(fu25, " incorrect eta range input on this i=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (et2 < 1 || et2 > kdimg(nbl)) {
                fprintf(fu25, " incorrect eta range input on this i=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (et2 <= et1) {
                fprintf(fu25, " incorrect eta range input on this i=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        }

        if (iindx(int_val, 2 * nfb + 3) / 10 == 2) {
            if (xi1 < 1 || xi1 > kdimg(nbl)) {
                fprintf(fu25, " incorrect xie range input on this j=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (xi2 < 1 || xi2 > kdimg(nbl)) {
                fprintf(fu25, " incorrect xie range input on this j=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (xi2 <= xi1) {
                fprintf(fu25, " incorrect xie range input on this j=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (et1 < 1 || et1 > idimg(nbl)) {
                fprintf(fu25, " incorrect eta range input on this j=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (et2 < 1 || et2 > idimg(nbl)) {
                fprintf(fu25, " incorrect eta range input on this j=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (et2 <= et1) {
                fprintf(fu25, " incorrect eta range input on this j=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        }

        if (iindx(int_val, 2 * nfb + 3) / 10 == 3) {
            if (xi1 < 1 || xi1 > jdimg(nbl)) {
                fprintf(fu25, " incorrect xie range input on this k=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (xi2 < 1 || xi2 > jdimg(nbl)) {
                fprintf(fu25, " incorrect xie range input on this k=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (xi2 <= xi1) {
                fprintf(fu25, " incorrect xie range input on this k=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (et1 < 1 || et1 > idimg(nbl)) {
                fprintf(fu25, " incorrect eta range input on this k=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (et2 < 1 || et2 > idimg(nbl)) {
                fprintf(fu25, " incorrect eta range input on this k=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (et2 <= et1) {
                fprintf(fu25, " incorrect eta range input on this k=constant \"to\" face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        }



        // check input search ranges of xie and eta - "from" block(s)
        // do 1716 l=1,nfb
        for (l = 1; l <= nfb; l++) {
            mbl = nblg(iindx(int_val, l + 1));
            if (iindx(int_val, l + nfb + 2) / 10 == 1) {
                if (xif1(l) < 1 || xif1(l) > jdimg(mbl)) {
                    fprintf(fu25, " incorrect xie search range input on the i=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (xif2(l) < 1 || xif2(l) > jdimg(mbl)) {
                    fprintf(fu25, " incorrect xie search range input on the i=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (xif2(l) <= xif1(l)) {
                    fprintf(fu25, " incorrect xie search range input on the i=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (etf1(l) < 1 || etf1(l) > kdimg(mbl)) {
                    fprintf(fu25, " incorrect eta search range input on the i=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (etf2(l) < 1 || etf2(l) > kdimg(mbl)) {
                    fprintf(fu25, " incorrect eta search range input on the i=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (etf2(l) <= etf1(l)) {
                    fprintf(fu25, " incorrect eta search range input on the i=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }

            if (iindx(int_val, l + nfb + 2) / 10 == 2) {
                if (xif1(l) < 1 || xif1(l) > kdimg(mbl)) {
                    fprintf(fu25, " incorrect xie search range input on the j=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (xif2(l) < 1 || xif2(l) > kdimg(mbl)) {
                    fprintf(fu25, " incorrect xie search range input on the j=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (xif2(l) <= xif1(l)) {
                    fprintf(fu25, " incorrect xie search range input on the j=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (etf1(l) < 1 || etf1(l) > idimg(mbl)) {
                    fprintf(fu25, " incorrect eta search range input on the j=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (etf2(l) < 1 || etf2(l) > idimg(mbl)) {
                    fprintf(fu25, " incorrect eta search range input on the j=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (etf2(l) <= etf1(l)) {
                    fprintf(fu25, " incorrect eta search range input on the j=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }

            if (iindx(int_val, l + nfb + 2) / 10 == 3) {
                if (xif1(l) < 1 || xif1(l) > jdimg(mbl)) {
                    fprintf(fu25, " incorrect xie search range input on the k=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (xif2(l) < 1 || xif2(l) > jdimg(mbl)) {
                    fprintf(fu25, " incorrect xie search range input on the k=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (xif2(l) <= xif1(l)) {
                    fprintf(fu25, " incorrect xie search range input on the k=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (etf1(l) < 1 || etf1(l) > idimg(mbl)) {
                    fprintf(fu25, " incorrect eta search range input on the k=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (etf2(l) < 1 || etf2(l) > idimg(mbl)) {
                    fprintf(fu25, " incorrect eta search range input on the k=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (etf2(l) <= etf1(l)) {
                    fprintf(fu25, " incorrect eta search range input on the k=constant \"from\" face%d\n", l);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }
        } // end do 1716

        iindx(int_val, 2 * nfb + 6) = xi1;
        iindx(int_val, 2 * nfb + 7) = xi2;
        iindx(int_val, 2 * nfb + 8) = et1;
        iindx(int_val, 2 * nfb + 9) = et2;
        // do 1707 l=1,nfb
        for (l = 1; l <= nfb; l++) {
            iindx(int_val, 2 * nfb + 9 + l) = xif1(l);
            iindx(int_val, 3 * nfb + 9 + l) = xif2(l);
            iindx(int_val, 4 * nfb + 9 + l) = etf1(l);
            iindx(int_val, 5 * nfb + 9 + l) = etf2(l);
        }

        iindx(int_val, 2 * nfb + 4) = (xi2 - xi1) * (et2 - et1);
        if (int_val == 1) {
            iindx(int_val, 2 * nfb + 5) = 1;
        } else {
            nfb1 = iindx(int_val - 1, 1);
            iindx(int_val, 2 * nfb + 5) = iindx(int_val - 1, 2 * nfb1 + 5)
                                         + iindx(int_val - 1, 2 * nfb1 + 4);
        }
    } // end do 1703



    // Lines 685-765: imode==0 summary by grids/levels
    if (imode == 0) {
        fprintf(fu25, "\n %s\n", "SUMMARY BY GRIDS");
        // format 117: 3x,4hgrid,2x,5hlevel,2x,5hblock,3x,4hjdim,3x,4hkdim,3x,4hidim,3x,9hgrid pts.
        fprintf(fu25, "   grid  level  block   jdim   kdim   idim   grid pts.\n");

        lfem        = 0;
        lfgm_local  = ncgmax + 1;
        lcgm        = 1;
        igptot      = 0;
        // do 105 igrid=1,ngrid
        for (igrid = 1; igrid <= ngrid; igrid++) {
            int iem  = iemg(igrid);
            nbl      = nblg(igrid);
            ncg      = ncgg(igrid);
            igridg(nbl) = igrid;
            levelg(nbl) = iem + lfgm_local;
            if (iem == 0) lfgm_local = std::max(lfgm_local, levelg(nbl));
            if (iem > 0)  lfem       = std::max(lfem, levelg(nbl));
            igpts  = jdimg(nbl) * kdimg(nbl) * idimg(nbl);
            igptot = igptot + igpts;
            // format 7: 6i7,i12
            fprintf(fu25, "%7d%7d%7d%7d%7d%7d%12d\n",
                    igrid, levelg(nbl), nbl,
                    jdimg(nbl), kdimg(nbl), idimg(nbl), igpts);
            if (ncg > 0 && iem == 0) {
                // do 104 n=1,ncg
                for (int nn = 1; nn <= ncg; nn++) {
                    nbl = nbl + 1;
                    igridg(nbl) = igrid;
                    levelg(nbl) = levelg(nbl - 1) - 1;
                    igpts  = jdimg(nbl) * kdimg(nbl) * idimg(nbl);
                    igptot = igptot + igpts;
                    fprintf(fu25, "%7d%7d%7d%7d%7d%7d%12d\n",
                            igrid, levelg(nbl), nbl,
                            jdimg(nbl), kdimg(nbl), idimg(nbl), igpts);
                }
            }
        }
        // format 2021: /,37x,5hTOTAL,i12
        fprintf(fu25, "\n                                     TOTAL%12d\n", igptot);

        fprintf(fu25, "\n %s\n", "SUMMARY BY LEVELS");
        // format 110: 2x,5hlevel,3x,4hgrid,2x,5hblock
        fprintf(fu25, "  level   grid  block\n");
        lf = lfem;
        if (lfem == 0) lf = lfgm_local;
        // do 115 levelc=lf,lcgm,-1
        for (levelc = lf; levelc >= lcgm; levelc--) {
            // do 114 nbl=1,nblock
            for (nbl_local = 1; nbl_local <= nblock; nbl_local++) {
                if (levelc != levelg(nbl_local)) continue;
                igrid = igridg(nbl_local);
                // format 7: 6i7,i12 (but only 3 values here)
                fprintf(fu25, "%7d%7d%7d\n", levelg(nbl_local), igrid, nbl_local);
            }
        }

        if (mseq > (ncgmax + 1)) {
            // format 153: 1x,12hmseq,ncgmax=,2i4
            fprintf(fu25, " mseq,ncgmax=%4d%4d\n", mseq, ncgmax);
            // format 155: 1x,28herror in input, mseq, ncgmax,2i5
            fprintf(fu25, " error in input, mseq, ncgmax%5d%5d\n", mseq, ncgmax);
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        // format 4008: /
        fprintf(fu25, "\n");
        if (lfem != 0) {
            // format 4009: 31h level of finest embedded mesh=,i3
            fprintf(fu25, " level of finest embedded mesh=%3d\n", lfem);
        }
        // format 4019: 31h level of coarsest global mesh=,i3
        fprintf(fu25, " level of coarsest global mesh=%3d\n", lcgm);
        // format 4029: 31h level of finest global mesh  =,i3
        fprintf(fu25, " level of finest global mesh  =%3d\n", lfgm_local);

        fprintf(fu25, "\n %s\n", "SUMMARY OF GRID SEQUENCES");
        // format 9: 1x,8hsequence,2x,14hstarting level,4x,12hending level
        fprintf(fu25, " sequence  starting level    ending level\n");
        // do 19 m=1,mseq
        for (m = 1; m <= mseq; m++) {
            // format 29: 2x,i7,2i16
            fprintf(fu25, "  %7d%16d%16d\n", m, levelt[m - 1], levelb[m - 1]);
        }
    } // end if (imode==0)



    // Lines 769-946: generate additional patch data for coarser blocks
    ntemp  = nintr;
    ifirst = 1;
    // do 8400 n=1,nintr
    for (n = 1; n <= nintr; n++) {
        nfb = iindx(n, 1);
        ifiner(n) = 0;
        if (ncgg(iindx(n, nfb + 2)) > 0) {
            if (ifirst == 1) {
                if (ioflag == 2) {
                    fprintf(fu25, "\n REARRANGED BLOCK NUMBERS AND DYNAMIC PATCH DATA GENERATED FOR COARSER MESHES\n");
                    fprintf(fu25, "           block references to CFL3D internal ordering\n");
                    fprintf(fu25, "                * ---> coarser level interpolation\n\n");
                } else {
                    fprintf(fu25, "\n REARRANGED BLOCK NUMBERS AND PATCH DATA GENERATED FOR COARSER MESHES\n");
                    fprintf(fu25, "           block references to CFL3D internal ordering\n");
                    fprintf(fu25, "                * ---> coarser level interpolation\n\n");
                }
                ifirst = 0;
            }
            ito    = (iindx(n, nfb + 2) * 100) + iindx(n, 2 * nfb + 3);
            itonew = nblg(iindx(n, nfb + 2)) * 100 + iindx(n, 2 * nfb + 3);
            xi1 = iindx(n, 2 * nfb + 6);
            xi2 = iindx(n, 2 * nfb + 7);
            et1 = iindx(n, 2 * nfb + 8);
            et2 = iindx(n, 2 * nfb + 9);
            for (i = 1; i <= nfb; i++) {
                ifrom(i) = nblg(iindx(n, i + 1)) * 100 + iindx(n, i + nfb + 2);
            }
            fprintf(fu25, "   int    to     xie1     xie2     eta1     eta2      nfb\n");
            fprintf(fu25, "%6d%6d%9d%9d%9d%9d%9d\n", n, itonew, xi1, xi2, et1, et2, nfb);
            for (l = 1; l <= nfb; l++) {
                fprintf(fu25, "        from     xie1     xie2     eta1     eta2\n");
                fprintf(fu25, "      %6d%9d%9d%9d%9d\n",
                        ifrom(l), xif1(l), xif2(l), etf1(l), etf2(l));
                if (ioflag == 2) {
                    fprintf(fu25, "        dx       dy       dz   dthetx   dthety   dthetz\n");
                    fprintf(fu25, "   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
                            (float)dx(n, l), (float)dy(n, l), (float)dz(n, l),
                            (float)dthetx(n, l) * (float)radtodeg,
                            (float)dthety(n, l) * (float)radtodeg,
                            (float)dthetz(n, l) * (float)radtodeg);
                }
            }

            // do 8200 m=1,ncgg(iindx(n,nfb+2))
            for (m = 1; m <= ncgg(iindx(n, nfb + 2)); m++) {
                ntemp = ntemp + 1;
                lprev = ntemp - 1;
                if (m == 1) {
                    ifiner(ntemp) = n;
                } else {
                    ifiner(ntemp) = ntemp - 1;
                }
                iifit(ntemp)  = iifit(n);
                iic0(ntemp)   = iic0(n);
                iiorph(ntemp) = iiorph(n);
                iitoss(ntemp) = iitoss(n);
                llimit(ntemp) = llimit(n);
                iitmax(ntemp) = iitmax(n);
                mmcxie(ntemp) = mmcxie(n);
                mmceta(ntemp) = mmceta(n);
                for (l = 1; l <= nfb; l++) {
                    dx(ntemp, l)     = dx(n, l);
                    dy(ntemp, l)     = dy(n, l);
                    dz(ntemp, l)     = dz(n, l);
                    dthetx(ntemp, l) = dthetx(n, l);
                    dthety(ntemp, l) = dthety(n, l);
                    dthetz(ntemp, l) = dthetz(n, l);
                }
                iindx(ntemp, 1)           = nfb;
                iindx(ntemp, nfb + 2)     = nblg(iindx(n, nfb + 2)) + m;
                iindx(ntemp, nfb + 3)     = iindx(n, nfb + 3);
                iindx(ntemp, 2 * nfb + 3) = iindx(n, 2 * nfb + 3);
                nfb1      = iindx(lprev, 1);
                nbl_local = iindx(ntemp, nfb + 2);

                xi1 = iindx(n, 2 * nfb + 6);
                xi2 = iindx(n, 2 * nfb + 7);
                et1 = iindx(n, 2 * nfb + 8);
                et2 = iindx(n, 2 * nfb + 9);

                iindx(ntemp, 2 * nfb + 6) = xi1 / (2 * m) + 1;
                iindx(ntemp, 2 * nfb + 7) = xi2 / (2 * m) + 1;
                iindx(ntemp, 2 * nfb + 8) = et1 / (2 * m) + 1;
                iindx(ntemp, 2 * nfb + 9) = et2 / (2 * m) + 1;

                if (iindx(n, 2 * nfb + 6) <= 2) iindx(ntemp, 2 * nfb + 6) = iindx(n, 2 * nfb + 6);
                if (iindx(n, 2 * nfb + 7) <= 2) iindx(ntemp, 2 * nfb + 7) = iindx(n, 2 * nfb + 7);
                if (iindx(n, 2 * nfb + 8) <= 2) iindx(ntemp, 2 * nfb + 8) = iindx(n, 2 * nfb + 8);
                if (iindx(n, 2 * nfb + 9) <= 2) iindx(ntemp, 2 * nfb + 9) = iindx(n, 2 * nfb + 9);

                xi1 = iindx(ntemp, 2 * nfb + 6);
                xi2 = iindx(ntemp, 2 * nfb + 7);
                et1 = iindx(ntemp, 2 * nfb + 8);
                et2 = iindx(ntemp, 2 * nfb + 9);



                // do 8202 l=1,nfb
                for (l = 1; l <= nfb; l++) {
                    xif11 = iindx(n, 2 * nfb + 9 + l);
                    xif22 = iindx(n, 3 * nfb + 9 + l);
                    etf11 = iindx(n, 4 * nfb + 9 + l);
                    etf22 = iindx(n, 5 * nfb + 9 + l);

                    iindx(ntemp, 2 * nfb + 9 + l) = xif11 / (2 * m) + 1;
                    iindx(ntemp, 3 * nfb + 9 + l) = xif22 / (2 * m) + 1;
                    iindx(ntemp, 4 * nfb + 9 + l) = etf11 / (2 * m) + 1;
                    iindx(ntemp, 5 * nfb + 9 + l) = etf22 / (2 * m) + 1;

                    if (iindx(n, 2 * nfb + 9 + l) <= 2)
                        iindx(ntemp, 2 * nfb + 9 + l) = iindx(n, 2 * nfb + 9 + l);
                    if (iindx(n, 3 * nfb + 9 + l) <= 2)
                        iindx(ntemp, 3 * nfb + 9 + l) = iindx(n, 3 * nfb + 9 + l);
                    if (iindx(n, 4 * nfb + 9 + l) <= 2)
                        iindx(ntemp, 4 * nfb + 9 + l) = iindx(n, 4 * nfb + 9 + l);
                    if (iindx(n, 5 * nfb + 9 + l) <= 2)
                        iindx(ntemp, 5 * nfb + 9 + l) = iindx(n, 5 * nfb + 9 + l);
                }

                iindx(ntemp, 2 * nfb + 4) = (xi2 - xi1) * (et2 - et1);
                iindx(ntemp, 2 * nfb + 5) = iindx(lprev, 2 * nfb1 + 4)
                                           + iindx(lprev, 2 * nfb1 + 5);
                mxxchk = iindx(ntemp, 2 * nfb + 4) + iindx(ntemp, 2 * nfb + 5) - 1;

                if (icall > 0) {
                    if (mxxchk > mxxe) {
                        char grdmov_buf[61];
                        std::memcpy(grdmov_buf, grdmov, 60);
                        grdmov_buf[60] = '\0';
                        fprintf(fu11, " program terminated in dynamic patching routines - see file %60s\n", grdmov_buf);
                        fprintf(fu25, "  stopping...parameter mxxe too small\n");
                        fprintf(fu25, "  must make  mxxe at least %d\n", mxxchk);
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                }

                // do 8100 i=1,nfb
                for (i = 1; i <= nfb; i++) {
                    iindx(ntemp, i + 1)       = nblg(iindx(n, i + 1)) + m;
                    iindx(ntemp, nfb + i + 2) = iindx(n, nfb + i + 2);
                    if (ncgg(iindx(n, i + 1)) < m)
                        iindx(ntemp, i + 1) = nblg(iindx(n, i + 1))
                                            + ncgg(iindx(n, i + 1));
                    ifrom(i) = (iindx(ntemp, i + 1) * 100) + iindx(ntemp, i + nfb + 2);
                }
                ito = (iindx(ntemp, nfb + 2) * 100) + iindx(ntemp, 2 * nfb + 3);
                // write(iunit25,8147)
                fprintf(fu25, "   int    to     xie1     xie2     eta1     eta2      nfb\n");
                // write(iunit25,37) ntemp,ito,xi1,xi2,et1,et2,nfb  -- format 37: 2h *,i4,i6,10i9
                fprintf(fu25, " *%4d%6d%9d%9d%9d%9d%9d\n",
                        ntemp, ito, xi1, xi2, et1, et2, nfb);
                // do 9302 l=1,nfb
                for (l = 1; l <= nfb; l++) {
                    fprintf(fu25, "        from     xie1     xie2     eta1     eta2\n");
                    fprintf(fu25, "      %6d%9d%9d%9d%9d\n",
                            ifrom(l),
                            iindx(ntemp, 2 * nfb + 9 + l),
                            iindx(ntemp, 3 * nfb + 9 + l),
                            iindx(ntemp, 4 * nfb + 9 + l),
                            iindx(ntemp, 5 * nfb + 9 + l));
                    if (ioflag == 2) {
                        fprintf(fu25, "        dx       dy       dz   dthetx   dthety   dthetz\n");
                        fprintf(fu25, "   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
                                (float)dx(n, l), (float)dy(n, l), (float)dz(n, l),
                                (float)dthetx(n, l) * (float)radtodeg,
                                (float)dthety(n, l) * (float)radtodeg,
                                (float)dthetz(n, l) * (float)radtodeg);
                    }
                }
            } // end do 8200
        } // end if ncgg > 0

        // renumber "to" and "from" block numbers
        iindx(n, nfb + 2) = nblg(iindx(n, nfb + 2));
        for (i = 1; i <= nfb; i++) {
            iindx(n, i + 1) = nblg(iindx(n, i + 1));
        }
    } // end do 8400

    nintr = ntemp;



    // Lines 950-1068: check parameter sizes, multigridability
    mpa = 0;
    // do 1713 n=1,nintr
    for (n = 1; n <= nintr; n++) {
        if (n > 1) nfb1 = iindx(n - 1, 1);
        nfb = iindx(n, 1);
        nbl_local = iindx(n, nfb + 2);
        xi1 = iindx(n, 2 * nfb + 6);
        xi2 = iindx(n, 2 * nfb + 7);
        et1 = iindx(n, 2 * nfb + 8);
        et2 = iindx(n, 2 * nfb + 9);
        int id2_v  = idimg(nbl_local);
        int jd2_v  = jdimg(nbl_local);
        int kd2_v  = kdimg(nbl_local);
        itop = iindx(n, 2 * nfb + 3) / 10;
        if (itop == 1) {
            mpa = std::max({jd2_v, kd2_v, mpa});
        }
        if (itop == 2) {
            mpa = std::max({id2_v, kd2_v, mpa});
        }
        if (itop == 3) {
            mpa = std::max({id2_v, jd2_v, mpa});
        }
        iindx(n, 2 * nfb + 4) = (xi2 - xi1) * (et2 - et1);
        // do 1714 l=1,nfb
        for (l = 1; l <= nfb; l++) {
            int mbl_v  = iindx(n, l + 1);
            mbl_v      = iindx(n, l + 1);
            int id1_v  = idimg(mbl_v);
            int jd1_v  = jdimg(mbl_v);
            int kd1_v  = kdimg(mbl_v);
            itop = (iindx(n, l + nfb + 2) - iindx(n, l + nfb + 2) / 100 * 100) / 10;
            if (itop == 1) {
                mpa = std::max({jd1_v, kd1_v, mpa});
            }
            if (itop == 2) {
                mpa = std::max({id1_v, kd1_v, mpa});
            }
            if (itop == 3) {
                mpa = std::max({id1_v, jd1_v, mpa});
            }

            lmptch = std::max((int)lmptch, mpa);
        } // end do 1714

        if (icall > 0) {
            if (mpa > mptch) {
                // format 1929: 1x,31hmptch is too small; should be: ,i4
                fprintf(fu25, " mptch is too small; should be: %4d\n", mpa);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        }

        if (n == 1) {
            iindx(n, 2 * nfb + 5) = 1;
        } else {
            iindx(n, 2 * nfb + 5) = iindx(n - 1, 2 * nfb1 + 4) + iindx(n - 1, 2 * nfb1 + 5);
        }
        mxxchk = iindx(nintr, 2 * nfb + 4) + iindx(nintr, 2 * nfb + 5) - 1;
    } // end do 1713

    if (icall > 0) {
        if (mxxchk > mxxe) {
            char grdmov_buf[61];
            std::memcpy(grdmov_buf, grdmov, 60);
            grdmov_buf[60] = '\0';
            fprintf(fu11, " program terminated in dynamic patching routines - see file %60s\n", grdmov_buf);
            fprintf(fu25, "  stopping...parameter mxxe too small\n");
            fprintf(fu25, "  must make mxxe at least %d\n", mxxchk);
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }

    // check multigridability - only check "to" side
    istop = 0;
    for (n = 1; n <= nintr; n++) {
        nfb = iindx(n, 1);
        nbl_local = iindx(n, nfb + 2);
        xi1 = iindx(n, 2 * nfb + 6);
        xi2 = iindx(n, 2 * nfb + 7);
        et1 = iindx(n, 2 * nfb + 8);
        et2 = iindx(n, 2 * nfb + 9);
        if (levelg(nbl_local) > 1) {
            if ((float)(xi1 / 2) == (float)xi1 / 2.f) {
                fprintf(fu25, " begining \"to\" xie index not multigridable for interpolation %4d\n", n);
                istop = 1;
            }
            if ((float)(xi2 / 2) == (float)xi2 / 2.f) {
                fprintf(fu25, " ending \"to\" xie index not multigridable for interpolation %4d\n", n);
                istop = 1;
            }
            // skip eta multigrid checks if 2d
            if (idimg(nbl_local) == 2) {
                if (et1 != 1 || et2 != 2) {
                    fprintf(fu25, " for 2D cases must have eta1 = 1 and eta2 = 2\n");
                    istop = 1;
                }
            } else {
                if ((float)(et1 / 2) == (float)et1 / 2.f) {
                    fprintf(fu25, " begining \"to\" eta index not multigridable for interpolation %4d\n", n);
                    istop = 1;
                }
                if ((float)(et2 / 2) == (float)et2 / 2.f) {
                    fprintf(fu25, " ending \"to\" eta index not multigridable for interpolation %4d\n", n);
                    istop = 1;
                }
            }
        }
    }
    if (istop > 0) {
        fprintf(fu25, " stopping\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }



    // Lines 1070-1100: debug iindx check (ichk=0 so never executes)
    ichk = 0;
    if (ichk > 0) {
        fprintf(fu8, "nintr: %d\n", nintr);
        fprintf(fu8, "\n");
        // do 6969 n=1,nintr
        for (n = 1; n <= nintr; n++) {
            fprintf(fu8, "\n");
            fprintf(fu8, "\n");
            fprintf(fu8, " interp. no.: %d\n", n);
            fprintf(fu8, "\n");
            fprintf(fu8, " nfb: %d\n", iindx(n, 1));
            nfb = iindx(n, 1);
            fprintf(fu8, " from block(s):");
            for (l = 1; l <= nfb; l++) {
                fprintf(fu8, " %d", iindx(n, l + 1));
            }
            fprintf(fu8, "\n");
            fprintf(fu8, "  to: %d\n", iindx(n, nfb + 2));
            fprintf(fu8, " topology (from):");
            for (l = 1; l <= nfb; l++) {
                fprintf(fu8, " %d", iindx(n, l + nfb + 2));
            }
            fprintf(fu8, "\n");
            fprintf(fu8, " topology (  to): %d\n", iindx(n, 2 * nfb + 3));
            fprintf(fu8, " number of points: %d\n", iindx(n, 2 * nfb + 4));
            fprintf(fu8, " starting index: %d\n", iindx(n, 2 * nfb + 5));
            fprintf(fu8, " xie range: %d %d\n", iindx(n, 2 * nfb + 6), iindx(n, 2 * nfb + 7));
            fprintf(fu8, " eta range: %d %d\n", iindx(n, 2 * nfb + 8), iindx(n, 2 * nfb + 9));
            nfb = iindx(n, 1);
            // do 6968 l=1,nfb
            for (l = 1; l <= nfb; l++) {
                fprintf(fu8, " xie search range in from block %d : %d %d\n",
                        l, iindx(n, 2 * nfb + 9 + l), iindx(n, 3 * nfb + 9 + l));
                fprintf(fu8, " eta search range, from block %d : %d %d\n",
                        l, iindx(n, 4 * nfb + 9 + l), iindx(n, 5 * nfb + 9 + l));
            }
        }
    }

    // Lines 1105-1127: adjust patch parameters (icall==0 branch)
    if (icall == 0) {
        lmsub1 = nfbmax;
        lnsub1 = std::max((int)lnsub1, (int)lmsub1);
        lintmx = nintr;
        if (lintmax > 1) {
            lintmax = lintmax + lintmx;
        } else {
            lintmax = lintmx;
        }
        lmxxe = 0;
        ntot = 0;
        for (int int_i = 1; int_i <= nintr; int_i++) {
            nfbl = iindx(int_i, 1);
            npts = iindx(int_i, 2 * nfbl + 4);
            nst  = iindx(int_i, 2 * nfbl + 5);
            ntot = nst + npts - 1;
        }
        lmxxe = ntot;
        if (lmaxxe > 1) {
            lmaxxe = lmaxxe + lmxxe;
        } else {
            lmaxxe = lmxxe;
        }
    }

    // 999 continue
    return;
}

} // namespace global2_ns
