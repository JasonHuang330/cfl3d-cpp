// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "global0.h"
#include "global.h"
#include "readkey.h"
#include "termn8.h"
#include "my_flush.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include "runtime/fortran_array.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace global0_ns {

// Forward declaration with int irw by value (allows passing integer literals)
static void echoinp(int& iunit5, int& iunit11, int irw);


// Helper: skip one Fortran unformatted record from file f
static void skip_unformatted_record(FILE* f) {
    if (!f) return;
    int32_t reclen = 0;
    if (fread(&reclen, sizeof(int32_t), 1, f) != 1) return;
    if (reclen < 0) reclen = -reclen;
    fseek(f, (long)reclen, SEEK_CUR);
    int32_t reclen2 = 0;
    fread(&reclen2, sizeof(int32_t), 1, f);
}

// Helper: read one Fortran unformatted record into buf (up to bufsize bytes)
static int read_unformatted_record(FILE* f, void* buf, int bufsize) {
    if (!f) return 0;
    int32_t reclen = 0;
    if (fread(&reclen, sizeof(int32_t), 1, f) != 1) return 0;
    if (reclen < 0) reclen = -reclen;
    int toread = (reclen < bufsize) ? reclen : bufsize;
    int nr = (int)fread(buf, 1, toread, f);
    if (reclen > bufsize) fseek(f, (long)(reclen - bufsize), SEEK_CUR);
    int32_t reclen2 = 0;
    fread(&reclen2, sizeof(int32_t), 1, f);
    return nr;
}

void global(int& myid, int& maxbl, int& maxgr, int& maxseg, int& maxcs, int& nplots, int& mxbli, FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, FortranArray2DRef<int> nblk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, FortranArray2DRef<double> rkap0g, int& nblock, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<int> idiagg, FortranArray1DRef<int> nblcg, FortranArray2DRef<int> idegg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray2DRef<int> mit, FortranArray1DRef<int> ilamlog, FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamlog, FortranArray1DRef<int> klamhig, FortranArray2DRef<int> iwfg, FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx, FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz, FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr, FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg, FortranArray1DRef<int> iovrlp, int& ninter, int& nplot3d, FortranArray2DRef<int> inpl3d, int& ip3dsurf, int& nprint, FortranArray2DRef<int> inpr, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> iforce, int& lfgm, int& ncs, FortranArray2DRef<int> icsinfo, int& ihstry, int& ncycmax, FortranArray1DRef<int> iv, FortranArray1DRef<double> time2, FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl, int& intmax, int& nsub1, FortranArray2DRef<int> iindex, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, int& icall, int& iunit11, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& ibufdim, int& nbuf, FortranArray1DRef<int> mglevg, FortranArray1DRef<int> nemgl, FortranArray2DRef<int> ipl3dtmp, int& ntr, FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae, FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae, FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass, FortranArray2DRef<double> damp, FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0, int& nmds, int& maxaes, FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg, FortranArray2DRef<double> xorgae0, FortranArray2DRef<double> yorgae0, FortranArray2DRef<double> zorgae0, FortranArray2DRef<int> icouple, int& iprnsurf)
{
    global_ns::global(myid, maxbl, maxgr, maxseg, maxcs, nplots, mxbli,
        bcvali, bcvalj, bcvalk, nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
        ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek, nblk, nbli,
        limblk, isva, nblon, rkap0g, nblock, levelg, igridg, iflimg, ifdsg,
        iviscg, jdimg, kdimg, idimg, idiagg, nblcg, idegg, jsg, ksg, isg,
        jeg, keg, ieg, mit, ilamlog, ilamhig, jlamlog, jlamhig, klamlog,
        klamhig, iwfg, utrans, vtrans, wtrans, omegax, omegay, omegaz,
        xorig, yorig, zorig, dxmx, dymx, dzmx, dthxmx, dthymx, dthzmx,
        thetax, thetay, thetaz, rfreqt, rfreqr, xorig0, yorig0, zorig0,
        itrans, irotat, idefrm, ngrid, ncgg, nblg, iemg, inewgg, iovrlp,
        ninter, nplot3d, inpl3d, ip3dsurf, nprint, inpr, iadvance, iforce,
        lfgm, ncs, icsinfo, ihstry, ncycmax, iv, time2, thetaxl, thetayl,
        thetazl, intmax, nsub1, iindex, lig, lbg, ibpntsg, iipntsg, icall,
        iunit11, nou, bou, ibufdim, nbuf, mglevg, nemgl, ipl3dtmp, ntr,
        bcfiles, mxbcfil, utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae,
        xorgae, yorgae, zorgae, thtxae, thtyae, thtzae, rfrqtae, rfrqrae,
        icsi, icsf, jcsi, jcsf, kcsi, kcsf, freq, gmass, damp, x0, gf0,
        nmds, maxaes, aesrfdat, perturb, iskip, jskip, kskip, nsegdfrm,
        idfrmseg, iaesurf, maxsegdg, xorgae0, yorgae0, zorgae0, icouple,
        iprnsurf);
}

void global0(int& nplots0, int& maxnode0, int& mxbli0, int& lbcprd0, int& lbcemb0,
             int& lbcrad0, int& maxbl0, int& maxgr0, int& maxseg0, int& maxcs0,
             int& ncycmax0, int& intmax0, int& nsub10, int& intmx0, int& mxxe0,
             int& mptch0, int& msub10, int& ibufdim0, int& nbuf0, int& mxbcfil0,
             int& nmds0, int& maxaes0, int& maxsegdg0, int& ntr, int& nnodes,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& iunit11, int& myid, int& idm0, int& jdm0, int& kdm0)
{
    // COMMON block aliases
    int& iunit5  = cmn_unit5.iunit5;
    int& icgns   = cmn_cgns.icgns;
    int& idef_ss = cmn_elastic_ss.idef_ss;

    // local variables
    float realval[21]; // 1-based indices 1..20
    int titlw[21];     // 1-based indices 1..20

    int nread, nn, n;
    int ncg, ncgmax;
    int idm, jdm, kdm;
    int mdum, mbci0, mbcidim, mbcj0, mbcjdim, mbck0, mbckdim;
    int igrid, idum1, idum2, idum3, idum4, idum5, idum6, idum7, ibctyp, ndat;
    int iglast, nprd, nrad, nsol, nsolgd;
    int mseq, ncyc, ncyctot, ntstep;
    int nbli0, nint0;
    int movie;
    int ncs;
    int ntrans, nrotat;
    int ndefrm, naesrf, naes, ngd, nmodes, iskyhk, iaes, nm, ng;
    int nskip, nskpi, nskpj, nskpk, igr;
    int ncoupl;
    int nint1;
    int iis, iie, jj;
    int iunst, irest, iflagts, ita, ichkd, i2d;
    int nplot3d, nprint, nwrest;
    int ititr;
    int iprint, ichek, ierrflg;
    double dt;
    float xmachw, alphw, reuew, time_val;
    int jt, kt, it_val;
    int nsb10;
    int nfb, int_val, ito, i1, i2, j1, j2;
    char string[81];

    // Initialize nou array
    for (nn = 1; nn <= nbuf0; nn++) {
        nou(nn) = 0;
    }


    // read i/o file names (14 lines)
    nread = 14;
    for (n = 1; n <= nread; n++) {
        echoinp(iunit5, iunit11, 1);
    }

    // read keyword-driven input, if any
    { int _rk_iunit11 = 0, _rk_ierrflg = -99; readkey_ns::readkey(ititr, myid, ibufdim0, nbuf0, bou, nou, _rk_iunit11, _rk_ierrflg); }

    // read title
    if (ititr == 0) {
        echoinp(iunit5, iunit11, 1);
    }

    // read input file up to irest (5 lines)
    nread = 5;
    for (n = 1; n <= nread; n++) {
        echoinp(iunit5, iunit11, 1);
    }

    // read(iunit5,*) realval(1),irest,iflagts,realval(2),iunst,realval(3)
    fortran_read_list(iunit5, &realval[1], &irest, &iflagts, &realval[2], &iunst, &realval[3]);
    fortran_write_unit(iunit11, "%10.5f%10d%10d%10.5f%10d%10.5f\n",
        (double)realval[1], irest, iflagts, (double)realval[2], iunst, (double)realval[3]);
    echoinp(iunit5, iunit11, 0);
    dt = (double)realval[1];
    if ((float)dt < 0.0f) iunst = 0;

    echoinp(iunit5, iunit11, 1);
    // read(iunit5,*) maxgr0,nplot3d,nprint,nwrest,ichkd,i2d,ntstep,ita
    fortran_read_list(iunit5, &maxgr0, &nplot3d, &nprint, &nwrest, &ichkd, &i2d, &ntstep, &ita);
    fortran_write_unit(iunit11, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
        maxgr0, nplot3d, nprint, nwrest, ichkd, i2d, ntstep, ita);
    echoinp(iunit5, iunit11, 0);

    maxgr0 = std::abs(maxgr0);

    echoinp(iunit5, iunit11, 1);

    // read ncg data
    ncgmax = 0;
    for (n = 1; n <= maxgr0; n++) {
        fortran_read_list(iunit5, &ncg, &idum1, &idum2, &idum3, &idum4, &idum5, &idum6);
        fortran_write_unit(iunit11, "%10d%10d%10d%10d%10d%10d%10d\n",
            ncg, idum1, idum2, idum3, idum4, idum5, idum6);
        echoinp(iunit5, iunit11, 0);
        ncgmax = std::max(ncg, ncgmax);
    }

    maxbl0   = (ncgmax + 1) * maxgr0;
    maxnode0 = std::max(maxgr0, nnodes);

    idm0 = 0;
    jdm0 = 0;
    kdm0 = 0;
    // read(iunit5,*)  -- blank line
    {
        FILE* f5 = fortran_get_unit(iunit5);
        if (f5) { char linebuf[256]; fgets(linebuf, sizeof(linebuf), f5); }
    }
    nread = maxgr0;
    for (n = 1; n <= nread; n++) {
        fortran_read_list(iunit5, &idm, &jdm, &kdm);
        if (idm > idm0) idm0 = idm;
        if (jdm > jdm0) jdm0 = jdm;
        if (kdm > kdm0) kdm0 = kdm;
    }
    nread = maxgr0 + 1;

    // read laminar regions section
    for (n = 1; n <= nread; n++) {
        echoinp(iunit5, iunit11, 1);
    }
    // read embedded grid section
    for (n = 1; n <= nread; n++) {
        echoinp(iunit5, iunit11, 1);
    }
    // read idiag/iflim section
    for (n = 1; n <= nread; n++) {
        echoinp(iunit5, iunit11, 1);
    }
    // read ifds/rkap0 section
    for (n = 1; n <= nread; n++) {
        echoinp(iunit5, iunit11, 1);
    }

    // read no. of bc segments section
    maxseg0  = 1;
    int nreadi0 = 0;
    int nreadid = 0;
    int nreadj0 = 0;
    int nreadjd = 0;
    int nreadk0 = 0;
    int nreadkd = 0;

    echoinp(iunit5, iunit11, 1);
    nread = maxgr0;
    for (n = 1; n <= nread; n++) {
        fortran_read_list(iunit5, &mdum, &mbci0, &mbcidim, &mbcj0, &mbcjdim, &mbck0, &mbckdim);
        fortran_write_unit(iunit11, "%10d%10d%10d%10d%10d%10d%10d\n",
            mdum, mbci0, mbcidim, mbcj0, mbcjdim, mbck0, mbckdim);
        echoinp(iunit5, iunit11, 0);
        maxseg0 = std::max({maxseg0, mbci0, mbcidim, mbcj0, mbcjdim, mbck0, mbckdim});
        nreadi0 = nreadi0 + mbci0;
        nreadid = nreadid + mbcidim;
        nreadj0 = nreadj0 + mbcj0;
        nreadjd = nreadjd + mbcjdim;
        nreadk0 = nreadk0 + mbck0;
        nreadkd = nreadkd + mbckdim;
    }


    // read through bc section, counting solid surface, periodic, radial equilibrium bcs
    nprd    = 0;
    nrad    = 0;
    nsol    = 0;
    nsolgd  = 0;
    iglast  = 0;
    lbcrad0 = 1;
    lbcprd0 = 1;

    // i0 face bcs
    echoinp(iunit5, iunit11, 1);
    for (n = 1; n <= nreadi0; n++) {
        fortran_read_list(iunit5, &igrid, &idum2, &ibctyp, &idum4, &idum5, &idum6, &idum7, &ndat);
        fortran_write_unit(iunit11, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
            igrid, idum2, ibctyp, idum4, idum5, idum6, idum7, ndat);
        echoinp(iunit5, iunit11, 0);
        if (igrid != iglast) nsolgd = 0;
        if (std::abs(ibctyp) == 2005) nprd = nprd + 1;
        if (std::abs(ibctyp) == 2006) nrad = nrad + 1;
        if (std::abs(ibctyp) == 2004 || std::abs(ibctyp) == 1005 ||
            std::abs(ibctyp) == 1006 || std::abs(ibctyp) == 2014 ||
            std::abs(ibctyp) == 2024 || std::abs(ibctyp) == 2034 ||
            std::abs(ibctyp) == 2016) {
            nsol   = nsol + 1;
            nsolgd = nsolgd + 1;
        }
        if (std::abs(ndat) > 0) {
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
        }
        iglast = igrid;
    }

    // idim face bcs
    echoinp(iunit5, iunit11, 1);
    for (n = 1; n <= nreadid; n++) {
        fortran_read_list(iunit5, &igrid, &idum2, &ibctyp, &idum4, &idum5, &idum6, &idum7, &ndat);
        fortran_write_unit(iunit11, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
            igrid, idum2, ibctyp, idum4, idum5, idum6, idum7, ndat);
        echoinp(iunit5, iunit11, 0);
        if (igrid != iglast) nsolgd = 0;
        if (std::abs(ibctyp) == 2005) nprd = nprd + 1;
        if (std::abs(ibctyp) == 2006) nrad = nrad + 1;
        if (std::abs(ibctyp) == 2004 || std::abs(ibctyp) == 1005 ||
            std::abs(ibctyp) == 1006 || std::abs(ibctyp) == 2014 ||
            std::abs(ibctyp) == 2024 || std::abs(ibctyp) == 2034 ||
            std::abs(ibctyp) == 2016) {
            nsol   = nsol + 1;
            nsolgd = nsolgd + 1;
        }
        if (std::abs(ndat) > 0) {
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
        }
        iglast = igrid;
    }

    // j0 face bcs
    echoinp(iunit5, iunit11, 1);
    for (n = 1; n <= nreadj0; n++) {
        fortran_read_list(iunit5, &igrid, &idum2, &ibctyp, &idum4, &idum5, &idum6, &idum7, &ndat);
        fortran_write_unit(iunit11, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
            igrid, idum2, ibctyp, idum4, idum5, idum6, idum7, ndat);
        echoinp(iunit5, iunit11, 0);
        if (igrid != iglast) nsolgd = 0;
        if (std::abs(ibctyp) == 2005) nprd = nprd + 1;
        if (std::abs(ibctyp) == 2006) nrad = nrad + 1;
        if (std::abs(ibctyp) == 2004 || std::abs(ibctyp) == 1005 ||
            std::abs(ibctyp) == 1006 || std::abs(ibctyp) == 2014 ||
            std::abs(ibctyp) == 2024 || std::abs(ibctyp) == 2034 ||
            std::abs(ibctyp) == 2016) {
            nsol   = nsol + 1;
            nsolgd = nsolgd + 1;
        }
        if (std::abs(ndat) > 0) {
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
        }
        iglast = igrid;
    }

    // jdim face bcs
    echoinp(iunit5, iunit11, 1);
    for (n = 1; n <= nreadjd; n++) {
        fortran_read_list(iunit5, &igrid, &idum2, &ibctyp, &idum4, &idum5, &idum6, &idum7, &ndat);
        fortran_write_unit(iunit11, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
            igrid, idum2, ibctyp, idum4, idum5, idum6, idum7, ndat);
        echoinp(iunit5, iunit11, 0);
        if (igrid != iglast) nsolgd = 0;
        if (std::abs(ibctyp) == 2005) nprd = nprd + 1;
        if (std::abs(ibctyp) == 2006) nrad = nrad + 1;
        if (std::abs(ibctyp) == 2004 || std::abs(ibctyp) == 1005 ||
            std::abs(ibctyp) == 1006 || std::abs(ibctyp) == 2014 ||
            std::abs(ibctyp) == 2024 || std::abs(ibctyp) == 2034 ||
            std::abs(ibctyp) == 2016) {
            nsol   = nsol + 1;
            nsolgd = nsolgd + 1;
        }
        if (std::abs(ndat) > 0) {
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
        }
        iglast = igrid;
    }


    // k0 face bcs
    echoinp(iunit5, iunit11, 1);
    for (n = 1; n <= nreadk0; n++) {
        fortran_read_list(iunit5, &igrid, &idum2, &ibctyp, &idum4, &idum5, &idum6, &idum7, &ndat);
        fortran_write_unit(iunit11, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
            igrid, idum2, ibctyp, idum4, idum5, idum6, idum7, ndat);
        echoinp(iunit5, iunit11, 0);
        if (igrid != iglast) nsolgd = 0;
        if (std::abs(ibctyp) == 2005) nprd = nprd + 1;
        if (std::abs(ibctyp) == 2006) nrad = nrad + 1;
        if (std::abs(ibctyp) == 2004 || std::abs(ibctyp) == 1005 ||
            std::abs(ibctyp) == 1006 || std::abs(ibctyp) == 2014 ||
            std::abs(ibctyp) == 2024 || std::abs(ibctyp) == 2034 ||
            std::abs(ibctyp) == 2016) {
            nsol   = nsol + 1;
            nsolgd = nsolgd + 1;
        }
        if (std::abs(ndat) > 0) {
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
        }
        iglast = igrid;
    }

    // kdim face bcs
    echoinp(iunit5, iunit11, 1);
    for (n = 1; n <= nreadkd; n++) {
        fortran_read_list(iunit5, &igrid, &idum2, &ibctyp, &idum4, &idum5, &idum6, &idum7, &ndat);
        fortran_write_unit(iunit11, "%10d%10d%10d%10d%10d%10d%10d%10d\n",
            igrid, idum2, ibctyp, idum4, idum5, idum6, idum7, ndat);
        echoinp(iunit5, iunit11, 0);
        if (igrid != iglast) nsolgd = 0;
        if (std::abs(ibctyp) == 2005) nprd = nprd + 1;
        if (std::abs(ibctyp) == 2006) nrad = nrad + 1;
        if (std::abs(ibctyp) == 2004 || std::abs(ibctyp) == 1005 ||
            std::abs(ibctyp) == 1006 || std::abs(ibctyp) == 2014 ||
            std::abs(ibctyp) == 2024 || std::abs(ibctyp) == 2034 ||
            std::abs(ibctyp) == 2016) {
            nsol   = nsol + 1;
            nsolgd = nsolgd + 1;
        }
        if (std::abs(ndat) > 0) {
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
        }
        iglast = igrid;
    }

    lbcprd0 = std::max(lbcprd0, nprd * (ncgmax + 1));
    lbcrad0 = std::max(lbcrad0, nrad * (ncgmax + 1));

    // read mseq data
    echoinp(iunit5, iunit11, 1);
    fortran_read_list(iunit5, &mseq, &idum1, &idum2, &idum3, &idum4);
    fortran_write_unit(iunit11, "%10d%10d%10d%10d%10d\n", mseq, idum1, idum2, idum3, idum4);
    echoinp(iunit5, iunit11, 0);
    echoinp(iunit5, iunit11, 1);
    echoinp(iunit5, iunit11, 1);
    echoinp(iunit5, iunit11, 1);
    ncyctot = 0;
    for (n = 1; n <= mseq; n++) {
        fortran_read_list(iunit5, &ncyc, &idum1, &idum2, &idum3);
        fortran_write_unit(iunit11, "%10d%10d%10d%10d\n", ncyc, idum1, idum2, idum3);
        echoinp(iunit5, iunit11, 0);
        if ((float)dt < 0.0e0f) {
            ncyctot = ncyctot + ncyc;
        } else {
            if (ncyc > 0) {
                ncyctot = ncyctot + ntstep;
            }
        }
    }
    echoinp(iunit5, iunit11, 1);
    for (n = 1; n <= mseq; n++) {
        echoinp(iunit5, iunit11, 1);
    }

    // read through 1-1 interface data
    echoinp(iunit5, iunit11, 1);
    echoinp(iunit5, iunit11, 1);
    fortran_read_list(iunit5, &nbli0);
    fortran_write_unit(iunit11, "%10d\n", nbli0);
    echoinp(iunit5, iunit11, 0);
    nread = 2 * nbli0 + 2;
    for (n = 1; n <= nread; n++) {
        echoinp(iunit5, iunit11, 1);
    }

    mxbli0 = std::max(1, nbli0 * (ncgmax + 1));

    // read through patch data file
    echoinp(iunit5, iunit11, 1);
    echoinp(iunit5, iunit11, 1);
    fortran_read_list(iunit5, &nint0);
    fortran_write_unit(iunit11, "%10d\n", nint0);
    echoinp(iunit5, iunit11, 0);
    if (nint0 == 0) {
        intmax0 = 1;
        nsub10  = 1;
    } else {
        nsub10 = 0;
        // read(22) nint0
        {
            FILE* f22 = fortran_get_unit(22);
            // Fortran: read(22) nint0 — OVERWRITES nint0 with the patch-file value
            // (the deck NINTER was just a flag; the real interface count comes from
            // unit 22). Reading into a throwaway temp left nint0 at the deck value
            // (-1) → intmax0 = nint0*(ncgmax+1) < 0 → bogus array sizing.
            read_unformatted_record(f22, &nint0, sizeof(int));
            echoinp(iunit5, iunit11, 0);
            for (n = 1; n <= std::abs(nint0); n++) {
                // read(22) nsb10
                read_unformatted_record(f22, &nsb10, sizeof(int));
                echoinp(iunit5, iunit11, 0);
                nsub10 = std::max(nsub10, nsb10);
                for (nn = 1; nn <= 12; nn++) {
                    // read(22)
                    skip_unformatted_record(f22);
                }
            }
            rewind(f22);
        }
        intmax0 = nint0 * (ncgmax + 1);
    }


    // plot3d, printout, and control surface data
    nplots0 = 1;
    if (nplot3d > 0) {
        nplots0 = nplot3d;
    } else if (nplot3d < 0) {
        nplots0 = std::max(nsol, std::abs(nplot3d));
        if (std::abs(i2d) > 0) nplots0 = std::max(maxgr0, std::abs(nplot3d));
    }
    if (nprint > 0) {
        nplots0 = std::max(nplots0, nprint);
    } else if (nprint < 0) {
        nplots0 = std::max(nplots0, nsol);
    }

    echoinp(iunit5, iunit11, 1);
    echoinp(iunit5, iunit11, 1);
    for (n = 1; n <= std::abs(nplot3d); n++) {
        echoinp(iunit5, iunit11, 1);
    }

    echoinp(iunit5, iunit11, 1);
    fortran_read_list(iunit5, &movie);
    fortran_write_unit(iunit11, "%10d\n", movie);
    echoinp(iunit5, iunit11, 0);

    echoinp(iunit5, iunit11, 1);
    echoinp(iunit5, iunit11, 1);
    for (n = 1; n <= std::abs(nprint); n++) {
        echoinp(iunit5, iunit11, 1);
    }

    echoinp(iunit5, iunit11, 1);
    echoinp(iunit5, iunit11, 1);
    fortran_read_list(iunit5, &ncs);
    fortran_write_unit(iunit11, "%10d\n", ncs);
    echoinp(iunit5, iunit11, 0);
    echoinp(iunit5, iunit11, 1);
    for (n = 1; n <= ncs; n++) {
        echoinp(iunit5, iunit11, 1);
    }

    maxcs0 = std::max(1, ncs);

    // rigid grid motion data
    if (iunst == 1 || iunst == 3) {
        // translation data
        echoinp(iunit5, iunit11, 1);
        echoinp(iunit5, iunit11, 1);
        fortran_read_list(iunit5, &ntrans);
        fortran_write_unit(iunit11, "%10d\n", ntrans);
        echoinp(iunit5, iunit11, 0);
        echoinp(iunit5, iunit11, 1);
        if (ntrans > 0) {
            echoinp(iunit5, iunit11, 1);
        }
        echoinp(iunit5, iunit11, 1);
        for (n = 1; n <= ntrans; n++) {
            echoinp(iunit5, iunit11, 1);
        }
        echoinp(iunit5, iunit11, 1);
        for (n = 1; n <= ntrans; n++) {
            echoinp(iunit5, iunit11, 1);
        }
        // rotation data
        echoinp(iunit5, iunit11, 1);
        echoinp(iunit5, iunit11, 1);
        fortran_read_list(iunit5, &nrotat);
        fortran_write_unit(iunit11, "%10d\n", nrotat);
        echoinp(iunit5, iunit11, 0);
        echoinp(iunit5, iunit11, 1);
        if (nrotat > 0) {
            echoinp(iunit5, iunit11, 1);
        }
        echoinp(iunit5, iunit11, 1);
        for (n = 1; n <= nrotat; n++) {
            echoinp(iunit5, iunit11, 1);
        }
        echoinp(iunit5, iunit11, 1);
        for (n = 1; n <= nrotat; n++) {
            echoinp(iunit5, iunit11, 1);
        }
    }

    // deforming grid motion data
    nmds0     = 1;
    maxaes0   = 1;
    maxsegdg0 = 1;

    if (iunst > 1 || idef_ss > 0) {
        // deforming mesh data
        echoinp(iunit5, iunit11, 1);
        echoinp(iunit5, iunit11, 1);
        fortran_read_list(iunit5, &ndefrm);
        fortran_write_unit(iunit11, "%10d\n", ndefrm);
        echoinp(iunit5, iunit11, 0);
        if (std::abs(ndefrm) > 0) {
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
            for (int ndef = 1; ndef <= std::abs(ndefrm); ndef++) {
                echoinp(iunit5, iunit11, 1);
            }
            echoinp(iunit5, iunit11, 1);
            for (int ndef = 1; ndef <= std::abs(ndefrm); ndef++) {
                echoinp(iunit5, iunit11, 1);
            }
        } else {
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
        }
        if (ndefrm > 0) {
            maxsegdg0 = nsol + ndefrm;
        } else {
            maxsegdg0 = nsol;
        }


        // aeroelastic mesh data
        echoinp(iunit5, iunit11, 1);
        echoinp(iunit5, iunit11, 1);
        fortran_read_list(iunit5, &naesrf);
        fortran_write_unit(iunit11, "%10d\n", naesrf);
        echoinp(iunit5, iunit11, 0);
        maxaes0 = std::max(maxaes0, naesrf);
        if (naesrf > 0) {
            for (naes = 1; naes <= naesrf; naes++) {
                echoinp(iunit5, iunit11, 1);
                fortran_read_list(iunit5, &iaes, &ngd,
                    &realval[1], &realval[2], &realval[3], &nmodes, &iskyhk);
                fortran_write_unit(iunit11, "%8d%9d%9.4f%9.4f%9.4f%9d%9d\n",
                    iaes, ngd, (double)realval[1], (double)realval[2],
                    (double)realval[3], nmodes, iskyhk);
                echoinp(iunit5, iunit11, 0);
                nmds0 = std::max(nmds0, nmodes);
                echoinp(iunit5, iunit11, 1);
                for (nm = 1; nm <= nmodes; nm++) {
                    echoinp(iunit5, iunit11, 1);
                }
                echoinp(iunit5, iunit11, 1);
                for (nm = 1; nm <= nmodes; nm++) {
                    echoinp(iunit5, iunit11, 1);
                }
                echoinp(iunit5, iunit11, 1);
                for (ng = 1; ng <= std::abs(ngd); ng++) {
                    echoinp(iunit5, iunit11, 1);
                }
                if (ngd >= 0) {
                    maxsegdg0 = maxsegdg0 + ngd;
                }
            }
        } else {
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
            echoinp(iunit5, iunit11, 1);
        }

        // offbody/multiblock mesh deformation data
        echoinp(iunit5, iunit11, 1);
        echoinp(iunit5, iunit11, 1);
        fortran_read_list(iunit5, &nskip, &idum1,
            &realval[1], &realval[2], &realval[3], &realval[4], &idum2);
        fortran_write_unit(iunit11, "%8d%8d%9.6f%9.6f%9.6f%9.6f%9d\n",
            nskip, idum1, (double)realval[1], (double)realval[2],
            (double)realval[3], (double)realval[4], idum2);
        if (std::abs(idum1) == 1) {
            echoinp(iunit5, iunit11, 1);
            if (std::abs(nskip) > 0) {
                for (ng = 1; ng <= std::abs(nskip); ng++) {
                    echoinp(iunit5, iunit11, 1);
                }
            }
        } else {
            echoinp(iunit5, iunit11, 1);
            if (std::abs(nskip) > 0) {
                for (ng = 1; ng <= std::abs(nskip); ng++) {
                    echoinp(iunit5, iunit11, 1);
                    fortran_read_list(iunit5, &igr, &nskpi, &nskpj, &nskpk);
                    echoinp(iunit5, iunit11, 1);
                    iis = -9;
                    iie =  0;
                    for (jj = 1; jj <= 50; jj++) {
                        iis = iis + 10;
                        iie = iie + 10;
                        if (iie > nskpi) iie = nskpi;
                        echoinp(iunit5, iunit11, 1);
                        if (iie == nskpi) goto label325;
                    }
                    label325:;
                    echoinp(iunit5, iunit11, 1);
                    iis = -9;
                    iie =  0;
                    for (jj = 1; jj <= 50; jj++) {
                        iis = iis + 10;
                        iie = iie + 10;
                        if (iie > nskpj) iie = nskpj;
                        echoinp(iunit5, iunit11, 1);
                        if (iie == nskpj) goto label350;
                    }
                    label350:;
                    echoinp(iunit5, iunit11, 1);
                    iis = -9;
                    iie =  0;
                    for (jj = 1; jj <= 50; jj++) {
                        iis = iis + 10;
                        iie = iie + 10;
                        if (iie > nskpk) iie = nskpk;
                        echoinp(iunit5, iunit11, 1);
                        if (iie == nskpk) goto label375;
                    }
                    label375:;
                }
            }
        }

        // multi-motion coupling
        echoinp(iunit5, iunit11, 1);
        echoinp(iunit5, iunit11, 1);
        fortran_read_list(iunit5, &ncoupl);
        fortran_write_unit(iunit11, "%8d\n", ncoupl);
        echoinp(iunit5, iunit11, 0);
        echoinp(iunit5, iunit11, 1);
        if (std::abs(ncoupl) > 0) {
            for (ng = 1; ng <= std::abs(ncoupl); ng++) {
                echoinp(iunit5, iunit11, 1);
            }
        }

        maxsegdg0 = std::max(1, maxsegdg0);

    } // end if (iunst > 1 || idef_ss > 0)


    // dynamic patch data
    intmx0 = 1;
    msub10 = 1;
    nint1  = 0;

    if (iunst > 0) {
        bool at_eof = false;
        FILE* f5 = fortran_get_unit(iunit5);
        // read(iunit5,'(a80)',end=888) string
        if (f5 && !at_eof) {
            if (fgets(string, 81, f5) == nullptr) {
                at_eof = true;
            } else {
                int slen = (int)strlen(string);
                // a80 record: if buffer filled without the newline, drain rest of line
                if (slen > 0 && string[slen-1] != '\n') { int c; while ((c=fgetc(f5))!=EOF && c!='\n') {} }
                if (slen > 0 && string[slen-1] == '\n') { string[slen-1] = '\0'; slen--; }
                fortran_write_unit(iunit11, "%-80s\n", string);
                echoinp(iunit5, iunit11, 0);
            }
        }
        // read(iunit5,'(a80)',end=888) string
        if (f5 && !at_eof) {
            if (fgets(string, 81, f5) == nullptr) {
                at_eof = true;
            } else {
                int slen = (int)strlen(string);
                if (slen > 0 && string[slen-1] != '\n') { int c; while ((c=fgetc(f5))!=EOF && c!='\n') {} }
                if (slen > 0 && string[slen-1] == '\n') { string[slen-1] = '\0'; slen--; }
                fortran_write_unit(iunit11, "%-80s\n", string);
                echoinp(iunit5, iunit11, 0);
            }
        }
        // read(iunit5,*,end=888) nint1
        if (f5 && !at_eof) {
            char linebuf[256];
            if (fgets(linebuf, sizeof(linebuf), f5) == nullptr) {
                at_eof = true;
            } else {
                sscanf(linebuf, "%d", &nint1);
                fortran_write_unit(iunit11, "%10d\n", nint1);
                echoinp(iunit5, iunit11, 0);
            }
        }
        // label 888: at_eof == true means we jumped here

        if (nint1 > 0) {
            intmx0 = nint1 * (ncgmax + 1);

            echoinp(iunit5, iunit11, 1);
            for (n = 1; n <= nint1; n++) {
                echoinp(iunit5, iunit11, 1);
            }
            for (n = 1; n <= nint1; n++) {
                echoinp(iunit5, iunit11, 1);
                fortran_read_list(iunit5, &int_val, &ito, &i1, &i2, &j1, &j2, &nfb);
                fortran_write_unit(iunit11, "%10d%10d%10d%10d%10d%10d%10d\n",
                    int_val, ito, i1, i2, j1, j2, nfb);
                echoinp(iunit5, iunit11, 0);
                for (nn = 1; nn <= nfb; nn++) {
                    echoinp(iunit5, iunit11, 1);
                    echoinp(iunit5, iunit11, 1);
                    echoinp(iunit5, iunit11, 1);
                    echoinp(iunit5, iunit11, 1);
                }
                msub10 = std::max(msub10, nfb);
            }
        }

    } // end if (iunst > 0)

    nsub10 = std::max(nsub10, msub10);

    // if a restart, read restart file to determine how many cycles done
    ntr = 0;
    if (irest != 0) {
        if (icgns != 1) {
            FILE* f2 = fortran_get_unit(2);
            if (f2) {
                // read(2,end=999) titlw,xmachw,jt,kt,it,alphw,reuew,ntr,time
                // This is an unformatted record containing:
                // 20 ints (titlw), float xmachw, int jt, int kt, int it,
                // float alphw, float reuew, int ntr, float time
                struct { int titlw[20]; float xmachw; int jt; int kt; int it_v;
                         float alphw; float reuew; int ntr_v; float time_v; } rec;
                int nr = read_unformatted_record(f2, &rec, (int)sizeof(rec));
                if (nr > 0) {
                    ntr = rec.ntr_v;
                }
                // label 999: continue
                rewind(f2);
                fortran_close_unit(2);
            }
        }
    }
    ncycmax0 = std::max(ncyctot + ntr, 1);

    // set remaining parameters needed by sizer (precfl3d)
    mxxe0   = 1;
    mptch0  = 1;
    lbcemb0 = 6 * maxbl0;

    // print out current parameters if desired or if one is zero (error)
    iprint = 0;
    ichek  = 1;

    if (nplots0   <= 0 ||
        maxnode0  <= 0 ||
        mxbli0    <= 0 ||
        lbcprd0   <= 0 ||
        lbcemb0   <= 0 ||
        lbcrad0   <= 0 ||
        maxbl0    <= 0 ||
        maxgr0    <= 0 ||
        maxseg0   <= 0 ||
        maxcs0    <= 0 ||
        ncycmax0  <= 0 ||
        intmax0   <= 0 ||
        nsub10    <= 0 ||
        intmx0    <= 0 ||
        mxxe0     <= 0 ||
        mptch0    <= 0 ||
        msub10    <= 0 ||
        ibufdim0  <= 0 ||
        nbuf0     <= 0 ||
        maxsegdg0 <= 0 ||
        mxbcfil0  <= 0) ichek = 0;

    if (ichek == 0) iprint = 1;

    if (iprint > 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "nplots0   = %d", nplots0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "maxnode0  = %d", maxnode0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "mxbli0    = %d", mxbli0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "lbcprd0   = %d", lbcprd0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "lbcemb0   = %d", lbcemb0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "lbcrad0   = %d", lbcrad0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "maxbl0    = %d", maxbl0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "maxgr0    = %d", maxgr0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "maxseg0   = %d", maxseg0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "maxcs0    = %d", maxcs0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "ncycmax0  = %d", ncycmax0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "intmax0   = %d", intmax0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "nsub10    = %d", nsub10);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "intmx0    = %d", intmx0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "mxxe0     = %d", mxxe0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "mptch0    = %d", mptch0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "msub10    = %d", msub10);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "ibufdim0  = %d", ibufdim0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "nbuf0     = %d", nbuf0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "mxbcfil0  = %d", mxbcfil0);
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120, "maxsegdg0 = %d", maxsegdg0);
    }

    if (ichek == 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim0);
        std::snprintf(bou(nou(1), 1), 120,
            "error in routine global0 - one of the parameters listed above is zero");
        ierrflg = -99;
        termn8_ns::termn8(myid, ierrflg, ibufdim0, nbuf0, bou, nou);
    }

    {
        FILE* f11 = fortran_get_unit(iunit11);
        if (f11) rewind(f11);
    }

    return;
}

static void echoinp(int& iunit5, int& iunit11, int irw)
{
    // Fortran: read(iunit5,'(a80)') string  — reads one full RECORD, advancing
    // past the record terminator regardless of line length. A too-small fgets
    // buffer (e.g. 81 for an 80-char + newline record) would read only 80 chars
    // and leave the newline, desyncing every subsequent read. Consume the full
    // line (through the newline), using a buffer large enough for the record.
    char string[512];
    if (irw > 0) {
        FILE* f5 = fortran_get_unit(iunit5);
        if (f5) {
            if (fgets(string, (int)sizeof(string), f5) != nullptr) {
                int slen = (int)strlen(string);
                // If the record is longer than the buffer, drain the rest of the
                // line so the file is positioned at the start of the next record.
                while (slen > 0 && string[slen-1] != '\n' && !feof(f5)) {
                    int c;
                    while ((c = fgetc(f5)) != EOF && c != '\n') { }
                    break;
                }
                if (slen > 0 && string[slen-1] == '\n') { string[slen-1] = '\0'; slen--; }
                fortran_write_unit(iunit11, "%-80s\n", string);
            }
        }
    }
    my_flush_ns::my_flush(iunit11);
}

} // namespace global0_ns
