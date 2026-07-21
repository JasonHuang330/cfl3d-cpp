// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "mgbl.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include "module_kwstm.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

#include "cputim.h"
#include "termn8.h"
#include "global.h"
#include "pointers.h"
#include "writ_buf.h"
#include "setup.h"
#include "qinter.h"
#include "mgblk.h"
#include "qout.h"
#include "csout.h"
#include "forceout.h"
#include "yplusout.h"
#include "histout.h"
#include "umalloc.h"

// Parameterized cputim overload (not in cputim.h but defined in cputim.cpp)
namespace cputim_ns { void cputim(int& mode, int& nnodes, char* string, int& myhost, int& myid, int& mycomm, int& iunit); }
namespace mgbl_ns {

void mgbl(int& mwork, int& mworki, int& nplots, int& minnode, int& iitot, int& intmax, int& maxxe, int& mxbli, int& nsub1, int& lbcprd, int& lbcemb, int& lbcrad, int& maxbl, int& maxgr, int& maxseg, int& maxcs, int& ncycmax, int& intmx, int& mxxe, int& mptch, int& msub1, int& ibufdim, int& nbuf, int& mxbcfil, int& istat_size, FortranArray1DRef<double> work, FortranArray1DRef<int> iwork, FortranArray3DRef<int> lwdat, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, FortranArray2DRef<int> istat2_bl, FortranArray2DRef<int> istat2_pa, FortranArray2DRef<int> istat2_em, FortranArray2DRef<int> istat2_pe, FortranArray2DRef<int> nblk, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, FortranArray1DRef<double> resmx, FortranArray1DRef<int> imx, FortranArray1DRef<int> jmx, FortranArray1DRef<int> kmx, FortranArray1DRef<double> vormax, FortranArray1DRef<int> ivmax, FortranArray1DRef<int> jvmax, FortranArray1DRef<int> kvmax, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray1DRef<int> iovrlp, FortranArray3DRef<double> qb, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg, FortranArray1DRef<double> dxintg, FortranArray1DRef<double> dyintg, FortranArray1DRef<double> dzintg, FortranArray1DRef<int> iiig, FortranArray1DRef<int> jjig, FortranArray1DRef<int> kkig, FortranArray2DRef<double> rkap0g, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<int> idiagg, FortranArray1DRef<int> nblcg, FortranArray2DRef<int> idegg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray2DRef<int> mit, FortranArray1DRef<int> ilamlog, FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamlog, FortranArray1DRef<int> klamhig, FortranArray2DRef<int> iwfg, FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx, FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz, FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr, FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0, FortranArray1DRef<double> time2, FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg, FortranArray2DRef<int> inpl3d, FortranArray2DRef<int> inpr, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> iforce, FortranArray1DRef<double> rms, FortranArray1DRef<double> clw, FortranArray1DRef<double> cdw, FortranArray1DRef<double> cdpw, FortranArray1DRef<double> cdvw, FortranArray1DRef<double> cxw, FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw, FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw, FortranArray1DRef<double> cmzw, int& n_clcd, FortranArray3DRef<double> clcd, int& nblocks_clcd, FortranArray2DRef<int> blocks_clcd, FortranArray1DRef<double> chdw, FortranArray1DRef<double> swetw, FortranArray1DRef<double> fmdotw, FortranArray1DRef<double> cfttotw, FortranArray1DRef<double> cftmomw, FortranArray1DRef<double> cftpw, FortranArray1DRef<double> cftvw, FortranArray1DRef<double> swett, FortranArray1DRef<double> clt, FortranArray1DRef<double> cdt, FortranArray1DRef<double> cxt, FortranArray1DRef<double> cyt, FortranArray1DRef<double> czt, FortranArray1DRef<double> cmxt, FortranArray1DRef<double> cmyt, FortranArray1DRef<double> cmzt, FortranArray1DRef<double> cdpt, FortranArray1DRef<double> cdvt, FortranArray1DRef<double> sx, FortranArray1DRef<double> sy, FortranArray1DRef<double> sz, FortranArray1DRef<double> stot, FortranArray1DRef<double> pav, FortranArray1DRef<double> ptav, FortranArray1DRef<double> tav, FortranArray1DRef<double> ttav, FortranArray1DRef<double> xmav, FortranArray1DRef<double> fmdot, FortranArray1DRef<double> cfxp, FortranArray1DRef<double> cfyp, FortranArray1DRef<double> cfzp, FortranArray1DRef<double> cfdp, FortranArray1DRef<double> cflp, FortranArray1DRef<double> cftp, FortranArray1DRef<double> cfxv, FortranArray1DRef<double> cfyv, FortranArray1DRef<double> cfzv, FortranArray1DRef<double> cfdv, FortranArray1DRef<double> cflv, FortranArray1DRef<double> cftv, FortranArray1DRef<double> cfxmom, FortranArray1DRef<double> cfymom, FortranArray1DRef<double> cfzmom, FortranArray1DRef<double> cfdmom, FortranArray1DRef<double> cflmom, FortranArray1DRef<double> cftmom, FortranArray1DRef<double> cfxtot, FortranArray1DRef<double> cfytot, FortranArray1DRef<double> cfztot, FortranArray1DRef<double> cfdtot, FortranArray1DRef<double> cfltot, FortranArray1DRef<double> cfttot, FortranArray2DRef<int> icsinfo, FortranArray2DRef<double> windex, FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt, FortranArray2DRef<double> windx, FortranArray2DRef<int> iindx, FortranArray1DRef<int> llimit, FortranArray1DRef<int> iitmax, FortranArray1DRef<int> mmcxie, FortranArray1DRef<int> mmceta, FortranArray1DRef<int> ncheck, FortranArray1DRef<int> iifit, FortranArray1DRef<int> mblkpt, FortranArray1DRef<int> iic0, FortranArray1DRef<int> iiorph, FortranArray1DRef<int> iitoss, FortranArray1DRef<int> ifiner, FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz, FortranArray2DRef<double> dx, FortranArray2DRef<double> dy, FortranArray2DRef<double> dz, FortranArray2DRef<double> dthetx, FortranArray2DRef<double> dthety, FortranArray2DRef<double> dthetz, FortranArray1DRef<int> lout, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, FortranArray1DRef<int> iiint1, FortranArray1DRef<int> iiint2, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, FortranArray1DRef<int> jte, FortranArray1DRef<int> kte, FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm, FortranArray1DRef<int> nblk1, FortranArray1DRef<int> nblk2, FortranArray3DRef<double> xte, FortranArray3DRef<double> yte, FortranArray3DRef<double> zte, FortranArray3DRef<double> xmi, FortranArray3DRef<double> ymi, FortranArray3DRef<double> zmi, FortranArray3DRef<double> xmie, FortranArray3DRef<double> ymie, FortranArray3DRef<double> zmie, FortranArray3DRef<double> sxie, FortranArray3DRef<double> seta, FortranArray2DRef<double> sxie2, FortranArray2DRef<double> seta2, FortranArray2DRef<double> xie2s, FortranArray2DRef<double> eta2s, FortranArray1DRef<double> temp, FortranArray2DRef<double> x2, FortranArray2DRef<double> y2, FortranArray2DRef<double> z2, FortranArray2DRef<double> x1, FortranArray2DRef<double> y1, FortranArray2DRef<double> z1, FortranArray2DRef<double> factjlo, FortranArray2DRef<double> factjhi, FortranArray2DRef<double> factklo, FortranArray2DRef<double> factkhi, FortranArray1DRef<int> ifrom, FortranArray1DRef<double> geom_miss, FortranArray1DRef<double> period_miss, FortranArray2DRef<int> isav_blk, FortranArray2DRef<int> isav_prd, FortranArray2DRef<int> isav_pat, FortranArray3DRef<int> isav_pat_b, FortranArray2DRef<int> isav_dpat, FortranArray3DRef<int> isav_dpat_b, FortranArray2DRef<int> isav_emb, FortranArray1DRef<int> mblk2nd, int& ntr, FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae, FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae, FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass, FortranArray2DRef<double> damp, FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0, int& nmds, int& maxaes, FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb, int& memuse, FortranArray1DRef<char[80]> bcfiles, FortranArray3DRef<int> islavept, int& nslave, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray3DRef<double> bmat, FortranArray3DRef<double> stm, FortranArray3DRef<double> stmi, FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn, FortranArray2DRef<double> gforcn, FortranArray2DRef<double> gforcnm, FortranArray2DRef<double> gforcs, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg, int& nmaster, FortranArray4DRef<double> aehist, FortranArray1DRef<double> timekeep, FortranArray2DRef<double> xorgae0, FortranArray2DRef<double> yorgae0, FortranArray2DRef<double> zorgae0, FortranArray2DRef<int> icouple, FortranArray2DRef<int> nblelst, FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax, FortranArray1DRef<int> kskmax, FortranArray1DRef<double> ue, FortranArray1DRef<int> irdrea)
{

    // COMMON block aliases
    int32_t& ndefrm     = cmn_elastic.ndefrm;
    int32_t& naesrf     = cmn_elastic.naesrf;
    int32_t& ialphit    = cmn_alphait.ialphit;
    float&   cltarg     = cmn_alphait.cltarg;
    float&   rlxalph    = cmn_alphait.rlxalph;
    float&   dalim      = cmn_alphait.dalim;
    float&   dalpha     = cmn_alphait.dalpha;
    int32_t& icycupdt   = cmn_alphait.icycupdt;
    float*   rate       = cmn_cpurate.rate;    // rate[5], 0-based
    float*   ratesub    = cmn_cpurate.ratesub; // ratesub[5], 0-based
    int32_t* ncell      = cmn_cpurate.ncell;   // ncell[20], 0-based
    float&   dt0        = cmn_cfl.dt0;
    float&   dtold      = cmn_cfl.dtold;
    float&   gamma      = cmn_fluid.gamma;
    float&   gm1        = cmn_fluid.gm1;
    float&   gp1        = cmn_fluid.gp1;
    float&   gm1g       = cmn_fluid.gm1g;
    float&   gp1g       = cmn_fluid.gp1g;
    float&   ggm1       = cmn_fluid.ggm1;
    float&   xmach      = cmn_info.xmach;
    float&   alpha_info = cmn_info.alpha;
    float&   beta_info  = cmn_info.beta;
    float&   dt         = cmn_info.dt;
    float&   fmax       = cmn_info.fmax;
    int32_t& nit        = cmn_info.nit;
    int32_t& ntt        = cmn_info.ntt;
    int32_t& nitfo      = cmn_info.nitfo;
    int32_t& iflagts    = cmn_info.iflagts;
    int32_t& nres       = cmn_info.nres;
    int32_t& mgflag     = cmn_info.mgflag;
    int32_t& iconsf     = cmn_info.iconsf;
    int32_t& mseq       = cmn_info.mseq;
    int32_t& ngam       = cmn_info.ngam;
    int32_t& iipv       = cmn_info.iipv;
    int32_t* ncyc1      = cmn_info.ncyc1;   // ncyc1[5], 0-based
    int32_t* levelt     = cmn_info.levelt;  // levelt[5], 0-based
    int32_t* nitfo1     = cmn_info.nitfo1;  // nitfo1[5], 0-based
    int32_t* nsm        = cmn_info.nsm;     // nsm[5], 0-based
    int32_t* levelb     = cmn_info.levelb;  // levelb[5], 0-based
    int32_t* idiag      = cmn_info.idiag;   // idiag[3], 0-based
    int32_t* iflim      = cmn_info.iflim;   // iflim[3], 0-based
    int32_t& ivmx       = cmn_maxiv.ivmx;
    int32_t& isklton    = cmn_sklton.isklton;
    float&   time_unst  = cmn_unst.time;
    float&   cfltau     = cmn_unst.cfltau;
    int32_t& ntstep     = cmn_unst.ntstep;
    int32_t& ita        = cmn_unst.ita;
    int32_t& iunst      = cmn_unst.iunst;
    float&   cfltau0    = cmn_unst.cfltau0;
    float&   cfltauMax  = cmn_unst.cfltaumax;
    int32_t& movie      = cmn_moov.movie;
    int32_t& nframes    = cmn_moov.nframes;
    int32_t& icall1     = cmn_moov.icall1;
    int32_t& lhdr       = cmn_moov.lhdr;
    int32_t& icoarsemovie = cmn_moov.icoarsemovie;
    int32_t& i2dmovie   = cmn_moov.i2dmovie;
    char*    grid       = cmn_filenam.grid;
    char*    plt3dg     = cmn_filenam.plt3dg;
    char*    plt3dq     = cmn_filenam.plt3dq;
    char*    output_fn  = cmn_filenam.output;
    char*    residual   = cmn_filenam.residual;
    char*    turbres    = cmn_filenam.turbres;
    char*    blomx      = cmn_filenam.blomx;
    char*    output2    = cmn_filenam.output2;
    char*    printout   = cmn_filenam.printout;
    char*    pplunge    = cmn_filenam.pplunge;
    char*    ovrlap     = cmn_filenam.ovrlap;
    char*    patch      = cmn_filenam.patch;
    char*    restrt     = cmn_filenam.restrt;
    char*    subres     = cmn_filenam.subres;
    char*    subtur     = cmn_filenam.subtur;
    char*    grdmov     = cmn_filenam.grdmov;
    char*    alphahist  = cmn_filenam.alphahist;
    char*    errfile    = cmn_filenam.errfile;
    char*    preout     = cmn_filenam.preout;
    char*    aeinp      = cmn_filenam.aeinp;
    char*    aeout      = cmn_filenam.aeout;
    char*    sdhist     = cmn_filenam.sdhist;
    char*    avgg       = cmn_filenam.avgg;
    char*    avgq       = cmn_filenam.avgq;
    int32_t& iunit5     = cmn_unit5.iunit5;
    float&   radtodeg   = cmn_conversion.radtodeg;
    int32_t& nnodes     = cmn_mydist2.nnodes;
    int32_t& myhost     = cmn_mydist2.myhost;
    int32_t& myid       = cmn_mydist2.myid;
    int32_t& mycomm     = cmn_mydist2.mycomm;
    int32_t& numprocs   = cmn_proces.numprocs;
    // tim(3,3) stored column-major in tim[9]: tim(i,j) = tim[(i-1) + (j-1)*3]
    float*   tim_arr    = cmn_time1.tim;
    float*   tm         = cmn_time1.tm;
    int32_t& ibin       = cmn_bin.ibin;
    int32_t& iblnk      = cmn_bin.iblnk;
    int32_t& iblnkfr    = cmn_bin.iblnkfr;
    int32_t& ip3dgrad   = cmn_bin.ip3dgrad;
    float&   beta1      = cmn_deformz.beta1;
    float&   beta2      = cmn_deformz.beta2;
    float&   alpha1     = cmn_deformz.alpha1;
    float&   alpha2     = cmn_deformz.alpha2;
    int32_t& isktyp     = cmn_deformz.isktyp;
    int32_t& negvol     = cmn_deformz.negvol;
    int32_t& meshdef    = cmn_deformz.meshdef;
    int32_t& nsprgit    = cmn_deformz.nsprgit;
    int32_t& ndgrd      = cmn_deformz.ndgrd;
    int32_t& ndwrt      = cmn_deformz.ndwrt;

    // Local variables
    char     string[51];
    float    walltime, totaltime;
    int      ndays;
    float    pi;
    int      iwk1, iwk2, iwk3, iwk4, iwk5, mworki1;
    int      nwork1;
    int      maxl, lembed, nstart, nwork;
    int      iunit, icallgl, icallpt, mpihost;
    int      ll, lou, kou, nf, iii, n, nbl, m;
    int      iseq, lglobal, ncells, nttuse;
    int      iflag;
    int      nummem;
    int      stats = 0;
    double   epsc0, epsrot;
    int      nblock, nbli, ngrid, ninter, nplot3d, ip3dsurf, nprint;
    int      lfgm, ncs, ihstry;
    int      nintr;
    int      iprnsurf;

    // Allocatable arrays (translated from Fortran allocatable)
    FortranArray2D<double> rmstr_arr;
    FortranArray2D<int>    nneg_arr;


    // ***********************************************************************
    // initialization
    // ***********************************************************************
    //
    // initialize timing utility
    //
    std::memset(string, ' ', 50);
    string[50] = '\0';
    {
        int _mode0 = 0, _nn1 = 1, _iunit11 = 11;
        cputim_ns::cputim(_mode0, _nn1, string, myhost, myid, mycomm, _iunit11);
    }

    //
    // initialize data for grids stored in array w
    //
    for (ll = 1; ll <= mwork; ll++) {
        work(ll) = 0.0;
    }
    //
    for (ll = 1; ll <= mworki; ll++) {
        iwork(ll) = 0;
    }
    //
    // initialize internal write buffer
    //
    for (lou = 1; lou <= nbuf; lou++) {
        for (kou = 1; kou <= ibufdim; kou++) {
            std::memset(bou(kou, lou), ' ', 120);
        }
        nou(lou) = 0;
    }
    //
    // initialize auxiliary file names for 2000 series bc's
    //
    for (nf = 1; nf <= mxbcfil; nf++) {
        std::memset(bcfiles(nf), ' ', 80);
        std::memcpy(bcfiles(nf), "null", 4);
    }
    //
    // initialize radtodeg
    //
    pi       = 4.0f * std::atan(1.0f);
    radtodeg = 180.0e0f / pi;
    //
    if (myid == myhost) {
        //
        // ***********************************************************************
        //    write version information, etc. to main output file
        // ***********************************************************************
        //
        // output banner
        //
        // format 83: 35(2h *) = "* " repeated 35 times = 70 chars
        // format 87: 2(2h *),62x,2(2h *) = "* * " + 62 spaces + "* *"
        // write(11,83) twice
        fortran_write_unit(11, "%s\n", "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
        fortran_write_unit(11, "%s\n", "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
        // write(11,87)
        fortran_write_unit(11, "%s\n", "*                                                                               *");
        // write(11,9900)
        fortran_write_unit(11, "%s\n", "* *   CFL3D :  CHARACTERISTIC-BASED SCHEME FOR STEADY/UNSTEADY   * *");
        fortran_write_unit(11, "%s\n", "* *            SOLUTIONS TO THE EULER/NAVIER-STOKES EQUATIONS        * *");
        // write(11,87)
        fortran_write_unit(11, "%s\n", "*                                                                               *");
        // write(11,9990)
        fortran_write_unit(11, "%s\n", "* *   VERSION 6.7 :  Computational Fluids Lab,  Mail Stop 128,    * *");
        fortran_write_unit(11, "%s\n", "* *                  NASA Langley Research Center, Hampton, VA   * *");
        fortran_write_unit(11, "%s\n", "* *                  Release Date:  February  1, 2017.           * *");
        // write(11,87)
        fortran_write_unit(11, "%s\n", "*                                                                               *");
        // write(11,83) twice
        fortran_write_unit(11, "%s\n", "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
        fortran_write_unit(11, "%s\n", "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");

        fortran_write_unit(11, "%s\n", "CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged");
        fortran_write_unit(11, "%s\n", "Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones");
        fortran_write_unit(11, "%s\n", "with point-matched, patched, overset, or embedded connectivities. Both");
        fortran_write_unit(11, "%s\n", "multigrid and mesh sequencing are available in time-accurate or");
        fortran_write_unit(11, "%s\n", "steady-state modes.");
        fortran_write_unit(11, "%s\n", " ");
        fortran_write_unit(11, "%s\n", "-----------------------------------------------------------------------");
        fortran_write_unit(11, "%s\n", "Copyright 2001 United States Government as represented by the Administrator");
        fortran_write_unit(11, "%s\n", "of the National Aeronautics and Space Administration. All Rights Reserved.");
        fortran_write_unit(11, "%s\n", " ");
        fortran_write_unit(11, "%s\n", "The CFL3D platform is licensed under the Apache License, Version 2.0");
        fortran_write_unit(11, "%s\n", "(the \"License\"); you may not use this file except in compliance with the");
        fortran_write_unit(11, "%s\n", "License. You may obtain a copy of the License at");
        fortran_write_unit(11, "%s\n", "http://www.apache.org/licenses/LICENSE-2.0.");
        fortran_write_unit(11, "%s\n", " ");
        fortran_write_unit(11, "%s\n", "Unless required by applicable law or agreed to in writing, software");
        fortran_write_unit(11, "%s\n", "distributed under the License is distributed on an \"AS IS\" BASIS, WITHOUT");
        fortran_write_unit(11, "%s\n", "WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the");
        fortran_write_unit(11, "%s\n", "License for the specific language governing permissions and limitations");
        fortran_write_unit(11, "%s\n", "under the License.");
        fortran_write_unit(11, "%s\n", "-----------------------------------------------------------------------");
        //
        // write(11,13) real(float(memuse))/1.e6
        // format 13: /,' memory allocation: ',f12.6,' Mbytes',' single precision'
        fortran_write_unit(11, "\n memory allocation: %12.6f Mbytes single precision\n", (float)memuse / 1.0e6f);
        //
        // write(11,88) - format 88: /19hinput/output files:
        fortran_write_unit(11, "\ninput/output files:\n");
        //
        // write(11,'(''  '',a60)') for each file name
        fortran_write_unit(11, "  %-60.60s\n", grid);
        fortran_write_unit(11, "  %-60.60s\n", plt3dg);
        fortran_write_unit(11, "  %-60.60s\n", plt3dq);
        fortran_write_unit(11, "  %-60.60s\n", output_fn);
        fortran_write_unit(11, "  %-60.60s\n", residual);
        fortran_write_unit(11, "  %-60.60s\n", turbres);
        fortran_write_unit(11, "  %-60.60s\n", blomx);
        fortran_write_unit(11, "  %-60.60s\n", output2);
        fortran_write_unit(11, "  %-60.60s\n", printout);
        fortran_write_unit(11, "  %-60.60s\n", pplunge);
        fortran_write_unit(11, "  %-60.60s\n", ovrlap);
        fortran_write_unit(11, "  %-60.60s\n", patch);
        fortran_write_unit(11, "  %-60.60s\n", restrt);
        //
        // ***********************************************************************
        //    read in global information
        // ***********************************************************************
        //
        iwk1    = 1;
        iwk2    = iwk1 + maxbl;
        iwk3    = iwk2 + maxbl;
        iwk4    = iwk3 + maxbl;
        iwk5    = iwk4 + 11 * nplots;
        mworki1 = mworki - iwk5 + 1;
        if (mworki1 <= 0) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120, "stopping...not enough integer work space for subroutine global");
            { int _m1 = -1; termn8_ns::termn8(myid, _m1, ibufdim, nbuf, bou, nou); }
        }
        //
        // read input file up to title line
        {
            int nread = 14;
            for (n = 1; n <= nread; n++) {
                // read(iunit5,*)  -- skip a line
                FILE* f5 = fortran_get_unit(iunit5);
                char linebuf[256];
                if (f5) fgets(linebuf, sizeof(linebuf), f5);
            }
        }
        //
        iunit   = 11;
        icallgl = 1;
        {
            // iwork(iwk1) slice starting at iwk1
            FortranArray1DRef<int> iwork_iwk1(iwork.data() + (iwk1 - 1), iwork.size(1) - iwk1 + 1);
            FortranArray1DRef<int> iwork_iwk2(iwork.data() + (iwk2 - 1), iwork.size(1) - iwk2 + 1);
            FortranArray1DRef<int> iwork_iwk3(iwork.data() + (iwk3 - 1), iwork.size(1) - iwk3 + 1);
            FortranArray2DRef<int> iwork_iwk4(&iwork(iwk4), 11, nplots);
            global_ns::global(myid, maxbl, maxgr, maxseg, maxcs, nplots, mxbli,
                              bcvali, bcvalj, bcvalk, nbci0, nbcj0, nbck0,
                              nbcidim, nbcjdim, nbckdim, ibcinfo, jbcinfo,
                              kbcinfo, bcfilei, bcfilej, bcfilek, nblk, nbli,
                              limblk, isva, nblon, rkap0g, nblock, levelg,
                              igridg, iflimg, ifdsg, iviscg, jdimg, kdimg,
                              idimg, idiagg, nblcg, idegg, jsg, ksg, isg, jeg,
                              keg, ieg, mit, ilamlog, ilamhig, jlamlog,
                              jlamhig, klamlog, klamhig, iwfg, utrans, vtrans,
                              wtrans, omegax, omegay, omegaz, xorig, yorig,
                              zorig, dxmx, dymx, dzmx, dthxmx, dthymx,
                              dthzmx, thetax, thetay, thetaz, rfreqt, rfreqr,
                              xorig0, yorig0, zorig0, itrans, irotat, idefrm,
                              ngrid, ncgg, nblg, iemg, inewgg, iovrlp, ninter,
                              nplot3d, inpl3d, ip3dsurf, nprint, inpr,
                              iadvance, iforce, lfgm, ncs, icsinfo, ihstry,
                              ncycmax, iwork_iwk1, time2, thetaxl, thetayl, thetazl,
                              intmax, nsub1, iindex, lig, lbg, ibpntsg,
                              iipntsg, icallgl, iunit, nou, bou, ibufdim, nbuf,
                              iwork_iwk2, iwork_iwk3, iwork_iwk4, ntr, bcfiles,
                              mxbcfil, utrnsae, vtrnsae, wtrnsae, omgxae,
                              omgyae, omgzae, xorgae, yorgae, zorgae, thtxae,
                              thtyae, thtzae, rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf,
                              kcsi, kcsf, freq, gmass, damp, x0, gf0, nmds, maxaes,
                              aesrfdat, perturb, iskip, jskip, kskip, nsegdfrm,
                              idfrmseg, iaesurf, maxsegdg, xorgae0, yorgae0, zorgae0,
                              icouple, iprnsurf);
        }
        for (iii = 1; iii <= iwk5; iii++) {
            iwork(iii) = 0;
        }
        //
        // ***********************************************************************
        //    set up block-to-node mapping array
        // ***********************************************************************
        //
        for (int ibloc = 1; ibloc <= maxbl; ibloc++) {
            mblk2nd(ibloc) = 0;
        }
        //
    } // end if (myid == myhost)

    //
    // ***********************************************************************
    //  compute starting locations for arrays and print out block to
    //  node mapping for mpi case
    // ***********************************************************************
    //
    if (mworki < 2 * maxbl) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, "stopping...not enough integer work space for subroutine pointers");
        { int _m1 = -1; termn8_ns::termn8(myid, _m1, ibufdim, nbuf, bou, nou); }
    }
    icallpt = 1;
    mpihost = 0;
    {
        FortranArray1DRef<int> iwork_1(iwork.data(), iwork.size(1));
        FortranArray1DRef<int> iwork_maxbl1(iwork.data() + maxbl, iwork.size(1) - maxbl);
        pointers_ns::pointers(lw, lw2, maxl, lembed, nstart, nwork, mwork, maxbl,
                              maxgr, maxseg, lwdat, levelg, igridg, iviscg,
                              idimg, jdimg, kdimg, nblcg, itrans, irotat, idefrm,
                              nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
                              ibcinfo, jbcinfo, kbcinfo, ngrid, ncgg, nblg,
                              iemg, nblock, myhost, myid, mblk2nd, nou, bou, nbuf,
                              ibufdim, iwork_1, ilamlog, jlamlog,
                              klamlog, ilamhig, jlamhig, klamhig, idegg, iwfg,
                              idiagg, iflimg, ifdsg, rkap0g, jsg, ksg, isg, jeg,
                              keg, ieg, iwork_maxbl1, icallpt, nmds, maxaes, mpihost);
    }
    //
    for (iii = 1; iii <= mworki; iii++) {
        iwork(iii) = 0;
    }
    //
    for (m = 1; m <= nnodes; m++) {
        for (nbl = 1; nbl <= nblock; nbl++) {
            if (mblk2nd(nbl) == m) {
                { int _iunit11 = 11; writ_buf_ns::writ_buf(nbl, _iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
                goto label_777;
            }
        }
        label_777:;
    }
    //
    // ***********************************************************************
    //   Read grid data and perform preliminary calculations.
    // ***********************************************************************
    //
    isklton = 1;
    //
    iwk1    = 1;
    iwk2    = iwk1 + maxgr;
    iwk3    = iwk2 + maxgr;
    iwk4    = iwk3 + maxgr;
    iwk5    = iwk4 + maxbl * 8;
    mworki1 = mworki - iwk5 + 1;
    if (mworki1 <= 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, "stopping...not enough integer work space for subroutine setup");
        { int _m1 = -1; termn8_ns::termn8(myid, _m1, ibufdim, nbuf, bou, nou); }
    }
    //
    // constants derived from gamma
    //
    gm1  = gamma - 1.0e0f;
    gp1  = gamma + 1.0e0f;
    gm1g = gm1 / gamma;
    gp1g = gp1 / gamma;
    ggm1 = gamma * gm1;
    // set nummem (usually 2 for 1 or 2 turbulence equations)
    // need to also set memory appropriately in pointers
    nummem = 2;
    if (ivmx == 30) nummem = 3;
    if (ivmx == 40) nummem = 4;
    if (ivmx == 72) nummem = (int)module_kwstm_ns::kws_get_nummem();
    rmstr_arr.allocate(ncycmax, nummem);
    { int _n = ncycmax * nummem, _itype = 0; char _txt[1]={'r'}; umalloc_ns::umalloc(_n, _itype, _txt, memuse, stats); }
    nneg_arr.allocate(ncycmax, nummem);
    { int _n = ncycmax * nummem, _itype = 1; char _txt[1]={'n'}; umalloc_ns::umalloc(_n, _itype, _txt, memuse, stats); }
    //
    {
        FortranArray1DRef<double> work_nstart1(work.data() + nstart, work.size(1) - nstart);
        FortranArray1DRef<int>    iwork_iwk5(iwork.data() + (iwk5 - 1), iwork.size(1) - iwk5 + 1);
        FortranArray1DRef<int>    iwork_iwk1(iwork.data() + (iwk1 - 1), iwork.size(1) - iwk1 + 1);
        FortranArray1DRef<int>    iwork_iwk2(iwork.data() + (iwk2 - 1), iwork.size(1) - iwk2 + 1);
        FortranArray1DRef<int>    iwork_iwk3(iwork.data() + (iwk3 - 1), iwork.size(1) - iwk3 + 1);
        FortranArray1DRef<int>    iwork_iwk4(iwork.data() + (iwk4 - 1), iwork.size(1) - iwk4 + 1);
        FortranArray2DRef<double> rmstr_ref = rmstr_arr.ref();
        FortranArray2DRef<int>    nneg_ref  = nneg_arr.ref();
        setup_ns::setup(lw, lw2, work, nstart, work_nstart1, nwork, iwork_iwk5,
                        mworki1, iwork_iwk1, iwork_iwk2, iwork_iwk3,
                        maxbl, mxbli, maxgr, maxseg, nsub1, maxxe, intmax,
                        iitot, ncycmax, lwdat, lig, lbg, iovrlp,
                        qb, nblock, iviscg, jdimg, kdimg, idimg, utrans,
                        vtrans, wtrans, omegax, omegay, omegaz, xorig,
                        yorig, zorig, dxmx, dymx, dzmx, dthxmx, dthymx,
                        dthzmx, thetax, thetay, thetaz, rfreqt, rfreqr,
                        xorig0, yorig0, zorig0, time2, thetaxl, thetayl,
                        thetazl, itrans, irotat, idefrm, bcvali, bcvalj,
                        bcvalk, nbci0, nbcidim, nbcj0, nbcjdim,
                        nbck0, nbckdim, ibcinfo, jbcinfo, kbcinfo, bcfilei,
                        bcfilej, bcfilek, ngrid, ncgg, nblg, iemg, inewgg,
                        rms, clw, cdw, cdpw, cdvw, cxw, cyw, czw, cmxw, cmyw,
                        cmzw, n_clcd, clcd, nblocks_clcd, blocks_clcd,
                        chdw, swetw, fmdotw, cfttotw, cftmomw, cftpw,
                        cftvw, rmstr_ref, nneg_ref, ntr, windex,
                        ninter, iindex, nblkpt, dthetxx, dthetyy, dthetzz,
                        iibg, kkbg, jjbg, ibcg, dxintg, dyintg, dzintg, iiig,
                        jjig, kkig, ibpntsg, iipntsg, mblk2nd, nou, bou, nbuf,
                        ibufdim, iwork_iwk4, igridg, bcfiles, mxbcfil,
                        utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae,
                        xorgae, yorgae, zorgae, thtxae, thtyae, thtzae,
                        rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf, kcsi, kcsf,
                        freq, gmass, damp, x0, gf0, nmds, maxaes, aesrfdat, perturb,
                        islavept, nslave, iskip, jskip, kskip, bmat, stm, stmi,
                        gforcn, gforcnm, xxn, nsegdfrm, idfrmseg, iaesurf,
                        maxsegdg, nmaster, aehist, timekeep, inpl3d, nplots,
                        nplot3d, levelg, iadvance, xs, gforcs, xorgae0, yorgae0,
                        zorgae0, icouple, lfgm, nblk, limblk, isva, nblelst, iskmax,
                        jskmax, kskmax, ue, irdrea, nbli, nummem);
    }
    //
    for (iii = 1; iii <= mworki; iii++) {
        iwork(iii) = 0;
    }
    for (iii = nstart + 1; iii <= mwork; iii++) {
        work(iii) = 0.0;
    }
    //
    // write header to cfl3d.alpha file if specified Cl option is used
    //
    if (ialphit != 0 && myid == myhost) {
        fortran_write_unit(27, "%s\n", "       it    log(res)           cl            alpha           dalpha");
    }
    //
    // open plot3d/printout files
    //
    if (myid == myhost) {
        char tmp3[81], tmp4[81], tmp17[81];
        // trim plt3dg, plt3dq, printout to null-terminated
        {
            int len = 80;
            while (len > 0 && plt3dg[len-1] == ' ') len--;
            std::memcpy(tmp3, plt3dg, len); tmp3[len] = '\0';
        }
        {
            int len = 80;
            while (len > 0 && plt3dq[len-1] == ' ') len--;
            std::memcpy(tmp4, plt3dq, len); tmp4[len] = '\0';
        }
        {
            int len = 80;
            while (len > 0 && printout[len-1] == ' ') len--;
            std::memcpy(tmp17, printout, len); tmp17[len] = '\0';
        }
        if (ibin == 1) {
            fortran_open_unit(3, tmp3, "wb");
            fortran_open_unit(4, tmp4, "wb");
        } else {
            fortran_open_unit(3, tmp3, "w");
            fortran_open_unit(4, tmp4, "w");
        }
        fortran_open_unit(17, tmp17, "w");
    }

    //
    // ***********************************************************************
    //  mesh sequencing - cycle through mesh sequences, coarser to finer
    // ***********************************************************************
    //
    dt0    = dt;
    nframes = 0;
    icall1  = 0;
    //
    for (iseq = 1; iseq <= mseq; iseq++) {
        //
        // skip mesh sequence levels with ncyc = 0
        //
        if (ncyc1[iseq - 1] > 0) {
            //
            // interpolate from coarser solution
            //
            if (iseq > 1 || iseq == mseq) {
                isklton = 1;
                qinter_ns::qinter(iseq, lembed, maxl, lw, lw2, work, nstart,
                                  FortranArray1DRef<double>(work.data() + nstart, work.size(1) - nstart),
                                  nwork, maxbl, maxgr,
                                  levelg, igridg, idimg, jdimg, kdimg, ngrid,
                                  nblg, iemg, inewgg, itrans, irotat, idefrm, xorig,
                                  yorig, zorig, xorig0, yorig0, zorig0, thetax, thetay,
                                  thetaz, time2, nou, bou, nbuf, ibufdim, myid, myhost,
                                  mycomm, mblk2nd, nsegdfrm, idfrmseg, xorgae,
                                  yorgae, zorgae, thtxae, thtyae, thtzae, maxsegdg,
                                  nummem);
                //
                for (iii = nstart + 1; iii <= mwork; iii++) {
                    work(iii) = 0.0;
                }
                //
            }
            //
            // ***********************************************************************
            //    multigrid/time-step cycles
            // ***********************************************************************
            //
            {
                std::memset(string, ' ', 50);
                string[50] = '\0';
                int _pp1 = +1, _iunit11 = 11;
                cputim_ns::cputim(_pp1, nnodes, string, myhost, myid, mycomm, _iunit11);
            }
            //
            nwork1 = nwork - maxbl;
            if (nwork1 <= 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, "stopping...not enough work space for subroutine mgblk");
                { int _m1 = -1; termn8_ns::termn8(myid, _m1, ibufdim, nbuf, bou, nou); }
            }
            iwk1 = 1;
            iwk2 = iwk1 + maxbl * 7 * 3;
            iwk3 = iwk2 + maxbl * 8;
            mworki1 = mworki - iwk3 + 1;
            if (mworki1 <= 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, "stopping...not enough integer work space for subroutine mgblk");
                { int _m1 = -1; termn8_ns::termn8(myid, _m1, ibufdim, nbuf, bou, nou); }
            }

            {
                // call mgblk(iseq,lw,lw2,work,nstart,work(nstart+1),
                //            work(nstart+1+maxbl),nwork1,
                //            iwork(iwk1),iwork(iwk2),iwork(iwk3),mworki1, ...)
                FortranArray1DRef<double> work_nstart1(work.data() + nstart, work.size(1) - nstart);
                FortranArray1DRef<double> work_nstart1_maxbl(work.data() + nstart + maxbl, work.size(1) - nstart - maxbl);
                FortranArray1DRef<int>    iwork_iwk1(iwork.data() + (iwk1 - 1), iwork.size(1) - iwk1 + 1);
                FortranArray1DRef<int>    iwork_iwk2(iwork.data() + (iwk2 - 1), iwork.size(1) - iwk2 + 1);
                FortranArray1DRef<int>    iwork_iwk3(iwork.data() + (iwk3 - 1), iwork.size(1) - iwk3 + 1);
                FortranArray2DRef<double> rmstr_ref = rmstr_arr.ref();
                FortranArray2DRef<int>    nneg_ref  = nneg_arr.ref();
                mgblk_ns::mgblk(iseq, lw, lw2, work, nstart, work_nstart1,
                                work_nstart1_maxbl, nwork1,
                                iwork_iwk1, iwork_iwk2, iwork_iwk3, mworki1,
                                maxbl, maxgr, maxseg, iitot, intmax, nsub1,
                                maxxe, ncycmax, iovrlp, lig, lbg, ibpntsg, iipntsg,
                                iibg, kkbg, jjbg, ibcg, dxintg, dyintg, dzintg, iiig,
                                jjig, kkig, qb, lwdat, nblk, nbli, limblk, isva, nblon,
                                nblock, levelg, igridg, iviscg, idimg, jdimg, kdimg,
                                jsg, ksg, isg, jeg, keg, ieg, nblcg, mit, bcvali, bcvalj,
                                bcvalk, nbci0, nbcidim, nbcj0, nbcjdim, nbck0,
                                nbckdim, ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej,
                                bcfilek, utrans, vtrans, wtrans, omegax, omegay,
                                omegaz, xorig, yorig, zorig, dxmx, dymx, dzmx, dthxmx,
                                dthymx, dthzmx, thetax, thetay, thetaz, rfreqt,
                                rfreqr, xorig0, yorig0, zorig0, time2, thetaxl,
                                thetayl, thetazl, itrans, irotat, idefrm,
                                rms, clw, cdw, cdpw, cdvw, cxw, cyw, czw, cmxw,
                                cmyw, cmzw, n_clcd, clcd, nblocks_clcd, blocks_clcd,
                                chdw, swetw, fmdotw, cfttotw, cftmomw,
                                cftpw, cftvw, rmstr_ref, nneg_ref, ntr,
                                ihstry, iadvance, iforce, lfgm, resmx, imx, jmx, kmx,
                                vormax, ivmax, jvmax, kvmax, sx, sy, sz, stot, pav,
                                ptav, tav, ttav, xmav, fmdot, cfxp, cfyp, cfzp, cfdp,
                                cflp, cftp, cfxv, cfyv, cfzv, cfdv, cflv, cftv, cfxmom,
                                cfymom, cfzmom, cfdmom, cflmom, cftmom, cfxtot,
                                cfytot, cfztot, cfdtot, cfltot, cfttot, icsinfo, ncs,
                                windex, ninter, iindex, nblkpt, windx, nintr, iindx,
                                llimit, iitmax, mmcxie, mmceta, ncheck, iifit, mblkpt,
                                iic0, iiorph, iitoss, ifiner, msub1, dthetxx, dthetyy,
                                dthetzz, swett, clt, cdt, cxt, cyt, czt, cmxt, cmyt,
                                cmzt, cdpt, cdvt, mblk2nd, geom_miss, epsc0,
                                period_miss, epsrot, isav_blk, isav_prd, lbcprd,
                                isav_pat, isav_pat_b, isav_emb, lbcemb, mxbli,
                                maxcs, intmx, mxxe, mptch, ncgg, nblg, iemg, ngrid,
                                dx, dy, dz, dthetx, dthety, dthetz, isav_dpat,
                                isav_dpat_b, lout, ifrom, xif1, xif2, etf1,
                                etf2, jjmax1, kkmax1, iiint1, iiint2, nblk1,
                                nblk2, jimage, kimage, jte, kte, jmm, kmm,
                                xte, yte, zte, xmi, ymi, zmi, xmie, ymie,
                                zmie, sxie, seta, sxie2, seta2, xie2s,
                                eta2s, temp, x2, y2, z2, x1, y1, z1, ip3dsurf,
                                factjlo, factjhi, factklo, factkhi, nplots,
                                nou, bou, nbuf, ibufdim, istat2_bl,
                                istat2_pa, istat2_pe, istat2_em, istat_size,
                                nplot3d, nprint, inpl3d, inpr, bcfiles, mxbcfil,
                                utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae,
                                xorgae, yorgae, zorgae, thtxae, thtyae, thtzae,
                                rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf, kcsi, kcsf,
                                freq, gmass, damp, x0, gf0, nmds, maxaes, aesrfdat, perturb,
                                islavept, nslave, iskip, jskip, kskip, bmat, stm, stmi,
                                xs, xxn, gforcn, gforcnm, gforcs, nsegdfrm, idfrmseg,
                                iaesurf, maxsegdg, nmaster, aehist, timekeep, xorgae0,
                                yorgae0, zorgae0, icouple, iprnsurf, nblelst, iskmax,
                                jskmax, kskmax, ue, nummem);
            }
            //
            for (iii = 1; iii <= mworki; iii++) {
                iwork(iii) = 0;
            }
            for (iii = nstart + 1; iii <= mwork; iii++) {
                work(iii) = 0.0;
            }
            //
            {
                std::memset(string, ' ', 50);
                string[50] = '\0';
                int _pp1 = +1, _iunit11 = 11;
                cputim_ns::cputim(_pp1, nnodes, string, myhost, myid, mycomm, _iunit11);
            }

            //
            // determine computational rate (note: ignores start-up overhead)
            // tim(3,2) contains the wall time between calls to cputim
            //
            // need a representative value for ncells for each sequence:
            // for non-embeded cases, include only those cells on the finest
            // global level; for embeded cases, include the cells on the finest
            // global level *and* the finest embeded level
            //
            lglobal = lfgm - (mseq - iseq);
            if (levelt[iseq - 1] == lglobal) {
                ncells = ncell[lglobal - 1];
            } else {
                ncells = ncell[lglobal - 1] + ncell[levelt[iseq - 1] - 1];
            }
            //
            // best attempt at correct wall time for > 24 hour runs:
            // tim(3,1) = tim_arr[2 + 0*3] = tim_arr[2]
            // tim(1,1) = tim_arr[0]
            // tim(2,1) = tim_arr[1]
            //
            walltime  = tim_arr[2];  // tim(3,1)
            totaltime = tim_arr[0] + tim_arr[1];  // tim(1,1) + tim(2,1)
            ndays = (int)((float)totaltime) / 86400;
            if (ndays > 0) {
                walltime = walltime + (float)(ndays * 86400);
            } else {
                walltime = std::max(walltime, totaltime);
            }
            if ((float)dt < 0.0f) {
                rate[iseq - 1] = 1.0e6f * walltime /
                                 (float)ncells / (float)ncyc1[iseq - 1];
            } else {
                rate[iseq - 1] = 1.0e6f * walltime /
                                 (float)ncells / (float)ntstep;
                ratesub[iseq - 1] = rate[iseq - 1] / (float)ncyc1[iseq - 1];
            }
            //
            // meshdef = 1 if flow solution is to be bipassed, and mesh deformation only:
            //
            if (meshdef == 1) goto label_8100;
            //
            // ***********************************************************************
            //    output for the current sequence level
            // ***********************************************************************
            //
            // output solution (plot3d and printout files)
            //
            if ((float)dt < 0.0f || ((float)dt > 0.0f && movie == 0)) {
                iflag = 0;
                if (iseq == mseq) iflag = 1;
                if (iseq < mseq && ncyc1[iseq] == 0) iflag = 1;
                if (iflag > 0) {
                    isklton = 0;
                    lhdr    = 1;
                    iwk1    = 1;
                    iwk2    = iwk1 + 3 * nplots;
                    iwk3    = iwk2 + maxbl;
                    mworki1 = mworki - iwk3 + 1;
                    if (mworki1 <= 0) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120, "stopping...not enough integer work space for subroutine qout");
                        { int _m1 = -1; termn8_ns::termn8(myid, _m1, ibufdim, nbuf, bou, nou); }
                    }
                    {
                        FortranArray1DRef<double> work_nstart1(work.data() + nstart, work.size(1) - nstart);
                        FortranArray1DRef<int>    iwork_iwk1(iwork.data() + (iwk1 - 1), iwork.size(1) - iwk1 + 1);
                        FortranArray1DRef<int>    iwork_iwk2(iwork.data() + (iwk2 - 1), iwork.size(1) - iwk2 + 1);
                        FortranArray1DRef<int>    iwork_iwk3(iwork.data() + (iwk3 - 1), iwork.size(1) - iwk3 + 1);
                        // ip3ddim is a 2D array (3, nplots) stored in iwork(iwk1)
                        FortranArray2DRef<int>    ip3ddim_ref(iwork.data() + (iwk1 - 1), 3, nplots);
                        qout_ns::qout(iseq, lw, lw2, work, nstart, work_nstart1, nwork,
                                      nplots, iovrlp, iibg, kkbg, jjbg, ibcg, lbg,
                                      ibpntsg, qb, lwdat, nbci0, nbcj0, nbck0,
                                      nbcidim, nbcjdim, nbckdim, jbcinfo, kbcinfo,
                                      ibcinfo, bcfilei, bcfilej, bcfilek, itrans,
                                      irotat, idefrm, nblock, levelg, igridg, iviscg,
                                      jdimg, kdimg, idimg, nblg, clw, ncycmax, nplot3d,
                                      inpl3d, ip3dsurf, nprint, inpr, iadvance, mycomm,
                                      myid, myhost, mblk2nd, nou, bou, nbuf, ibufdim, maxbl,
                                      maxgr, maxseg, iitot, jsg, ksg, isg, jeg, keg, ieg,
                                      ninter, windex, iindex, nblkpt, intmax, nsub1, maxxe,
                                      nblk, nbli, limblk, isva, nblon, mxbli, thetay,
                                      ip3ddim_ref, iwork_iwk2, iwork_iwk3, mworki1,
                                      xorig, yorig, zorig, period_miss, geom_miss,
                                      epsc0, epsrot, isav_blk, isav_pat, isav_pat_b,
                                      isav_emb, isav_prd, lbcprd,
                                      lbcemb, dthetxx, dthetyy, dthetzz, nblcg,
                                      lfgm, istat2_bl, istat2_pa, istat2_pe, istat2_em,
                                      istat_size, vormax, ivmax, jvmax, kvmax, bcfiles,
                                      mxbcfil, iprnsurf, nummem);
                    }
                    //
                    for (iii = 1; iii <= mworki; iii++) {
                        iwork(iii) = 0;
                    }
                    for (iii = nstart + 1; iii <= mwork; iii++) {
                        work(iii) = 0.0;
                    }
                    //
                } // end if (iflag > 0)
            } // end if dt < 0 or (dt > 0 and movie == 0)

            //
            // print out control surface data
            //
            iflag = 0;
            if (iseq == mseq) iflag = 1;
            if (iseq < mseq && ncyc1[iseq] == 0) iflag = 1;
            if (iflag > 0) {
                if (myid == myhost) {
                    csout_ns::csout(iseq, maxbl, maxcs, igridg, levelg, ncs, sx, sy, sz,
                                    stot, pav, ptav, tav, ttav, xmav, fmdot, cfxp, cfyp,
                                    cfzp, cfdp, cflp, cftp, cfxv, cfyv, cfzv, cfdv, cflv,
                                    cftv, cfxmom, cfymom, cfzmom, cfdmom, cflmom,
                                    cftmom, cfxtot, cfytot, cfztot, cfdtot, cfltot,
                                    cfttot, icsinfo);
                }
            }
            //
            // print out force and moment data
            //
            forceout_ns::forceout(iseq, maxbl, maxgr, maxseg, nblock, iforce,
                                  igridg, nbci0, nbcj0, nbck0, nbcidim, nbcjdim,
                                  nbckdim, levelg, ibcinfo, jbcinfo, kbcinfo,
                                  swett, clt, cdt, cxt, cyt, czt, cmxt, cmyt, cmzt,
                                  cdpt, cdvt, swetw, clw, cdw, cxw, cyw, czw, cmxw,
                                  cmyw, cmzw, cdpw, cdvw, ncycmax, myhost, myid,
                                  mycomm, mblk2nd);
            //
            // calculate and print out yplus statistics for turbulent flows
            //
            iflag = 0;
            if (iseq == mseq) iflag = 1;
            if (iseq < mseq && ncyc1[iseq] == 0) iflag = 1;
            if (ivmx < 2) iflag = 0;
            if (iflag > 0) {
                nttuse = std::max(ntt - 1, 1);
                iwk1    = 1;
                iwk2    = iwk1 + 6 * maxbl;
                iwk3    = iwk2 + 6 * maxbl;
                mworki1 = mworki - iwk3;
                if (mworki1 < 0) {
                    nou(1) = std::min(nou(1) + 1, ibufdim);
                    std::snprintf(bou(nou(1), 1), 120, "stopping...not enough integer work space for subroutine yplusout");
                    { int _m1 = -1; termn8_ns::termn8(myid, _m1, ibufdim, nbuf, bou, nou); }
                }
                {
                    FortranArray1DRef<double> work_nstart1(work.data() + nstart, work.size(1) - nstart);
                    FortranArray1DRef<int>    iwork_iwk1(iwork.data() + (iwk1 - 1), iwork.size(1) - iwk1 + 1);
                    FortranArray1DRef<int>    iwork_iwk2(iwork.data() + (iwk2 - 1), iwork.size(1) - iwk2 + 1);
                    FortranArray1DRef<int>    iwork_iwk3(iwork.data() + (iwk3 - 1), iwork.size(1) - iwk3 + 1);
                    double clw_nttuse = (double)clw(nttuse);
                    yplusout_ns::yplusout(iseq, lw, lw2, work, nstart, work_nstart1, nwork,
                                         clw_nttuse, maxbl, maxgr, maxseg, nblock,
                                         lwdat, levelg, igridg, jdimg, kdimg, idimg,
                                         nbci0, nbcj0, nbck0, nbcidim, nbcjdim,
                                         nbckdim, bcfilei, bcfilej, bcfilek,
                                         itrans, irotat, idefrm, nblg, ibcinfo, jbcinfo,
                                         kbcinfo, iadvance, iovrlp, myid, myhost,
                                         mycomm, mblk2nd, iwork_iwk1, iwork_iwk2,
                                         iwork_iwk3, nou, bou, nbuf, ibufdim, bcfiles,
                                         mxbcfil, nummem);
                }
                //
                for (iii = 1; iii <= mworki; iii++) {
                    iwork(iii) = 0;
                }
                for (iii = nstart + 1; iii <= mwork; iii++) {
                    work(iii) = 0.0;
                }
                //
            } // end if (iflag > 0) for yplusout
            //
        } // end if (ncyc1(iseq) > 0)
        //
    } // end do 8000 iseq=1,mseq
    label_8100:;
    //
    // end of mesh sequencing
    //
    // ***********************************************************************
    //  output convergence history
    // ***********************************************************************
    //
    // meshdef = 1 if flow solution is to be bipassed, and mesh deformation only:
    //
    if (meshdef == 1) goto label_8110;
    //
    {
        FortranArray2DRef<double> rmstr_ref = rmstr_arr.ref();
        FortranArray2DRef<int>    nneg_ref  = nneg_arr.ref();
        histout_ns::histout(ihstry, rms, clw, cdw, cdpw, cdvw, cxw, cyw, czw,
                            cmxw, cmyw, cmzw, n_clcd, clcd, nblocks_clcd, blocks_clcd,
                            chdw, swetw, fmdotw, cfttotw,
                            cftmomw, cftpw, cftvw, rmstr_ref, nneg_ref,
                            ncycmax, aehist, aesrfdat, nmds, maxaes,
                            timekeep, nummem);
    }
    //
    if (std::abs(movie) > 0) {
        // write(11,161) nframes
        // format 161: 7h output,i4,36h frames of plot3d data for animation
        fortran_write_unit(11, " output%4d frames of plot3d data for animation\n", nframes);
    }
    label_8110:;
    //
    // ***********************************************************************
    //  output computational rates based on *wall time*
    // ***********************************************************************
    //
    // write(11,9994)
    // format 9994: /1x,36hcomputational rate by mesh sequence ,21h(based on wall time):
    fortran_write_unit(11, "\n computational rate by mesh sequence (based on wall time):\n");
    for (iseq = 1; iseq <= mseq; iseq++) {
        if ((float)dt < 0.0f) {
            // format 9995: 1x,4hiseq,i2,1x,f7.2,28h microseconds/cell/iteration
            fortran_write_unit(11, " iseq%2d %7.2f microseconds/cell/iteration\n", iseq, (float)rate[iseq - 1]);
        } else {
            // format 9996: 1x,4hiseq,i2,1x,f7.2,28h microseconds/cell/time step
            fortran_write_unit(11, " iseq%2d %7.2f microseconds/cell/time step\n", iseq, (float)rate[iseq - 1]);
            // format 9997: 8x,f7.2,31h microseconds/cell/subiteration
            fortran_write_unit(11, "        %7.2f microseconds/cell/subiteration\n", (float)ratesub[iseq - 1]);
        }
    }
    //
    // ***********************************************************************
    //  output final timings for this run
    // ***********************************************************************
    //
    {
        char timing_str[51];
        std::memset(timing_str, ' ', 50);
        const char* ts = "    timing for complete run - time in seconds     ";
        std::memcpy(timing_str, ts, std::min((int)std::strlen(ts), 50));
        timing_str[50] = '\0';
        int _mm1 = -1, _iunit11 = 11;
        cputim_ns::cputim(_mm1, nnodes, timing_str, myhost, myid, mycomm, _iunit11);
    }
    //
    rmstr_arr.deallocate();
    nneg_arr.deallocate();
    return;
}

} // namespace mgbl_ns
