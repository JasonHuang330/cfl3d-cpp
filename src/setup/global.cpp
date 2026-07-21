// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "global.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include "bc.h"
#include "blocki.h"
#include "readkey.h"
#include "termn8.h"
#include "getibk0.h"
#include "getdhdr.h"
#include "rpatch0.h"
#include "ccomplex.h"
#include "setseg.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace global_ns {

// -----------------------------------------------------------------------
// bc: thin wrapper — actual implementation is in bc_ns
// -----------------------------------------------------------------------
void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork,
        double& cl, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
        int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg,
        FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
        FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg,
        FortranArray1DRef<int> nblg, FortranArray1DRef<int> nbci0,
        FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
        FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim,
        FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo,
        FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
        FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej,
        FortranArray3DRef<int> bcfilek, FortranArray3DRef<int> lwdat,
        int& myid, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
        FortranArray1DRef<int> kdimg, FortranArray1DRef<char[80]> bcfiles,
        int& mxbcfil, int& nummem)
{
    bc_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, cl, nou, bou, nbuf,
              ibufdim, maxbl, maxgr, maxseg, itrans, irotat, idefrm, igridg,
              nblg, nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
              ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek, lwdat,
              myid, idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
}

// -----------------------------------------------------------------------
// blocki: thin wrapper — actual implementation is in blocki_ns
// -----------------------------------------------------------------------
void blocki(FortranArray4DRef<double> q, FortranArray4DRef<double> qi0,
            int& idimr, int& jdimr, int& kdimr, int& jdimt, int& kdimt,
            FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva,
            int& it, int& ir, int& ldim, FortranArray3DRef<double> bci,
            int& iedge, int& ivolflg)
{
    blocki_ns::blocki(q, qi0, idimr, jdimr, kdimr, jdimt, kdimt,
                      limblk, isva, it, ir, ldim, bci, iedge, ivolflg);
}


// -----------------------------------------------------------------------
// global: main input reader
// -----------------------------------------------------------------------
void global(int& myid, int& maxbl, int& maxgr, int& maxseg, int& maxcs,
            int& nplots, int& mxbli,
            FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj,
            FortranArray4DRef<double> bcvalk,
            FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0,
            FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim,
            FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim,
            FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo,
            FortranArray4DRef<int> kbcinfo,
            FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej,
            FortranArray3DRef<int> bcfilek,
            FortranArray2DRef<int> nblk, int& nbli,
            FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva,
            FortranArray1DRef<int> nblon, FortranArray2DRef<double> rkap0g,
            int& nblock, FortranArray1DRef<int> levelg,
            FortranArray1DRef<int> igridg, FortranArray2DRef<int> iflimg,
            FortranArray2DRef<int> ifdsg, FortranArray2DRef<int> iviscg,
            FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg,
            FortranArray1DRef<int> idimg, FortranArray2DRef<int> idiagg,
            FortranArray1DRef<int> nblcg, FortranArray2DRef<int> idegg,
            FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg,
            FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg,
            FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg,
            FortranArray2DRef<int> mit, FortranArray1DRef<int> ilamlog,
            FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamlog,
            FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamlog,
            FortranArray1DRef<int> klamhig, FortranArray2DRef<int> iwfg,
            FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans,
            FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax,
            FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz,
            FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig,
            FortranArray1DRef<double> zorig, FortranArray1DRef<double> dxmx,
            FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx,
            FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx,
            FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> thetax,
            FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz,
            FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr,
            FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0,
            FortranArray1DRef<double> zorig0,
            FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
            FortranArray1DRef<int> idefrm,
            int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg,
            FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg,
            FortranArray1DRef<int> iovrlp, int& ninter, int& nplot3d,
            FortranArray2DRef<int> inpl3d, int& ip3dsurf, int& nprint,
            FortranArray2DRef<int> inpr, FortranArray1DRef<int> iadvance,
            FortranArray1DRef<int> iforce, int& lfgm, int& ncs,
            FortranArray2DRef<int> icsinfo, int& ihstry, int& ncycmax,
            FortranArray1DRef<int> iv, FortranArray1DRef<double> time2,
            FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl,
            FortranArray1DRef<double> thetazl, int& intmax, int& nsub1,
            FortranArray2DRef<int> iindex, FortranArray1DRef<int> lig,
            FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg,
            FortranArray1DRef<int> iipntsg, int& icall, int& iunit11,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& ibufdim, int& nbuf, FortranArray1DRef<int> mglevg,
            FortranArray1DRef<int> nemgl, FortranArray2DRef<int> ipl3dtmp,
            int& ntr, FortranArray1DRef<char[80]> bcfiles, int& mxbcfil,
            FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae,
            FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae,
            FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae,
            FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae,
            FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae,
            FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae,
            FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae,
            FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
            FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
            FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
            FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass,
            FortranArray2DRef<double> damp, FortranArray2DRef<double> x0,
            FortranArray2DRef<double> gf0, int& nmds, int& maxaes,
            FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb,
            FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip,
            FortranArray2DRef<int> kskip, FortranArray1DRef<int> nsegdfrm,
            FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf,
            int& maxsegdg, FortranArray2DRef<double> xorgae0,
            FortranArray2DRef<double> yorgae0, FortranArray2DRef<double> zorgae0,
            FortranArray2DRef<int> icouple, int& iprnsurf)
{

    // COMMON block aliases
    float& xnumavg    = cmn_avgdata.xnumavg;
    int&   iteravg    = cmn_avgdata.iteravg;
    float& xnumavg2   = cmn_avgdata.xnumavg2;
    int&   ipertavg   = cmn_avgdata.ipertavg;
    int&   iclcd      = cmn_avgdata.iclcd;
    int&   isubit_r   = cmn_avgdata.isubit_r;
    float& dt0        = cmn_cfl.dt0;
    float& dtold      = cmn_cfl.dtold;
    float& xmach_img  = cmn_complx.xmach_img;
    float& alpha_img  = cmn_complx.alpha_img;
    float& beta_img   = cmn_complx.beta_img;
    float& reue_img   = cmn_complx.reue_img;
    float& tinf_img   = cmn_complx.tinf_img;
    float& geom_img   = cmn_complx.geom_img;
    float& surf_img   = cmn_complx.surf_img;
    float& xrotrate_img = cmn_complx.xrotrate_img;
    float& yrotrate_img = cmn_complx.yrotrate_img;
    float& zrotrate_img = cmn_complx.zrotrate_img;
    float& cdes       = cmn_des.cdes;
    int&   ides       = cmn_des.ides;
    float& cddes      = cmn_des.cddes;
    int&   ndefrm     = cmn_elastic.ndefrm;
    int&   naesrf     = cmn_elastic.naesrf;
    int&   idef_ss    = cmn_elastic_ss.idef_ss;
    int&   lmaxgr     = cmn_params.lmaxgr;
    int&   lmaxbl     = cmn_params.lmaxbl;
    int&   lmxseg     = cmn_params.lmxseg;
    int&   lmaxcs     = cmn_params.lmaxcs;
    int&   lnplts     = cmn_params.lnplts;
    int&   lmxbli     = cmn_params.lmxbli;
    int&   lmaxxe     = cmn_params.lmaxxe;
    int&   lnsub1     = cmn_params.lnsub1;
    int&   lintmx     = cmn_params.lintmx;
    int&   lmxxe      = cmn_params.lmxxe;
    int&   liitot     = cmn_params.liitot;
    int&   isum       = cmn_params.isum;
    int&   lncycm     = cmn_params.lncycm;
    int&   isum_n     = cmn_params.isum_n;
    int&   lminnode   = cmn_params.lminnode;
    int&   isumi      = cmn_params.isumi;
    int&   isumi_n    = cmn_params.isumi_n;
    int&   lmptch     = cmn_params.lmptch;
    int&   lmsub1     = cmn_params.lmsub1;
    int&   lintmax    = cmn_params.lintmax;
    int&   libufdim   = cmn_params.libufdim;
    int&   lnbuf      = cmn_params.lnbuf;
    int&   llbcprd    = cmn_params.llbcprd;
    int&   llbcemb    = cmn_params.llbcemb;
    int&   llbcrad    = cmn_params.llbcrad;
    int&   lnmds      = cmn_params.lnmds;
    int&   lmaxaes    = cmn_params.lmaxaes;
    int&   lnslave    = cmn_params.lnslave;
    int&   lmxsegdg   = cmn_params.lmxsegdg;
    int&   lnmaster   = cmn_params.lnmaster;
    float& cprec      = cmn_precond.cprec;
    float& uref       = cmn_precond.uref;
    float& avn        = cmn_precond.avn;
    float* rate       = cmn_cpurate.rate;
    float* ratesub    = cmn_cpurate.ratesub;
    int*   ncell      = cmn_cpurate.ncell;
    int&   ichk       = cmn_chk.ichk;
    int*   ideg       = cmn_degshf.ideg;
    float& sref       = cmn_fsum.sref;
    float& cref       = cmn_fsum.cref;
    float& bref       = cmn_fsum.bref;
    float& xmc        = cmn_fsum.xmc;
    float& ymc        = cmn_fsum.ymc;
    float& zmc        = cmn_fsum.zmc;
    float* rkap0      = cmn_fvfds.rkap0;
    int*   ifds       = cmn_fvfds.ifds;
    float* title_arr  = cmn_info.title;
    float* rkap       = cmn_info.rkap;
    float& xmach      = cmn_info.xmach;
    float& alpha      = cmn_info.alpha;
    float& beta       = cmn_info.beta;
    float& dt         = cmn_info.dt;
    float& fmax       = cmn_info.fmax;
    int&   nit        = cmn_info.nit;
    int&   ntt        = cmn_info.ntt;
    int*   idiag      = cmn_info.idiag;
    int&   nitfo      = cmn_info.nitfo;
    int&   iflagts    = cmn_info.iflagts;
    int*   iflim      = cmn_info.iflim;
    int&   nres       = cmn_info.nres;
    int*   levelb     = cmn_info.levelb;
    int&   mgflag     = cmn_info.mgflag;
    int&   iconsf     = cmn_info.iconsf;
    int&   mseq       = cmn_info.mseq;
    int*   ncyc1      = cmn_info.ncyc1;
    int*   levelt     = cmn_info.levelt;
    int*   nitfo1     = cmn_info.nitfo1;
    int&   ngam       = cmn_info.ngam;
    int*   nsm        = cmn_info.nsm;
    int&   iipv       = cmn_info.iipv;
    int&   ilamlo     = cmn_lam.ilamlo;
    int&   ilamhi     = cmn_lam.ilamhi;
    int&   jlamlo     = cmn_lam.jlamlo;
    int&   jlamhi     = cmn_lam.jlamhi;
    int&   klamlo     = cmn_lam.klamlo;
    int&   klamhi     = cmn_lam.klamhi;
    int&   i_lam_forcezero = cmn_lam.i_lam_forcezero;
    int&   ivmx       = cmn_maxiv.ivmx;
    int&   levt       = cmn_mgrd.levt;
    int&   kode       = cmn_mgrd.kode;
    int&   mode       = cmn_mgrd.mode;
    int&   ncyc       = cmn_mgrd.ncyc;
    int&   mtt        = cmn_mgrd.mtt;
    int&   icyc       = cmn_mgrd.icyc;
    int&   level      = cmn_mgrd.level;
    int&   lglobal    = cmn_mgrd.lglobal;
    float* epsssc     = cmn_mgv.epsssc;
    float* epsssr     = cmn_mgv.epsssr;
    int&   issc       = cmn_mgv.issc;
    int&   issr       = cmn_mgv.issr;
    int&   ncyctot    = cmn_ncyct.ncyctot;
    int&   ialphit    = cmn_alphait.ialphit;
    float& cltarg     = cmn_alphait.cltarg;
    float& rlxalph    = cmn_alphait.rlxalph;
    float& dalim      = cmn_alphait.dalim;
    float& dalpha     = cmn_alphait.dalpha;
    int&   icycupdt   = cmn_alphait.icycupdt;
    float& reue       = cmn_reyue.reue;
    float& tinf       = cmn_reyue.tinf;
    int*   ivisc      = cmn_reyue.ivisc;
    int&   isminc     = cmn_sminn.isminc;
    int&   ismincforce = cmn_sminn.ismincforce;
    int&   i2d        = cmn_twod.i2d;
    float& time       = cmn_unst.time;
    float& cfltau     = cmn_unst.cfltau;
    int&   ntstep     = cmn_unst.ntstep;
    int&   ita        = cmn_unst.ita;
    int&   iunst      = cmn_unst.iunst;
    float& cfltau0    = cmn_unst.cfltau0;
    float& cfltauMax  = cmn_unst.cfltaumax;
    int*   iwf        = cmn_wallfun.iwf;
    int&   nwrest     = cmn_wrbl.nwrest;
    int&   irest      = cmn_wrestq.irest;
    int&   irest2     = cmn_wrestq.irest2;
    int&   movie      = cmn_moov.movie;
    int&   nframes    = cmn_moov.nframes;
    int&   icall1     = cmn_moov.icall1;
    int&   lhdr       = cmn_moov.lhdr;
    int&   icoarsemovie = cmn_moov.icoarsemovie;
    int&   i2dmovie   = cmn_moov.i2dmovie;
    int&   ip3dgrd    = cmn_igrdtyp.ip3dgrd;
    int&   ialph      = cmn_igrdtyp.ialph;
    float& radtodeg   = cmn_conversion.radtodeg;
    int&   iunit5     = cmn_unit5.iunit5;

    float& xmc0       = cmn_motionmc.xmc0;
    float& ymc0       = cmn_motionmc.ymc0;
    float& zmc0       = cmn_motionmc.zmc0;
    float& utransmc   = cmn_motionmc.utransmc;
    float& vtransmc   = cmn_motionmc.vtransmc;
    float& wtransmc   = cmn_motionmc.wtransmc;
    float& omegaxmc   = cmn_motionmc.omegaxmc;
    float& omegaymc   = cmn_motionmc.omegaymc;
    float& omegazmc   = cmn_motionmc.omegazmc;
    float& xorigmc    = cmn_motionmc.xorigmc;
    float& yorigmc    = cmn_motionmc.yorigmc;
    float& zorigmc    = cmn_motionmc.zorigmc;
    float& xorig0mc   = cmn_motionmc.xorig0mc;
    float& yorig0mc   = cmn_motionmc.yorig0mc;
    float& zorig0mc   = cmn_motionmc.zorig0mc;
    float& thetaxmc   = cmn_motionmc.thetaxmc;
    float& thetaymc   = cmn_motionmc.thetaymc;
    float& thetazmc   = cmn_motionmc.thetazmc;
    float& dxmxmc     = cmn_motionmc.dxmxmc;
    float& dymxmc     = cmn_motionmc.dymxmc;
    float& dzmxmc     = cmn_motionmc.dzmxmc;
    float& dthxmxmc   = cmn_motionmc.dthxmxmc;
    float& dthymxmc   = cmn_motionmc.dthymxmc;
    float& dthzmxmc   = cmn_motionmc.dthzmxmc;
    float& rfreqtmc   = cmn_motionmc.rfreqtmc;
    float& rfreqrmc   = cmn_motionmc.rfreqrmc;
    int&   itransmc   = cmn_motionmc.itransmc;
    int&   irotatmc   = cmn_motionmc.irotatmc;
    float& time2mc    = cmn_motionmc.time2mc;
    float& beta1      = cmn_deformz.beta1;
    float& beta2      = cmn_deformz.beta2;
    float& alpha1     = cmn_deformz.alpha1;
    float& alpha2     = cmn_deformz.alpha2;
    int&   isktyp     = cmn_deformz.isktyp;
    int&   negvol     = cmn_deformz.negvol;
    int&   meshdef    = cmn_deformz.meshdef;
    int&   nsprgit    = cmn_deformz.nsprgit;
    int&   ndgrd      = cmn_deformz.ndgrd;
    int&   ndwrt      = cmn_deformz.ndwrt;
    float* cflturb    = cmn_turbconv.cflturb;
    float& edvislim   = cmn_turbconv.edvislim;
    int&   iturbprod  = cmn_turbconv.iturbprod;
    int&   nsubturb   = cmn_turbconv.nsubturb;
    int&   nfreeze    = cmn_turbconv.nfreeze;
    int&   iwarneddy  = cmn_turbconv.iwarneddy;
    int&   itime2read = cmn_turbconv.itime2read;
    int&   itaturb    = cmn_turbconv.itaturb;
    float& tur1cut    = cmn_turbconv.tur1cut;
    float& tur2cut    = cmn_turbconv.tur2cut;
    int&   iturbord   = cmn_turbconv.iturbord;
    float& tur1cutlev = cmn_turbconv.tur1cutlev;
    float& tur2cutlev = cmn_turbconv.tur2cutlev;
    int&   ikoprod    = cmn_konew.ikoprod;
    int&   isstdenom  = cmn_konew.isstdenom;
    float& pklimterm  = cmn_konew.pklimterm;
    int&   ibeta8kzeta = cmn_konew.ibeta8kzeta;
    int&   i_bsl      = cmn_konew.i_bsl;
    int&   keepambient = cmn_konew.keepambient;
    float& re_thetat0 = cmn_konew.re_thetat0;
    int&   i_wilcox06 = cmn_konew.i_wilcox06;
    int&   i_wilcox06_chiw = cmn_konew.i_wilcox06_chiw;
    int&   i_turbprod_kterm = cmn_konew.i_turbprod_kterm;
    int&   i_catris_kw = cmn_konew.i_catris_kw;
    float& prod2d3dtrace = cmn_konew.prod2d3dtrace;
    int&   i_compress_correct = cmn_konew.i_compress_correct;
    int&   isstsf     = cmn_konew.isstsf;
    int&   i_wilcox98 = cmn_konew.i_wilcox98;
    int&   i_wilcox98_chiw = cmn_konew.i_wilcox98_chiw;
    int&   isst2003   = cmn_konew.isst2003;
    int&   iaxi2plane = cmn_axisym.iaxi2plane;
    int&   iaxi2planeturb = cmn_axisym.iaxi2planeturb;
    int&   istrongturbdis = cmn_axisym.istrongturbdis;
    int&   iforcev0   = cmn_axisym.iforcev0;
    int&   ifullns    = cmn_fullns.ifullns;
    int&   isarc2d    = cmn_curvat.isarc2d;
    float& sarccr3    = cmn_curvat.sarccr3;
    int&   ieasmcc2d  = cmn_curvat.ieasmcc2d;
    int&   isstrc     = cmn_curvat.isstrc;
    float& sstrc_crc  = cmn_curvat.sstrc_crc;
    int&   isar       = cmn_curvat.isar;
    float& crot       = cmn_curvat.crot;
    int&   isarc3d    = cmn_curvat.isarc3d;
    int&   iexact_trunc = cmn_mms.iexact_trunc;
    int&   iexact_disc  = cmn_mms.iexact_disc;
    int&   iexact_ring  = cmn_mms.iexact_ring;
    int&   i_nonlin   = cmn_constit.i_nonlin;
    float& c_nonlin   = cmn_constit.c_nonlin;
    float& snonlin_lim = cmn_constit.snonlin_lim;
    int&   i_tauijs   = cmn_constit.i_tauijs;
    int&   i_qcr2000  = cmn_constit.i_qcr2000;
    int&   i_qcr2013  = cmn_constit.i_qcr2013;
    int&   i_qcr2013v = cmn_constit.i_qcr2013v;
    float& cmulim     = cmn_easmlim.cmulim;
    float& p0         = cmn_ivals.p0;
    float& rho0       = cmn_ivals.rho0;
    float& c0         = cmn_ivals.c0;
    float& u0         = cmn_ivals.u0;
    float& v0         = cmn_ivals.v0;
    float& w0         = cmn_ivals.w0;
    float& et0        = cmn_ivals.et0;
    float& h0         = cmn_ivals.h0;
    float& pt0        = cmn_ivals.pt0;
    float& rhot0      = cmn_ivals.rhot0;
    float* qiv        = cmn_ivals.qiv;
    float* tur10      = cmn_ivals.tur10;
    int&   icgns      = cmn_cgns.icgns;
    int&   iccg       = cmn_cgns.iccg;
    int&   ibase      = cmn_cgns.ibase;
    int&   nzones     = cmn_cgns.nzones;
    int&   nsoluse    = cmn_cgns.nsoluse;
    int&   irind      = cmn_cgns.irind;
    int&   jrind      = cmn_cgns.jrind;
    int&   krind      = cmn_cgns.krind;
    int&   ifort50write = cmn_writestuff.ifort50write;
    int&   j_ifort50write = cmn_writestuff.j_ifort50write;
    int&   i_ifort50write = cmn_writestuff.i_ifort50write;
    int&   i_saneg    = cmn_sa_options.i_saneg;
    int&   i_sanoft2  = cmn_sa_options.i_sanoft2;
    float& sa_cw2     = cmn_sa_options.sa_cw2;
    float& sa_cw3     = cmn_sa_options.sa_cw3;
    float& sa_cv1     = cmn_sa_options.sa_cv1;
    float& sa_ct3     = cmn_sa_options.sa_ct3;
    float& sa_ct4     = cmn_sa_options.sa_ct4;
    float& sa_cb1     = cmn_sa_options.sa_cb1;
    float& sa_cb2     = cmn_sa_options.sa_cb2;
    float& sa_sigma   = cmn_sa_options.sa_sigma;
    float& sa_karman  = cmn_sa_options.sa_karman;
    int&   i_specialtop_kmax1001 = cmn_specialtop_kmax1001.i_specialtop_kmax1001;
    float& a_specialtop_kmax1001 = cmn_specialtop_kmax1001.a_specialtop_kmax1001;
    float& xc_specialtop_kmax1001 = cmn_specialtop_kmax1001.xc_specialtop_kmax1001;
    float& sig_specialtop_kmax1001 = cmn_specialtop_kmax1001.sig_specialtop_kmax1001;
    float& vtp_specialtop_kmax1001 = cmn_specialtop_kmax1001.vtp_specialtop_kmax1001;
    float& wc_specialtop_kmax1001 = cmn_specialtop_kmax1001.wc_specialtop_kmax1001;
    float& fac_specialtop_kmax1001 = cmn_specialtop_kmax1001.fac_specialtop_kmax1001;
    float& cc_specialtop_kmax1001 = cmn_specialtop_kmax1001.cc_specialtop_kmax1001;
    float& xerf_specialtop_kmax1001 = cmn_specialtop_kmax1001.xerf_specialtop_kmax1001;
    float& sigerf_specialtop_kmax1001 = cmn_specialtop_kmax1001.sigerf_specialtop_kmax1001;
    int&   issglrrw2012 = cmn_reystressmodel.issglrrw2012;
    int&   i_sas_rsm  = cmn_reystressmodel.i_sas_rsm;
    int&   i_yapterm  = cmn_reystressmodel.i_yapterm;
    int&   iupdatemean = cmn_iupdate.iupdatemean;


    // Local variables
    float pi;
    int itoti, itotb, nfiles, ierrflg;
    float realval[21]; // 1-based: use [1..20]
    int im[12];        // 1-based: use [1..11]
    // local ivisc array (separate from common /reyue/ ivisc)
    int ivisc_loc[4];  // 1-based: use [1..3]
    int idiag_loc[4];  // 1-based: use [1..3]
    int iflim_loc[4];  // 1-based: use [1..3]
    int ideg_loc[4];   // 1-based: use [1..3]
    int iwf_loc[4];    // 1-based: use [1..3]
    // allocatable arrays
    FortranArray2D<float> riskp, rjskp, rkskp;
    FortranArray2D<int>   iskipt, jskipt, kskipt;
    // local integer arrays (1-based)
    FortranArray1D<int> nskpi1(maxbl), nskpj1(maxbl), nskpk1(maxbl);

    // local scalars
    int ititr, nchk, ncg, iem, iad, ifor;
    int igrid, nbl, n, m, l, ib, ind;
    int ncgmax, iemtot;
    float uinf, vinf, winf, temp;
    int ntest, ncs2, msegment, mxmds, mxaes, mxdefseg;
    int nintr;
    int lf, lfem, lcgm, levelc;
    int igptot, igpts, ncellseq;
    int nblict, ntblk, iflag;
    int i1b1, i2b1, i1b2, i2b2;
    int idif1, idif2, idif3, idif4;
    int lpoint;
    int itest1, itest2;
    float lreftra, lrefrot, lrefdef;
    int iseg, nseg;
    int ibl;
    int n1, n2;
    int nblk1, nblk2;
    int lmblk11, lmblk12, lmblk21, lmblk22;
    int icrd1, icrd2;
    int ijkdim1, ijkdim2;
    int nskp1, nskp2, nskp3;
    int nskipt1, nskipt2, nskipt3;
    int iptype;
    int idim, jdim, kdim;
    int inewg, igridc, is, js, ks, ie, je, ke;
    int ig, iover;
    int ibctyp, ndata;
    int itemp1, itemp2;
    int ifoseg, mfoseg;
    int ifo, jfo, kfo;
    int istop;
    int ii, jj, kk;
    int nbl2, nbl3;
    int ntr_loc;
    int ncyctot_loc;
    int lfgm_loc, lfem_loc, lcgm_loc;
    int ncg_loc;
    int iemtot_loc;
    int nblock_loc;
    int nbl_loc;
    int igrid_loc;
    int ncgmax_loc;
    int ivmx_loc;
    int nblict_loc;
    int ntblk_loc;
    int iflag_loc;
    int i1b1_loc, i2b1_loc, i1b2_loc, i2b2_loc;
    int idif1_loc, idif2_loc, idif3_loc, idif4_loc;
    int lpoint_loc;
    int itest1_loc, itest2_loc;
    int igptot_loc, igpts_loc;
    int ncellseq_loc;
    int levelc_loc;
    int lf_loc;
    int n1_loc, n2_loc;
    int nblk1_loc, nblk2_loc;
    int lmblk11_loc, lmblk12_loc, lmblk21_loc, lmblk22_loc;
    int icrd1_loc, icrd2_loc;
    int ijkdim1_loc, ijkdim2_loc;
    int nskp1_loc, nskp2_loc, nskp3_loc;
    int nskipt1_loc, nskipt2_loc, nskipt3_loc;
    int ibl_loc;
    int iseg_loc, nseg_loc;
    int iptype_loc;
    int iunit5_loc;
    int ncs2_loc, mxmds_loc, mxaes_loc, mxdefseg_loc;
    int msegment_loc, ntest_loc;
    int nbl_tmp, igrid_tmp;
    int ncg_tmp;
    int iemtot_tmp;
    int nblock_tmp;
    int ncgmax_tmp;
    int ivmx_tmp;
    int nblict_tmp;
    int ntblk_tmp;
    int iflag_tmp;
    int i1b1_tmp, i2b1_tmp, i1b2_tmp, i2b2_tmp;
    int idif1_tmp, idif2_tmp, idif3_tmp, idif4_tmp;
    int lpoint_tmp;
    int itest1_tmp, itest2_tmp;
    int igptot_tmp, igpts_tmp;
    int ncellseq_tmp;
    int levelc_tmp;
    int lf_tmp;
    int lfem_tmp, lcgm_tmp;
    int n1_tmp, n2_tmp;
    int nblk1_tmp, nblk2_tmp;
    int lmblk11_tmp, lmblk12_tmp, lmblk21_tmp, lmblk22_tmp;
    int icrd1_tmp, icrd2_tmp;
    int ijkdim1_tmp, ijkdim2_tmp;
    int nskp1_tmp, nskp2_tmp, nskp3_tmp;
    int nskipt1_tmp, nskipt2_tmp, nskipt3_tmp;
    int ibl_tmp;
    int iseg_tmp, nseg_tmp;
    int iptype_tmp;
    int nfiles_tmp;
    int iunit5_tmp;
    int ntr_tmp;
    int ncyctot_tmp;
    int lfgm_tmp, lfem_tmp2, lcgm_tmp2;
    int nblt, mm, nfl;
    char datahdr[11][10]; // 1-based [1..10], each 10 chars (Fortran character*10)
    int stats = 0;
    int ncs2_tmp, mxmds_tmp, mxaes_tmp, mxdefseg_tmp;
    int msegment_tmp, ntest_tmp;
    int ntr_loc2, ncyctot_loc2;
    int lfgm_loc2, lfem_loc2, lcgm_loc2;
    int ncg_loc2, iemtot_loc2, nblock_loc2, nbl_loc2, igrid_loc2;
    int ncgmax_loc2, ivmx_loc2, nblict_loc2, ntblk_loc2, iflag_loc2;
    int i1b1_loc2, i2b1_loc2, i1b2_loc2, i2b2_loc2;
    int idif1_loc2, idif2_loc2, idif3_loc2, idif4_loc2;
    int lpoint_loc2, itest1_loc2, itest2_loc2;
    int igptot_loc2, igpts_loc2, ncellseq_loc2, levelc_loc2, lf_loc2;
    int n1_loc2, n2_loc2, nblk1_loc2, nblk2_loc2;
    int lmblk11_loc2, lmblk12_loc2, lmblk21_loc2, lmblk22_loc2;
    int icrd1_loc2, icrd2_loc2, ijkdim1_loc2, ijkdim2_loc2;
    int nskp1_loc2, nskp2_loc2, nskp3_loc2;
    int nskipt1_loc2, nskipt2_loc2, nskipt3_loc2;
    int ibl_loc2, iseg_loc2, nseg_loc2, iptype_loc2;
    int iunit5_loc2, ncs2_loc2, mxmds_loc2, mxaes_loc2, mxdefseg_loc2;
    int msegment_loc2, ntest_loc2;
    int nbl_tmp2, igrid_tmp2, ncg_tmp2, iemtot_tmp2, nblock_tmp2, ncgmax_tmp2;
    int ivmx_tmp2, nblict_tmp2, ntblk_tmp2, iflag_tmp2;
    int i1b1_tmp2, i2b1_tmp2, i1b2_tmp2, i2b2_tmp2;
    int idif1_tmp2, idif2_tmp2, idif3_tmp2, idif4_tmp2;
    int lpoint_tmp2, itest1_tmp2, itest2_tmp2;
    int n1_tmp2, n2_tmp2, nblk1_tmp2, nblk2_tmp2;

    int lmblk11_tmp2, lmblk12_tmp2, lmblk21_tmp2, lmblk22_tmp2;
    int icrd1_tmp2, icrd2_tmp2, ijkdim1_tmp2, ijkdim2_tmp2;
    int nskp1_tmp2, nskp2_tmp2, nskp3_tmp2;
    int nskipt1_tmp2, nskipt2_tmp2, nskipt3_tmp2;
    int ibl_tmp2, iseg_tmp2, nseg_tmp2, iptype_tmp2;
    int nfiles_tmp2, iunit5_tmp2, ntr_tmp2, ncyctot_tmp2;
    int lfgm_tmp2, lfem_tmp3, lcgm_tmp3;
    int nbl_tmp3, igrid_tmp3, ncg_tmp3, iemtot_tmp3, nblock_tmp3, ncgmax_tmp3;
    int ivmx_tmp3, nblict_tmp3, ntblk_tmp3, iflag_tmp3;
    int i1b1_tmp3, i2b1_tmp3, i1b2_tmp3, i2b2_tmp3;
    int idif1_tmp3, idif2_tmp3, idif3_tmp3, idif4_tmp3;
    int lpoint_tmp3, itest1_tmp3, itest2_tmp3;
    int igptot_tmp3, igpts_tmp3, ncellseq_tmp3, levelc_tmp3, lf_tmp3;
    int lfem_tmp4, lcgm_tmp4;
    int n1_tmp3, n2_tmp3, nblk1_tmp3, nblk2_tmp3;
    int lmblk11_tmp3, lmblk12_tmp3, lmblk21_tmp3, lmblk22_tmp3;
    int icrd1_tmp3, icrd2_tmp3, ijkdim1_tmp3, ijkdim2_tmp3;
    int nskp1_tmp3, nskp2_tmp3, nskp3_tmp3;
    int nskipt1_tmp3, nskipt2_tmp3, nskipt3_tmp3;
    int ibl_tmp3, iseg_tmp3, nseg_tmp3, iptype_tmp3;
    int nfiles_tmp3, iunit5_tmp3, ntr_tmp3, ncyctot_tmp3;
    int lfgm_tmp3, lfem_tmp5, lcgm_tmp5;
    int nbl_tmp4, igrid_tmp4, ncg_tmp4, iemtot_tmp4, nblock_tmp4, ncgmax_tmp4;
    int ivmx_tmp4, nblict_tmp4, ntblk_tmp4, iflag_tmp4;
    int i1b1_tmp4, i2b1_tmp4, i1b2_tmp4, i2b2_tmp4;
    int idif1_tmp4, idif2_tmp4, idif3_tmp4, idif4_tmp4;
    int lpoint_tmp4, itest1_tmp4, itest2_tmp4;
    int igptot_tmp4, igpts_tmp4, ncellseq_tmp4, levelc_tmp4, lf_tmp4;
    int lfem_tmp6, lcgm_tmp6;
    int n1_tmp4, n2_tmp4, nblk1_tmp4, nblk2_tmp4;
    int lmblk11_tmp4, lmblk12_tmp4, lmblk21_tmp4, lmblk22_tmp4;
    int icrd1_tmp4, icrd2_tmp4, ijkdim1_tmp4, ijkdim2_tmp4;
    int nskp1_tmp4, nskp2_tmp4, nskp3_tmp4;
    int nskipt1_tmp4, nskipt2_tmp4, nskipt3_tmp4;
    int ibl_tmp4, iseg_tmp4, nseg_tmp4, iptype_tmp4;
    int nfiles_tmp4, iunit5_tmp4, ntr_tmp4, ncyctot_tmp4;
    int lfgm_tmp4, lfem_tmp7, lcgm_tmp7;
    int nbl_tmp5, igrid_tmp5, ncg_tmp5, iemtot_tmp5, nblock_tmp5, ncgmax_tmp5;
    int ivmx_tmp5, nblict_tmp5, ntblk_tmp5, iflag_tmp5;
    int i1b1_tmp5, i2b1_tmp5, i1b2_tmp5, i2b2_tmp5;
    int idif1_tmp5, idif2_tmp5, idif3_tmp5, idif4_tmp5;
    int lpoint_tmp5, itest1_tmp5, itest2_tmp5;
    int igptot_tmp5, igpts_tmp5, ncellseq_tmp5, levelc_tmp5, lf_tmp5;
    int lfem_tmp8, lcgm_tmp8;
    int n1_tmp5, n2_tmp5, nblk1_tmp5, nblk2_tmp5;
    int lmblk11_tmp5, lmblk12_tmp5, lmblk21_tmp5, lmblk22_tmp5;
    int icrd1_tmp5, icrd2_tmp5, ijkdim1_tmp5, ijkdim2_tmp5;
    int nskp1_tmp5, nskp2_tmp5, nskp3_tmp5;
    int nskipt1_tmp5, nskipt2_tmp5, nskipt3_tmp5;
    int ibl_tmp5, iseg_tmp5, nseg_tmp5, iptype_tmp5;
    int nfiles_tmp5, iunit5_tmp5, ntr_tmp5, ncyctot_tmp5;
    int lfgm_tmp5, lfem_tmp9, lcgm_tmp9;
    // local character arrays
    char grid[81], plt3dg[81], plt3dq[81], output[81], residual[81];
    char turbres[81], blomx[81], output2[81], printout[81], pplunge[81];
    char ovrlap[81], patch[81], restrt[81], subres[81], subtur[81];
    char grdmov[81], alphahist[81], errfile[81];

    // ---- executable body ----
    pi = 4.0f * std::atan(1.0f);
    itoti  = 0;
    itotb  = 0;
    nfiles = 1;

    // set error flag to -99 if this routine is called during array sizing step
    if (icall == 0) {
        ierrflg = -99;
    } else {
        ierrflg = -1;
    }

    // read keyword-driven input, if any
    readkey_ns::readkey(ititr, myid, ibufdim, nbuf, bou, nou, iunit11, ierrflg);

    // free-form alphanumeric title
    if (ititr == 0) {
        // read(iunit5,10)(realval(i),i=1,20)  format 10 = 20a4
        // This is a FORMATTED character read of exactly ONE line (20 fields of
        // 4 chars = 80 chars), NOT a list-directed read.  Reading it list-
        // directed would gobble following lines until 20 tokens were collected,
        // desynchronising every subsequent read.  Read one line and pack it as
        // 20 four-character groups into the title (float-aliased char storage).
        char _titleline[210];
        std::memset(_titleline, ' ', sizeof(_titleline));
        if (fgets(_titleline, 205, fortran_get_unit(iunit5)) != nullptr) {
            int _len = (int)std::strlen(_titleline);
            if (_len > 0 && _titleline[_len-1] == '\n') { _titleline[_len-1] = ' '; _len--; }
            for (int ii = _len; ii < 200; ii++) _titleline[ii] = ' ';
        }
        for (int i = 1; i <= 20; i++) {
            std::memcpy(&title_arr[i-1], &_titleline[(i-1)*4], 4);
        }
    }
    // write(iunit11,11)(real(title(i)),i=1,20)  format 11 = /1h ,20a4
    {
        FILE* f11 = fortran_get_unit(iunit11);
        fprintf(f11, "\n ");
        for (int i = 1; i <= 20; i++) {
            // title is stored as float (4-byte), interpret as 4 chars
            char buf[5];
            std::memcpy(buf, &title_arr[i-1], 4);
            buf[4] = '\0';
            fprintf(f11, "%.4s", buf);
        }
        fprintf(f11, "\n");
    }

    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    // read(iunit5,*) (realval(i),i=1,5),ialph,ihstry
    fortran_read_list(iunit5, &realval[1], &realval[2], &realval[3],
                      &realval[4], &realval[5], &ialph, &ihstry);
    xmach = realval[1];
    alpha = realval[2];
    beta  = realval[3];
    reue  = realval[4];
    tinf  = realval[5];
    // write(iunit11,24)
    fprintf(fortran_get_unit(iunit11),
        "      Mach     alpha      beta      ReUe   Tinf,dR     ialph    ihstry\n");

    if ((float)reue <= 0.0f) reue = 1.0f;
    reue = reue * 1.0e+06f;
    if ((float)tinf <= 0.0f) tinf = 460.0f;

    // write(iunit11,20) real(xmach),real(alpha),real(beta),real(reue),real(tinf),ialph,ihstry
    // format 20: 3f10.5,e10.3,f10.5,i10,2i10
    fprintf(fortran_get_unit(iunit11), "%10.5f%10.5f%10.5f%10.3e%10.5f%10d%10d\n",
            (float)xmach, (float)alpha, (float)beta, (float)reue, (float)tinf, ialph, ihstry);
    if (cltarg != 99999.0f) {
        fprintf(fortran_get_unit(iunit11),
            " Note: will adjust alpha to reach cltarg =%12.7f\n", (float)cltarg);
    }

    alpha = alpha / radtodeg;
    beta  = beta  / radtodeg;

    // write(iunit11,16)
    fprintf(fortran_get_unit(iunit11),
        "\n note: isnd and c2spe no longer used -see notes on bc2004\n\n");

    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    // read(iunit5,*) (realval(i),i=1,6)
    fortran_read_list(iunit5, &realval[1], &realval[2], &realval[3],
                      &realval[4], &realval[5], &realval[6]);
    sref = realval[1];
    cref = realval[2];
    bref = realval[3];
    xmc  = realval[4];
    ymc  = realval[5];
    zmc  = realval[6];

    if ((float)sref == 0.0f) sref = 1.0f;
    if ((float)cref <= 0.0f) cref = 1.0f;
    if ((float)bref <= 0.0f) bref = 1.0f;

    // write(iunit11,23)
    fprintf(fortran_get_unit(iunit11),
        "      sref      cref      bref       xmc       ymc       zmc\n");
    // write(iunit11,22) real(sref),real(cref),real(bref),real(xmc),real(ymc),real(zmc)
    // format 22: f10.3,5f10.5
    fprintf(fortran_get_unit(iunit11), "%10.3f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
            (float)sref, (float)cref, (float)bref, (float)xmc, (float)ymc, (float)zmc);
    if ((float)sref < 0.0f) {
        fprintf(fortran_get_unit(iunit11), "WARNING: input sref<0, taking abs(sref)\n");
        fprintf(fortran_get_unit(iunit11),
            "  Note: fixed Cl option is no longer triggered with sref<0...use keyword input\n");
        { double _tmp = sref; sref = (float)ccomplex_ns::ccabs(_tmp); }
    }

    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    // read(iunit5,*) realval(1),irest,iflagts,realval(2),iunst,realval(3)
    fortran_read_list(iunit5, &realval[1], &irest, &iflagts,
                      &realval[2], &iunst, &realval[3]);
    dt     = realval[1];
    dtold  = dt;
    fmax   = realval[2];
    cfltau = realval[3];
    if ((float)dt < 0.0f) iunst = 0;
    if ((float)dt < 0.0f && ides >= 1) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: DES/DDES SHOULD be run time-accurately (dt>0)\n");
    }
    if ((float)dt > 0.0f && cltarg != 99999.0f) {
        fprintf(fortran_get_unit(iunit11), " Cannot use cltarg with time-accurate\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    // write(iunit11,31)
    fprintf(fortran_get_unit(iunit11),
        "        dt     irest   iflagts      fmax     iunst   cfl_tau\n");
    // write(iunit11,30) real(dt),irest,iflagts,real(fmax),iunst,real(cfltau)
    // format 30: f10.5,2i10,f10.5,i10,f10.5
    fprintf(fortran_get_unit(iunit11), "%10.5f%10d%10d%10.5f%10d%10.5f\n",
            (float)dt, irest, iflagts, (float)fmax, iunst, (float)cfltau);

    if (irest < 0) {
        irest2 = 1;
        irest  = -irest;
    } else {
        irest2 = 0;
    }
    if (irest == 0 && iteravg == 2) {
        fprintf(fortran_get_unit(iunit11),
            "\n when irest=0, iteravg=2 starts running-average from scratch (as if iteravg=1)\n\n");
    }

    { double _tmp = fmax; fmax = (float)ccomplex_ns::ccabs(_tmp); }


    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    // read(iunit5,*) ngrid,nplot3d,nprint,nwrest,ichk,i2d,ntstep,ita
    fortran_read_list(iunit5, &ngrid, &nplot3d, &nprint, &nwrest, &ichk, &i2d, &ntstep, &ita);
    // write(iunit11,1639)
    fprintf(fortran_get_unit(iunit11),
        "     ngrid   nplot3d    nprint    nwrest      ichk       i2d    ntstep       ita\n");
    if ((float)dt < 0.0f) ita = 1;
    if (ita == 0 && (float)dt > 0.0f) {
        fprintf(fortran_get_unit(iunit11),
            " stopping...must chose abs(ita) > 0 for time-accurate computations\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (i2d < -1 || i2d > 1) i2d = 0;
    if ((float)dt < 0.0f) ntstep = 1;
    if (nwrest <= 0) {
        fprintf(fortran_get_unit(iunit11),
            " NWREST set incorrectly.  Resetting to 2000\n");
        nwrest = 2000;
    }
    if (nwrest < 20) {
        fprintf(fortran_get_unit(iunit11),
            " NWREST is very small... be sure you want this\n");
    }
    // write(iunit11,36) ngrid,nplot3d,nprint,nwrest,ichk,i2d,ntstep,ita
    // format 36: 8i10
    fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%10d%10d\n",
            ngrid, nplot3d, nprint, nwrest, ichk, i2d, ntstep, ita);
    iipv = 0;
    if (i2d == -1) {
        iipv = 1;
        i2d  = 1;
    }
    // DES must be run in 3-D
    if (ides >= 1 && i2d == 1) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: DES/DDES SHOULD be run 3-D (i2d=0)\n");
    }
    // Check for SARC & EASMCC models
    if ((isarc2d == 1 || ieasmcc2d == 1) && i2d == 0) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: SARC & EASMCC currently only coded for 2-D curvature in x-z or x-y planes;\n");
        fprintf(fortran_get_unit(iunit11),
            "          (i-index MUST be in spanwise direction; i-index derivatives are not accounted for)...\n");
        fprintf(fortran_get_unit(iunit11),
            "          usage in 3-D cases is discouraged\n");
    }
    if (isarc2d == 1 && isar == 1) {
        fprintf(fortran_get_unit(iunit11), " ERROR: cannot have isarc2d and isar=1\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (isarc3d == 1 && isar == 1) {
        fprintf(fortran_get_unit(iunit11), " ERROR: cannot have isarc3d and isar=1\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (isarc3d == 1 && isarc2d == 1) {
        fprintf(fortran_get_unit(iunit11), " ERROR: cannot have isarc3d and isarc2d=1\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    // for cfl3d type grids, set ialph = 0 regardless of input value
    if (ngrid > 0 && ialph != 0) {
        ialph = 0;
        fprintf(fortran_get_unit(iunit11),
            "\n WARNING: IALPH currently hardwired to 0 for NGRID>0 - type grids.  See manual.\n");
    }

    ip3dgrd = 0;
    if (ngrid < 0) {
        ip3dgrd = 1;
        ngrid   = std::abs(ngrid);
    }

    if (ialph == 0) {
        uinf = xmach * std::cos(alpha) * std::cos(beta);
        winf = xmach * std::sin(alpha) * std::cos(beta);
        vinf = -xmach * std::sin(beta);
    } else {
        uinf = xmach * std::cos(alpha) * std::cos(beta);
        winf = xmach * std::sin(beta);
        vinf = xmach * std::sin(alpha) * std::cos(beta);
        temp = ymc;
        ymc  = -zmc;
        zmc  = temp;
    }
    // write(iunit11,18)
    fprintf(fortran_get_unit(iunit11),
        "\n the input values of Mach, alpha, beta and ialph yield a freestream velocity\n"
        " with x,y,z components:  u0 = %8.4f,   v0 = %8.4f,   w0 = %8.4f\n"
        "           ** make sure this is consistent with your grid **\n\n",
        (float)uinf, (float)vinf, (float)winf);

    nchk = maxgr - ngrid;

    // check maximum number of grids
    if (nchk < 0) {
        fprintf(fortran_get_unit(iunit11),
            " stopping - insufficient maximum number of grids(maxgr)\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    // check maximum number of output sets
    if (nplot3d > nplots || nprint > nplots) {
        fprintf(fortran_get_unit(iunit11), " nplots too small\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl    = 0;
    nblock = 0;
    // write(iunit11,61)
    fprintf(fortran_get_unit(iunit11),
        "       ncg       iem  iadvance    iforce  ivisc(i)  ivisc(j)  ivisc(k)\n");
    iemtot = 0;
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        // read(iunit5,*) ncg,iem,iad,ifor,ivisc
        // ivisc is the common block array ivisc[3] (1-based)
        fortran_read_list(iunit5, &ncg, &iem, &iad, &ifor,
                          &ivisc[0], &ivisc[1], &ivisc[2]);
        // write(iunit11,36) ncg,iem,iad,ifor,ivisc
        fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%10d\n",
                ncg, iem, iad, ifor, ivisc[0], ivisc[1], ivisc[2]);
        nblock = nblock + ncg + 1;
        if (iem > 0 && ncg > 0) {
            fprintf(fortran_get_unit(iunit11),
                " Embedded grid must have ncg=0 (it uses grid\n");
            fprintf(fortran_get_unit(iunit11),
                " in which it is embedded for coarser levels)\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (ncg > 4) {
            fprintf(fortran_get_unit(iunit11), " known bug: ncg must not exceed 4\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }


        if (ivisc[0] == -2 || ivisc[1] == -2 || ivisc[2] == -2) {
            fprintf(fortran_get_unit(iunit11),
                " NOTICE:  ivisc=-2 is no longer B-L with Degani-Schiff\n");
            fprintf(fortran_get_unit(iunit11),
                "     ivisc < 0 now indicates Wall Function\n");
            fprintf(fortran_get_unit(iunit11),
                "     Use ivisc=3 for B-L with Degani-Schiff\n");
        }

        if (ivisc[0] == -16 || ivisc[1] == -16 || ivisc[2] == -16 ||
            ivisc[0] == -30 || ivisc[1] == -30 || ivisc[2] == -30 ||
            ivisc[0] == -40 || ivisc[1] == -40 || ivisc[2] == -40) {
            fprintf(fortran_get_unit(iunit11),
                " cannot use wall fns with ivisc 16,30,40\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        // Set wall iwf, wall function parameter:
        iwf[0] = 0; iwf[1] = 0; iwf[2] = 0;
        if (ivisc[0] < 0) { iwf[0] = 1; ivisc[0] = std::abs(ivisc[0]); }
        if (ivisc[1] < 0) { iwf[1] = 1; ivisc[1] = std::abs(ivisc[1]); }
        if (ivisc[2] < 0) { iwf[2] = 1; ivisc[2] = std::abs(ivisc[2]); }

        // DES currently only works with models 5,6,7
        if (ides == 1 || ides == 2 || ides == 3) {
            if (i2d != 1) {
                if (!((ivisc[0]==5||ivisc[0]==6||ivisc[0]==7||ivisc[0]==40) &&
                      (ivisc[1]==5||ivisc[1]==6||ivisc[1]==7||ivisc[1]==40) &&
                      (ivisc[2]==5||ivisc[2]==6||ivisc[2]==7||ivisc[2]==40))) {
                    fprintf(fortran_get_unit(iunit11),
                        " ides=1,2,3 currently works with ivisc 5,6,7,40 only, and must be on in all directions\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            } else {
                if (!((ivisc[1]==5||ivisc[1]==6||ivisc[1]==7||ivisc[1]==40) &&
                      (ivisc[2]==5||ivisc[2]==6||ivisc[2]==7||ivisc[2]==40))) {
                    fprintf(fortran_get_unit(iunit11),
                        " ides=1,2,3 currently works with ivisc 5,6,7,40 only, and must be on in all directions\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }
        }
        if (ides > 3) {
            if (i2d != 1) {
                if ((ivisc[0]<5||ivisc[0]>5) && (ivisc[1]<5||ivisc[1]>5) && (ivisc[2]<5||ivisc[2]>5)) {
                    fprintf(fortran_get_unit(iunit11),
                        " ides>3 currently works with ivisc 5 only, and must be on in all directions\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            } else {
                if ((ivisc[1]<5||ivisc[1]>5) && (ivisc[2]<5||ivisc[2]>5)) {
                    fprintf(fortran_get_unit(iunit11),
                        " ides>3 currently works with ivisc 5 only, and must be on in all directions\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }
        }
        // for Degani schiff, set ideg=1 and reset ivisc=2 within code:
        ideg[0] = 0; ideg[1] = 0; ideg[2] = 0;
        if (ivisc[0] == 3) { ideg[0] = 1; ivisc[0] = 2; }
        if (ivisc[1] == 3) { ideg[1] = 1; ivisc[1] = 2; }
        if (ivisc[2] == 3) { ideg[2] = 1; ivisc[2] = 2; }
        iv(igrid) = std::max(ivisc[0], ivisc[1]);
        iv(igrid) = std::max(ivisc[2], iv(igrid));
        if (ivisc[0]>16 || ivisc[1]>16 || ivisc[2]>16) {
            if (ivisc[0]!=25 && ivisc[1]!=25 && ivisc[2]!=25 &&
                ivisc[0]!=30 && ivisc[1]!=30 && ivisc[2]!=30 &&
                ivisc[0]!=40 && ivisc[1]!=40 && ivisc[2]!=40 &&
                ivisc[0]!=72 && ivisc[1]!=72 && ivisc[2]!=72) {
                fprintf(fortran_get_unit(iunit11), " Unknown turb model type chosen\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        }
        if (ivisc[0]>=4 || ivisc[1]>=4 || ivisc[2]>=4) {
            isminc = 1;
        } else {
            isminc = 0;
        }
        if (ivisc[0]>1 && ivisc[1]>1 && ivisc[0]!=ivisc[1]) {
            fprintf(fortran_get_unit(iunit11), " cannot mix turbulence model types\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (ivisc[0]>1 && ivisc[2]>1 && ivisc[0]!=ivisc[2]) {
            fprintf(fortran_get_unit(iunit11), " cannot mix turbulence model types\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (ivisc[1]>1 && ivisc[2]>1 && ivisc[1]!=ivisc[2]) {
            fprintf(fortran_get_unit(iunit11), " cannot mix turbulence model types\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        if (ivisc[1] == 2 && ivisc[0] == 2) {
            fprintf(fortran_get_unit(iunit11),
                " B-L progr. only for turbulent flow in two directions either K-J or K-I directions only\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        iemtot         = iemtot + iem;
        ncgg(igrid)    = ncg;
        if (igrid == 1) {
            ncgmax = ncg;
        } else {
            ncgmax = std::max(ncgmax, ncg);
        }

        iemg(igrid)    = iem;
        iadvance(nbl)  = iad;
        iforce(nbl)    = ifor;
        iviscg(nbl, 1) = ivisc[0];
        iviscg(nbl, 2) = ivisc[1];
        iviscg(nbl, 3) = ivisc[2];
        idegg(nbl, 1)  = ideg[0];
        idegg(nbl, 2)  = ideg[1];
        idegg(nbl, 3)  = ideg[2];
        iwfg(nbl, 1)   = iwf[0];
        iwfg(nbl, 2)   = iwf[1];
        iwfg(nbl, 3)   = iwf[2];
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl              = nbl + 1;
                iadvance(nbl)    = iadvance(nbl-1);
                iforce(nbl)      = iforce(nbl-1);
                iviscg(nbl, 1)   = ivisc[0];
                iviscg(nbl, 2)   = ivisc[1];
                iviscg(nbl, 3)   = ivisc[2];
                idegg(nbl, 1)    = ideg[0];
                idegg(nbl, 2)    = ideg[1];
                idegg(nbl, 3)    = ideg[2];
                iwfg(nbl, 1)     = iwf[0];
                iwfg(nbl, 2)     = iwf[1];
                iwfg(nbl, 3)     = iwf[2];
            }
        }
    } // end do 71 igrid=1,ngrid


    // check maximum number of blocks
    nchk = maxbl - nblock;
    if (nchk < 0) {
        fprintf(fortran_get_unit(iunit11),
            " stopping - insufficient maximum number of blocks(maxbl)\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    ivmx = iv(1);
    if (ngrid > 1) {
        for (igrid = 2; igrid <= ngrid; igrid++) {
            if (iv(igrid) > ivmx) ivmx = iv(igrid);
        }
        if (ivmx > 3) {
            for (igrid = 2; igrid <= ngrid; igrid++) {
                if (iv(igrid) != iv(igrid-1)) {
                    fprintf(fortran_get_unit(iunit11),
                        " for turb models > 3, need same number on ALL blocks\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }
        }
    }

    if ((isarc2d == 1 || isarc3d == 1) && ivmx != 5) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: isarc2d/3d not active unless SA model (ivisc=5) used\n");
    }
    if (ieasmcc2d == 1 && (ivmx!=8 && ivmx!=9 && ivmx!=11 && ivmx!=12 && ivmx!=13 && ivmx!=14)) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: ieasmcc2d not active unless EASM model (ivisc=8,9,11,12,13, or 14) used\n");
    }

    if (i_wilcox06 == 1 && ivmx != 6) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: i_wilcox06 not active unless Wilcox k-o model (ivisc=6) used\n");
    }
    if (i_wilcox98 == 1 && ivmx != 6) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: i_wilcox98 not active unless Wilcox k-o model (ivisc=6) used\n");
    }
    if (i_wilcox98 == 1 && i_wilcox06 == 1) {
        fprintf(fortran_get_unit(iunit11),
            " Error: cannot have both i_wilcox98=1 and i_wilcox06=1\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    if (i_turbprod_kterm == 1 && (ivmx != 6 && ivmx != 7)) {
        if (ikoprod != 1) {
            fprintf(fortran_get_unit(iunit11),
                " WARNING: i_turbprod_kterm not active unless ivisc=6 or 7 used AND ikoprod=1\n");
        }
    }

    if (i_catris_kw == 1 && (ivmx!=6 && ivmx!=7 && ivmx!=8 && ivmx!=12 && ivmx!=14)) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: i_catris_kw not active unless ivisc=6,7,8,12, or 14 used\n");
    }

    if (isstrc > 0 && (ivmx != 6 && ivmx != 7)) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: isstrc not active unless ivisc=6 or 7\n");
    }
    if (isstsf == 1 && (ivmx != 6 && ivmx != 7)) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: isstsf not active unless ivisc=6 or 7\n");
    }

    if (isst2003 == 1 && ivmx != 7) {
        fprintf(fortran_get_unit(iunit11), " ERROR: isst2003 must have ivisc=7\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    if (std::abs(prod2d3dtrace - 0.5f) < 0.01f) {
        prod2d3dtrace = 0.5f;
    } else if (std::abs(prod2d3dtrace - 0.33333333f) < 0.01f) {
        prod2d3dtrace = 0.33333333f;
    } else {
        prod2d3dtrace = 0.0f;
    }

    if (ikoprod > 0 && (ivmx!=6 && ivmx!=7 && ivmx!=10 && ivmx!=15 && ivmx!=16 && ivmx!=30 && ivmx!=40)) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: ikoprod not active unless ivisc=6,7,10,15,16,30, or 40 used\n");
    }

    if (iturbprod > 0 && (ivmx!=8 && ivmx!=9 && ivmx!=13 && ivmx!=14)) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: iturbprod not active unless ivisc=8,9,13, or 14 used\n");
    }

    if (isstdenom != 0 && ivmx != 7) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: isstdenom not active unless ivisc=7 used\n");
    } else if (isstdenom != 0 && i_bsl != 0) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: isstdenom not active when i_bsl=1 (BSL option in effect)\n");
    }

    if (ibeta8kzeta != 0 && ivmx != 15) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: ibeta8kzeta not active unless ivisc=15 used\n");
    }

    if (i_bsl != 0 && ivmx != 7) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: i_bsl not active unless ivisc=7 used\n");
    }

    if (keepambient != 0 && (ivmx!=6 && ivmx!=7 && ivmx!=8 && ivmx!=9 && ivmx!=10 &&
        ivmx!=11 && ivmx!=12 && ivmx!=13 && ivmx!=14 && ivmx!=30 && ivmx!=40)) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: keepambient not active unless ivisc=6,7,8,9,10,11,12,13,14,30, or 40 used\n");
    }

    if (ivmx==25 || ivmx==30 || ivmx==40 || ivmx==72 || ivmx==16) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: ivisc=%4d still under development... use at your own risk!\n", ivmx);
    }

    // Need to check for keywords that do not work for ivmx=72:
    if (ivmx == 72) {
        fprintf(fortran_get_unit(iunit11),
            " WARNING: keyword ifullns has no effect for ivisc=72 - it is ALWAYS full NS\n");
        fprintf(fortran_get_unit(iunit11),
            "          (but diffusion terms in turb eqns are still thin-layer)\n");
        if (ifullns != 0) {
            fprintf(fortran_get_unit(iunit11), " ERROR: ifullns must be 0 for ivisc=72\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (edvislim != 1.0e10f) {
            fprintf(fortran_get_unit(iunit11),
                " WARNING: keyword edvislim has influence only in turb heat flux in energy eqn for ivisc=72\n");
        }
        if (icgns != 0) {
            fprintf(fortran_get_unit(iunit11), " ERROR: CGNS not working yet for ivisc=72\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (itaturb != 1) {
            fprintf(fortran_get_unit(iunit11), " ERROR: itaturb not implemented for ivisc=72\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (iturbord != 1) {
            fprintf(fortran_get_unit(iunit11),
                " ERROR: iturbord not implemented for ivisc=72 (default is 1st order)\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (istrongturbdis != 0) {
            fprintf(fortran_get_unit(iunit11),
                " ERROR: istrongturbdis not implemented for ivisc=72\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (tur1cut != 1.0e-20f || tur2cut != 1.0e-20f) {
            fprintf(fortran_get_unit(iunit11),
                " ERROR: tur1cut, tur2cut not implemented for ivisc=72\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (tur1cutlev != 0.0f || tur2cutlev != 0.0f) {
            fprintf(fortran_get_unit(iunit11),
                " ERROR: tur1cutlev, tur2cutlev not implemented for ivisc=72\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (i_lam_forcezero != 0) {
            fprintf(fortran_get_unit(iunit11),
                " ERROR: i_lam_forcezero not implemented for ivisc=72\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }


    if (i_nonlin != 0 && (ivmx < 2 || ivmx==8 || ivmx==9 || ivmx==11 || ivmx==12 ||
        ivmx==13 || ivmx==14 || ivmx > 15)) {
        fprintf(fortran_get_unit(iunit11), " i_nonlin cannot be used with this ivisc\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    if (i_tauijs == 1) {
        if (ivmx < 70 && ngrid != 1) {
            fprintf(fortran_get_unit(iunit11),
                " ERROR: across-boundary BCs not correct for i_tauijs=1, ivmx<70; cannot use multiple zones currently\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (ifullns != 0) {
            fprintf(fortran_get_unit(iunit11), " ERROR: ifullns must be 0 when i_tauijs=1\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (ivmx==11 || ivmx==12 || ivmx==13 || ivmx==14) {
            fprintf(fortran_get_unit(iunit11),
                " i_tauijs=1 cannot be used with ivisc=11,12,13,14\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (i_nonlin != 0) {
            fprintf(fortran_get_unit(iunit11), " i_tauijs=1 cannot be used with i_nonlin\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }

    if (i_qcr2000 == 1 && i_tauijs != 1) {
        fprintf(fortran_get_unit(iunit11), " i_qcr2000 must be used with i_tauijs=1\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (i_qcr2013 == 1 && i_tauijs != 1) {
        fprintf(fortran_get_unit(iunit11), " i_qcr2013 must be used with i_tauijs=1\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (i_qcr2013v == 1 && i_tauijs != 1) {
        fprintf(fortran_get_unit(iunit11), " i_qcr2013v must be used with i_tauijs=1\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (i_qcr2000 == 1 && i_qcr2013 == 1) {
        fprintf(fortran_get_unit(iunit11), " cannot have both i_qcr2000=1 and i_qcr2013=1\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (i_qcr2000 == 1 && i_qcr2013v == 1) {
        fprintf(fortran_get_unit(iunit11), " cannot have both i_qcr2000=1 and i_qcr2013v=1\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (i_qcr2013 == 1 && i_qcr2013v == 1) {
        fprintf(fortran_get_unit(iunit11), " cannot have both i_qcr2013=1 and i_qcr2013v=1\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (i_qcr2013 == 1 && ivmx != 5) {
        fprintf(fortran_get_unit(iunit11), " i_qcr2013 is only usable with SA model\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (i_qcr2013v == 1 && ivmx != 5) {
        fprintf(fortran_get_unit(iunit11), " i_qcr2013v is only usable with SA model\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    if ((ivmx==8 || ivmx==9 || ivmx==13 || ivmx==14) && (float)cmulim > 0.04f) {
        fprintf(fortran_get_unit(iunit11), " error... cmulim set too high\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    // grid dimensions
    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    // write(iunit11,1631)
    fprintf(fortran_get_unit(iunit11), "      idim      jdim      kdim\n");
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        // read(iunit5,*) idim,jdim,kdim
        fortran_read_list(iunit5, &idim, &jdim, &kdim);
        // write(iunit11,336) idim,jdim,kdim
        fprintf(fortran_get_unit(iunit11), "%10d%10d%10d\n", idim, jdim, kdim);
        if (i2d > 0 && idim != 2) {
            fprintf(fortran_get_unit(iunit11),
                " stopping...must have idim=2 for 2d cases!\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (i2d == 0 && idim == 2) {
            fprintf(fortran_get_unit(iunit11),
                " WARNING: idim=2 for i2d=0 may not be a good idea\n");
        }
        if (iaxi2plane == 1 && (i2d != 0 || idim != 2)) {
            fprintf(fortran_get_unit(iunit11), " WARNING: iaxi2plane set back to 0\n");
            fprintf(fortran_get_unit(iunit11), "    (i2d must=0 and idim must=2)\n");
            iaxi2plane = 0;
        }
        if (iaxi2planeturb == 1 && (i2d != 0 || idim != 2)) {
            fprintf(fortran_get_unit(iunit11), " WARNING: iaxi2planeturb set back to 0\n");
            fprintf(fortran_get_unit(iunit11), "    (i2d must=0 and idim must=2)\n");
            iaxi2planeturb = 0;
        }
        // check usage of full N-S
        if (ifullns != 0) {
            fprintf(fortran_get_unit(iunit11),
                " Note that even with ifullns=1, diffusion terms in turb eqns are still thin-layer)\n");
            if (i2d == 0 && idim > 2 && (ivisc[0]==0 || ivisc[1]==0 || ivisc[2]==0)) {
                fprintf(fortran_get_unit(iunit11),
                    " WARNING: When using full N-S in 3-D, viscous terms (ivisc) usually \n");
                fprintf(fortran_get_unit(iunit11), "must be ON in all 3 directions\n");
            }
            if (i2d == 0 && idim == 2 && (ivisc[1]==0 || ivisc[2]==0)) {
                fprintf(fortran_get_unit(iunit11),
                    " ERROR: When using full N-S in 3-D, with idim=2, viscous terms (ivisc) must be\n");
                fprintf(fortran_get_unit(iunit11), " ON in both J and K directions\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (i2d != 0 && (ivisc[1]==0 || ivisc[2]==0)) {
                fprintf(fortran_get_unit(iunit11),
                    " ERROR: When using full N-S in 2-D, viscous terms (ivisc) must be\n");
                fprintf(fortran_get_unit(iunit11), " ON in both J and K directions\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (idim == 2 && ivisc[0] != 0) {
                fprintf(fortran_get_unit(iunit11),
                    " ERROR: Having ifullns .ne. 0 with idim=2 and ivisc(1) .ne. 0\n");
                fprintf(fortran_get_unit(iunit11),
                    " can cause major problems; recommend setting ivisc(1)=0\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        }
        if (ides >= 1 && idim == 2) {
            fprintf(fortran_get_unit(iunit11),
                " WARNING: idim=2 for DES/DDES may not be a good idea\n");
        }
        nblg(igrid) = nbl;
        idimg(nbl)  = idim;
        jdimg(nbl)  = jdim;
        kdimg(nbl)  = kdim;
        if (jdim == 2 || kdim == 2) {
            fprintf(fortran_get_unit(iunit11), " Error.  Neither JDIM nor KDIM can be 2.\n");
            fprintf(fortran_get_unit(iunit11),
                " (can cause problems when need 2 ghost BCs)\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl        = nbl + 1;
                idimg(nbl) = idimg(nbl-1)/2 + 1;
                jdimg(nbl) = jdimg(nbl-1)/2 + 1;
                kdimg(nbl) = kdimg(nbl-1)/2 + 1;
                if (idim == 2) {
                    idimg(nbl) = 2;
                }
                if (jdimg(nbl) == 2 || kdimg(nbl) == 2) {
                    fprintf(fortran_get_unit(iunit11),
                        " Error.  Coarser levels of JDIM and KDIM cannot be 2.\n");
                    fprintf(fortran_get_unit(iunit11),
                        " (can cause problems when need 2 ghost BCs)\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                istop = 0;
                if ((float)(idimg(nbl-1)/2) == (float)idimg(nbl-1)/2.0f && idim > 2) {
                    fprintf(fortran_get_unit(iunit11),
                        " Cannot create coarser level for idim past%6d\n", idimg(nbl-1));
                    istop = 1;
                }
                if ((float)(jdimg(nbl-1)/2) == (float)jdimg(nbl-1)/2.0f) {
                    fprintf(fortran_get_unit(iunit11),
                        " Cannot create coarser level for jdim past%6d\n", jdimg(nbl-1));
                    istop = 1;
                }
                if ((float)(kdimg(nbl-1)/2) == (float)kdimg(nbl-1)/2.0f) {
                    fprintf(fortran_get_unit(iunit11),
                        " Cannot create coarser level for kdim past%6d\n", kdimg(nbl-1));
                    istop = 1;
                }
                if (istop == 1) {
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }
        }
    } // end do 7001 igrid=1,ngrid


    // laminar region indices
    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    // write(iunit11,1939)
    fprintf(fortran_get_unit(iunit11),
        "    ilamlo    ilamhi    jlamlo    jlamhi    klamlo    klamhi\n");
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        // read(iunit5,*) ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi
        fortran_read_list(iunit5, &ilamlo, &ilamhi, &jlamlo, &jlamhi, &klamlo, &klamhi);
        if (ilamhi > idimg(nbl)) ilamhi = idimg(nbl);
        if (jlamhi > jdimg(nbl)) jlamhi = jdimg(nbl);
        if (klamhi > kdimg(nbl)) klamhi = kdimg(nbl);
        if (ilamlo <= 0) { ilamlo = 0; ilamhi = 0; }
        if (ilamhi <= 0) { ilamlo = 0; ilamhi = 0; }
        if (jlamlo <= 0) { jlamlo = 0; jlamhi = 0; }
        if (jlamhi <= 0) { jlamlo = 0; jlamhi = 0; }
        if (klamlo <= 0) { klamlo = 0; klamhi = 0; }
        if (klamhi <= 0) { klamlo = 0; klamhi = 0; }
        // write(iunit11,36) ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi
        fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d\n",
                ilamlo, ilamhi, jlamlo, jlamhi, klamlo, klamhi);
        ilamlog(nbl) = ilamlo;
        ilamhig(nbl) = ilamhi;
        jlamlog(nbl) = jlamlo;
        jlamhig(nbl) = jlamhi;
        klamlog(nbl) = klamlo;
        klamhig(nbl) = klamhi;

        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                ilamlog(nbl) = ilamlog(nbl-1)/2 + 1;
                ilamhig(nbl) = ilamhig(nbl-1)/2 + 1;
                jlamlog(nbl) = jlamlog(nbl-1)/2 + 1;
                jlamhig(nbl) = jlamhig(nbl-1)/2 + 1;
                klamlog(nbl) = klamlog(nbl-1)/2 + 1;
                klamhig(nbl) = klamhig(nbl-1)/2 + 1;
                if (i2d == 1) {
                    ilamlog(nbl) = ilamlog(nbl-1);
                    ilamhig(nbl) = ilamhig(nbl-1);
                }
                if (ilamlog(nbl-1) == 0) { ilamlog(nbl) = 0; ilamhig(nbl) = 0; }
                if (jlamlog(nbl-1) == 0) { jlamlog(nbl) = 0; jlamhig(nbl) = 0; }
                if (klamlog(nbl-1) == 0) { klamlog(nbl) = 0; klamhig(nbl) = 0; }
            }
        }
    } // end do 7201

    // inewg, igridc, is, js, ks, ie, je, ke
    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    // write(iunit11,1931)
    fprintf(fortran_get_unit(iunit11),
        "     inewg    igridc        is        js        ks        ie        je        ke\n");
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        // read(iunit5,*) inewg,igridc,is,js,ks,ie,je,ke
        fortran_read_list(iunit5, &inewg, &igridc, &is, &js, &ks, &ie, &je, &ke);
        // write(iunit11,36) inewg,igridc,is,js,ks,ie,je,ke
        fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                inewg, igridc, is, js, ks, ie, je, ke);
        if (igridc > 0) {
            nblcg(nbl) = nblg(igridc);
        } else {
            nblcg(nbl) = nbl;
        }
        inewgg(igrid) = inewg;
        jsg(nbl) = js;
        ksg(nbl) = ks;
        isg(nbl) = is;
        jeg(nbl) = je;
        keg(nbl) = ke;
        ieg(nbl) = ie;

        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl        = nbl + 1;
                nblcg(nbl) = nblcg(nbl-1) + 1;
                jsg(nbl)   = jsg(nbl-1)/2 + 1;
                ksg(nbl)   = ksg(nbl-1)/2 + 1;
                isg(nbl)   = isg(nbl-1)/2 + 1;
                jeg(nbl)   = jeg(nbl-1)/2 + 1;
                keg(nbl)   = keg(nbl-1)/2 + 1;
                ieg(nbl)   = ieg(nbl-1)/2 + 1;
            }
        }
    } // end do 7002

    // idiag, iflim
    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    // write(iunit11,9004)
    fprintf(fortran_get_unit(iunit11),
        "  idiag(i)  idiag(j)  idiag(k)  iflim(i)  iflim(j)  iflim(k)\n");
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        // read(iunit5,*) idiag,iflim
        fortran_read_list(iunit5, &idiag_loc[1], &idiag_loc[2], &idiag_loc[3],
                          &iflim_loc[1], &iflim_loc[2], &iflim_loc[3]);
        // write(iunit11,36) idiag,iflim
        fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d\n",
                idiag_loc[1], idiag_loc[2], idiag_loc[3],
                iflim_loc[1], iflim_loc[2], iflim_loc[3]);
        // limit on ncg when iflim=4
        if ((iflim_loc[1]==4 || iflim_loc[2]==4 || iflim_loc[3]==4) && ncg >= 10) {
            fprintf(fortran_get_unit(iunit11),
                " currently ncg limited < 10 for iflim=4\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        idiagg(nbl, 1) = idiag_loc[1];
        idiagg(nbl, 2) = idiag_loc[2];
        idiagg(nbl, 3) = idiag_loc[3];
        iflimg(nbl, 1) = iflim_loc[1];
        iflimg(nbl, 2) = iflim_loc[2];
        iflimg(nbl, 3) = iflim_loc[3];

        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                idiagg(nbl, 1) = idiag_loc[1];
                idiagg(nbl, 2) = idiag_loc[2];
                idiagg(nbl, 3) = idiag_loc[3];
                iflimg(nbl, 1) = iflim_loc[1];
                iflimg(nbl, 2) = iflim_loc[2];
                iflimg(nbl, 3) = iflim_loc[3];
            }
        }
    } // end do 9002


    // ifds, rkap0
    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    // write(iunit11,1951)
    fprintf(fortran_get_unit(iunit11),
        "   ifds(i)   ifds(j)   ifds(k)  rkap0(i)  rkap0(j)  rkap0(k)\n");
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        // read(iunit5,*) ifds(1),ifds(2),ifds(3),realval(1),realval(2),realval(3)
        fortran_read_list(iunit5, &ifds[0], &ifds[1], &ifds[2],
                          &realval[1], &realval[2], &realval[3]);
        rkap0[0] = realval[1];
        rkap0[1] = realval[2];
        rkap0[2] = realval[3];
        // write(iunit11,1950) ifds(1),ifds(2),ifds(3),real(rkap0(1)),real(rkap0(2)),real(rkap0(3))
        // format 1950: 3i10,3f10.4
        fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10.4f%10.4f%10.4f\n",
                ifds[0], ifds[1], ifds[2], (float)rkap0[0], (float)rkap0[1], (float)rkap0[2]);
        if (iaxi2plane == 1 && (ifds[0]==1 || ifds[1]==1 || ifds[2]==1)) {
            fprintf(fortran_get_unit(iunit11),
                " WARNING: if problems occur when running two-plane axisymmetric case with\n");
            fprintf(fortran_get_unit(iunit11),
                "   singular or near-singular axis, the use of IFDS=0 may help\n");
        }
        ifdsg(nbl, 1)  = ifds[0];
        ifdsg(nbl, 2)  = ifds[1];
        ifdsg(nbl, 3)  = ifds[2];
        rkap0g(nbl, 1) = rkap0[0];
        rkap0g(nbl, 2) = rkap0[1];
        rkap0g(nbl, 3) = rkap0[2];

        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                ifdsg(nbl, 1)  = ifds[0];
                ifdsg(nbl, 2)  = ifds[1];
                ifdsg(nbl, 3)  = ifds[2];
                rkap0g(nbl, 1) = rkap0[0];
                rkap0g(nbl, 2) = rkap0[1];
                rkap0g(nbl, 3) = rkap0[2];
            }
        }
    } // end do 7012

    // grid, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim, iovrlp
    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    // write(iunit11,9009)
    fprintf(fortran_get_unit(iunit11),
        "      grid     nbci0   nbcidim     nbcj0   nbcjdim     nbck0   nbckdim    iovrlp\n");
    msegment = 0;
    lig(1) = 1;
    lbg(1) = 1;
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        // read(iunit5,*) ig,nbci0(nbl),nbcidim(nbl),nbcj0(nbl),nbcjdim(nbl),nbck0(nbl),nbckdim(nbl),iover
        fortran_read_list(iunit5, &ig, &nbci0(nbl), &nbcidim(nbl), &nbcj0(nbl),
                          &nbcjdim(nbl), &nbck0(nbl), &nbckdim(nbl), &iover);
        if (ig != igrid) {
            fprintf(fortran_get_unit(iunit11), " you must put these lines in order by grid!\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        // write(iunit11,36) ig,nbci0(nbl),nbcidim(nbl),nbcj0(nbl),nbcjdim(nbl),nbck0(nbl),nbckdim(nbl),iover
        fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                ig, nbci0(nbl), nbcidim(nbl), nbcj0(nbl),
                nbcjdim(nbl), nbck0(nbl), nbckdim(nbl), iover);
        if (nbci0(nbl)<1 || nbcidim(nbl)<1 || nbcj0(nbl)<1 ||
            nbcjdim(nbl)<1 || nbck0(nbl)<1 || nbckdim(nbl)<1) {
            fprintf(fortran_get_unit(iunit11), " Error... nbci0 etc must be at least 1\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        iovrlp(nbl) = iover;
        // determine size requirements for overset applications
        if (icall == 0 && iovrlp(nbl) == 1) {
            if (nbl > 1) {
                lig(nbl) = lig(nbl - ncgg(igrid));
                lbg(nbl) = lbg(nbl - ncgg(igrid));
            }
            getibk0_ns::getibk0(jdimg(nbl), kdimg(nbl), idimg(nbl), nbl, itotb,
                                 itoti, maxbl, lig, lbg, ibpntsg, iipntsg, nou, bou,
                                 ibufdim, nbuf, ierrflg, myid);
        }
        msegment = std::max(msegment, nbci0(nbl));
        msegment = std::max(msegment, nbcidim(nbl));
        msegment = std::max(msegment, nbcj0(nbl));
        msegment = std::max(msegment, nbcjdim(nbl));
        msegment = std::max(msegment, nbck0(nbl));
        msegment = std::max(msegment, nbckdim(nbl));

        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                iovrlp(nbl) = 0;
                if (iover > 0) iovrlp(nbl) = -1;
                nbci0(nbl)   = nbci0(nbl-1);
                nbcidim(nbl) = nbcidim(nbl-1);
                nbcj0(nbl)   = nbcj0(nbl-1);
                nbcjdim(nbl) = nbcjdim(nbl-1);
                nbck0(nbl)   = nbck0(nbl-1);
                nbckdim(nbl) = nbckdim(nbl-1);
            }
        }
    } // end do 9007

    if (icall == 0) {
        rewind(fortran_get_unit(21));
    }

    if (msegment > maxseg) {
        fprintf(fortran_get_unit(iunit11),
            " increase maxseg parameter; should be: %d\n", msegment);
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }


    // ---- I0 Boundary ----
    // read(iunit5,10) - skip a line
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    fprintf(fortran_get_unit(iunit11),
        "i0:   grid   segment    bctype      jsta      jend      ksta      kend     ndata\n");
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        ifor = iforce(nbl);
        ifo  = ifor / 100;
        jfo  = (ifor - ifo*100) / 10;
        kfo  = (ifor - ifo*100 - jfo*10);
        for (iseg = 1; iseg <= nbci0(nbl); iseg++) {
            ibcinfo(nbl, iseg, 1, 1) = -99;
            ibcinfo(nbl, iseg, 6, 1) = 0;
            if (ifo==1 || ifo==3) ibcinfo(nbl, iseg, 6, 1) = 1;
            bcfilei(nbl, iseg, 1) = 1;
            for (mm = 1; mm <= 12; mm++) bcvali(nbl, iseg, mm, 1) = -1.0e15;
        }
        for (iseg = 1; iseg <= nbci0(nbl); iseg++) {
            { int _v2,_v3,_v4,_v5;
              fortran_read_list(iunit5, &ig, &nseg, &ibctyp, &_v2, &_v3, &_v4, &_v5, &ndata);
              ibcinfo(nbl, std::abs(nseg), 2, 1) = _v2;
              ibcinfo(nbl, std::abs(nseg), 3, 1) = _v3;
              ibcinfo(nbl, std::abs(nseg), 4, 1) = _v4;
              ibcinfo(nbl, std::abs(nseg), 5, 1) = _v5; }
            if (ibctyp == 1004) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...bc1004 no longer available...use 2004 instead\n see bc.f for usage notes\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if ((std::abs(ibctyp)==2004 || std::abs(ibctyp)==2014 ||
                 std::abs(ibctyp)==2024 || std::abs(ibctyp)==2034 ||
                 std::abs(ibctyp)==2016) && iviscg(nbl,1)==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...must have ivisc(i) > 0 if bc2004/14/24/34/16 is being used on an i-face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(ibctyp)==2014 && ivmx<=3) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2014 can only be used for turb models ivisc>3\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(ibctyp)==2024 && ivmx!=30) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2024 can only be used for turb models ivisc=30\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibctyp>=2000 && ibctyp<3000 && ndata==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2000 series bc's require abs(ndata) > 0\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibctyp==9999 && iexact_trunc==0 && iexact_disc==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...9999 bc requires keyword iexact_trunc or iexact_disc\n");
            }
            itemp1 = ibcinfo(nbl, std::abs(nseg), 2, 1);
            itemp2 = ibcinfo(nbl, std::abs(nseg), 3, 1);
            if (itemp1 > itemp2) {
                ibcinfo(nbl, std::abs(nseg), 2, 1) = itemp2;
                ibcinfo(nbl, std::abs(nseg), 3, 1) = itemp1;
            }
            itemp1 = ibcinfo(nbl, std::abs(nseg), 4, 1);
            itemp2 = ibcinfo(nbl, std::abs(nseg), 5, 1);
            if (itemp1 > itemp2) {
                ibcinfo(nbl, std::abs(nseg), 4, 1) = itemp2;
                ibcinfo(nbl, std::abs(nseg), 5, 1) = itemp1;
            }
            ifoseg = 1;
            if (nseg < 0) ifoseg = 0;
            mfoseg = 0;
            if (std::abs(ibctyp)==2004 || std::abs(ibctyp)==1005 ||
                std::abs(ibctyp)==1006 || std::abs(ibctyp)==2014 ||
                std::abs(ibctyp)==2024 || std::abs(ibctyp)==2034 ||
                std::abs(ibctyp)==2016) mfoseg = 1;
            ifoseg = ifoseg * mfoseg;
            nseg = std::abs(nseg);
            if (ndata > 0) {
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                fortran_read_list_array(iunit5, realval+1, ndata);
                for (mm = 1; mm <= ndata; mm++) bcvali(nbl, nseg, mm, 1) = realval[mm];
            }
            if (ndata < 0) {
                nfiles = nfiles + 1;
                if (nfiles > mxbcfil) {
                    fprintf(fortran_get_unit(iunit11),
                        " too many bc files specified...increase parameter mxbcfil\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                bcfilei(nbl, nseg, 1) = nfiles;
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                {
                    char _line[256]; fgets(_line, sizeof(_line), fortran_get_unit(iunit5));
                    int _len = (int)strlen(_line);
                    if (_len > 0 && _line[_len-1] == '\n') _line[--_len] = '\0';
                    std::memset(bcfiles(nfiles), ' ', 80);
                    std::memcpy(bcfiles(nfiles), _line, std::min(_len, 60));
                }
            }
            if (ig != igrid) {
                fprintf(fortran_get_unit(iunit11), " you must put these lines in order by grid!\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibcinfo(nbl, nseg, 1, 1) != -99) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...attempting to set data for segment %d more than once\n", nseg);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            ibcinfo(nbl, nseg, 1, 1) = ibctyp;
            if (ibcinfo(nbl,nseg,2,1)==0 && ibcinfo(nbl,nseg,3,1)==0) {
                ibcinfo(nbl, nseg, 2, 1) = 1; ibcinfo(nbl, nseg, 3, 1) = jdimg(nbl);
            }
            if (ibcinfo(nbl,nseg,4,1)==0 && ibcinfo(nbl,nseg,5,1)==0) {
                ibcinfo(nbl, nseg, 4, 1) = 1; ibcinfo(nbl, nseg, 5, 1) = kdimg(nbl);
            }
            if (ifoseg==0 || ibctyp==0) ibcinfo(nbl, nseg, 6, 1) = 0;
            ibcinfo(nbl, nseg, 7, 1) = ndata;
            fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    ig, nseg, ibcinfo(nbl,nseg,1,1), ibcinfo(nbl,nseg,2,1),
                    ibcinfo(nbl,nseg,3,1), ibcinfo(nbl,nseg,4,1),
                    ibcinfo(nbl,nseg,5,1), ibcinfo(nbl,nseg,7,1));
            if (ndata > 0) {
                { FortranArray1DRef<char[10]> datahdr_ref(datahdr, 10); getdhdr_ns::getdhdr(datahdr_ref, ibctyp, ndata); }
                for (mm = 1; mm <= ndata; mm++) fprintf(fortran_get_unit(iunit11), "%-10.10s", datahdr[mm]);
                fprintf(fortran_get_unit(iunit11), "\n");
                for (mm = 1; mm <= std::min(ndata,4); mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.5f", (float)bcvali(nbl,nseg,mm,1));
                for (mm = 5; mm <= ndata; mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.3e", (float)bcvali(nbl,nseg,mm,1));
                if (ndata > 0) fprintf(fortran_get_unit(iunit11), "\n");
            }
            if (ndata < 0) {
                fprintf(fortran_get_unit(iunit11), "     bc data file\n");
                nfl = bcfilei(nbl, nseg, 1);
                fprintf(fortran_get_unit(iunit11), "     %-60.60s\n", bcfiles(nfl));
            }
            if (ibcinfo(nbl,nseg,2,1)<1 || ibcinfo(nbl,nseg,2,1)>jdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...jsta is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibcinfo(nbl,nseg,3,1)<1 || ibcinfo(nbl,nseg,3,1)>jdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...jend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibcinfo(nbl,nseg,4,1)<1 || ibcinfo(nbl,nseg,4,1)>kdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...ksta is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibcinfo(nbl,nseg,5,1)<1 || ibcinfo(nbl,nseg,5,1)>kdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...kend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        } // end do 1201 iseg


        // set boundary condition type array to 21 for embedded mesh
        if (iemg(igrid) > 0 && isg(nbl) != 1) {
            if (nbci0(nbl) > 1) {
                fprintf(fortran_get_unit(iunit11),
                    " error: embedded mesh boundary at i=1 must extend over entire block face\n");
                fprintf(fortran_get_unit(iunit11), "        segmentation not allowed\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            } else {
                ibcinfo(nbl, 1, 1, 1) = 21;
                ibcinfo(nbl, 1, 2, 1) = 1;
                ibcinfo(nbl, 1, 3, 1) = jdimg(nbl);
                ibcinfo(nbl, 1, 4, 1) = 1;
                ibcinfo(nbl, 1, 5, 1) = kdimg(nbl);
            }
        }
        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
                    ibcinfo(nbl, nseg, 1, 1) = ibcinfo(nbl-1, nseg, 1, 1);
                    if (iovrlp(nbl) < 0) {
                        if (ibcinfo(nbl, nseg, 1, 1) == 0) ibcinfo(nbl, nseg, 1, 1) = 1002;
                    }
                    ibcinfo(nbl, nseg, 2, 1) = ibcinfo(nbl-1, nseg, 2, 1)/2 + 1;
                    ibcinfo(nbl, nseg, 3, 1) = ibcinfo(nbl-1, nseg, 3, 1)/2 + 1;
                    ibcinfo(nbl, nseg, 4, 1) = ibcinfo(nbl-1, nseg, 4, 1)/2 + 1;
                    ibcinfo(nbl, nseg, 5, 1) = ibcinfo(nbl-1, nseg, 5, 1)/2 + 1;
                    ibcinfo(nbl, nseg, 6, 1) = ibcinfo(nbl-1, nseg, 6, 1);
                    ibcinfo(nbl, nseg, 7, 1) = ibcinfo(nbl-1, nseg, 7, 1);
                    bcfilei(nbl, nseg, 1)    = bcfilei(nbl-1, nseg, 1);
                    for (l = 1; l <= 12; l++) bcvali(nbl, nseg, l, 1) = bcvali(nbl-1, nseg, l, 1);
                }
            }
        }
        nblt = (igrid-1) * (ncgg(igrid)+1) + 1;
        isum = 0;
        for (nseg = 1; nseg <= nbci0(nblt); nseg++) {
            isum += (ibcinfo(nblt,nseg,3,1) - ibcinfo(nblt,nseg,2,1)) *
                    (ibcinfo(nblt,nseg,5,1) - ibcinfo(nblt,nseg,4,1));
        }
        if (isum < (jdimg(nblt)-1)*(kdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  I0 BCs do not span the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in j and 1-%4d in k\n", jdimg(nblt), kdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (isum > (jdimg(nblt)-1)*(kdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  I0 BCs overspan the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in j and 1-%4d in k\n", jdimg(nblt), kdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    } // end do 1207 igrid


    // ---- IDIM Boundary ----
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    fprintf(fortran_get_unit(iunit11),
        "idim:  grid   segment    bctype      jsta      jend      ksta      kend     ndata\n");
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        ifor = iforce(nbl);
        ifo  = ifor / 100;
        jfo  = (ifor - ifo*100) / 10;
        kfo  = (ifor - ifo*100 - jfo*10);
        for (iseg = 1; iseg <= nbcidim(nbl); iseg++) {
            ibcinfo(nbl, iseg, 1, 2) = -99;
            ibcinfo(nbl, iseg, 6, 2) = 0;
            if (ifo==2 || ifo==3) ibcinfo(nbl, iseg, 6, 2) = 1;
            bcfilei(nbl, iseg, 2) = 1;
            for (mm = 1; mm <= 12; mm++) bcvali(nbl, iseg, mm, 2) = -1.0e15;
        }
        for (iseg = 1; iseg <= nbcidim(nbl); iseg++) {
            { int _v2,_v3,_v4,_v5;
              fortran_read_list(iunit5, &ig, &nseg, &ibctyp, &_v2, &_v3, &_v4, &_v5, &ndata);
              ibcinfo(nbl, std::abs(nseg), 2, 2) = _v2;
              ibcinfo(nbl, std::abs(nseg), 3, 2) = _v3;
              ibcinfo(nbl, std::abs(nseg), 4, 2) = _v4;
              ibcinfo(nbl, std::abs(nseg), 5, 2) = _v5; }
            if (ibctyp == 1004) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...bc1004 no longer available...use 2004 instead\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if ((std::abs(ibctyp)==2004 || std::abs(ibctyp)==2014 ||
                 std::abs(ibctyp)==2024 || std::abs(ibctyp)==2034 ||
                 std::abs(ibctyp)==2016) && iviscg(nbl,1)==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...must have ivisc(i) > 0 if bc2004/14/24/34/16 is being used on an i-face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(ibctyp)==2014 && ivmx<=3) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2014 can only be used for turb models ivisc>3\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(ibctyp)==2024 && ivmx!=30) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2024 can only be used for turb models ivisc=30\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibctyp>=2000 && ibctyp<3000 && ndata==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2000 series bc's require abs(ndata) > 0\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibctyp==9999 && iexact_trunc==0 && iexact_disc==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...9999 bc requires keyword iexact_trunc or iexact_disc\n");
            }
            itemp1 = ibcinfo(nbl, std::abs(nseg), 2, 2);
            itemp2 = ibcinfo(nbl, std::abs(nseg), 3, 2);
            if (itemp1 > itemp2) {
                ibcinfo(nbl, std::abs(nseg), 2, 2) = itemp2;
                ibcinfo(nbl, std::abs(nseg), 3, 2) = itemp1;
            }
            itemp1 = ibcinfo(nbl, std::abs(nseg), 4, 2);
            itemp2 = ibcinfo(nbl, std::abs(nseg), 5, 2);
            if (itemp1 > itemp2) {
                ibcinfo(nbl, std::abs(nseg), 4, 2) = itemp2;
                ibcinfo(nbl, std::abs(nseg), 5, 2) = itemp1;
            }
            ifoseg = 1;
            if (nseg < 0) ifoseg = 0;
            mfoseg = 0;
            if (std::abs(ibctyp)==2004 || std::abs(ibctyp)==1005 ||
                std::abs(ibctyp)==1006 || std::abs(ibctyp)==2014 ||
                std::abs(ibctyp)==2024 || std::abs(ibctyp)==2034 ||
                std::abs(ibctyp)==2016) mfoseg = 1;
            ifoseg = ifoseg * mfoseg;
            nseg = std::abs(nseg);
            if (ndata > 0) {
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                fortran_read_list_array(iunit5, realval+1, ndata);
                for (mm = 1; mm <= ndata; mm++) bcvali(nbl, nseg, mm, 2) = realval[mm];
            }
            if (ndata < 0) {
                nfiles = nfiles + 1;
                if (nfiles > mxbcfil) {
                    fprintf(fortran_get_unit(iunit11),
                        " too many bc files specified...increase parameter mxbcfil\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                bcfilei(nbl, nseg, 2) = nfiles;
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                {
                    char _line[256]; fgets(_line, sizeof(_line), fortran_get_unit(iunit5));
                    int _len = (int)strlen(_line);
                    if (_len > 0 && _line[_len-1] == '\n') _line[--_len] = '\0';
                    std::memset(bcfiles(nfiles), ' ', 80);
                    std::memcpy(bcfiles(nfiles), _line, std::min(_len, 60));
                }
            }
            if (ig != igrid) {
                fprintf(fortran_get_unit(iunit11), " you must put these lines in order by grid!\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibcinfo(nbl, nseg, 1, 2) != -99) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...attempting to set data for segment %d more than once\n", nseg);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            ibcinfo(nbl, nseg, 1, 2) = ibctyp;
            if (ibcinfo(nbl,nseg,2,2)==0 && ibcinfo(nbl,nseg,3,2)==0) {
                ibcinfo(nbl, nseg, 2, 2) = 1; ibcinfo(nbl, nseg, 3, 2) = jdimg(nbl);
            }
            if (ibcinfo(nbl,nseg,4,2)==0 && ibcinfo(nbl,nseg,5,2)==0) {
                ibcinfo(nbl, nseg, 4, 2) = 1; ibcinfo(nbl, nseg, 5, 2) = kdimg(nbl);
            }
            if (ifoseg==0 || ibctyp==0) ibcinfo(nbl, nseg, 6, 2) = 0;
            ibcinfo(nbl, nseg, 7, 2) = ndata;
            fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    ig, nseg, ibcinfo(nbl,nseg,1,2), ibcinfo(nbl,nseg,2,2),
                    ibcinfo(nbl,nseg,3,2), ibcinfo(nbl,nseg,4,2),
                    ibcinfo(nbl,nseg,5,2), ibcinfo(nbl,nseg,7,2));
            if (ndata > 0) {
                { FortranArray1DRef<char[10]> datahdr_ref(datahdr, 10); getdhdr_ns::getdhdr(datahdr_ref, ibctyp, ndata); }
                for (mm = 1; mm <= ndata; mm++) fprintf(fortran_get_unit(iunit11), "%-10.10s", datahdr[mm]);
                fprintf(fortran_get_unit(iunit11), "\n");
                for (mm = 1; mm <= std::min(ndata,4); mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.5f", (float)bcvali(nbl,nseg,mm,2));
                for (mm = 5; mm <= ndata; mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.3e", (float)bcvali(nbl,nseg,mm,2));
                if (ndata > 0) fprintf(fortran_get_unit(iunit11), "\n");
            }
            if (ndata < 0) {
                fprintf(fortran_get_unit(iunit11), "     bc data file\n");
                nfl = bcfilei(nbl, nseg, 2);
                fprintf(fortran_get_unit(iunit11), "     %-60.60s\n", bcfiles(nfl));
            }
            if (ibcinfo(nbl,nseg,2,2)<1 || ibcinfo(nbl,nseg,2,2)>jdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...jsta is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibcinfo(nbl,nseg,3,2)<1 || ibcinfo(nbl,nseg,3,2)>jdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...jend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibcinfo(nbl,nseg,4,2)<1 || ibcinfo(nbl,nseg,4,2)>kdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...ksta is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ibcinfo(nbl,nseg,5,2)<1 || ibcinfo(nbl,nseg,5,2)>kdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...kend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        } // end do 1301 iseg


        if (iemg(igrid) > 0 && ieg(nbl) != idimg(nblcg(nbl))) {
            if (nbcidim(nbl) > 1) {
                fprintf(fortran_get_unit(iunit11),
                    " error: embedded mesh boundary at i=idim must extend over entire block face\n");
                fprintf(fortran_get_unit(iunit11), "        segmentation not allowed\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            } else {
                ibcinfo(nbl, 1, 1, 2) = 21;
                ibcinfo(nbl, 1, 2, 2) = 1;
                ibcinfo(nbl, 1, 3, 2) = jdimg(nbl);
                ibcinfo(nbl, 1, 4, 2) = 1;
                ibcinfo(nbl, 1, 5, 2) = kdimg(nbl);
            }
        }
        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
                    ibcinfo(nbl, nseg, 1, 2) = ibcinfo(nbl-1, nseg, 1, 2);
                    if (iovrlp(nbl) < 0) {
                        if (ibcinfo(nbl, nseg, 1, 2) == 0) ibcinfo(nbl, nseg, 1, 2) = 1002;
                    }
                    ibcinfo(nbl, nseg, 2, 2) = ibcinfo(nbl-1, nseg, 2, 2)/2 + 1;
                    ibcinfo(nbl, nseg, 3, 2) = ibcinfo(nbl-1, nseg, 3, 2)/2 + 1;
                    ibcinfo(nbl, nseg, 4, 2) = ibcinfo(nbl-1, nseg, 4, 2)/2 + 1;
                    ibcinfo(nbl, nseg, 5, 2) = ibcinfo(nbl-1, nseg, 5, 2)/2 + 1;
                    ibcinfo(nbl, nseg, 6, 2) = ibcinfo(nbl-1, nseg, 6, 2);
                    ibcinfo(nbl, nseg, 7, 2) = ibcinfo(nbl-1, nseg, 7, 2);
                    bcfilei(nbl, nseg, 2)    = bcfilei(nbl-1, nseg, 2);
                    for (l = 1; l <= 12; l++) bcvali(nbl, nseg, l, 2) = bcvali(nbl-1, nseg, l, 2);
                }
            }
        }
        nblt = (igrid-1) * (ncgg(igrid)+1) + 1;
        isum = 0;
        for (nseg = 1; nseg <= nbcidim(nblt); nseg++) {
            isum += (ibcinfo(nblt,nseg,3,2) - ibcinfo(nblt,nseg,2,2)) *
                    (ibcinfo(nblt,nseg,5,2) - ibcinfo(nblt,nseg,4,2));
        }
        if (isum < (jdimg(nblt)-1)*(kdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  IDIM BCs do not span the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in j and 1-%4d in k\n", jdimg(nblt), kdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (isum > (jdimg(nblt)-1)*(kdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  IDIM BCs overspan the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in j and 1-%4d in k\n", jdimg(nblt), kdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    } // end do 1307 igrid


    // ---- J0 Boundary ----
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    fprintf(fortran_get_unit(iunit11),
        "j0:   grid   segment    bctype      ista      iend      ksta      kend     ndata\n");
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        ifor = iforce(nbl);
        ifo  = ifor / 100;
        jfo  = (ifor - ifo*100) / 10;
        kfo  = (ifor - ifo*100 - jfo*10);
        for (iseg = 1; iseg <= nbcj0(nbl); iseg++) {
            jbcinfo(nbl, iseg, 1, 1) = -99;
            jbcinfo(nbl, iseg, 6, 1) = 0;
            if (jfo==1 || jfo==3) jbcinfo(nbl, iseg, 6, 1) = 1;
            bcfilej(nbl, iseg, 1) = 1;
            for (mm = 1; mm <= 12; mm++) bcvalj(nbl, iseg, mm, 1) = -1.0e15;
        }
        for (iseg = 1; iseg <= nbcj0(nbl); iseg++) {
            int jbctyp, jfoseg;
            { int _v2,_v3,_v4,_v5;
              fortran_read_list(iunit5, &ig, &nseg, &jbctyp, &_v2, &_v3, &_v4, &_v5, &ndata);
              jbcinfo(nbl, std::abs(nseg), 2, 1) = _v2;
              jbcinfo(nbl, std::abs(nseg), 3, 1) = _v3;
              jbcinfo(nbl, std::abs(nseg), 4, 1) = _v4;
              jbcinfo(nbl, std::abs(nseg), 5, 1) = _v5; }
            if (jbctyp == 1004) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...bc1004 no longer available...use 2004 instead\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if ((std::abs(jbctyp)==2004 || std::abs(jbctyp)==2014 ||
                 std::abs(jbctyp)==2024 || std::abs(jbctyp)==2034 ||
                 std::abs(jbctyp)==2016) && iviscg(nbl,2)==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...must have ivisc(j) > 0 if bc2004/14/24/34/16 is being used on a j-face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(jbctyp)==2014 && ivmx<=3) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2014 can only be used for turb models ivisc>3\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(jbctyp)==2024 && ivmx!=30) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2024 can only be used for turb models ivisc=30\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbctyp>=2000 && jbctyp<3000 && ndata==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2000 series bc's require abs(ndata) > 0\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbctyp==9999 && iexact_trunc==0 && iexact_disc==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...9999 bc requires keyword iexact_trunc or iexact_disc\n");
            }
            itemp1 = jbcinfo(nbl, std::abs(nseg), 2, 1);
            itemp2 = jbcinfo(nbl, std::abs(nseg), 3, 1);
            if (itemp1 > itemp2) {
                jbcinfo(nbl, std::abs(nseg), 2, 1) = itemp2;
                jbcinfo(nbl, std::abs(nseg), 3, 1) = itemp1;
            }
            itemp1 = jbcinfo(nbl, std::abs(nseg), 4, 1);
            itemp2 = jbcinfo(nbl, std::abs(nseg), 5, 1);
            if (itemp1 > itemp2) {
                jbcinfo(nbl, std::abs(nseg), 4, 1) = itemp2;
                jbcinfo(nbl, std::abs(nseg), 5, 1) = itemp1;
            }
            jfoseg = 1;
            if (nseg < 0) jfoseg = 0;
            mfoseg = 0;
            if (std::abs(jbctyp)==2004 || std::abs(jbctyp)==1005 ||
                std::abs(jbctyp)==1006 || std::abs(jbctyp)==2014 ||
                std::abs(jbctyp)==2024 || std::abs(jbctyp)==2034 ||
                std::abs(jbctyp)==2016) mfoseg = 1;
            jfoseg = jfoseg * mfoseg;
            nseg = std::abs(nseg);
            if (ndata > 0) {
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                fortran_read_list_array(iunit5, realval+1, ndata);
                for (mm = 1; mm <= ndata; mm++) bcvalj(nbl, nseg, mm, 1) = realval[mm];
            }
            if (ndata < 0) {
                nfiles = nfiles + 1;
                if (nfiles > mxbcfil) {
                    fprintf(fortran_get_unit(iunit11),
                        " too many bc files specified...increase parameter mxbcfil\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                bcfilej(nbl, nseg, 1) = nfiles;
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                {
                    char _line[256]; fgets(_line, sizeof(_line), fortran_get_unit(iunit5));
                    int _len = (int)strlen(_line);
                    if (_len > 0 && _line[_len-1] == '\n') _line[--_len] = '\0';
                    std::memset(bcfiles(nfiles), ' ', 80);
                    std::memcpy(bcfiles(nfiles), _line, std::min(_len, 60));
                }
            }
            if (ig != igrid) {
                fprintf(fortran_get_unit(iunit11), " you must put these lines in order by grid!\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbcinfo(nbl, nseg, 1, 1) != -99) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...attempting to set data for segment %d more than once\n", nseg);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            jbcinfo(nbl, nseg, 1, 1) = jbctyp;
            if (jbcinfo(nbl,nseg,2,1)==0 && jbcinfo(nbl,nseg,3,1)==0) {
                jbcinfo(nbl, nseg, 2, 1) = 1; jbcinfo(nbl, nseg, 3, 1) = idimg(nbl);
            }
            if (jbcinfo(nbl,nseg,4,1)==0 && jbcinfo(nbl,nseg,5,1)==0) {
                jbcinfo(nbl, nseg, 4, 1) = 1; jbcinfo(nbl, nseg, 5, 1) = kdimg(nbl);
            }
            if (jfoseg==0 || jbctyp==0) jbcinfo(nbl, nseg, 6, 1) = 0;
            jbcinfo(nbl, nseg, 7, 1) = ndata;
            fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    ig, nseg, jbcinfo(nbl,nseg,1,1), jbcinfo(nbl,nseg,2,1),
                    jbcinfo(nbl,nseg,3,1), jbcinfo(nbl,nseg,4,1),
                    jbcinfo(nbl,nseg,5,1), jbcinfo(nbl,nseg,7,1));
            if (ndata > 0) {
                { FortranArray1DRef<char[10]> datahdr_ref(datahdr, 10); getdhdr_ns::getdhdr(datahdr_ref, jbctyp, ndata); }
                for (mm = 1; mm <= ndata; mm++) fprintf(fortran_get_unit(iunit11), "%-10.10s", datahdr[mm]);
                fprintf(fortran_get_unit(iunit11), "\n");
                for (mm = 1; mm <= std::min(ndata,4); mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.5f", (float)bcvalj(nbl,nseg,mm,1));
                for (mm = 5; mm <= ndata; mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.3e", (float)bcvalj(nbl,nseg,mm,1));
                if (ndata > 0) fprintf(fortran_get_unit(iunit11), "\n");
            }
            if (ndata < 0) {
                fprintf(fortran_get_unit(iunit11), "     bc data file\n");
                nfl = bcfilej(nbl, nseg, 1);
                fprintf(fortran_get_unit(iunit11), "     %-60.60s\n", bcfiles(nfl));
            }
            if (jbcinfo(nbl,nseg,2,1)<1 || jbcinfo(nbl,nseg,2,1)>idimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...ista is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbcinfo(nbl,nseg,3,1)<1 || jbcinfo(nbl,nseg,3,1)>idimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...iend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbcinfo(nbl,nseg,4,1)<1 || jbcinfo(nbl,nseg,4,1)>kdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...ksta is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbcinfo(nbl,nseg,5,1)<1 || jbcinfo(nbl,nseg,5,1)>kdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...kend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        } // end do 1401 iseg


        if (iemg(igrid) > 0 && jsg(nbl) != 1) {
            if (nbcj0(nbl) > 1) {
                fprintf(fortran_get_unit(iunit11),
                    " error: embedded mesh boundary at j=1 must extend over entire block face\n");
                fprintf(fortran_get_unit(iunit11), "        segmentation not allowed\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            } else {
                jbcinfo(nbl, 1, 1, 1) = 21;
                jbcinfo(nbl, 1, 2, 1) = 1;
                jbcinfo(nbl, 1, 3, 1) = idimg(nbl);
                jbcinfo(nbl, 1, 4, 1) = 1;
                jbcinfo(nbl, 1, 5, 1) = kdimg(nbl);
            }
        }
        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
                    jbcinfo(nbl, nseg, 1, 1) = jbcinfo(nbl-1, nseg, 1, 1);
                    if (iovrlp(nbl) < 0) {
                        if (jbcinfo(nbl, nseg, 1, 1) == 0) jbcinfo(nbl, nseg, 1, 1) = 1002;
                    }
                    jbcinfo(nbl, nseg, 2, 1) = jbcinfo(nbl-1, nseg, 2, 1)/2 + 1;
                    jbcinfo(nbl, nseg, 3, 1) = jbcinfo(nbl-1, nseg, 3, 1)/2 + 1;
                    jbcinfo(nbl, nseg, 4, 1) = jbcinfo(nbl-1, nseg, 4, 1)/2 + 1;
                    jbcinfo(nbl, nseg, 5, 1) = jbcinfo(nbl-1, nseg, 5, 1)/2 + 1;
                    jbcinfo(nbl, nseg, 6, 1) = jbcinfo(nbl-1, nseg, 6, 1);
                    jbcinfo(nbl, nseg, 7, 1) = jbcinfo(nbl-1, nseg, 7, 1);
                    bcfilej(nbl, nseg, 1)    = bcfilej(nbl-1, nseg, 1);
                    for (l = 1; l <= 12; l++) bcvalj(nbl, nseg, l, 1) = bcvalj(nbl-1, nseg, l, 1);
                }
            }
        }
        nblt = (igrid-1) * (ncgg(igrid)+1) + 1;
        isum = 0;
        for (nseg = 1; nseg <= nbcj0(nblt); nseg++) {
            isum += (jbcinfo(nblt,nseg,3,1) - jbcinfo(nblt,nseg,2,1)) *
                    (jbcinfo(nblt,nseg,5,1) - jbcinfo(nblt,nseg,4,1));
        }
        if (isum < (idimg(nblt)-1)*(kdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  J0 BCs do not span the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in i and 1-%4d in k\n", idimg(nblt), kdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (isum > (idimg(nblt)-1)*(kdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  J0 BCs overspan the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in i and 1-%4d in k\n", idimg(nblt), kdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    } // end do 1407 igrid


    // ---- JDIM Boundary ----
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    fprintf(fortran_get_unit(iunit11),
        "jdim:  grid   segment    bctype      ista      iend      ksta      kend     ndata\n");
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        ifor = iforce(nbl);
        ifo  = ifor / 100;
        jfo  = (ifor - ifo*100) / 10;
        kfo  = (ifor - ifo*100 - jfo*10);
        for (iseg = 1; iseg <= nbcjdim(nbl); iseg++) {
            jbcinfo(nbl, iseg, 1, 2) = -99;
            jbcinfo(nbl, iseg, 6, 2) = 0;
            if (jfo==2 || jfo==3) jbcinfo(nbl, iseg, 6, 2) = 1;
            bcfilej(nbl, iseg, 2) = 1;
            for (mm = 1; mm <= 12; mm++) bcvalj(nbl, iseg, mm, 2) = -1.0e15;
        }
        for (iseg = 1; iseg <= nbcjdim(nbl); iseg++) {
            int jbctyp, jfoseg;
            { int _v2,_v3,_v4,_v5;
              fortran_read_list(iunit5, &ig, &nseg, &jbctyp, &_v2, &_v3, &_v4, &_v5, &ndata);
              jbcinfo(nbl, std::abs(nseg), 2, 2) = _v2;
              jbcinfo(nbl, std::abs(nseg), 3, 2) = _v3;
              jbcinfo(nbl, std::abs(nseg), 4, 2) = _v4;
              jbcinfo(nbl, std::abs(nseg), 5, 2) = _v5; }
            if (jbctyp == 1004) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...bc1004 no longer available...use 2004 instead\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if ((std::abs(jbctyp)==2004 || std::abs(jbctyp)==2014 ||
                 std::abs(jbctyp)==2024 || std::abs(jbctyp)==2034 ||
                 std::abs(jbctyp)==2016) && iviscg(nbl,2)==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...must have ivisc(j) > 0 if bc2004/14/24/34/16 is being used on a j-face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(jbctyp)==2014 && ivmx<=3) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2014 can only be used for turb models ivisc>3\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(jbctyp)==2024 && ivmx!=30) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2024 can only be used for turb models ivisc=30\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbctyp>=2000 && jbctyp<3000 && ndata==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2000 series bc's require abs(ndata) > 0\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbctyp==9999 && iexact_trunc==0 && iexact_disc==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...9999 bc requires keyword iexact_trunc or iexact_disc\n");
            }
            itemp1 = jbcinfo(nbl, std::abs(nseg), 2, 2);
            itemp2 = jbcinfo(nbl, std::abs(nseg), 3, 2);
            if (itemp1 > itemp2) {
                jbcinfo(nbl, std::abs(nseg), 2, 2) = itemp2;
                jbcinfo(nbl, std::abs(nseg), 3, 2) = itemp1;
            }
            itemp1 = jbcinfo(nbl, std::abs(nseg), 4, 2);
            itemp2 = jbcinfo(nbl, std::abs(nseg), 5, 2);
            if (itemp1 > itemp2) {
                jbcinfo(nbl, std::abs(nseg), 4, 2) = itemp2;
                jbcinfo(nbl, std::abs(nseg), 5, 2) = itemp1;
            }
            jfoseg = 1;
            if (nseg < 0) jfoseg = 0;
            mfoseg = 0;
            if (std::abs(jbctyp)==2004 || std::abs(jbctyp)==1005 ||
                std::abs(jbctyp)==1006 || std::abs(jbctyp)==2014 ||
                std::abs(jbctyp)==2024 || std::abs(jbctyp)==2034 ||
                std::abs(jbctyp)==2016) mfoseg = 1;
            jfoseg = jfoseg * mfoseg;
            nseg = std::abs(nseg);
            if (ndata > 0) {
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                fortran_read_list_array(iunit5, realval+1, ndata);
                for (mm = 1; mm <= ndata; mm++) bcvalj(nbl, nseg, mm, 2) = realval[mm];
            }
            if (ndata < 0) {
                nfiles = nfiles + 1;
                if (nfiles > mxbcfil) {
                    fprintf(fortran_get_unit(iunit11),
                        " too many bc files specified...increase parameter mxbcfil\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                bcfilej(nbl, nseg, 2) = nfiles;
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                {
                    char _line[256]; fgets(_line, sizeof(_line), fortran_get_unit(iunit5));
                    int _len = (int)strlen(_line);
                    if (_len > 0 && _line[_len-1] == '\n') _line[--_len] = '\0';
                    std::memset(bcfiles(nfiles), ' ', 80);
                    std::memcpy(bcfiles(nfiles), _line, std::min(_len, 60));
                }
            }
            if (ig != igrid) {
                fprintf(fortran_get_unit(iunit11), " you must put these lines in order by grid!\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbcinfo(nbl, nseg, 1, 2) != -99) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...attempting to set data for segment %d more than once\n", nseg);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            jbcinfo(nbl, nseg, 1, 2) = jbctyp;
            if (jbcinfo(nbl,nseg,2,2)==0 && jbcinfo(nbl,nseg,3,2)==0) {
                jbcinfo(nbl, nseg, 2, 2) = 1; jbcinfo(nbl, nseg, 3, 2) = idimg(nbl);
            }
            if (jbcinfo(nbl,nseg,4,2)==0 && jbcinfo(nbl,nseg,5,2)==0) {
                jbcinfo(nbl, nseg, 4, 2) = 1; jbcinfo(nbl, nseg, 5, 2) = kdimg(nbl);
            }
            if (jfoseg==0 || jbctyp==0) jbcinfo(nbl, nseg, 6, 2) = 0;
            jbcinfo(nbl, nseg, 7, 2) = ndata;
            fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    ig, nseg, jbcinfo(nbl,nseg,1,2), jbcinfo(nbl,nseg,2,2),
                    jbcinfo(nbl,nseg,3,2), jbcinfo(nbl,nseg,4,2),
                    jbcinfo(nbl,nseg,5,2), jbcinfo(nbl,nseg,7,2));
            if (ndata > 0) {
                { FortranArray1DRef<char[10]> datahdr_ref(datahdr, 10); getdhdr_ns::getdhdr(datahdr_ref, jbctyp, ndata); }
                for (mm = 1; mm <= ndata; mm++) fprintf(fortran_get_unit(iunit11), "%-10.10s", datahdr[mm]);
                fprintf(fortran_get_unit(iunit11), "\n");
                for (mm = 1; mm <= std::min(ndata,4); mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.5f", (float)bcvalj(nbl,nseg,mm,2));
                for (mm = 5; mm <= ndata; mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.3e", (float)bcvalj(nbl,nseg,mm,2));
                if (ndata > 0) fprintf(fortran_get_unit(iunit11), "\n");
            }
            if (ndata < 0) {
                fprintf(fortran_get_unit(iunit11), "     bc data file\n");
                nfl = bcfilej(nbl, nseg, 2);
                fprintf(fortran_get_unit(iunit11), "     %-60.60s\n", bcfiles(nfl));
            }
            if (jbcinfo(nbl,nseg,2,2)<1 || jbcinfo(nbl,nseg,2,2)>idimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...ista is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbcinfo(nbl,nseg,3,2)<1 || jbcinfo(nbl,nseg,3,2)>idimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...iend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbcinfo(nbl,nseg,4,2)<1 || jbcinfo(nbl,nseg,4,2)>kdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...ksta is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (jbcinfo(nbl,nseg,5,2)<1 || jbcinfo(nbl,nseg,5,2)>kdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...kend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        } // end do 1501 iseg


        if (iemg(igrid) > 0 && jeg(nbl) != jdimg(nblcg(nbl))) {
            if (nbcjdim(nbl) > 1) {
                fprintf(fortran_get_unit(iunit11),
                    " error: embedded mesh boundary at j=jdim must extend over entire block face\n");
                fprintf(fortran_get_unit(iunit11), "        segmentation not allowed\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            } else {
                jbcinfo(nbl, 1, 1, 2) = 21;
                jbcinfo(nbl, 1, 2, 2) = 1;
                jbcinfo(nbl, 1, 3, 2) = idimg(nbl);
                jbcinfo(nbl, 1, 4, 2) = 1;
                jbcinfo(nbl, 1, 5, 2) = kdimg(nbl);
            }
        }
        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
                    jbcinfo(nbl, nseg, 1, 2) = jbcinfo(nbl-1, nseg, 1, 2);
                    if (iovrlp(nbl) < 0) {
                        if (jbcinfo(nbl, nseg, 1, 2) == 0) jbcinfo(nbl, nseg, 1, 2) = 1002;
                    }
                    jbcinfo(nbl, nseg, 2, 2) = jbcinfo(nbl-1, nseg, 2, 2)/2 + 1;
                    jbcinfo(nbl, nseg, 3, 2) = jbcinfo(nbl-1, nseg, 3, 2)/2 + 1;
                    jbcinfo(nbl, nseg, 4, 2) = jbcinfo(nbl-1, nseg, 4, 2)/2 + 1;
                    jbcinfo(nbl, nseg, 5, 2) = jbcinfo(nbl-1, nseg, 5, 2)/2 + 1;
                    jbcinfo(nbl, nseg, 6, 2) = jbcinfo(nbl-1, nseg, 6, 2);
                    jbcinfo(nbl, nseg, 7, 2) = jbcinfo(nbl-1, nseg, 7, 2);
                    bcfilej(nbl, nseg, 2)    = bcfilej(nbl-1, nseg, 2);
                    for (l = 1; l <= 12; l++) bcvalj(nbl, nseg, l, 2) = bcvalj(nbl-1, nseg, l, 2);
                }
            }
        }
        nblt = (igrid-1) * (ncgg(igrid)+1) + 1;
        isum = 0;
        for (nseg = 1; nseg <= nbcjdim(nblt); nseg++) {
            isum += (jbcinfo(nblt,nseg,3,2) - jbcinfo(nblt,nseg,2,2)) *
                    (jbcinfo(nblt,nseg,5,2) - jbcinfo(nblt,nseg,4,2));
        }
        if (isum < (idimg(nblt)-1)*(kdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  JDIM BCs do not span the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in i and 1-%4d in k\n", idimg(nblt), kdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (isum > (idimg(nblt)-1)*(kdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  JDIM BCs overspan the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in i and 1-%4d in k\n", idimg(nblt), kdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    } // end do 1507 igrid


    // ---- K0 Boundary ----
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    fprintf(fortran_get_unit(iunit11),
        "k0:   grid   segment    bctype      ista      iend      jsta      jend     ndata\n");
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        ifor = iforce(nbl);
        ifo  = ifor / 100;
        jfo  = (ifor - ifo*100) / 10;
        kfo  = (ifor - ifo*100 - jfo*10);
        for (iseg = 1; iseg <= nbck0(nbl); iseg++) {
            kbcinfo(nbl, iseg, 1, 1) = -99;
            kbcinfo(nbl, iseg, 6, 1) = 0;
            if (kfo==1 || kfo==3) kbcinfo(nbl, iseg, 6, 1) = 1;
            bcfilek(nbl, iseg, 1) = 1;
            for (mm = 1; mm <= 12; mm++) bcvalk(nbl, iseg, mm, 1) = -1.0e15;
        }
        for (iseg = 1; iseg <= nbck0(nbl); iseg++) {
            int kbctyp, kfoseg;
            { int _v2,_v3,_v4,_v5;
              fortran_read_list(iunit5, &ig, &nseg, &kbctyp, &_v2, &_v3, &_v4, &_v5, &ndata);
              kbcinfo(nbl, std::abs(nseg), 2, 1) = _v2;
              kbcinfo(nbl, std::abs(nseg), 3, 1) = _v3;
              kbcinfo(nbl, std::abs(nseg), 4, 1) = _v4;
              kbcinfo(nbl, std::abs(nseg), 5, 1) = _v5; }
            if (kbctyp == 1004) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...bc1004 no longer available...use 2004 instead\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if ((std::abs(kbctyp)==2004 || std::abs(kbctyp)==2014 ||
                 std::abs(kbctyp)==2024 || std::abs(kbctyp)==2034 ||
                 std::abs(kbctyp)==2016) && iviscg(nbl,3)==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...must have ivisc(k) > 0 if bc2004/14/24/34/16 is being used on a k-face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(kbctyp)==2014 && ivmx<=3) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2014 can only be used for turb models ivisc>3\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(kbctyp)==2024 && ivmx!=30) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2024 can only be used for turb models ivisc=30\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbctyp>=2000 && kbctyp<3000 && ndata==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2000 series bc's require abs(ndata) > 0\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbctyp==9999 && iexact_trunc==0 && iexact_disc==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...9999 bc requires keyword iexact_trunc or iexact_disc\n");
            }
            itemp1 = kbcinfo(nbl, std::abs(nseg), 2, 1);
            itemp2 = kbcinfo(nbl, std::abs(nseg), 3, 1);
            if (itemp1 > itemp2) {
                kbcinfo(nbl, std::abs(nseg), 2, 1) = itemp2;
                kbcinfo(nbl, std::abs(nseg), 3, 1) = itemp1;
            }
            itemp1 = kbcinfo(nbl, std::abs(nseg), 4, 1);
            itemp2 = kbcinfo(nbl, std::abs(nseg), 5, 1);
            if (itemp1 > itemp2) {
                kbcinfo(nbl, std::abs(nseg), 4, 1) = itemp2;
                kbcinfo(nbl, std::abs(nseg), 5, 1) = itemp1;
            }
            kfoseg = 1;
            if (nseg < 0) kfoseg = 0;
            mfoseg = 0;
            if (std::abs(kbctyp)==2004 || std::abs(kbctyp)==1005 ||
                std::abs(kbctyp)==1006 || std::abs(kbctyp)==2014 ||
                std::abs(kbctyp)==2024 || std::abs(kbctyp)==2034 ||
                std::abs(kbctyp)==2016) mfoseg = 1;
            kfoseg = kfoseg * mfoseg;
            nseg = std::abs(nseg);
            if (ndata > 0) {
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                fortran_read_list_array(iunit5, realval+1, ndata);
                for (mm = 1; mm <= ndata; mm++) bcvalk(nbl, nseg, mm, 1) = realval[mm];
            }
            if (ndata < 0) {
                nfiles = nfiles + 1;
                if (nfiles > mxbcfil) {
                    fprintf(fortran_get_unit(iunit11),
                        " too many bc files specified...increase parameter mxbcfil\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                bcfilek(nbl, nseg, 1) = nfiles;
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                {
                    char _line[256]; fgets(_line, sizeof(_line), fortran_get_unit(iunit5));
                    int _len = (int)strlen(_line);
                    if (_len > 0 && _line[_len-1] == '\n') _line[--_len] = '\0';
                    std::memset(bcfiles(nfiles), ' ', 80);
                    std::memcpy(bcfiles(nfiles), _line, std::min(_len, 60));
                }
            }
            if (ig != igrid) {
                fprintf(fortran_get_unit(iunit11), " you must put these lines in order by grid!\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbcinfo(nbl, nseg, 1, 1) != -99) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...attempting to set data for segment %d more than once\n", nseg);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            kbcinfo(nbl, nseg, 1, 1) = kbctyp;
            if (kbcinfo(nbl,nseg,2,1)==0 && kbcinfo(nbl,nseg,3,1)==0) {
                kbcinfo(nbl, nseg, 2, 1) = 1; kbcinfo(nbl, nseg, 3, 1) = idimg(nbl);
            }
            if (kbcinfo(nbl,nseg,4,1)==0 && kbcinfo(nbl,nseg,5,1)==0) {
                kbcinfo(nbl, nseg, 4, 1) = 1; kbcinfo(nbl, nseg, 5, 1) = jdimg(nbl);
            }
            if (kfoseg==0 || kbctyp==0) kbcinfo(nbl, nseg, 6, 1) = 0;
            kbcinfo(nbl, nseg, 7, 1) = ndata;
            fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    ig, nseg, kbcinfo(nbl,nseg,1,1), kbcinfo(nbl,nseg,2,1),
                    kbcinfo(nbl,nseg,3,1), kbcinfo(nbl,nseg,4,1),
                    kbcinfo(nbl,nseg,5,1), kbcinfo(nbl,nseg,7,1));
            if (ndata > 0) {
                { FortranArray1DRef<char[10]> datahdr_ref(datahdr, 10); getdhdr_ns::getdhdr(datahdr_ref, kbctyp, ndata); }
                for (mm = 1; mm <= ndata; mm++) fprintf(fortran_get_unit(iunit11), "%-10.10s", datahdr[mm]);
                fprintf(fortran_get_unit(iunit11), "\n");
                for (mm = 1; mm <= std::min(ndata,4); mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.5f", (float)bcvalk(nbl,nseg,mm,1));
                for (mm = 5; mm <= ndata; mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.3e", (float)bcvalk(nbl,nseg,mm,1));
                if (ndata > 0) fprintf(fortran_get_unit(iunit11), "\n");
            }
            if (ndata < 0) {
                fprintf(fortran_get_unit(iunit11), "     bc data file\n");
                nfl = bcfilek(nbl, nseg, 1);
                fprintf(fortran_get_unit(iunit11), "     %-60.60s\n", bcfiles(nfl));
            }
            if (kbcinfo(nbl,nseg,2,1)<1 || kbcinfo(nbl,nseg,2,1)>idimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...ista is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbcinfo(nbl,nseg,3,1)<1 || kbcinfo(nbl,nseg,3,1)>idimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...iend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbcinfo(nbl,nseg,4,1)<1 || kbcinfo(nbl,nseg,4,1)>jdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...jsta is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbcinfo(nbl,nseg,5,1)<1 || kbcinfo(nbl,nseg,5,1)>jdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...jend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        } // end do 1601 iseg


        if (iemg(igrid) > 0 && ksg(nbl) != 1) {
            if (nbck0(nbl) > 1) {
                fprintf(fortran_get_unit(iunit11),
                    " error: embedded mesh boundary at k=1 must extend over entire block face\n");
                fprintf(fortran_get_unit(iunit11), "        segmentation not allowed\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            } else {
                kbcinfo(nbl, 1, 1, 1) = 21;
                kbcinfo(nbl, 1, 2, 1) = 1;
                kbcinfo(nbl, 1, 3, 1) = idimg(nbl);
                kbcinfo(nbl, 1, 4, 1) = 1;
                kbcinfo(nbl, 1, 5, 1) = jdimg(nbl);
            }
        }
        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
                    kbcinfo(nbl, nseg, 1, 1) = kbcinfo(nbl-1, nseg, 1, 1);
                    if (iovrlp(nbl) < 0) {
                        if (kbcinfo(nbl, nseg, 1, 1) == 0) kbcinfo(nbl, nseg, 1, 1) = 1002;
                    }
                    kbcinfo(nbl, nseg, 2, 1) = kbcinfo(nbl-1, nseg, 2, 1)/2 + 1;
                    kbcinfo(nbl, nseg, 3, 1) = kbcinfo(nbl-1, nseg, 3, 1)/2 + 1;
                    kbcinfo(nbl, nseg, 4, 1) = kbcinfo(nbl-1, nseg, 4, 1)/2 + 1;
                    kbcinfo(nbl, nseg, 5, 1) = kbcinfo(nbl-1, nseg, 5, 1)/2 + 1;
                    kbcinfo(nbl, nseg, 6, 1) = kbcinfo(nbl-1, nseg, 6, 1);
                    kbcinfo(nbl, nseg, 7, 1) = kbcinfo(nbl-1, nseg, 7, 1);
                    bcfilek(nbl, nseg, 1)    = bcfilek(nbl-1, nseg, 1);
                    for (l = 1; l <= 12; l++) bcvalk(nbl, nseg, l, 1) = bcvalk(nbl-1, nseg, l, 1);
                }
            }
        }
        nblt = (igrid-1) * (ncgg(igrid)+1) + 1;
        isum = 0;
        for (nseg = 1; nseg <= nbck0(nblt); nseg++) {
            isum += (kbcinfo(nblt,nseg,3,1) - kbcinfo(nblt,nseg,2,1)) *
                    (kbcinfo(nblt,nseg,5,1) - kbcinfo(nblt,nseg,4,1));
        }
        if (isum < (idimg(nblt)-1)*(jdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  K0 BCs do not span the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in i and 1-%4d in j\n", idimg(nblt), jdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (isum > (idimg(nblt)-1)*(jdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  K0 BCs overspan the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in i and 1-%4d in j\n", idimg(nblt), jdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    } // end do 1607 igrid


    // ---- KDIM Boundary ----
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    nbl = 0;
    fprintf(fortran_get_unit(iunit11),
        "kdim:  grid   segment    bctype      ista      iend      jsta      jend     ndata\n");
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        ifor = iforce(nbl);
        ifo  = ifor / 100;
        jfo  = (ifor - ifo*100) / 10;
        kfo  = (ifor - ifo*100 - jfo*10);
        for (iseg = 1; iseg <= nbckdim(nbl); iseg++) {
            kbcinfo(nbl, iseg, 1, 2) = -99;
            kbcinfo(nbl, iseg, 6, 2) = 0;
            if (kfo==2 || kfo==3) kbcinfo(nbl, iseg, 6, 2) = 1;
            bcfilek(nbl, iseg, 2) = 1;
            for (mm = 1; mm <= 12; mm++) bcvalk(nbl, iseg, mm, 2) = -1.0e15;
        }
        for (iseg = 1; iseg <= nbckdim(nbl); iseg++) {
            int kbctyp, kfoseg;
            { int _v2,_v3,_v4,_v5;
              fortran_read_list(iunit5, &ig, &nseg, &kbctyp, &_v2, &_v3, &_v4, &_v5, &ndata);
              kbcinfo(nbl, std::abs(nseg), 2, 2) = _v2;
              kbcinfo(nbl, std::abs(nseg), 3, 2) = _v3;
              kbcinfo(nbl, std::abs(nseg), 4, 2) = _v4;
              kbcinfo(nbl, std::abs(nseg), 5, 2) = _v5; }
            if (kbctyp == 1004) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...bc1004 no longer available...use 2004 instead\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if ((std::abs(kbctyp)==2004 || std::abs(kbctyp)==2014 ||
                 std::abs(kbctyp)==2024 || std::abs(kbctyp)==2034 ||
                 std::abs(kbctyp)==2016) && iviscg(nbl,3)==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...must have ivisc(k) > 0 if bc2004/14/24/34/16 is being used on a k-face\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(kbctyp)==2014 && ivmx<=3) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2014 can only be used for turb models ivisc>3\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (std::abs(kbctyp)==2024 && ivmx!=30) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2024 can only be used for turb models ivisc=30\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbctyp>=2000 && kbctyp<3000 && ndata==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...2000 series bc's require abs(ndata) > 0\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbctyp==9999 && iexact_trunc==0 && iexact_disc==0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...9999 bc requires keyword iexact_trunc or iexact_disc\n");
            }
            itemp1 = kbcinfo(nbl, std::abs(nseg), 2, 2);
            itemp2 = kbcinfo(nbl, std::abs(nseg), 3, 2);
            if (itemp1 > itemp2) {
                kbcinfo(nbl, std::abs(nseg), 2, 2) = itemp2;
                kbcinfo(nbl, std::abs(nseg), 3, 2) = itemp1;
            }
            itemp1 = kbcinfo(nbl, std::abs(nseg), 4, 2);
            itemp2 = kbcinfo(nbl, std::abs(nseg), 5, 2);
            if (itemp1 > itemp2) {
                kbcinfo(nbl, std::abs(nseg), 4, 2) = itemp2;
                kbcinfo(nbl, std::abs(nseg), 5, 2) = itemp1;
            }
            kfoseg = 1;
            if (nseg < 0) kfoseg = 0;
            mfoseg = 0;
            if (std::abs(kbctyp)==2004 || std::abs(kbctyp)==1005 ||
                std::abs(kbctyp)==1006 || std::abs(kbctyp)==2014 ||
                std::abs(kbctyp)==2024 || std::abs(kbctyp)==2034 ||
                std::abs(kbctyp)==2016) mfoseg = 1;
            kfoseg = kfoseg * mfoseg;
            nseg = std::abs(nseg);
            if (ndata > 0) {
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                fortran_read_list_array(iunit5, realval+1, ndata);
                for (mm = 1; mm <= ndata; mm++) bcvalk(nbl, nseg, mm, 2) = realval[mm];
            }
            if (ndata < 0) {
                nfiles = nfiles + 1;
                if (nfiles > mxbcfil) {
                    fprintf(fortran_get_unit(iunit11),
                        " too many bc files specified...increase parameter mxbcfil\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                bcfilek(nbl, nseg, 2) = nfiles;
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                {
                    char _line[256]; fgets(_line, sizeof(_line), fortran_get_unit(iunit5));
                    int _len = (int)strlen(_line);
                    if (_len > 0 && _line[_len-1] == '\n') _line[--_len] = '\0';
                    std::memset(bcfiles(nfiles), ' ', 80);
                    std::memcpy(bcfiles(nfiles), _line, std::min(_len, 60));
                }
            }
            if (ig != igrid) {
                fprintf(fortran_get_unit(iunit11), " you must put these lines in order by grid!\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbcinfo(nbl, nseg, 1, 2) != -99) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...attempting to set data for segment %d more than once\n", nseg);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            kbcinfo(nbl, nseg, 1, 2) = kbctyp;
            if (kbcinfo(nbl,nseg,2,2)==0 && kbcinfo(nbl,nseg,3,2)==0) {
                kbcinfo(nbl, nseg, 2, 2) = 1; kbcinfo(nbl, nseg, 3, 2) = idimg(nbl);
            }
            if (kbcinfo(nbl,nseg,4,2)==0 && kbcinfo(nbl,nseg,5,2)==0) {
                kbcinfo(nbl, nseg, 4, 2) = 1; kbcinfo(nbl, nseg, 5, 2) = jdimg(nbl);
            }
            if (kfoseg==0 || kbctyp==0) kbcinfo(nbl, nseg, 6, 2) = 0;
            kbcinfo(nbl, nseg, 7, 2) = ndata;
            fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    ig, nseg, kbcinfo(nbl,nseg,1,2), kbcinfo(nbl,nseg,2,2),
                    kbcinfo(nbl,nseg,3,2), kbcinfo(nbl,nseg,4,2),
                    kbcinfo(nbl,nseg,5,2), kbcinfo(nbl,nseg,7,2));
            if (ndata > 0) {
                { FortranArray1DRef<char[10]> datahdr_ref(datahdr, 10); getdhdr_ns::getdhdr(datahdr_ref, kbctyp, ndata); }
                for (mm = 1; mm <= ndata; mm++) fprintf(fortran_get_unit(iunit11), "%-10.10s", datahdr[mm]);
                fprintf(fortran_get_unit(iunit11), "\n");
                for (mm = 1; mm <= std::min(ndata,4); mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.5f", (float)bcvalk(nbl,nseg,mm,2));
                for (mm = 5; mm <= ndata; mm++)
                    fprintf(fortran_get_unit(iunit11), "%10.3e", (float)bcvalk(nbl,nseg,mm,2));
                if (ndata > 0) fprintf(fortran_get_unit(iunit11), "\n");
            }
            if (ndata < 0) {
                fprintf(fortran_get_unit(iunit11), "     bc data file\n");
                nfl = bcfilek(nbl, nseg, 2);
                fprintf(fortran_get_unit(iunit11), "     %-60.60s\n", bcfiles(nfl));
            }
            if (kbcinfo(nbl,nseg,2,2)<1 || kbcinfo(nbl,nseg,2,2)>idimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...ista is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbcinfo(nbl,nseg,3,2)<1 || kbcinfo(nbl,nseg,3,2)>idimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...iend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbcinfo(nbl,nseg,4,2)<1 || kbcinfo(nbl,nseg,4,2)>jdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...jsta is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (kbcinfo(nbl,nseg,5,2)<1 || kbcinfo(nbl,nseg,5,2)>jdimg(nbl)) {
                fprintf(fortran_get_unit(iunit11), " stopping...jend is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        } // end do 1701 iseg


        if (iemg(igrid) > 0 && keg(nbl) != kdimg(nblcg(nbl))) {
            if (nbckdim(nbl) > 1) {
                fprintf(fortran_get_unit(iunit11),
                    " error: embedded mesh boundary at k=kdim must extend over entire block face\n");
                fprintf(fortran_get_unit(iunit11), "        segmentation not allowed\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            } else {
                kbcinfo(nbl, 1, 1, 2) = 21;
                kbcinfo(nbl, 1, 2, 2) = 1;
                kbcinfo(nbl, 1, 3, 2) = idimg(nbl);
                kbcinfo(nbl, 1, 4, 2) = 1;
                kbcinfo(nbl, 1, 5, 2) = jdimg(nbl);
            }
        }
        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
                    kbcinfo(nbl, nseg, 1, 2) = kbcinfo(nbl-1, nseg, 1, 2);
                    if (iovrlp(nbl) < 0) {
                        if (kbcinfo(nbl, nseg, 1, 2) == 0) kbcinfo(nbl, nseg, 1, 2) = 1002;
                    }
                    kbcinfo(nbl, nseg, 2, 2) = kbcinfo(nbl-1, nseg, 2, 2)/2 + 1;
                    kbcinfo(nbl, nseg, 3, 2) = kbcinfo(nbl-1, nseg, 3, 2)/2 + 1;
                    kbcinfo(nbl, nseg, 4, 2) = kbcinfo(nbl-1, nseg, 4, 2)/2 + 1;
                    kbcinfo(nbl, nseg, 5, 2) = kbcinfo(nbl-1, nseg, 5, 2)/2 + 1;
                    kbcinfo(nbl, nseg, 6, 2) = kbcinfo(nbl-1, nseg, 6, 2);
                    kbcinfo(nbl, nseg, 7, 2) = kbcinfo(nbl-1, nseg, 7, 2);
                    bcfilek(nbl, nseg, 2)    = bcfilek(nbl-1, nseg, 2);
                    for (l = 1; l <= 12; l++) bcvalk(nbl, nseg, l, 2) = bcvalk(nbl-1, nseg, l, 2);
                }
                // reset iovrlp flag for coarser meshes (Fortran global.F:3140)
                if (iovrlp(nbl) < 0) iovrlp(nbl) = 0;
            }
        }
        nblt = (igrid-1) * (ncgg(igrid)+1) + 1;
        isum = 0;
        for (nseg = 1; nseg <= nbckdim(nblt); nseg++) {
            isum += (kbcinfo(nblt,nseg,3,2) - kbcinfo(nblt,nseg,2,2)) *
                    (kbcinfo(nblt,nseg,5,2) - kbcinfo(nblt,nseg,4,2));
        }
        if (isum < (idimg(nblt)-1)*(jdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  KDIM BCs do not span the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in i and 1-%4d in j\n", idimg(nblt), jdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (isum > (idimg(nblt)-1)*(jdimg(nblt)-1)) {
            fprintf(fortran_get_unit(iunit11),
                " Error.  KDIM BCs overspan the space for grid %5d\n", igrid);
            fprintf(fortran_get_unit(iunit11),
                "   Must exactly cover 1-%4d in i and 1-%4d in j\n", idimg(nblt), jdimg(nblt));
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    } // end do 1707 igrid


    // if nplot3d < 0, set plot3d output to contain all solid surfaces
    if (nplot3d < 0 || nprint < 0) {
        iptype = 0;
        int nplt3dtmp = 0;
        for (int ng = 1; ng <= ngrid; ng++) {
            nbl  = nblg(ng);
            idim = idimg(nbl);
            jdim = jdimg(nbl);
            kdim = kdimg(nbl);
            if (i2d == 0) {
                // i0 face
                for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
                    if (std::abs(ibcinfo(nbl,nseg,1,1))==2004 ||
                        std::abs(ibcinfo(nbl,nseg,1,1))==1005 ||
                        std::abs(ibcinfo(nbl,nseg,1,1))==1006 ||
                        std::abs(ibcinfo(nbl,nseg,1,1))==2014 ||
                        std::abs(ibcinfo(nbl,nseg,1,1))==2024 ||
                        std::abs(ibcinfo(nbl,nseg,1,1))==2034 ||
                        std::abs(ibcinfo(nbl,nseg,1,1))==2016) {
                        ii = 1; nplt3dtmp++;
                        ipl3dtmp(1,nplt3dtmp)=ng; ipl3dtmp(2,nplt3dtmp)=iptype;
                        ipl3dtmp(3,nplt3dtmp)=ii; ipl3dtmp(4,nplt3dtmp)=ii;
                        ipl3dtmp(5,nplt3dtmp)=1;
                        ipl3dtmp(6,nplt3dtmp)=ibcinfo(nbl,nseg,2,1);
                        ipl3dtmp(7,nplt3dtmp)=ibcinfo(nbl,nseg,3,1);
                        ipl3dtmp(8,nplt3dtmp)=1;
                        ipl3dtmp(9,nplt3dtmp)=ibcinfo(nbl,nseg,4,1);
                        ipl3dtmp(10,nplt3dtmp)=ibcinfo(nbl,nseg,5,1);
                        ipl3dtmp(11,nplt3dtmp)=1;
                    }
                }
                // idim face
                for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
                    if (std::abs(ibcinfo(nbl,nseg,1,2))==2004 ||
                        std::abs(ibcinfo(nbl,nseg,1,2))==1005 ||
                        std::abs(ibcinfo(nbl,nseg,1,2))==1006 ||
                        std::abs(ibcinfo(nbl,nseg,1,2))==2014 ||
                        std::abs(ibcinfo(nbl,nseg,1,2))==2024 ||
                        std::abs(ibcinfo(nbl,nseg,1,2))==2034 ||
                        std::abs(ibcinfo(nbl,nseg,1,2))==2016) {
                        ii = idim; nplt3dtmp++;
                        ipl3dtmp(1,nplt3dtmp)=ng; ipl3dtmp(2,nplt3dtmp)=iptype;
                        ipl3dtmp(3,nplt3dtmp)=ii; ipl3dtmp(4,nplt3dtmp)=ii;
                        ipl3dtmp(5,nplt3dtmp)=1;
                        ipl3dtmp(6,nplt3dtmp)=ibcinfo(nbl,nseg,2,2);
                        ipl3dtmp(7,nplt3dtmp)=ibcinfo(nbl,nseg,3,2);
                        ipl3dtmp(8,nplt3dtmp)=1;
                        ipl3dtmp(9,nplt3dtmp)=ibcinfo(nbl,nseg,4,2);
                        ipl3dtmp(10,nplt3dtmp)=ibcinfo(nbl,nseg,5,2);
                        ipl3dtmp(11,nplt3dtmp)=1;
                    }
                }
                // j0 face
                for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
                    if (std::abs(jbcinfo(nbl,nseg,1,1))==2004 ||
                        std::abs(jbcinfo(nbl,nseg,1,1))==1005 ||
                        std::abs(jbcinfo(nbl,nseg,1,1))==1006 ||
                        std::abs(jbcinfo(nbl,nseg,1,1))==2014 ||
                        std::abs(jbcinfo(nbl,nseg,1,1))==2024 ||
                        std::abs(jbcinfo(nbl,nseg,1,1))==2034 ||
                        std::abs(jbcinfo(nbl,nseg,1,1))==2016) {
                        jj = 1; nplt3dtmp++;
                        ipl3dtmp(1,nplt3dtmp)=ng; ipl3dtmp(2,nplt3dtmp)=iptype;
                        ipl3dtmp(3,nplt3dtmp)=jbcinfo(nbl,nseg,2,1);
                        ipl3dtmp(4,nplt3dtmp)=jbcinfo(nbl,nseg,3,1);
                        ipl3dtmp(5,nplt3dtmp)=1;
                        ipl3dtmp(6,nplt3dtmp)=jj; ipl3dtmp(7,nplt3dtmp)=jj;
                        ipl3dtmp(8,nplt3dtmp)=1;
                        ipl3dtmp(9,nplt3dtmp)=jbcinfo(nbl,nseg,4,1);
                        ipl3dtmp(10,nplt3dtmp)=jbcinfo(nbl,nseg,5,1);
                        ipl3dtmp(11,nplt3dtmp)=1;
                    }
                }
                // jdim face
                for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
                    if (std::abs(jbcinfo(nbl,nseg,1,2))==2004 ||
                        std::abs(jbcinfo(nbl,nseg,1,2))==1005 ||
                        std::abs(jbcinfo(nbl,nseg,1,2))==1006 ||
                        std::abs(jbcinfo(nbl,nseg,1,2))==2014 ||
                        std::abs(jbcinfo(nbl,nseg,1,2))==2024 ||
                        std::abs(jbcinfo(nbl,nseg,1,2))==2034 ||
                        std::abs(jbcinfo(nbl,nseg,1,2))==2016) {
                        jj = jdim; nplt3dtmp++;
                        ipl3dtmp(1,nplt3dtmp)=ng; ipl3dtmp(2,nplt3dtmp)=iptype;
                        ipl3dtmp(3,nplt3dtmp)=jbcinfo(nbl,nseg,2,2);
                        ipl3dtmp(4,nplt3dtmp)=jbcinfo(nbl,nseg,3,2);
                        ipl3dtmp(5,nplt3dtmp)=1;
                        ipl3dtmp(6,nplt3dtmp)=jj; ipl3dtmp(7,nplt3dtmp)=jj;
                        ipl3dtmp(8,nplt3dtmp)=1;
                        ipl3dtmp(9,nplt3dtmp)=jbcinfo(nbl,nseg,4,2);
                        ipl3dtmp(10,nplt3dtmp)=jbcinfo(nbl,nseg,5,2);
                        ipl3dtmp(11,nplt3dtmp)=1;
                    }
                }
                // k0 face
                for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
                    if (std::abs(kbcinfo(nbl,nseg,1,1))==2004 ||
                        std::abs(kbcinfo(nbl,nseg,1,1))==1005 ||
                        std::abs(kbcinfo(nbl,nseg,1,1))==1006 ||
                        std::abs(kbcinfo(nbl,nseg,1,1))==2014 ||
                        std::abs(kbcinfo(nbl,nseg,1,1))==2024 ||
                        std::abs(kbcinfo(nbl,nseg,1,1))==2034 ||
                        std::abs(kbcinfo(nbl,nseg,1,1))==2016) {
                        kk = 1; nplt3dtmp++;
                        ipl3dtmp(1,nplt3dtmp)=ng; ipl3dtmp(2,nplt3dtmp)=iptype;
                        ipl3dtmp(3,nplt3dtmp)=kbcinfo(nbl,nseg,2,1);
                        ipl3dtmp(4,nplt3dtmp)=kbcinfo(nbl,nseg,3,1);
                        ipl3dtmp(5,nplt3dtmp)=1;
                        ipl3dtmp(6,nplt3dtmp)=kbcinfo(nbl,nseg,4,1);
                        ipl3dtmp(7,nplt3dtmp)=kbcinfo(nbl,nseg,5,1);
                        ipl3dtmp(8,nplt3dtmp)=1;
                        ipl3dtmp(9,nplt3dtmp)=kk; ipl3dtmp(10,nplt3dtmp)=kk;
                        ipl3dtmp(11,nplt3dtmp)=1;
                    }
                }
                // kdim face
                for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
                    if (std::abs(kbcinfo(nbl,nseg,1,2))==2004 ||
                        std::abs(kbcinfo(nbl,nseg,1,2))==1005 ||
                        std::abs(kbcinfo(nbl,nseg,1,2))==1006 ||
                        std::abs(kbcinfo(nbl,nseg,1,2))==2014 ||
                        std::abs(kbcinfo(nbl,nseg,1,2))==2024 ||
                        std::abs(kbcinfo(nbl,nseg,1,2))==2034 ||
                        std::abs(kbcinfo(nbl,nseg,1,2))==2016) {
                        kk = kdim; nplt3dtmp++;
                        ipl3dtmp(1,nplt3dtmp)=ng; ipl3dtmp(2,nplt3dtmp)=iptype;
                        ipl3dtmp(3,nplt3dtmp)=kbcinfo(nbl,nseg,2,2);
                        ipl3dtmp(4,nplt3dtmp)=kbcinfo(nbl,nseg,3,2);
                        ipl3dtmp(5,nplt3dtmp)=1;
                        ipl3dtmp(6,nplt3dtmp)=kbcinfo(nbl,nseg,4,2);
                        ipl3dtmp(7,nplt3dtmp)=kbcinfo(nbl,nseg,5,2);
                        ipl3dtmp(8,nplt3dtmp)=1;
                        ipl3dtmp(9,nplt3dtmp)=kk; ipl3dtmp(10,nplt3dtmp)=kk;
                        ipl3dtmp(11,nplt3dtmp)=1;
                    }
                }
            } else {
                // 2D case: add all zones
                nplt3dtmp++;
                ipl3dtmp(1,nplt3dtmp)=ng; ipl3dtmp(2,nplt3dtmp)=iptype;
                ipl3dtmp(3,nplt3dtmp)=1; ipl3dtmp(4,nplt3dtmp)=1;
                ipl3dtmp(5,nplt3dtmp)=1;
                ipl3dtmp(6,nplt3dtmp)=1; ipl3dtmp(7,nplt3dtmp)=jdim;
                ipl3dtmp(8,nplt3dtmp)=1;
                ipl3dtmp(9,nplt3dtmp)=1; ipl3dtmp(10,nplt3dtmp)=kdim;
                ipl3dtmp(11,nplt3dtmp)=1;
            }
        } // end do 4800 ng


        if (nplt3dtmp > nplots) {
            fprintf(fortran_get_unit(iunit11), " nplots too small for auto-surface detection\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    } // end if (nplot3d < 0 || nprint < 0)

    // mseq, mgflag, iconsf, mtt, ngam
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    fortran_read_list(iunit5, &mseq, &mgflag, &iconsf, &mtt, &ngam);
    fprintf(fortran_get_unit(iunit11),
        "      mseq    mgflag    iconsf       mtt      ngam\n");
    fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d\n",
            mseq, mgflag, iconsf, mtt, ngam);

    if (mseq > 5) {
        fprintf(fortran_get_unit(iunit11), " Stopping... mseq cannot exceed 5\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (mseq > 1) {
        nbl = 0;
        for (igrid = 1; igrid <= ngrid; igrid++) {
            nbl = nbl + 1;
            if (iovrlp(nbl) > 0) {
                fprintf(fortran_get_unit(iunit11),
                    " Stopping!  Mesh sequencing with Chimera not allowed!\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            ncg = ncgg(igrid);
            if (ncg > 0) nbl = nbl + ncg;
        }
    }
    if (mseq > 1 && iemtot != 0 && iteravg != 0) {
        fprintf(fortran_get_unit(iunit11),
            " Stopping!  Cannot use iteravg>0 with mseq>1 and embedded grids\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (iexact_trunc != 0 && iexact_disc != 0) {
        fprintf(fortran_get_unit(iunit11),
            " Stopping!  cannot have both iexact_trunc and iexact_disc .ne. 0\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (iexact_trunc != 0) {
        if (ncg != 0 || mgflag != 0 || ifullns != 1 || i2d != 1) {
            fprintf(fortran_get_unit(iunit11),
                " Stopping!  must have ncg=0,mgflag=0,ifullns=1,i2d=1 when iexact_trunc.ne.0\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (ivmx != 5) {
            fprintf(fortran_get_unit(iunit11),
                " Stopping!  currently must use ivisc=5 for iexact_trunc.ne.0\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (iexact_trunc != 1 && iexact_trunc != 2 && iexact_trunc != 4) {
            fprintf(fortran_get_unit(iunit11),
                " Stopping!  currently must have iexact_trunc=1,2, or 4\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }
    if (iexact_disc != 0) {
        if (ifullns != 1 || i2d != 1) {
            fprintf(fortran_get_unit(iunit11),
                " Stopping!  must have ifullns=1,i2d=1 when iexact_disc.ne.0\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (ivmx != 5) {
            fprintf(fortran_get_unit(iunit11),
                " Stopping!  currently must use ivisc=5 for iexact_disc.ne.0\n");
            if (iexact_disc != 1 && iexact_disc != 2 && iexact_disc != 4) {
                fprintf(fortran_get_unit(iunit11),
                    " Stopping!  currently must have iexact_disc=1,2, or 4\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }
    if (iexact_ring != 0) {
        if (iexact_trunc == 0 && iexact_disc == 0) {
            fprintf(fortran_get_unit(iunit11),
                " iexact_ring must be used with iexact_trunc=1 or iexact_disc=1\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }

    // issc, epsssc, issr, epsssr
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    fortran_read_list(iunit5, &issc, &realval[1], &realval[2], &realval[3],
                      &issr, &realval[4], &realval[5], &realval[6]);
    for (n = 1; n <= 3; n++) {
        epsssc[n-1] = realval[n];
        epsssr[n-1] = realval[n+3];
    }
    fprintf(fortran_get_unit(iunit11),
        "      issc epsssc(1) epsssc(2) epsssc(3)      issr epsssr(1) epsssr(2) epsssr(3)\n");
    // format 360: 2(i10,3f10.4)
    fprintf(fortran_get_unit(iunit11), "%10d%10.4f%10.4f%10.4f%10d%10.4f%10.4f%10.4f\n",
            issc, (float)epsssc[0], (float)epsssc[1], (float)epsssc[2],
            issr, (float)epsssr[0], (float)epsssr[1], (float)epsssr[2]);

    // ncyc1, mglevg, nemgl, nitfo1 for each sequence
    ncyctot = 0;
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    fprintf(fortran_get_unit(iunit11),
        "      ncyc    mglevg     nemgl     nitfo\n");
    for (m = 1; m <= mseq; m++) {
        int idum;
        fortran_read_list(iunit5, &ncyc1[m-1], &mglevg(m), &nemgl(m), &nitfo1[m-1]);
        idum = mglevg(m) + nemgl(m);
        if (m > 1 && ncyc1[m-1] != 0 && nfreeze != 0) {
            fprintf(fortran_get_unit(iunit11),
                " not allowed to freeze turb model immediately upon mesh-sequencing-up to next level\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (m > 1 && ncyc1[m-1] != 0 && iteravg > 0) {
            fprintf(fortran_get_unit(iunit11),
                " not allowed to use a running-average file (iteravg>0) with mesh sequencing\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (mglevg(m) <= 0 || nemgl(m) < 0 || idum > maxbl) {
            fprintf(fortran_get_unit(iunit11),
                " incorrect input value for mglevg and/or nemgl\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (nemgl(m) > 0 && iemtot == 0) {
            fprintf(fortran_get_unit(iunit11),
                " inconsistent values assigned to nemgl and iem\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }


        if (mglevg(m) > 1 && mgflag < 1) {
            fprintf(fortran_get_unit(iunit11),
                " inconsistent values assigned to mglevg and mgflag\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (m < mseq && nemgl(m) > 0) {
            fprintf(fortran_get_unit(iunit11),
                " embedded grids only allowed on finest global grid\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if ((float)dt < 0.0f) {
            ncyctot = ncyctot + ncyc1[m-1];
        } else {
            if (ncyc1[m-1] > 0) ncyctot = ncyctot + ntstep;
        }
        if ((float)dt > 0.0f && mgflag != 0 && ncyc1[m-1] == 1) {
            fprintf(fortran_get_unit(iunit11),
                " must do subits (ncyc>1) with MG and time-accurate\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    } // end do 12 m=1,mseq

    if (mseq > 1 && ivmx >= 4) {
        int iset = 0;
        for (int iseq = 1; iseq <= mseq; iseq++) {
            if (ncyc1[iseq-1] > 0) iset++;
        }
        if (iset > 1 && irest != 0) isminc = 2;
    }
    for (igrid = 1; igrid <= ngrid; igrid++) {
        if (inewgg(igrid) > 0 && ivmx >= 4) isminc = 2;
    }

    // determine levelt and levelb
    fprintf(fortran_get_unit(iunit11),
        "      ncyc    mglevg     nemgl     nitfo\n");
    for (m = 1; m <= mseq; m++) {
        levelt[m-1] = ncgmax - (mseq - m) + nemgl(m) + 1;
        levelb[m-1] = levelt[m-1] - (mglevg(m) - 1) - nemgl(m);
        if (levelb[m-1] < 1) {
            fprintf(fortran_get_unit(iunit11),
                " error in input, m, levelt, levelb, ncgmax=%5d%5d%5d%5d\n",
                m, levelt[m-1], levelb[m-1], ncgmax);
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d\n",
                ncyc1[m-1], mglevg(m), nemgl(m), nitfo1[m-1]);
        if (std::abs(ita) > 1 && ncyc1[m-1] == 1) {
            fprintf(fortran_get_unit(iunit11),
                " you MUST use subiterations (ncyc>1) to get second-order accuracy in time!\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }

    if (icall > 0) {
        if (ncyctot > ncycmax) {
            fprintf(fortran_get_unit(iunit11),
                " number of cycles (%5d) exceeds ncycmax (%5d)\n", ncyctot, ncycmax);
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }
    if (ncyc1[0] <= 0) {
        fprintf(fortran_get_unit(iunit11),
            " you MUST do at least one cycle on the starting level!\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    // mit - iterations on each level
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    fprintf(fortran_get_unit(iunit11),
        "      mit1      mit2      mit3      mit4      mit5  ...\n");
    for (m = 1; m <= mseq; m++) {
        int nlev = mglevg(m) + nemgl(m);
        // Fortran: read(iunit5,*) (mit(i,m),i=1,nlev) -- all nlev values in ONE record
        std::vector<int> _mitbuf(nlev);
        fortran_read_list_array(iunit5, _mitbuf.data(), nlev);
        for (int i = 1; i <= nlev; i++) {
            mit(i, m) = _mitbuf[i-1];
        }
        for (int i = 1; i <= nlev; i++) {
            fprintf(fortran_get_unit(iunit11), "%10d", mit(i, m));
        }
        fprintf(fortran_get_unit(iunit11), "\n");
    }

    // nbli - number of 1:1 block interfaces
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    fprintf(fortran_get_unit(iunit11), " 1-1 blocking data:\n");
    fprintf(fortran_get_unit(iunit11), "      nbli\n");
    fortran_read_list(iunit5, &nbli);
    ntest = std::abs(nbli);
    ntest = ntest + ncgmax * ntest;
    if (ntest > mxbli) {
        fprintf(fortran_get_unit(iunit11),
            " too many block boundaries specified - increase mxbli\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    fprintf(fortran_get_unit(iunit11), "%10d\n", nbli);
    if (std::abs(nbli) > 0) {
        for (n = 1; n <= std::abs(nbli); n++) nblon(n) = 0;
    }

    // first block of 1:1 interface data
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    fprintf(fortran_get_unit(iunit11),
        "  number    grid    ista    jsta    ksta    iend    jend    kend   isva1   isva2\n");
    if (std::abs(nbli) > 0) {
        for (n = 1; n <= std::abs(nbli); n++) {
            int mdum;
            fortran_read_list(iunit5, &mdum, &nblk(1,n),
                              &limblk(1,1,n), &limblk(1,2,n), &limblk(1,3,n),
                              &limblk(1,4,n), &limblk(1,5,n), &limblk(1,6,n),
                              &isva(1,1,n), &isva(1,2,n));
            fprintf(fortran_get_unit(iunit11), "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
                    mdum, nblk(1,n),
                    limblk(1,1,n), limblk(1,2,n), limblk(1,3,n),
                    limblk(1,4,n), limblk(1,5,n), limblk(1,6,n),
                    isva(1,1,n), isva(1,2,n));
            // preliminary 1-1 interface check
            for (m = 1; m <= 3; m++) {
                if (m + isva(1,1,n) + isva(1,2,n) == 6) {
                    if (limblk(1,m,n) != limblk(1,m+3,n)) {
                        const char* dir = (m==1)?"i":(m==2)?"j":"k";
                        fprintf(fortran_get_unit(iunit11),
                            " 1-1 interface error: %s-index not constant on interface\n", dir);
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }
    }

    // second block of 1:1 interface data
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    fprintf(fortran_get_unit(iunit11),
        "  number    grid    ista    jsta    ksta    iend    jend    kend   isva1   isva2\n");
    if (std::abs(nbli) > 0) {
        for (n = 1; n <= std::abs(nbli); n++) {
            int mdum;
            fortran_read_list(iunit5, &mdum, &nblk(2,n),
                              &limblk(2,1,n), &limblk(2,2,n), &limblk(2,3,n),
                              &limblk(2,4,n), &limblk(2,5,n), &limblk(2,6,n),
                              &isva(2,1,n), &isva(2,2,n));
            fprintf(fortran_get_unit(iunit11), "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
                    mdum, nblk(2,n),
                    limblk(2,1,n), limblk(2,2,n), limblk(2,3,n),
                    limblk(2,4,n), limblk(2,5,n), limblk(2,6,n),
                    isva(2,1,n), isva(2,2,n));
            // preliminary 1-1 interface check
            for (m = 1; m <= 3; m++) {
                if (m + isva(2,1,n) + isva(2,2,n) == 6) {
                    if (limblk(2,m,n) != limblk(2,m+3,n)) {
                        const char* dir = (m==1)?"i":(m==2)?"j":"k";
                        fprintf(fortran_get_unit(iunit11),
                            " 1-1 interface error: %s-index not constant on interface\n", dir);
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }
    }


    // ninter - patch interface data
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    fprintf(fortran_get_unit(iunit11), " patch interface data:\n");
    fprintf(fortran_get_unit(iunit11), "    ninter\n");
    fortran_read_list(iunit5, &ninter);
    fprintf(fortran_get_unit(iunit11), "%10d\n", ninter);

    if (iunst == 0 && ninter > 0) {
        fprintf(fortran_get_unit(iunit11),
            " coefficients for patching must be precalculated and available in file patch.bin\n");
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    if (ninter == 0) {
        lnsub1  = 1;
        lintmax = 1;
        lmaxxe  = 1;
    } else {
        rpatch0_ns::rpatch0(intmax, nsub1, iindex, ninter);
        rewind(fortran_get_unit(22));
    }

    // set level and grid for each block
    lfem = 0;
    lfgm = ncgmax + 1;
    lcgm = 1;
    for (igrid = 1; igrid <= ngrid; igrid++) {
        iem         = iemg(igrid);
        nbl         = nblg(igrid);
        ncg         = ncgg(igrid);
        igridg(nbl) = igrid;
        levelg(nbl) = iem + lfgm;
        if (iem == 0) lfgm = std::max(lfgm, levelg(nbl));
        if (iem > 0)  lfem = std::max(lfem, levelg(nbl));
        if (ncg > 0 && iem == 0) {
            for (n = 1; n <= ncg; n++) {
                nbl         = nbl + 1;
                igridg(nbl) = igrid;
                levelg(nbl) = levelg(nbl-1) - 1;
            }
        }
    }

    // determine highest levels
    int levtt = levelt[mseq-1];
    int levgg = levtt;
    if (nemgl(mseq) > 0) levgg = levtt - nemgl(mseq);
    int msmax = 0;
    for (int ms = 1; ms <= mseq; ms++) {
        if (ncyc1[ms-1] > 0) msmax = ms;
    }
    int levtop = levelt[msmax-1];
    int levglb = levtop;
    if (nemgl(msmax) > 0) levglb = levtop - nemgl(msmax);

    int nplt3dtmp_local = 0;
    // plot3d output data
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }

    if (nplot3d != 0) {
        iflag    = 0;
        ip3dsurf = 0;
        if (nplot3d < 0) {
            iflag = 1;
            if (i2d == 0) ip3dsurf = 1;
        }
        fprintf(fortran_get_unit(iunit11), " plot3d output:\n");
        fprintf(fortran_get_unit(iunit11),
            "  grid block iptyp  ista  iend  iinc  jsta  jend  jinc  ksta  kend  kinc\n");

        int np3dhold = std::abs(nplot3d);
        int nplt3dtmp_local = 0;
        if (nplot3d < 0) {
            // count from ipl3dtmp
            for (int nn = 1; nn <= np3dhold; nn++) {
                if (ipl3dtmp(1,nn) > 0) nplt3dtmp_local++;
            }
            nplot3d = nplt3dtmp_local;
        }
        int np3dread = std::max(nplot3d, np3dhold);
        n1 = 0;
        int iptyp = 0;
        for (n = 1; n <= np3dread; n++) {
            if (iflag == 0) {
                fortran_read_list(iunit5, &im[1], &im[2], &im[3], &im[4], &im[5],
                                  &im[6], &im[7], &im[8], &im[9], &im[10], &im[11]);
                if (im[2] == -1) {
                    fprintf(fortran_get_unit(iunit11),
                        " stopping...iptyp should not be -1 unless nplot3d=-1 (surface only)\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (im[2] == 7) {
                    fprintf(fortran_get_unit(iunit11),
                        " stopping...iptyp should not be 7 unless nplot3d=-1 (surface only)\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (ivmx == 72 && im[2] == 2) {
                    fprintf(fortran_get_unit(iunit11),
                        " stopping...iptyp=2 not implemented for ivisc=72\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            } else {
                if (n <= np3dhold) {
                    int idum_tmp;
                    fortran_read_list(iunit5, &idum_tmp, &iptyp,
                                      &idum_tmp, &idum_tmp, &idum_tmp,
                                      &idum_tmp, &idum_tmp, &idum_tmp,
                                      &idum_tmp, &idum_tmp, &idum_tmp);
                }
                if (iptyp == 7 && i2d != 0) {
                    fprintf(fortran_get_unit(iunit11),
                        " stopping...iptyp 7 does not work for 2D\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (ivmx == 72 && iptyp == 2) {
                    fprintf(fortran_get_unit(iunit11),
                        " stopping...iptyp=2 not implemented for ivisc=72\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }
            if (n <= nplot3d) {
                if (iflag != 0) {
                    im[1] = ipl3dtmp(1, n);
                    im[2] = iptyp;
                    for (l = 3; l <= 11; l++) im[l] = ipl3dtmp(l, n);
                    // fixup range if cell center
                    if (iptyp > 0) {
                        int idimm = idimg(nblg(im[1]));
                        int jdimm = jdimg(nblg(im[1]));
                        int kdimm = kdimg(nblg(im[1]));
                        if (im[4] > im[3]) {
                            im[4] = im[4] - 1;
                        } else if (im[3] == idimm) {
                            im[3] = idimm - 1; im[4] = idimm - 1;
                        }
                        if (im[7] > im[6]) {
                            im[7] = im[7] - 1;
                        } else if (im[6] == jdimm) {
                            im[6] = jdimm - 1; im[7] = jdimm - 1;
                        }
                        if (im[10] > im[9]) {
                            im[10] = im[10] - 1;
                        } else if (im[9] == kdimm) {
                            im[9] = kdimm - 1; im[10] = kdimm - 1;
                        }
                    }
                }


                // check for inappropriate iptyp
                iptyp = im[2];
                if ((iptyp-0)*(iptyp-1)*(iptyp-1)*(iptyp-3)*
                    (iptyp-4)*(iptyp-5)*(iptyp+5)*(iptyp-6)*
                    (iptyp+6)*(iptyp-2)*(iptyp+4)*(iptyp+1)*
                    (iptyp-7) != 0) {
                    fprintf(fortran_get_unit(iunit11),
                        " stopping...inappropriate value for iptyp for plot3d data\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (ifort50write != 0 && (ngrid > 1 || iptyp != 2)) {
                    fprintf(fortran_get_unit(iunit11),
                        " stopping...ifort50write must have ngrid=+-1 and iptyp=2\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                igrid = im[1];
                nbl   = nblg(igrid);
                int jd = jdimg(nbl);
                int kd = kdimg(nbl);
                int id = idimg(nbl);
                // set entire index range if all zeros
                if (im[3] == 0 && im[4] == 0) { im[3]=1; im[4]=id; im[5]=1; }
                if (im[6] == 0 && im[7] == 0) { im[6]=1; im[7]=jd; im[8]=1; }
                if (im[9] == 0 && im[10] == 0) { im[9]=1; im[10]=kd; im[11]=1; }
                // find block number at highest level
                nbl = 0;
                int iflg = 0;
                for (int nnn = 1; nnn <= nblock; nnn++) {
                    if (levelg(nnn) == levglb) {
                        if (igridg(nnn) == igrid) nbl = nnn;
                    }
                    if (iemg(igrid) > 0) {
                        if (levelg(nnn) > levgg) {
                            if (igridg(nnn) == igrid) nbl = nnn;
                            if (levtop < levtt) iflg = 1;
                        }
                    }
                }
                if (nbl == 0) {
                    fprintf(fortran_get_unit(iunit11),
                        "   stopping: no block corresponds to grid %3d at the highest level computed on\n", igrid);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (iflg == 0) {
                    n1++;
                    for (l = 1; l <= 11; l++) inpl3d(n1, l) = im[l];
                    inpl3d(n1, 1) = nbl;
                    jd = jdimg(nbl); kd = kdimg(nbl); id = idimg(nbl);
                    if (levglb < levgg) {
                        for (int mm = 3; mm <= 11; mm++) {
                            if (mm == 4 && inpl3d(n1,mm) == 2) continue;
                            inpl3d(n1,mm) = (inpl3d(n1,mm)-1) / (1 << (levgg-levglb)) + 1;
                        }
                    }
                    if (inpl3d(n1,3)  <= 0)  inpl3d(n1,3)  = 1;
                    if (inpl3d(n1,3)  > id)  inpl3d(n1,3)  = id;
                    if (inpl3d(n1,4)  <= 0)  inpl3d(n1,4)  = id;
                    if (inpl3d(n1,4)  > id)  inpl3d(n1,4)  = id;
                    if (inpl3d(n1,5)  == 0)  inpl3d(n1,5)  = 1;
                    if (inpl3d(n1,6)  <= 0)  inpl3d(n1,6)  = 1;
                    if (inpl3d(n1,6)  > jd)  inpl3d(n1,6)  = jd;
                    if (inpl3d(n1,7)  <= 0)  inpl3d(n1,7)  = jd;
                    if (inpl3d(n1,7)  > jd)  inpl3d(n1,7)  = jd;
                    if (inpl3d(n1,8)  == 0)  inpl3d(n1,8)  = 1;
                    if (inpl3d(n1,9)  <= 0)  inpl3d(n1,9)  = 1;
                    if (inpl3d(n1,9)  > kd)  inpl3d(n1,9)  = kd;
                    if (inpl3d(n1,10) <= 0)  inpl3d(n1,10) = kd;
                    if (inpl3d(n1,10) > kd)  inpl3d(n1,10) = kd;
                    if (inpl3d(n1,11) == 0)  inpl3d(n1,11) = 1;
                    if (i2d == 1) {
                        inpl3d(n1,3) = 1; inpl3d(n1,4) = 1; inpl3d(n1,5) = 1;
                    }
                    fprintf(fortran_get_unit(iunit11), "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
                            igrid, inpl3d(n1,1), inpl3d(n1,2), inpl3d(n1,3),
                            inpl3d(n1,4), inpl3d(n1,5), inpl3d(n1,6), inpl3d(n1,7),
                            inpl3d(n1,8), inpl3d(n1,9), inpl3d(n1,10), inpl3d(n1,11));
                } else {
                    fprintf(fortran_get_unit(iunit11),
                        "    stopping: embedded grids are computed only on finest level\n"
                        "    the way you are mesh sequencing, can't output solution on grid%3d\n", igrid);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            } // if n <= nplot3d
        } // end do 136

        // check all zones have same iptype if any > 2
        int ifunc = 0;
        if (nplot3d > 0) {
            for (int nnn = 1; nnn <= nplot3d; nnn++) {
                if (inpl3d(nnn,2) > 2) ifunc = inpl3d(nnn,2);
            }
            if (ifunc > 0) {
                for (int nnn = 1; nnn <= nplot3d; nnn++) {
                    if (inpl3d(nnn,2) != ifunc) {
                        fprintf(fortran_get_unit(iunit11),
                            " stopping...must have iptype same for all zones if any zone has iptype > 2\n");
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }
    } // end if nplot3d != 0

    // movie
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    fortran_read_list(iunit5, &movie);
    fprintf(fortran_get_unit(iunit11), " movie\n");
    fprintf(fortran_get_unit(iunit11), "%6d\n", movie);


    // printed output data
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }

    if (nprint != 0) {
        iflag    = 0;
        iprnsurf = 0;
        if (nprint < 0) {
            iflag = 1;
            if (i2d == 0) iprnsurf = 1;
        }
        fprintf(fortran_get_unit(iunit11), " print out:\n");
        fprintf(fortran_get_unit(iunit11),
            "  grid block iptyp  ista  iend  iinc  jsta  jend  jinc  ksta  kend  kinc\n");

        int nprnthold = 0;
        if (nprint < 0) {
            nprnthold = std::abs(nprint);
            nprint = nplt3dtmp_local;
        }
        n1 = 0;
        int iptyp_pr = 0;
        for (n = 1; n <= nprint; n++) {
            if (iflag == 0) {
                fortran_read_list(iunit5, &im[1], &im[2], &im[3], &im[4], &im[5],
                                  &im[6], &im[7], &im[8], &im[9], &im[10], &im[11]);
            } else {
                if (n <= nprnthold) {
                    int idum_tmp;
                    fortran_read_list(iunit5, &idum_tmp, &iptyp_pr,
                                      &idum_tmp, &idum_tmp, &idum_tmp,
                                      &idum_tmp, &idum_tmp, &idum_tmp,
                                      &idum_tmp, &idum_tmp, &idum_tmp);
                }
            }
            if (im[5] < 0 || im[8] < 0 || im[11] < 0) {
                fprintf(fortran_get_unit(iunit11), " input error: inc < 0\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (im[3] > im[4] || im[6] > im[7] || im[9] > im[10]) {
                fprintf(fortran_get_unit(iunit11), " input error: end > start\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (iflag != 0) {
                im[1] = ipl3dtmp(1, n);
                im[2] = iptyp_pr;
                for (l = 3; l <= 11; l++) im[l] = ipl3dtmp(l, n);
                if (iptyp_pr > 0) {
                    int idimm = idimg(nblg(im[1]));
                    int jdimm = jdimg(nblg(im[1]));
                    int kdimm = kdimg(nblg(im[1]));
                    if (im[4] > im[3]) { im[4] = im[4]-1; }
                    else if (im[3] == idimm) { im[3] = idimm-1; im[4] = idimm-1; }
                    if (im[7] > im[6]) { im[7] = im[7]-1; }
                    else if (im[6] == jdimm) { im[6] = jdimm-1; im[7] = jdimm-1; }
                    if (im[10] > im[9]) { im[10] = im[10]-1; }
                    else if (im[9] == kdimm) { im[9] = kdimm-1; im[10] = kdimm-1; }
                }
            }
            iptyp_pr = im[2];
            if ((iptyp_pr-0)*(iptyp_pr-1) != 0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...inappropriate value for iptyp for printout data\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            igrid = im[1];
            nbl   = nblg(igrid);
            int jd = jdimg(nbl); int kd = kdimg(nbl); int id = idimg(nbl);
            if (im[3] == 0 && im[4] == 0) { im[3]=1; im[4]=id; im[5]=1; }
            if (im[6] == 0 && im[7] == 0) { im[6]=1; im[7]=jd; im[8]=1; }
            if (im[9] == 0 && im[10] == 0) { im[9]=1; im[10]=kd; im[11]=1; }
            nbl = 0;
            int iflg = 0;
            for (int nnn = 1; nnn <= nblock; nnn++) {
                if (levelg(nnn) == levglb) {
                    if (igridg(nnn) == igrid) nbl = nnn;
                }
                if (iemg(igrid) > 0) {
                    if (levelg(nnn) > levgg) {
                        if (igridg(nnn) == igrid) nbl = nnn;
                        if (levtop < levtt) iflg = 1;
                    }
                }
            }
            if (nbl == 0) {
                fprintf(fortran_get_unit(iunit11),
                    "   stopping: no block corresponds to grid %3d at the highest level computed on\n", igrid);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (iflg == 0) {
                n1++;
                for (l = 1; l <= 11; l++) inpr(n1, l) = im[l];
                inpr(n1, 1) = nbl;
                jd = jdimg(nbl); kd = kdimg(nbl); id = idimg(nbl);
                if (nemgl(mseq) > 0) levgg = levtt - nemgl(mseq);
                if (levglb < levgg) {
                    for (int mm = 3; mm <= 11; mm++) {
                        if (mm == 4 && inpr(n1,mm) == 2) continue;
                        inpr(n1,mm) = (inpr(n1,mm)-1) / (1 << (levgg-levglb)) + 1;
                    }
                }
                if (inpr(n1,3)  <= 0)  inpr(n1,3)  = 1;
                if (inpr(n1,3)  > id)  inpr(n1,3)  = id;
                if (inpr(n1,4)  <= 0)  inpr(n1,4)  = id;
                if (inpr(n1,4)  > id)  inpr(n1,4)  = id;
                if (inpr(n1,5)  == 0)  inpr(n1,5)  = 1;
                if (inpr(n1,6)  <= 0)  inpr(n1,6)  = 1;
                if (inpr(n1,6)  > jd)  inpr(n1,6)  = jd;
                if (inpr(n1,7)  <= 0)  inpr(n1,7)  = jd;
                if (inpr(n1,7)  > jd)  inpr(n1,7)  = jd;
                if (inpr(n1,8)  == 0)  inpr(n1,8)  = 1;
                if (inpr(n1,9)  <= 0)  inpr(n1,9)  = 1;
                if (inpr(n1,9)  > kd)  inpr(n1,9)  = kd;
                if (inpr(n1,10) <= 0)  inpr(n1,10) = kd;
                if (inpr(n1,10) > kd)  inpr(n1,10) = kd;
                if (inpr(n1,11) == 0)  inpr(n1,11) = 1;
                if (i2d == 1) { inpr(n1,3)=1; inpr(n1,4)=1; inpr(n1,5)=1; }
                fprintf(fortran_get_unit(iunit11), "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
                        igrid, inpr(n1,1), inpr(n1,2), inpr(n1,3),
                        inpr(n1,4), inpr(n1,5), inpr(n1,6), inpr(n1,7),
                        inpr(n1,8), inpr(n1,9), inpr(n1,10), inpr(n1,11));
            } else {
                fprintf(fortran_get_unit(iunit11),
                    "    stopping: embedded grids are computed only on finest level\n"
                    "    the way you are mesh sequencing, can't output solution on grid%3d\n", igrid);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        } // end do 236
    } // end if nprint != 0


    // control surface data
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
    fortran_read_list(iunit5, &ncs);
    { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }

    ncs2 = 0;
    if (ncs != 0) {
        fprintf(fortran_get_unit(iunit11), " control surfaces:\n");
        fprintf(fortran_get_unit(iunit11), "   ncs\n");
        fprintf(fortran_get_unit(iunit11), "%6d\n", ncs);
        fprintf(fortran_get_unit(iunit11),
            "  grid block  ista  iend  jsta  jend  ksta  kend iwall inorm\n");

        int ncs1 = ncs;
        ncs = 0;
        int ntestcs = 0;
        for (n = 1; n <= ncs1; n++) {
            fortran_read_list(iunit5, &im[1], &im[2], &im[3], &im[4], &im[5],
                              &im[6], &im[7], &im[8], &im[9]);
            igrid = im[1];
            nbl   = nblg(igrid);
            int jd = jdimg(nbl); int kd = kdimg(nbl); int id = idimg(nbl);
            if (im[2] == 0 && im[3] == 0) { im[2]=1; im[3]=id; }
            if (im[4] == 0 && im[5] == 0) { im[4]=1; im[5]=jd; }
            if (im[6] == 0 && im[7] == 0) { im[6]=1; im[7]=kd; }
            ntestcs = ntestcs + (ncgg(igrid) + 1);
            if (icall > 0) {
                if (ntestcs > maxcs - 1) {
                    fprintf(fortran_get_unit(iunit11),
                        " stopping...to accomodate fine and coarse levels, for an input ncs = %d\n", ncs1);
                    fprintf(fortran_get_unit(iunit11),
                        " the parameter maxcs must be at least %d\n",
                        ncs1*(ncgg(igrid)+1) + 1);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }
            for (ncg = 1; ncg <= ncgg(igrid)+1; ncg++) {
                nbl  = nblg(igrid) + ncg - 1;
                ncs  = ncs + 1;
                ncs2 = ncs2 + 1;
                if (icall == 0) ncs = std::min(maxcs, ncs);
                icsinfo(ncs, 1) = nbl;
                for (int mm = 2; mm <= 7; mm++) {
                    icsinfo(ncs, mm) = (im[mm]-1) / (1 << (ncg-1)) + 1;
                }
                icsinfo(ncs, 8) = im[8];
                icsinfo(ncs, 9) = im[9];
                jd = jdimg(nbl); kd = kdimg(nbl); id = idimg(nbl);
                if (icsinfo(ncs,2) <= 0)  icsinfo(ncs,2) = 1;
                if (icsinfo(ncs,2) > id)  icsinfo(ncs,2) = id;
                if (icsinfo(ncs,3) <= 0)  icsinfo(ncs,3) = id;
                if (icsinfo(ncs,3) > id)  icsinfo(ncs,3) = id;
                if (icsinfo(ncs,4) <= 0)  icsinfo(ncs,4) = 1;
                if (icsinfo(ncs,4) > jd)  icsinfo(ncs,4) = jd;
                if (icsinfo(ncs,5) <= 0)  icsinfo(ncs,5) = jd;
                if (icsinfo(ncs,5) > jd)  icsinfo(ncs,5) = jd;
                if (icsinfo(ncs,6) <= 0)  icsinfo(ncs,6) = 1;
                if (icsinfo(ncs,6) > kd)  icsinfo(ncs,6) = kd;
                if (icsinfo(ncs,7) <= 0)  icsinfo(ncs,7) = kd;
                if (icsinfo(ncs,7) > kd)  icsinfo(ncs,7) = kd;
                if (i2d == 1) { icsinfo(ncs,2)=1; icsinfo(ncs,3)=2; }
                if (ncg == 1) {
                    fprintf(fortran_get_unit(iunit11), "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
                            igrid, icsinfo(ncs,1), icsinfo(ncs,2), icsinfo(ncs,3),
                            icsinfo(ncs,4), icsinfo(ncs,5), icsinfo(ncs,6),
                            icsinfo(ncs,7), icsinfo(ncs,8), icsinfo(ncs,9));
                }
            }
        }
    } // end if ncs != 0


    // rigid grid motion - set defaults
    itransmc = 0; irotatmc = 0;
    rfreqtmc = 0.0f; rfreqrmc = 0.0f;
    utransmc = 0.0f; vtransmc = 0.0f; wtransmc = 0.0f;
    omegaxmc = 0.0f; omegaymc = 0.0f; omegazmc = 0.0f;
    thetaxmc = 0.0f; thetaymc = 0.0f; thetazmc = 0.0f;
    xorigmc  = 0.0f; yorigmc  = 0.0f; zorigmc  = 0.0f;
    xorig0mc = 0.0f; yorig0mc = 0.0f; zorig0mc = 0.0f;
    dxmxmc   = 0.0f; dymxmc   = 0.0f; dzmxmc   = 0.0f;
    dthxmxmc = 0.0f; dthymxmc = 0.0f; dthzmxmc = 0.0f;
    time2mc  = 0.0f;
    xmc0 = 0.0f; ymc0 = 0.0f; zmc0 = 0.0f;

    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nblg(igrid);
        ncg = ncgg(igrid);
        itrans(nbl)  = 0;  utrans(nbl)  = 0.; vtrans(nbl)  = 0.; wtrans(nbl)  = 0.;
        dxmx(nbl)    = 0.; dymx(nbl)    = 0.; dzmx(nbl)    = 0.;
        irotat(nbl)  = 0;  omegax(nbl)  = 0.; omegay(nbl)  = 0.; omegaz(nbl)  = 0.;
        xorig(nbl)   = 0.; yorig(nbl)   = 0.; zorig(nbl)   = 0.;
        xorig0(nbl)  = 0.; yorig0(nbl)  = 0.; zorig0(nbl)  = 0.;
        thetax(nbl)  = 0.; thetay(nbl)  = 0.; thetaz(nbl)  = 0.;
        dthxmx(nbl)  = 0.; dthymx(nbl)  = 0.; dthzmx(nbl)  = 0.;
        rfreqt(nbl)  = 0.; rfreqr(nbl)  = 0.; time2(nbl)   = 0.;
        thetaxl(nbl) = 0.; thetayl(nbl) = 0.; thetazl(nbl) = 0.;
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl++;
                itrans(nbl)  = itrans(nbl-1);  utrans(nbl)  = utrans(nbl-1);
                vtrans(nbl)  = vtrans(nbl-1);   wtrans(nbl)  = wtrans(nbl-1);
                irotat(nbl)  = irotat(nbl-1);   omegax(nbl)  = omegax(nbl-1);
                omegay(nbl)  = omegay(nbl-1);   omegaz(nbl)  = omegaz(nbl-1);
                xorig(nbl)   = xorig(nbl-1);    yorig(nbl)   = yorig(nbl-1);
                zorig(nbl)   = zorig(nbl-1);    xorig0(nbl)  = xorig0(nbl-1);
                yorig0(nbl)  = yorig0(nbl-1);   zorig0(nbl)  = zorig0(nbl-1);
                thetax(nbl)  = thetax(nbl-1);   thetay(nbl)  = thetay(nbl-1);
                thetaz(nbl)  = thetaz(nbl-1);   rfreqt(nbl)  = rfreqt(nbl-1);
                rfreqr(nbl)  = rfreqr(nbl-1);   dxmx(nbl)    = dxmx(nbl-1);
                dymx(nbl)    = dymx(nbl-1);     dzmx(nbl)    = dzmx(nbl-1);
                dthxmx(nbl)  = dthxmx(nbl-1);  dthymx(nbl)  = dthymx(nbl-1);
                dthzmx(nbl)  = dthzmx(nbl-1);  time2(nbl)   = time2(nbl-1);
                thetaxl(nbl) = thetaxl(nbl-1); thetayl(nbl) = thetayl(nbl-1);
                thetazl(nbl) = thetazl(nbl-1);
            }
        }
    }

    if (iunst == 1 || iunst == 3) {
        // translation data
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        int ntrans;
        fortran_read_list(iunit5, &ntrans);
        fprintf(fortran_get_unit(iunit11), " translation data:\n");
        fprintf(fortran_get_unit(iunit11), " ntrans\n");
        fprintf(fortran_get_unit(iunit11), "%6d\n", ntrans);
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        fprintf(fortran_get_unit(iunit11), " lreftra\n");
        double lreftra = 1.0;
        if (ntrans > 0) {
            fortran_read_list(iunit5, &realval[1]);
            lreftra = realval[1];
            fprintf(fortran_get_unit(iunit11), "%g\n", (float)lreftra);
        }
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        fprintf(fortran_get_unit(iunit11), " grid itrans rfreq utrans vtrans wtrans\n");
        if (ntrans > 0) {
            int ichk1 = 0;
            if (ntrans > ngrid) {
                if (ntrans == ngrid + 1) {
                    ichk1 = 1;
                } else {
                    fprintf(fortran_get_unit(iunit11),
                        " stopping...must have ntrans <= ngrid + 1\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }
            for (igrid = 1; igrid <= ntrans; igrid++) {
                int igr, itran;
                fortran_read_list(iunit5, &igr, &itran,
                                  &realval[1], &realval[2], &realval[3], &realval[4]);
                double rfreq = realval[1], utran = realval[2];
                double vtran = realval[3], wtran = realval[4];
                if (igr > 0) {
                    fprintf(fortran_get_unit(iunit11), "%6d%6d%12.5f%12.5f%12.5f%12.5f\n",
                            igr, itran, (float)rfreq, (float)utran, (float)vtran, (float)wtran);
                    nbl = nblg(igr); ncg = ncgg(igr);
                    itrans(nbl) = itran;
                    utrans(nbl) = utran; vtrans(nbl) = vtran; wtrans(nbl) = wtran;
                    if (itran == 1) lreftra = 1.0;
                    rfreqt(nbl) = 2.0 * pi * rfreq / lreftra;
                    if (ncg > 0) {
                        for (n = 1; n <= ncg; n++) {
                            nbl++;
                            itrans(nbl) = itrans(nbl-1);
                            utrans(nbl) = utrans(nbl-1); vtrans(nbl) = vtrans(nbl-1);
                            wtrans(nbl) = wtrans(nbl-1); rfreqt(nbl) = rfreqt(nbl-1);
                        }
                    }
                }
                if (igr == 0) {
                    fprintf(fortran_get_unit(iunit11), " mc: %6d%12.5f%12.5f%12.5f%12.5f\n",
                            itran, (float)rfreq, (float)utran, (float)vtran, (float)wtran);
                    itransmc = itran;
                    utransmc = (float)utran; vtransmc = (float)vtran; wtransmc = (float)wtran;
                    rfreqtmc = (float)(2.0 * pi * rfreq / lreftra);
                    ichk1 = 0;
                }
            }
            if (ichk1 != 0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...since ntrans = ngrid + 1, must have grid = 0 (moment center) in\n"
                    " one input line of translation data\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        }


        // dxmx/dymx/dzmx data
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        fprintf(fortran_get_unit(iunit11), "    grid    dxmx    dymx    dzmx\n");
        if (ntrans > 0) {
            for (igrid = 1; igrid <= ntrans; igrid++) {
                int igr;
                fortran_read_list(iunit5, &igr, &realval[1], &realval[2], &realval[3]);
                double dxmax = realval[1], dymax = realval[2], dzmax = realval[3];
                if (igr > 0) {
                    fprintf(fortran_get_unit(iunit11), "     %3d %8.4f %8.4f %8.4f\n",
                            igr, (float)dxmax, (float)dymax, (float)dzmax);
                    nbl = nblg(igr); ncg = ncgg(igr);
                    dxmx(nbl) = dxmax; dymx(nbl) = dymax; dzmx(nbl) = dzmax;
                    if (ncg > 0) {
                        for (n = 1; n <= ncg; n++) {
                            nbl++;
                            dxmx(nbl) = dxmx(nbl-1); dymx(nbl) = dymx(nbl-1); dzmx(nbl) = dzmx(nbl-1);
                        }
                    }
                }
                if (igr == 0) {
                    fprintf(fortran_get_unit(iunit11), "      MC %8.4f %8.4f %8.4f\n",
                            (float)dxmax, (float)dymax, (float)dzmax);
                    dxmxmc = (float)dxmax; dymxmc = (float)dymax; dzmxmc = (float)dzmax;
                }
            }
        }

        // rotation data
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        int nrotat;
        fortran_read_list(iunit5, &nrotat);
        fprintf(fortran_get_unit(iunit11), " dynamic grid data - rotation:\n");
        fprintf(fortran_get_unit(iunit11), "   nrotat\n");
        fprintf(fortran_get_unit(iunit11), "     %3d\n", nrotat);
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        fprintf(fortran_get_unit(iunit11), "     lref\n");
        double lrefrot = 1.0;
        if (nrotat > 0) {
            fortran_read_list(iunit5, &realval[1]);
            lrefrot = realval[1];
            fprintf(fortran_get_unit(iunit11), "%8.4f\n", (float)lrefrot);
        }
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        fprintf(fortran_get_unit(iunit11),
            "    grid  irotat  rfreqr  omegax  omegay  omegaz   xorig   yorig   zorig\n");
        if (nrotat > 0) {
            int ichk1 = 0;
            if (nrotat > ngrid) {
                if (nrotat == ngrid + 1) { ichk1 = 1; }
                else {
                    fprintf(fortran_get_unit(iunit11),
                        " stopping...must have nrotat <= ngrid + 1\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }
            for (igrid = 1; igrid <= nrotat; igrid++) {
                int igr, irota;
                fortran_read_list(iunit5, &igr, &irota,
                                  &realval[1], &realval[2], &realval[3], &realval[4],
                                  &realval[5], &realval[6], &realval[7]);
                double rfreq = realval[1], omegx = realval[2], omegy = realval[3];
                double omegz = realval[4], xorg = realval[5], yorg = realval[6], zorg = realval[7];
                if (igr > 0) {
                    fprintf(fortran_get_unit(iunit11),
                            "     %3d       %2d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                            igr, irota, (float)rfreq, (float)omegx, (float)omegy,
                            (float)omegz, (float)xorg, (float)yorg, (float)zorg);
                    nbl = nblg(igr); ncg = ncgg(igr);
                    irotat(nbl) = irota;
                    xorig(nbl) = xorg; yorig(nbl) = yorg; zorig(nbl) = zorg;
                    if (irota == 1) {
                        omegax(nbl) = 2.0*pi*omegx/lrefrot;
                        omegay(nbl) = 2.0*pi*omegy/lrefrot;
                        omegaz(nbl) = 2.0*pi*omegz/lrefrot;
                    }
                    if (irota == 2 || irota == 3) {
                        omegax(nbl) = omegx/radtodeg;
                        omegay(nbl) = omegy/radtodeg;
                        omegaz(nbl) = omegz/radtodeg;
                    }
                    rfreqr(nbl) = 2.0*pi*rfreq/lrefrot;
                    if (ncg > 0) {
                        for (n = 1; n <= ncg; n++) {
                            nbl++;
                            irotat(nbl) = irotat(nbl-1); omegax(nbl) = omegax(nbl-1);
                            omegay(nbl) = omegay(nbl-1); omegaz(nbl) = omegaz(nbl-1);
                            xorig(nbl)  = xorig(nbl-1);  yorig(nbl)  = yorig(nbl-1);
                            zorig(nbl)  = zorig(nbl-1);  thetax(nbl) = thetax(nbl-1);
                            thetay(nbl) = thetay(nbl-1); thetaz(nbl) = thetaz(nbl-1);
                            rfreqr(nbl) = rfreqr(nbl-1);
                        }
                    }
                }
                if (igr == 0) {
                    fprintf(fortran_get_unit(iunit11),
                            "      MC       %2d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                            irota, (float)rfreq, (float)omegx, (float)omegy,
                            (float)omegz, (float)xorg, (float)yorg, (float)zorg);
                    irotatmc = irota;
                    xorigmc = (float)xorg; yorigmc = (float)yorg; zorigmc = (float)zorg;
                    if (irota == 1) {
                        omegaxmc = (float)(2.0*pi*omegx/lrefrot);
                        omegaymc = (float)(2.0*pi*omegy/lrefrot);
                        omegazmc = (float)(2.0*pi*omegz/lrefrot);
                    }
                    if (irota == 2 || irota == 3) {
                        omegaxmc = (float)(omegx/radtodeg);
                        omegaymc = (float)(omegy/radtodeg);
                        omegazmc = (float)(omegz/radtodeg);
                    }
                    rfreqrmc = (float)(2.0*pi*rfreq/lrefrot);
                    ichk1 = 0;
                }
            }
            if (ichk1 != 0) {
                fprintf(fortran_get_unit(iunit11),
                    " stopping...since nrotat = ngrid + 1, must have grid = 0 (moment center) in\n"
                    " one input line of rotation data\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        }
        // dthxmx/dthymx/dthzmx data
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        fprintf(fortran_get_unit(iunit11), "    grid  dthxmx  dthymx  dthzmx\n");
        if (nrotat > 0) {
            for (igrid = 1; igrid <= nrotat; igrid++) {
                int igr;
                fortran_read_list(iunit5, &igr, &realval[1], &realval[2], &realval[3]);
                double dtxmax = realval[1], dtymax = realval[2], dtzmax = realval[3];
                if (igr > 0) {
                    fprintf(fortran_get_unit(iunit11), "     %3d %8.4f %8.4f %8.4f\n",
                            igr, (float)dtxmax, (float)dtymax, (float)dtzmax);
                    nbl = nblg(igr); ncg = ncgg(igr);
                    dthxmx(nbl) = dtxmax/radtodeg;
                    dthymx(nbl) = dtymax/radtodeg;
                    dthzmx(nbl) = dtzmax/radtodeg;
                    if (ncg > 0) {
                        for (n = 1; n <= ncg; n++) {
                            nbl++;
                            dthxmx(nbl) = dthxmx(nbl-1);
                            dthymx(nbl) = dthymx(nbl-1);
                            dthzmx(nbl) = dthzmx(nbl-1);
                        }
                    }
                }
                if (igr == 0) {
                    fprintf(fortran_get_unit(iunit11), "      MC %8.4f %8.4f %8.4f\n",
                            (float)dtxmax, (float)dtymax, (float)dtzmax);
                    dthxmxmc = (float)(dtxmax/radtodeg);
                    dthymxmc = (float)(dtymax/radtodeg);
                    dthzmxmc = (float)(dtzmax/radtodeg);
                }
            }
        }
        // store t=0 values
        xmc0 = (float)xmc; ymc0 = (float)ymc; zmc0 = (float)zmc;
        xorig0mc = xorigmc; yorig0mc = yorigmc; zorig0mc = zorigmc;
        for (n = 1; n <= nblock; n++) {
            xorig0(n) = xorig(n); yorig0(n) = yorig(n); zorig0(n) = zorig(n);
        }
    } // end if iunst == 1 || iunst == 3


    // elastic grid motion defaults
    mxmds    = 1;
    mxaes    = 1;
    mxdefseg = 1;
    nsprgit  = 0;
    naesrf   = 0;
    ndefrm   = 0;

    for (int iaes = 1; iaes <= maxaes; iaes++) {
        for (int ll = 1; ll <= 5; ll++) aesrfdat(ll, iaes) = 0.0;
        for (n = 1; n <= nmds; n++) {
            int nm1 = 2*n-1, nm2 = 2*n;
            freq(n, iaes)   = 0.; gmass(n, iaes) = 0.; damp(n, iaes) = 0.;
            x0(nm1, iaes)   = 0.; x0(nm2, iaes)  = 0.;
            gf0(nm1, iaes)  = 0.; gf0(nm2, iaes)  = 0.;
            for (int ll = 1; ll <= 4; ll++) perturb(n, iaes, ll) = 0.;
        }
    }

    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nblg(igrid);
        ncg = ncgg(igrid);
        idefrm(nbl)   = 0;
        nsegdfrm(nbl) = 0;
        for (int iseg = 1; iseg <= maxsegdg; iseg++) {
            utrnsae(nbl,iseg)=0.; vtrnsae(nbl,iseg)=0.; wtrnsae(nbl,iseg)=0.;
            omgxae(nbl,iseg)=0.;  omgyae(nbl,iseg)=0.;  omgzae(nbl,iseg)=0.;
            xorgae(nbl,iseg)=0.;  yorgae(nbl,iseg)=0.;  zorgae(nbl,iseg)=0.;
            xorgae0(nbl,iseg)=0.; yorgae0(nbl,iseg)=0.; zorgae0(nbl,iseg)=0.;
            icouple(nbl,iseg)=0;
            thtxae(nbl,iseg)=0.;  thtyae(nbl,iseg)=0.;  thtzae(nbl,iseg)=0.;
            rfrqtae(nbl,iseg)=0.; rfrqrae(nbl,iseg)=0.;
            icsi(nbl,iseg)=0; icsf(nbl,iseg)=0;
            jcsi(nbl,iseg)=0; jcsf(nbl,iseg)=0;
            kcsi(nbl,iseg)=0; kcsf(nbl,iseg)=0;
            idfrmseg(nbl,iseg)=0; iaesurf(nbl,iseg)=0;
        }
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl++;
                idefrm(nbl)   = idefrm(nbl-1);
                nsegdfrm(nbl) = nsegdfrm(nbl-1);
                for (int iseg = 1; iseg <= maxsegdg; iseg++) {
                    utrnsae(nbl,iseg)=utrnsae(nbl-1,iseg); vtrnsae(nbl,iseg)=vtrnsae(nbl-1,iseg);
                    wtrnsae(nbl,iseg)=wtrnsae(nbl-1,iseg); omgxae(nbl,iseg)=omgxae(nbl-1,iseg);
                    omgyae(nbl,iseg)=omgyae(nbl-1,iseg);   omgzae(nbl,iseg)=omgzae(nbl-1,iseg);
                    xorgae(nbl,iseg)=xorgae(nbl-1,iseg);   yorgae(nbl,iseg)=yorgae(nbl-1,iseg);
                    zorgae(nbl,iseg)=zorgae(nbl-1,iseg);   xorgae0(nbl,iseg)=xorgae0(nbl-1,iseg);
                    yorgae0(nbl,iseg)=yorgae0(nbl-1,iseg); zorgae0(nbl,iseg)=zorgae0(nbl-1,iseg);
                    icouple(nbl,iseg)=icouple(nbl-1,iseg);
                    thtxae(nbl,iseg)=thtxae(nbl-1,iseg);   thtyae(nbl,iseg)=thtyae(nbl-1,iseg);
                    thtzae(nbl,iseg)=thtzae(nbl-1,iseg);   rfrqtae(nbl,iseg)=rfrqtae(nbl-1,iseg);
                    rfrqrae(nbl,iseg)=rfrqrae(nbl-1,iseg);
                    icsi(nbl,iseg)=icsi(nbl-1,iseg);       icsf(nbl,iseg)=icsf(nbl-1,iseg);
                    jcsi(nbl,iseg)=jcsi(nbl-1,iseg);       jcsf(nbl,iseg)=jcsf(nbl-1,iseg);
                    kcsi(nbl,iseg)=kcsi(nbl-1,iseg);       kcsf(nbl,iseg)=kcsf(nbl-1,iseg);
                    idfrmseg(nbl,iseg)=idfrmseg(nbl-1,iseg); iaesurf(nbl,iseg)=iaesurf(nbl-1,iseg);
                }
            }
        }
    }

    if (iunst > 1 || idef_ss > 0) {
        riskp.allocate(maxbl, 500);
        rjskp.allocate(maxbl, 500);
        rkskp.allocate(maxbl, 500);
        iskipt.allocate(maxbl, 500);
        jskipt.allocate(maxbl, 500);
        kskipt.allocate(maxbl, 500);

        // set default skip values
        for (int nb = 1; nb <= maxbl; nb++) {
            for (int ii = 1; ii <= 500; ii++) {
                iskip(nb,ii) = 0; jskip(nb,ii) = 0; kskip(nb,ii) = 0;
            }
        }
        for (int ig = 1; ig <= ngrid; ig++) {
            nbl = nblg(ig);
            int jdim1 = jdimg(nbl)-1, kdim1 = kdimg(nbl)-1, idim1 = idimg(nbl)-1;
            int ntimes = 10;
            int jskip_i0=1, kskip_i0=1, jskip_id=1, kskip_id=1;
            int jrngmin = jdim1, krngmin = kdim1;
            for (int iseg = 1; iseg <= nbci0(nbl); iseg++) {
                int jrange = std::abs(ibcinfo(nbl,iseg,3,1)-ibcinfo(nbl,iseg,2,1));
                int krange = std::abs(ibcinfo(nbl,iseg,5,1)-ibcinfo(nbl,iseg,4,1));
                if (jrange < jrngmin) jrngmin = jrange;
                if (krange < krngmin) krngmin = krange;
                for (int jexp=1; jexp<=ntimes; jexp++) {
                    int power2 = 1<<(jexp-1);
                    float fact = (float)jrange/(float)power2;
                    if ((float)(int)fact == fact) jskip_i0 = power2;
                }
                for (int kexp=1; kexp<=ntimes; kexp++) {
                    int power2 = 1<<(kexp-1);
                    float fact = (float)krange/(float)power2;
                    if ((float)(int)fact == fact) kskip_i0 = power2;
                }
            }
            for (int iseg = 1; iseg <= nbcidim(nbl); iseg++) {
                int jrange = std::abs(ibcinfo(nbl,iseg,3,2)-ibcinfo(nbl,iseg,2,2));
                int krange = std::abs(ibcinfo(nbl,iseg,5,2)-ibcinfo(nbl,iseg,4,2));
                if (jrange < jrngmin) jrngmin = jrange;
                if (krange < krngmin) krngmin = krange;
                for (int jexp=1; jexp<=ntimes; jexp++) {
                    int power2 = 1<<(jexp-1);
                    float fact = (float)jrange/(float)power2;
                    if ((float)(int)fact == fact) jskip_id = power2;
                }
                for (int kexp=1; kexp<=ntimes; kexp++) {
                    int power2 = 1<<(kexp-1);
                    float fact = (float)krange/(float)power2;
                    if ((float)(int)fact == fact) kskip_id = power2;
                }
            }
            int jskip_i, kskip_i;
            if (jdim1/jrngmin*jrngmin == jdim1) jskip_i = jrngmin;
            else jskip_i = std::min(jskip_i0, jskip_id);
            if (kdim1/krngmin*krngmin == kdim1) kskip_i = krngmin;
            else kskip_i = std::min(kskip_i0, kskip_id);


            int iskip_j0=1, kskip_j0=1, iskip_jd=1, kskip_jd=1;
            int irngmin = idim1;
            krngmin = kdim1;
            for (int jseg = 1; jseg <= nbcj0(nbl); jseg++) {
                int irange = std::abs(jbcinfo(nbl,jseg,3,1)-jbcinfo(nbl,jseg,2,1));
                int krange = std::abs(jbcinfo(nbl,jseg,5,1)-jbcinfo(nbl,jseg,4,1));
                if (irange < irngmin) irngmin = irange;
                if (krange < krngmin) krngmin = krange;
                for (int iexp=1; iexp<=ntimes; iexp++) {
                    int power2 = 1<<(iexp-1);
                    float fact = (float)irange/(float)power2;
                    if ((float)(int)fact == fact) iskip_j0 = power2;
                }
                for (int kexp=1; kexp<=ntimes; kexp++) {
                    int power2 = 1<<(kexp-1);
                    float fact = (float)krange/(float)power2;
                    if ((float)(int)fact == fact) kskip_j0 = power2;
                }
            }
            for (int jseg = 1; jseg <= nbcjdim(nbl); jseg++) {
                int irange = std::abs(jbcinfo(nbl,jseg,3,2)-jbcinfo(nbl,jseg,2,2));
                int krange = std::abs(jbcinfo(nbl,jseg,5,2)-jbcinfo(nbl,jseg,4,2));
                if (irange < irngmin) irngmin = irange;
                if (krange < krngmin) krngmin = krange;
                for (int iexp=1; iexp<=ntimes; iexp++) {
                    int power2 = 1<<(iexp-1);
                    float fact = (float)irange/(float)power2;
                    if ((float)(int)fact == fact) iskip_jd = power2;
                }
                for (int kexp=1; kexp<=ntimes; kexp++) {
                    int power2 = 1<<(kexp-1);
                    float fact = (float)krange/(float)power2;
                    if ((float)(int)fact == fact) kskip_jd = power2;
                }
            }
            int iskip_j, kskip_j;
            if (idim1/irngmin*irngmin == idim1) iskip_j = irngmin;
            else iskip_j = std::min(iskip_j0, iskip_jd);
            if (kdim1/krngmin*krngmin == kdim1) kskip_j = krngmin;
            else kskip_j = std::min(kskip_j0, kskip_jd);

            int iskip_k0=1, jskip_k0=1, iskip_kd=1, jskip_kd=1;
            irngmin = idim1; jrngmin = jdim1;
            for (int kseg = 1; kseg <= nbck0(nbl); kseg++) {
                int irange = std::abs(kbcinfo(nbl,kseg,3,1)-kbcinfo(nbl,kseg,2,1));
                int jrange = std::abs(kbcinfo(nbl,kseg,5,1)-kbcinfo(nbl,kseg,4,1));
                if (irange < irngmin) irngmin = irange;
                if (jrange < jrngmin) jrngmin = jrange;
                for (int iexp=1; iexp<=ntimes; iexp++) {
                    int power2 = 1<<(iexp-1);
                    float fact = (float)irange/(float)power2;
                    if ((float)(int)fact == fact) iskip_k0 = power2;
                }
                for (int jexp=1; jexp<=ntimes; jexp++) {
                    int power2 = 1<<(jexp-1);
                    float fact = (float)jrange/(float)power2;
                    if ((float)(int)fact == fact) jskip_k0 = power2;
                }
            }
            for (int kseg = 1; kseg <= nbckdim(nbl); kseg++) {
                int irange = std::abs(kbcinfo(nbl,kseg,3,2)-kbcinfo(nbl,kseg,2,2));
                int jrange = std::abs(kbcinfo(nbl,kseg,5,2)-kbcinfo(nbl,kseg,4,2));
                if (irange < irngmin) irngmin = irange;
                if (jrange < jrngmin) jrngmin = jrange;
                for (int iexp=1; iexp<=ntimes; iexp++) {
                    int power2 = 1<<(iexp-1);
                    float fact = (float)irange/(float)power2;
                    if ((float)(int)fact == fact) iskip_kd = power2;
                }
                for (int jexp=1; jexp<=ntimes; jexp++) {
                    int power2 = 1<<(jexp-1);
                    float fact = (float)jrange/(float)power2;
                    if ((float)(int)fact == fact) jskip_kd = power2;
                }
            }
            int iskip_k, jskip_k;
            if (idim1/irngmin*irngmin == idim1) iskip_k = irngmin;
            else iskip_k = std::min(iskip_k0, iskip_kd);
            if (jdim1/jrngmin*jrngmin == jdim1) jskip_k = jrngmin;
            else jskip_k = std::min(jskip_k0, jskip_kd);

            iskip(nbl,1) = std::min(iskip_j, iskip_k);
            jskip(nbl,1) = std::min(jskip_i, jskip_k);
            kskip(nbl,1) = std::min(kskip_i, kskip_j);
        } // end for ig

        // account for 1-1 interfaces
        if (std::abs(nbli) > 0) {
            for (int it = 1; it <= ngrid; it++) {
                int nskpdif = 0;
                for (n = 1; n <= std::abs(nbli); n++) {
                    int igr1 = nblk(1,n), igr2 = nblk(2,n);
                    int nbl1 = nblg(igr1), nbl2 = nblg(igr2);
                    iskipt(nbl1,1)=iskip(nbl1,1); iskipt(nbl2,1)=iskip(nbl2,1);
                    jskipt(nbl1,1)=jskip(nbl1,1); jskipt(nbl2,1)=jskip(nbl2,1);
                    kskipt(nbl1,1)=kskip(nbl1,1); kskipt(nbl2,1)=kskip(nbl2,1);
                    // isva(1,1,n) vs isva(2,1,n)
                    if (isva(1,1,n)==1) {
                        if (isva(2,1,n)==1) blocking_skip(iskip(nbl1,1),iskip(nbl2,1),iskipt(nbl1,1),iskipt(nbl2,1));
                        else if (isva(2,1,n)==2) blocking_skip(iskip(nbl1,1),jskip(nbl2,1),iskipt(nbl1,1),jskipt(nbl2,1));
                        else if (isva(2,1,n)==3) blocking_skip(iskip(nbl1,1),kskip(nbl2,1),iskipt(nbl1,1),kskipt(nbl2,1));
                    } else if (isva(1,1,n)==2) {
                        if (isva(2,1,n)==1) blocking_skip(jskip(nbl1,1),iskip(nbl2,1),jskipt(nbl1,1),iskipt(nbl2,1));
                        else if (isva(2,1,n)==2) blocking_skip(jskip(nbl1,1),jskip(nbl2,1),jskipt(nbl1,1),jskipt(nbl2,1));
                        else if (isva(2,1,n)==3) blocking_skip(jskip(nbl1,1),kskip(nbl2,1),jskipt(nbl1,1),kskipt(nbl2,1));
                    } else if (isva(1,1,n)==3) {
                        if (isva(2,1,n)==1) blocking_skip(kskip(nbl1,1),iskip(nbl2,1),kskipt(nbl1,1),iskipt(nbl2,1));
                        else if (isva(2,1,n)==2) blocking_skip(kskip(nbl1,1),jskip(nbl2,1),kskipt(nbl1,1),jskipt(nbl2,1));
                        else if (isva(2,1,n)==3) blocking_skip(kskip(nbl1,1),kskip(nbl2,1),kskipt(nbl1,1),kskipt(nbl2,1));
                    }
                    // isva(1,2,n) vs isva(2,2,n)
                    if (isva(1,2,n)==1) {
                        if (isva(2,2,n)==1) blocking_skip(iskip(nbl1,1),iskip(nbl2,1),iskipt(nbl1,1),iskipt(nbl2,1));
                        else if (isva(2,2,n)==2) blocking_skip(iskip(nbl1,1),jskip(nbl2,1),iskipt(nbl1,1),jskipt(nbl2,1));
                        else if (isva(2,2,n)==3) blocking_skip(iskip(nbl1,1),kskip(nbl2,1),iskipt(nbl1,1),kskipt(nbl2,1));
                    } else if (isva(1,2,n)==2) {
                        if (isva(2,2,n)==1) blocking_skip(jskip(nbl1,1),iskip(nbl2,1),jskipt(nbl1,1),iskipt(nbl2,1));
                        else if (isva(2,2,n)==2) blocking_skip(jskip(nbl1,1),jskip(nbl2,1),jskipt(nbl1,1),jskipt(nbl2,1));
                        else if (isva(2,2,n)==3) blocking_skip(jskip(nbl1,1),kskip(nbl2,1),jskipt(nbl1,1),kskipt(nbl2,1));
                    } else if (isva(1,2,n)==3) {
                        if (isva(2,2,n)==1) blocking_skip(kskip(nbl1,1),iskip(nbl2,1),kskipt(nbl1,1),iskipt(nbl2,1));
                        else if (isva(2,2,n)==2) blocking_skip(kskip(nbl1,1),jskip(nbl2,1),kskipt(nbl1,1),jskipt(nbl2,1));
                        else if (isva(2,2,n)==3) blocking_skip(kskip(nbl1,1),kskip(nbl2,1),kskipt(nbl1,1),kskipt(nbl2,1));
                    }
                    nskpdif += std::abs(iskip(nbl1,1)-iskipt(nbl1,1)) + std::abs(iskip(nbl2,1)-iskipt(nbl2,1));
                    nskpdif += std::abs(jskip(nbl1,1)-jskipt(nbl1,1)) + std::abs(jskip(nbl2,1)-jskipt(nbl2,1));
                    nskpdif += std::abs(kskip(nbl1,1)-kskipt(nbl1,1)) + std::abs(kskip(nbl2,1)-kskipt(nbl2,1));
                    iskip(nbl1,1)=iskipt(nbl1,1); iskip(nbl2,1)=iskipt(nbl2,1);
                    jskip(nbl1,1)=jskipt(nbl1,1); jskip(nbl2,1)=jskipt(nbl2,1);
                    kskip(nbl1,1)=kskipt(nbl1,1); kskip(nbl2,1)=kskipt(nbl2,1);
                }
                if (nskpdif == 0) break; // goto 850
            }
        }
        // label 850:

        // deforming mesh data
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        fortran_read_list(iunit5, &ndefrm);
        fprintf(fortran_get_unit(iunit11),
            " moving grid data - deforming surface (forced motion):\n");
        fprintf(fortran_get_unit(iunit11), "  ndefrm\n");
        if (ndefrm >= 0) fprintf(fortran_get_unit(iunit11), "%8d\n", ndefrm);


        if (ndefrm > 0) {
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            fortran_read_list(iunit5, &realval[1]);
            double lrefdef = realval[1];
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            fprintf(fortran_get_unit(iunit11), "    lref\n");
            fprintf(fortran_get_unit(iunit11), "%8.4f\n", (float)lrefdef);
            fprintf(fortran_get_unit(iunit11),
                "    grid   idefrm    rfreq u/omegax v/omegay w/omegaz    xorig    yorig    zorig\n");
            for (int ndef = 1; ndef <= ndefrm; ndef++) {
                int igr, idef;
                fortran_read_list(iunit5, &igr, &idef,
                                  &realval[1], &realval[2], &realval[3], &realval[4],
                                  &realval[5], &realval[6], &realval[7]);
                double rfr=realval[1], omgx=realval[2], omgy=realval[3], omgz=realval[4];
                double origx=realval[5], origy=realval[6], origz=realval[7];
                nbl = nblg(igr); ncg = ncgg(igr);
                idefrm(nbl) = 1;
                nsegdfrm(nbl) = nsegdfrm(nbl) + 1;
                int iseg = nsegdfrm(nbl);
                idfrmseg(nbl,iseg) = idef;
                if (idef == 2 || idef == 3) {
                    rfrqrae(nbl,iseg) = (float)(2.*pi*rfr/lrefdef);
                    omgxae(nbl,iseg)  = (float)(omgx/radtodeg);
                    omgyae(nbl,iseg)  = (float)(omgy/radtodeg);
                    omgzae(nbl,iseg)  = (float)(omgz/radtodeg);
                    xorgae(nbl,iseg)  = (float)origx; yorgae(nbl,iseg)  = (float)origy; zorgae(nbl,iseg)  = (float)origz;
                    xorgae0(nbl,iseg) = (float)origx; yorgae0(nbl,iseg) = (float)origy; zorgae0(nbl,iseg) = (float)origz;
                    fprintf(fortran_get_unit(iunit11), "%8d%9d%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
                            igr, idef, (float)rfr, (float)omgx, (float)omgy, (float)omgz,
                            (float)origx, (float)origy, (float)origz);
                } else if (idef == 1) {
                    rfrqtae(nbl,iseg) = (float)(2.*pi*rfr/lrefdef);
                    utrnsae(nbl,iseg) = (float)omgx; vtrnsae(nbl,iseg) = (float)omgy; wtrnsae(nbl,iseg) = (float)omgz;
                    xorgae(nbl,iseg)  = (float)origx; yorgae(nbl,iseg)  = (float)origy; zorgae(nbl,iseg)  = (float)origz;
                    xorgae0(nbl,iseg) = (float)origx; yorgae0(nbl,iseg) = (float)origy; zorgae0(nbl,iseg) = (float)origz;
                    fprintf(fortran_get_unit(iunit11), "%8d%9d%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
                            igr, idef, (float)rfr, (float)omgx, (float)omgy, (float)omgz,
                            (float)origx, (float)origy, (float)origz);
                } else {
                    fprintf(fortran_get_unit(iunit11), "%8d%9d%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
                            igr, idfrmseg(nbl,iseg),
                            rfrqrae(nbl,iseg), omgxae(nbl,iseg), omgyae(nbl,iseg), omgzae(nbl,iseg),
                            xorgae(nbl,iseg), yorgae(nbl,iseg), zorgae(nbl,iseg));
                }
                if (ncg > 0) {
                    for (n = 1; n <= ncg; n++) {
                        nbl++;
                        idefrm(nbl)=idefrm(nbl-1); nsegdfrm(nbl)=nsegdfrm(nbl-1);
                        idfrmseg(nbl,iseg)=idfrmseg(nbl-1,iseg);
                        rfrqrae(nbl,iseg)=rfrqrae(nbl-1,iseg); omgxae(nbl,iseg)=omgxae(nbl-1,iseg);
                        omgyae(nbl,iseg)=omgyae(nbl-1,iseg);   omgzae(nbl,iseg)=omgzae(nbl-1,iseg);
                        rfrqtae(nbl,iseg)=rfrqtae(nbl-1,iseg); utrnsae(nbl,iseg)=utrnsae(nbl-1,iseg);
                        vtrnsae(nbl,iseg)=vtrnsae(nbl-1,iseg); wtrnsae(nbl,iseg)=wtrnsae(nbl-1,iseg);
                        xorgae(nbl,iseg)=xorgae(nbl-1,iseg);   yorgae(nbl,iseg)=yorgae(nbl-1,iseg);
                        zorgae(nbl,iseg)=zorgae(nbl-1,iseg);   xorgae0(nbl,iseg)=xorgae0(nbl-1,iseg);
                        yorgae0(nbl,iseg)=yorgae0(nbl-1,iseg); zorgae0(nbl,iseg)=zorgae0(nbl-1,iseg);
                    }
                }
            }
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            fprintf(fortran_get_unit(iunit11),
                "    grid     icsi     icsf     jcsi     jcsf     kcsi     kcsf\n");
            // reset nsegdfrm
            for (nbl = 1; nbl <= nblock; nbl++) nsegdfrm(nbl) = 0;
            for (int ndef = 1; ndef <= ndefrm; ndef++) {
                int igr, is, ie, js, je, ks, ke;
                fortran_read_list(iunit5, &igr, &is, &ie, &js, &je, &ks, &ke);
                fprintf(fortran_get_unit(iunit11), "%8d%9d%9d%9d%9d%9d%9d\n",
                        igr, is, ie, js, je, ks, ke);
                nbl = nblg(igr); ncg = ncgg(igr);
                nsegdfrm(nbl) = nsegdfrm(nbl) + 1;
                int iseg = nsegdfrm(nbl);
                icsi(nbl,iseg)=is; icsf(nbl,iseg)=ie;
                jcsi(nbl,iseg)=js; jcsf(nbl,iseg)=je;
                kcsi(nbl,iseg)=ks; kcsf(nbl,iseg)=ke;
                if (ncg > 0) {
                    for (n = 1; n <= ncg; n++) {
                        nbl++;
                        nsegdfrm(nbl)=nsegdfrm(nbl-1);
                        icsi(nbl,iseg)=icsi(nbl-1,iseg)/2+1; icsf(nbl,iseg)=icsf(nbl-1,iseg)/2+1;
                        jcsi(nbl,iseg)=jcsi(nbl-1,iseg)/2+1; jcsf(nbl,iseg)=jcsf(nbl-1,iseg)/2+1;
                        kcsi(nbl,iseg)=kcsi(nbl-1,iseg)/2+1; kcsf(nbl,iseg)=kcsf(nbl-1,iseg)/2+1;
                    }
                }
            }
        } else if (ndefrm < 0) {
            int ndefrm0 = 0;
            {
                FortranArray1DRef<int> nblg_ref = nblg; FortranArray1DRef<int> ncgg_ref = ncgg;
                FortranArray1DRef<int> idimg_ref = idimg; FortranArray1DRef<int> jdimg_ref = jdimg;
                FortranArray1DRef<int> kdimg_ref = kdimg;
                FortranArray4DRef<int> ibcinfo_ref = ibcinfo; FortranArray4DRef<int> jbcinfo_ref = jbcinfo;
                FortranArray4DRef<int> kbcinfo_ref = kbcinfo;
                FortranArray1DRef<int> nbci0_ref = nbci0; FortranArray1DRef<int> nbcidim_ref = nbcidim;
                FortranArray1DRef<int> nbcj0_ref = nbcj0; FortranArray1DRef<int> nbcjdim_ref = nbcjdim;
                FortranArray1DRef<int> nbck0_ref = nbck0; FortranArray1DRef<int> nbckdim_ref = nbckdim;
                FortranArray1DRef<int> idefrm_ref = idefrm; FortranArray1DRef<int> nsegdfrm_ref = nsegdfrm;
                FortranArray2DRef<int> icsi_ref = icsi; FortranArray2DRef<int> icsf_ref = icsf;
                FortranArray2DRef<int> jcsi_ref = jcsi; FortranArray2DRef<int> jcsf_ref = jcsf;
                FortranArray2DRef<int> kcsi_ref = kcsi; FortranArray2DRef<int> kcsf_ref = kcsf;
                setseg_ns::setseg(maxgr, maxbl, maxseg, nblg_ref, ncgg_ref, idimg_ref, jdimg_ref,
                                  kdimg_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, nbci0_ref,
                                  nbcidim_ref, nbcj0_ref, nbcjdim_ref, nbck0_ref, nbckdim_ref,
                                  ndefrm0, idefrm_ref, nsegdfrm_ref, icsi_ref, icsf_ref,
                                  jcsi_ref, jcsf_ref, kcsi_ref, kcsf_ref, maxsegdg, ngrid);
            }
            fprintf(fortran_get_unit(iunit11), "%8d\n", ndefrm0);
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            fortran_read_list(iunit5, &realval[1]);
            double lrefdef = realval[1];
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            fprintf(fortran_get_unit(iunit11), "    lref\n");
            fprintf(fortran_get_unit(iunit11), "%8.4f\n", (float)lrefdef);
            fprintf(fortran_get_unit(iunit11),
                "    grid   idefrm    rfreq u/omegax v/omegay w/omegaz    xorig    yorig    zorig\n");
            int idef_last = 0;
            double rfr_last=0., omgx_last=0., omgy_last=0., omgz_last=0.;
            double origx_last=0., origy_last=0., origz_last=0.;
            for (int ndef = 1; ndef <= std::abs(ndefrm); ndef++) {
                int igr, idef;
                fortran_read_list(iunit5, &igr, &idef,
                                  &realval[1], &realval[2], &realval[3], &realval[4],
                                  &realval[5], &realval[6], &realval[7]);
                idef_last=idef; rfr_last=realval[1]; omgx_last=realval[2]; omgy_last=realval[3];
                omgz_last=realval[4]; origx_last=realval[5]; origy_last=realval[6]; origz_last=realval[7];
            }
            for (igrid = 1; igrid <= ngrid; igrid++) {
                nbl = nblg(igrid); ncg = ncgg(igrid);
                for (int iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {
                    idfrmseg(nbl,iseg) = idef_last;
                    if (idef_last == 2 || idef_last == 3) {
                        rfrqrae(nbl,iseg)=(float)(2.*pi*rfr_last/lrefdef);
                        omgxae(nbl,iseg)=(float)(omgx_last/radtodeg); omgyae(nbl,iseg)=(float)(omgy_last/radtodeg);
                        omgzae(nbl,iseg)=(float)(omgz_last/radtodeg);
                        xorgae(nbl,iseg)=(float)origx_last; yorgae(nbl,iseg)=(float)origy_last; zorgae(nbl,iseg)=(float)origz_last;
                        xorgae0(nbl,iseg)=(float)origx_last; yorgae0(nbl,iseg)=(float)origy_last; zorgae0(nbl,iseg)=(float)origz_last;
                        fprintf(fortran_get_unit(iunit11), "%8d%9d%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
                                igrid, idef_last, (float)rfr_last, (float)omgx_last, (float)omgy_last, (float)omgz_last,
                                (float)origx_last, (float)origy_last, (float)origz_last);
                    } else if (idef_last == 1) {
                        rfrqtae(nbl,iseg)=(float)(2.*pi*rfr_last/lrefdef);
                        utrnsae(nbl,iseg)=(float)omgx_last; vtrnsae(nbl,iseg)=(float)omgy_last; wtrnsae(nbl,iseg)=(float)omgz_last;
                        xorgae(nbl,iseg)=(float)origx_last; yorgae(nbl,iseg)=(float)origy_last; zorgae(nbl,iseg)=(float)origz_last;
                        xorgae0(nbl,iseg)=(float)origx_last; yorgae0(nbl,iseg)=(float)origy_last; zorgae0(nbl,iseg)=(float)origz_last;
                        fprintf(fortran_get_unit(iunit11), "%8d%9d%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
                                igrid, idef_last, (float)rfr_last, (float)omgx_last, (float)omgy_last, (float)omgz_last,
                                (float)origx_last, (float)origy_last, (float)origz_last);
                    } else {
                        fprintf(fortran_get_unit(iunit11), "%8d%9d%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
                                igrid, idfrmseg(nbl,iseg),
                                rfrqrae(nbl,iseg), omgxae(nbl,iseg), omgyae(nbl,iseg), omgzae(nbl,iseg),
                                xorgae(nbl,iseg), yorgae(nbl,iseg), zorgae(nbl,iseg));
                    }
                }
                if (ncg > 0) {
                    for (n = 1; n <= ncg; n++) {
                        nbl++;
                        idefrm(nbl)=idefrm(nbl-1); nsegdfrm(nbl)=nsegdfrm(nbl-1);
                        for (int iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {
                            idfrmseg(nbl,iseg)=idfrmseg(nbl-1,iseg);
                            rfrqrae(nbl,iseg)=rfrqrae(nbl-1,iseg); omgxae(nbl,iseg)=omgxae(nbl-1,iseg);
                            omgyae(nbl,iseg)=omgyae(nbl-1,iseg);   omgzae(nbl,iseg)=omgzae(nbl-1,iseg);
                            rfrqtae(nbl,iseg)=rfrqtae(nbl-1,iseg); utrnsae(nbl,iseg)=utrnsae(nbl-1,iseg);
                            vtrnsae(nbl,iseg)=vtrnsae(nbl-1,iseg); wtrnsae(nbl,iseg)=wtrnsae(nbl-1,iseg);
                            xorgae(nbl,iseg)=xorgae(nbl-1,iseg);   yorgae(nbl,iseg)=yorgae(nbl-1,iseg);
                            zorgae(nbl,iseg)=zorgae(nbl-1,iseg);   xorgae0(nbl,iseg)=xorgae0(nbl-1,iseg);
                            yorgae0(nbl,iseg)=yorgae0(nbl-1,iseg); zorgae0(nbl,iseg)=zorgae0(nbl-1,iseg);
                        }
                    }
                }
            }
        } // end if ndefrm > 0 / < 0
    } // end if iunst > 1 || idef_ss > 0


    // ndefrm < 0 continuation: read index ranges and output
    // (already handled above in the ndefrm < 0 branch)
    // ndefrm == 0 case: skip headers
    if (iunst > 1 || idef_ss > 0) {
        // naesrf - aeroelastic surface data
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        fortran_read_list(iunit5, &naesrf);
        fprintf(fortran_get_unit(iunit11),
            " moving grid data - aeroelastic surface (aeroelastic motion):\n");
        fprintf(fortran_get_unit(iunit11), "  naesrf\n");
        fprintf(fortran_get_unit(iunit11), "%8d\n", naesrf);
        if (naesrf > 0) {
            if (naesrf > maxaes) {
                fprintf(fortran_get_unit(iunit11), " input error: naesrf > maxaes\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            for (int naes = 1; naes <= naesrf; naes++) {
                { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
                fprintf(fortran_get_unit(iunit11),
                    "  iaesrf    ngrid    grefl     uinf     qinf   nmodes iskyhook\n");
                int iaes, ngd, nmodes, iskyhk;
                fortran_read_list(iunit5, &iaes, &ngd,
                                  &realval[1], &realval[2], &realval[3], &nmodes, &iskyhk);
                double grefl=realval[1], uinf_ae=realval[2], qinf_ae=realval[3];
                aesrfdat(1,iaes)=(float)iskyhk; aesrfdat(2,iaes)=(float)grefl;
                aesrfdat(3,iaes)=(float)uinf_ae; aesrfdat(4,iaes)=(float)qinf_ae;
                aesrfdat(5,iaes)=(float)nmodes;
                fprintf(fortran_get_unit(iunit11), "%8d%9d%9.4f%9.4f%9.4f%9d%9d\n",
                        iaes, ngd, (float)grefl, (float)uinf_ae, (float)qinf_ae, nmodes, iskyhk);
                if ((float)grefl <= 0.0f) {
                    fprintf(fortran_get_unit(iunit11), " input error: invalid grefl\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (nmodes > nmds) {
                    fprintf(fortran_get_unit(iunit11), " input error: number of modes exceeds nmds\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
                fprintf(fortran_get_unit(iunit11),
                    "    freq    gmass     damp x0(2n-1)   x0(2n)  gf0(2n)\n");
                for (int nm = 1; nm <= nmodes; nm++) {
                    int nm1 = 2*nm-1, nm2 = 2*nm;
                    fortran_read_list(iunit5, &realval[1], &realval[2], &realval[3],
                                      &realval[4], &realval[5], &realval[6]);
                    double freq1=realval[1], gmass1=realval[2], damp1=realval[3];
                    double x01=realval[4], x02=realval[5], gf02=realval[6];
                    freq(nm,iaes)=(float)freq1; gmass(nm,iaes)=(float)gmass1; damp(nm,iaes)=(float)damp1;
                    x0(nm1,iaes)=(float)x01; x0(nm2,iaes)=(float)x02;
                    gf0(nm1,iaes)=(float)gf02; gf0(nm2,iaes)=(float)gf02;
                    fprintf(fortran_get_unit(iunit11), "%8.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
                            (float)freq1, (float)gmass1, (float)damp1, (float)x01, (float)x02, (float)gf02);
                    if ((float)freq1 <= 0.0f) {
                        fprintf(fortran_get_unit(iunit11), " input error: invalid frequency\n");
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                    if ((float)damp1 >= 1.0f) {
                        fprintf(fortran_get_unit(iunit11), " input error: invalid  damping\n");
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                }
                { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
                fprintf(fortran_get_unit(iunit11), "  moddfl      amp     freq       t0\n");
                for (int nm = 1; nm <= nmodes; nm++) {
                    int moddfl;
                    fortran_read_list(iunit5, &moddfl, &realval[1], &realval[2], &realval[3]);
                    double amp=realval[1], freqp=realval[2], t0=realval[3];
                    perturb(nm,iaes,1)=(float)moddfl; perturb(nm,iaes,2)=(float)amp;
                    perturb(nm,iaes,3)=(float)freqp;  perturb(nm,iaes,4)=(float)t0;
                    fprintf(fortran_get_unit(iunit11), "%8d%9.4f%9.4f%9.4f\n",
                            moddfl, (float)amp, (float)freqp, (float)t0);
                }
                { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
                fprintf(fortran_get_unit(iunit11),
                    "    grid     iaei     iaef     jaei     jaef     kaei     kaef\n");
                if (ngd > 0) {
                    for (int ng = 1; ng <= ngd; ng++) {
                        int igr, is, ie, js, je, ks, ke;
                        fortran_read_list(iunit5, &igr, &is, &ie, &js, &je, &ks, &ke);
                        fprintf(fortran_get_unit(iunit11), "%8d%9d%9d%9d%9d%9d%9d\n",
                                igr, is, ie, js, je, ks, ke);
                        nbl = nblg(igr); ncg = ncgg(igr);
                        nsegdfrm(nbl) = nsegdfrm(nbl) + 1;
                        int iseg = nsegdfrm(nbl);
                        icsi(nbl,iseg)=is; icsf(nbl,iseg)=ie;
                        jcsi(nbl,iseg)=js; jcsf(nbl,iseg)=je;
                        kcsi(nbl,iseg)=ks; kcsf(nbl,iseg)=ke;
                        idfrmseg(nbl,iseg)=99; iaesurf(nbl,iseg)=iaes; idefrm(nbl)=1;
                        if (ncg > 0) {
                            for (n = 1; n <= ncg; n++) {
                                nbl++;
                                idefrm(nbl)=idefrm(nbl-1); nsegdfrm(nbl)=nsegdfrm(nbl-1);
                                idfrmseg(nbl,iseg)=idfrmseg(nbl-1,iseg);
                                iaesurf(nbl,iseg)=iaesurf(nbl-1,iseg);
                                icsi(nbl,iseg)=icsi(nbl-1,iseg)/2+1; icsf(nbl,iseg)=icsf(nbl-1,iseg)/2+1;
                                jcsi(nbl,iseg)=jcsi(nbl-1,iseg)/2+1; jcsf(nbl,iseg)=jcsf(nbl-1,iseg)/2+1;
                                kcsi(nbl,iseg)=kcsi(nbl-1,iseg)/2+1; kcsf(nbl,iseg)=kcsf(nbl-1,iseg)/2+1;
                            }
                        }
                    }
                } else if (ngd < 0) {
                    int ndefrm0 = 0;
                    {
                        FortranArray1DRef<int> nblg_r=nblg, ncgg_r=ncgg, idimg_r=idimg, jdimg_r=jdimg, kdimg_r=kdimg;
                        FortranArray4DRef<int> ibcinfo_r=ibcinfo, jbcinfo_r=jbcinfo, kbcinfo_r=kbcinfo;
                        FortranArray1DRef<int> nbci0_r=nbci0, nbcidim_r=nbcidim, nbcj0_r=nbcj0, nbcjdim_r=nbcjdim;
                        FortranArray1DRef<int> nbck0_r=nbck0, nbckdim_r=nbckdim;
                        FortranArray1DRef<int> idefrm_r=idefrm, nsegdfrm_r=nsegdfrm;
                        FortranArray2DRef<int> icsi_r=icsi, icsf_r=icsf, jcsi_r=jcsi, jcsf_r=jcsf, kcsi_r=kcsi, kcsf_r=kcsf;
                        setseg_ns::setseg(maxgr, maxbl, maxseg, nblg_r, ncgg_r, idimg_r, jdimg_r,
                                          kdimg_r, ibcinfo_r, jbcinfo_r, kbcinfo_r, nbci0_r,
                                          nbcidim_r, nbcj0_r, nbcjdim_r, nbck0_r, nbckdim_r,
                                          ndefrm0, idefrm_r, nsegdfrm_r, icsi_r, icsf_r,
                                          jcsi_r, jcsf_r, kcsi_r, kcsf_r, maxsegdg, ngrid);
                    }
                    for (int ng = 1; ng <= std::abs(ngd); ng++) {
                        int igr, is, ie, js, je, ks, ke;
                        fortran_read_list(iunit5, &igr, &is, &ie, &js, &je, &ks, &ke);
                    }
                    for (igrid = 1; igrid <= ngrid; igrid++) {
                        nbl = nblg(igrid); ncg = ncgg(igrid);
                        for (int iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {
                            idfrmseg(nbl,iseg)=99; iaesurf(nbl,iseg)=iaes;
                            fprintf(fortran_get_unit(iunit11), "%8d%9d%9d%9d%9d%9d%9d\n",
                                    igrid, icsi(nbl,iseg), icsf(nbl,iseg),
                                    jcsi(nbl,iseg), jcsf(nbl,iseg), kcsi(nbl,iseg), kcsf(nbl,iseg));
                        }
                        if (ncg > 0) {
                            for (n = 1; n <= ncg; n++) {
                                nbl++;
                                for (int iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {
                                    idfrmseg(nbl,iseg)=idfrmseg(nbl-1,iseg);
                                    iaesurf(nbl,iseg)=iaesurf(nbl-1,iseg);
                                }
                            }
                        }
                    }
                }
            }
            mxmds = std::max(mxmds, nmds);
        } else {
            // naesrf == 0: skip 4 header lines
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            fprintf(fortran_get_unit(iunit11),
                "  iaesrf    ngrid    grefl     uinf     qinf   nmodes iskyhook\n");
            fprintf(fortran_get_unit(iunit11),
                "    freq    gmass     damp x0(2n-1)   x0(2n)  gf0(2n)\n");
            fprintf(fortran_get_unit(iunit11), "  moddfl      amp     freq       t0\n");
            fprintf(fortran_get_unit(iunit11),
                "    grid     iaei     iaef     jaei     jaef     kaei     kaef\n");
        }
        mxaes = std::max(mxaes, naesrf);


        // deforming mesh skip data
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        int nskip;
        fortran_read_list(iunit5, &nskip, &isktyp,
                          &realval[1], &realval[2], &realval[3], &realval[4], &nsprgit);
        beta1  = (float)realval[1]; alpha1 = (float)realval[2];
        beta2  = (float)realval[3]; alpha2 = (float)realval[4];
        fprintf(fortran_get_unit(iunit11),
            " moving grid data - data for field/multiblock mesh movement\n");
        fprintf(fortran_get_unit(iunit11),
            "   nskip   isktyp    beta1   alpha1    beta2   alpha2 nsprngit\n");
        int itcpadd = 0;
        if (std::abs(isktyp) < 1 || std::abs(isktyp) > 2) {
            fprintf(fortran_get_unit(iunit11), "%8d%9d%9.6f%9.6f%9.6f%9.6f%9d\n",
                    nskip, isktyp, beta1, alpha1, beta2, alpha2, nsprgit);
            fprintf(fortran_get_unit(iunit11), " input error: Invalid isktyp value\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        if (std::abs(isktyp) == 1) {
            if (nskip > 0) {
                if (nskip > ngrid) {
                    fprintf(fortran_get_unit(iunit11), "%8d%9d%9.6f%9.6f%9.6f%9.6f%9d\n",
                            nskip, isktyp, beta1, alpha1, beta2, alpha2, nsprgit);
                    fprintf(fortran_get_unit(iunit11), " input error: nskip > ngrid\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                fprintf(fortran_get_unit(iunit11), "%8d%9d%9.6f%9.6f%9.6f%9.6f%9d\n",
                        nskip, isktyp, beta1, alpha1, beta2, alpha2, nsprgit);
                if (beta1 < 0.2f) {
                    fprintf(fortran_get_unit(iunit11),
                        "\n WARNING: beta this small may severely compromise fidelity of surface geometry\n\n");
                }
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                for (int ng = 1; ng <= nskip; ng++) {
                    int igr, iskp, jskp, kskp;
                    fortran_read_list(iunit5, &igr, &iskp, &jskp, &kskp);
                    nbl = nblg(igr); ncg = ncgg(igr);
                    if (iskp == 0) iskp = idimg(nbl)-1;
                    if (jskp == 0) jskp = jdimg(nbl)-1;
                    if (kskp == 0) kskp = kdimg(nbl)-1;
                    if (i2d != 0) iskp = idimg(nbl)-1;
                    int istop = 0;
                    if ((idimg(nbl)-1)/iskp*iskp != (idimg(nbl)-1)) {
                        istop=1;
                        fprintf(fortran_get_unit(iunit11),
                            " stopping...iskip must divide evenly into idim-1 for grid%4d\n", igr);
                    }
                    if ((jdimg(nbl)-1)/jskp*jskp != (jdimg(nbl)-1)) {
                        istop=1;
                        fprintf(fortran_get_unit(iunit11),
                            " stopping...jskip must divide evenly into jdim-1 for grid%4d\n", igr);
                    }
                    if ((kdimg(nbl)-1)/kskp*kskp != (kdimg(nbl)-1)) {
                        istop=1;
                        fprintf(fortran_get_unit(iunit11),
                            " stopping...kskip must divide evenly into kdim-1 for grid%4d\n", igr);
                    }
                    if (istop == 1) termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    iskip(nbl,1) = iskp; jskip(nbl,1) = jskp; kskip(nbl,1) = kskp;
                }
            } else if (nskip < 0) {
                fprintf(fortran_get_unit(iunit11), "%8d%9d%9.6f%9.6f%9.6f%9.6f%9d\n",
                        ngrid, isktyp, beta1, alpha1, beta2, alpha2, nsprgit);
                if (beta1 < 0.2f) {
                    fprintf(fortran_get_unit(iunit11),
                        "\n WARNING: beta this small may severely compromise fidelity of surface geometry\n\n");
                }
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
                int iskp0=0, jskp0=0, kskp0=0;
                for (int ng = 1; ng <= std::abs(nskip); ng++) {
                    int igr;
                    fortran_read_list(iunit5, &igr, &iskp0, &jskp0, &kskp0);
                }
                for (int igr = 1; igr <= ngrid; igr++) {
                    nbl = nblg(igr); ncg = ncgg(igr);
                    int iskp = iskp0, jskp = jskp0, kskp = kskp0;
                    if (iskp0 == 0) iskp = idimg(nbl)-1;
                    if (jskp0 == 0) jskp = jdimg(nbl)-1;
                    if (kskp0 == 0) kskp = kdimg(nbl)-1;
                    int istop = 0;
                    if ((idimg(nbl)-1)/iskp*iskp != (idimg(nbl)-1)) {
                        istop=1;
                        fprintf(fortran_get_unit(iunit11),
                            " stopping...iskip must divide evenly into idim-1 for grid%4d\n", igr);
                    }
                    if ((jdimg(nbl)-1)/jskp*jskp != (jdimg(nbl)-1)) {
                        istop=1;
                        fprintf(fortran_get_unit(iunit11),
                            " stopping...jskip must divide evenly into jdim-1 for grid%4d\n", igr);
                    }
                    if ((kdimg(nbl)-1)/kskp*kskp != (kdimg(nbl)-1)) {
                        istop=1;
                        fprintf(fortran_get_unit(iunit11),
                            " stopping...kskip must divide evenly into kdim-1 for grid%4d\n", igr);
                    }
                    if (istop == 1) termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    iskip(nbl,1) = iskp; jskip(nbl,1) = jskp; kskip(nbl,1) = kskp;
                }
            } else {
                // nskip == 0
                fprintf(fortran_get_unit(iunit11), "%8d%9d%9.6f%9.6f%9.6f%9.6f%9d\n",
                        ngrid, isktyp, beta1, alpha1, beta2, alpha2, nsprgit);
                if (beta1 < 0.2f) {
                    fprintf(fortran_get_unit(iunit11),
                        "\n WARNING: beta this small may severely compromise fidelity of surface geometry\n\n");
                }
                { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }
            }
        } else {
            // isktyp == 2
            if (nskip != ngrid && nskip != 0) {
                fprintf(fortran_get_unit(iunit11), "%8d%9d%9.6f%9.6f%9.6f%9.6f%9d\n",
                        nskip, isktyp, beta1, alpha1, beta2, alpha2, nsprgit);
                fprintf(fortran_get_unit(iunit11),
                    " input error: nskip .ne. ngrid. May mean improper isktyp value.\n\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            fprintf(fortran_get_unit(iunit11), "%8d%9d%9.6f%9.6f%9.6f%9.6f%9d\n",
                    ngrid, isktyp, beta1, alpha1, beta2, alpha2, nsprgit);
            if (beta1 < 0.2f) {
                fprintf(fortran_get_unit(iunit11),
                    "\n WARNING: beta this small may severely compromise fidelity of surface geometry\n\n");
            }
            { char _skip2[256]; fgets(_skip2, sizeof(_skip2), fortran_get_unit(iunit5)); }


            if (nskip == ngrid) {
                for (int ng = 1; ng <= nskip; ng++) {
                    { char _skip3[256]; fgets(_skip3, sizeof(_skip3), fortran_get_unit(iunit5)); }
                    int igr;
                    fortran_read_list(iunit5, &igr, &nskpi1(ng), &nskpj1(ng), &nskpk1(ng));
                    if (igr != ng) {
                        fprintf(fortran_get_unit(iunit11),
                            "\n ERROR: control point input must include all blocks, and be entered in \nascending block order\n\n");
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                    if (nskpi1(ng) > 500 || nskpj1(ng) > 500 || nskpk1(ng) > 500) {
                        fprintf(fortran_get_unit(iunit11),
                            "\n ERROR: control point input number of control points must be <= 500\n\n");
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                    nbl = nblg(igr);
                    { char _skip3[256]; fgets(_skip3, sizeof(_skip3), fortran_get_unit(iunit5)); }
                    int iis = -9, iie = 0;
                    for (int jj = 1; jj <= 50; jj++) {
                        iis += 10; iie += 10;
                        if (iie > nskpi1(ng)) iie = nskpi1(ng);
                        {
                            int _cnt = iie - iis + 1;
                            int _buf[10];
                            fortran_read_list_array(iunit5, _buf, _cnt);
                            for (int _m = 0; _m < _cnt; _m++) iskip(nbl, iis + _m) = _buf[_m];
                        }
                        if (iie == nskpi1(ng)) break;
                    }
                    { char _skip3[256]; fgets(_skip3, sizeof(_skip3), fortran_get_unit(iunit5)); }
                    iis = -9; iie = 0;
                    for (int jj = 1; jj <= 50; jj++) {
                        iis += 10; iie += 10;
                        if (iie > nskpj1(ng)) iie = nskpj1(ng);
                        {
                            int _cnt = iie - iis + 1;
                            int _buf[10];
                            fortran_read_list_array(iunit5, _buf, _cnt);
                            for (int _m = 0; _m < _cnt; _m++) jskip(nbl, iis + _m) = _buf[_m];
                        }
                        if (iie == nskpj1(ng)) break;
                    }
                    { char _skip3[256]; fgets(_skip3, sizeof(_skip3), fortran_get_unit(iunit5)); }
                    iis = -9; iie = 0;
                    for (int jj = 1; jj <= 50; jj++) {
                        iis += 10; iie += 10;
                        if (iie > nskpk1(ng)) iie = nskpk1(ng);
                        {
                            int _cnt = iie - iis + 1;
                            int _buf[10];
                            fortran_read_list_array(iunit5, _buf, _cnt);
                            for (int _m = 0; _m < _cnt; _m++) kskip(nbl, iis + _m) = _buf[_m];
                        }
                        if (iie == nskpk1(ng)) break;
                    }
                    int idim = idimg(nbl), jdim = jdimg(nbl), kdim = kdimg(nbl);
                    if (iskip(nbl,nskpi1(ng)) != idim) {
                        fprintf(fortran_get_unit(iunit11), "\n STOPPING: ic(iskip) ne idim\n\n"); std::exit(1);
                    }
                    if (jskip(nbl,nskpj1(ng)) != jdim) {
                        fprintf(fortran_get_unit(iunit11), "\n STOPPING: ic(jskip) ne jdim\n\n"); std::exit(1);
                    }
                    if (kskip(nbl,nskpk1(ng)) != kdim) {
                        fprintf(fortran_get_unit(iunit11), "\n STOPPING: ic(kskip) ne kdim\n\n"); std::exit(1);
                    }
                    for (int jj = 1; jj <= nskpi1(ng); jj++) iskipt(nbl,jj) = iskip(nbl,jj);
                    for (int jj = 1; jj <= nskpj1(ng); jj++) jskipt(nbl,jj) = jskip(nbl,jj);
                    for (int jj = 1; jj <= nskpk1(ng); jj++) kskipt(nbl,jj) = kskip(nbl,jj);
                }
            }
            if (nskip == 0) {
                for (int ng = 1; ng <= ngrid; ng++) {
                    nbl = nblg(ng);
                    int idim = idimg(nbl), jdim = jdimg(nbl), kdim = kdimg(nbl);
                    nskpi1(ng)=2; nskpj1(ng)=2; nskpk1(ng)=2;
                    iskip(nbl,1)=1; iskip(nbl,2)=idim;
                    jskip(nbl,1)=1; jskip(nbl,2)=jdim;
                    kskip(nbl,1)=1; kskip(nbl,2)=kdim;
                }
                for (n = 1; n <= std::abs(nbli); n++) {
                    for (int ii = 1; ii <= 2; ii++) {
                        nbl = nblg(nblk(ii,n));
                        int idim = idimg(nbl), jdim = jdimg(nbl), kdim = kdimg(nbl);
                        for (int jj = 1; jj <= 2; jj++) {
                            if (isva(ii,jj,n) == 1) {
                                if (limblk(ii,1,n)>1 && limblk(ii,1,n)<idim)
                                    skordr(nbl, nskpi1(nblk(ii,n)), limblk(ii,1,n), iskip, maxbl);
                                if (limblk(ii,4,n)>1 && limblk(ii,4,n)<idim)
                                    skordr(nbl, nskpi1(nblk(ii,n)), limblk(ii,4,n), iskip, maxbl);
                            } else if (isva(ii,jj,n) == 2) {
                                if (limblk(ii,2,n)>1 && limblk(ii,2,n)<jdim)
                                    skordr(nbl, nskpj1(nblk(ii,n)), limblk(ii,2,n), jskip, maxbl);
                                if (limblk(ii,5,n)>1 && limblk(ii,5,n)<jdim)
                                    skordr(nbl, nskpj1(nblk(ii,n)), limblk(ii,5,n), jskip, maxbl);
                            } else if (isva(ii,jj,n) == 3) {
                                if (limblk(ii,3,n)>1 && limblk(ii,3,n)<kdim)
                                    skordr(nbl, nskpk1(nblk(ii,n)), limblk(ii,3,n), kskip, maxbl);
                                if (limblk(ii,6,n)>1 && limblk(ii,6,n)<kdim)
                                    skordr(nbl, nskpk1(nblk(ii,n)), limblk(ii,6,n), kskip, maxbl);
                            }
                        }
                    }
                }
            }
        } // end isktyp == 2
        // NOTE: the iunst>1||idef_ss>0 block continues below; it closes at line ~5252.


        // print skip data
        if (std::abs(isktyp) == 1) {
            fprintf(fortran_get_unit(iunit11), "    grid    iskip    jskip    kskip\n");
            for (int ng = 1; ng <= ngrid; ng++) {
                nbl = nblg(ng); ncg = ncgg(ng);
                fprintf(fortran_get_unit(iunit11), "%8d%9d%9d%9d\n",
                        ng, iskip(nbl,1), jskip(nbl,1), kskip(nbl,1));
                if (ncg > 0) {
                    for (n = 1; n <= ncg; n++) {
                        nbl++;
                        iskip(nbl,1) = iskip(nbl-1,1);
                        jskip(nbl,1) = jskip(nbl-1,1);
                        kskip(nbl,1) = kskip(nbl-1,1);
                    }
                }
            }
        }
        if (std::abs(isktyp) == 2) {
            fprintf(fortran_get_unit(iunit11), "    grid  nskpi1  nskpj1  nskpk1\n");
            for (int ng = 1; ng <= ngrid; ng++) {
                nbl = nblg(ng);
                fprintf(fortran_get_unit(iunit11), "%8d%8d%8d%8d\n",
                        ng, nskpi1(ng), nskpj1(ng), nskpk1(ng));
            }
        }
        // set mxdefseg
        for (nbl = 1; nbl <= nblock; nbl++) {
            mxdefseg = std::max(mxdefseg, nsegdfrm(nbl));
        }


        // For isktyp==1: set coarse grid skip values and build iskipt arrays
        if (std::abs(isktyp) == 1) {
            for (int ng = 1; ng <= ngrid; ng++) {
                nbl = nblg(ng); ncg = ncgg(ng);
                int idim = idimg(nbl), jdim = jdimg(nbl), kdim = kdimg(nbl);
                if (ncg > 0) {
                    for (n = 1; n <= ncg; n++) {
                        nbl++;
                        if (idim > 2) iskip(nbl,1) = iskip(nbl-1,1)/2;
                        else          iskip(nbl,1) = iskip(nbl-1,1);
                        jskip(nbl,1) = jskip(nbl-1,1)/2;
                        kskip(nbl,1) = kskip(nbl-1,1)/2;
                    }
                    nbl = nblg(ng);
                }
                int ncgt = ncg;
                if (idim == 2) nskpi1(ng) = 2;
                else           nskpi1(ng) = (idimg(nbl)-1)/iskip(nbl,1) + 1;
                nskpj1(ng) = (jdimg(nbl)-1)/jskip(nbl,1) + 1;
                nskpk1(ng) = (kdimg(nbl)-1)/kskip(nbl,1) + 1;
                iskipt(nbl,1) = 1; jskipt(nbl,1) = 1; kskipt(nbl,1) = 1;
                if (nskpi1(ng) > 500) {
                    fprintf(fortran_get_unit(iunit11),
                        "\n ERROR: nskpi1 > 500  Too many i-control points (iskip too small)\n\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (nskpj1(ng) > 500) {
                    fprintf(fortran_get_unit(iunit11),
                        "\n ERROR: nskpj1 > 500  Too many j-control points (jskip too small)\n\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                if (nskpk1(ng) > 500) {
                    fprintf(fortran_get_unit(iunit11),
                        "\n ERROR: nskpk1 > 500  Too many k-control points (kskip too small)\n\n");
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
                for (int jj = 2; jj <= nskpi1(ng); jj++) iskipt(nbl,jj) = iskipt(nbl,jj-1)+iskip(nbl,1);
                for (int jj = 2; jj <= nskpj1(ng); jj++) jskipt(nbl,jj) = jskipt(nbl,jj-1)+jskip(nbl,1);
                for (int jj = 2; jj <= nskpk1(ng); jj++) kskipt(nbl,jj) = kskipt(nbl,jj-1)+kskip(nbl,1);
            }
        } else {
            // isktyp==2: propagate to coarse grids
            for (int ng = 1; ng <= ngrid; ng++) {
                nbl = nblg(ng); ncg = ncgg(ng);
                int idim = idimg(nbl);
                if (ncg > 0) {
                    nbl = nblg(ng);
                    for (int ii = 1; ii <= nskpi1(ng); ii++) iskipt(nbl,ii) = iskip(nbl,ii);
                    if (idim > 2) {
                        for (n = 1; n <= ncg; n++) {
                            nbl++;
                            for (int ii = 1; ii <= nskpi1(ng); ii++) {
                                iskip(nbl,ii) = (iskip(nbl-1,ii)-1)/2 + 1;
                                iskipt(nbl,ii) = iskip(nbl,ii);
                                riskp(nbl,ii) = 0.5f*(float)(iskip(nbl-1,ii)-1) + 1.0f;
                            }
                            for (int ii = 2; ii <= nskpi1(ng)-1; ii++) {
                                if (iskip(nbl,ii) >= iskip(nbl,ii+1) && (float)iskip(nbl,ii) != riskp(nbl,ii))
                                    iskip(nbl,ii)--;
                                if (iskip(nbl,ii) <= iskip(nbl,ii-1) && (float)iskip(nbl,ii) != riskp(nbl,ii))
                                    iskip(nbl,ii)++;
                            }
                            for (int ii = 1; ii <= nskpi1(ng); ii++) iskipt(nbl,ii) = iskip(nbl,ii);
                        }
                    } else {
                        for (n = 1; n <= ncg; n++) {
                            nbl++;
                            for (int ii = 1; ii <= nskpi1(ng); ii++) {
                                iskip(nbl,ii) = iskip(nbl-1,ii);
                                iskipt(nbl,ii) = iskip(nbl,ii);
                            }
                        }
                    }
                    nbl = nblg(ng);
                    for (int ii = 1; ii <= nskpj1(ng); ii++) jskipt(nbl,ii) = jskip(nbl,ii);
                    for (n = 1; n <= ncg; n++) {
                        nbl++;
                        for (int ii = 1; ii <= nskpj1(ng); ii++) {
                            jskip(nbl,ii) = (jskip(nbl-1,ii)-1)/2 + 1;
                            rjskp(nbl,ii) = 0.5f*(float)(jskip(nbl-1,ii)-1) + 1.0f;
                        }
                        for (int ii = 2; ii <= nskpj1(ng)-1; ii++) {
                            if (jskip(nbl,ii) >= jskip(nbl,ii+1) && (float)jskip(nbl,ii) != rjskp(nbl,ii))
                                jskip(nbl,ii)--;
                            if (jskip(nbl,ii) <= jskip(nbl,ii-1) && (float)jskip(nbl,ii) != rjskp(nbl,ii))
                                jskip(nbl,ii)++;
                        }
                        for (int ii = 1; ii <= nskpj1(ng); ii++) jskipt(nbl,ii) = jskip(nbl,ii);
                    }
                    nbl = nblg(ng);
                    for (int ii = 1; ii <= nskpk1(ng); ii++) kskipt(nbl,ii) = kskip(nbl,ii);
                    for (n = 1; n <= ncg; n++) {
                        nbl++;
                        for (int ii = 1; ii <= nskpk1(ng); ii++) {
                            kskip(nbl,ii) = (kskip(nbl-1,ii)-1)/2 + 1;
                            rkskp(nbl,ii) = 0.5f*(float)(kskip(nbl-1,ii)-1) + 1.0f;
                        }
                        for (int ii = 2; ii <= nskpk1(ng)-1; ii++) {
                            if (kskip(nbl,ii) >= kskip(nbl,ii+1) && (float)kskip(nbl,ii) != rkskp(nbl,ii))
                                kskip(nbl,ii)--;
                            if (kskip(nbl,ii) <= kskip(nbl,ii-1) && (float)kskip(nbl,ii) != rkskp(nbl,ii))
                                kskip(nbl,ii)++;
                        }
                        for (int ii = 1; ii <= nskpk1(ng); ii++) kskipt(nbl,ii) = kskip(nbl,ii);
                    }
                }
            }
        }
        // write meshdef.inp if needed
        if ((std::abs(isktyp)==1 || (std::abs(isktyp)==2 && nskip==0) || (itcpadd>1)) && meshdef==1) {
            fortran_open_unit(196, "meshdef.inp", "w");
            for (int ng = 1; ng <= ngrid; ng++) {
                nbl = nblg(ng);
                fprintf(fortran_get_unit(196), "   GRID   NIND   NJND   NKND \n");
                fprintf(fortran_get_unit(196), "%7d%7d%7d%7d\n", ng, nskpi1(ng), nskpj1(ng), nskpk1(ng));
                fprintf(fortran_get_unit(196), "************************** I NODE INDICES ****************************\n");
                int iis = -9, iie = 0;
                for (int jj = 1; jj <= 50; jj++) {
                    iis += 10; iie += 10;
                    if (iie > nskpi1(ng)) iie = nskpi1(ng);
                    for (int ii = iis; ii <= iie; ii++) fprintf(fortran_get_unit(196), "%7d", iskipt(nbl,ii));
                    fprintf(fortran_get_unit(196), "\n");
                    if (iie == nskpi1(ng)) break;
                }
                fprintf(fortran_get_unit(196), "************************** J NODE INDICES ****************************\n");
                iis = -9; iie = 0;
                for (int jj = 1; jj <= 50; jj++) {
                    iis += 10; iie += 10;
                    if (iie > nskpj1(ng)) iie = nskpj1(ng);
                    for (int ii = iis; ii <= iie; ii++) fprintf(fortran_get_unit(196), "%7d", jskipt(nbl,ii));
                    fprintf(fortran_get_unit(196), "\n");
                    if (iie == nskpj1(ng)) break;
                }
                fprintf(fortran_get_unit(196), "************************** K NODE INDICES ****************************\n");
                iis = -9; iie = 0;
                for (int jj = 1; jj <= 50; jj++) {
                    iis += 10; iie += 10;
                    if (iie > nskpk1(ng)) iie = nskpk1(ng);
                    for (int ii = iis; ii <= iie; ii++) fprintf(fortran_get_unit(196), "%7d", kskipt(nbl,ii));
                    fprintf(fortran_get_unit(196), "\n");
                    if (iie == nskpk1(ng)) break;
                }
            }
            fortran_close_unit(196);
        }


        // print isktyp==2 control point indices
        if (std::abs(isktyp) == 2) {
            for (int ng = 1; ng <= ngrid; ng++) {
                nbl = nblg(ng); ncg = ncgg(ng);
                int ncgt = ncg;
                fprintf(fortran_get_unit(iunit11), "      ng     nipt     njpt     nkpt  \n");
                fprintf(fortran_get_unit(iunit11), "%8d%9d%9d%9d\n", ng, nskpi1(ng), nskpj1(ng), nskpk1(ng));
                nbl = nblg(ng);
                fprintf(fortran_get_unit(iunit11), "    control point i-indices for grid levels  ");
                for (int i = 1; i <= ncgt+1; i++) fprintf(fortran_get_unit(iunit11), "%4d", nbl+i-1);
                fprintf(fortran_get_unit(iunit11), "\n");
                for (int i = 1; i <= nskpi1(ng); i++) {
                    for (int ii = 1; ii <= ncgt+1; ii++) fprintf(fortran_get_unit(iunit11), "%9d", iskipt(nbl+ii-1,i));
                    fprintf(fortran_get_unit(iunit11), "\n");
                }
                nbl = nblg(ng);
                fprintf(fortran_get_unit(iunit11), "    control point j-indices for grid levels  ");
                for (int i = 1; i <= ncgt+1; i++) fprintf(fortran_get_unit(iunit11), "%4d", nbl+i-1);
                fprintf(fortran_get_unit(iunit11), "\n");
                for (int i = 1; i <= nskpj1(ng); i++) {
                    for (int ii = 1; ii <= ncgt+1; ii++) fprintf(fortran_get_unit(iunit11), "%9d", jskipt(nbl+ii-1,i));
                    fprintf(fortran_get_unit(iunit11), "\n");
                }
                nbl = nblg(ng);
                fprintf(fortran_get_unit(iunit11), "    control point k-indices for grid levels  ");
                for (int i = 1; i <= ncgt+1; i++) fprintf(fortran_get_unit(iunit11), "%4d", nbl+i-1);
                fprintf(fortran_get_unit(iunit11), "\n");
                for (int i = 1; i <= nskpk1(ng); i++) {
                    for (int ii = 1; ii <= ncgt+1; ii++) fprintf(fortran_get_unit(iunit11), "%9d", kskipt(nbl+ii-1,i));
                    fprintf(fortran_get_unit(iunit11), "\n");
                }
            }
        }
        // set idefrm=999 for blocks not yet identified as deforming
        for (nbl = 1; nbl <= nblock; nbl++) {
            if (idefrm(nbl) == 0) {
                idefrm(nbl) = 999;
                nsegdfrm(nbl) = 1;
                for (int iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {
                    icsi(nbl,iseg)=1; icsf(nbl,iseg)=idimg(nbl);
                    jcsi(nbl,iseg)=1; jcsf(nbl,iseg)=jdimg(nbl);
                    kcsi(nbl,iseg)=1; kcsf(nbl,iseg)=kdimg(nbl);
                }
            }
        }
        // multi-motion coupling data
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
        int ncoupl;
        fortran_read_list(iunit5, &ncoupl);
        fprintf(fortran_get_unit(iunit11), " moving grid data - multi-motion coupling\n");
        fprintf(fortran_get_unit(iunit11), "  ncoupl\n");
        if (ncoupl > 0) {
            fprintf(fortran_get_unit(iunit11), "%8d\n", ncoupl);
            if (ncoupl > ngrid) {
                fprintf(fortran_get_unit(iunit11), " input error: ncoupl > ngrid\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            fprintf(fortran_get_unit(iunit11), "  slave   master   xorig   yorig   zorig\n");
            for (int nn = 1; nn <= ncoupl; nn++) {
                int igslav, igmast;
                fortran_read_list(iunit5, &igslav, &igmast, &realval[1], &realval[2], &realval[3]);
                double xorg=realval[1], yorg=realval[2], zorg=realval[3];
                fprintf(fortran_get_unit(iunit11), "%8d%9d%9.4f%9.4f%9.4f\n",
                        igslav, igmast, (float)xorg, (float)yorg, (float)zorg);
                int nblslav = nblg(igslav);
                int nblmast = 0;
                if (igmast > 0) {
                    nblmast = nblg(igmast);
                    xorig(nblmast)=(float)xorg; yorig(nblmast)=(float)yorg; zorig(nblmast)=(float)zorg;
                    xorig0(nblmast)=(float)xorg; yorig0(nblmast)=(float)yorg; zorig0(nblmast)=(float)zorg;
                }
                ncg = ncgg(igslav);
                for (int iseg = 1; iseg <= nsegdfrm(nblslav); iseg++) icouple(nblslav,iseg) = nblmast;
                if (ncg > 0) {
                    for (n = 1; n <= ncg; n++) {
                        nblslav++;
                        if (nblmast > 0) {
                            nblmast++;
                            xorig(nblmast)=xorig(nblmast-1); yorig(nblmast)=yorig(nblmast-1); zorig(nblmast)=zorig(nblmast-1);
                            xorig0(nblmast)=xorig0(nblmast-1); yorig0(nblmast)=yorig0(nblmast-1); zorig0(nblmast)=zorig0(nblmast-1);
                        }
                        for (int iseg = 1; iseg <= nsegdfrm(nblslav); iseg++) icouple(nblslav,iseg) = nblmast;
                    }
                }
            }
        } else if (ncoupl < 0) {
            int igslav0, igmast0;
            fortran_read_list(iunit5, &igslav0, &igmast0, &realval[1], &realval[2], &realval[3]);
            double xorg=realval[1], yorg=realval[2], zorg=realval[3];
            ncoupl = ngrid;
            fprintf(fortran_get_unit(iunit11), "%8d\n", ncoupl);
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            fprintf(fortran_get_unit(iunit11), "  slave   master   xorig   yorig   zorig\n");
            for (int nn = 1; nn <= ncoupl; nn++) {
                int igslav = nn, igmast = igmast0;
                fprintf(fortran_get_unit(iunit11), "%8d%9d%9.4f%9.4f%9.4f\n",
                        igslav, igmast, (float)xorg, (float)yorg, (float)zorg);
                int nblslav = nblg(igslav);
                int nblmast = 0;
                if (igmast > 0) {
                    nblmast = nblg(igmast);
                    xorig(nblmast)=(float)xorg; yorig(nblmast)=(float)yorg; zorig(nblmast)=(float)zorg;
                    xorig0(nblmast)=(float)xorg; yorig0(nblmast)=(float)yorg; zorig0(nblmast)=(float)zorg;
                }
                ncg = ncgg(igslav);
                for (int iseg = 1; iseg <= nsegdfrm(nblslav); iseg++) icouple(nblslav,iseg) = nblmast;
                if (ncg > 0) {
                    for (n = 1; n <= ncg; n++) {
                        nblslav++;
                        if (nblmast > 0) {
                            nblmast++;
                            xorig(nblmast)=xorig(nblmast-1); yorig(nblmast)=yorig(nblmast-1); zorig(nblmast)=zorig(nblmast-1);
                            xorig0(nblmast)=xorig0(nblmast-1); yorig0(nblmast)=yorig0(nblmast-1); zorig0(nblmast)=zorig0(nblmast-1);
                        }
                        for (int iseg = 1; iseg <= nsegdfrm(nblslav); iseg++) icouple(nblslav,iseg) = nblmast;
                    }
                }
            }
        } else {
            fprintf(fortran_get_unit(iunit11), "%8d\n", ncoupl);
            { char _skip[256]; fgets(_skip, sizeof(_skip), fortran_get_unit(iunit5)); }
            fprintf(fortran_get_unit(iunit11), "  slave   master   xorig   yorig   zorig\n");
        }
        // deallocate temporary arrays
        riskp.deallocate(); rjskp.deallocate(); rkskp.deallocate();
        iskipt.deallocate(); jskipt.deallocate(); kskipt.deallocate();
    } // end if iunst > 1 || idef_ss > 0

    // final mxdefseg
    for (nbl = 1; nbl <= nblock; nbl++) mxdefseg = std::max(mxdefseg, nsegdfrm(nbl));


    // check deforming surface segments lie on block faces
    for (nbl = 1; nbl <= nblock; nbl++) {
        int jdim = jdimg(nbl), kdim = kdimg(nbl), idim = idimg(nbl);
        if (idefrm(nbl) < 99) {
            for (int iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {
                int ist=icsi(nbl,iseg), ifn=icsf(nbl,iseg);
                int jst=jcsi(nbl,iseg), jfn=jcsf(nbl,iseg);
                int kst=kcsi(nbl,iseg), kfn=kcsf(nbl,iseg);
                int iok=1, jok=1, kok=1, lconst=1;
                if (ist==ifn) {
                    if (ist==1 || ist==idim) iok=1; else iok=0;
                } else if (jst==jfn) {
                    if (jst==1 || jst==jdim) jok=1; else jok=0;
                } else if (kst==kfn) {
                    if (kst==1 || kst==kdim) kok=1; else kok=0;
                } else {
                    lconst=0;
                }
                if (iok==0 || jok==0 || kok==0)
                    fprintf(fortran_get_unit(iunit11), " stopping...deforming surface does not lie on a block face\n");
                if (lconst==0)
                    fprintf(fortran_get_unit(iunit11), " stopping...at least one index must be constant on a deforming surface\n");
                int isum = iok+jok+kok+lconst;
                if (isum != 4) {
                    fprintf(fortran_get_unit(iunit11), "   segment%3d block %4d:\n", iseg, nbl);
                    fprintf(fortran_get_unit(iunit11), "   icsi,icsf,jcsi,jcsf,kcsi,kcsf %4d%4d%4d%4d%4d%4d\n",
                            ist, ifn, jst, jfn, kst, kfn);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }
        }
    }

    // swap y and z grid motion parameters for ialph > 0
    if (ialph > 0) {
        float temp;
        temp=vtransmc; vtransmc=-wtransmc; wtransmc=temp;
        temp=omegaymc; omegaymc=-omegazmc; omegazmc=temp;
        temp=thetaymc; thetaymc=-thetazmc; thetazmc=temp;
        temp=yorigmc;  yorigmc=-zorigmc;   zorigmc=temp;
        temp=yorig0mc; yorig0mc=-zorig0mc; zorig0mc=temp;
        temp=dymxmc;   dymxmc=-dzmxmc;    dzmxmc=temp;
        temp=dthymxmc; dthymxmc=-dthzmxmc; dthzmxmc=temp;
        temp=ymc0;     ymc0=-zmc0;        zmc0=temp;
        for (nbl = 1; nbl <= nblock; nbl++) {
            double tmp;
            tmp=vtrans(nbl); vtrans(nbl)=-wtrans(nbl); wtrans(nbl)=tmp;
            tmp=dymx(nbl);   dymx(nbl)=-dzmx(nbl);    dzmx(nbl)=tmp;
            tmp=omegay(nbl); omegay(nbl)=-omegaz(nbl); omegaz(nbl)=tmp;
            tmp=yorig(nbl);  yorig(nbl)=-zorig(nbl);   zorig(nbl)=tmp;
            tmp=yorig0(nbl); yorig0(nbl)=-zorig0(nbl); zorig0(nbl)=tmp;
            tmp=thetay(nbl); thetay(nbl)=-thetaz(nbl); thetaz(nbl)=tmp;
            tmp=dthymx(nbl); dthymx(nbl)=-dthzmx(nbl); dthzmx(nbl)=tmp;
            tmp=thetayl(nbl); thetayl(nbl)=-thetazl(nbl); thetazl(nbl)=tmp;
        }
        for (nbl = 1; nbl <= nblock; nbl++) {
            for (int iseg = 1; iseg <= maxsegdg; iseg++) {
                double tmp;
                tmp=vtrnsae(nbl,iseg); vtrnsae(nbl,iseg)=-wtrnsae(nbl,iseg); wtrnsae(nbl,iseg)=tmp;
                tmp=omgyae(nbl,iseg);  omgyae(nbl,iseg)=-omgzae(nbl,iseg);   omgzae(nbl,iseg)=tmp;
                tmp=yorgae(nbl,iseg);  yorgae(nbl,iseg)=zorgae(nbl,iseg);    zorgae(nbl,iseg)=tmp;
                tmp=yorgae0(nbl,iseg); yorgae0(nbl,iseg)=-zorgae0(nbl,iseg); zorgae0(nbl,iseg)=tmp;
                tmp=thtyae(nbl,iseg);  thtyae(nbl,iseg)=-thtyae(nbl,iseg);   thtzae(nbl,iseg)=tmp;
            }
        }
    }

    // print summary by grids
    fprintf(fortran_get_unit(iunit11), "\n SUMMARY BY GRIDS\n");
    fprintf(fortran_get_unit(iunit11), "      grid     level     block       jdim       kdim       idim   grid pts.\n");
    igptot = 0;
    for (igrid = 1; igrid <= ngrid; igrid++) {
        iem = iemg(igrid); nbl = nblg(igrid); ncg = ncgg(igrid);
        igpts = jdimg(nbl)*kdimg(nbl)*idimg(nbl);
        igptot += igpts;
        fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%12d\n",
                igrid, levelg(nbl), nbl, jdimg(nbl), kdimg(nbl), idimg(nbl), igpts);
        if (ncg > 0 && iem == 0) {
            for (n = 1; n <= ncg; n++) {
                nbl++;
                igpts = jdimg(nbl)*kdimg(nbl)*idimg(nbl);
                igptot += igpts;
                fprintf(fortran_get_unit(iunit11), "%10d%10d%10d%10d%10d%10d%12d\n",
                        igrid, levelg(nbl), nbl, jdimg(nbl), kdimg(nbl), idimg(nbl), igpts);
            }
        }
    }
    fprintf(fortran_get_unit(iunit11), "\n                                                       TOTAL%12d\n", igptot);

    // print summary by levels
    fprintf(fortran_get_unit(iunit11), "\n SUMMARY BY LEVELS\n");
    fprintf(fortran_get_unit(iunit11), "     level      grid     block\n");
    lf = lfem;
    if (lfem == 0) lf = lfgm;
    for (int levelc = lf; levelc >= lcgm; levelc--) {
        ncell[levelc-1] = 0;
        for (nbl = 1; nbl <= nblock; nbl++) {
            if (levelc != levelg(nbl)) continue;
            igrid = igridg(nbl);
            ncell[levelc-1] += (jdimg(nbl)-1)*(kdimg(nbl)-1)*(idimg(nbl)-1);
            fprintf(fortran_get_unit(iunit11), "%10d%10d%10d\n", levelg(nbl), igrid, nbl);
        }
    }
    fprintf(fortran_get_unit(iunit11), "\n     level     cells\n");
    for (int levelc = lf; levelc >= lcgm; levelc--) {
        fprintf(fortran_get_unit(iunit11), "%10d%10d\n", levelc, ncell[levelc-1]);
    }
    if (mseq > (ncgmax+1)) {
        fprintf(fortran_get_unit(iunit11), " mseq,ncgmax=%4d%4d\n", mseq, ncgmax);
        fprintf(fortran_get_unit(iunit11), " error in input, mseq, ncgmax%5d%5d\n", mseq, ncgmax);
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }
    if (lfem != 0) fprintf(fortran_get_unit(iunit11), "\n level of finest embedded mesh=%3d\n", lfem);
    fprintf(fortran_get_unit(iunit11), "\n level of coarsest global mesh=%3d\n", lcgm);
    fprintf(fortran_get_unit(iunit11), "\n level of finest global mesh  =%3d\n", lfgm);


    // print summary of grid sequences
    fprintf(fortran_get_unit(iunit11), "\n SUMMARY OF GRID SEQUENCES\n");
    fprintf(fortran_get_unit(iunit11), "  sequence   starting level     ending level  cells(finest global+embeded)\n");
    for (m = 1; m <= mseq; m++) {
        int ncellseq;
        if (nemgl(m) == 0) ncellseq = ncell[levelt[m-1]-1];
        else ncellseq = ncell[levelt[m-1]-nemgl(m)-1] + ncell[levelt[m-1]-1];
        fprintf(fortran_get_unit(iunit11), "%10d%17d%17d%17d\n", m, levelt[m-1], levelb[m-1], ncellseq);
    }

    // validate and expand 1:1 blocking
    if (std::abs(nbli) > 0) {
        nblict = std::abs(nbli);
        for (n = 1; n <= std::abs(nbli); n++) {
            if (nblk(1,n) < 1 || nblk(1,n) > ngrid) {
                fprintf(fortran_get_unit(iunit11), " side number 1 is out of range for 1:1 plane%3d\n", n);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (nblk(2,n) < 1 || nblk(2,n) > ngrid) {
                fprintf(fortran_get_unit(iunit11), " side number 2 is out of range for 1:1 plane%3d\n", n);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (isva(1,1,n) < 1 || isva(1,1,n) > 3) {
                fprintf(fortran_get_unit(iunit11), " ind1 for side 1 is out of range for 1:1 plane%3d\n", n);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (isva(1,2,n) < 1 || isva(1,2,n) > 3) {
                fprintf(fortran_get_unit(iunit11), " ind2 for side 1 is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (isva(1,1,n) == isva(1,2,n)) {
                fprintf(fortran_get_unit(iunit11), " ind1 = ind2 for block 1 for 1:1 plane%3d\n", n);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (isva(2,1,n) < 1 || isva(2,1,n) > 3) {
                fprintf(fortran_get_unit(iunit11), " ind1 for side 2 is out of range for 1:1 plane%3d\n", n);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (isva(2,2,n) < 1 || isva(2,2,n) > 3) {
                fprintf(fortran_get_unit(iunit11), " ind2 for side 2 is out of range\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (isva(2,1,n) == isva(2,2,n)) {
                fprintf(fortran_get_unit(iunit11), " ind1 = ind2 for block 2 for 1:1 plane%3d\n", n);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            for (int ib = 1; ib <= 2; ib++) {
                ntblk = nblg(nblk(ib,n));
                iflag = 0;
                if (limblk(ib,1,n)<1 || limblk(ib,1,n)>idimg(ntblk)) iflag=1;
                if (limblk(ib,2,n)<1 || limblk(ib,2,n)>jdimg(ntblk)) iflag=1;
                if (limblk(ib,3,n)<1 || limblk(ib,3,n)>kdimg(ntblk)) iflag=1;
                if (limblk(ib,4,n)<1 || limblk(ib,4,n)>idimg(ntblk)) iflag=1;
                if (limblk(ib,5,n)<1 || limblk(ib,5,n)>jdimg(ntblk)) iflag=1;
                if (limblk(ib,6,n)<1 || limblk(ib,6,n)>kdimg(ntblk)) iflag=1;
                if (iflag == 1) {
                    fprintf(fortran_get_unit(iunit11), " limits out of range for side%4d for 1:1 plane%3d\n", ib, n);
                    termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                }
            }
            i1b1 = isva(1,1,n); i2b1 = isva(1,2,n);
            i1b2 = isva(2,1,n); i2b2 = isva(2,2,n);
            idif1 = std::abs(limblk(1,i1b1,n) - limblk(1,i1b1+3,n));
            idif2 = std::abs(limblk(1,i2b1,n) - limblk(1,i2b1+3,n));
            idif3 = std::abs(limblk(2,i1b2,n) - limblk(2,i1b2+3,n));
            idif4 = std::abs(limblk(2,i2b2,n) - limblk(2,i2b2+3,n));
            if (idif1 != idif3) {
                fprintf(fortran_get_unit(iunit11), " the limits of ind1 are not the same for both sides for 1:1 plane%3d\n", n);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (idif2 != idif4) {
                fprintf(fortran_get_unit(iunit11), " the limits of ind2 are not the same for both sides for 1:1 plane%3d\n", n);
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ncgg(nblk(1,n)) != ncgg(nblk(2,n))) {
                fprintf(fortran_get_unit(iunit11), " both sides are not at the same level of coarseness\n");
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
            if (ncgg(nblk(1,n)) > 0) {
                for (m = 1; m <= ncgg(nblk(1,n)); m++) {
                    nblict++;
                    nblon(nblict) = nblon(n);
                    for (int ib = 1; ib <= 2; ib++) {
                        nblk(ib,nblict) = nblg(nblk(ib,n)) + m;
                        int lpoint = n;
                        if (m > 1) lpoint = nblict - 1;
                        for (int ind = 1; ind <= 2; ind++) isva(ib,ind,nblict) = isva(ib,ind,n);
                        for (l = 1; l <= 6; l++) {
                            itest1 = (limblk(ib,l,lpoint)-1)/2 + 1;
                            itest2 = limblk(ib,l,lpoint)/2 + 1;
                            if (itest1 != itest2 && itest1 != 1) {
                                fprintf(fortran_get_unit(iunit11),
                                    " one of the points for blocking cannot be mapped to a coarser grid in block %d\n", ib);
                                fprintf(fortran_get_unit(iunit11), " %d %d\n", itest1, itest2);
                                fprintf(fortran_get_unit(iunit11), " %d\n", limblk(ib,l,lpoint));
                                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                            }
                            if (itest1 == 1) itest1 = limblk(ib,l,lpoint);
                            limblk(ib,l,nblict) = itest1;
                        }
                    }
                }
            }
            nblk(1,n) = nblg(nblk(1,n));
            nblk(2,n) = nblg(nblk(2,n));
        } // end do 5000
        if (nbli > 0) nbli =  nblict;
        if (nbli < 0) nbli = -nblict;
    }


    // write summary of 1:1 blocking indices (reset limits to cell center)
    for (n = 1; n <= std::abs(nbli); n++) {
        for (int ib = 1; ib <= 2; ib++) {
            for (m = 1; m <= 3; m++) {
                if (m + isva(ib,1,n) + isva(ib,2,n) != 6) {
                    if (limblk(ib,m,n) < limblk(ib,m+3,n))
                        limblk(ib,m+3,n) = limblk(ib,m+3,n) - 1;
                    else
                        limblk(ib,m,n) = limblk(ib,m,n) - 1;
                }
            }
        }
    }

    // determine more array size requirements
    if (icall == 0) {
        lmaxgr   = std::max(ngrid, 1);
        lmaxbl   = std::max(nblock, 1);
        lmxseg   = std::max(msegment, 1);
        lnplts   = std::max(std::max(nplot3d, nprint), 1);
        lmaxcs   = std::max(ncs2+1, 1);
        lmxbli   = std::max(ntest, 1);
        lncycm   = std::max(ncyctot+ntr, 1);
        nintr    = 0;
        lmptch   = 1;
        lintmx   = 1;
        lmxxe    = 1;
        lmsub1   = 1;
        lnmds    = mxmds;
        lmaxaes  = mxaes;
        lmxsegdg = mxdefseg;
        lnmaster = 21;
    }
} // end global()

void blocking_skip(int& nskip1, int& nskip2, int& nskipt1, int& nskipt2) {
    if (nskip1 > nskip2) {
        if ((float)nskip1/(float)nskip2 == (float)(nskip1/nskip2)) {
            nskipt1 = nskip2; nskipt2 = nskip2;
        } else {
            for (int i = 1; i <= 8; i++) {
                nskip2 = nskip2/2;
                if ((float)nskip1/(float)nskip2 == (float)(nskip1/nskip2)) {
                    nskipt1 = nskip2; nskipt2 = nskip2; return;
                }
            }
        }
    } else {
        if ((float)nskip2/(float)nskip1 == (float)(nskip2/nskip1)) {
            nskipt1 = nskip1; nskipt2 = nskip1;
        } else {
            for (int i = 1; i <= 8; i++) {
                nskip1 = nskip1/2;
                if ((float)nskip2/(float)nskip1 == (float)(nskip2/nskip1)) {
                    nskipt1 = nskip1; nskipt2 = nskip1; return;
                }
            }
        }
    }
}

void skordr(int& nbl, int& nskp, int& i1, FortranArray2DRef<int> ijkskip, int& maxbl) {
    for (int ii1 = nskp; ii1 >= 1; ii1--) {
        if (ijkskip(nbl, ii1) == i1) return;
    }
    nskp++;
    ijkskip(nbl, nskp) = ijkskip(nbl, nskp-1);
    for (int ii1 = nskp-1; ii1 >= 2; ii1--) {
        if (ijkskip(nbl, ii1-1) < i1) {
            ijkskip(nbl, ii1) = i1;
            return;
        } else {
            ijkskip(nbl, ii1) = ijkskip(nbl, ii1-1);
        }
    }
}

void bndchk(int& n, int& nbl, int& iseg, int& n1, int& n2, int& n3,
            int& maxbl, int& maxsegdg, int& maxseg,
            FortranArray2DRef<int> ijkskip, FortranArray1DRef<int> nskpijk1,
            FortranArray4DRef<int> nbcinfo, FortranArray1DRef<int> nbc,
            FortranArray2DRef<int> ncs) {
    for (int nseg = 1; nseg <= nbc(nbl); nseg++) {
        if (nbcinfo(nbl,nseg,1,n3)==2004 || nbcinfo(nbl,nseg,1,n3)==1005 ||
            nbcinfo(nbl,nseg,1,n3)==1006 || nbcinfo(nbl,nseg,1,n3)==2014 ||
            nbcinfo(nbl,nseg,1,n3)==2024 || nbcinfo(nbl,nseg,1,n3)==2034 ||
            nbcinfo(nbl,nseg,1,n3)==2016) {
            if (nbcinfo(nbl,nseg,n1,n3) < ncs(nbl,iseg) &&
                nbcinfo(nbl,nseg,n2,n3) >= ncs(nbl,iseg)) {
                int val = ncs(nbl,iseg)-1;
                skordr(nbl, nskpijk1(n), val, ijkskip, maxbl);
            }
            if (nbcinfo(nbl,nseg,n1,n3) <= ncs(nbl,iseg) &&
                nbcinfo(nbl,nseg,n2,n3) > ncs(nbl,iseg)) {
                int val = ncs(nbl,iseg)+1;
                skordr(nbl, nskpijk1(n), val, ijkskip, maxbl);
            }
        }
    }
}


void cpadd(int& nbl1, int& nbl2, int& nblk1, int& nblk2,
           int& lmblk11, int& lmblk12, int& lmblk21, int& lmblk22,
           int& nskp1, FortranArray1DRef<int> ijkskip1,
           int& nskp2, FortranArray1DRef<int> ijkskip2,
           int& maxbl, int& ijkdim1, int& ijkdim2, int& icrd1, int& icrd2) {
    int lmskp11=0, lmskp12=0, lmskp21=0, lmskp22=0;
    FortranArray2D<int> ijkskpt1(maxbl, 500);
    FortranArray2D<int> ijkskpt2(maxbl, 500);
    for (int nb = 1; nb <= maxbl; nb++)
        for (int ii = 1; ii <= 500; ii++) { ijkskpt1(nb,ii)=0; ijkskpt2(nb,ii)=0; }

    for (int nn = 1; nn <= nskp1; nn++) {
        ijkskpt1(nbl1, nn) = ijkskip1(nn);
        if (ijkskpt1(nbl1,nn) == lmblk11) lmskp11 = nn;
        if (ijkskpt1(nbl1,nn) == lmblk12) lmskp12 = nn;
    }
    for (int nn = 1; nn <= nskp2; nn++) {
        ijkskpt2(nbl2, nn) = ijkskip2(nn);
        if (ijkskpt2(nbl2,nn) == lmblk21) lmskp21 = nn;
        if (ijkskpt2(nbl2,nn) == lmblk22) lmskp22 = nn;
    }

    FortranArray2DRef<int> ijkskpt1_ref = ijkskpt1.ref();
    FortranArray2DRef<int> ijkskpt2_ref = ijkskpt2.ref();

    if (nblk1 == nblk2 && icrd1 == icrd2) {
        if (lmblk11 < lmblk12) {
            if (lmblk21 < lmblk22) {
                for (int n1 = lmskp11; n1 <= lmskp12; n1++) {
                    int itst1 = ijkskip1(n1) - lmblk11 + lmblk21;
                    bool found = false;
                    for (int n2 = lmskp21; n2 <= lmskp22; n2++) { if (ijkskip2(n2)==itst1) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
                for (int n2 = lmskp21; n2 <= lmskp22; n2++) {
                    int itst1 = ijkskip2(n2) - lmblk21 + lmblk11;
                    bool found = false;
                    for (int n1 = lmskp11; n1 <= lmskp12; n1++) { if (itst1==ijkskip1(n1)) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
            } else {
                for (int n1 = lmskp11; n1 <= lmskp12; n1++) {
                    int itst1 = lmblk12 - ijkskip1(n1) + lmblk22;
                    bool found = false;
                    for (int n2 = lmskp22; n2 <= lmskp21; n2++) { if (ijkskip2(n2)==itst1) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
                for (int n2 = lmskp22; n2 <= lmskp21; n2++) {
                    int itst1 = lmblk21 - ijkskip2(n2) + lmblk11;
                    bool found = false;
                    for (int n1 = lmskp11; n1 <= lmskp12; n1++) { if (itst1==ijkskip1(n1)) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
            }
        } else {
            if (lmblk21 < lmblk22) {
                for (int n1 = lmskp12; n1 <= lmskp11; n1++) {
                    int itst1 = lmblk11 - ijkskip1(n1) + lmblk21;
                    bool found = false;
                    for (int n2 = lmskp21; n2 <= lmskp22; n2++) { if (ijkskip2(n2)==itst1) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
                for (int n2 = lmskp21; n2 <= lmskp22; n2++) {
                    int itst1 = lmblk22 - ijkskip2(n2) + lmblk12;
                    bool found = false;
                    for (int n1 = lmskp12; n1 <= lmskp11; n1++) { if (itst1==ijkskip1(n1)) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
            } else {
                for (int n1 = lmskp12; n1 <= lmskp11; n1++) {
                    int itst1 = ijkskip1(n1) - lmblk12 + lmblk22;
                    bool found = false;
                    for (int n2 = lmskp22; n2 <= lmskp21; n2++) { if (ijkskip2(n2)==itst1) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
                for (int n2 = lmskp22; n2 <= lmskp21; n2++) {
                    int itst1 = ijkskip2(n2) - lmblk22 + lmblk12;
                    bool found = false;
                    for (int n1 = lmskp12; n1 <= lmskp11; n1++) { if (itst1==ijkskip1(n1)) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
            }
        }
        nskp1 = nskp2;
        for (int nn = 1; nn <= nskp1; nn++) ijkskip1(nn) = ijkskpt2(nbl2, nn);
        for (int nn = 1; nn <= nskp2; nn++) ijkskip2(nn) = ijkskpt2(nbl2, nn);
    } else {
        if (lmblk11 < lmblk12) {
            if (lmblk21 < lmblk22) {
                for (int n1 = lmskp11; n1 <= lmskp12; n1++) {
                    int itst1 = ijkskip1(n1) - lmblk11 + lmblk21;
                    bool found = false;
                    for (int n2 = lmskp21; n2 <= lmskp22; n2++) { if (ijkskip2(n2)==itst1) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
                for (int n2 = lmskp21; n2 <= lmskp22; n2++) {
                    int itst1 = ijkskip2(n2) - lmblk21 + lmblk11;
                    bool found = false;
                    for (int n1 = lmskp11; n1 <= lmskp12; n1++) { if (itst1==ijkskip1(n1)) { found=true; break; } }
                    if (!found) skordr(nbl1, nskp1, itst1, ijkskpt1_ref, maxbl);
                }
            } else {
                for (int n1 = lmskp11; n1 <= lmskp12; n1++) {
                    int itst1 = lmblk12 - ijkskip1(n1) + lmblk22;
                    bool found = false;
                    for (int n2 = lmskp22; n2 <= lmskp21; n2++) { if (ijkskip2(n2)==itst1) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
                for (int n2 = lmskp22; n2 <= lmskp21; n2++) {
                    int itst1 = lmblk21 - ijkskip2(n2) + lmblk11;
                    bool found = false;
                    for (int n1 = lmskp11; n1 <= lmskp12; n1++) { if (itst1==ijkskip1(n1)) { found=true; break; } }
                    if (!found) skordr(nbl1, nskp1, itst1, ijkskpt1_ref, maxbl);
                }
            }
        } else {
            if (lmblk21 < lmblk22) {
                for (int n1 = lmskp12; n1 <= lmskp11; n1++) {
                    int itst1 = lmblk11 - ijkskip1(n1) + lmblk21;
                    bool found = false;
                    for (int n2 = lmskp21; n2 <= lmskp22; n2++) { if (ijkskip2(n2)==itst1) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
                for (int n2 = lmskp21; n2 <= lmskp22; n2++) {
                    int itst1 = lmblk22 - ijkskip2(n2) + lmblk12;
                    bool found = false;
                    for (int n1 = lmskp12; n1 <= lmskp11; n1++) { if (itst1==ijkskip1(n1)) { found=true; break; } }
                    if (!found) skordr(nbl1, nskp1, itst1, ijkskpt1_ref, maxbl);
                }
            } else {
                for (int n1 = lmskp12; n1 <= lmskp11; n1++) {
                    int itst1 = ijkskip1(n1) - lmblk12 + lmblk22;
                    bool found = false;
                    for (int n2 = lmskp22; n2 <= lmskp21; n2++) { if (ijkskip2(n2)==itst1) { found=true; break; } }
                    if (!found) skordr(nbl2, nskp2, itst1, ijkskpt2_ref, maxbl);
                }
                for (int n2 = lmskp22; n2 <= lmskp21; n2++) {
                    int itst1 = ijkskip2(n2) - lmblk22 + lmblk12;
                    bool found = false;
                    for (int n1 = lmskp12; n1 <= lmskp11; n1++) { if (itst1==ijkskip1(n1)) { found=true; break; } }
                    if (!found) skordr(nbl1, nskp1, itst1, ijkskpt1_ref, maxbl);
                }
            }
        }
        for (int nn = 1; nn <= nskp1; nn++) ijkskip1(nn) = ijkskpt1(nbl1, nn);
        for (int nn = 1; nn <= nskp2; nn++) ijkskip2(nn) = ijkskpt2(nbl2, nn);
    }
}

} // namespace global_ns
