// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "sizer.h"
#include <cstdio>
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include "umalloc.h"
#include "lead.h"
#include "pointers.h"
#include "termn8.h"
#include "global.h"
#include "global2.h"
#include "compg2n.h"
#include "cntsurf.h"
#include "pre_patch.h"
#include "pre_blockbc.h"
#include "pre_period.h"
#include "pre_embed.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace sizer_ns {

void sizer(int& mwork, int& mworki, int& nplots, int& minnode, int& iitot,
           int& intmax, int& maxxe, int& mxbli, int& nsub1, int& lbcprd,
           int& lbcemb, int& lbcrad, int& maxbl, int& maxgr, int& maxseg,
           int& maxcs, int& ncycmax, int& intmx, int& mxxe, int& mptch,
           int& msub1, int& nmds, int& maxaes, int& maxsegdg, int& nnodes,
           int& nslave, int& nmaster, int& myhost, int& myid, int& mycomm,
           int& nplots0, int& maxnode0, int& mxbli0, int& lbcprd0,
           int& lbcemb0, int& lbcrad0, int& maxbl0, int& maxgr0,
           int& maxseg0, int& maxcs0, int& ncycmax0, int& intmax0,
           int& nsub10, int& intmx0, int& mxxe0, int& mptch0, int& msub10,
           int& ibufdim0, int& nbuf0, int& mxbcfil0, int& nmds0,
           int& maxaes0, int& maxsegdg0, int& imode, int& ntr,
           FortranArray1DRef<char[80]> bcfiles,
           FortranArray2DRef<char[120]> bou,
           FortranArray1DRef<int> nou)
{
    // COMMON block aliases
    float& gamma     = cmn_fluid.gamma;
    float& gm1       = cmn_fluid.gm1;
    float& gp1       = cmn_fluid.gp1;
    float& gm1g      = cmn_fluid.gm1g;
    float& gp1g      = cmn_fluid.gp1g;
    float& ggm1      = cmn_fluid.ggm1;
    int32_t& jdim    = cmn_ginfo.jdim;
    int32_t& kdim    = cmn_ginfo.kdim;
    int32_t& idim    = cmn_ginfo.idim;
    int32_t& jj2     = cmn_ginfo.jj2;
    int32_t& kk2     = cmn_ginfo.kk2;
    int32_t& ii2     = cmn_ginfo.ii2;
    int32_t& nblc    = cmn_ginfo.nblc;
    int32_t& nit     = cmn_info.nit;
    int32_t* idiag   = cmn_info.idiag;
    int32_t* levelb  = cmn_info.levelb;
    int32_t& mgflag  = cmn_info.mgflag;
    int32_t& iconsf  = cmn_info.iconsf;
    int32_t& mseq    = cmn_info.mseq;
    int32_t* ncyc1   = cmn_info.ncyc1;
    int32_t* levelt  = cmn_info.levelt;
    int32_t& ngam    = cmn_info.ngam;
    int32_t* nsm     = cmn_info.nsm;
    int32_t& ivmx    = cmn_maxiv.ivmx;
    int32_t& levt    = cmn_mgrd.levt;
    int32_t& mtt     = cmn_mgrd.mtt;
    int32_t& level   = cmn_mgrd.level;
    int32_t& lglobal = cmn_mgrd.lglobal;
    float&   radtodeg = cmn_conversion.radtodeg;
    int32_t& iunit5  = cmn_unit5.iunit5;
    int32_t& iunst   = cmn_unst.iunst;
    int32_t& i2d     = cmn_twod.i2d;
    int32_t* is_blk  = cmn_is_blockbc.is_blk;
    int32_t* ie_blk  = cmn_is_blockbc.ie_blk;
    int32_t* is_prd  = cmn_is_perbc.is_prd;
    int32_t* ie_prd  = cmn_is_perbc.ie_prd;
    int32_t& nbcprd  = cmn_is_perbc.nbcprd;
    int32_t* is_pat  = cmn_is_patch.is_pat;
    int32_t* ie_pat  = cmn_is_patch.ie_pat;
    int32_t* is_emb  = cmn_is_embedbc.is_emb;
    int32_t* ie_emb  = cmn_is_embedbc.ie_emb;
    int32_t& nbcemb  = cmn_is_embedbc.nbcemb;
    char* grid       = cmn_filenam.grid;
    char* plt3dg     = cmn_filenam.plt3dg;
    char* plt3dq     = cmn_filenam.plt3dq;
    char* output     = cmn_filenam.output;
    char* residual   = cmn_filenam.residual;
    char* turbres    = cmn_filenam.turbres;
    char* blomx      = cmn_filenam.blomx;
    char* output2    = cmn_filenam.output2;
    char* printout   = cmn_filenam.printout;
    char* pplunge    = cmn_filenam.pplunge;
    char* ovrlap     = cmn_filenam.ovrlap;
    char* patch      = cmn_filenam.patch;
    char* restrt     = cmn_filenam.restrt;
    int32_t& movie   = cmn_moov.movie;
    int32_t& idef_ss = cmn_elastic_ss.idef_ss;
    int32_t& memadd  = cmn_memory.memadd;
    int32_t& memaddi = cmn_memory.memaddi;
    int32_t& iturbord = cmn_turbconv.iturbord;
    int32_t& isktyp  = cmn_deformz.isktyp;
    int32_t& isarc2d  = cmn_curvat.isarc2d;
    int32_t& ieasmcc2d = cmn_curvat.ieasmcc2d;
    int32_t& isstrc   = cmn_curvat.isstrc;
    int32_t& isarc3d  = cmn_curvat.isarc3d;
    int32_t& iteravg  = cmn_avgdata.iteravg;
    int32_t& ides     = cmn_des.ides;
    int32_t* ivisc    = cmn_reyue.ivisc;
    // /params/ aliases
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
    // Local scalar variables
    int icall, memuse = 0, stats = 0, ierrflg;  // stats: C++ allocate throws on OOM, never sets a code → always 0
    int isum_h;
    int iunit, icallgl, icallpt;
    int ngrid, nblock, nbli, ninter, nplot3d, ip3dsurf, nprint, lfgm, ncs;
    int ihstry, iprnsurf;
    int myiduse, mpihost;
    int nstart_hseq, nstart_hmpi, nstart, nwork, maxl, lembed;
    int nodel, nfiner, iaug, npoints;
    int npts_max, mem_w_max, mem_wk_max, tot_nod;
    int nptsum, npt, nptt, maxpt;
    int need, needi, needi_node;
    int needi_glo, needi_trn, needi_pnt, needi_set, needi_mgblk, needi_qout, needi_yout;
    int itemp, itempt, itemp1, itemp_node;
    int nroomf, nroomi, nroomf0, nroomf_node, nroomi_node;
    int nadd, nadd_node;
    int nworkxs, nworkixs, nwklsmin, nwkireq1, nwkireq2;
    int nroom, nroom_node;
    int ns2004, ns2014, ns2024, ns2034, ns2016, ns1005, ns1006, nsurf;
    int minbox, nbb;
    int nv, maxmem, mem_updatedg;
    int nwfa;
    int icount_pat, icount_blk, icount_prd, icount_emb, icount_dpat, icount_rad;
    int icount_pat1, icount_blk1, icount_prd1, icount_emb1;
    int levl, levb;
    int ivert, itot, jtot, ktot;
    int ioflag, imode2;
    int ifree;
    int ichek;
    int n14_plus;
    int imovie;
    int itempi;
    int iskipz;
    int ibufdim, nbuf;
    double pi;
    double ftot_seq, fmemsum;
    float speedopt;
    int i, j, k, n, nbl, nblf, nblz, nblcc, igrid, igridc;
    int iprint = 0;
    int ii, jj, ll, mbl, nf, ind, nnn, nod;
    int iseq, ntime, nit1, ilc;
    int kxpand;
    int lcnt, intr, myid0;
    int lmax_val, itb, ifb, ityp;
    int ic_blk, in_blk, nd_dest, nd_srce, nd_recv;
    int jface, idimn, jdimn, kdimn, maxdims;
    int jmax1, kmax1;
    int nbll, nblp, nseg, nface, ldata;
    int ista, iend, jsta, jend, ksta, kend;
    int idimp, jdimp, kdimp;
    int idimc, jdimc, kdimc;
    int jdw, kdw, idw;
    int i1, i2, i3, j1, j2, j3, k1, k2, k3;
    int ipw;
    int iskp, jskp, kskp_v;
    int ihead, inmax, inmax1;
    int kloop, jloop, iloop;
    int in;
    int jc, kc, ic;
    int iwk1, iwk2, iwk3, iwk4, iwk5, iwk6, iwk7, iwk8, iwk9, iwk10;
    int iwk11, iwk12, iwk13, iwk14, iwk15, iwk16, iwk17, iwk18, iwk19, iwk20;
    int iwk21, iwk22, iwk23, iwk24, iwk25, iwk26, iwk27, iwk28, iwk29, iwk30;
    int iwk31, iwk32, iwk33, iwk34, iwk35, iwk36, iwk37, iwk38, iwk39, iwk40;
    int iwk41, iwk42, iwk43, iwk44, iwk45, iwk46, iwk47, iwk48, iwk49, iwk50;
    int iwk51, iwk52, iwk53, iwk54, iwk55, iwk56, iwk57, iwk58, iwk59, iwk60;
    int iwk61, iwk62, iwk63, iwk64, iwk65, iwk66, iwk67, iwk68, iwk69, iwk70;
    int iwk71, iwk72, iwk73, iwk74, iwk75, iwk76, iwk77, iwk78, iwk79, iwk80;
    int iwk81, iwk82;
    int inmx, ndim;
    int iex, iex2, iex3;
    int iuns;

    // Allocatable arrays - real (double)
    FortranArray2D<double> aesrfdat;
    FortranArray4D<double> bcvali, bcvalj, bcvalk;
    FortranArray2D<double> damp;
    FortranArray2D<double> dthetx, dthety, dthetz;
    FortranArray1D<double> dthxmx, dthymx, dthzmx;
    FortranArray2D<double> dx, dy, dz;
    FortranArray1D<double> dxmx, dymx, dzmx;
    FortranArray1D<double> etf1_d, etf2_d;
    FortranArray2D<double> freq;
    FortranArray2D<double> gf0;
    FortranArray2D<double> gmass;
    FortranArray2D<double> omgxae, omgyae, omgzae;
    FortranArray1D<double> omegax, omegay, omegaz;
    FortranArray3D<double> perturb;
    FortranArray1D<double> rfreqr, rfreqt;
    FortranArray2D<double> rfrqrae, rfrqtae;
    FortranArray2D<double> rkap0g;
    FortranArray1D<double> thetax, thetaxl, thetay, thetayl, thetaz, thetazl;
    FortranArray2D<double> thtxae, thtyae, thtzae;
    FortranArray1D<double> time2;
    FortranArray1D<double> utrans, vtrans, wtrans;
    FortranArray2D<double> utrnsae, vtrnsae, wtrnsae;
    FortranArray2D<double> x0;
    FortranArray2D<double> xorgae, xorgae0, yorgae, yorgae0, zorgae, zorgae0;
    FortranArray1D<double> xorig, xorig0, yorig, yorig0, zorig, zorig0;
    // Allocatable arrays - integer
    FortranArray3D<int> bcfilei, bcfilej, bcfilek;
    FortranArray1D<int> iadvance;
    FortranArray2D<int> iaesurf;
    FortranArray4D<int> ibcinfo, jbcinfo, kbcinfo;
    FortranArray2D<int> ibpntsg;
    FortranArray2D<int> icouple;
    FortranArray2D<int> icsf, icsi;
    FortranArray2D<int> icsinfo;
    FortranArray1D<int> idefrm;
    FortranArray2D<int> idegg;
    FortranArray2D<int> idfrmseg;
    FortranArray2D<int> idiagg;
    FortranArray1D<int> idimg;
    FortranArray1D<int> ieg;
    FortranArray1D<int> iemg;
    FortranArray2D<int> ifdsg;
    FortranArray1D<int> ifiner;
    FortranArray2D<int> iflimg;
    FortranArray1D<int> iforce;
    FortranArray1D<int> ifrom;
    FortranArray1D<int> igridg;
    FortranArray1D<int> iic0;
    FortranArray1D<int> iifit;
    FortranArray1D<int> iiint1, iiint2;
    FortranArray2D<int> iindex;
    FortranArray2D<int> iindx;
    FortranArray1D<int> iiorph;
    FortranArray1D<int> iipntsg;
    FortranArray1D<int> iitmax;
    FortranArray1D<int> iitoss;
    FortranArray1D<int> ilamhig, ilamlog;
    FortranArray1D<int> inewgg;
    FortranArray2D<int> inpl3d;
    FortranArray2D<int> inpr;
    FortranArray1D<int> iovrlp;
    FortranArray2D<int> ipl3dtmp;
    FortranArray1D<int> irotat;
    FortranArray2D<int> isav_blk;
    FortranArray2D<int> isav_dpat;
    FortranArray3D<int> isav_dpat_b;
    FortranArray2D<int> isav_emb;
    FortranArray2D<int> isav_pat;
    FortranArray3D<int> isav_pat_b;
    FortranArray2D<int> isav_prd;
    FortranArray1D<int> isg;
    FortranArray2D<int> iskip;
    FortranArray3D<int> isva;
    FortranArray1D<int> itrans;
    FortranArray1D<int> iv;
    FortranArray2D<int> iviscg;
    FortranArray2D<int> iwfg;
    FortranArray2D<int> iwrk;
    FortranArray2D<int> jcsf, jcsi;
    FortranArray1D<int> jdimg;
    FortranArray1D<int> jeg;
    FortranArray1D<int> jjmax1;
    FortranArray1D<int> jlamhig, jlamlog;
    FortranArray1D<int> jsg;
    FortranArray2D<int> jskip;
    FortranArray2D<int> kcsf, kcsi;
    FortranArray1D<int> kdimg;
    FortranArray1D<int> keg;
    FortranArray1D<int> kkmax1;
    FortranArray1D<int> klamhig, klamlog;
    FortranArray1D<int> ksg;
    FortranArray2D<int> kskip;
    FortranArray1D<int> lbg;
    FortranArray1D<int> levelg;
    FortranArray1D<int> lig;
    FortranArray3D<int> limblk;
    FortranArray1D<int> llimit;
    FortranArray1D<int> ltot;
    FortranArray2D<int> lw;
    FortranArray2D<int> lw2;
    FortranArray3D<int> lwdat;
    FortranArray1D<int> mblk2nd;
    FortranArray1D<int> mblkpt;
    FortranArray1D<int> mem_req_node;
    FortranArray1D<int> memblock;
    FortranArray1D<int> mglevg;
    FortranArray2D<int> mit;
    FortranArray1D<int> mmceta, mmcxie;
    FortranArray1D<int> n14;
    FortranArray1D<int> nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim;
    FortranArray1D<int> nblcg;
    FortranArray1D<int> nblfine;
    FortranArray1D<int> nblg;
    FortranArray2D<int> nblk;
    FortranArray1D<int> nblon;
    FortranArray1D<int> ncgg;
    FortranArray1D<int> ncheck;
    FortranArray1D<int> nemgl;
    FortranArray1D<int> no_of_points;
    FortranArray1D<int> nsegdfrm;
    FortranArray1D<int> xif1, xif2, etf1, etf2;

    // Helper lambda for umalloc calls
    auto do_umalloc = [&](int n_words, int intflag, int& mem, int& st) {
        char txt[1] = {0};
        int nw = n_words, iflag = intflag;
        umalloc_ns::umalloc(nw, iflag, txt, mem, st);
    };


    // Allocate all arrays
    icall  = 0;
    memuse = 0;

    aesrfdat.allocate(5, maxaes0);
    do_umalloc(5*maxaes0, 0, memuse, stats);
    bcfilei.allocate(maxbl0, maxseg0, 2);
    do_umalloc(maxbl0*maxseg0*2, 1, memuse, stats);
    bcfilej.allocate(maxbl0, maxseg0, 2);
    do_umalloc(maxbl0*maxseg0*2, 1, memuse, stats);
    bcfilek.allocate(maxbl0, maxseg0, 2);
    do_umalloc(maxbl0*maxseg0*2, 1, memuse, stats);
    bcvali.allocate(maxbl0, maxseg0, 12, 2);
    do_umalloc(maxbl0*maxseg0*12*2, 0, memuse, stats);
    bcvalj.allocate(maxbl0, maxseg0, 12, 2);
    do_umalloc(maxbl0*maxseg0*12*2, 0, memuse, stats);
    bcvalk.allocate(maxbl0, maxseg0, 12, 2);
    do_umalloc(maxbl0*maxseg0*12*2, 0, memuse, stats);
    damp.allocate(nmds0, maxaes0);
    do_umalloc(nmds0*maxaes0, 0, memuse, stats);
    dthetx.allocate(intmx0, msub10);
    do_umalloc(intmx0*msub10, 0, memuse, stats);
    dthety.allocate(intmx0, msub10);
    do_umalloc(intmx0*msub10, 0, memuse, stats);
    dthetz.allocate(intmx0, msub10);
    do_umalloc(intmx0*msub10, 0, memuse, stats);
    dthxmx.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    dthymx.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    dthzmx.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    dx.allocate(intmx0, msub10);
    do_umalloc(intmx0*msub10, 0, memuse, stats);
    dxmx.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    dy.allocate(intmx0, msub10);
    do_umalloc(intmx0*msub10, 0, memuse, stats);
    dymx.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    dz.allocate(intmx0, msub10);
    do_umalloc(intmx0*msub10, 0, memuse, stats);
    dzmx.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    etf1.allocate(msub10);
    do_umalloc(msub10, 1, memuse, stats);
    etf2.allocate(msub10);
    do_umalloc(msub10, 1, memuse, stats);
    freq.allocate(nmds0, maxaes0);
    do_umalloc(nmds0*maxaes0, 0, memuse, stats);
    gf0.allocate(2*nmds0, maxaes0);
    do_umalloc(2*nmds0*maxaes0, 0, memuse, stats);
    gmass.allocate(nmds0, maxaes0);
    do_umalloc(nmds0*maxaes0, 0, memuse, stats);
    iadvance.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    iaesurf.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 1, memuse, stats);
    ibcinfo.allocate(maxbl0, maxseg0, 7, 2);
    do_umalloc(maxbl0*maxseg0*7*2, 1, memuse, stats);
    ibpntsg.allocate(maxbl0, 4);
    do_umalloc(maxbl0*4, 1, memuse, stats);
    icouple.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 1, memuse, stats);
    icsf.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 1, memuse, stats);
    icsi.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 1, memuse, stats);
    icsinfo.allocate(maxcs0, 9);
    do_umalloc(maxcs0*9, 1, memuse, stats);
    idefrm.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    idegg.allocate(maxbl0, 3);
    do_umalloc(maxbl0*3, 1, memuse, stats);
    idfrmseg.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 1, memuse, stats);
    idiagg.allocate(maxbl0, 3);
    do_umalloc(maxbl0*3, 1, memuse, stats);
    idimg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    ieg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    iemg.allocate(maxgr0);
    do_umalloc(maxgr0, 1, memuse, stats);
    ifdsg.allocate(maxbl0, 3);
    do_umalloc(maxbl0*3, 1, memuse, stats);
    ifiner.allocate(intmx0);
    do_umalloc(intmx0, 1, memuse, stats);
    iflimg.allocate(maxbl0, 3);
    do_umalloc(maxbl0*3, 1, memuse, stats);
    iforce.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    ifrom.allocate(msub10);
    do_umalloc(msub10, 1, memuse, stats);
    igridg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    iic0.allocate(intmx0);
    do_umalloc(intmx0, 1, memuse, stats);
    iifit.allocate(intmx0);
    do_umalloc(intmx0, 1, memuse, stats);
    iiint1.allocate(nsub10);
    do_umalloc(nsub10, 1, memuse, stats);
    iiint2.allocate(nsub10);
    do_umalloc(nsub10, 1, memuse, stats);
    iindex.allocate(intmax0, 6*nsub10+9);
    do_umalloc(intmax0*(6*nsub10+9), 1, memuse, stats);
    iindx.allocate(intmx0, 6*msub10+9);
    do_umalloc(intmx0*(6*msub10+9), 1, memuse, stats);
    iiorph.allocate(intmx0);
    do_umalloc(intmx0, 1, memuse, stats);
    iipntsg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    iitmax.allocate(intmx0);
    do_umalloc(intmx0, 1, memuse, stats);
    iitoss.allocate(intmx0);
    do_umalloc(intmx0, 1, memuse, stats);
    ilamhig.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    ilamlog.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    inewgg.allocate(maxgr0);
    do_umalloc(maxgr0, 1, memuse, stats);
    inpl3d.allocate(nplots0, 11);
    do_umalloc(nplots0*11, 1, memuse, stats);
    inpr.allocate(nplots0, 11);
    do_umalloc(nplots0*11, 1, memuse, stats);
    iovrlp.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    ipl3dtmp.allocate(11, nplots0);
    do_umalloc(11*nplots0, 1, memuse, stats);
    irotat.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);


    isav_blk.allocate(2*mxbli0, 17);
    do_umalloc(2*mxbli0*17, 1, memuse, stats);
    isav_dpat.allocate(intmx0, 17);
    do_umalloc(intmx0*17, 1, memuse, stats);
    isav_dpat_b.allocate(intmx0, msub10, 6);
    do_umalloc(intmx0*msub10*6, 1, memuse, stats);
    isav_emb.allocate(lbcemb0, 12);
    do_umalloc(lbcemb0*12, 1, memuse, stats);
    isav_pat.allocate(intmax0, 17);
    do_umalloc(intmax0*17, 1, memuse, stats);
    isav_pat_b.allocate(intmax0, nsub10, 6);
    do_umalloc(intmax0*nsub10*6, 1, memuse, stats);
    isav_prd.allocate(lbcprd0, 12);
    do_umalloc(lbcprd0*12, 1, memuse, stats);
    isg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    iskip.allocate(maxbl0, 500);
    do_umalloc(500*maxbl0, 1, memuse, stats);
    isva.allocate(2, 2, mxbli0);
    do_umalloc(2*2*mxbli0, 1, memuse, stats);
    itrans.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    iv.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    iviscg.allocate(maxbl0, 3);
    do_umalloc(maxbl0*3, 1, memuse, stats);
    iwfg.allocate(maxbl0, 3);
    do_umalloc(maxbl0*3, 1, memuse, stats);
    iwrk.allocate(maxbl0, 3);
    do_umalloc(maxbl0*3, 1, memuse, stats);
    jbcinfo.allocate(maxbl0, maxseg0, 7, 2);
    do_umalloc(maxbl0*maxseg0*7*2, 1, memuse, stats);
    jcsf.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 1, memuse, stats);
    jcsi.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 1, memuse, stats);
    jdimg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    jeg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    jjmax1.allocate(nsub10);
    do_umalloc(nsub10, 1, memuse, stats);
    jlamhig.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    jlamlog.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    jsg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    jskip.allocate(maxbl0, 500);
    do_umalloc(500*maxbl0, 1, memuse, stats);
    kbcinfo.allocate(maxbl0, maxseg0, 7, 2);
    do_umalloc(maxbl0*maxseg0*7*2, 1, memuse, stats);
    kcsf.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 1, memuse, stats);
    kcsi.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 1, memuse, stats);
    kdimg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    keg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    kkmax1.allocate(nsub10);
    do_umalloc(nsub10, 1, memuse, stats);
    klamhig.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    klamlog.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    ksg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    kskip.allocate(maxbl0, 500);
    do_umalloc(500*maxbl0, 1, memuse, stats);
    lbg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    levelg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    lig.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    limblk.allocate(2, 6, mxbli0);
    do_umalloc(2*6*mxbli0, 1, memuse, stats);
    llimit.allocate(intmx0);
    do_umalloc(intmx0, 1, memuse, stats);
    ltot.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    lw.allocate(65, maxbl0);
    do_umalloc(65*maxbl0, 1, memuse, stats);
    lw2.allocate(43, maxbl0);
    do_umalloc(43*maxbl0, 1, memuse, stats);
    lwdat.allocate(maxbl0, maxseg0, 6);
    do_umalloc(maxbl0*maxseg0*6, 1, memuse, stats);
    mblk2nd.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    mblkpt.allocate(mxxe0);
    do_umalloc(mxxe0, 1, memuse, stats);
    mem_req_node.allocate(maxnode0);
    do_umalloc(maxnode0, 1, memuse, stats);
    memblock.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    mglevg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    mit.allocate(5, maxbl0);
    do_umalloc(5*maxbl0, 1, memuse, stats);
    mmceta.allocate(intmx0);
    do_umalloc(intmx0, 1, memuse, stats);
    mmcxie.allocate(intmx0);
    do_umalloc(intmx0, 1, memuse, stats);
    n14.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    nbci0.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    nbcidim.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    nbcj0.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    nbcjdim.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    nbck0.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    nbckdim.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    nblcg.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    nblfine.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    nblg.allocate(maxgr0);
    do_umalloc(maxgr0, 1, memuse, stats);
    nblk.allocate(2, mxbli0);
    do_umalloc(2*mxbli0, 1, memuse, stats);
    nblon.allocate(mxbli0);
    do_umalloc(mxbli0, 1, memuse, stats);
    ncgg.allocate(maxgr0);
    do_umalloc(maxgr0, 1, memuse, stats);
    ncheck.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    nemgl.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);
    no_of_points.allocate(maxnode0);
    do_umalloc(maxnode0, 1, memuse, stats);
    nsegdfrm.allocate(maxbl0);
    do_umalloc(maxbl0, 1, memuse, stats);


    omegax.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    omegay.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    omegaz.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    omgxae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    omgyae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    omgzae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    perturb.allocate(nmds0, maxaes0, 4);
    do_umalloc(nmds0*maxaes0*4, 0, memuse, stats);
    rfreqr.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    rfreqt.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    rfrqrae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    rfrqtae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    rkap0g.allocate(maxbl0, 3);
    do_umalloc(maxbl0*3, 0, memuse, stats);
    thetax.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    thetaxl.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    thetay.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    thetayl.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    thetaz.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    thetazl.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    thtxae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    thtyae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    thtzae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    time2.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    utrans.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    utrnsae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    vtrans.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    vtrnsae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    wtrans.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    wtrnsae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    x0.allocate(2*nmds0, maxaes0);
    do_umalloc(2*nmds0*maxaes0, 0, memuse, stats);
    xif1.allocate(msub10);
    do_umalloc(msub10, 1, memuse, stats);
    xif2.allocate(msub10);
    do_umalloc(msub10, 1, memuse, stats);
    xorgae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    xorgae0.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    xorig.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    xorig0.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    yorgae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    yorgae0.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    yorig.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    yorig0.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    zorgae.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    zorgae0.allocate(maxbl0, maxsegdg0);
    do_umalloc(maxbl0*maxsegdg0, 0, memuse, stats);
    zorig.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);
    zorig0.allocate(maxbl0);
    do_umalloc(maxbl0, 0, memuse, stats);

    // Now read the Fortran source lines 668 onwards
    // pi constant
    pi = 4.0 * std::atan(1.0);

    // ierrflg = -99 to indicate errors occur during memory allocation phase
    ierrflg = -99;

    // default to dimensions of 1
    lmaxgr   = 1;
    lmaxbl   = 1;
    lmxseg   = 1;
    lmaxcs   = 1;
    lnplts   = 1;
    lmxbli   = 1;
    lmaxxe   = 1;
    lnsub1   = 1;
    lintmx   = 1;
    lmxxe    = 1;
    liitot   = 1;
    isum     = 1;
    lncycm   = 1;
    isum_h   = 1;
    isum_n   = 1;
    lminnode = 1;
    isumi    = 1;
    isumi_n  = 1;
    lmptch   = 1;
    lmsub1   = 1;
    lintmax  = 1;
    libufdim = 1;
    lnbuf    = 1;
    llbcprd  = 1;
    llbcemb  = 1;
    llbcrad  = 1;
    lmaxaes  = 1;
    lnmds    = 1;
    lnslave  = 1;
    lnmaster = 1;
    lmxsegdg = 1;


    // Read filenames from iunit5
    {
        FILE* fu = fortran_get_unit(iunit5);
        char linebuf[256];  // >81: records are 80 chars + newline; fgets must consume the newline
        // skip one line
        fgets(linebuf, sizeof(linebuf), fu);
        // read grid
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(grid,' ',80); std::memcpy(grid,linebuf,std::min(len,60)); }
        // read plt3dg
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(plt3dg,' ',80); std::memcpy(plt3dg,linebuf,std::min(len,60)); }
        // read plt3dq
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(plt3dq,' ',80); std::memcpy(plt3dq,linebuf,std::min(len,60)); }
        // read output
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(output,' ',80); std::memcpy(output,linebuf,std::min(len,60)); }
        // read residual
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(residual,' ',80); std::memcpy(residual,linebuf,std::min(len,60)); }
        // read turbres
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(turbres,' ',80); std::memcpy(turbres,linebuf,std::min(len,60)); }
        // read blomx
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(blomx,' ',80); std::memcpy(blomx,linebuf,std::min(len,60)); }
        // read output2
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(output2,' ',80); std::memcpy(output2,linebuf,std::min(len,60)); }
        // read printout
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(printout,' ',80); std::memcpy(printout,linebuf,std::min(len,60)); }
        // read pplunge
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(pplunge,' ',80); std::memcpy(pplunge,linebuf,std::min(len,60)); }
        // read ovrlap
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(ovrlap,' ',80); std::memcpy(ovrlap,linebuf,std::min(len,60)); }
        // read patch
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(patch,' ',80); std::memcpy(patch,linebuf,std::min(len,60)); }
        // read restrt
        fgets(linebuf, sizeof(linebuf), fu);
        { int len=(int)strlen(linebuf); if(len>0&&linebuf[len-1]=='\n'){linebuf[len-1]=' ';len--;} std::memset(restrt,' ',80); std::memcpy(restrt,linebuf,std::min(len,60)); }
    }

    // Initialize fluid constants
    gamma    = 1.4f;
    gm1      = gamma - 1.0f;
    gp1      = gamma + 1.0f;
    gm1g     = gm1 / gamma;
    gp1g     = gp1 / gamma;
    ggm1     = gamma * gm1;
    pi       = 4.0 * std::atan(1.0);
    radtodeg = (float)(180.0 / pi);

    // Output banner
    FILE* f66 = fortran_get_unit(66);
    fprintf(f66, " ******************************************************************************\n");
    fprintf(f66, " ******************************************************************************\n");
    fprintf(f66, " **                                                                          **\n");
    fprintf(f66, " **                  PRECFL3D - CFL3D PREPROCESSOR           **\n");
    fprintf(f66, " **                                                                          **\n");
    fprintf(f66, " **   VERSION 6.7 :  Computational Fluids Lab, Mail Stop 128,    **\n");
    fprintf(f66, " **                  NASA Langley Research Center, Hampton, VA   **\n");
    fprintf(f66, " **                  Release Date:  February  1, 2017.           **\n");
    fprintf(f66, " **                                                                          **\n");
    fprintf(f66, " ******************************************************************************\n");
    fprintf(f66, " ******************************************************************************\n");

    fprintf(f66, "\n memory allocation: %12.6f Mbytes, single precision\n",
            (float)memuse / 1.e6f);

    fprintf(f66, "\ninput/output files:\n");
    fprintf(f66, "  %.60s\n", grid);
    fprintf(f66, "  %.60s\n", plt3dg);
    fprintf(f66, "  %.60s\n", plt3dq);
    fprintf(f66, "  %.60s\n", output);
    fprintf(f66, "  %.60s\n", residual);
    fprintf(f66, "  %.60s\n", turbres);
    fprintf(f66, "  %.60s\n", blomx);
    fprintf(f66, "  %.60s\n", output2);
    fprintf(f66, "  %.60s\n", printout);
    fprintf(f66, "  %.60s\n", pplunge);
    fprintf(f66, "  %.60s\n", ovrlap);
    fprintf(f66, "  %.60s\n", patch);
    fprintf(f66, "  %.60s\n", restrt);


    // Read global information
    iunit = 66;
    for (ll = 1; ll <= nbuf0; ll++) {
        nou(ll) = 0;
    }
    icallgl = 0;
    {
        FortranArray1DRef<int> nbci0_ref = nbci0.ref();
        FortranArray1DRef<int> nbcj0_ref = nbcj0.ref();
        FortranArray1DRef<int> nbck0_ref = nbck0.ref();
        FortranArray1DRef<int> nbcidim_ref = nbcidim.ref();
        FortranArray1DRef<int> nbcjdim_ref = nbcjdim.ref();
        FortranArray1DRef<int> nbckdim_ref = nbckdim.ref();
        FortranArray4DRef<int> ibcinfo_ref = ibcinfo.ref();
        FortranArray4DRef<int> jbcinfo_ref = jbcinfo.ref();
        FortranArray4DRef<int> kbcinfo_ref = kbcinfo.ref();
        FortranArray3DRef<int> bcfilei_ref = bcfilei.ref();
        FortranArray3DRef<int> bcfilej_ref = bcfilej.ref();
        FortranArray3DRef<int> bcfilek_ref = bcfilek.ref();
        FortranArray2DRef<int> nblk_ref = nblk.ref();
        FortranArray3DRef<int> limblk_ref = limblk.ref();
        FortranArray3DRef<int> isva_ref = isva.ref();
        FortranArray1DRef<int> nblon_ref = nblon.ref();
        FortranArray2DRef<double> rkap0g_ref = rkap0g.ref();
        FortranArray1DRef<int> levelg_ref = levelg.ref();
        FortranArray1DRef<int> igridg_ref = igridg.ref();
        FortranArray2DRef<int> iflimg_ref = iflimg.ref();
        FortranArray2DRef<int> ifdsg_ref = ifdsg.ref();
        FortranArray2DRef<int> iviscg_ref = iviscg.ref();
        FortranArray1DRef<int> jdimg_ref = jdimg.ref();
        FortranArray1DRef<int> kdimg_ref = kdimg.ref();
        FortranArray1DRef<int> idimg_ref = idimg.ref();
        FortranArray2DRef<int> idiagg_ref = idiagg.ref();
        FortranArray1DRef<int> nblcg_ref = nblcg.ref();
        FortranArray2DRef<int> idegg_ref = idegg.ref();
        FortranArray1DRef<int> jsg_ref = jsg.ref();
        FortranArray1DRef<int> ksg_ref = ksg.ref();
        FortranArray1DRef<int> isg_ref = isg.ref();
        FortranArray1DRef<int> jeg_ref = jeg.ref();
        FortranArray1DRef<int> keg_ref = keg.ref();
        FortranArray1DRef<int> ieg_ref = ieg.ref();
        FortranArray2DRef<int> mit_ref = mit.ref();
        FortranArray1DRef<int> ilamlog_ref = ilamlog.ref();
        FortranArray1DRef<int> ilamhig_ref = ilamhig.ref();
        FortranArray1DRef<int> jlamlog_ref = jlamlog.ref();
        FortranArray1DRef<int> jlamhig_ref = jlamhig.ref();
        FortranArray1DRef<int> klamlog_ref = klamlog.ref();
        FortranArray1DRef<int> klamhig_ref = klamhig.ref();
        FortranArray2DRef<int> iwfg_ref = iwfg.ref();
        FortranArray1DRef<double> utrans_ref = utrans.ref();
        FortranArray1DRef<double> vtrans_ref = vtrans.ref();
        FortranArray1DRef<double> wtrans_ref = wtrans.ref();
        FortranArray1DRef<double> omegax_ref = omegax.ref();
        FortranArray1DRef<double> omegay_ref = omegay.ref();
        FortranArray1DRef<double> omegaz_ref = omegaz.ref();
        FortranArray1DRef<double> xorig_ref = xorig.ref();
        FortranArray1DRef<double> yorig_ref = yorig.ref();
        FortranArray1DRef<double> zorig_ref = zorig.ref();
        FortranArray1DRef<double> dxmx_ref = dxmx.ref();
        FortranArray1DRef<double> dymx_ref = dymx.ref();
        FortranArray1DRef<double> dzmx_ref = dzmx.ref();
        FortranArray1DRef<double> dthxmx_ref = dthxmx.ref();
        FortranArray1DRef<double> dthymx_ref = dthymx.ref();
        FortranArray1DRef<double> dthzmx_ref = dthzmx.ref();
        FortranArray1DRef<double> thetax_ref = thetax.ref();
        FortranArray1DRef<double> thetay_ref = thetay.ref();
        FortranArray1DRef<double> thetaz_ref = thetaz.ref();
        FortranArray1DRef<double> rfreqt_ref = rfreqt.ref();
        FortranArray1DRef<double> rfreqr_ref = rfreqr.ref();
        FortranArray1DRef<double> xorig0_ref = xorig0.ref();
        FortranArray1DRef<double> yorig0_ref = yorig0.ref();
        FortranArray1DRef<double> zorig0_ref = zorig0.ref();
        FortranArray1DRef<int> itrans_ref = itrans.ref();
        FortranArray1DRef<int> irotat_ref = irotat.ref();
        FortranArray1DRef<int> idefrm_ref = idefrm.ref();
        FortranArray1DRef<int> ncgg_ref = ncgg.ref();
        FortranArray1DRef<int> nblg_ref = nblg.ref();
        FortranArray1DRef<int> iemg_ref = iemg.ref();
        FortranArray1DRef<int> inewgg_ref = inewgg.ref();
        FortranArray1DRef<int> iovrlp_ref = iovrlp.ref();
        FortranArray2DRef<int> inpl3d_ref = inpl3d.ref();
        FortranArray2DRef<int> inpr_ref = inpr.ref();
        FortranArray1DRef<int> iadvance_ref = iadvance.ref();
        FortranArray1DRef<int> iforce_ref = iforce.ref();
        FortranArray2DRef<int> icsinfo_ref = icsinfo.ref();
        FortranArray1DRef<int> iv_ref = iv.ref();
        FortranArray1DRef<double> time2_ref = time2.ref();
        FortranArray1DRef<double> thetaxl_ref = thetaxl.ref();
        FortranArray1DRef<double> thetayl_ref = thetayl.ref();
        FortranArray1DRef<double> thetazl_ref = thetazl.ref();
        FortranArray2DRef<int> iindex_ref = iindex.ref();
        FortranArray1DRef<int> lig_ref = lig.ref();
        FortranArray1DRef<int> lbg_ref = lbg.ref();
        FortranArray2DRef<int> ibpntsg_ref = ibpntsg.ref();
        FortranArray1DRef<int> iipntsg_ref = iipntsg.ref();
        FortranArray1DRef<int> mglevg_ref = mglevg.ref();
        FortranArray1DRef<int> nemgl_ref = nemgl.ref();
        FortranArray2DRef<int> ipl3dtmp_ref = ipl3dtmp.ref();
        FortranArray2DRef<double> utrnsae_ref = utrnsae.ref();
        FortranArray2DRef<double> vtrnsae_ref = vtrnsae.ref();
        FortranArray2DRef<double> wtrnsae_ref = wtrnsae.ref();
        FortranArray2DRef<double> omgxae_ref = omgxae.ref();
        FortranArray2DRef<double> omgyae_ref = omgyae.ref();
        FortranArray2DRef<double> omgzae_ref = omgzae.ref();
        FortranArray2DRef<double> xorgae_ref = xorgae.ref();
        FortranArray2DRef<double> yorgae_ref = yorgae.ref();
        FortranArray2DRef<double> zorgae_ref = zorgae.ref();
        FortranArray2DRef<double> thtxae_ref = thtxae.ref();
        FortranArray2DRef<double> thtyae_ref = thtyae.ref();
        FortranArray2DRef<double> thtzae_ref = thtzae.ref();
        FortranArray2DRef<double> rfrqtae_ref = rfrqtae.ref();
        FortranArray2DRef<double> rfrqrae_ref = rfrqrae.ref();
        FortranArray2DRef<int> icsi_ref = icsi.ref();
        FortranArray2DRef<int> icsf_ref = icsf.ref();
        FortranArray2DRef<int> jcsi_ref = jcsi.ref();
        FortranArray2DRef<int> jcsf_ref = jcsf.ref();
        FortranArray2DRef<int> kcsi_ref = kcsi.ref();
        FortranArray2DRef<int> kcsf_ref = kcsf.ref();
        FortranArray2DRef<double> freq_ref = freq.ref();
        FortranArray2DRef<double> gmass_ref = gmass.ref();
        FortranArray2DRef<double> damp_ref = damp.ref();
        FortranArray2DRef<double> x0_ref = x0.ref();
        FortranArray2DRef<double> gf0_ref = gf0.ref();
        FortranArray2DRef<double> aesrfdat_ref = aesrfdat.ref();
        FortranArray3DRef<double> perturb_ref = perturb.ref();
        FortranArray2DRef<int> iskip_ref = iskip.ref();
        FortranArray2DRef<int> jskip_ref = jskip.ref();
        FortranArray2DRef<int> kskip_ref = kskip.ref();
        FortranArray1DRef<int> nsegdfrm_ref = nsegdfrm.ref();
        FortranArray2DRef<int> idfrmseg_ref = idfrmseg.ref();
        FortranArray2DRef<int> iaesurf_ref = iaesurf.ref();
        FortranArray2DRef<double> xorgae0_ref = xorgae0.ref();
        FortranArray2DRef<double> yorgae0_ref = yorgae0.ref();
        FortranArray2DRef<double> zorgae0_ref = zorgae0.ref();
        FortranArray2DRef<int> icouple_ref = icouple.ref();
        FortranArray4DRef<double> bcvali_ref = bcvali.ref();
        FortranArray4DRef<double> bcvalj_ref = bcvalj.ref();
        FortranArray4DRef<double> bcvalk_ref = bcvalk.ref();
        int nintr = 0;
        global_ns::global(myid, maxbl0, maxgr0, maxseg0, maxcs0, nplots0, mxbli0,
                          bcvali_ref, bcvalj_ref, bcvalk_ref, nbci0_ref, nbcj0_ref, nbck0_ref,
                          nbcidim_ref, nbcjdim_ref, nbckdim_ref, ibcinfo_ref, jbcinfo_ref,
                          kbcinfo_ref, bcfilei_ref, bcfilej_ref, bcfilek_ref, nblk_ref, nbli,
                          limblk_ref, isva_ref, nblon_ref, rkap0g_ref, nblock, levelg_ref,
                          igridg_ref, iflimg_ref, ifdsg_ref, iviscg_ref, jdimg_ref, kdimg_ref,
                          idimg_ref, idiagg_ref, nblcg_ref, idegg_ref, jsg_ref, ksg_ref, isg_ref,
                          jeg_ref, keg_ref, ieg_ref, mit_ref, ilamlog_ref, ilamhig_ref,
                          jlamlog_ref, jlamhig_ref, klamlog_ref, klamhig_ref, iwfg_ref,
                          utrans_ref, vtrans_ref, wtrans_ref, omegax_ref, omegay_ref, omegaz_ref,
                          xorig_ref, yorig_ref, zorig_ref, dxmx_ref, dymx_ref, dzmx_ref,
                          dthxmx_ref, dthymx_ref, dthzmx_ref, thetax_ref, thetay_ref, thetaz_ref,
                          rfreqt_ref, rfreqr_ref, xorig0_ref, yorig0_ref, zorig0_ref,
                          itrans_ref, irotat_ref, idefrm_ref, ngrid, ncgg_ref, nblg_ref,
                          iemg_ref, inewgg_ref, iovrlp_ref, ninter, nplot3d, inpl3d_ref,
                          ip3dsurf, nprint, inpr_ref, iadvance_ref, iforce_ref, lfgm, ncs,
                          icsinfo_ref, ihstry, ncycmax0, iv_ref, time2_ref, thetaxl_ref,
                          thetayl_ref, thetazl_ref, intmax0, nsub10, iindex_ref, lig_ref,
                          lbg_ref, ibpntsg_ref, iipntsg_ref, icallgl, iunit, nou, bou,
                          ibufdim0, nbuf0, mglevg_ref, nemgl_ref, ipl3dtmp_ref, ntr, bcfiles,
                          mxbcfil0, utrnsae_ref, vtrnsae_ref, wtrnsae_ref, omgxae_ref,
                          omgyae_ref, omgzae_ref, xorgae_ref, yorgae_ref, zorgae_ref,
                          thtxae_ref, thtyae_ref, thtzae_ref, rfrqtae_ref, rfrqrae_ref,
                          icsi_ref, icsf_ref, jcsi_ref, jcsf_ref, kcsi_ref, kcsf_ref,
                          freq_ref, gmass_ref, damp_ref, x0_ref, gf0_ref, nmds0, maxaes0,
                          aesrfdat_ref, perturb_ref, iskip_ref, jskip_ref, kskip_ref,
                          nsegdfrm_ref, idfrmseg_ref, iaesurf_ref, maxsegdg0,
                          xorgae0_ref, yorgae0_ref, zorgae0_ref, icouple_ref, iprnsurf);
    }
    // For nplot3d<0 / nprint<0 the plot3d/print surfaces are auto-detected;
    // global() returns the detected counts in ip3dsurf/iprnsurf. Size the
    // surface arrays (lnplts -> nplots) to hold them, else the real run aborts
    // "nplots too small for auto-surface detection" (lnplts stayed 1).
    if (ip3dsurf > lnplts) lnplts = ip3dsurf;
    if (iprnsurf > lnplts) lnplts = iprnsurf;


    // Read dynamic patch input data
    if (iunst > 0) {
        ioflag = 2;
        imode2 = 1;
        int nintr = 0;
        FortranArray1DRef<int> idimg_ref = idimg.ref();
        FortranArray1DRef<int> jdimg_ref = jdimg.ref();
        FortranArray1DRef<int> kdimg_ref = kdimg.ref();
        FortranArray1DRef<int> levelg_ref = levelg.ref();
        FortranArray1DRef<int> ncgg_ref = ncgg.ref();
        FortranArray1DRef<int> nblg_ref = nblg.ref();
        FortranArray2DRef<int> iindx_ref = iindx.ref();
        FortranArray1DRef<int> llimit_ref = llimit.ref();
        FortranArray1DRef<int> iitmax_ref = iitmax.ref();
        FortranArray1DRef<int> mmcxie_ref = mmcxie.ref();
        FortranArray1DRef<int> mmceta_ref = mmceta.ref();
        FortranArray1DRef<int> ncheck_ref = ncheck.ref();
        FortranArray1DRef<int> iifit_ref = iifit.ref();
        FortranArray1DRef<int> iic0_ref = iic0.ref();
        FortranArray1DRef<int> iiorph_ref = iiorph.ref();
        FortranArray1DRef<int> iitoss_ref = iitoss.ref();
        FortranArray1DRef<int> ifiner_ref = ifiner.ref();
        FortranArray2DRef<double> dx_ref = dx.ref();
        FortranArray2DRef<double> dy_ref = dy.ref();
        FortranArray2DRef<double> dz_ref = dz.ref();
        FortranArray2DRef<double> dthetx_ref = dthetx.ref();
        FortranArray2DRef<double> dthety_ref = dthety.ref();
        FortranArray2DRef<double> dthetz_ref = dthetz.ref();
        FortranArray1DRef<int> ifrom_ref = ifrom.ref();
        FortranArray1DRef<int> xif1_ref = xif1.ref();
        FortranArray1DRef<int> etf1_ref = etf1.ref();
        FortranArray1DRef<int> xif2_ref = xif2.ref();
        FortranArray1DRef<int> etf2_ref = etf2.ref();
        FortranArray1DRef<int> igridg_ref = igridg.ref();
        FortranArray1DRef<int> iemg_ref = iemg.ref();
        global2_ns::global2(maxbl0, maxgr0, msub10, nintr, intmx0, ngrid,
                            idimg_ref, jdimg_ref, kdimg_ref, levelg_ref,
                            ncgg_ref, nblg_ref, iindx_ref, llimit_ref,
                            iitmax_ref, mmcxie_ref, mmceta_ref, ncheck_ref,
                            iifit_ref, iic0_ref, iiorph_ref, iitoss_ref,
                            ifiner_ref, dx_ref, dy_ref, dz_ref, dthetx_ref,
                            dthety_ref, dthetz_ref, myid, mptch0, mxxe0,
                            icallgl, iunit, nou, bou, ibufdim0, nbuf0,
                            ifrom_ref, xif1_ref, etf1_ref, xif2_ref, etf2_ref,
                            igridg_ref, iemg_ref, nblock, ioflag, imode2);
    }

    // Allocate blocks to nodes
    {
        FortranArray1DRef<int> ncgg_ref = ncgg.ref();
        FortranArray1DRef<int> nblg_ref = nblg.ref();
        FortranArray1DRef<int> idimg_ref = idimg.ref();
        FortranArray1DRef<int> jdimg_ref = jdimg.ref();
        FortranArray1DRef<int> kdimg_ref = kdimg.ref();
        FortranArray1DRef<int> nblcg_ref = nblcg.ref();
        FortranArray2DRef<int> iwrk_ref = iwrk.ref();
        FortranArray1DRef<int> mblk2nd_ref = mblk2nd.ref();
        compg2n_ns::compg2n(nblock, ngrid, ncgg_ref, nblg_ref, idimg_ref,
                            jdimg_ref, kdimg_ref, nblcg_ref, nnodes, iwrk_ref,
                            myid, myhost, mblk2nd_ref, mycomm, maxgr0, maxbl0,
                            ierrflg, ibufdim0, nbuf0, bou, nou);
    }

    fprintf(f66, "\n BLOCK TO NODE MAPPING\n");
    fprintf(f66, " no. of blocks = %4d\n no. of  nodes = %4d\n block    node\n", nblock, nnodes);
    for (i = 1; i <= nblock; i++) {
        fprintf(f66, "%5d   %5d\n", i, mblk2nd(i));
    }

    // Check total number of nodes actually used
    nodel = 1;
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (mblk2nd(nbl) > nodel) nodel = mblk2nd(nbl);
    }
    lminnode = nodel;

    // Set up lw, lw2 pointer arrays
    icallpt = 0;
    for (nbl = 1; nbl <= nblock; nbl++) {
        memblock(nbl) = 0;
    }

    // Call pointers with myiduse=myhost and mpihost=0 for sequential host
    myiduse = myhost;
    mpihost = 0;
    {
        FortranArray2DRef<int> lw_ref = lw.ref();
        FortranArray2DRef<int> lw2_ref = lw2.ref();
        FortranArray3DRef<int> lwdat_ref = lwdat.ref();
        FortranArray1DRef<int> levelg_ref = levelg.ref();
        FortranArray1DRef<int> igridg_ref = igridg.ref();
        FortranArray2DRef<int> iviscg_ref = iviscg.ref();
        FortranArray1DRef<int> idimg_ref = idimg.ref();
        FortranArray1DRef<int> jdimg_ref = jdimg.ref();
        FortranArray1DRef<int> kdimg_ref = kdimg.ref();
        FortranArray1DRef<int> nblcg_ref = nblcg.ref();
        FortranArray1DRef<int> itrans_ref = itrans.ref();
        FortranArray1DRef<int> irotat_ref = irotat.ref();
        FortranArray1DRef<int> idefrm_ref = idefrm.ref();
        FortranArray1DRef<int> nbci0_ref = nbci0.ref();
        FortranArray1DRef<int> nbcj0_ref = nbcj0.ref();
        FortranArray1DRef<int> nbck0_ref = nbck0.ref();
        FortranArray1DRef<int> nbcidim_ref = nbcidim.ref();
        FortranArray1DRef<int> nbcjdim_ref = nbcjdim.ref();
        FortranArray1DRef<int> nbckdim_ref = nbckdim.ref();
        FortranArray4DRef<int> ibcinfo_ref = ibcinfo.ref();
        FortranArray4DRef<int> jbcinfo_ref = jbcinfo.ref();
        FortranArray4DRef<int> kbcinfo_ref = kbcinfo.ref();
        FortranArray1DRef<int> ncgg_ref = ncgg.ref();
        FortranArray1DRef<int> nblg_ref = nblg.ref();
        FortranArray1DRef<int> iemg_ref = iemg.ref();
        FortranArray1DRef<int> mblk2nd_ref = mblk2nd.ref();
        FortranArray1DRef<int> nblfine_ref = nblfine.ref();
        FortranArray1DRef<int> ilamlog_ref = ilamlog.ref();
        FortranArray1DRef<int> jlamlog_ref = jlamlog.ref();
        FortranArray1DRef<int> klamlog_ref = klamlog.ref();
        FortranArray1DRef<int> ilamhig_ref = ilamhig.ref();
        FortranArray1DRef<int> jlamhig_ref = jlamhig.ref();
        FortranArray1DRef<int> klamhig_ref = klamhig.ref();
        FortranArray2DRef<int> idegg_ref = idegg.ref();
        FortranArray2DRef<int> iwfg_ref = iwfg.ref();
        FortranArray2DRef<int> idiagg_ref = idiagg.ref();
        FortranArray2DRef<int> iflimg_ref = iflimg.ref();
        FortranArray2DRef<int> ifdsg_ref = ifdsg.ref();
        FortranArray2DRef<double> rkap0g_ref = rkap0g.ref();
        FortranArray1DRef<int> jsg_ref = jsg.ref();
        FortranArray1DRef<int> ksg_ref = ksg.ref();
        FortranArray1DRef<int> isg_ref = isg.ref();
        FortranArray1DRef<int> jeg_ref = jeg.ref();
        FortranArray1DRef<int> keg_ref = keg.ref();
        FortranArray1DRef<int> ieg_ref = ieg.ref();
        FortranArray1DRef<int> memblock_ref = memblock.ref();
        pointers_ns::pointers(lw_ref, lw2_ref, maxl, lembed, nstart_hseq, nwork, mwork,
                          maxbl0, maxgr0, maxseg0, lwdat_ref, levelg_ref, igridg_ref,
                          iviscg_ref, idimg_ref, jdimg_ref, kdimg_ref, nblcg_ref,
                          itrans_ref, irotat_ref, idefrm_ref, nbci0_ref, nbcj0_ref,
                          nbck0_ref, nbcidim_ref, nbcjdim_ref, nbckdim_ref,
                          ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ngrid, ncgg_ref,
                          nblg_ref, iemg_ref, nblock, myhost, myiduse, mblk2nd_ref,
                          nou, bou, nbuf0, ibufdim0, nblfine_ref, ilamlog_ref,
                          jlamlog_ref, klamlog_ref, ilamhig_ref, jlamhig_ref,
                          klamhig_ref, idegg_ref, iwfg_ref, idiagg_ref, iflimg_ref,
                          ifdsg_ref, rkap0g_ref, jsg_ref, ksg_ref, isg_ref, jeg_ref,
                          keg_ref, ieg_ref, memblock_ref, icallpt, nmds0, maxaes0, mpihost);
    }


    // Output contents of lw() for sequential host
    fprintf(f66, "  \n");
    fprintf(f66, "***********************************************\n");
    fprintf(f66, "\n");
    fprintf(f66, "    PERMANENT STORAGE REQUIREMENTS - W ARRAY\n");
    fprintf(f66, "                SEQUENTIAL BUILD\n");
    fprintf(f66, "\n");
    fprintf(f66, "***********************************************\n");
    fprintf(f66, "\n summary of starting locations for block information on single node\n");

    for (nbl = 1; nbl <= nblock; nbl++) {
        igrid = igridg(nbl);
        j = jdimg(nbl);
        k = kdimg(nbl);
        i = idimg(nbl);
        fprintf(f66, "\n igrid block  jdim  kdim  idim\n");
        fprintf(f66, "%5d%5d%5d%5d%5d\n", igrid, nbl, j, k, i);
        for (ii = 1; ii <= 65; ii++) {
            fprintf(f66, "   lw(%2d,%3d)=       %10d\n", ii, nbl, lw(ii, nbl));
        }
        for (jj = 1; jj <= 6; jj++) {
            if (jj == 1) nseg = nbci0(nbl);
            if (jj == 2) nseg = nbcidim(nbl);
            if (jj == 3) nseg = nbcj0(nbl);
            if (jj == 4) nseg = nbcjdim(nbl);
            if (jj == 5) nseg = nbck0(nbl);
            if (jj == 6) nseg = nbckdim(nbl);
            for (ii = 1; ii <= nseg; ii++) {
                fprintf(f66, "   lwdat(%3d,%2d,%2d)= %10d\n", nbl, ii, jj, lwdat(nbl, ii, jj));
            }
        }
    }

    // Call pointers with myiduse=myhost and mpihost=1 for parallel host
    myiduse = myhost;
    mpihost = 1;
    {
        FortranArray2DRef<int> lw_ref = lw.ref();
        FortranArray2DRef<int> lw2_ref = lw2.ref();
        FortranArray3DRef<int> lwdat_ref = lwdat.ref();
        FortranArray1DRef<int> levelg_ref = levelg.ref();
        FortranArray1DRef<int> igridg_ref = igridg.ref();
        FortranArray2DRef<int> iviscg_ref = iviscg.ref();
        FortranArray1DRef<int> idimg_ref = idimg.ref();
        FortranArray1DRef<int> jdimg_ref = jdimg.ref();
        FortranArray1DRef<int> kdimg_ref = kdimg.ref();
        FortranArray1DRef<int> nblcg_ref = nblcg.ref();
        FortranArray1DRef<int> itrans_ref = itrans.ref();
        FortranArray1DRef<int> irotat_ref = irotat.ref();
        FortranArray1DRef<int> idefrm_ref = idefrm.ref();
        FortranArray1DRef<int> nbci0_ref = nbci0.ref();
        FortranArray1DRef<int> nbcj0_ref = nbcj0.ref();
        FortranArray1DRef<int> nbck0_ref = nbck0.ref();
        FortranArray1DRef<int> nbcidim_ref = nbcidim.ref();
        FortranArray1DRef<int> nbcjdim_ref = nbcjdim.ref();
        FortranArray1DRef<int> nbckdim_ref = nbckdim.ref();
        FortranArray4DRef<int> ibcinfo_ref = ibcinfo.ref();
        FortranArray4DRef<int> jbcinfo_ref = jbcinfo.ref();
        FortranArray4DRef<int> kbcinfo_ref = kbcinfo.ref();
        FortranArray1DRef<int> ncgg_ref = ncgg.ref();
        FortranArray1DRef<int> nblg_ref = nblg.ref();
        FortranArray1DRef<int> iemg_ref = iemg.ref();
        FortranArray1DRef<int> mblk2nd_ref = mblk2nd.ref();
        FortranArray1DRef<int> nblfine_ref = nblfine.ref();
        FortranArray1DRef<int> ilamlog_ref = ilamlog.ref();
        FortranArray1DRef<int> jlamlog_ref = jlamlog.ref();
        FortranArray1DRef<int> klamlog_ref = klamlog.ref();
        FortranArray1DRef<int> ilamhig_ref = ilamhig.ref();
        FortranArray1DRef<int> jlamhig_ref = jlamhig.ref();
        FortranArray1DRef<int> klamhig_ref = klamhig.ref();
        FortranArray2DRef<int> idegg_ref = idegg.ref();
        FortranArray2DRef<int> iwfg_ref = iwfg.ref();
        FortranArray2DRef<int> idiagg_ref = idiagg.ref();
        FortranArray2DRef<int> iflimg_ref = iflimg.ref();
        FortranArray2DRef<int> ifdsg_ref = ifdsg.ref();
        FortranArray2DRef<double> rkap0g_ref = rkap0g.ref();
        FortranArray1DRef<int> jsg_ref = jsg.ref();
        FortranArray1DRef<int> ksg_ref = ksg.ref();
        FortranArray1DRef<int> isg_ref = isg.ref();
        FortranArray1DRef<int> jeg_ref = jeg.ref();
        FortranArray1DRef<int> keg_ref = keg.ref();
        FortranArray1DRef<int> ieg_ref = ieg.ref();
        FortranArray1DRef<int> memblock_ref = memblock.ref();
        pointers_ns::pointers(lw_ref, lw2_ref, maxl, lembed, nstart_hmpi, nwork, mwork,
                          maxbl0, maxgr0, maxseg0, lwdat_ref, levelg_ref, igridg_ref,
                          iviscg_ref, idimg_ref, jdimg_ref, kdimg_ref, nblcg_ref,
                          itrans_ref, irotat_ref, idefrm_ref, nbci0_ref, nbcj0_ref,
                          nbck0_ref, nbcidim_ref, nbcjdim_ref, nbckdim_ref,
                          ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ngrid, ncgg_ref,
                          nblg_ref, iemg_ref, nblock, myhost, myiduse, mblk2nd_ref,
                          nou, bou, nbuf0, ibufdim0, nblfine_ref, ilamlog_ref,
                          jlamlog_ref, klamlog_ref, ilamhig_ref, jlamhig_ref,
                          klamhig_ref, idegg_ref, iwfg_ref, idiagg_ref, iflimg_ref,
                          ifdsg_ref, rkap0g_ref, jsg_ref, ksg_ref, isg_ref, jeg_ref,
                          keg_ref, ieg_ref, memblock_ref, icallpt, nmds0, maxaes0, mpihost);
    }

    // Loop over parallel nodes
    for (myiduse = 1; myiduse <= nnodes; myiduse++) {
        FortranArray2DRef<int> lw_ref = lw.ref();
        FortranArray2DRef<int> lw2_ref = lw2.ref();
        FortranArray3DRef<int> lwdat_ref = lwdat.ref();
        FortranArray1DRef<int> levelg_ref = levelg.ref();
        FortranArray1DRef<int> igridg_ref = igridg.ref();
        FortranArray2DRef<int> iviscg_ref = iviscg.ref();
        FortranArray1DRef<int> idimg_ref = idimg.ref();
        FortranArray1DRef<int> jdimg_ref = jdimg.ref();
        FortranArray1DRef<int> kdimg_ref = kdimg.ref();
        FortranArray1DRef<int> nblcg_ref = nblcg.ref();
        FortranArray1DRef<int> itrans_ref = itrans.ref();
        FortranArray1DRef<int> irotat_ref = irotat.ref();
        FortranArray1DRef<int> idefrm_ref = idefrm.ref();
        FortranArray1DRef<int> nbci0_ref = nbci0.ref();
        FortranArray1DRef<int> nbcj0_ref = nbcj0.ref();
        FortranArray1DRef<int> nbck0_ref = nbck0.ref();
        FortranArray1DRef<int> nbcidim_ref = nbcidim.ref();
        FortranArray1DRef<int> nbcjdim_ref = nbcjdim.ref();
        FortranArray1DRef<int> nbckdim_ref = nbckdim.ref();
        FortranArray4DRef<int> ibcinfo_ref = ibcinfo.ref();
        FortranArray4DRef<int> jbcinfo_ref = jbcinfo.ref();
        FortranArray4DRef<int> kbcinfo_ref = kbcinfo.ref();
        FortranArray1DRef<int> ncgg_ref = ncgg.ref();
        FortranArray1DRef<int> nblg_ref = nblg.ref();
        FortranArray1DRef<int> iemg_ref = iemg.ref();
        FortranArray1DRef<int> mblk2nd_ref = mblk2nd.ref();
        FortranArray1DRef<int> nblfine_ref = nblfine.ref();
        FortranArray1DRef<int> ilamlog_ref = ilamlog.ref();
        FortranArray1DRef<int> jlamlog_ref = jlamlog.ref();
        FortranArray1DRef<int> klamlog_ref = klamlog.ref();
        FortranArray1DRef<int> ilamhig_ref = ilamhig.ref();
        FortranArray1DRef<int> jlamhig_ref = jlamhig.ref();
        FortranArray1DRef<int> klamhig_ref = klamhig.ref();
        FortranArray2DRef<int> idegg_ref = idegg.ref();
        FortranArray2DRef<int> iwfg_ref = iwfg.ref();
        FortranArray2DRef<int> idiagg_ref = idiagg.ref();
        FortranArray2DRef<int> iflimg_ref = iflimg.ref();
        FortranArray2DRef<int> ifdsg_ref = ifdsg.ref();
        FortranArray2DRef<double> rkap0g_ref = rkap0g.ref();
        FortranArray1DRef<int> jsg_ref = jsg.ref();
        FortranArray1DRef<int> ksg_ref = ksg.ref();
        FortranArray1DRef<int> isg_ref = isg.ref();
        FortranArray1DRef<int> jeg_ref = jeg.ref();
        FortranArray1DRef<int> keg_ref = keg.ref();
        FortranArray1DRef<int> ieg_ref = ieg.ref();
        FortranArray1DRef<int> memblock_ref = memblock.ref();
        pointers_ns::pointers(lw_ref, lw2_ref, maxl, lembed, nstart, nwork, mwork,
                          maxbl0, maxgr0, maxseg0, lwdat_ref, levelg_ref, igridg_ref,
                          iviscg_ref, idimg_ref, jdimg_ref, kdimg_ref, nblcg_ref,
                          itrans_ref, irotat_ref, idefrm_ref, nbci0_ref, nbcj0_ref,
                          nbck0_ref, nbcidim_ref, nbcjdim_ref, nbckdim_ref,
                          ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ngrid, ncgg_ref,
                          nblg_ref, iemg_ref, nblock, myhost, myiduse, mblk2nd_ref,
                          nou, bou, nbuf0, ibufdim0, nblfine_ref, ilamlog_ref,
                          jlamlog_ref, klamlog_ref, ilamhig_ref, jlamhig_ref,
                          klamhig_ref, idegg_ref, iwfg_ref, idiagg_ref, iflimg_ref,
                          ifdsg_ref, rkap0g_ref, jsg_ref, ksg_ref, isg_ref, jeg_ref,
                          keg_ref, ieg_ref, memblock_ref, icallpt, nmds0, maxaes0, mpihost);
    }


    // Determine permanent w array requirement for nodes
    for (ii = 1; ii <= nodel; ii++) {
        mem_req_node(ii) = 0;
        no_of_points(ii) = 0;
        for (nbl = 1; nbl <= nblock; nbl++) {
            nfiner = 0;
            for (mbl = 1; mbl <= nblock; mbl++) {
                if (levelg(mbl) > lglobal) {
                    if (nblcg(mbl) == nbl) {
                        nfiner = nfiner + 1;
                        nblfine(nfiner) = mbl;
                    }
                }
            }
            nblc = nblcg(nbl);
            iaug = 0;
            if (ii == mblk2nd(nbl)) iaug = 1;
            if (ii == mblk2nd(nblc)) iaug = 1;
            if (nfiner > 0) {
                for (nf = 1; nf <= nfiner; nf++) {
                    if (ii == mblk2nd(nblfine(nf))) iaug = 1;
                }
            }
            if (iaug > 0) {
                npoints = (idimg(nbl)-1)*(jdimg(nbl)-1)*(kdimg(nbl)-1);
                mem_req_node(ii) = mem_req_node(ii) + memblock(nbl);
                no_of_points(ii) = no_of_points(ii) + npoints;
            }
        }
    }

    npts_max  = 0;
    mem_w_max = nstart_hmpi;
    for (i = 1; i <= nodel; i++) {
        if (mem_req_node(i) > mem_w_max) mem_w_max = mem_req_node(i);
        if (no_of_points(i) > npts_max)  npts_max  = no_of_points(i);
    }

    fprintf(f66, "\n");
    fprintf(f66, "***********************************************\n");
    fprintf(f66, "\n");
    fprintf(f66, "   PERMANENT STORAGE REQUIREMENTS - W ARRAY    \n");
    fprintf(f66, "            SUMMARY FOR ALL BUILDS\n");
    fprintf(f66, "\n");
    fprintf(f66, "***********************************************\n");

    ftot_seq = (double)nstart_hseq / 1.e6;
    fprintf(f66, "\n memory (mw) for w storage (sequential) = %9.4f\n", (float)ftot_seq);

    fprintf(f66, "\n memory (mw) for w storage (nodes)\n");
    fprintf(f66, "\n node   total points  memory (mw)\n");
    fmemsum = nstart_hmpi / 1.e6;
    nptsum  = 0;
    fprintf(f66, "%6d       %8d    %9.4f\n", myhost, npts_max, (float)fmemsum);
    for (i = 1; i <= nodel; i++) {
        fprintf(f66, "%6d       %8d    %9.4f\n", i, no_of_points(i), (float)(mem_req_node(i)/1.e6));
        fmemsum = fmemsum + mem_req_node(i)/1.e6;
        nptsum  = nptsum + no_of_points(i);
        mem_req_node(i) = 0;
    }
    fprintf(f66, "\n     total =  %8d    %9.4f\n", nptsum, (float)fmemsum);
    fprintf(f66, "\n   maximum =  %8d    %9.4f\n", npts_max, (float)((float)mem_w_max/1.e6f));
    fprintf(f66, "\n");

    ftot_seq = (double)nstart_hseq;

    // Evaluate requirements for temporary real array size (work/wk array)
    fprintf(f66, "\n");
    fprintf(f66, "***********************************************\n");
    fprintf(f66, "\n");
    fprintf(f66, "   TEMPORARY STORAGE REQUIREMENTS - WK ARRAY\n");
    fprintf(f66, "        SEQUENTIAL AND PARALLEL BUILDS\n");
    fprintf(f66, "\n");
    fprintf(f66, "***********************************************\n");

    for (n = 1; n <= maxbl0; n++) {
        ltot(n) = 0;
    }

    need  = lmaxbl;
    needi = 1;
    needi_node = 1;


    // Requirements of subroutine addx
    // Outer loop over iseq (mesh sequencing)
    for (iseq = 1; iseq <= mseq; iseq++) {
        if (ncyc1[iseq-1] <= 0) goto label_8500;
        if (iseq > 1 || iseq == mseq) {
            if (iseq <= mseq && iseq != 1) {
                for (igrid = 1; igrid <= ngrid; igrid++) {
                    int iem = iemg(igrid);
                    if (iem > 0) continue;
                    nblf  = nblg(igrid) + (mseq - iseq);
                    nblz  = nblf + 1;
                    {
                        FortranArray2DRef<int> lw_ref = lw.ref();
                        FortranArray2DRef<int> lw2_ref = lw2.ref();
                        lead_ns::lead(nblf, lw_ref, lw2_ref, maxbl0);
                    }
                    iwk1  = 1;
                    iwk2  = iwk1 + jdim*kdim*idim*7;
                    iwk3  = iwk2 + jj2*kk2*ii2*7;
                    itemp = iwk3 + jdim*kk2*ii2*7;
                    need  = std::max(need, itemp);
                    ii = mblk2nd(nblf);
                    mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
                    fprintf(f66, " 1(addx)         itemp, need = %d %d\n", itemp, need);
                    fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
                }
            }
            // Requirements of subroutine add2x
            if (lembed == 0 || iseq != mseq) goto label_333;
            for (level = mseq+1; level <= maxl; level++) {
                for (igrid = 1; igrid <= ngrid; igrid++) {
                    nbl   = nblg(igrid);
                    int iem   = iemg(igrid);
                    int inewg = inewgg(igrid);
                    if (iem == 0) continue;
                    if (iseq == 1 && inewg == 0) continue;
                    if (levelg(nbl) != level) continue;
                    {
                        FortranArray2DRef<int> lw_ref = lw.ref();
                        FortranArray2DRef<int> lw2_ref = lw2.ref();
                        lead_ns::lead(nbl, lw_ref, lw2_ref, maxbl0);
                    }
                    jj2   = jdimg(nblc);
                    kk2   = kdimg(nblc);
                    ii2   = idimg(nblc);
                    iwk1  = 1;
                    iwk2  = iwk1 + jdim*kdim*idim*7;
                    iwk3  = iwk2 + jj2*kk2*ii2*7;
                    iwk4  = iwk3 + jdim*kk2*ii2;
                    itemp = iwk4 + jdim*kdim*ii2;
                    need  = std::max(need, itemp);
                    ii = mblk2nd(nbl);
                    mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
                    fprintf(f66, " 2(add2x)        itemp, need = %d %d\n", itemp, need);
                    fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
                }
            }
        }
        label_333:
        levt  = levelt[iseq-1];
        levb  = levelb[iseq-1];
        level = levelt[iseq-1];
        lglobal = lfgm - (mseq - iseq);
        label_7000:
        for (i = 1; i <= level-levb+1; i++) {
            nsm[i-1] = ngam;
        }
        nit    = mit(level-levb+1, iseq);
        kxpand = 1;
        label_9000:
        ilc    = 0;

        for (nbl = 1; nbl <= nblock; nbl++) {
            if (level != levelg(nbl)) continue;
            {
                FortranArray2DRef<int> lw_ref = lw.ref();
                FortranArray2DRef<int> lw2_ref = lw2.ref();
                lead_ns::lead(nbl, lw_ref, lw2_ref, maxbl0);
            }
            nit1 = nit + std::min(1, mgflag);
            if (level == levb) nit1 = nit;
            if (level > lglobal) nit1 = nit + 1;
            if (level == lglobal) nit1 = nit;
            if (nit1 > 0) {
                for (ntime = 1; ntime <= nit1; ntime++) {
                    if (ilc == 0) {
                        nsm[level-levb] = nsm[level-levb] - 1;
                        if (nsm[level-levb] < 0) nsm[level-levb] = ngam - 1;
                        ilc = 1;
                    }
                    igrid = igridg(nbl);
                    int iem   = iemg(igrid);
                    int ifluxa = 0;
                    if (level >= lglobal && level != levt) ifluxa = std::min(1, iconsf);
                    if (ifluxa > 0) {
                        if (ntime == 1) {
                            FortranArray1DRef<int> nblcg_ref = nblcg.ref();
                            FortranArray1DRef<int> ieg_ref = ieg.ref();
                            FortranArray1DRef<int> isg_ref = isg.ref();
                            FortranArray1DRef<int> jdimg_ref = jdimg.ref();
                            FortranArray1DRef<int> kdimg_ref = kdimg.ref();
                            FortranArray1DRef<int> idimg_ref = idimg.ref();
                            FortranArray1DRef<int> nbcj0_ref = nbcj0.ref();
                            FortranArray1DRef<int> nbck0_ref = nbck0.ref();
                            FortranArray1DRef<int> nbci0_ref = nbci0.ref();
                            FortranArray1DRef<int> nbcjdim_ref = nbcjdim.ref();
                            FortranArray1DRef<int> nbckdim_ref = nbckdim.ref();
                            FortranArray1DRef<int> nbcidim_ref = nbcidim.ref();
                            FortranArray4DRef<int> jbcinfo_ref = jbcinfo.ref();
                            FortranArray4DRef<int> kbcinfo_ref = kbcinfo.ref();
                            FortranArray4DRef<int> ibcinfo_ref = ibcinfo.ref();
                            cntfa(nbl, nwfa, maxbl0, maxseg0, nblcg_ref, ieg_ref, isg_ref,
                                  jdimg_ref, kdimg_ref, idimg_ref, nbcj0_ref, nbck0_ref,
                                  nbci0_ref, nbcjdim_ref, nbckdim_ref, nbcidim_ref,
                                  jbcinfo_ref, kbcinfo_ref, ibcinfo_ref, nblock);
                        }
                    } else {
                        nwfa = 0;
                    }


                    int lres      = nwfa  + 1;
                    int lwj0      = lres  + jdim*kdim*(idim-1)*5;
                    int lwk0      = lwj0  + kdim*idim*22;
                    int lwi0      = lwk0  + jdim*idim*22;
                    int lvmuk     = lwi0  + kdim*jdim*22;
                    int lvmuj     = lvmuk + 2*(jdim-1)*(idim-1);
                    int lvmui     = lvmuj + 2*(kdim-1)*(idim-1);
                    ltot(nbl) = lvmui + 2*(kdim-1)*(jdim-1) + maxbl0;
                    itemp = 0;

                    if (ivmx >= 2) {
                        if (ivmx == 2) {
                            // Baldwin-Lomax
                            inmx = std::max({jdim, kdim, idim});
                            iwk1 = 1; iwk2 = iwk1+inmx; iwk3 = iwk2+inmx;
                            iwk4 = iwk3+inmx; iwk5 = iwk4+inmx; iwk6 = iwk5+inmx;
                            iwk7 = iwk6+inmx; iwk8 = iwk7+inmx; iwk9 = iwk8+inmx;
                            itemp= iwk9+inmx;
                        } else if (ivmx == 4) {
                            // Baldwin-Barth
                            iwk1 = 1;
                            if (iturbord == 1) iwk5 = iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                            else               iwk5 = iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d));
                            inmx=(kdim-1)*(jdim-1); iwk6=iwk5+inmx; iwk7=iwk6+inmx; iwk8=iwk7+inmx; iwk9=iwk8+inmx; iwk10=iwk9+inmx;
                            inmx=(jdim-1)*(kdim-1); iwk11=iwk10+inmx; iwk12=iwk11+inmx; iwk13=iwk12+inmx; iwk14=iwk13+inmx; iwk15=iwk14+inmx;
                            inmx=(kdim-1)*(idim-1); iwk16=iwk15+inmx; iwk17=iwk16+inmx; iwk18=iwk17+inmx; iwk19=iwk18+inmx; itemp=iwk19+inmx;
                        } else if (ivmx == 5) {
                            // Spalart
                            iwk1 = 1;
                            if (iturbord == 1) iwk3 = iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                            else               iwk3 = iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d));
                            iwk4 = iwk3+(jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                            inmx=(kdim-1)*(jdim-1); iwk5=iwk4+inmx; iwk6=iwk5+inmx; iwk7=iwk6+inmx; iwk8=iwk7+inmx; iwk9=iwk8+inmx;
                            inmx=(jdim-1)*(kdim-1); iwk10=iwk9+inmx; iwk11=iwk10+inmx; iwk12=iwk11+inmx; iwk13=iwk12+inmx; iwk14=iwk13+inmx;
                            inmx=(kdim-1)*(idim-1); iwk15=iwk14+inmx; iwk16=iwk15+inmx; iwk17=iwk16+inmx; iwk18=iwk17+inmx; iwk19=iwk18+inmx;
                            if (isarc2d == 1)      inmx=(jdim+1)*(kdim+1)*(idim-1)*4;
                            else if (isarc3d == 1) inmx=(jdim+1)*(kdim+1)*(idim-1)*6;
                            else                   inmx=0;
                            itemp=iwk19+inmx;
                        } else if (ivmx >= 6 && ivmx <= 16) {
                            // Two-equation turb models
                            iwk1 = 1;
                            if (iturbord == 1) iwk4=iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))*2;
                            else               iwk4=iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d))*2;
                            iwk5=iwk4+(jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                            inmx=(kdim-1)*(jdim-1);
                            iwk6=iwk5+inmx; iwk7=iwk6+inmx; iwk8=iwk7+inmx; iwk9=iwk8+inmx; iwk10=iwk9+inmx;
                            iwk11=iwk10+inmx; iwk12=iwk11+inmx; iwk13=iwk12+inmx; iwk14=iwk13+inmx;
                            inmx=(jdim-1)*(kdim-1);
                            iwk15=iwk14+inmx; iwk16=iwk15+inmx; iwk17=iwk16+inmx; iwk18=iwk17+inmx; iwk19=iwk18+inmx;
                            iwk20=iwk19+inmx; iwk21=iwk20+inmx; iwk22=iwk21+inmx; iwk23=iwk22+inmx;
                            inmx=(kdim-1)*(idim-1);
                            iwk24=iwk23+inmx; iwk25=iwk24+inmx; iwk26=iwk25+inmx; iwk27=iwk26+inmx; iwk28=iwk27+inmx;
                            iwk29=iwk28+inmx; iwk30=iwk29+inmx; iwk31=iwk30+inmx; iwk32=iwk31+inmx;
                            inmx=(jdim-1)*(kdim-1)*(idim-1)*2;
                            iwk33=iwk32+inmx;
                            inmx=(jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                            // allocate for gradients of vorticity and turbulence - k-zeta only:
                            iwk34=iwk33+inmx;
                            if (ivmx == 15) inmx=jdim; else inmx=0;
                            iwk35=iwk34+inmx; iwk36=iwk35+inmx; iwk37=iwk36+inmx; iwk38=iwk37+inmx;
                            iwk39=iwk38+inmx; iwk40=iwk39+inmx; iwk41=iwk40+inmx; iwk42=iwk41+inmx;
                            iwk43=iwk42+inmx; iwk44=iwk43+inmx; iwk45=iwk44+inmx; iwk46=iwk45+inmx;
                            iwk47=iwk46+inmx; iwk48=iwk47+inmx; iwk49=iwk48+inmx; iwk50=iwk49+inmx;
                            iwk51=iwk50+inmx; iwk52=iwk51+inmx; iwk53=iwk52+inmx; iwk54=iwk53+inmx;
                            iwk55=iwk54+inmx; iwk56=iwk55+inmx; iwk57=iwk56+inmx; iwk58=iwk57+inmx;
                            iwk59=iwk58+inmx; iwk60=iwk59+inmx; iwk61=iwk60+inmx; iwk62=iwk61+inmx;
                            iwk63=iwk62+inmx; iwk64=iwk63+inmx; iwk65=iwk64+inmx; iwk66=iwk65+inmx;
                            iwk67=iwk66+inmx; iwk68=iwk67+inmx; iwk69=iwk68+inmx; iwk70=iwk69+inmx;
                            iwk71=iwk70+inmx; iwk72=iwk71+inmx; iwk73=iwk72+inmx; iwk74=iwk73+inmx;
                            iwk75=iwk74+inmx; iwk76=iwk75+inmx;
                            if (ivmx == 15) inmx=(jdim+1)*(kdim+1)*(idim+1); else inmx=0;
                            iwk77=iwk76+inmx; iwk78=iwk77+inmx; iwk79=iwk78+inmx;
                            if (ivmx == 15) inmx=(jdim+1)*(kdim+1)*(idim+1)*2; else inmx=0;
                            iwk80=iwk79+inmx;
                            // allocate for easmcc2d
                            if (ieasmcc2d == 1)      inmx=(jdim+1)*(kdim+1)*(idim-1)*4;
                            else if (isstrc == 2)     inmx=(jdim+1)*(kdim+1)*(idim-1)*6;
                            else if (ivmx == 16)      inmx=(jdim+1)*(kdim+1)*(idim-1)*4;
                            else                      inmx=0;
                            iwk81=iwk80+inmx;
                            // allocate for DES and DDES
                            if (ides != 0) inmx=(jdim-1)*(kdim-1)*(idim-1); else inmx=0;
                            iwk82=iwk81+inmx;
                            if (ides == 3) inmx=(jdim-1)*(kdim-1)*(idim-1); else inmx=0;
                            itemp=iwk82+inmx;
                        } else if (ivmx == 30) {
                            // three-equation turb model
                            iwk1 = 1;
                            if (iturbord == 1) { iwk4=iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))*3; iex=0; iex2=-i2d; }
                            else               { iwk4=iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d))*3; iex=1; iex2=1-(2*i2d); }
                            iwk5=iwk4+(jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                            iex3=-i2d;
                            inmx=(kdim-1)*(jdim-1);
                            iwk6=iwk5+inmx; iwk7=iwk6+inmx; iwk8=iwk7+inmx; iwk9=iwk8+inmx; iwk10=iwk9+inmx;
                            iwk11=iwk10+inmx; iwk12=iwk11+inmx; iwk13=iwk12+inmx; iwk14=iwk13+inmx;
                            inmx=(jdim-1)*(kdim-1);
                            iwk15=iwk14+inmx; iwk16=iwk15+inmx; iwk17=iwk16+inmx; iwk18=iwk17+inmx; iwk19=iwk18+inmx;
                            iwk20=iwk19+inmx; iwk21=iwk20+inmx; iwk22=iwk21+inmx; iwk23=iwk22+inmx;
                            inmx=(kdim-1)*(idim-1);
                            iwk24=iwk23+inmx; iwk25=iwk24+inmx; iwk26=iwk25+inmx; iwk27=iwk26+inmx; iwk28=iwk27+inmx;
                            iwk29=iwk28+inmx; iwk30=iwk29+inmx; iwk31=iwk30+inmx; iwk32=iwk31+inmx;
                            inmx=(jdim-1)*(kdim-1)*(idim-1)*3;
                            iwk33=iwk32+inmx;
                            inmx=(jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                            iwk34=iwk33+inmx;
                            inmx=(kdim-1)*(jdim-1);
                            iwk35=iwk34+inmx; iwk36=iwk35+inmx; iwk37=iwk36+inmx; iwk38=iwk37+inmx;
                            inmx=(jdim-1)*(kdim-1);
                            iwk39=iwk38+inmx; iwk40=iwk39+inmx; iwk41=iwk40+inmx; iwk42=iwk41+inmx;
                            inmx=(kdim-1)*(idim-1);
                            iwk43=iwk42+inmx; iwk44=iwk43+inmx; iwk45=iwk44+inmx; itemp=iwk45+inmx;


                        } else if (ivmx == 40) {
                            // four-equation turb model
                            iwk1 = 1;
                            if (iturbord == 1) { iwk4=iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))*4; iex=0; iex2=-i2d; }
                            else               { iwk4=iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d))*4; iex=1; iex2=1-(2*i2d); }
                            iwk5=iwk4+(jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                            iex3=-i2d;
                            inmx=(kdim-1)*(jdim-1);
                            iwk6=iwk5+inmx; iwk7=iwk6+inmx; iwk8=iwk7+inmx; iwk9=iwk8+inmx; iwk10=iwk9+inmx;
                            iwk11=iwk10+inmx; iwk12=iwk11+inmx; iwk13=iwk12+inmx; iwk14=iwk13+inmx;
                            inmx=(jdim-1)*(kdim-1);
                            iwk15=iwk14+inmx; iwk16=iwk15+inmx; iwk17=iwk16+inmx; iwk18=iwk17+inmx; iwk19=iwk18+inmx;
                            iwk20=iwk19+inmx; iwk21=iwk20+inmx; iwk22=iwk21+inmx; iwk23=iwk22+inmx;
                            inmx=(kdim-1)*(idim-1);
                            iwk24=iwk23+inmx; iwk25=iwk24+inmx; iwk26=iwk25+inmx; iwk27=iwk26+inmx; iwk28=iwk27+inmx;
                            iwk29=iwk28+inmx; iwk30=iwk29+inmx; iwk31=iwk30+inmx; iwk32=iwk31+inmx;
                            inmx=(jdim-1)*(kdim-1)*(idim-1)*4;
                            iwk33=iwk32+inmx;
                            inmx=(jdim+1)*(kdim+1)*(idim+1-(2*i2d));
                            iwk34=iwk33+inmx;
                            inmx=(kdim-1)*(jdim-1);
                            iwk35=iwk34+inmx; iwk36=iwk35+inmx; iwk37=iwk36+inmx; iwk38=iwk37+inmx;
                            inmx=(jdim-1)*(kdim-1);
                            iwk39=iwk38+inmx; iwk40=iwk39+inmx; iwk41=iwk40+inmx; iwk42=iwk41+inmx;
                            inmx=(kdim-1)*(idim-1);
                            iwk43=iwk42+inmx; iwk44=iwk43+inmx; iwk45=iwk44+inmx; iwk46=iwk45+inmx;
                            inmx=(kdim-1)*(jdim-1);
                            iwk47=iwk46+inmx; iwk48=iwk47+inmx; iwk49=iwk48+inmx; iwk50=iwk49+inmx;
                            inmx=(jdim-1)*(kdim-1);
                            iwk51=iwk50+inmx; iwk52=iwk51+inmx; iwk53=iwk52+inmx; iwk54=iwk53+inmx;
                            inmx=(kdim-1)*(idim-1);
                            iwk55=iwk54+inmx; iwk56=iwk55+inmx; iwk57=iwk56+inmx; iwk58=iwk57+inmx;
                            // allocate for DES and DDES
                            if (ides != 0) inmx=(jdim-1)*(kdim-1)*(idim-1); else inmx=0;
                            iwk59=iwk58+inmx;
                            if (ides == 3) inmx=(jdim-1)*(kdim-1)*(idim-1); else inmx=0;
                            itemp=iwk59+inmx;
                        } else if (ivmx == 72) {
                            itemp = 0;
                        }
                    } // end if ivmx >= 2

                    need = std::max(need, itemp);
                    ii = mblk2nd(nbl);
                    mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
                    fprintf(f66, " 3(turbs)        itemp, need = %d %d\n", itemp, need);
                    fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
                    itempt = itemp;

                    nroomf = 0;
                    nroomi = 0;

                    if (ivisc[0] >= 4 || ivisc[1] >= 4 || ivisc[2] >= 4) {
                        FortranArray1DRef<int> nblg_ref = nblg.ref();
                        FortranArray1DRef<int> nbci0_ref = nbci0.ref();
                        FortranArray1DRef<int> nbcj0_ref = nbcj0.ref();
                        FortranArray1DRef<int> nbck0_ref = nbck0.ref();
                        FortranArray1DRef<int> nbcidim_ref = nbcidim.ref();
                        FortranArray1DRef<int> nbcjdim_ref = nbcjdim.ref();
                        FortranArray1DRef<int> nbckdim_ref = nbckdim.ref();
                        FortranArray4DRef<int> ibcinfo_ref = ibcinfo.ref();
                        FortranArray4DRef<int> jbcinfo_ref = jbcinfo.ref();
                        FortranArray4DRef<int> kbcinfo_ref = kbcinfo.ref();
                        int ibctyp2004 = 2004, ibctyp2014 = 2014, ibctyp2024 = 2024;
                        int ibctyp2034 = 2034, ibctyp2016 = 2016;
                        cntsurf_ns::cntsurf(ns2004, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref,
                                            nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref,
                                            nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp2004);
                        cntsurf_ns::cntsurf(ns2014, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref,
                                            nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref,
                                            nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp2014);
                        cntsurf_ns::cntsurf(ns2024, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref,
                                            nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref,
                                            nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp2024);
                        cntsurf_ns::cntsurf(ns2034, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref,
                                            nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref,
                                            nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp2034);
                        cntsurf_ns::cntsurf(ns2016, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref,
                                            nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref,
                                            nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp2016);
                        nsurf = ns2004+ns2014+ns2024+ns2034+ns2016;
                        minbox = (int)std::sqrt((float)nsurf);
                        minbox = std::max(minbox, 50);
                        nbb = 3*nsurf/minbox;
                        for (int nbll2 = 1; nbll2 <= ngrid; nbll2++) {
                            int nbll = nblg(nbll2);
                            int ntempf = 9*nsurf + 7*nbb + 4*jdimg(nbll)*kdimg(nbll)*idimg(nbll);
                            int ntempi;
                            if (ivmx == 4 || ivmx == 25) ntempi = 15*nsurf + 2*nbb + jdimg(nbll)*kdimg(nbll)*idimg(nbll);
                            else                          ntempi = 11*nsurf + 2*nbb + jdimg(nbll)*kdimg(nbll)*idimg(nbll);
                            nroomf = std::max(nroomf, ntempf);
                            nroomi = std::max(nroomi, ntempi);
                        }
                        nadd_node = 0;
                        for (int nbll2 = 1; nbll2 <= nblock; nbll2++) {
                            if (mblk2nd(nbll2) == mblk2nd(nbl)) nadd_node += jdimg(nbll2)*kdimg(nbll2)*idimg(nbll2);
                        }
                        nadd = 0;
                        for (int nbll2 = 1; nbll2 <= nblock; nbll2++) {
                            nadd += jdimg(nbll2)*kdimg(nbll2)*idimg(nbll2);
                        }
                        nroomf0 = nroomf;
                        nroomf = nroomf0 + nadd;
                        nroomf_node = nroomf0 + nadd_node;
                        nworkxs = 4*nsurf;
                        if (ivmx == 4 || ivmx == 25) nworkixs = 4*nsurf; else nworkixs = 0;
                        nwklsmin = maxbl0;
                        nwkireq1 = maxbl0*maxseg0*6;
                        nwkireq2 = maxbl0*maxseg0*6;
                        nroom      = nroomf + nworkxs;
                        nroom_node = nroomf_node + nworkxs;
                        nroomi     = nroomi + nworkixs + nwklsmin + nwkireq1 + nwkireq2;
                        nroomi_node= nroomi + nworkixs + nwklsmin + nwkireq1 + nwkireq2;
                        itemp      = nroom;
                        itemp_node = nroom_node;
                        need       = std::max(need, itemp);
                        needi      = std::max(needi, nroomi);
                        needi_node = std::max(needi_node, nroomi_node);
                        ii = mblk2nd(nbl);
                        mem_req_node(ii) = std::max(mem_req_node(ii), itemp_node);
                        fprintf(f66, " 4(findmin_new)  itemp, need = %d %d\n", itemp, need);
                        fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp_node, mem_req_node(ii));
                    }


                    nv = 35;
                    maxmem = std::max(jdim*kdim*nv, jdim*idim*nv);
                    if (idefrm(nbl) > 0) {
                        nv = 41;
                        maxmem = std::max(jdim*kdim*nv, jdim*idim*nv);
                        mem_updatedg = 16*idim*jdim*kdim + 30*jdim*kdim + 24*jdim*idim + 24*kdim*idim;
                        maxmem = maxmem + mem_updatedg;
                    }
                    itemp  = maxmem + ltot(nbl);
                    itempt = itempt + ltot(nbl);
                    itemp  = std::max(itemp, itempt);
                    need   = std::max(need, itemp);
                    ii = mblk2nd(nbl);
                    mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
                    fprintf(f66, " 6(resid)        itemp, need = %d %d\n", itemp, need);
                    fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));

                    if (ntime == nit1 && level != levb) {
                        if (level <= lglobal) {
                            itemp = lwj0 + jdim*kdim*idim*5;
                            need  = std::max(need, itemp);
                            ii = mblk2nd(nbl);
                            mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
                            fprintf(f66, " 7(collq)        itemp, need = %d %d\n", itemp, need);
                            fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
                        } else {
                            itemp = lwj0 + jdim*kdim*idim*5;
                            need  = std::max(need, itemp);
                            ii = mblk2nd(nbl);
                            mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
                            fprintf(f66, " 8(coll2q)       itemp, need = %d %d\n", itemp, need);
                            fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
                        }
                    }
                } // end ntime loop
            } // end if nit1 > 0
        } // end nbl loop

        if (level == lglobal) goto label_671;
        label_2222:
        if (kxpand == -1) goto label_2012;
        if (level == levb) goto label_7000;

        level  = level - 1;
        ntime  = 0;
        kxpand = 1;
        goto label_9000;

        label_2012:
        if (level == levt) goto label_8500;
        level  = level + 1;
        kxpand = -1;
        goto label_9000;

        label_671:
        // Requirements of subroutine blomax
        if (ivmx == 2) {
            for (nbl = 1; nbl <= nblock; nbl++) {
                if (level != levelg(nbl)) continue;
                {
                    FortranArray2DRef<int> lw_ref = lw.ref();
                    FortranArray2DRef<int> lw2_ref = lw2.ref();
                    lead_ns::lead(nbl, lw_ref, lw2_ref, maxbl0);
                }
                FortranArray1DRef<int> nbcj0_ref = nbcj0.ref();
                FortranArray1DRef<int> nbck0_ref = nbck0.ref();
                FortranArray1DRef<int> nbci0_ref = nbci0.ref();
                FortranArray1DRef<int> nbcjdim_ref = nbcjdim.ref();
                FortranArray1DRef<int> nbckdim_ref = nbckdim.ref();
                FortranArray1DRef<int> nbcidim_ref = nbcidim.ref();
                FortranArray4DRef<int> jbcinfo_ref = jbcinfo.ref();
                FortranArray4DRef<int> kbcinfo_ref = kbcinfo.ref();
                FortranArray4DRef<int> ibcinfo_ref = ibcinfo.ref();
                cntblmx(nbl, jdim, kdim, idim, iovrlp(nbl), iprint, n14(nbl), maxbl0, maxseg0,
                         nbcj0_ref, nbck0_ref, nbci0_ref, nbcjdim_ref, nbckdim_ref, nbcidim_ref,
                         jbcinfo_ref, kbcinfo_ref, ibcinfo_ref);
                n14_plus = n14(nbl) + 1;
                inmx = std::max({jdim, kdim, idim});
                iwk1 = 1;
                iwk2 = iwk1 + inmx;
                iwk3 = iwk2 + inmx;
                iwk4 = iwk3 + inmx;
                iwk5 = iwk4 + inmx;
                iwk6 = iwk5 + inmx;
                iwk7 = iwk6 + inmx;
                iwk8 = iwk7 + inmx;
                iwk9 = iwk8 + inmx;
                iwk10= iwk9 + inmx;
                iwk11= iwk10+ n14_plus;
                iwk12= iwk11+ n14_plus;
                iwk13= iwk12+ n14_plus;
                iwk14= iwk13+ n14_plus;
                iwk15= iwk14+ n14_plus;
                iwk16= iwk15+ n14_plus;
                iwk17= iwk16+ n14_plus;
                iwk18= iwk17+ n14_plus;
                iwk19= iwk18+ n14_plus;
                iwk20= iwk19+ n14_plus;
                iwk21= iwk20+ n14_plus;
                iwk22= iwk21+ n14_plus;
                iwk23= iwk22+ n14_plus;
                iwk24= iwk23+ n14_plus;
                iwk25= iwk24+ n14_plus;
                iwk26= iwk25+ n14_plus;
                iwk27= iwk26+ n14_plus;
                iwk28= iwk27+ n14_plus;
                iwk29= iwk28+ n14_plus;
                iwk30= iwk29+ n14_plus;
                iwk31= iwk30+ n14_plus;
                iwk32= iwk31+ n14_plus;
                iwk33= iwk32+ n14_plus;
                iwk34= iwk33+ n14_plus;
                iwk35= iwk34+ n14_plus;
                iwk36= iwk35+ n14_plus;
                iwk37= iwk36+ n14_plus;
                iwk38= iwk37+ n14_plus;
                iwk39= iwk38+ n14_plus;
                iwk40= iwk39+ n14_plus;
                iwk41= iwk40+ n14_plus;
                iwk42= iwk41+ n14_plus;
                iwk43= iwk42+ n14_plus;
                iwk44= iwk43+ n14_plus;
                iwk45= iwk44+ n14_plus;
                iwk46= iwk45+ n14_plus;
                iwk47= iwk46+ n14_plus;
                iwk48= iwk47+ n14_plus;
                iwk49= iwk48+ n14_plus;
                iwk50= iwk49+ n14_plus;
                iwk51= iwk50+ n14_plus;
                iwk52= iwk51+ n14_plus;
                iwk53= iwk52+ n14_plus;
                iwk54= iwk53+ n14_plus;
                iwk55= iwk54+ n14_plus;
                iwk56= iwk55+ n14_plus;
                iwk57= iwk56+ n14_plus;
                iwk58= iwk57+ n14_plus;
                iwk59= iwk58+ n14_plus;
                iwk60= iwk59+ n14_plus;
                iwk61= iwk60+ n14_plus;
                iwk62= iwk61+ n14_plus;
                iwk63= iwk62+ n14_plus;
                iwk64= iwk63+ n14_plus;
                iwk65= iwk64+ n14_plus;
                iwk66= iwk65+ n14_plus;
                iwk67= iwk66+ n14_plus;
                iwk68= iwk67+ n14_plus;
                iwk69= iwk68+ n14_plus;
                iwk70= iwk69+ n14_plus;
                iwk71= iwk70+ n14_plus;
                iwk72= iwk71+ n14_plus;
                iwk73= iwk72+ n14_plus;
                iwk74= iwk73+ n14_plus;
                iwk75= iwk74+ n14_plus;
                iwk76= iwk75+ n14_plus;
                iwk77= iwk76+ n14_plus;
                iwk78= iwk77+ n14_plus;
                iwk79= iwk78+ n14_plus;
                iwk80= iwk79+ n14_plus;
                iwk81= iwk80+ n14_plus;
                iwk82= iwk81+ n14_plus;
                itemp= iwk82+ n14_plus;
                need = std::max(need, itemp);
                ii = mblk2nd(nbl);
                mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
                fprintf(f66, " 9(blomax)       itemp, need = %d %d\n", itemp, need);
                fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
            }
        }

        label_8500:
        (void)0; // label target
    } // end iseq loop


    // Requirements for diagonal solvers and metric
    for (nbl = 1; nbl <= nblock; nbl++) {
        {
            FortranArray2DRef<int> lw_ref = lw.ref();
            FortranArray2DRef<int> lw2_ref = lw2.ref();
            lead_ns::lead(nbl, lw_ref, lw2_ref, maxbl0);
        }
        ndim = std::max(kdim, idim);
        if (ifdsg(nbl,2) == 0 || idiag[1] == 0) {
            itemp = ltot(nbl) + jdim*ndim*145;
            need  = std::max(need, itemp);
            ii = mblk2nd(nbl);
            mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
            fprintf(f66, " 13(5x5 J)       itemp, need = %d %d\n", itemp, need);
            fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
        } else {
            itemp = ltot(nbl) + jdim*ndim*35;
            need  = std::max(need, itemp);
            ii = mblk2nd(nbl);
            mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
            fprintf(f66, " 14(diagonal J)  itemp, need = %d %d\n", itemp, need);
            fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
        }
        if (ifdsg(nbl,3) == 0 || idiag[2] == 0) {
            itemp = ltot(nbl) + jdim*ndim*145;
            need  = std::max(need, itemp);
            ii = mblk2nd(nbl);
            mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
            fprintf(f66, " 15(5x5 K)       itemp, need = %d %d\n", itemp, need);
            fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
        } else {
            itemp = ltot(nbl) + jdim*ndim*35;
            need  = std::max(need, itemp);
            ii = mblk2nd(nbl);
            mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
            fprintf(f66, " 16(diagonal K)  itemp, need = %d %d\n", itemp, need);
            fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
        }
        if (std::abs(i2d) != 1) {
            if (ifdsg(nbl,1) == 0 || idiag[0] == 0) {
                itemp = ltot(nbl) + jdim*ndim*145;
                need  = std::max(need, itemp);
                ii = mblk2nd(nbl);
                mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
                fprintf(f66, " 17(5x5 I)       itemp, need = %d %d\n", itemp, need);
                fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
            } else {
                itemp = ltot(nbl) + jdim*ndim*35;
                need  = std::max(need, itemp);
                ii = mblk2nd(nbl);
                mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
                fprintf(f66, " 18(diagonal I)  itemp, need = %d %d\n", itemp, need);
                fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
            }
        }
        itemp = ltot(nbl) + jdim*kdim*6 + jdim*kdim*idim*5;
        need  = std::max(need, itemp);
        ii = mblk2nd(nbl);
        mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
        fprintf(f66, " 19(metric)      itemp, need = %d %d\n", itemp, need);
        fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));

        iuns = std::max({itrans(nbl), irotat(nbl), idefrm(nbl)});
        if (iuns > 0) {
            itemp = ltot(nbl) + jdim*kdim*idim*3 + jdim*kdim*3*2 + kdim*idim*3*2 + jdim*idim*3*2;
            need  = std::max(need, itemp);
            ii = mblk2nd(nbl);
            mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
            fprintf(f66, " 20(trans)       itemp, need = %d %d\n", itemp, need);
            fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));

            itemp = ltot(nbl) + 2*jdim*kdim*idim*3 + jdim*kdim*3*2 + kdim*idim*3*2 + jdim*idim*3*2;
            need  = std::max(need, itemp);
            ii = mblk2nd(nbl);
            mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
            fprintf(f66, " 21(rotate)      itemp, need = %d %d\n", itemp, need);
            fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));

            itemp = ltot(nbl) + 2*jdim*kdim*idim*3 + jdim*kdim*6 + jdim*kdim*idim*5 + jdim*kdim*3*2 + kdim*idim*3*2 + jdim*idim*3*2;
            need  = std::max(need, itemp);
            ii = mblk2nd(nbl);
            mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
            fprintf(f66, " 22(tmetric)     itemp, need = %d %d\n", itemp, need);
            fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));

            itemp = ltot(nbl) + jdim*kdim*idim*3 + jdim*kdim*3*2 + kdim*idim*3*2 + jdim*idim*3*2;
            need  = std::max(need, itemp);
            ii = mblk2nd(nbl);
            mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
            fprintf(f66, " 23(xtbatb)      itemp, need = %d %d\n", itemp, need);
            fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
        }
    } // end nbl loop for diagonal/metric

    // Requirements for plot3d output
    for (n = 1; n <= nplot3d; n++) {
        nbl = inpl3d(n, 1);
        i1 = inpl3d(n, 3); i2 = inpl3d(n, 4); i3 = inpl3d(n, 5);
        j1 = inpl3d(n, 6); j2 = inpl3d(n, 7); j3 = inpl3d(n, 8);
        k1 = inpl3d(n, 9); k2 = inpl3d(n,10); k3 = inpl3d(n,11);
        {
            FortranArray2DRef<int> lw_ref = lw.ref();
            FortranArray2DRef<int> lw2_ref = lw2.ref();
            lead_ns::lead(nbl, lw_ref, lw2_ref, maxbl0);
        }
        if (inpl3d(n,2) <= 0) {
            jdw = (j2-j1)/j3 + 1; kdw = (k2-k1)/k3 + 1; idw = (i2-i1)/i3 + 1;
            itemp = jdw*kdw*idw*14 + jdim*kdim*idim*2;
        } else if (inpl3d(n,2) == 1 || inpl3d(n,2) > 2) {
            i2 = std::min(idim-1,i2); j2 = std::min(jdim-1,j2); k2 = std::min(kdim-1,k2);
            i1 = std::min(idim-1,i1); j1 = std::min(jdim-1,j1); k1 = std::min(kdim-1,k1);
            jdw = (j2-j1)/j3 + 1; kdw = (k2-k1)/k3 + 1; idw = (i2-i1)/i3 + 1;
            itemp = jdw*kdw*idw*9 + jdim*kdim*idim;
        } else {
            if (ivmx > 1) {
                i2 = std::min(idim-1,i2); j2 = std::min(jdim-1,j2); k2 = std::min(kdim-1,k2);
                i1 = std::min(idim-1,i1); j1 = std::min(jdim-1,j1); k1 = std::min(kdim-1,k1);
                jdw = (j2-j1)/j3 + 1; kdw = (k2-k1)/k3 + 1; idw = (i2-i1)/i3 + 1;
                if (ivmx == 8 || ivmx == 9 || ivmx >= 11)
                    itemp = jdw*kdw*idw*9 + jdim*kdim*idim + jdim*kdim*9;
                else
                    itemp = jdw*kdw*idw*9 + jdim*kdim*idim + (jdim-1)*(kdim-1)*(idim-1)*9 + jdim*kdim*9;
            } else {
                itemp = 0;
            }
        }
        need  = std::max(need, itemp);
        ii = mblk2nd(nbl);
        mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
        fprintf(f66, " 23(plot3d)      itemp, need = %d %d\n", itemp, need);
        fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
    }

    // Requirements for qoutavg (plot3davg) output for iteravg>0
    if (iteravg > 0) {
        for (n = 1; n <= nplot3d; n++) {
            nbl = inpl3d(n, 1);
            {
                FortranArray2DRef<int> lw_ref = lw.ref();
                FortranArray2DRef<int> lw2_ref = lw2.ref();
                lead_ns::lead(nbl, lw_ref, lw2_ref, maxbl0);
            }
            i1=1; i2=idim; i3=1; j1=1; j2=jdim; j3=1; k1=1; k2=kdim; k3=1;
            jdw = (j2-j1)/j3 + 1; kdw = (k2-k1)/k3 + 1; idw = (i2-i1)/i3 + 1;
            itemp = jdw*kdw*idw*14 + jdim*kdim*idim*2;
            need  = std::max(need, itemp);
            ii = mblk2nd(nbl);
            mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
            fprintf(f66, " 23(plot3davg)   itemp, need = %d %d\n", itemp, need);
            fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
        }
    }

    // Requirements for print output
    for (n = 1; n <= nprint; n++) {
        nbl = inpr(n, 1);
        i1 = inpr(n,3); i2 = inpr(n,4); i3 = inpr(n,5);
        j1 = inpr(n,6); j2 = inpr(n,7); j3 = inpr(n,8);
        k1 = inpr(n,9); k2 = inpr(n,10); k3 = inpr(n,11);
        {
            FortranArray2DRef<int> lw_ref = lw.ref();
            FortranArray2DRef<int> lw2_ref = lw2.ref();
            lead_ns::lead(nbl, lw_ref, lw2_ref, maxbl0);
        }
        if (inpr(n,2) == 0) {
            jdw = (j2-j1)/j3 + 1; kdw = (k2-k1)/k3 + 1; idw = (i2-i1)/i3 + 1;
            itemp = jdw*kdw*idw*17 + jdim*kdim*idim*2;
        } else {
            i2 = std::min(idim-1,i2); j2 = std::min(jdim-1,j2); k2 = std::min(kdim-1,k2);
            i1 = std::min(idim-1,i1); j1 = std::min(jdim-1,j1); k1 = std::min(kdim-1,k1);
            jdw = (j2-j1)/j3 + 1; kdw = (k2-k1)/k3 + 1; idw = (i2-i1)/i3 + 1;
            itemp = jdw*kdw*idw*9 + jdim*kdim*idim;
        }
        need  = std::max(need, itemp);
        ii = mblk2nd(nbl);
        mem_req_node(ii) = std::max(mem_req_node(ii), itemp);
        fprintf(f66, " 23(print)       itemp, need = %d %d\n", itemp, need);
        fprintf(f66, "       parallel: itemp, need = %d %d\n", itemp, mem_req_node(ii));
    }

    // Determine i/o buffer size
    libufdim = ibufdim0;


    // Determine i/o buffer size for Baldwin-Lomax
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nblg(igrid);
        {
            FortranArray2DRef<int> lw_ref = lw.ref();
            FortranArray2DRef<int> lw2_ref = lw2.ref();
            lead_ns::lead(nbl, lw_ref, lw2_ref, maxbl0);
        }
        n14(nbl) = 0;
        ivisc[0] = iviscg(nbl, 1);
        ivisc[1] = iviscg(nbl, 2);
        ivisc[2] = iviscg(nbl, 3);
        if (ivisc[2] == 2 || ivisc[2] == 3 || ivisc[1] == 2 || ivisc[1] == 3 || ivisc[0] == 2 || ivisc[0] == 3) {
            ipw = 2;
            FortranArray1DRef<int> nbcj0_ref = nbcj0.ref();
            FortranArray1DRef<int> nbck0_ref = nbck0.ref();
            FortranArray1DRef<int> nbci0_ref = nbci0.ref();
            FortranArray1DRef<int> nbcjdim_ref = nbcjdim.ref();
            FortranArray1DRef<int> nbckdim_ref = nbckdim.ref();
            FortranArray1DRef<int> nbcidim_ref = nbcidim.ref();
            FortranArray4DRef<int> jbcinfo_ref = jbcinfo.ref();
            FortranArray4DRef<int> kbcinfo_ref = kbcinfo.ref();
            FortranArray4DRef<int> ibcinfo_ref = ibcinfo.ref();
            cntblmx(nbl, jdim, kdim, idim, iovrlp(nbl), ipw, n14(nbl), maxbl0, maxseg0,
                     nbcj0_ref, nbck0_ref, nbci0_ref, nbcjdim_ref, nbckdim_ref, nbcidim_ref,
                     jbcinfo_ref, kbcinfo_ref, ibcinfo_ref);
            n14_plus = n14(nbl) + 1;
            libufdim = std::max(n14_plus, libufdim);
        }
    }

    // Async message passing requirements
    icount_pat  = 0;
    icount_blk  = 0;
    icount_prd  = 0;
    icount_emb  = 0;
    icount_dpat = 0;
    icount_rad  = 0;
    for (levl = 1; levl <= levt; levl++) {
        is_blk[levl-1] = icount_blk + 1;
        is_pat[levl-1] = icount_pat + 1;
        is_prd[levl-1] = icount_prd + 1;
        is_emb[levl-1] = icount_emb + 1;
        ie_blk[levl-1] = is_blk[levl-1] - 1;
        ie_pat[levl-1] = is_pat[levl-1] - 1;
        ie_prd[levl-1] = is_prd[levl-1] - 1;
        ie_emb[levl-1] = is_emb[levl-1] - 1;
        for (nbl = 1; nbl <= nblock; nbl++) {
            if (levl != levelg(nbl)) continue;
            icount_pat1 = icount_pat;
            {
                FortranArray2DRef<int> lw_ref = lw.ref();
                FortranArray2DRef<int> iindex_ref = iindex.ref();
                FortranArray2DRef<int> isav_pat_ref = isav_pat.ref();
                FortranArray3DRef<int> isav_pat_b_ref = isav_pat_b.ref();
                FortranArray1DRef<int> jjmax1_ref = jjmax1.ref();
                FortranArray1DRef<int> kkmax1_ref = kkmax1.ref();
                FortranArray1DRef<int> iiint1_ref = iiint1.ref();
                FortranArray1DRef<int> iiint2_ref = iiint2.ref();
                FortranArray1DRef<int> jdimg_ref = jdimg.ref();
                FortranArray1DRef<int> kdimg_ref = kdimg.ref();
                FortranArray1DRef<int> idimg_ref = idimg.ref();
                int nintr2 = ninter;
                pre_patch_ns::pre_patch(nbl, lw_ref, icount_pat, nintr2,
                                        iindex_ref, intmax0, nsub10, isav_pat_ref,
                                        isav_pat_b_ref, jjmax1_ref, kkmax1_ref,
                                        iiint1_ref, iiint2_ref, maxbl0, jdimg_ref,
                                        kdimg_ref, idimg_ref, ierrflg);
            }
            icount_blk1 = icount_blk;
            {
                FortranArray2DRef<int> lw_ref = lw.ref();
                FortranArray2DRef<int> isav_blk_ref = isav_blk.ref();
                FortranArray2DRef<int> nblk_ref = nblk.ref();
                FortranArray3DRef<int> limblk_ref = limblk.ref();
                FortranArray3DRef<int> isva_ref = isva.ref();
                FortranArray1DRef<int> nblon_ref = nblon.ref();
                FortranArray1DRef<int> idimg_ref = idimg.ref();
                FortranArray1DRef<int> jdimg_ref = jdimg.ref();
                FortranArray1DRef<int> kdimg_ref = kdimg.ref();
                pre_blockbc_ns::pre_blockbc(nbl, lw_ref, icount_blk, idimg_ref,
                                            jdimg_ref, kdimg_ref, isav_blk_ref, nblk_ref, nbli,
                                            limblk_ref, isva_ref, nblon_ref, mxbli0, nou, bou, nbuf0,
                                            ibufdim0, myid, maxbl0, ierrflg);
            }
            icount_prd1 = icount_prd;
            {
                FortranArray2DRef<int> lw_ref = lw.ref();
                FortranArray2DRef<int> lw2_ref = lw2.ref();
                FortranArray1DRef<int> nbcj0_ref = nbcj0.ref();
                FortranArray1DRef<int> nbck0_ref = nbck0.ref();
                FortranArray1DRef<int> nbci0_ref = nbci0.ref();
                FortranArray1DRef<int> nbcjdim_ref = nbcjdim.ref();
                FortranArray1DRef<int> nbckdim_ref = nbckdim.ref();
                FortranArray1DRef<int> nbcidim_ref = nbcidim.ref();
                FortranArray4DRef<int> jbcinfo_ref = jbcinfo.ref();
                FortranArray4DRef<int> kbcinfo_ref = kbcinfo.ref();
                FortranArray4DRef<int> ibcinfo_ref = ibcinfo.ref();
                FortranArray1DRef<int> igridg_ref = igridg.ref();
                FortranArray1DRef<int> jdimg_ref = jdimg.ref();
                FortranArray1DRef<int> kdimg_ref = kdimg.ref();
                FortranArray1DRef<int> idimg_ref = idimg.ref();
                FortranArray2DRef<int> isav_prd_ref = isav_prd.ref();
                FortranArray4DRef<double> bcvali_ref = bcvali.ref();
                FortranArray4DRef<double> bcvalj_ref = bcvalj.ref();
                FortranArray4DRef<double> bcvalk_ref = bcvalk.ref();
                FortranArray1DRef<int> nblg_ref = nblg.ref();
                pre_period_ns::pre_period(nbl, lw_ref, lw2_ref, icount_prd,
                                          maxbl0, maxseg0, lbcprd0,
                                          nbcj0_ref, nbck0_ref, nbci0_ref, nbcjdim_ref, nbckdim_ref,
                                          nbcidim_ref, jbcinfo_ref, kbcinfo_ref, ibcinfo_ref,
                                          igridg_ref, jdimg_ref, kdimg_ref, idimg_ref, isav_prd_ref,
                                          FortranArray1DRef<int>(is_prd, 5), FortranArray1DRef<int>(ie_prd, 5), nbcprd, nou, bou, nbuf0, ibufdim0,
                                          bcvali_ref, bcvalj_ref, bcvalk_ref, myid, nblg_ref, maxgr0,
                                          ierrflg);
            }
            icount_emb1 = icount_emb;
            {
                FortranArray2DRef<int> lw_ref = lw.ref();
                FortranArray2DRef<int> lw2_ref = lw2.ref();
                FortranArray1DRef<int> nbcj0_ref = nbcj0.ref();
                FortranArray1DRef<int> nbck0_ref = nbck0.ref();
                FortranArray1DRef<int> nbci0_ref = nbci0.ref();
                FortranArray1DRef<int> nbcjdim_ref = nbcjdim.ref();
                FortranArray1DRef<int> nbckdim_ref = nbckdim.ref();
                FortranArray1DRef<int> nbcidim_ref = nbcidim.ref();
                FortranArray4DRef<int> jbcinfo_ref = jbcinfo.ref();
                FortranArray4DRef<int> kbcinfo_ref = kbcinfo.ref();
                FortranArray4DRef<int> ibcinfo_ref = ibcinfo.ref();
                FortranArray1DRef<int> iemg_ref = iemg.ref();
                FortranArray1DRef<int> igridg_ref = igridg.ref();
                FortranArray1DRef<int> jdimg_ref = jdimg.ref();
                FortranArray1DRef<int> kdimg_ref = kdimg.ref();
                FortranArray1DRef<int> idimg_ref = idimg.ref();
                FortranArray2DRef<int> isav_emb_ref = isav_emb.ref();
                FortranArray1DRef<int> nblg_ref = nblg.ref();
                pre_embed_ns::pre_embed(nbl, lw_ref, lw2_ref, icount_emb,
                                        maxbl0, maxseg0, lbcemb0,
                                        nbcj0_ref, nbck0_ref, nbci0_ref, nbcjdim_ref, nbckdim_ref,
                                        nbcidim_ref, jbcinfo_ref, kbcinfo_ref, ibcinfo_ref, iemg_ref,
                                        igridg_ref, jdimg_ref, kdimg_ref, idimg_ref, isav_emb_ref,
                                        FortranArray1DRef<int>(is_emb, 5), FortranArray1DRef<int>(ie_emb, 5), nbcemb, nou, bou, nbuf0, ibufdim0,
                                        myid, maxgr0, ierrflg);
            }
            if (iunst > 0) {
                FortranArray2DRef<int> lw_ref = lw.ref();
                FortranArray2DRef<int> iindx_ref = iindx.ref();
                FortranArray2DRef<int> isav_dpat_ref = isav_dpat.ref();
                FortranArray3DRef<int> isav_dpat_b_ref = isav_dpat_b.ref();
                FortranArray1DRef<int> jjmax1_ref = jjmax1.ref();
                FortranArray1DRef<int> kkmax1_ref = kkmax1.ref();
                FortranArray1DRef<int> iiint1_ref = iiint1.ref();
                FortranArray1DRef<int> iiint2_ref = iiint2.ref();
                FortranArray1DRef<int> jdimg_ref = jdimg.ref();
                FortranArray1DRef<int> kdimg_ref = kdimg.ref();
                FortranArray1DRef<int> idimg_ref = idimg.ref();
                int nintr2 = 0;
                pre_patch_ns::pre_patch(nbl, lw_ref, icount_dpat, nintr2,
                                        iindx_ref, intmx0, msub10, isav_dpat_ref,
                                        isav_dpat_b_ref, jjmax1_ref, kkmax1_ref,
                                        iiint1_ref, iiint2_ref, maxbl0, jdimg_ref,
                                        kdimg_ref, idimg_ref, ierrflg);
            }
            if (icount_blk > icount_blk1) ie_blk[levl-1] = icount_blk;
            if (icount_pat > icount_pat1) ie_pat[levl-1] = icount_pat;
            if (icount_prd > icount_prd1) ie_prd[levl-1] = icount_prd;
            if (icount_emb > icount_emb1) ie_emb[levl-1] = icount_emb;
        }
    }

    llbcprd = std::max(icount_prd, 1);
    llbcemb = std::max(icount_emb, 1);
    llbcrad = std::max(icount_rad, 1);

    if (icount_blk > 0 || icount_pat > 0 || icount_prd > 0 || icount_emb > 0 || icount_dpat > 0 || icount_rad > 0) {
        fprintf(f66, " -async are requirements for asynchronous message passing, by node\n");
    }


    // For 1-1 interfaces
    if (icount_blk > 0) {
        for (myid0 = 1; myid0 <= nodel; myid0++) {
            itemp = 1;
            // for wk(kqintl) allocation
            for (lcnt = 1; lcnt <= icount_blk; lcnt++) {
                ic_blk  = isav_blk(lcnt, 4);
                in_blk  = isav_blk(lcnt, 5);
                nd_dest = mblk2nd(ic_blk);
                nd_srce = mblk2nd(in_blk);
                if (nd_srce != myid0) {
                    if (nd_dest == myid0) {
                        jface = isav_blk(lcnt, 6);
                        idimn = idimg(in_blk); jdimn = jdimg(in_blk); kdimn = kdimg(in_blk);
                        if (jface == 1) maxdims = jdimn*kdimn;
                        if (jface == 2) maxdims = kdimn*idimn;
                        if (jface == 3) maxdims = jdimn*idimn;
                        itemp = itemp + maxdims*13;
                        if (ivmx >= 2) itemp = itemp + maxdims*2;
                        if (ivmx >= 4) itemp = itemp + maxdims*7*2;
                        if (ivmx >= 1) itemp = itemp + maxdims*2;
                    }
                }
            }
            // for wk(ktl) allocation
            iskipz = 0;
            if (iskipz == 0) {
                for (lcnt = 1; lcnt <= icount_blk; lcnt++) {
                    ic_blk  = isav_blk(lcnt, 4);
                    in_blk  = isav_blk(lcnt, 5);
                    nd_dest = mblk2nd(ic_blk);
                    nd_srce = mblk2nd(in_blk);
                    if (nd_dest != myid0) {
                        if (nd_srce == myid0) {
                            jface = isav_blk(lcnt, 6);
                            idimn = idimg(in_blk); jdimn = jdimg(in_blk); kdimn = kdimg(in_blk);
                            if (jface == 1) maxdims = jdimn*kdimn;
                            if (jface == 2) maxdims = kdimn*idimn;
                            if (jface == 3) maxdims = jdimn*idimn;
                            itemp = itemp + maxdims*13;
                            if (ivmx >= 2) itemp = itemp + maxdims*2;
                            if (ivmx >= 4) itemp = itemp + maxdims*7*2;
                            if (ivmx >= 1) itemp = itemp + maxdims*2;
                        }
                    }
                }
            }
            mem_req_node(myid0) = std::max(mem_req_node(myid0), itemp);
            fprintf(f66, " 22(1:1-async)   itemp, need = %d %d\n", itemp, mem_req_node(myid0));
        }
    }

    // For patched interfaces
    if (icount_pat > 0) {
        for (myid0 = 1; myid0 <= nodel; myid0++) {
            itemp = 1;
            // for wk(kqintl) allocation
            for (intr = 1; intr <= icount_pat; intr++) {
                lmax_val = isav_pat(intr, 2);
                for (nf = 1; nf <= lmax_val; nf++) {
                    itb  = isav_pat(intr, 1);
                    ifb  = isav_pat_b(intr, nf, 1);
                    ityp = isav_pat_b(intr, nf, 2);
                    if (mblk2nd(ifb) != myid0) {
                        if (mblk2nd(itb) == myid0) {
                            if (ityp == 1) { jmax1 = jdimg(ifb); kmax1 = kdimg(ifb); }
                            else if (ityp == 2) { jmax1 = kdimg(ifb); kmax1 = idimg(ifb); }
                            else if (ityp == 3) { jmax1 = jdimg(ifb); kmax1 = idimg(ifb); }
                            itemp = itemp + jmax1*kmax1*5*2;
                            if (ivmx >= 2) itemp = itemp + jmax1*kmax1*1*2;
                            if (ivmx >= 4) itemp = itemp + jmax1*kmax1*7*2;
                        }
                    }
                }
            }
            // for wk(ktl) allocation
            for (intr = 1; intr <= icount_pat; intr++) {
                lmax_val = isav_pat(intr, 2);
                for (nf = 1; nf <= lmax_val; nf++) {
                    itb  = isav_pat(intr, 1);
                    ifb  = isav_pat_b(intr, nf, 1);
                    ityp = isav_pat_b(intr, nf, 2);
                    if (mblk2nd(ifb) == myid0) {
                        if (mblk2nd(itb) != myid0) {
                            if (ityp == 1) { jmax1 = jdimg(ifb); kmax1 = kdimg(ifb); }
                            else if (ityp == 2) { jmax1 = kdimg(ifb); kmax1 = idimg(ifb); }
                            else if (ityp == 3) { jmax1 = jdimg(ifb); kmax1 = idimg(ifb); }
                            itemp = itemp + jmax1*kmax1*5*2;
                            if (ivmx >= 2) itemp = itemp + jmax1*kmax1*1*2;
                            if (ivmx >= 4) itemp = itemp + jmax1*kmax1*7*2;
                        }
                    }
                }
            }
            // for subroutine int2
            itemp1 = 0;
            for (intr = 1; intr <= icount_pat; intr++) {
                lmax_val = isav_pat(intr, 2);
                for (nf = 1; nf <= lmax_val; nf++) {
                    itb  = isav_pat(intr, 1);
                    ifb  = isav_pat_b(intr, nf, 1);
                    ityp = isav_pat_b(intr, nf, 2);
                    if (ityp == 1) { jmax1 = jdimg(ifb); kmax1 = kdimg(ifb); }
                    else if (ityp == 2) { jmax1 = kdimg(ifb); kmax1 = idimg(ifb); }
                    else if (ityp == 3) { jmax1 = jdimg(ifb); kmax1 = idimg(ifb); }
                    itemp1 = std::max(itemp1, 4*jmax1*kmax1);
                }
            }
            itemp = itemp + itemp1;
            mem_req_node(myid0) = std::max(mem_req_node(myid0), itemp);
            fprintf(f66, " 23(patch-async) itemp, need = %d %d\n", itemp, mem_req_node(myid0));
        }
    }


    // For periodic interfaces
    if (icount_prd > 0) {
        for (myid0 = 1; myid0 <= nodel; myid0++) {
            itemp = 1;
            for (lcnt = 1; lcnt <= icount_prd; lcnt++) {
                int nbll2 = isav_prd(lcnt, 1);
                nseg  = isav_prd(lcnt, 11);
                nface = isav_prd(lcnt, 2);
                int nblp2 = isav_prd(lcnt, 12);
                nd_recv = mblk2nd(nbll2);
                nd_srce = mblk2nd(nblp2);
                if (nd_recv == myid0) {
                    if (nd_srce != myid0) {
                        idimp = idimg(nblp2); jdimp = jdimg(nblp2); kdimp = kdimg(nblp2);
                        if (nface == 1 || nface == 2) maxdims = jdimp*kdimp;
                        if (nface == 3 || nface == 4) maxdims = kdimp*idimp;
                        if (nface == 5 || nface == 6) maxdims = jdimp*idimp;
                        itemp = itemp + maxdims*13;
                        if (ivmx >= 2) itemp = itemp + maxdims*2;
                        if (ivmx >= 4) itemp = itemp + maxdims*7*2;
                    }
                }
            }
            for (lcnt = 1; lcnt <= icount_prd; lcnt++) {
                int nbll2 = isav_prd(lcnt, 1);
                int nblp2 = isav_prd(lcnt, 12);
                nface = isav_prd(lcnt, 2);
                nd_recv = mblk2nd(nbll2);
                nd_srce = mblk2nd(nblp2);
                if (nd_srce == myid0) {
                    if (nd_recv != myid0) {
                        idimp = idimg(nblp2); jdimp = jdimg(nblp2); kdimp = kdimg(nblp2);
                        if (nface == 1 || nface == 2) maxdims = jdimp*kdimp;
                        if (nface == 3 || nface == 4) maxdims = kdimp*idimp;
                        if (nface == 5 || nface == 6) maxdims = jdimp*idimp;
                        itemp = itemp + maxdims*13;
                        if (ivmx >= 2) itemp = itemp + maxdims*2;
                        if (ivmx >= 4) itemp = itemp + maxdims*7*2;
                    }
                }
            }
            mem_req_node(myid0) = std::max(mem_req_node(myid0), itemp);
            fprintf(f66, " 24(period-async)itemp, need = %d %d\n", itemp, mem_req_node(myid0));
        }
    }

    // For embedded interfaces
    if (icount_emb > 0) {
        for (myid0 = 1; myid0 <= nodel; myid0++) {
            itemp = 1;
            for (lcnt = 1; lcnt <= icount_emb; lcnt++) {
                int nblf2 = isav_emb(lcnt, 1);
                nface = isav_emb(lcnt, 2);
                int nblc2 = isav_emb(lcnt, 9);
                nd_recv = mblk2nd(nblf2);
                nd_srce = mblk2nd(nblc2);
                if (nd_recv == myid0) {
                    if (nd_srce != myid0) {
                        idimc = idimg(nblc2); jdimc = jdimg(nblc2); kdimc = kdimg(nblc2);
                        if (nface == 1 || nface == 2) maxdims = jdimc*kdimc;
                        if (nface == 3 || nface == 4) maxdims = kdimc*idimc;
                        if (nface == 5 || nface == 6) maxdims = jdimc*idimc;
                        itemp = itemp + maxdims*15;
                        if (ivmx >= 2) itemp = itemp + maxdims*3;
                        if (ivmx >= 4) itemp = itemp + maxdims*7*3;
                    }
                }
            }
            for (lcnt = 1; lcnt <= icount_emb; lcnt++) {
                int nblf2 = isav_emb(lcnt, 1);
                int nblc2 = isav_emb(lcnt, 9);
                nface = isav_emb(lcnt, 2);
                nd_recv = mblk2nd(nblf2);
                nd_srce = mblk2nd(nblc2);
                if (nd_srce == myid0) {
                    if (nd_recv != myid0) {
                        idimc = idimg(nblc2); jdimc = jdimg(nblc2); kdimc = kdimg(nblc2);
                        if (nface == 1 || nface == 2) maxdims = jdimc*kdimc;
                        if (nface == 3 || nface == 4) maxdims = kdimc*idimc;
                        if (nface == 5 || nface == 6) maxdims = jdimc*idimc;
                        itemp = itemp + maxdims*15;
                        if (ivmx >= 2) itemp = itemp + maxdims*3;
                        if (ivmx >= 4) itemp = itemp + maxdims*7*3;
                    }
                }
            }
            mem_req_node(myid0) = std::max(mem_req_node(myid0), itemp);
            fprintf(f66, " 25(embed-async) itemp, need = %d %d\n", itemp, mem_req_node(myid0));
        }
    }

    // For dynamic grids
    if (iunst > 0) {
        for (myid0 = 1; myid0 <= nodel; myid0++) {
            itemp = 1;
            for (intr = 1; intr <= icount_dpat; intr++) {
                lmax_val = isav_dpat(intr, 2);
                for (nf = 1; nf <= lmax_val; nf++) {
                    itb = isav_dpat(intr, 1);
                    ifb = isav_dpat_b(intr, nf, 1);
                    int mdim1 = isav_dpat_b(intr, nf, 3);
                    int ndim1 = isav_dpat_b(intr, nf, 4);
                    if (mblk2nd(ifb) != myid0) {
                        if (mblk2nd(itb) == myid0) {
                            itemp = itemp + 3*mdim1*ndim1;
                        }
                    }
                }
            }
            mem_req_node(myid0) = std::max(mem_req_node(myid0), itemp);
            fprintf(f66, " 26(dyn-async)   itemp, need = %d %d\n", itemp, mem_req_node(myid0));
        }
    }

    // Also need space for int2 in sequential case
    if (ninter > 0) {
        need = std::max(need, itemp1);
        fprintf(f66, " 27(int2)        itemp, need = %d %d\n", itemp1, need);
    }

    // Add space for timesave array in mgblk
    need = need + lmaxbl;
    for (ii = 1; ii <= nodel; ii++) {
        mem_req_node(ii) = mem_req_node(ii) + lmaxbl;
    }

    mem_wk_max = mem_req_node(1);
    for (myid0 = 2; myid0 <= nodel; myid0++) {
        if (mem_req_node(myid0) > mem_wk_max) mem_wk_max = mem_req_node(myid0);
    }
    tot_nod = mem_w_max + mem_wk_max;

    fprintf(f66, "\n");
    fprintf(f66, "***********************************************\n");
    fprintf(f66, "\n");
    fprintf(f66, "      TEMPORARY STORAGE REQUIREMENTS - WK\n");
    fprintf(f66, "            SUMMARY FOR ALL BUILDS\n");
    fprintf(f66, "\n");
    fprintf(f66, "***********************************************\n");
    fprintf(f66, "\n memory (mw) for wk storage (sequential) = %9.4f\n", (float)need/1.e6f);
    fprintf(f66, "\n memory (mw) for wk storage (nodes)\n");
    fprintf(f66, "\n node   total points  memory (mw)\n");
    fmemsum = 0.0;
    nptsum  = 0;
    for (i = 1; i <= nodel; i++) {
        fprintf(f66, "%6d       %8d    %9.4f\n", i, no_of_points(i), (float)(mem_req_node(i)/1.e6));
        fmemsum = fmemsum + mem_req_node(i)/1.e6;
        nptsum  = nptsum + no_of_points(i);
        mem_req_node(i) = 0;
    }
    fprintf(f66, "\n   maximum =  %8d    %9.4f\n", npts_max, (float)mem_wk_max/1.e6f);


    // Evaluate requirements for temporary integer array size (iwork/iwk array)
    // in mgbl, to call global
    needi_glo = lmaxbl*3 + 11*lnplts;
    needi = std::max(needi, needi_glo);
    // in mgbl, to call trnsfr_vals
    needi_trn = 1550*lmaxbl + 119 + 9*lmaxbl*lmxsegdg
              + 48*lmaxbl*lmxseg + 19*lmxbli + 22*lnplts
              + 4*lmaxgr + 9*lmaxcs;
    needi = std::max(needi, needi_trn);
    // in mgbl, to call pointers
    needi_pnt = lmaxbl*2;
    needi = std::max(needi, needi_pnt);
    // in mgbl, to call setup
    needi_set = 3*lmaxgr + lmaxbl*8;
    needi = std::max(needi, needi_set);
    // in mgbl, to call mgblk
    needi_mgblk = lmaxbl*7*3 + lmaxbl*8;
    needi = std::max(needi, needi_mgblk);
    // in mgbl, to call qout
    needi_qout = 3*lnplts + lmaxbl;
    needi = std::max(needi, needi_qout);
    // in mgbl, to call yplusout
    needi_yout = 12*lmaxbl;
    needi = std::max(needi, needi_yout);
    // in setup, to call findmin_new
    itempi = nroomi + needi_set;
    needi = std::max(needi, itempi);
    // in setslave, to call setcorner
    if (iunst > 1 || idef_ss > 0) {
        FortranArray1DRef<int> nblg_ref = nblg.ref();
        FortranArray1DRef<int> nbci0_ref = nbci0.ref();
        FortranArray1DRef<int> nbcj0_ref = nbcj0.ref();
        FortranArray1DRef<int> nbck0_ref = nbck0.ref();
        FortranArray1DRef<int> nbcidim_ref = nbcidim.ref();
        FortranArray1DRef<int> nbcjdim_ref = nbcjdim.ref();
        FortranArray1DRef<int> nbckdim_ref = nbckdim.ref();
        FortranArray4DRef<int> ibcinfo_ref = ibcinfo.ref();
        FortranArray4DRef<int> jbcinfo_ref = jbcinfo.ref();
        FortranArray4DRef<int> kbcinfo_ref = kbcinfo.ref();
        int ibctyp2004 = 2004, ibctyp2014 = 2014, ibctyp2024 = 2024;
        int ibctyp2034 = 2034, ibctyp1005 = 1005, ibctyp1006 = 1006, ibctyp2016 = 2016;
        cntsurf_ns::cntsurf(ns2004, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref, nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref, nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp2004);
        cntsurf_ns::cntsurf(ns2014, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref, nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref, nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp2014);
        cntsurf_ns::cntsurf(ns2024, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref, nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref, nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp2024);
        cntsurf_ns::cntsurf(ns2034, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref, nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref, nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp2034);
        cntsurf_ns::cntsurf(ns1005, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref, nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref, nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp1005);
        cntsurf_ns::cntsurf(ns1006, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref, nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref, nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp1006);
        cntsurf_ns::cntsurf(ns2016, maxbl0, maxgr0, maxseg0, ngrid, nblg_ref, nbci0_ref, nbcj0_ref, nbck0_ref, nbcidim_ref, nbcjdim_ref, nbckdim_ref, ibcinfo_ref, jbcinfo_ref, kbcinfo_ref, ibctyp2016);
        nsurf = ns2004 + ns2014 + ns2024 + ns2034 + ns1005 + ns1006 + ns2016;
        itempi = nsurf + needi_set;
        needi = std::max(needi, itempi);
    }
    // in mgblk, to output movie data
    imovie = 0;
    if (std::abs(movie) > 0) {
        imovie = lmaxbl + 3*lnplts;
    }
    // in mgblk, to call dynptch
    if (iunst > 0) {
        itempi = 65*lmaxbl + 3*lintmx + needi_mgblk + imovie;
        needi = std::max(needi, itempi);
    }
    // in mgblk, to call pre_bc
    itempi = 65*lmaxbl + 4*lnsub1 + needi_mgblk + imovie;
    needi = std::max(needi, itempi);
    // in mgblk, to call qout
    itempi = lmaxbl + 3*lnplts + needi_mgblk + imovie;
    needi = std::max(needi, itempi);
    // in mgblk, to call resetg
    itempi = lmaxbl + needi_mgblk + imovie;
    needi = std::max(needi, itempi);
    // in mgblk, to call bc_period
    itempi = 30*llbcprd + needi_mgblk + imovie;
    needi = std::max(needi, itempi);
    // in mgblk, to call bc_embed
    itempi = 18*llbcemb + needi_mgblk + imovie;
    needi = std::max(needi, itempi);
    // in mgblk, to call bc_blkint
    itempi = 30*lmxbli + needi_mgblk + imovie;
    needi = std::max(needi, itempi);
    // in mgblk, to call bc_patch
    itempi = 21*lintmax*lnsub1 + needi_mgblk + imovie;
    needi = std::max(needi, itempi);
    // in qout, to call bc_period
    itempi = 30*llbcprd + needi_qout;
    needi = std::max(needi, itempi);
    // in qout, to call bc_embed
    itempi = 18*llbcemb + needi_qout;
    needi = std::max(needi, itempi);
    // in qout, to call bc_blkint
    itempi = 30*lmxbli + needi_qout;
    needi = std::max(needi, itempi);
    // in qout, to call bc_patch
    itempi = 21*lintmax*lnsub1 + needi_qout;
    needi = std::max(needi, itempi);
    // need to add 1 for safety
    needi = needi + 1;
    // currently, mpi code requires same integer work space as sequential code
    needi_node = needi;

    // Set various work array sizes
    isum   = (int)ftot_seq + need;
    isumi  = std::max(needi, 1);
    // parallel build real and integer work array sizes (per compute node)
    isum_n  = tot_nod;
    isumi_n = std::max(needi_node, 1);

    // Determine number of slave points for deforming mesh case
    ivert = 0;
    if (iunst > 1 || idef_ss > 0) {
        for (nbl = 1; nbl <= nblock; nbl++) {
            if (levelg(nbl) >= lglobal) {
                {
                    FortranArray2DRef<int> lw_ref = lw.ref();
                    FortranArray2DRef<int> lw2_ref = lw2.ref();
                    lead_ns::lead(nbl, lw_ref, lw2_ref, maxbl0);
                }
                if (std::abs(isktyp) == 1) {
                    iskp = iskip(nbl, 1);
                    jskp = jskip(nbl, 1);
                    kskp_v = kskip(nbl, 1);
                    for (i = 1; i <= idim; i += iskp) {
                        for (j = 1; j <= jdim; j += jskp) {
                            for (k = 1; k <= kdim; k += kskp_v) {
                                ivert = ivert + 1;
                            }
                        }
                    }
                } else if (std::abs(isktyp) == 2) {
                    itot = 0; jtot = 0; ktot = 0;
                    for (ii = 1; ii <= 500; ii++) {
                        if (iskip(nbl, ii) == 0) break;
                        itot = itot + 1;
                    }
                    for (ii = 1; ii <= 500; ii++) {
                        if (jskip(nbl, ii) == 0) break;
                        jtot = jtot + 1;
                    }
                    for (ii = 1; ii <= 500; ii++) {
                        if (kskip(nbl, ii) == 0) break;
                        ktot = ktot + 1;
                    }
                    ivert = ivert + itot*jtot*ktot;
                }
            }
        }
    }
    lnslave = std::max(ivert, lnslave);

    // Add any extra memory input by the user
    if (nnodes == 1) {
        isum  = isum  + memadd;
        isumi = isumi + memaddi;
    } else {
        isum_n  = isum_n  + memadd;
        isumi_n = isumi_n + memaddi;
    }

    // Set the auxiliary array sizes
    if (nnodes == 1) {
        mwork  = isum;
        mworki = isumi;
    } else {
        mwork  = isum_n;
        mworki = isumi_n;
    }
    nplots   = lnplts;
    minnode  = lminnode;
    iitot    = liitot;
    intmax   = lintmax;
    maxxe    = lmaxxe;
    mxbli    = lmxbli;
    nsub1    = lnsub1;
    lbcprd   = llbcprd;
    lbcemb   = llbcemb;
    lbcrad   = llbcrad;
    maxbl    = lmaxbl;
    maxgr    = lmaxgr;
    maxseg   = lmxseg;
    maxcs    = lmaxcs;
    ncycmax  = lncycm;
    intmx    = lintmx;
    mxxe     = lmxxe;
    mptch    = lmptch;
    msub1    = lmsub1;
    ibufdim  = libufdim;
    nbuf     = nbuf0;
    nmds     = lnmds;
    maxaes   = lmaxaes;
    nslave   = lnslave;
    nmaster  = lnmaster;
    maxsegdg = lmxsegdg;


    // Print parameter sizes
    fprintf(f66, "\n");
    fprintf(f66, "***********************************************\n");
    fprintf(f66, "\n");
    fprintf(f66, "    PARAMETER SIZES REQUIRED FOR THIS CASE:\n");
    fprintf(f66, "\n");
    fprintf(f66, "***********************************************\n");
    fprintf(f66, "\n");
    if (imode == 0) {
        fprintf(f66, " mwork    = %10d (sequential)\n", isum);
        fprintf(f66, " mworki   = %10d (sequential)\n", isumi);
        fprintf(f66, " mwork    = %10d (parallel)\n", isum_n);
        fprintf(f66, " mworki   = %10d (parallel)\n", isumi_n);
    } else {
        fprintf(f66, " mwork    = %10d\n", mwork);
        fprintf(f66, " mworki   = %10d\n", mworki);
    }
    fprintf(f66, " nplots   = %10d\n", nplots);
    fprintf(f66, " minnode  = %10d\n", minnode);
    fprintf(f66, " iitot    = %10d\n", iitot);
    fprintf(f66, " intmax   = %10d\n", intmax);
    fprintf(f66, " maxxe    = %10d\n", maxxe);
    fprintf(f66, " mxbli    = %10d\n", mxbli);
    fprintf(f66, " nsub1    = %10d\n", nsub1);
    fprintf(f66, " lbcprd   = %10d\n", lbcprd);
    fprintf(f66, " lbcemb   = %10d\n", lbcemb);
    fprintf(f66, " lbcrad   = %10d\n", lbcrad);
    fprintf(f66, " maxbl    = %10d\n", maxbl);
    fprintf(f66, " maxgr    = %10d\n", maxgr);
    fprintf(f66, " maxseg   = %10d\n", maxseg);
    fprintf(f66, " maxcs    = %10d\n", maxcs);
    fprintf(f66, " ncycmax  = %10d\n", ncycmax);
    fprintf(f66, " intmx    = %10d\n", intmx);
    fprintf(f66, " mxxe     = %10d\n", mxxe);
    fprintf(f66, " mptch    = %10d\n", mptch);
    fprintf(f66, " msub1    = %10d\n", msub1);
    fprintf(f66, " nmds     = %10d\n", nmds);
    fprintf(f66, " maxaes   = %10d\n", maxaes);
    fprintf(f66, " ibufdim  = %10d\n", ibufdim);
    fprintf(f66, " nbuf     = %10d\n", nbuf);
    fprintf(f66, " nslave   = %10d\n", nslave);
    fprintf(f66, " nmaster  = %10d\n", nmaster);
    fprintf(f66, " maxsegdg = %10d\n", maxsegdg);
    fprintf(f66, "\n");

    // Stop if any of the parameters are non-zero
    ichek = 1;
    if (mwork <= 0 || mworki <= 0 || nplots <= 0 || minnode <= 0 || iitot <= 0 ||
        intmax <= 0 || maxxe <= 0 || mxbli <= 0 || nsub1 <= 0 || lbcprd <= 0 ||
        lbcemb <= 0 || lbcrad <= 0 || maxbl <= 0 || maxgr <= 0 || maxseg <= 0 ||
        maxcs <= 0 || ncycmax <= 0 || intmx <= 0 || mxxe <= 0 || mptch <= 0 ||
        msub1 <= 0 || ibufdim <= 0 || nbuf <= 0 || nmds <= 0 || maxaes <= 0 ||
        nslave <= 0 || nmaster <= 0 || maxsegdg <= 0) ichek = 0;
    if (ichek == 0) {
        nou(1) = std::min(nou(1)+1, ibufdim0);
        std::memset(bou(nou(1), 1), ' ', 120);
        const char* msg = "error in routine sizer - one of the parameters listed above is zero";
        std::memcpy(bou(nou(1), 1), msg, std::min((int)strlen(msg), 120));
        ierrflg = -99;
        termn8_ns::termn8(myid, ierrflg, ibufdim0, nbuf0, bou, nou);
    }

    // Output storage requirement summary
    fprintf(f66, "\n");
    fprintf(f66, "\n");
    fprintf(f66, "**************************************************************\n");
    fprintf(f66, "\n");
    fprintf(f66, "   SUMMARY OF STORAGE REQUIREMENTS - W + WK ARRAYS\n");
    fprintf(f66, "\n");
    fprintf(f66, " sequential version:\n");
    fprintf(f66, "\n");
    fprintf(f66, "          permanent array w   requires %10d (words)\n", (int)ftot_seq);
    fprintf(f66, "          temporary array wk  requires %10d (words)\n", need);
    fprintf(f66, "          temporary array iwk requires %10d (words)\n", needi);
    fprintf(f66, "\n");
    fprintf(f66, " parallel version, per node:\n");
    fprintf(f66, "\n");
    fprintf(f66, "          permanent array w   requires %10d (words)\n", mem_w_max);
    fprintf(f66, "          temporary array wk  requires %10d (words)\n", mem_wk_max);
    fprintf(f66, "          temporary array iwk requires %10d (words)\n", needi_node);
    fprintf(f66, "  \n");
    fprintf(f66, "  \n");
    fprintf(f66, " >>> Estimate for mwork      (sequential)     = %10d <<<\n", isum);
    fprintf(f66, "\n >>> Estimate for mworki     (sequential)     = %10d <<<\n", isumi);
    fprintf(f66, "\n >>> Estimate for mwork  (per node, parallel) = %10d <<<\n", isum_n);
    fprintf(f66, "\n >>> Estimate for mworki (per node, parallel) = %10d <<<\n", isumi_n);
    fprintf(f66, "\n >>> Parallel code sized for%4d nodes, min. (+host)       <<<\n", nodel);

    // Calculate best possible speedup
    maxpt = 0;
    nptt  = 0;
    for (ind = 1; ind <= nodel; ind++) {
        npt  = no_of_points(ind);
        nptt = nptt + npt;
        if (npt > maxpt) maxpt = npt;
    }

    fprintf(f66, "\n");
    fprintf(f66, "\n");
    fprintf(f66, "**************************************************************\n");
    fprintf(f66, "\n");

    // Output ideal speedup vs no. of nodes
    if (imode == 0) {
        char speedup_file[] = "ideal_speedup.dat";
        fortran_open_unit(67, speedup_file, "w");
        FILE* f67 = fortran_get_unit(67);
        fprintf(f67, " compute nodes  ideal speedup\n");
        for (nnn = 1; nnn <= ngrid; nnn++) {
            {
                FortranArray1DRef<int> ncgg_ref = ncgg.ref();
                FortranArray1DRef<int> nblg_ref = nblg.ref();
                FortranArray1DRef<int> idimg_ref = idimg.ref();
                FortranArray1DRef<int> jdimg_ref = jdimg.ref();
                FortranArray1DRef<int> kdimg_ref = kdimg.ref();
                FortranArray1DRef<int> nblcg_ref = nblcg.ref();
                FortranArray2DRef<int> iwrk_ref = iwrk.ref();
                FortranArray1DRef<int> mblk2nd_ref = mblk2nd.ref();
                compg2n_ns::compg2n(nblock, ngrid, ncgg_ref, nblg_ref, idimg_ref,
                                    jdimg_ref, kdimg_ref, nblcg_ref, nnn, iwrk_ref,
                                    myid, myhost, mblk2nd_ref, mycomm, maxgr0, maxbl0,
                                    ierrflg, ibufdim, nbuf, bou, nou);
            }
            for (ii = 1; ii <= nnn; ii++) {
                mem_req_node(ii) = 0;
                no_of_points(ii) = 0;
            }
            for (i = 1; i <= nblock; i++) {
                npoints = (idimg(i)-1)*(jdimg(i)-1)*(kdimg(i)-1);
                nod = mblk2nd(i);
                mem_req_node(nod) = mem_req_node(nod) + memblock(i);
                no_of_points(nod) = no_of_points(nod) + npoints;
            }
            npts_max  = no_of_points(1);
            mem_w_max = mem_req_node(1);
            for (i = 2; i <= nnn; i++) {
                if (mem_req_node(i) > mem_w_max) mem_w_max = mem_req_node(i);
                if (no_of_points(i) > npts_max)  npts_max  = no_of_points(i);
            }
            maxpt = 0; nptt = 0;
            for (ind = 1; ind <= nnn; ind++) {
                npt  = no_of_points(ind);
                nptt = nptt + npt;
                if (npt > maxpt) maxpt = npt;
            }
            speedopt = (float)nptt / (float)maxpt;
            fprintf(f67, "          %4d        %7.2f\n", nnn, (float)speedopt);
        }
        fprintf(f66, "\n an estimate of ideal parallel speedup has been put in file ideal_speedup.dat\n\n");
    }

    if (nodel < nnodes) {
        fprintf(f66, "\n NOTE: the value of nodes on the command line was overspecified\n for the number of global-level grids\n");
        FILE* f11 = fortran_get_unit(11);
        if (f11) fprintf(f11, "\n NOTE: the value of nodes on the command line was overspecified\n for the number of global-level grids\n");
    }

    // Rewind iunit5
    {
        FILE* fu = fortran_get_unit(iunit5);
        if (fu) rewind(fu);
    }


    // Free up memory used by sizing routine
    ifree = 1;
    if (ifree > 0) {
        ltot.deallocate();
        jjmax1.deallocate();
        kkmax1.deallocate();
        iiint1.deallocate();
        iiint2.deallocate();
        n14.deallocate();
        iwrk.deallocate();
        icsinfo.deallocate();
        iv.deallocate();
        memblock.deallocate();
        lwdat.deallocate();
        nblfine.deallocate();
        mem_req_node.deallocate();
        no_of_points.deallocate();
        lw.deallocate();
        lw2.deallocate();
        nblk.deallocate();
        limblk.deallocate();
        isva.deallocate();
        nblon.deallocate();
        lig.deallocate();
        lbg.deallocate();
        iovrlp.deallocate();
        ibpntsg.deallocate();
        iipntsg.deallocate();
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
        jlamlog.deallocate();
        klamlog.deallocate();
        ilamlog.deallocate();
        jlamhig.deallocate();
        klamhig.deallocate();
        ilamhig.deallocate();
        iwfg.deallocate();
        utrans.deallocate();
        vtrans.deallocate();
        wtrans.deallocate();
        omegax.deallocate();
        omegay.deallocate();
        omegaz.deallocate();
        xorig.deallocate();
        yorig.deallocate();
        zorig.deallocate();
        dxmx.deallocate();
        dymx.deallocate();
        dzmx.deallocate();
        dthxmx.deallocate();
        dthymx.deallocate();
        dthzmx.deallocate();
        thetax.deallocate();
        thetay.deallocate();
        thetaz.deallocate();
        rfreqt.deallocate();
        rfreqr.deallocate();
        xorig0.deallocate();
        yorig0.deallocate();
        zorig0.deallocate();
        time2.deallocate();
        thetaxl.deallocate();
        thetayl.deallocate();
        thetazl.deallocate();
        itrans.deallocate();
        irotat.deallocate();
        idefrm.deallocate();
        bcvali.deallocate();
        bcvalj.deallocate();
        bcvalk.deallocate();
        nbci0.deallocate();
        nbcj0.deallocate();
        nbck0.deallocate();
        nbcidim.deallocate();
        nbcjdim.deallocate();
        nbckdim.deallocate();
        ibcinfo.deallocate();
        jbcinfo.deallocate();
        kbcinfo.deallocate();
        bcfilei.deallocate();
        bcfilej.deallocate();
        bcfilek.deallocate();
        ncgg.deallocate();
        nblg.deallocate();
        iemg.deallocate();
        inewgg.deallocate();
        inpl3d.deallocate();
        inpr.deallocate();
        iadvance.deallocate();
        iforce.deallocate();
        iindex.deallocate();
        iindx.deallocate();
        llimit.deallocate();
        iitmax.deallocate();
        mmcxie.deallocate();
        mmceta.deallocate();
        ncheck.deallocate();
        iifit.deallocate();
        mblkpt.deallocate();
        iic0.deallocate();
        iiorph.deallocate();
        iitoss.deallocate();
        ifiner.deallocate();
        dx.deallocate();
        dy.deallocate();
        dz.deallocate();
        dthetx.deallocate();
        dthety.deallocate();
        dthetz.deallocate();
        isav_blk.deallocate();
        isav_prd.deallocate();
        isav_pat.deallocate();
        isav_pat_b.deallocate();
        isav_dpat.deallocate();
        isav_dpat_b.deallocate();
        isav_emb.deallocate();
        mblk2nd.deallocate();
        mglevg.deallocate();
        nemgl.deallocate();
        ipl3dtmp.deallocate();
        ifrom.deallocate();
        xif1.deallocate();
        etf1.deallocate();
        xif2.deallocate();
        etf2.deallocate();
        utrnsae.deallocate();
        vtrnsae.deallocate();
        wtrnsae.deallocate();
        omgxae.deallocate();
        omgyae.deallocate();
        omgzae.deallocate();
        xorgae.deallocate();
        yorgae.deallocate();
        zorgae.deallocate();
        xorgae0.deallocate();
        yorgae0.deallocate();
        zorgae0.deallocate();
        icouple.deallocate();
        thtxae.deallocate();
        thtyae.deallocate();
        thtzae.deallocate();
        rfrqtae.deallocate();
        rfrqrae.deallocate();
        icsi.deallocate();
        icsf.deallocate();
        jcsi.deallocate();
        jcsf.deallocate();
        kcsi.deallocate();
        kcsf.deallocate();
        idfrmseg.deallocate();
        iaesurf.deallocate();
        nsegdfrm.deallocate();
        freq.deallocate();
        gmass.deallocate();
        x0.deallocate();
        gf0.deallocate();
        damp.deallocate();
        perturb.deallocate();
        aesrfdat.deallocate();
        iskip.deallocate();
        jskip.deallocate();
        kskip.deallocate();
    }

    if (imode == 1) {
        fprintf(f66, "\n memory for precfl3d has been deallocated\n");
    } else {
        printf("\nprecfl3d has completed successfully\n");
        printf("precfl3d information has been put in file precfl3d.out\n\n");
    }

    fortran_close_unit(66);

    return;
} // end sizer


void cntfa(int& nbl, int& nwfa, int& maxbl, int& maxseg,
           FortranArray1DRef<int> nblcg, FortranArray1DRef<int> ieg,
           FortranArray1DRef<int> isg, FortranArray1DRef<int> jdimg,
           FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg,
           FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
           FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcjdim,
           FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbcidim,
           FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
           FortranArray4DRef<int> ibcinfo, int& nblock)
{
    int nblc, nblcc, jfm1, kfm1, ifm1, nsi, nseg, l, i, j, k;
    nwfa = 0;
    for (nblc = 1; nblc <= nblock; nblc++) {
        if (nbl == nblc) continue;
        nblcc = nblcg(nblc);
        if (nblcc == nbl) {
            jfm1 = jdimg(nblc) - 1;
            kfm1 = kdimg(nblc) - 1;
            ifm1 = idimg(nblc) - 1;
            nsi  = ifm1 / (ieg(nblc) - isg(nblc));
            for (nseg = 1; nseg <= nbcj0(nblc); nseg++) {
                if (jbcinfo(nblc, nseg, 1, 1) == 21) {
                    for (l = 1; l <= 5; l++)
                        for (i = 1; i <= ifm1; i += nsi)
                            for (k = 1; k <= kfm1; k += 2)
                                nwfa = nwfa + 1;
                }
            }
            for (nseg = 1; nseg <= nbcjdim(nblc); nseg++) {
                if (jbcinfo(nblc, nseg, 1, 2) == 21) {
                    for (l = 1; l <= 5; l++)
                        for (i = 1; i <= ifm1; i += nsi)
                            for (k = 1; k <= kfm1; k += 2)
                                nwfa = nwfa + 1;
                }
            }
            for (nseg = 1; nseg <= nbck0(nblc); nseg++) {
                if (kbcinfo(nblc, nseg, 1, 1) == 21) {
                    for (l = 1; l <= 5; l++)
                        for (i = 1; i <= ifm1; i += nsi)
                            for (j = 1; j <= jfm1; j += 2)
                                nwfa = nwfa + 1;
                }
            }
            for (nseg = 1; nseg <= nbckdim(nblc); nseg++) {
                if (kbcinfo(nblc, nseg, 1, 2) == 21) {
                    for (l = 1; l <= 5; l++)
                        for (i = 1; i <= ifm1; i += nsi)
                            for (j = 1; j <= jfm1; j += 2)
                                nwfa = nwfa + 1;
                }
            }
            for (nseg = 1; nseg <= nbci0(nblc); nseg++) {
                if (ibcinfo(nblc, nseg, 1, 1) == 21) {
                    for (l = 1; l <= 5; l++)
                        for (k = 1; k <= kfm1; k += nsi)
                            for (j = 1; j <= jfm1; j += 2)
                                nwfa = nwfa + 1;
                }
            }
            for (nseg = 1; nseg <= nbcidim(nblc); nseg++) {
                if (ibcinfo(nblc, nseg, 1, 2) == 21) {
                    for (l = 1; l <= 5; l++)
                        for (k = 1; k <= kfm1; k += nsi)
                            for (j = 1; j <= jfm1; j += 2)
                                nwfa = nwfa + 1;
                }
            }
        }
    }
}

void cntblmx(int& nbl, int& jdim, int& kdim, int& idim, int& iover,
             int& iprint, int& n14, int& maxbl, int& maxseg,
             FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
             FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcjdim,
             FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbcidim,
             FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
             FortranArray4DRef<int> ibcinfo)
{
    int32_t* ivisc = cmn_reyue.ivisc;
    int ilfreq1 = 1, ilfreq2 = 5, ipfreq1 = 10, ipfreq2 = 10;
    int jdim1, kdim1, idim1;
    int ibcjmin, ibcjmax, ibckmin, ibckmax, ibcimin, ibcimax;
    int mseg, ihead, i, j, k, in;
    int kloop, jloop, iloop, inmax, inmax1;

    if (idim == 2 || jdim == 2 || kdim == 2) { ilfreq2 = 1; ipfreq2 = 1; }
    jdim1 = jdim - 1; kdim1 = kdim - 1; idim1 = idim - 1;

    ibcjmin = 0; ibcjmax = 0;
    for (mseg = 1; mseg <= nbcj0(nbl); mseg++) {
        int v = std::abs(jbcinfo(nbl, mseg, 1, 1));
        if (v==2004||v==2014||v==2024||v==2034||v==2016) ibcjmin=1;
    }
    for (mseg = 1; mseg <= nbcjdim(nbl); mseg++) {
        int v = std::abs(jbcinfo(nbl, mseg, 1, 2));
        if (v==2004||v==2014||v==2024||v==2034||v==2016) ibcjmax=1;
    }
    ibckmin = 0; ibckmax = 0;
    for (mseg = 1; mseg <= nbck0(nbl); mseg++) {
        int v = std::abs(kbcinfo(nbl, mseg, 1, 1));
        if (v==2004||v==2014||v==2024||v==2034||v==2016) ibckmin=1;
    }
    for (mseg = 1; mseg <= nbckdim(nbl); mseg++) {
        int v = std::abs(kbcinfo(nbl, mseg, 1, 2));
        if (v==2004||v==2014||v==2024||v==2034||v==2016) ibckmax=1;
    }
    ibcimin = 0; ibcimax = 0;
    for (mseg = 1; mseg <= nbci0(nbl); mseg++) {
        int v = std::abs(ibcinfo(nbl, mseg, 1, 1));
        if (v==2004||v==2014||v==2024||v==2034||v==2016) ibcimin=1;
    }
    for (mseg = 1; mseg <= nbcidim(nbl); mseg++) {
        int v = std::abs(ibcinfo(nbl, mseg, 1, 2));
        if (v==2004||v==2014||v==2024||v==2034||v==2016) ibcimax=1;
    }


    if (ivisc[2] > 1) {
        if (ibckmin == 1 || ibckmax == 0) {
            ihead = 0;
            for (i = 1; i <= idim1; i++) {
                for (j = 1; j <= jdim1; j++) {
                    kloop = (int)(0.80 * kdim);
                    inmax = kloop;
                    if (iover==1 && iprint>=1 && j==j/ilfreq1*ilfreq1 && i==i/ilfreq2*ilfreq2) {
                        if (ihead==0) n14++;
                        n14++; ihead++;
                    }
                    if (iprint>=1 && j==(j/ilfreq1)*ilfreq1 && i==(i/ilfreq2)*ilfreq2) {
                        if (ihead==0) n14++;
                        ihead++; n14++;
                    }
                    if (iprint>1 && j==(j/ipfreq1)*ipfreq1 && i==(i/ipfreq2)*ipfreq2) {
                        ihead=0; n14++;
                        for (in=1; in<=inmax; in++) n14++;
                    }
                }
            }
        }
        if (ibckmax == 1) {
            ihead = 0;
            for (i = 1; i <= idim1; i++) {
                for (j = 1; j <= jdim1; j++) {
                    kloop = (int)(0.80 * kdim);
                    kloop = kdim1 - kloop + 1;
                    kloop = std::max(1, kloop);
                    inmax = kloop; inmax1 = inmax + 1;
                    if (iover==1 && iprint>=1 && j==j/ilfreq1*ilfreq1 && i==i/ilfreq2*ilfreq2) {
                        if (ihead==0) n14++;
                        n14++; ihead++;
                    }
                    if (iprint>=1 && j==(j/ilfreq1)*ilfreq1 && i==(i/ilfreq2)*ilfreq2) {
                        if (ihead==0) n14++;
                        ihead++; n14++;
                    }
                    if (iprint>1 && j==(j/ipfreq1)*ipfreq1 && i==(i/ipfreq2)*ipfreq2) {
                        ihead=0; n14++;
                        for (in=kdim1; in>=inmax1; in--) n14++;
                    }
                }
            }
        }
    }
    if (ivisc[1] > 1) {
        if (ibcjmin == 1 || ibcjmax == 0) {
            ihead = 0;
            for (i = 1; i <= idim1; i++) {
                for (k = 1; k <= kdim1; k++) {
                    jloop = (int)(0.80 * jdim);
                    inmax = jloop;
                    if (iover==1 && iprint>=1 && k==k/ilfreq1*ilfreq1 && i==i/ilfreq2*ilfreq2) {
                        if (ihead==0) n14++;
                        n14++; ihead++;
                    }
                    if (iprint>=1 && k==(k/ilfreq1)*ilfreq1 && i==(i/ilfreq2)*ilfreq2) {
                        if (ihead==0) n14++;
                        ihead++; n14++;
                    }
                    if (iprint>1 && k==(k/ipfreq1)*ipfreq1 && i==(i/ipfreq2)*ipfreq2) {
                        ihead=0; n14++;
                        for (in=1; in<=inmax; in++) n14++;
                    }
                }
            }
        }
        if (ibcjmax == 1) {
            ihead = 0;
            for (i = 1; i <= idim1; i++) {
                for (k = 1; k <= kdim1; k++) {
                    jloop = (int)(0.80 * jdim);
                    jloop = jdim1 - jloop + 1;
                    jloop = std::max(1, jloop);
                    inmax = jloop; inmax1 = inmax + 1;
                    if (iover==1 && iprint>=1 && k==k/ilfreq1*ilfreq1 && i==i/ilfreq2*ilfreq2) {
                        if (ihead==0) n14++;
                        n14++; ihead++;
                    }
                    if (iprint>=1 && k==(k/ilfreq1)*ilfreq1 && i==(i/ilfreq2)*ilfreq2) {
                        if (ihead==0) n14++;
                        ihead++; n14++;
                    }
                    if (iprint>1 && k==(k/ipfreq1)*ipfreq1 && i==(i/ipfreq2)*ipfreq2) {
                        ihead=0; n14++;
                        for (in=jdim1; in>=inmax1; in--) n14++;
                    }
                }
            }
        }
    }
    if (ivisc[0] > 1) {
        if (ibcimin == 1 || ibcimax == 0) {
            ihead = 0;
            for (k = 1; k <= kdim1; k++) {
                for (j = 1; j <= jdim1; j++) {
                    iloop = (int)(0.80 * idim);
                    inmax = iloop;
                    if (iover==1 && iprint>=1 && j==j/ilfreq1*ilfreq1 && k==k/ilfreq2*ilfreq2) {
                        if (ihead==0) n14++;
                        n14++; ihead++;
                    }
                    if (iprint>=1 && j==(j/ilfreq1)*ilfreq1 && k==(k/ilfreq2)*ilfreq2) {
                        if (ihead==0) n14++;
                        ihead++; n14++;
                    }
                    if (iprint>1 && j==(j/ipfreq1)*ipfreq1 && k==(k/ipfreq2)*ipfreq2) {
                        ihead=0; n14++;
                        for (in=1; in<=inmax; in++) n14++;
                    }
                }
            }
        }
        if (ibcimax == 1) {
            ihead = 0;
            for (k = 1; k <= kdim1; k++) {
                for (j = 1; j <= jdim1; j++) {
                    iloop = (int)(0.80 * idim);
                    iloop = idim1 - iloop + 1;
                    iloop = std::max(1, iloop);
                    inmax = iloop; inmax1 = inmax + 1;
                    if (iover==1 && iprint>=1 && j==j/ilfreq1*ilfreq1 && k==k/ilfreq2*ilfreq2) {
                        if (ihead==0) n14++;
                        n14++; ihead++;
                    }
                    if (iprint>=1 && j==(j/ilfreq1)*ilfreq1 && k==(k/ilfreq2)*ilfreq2) {
                        if (ihead==0) n14++;
                        ihead++; n14++;
                    }
                    if (iprint>1 && j==(j/ipfreq1)*ipfreq1 && k==(k/ipfreq2)*ipfreq2) {
                        ihead=0; n14++;
                        for (in=idim1; in>=inmax1; in--) n14++;
                    }
                }
            }
        }
    }
}


} // namespace sizer_ns
