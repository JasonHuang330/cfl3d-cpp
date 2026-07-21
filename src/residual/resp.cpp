// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "resp.h"
#include <cstdio>
#include "runtime/fortran_io.h"
#include "l2norm.h"
#include "l2norm2.h"
#include "termn8.h"
#include "force.h"
#include "csurf.h"
#include "blkmax.h"
#include "prntcp.h"
#include "my_flush.h"
#include <cmath>
#include <algorithm>

namespace resp_ns {

void resp(int& nbl, int& ntime, int& ntt0, int& jdim, int& kdim, int& idim,
          FortranArray2DRef<double> q, FortranArray2DRef<double> sj,
          FortranArray2DRef<double> sk, FortranArray2DRef<double> si,
          FortranArray1DRef<double> vol, FortranArray1DRef<double> x,
          FortranArray1DRef<double> y, FortranArray1DRef<double> z,
          FortranArray2DRef<double> res, FortranArray1DRef<double> wk0,
          FortranArray3DRef<double> vmuk, FortranArray3DRef<double> vmuj,
          FortranArray3DRef<double> vmui, FortranArray1DRef<double> wk,
          int& nwork, int& nblstag, int& nblendg, int& nblstat, int& nblendt,
          int& nptsr, int& nptsrb, FortranArray1DRef<int> nnegb,
          FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
          FortranArray3DRef<double> bci, FortranArray3DRef<double> blank,
          int& nt, FortranArray1DRef<double> sumn, FortranArray1DRef<int> negn,
          double& rmst, FortranArray4DRef<double> xtbj,
          FortranArray4DRef<double> xtbk, FortranArray4DRef<double> xtbi,
          FortranArray4DRef<double> qc0, FortranArray4DRef<double> dqc0,
          int& iseq, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
          FortranArray4DRef<double> qi0, int& maxbl, int& maxgr, int& maxseg,
          FortranArray1DRef<double> rms, double& rmsb,
          FortranArray1DRef<double> rmstb, FortranArray1DRef<double> clw,
          FortranArray1DRef<double> cdw, FortranArray1DRef<double> cdpw,
          FortranArray1DRef<double> cdvw, FortranArray1DRef<double> cxw,
          FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw,
          FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw,
          FortranArray1DRef<double> cmzw, int& n_clcd,
          FortranArray3DRef<double> clcd, int& nblocks_clcd,
          FortranArray2DRef<int> blocks_clcd, FortranArray1DRef<double> chdw,
          FortranArray1DRef<double> swetw, FortranArray1DRef<double> fmdotw,
          FortranArray1DRef<double> cfttotw, FortranArray1DRef<double> cftmomw,
          FortranArray1DRef<double> cftpw, FortranArray1DRef<double> cftvw,
          FortranArray2DRef<double> rmstr, FortranArray2DRef<int> nneg,
          int& ihstry, int& ngrid, FortranArray1DRef<int> nblg,
          FortranArray1DRef<int> iemg, FortranArray1DRef<int> levelg,
          FortranArray2DRef<int> iviscg, FortranArray1DRef<int> itrans,
          FortranArray1DRef<int> irotat, FortranArray1DRef<int> iforce,
          FortranArray1DRef<double> swett, FortranArray1DRef<double> clt,
          FortranArray1DRef<double> cdt, FortranArray1DRef<double> cxt,
          FortranArray1DRef<double> cyt, FortranArray1DRef<double> czt,
          FortranArray1DRef<double> cmxt, FortranArray1DRef<double> cmyt,
          FortranArray1DRef<double> cmzt, FortranArray1DRef<double> cdpt,
          FortranArray1DRef<double> cdvt, FortranArray1DRef<int> nbci0,
          FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
          FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim,
          FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo,
          FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
          FortranArray1DRef<double> resmx, FortranArray1DRef<int> imx,
          FortranArray1DRef<int> jmx, FortranArray1DRef<int> kmx,
          FortranArray1DRef<double> vormax, FortranArray1DRef<int> ivmax,
          FortranArray1DRef<int> jvmax, FortranArray1DRef<int> kvmax,
          FortranArray1DRef<double> sx, FortranArray1DRef<double> sy,
          FortranArray1DRef<double> sz, FortranArray1DRef<double> stot,
          FortranArray1DRef<double> pav, FortranArray1DRef<double> ptav,
          FortranArray1DRef<double> tav, FortranArray1DRef<double> ttav,
          FortranArray1DRef<double> xmav, FortranArray1DRef<double> fmdot,
          FortranArray1DRef<double> cfxp, FortranArray1DRef<double> cfyp,
          FortranArray1DRef<double> cfzp, FortranArray1DRef<double> cfdp,
          FortranArray1DRef<double> cflp, FortranArray1DRef<double> cftp,
          FortranArray1DRef<double> cfxv, FortranArray1DRef<double> cfyv,
          FortranArray1DRef<double> cfzv, FortranArray1DRef<double> cfdv,
          FortranArray1DRef<double> cflv, FortranArray1DRef<double> cftv,
          FortranArray1DRef<double> cfxmom, FortranArray1DRef<double> cfymom,
          FortranArray1DRef<double> cfzmom, FortranArray1DRef<double> cfdmom,
          FortranArray1DRef<double> cflmom, FortranArray1DRef<double> cftmom,
          FortranArray1DRef<double> cfxtot, FortranArray1DRef<double> cfytot,
          FortranArray1DRef<double> cfztot, FortranArray1DRef<double> cfdtot,
          FortranArray1DRef<double> cfltot, FortranArray1DRef<double> cfttot,
          int& ncs, FortranArray2DRef<int> icsinfo, int& myid, int& myhost,
          int& mycomm, FortranArray1DRef<int> mblk2nd,
          FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
          int& nbuf, int& ibufdim, int& ncycmax, int& maxcs,
          FortranArray1DRef<double> thetay, FortranArray1DRef<int> iadvance,
          FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg,
          int& nummem)
{
    // COMMON block references
    int& i2d        = cmn_twod.i2d;
    int& lq2avg     = cmn_ginfo2.lq2avg;
    int& iskip_blocks = cmn_ginfo2.iskip_blocks;
    float& xnumavg  = cmn_avgdata.xnumavg;
    int& iteravg    = cmn_avgdata.iteravg;
    float& xnumavg2 = cmn_avgdata.xnumavg2;
    int& ipertavg   = cmn_avgdata.ipertavg;
    int& iclcd      = cmn_avgdata.iclcd;
    int& isubit_r   = cmn_avgdata.isubit_r;
    int& iaccnt     = cmn_account.iaccnt;
    int& ioutsub    = cmn_account.ioutsub;
    float& sref     = cmn_fsum.sref;
    float& cref     = cmn_fsum.cref;
    float& bref     = cmn_fsum.bref;
    float& xmc      = cmn_fsum.xmc;
    float& ymc      = cmn_fsum.ymc;
    float& zmc      = cmn_fsum.zmc;
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
    int& levt       = cmn_mgrd.levt;
    int& kode       = cmn_mgrd.kode;
    int& mode       = cmn_mgrd.mode;
    int& ncyc       = cmn_mgrd.ncyc;
    int& mtt        = cmn_mgrd.mtt;
    int& icyc       = cmn_mgrd.icyc;
    int& level      = cmn_mgrd.level;
    int& lglobal    = cmn_mgrd.lglobal;
    float& time     = cmn_unst.time;
    float& cfltau   = cmn_unst.cfltau;
    int& ntstep     = cmn_unst.ntstep;
    int& ita        = cmn_unst.ita;
    int& iunst      = cmn_unst.iunst;
    float& cfltau0  = cmn_unst.cfltau0;
    float& cfltauMax = cmn_unst.cfltaumax;
    int& nwrest     = cmn_wrbl.nwrest;
    float& cdv      = cmn_drag.cdv;
    float& cdp      = cmn_drag.cdp;
    float& radtodeg = cmn_conversion.radtodeg;
    float& cltsub   = cmn_subit.cltsub;
    float& cdtsub   = cmn_subit.cdtsub;
    float& cxtsub   = cmn_subit.cxtsub;
    float& cytsub   = cmn_subit.cytsub;
    float& cztsub   = cmn_subit.cztsub;
    float& cmxtsub  = cmn_subit.cmxtsub;
    float& cmytsub  = cmn_subit.cmytsub;
    float& cmztsub  = cmn_subit.cmztsub;
    float& cdptsub  = cmn_subit.cdptsub;
    float& cdvtsub  = cmn_subit.cdvtsub;
    float& sxsub    = cmn_subit.sxsub;
    float& sysub    = cmn_subit.sysub;
    float& szsub    = cmn_subit.szsub;
    float& stotsub  = cmn_subit.stotsub;
    float& fmdotsub = cmn_subit.fmdotsub;
    float& cfxpsub  = cmn_subit.cfxpsub;
    float& cfypsub  = cmn_subit.cfypsub;
    float& cfzpsub  = cmn_subit.cfzpsub;
    float& cflpsub  = cmn_subit.cflpsub;
    float& cfdpsub  = cmn_subit.cfdpsub;
    float& cftpsub  = cmn_subit.cftpsub;
    float& cfxvsub  = cmn_subit.cfxvsub;
    float& cfyvsub  = cmn_subit.cfyvsub;
    float& cfzvsub  = cmn_subit.cfzvsub;
    float& cflvsub  = cmn_subit.cflvsub;
    float& cfdvsub  = cmn_subit.cfdvsub;
    float& cftvsub  = cmn_subit.cftvsub;
    float& cfxmomsub = cmn_subit.cfxmomsub;
    float& cfymomsub = cmn_subit.cfymomsub;
    float& cfzmomsub = cmn_subit.cfzmomsub;
    float& cflmomsub = cmn_subit.cflmomsub;
    float& cfdmomsub = cmn_subit.cfdmomsub;
    float& cftmomsub = cmn_subit.cftmomsub;
    float& cfxtotsub = cmn_subit.cfxtotsub;
    float& cfytotsub = cmn_subit.cfytotsub;
    float& cfztotsub = cmn_subit.cfztotsub;
    float& cfdtotsub = cmn_subit.cfdtotsub;
    float& cfltotsub = cmn_subit.cfltotsub;
    float& cfttotsub = cmn_subit.cfttotsub;
    int& ip3dgrd    = cmn_igrdtyp.ip3dgrd;
    int& ialph      = cmn_igrdtyp.ialph;

    // Local variables
    int jdim1, kdim1, idim1;
    int iflag, ifor, ifo, jfo, kfo, icall, iuns;
    int l, ics, igrid, iskip, nbl1, iem, ivar;
    int jm, km, im, nb, nn, nout, iwrit3, irite, iwk1, iwk2, iwk3;
    double resd, rmssum, t1, rmstx, rmsts1, fnorm, alot;
    double cl, cd, czz, cyy, cxx, cmy, cmx, cmz, chd, swet;
    double yref, zref, xref;
    double sxi, syi, szi, stoti, pavi, ptavi, tavi, ttavi, xmavi, fmdoti;
    double cfxpi, cfypi, cfzpi, cfdpi, cflpi, cftpi;
    double cfxvi, cfyvi, cfzvi, cfdvi, cflvi, cftvi;
    double cfxmomi, cfymomi, cfzmomi, cfdmomi, cflmomi, cftmomi;
    double cfxtoti, cfytoti, cfztoti, cfdtoti, cfltoti, cfttoti;
    double resmax;
    double resa[6];   // 1-based, indices 1..5
    double rmsba[6];  // 1-based, indices 1..5



    // Build multi-dimensional views for array arguments
    // res is passed as (jdim*kdim*(idim-1), 5) → reshape to 4D for l2norm
    FortranArray4DRef<double> res4d(res.data(), jdim, kdim, idim-1, 5);
    // vol is passed as 1D (jdim*kdim*(idim-1)) → reshape to 3D
    FortranArray3DRef<double> vol3d(vol.data(), jdim, kdim, idim-1);
    // x,y,z are 1D (jdim*kdim*idim) → reshape to 3D
    FortranArray3DRef<double> x3d(x.data(), jdim, kdim, idim);
    FortranArray3DRef<double> y3d(y.data(), jdim, kdim, idim);
    FortranArray3DRef<double> z3d(z.data(), jdim, kdim, idim);
    // q is (jdim*kdim*idim, 5) → reshape to 4D
    FortranArray4DRef<double> q4d(q.data(), jdim, kdim, idim, 5);
    // sk,sj,si are (jdim*kdim*(idim-1), 5) → reshape to 4D
    FortranArray4DRef<double> sk4d(sk.data(), jdim, kdim, idim-1, 5);
    FortranArray4DRef<double> sj4d(sj.data(), jdim, kdim, idim-1, 5);
    FortranArray4DRef<double> si4d(si.data(), jdim, kdim, idim, 5);
    // q(1,2), q(1,3), q(1,4) → 3D views of columns 2,3,4
    FortranArray3DRef<double> q_ub(&q(1,2), jdim, kdim, idim);
    FortranArray3DRef<double> q_vb(&q(1,3), jdim, kdim, idim);
    FortranArray3DRef<double> q_wb(&q(1,4), jdim, kdim, idim);

    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;
    // global iteration convergence
    nres = ntt;
    iflag = 0;

    if ((float)dt < 0.f && nt == 1) {
        iflag = 1;
        resd = 0.;
        if (myid == mblk2nd(nbl) && iadvance(nbl) >= 0) {
            int irdq0 = 0;
            l2norm_ns::l2norm(nbl, icyc, resd, irdq0, jdim, kdim, idim, res4d, vol3d);
            if (((float)(resd+resd) == (float)(resd) &&
                 std::abs((float)(resd)) > 1.e-20f) ||
                !((float)(resd) < std::abs((float)(resd)) ||
                  (float)(resd) >= 0.e0f)) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                snprintf(bou(nou(1),1), 120,
                    " NaN detected after residual evaluation, block%4d cycle%5d",
                    nbl, icyc);
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
        }
    }

    if ((float)dt > 0.f && icyc == ioutsub) {
        iflag = 1;
        // l2-norm of density residual EXCLUDING time terms
        resd = 0.;
        if (myid == mblk2nd(nbl) && iadvance(nbl) >= 0) {
            int irdq0 = 0;
            l2norm2_ns::l2norm2(nbl, icyc, resd, irdq0, jdim, kdim, idim,
                                res4d, vol3d, qc0, dqc0, q4d, blank);
            if (((float)(resd+resd) == (float)(resd) &&
                 std::abs((float)(resd)) > 1.e-20f) ||
                !((float)(resd) < std::abs((float)(resd)) ||
                  (float)(resd) >= 0.e0f)) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                snprintf(bou(nou(1),1), 120,
                    "NaN detected after residual evaluation, block%4d time step%5d",
                    nbl, nt);
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
        }
    }
    if (iflag == 1) {
        // density residual
        if (myid == myhost && iadvance(nbl) >= 0) {
            rmssum = 0.;
            if (nptsr > 0)
                rmssum = rms(ntt) * rms(ntt) * (float)nptsr;
            t1    = (float)((jdim-1)*(kdim-1)*(idim-1));
            rmst  = std::sqrt(resd/t1);
            nptsr = nptsr + (jdim-1)*(kdim-1)*(idim-1);
            rmssum       = rmssum + resd;
            rms(ntt)    = std::sqrt(rmssum/(float)nptsr);
        }

        // turbulence residual
        if ((iviscg(nbl,1) > 2 || iviscg(nbl,2) > 2
             || iviscg(nbl,3) > 2) && iadvance(nbl) >= 0) {
            if (myid == myhost) {
                for (l = 1; l <= nummem; l++) {
                    rmstx = 0.;
                    if (nptsr > 0) {
                        nptsr = nptsr - jdim1*kdim1*idim1;
                        rmstx = rmstr(nres,l)*rmstr(nres,l)*(float)nptsr;
                    }
                    nptsr        = nptsr + jdim1*kdim1*idim1;
                    rmstx        = rmstx
                                 + sumn(l)*sumn(l)*(float)(jdim1*kdim1*idim1);
                    rmstr(nres,l) = std::sqrt(rmstx / (float)nptsr);
                    nneg(nres,l)  = nneg(nres,l) + negn(l);
                }
            }
        }

        // forces and moments
        if (iforce(nbl) > 0) {
            ifor  = iforce(nbl);
            ifo   = ifor/100;
            jfo   = (ifor - ifo*100)/10;
            kfo   = (ifor - ifo*100 - jfo*10);
            icall = 0;
            iuns  = std::max({itrans(nbl), irotat(nbl), idefrm(nbl)});

            if (myid == mblk2nd(nbl)) {
                force_ns::force(jdim, kdim, idim, x3d, y3d, z3d,
                                sk4d, sj4d, si4d,
                                cl, cd, czz, cyy, cxx, cmy, cmx, cmz,
                                chd, swet, icall,
                                q_ub, q_vb, q_wb,
                                vmuk, vmuj, vmui, vol3d,
                                ifo, jfo, kfo, bcj, bck, bci, blank, nbl,
                                xtbj, xtbk, xtbi, iuns,
                                qj0, qk0, qi0,
                                nbci0, nbcj0, nbck0,
                                nbcidim, nbcjdim, nbckdim,
                                ibcinfo, jbcinfo, kbcinfo, maxbl, maxseg);
            }


            if (myid == mblk2nd(nbl)) {
                if (ialph == 0) {
                    yref = cref;
                    zref = bref;
                    xref = bref;
                } else {
                    yref = bref;
                    zref = cref;
                    xref = bref;
                }
                clt(nbl)   = cl/sref;
                if(!std::isfinite((double)clt(nbl))){static int nb=0;nb++;if(nb<=3)std::fprintf(stderr,"R_CLT nbl=%d cl=%g sref=%g clt=%g\n",nbl,(double)cl,(double)sref,(double)clt(nbl));}
                cmyt(nbl)  = cmy/(sref*yref);
                cmxt(nbl)  = cmx/(sref*xref);
                cmzt(nbl)  = cmz/(sref*zref);
                cdt(nbl)   = cd/sref;
                czt(nbl)   = czz/sref;
                cyt(nbl)   = cyy/sref;
                cxt(nbl)   = cxx/sref;
                cdpt(nbl)  = cdp/sref;
                cdvt(nbl)  = cdv/sref;
                swett(nbl) = swet;
            }

            // sum contributions of blocks on global level
            if (myid == myhost) {
                if (level == lglobal) {
                    chdw(nres)  = chdw(nres)  + chd;
                    swetw(nres) = swetw(nres) + swett(nbl);
                    clw(nres)   = clw(nres)   + clt(nbl);
                    cmyw(nres)  = cmyw(nres)  + cmyt(nbl);
                    cdw(nres)   = cdw(nres)   + cdt(nbl);
                    cmxw(nres)  = cmxw(nres)  + cmxt(nbl);
                    cmzw(nres)  = cmzw(nres)  + cmzt(nbl);
                    czw(nres)   = czw(nres)   + czt(nbl);
                    cyw(nres)   = cyw(nres)   + cyt(nbl);
                    cxw(nres)   = cxw(nres)   + cxt(nbl);
                    cdpw(nres)  = cdpw(nres)  + cdpt(nbl);
                    cdvw(nres)  = cdvw(nres)  + cdvt(nbl);

                    if (iclcd == 1 || iclcd == 2) {
                        for (nb = 1; nb <= nblocks_clcd; nb++) {
                            if (igridg(nbl) == blocks_clcd(2,nb)) {
                                nn = blocks_clcd(1,nb);
                                clcd(1,nn,nres) = clcd(1,nn,nres) + clt(nbl);
                                clcd(2,nn,nres) = clcd(2,nn,nres) + cdt(nbl);
                            }
                        }
                    }
                }
            }
        } // end if (iforce(nbl) > 0)

        // control surface data
        for (ics = 1; ics <= ncs; ics++) {
            if (nbl == icsinfo(ics,1) && levelg(nbl) >= lglobal) {
                iuns = std::max({itrans(nbl), irotat(nbl), idefrm(nbl)});
                if (myid == mblk2nd(nbl)) {
                    csurf_ns::csurf(jdim, kdim, idim, x3d, y3d, z3d,
                                    sk4d, sj4d, si4d, q4d, ics,
                                    q_ub, q_vb, q_wb,
                                    vmuk, vmuj, vmui, vol3d,
                                    bcj, bck, bci, blank,
                                    xtbj, xtbk, xtbi, iuns, ncs, icsinfo,
                                    sxi, syi, szi, stoti,
                                    pavi, ptavi, tavi, ttavi, xmavi, fmdoti,
                                    cfxpi, cfypi, cfzpi, cfdpi, cflpi, cftpi,
                                    cfxvi, cfyvi, cfzvi, cfdvi, cflvi, cftvi,
                                    cfxmomi, cfymomi, cfzmomi, cfdmomi, cflmomi,
                                    cftmomi, cfxtoti, cfytoti, cfztoti, cfdtoti,
                                    cfltoti, cfttoti, maxcs, qj0, qk0, qi0);
                    sx(ics)     = sxi;
                    sy(ics)     = syi;
                    sz(ics)     = szi;
                    stot(ics)   = stoti;
                    pav(ics)    = pavi;
                    ptav(ics)   = ptavi;
                    tav(ics)    = tavi;
                    ttav(ics)   = ttavi;
                    xmav(ics)   = xmavi;
                    fmdot(ics)  = fmdoti;
                    cfxp(ics)   = cfxpi;
                    cfyp(ics)   = cfypi;
                    cfzp(ics)   = cfzpi;
                    cflp(ics)   = cflpi;
                    cfdp(ics)   = cfdpi;
                    cftp(ics)   = cftpi;
                    cfxv(ics)   = cfxvi;
                    cfyv(ics)   = cfyvi;
                    cfzv(ics)   = cfzvi;
                    cflv(ics)   = cflvi;
                    cfdv(ics)   = cfdvi;
                    cftv(ics)   = cftvi;
                    cfxmom(ics) = cfxmomi;
                    cfymom(ics) = cfymomi;
                    cfzmom(ics) = cfzmomi;
                    cflmom(ics) = cflmomi;
                    cfdmom(ics) = cfdmomi;
                    cftmom(ics) = cftmomi;
                    cfxtot(ics) = cfxtoti;
                    cfytot(ics) = cfytoti;
                    cfztot(ics) = cfztoti;
                    cfdtot(ics) = cfdtoti;
                    cfltot(ics) = cfltoti;
                    cfttot(ics) = cfttoti;
                }
            }
        }



        // sum contributions over all control surfaces with inorm .ne. 0
        // (global level only)
        if (myid == myhost) {
            if (nbl == nblendg && levelg(nbl) == lglobal) {
                sx(ncs+1)     = 0.e0;
                sy(ncs+1)     = 0.e0;
                sz(ncs+1)     = 0.e0;
                stot(ncs+1)   = 0.e0;
                fmdot(ncs+1)  = 0.e0;
                cfxp(ncs+1)   = 0.e0;
                cfyp(ncs+1)   = 0.e0;
                cfzp(ncs+1)   = 0.e0;
                cflp(ncs+1)   = 0.e0;
                cfdp(ncs+1)   = 0.e0;
                cfxv(ncs+1)   = 0.e0;
                cfyv(ncs+1)   = 0.e0;
                cfzv(ncs+1)   = 0.e0;
                cflv(ncs+1)   = 0.e0;
                cfdv(ncs+1)   = 0.e0;
                cfxmom(ncs+1) = 0.e0;
                cfymom(ncs+1) = 0.e0;
                cfzmom(ncs+1) = 0.e0;
                cflmom(ncs+1) = 0.e0;
                cfdmom(ncs+1) = 0.e0;
                cfxtot(ncs+1) = 0.e0;
                cfytot(ncs+1) = 0.e0;
                cfztot(ncs+1) = 0.e0;
                cfdtot(ncs+1) = 0.e0;
                cfltot(ncs+1) = 0.e0;
                for (ics = 1; ics <= ncs; ics++) {
                    fnorm = (float)(icsinfo(ics,9));
                    int nbl1_cs = icsinfo(ics,1);
                    if ((float)fnorm != 0.e0f &&
                        levelg(nbl1_cs) == lglobal) {
                        sx(ncs+1)     = sx(ncs+1)     + sx(ics);
                        sy(ncs+1)     = sy(ncs+1)     + sy(ics);
                        sz(ncs+1)     = sz(ncs+1)     + sz(ics);
                        stot(ncs+1)   = stot(ncs+1)   + stot(ics);
                        fmdot(ncs+1)  = fmdot(ncs+1)  + fmdot(ics);
                        cfxp(ncs+1)   = cfxp(ncs+1)   + cfxp(ics);
                        cfyp(ncs+1)   = cfyp(ncs+1)   + cfyp(ics);
                        cfzp(ncs+1)   = cfzp(ncs+1)   + cfzp(ics);
                        cflp(ncs+1)   = cflp(ncs+1)   + cflp(ics);
                        cfdp(ncs+1)   = cfdp(ncs+1)   + cfdp(ics);
                        cfxv(ncs+1)   = cfxv(ncs+1)   + cfxv(ics);
                        cfyv(ncs+1)   = cfyv(ncs+1)   + cfyv(ics);
                        cfzv(ncs+1)   = cfzv(ncs+1)   + cfzv(ics);
                        cflv(ncs+1)   = cflv(ncs+1)   + cflv(ics);
                        cfdv(ncs+1)   = cfdv(ncs+1)   + cfdv(ics);
                        cfxmom(ncs+1) = cfxmom(ncs+1) + cfxmom(ics);
                        cfymom(ncs+1) = cfymom(ncs+1) + cfymom(ics);
                        cfzmom(ncs+1) = cfzmom(ncs+1) + cfzmom(ics);
                        cflmom(ncs+1) = cflmom(ncs+1) + cflmom(ics);
                        cfdmom(ncs+1) = cfdmom(ncs+1) + cfdmom(ics);
                        cfxtot(ncs+1) = cfxtot(ncs+1) + cfxtot(ics);
                        cfytot(ncs+1) = cfytot(ncs+1) + cfytot(ics);
                        cfztot(ncs+1) = cfztot(ncs+1) + cfztot(ics);
                        cfltot(ncs+1) = cfltot(ncs+1) + cfltot(ics);
                        cfdtot(ncs+1) = cfdtot(ncs+1) + cfdtot(ics);

                        cftp(ncs+1)   = std::sqrt(cfxp(ncs+1)*cfxp(ncs+1)
                                                 + cfyp(ncs+1)*cfyp(ncs+1)
                                                 + cfzp(ncs+1)*cfzp(ncs+1));
                        cftv(ncs+1)   = std::sqrt(cfxv(ncs+1)*cfxv(ncs+1)
                                                 + cfyv(ncs+1)*cfyv(ncs+1)
                                                 + cfzv(ncs+1)*cfzv(ncs+1));
                        cftmom(ncs+1) = std::sqrt(cfxmom(ncs+1)*cfxmom(ncs+1)
                                                 + cfymom(ncs+1)*cfymom(ncs+1)
                                                 + cfzmom(ncs+1)*cfzmom(ncs+1));
                        cfttot(ncs+1) = std::sqrt(cfxtot(ncs+1)*cfxtot(ncs+1)
                                                 + cfytot(ncs+1)*cfytot(ncs+1)
                                                 + cfztot(ncs+1)*cfztot(ncs+1));

                        // store selected totals for convergence history
                        fmdotw(nres)  = fmdot(ncs+1);
                        cftmomw(nres) = cftmom(ncs+1);
                        cftpw(nres)   = cftp(ncs+1);
                        cftvw(nres)   = cftv(ncs+1);
                        cfttotw(nres) = cfttot(ncs+1);
                    }
                }
            }
        }

    } // end if (iflag == 1)



    // print residual and force coefficients
    if (myid == myhost) {
        if ((float)dt < 0.f && nt == 1) {
            iflag = 0;
            if (icyc == 1) iflag = 1;
            if (icyc == 2 && nbl == nblendg) iflag = 1;
            if (icyc == ncyc && nbl == nblstat)  iflag = 1;
            if (icyc/nwrest*nwrest == icyc && nbl == nblstat)
                iflag = 1;
            if ((icyc-1)/nwrest*nwrest == (icyc-1) &&
                nbl == nblstag) iflag = 1;
            if (iflag == 1) {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, "\n level    blk   iter   residual   total res.     lift        drag     side force\n");
            }
            iflag = 0;
            if (icyc == 1) iflag = 1;
            if (icyc >= 2 && nbl == nblendg) iflag = 1;
            if (icyc == ncyc && nbl >= nblstag) iflag = 1;
            if (icyc/nwrest*nwrest == icyc && nbl >= nblstag)
                iflag = 1;
            if (iflag > 0) {
                FILE* f11 = fortran_get_unit(11);
                if (ialph != 0) {
                    fprintf(f11, "     %1d %6d %6d %11.4e %11.4e %11.4e %11.4e %11.4e\n",
                            level, nbl, nres, (float)rmst,
                            (float)rms(nres), (float)clw(nres),
                            (float)cdw(nres), (float)czw(nres));
                } else {
                    fprintf(f11, "     %1d %6d %6d %11.4e %11.4e %11.4e %11.4e %11.4e\n",
                            level, nbl, nres, (float)rmst,
                            (float)rms(nres), (float)clw(nres),
                            (float)cdw(nres), (float)cyw(nres));
                }
            }
        }

        if ((float)dt > 0.f && icyc == ioutsub) {
            iflag = 0;
            if (ioutsub == 1) {
                if (nt == 1) iflag = 1;
                if (nt == 2 && nbl == nblendg) iflag = 1;
                if (nt == ntstep && nbl == nblstat) iflag = 1;
                if (nt/nwrest*nwrest == nt && nbl == nblstat)
                    iflag = 1;
                if (nt > 1 && (nt-1)/nwrest*nwrest == (nt-1) &&
                    nbl == nblstag) iflag = 1;
            }
            if (ioutsub == ncyc) {
                if (nt == 1 && nbl == nblstat) iflag = 1;
                if (nt == 2 && nbl == nblendg) iflag = 1;
                if (nt == ntstep && nbl == nblstat) iflag = 1;
                if (nt/nwrest*nwrest == nt && nbl == nblstat)
                    iflag = 1;
                if (nt > 1 && (nt-1)/nwrest*nwrest == (nt-1) &&
                    nbl == nblstag) iflag = 1;
            }
            if (iflag > 0) {
                FILE* f11 = fortran_get_unit(11);
                if (iunst == 0) {
                    fprintf(f11, "\n level    blk  ntime   residual   total res.     lift        drag     side force\n");
                } else {
                    fprintf(f11, "\n level    blk  ntime   residual   total res.     lift        drag     side force     time       alpha\n");
                }
            }
            iflag = 0;
            if (ioutsub == 1) {
                if (nt == 1) iflag = 1;
                if (nt >= 2 && nbl == nblendg) iflag = 1;
                if (nt == ntstep && nbl >= nblstag) iflag = 1;
                if (nt/nwrest*nwrest == nt && nbl >= nblstag)
                    iflag = 1;
            }
            if (ioutsub == ncyc) {
                if (nt == 1 && nbl >= nblstag) iflag = 1;
                if (nt >= 2 && nbl == nblendg) iflag = 1;
                if (nt == ntstep && nbl >= nblstag) iflag = 1;
                if (nt/nwrest*nwrest == nt && nbl >= nblstag)
                    iflag = 1;
            }
            if (iflag > 0) {
                FILE* f11 = fortran_get_unit(11);
                if (iunst == 0) {
                    if (ialph != 0) {
                        fprintf(f11, "     %1d %6d %6d %11.4e %11.4e %11.4e %11.4e %11.4e\n",
                                level, nbl, nres, (float)rmst,
                                (float)rms(nres), (float)clw(nres),
                                (float)cdw(nres), (float)czw(nres));
                    } else {
                        fprintf(f11, "     %1d %6d %6d %11.4e %11.4e %11.4e %11.4e %11.4e\n",
                                level, nbl, nres, (float)rmst,
                                (float)rms(nres), (float)clw(nres),
                                (float)cdw(nres), (float)cyw(nres));
                    }
                } else {
                    alot = ((double)alpha + (double)thetay(nbl)) * (double)radtodeg;
                    if (ialph != 0) {
                        fprintf(f11, "     %1d %6d %6d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",
                                level, nbl, nres, (float)rmst,
                                (float)rms(nres), (float)clw(nres),
                                (float)cdw(nres), (float)czw(nres),
                                (float)time, (float)alot);
                    } else {
                        fprintf(f11, "     %1d %6d %6d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",
                                level, nbl, nres, (float)rmst,
                                (float)rms(nres), (float)clw(nres),
                                (float)cdw(nres), (float)cyw(nres),
                                (float)time, (float)alot);
                    }
                }
            }
        }
    } // end if (myid == myhost) for print section



    // print maximum residual/location and, if turbulent,
    // maximum vorticity/location - print whenever restart
    // file is written and at end of computation
    iflag = 0;
    if ((float)dt <= 0.f && nt == 1) {
        if (icyc/nwrest*nwrest == icyc || icyc == ncyc) iflag = 1;
    }
    if ((float)dt > 0.f && icyc == ioutsub) {
        if (nt/nwrest*nwrest == nt || nt == ntstep) iflag = 1;
    }

    if (iflag > 0) {
        if (myid == mblk2nd(nbl)) {
            blkmax_ns::blkmax(jdim, kdim, idim, res4d, resmax, jm, km, im);
            resmx(nbl) = resmax;
            imx(nbl)   = im;
            jmx(nbl)   = jm;
            kmx(nbl)   = km;
        }
    }

    if (myid == myhost) {
        iflag = 0;
        if ((float)dt <= 0.f && nt == 1) {
            if (level == lglobal && nbl == nblendg) {
                if (icyc/nwrest*nwrest == icyc || icyc == ncyc) iflag = 1;
            }
        }
        if ((float)dt > 0.f && icyc == ioutsub) {
            if (level == lglobal && nbl == nblendg) {
                if (nt/nwrest*nwrest == nt || nt == ntstep) iflag = 1;
            }
        }

        if (iflag > 0) {
            FILE* f11 = fortran_get_unit(11);
            for (igrid = 1; igrid <= ngrid; igrid++) {
                iskip = 0;
                if (igrid == 1) iskip = 1;
                nbl1 = nblg(igrid);
                iem  = iemg(igrid);
                if (iseq != mseq && iem > 0) continue;
                if (iem == 0) nbl1 = nbl1 + (mseq - iseq);
                if (iskip > 0) {
                    fprintf(f11,
                        "\n max residual on block  %6d (grid %6d) is %12.5e at j,k,i = %4d,%4d,%4d.\n",
                        nbl1, igrid, (float)resmx(nbl1),
                        jmx(nbl1), kmx(nbl1), imx(nbl1));
                } else {
                    fprintf(f11,
                        " max residual on block  %6d (grid %6d) is %12.5e at j,k,i = %4d,%4d,%4d.\n",
                        nbl1, igrid, (float)resmx(nbl1),
                        jmx(nbl1), kmx(nbl1), imx(nbl1));
                }
            }
            for (igrid = 1; igrid <= ngrid; igrid++) {
                iskip = 0;
                if (igrid == 1) iskip = 1;
                nbl1 = nblg(igrid);
                iem  = iemg(igrid);
                if (iseq != mseq && iem > 0) continue;
                if (iem == 0) nbl1 = nbl1 + (mseq - iseq);
                if ((iviscg(nbl1,3) > 1 || iviscg(nbl1,2) > 1
                     || iviscg(nbl1,1) > 1) && iadvance(nbl) >= 0) {
                    if (iskip > 0) {
                        fprintf(f11,
                            "\n max vorticity on block %6d (grid %6d) is %12.5e at j,k,i = %4d,%4d,%4d.\n",
                            nbl1, igrid, (float)vormax(nbl1),
                            jvmax(nbl1), kvmax(nbl1), ivmax(nbl1));
                    } else {
                        fprintf(f11,
                            " max vorticity on block %6d (grid %6d) is %12.5e at j,k,i = %4d,%4d,%4d.\n",
                            nbl1, igrid, (float)vormax(nbl1),
                            jvmax(nbl1), kvmax(nbl1), ivmax(nbl1));
                    }
                }
            }
        }
    }



    // output unsteady cp data for dynamic mesh cases
    // irite = 0 don't output data
    //         1 output data with old headers
    //         2 output data with new headers
    irite = 0;
    iwk1  = 1;
    iwk2  = iwk1 + jdim*kdim;
    iwk3  = iwk2 + kdim*idim;
    if (irite > 0 && iunst > 0) {
        FortranArray1DRef<double> wk_iwk1(&wk(iwk1), jdim*kdim);
        FortranArray1DRef<double> wk_iwk2(&wk(iwk2), kdim*idim);
        FortranArray1DRef<double> wk_iwk3(&wk(iwk3), jdim*idim);
        prntcp_ns::prntcp(jdim, kdim, idim,
                          wk_iwk1, wk_iwk2, wk_iwk3, q4d,
                          nbl, maxbl, maxseg,
                          ibcinfo, jbcinfo, kbcinfo,
                          nbci0, nbcj0, nbck0,
                          nbcidim, nbcjdim, nbckdim,
                          thetay, mblk2nd, myid, myhost, mycomm, irite);
    }

    // sub-iteration convergence - monitor density residual,
    // force/moment coefficients (and turbulence residuals,
    // if appropriate), and control surface data
    if ((float)dt >= 0.f && ncyc > 1) {

        // l2-norm of full residual (INCLUDING time terms)
        resd = 0.;
        if (myid == mblk2nd(nbl) && iadvance(nbl) >= 0) {
            int irdq0 = 0;
            l2norm_ns::l2norm(nbl, icyc, resd, irdq0, jdim, kdim, idim, res4d, vol3d);
        }

        if (myid == myhost && iadvance(nbl) >= 0) {
            rmssum = 0.;
            if (nptsrb > 0) rmssum = rmsb*rmsb*(float)nptsrb;
            t1     = (float)(1*jdim1*kdim1*idim1);
            // rmstb  = sqrt(resd/t1)
            nptsrb = nptsrb + jdim1*kdim1*idim1;
            rmssum = rmssum + resd;
            rmsb   = std::sqrt(rmssum / (float)nptsrb);
        }

        if (isubit_r != 0) {
            // l2-norm of full residual (INCLUDING time terms)
            for (int ii = 1; ii <= 5; ii++) resa[ii] = 0.;
            FortranArray1DRef<double> resa_ref(&resa[1], 5);
            if (myid == mblk2nd(nbl) && iadvance(nbl) >= 0) {
                int irdq0 = 0;
                l2norm_ns::l2normall(nbl, icyc, resa_ref, irdq0, jdim, kdim, idim, res4d, vol3d);
            }

            if (myid == myhost && iadvance(nbl) >= 0) {
                for (ivar = 1; ivar <= 5; ivar++) {
                    rmssum = 0.;
                    if (nptsrb > jdim1*kdim1*idim1) {
                        rmssum = rmsba[ivar]*rmsba[ivar]
                               * (float)(nptsrb - jdim1*kdim1*idim1);
                    }
                    rmssum = rmssum + resa[ivar];
                    rmsba[ivar] = std::sqrt(rmssum / (float)nptsrb);
                }
            }
        }

        // turbulence residual
        if ((iviscg(nbl,1) > 2 || iviscg(nbl,2) > 2 ||
             iviscg(nbl,3) > 2) && iadvance(nbl) >= 0) {
            if (myid == myhost) {
                for (l = 1; l <= nummem; l++) {
                    rmsts1 = 0.;
                    if (nptsrb > 0) {
                        nptsrb = nptsrb - jdim1*kdim1*idim1;
                        // if (nptsrb.gt.0) then
                        rmsts1 = rmstb(l)*rmstb(l)*(float)nptsrb;
                        // end if
                    }
                    nptsrb = nptsrb + jdim1*kdim1*idim1;
                    rmsts1 = rmsts1
                           + sumn(l)*sumn(l)*(float)(jdim1*kdim1*idim1);
                    rmstb(l) = std::sqrt(rmsts1 / (float)nptsrb);
                    nnegb(l) = nnegb(l) + negn(l);
                }
            }
        }



        // compute forces and moments (global level only)
        if (level == lglobal) {
            if (myid == myhost) {
                // zero out values if this is the first block on the global level
                if (nbl == nblstag) {
                    cltsub   = 0.f;
                    cmytsub  = 0.f;
                    cmxtsub  = 0.f;
                    cmztsub  = 0.f;
                    cdtsub   = 0.f;
                    cztsub   = 0.f;
                    cytsub   = 0.f;
                    cxtsub   = 0.f;
                    cdptsub  = 0.f;
                    cdvtsub  = 0.f;
                }
            }

            if (iforce(nbl) > 0) {
                ifor  = iforce(nbl);
                ifo   = ifor/100;
                jfo   = (ifor - ifo*100)/10;
                kfo   = (ifor - ifo*100 - jfo*10);
                icall = 0;
                iuns  = std::max({itrans(nbl), irotat(nbl), idefrm(nbl)});
                if (myid == mblk2nd(nbl)) {
                    force_ns::force(jdim, kdim, idim, x3d, y3d, z3d,
                                    sk4d, sj4d, si4d,
                                    cl, cd, czz, cyy, cxx, cmy, cmx, cmz,
                                    chd, swet, icall,
                                    q_ub, q_vb, q_wb,
                                    vmuk, vmuj, vmui, vol3d,
                                    ifo, jfo, kfo, bcj, bck, bci, blank, nbl,
                                    xtbj, xtbk, xtbi, iuns,
                                    qj0, qk0, qi0,
                                    nbci0, nbcj0, nbck0,
                                    nbcidim, nbcjdim, nbckdim,
                                    ibcinfo, jbcinfo, kbcinfo, maxbl, maxseg);
                }
                if (myid == myhost) {
                    // sum contributions of blocks on global level
                    if (ialph == 0) {
                        yref = cref;
                        zref = bref;
                        xref = bref;
                    } else {
                        yref = bref;
                        zref = cref;
                        xref = bref;
                    }
                    cltsub   = cltsub   + (float)(cl/sref);
                    cmytsub  = cmytsub  + (float)(cmy/(sref*yref));
                    cmxtsub  = cmxtsub  + (float)(cmx/(sref*xref));
                    cmztsub  = cmztsub  + (float)(cmz/(sref*zref));
                    cdtsub   = cdtsub   + (float)(cd/sref);
                    cztsub   = cztsub   + (float)(czz/sref);
                    cytsub   = cytsub   + (float)(cyy/sref);
                    cxtsub   = cxtsub   + (float)(cxx/sref);
                    cdptsub  = cdptsub  + (float)(cdp/sref);
                    cdvtsub  = cdvtsub  + (float)(cdv/sref);
                }
            }
        } // end if (level == lglobal) for forces

        // control surface data (global level only)
        if (level == lglobal) {
            if (myid == myhost) {
                // zero out values if this is the first block on the global level
                if (nbl == nblstag) {
                    sxsub     = 0.e0f;
                    sysub     = 0.e0f;
                    szsub     = 0.e0f;
                    stotsub   = 0.e0f;
                    fmdotsub  = 0.e0f;
                    cfxpsub   = 0.e0f;
                    cfypsub   = 0.e0f;
                    cfzpsub   = 0.e0f;
                    cflpsub   = 0.e0f;
                    cfdpsub   = 0.e0f;
                    cfxvsub   = 0.e0f;
                    cfyvsub   = 0.e0f;
                    cfzvsub   = 0.e0f;
                    cflvsub   = 0.e0f;
                    cfdvsub   = 0.e0f;
                    cfxmomsub = 0.e0f;
                    cfymomsub = 0.e0f;
                    cfzmomsub = 0.e0f;
                    cflmomsub = 0.e0f;
                    cfdmomsub = 0.e0f;
                    cfxtotsub = 0.e0f;
                    cfytotsub = 0.e0f;
                    cfztotsub = 0.e0f;
                    cfdtotsub = 0.e0f;
                    cfltotsub = 0.e0f;
                }
            }

            for (ics = 1; ics <= ncs; ics++) {
                fnorm = (float)(icsinfo(ics,9));
                if (nbl == icsinfo(ics,1) &&
                    (float)fnorm != 0.e0f) {
                    iuns = std::max({itrans(nbl), irotat(nbl), idefrm(nbl)});
                    if (myid == mblk2nd(nbl)) {
                        csurf_ns::csurf(jdim, kdim, idim, x3d, y3d, z3d,
                                        sk4d, sj4d, si4d, q4d, ics,
                                        q_ub, q_vb, q_wb,
                                        vmuk, vmuj, vmui, vol3d,
                                        bcj, bck, bci, blank,
                                        xtbj, xtbk, xtbi, iuns, ncs, icsinfo,
                                        sxi, syi, szi, stoti,
                                        pavi, ptavi, tavi, ttavi, xmavi, fmdoti,
                                        cfxpi, cfypi, cfzpi, cfdpi, cflpi, cftpi,
                                        cfxvi, cfyvi, cfzvi, cfdvi, cflvi, cftvi,
                                        cfxmomi, cfymomi, cfzmomi, cfdmomi, cflmomi,
                                        cftmomi, cfxtoti, cfytoti, cfztoti, cfdtoti,
                                        cfltoti, cfttoti, maxcs, qj0, qk0, qi0);
                    }
                    if (myid == myhost) {
                        sxsub     = sxsub     + (float)sxi;
                        sysub     = sysub     + (float)syi;
                        szsub     = szsub     + (float)szi;
                        stotsub   = stotsub   + (float)stoti;
                        fmdotsub  = fmdotsub  + (float)fmdoti;
                        cfxpsub   = cfxpsub   + (float)cfxpi;
                        cfypsub   = cfypsub   + (float)cfypi;
                        cfzpsub   = cfzpsub   + (float)cfzpi;
                        cflpsub   = cflpsub   + (float)cflpi;
                        cfdpsub   = cfdpsub   + (float)cfdpi;
                        cfxvsub   = cfxvsub   + (float)cfxvi;
                        cfyvsub   = cfyvsub   + (float)cfyvi;
                        cfzvsub   = cfzvsub   + (float)cfzvi;
                        cflvsub   = cflvsub   + (float)cflvi;
                        cfdvsub   = cfdvsub   + (float)cfdvi;
                        cfxmomsub = cfxmomsub + (float)cfxmomi;
                        cfymomsub = cfymomsub + (float)cfymomi;
                        cfzmomsub = cfzmomsub + (float)cfzmomi;
                        cflmomsub = cflmomsub + (float)cflmomi;
                        cfdmomsub = cfdmomsub + (float)cfdmomi;
                        cfxtotsub = cfxtotsub + (float)cfxtoti;
                        cfytotsub = cfytotsub + (float)cfytoti;
                        cfztotsub = cfztotsub + (float)cfztoti;
                        cfdtotsub = cfdtotsub + (float)cfdtoti;
                        cfltotsub = cfltotsub + (float)cfltoti;

                        cftpsub   = std::sqrt((double)cfxpsub*(double)cfxpsub
                                            + (double)cfypsub*(double)cfypsub
                                            + (double)cfzpsub*(double)cfzpsub);
                        cftvsub   = std::sqrt((double)cfxvsub*(double)cfxvsub
                                            + (double)cfyvsub*(double)cfyvsub
                                            + (double)cfzvsub*(double)cfzvsub);
                        cftmomsub = std::sqrt((double)cfxmomsub*(double)cfxmomsub
                                            + (double)cfymomsub*(double)cfymomsub
                                            + (double)cfzmomsub*(double)cfzmomsub);
                        cfttotsub = std::sqrt((double)cfxtotsub*(double)cfxtotsub
                                            + (double)cfytotsub*(double)cfytotsub
                                            + (double)cfztotsub*(double)cfztotsub);
                    }
                }
            }
        } // end if (level == lglobal) for control surfaces



        // output subiteration convergence if iwrit3 > 0
        // iwrit3 > 0...density residual and forces (or cont. surf. data) only
        // iwrit3 > 1...density residual and forces (or cont. surf. data) AND
        //              turbulence residual(s) (field eqn. turb. models only)
        iwrit3 = 2;

        if (myid == myhost) {
            if (ncyc > 1 && nbl == nblendg) {
                nout = (nt-1)*ncyc + icyc;
                if (iwrit3 > 0) {
                    FILE* f23 = fortran_get_unit(23);
                    if (ihstry == 0) {
                        if (ialph == 0) {
                            if (nout == 1) {
                                fprintf(f23, "    subit   log(subres)       cl"
                                             "            cd            cy"
                                             "           cmy\n");
                            }
                            fprintf(f23, "   %6d%14.5e%14.5e%14.5e%14.5e%14.5e\n",
                                    nout, std::log10((float)rmsb),
                                    (float)cltsub, (float)cdtsub,
                                    (float)cytsub, (float)cmytsub);
                        } else {
                            if (nout == 1) {
                                fprintf(f23, "    subit   log(subres)       cl"
                                             "            cd            cz"
                                             "           cmz\n");
                            }
                            fprintf(f23, "   %6d%14.5e%14.5e%14.5e%14.5e%14.5e\n",
                                    nout, std::log10((float)rmsb),
                                    (float)cltsub, (float)cdtsub,
                                    (float)cztsub, (float)cmztsub);
                        }
                        if (isubit_r != 0) {
                            FILE* f111 = fortran_get_unit(111);
                            fprintf(f111, "%6d %14.5e%14.5e%14.5e%14.5e%14.5e\n",
                                    nout,
                                    std::log10((float)rmsba[1]),
                                    std::log10((float)rmsba[2]),
                                    std::log10((float)rmsba[3]),
                                    std::log10((float)rmsba[4]),
                                    std::log10((float)rmsba[5]));
                        }
                    } else if (ihstry == 1) {
                        if (nout == 1) {
                            fprintf(f23, "    subit   log(subres)   mass_flow"
                                         "        cftp      cftv    cftmom\n");
                        }
                        fprintf(f23, "   %6d%14.5e%14.5e%14.5e%14.5e%14.5e\n",
                                nout, std::log10((float)rmsb),
                                (float)fmdotsub, (float)cftpsub,
                                (float)cftvsub, (float)cftmomsub);
                    } else {
                        if (nout == 1) {
                            fprintf(f23, "    subit   log(subres)       cl"
                                         "            cd            cx"
                                         "            cy            cz"
                                         "           cmx           cmy"
                                         "           cmz\n");
                        }
                        fprintf(f23, "   %6d%14.5e%14.5e%14.5e%14.5e%14.5e%14.5e%14.5e%14.5e%14.5e\n",
                                nout, std::log10((float)rmsb),
                                (float)cltsub, (float)cdtsub,
                                (float)cxtsub, (float)cytsub,
                                (float)cztsub, (float)cmxtsub,
                                (float)cmytsub, (float)cmztsub);
                        if (isubit_r != 0) {
                            FILE* f111 = fortran_get_unit(111);
                            fprintf(f111, "%6d %14.5e%14.5e%14.5e%14.5e%14.5e\n",
                                    nout,
                                    std::log10((float)rmsba[1]),
                                    std::log10((float)rmsba[2]),
                                    std::log10((float)rmsba[3]),
                                    std::log10((float)rmsba[4]),
                                    std::log10((float)rmsba[5]));
                        }
                    }
                }
                if (iwrit3 > 1) {
                    if (iviscg(nbl,1) > 2 || iviscg(nbl,2) > 2
                        || iviscg(nbl,3) > 2) {
                        FILE* f24 = fortran_get_unit(24);
                        if (nummem == 2) {
                            if (nout == 1) {
                                fprintf(f24, "    subit  log(turres1)  log(turres2)"
                                             "  nneg1  nneg2\n");
                            }
                            fprintf(f24, "   %6d%14.5e%14.5e %6d %6d\n",
                                    nout,
                                    std::log10((float)rmstb(1)),
                                    std::log10((float)rmstb(2)),
                                    nnegb(1), nnegb(2));
                        } else if (nummem == 3) {
                            if (nout == 1) {
                                fprintf(f24, "    subit  log(turres1)  log(turres2)"
                                             "  log(turres3)  nneg1  nneg2  nneg3\n");
                            }
                            fprintf(f24, "   %6d%14.5e%14.5e%14.5e %6d %6d %6d\n",
                                    nout,
                                    std::log10((float)rmstb(1)),
                                    std::log10((float)rmstb(2)),
                                    std::log10((float)rmstb(3)),
                                    nnegb(1), nnegb(2), nnegb(3));
                        } else {
                            if (nout == 1) {
                                fprintf(f24, "    subit  log(turres1)  log(turres2)"
                                             "  log(turres3)  log(turres4)  log(turres5)"
                                             "  log(turres6)  log(turres7)  nneg1"
                                             "  nneg2  nneg3  nneg4  nneg5  nneg6  nneg7\n");
                            }
                            fprintf(f24, "   %6d%14.5e%14.5e%14.5e%14.5e%14.5e%14.5e%14.5e"
                                         " %6d %6d %6d %6d %6d %6d %6d\n",
                                    nout,
                                    std::log10((float)rmstb(1)),
                                    std::log10((float)rmstb(2)),
                                    std::log10((float)rmstb(3)),
                                    std::log10((float)rmstb(4)),
                                    std::log10((float)rmstb(5)),
                                    std::log10((float)rmstb(6)),
                                    std::log10((float)rmstb(7)),
                                    nnegb(1), nnegb(2), nnegb(3), nnegb(4),
                                    nnegb(5), nnegb(6), nnegb(7));
                        }
                    }
                }
            }
        }

    } // end if ((float)dt >= 0.f && ncyc > 1)

    if (myid == myhost) {
        int iunit11 = 11;
        my_flush_ns::my_flush(iunit11);
    }
    if ((float)dt >= 0.f) {
        if (myid == myhost) {
            int iunit23 = 23;
            my_flush_ns::my_flush(iunit23);
        }
        if (myid == myhost) {
            int iunit24 = 24;
            my_flush_ns::my_flush(iunit24);
        }
    }

} // end resp

} // namespace resp_ns
