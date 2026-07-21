// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "twoeqn.h"
#include "runtime/fortran_io.h"
#include "ccomplex.h"
#include "triv.h"
#include "sijrate2d.h"
#include "sijrate3d.h"
#include "u_doubleprime.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace twoeqn_ns {

using namespace ccomplex_ns;
using namespace triv_ns;
using namespace sijrate2d_ns;
using namespace sijrate3d_ns;
using namespace u_doubleprime_ns;

// Local by-value wrappers so rvalue expressions can be passed to cc* and triv functions
static inline double ccmin(double a1, double a2) { return ccomplex_ns::ccmin(a1, a2); }
static inline double ccmax(double a1, double a2) { return ccomplex_ns::ccmax(a1, a2); }
static inline double ccmincr(double a1, double a2) { return ccomplex_ns::ccmincr(a1, a2); }
static inline double ccmaxcr(double a1, double a2) { return ccomplex_ns::ccmaxcr(a1, a2); }
static inline double ccmaxrc(double a1, double a2) { return ccomplex_ns::ccmaxrc(a1, a2); }
static inline double ccminrc(double a1, double a2) { return ccomplex_ns::ccminrc(a1, a2); }
static inline double ccabs(double a) { return ccomplex_ns::ccabs(a); }
static inline double ccsign(double a, double b) { return ccomplex_ns::ccsign(a, b); }
static inline double ccsignrc(double a, double b) { return ccomplex_ns::ccsignrc(a, b); }
static inline double cctanh(double a) { return ccomplex_ns::cctanh(a); }
static inline double ccacos(double z) { return ccomplex_ns::ccacos(z); }
static inline double cclog10(double z) { return ccomplex_ns::cclog10(z); }
static inline double ccdim(double a1, double a2) { return ccomplex_ns::ccdim(a1, a2); }
static inline void triv(int jdim_, int kdim_, int jl, int ju, int kl, int ku,
    FortranArray2DRef<double> x, FortranArray2DRef<double> a, FortranArray2DRef<double> b,
    FortranArray2DRef<double> c, FortranArray2DRef<double> f) {
    triv_ns::triv(jdim_, kdim_, jl, ju, kl, ku, x, a, b, c, f);
}

void twoeqn(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, FortranArray3DRef<double> dtj, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vor, FortranArray3DRef<double> smin, FortranArray4DRef<double> zksav, FortranArray4DRef<double> turre, FortranArray3DRef<double> damp1, FortranArray3DRef<double> blend, FortranArray3DRef<double> timestp, FortranArray3DRef<double> fnu, FortranArray2DRef<double> bx, FortranArray2DRef<double> bx2, FortranArray2DRef<double> cx, FortranArray2DRef<double> cx2, FortranArray2DRef<double> dx, FortranArray2DRef<double> dx2, FortranArray2DRef<double> fx, FortranArray2DRef<double> fx2, FortranArray2DRef<double> workx, FortranArray2DRef<double> by, FortranArray2DRef<double> by2, FortranArray2DRef<double> cy, FortranArray2DRef<double> cy2, FortranArray2DRef<double> dy, FortranArray2DRef<double> dy2, FortranArray2DRef<double> fy, FortranArray2DRef<double> fy2, FortranArray2DRef<double> worky, FortranArray2DRef<double> bz, FortranArray2DRef<double> bz2, FortranArray2DRef<double> cz, FortranArray2DRef<double> cz2, FortranArray2DRef<double> dz, FortranArray2DRef<double> dz2, FortranArray2DRef<double> fz, FortranArray2DRef<double> fz2, FortranArray2DRef<double> workz, int& ntime, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, int& nbl, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, FortranArray3DRef<double> blank, int& iover, double& sumn1, double& sumn2, int& negn1, int& negn2, FortranArray4DRef<double> ux, FortranArray4DRef<double> rhside, FortranArray4DRef<double> zksav2, FortranArray3DRef<double> v3dtmp, FortranArray3DRef<double> cmuv, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& maxbl, int& maxseg, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iex, int& iex2, int& iex3, FortranArray1DRef<double> dkdx, FortranArray1DRef<double> dkdy, FortranArray1DRef<double> dkdz, FortranArray1DRef<double> dzdx, FortranArray1DRef<double> dzdy, FortranArray1DRef<double> dzdz, FortranArray1DRef<double> dw1dx, FortranArray1DRef<double> dw1dy, FortranArray1DRef<double> dw1dz, FortranArray1DRef<double> dw2dx, FortranArray1DRef<double> dw2dy, FortranArray1DRef<double> dw2dz, FortranArray1DRef<double> dw3dx, FortranArray1DRef<double> dw3dy, FortranArray1DRef<double> dw3dz, FortranArray1DRef<double> drdx, FortranArray1DRef<double> drdy, FortranArray1DRef<double> drdz, FortranArray1DRef<double> dpdx, FortranArray1DRef<double> dpdy, FortranArray1DRef<double> dpdz, FortranArray1DRef<double> dkdj, FortranArray1DRef<double> dkdk, FortranArray1DRef<double> dkdi, FortranArray1DRef<double> dzdj, FortranArray1DRef<double> dzdk, FortranArray1DRef<double> dzdi, FortranArray1DRef<double> dw1dj, FortranArray1DRef<double> dw1dk, FortranArray1DRef<double> dw1di, FortranArray1DRef<double> dw2dj, FortranArray1DRef<double> dw2dk, FortranArray1DRef<double> dw2di, FortranArray1DRef<double> dw3dj, FortranArray1DRef<double> dw3dk, FortranArray1DRef<double> dw3di, FortranArray1DRef<double> drdj, FortranArray1DRef<double> drdk, FortranArray1DRef<double> drdi, FortranArray1DRef<double> dpdj, FortranArray1DRef<double> dpdk, FortranArray1DRef<double> dpdi, FortranArray3DRef<double> w1, FortranArray3DRef<double> w2, FortranArray3DRef<double> w3, FortranArray4DRef<double> srce, FortranArray4DRef<double> vx, FortranArray3DRef<double> xlscale, FortranArray3DRef<double> fdsav, int& nummem, int& iccnum)
{
    // COMMON block aliases
    float& cdes = cmn_des.cdes;
    int& ides = cmn_des.ides;
    float& cddes = cmn_des.cddes;
    float& xmach = cmn_info.xmach;
    float& dt = cmn_info.dt;
    int* ncyc1 = cmn_info.ncyc1;
    int& iconsf = cmn_info.iconsf;
    float* tur10 = cmn_ivals.tur10;
    float& c10 = cmn_easmv.c10;
    float& c11 = cmn_easmv.c11;
    float& c2 = cmn_easmv.c2;
    float& c3 = cmn_easmv.c3;
    float& c4 = cmn_easmv.c4;
    float& c5 = cmn_easmv.c5;
    int& ieasm_type = cmn_easmv.ieasm_type;
    float& gamma = cmn_fluid.gamma;
    float& cbar = cmn_fluid2.cbar;
    int& ilamlo = cmn_lam.ilamlo;
    int& ilamhi = cmn_lam.ilamhi;
    int& jlamlo = cmn_lam.jlamlo;
    int& jlamhi = cmn_lam.jlamhi;
    int& klamlo = cmn_lam.klamlo;
    int& klamhi = cmn_lam.klamhi;
    int& i_lam_forcezero = cmn_lam.i_lam_forcezero;
    int& ivmx = cmn_maxiv.ivmx;
    int& ncyc = cmn_mgrd.ncyc;
    int& icyc = cmn_mgrd.icyc;
    float& reue = cmn_reyue.reue;
    float& tinf = cmn_reyue.tinf;
    int& isklton = cmn_sklton.isklton;
    int& i2d = cmn_twod.i2d;
    int& iexp = cmn_zero.iexp;
    int* iwf = cmn_wallfun.iwf;
    float* cflturb = cmn_turbconv.cflturb;
    float& edvislim = cmn_turbconv.edvislim;
    int& iturbprod = cmn_turbconv.iturbprod;
    int& nsubturb = cmn_turbconv.nsubturb;
    int& nfreeze = cmn_turbconv.nfreeze;
    int& itaturb = cmn_turbconv.itaturb;
    float& tur1cut = cmn_turbconv.tur1cut;
    float& tur2cut = cmn_turbconv.tur2cut;
    int& iturbord = cmn_turbconv.iturbord;
    float& tur1cutlev = cmn_turbconv.tur1cutlev;
    float& tur2cutlev = cmn_turbconv.tur2cutlev;
    int& ita = cmn_unst.ita;
    int& itaturb_unst = cmn_unst.ita; // same as ita
    int& ikoprod = cmn_konew.ikoprod;
    int& isstdenom = cmn_konew.isstdenom;
    float& pklimterm = cmn_konew.pklimterm;
    int& ibeta8kzeta = cmn_konew.ibeta8kzeta;
    int& i_bsl = cmn_konew.i_bsl;
    int& keepambient = cmn_konew.keepambient;
    int& i_wilcox06 = cmn_konew.i_wilcox06;
    int& i_wilcox06_chiw = cmn_konew.i_wilcox06_chiw;
    int& i_turbprod_kterm = cmn_konew.i_turbprod_kterm;
    int& i_catris_kw = cmn_konew.i_catris_kw;
    float& prod2d3dtrace = cmn_konew.prod2d3dtrace;
    int& i_compress_correct = cmn_konew.i_compress_correct;
    int& isstsf = cmn_konew.isstsf;
    int& i_wilcox98 = cmn_konew.i_wilcox98;
    int& i_wilcox98_chiw = cmn_konew.i_wilcox98_chiw;
    int& isst2003 = cmn_konew.isst2003;
    int& isarc2d = cmn_curvat.isarc2d;
    float& sarccr3 = cmn_curvat.sarccr3;
    int& ieasmcc2d = cmn_curvat.ieasmcc2d;
    int& isstrc = cmn_curvat.isstrc;
    float& sstrc_crc = cmn_curvat.sstrc_crc;
    int& isar = cmn_curvat.isar;
    float& crot = cmn_curvat.crot;
    int& isarc3d = cmn_curvat.isarc3d;
    float& cmulim = cmn_easmlim.cmulim;
    int& iaxi2planeturb = cmn_axisym.iaxi2planeturb;
    int& istrongturbdis = cmn_axisym.istrongturbdis;

    // Local scalar variables
    double xminn, catris_kw, turb_mach_0, xsi_star, gam_comp;
    int nsubit;
    double factor, factor1, factor2;
    double vk, a1, cmuc1, cmuc2;
    double beta1, beta2, sigo1, sigo2, sigk1, sigk2, alp1, alp2, sigkmu;
    double sigd0, phi, re, c2b, c2bp;
    int jd2, iwrite;
    double css, gg;
    int idurbinlim;
    double c4new;
    double al10, al1, al2, al3, al4, pi, tpi3;
    double zeta1_kkl, zeta2_kkl, zeta3_kkl, c11_kkl, c12_kkl, cd1_kkl;
    double sumno, sumnk;
    int negno, negnk;
    double xa, ya, za, xp, yp, zp, xm, ym, zm, xc, yc, zc, tc;
    double ttp, ttm, ttpo, ttmo, ttpn, ttmn;
    double anutp, anutm, fnup, fnum, cdp, cdm;
    double rhop, rhom;
    double byy, cyy, dyy, bxx, cxx, dxx, bzz, czz, dzz;
    double sigkp, sigkm, sigop, sigom, dfacep, dfacem;
    double uu, sgnu, app, apm;
    double volku, volkl, volju, voljl, voliu, volil;
    double s11, s22, s33, s12, s13, s23, w12, w13, w23, xis, wis;
    double eta, squig, cmu, cmuj, cmub, cmuc_var, alfa1, alfa2, alfa3;
    double eta1, eta2, eta11, ap, ar, aq, aa, ab, ad, raat, rbbt, raa, rbb;
    double coss, theta, facy;
    double ret, rek, fmu, denom, denom1, denom2, arg2, f2;
    double tracepart, s11t, s22t, s33t, xxx;
    double tau, taulim, fff, eta1_girimaji;
    double omegatemp, zktemp, utau, tauw, dudy;
    int nss;
    double constt;
    double tt, ca;
    double cutoff;
    double pk, dk, dk1, dk2;
    double t11, t22, t33, t12, t13, t23;
    double alpa1, alpa2;
    double factre;
    double fbeta;
    double cyadd, cy2add;
    double betax;
    double term1, term2;
    double uprime, uprimeprime;
    double zlvk, xmin_kkl, xmax_kkl, c_phi1, f_phi, pw, dw;
    // k-enstrophy (ivmx=15) variables
    double cmu_ke, rcmu, alpha3, beta4, beta5, beta6, beta7, beta8;
    double sigr, c1_ke, czeta1, ck_ke, sigp, sigrho, fnuef, delta_ke;
    double wt1, wt2, wlim, slim;
    double sj4, sj4p, sk4, sk4p, si4, si4p;
    double xjp, yjp, zjp, xkp, ykp, zkp, xip, yip, zip;
    double rvol;
    double fnut, divv, omemag;
    double p3, p4p6, p5, p7, p8;
    double wiwjtauij;
    double ct1, ct2, ct3, ct4, ct5;
    double DpDt, termn, termd;
    double xnu, rt, ft, d5;
    double rtaurho, rtauk;
    double pklim;
    // loop indices
    int i, j, k, il, iu, jl, ju, kl, ku;
    int not_iter;
    int kk, jj, ii;
    int kstop, jstop, istop;
    int kset, jset_bc, iset;
    int ibeg, iend, jbeg, jend, kbeg, kend;
    int ip, im, kp, km;
    // misc
    double sigd;
    double chiw, fw;
    double mt2, mt02, f_mt;
    double fr1, fr2, fr3;
    double r_val;
    double sstrc_r;
    double sstsf_r;
    double clim_w06;
    double chiw98, fw98;
    double lturb, ldes;
    double deltaj, deltak, deltai, delta, ell, dist, velterm, rd, fd, term;
    double ti_val, tk_val;
    double srce1, srce2;
    double uu_wf;
    double qset;
    int jset_out, nnumb;
    double ypl, uplus, zkplus, eplus, uvplus;
    int iwriteaux;
    double dkdomega_dot;
    double xsi_w06;
    double betastar_w06;
    double sij_sij, wij_wij;
    double ds11, ds22, ds33, ds12, ds13, ds23;
    double dw12, dw13, dw23;
    double dsdt11, dsdt22, dsdt33, dsdt12, dsdt13, dsdt23;
    double dw12dt, dw13dt, dw23dt;
    double sij_wjk_ski, wij_sjk_ski;
    double sij_sjk_ski;
    double sij_wjk_wki, wij_wjk_ski;
    double sij_sjk_wki;
    double sstrc_fr1;
    double sstsf_fr;
    double alp1_sst2003;
    double chiw98_val;
    double fw98_val;
    double sstrc_r_val;
    double sstsf_r_val;
    double mt2_val;
    double uprimeprime_val;
    double vx_curv;
    double lturb_val, ldes_val;
    double ip_f, im_f, kp_f, km_f;
    double pklimtermset;


    // isklton output messages - part 1
    if (isklton > 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     Computing turbulent viscosity using 2-eqns, block=%5d", nbl);
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "     Freestream tur10,tur20 = %19.8E%19.8E", (float)tur10[0], (float)tur10[1]);
        if (iturbord == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     1st order advection on RHS");
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     2nd order advection on RHS");
        }
    }
    if (isklton > 0) {
        if (ivmx == 6) {
            if (i_wilcox06 == 1) {
                if (i_wilcox06_chiw == 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "     k-omega (Wilcox 06), with vortex stretching term");
                } else {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "     k-omega (Wilcox 06), without vortex stretching term");
                }
            } else if (i_wilcox98 == 1) {
                if (i_wilcox98_chiw == 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "     k-omega (Wilcox 98), with vortex stretching term");
                } else {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "     k-omega (Wilcox 98), without vortex stretching term");
                }
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     k-omega (Wilcox 88)");
            }
            if (ikoprod == 1) {
                if (i_turbprod_kterm != 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "     strain-based production term");
                } else {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "     strain-based production term with 2/3rho*k subtracted");
                }
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     approx (vort) production term");
            }
        }
        if (ivmx == 7 && isst2003 == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-omega SST-2003 (Menter)");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "       ignores ikoprod,isstdenom (these will ALWAYS be on)");
        }
        if (isst2003 != 1) {
            if (ivmx == 7 && i_bsl != 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     k-omega SST (Menter)");
                if (ikoprod == 1) {
                    if (i_turbprod_kterm != 1) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "     strain-based production term");
                    } else {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "     strain-based production term with 2/3rho*k subtracted");
                    }
                } else {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "     approx (vort) production term");
                }
                if (isstdenom == 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "     vort in denom of mut term");
                } else {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "     S in denom of mut term");
                }
            } else if (ivmx == 7 && i_bsl == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     k-omega BSL (Menter)");
                if (ikoprod == 1) {
                    if (i_turbprod_kterm != 1) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "     strain-based production term");
                    } else {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "     strain-based production term with 2/3rho*k subtracted");
                    }
                } else {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "     approx (vort) production term");
                }
            }
        }
        if (ivmx == 8) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-omega (EASM, linear)");
            if (iturbprod == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     strain-based production term");
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     exact production term");
            }
            if (ieasmcc2d == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     curvature correction ON");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "       WARNING: curv terms are active in 2-D sense only!!!");
            }
        }
        if (ivmx == 9) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-epsilon (EASM,linear)");
            if (iturbprod == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     strain-based production term");
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     exact production term");
            }
            if (ieasmcc2d == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     curvature correction ON");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "       WARNING: curv terms are active in 2-D sense only!!!");
            }
        }
    } // end isklton > 0 part 1

    // isklton output messages - part 2
    if (isklton > 0) {
        if (ivmx == 10) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-epsilon (Abid)");
            if (ikoprod == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     strain-based production term");
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     approx (vort) production term");
            }
        }
        if (ivmx == 11) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-epsilon (G-S EASM, nonlinear)");
            if (ieasmcc2d == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     curvature correction ON");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "       WARNING: curv terms are active in 2-D sense only!!!");
            }
        }
        if (ivmx == 12) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-omega (G-S EASM, nonlinear)");
            if (ieasmcc2d == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     curvature correction ON");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "       WARNING: curv terms are active in 2-D sense only!!!");
            }
        }
        if (ivmx == 13) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-epsilon (EASM,nonlinear)");
            if (iturbprod == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     strain-based production term");
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     exact production term");
            }
            if (ieasmcc2d == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     curvature correction ON");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "       WARNING: curv terms are active in 2-D sense only!!!");
            }
        }
        if (ivmx == 14) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-omega (EASM,nonlinear)");
            if (iturbprod == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     strain-based production term");
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     exact production term");
            }
            if (ieasmcc2d == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     curvature correction ON");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "       WARNING: curv terms are active in 2-D sense only!!!");
            }
        }
        if (ivmx == 15) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-enstrophy (Note: default tur1cut=-1 if no keyword)");
            if (ikoprod == 2) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     exact production term");
            } else if (ikoprod == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     strain-based production term");
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     approx (vort) production term");
            }
        }
        if (ivmx == 16) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     k-kL-MEAH2015");
            if (ikoprod == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     strain-based production term");
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     approx (vort) production term");
            }
        }
        if (ides == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   using model in conjunction with DES, cdes=%7.3f", (float)cdes);
        } else if (ides == 2) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   using model in conjunction with DDES, cdes=%7.3f", (float)cdes);
        } else if (ides == 3) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   using model in conjunction with MDDES, cdes=%7.3f, cddes=%7.3f", (float)cdes, (float)cddes);
        }
        if (ikoprod == 1 && (ivmx==6||ivmx==7||ivmx==10)) {
            if (std::abs((float)prod2d3dtrace-0.5f) < 0.001f) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     Sij used in 2SijSij prod term forced to be traceless in 2-D sense");
            } else if (std::abs((float)prod2d3dtrace-0.33333333f) < 0.001f) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     Sij used in 2SijSij prod term forced to be traceless in 3-D sense");
            }
        }
        if (ivmx==6 && i_wilcox06==1) {
            if (std::abs((float)prod2d3dtrace-0.5f) < 0.001f) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     Sij used in 2SijSij in W06 stress-limiter forced to be traceless in 2-D sense");
            } else if (std::abs((float)prod2d3dtrace-0.33333333f) < 0.001f) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     Sij used in 2SijSij in W06 stress-limiter forced to be traceless in 3-D sense");
            }
        }
        if (ivmx==9||ivmx==10||ivmx==11||ivmx==13||ivmx==15) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     WARNING: k-epsilon and k-enstrophy models sometimes fail to go turbulent!");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "         We recommend either restarting from a different converged model,");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "         or, if starting from scratch, freezing the model to its initialzed");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "         levels (using keyword NFREEZE) until the flowfield is converged enough");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "         for turbulence to sustain itself");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     ALWAYS check vist3d levels to insure turbulence has tripped!");
            if (ikoprod != 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     Also check that there is not excessive turbulence at stagnation regions.");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "         Try ikoprod=0 if there is.");
            }
        }
        if (ivmx==9||ivmx==10||ivmx==11||ivmx==13) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     Also, k-e models are NOT RECOMMENDED for");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "         adverse-pressure-gradient wall-bounded flows.");
            if (iturbord == 2) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     WARNING: 2nd order advection (iturbord=2)");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "          can converge poorly for k-epsilon models");
            }
        }
        if (ivmx==11||ivmx==12||ivmx==13||ivmx==14) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     WARNING: when using nonlinear models, the grid");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "         needs to have min y+ of order 1.  If it is much larger,");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "         you may need WallFunctions to keep from blowing up!");
        }
        if (ivmx==6||ivmx==7||ivmx==8||ivmx==12||ivmx==14) {
            if (i_catris_kw == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     Catris compressible correction for turb diffusion terms being used");
            }
        }
        if (ivmx==6||ivmx==7) {
            if (i_compress_correct == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     Wilcox-type dilatation-dissipation compressibility correction employed");
            } else if (i_compress_correct == 2) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     Zeman-BL-type dilatation-dissipation compressibility correction employed");
            }
            if (isstrc == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     SSTRC-type curvature correction employed (AIAA 98-2554), sstrc_crc=%5.2f", (float)sstrc_crc);
            }
            if (isstrc == 2) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     SSTRC-type curvature correction employed (Smirnov & Menter)");
            }
            if (isstsf == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "     SST-sf separation fix correction employed");
            }
        }
        if (keepambient == 1 && ivmx != 15) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     2-eqn ambient turbulence levels not allowed to decay");
        }
        if (iaxi2planeturb == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     2-eqn model ignoring i-dir");
        }
        if (istrongturbdis == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "     strong conserv - diss terms");
        }
    } // end isklton > 0 part 2


    // Setup constants
    xminn = std::pow(10.0, (double)(-iexp+1));
    catris_kw = 0.0;
    if (i_catris_kw == 1) catris_kw = 1.0;
    if (i_compress_correct == 1) {
        // Wilcox-type (Turb Modeling for CFD, Ed 3, 2006, p. 258)
        turb_mach_0 = 0.25;
        xsi_star = 2.0;
        gam_comp = 0.0;
    } else if (i_compress_correct == 2) {
        // Zeman B.L.-type (AIAA 93-0897)
        turb_mach_0 = 0.2;
        xsi_star = 0.75;
        gam_comp = 0.66;
    } else {
        turb_mach_0 = 0.0; xsi_star = 0.0; gam_comp = 0.0;
    }
    nsubit = nsubturb;
    // Set CFL factors
    if (ivmx == 6) factor = 2.;
    if (ivmx == 7) factor = 10.;
    if (ivmx == 8 || ivmx == 12 || ivmx == 14) factor = 5.;
    if (ivmx == 10 || ivmx == 11) factor = 2.;
    if (ivmx == 9 || ivmx == 13) factor = 5.;
    if (ivmx == 15) factor = 10.;
    if (ivmx == 16) factor = 10.;
    factor1 = factor;
    factor2 = factor;
    if ((float)cflturb[0] != 0.) factor1 = cflturb[0];
    if ((float)cflturb[1] != 0.) factor2 = cflturb[1];
    factor2 = factor2 / factor1;
    // Timestep for turb model
    if ((float)dt < 0) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            timestp(j,k,i) = factor1 * vol(j,k,i) / dtj(j,k,i);
            timestp(j,k,i) = ccmincr(timestp(j,k,i), 100.);
        }
    } else {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            timestp(j,k,i) = dt;
            factor2 = 1.;
        }
    }
    // Set up constants
    vk = .41;
    a1 = .31;
    if (i_bsl == 1) a1 = 1.e8;
    cmuc1 = 0.09;
    if (ivmx==8||ivmx==11||ivmx==12||ivmx==9||ivmx==13||ivmx==14) {
        cmuc1 = 0.081;
        if (ivmx==9||ivmx==13) cmuc1 = 0.0885;
        if (ivmx==8||ivmx==14) cmuc1 = 0.0895;
        if (ieasm_type == 1) {
            c10 = 3.6; c11 = 0.; c2 = 0.8;
            constt = 5./9.;
            c3 = 6./11.*(2.+(3.*constt));
            c4 = 2./11.*(10.-(7.*constt));
        } else if (ieasm_type == 2) {
            c10 = 3.0; c11 = 0.; c2 = 0.8;
            constt = 0.4;
            c3 = 6./11.*(2.+(3.*constt));
            c4 = 2./11.*(10.-(7.*constt));
        } else if (ieasm_type == 3 || ieasm_type == 4) {
            c10 = 3.4; c11 = 1.8; c2 = 1.2; c3 = 1.25; c4 = 0.40;
            css = 0.84;
        } else {
            c10 = 3.4; c11 = 1.8; c2 = 0.36; c3 = 1.25; c4 = 0.40; c5 = 1.88;
            css = 0.0;
            gg = 1./(c10+c5-1.);
        }
        idurbinlim = 0;
    } else {
        css = 0.0; gg = 0.0; idurbinlim = 0;
        constt = 0.0;
    }
    // beta1
    if (ivmx==6 && i_wilcox98==1) beta1 = 0.072;
    if (ivmx==6 && i_wilcox98!=1) beta1 = 0.075;
    if (ivmx==7) beta1 = 0.075;
    if (ivmx==12) beta1 = 0.83;
    if (ivmx==8||ivmx==14) beta1 = 0.83;
    if (ivmx==10||ivmx==11||ivmx==9||ivmx==13) beta1 = 1.83;
    // sigo1
    if (ivmx==6||ivmx==7) sigo1 = 0.5;
    if (ivmx==12) sigo1 = 1./2.20;
    if (ivmx==10) sigo1 = 0.71429;
    if (ivmx==11) sigo1 = 1./1.3;
    if (ivmx==15) sigo1 = 1.46;
    // sigk1
    if (ivmx==6) sigk1 = 0.5;
    if (ivmx==7 && i_bsl!=1) sigk1 = 0.85;
    if (ivmx==7 && i_bsl==1) sigk1 = 0.5;
    if (ivmx==12) sigk1 = 1./1.4;
    if (ivmx==10||ivmx==11||ivmx==9||ivmx==13) sigk1 = 1.00;
    if (ivmx==8||ivmx==14) sigk1 = 1./1.0;
    if (ivmx==15) sigk1 = 1.8;
    // alp1
    if (ivmx==6 && i_wilcox98==1) alp1 = 13./25.;
    if (ivmx==6 && i_wilcox98!=1) alp1 = beta1/cmuc1 - sigo1*vk*vk/std::sqrt(cmuc1);
    if (ivmx==7 && isst2003!=1) alp1 = beta1/cmuc1 - sigo1*vk*vk/std::sqrt(cmuc1);
    if (ivmx==7 && isst2003==1) alp1 = 5./9.;
    if (ivmx==10) alp1 = 1.45;
    if (ivmx==11||ivmx==12) alp1 = beta1 - vk*vk*sigo1/std::sqrt(cmuc1);
    if (ivmx==9||ivmx==13) alp1 = 1.44;
    if (ivmx==8||ivmx==14) alp1 = 0.53;
    if (ivmx==9||ivmx==13) sigo1 = std::sqrt(cmuc1)*(beta1-alp1)/(vk*vk);
    if (ivmx==8||ivmx==14) sigo1 = std::sqrt(cmuc1)*(beta1-alp1)/(vk*vk);
    if ((ivmx==8||ivmx==14) && (ieasm_type==3||ieasm_type==4)) sigo1 = 0.65;
    // sigkmu
    sigkmu = 1.0;
    if (ivmx==15) sigkmu = 1./3.;
    // Constants for Set 2
    if (ivmx==6 && i_wilcox06==1) {
        cmuc2 = 0.09; beta2 = 0.0708; sigo2 = 0.5; sigk2 = 0.6;
        alp2 = 13./25.; sigd0 = 0.125;
    } else {
        cmuc2 = 0.09; beta2 = 0.0828; sigo2 = 0.856; sigk2 = 1.00;
        alp2 = beta2/cmuc2 - sigo2*vk*vk/std::sqrt(cmuc2);
        if (isst2003 == 1) alp2 = 0.44;
        sigd0 = 0.0;
    }
    // k-kL-MEAH2015 constants
    if (ivmx == 16) {
        sigk1 = 1.0; sigo1 = 1.0;
        zeta1_kkl = 1.2; zeta2_kkl = 0.97; zeta3_kkl = 0.13;
        c11_kkl = 10.; c12_kkl = 1.3; cd1_kkl = 4.7;
    } else {
        zeta1_kkl = 0.; zeta2_kkl = 0.; zeta3_kkl = 0.;
        c11_kkl = 0.; c12_kkl = 0.; cd1_kkl = 0.;
    }
    jd2 = (jdim-1)/2;
    re = (double)reue / (double)xmach;
    c2b = (double)cbar / (double)tinf;
    c2bp = c2b + 1.0;
    iwrite = 0;
    // Check nfreeze
    if (icyc <= nfreeze) {
        if (isklton > 0) {
            nss = std::min(ncyc, nfreeze);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " ");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " turbulence model is frozen for %5d iterations or subits", nss);
        }
        sumn1 = 0.;
        sumn2 = 0.;
        negn1 = 0;
        negn2 = 0;
        return;
    }


    // phi setup
    phi = 0.;
    if ((float)dt > 0.) {
        if (std::abs(ita) == 2) {
            phi = 0.5;
        } else {
            phi = 0.;
        }
        if (itaturb == 0) {
            phi = 0.;
            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   turb model is 1st order in time");
            }
        } else {
            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   turb model is same order in time as mean flow eqns");
            }
        }
    }
    // Get laminar viscosity at cell centers
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++)
    for (j = 1; j <= jdim-1; j++) {
        tt = (double)gamma * q(j,k,i,5) / q(j,k,i,1);
        fnu(j,k,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
    }
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++) {
        tt = (double)gamma * qj0(k,i,5,1) / qj0(k,i,1,1);
        fnu(0,k,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
        tt = (double)gamma * qj0(k,i,5,3) / qj0(k,i,1,3);
        fnu(jdim,k,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
    }
    for (i = 1; i <= idim-1; i++)
    for (j = 1; j <= jdim-1; j++) {
        tt = (double)gamma * qk0(j,i,5,1) / qk0(j,i,1,1);
        fnu(j,0,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
        tt = (double)gamma * qk0(j,i,5,3) / qk0(j,i,1,3);
        fnu(j,kdim,i) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
    }
    if (i2d != 1 && iaxi2planeturb != 1) {
        for (j = 1; j <= jdim-1; j++)
        for (k = 1; k <= kdim-1; k++) {
            tt = (double)gamma * qi0(j,k,5,1) / qi0(j,k,1,1);
            fnu(j,k,0) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
            tt = (double)gamma * qi0(j,k,5,3) / qi0(j,k,1,3);
            fnu(j,k,idim) = c2bp*tt*std::sqrt(tt)/(c2b+tt);
        }
    }
    // Load appropriate turb viscosity at cell centers
    if (ivmx == 12) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            v3dtmp(j,k,i) = cmuc1*q(j,k,i,1)*zksav(j,k,i,2)/zksav(j,k,i,1);
            v3dtmp(j,k,i) = ccmin(v3dtmp(j,k,i), (double)edvislim);
        }
    } else if (ivmx==11||ivmx==9||ivmx==13) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            v3dtmp(j,k,i) = cmuc1*q(j,k,i,1)*zksav(j,k,i,2)*zksav(j,k,i,2)/zksav(j,k,i,1);
            v3dtmp(j,k,i) = ccmin(v3dtmp(j,k,i), (double)edvislim);
        }
    } else if (ivmx==6 && i_wilcox06==1) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            v3dtmp(j,k,i) = q(j,k,i,1)*zksav(j,k,i,2)/zksav(j,k,i,1);
            v3dtmp(j,k,i) = ccmin(v3dtmp(j,k,i), (double)edvislim);
        }
    } else {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            v3dtmp(j,k,i) = vist3d(j,k,i);
        }
    }
    // Load vist3d value into ghost cells
    for (i = 1; i <= idim-1; i++)
    for (j = 1; j <= jdim-1; j++) {
        v3dtmp(j,0,i) = bck(j,i,1)*(iwf[2]*v3dtmp(j,1,i) +
            (1-iwf[2])*2.*vk0(j,i,1,1)-v3dtmp(j,1,i)) +
            (1.-bck(j,i,1))*vk0(j,i,1,1);
        v3dtmp(j,kdim,i) = bck(j,i,2)*(iwf[2]*v3dtmp(j,kdim-1,i) +
            (1-iwf[2])*2.*vk0(j,i,1,3)-v3dtmp(j,kdim-1,i)) +
            (1.-bck(j,i,2))*vk0(j,i,1,3);
    }
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++) {
        v3dtmp(0,k,i) = bcj(k,i,1)*(iwf[1]*v3dtmp(1,k,i) +
            (1-iwf[1])*2.*vj0(k,i,1,1)-v3dtmp(1,k,i)) +
            (1.-bcj(k,i,1))*vj0(k,i,1,1);
        v3dtmp(jdim,k,i) = bcj(k,i,2)*(iwf[1]*v3dtmp(jdim-1,k,i) +
            (1-iwf[1])*2.*vj0(k,i,1,3)-v3dtmp(jdim-1,k,i)) +
            (1.-bcj(k,i,2))*vj0(k,i,1,3);
    }
    if (i2d != 1 && iaxi2planeturb != 1) {
        for (j = 1; j <= jdim-1; j++)
        for (k = 1; k <= kdim-1; k++) {
            v3dtmp(j,k,0) = bci(j,k,1)*(iwf[0]*v3dtmp(j,k,1) +
                (1-iwf[0])*2.*vi0(j,k,1,1)-v3dtmp(j,k,1)) +
                (1.-bci(j,k,1))*vi0(j,k,1,1);
            v3dtmp(j,k,idim) = bci(j,k,2)*(iwf[0]*v3dtmp(j,k,idim-1) +
                (1-iwf[0])*2.*vi0(j,k,1,3)-v3dtmp(j,k,idim-1)) +
                (1.-bci(j,k,2))*vi0(j,k,1,3);
        }
    }
    // Save zksav at time step n for time-accurate computation
    if ((float)dt > 0. && icyc == 1) {
        if (std::abs(ita) == 2) {
            if ((float)zksav2(1,1,1,1) == 0.) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    zksav2(j,k,i,3) = 0.;
                    zksav2(j,k,i,4) = 0.;
                }
            } else {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    zksav2(j,k,i,3) = zksav(j,k,i,1) - zksav2(j,k,i,1);
                    zksav2(j,k,i,4) = zksav(j,k,i,2) - zksav2(j,k,i,2);
                }
            }
        }
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            zksav2(j,k,i,1) = zksav(j,k,i,1);
            zksav2(j,k,i,2) = zksav(j,k,i,2);
        }
    }
    // Get TURRE values
    for (i = 1; i <= idim-1; i++)
    for (k = 1; k <= kdim-1; k++)
    for (j = 1; j <= jdim-1; j++) {
        turre(j,k,i,1) = zksav(j,k,i,1);
        turre(j,k,i,2) = zksav(j,k,i,2);
    }


    // Iterate to solve the equations
    for (not_iter = 1; not_iter <= nsubit; not_iter++) {
        // Set up boundary conditions (ghost cells)
        // (1) k=0 boundary:
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            turre(j,0,i,1) = tk0(j,i,1,1);
            turre(j,0,i,2) = tk0(j,i,2,1);
        }
        // (2) k=kdim boundary:
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            turre(j,kdim,i,1) = tk0(j,i,1,3);
            turre(j,kdim,i,2) = tk0(j,i,2,3);
        }
        // (3) j=0 boundary:
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            turre(0,k,i,1) = tj0(k,i,1,1);
            turre(0,k,i,2) = tj0(k,i,2,1);
        }
        // (4) j=jdim boundary:
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++) {
            turre(jdim,k,i,1) = tj0(k,i,1,3);
            turre(jdim,k,i,2) = tj0(k,i,2,3);
        }
        if (i2d != 1 && iaxi2planeturb != 1) {
            // (5) i=0 boundary:
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,k,0,1) = ti0(j,k,1,1);
                turre(j,k,0,2) = ti0(j,k,2,1);
            }
            // (6) i=idim boundary:
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,k,idim,1) = ti0(j,k,1,3);
                turre(j,k,idim,2) = ti0(j,k,2,3);
            }
        }
        if (iturbord != 1) {
            // (1) k=0 boundary:
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,-1,i,1) = tk0(j,i,1,2);
                turre(j,-1,i,2) = tk0(j,i,2,2);
            }
            // (2) k=kdim boundary:
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                turre(j,kdim+1,i,1) = tk0(j,i,1,4);
                turre(j,kdim+1,i,2) = tk0(j,i,2,4);
            }
            // (3) j=0 boundary:
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                turre(-1,k,i,1) = tj0(k,i,1,2);
                turre(-1,k,i,2) = tj0(k,i,2,2);
            }
            // (4) j=jdim boundary:
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                turre(jdim+1,k,i,1) = tj0(k,i,1,4);
                turre(jdim+1,k,i,2) = tj0(k,i,2,4);
            }
            if (i2d != 1 && iaxi2planeturb != 1) {
                // (5) i=0 boundary:
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    turre(j,k,-1,1) = ti0(j,k,1,2);
                    turre(j,k,-1,2) = ti0(j,k,2,2);
                }
                // (6) i=idim boundary:
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    turre(j,k,idim+1,1) = ti0(j,k,1,4);
                    turre(j,k,idim+1,2) = ti0(j,k,2,4);
                }
            }
        }
        // Set default for damp1 (when not used)
        for (k = 1; k <= kdim-1; k++)
        for (i = 1; i <= idim-1; i++)
        for (j = 1; j <= jdim-1; j++) {
            damp1(j,k,i) = 0.0;
        }

        // Get damp1 = CD = cross derivative term for SST or EASM k-o models only:
        if (ivmx==7||ivmx==8||ivmx==14) {
            double temp, xderiv, arg1, arga, argb, arg;
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tt = xa*xa+ya*ya+za*za;
                ca = 2.*sigo2*tt/(turre(j,k,i,1)*re);
                damp1(j,k,i) = 0.25*ca*(turre(j,k+1,i,1)-turre(j,k-1,i,1))*(turre(j,k+1,i,2)-turre(j,k-1,i,2));
            }
            for (j = 1; j <= jdim-1; j++)
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tt = xa*xa+ya*ya+za*za;
                ca = 2.*sigo2*tt/(turre(j,k,i,1)*re);
                damp1(j,k,i) += 0.25*ca*(turre(j+1,k,i,1)-turre(j-1,k,i,1))*(turre(j+1,k,i,2)-turre(j-1,k,i,2));
            }
            if (i2d != 1 && iaxi2planeturb != 1) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tt = xa*xa+ya*ya+za*za;
                    ca = 2.*sigo2*tt/(turre(j,k,i,1)*re);
                    damp1(j,k,i) += 0.25*ca*(turre(j,k,i+1,1)-turre(j,k,i-1,1))*(turre(j,k,i+1,2)-turre(j,k,i-1,2));
                }
            }
            if (ivmx==8||ivmx==14) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    temp = cmuc1*cmuc1*damp1(j,k,i)/(2.*sigo2*re*turre(j,k,i,1)*turre(j,k,i,1));
                    temp = ccmaxcr(temp, 0.);
                    damp1(j,k,i) = (1.+680.*temp*temp)/(1.+400.*temp*temp);
                }
            }
        }
        // Get damp1 = CD = cross derivative term for Wilcox06 model only:
        if (ivmx==6 && i_wilcox06==1) {
            double xderiv;
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tt = xa*xa+ya*ya+za*za;
                xderiv = 0.25*(turre(j,k+1,i,1)-turre(j,k-1,i,1))*(turre(j,k+1,i,2)-turre(j,k-1,i,2));
                damp1(j,k,i) = xderiv*sigd0*tt/(turre(j,k,i,1)*re);
            }
            for (j = 1; j <= jdim-1; j++)
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tt = xa*xa+ya*ya+za*za;
                xderiv = 0.25*(turre(j+1,k,i,1)-turre(j-1,k,i,1))*(turre(j+1,k,i,2)-turre(j-1,k,i,2));
                damp1(j,k,i) += xderiv*sigd0*tt/(turre(j,k,i,1)*re);
            }
            if (i2d != 1 && iaxi2planeturb != 1) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tt = xa*xa+ya*ya+za*za;
                    xderiv = 0.25*(turre(j,k,i+1,1)-turre(j,k,i-1,1))*(turre(j,k,i+1,2)-turre(j,k,i-1,2));
                    damp1(j,k,i) += xderiv*sigd0*tt/(turre(j,k,i,1)*re);
                }
            }
            // for Wilcox06, damp1 only active if it is positive
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                damp1(j,k,i) = ccmaxcr(damp1(j,k,i), 0.);
            }
        }
        // Get damp1 = CD = cross derivative term for Wilcox98 model only:
        if (ivmx==6 && i_wilcox98==1) {
            double xderiv;
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tt = xa*xa+ya*ya+za*za;
                xderiv = 0.25*(turre(j,k+1,i,1)-turre(j,k-1,i,1))*(turre(j,k+1,i,2)-turre(j,k-1,i,2));
                damp1(j,k,i) = xderiv*tt/(turre(j,k,i,1)*turre(j,k,i,1)*turre(j,k,i,1)*re*re);
            }
            for (j = 1; j <= jdim-1; j++)
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tt = xa*xa+ya*ya+za*za;
                xderiv = 0.25*(turre(j+1,k,i,1)-turre(j-1,k,i,1))*(turre(j+1,k,i,2)-turre(j-1,k,i,2));
                damp1(j,k,i) += xderiv*tt/(turre(j,k,i,1)*turre(j,k,i,1)*turre(j,k,i,1)*re*re);
            }
            if (i2d != 1 && iaxi2planeturb != 1) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tt = xa*xa+ya*ya+za*za;
                    xderiv = 0.25*(turre(j,k,i+1,1)-turre(j,k,i-1,1))*(turre(j,k,i+1,2)-turre(j,k,i-1,2));
                    damp1(j,k,i) += xderiv*tt/(turre(j,k,i,1)*turre(j,k,i,1)*turre(j,k,i,1)*re*re);
                }
            }
            // for Wilcox98, damp1 only active if it is positive
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                damp1(j,k,i) = ccmaxcr(damp1(j,k,i), 0.);
            }
        }

        // Get blend = F1 factor
        if (ivmx == 7) {
            double small_val = 1.e-20;
            if (isst2003 == 1) small_val = 1.e-10;
            double arg1, arga, argb, arg;
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                arg1 = std::sqrt(turre(j,k,i,2))/(.09*re*turre(j,k,i,1)*ccabs(smin(j,k,i)));
                arg2 = 500.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*re*re*smin(j,k,i)*turre(j,k,i,1));
                arga = ccmax(arg1, arg2);
                double temp_val = ccmaxcr(damp1(j,k,i)*re, small_val);
                argb = 4.*sigo2*turre(j,k,i,2)/(temp_val*smin(j,k,i)*smin(j,k,i));
                arg = ccmin(arga, argb);
                blend(j,k,i) = cctanh(arg*arg*arg*arg);
            }
        } else if (ivmx==6 && i_wilcox06==1) {
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                blend(j,k,i) = 0.0;
            }
        } else {
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                blend(j,k,i) = 1.0;
            }
        }

        // F_eta_eta viscous terms - Interior points
        for (k = 2; k <= kdim-2; k++) {
            kl = k-1; ku = k+1;
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
                dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
                volku = vol(j,ku,i);
                xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                volkl = vol(j,kl,i);
                xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                fnum = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                rhop = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j,k+1,i,1)+q(j,k,i,1)));
                rhom = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j,k-1,i,1)+q(j,k,i,1)));
                cdp = (fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,1) = (1.0-catris_kw)*(-byy*turre(j,k-1,i,1)-cyy*turre(j,k,i,1)-dyy*turre(j,k+1,i,1)) +
                    catris_kw*(-byy*std::sqrt(q(j,k-1,i,1))*turre(j,k-1,i,1)-cyy*std::sqrt(q(j,k,i,1))*turre(j,k,i,1)-dyy*std::sqrt(q(j,k+1,i,1))*turre(j,k+1,i,1));
                rhop = (1.0-catris_kw)+catris_kw*.5*(q(j,k+1,i,1)+q(j,k,i,1));
                rhom = (1.0-catris_kw)+catris_kw*.5*(q(j,k-1,i,1)+q(j,k,i,1));
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,2) = (1.0-catris_kw)*(-byy*turre(j,k-1,i,2)-cyy*turre(j,k,i,2)-dyy*turre(j,k+1,i,2)) +
                    catris_kw*(-byy*q(j,k-1,i,1)*turre(j,k-1,i,2)-cyy*q(j,k,i,1)*turre(j,k,i,2)-dyy*q(j,k+1,i,1)*turre(j,k+1,i,2));
            }
        }
        // F_eta_eta K0 boundary points
        {
            k = 1; kl = 1; ku = std::min(2, kdim-1);
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
                dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
                volku = vol(j,ku,i);
                xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                volkl = volk0(j,i,1);
                xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                fnum = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                rhop = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j,k+1,i,1)+q(j,k,i,1)));
                rhom = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(qk0(j,i,1,1)+q(j,k,i,1)));
                cdp = (fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,1) = (1.0-catris_kw)*(-byy*turre(j,k-1,i,1)-cyy*turre(j,k,i,1)-dyy*turre(j,k+1,i,1)) +
                    catris_kw*(-byy*std::sqrt(qk0(j,i,1,1))*turre(j,k-1,i,1)-cyy*std::sqrt(q(j,k,i,1))*turre(j,k,i,1)-dyy*std::sqrt(q(j,k+1,i,1))*turre(j,k+1,i,1));
                rhop = (1.0-catris_kw)+catris_kw*.5*(q(j,k+1,i,1)+q(j,k,i,1));
                rhom = (1.0-catris_kw)+catris_kw*.5*(qk0(j,i,1,1)+q(j,k,i,1));
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,2) = (1.0-catris_kw)*(-byy*turre(j,k-1,i,2)-cyy*turre(j,k,i,2)-dyy*turre(j,k+1,i,2)) +
                    catris_kw*(-byy*qk0(j,i,1,1)*turre(j,k-1,i,2)-cyy*q(j,k,i,1)*turre(j,k,i,2)-dyy*q(j,k+1,i,1)*turre(j,k+1,i,2));
            }
        }

        // F_eta_eta KDIM boundary points
        {
            k = kdim-1; kl = kdim-2; ku = kdim-1;
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                dfacep = 0.5*(blend(j,k,i)+blend(j,ku,i));
                dfacem = 0.5*(blend(j,k,i)+blend(j,kl,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
                volku = volk0(j,i,3);
                xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                yp = sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                volkl = vol(j,kl,i);
                xm = sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                ym = sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                zm = sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                xa = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                fnup = .5*(fnu(j,k+1,i)+fnu(j,k,i));
                fnum = .5*(fnu(j,k-1,i)+fnu(j,k,i));
                rhop = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(qk0(j,i,1,3)+q(j,k,i,1)));
                rhom = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j,k-1,i,1)+q(j,k,i,1)));
                cdp = (fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,1) = (1.0-catris_kw)*(-byy*turre(j,k-1,i,1)-cyy*turre(j,k,i,1)-dyy*turre(j,k+1,i,1)) +
                    catris_kw*(-byy*std::sqrt(q(j,k-1,i,1))*turre(j,k-1,i,1)-cyy*std::sqrt(q(j,k,i,1))*turre(j,k,i,1)-dyy*std::sqrt(qk0(j,i,1,3))*turre(j,k+1,i,1));
                rhop = (1.0-catris_kw)+catris_kw*.5*(qk0(j,i,1,3)+q(j,k,i,1));
                rhom = (1.0-catris_kw)+catris_kw*.5*(q(j,k-1,i,1)+q(j,k,i,1));
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                byy = -cdm; cyy = cdp+cdm; dyy = -cdp;
                rhside(j,k,i,2) = (1.0-catris_kw)*(-byy*turre(j,k-1,i,2)-cyy*turre(j,k,i,2)-dyy*turre(j,k+1,i,2)) +
                    catris_kw*(-byy*q(j,k-1,i,1)*turre(j,k-1,i,1)-cyy*q(j,k,i,1)*turre(j,k,i,1)-dyy*qk0(j,i,1,3)*turre(j,k+1,i,1));
            }
        }
        // Advective terms in eta
        if (iturbord == 1) {
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                sgnu = ccsignrc(1.,uu);
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                rhside(j,k,i,1) -= uu*(app*(turre(j,k,i,1)-turre(j,k-1,i,1))+apm*(turre(j,k+1,i,1)-turre(j,k,i,1)));
                rhside(j,k,i,2) -= uu*(app*(turre(j,k,i,2)-turre(j,k-1,i,2))+apm*(turre(j,k+1,i,2)-turre(j,k,i,2)));
            }
        } else {
            // 2nd order upwind; LHS remains 1st order everywhere
            for (k = 1; k <= kdim-1; k++)
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                sgnu = ccsignrc(1.,uu);
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                rhside(j,k,i,1) += -0.5*uu*app*turre(j,k-2,i,1)+2.*uu*app*turre(j,k-1,i,1)-1.5*uu*app*turre(j,k,i,1)
                    +1.5*uu*apm*turre(j,k,i,1)-2.*uu*apm*turre(j,k+1,i,1)+0.5*uu*apm*turre(j,k+2,i,1);
                rhside(j,k,i,2) += -0.5*uu*app*turre(j,k-2,i,2)+2.*uu*app*turre(j,k-1,i,2)-1.5*uu*app*turre(j,k,i,2)
                    +1.5*uu*apm*turre(j,k,i,2)-2.*uu*apm*turre(j,k+1,i,2)+0.5*uu*apm*turre(j,k+2,i,2);
            }
        }

        // F_xi_xi viscous terms - Interior points
        for (j = 2; j <= jdim-2; j++) {
            jl = j-1; ju = j+1;
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
                dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
                volju = vol(ju,k,i);
                xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                voljl = vol(jl,k,i);
                xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                fnum = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                rhop = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j+1,k,i,1)+q(j,k,i,1)));
                rhom = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j-1,k,i,1)+q(j,k,i,1)));
                cdp = (fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,1) += (1.0-catris_kw)*(-bxx*turre(j-1,k,i,1)-cxx*turre(j,k,i,1)-dxx*turre(j+1,k,i,1)) +
                    catris_kw*(-bxx*std::sqrt(q(j-1,k,i,1))*turre(j-1,k,i,1)-cxx*std::sqrt(q(j,k,i,1))*turre(j,k,i,1)-dxx*std::sqrt(q(j+1,k,i,1))*turre(j+1,k,i,1));
                rhop = (1.0-catris_kw)+catris_kw*.5*(q(j+1,k,i,1)+q(j,k,i,1));
                rhom = (1.0-catris_kw)+catris_kw*.5*(q(j-1,k,i,1)+q(j,k,i,1));
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,2) += (1.0-catris_kw)*(-bxx*turre(j-1,k,i,2)-cxx*turre(j,k,i,2)-dxx*turre(j+1,k,i,2)) +
                    catris_kw*(-bxx*q(j-1,k,i,1)*turre(j-1,k,i,2)-cxx*q(j,k,i,1)*turre(j,k,i,2)-dxx*q(j+1,k,i,1)*turre(j+1,k,i,2));
            }
        }
        // F_xi_xi J0 boundary points
        {
            j = 1; jl = 1; ju = std::min(2, jdim-1);
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
                dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
                volju = vol(ju,k,i);
                xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                voljl = volj0(k,i,1);
                xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                fnum = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                rhop = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j+1,k,i,1)+q(j,k,i,1)));
                rhom = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(qj0(k,i,1,1)+q(j,k,i,1)));
                cdp = (fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,1) += (1.0-catris_kw)*(-bxx*turre(j-1,k,i,1)-cxx*turre(j,k,i,1)-dxx*turre(j+1,k,i,1)) +
                    catris_kw*(-bxx*std::sqrt(qj0(k,i,1,1))*turre(j-1,k,i,1)-cxx*std::sqrt(q(j,k,i,1))*turre(j,k,i,1)-dxx*std::sqrt(q(j+1,k,i,1))*turre(j+1,k,i,1));
                rhop = (1.0-catris_kw)+catris_kw*.5*(q(j+1,k,i,1)+q(j,k,i,1));
                rhom = (1.0-catris_kw)+catris_kw*.5*(qj0(k,i,1,1)+q(j,k,i,1));
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,2) += (1.0-catris_kw)*(-bxx*turre(j-1,k,i,2)-cxx*turre(j,k,i,2)-dxx*turre(j+1,k,i,2)) +
                    catris_kw*(-bxx*qj0(k,i,1,1)*turre(j-1,k,i,1)-cxx*q(j,k,i,1)*turre(j,k,i,1)-dxx*q(j+1,k,i,1)*turre(j+1,k,i,1));
            }
        }

        // F_xi_xi JDIM boundary points
        {
            j = jdim-1; jl = jdim-2; ju = jdim-1;
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++) {
                dfacep = 0.5*(blend(j,k,i)+blend(ju,k,i));
                dfacem = 0.5*(blend(j,k,i)+blend(jl,k,i));
                sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
                volju = volj0(k,i,3);
                xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                yp = sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                voljl = vol(jl,k,i);
                xm = sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                ym = sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                zm = sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                xa = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                ya = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                za = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                ttpo = xp*xa+yp*ya+zp*za;
                ttmo = xm*xa+ym*ya+zm*za;
                ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                fnup = .5*(fnu(j+1,k,i)+fnu(j,k,i));
                fnum = .5*(fnu(j-1,k,i)+fnu(j,k,i));
                rhop = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(qj0(k,i,1,3)+q(j,k,i,1)));
                rhom = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j-1,k,i,1)+q(j,k,i,1)));
                cdp = (fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,1) += (1.0-catris_kw)*(-bxx*turre(j-1,k,i,1)-cxx*turre(j,k,i,1)-dxx*turre(j+1,k,i,1)) +
                    catris_kw*(-bxx*std::sqrt(q(j-1,k,i,1))*turre(j-1,k,i,1)-cxx*std::sqrt(q(j,k,i,1))*turre(j,k,i,1)-dxx*std::sqrt(qj0(k,i,1,3))*turre(j+1,k,i,1));
                rhop = (1.0-catris_kw)+catris_kw*.5*(qj0(k,i,1,3)+q(j,k,i,1));
                rhom = (1.0-catris_kw)+catris_kw*.5*(q(j-1,k,i,1)+q(j,k,i,1));
                cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                bxx = -cdm; cxx = cdp+cdm; dxx = -cdp;
                rhside(j,k,i,2) += (1.0-catris_kw)*(-bxx*turre(j-1,k,i,2)-cxx*turre(j,k,i,2)-dxx*turre(j+1,k,i,2)) +
                    catris_kw*(-bxx*q(j-1,k,i,1)*turre(j-1,k,i,1)-cxx*q(j,k,i,1)*turre(j,k,i,1)-dxx*qj0(k,i,1,3)*turre(j+1,k,i,1));
            }
        }
        // Advective terms in xi
        if (iturbord == 1) {
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                sgnu = ccsignrc(1.,uu);
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                rhside(j,k,i,1) -= uu*(app*(turre(j,k,i,1)-turre(j-1,k,i,1))+apm*(turre(j+1,k,i,1)-turre(j,k,i,1)));
                rhside(j,k,i,2) -= uu*(app*(turre(j,k,i,2)-turre(j-1,k,i,2))+apm*(turre(j+1,k,i,2)-turre(j,k,i,2)));
            }
        } else {
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                yc = 0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                sgnu = ccsignrc(1.,uu);
                app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                rhside(j,k,i,1) += -0.5*uu*app*turre(j-2,k,i,1)+2.*uu*app*turre(j-1,k,i,1)-1.5*uu*app*turre(j,k,i,1)
                    +1.5*uu*apm*turre(j,k,i,1)-2.*uu*apm*turre(j+1,k,i,1)+0.5*uu*apm*turre(j+2,k,i,1);
                rhside(j,k,i,2) += -0.5*uu*app*turre(j-2,k,i,2)+2.*uu*app*turre(j-1,k,i,2)-1.5*uu*app*turre(j,k,i,2)
                    +1.5*uu*apm*turre(j,k,i,2)-2.*uu*apm*turre(j+1,k,i,2)+0.5*uu*apm*turre(j+2,k,i,2);
            }
        }

        // F_zeta_zeta viscous terms
        if (i2d != 1 && iaxi2planeturb != 1) {
            // Interior points
            for (i = 2; i <= idim-2; i++) {
                il = i-1; iu = i+1;
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                    dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                    sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
                    voliu = vol(j,k,iu);
                    xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    volil = vol(j,k,il);
                    xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                    anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                    fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                    fnum = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                    rhop = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j,k,i+1,1)+q(j,k,i,1)));
                    rhom = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j,k,i-1,1)+q(j,k,i,1)));
                    cdp = (fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,1) += (1.0-catris_kw)*(-bzz*turre(j,k,i-1,1)-czz*turre(j,k,i,1)-dzz*turre(j,k,i+1,1)) +
                        catris_kw*(-bzz*std::sqrt(q(j,k,i-1,1))*turre(j,k,i-1,1)-czz*std::sqrt(q(j,k,i,1))*turre(j,k,i,1)-dzz*std::sqrt(q(j,k,i+1,1))*turre(j,k,i+1,1));
                    rhop = (1.0-catris_kw)+catris_kw*.5*(q(j,k,i+1,1)+q(j,k,i,1));
                    rhom = (1.0-catris_kw)+catris_kw*.5*(q(j,k,i-1,1)+q(j,k,i,1));
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,2) += (1.0-catris_kw)*(-bzz*turre(j,k,i-1,2)-czz*turre(j,k,i,2)-dzz*turre(j,k,i+1,2)) +
                        catris_kw*(-bzz*q(j,k,i-1,1)*turre(j,k,i-1,1)-czz*q(j,k,i,1)*turre(j,k,i,1)-dzz*q(j,k,i+1,1)*turre(j,k,i+1,1));
                }
            }
            // I0 boundary points
            {
                i = 1; il = 1; iu = std::min(2, idim-1);
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                    dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                    sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
                    voliu = vol(j,k,iu);
                    xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    volil = voli0(j,k,1);
                    xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                    anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                    fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                    fnum = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                    rhop = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j,k,iu,1)+q(j,k,i,1)));
                    rhom = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(qi0(j,k,1,1)+q(j,k,i,1)));
                    cdp = (fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,1) += (1.0-catris_kw)*(-bzz*turre(j,k,i-1,1)-czz*turre(j,k,i,1)-dzz*turre(j,k,i+1,1)) +
                        catris_kw*(-bzz*std::sqrt(qi0(j,k,1,1))*turre(j,k,i-1,1)-czz*std::sqrt(q(j,k,i,1))*turre(j,k,i,1)-dzz*std::sqrt(q(j,k,iu,1))*turre(j,k,i+1,1));
                    rhop = (1.0-catris_kw)+catris_kw*.5*(q(j,k,iu,1)+q(j,k,i,1));
                    rhom = (1.0-catris_kw)+catris_kw*.5*(qi0(j,k,1,1)+q(j,k,i,1));
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,2) += (1.0-catris_kw)*(-bzz*turre(j,k,i-1,2)-czz*turre(j,k,i,2)-dzz*turre(j,k,i+1,2)) +
                        catris_kw*(-bzz*qi0(j,k,1,1)*turre(j,k,i-1,1)-czz*q(j,k,i,1)*turre(j,k,i,1)-dzz*q(j,k,iu,1)*turre(j,k,i+1,1));
                }
            }

            // IDIM boundary points
            {
                i = idim-1; il = std::max(1, idim-2); iu = idim-1;
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    dfacep = 0.5*(blend(j,k,i)+blend(j,k,iu));
                    dfacem = 0.5*(blend(j,k,i)+blend(j,k,il));
                    sigkp = dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm = dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop = dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom = dfacem*sigo1+(1.-dfacem)*sigo2;
                    voliu = voli0(j,k,3);
                    xp = si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    yp = si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    zp = si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                    volil = vol(j,k,il);
                    xm = si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    ym = si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    zm = si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                    xa = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    ya = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    za = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    ttpo = xp*xa+yp*ya+zp*za;
                    ttmo = xm*xa+ym*ya+zm*za;
                    ttpn = (xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                    ttmn = (xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                    ttp = ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm = ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                    anutm = .5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                    fnup = .5*(fnu(j,k,i+1)+fnu(j,k,i));
                    fnum = .5*(fnu(j,k,i-1)+fnu(j,k,i));
                    rhop = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(qi0(j,k,1,3)+q(j,k,i,1)));
                    rhom = (1.0-catris_kw)+catris_kw*std::sqrt(.5*(q(j,k,il,1)+q(j,k,i,1)));
                    cdp = (fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                    cdm = (fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,1) += (1.0-catris_kw)*(-bzz*turre(j,k,i-1,1)-czz*turre(j,k,i,1)-dzz*turre(j,k,i+1,1)) +
                        catris_kw*(-bzz*std::sqrt(q(j,k,il,1))*turre(j,k,i-1,1)-czz*std::sqrt(q(j,k,i,1))*turre(j,k,i,1)-dzz*std::sqrt(qi0(j,k,1,3))*turre(j,k,i+1,1));
                    rhop = (1.0-catris_kw)+catris_kw*.5*(qi0(j,k,1,3)+q(j,k,i,1));
                    rhom = (1.0-catris_kw)+catris_kw*.5*(q(j,k,il,1)+q(j,k,i,1));
                    cdp = (sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re);
                    cdm = (sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re);
                    bzz = -cdm; czz = cdp+cdm; dzz = -cdp;
                    rhside(j,k,i,2) += (1.0-catris_kw)*(-bzz*turre(j,k,i-1,2)-czz*turre(j,k,i,2)-dzz*turre(j,k,i+1,2)) +
                        catris_kw*(-bzz*q(j,k,il,1)*turre(j,k,i-1,1)-czz*q(j,k,i,1)*turre(j,k,i,1)-dzz*qi0(j,k,1,3)*turre(j,k,i+1,1));
                }
            }
            // Advective terms in zeta
            if (iturbord == 1) {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    sgnu = ccsignrc(1.,uu);
                    app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                    rhside(j,k,i,1) -= uu*(app*(turre(j,k,i,1)-turre(j,k,i-1,1))+apm*(turre(j,k,i+1,1)-turre(j,k,i,1)));
                    rhside(j,k,i,2) -= uu*(app*(turre(j,k,i,2)-turre(j,k,i-1,2))+apm*(turre(j,k,i+1,2)-turre(j,k,i,2)));
                }
            } else {
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                    yc = 0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                    zc = 0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                    tc = 0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                    uu = xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    sgnu = ccsignrc(1.,uu);
                    app = 0.5*(1.+sgnu); apm = 0.5*(1.-sgnu);
                    rhside(j,k,i,1) += -0.5*uu*app*turre(j,k,i-2,1)+2.*uu*app*turre(j,k,i-1,1)-1.5*uu*app*turre(j,k,i,1)
                        +1.5*uu*apm*turre(j,k,i,1)-2.*uu*apm*turre(j,k,i+1,1)+0.5*uu*apm*turre(j,k,i+2,1);
                    rhside(j,k,i,2) += -0.5*uu*app*turre(j,k,i-2,2)+2.*uu*app*turre(j,k,i-1,2)-1.5*uu*app*turre(j,k,i,2)
                        +1.5*uu*apm*turre(j,k,i,2)-2.*uu*apm*turre(j,k,i+1,2)+0.5*uu*apm*turre(j,k,i+2,2);
                }
            }
        } // end i2d != 1 && iaxi2planeturb != 1

        // Curvature terms for EASMCC
        if ((ivmx==8||ivmx==9||ivmx==11||ivmx==12||ivmx==13||ivmx==14) && ieasmcc2d==1) {
            sijrate2d(idim, jdim, kdim, q, qj0, qk0, bcj, bck, vol, sj, sk, vx);
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                s11 = ux(j,k,i,1) - (ux(j,k,i,1)+ux(j,k,i,9))/2.;
                s13 = 0.5*(ux(j,k,i,3)+ux(j,k,i,7));
                denom = 2.*(s11*s11+s13*s13);
                denom = ccmax(denom, xminn);
                vx(j,k,i,1) = (s11*vx(j,k,i,4)-s13*vx(j,k,i,3))/denom;
            }
        }
        // ADD SOURCE TERMS TO RHS
        if (isklton > 0) {
            if (ilamlo==0||jlamlo==0||klamlo==0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " Block #%5d in 2-eqn turb model has no laminar regions", nbl);
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " Block #%5d in 2-eqn turb model - laminar region is:", nbl);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " i=%5d to%5d, j=%5d to%5d, k=%5d to%5d", ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi);
                if (i_lam_forcezero == 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "    ...forcing vist3d=0");
                }
            }
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " NOTE:  This particular model <<transitions>> on its own, but there is");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " no guarantee that it will transition at all.  Check vist3d levels if unsure.");
        }
        // DES implementation
        if (ides==1 && (ivmx==6||ivmx==7)) {
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4));
                deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4));
                deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4));
                delta = ccmax(deltaj, deltak);
                if (i2d != 1 && iaxi2planeturb != 1) delta = ccmax(delta, deltai);
                ell = std::sqrt(turre(j,k,i,2))/(cmuc1*turre(j,k,i,1)*re);
                xlscale(j,k,i) = ccmin(ell, (double)cdes*delta);
            }
        } else if (ides >= 2 && (ivmx==6||ivmx==7)) {
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4));
                deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4));
                deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4));
                delta = ccmax(deltaj, deltak);
                if (i2d != 1 && iaxi2planeturb != 1) delta = ccmax(delta, deltai);
                ell = std::sqrt(turre(j,k,i,2))/(cmuc1*turre(j,k,i,1)*re);
                dist = ccabs(smin(j,k,i));
                velterm = ux(j,k,i,1)*ux(j,k,i,1)+ux(j,k,i,2)*ux(j,k,i,2)+ux(j,k,i,3)*ux(j,k,i,3)+
                    ux(j,k,i,4)*ux(j,k,i,4)+ux(j,k,i,5)*ux(j,k,i,5)+ux(j,k,i,6)*ux(j,k,i,6)+
                    ux(j,k,i,7)*ux(j,k,i,7)+ux(j,k,i,8)*ux(j,k,i,8)+ux(j,k,i,9)*ux(j,k,i,9);
                rd = (vist3d(j,k,i)+fnu(j,k,i))/(q(j,k,i,1)*std::sqrt(velterm)*vk*vk*dist*dist*re);
                fd = 1.0-cctanh((8.0*rd)*(8.0*rd)*(8.0*rd));
                term = ccmaxrc(0.0, ell-(double)cdes*delta);
                xlscale(j,k,i) = ell - fd*term;
                if (ides == 3) fdsav(j,k,i) = fd;
            }
        }

        // Source terms for k-omega (Wilcox & SST):
        if (ivmx == 6 || ivmx == 7) {
            pklimtermset = pklimterm;
            if (isst2003 == 1) pklimtermset = 10.;
            if (isstrc == 2) {
                sijrate3d(idim, jdim, kdim, q, ux, vol, si, sj, sk, vx);
            }
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                if ((i>=ilamlo && i<ilamhi && j>=jlamlo && j<jlamhi && k>=klamlo && k<klamhi) ||
                    (float)smin(j,k,i) < 0.) {
                    cutoff = 0.;
                } else if (ides == 3 && (float)fdsav(j,k,i) > (float)cddes) {
                    cutoff = (1.0 - fdsav(j,k,i))/(1.0-(double)cddes);
                } else {
                    cutoff = 1.;
                }
                betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                double cmuc_loc = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                double alp = blend(j,k,i)*alp1+(1.-blend(j,k,i))*alp2;
                if ((ivmx==6 && i_wilcox06==1) || (ivmx==6 && i_wilcox98==1) ||
                    ikoprod==1 || isstrc==1 || isstrc==2 || isst2003==1) {
                    s11 = ux(j,k,i,1); s22 = ux(j,k,i,5); s33 = ux(j,k,i,9);
                    s12 = 0.5*(ux(j,k,i,2)+ux(j,k,i,4));
                    s13 = 0.5*(ux(j,k,i,3)+ux(j,k,i,7));
                    s23 = 0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                    tracepart = (s11+s22+s33)*prod2d3dtrace;
                    s11t = s11-tracepart; s22t = s22-tracepart; s33t = s33-tracepart;
                    xis = s11t*s11t+s22t*s22t+s33t*s33t+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                    w12 = 0.5*(ux(j,k,i,2)-ux(j,k,i,4));
                    w13 = 0.5*(ux(j,k,i,3)-ux(j,k,i,7));
                    w23 = 0.5*(ux(j,k,i,6)-ux(j,k,i,8));
                    wis = 2.*w12*w12+2.*w13*w13+2.*w23*w23;
                }
                double fr1_loc = 1.0;
                if (isstrc == 2) {
                    double ss = std::sqrt(2.*wis); ss = ccmax(ss, xminn);
                    double xisabs = std::sqrt(2.*xis);
                    double rstar = xisabs/ss;
                    double ddd = ccmax(xisabs, 0.3*turre(j,k,i,1)*re);
                    double rtilde = 2./(ddd*ddd*ddd*ss)*
                        (-w12*vx(j,k,i,2)*(s11-s22)-w13*vx(j,k,i,3)*(s11-s33)
                         -w23*vx(j,k,i,5)*(s22-s33)
                         +s12*(-w12*(vx(j,k,i,4)-vx(j,k,i,1))-w13*vx(j,k,i,5)-w23*vx(j,k,i,3))
                         +s13*(-w13*(vx(j,k,i,6)-vx(j,k,i,1))-w12*vx(j,k,i,5)+w23*vx(j,k,i,2))
                         +s23*(-w23*(vx(j,k,i,6)-vx(j,k,i,4))+w12*vx(j,k,i,3)+w13*vx(j,k,i,2)));
                    fr1_loc = 4.*rstar/(1.+rstar)*(1.-ccacos(2.*rtilde))-1.;
                    fr1_loc = ccmincr(fr1_loc, 1.25);
                    fr1_loc = ccmaxcr(fr1_loc, 0.0);
                }
                double w06mult = 1.0;
                if (ivmx==6 && i_wilcox06==1) {
                    double wtilde = ccmax(turre(j,k,i,1), std::sqrt(2.*xis/cmuc_loc)*0.875/re);
                    w06mult = turre(j,k,i,1)/wtilde;
                }
                if (ivmx==6 && (i_wilcox06==1 || i_wilcox98==1)) {
                    double tracepart2 = 0.5*(s11+s22+s33);
                    double xnum = -(w12*w12+w13*w13)*(s11-tracepart2)
                        -(w12*w12+w23*w23)*(s22-tracepart2)
                        -(w13*w13+w23*w23)*(s33-tracepart2)-2.*w13*w23*s12
                        +2.*w12*w23*s13-2.*w12*w13*s23;
                    chiw = std::abs(xnum/std::pow(re*cmuc_loc*turre(j,k,i,1),3.0));
                    if (i_wilcox06_chiw == 1) {
                        fbeta = (1.+85.*chiw)/(1.+100.*chiw);
                    } else if (i_wilcox98_chiw == 1) {
                        fbeta = (1.+70.*chiw)/(1.+80.*chiw);
                    } else {
                        fbeta = 1.;
                    }
                    betax = betax*fbeta;
                }
                if (ivmx==6 && i_wilcox98==1) {
                    double fbetastar = (1.+680.*damp1(j,k,i)*damp1(j,k,i))/(1.+400.*damp1(j,k,i)*damp1(j,k,i));
                    cmuc_loc = cmuc_loc*fbetastar;
                }
                if (i_compress_correct == 1) {
                    double turb_mach2 = 2.*q(j,k,i,1)*turre(j,k,i,2)/(gamma*q(j,k,i,5));
                    double f_turb_mach = ccmaxcr(turb_mach2-turb_mach_0*turb_mach_0, 0.0);
                    betax = betax-cmuc_loc*xsi_star*f_turb_mach;
                    cmuc_loc = cmuc_loc*(1.+xsi_star*f_turb_mach);
                } else if (i_compress_correct == 2) {
                    double turb_mach = std::sqrt(2.*q(j,k,i,1)*turre(j,k,i,2)/(gamma*q(j,k,i,5)));
                    double zeman_lag = ccmaxcr(turb_mach-turb_mach_0, 0.0);
                    double f_turb_mach = 1.0-(std::exp(-std::pow(zeman_lag/gam_comp,2.0)));
                    betax = betax-cmuc_loc*xsi_star*f_turb_mach;
                    cmuc_loc = cmuc_loc*(1.+xsi_star*f_turb_mach);
                }
                double f4 = 1.0;
                if (isstrc == 1) {
                    double sij = std::sqrt(2.0*xis)+1.e-20;
                    double ri = (vor(j,k,i)/sij)*(vor(j,k,i)/sij-1.0);
                    f4 = 1.0/(1.0+(double)sstrc_crc*ri);
                }
                double f5 = 1.0;
                if (isstsf == 1) {
                    double povere = vor(j,k,i)*vor(j,k,i)*vist3d(j,k,i)/(0.09*q(j,k,i,1)*turre(j,k,i,1)*turre(j,k,i,2)*re*re);
                    f5 = 4.*povere-5.0;
                    f5 = ccmaxrc(1.0, f5);
                    f5 = ccminrc(12.0, f5);
                    velterm = ux(j,k,i,1)*ux(j,k,i,1)+ux(j,k,i,2)*ux(j,k,i,2)+ux(j,k,i,3)*ux(j,k,i,3)+
                        ux(j,k,i,4)*ux(j,k,i,4)+ux(j,k,i,5)*ux(j,k,i,5)+ux(j,k,i,6)*ux(j,k,i,6)+
                        ux(j,k,i,7)*ux(j,k,i,7)+ux(j,k,i,8)*ux(j,k,i,8)+ux(j,k,i,9)*ux(j,k,i,9);
                    rd = (vist3d(j,k,i)+fnu(j,k,i))/(q(j,k,i,1)*std::sqrt(velterm)*0.41*0.41*smin(j,k,i)*smin(j,k,i)*re);
                    fd = 1.0-cctanh((8.0*rd)*(8.0*rd)*(8.0*rd));
                    f5 = f5*fd+(1.0-fd);
                }
                if (ikoprod==1 || isst2003==1) {
                    pk = vist3d(j,k,i)/(q(j,k,i,1)*re)*2.*xis;
                    if (i_turbprod_kterm == 1) pk = pk-(2./3.*turre(j,k,i,2)*(s11+s22+s33));
                } else {
                    pk = vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)*vor(j,k,i);
                }
                pk = pk*fr1_loc;
                if (ides != 0) {
                    dk = std::pow(turre(j,k,i,2),1.5)/xlscale(j,k,i);
                } else {
                    dk = re*cmuc_loc*turre(j,k,i,1)*turre(j,k,i,2);
                }
                pklim = ccmin(pk, pklimtermset*dk);
                if (isst2003 == 1) {
                    pw = alp*pklim*q(j,k,i,1)/vist3d(j,k,i)*w06mult;
                } else {
                    pw = alp*pk*q(j,k,i,1)/vist3d(j,k,i)*w06mult;
                }
                dw = f5*f4*re*betax*turre(j,k,i,1)*turre(j,k,i,1);
                rhside(j,k,i,1) = rhside(j,k,i,1)+cutoff*pw-dw+(1.-blend(j,k,i))*damp1(j,k,i);
                rhside(j,k,i,2) = rhside(j,k,i,2)+cutoff*pklim-dk;
                rhside(j,k,i,1) = rhside(j,k,i,1)+keepambient*re*betax*(double)tur10[0]*(double)tur10[0];
                rhside(j,k,i,2) = rhside(j,k,i,2)+keepambient*re*cmuc_loc*(double)tur10[0]*(double)tur10[1];
            }
        } // end ivmx==6||ivmx==7

        // Source terms for nonlinear EASM k-omega (ivmx==12):
        else if (ivmx == 12) {
            alpa1 = (2.-c4)/2.*gg;
            alpa2 = (2.-c3)*gg;
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                if ((i>=ilamlo && i<ilamhi && j>=jlamlo && j<jlamhi && k>=klamlo && k<klamhi) ||
                    (float)smin(j,k,i) < 0.) { cutoff = 0.; } else { cutoff = 1.; }
                s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
                s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                w12=0.5*(ux(j,k,i,2)-ux(j,k,i,4)); w13=0.5*(ux(j,k,i,3)-ux(j,k,i,7)); w23=0.5*(ux(j,k,i,6)-ux(j,k,i,8));
                xis=s11*s11+s22*s22+s33*s33+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                wis=2.*w12*w12+2.*w13*w13+2.*w23*w23;
                omegatemp = ccmax(turre(j,k,i,1), (double)tur10[0]);
                eta=(2.-c3)*(2.-c3)*(gg*gg/4.)*xis/(omegatemp*re)/(omegatemp*re);
                squig=(2.-c4)*(2.-c4)*(gg*gg/4.)*wis/(omegatemp*re)/(omegatemp*re);
                eta=ccmincr(eta,10.); squig=ccmincr(squig,10.);
                factre=(3.*(1.+eta)+.2e-8*(eta*eta*eta+squig*squig*squig))/(3.*(1.+eta)+.2*(eta*eta*eta+squig*squig*squig));
                t11=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s11-0.33333*(s11+s22+s33))
                    -4.*alpa1*vist3d(j,k,i)*factre/omegatemp*(-s12*w12-s13*w13)/re
                    +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                t11=ccmaxcr(t11,0.);
                t22=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s22-0.33333*(s11+s22+s33))
                    -4.*alpa1*vist3d(j,k,i)*factre/omegatemp*(s12*w12-s23*w23)/re
                    +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                t22=ccmaxcr(t22,0.);
                t33=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s33-0.33333*(s11+s22+s33))
                    -4.*alpa1*vist3d(j,k,i)*factre/omegatemp*(s23*w23+s13*w13)/re
                    +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                t33=ccmaxcr(t33,0.);
                t12=-2.*vist3d(j,k,i)*s12-2.*alpa1*vist3d(j,k,i)*factre/omegatemp*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                    +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*(s11*s12+s12*s22+s13*s23)/re;
                t13=-2.*vist3d(j,k,i)*s13-2.*alpa1*vist3d(j,k,i)*factre/omegatemp*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                    +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*(s11*s13+s12*s23+s13*s33)/re;
                t23=-2.*vist3d(j,k,i)*s23-2.*alpa1*vist3d(j,k,i)*factre/omegatemp*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                    +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*(s12*s13+s22*s23+s23*s33)/re;
                pk=-(t11*ux(j,k,i,1)+t22*ux(j,k,i,5)+t33*ux(j,k,i,9)+t12*(ux(j,k,i,2)+ux(j,k,i,4))
                    +t13*(ux(j,k,i,3)+ux(j,k,i,7))+t23*(ux(j,k,i,6)+ux(j,k,i,8)))/(q(j,k,i,1)*re);
                pk=ccabs(pk);
                rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*pk*turre(j,k,i,1)/turre(j,k,i,2)-re*beta1*turre(j,k,i,1)*turre(j,k,i,1);
                dk=re*turre(j,k,i,1)*turre(j,k,i,2);
                pk=ccmin(pk,pklimterm*dk);
                rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk-dk;
                rhside(j,k,i,1)=rhside(j,k,i,1)+keepambient*re*beta1*(double)tur10[0]*(double)tur10[0];
                rhside(j,k,i,2)=rhside(j,k,i,2)+keepambient*re*(double)tur10[0]*(double)tur10[1];
            }
        } // end ivmx==12

        // Source terms for k-epsilon (Abid) ivmx==10:
        else if (ivmx == 10) {
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                if ((i>=ilamlo && i<ilamhi && j>=jlamlo && j<jlamhi && k>=klamlo && k<klamhi) ||
                    (float)smin(j,k,i) < 0.) { cutoff = 0.; } else { cutoff = 1.; }
                rek = q(j,k,i,1)*std::sqrt(turre(j,k,i,2))*re*ccabs(smin(j,k,i))/fnu(j,k,i);
                f2 = (1.-std::exp(-rek/12.));
                damp1(j,k,i) = 1.5*re*beta1*turre(j,k,i,1)*f2/turre(j,k,i,2);
                if (ikoprod == 1) {
                    s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
                    s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                    tracepart=(s11+s22+s33)*prod2d3dtrace;
                    s11t=s11-tracepart; s22t=s22-tracepart; s33t=s33-tracepart;
                    xis=s11t*s11t+s22t*s22t+s33t*s33t+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                    rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*turre(j,k,i,1)*vist3d(j,k,i)*2.*xis/(re*turre(j,k,i,2)*q(j,k,i,1))
                        -re*beta1*turre(j,k,i,1)*turre(j,k,i,1)*f2/turre(j,k,i,2);
                    pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*2.*xis;
                } else {
                    rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*turre(j,k,i,1)*vist3d(j,k,i)*vor(j,k,i)*vor(j,k,i)/(re*turre(j,k,i,2)*q(j,k,i,1))
                        -re*beta1*turre(j,k,i,1)*turre(j,k,i,1)*f2/turre(j,k,i,2);
                    pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)*vor(j,k,i);
                }
                dk=re*turre(j,k,i,1);
                pk=ccmin(pk,pklimterm*dk);
                rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk-dk;
                rhside(j,k,i,1)=rhside(j,k,i,1)+keepambient*re*beta1*(double)tur10[0]*(double)tur10[0]/(double)tur10[1];
                rhside(j,k,i,2)=rhside(j,k,i,2)+keepambient*re*(double)tur10[0];
            }
        } // end ivmx==10
        // Source terms for nonlinear EASM k-epsilon (ivmx==11):
        else if (ivmx == 11) {
            alpa1 = (2.-c4)/2.*gg;
            alpa2 = (2.-c3)*gg;
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                if ((i>=ilamlo && i<ilamhi && j>=jlamlo && j<jlamhi && k>=klamlo && k<klamhi) ||
                    (float)smin(j,k,i) < 0.) { cutoff = 0.; } else { cutoff = 1.; }
                rek=q(j,k,i,1)*std::sqrt(turre(j,k,i,2))*re*ccabs(smin(j,k,i))/fnu(j,k,i);
                f2=(1.-std::exp(-rek/12.));
                damp1(j,k,i)=1.5*re*beta1*turre(j,k,i,1)*f2/turre(j,k,i,2);
                s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
                s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                w12=0.5*(ux(j,k,i,2)-ux(j,k,i,4)); w13=0.5*(ux(j,k,i,3)-ux(j,k,i,7)); w23=0.5*(ux(j,k,i,6)-ux(j,k,i,8));
                xis=s11*s11+s22*s22+s33*s33+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                wis=2.*w12*w12+2.*w13*w13+2.*w23*w23;
                eta=(2.-c3)*(2.-c3)*(gg*gg/4.)*xis*turre(j,k,i,2)*turre(j,k,i,2)/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
                squig=(2.-c4)*(2.-c4)*(gg*gg/4.)*wis*turre(j,k,i,2)*turre(j,k,i,2)/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
                eta=ccmincr(eta,10.); squig=ccmincr(squig,10.);
                factre=(3.*(1.+eta)+.2e-8*(eta*eta*eta+squig*squig*squig))/(3.*(1.+eta)+.2*(eta*eta*eta+squig*squig*squig));
                t11=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s11-0.33333*(s11+s22+s33))
                    -4.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(-s12*w12-s13*w13)/re
                    +2.*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                t11=ccmaxcr(t11,0.);
                t22=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s22-0.33333*(s11+s22+s33))
                    -4.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s12*w12-s23*w23)/re
                    +2.*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                t22=ccmaxcr(t22,0.);
                t33=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s33-0.33333*(s11+s22+s33))
                    -4.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s23*w23+s13*w13)/re
                    +2.*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                t33=ccmaxcr(t33,0.);
                t12=-2.*vist3d(j,k,i)*s12-2.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                    +2.*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s12+s12*s22+s13*s23)/re;
                t13=-2.*vist3d(j,k,i)*s13-2.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                    +2.*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s13+s12*s23+s13*s33)/re;
                t23=-2.*vist3d(j,k,i)*s23-2.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                    +2.*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s12*s13+s22*s23+s23*s33)/re;
                pk=-(t11*ux(j,k,i,1)+t22*ux(j,k,i,5)+t33*ux(j,k,i,9)+t12*(ux(j,k,i,2)+ux(j,k,i,4))
                    +t13*(ux(j,k,i,3)+ux(j,k,i,7))+t23*(ux(j,k,i,6)+ux(j,k,i,8)))/(q(j,k,i,1)*re);
                pk=ccabs(pk);
                rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*pk*turre(j,k,i,1)/turre(j,k,i,2)-re*beta1*turre(j,k,i,1)*turre(j,k,i,1)*f2/turre(j,k,i,2);
                dk=re*turre(j,k,i,1);
                pk=ccmin(pk,pklimterm*dk);
                rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk-dk;
                rhside(j,k,i,1)=rhside(j,k,i,1)+keepambient*re*beta1*(double)tur10[0]*(double)tur10[0]/(double)tur10[1];
                rhside(j,k,i,2)=rhside(j,k,i,2)+keepambient*re*(double)tur10[0];
            }
        } // end ivmx==11

        // Source terms for EASM k-epsilon (ivmx==9 or ivmx==13):
        else if (ivmx == 9 || ivmx == 13) {
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                if ((i>=ilamlo && i<ilamhi && j>=jlamlo && j<jlamhi && k>=klamlo && k<klamhi) ||
                    (float)smin(j,k,i) < 0.) { cutoff = 0.; } else { cutoff = 1.; }
                rek=q(j,k,i,1)*std::sqrt(turre(j,k,i,2))*re*ccabs(smin(j,k,i))/fnu(j,k,i);
                f2=(1.-std::exp(-rek/10.8));
                damp1(j,k,i)=1.5*re*beta1*turre(j,k,i,1)*f2/turre(j,k,i,2);
                s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
                s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                xis=s11*s11+s22*s22+s33*s33+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                if (iturbprod == 0) {
                    pk=2.*vist3d(j,k,i)*xis/(q(j,k,i,1)*re);
                } else {
                    al10=0.5*c10-1.; al1=2.*(0.5*c11+1.);
                    if (ieasm_type == 0) { al10=al10+1.8864; al1=al1-2.; }
                    al2=0.5*c2-2./3.; al3=0.5*c3-1.; al4=0.5*c4-1.;
                    w12=0.5*(ux(j,k,i,2)-ux(j,k,i,4)); w13=0.5*(ux(j,k,i,3)-ux(j,k,i,7)); w23=0.5*(ux(j,k,i,6)-ux(j,k,i,8));
                    eta1=xis*turre(j,k,i,2)*turre(j,k,i,2)/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
                    alpa1=-al4/(al10-eta1*al1*cmuv(j,k,i));
                    alpa2=-2.*al3/(al10-eta1*cmuv(j,k,i)*al1);
                    if (ivmx == 9) { alpa1=0.; alpa2=0.; }
                    t11=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s11-0.33333*(s11+s22+s33))
                        -4.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(-s12*w12-s13*w13)/re
                        +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                    t11=ccmaxcr(t11,0.);
                    t22=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s22-0.33333*(s11+s22+s33))
                        -4.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s12*w12-s23*w23)/re
                        +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                    t22=ccmaxcr(t22,0.);
                    t33=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s33-0.33333*(s11+s22+s33))
                        -4.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s23*w23+s13*w13)/re
                        +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                    t33=ccmaxcr(t33,0.);
                    t12=-2.*vist3d(j,k,i)*s12-2.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                        +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s12+s12*s22+s13*s23)/re;
                    t13=-2.*vist3d(j,k,i)*s13-2.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                        +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s13+s12*s23+s13*s33)/re;
                    t23=-2.*vist3d(j,k,i)*s23-2.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                        +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s12*s13+s22*s23+s23*s33)/re;
                    pk=-(t11*ux(j,k,i,1)+t22*ux(j,k,i,5)+t33*ux(j,k,i,9)+t12*(ux(j,k,i,2)+ux(j,k,i,4))
                        +t13*(ux(j,k,i,3)+ux(j,k,i,7))+t23*(ux(j,k,i,6)+ux(j,k,i,8)))/(q(j,k,i,1)*re);
                    pk=ccabs(pk);
                }
                rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*pk*turre(j,k,i,1)/turre(j,k,i,2)-re*beta1*turre(j,k,i,1)*turre(j,k,i,1)*f2/turre(j,k,i,2);
                dk=re*turre(j,k,i,1);
                pk=ccmin(pk,pklimterm*dk);
                rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk-dk;
                rhside(j,k,i,1)=rhside(j,k,i,1)+keepambient*re*beta1*(double)tur10[0]*(double)tur10[0]/(double)tur10[1];
                rhside(j,k,i,2)=rhside(j,k,i,2)+keepambient*re*(double)tur10[0];
            }
        } // end ivmx==9||ivmx==13

        // Source terms for EASM k-omega (ivmx==8 or ivmx==14):
        else if (ivmx == 8 || ivmx == 14) {
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                if ((i>=ilamlo && i<ilamhi && j>=jlamlo && j<jlamhi && k>=klamlo && k<klamhi) ||
                    (float)smin(j,k,i) < 0.) { cutoff = 0.; } else { cutoff = 1.; }
                s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
                s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                xis=s11*s11+s22*s22+s33*s33+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                if (iturbprod == 0) {
                    pk=2.*vist3d(j,k,i)*xis/(q(j,k,i,1)*re);
                } else {
                    al10=0.5*c10-1.; al1=2.*(0.5*c11+1.);
                    if (ieasm_type==0||ieasm_type==3||ieasm_type==4) { al10=al10+1.8864; al1=al1-2.; }
                    al2=0.5*c2-2./3.; al3=0.5*c3-1.; al4=0.5*c4-1.;
                    w12=0.5*(ux(j,k,i,2)-ux(j,k,i,4)); w13=0.5*(ux(j,k,i,3)-ux(j,k,i,7)); w23=0.5*(ux(j,k,i,6)-ux(j,k,i,8));
                    if (ieasm_type == 4) {
                        wis=2.*w12*w12+2.*w13*w13+2.*w23*w23;
                        eta1_girimaji=xis/(xis+wis);
                        if ((float)eta1_girimaji < 0.5) {
                            c4new=2.0-((2.0-c4)*std::pow(eta1_girimaji/(1.-eta1_girimaji),0.75));
                        } else { c4new=c4; }
                        al4=0.5*c4new-1.;
                    }
                    eta1=xis/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
                    if (idurbinlim != 0 && (ieasm_type==3||ieasm_type==4)) {
                        tau=1./turre(j,k,i,1);
                        taulim=6.*std::sqrt(fnu(j,k,i)/(q(j,k,i,1)*turre(j,k,i,1)*turre(j,k,i,2)));
                        tau=ccmax(tau,taulim);
                        eta1=xis*(tau/re)*(tau/re);
                    }
                    alpa1=-al4/(al10-eta1*al1*cmuv(j,k,i));
                    alpa2=-2.*al3/(al10-eta1*cmuv(j,k,i)*al1);
                    if (ivmx == 8) { alpa1=0.; alpa2=0.; }
                    t11=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s11-0.33333*(s11+s22+s33))
                        -4.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(-s12*w12-s13*w13)/re
                        +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s11*s11+s12*s12+s13*s13-0.33333*xis)/re;
                    t11=ccmaxcr(t11,0.);
                    t22=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s22-0.33333*(s11+s22+s33))
                        -4.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(s12*w12-s23*w23)/re
                        +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s22*s22+s12*s12+s23*s23-0.33333*xis)/re;
                    t22=ccmaxcr(t22,0.);
                    t33=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*vist3d(j,k,i)*(s33-0.33333*(s11+s22+s33))
                        -4.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(s23*w23+s13*w13)/re
                        +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s33*s33+s23*s23+s13*s13-0.33333*xis)/re;
                    t33=ccmaxcr(t33,0.);
                    t12=-2.*vist3d(j,k,i)*s12-2.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(s11*w12-s22*w12-s13*w23-s23*w13)/re
                        +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s11*s12+s12*s22+s13*s23)/re;
                    t13=-2.*vist3d(j,k,i)*s13-2.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(s11*w13+s12*w23-s23*w12-s33*w13)/re
                        +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s11*s13+s12*s23+s13*s33)/re;
                    t23=-2.*vist3d(j,k,i)*s23-2.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(s12*w13+s22*w23+s13*w12-s33*w23)/re
                        +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s12*s13+s22*s23+s23*s33)/re;
                    pk=-(t11*ux(j,k,i,1)+t22*ux(j,k,i,5)+t33*ux(j,k,i,9)+t12*(ux(j,k,i,2)+ux(j,k,i,4))
                        +t13*(ux(j,k,i,3)+ux(j,k,i,7))+t23*(ux(j,k,i,6)+ux(j,k,i,8)))/(q(j,k,i,1)*re);
                    pk=ccabs(pk);
                }
                fbeta=1.;
                rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*pk*turre(j,k,i,1)/turre(j,k,i,2)-fbeta*re*beta1*turre(j,k,i,1)*turre(j,k,i,1);
                dk=damp1(j,k,i)*re*turre(j,k,i,1)*turre(j,k,i,2);
                pk=ccmin(pk,pklimterm*dk);
                rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk-dk;
                rhside(j,k,i,1)=rhside(j,k,i,1)+keepambient*re*beta1*(double)tur10[0]*(double)tur10[0];
                rhside(j,k,i,2)=rhside(j,k,i,2)+keepambient*re*(double)tur10[0]*(double)tur10[1];
            }
        } // end ivmx==8||ivmx==14

        // Source terms for k-enstrophy (ivmx==15):
        else if (ivmx == 15) {
            cmu_ke = 0.09; rcmu = 1./cmu_ke; alpha3 = 0.35;
            beta4 = 0.42; beta5 = 2.37; beta6 = 0.10; beta7 = 1.50;
            if (ibeta8kzeta != 0) beta8 = 2.3; else beta8 = 0.0;
            sigr = 0.07; c1_ke = 0.6; czeta1 = 2.1; ck_ke = 2.0;
            sigp = 0.13; sigrho = 91.9; fnuef = 198.6/(double)tinf;
            delta_ke = 0.1;
            wt1 = 1.0; wt2 = 1.0-wt1;
            wlim = 1.e-6; slim = 1.e-6;
            for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                w1(j,k,i) = ux(j,k,i,8)-ux(j,k,i,6);
                w2(j,k,i) = ux(j,k,i,3)-ux(j,k,i,7);
                w3(j,k,i) = ux(j,k,i,4)-ux(j,k,i,2);
            }
            for (i = 1; i <= idim-1; i++) {
                int ip_i = 1, im_i = 1;
                if (i == 1) im_i = 0;
                if (i == idim-1) ip_i = 0;
                double ti_val;
                if (idim-1 == 1) ti_val = 1.; else ti_val = 1./(double)(ip_i+im_i);
                for (k = 1; k <= kdim-1; k++) {
                    int kp_i = 1, km_i = 1;
                    if (k == 1) km_i = 0;
                    if (k == kdim-1) kp_i = 0;
                    double tk_val;
                    if (kdim-1 == 1) tk_val = 1.; else tk_val = 1./(double)(kp_i+km_i);
                    // compute gradients of vorticity and turbulence variables
                    j = 1;
                    dkdj(j)=(turre(j+1,k,i,2)-turre(j,k,i,2));
                    dzdj(j)=(turre(j+1,k,i,1)-turre(j,k,i,1));
                    dw1dj(j)=(w1(j+1,k,i)-w1(j,k,i));
                    dw2dj(j)=(w2(j+1,k,i)-w2(j,k,i));
                    dw3dj(j)=(w3(j+1,k,i)-w3(j,k,i));
                    dpdj(j)=(q(j+1,k,i,5)-q(j,k,i,5));
                    drdj(j)=(q(j+1,k,i,1)-q(j,k,i,1));
                    for (j = 2; j <= jdim-2; j++) {
                        dkdj(j)=0.5*(turre(j+1,k,i,2)-turre(j-1,k,i,2));
                        dzdj(j)=0.5*(turre(j+1,k,i,1)-turre(j-1,k,i,1));
                        dw1dj(j)=0.5*(w1(j+1,k,i)-w1(j-1,k,i));
                        dw2dj(j)=0.5*(w2(j+1,k,i)-w2(j-1,k,i));
                        dw3dj(j)=0.5*(w3(j+1,k,i)-w3(j-1,k,i));
                        dpdj(j)=0.5*(q(j+1,k,i,5)-q(j-1,k,i,5));
                        drdj(j)=0.5*(q(j+1,k,i,1)-q(j-1,k,i,1));
                    }
                    j = jdim-1;
                    dkdj(j)=(turre(j,k,i,2)-turre(j-1,k,i,2));
                    dzdj(j)=(turre(j,k,i,1)-turre(j-1,k,i,1));
                    dw1dj(j)=(w1(j,k,i)-w1(j-1,k,i));
                    dw2dj(j)=(w2(j,k,i)-w2(j-1,k,i));
                    dw3dj(j)=(w3(j,k,i)-w3(j-1,k,i));
                    dpdj(j)=(q(j,k,i,5)-q(j-1,k,i,5));
                    drdj(j)=(q(j,k,i,1)-q(j-1,k,i,1));
                    for (j = 1; j <= jdim-1; j++) {
                        dkdk(j)=tk_val*(turre(j,k+kp_i,i,2)-turre(j,k-km_i,i,2));
                        dkdi(j)=ti_val*(turre(j,k,i+ip_i,2)-turre(j,k,i-im_i,2));
                        dzdk(j)=tk_val*(turre(j,k+kp_i,i,1)-turre(j,k-km_i,i,1));
                        dzdi(j)=ti_val*(turre(j,k,i+ip_i,1)-turre(j,k,i-im_i,1));
                        dw1dk(j)=tk_val*(w1(j,k+kp_i,i)-w1(j,k-km_i,i));
                        dw1di(j)=ti_val*(w1(j,k,i+ip_i)-w1(j,k,i-im_i));
                        dw2dk(j)=tk_val*(w2(j,k+kp_i,i)-w2(j,k-km_i,i));
                        dw2di(j)=ti_val*(w2(j,k,i+ip_i)-w2(j,k,i-im_i));
                        dw3dk(j)=tk_val*(w3(j,k+kp_i,i)-w3(j,k-km_i,i));
                        dw3di(j)=ti_val*(w3(j,k,i+ip_i)-w3(j,k,i-im_i));
                        dpdk(j)=tk_val*(q(j,k+kp_i,i,5)-q(j,k-km_i,i,5));
                        dpdi(j)=ti_val*(q(j,k,i+ip_i,5)-q(j,k,i-im_i,5));
                        drdk(j)=tk_val*(q(j,k+kp_i,i,1)-q(j,k-km_i,i,1));
                        drdi(j)=ti_val*(q(j,k,i+ip_i,1)-q(j,k,i-im_i,1));
                    }
                    for (j = 1; j <= jdim-1; j++) {
                        rvol = 1.0/vol(j,k,i);
                        sj4=sj(j,k,i,4); sj4p=sj(j+1,k,i,4);
                        xjp=0.5*(sj(j,k,i,1)*sj4+sj(j+1,k,i,1)*sj4p);
                        yjp=0.5*(sj(j,k,i,2)*sj4+sj(j+1,k,i,2)*sj4p);
                        zjp=0.5*(sj(j,k,i,3)*sj4+sj(j+1,k,i,3)*sj4p);
                        sk4=sk(j,k,i,4); sk4p=sk(j,k+1,i,4);
                        xkp=0.5*(sk(j,k,i,1)*sk4+sk(j,k+1,i,1)*sk4p);
                        ykp=0.5*(sk(j,k,i,2)*sk4+sk(j,k+1,i,2)*sk4p);
                        zkp=0.5*(sk(j,k,i,3)*sk4+sk(j,k+1,i,3)*sk4p);
                        si4=si(j,k,i,4); si4p=si(j,k,i+1,4);
                        xip=0.5*(si(j,k,i,1)*si4+si(j,k,i+1,1)*si4p);
                        yip=0.5*(si(j,k,i,2)*si4+si(j,k,i+1,2)*si4p);
                        zip=0.5*(si(j,k,i,3)*si4+si(j,k,i+1,3)*si4p);
                        dkdx(j)=(xjp*dkdj(j)+xkp*dkdk(j)+xip*dkdi(j))*rvol;
                        dkdy(j)=(yjp*dkdj(j)+ykp*dkdk(j)+yip*dkdi(j))*rvol;
                        dkdz(j)=(zjp*dkdj(j)+zkp*dkdk(j)+zip*dkdi(j))*rvol;
                        dzdx(j)=(xjp*dzdj(j)+xkp*dzdk(j)+xip*dzdi(j))*rvol;
                        dzdy(j)=(yjp*dzdj(j)+ykp*dzdk(j)+yip*dzdi(j))*rvol;
                        dzdz(j)=(zjp*dzdj(j)+zkp*dzdk(j)+zip*dzdi(j))*rvol;
                        dw1dx(j)=(xjp*dw1dj(j)+xkp*dw1dk(j)+xip*dw1di(j))*rvol;
                        dw1dy(j)=(yjp*dw1dj(j)+ykp*dw1dk(j)+yip*dw1di(j))*rvol;
                        dw1dz(j)=(zjp*dw1dj(j)+zkp*dw1dk(j)+zip*dw1di(j))*rvol;
                        dw2dx(j)=(xjp*dw2dj(j)+xkp*dw2dk(j)+xip*dw2di(j))*rvol;
                        dw2dy(j)=(yjp*dw2dj(j)+ykp*dw2dk(j)+yip*dw2di(j))*rvol;
                        dw2dz(j)=(zjp*dw2dj(j)+zkp*dw2dk(j)+zip*dw2di(j))*rvol;
                        dw3dx(j)=(xjp*dw3dj(j)+xkp*dw3dk(j)+xip*dw3di(j))*rvol;
                        dw3dy(j)=(yjp*dw3dj(j)+ykp*dw3dk(j)+yip*dw3di(j))*rvol;
                        dw3dz(j)=(zjp*dw3dj(j)+zkp*dw3dk(j)+zip*dw3di(j))*rvol;
                        dpdx(j)=(xjp*dpdj(j)+xkp*dpdk(j)+xip*dpdi(j))*rvol;
                        dpdy(j)=(yjp*dpdj(j)+ykp*dpdk(j)+yip*dpdi(j))*rvol;
                        dpdz(j)=(zjp*dpdj(j)+zkp*dpdk(j)+zip*dpdi(j))*rvol;
                        drdx(j)=(xjp*drdj(j)+xkp*drdk(j)+xip*drdi(j))*rvol;
                        drdy(j)=(yjp*drdj(j)+ykp*drdk(j)+yip*drdi(j))*rvol;
                        drdz(j)=(zjp*drdj(j)+zkp*drdk(j)+zip*drdi(j))*rvol;
                    }
                    for (j = 1; j <= jdim-1; j++) {
                        fnut=vist3d(j,k,i);
                        s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
                        divv=s11+s22+s33;
                        s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                        wis=w1(j,k,i)*w1(j,k,i)+w2(j,k,i)*w2(j,k,i)+w3(j,k,i)*w3(j,k,i);
                        omemag=std::sqrt(wis);
                        xis=s11*s11+s22*s22+s33*s33+2.*s12*s12+2.*s13*s13+2.*s23*s23+slim;
                        t11=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*fnut*(s11-0.33333*divv); t11=ccmaxcr(t11,0.);
                        t22=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*fnut*(s22-0.33333*divv); t22=ccmaxcr(t22,0.);
                        t33=2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.-2.*fnut*(s33-0.33333*divv); t33=ccmaxcr(t33,0.);
                        t12=-2.*fnut*s12; t13=-2.*fnut*s13; t23=-2.*fnut*s23;
                        p3=2.0*fnut*alpha3*turre(j,k,i,1)*(xis-slim)/turre(j,k,i,2)+2.0*q(j,k,i,1)*re*turre(j,k,i,1)*divv/3.;
                        p8=-2.0*beta8/turre(j,k,i,2)/(xis+0.5*wis)*
                            ((w1(j,k,i)*t13+w2(j,k,i)*t23+w3(j,k,i)*t33)*(dkdx(j)*dzdy(j)-dkdy(j)*dzdx(j))+
                             (w1(j,k,i)*t11+w2(j,k,i)*t12+w3(j,k,i)*t13)*(dkdy(j)*dzdz(j)-dkdz(j)*dzdy(j))+
                             (w1(j,k,i)*t12+w2(j,k,i)*t22+w3(j,k,i)*t23)*(dkdz(j)*dzdx(j)-dkdx(j)*dzdz(j)));
                        wiwjtauij=w1(j,k,i)*w1(j,k,i)*t11+w2(j,k,i)*w2(j,k,i)*t22+w3(j,k,i)*w3(j,k,i)*t33
                            +2.*(w1(j,k,i)*w2(j,k,i)*t12+w1(j,k,i)*w3(j,k,i)*t13+w2(j,k,i)*w3(j,k,i)*t23);
                        p4p6=(beta4*turre(j,k,i,1)/(omemag+wlim)+2.0*beta6*fnut*omemag/(fnu(j,k,i)*re*re))*wiwjtauij/turre(j,k,i,2);
                        p7=beta7*re*q(j,k,i,1)*turre(j,k,i,1)*(w1(j,k,i)*w1(j,k,i)*s11+w2(j,k,i)*w2(j,k,i)*s22+w3(j,k,i)*w3(j,k,i)*s33
                            +2.*(w1(j,k,i)*w2(j,k,i)*s12+w1(j,k,i)*w3(j,k,i)*s13+w2(j,k,i)*w3(j,k,i)*s23))/(wis+wlim);
                        p5=fnut/(sigr*re*re)*(2.*(dw1dx(j)*dw1dx(j)+dw2dy(j)*dw2dy(j)+dw3dz(j)*dw3dz(j))+
                            (dw2dx(j)+dw1dy(j))*(dw2dx(j)+dw1dy(j))+(dw3dx(j)+dw1dz(j))*(dw3dx(j)+dw1dz(j))+
                            (dw3dy(j)+dw2dz(j))*(dw3dy(j)+dw2dz(j)));
                        xnu=fnu(j,k,i)/q(j,k,i,1);
                        rt=(turre(j,k,i,2)/xnu)*(turre(j,k,i,2)/xnu)/turre(j,k,i,1);
                        ft=1.0; // ft=1.0 (overrides the formula above)
                        d5=re*beta5*ft*std::pow(turre(j,k,i,1),1.5)/(std::sqrt(rt)+delta_ke);
                        ct1=-fnut*(drdx(j)*dpdx(j)+drdy(j)*dpdy(j)+drdz(j)*dpdz(j))/(re*ck_ke*q(j,k,i,1)*q(j,k,i,1)*q(j,k,i,1));
                        rtaurho=std::sqrt(turre(j,k,i,2)*(drdx(j)*drdx(j)+drdy(j)*drdy(j)+drdz(j)*drdz(j)))/q(j,k,i,1);
                        ct2=-c1_ke*turre(j,k,i,2)*rtaurho;
                        DpDt=q(j,k,i,2)*dpdx(j)+q(j,k,i,3)*dpdy(j)+q(j,k,i,4)*dpdz(j);
                        termn=(q(j,k,i,1)*turre(j,k,i,2)*omemag)/(xnu*q(j,k,i,5)*sigp);
                        termd=1.+sigrho*rtaurho*std::sqrt(rt/turre(j,k,i,1))/re;
                        ct3=ccmaxcr(termn*DpDt/termd,0.);
                        ct4=-czeta1*fnut*omemag*turre(j,k,i,1)*rtaurho/turre(j,k,i,2);
                        ct5=-2.0*re*q(j,k,i,1)*turre(j,k,i,1)*divv;
                        if (ikoprod == 2) {
                            pk=-(t11*ux(j,k,i,1)+t22*ux(j,k,i,5)+t33*ux(j,k,i,9)+t12*(ux(j,k,i,2)+ux(j,k,i,4))
                                +t13*(ux(j,k,i,3)+ux(j,k,i,7))+t23*(ux(j,k,i,6)+ux(j,k,i,8)))/(q(j,k,i,1)*re);
                        } else if (ikoprod == 1) {
                            pk=2.*vist3d(j,k,i)*xis/(q(j,k,i,1)*re);
                        } else {
                            pk=vist3d(j,k,i)*vor(j,k,i)*vor(j,k,i)/(q(j,k,i,1)*re);
                        }
                        rtauk=re*fnu(j,k,i)*turre(j,k,i,1)/(turre(j,k,i,2)*q(j,k,i,1));
                        dk=turre(j,k,i,2)*rtauk;
                        pklim=ccmin(pk,pklimterm*dk);
                        pk=pklim;
                        srce(j,k,i,2)=pk+ct1+ct2-dk;
                        rhside(j,k,i,2)=rhside(j,k,i,2)+srce(j,k,i,2);
                        srce(j,k,i,1)=(p3+p4p6+p5+p7+p8+ct3+ct4+ct5)/(q(j,k,i,1)*re)-d5;
                        rhside(j,k,i,1)=rhside(j,k,i,1)+srce(j,k,i,1);
                    }
                }
            }
        } // end ivmx==15

        // Source terms for k-kL-MEAH2015 (ivmx==16):
        else if (ivmx == 16) {
            u_doubleprime(idim, jdim, kdim, q, ux, vol, si, sj, sk, vx);
            for (i = 1; i <= idim-1; i++)
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                if ((i>=ilamlo && i<ilamhi && j>=jlamlo && j<jlamhi && k>=klamlo && k<klamhi) ||
                    (float)smin(j,k,i) < 0.) { cutoff = 0.; } else { cutoff = 1.; }
                s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
                s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                tracepart=(s11+s22+s33)*prod2d3dtrace;
                s11t=s11-tracepart; s22t=s22-tracepart; s33t=s33-tracepart;
                xis=s11t*s11t+s22t*s22t+s33t*s33t+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                uprime=std::sqrt(2.*xis);
                uprimeprime=vx(j,k,i,1);
                uprimeprime=ccabs(uprimeprime);
                uprimeprime=ccmaxcr(uprimeprime,1.e-20);
                if (ikoprod == 1) {
                    pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*2.*xis;
                } else {
                    pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)*vor(j,k,i);
                }
                dk1=std::pow(cmuc1,0.75)*std::pow(turre(j,k,i,2),2.5)/turre(j,k,i,1)*re;
                dk2=2.*fnu(j,k,i)/q(j,k,i,1)*turre(j,k,i,2)/(smin(j,k,i)*smin(j,k,i))/re;
                dk=dk1+dk2;
                pk=ccmin(pk,20.*dk1);
                zlvk=vk*(uprime/uprimeprime);
                xmin_kkl=turre(j,k,i,1)/(turre(j,k,i,2)*c11_kkl)/re;
                xmax_kkl=pk*turre(j,k,i,1)/((std::pow(cmuc1,0.75))*std::pow(turre(j,k,i,2),2.5)*re);
                xmax_kkl=ccmax(xmax_kkl,0.5);
                xmax_kkl=ccmin(xmax_kkl,1.0);
                xmax_kkl=c12_kkl*vk*smin(j,k,i)*xmax_kkl;
                zlvk=ccmax(zlvk,xmin_kkl);
                zlvk=ccmin(zlvk,xmax_kkl);
                c_phi1=zeta1_kkl-(zeta2_kkl*(turre(j,k,i,1)/(turre(j,k,i,2)*zlvk*re))*(turre(j,k,i,1)/(turre(j,k,i,2)*zlvk*re)));
                squig=q(j,k,i,1)*smin(j,k,i)*std::sqrt(0.3*turre(j,k,i,2))/(20.*fnu(j,k,i))*re;
                f_phi=(1.+(cd1_kkl*squig))/(1.+squig*squig*squig*squig);
                damp1(j,k,i)=f_phi;
                pw=c_phi1*turre(j,k,i,1)/turre(j,k,i,2)*pk;
                dw=zeta3_kkl*std::pow(turre(j,k,i,2),1.5)*re+6.*f_phi*fnu(j,k,i)/q(j,k,i,1)*turre(j,k,i,1)/(smin(j,k,i)*smin(j,k,i))/re;
                rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*pw-dw;
                rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk-dk;
            }
        } // end ivmx==16

        // Implicit F_eta_eta viscous terms (LHS) - do over all i's
        for (i = 1; i <= idim-1; i++) {
            // Helper lambda for source term diagonal additions
            auto get_cyadd = [&](int j, int k) -> std::pair<double,double> {
                double cyadd_v, cy2add_v;
                if (ivmx == 6) {
                    cyadd_v = 2.*re*beta1*turre(j,k,i,1);
                    cy2add_v = re*cmuc1*turre(j,k,i,1);
                } else if (ivmx == 7) {
                    cmuc_var = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                    betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                    cyadd_v = 2.*re*betax*turre(j,k,i,1)+ccabs(damp1(j,k,i))/turre(j,k,i,1);
                    cy2add_v = re*cmuc_var*turre(j,k,i,1);
                } else if (ivmx == 12) {
                    cyadd_v = 2.*re*beta1*turre(j,k,i,1);
                    cy2add_v = re*turre(j,k,i,1);
                } else if (ivmx==8||ivmx==14) {
                    cyadd_v = 2.*re*beta1*turre(j,k,i,1);
                    cy2add_v = damp1(j,k,i)*re*turre(j,k,i,1);
                } else if (ivmx==9||ivmx==10||ivmx==11||ivmx==13) {
                    cyadd_v = damp1(j,k,i);
                    cy2add_v = 2.*re*turre(j,k,i,1)/turre(j,k,i,2);
                } else if (ivmx == 15) {
                    term1 = -0.5*(srce(j,k,i,1)-ccabs(srce(j,k,i,1)));
                    term2 = -0.5*(srce(j,k,i,2)-ccabs(srce(j,k,i,2)));
                    cyadd_v = term1/turre(j,k,i,1);
                    cy2add_v = term2/turre(j,k,i,2);
                } else if (ivmx == 16) {
                    cyadd_v = 6.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*smin(j,k,i))*damp1(j,k,i)/re;
                    cy2add_v = 2.5*std::pow(cmuc1,0.75)*std::pow(turre(j,k,i,2),1.5)/turre(j,k,i,1)*re
                        +2.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*smin(j,k,i))/re;
                } else {
                    cyadd_v = 0.; cy2add_v = 0.;
                }
                return {cyadd_v, cy2add_v};
            };
            // Interior points
            for (k = 2; k <= kdim-2; k++) {
                kl = k-1; ku = k+1;
                for (j = 1; j <= jdim-1; j++) {
                    dfacep=0.5*(blend(j,k,i)+blend(j,ku,i));
                    dfacem=0.5*(blend(j,k,i)+blend(j,kl,i));
                    sigkp=dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm=dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop=dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom=dfacem*sigo1+(1.-dfacem)*sigo2;
                    volku=vol(j,ku,i);
                    xp=sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    yp=sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    zp=sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    volkl=vol(j,kl,i);
                    xm=sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    ym=sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    zm=sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    xa=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    ya=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    za=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    ttpo=xp*xa+yp*ya+zp*za;
                    ttmo=xm*xa+ym*ya+zm*za;
                    ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                    ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                    ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                    anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                    fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i));
                    fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i));
                    cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    by(j,k)=-cdm; cy(j,k)=cdp+cdm; dy(j,k)=-cdp;
                    cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    by2(j,k)=-cdm; cy2(j,k)=cdp+cdm; dy2(j,k)=-cdp;
                }
                for (j = 1; j <= jdim-1; j++) {
                    xc=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    yc=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    zc=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    tc=0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                    uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    sgnu=ccsignrc(1.,uu); app=0.5*(1.+sgnu); apm=0.5*(1.-sgnu);
                    by(j,k)=by(j,k)-uu*app; by2(j,k)=by2(j,k)-uu*app;
                    cy(j,k)=cy(j,k)+uu*(app-apm); cy2(j,k)=cy2(j,k)+uu*(app-apm);
                    dy(j,k)=dy(j,k)+uu*apm; dy2(j,k)=dy2(j,k)+uu*apm;
                }
                for (j = 1; j <= jdim-1; j++) {
                    auto [cyadd_v, cy2add_v] = get_cyadd(j,k);
                    cy(j,k)=cy(j,k)+cyadd_v; cy2(j,k)=cy2(j,k)+cy2add_v;
                }
                for (j = 1; j <= jdim-1; j++) {
                    by(j,k)=by(j,k)*timestp(j,k,i);
                    by2(j,k)=by2(j,k)*timestp(j,k,i)*factor2;
                    cy(j,k)=cy(j,k)*timestp(j,k,i)+1.0*(1.+phi);
                    cy2(j,k)=cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                    dy(j,k)=dy(j,k)*timestp(j,k,i);
                    dy2(j,k)=dy2(j,k)*timestp(j,k,i)*factor2;
                    fy(j,k)=rhside(j,k,i,1)*timestp(j,k,i);
                    fy2(j,k)=rhside(j,k,i,2)*timestp(j,k,i)*factor2;
                }
                if ((float)dt > 0.) {
                    for (j = 1; j <= jdim-1; j++) {
                        fy(j,k)=fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-turre(j,k,i,1))+phi*zksav2(j,k,i,3);
                        fy2(j,k)=fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-turre(j,k,i,2))+phi*zksav2(j,k,i,4);
                    }
                }
            } // end interior k loop

            // K0 boundary points
            {
                k = 1; kl = 1; ku = std::min(2, kdim-1);
                for (j = 1; j <= jdim-1; j++) {
                    dfacep=0.5*(blend(j,k,i)+blend(j,ku,i));
                    dfacem=0.5*(blend(j,k,i)+blend(j,kl,i));
                    sigkp=dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm=dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop=dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom=dfacem*sigo1+(1.-dfacem)*sigo2;
                    volku=vol(j,ku,i);
                    xp=sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    yp=sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    zp=sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    volkl=volk0(j,i,1);
                    xm=sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    ym=sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    zm=sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    xa=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    ya=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    za=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    ttpo=xp*xa+yp*ya+zp*za;
                    ttmo=xm*xa+ym*ya+zm*za;
                    ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                    ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                    ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                    anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                    fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i));
                    fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i));
                    cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    by(j,k)=-cdm; cy(j,k)=cdp+cdm; dy(j,k)=-cdp;
                    cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    by2(j,k)=-cdm; cy2(j,k)=cdp+cdm; dy2(j,k)=-cdp;
                }
                for (j = 1; j <= jdim-1; j++) {
                    xc=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    yc=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    zc=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    tc=0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                    uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    sgnu=ccsignrc(1.,uu); app=0.5*(1.+sgnu); apm=0.5*(1.-sgnu);
                    by(j,k)=by(j,k)-uu*app; by2(j,k)=by2(j,k)-uu*app;
                    cy(j,k)=cy(j,k)+uu*(app-apm); cy2(j,k)=cy2(j,k)+uu*(app-apm);
                    dy(j,k)=dy(j,k)+uu*apm; dy2(j,k)=dy2(j,k)+uu*apm;
                }
                for (j = 1; j <= jdim-1; j++) {
                    auto [cyadd_v, cy2add_v] = get_cyadd(j,k);
                    cy(j,k)=cy(j,k)+cyadd_v; cy2(j,k)=cy2(j,k)+cy2add_v;
                }
                for (j = 1; j <= jdim-1; j++) {
                    by(j,k)=by(j,k)*timestp(j,k,i);
                    by2(j,k)=by2(j,k)*timestp(j,k,i)*factor2;
                    cy(j,k)=cy(j,k)*timestp(j,k,i)+1.0*(1.+phi);
                    cy2(j,k)=cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                    dy(j,k)=dy(j,k)*timestp(j,k,i);
                    dy2(j,k)=dy2(j,k)*timestp(j,k,i)*factor2;
                    fy(j,k)=rhside(j,k,i,1)*timestp(j,k,i);
                    fy2(j,k)=rhside(j,k,i,2)*timestp(j,k,i)*factor2;
                }
                if ((float)dt > 0.) {
                    for (j = 1; j <= jdim-1; j++) {
                        fy(j,k)=fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-turre(j,k,i,1))+phi*zksav2(j,k,i,3);
                        fy2(j,k)=fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-turre(j,k,i,2))+phi*zksav2(j,k,i,4);
                    }
                }
            }
            // KDIM boundary points
            {
                k = kdim-1; kl = kdim-2; ku = kdim-1;
                for (j = 1; j <= jdim-1; j++) {
                    dfacep=0.5*(blend(j,k,i)+blend(j,ku,i));
                    dfacem=0.5*(blend(j,k,i)+blend(j,kl,i));
                    sigkp=dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm=dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop=dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom=dfacem*sigo1+(1.-dfacem)*sigo2;
                    volku=volk0(j,i,3);
                    xp=sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    yp=sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    zp=sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volku));
                    volkl=vol(j,kl,i);
                    xm=sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    ym=sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    zm=sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)+volkl));
                    xa=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    ya=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    za=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    ttpo=xp*xa+yp*ya+zp*za;
                    ttmo=xm*xa+ym*ya+zm*za;
                    ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/vol(j,k,i);
                    ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/vol(j,k,i);
                    ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i));
                    anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i));
                    fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i));
                    fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i));
                    cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    by(j,k)=-cdm; cy(j,k)=cdp+cdm; dy(j,k)=-cdp;
                    cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    by2(j,k)=-cdm; cy2(j,k)=cdp+cdm; dy2(j,k)=-cdp;
                }
                for (j = 1; j <= jdim-1; j++) {
                    xc=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i);
                    yc=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i);
                    zc=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i);
                    tc=0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+sk(j,k,i,5)*sk(j,k,i,4))/vol(j,k,i);
                    uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    sgnu=ccsignrc(1.,uu); app=0.5*(1.+sgnu); apm=0.5*(1.-sgnu);
                    by(j,k)=by(j,k)-uu*app; by2(j,k)=by2(j,k)-uu*app;
                    cy(j,k)=cy(j,k)+uu*(app-apm); cy2(j,k)=cy2(j,k)+uu*(app-apm);
                    dy(j,k)=dy(j,k)+uu*apm; dy2(j,k)=dy2(j,k)+uu*apm;
                }
                for (j = 1; j <= jdim-1; j++) {
                    auto [cyadd_v, cy2add_v] = get_cyadd(j,k);
                    cy(j,k)=cy(j,k)+cyadd_v; cy2(j,k)=cy2(j,k)+cy2add_v;
                }
                for (j = 1; j <= jdim-1; j++) {
                    by(j,k)=by(j,k)*timestp(j,k,i);
                    by2(j,k)=by2(j,k)*timestp(j,k,i)*factor2;
                    cy(j,k)=cy(j,k)*timestp(j,k,i)+1.0*(1.+phi);
                    cy2(j,k)=cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                    dy(j,k)=dy(j,k)*timestp(j,k,i);
                    dy2(j,k)=dy2(j,k)*timestp(j,k,i)*factor2;
                    fy(j,k)=rhside(j,k,i,1)*timestp(j,k,i);
                    fy2(j,k)=rhside(j,k,i,2)*timestp(j,k,i)*factor2;
                }
                if ((float)dt > 0.) {
                    for (j = 1; j <= jdim-1; j++) {
                        fy(j,k)=fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-turre(j,k,i,1))+phi*zksav2(j,k,i,3);
                        fy2(j,k)=fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-turre(j,k,i,2))+phi*zksav2(j,k,i,4);
                    }
                }
            }
            // Apply blank array for overset
            if (iover == 1) {
                for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++) {
                    fy(j,k)=fy(j,k)*blank(j,k,i);
                    by(j,k)=by(j,k)*blank(j,k,i);
                    dy(j,k)=dy(j,k)*blank(j,k,i);
                    cy(j,k)=cy(j,k)*blank(j,k,i)+(1.-blank(j,k,i));
                    fy2(j,k)=fy2(j,k)*blank(j,k,i);
                    by2(j,k)=by2(j,k)*blank(j,k,i);
                    dy2(j,k)=dy2(j,k)*blank(j,k,i);
                    cy2(j,k)=cy2(j,k)*blank(j,k,i)+(1.-blank(j,k,i));
                }
            }
            // Solve tridiagonal systems
            triv(jdim-1, kdim-1, 1, jdim-1, 1, kdim-1, worky, by, cy, dy, fy);
            triv(jdim-1, kdim-1, 1, jdim-1, 1, kdim-1, worky, by2, cy2, dy2, fy2);

            // Store results back
            for (k = 1; k <= kdim-1; k++)
            for (j = 1; j <= jdim-1; j++) {
                rhside(j,k,i,1) = fy(j,k);
                rhside(j,k,i,2) = fy2(j,k);
            }
        } // end i loop for eta sweep

        // Implicit F_xi_xi viscous terms. Do over all i's
        for (i = 1; i <= idim-1; i++) {
            auto get_cxadd = [&](int j, int k) -> std::pair<double,double> {
                double cxadd_v, cx2add_v;
                if (ivmx == 6) {
                    cxadd_v = 2.*re*beta1*turre(j,k,i,1);
                    cx2add_v = re*cmuc1*turre(j,k,i,1);
                } else if (ivmx == 7) {
                    cmuc_var = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                    betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                    cxadd_v = 2.*re*betax*turre(j,k,i,1)+ccabs(damp1(j,k,i))/turre(j,k,i,1);
                    cx2add_v = re*cmuc_var*turre(j,k,i,1);
                } else if (ivmx == 12) {
                    cxadd_v = 2.*re*beta1*turre(j,k,i,1);
                    cx2add_v = re*turre(j,k,i,1);
                } else if (ivmx==8||ivmx==14) {
                    cxadd_v = 2.*re*beta1*turre(j,k,i,1);
                    cx2add_v = damp1(j,k,i)*re*turre(j,k,i,1);
                } else if (ivmx==9||ivmx==10||ivmx==11||ivmx==13) {
                    cxadd_v = damp1(j,k,i);
                    cx2add_v = 2.*re*turre(j,k,i,1)/turre(j,k,i,2);
                } else if (ivmx == 15) {
                    term1 = -0.5*(srce(j,k,i,1)-ccabs(srce(j,k,i,1)));
                    term2 = -0.5*(srce(j,k,i,2)-ccabs(srce(j,k,i,2)));
                    cxadd_v = term1/turre(j,k,i,1);
                    cx2add_v = term2/turre(j,k,i,2);
                } else if (ivmx == 16) {
                    cxadd_v = 6.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*smin(j,k,i))*damp1(j,k,i)/re;
                    cx2add_v = 2.5*std::pow(cmuc1,0.75)*std::pow(turre(j,k,i,2),1.5)/turre(j,k,i,1)*re
                        +2.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*smin(j,k,i))/re;
                } else {
                    cxadd_v = 0.; cx2add_v = 0.;
                }
                return {cxadd_v, cx2add_v};
            };
            // Interior points
            for (j = 2; j <= jdim-2; j++) {
                jl = j-1; ju = j+1;
                for (k = 1; k <= kdim-1; k++) {
                    dfacep=0.5*(blend(j,k,i)+blend(ju,k,i));
                    dfacem=0.5*(blend(j,k,i)+blend(jl,k,i));
                    sigkp=dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm=dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop=dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom=dfacem*sigo1+(1.-dfacem)*sigo2;
                    volju=vol(ju,k,i);
                    xp=sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    yp=sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    zp=sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    voljl=vol(jl,k,i);
                    xm=sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    ym=sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    zm=sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    xa=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    ya=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    za=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    ttpo=xp*xa+yp*ya+zp*za;
                    ttmo=xm*xa+ym*ya+zm*za;
                    ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                    ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                    ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                    anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                    fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i));
                    fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i));
                    cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    bx(k,j)=-cdm; cx(k,j)=cdp+cdm; dx(k,j)=-cdp;
                    cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    bx2(k,j)=-cdm; cx2(k,j)=cdp+cdm; dx2(k,j)=-cdp;
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    yc=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    zc=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    tc=0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                    uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    sgnu=ccsignrc(1.,uu); app=0.5*(1.+sgnu); apm=0.5*(1.-sgnu);
                    bx(k,j)=bx(k,j)-uu*app; bx2(k,j)=bx2(k,j)-uu*app;
                    cx(k,j)=cx(k,j)+uu*(app-apm); cx2(k,j)=cx2(k,j)+uu*(app-apm);
                    dx(k,j)=dx(k,j)+uu*apm; dx2(k,j)=dx2(k,j)+uu*apm;
                }
                // NOTE: the source-diagonal (cxadd) belongs ONLY to the first (eta)
                // ADI sweep in twoeqn.f; the xi sweep must NOT re-add it (doing so
                // makes the diagonal ~1e6x too large and kills the implicit update).
                for (k = 1; k <= kdim-1; k++) {
                    bx(k,j)=bx(k,j)*timestp(j,k,i);
                    bx2(k,j)=bx2(k,j)*timestp(j,k,i)*factor2;
                    cx(k,j)=cx(k,j)*timestp(j,k,i)+1.0*(1.+phi);
                    cx2(k,j)=cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                    dx(k,j)=dx(k,j)*timestp(j,k,i);
                    dx2(k,j)=dx2(k,j)*timestp(j,k,i)*factor2;
                    fx(k,j)=rhside(j,k,i,1)*(1.+phi);
                    fx2(k,j)=rhside(j,k,i,2)*(1.+phi);
                }
            } // end interior j loop

            // J0 boundary points
            {
                j = 1; jl = 1; ju = std::min(2, jdim-1);
                for (k = 1; k <= kdim-1; k++) {
                    dfacep=0.5*(blend(j,k,i)+blend(ju,k,i));
                    dfacem=0.5*(blend(j,k,i)+blend(jl,k,i));
                    sigkp=dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm=dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop=dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom=dfacem*sigo1+(1.-dfacem)*sigo2;
                    volju=vol(ju,k,i);
                    xp=sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    yp=sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    zp=sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    voljl=volj0(k,i,1);
                    xm=sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    ym=sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    zm=sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    xa=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    ya=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    za=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    ttpo=xp*xa+yp*ya+zp*za;
                    ttmo=xm*xa+ym*ya+zm*za;
                    ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                    ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                    ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                    anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                    fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i));
                    fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i));
                    cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    bx(k,j)=-cdm; cx(k,j)=cdp+cdm; dx(k,j)=-cdp;
                    cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    bx2(k,j)=-cdm; cx2(k,j)=cdp+cdm; dx2(k,j)=-cdp;
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    yc=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    zc=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    tc=0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                    uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    sgnu=ccsignrc(1.,uu); app=0.5*(1.+sgnu); apm=0.5*(1.-sgnu);
                    bx(k,j)=bx(k,j)-uu*app; bx2(k,j)=bx2(k,j)-uu*app;
                    cx(k,j)=cx(k,j)+uu*(app-apm); cx2(k,j)=cx2(k,j)+uu*(app-apm);
                    dx(k,j)=dx(k,j)+uu*apm; dx2(k,j)=dx2(k,j)+uu*apm;
                }
                // NOTE: the source-diagonal (cxadd) belongs ONLY to the first (eta)
                // ADI sweep in twoeqn.f; the xi sweep must NOT re-add it (doing so
                // makes the diagonal ~1e6x too large and kills the implicit update).
                for (k = 1; k <= kdim-1; k++) {
                    bx(k,j)=bx(k,j)*timestp(j,k,i);
                    bx2(k,j)=bx2(k,j)*timestp(j,k,i)*factor2;
                    cx(k,j)=cx(k,j)*timestp(j,k,i)+1.0*(1.+phi);
                    cx2(k,j)=cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                    dx(k,j)=dx(k,j)*timestp(j,k,i);
                    dx2(k,j)=dx2(k,j)*timestp(j,k,i)*factor2;
                    fx(k,j)=rhside(j,k,i,1)*(1.+phi);
                    fx2(k,j)=rhside(j,k,i,2)*(1.+phi);
                }
            }
            // JDIM boundary points
            {
                j = jdim-1; jl = jdim-2; ju = jdim-1;
                for (k = 1; k <= kdim-1; k++) {
                    dfacep=0.5*(blend(j,k,i)+blend(ju,k,i));
                    dfacem=0.5*(blend(j,k,i)+blend(jl,k,i));
                    sigkp=dfacep*sigk1+(1.-dfacep)*sigk2;
                    sigkm=dfacem*sigk1+(1.-dfacem)*sigk2;
                    sigop=dfacep*sigo1+(1.-dfacep)*sigo2;
                    sigom=dfacem*sigo1+(1.-dfacem)*sigo2;
                    volju=volj0(k,i,3);
                    xp=sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    yp=sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    zp=sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volju));
                    voljl=vol(jl,k,i);
                    xm=sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    ym=sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    zm=sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)+voljl));
                    xa=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    ya=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    za=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    ttpo=xp*xa+yp*ya+zp*za;
                    ttmo=xm*xa+ym*ya+zm*za;
                    ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/vol(j,k,i);
                    ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/vol(j,k,i);
                    ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                    ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                    anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i));
                    anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i));
                    fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i));
                    fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i));
                    cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                    bx(k,j)=-cdm; cx(k,j)=cdp+cdm; dx(k,j)=-cdp;
                    cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                    cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                    bx2(k,j)=-cdm; cx2(k,j)=cdp+cdm; dx2(k,j)=-cdp;
                }
                for (k = 1; k <= kdim-1; k++) {
                    xc=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i);
                    yc=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i);
                    zc=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i);
                    tc=0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+sj(j,k,i,5)*sj(j,k,i,4))/vol(j,k,i);
                    uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                    sgnu=ccsignrc(1.,uu); app=0.5*(1.+sgnu); apm=0.5*(1.-sgnu);
                    bx(k,j)=bx(k,j)-uu*app; bx2(k,j)=bx2(k,j)-uu*app;
                    cx(k,j)=cx(k,j)+uu*(app-apm); cx2(k,j)=cx2(k,j)+uu*(app-apm);
                    dx(k,j)=dx(k,j)+uu*apm; dx2(k,j)=dx2(k,j)+uu*apm;
                }
                // NOTE: the source-diagonal (cxadd) belongs ONLY to the first (eta)
                // ADI sweep in twoeqn.f; the xi sweep must NOT re-add it (doing so
                // makes the diagonal ~1e6x too large and kills the implicit update).
                for (k = 1; k <= kdim-1; k++) {
                    bx(k,j)=bx(k,j)*timestp(j,k,i);
                    bx2(k,j)=bx2(k,j)*timestp(j,k,i)*factor2;
                    cx(k,j)=cx(k,j)*timestp(j,k,i)+1.0*(1.+phi);
                    cx2(k,j)=cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                    dx(k,j)=dx(k,j)*timestp(j,k,i);
                    dx2(k,j)=dx2(k,j)*timestp(j,k,i)*factor2;
                    fx(k,j)=rhside(j,k,i,1)*(1.+phi);
                    fx2(k,j)=rhside(j,k,i,2)*(1.+phi);
                }
            }
            // Apply blank array for overset
            if (iover == 1) {
                for (j = 1; j <= jdim-1; j++)
                for (k = 1; k <= kdim-1; k++) {
                    fx(k,j)=fx(k,j)*blank(j,k,i);
                    bx(k,j)=bx(k,j)*blank(j,k,i);
                    dx(k,j)=dx(k,j)*blank(j,k,i);
                    cx(k,j)=cx(k,j)*blank(j,k,i)+(1.-blank(j,k,i));
                    fx2(k,j)=fx2(k,j)*blank(j,k,i);
                    bx2(k,j)=bx2(k,j)*blank(j,k,i);
                    dx2(k,j)=dx2(k,j)*blank(j,k,i);
                    cx2(k,j)=cx2(k,j)*blank(j,k,i)+(1.-blank(j,k,i));
                }
            }
            // Solve tridiagonal systems
            triv(kdim-1, jdim-1, 1, kdim-1, 1, jdim-1, workx, bx, cx, dx, fx);
            triv(kdim-1, jdim-1, 1, kdim-1, 1, jdim-1, workx, bx2, cx2, dx2, fx2);
            // Store results back
            for (j = 1; j <= jdim-1; j++)
            for (k = 1; k <= kdim-1; k++) {
                rhside(j,k,i,1) = fx(k,j);
                rhside(j,k,i,2) = fx2(k,j);
            }
        } // end i loop for xi sweep

        // Implicit F_zeta_zeta viscous terms. Do over all j's
        if (i2d != 1 && iaxi2planeturb != 1) {
            for (j = 1; j <= jdim-1; j++) {
                auto get_czadd = [&](int k, int i) -> std::pair<double,double> {
                    double czadd_v, cz2add_v;
                    if (ivmx == 6) {
                        czadd_v = 2.*re*beta1*turre(j,k,i,1);
                        cz2add_v = re*cmuc1*turre(j,k,i,1);
                    } else if (ivmx == 7) {
                        cmuc_var = blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2;
                        betax = blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2;
                        czadd_v = 2.*re*betax*turre(j,k,i,1)+ccabs(damp1(j,k,i))/turre(j,k,i,1);
                        cz2add_v = re*cmuc_var*turre(j,k,i,1);
                    } else if (ivmx == 12) {
                        czadd_v = 2.*re*beta1*turre(j,k,i,1);
                        cz2add_v = re*turre(j,k,i,1);
                    } else if (ivmx==8||ivmx==14) {
                        czadd_v = 2.*re*beta1*turre(j,k,i,1);
                        cz2add_v = damp1(j,k,i)*re*turre(j,k,i,1);
                    } else if (ivmx==9||ivmx==10||ivmx==11||ivmx==13) {
                        czadd_v = damp1(j,k,i);
                        cz2add_v = 2.*re*turre(j,k,i,1)/turre(j,k,i,2);
                    } else if (ivmx == 15) {
                        term1 = -0.5*(srce(j,k,i,1)-ccabs(srce(j,k,i,1)));
                        term2 = -0.5*(srce(j,k,i,2)-ccabs(srce(j,k,i,2)));
                        czadd_v = term1/turre(j,k,i,1);
                        cz2add_v = term2/turre(j,k,i,2);
                    } else if (ivmx == 16) {
                        czadd_v = 6.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*smin(j,k,i))*damp1(j,k,i)/re;
                        cz2add_v = 2.5*std::pow(cmuc1,0.75)*std::pow(turre(j,k,i,2),1.5)/turre(j,k,i,1)*re
                            +2.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*smin(j,k,i))/re;
                    } else {
                        czadd_v = 0.; cz2add_v = 0.;
                    }
                    return {czadd_v, cz2add_v};
                };
                // Interior points
                for (i = 2; i <= idim-2; i++) {
                    il = i-1; iu = i+1;
                    for (k = 1; k <= kdim-1; k++) {
                        dfacep=0.5*(blend(j,k,i)+blend(j,k,iu));
                        dfacem=0.5*(blend(j,k,i)+blend(j,k,il));
                        sigkp=dfacep*sigk1+(1.-dfacep)*sigk2;
                        sigkm=dfacem*sigk1+(1.-dfacem)*sigk2;
                        sigop=dfacep*sigo1+(1.-dfacep)*sigo2;
                        sigom=dfacem*sigo1+(1.-dfacem)*sigo2;
                        voliu=vol(j,k,iu);
                        xp=si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        yp=si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        zp=si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        volil=vol(j,k,il);
                        xm=si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        ym=si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        zm=si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        xa=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        ya=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        za=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        ttpo=xp*xa+yp*ya+zp*za;
                        ttmo=xm*xa+ym*ya+zm*za;
                        ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                        ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                        ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                        ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                        anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                        anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                        fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i));
                        fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i));
                        cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                        bz(k,i)=-cdm; cz(k,i)=cdp+cdm; dz(k,i)=-cdp;
                        cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                        bz2(k,i)=-cdm; cz2(k,i)=cdp+cdm; dz2(k,i)=-cdp;
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        xc=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        yc=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        zc=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        tc=0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                        uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                        sgnu=ccsignrc(1.,uu); app=0.5*(1.+sgnu); apm=0.5*(1.-sgnu);
                        bz(k,i)=bz(k,i)-uu*app; bz2(k,i)=bz2(k,i)-uu*app;
                        cz(k,i)=cz(k,i)+uu*(app-apm); cz2(k,i)=cz2(k,i)+uu*(app-apm);
                        dz(k,i)=dz(k,i)+uu*apm; dz2(k,i)=dz2(k,i)+uu*apm;
                    }
                    // NOTE: the source-term diagonal (czadd) is deliberately NOT
                    // added in the zeta (i) sweep. In the factored-ADI scheme the
                    // source is folded into the eta (j) sweep diagonal ONLY; adding
                    // it again here (as the original translation did) inflated the
                    // i-sweep diagonal ~500x, collapsing the turbulence update and
                    // diverging the solution. Fortran twoeqn's zeta sweep has no czadd.
                    for (k = 1; k <= kdim-1; k++) {
                        bz(k,i)=bz(k,i)*timestp(j,k,i);
                        bz2(k,i)=bz2(k,i)*timestp(j,k,i)*factor2;
                        cz(k,i)=cz(k,i)*timestp(j,k,i)+1.0*(1.+phi);
                        cz2(k,i)=cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                        dz(k,i)=dz(k,i)*timestp(j,k,i);
                        dz2(k,i)=dz2(k,i)*timestp(j,k,i)*factor2;
                        fz(k,i)=rhside(j,k,i,1)*(1.+phi);
                        fz2(k,i)=rhside(j,k,i,2)*(1.+phi);
                    }
                } // end interior i loop

                // I0 boundary points
                {
                    i = 1; il = 1; iu = std::min(2, idim-1);
                    for (k = 1; k <= kdim-1; k++) {
                        dfacep=0.5*(blend(j,k,i)+blend(j,k,iu));
                        dfacem=0.5*(blend(j,k,i)+blend(j,k,il));
                        sigkp=dfacep*sigk1+(1.-dfacep)*sigk2;
                        sigkm=dfacem*sigk1+(1.-dfacem)*sigk2;
                        sigop=dfacep*sigo1+(1.-dfacep)*sigo2;
                        sigom=dfacem*sigo1+(1.-dfacem)*sigo2;
                        voliu=vol(j,k,iu);
                        xp=si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        yp=si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        zp=si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        volil=voli0(j,k,1);
                        xm=si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        ym=si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        zm=si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        xa=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        ya=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        za=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        ttpo=xp*xa+yp*ya+zp*za;
                        ttmo=xm*xa+ym*ya+zm*za;
                        ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                        ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                        ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                        ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                        anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                        anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                        fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i));
                        fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i));
                        cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                        bz(k,i)=-cdm; cz(k,i)=cdp+cdm; dz(k,i)=-cdp;
                        cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                        bz2(k,i)=-cdm; cz2(k,i)=cdp+cdm; dz2(k,i)=-cdp;
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        xc=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        yc=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        zc=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        tc=0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                        uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                        sgnu=ccsignrc(1.,uu); app=0.5*(1.+sgnu); apm=0.5*(1.-sgnu);
                        bz(k,i)=bz(k,i)-uu*app; bz2(k,i)=bz2(k,i)-uu*app;
                        cz(k,i)=cz(k,i)+uu*(app-apm); cz2(k,i)=cz2(k,i)+uu*(app-apm);
                        dz(k,i)=dz(k,i)+uu*apm; dz2(k,i)=dz2(k,i)+uu*apm;
                    }
                    // NOTE: the source-term diagonal (czadd) is deliberately NOT
                    // added in the zeta (i) sweep. In the factored-ADI scheme the
                    // source is folded into the eta (j) sweep diagonal ONLY; adding
                    // it again here (as the original translation did) inflated the
                    // i-sweep diagonal ~500x, collapsing the turbulence update and
                    // diverging the solution. Fortran twoeqn's zeta sweep has no czadd.
                    for (k = 1; k <= kdim-1; k++) {
                        bz(k,i)=bz(k,i)*timestp(j,k,i);
                        bz2(k,i)=bz2(k,i)*timestp(j,k,i)*factor2;
                        cz(k,i)=cz(k,i)*timestp(j,k,i)+1.0*(1.+phi);
                        cz2(k,i)=cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                        dz(k,i)=dz(k,i)*timestp(j,k,i);
                        dz2(k,i)=dz2(k,i)*timestp(j,k,i)*factor2;
                        fz(k,i)=rhside(j,k,i,1)*(1.+phi);
                        fz2(k,i)=rhside(j,k,i,2)*(1.+phi);
                    }
                }
                // IDIM boundary points
                {
                    i = idim-1; il = idim-2; iu = idim-1;
                    for (k = 1; k <= kdim-1; k++) {
                        dfacep=0.5*(blend(j,k,i)+blend(j,k,iu));
                        dfacem=0.5*(blend(j,k,i)+blend(j,k,il));
                        sigkp=dfacep*sigk1+(1.-dfacep)*sigk2;
                        sigkm=dfacem*sigk1+(1.-dfacem)*sigk2;
                        sigop=dfacep*sigo1+(1.-dfacep)*sigo2;
                        sigom=dfacem*sigo1+(1.-dfacem)*sigo2;
                        voliu=voli0(j,k,3);
                        xp=si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        yp=si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        zp=si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)+voliu));
                        volil=vol(j,k,il);
                        xm=si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        ym=si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        zm=si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)+volil));
                        xa=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        ya=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        za=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        ttpo=xp*xa+yp*ya+zp*za;
                        ttmo=xm*xa+ym*ya+zm*za;
                        ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/vol(j,k,i);
                        ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/vol(j,k,i);
                        ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis;
                        ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis;
                        anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1));
                        anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1));
                        fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i));
                        fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i));
                        cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re);
                        bz(k,i)=-cdm; cz(k,i)=cdp+cdm; dz(k,i)=-cdp;
                        cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re);
                        cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re);
                        bz2(k,i)=-cdm; cz2(k,i)=cdp+cdm; dz2(k,i)=-cdp;
                    }
                    for (k = 1; k <= kdim-1; k++) {
                        xc=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i);
                        yc=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i);
                        zc=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i);
                        tc=0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+si(j,k,i,5)*si(j,k,i,4))/vol(j,k,i);
                        uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc;
                        sgnu=ccsignrc(1.,uu); app=0.5*(1.+sgnu); apm=0.5*(1.-sgnu);
                        bz(k,i)=bz(k,i)-uu*app; bz2(k,i)=bz2(k,i)-uu*app;
                        cz(k,i)=cz(k,i)+uu*(app-apm); cz2(k,i)=cz2(k,i)+uu*(app-apm);
                        dz(k,i)=dz(k,i)+uu*apm; dz2(k,i)=dz2(k,i)+uu*apm;
                    }
                    // NOTE: the source-term diagonal (czadd) is deliberately NOT
                    // added in the zeta (i) sweep. In the factored-ADI scheme the
                    // source is folded into the eta (j) sweep diagonal ONLY; adding
                    // it again here (as the original translation did) inflated the
                    // i-sweep diagonal ~500x, collapsing the turbulence update and
                    // diverging the solution. Fortran twoeqn's zeta sweep has no czadd.
                    for (k = 1; k <= kdim-1; k++) {
                        bz(k,i)=bz(k,i)*timestp(j,k,i);
                        bz2(k,i)=bz2(k,i)*timestp(j,k,i)*factor2;
                        cz(k,i)=cz(k,i)*timestp(j,k,i)+1.0*(1.+phi);
                        cz2(k,i)=cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.+phi);
                        dz(k,i)=dz(k,i)*timestp(j,k,i);
                        dz2(k,i)=dz2(k,i)*timestp(j,k,i)*factor2;
                        fz(k,i)=rhside(j,k,i,1)*(1.+phi);
                        fz2(k,i)=rhside(j,k,i,2)*(1.+phi);
                    }
                }
                // Apply blank array for overset
                if (iover == 1) {
                    for (i = 1; i <= idim-1; i++)
                    for (k = 1; k <= kdim-1; k++) {
                        fz(k,i)=fz(k,i)*blank(j,k,i);
                        bz(k,i)=bz(k,i)*blank(j,k,i);
                        dz(k,i)=dz(k,i)*blank(j,k,i);
                        cz(k,i)=cz(k,i)*blank(j,k,i)+(1.-blank(j,k,i));
                        fz2(k,i)=fz2(k,i)*blank(j,k,i);
                        bz2(k,i)=bz2(k,i)*blank(j,k,i);
                        dz2(k,i)=dz2(k,i)*blank(j,k,i);
                        cz2(k,i)=cz2(k,i)*blank(j,k,i)+(1.-blank(j,k,i));
                    }
                }
                // Solve tridiagonal systems
                triv(kdim-1, idim-1, 1, kdim-1, 1, idim-1, workz, bz, cz, dz, fz);
                triv(kdim-1, idim-1, 1, kdim-1, 1, idim-1, workz, bz2, cz2, dz2, fz2);
                // Store results back
                for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++) {
                    rhside(j,k,i,1) = fz(k,i);
                    rhside(j,k,i,2) = fz2(k,i);
                }
            } // end j loop for zeta sweep
        } // end i2d != 1 && iaxi2planeturb != 1

        // Update TURRE
        sumno = 0.;
        sumnk = 0.;
        negno = 0;
        negnk = 0;
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            sumno = sumno + rhside(j,k,i,1)*rhside(j,k,i,1);
            if ((float)(turre(j,k,i,1)+rhside(j,k,i,1)) <= (float)tur1cutlev) {
                negno = negno + 1;
                if ((float)tur1cut > 0.) turre(j,k,i,1) = tur1cut;
            } else {
                turre(j,k,i,1) = turre(j,k,i,1) + rhside(j,k,i,1);
            }
            sumnk = sumnk + rhside(j,k,i,2)*rhside(j,k,i,2);
            if ((float)(turre(j,k,i,2)+rhside(j,k,i,2)) <= (float)tur2cutlev) {
                negnk = negnk + 1;
                if ((float)tur2cut > 0.) turre(j,k,i,2) = tur2cut;
            } else {
                turre(j,k,i,2) = turre(j,k,i,2) + rhside(j,k,i,2);
            }
        }
        sumno = std::sqrt(sumno) / (float)((kdim-1)*(jdim-1)*(idim-1));
        sumnk = std::sqrt(sumnk) / (float)((kdim-1)*(jdim-1)*(idim-1));

    } // end do 500 not=1,nsubit

    sumn1 = sumno;
    sumn2 = sumnk;
    negn1 = negno;
    negn2 = negnk;

    // Update VIST3D and save omega and k values
    if (ivmx == 6) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            if (i_wilcox06 == 1) {
                s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
                s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                tracepart=(s11+s22+s33)*prod2d3dtrace;
                s11t=s11-tracepart; s22t=s22-tracepart; s33t=s33-tracepart;
                xis=s11t*s11t+s22t*s22t+s33t*s33t+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                xxx=std::sqrt(2.*xis);
                denom=ccmax(turre(j,k,i,1),xxx*0.875/(re*std::sqrt(cmuc2)));
                vist3d(j,k,i)=q(j,k,i,1)*turre(j,k,i,2)/denom;
            } else {
                vist3d(j,k,i)=q(j,k,i,1)*turre(j,k,i,2)/turre(j,k,i,1);
            }
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),(double)edvislim);
            zksav(j,k,i,1)=turre(j,k,i,1); zksav(j,k,i,2)=turre(j,k,i,2);
        }
    } else if (ivmx == 7) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            denom1=2./re*std::sqrt(turre(j,k,i,2))/(.09*turre(j,k,i,1)*ccabs(smin(j,k,i)));
            denom2=500.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*re*re*smin(j,k,i)*turre(j,k,i,1));
            arg2=ccmax(denom1,denom2);
            f2=cctanh(arg2*arg2);
            if (isstdenom==1 || isst2003==1) {
                s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
                s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
                xis=s11*s11+s22*s22+s33*s33+2.*s12*s12+2.*s13*s13+2.*s23*s23;
                xxx=std::sqrt(2.*xis);
                denom=ccmax(a1*turre(j,k,i,1),xxx*f2/re);
            } else {
                denom=ccmax(a1*turre(j,k,i,1),vor(j,k,i)*f2/re);
            }
            vist3d(j,k,i)=a1*q(j,k,i,1)*turre(j,k,i,2)/denom;
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),(double)edvislim);
            zksav(j,k,i,1)=turre(j,k,i,1); zksav(j,k,i,2)=turre(j,k,i,2);
        }
    } else if (ivmx == 12) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            alfa1=(1.333333-c2)*gg/2.;
            alfa2=(2.-c3)*(2.-c3)*(gg*gg/4.);
            alfa3=(2.-c4)*(2.-c4)*(gg*gg/4.);
            s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
            s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
            w12=0.5*(ux(j,k,i,2)-ux(j,k,i,4)); w13=0.5*(ux(j,k,i,3)-ux(j,k,i,7));
            if (ieasmcc2d==1) w13=w13+(2./(c4-2.))*vx(j,k,i,1);
            w23=0.5*(ux(j,k,i,6)-ux(j,k,i,8));
            xis=s11*s11+s22*s22+s33*s33+2.*s12*s12+2.*s13*s13+2.*s23*s23;
            wis=2.*w12*w12+2.*w13*w13+2.*w23*w23;
            eta=alfa2*xis/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
            squig=alfa3*wis/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
            eta=ccmincr(eta,10.); squig=ccmincr(squig,10.);
            cmu=alfa1*((3.*(1.+eta)+0.2*eta*eta*eta+0.2*squig*squig*squig)/
                (3.+eta+6.*eta*squig+6.*squig+eta*eta*eta+squig*squig*squig));
            cmu=ccmincr(cmu,.187); cmu=ccmaxcr(cmu,.005);
            vist3d(j,k,i)=cmu*q(j,k,i,1)*turre(j,k,i,2)/turre(j,k,i,1);
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),(double)edvislim);
            zksav(j,k,i,1)=turre(j,k,i,1); zksav(j,k,i,2)=turre(j,k,i,2);
        }
    } else if (ivmx == 10) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            ret=q(j,k,i,1)*turre(j,k,i,2)*turre(j,k,i,2)/(fnu(j,k,i)*turre(j,k,i,1));
            rek=q(j,k,i,1)*std::sqrt(turre(j,k,i,2))*re*ccabs(smin(j,k,i))/fnu(j,k,i);
            fmu=(1.+4.*std::pow(ret,-0.75))*cctanh(.008*rek);
            fmu=ccmincr(fmu,1.0);
            vist3d(j,k,i)=cmuc1*fmu*q(j,k,i,1)*turre(j,k,i,2)*turre(j,k,i,2)/turre(j,k,i,1);
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),(double)edvislim);
            zksav(j,k,i,1)=turre(j,k,i,1); zksav(j,k,i,2)=turre(j,k,i,2);
        }
    } else if (ivmx == 11) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            alfa1=(1.333333-c2)*gg/2.;
            alfa2=(2.-c3)*(2.-c3)*(gg*gg/4.);
            alfa3=(2.-c4)*(2.-c4)*(gg*gg/4.);
            s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
            s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
            w12=0.5*(ux(j,k,i,2)-ux(j,k,i,4)); w13=0.5*(ux(j,k,i,3)-ux(j,k,i,7));
            if (ieasmcc2d==1) w13=w13+(2./(c4-2.))*vx(j,k,i,1);
            w23=0.5*(ux(j,k,i,6)-ux(j,k,i,8));
            xis=s11*s11+s22*s22+s33*s33+2.*s12*s12+2.*s13*s13+2.*s23*s23;
            wis=2.*w12*w12+2.*w13*w13+2.*w23*w23;
            eta=alfa2*xis*turre(j,k,i,2)*turre(j,k,i,2)/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
            squig=alfa3*wis*turre(j,k,i,2)*turre(j,k,i,2)/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
            eta=ccmincr(eta,10.); squig=ccmincr(squig,10.);
            cmu=alfa1*((3.*(1.+eta)+0.2*eta*eta*eta+0.2*squig*squig*squig)/
                (3.+eta+6.*eta*squig+6.*squig+eta*eta*eta+squig*squig*squig));
            cmu=ccmincr(cmu,.187); cmu=ccmaxcr(cmu,.005);
            vist3d(j,k,i)=cmu*q(j,k,i,1)*turre(j,k,i,2)*turre(j,k,i,2)/turre(j,k,i,1);
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),(double)edvislim);
            zksav(j,k,i,1)=turre(j,k,i,1); zksav(j,k,i,2)=turre(j,k,i,2);
        }
    }

    else if (ivmx==9 || ivmx==13) {
        pi=std::acos(-1.); tpi3=2.*pi/3.;
        al10=0.5*c10-1.; al1=2.*(0.5*c11+1.);
        if (ieasm_type==0) { al10=al10+1.8864; al1=al1-2.; }
        al2=0.5*c2-2./3.; al3=0.5*c3-1.; al4=0.5*c4-1.;
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            fmu=1.0;
            s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
            s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
            w12=0.5*(ux(j,k,i,2)-ux(j,k,i,4)); w13=0.5*(ux(j,k,i,3)-ux(j,k,i,7));
            if (ieasmcc2d==1) w13=w13+(2./(c4-2.))*vx(j,k,i,1);
            w23=0.5*(ux(j,k,i,6)-ux(j,k,i,8));
            xis=s11*s11+s22*s22+s33*s33+2.*s12*s12+2.*s13*s13+2.*s23*s23;
            wis=2.*w12*w12+2.*w13*w13+2.*w23*w23;
            eta1=xis*turre(j,k,i,2)*turre(j,k,i,2)/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
            eta2=wis*turre(j,k,i,2)*turre(j,k,i,2)/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
            eta1=ccmincr(eta1,1200.); eta2=ccmincr(eta2,1200.);
            if ((float)eta1 <= 1.e-6) {
                cmuj=al10*al2/(al10*al10+2.*eta2*al4*al4);
            } else {
                eta11=eta1*al1;
                ap=-2.*al10/eta11;
                ar=-al10*al2/(eta11*eta11);
                aq=al10*al10+eta11*al2-0.666667*eta1*al3*al3+2.*eta2*al4*al4;
                aq=aq/(eta11*eta11);
                aa=(aq-ap*ap/3.);
                ab=(2.*ap*ap*ap-9.*ap*aq+27.*ar)/27.;
                ad=(ab*ab/4.)+(aa*aa*aa)/27.;
                if ((float)ad > 0.) {
                    raat=-0.5*ab+std::sqrt(ad);
                    rbbt=-0.5*ab-std::sqrt(ad);
                    raa=std::pow(ccabs(raat),1./3.); raa=ccsign(raa,raat);
                    rbb=std::pow(ccabs(rbbt),1./3.); rbb=ccsign(rbb,rbbt);
                    cmuj=-ap/3.+raa+rbb;
                    cmub=-ap/3.-.5*raa-.5*rbb;
                    cmuj=ccmin(cmuj,cmub);
                } else {
                    coss=-ab/2./std::sqrt(-aa*aa*aa/27.);
                    theta=ccacos(coss);
                    cmuj=-ap/3.+2.*std::sqrt(-aa/3.)*std::cos(theta/3.);
                    cmub=-ap/3.+2.*std::sqrt(-aa/3.)*std::cos(tpi3+theta/3.);
                    cmuc_var=-ap/3.+2.*std::sqrt(-aa/3.)*std::cos(2.*tpi3+theta/3.);
                    cmuj=ccmin(cmuj,cmub); cmuj=ccmin(cmuj,cmuc_var);
                }
            }
            cmuj=ccmincr(cmuj,-(double)cmulim);
            cmu=-fmu*cmuj;
            vist3d(j,k,i)=cmu*q(j,k,i,1)*turre(j,k,i,2)*turre(j,k,i,2)/turre(j,k,i,1);
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),(double)edvislim);
            zksav(j,k,i,1)=turre(j,k,i,1); zksav(j,k,i,2)=turre(j,k,i,2);
            cmuv(j,k,i)=cmuj;
        }
    } else if (ivmx==8 || ivmx==14) {
        pi=std::acos(-1.); tpi3=2.*pi/3.;
        al10=0.5*c10-1.; al1=2.*(0.5*c11+1.);
        if (ieasm_type==0||ieasm_type==3||ieasm_type==4) { al10=al10+1.8864; al1=al1-2.; }
        al2=0.5*c2-2./3.; al3=0.5*c3-1.; al4=0.5*c4-1.;
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            s11=ux(j,k,i,1); s22=ux(j,k,i,5); s33=ux(j,k,i,9);
            s12=0.5*(ux(j,k,i,2)+ux(j,k,i,4)); s13=0.5*(ux(j,k,i,3)+ux(j,k,i,7)); s23=0.5*(ux(j,k,i,6)+ux(j,k,i,8));
            w12=0.5*(ux(j,k,i,2)-ux(j,k,i,4)); w13=0.5*(ux(j,k,i,3)-ux(j,k,i,7));
            if (ieasmcc2d==1) w13=w13+(2./(c4-2.))*vx(j,k,i,1);
            w23=0.5*(ux(j,k,i,6)-ux(j,k,i,8));
            xis=s11*s11+s22*s22+s33*s33+2.*s12*s12+2.*s13*s13+2.*s23*s23;
            wis=2.*w12*w12+2.*w13*w13+2.*w23*w23;
            eta1=xis/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
            eta2=wis/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
            if (idurbinlim!=0 && (ieasm_type==3||ieasm_type==4)) {
                tau=1./turre(j,k,i,1);
                taulim=6.*std::sqrt(fnu(j,k,i)/(q(j,k,i,1)*turre(j,k,i,1)*turre(j,k,i,2)));
                tau=ccmax(tau,taulim);
                eta1=xis*(tau/re)*(tau/re);
                eta2=wis*(tau/re)*(tau/re);
            }
            fff=0.5*(1.+cctanh(4.0/1.414214*std::sqrt(eta1)-10.0));
            al2=-0.5*(4./3.-c2+(c2-0.36)*fff);
            if (ieasm_type==4) {
                eta1_girimaji=xis/(xis+wis);
                if ((float)eta1_girimaji < 0.5) {
                    c4new=2.0-((2.0-c4)*std::pow(eta1_girimaji/(1.-eta1_girimaji),0.75));
                } else { c4new=c4; }
                al4=0.5*c4new-1.;
            }
            if ((float)eta1 <= 1.e-6) {
                cmuj=al10*al2/(al10*al10+2.*eta2*al4*al4);
            } else {
                eta11=eta1*al1;
                ap=-2.*al10/eta11;
                ar=-al10*al2/(eta11*eta11);
                aq=al10*al10+eta11*al2-0.666667*eta1*al3*al3+2.*eta2*al4*al4;
                if (ieasm_type==3||ieasm_type==4) {
                    ap=ap+css*(1.-fff)/al1;
                    aq=aq-al10*css*(1.-fff)*eta1;
                }
                aq=aq/(eta11*eta11);
                aa=(aq-ap*ap/3.);
                ab=(2.*ap*ap*ap-9.*ap*aq+27.*ar)/27.;
                ad=(ab*ab/4.)+(aa*aa*aa)/27.;
                if ((float)ad > 0.) {
                    raat=-0.5*ab+std::sqrt(ad);
                    rbbt=-0.5*ab-std::sqrt(ad);
                    raa=std::pow(ccabs(raat),1./3.); raa=ccsign(raa,raat);
                    rbb=std::pow(ccabs(rbbt),1./3.); rbb=ccsign(rbb,rbbt);
                    cmuj=-ap/3.+raa+rbb;
                    cmub=-ap/3.-.5*raa-.5*rbb;
                    cmuj=ccmin(cmuj,cmub);
                } else {
                    coss=-ab/2./std::sqrt(-aa*aa*aa/27.);
                    theta=ccacos(coss);
                    cmuj=-ap/3.+2.*std::sqrt(-aa/3.)*std::cos(theta/3.);
                    cmub=-ap/3.+2.*std::sqrt(-aa/3.)*std::cos(tpi3+theta/3.);
                    cmuc_var=-ap/3.+2.*std::sqrt(-aa/3.)*std::cos(2.*tpi3+theta/3.);
                    cmuj=ccmin(cmuj,cmub); cmuj=ccmin(cmuj,cmuc_var);
                }
            }
            cmuj=ccmincr(cmuj,-(double)cmulim);
            if (ieasm_type==3||ieasm_type==4) {
                facy=-1./(std::sqrt(6.*eta1));
                cmuj=ccmaxcr(cmuj,facy);
            }
            cmu=-cmuj;
            if (idurbinlim!=0 && (ieasm_type==3||ieasm_type==4)) {
                vist3d(j,k,i)=cmu*q(j,k,i,1)*turre(j,k,i,2)*tau;
            } else {
                vist3d(j,k,i)=cmu*q(j,k,i,1)*turre(j,k,i,2)/turre(j,k,i,1);
            }
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),(double)edvislim);
            zksav(j,k,i,1)=turre(j,k,i,1); zksav(j,k,i,2)=turre(j,k,i,2);
            cmuv(j,k,i)=cmuj;
        }
    } else if (ivmx == 15) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            vist3d(j,k,i)=q(j,k,i,1)*q(j,k,i,1)*cmu_ke*turre(j,k,i,2)*turre(j,k,i,2)/(fnu(j,k,i)*turre(j,k,i,1));
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),(double)edvislim);
            zksav(j,k,i,1)=turre(j,k,i,1); zksav(j,k,i,2)=turre(j,k,i,2);
        }
    } else if (ivmx == 16) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            vist3d(j,k,i)=std::pow(cmuc1,0.25)*q(j,k,i,1)*turre(j,k,i,1)/std::sqrt(turre(j,k,i,2));
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),(double)edvislim);
            zksav(j,k,i,1)=turre(j,k,i,1); zksav(j,k,i,2)=turre(j,k,i,2);
        }
    }

    // Wall functions - k-direction
    if (iwf[2] == 1) {
        for (int kk = 1; kk <= 2; kk++) {
            if (kk == 1) { k = 1; kstop = nbck0(nbl); }
            else { k = kdim-1; kstop = nbckdim(nbl); }
            for (kset = 1; kset <= kstop; kset++) {
                int bc = std::abs(kbcinfo(nbl,kset,1,kk));
                if (bc==2004||bc==2024||bc==2034||bc==2016) {
                    ibeg=kbcinfo(nbl,kset,2,kk); iend=kbcinfo(nbl,kset,3,kk)-1;
                    jbeg=kbcinfo(nbl,kset,4,kk); jend=kbcinfo(nbl,kset,5,kk)-1;
                    for (i = ibeg; i <= iend; i++)
                    for (j = jbeg; j <= jend; j++) {
                        if ((i>=ilamlo&&i<ilamhi&&j>=jlamlo&&j<jlamhi&&k>=klamlo&&k<klamhi)||(float)smin(j,k,i)<0.) continue;
                        uu=std::sqrt((q(j,k,i,2)-qk0(j,i,2,kk+kk-1))*(q(j,k,i,2)-qk0(j,i,2,kk+kk-1))+
                            (q(j,k,i,3)-qk0(j,i,3,kk+kk-1))*(q(j,k,i,3)-qk0(j,i,3,kk+kk-1))+
                            (q(j,k,i,4)-qk0(j,i,4,kk+kk-1))*(q(j,k,i,4)-qk0(j,i,4,kk+kk-1)));
                        dudy=uu/ccabs(smin(j,k,i));
                        tauw=(fnu(j,k,i)+vk0(j,i,1,kk+kk-1))*dudy;
                        utau=std::sqrt(tauw/q(j,k,i,1)/re);
                        omegatemp=utau/(std::sqrt(cmuc1)*vk*ccabs(smin(j,k,i))*re);
                        omegatemp=ccmax(omegatemp,(double)tur10[0]);
                        zktemp=utau*utau/std::sqrt(cmuc1);
                        vist3d(j,k,i)=q(j,k,i,1)*zktemp/omegatemp;
                        if (ivmx==6||ivmx==7) { zksav(j,k,i,1)=omegatemp; zksav(j,k,i,2)=zktemp; }
                        else if (ivmx==8||ivmx==12||ivmx==14) { zksav(j,k,i,1)=omegatemp*cmuc1; zksav(j,k,i,2)=zktemp; }
                        else if (ivmx==15) { zksav(j,k,i,1)=omegatemp*cmuc1*zktemp*q(j,k,i,1)/fnu(j,k,i); zksav(j,k,i,2)=zktemp; }
                        else { zksav(j,k,i,1)=omegatemp*cmuc1*zktemp; zksav(j,k,i,2)=zktemp; }
                    }
                }
            }
        }
    }
    // Wall functions - j-direction
    if (iwf[1] == 1) {
        for (int jj = 1; jj <= 2; jj++) {
            if (jj == 1) { j = 1; jstop = nbcj0(nbl); }
            else { j = jdim-1; jstop = nbcjdim(nbl); }
            for (jset_bc = 1; jset_bc <= jstop; jset_bc++) {
                int bc = std::abs(jbcinfo(nbl,jset_bc,1,jj));
                if (bc==2004||bc==2024||bc==2034||bc==2016) {
                    ibeg=jbcinfo(nbl,jset_bc,2,jj); iend=jbcinfo(nbl,jset_bc,3,jj)-1;
                    kbeg=jbcinfo(nbl,jset_bc,4,jj); kend=jbcinfo(nbl,jset_bc,5,jj)-1;
                    for (i = ibeg; i <= iend; i++)
                    for (k = kbeg; k <= kend; k++) {
                        if ((i>=ilamlo&&i<ilamhi&&j>=jlamlo&&j<jlamhi&&k>=klamlo&&k<klamhi)||(float)smin(j,k,i)<0.) continue;
                        uu=std::sqrt((q(j,k,i,2)-qj0(k,i,2,jj+jj-1))*(q(j,k,i,2)-qj0(k,i,2,jj+jj-1))+
                            (q(j,k,i,3)-qj0(k,i,3,jj+jj-1))*(q(j,k,i,3)-qj0(k,i,3,jj+jj-1))+
                            (q(j,k,i,4)-qj0(k,i,4,jj+jj-1))*(q(j,k,i,4)-qj0(k,i,4,jj+jj-1)));
                        dudy=uu/ccabs(smin(j,k,i));
                        tauw=(fnu(j,k,i)+vj0(k,i,1,jj+jj-1))*dudy;
                        utau=std::sqrt(tauw/q(j,k,i,1)/re);
                        omegatemp=utau/(std::sqrt(cmuc1)*vk*ccabs(smin(j,k,i))*re);
                        omegatemp=ccmax(omegatemp,(double)tur10[0]);
                        zktemp=utau*utau/std::sqrt(cmuc1);
                        vist3d(j,k,i)=q(j,k,i,1)*zktemp/omegatemp;
                        if (ivmx==6||ivmx==7) { zksav(j,k,i,1)=omegatemp; zksav(j,k,i,2)=zktemp; }
                        else if (ivmx==8||ivmx==12||ivmx==14) { zksav(j,k,i,1)=omegatemp*cmuc1; zksav(j,k,i,2)=zktemp; }
                        else if (ivmx==15) { zksav(j,k,i,1)=omegatemp*cmuc1*zktemp*q(j,k,i,1)/fnu(j,k,i); zksav(j,k,i,2)=zktemp; }
                        else { zksav(j,k,i,1)=omegatemp*cmuc1*zktemp; zksav(j,k,i,2)=zktemp; }
                    }
                }
            }
        }
    }
    // Wall functions - i-direction
    if (i2d != 1 && iwf[0] == 1 && iaxi2planeturb != 1) {
        for (int ii = 1; ii <= 2; ii++) {
            if (ii == 1) { i = 1; istop = nbci0(nbl); }
            else { i = idim-1; istop = nbcidim(nbl); }
            for (iset = 1; iset <= istop; iset++) {
                int bc = std::abs(ibcinfo(nbl,iset,1,ii));
                if (bc==2004||bc==2024||bc==2034||bc==2016) {
                    jbeg=ibcinfo(nbl,iset,2,ii); jend=ibcinfo(nbl,iset,3,ii)-1;
                    kbeg=ibcinfo(nbl,iset,4,ii); kend=ibcinfo(nbl,iset,5,ii)-1;
                    for (j = jbeg; j <= jend; j++)
                    for (k = kbeg; k <= kend; k++) {
                        if ((i>=ilamlo&&i<ilamhi&&j>=jlamlo&&j<jlamhi&&k>=klamlo&&k<klamhi)||(float)smin(j,k,i)<0.) continue;
                        uu=std::sqrt((q(j,k,i,2)-qi0(k,i,2,ii+ii-1))*(q(j,k,i,2)-qi0(k,i,2,ii+ii-1))+
                            (q(j,k,i,3)-qi0(k,i,3,ii+ii-1))*(q(j,k,i,3)-qi0(k,i,3,ii+ii-1))+
                            (q(j,k,i,4)-qi0(k,i,4,ii+ii-1))*(q(j,k,i,4)-qi0(k,i,4,ii+ii-1)));
                        dudy=uu/ccabs(smin(j,k,i));
                        tauw=(fnu(j,k,i)+vi0(k,i,1,ii+ii-1))*dudy;
                        utau=std::sqrt(tauw/q(j,k,i,1)/re);
                        omegatemp=utau/(std::sqrt(cmuc1)*vk*ccabs(smin(j,k,i))*re);
                        omegatemp=ccmax(omegatemp,(double)tur10[0]);
                        zktemp=utau*utau/std::sqrt(cmuc1);
                        vist3d(j,k,i)=q(j,k,i,1)*zktemp/omegatemp;
                        if (ivmx==6||ivmx==7) { zksav(j,k,i,1)=omegatemp; zksav(j,k,i,2)=zktemp; }
                        else if (ivmx==8||ivmx==12||ivmx==14) { zksav(j,k,i,1)=omegatemp*cmuc1; zksav(j,k,i,2)=zktemp; }
                        else if (ivmx==15) { zksav(j,k,i,1)=omegatemp*cmuc1*zktemp*q(j,k,i,1)/fnu(j,k,i); zksav(j,k,i,2)=zktemp; }
                        else { zksav(j,k,i,1)=omegatemp*cmuc1*zktemp; zksav(j,k,i,2)=zktemp; }
                    }
                }
            }
        }
    }
    // Force vist3d=0 in laminar region
    if (i_lam_forcezero == 1) {
        for (i = 1; i <= idim-1; i++)
        for (k = 1; k <= kdim-1; k++)
        for (j = 1; j <= jdim-1; j++) {
            if ((i>=ilamlo&&i<ilamhi&&j>=jlamlo&&j<jlamhi&&k>=klamlo&&k<klamhi)||(float)smin(j,k,i)<0.) {
                vist3d(j,k,i)=0.;
            }
        }
    }
    // Hardwired output section (iwriteaux=0 always, so this is a no-op)
    iwriteaux = 0;
    // (iwriteaux is always 0 in the Fortran code, so no output is written)
} // end twoeqn

} // namespace twoeqn_ns
