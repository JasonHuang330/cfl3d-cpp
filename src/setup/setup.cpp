// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "setup.h"
#include "intrbc.h"
#include "metric.h"
#include "cellvol.h"
#include "dird.h"
#include "findmin_new.h"
#include "runtime/fortran_io.h"
#include "module_kwstm.h"
#include "cgnstools.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <algorithm>

// External subroutine declarations (called from setup())
// These are defined in other translation units
namespace lead_ns { void lead(int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, int& maxbl); }
namespace rcfl_ns { void rcfl(FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& jdim, int& kdim, int& idim, int& igrid, int& ialph, int& igeom_img, int& irr); }
namespace rp3d_ns { void rp3d(FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& jdim, int& kdim, int& idim, int& igrid, int& ialph, int& igeom_img, int& irr); }
namespace grdmove_ns { void grdmove(int& nbl, int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, double& xorig, double& yorig, double& zorig, double& xorg, double& yorg, double& zorg, double& thetax, double& thetay, double& thetaz); }
namespace collx_ns { void collx(FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> xx, FortranArray3DRef<double> yy, FortranArray3DRef<double> zz, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2); }
namespace writ_buf_ns { void writ_buf(int& nbl, int& iunit, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myhost, int& myid, int& mycomm, FortranArray1DRef<int> mblk2nd, int& maxbl); }
namespace collv_ns { void collv(FortranArray3DRef<double> vol, FortranArray3DRef<double> volc, int& nj, int& nk, int& ni, int& jj2, int& kk2, int& ii2); }
namespace readdat_ns { void readdat(FortranArray4DRef<double> bcdata, int& mdim, int& ndim, char (&filname)[80], FortranArray1DRef<double> bcdat, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl); }
namespace colldat_ns { void colldat(FortranArray4DRef<double> bcdata, int& mdim, int& ndim, FortranArray4DRef<double> bcdatac, int& mdimc, int& ndimc); }
namespace getibk_ns { void getibk(FortranArray3DRef<double> blk, int& jdim, int& kdim, int& idim, int& nbl, FortranArray1DRef<int> intpts, int& nblpts, int& ibpnts, int& iipnts, int& iitot, int& maxbl, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray1DRef<double> dxintg, FortranArray1DRef<double> dyintg, FortranArray1DRef<double> dzintg, FortranArray1DRef<int> iiig, FortranArray1DRef<int> jjig, FortranArray1DRef<int> kkig, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, int& myid, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou); }
namespace setblk_ns { void setblk(FortranArray3DRef<double> blk, int& jdim, int& kdim, int& idim, int& nbl); }
namespace rpatch_ns { void rpatch(int& maxbl, int& maxxe, int& intmax, int& nsub1, FortranArray2DRef<double> windex, int& ninter, FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt, FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid); }
namespace init_ns { void init(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, FortranArray3DRef<double> vol, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, int& nummem, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iflagprnt); }
namespace initnonin_ns { void initnonin(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray3DRef<double> vol, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z); }
namespace bc_ns { void bc(int& idum, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, double& clw, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, FortranArray3DRef<int> lwdat, int& myid, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, int& nummem); }
namespace init_rb_ns { void init_rb(int& maxaes, int& maxbl, FortranArray1DRef<double> zorig, FortranArray2DRef<double> aesrfdat); }
namespace init_trim_ns { void init_trim(); }
namespace init_ae_ns { void init_ae(int& nmds, int& maxaes, FortranArray2DRef<double> x0, FortranArray2DRef<double> xxn, FortranArray2DRef<double> temp, FortranArray3DRef<double> bmat, FortranArray3DRef<double> stm, FortranArray3DRef<double> stmi, FortranArray2DRef<double> gforcn, FortranArray2DRef<double> gforcnm, FortranArray2DRef<double> freq, FortranArray2DRef<double> damp, FortranArray2DRef<double> gmass, FortranArray2DRef<double> aesrfdat, FortranArray2DRef<double> xs, FortranArray2DRef<double> gforcs); }
namespace pltmode_ns { void pltmode(int& nm, int& iaes, FortranArray1DRef<int> nblg, int& ngrid, int& maxgr, int& maxbl, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> iaesurf, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, int& maxseg, int& maxsegdg, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray5DRef<double> xmdj, FortranArray5DRef<double> xmdk, FortranArray5DRef<double> xmdi, int& jdim1, int& kdim1, int& idim1, int& maxaes, int& nmds, int& nbl1, int& iunitw, int& iflag); }
namespace modread_ns { void modread(int& idim, int& jdim, int& kdim, int& n, int& nbl, int& iunit, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, int& maxbl, int& maxseg, int& nmds, FortranArray5DRef<double> xmdj, FortranArray5DRef<double> xmdk, FortranArray5DRef<double> xmdi, int& iaes, int& maxaes); }
namespace collmod_ns { void collmod(FortranArray5DRef<double> xmdj, FortranArray5DRef<double> xmdk, FortranArray5DRef<double> xmdi, FortranArray5DRef<double> xmdjc, FortranArray5DRef<double> xmdkc, FortranArray5DRef<double> xmdic, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, int& nm, int& nmds, int& iaes, int& maxaes); }
namespace rrest_ns { void rrest(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, int& ncycmax, int& ntr, FortranArray1DRef<double> rms, FortranArray1DRef<double> clw, FortranArray1DRef<double> cdw, FortranArray1DRef<double> cdpw, FortranArray1DRef<double> cdvw, FortranArray1DRef<double> cxw, FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw, FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw, FortranArray1DRef<double> cmzw, int& n_clcd, FortranArray3DRef<double> clcd, int& nblocks_clcd, FortranArray2DRef<int> blocks_clcd, FortranArray1DRef<double> fmdotw, FortranArray1DRef<double> cftmomw, FortranArray1DRef<double> cftpw, FortranArray1DRef<double> cftvw, FortranArray1DRef<double> cfttotw, FortranArray2DRef<double> rmstr, FortranArray2DRef<int> nneg, int& iskip, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> tursav, FortranArray3DRef<double> smin, FortranArray3DRef<double> xjb, FortranArray3DRef<double> xkb, FortranArray3DRef<double> blnum, FortranArray3DRef<double> cmuv, int& maxbl, FortranArray1DRef<int> mblk2nd, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& igrid, FortranArray3DRef<double> wk, int& idima, int& jdima, int& kdima, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0, FortranArray4DRef<double> qavg, FortranArray4DRef<double> q2avg, int& nummem); }
namespace rrestg_ns { void rrestg(int& nbl, int& igrid, int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> xnm2, FortranArray3DRef<double> ynm2, FortranArray3DRef<double> znm2, FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, FortranArray4DRef<double> delti, FortranArray4DRef<double> qc0, int& igrdflag, int& iuns, FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx, FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz, FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr, FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0, FortranArray1DRef<double> time2, FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae, FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae, FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass, FortranArray2DRef<double> damp, FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0, int& nmds, int& maxaes, FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb, int& myhost, int& myid, int& mycomm, FortranArray1DRef<int> mblk2nd, int& maxbl, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg, FortranArray1DRef<double> wk, int& nwork, int& idima, int& jdima, int& kdima, FortranArray4DRef<double> tursav2, int& nummem); }
namespace setslave_ns { void setslave(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& maxbl, int& mxbli, int& maxgr, int& maxseg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& nblock, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray3DRef<int> islavept, int& nslave, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, int& maxsegdg, FortranArray1DRef<int> iwk, int& iwork, int& nmaster, int& ngrid, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray2DRef<int> iskip, FortranArray1DRef<int> nblg, FortranArray1DRef<int> levelg, int& lfgm, FortranArray2DRef<int> nblk, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray2DRef<int> nblelst, int& nnodes, FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax, FortranArray1DRef<int> kskmax, int& nbli); }
namespace rsurf_ns { void rsurf(int& maxbl, int& maxsegdg, int& idim, int& jdim, int& kdim, FortranArray4DRef<double> delti, FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, int& nbl, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, int& iseg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iunit); }
namespace getdelt_ns { void getdelt(int& maxbl, int& maxsegdg, int& idim, int& jdim, int& kdim, FortranArray4DRef<double> delti, FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& nbl, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, int& iseg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray3DRef<double> wkj, FortranArray3DRef<double> wkk, FortranArray3DRef<double> wki); }
namespace cntsurf_ns { void cntsurf(int& ns, int& maxbl, int& maxgr, int& maxseg, int& ngrid, FortranArray1DRef<int> nblg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& bctype); }
namespace initvist_ns { void initvist(int& nbl, int& jdim, int& kdim, int& idim, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> tursav, FortranArray3DRef<double> smin, FortranArray3DRef<double> cmuv, int& nummem, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z); }
namespace lamfix_ns { void lamfix(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> snk0, int& mdim, int& ndim, FortranArray4DRef<double> ldata, int& ibeg, int& iend, int& jbeg, int& jend, int& kbeg, int& kend, int& iface, int& llev); }
namespace termn8_ns { void termn8(int& myid, int& ierr, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou); }
namespace ccabs_ns { float ccabs(float x); }

namespace setup_ns {

void setup(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, FortranArray1DRef<int> iwk, int& iwork, FortranArray1DRef<int> itest, FortranArray1DRef<int> jtest, FortranArray1DRef<int> ktest, int& maxbl, int& mxbli, int& maxgr, int& maxseg, int& nsub1, int& maxxe, int& intmax, int& iitot, int& ncycmax, FortranArray3DRef<int> lwdat, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray1DRef<int> iovrlp, FortranArray3DRef<double> qb, int& nblock, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx, FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz, FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr, FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0, FortranArray1DRef<double> time2, FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg, FortranArray1DRef<double> rms, FortranArray1DRef<double> clw, FortranArray1DRef<double> cdw, FortranArray1DRef<double> cdpw, FortranArray1DRef<double> cdvw, FortranArray1DRef<double> cxw, FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw, FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw, FortranArray1DRef<double> cmzw, int& n_clcd, FortranArray3DRef<double> clcd, int& nblocks_clcd, FortranArray2DRef<int> blocks_clcd, FortranArray1DRef<double> chdw, FortranArray1DRef<double> swetw, FortranArray1DRef<double> fmdotw, FortranArray1DRef<double> cfttotw, FortranArray1DRef<double> cftmomw, FortranArray1DRef<double> cftpw, FortranArray1DRef<double> cftvw, FortranArray2DRef<double> rmstr, FortranArray2DRef<int> nneg, int& ntr, FortranArray2DRef<double> windex, int& ninter, FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt, FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy, FortranArray2DRef<double> dthetzz, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg, FortranArray1DRef<double> dxintg, FortranArray1DRef<double> dyintg, FortranArray1DRef<double> dzintg, FortranArray1DRef<int> iiig, FortranArray1DRef<int> jjig, FortranArray1DRef<int> kkig, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray1DRef<int> ireq_qb, FortranArray1DRef<int> igridg, FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae, FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae, FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass, FortranArray2DRef<double> damp, FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0, int& nmds, int& maxaes, FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb, FortranArray3DRef<int> islavept, int& nslave, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray3DRef<double> bmat, FortranArray3DRef<double> stm, FortranArray3DRef<double> stmi, FortranArray2DRef<double> gforcn, FortranArray2DRef<double> gforcnm, FortranArray2DRef<double> xxn, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg, int& nmaster, FortranArray4DRef<double> aehist, FortranArray1DRef<double> timekeep, FortranArray2DRef<int> inpl3d, int& nplots, int& nplot3d, FortranArray1DRef<int> levelg, FortranArray1DRef<int> iadvance, FortranArray2DRef<double> xs, FortranArray2DRef<double> gforcs, FortranArray2DRef<double> xorgae0, FortranArray2DRef<double> yorgae0, FortranArray2DRef<double> zorgae0, FortranArray2DRef<int> icouple, int& lfgm, FortranArray2DRef<int> nblk, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray2DRef<int> nblelst, FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax, FortranArray1DRef<int> kskmax, FortranArray1DRef<double> ue, FortranArray1DRef<int> irdrea, int& nbli, int& nummem)
{
    // COMMON block references
    int32_t& icgns     = cmn_cgns.icgns;
    int32_t& nzones    = cmn_cgns.nzones;
    int32_t& iccg      = cmn_cgns.iccg;
    int32_t& ibase     = cmn_cgns.ibase;
    int32_t& jdim      = cmn_ginfo.jdim;
    int32_t& kdim      = cmn_ginfo.kdim;
    int32_t& idim      = cmn_ginfo.idim;
    int32_t& jj2       = cmn_ginfo.jj2;
    int32_t& kk2       = cmn_ginfo.kk2;
    int32_t& ii2       = cmn_ginfo.ii2;
    int32_t& lq        = cmn_ginfo.lq;
    int32_t& lqj0      = cmn_ginfo.lqj0;
    int32_t& lqk0      = cmn_ginfo.lqk0;
    int32_t& lqi0      = cmn_ginfo.lqi0;
    int32_t& lsj       = cmn_ginfo.lsj;
    int32_t& lsk       = cmn_ginfo.lsk;
    int32_t& lsi       = cmn_ginfo.lsi;
    int32_t& lvol      = cmn_ginfo.lvol;
    int32_t& ldtj      = cmn_ginfo.ldtj;
    int32_t& lx        = cmn_ginfo.lx;
    int32_t& ly        = cmn_ginfo.ly;
    int32_t& lz        = cmn_ginfo.lz;
    int32_t& lvis      = cmn_ginfo.lvis;
    int32_t& lsnk0     = cmn_ginfo.lsnk0;
    int32_t& lsni0     = cmn_ginfo.lsni0;
    int32_t& lq1       = cmn_ginfo.lq1;
    int32_t& lqr       = cmn_ginfo.lqr;
    int32_t& lblk      = cmn_ginfo.lblk;
    int32_t& lxib      = cmn_ginfo.lxib;
    int32_t& lsig      = cmn_ginfo.lsig;
    int32_t& lsqtq     = cmn_ginfo.lsqtq;
    int32_t& lg        = cmn_ginfo.lg;
    int32_t& ltj0      = cmn_ginfo.ltj0;
    int32_t& ltk0      = cmn_ginfo.ltk0;
    int32_t& lti0      = cmn_ginfo.lti0;
    int32_t& lxkb      = cmn_ginfo.lxkb;
    int32_t& lnbl      = cmn_ginfo.lnbl;
    int32_t& lvj0      = cmn_ginfo.lvj0;
    int32_t& lvk0      = cmn_ginfo.lvk0;
    int32_t& lvi0      = cmn_ginfo.lvi0;
    int32_t& lbcj      = cmn_ginfo.lbcj;
    int32_t& lbck      = cmn_ginfo.lbck;
    int32_t& lbci      = cmn_ginfo.lbci;
    int32_t& lqc0      = cmn_ginfo.lqc0;
    int32_t& ldqc0     = cmn_ginfo.ldqc0;
    int32_t& lxtbi     = cmn_ginfo.lxtbi;
    int32_t& lxtbj     = cmn_ginfo.lxtbj;
    int32_t& lxtbk     = cmn_ginfo.lxtbk;
    int32_t& latbi     = cmn_ginfo.latbi;
    int32_t& latbj     = cmn_ginfo.latbj;
    int32_t& latbk     = cmn_ginfo.latbk;
    int32_t& lbcdj     = cmn_ginfo.lbcdj;
    int32_t& lbcdk     = cmn_ginfo.lbcdk;
    int32_t& lbcdi     = cmn_ginfo.lbcdi;
    int32_t& lxib2     = cmn_ginfo.lxib2;
    int32_t& lux       = cmn_ginfo.lux;
    int32_t& lcmuv     = cmn_ginfo.lcmuv;
    int32_t& lvolj0    = cmn_ginfo.lvolj0;
    int32_t& lvolk0    = cmn_ginfo.lvolk0;
    int32_t& lvoli0    = cmn_ginfo.lvoli0;
    int32_t& lxmdj     = cmn_ginfo.lxmdj;
    int32_t& lxmdk     = cmn_ginfo.lxmdk;
    int32_t& lxmdi     = cmn_ginfo.lxmdi;
    int32_t& lvelg     = cmn_ginfo.lvelg;
    int32_t& ldeltj    = cmn_ginfo.ldeltj;
    int32_t& ldeltk    = cmn_ginfo.ldeltk;
    int32_t& ldelti    = cmn_ginfo.ldelti;
    int32_t& lxnm2     = cmn_ginfo.lxnm2;
    int32_t& lynm2     = cmn_ginfo.lynm2;
    int32_t& lznm2     = cmn_ginfo.lznm2;
    int32_t& lxnm1     = cmn_ginfo.lxnm1;
    int32_t& lynm1     = cmn_ginfo.lynm1;
    int32_t& lznm1     = cmn_ginfo.lznm1;
    int32_t& lqavg     = cmn_ginfo.lqavg;
    int32_t& lq2avg    = cmn_ginfo2.lq2avg;
    float&   xmach     = cmn_info.xmach;
    float&   alpha     = cmn_info.alpha;
    float&   beta      = cmn_info.beta;
    float&   dt        = cmn_info.dt;
    int32_t& nit       = cmn_info.nit;
    int32_t& ntt       = cmn_info.ntt;
    int32_t* idiag     = cmn_info.idiag;
    int32_t& nitfo     = cmn_info.nitfo;
    int32_t& iflagts   = cmn_info.iflagts;
    int32_t* iflim     = cmn_info.iflim;
    int32_t& nres      = cmn_info.nres;
    int32_t* levelb    = cmn_info.levelb;
    int32_t& mgflag    = cmn_info.mgflag;
    int32_t& iconsf    = cmn_info.iconsf;
    int32_t& mseq      = cmn_info.mseq;
    int32_t* ncyc1     = cmn_info.ncyc1;
    int32_t* levelt    = cmn_info.levelt;
    int32_t* nitfo1    = cmn_info.nitfo1;
    int32_t& ngam      = cmn_info.ngam;
    int32_t* nsm       = cmn_info.nsm;
    int32_t& iipv      = cmn_info.iipv;
    int32_t& ivmx      = cmn_maxiv.ivmx;
    int32_t& ncyctot   = cmn_ncyct.ncyctot;
    float&   reue      = cmn_reyue.reue;
    float&   tinf      = cmn_reyue.tinf;
    int32_t* ivisc     = cmn_reyue.ivisc;
    int32_t& isklton   = cmn_sklton.isklton;
    int32_t& isminc    = cmn_sminn.isminc;
    int32_t& ismincforce = cmn_sminn.ismincforce;
    float&   time      = cmn_unst.time;
    float&   cfltau    = cmn_unst.cfltau;
    int32_t& ntstep    = cmn_unst.ntstep;
    int32_t& ita       = cmn_unst.ita;
    int32_t& iunst     = cmn_unst.iunst;
    int32_t& irest     = cmn_wrestq.irest;
    int32_t& irest2    = cmn_wrestq.irest2;
    char*    grid      = cmn_filenam.grid;
    char*    restrt    = cmn_filenam.restrt;
    int32_t& ip3dgrd   = cmn_igrdtyp.ip3dgrd;
    int32_t& ialph     = cmn_igrdtyp.ialph;
    int32_t& nnodes    = cmn_mydist2.nnodes;
    int32_t& myhost    = cmn_mydist2.myhost;
    int32_t& myid      = cmn_mydist2.myid;
    int32_t& mycomm    = cmn_mydist2.mycomm;
    int32_t& ndefrm    = cmn_elastic.ndefrm;
    int32_t& naesrf    = cmn_elastic.naesrf;
    int32_t& idef_ss   = cmn_elastic_ss.idef_ss;
    int32_t& negvol    = cmn_deformz.negvol;
    int32_t& ndgrd     = cmn_deformz.ndgrd;
    int32_t& irghost   = cmn_ghost.irghost;
    int32_t& iwghost   = cmn_ghost.iwghost;
    int32_t& irigb     = cmn_rigidbody.irigb;
    int32_t& irbtrim   = cmn_rigidbody.irbtrim;
    float&   xcentrot  = cmn_noninertial.xcentrot;
    float&   ycentrot  = cmn_noninertial.ycentrot;
    float&   zcentrot  = cmn_noninertial.zcentrot;
    int32_t& noninflag = cmn_noninertial.noninflag;
    float&   xnumavg   = cmn_avgdata.xnumavg;
    float&   xnumavg2  = cmn_avgdata.xnumavg2;
    int32_t& iteravg   = cmn_avgdata.iteravg;
    int32_t& ipertavg  = cmn_avgdata.ipertavg;
    int32_t& iclcd     = cmn_avgdata.iclcd;
    float&   roll_angle = cmn_gridtrans.roll_angle;
    char*    clcds     = cmn_filenam2.clcds;

    // Local variables
    double bcdat[12];
    int intpts[4];
    int isndrea[2];
    int mytag_qb[8];
    float time2mc = 0.0f;
    int ntq = 0;
    int igrid, nbl, ncg, iem;
    int iflag, icnt, iflagv;
    int imin, imax, jmin, jmax, kmin, kmax;
    int igeom_img;
    int irr, jdum, kdum, idum, ndum, ll;
    int lxc, lyc, lzc, lvolc;
    int lbcdjc, lbcdkc, lbcdic;
    int ldata, mdim, ndim, ldatac, mdimc, ndimc;
    char filname[80];
    int ibpnts, iipnts, nblpts;
    int inewg, iskipz, ncycchk, iflagg, iuns, iunsn;
    int idima, jdima, kdima;
    int ngridt, idimt, jdimt, kdimt;
    int iflg, ldim, need;
    int isksav, nttuse;
    int j1, j1sav, itry, iset1, iset2, ntot, izz;
    int nworkf, nworki, nworkxs, nworkixs, nwklsmin, nwkireq, nroomf, nroomi;
    int nadd, minbox, nbb, ntempf, ntempi;
    int ns2004, ns2014, ns2024, ns2034, ns2016, nsurf;
    int lwk1, lwk2, lwk3, lwk4, lwk5, lwk6;
    int iwk1, iwk2, iwk3, iwk4, iwk5, iwk6;
    int iwork1;
    int iunitr, iunitw, iunit;
    int iaes, nmodes, iaesrf, nm, is, m, n, nc, ns;
    int lll;
    int nblout, isk1, lxmdjc, lxmdkc, lxmdic, jst, kst, ist;
    int mxaedum, nmddum, nbldum;
    int iter, nn;
    int iset, nseg, ibeg, iend, jbeg, jend, kbeg, kend, llev;
    int iseqr, iupdat, iwk1i;
    int iseg;
    int nblc_local;
    int iuns_local;
    int mm;
    int basedesired[32]; // character*32 - not actually used as string in logic
    // termn8 error code
    int termn8_err = -1;


    // if (myid.eq.myhost) then
    if (myid == myhost) {
        if (icgns != 1) {
            nzones = ngrid;
            fortran_open_unit(2, restrt, "r+b");
            if (iclcd == 1 || iclcd == 2) {
                char clcds_tmp[242];
                int clen = 241;
                while (clen > 0 && clcds[clen-1] == ' ') clen--;
                std::memcpy(clcds_tmp, clcds, clen);
                clcds_tmp[clen] = '\0';
                fortran_open_unit(102, clcds_tmp, "r+b");
            }
        }
    }

    ntt      = 0;
    ntq      = 0;
    ntr      = 0;
    time     = 0.0f;
    time2mc  = 0.0f;
    xnumavg  = 0.0f;
    xnumavg2 = 0.0f;
    for (nbl = 1; nbl <= nblock; nbl++) {
        time2(nbl) = 0.0;
    }
    for (int nc2 = 1; nc2 <= ncycmax; nc2++) {
        timekeep(nc2) = 0.0;
        if (naesrf > 0) {
            for (int iaes2 = 1; iaes2 <= naesrf; iaes2++) {
                nmodes = (int)aesrfdat(5, iaes2);
                for (int nm2 = 1; nm2 <= nmodes; nm2++) {
                    aehist(nc2, 1, nm2, iaes2) = 0.0;
                    aehist(nc2, 2, nm2, iaes2) = 0.0;
                    aehist(nc2, 3, nm2, iaes2) = 0.0;
                }
            }
        }
    }
    if (iunst > 1) {
        for (int ns2 = 1; ns2 <= 3*nslave; ns2++) {
            ue(ns2) = 0.0;
        }
    }

    lig(1) = 1;
    lbg(1) = 1;

    if (myid == mblk2nd(1)) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "***** BEGINNING INITIALIZATION *****");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
    }

    igeom_img = 0;

    if (icgns == 1 && myid == myhost) {
        // open the CGNS grid database (read-only); sets iccg, ibase, nzones.
        // 'grid' is an 80-char space-padded Fortran field -> trim + NUL-terminate.
        char cgfname[81];
        std::memcpy(cgfname, grid, 80); cgfname[80] = '\0';
        for (int t = 79; t >= 0 && (cgfname[t] == ' ' || cgfname[t] == '\0'); t--) cgfname[t] = '\0';
        int idimdesired = 3;   // 3-D structured index dimension
        int rc = cgnstools_ns::ropencgns_c(cgfname, idimdesired, iccg, ibase, nzones);
        if (rc != 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " error opening cgns grid file %s", grid);
            termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
        }
        if (nzones != ngrid) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                " inconsistency between cgns # of zones & input. nzones,ngrid=%6d%6d",
                nzones, ngrid);
            termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
        }
    }

    // do 1000 igrid=1,ngrid
    for (igrid = 1; igrid <= ngrid; igrid++) {
        iem = iemg(igrid);
        nbl = nblg(igrid);
        lead_ns::lead(nbl, lw, lw2, maxbl);

        if (icgns == 1) {
            // read this zone's grid coordinates from the CGNS database
            if (myid == myhost) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " reading grid%6d of dimensions (I/J/K) :%6d%6d%6d",
                    igrid, idim, jdim, kdim);
                int rc = cgnstools_ns::getgrid_c(iccg, ibase, igrid, idim, jdim, kdim, ialph,
                    &w(lx), &w(ly), &w(lz));
                if (rc != 0) {
                    FILE* f11 = fortran_get_unit(11);
                    std::fprintf(f11, " Stopping... error reading cgns grid for zone%6d\n", igrid);
                    termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                }
            }
        } else {
            if (ip3dgrd == 0) {
                // cfl3d type
                if (myid == myhost) {
                    FILE* f1 = fortran_get_unit(1);
                    // read(1) jdum,kdum,idum  (unformatted sequential record)
                    FortranRecordReader rec(f1);
                    rec.read_scalar(jdum);
                    rec.read_scalar(kdum);
                    rec.read_scalar(idum);
                    rec.finish();
                }
                if (myid == mblk2nd(nbl)) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " reading grid%6d of dimensions (I/J/K) :%6d%6d%6d",
                        igrid, idim, jdim, kdim);
                }
                if (myid == myhost) {
                    if (jdum != jdim || kdum != kdim || idum != idim) {
                        FILE* f11 = fortran_get_unit(11);
                        std::fprintf(f11, " stopping....inconsistency in\n");
                        std::fprintf(f11, " grid data file\n");
                        std::fprintf(f11, " for grid%6d\n", igrid);
                        std::fprintf(f11, " input file: idim,jdim,kdim%6d%6d%6d\n", idim, jdim, kdim);
                        std::fprintf(f11, " grid  file: idim,jdim,kdim%6d%6d%6d\n", idum, jdum, kdum);
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    if (ialph != 0 && iipv == 1) {
                        FILE* f11 = fortran_get_unit(11);
                        std::fprintf(f11, " To use farfield point vortex, grid must be in x-z plane\n");
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    irr = 0;
                    rcfl_ns::rcfl(FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                                  jdim, kdim, idim, igrid, ialph, igeom_img, irr);
                    if (irr != 0) {
                        FILE* f11 = fortran_get_unit(11);
                        std::fprintf(f11, " Stopping... error reading grid...\n");
                        std::fprintf(f11, " (Common error:  grid not written in same precision that CFL3D was compiled)\n");
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                }
            } else {
                // plot3d type
                if (myid == myhost) {
                    if (igrid == 1) {
                        FILE* f1 = fortran_get_unit(1);
                        // read(1) ndum  — a Fortran UNFORMATTED record: the int
                        // payload is wrapped by 4-byte record-length markers, so
                        // it must be read marker-aware (a bare fread would return
                        // the leading marker instead of ndum).
                        { FortranRecordReader rr(f1); rr.read_scalar(ndum); rr.finish(); }
                        if (ndum != ngrid) {
                            FILE* f11 = fortran_get_unit(11);
                            std::fprintf(f11, " stopping....ngrid = %6d but grid file contains %6d grids\n", ngrid, ndum);
                            termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                        }
                        // read(1) (itest(ll),jtest(ll),ktest(ll),ll=1,ngrid) — one record
                        { FortranRecordReader rr(f1);
                          for (ll = 1; ll <= ngrid; ll++) {
                              rr.read_scalar(itest(ll));
                              rr.read_scalar(jtest(ll));
                              rr.read_scalar(ktest(ll));
                          }
                          rr.finish(); }
                    }
                }
                if (myid == mblk2nd(nbl)) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " reading grid%6d of dimensions (I/J/K) :%6d%6d%6d",
                        igrid, idim, jdim, kdim);
                }

                if (myid == myhost) {
                    if (jtest(igrid) != jdim || ktest(igrid) != kdim || itest(igrid) != idim) {
                        FILE* f11 = fortran_get_unit(11);
                        std::fprintf(f11, " stopping....inconsistency in grid data file\n");
                        std::fprintf(f11, " for grid%6d\n", igrid);
                        std::fprintf(f11, " input file: idim,jdim,kdim%6d%6d%6d\n", idim, jdim, kdim);
                        std::fprintf(f11, " grid  file: idim,jdim,kdim%6d%6d%6d\n",
                                     itest(igrid), jtest(igrid), ktest(igrid));
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    irr = 0;
                    rp3d_ns::rp3d(FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                                  jdim, kdim, idim, igrid, ialph, igeom_img, irr);
                    if (irr != 0) {
                        FILE* f11 = fortran_get_unit(11);
                        std::fprintf(f11, " Stopping... error reading grid...\n");
                        std::fprintf(f11, " (Common error:  grid not written in same precision that CFL3D was compiled)\n");
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                }
            } // end ip3dgrd else
        } // end icgns else

        if (myid == myhost) {
            if (igrid == 1) {
                FILE* f11 = fortran_get_unit(11);
                std::fprintf(f11, "\n");
                std::fprintf(f11, "input roll angle = %9.4f degrees\n", (float)roll_angle);
                roll_angle = roll_angle * (float)(std::atan(1.0)/45.0);
            }
            if (std::fabs((float)roll_angle) > 0.0f) {
                nbl = nblg(igrid);
                { double gm_xc1 = xcentrot, gm_yc1 = ycentrot, gm_zc1 = zcentrot;
                  double gm_xc2 = xcentrot, gm_yc2 = ycentrot, gm_zc2 = zcentrot;
                  double gm_rx = roll_angle, gm_ry = 0.0, gm_rz = 0.0;
                grdmove_ns::grdmove(nbl, jdim, kdim, idim,
                    FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                    gm_xc1, gm_yc1, gm_zc1,
                    gm_xc2, gm_yc2, gm_zc2,
                    gm_rx, gm_ry, gm_rz); }
            }
        }

        // store off original grid if mesh deforms
        if (idefrm(nbl) > 0) {
            if (myid == mblk2nd(nbl)) {
                for (lll = lx; lll <= lvis-1; lll++) {
                    w(lxnm1 + lll - lx) = w(lll);
                }
            }
        }

        // collocate grid points on coarser levels
        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (m = 1; m <= ncg; m++) {
                nbl = nbl + 1;
                if (myid == mblk2nd(nbl)) {
                    lxc = lw(10, nbl);
                    lyc = lw(11, nbl);
                    lzc = lw(12, nbl);
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "   creating coarser block%6d of dimensions (I/J/K) :%6d%6d%6d",
                        nbl, ii2, jj2, kk2);
                    collx_ns::collx(FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                                    FortranArray3DRef<double>(&w(lxc), jj2, kk2, ii2), FortranArray3DRef<double>(&w(lyc), jj2, kk2, ii2), FortranArray3DRef<double>(&w(lzc), jj2, kk2, ii2),
                                    jdim, kdim, idim, jj2, kk2, ii2);
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    if (idefrm(nbl) > 0) {
                        for (lll = lx; lll <= lvis-1; lll++) {
                            w(lxnm1 + lll - lx) = w(lll);
                        }
                    }
                }
            }
        }

        { int iunit11 = 11; writ_buf_ns::writ_buf(nblout, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
    } // end do 1000

    // initialize slavept data for mesh deformation
    if (iunst > 1 || idef_ss > 0) {
        setslave_ns::setslave(lw, lw2, w, mgwk, wk, nwork, maxbl, mxbli, maxgr, maxseg,
                              nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
                              ibcinfo, jbcinfo, kbcinfo, nblock, idefrm,
                              iadvance, nou, bou, nbuf, ibufdim, myid, myhost,
                              mycomm, mblk2nd, icsi, icsf, jcsi, jcsf, kcsi,
                              kcsf, islavept, nslave, nsegdfrm, idfrmseg,
                              maxsegdg, iwk, iwork, nmaster, ngrid, jskip, kskip,
                              iskip, nblg, levelg, lfgm, nblk, limblk, isva, nblelst,
                              nnodes, iskmax, jskmax, kskmax, nbli);
    }


    if (idef_ss > 0) {
        // read in new surface grid and deform mesh to reflect the new surface
        iunitr = 39;
        fortran_open_unit(iunitr, "newsurf.p3d", "r");
        // read(iunitr,*) ndum
        fortran_read_list(iunitr, &ndum);
        // read(iunitr,*) (idum,jdum,kdum,ll=1,ndum) - reads 3*ndum values from one record
        // Use a temporary array to read all values at once
        {
            FortranArray1D<int> tmp_dims(3*ndum);
            fortran_read_list_array(iunitr, tmp_dims, 1, 3*ndum);
            for (ll = 1; ll <= ndum; ll++) {
                idum = tmp_dims(3*(ll-1)+1);
                jdum = tmp_dims(3*(ll-1)+2);
                kdum = tmp_dims(3*(ll-1)+3);
            }
        }

        // read in the new surface data and temporarily store in deltj/deltk/delti arrays
        for (igrid = 1; igrid <= ngrid; igrid++) {
            nbl = nblg(igrid);
            if (idefrm(nbl) == 1) {
                lead_ns::lead(nbl, lw, lw2, maxbl);
                for (iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {
                    rsurf_ns::rsurf(maxbl, maxsegdg, idim, jdim, kdim,
                                    FortranArray4DRef<double>(&w(ldelti), jdim, kdim, 3, 2), FortranArray4DRef<double>(&w(ldeltj), kdim, idim, 3, 2), FortranArray4DRef<double>(&w(ldeltk), jdim, idim, 3, 2),
                                    nbl, icsi, icsf, jcsi, jcsf, kcsi, kcsf,
                                    iseg, nou, bou, nbuf, ibufdim, iunitr);
                }
                // calculate delta displacement between new surface and current surface
                lwk1 = 1;
                lwk2 = lwk1 + kdim*idim*2;
                lwk3 = lwk2 + jdim*idim*2;
                lwk4 = lwk3 + jdim*kdim*2;
                if (nwork < lwk4) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " not enough work space for subroutine getdelt");
                    termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                }
                for (lll = lwk1; lll <= lwk4; lll++) {
                    wk(lll) = 1.0;
                }
                for (iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {
                    getdelt_ns::getdelt(maxbl, maxsegdg, idim, jdim, kdim,
                                        FortranArray4DRef<double>(&w(ldelti), jdim, kdim, 3, 2), FortranArray4DRef<double>(&w(ldeltj), kdim, idim, 3, 2), FortranArray4DRef<double>(&w(ldeltk), jdim, idim, 3, 2),
                                        FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                                        nbl, icsi, icsf, jcsi, jcsf, kcsi, kcsf,
                                        iseg, nou, bou, nbuf, ibufdim,
                                        FortranArray3DRef<double>(&wk(lwk1), kdim, idim, 2), FortranArray3DRef<double>(&wk(lwk2), jdim, idim, 2), FortranArray3DRef<double>(&wk(lwk3), jdim, kdim, 2));
                }
            }
        }

        // deform volume grid to fit new surface
        iwk1i = 2*maxbl;
        iwork1 = iwork - iwk1i;
        if (iwork1 < 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "stopping...not enough integer work space for subroutine updatedg");
            termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
        }
        for (int iii = 1; iii <= iwk1i; iii++) {
            iwk(iii) = 0;
        }
        iseqr  = mseq;
        iupdat = 1;
        {
            int nt_local = 1;
            FortranArray1DRef<int> iwk_slice = FortranArray1DRef<int>(&iwk(maxbl), iwk.size() - (maxbl) + 1);
            updatedg(lw, lw2, w, mgwk, wk, nwork, iupdat, iseqr, maxbl,
                     maxgr, maxseg, nbci0, nbcj0, nbck0, nbcidim,
                     nbcjdim, nbckdim, ibcinfo, jbcinfo, kbcinfo,
                     nblock, levelg, igridg, idefrm, ncgg, iadvance, nou,
                     bou, nbuf, ibufdim, myid, myhost, mycomm, mblk2nd,
                     utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae,
                     xorgae, yorgae, zorgae, thtxae, thtyae, thtzae,
                     rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf,
                     kcsi, kcsf, freq, gmass, damp, x0, gf0, nmds, maxaes,
                     aesrfdat, perturb, itrans, irotat, islavept, nslave,
                     iskip, jskip, kskip, xs, xxn, nsegdfrm, idfrmseg,
                     iaesurf, maxsegdg, iwk, nmaster, nt_local, xorig, yorig,
                     zorig, xorgae0, yorgae0, zorgae0, icouple,
                     iwk_slice, nnodes, nblelst, iskmax, jskmax, kskmax, ue);
        }

        // collocate deformed grid points on coarser levels
        for (igrid = 1; igrid <= ngrid; igrid++) {
            nbl = nblg(igrid);
            lead_ns::lead(nbl, lw, lw2, maxbl);
            ncg = ncgg(igrid);
            if (ncg > 0) {
                for (m = 1; m <= ncg; m++) {
                    nbl = nbl + 1;
                    if (myid == mblk2nd(nbl)) {
                        lxc = lw(10, nbl);
                        lyc = lw(11, nbl);
                        lzc = lw(12, nbl);
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "   creating coarser block%6d of dimensions (I/J/K) :%6d%6d%6d",
                            nbl, ii2, jj2, kk2);
                        collx_ns::collx(FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                                        FortranArray3DRef<double>(&w(lxc), jj2, kk2, ii2), FortranArray3DRef<double>(&w(lyc), jj2, kk2, ii2), FortranArray3DRef<double>(&w(lzc), jj2, kk2, ii2),
                                        jdim, kdim, idim, jj2, kk2, ii2);
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                    }
                }
            }
        }
    } // end idef_ss > 0


    // compute grid metrics
    if (myid == mblk2nd(1)) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
    }
    iflag = 0;
    icnt  = 0;
    // do 1200 igrid=1,ngrid
    for (igrid = 1; igrid <= ngrid; igrid++) {
        iem = iemg(igrid);
        nbl = nblg(igrid);
        if (myid == mblk2nd(nbl)) {
            lead_ns::lead(nbl, lw, lw2, maxbl);
            if (nwork < jdim*kdim*6 + jdim*kdim*idim*5) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " not enough work space for subroutine metric (must increase mwork)");
                termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
            }
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                " computing metrics for block%6d  (grid%6d)", nbl, igrid);
            metric_ns::metric(jdim, kdim, idim, FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                   FortranArray3DRef<double>(&w(lsj), jdim*kdim, idim-1, 5), FortranArray3DRef<double>(&w(lsk), jdim*kdim, idim-1, 5), FortranArray3DRef<double>(&w(lsi), jdim*kdim, idim, 5),
                   FortranArray2DRef<double>(&wk(1), jdim*kdim, 6), FortranArray3DRef<double>(&wk(jdim*kdim*6+1), jdim*kdim, idim, 5),
                   nbl, iflag, icnt,
                   nbci0, nbcj0, nbck0, nbcidim, nbcjdim,
                   nbckdim, ibcinfo, jbcinfo, kbcinfo, maxbl, maxseg,
                   nou, bou, nbuf, ibufdim, myid, mblk2nd);
        }

        // compute grid metrics on coarser levels
        ncg = ncgg(igrid);
        iem = iemg(igrid);
        if (ncg > 0) {
            isk1    = isklton;
            isklton = 2;
            for (m = 1; m <= ncg; m++) {
                nbl = nbl + 1;
                if (mblk2nd(nbl) == myid) {
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "   computing metrics for coarser block%6d", nbl);
                    metric_ns::metric(jdim, kdim, idim, FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                           FortranArray3DRef<double>(&w(lsj), jdim*kdim, idim-1, 5), FortranArray3DRef<double>(&w(lsk), jdim*kdim, idim-1, 5), FortranArray3DRef<double>(&w(lsi), jdim*kdim, idim, 5),
                           FortranArray2DRef<double>(&wk(1), jdim*kdim, 6), FortranArray3DRef<double>(&wk(jdim*kdim*6+1), jdim*kdim, idim, 5),
                           nbl, iflag, icnt,
                           nbci0, nbcj0, nbck0, nbcidim, nbcjdim,
                           nbckdim, ibcinfo, jbcinfo, kbcinfo, maxbl, maxseg,
                           nou, bou, nbuf, ibufdim, myid, mblk2nd);
                }
            }
            isklton = isk1;
        }

        nblout = nblg(igrid);
        { int iunit11 = 11; writ_buf_ns::writ_buf(nblout, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
    } // end do 1200

    if (iflag > 0 && myid == myhost) {
        FILE* f11 = fortran_get_unit(11);
        std::fprintf(f11, " Fatal error(s) uncovered in grid metrics\n");
        if (negvol == 0) {
            termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
        }
    }

    // compute cell volumes
    if (myid == mblk2nd(1)) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
    }

    // do 1300 igrid=1,ngrid
    for (igrid = 1; igrid <= ngrid; igrid++) {
        iem = iemg(igrid);
        nbl = nblg(igrid);
        if (myid == mblk2nd(nbl)) {
            lead_ns::lead(nbl, lw, lw2, maxbl);
            if (nwork < jdim*kdim*15) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " not enough work space for subroutine cellvol (must increase mwork)");
                termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
            }
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                " computing cell volumes for block%6d  (grid%6d)", nbl, igrid);
            {
                int iflagv1_local = 0;
                cellvol_ns::cellvol(jdim, kdim, idim, FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                        FortranArray3DRef<double>(&w(lsj), jdim*kdim, idim-1, 5), FortranArray3DRef<double>(&w(lsk), jdim*kdim, idim-1, 5), FortranArray3DRef<double>(&w(lsi), jdim*kdim, idim, 5), FortranArray3DRef<double>(&w(lvol), jdim, kdim, idim-1),
                        FortranArray2DRef<double>(&wk(1), jdim*kdim, 15), nou, bou, nbuf, ibufdim,
                        myid, mblk2nd, maxbl, nbl, iflagv1_local, iflagv,
                        imin, imax, jmin, jmax, kmin, kmax);
            }
        }

        // collocate cell volumes on coarser levels
        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (m = 1; m <= ncg; m++) {
                nbl = nbl + 1;
                if (myid == mblk2nd(nbl)) {
                    lvolc = lw(8, nbl);
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "   computing cell volumes for coarser block%6d", nbl);
                    collv_ns::collv(FortranArray3DRef<double>(&w(lvol), jdim, kdim, idim-1), FortranArray3DRef<double>(&w(lvolc), jj2, kk2, ii2-1),
                                    jdim, kdim, idim, jj2, kk2, ii2);
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                }
            }
        }

        nblout = nblg(igrid);
        { int iunit11 = 11; writ_buf_ns::writ_buf(nblout, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
    } // end do 1300


    // read/set auxilliary data arrays for 2000 series boundary conditions
    if (myid == mblk2nd(1)) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
    }

    // do 1400 igrid=1,ngrid
    for (igrid = 1; igrid <= ngrid; igrid++) {
        iem = iemg(igrid);
        nbl = nblg(igrid);
        lead_ns::lead(nbl, lw, lw2, maxbl);

        // do 1420 nseg=1,nbci0(nbl)
        for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
            mm = 1;
            if (std::abs(ibcinfo(nbl,nseg,1,mm)) >= 2000 &&
                std::abs(ibcinfo(nbl,nseg,1,mm)) < 3000) {
                ldata = lwdat(nbl, nseg, 1);
                mdim = ibcinfo(nbl,nseg,3,mm) - ibcinfo(nbl,nseg,2,mm);
                ndim = ibcinfo(nbl,nseg,5,mm) - ibcinfo(nbl,nseg,4,mm);
                std::memcpy(filname, bcfiles(bcfilei(nbl,nseg,mm)), 80);
                for (ll = 1; ll <= 12; ll++) bcdat[ll-1] = bcvali(nbl,nseg,ll,mm);
                readdat_ns::readdat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim, filname, FortranArray1DRef<double>(bcdat, 12),
                                    nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl);
            }
        }

        // do 1430 nseg=1,nbcidim(nbl)
        for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
            mm = 2;
            if (std::abs(ibcinfo(nbl,nseg,1,mm)) >= 2000 &&
                std::abs(ibcinfo(nbl,nseg,1,mm)) < 3000) {
                ldata = lwdat(nbl, nseg, 2);
                mdim = ibcinfo(nbl,nseg,3,mm) - ibcinfo(nbl,nseg,2,mm);
                ndim = ibcinfo(nbl,nseg,5,mm) - ibcinfo(nbl,nseg,4,mm);
                std::memcpy(filname, bcfiles(bcfilei(nbl,nseg,mm)), 80);
                for (ll = 1; ll <= 12; ll++) bcdat[ll-1] = bcvali(nbl,nseg,ll,mm);
                readdat_ns::readdat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim, filname, FortranArray1DRef<double>(bcdat, 12),
                                    nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl);
            }
        }

        // do 1440 nseg=1,nbcj0(nbl)
        for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
            mm = 1;
            if (std::abs(jbcinfo(nbl,nseg,1,mm)) >= 2000 &&
                std::abs(jbcinfo(nbl,nseg,1,mm)) < 3000) {
                ldata = lwdat(nbl, nseg, 3);
                mdim = jbcinfo(nbl,nseg,5,mm) - jbcinfo(nbl,nseg,4,mm);
                ndim = jbcinfo(nbl,nseg,3,mm) - jbcinfo(nbl,nseg,2,mm);
                std::memcpy(filname, bcfiles(bcfilej(nbl,nseg,mm)), 80);
                for (ll = 1; ll <= 12; ll++) bcdat[ll-1] = bcvalj(nbl,nseg,ll,mm);
                readdat_ns::readdat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim, filname, FortranArray1DRef<double>(bcdat, 12),
                                    nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl);
            }
        }

        // do 1450 nseg=1,nbcjdim(nbl)
        for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
            mm = 2;
            if (std::abs(jbcinfo(nbl,nseg,1,mm)) >= 2000 &&
                std::abs(jbcinfo(nbl,nseg,1,mm)) < 3000) {
                ldata = lwdat(nbl, nseg, 4);
                mdim = jbcinfo(nbl,nseg,5,mm) - jbcinfo(nbl,nseg,4,mm);
                ndim = jbcinfo(nbl,nseg,3,mm) - jbcinfo(nbl,nseg,2,mm);
                std::memcpy(filname, bcfiles(bcfilej(nbl,nseg,mm)), 80);
                for (ll = 1; ll <= 12; ll++) bcdat[ll-1] = bcvalj(nbl,nseg,ll,mm);
                readdat_ns::readdat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim, filname, FortranArray1DRef<double>(bcdat, 12),
                                    nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl);
            }
        }

        // do 1460 nseg=1,nbck0(nbl)
        for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
            mm = 1;
            if (std::abs(kbcinfo(nbl,nseg,1,mm)) >= 2000 &&
                std::abs(kbcinfo(nbl,nseg,1,mm)) < 3000) {
                ldata = lwdat(nbl, nseg, 5);
                mdim = kbcinfo(nbl,nseg,5,mm) - kbcinfo(nbl,nseg,4,mm);
                ndim = kbcinfo(nbl,nseg,3,mm) - kbcinfo(nbl,nseg,2,mm);
                std::memcpy(filname, bcfiles(bcfilek(nbl,nseg,mm)), 80);
                for (ll = 1; ll <= 12; ll++) bcdat[ll-1] = bcvalk(nbl,nseg,ll,mm);
                readdat_ns::readdat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim, filname, FortranArray1DRef<double>(bcdat, 12),
                                    nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl);
            }
        }

        // do 1470 nseg=1,nbckdim(nbl)
        for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
            mm = 2;
            if (std::abs(kbcinfo(nbl,nseg,1,mm)) >= 2000 &&
                std::abs(kbcinfo(nbl,nseg,1,mm)) < 3000) {
                ldata = lwdat(nbl, nseg, 6);
                mdim = kbcinfo(nbl,nseg,5,mm) - kbcinfo(nbl,nseg,4,mm);
                ndim = kbcinfo(nbl,nseg,3,mm) - kbcinfo(nbl,nseg,2,mm);
                std::memcpy(filname, bcfiles(bcfilek(nbl,nseg,mm)), 80);
                for (ll = 1; ll <= 12; ll++) bcdat[ll-1] = bcvalk(nbl,nseg,ll,mm);
                readdat_ns::readdat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim, filname, FortranArray1DRef<double>(bcdat, 12),
                                    nou, bou, nbuf, ibufdim, myid, mblk2nd, maxbl);
            }
        }


        // collocate auxilliary data on coarser grids
        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (m = 1; m <= ncg; m++) {
                nbl = nbl + 1;
                if (mblk2nd(nbl) == myid) {
                    lbcdjc = lw(42, nbl);
                    lbcdkc = lw(43, nbl);
                    lbcdic = lw(44, nbl);

                    for (nseg = 1; nseg <= nbci0(nbl-1); nseg++) {
                        mm = 1;
                        if (std::abs(ibcinfo(nbl,nseg,1,mm)) >= 2000 &&
                            std::abs(ibcinfo(nbl,nseg,1,mm)) < 3000) {
                            ldata  = lwdat(nbl-1, nseg, 1);
                            mdim   = ibcinfo(nbl-1,nseg,3,mm) - ibcinfo(nbl-1,nseg,2,mm);
                            ndim   = ibcinfo(nbl-1,nseg,5,mm) - ibcinfo(nbl-1,nseg,4,mm);
                            ldatac = lwdat(nbl, nseg, 1);
                            mdimc  = ibcinfo(nbl,nseg,3,mm) - ibcinfo(nbl,nseg,2,mm);
                            ndimc  = ibcinfo(nbl,nseg,5,mm) - ibcinfo(nbl,nseg,4,mm);
                            colldat_ns::colldat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim,
                                                FortranArray4DRef<double>(&w(ldatac), mdimc, ndimc, 2, 12), mdimc, ndimc);
                        }
                    }

                    for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
                        mm = 2;
                        if (std::abs(ibcinfo(nbl,nseg,1,mm)) >= 2000 &&
                            std::abs(ibcinfo(nbl,nseg,1,mm)) < 3000) {
                            ldata  = lwdat(nbl-1, nseg, 2);
                            mdim   = ibcinfo(nbl-1,nseg,3,mm) - ibcinfo(nbl-1,nseg,2,mm);
                            ndim   = ibcinfo(nbl-1,nseg,5,mm) - ibcinfo(nbl-1,nseg,4,mm);
                            ldatac = lwdat(nbl, nseg, 2);
                            mdimc  = ibcinfo(nbl,nseg,3,mm) - ibcinfo(nbl,nseg,2,mm);
                            ndimc  = ibcinfo(nbl,nseg,5,mm) - ibcinfo(nbl,nseg,4,mm);
                            colldat_ns::colldat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim,
                                                FortranArray4DRef<double>(&w(ldatac), mdimc, ndimc, 2, 12), mdimc, ndimc);
                        }
                    }

                    for (nseg = 1; nseg <= nbcj0(nbl-1); nseg++) {
                        mm = 1;
                        if (std::abs(jbcinfo(nbl,nseg,1,mm)) >= 2000 &&
                            std::abs(jbcinfo(nbl,nseg,1,mm)) < 3000) {
                            ldata  = lwdat(nbl-1, nseg, 3);
                            mdim   = jbcinfo(nbl-1,nseg,5,mm) - jbcinfo(nbl-1,nseg,4,mm);
                            ndim   = jbcinfo(nbl-1,nseg,3,mm) - jbcinfo(nbl-1,nseg,2,mm);
                            ldatac = lwdat(nbl, nseg, 3);
                            mdimc  = jbcinfo(nbl,nseg,5,mm) - jbcinfo(nbl,nseg,4,mm);
                            ndimc  = jbcinfo(nbl,nseg,3,mm) - jbcinfo(nbl,nseg,2,mm);
                            colldat_ns::colldat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim,
                                                FortranArray4DRef<double>(&w(ldatac), mdimc, ndimc, 2, 12), mdimc, ndimc);
                        }
                    }

                    for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
                        mm = 2;
                        if (std::abs(jbcinfo(nbl,nseg,1,mm)) >= 2000 &&
                            std::abs(jbcinfo(nbl,nseg,1,mm)) < 3000) {
                            ldata  = lwdat(nbl-1, nseg, 4);
                            mdim   = jbcinfo(nbl-1,nseg,5,mm) - jbcinfo(nbl-1,nseg,4,mm);
                            ndim   = jbcinfo(nbl-1,nseg,3,mm) - jbcinfo(nbl-1,nseg,2,mm);
                            ldatac = lwdat(nbl, nseg, 4);
                            mdimc  = jbcinfo(nbl,nseg,5,mm) - jbcinfo(nbl,nseg,4,mm);
                            ndimc  = jbcinfo(nbl,nseg,3,mm) - jbcinfo(nbl,nseg,2,mm);
                            colldat_ns::colldat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim,
                                                FortranArray4DRef<double>(&w(ldatac), mdimc, ndimc, 2, 12), mdimc, ndimc);
                        }
                    }

                    for (nseg = 1; nseg <= nbck0(nbl-1); nseg++) {
                        mm = 1;
                        if (std::abs(kbcinfo(nbl,nseg,1,mm)) >= 2000 &&
                            std::abs(kbcinfo(nbl,nseg,1,mm)) < 3000) {
                            ldata  = lwdat(nbl-1, nseg, 5);
                            mdim   = kbcinfo(nbl-1,nseg,5,mm) - kbcinfo(nbl-1,nseg,4,mm);
                            ndim   = kbcinfo(nbl-1,nseg,3,mm) - kbcinfo(nbl-1,nseg,2,mm);
                            ldatac = lwdat(nbl, nseg, 5);
                            mdimc  = kbcinfo(nbl,nseg,5,mm) - kbcinfo(nbl,nseg,4,mm);
                            ndimc  = kbcinfo(nbl,nseg,3,mm) - kbcinfo(nbl,nseg,2,mm);
                            colldat_ns::colldat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim,
                                                FortranArray4DRef<double>(&w(ldatac), mdimc, ndimc, 2, 12), mdimc, ndimc);
                        }
                    }

                    for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
                        mm = 2;
                        if (std::abs(kbcinfo(nbl,nseg,1,mm)) >= 2000 &&
                            std::abs(kbcinfo(nbl,nseg,1,mm)) < 3000) {
                            ldata  = lwdat(nbl-1, nseg, 6);
                            mdim   = kbcinfo(nbl-1,nseg,5,mm) - kbcinfo(nbl-1,nseg,4,mm);
                            ndim   = kbcinfo(nbl-1,nseg,3,mm) - kbcinfo(nbl-1,nseg,2,mm);
                            ldatac = lwdat(nbl, nseg, 6);
                            mdimc  = kbcinfo(nbl,nseg,5,mm) - kbcinfo(nbl,nseg,4,mm);
                            ndimc  = kbcinfo(nbl,nseg,3,mm) - kbcinfo(nbl,nseg,2,mm);
                            colldat_ns::colldat(FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12), mdim, ndim,
                                                FortranArray4DRef<double>(&w(ldatac), mdimc, ndimc, 2, 12), mdimc, ndimc);
                        }
                    }
                }
            }
        }
    } // end do 1400


    // read overlapped grid interpolation data and set iblank arrays
    // do 1600 igrid=1,ngrid
    for (igrid = 1; igrid <= ngrid; igrid++) {
        iem = iemg(igrid);
        nbl = nblg(igrid);
        lead_ns::lead(nbl, lw, lw2, maxbl);
        if (iovrlp(nbl) == 1) {
            if (myid == mblk2nd(nbl)) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " reading overlap information for block%6d  (grid%6d)", nbl, igrid);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " reading overlap information for block %6d", nbl);
            }
            if (myid == myhost) {
                if (nbl > 1) {
                    if (iovrlp(nbl-1) == 0) {
                        lig(nbl) = lig(nbl-1);
                        lbg(nbl) = lbg(nbl-1);
                    }
                }
                getibk_ns::getibk(FortranArray3DRef<double>(&w(lblk), jdim, kdim, idim), jdim, kdim, idim, nbl, FortranArray1DRef<int>(intpts, 4),
                                  nblpts, ibpnts, iipnts, iitot, maxbl,
                                  iibg, kkbg, jjbg, ibcg, lig, lbg,
                                  dxintg, dyintg, dzintg, iiig, jjig, kkig,
                                  ibpntsg, iipntsg, myid, ibufdim, nbuf, bou, nou);
            }
            if (myid == mblk2nd(nbl)) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "  ibpnts,intpts,iipnts = %8d%8d%8d%8d%8d%8d",
                    ibpnts, intpts[0], intpts[1], intpts[2], intpts[3], iipnts);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "  counting iblank values....");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "  number of blank values .ne. 1 are %8d", nblpts);
            }
        } else {
            if (myid == mblk2nd(nbl)) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " setting blank array to 1.0 for block %6d  (grid%6d)", nbl, igrid);
            }
            if (myid == myhost) {
                setblk_ns::setblk(FortranArray3DRef<double>(&w(lblk), jdim, kdim, idim), jdim, kdim, idim, nbl);
                if (nbl > 1) {
                    lig(nbl) = lig(nbl-1);
                    lbg(nbl) = lbg(nbl-1);
                }
            }
        }

        // set iblank for coarser grids
        ncg = ncgg(igrid);
        if (ncg > 0) {
            for (m = 1; m <= ncg; m++) {
                nbl = nbl + 1;
                int lblkc = lw(18, nbl);
                if (iovrlp(nbl-1) == 0) {
                    lig(nbl) = lig(nbl-1);
                    lbg(nbl) = lbg(nbl-1);
                }
                if (mblk2nd(nbl) == myid) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "   setting blank array to 1.0 for coarser block%6d", nbl);
                    setblk_ns::setblk(FortranArray3DRef<double>(&w(lblkc), jj2, kk2, ii2), jj2, kk2, ii2, nbl);
                }
                lead_ns::lead(nbl, lw, lw2, maxbl);
            }
        }

        nblout = nblg(igrid);
        { int iunit11 = 11; writ_buf_ns::writ_buf(nblout, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
    } // end do 1600

    // read patched-grid interpolation data from file
    if (ninter < 0) {
        if (myid == myhost) {
            rpatch_ns::rpatch(maxbl, maxxe, intmax, nsub1, windex, ninter, iindex,
                              nblkpt, dthetxx, dthetyy, dthetzz, jdimg, kdimg,
                              idimg, nou, bou, nbuf, ibufdim, myid);
        }
    }


    // set initial conditions as freestream on all meshes
    if (myid == mblk2nd(1)) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
    }

    // Get some settings for RSM, and print some info:
    if (ivmx == 72) {
        module_kwstm_ns::kws_init(nou, bou, nbuf, ibufdim);
    }

    // do 1100 igrid=1,ngrid
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nblg(igrid) - 1;
        ncg = ncgg(igrid);
        // do 1110 n=1,ncg+1
        for (n = 1; n <= ncg+1; n++) {
            iem = iemg(igrid);
            nbl = nbl + 1;
            if (myid == mblk2nd(nbl)) {
                lead_ns::lead(nbl, lw, lw2, maxbl);
                isksav = isklton;
                isklton = 0;
                idum = 0;
                nttuse = std::max(ntt-1, 1);
                // first call to init for safety
                init_ns::init(nbl, jdim, kdim, idim,
                              FortranArray4DRef<double>(&w(lq), jdim, kdim, idim, 5), FortranArray4DRef<double>(&w(lqj0), kdim, idim-1, 5, 4), FortranArray4DRef<double>(&w(lqk0), jdim, idim-1, 5, 4), FortranArray4DRef<double>(&w(lqi0), jdim, kdim, 5, 4),
                              FortranArray4DRef<double>(&w(ltj0), kdim, idim-1, nummem, 4), FortranArray4DRef<double>(&w(ltk0), jdim, idim-1, nummem, 4), FortranArray4DRef<double>(&w(lti0), jdim, kdim, nummem, 4),
                              FortranArray3DRef<double>(&w(lvol), jdim, kdim, idim-1), FortranArray3DRef<double>(&w(lvolj0), kdim, idim-1, 4), FortranArray3DRef<double>(&w(lvolk0), jdim, idim-1, 4), FortranArray3DRef<double>(&w(lvoli0), jdim, kdim, 4),
                              nummem, FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                              nou, bou, nbuf, ibufdim, idum);
                // call the noninertial initialization to speed convergence
                if (noninflag > 0) {
                    initnonin_ns::initnonin(nbl, jdim, kdim, idim,
                                            FortranArray4DRef<double>(&w(lq), jdim, kdim, idim, 5), FortranArray4DRef<double>(&w(lqj0), kdim, idim-1, 5, 4), FortranArray4DRef<double>(&w(lqk0), jdim, idim-1, 5, 4), FortranArray4DRef<double>(&w(lqi0), jdim, kdim, 5, 4),
                                            FortranArray3DRef<double>(&w(lvol), jdim, kdim, idim-1), FortranArray3DRef<double>(&w(lvolj0), kdim, idim-1, 4), FortranArray3DRef<double>(&w(lvolk0), jdim, idim-1, 4), FortranArray3DRef<double>(&w(lvoli0), jdim, kdim, 4),
                                            FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim));
                }
                bc_ns::bc(idum, nbl, lw, lw2, w, mgwk, wk, nwork, clw(nttuse),
                          nou, bou, nbuf, ibufdim, maxbl, maxgr, maxseg, itrans,
                          irotat, idefrm, igridg, nblg, nbci0, nbcj0, nbck0,
                          nbcidim, nbcjdim, nbckdim, ibcinfo, jbcinfo,
                          kbcinfo, bcfilei, bcfilej, bcfilek, lwdat, myid,
                          idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
                isklton = isksav;
                if (n == 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " setting initial conditions as freestream for block%6d  (grid%6d)",
                        nbl, igrid);
                    if (noninflag > 0) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "  correcting initial conditions for NONINERTIAL reference frame");
                    }
                } else {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "   setting initial conditions as freestream for coarser block%6d", nbl);
                    if (noninflag > 0) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "    correcting initial conditions for NONINERTIAL reference frame");
                    }
                }
                {
                    int iflag2 = 1;
                    init_ns::init(nbl, jdim, kdim, idim,
                                  FortranArray4DRef<double>(&w(lq), jdim, kdim, idim, 5), FortranArray4DRef<double>(&w(lqj0), kdim, idim-1, 5, 4), FortranArray4DRef<double>(&w(lqk0), jdim, idim-1, 5, 4), FortranArray4DRef<double>(&w(lqi0), jdim, kdim, 5, 4),
                                  FortranArray4DRef<double>(&w(ltj0), kdim, idim-1, nummem, 4), FortranArray4DRef<double>(&w(ltk0), jdim, idim-1, nummem, 4), FortranArray4DRef<double>(&w(lti0), jdim, kdim, nummem, 4),
                                  FortranArray3DRef<double>(&w(lvol), jdim, kdim, idim-1), FortranArray3DRef<double>(&w(lvolj0), kdim, idim-1, 4), FortranArray3DRef<double>(&w(lvolk0), jdim, idim-1, 4), FortranArray3DRef<double>(&w(lvoli0), jdim, kdim, 4),
                                  nummem, FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                                  nou, bou, nbuf, ibufdim, iflag2);
                }
                if (noninflag > 0) {
                    initnonin_ns::initnonin(nbl, jdim, kdim, idim,
                                            FortranArray4DRef<double>(&w(lq), jdim, kdim, idim, 5), FortranArray4DRef<double>(&w(lqj0), kdim, idim-1, 5, 4), FortranArray4DRef<double>(&w(lqk0), jdim, idim-1, 5, 4), FortranArray4DRef<double>(&w(lqi0), jdim, kdim, 5, 4),
                                            FortranArray3DRef<double>(&w(lvol), jdim, kdim, idim-1), FortranArray3DRef<double>(&w(lvolj0), kdim, idim-1, 4), FortranArray3DRef<double>(&w(lvolk0), jdim, idim-1, 4), FortranArray3DRef<double>(&w(lvoli0), jdim, kdim, 4),
                                            FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim));
                }
            }
        } // end do 1110

        nblout = nblg(igrid);
        { int iunit11 = 11; writ_buf_ns::writ_buf(nblout, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
    } // end do 1100

    // initialize rigid body dynamics arrays/data
    if (irigb > 0 || irbtrim > 0) {
        init_rb_ns::init_rb(maxaes, maxbl, zorig, aesrfdat);
    }
    init_trim_ns::init_trim();


    // initialize aeroelastic arrays/data and read modal surface definitions
    if (naesrf > 0) {
        if (myid == myhost) {
            init_ae_ns::init_ae(nmds, maxaes, x0, xxn, FortranArray2DRef<double>(&wk(1), 2*nmds, 2*nmds), bmat, stm, stmi,
                                gforcn, gforcnm, freq, damp, gmass, aesrfdat, xs, gforcs);
        }

        // read in modal surface definitions
        for (iaes = 1; iaes <= naesrf; iaes++) {
            nmodes = (int)aesrfdat(5, iaes);

            // setup for plot3d output of modal surface definitions
            if (myid == myhost) {
                idum    = 1;
                jdum    = 1;
                kdum    = 1;
                mxaedum = 1;
                nmddum  = 1;
                nbldum  = 1;
                lwk1    = 1;
                lwk2    = lwk1 + idum*jdum*kdum;
                lwk3    = lwk2 + idum*jdum*kdum;
                lwk4    = lwk3 + kdum*idum*6*nmddum*mxaedum;
                lwk5    = lwk4 + kdum*idum*6*nmddum*mxaedum;
                lwk6    = lwk5 + kdum*idum*6*nmddum*mxaedum;
                for (nm = 1; nm <= nmodes; nm++) {
                    iunitw = 200 + nm;
                    { int iflag_pm = 0;
                    pltmode_ns::pltmode(nm, iaes, nblg, ngrid, maxgr, maxbl,
                                        nsegdfrm, iaesurf, jbcinfo, kbcinfo,
                                        ibcinfo, nbcj0, nbcjdim, nbck0, nbckdim,
                                        nbci0, nbcidim, maxseg, maxsegdg, lw, lw2,
                                        FortranArray3DRef<double>(&wk(lwk1), jdum, kdum, idum), FortranArray3DRef<double>(&wk(lwk2), jdum, kdum, idum), FortranArray3DRef<double>(&wk(lwk3), jdum, kdum, idum),
                                        FortranArray5DRef<double>(&wk(lwk4), kdum, idum, 6, nmddum, mxaedum), FortranArray5DRef<double>(&wk(lwk5), jdum, idum, 6, nmddum, mxaedum), FortranArray5DRef<double>(&wk(lwk6), jdum, kdum, 6, nmddum, mxaedum),
                                        jdum, kdum, idum, mxaedum, nmddum, nbldum, iunitw, iflag_pm); }
                }
            }

            for (igrid = 1; igrid <= ngrid; igrid++) {
                nbl = nblg(igrid);
                ncg = ncgg(igrid);
                lead_ns::lead(nbl, lw, lw2, maxbl);
                // initialize modal data
                for (lll = lxmdj; lll <= lvelg-1; lll++) {
                    w(lll) = 0.0;
                }
            }
            for (igrid = 1; igrid <= ngrid; igrid++) {
                nbl = nblg(igrid);
                ncg = ncgg(igrid);
                lead_ns::lead(nbl, lw, lw2, maxbl);
                iaesrf = 0;
                for (is = 1; is <= nsegdfrm(nbl); is++) {
                    iaesrf = iaesrf + iaesurf(nbl, is);
                }
                if (iaesrf != 0 && myid == myhost) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " ");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " reading modal surface for block%6d of dimensions (I/J/K) :%6d%6d%6d",
                        nbl, idim, jdim, kdim);
                    iunit = 33;
                    jst = lxmdj;
                    kst = lxmdk;
                    ist = lxmdi;
                    for (nm = 1; nm <= nmodes; nm++) {
                        modread_ns::modread(idim, jdim, kdim, nm, nbl, iunit,
                                            jbcinfo, kbcinfo, ibcinfo,
                                            nbcj0, nbcjdim, nbck0, nbckdim,
                                            nbci0, nbcidim, maxbl, maxseg,
                                            nmds, FortranArray5DRef<double>(&w(jst), kdim, idim, 6, nmds, maxaes), FortranArray5DRef<double>(&w(kst), jdim, idim, 6, nmds, maxaes), FortranArray5DRef<double>(&w(ist), jdim, kdim, 6, nmds, maxaes),
                                            iaes, maxaes);
                        iunitw = 200 + nm;
                        { int iflag_pm = 1;
                        pltmode_ns::pltmode(nm, iaes, nblg, ngrid, maxgr, maxbl,
                                            nsegdfrm, iaesurf, jbcinfo, kbcinfo,
                                            ibcinfo, nbcj0, nbcjdim, nbck0, nbckdim,
                                            nbci0, nbcidim, maxseg, maxsegdg, lw, lw2,
                                            FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                                            FortranArray5DRef<double>(&w(jst), kdim, idim, 6, nmds, maxaes), FortranArray5DRef<double>(&w(kst), jdim, idim, 6, nmds, maxaes), FortranArray5DRef<double>(&w(ist), jdim, kdim, 6, nmds, maxaes),
                                            jdim, kdim, idim, maxaes, nmds, nbl, iunitw, iflag_pm); }
                    }
                }

                // create modal surface definitions for coarser blocks
                if (ncg > 0) {
                    for (m = 1; m <= ncg; m++) {
                        nbl = nbl + 1;
                        if (iaesrf != 0 && myid == mblk2nd(nbl)) {
                            lxmdjc = lw(52, nbl);
                            lxmdkc = lw(53, nbl);
                            lxmdic = lw(54, nbl);
                            nou(1) = std::min(nou(1)+1, ibufdim);
                            std::snprintf(bou(nou(1),1), 120,
                                "   creating modal surface for coarser block%6d", nbl);
                            jst = lxmdj;
                            kst = lxmdk;
                            ist = lxmdi;
                            collmod_ns::collmod(FortranArray5DRef<double>(&w(lxmdj), kdim, idim, 6, nmds, maxaes), FortranArray5DRef<double>(&w(lxmdk), jdim, idim, 6, nmds, maxaes), FortranArray5DRef<double>(&w(lxmdi), jdim, kdim, 6, nmds, maxaes),
                                                FortranArray5DRef<double>(&w(lxmdjc), kk2, ii2, 6, nmds, maxaes), FortranArray5DRef<double>(&w(lxmdkc), jj2, ii2, 6, nmds, maxaes), FortranArray5DRef<double>(&w(lxmdic), jj2, kk2, 6, nmds, maxaes),
                                                jdim, kdim, idim, jj2, kk2, ii2,
                                                nm, nmds, iaes, maxaes);
                            lead_ns::lead(nbl, lw, lw2, maxbl);
                        }
                    }
                }

                nblout = nblg(igrid);
                { int iunit11 = 11; writ_buf_ns::writ_buf(nblout, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
            } // end igrid loop

            for (nm = 1; nm <= nmodes; nm++) {
                fortran_close_unit(200 + nm);
            }
        } // end iaes loop
    } // end naesrf > 0


    // read restart data on mseq-1 level from finest level
    if (irest > 0) {
        if (myid == mblk2nd(1)) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " ");
        }

        if (iteravg == 2) {
            if (myid == myhost) {
                FILE* f97 = fortran_get_unit(97);
                // read(97,end=1011,err=1011) ngx
                int ngx = 0;
                bool ok = true;
                if (std::fread(&ngx, sizeof(int), 1, f97) != 1) ok = false;
                if (ok && ngx != ngrid) ok = false;
                if (ok) {
                    for (igrid = 1; igrid <= ngx; igrid++) {
                        if (std::fread(&itest(igrid), sizeof(int), 1, f97) != 1) { ok = false; break; }
                        if (std::fread(&jtest(igrid), sizeof(int), 1, f97) != 1) { ok = false; break; }
                        if (std::fread(&ktest(igrid), sizeof(int), 1, f97) != 1) { ok = false; break; }
                    }
                }
                if (ok && ipertavg != 0) {
                    for (igrid = 1; igrid <= ngrid; igrid++) {
                        if (itest(igrid) != idimg(nblg(igrid)+mseq-1) ||
                            jtest(igrid) != jdimg(nblg(igrid)+mseq-1) ||
                            ktest(igrid) != kdimg(nblg(igrid)+mseq-1)) { ok = false; break; }
                    }
                    if (ok) {
                        FILE* f98 = fortran_get_unit(98);
                        int ngx2 = 0;
                        if (std::fread(&ngx2, sizeof(int), 1, f98) != 1) ok = false;
                        if (ok && ngx2 != ngrid) ok = false;
                        if (ok) {
                            for (igrid = 1; igrid <= ngx2; igrid++) {
                                if (std::fread(&itest(igrid), sizeof(int), 1, f98) != 1) { ok = false; break; }
                                if (std::fread(&jtest(igrid), sizeof(int), 1, f98) != 1) { ok = false; break; }
                                if (std::fread(&ktest(igrid), sizeof(int), 1, f98) != 1) { ok = false; break; }
                            }
                        }
                        if (ok) {
                            for (igrid = 1; igrid <= ngrid; igrid++) {
                                if (itest(igrid) != idimg(nblg(igrid)+mseq-1) ||
                                    jtest(igrid) != jdimg(nblg(igrid)+mseq-1) ||
                                    ktest(igrid) != kdimg(nblg(igrid)+mseq-1)) { ok = false; break; }
                            }
                        }
                    }
                } else if (ok) {
                    for (igrid = 1; igrid <= ngrid; igrid++) {
                        if (itest(igrid) != idimg(nblg(igrid)+mseq-1)-1 ||
                            jtest(igrid) != jdimg(nblg(igrid)+mseq-1)-1 ||
                            ktest(igrid) != kdimg(nblg(igrid)+mseq-1)-1) { ok = false; break; }
                    }
                }
                if (!ok) {
                    FILE* f11 = fortran_get_unit(11);
                    std::fprintf(f11, "\n stopping... flag iteravg=2, but cannot read running-average Q file\n");
                    std::fprintf(f11, " ... it either does not exist or it is an incorrect file\n");
                    termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                }
            }
        }

        // do 1700 igrid=1,ngrid
        for (igrid = 1; igrid <= ngrid; igrid++) {
            iskipz = 0;
            if (igrid == 1) iskipz = 1;
            iem   = iemg(igrid);
            nbl   = nblg(igrid);
            inewg = inewgg(igrid);
            if (inewg > 0) goto label_1701;
            if (iem >= 1 && mseq > 1) goto label_1701;
            if (icgns == 1) {
                lead_ns::lead(nbl, lw, lw2, maxbl);
                idima = idim;
                jdima = jdim;
                kdima = kdim;
            }
            nbl = nbl + (mseq - 1);
            lead_ns::lead(nbl, lw, lw2, maxbl);
            if (myid == mblk2nd(nbl)) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " reading restart file for block%6d  (grid%6d)", nbl, igrid);
                if (iteravg == 1) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " running-average to be started from scratch for block%6d", nbl);
                }
                if (iteravg == 2) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " reading running-average file for block%6d", nbl);
                }
                if (icgns == 1) {
                    if (nwork < (idima+1)*(jdima+1)*(kdima+1)) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " not enough memory for cgns Q read.");
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " nwork in wk=%6d.  Needed = %6d",
                            nwork, (idima+1)*(jdima+1)*(kdima+1));
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                } else {
                    idima = idim;
                    jdima = jdim;
                    kdima = kdim;
                }
            }
            if (mblk2nd(nbl) == myid || myid == myhost) {
                rrest_ns::rrest(nbl, jdim, kdim, idim,
                                FortranArray4DRef<double>(&w(lq), jdim, kdim, idim, 5), FortranArray4DRef<double>(&w(lqj0), kdim, idim-1, 5, 4), FortranArray4DRef<double>(&w(lqk0), jdim, idim-1, 5, 4), FortranArray4DRef<double>(&w(lqi0), jdim, kdim, 5, 4),
                                ncycmax, ntr, rms, clw, cdw, cdpw, cdvw, cxw, cyw, czw,
                                cmxw, cmyw, cmzw,
                                n_clcd, clcd, nblocks_clcd, blocks_clcd,
                                fmdotw, cftmomw, cftpw, cftvw,
                                cfttotw, rmstr, nneg, iskipz,
                                FortranArray3DRef<double>(&w(lvis), jdim, kdim, idim), FortranArray4DRef<double>(&w(lxib), jdim, kdim, idim, nummem), FortranArray3DRef<double>(&w(lsnk0), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lsni0), jdim-1, kdim-1, idim-1),
                                FortranArray3DRef<double>(&w(lxkb), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lnbl), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lcmuv), jdim-1, kdim-1, idim-1),
                                maxbl, mblk2nd, myid, myhost, mycomm,
                                nou, bou, nbuf, ibufdim, igrid,
                                FortranArray3DRef<double>(&wk(1), idima+1, jdima+1, kdima+1), idima, jdima, kdima,
                                FortranArray4DRef<double>(&w(lvj0), kdim, idim-1, 1, 4), FortranArray4DRef<double>(&w(lvk0), jdim, idim-1, 1, 4), FortranArray4DRef<double>(&w(lvi0), jdim, kdim, 1, 4),
                                FortranArray4DRef<double>(&w(ltj0), kdim, idim-1, nummem, 4), FortranArray4DRef<double>(&w(ltk0), jdim, idim-1, nummem, 4), FortranArray4DRef<double>(&w(lti0), jdim, kdim, nummem, 4),
                                FortranArray4DRef<double>(&w(lqavg), jdim, kdim, idim, 5), FortranArray4DRef<double>(&w(lq2avg), jdim, kdim, idim, 5), nummem);
                ntq = ntt;
                ncycchk = ncyctot + ntt;
                if (myid == myhost) {
                    if (ncycchk > ncycmax) {
                        FILE* f11 = fortran_get_unit(11);
                        std::fprintf(f11, "\n");
                        std::fprintf(f11, " stopping...ncycmax too small to store prior convergence history AND current history \n");
                        std::fprintf(f11, " increase value of ncycmax to at LEAST %6d%6d%6d%6d\n",
                                     ncycchk, ncyctot, ncycchk, ntt);
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                }
            }
            label_1701:;
            { int iunit11 = 11; writ_buf_ns::writ_buf(nbl, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
        } // end do 1700


        iflagg = 0;
        if (icgns != 1) {
            FILE* f2 = fortran_get_unit(2);
            if (std::fread(&iflagg, sizeof(int), 1, f2) != 1) {
                iflagg = 0;
                goto label_8989;
            }
        }

        // set iflagg = 0 if steady state (dt < 0)
        if ((float)dt < 0.0f) iflagg = 0;
        if (iflagg == 0) goto label_8989;

        // read data for 2nd order time advancement
        if (iflagg == 1 || iflagg == 3) {
            // do 1810 igrid=1,ngrid
            for (igrid = 1; igrid <= ngrid; igrid++) {
                iem   = iemg(igrid);
                nbl   = nblg(igrid);
                inewg = inewgg(igrid);
                if (irest == 0 || inewg > 0) goto label_1811;
                if (iem >= 1 && mseq > 1) goto label_1811;
                if (icgns == 1) {
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    idima = idim;
                    jdima = jdim;
                    kdima = kdim;
                }
                nbl = nbl + (mseq - 1);
                if (mblk2nd(nbl) == myid || myid == myhost) {
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    if (icgns == 1) {
                        if (nwork < (idima+1)*(jdima+1)*(kdima+1)) {
                            FILE* f11 = fortran_get_unit(11);
                            std::fprintf(f11, " not enough memory for cgns unsteady Q read.\n");
                            std::fprintf(f11, " nwork in wk=%6d.  Needed = %6d\n",
                                         nwork, (idima+1)*(jdima+1)*(kdima+1));
                            termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                        }
                    } else {
                        idima = idim;
                        jdima = jdim;
                        kdima = kdim;
                    }
                    {
                        int igrdflag_local = 0;
                        rrestg_ns::rrestg(nbl, igrid, jdim, kdim, idim,
                                          FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                                          FortranArray3DRef<double>(&w(lxnm2), jdim, kdim, idim), FortranArray3DRef<double>(&w(lynm2), jdim, kdim, idim), FortranArray3DRef<double>(&w(lznm2), jdim, kdim, idim),
                                          FortranArray4DRef<double>(&w(ldeltj), kdim, idim, 3, 2), FortranArray4DRef<double>(&w(ldeltk), jdim, idim, 3, 2), FortranArray4DRef<double>(&w(ldelti), jdim, kdim, 3, 2),
                                          FortranArray4DRef<double>(&w(lqc0), jdim, kdim, idim-1, 5), igrdflag_local, iuns,
                                          utrans, vtrans, wtrans,
                                          omegax, omegay, omegaz, xorig, yorig, zorig,
                                          dxmx, dymx, dzmx, dthxmx, dthymx,
                                          dthzmx, thetax, thetay, thetaz, rfreqt,
                                          rfreqr, xorig0, yorig0, zorig0, time2,
                                          thetaxl, thetayl, thetazl, itrans, irotat, idefrm,
                                          utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae,
                                          xorgae, yorgae, zorgae, thtxae, thtyae, thtzae,
                                          rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf,
                                          kcsi, kcsf, freq, gmass, damp, x0, gf0, nmds, maxaes,
                                          aesrfdat, perturb, myhost, myid, mycomm, mblk2nd,
                                          maxbl, ibufdim, nbuf, bou, nou, nsegdfrm, idfrmseg,
                                          iaesurf, maxsegdg, wk, nwork, idima, jdima, kdima,
                                          FortranArray4DRef<double>(&w(lxib2), jdim, kdim, idim, 2*nummem), nummem);
                    }
                }
                label_1811:;
            } // end do 1810
        } else {
            if (myid == mblk2nd(1)) {
                if (time != 0.0f) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " WARNING:  You are restarting a 2nd order in time solution with no");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "           qc0 data saved in restart. You may see a glitch in the");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "           restart as a result.");
                }
            }
        }


        // read dynamic mesh data on restart level and reconstruct last grid position
        if ((float)dt > 0.0f && iunst > 0) {
            if (iflagg == 2 || iflagg == 3) {
                if (ndgrd != 0 && iunst > 1) {
                    if (myid == myhost) {
                        fortran_open_unit(98, "dgplot3d.bin", "r");
                        // read(98,*) ngridt
                        fortran_read_list(98, &ngridt);
                        // read(98,*) (idimt,jdimt,kdimt,igrid=1,ngrid) - single record
                        {
                            FortranArray1D<int> tmp_dims(3*ngrid);
                            fortran_read_list_array(98, tmp_dims, 1, 3*ngrid);
                            for (int ig = 1; ig <= ngrid; ig++) {
                                idimt = tmp_dims(3*(ig-1)+1);
                                jdimt = tmp_dims(3*(ig-1)+2);
                                kdimt = tmp_dims(3*(ig-1)+3);
                            }
                        }
                    } // end if (myid == myhost)
                } // end if (ndgrd != 0 && iunst > 1)
                // do 1800 igrid=1,ngrid
                for (igrid = 1; igrid <= ngrid; igrid++) {
                    iem   = iemg(igrid);
                    nbl   = nblg(igrid);
                    inewg = inewgg(igrid);
                    if (irest == 0 || inewg > 0) goto label_1801;
                    if (iem >= 1 && mseq > 1) goto label_1801;
                    if (icgns == 1) {
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        idima = idim;
                        jdima = jdim;
                        kdima = kdim;
                    }
                    nbl = nbl + (mseq - 1);
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    if (icgns == 1) {
                        if (nwork < (idima+1)*(jdima+1)*(kdima+1)) {
                            FILE* f901 = fortran_get_unit(901);
                            std::fprintf(f901, " not enough memory for cgns unsteady Q read.\n");
                            std::fprintf(f901, " nwork in wk=%6d.  Needed = %6d\n",
                                         nwork, (idima+1)*(jdima+1)*(kdima+1));
                            termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                        }
                    } else {
                        idima = idim;
                        jdima = jdim;
                        kdima = kdim;
                    }
                    if (mblk2nd(nbl) == myid || myid == myhost) {
                        iunsn = std::max({itrans(nbl), irotat(nbl), idefrm(nbl)});
                        {
                            int igrdflag_local = 1;
                            rrestg_ns::rrestg(nbl, igrid, jdim, kdim, idim,
                                              FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                                              FortranArray3DRef<double>(&w(lxnm2), jdim, kdim, idim), FortranArray3DRef<double>(&w(lynm2), jdim, kdim, idim), FortranArray3DRef<double>(&w(lznm2), jdim, kdim, idim),
                                              FortranArray4DRef<double>(&w(ldeltj), kdim, idim, 3, 2), FortranArray4DRef<double>(&w(ldeltk), jdim, idim, 3, 2), FortranArray4DRef<double>(&w(ldelti), jdim, kdim, 3, 2),
                                              FortranArray4DRef<double>(&w(lqc0), jdim, kdim, idim-1, 5), igrdflag_local, iuns,
                                              utrans, vtrans, wtrans,
                                              omegax, omegay, omegaz, xorig, yorig, zorig,
                                              dxmx, dymx, dzmx, dthxmx, dthymx,
                                              dthzmx, thetax, thetay, thetaz, rfreqt,
                                              rfreqr, xorig0, yorig0, zorig0, time2,
                                              thetaxl, thetayl, thetazl, itrans, irotat, idefrm,
                                              utrnsae, vtrnsae, wtrnsae, omgxae, omgyae, omgzae,
                                              xorgae, yorgae, zorgae, thtxae, thtyae, thtzae,
                                              rfrqtae, rfrqrae, icsi, icsf, jcsi, jcsf,
                                              kcsi, kcsf, freq, gmass, damp, x0, gf0, nmds, maxaes,
                                              aesrfdat, perturb, myhost, myid, mycomm, mblk2nd,
                                              maxbl, ibufdim, nbuf, bou, nou, nsegdfrm, idfrmseg,
                                              iaesurf, maxsegdg, wk, nwork, idima, jdima, kdima,
                                              FortranArray4DRef<double>(&w(lxib2), jdim, kdim, idim, 2*nummem), nummem);
                        }
                        if (myid == myhost) {
                            if (iuns != iunsn) {
                                FILE* f11 = fortran_get_unit(11);
                                std::fprintf(f11, " Stopping: cannot alter type of grid motion between restarts\n");
                                termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                            }
                        }
                    }
                    label_1801:;
                    { int iunit11 = 11; writ_buf_ns::writ_buf(nbl, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
                } // end do 1800

                // reconstruct last grid position for purely rigid grid case
                if (iunst > 0) {
                    for (igrid = 1; igrid <= ngrid; igrid++) {
                        nbl = nblg(igrid);
                        nbl = nbl + (mseq - 1);
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        iuns_local = std::max(itrans(nbl), irotat(nbl));
                        if (mblk2nd(nbl) == myid) {
                            if (idefrm(nbl) == 0 && iuns_local > 0) {
                                nou(1) = std::min(nou(1)+1, ibufdim);
                                std::snprintf(bou(nou(1),1), 120,
                                    " reconstructing last grid position from data for block%6d", nbl);
                                { double gm_xc1 = xorig0(nbl), gm_yc1 = yorig0(nbl), gm_zc1 = zorig0(nbl);
                                  double gm_xc2 = xorig(nbl), gm_yc2 = yorig(nbl), gm_zc2 = zorig(nbl);
                                  double gm_rx = thetax(nbl), gm_ry = thetay(nbl), gm_rz = thetaz(nbl);
                                grdmove_ns::grdmove(nbl, jdim, kdim, idim,
                                    FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                                    gm_xc1, gm_yc1, gm_zc1,
                                    gm_xc2, gm_yc2, gm_zc2,
                                    gm_rx, gm_ry, gm_rz); }
                            }
                        }
                        { int iunit11 = 11; writ_buf_ns::writ_buf(nbl, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
                    }
                }
            }
        }

        goto label_8990;

        label_8989:;
        if (myid == mblk2nd(1)) {
            if ((float)dt > 0.0f && std::abs(ita) == 2) {
                if (time != 0.0f) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " WARNING:  You are restarting a 2nd order in time solution with no");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "           qc0 data saved in restart. You may see a glitch in the");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "           restart as a result.");
                }
            }
        }

        label_8990:;


        // read modal data if aeroelastic
        if (myid == myhost) {
            if (naesrf > 0) {
                if (icgns != 1) {
                    FILE* f2 = fortran_get_unit(2);
                    // read(2,end=8995) (timekeep(nn),nn=1,ntt)
                    bool ok8995 = true;
                    for (nn = 1; nn <= ntt && ok8995; nn++) {
                        double tmp;
                        if (std::fread(&tmp, sizeof(double), 1, f2) != 1) { ok8995 = false; break; }
                        timekeep(nn) = tmp;
                    }
                    if (ok8995) {
                        for (iaes = 1; iaes <= naesrf && ok8995; iaes++) {
                            nmodes = (int)aesrfdat(5, iaes);
                            // read(2,end=8995) (xxn(n,iaes),n=1,2*nmodes)
                            for (n = 1; n <= 2*nmodes && ok8995; n++) {
                                double tmp;
                                if (std::fread(&tmp, sizeof(double), 1, f2) != 1) { ok8995 = false; break; }
                                xxn(n, iaes) = tmp;
                            }
                            // read(2,end=8995) (gforcn(n,iaes),n=1,2*nmodes)
                            for (n = 1; n <= 2*nmodes && ok8995; n++) {
                                double tmp;
                                if (std::fread(&tmp, sizeof(double), 1, f2) != 1) { ok8995 = false; break; }
                                gforcn(n, iaes) = tmp;
                            }
                            // read(2,end=8995) (gforcnm(n,iaes),n=1,2*nmodes)
                            for (n = 1; n <= 2*nmodes && ok8995; n++) {
                                double tmp;
                                if (std::fread(&tmp, sizeof(double), 1, f2) != 1) { ok8995 = false; break; }
                                gforcnm(n, iaes) = tmp;
                            }
                            // read(2,end=8995) (((aehist(nn,ll,n,iaes),nn=1,ntt),ll=1,3),n=1,nmodes)
                            for (n = 1; n <= nmodes && ok8995; n++) {
                                for (ll = 1; ll <= 3 && ok8995; ll++) {
                                    for (nn = 1; nn <= ntt && ok8995; nn++) {
                                        double tmp;
                                        if (std::fread(&tmp, sizeof(double), 1, f2) != 1) { ok8995 = false; break; }
                                        aehist(nn, ll, n, iaes) = tmp;
                                    }
                                }
                            }
                            nou(1) = std::min(nou(1)+1, ibufdim);
                            std::snprintf(bou(nou(1),1), 120,
                                " reading generalized force/displacement data for aeroelastic surface %6d",
                                iaes);
                            // if x0 data is non-zero in input file, overwrite xxn with x0
                            for (n = 1; n <= 2*nmodes; n++) {
                                if (std::fabs((float)x0(n, iaes)) != 0.0f) {
                                    xxn(n, iaes) = x0(n, iaes);
                                }
                            }
                        }
                    }
                    // label_8995:
                }
            }
        }
    } // end irest > 0

    if (myid == myhost) {
        if (icgns != 1) {
            fortran_close_unit(2);
            fortran_close_unit(102);
            if (ndgrd != 0 && iunst > 1) {
                fortran_close_unit(98);
            }
        }
    }


    // compute directed distance function for Baldwin-Lomax turbulence model
    if (myid == mblk2nd(1)) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
    }

    // do 1900 igrid=1,ngrid
    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nblg(igrid);
        lead_ns::lead(nbl, lw, lw2, maxbl);
        if (myid == mblk2nd(nbl)) {
            if (ivisc[0] == 2 || ivisc[1] == 2 || ivisc[2] == 2 ||
                ivisc[0] == 3 || ivisc[1] == 3 || ivisc[2] == 3) {
                inewg = inewgg(igrid);
                int numchk = (jdim-1)*(kdim-1)*(idim-1)*4;
                if (nwork < numchk) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " not enough work space for subroutine dird (must increase mwork)");
                    termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                }
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " computing directed distances for block%6d  (grid%6d)", nbl, igrid);
                {
                    FortranArray1DRef<int> ivisc_ref(ivisc, 1, 3);
                    dird_ns::dird(jdim, kdim, idim, FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                         FortranArray4DRef<double>(&w(lsj), jdim, kdim, idim-1, 5), FortranArray4DRef<double>(&w(lsk), jdim, kdim, idim-1, 5), FortranArray4DRef<double>(&w(lsi), jdim, kdim, idim, 5),
                         FortranArray3DRef<double>(&w(lsni0), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lsnk0), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lsni0), jdim-1, kdim-1, idim-1),
                         FortranArray3DRef<double>(&w(lnbl), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lxkb), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lnbl), jdim-1, kdim-1, idim-1),
                         FortranArray4DRef<double>(&wk(1), jdim-1, kdim-1, idim-1, 4), ivisc_ref, nou, bou, nbuf, ibufdim);
                }

                // compute directed distances on coarser levels
                ncg = ncgg(igrid);
                iem = iemg(igrid);
                if (ncg > 0) {
                    for (m = 1; m <= ncg; m++) {
                        nbl = nbl + 1;
                        lead_ns::lead(nbl, lw, lw2, maxbl);
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "   computing directed distances for coarser block%6d", nbl);
                        FortranArray1DRef<int> ivisc_ref(ivisc, 1, 3);
                        dird_ns::dird(jdim, kdim, idim, FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim),
                             FortranArray4DRef<double>(&w(lsj), jdim, kdim, idim-1, 5), FortranArray4DRef<double>(&w(lsk), jdim, kdim, idim-1, 5), FortranArray4DRef<double>(&w(lsi), jdim, kdim, idim, 5),
                             FortranArray3DRef<double>(&w(lsni0), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lsnk0), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lsni0), jdim-1, kdim-1, idim-1),
                             FortranArray3DRef<double>(&w(lnbl), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lxkb), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lnbl), jdim-1, kdim-1, idim-1),
                             FortranArray4DRef<double>(&wk(1), jdim-1, kdim-1, idim-1, 4), ivisc_ref, nou, bou, nbuf, ibufdim);
                    }
                }
            }
        }

        nblout = nblg(igrid);
        { int iunit11 = 11; writ_buf_ns::writ_buf(nblout, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
    } // end do 1900


    // minimum distance function for field equation turbulence models
    if (isminc >= 1) {
        if (isminc == 1 || isminc == 2) {
            // avoid computing min dist function on grid levels where it's not needed
            if (mseq == 2 && ncyc1[1] == 0) {
                j1 = 1;
            } else if (mseq == 3 && ncyc1[1] == 0) {
                j1 = 2;
            } else if (mseq == 3 && ncyc1[2] == 0) {
                j1 = 1;
            } else if (mseq == 4 && ncyc1[1] == 0) {
                j1 = 3;
            } else if (mseq == 4 && ncyc1[2] == 0) {
                j1 = 2;
            } else if (mseq == 4 && ncyc1[3] == 0) {
                j1 = 1;
            } else if (mseq == 5 && ncyc1[1] == 0) {
                j1 = 4;
            } else if (mseq == 5 && ncyc1[2] == 0) {
                j1 = 3;
            } else if (mseq == 5 && ncyc1[3] == 0) {
                j1 = 2;
            } else if (mseq == 5 && ncyc1[4] == 0) {
                j1 = 1;
            } else {
                j1 = 0;
            }

            // for recursive box algorithm, need some bc info on finest level
            j1sav = j1;
            j1 = 0;

            // Initialize smin to large value over all blocks
            itry = 0;
            for (n = 1; n <= ngrid; n++) {
                iset1 = nblg(n) + j1;
                iset2 = nblg(n) + ncgg(n);
                for (nbl = iset1; nbl <= iset2; nbl++) {
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    if (ivisc[0] >= 4 || ivisc[1] >= 4 || ivisc[2] >= 4) {
                        itry = itry + 1;
                        if (myid == mblk2nd(nbl)) {
                            ntot = (jdim-1)*(kdim-1);
                            izz = 0;
                            for (int i2 = 1; i2 <= idim-1; i2++) {
                                for (int m2 = 1; m2 <= ntot; m2++) {
                                    izz = izz + 1;
                                    w(lsnk0 + izz - 1) = 1.0e20;
                                }
                            }
                        }
                    }
                }
            }

            // If no walls at all, don't do min distance calculation
            for (n = 1; n <= ngrid; n++) {
                nbl = nblg(n);
                bool found_wall = false;
                for (iset = 1; iset <= nbci0(nbl) && !found_wall; iset++) {
                    int bc1 = std::abs(ibcinfo(nbl,iset,1,1));
                    if (bc1==2004||bc1==2014||bc1==2024||bc1==2034||bc1==2016) found_wall=true;
                }
                for (iset = 1; iset <= nbcidim(nbl) && !found_wall; iset++) {
                    int bc1 = std::abs(ibcinfo(nbl,iset,1,2));
                    if (bc1==2004||bc1==2014||bc1==2024||bc1==2034||bc1==2016) found_wall=true;
                }
                for (iset = 1; iset <= nbcj0(nbl) && !found_wall; iset++) {
                    int bc1 = std::abs(jbcinfo(nbl,iset,1,1));
                    if (bc1==2004||bc1==2014||bc1==2024||bc1==2034||bc1==2016) found_wall=true;
                }
                for (iset = 1; iset <= nbcjdim(nbl) && !found_wall; iset++) {
                    int bc1 = std::abs(jbcinfo(nbl,iset,1,2));
                    if (bc1==2004||bc1==2014||bc1==2024||bc1==2034||bc1==2016) found_wall=true;
                }
                for (iset = 1; iset <= nbck0(nbl) && !found_wall; iset++) {
                    int bc1 = std::abs(kbcinfo(nbl,iset,1,1));
                    if (bc1==2004||bc1==2014||bc1==2024||bc1==2034||bc1==2016) found_wall=true;
                }
                for (iset = 1; iset <= nbckdim(nbl) && !found_wall; iset++) {
                    int bc1 = std::abs(kbcinfo(nbl,iset,1,2));
                    if (bc1==2004||bc1==2014||bc1==2024||bc1==2034||bc1==2016) found_wall=true;
                }
                if (found_wall) goto label_5860;
            }
            itry = 0;
            label_5860:;


            if (itry != 0) {
                // calculate minimum distance via recursive-box algorithm
                j1 = j1sav;

                // check storage availability
                nworkf = 0;
                nworki = 0;
                if (myid == myhost) {
                    { int bctype = 2004; cntsurf_ns::cntsurf(ns2004, maxbl, maxgr, maxseg, ngrid, nblg,
                                        nbci0, nbcj0, nbck0, nbcidim, nbcjdim,
                                        nbckdim, ibcinfo, jbcinfo, kbcinfo, bctype); }
                    { int bctype = 2014; cntsurf_ns::cntsurf(ns2014, maxbl, maxgr, maxseg, ngrid, nblg,
                                        nbci0, nbcj0, nbck0, nbcidim, nbcjdim,
                                        nbckdim, ibcinfo, jbcinfo, kbcinfo, bctype); }
                    { int bctype = 2024; cntsurf_ns::cntsurf(ns2024, maxbl, maxgr, maxseg, ngrid, nblg,
                                        nbci0, nbcj0, nbck0, nbcidim, nbcjdim,
                                        nbckdim, ibcinfo, jbcinfo, kbcinfo, bctype); }
                    { int bctype = 2034; cntsurf_ns::cntsurf(ns2034, maxbl, maxgr, maxseg, ngrid, nblg,
                                        nbci0, nbcj0, nbck0, nbcidim, nbcjdim,
                                        nbckdim, ibcinfo, jbcinfo, kbcinfo, bctype); }
                    { int bctype = 2016; cntsurf_ns::cntsurf(ns2016, maxbl, maxgr, maxseg, ngrid, nblg,
                                        nbci0, nbcj0, nbck0, nbcidim, nbcjdim,
                                        nbckdim, ibcinfo, jbcinfo, kbcinfo, bctype); }
                    nsurf = ns2004 + ns2014 + ns2024 + ns2034 + ns2016;
                }
                minbox = (int)std::sqrt((float)nsurf);
                minbox = std::max(minbox, 50);
                nbb = 3*nsurf/minbox;
                for (igrid = 1; igrid <= ngrid; igrid++) {
                    nbl = nblg(igrid);
                    ntempf = 9*nsurf + 7*nbb
                           + 4*jdimg(nbl)*kdimg(nbl)*idimg(nbl);
                    if (ivmx == 4 || ivmx == 25) {
                        ntempi = 15*nsurf + 2*nbb
                               + jdimg(nbl)*kdimg(nbl)*idimg(nbl);
                    } else {
                        ntempi = 11*nsurf + 2*nbb
                               + jdimg(nbl)*kdimg(nbl)*idimg(nbl);
                    }
                    nworkf = std::max(nworkf, ntempf);
                    nworki = std::max(nworki, ntempi);
                }
                // additional memory for temporary storage of smin
                nadd = 0;
                for (nbl = 1; nbl <= nblock; nbl++) {
                    if (mblk2nd(nbl) == myid) {
                        nadd = nadd + jdimg(nbl)*kdimg(nbl)*idimg(nbl);
                    }
                }
                nworkf = nworkf + nadd;

                nworkxs  = 4*nsurf;
                if (ivmx == 4 || ivmx == 25) {
                    nworkixs = 4*nsurf;
                } else {
                    nworkixs = 0;
                }

                nwklsmin = maxbl;
                nwkireq  = maxbl*maxseg*6;

                nroomf = nwork  - nworkf  - nworkxs;
                nroomi = iwork  - nworki  - nworkixs
                       - nwklsmin - 2*nwkireq;
                if (nroomf < 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " not enough work space for subroutine findmin_new");
                    termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                }
                if (nroomi < 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " not enough integer work space for subroutine findmin_new");
                    termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                }

                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " ");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " setting up for minimum distance calculation using recursive-box algorithm");
                lwk1 = 1;
                lwk2 = lwk1 + nworkf;
                iwk1 = 1;
                iwk2 = iwk1 + nworki;
                iwk3 = iwk2 + nworkixs;
                iwk4 = iwk3 + nwklsmin;
                iwk5 = iwk4 + nwkireq;
                iwk6 = iwk5 + nwkireq;
                for (int iii = 1; iii <= iwk6; iii++) {
                    iwk(iii) = 0;
                }
                findmin_new_ns::findmin_new(lw, lw2, w, mgwk, FortranArray1DRef<double>(&wk(lwk1), wk.size() - (lwk1) + 1), nworkf,
                            FortranArray1DRef<int>(&iwk(iwk1), iwk.size() - (iwk1) + 1), nworki, nsurf, j1,
                            FortranArray1DRef<double>(&wk(lwk2), wk.size() - (lwk2) + 1), FortranArray1DRef<int>(&iwk(iwk2), iwk.size() - (iwk2) + 1), FortranArray1DRef<int>(&iwk(iwk3), iwk.size() - (iwk3) + 1),
                            FortranArray1DRef<int>(&iwk(iwk4), iwk.size() - (iwk4) + 1), FortranArray1DRef<int>(&iwk(iwk5), iwk.size() - (iwk5) + 1), ngrid, ncgg,
                            nbci0, nbcj0, nbck0, nbcidim, nbcjdim,
                            nbckdim, jbcinfo, kbcinfo, ibcinfo,
                            nblg, nou, bou, nbuf, ibufdim, maxbl,
                            maxgr, maxseg, mblk2nd);

                for (nbl = 1; nbl <= nblock; nbl++) {
                    { int iunit11 = 11; writ_buf_ns::writ_buf(nbl, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
                }
            } // end itry != 0
        } // end isminc == 1 || isminc == 2


        // need to initialize turb quantities if there is turbulence on any block
        if (isminc == 1 || isminc == 3) {
            for (nbl = 1; nbl <= nblock; nbl++) {
                if (mblk2nd(nbl) == myid) {
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    if (nbl == 1) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " ");
                    }
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " setting turbulent initial conditions, block%6d", nbl);
                    initvist_ns::initvist(nbl, jdim, kdim, idim,
                                          FortranArray3DRef<double>(&w(lvis), jdim, kdim, idim), FortranArray4DRef<double>(&w(lxib), jdim, kdim, idim, nummem),
                                          FortranArray3DRef<double>(&w(lsnk0), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lcmuv), jdim-1, kdim-1, idim-1),
                                          nummem, FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim));
                }
                { int iunit11 = 11; writ_buf_ns::writ_buf(nbl, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
            }
        }
    } // end isminc >= 1

    // make smin negative in regions where want laminar flow (field equation models only)
    if (ivmx >= 4) {
        for (n = 1; n <= ngrid; n++) {
            llev = 0;
            for (nbl = nblg(n); nbl <= nblg(n)+ncgg(n); nbl++) {
                llev = llev + 1;
                if (mblk2nd(nbl) == myid) {
                    lead_ns::lead(nbl, lw, lw2, maxbl);
                    // J-dir:
                    for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
                        if (std::abs(jbcinfo(nbl,nseg,1,1)) == 2014) {
                            ibeg = jbcinfo(nbl,nseg,2,1);
                            iend = jbcinfo(nbl,nseg,3,1);
                            jbeg = 1;
                            jend = 1;
                            kbeg = jbcinfo(nbl,nseg,4,1);
                            kend = jbcinfo(nbl,nseg,5,1);
                            mdim = kend - kbeg;
                            ndim = iend - ibeg;
                            ldata = lwdat(nbl, nseg, 3);
                            { int iface_lf = 3;
                            lamfix_ns::lamfix(jdim, kdim, idim, FortranArray3DRef<double>(&w(lsnk0), jdim-1, kdim-1, idim-1),
                                              mdim, ndim, FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12),
                                              ibeg, iend, jbeg, jend, kbeg, kend, iface_lf, llev); }
                        }
                    }
                    for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
                        if (std::abs(jbcinfo(nbl,nseg,1,2)) == 2014) {
                            ibeg = jbcinfo(nbl,nseg,2,2);
                            iend = jbcinfo(nbl,nseg,3,2);
                            jbeg = jdim;
                            jend = jdim;
                            kbeg = jbcinfo(nbl,nseg,4,2);
                            kend = jbcinfo(nbl,nseg,5,2);
                            mdim = kend - kbeg;
                            ndim = iend - ibeg;
                            ldata = lwdat(nbl, nseg, 4);
                            { int iface_lf = 4;
                            lamfix_ns::lamfix(jdim, kdim, idim, FortranArray3DRef<double>(&w(lsnk0), jdim-1, kdim-1, idim-1),
                                              mdim, ndim, FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12),
                                              ibeg, iend, jbeg, jend, kbeg, kend, iface_lf, llev); }
                        }
                    }
                    // K-dir:
                    for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
                        if (std::abs(kbcinfo(nbl,nseg,1,1)) == 2014) {
                            ibeg = kbcinfo(nbl,nseg,2,1);
                            iend = kbcinfo(nbl,nseg,3,1);
                            jbeg = kbcinfo(nbl,nseg,4,1);
                            jend = kbcinfo(nbl,nseg,5,1);
                            kbeg = 1;
                            kend = 1;
                            mdim = jend - jbeg;
                            ndim = iend - ibeg;
                            ldata = lwdat(nbl, nseg, 5);
                            { int iface_lf = 5;
                            lamfix_ns::lamfix(jdim, kdim, idim, FortranArray3DRef<double>(&w(lsnk0), jdim-1, kdim-1, idim-1),
                                              mdim, ndim, FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12),
                                              ibeg, iend, jbeg, jend, kbeg, kend, iface_lf, llev); }
                        }
                    }
                    for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
                        if (std::abs(kbcinfo(nbl,nseg,1,2)) == 2014) {
                            ibeg = kbcinfo(nbl,nseg,2,2);
                            iend = kbcinfo(nbl,nseg,3,2);
                            jbeg = kbcinfo(nbl,nseg,4,2);
                            jend = kbcinfo(nbl,nseg,5,2);
                            kbeg = kdim;
                            kend = kdim;
                            mdim = jend - jbeg;
                            ndim = iend - ibeg;
                            ldata = lwdat(nbl, nseg, 6);
                            { int iface_lf = 6;
                            lamfix_ns::lamfix(jdim, kdim, idim, FortranArray3DRef<double>(&w(lsnk0), jdim-1, kdim-1, idim-1),
                                              mdim, ndim, FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12),
                                              ibeg, iend, jbeg, jend, kbeg, kend, iface_lf, llev); }
                        }
                    }
                    // I-dir:
                    for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
                        if (std::abs(ibcinfo(nbl,nseg,1,1)) == 2014) {
                            ibeg = 1;
                            iend = 1;
                            jbeg = ibcinfo(nbl,nseg,2,1);
                            jend = ibcinfo(nbl,nseg,3,1);
                            kbeg = ibcinfo(nbl,nseg,4,1);
                            kend = ibcinfo(nbl,nseg,5,1);
                            mdim = jend - jbeg;
                            ndim = kend - kbeg;
                            ldata = lwdat(nbl, nseg, 1);
                            { int iface_lf = 1;
                            lamfix_ns::lamfix(jdim, kdim, idim, FortranArray3DRef<double>(&w(lsnk0), jdim-1, kdim-1, idim-1),
                                              mdim, ndim, FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12),
                                              ibeg, iend, jbeg, jend, kbeg, kend, iface_lf, llev); }
                        }
                    }
                    for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
                        if (std::abs(ibcinfo(nbl,nseg,1,2)) == 2014) {
                            ibeg = idim;
                            iend = idim;
                            jbeg = ibcinfo(nbl,nseg,2,2);
                            jend = ibcinfo(nbl,nseg,3,2);
                            kbeg = ibcinfo(nbl,nseg,4,2);
                            kend = ibcinfo(nbl,nseg,5,2);
                            mdim = jend - jbeg;
                            ndim = kend - kbeg;
                            ldata = lwdat(nbl, nseg, 2);
                            { int iface_lf = 2;
                            lamfix_ns::lamfix(jdim, kdim, idim, FortranArray3DRef<double>(&w(lsnk0), jdim-1, kdim-1, idim-1),
                                              mdim, ndim, FortranArray4DRef<double>(&w(ldata), mdim, ndim, 2, 12),
                                              ibeg, iend, jbeg, jend, kbeg, kend, iface_lf, llev); }
                        }
                    }
                }
            }
        }
    } // end ivmx >= 4


    // initialize viscosity arrays
    for (nbl = 1; nbl <= nblock; nbl++) {
        if (mblk2nd(nbl) == myid) {
            lead_ns::lead(nbl, lw, lw2, maxbl);
            initvist_ns::initvist(nbl, jdim, kdim, idim,
                                  FortranArray3DRef<double>(&w(lvis), jdim, kdim, idim), FortranArray4DRef<double>(&w(lxib), jdim, kdim, idim, nummem),
                                  FortranArray3DRef<double>(&w(lsnk0), jdim-1, kdim-1, idim-1), FortranArray3DRef<double>(&w(lcmuv), jdim-1, kdim-1, idim-1),
                                  nummem, FortranArray3DRef<double>(&w(lx), jdim, kdim, idim), FortranArray3DRef<double>(&w(ly), jdim, kdim, idim), FortranArray3DRef<double>(&w(lz), jdim, kdim, idim));
        }
    }

    // initialize RSM turbulence model if needed
    if (ivmx == 72) {
        for (nbl = 1; nbl <= nblock; nbl++) {
            if (mblk2nd(nbl) == myid) {
                lead_ns::lead(nbl, lw, lw2, maxbl);
                {
                    FortranArray4DRef<double> cij_ref = FortranArray4DRef<double>(&wk(lux), jdim, kdim, idim, nummem);
                    module_kwstm_ns::kws_initvist(jdim, kdim, idim, nummem, cij_ref, FortranArray3DRef<double>(&w(lvis), jdim, kdim, idim));
                }
            }
        }
        module_kwstm_ns::kws_read_parm(nou, bou, nbuf, ibufdim);
    }

    // set up interpolation for patched grids
    if (ninter > 0) {
        if (myid == mblk2nd(1)) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " ");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " setting up patched-grid interpolation");
        }
        for (igrid = 1; igrid <= ngrid; igrid++) {
            nbl = nblg(igrid);
            lead_ns::lead(nbl, lw, lw2, maxbl);
            if (myid == mblk2nd(nbl)) {
                { int intrbc_icorr = 0; int intrbc_iset = 0;
                intrbc(FortranArray4DRef<double>(&w(lq), jdim, kdim, idim, ldim), jdim, kdim, idim, nbl, ldim, maxbl, iitot,
                       lig, iipntsg, dxintg, dyintg, dzintg, iiig, jjig, kkig,
                       qb, FortranArray4DRef<double>(&w(lqj0), kdim, idim-1, ldim, 4), FortranArray4DRef<double>(&w(lqk0), jdim, idim-1, ldim, 4), FortranArray4DRef<double>(&w(lqi0), jdim, kdim, ldim, 4),
                       FortranArray4DRef<double>(&wk(1), jdim+1, kdim+1, idim+1, ldim), FortranArray3DRef<double>(&w(lbcj), kdim, idim-1, 2), FortranArray3DRef<double>(&w(lbck), jdim, idim-1, 2), FortranArray3DRef<double>(&w(lbci), jdim, kdim, 2),
                       nou, bou, nbuf, ibufdim, intrbc_icorr, intrbc_iset); }
            }
        }
    }

    // use overlapped grid interpolation data to set conditions in overlapped
    // cells (chimera). This setup pass fills the qb buffer so the first
    // bc_xmera/xupdt scatters interpolated donor values (not 0) into the fringe
    // cells before the initial bcchk. (Fortran setup.F do-2100 loop; was
    // MISSING in the C++ translation, leaving qb=0 -> fringe density zeroed ->
    // "negative density" abort at cycle 1 on overset grids.)
    for (igrid = 1; igrid <= ngrid; igrid++) {
        ncg = ncgg(igrid);
        nbl = nblg(igrid) - 1;
        for (int icg = 1; icg <= ncg + 1; icg++) {
            nbl = nbl + 1;
            if (iovrlp(nbl) == 1 && myid == mblk2nd(nbl)) {
                lead_ns::lead(nbl, lw, lw2, maxbl);
                ldim = 5;
                { int ic = irghost, is = 1;
                intrbc(FortranArray4DRef<double>(&w(lq), jdim, kdim, idim, ldim), jdim, kdim, idim, nbl, ldim, maxbl, iitot,
                       lig, iipntsg, dxintg, dyintg, dzintg, iiig, jjig, kkig,
                       qb, FortranArray4DRef<double>(&w(lqj0), kdim, idim-1, ldim, 4), FortranArray4DRef<double>(&w(lqk0), jdim, idim-1, ldim, 4), FortranArray4DRef<double>(&w(lqi0), jdim, kdim, ldim, 4),
                       FortranArray4DRef<double>(&wk(1), jdim+1, kdim+1, idim+1, ldim), FortranArray3DRef<double>(&w(lbcj), kdim, idim-1, 2), FortranArray3DRef<double>(&w(lbck), jdim, idim-1, 2), FortranArray3DRef<double>(&w(lbci), jdim, kdim, 2),
                       nou, bou, nbuf, ibufdim, ic, is); }
                if (iviscg(nbl,1) >= 2 || iviscg(nbl,2) >= 2 || iviscg(nbl,3) >= 2) {
                    ldim = 1;
                    { int ic = irghost, is = 2;
                    intrbc(FortranArray4DRef<double>(&w(lvis), jdim, kdim, idim, ldim), jdim, kdim, idim, nbl, ldim, maxbl, iitot,
                           lig, iipntsg, dxintg, dyintg, dzintg, iiig, jjig, kkig,
                           qb, FortranArray4DRef<double>(&w(lvj0), kdim, idim-1, ldim, 4), FortranArray4DRef<double>(&w(lvk0), jdim, idim-1, ldim, 4), FortranArray4DRef<double>(&w(lvi0), jdim, kdim, ldim, 4),
                           FortranArray4DRef<double>(&wk(1), jdim+1, kdim+1, idim+1, ldim), FortranArray3DRef<double>(&w(lbcj), kdim, idim-1, 2), FortranArray3DRef<double>(&w(lbck), jdim, idim-1, 2), FortranArray3DRef<double>(&w(lbci), jdim, kdim, 2),
                           nou, bou, nbuf, ibufdim, ic, is); }
                }
                if (iviscg(nbl,1) >= 4 || iviscg(nbl,2) >= 4 || iviscg(nbl,3) >= 4) {
                    ldim = nummem;
                    { int ic = irghost, is = 3;
                    intrbc(FortranArray4DRef<double>(&w(lxib), jdim, kdim, idim, ldim), jdim, kdim, idim, nbl, ldim, maxbl, iitot,
                           lig, iipntsg, dxintg, dyintg, dzintg, iiig, jjig, kkig,
                           qb, FortranArray4DRef<double>(&w(ltj0), kdim, idim-1, ldim, 4), FortranArray4DRef<double>(&w(ltk0), jdim, idim-1, ldim, 4), FortranArray4DRef<double>(&w(lti0), jdim, kdim, ldim, 4),
                           FortranArray4DRef<double>(&wk(1), jdim+1, kdim+1, idim+1, ldim), FortranArray3DRef<double>(&w(lbcj), kdim, idim-1, 2), FortranArray3DRef<double>(&w(lbck), jdim, idim-1, 2), FortranArray3DRef<double>(&w(lbci), jdim, kdim, 2),
                           nou, bou, nbuf, ibufdim, ic, is); }
                }
            }
        }
    }

    // apply boundary conditions to initialize ghost cells
    if (myid == mblk2nd(1)) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " applying boundary conditions");
    }

    for (igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nblg(igrid);
        ncg = ncgg(igrid);
        for (n = 1; n <= ncg+1; n++) {
            if (myid == mblk2nd(nbl)) {
                lead_ns::lead(nbl, lw, lw2, maxbl);
                idum = 0;
                nttuse = std::max(ntt-1, 1);
                bc_ns::bc(idum, nbl, lw, lw2, w, mgwk, wk, nwork, clw(nttuse),
                          nou, bou, nbuf, ibufdim, maxbl, maxgr, maxseg, itrans,
                          irotat, idefrm, igridg, nblg, nbci0, nbcj0, nbck0,
                          nbcidim, nbcjdim, nbckdim, ibcinfo, jbcinfo,
                          kbcinfo, bcfilei, bcfilej, bcfilek, lwdat, myid,
                          idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
            }
            nbl = nbl + 1;
        }
        nblout = nblg(igrid);
        { int iunit11 = 11; writ_buf_ns::writ_buf(nblout, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }
    }

    // print final initialization message
    if (myid == mblk2nd(1)) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "***** INITIALIZATION COMPLETE *****");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
    }

    { int iunit11 = 11; int nbl1 = 1; writ_buf_ns::writ_buf(nbl1, iunit11, nou, bou, nbuf, ibufdim, myhost, myid, mycomm, mblk2nd, maxbl); }

} // end setup()

// Forward implementations for procedures defined in other source files
// These are declared in setup.h but implemented in separate Fortran files

void cellvol(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> sj, FortranArray3DRef<double> sk, FortranArray3DRef<double> si, FortranArray3DRef<double> vol, FortranArray2DRef<double> t, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, int& nbl, int& iflagv1, int& iflagv, int& imin, int& imax, int& jmin, int& jmax, int& kmin, int& kmax)
{
    // Implemented in cellvol.f / cellvol.cpp
}

void dird(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> snj0, FortranArray3DRef<double> snk0, FortranArray3DRef<double> sni0, FortranArray3DRef<double> snjm, FortranArray3DRef<double> snkm, FortranArray3DRef<double> snim, FortranArray4DRef<double> w, FortranArray1DRef<int> ivisc, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    // Implemented in dird.f / dird.cpp
}

void findmin_new(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, FortranArray1DRef<int> iwk, int& iwork, int& nsurf, int& j1, FortranArray1DRef<double> xs, FortranArray1DRef<int> ixs, FortranArray1DRef<int> lsminn, FortranArray1DRef<int> ireq_xs, FortranArray1DRef<int> ireq_bb, int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, FortranArray1DRef<int> nblg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg, FortranArray1DRef<int> mblk2nd)
{
    // Implemented in findmin_new.f / findmin_new.cpp
}

void global(int& myid, int& maxbl, int& maxgr, int& maxseg, int& maxcs, int& nplots, int& mxbli, FortranArray4DRef<double> bcvali, FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, FortranArray2DRef<int> nblk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, FortranArray2DRef<double> rkap0g, int& nblock, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<int> idiagg, FortranArray1DRef<int> nblcg, FortranArray2DRef<int> idegg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray2DRef<int> mit, FortranArray1DRef<int> ilamlog, FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamlog, FortranArray1DRef<int> klamhig, FortranArray2DRef<int> iwfg, FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx, FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz, FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr, FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg, FortranArray1DRef<int> iovrlp, int& ninter, int& nplot3d, FortranArray2DRef<int> inpl3d, int& ip3dsurf, int& nprint, FortranArray2DRef<int> inpr, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> iforce, int& lfgm, int& ncs, FortranArray2DRef<int> icsinfo, int& ihstry, int& ncycmax, FortranArray1DRef<int> iv, FortranArray1DRef<double> time2, FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl, int& intmax, int& nsub1, FortranArray2DRef<int> iindex, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, int& icall, int& iunit11, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& ibufdim, int& nbuf, FortranArray1DRef<int> mglevg, FortranArray1DRef<int> nemgl, FortranArray2DRef<int> ipl3dtmp, int& ntr, FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae, FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae, FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass, FortranArray2DRef<double> damp, FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0, int& nmds, int& maxaes, FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg, FortranArray2DRef<double> xorgae0, FortranArray2DRef<double> yorgae0, FortranArray2DRef<double> zorgae0, FortranArray2DRef<int> icouple, int& iprnsurf)
{
    // Implemented in global.f / global.cpp
}

void intrbc(FortranArray4DRef<double> q, int& jdim, int& kdim, int& idim, int& nbl, int& ldim, int& maxbl, int& iitot, FortranArray1DRef<int> lig, FortranArray1DRef<int> iipntsg, FortranArray1DRef<double> dxintg, FortranArray1DRef<double> dyintg, FortranArray1DRef<double> dzintg, FortranArray1DRef<int> iiig, FortranArray1DRef<int> jjig, FortranArray1DRef<int> kkig, FortranArray3DRef<double> qb, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qq, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& icorr, int& iset)
{
    // Forward to the real implementation. Previously an EMPTY stub, which C++
    // unqualified lookup resolved for the setup-stage intrbc calls (patched-
    // grid + overset chimera qb fill), so qb was never filled -> overset
    // fringe cells zeroed -> "negative density" abort. (bug class #7)
    intrbc_ns::intrbc(q, jdim, kdim, idim, nbl, ldim, maxbl, iitot, lig, iipntsg,
                      dxintg, dyintg, dzintg, iiig, jjig, kkig, qb, qj0, qk0, qi0,
                      qq, bcj, bck, bci, nou, bou, nbuf, ibufdim, icorr, iset);
}

void metric(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> sj, FortranArray3DRef<double> sk, FortranArray3DRef<double> si, FortranArray2DRef<double> t, FortranArray3DRef<double> t1, int& nbl, int& iflag, int& icnt, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& maxbl, int& maxseg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd)
{
    // Implemented in metric.f / metric.cpp
}

void update(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray2DRef<double> qj0, FortranArray2DRef<double> qk0, FortranArray2DRef<double> qi0, FortranArray2DRef<double> sj, FortranArray2DRef<double> sk, FortranArray2DRef<double> si, FortranArray1DRef<double> vol, FortranArray1DRef<double> dtj, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> blank, FortranArray1DRef<double> x, FortranArray1DRef<double> y, FortranArray1DRef<double> z, FortranArray4DRef<double> res, FortranArray1DRef<double> wk0, FortranArray3DRef<double> vmuk, FortranArray3DRef<double> vmuj, FortranArray3DRef<double> vmui, FortranArray1DRef<double> wk, int& nwork, int& nbl, int& iover, FortranArray4DRef<double> vk0, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, FortranArray3DRef<double> volk0, FortranArray4DRef<double> tursav, FortranArray4DRef<double> tk0, FortranArray3DRef<double> cmuv, FortranArray1DRef<int> iadvance, int& nummem, FortranArray4DRef<double> ux)
{
    // Implemented in update.f / update.cpp
}

void updatedg(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& iupdat, int& iseqr, int& maxbl, int& maxgr, int& maxseg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& nblock, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae, FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae, FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass, FortranArray2DRef<double> damp, FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0, int& nmds, int& maxaes, FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray3DRef<int> islavept, int& nslave, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf, int& maxsegdg, FortranArray1DRef<int> iwk, int& nmaster, int& nt, FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig, FortranArray2DRef<double> xorgae0, FortranArray2DRef<double> yorgae0, FortranArray2DRef<double> zorgae0, FortranArray2DRef<int> icouple, FortranArray1DRef<int> ireq, int& nnodes, FortranArray2DRef<int> nblelst, FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax, FortranArray1DRef<int> kskmax, FortranArray1DRef<double> ue)
{
    // Implemented in updatedg.f / updatedg.cpp
}

} // namespace setup_ns
