// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "pre_blockbc.h"
#include "termn8.h"
#include <cstring>
#include <cstdio>
#include <algorithm>

namespace pre_blockbc_ns {

void pre_blockbc(int& nbl, FortranArray2DRef<int> lw, int& icount,
                 FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
                 FortranArray1DRef<int> kdimg, FortranArray2DRef<int> isav_blk,
                 FortranArray2DRef<int> nblk, int& nbli, FortranArray3DRef<int> limblk,
                 FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, int& mxbli,
                 FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
                 int& nbuf, int& ibufdim, int& myid, int& maxbl, int& ierrflg)
{
    // local variables
    int n, it, ir, itime, iti;
    int ic_blk, in_blk;
    int idimn, jdimn, kdimn, idimc, jdimc, kdimc;
    int jface, jedge, iedge, lwt;
    int iss, ise, jss, jse, kss, kse;

    if (std::abs(nbli) > 0) {
//
//     loop over all 1-1 interfaces
//
        for (n = 1; n <= std::abs(nbli); n++) {
//
//           check to see if block interface n is turned on
//
            if (nblon(n) >= 0) {
//
//               if block interface n is a boundary of current block then
//               set up blocking parameters
//
                if (nbl == nblk(1,n) || nbl == nblk(2,n)) {
//
//                  it = 1  if nblk(1,n) is the block being advanced
//                  it = 2  if nblk(2,n) is the block being advanced
//                  ir = 1  if nblk(1,n) is the neighboring block
//                  ir = 2  if nblk(2,n) is the neighboring block
//
                    it = 1;
                    ir = 2;
                    if (nbl == nblk(2,n)) it = 2;
                    if (nbl == nblk(2,n)) ir = 1;
//
//                  allow for 1-1 blocking in same grid
//
                    itime = 1;
                    if (nblk(1,n) == nblk(2,n)) itime = 2;
                    for (iti = 1; iti <= itime; iti++) {
                        if (iti > 1) {
                            it = 1;
                            ir = 2;
                        }
//
//                       define current and neighbor blocks
//
                        ic_blk = nblk(it,n);
                        in_blk = nblk(ir,n);
//
//                       set dimensions of current blocks involved
//
                        idimn = idimg(in_blk);
                        jdimn = jdimg(in_blk);
                        kdimn = kdimg(in_blk);
                        idimc = idimg(ic_blk);
                        jdimc = jdimg(ic_blk);
                        kdimc = kdimg(ic_blk);
//
                        if (isva(ir,1,n) + isva(ir,2,n) == 3) {
                            jface = 3;
                            jedge = 1;
                            if (limblk(ir,3,n) != 1) jedge = kdimg(in_blk);
                        }
                        if (isva(ir,1,n) + isva(ir,2,n) == 4) {
                            jface = 2;
                            jedge = 1;
                            if (limblk(ir,2,n) != 1) jedge = jdimg(in_blk);
                        }
                        if (isva(ir,1,n) + isva(ir,2,n) == 5) {
                            jface = 1;
                            jedge = 1;
                            if (limblk(ir,1,n) != 1) jedge = idimg(in_blk);
                        }
//
//                       k = constant interface
//
                        if (isva(it,1,n) + isva(it,2,n) == 3) {
                            lwt   = lw(3,ic_blk);
                            iedge = 1;
                            if (limblk(it,3,n) != 1) {
                                lwt   = lwt + jdimc*(idimc-1)*5*2;
                                iedge = 2;
                            }
                            iss = limblk(it,1,n);
                            ise = limblk(it,4,n);
                            jss = limblk(it,2,n);
                            jse = limblk(it,5,n);
                            kss = limblk(it,3,n);
                            kse = limblk(it,6,n);
                            if (iss > ise) {
                                iss = iss + 1;
                            } else {
                                ise = ise + 1;
                            }
                            if (jss > jse) {
                                jss = jss + 1;
                            } else {
                                jse = jse + 1;
                            }
//
//                       j = constant interface
//
                        } else if (isva(it,1,n) + isva(it,2,n) == 4) {
                            lwt   = lw(2,ic_blk);
                            iedge = 1;
                            if (limblk(it,2,n) != 1) {
                                lwt   = lwt + kdimc*(idimc-1)*5*2;
                                iedge = 2;
                            }
                            iss = limblk(it,1,n);
                            ise = limblk(it,4,n);
                            jss = limblk(it,2,n);
                            jse = limblk(it,5,n);
                            kss = limblk(it,3,n);
                            kse = limblk(it,6,n);
                            if (iss > ise) {
                                iss = iss + 1;
                            } else {
                                ise = ise + 1;
                            }
                            if (kss > kse) {
                                kss = kss + 1;
                            } else {
                                kse = kse + 1;
                            }
//
//                       i = constant interface
//
                        } else if (isva(it,1,n) + isva(it,2,n) == 5) {
                            lwt   = lw(4,ic_blk);
                            iedge = 1;
                            if (limblk(it,1,n) != 1) {
                                lwt   = lwt + jdimc*kdimc*5*2;
                                iedge = 2;
                            }
                            iss = limblk(it,1,n);
                            ise = limblk(it,4,n);
                            kss = limblk(it,3,n);
                            kse = limblk(it,6,n);
                            jss = limblk(it,2,n);
                            jse = limblk(it,5,n);
                            if (kss > kse) {
                                kss = kss + 1;
                            } else {
                                kse = kse + 1;
                            }
                            if (jss > jse) {
                                jss = jss + 1;
                            } else {
                                jse = jse + 1;
                            }
                        }
//
//                       put appropriate values into isav_blk array
//
                        icount = icount + 1;
                        if (icount > 2*mxbli) {
                            nou(1) = std::min(nou(1)+1, ibufdim);
                            std::snprintf(bou(nou(1),1), 120,
                                " err: icount in pre_bcblk = %d", icount);
                            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                        }
                        isav_blk(icount,1)  = n;
                        isav_blk(icount,2)  = it;
                        isav_blk(icount,3)  = ir;
                        isav_blk(icount,4)  = ic_blk;
                        isav_blk(icount,5)  = in_blk;
                        isav_blk(icount,6)  = jface;
                        isav_blk(icount,7)  = jedge;
                        isav_blk(icount,8)  = isva(it,1,n) + isva(it,2,n);
                        isav_blk(icount,9)  = lwt;
                        isav_blk(icount,10) = iedge;
                        isav_blk(icount,11) = iss;
                        isav_blk(icount,12) = ise;
                        isav_blk(icount,13) = jss;
                        isav_blk(icount,14) = jse;
                        isav_blk(icount,15) = kss;
                        isav_blk(icount,16) = kse;
                        isav_blk(icount,17) = iti;

                    } // end do 101 iti
                } // end if nbl == nblk
            } // end if nblon(n) >= 0
        } // end do 100 n

    } // end if abs(nbli) > 0

    return;
}

} // namespace pre_blockbc_ns
