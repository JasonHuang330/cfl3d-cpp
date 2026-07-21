// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "compg2n.h"
#include "termn8.h"
#include "runtime/fortran_io.h"

namespace compg2n_ns {

void compg2n(int& nblock, int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> nblcg, int& nodes, FortranArray2DRef<int> iwrk, int& myid, int& myhost, FortranArray1DRef<int> mblk2nd, int& mycomm, int& maxgr, int& maxbl, int& ierrflg, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou)
{
    //
    // Purpose: To assign blocks to nodes based on size considerations
    // and multigrid structure
    //

    // simple case: 1 node
    if (nodes == 1) {
        for (int ibloc = 1; ibloc <= nblock; ibloc++) {
            mblk2nd(ibloc) = 1;
        }
        return;
    }

    // zero out block assignments
    for (int ibloc = 1; ibloc <= nblock; ibloc++) {
        mblk2nd(ibloc) = 0;
        iwrk(ibloc, 3) = 0;
    }

    // find parent grid sizes and blocks associated with coarser grids
    // sizes are in terms of number of cells
    for (int igrid = 1; igrid <= ngrid; igrid++) {
        int ibloc = nblg(igrid);
        iwrk(igrid, 1) = igrid;
        iwrk(igrid, 2) =
            (idimg(ibloc) - 1) * (jdimg(ibloc) - 1) * (kdimg(ibloc) - 1);
    }

    // sort grids from largest to smallest
    for (int igrid = 1; igrid <= ngrid - 1; igrid++) {
        for (int jgrid = igrid + 1; jgrid <= ngrid; jgrid++) {
            if (iwrk(igrid, 2) < iwrk(jgrid, 2)) {
                int itmp        = iwrk(igrid, 1);
                iwrk(igrid, 1)  = iwrk(jgrid, 1);
                iwrk(jgrid, 1)  = itmp;
                itmp            = iwrk(igrid, 2);
                iwrk(igrid, 2)  = iwrk(jgrid, 2);
                iwrk(jgrid, 2)  = itmp;
            }
        }
    }

    // option 1: allow embedded blocks to reside on any processor in
    //           order to load balance (some cases may not work
    //           correctly, and is just about as slow as option 2)
    //
    // option 2: restrict embedded blocks to the same processor as
    //           the block in which they are embedded (many cases
    //           won't load balance very well).
    //
    // option 2 is recommended at this point, 1 is only for testing

    int ioption = 2;

    if (ioption == 1) {

        // assign global-level grids in order (largest to smallest) to
        // the nodes with the fewest points
        for (int jgrid = 1; jgrid <= ngrid; jgrid++) {
            int igrid = iwrk(jgrid, 1);
            int ibloc = nblg(igrid);
            if (mblk2nd(ibloc) == 0) {
                int imin = 1;
                for (int inode = 2; inode <= nodes; inode++) {
                    if (iwrk(inode, 3) < iwrk(imin, 3)) imin = inode;
                }
                mblk2nd(ibloc) = imin;
                iwrk(imin, 3) = iwrk(imin, 3) + iwrk(jgrid, 2);
            }
        }

    } else {

        // assign global-level grids in order (largest to smallest) to
        // the nodes with the fewest points, skipping embedded blocks
        for (int jgrid = 1; jgrid <= ngrid; jgrid++) {
            int igrid = iwrk(jgrid, 1);
            int ibloc = nblg(igrid);
            if (nblcg(ibloc) == ibloc) {
                // block is not embedded
                if (mblk2nd(ibloc) == 0) {
                    int imin = 1;
                    for (int inode = 2; inode <= nodes; inode++) {
                        if (iwrk(inode, 3) < iwrk(imin, 3)) imin = inode;
                    }
                    mblk2nd(ibloc) = imin;
                    iwrk(imin, 3) = iwrk(imin, 3) + iwrk(jgrid, 2);
                }
            }
        }

        // assign each embedded block to the same node as the block
        // into which it is embedded
        for (int jgrid = 1; jgrid <= ngrid; jgrid++) {
            int igrid = iwrk(jgrid, 1);
            int ibloc = nblg(igrid);
            if (nblcg(ibloc) != ibloc) {
                // block is embedded
                mblk2nd(ibloc) = mblk2nd(nblcg(ibloc));
            }
        }

    }

    // assign each coarser-level grid to the same node as its parent
    int nbl = 0;
    for (int igrid = 1; igrid <= ngrid; igrid++) {
        nbl = nbl + 1;
        int ibloc = nblg(igrid);
        if (nbl != ibloc) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, " block numbers mismatch, in compg2n\n");
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
        int nodel = mblk2nd(ibloc);
        int ncg   = ncgg(igrid);
        if (ncg > 0) {
            for (int n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                mblk2nd(nbl) = nodel;
            }
        }
    }

    return;
}

} // namespace compg2n_ns
