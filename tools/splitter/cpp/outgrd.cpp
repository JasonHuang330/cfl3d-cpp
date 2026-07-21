// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "outgrd.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>

namespace outgrd_ns {

void outgrd(int& nf, int& ib, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& ni, int& nj, int& nk, int& i0, int& i1, int& isk, int& j0, int& j1, int& jsk, int& k0, int& k1, int& ksk, int& ip3dgrd)
{
    if (ip3dgrd > 0) {
        if (ib == 0) {
            // write(nf,*) (((x(i,j,k),i=i0,i1,isk),j=j0,j1,jsk),k=k0,k1,ksk),
            //             (((y(i,j,k),i=i0,i1,isk),j=j0,j1,jsk),k=k0,k1,ksk),
            //             (((z(i,j,k),i=i0,i1,isk),j=j0,j1,jsk),k=k0,k1,ksk)
            FILE* fp = fortran_get_unit(nf);
            for (int k = k0; k <= k1; k += ksk)
                for (int j = j0; j <= j1; j += jsk)
                    for (int i = i0; i <= i1; i += isk)
                        fprintf(fp, " %g", x(i, j, k));
            for (int k = k0; k <= k1; k += ksk)
                for (int j = j0; j <= j1; j += jsk)
                    for (int i = i0; i <= i1; i += isk)
                        fprintf(fp, " %g", y(i, j, k));
            for (int k = k0; k <= k1; k += ksk)
                for (int j = j0; j <= j1; j += jsk)
                    for (int i = i0; i <= i1; i += isk)
                        fprintf(fp, " %g", z(i, j, k));
            fprintf(fp, "\n");
        } else {
            // write(nf) (((x(i,j,k),i=i0,i1,isk),j=j0,j1,jsk),k=k0,k1,ksk),
            //           (((y(i,j,k),i=i0,i1,isk),j=j0,j1,jsk),k=k0,k1,ksk),
            //           (((z(i,j,k),i=i0,i1,isk),j=j0,j1,jsk),k=k0,k1,ksk)
            FILE* fp = fortran_get_unit(nf);
            // Count elements
            int ni_out = 0, nj_out = 0, nk_out = 0;
            for (int i = i0; i <= i1; i += isk) ni_out++;
            for (int j = j0; j <= j1; j += jsk) nj_out++;
            for (int k = k0; k <= k1; k += ksk) nk_out++;
            int nelem = ni_out * nj_out * nk_out;
            int reclen = 3 * nelem * (int)sizeof(double);
            fwrite(&reclen, sizeof(int), 1, fp);
            for (int k = k0; k <= k1; k += ksk)
                for (int j = j0; j <= j1; j += jsk)
                    for (int i = i0; i <= i1; i += isk)
                        fwrite(&x(i, j, k), sizeof(double), 1, fp);
            for (int k = k0; k <= k1; k += ksk)
                for (int j = j0; j <= j1; j += jsk)
                    for (int i = i0; i <= i1; i += isk)
                        fwrite(&y(i, j, k), sizeof(double), 1, fp);
            for (int k = k0; k <= k1; k += ksk)
                for (int j = j0; j <= j1; j += jsk)
                    for (int i = i0; i <= i1; i += isk)
                        fwrite(&z(i, j, k), sizeof(double), 1, fp);
            fwrite(&reclen, sizeof(int), 1, fp);
        }
    } else {
        if (ib == 0) {
            // write(nf,*) (((x(i,j,k),j=j0,j1,jsk),k=k0,k1,ksk),i=i0,i1,isk),
            //             (((y(i,j,k),j=j0,j1,jsk),k=k0,k1,ksk),i=i0,i1,isk),
            //             (((z(i,j,k),j=j0,j1,jsk),k=k0,k1,ksk),i=i0,i1,isk)
            FILE* fp = fortran_get_unit(nf);
            for (int i = i0; i <= i1; i += isk)
                for (int k = k0; k <= k1; k += ksk)
                    for (int j = j0; j <= j1; j += jsk)
                        fprintf(fp, " %g", x(i, j, k));
            for (int i = i0; i <= i1; i += isk)
                for (int k = k0; k <= k1; k += ksk)
                    for (int j = j0; j <= j1; j += jsk)
                        fprintf(fp, " %g", y(i, j, k));
            for (int i = i0; i <= i1; i += isk)
                for (int k = k0; k <= k1; k += ksk)
                    for (int j = j0; j <= j1; j += jsk)
                        fprintf(fp, " %g", z(i, j, k));
            fprintf(fp, "\n");
        } else {
            // write(nf) (((x(i,j,k),j=j0,j1,jsk),k=k0,k1,ksk),i=i0,i1,isk),
            //           (((y(i,j,k),j=j0,j1,jsk),k=k0,k1,ksk),i=i0,i1,isk),
            //           (((z(i,j,k),j=j0,j1,jsk),k=k0,k1,ksk),i=i0,i1,isk)
            FILE* fp = fortran_get_unit(nf);
            // Count elements
            int ni_out = 0, nj_out = 0, nk_out = 0;
            for (int i = i0; i <= i1; i += isk) ni_out++;
            for (int j = j0; j <= j1; j += jsk) nj_out++;
            for (int k = k0; k <= k1; k += ksk) nk_out++;
            int nelem = ni_out * nj_out * nk_out;
            int reclen = 3 * nelem * (int)sizeof(double);
            fwrite(&reclen, sizeof(int), 1, fp);
            for (int i = i0; i <= i1; i += isk)
                for (int k = k0; k <= k1; k += ksk)
                    for (int j = j0; j <= j1; j += jsk)
                        fwrite(&x(i, j, k), sizeof(double), 1, fp);
            for (int i = i0; i <= i1; i += isk)
                for (int k = k0; k <= k1; k += ksk)
                    for (int j = j0; j <= j1; j += jsk)
                        fwrite(&y(i, j, k), sizeof(double), 1, fp);
            for (int i = i0; i <= i1; i += isk)
                for (int k = k0; k <= k1; k += ksk)
                    for (int j = j0; j <= j1; j += jsk)
                        fwrite(&z(i, j, k), sizeof(double), 1, fp);
            fwrite(&reclen, sizeof(int), 1, fp);
        }
    }
}

} // namespace outgrd_ns
