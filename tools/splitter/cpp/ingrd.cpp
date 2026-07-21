// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "ingrd.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstdlib>

namespace ingrd_ns {

// ingrd: reads grid coordinates (x, y, z) from file unit nf.
//   ib == 0  => formatted (list-directed) read
//   ib != 0  => unformatted (binary) read
//   ip3dgrd > 0  => read order: i varies fastest (i=1..ni, j=1..nj, k=1..nk)
//   ip3dgrd <= 0 => read order: j varies fastest (j=1..nj, k=1..nk, i=1..ni)
void ingrd(int& nf, int& ib, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& ni, int& nj, int& nk, int& ip3dgrd)
{
    int total = ni * nj * nk;

    if (ip3dgrd > 0) {
        if (ib == 0) {
            // Formatted list-directed read, order: (i=1..ni, j=1..nj, k=1..nk)
            // Read all x values, then all y values, then all z values.
            // FortranArray3D stores in column-major (i fastest), matching this order.
            // Use a flat temp buffer and fortran_read_list_array.
            FortranArray1D<double> tmp(total);
            fortran_read_list_array(nf, tmp, 1, total);
            int idx = 1;
            for (int k = 1; k <= nk; k++)
                for (int j = 1; j <= nj; j++)
                    for (int i = 1; i <= ni; i++)
                        x(i, j, k) = tmp(idx++);

            fortran_read_list_array(nf, tmp, 1, total);
            idx = 1;
            for (int k = 1; k <= nk; k++)
                for (int j = 1; j <= nj; j++)
                    for (int i = 1; i <= ni; i++)
                        y(i, j, k) = tmp(idx++);

            fortran_read_list_array(nf, tmp, 1, total);
            idx = 1;
            for (int k = 1; k <= nk; k++)
                for (int j = 1; j <= nj; j++)
                    for (int i = 1; i <= ni; i++)
                        z(i, j, k) = tmp(idx++);
        } else {
            // Unformatted (binary) read, order: (i=1..ni, j=1..nj, k=1..nk)
            // Fortran unformatted record: [4-byte len][data][4-byte len]
            FILE* fp = fortran_get_unit(nf);
            int reclen;
            // Read record: x, y, z all in one unformatted record
            fread(&reclen, sizeof(int), 1, fp);
            for (int k = 1; k <= nk; k++)
                for (int j = 1; j <= nj; j++)
                    for (int i = 1; i <= ni; i++)
                        fread(&x(i, j, k), sizeof(double), 1, fp);
            for (int k = 1; k <= nk; k++)
                for (int j = 1; j <= nj; j++)
                    for (int i = 1; i <= ni; i++)
                        fread(&y(i, j, k), sizeof(double), 1, fp);
            for (int k = 1; k <= nk; k++)
                for (int j = 1; j <= nj; j++)
                    for (int i = 1; i <= ni; i++)
                        fread(&z(i, j, k), sizeof(double), 1, fp);
            fread(&reclen, sizeof(int), 1, fp);
        }
    } else {
        if (ib == 0) {
            // Formatted list-directed read, order: (j=1..nj, k=1..nk, i=1..ni)
            FortranArray1D<double> tmp(total);
            fortran_read_list_array(nf, tmp, 1, total);
            int idx = 1;
            for (int i = 1; i <= ni; i++)
                for (int k = 1; k <= nk; k++)
                    for (int j = 1; j <= nj; j++)
                        x(i, j, k) = tmp(idx++);

            fortran_read_list_array(nf, tmp, 1, total);
            idx = 1;
            for (int i = 1; i <= ni; i++)
                for (int k = 1; k <= nk; k++)
                    for (int j = 1; j <= nj; j++)
                        y(i, j, k) = tmp(idx++);

            fortran_read_list_array(nf, tmp, 1, total);
            idx = 1;
            for (int i = 1; i <= ni; i++)
                for (int k = 1; k <= nk; k++)
                    for (int j = 1; j <= nj; j++)
                        z(i, j, k) = tmp(idx++);
        } else {
            // Unformatted (binary) read, order: (j=1..nj, k=1..nk, i=1..ni)
            FILE* fp = fortran_get_unit(nf);
            int reclen;
            fread(&reclen, sizeof(int), 1, fp);
            for (int i = 1; i <= ni; i++)
                for (int k = 1; k <= nk; k++)
                    for (int j = 1; j <= nj; j++)
                        fread(&x(i, j, k), sizeof(double), 1, fp);
            for (int i = 1; i <= ni; i++)
                for (int k = 1; k <= nk; k++)
                    for (int j = 1; j <= nj; j++)
                        fread(&y(i, j, k), sizeof(double), 1, fp);
            for (int i = 1; i <= ni; i++)
                for (int k = 1; k <= nk; k++)
                    for (int j = 1; j <= nj; j++)
                        fread(&z(i, j, k), sizeof(double), 1, fp);
            fread(&reclen, sizeof(int), 1, fp);
        }
    }
}

} // namespace ingrd_ns
