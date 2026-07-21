// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// writ_buf.cpp — write the internal buffer for block nbl to unit iunit.
// Sequential (non-FASTIO, non-DIST_MPI) build. bou(ibufdim,nbuf), nou(nbuf),
// mblk2nd(maxbl). Only the single-node path remains.
#include "ron_common.h"

void writ_buf(int /*nbl*/,int iunit,int* nou,FStr* bou,int nbuf,int ibufdim,
              int /*myhost*/,int /*myid*/,int /*mycomm*/,int* /*mblk2nd*/,
              int /*maxbl*/){
    int iou=1;
    if (iunit==11) iou=1;
    if (iunit== 9) iou=2;
    if (iunit==14) iou=3;
    if (iunit==25) iou=4;

    if (nou[iou-1]>0){
        for(int kou=1;kou<=nou[iou-1];++kou)
            outbuf(bou[(kou-1)+(iou-1)*ibufdim],iunit);

        if (nou[iou-1]>=ibufdim){
            FILE* f=FUNIT(iunit);
            if(f){
                std::fprintf(f,"WARNING: internal buffer length exceeded -"
                               "  make parameter ibufdim > %6d\n",nou[iou-1]);
                std::fprintf(f,"continuing, but you will be missing"
                               " some output data \n");
            }
        }
        for(int kou=1;kou<=nou[iou-1];++kou)
            bou[(kou-1)+(iou-1)*ibufdim]=FStr(120);
        nou[iou-1]=0;
    }
    (void)nbuf;
}
