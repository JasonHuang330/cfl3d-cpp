// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// termn8.cpp — flush/close output, write the error flag to unit 99, terminate.
// Sequential build. bou is character*120 bou(ibufdim,nbuf); nou(nbuf).
#include "ron_common.h"

void termn8(int /*myid_stop*/,int ierrflg,int ibufdim,int nbuf,FStr* bou,int* nou){
    FILE* f99=FUNIT(99);
    // format 99: ' error code:',/,i4
    if(f99) std::fprintf(f99," error code:\n%s\n",f_iw(ierrflg,4).c_str());

    if (ierrflg==0){
        if(f99) std::fprintf(f99,"\n execution terminated normally\n");
    } else {
        if (ierrflg!=-999){
            if (ierrflg==-99){
                if (nou[0]>0){
                    if(f99) std::fprintf(f99,
                      "\n abnormal termination during array sizing (preronnie)\n"
                      " (error message follows)\n\n");
                } else {
                    if(f99) std::fprintf(f99,
                      "\n abnormal termination during array sizing (preronnie)\n"
                      " (see preronnie.out for error message)\n\n");
                }
            } else {
                int iflag=0;
                for(int nn=0;nn<nbuf;++nn) if(nou[nn]>0) iflag=1;
                if (iflag>0){
                    if(f99) std::fprintf(f99,
                      "\n abnormal termination due to cfl3d error check\n"
                      " (error message follows)\n\n");
                } else {
                    if(f99) std::fprintf(f99,
                      "\n abnormal termination due to cfl3d error check\n"
                      " (see main output file for error message)\n"
                      " (also check dynamic patch output file if applicable)\n\n");
                }
            }
        } else {
            if(f99) std::fprintf(f99,
              "\n abnormal termination due to receipt of system signal\n"
              " (kill, floating pt. exception, segmentation fault, etc.)\n\n");
        }
        // dump all internal buffers (list-directed write(99,*) => leading blank)
        for(int nn=1;nn<=nbuf;++nn){
            if (nou[nn-1]>0){
                if(f99) std::fprintf(f99," \n");
                if (nn==1){
                    if(f99) std::fprintf(f99," dump of unit 11 (main output) buffer:\n");
                } else if (nn==2){
                    if(f99) std::fprintf(f99," dump of unit 09 (fort.9, dyn. patch) buffer:\n");
                } else if (nn==3){
                    if(f99) std::fprintf(f99," dump of unit 14 (baldwin-lomax) buffer:\n");
                } else if (nn==4){
                    if(f99) std::fprintf(f99," dump of unit 25 (cfl3d.dynamic_patch) buffer:\n");
                } else {
                    if(f99) std::fprintf(f99," need a message for buffer          %d\n",nn);
                }
                if(f99) std::fprintf(f99," \n");
                for(int kou=1;kou<=nou[nn-1];++kou)
                    outbuf(bou[(kou-1)+(nn-1)*ibufdim],99);
            }
        }
    }

    my_flush(25);
    my_flush(99);
    if(FUNIT(25)){ std::fclose(FUNIT(25)); FUNIT(25)=nullptr; }
    if(FUNIT(99)){ std::fclose(FUNIT(99)); FUNIT(99)=nullptr; }

    if (ierrflg==0){
        std::exit(0);
    } else if (ierrflg==-99){
        std::fprintf(stdout," program termination due to a preronnie error check\n");
        std::fprintf(stdout," see file preronnie.error\n");
        std::exit(1);
    } else if (ierrflg>-99){
        std::fprintf(stdout," program termination due to a ronnie error check\n");
        std::fprintf(stdout," see file ronnie.error\n");
        std::exit(1);
    } else { // ierrflg == -999
        std::fprintf(stdout," program termination due to receipt of a system signal\n");
        std::fprintf(stdout," see files preronnie.error and ronnie.error\n");
        std::exit(1);
    }
}
