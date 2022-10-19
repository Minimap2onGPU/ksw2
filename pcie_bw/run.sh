#!/bin/bash

niter=20
nthread=4
async=1
pinned=1

#for pinned in {0,1}
#do
    nbuf=1
    for bufsize in {64,128,256,512,1024,2048,4096,8192,16000,32000,48000,60000}
    do
        cmd="OMP_NUM_THREADS=$nthread ./test $nbuf $bufsize $pinned $niter $async"
        eval $cmd
    done

    nbuf=16
    for bufsize in {1,2,4,8,10,12,14,16,32,48,64,128,256,512,1024,2048}
    do
        cmd="OMP_NUM_THREADS=$nthread ./test $nbuf $bufsize $pinned $niter $async"
        eval $cmd
    done
#done

