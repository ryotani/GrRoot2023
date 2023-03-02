#!/bin/bash

export GRROOT_DIR=/shared/storage/physnp/NSCLData/GrROOT6/GrROOT

if [ `file /bin/bash | awk '{print $3}'` = "64-bit" ]; then
    export GRROOT_BINDIR=$GRROOT_DIR/bin/bin64
    export GRROOT_LIBDIR=$GRROOT_DIR/lib/lib64
else
    export GRROOT_BINDIR=$GRROOT_DIR/bin/bin32
    export GRROOT_LIBDIR=$GRROOT_DIR/lib/lib32
fi
export PATH=$PATH:$GRROOT_BINDIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GRROOT_LIBDIR
