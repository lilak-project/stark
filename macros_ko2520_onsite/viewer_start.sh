#!/usr/bin/env bash
cd "$(dirname "$0")" || exit
source ../../lilak.sh

# Pin the JSROOT server to the same ROOT release used by reconstruction and
# the 9093 Rutherford server.  An older ROOT on the inherited PATH can crash
# while serializing the newer canvases.
export ROOTSYS=/root/opt/root-6.40.04
export PATH="$ROOTSYS/bin:/usr/bin:/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:/root/lilak/build:/usr/local/gru/lib:/usr/local/get/lib"

lilak get_viewer -K -I 0.0.0.0:8877 
lilak get_viewer -K -I 0.0.0.0:8878 
lilak get_viewer -I 0.0.0.0:8877 -M mapping_ko2520_0904 -B data_raw
lilak get_viewer -I 0.0.0.0:8878 -M mapping_ko2520_0904 -B data_reco
