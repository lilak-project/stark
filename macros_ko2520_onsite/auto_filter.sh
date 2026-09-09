#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
./auto_run_stage.sh filter "$@"
[[ "${1:-}" == --dry-run ]] && { echo "RUN=${RUN:-<default_run>} root draw_sihit.C -q -b"; exit; }
source ../../lilak.sh
export ROOTSYS=/root/opt/root-6.40.04
export PATH="$ROOTSYS/bin:/usr/bin:/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:$LILAK_PATH/build:/usr/local/gru/lib:/usr/local/get/lib"
#root draw_sihit.C -q -b
