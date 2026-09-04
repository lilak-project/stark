#!/usr/bin/env bash
cd "$(dirname "$0")" || exit
source ../../lilak.sh

# Keep one background watcher alive.  It waits until an incoming reco file is
# stable and writes the complete fit canvas into rutherford_sidecars/.
mkdir -p data_reco3/rutherford_sidecars
setsid -f ./watch_rutherford.sh data_reco3 data_reco3/rutherford_beam.C \
  >>data_reco3/.rutherford_watch.log 2>&1 </dev/null

# Forward only reconstruction files created after the bridge's first launch.
# Its persistent baseline keeps every older compass reco file excluded.
compass_bridge=/root/lilak/stark/macros_ko2520_compass/watch_new_reco_for_9093.sh
compass_reco=/root/lilak/stark/macros_ko2520_compass/data_reco
compass_bridge_log=/root/lilak/stark/macros_ko2520_compass/.reco_9093_bridge/watch.log
mkdir -p "$(dirname "$compass_bridge_log")"
setsid -f "$compass_bridge" "$compass_reco" "$(readlink -f data_reco3)" \
  >>"$compass_bridge_log" 2>&1 </dev/null

# The generic directory timer hangs in the freshly started ROOT version used
# here.  This dedicated server loads each original reco and Rutherford file
# once, keeps legacy canvases alive, and exposes the requested hierarchy.
server_pid_file=/tmp/rutherford_js_9093.pid
if [ -r "$server_pid_file" ]; then
  read -r old_server_pid <"$server_pid_file"
  if [[ "$old_server_pid" =~ ^[0-9]+$ ]] && kill -0 "$old_server_pid" 2>/dev/null; then
    kill "$old_server_pid"
  fi
fi

server_directory=$(readlink -f data_reco3)
fit_directory=$(readlink -f data_reco3/rutherford_sidecars)
server_macro=$(readlink -f ./serve_rutherford_9093.C)
server_log=/tmp/rutherford_js_9093.log
root_bin=/root/opt/root-6.40.04/bin/root
root_library_path=/root/opt/root-6.40.04/lib:/root/lilak/build:/usr/local/gru/lib:/usr/local/get/lib

setsid -f env \
  RUTHERFORD_JS_ADDRESS=192.168.1.102:9093 \
  RUTHERFORD_JS_DIRECTORY="$server_directory" \
  RUTHERFORD_JS_FIT_DIRECTORY="$fit_directory" \
  RUTHERFORD_JS_PIDFILE="$server_pid_file" \
  ROOTSYS=/root/opt/root-6.40.04 \
  LD_LIBRARY_PATH="$root_library_path" \
  "$root_bin" -n -l -b -q "$server_macro" \
  >>"$server_log" 2>&1 </dev/null

echo "Rutherford JSROOT URL: http://192.168.1.102:9093"
echo "Watched directory: $server_directory"
