#!/usr/bin/env bash

# Add the Rutherford-fit canvas whenever a completed compass reco file arrives.
# A per-file stamp prevents repeated processing, and flock prevents duplicate
# watchers when start_js_3.sh is run more than once.

set -u

watch_dir=${1:-data_reco3}
macro=${2:-${watch_dir}/rutherford_beam.C}
publish_dir=${3:-}
sanitizer=${4:-}
poll_seconds=${RUTHERFORD_POLL_SECONDS:-5}
state_dir="${watch_dir}/.rutherford_state"
lock_file="${state_dir}/watch.lock"
root_bin=/root/opt/root-6.40.04/bin/root
root_library_path=/root/opt/root-6.40.04/lib:/root/lilak/build:/usr/local/gru/lib:/usr/local/get/lib

mkdir -p "$state_dir"
[ -z "$publish_dir" ] || mkdir -p "$publish_dir"
exec 9>"$lock_file"
if ! flock -n 9; then
  exit 0
fi

declare -A previous_stamp

while true; do
  for reco_file in "$watch_dir"/compass_0*.reco.root; do
    [ -f "$reco_file" ] || continue

    stamp=$(stat -c '%s:%Y' "$reco_file" 2>/dev/null) || continue
    if [ "${previous_stamp[$reco_file]:-}" != "$stamp" ]; then
      previous_stamp[$reco_file]=$stamp
      continue
    fi

    base_name=${reco_file##*/}
    stamp_file="${state_dir}/${base_name}.stamp"
    if [ -f "$stamp_file" ] && [ "$(<"$stamp_file")" = "$stamp" ]; then
      continue
    fi

    printf 'LILAK_JS_STATUS Rutherford fit: %s\n' "$base_name"
    root_output=$(ROOTSYS=/root/opt/root-6.40.04 \
      LD_LIBRARY_PATH="$root_library_path" "$root_bin" -l -b -q \
      "${macro}(\"${reco_file}\",8.0,2,185.,-1.,true)" 2>&1)
    root_status=$?
    printf '%s\n' "$root_output"

    if [ "$root_status" -eq 0 ] && grep -q 'canvas written to' <<<"$root_output"; then
      run_stem=${base_name%.reco.root}
      sidecar="${watch_dir}/${run_stem}.rutherford.root"
      run_number=${run_stem#compass_}

      # Starting with run 0447, append fixed five-minute Rutherford fits and
      # a five-minute beam-rate trend.  Stop JSROOT while the sidecar is being
      # updated, then restart it after the file is complete.
      if (( 10#$run_number >= 447 )); then
        server_pid_file=/tmp/rutherford_js_9093.pid
        if [ -r "$server_pid_file" ]; then
          read -r server_pid <"$server_pid_file"
          if [[ "$server_pid" =~ ^[0-9]+$ ]] && kill -0 "$server_pid" 2>/dev/null; then
            kill "$server_pid"
          fi
        fi

        window_output=$(ROOTSYS=/root/opt/root-6.40.04 \
          LD_LIBRARY_PATH="$root_library_path" "$root_bin" -l -b -q \
          -e ".L ${macro}" \
          -e "rutherford_beam_windows(\"${reco_file}\",5.0,0,8.0,2,185.0,-1.0)" 2>&1)
        window_status=$?
        printf '%s\n' "$window_output"
        /root/lilak/stark/macros_ko2520_onsite/start_js_3.sh
        if [ "$window_status" -ne 0 ] || ! grep -q 'beam-rate trend written' <<<"$window_output"; then
          printf 'Five-minute Rutherford fits failed; will retry: %s\n' "$base_name" >&2
          continue
        fi
      fi

      if [ -n "$publish_dir" ] && [ -f "$sidecar" ]; then
        # Publish atomically so JSROOT never opens a half-copied file.
        if [ -n "$sanitizer" ]; then
          ROOTSYS=/root/opt/root-6.40.04 \
            LD_LIBRARY_PATH="$root_library_path" "$root_bin" -n -l -b -q \
            "${sanitizer}(\"${sidecar}\",\"${publish_dir}/.${base_name}.tmp\")"
        else
          cp "$sidecar" "${publish_dir}/.${base_name}.tmp"
        fi
        mv "${publish_dir}/.${base_name}.tmp" "${publish_dir}/${base_name}"
      fi
      # The macro changes the ROOT file, so record its post-write stamp.
      stamp=$(stat -c '%s:%Y' "$reco_file" 2>/dev/null) || continue
      printf '%s\n' "$stamp" >"${stamp_file}.tmp"
      mv "${stamp_file}.tmp" "$stamp_file"
      previous_stamp[$reco_file]=$stamp
      printf 'Rutherford fit complete: %s\n' "$base_name"
    else
      printf 'Rutherford fit failed; will retry: %s\n' "$base_name" >&2
    fi
  done
  sleep "$poll_seconds"
done
