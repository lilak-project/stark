#!/usr/bin/env bash
set -eo pipefail

cd "$(dirname "$0")" || exit
source ../../lilak.sh
set -u

# Pin the JSROOT server to the same ROOT release used by reconstruction and
# the 9093 Rutherford server.  An older ROOT on the inherited PATH can crash
# while serializing the newer canvases.
export ROOTSYS=/root/opt/root-6.40.04
export PATH="$ROOTSYS/bin:/usr/bin:/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:/root/lilak/build:/usr/local/gru/lib:/usr/local/get/lib"

expected_root="$ROOTSYS/bin/root"
actual_root="$(command -v root)"
actual_root_version="$(root-config --version)"
if [[ "$actual_root" != "$expected_root" || "$actual_root_version" != "6.40.04" ]]; then
  printf 'ROOT environment mismatch: expected %s (6.40.04), got %s (%s)\n' \
    "$expected_root" "$actual_root" "$actual_root_version" >&2
  exit 1
fi
printf 'Using ROOT %s from %s\n' "$actual_root_version" "$actual_root"

production_runs=(
  83 84 85 86 87 88 89 90
  97 98 99 100 101 102 105
  112 113 114
  134 135 136 139 141 142 146 147 148
)

start_run="${START_RUN:-${1:-${production_runs[0]}}}"
start_found=false
processed_count=0
for run in "${production_runs[@]}"; do
  if [[ "$run" == "$start_run" ]]; then
    start_found=true
  fi
  if [[ "$start_found" != true ]]; then
    continue
  fi
  ((processed_count += 1))

  reco_macro=get_reco.mac
  recovered_pattern="data_raw/recovered/run_$(printf '%04d' "$run").*"
  if compgen -G "$recovered_pattern" >/dev/null; then
    reco_macro=get_reco_recovered.mac
    printf '[%s] Using sanitized raw input for production run %04d\n' \
      "$(date -u +%FT%TZ)" "$run"
  fi

  #printf '\n[%s] Reconstructing production run %04d with mapping_ko2520_0904\n' \
  #  "$(date -u +%FT%TZ)" "$run"
  #RUN="$run" lilak run "$reco_macro"

  printf '\n[%s] Building SiHit for production run %04d\n' \
    "$(date -u +%FT%TZ)" "$run"
  RUN="$run" lilak run get_sihit.mac
done

if [[ "$start_found" != true ]]; then
  printf 'START_RUN=%s is not in the production run list.\n' "$start_run" >&2
  exit 2
fi

printf '\n[%s] Finished all %d production runs.\n' \
  "$(date -u +%FT%TZ)" "$processed_count"
