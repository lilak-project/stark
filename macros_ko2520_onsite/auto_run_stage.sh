#!/usr/bin/env bash
set -eo pipefail
stage=${1:?Usage: auto_run_stage.sh STAGE_OR_MACRO.C [--dry-run]}; shift
case "${1:-}" in
  ''|--dry-run) ;;
  -h|--help) echo "RUN=8000 [START_RUN=97] ./auto_run_stage.sh $stage [--dry-run]"; exit 0 ;;
  *) exit 2 ;;
esac
(( $# <= 1 )) || exit 2
[[ -z "${RUN:-}" || "$RUN" =~ ^[0-9]+$ ]] || { echo 'Invalid RUN' >&2; exit 2; }

cd "$(dirname "$0")"
input_macro="get_${stage}.mac"
[[ "$stage" == *.C ]] && input_macro="$stage"
[[ -f "$input_macro" ]] || { echo "Missing input: $input_macro" >&2; exit 2; }
source ../../lilak.sh
export ROOTSYS=/root/opt/root-6.40.04
export PATH="$ROOTSYS/bin:/usr/bin:/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:$LILAK_PATH/build:/usr/local/gru/lib:/usr/local/get/lib"
plan_dir=$(mktemp -d "${TMPDIR:-/tmp}/auto_run.XXXXXX")
trap 'rm -rf -- "$plan_dir"' EXIT
export LILAK_ANALYSIS_PLAN="$plan_dir/runs.txt"
if ! root -l -b -q resolve_analysis_runs.C > "$plan_dir/log" 2>&1 || [[ ! -s "$LILAK_ANALYSIS_PLAN" ]]; then
  cat "$plan_dir/log" >&2
  exit 1
fi
mapfile -t plan < "$LILAK_ANALYSIS_PLAN"
runs=("${plan[@]:2}")
(( ${#runs[@]} > 0 )) || { echo 'No runs selected' >&2; exit 2; }

if [[ -n "${START_RUN:-}" ]]; then
  i=0
  while (( i < ${#runs[@]} )) && [[ "${runs[i]}" != "$START_RUN" ]]; do ((i+=1)); done
  (( i < ${#runs[@]} )) || { echo "START_RUN=$START_RUN is not in the list" >&2; exit 2; }
  runs=("${runs[@]:i}")
fi
macros=()
for run in "${runs[@]}"; do
  macro="$input_macro"
  if [[ "$stage" == reco ]] && compgen -G "data_raw/recovered/run_$(printf '%04d' "$run").*" >/dev/null; then
    # Replace InputPath instead of adding a second path that would read both datasets.
    macro="$plan_dir/get_reco_recovered.mac"
    if [[ ! -f "$macro" ]]; then
      awk '/^LKRun\/InputPath[ \t]/ { $0="LKRun/InputPath data_raw/recovered"; count++ }
           { print } END { if (count != 1) exit 2 }' get_reco.mac > "$macro"
    fi
  fi
  [[ -f "$macro" ]] || { echo "Missing configuration: $macro" >&2; exit 2; }
  macros+=("$macro")
done

printf 'Stage: %s | selection: %s | analysis: %s | runs: %d\n' "$stage" "${plan[0]}" "${plan[1]}" "${#runs[@]}"
for i in "${!runs[@]}"; do
  command=(lilak run "${macros[i]}")
  [[ "$stage" == *.C ]] && command=(root "${macros[i]}" -q -b)
  printf 'RUN=%s %s\n' "${runs[i]}" "${command[*]}"
  if [[ "${1:-}" != --dry-run ]]; then RUN="${runs[i]}" "${command[@]}"; fi
done
