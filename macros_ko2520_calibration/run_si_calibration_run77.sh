#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

# libLILAK in this workspace was built against ROOT 6.40.04.
source /root/opt/root-6.40.04/bin/thisroot.sh

if ! command -v lilak >/dev/null 2>&1; then
    # The LILAK command is defined by this setup script as a shell function.
    set +u
    source /root/lilak/lilak.sh
    set -u
fi

start_stage="${START_STAGE:-0}"
if [[ ! "$start_stage" =~ ^[0-6]$ ]]; then
    echo "START_STAGE must be an integer from 0 through 6" >&2
    exit 2
fi

for ((stage = start_stage; stage <= 6; stage++)); do
    echo "[run 77 calibration] starting stage ${stage}"
    set +e
    STAGE="$stage" lilak run si_calibration_run77.mac
    status=$?
    set -e

    # ROOT 6.40 can segfault while destroying canvases after LILAK has already
    # completed Terminate and saved them.  Keep advancing only for that exact
    # post-run signal; all other failures remain fatal.
    if ((status != 0 && status != 11 && status != 139)); then
        echo "[run 77 calibration] stage ${stage} failed with status ${status}" >&2
        exit "$status"
    fi
    if ((status == 11 || status == 139)); then
        echo "[run 77 calibration] stage ${stage}: ignoring ROOT post-terminate SIGSEGV (status ${status})"
    fi
    echo "[run 77 calibration] completed stage ${stage}"
done

echo "[run 77 calibration] all stages completed"
