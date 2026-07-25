#!/usr/bin/env bash
#SBATCH --job-name=damascus-sun
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.out

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
BUILD_DIR="${BUILD_DIR:-$PROJECT_ROOT/build}"
ORIGINAL_ARGS=("$@")
REBUILD=0

usage() {
	echo "Usage: submit.sh [--rebuild] [--build-dir DIR] <config.cfg>"
	echo "Run directly to submit, or use: sbatch [Slurm options] bin/submit.sh ..."
}

while [[ $# -gt 0 ]]; do
	case "$1" in
		--rebuild) REBUILD=1; shift ;;
		--build-dir) BUILD_DIR="${2:?--build-dir requires a path}"; shift 2 ;;
		--help|-h) usage; exit 0 ;;
		--*) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
		*) break ;;
	esac
done

[[ $# -eq 1 ]] || { usage >&2; exit 2; }
CONFIG_FILE="$1"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
	command -v sbatch >/dev/null 2>&1 || {
		echo "sbatch is unavailable; run this script on a Slurm login node." >&2
		exit 2
	}
	exec sbatch "$SCRIPT_DIR/submit.sh" "${ORIGINAL_ARGS[@]}"
fi

[[ "$BUILD_DIR" == /* ]] || BUILD_DIR="$PROJECT_ROOT/$BUILD_DIR"
[[ "$CONFIG_FILE" == /* ]] || CONFIG_FILE="${SLURM_SUBMIT_DIR:-$PWD}/$CONFIG_FILE"
[[ -f "$CONFIG_FILE" ]] || {
	echo "Configuration file not found: $CONFIG_FILE" >&2
	exit 2
}

COMPILE_ARGS=(--build-dir "$BUILD_DIR")
[[ "$REBUILD" -eq 0 ]] || COMPILE_ARGS+=(--rebuild)
"$SCRIPT_DIR/compile.sh" "${COMPILE_ARGS[@]}"

BINARY="$BUILD_DIR/src/DaMaSCUS-SUN"
[[ -x "$BINARY" ]] || {
	echo "Executable not found after build: $BINARY" >&2
	exit 2
}

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
exec srun "$BINARY" "$CONFIG_FILE"
