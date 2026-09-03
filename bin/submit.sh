#!/usr/bin/env bash
#SBATCH --job-name=damascus-sun
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --output=/lustre/project/kennyng/damascus-earth/DaMaSCUS-NewEarth-EVAP/bin/job-%x-%j.out
#SBATCH --error=/lustre/project/kennyng/damascus-earth/DaMaSCUS-NewEarth-EVAP/bin/job-%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=1155234222@link.cuhk.edu.hk

set -Eeuo pipefail

# 硬编码项目根目录，不依赖 git / BASH_SOURCE（SLURM spool 会破坏它们）
PROJECT_ROOT="${PROJECT_ROOT:-/lustre/project/kennyng/damascus-earth/DaMaSCUS-NewEarth-EVAP}"
export PROJECT_ROOT
SCRIPT_DIR="$PROJECT_ROOT/bin"
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

# ---- 进入 SLURM 作业后：加载编译/运行环境 ----
module purge
module load gcc/10.2.0
module load boost/1.77.0
module load openmpi/gcc/4.0.5
module load mpfr/3.1.4
export PATH="$HOME/cmake-3.12.4-Linux-x86_64/bin:$PATH"
export CMAKE_BIN="$HOME/cmake-3.12.4-Linux-x86_64/bin/cmake"
export LIBCONFIG_PREFIX="$HOME/libconfig-1.7.3"
export LD_LIBRARY_PATH="$HOME/libconfig-1.7.3/lib:${LD_LIBRARY_PATH:-}"

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
exec mpirun --mca btl ^openib "$BINARY" "$CONFIG_FILE"
