#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CONFIG_FILE="${1:-config_Lingyu.cfg}"
if [[ "$CONFIG_FILE" != /* ]]; then
	CONFIG_FILE="$SCRIPT_DIR/$CONFIG_FILE"
fi

config_value() {
	local key="$1"
	awk -v key="$key" '
	{
		sub(/\/\/.*/, "", $0)
		if ($0 ~ "^[[:space:]]*" key "[[:space:]]*=") {
			value = $0
			sub(/^[^=]*=[[:space:]]*/, "", value)
			sub(/[[:space:]]*;.*/, "", value)
			gsub(/"/, "", value)
			gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
			print value
			exit
		}
	}' "$CONFIG_FILE"
}

normalize_number() {
	awk -v value="$1" 'BEGIN { printf "%.12g", value + 0 }' \
		| sed -e 's/e-0*/e-/' -e 's/e+0*/e/'
}

safe_name() {
	printf '%s' "$1" | sed -e 's/[[:space:]]//g' -e 's/+//g' -e 's/[^A-Za-z0-9._-]/_/g'
}

DM_MASS_RAW="$(config_value DM_mass)"
DM_CROSS_SECTION_RAW="$(config_value DM_cross_section_nucleon)"

if [[ -z "$DM_MASS_RAW" || -z "$DM_CROSS_SECTION_RAW" ]]; then
	echo "Could not read DM_mass or DM_cross_section_nucleon from $CONFIG_FILE" >&2
	exit 1
fi

DM_MASS="$(normalize_number "$DM_MASS_RAW")"
DM_CROSS_SECTION="$(normalize_number "$DM_CROSS_SECTION_RAW")"
RUN_TAG="$(safe_name "m${DM_MASS}GeV_sigmaN${DM_CROSS_SECTION}")"

echo "Submitting job_Lingyu.sh"
echo "Config: $CONFIG_FILE"
echo "Job name: $RUN_TAG"
echo "Output: job-${RUN_TAG}-%j"

sbatch --job-name="$RUN_TAG" --output="job-${RUN_TAG}-%j" job_Lingyu.sh "$CONFIG_FILE"