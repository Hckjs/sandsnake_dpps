#!/usr/bin/env bash
set -Eeuo pipefail

readonly BASE_DIR="/cephfs/projects/CTA/simulations/v0.23.1/PROD6/LaPalma"

readonly MCCAMPAIGN="PROD6"
readonly SITE="LaPalma"
readonly ARRAY_LAYOUT="Alpha"
readonly ANALYSIS_PROG="ctapipe-merge"
readonly MERGED="0"
readonly DATA_LEVEL="2"
readonly THETAP="20"
readonly PHIP="0"
readonly AZ="180"

readonly DOWNLOAD_JOBS=8

readonly FILE_DIR="dl${DATA_LEVEL}/zen_${THETAP}/az_${AZ}/merged_all"
readonly DIR="${BASE_DIR}/${FILE_DIR}"
readonly LOGS_DIR="${DIR}/logs"

readonly FILELIST="${DIR}/filelist.txt"
readonly LOGLIST="${LOGS_DIR}/loglist.txt"

find_data_lfns() {
  local output_file=$1

  dirac-dms-find-lfns \
    MCCampaign="$MCCAMPAIGN" \
    thetaP="$THETAP" \
    phiP="$PHIP" \
    site="$SITE" \
    array_layout="$ARRAY_LAYOUT" \
    analysis_prog="$ANALYSIS_PROG" \
    merged="$MERGED" \
    data_level="$DATA_LEVEL" \
    outputType=Data \
    >"$output_file"
}

find_log_lfns() {
  local output_file=$1

  dirac-dms-find-lfns \
    MCCampaign="$MCCAMPAIGN" \
    thetaP="$THETAP" \
    phiP="$PHIP" \
    site="$SITE" \
    array_layout="$ARRAY_LAYOUT" \
    outputType=Log \
    >"$output_file"
}

download_filelist() {
  local filelist=$1
  local workdir=$2
  local label=$3

  if [[ -s "$filelist" ]]; then
    (
      cd "$workdir" || exit 1
      echo "Downloading $label from $filelist"
      cta-prod-get-file -ddd -k "$DOWNLOAD_JOBS" "$filelist"
    )
  else
    echo "No $label found in $filelist"
  fi
}

mkdir -p "$DIR" "$LOGS_DIR"

find_data_lfns "$FILELIST"
find_log_lfns "$LOGLIST"

download_filelist "$FILELIST" "$DIR" "data files"
download_filelist "$LOGLIST" "$LOGS_DIR" "log files"

chmod -R ug+rwX,o-rwx "$DIR"

echo "Finished: $DIR"
echo
