#!/usr/bin/bash

# source diracos/diracosrc
# dirac-proxy-init
#
# tmux new -s download_mcs 'bash download_files_dl1_images.sh | tee -a ~/logs/download_files_dl1_images.log'
# tmux attach -t download_mcs
# Ctrl B, D
#
# For Updating crls:
# download 'igtf-preinstalled-bundle-classic.tar.gz' from
# https://dist.eugridpma.info/distribution/igtf/current/accredited/
# substitute into dirac/diracos/etc/grid-security/certificates/
# install fetch-crl
# fetch-crl --verbose --infodir dirac/diracos/etc/grid-security/certificates/

set -Eeuo pipefail
shopt -s nullglob

readonly BASE_DIR="/cephfs/projects/CTA/simulations/v0.23.1/PROD6/LaPalma"
readonly SPLIT_LINES=8
readonly DOWNLOAD_JOBS=8

readonly MCCAMPAIGN="PROD6"
readonly SITE="LaPalma"
readonly ARRAY_LAYOUT="Alpha"
readonly ANALYSIS_PROG="ctapipe-process"
readonly DATA_LEVEL="1"
readonly THETAP="20"
readonly PHIP="180"
readonly AZ="0"

PARTICLES=(
  "gamma"
  "gamma-diffuse"
  "electron"
  "proton"
)

readonly FILE_DIR="dl${DATA_LEVEL}/zen_${THETAP}/az_${AZ}"

find_data_lfns() {
  local output_file=$1
  local particle=$2

  dirac-dms-find-lfns \
    MCCampaign="$MCCAMPAIGN" \
    thetaP="$THETAP" \
    phiP="$PHIP" \
    site="$SITE" \
    array_layout="$ARRAY_LAYOUT" \
    analysis_prog="$ANALYSIS_PROG" \
    data_level="$DATA_LEVEL" \
    particle="$particle" \
    outputType=Data \
    >"$output_file"
}

find_log_lfns() {
  local output_file=$1
  local particle=$2

  dirac-dms-find-lfns \
    MCCampaign="$MCCAMPAIGN" \
    thetaP="$THETAP" \
    phiP="$PHIP" \
    site="$SITE" \
    array_layout="$ARRAY_LAYOUT" \
    data_level="$DATA_LEVEL" \
    particle="$particle" \
    outputType=Log \
    >"$output_file"
}

download_batches() {
  local batch_glob_dir=$1
  local workdir=$2
  local prefix=$3

  rm -f "$batch_glob_dir"/"${prefix}"_batch_*
  split -l "$SPLIT_LINES" "$batch_glob_dir"/"${prefix}"list.txt "$batch_glob_dir"/"${prefix}"_batch_

  (
    cd "$workdir" || exit 1
    for f in "$batch_glob_dir"/"${prefix}"_batch_*; do
      [[ -e "$f" ]] || continue
      echo "Processing file: $f"
      cta-prod-get-file -ddd -n -k "$DOWNLOAD_JOBS" "$f"
    done
  )
}

for particle in "${PARTICLES[@]}"; do
  echo "Particle: $particle"

  dir="${BASE_DIR}/${FILE_DIR}/${particle}"
  logs_dir="${dir}/logs"
  filelist_dir="${dir}/filelists"

  filelist="${filelist_dir}/filelist.txt"
  loglist="${logs_dir}/loglist.txt"

  mkdir -p "$filelist_dir" "$logs_dir"

  find_data_lfns "$filelist" "$particle"
  find_log_lfns "$loglist" "$particle"

  if [[ -s "$filelist" ]]; then
    download_batches "$filelist_dir" "$dir" "file"
  else
    echo "No files found in $filelist"
  fi

  if [[ -s "$loglist" ]]; then
    download_batches "$logs_dir" "$logs_dir" "log"
  else
    echo "No log files found in $loglist"
  fi

  chmod -R ug+rwX,o-rwx "$dir"

  echo "Finished: $dir"
  echo
done
