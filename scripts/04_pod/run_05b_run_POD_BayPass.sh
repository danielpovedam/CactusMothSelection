#!/usr/bin/env bash

set -euo pipefail


############################################################
## 49B_run_POD_BayPass.sh
##
## Run BayPass v3.1 on neutral POD data.
##
## Produces:
##   - XtXst null distribution
##   - C2 host-contrast null distribution
##
## Empirical Omega is supplied explicitly.
############################################################


BASE="/home/dpoveda/micetral/Population_snps/baypass"

PODDIR="${BASE}/POD_analysis/01_simulation"

RUNDIR="${BASE}/POD_analysis/02_baypass_run"

LOGDIR="${BASE}/POD_analysis/logs"

mkdir -p "${RUNDIR}" "${LOGDIR}"


POD_FILE="${PODDIR}/G.POD_neutral_100k"

OMEGA_FILE="${BASE}/runs_10x/all_snps/seed_12345/BP_all_snps_seed12345_mat_omega.out"

CONTRAST_FILE="${BASE}/covariates/Host_C2.cov"

PREFIX="${RUNDIR}/BP_POD_neutral100k_seed24680"


############################################################
## Locate BayPass
############################################################

if command -v baypass >/dev/null 2>&1; then

    BAYPASS_BIN="$(command -v baypass)"

else

    BAYPASS_BIN="$(
        find /home/dpoveda/micetral/selection/baypass \
        -type f \
        -name "baypass" \
        -perm -u+x \
        -print \
        -quit
    )"

fi


if [[ -z "${BAYPASS_BIN:-}" || ! -x "${BAYPASS_BIN}" ]]; then

    echo "ERROR: BayPass executable not found."
    echo
    echo "Locate it with:"
    echo
    echo "find /home/dpoveda -type f -name baypass -perm -u+x 2>/dev/null"
    exit 1

fi


echo "============================================================"
echo "49B — BAYPASS POD RUN"
echo "============================================================"
echo
echo "BayPass binary:"
echo "${BAYPASS_BIN}"
echo


############################################################
## Input QC
############################################################

for f in \
    "${POD_FILE}" \
    "${OMEGA_FILE}" \
    "${CONTRAST_FILE}"
do

    if [[ ! -s "${f}" ]]; then
        echo "ERROR: Missing/empty input: ${f}"
        exit 1
    fi

done


N_POD="$(wc -l < "${POD_FILE}")"

if [[ "${N_POD}" -ne 100000 ]]; then
    echo "ERROR: Expected 100000 POD SNPs; found ${N_POD}"
    exit 1
fi


N_COLS="$(awk 'NR==1{print NF}' "${POD_FILE}")"

if [[ "${N_COLS}" -ne 56 ]]; then
    echo "ERROR: Expected 56 POD columns; found ${N_COLS}"
    exit 1
fi


N_CONTRAST_COLS="$(awk 'NR==1{print NF}' "${CONTRAST_FILE}")"

if [[ "${N_CONTRAST_COLS}" -ne 28 ]]; then
    echo "ERROR: Host contrast does not contain 28 populations."
    exit 1
fi


echo "POD SNPs: ${N_POD}"
echo "POD columns: ${N_COLS}"
echo "Contrast columns: ${N_CONTRAST_COLS}"
echo


############################################################
## Run BayPass
##
## Important:
## - empirical Omega supplied through -omegafile
## - same host contrast as empirical C2 analysis
## - same MCMC configuration as empirical analysis
############################################################

"${BAYPASS_BIN}" \
    -countdatafile "${POD_FILE}" \
    -npop 28 \
    -omegafile "${OMEGA_FILE}" \
    -contrastfile "${CONTRAST_FILE}" \
    -npilot 20 \
    -pilotlength 500 \
    -burnin 10000 \
    -nval 2500 \
    -thin 20 \
    -seed 24680 \
    -nthreads 8 \
    -outprefix "${PREFIX}" \
    > "${LOGDIR}/49B_POD_BayPass.stdout.log" \
    2> "${LOGDIR}/49B_POD_BayPass.stderr.log"


############################################################
## Expected output QC
############################################################

XTX_FILE="${PREFIX}_summary_pi_xtx.out"

C2_FILE="${PREFIX}_summary_contrast.out"


if [[ ! -s "${XTX_FILE}" ]]; then
    echo "ERROR: XtX POD output was not created:"
    echo "${XTX_FILE}"
    exit 1
fi


if [[ ! -s "${C2_FILE}" ]]; then
    echo "ERROR: C2 POD output was not created:"
    echo "${C2_FILE}"
    exit 1
fi


N_XTX="$(( $(wc -l < "${XTX_FILE}") - 1 ))"

N_C2="$(( $(wc -l < "${C2_FILE}") - 1 ))"


echo
echo "============================================================"
echo "49B COMPLETE"
echo "============================================================"
echo
echo "XtX POD rows: ${N_XTX}"
echo "C2 POD rows:  ${N_C2}"
echo
echo "XtX file:"
echo "${XTX_FILE}"
echo
echo "C2 file:"
echo "${C2_FILE}"
echo
echo "DONE"
