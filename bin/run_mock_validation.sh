#!/bin/sh
# Realistic mock validation (objective 2, item 4): two truths built from the counts-based best fit,
# each with the published kinematic noise and Poisson counts (seeded), fitted with a free cored halo
# from the truth's own parameters plus one Latin start.
#   mockA: the fitted model itself (rho20 = 0)   -> does the fit find rho20 ~ 0?
#   mockB: the same stars with an injected halo    -> is the injected rho20 recovered?
# Usage: sh bin/run_mock_validation.sh [RHO20_INJECTED (default 1.0)] [TRUTH_JOB (default observed_no_halo_start0)]
set -eu
PY=/Users/vasilybelokurov/Work/venvs/.venv/bin/python
DRV=bin/run_compact_df_recovery.py
TAG=20260923
RHO=${1:-1.0}
JOB=${2:-observed_no_halo_start0}
TRUTH="results/df/dftwo_counts_${TAG}::${JOB}"
COMMON="--two-transition --branches free_halo --hand-starts 0 --n-starts 1 --start-from $TRUTH \
 --free-distance 4.8:6.0 --distance-prior 5.43:0.05 \
 --stop-delta 0.01 --jacobian-seed base --iteration-tolerance 5e-5 --probe-workers 5 \
 --bound stellar.J0=30:800 --bound stellar.J_a=5:500 --bound stellar.b_out=-2:2"
A=results/df/mockA_nohalo_${TAG}
B=results/df/mockB_rho${RHO}_${TAG}
[ -f $A/batch.json ] || $PY $DRV prepare-mock --out $A --truth "$TRUTH" --seed 101 $COMMON
[ -f $B/batch.json ] || $PY $DRV prepare-mock --out $B --truth "$TRUTH" --inject matter.rho20=$RHO --inject matter.r_s=30 --seed 202 $COMMON
for name in $A $B; do
    echo "$(date -u +%H:%M:%S) start $name"
    $PY $DRV run --out $name --workers 2 --max-calls 600 > $name/controller.log 2>&1
    echo "$(date -u +%H:%M:%S) done $name"
done
