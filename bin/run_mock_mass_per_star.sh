#!/bin/sh
# Mock test of the pinned stellar mass per counted star: the coverage mocks (same seeds, hence identical
# mock data to results/df/cov_*), refitted with the HST count amplitude pinned at M/N19 = M. M = truth
# pins the mock's own value (10.67 Msun per F625W<19 star): does the rho20 recovery scatter shrink?
# A wrong M measures the bias a wrong external stellar-mass assumption induces in rho20.
# Truth A: rho20 = 0; truth B: injected rho20 = RHO, r_s = 30 pc. One generic hand-picked start per
# mock (never the truth), A and B of a seed run concurrently (12 cores).
# Usage: sh bin/run_mock_mass_per_star.sh [M=truth] [SEEDS="11 12 13 14 15"] [RHO=1.0]
set -eu
PY=/Users/vasilybelokurov/Work/venvs/.venv/bin/python
DRV=bin/run_compact_df_recovery.py
TAG=20260923
M=${1:-truth}
SEEDS=${2:-"11 12 13 14 15"}
RHO=${3:-1.0}
TRUTH="results/df/dftwo_counts_${TAG}::observed_no_halo_start0"
COMMON="--two-transition --branches free_halo --hand-starts 1 --n-starts 0 \
 --free-distance 4.8:6.0 --distance-prior 5.43:0.05 \
 --stop-delta 0.01 --jacobian-seed base --iteration-tolerance 5e-5 --probe-workers 5 \
 --bound stellar.J0=30:800 --bound stellar.J_a=5:500 --bound stellar.b_out=-2:2 \
 --mass-per-star ocen_counts_hst_f625w19=$M"
for S in $SEEDS; do
    A=results/df/mpsmock_${M}_A_s${S}_${TAG}; B=results/df/mpsmock_${M}_B_rho${RHO}_s${S}_${TAG}
    [ -f $A/batch.json ] || $PY $DRV prepare-mock --out $A --truth "$TRUTH" --seed $S $COMMON
    [ -f $B/batch.json ] || $PY $DRV prepare-mock --out $B --truth "$TRUTH" --inject matter.rho20=$RHO --inject matter.r_s=30 --seed $((S+1000)) $COMMON
    echo "$(date -u +%H:%M:%S) start M=$M seed $S"
    $PY $DRV run --out $A --workers 1 --max-calls 600 > $A/controller.log 2>&1 &
    $PY $DRV run --out $B --workers 1 --max-calls 600 > $B/controller.log 2>&1 &
    wait
    echo "$(date -u +%H:%M:%S) done M=$M seed $S"
done
