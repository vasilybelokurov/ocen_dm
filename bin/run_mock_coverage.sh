#!/bin/sh
# Coverage test of the rho20 recovery: N noise seeds for two truths built from the counts-based
# best fit (A: rho20 = 0; B: injected rho20 = RHO, r_s = 30 pc), each fitted with a free halo from one
# generic hand-picked start (never from the truth). A and B of each seed run concurrently
# (2 x (5 probe workers + 1) = 12 cores).
# Usage: sh bin/run_mock_coverage.sh [RHO=1.0] [SEEDS="11 12 13 14 15"]
set -eu
PY=/Users/vasilybelokurov/Work/venvs/.venv/bin/python
DRV=bin/run_compact_df_recovery.py
TAG=20260923
RHO=${1:-1.0}
SEEDS=${2:-"11 12 13 14 15"}
TRUTH="results/df/dftwo_counts_${TAG}::observed_no_halo_start0"
COMMON="--two-transition --branches free_halo --hand-starts 1 --n-starts 0 \
 --free-distance 4.8:6.0 --distance-prior 5.43:0.05 \
 --stop-delta 0.01 --jacobian-seed base --iteration-tolerance 5e-5 --probe-workers 5 \
 --bound stellar.J0=30:800 --bound stellar.J_a=5:500 --bound stellar.b_out=-2:2"
for S in $SEEDS; do
    A=results/df/cov_A_s${S}_${TAG}; B=results/df/cov_B_rho${RHO}_s${S}_${TAG}
    [ -f $A/batch.json ] || $PY $DRV prepare-mock --out $A --truth "$TRUTH" --seed $S $COMMON
    [ -f $B/batch.json ] || $PY $DRV prepare-mock --out $B --truth "$TRUTH" --inject matter.rho20=$RHO --inject matter.r_s=30 --seed $((S+1000)) $COMMON
    echo "$(date -u +%H:%M:%S) start seed $S"
    $PY $DRV run --out $A --workers 1 --max-calls 600 > $A/controller.log 2>&1 &
    $PY $DRV run --out $B --workers 1 --max-calls 600 > $B/controller.log 2>&1 &
    wait
    echo "$(date -u +%H:%M:%S) done seed $S"
done
