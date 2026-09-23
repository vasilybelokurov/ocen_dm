#!/bin/sh
# Profile likelihood in rho20 for the two-transition DF with Poisson counts and the distance prior
# (docs/DF_ANISOTROPY_DIAGNOSTICS.md, objective 2). For each fixed rho20 all other coordinates are
# refitted (r_s free, 5-500 pc); starts: the counts-based no-halo best fit plus one Latin start.
# Batches run one after another, 2 jobs x (5 probe workers + 1) = 12 cores.
set -eu
PY=/Users/vasilybelokurov/Work/venvs/.venv/bin/python
DRV=bin/run_compact_df_recovery.py
TAG=20260923
SEED="results/df/dftwo_counts_${TAG}::observed_no_halo_start0"
COMMON="--two-transition --branches free_halo --hand-starts 0 --n-starts 1 --start-from $SEED \
 --counts ocen_counts_hst_f625w19,ocen_counts_gaia_g17 --free-distance 4.8:6.0 --distance-prior 5.43:0.05 \
 --stop-delta 0.01 --jacobian-seed base --iteration-tolerance 5e-5 --probe-workers 5 \
 --bound stellar.J0=30:800 --bound stellar.J_a=5:500 --bound stellar.b_out=-2:2"
for RHO in 0.25 0.5 1 2 4; do
    name=rho20scan_${RHO}_${TAG}
    [ -f results/df/$name/batch.json ] || $PY $DRV prepare-real --out results/df/$name $COMMON --fix matter.rho20=$RHO
    echo "$(date -u +%H:%M:%S) start $name"
    $PY $DRV run --out results/df/$name --workers 2 --max-calls 600 > results/df/$name/controller.log 2>&1
    echo "$(date -u +%H:%M:%S) done $name"
done
