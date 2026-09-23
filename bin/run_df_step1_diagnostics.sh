#!/bin/sh
# Step-1 diagnostics for the observed-data DF failure (docs/OBSERVED_DF_FLEXIBILITY.md).
# Current DF, halo frozen (no-halo branch), two hand-picked starts, wide bounds.
# Batches run one after another; each uses 2 jobs x (5 probe workers + 1) = 12 cores.
set -eu
PY=/Users/vasilybelokurov/Work/venvs/.venv/bin/python
DRV=bin/run_compact_df_recovery.py
TAG=20260923
COMMON="--stop-delta 0.01 --jacobian-seed base --iteration-tolerance 5e-5 --probe-workers 5 \
 --bound stellar.J0=30:800 --bound stellar.J_a=5:3000 --bound matter.r_s=5:500 --branches no_halo"

prepare() {  # name, extra options
    name=$1; shift
    [ -f results/df/$name/batch.json ] || $PY $DRV prepare-real --out results/df/$name $COMMON "$@"
}
prepare dfdiag_hst_muse_$TAG --datasets hst_pm_radial_ours,hst_pm_tangential_ours,muse_los_dispersion
prepare dfdiag_gaia_$TAG --datasets gaia_edr3_ours_radial,gaia_edr3_ours_tangential
prepare dfdiag_gaiafloor01_$TAG --gaia-error-floor 0.01
prepare dfdiag_gaiafloor02_$TAG --gaia-error-floor 0.02
prepare dfdiag_freedist_$TAG --free-distance 4.8:6.0

for name in dfdiag_hst_muse_$TAG dfdiag_gaia_$TAG dfdiag_gaiafloor01_$TAG dfdiag_gaiafloor02_$TAG dfdiag_freedist_$TAG; do
    echo "$(date -u +%H:%M:%S) start $name"
    $PY $DRV run --out results/df/$name --workers 2 --max-calls 600 > results/df/$name/controller.log 2>&1
    echo "$(date -u +%H:%M:%S) done $name"
done
