#!/bin/sh
# rho20 as a function of the assumed stellar mass per counted star (objective 2, external stellar-mass
# constraint). For each value of M/N19 (Msun of mass-that-follows-light per star with F625W < 19, the
# HST count product's selection) the HST count amplitude is pinned, so the counts fix the stellar mass
# normalization, and everything else (DF shape, remnants, cored halo rho20 + r_s, distance with prior)
# is refitted. The photometric mass function gives M/N19 = 8.4-9.5 in the zones dominating the
# kinematics (JOURNAL 2026-09-23); the unconstrained counts fit returns 10.67. Starts: the counts-based
# no-halo best fit plus one Latin start; 2 jobs x (5 probe workers + 1) = 12 cores per batch.
# Usage: sh bin/run_mass_per_star_scan.sh [M/N19 values, default "7 8 9 10 11 12"]
set -eu
PY=/Users/vasilybelokurov/Work/venvs/.venv/bin/python
DRV=bin/run_compact_df_recovery.py
TAG=20260923
SEED="results/df/dftwo_counts_${TAG}::observed_no_halo_start0"
COMMON="--two-transition --branches free_halo --hand-starts 0 --n-starts 1 --start-from $SEED \
 --counts ocen_counts_hst_f625w19,ocen_counts_gaia_g17 --free-distance 4.8:6.0 --distance-prior 5.43:0.05 \
 --stop-delta 0.01 --jacobian-seed base --iteration-tolerance 5e-5 --probe-workers 5 \
 --bound stellar.J0=30:800 --bound stellar.J_a=5:500 --bound stellar.b_out=-2:2"
for M in ${*:-7 8 9 10 11 12}; do
    name=mpsscan_${M}_${TAG}
    [ -f results/df/$name/batch.json ] || $PY $DRV prepare-real --out results/df/$name $COMMON \
        --mass-per-star ocen_counts_hst_f625w19=$M
    echo "$(date -u +%H:%M:%S) start $name"
    $PY $DRV run --out results/df/$name --workers 2 --max-calls 600 > results/df/$name/controller.log 2>&1
    echo "$(date -u +%H:%M:%S) done $name"
done
