# Source this to run gyrfalcON (NEMO build in ~/Work/src/nemo_satellite_experiment) with the
# AGAMA external-potential plug-in on this Mac. Findings of 2026-09-20 (docs/NBODY_CODE_ASSESSMENT.md):
#  * the falcON binaries carry relative install names, so both library dirs must be on DYLD_LIBRARY_PATH;
#  * never launch through /usr/bin/time or other SIP-protected binaries (DYLD_* is stripped);
#  * agama.so links LLVM libomp, whose initialisation leaves a stale dlerror() string that NEMO's
#    loadobj misreads as a failed dlopen -> preload a stub defining the looked-up symbols.
source ~/Work/src/nemo_satellite_experiment/nemo_start.sh > /dev/null 2>&1
export DYLD_LIBRARY_PATH=$FALCON/lib:$FALCON/utils/lib:$NEMOLIB
_stub_dir=${OCEN_NBODY_SCRATCH:-$HOME/.cache/ocen_dm}; mkdir -p "$_stub_dir"
if [ ! -f "$_stub_dir/libompstub.dylib" ]; then
  printf 'void llvm_omp_target_unlock_mem(void){}\nvoid llvm_omp_target_lock_mem(void){}\nvoid llvm_omp_target_alloc_host(void){}\nvoid llvm_omp_target_free_host(void){}\n' > "$_stub_dir/ompstub.c"
  cc -dynamiclib -o "$_stub_dir/libompstub.dylib" "$_stub_dir/ompstub.c"
fi
export DYLD_INSERT_LIBRARIES="$_stub_dir/libompstub.dylib"
[ -e "$NEMOOBJ/acc/agama.so" ] || ln -s "$HOME/Work/venvs/.venv/lib/python3.13/site-packages/agama/agama.so" "$NEMOOBJ/acc/agama.so"
# usage: gyrfalcON in=ic.nemo out=run.nemo eps=... kmax=... Grav=4.30091e-6 accname=agama accpars=0,4.30091e-6 accfile=pot.ini
