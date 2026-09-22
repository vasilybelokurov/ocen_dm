#!/bin/zsh
# Called inside the selected NEMO environment. All output remains in the run.
set -eu
ocen_compiler=${OCEN_NEMO_CXX:-g++-15}
ocen_precision=${OCEN_NEMO_PRECISION:-SINGLE}
case "$ocen_precision" in SINGLE|DOUBLE) ;; *) exit 2 ;; esac
exec "$ocen_compiler" -std=c++11 -O2 -fPIC -shared \
  -DfalcON_NEMO "-DfalcON_$ocen_precision" \
  -I"${0:A:h}/nbody_compat" \
  -I"$FALCON/inc" -I"$FALCON/inc/utils" -I"$NEMOINC" \
  "$1" -o "$2" \
  -L"$FALCON/lib" -lfalcON -L"$FALCON/utils/lib" -lWDutils \
  -L"$NEMOLIB" -lnemo
