#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

nvg index -x hla.xg msgas/hla_v.vg
nvg dotplot -x hla.xg >/dev/null

is "${?}" "0" "nvg dotplot runs successfully"

rm -f hla.xg


