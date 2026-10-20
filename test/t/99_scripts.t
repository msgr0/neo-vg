#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 2

nvg construct -a -r small/x.fa -v small/x.vcf.gz >x.vg
nvg gbwt -o x-haps.gbwt -x x.vg -v small/x.vcf.gz
nvg gbwt -o x-paths.gbwt -x x.vg --index-paths
nvg gbwt -o x-merged.gbwt -m x-haps.gbwt x-paths.gbwt
nvg gbwt -g x.gbz --gbz-format --augment-gbwt -x x.vg x-merged.gbwt
nvg index -j x.dist x.vg
nvg minimizer -k 29 -w 11 -d x.dist -o x.shortread.withzip.min -z x.shortread.zipcodes x.gbz
nvg giraffe --track-provenance -Z x.gbz -m x.shortread.withzip.min -z x.shortread.zipcodes -d x.dist -f reads/small.middle.ref.fq > mapped1.gam

../scripts/giraffe-facts.py mapped1.gam facts-out >facts.txt

is "$?" "0" "giraffe-facts.py processes Giraffe output successfully"
is "$(grep cluster-coverage facts.txt | wc -l)" "1" "giraffe-facts.py lists filters used when provenance is tracked"
rm -f x.vg x-haps.gbwt x-paths.gbwt x-merged.gbwt x.gbz x.giraffe.gbz x.dist x.shortread.withzip.min x.shortread.zipcodes facts.txt mapped1.gam
rm -Rf facts-out



