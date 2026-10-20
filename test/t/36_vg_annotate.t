#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 9

nvg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >t.vg

nvg mod -N t.vg >t.ref.vg
nvg index -x t.xg t.vg
nvg index -x t.ref.xg t.ref.vg

is "$(nvg annotate -n -x t.ref.xg -a tiny/tiny-s7331-n10-l50.gam | awk '{ if ($5 < 50) print }' | wc -l)" "10" "we can detect when reads contain non-reference variation"

nvg annotate -b tiny/tiny.bed -x t.ref.xg -a tiny/tiny-s543-n30-l10.gam > annotated.gam
is "$(nvg view -aj annotated.gam | jq -c '.annotation.features' | grep feat1 | wc -l)" 3 "nvg annotate finds the right number of reads overlapping a feature"
is "$(nvg view -aj annotated.gam | grep feat1 | grep -e '"node_id": *"1"' | wc -l)" "$(nvg view -aj annotated.gam | grep feat1 | wc -l)" "all reads overlapping a feature fall on its node"
is "$(nvg view -aj annotated.gam | jq -c '.annotation.features' | grep feat1 | grep feat2 | wc -l)" 0 "nvg annotate finds no reads touching both of two distant features"
is "$(nvg view -aj annotated.gam | jq -c '.annotation.features' | grep feat2 | grep feat3 | wc -l)" 2 "nvg annotate shows reads having to go through one feature to get to another at the end"
is "$(nvg view -aj annotated.gam | jq -c '.annotation.features' | grep featAll | wc -l)" 30 "nvg annotate shows all reads overlapping a whole-reference-covering feature"

rm -f t.vg t.ref.vg t.xg t.ref.xg annotated.gam

nvg construct -r small/x.fa -v small/x.vcf.gz >x.vg
nvg index -x x.xg x.vg

nvg annotate -p -x x.xg -a small/x-s1337-n1.gam > annot.gam

is "$(nvg annotate -p -x x.xg -a small/x-s1337-n1.gam | nvg view -aj - | jq -c '.refpos[]' | wc -l)" "1" "annotating with earliest path position produces one position"
is "$(nvg annotate -m -x x.xg -a small/x-s1337-n1.gam | nvg view -aj - | jq -c '.refpos[]' | wc -l)" "13" "annotating with multiple path positions produces multiple positions"

rm -f x.xg x.vg annot.gam

nvg view -Jv cyclic/circular_path.json > circular_path.vg
nvg index -x circular_path.xg circular_path.vg
nvg annotate -p -x circular_path.xg -b cyclic/circular_path_origin.bed > circular_path_origin.gam

is "$(nvg view -aj circular_path_origin.gam | jq -c '.path.mapping[].position')" "$(printf '{"node_id":"1","offset":"5"}\n{"node_id":"1"}')" "annotations derived from BED can span the end/start joins of circular paths"

rm -f circular_path.vg circular_path.xg circular_path_origin.gam
