#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 18

nvg construct -r small/x.fa -v small/x.vcf.gz >x.vg
nvg index -x x.xg  x.vg
nvg sim -n 1000 -l 100 -e 0.01 -i 0.005 -x x.xg -a >x.gam


# GAM
nvg gamsort x.gam >x.sorted.gam

nvg view -aj x.sorted.gam | jq -r '.path.mapping | ([.[] | .position.node_id | tonumber] | min)' >min_ids.gamsorted.txt
nvg view -aj x.gam | jq -r '.path.mapping | ([.[] | .position.node_id | tonumber] | min)' | sort -n >min_ids.sorted.txt

is "$(md5sum <min_ids.gamsorted.txt)" "$(md5sum <min_ids.sorted.txt)" "Sorting a GAM orders the alignments by min node ID"

nvg gamsort x.gam -i x.sorted.gam.gai >x.sorted.gam
is "$?" "0" "sorted GAMs can be indexed during the sort"

nvg gamsort --shuffle x.sorted.gam >x.shuffled.gam
is "$?" "0" "GAMs can be shuffled"
is "$(nvg stats -a x.shuffled.gam)" "$(nvg stats -a x.sorted.gam)" "Shuffling preserves read data"

rm -f x.sorted.gam x.sorted.2.gam x.shuffled.gam min_ids.gamsorted.txt min_ids.sorted.txt x.sorted.gam.gai x.sorted.2.gam.gai


# GAF. Correctness is tested in unit tests, so we just check that the commands work.
nvg convert -G x.gam x.xg > x.gaf
grep "^@" x.gaf > x.gaf.header
grep -v "^@" x.gaf > x.gaf.no_header
sort x.gaf.no_header > x.gaf.lexicographic
rm -f x.gbwt

nvg gamsort -G x.gaf > x.sorted.gaf
is "$?" "0" "GAFs can be sorted"
is "$(ls x.gbwt 2> /dev/null)" "" "GBWT index is not created by default"
grep "^@" x.gaf > x.sorted.gaf.header
cmp x.gaf.header x.sorted.gaf.header
is "$?" "0" "Sorting a GAF preserves header"
grep -v "^@" x.sorted.gaf | sort > x.sorted.gaf.lexicographic
cmp x.gaf.lexicographic x.sorted.gaf.lexicographic
is "$?" "0" "Sorting a GAF preserves read data"

nvg gamsort -G x.gaf.no_header > x.sorted.gaf
is "$?" "0" "GAFs without headers can be sorted"
sort x.sorted.gaf > x.sorted.gaf.lexicographic
cmp x.gaf.lexicographic x.sorted.gaf.lexicographic
is "$?" "0" "Sorting a GAF without header preserves read data"

nvg gamsort -G x.gaf -g x.gbwt > x.sorted.gaf
is "$?" "0" "GAFs can be sorted while creating a GBWT index"
is "$(ls x.gbwt 2> /dev/null)" "x.gbwt" "GBWT index can be created"
grep -v "^@" x.sorted.gaf | sort > x.sorted.gaf.lexicographic
cmp x.gaf.lexicographic x.sorted.gaf.lexicographic
is "$?" "0" "Sorting a GAF preserves read data with GBWT index"

nvg gamsort -G x.gaf -g x.gbwt --bidirectional > x.sorted.gaf
is "$?" "0" "GAFs can be sorted while creating a bidirectional GBWT index"
is "$(ls x.gbwt 2> /dev/null)" "x.gbwt" "Bidirectional GBWT index can be created"
grep -v "^@" x.sorted.gaf | sort > x.sorted.gaf.lexicographic
cmp x.gaf.lexicographic x.sorted.gaf.lexicographic
is "$?" "0" "Sorting a GAF preserves read data with bidirectional GBWT index"

nvg gamsort -G --shuffle x.gaf > x.shuffled.gaf
is "$?" "0" "GAFs can be shuffled"
grep -v "^@" x.shuffled.gaf | sort > x.shuffled.gaf.lexicographic
cmp x.gaf.lexicographic x.shuffled.gaf.lexicographic
is "$?" "0" "Shuffling a GAF preserves read data"

rm -f x.gaf x.gaf.header x.gaf.no_header x.gaf.lexicographic
rm -f x.sorted.gaf x.sorted.gaf.header x.sorted.gaf.lexicographic
rm -f x.shuffled.gaf x.shuffled.gaf.lexicographic
rm -f x.gbwt


# Cleanup
rm -f x.vg x.xg x.gam
