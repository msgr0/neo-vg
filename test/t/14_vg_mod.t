#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order

plan tests 32

is $(nvg construct -r small/x.fa -v small/x.vcf.gz | nvg mod -k x - | nvg view - | grep "^P" | cut -f 3 | grep -o "[0-9]\+" |  wc -l) \
    $(nvg construct -r small/x.fa -v small/x.vcf.gz | nvg mod -k x - | nvg view - | grep "^S" | wc -l) \
    "nvg mod yields a graph with only a particular path"

is $(nvg mod graphs/orphans.vg | nvg view - | wc -l) 8 "orphan edge removal is automatic"

is $(nvg construct -m 1000 -r small/x.fa -v small/x.vcf.gz | nvg mod -pl 10 -e 3 - | nvg stats -E - ) 285 "graph complexity reduction works as expected"

is $(nvg construct -r small/x.fa -v small/x.vcf.gz | nvg mod -pl 10 -e 3 -t 16 - | nvg mod -S -l 200 - | nvg stats -l - | cut -f 2) 983 "short subgraph pruning works"

is $(nvg mod -U 10 graphs/q_redundant.vg | nvg view - | grep ^S | wc -l) 4 "normalization produces the correct number of nodes"

is $(nvg view -vF graphs/redundant-snp.gfa | nvg mod -n - | nvg view - | grep ^S | wc -l) 4 "normalization removes redundant SNP alleles"

nvg mod -n graphs/q_redundant.vg | nvg validate -
is $? 0 "normalization produces a valid graph"

nvg mod -U 10 graphs/q_redundant.vg | nvg validate -
is $? 0 "looped normalization produces a valid graph"

nvg mod -u graphs/q_redundant.vg | nvg validate -
is $? 0 "unchop produces a valid graph"

is $(nvg mod -U 10 graphs/q_redundant.vg | nvg stats -l - | cut -f 2) 154 "normalization removes redundant sequence in the graph"

is $(nvg view -Fv graphs/normalize_me.gfa | nvg mod -U 10 - | nvg view - | sort | md5sum | cut -f 1 -d\ ) $(md5sum graphs/normalize_me.norm.gfa | cut -f 1 -d\ ) "normalization doesn't introduce cycles and does remove redundancy in bubbles"

# shows that after mod we have == numbers of path annotations and nodes
# in this one-path graph
is $(nvg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | nvg mod -N - | nvg view - | grep "^P" | cut -f 3 | grep -o "[0-9]\+" |wc -l) \
   $(nvg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | nvg mod -N - | nvg view - | grep "^S" | wc -l) \
   "nvg mod removes non-path nodes and edge"

set -o pipefail
nvg view -Jv reversing/reversing_path.json | nvg mod -X 3 - | nvg validate -
is "$?" "0" "chopping a graph works correctly with reverse mappings"
set +o pipefail

is $(nvg mod -X 5 msgas/s.vg | nvg mod -u - | nvg validate - && nvg mod -X 5 msgas/s.vg | nvg mod -u - | nvg paths -v - -X | nvg view -a - | jq '.sequence' | sort | md5sum | cut -f 1 -d\ ) 2f785068c91dbe84177c1fd679b6f133 "unchop correctly handles paths"

is "$(nvg view -Jv graphs/inv-mess.json | nvg mod -u - | nvg validate - && nvg view -Jv graphs/inv-mess.json | nvg mod -u - | nvg view - | sort | diff - correct/14_vg_mod/inv-mess-unchopped.gfa)" "" "unchop correctly handles a graph with an inversion"

is "$(nvg view -Jv reversing/double_reversing.json | nvg mod -u - | nvg stats -z - | grep "nodes" | cut -f2)" "1" "unchop handles doubly-reversing edges"

is "$(nvg view -Jv graphs/inv-mess.json | nvg mod -U 10 - | nvg validate - && nvg view -Jv graphs/inv-mess.json | nvg mod -U 10 - | nvg view - | sort | diff - correct/14_vg_mod/inv-mess-normalized.gfa)" "" "normalization works on a graph with an inversion"

nvg mod -u msgas/s2.vg >/dev/null
is $? 0 "mod successfully unchops a difficult graph"

nvg mod -X 8 msgas/l.vg | nvg validate -
is $? 0 "chopping self-cycling nodes retains the cycle"

nvg mod -U 3 graphs/atgclinv2.vg | nvg validate -
is $? 0 "unrolling works and produces a valid graph"

nvg mod -U 3 graphs/atgclinv2.vg | nvg mod -f 50 - | nvg validate -
is $? 0 "unfolding works and produces a valid graph"

is $(comm -12 <(nvg mod -d 7 graphs/atgclinv2.vg | nvg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) <(cat graphs/atgclinv2.vg| nvg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) | wc -l) $(cat graphs/atgclinv2.vg| nvg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq | wc -l) "dagify-unroll produces a graph with the same kmers as the original graph"

is $(comm -12 <(nvg mod -f 7 graphs/atgclinv2.vg | nvg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) <(cat graphs/atgclinv2.vg| nvg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) | wc -l) $(cat graphs/atgclinv2.vg| nvg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq | wc -l) "unfold produces a graph with the same kmers as the original graph"

nvg mod -f 6 graphs/atgclinv2.vg | nvg mod -d 5 - | nvg validate -
is $? 0 "dagify handles a graph with two strongly connected components"

is $(nvg mod -f 6 graphs/atgclinv2.vg | nvg mod -d 5 - | nvg stats -c - | wc -l) $(nvg mod -f 6 graphs/atgclinv2.vg | nvg mod -d 5 - | nvg stats -N -) "unfold followed by dagify produces a graph with no cycles"

is $(comm -12 <(nvg mod -f 7 graphs/atgclinv2.vg | nvg mod -d 7 - | nvg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) <(cat graphs/atgclinv2.vg| nvg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) | wc -l) $(cat graphs/atgclinv2.vg| nvg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq | wc -l) "unfold followed by dagify produces a graph with the same kmers as the original graph"

nvg mod -w 100 graphs/ununrollable.vg | nvg validate -
is $? 0 "dagify unrolls the un-unrollable graph"

nvg mod -s graphs/not-simple.vg | nvg validate -
is $? 0 "sibling simplification does not disrupt paths"

nvg mod -N msgas/c2.vg | nvg paths -d -v - | nvg mod -U 10 - >c.vg
is $(cat c.vg | nvg mod -w 100 - | nvg stats -N -) 4 "dagify correctly calculates the minimum distance through the unrolled component"
is $(cat c.vg | nvg mod -w 100 - | nvg stats -l - | cut -f 2) 200 "dagify produces a graph of the correct size"
rm -f c.vg

# Note the math (and subsetting) only works out on a flat alleles graph
nvg construct -r small/x.fa -a -f -v small/x.vcf.gz >x.vg
nvg mod -v small/x.vcf.gz x.vg >x.sample.vg
is "$(nvg stats x.sample.vg -l | cut -f 2)" "1041" "subsetting a flat-alleles graph to a sample graph works"
rm -f x.vg x.sample.vg 

is $(nvg mod -M 5 jumble/j.vg|  nvg stats -s - | wc -l) 7 "removal of high-degree nodes results in the expected number of subgraphs"
