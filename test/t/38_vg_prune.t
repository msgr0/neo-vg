#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 21


# Build a graph with one path and two threads
nvg construct -m 32 -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
nvg gbwt -o x.gbwt -v small/xy2.vcf.gz -x x.vg

# Basic pruning: 5 components, 51 nodes, 51 edges
nvg prune -e 1 x.vg > y.vg
is $(nvg stats -s y.vg | wc -l) 5 "pruning produces the correct number of components"
is $(nvg stats -N y.vg) 51 "pruning leaves the correct number of nodes"
is $(nvg stats -E y.vg) 51 "pruning leaves the correct number of edges"
rm -f y.vg

# Remove high-degree nodes: 6 components, 50 nodes, 47 edges
nvg prune -e 1 -M 3 x.vg > y.vg
is $(nvg stats -s y.vg | wc -l) 6 "pruning without high-degree nodes produces the correct number of components"
is $(nvg stats -N y.vg) 50 "pruning without high-degree nodes leaves the correct number of nodes"
is $(nvg stats -E y.vg) 47 "pruning without high-degree nodes leaves the correct number of edges"
rm -f y.vg

# Restore paths: 1 component, 64 nodes, 68 edges
nvg prune -r -e 1 x.vg > y.vg
is $(nvg stats -s y.vg | wc -l) 1 "pruning with path restoring produces the correct number of components"
is $(nvg stats -N y.vg) 64 "pruning with path restoring leaves the correct number of nodes"
is $(nvg stats -E y.vg) 68 "pruning with path restoring leaves the correct number of edges"
rm -f y.vg

# Unfold paths and threads: 1 component, 80 nodes, 92 edges
nvg prune -u -m x.mapping -g x.gbwt -e 1 x.vg > y.vg
is $(nvg stats -s y.vg | wc -l) 1 "pruning with path/thread unfolding produces the correct number of components"
is $(nvg stats -N y.vg) 80 "pruning with path/thread unfolding produces the correct number of nodes"
is $(nvg stats -E y.vg) 92 "pruning with path/thread unfolding produces the correct number of edges"
rm -f x.mapping y.vg

# Unfold only paths: 1 component, 64 nodes, 68 edges
nvg prune -u -m x.mapping -e 1 x.vg > y.vg
is $(nvg stats -s y.vg | wc -l) 1 "pruning with path unfolding produces the correct number of components"
is $(nvg stats -N y.vg) 64 "pruning with path unfolding produces the correct number of nodes"
is $(nvg stats -E y.vg) 68 "pruning with path unfolding produces the correct number of edges"
rm -f x.mapping y.vg

rm -f x.vg x.gbwt


# Build vg graphs for two chromosomes
nvg construct -m 32 -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
nvg construct -m 32 -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
nvg ids -j -m xy.mapping x.vg y.vg
nvg gbwt -o x.gbwt -v small/xy2.vcf.gz -x x.vg
nvg gbwt -o y.gbwt -v small/xy2.vcf.gz -x y.vg

# Prune a single-chromosome graph using multi-chromosome GBWT
nvg gbwt -m -o xy.gbwt x.gbwt y.gbwt
nvg prune -u -m x.mapping -g xy.gbwt -e 1 x.vg > pruned.vg
is $(nvg stats -s pruned.vg | wc -l) 1 "unfolding with multi-chromosome GBWT produces the correct number of components"
is $(nvg stats -N pruned.vg) 80 "unfolding with multi-chromosome GBWT produces the correct number of nodes"
is $(nvg stats -E pruned.vg) 92 "unfolding with multi-chromosome GBWT produces the correct number of edges"

# Check the range in the node mapping at each point
is $(head -c 16 xy.mapping | hexdump -v -e '4/4 "%u_"') 139_0_139_0_ "empty mapping starts from the right node id"
nvg prune -u -g x.gbwt -e 1 -a -m xy.mapping x.vg > /dev/null
is $(head -c 16 xy.mapping | hexdump -v -e '4/4 "%u_"') 139_0_168_0_ "the first unfolded graph adds the correct number of nodes to the mapping"
nvg prune -u -g y.gbwt -e 1 -a -m xy.mapping y.vg > /dev/null
is $(head -c 16 xy.mapping | hexdump -v -e '4/4 "%u_"') 139_0_196_0_ "the second unfolded graph adds the correct number of nodes to the mapping"

rm -f x.vg y.vg x.gbwt y.gbwt
rm -f xy.gbwt x.mapping pruned.vg
rm -f xy.mapping
