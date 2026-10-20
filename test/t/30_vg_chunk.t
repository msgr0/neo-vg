#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 46

# Construct a graph with alt paths so we can make a GBWT and a GBZ
nvg construct -m 1000 -r small/x.fa -v small/x.vcf.gz -a >x.vg
nvg index -x x.xg x.vg
nvg gbwt -x x.vg -v small/x.vcf.gz -o x.haps.gbwt
nvg gbwt -x x.vg -E -o x.paths.gbwt
nvg gbwt -m x.haps.gbwt x.paths.gbwt -o x.gbwt
nvg gbwt -x x.vg x.gbwt --gbz-format -g x.gbz

nvg view -a small/x-l100-n1000-s10-e0.01-i0.01.gam > x.gam.json

# sanity check: does passing no options preserve input
is $(nvg chunk -x x.xg -p x -c 10| nvg stats - -N) 210 "nvg chunk with no options preserves nodes"
is $(nvg chunk -x x.xg -p x -c 10| nvg stats - -E) 291 "nvg chunk with no options preserves edges"

# check a small chunk
is $(nvg chunk -x x.xg -p x:20-30 -c 0 | nvg view - -j | jq -c '.path[0].mapping[].position' | jq 'select ((.node_id == "9"))' | grep node | sed s/,// | sort | uniq | wc -l) 1 "chunk has path going through node 9"
is $(nvg chunk -x x.gbz -p x:20-30 -c 0 | nvg view - -j | jq -c '.path[0].mapping[].position' | jq 'select ((.node_id == "9"))' | grep node | sed s/,// | sort | uniq | wc -l) 1 "chunk can be found from GBZ"

# check snarl chunking
nvg snarls x.xg > x.snarls
is $(nvg chunk -x x.xg -S x.snarls -p x:10-20 | nvg view - | grep ^S | awk '{print $2}' | sort -n | awk 'BEGIN { ORS = "" } {print}') 6789 "snarl chunk works on simple example"
rm -f x.snarls

# check a small chunk, but using vg input and packed graph output
is $(nvg chunk -x x.vg -p x:20-30 -c 0 -O pg | nvg convert -v - | nvg view - -j | jq -c '.path[0].mapping[].position' | jq 'select ((.node_id == "9"))' | grep node | sed s/,// | sort | uniq | wc -l) 1 "chunk has path going through node 9"

# check no crash when using chunk_size, and filenames deterministic
rm -f _chunk_test*
nvg chunk -x x.xg -p x -s 233 -o 50 -b _chunk_test -c 0 -t 2
nvg chunk -x x.xg -p x -s 233 -o 50 -b _chunk_test -c 0 -t 1
is $(ls -l _chunk_test*.vg | wc -l) 6 "-s produces correct number of chunks"
rm -f _chunk_test*

#check that gam chunker runs through without crashing
nvg gamsort small/x-l100-n1000-s10-e0.01-i0.01.gam -i x.sorted.gam.gai > x.sorted.gam
printf "x\t2\t200\nx\t500\t600\n" > _chunk_test_bed.bed
nvg chunk -x x.xg -a x.sorted.gam -g -b _chunk_test -e _chunk_test_bed.bed -E _chunk_test_out.bed -c 0
is $(ls -l _chunk_test*.vg | wc -l) 2 "gam chunker produces correct number of graphs"
is $(ls -l _chunk_test*.gam | wc -l) 2 "gam chunker produces correct number of gams"
is $(grep x _chunk_test_out.bed | wc -l) 2 "gam chunker produces bed with correct number of chunks"
is "$(nvg view -aj _chunk_test_0_x_0_199.gam | wc -l)" "$(nvg view -aj _chunk_test_0_x_0_199.gam | sort | uniq | wc -l)" "gam chunker emits each matching read at most once"
is "$(nvg view -aj _chunk_test_1_x_500_627.gam | wc -l)" "225" "chunk contains the expected number of alignments"
rm -f _chunk_test*

#check that we can chunk by read count
nvg chunk -a small/x-l100-n1000-s10-e0.01-i0.01.gam -m 100 -b _chunk_test
is $(ls -l _chunk_test*.gam | wc -l) 10 "simple gam chunker produces correct number of gams"
is "$(nvg view -aj _chunk_test000005.gam | wc -l)" "100" "simple chunk contains the expected number of alignments"

#check that id ranges work
is $(nvg chunk -x x.xg -r 1:3 -c 0 | nvg view - -j | jq .node | grep id |  wc -l) 3 "id chunker produces correct chunk size"
is $(nvg chunk -x x.xg -r 1 -c 0 | nvg view - -j | jq .node | grep id | wc -l) 1 "id chunker produces correct single chunk"

# Check that traces work without any haplotypes
# TODO: Make the tube map stop doing that?
is $(nvg chunk -x x.xg -r 1:1 -c 2 -T | nvg view - -j | jq .node | grep id | wc -l) 5 "id chunker traces correct chunk size without haplotypes"

# Check that traces work on a GBWT
is $(nvg chunk -x x.xg -G x.gbwt -r 1:1 -c 2 -T | nvg view - -j | jq .node | grep id | wc -l) 5 "id chunker traces correct chunk size"
is "$(nvg chunk -x x.xg -r 1:1 -c 2 -T | nvg view - -j | jq -c '.path[] | select(.name | startswith("x[0") | not)' | wc -l)" 0 "chunker extracts no threads from an empty gPBWT"
is "$(nvg chunk -x x.xg -G x.haps.gbwt -r 1:1 -c 2 -T | nvg view - -j | jq -c '.path[] | select(.name | startswith("x[0") | not)' | wc -l)" 2 "chunker extracts 2 local threads from a gBWT with 2 locally distinct threads in it"
is "$(nvg chunk -x x.xg -G x.gbwt -r 1:1 -c 2 -T | nvg view - -j | jq -r '.path[] | select(.name == "thread_0") | .mapping | length')" 3 "chunker can extract a partial haplotype from a GBWT"
is "$(nvg chunk -x x.gbz -r 1:1 -c 2 -T | nvg view - -j | jq -r '.path[] | select(.name == "thread_0") | .mapping | length')" 3 "chunker can extract a partial haplotype from a GBZ"
is "$(nvg chunk -x x.gbz -r 1:1 -c 2 -T --no-embedded-haplotypes | nvg view - -j | jq -r '.path[] | select(.name == "thread_0") | .mapping | length')" "" "chunker doesn't see haplotypes in the GBZ if asked not to"

nvg chunk -x x.xg -G x.gbwt -p x:50-100 -c 0 -T | nvg view - | grep ^S | sort > cnodes
nvg find -x x.xg -p x:50-100 -c 0 | nvg view - | grep ^S | sort > fnodes
diff cnodes fnodes
is "$?" 0 "trace consistent on path coordinates"
rm -f cnodes fnodes

#check that n-chunking works
# We know that it will drop _alt paths so we remake the graph without them for comparison.
nvg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
mkdir -p x.chunk
nvg chunk -x x.xg -n 5 -b x.chunk/ -O vg
# nvg chunk no longer keeps ranks, so there's no way to stitch the path back together anymore!
# we grep it out of the comparison
is $(cat x.chunk/*vg | nvg view -V - | grep -v P 2>/dev/null | sort |  md5sum | cut -f 1 -d\ ) $(nvg view x.vg | grep -v P | sort  | md5sum | cut -f 1 -d\ ) "n-chunking works and chunks over the full graph"

rm -rf x.sorted.gam x.sorted.gam.gai _chunk_test_bed.bed _chunk_test* x.chunk
rm -f x.vg x.xg x.gbwt x.gbz x.gam.json filter_chunk*.gam chunks.bed
rm -f chunk_*.annotate.txt

nvg construct -r small/xy.fa -v small/xy.vcf.gz > xy.vg
nvg construct -r small/xy.fa -v small/xy.vcf.gz -R x > x.vg
nvg construct -r small/xy.fa -v small/xy.vcf.gz -R y > y.vg
nvg ids -j x.vg y.vg
nvg sim -x x.vg -n 50 -a > x.gam
nvg sim -x y.vg -n 100 -a > y.gam
cat x.gam y.gam > xy.gam
# test that exploding into components works
nvg chunk -x xy.vg -M -b path_chunk -O hg -a xy.gam -g
nvg view x.vg | grep "^S" | awk '{print $3}' | sort > x_nodes.txt
nvg view y.vg | grep "^S" | awk '{print $3}' | sort > y_nodes.txt
nvg convert path_chunk_x.vg -v | nvg view - | grep "^S" | awk '{print $3}' | sort > pc_x_nodes.txt
nvg convert path_chunk_y.vg -v | nvg view - | grep "^S" | awk '{print $3}' | sort > pc_y_nodes.txt
diff x_nodes.txt pc_x_nodes.txt && diff y_nodes.txt pc_y_nodes.txt
is "$?" 0 "path-based components finds subgraphs"
is $(nvg view -a path_chunk_x.gam | wc -l) $(nvg view -a x.gam | wc -l) "x gam chunk has correct number of reads"
is $(nvg view -a path_chunk_y.gam | wc -l) $(nvg view -a y.gam | wc -l) "y gam chunk has correct number of reads"
nvg chunk -x xy.vg -C -p x -b path_chunk_ind -O hg -a xy.gam -g > /dev/null
nvg chunk -x xy.vg -C -p y -b path_chunk_ind -O hg -a xy.gam -g > /dev/null
is $(nvg view -a path_chunk_ind_x.gam | wc -l) $(nvg view -a x.gam | wc -l) "x gam chunk has correct number of reads with -p path"
is $(nvg view -a path_chunk_ind_y.gam | wc -l) $(nvg view -a y.gam | wc -l) "y gam chunk has correct number of reads with -p path"
nvg paths -v x.vg -E > x_paths.txt
nvg paths -v path_chunk_x.vg -E > pc_x_paths.txt
diff pc_x_paths.txt x_paths.txt
is "$?" 0 "path-based component contains correct path length"
nvg chunk -x xy.vg -C -b components_chunk
nvg view components_chunk_0.vg | grep "^S" | awk '{print $3}' > comp_0_nodes.txt
nvg view components_chunk_1.vg | grep "^S" | awk '{print $3}' > comp_1_nodes.txt
cat comp_0_nodes.txt comp_1_nodes.txt | sort > comp_nodes.txt
cat x_nodes.txt y_nodes.txt | sort > nodes.txt
diff comp_nodes.txt nodes.txt
is "$?" 0 "components finds subgraphs"

rm -f xy.vg x.vg y.vg x_nodes.txt y_nodes.txt convert path_chunk_x.vg  convert path_chunk_y.vg pc_x_nodes.txt pc_y_nodes.txt x_paths.txt pc_x_paths.txt components_chunk_0.vg components_chunk_1.vg comp_0_nodes.txt comp_1_nodes.txt comp_nodes.txt nodes.txt x.gam y.gam xy.gam path_chunk_x.gam path_chunk_y.gam

nvg gbwt --gbz-format --graph-name graph.gbz --gfa-input graphs/gfa_with_reference.gfa
nvg chunk -x graph.gbz -p sample1#1#chr1#0:1-2 -c 1 >part.vg 2>log.txt
grep "out_of_range" log.txt
is "$?" "1" "chunking on a haplotype path does not produce an out of range error"
grep "not found in" log.txt
is "$?" "1" "chunking on a haplotype path does not produce a path not found error"
is "$(nvg stats -z part.vg | grep nodes | cut -f2)" "4" "chunking on a haplotype produces the correct size graph"

nvg chunk -x graph.gbz -p GRCh38#0#chr1:5-7 -c 1 >part.vg
nvg chunk -x part.vg -p GRCh38#0#chr1:5-7 -c 1 >subpart.vg
is "$?" "0" "chunking the same base path range out of a chunk works"
is "$(nvg stats -l part.vg)" "$(nvg stats -l subpart.vg)" "chunking the same base path range out of a chunk gets the same nodes"
is "$(nvg paths -L -H -x part.vg | wc -l)" "0" "chunking a GBZ leaves out haplotype paths by default"
nvg chunk -x graph.gbz -p GRCh38#0#chr1:5-7 -c 1 -H >part.vg
is "$(nvg paths -L -H -x part.vg | wc -l)" "3" "chunking a GBZ with -H keeps the haplotype paths"

nvg convert -g graphs/gfa_with_reference.gfa >graph.vg
nvg chunk -x graph.vg -p GRCh38#0#chr1:5-7 -c 1 >part.vg
is "$(nvg paths -L -H -x part.vg | wc -l)" "0" "chunking a non-GBZ graph leaves out haplotype paths by default"
nvg chunk -x graph.vg -p GRCh38#0#chr1:5-7 -c 1 -H >part.vg
is "$(nvg paths -L -H -x part.vg | wc -l)" "3" "chunking a non-GBZ graph with -H keeps the haplotype paths"

rm -f graph.gbz graph.vg part.vg subpart.vg log.txt

# GBZ chunking
nvg gbwt --graph-name graph.gbz --set-reference sample --gfa-input graphs/components_walks.gfa
nvg chunk -x graph.gbz --gbz
is "$?" "0" "chunking a GBZ with option --gbz works"
is "$(ls -l chunk_*.gbz | wc -l)" "2" "GBZ chunking produces correct number of chunks"
nvg chunk -x graph.gbz --gbz --contig A --prefix single
is "$?" "0" "chunking a GBZ with options --gbz and --contig works"
is "$(ls -l single_*.gbz | wc -l)" "1" "GBZ chunking with --contig produces correct number of chunks"
cmp chunk_0_A.gbz single_0_A.gbz
is "$?" "0" "GBZ chunking with --contig produces the expected chunk"

rm -f graph.gbz chunk_*.gbz single_*.gbz
