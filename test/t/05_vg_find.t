#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 31

nvg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
is $? 0 "construction"

nvg index -x x.xg x.vg 2>/dev/null
is $(nvg find -x x.xg -p x:200-300 -c 2 | nvg view - | grep CTACTGACAGCAGA | cut -f 2) 72 "a path can be queried from the xg index"
# todo: I complain because context expansion no longer puts ranks in paths, so edge check does't see that path is discontinuous
is $(nvg find -x x.xg -n 203 -c 1 | nvg view - | grep CTACCCAGGCCATTTTAAGTTTCCTGT | wc -l) 1 "a node near another can be obtained using context from the xg index"

nvg index -x x.xg -g x.gcsa -k 16 x.vg
is $(( for seq in $(nvg sim -l 50 -n 100 -x x.xg); do nvg find -M $seq -g x.gcsa; done ) | jq length | grep ^1$ | wc -l) 100 "each perfect read contains one maximal exact match"

nvg index -x x.xg -g x.gcsa -k 16 x.vg
is $(nvg find -n 1 -n 3 -D -x x.xg ) 8 "nvg find -D finds approximate distance between 2 adjacent node starts"
is $(nvg find -n 1 -n 2 -D -x x.xg ) 8 "nvg find -D finds approximate distance between node start and adjacent snp"
is $(nvg find -n 17 -n 20 -D -x x.xg ) 7 "nvg find -D jumps deletion"
is $(nvg find -n 16 -n 20 -D -x x.xg ) 7 "nvg find -D jumps deletion from other allele in snp"

is $(nvg find -n 2 -n 3 -c 1 -L -x x.xg | nvg view -g - | grep "^S" | wc -l) 5 "nvg find -L finds same number of nodes (with -c 1)"

is $(nvg find -r 6:5 -L -x x.xg | nvg view -g - | grep S | wc -l) 3 "nvg find -L works with -r.  it scans from start position of first node in range "

rm -f x.idx x.xg x.gcsa x.gcsa.lcp x.vg

nvg index -x m.xg inverting/m.vg
is $(nvg find -n 2308 -c 10 -L -x m.xg | nvg view -g - | grep S | wc -l) 5 "nvg find -L tracks length from start position of input node"
is $(nvg find -n 2315 -n 183 -n 176 -c 1 -L -x m.xg | nvg view -g - | grep S | wc -l) 7 "nvg find -L works with more than one input node"
rm m.xg

nvg construct -m 1000 -rmem/h.fa >h.vg
nvg index -g h.gcsa -k 16 h.vg
is $(nvg find -M ACCGTTAGAGTCAG -g h.gcsa) '[["ACC",["1:-32"]],["CCGTTAG",["1:5"]],["GTTAGAGT",["1:19"]],["TAGAGTCAG",["1:40"]]]' "we find the 4 canonical SMEMs from @lh3's bwa mem poster"
rm -f h.gcsa h.gcsa.lcp h.vg

nvg construct -m 1000 -r minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz -m 64 >giab.vg
nvg index -x giab.xg -g giab.gcsa -k 11 giab.vg
is $(nvg find -M ATTCATNNNNAGTTAA -g giab.gcsa | md5sum | cut -f -1 -d\ ) $(md5sum correct/05_vg_find/28.txt | cut -f -1 -d\ ) "we can find the right MEMs for a sequence with Ns"
is $(nvg find -M ATTCATNNNNAGTTAA -g giab.gcsa | md5sum | cut -f -1 -d\ ) $(nvg find -M ATTCATNNNNNNNNAGTTAA -g giab.gcsa | md5sum | cut -f -1 -d\ ) "we find the same MEMs sequences with different lengths of Ns"
rm -f giab.vg giab.xg giab.gcsa{,.lcp}

#nvg construct -r mem/r.fa > r.vg
#nvg index -x r.xg -g r.gcsa -k 16 r.vg
#is $(nvg find -M GAGGCAGTGAAGAGATCGTGGGAGGGAC -R 10 -Z 10 -g r.gcsa -x r.xg) '[["GAGGCAGTGAAGAGATCGTGGGAG",["1:20"]],["GCAGTGAAGAGATCGTGGGAGGGAC",["1:43"]],["CAGTGAAGAGATCGTGGGAGG",["1:0"]]]' "we can find sub-MEMs and only return hits that are not within containing MEMs"
#is $(nvg find -M GAGGCAGTGAAGAGATCGTGGGAGGGAC -R 10 -Z 10 -f -g r.gcsa -x r.xg) '[["GAGGCAGTGAAGAGATCGTGGGAG",["1:20"]],["GCAGTGAAGAGATCGTGGGAGGGAC",["1:43"]],["CAGTGAAGAGATCGTGGGAGG",["1:0"]]]' "we get the same (sufficiently long) MEMs with the fast sub-MEM option"
#rm -rf r.vg r.xg r.gcsa* mem/r.fa.fai

#nvg view -J mem/s.json -v > s.vg
#nvg index -x s.xg -g s.gcsa -k 16 s.vg
#is $(nvg find -M ACGTGCCGTTAGCCAGTGGGTTAG -R 10 -Z 10 -x s.xg -g s.gcsa) '[["ACGTGCCGTTAGCCAGTGGGTTAG",["3:11"]],["AGCCAGTGGGTTA",["1:0","2:0"]]]' "we can find sub-MEMs in a hard case of de-duplication"
#is $(nvg find -M ACGTGCCGTTAGCCAGTGGGTTAG -R 10 -Z 10 -f -x s.xg -g s.gcsa) '[["ACGTGCCGTTAGCCAGTGGGTTAG",["3:11"]],["AGCCAGTGGGTTA",["1:0","2:0"]]]' "we can find the same (sufficiently long) sub-MEMs the hard de-duplication case with the fast sub-MEM option"
#rm -rf s.vg s.xg s.gcsa*

nvg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
nvg index -x x.xg -g x.gcsa -k 11 x.vg
nvg map -x x.xg -g x.gcsa -T small/x-s1337-n100.reads >x.gam

nvg gamsort -i x.sorted.gam.gai x.gam > x.sorted.gam
is $(nvg find -o 127 --sorted-gam x.sorted.gam | nvg view -a - | wc -l) 6 "the GAM index can return the set of alignments mapping to a node"
is $(nvg find -A <(nvg find -N <(seq 37 52 ) -x x.xg ) --sorted-gam x.sorted.gam | nvg view -a - | wc -l) 15 "a subgraph query may be used to obtain a particular subset of alignments from a sorted GAM"

rm -rf x.gam x.sorted.gam x.sorted.gam.gai

is "$(nvg find -G small/x-s1337-n1.gam -x x.xg | nvg view - | grep ATTAGCCATGTGACTTTGAACAAGTTAGTTAATCTCTCTGAACTTCAGTT | wc -l)" "1" "the index can be queried using GAM alignments"

is "$(nvg find -G small/x-s1337-n1.gam -x x.xg | nvg paths --list -x -)" "x[121-272]" "querying the index with GAM alignments finds relevant subpaths"

rm -rf x.vg x.xg x.gcsa{,.lcp}

nvg construct -m 1000 -r tiny/tiny.fa -v tiny/tiny.vcf.gz >tiny.vg
nvg index -x tiny.xg tiny.vg 
is $(nvg find -x tiny.xg -n 12 -n 13 -n 14 -n 15 | nvg view - | grep ^L | wc -l) 4 "find gets connected edges between queried nodes by default"
echo 12 13 >get.nodes
echo 14 >>get.nodes
echo 15 >>get.nodes
is $(nvg find -x tiny.xg -N get.nodes | nvg view - | grep ^S | wc -l) 4 "find gets nodes provided in a node file list"
rm -rf tiny.xg tiny.vg get.nodes

echo '{"node": [{"id": 1, "sequence": "A"}, {"id": 2, "sequence": "A"}], "edge": [{"from": 1, "to": 2}], "path": [{"name": "ref", "mapping": [{"position": {"node_id": 1}}]}]}' | nvg view -Jv - >test.vg
nvg index -x test.xg test.vg
is "$(nvg find -x test.xg -p ref:0 -c 10 | nvg view -j - | jq '.edge | length')" "1" "extracting by path adds no extra edges"

rm -f test.vg test.xg

nvg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >t.vg
nvg index -x t.xg t.vg
is $(nvg find -x t.xg -E -p x:30-35 | nvg view - | grep ^S | wc -l) 4 "path DAG range query works"

nvg find -x t.xg -E -p x:30-35 -p x:10-20 -W t.
is $((nvg view t.x:30:35.vg; nvg view t.x:10:20.vg) | wc -l) 20 "we can extract a set of targets to separate files"

echo x 30 36 | tr ' ' '\t' >t.bed
echo x 10 21 | tr ' ' '\t' >>t.bed
nvg find -x t.xg -E -R t.bed -W q.
is $((nvg view q.x:10:20.vg; nvg view q.x:30:35.vg) | md5sum | cut -f 1 -d\ ) $((nvg view t.x:10:20.vg ; nvg view t.x:30:35.vg)| md5sum | cut -f 1 -d\ ) "the same extraction can be made using BED input"

is $(nvg find -x t.xg -E -p x:30-35 -p x:10-20 -K 5 | wc -l) 36 "we see the expected number of kmers in the given targets"

rm -f t.xg t.vg t.x:30:35.vg t.x:10:20.vg q.x:30:35.vg q.x:10:20.vg t.bed

nvg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
nvg index -x x.xg x.vg
nvg gbwt -v small/xy2.vcf.gz -o x.gbwt -x x.vg
is $(nvg find -p x -x x.xg -K 16 -H x.gbwt | cut -f 5 | sort | uniq -c  | tail -n 1 | awk '{ print $1 }') 1510 "we find the expected number of kmers with haplotype frequency equal to 2"
rm -f x.vg x.xg x.gbwt

nvg convert -gp tiny/tiny.gfa | nvg find -x - -n 1 -c 2 | nvg convert -f - | nvg ids -s - | sort > found1.gfa
nvg find -x tiny/tiny.gfa -n 1 -c 2| nvg ids -s - | sort > found2.gfa
diff found1.gfa found2.gfa
is $? 0 "GFA i/o for find -n consistent with converting both ways"
nvg find -x tiny/tiny.gfaz -n 1 -c 2 >/dev/null
is $? 0 "find accepts GFAZ graph input"


# Find nodes that map to the provided ids
nvg construct -m 32 -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
nvg gbwt -v small/xy2.vcf.gz -o x.gbwt -x x.vg
nvg prune -u -m x.mapping -g x.gbwt -e 1 x.vg > x.unfolded.vg

rm -f expected.gfa
printf "H\tVN:Z:1.1\n" >> expected.gfa
printf "S\t72\tTGGAGTTCTATTATATTCC\n" >> expected.gfa
printf "S\t76\tTGGAGTTCTATTATATTCC\n" >> expected.gfa
printf "S\t97\tTCT\n" >> expected.gfa
printf "S\t98\tTCT\n" >> expected.gfa
printf "S\t82\tTGGAGTTCTATTATATTCC\n" >> expected.gfa
nvg find -n 5 -n 23 --mapping x.mapping -x x.unfolded.vg 2> /dev/null | nvg view -g - > found.gfa
cmp found.gfa expected.gfa
is $? 0 "find nodes that map to the provided node ids"

rm -f x.vg x.gbwt x.mapping x.unfolded.vg
rm -f expected.gfa found.gfa

# We wish we could test specifically for which paths are indexed, but we can't.
# So we test to make sure at least paths that shouldn't normally be indexed are
# indexed when asked about.
is "$(nvg find -n 5 -P sample1#1#chr1#0 -x graphs/gfa_with_reference.gfa | cut -f2)" "4" "nvg find can find positions along non-reference paths asked about specifically"

