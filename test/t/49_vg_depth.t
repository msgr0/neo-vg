#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 9

nvg construct -m 10 -r tiny/tiny.fa >flat.vg
nvg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | nvg view -Fv - >2snp.vg
nvg sim -l 30 -x 2snp.vg -n 30 -a -s 1 >2snp.sim
nvg index -x flat.xg -g flat.gcsa -k 16 flat.vg
nvg map -g flat.gcsa -x flat.xg -G 2snp.sim -k 8 >2snp.gam
nvg pack -x flat.xg -o 2snp.gam.cx -g 2snp.gam
# total read bases (30 * 30) / total graph bases 50 = 18
is $(nvg depth flat.vg -g 2snp.gam | awk '{print $1}') 18 "nvg depth gets correct depth from gam"
is $(nvg depth flat.xg -k 2snp.gam.cx -b 100000 | awk '{print int($4)}') 18 "nvg depth gets correct depth from pack"
is $(nvg depth flat.xg -k 2snp.gam.cx -b 10 | wc -l) 5 "nvg depth gets correct number of bins"
nvg convert flat.vg -G 2snp.gam | gzip > 2snp.gaf.gz
is $(nvg depth flat.vg -a 2snp.gaf.gz | awk '{print $1}') 18 "nvg depth gets correct depth from gaf"
nvg augment flat.vg 2snp.gam -i > flat-aug.vg
is $(nvg depth flat-aug.vg | awk '{print $1}' | uniq | wc -l) $(nvg paths -Lv flat-aug.vg | wc -l) "nvg depth of paths reports all paths"
is $(nvg depth flat-aug.vg -P x | awk '{print $1}' | uniq | wc -l) 1 "nvg depth of paths reports just path with selected prefix"
rm -f flat.vg flat.gcsa flat.xg 2snp.vg 2snp.sim 2snp.gam 2snp.gam.cx 2snp.gaf.gz flat-aug.vg

# Build a GBZ with a generic reference 'x' and two haplotype paths from x.vcf samples
nvg autoindex -r small/x.fa -v small/x.vcf.gz -w giraffe -p depth_x 2>/dev/null
mv depth_x.giraffe.gbz depth_x.gbz

# Path coverage on a GBZ counts haplotype paths.  'x' is the only generic path,
# so without haplotype counting max coverage on x would be 0; with the two
# haplotypes from the VCF, positions shared by both report coverage 2.
is $(nvg depth -p x depth_x.gbz | awk '{print $3}' | sort -nr | head -1) 2 "nvg depth path coverage on a GBZ counts haplotype paths"

# Pack coverage on a GBZ uses pack values, not path counts.  Sim 100 30bp reads
# from x giving ~3x mean coverage; a path-count-based result would be 2 instead.
nvg convert depth_x.gbz -p > depth_x.pg
nvg sim -x depth_x.pg -l 30 -n 100 -a -s 1 -P x > depth_x.sim.gam
nvg pack -x depth_x.pg -o depth_x.pack -g depth_x.sim.gam
is $(nvg depth depth_x.gbz -k depth_x.pack -b 100000 -p x | awk '{print int($4)}') 3 "nvg depth pack coverage on a GBZ is unaffected by haplotype paths"

# Parallel path depth must be deterministic: -t 1 and -t 4 should produce
# byte-identical output (the outer ordered parallel for preserves path order).
nvg depth -t 1 depth_x.gbz > depth_x.t1.out
nvg depth -t 4 depth_x.gbz > depth_x.t4.out
cmp depth_x.t1.out depth_x.t4.out
is $? 0 "nvg depth path coverage is deterministic across thread counts"

rm -f depth_x.gbz depth_x.pg depth_x.sim.gam depth_x.pack depth_x.dist depth_x.shortread.withzip.min depth_x.shortread.zipcodes depth_x.t1.out depth_x.t4.out
