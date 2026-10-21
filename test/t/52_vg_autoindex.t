#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 66

rm auto.*

nvg autoindex -p auto -w map -r tiny/tiny.fa -v tiny/tiny.vcf.gz --force-unphased
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg map with basic input"
is $(ls auto.* | wc -l) 3 "autoindexing makes 3 outputs for nvg map" 
is $(ls auto.xg | wc -l) 1 "autoindexing makes an XG for nvg map"
is $(ls auto.gcsa* | wc -l) 2 "autoindexing makes a GCSA2/LCP pair for nvg map"
nvg sim -x auto.xg -n 20 -a -l 10 | nvg map -d auto -t 1 -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by nvg map"

# Backdate the index so a rewrite is visible even with 1-second mtime resolution
touch -t 200001010000 auto.xg auto.xg.stamp
nvg autoindex -p auto -w map -r tiny/tiny.fa -v tiny/tiny.vcf.gz --force-unphased --no-guessing
same_time="yes"
if [ auto.xg -nt auto.xg.stamp ]; then
    same_time="no"
fi
is "$same_time" "no" "autoindexing with --no-guessing overwrites existing files"

rm auto.*

nvg autoindex -p auto -w map -r small/x.fa -v small/x.vcf.gz -r small/y.fa -v small/y.vcf.gz 
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg map with chunked input"
nvg sim -x auto.xg -n 20 -a -l 10 | nvg map -d auto -t 1 -G - > /dev/null
is $(echo $?) 0 "chunked autoindexing results can be used by nvg map"

rm auto.*

nvg autoindex -p auto -w map -r small/xy.fa -v small/xy2.vcf.gz
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg map with phased input"
nvg sim -x auto.xg -n 20 -a -l 10 | nvg map -d auto -t 1 -G - > /dev/null
is $(echo $?) 0 "phased autoindexing results can be used by nvg map"

rm auto.*

nvg autoindex -p auto -w mpmap -w rpvg -r tiny/tiny.fa -v tiny/tiny.vcf.gz -x tiny/tiny.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg mpmap with unchunked input"
is $(ls auto.* | wc -l) 6 "autoindexing creates 6 files for mpmap/rpvg"
nvg sim -x auto.spliced.xg -n 20 -a -l 10 | nvg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by nvg mpmap"
is $(nvg paths -g auto.haplotx.gbwt -L | wc -l) 6 "haplotype transcript GBWT made by autoindex is valid"
is $(cat auto.txorigin.tsv | wc -l) 7 "transcript origin table has expected number of rows" 

rm auto.*

nvg autoindex -p auto -w mpmap  -r tiny/tiny.fa -v tiny/tiny.vcf.gz -x tiny/tiny.gtf --force-unphased
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg mpmap with unchunked, unphased input"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 files for mpmap/rpvg"
nvg sim -x auto.spliced.xg -n 20 -a -l 10 | nvg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "basic unphased autoindexing results can be used by nvg mpmap"

rm auto.*

nvg autoindex -p auto -w mpmap -r tiny/tiny.fa -x tiny/tiny.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg mpmap without variants"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 files for mpmap without variants"
nvg sim -x auto.spliced.xg -n 20 -a -l 10 | nvg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "autoindexing results with no variants can be used by nvg mpmap"

rm auto.*

nvg autoindex -p auto -w mpmap -w rpvg -r small/x.fa -r small/y.fa -v small/x.vcf.gz -v small/y.vcf.gz -x small/x.gtf -x small/y.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg mpmap with chunked input"
is $(ls auto.* | wc -l) 6 "autoindexing creates 6 files for mpmap/rpvg with chunked input"

rm auto.*

nvg autoindex -p auto -w mpmap -r small/x.fa -r small/y.fa -v small/x.vcf.gz -v small/y.vcf.gz -x small/x.gtf -x small/y.gtf --force-unphased
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg mpmap with unphased chunked input"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 files for mpmap/rpvg with chunked input"
nvg sim -x auto.spliced.xg -n 20 -a -l 10 | nvg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "autoindexing results with chunked unphased input can be used by nvg mpmap"

rm auto.*

nvg autoindex -p auto -w mpmap -w rpvg -r small/x.fa -r small/y.fa -v small/x.vcf.gz -v small/y.vcf.gz -x small/xy.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg mpmap with partially chunked input"

rm auto.*

nvg autoindex -p auto -w mpmap -w rpvg -r small/xy.fa -v small/x.vcf.gz -v small/y.vcf.gz -x small/xy.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg mpmap with another partially chunked input"

rm auto.*

nvg autoindex -p auto -w mpmap -g graphs/gfa_with_w_lines.gfa -x graphs/gfa_with_w_lines.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg mpmap with GFA input"
nvg sim -x auto.spliced.xg -n 20 -a -l 5 | nvg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "autoindexing results with GFA input can be used by nvg mpmap"

rm auto.*

nvg autoindex -p auto -w mpmap -g graphs/gfa_with_w_lines.gfa -x graphs/gfa_with_w_lines.gtf -H graphs/gfa_with_w_lines_haplo.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg mpmap with GFA input and a haplotype GTF"
nvg sim -x auto.spliced.xg -n 20 -a -l 5 | nvg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "autoindexing results with GFA input and a haplotype GTF can be used by nvg mpmap"

rm auto.*

nvg autoindex -p auto -w rpvg -w mpmap -g graphs/gfa_with_w_lines.gfa -x graphs/gfa_with_w_lines.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg mpmap and rpvg with GFA input"

rm auto.*

nvg autoindex -p auto -w rpvg -w mpmap -g graphs/gfa_with_w_lines.gfa -x graphs/gfa_with_w_lines.gtf -H graphs/gfa_with_w_lines_haplo.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg mpmap and rpvg with GFA input and a haplotype GTF"

rm auto.*

nvg autoindex -p auto -w giraffe -r tiny/tiny.fa -v tiny/tiny.vcf.gz 
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg giraffe with unchunked input"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 inputs for nvg giraffe"
is "$(nvg describe auto.giraffe.gbz | grep -c 'pggname =')" 1 "GBZ has a graph name"
is "$(nvg describe auto.shortread.withzip.min | grep 'pggname =')" "$(nvg describe auto.giraffe.gbz | grep 'pggname =')" "graph name was copied to minimizer index"
nvg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > t.vg
nvg index -x t.xg t.vg
nvg sim -x t.xg -n 20 -a -l 10 | nvg giraffe -Z auto.giraffe.gbz -m auto.shortread.withzip.min -z auto.shortread.zipcodes -d auto.dist -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by nvg giraffe"

rm auto.*
rm t.*

nvg autoindex -p auto -w sr-giraffe -r tiny/tiny.fa -v tiny/tiny.vcf.gz 
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg giraffe with unchunked input"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 inputs for short read nvg giraffe"
nvg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > t.vg
nvg index -x t.xg t.vg
nvg sim -x t.xg -n 20 -a -l 10 | nvg giraffe -Z auto.giraffe.gbz -m auto.shortread.withzip.min -z auto.shortread.zipcodes -d auto.dist -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by nvg giraffe"

rm auto.*
rm t.*

nvg autoindex -p auto -w lr-giraffe -r tiny/tiny.fa -v tiny/tiny.vcf.gz 
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg giraffe with unchunked input"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 inputs for long read nvg giraffe"
nvg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > t.vg
nvg index -x t.xg t.vg
nvg sim -x t.xg -n 20 -a -l 10 | nvg giraffe -Z auto.giraffe.gbz -m auto.longread.withzip.min -z auto.longread.zipcodes -d auto.dist -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by nvg giraffe"

rm auto.*
rm t.*

nvg autoindex -p auto -w giraffe -g graphs/gfa_with_w_lines.gfa 
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg giraffe with GFA input with W-lines"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 inputs for nvg giraffe from GFA input"
is "$(nvg describe auto.giraffe.gbz | grep -c 'pggname =')" 1 "GBZ has a graph name"
is "$(nvg describe auto.shortread.withzip.min | grep 'pggname =')" "$(nvg describe auto.giraffe.gbz | grep 'pggname =')" "graph name was copied to minimizer index"
nvg convert -g graphs/gfa_with_w_lines.gfa -x > g.xg
nvg sim -x g.xg -n 20 -a -l 4 | nvg giraffe -Z auto.giraffe.gbz -m auto.shortread.withzip.min -z auto.shortread.zipcodes -d auto.dist -G - > /dev/null
is $(echo $?) 0 "GFA autoindexing results can be used by nvg giraffe"

rm auto.*
rm g.xg

nvg autoindex -p auto -w giraffe -g graphs/components_walks_compressed.gfa
is $(echo $?) 0 "autoindexing successfully completes indexing for nvg giraffe with grammar-compressed GFA input"
nvg gbwt -G graphs/components_walks_compressed.gfa -g expected.gbz
cmp expected.gbz auto.giraffe.gbz
is $(echo $?) 0 "nvg autoindex handles grammar-compressed GFA input correctly"

rm auto.*
rm expected.gbz

nvg autoindex -p auto -w giraffe -g graphs/named_with_walk.gfa 
is $(echo $?) 0 "autoindexing successfully completes on a GFA with named segments and W-lines"
printf '@read\nGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGATTACACATTAGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n' > read.fq
nvg giraffe -Z auto.giraffe.gbz -m auto.shortread.withzip.min -d auto.dist -f read.fq --named-coordinates > read.gam
is "$(nvg view -aj read.gam | jq -r '.path.mapping[].position.name')" "Ishmael" "GFA segment names are available in output GAM when a walk exists"

rm auto.*
rm read.fq read.gam

nvg autoindex -p auto -w giraffe -g graphs/named.gfa 
is $(echo $?) 0 "autoindexing successfully completes on a GFA with named segments and no W-lines"
printf '@read\nGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGATTACACATTAGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n' > read.fq
nvg giraffe -Z auto.giraffe.gbz -m auto.shortread.withzip.min -d auto.dist -f read.fq --named-coordinates > read.gam
is "$(nvg view -aj read.gam | jq -r '.path.mapping[].position.name')" "Ishmael" "GFA segment names are available in output GAM when no walks exist"

# reduce disk limit to 1MB to trigger it during k-mer generation
is "$(nvg autoindex -p auto -w map --gcsa-size-limit 1000000 -g graphs/linked_cycles.gfa 2>&1 | grep Rewind | wc -l | tr -d " ")" 1 "Running out of room during k-mer enumeration triggers a rewind"
is "$(echo $?)" 0 "Indexing is successful after rewinding from k-mer generation"

rm -f auto.gcsa auto.gcsa.lcp

# reduce disk limit to 2MB to trigger it during doubling steps
is "$(nvg autoindex -p auto -w map --gcsa-size-limit 2000000 -g graphs/linked_cycles.gfa 2>&1 | grep Rewind | wc -l | tr -d " ")" 1 "Running out of room during GCSA2 indexing triggers a rewind"
is "$(echo $?)" 0 "Indexing is successful after rewinding from GCSA2 indexing"

rm -f auto.gcsa auto.gcsa.lcp

# use the memory limit to trigger a rewide
is "$(nvg autoindex -p auto -w map -M 512M -g graphs/linked_cycles.gfa 2>&1 | grep Rewind | wc -l | tr -d " ")" 1 "Running out of memory during GCSA2 indexing triggers a rewind"
is "$(echo $?)" 0 "Indexing is successful after rewinding from GCSA2 indexing"

rm auto.*

# index a graph with oversized snarls
nvg index -j auto.dist --snarl-limit 5 graphs/chain-clip.gfa 2> log.txt
is "$(grep 'oversized snarls' log.txt | wc -l | tr -d " ")" 1 "graph has oversized snarls"
nvg autoindex -p auto -w lr-giraffe --gfa graphs/chain-clip.gfa
is "$(echo $?)" 0 "autoindexing successfully completes indexing of graph with oversized snarls"

rm auto.* log.txt

nvg autoindex -p auto -w sampling --gfa graphs/gfa_with_reference.gfa
is "$(echo $?)" 0 "autoindexing successfully completes indexing for haplotype sampling"
nvg haplotypes -i auto.hapl -k haplotype-sampling/HG003.kff -g sampled.gbz auto.gbz
is "$(echo $?)" 0 "haplotype sampling indexes can be used"

rm auto.* sampled.gbz
rm read.fq read.gam
