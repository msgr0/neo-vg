#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 7

nvg construct -r small/x.fa -v small/x.vcf.gz >s.vg
nvg index -x s.xg -g s.gcsa s.vg
nvg snarls -T s.vg > s.snarls
nvg index -j s.dist s.vg
nvg sim -n 1000 -l 100 -e 0.01 -i 0.005 -x s.xg -a >s.sim

is $(nvg map -x s.xg -g s.gcsa -G s.sim --surject-to sam | nvg inject -x s.xg - | nvg gamcompare - s.sim | nvg view -a - | wc -l) 1000 "gamcompare completes"

is $(nvg gamcompare --range 10 s.sim s.sim | nvg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l) 1000 "gamcompare says the truth is correctly mapped"

# Map a couple adjacent reads with multi-positioning
nvg map -x s.xg  -g s.gcsa -s "AATCTCTCTGAACTTCAGTTTAATTATC" > read1.gam
nvg annotate -a read1.gam -p -x s.xg > read1.single.gam
nvg annotate -a read1.gam -m -x s.xg > read1.multi.gam
nvg map -x s.xg  -g s.gcsa -s "TCTAATATGGAGATGATACTACTGACAG" > read2.gam
nvg annotate -a read2.gam -p -x s.xg > read2.single.gam
nvg annotate -a read2.gam -m -x s.xg > read2.multi.gam

is "$(nvg gamcompare -r 30 read1.single.gam read2.single.gam | nvg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l | tr -d " ")" "1" "Reads annotated with leftmost position only are close enough at a long distance"
is "$(nvg gamcompare -r 30 -d s.dist read1.gam read2.gam | nvg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l | tr -d " ")" "1" "Reads are close enough at a long distance with a distance index"
is "$(nvg gamcompare -r 10 read1.single.gam read2.single.gam | nvg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l | tr -d " ")" "0" "Reads annotated with leftmost position only are too far apart at a short distance"
is "$(nvg gamcompare -r 10 -d s.dist read1.gam read2.gam | nvg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l | tr -d " ")" "0" "Reads are too far apart at a short distance with a distance index"
is "$(nvg gamcompare -r 10 read1.multi.gam read2.multi.gam | nvg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l | tr -d " ")" "1" "Reads annotated with multiple positions only are close enough at a short distance"

rm -f read1.gam read1.single.gam read1.multi.gam read2.gam read2.single.gam read2.multi.gam 

rm -f s.vg s.xg s.gcsa s.gcsa.lcp s.sim s.snarls s.dist
