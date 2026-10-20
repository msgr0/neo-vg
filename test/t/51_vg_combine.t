#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 10

nvg construct -r small/x.fa -v small/x.vcf.gz  >  x.vg
nvg construct -r small/x.fa -v small/x.vcf.gz | nvg paths -d -v - > y.vg
nvg construct -r small/x.fa -v small/x.vcf.gz | nvg paths -d -v - > z.vg

nvg combine x.vg y.vg z.vg > xyz.vg

is $(nvg stats -z xyz.vg | grep nodes | awk '{print $2}') 645 "combined graph as correct total nodes"
is $(nvg stats -z xyz.vg | grep edges | awk '{print $2}') 888 "combined graph as correct total edges"
nvg paths -Ev x.vg > x.paths
nvg paths -Ev xyz.vg > xyz.paths
diff x.paths xyz.paths
is $? 0  "combined graph has same input path"

rm -f x.vg y.vg z.vg xyz.vg x.paths xyz.paths

nvg construct -r small/x.fa -v small/x.vcf.gz | nvg convert -p - >  x.vg
nvg construct -r small/x.fa -v small/x.vcf.gz | nvg paths -d -v - | nvg convert -a - > y.vg
nvg construct -r small/x.fa -v small/x.vcf.gz | nvg paths -d -v - | nvg convert -p - > z.vg

nvg combine x.vg y.vg z.vg > xyz.vg

is $(nvg stats -z xyz.vg | grep nodes | awk '{print $2}') 645 "combined handle graph as correct total nodes"
is $(nvg stats -z xyz.vg | grep edges | awk '{print $2}') 888 "combined handle graph as correct total edges"
nvg paths -Ev x.vg > x.paths
nvg paths -Ev xyz.vg > xyz.paths
diff x.paths xyz.paths
is $? 0  "combined handle graph has same input path"

rm -f x.vg y.vg z.vg xyz.vg x.paths xyz.paths

nvg construct -r small/x.fa -v small/x.vcf.gz | nvg convert -p - >  x.vg
nvg construct -r small/x.fa -v small/x.vcf.gz | nvg convert -a - > y.vg
nvg construct -r small/x.fa -v small/x.vcf.gz | nvg convert -v - > z.vg

nvg combine -p x.vg y.vg z.vg > xyz.vg

is $(nvg stats -z xyz.vg | grep nodes | awk '{print $2}') 645 "path combined graph as correct total nodes"
is $(nvg stats -z xyz.vg | grep edges | awk '{print $2}') 890 "path combined path graph as correct total edges"
is $(nvg paths -Ev xyz.vg | wc -l) 1 "path combined graph as correct number of paths"
is $(nvg paths -Ev xyz.vg | awk '{print $2}') 3003 "combined path has correct size"

rm -f x.vg y.vg z.vg xyz.vg


