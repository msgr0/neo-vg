#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 33

is $(nvg construct -m 1000 -r small/x.fa -v small/x.vcf.gz | nvg stats -z - | grep nodes | cut -f 2) 210 "construction produces the right number of nodes"

is $(nvg construct -m 1000 -r small/x.fa -v small/x.vcf.gz | nvg stats -z - | grep edges | cut -f 2) 291 "construction produces the right number of edges"

is $(nvg construct -r small/x.fa --rename chrX=x -R chrX:1-2 | nvg stats -l - | cut -f 2 ) 2 "construction obeys rename and region options"

nvg construct -m 1000 -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >z.vg
is $? 0 "construction of a 1 megabase graph from the 1000 Genomes succeeds"

nodes=$(nvg stats -z z.vg | head -1 | cut -f 2)
is $nodes 84559 "the 1mb graph has the expected number of nodes"

edges=$(nvg stats -z z.vg | tail -1 | cut -f 2)
is $edges 115375 "the 1mb graph has the expected number of edges"

# Our 1000 Genomes subset file has changed variant start positions but not stored end positions.
nvg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz -S >/dev/null
is "${?}" "0" "nvg construct can skip self-inconsistent structural variants"

rm -f z.vg

is $(nvg construct -r 1mb1kgp/z.fa | nvg view -j - | jq -c '.node[] | select((.sequence | length) >= 1024)' | wc -l) 0 "node size is manageable by default"

nvg construct -m 1000 -r complex/c.fa -v complex/c.vcf.gz >c.vg
is $? 0 "construction of a very complex region succeeds"

nodes=$(nvg stats -z c.vg | head -1 | cut -f 2)
is $nodes 71 "the complex graph has the expected number of nodes"

edges=$(nvg stats -z c.vg | tail -1 | cut -f 2)
is $edges 116 "the complex graph has the expected number of edges"

rm -f c.vg

order_a=$(nvg construct -r order/n.fa -v order/x.vcf.gz | md5sum | cut -f 1 -d\ )
order_b=$(nvg construct -r order/n.fa -v order/y.vcf.gz | md5sum | cut -f 1 -d\ )

is $order_a $order_b "the ordering of variants at the same position has no effect on the resulting graph"

nvg construct -r order/n.fa -v order/z.vcf.gz -R n:47-73 >/dev/null
is $? 0 "construction does not fail when the first position in the VCF is repeated and has an indel"

x1=$(for i in $(seq 100); do size=$(shuf -i 1-100 -n 1); threads=1; nvg construct -r small/x.fa -v small/x.vcf.gz -z $size -t $threads | nvg view -g - | sort -n -k 2 | md5sum; done | sort | uniq | wc -l)

is $x1 1 "the size of the regions used in construction has no effect on the graph"

x2=$(for i in $(seq 100); do
    size=10;
    threads=$(shuf -i 1-$(nproc) -n 1);
    nvg construct -r small/x.fa -v small/x.vcf.gz -z $size -t $threads | nvg view -g - | sort -n -k 2 | md5sum;
    done | sort | uniq | wc -l)

is $x2 1 "the number of threads used in construction has no effect on the graph"

x3=$(for i in $(seq 100); do
    size=$(shuf -i 1-$(nproc) -n 1);
    threads=$(shuf -i 1-100 -n 1);
    nvg construct -r small/x.fa -v small/x.vcf.gz -z $size -t $threads | nvg view -g - | sort -n -k 2 | md5sum;
    done | sort | uniq | wc -l)

is $x3 1 "the number of threads and regions used in construction has no effect on the graph"

nvg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz -R z:10-20 >/dev/null
is $? 0 "construction of a graph with two head nodes succeeds"

# in case there were failures in topological sort
rm -f fail.vg

# check that we produce a full graph

refbp=$(tail -n+2 small/x.fa | tr -d '\n' | wc -c)
variantbp=$(zcat < small/x.vcf.gz | grep -v "#" | cut -f 5,4 \
    | awk '{ x=length($2)-length($1); if (x > 0) { print x; } else if (x == 0) { print length($2); } }' \
        | awk '{ sum += $1 } END { print sum }')

graphbp=$(nvg construct -r small/x.fa -v small/x.vcf.gz | nvg stats -l - | cut -f 2)

is "$graphbp" $(echo "$refbp + $variantbp" | bc) "the graph contains all the sequence in the reference and VCF"

is $(for i in $(seq 100); do nvg construct -r small/x.fa -v small/x.vcf.gz -m $i | nvg stats -l - | cut -f 2; done | sort | uniq | wc -l) 1 "varying the max node size does not affect graph length"

max_node_size=$(nvg construct -r small/x.fa -v small/x.vcf.gz -m 12 | nvg view -g - | grep ^S | cut -f 3 | awk '{ print length($1) }' | sort -n | tail -1)

is $max_node_size 12 "nodes are correctly capped in size"

## Check the length of the longest node
is $(nvg construct -m 1000 -R z:10000-20000 -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz | nvg view - | grep "^S" | awk '{ print length($3); }' | sort -n | tail -1) 241 "-R --region flag is respected"

nvg construct -r small/x.fa >/dev/null
is $? 0 "nvg construct does not require a vcf"

(( short_enough=$(nvg construct -r small/x.fa -m 50 | nvg view - | grep "^S" | cut -f3 | wc -l) < 51 ? 1 : 0 ))
is $short_enough 1 "nvg construct respects node size limit"

is $(nvg construct -CR 'gi|568815592:29791752-29792749' -r GRCh38_alts/FASTA/HLA/V-352962.fa | nvg view - | grep TCTAGAAGAGTCCACGGGGACAGGTAAG | wc -l) 1 "--region can be interpreted to be a reference sequence (and not parsed as a region spec)"

is "$(nvg construct -r sv/x.fa -v sv/x.inv.vcf -S | nvg view - | sort | md5sum | cut -f 1 -d\ )" "$(cat sv/x.inv.gfa | sort | md5sum | cut -f1 -d\ )" "vg constructs the correct graph for inversions"

is "$(nvg construct -r tiny/tiny.fa -v tiny/dots.vcf | nvg mod -u - | nvg stats -z - | grep nodes | cut -f2)" "1" "nvg construct skips variants with . ALTs"

nvg construct -r small/x.fa -v small/x.v4.3.vcf.gz >x.vg 2>log.txt
is $? 0  "VCF with * alleles is accepted"
grep "skipping variant" log.txt >/dev/null
is $? 0  "VCF record with * allele is rejected"
rm -f x.vg log.txt

nvg construct -r tiny/ambiguous.fa >tiny.vg
is $? 0 "Reference with ambiguity codes has them coerced to Ns"
is "$(nvg view -j tiny.vg | jq -r '.node[].sequence' | tr -d 'ACGT\n' | wc -c)" "10" "Expected number of Ns are created"
rm -f tiny.vg

nvg construct -r small/x.fa -r small/x.fa > /dev/null
is $? 1 "Names may not be duplicated across files"
sed "s/y/x/" small/xy.fa > xx.fa
samtools faidx xx.fa
# Samtools drops duplicate entries so we need to manually make an index with duplicates.
cat xx.fa.fai xx.fa.fai >xx.fa.fai2
mv xx.fa.fai2 xx.fa.fai
nvg construct -r xx.fa > /dev/null
is $? 1 "Names may not be duplicated within an index"
rm -f xx.fa.fai
nvg construct -r xx.fa > /dev/null
is $? 1 "Names may not be duplicated within an unindexed file"
rm -f xx.fa xx.fa.fai
