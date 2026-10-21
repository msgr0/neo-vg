#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 33

nvg construct -r small/x.fa > j.vg
nvg index -x j.xg j.vg
nvg construct -r small/x.fa -v small/x.vcf.gz > x.vg
nvg index -k 11 -g x.gcsa -x x.xg x.vg

#samtools view -f 0 small/x.bam | awk '$6 ~ /100M/' > p.bam
#samtools view -f 0 small/x.bam | awk '$6 !~ /100M/' > pm.bam
#samtools view -O bam -f 16 small/x.bam > small/i.bam
#samtools view -f 16 small/x.bam | awk '$6 ~ /100M/' > i.bam
#samtools view -f 16 small/x.bam | awk '$6 !~ /100M/' > im.bam

is $(nvg inject -x x.xg small/x.bam | nvg view -a - | wc -l) \
    1000 "reads are generated"

is $(nvg inject -x x.xg small/x.bam | nvg surject -x x.xg -t 1 - | nvg view -a - | wc -l) \
    1000 "nvg inject works for all reads included in the original bam"

is "$(samtools view small/x.bam | cut -f 1 | sort)" "$(nvg inject -x x.xg small/x.bam | nvg surject -p x -x x.xg -t 1 -s - | samtools view -S - | cut -f 1 | sort)" \
    "nvg inject retains read names"

is "$(samtools view small/x.bam | cut -f 4)" "$(nvg inject -x j.xg small/x.bam | nvg surject -p x -x j.xg -t 1 -s - | samtools view -S - | cut -f 4 | sort -n)" \
    "nvg inject works perfectly for the position of alignment in the simple graph"

is "$(samtools view small/x.bam | cut -f 4)" "$(nvg inject -x x.xg small/x.bam | nvg surject -p x -x x.xg -t 1 -s - | samtools view -S - | cut -f 4 | sort -n)" \
    "nvg inject works perfectly for the position of alignment"

is "$(samtools view small/x.bam | cut -f 5 | grep 60 | wc -l | tr -d " ")" "$(nvg inject -x x.xg small/x.bam | nvg surject -p x -x x.xg -t 1 -s - | samtools view -S - | cut -f 5 | grep 60 | wc -l | tr -d " ")" \
    "nvg inject works perfectly for the mapping quality of alignment"

is "$(samtools view small/x.bam | sort -k 1 | cut -f 6)" "$(nvg inject -x j.xg small/x.bam | nvg surject -p x -x j.xg -t 1 -s - | samtools view -S - | sort -k 1 | cut -f 6)" \
    "nvg inject works perfectly for the cigar of alignment in the simple graph"

is "$(samtools view small/x.bam | sort -k 1 | cut -f 6)" "$(nvg inject -x x.xg small/x.bam | nvg surject -p x -x x.xg -t 1 -s - | samtools view -S - | sort -k 1 | cut -f 6)" \
    "nvg inject works perfectly for the cigar of alignment"

is "$(nvg inject -x x.xg small/i.bam | nvg view -a - | wc -l | tr -d " ")" 470 "nvg inject preserves all reads"

is "$(nvg inject -x x.xg small/i.bam | nvg view -a - | jq .path.mapping[0].position.is_reverse | grep -v true | wc -l | tr -d " ")" \
    0 "nvg inject works perfectly for the reads flagged as is_reverse"

cat <(samtools view -H small/x.bam) <(printf "name\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTACGT\tHHHHHHHHHHHH\n") > unmapped.sam
is "$(nvg inject -x x.xg unmapped.sam | nvg view -aj - | grep "path" | wc -l | tr -d " ")" 0 "nvg inject does not make an alignment for an umapped read"
is "$(echo $?)" 0 "nvg inject does not crash on an unmapped read"

is $(nvg inject -x x.xg small/x.bam -o GAF | grep -v "^@" | wc -l) \
    1000 "nvg inject supports GAF output"

samtools cat small/x.bam small/i.bam > all.bam
is "$(nvg inject -x x.xg all.bam | nvg validate -A -c -a - x.xg 2>&1 | grep invalid | wc -l | tr -d " ")" 0 "nvg inject produces valid alignments of the read and reference"

nvg inject -x x.xg small/atstart.sam >/dev/null 2>log.txt
is "${?}" 0 "nvg inject can inject a read abutting the start of the contig"
is "$(cat log.txt | wc -l | tr -d " ")" "0" "nvg inject does not warn about a read abutting the start of the contig"

# TODO: We can't see a 0-position-in-the-SAM past-start mapped read properly,
# because htslib prints a warning and fixes it to unmapped for us.
nvg inject -x x.xg small/paststart.sam >/dev/null 2>log.txt
is "$(grep '\[W' log.txt | wc -l | tr -d " ")" "1" "nvg inject warns about a read mapping extending past the start of the contig"

nvg inject -x x.xg small/atend.sam >/dev/null 2>log.txt
is "${?}" 0 "nvg inject can inject a read abutting the end of the contig"
is "$(cat log.txt | wc -l | tr -d " ")" "0" "nvg inject does not warn about a read abutting the end of the contig"

nvg inject -x x.xg small/pastend.sam >/dev/null 2>log.txt
is "${?}" 1 "nvg inject aborts when given a read extending past the end of the contig"
is "$(grep error log.txt | grep 1001 | grep 1002 | wc -l | tr -d " ")" "1" "nvg inject reports a useful error message about a read extending past the end of the contig"

nvg inject -x x.xg small/unmapped.sam | nvg stats -a - >stats.txt
is "$(grep "Total alignments: 4" stats.txt | wc -l | tr -d " ")" "1" "Injecting reads flagged as unmapped produces alignment records"
is "$(grep "Total aligned: 0" stats.txt | wc -l | tr -d " ")" "1" "Injecting reads flagged as unmapped produces unaligned alignment records"

nvg inject -x x.xg small/bad_contig.sam >/dev/null 2>log.txt
is "${?}" 1 "nvg inject aborts when given a read on a contig not in the graph"
is "$(grep error log.txt | grep 'not present in graph' | wc -l | tr -d " ")" "1" \
"nvg inject reports a useful error message about a read on a contig not in the graph"

nvg inject -t 2 -x x.xg small/bad_contig.sam >/dev/null 2>log.txt
is "${?}" 1 "nvg inject aborts when given a read on a contig not in the graph (parallel)"
is "$(grep error log.txt | grep 'not present in graph' | wc -l | tr -d " ")" "1" \
"nvg inject reports a useful error message about a read on a contig not in the graph (parallel)"

nvg inject -a -x x.xg small/bad_contig.sam > bad_contig.gam
is "${?}" 0 "nvg inject -a allows a read on a contig not in the graph"
nvg stats -a bad_contig.gam >stats.txt
is "$(grep "Total alignments: 2" stats.txt | wc -l | tr -d " ")" "1" "Injecting reads on missing contigs still produces alignment records"
is "$(grep "Total aligned: 1" stats.txt | wc -l | tr -d " ")" "1" "Injecting reads on missing contigs produces unmapped alignment records"

nvg inject -t 2 -a -x x.xg small/bad_contig.sam > bad_contig.gam
is "${?}" 0 "nvg inject -a allows a read on a contig not in the graph (parallel)"
nvg stats -a bad_contig.gam >stats.txt
is "$(grep "Total alignments: 2" stats.txt | wc -l | tr -d " ")" "1" "Injecting reads on missing contigs still produces alignment records (parallel)"
is "$(grep "Total aligned: 1" stats.txt | wc -l | tr -d " ")" "1" "Injecting reads on missing contigs produces unmapped alignment records (parallel)"

rm j.vg j.xg x.vg x.gcsa x.gcsa.lcp x.xg unmapped.sam all.bam log.txt stats.txt bad_contig.gam
