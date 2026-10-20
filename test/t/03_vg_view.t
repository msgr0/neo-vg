#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 24

is $(nvg construct -m 1000 -r small/x.fa -v small/x.vcf.gz | nvg view -d - | wc -l) 505 "view produces the expected number of lines of dot output"
is $(nvg construct -m 1000 -r small/x.fa -v small/x.vcf.gz | nvg view -g - | wc -l) 503 "view produces the expected number of lines of GFA output"
# This one may throw warnings related to dangling edges because we just take an arbitrary subset of the GFA
is $(nvg construct -r small/x.fa -v small/x.vcf.gz | nvg view - | head | nvg view -Fv - | nvg view - | wc -l) 10 "view converts back and forth between GFA and vg format"

is $(samtools view -u minigiab/NA12878.chr22.tiny.bam | nvg view -bG - | nvg view -a - | wc -l) $(samtools view -u minigiab/NA12878.chr22.tiny.bam | samtools view - | wc -l) "view can convert BAM to GAM"

is "$(samtools view -u minigiab/NA12878.chr22.tiny.bam | nvg view -bG - | nvg view -aj - | jq -c --sort-keys . | sort | md5sum)" "$(samtools view -u minigiab/NA12878.chr22.tiny.bam | nvg view -bG - | nvg view -aj - | nvg view -JGa - | nvg view -aj - | jq -c --sort-keys . | sort | md5sum)" "view can round-trip JSON and GAM"

# Look at only forward strand reads here to avoid having to reverse qualities.
is "$(samtools view -F 16 minigiab/NA12878.chr22.tiny.bam | head -n1 | cut -f11)" "$(samtools view -F 16 -u minigiab/NA12878.chr22.tiny.bam | nvg view -bG - | nvg view -aj - | nvg view -JGa - | nvg view -aj - | head -n1 | jq -r '.quality' | base64 -d | tr '\000-\077' '\041-\150')" "view can round-trip qualities and encodes tham as PHRED+0 base-64-encoded strings in JSON"

# We need to run through GFA because nvg construct doesn't necessarily chunk the
# graph the way nvg view wants to. We also need to treat as vg::VG to preserve ordering.
nvg construct -r small/x.fa -v small/x.vcf.gz | nvg view -Vg - | nvg view -Fv - >x.vg
nvg view -Vj x.vg | jq . | nvg view -Jv - | diff x.vg -
is $? 0 "view can reconstruct a VG graph from JSON"

nvg view -vV x.vg | cmp -s - x.vg
is $? 0 "view can pass through VG when loading as vg::VG"

nvg view x.vg | sort > x.gfa.sorted
nvg view -v x.vg | nvg view - | sort | cmp -s - x.gfa.sorted
is $? 0 "view can pass through semantically identical VG normally"

rm -f x.vg x.gfa.sorted

is $(samtools view -u minigiab/NA12878.chr22.tiny.bam | nvg view -bG - | nvg view -a - | jq .sample_name | grep -v '^"1"$' | wc -l ) 0 "view parses sample names"

is $(nvg view -f ./small/x.fa_1.fastq  ./small/x.fa_2.fastq | nvg view -a - | wc -l) 2000 "view can handle fastq input"

is $(nvg view -Jv ./cyclic/two_node.json | nvg view -j - | jq ".edge | length") 4 "view can translate graphs with 2-node cycles"

is $(nvg view -g ./cyclic/all.vg | tr '\t' ' ' | grep "4 + 4 -" | wc -l) 1 "view outputs properly oriented GFA"

is $(nvg view -d ./cyclic/all.vg | wc -l) 23 "view produces the expected number of lines of dot output from a cyclic graph"

# We need to make a single-chunk graph
nvg construct -r small/x.fa -v small/x.vcf.gz | nvg view -v - >x.vg
is $(cat x.vg x.vg x.vg x.vg | nvg view -c - | wc -l) 4 "streaming JSON output produces the expected number of chunks"

is "$(cat x.vg x.vg | nvg view -vVD - 2>&1 > /dev/null | grep warning | grep -v deprecated | wc -l)" 0 "duplicate warnings can be suppressed when loading as vg::VG"

rm x.vg

nvg view -Fv overlaps/two_snvs_assembly1.gfa >/dev/null 2>errors.txt
is "${?}" "1" "gfa graphs with overlaps are rejected"
is "$(cat errors.txt | fgrep "Input GFA is not acceptable" | wc -l)" "1" "GFA import produces a concise error message when overlaps are present"

nvg view -Fv overlaps/incorrect_overlap.gfa >/dev/null 2>errors.txt
is "$?" "1" "GFA import rejects a GFA file with an overlap that goes beyond its sequences"
is "$(cat errors.txt | fgrep "Input GFA is not acceptable" | wc -l)" "1" "GFA import produces a concise error message in that case"

rm -f errors.txt

nvg paths -M -x graphs/rgfa_with_reference.rgfa | sort > paths.truth.txt
nvg view graphs/rgfa_with_reference.rgfa | nvg paths -M -x - | sort > paths.test.txt
cmp paths.test.txt paths.truth.txt
is "${?}" "0" "nvg view preserves path metadata of rGFA file"

nvg paths -M -x graphs/gfa_with_reference.gfa | sort > paths.truth.txt
nvg view graphs/gfa_with_reference.gfa | nvg paths -M -x - | sort > paths.test.txt
cmp paths.test.txt paths.truth.txt
is "${?}" "0" "nvg view preserves path metadata of GFA file"

rm -f paths.truth.txt paths.test.txt

nvg view graphs/hap_subrange.gfa >/dev/null
is "${?}" "0" "nvg view can emit a GFA file with subranges on a haplotype without explicit phase blocks"

diff <(nvg view tiny/tiny.gfa | grep '^S' | cut -f3 | sort) <(nvg view tiny/tiny.gfaz | grep '^S' | cut -f3 | sort)
is "${?}" "0" "nvg view preserves segment sequences for matching GFA and GFAZ inputs"
