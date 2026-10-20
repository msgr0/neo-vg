#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 195


# Build vg graphs for two chromosomes
nvg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
nvg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
nvg ids -j x.vg y.vg

nvg index -x x.xg x.vg
nvg index -x xy.xg x.vg y.vg
nvg index -x xy-alt.xg -L x.vg y.vg


# Single chromosome: haplotypes
nvg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
is $? 0 "chromosome X GBWT with nvg gbwt"
nvg gbwt -x x.vg -o parse --parse-only -v small/xy2.vcf.gz
is $? 0 "chromosome X VCF parse"
../deps/gbwt/bin/build_gbwt -p -r parse_x > /dev/null 2> /dev/null
is $? 0 "chromosome X GBWT from VCF parse"
cmp x.gbwt parse_x.gbwt
is $? 0 "identical construction results with nvg gbwt and from VCF parse"

# Single chromosome: metadata for haplotypes
is $(nvg gbwt -c x.gbwt) 2 "chromosome X: 2 paths"
is $(nvg gbwt -C x.gbwt) 1 "chromosome X: 1 contig"
is $(nvg gbwt -H x.gbwt) 2 "chromosome X: 2 haplotypes"
is $(nvg gbwt -S x.gbwt) 1 "chromosome X: 1 sample"
is $(nvg gbwt -T x.gbwt | wc -l) 2 "chromosome X: 2 path names"
is $(nvg gbwt -C -L x.gbwt | wc -l) 1 "chromosome X: 1 contig name"
is $(nvg gbwt -S -L x.gbwt | wc -l) 1 "chromosome X: 1 sample name"

rm -f x.gbwt parse parse_x.gbwt
rm -f parse_x parse_x_0_1


# Single chromosome: paths
nvg gbwt -E -o x.ref.gbwt -x x.vg
is $? 0 "chromosome X reference GBWT with nvg gbwt"
is $(nvg gbwt -c x.ref.gbwt) 1 "chromosome X reference: 1 path"

rm -f x.ref.gbwt


# Single chromosome: alignments
nvg paths -x x.vg -X -Q _alt > x.alts.gam
nvg convert -G x.alts.gam x.vg > x.alts.gaf
nvg gbwt -A --num-jobs 1 -o x.alts.gaf.gbwt -x x.vg x.alts.gaf
is $? 0 "chromosome X GAF with nvg gbwt"
nvg gbwt -A --num-jobs 1 --gam-format -o x.alts.gam.gbwt -x x.vg x.alts.gam
is $? 0 "chromosome X GAM with nvg gbwt"
cmp x.alts.gaf.gbwt x.alts.gaf.gbwt
is $? 0 "identical construction results from GAF and GAM"

rm -f x.alts.gam x.alts.gaf
rm -f x.alts.gaf.gbwt x.alts.gam.gbwt


# Graph region: haplotypes
nvg construct -r small/x.fa -v small/x.vcf.gz -a --region x:100-200 > x.part.vg
nvg gbwt -x x.part.vg -o x.part.gbwt --vcf-region x:100-200 -v small/x.vcf.gz 2> log.txt
is $? 0 "chromosome X subgraph GBWT with nvg gbwt"
is "$(cat log.txt | wc -c)" 0 "no warnings about missing variants"

rm -f x.part.vg x.part.gbwt log.txt


# Multiple chromosomes: haplotypes
nvg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
nvg gbwt -x y.vg -o y.gbwt -v small/xy2.vcf.gz
nvg gbwt -m -o xy.merge.gbwt x.gbwt y.gbwt
is $? 0 "insertion merging with multiple chromosomes"
nvg gbwt -f -o xy.fast.gbwt x.gbwt y.gbwt
is $? 0 "fast merging with multiple chromosomes"
cmp xy.merge.gbwt xy.fast.gbwt
is $? 0 "identical merging results with the insertion and fast algorithms"
nvg gbwt -b --merge-jobs 1 -o xy.parallel.gbwt x.gbwt y.gbwt
is $? 0 "parallel merging with multiple chromosomes"
cmp xy.merge.gbwt xy.parallel.gbwt
is $? 0 "identical merging results with the insertion and parallel algorithms"
nvg gbwt -x xy-alt.xg -o xy.direct.gbwt -v small/xy2.vcf.gz
is $? 0 "direct construction with multiple chromosomes and a single VCF"
cmp xy.direct.gbwt xy.merge.gbwt
is $? 0 "identical results with direct construction and merging"
nvg gbwt -x xy-alt.xg -o xy.multi.gbwt -v --inputs-as-jobs small/xy2_x.vcf.gz small/xy2_y.vcf.gz
is $? 0 "direct construction with multiple chromosomes and multiple VCFs"
cmp xy.direct.gbwt xy.multi.gbwt
is $? 0 "identical construction results with a single VCF and multiple VCFs"

# Multiple chromosomes: haplotypes with presets
nvg gbwt -x xy-alt.xg -o xy.1000gp.gbwt --preset 1000gp -v small/xy2.vcf.gz
is $? 0 "construction preset: 1000gp"

# Multiple chromosomes: metadata for haplotypes
is $(nvg gbwt -c xy.merge.gbwt) 4 "multiple chromosomes: 4 paths"
is $(nvg gbwt -C xy.merge.gbwt) 2 "multiple chromosomes: 2 contigs"
is $(nvg gbwt -H xy.merge.gbwt) 2 "multiple chromosomes: 2 haplotypes"
is $(nvg gbwt -S xy.merge.gbwt) 1 "multiple chromosomes: 1 sample"

rm -f x.gbwt y.gbwt xy.merge.gbwt xy.fast.gbwt xy.parallel.gbwt xy.direct.gbwt xy.multi.gbwt
rm -f xy.1000gp.gbwt


# Multiple chromosomes: alignments
nvg paths -x xy-alt.xg -X -Q _alt > xy.alts.gam
nvg convert -G xy.alts.gam xy.xg > xy.alts.gaf
nvg gbwt -A --num-jobs 1 -o xy.alts.gaf.gbwt -x xy.xg xy.alts.gaf
is $? 0 "multi-chromosome GAF with nvg gbwt"
nvg gbwt -A --num-jobs 1 --gam-format -o xy.alts.gam.gbwt -x xy.xg xy.alts.gam
is $? 0 "multi-chromosome GAM with nvg gbwt"
cmp xy.alts.gaf.gbwt xy.alts.gaf.gbwt
is $? 0 "identical construction results from GAF and GAM"

nvg gbwt -A --num-jobs 2 -o multi.gbwt -x xy.xg xy.alts.gaf
is $? 0 "multi-chromosome GAF with nvg gbwt using multiple jobs"
is $(nvg gbwt -c xy.alts.gaf.gbwt) 58 "single job: 58 paths"
is $(nvg gbwt -c multi.gbwt) 58 "multiple jobs: 58 paths"

rm -f xy.alts.gam xy.alts.gaf
rm -f xy.alts.gaf.gbwt xy.alts.gam.gbwt multi.gbwt


# Multiple chromosomes: paths as contigs
nvg gbwt -E -o xy.contigs.gbwt -x xy.xg
is $? 0 "paths as contigs with nvg gbwt"
is $(nvg gbwt -c xy.contigs.gbwt) 2 "paths as contigs: 2 paths"
is $(nvg gbwt -C xy.contigs.gbwt) 2 "paths as contigs: 2 contigs"
is $(nvg gbwt -H xy.contigs.gbwt) 1 "paths as contigs: 1 haplotype"
is $(nvg gbwt -S xy.contigs.gbwt) 1 "paths as contigs: 1 sample"

rm -f xy.contigs.gbwt


# Build an r-index
# TODO: More tests once we can use the r-index for something
nvg gbwt -x xy-alt.xg -o xy.gbwt -r xy.ri --num-threads 1 -v small/xy2.vcf.gz
is $? 0 "r-index construction"

rm -f xy.gbwt xy.ri


# Non-empty and empty GBWTs
nvg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
../deps/gbwt/bin/build_gbwt -e empty > /dev/null

# Normal merging, x + empty
nvg gbwt -m -o x2.gbwt x.gbwt empty.gbwt
is $? 0 "normal merging: non-empty + empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

# Normal merging, empty + x; silence the warning about the merging order
nvg gbwt -m -o x2.gbwt empty.gbwt x.gbwt 2> /dev/null
is $? 0 "normal merging: empty + non-empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

# Fast merging, x + empty
nvg gbwt -f -o x2.gbwt x.gbwt empty.gbwt
is $? 0 "fast merging: non-empty + empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

# Fast merging, empty + x
nvg gbwt -f -o x2.gbwt empty.gbwt x.gbwt
is $? 0 "fast merging: empty + non-empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

# Parallel merging, x + empty
nvg gbwt -b --merge-jobs 1 -o x2.gbwt x.gbwt empty.gbwt
is $? 0 "parallel merging: non-empty + empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

# Parallel merging, empty + x; silence the warning about the merging order
nvg gbwt -b --merge-jobs 1 -o x2.gbwt empty.gbwt x.gbwt 2> /dev/null
is $? 0 "parallel merging: empty + non-empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

rm -f x.gbwt empty.gbwt x2.gbwt


# Build a GBWT with both haplotypes and embedded paths
nvg gbwt -x xy-alt.xg -o xy.gbwt -v small/xy2.vcf.gz
nvg gbwt -E -o xy.ref.gbwt -x xy.xg
nvg gbwt -m -o xy.both.gbwt xy.gbwt xy.ref.gbwt
is $(nvg gbwt -c xy.both.gbwt) 6 "haplotypes and paths: 6 paths"

# Remove the reference sample that GBWTs use for paths
nvg gbwt -R _gbwt_ref -o xy.removed.gbwt xy.both.gbwt
is $? 0 "samples can be removed from a GBWT index"
is $(nvg gbwt -c xy.removed.gbwt) 4 "haplotypes only: 4 paths"

rm -f xy.gbwt xy.ref.gbwt xy.both.gbwt xy.removed.gbwt


# Build a three-sample GBWT from a simple GFA
nvg gbwt -o all.gbwt -G graphs/three_samples.gfa
is $(nvg gbwt -c all.gbwt) 12 "all samples: 12 paths"
is $(nvg gbwt -H all.gbwt) 6 "all samples: 6 haplotypes"

# Remove samples 1 and 3
nvg gbwt -R sample1 -R sample3 -o removed.gbwt all.gbwt
is $? 0 "multiple samples can be removed from a GBWT index"
is $(nvg gbwt -c removed.gbwt) 4 "sample 2: 4 paths"
is $(nvg gbwt -H removed.gbwt) 2 "sample 2: 2 haplotypes"

rm -f all.gbwt removed.gbwt


# Extract paths from GBWT
nvg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
nvg gbwt -e x.extract x.gbwt
is $? 0 "paths can be extracted from GBWT"
is $(cat x.extract | wc -c) 121 "correct size for the paths file"

rm -f x.extract


# Build and serialize GBZ from an existing GBWT
nvg gbwt -x x.vg -g x.gbz x.gbwt
is $? 0 "GBZ construction from GBWT"
is $(md5sum x.gbz | cut -f 1 -d\ ) 15bdc4856e19145832885e32e8413321 "GBZ was serialized correctly"
is "$(nvg describe x.gbz | grep -c 'pggname =')" 1 "GBZ contains graph name"

# Try unsetting and setting the name
nvg gbwt --unset-pggname -Z x.gbz -g x.gbz
is $? 0 "Graph name can be unset"
is "$(nvg describe x.gbz | grep -c 'pggname =')" 0 "GBZ does not contain graph name after unsetting"
nvg gbwt --set-pggname -Z x.gbz -g x.gbz
is $? 0 "Graph name can be set"
is "$(nvg describe x.gbz | grep -c 'pggname =')" 1 "GBZ contains graph name after setting"

# Set graph name and a subgraph relationship
nvg gbwt -g supergraph.gbz -G graphs/components_walks.gfa
is $? 0 "GBZ construction for supergraph"
nvg chunk --gbz --contig A -x supergraph.gbz
is $? 0 "Graph component extraction"
nvg gbwt --subgraph-of supergraph.gbz -g subgraph.gbz -Z chunk_0_component_0.gbz
is $? 0 "Subgraph relationship can be set"
is "$(nvg describe subgraph.gbz | grep -c 'subgraph =')" 1 "GBZ contains subgraph relationship"
rm -f supergraph.gbz chunk_0_component_0.gbz subgraph.gbz

# Build and serialize GBZ from VCF
nvg gbwt -x x.vg -g x2.gbz -v small/xy2.vcf.gz
is $? 0 "GBZ construction from VCF"
cmp x.gbz x2.gbz
is $? 0 "Identical construction results from GBWT and VCF"

# Extract GBWT from GBZ
nvg gbwt -o extracted.gbwt -Z x.gbz
is $? 0 "GBWT extraction from GBZ"
cmp x.gbwt extracted.gbwt
is $? 0 "Identical GBWT indexes"

rm -f x.gbwt x.gbz x2.gbz
rm -f extracted.gbwt


# Build GBZ from a 16-path cover
nvg gbwt -P -n 16 -x xy.xg -g xy.cover.gbz
is $? 0 "Path cover GBZ construction"
is $(md5sum xy.cover.gbz | cut -f 1 -d\ ) 9a22d6eb0e6dbdc308279aeb10f655e0 "GBZ was serialized correctly"
is $(nvg gbwt -c -Z xy.cover.gbz) 32 "path cover: 32 paths"
is $(nvg gbwt -C -Z xy.cover.gbz) 2 "path cover: 2 contigs"
is $(nvg gbwt -H -Z xy.cover.gbz) 16 "path cover: 16 haplotypes"
is $(nvg gbwt -S -Z xy.cover.gbz) 16 "path cover: 16 samples"
is "$(nvg describe xy.cover.gbz | grep -c 'pggname =')" 1 "GBZ contains graph name"

rm -f xy.cover.gbz

# Build GBZ from a 16-path cover, passing through named paths
nvg gbwt -P -n 16 -x xy.xg -g xy.cover.gbz --pass-paths
is $? 0 "Path cover GBZ construction"
is $(md5sum xy.cover.gbz | cut -f 1 -d\ ) e14c5b97d84ab28c2005b8ac9c9a0fad "GBZ was serialized correctly"
is $(nvg gbwt -c -Z xy.cover.gbz) 34 "path cover w/ reference paths: 34 paths"
is $(nvg gbwt -C -Z xy.cover.gbz) 2 "path cover w/ reference paths: 2 contigs"
is $(nvg gbwt -H -Z xy.cover.gbz) 17 "path cover w/ reference paths: 17 haplotypes"
is $(nvg gbwt -S -Z xy.cover.gbz) 17 "path cover w/ reference paths: 17 samples"

rm -f xy.cover.gbz

# Build GBZ from 16 paths of local haplotypes with a single job
nvg gbwt -x xy-alt.xg -g xy.local.gbz -l -n 16 --num-jobs 1 -v small/xy2.vcf.gz
is $? 0 "Local haplotypes GBZ construction"
is $(md5sum xy.local.gbz | cut -f 1 -d\ ) b8f55c89391700c27be055cbd5fe96f5 "GBZ was serialized correctly"
is $(nvg gbwt -c -Z xy.local.gbz) 32 "local haplotypes: 32 paths"
is $(nvg gbwt -C -Z xy.local.gbz) 2 "local haplotypes: 2 contigs"
is $(nvg gbwt -H -Z xy.local.gbz) 16 "local haplotypes: 16 haplotypes"
is $(nvg gbwt -S -Z xy.local.gbz) 16 "local haplotypes: 16 samples"

rm -f xy.local.gbz

# As above, but with two parallel jobs
nvg gbwt -x xy-alt.xg -g xy.local.gbz -l -n 16 --num-jobs 2 -v small/xy2.vcf.gz
is $? 0 "Local haplotypes GBZ construction (two jobs)"
is $(md5sum xy.local.gbz | cut -f 1 -d\ ) b8f55c89391700c27be055cbd5fe96f5 "GBZ was serialized correctly"

rm -f xy.local.gbz

# Build GBZ from 16 paths of local haplotypes, passing through named paths
nvg gbwt -x xy-alt.xg -g xy.local.gbz -l -n 16 -v small/xy2.vcf.gz --pass-paths
is $? 0 "Local haplotypes GBZ construction"
is $(md5sum xy.local.gbz | cut -f 1 -d\ ) cfe5a617ce959b97f7b52b4f5fe9dadc "GBZ was serialized correctly"
is $(nvg gbwt -c -Z xy.local.gbz) 34 "local haplotypes w/ reference paths: 34 paths"
is $(nvg gbwt -C -Z xy.local.gbz) 2 "local haplotypes w/ reference paths: 2 contigs"
is $(nvg gbwt -H -Z xy.local.gbz) 17 "local haplotypes w/ reference paths: 17 haplotypes"
is $(nvg gbwt -S -Z xy.local.gbz) 17 "local haplotypes w/ reference paths: 17 samples"

rm -f xy.local.gbz

# Build GBZ from a GFA and then build a local haplotype cover with reference paths from the GBZ
nvg gbwt -G haplotype-sampling/micb-kir3dl1.gfa -g large.gbz
nvg gbwt -Z large.gbz -l -n 16 --pass-paths -o large.local.gbwt
is $? 0 "Local haplotypes with reference paths from a larger GBZ"
is $(nvg gbwt -c large.local.gbwt) 36 "local haplotypes w/ reference paths: 36 paths"
is $(nvg gbwt -C large.local.gbwt) 2 "local haplotypes w/ reference paths: 2 contigs"
is $(nvg gbwt -H large.local.gbwt) 18 "local haplotypes w/ reference paths: 18 haplotypes"
is $(nvg gbwt -S large.local.gbwt) 18 "local haplotypes w/ reference paths: 18 samples"
is $(nvg gbwt --tags large.local.gbwt | grep -c reference_samples) 1 "local haplotypes w/ reference paths: reference_samples set"

# Now the same but with a specified reference sample
nvg gbwt -Z large.gbz -l -n 16 --pass-paths --set-reference GRCh38 -o large.local.GCRh38.gbwt
is $? 0 "Local haplotypes with a specified reference sample"
is $(nvg gbwt -c large.local.GCRh38.gbwt) 34 "local haplotypes w/ reference paths: 34 paths"
is $(nvg gbwt -C large.local.GCRh38.gbwt) 2 "local haplotypes w/ reference paths: 2 contigs"
is $(nvg gbwt -H large.local.GCRh38.gbwt) 17 "local haplotypes w/ reference paths: 17 haplotypes"
is $(nvg gbwt -S large.local.GCRh38.gbwt) 17 "local haplotypes w/ reference paths: 17 samples"
is $(nvg gbwt --tags large.local.GCRh38.gbwt | grep -c reference_samples) 1 "local haplotypes w/ reference paths: reference_samples set"

rm -f large.gbz large.local.gbwt large.local.GCRh38.gbwt


# Build GBZ from an augmented GBWT
nvg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
nvg gbwt -a -n 16 -x xy.xg -g augmented.gbz x.gbwt
is $? 0 "Augmented GBZ construction"
is $(md5sum augmented.gbz | cut -f 1 -d\ ) bcfb39767a4605320d6878d2551c1895 "GBZ was serialized correctly"
is $(nvg gbwt -c -Z augmented.gbz) 18 "augmented: 18 paths"
is $(nvg gbwt -C -Z augmented.gbz) 2 "augmented: 2 contigs"
is $(nvg gbwt -H -Z augmented.gbz) 18 "augmented: 18 haplotypes"
is $(nvg gbwt -S -Z augmented.gbz) 17 "augmented: 17 samples"

rm -f x.gbwt augmented.gbz


# Remove the graphs
rm -f x.vg y.vg x.xg xy.xg xy-alt.xg


# Build GBWT from GFA
nvg gbwt -o gfa.gbwt -G graphs/components_walks.gfa
is $? 0 "GBWT construction from GFA"
is $(md5sum gfa.gbwt | cut -f 1 -d\ ) 44c27c37c7af6911c26aea2a41008460 "GBWT was serialized correctly"
is $(nvg gbwt -c gfa.gbwt) 4 "gfa: 4 paths"
is $(nvg gbwt -C gfa.gbwt) 2 "gfa: 2 contigs"
is $(nvg gbwt -H gfa.gbwt) 2 "gfa: 2 haplotypes"
is $(nvg gbwt -S gfa.gbwt) 1 "gfa: 1 sample"

# Build GBZ from GFA
nvg gbwt -g gfa2.gbz -G graphs/components_walks.gfa
is $? 0 "GBZ construction from GFA"
is $(md5sum gfa2.gbz | cut -f 1 -d\ ) 74d119e982352d10a3bd29f23db0ccc3 "GBZ was serialized correctly"
nvg gbwt -Z gfa2.gbz -o /dev/null --translation gfa2.trans
is $(wc -l < gfa2.trans) 0 "no chopping: 0 translations"
is "$(nvg describe gfa2.gbz | grep -c 'pggname =')" 1 "GBZ contains graph name"

# Build GBZ from grammar-compressed GFA
nvg gbwt -g gfa_grammar.gbz -G graphs/components_walks_compressed.gfa
is $? 0 "GBZ construction from grammar-compressed GFA"
cmp gfa_grammar.gbz gfa2.gbz
is $? 0 "grammar-compressed GFA produces identical GBZ"

# Build GBZ from GFA with node chopping
nvg gbwt -g chopping.gbz --max-node 2 --translation chopping.trans -G graphs/chopping_walks.gfa
is $? 0 "GBZ construction from GFA with chopping"
is $(nvg gbwt -c -Z chopping.gbz) 3 "chopping: 3 paths"
is $(nvg gbwt -C -Z chopping.gbz) 1 "chopping: 1 contig"
is $(nvg gbwt -H -Z chopping.gbz) 3 "chopping: 3 haplotypes"
is $(nvg gbwt -S -Z chopping.gbz) 2 "chopping: 2 samples"
is $(wc -l < chopping.trans) 9 "chopping: 9 translations"
nvg gbwt -Z chopping.gbz -o /dev/null --translation from_gbz.trans
is $? 0 "Translation can be extracted from GBZ"
is $(wc -l < from_gbz.trans) 8 "from GBZ: 8 translations"
is "$(nvg describe chopping.gbz | grep -c 'pggname =')" 1 "GBZ contains graph name"
is "$(nvg describe chopping.gbz | grep -c 'translation =')" 1 "GBZ contains a translation relationship"

# Build GBZ from GFA with both paths and walks
nvg gbwt -g ref_paths.gbz --translation ref_paths.trans -G graphs/components_paths_walks.gfa
is $? 0 "GBZ construction from GFA with reference paths"
is $(nvg gbwt -c -Z ref_paths.gbz) 6 "ref paths: 6 paths"
is $(nvg gbwt -C -Z ref_paths.gbz) 2 "ref paths: 2 contigs"
is $(nvg gbwt -H -Z ref_paths.gbz) 3 "ref paths: 3 haplotypes"
is $(nvg gbwt -S -Z ref_paths.gbz) 2 "ref paths: 2 samples"
is $(wc -l < ref_paths.trans) 0 "ref paths: 0 translations"

rm -f gfa.gbwt
rm -f gfa2.trans gfa2.gbz gfa_grammar.gbz
rm -f ref_paths.gbz ref_paths.trans
rm -f chopping.gbz chopping.trans from_gbz.trans

# Build a GBZ from a graph with a reference
nvg gbwt -g gfa.gbz -G graphs/gfa_with_reference.gfa
is $? 0 "GBZ construction from GFA with reference"
nvg gbwt -Z gfa.gbz --tags >tags.tsv
is $? 0 "GBZ GBWT tag extraction works"
is "$(grep reference_samples tags.tsv | cut -f2 | tr ' ' '\\n' | sort | tr '\\n' ' ')" "GRCh37 GRCh38" "GBWT tags contain the correct reference samples"
nvg gbwt -g gfa2.gbz -Z gfa.gbz --set-tag "reference_samples=GRCh38 CHM13"
is $? 0 "GBZ GBWT tag modification works"
is "$(nvg paths -M -S GRCh37 -x gfa2.gbz | grep -v "^#" | cut -f2 | grep HAPLOTYPE | wc -l)" "1" "Changing reference_samples tag can make a reference a haplotype"
is "$(nvg paths -M -S CHM13 -x gfa2.gbz | grep -v "^#" | cut -f2 | grep REFERENCE | wc -l)" "1" "Changing reference_samples tag can make a haplotype a reference"
nvg gbwt -g gfa2.gbz -Z gfa.gbz --set-tag "reference_samples=GRCh38#1 CHM13" 2>/dev/null
is $? 1 "GBZ GBWT tag modification validation works"
nvg gbwt -g gfa3.gbz --set-reference GRCh37 --set-reference CHM13 -Z gfa2.gbz
is $? 0 "Samples can be direcly set as references"
is "$(nvg gbwt --tags -Z gfa3.gbz | grep reference_samples | cut -f 2)" "CHM13 GRCh37" "Direct reference assignment works"

# Try writing and reading GBZ v1
nvg gbwt -Z gfa.gbz -g gfa_version.gbz --gbz-version 1
is $? 0 "GBZ version 1 output"
is "$(nvg describe gfa_version.gbz | grep -A 1 '^GBZ header' | grep -c 'Version 1')" "1" "Correct version for GBZ v1 output"
nvg gbwt -Z gfa_version.gbz -g gfa_default.gbz
is $? 0 "GBZ v1 can be read and rewritten as the default GBZ version"
cmp gfa.gbz gfa_default.gbz
is $? 0 "GBZ v1 was converted to the default GBZ version without changes"

# Try writing and reading GBZ v2
nvg gbwt -Z gfa.gbz -g gfa_version.gbz --gbz-version 2
is $? 0 "GBZ version 2 output"
is "$(nvg describe gfa_version.gbz | grep -A 1 '^GBZ header' | grep -c 'Version 2')" "1" "Correct version for GBZ v2 output"
nvg gbwt -Z gfa_version.gbz -g gfa_default.gbz
is $? 0 "GBZ v2 can be read and rewritten as the default GBZ version"
cmp gfa.gbz gfa_default.gbz
is $? 0 "GBZ v2 was converted to the default GBZ version without changes"

# Try writing and reading GBZ v3
nvg gbwt -Z gfa.gbz -g gfa_version.gbz --gbz-version 3
is $? 0 "GBZ version 3 output"
is "$(nvg describe gfa_version.gbz | grep -A 1 '^GBZ header' | grep -c 'Version 3')" "1" "Correct version for GBZ v3 output"
nvg gbwt -Z gfa_version.gbz -g gfa_default.gbz
is $? 0 "GBZ v3 can be read and rewritten as the default GBZ version"
cmp gfa.gbz gfa_default.gbz
is $? 0 "GBZ v3 was converted to the default GBZ version without changes"

rm -f gfa.gbz gfa2.gbz gfa3.gbz tags.tsv gfa_version.gbz gfa_default.gbz

# Build a GBZ from a graph with a reference but no haplotype phase number
nvg gbwt -g gfa.gbz -G graphs/gfa_two_part_reference.gfa
is "$(nvg paths -M --reference-paths -x gfa.gbz | grep -v "^#" | cut -f4 | grep NO_HAPLOTYPE | wc -l)" "2" "GBZ can represent reference paths without haplotype numbers"

rm -f gfa.gbz

# Wrapping contigs: append each haplotype's origin fragment onto its last
# fragment, restoring the end-to-start adjacency of a circular sequence.
# components_walks.gfa has two contigs (A, B), each with a single unfragmented
# count==0 haplotype walk, so the origin and the last fragment coincide and the
# path is exactly doubled.
nvg gbwt -g wrap_src.gbz -G graphs/components_walks.gfa
is $? 0 "built a GBZ with haplotype walks for wrapping"

nvg gbwt -Z wrap_src.gbz --wrap-contig A -g wrap_gbz.gbz
is $? 0 "wrapping a contig in a GBZ succeeds"
is "$(nvg paths -E -x wrap_gbz.gbz | awk '$1 == "sample#1#A#0" {print $2}')" "10" "the wrapped haplotype on contig A is doubled (5 -> 10)"
is "$(nvg paths -E -x wrap_gbz.gbz | awk '$1 == "sample#1#B#0" {print $2}')" "5" "haplotypes on the unwrapped contig B are unchanged"
is "$(nvg view -g wrap_gbz.gbz | grep -c "$(printf '^L\t17\t+\t11\t+')")" "1" "the GBZ graph gains the wrap edge from the end of A back to its start"

# The bare GBWT can be wrapped without a graph.
nvg gbwt -o wrap_src.gbwt -G graphs/components_walks.gfa
nvg gbwt wrap_src.gbwt --wrap-contig A -o wrap_bare.gbwt
is $? 0 "wrapping a contig in a bare GBWT succeeds"
is "$(nvg gbwt -M wrap_bare.gbwt)" "$(nvg gbwt -M wrap_src.gbwt)" "wrapping preserves the path/sample/haplotype/contig metadata"

# Wrapping multiple contigs at once.
nvg gbwt -Z wrap_src.gbz --wrap-contig A --wrap-contig B -g wrap_both.gbz
is $? 0 "wrapping multiple contigs succeeds"
is "$(nvg paths -E -x wrap_both.gbz | awk '$1 == "sample#1#A#0" {print $2}')" "10" "contig A is doubled when wrapping both"
is "$(nvg paths -E -x wrap_both.gbz | awk '$1 == "sample#1#B#0" {print $2}')" "10" "contig B is doubled when wrapping both"

# Wrapping requires an output GBWT.
nvg gbwt -Z wrap_src.gbz --wrap-contig A 2> wrap_err.txt
is $? 1 "wrapping without an output fails"
is "$(grep -c "wrapping contigs requires" wrap_err.txt)" "1" "an appropriate error message is printed"

# A fragmented haplotype is wrapped by appending its origin (first) fragment
# onto its last fragment, not by doubling the origin. wrap_fragmented.gfa has a
# single haplotype on contig A split into two fragments (count 0: nodes 11,12,13;
# count 3: nodes 14,15,16). The origin fragment must stay length 3, and the last
# fragment must grow to 6 (its 3 nodes plus the 3 origin nodes).
nvg gbwt -g wrap_frag.gbz -G graphs/wrap_fragmented.gfa
is $? 0 "built a GBZ with a fragmented haplotype for wrapping"
nvg gbwt -Z wrap_frag.gbz --wrap-contig A -g wrap_frag_out.gbz
is $? 0 "wrapping a fragmented haplotype succeeds"
is "$(nvg paths -E -x wrap_frag_out.gbz | awk '$1 == "sample#1#A#0" {print $2}')" "3" "the origin fragment is left unchanged (length 3)"
is "$(nvg paths -E -x wrap_frag_out.gbz | awk '$1 == "sample#1#A#3" {print $2}')" "6" "the origin is appended onto the last fragment (3 -> 6)"
is "$(nvg view -g wrap_frag_out.gbz | grep -c "$(printf '^L\t16\t+\t11\t+')")" "1" "the wrap edge joins the end of the last fragment back to the origin start"

# A haplotype whose origin (count 0) fragment is absent cannot be wrapped.
nvg gbwt -g wrap_nz.gbz -G graphs/wrap_nonzero_offset.gfa
nvg gbwt -Z wrap_nz.gbz --wrap-contig A -g /dev/null 2> wrap_nz_err.txt
is $? 1 "wrapping a contig with no origin fragment fails"
is "$(grep -c "has no origin fragment" wrap_nz_err.txt)" "1" "an appropriate error message is printed"

rm -f wrap_src.gbz wrap_gbz.gbz wrap_both.gbz wrap_src.gbwt wrap_bare.gbwt wrap_nz.gbz wrap_err.txt wrap_nz_err.txt wrap_frag.gbz wrap_frag_out.gbz
