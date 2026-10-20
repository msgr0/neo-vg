# file-info

This file contains extra text that will be added to the man pages generated with doc/nvgmanmd.py
The `# description` section is added to the top of the page, and each `# subcommand` section will be added to given subcommand
When adding a new subcommand, add it to the appropriate section(s) in the description

# description

[vg](https://github.com/vgteam/vg) is a toolkit for variation graph data structures, interchange formats, alignment, genotyping, and variant calling methods.

For more in-depth explanations of tools and workflows, see the [general wiki page](https://github.com/vgteam/vg/wiki).

# synopsis

- **Graph construction and indexing**
    See the [wiki page](https://github.com/vgteam/vg/wiki/File-Types) for an overview of file types.
    - [`nvg autoindex`](#autoindex): automatically construct a graph and indexes for a specific workflow (e.g. giraffe, rpvg). [wiki page](https://github.com/vgteam/vg/wiki/Automatic-indexing-for-read-mapping-and-downstream-inference)
    - [`nvg construct`](#construct): manually construct a graph from a reference and variants. [wiki page](https://github.com/vgteam/vg/wiki/Construction)
    - [`nvg gbwt`](#gbwt): manually build and manipulate GBWTs and indexes (GBWTgraph, GBZ, r-index). [wiki page](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand)
    - [`nvg haplotypes`](#haplotypes): haplotype sample a graph. Recommended for mapping with giraffe. [wiki page](https://github.com/vgteam/vg/wiki/Haplotype-Sampling)
    - [`nvg index`](#index): manually build individual indexes (xg, distance, GCSA, etc). [wiki page](https://github.com/vgteam/vg/wiki/Index-Construction) 
    - [`nvg minimizer`](#minimizer): manually build a minimizer index for mapping. 
    - [`nvg rna`](#rna): construct splicing graphs and pantranscriptomes. [wiki page](https://github.com/vgteam/vg/wiki/Transcriptomic-analyses). Also see [rpvg](https://github.com/jonassibbesen/rpvg) 
- **Read mapping**
    - [`nvg giraffe`](#giraffe): fast haplotype-aware short read alignment. [wiki page](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe)
    - [`nvg map`](#map): MEM-based read alignment.
    - [`nvg mpmap`](#mpmap): splice-aware multipath alignment of short reads. [wiki page](https://github.com/vgteam/vg/wiki/Multipath-alignments-and-vg-mpmap)
- **Downstream analyses**
    - [`nvg call`](#call): call or genotype VCF variants. Uses [nvg pack](#pack). [wiki page](https://github.com/vgteam/vg/wiki/SV-Genotyping-and-variant-calling)
    - [`nvg deconstruct`](#deconstruct): create a VCF from variation in the graph. [wiki page](https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct)
    - [`nvg depth`](#depth): calculate path coverage depth or aggregate alignment coverage Also see [nvg pack](#pack)
    - [`nvg pack`](#pack): convert alignments to a compact coverage index. Used with [nvg call](#call)
    - [`nvg primers`](#primers): filter primers based on variation in a graph. [wiki page](https://github.com/vgteam/vg/wiki/Primer-Filter)
- **Working with read alignments**
    - [`nvg augment`](#augment): embed alignments into a graph.
    - [`nvg filter`](#filter): filter alignments by properties.
    - [`nvg gamcompare`](#gamcompare): compare GAM-format alignments to a truth.
    - [`nvg gampcompare`](#gampcompare): compare [GAMP-format alignments](https://github.com/vgteam/vg/wiki/Multipath-alignments-and-the-GAMP-format) to a truth.
    - [`nvg gamsort`](#gamsort): sort a GAM/GAF file or index a sorted GAM file.
    - [`nvg sim`](#sim): simulate reads from a graph. [wiki page](https://github.com/vgteam/vg/wiki/Simulating-reads-with-vg-sim)
- **Graph and read statistics**
    - [`nvg chains`](#chains): get top-level chains from graph. [wiki page](https://github.com/vgteam/vg/wiki/Snarls-and-chains)
    - [`nvg dotplot`](#dotplot): get dotplot matrix for paths in graph.
    - [`nvg filter`](#filter): get stats about alignments (use `--tsv-out`).
    - [`nvg gbwt`](#gbwt): get stats about a GBWT. [wiki page](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand)
    - [`nvg paths`](#paths): get stats about the paths. [wiki page](https://github.com/vgteam/vg/wiki/Path-Metadata-Model)
    - [`nvg snarls`](#snarls): get a list of snarls within the graph. [wiki page](https://github.com/vgteam/vg/wiki/Snarls-and-chains)
    - [`nvg stats`](#stats): get stats about the graph.
- **Manipulating a graph**
    - [`nvg annotate`](#annotate): annotate a graph or alignments.
    - [`nvg circularize`](#circularize): connect head and tail nodes to circularize paths.
    - [`nvg clip`](#clip): remove variation from a graph.
    - [`nvg combine`](#combine): merge graphs into a combined graph. Useful for per-chromosome graphs.
    - [`nvg mask`](#mask): N-mask regions of a graph.
    - [`nvg mod`](#mod): filter, transform, and edit the graph.
    - [`nvg gbwt`](#gbwt): manipulate GBWTs and associated indexes. [wiki page](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand)
    - [`nvg ids`](#ids): manipulate graph node ids.
    - [`nvg paths`](#paths): manipulate paths in a graph. [wiki page](https://github.com/vgteam/vg/wiki/Path-Metadata-Model)
    - [`nvg prune`](#prune): prune the graph for GCSA2 indexing.
    - [`nvg simplify`](#simplify): simplify graph by removing small leaf sites or rare variants.
- **Conversion between formats**
    - [`nvg convert`](#convert): convert between handle graph formats and GFA, and between alignment formats.
    - [`nvg inject`](#inject): project alignments on a linear reference onto a graph (bam/sam/cram->gam/gaf).
    - [`nvg paths`](#paths): extract a fasta from a graph. [wiki page](https://github.com/vgteam/vg/wiki/Linear-references-and-vg)
    - [`nvg surject`](#surject): project alignments on a graph onto a linear reference (gam/gaf->bam/sam/cram).
    - [`nvg vectorize`](#vectorize): export GAM alignments to other vector formats.
    - [`nvg view`](#view): convert between non-handle graph formats and alignment formats (dot, json, turtle...).
- **Subgraph extraction**
    - [`nvg chunk`](#chunk): split a graph and/or alignment into smaller chunks.
    - [`nvg find`](#find): use an index to find nodes, edges, kmers, paths, or positions.
    - [`nvg paths`](#paths): subset a graph to keep or remove particular paths.
    - [`nvg trace`](#trace): extract (sub-)haplotypes from a graph.
- **Extremely specific analyses** (i.e. don't use unless you know what you're doing)
    - [`nvg align`](#align): run highly configurable alignment of a single sequence.
- **Developer tools**
    - [`nvg describe`](#describe): determine file format of input. [wiki page](https://github.com/vgteam/vg/wiki/File-Types)
    - [`nvg help`](#help): display all subcommands.
    - [`nvg test`](#test): run unit tests.
    - [`nvg validate`](#validate): check if a graph is well-formed.
    - [`nvg version`](#version): display current version number.

# bugs

Bugs can be reported at: https://github.com/vgteam/vg/issues

For technical support, please visit: https://www.biostars.org/tag/vg/

