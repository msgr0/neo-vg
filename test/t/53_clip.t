#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 22

nvg mod -U 10 msgas/hla_v.vg | nvg mod -c - > hla_v.vg

#flatten it to one path
printf "gi|568815551:1054737-1055734\t0\t1000\n" > region.bed
nvg clip hla_v.vg -b region.bed > clip_flat.vg
nvg validate clip_flat.vg
is "$?" 0 "clipped graph is valid"
for i in `nvg view clip_flat.vg | grep '^S' | awk '{print $1}'` ; do nvg view clip_flat.vg -n $i | grep ^P | grep "gi|568815551:1054737-1055734"; done | wc -l > step_count
nvg view clip_flat.vg | grep ^S | wc -l > node_count
diff step_count node_count
is "$?" 0 "every step in clipped graph belongs to reference path"
is "$(nvg paths -Ev hla_v.vg -Q "gi|568815551:1054737-1055734" | awk '{ print $2 }')" "$(nvg stats -l clip_flat.vg | awk '{ print $2 }')" "clipped graph has same length as ref path"

rm -f region.bed step_count node_count

#clip out one snarl
printf "gi|157734152:29563108-29564082\t90\t92\n" > region.bed
nvg clip hla_v.vg -b region.bed > clip.vg
nvg validate clip.vg
is "$?" 0 "clipped graph is valid"
is "$(nvg view clip.vg | grep ^S | wc -l)" "49" "Just one node filtered"

rm -f region.bed clip.vg

# clip out one edge
printf "gi|568815564:1054403-1055400\t150\t153\n" > region.bed
nvg clip hla_v.vg -b region.bed > clip.vg
nvg validate clip.vg
is "$?" 0 "clipped graph is valid"
is "$(nvg view clip.vg | grep ^L | wc -l)" "65" "Just one edge filtered"

rm -f region.bed clip.vg

# clip out low coverage node
nvg clip hla_v.vg -d 4 -P "gi|568815551:1054737-1055734" > clip.vg
nvg validate clip.vg
is "$?" 0 "clipped graph is valid"
is "$(nvg view clip.vg | grep ^S | wc -l)" "49" "Just one node filtered"

rm -f clip.vg

# clip out out-of-bounds low coverage node
printf "gi|568815551:1054737-1055734\t5\t25\n" > region.bed
nvg clip hla_v.vg -b region.bed -d 4 > clip.vg
nvg validate clip.vg
is "$?" 0 "clipped graph is valid"
nvg view hla_v.vg | sort > hla.gfa
nvg view clip.vg | sort > clip.gfa
diff hla.gfa clip.gfa
is "$?" 0 "clipping bad region changes nothing"

rm -f clip.vg hla.gfa clip.gfa

# clip out in-bounds low coverage node
printf "gi|568815551:1054737-1055734\t600\t650\n" > region.bed
nvg clip hla_v.vg -b region.bed -d 4 > clip.vg
nvg validate clip.vg
is "$?" 0 "clipped graph is valid"
is "$(nvg view clip.vg | grep ^S | wc -l)" "49" "Just one node filtered"

rm -f region.bed clip.vg

rm -f hla_v.vg

nvg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | nvg view - | sort  >  tiny.gfa
cp tiny.gfa tiny-stubs.gfa
printf "S\t0\tA\n" >> tiny-stubs.gfa
printf "L\t0\t+\t1\t+\t0M\n" >> tiny-stubs.gfa
printf "S\t100\tA\n" >> tiny-stubs.gfa
printf "L\t0\t+\t100\t+\t0M\n" >> tiny-stubs.gfa
printf "S\t200\tA\n" >> tiny-stubs.gfa
printf "L\t5\t+\t200\t+\t0M\n" >> tiny-stubs.gfa
printf "S\t300\tA\n" >> tiny-stubs.gfa
printf "L\t200\t+\t300\t+\t0M\n" >> tiny-stubs.gfa
nvg clip tiny.gfa -s -P x | sort > tiny-nostubs.gfa
diff tiny.gfa tiny-nostubs.gfa
is "$?" 0 "stub clipping removed all stubs"
diff <(nvg clip tiny/tiny.gfa -s -P x | nvg paths -v - -E | sort) <(nvg clip tiny/tiny.gfaz -s -P x | nvg paths -v - -E | sort)
is "$?" 0 "clip preserves path sequences for matching GFA and GFAZ inputs"

printf "x\t5\t25\n" > region.bed
is "$(nvg clip tiny-stubs.gfa -s -b region.bed | nvg stats -N -)" "17" "region clipping filtered out only 2 / 4 stub nodes"

printf "L\t100\t+\t2\t-\t0M\n" >> tiny-stubs.gfa
printf "L\t15\t+\t13\t-\t0M\n" >> tiny-stubs.gfa
is "$(nvg clip tiny-stubs.gfa -sS -P x | nvg stats -HT - | sort -nk 2 | awk '{print $2}' | head -1)" "1" "Correct head after path stubbification"
is "$(nvg clip tiny-stubs.gfa -sS -P x | nvg stats -HT - | sort -nk 2 | awk '{print $2}' | tail -1)" "15" "Correct tail after path stubbification"


rm -f tiny.gfa tiny-stubs.gfa region.bed tiny-nostubs.gfa

nvg clip graphs/snarl-clip.gfa -A 2 -d 2 -P x > sc-A2d2.gfa
is "$(nvg find -x sc-A2d2.gfa -n 4 | wc -l)" "0" "Node 4 correctly clipped with snarl length and depth filter"
is "$(nvg stats sc-A2d2.gfa -N)" "14" "No other nodes were clipped with snarl length and depth filter"

rm -f sc-A2d2.gfa

nvg clip graphs/chain-clip.gfa -A 2 -d1 -P x -g > sc-A2d1g.gfa
is "$(nvg stats -E sc-A2d1g.gfa)" "30" "One net edges clipped"
is "$(nvg stats -N sc-A2d1g.gfa)" "22" "No other nodes were clipped with net edge filter"

rm -f sc-A2d1g.gfa
