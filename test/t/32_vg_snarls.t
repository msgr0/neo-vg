#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 42

nvg view -J -v snarls/snarls.json > snarls.vg
nvg snarls -t 1 snarls.vg -r st.pb > snarls.pb
is $(nvg view -R snarls.pb | wc -l) 3 "nvg snarls made right number of protobuf Snarls"
is $(nvg view -E st.pb | wc -l) 6 "nvg snarls made right number of protobuf SnarlTraversals"
is $(nvg view -R snarls.pb | jq -r '[(.start.node_id | tonumber), (.end.node_id | tonumber)] | min' | tr '\n' ',') "1,3,7," "nvg snarls made snarls in the right order"

nvg snarls -l -t 1 snarls.vg > /dev/null 2> error.txt
is $? 1 "nvg snarls does not allow --leaf-only to be used without --traversals"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"
nvg snarls -l -t 1 snarls.vg -r /dev/null > leaf_only.pb
diff leaf_only.pb snarls.pb
is $? 0 "nvg snarls --leaf-only matches full snarls"

nvg snarls -o -t 1 snarls.vg > /dev/null 2> error.txt
is $? 1 "nvg snarls does not allow --top-level to be used without --traversals"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"
nvg snarls -o -t 1 snarls.vg -r /dev/null > top_level.pb
diff top_level.pb snarls.pb
is $? 0 "nvg snarls --top-level matches full snarls"

nvg snarls -a -t 1 snarls.vg > /dev/null 2> error.txt
is $? 1 "nvg snarls does not allow --any-snarl-type to be used without --traversals"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"
nvg snarls -a -t 1 snarls.vg -r /dev/null > any.pb
diff any.pb snarls.pb
is $? 0 "nvg snarls --any-snarl-type matches full snarls"

nvg snarls -m 5 -t 1 snarls.vg > /dev/null 2> error.txt
is $? 1 "nvg snarls does not allow --max-nodes to be used without --traversals"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"
nvg snarls -m 5 -t 1 snarls.vg -r /dev/null > small.pb
diff small.pb snarls.pb
is $? 0 "nvg snarls with --max-nodes changed matches full snarls"

rm -f snarls.pb st.pb leaf_only.pb top_level.pb any.pb small.pb error.txt

nvg index snarls.vg -x snarls.xg
is $(nvg snarls snarls.xg -r st.pb | nvg view -R - | wc -l) 3 "nvg snarls on xg made right number of protobuf Snarls"
is $(nvg view -E st.pb | wc -l) 6 "nvg snarls on xg made right number of protobuf SnarlTraversals"

rm -f snarls.vg snarls.xg st.pb

# vcf alt traversals in tiny graph
nvg construct -Saf -v tiny/tiny.vcf.gz -r tiny/tiny.fa > tiny.vg
nvg snarls tiny.vg -r tiny.exhaustive.trav > /dev/null
nvg snarls tiny.vg -r tiny.vcf.trav -v tiny/tiny.vcf.gz > /dev/null
nvg view -E tiny.exhaustive.trav | sort > tiny.exhaustive.trav.sort
nvg view -E tiny.vcf.trav | sort > tiny.vcf.trav.sort
diff tiny.exhaustive.trav.sort tiny.vcf.trav.sort
is $? 0 "vcf traversals are the same as exhaustive traversals for tiny graph"

nvg snarls tiny.vg -v tiny/tiny.vcf.gz > /dev/null 2> error.txt
is $? 1 "nvg snarls requires --traversals to be given if --vcf is"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"

nvg snarls tiny.vg -r tiny.vcf.trav -v tiny/tiny.vcf.gz -e > /dev/null 2> error.txt
is $? 1 "nvg snarls does not allow both --vcf and --path-traversals"
grep "cannot be used with" error.txt
is $? 0 "the problem is explained"

nvg snarls tiny.vg -e > /dev/null 2> error.txt
is $? 1 "nvg snarls requires --traversals to be given if --path-traversals is"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"

nvg snarls tiny.vg -e -m 5 -r /dev/null > /dev/null 2> warning.txt
grep "ignored" warning.txt
is $? 0 "warning about ignoring is given"

rm -f tiny.vg tiny.exhaustive.trav.sort tiny.exhaustive.trav tiny.vcf.trav.sort tiny.vcf.trav error.txt warning.txt

# vcf alt traversals on an inversion
nvg construct -Saf -v sv/x.inv.vcf -r sv/x.fa > x.inv.vg
nvg snarls x.inv.vg -a -r x.inv.exhaustive.trav > /dev/null
nvg snarls x.inv.vg -a -r x.inv.vcf.trav -v sv/x.inv.vcf -f sv/x.fa > /dev/null
nvg view -E x.inv.exhaustive.trav | sort > x.inv.exhaustive.trav.sort
nvg view -E x.inv.vcf.trav | sort > x.inv.vcf.trav.sort
diff x.inv.exhaustive.trav.sort x.inv.vcf.trav.sort
is $? 0 "vcf traversals are the same as exhaustive traversals for inversion"

nvg snarls x.inv.vg -a -r x.inv.vcf.trav -f sv/x.fa > /dev/null 2> error.txt
is $? 1 "nvg snarls requires --vcf to be given if --fasta is"
grep "requires --vcf" error.txt
is $? 0 "the problem is explained"

nvg snarls x.inv.vg -a -r x.inv.vcf.trav -i sv/x.fa > /dev/null 2> error.txt
is $? 1 "nvg snarls requires --vcf to be given if --ins-fasta is"
grep "requires --vcf" error.txt
is $? 0 "the problem is explained"

rm -f x.inv.vg x.inv.exhaustive.trav.sort x.inv.exhaustive.trav x.inv.vcf.trav.sort x.inv.vcf.trav error.txt

# vcf alt traversals on ins_and_del
nvg construct -Saf -v ins_and_del/ins_and_del.vcf.gz -r ins_and_del/ins_and_del.fa > ins_and_del.vg
nvg snarls ins_and_del.vg -r ins_and_del.exhaustive.trav > /dev/null
nvg snarls ins_and_del.vg -r ins_and_del.vcf.trav -v ins_and_del/ins_and_del.vcf.gz > /dev/null
nvg view -E ins_and_del.exhaustive.trav | sort > ins_and_del.exhaustive.trav.sort
nvg view -E ins_and_del.vcf.trav | sort > ins_and_del.vcf.trav.sort
diff ins_and_del.exhaustive.trav.sort ins_and_del.vcf.trav.sort
is $? 0 "vcf traversals are the same as exhaustive traversals for ins_and_del graph"

rm -f ins_and_del.vg ins_and_del.exhaustive.trav.sort ins_and_del.exhaustive.trav ins_and_del.vcf.trav.sort ins_and_del.vcf.trav

# parallelizing on components (which is deactivated with -t 1)
nvg construct -r small/xy.fa -v small/xy.vcf.gz > xy.vg
nvg construct -r small/xy.fa -v small/xy.vcf.gz -R x > x.vg
nvg construct -r small/xy.fa -v small/xy.vcf.gz -R y > y.vg
nvg snarls x.vg > xy.snarls
nvg snarls y.vg >> xy.snarls
is $(nvg snarls xy.vg | nvg view -R - | wc -l) 35 "correct number of snarls when parallelizing on compoents"
is $(nvg snarls xy.vg | nvg view -R - | wc -l) $(nvg view -R xy.snarls | wc -l) "same number of snarls when parallelizing on components"
rm -f xy.vg xy.snarls x.vg y.vg

# Find trivial snarls from a GBZ
nvg gbwt -g graph.gbz --gbz-format -G graphs/components_walks.gfa
nvg snarls -T graph.gbz > graph.snarls
is $? 0 "GBZ graphs can be used for finding snarls"
is $(nvg view -R graph.snarls | wc -l) 5 "correct number of snarls in the GFA W-line example"
rm -f graph.gbz graph.snarls

# Check snarl output order in a more complex case.
# Snarls need to come out in order along chains, recursing down and then coming back up to the right place in the chain.
nvg view -J -v snarls/nested.json > nested.vg
nvg snarls -t 1 -T nested.vg > nested.snarls
# Could go either way first in the first snarl
echo "1,2,5,5,7,3,9," > possibilities.txt
echo "1,3,2,5,5,7,9," >> possibilities.txt
GOT="$(nvg view -R nested.snarls | jq -r '[(.start.node_id | tonumber), (.end.node_id | tonumber)] | min' | tr '\n' ',')"
is $(echo $GOT | wc -c) 15 "nvg snarls made the right snarls when dealing with nested chains"
grep $GOT possibilities.txt >/dev/null
is $? 0 "nvg snarls made snarls in the right order when dealing with nested chains"
rm -f nested.vg nested.snarls possibilities.txt

# Find snarls from a GFA in GFA space
is "$(nvg snarls --named-coordinates graphs/components_walks_named.gfa | nvg view -Rj - | jq -c '[.start.name, .start.backward, .end.name, .end.backward]' | sort)" \
   "$(printf '["cats",null,"goldfish",null]\n["goldfish",null,"sparrows",null]\n["pigs",true,"squirrels",null]\n["squirrels",null,"rabbits",true]\n')" \
   "Correct snarls are found in GFA space"

# Find traversals from a chopped GBZ-ified GFA in node ID space
nvg autoindex -p index -g graphs/big_snarl_named.gfa -w giraffe
nvg snarls index.giraffe.gbz --traversals traversals.dat >/dev/null
is "$(nvg view -Ej traversals.dat | jq -c 'select(.visit | length > 2) | .visit[] | select(.node_id // .name)' | wc -l | tr -d " ")" "7" "In node ID space traversal visits 7 nodes"
is "$(nvg view -Ej traversals.dat | jq -c 'select(.visit | length > 2) | .visit[] | select(.snarl)' | wc -l | tr -d " ")" "4" "In node ID space traversal visits 4 trivial snarls"

# And in GFA space
nvg snarls index.giraffe.gbz --traversals traversals.dat --named-coordinates >/dev/null
is "$(nvg view -Ej traversals.dat | jq -c 'select(.visit | length > 2) | .visit[] | select(.node_id // .name)' | wc -l | tr -d " ")" "3" "In segment name space traversal visits 3 segments"
is "$(nvg view -Ej traversals.dat | jq -c 'select(.visit | length > 2) | .visit[] | select(.snarl)' | wc -l | tr -d " ")" "0" "In segment name space traversal visits 0 trivial snarls"

rm -f index.* traversals.dat
