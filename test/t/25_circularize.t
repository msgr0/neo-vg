#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 7

nvg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > linear.vg
nvg circularize -p x linear.vg > circular.vg

is $(nvg view -j circular.vg | jq -c '.path[] | select(.is_circular)' | wc -l) 1 "a path may be circularized"

nvg index -x circular.xg circular.vg
nvg convert circular.xg -v  > extracted.vg

is $(nvg view -j extracted.vg | jq -c '.path[] | select(.is_circular)' | wc -l) 1 "a circular path survives a round trip to/from xg"

nvg circularize -p y linear.vg
is $? 1 "Not allowed to circularize a nonexistent path (--path)"

echo "y" > paths.txt
nvg circularize -P paths.txt linear.vg
is $? 1 "Not allowed to circularize a nonexistent path (--pathfile)"

nvg circularize -a 2 -z 1 linear.vg
is $? 1 "Not allowed to have tail ID smaller than head ID"

nvg circularize -a 1 linear.vg
is $? 1 "Not allowed to have only a head ID"

nvg circularize -z 2 linear.vg
is $? 1 "Not allowed to have only a tail ID"

rm -f circular.vg circular.xg extracted.vg linear.vg paths.txt

