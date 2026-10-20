#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 2

is $(nvg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | nvg view -t -r 'http://example.org' - | wc -l) 90 "nvg view produces the expected number of lines of turtle"
is $(nvg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | nvg view -tC -r 'http://example.org' - | wc -l) 5 "nvg view produces the expected number of lines of turtle"
