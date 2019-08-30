#!/bin/tcsh
if( $#argv != 4 ) then
    cat <<EOF
    ./listabcd.csh  faiinex R .latfileDir outfile
    faiindex: one of 1~12; fai region
    R:    Moliere unit where 2piRrho(R) is compute
    .latfileDir: where .lat  files exist
    outfile:   output files

    outfile will containe
    title
    term captions
    code region a b c pw  age cog maxdiff R 2pirrho
    if R is not is teh "regon", 2pirrho=0.
EOF
    exit 1
endif
awk -f listabcd.awk  fai=$1 R=$2 $3/*.lat  > $4

