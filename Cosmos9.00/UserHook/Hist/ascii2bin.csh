#!/bin/tcsh
# set ARCH
if ( $# == 0 ) then
    echo "Usage: ./ascii2bin.csh  files"
    echo "  files are ascii .hist files"
    echo "  They are renamed to .ahist "
    echo "  and new binary .hist is created"
    exit 1
endif

make -f ascii2bin.mk
if ( $status != 0 ) then
    echo ascii2bin failed
    exit 1
endif

source  $COSMOSTOP/Scrpt/setarch
echo $ARCH
cat <<EOF
Each xxx.hist ascii file is renamed to xxx.ahist
and converted binary file is put as xxx.hist in the same
directory.
EOF

foreach f($*)
    echo $f
 set   ext=$f:e
 set  other=$f:r
 mv $f ${other}.ahist
 ./ascii2bin$ARCH  ${other}.ahist
end
