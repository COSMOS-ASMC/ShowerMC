NR==1{print; pcode=$1; psubcode=$2; pchg=$3; ppx=$4; ppy=$5; ppz=$6;next}
NR==2{print; tcode=$1; tsubcode=$2; tchg=$3; tpx=$4; tpy=$5; tpz=$6;next}
END {print "end"}

