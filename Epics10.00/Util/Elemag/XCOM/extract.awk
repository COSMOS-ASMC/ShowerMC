$1=="#" {print ;next}
$1 ~ "#-------" {print; next}
NF==8 && ($1 !~ "MeV"  && $1 !~ "ENERGY") {$1=$1;printf("  "); print;next}
NF==10 {$1=""; $2="";print;next}
NF==9 {$1=""; print; next}
