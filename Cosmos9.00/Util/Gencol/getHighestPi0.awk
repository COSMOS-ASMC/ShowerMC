BEGIN {Emax=0}
$1=="h" {if(Emax > 0) print Emax, code, chg; Emax=0; next}
$1==4 {if($7>Emax ) { Emax=$7; code=$1; chg=$3 }}
END {if(Emax > 0) print Emax, code, chg; Emax=0}

    
