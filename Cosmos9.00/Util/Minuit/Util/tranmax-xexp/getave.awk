{gsum1+=$1; gsum2+=$2; esum1+=$3; esum2+=$4; Esum1+=$5; Esum2+=$6; hsum1+=$7; hsum2+=$8; n++}
END {print "Ng(", gsum1/n, gsum2/n,") Ne(", esum1/n, esum2/n,") dE(", Esum1/n, Esum2/n,") hNe (",hsum1/n, hsum2/n,")"}


