#  awk -f ./convSrim.awk Zin=xx Z4-SCIN.srim > outputfile
$6=="Mass" {A=$8; i=index($4,"]"); Z=substr($4,2,i-2);  \
    Zin=Zin*1; if(Z != Zin) {print "error: Zin="Zin"!=internal Z="Z > "/dev/stderr"; exit};   next}
$1=="-----------" {ok=1;next};ok==1 && NF==1 {exit}; 
ok==1 {if($2=="eV") $1=$1/1e9;
    else if($2=="keV") $1=$1/1e6;
    else if($2=="MeV")  $1=$1/1e3;
    else if($2=="GeV")  $1=$1;
    else {print "err"  > "/dev/stderr"; exit}
    print $1,$3+$4
}
