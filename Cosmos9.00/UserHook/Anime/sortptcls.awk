$4==2 {print $1, $2, $3 >> "electron.dat"; next};
$4==3 {print $1, $2, $3 >> "muon.dat" ;next};
$4>3 {print $1, $2, $3 >> "other.dat" ;next};

