BEGIN{Ne1=3.e10; Ne2=7.e10}
$1=="h" {Ne=0;next}
$1=="t" {if($8>Ne) Ne=$8}
END {print "Nemax=   " Ne;
 if(Ne > Ne2 || Ne <Ne1) {
     print "Check transition !!!***********************************";
 }
}

