BEGIN {put="no"}
$1=="#hist1"  { if($2==evn) {put="yes"; print; next}
		else {
		  if ( $2< evn) { put="no"; next}
		  else exit;
		} 
}
$1=="#hist2"  { if($2==evn) {put="yes"; print; next}
		else {
		  if ( $2< evn) { put="no"; next}
		  else exit;
		} 
}
$1=="#hist3"  { if($2==evn) {put="yes"; print; next}
		else {
		  if ( $2< evn) { put="no"; next}
		  else exit;
		} 
}

put=="yes" {print;next}


