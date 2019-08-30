BEGIN{level[0]=""}
$3==level[$2] {next}
{level[$2]=$3; level[$2+1]="";print}


