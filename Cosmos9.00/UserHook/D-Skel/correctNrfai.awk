BEGIN{modify=0}
$1 != "dE/dx" &&  modify == 0 {print;next}
$1 == "dE/dx" {print; modify=1;next}
modify==1 && NF>0 {for(i=1; i<=NF; i++) printf("%11.3E", $i=$i/enhance); printf("\n");next}
  {print}

