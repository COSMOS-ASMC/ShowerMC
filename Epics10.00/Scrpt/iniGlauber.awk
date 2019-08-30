BEGIN {nel=0; n=0; end1 = 0; loc=0}
!/^#/ && NF==7 && nel == 0 {nel = $1;  next}
!/^#/ && nel>0 && end1 == 0 {++n; z[n]=$1; a[n]=$2; w[n]=$3}
n> 0 && nel == n && end1==0 {
 print "* dpmjet input file for preinitialization of Glauber calculaiton";
 print "* This is for " BaseM;
 end1=1;
 next}
end1 == 1 && $2 == "EMULSION" {
   for(j=1; j<=n; j++) {
     printf("%s%10d%10d%10.3f\n","EMULSION  ",int(a[j]+0.5),z[j],w[j]);
   }
 next}
end1 == 1 && $6 == "dpmjet" {$6=filename;
  printf("%-10s%10.1f%10.1f%10.1f%10.1f%-20s%-10s\n", $1, $2, $3, $4, $5, \
  "                  ", $6);next
  }
#end1 ==1 {loc=index(FILENAME,"GlaubIniTemplate.inp")};
end1 ==1 {loc=index(FILENAME,"GlaubIniTemplate")}
end1 ==1 && NF>0 && loc>0 {print; next }
end1 ==1 && NF>0 && loc==0 {print "*", $0 }
