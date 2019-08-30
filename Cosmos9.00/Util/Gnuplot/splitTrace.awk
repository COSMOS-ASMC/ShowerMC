BEGIN {ok=0; nc=0; Emin=10e-3;
      fileN[2,-1]="trace_e-"; fileN[2,1]="trace_e+";
      fileN[3,-1]="trace_mu"; fileN[3,1]="trace_mu";
      fileN[4,-1]="trace_pi"; fileN[4,1]="trace_pi";
      fileN[5,-1]="trace_K";  fileN[5,1]="trace_K";
      fileN[6,-1]="trace_p";  fileN[6,1]="trace_p";
      fileN[9, 1]="trace_H";
      }
NF > 0 && $4 < 7 && $6 != 0 && ( $5>Emin  || nc==1) {
    nc++; ok=1;  out=Work"/"fileN[$4,$6];
    print $1, $2, $3 > out;     next;
};
#(NF > 0 && $4==9 || $6 !=0 )  {ok=1;
 NF > 0 && $4==9 && ( $5>Emin  || nc==1) {
     nc++; ok=1;
    out=Work"/"fileN[9,1];
    print $1, $2, $3  > out; next;
}
 NF==0 && ok==1  {print ""> out; print "" >out; ok=2;nc=0; next;
}
