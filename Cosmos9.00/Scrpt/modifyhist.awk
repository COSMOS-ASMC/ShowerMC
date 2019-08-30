#   modify old ascii hist data to the one which can be processed by splitHisto.sh
# 
BEGIN{mode=0}
$1=="#hist1" {print;mode=1; xold=0;print "#c "$13; 
	      dir="#d "$13"/"$12;\
	      print "#k " $10, $11, $12; \
	      print "#t "$13,$12 ;\
	      print "#pw  2"; \
	      norm=0; \
	      next}
$1=="#hist2" {print; hist2=$0; nd=$11; mode=2; categ=$15; \
	      print "#c "categ; dir1=$15"/"$14; dir2="/d"nd; dir=dir1""dir2; \

	      key=$12" "$13" "$14" "$15; print "#k "key; \
	      if($15 == "re") {logx=1;pw=2} \
	      if($15 == "rt") {logx=1; pw=0} \
	      if($15 == "ef") {logx=0; pw=0} \
	      if($15 == "rf") {logx=0; pw=0}\
	      if($15 == "rt2"){ logx=0;pw=0}\
	      if($15 == "rz") {logx=0; pw=0}\
	      if($15 == "zf") {logx=0; pw=0}\
	      title=$14" "$15; \
	      print "#t "title; \
	      print "#pw "pw; \
	      print "#l "logx"  1";\
	      nc = 0;  \
	      next}
$1=="#hist3" {exit}
mode==2 &&  $2==0 && $3==0 { sum=0; $0=a[1]; y1=$4; $0=a[2];y2=$4;\
			     if(logx==1) {
			       bin=log(y2/y1)*0.4343; y1u=y1*10.**(bin/2);\
  			       y1l=y1/10**(bin/2); dy=y1u-y1l; dd=10.0**bin;} 
			     else dy=y2-y1;  \
			     for(i=1; i<=nc; i++) {$0=a[i]; sum+=$5*dy; \
  			       if(logx==1){y1u=y1u*dd;y1l=y1l*dd; dy=y1u-y1l}
			     }
			     print "#n  "sum" "sum ; \
			     print "#k "key" y="$3; xold=0;\
			     print "#d "dir; \
			     for(i=1; i<=nc;i++) {$0=a[i]; $5=$5/sum;\
			     print $0,$4-xold; xold=$4 } \
			     print 0,0,0,0,0,0 ; print hist2; print "#c "categ; \
			    print "#d "dir; \
			    print "#t "title ;\
			    xold=-99999;  \
			    print "#pw "pw; 
			    print "#l "logx"  1";\
			     nc=0;\
			    next}
mode==1 {if(norm == 0) {norm=$5; sum=norm;  print "#n "sum" "norm ;print dir;} \
	 print $0,  $2-xold; xold=$2; next}
mode==2 { nc++; a[nc]=$0; next}


