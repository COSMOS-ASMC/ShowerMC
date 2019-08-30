#   split an ascii hist file (may be  made by bin2ascii)
#   into individual
#   histogram files. with a directory hierarchy 
#   change style as you like
BEGIN{fc=0; mode=0; dirsave=""; nc=0;  dq = "\""; bs = ",\\"; hlab ="x";
    style="w his lw 3";

};


$1=="#hist1" {mode=1; key=""; normf=1; close(file); file="";  next};
$1=="#hist2" {mode=2; key=""; normf=1; close(file); file=""; next};
$1=="#hist3" {mode=3; key=""; normf=1; close(file); file=""; next};

$1=="0" && $2=="0" && $3=="0" {next};  # end of 1 hist

$1=="#t" {$1=""; title=$0; next};
$1=="#x" {xlab=$2; xunit=$3; if(mode==1 && xlab !="") { hlab=xlab; hunit=xunit}; next};
$1=="#y" {ylab=$2; yunit=$3; if(mode==2 && ylab !="") { hlab=ylab; hunit=yunit}; next};
$1=="#z" {zlab=$2; zunit=$3; if(mode==3 && zlab !="") { hlab=zlab; hunit=zunit}; next};
$1=="#pw" {pw=$2;next};
$1=="#dN" {dNu=$2;next};
$1=="#k" {$1=""; key=$0; next};
$1=="#l" {logx=$2; logy=$3;next};
$1=="#c" {categ=$2; next}
#     normf == 1 means dN/dx is not normalized. nf=1./sum;  $3=dN/dx*nf
#           != 1 means dN/dx is  normalized.  nf=sum; $3=dN/dx*nf
$1=="#n" {sum=$2;normf=$3;
	  if( normf == 1.0) nf=1./sum;
	  else  nf = sum; 
	  next;}

$1=="#o" {imin=$2; imax=$3; inc=$4;next}

$1=="#d" {dirspec(); next}
# only next may come after "#d"
$1=="#f" { filespec(); next}
file=="" { filespec2();next}

fc==0 { fc++; mktab();	next}

      { mktab(); next}

END { rmlastbs(); }
#####################
function mkdir(){
   if(dir == "") {
     if( categ == "") dir=maindir"/";
     else dir=maindir"/"categ"/";
   }
   else {
     dir=maindir"/"dir"/";
   }
   system("mkdir -p " dir);
   return dir;
 }
function putcomment(){
  if(mode==1) print "# tab is: x dn/dx dn/dx* dn  n(>x) <x> dx" >> gfile;
  if(mode==2) print "# tab is: y dn/dx/dy dn/dx/dy* dn  dy x" >> gfile;
  if(mode==2) print "# tab is: z dn/dx/dy/z dn/dx/dy/dz* dn dz x y" >> gfile;
  if( normf == 1.0 ) {
    print "# graph(dn/dx..) is not normalzied." >> gfile;
    print "# To show normalzied one(dn/dx..*); change $2-->$3 in " \
      " the last call command" >> gfile;
  }
  else {
    print "# graph(dn/dx..) is normalzied." >> gfile;
    print "# To show unnormalzied one(dn/dx..*); change $2-->$3 in" \
      " the last call command" >> gfile;
  };
#
  print "# To change the line style etc, modify the last "dq"w his"dq >> gfile;
}

#####################
function dirspec() {
  dir=$2; dir=mkdir();
  if( dir == dirsave ) {nc++;
##    print "   "  dq""nc".dat"dq	\
##	" using 1:($1**pw*$2) title " dq""key""dq " " style  bs  >> vfile;
    }
  else {
    close(gfile); close(vfile);
#    rm last ,\ in vfile
    if(dirsave != "" ) rmlastbs();

    nc=1; dirsave=dir ;
    
    gfile=dir"plot.gp";
    vfile=dir"var.gp";
    print "# variable part of gnuplot commands" > vfile;
##   print "plot " dq""nc".dat"dq	\
##      " using 1:($1**pw*$2) title " dq""key""dq " " style   bs  >> vfile;

    print "# gnuplot commands" > gfile;

##    print "style="dq""style""dq >> gfile;
    print "set title "dq""title""dq \
      " font "dq"Times,20"dq >> gfile;

    mkxlab();
    print "set xlabel " dq""xlabel""dq  \
      " offset  0,0.5 font "dq"Times-Italic,22"dq >> gfile;

    print "pw ="dq""pw""dq >> gfile;

    mkylab();
    print "set ylabel " dq""ylabel""dq" offset  2.5,0  font" \
      " "dq"Times-Italic,22"dq >> gfile;

    print "set grid mxtic xtic  mytic ytic" >> gfile;
    if( logx == 1 ) print "set log x" >> gfile;
    if( logy == 1 ) print "set log y" >> gfile;
    print "#  set key x,y" >> gfile;
    print "#  set xrange [1:10]" >> gfile;
    print "#  set yrange [1:10]" >> gfile;
#   print "#  sum="sum " norm="normf >> vfile;
    putcomment();
    print "call "dq"var.gp"dq   >> gfile;
  };
};
#################
function filespec() {
  if(dir =="") dir=mkdir(); 
  file=dir$2; fc=0;
  if(nc ==1) print "plot " dq""$2""dq \
	       " using 1:($1**pw*$2) title "dq""key""dq " " style bs  >> vfile;
  else {if (nc>1) print "  "  dq""$2""dq \
		    " using 1:($1**pw*$2) title "dq""key""dq " " style  bs  >> vfile};
}
##############
function  filespec2(){
  if(dir =="") dir=mkdir();
  file=dir""nc".dat"; fc=1;
  if( nc == 1) {
     print "plot " dq""nc".dat"dq\
      " using 1:($1**pw*$2) title " dq""key""dq " " style   bs  >> vfile }
  else
  { print "   "  dq""nc".dat"dq						\
   " using 1:($1**pw*$2) title " dq""key""dq " " style  bs  >> vfile
  }
#   mktab();

  
#  if( mode == 1) {print $2,$3,$3*nf, $4,$5,$6,$7 > file };
#  if( mode == 2) {print $4,$5,$5*nf, $6,$7,$3 > file};
#  if( mode == 3) {print $6,$7,$7*nf, $8,$9,$4,$5 > file};
}

##############remove last ,\
function rmlastbs() {
      system("mv "vfile"   temp.gp");
      system("sed '$,$s/,\\\\//' temp.gp  > "vfile);
      close("sed '$,$s/,\\\\//' temp.gp  > "vfile);
###??      close("mv "vfile"   temp.gp");
      system("rm -f temp.gp");
}
#################
function mktab(){
  if( mode == 1) {print  $2,$3,$3*nf, $4,$5,$6,$7 > file };
  if( mode == 2) {print  $4,$5,$5*nf, $6,$7,$3 > file};
  if( mode == 3) {print  $6,$7,$7*nf, $8,$9,$4,$5 > file};
}
#############
function mkxlab() {
  if(hunit != "") {
    xlabel=hlab"("hunit")";}
  else {
    xlabel=hlab;
  } 
}
##############
function mkylab(){
#               (x=hlab)
#           normalized: pw==1-->ylabel= xdN/Ndx           
#                       pw==0-->ylabel= dN/Ndx (/hunit) 
#                       other  pw  yhlabel= x^pwdN/Ndx (hunit^pw-1)
#           unnormalized  pw=1--> ylabel =xdN/dx (dNu)
#                         pw=0--> ylabel = dN/dx (dN/hunit)
#                       otgher pw ylabel = x^pwdN/dx (dNu*hunit^pw-1)
#  
  if(normf != 1.0) {
#      normalized
    if(pw == 0.) { 
      if(hunit != "") {
	ylabel="dN/Nd"hlab"(/"hunit")";}
      else {
	ylabel="dN/Nd"hlab; }
      return;
    }
    if(pw == 1.) {ylabel=hlab"dN/Nd"hlab; return};
###    other pw 
    pwm=pw-1; 
    if(hunit != "") {
      ylabel=hlab"^"pw"dN/Nd"hlab"("hunit"^"pwm")";}
    else {
      ylabel=hlab"^"pw"dN/Nd"hlab; }
    return;
  }
  else {
#    unnormalized
    if(pw == 0.)  {
      if(hunit != "" ) {
	ylabel="dN/d"hlab"("dNu"/"hunit")"}
      else {
	ylabel="dN/d"hlab"("dNu")";}
      return;
    }
    if( pw == 1.) {
      if(dNu != "") {
	ylabel=hlab"dN/d"hlab" ("dNu")"}
      else {
	ylabel=hlab"dN/d"hlab}
      return;
    }
## other pw
    pwm=pw-1; 
    if( hunit != "" ) {
      ylabel=hlab"^"pw"dN/d"hlab"("dNu"/"hunit"^"pwm")"}
    else { 
      if(dNu != "" ) {
	ylabel=hlab"^"pw"dN/d"hlab"("dNu")";}
      else ylabel=hlab"^"pw"dN/d"hlab;
    }
  }
}
