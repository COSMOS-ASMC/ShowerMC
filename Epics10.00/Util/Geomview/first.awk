#
#  this converts trace data into the main vertex data to be manipulated
#  by Geomview.
#
#  usage:  awk -f first.awk  code=xxx tracedata
#  where xxx is the particle code to be extracted.
# 
NF==0 && n>0  { for(i=1; i<=n; i++) print x[i], y[i], z[i];
                print ""; print ""; n=0; next};

NF>0 && ( chg == "p" ? $6 > 0: (chg == "n" ? $6 < 0: $6 == 0)) && $4 == code  \
      {n++; x[n]=$1; y[n]=$2; z[n]=$3; next}

NF>0 && ( chg == "p" ? $6 > 0 : (chg == "n" ? $6 < 0: $6 == 0) ) && code == 9 && $4 > code   \
       {n++; x[n]=$1; y[n]=$2; z[n]=$3}

END { for(i=1; i<=n; i++){ print x[i], y[i], z[i] };
                          print ""; print "";
   }

