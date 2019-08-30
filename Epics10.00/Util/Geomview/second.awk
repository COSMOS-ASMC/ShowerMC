# this is used to get header from output by first.awk
# Usage:   awk -f second.awk  outputfromfirst
#
BEGIN {prev="0"; tv=0; nv=0}
NF > 1 && prev=="0" {nl++; prev="1"};
NF > 1 {nv++; tv++}
NF == 0 && prev=="1" {v[nl]=nv; 
       if( nl == 1 ) c[nl]=1;
       else  c[nl] = 0;
       prev="0"; nv=0}
#        last 1 is 1 color for this particle code.
END { if( prev=="1" )  {v[nl]=nv; 
			 if( nl == 1 ) c[nl]=1;
			 else  c[nl] = 0;
		       }
#        last 1 is 1 color for this particle code.
	print nl, tv, 1;
        for(i=1; i<=nl; i++) print  v[i];
        print ""; 
	for(i=1; i<=nl; i++) print  c[i];
        print ""
    }

