BEGIN{xmin=100.; xmax=0.}
{if(xmin > $i)  {xmin=$i;fmin=$3};
 if(xmax<$i) { xmax=$i; fmax=$3}
}
END { if( i == 1)  msg="age"; else if (i == 2) msg="cog";
     {print "max " msg"=", xmax, " file=",fmax;
     print "min " msg"=", xmin, " file=",fmin}
}

        
