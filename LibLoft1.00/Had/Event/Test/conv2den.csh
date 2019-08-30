#/bin/tcsh
@ i = 1
while ( $i < 12 )
 set  norm = `awk 'END{print 3.1415*(r**2 - (r-1)**2)*1.e5}' r=$i /dev/null`
 set  nm = `awk 'END {print (r*2-1)/2}' r=$i /dev/null`
 echo $norm
 echo $nm
 echo "histn_"$nm
# awk '$2 > r-1 && $2 < r {print $1}' r=$i $1 | histo -l 80 0.025  $norm > histn_$nm
 awk '$2 > r-1 && $2 < r {print $1}' r=$i $1 | histo 1000  50  $norm > histn_$nm
 @ i++
 echo $i
end
