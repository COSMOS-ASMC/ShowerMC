#  This is to make transformation matrix to form an anray.
#  one unit of detector is copied to various positions.
#  input is x,y,z of each detector
#  
{print "1 0 0 0"; print "0 1 0 0"; print "0 0 1 0"; 
 print $1, $2, $3, " 1"; print " "}



