# extract data with a given fai index and given code from the output from
# procTime. #  (default layer is 1)
# 
# Usage: e.g
# awk -f getTimeAtFai.awk fai=4 code=2  tf.data
#
BEGIN {ly=1}
$1==fai && $2==code && $3==ly {ok=1;next}; $1==0 && ok==1{exit};ok==1{print}


# to gnuplot  data from this output ( assume temp.data)
# plot "temp.data" u (0.01*10**(($1-1)*0.1)):2
# gnuplot> replot "temp.data" u (0.01*10**(($1-1)*0.1)):4
#  etc 

