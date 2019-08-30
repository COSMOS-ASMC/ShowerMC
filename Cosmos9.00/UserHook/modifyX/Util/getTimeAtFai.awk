# extract data with a given fai index and given code from the output which
# in turn from  procTime.
#  (default layer is 1)
# 
# Usage: e.g
# awk -f getTimeAtFai.awk fai=4 code=2  tf.data
#
BEGIN {ly=1}
$1==fai && $2==code && $3==ly {ok=1;next}; $1==0 && ok==1{exit};ok==1{print}
