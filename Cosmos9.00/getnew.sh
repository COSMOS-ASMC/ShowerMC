#/bin/bash
x=(`date`)
echo "cosnew${x[2]}"
#find . -newer Version/TimeStamp -print > /tmp/cosnew${x[2]}
find . -newer Version/TimeStamp -print > /tmp/cosnew${x[2]}
