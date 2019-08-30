#/bin/bash
x=(`date`)
echo "epicsnew${x[2]}"
#find . -newer Version/TimeStamp -print > /tmp/epicsnew${x[2]}
find . -newer TimeStamp -print > /tmp/epicsnew${x[2]}
