# use this:
#  awk -f getThickness.awk  output  > data_to_be_shown
# output: the one from epics program run.
# data_to_be_shown :   x, y, thicness (in r.l) list



#  this is for usual  situation: such as InputP='u+z' 'u+y' 
$1=="h" {nl=$3; x0=$4; y0=$5; z0=$6}
$1==nl && $NF >0 { print x0,z0,$NF}

#  this is for InputP='u->sph2' (or 'g->sph2')
#  if position info on the projected plane; use the one in "projp" line
#$1=="projp" {x0=$2;y0=$3; z0=$4}
#$1=="h" {nl=$3}
#$1==nl && $NF >0 { print x0,y0,$NF}

