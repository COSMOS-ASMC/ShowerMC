#----------------------  start script
#
BEGIN{yes=0}
# lines with #include (from 1st colmumn) is printed and next line
# is processed
/^#include/ {print; next}
# f90 module may have "contains" assume not comment for old fort.
/^contains/{print; next} 
#  classic comment lines.  c or C is replaced by !
/^(c|C)/{sub("^(c|C)","!"); print;next}

#  fortran intrinsic include
/^[ \t]*include/ {print;next}

#  lines starting from  !
/^[ \t]*!/{print;next}

#    record /track/  xxxx  ==> type(track):: xxxx
#    record/track/xxxx     ==> type(track):: xxxx
#     record /epMedia/  media
/^[ \t]*record/ { print gensub(/^[ \t]*record *\/ *(\w*) *\/([^\/]*)/, "       type(\\1):: \\2",1); next}

    
#     avoid applying thie to Makefile
#          next one hit  cTrack.p.fm.p(4)    
# $1 ~ /\w*\.(f90|f|o)/ {print;next}
#  $1 ~ /\w*\.(f\,|o\,)/ {print;next}
#         space not work
# $1 ~ /\w*\.(f90|f |o\ )/ {print;next}


#           structure/element/    
#           structure /element/   both ok by the next
#                  --> type  element
/^[ \t]**structure/ { \
    $0=gensub(/^[ \t]*structure *\/ *(.*) *\//, "       type \\1 ", 1); \
    name=$2; print;  print "       sequence"; next} 
##########

/^      *end *structure/{ \
    print "       end type " name; yes=0;next}


#  track.p.fm.p(4) ==>   track%p%fm%p(4)

#    number part   12.5   0.5   .5  +.5 =.5   /.5 *.5  etc  . --> %%
#  {print gensub(/([ 0-9+-=\/\*])\.([ 0-9+-ed])/, "\\1%%\\2", "g"); next}
#   if next misunderstande number part, do above first and then next
#    and then make  @|-->.
#
#     ab.abc.xyz     P(4).fm.p(3)  ---> . -> %
{l=index($0, "!");
    if(l == 0 ) {top=$0; tail=""}
    else {top=substr($0,1,l-1); tail=substr($0,l)};
   top=gensub(/([ 0-9\+\-=\/\*])\.([ 0-9\+\-ed])/, "\\1%%\\2", "g",top);
top=gensub(/([[:alnum:]_\)]+)\.([[:alpha:]])/, "\\1%\\2", "g",top); 
top=gensub(/([[:alnum:]_\)]+)\.([[:alpha:]])/, "\\1%\\2", "g",top);
top=gensub(/([[:alpha:]\)])%((eq|le|lt|ge|gt|ne|or|not). )/, "\\1 .\\2", "g",top);
top=gensub(/( .(eq|le|lt|ge|gt|ne|or|not))%([\([:alpha:]])/, "\\1. \\3", "g",top);
top=gensub(/%%/, "\\.", "g", top); print top tail;
next
}
#{$0=gensub(/([[:alpha:]_]\w*\)*)\.([[:alpha:]])/, "\\1%\\2", "g");\
# $0=gensub(/([[:alpha:]_]\w*\)*)\.([[:alpha:]])/, "\\1%\\2", "g"); \
# $0=gensub(/([[:alpha:]])%((eq|le|lt|ge|gt|ne). )/, "\\1 .\\2", "g"); \
# print gensub(/( .(eq|le|lt|ge|gt|ne))%([[:alpha:]])/, "\\1. \\3", "g"); \
#    next}
