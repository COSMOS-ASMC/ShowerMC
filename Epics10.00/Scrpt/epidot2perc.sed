#  This is a sed script to transform
#     track.p.fm.p(4)
#     xyz(where).height etc
#  in a file  to 
#     track%p%fm%p(4)
#     xyz(where).height etc
#  while
#  #include "..."   (for cpp include)
#   or   
#      include "..."  (for fortran intrinsic iclude)
#  lines are untouched.
# Also comment lines are untouched even they have track.p.mass like
#  lines.  (If comment lines start from C or c in the first column
#  they will be changed to those  starting  with !)
# The converted result will be put on stdout.
#
# Usage:
#  for a single file, say, abc.f,
#   sed -f $COSMOSTOP/Scrpt/dot2perc.sed agc.f > temp.f
#   diff agc.f temp.f
#  For an actul  application, use Cosmos/Scrpt/dot2perc.sh
#  as   dot2perc.sh  inputfile
#----------------------  start script
#
# lines with #include (from 1st colmumn) is printed and next line
# is processed from the top of sed.
#    /^#include/p;d
# expected to work according to the net info.,
# but it's behaviour is very strange. So we always use next
# format.

#  classic comment lines.  c or C is replaced by !
/^[cC]/{
s/^[cC]/!/
p
d
}

 /^#include/{
 p
 d
 }
# the same for fortran intrinsic include
/ *include/{
 p
 d
}
#  lines starting from  !
/^ *!/{
p
d
}
/ *record/{
s/\( *\)record *\/\(.*\)\//\1type(\2)::/
p
d
}

#  main target ; for all lines which come here check next
{
s/\([a-zA-Z)][0-9]*\)\.\([a-zA-Z]\)/\1%\2/g
#   by the above command,
#  track.p.fm.p(4) etc will be converted to
#  track%p.fm%p(4) but p.fm part is not. this is because
#  track.p matches first, and 2nd matching search is not done for the
#  p.fm part but for .fm.p(4) so that p.fm is overlooked.
#  So next command changes such parts
s/\([a-zA-Z]\)\.\([a-zA-Z]\)/\1%\2/g
#  p1.mass type : this is included 1st s/../
# s/\([a-zA-Z][0-9]\)\.\([a-zA-Z]\)/\1%\2/g
#
#     ctrack.f  or ctrack.f90 or '...ctrack.f' might have
#   been changed to ctrack%.f etc.  restore them
s/\([a-zA-Z]\)%f90/\1.f90/g
s/\([a-zA-Z]\)%f /\1.f /g
s/\([a-zA-Z]\)%f\./\1.f\./g
s/\([a-zA-Z]\)%f'/\1.f'/g
s/\([a-zA-Z]\)%f"/\1.f"/g
p
d
  }
# comment lines starting middle
# Though looks simple, this treatment is rather difficult
# so disregards

