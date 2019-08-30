
#   We treat two types of bianry .hist  data and
#   ascii .hist data.
#  *) .hist data made by Assembiling procedure is incomplete
#      since it dose not include Assembled .hyb data info. yet.
#      (If you organized not to use hybrid AS data, don't worry 
#      about this. .hist is aready complete). 
#      It becomes a complete binary .hist data after being
#      modified with Assembled .hyb data.
#  *) binary .hist data can be converted to ascii .hist data
#  *) ascii  .hist data is need when converting it to 
#     files for plotting.
#   For conveniencee, we express
#      Incomplete binary .hist by .hist
#      Complete  binary  .hist by .chist
#      Complete  ascii   .hist by .achist
#      Incomplete ascii  .hist by .ahist
# 
#   Following type of jobs are supposed.
#  1)  .hist + .hyb ---> .chist  
#  2)  .hist + .hyb ---> .achist 
#  3)  .chist       ---> .achist 
#  4)  .achist      --->  plotting files
#   
#  5)  .hist        ---> .ahist
#  6)  .ahist       --->  plotting files: graph itself is corrcect but
#                                        key comment becoems wrong.
#  7)  .chist + .hyb ---> .achist (possilble but redundant)
#    
#  *)  2+4) is the shortest way, but it will be better to keep
#  complete binary .chist.  So 1+3+4) is the normal way. 
#
#  Others such as  .hist--> .hist is possible but  nonsense.
#
#  Job type
JOBTYPE=" 1  3 4"
#   next is input .hist file. needed when 
#    1), 2), 3)  [5,7]
HISTFILE0=../DisParaForTA/Assemble/$EXECID.hist
export HISTFILE0

#   next is input .hyb file.  needed when
#    1), 2)   [7]
HYBFILE0=../DisParaForTA/Assemble/$EXECID.hyb
export HYBFILE0

#  next is output .hist file. needed when
#  1)
HISTFILE1=./$EXECID.chist
export HISTFILE1

#  next is output/input ascii .achist file. needed when
#  2,3,4, [5,6] (for 5,6, .achist is not appropriate)
AHISTFILE=./$EXECID.achist
#      *******************************************
#####  ascii output is by redirection
#      *******************************************
#   next is used when
#     4) [6]. 
#   diectory below which all data and command for plotting
#   are stored.  The directory will be created.
PLOTDIR=./$EXECID-Plot
export PLOTDIR


