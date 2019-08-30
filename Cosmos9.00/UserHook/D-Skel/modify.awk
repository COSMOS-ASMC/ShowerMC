toupper($0) ~ "USERHOOKC" ,  toupper($0) ~"USERHOOKI" { if( toupper($0) ~ "USERHOOKI") { print; next}; if(yet1==0) { print " UserHookc = " UserHookc ; yet1=1}; next}; 

toupper($0) ~ "SKELETONFILE" ,  toupper($0) ~ "SOURCEDEC" { if( toupper($0) ~ "SOURCEDEC") { print; next}; if(yet2==0) { print " SkeletonFile = " paramfile ; yet2=1}; next}; 

{print;next};





