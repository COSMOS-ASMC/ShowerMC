#   This is used to set EPICSTOP and EPICSINC at
#   installation time (make, make install)
#   Use this at Top directory as 
#   source  Scrpt/sevi

	export  EPICSTOP=`pwd`
	export  EPICSINC=$EPICSTOP/epics
        export  PATH=$EPICSTOP/Scrpt:$PATH







