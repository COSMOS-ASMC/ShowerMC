#   This is used to set EPICSTOP and EPICSINC at
#   installation time (make, make install)
#   Use this at Top directory as 
#   source  Scrpt/sevi

	setenv EPICSTOP `pwd`
	setenv EPICSINC $EPICSTOP/epics
        set path = ($EPICSTOP/Scrpt $path)






