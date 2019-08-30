#!/bin/bash
#   This is used to set COSMOSTOP and COSMOSINC at
#   installation time (make, make install)
#   Use this at Top directory as 
#   source  Scrpt/SetEnv

	export COSMOSTOP=`pwd`
#	export COSMOSINC=${COSMOSTOP}/cosmos
	export PATH=${COSMOSTOP}/Scrpt:$PATH


