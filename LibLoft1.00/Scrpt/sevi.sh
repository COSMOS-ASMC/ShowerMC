#!/bin/bash
#   This is used to set LIBLOFT  at
#   installation time (make)
#   Use this at Top directory (LibLoft/)
#   source  Scrpt/sevi.sh

	export LIBLOFT=`pwd`
	export LIBLOFTINC=${LIBLOFT}/Header
#	export COSMOSINC=${COSMOSTOP}/cosmos
	export PATH=${LIBLOFT}/Scrpt:$PATH


