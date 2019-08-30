#!/bin/bash
if [ $# -ne 1 ]; then
    echo "Usage: ./setupenv.sh [cls|mpi]"
    echo "cls: for cluster env. job" 
    echo "mpi: for MPI env. job" 
    exit
fi

if [ $1 == "cls" ]; then
    ./setupenvcls.sh
elif [ $1 == "mpi" ]; then
    ./setupenvmpi.sh
else
    echo "usage: ./setupenv.sh [cls|mpi]"
fi
