#!/bin/bash
if [ $# -ne 1 ]; then
    echo "Usage: ./setupenv.sh  [cls|mpi]"
    echo "cls: for cluster environment"
    echo "mpi: for mpi environment"
    exit
elif [ $1 == "cls" ]; then
    ./setupenvcls.sh
elif [ $1 == "mpi" ]; then
    ./setupenvmpi.sh
else
    echo "Usage: ./setupenv.sh  [cls|mpi]"
    echo "cls: for cluster environment"
    echo "mpi: for mpi environment"
    exit
fi
