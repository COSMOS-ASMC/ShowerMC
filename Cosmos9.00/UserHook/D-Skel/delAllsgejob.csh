qdel `qstat | awk  '$4==user {print $1}' user=$USER`
