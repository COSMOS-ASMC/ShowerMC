awk '$1=="ARCH" && $2=="=" {print $3}' $COSMOSTOP/site.config
