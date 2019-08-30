BEGIN{nrsave=0};
$0~"F1/" {nrsave=NR}
nrsave>0 && NR<nrsave+10 {print ; next}
nrsave==0 {next}
{exit}

