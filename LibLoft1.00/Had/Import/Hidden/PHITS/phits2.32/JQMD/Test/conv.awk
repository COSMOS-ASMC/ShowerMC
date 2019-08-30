
$1!="*" && $1!="c" && $1!="C" && $1=="function" {endflag="function"}
$1!="*" && $1!="c" && $1!="C" && $1=="subroutine" {endflag="subroutine"}
$1!="*" && $1!="c" && $1!="C" && $1=="block" && $2=="data" {endflag="block data"}

$1=="end" && NF==1 {print "       end " endflag;next}; {print }


