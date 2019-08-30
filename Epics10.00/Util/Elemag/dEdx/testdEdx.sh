#!/bin/bash
if [ $# != 3  ] && [ $# != 4 ] ; then
    cat <<EOF
Usage: ./testdEdx.sh  code subcode charge [outputfile]
     where 
           code:  incident particle code (e-->2, mu->3, pi->4, K->5,p->6
                  heavy 9)  
        subcode:  anti particle-->1  regular particle -1 
                  heavy--> A
         charge:  Z
      ouputfile:  file path where  the result is stored; if not given, 
                  stdout 
      *** other implicit input:
      ./epicsfile:  StoppingPw is used
                    SrimEmax is used
                 non-zero RecoilKeMin is used if AutoEmin=0
        ./config:  Media is obtained from here
EOF
exit
fi


if [ $3 -eq 0 ]; then
    echo charge must not be 0
    exit
fi

make  -f dEdx.mk clean 
make -f dEdx.mk


if [ $# -eq 4 ] ; then
  echo  ./epicsfile ./config $1 $2   $3   | ./a.out > $4
  echo output has been directed to  $4
else
  echo  ./epicsfile ./config $1 $2   $3  | ./a.out 
fi
