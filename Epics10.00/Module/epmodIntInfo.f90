module  modIntInfo
  integer,parameter::maxcodeForInt=10
  integer,save::codeAforInt(-1:maxcodeForInt)=0
  integer,save::IntInfo1, IntInfo2 ! first and last stack position
                     !  of interaction products 
  integer,save::kintInfo ! acutal ptcle code of interacting ptcl
                      ! if >maxcodeForInt, maxcodeForInt is used
end module modIntInfo



