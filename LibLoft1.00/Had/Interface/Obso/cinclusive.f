!         hadron Air collision for inclusive treatment
        subroutine cinclusive(pj, a, np)
        implicit none
#include  "Zmanagerp.h"
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnp.h"
!
        type(ptcl):: pj, a(*)
        integer np
        integer maxinclusive
        integer chg, code, subcode, nchild
        parameter (maxinclusive=100)
        integer kchild(maxinclusive), chgchild(maxinclusive)
        integer subchild(maxinclusive)
        real*8  echild(maxinclusive)
        integer i
        logical first/.true./
        save    first

        if(first) then  
!               read inclusive data table.
           call rdtbl(TempDev, InclusiveFile)
           first = .false.
        endif
!
        code = pj.code
        chg = pj.charge
        subcode = pj.subcode

        call ptlint(code, chg, subcode, pj.fm.p(4),
     *   nchild, kchild, chgchild, subchild, echild)

        if(nchild .gt. maxinclusive) then
           call cerrorMsg(
     *      '# of ptcls by inclusive prod. exceeded limit', 1)
           call cerrorMsg(
     *       'enlarge maxinclusive in cinclusive.f', 0)
        endif
!
        do i = 1, nchild
           call cmkptc(kchild(i), subchild(i), chgchild(i),
     *         a(i))
           a(i).fm.p(4) = echild(i) + a(i).mass
           a(i).fm.p(1) = 0.  
           a(i).fm.p(2) = 0.
           a(i).fm.p(3) = sqrt( echild(i) * (echild(i)+2*a(i).mass))
        enddo
        np = nchild
        end


