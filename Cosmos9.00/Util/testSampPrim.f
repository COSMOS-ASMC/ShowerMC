!              test primary sampling routines
#include  "BlockData/cblkGene.h"

      implicit none
#include  "Zmanagerp.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zprimary.h"
#include  "Zcode.h"
      
      type (ptcl):: p
      type (track)::incident
      type (coord):: angle

      real*8 eorporrig
      integer label, icon, i
      real*8 mom, kepn, kepA
      integer fin

      call creadParam(5)
      call cbeginRun
      call cprintPrim(ErrorOut)
      write(0,*) 'id code subc chg  orgE TKE TE P KE/n'
      do i = 1, abs(DestEventNo(1))
         icon = 1
         do while (icon .ne. 0)
            call csampPrimary(incident%p, fin)
            if(fin .ne. 0) goto 10 
            call csPrimAng(angle)
            call cmkInc(incident, angle)
            if(CutOffFile .ne. ' ') then
               call cifCutOff(icon)
            else
              icon =0
            endif
         enddo
 10      continue
         call cqPrimE(eorporrig)
         call cqPrimLabel(label)
         p = incident%p
         kepn =( p%fm%p(4)-p%mass) 

         kepA = kepn
         if(p%code .eq. kgnuc) then
            kepn = kepn / p%subcode
         endif

         mom =  sqrt(p%fm%p(4)**2 - p%mass**2)
         write(*, '(i2, 3i4,1p, 5g14.4)')
     *      label, 
     *      p%code, p%subcode, p%charge,
     *      eorporrig, kepA, p%fm%p(4),  mom, kepn
      enddo
      end
      subroutine chookTrace
      end
      subroutine chookCeren
      end
      subroutine chookCerenS
      end
      subroutine chookCerenE
      end
      subroutine chookBgRun
      end
