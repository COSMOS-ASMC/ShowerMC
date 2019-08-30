      function varimlt(iin,iout,ie)
      implicit real*8 (a-h,o-z)
      common /variparm3/varim(16,10,9)
      varimlt = varim(ie,iout,iin)
      end

      function variamp(iin,iout,ie)
      implicit real*8 (a-h,o-z)
      common /variparm2/varia(16,10,9)
      variamp = varia(ie,iout,iin)
      end

      function varifunc(x0, x)
      implicit real*8 (a-h,o-z)
      parameter (p0=1.d0)
      x1 = (x-x0)
      varifunc =  x1/sqrt(x1**2 + p0**2)
      end

      subroutine initvarip(filn, outputfile)
      implicit real*8 (a-h,o-z)
      common /variparm2/varia(16,10,9)
      common /variparm3/varim(16,10,9)

      dimension dae(16,10,9)
      logical rescale, logfilexist
      character line*80, filn*80, logfil*80, cptl(9)*3, outputfile*80
      data cptl/'pi+','pi-','k+ ','k- ','k0 ','p+ ','p- ','n+ ','n- '/

      do 101 i=1,9
         do 102 j=1,10
            do 103 k=1,16
               dae  (k,j,i) = 0.
               varim(k,j,i) = 0.
               varia(k,j,i) = 0.
 103        continue
 102     continue
 101  continue

      open (8,file=filn,status='old',err=2999)
      do 1000 iii = 1, 1000
         read(8,'(a)',end=1999) line
         do 1001 ip = 1, 9
            if(line(1:3).eq.cptl(ip)) then
               iin = ip
               go to 1002
            end if
 1001    continue
         go to 1000
 1002    continue

         if(line(4:4).eq.'m')then
            do 1003 ie=1, 16
               read(8,*) ie1, (varim(ie,iout,iin), iout=1,10)
               if(ie1.ne.ie) then
                  stop ': Interaction variation file error'
               end if
 1003       continue
         else
            do 1004 ie=1, 16
               read(8,*) ie1, (dae(ie,iout,iin), iout=1,10)
               if(ie1.ne.ie) then
                  stop ': Interaction variation file error'
               end if
 1004       continue
         end if
            
 1000 continue
 1999 continue
      close(8)

      ll=lchend(filn, 80)
      write(0,*) 'Modification of interaction from ', filn(1:ll)
      do 1100 iin = 1, 9
         do 1101 ie = 1, 16
c                                        Check the energy balance 
            detot  = 0.
            acomp  = 0.
            do 1102 iout = 1, 10
               amult = 10.d0**hlmulti(iin, iout, ie)
               if(iout.eq.7 .or. iout.eq.9) then
c                                  Compensation factor by nucleon and proton
                  call getvarip(iin,iout,ie,xl0, xint, vxint)
c                  acomp = acomp + amult * vxint
                  acomp = acomp + amult * xint
               else
                  call getvarip(iin,iout,ie,xl0,xint,vxint)
                  if(abs(vxint).lt.1.d-3) then
                     dae(ie,iout,iin) = 0.
                     varia(ie,iout,iin) = 0.
                  else
                     dae(ie,iout,iin) = dae(ie,iout,iin)/100
                     varia(ie,iout,iin) = dae(ie,iout,iin)*xint/vxint
                  end if
c                                        |Variation parameter| < 0.9
                  if(varia(ie,iout,iin).gt.0.9d0) then
                     dae(ie,iout,iin) = 0.9d0/abs(varia(ie,iout,iin))
     &                    * dae(ie,iout,iin)
                     varia(ie,iout,iin) = sign(0.9d0,varia(ie,iout,iin))
                  end if
                  varim(ie,iout,iin)=varim(ie,iout,iin)/100
                  detot = detot + amult * xint * 
     &             ( (1.d0+dae(ie,iout,iin))*(1.d0+varim(ie,iout,iin))
     &                - 1.0d0 )
               end if
 1102       continue

            rcomp = -detot/acomp
            rescale = (abs(detot).gt.acomp)
            if(rescale) then
               rfact = abs(0.9d0/rcomp)
               rcomp = -sign(0.9d0,rcomp)
            end if

            do 1103 iout=1, 10
               if(iout.eq.7.or.iout.eq.9) then
                  call getvarip(iin,iout,ie,xl0, xint, vxint)
                  varia(ie,iout,iin) = 0.0
c                  varia(ie,iout,iin) = rcomp
                  
                  varim(ie,iout,iin) = rcomp
                  dae  (ie,iout,iin) = 0.0
               else if(rescale) then
                  a = dae  (ie,iout,iin)
                  b = varim(ie,iout,iin)
                  if(abs(a*b).gt.1e-6) then
                     rf = (sqrt((a+b)**2+4*rfact*(a+b+a*b)*a*b)-(a+b))
     &                    /(2*a*b)
                     dae  (ie,iout,iin) = rf * a 
                     varia(ie,iout,iin) = rf * varia(ie,iout,iin)
                     varim(ie,iout,iin) = rf * b
                  else if(abs(a).gt.1e-6) then
                     dae  (ie,iout,iin) = rfact * a 
                     varia(ie,iout,iin) = rfact * varia(ie,iout,iin)
                     varim(ie,iout,iin) = 0.
                  else if(abs(b).gt.1e-6) then
                     dae  (ie,iout,iin) = 0.
                     varia(ie,iout,iin) = 0.
                     varim(ie,iout,iin) = rfact * b
                  else
                     dae  (ie,iout,iin) = 0.
                     varia(ie,iout,iin) = 0.
                     varim(ie,iout,iin) = 0.
                  end if
               end if
 1103       continue
 1101    continue
 1100 continue

      ll = lchend(outputfile,80)
      ls = ll
      do 10 while(outputfile(ls:ls).ne.'.')
         ls = ls - 1
 10   continue
      logfil = outputfile(1:ls-1)//'.xlog'

      inquire(file=logfil, EXIST=logfilexist)
      if(logfilexist) then
         return
      else
         open(10, file=logfil, status='unknown')
         write(10,*) 'Modification of interaction from ', filn(1:ll)

         do 1700 iin = 1,9
            do 1701 ie = 1,16
               write(10,'(2i3,10F7.3,a)')iin,ie,
     &              (dae  (ie,iout,iin),iout=1,10),' dx'
               write(10,'(2i3,10F7.3,a)')iin,ie,
     &              (varia(ie,iout,iin),iout=1,10),' da'
               write(10,'(2i3,10F7.3,a)') iin,ie,
     &              (varim(ie,iout,iin),iout=1,10),' dm'
 1701       continue
 1700    continue
         close(10)
         return
      end if
      
 2999 continue
      write(0,*) 'No modification of interaction doe to file error'
      end
