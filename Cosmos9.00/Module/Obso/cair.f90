  module modAir
    implicit none
    integer,parameter:: noOfElem=3
    real(8),save:: elemA(noOfElem) = (/14.0d0, 16.0d0, 40.0d0/)
    real(8),save:: elemZ(noOfElem) = (/7.0d0, 8.0d0, 18.0d0/)
    real(8),save:: No(noOfElem) =  ( /1.56d0, 0.42d0, 0.01d0/)
    real(8),save:: frac(noOfElem) 
    real(8),save:: TargetMassN, TargetAtomicN, TargetZ2, & 
                   TargetZ2_3rd, TargetZ1_3rd
!       (A,Z) of  selected target nucleus
    real(8),save:: TargetNucleonNo, TargetProtonNo 


    contains
      subroutine cinimodAir
      implicit none
      real(8):: sumN
      integer:: i, j

      sumN = sum(No(:))
      TargetMassN =  sum(No(:)*elemA(:))/sumN
      TargetAtomicN =sum(No(:)*elemZ(:))/sumN
      TargetZ2 = sum(No(:)*elemZ(:)**2)/sumN
!              <Z>^0.666 = 3.7464 <Z^0.666>=3.7457
      TargetZ2_3rd = TargetAtomicN**0.66666
!              <Z>^0.333 =1.9367 <Z^0.3333>=1.9347
      TargetZ1_3rd = TargetAtomicN**0.33333
      frac(:) = No(:)/sumN
!//////////
      write(0,*) "<A>=", TargetMassN
      write(0,*) "<Z>=", TargetAtomicN 
      write(0,*) "<Z^2>=", TargetZ2
      write(0,*) "<Z>^2=", TargetAtomicN**2
      write(0,*) "<Z>^(2/3)=", TargetZ2_3rd
      write(0,*) "<Z>^(1/3)=", TargetZ1_3rd
!//////////// 
      TargetNucleonNo = 14   !  these two will be updated 
      TargetProtonNo =7      ! for each collision


      do j = 1, noOfElem
         shp = xsecmin          ! (say 10mb)
         do i = 1, nxsec
            call cxp2xAXsec(elemA(j), shp, elemS(j,i))
            shp = shp + xsecstep
         enddo
      enddo


      do i = 1, nxsec
         sigma(i) = sum(No(:)*elemS(:,i))
      enddo

      do i = 1, nxsec
         do j = 1, noOfElem
            cumsigma(j, i) =  No(j)*elemS(j,i)/sigma(i)
         enddo
      enddo

      do i = 1, nxsec
         do j = 2, noOfElem
            cumsigma(j, i) = 
     *         cumsigma(j, i) + cumsigma(j-1, i)
         enddo
c          for safety
         cumsigma(noOfElem, i) = 1.0
      enddo
    end subroutine cinimodAir
    subroutine  cfixTarget(xs)
      implicit none
      integer:: j
c             xs   xs                        xs
c      jcon = 1    0                         0
c       i          1                         nxsec
      call kfrge(sigma, 1, nxsec,  xs,  i, jcon)
      if(jcon .ne. 0 ) then
         mbindex = 1
      elseif( i .lt. nxsec - 1) then
         mbindex = i + 1
      else
         mbindex = nxsec
      endif

      call rndc(u)
      do j = 1, noOfElem
         if(u .le. cumsigma(j, mbindex))  goto 10
      enddo
 10   continue
      TargetNucleonNo = elemA(j) 
      TargetProtonNo =elemZ(j)
      if( ActiveMdl  == 'jam' .and. 
     *    TrackBefMove.p.code /= kmuon ) then
c          we need xs on the target nucleon for bmax; for other models
c          need not do this
         tgA =  TargetNucleonNo 
         tgZ =  TargetProtonNo
c         if(TrackBefMove.p.code <= knuc  .and. 
c     *     (TrackBefMove.p.fm.p(4)-TrackBefMove.p.mass)  < 10.d0) then
         if(TrackBefMove.p.code <= knuc ) then
            if( JamXs == 1 ) then
               call ctotx(TrackBefMove.p, tgA,  TargetXs)
            elseif( JamXs == 0 ) then
               call cinelx(TrackBefMove.p,  tgA, tgZ, TargetXs)
            else
               write(0,*) ' JamXs =',JamXs,' not usable in cixsec'
               stop
            endif
         else
            call cinelx(TrackBefMove.p, tgA, tgZ, TargetXs)
         endif
      endif
c///////////////////
c      write(0,*) ' A, Z=', TargetNucleonNo, TargetProtonNo 
c///////////////
end subroutine cfixTarget
end module modAir
