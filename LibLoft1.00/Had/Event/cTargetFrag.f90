module modTargetFrag
!         to be used to get free N, heavy fragment,  remnant fragment
!         for a given projectile, target, and spectator Nucleons from
!         the target.  Target mass # > 56
!    The user may  use only csampFragNHR.
!    test prog. is at the last part
    integer,parameter:: SpecN=21
    integer,parameter:: RestN = 21
    integer,parameter:: dpjA=9 ! p,He,B,C,O,Ne,Si,Ca,Fe
    integer,parameter:: dtgA=8 ! Fe,Cu,Zr,Sn,Xe,Gd,W,Pb
    integer,save:: dpjA2pjA(dpjA)=(/1, 4, 11, 12, 16, 20, 29, 40, 56/)
    integer,save:: dtgA2tgA(dtgA)=(/56, 64, 91, 119, 131, 156, 184,207/)
    real(4):: SpecN2NucAv(SpecN,  dpjA, dtgA)
    real(4):: SpecN2NucRms(SpecN, dpjA, dtgA)
    real(4):: RestN2HeavyNAv(RestN, dpjA, dtgA)
    real(8):: ava(2, 2, 2)  ! average  # of free N or heavy fragment　for two
                  !  sets of spectator nucleon numbers, projectile mass # and target mass #
                  !  in-between which the input values are found.
    real(8):: rmsa(2, 2, 2)  ! same for rms of free N.  not used for heavy fragmnt


    real(8),parameter:: portionOfp= 0.4  ! average proton fraction for A>56
  contains
    subroutine csetFrag
      integer i, j, k
      integer pA, tA
      character(60)::file
      integer:: icon, sN, sRN, dummy
      real(4)::av, rms

      SpecN2NucAv(:,:,:)=0.
      SpecN2NucRms(:,:,:)=0.
      RestN2HeavyNAv(:,:,:)= 0.

      do i = 1, dpjA   ! diff. # of  proj. A
         pA = dpjA2pjA(i)    ! index to progj. A
         do j = 1, dtgA  ! diff # of target A
            tA =  dtgA2tgA(j)  ! index to target. A
!            file=" " ! maka file name
            file = '$COSMOSTOP/Data/Fragment/' ! for implementaion add next after this
            write(file(len(trim(file))+1:), '(a, i0, a, i0,a)') "N-",pA,"-",tA,".parm"

!            write(0,*) file
            call copenf(11,file,icon)
            if(icon /= 0 )  then
               write(0,*) ' could not be opened'
               cycle
            endif
            do while (.true.)
               read(11,*,end=100) av, rms, dummy, sN
               k = (sN-5)/10 + 1
               SpecN2NucAv(k, i, j) = av
               SpecN2NucRms(k, i, j) = rms
            enddo
100         continue
            close(11)
!                read heavy fragment data (only average # is used)
!            file = " "
            file = '$COSMOSTOP/Data/Fragment/' ! for implementaion add next after this
            write(file(len(trim(file))+1:), '(a, i0, a, i0,a)') "H-",pA,"-",tA,".parm"
!            write(0,*) file
            call copenf(11,file,icon)
            if(icon /= 0 )  then
               write(0,*) ' could not be opened'
               cycle
            endif
            do while (.true.)
               read(11,*,end=200) av, rms, dummy, sRN
               k = (sRN-5)/10 + 1
               RestN2HeavyNAv(k,i,j)= av
            enddo
200         continue
            close(11)
         enddo
      enddo
    end subroutine csetFrag

    subroutine cwriteFrag
      implicit none
      integer i, j, k
      do i = 1, dpjA
         do j = 1, dtgA
            write(0,*) 'Nucleons pjA=',dpjA2pjA(i),  ' tgA=', dtgA2tgA(j)
            do k = 1, SpecN
               write(0,'(1p,2g12.3,i4)')  &
                    & SpecN2NucAv(k,i,j), SpecN2NucRms(k,i,j), (k-1)*10 + 5
            enddo
            write(0,*)
         enddo
      enddo

      do i = 1, dpjA
         do j = 1, dtgA
            write(0,*) 'heavy frag; pjA=',dpjA2pjA(i),  ' tgA=', dtgA2tgA(j)
            do k = 1, RestN
               write(0,'(1p,g12.3,i4)')  &
                    & RestN2HeavyNAv(k, i,j), (k-1)*10 + 5
            enddo
            write(0,*)
         enddo
      enddo
    end subroutine cwriteFrag


    subroutine cSampSpecN(N)
!       some very rough example of # of spectator nucleons
!       at target side   for C-W  collisions
      implicit none
      integer,intent(out):: N
      
      real(8):: u, x
      call rndc(u)
      x= log(u)
      if(u > 0.5) then
         N= ( -13.73*x+ 16.12)*x +  184.415
      elseif(u > 0.1 ) then
         N= ( 4.36073*x+ 40.37)*x +  193.026
      elseif( u> 0.01) then
         N = ( 2.3922*x+ 28.462)*x + 175.84
      elseif( u> 0.0001) then
         N = (0.5739*x+ 13.188)*x + 143.53
      elseif( u> 1.e-5) then
         N= (0.9418*x+ 23.482)*x+ 206.60
      else
         N = 60
      endif
    end subroutine cSampSpecN

   subroutine cSampFragNHR(pA, pZ, tA, tZ, sn,  nn, hn, rn, oA, oZ)
!   this is to generate nucleons,light heavy ions and remnant heavy ions
!       for  a given unmber of target spectator nucleon number "sn"
!    Target nuleus mass #, tA, must be > 56 (though 56 can be used)
!    
      implicit none
      integer,intent(in):: pA  ! projectile A
      integer,intent(in):: pZ  ! projectile Z (not used)
      integer,intent(in):: tA  ! target A (>= 56)
      integer,intent(in):: tZ  ! target Z
      integer,intent(in):: sn !  total number of Spectaor nucleons
      integer,intent(out)::nn !  number of outgoing  nucleons among sn
      integer,intent(out)::hn !  number of outgoing  heavy ions (Z=1, 2) among sn
      integer,intent(out)::rn !  number of outgoing remnant ions (high Z) among sn
                              !  this inludes heavy ions Z>2 
      integer,intent(out)::oA(*) ! size must be >= nn+hn+rn. oA(1:nn)  mass # of nucleons =1 
                             !           oA(nn+1:nn+hn) heavy ion mass #
                             !           oA(nn+hn+1:nn+hn+rn): remnant mass #
      integer,intent(out)::oZ(*) ! same as oA for charge.

      logical,save:: first=.true.
      if( first ) then  ! read fragmentation data from Cosmos/DaAta/Fragment/...
         call csetFrag
         first = .false.
      endif

      call cSampFragN(pA, pZ, tA, tZ, sn,  nn,  oA, oZ)
      if(nn < sn) then
!                         nn may be chaged (increased)
         call cSampFragHR(pA, pZ, tA, tZ, sn,  nn, hn, rn, oA, oZ)
      else
         hn = 0
         rn = 0
      endif
    end subroutine cSampFragNHR

    subroutine cSampFragN(pA, pZ, tA, tZ, sn,  nn, oA, oZ)
!   this is to generate nucleons
!       for  a given unmber of target spectator nucleon number "sn"
!    Target nuleus mass #, tA, must be > 56 (though 56 can be used)
!    
      implicit none
      integer,intent(in):: pA  ! projectile A
      integer,intent(in):: pZ  ! projectile Z (not used)
      integer,intent(in):: tA  ! target A (>= 56)
      integer,intent(in):: tZ  ! target Z
      integer,intent(in):: sn !  total number of Spectaor nucleons
      integer,intent(out)::nn !  number of outgoing  nucleons among sn
      integer,intent(out)::oA(*) ! size must be >= nn+hn+rn. oA(1:nn)  mass # of nucleons =1 
                             !           oA(nn+1:nn+hn) heavy ion mass #
                             !           oA(nn+hn+1:nn+hn+rn): remnant mass #
      integer,intent(out)::oZ(*) ! same as oA for charge.
                       !         
      
      integer::i, j, k
      integer:: snindex
      integer:: indexpA1, indextA1, snindex1
      integer:: indexpA2, indextA2, snindex2
      integer:: qsn
      integer:: pAi(2)  !to store two indexes pjA(pAi(1))<= pA < pjA(pAi(2))
                        !       
      integer:: tAi(2)  !to store two indexestgA(tAi(1)) <= tA < tgA(tAi(2))
      integer:: sni(2)  ! to stoe //         for  SpecN2NucRms(snindex, 

      real(8):: ansav(2), ansrms(2)  ! interpolated av and rms
      real(8):: avN, rmsN             ! final answer
      logical::exactpA, exacttA  ! becomes t, if exactly the same mass #  as input is 
                       ! found  in the table.
      real(8):: x0, y0, hx, hy, x, y, temp

      indexpA1 = 0  !             

      do i = 1, dpjA  ! # of diff. proj. A
         if( pA <= dpjA2pjA(i) ) then   
            indexpA1 = i  !  find dpjA2pjA(i) <= pA < dpjA2pjA(i+1)
            exactpA = ( dpjA2pja(i) == pA )
            exit
         endif
      enddo
      if( indexpA1 == 0 ) then
!         write(0,*) ' warning: input pA=',pA, 'too large for ', &
!          & ' cSampFragNHR'
         indexpA1 = dpjA
         exactpA = .true.
      endif
      if(exactpA) then
         indexpA2 = indexpA1
      else
         indexpA2 = indexpA1
         indexpA1 = indexpA1 - 1
         if( indexpA1 < 1) then
            write(0,*) ' indexpA1 =', indexpA1
            write(0,*) 'pA, pZ, tA, tZ, sn'
            write(0,*) pA, pZ, tA, tZ, sn
         endif
      endif

      pAi(1) = indexpA1
      pAi(2) = indexpA2
!/////////////
!      write(0,*) ' pAi=',pAi
!//////////////

!          same for target
      indextA1 = 0
      do i = 1, dtgA
         if( tA <= dtgA2tgA(i) ) then
            indextA1 = i
            exacttA = ( dtgA2tgA(i) == tA )
            exit
         endif
      enddo
      if( indextA1 == 0 ) then
!         write(0,*) ' warning: input pA=',tA, 'too large for ', &
!          & ' cSampFragNHR'
         indextA1= dtgA
         exacttA= .true.
      endif
      if( exacttA) then
         indextA2 = indextA1
      else
         indextA2 = indextA1
         indextA1 = indextA1 -1
         if( indextA1 < 1) then
            write(0,*) ' indextA1 =', indextA1
            write(0,*) 'pA, pZ, tA, tZ, sn'
            write(0,*) pA, pZ, tA, tZ, sn
         endif

      endif
      tAi(1) = indextA1
      tAi(2) = indextA2
!//////////
!      write(0,*) ' tAi=',tAi
!////////////
!      (    nearest # in 5, 15, 25 .. from sn )
      qsn = ((sn-1)/10)*10 + 5   ! quantized sn (5 15 ....)
!                 /5 = 1,2,...
      snindex =(qsn+5)/10   ! index for   SpecN2NucAv(snindex, .) and SpecN2NucRms(snindex,..)
      if(snindex > SpecN ) then
         write(0,*) 'Warning: input sn to cSampFragNHR=',sn,' > SpecN*10=',SpecN*10
      endif
!///////////
!      write(0,*) ' qsn= ',qsn, 'snindex=', snindex
!/////////////
      snindex1 = min(snindex, SpecN)
      if( sn > qsn ) then
         snindex2 = min(snindex+1, SpecN)
      else
         snindex2 = max(1, snindex-1)
      endif
      sni(1) = snindex1
      sni(2) = snindex2
!//////////////////
!      write(0,*) ' sni=',sni
!////////////////
!          examin following combinations. 
!       snindex    indexpA   indextA
!         2            2         2
! 
      do i = 1, 2
         do j = 1, 2
            do k = 1, 2
               ava(k,j,i) = SpecN2NucAv(sni(i), pAi(j), tAi(k))
               rmsa(k,j,i) = SpecN2NucRms(sni(i), pAi(j), tAi(k))
            enddo
         enddo
         x0 = dtgA2tgA(tAi(1))
         y0 = dpjA2pjA(pAi(1))
         hx = dtgA2tgA(tAi(2))- x0
         hy = dpjA2pjA(pAi(2))- y0
         x = min(tA, dtgA2tgA(dtgA))
         y = min(pA, dpjA2pjA(dpjA))
         call c4ptdi(ava(1,1,i),  x0, y0, hx, hy, x, y, ansav(i))
         call c4ptdi(rmsa(1,1,i),  x0, y0, hx, hy, x, y, ansrms(i))
!/////////////
!         write(0,*) 'x0,y0=',x0,y0
!         write(0,*) 'hx,hy=',hx, hy
!         write(0,*) 'x, y=', x,y
!//////////////////
      enddo
!//////////
!      write(0,*)   ' ansav=',ansav
!      write(0,*) ' ansrms=',ansrms
!/////////////
      if( ansav(1) > 0. .and. ansav(2) > 0.) then
         if(ansav(1) /=  ansav(2) ) then
            avN = (ansav(2)- ansav(1))/(sni(2)-sni(1))/10.0 * (sn- qsn) + ansav(1)
            rmsN = (ansrms(2)- ansrms(1))/(sni(2)-sni(1))/10.0 *(sn- qsn) + ansrms(1)
         else
            avN = ansav(1)
            rmsN = ansrms(1)
         endif
      else
         avN = max( ansav(1), ansav(2))
         rmsN =max( ansrms(1), ansrms(1))
      endif
!/////////
!         write(0,*) ' avN=', avN, ' rmsN=',rmsN
!/////////////
      if( avN == 0.) then
         !                  assume all spectator nucleons go into independent nucleons 
         !  nn !  number of outgoing  nucleons among sn
         !  oA(*) ! size must be >= nn+hn+rn. oA(1:nn)  mass # of nucleons =1 
         !           oA(nn+1:nn+hn) heavy ion mass #
         !           oA(nn+hn+1:nn+hn+rn): remnant mass #
         ! oZ(*) ! same as oA for charge.
         nn = sn
      else
         call kgauss(avN, rmsN, temp)
         nn =max( temp, 0.)
         nn = min (nn, sn)
      endif
      if(sn - nn  == 1) then   ! no heavy frag. possible, so all go into free N
         nn = nn + 1
      endif
      call cMkNucFrag(nn, oA, oZ)   ! make nn nucleons
    end subroutine cSampFragN

    subroutine cSampFragHR(pA, pZ, tA, tZ, sn,  nn, hn, rn, oA, oZ)
!   this is to generate light heavy fragments
!       for  a given unmber of non-free  target spectator nucleon number "nn"
!    Target nucleus mass #, tA, must be > 56 (though 56 can be used)
!    
      implicit none
      integer,intent(in):: pA  ! projectile A
      integer,intent(in):: pZ  ! projectile Z (not used)
      integer,intent(in):: tA  ! target A (>= 56)
      integer,intent(in):: tZ  ! target Z
      integer,intent(in):: sn !  total number of Spectaor nucleons
      integer,intent(inout)::nn !  number of outgoing  nucleons among sn
                         ! input vlaues may be changed (increased) to  adjust
                         ! the total number.
      integer,intent(out)::hn !  number of outgoing  heavy ions (Z=1, 2) among sn
      integer,intent(out)::rn !  number of outgoing remnant ions (high Z) among sn
                              !  this inludes heavy ions Z>2 
      integer,intent(inout)::oA(*) ! size must be >= nn+hn+rn. oA(1:nn)  mass # of nucleons =1 
                             !           oA(nn+1:nn+hn) heavy ion mass #
                             !           oA(nn+hn+1:nn+hn+rn): remnant mass #
      integer,intent(inout)::oZ(*) ! same as oA for charge.



      integer:: nfn !  total number of (spectaor nucleons - spectator free nucleons)      
      integer:: pAi(2)  !to store two indexes pjA(pAi(1))<= pA < pjA(pAi(2))
                        !       
      integer:: tAi(2)  !to store two indexestgA(tAi(1)) <= tA < tgA(tAi(2))
      integer:: sni(2)  ! to stoe //         for  SpecN2NucRms(snindex, 
      integer:: nfni(2) ! // for non free nucleon index
      real(8):: ansav(2) ! interpolated av and rms
      real(8):: avN            ! final answer
      logical::exactpA, exacttA  ! becomes t, if exactly the same mass #  as input is 
                       ! found  in the table.

      integer::i, j, k
      real(8):: x0, y0, hx, hy, x, y
      integer:: qnfn, nfnindex, snindex
      integer:: indexpA1, indextA1, snindex1, nfnindex1
      integer:: indexpA2, indextA2, snindex2, nfnindex2
      real(8):: u
      rn = 0
      nfn = sn - nn
      if( nfn <= 1 ) then !  This case will not happen
         write(0,*) ' # of non free nucleons <=', nfn, 'strange in cSampFragHR'
         write(0,*) 'sn =',sn, ' nn =',nn 
         stop
      endif

      if( nfn == 2 ) then
         hn = 1   ! deuteron
         oA(nn+1) = 2
         oZ(nn+1) = 1
         return ! ************
      elseif( nfn == 3 ) then
         hn= 1   !  He3
         oA(nn+1) = 3
         oZ(nn+1) = 2
         return  !   ****
      elseif( nfn == 4) then
         call rndc(u)   
         if( u < 0.25) then
            !  2* deuteron
            do i = 1, 2
               oA(nn+i) = 2
               oZ(nn+i) = 1
            enddo
            hn = 2
            return !   ****
         else
            !  He4
            hn = 1
            oA(nn+1) = 4
            oZ(nn+1) = 2
            return !  ****
         endif
      endif


      indexpA1 = 0
      do i = 1, dpjA
         if( pA <= dpjA2pjA(i) ) then
            indexpA1 = i
            exactpA = ( dpjA2pjA(i) == pA )
            exit
         endif
      enddo
      if( indexpA1 == 0 ) then
!         write(0,*) ' warning: input pA=',pA, 'too large for ', &
!          & ' cSampFragNHR'
         indexpA1 = dpjA
         exactpA = .true.
      endif
      if(exactpA) then
         indexpA2 = indexpA1
      else
         indexpA2 = indexpA1
         indexpA1 = indexpA1 -1
      endif
      pAi(1) = indexpA1
      pAi(2) = indexpA2
!          same for target
      indextA1 = 0
      do i = 1, dtgA
         if( tA <= dtgA2tga(i) ) then
            indextA1 = i
            exacttA = ( dtgA2tgA(i) == tA )
            exit
         endif
      enddo
      if( indextA1 == 0 ) then
!         write(0,*) ' warning: input pA=',tA, 'too large for ', &
!          & ' cSampFragNHR'
         indextA1= dtgA
         exacttA= .true.
      endif
      if( exacttA) then
         indextA2 = indextA1
      else
         indextA2 = indextA1
         indextA1 = indextA1- 1
      endif
      tAi(1) = indextA1
      tAi(2) = indextA2

!      (    nearest # in 5, 15, 25 .. from sn )
      qnfn = ((nfn-1)/10)*10 + 5   ! quantized sn (5 15 ....)
!                 /5 = 1,2,...
      nfnindex =(qnfn+5)/10   ! index for;  RestN2HeavyNAv(RestN, dpjA, dtgA)
      if(nfnindex > RestN ) then
         write(0,*) 'Warning: input nfn to cSampFragH=',nfn,' > RestN*10=',RestN*10
      endif
      nfnindex1 = min(nfnindex, RestN)
      if( nfn > qnfn ) then
         nfnindex2 = min(nfnindex+1, RestN)
      else
         nfnindex2 = max(1, nfnindex-1)
      endif
      nfni(1) = nfnindex1
      nfni(2) = nfnindex2
      
!          examin following combinations. 
!       nfnindex    indexpA   indextA
!         2            2         2
! 
      do i = 1, 2   ! nfn
         do j = 1, 2   ! pj
            do k = 1, 2  ! tg
               ava(k,j,i) =  RestN2HeavyNAv( nfni(i), pAi(j), tAi(k) )
            enddo
         enddo
         x0 = dtgA2tgA(tAi(1))
         y0 = dpjA2pjA(pAi(1))
         hx = dtgA2tgA(tAi(2))- x0
         hy = dpjA2pjA(pAi(2))- y0
         x = min(tA, dtgA2tgA(dtgA))
         y = min(pA, dpjA2pjA(dpjA))
         call c4ptdi(ava(1,1,1),  x0, y0, hx, hy, x, y, ansav(i))
      enddo

      if( ansav(1) > 0. .and. ansav(2) > 0.) then
         if(ansav(1) /=  ansav(2) ) then
            avN = (ansav(2)- ansav(1))/(nfni(2) - nfni(1))/10.  * (nfn-qnfn) + ansav(1)
         else
            avN = ansav(1)
         endif
      else
         avN = max(ansav(1), ansav(2))
      endif
      if( avN == 0.) then
         hn = 0
         oA(nn+1) = nfn ! > 4
         oZ(nn+1) = oA(nn+1)*portionOfp
         rn = 1
      else
         call kpoisn(avN, hn)
         call cMkHFrag(tA, tZ, sn, nn, hn, rn, oA, oZ)
      endif 
    end subroutine cSampFragHR


    subroutine cMkNucFrag(nnuc, A, Z)
      implicit none
      integer,intent(in):: nnuc
      integer,intent(out)::A(*) ! always 1
      integer,intent(out)::Z(*) ! 0 or 1
     

      real(8):: u
      integer:: i
      do i = 1, nnuc
         A(i) = 1
         call rndc(u)
         if( u < portionOfp ) then 
            Z(i) = 1
         else
            Z(i) = 0
         endif
      enddo
    end subroutine cMkNucFrag

    subroutine cMkHFrag(tA, tZ, sn, nn, hn, rn, oA, oZ)
      implicit none
      integer,intent(in)::tA  ! target mass #
      integer,intent(in)::tZ  ! target charge
      integer,intent(in)::sn  ! # of spectator nucleons
      integer,intent(inout)::nn  ! # of  free nucleons 
      integer,intent(inout)::hn  ! # of hevy fragment
      integer,intent(out)::rn  ! # of  remnant hevy fragment
      integer,intent(inout)::oA(*)  !
      integer,intent(inout)::oZ(*) 
     
      real(8),parameter:: portionOfd= 0.25d0
      real(8),parameter:: portionOfHe3= 0.25d0 + portionOfd
      real(8),parameter:: portionOfHe4= 0.5d0 + portionOfHe3
      real(8),parameter:: portionOfp = 0.4
      real(8):: u
      integer:: i, nad, rnA, hA
      integer:: sumA
      integer:: nfn  ! # of non free nucleons
      integer(2):: tempA(100), tempZ(100)


      nfn = sn - nn
      sumA = 0
      do i = 1, hn
         call rndc(u)
         if(u < portionOfd ) then
            sumA = sumA + 2
            tempA(i) = 2
            tempZ(i) = 1
         elseif( u < portionOfHe3) then
            sumA = sumA + 3
            tempA(i) = 3
            tempZ(i) = 2
         else
            sumA = sumA + 4
            tempA(i) = 4
            tempZ(i) = 2
         endif
         if( sumA >= nfn ) then
            hn = i
            exit
         endif
      enddo
      nad = sumA- nfn

      if( nad > 0 ) then
         hA = tempA(hn) - nad  ! # of extra nucleons
         if( hA == 1 ) then
            hn = i-1
            nn = nn + 1
            oA(nn) = 1
            call rndc(u)
            if( u < portionOfp ) then
               oZ(nn)  = 1
            else
               oZ(nn) = 0
            endif
         else
            tempA(hn) = hA
            if( hA== 2 ) then
               tempZ(hn) = 1
            else
               tempZ(hn) = 2
            endif
         endif
         sumA = sumA - hA
      endif
      rnA = max(0, nfn - sumA)    ! remnant A
      if( rnA <= 4 ) then    ! make them nucleons
         do i = 1, rnA
            nn = nn + 1
            oA(nn) = 1
            call rndc(u)
            if( u < portionOfp ) then
               oZ(nn)  = 1
            else
               oZ(nn) = 0
            endif
         enddo
      endif
      do i = 1, hn              ! move heavy fragment into oA,oZ
         oA(nn+i) = tempA(i)
         oZ(nn+i) = tempZ(i)
      enddo
      if( rnA > 4 ) then       ! make 1  remnantt
         oA(nn + hn + 1) = rnA
         oz(nn + hn + 1) = rnA*portionOfp
         rn = 1
      else
         rn = 0
      endif
    end subroutine cMkHFrag
!     ****************************************************************
!     *                                                              *
!     * c4ptdi:  4-point two dimensional interpolation              *
!     *   special for cSampFragNHR
!     *          each of two points maybe the same points
!     *   if function vlaue is 0, not used.
!     ****************************************************************
!
!   /usage/
!
!       call c4ptdi(f, x0, y0, hx, hy, x, y, ans)
!
!     f containes the function values at (x0,y0), (x0+hx, y0+hy)
!     ans gets the value of the funtion at (x,y)
!     linear interpolation is tried.
!     hx can be 0, if x is the same as x0
!     hy can be 0, if y is the same as y0
!     if f is 0, only non zero value is used.
!     if all f is 0, ans=0.
      subroutine c4ptdi(f,  x0, y0, hx, hy, x, y, ans)
      implicit none
      real*8 f(2,2),  x0, y0, hx, hy, x, y, ans
!
      integer i, j, i1, j1
      real(8):: a, b, pL, pH, q, pL1, pH1, q1
      real(8):: fij, fi1j, fij1, fi1j1, fa, fb

      if(hx /= 0.0) then
         a=(x-x0)/hx
      else
         a = 0.
         if(x /= x0) then
            write(0,*) ' x=',x, ' != x0=',x0,' while hx=0'
            write(0,*) ' in c4ptdi'
            stop
         endif
      endif
      if( hy /= 0.0) then
         b = (y-y0)/hy
      else
         b = 0.
         if(y /= y0) then
            write(0,*) ' y=',y, ' != y0=',y0,' while hy=0'
            write(0,*) ' in c4ptdi'
            stop
         endif
      endif

      i=1
      j=1
      pL=a
      pH=a  
      pL1=1.-pL
      pH1=1.-pH
      q=b
      q1=1.-q
      if( hx == 0.) then
         i1 = i
      else
         i1 = i + 1
      endif
      if( hy  == 0.) then
         j1 = j
      else
         j1 = j + 1
      endif
      fij = f(i,j)
      fi1j= f(i1,j)
      fij1= f(i,j1)
      fi1j1= f(i1,j1)  

      if( fij == 0.) then
         pL = 1.
         pL1 =0.
      endif
      if( fi1j == 0.) then
         pL1 = 1
         pL =  0.
      endif
      if( fij1 == 0. ) then
         pH = 1.
         pH1 = 0.
      endif
      if( fi1j1 == 0.) then
         pH = 0.
         pH1 = 1.
      endif
      fa =  fij*pL1 + fi1j*pL
      fb =  fij1*pH1 + fi1j1*pH 

      if( fa == 0.) then
         q1 =0.
         q = 1.
      elseif(fb == 0.) then
         q1 = 1.
         q = 0.
      endif
      ans=  fa  * q1 +  fb* q
   end subroutine c4ptdi



  end module modTargetFrag
  
!program main
!  use  modTargetFrag
!  implicit none
!  integer N, i, j
!  integer:: oA(210), oZ(210)
!  integer:: nn, hn, rn 
!!  call csetFrag   ! this is now inside cSampFragNHR
!!  call cwriteFrag  ! so this cannot be here.
!  do i = 1, 100
!     call cSampSpecN(N)
!     call  cSampFragNHR(12, 6, 184, 74, N,  nn, hn, rn, oA, oZ)
!     write(0,*) '*****', N, nn, hn, rn
!     do j = 1, hn
!        write(0,*)  oA(nn+j), oZ(nn+j)
!     enddo
!     do j = 1, rn
!        write(0,*) oA(nn+hn+j), oZ(nn+hn+j)
!     enddo
!  enddo
!
!end program main
