
! Suppose a detector consisting of many components:
! a subdetector A contains subdetectors B each of which
! intern contains components C (we call this target;
! say SciFi; ) to measure  energy  deposit. 
! The detector contains  A('s).
! A's  may be contained by another larger component, 
! B('s) may be contained by A indirectory. C('s) also may 
! be contained indirectly by B. 
!*** But to specify target C's, the subdetector hierarchy
!    A B C should be  enough.
!    i.e., if there are other C's in the detector not
!    belonging to this hierarchy, they should not be mixed.
!    So, if only 1 A exists in the detecctor and if  B C
!    exist not contained by A, and such C is not a target,
!    A must be specified. If no such B C exist, A need not
!    be sepcified.
!*** At first, for simplicity, we treat the case where  C
!  is a simple component not  containing other components.
! 
! Let's suppose # of A, B, C is nA, nB, nC respectively.
! In this example, we can inform the name of array, say,
! defined as
!   integer:: abc(nA, nB, nC)
! to epGetIndex. After that, you can use this 'abc'
! to get the component number, cmpn,  of a particular C as
!    cmpn=abc(i, j, k)
! where i is the i-th A and j the j-th B  and k the k-th C. 
! By referring to cmpn you can get various quantities related
! to the component. ( i in (1,nA), j in (1,nB), k in (1,nC)).
! Maximum of hierarchy is 7  (normally, larger rank
! array is not permitted in Fortran).  In normal applications, 
! 2 or 3 is enough. 
!  
!       Generally nA, nB, nC can be known if one counts componens
!  in a  config file carefully but sometimes it is not so easy.
!  So a counter routine is prepared which is a little  bit
!  modified version of epGetIndex. A typical usage will be:
!  (in, say,  uiaev in ephook.f)
!
!     integer,allocatable:: abc(:,:,:)   
!     character(len=12)::  spec='A B C'
!     integer::n, n1, n2  ! n should become 3
!     integer::shape(3)
!     
!          get size of  A,B,C and fill them in shape
!     call epCountSubdTree(spec,  shape, n, n1, n2)
!          allocate abc with that size
!     allocate(abc(shape(1), shape(2), shape(3)))
!          fill abc
!     call  epGetIndex(spec, shape, abc, n1, n2)
!
!       After this, the user can know the compoent number
!     by abc(i, j, k) for i-th A, j-th B and k-th C.
!
      module modGetIndex
        implicit none
        integer,parameter:: maxf=7  ! max hierarchy
        integer,parameter:: maxOdd=2 ! max Odd subds
                      ! making maxOdd > 2 is very difficult
                      !  
        character(len=16):: splitSpec(maxf) ! name list of
                   !  sub detectors
         ! we permit A B(B') C type spec. One upper hierarchy
         ! than target  may be B' instead of B.
         !  To store A B' C next is used.
        integer,parameter:: maxbr=4 ! max possilbe top branches
        integer:: ntbr  ! # of top branches actully used
        character(len=16):: topbranch(maxbr)  ! branch names
!     　　topbranch specification is used to combine two or 
!        more (upto 4) hirarchy path, say,
!         A1   B   C
!         A2   B   C
!       if we get index for A1 B C,  array abc1(:,:,:)
!       is needed.      For A2  B  C,   array abc2(:,:,:)
!       If they are tightly connected, one may want to
!       use index array abc(2, :,:,:). abc(1,..) can play
!       a role of abc1 and   abc(2,...)  abc2(..).
        character(len=16):: OddSubdName(maxOdd) ! to store names in
!         e.g,  (B,B') 
        integer:: OddSubd  ! if B' exists > 1 else 1.
        integer,save:: nitem  ! actual number of items in splitSpecs
        integer,save:: filled(maxf)  ! # of componenets filled
                                    ! for each sub detector

        integer,save:: filledmax(maxf) ! to keep max of filled

        integer,save:: depthList(maxf) ! hierarchy depth list
        character(len=16),parameter::targetdef='simple'
        integer,parameter:: dEdef=-100
        integer,parameter:: IOdef=-100

        character(len=16),save::targetcpy=targetdef
        integer,save:: dEcpy=dEdef
        integer,save:: IOcpy=IOdef
        contains


      subroutine epGetIndex(spec, shape, idarray, 
     * num, target, dE, IO, judgeBy)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
       ! Normally upto idarray may be given.
      character(*),intent(in)::spec ! sub detector 
       ! hierarchy specification. From larg to small
       ! In the above example,
       !    A B C
       ! This is default format which assumes that target
       ! C is the  media name of a simple component (i.e,
       ! it does not contain other components).
       ! A B must be subdetector name. For non default 
       ! specification, see a separate document.
       ! Let n be the number of items in spec.
       ! In the above example, n=3 
      integer,intent(in):: shape(*)
       ! The size must be >=n. 1:n is used. 
       ! In the case of above example,
       ! shape(1): the number of A's 
       ! shape(2): the number of B's
       ! shape(3): the number of C's  ...
      integer,intent(out):: idarray(*)
         ! idarray(N1,N2,...,Nn)  where Ni=shape(i).
         ! In our example, this is  abc.
         ! It is treated as rank 1 array here so the
         ! user also may define it as a rank 1 array with
         ! the size s >= N1*N2*..Nn
         ! ie, integer:: idarray(s).
         ! In that case,  after calling this, the user must
         ! reshape the idarray.  The user must define an array
         !    integer:: somearray(N1, N2,...)
         ! and fill it by
         !     somearray=reshape(idarray, shape)
         ! Then, the user can get a component number via
         ! somearray(i,j...)
         ! 
         ! This "big to small" hierarchical approach generates
         ! C-language style indexing of array ingredients. 
         ! However, we are dealing with Fortran, so, 
         ! in our example, abc(i,j,k) and abc(i,j,k') (k~k') are
         ! probably located rather far distant point in the
         ! memory space while they would be accessed by changing
         ! k continuously. This will lead to slower memory access
         ! than they were located in a consecutive memory space.
         ! This overhead may be negligible but in some case, if
         ! the user access abc frequently, may not be desirable. 
         ! In such a case, the user  may reshape the array by
         ! using reshape function as
         !    integer:: Farray(nC,nB,nA)
         !    Farray = reshape( abc, (/nC,nB,nA/), order=(/3,2,1/))
         ! Then, the user can use Farray to get a component
         ! number as Farray(k, j, i).  Farray(k, j,i) (k=1,nC) 
         ! are located consecutively in memory space.
      integer,pointer,intent(out)::num(:,:)
              ! 1st dim. is N1 if |A1,A2..| spec is used
              !     else 1. 
              ! num(1,1)  # of C's for A B  C  hier.
              ! num(1,2)  # of C's for A B' C
              ! if |A1,A2,A3| B(B') C hier,
              ! num(i, 1) is # of C's for Ai B C
              ! num(i, 2) is # of C's for Ai B' C
      character(*),intent(in), optional::target
        ! This is optional and probably rarely needed.
        ! The user may use next format, if 'spec' is not enough
        ! to fix C in the default mode.
        !
        ! target=xxxx where xxxx is a character strig and one of
        !
        !  'subd'     : C is not medium name but subdetector name
        !
        !       next 3 are for C being medium name.
        !  'container': C is a medium name but it contains other
        !               component(s), i.e, a subdetector.
        !  'any'      : C is a medium name but it can be either 
        !               a simple component or container
        !  'simple'   : C is a medium name and it is a simple
        !               component. 
        !               This is default so it need not be given
        !               but acceptable.
      integer,intent(in),optional:: dE  ! target has specifcation
                    ! of this value for energy loss count.
                    ! dE={-2,-1,0, 1, 2}
      integer,intent(in),optional:: IO  ! target has specifcation
                    ! of this value for In/Out count.
                    ! IO={0,1,2,3}
      ! If both dE and IO are given, .or. is taken.
      logical, external,optional:: judgeBy ! If C cannot be fixed
             ! by using spec, target, dE, IO, the user may
             ! give a logical function as judgeBy=myfunc
          ! where myfunc represents a function name which
          ! gives logical T/F. If F, current C is not accepted.
          ! If T, current C may be accepted if other conditoins
          ! are OK.  
          ! myfunc (--name may be arbitrary--) should be specified
          ! as external in a program calling this subroutine.
          ! Its arguments are (
          ! and can use module "modGetIndex"
          ! It must be defined as,
          !  function myfunc( .. )  result(ans)
          !  use  ....
          !  implicit none
          !  logical:: ans
          ! ...
          !  ans = ...
          !  end function myfunc

      integer:: arraysize       ! min array size needed for idarray
      integer:: subdc    !  hierarchy counter
      integer:: prod   ! working variable to get some product
      integer:: lp, rp  ! used to see ( )
      
      integer:: ntbrc, n, OddSubdc 
      integer:: occupied

      
!         ntbr, splitSpec, topbranch, OddSubd, nitem are fixed
      call epProcSpecForIdx(spec)

      if( ntbr >= 2 ) then
         n = nitem + 1
      else
         n = nitem
      endif
      arraysize = product(shape(1:n))  ! not nitem
      ! clear idx area    
      idarray(1:arraysize) = 0

      if( PRESENT(target) ) then
         targetcpy=target
      endif
      if( PRESENT(dE) ) then
         dEcpy = dE
      endif
      if( PRESENT(IO) ) then
         IOcpy = IO
      endif
      if(associated(num)) deallocate( num )
      allocate( num(ntbr,OddSubd) )
      num(:,:) = 0  ! clear counter
      do ntbrc = 1, ntbr
         splitSpec(1) = topbranch(ntbrc)
         occupied = 0
         do OddSubdc = 1, OddSubd
            filledmax(1:nitem) = 0
            filled(1:nitem) = 0
            depthList(1:nitem) = 0
            if( nitem > 1 .and. OddSubd > 1 ) then
               splitSpec(nitem-1) = OddSubdName(OddSubdc)
            endif

            subdc = 1           ! search from the hierarchy top
            call epGetIndex0(0, Det%nct, subdc, shape, ntbrc,
     *      idarray(occupied+1), arraysize-occupied, judgeBy)

            num(ntbrc, OddSubdc) = filledmax(nitem) 
            occupied = num(ntbrc,1)*product(shape(1:nitem-1))
         enddo
      enddo

!        restore default
      targetcpy=targetdef
      dEcpy=dEdef
      IOcpy=IOdef
      
      end subroutine epGetIndex
      
      subroutine epProcSpecForIdx(specin)
      implicit none
      character(len=*),intent(in):: specin

      character(len=len(trim(specin)) )::spec  ! copy of specin

      integer:: i,  icon, i1, i2, nc
      character(len=16)::tempc(2)
      logical error

      spec = trim(specin)

      i1 = index(spec, '|')
      if(i1 > 1 .and. spec(1:i1-1) /= ' ') then
         write(0,*) 'Indexing input error: spec=',specin
         write(0,*) '| is used non-first items'
         stop
      endif
      if(i1 > 0 ) then
         i2 = index(spec(i1+2:),'|')
         if(i2==0 ) then
            write(0,*) ' error usage of | for Indexing '
            write(0,*) ' the input is :', specin
            stop
         endif
!                            i1             i1+i2+1
!           treat input like |abc, xyaz,xxxx|  and
!           extract abc xyaz xxx to store them in
!           topbranch  if | abc| or |abc,| ntbr=1
         call ksplit2(spec(i1+1:i1+i2), ',', topbranch, maxbr, ntbr,
     *   icon)
         if(icon /= 0 ) then
            write(0,*) 'For indexing, error from ksplit2 icon=',icon
            write(0,*) 'NG in |... | where ... is ', spec(i1+1:i1+i2)
            stop
         endif
         spec(i1:i1+i2+1) = topbranch(1) ! fill only first one in spec
      else
         ntbr = 1   ! ntbr = 1 should be treaed specially
      endif


      i1 = index(spec, '(')
      if(i1 > 0 ) then
         i2 = index(spec(i1+2:),')')
         if(i2==0 ) then
            write(0,*) ' error usage of ( ) for Indexing '
            write(0,*) ' the input is :', specin
            stop
         endif
         ! check () is items last but one; after ), there
         ! must be one blank and non blank.
         if( i1 + i2 +2  <= len(spec)-1 ) then
            call kgetField(spec(i1+i2+2:), tempc, 2, nc)
            error = nc /= 1
         else
            error = .true.
         endif
         if( error ) then
            write(0,*) 'Indexing spec. error; () position NG'
            write(0,*) ' input is :',  specin
            stop
         endif
               
         call ksplit2(spec(i1+1:i1+i2), ',', OddSubdName, maxOdd, 
     *   OddSubd,icon)
         if(icon /= 0 ) then
            write(0,*) 'For indexing, error from ksplit2 icon=',icon
            write(0,*) 'NG in (... ) where ... is ', spec(i1+1:i1+i2)
            stop
         endif
         spec(i1:i1+i2+1) = OddSubdName(1) ! fill only first one in spec
      else
         OddSubd = 1  !   should be treated  specially
      endif
      
      call kgetField(spec, splitSpec, maxf, nitem)
      if( nitem > 1 ) then
         OddSubdName(1) = splitSpec(nitem-1) ! if OddSubd==1, must fill
      endif
      topbranch(1) = splitSpec(1)
      end subroutine epProcSpecForIdx

      

 
      recursive subroutine epGetIndex0(mother, i,  subdc,  shape,
     *    ntbrc,  idarray, size, judgeBy)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer,intent(in)::i     ! comp. #
      integer,intent(in)::mother ! subd index of the 
                !  mother component of comp. # i
                !  if  i is  world, 0

              ! search subd or
              !  media specified by splitSpec(subdc)
      integer,intent(inout):: subdc 
              ! # of subd's in splitSpec(subdc)
      integer,intent(in):: shape(*) ! 1:nitem is used.
              !  idaaray size 
      integer,intent(in):: ntbrc   ! branch counter
      integer,intent(in):: size 
      integer,intent(out):: idarray(size) !
      logical,external,optional:: judgeBy

!      Below  'save' should not be used. If used,
!      say, flag becomes 0 which should be non 0.
       integer:: depth
       integer:: j, m, k, ii, l, prod
       character(len=16):: sname ! subd name
       integer :: msubdi, flag, ll, namec, pos
       integer :: tempc, nornf, fsubdc
        ! fsubdc was logical so far( June 17,2016 why was ok?)
       logical :: ok
       depth = Det%cmp(i)%level    ! hier depth
       j = Det%cmp(i)%NMatreska    ! # of contained comp.
       tempc = subdc   ! temporary counter
       do while( tempc > 1 ) 
          ! if the hier depth of i is higher than
          ! one more upper hier of current hier.
          ! we should search higher hier. 
          ! repeat it.
          if( depth <= depthList(tempc-1) ) then
             tempc = tempc - 1
             depthList(tempc:nitem) = 0
          else
             exit
          endif
       enddo
          ! see if depth is really higher
       if( subdc /= tempc) then
!         goto deeper hierarchy.
!         before that, see fill is complete
          fsubdc = subdc
          if( ntbr >= 2 ) fsubdc = subdc +1

          if( filled(subdc) /= shape(fsubdc)) then
               ! this should not happen in principl. 
               ! for safety check next
             if(OddSubd <= 1 .and.
     *            filled(subdc) < shape(fsubdc)) then
!                   ! warning
                write(0,*) '# of media or subd: ', splitSpec(subdc),
     *               '= ', filled(subdc), ' < ', shape(fsubdc)
             elseif(filled(subdc) > shape(fsubdc)) then
!                   !  error
                write(0,*) '# of media or subd: ', splitSpec(subdc),
     *               '= ', filled(subdc), ' > ', shape(fsubdc)
                stop
             endif
          endif
!!           clear lower  counter
          filled(subdc:nitem) = 0
          subdc = tempc   ! go to upper hier.
       endif

           ! get media name of i-th comp.
       call epqSubdMediaName(mother, i, sname, subdc, flag)
!              additional check by user
       if( subdc == nitem ) then
          call epGetIndexUserChk(mother, i, flag, judgeBy)
       endif

       if( flag == 0 .and. ( sname == splitSpec(subdc) )) then

!              found a subd. fill the counter

          filled(subdc) = filled(subdc) + 1
          filledmax(subdc) = max( filledmax(subdc), filled(subdc) )
          depthList(subdc) = depth
          if( subdc == nitem ) then  ! see if last node
 !           if so fix the location 'pos' in rank 1 idarray
 !          (to understand, think rank 2 or 3 case.)
             pos = 0
             prod = 1
             if( ntbr <= 1 ) then
                do k = 2, nitem
                   prod = prod * shape(k-1)
                   pos = pos + (filled(k)-1)*prod 
                enddo
                pos = pos + filled(1)
             else
                do k = 2, nitem+1
                   prod = prod * shape(k-1)
                   pos = pos + (filled(k-1)-1)* prod
                enddo
                pos = pos  + ntbrc
             endif
                  ! safety check
             if(pos <= 0 .or. pos > size) then
                write(0,*) 'In epGetIndex:'
                write(0,*) " filled(1:nitem)=", filled(1:nitem)
                write(0,*) ' pos is strange=', pos, ' i =',i
                stop
             else
                idarray(pos) = i
             endif
          else
!             not found yet;    see deeper subd.
             subdc = subdc + 1
          endif
       endif                          
       call epqSubdIdx(i, msubdi)
         ! msubdi is a  subdetector index to which         
         ! the component i belongs (if compn                
         ! is a simple comp.)                             
         ! If it is a subdetector, msubdi becomes            
         ! its subd index.                                
         ! if i is invalid, -1                        
         ! if i  does not belong to a subd, 0 
       do m = j, 1, -1
                 ! for partially contained one, CnArea.... < 0
          ii =abs( CnArea( Det%cmp(i)%ContainsR+m ) )
          call epGetIndex0(msubdi, ii, subdc, shape, ntbrc,
     *      idarray, size)
       enddo

       end subroutine epGetIndex0

       subroutine epGetIndexUserChk(mother, i, flag, judgeBy)
!              additional check by user       
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
       integer,intent(in):: mother, i
       integer,intent(inout):: flag
       logical,external,optional:: judgeBy

       logical::ok

          if(dEcpy /= dEdef .or. IOcpy /= IOdef ) then
             ok =(dEcpy ==  Det%cmp(i)%countDE)  .or.
     *            (IOcpy ==  Det%cmp(i)%countIO)
          else
             ok = .true.
          endif
          if( .not. ok ) then
             flag = 1 
          else
             if( PRESENT( judgeBy ) ) then
                ok = judgeBy(mother, i)
             endif
             if(.not. ok ) flag = 1   
          endif
       end subroutine epGetIndexUserChk

       subroutine epCountSubdTree(spec,  shape, n, num,
     *  target, dE, IO, judgeBy)

!!!      use modGetIndex
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      character(*),intent(in)::spec ! sub detector 
       ! hierarchy specification. From larger to
       ! smaller.  In the case of the example given in
       ! epGetIndex.f
       !    'A B C'
      integer,intent(out):: shape(*)
       ! The size must be >=n. 1:n is used. n is # of items in spec
       ! In the case of above example,
       ! shape(1): the number of A's 
       ! shape(2): the number of B's
       ! shape(3): the number of C's  ...
      integer,intent(out):: n  ! see above. so the user knows it
            ! but we return it. But if |A1,A2..| type spec is used,
            ! 1 is added.  
      integer,pointer,intent(out)::num(:,:)
              ! 1st dim. is N1 if |A1,A2..| spec is used
              !     else 1. 
              ! num(1,1)  # of C's for A B  C  hier.
              ! num(1,2)  # of C's for A B' C
              ! if |A1,A2,A3| B(B') C hier,
              ! num(i, 1) is # of C's for Ai B C
              ! num(i, 2) is # of C's for Ai B' C

      character(*),intent(in), optional::target
      integer,intent(in),optional:: dE  
      integer,intent(in),optional:: IO  
      logical, external,optional:: judgeBy 
      integer:: ntbrc, OddSubdc

      integer:: subdc    !  hierarchy counter
      integer:: prod   ! working variable to get some product
       ! split spec into nitem parts in splitSpec

      call epProcSpecForIdx(spec)

      if( ntbr >= 2 ) then
         n = nitem + 1
      else
         n = nitem
      endif

      if( PRESENT(target) ) then
         targetcpy=target
      endif
      if( PRESENT(dE) ) then
         dEcpy = dE
      endif
      if( PRESENT(IO) ) then
         IOcpy = IO
      endif
      if(associated(num)) deallocate( num )
      allocate( num(ntbr,OddSubd) )
      num(:,:) = 0  ! clear counter

      if( ntbr >= 2 ) then
         shape(2:n) = 0
         shape(1) = ntbr
      else
         shape(1:n) =0
      endif
      do ntbrc = 1, ntbr
         splitSpec(1) = topbranch(ntbrc)
         do OddSubdc = 1, OddSubd
            filledmax(1:nitem) = 0
            filled(1:nitem) = 0
            depthList(1:nitem) = 0
            subdc = 1           ! from the top 

            if( nitem > 1 .and. OddSubd > 1 ) then
               splitSpec(nitem-1) = OddSubdName(OddSubdc)
            endif

            call epCsubdTree0(0, Det%nct, subdc, judgeby)
            num(ntbrc, OddSubdc) = filledmax(nitem)
            if( ntbr == 1 ) then
               shape(1:nitem-1) =
     *           max(filledmax(1:nitem-1), shape(1:nitem-1))
               if( OddSubd >= 2 ) then
                  shape(nitem) = shape(nitem) + filledmax(nitem)
               else
                  shape(nitem) = max(shape(nitem), filledmax(nitem))
               endif
            else
               shape(2:n-1) =
     *           max( filledmax(1:nitem-1), shape(2:n-1))
               if( OddSubd >= 2 )  then
                  shape(n) = shape(n) + filledmax(nitem)
               else
                  shape(n) =  max( filledmax(nitem), shape(n))
               endif
            endif
         enddo
      enddo

 !          restore default
      targetcpy=targetdef
      dEcpy=dEdef
      IOcpy=IOdef


      end subroutine epCountSubdTree


      recursive subroutine epCsubdTree0(mother, i, 
     *       subdc, judgeBy)
!!!       use  modGetIndex
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
       integer,intent(in)::i  ! comp. #
       integer,intent(in)::mother ! subd index of the 
                !  mother component of comp. # i
                !  if  i is  world, 0

       integer,intent(inout):: subdc  ! find subd or
              !  media specified by splitSpec(subdc)
       logical,external,optional:: judgeBy

!      Below  save should not be used. some strange things happen
!       say, flag becomes 0 which should be non 0 and sname
!       becomes wrong. 
       integer:: depth
       integer:: j, m, k, ii, l, prod
       character(len=16):: sname ! subd name
       integer :: msubdi, flag, ll, namec, pos
       integer :: tempc

       depth = Det%cmp(i)%level
       j = Det%cmp(i)%NMatreska
       tempc = subdc
       do while( tempc > 1 ) 
          if( depth <= depthList(tempc-1) ) then
             tempc = tempc - 1
             depthList(tempc:nitem) = 0
          else
             exit
          endif
       enddo

       if( subdc /= tempc) then
!         goto upper hierarchy.
          filled(subdc:nitem) = 0
          subdc = tempc
       endif
       call epqSubdMediaName(mother, i, sname, subdc,flag)

!            additional check by user  
       if( subdc == nitem ) then
          call epGetIndexUserChk(mother, i, flag, judgeBy)
       endif

       if(flag == 0 .and. (sname == splitSpec(subdc) ) ) then
!              found a subd/media. fill the counter
          filled(subdc) = filled(subdc) + 1
          filledmax(subdc) = max( filled(subdc), filledmax(subdc))
          depthList(subdc) = depth
                  ! see if last node
          if( subdc < nitem ) then
!             see deeper subd.
             subdc = subdc + 1
          endif
       endif                          
       call epqSubdIdx(i, msubdi)
         ! msubdi is a  subdetector index to which         
         ! the component i belongs (if compn                
         ! is a simple comp.)                             
         ! If it is a subdetector, msubdi becomes            
         ! its subd index.                                
         ! if i is invalid, -1                        
         ! if i  does not belong to a subd, 0 
       do m = j, 1, -1
                 ! for partially contained one, CnArea.... < 0
          ii =abs( CnArea( Det%cmp(i)%ContainsR+m ) )
          call epCsubdTree0(msubdi, ii, subdc)
       enddo

       end subroutine epCsubdTree0

      

      end module modGetIndex

      subroutine epqSubdMediaName(mother, i, sname, subdc,icon)
      use modGetIndex
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer,intent(in) ::  mother
      integer,intent(in) ::  i
      character(len=*),intent(out):: sname  ! subd name or media name

      integer,intent(in):: subdc

      integer,intent(out) :: icon
!                    targetcpy    flag=0    flag = 1
!                                 sname     sname
!      spec(nitem)       simple      matter    no chg*
!                     container   no chg*   matter
!                     any         matter    matter
!                     subd        no chg*   no chg    
!         icon is made to be 0 except for
!         * case; in that case, icon = 1 results which means
!         this i does not match the requied condition
!         so need not compare the target and sname.
!         (If a medium name is used  as a subdetector name
!         (or vice versa), target  and sname could match;
!         so if icon is referred  such mistake can be avoided).
!   
!                           
      integer:: flag

      call epqSubdName2(mother, i, sname, flag)

      icon = 0 

      if(subdc /= nitem) return !!!!!!!!!!!

      if( flag == 0 ) then
         if( targetcpy == targetdef ) then
            sname = Det%cmp(i)%matter ! may be alias 
         elseif(targetcpy == 'container') then
            icon = 1
         elseif(targetcpy == 'any') then
            sname = Det%cmp(i)%matter ! may be alias
         elseif(targetcpy == 'subd') then
            icon = 1
         else
            write(0,*) ' target=',targetcpy, 
     *       ' invalid for epqSubdMediaName'
            stop
         endif
      else
         if( targetcpy == targetdef ) then
            icon = 1
         elseif(targetcpy == 'container') then
            sname = Det%cmp(i)%matter ! may be alias
         elseif(targetcpy == 'any') then
            sname = Det%cmp(i)%matter ! may be alias
         elseif(targetcpy == 'subd') then
         else
            write(0,*) ' target=',targetcpy,
     *       ' invalid for epqSubdMediaName'
            stop
         endif
      endif            
 
      end  subroutine epqSubdMediaName

