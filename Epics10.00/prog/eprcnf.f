#include  "ZepMaxdef.h"
      module modSubdWorld
      implicit none

!     %%%%%%%%
      integer,save:: bugbug=0
!     %%%%%%%

!        what follows is applied when a subdetector or 
!        main detecor is being proccessed. 
!        For each detector, counter etc are cleared
!                max # of subdetectors  that some comp.
!            can directly contains

      integer,parameter:: maxSubContainable = 100

      integer,parameter:: maxSubContainableNow = maxSubContainable
!              max # of comps that directly contains subdetecors
      integer,parameter:: maxSubContainers=40 
      integer,save::  maxSubContainersNow=maxSubContainers
!              such comp. counter
      integer,save:: containersCounter = 0  
!              counter for the # of subdetectors
!            that are contained in a given  container
      integer,save:: noOfContainedSub(maxSubContainers)
!                container's comp # List
      integer,save:: containersList(maxSubContainers)
!               list of contained sub detector comp. #
!                for each container
      integer,save::
     *  containedSubList(maxSubContainable, maxSubContainers)
!               material of each  container.  In some case
!             this matter will replace subdetercor's world matter 
      character(len=MAX_MEDIANAMELENG),save::
     *      containersMatter(maxSubContainers)


      end  module modSubdWorld

      module modNegativeC
      implicit none
      integer,save:: negativeABC =0 ! counter for boxes with negative egde length
      integer,save:: negativeCylH =0 ! counter for Cyl with negative height
      integer,save:: negativePipeH =0 ! counter for Pipe with negative height
      integer,save:: negativePrism =0 ! prism with negative param
      end  module modNegativeC


         ! if alias, replace it by true
       subroutine eprcnf(dsn)
!!!!        %%%%%%%%
       use modSubdWorld, only: bugbug
       use modNegativeC
!!!!        %%%%%%%%
       implicit none

#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
!
          character*(*) dsn
!
          character(MAX_STRUCCHR)::epparaphrase
          character(MAX_STRUCCHR)::tempph

          character(maxtotcha):: datax, dataout
          character*256 msg
          integer icon, i, j, nthline, klena
          integer incLevel
          character*120 incfile
          character*2 rep
          logical eattail
          character*1 HTab
!                       file # of include file starts from lvloffset+1
!                       max depth of inc is maxinc
          integer maxinc, lvloffset, maxparam
          parameter (maxinc=10, lvloffset=70, maxparam=9)
          integer ionum(maxinc), incparamN(maxinc)

          character*24 incparam(maxparam, maxinc)

          save eattail, incLevel, ionum
          integer::idpos
          logical,save:: parroting = .false. ! #parroting on/off will
                   ! make this t/f.  If t, every line is printed
          integer::uscl   ! under score location or max character pos
          character*2 bq/'\\'/  ! solaris cannot use '\'

!          next is for dpmjet file management. 
!          dpmfiles for Epcis is assumed to be in the
!          same directory as config file resides so
!          extract the directory where the config file is.
!          this must be called before cintModels is called.
          call cfixPrefix(dsn)
          call csetCosOrEpi('epics')

          eattail  =  .false.
          Incused = .false.
          nthline = 0
          incLevel = 0
          HTab = char(9)      
!           include file nestable upto maxinc; file # from 71 to ...
          do i = 1, maxinc
             ionum(i) = lvloffset + i
          enddo

          call copenf(iowk, dsn, icon)
          if(icon .ne. 0) then
             call cerrorMsg('file not exist', 1)
             call cerrorMsg(dsn, 0)
          endif
!              find separater
          call afsep(iowk)
          mode = 0
          Nequates = 0
          EqualF = -97531.d0
          EqualFshort = -6321
          PlusF =  -13579.d0
          MinusF = -3.1415d3
          OrgEqBy = 'o'     ! 'ebb'/EBB'/'BB';  'Oxyz' 'O'
          wrtcom = 1
!              clear counter
          call  epclCnfCnt
          call epiniparaph      ! clear paraphrase counter (only once)
          NsubD = 0
          cumsubdloc = 0
! ((((((((((((((((
!                these are cleared only once. 
!           
          AttrCounter = 0
          CnCounter = 0
          do i = 1, maxCnArea
             CnArea(i) = 0
          enddo
! ((((((((((((((((

          coment = 0
          icon=0
!          *** until loop*** 
          do while (icon .eq. 0)
             nthline = nthline + 1
             if(incLevel .eq. 0) then
!                      treatment of \ continuation.
!                      if \ is not in the '# type commen line', 
!                      it is regareded that the next line is the
!                      the continuation ; e.g
!                              0 0 1.0 \
!                              10  2  5
!                        ==>     0 0 1.0 10 2 5
!
                i = 1
                do while (i > 0) 
                   read(iowk, '(a)') datax(i:)
!                         <<<<  Mar.15.2014 for debug
                   if(i == 1) then
                      call epcheckDebugOnOff(datax, parroting)
                      if(parroting .or. MsgLevel >= 2) then
                         write(0,*) trim(datax(i:))
                      endif
                   endif
!                         >>>>
                   i = index(datax, bq(1:1))
                   if(i > 0 ) then   ! \  must not be at the col 1.
                      if( datax(1:2) == "# ") then
                         i = 0   ! if comment line, neglect "\"
                      else
                         i = i-1
                      endif
                   endif
                enddo
             else
                datax = ' '
                i = 1
                do while (i > 0) 
                   read(ionum(incLevel), '(a)', end=60) datax(i:)
!                      <<<<  Mar.15.2014 for debug                 
                   if(i == 1) then
                      call epcheckDebugOnOff(datax, parroting)
                      if(parroting .or. MsgLevel >= 2 ) then
                         write(0,*) trim(datax(i:))
                      endif
                   endif
!                        >>>>>>>> 
                   i = index(datax, bq(1:1))
                   if(i > 0 ) then   
                      if( datax(1:2) == "# ") then
                         i = 0   ! if comment line, neglect "\"
                      else
                         i = i-1
                      endif
                   endif
                enddo

                do i = 1, incparamN(incLevel)
                   rep = ' '
                   write(rep,'("$",i1)') i
                   call kgsub(rep, incparam(i, incLevel), datax,
     *               dataout)
                   datax = dataout
                enddo
             endif
!               decompse into  field(only first ~80 char) without
!               converting into lower case letter
             call kgetField(datax(1:), FieldAsItis, 50,  Nfield)
             call c2lowerCase(datax, confdata)
!               decompose into field(only first ~100 char), Nfield should
!               be the same as for FieldAsItis
             call kgetField(confdata(1:100), Field, 50, Nfield)
             if(Nfield .eq. 0) goto 30
             if(Field(1) .eq. '/*') then
                coment = coment  +1 
             elseif(Field(1) .eq. '*/') then
                if(coment .eq. 0) then
                   call cerrorMsg(
     *            'No /* corresponding  to */ in config', 1)
                   write(msg, *) ' The line number is ',nthline
                   call cerrorMsg(msg, 0)
                endif
                coment = coment - 1
                goto 30
             endif
             if(coment .gt. 0) then
                goto 30
             endif

             if(Nfield .ge. 2 .and. Field(1) .eq. '#subd')  then
                call eppDefSubd   ! #subd line
             elseif(Nfield .ge. 2 .and. Field(1) .eq. '#end') then
                call  eppEnd    !  #end line
             elseif(Nfield .ge. 3 .and.( Field(1) .eq. '#eq' .or.
     *        Field(1) .eq. '#feq' ) ) then
                call eppEquate
             elseif(Nfield .ge. 2 .and. Field(1) .eq. '#com') then
                read(Field(2), *) wrtcom
             elseif(Nfield == 2 .and. Field(1) == "#org=" ) then
                read(Field(2), '(a)') OrgEqBy
                call c2lowerCase(OrgEqBy, OrgEqBy)
                if( OrgEqBy /= "o" .and.  OrgEqBy /= "ebb" )  then
                   write(0,*) ' #org= ',OrgEqBy, ' invalid'
                   stop
                else
                   write(0,*)' now #org ',OrgEqBy, ' is active'
                endif
             elseif(Nfield .ge. 2 .and. Field(1) .eq. '#inc') then
                 incLevel = incLevel + 1
                 Incused = .true.
                 if(incLevel .gt. maxinc) then
                    call cerrorMsg(
     *              '#inc nest is too deep; must be <20 ', 0)
                 endif
                 call cformFullPath(FieldAsItis(2), incfile)
                 call copenf(ionum(incLevel), incfile, icon)
                 if(icon .ne. 0) then
                    call cerrorMsg('file not exist', 1)
                    call cerrorMsg(incfile, 0)
                 endif
!                 save parameter

                 if(  Nfield-2 .gt. maxparam) then
                    call cerrorMsg('too many parameters for #inc',0)
                 endif
                 incparamN(incLevel) = Nfield-2                    
                 do i = 1, incparamN(incLevel)
                    incparam(i, incLevel) = ' '
                    incparam(i, incLevel) = FieldAsItis(i+2)
                 enddo
                    
                 goto 90
             elseif(Nfield .ge. 3 .and. 
     *            Field(1) .eq.  '#news') then
                tempph = epparaphrase(Field(2))
                if(tempph .eq. Field(2) ) then
!                      new definition; depo
                   call epparadepo(Field(2), FieldAsItis(3))
                else
                   call cerrorMsg(Field(2), 1)
                   call cerrorMsg(' is defined more than once', 1)
                   write(msg, *) 'The line number is ',nthline
                   call cerrorMsg(msg, 1)
                endif
             elseif( Nfield == 3 .and. Field(1) == "#alias" ) then
                call epProcAlias(FieldAsItis)
             elseif(confdata(2:10) .eq. '---------') then
                call eppfsep   !     end of config  data
                icon = 1
             elseif(confdata(1:1) .eq. 'c'.or. confdata(1:1) .eq. '!'
     *               .or. confdata(1:1) .eq. '#' 
     *               .or. datax(1:5) .eq. '     '  ) then
                call epsvcom(datax)
             else
                if( Det%nct .ge. ncmax) then
                   write(msg, *) ' Too many components >',ncmax
                   call cerrorMsg(msg, 0)
                else
                   Det%nct = Det%nct + 1
                endif
!               &&&&&
!                Field(2) : keep capital lett.  
!                so we replace it by FieldAsItis
             
!///////////                Field(2) = FieldAsItis(2)   /////////////
                Det%cmp(Det%nct)%cn = Det%nct

                Det%cmp(Det%nct)%level = 0 ! (((((
!!!!/////            Det.cmp(Det.nct).subdidx = 0 !^^^^^^
                if( mode == 1 ) then
                   Det%cmp(Det%nct)%subdidx = NsubD+1
                else
                   Det%cmp(Det%nct)%subdidx = 0
                endif
!!!!////
                Det%cmp(Det%nct)%fsubdc = 0 !^^^^^^
!                  proces  for  one component def.
                do i = 1, 9
                   Det%cmp(Det%nct)%direc(i) = 0.
                enddo
                Det%cmp(Det%nct)%direc(1) = 1.d0
                Det%cmp(Det%nct)%direc(5) = 1.d0
                Det%cmp(Det%nct)%direc(9) = 1.d0
! (((((((((((((
!                do i = 1, maxattr
!                   Det.cmp(Det.nct).vol(i) = 0.
!                enddo
! ))))))))))))
                
                if(Field(2) == 'box' .or.
     *             Field(2) == 'box_w' ) then
                   call eprbox(Det%cmp(Det%nct))

                elseif(Field(2) == 'cyl' .or.
     *                 Field(2)(1:4) == 'cyl_') then
                   call eprcyl(Det%cmp(Det%nct))

                elseif(Field(2) == 'pipe' .or.
     *                 Field(2)(1:5) == 'pipe_') then
                   call eprpip(Det%cmp(Det%nct))

                elseif(Field(2) .eq. 'prism' .or.
     *                Field(2)(1:6) .eq. 'prism_') then 
                   call eprprs(Det%cmp(Det%nct))

                elseif(Field(2) == 'sphere' .or.
     *                 Field(2) == 'sphere_w' ) then
                   call eprsph(Det%cmp(Det%nct))

                elseif(index(confdata, 'new-') .ne. 0) then
                   tempph = epparaphrase(Field(2))
                   if(tempph .eq. Field(2) ) then
!                      this is new new-xx. new-xx is deposited
                      call epparadepo(Field(2), Field(2))
                   endif
                   call eprNew(Det%cmp(Det%nct),Field(2))
                else
                   call epseeUnderScore(Field(2), uscl)
                   tempph = epparaphrase(Field(2)(1:uscl))
                   if(tempph(1:4) .eq. 'new-')  then
                      call eprNew(Det%cmp(Det%nct), Field(2))
                   else
!                      see if the name is defined as a SubDetector 
                      do  j = 1, NsubD
                         if(SubDName(j) .eq. Field(2)) goto 20
                      enddo
!                           invalid comp. name
                      call cerrorMsg('component name  error', 1)
                      write(msg, *) 'It is in ',
     *                nthline, '-th line which is:'
                      call cerrorMsg(msg, 1)
                   !!   call cerrorMsg(confdata, 1) -->next v9.164
                      call cerrorMsg(datax, 1) ! original case
                      call cerrorMsg("Forgot #news...?",1)
                      if(index(confdata, HTab) .gt. 0) then
                         call cerrorMsg(
     *                   'If only blanks are seen above, "Tab"', 1)
                         call cerrorMsg(
     *                   'in the line may be the criminal',0)
                      endif
                      stop


!                    -----------------
 20                   continue

                      if(SubDUsed .lt. maxSubDRef) then
                         SubDUsed = SubDUsed + 1 ! SubDetctor counter
                      else
                         call cerrorMsg('too many ref. to sub det', 0)
                      endif
                      
                      call eprIncSubD( Det%cmp(Det%nct) )
                      
                      PosInDet(SubDUsed) = Det%nct
                      SubDNumb(SubDUsed) = j
                      if(SubDUsed .eq. 1) then
                         SumOffset(SubDUsed) =  SubD(j)%nct -1
                      else
                         SumOffset(SubDUsed) =  SumOffset(SubDUsed-1)
     *                        +  SubD(j)%nct -1
                      endif
                   endif
                endif
!                     contain-contained relation
                call epcontain(eattail)
!!                     v9.164.  all comp can be now world
!!                  so omit next if
!!!            if(index(Det.cmp(Det.nct).struc, '_w') .eq. 0) then
!                     judge rotation 
                call eprotation(Det%cmp(Det%nct))
!!!            endif
             endif
 30          continue
             goto 90
 60          continue
!               include end
             close(ionum(incLevel))
             incLevel = incLevel - 1 
!!!!!!!!!!! This is not enough to check cmp # seq. 
!             if(incLevel == 0 ) then
!                Incused = .false.
!             endif
!!!!!!!!!!!!!
                
 90          continue
          enddo
  100     continue
          
          close(iowk)
!              ouput expanded config/count each comp.
          call eprecount  ! many of business here has
                  ! been mvoed to  epcountStruc
!               reset world
          if(Det%nworld .gt. 0) then
             call epresetWorld
          else
             call epenvlpAll
          endif
!!           v9.14   change matter: 'sp2'->'sp'
!!                                  'world'-->"the one for world"
          call epresetMatter
!!!!
!               set form such as form='mix','box' etc
          call epsetform


          if(MsgLevel .ge. 1) then
             write(msg,*)
     *       ' configuration has been read from '
             call cerrorMsg(msg, 1)
             call cerrorMsg(dsn, 1)
             call epcountStruc
             call epwriteStruc
!             write(msg,*) ' :the # of components=',Det.nct
!             call cerrorMsg(msg, 1)
!             if(Det.nbox .gt. 0) then
!                write(msg,*) ' # of box=',Det.nbox
!                call cerrorMsg(msg, 1)
!             endif
!             if(Det.ncyl .gt. 0) then
!                write(msg,*) ' # of cylinder=',Det.ncyl
!                call cerrorMsg(msg, 1)
!             endif
!             if(Det.npip .gt. 0) then
!                write(msg,*) ' # of pipe=',Det.npip
!                call cerrorMsg(msg, 1)
!             endif
!             if(Det.nprs .gt. 0) then
!                write(msg,*) ' # of prism=',Det.nprs
!                call cerrorMsg(msg, 1)
!             endif
!             if(Det.nsph .gt. 0) then
!                write(msg, *) ' # of spheres=', Det.nsph
!                call cerrorMsg(msg, 1)
!             endif             
!/////////////

!             do i = 1, MaxNewStruc
!                if(Det.nnew(i) .gt. 0) then
!                   if(i .le. 9) then
!                      write(tempph , '("new-",i1)')  i
!                   else
!                      write(tempph , '("new-",i2)')  i
!                   endif
!                   write(msg, '(" # of ",a,"=", i5)')
!     *                  epparaphrase(tempph), Det.nnew(i)
!                   call cerrorMsg(msg, 1)
!                endif
!             enddo
             if( negativeABC > 0 ) then
                write(0,'(a,i5)')
     *          "# of boxes defined with negative edge length=", 
     *          negativeABC
             endif
             if( negativeCylH > 0 ) then
                write(0,'(a,i5)')
     *          "# of Cylinders defined with negative height=", 
     *          negativeCylH
             endif
             if( negativePipeH > 0 ) then
                write(0,'(a,i5)')
     *          "# of Pipe defined with negative height=",
     *           negativePipeH
             endif
             if( negativePrism > 0 ) then
                write(0,'(a,i5)') 
     *         "# of prism defined with negative attributes",
     *         negativePrism
             endif
             if( negativeABC > 0 .or.  negativeCylH > 0 .or.  
     *           negativePipeH > 0 .or. negativePrism >0 ) then
                write(0,*)
     *          '****** Be sure negative one is not mistake *****'
                write(0,*) '       see ****** Info above, if any '
                write(0,*)
     *          '************************************************'
             endif
             call epshowWorldInfo
             call epseeContainers ! may be comment out
          endif
!              see if given matter is ready for sampling
          call epgetMedia
          if(eattail) then
             msg='Self eating comp. exists; see *** line above'
             call cerrorMsg(msg, 0)
          endif
          if(MsgLevel .ge. 1) then
             call epwriteMedia
          endif
!             adjust max path by referring to the thickness of each
!             layer
          call epresetMaxPath
      end
      subroutine epshowWorldInfo
      implicit none
!  #include  "ZepMaxdef.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      real(8):: org(3), abc(3), vol(MAX_IATTR)
      integer:: na 
      if(Det%nworld .gt.  0)  then
         write(0, *) trim( Det%cmp(Det%nct)%struc ),
     *   " world defined"
         if(Det%cmp(Det%nct)%struc == "box_w") then
            call epqcnf(org, abc)
            write(0,'(a, 1p,3g15.6)')
     *      "world(bounding box) origin=", org(:)
            write(0,'(a, 1p, 3g15.6)')
     *       'box length=',abc(:)
         elseif(Det%cmp(Det%nct)%struc == "sphere_w") then
            call epqorg(Det%nct, org)
            call epqvolatr(Det%nct, na, vol)
            write(0,'(a, 1p, 3g15.6)')
     *      'world(sphere) origin=',org(:)
            write(0,*) 'sphere radius =', vol(1)
         elseif(Det%cmp(Det%nct)%struc == "cyl_w") then
            call epqorg(Det%nct, org)
            call epqvolatr(Det%nct, na, vol)
            write(0,'(a,1p, 3g15.6)')
     *       'world(cyl) origin=',org(:)
            write(0,*) 'radiu and height =', vol(1:2)
!               next   v9.164
         elseif( index( trim( Det%cmp(Det%nct)%struc), "_w") >0) then
            call epqorg(Det%nct, org)
            call epqvolatr(Det%nct, na, vol)
            write(0,'(a,a,a, 1p, 3g15.6)')
     *       'world by ', Det%cmp(Det%nct)%struc,' origin=',org(:)
            write(0,*) 'attributes=', vol(1:na)
         else
            write(0,*)
     *      Det%cmp(Det%nct)%struc, ' cannot be world'
            stop
         endif
      else
         call epqcnf(org, abc)
         write(0,*) 'no world is given'
         write(0,'(a, 1p, 3g15.6)')
     *      "but the bounding box origin=", org(:)
         write(0,'(a, 1p, 3g15.6)') 'box length=',abc(:)
      endif
      end  subroutine epshowWorldInfo

      subroutine epseeContainers
      implicit none
!   #include  "ZepMaxdef.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
      character(256)::msg
      integer:: i, j
      integer:: Npc
!////////////// 
      integer:: mother, mth
!////////////// 
      Npc = 0
      do i= 1, Det%nct
         if(Det%cmp(i)%NPContainer .gt. 0) then
            write(msg,*)
     *      'Comp. ', i, ' is partially contained by',
     *      ' comp. ',
     *      ( CnArea( j+Det%cmp(i)%PContained ),
     *      j=1, Det%cmp(i)%NPContainer )
            call cerrorMsg(msg, 1)
            Npc = Npc  + 1
         endif
      enddo
      write(msg,*) " # of partially contained comps.=", Npc
      call cerrorMsg(msg, 1)

      do i = 1, Det%nct-1
         if(Det%cmp(i)%NContainer == 1) then  ! if contained 1
!            write(msg, *)
!     *           '#r of components which contain comp ',i,
!     *           '=', Det.cmp(i).NContainer
!            call cerrorMsg(msg, 1)
         elseif(Det%nworld > 0 .and. Det%cmp(i)%NContainer == 0) then
!//////////////
            !  see if i is really not listed in matreska
!            write(0,*) '# ',i, ' really solitary ?'
            call epseeIfNotInMatreska(i, mother, mth)
            if( mother /= 0 ) then  !!!!
!               write(0, *) " it's mother is #", mother
!               write(0, *) mth,"-th baby" 
               Det%cmp(i)%NContainer = 1
            else
               write(msg,*)
     *           'There is no  comp. which contains comp ',i
               call cerrorMsg(msg, 1)
               write(0,*) ' Error ????'
            endif  !!!!
         endif
      enddo
      end


      subroutine epcheckDebugOnOff(line, onoff)
      implicit none
      character(len=*),intent(in):: line
      logical,intent(inout)::onoff

      integer::nf

      character(len=24)::field(2)
      call kgetField(line, field, 2,  nf)
      if( nf>0 ) then
         if( field(1) == "#parroting" ) then
            if( nf == 2 ) then
               if( field(2) == "on") then
                  onoff = .true.
               elseif( field(2) == "off") then
                  if(onoff) then
                     write(0,*) trim(line)
                  endif
                  onoff = .false.
               endif
            endif
         endif
      endif
      end  subroutine epcheckDebugOnOff

      subroutine epseeIfNotInMatreska(i, mother, mth)
      implicit none
!cc     #include  "ZepMaxdef.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
      integer,intent(in):: i !  component #
      integer,intent(out):: mother
      integer,intent(out):: mth

      integer::j, k, l
      
      do j = 1, Det%nct
         k = Det%cmp(j)%NMatreska
         do l = 1, k
            if( CnArea( Det%cmp(j)%ContainsR+l ) == i) then
               mother = j
               mth = l
               return
            endif
         enddo
      enddo
      mother = 0
      mth = 0
      end

!     ********************
      subroutine eppDefSubd
      use modSubdWorld
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!         process #subd  line


       character*120 msg 
       integer i
       type(ep3Vec)::   org, abc

       if(mode .ne. 0 ) then
          call cerrorMsg(confdata,  1)
          call cerrorMsg('#subd  is nested; not permited', 0)
       endif
       mode =  1
       do i = 1,  NsubD
          if(SubDName(i) .eq. Field(2) ) then
             msg = 'sub detector ' // Field(2) //
     *             ' doubly defined'
             call cerrorMsg(msg, 0)
          endif
       enddo
       do  i = 1, NPreDefName
          if(PreDefName(i) .eq. Field(2)) then
             msg = 'sub detector ' // Field(2) //
     *             ' has been used by Epics'
             call cerrorMsg(msg, 0)
          endif
       enddo

       if(NsubD .ge. maxSubD) then
          call cerrorMsg(
     *     ' MAX_SUBD in ZepMaxdef%h is too small',0)
       endif
       SubDName(NsubD +1) = Field(2)
       return
!      **************
       entry  eppEnd
!      **************
!         process #end line
       if(mode .ne. 1) then
          call cerrorMsg('#end without #subd', 0) 
       endif
       if(Field(2) .ne. SubDName(NsubD+1)) then
          call cerrorMsg(confdata,  1)
          msg = Field(2) // ' should be '  // 
     *       SubDName(NsubD+1)
          call cerrorMsg(msg, 0)
       endif
!        examine if a subdetector is contained by another  comp.
!        it must containe world at last.
       call epexamSubd
       if(Det%nworld .gt. 0) then
!           set flag so that world contains the
!           comp.
          call epworld
       else
!          msg = SubDName(NsubD+1) // ' subdetector has'
!     *      // ' no world component at the end'
!          call cerrorMsg(msg, 0)
       endif

!         +, = for count,de, maxpthl, attrib and  X,Y, Z pos. are processed
       call epAttrb  ! + for attrib is processed
       call epXpos
       call epYpos
       call epZpos


       if(SubDUsed .gt. 0) then
!          if the sub detector  includes other subdetctors,
!          we need expand them
          call epexpand(0)
       endif
       NsubD =  NsubD + 1
       if( NsubD .gt. maxSubD ) then
          call cerrorMsg('MAX_SUBD in ZepMaxdef%h is too small',0)
       endif
       SubD(NsubD)%loc =  cumsubdloc 
!          move the sub detector 

!         renumber
       do i = 1, Det%nct
          Det%cmp(i)%cn = i
       enddo

!         reset world, if any       
       if(Det%nworld .gt. 0) then
          call epresetWorld
!              get x,y,z thickness
          call epGetThick(Det%nct, 
     *    XYZthick(NsubD), ORGsubd(1,NsubD) )
       else
          call epenvlpAll
!          inquire the result
          call epqcnf(org, abc)
          XYZthick(NsubD) = abc
          ORGsubd(1,NsubD)  =org%x
          ORGsubd(2,NsubD)  =org%y
          ORGsubd(3,NsubD)  =org%z
       endif
!
       call epmvDet(Det,  SubdArea( SubD(NsubD)%loc+1 ), SubD(NsubD))

       mode = 0
       call epclCnfCnt          ! clear counter
       return
!      *************
       entry eppfsep
!      *************
       if(mode .ne.  0) then
          msg ='#subd for '// SubDName(NsubD+1)
     *         // ' is not completed'
          call  cerrorMsg(msg, 0)
       endif
!       examine subdetector is contained by anohter comp., if so
!       the last comp of the sub detector must be world.

       call epexamSubd

       call epadjPC
       if(Det%nworld .gt. 0) then
          call epworld
       endif
!         +, = count,DE,maxpathl, for X,Y,Z pos. are processed
       call  epAttrb
       call  epXpos
       call  epYpos
       call  epZpos
       if(SubDUsed .gt. 0) then
!          if the detector  includes other subdetctors,
!          we nee expand it
          call epexpand(1)
       endif
       end
!     ********************
      subroutine eppEquate
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!         process #eq or #feq line
       character*120 msg 
       character*16 name
       real*8 val
       integer ival
       integer i, icon, j
 
       if(Nequates .eq. maxEq) then
          call cerrorMsg('no more #eq can be used', 1)
          write(msg,*)
     *    ' Expand MAX_EQUATE=',maxEq, ' in ZepMaxdef%h'
          call cerrorMsg(msg, 0)
       endif
!         check double definition
       do i = 1, Nequates
          if(Field(2) .eq. EqName(i)) then
             if(Field(1) .eq. '#eq') then
                write(msg, *) Field(2), ' is doubly defined'
                call cerrorMsg(msg, 0)
             else
!                 should be #feq so overwrite previous value
                j = i
                goto 10
             endif
          endif
       enddo
       Nequates = Nequates + 1
       EqName(Nequates) = Field(2)
       j = Nequates 
 10    continue
       read(Field(3), *) EqValue(j)
       return
!     ***********
       entry epexamEq(name, val, icon)
!         name: input. if name if found in EqValue, icon =0 
!                      and val becomes correspoiding EqValue
!         val:  output.
!         icon: output
       do i = 1, Nequates
          if(EqName(i) .eq. name) then
             val = EqValue(i)
             icon = 0
             return     ! ******
          endif
       enddo
       icon = 1
       return
!     ***********
       entry epexamEqI(name, ival, icon)
!         name: input. if name is found in EqValue, icon =0 
!                      and val becomes correspoiding EqValue
!         ival:  output.
!         icon: output
       do i = 1, Nequates
          if(EqName(i) .eq. name) then
             ival = EqValue(i) + 0.00001
             icon = 0
             return     ! ******
          endif
       enddo
       icon = 1
       end
!      ************
       subroutine eperrEq(comp, fieldin)
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
       type(Component)::  comp 
       character*(*) fieldin

       character*128 msg
       write(msg, *) comp%cn,'-th component has undefined variable'
       call cerrorMsg(msg, 1)
       write(msg, *) 'which is: ', fieldin
       call cerrorMsg(msg, 0)
       end
!     ****************
      subroutine epexamSubd
       use  modSubdWorld
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!           see if a sub detector is contained by another
!         component. if so , the last component of the
!         sub detector must be a world. (may be implicit)

      character*150 msg

      integer i, j, epIsSubD, k, m
      integer jj, kk ! v9.14
      integer:: locus  ! v9.164
!!       v9.14
      ! clear counter for world business
      containersCounter = 0
      noOfContainedSub(:) = 0
      containersList(:) = 0
      containedSubList(:,:) = 0
      containersMatter(:) =" "

      do i = 1,  Det%nct

         do j = 1, Det%cmp(i)%NMatreska
!!!!    V9.164    check specified cmp # to be contained  >= Det.nct
            if(  CnArea( Det%cmp(i)%Contains+j)  >= Det%nct ) then
               if( SubDName(NsubD+1) /= " ") then
                  write(0,*) 'Error: in subdetector ', SubDName(NsubD+1)
               endif
               write(0,'(a, i0,a, i0, a,i0,a,i0)')
     *          ' "|cn #| of ', j,'-th daughter of component # ', i,
     *          '" is ',
     *          CnArea( Det%cmp(i)%Contains+j), ' and is >= ', Det%nct
               stop
            endif
!(((((((((((((((((
!           k = epIsSubD(  Det.cmp(i).Contains(j) )
            k = epIsSubD( CnArea( Det%cmp(i)%Contains+j) )
!)))))))))))))))))

            if( k .gt. 0) then
!               a subdetector k is contained by i-th detector
!               directly, so that k must have world at its end
               m = SubD(k)%nct

               if( index(SubdArea( SubD(k)%loc+m)%struc, '_w')
     *               .ne. 0 ) then
!                   world exists
!!!!!!!!!!!!!!!!!!! v9.14
                  kk = 0
                   ! examin each container  of a subd so far
                  do jj = 1, containersCounter
                      ! see if the container's comp # is i
                     if( containersList(jj) == i ) then
                        kk = jj   ! kk-th container 
                        exit
                     endif
                  enddo
                  if( kk == 0 ) then
                     if( containersCounter ==
     *                     maxSubContainersNow) then
                        write(0,*)
     *                  'too many containers of subdetectors'
                        write(0,*)
     *                   ' change maxSubContainers in eprcnf%f '
                        write(0,*) 
     *                   ' or implement dynamic memory alloc. '
                        stop
                     endif
                       ! update # of containers 
                     containersCounter = containersCounter + 1
                     kk =  containersCounter  ! kk-th container
                     containersList(kk) = i  ! is i
                            ! it's matter is memorized
                     containersMatter(kk) = Det%cmp(i)%matter
                  endif
               !   see if sub contained counter is full
                  if( noOfContainedSub(kk) == 
     *                  maxSubContainableNow ) then
                     write(0,*)
     *                'too many contained subdetectors>',
     *                 maxSubContainable
                     write(0,*) 
     *                'consider making subdeterctor(s) to ',
     *                'contain them  or '
                     write(0,*)
     *                ' change maxSubContainable in eprcnf%f '
                     stop
                  endif
                    !  update kk-th contained subd counter
                  noOfContainedSub(kk) = noOfContainedSub(kk) + 1
                    ! memorize j-th matreska's comp #
                  containedSubList(noOfContainedSub(kk), kk) =
     *                     CnArea( Det%cmp(i)%Contains+j )
!!!!!!!!!!!!!!!
               else
                 if(m .gt. 1 .and. 
!     *               SubD(k).cmp(m).NMatreska .eq. 0) then
     *               SubdArea( SubD(k)%loc+m )%NMatreska .eq. 0) then
                    write(msg, *)' Subdetector=', SubDName(k),
     *                ' has no implicit world in the last comp.'
                    call cerrorMsg(msg, 1)
                    call cerrorMsg(
     *              'Sub-detetctor to be contained by another'//
     *              ' comp. must have (implicit) world in the'//
     *              ' last component', 0)
                 elseif( m .gt. 1) then
!                      issue warning message
                    call cerrorMsg('********* Info *******',1)
                    write(msg, *)
     *              ' Now cheking  that  Subdetector=', SubDName(k),
     *              ' has an implict world in the last comp.'
                    call cerrorMsg(msg, 1)
!                       next v9.14
                    call episeqvWorld(
     *                k, SubD(k), SubdArea( SubD(k)%loc+1 ) )

                 endif
              endif
           endif   
        enddo
      enddo
      end
!      ***********************
       subroutine epresetMaxPath
       use epModify
       implicit none
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"
#include "Zmass.h"



       type(ep3Vec)::  vec

       integer i, mediaindex
       real(8):: mingram 
       integer modi
       logical usemodifier
       real(8):: dummy(3)

!         max path length allowable in r.l
       do i = 1, Det%nct
          mediaindex = Det%Cn2media(i)
          if(Det%cmp(i)%MaxPathL .eq.  0.) then
!               we set a negative value to indicate it is by system
             if(Det%cmp(i)%matter .eq. 'sp' .or. 
     *          Det%cmp(i)%matter .eq. 'sp2' .or. 
     *          Det%cmp(i)%matter .eq. 'hollow' ) then
                Det%cmp(i)%MaxPathL = -10000.  
             else
                call epGetThick(i, vec, dummy)
                if( min(vec%z,vec%x, vec%y)/Media(mediaindex)%X0 *
     *               Media(mediaindex)%rhoc  .lt. 5.d-2) then
!                        for very thin layer, we assume  such as
!                        emulsion. no effective scattering. 
!                        unless maxpath is specified.                  
                   Det%cmp(i)%MaxPathL= -1.e-2*Media(mediaindex)%X0 /
     *                 Media(mediaindex)%rhoc
                else
                  Det%cmp(i)%MaxPathL = 
     *             - max(0.05d0 * min(vec%x, vec%y, vec%z),
     *                    0.02d0*Media(mediaindex)%X0/
     *                 Media(mediaindex)%rhoc )
                endif
             endif
          else
!            use given value (cm) as it is
!
          endif
!             Emin business:  see if modification requested for Emin
          modi = Det%cmp(i)%modifier 
!             see if ModifyFile existed
          usemodifier = allocated(modify) .and. modi > 0
          if( usemodifier ) then
!               see if emin data was given
             if( modi <= maxModifyNum ) then
                usemodifier = IBITS(modify( modi )%kind, bitEmin, 1) > 0
             else
                write(0,*) ' modifier =', modi, ' in component #=',i
                write(0,*) ' > that (=',maxModifyNum,
     *             ') given in ModifyFile '
                write(0,*) ' To run the program '
                write(0,*) ' 1) correct data in modifier or ModifyFile'
                write(0,*) ' or '
                write(0,*) ' 2) give " " as ModifyFile to neglect all ',
     *                     ' modifiers'
                stop
             endif
          endif

          if( usemodifier ) then 
             Det%cmp(i)%EminG = modify(modi)%Em%Egmin
             Det%cmp(i)%EminE = modify(modi)%Em%Eemin
             Det%cmp(i)%RecoilE =modify(modi)%Em%Recoilmin
          else
             if(AutoEmin >= 1) then
!                          AutoEmin 
!                        fix Emin by looking at min size in g/cm2
!     
                mingram =  min(vec%x, vec%y,vec%z)*
     *               Media(mediaindex)%rho  *Det%cmp(i)%rhoc
                call epAutoEmin( mingram, Media(mediaindex), 
     *          Det%cmp(i)%EminG, Det%cmp(i)%EminE, Det%cmp(i)%RecoilE)
!/////////////////
!                write(0,*) ' mingram=', mingram
!                write(0,*) ' vec ',vec.x, vec.y, vec.z
!                write(0,*) ' min =',  
!     *          Det.cmp(i).EminG, Det.cmp(i).EminE, Det.cmp(i).RecoilE
!////////////////
             else
                Det%cmp(i)%EminG = EminGsave
                Det%cmp(i)%EminE = EminEsave
                Det%cmp(i)%RecoilE = RecoilEsave
!////////////////
!                write(0,*) 'Auto=',AutoEmin
!                write(0,*) 'min=',
!     *          Det.cmp(i).EminG, Det.cmp(i).EminE, Det.cmp(i).RecoilE
!/////////////////
             endif
          endif   
       enddo
       end

!      ***********************
       subroutine epresetWorld
       implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"


       type(epPos)::  org
       type(ep3Vec)::  abc
       integer i  ! ((()))
!
!          if Det( or Det to be a sub detector)  includes world
!        compute the actual world size here, and reset its size
!
!          cmopute envloping box
       call epenvlpAll
!          inquire the result
       call epqcnf(org, abc)
!
       if(Det%cmp(Det%nct)%struc .eq. 'box_w') then
          Det%cmp(Det%nct)%orgx = org%x
          Det%cmp(Det%nct)%orgy = org%y
          Det%cmp(Det%nct)%orgz = org%z
!(((((((((((
!          Det.cmp(Det.nct).vol(boxa) = abc.x
!          Det.cmp(Det.nct).vol(boxb) = abc.y
!          Det.cmp(Det.nct).vol(boxc) = abc.z
          Volat(Det%cmp(Det%nct)%vol+boxa) = abc%x
          Volat(Det%cmp(Det%nct)%vol+boxb) = abc%y
          Volat(Det%cmp(Det%nct)%vol+boxc) = abc%z
!)))))))))))

       elseif( Det%cmp(Det%nct)%struc .eq. 'sphere_w' ) then
          if(Volat( Det%cmp(Det%nct)%vol + 1) .eq. 0.) then
!                 world size not given. so set it here
             Det%cmp(Det%nct)%orgx = org%x + abc%x/2
             Det%cmp(Det%nct)%orgy = org%y + abc%y/2
             Det%cmp(Det%nct)%orgz = org%z + abc%z/2
!(((((((((((((9
!          Det.cmp(Det.nct).vol(sphr) =
!     *    sqrt( abc.x**2 +  abc.y**2 +  abc.z**2 )/2.d0 
             Volat( Det%cmp(Det%nct)%vol+sphr ) =
     *            sqrt( abc%x**2 +  abc%y**2 +  abc%z**2 )/2.d0 
          else
!           do nothing . respect original attribute
          endif
!))))))))))))
       elseif( Det%cmp(Det%nct)%struc(1:4) .eq. 'cyl_' .and.
     *        index(Det%cmp(Det%nct)%struc, '_w') > 0 ) then

          if( Volat( Det%cmp(Det%nct)%vol + 1) .eq. 0.) then

             if( Det%cmp(Det%nct)%struc(1:5) == 'cyl_w' .or.
     *           Det%cmp(Det%nct)%struc(1:7) == 'cyl_z_w' ) then
                Det%cmp(Det%nct)%orgx = org%x + abc%x/2
                Det%cmp(Det%nct)%orgy = org%y + abc%y/2
                Det%cmp(Det%nct)%orgz = org%z 
                Volat( Det%cmp(Det%nct)%vol+cylr ) =
     *            sqrt( abc%x**2 +  abc%y**2 )/2.d0 
                Volat( Det%cmp(Det%nct)%vol+cylh ) = abc%z
             elseif(Det%cmp(Det%nct)%struc(1:7) == 'cyl_y_w') then
                Det%cmp(Det%nct)%orgx = org%x + abc%x/2
                Det%cmp(Det%nct)%orgy = org%y 
                Det%cmp(Det%nct)%orgz = org%z + abc%z/2
                Volat( Det%cmp(Det%nct)%vol+cylr ) =
     *            sqrt( abc%x**2 +  abc%z**2 )/2.d0 
                Volat( Det%cmp(Det%nct)%vol+cylh ) = abc%y
             elseif(Det%cmp(Det%nct)%struc == 'cyl_x_w') then
                Det%cmp(Det%nct)%orgx = org%x 
                Det%cmp(Det%nct)%orgy = org%y + abc%y/2
                Det%cmp(Det%nct)%orgz = org%z + abc%z/2
                Volat( Det%cmp(Det%nct)%vol+cylr ) =
     *            sqrt( abc%y**2 +  abc%z**2 )/2.d0 
                Volat( Det%cmp(Det%nct)%vol+cylh ) = abc%x
             else
                write(0,*) 'cylinder: '
                write(0,*) Det%cmp(Det%nct)%struc
                write(0,*) ' is invalid '
                stop
             endif
!(((((((((((((
!          Det.cmp(Det.nct).vol(cylr) =
!     *    sqrt( abc.x**2 +  abc.y**2 )/2.d0 
!          Det.cmp(Det.nct).vol(cylh) = abc.z
          endif
!))))))))))))))
!     !       else   V9.164 replaced by next 6 lines
       elseif(index(trim(Det%cmp(Det%nct)%struc),'_w') /= 0 ) then
          if( Volat( Det%cmp(Det%nct)%vol + 1) .eq. 0.) then
             call cerrorMsg( Det%cmp(Det%nct)%struc, 1)
             call cerrorMsg(' is used as world:', 1)
             call cerrorMsg(' it must have volume attrib',0)
          endif
!          call cerrorMsg(Det.cmp(Det.nct).struc, 1)
!          call cerrorMsg(' is not suppored for world yet', 0)
       endif
!((((((((  1 level higher for comp. in  a world.
       do i =1 , Det%nct - 1
          Det%cmp(i)%level = Det%cmp(i)%level +1
       enddo
!))))))))
       end
!     *****************
      subroutine epXpos
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
!          X position is embeded. If = or + is
!        included, they are replaced by an actual value
       integer i, j
       real(8):: Borg(3), Babc(3)
       real(8):: Borgp(3), Babcp(3)


       real(8),external::epwBBmax, epwBBmin, eplBBmin
       integer,external:: epIsSubD


       do i =1, (Det%nct-Det%nworld)
          if(Det%cmp(i)%orgx .eq. EqualF) then
             if(i .eq. 1) then
                call cerrorMsg(
     *         '1st comp. must not have = in x-pos.',0)
             else
                if(epIsSubD(i-1) > 0 .or.  epIsSubD(i) > 0 .or. 
     *             OrgEqBy == "ebb" ) then
!                ELBBmin is mapped onto EWBBmin and offset is
!                added. (LEBBmin  Effective local B.B min;
!                WEBBmin  Effective world B.B min)
                   Det%cmp(i)%orgx = 
     *                  epwBBmin(i-1, 1) - eplBBmin(i, 1) 
     *                  + Det%cmp(i)%offsetx
                else  ! use simple mapping (old def).
                   Det%cmp(i)%orgx = Det%cmp(i-1)%orgx
     *                      + Det%cmp(i)%offsetx
                endif
             endif
          elseif(Det%cmp(i)%orgx .eq. PlusF) then
             if(i .eq. 1) then
                call cerrorMsg(
     *           '1st comp. must not have +',0)
             else
!                   ELBBmin is mapped onto EWBBmax and offset is
!                   added
                Det%cmp(i)%orgx = 
     *               epwBBmax(i-1, 1) - eplBBmin(i, 1) 
     *                      + Det%cmp(i)%offsetx
             endif
          elseif(Det%cmp(i)%orgx .eq. MinusF) then
             if(i .eq. 1) then
                call cerrorMsg(
     *           '1st comp. for origx must not have -',0)
             else
                call epqenvlper0(i, Borg, Babc)
                !  only for simple comp.
                call epqenvlper(i-1, Borgp, Babcp)
                Det%cmp(i)%orgx =Borgp(1) + Babcp(1) +
     *                             Borg(1)
             endif
          else ! simple numerical value is in orgx 
               ! mapping is automatic
          endif
       enddo
       end
!     *****************
      subroutine epYpos
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
!          Y position is embeded. If = or + is
!        included, they are replaced by an actual value
       integer i, j
       real(8):: Borg(3), Babc(3)
       real(8):: Borgp(3), Babcp(3)



       real(8),external::epwBBmax, epwBBmin, eplBBmin
       integer,external:: epIsSubD


       do i =1, (Det%nct-Det%nworld)
          if(Det%cmp(i)%orgy .eq. EqualF) then
             if(i .eq. 1) then
                call cerrorMsg(
     *         '1st comp. for origy must not have =',0)
             else
                if(epIsSubD(i-1) > 0 .or.  epIsSubD(i) > 0 .or. 
     *             OrgEqBy == "ebb" ) then
!                ELBBmin is mapped onto EWBBmin and offset is
!                added. (LEBBmin  Effective local B.B min;
!                 WEBBmin  Effective world B.B min)
                   Det%cmp(i)%orgy = 
     *                  epwBBmin(i-1, 2) - eplBBmin(i, 2) 
     *                  + Det%cmp(i)%offsety
                else  ! use simple mapping (old def).
                   Det%cmp(i)%orgy = Det%cmp(i-1)%orgy
     *                      + Det%cmp(i)%offsety
                endif
             endif
          elseif(Det%cmp(i)%orgy .eq. PlusF) then
             if(i .eq. 1) then
                call cerrorMsg(
     *           '1st comp. must not have +; origY ',0)
             else
                Det%cmp(i)%orgy = 
     *               epwBBmax(i-1, 2) - eplBBmin(i, 2) 
     *                      + Det%cmp(i)%offsety
             endif
          elseif(Det%cmp(i)%orgy .eq. MinusF) then
             if(i .eq. 1) then
                call cerrorMsg(
     *           '1st comp. for origy must not have -',0)
             else
                call epqenvlper0(i, Borg, Babc)
                !  only for simple comp.
                call epqenvlper(i-1, Borgp, Babcp)
                Det%cmp(i)%orgy =Borgp(2) + Babcp(2) +
     *                             Borg(2)
             endif
          else ! simple numerical value is in orgy
          endif
       enddo
       end
!     *****************
      subroutine epZpos
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
!          Z position is embeded. If = or + is
!        included, they are replaced by an actual value
       integer i, j
       real(8):: Borg(3), Babc(3)
       real(8):: Borgp(3), Babcp(3)


       real(8),external:: epwBBmax,epwBBmin, eplBBmin
       integer,external:: epIsSubD


       do i =1, (Det%nct-Det%nworld)
          if(Det%cmp(i)%orgz .eq. EqualF) then
             if(i .eq. 1) then
                call cerrorMsg(
     *         '1st comp. must not have = in z-pos.',0)
             else
                if(epIsSubD(i-1) > 0 .or.  epIsSubD(i) > 0 .or. 
     *             OrgEqBy == "ebb" ) then
                   Det%cmp(i)%orgz = 
     *                  epwBBmin(i-1, 3) - eplBBmin(i, 3) 
     *                  + Det%cmp(i)%offsetz
                else  ! use simple mapping (old def).
                   Det%cmp(i)%orgz = Det%cmp(i-1)%orgz
     *                      + Det%cmp(i)%offsetz
                endif
             endif
          elseif(Det%cmp(i)%orgz .eq. PlusF) then
             if(i .eq. 1) then
                call cerrorMsg(
     *           '1st comp. for origz must not have +',0)
             else
                Det%cmp(i)%orgz = 
     *               epwBBmax(i-1, 3) - eplBBmin(i, 3) 
     *                      + Det%cmp(i)%offsetz
             endif
          elseif(Det%cmp(i)%orgz .eq. MinusF) then
             if(i .eq. 1) then
                call cerrorMsg(
     *           '1st comp. for origz must not have -',0)
             else
                call epqenvlper0(i, Borg, Babc)
                !  only for simple comp.
                call epqenvlper(i-1, Borgp, Babcp)
                Det%cmp(i)%orgz =Borgp(3) + Babcp(3) +
     *                             Borg(3)
             endif
          else ! simple numerical value is in orgz 
          endif
       enddo
       end
      function epwBBmax(i, axis) result(ans)
      implicit none
      ! get effective world BB max coordinate value 
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer,intent(in):: i ! i-th componet
      integer,intent(in):: axis ! one of 1,2,3 for x, y, z
                     ! axis
      real(8):: ans
      real(8):: Borg(3), Babc(3)
      integer::j 
      integer,external:: epIsSubD

      j = epIsSubD(i) 
      if( j > 0 ) then
         if( axis == 1 ) then
            ans = Det%cmp(i)%orgx +  XYZthick(j)%x 
     *      + ORGsubd(1, j)
         elseif( axis == 2 ) then
            ans = Det%cmp(i)%orgy +  XYZthick(j)%y 
     *      + ORGsubd(2, j)
         elseif(  axis == 3 ) then
            ans = Det%cmp(i)%orgz +  XYZthick(j)%z 
     *      + ORGsubd(3, j)
         else
            write(0,*)
     *      'axis =',axis, ' invalid : epwBBmax 1'
            stop
         endif
      else
         call epqenvlper(i, Borg, Babc)
         if(axis >= 1  .and. axis <=3 ) then
            ans = Borg(axis) + Babc(axis)
         else
            write(0,*)
     *      'axis =',axis, ' invalid : epwBBmax 2'
            stop
         endif
      endif
      end    function epwBBmax
      function epwBBmin(i, axis) result(ans)
      implicit none
      ! get effective world BB min coordinate value 
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer,intent(in):: i ! i-th componet
      integer,intent(in):: axis ! one of 1,2,3 for x, y, z
                     ! axis
      real(8):: ans

      real(8):: Borg(3), Babc(3)
      integer::j 
      integer,external:: epIsSubD

      j = epIsSubD(i) 
      if( j > 0 ) then
         if( axis == 1 ) then
            ans = Det%cmp(i)%orgx + min( ORGsubd(1,j), 0.d0)
         elseif( axis == 2 ) then
            ans = Det%cmp(i)%orgy + min( ORGsubd(2,j), 0.d0)
         elseif(  axis == 3 ) then
            ans = Det%cmp(i)%orgz + min( ORGsubd(3,j), 0.d0)
         else
            write(0,*)
     *      'axis =',axis, ' invalid : epwBBmin 1'
            stop
         endif
      else
         call epqenvlper(i, Borg, Babc)
         if(axis >= 1  .and. axis <=3 ) then
            ans = Borg(axis) 
         else
            write(0,*)
     *      'axis =',axis, ' invalid : epwBBmin 2'
            stop
         endif
      endif
      end    function epwBBmin

      function eplBBmin(i, axis) result(ans)
      implicit none
      ! get effective local BB min coordinate value 
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer,intent(in):: i ! i-th componet
      integer,intent(in):: axis ! one of 1,2,3 for x, y, z
                     ! axis
      real(8):: ans


      real(8):: Borg(3), Babc(3)
      integer:: j
      integer,external:: epIsSubD
      if(i <  1) then
         write(0,*) ' i=',i, ' invalid for eplBBmin'
         stop
      endif

      j = epIsSubD(i) 
      if( j > 0 ) then
         if( axis >= 1 .and. axis <= 3 ) then
            ans = min(ORGsubd(axis,j), 0.d0)
         else
            write(0,*)
     *      'axis =',axis, ' invalid : epwBBmin 2'
            stop
         endif
      else
         call epqenvlper0(i, Borg, Babc)
         ans = Borg(axis)
      endif
      end   function eplBBmin

      function eporigDelta(i, axis) result(ans)
!          when = is specified, effective 
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer,intent(in):: i ! comp. #
      integer,intent(in)::axis !  1 2 or 3 for x, y, 
      real(8):: ans   ! for

      integer::j
      integer,external:: epIsSubD

      j = epIsSubD(i) 
      if( j > 0 ) then
         if( axis >= 1 .and. axis<= 3 ) then
            if( ORGsubd(axis,j) > 0.d0 ) then
               ans = 0.
            else
               ans = 0.
            endif
         else
            write(0,*) ' axis =',axis, ' invalid'
            write(0,*) ' for eporigDelta'
            stop
         endif
      else
         ans = 0.
      endif
      end    function eporigDelta
!     *****************
      subroutine epAttrb
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!       If = is in   attribute, replaced by actual value 
!    This has been done in epAttrb1 already. (dangerous though)
!     see Note above "call  epAttrb1"
!
       integer i, j, epIsSubD
       character*120 msg
       
       do i =1, (Det%nct-Det%nworld)
          if(Det%cmp(i)%CountIO .eq. EqualFshort) then
             if(i .eq. 1) then
                call cerrorMsg(
     *          '1st comp. must not have = in C field',0)
             endif
             Det%cmp(i)%CountIO = Det%cmp(i-1)%CountIO
          endif
          if(Det%cmp(i)%CountDE .eq. EqualFshort) then
             if(i .eq. 1) then
                call cerrorMsg(
     *          '1st comp. must not have = in DE field',0)
             endif
             Det%cmp(i)%CountDE = Det%cmp(i-1)%CountDE
          endif
          if(Det%cmp(i)%MaxPathL .eq. EqualF) then
             if(i .eq. 1) then
                call cerrorMsg(
     *          '1st comp. must not have = in MaxPathL field',0)
             endif
             Det%cmp(i)%MaxPathL = Det%cmp(i-1)%MaxPathL
          endif
          if(Det%cmp(i)%modifier .eq. EqualFshort ) then
             if(i .eq. 1) then
                call cerrorMsg(
     *          '1st comp. must not have = in modifere field',0)
             endif
             Det%cmp(i)%modifier = Det%cmp(i-1)%modifier
          endif
          if(Det%cmp(i)%chno .eq. EqualFshort ) then
             if(i .eq. 1) then
                call cerrorMsg(
     *          '1st comp. must not have = in chno field',0)
             endif
             Det%cmp(i)%chno = Det%cmp(i-1)%chno
          endif
!     
          do j = 1, Det%cmp(i)%Nattributes
!(((((((((((
             if(Volat( Det%cmp(i)%vol+j) .eq. EqualF) then
!)))))))))))
                if(i .eq. 1) then
                   call cerrorMsg(
     *            '1st comp. must not have = in vol. attrb.',0)
                elseif(epIsSubD(i-1) .gt. 0) then
!                     previous one is sub detector; avoid use of = 
                   write(msg, *)
     *             ' previous comp of ', i, '-th comp. is sub-D' // 
     *             ' avoid use of = in such case '                    
                   call cerrorMsg(msg, 0)
                else            ! previous is ordinary detector
!((((((((((
                  Volat( Det%cmp(i)%vol+j) =Volat(Det%cmp(i-1)%vol+j)
     *                 + VolatEq(Det%cmp(i)%vol+j)                 
!))))))))))
                endif
             endif
          enddo
       enddo
       end
!     *****************
      subroutine epAttrb1
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!       If = is in   attribut, replaced by actual value 
!    see Note above the "call epAttrb1"

       integer i, j, epIsSubD, k
       character*120 msg

      integer loc(maxiattr)
      logical got

      got = .false.

      do i = Det%nct, Det%nct-Det%nworld
         do j = 1, Det%cmp(i)%Nattributes
!(((((((((((
            if(Volat( Det%cmp(i)%vol+j) .eq. EqualF) then
!)))))))))))
               if(i .eq. 1) then
                  call cerrorMsg(
     *            '1st comp. must not have = in vol. attrb.',0)
               elseif(epIsSubD(i-1) .gt. 0) then
!                     previous one is sub detector; avoid use of = 
                  write(msg, *)
     *             ' previous comp of ', i, '-th comp. is sub-D' // 
     *             ' avoid use of = in such case '                    
                  call cerrorMsg(msg, 0)
               else             ! previous is ordinary detector
                  if(Det%cmp(i-1)%struc .ne. Det%cmp(i)%struc) then
!                     structure is diff. from previous. cannot use
!                     "=" notation.
                     write(msg,*)
     *               ' structure=', Det%cmp(i)%struc, 
     *               ' is diff. from previous one: ',
     *               Det%cmp(i-1)%struc
                     call cerrorMsg(msg, 1)
                     call cerrorMsg(
     *               'so = notation for the attrib cannot be used',0)
                  endif
!((((((((((
                  if(.not. got) then
!                         correspondence 
                     call epatloc(Det%cmp(Det%nct-1), loc)       
                     got = .true.
                  endif
                  k = loc(j)  ! jth original value is at k
                  Volat( Det%cmp(i)%vol+j) =Volat(Det%cmp(i-1)%vol+k)
     *                 + VolatEq(Det%cmp(i)%vol+j)                 
!))))))))))
               endif
            endif
         enddo
       enddo
       end

      integer function epIsSubD(i)
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!       Is the i-th component a sub detector ?
!       if yes, function gets the sub detector number
!       else gets 0

       integer i,  j

!       do j = 1, NsubD
       do j = 1, SubDUsed
          if(PosInDet(j) .eq. i) then
             epIsSubD = SubDNumb(j)
             goto 10
          endif
       enddo
       epIsSubD = 0
 10    continue
       end

      subroutine eprIncSubD(comp)
       implicit none
#include  "Zep3Vec.h"
#include "Zcnfig.h"
       type(Component)::  comp 
!            read a line including a sub detector 

!              read part of the structure
       call eprpst(comp,  0, 0, 1,  6)

       end
!      *****************
       subroutine ep2actn(n, actn, subDx)
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!         given component number n is
!         converted into actual component number
!       ( this  is needed when sub detectors are
!         included )
       integer n     ! input.  component number in Det.
       integer actn  ! output. component number when
                     !         referred sub detector is 
                     !         expanded.
       integer subDx ! output. sub detector number  if n is one 
                     !         of  PosInDet(1)
                     !         ~ PosIndet(SubDUsed)
                     !         else  0
       integer j
       

       do j = SubDUsed, 1, -1
!             see if n is >= j-th referred SubDector
!             pos in Det config.

          if(n  .gt. PosInDet(j)) then
             subDx = 0
             actn = n + SumOffset(j)
             goto 10
          elseif(n .eq. PosInDet(j)) then
             subDx = SubDNumb(j)
             actn = n + SumOffset(j)
             goto 10
          endif
       enddo
       subDx = 0
       actn = n
 10    continue
       end
!         move included sub Det into Vcomp as a compoent.
       subroutine  epexpand(final)
       use  modSubdWorld
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"

       integer final  ! input.  0--> within subdetector
                      !         1--> within detector
!            expand included sub detectors
!     Supose the current Det config. as follows
!
!     Det     Contains/PContains
!      1     /6
!      2     /3
!      3     /5 11
!      4     /10
!      5 sub1           sub1  1 ~ n1
!      6     
!      7     /12
!      8     /9
!      9
!     10  sub2          sub2  1 ~ n2 
!     11
!     ..
!    This should be expanded as follows
!
!    1       /6 + n1 -1 
!    2       /3
!    3       /5 + n1 -1, 11 + n1-1 + n2 -1
!    4       /10+ n1-1 + n2 -1
!    5      expanded sub1
!    6
!   ..
!    5+n1-1  last of sub1
!    6+n1-1  
!    7+n1-1  /12+n1-1 + n2-1
!    8+n1-1  /9 +n1-1
!    9+n1-1  /
!   10+n1-1  expanded sub2
!    ..
!   10+n1-1+n2-1  last of sub2
!   11+n1-1+n2-1 ..
!   ..
!          
!     aux  variables
!        
!        SubDUsed  PosInDet    SubDNumb       SumOffset
!          1           5          2  (sub1)      n1-1
!          2           10         1  (sub2)      n1-1 +n2-1
!
!          ..
!
       integer n, m, subDx, subDy, j, i

       integer cons, consr, pcon, nmat
       logical firstcons, firstconsr, firstpcon
       integer::idx 
!!!!!!!!!v9.14
       integer kk, jj, ll
       character(8):: mrepl 
!!!!!!
! ^^^^^^^^^^^^^^^
       integer  subdidxc
       save  subdidxc
       subdidxc = 0
!^^^^^^^^^^^^^^^^


       do j = Det%nct, 1, -1

          cons = Det%cmp(j)%Contains
          consr= Det%cmp(j)%ContainsR
          pcon = Det%cmp(j)%PContained
          nmat = Det%cmp(j)%NMatreska

          call ep2actn(j, n, subDx)

          firstcons = .true.
          firstpcon = .true.
          firstconsr = .true.

          if(subDx .ne. 0) then
!              sub detector subDx is being included to
!              Det.cmp(m); m= n - SubD(subDx).nct + 1
             if(final == 0 ) then
                idx = maxSubD +1
             else
                idx = NsubD + 1  ! NsubD is not yet updated so +1
             endif
             if(final .ne. 0) then
                subdidxc = subdidxc + 1   ! count subd in main body
             endif
!           ^^^^^^^^^^^
!!!!!!!!!!!! v9.14
             mrepl = " "
             do kk = 1, containersCounter
                do ll = 1, noOfContainedSub(kk)
                   if( j ==  containedSubList(ll, kk) ) then
                      jj =  containersList(kk)
                      mrepl = containersMatter(kk)
                      exit
                   endif
                enddo
                if( mrepl /= " ") exit
             enddo
!!!!!!!!!!!!!!!!!!!!!!!!! v9.14
             call epIncSubD(subdidxc, ! v9.14
     *                  j, subDx, n - SubD(subDx)%nct +1, mrepl)
!                move comment pos
             if(wrtcom .ne. 0 ) then
                call epresetcom(j, n - SubD(subDx)%nct +1)
             endif
          else
!              move component
             call epmvComp(j, n, Det%cmp(j), Det%cmp(n), 0, NotGiven)
!             call epprocIDmv(j, n)
!                  move comment pos.
             if(wrtcom .ne. 0 ) then
                call epresetcom(j, n)
             endif
!              other adjustment

!             do i = 1, Det.cmp(j).NMatreska
             do i = 1, nmat
!((((((((((((((((((((((
!                call ep2actn(Det.cmp(j).Contains(i), m, subDy)
                call ep2actn( CnArea( cons+i ), m, subDy)
!cc                Det.cmp(n).Contains(i) = m
                if(firstcons) then
                   Det%cmp(n)%Contains = CnCounter
                   firstcons = .false.
                endif
!cccc             CnArea( Det.cmp(n).Contains+i) = m
                call epCnArea( m )
!                Det.cmp(n).ContainsR(i) = 
!     *            sign(m, Det.cmp(j).ContainsR(i))  ! this is next job
!))))))))))))))))))))))
             enddo

!c             do i = 1, Det.cmp(j).NMatreska
               do i = 1, nmat
!((((((((((((((((((((((
!                call ep2actn(Det.cmp(j).Contains(i), m, subDy)
                call ep2actn(CnArea( cons+i ), m, subDy)
!                Det.cmp(n).ContainsR(i) = 
!     *            sign(m, Det.cmp(j).ContainsR(i))
                if(firstconsr) then
                   Det%cmp(n)%ContainsR = CnCounter
                   firstconsr = .false.
                endif
                call epCnArea( sign(m, 
     *              CnArea( consr+i) ) )
!))))))))))))))))))))))
             enddo

!
             do i = 1, Det%cmp(j)%NPContainer
!(((((((((((((((((((((((
                if(firstpcon) then
                   Det%cmp(n)%PContained = CnCounter
                   firstpcon = .false.
                endif
!                call ep2actn(Det.cmp(j).PContained(i), m, subDy)
!                Det.cmp(n).PContained(i) = m
                call ep2actn(
     *              CnArea( pcon+i ), m, subDy)
                call epCnArea(m)
!))))))))))))))))
             enddo
          endif
       enddo
!         reset Det.nct
       call  ep2actn(Det%nct, j, i)
       Det%nct = j
       end
!     **********8
      subroutine epresetcom(j, n)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer j,n    ! input comment of j-th component is moved to
                     !       as for n-th comp.
      integer i

      do i = comcounter, 1, -1
         if(comflag(i) .eq. j) then
            comflag(i) = n
         elseif(comflag(i) .lt. j) then
            goto 10
         endif
      enddo
 10   continue
      end
!     *******************************                   v9.14
      subroutine epIncSubD(subdidxc, orgPos, subDx, pos, mrepl)
!!!      subroutine epIncSubD(subdidxc, orgPos, subDx, pos)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!     ^^^^^^^^^^^^
      integer subdidxc ! input.  incremental counter for subdetector
                       !         in the final phase

      integer orgPos  ! input.  pos in Det where the sub detector
                      !         is referred.
      integer subDx   ! input.  sub detector number to be 
                      !         included
      integer  pos    ! input.  Det.cmp(pos) is the first
                      !         position where sub detector
                      !         is to be moved.
!!!!! v9.14
      character(*),intent(in):: mrepl ! if /= ' ', replace
                     ! world matter by this one
!!!!! v9.14
      integer offset, i, j, temp
       type(Component):: Vcomp


       offset = pos -1

!         move included sub Det into Vcomp as a compoent.
!
       Vcomp = Det%cmp(orgPos)
!ccc       call epmvComp(orgPos, 0, Det.cmp(orgPos), Vcomp, 0)
!
!         
       
       j = pos -1
!         expand sub detector
!///////////////
       do i = 1, SubD(subDx)%nct
          j = j + 1
!                move  one comp of sub det.
!          call epmvComp(i, 0, SubD(subDx).cmp(i), Det.cmp(j), offset)
          call epmvComp(i, 0, 
     *    SubdArea( SubD(subDx)%loc+i ),  Det%cmp(j), offset,
     *     Vcomp%chno )
!    ^^^^^^^^^^^
          if( subdidxc .gt.  0) then
!!!!////             Det.cmp(j).subdidx = subDx
             Det%cmp(j)%fsubdc = subdidxc
          endif
          if(Vcomp%matter .ne. ' ' .and. 
     *       Vcomp%matter .ne. 'sp' .and.
     *       Vcomp%matter .ne. 'sp2' ) then
!              replace all matter of the subdetector by this one
             Det%cmp(j)%matter = Vcomp%matter
          endif
!     ^^^^^^^^^^^
!((((((((((
!ccc          Det.cmp(j).level = Det.cmp(j).level  + 1
!))))))))))
!           adjust  displacement: R0+ M^-1 r0 --> r0
          call epnewDisp(Vcomp, Det%cmp(j))
!           adjust  rotation matrix; mM --> m 
!                         M             m
          call epnewRot(Vcomp%direc, Det%cmp(j)%direc)
          Det%cmp(j)%rotation =
     *     Det%cmp(j)%direc(1)  .ne. 1.d0
     * .or.
     *     Det%cmp(j)%direc(5)  .ne. 1.d0
     * .or.
     *     Det%cmp(j)%direc(9)  .ne. 1.d0

!c wrong          Det.cmp(j).rotation = Vcomp.rotation 
       enddo
!           if  the last comp is world; _w is
!         dropped now
!       if(Det.cmp(j).struc(1:5) .eq. 'box_w') then
!          Det.cmp(j).struc = 'box'
!       elseif(Det.cmp(j).struc(1:8) .eq. 'sphere_w') then
!          Det.cmp(j).struc = 'sphere'
       temp =index(trim(Det%cmp(j)%struc), '_w')  ! trim v9.14
       if(temp .gt. 0) then
          Det%cmp(j)%struc = Det%cmp(j)%struc(1:temp-1)
!!!!!!!  v9.14
          if( mrepl /= " " .and. Det%cmp(j)%matter == "sp") then
             Det%cmp(j)%matter = mrepl
          endif
!!!!!
       endif
       end
      subroutine epnewDisp(cmpA, cmpB)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!         M = cmpA's rot matrix
!         R0 = disp vec of cmpA
!         r0 = disp vec of cmpB
!      
!         r0 <--  M^-1 r0 + R0 is obtained

       type(Component)::  cmpA, cmpB

      real*8 x, y, z
      
      x =  cmpA%direc(1)* cmpB%orgx + 
     *     cmpA%direc(4)* cmpB%orgy +
     *     cmpA%direc(7)* cmpB%orgz

      y =  cmpA%direc(2)* cmpB%orgx + 
     *     cmpA%direc(5)* cmpB%orgy +
     *     cmpA%direc(8)* cmpB%orgz


      z =  cmpA%direc(3)* cmpB%orgx + 
     *     cmpA%direc(6)* cmpB%orgy +
     *     cmpA%direc(9)* cmpB%orgz

       cmpB%orgx = cmpA%orgx  + x 
       cmpB%orgy = cmpA%orgy  + y
       cmpB%orgz = cmpA%orgz  + z
      end
      subroutine epnewRot(Ml, ms)
      implicit none
!         Ml 
!         ms 
!         ms <-- msMl is obtained
!      
      real*8 Ml(9), ms(9)
      integer i
      real*8  v(9)
      call epmatmul(Ml, ms, v)     
      do i = 1, 9
         ms(i) = v(i)
      enddo
      end

      subroutine eprecount
      implicit none
!cc      #include  "ZepMaxdef.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!            

!!!xxx     character(MAX_STRUCCHR)::epparaphrase
!!!xxx      character(MAX_STRUCCHR)::tempph

      integer  i, j

!      Det.nbox=0
!      Det.ncyl=0
!      Det.npip=0
!      Det.nprs=0
!      Det.nsph=0
      Det%nworld = 0
!!!xxx      do  i = 1, MaxNewStruc
!!!xxx         Det.nnew(i) = 0
!!!xxx      enddo
      do  i = 1, Det%nct
         Det%cmp(i)%cn = i  ! renumber
      enddo
           ! trim v9.14
      if(index(trim(Det%cmp(Det%nct)%struc), '_w') .ne. 0) then
         Det%nworld = 1
      endif
      end

      subroutine epresetMatter  ! v9.14
      implicit none
!cc     #include  "ZepMaxdef.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!            
      integer  i, j
      character(8):: tempc
      logical,save:: first=.true.

      if( Det%nworld > 0 ) then
         tempc = Det%cmp(Det%nct)%matter
      else
         tempc = "sp"
      endif
      do  i = 1, Det%nct-Det%nworld
         if( Det%cmp(i)%matter == "world" ) then
            Det%cmp(i)%matter = tempc
            if(first .and. Det%nworld == 0 ) then
               ! issue warning (only once)
               write(0,*) '***********WARNING****************'
               write(0,*) 'As medium of component #=',i
               write(0,*) '"world" is used, but there is no world'
               write(0,*) 'defined in the main detector,'
               write(0,*) 'so it is assumed to be "sp". However,'
               write(0,*) 'this may result in error at run time,'
               write(0,*) 'if there is no "sp" in other part of'
               write(0,*) 'the config file'
               first =.false.
            endif
         elseif( Det%cmp(i)%matter == "sp2" ) then
            Det%cmp(i)%matter = "sp"
         endif
      enddo
      end  subroutine epresetMatter  

      subroutine epclCnfCnt
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!           clear counters for Det.
!         This is used at the init. time of program run
!         and also after subdetector components are 
!         moved to subdetector area.
!         (Each subdetector is read once in Det. area and
!          then, moved to subdetector area)
!         When the filnal main detector is put in Det.
!         there is no need to  move components in Det.
!         so also the counters are not cleared here.
      integer  i

      Det%nct=0
      comcounter = 0
      comloc = 1
!      Det.nbox=0
!      Det.ncyl=0
!      Det.npip=0
!      Det.nprs=0
!      Det.nsph=0
      Det%nworld = 0
!          clear the counter of parenet matreskas
!          which partially /fully contain i-th comp.
      do i = 1, ncmax
         Det%cmp(i)%NPContainer = 0
         Det%cmp(i)%NContainer = 0
      enddo
      SubDUsed = 0


      end
      subroutine epsvcom(cmnt)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      character*(*) cmnt
      if(wrtcom .ne.0) then
         if(comcounter .lt. maxcom) then
            comcounter = comcounter + 1
            comarea(comcounter) = cmnt(1:80)
            comflag(comcounter) = Det%nct + 1
         else
            call cerrorMsg('too many comment in config',1)
         endif
      endif
      end
!      
      subroutine epmvDet(det1, comp, subdet)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
       type(Detector)::  det1  ! input
       type(Component):: comp(*)  ! ouptut
       type(SubDetector)::  subdet  ! output

      integer i
!        move det1 into subdetector
      
      if(det1%nct + cumsubdloc .gt. ncmaxInSubD) then
         call cerrorMsg(
     *   'MAX_COMPONENT_IN_SUBD in ZepMaxdef%h is too small',0)
      endif


      do i = 1, det1%nct
!cc         call epmvComp(i, 0, det1.cmp(i), det2.cmp(i),  0)
!         det2.cmp(i) = det1.cmp(i)
         comp(i) = det1%cmp(i)
      enddo

!      subd.nbox = det1.nbox
!      subd.ncyl = det1.ncyl
!      subd.npip = det1.npip
!      subd.nprs = det1.nprs
!      subd.nsph = det1.nsph

!!!xxx      do i = 1, MaxNewStruc
!!!xxx         subdet.nnew(i) = det1.nnew(i)
!!!xxx      enddo
      subdet%nworld = det1%nworld
      subdet%nct = det1%nct
      cumsubdloc = cumsubdloc + subdet%nct
      end

      subroutine epmvComp(n1, n2, cmpA, cmpB, offset, chno )
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
       type(Component)::  cmpA, cmpB
      integer n1, n2  ! input. n1==n2 ==> *cmpA==*cmpB
                      !   n1!=n2==> *cmpA != *cmpB even if cmpA.cn ==cmpB.cn
      integer offset, j
      integer chno    ! input. if this is not NotGiven, copied to cmpB
      integer cons, consr, pcon

      if(n1 .ne. n2) then
         cmpB = cmpA 
         if(chno .ne. NotGiven) then
            cmpB%chno = chno
         endif
      endif
!         cmpB.orgx = cmpA.orgx
!         cmpB.orgy = cmpA.orgy
!         cmpB.orgz = cmpA.orgz
!         do j = 1, maxattr
!            cmpB.vol(j) = cmpA.vol(j)
!         enddo
!         do j = 1, 9
!            cmpB.direc(j) = cmpA.direc(j)
!         enddo
!         cmpB.NMatreska = cmpA.NMatreska
!         cmpB.NContainer = cmpA.NContainer
!((((((((((((((((((((
!
!          save cmpA index, so that even if cmpB==cmpA 
!          we can move cmpA safely.
!
         cons = cmpA%Contains
         consr = cmpA%ContainsR
         pcon = cmpA%PContained
!
         do j = 1, cmpA%NMatreska
!            cmpB.Contains(j) =
!     *         cmpA.Contains(j) + offset
!            cmpB.ContainsR(j) =  sign(cmpB.Contains(j),
!     *         cmpA.ContainsR(j)) 
            if(j .eq. 1) then
                cmpB%Contains = CnCounter
             endif
             call epCnArea( CnArea(cons +j ) + offset)
         enddo
         do j = 1, cmpA%NMatreska
!            cmpB.ContainsR(j) =  sign(cmpB.Contains(j),
!     *         cmpA.ContainsR(j)) 
            if(j .eq. 1) then
               cmpB%ContainsR =CnCounter
            endif
            call epCnArea(
     *           sign(CnArea(cmpB%Contains+j), CnArea(consr+j)))
         enddo
!))))))))))))))))))))
!         cmpB.NPContainer = cmpA.NPContainer
         do j = 1, cmpA%NPContainer
!(((((((((((((((((((
!            cmpB.PContained(j) =
!     *         cmpA.PContained(j) + offset
            if(j .eq. 1 ) then
               cmpB%PContained = CnCounter
            endif
            call epCnArea( CnArea( pcon+j ) + offset)
!)))))))))))))
         enddo
!         cmpB.struc = cmpA.struc
!         cmpB.strucNo = cmpA.strucNo
!         cmpB.matter = cmpA.matter
!         cmpB.rotation = cmpA.rotation
!         do j = 1, MaxNewStruc
!            cmpB.nnew(j) = cmpA.nnew(j)
!         enddo
      end

      subroutine epadjPC
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!            adjust partially contained relation;
!        if a component is contained by another one which is
!        contaied by still other one, ..., and the first
!        one is partially contained, all are assumed to be contained
!        parially by the parent:
      integer maxstack, maxpc
      parameter (maxstack = 1000, maxpc=10000)
      integer stack(maxstack), stackp
      integer k2i(maxpc, 2), k2ic
      integer used(maxpc), usedc
!
      integer i, j, k, m
!          reset NPContainer
      do i = 1, Det%nct
         Det%cmp(i)%NPContainer = 0
      enddo
      k2ic= 0
      do i = 1, Det%nct - Det%nworld
         do j = 1, Det%cmp(i)%NMatreska
!((((((((((((((((
!            if(Det.cmp(i).ContainsR(j) .lt. 0) then
!               k = - Det.cmp(i).ContainsR(j)
            if(CnArea( Det%cmp(i)%ContainsR+j ) .lt. 0) then
               k = - CnArea( Det%cmp(i)%ContainsR+j )
!))))))))))))
!                The following business could be done
!                by recursive call very easily but
!                we use stack method for wider compatibilty
!
               stackp = 1
!                    put k in stack
               stack(stackp) = k
!                extrack from stack
               do while(stackp .gt. 0)
                  k = stack(stackp)
                  stackp = stackp - 1
!                    set parent i  in k-th comp
                  Det%cmp(k)%NPContainer =
     *                  Det%cmp(k)%NPContainer + 1
!                  Det.cmp(k).PContained( Det.cmp(k).NPContainer )
!     *               = i
                  k2ic = k2ic + 1

                  if(k2ic .gt. maxpc) then
                     call cerrorMsg(
     *              ' too many partially contained cmp', 0)
                  endif
                  used(k2ic) =0
                  k2i(k2ic,1) = k
                  k2i(k2ic,2) = i
                  do m = 1, Det%cmp(k)%NMatreska
                     stackp = stackp + 1
!(((((((((((
!                     stack(stackp) = Det.cmp(k).ContainsR(m)
                     stack(stackp) =CnArea(Det%cmp(k)%ContainsR+m)
!)))))))))))))
                  enddo
               enddo
            endif
         enddo
      enddo
      usedc = 0
      do i = 1, k2ic
         k = k2i(i, 1)
         if(k .gt. 0) then
            do j = 1, usedc
               if(used(j) .eq.  k) goto 10
            enddo
!              k is not used yet
            usedc = usedc + 1
            used(usedc) = k
            Det%cmp(k)%PContained = CnCounter
 10         continue
            do j = i, k2ic
               if( k2i(j, 1) .eq. k ) then
                  call epCnArea( k2i(j, 2))
                  k2i(j, 1) = 0
               endif
            enddo
         endif
      enddo
            
      end
      subroutine epworld
       use  modSubdWorld, only:bugbug
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      
!          make the world to contain all  comp.      
!           which are not contained by other 
!           comp.
!          partially contained one is not added in the
!        list.
! 
      integer i, k, epIsSubD, j
      character*200 msg
      logical contained
      
!          trim v9.14
      if(index(trim(Det%cmp(Det%nct)%struc), '_w') .ne. 0 .and.
     *     Det%nworld .eq. 1) then
!             ok
      else
         call cerrorMsg(
     *    'world must be the last only 1 component', 0)
      endif
!!        now world has no inner comps.  V9.157 to permit world force contain 

      do i =  Det%nct-1, 1, -1  
         if(Det%cmp(i)%NContainer .eq. 0 ) then
!              partially contained one is not included
!              check if i-th compnent appears in some matreska
!          evenif NContainer =0, sometimes already it's matreshkar 

            call epCheckMat(i, contained)   ! v9.14

            if( .not. contained )  then   ! v9.14
               Det%cmp(Det%nct)%NMatreska =
     *              Det%cmp(Det%nct)%NMatreska + 1
               k = Det%cmp(Det%nct)%NMatreska
               if( Det%cmp(Det%nct)%Contains .eq. -1) then
!                  two can use  the same array point
                  Det%cmp(Det%nct)%Contains =CnCounter
                  Det%cmp(Det%nct)%ContainsR =CnCounter
               endif

               call epCnArea(i)
               j = epIsSubD(i)  ! subD number
               if( j .gt. 0 )  then
!              i is subD.
                  if( SubD(j)%nct .gt. 1) then ! more than one comp. in subD
                     if( SubD(j)%nworld .gt. 0) then
!                       o.k
!                 elseif(SubD(j).cmp(SubD(j).nct).NMatreska .eq. 0) then
                     elseif( 
     *                 SubdArea(SubD(j)%loc+ SubD(j)%nct )%NMatreska
     *               .eq. 0) then
!                         last one seems a simple component
                        msg='subdetetor( '//SubDName(j)//') must have '
     *      //'implicit or explicit world at its last since it is used '
     *      //'by another (sub)Detctor with a world'
                        call cerrorMsg(msg, 0)
                     else
!                      check if the last one is equivalent to world
                        call episeqvWorld(
     *                       j, SubD(j), SubdArea( SubD(j)%loc+1 ) )
                     endif
                  endif
               endif
            endif  !v9.14
!                 now i is contained by the world.or other comp set flag
            Det%cmp(i)%NContainer = 1 ! 
         endif 
      enddo
      end
!!!!!!!!!!!!  next sub v9.14
      subroutine epCheckMat(compn, contained)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer,intent(in):: compn ! component number which is checked
                          ! whether it is listed in matreska
      logical,intent(out):: contained ! if compn is not listed as mat
                             !   .false. else .true. is returned

      integer::i, j, matcompn
      integer::NoOfMat

      contained = .false.
!!!!      do i = Det.nct-1, 1, -1   !next  v9.157 (to exclude the one if already
!11                     listed in world)
      do i = Det%nct, 1, -1
         NoOfMat= Det%cmp(i)%NMatreska 
         do j = 1, NoOfMat
            ! component # of j-th matreska
            matcompn = CnArea(Det%cmp(i)%ContainsR+j)
            if( matcompn .eq. compn) then
                ! i is aleady listed as matreska
               contained =.true.
               return  !!!!!!!!!!!!!
            endif
         enddo
      enddo
      end


      subroutine episeqvWorld(n, subdin, comp)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer n  ! input. subdetector number
       type(SubDetector)::  subdin   ! input. a subDetector
       type(Component):: comp(*)
!
!          examines if the last component is equivalent to
!          a world
!     
      integer*2 numb(ncmaxInSubD)
      
      integer i, m, mt, j
      character*200 msg
      logical error
      do i = 1, subdin%nct
         numb(i) = 0
      enddo
      error = .false.
      do i =1, subdin%nct
         mt = comp(i)%NMatreska 
         if(mt .gt. 0) then
!            m = subdin.cmp(i).Contains
            m = comp(i)%Contains
            do j = 1, mt
                numb(CnArea( m+j)) =1
            enddo
         endif
      enddo
      do i = 1, subdin%nct-1
         if(numb(i) .eq. 0) then
            write(msg, *) ' subdetector=',SubDName(n), 
     *     ':  implicit world  dose not contain ', i, '-th',
     *     ' component directoy or indirectory'
           call  cerrorMsg(msg, 1)
           error = .true.
         endif
      enddo
      if(error) then
         call cerrorMsg(
     *    'Some implicit world in some subdetectors may not'
     *   //' contain all its child', 1)
         call cerrorMsg(' check the subdetector',0)
      else
         call cerrorMsg('All implicit worlds are O.K', 1)
      endif
      end
      subroutine eprotation(comp)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepDirec.h"

       type(Component)::  comp  ! input/output.
      character*150 msg

       type(epDirec)::  wx, wy, wz

      real*8 a, b, c, epgetSProd

       comp%rotation =
     *     comp%direc(1)  .ne. 1.d0
     * .or.
     *     comp%direc(5)  .ne. 1.d0
     * .or.
     *     comp%direc(9)  .ne. 1.d0

      if(comp%rotation ) then
!          check rotation matrix         
         wx%x = comp%direc(1)
         wx%y = comp%direc(2)
         wx%z = comp%direc(3)

         wy%x = comp%direc(4)
         wy%y = comp%direc(5)
         wy%z = comp%direc(6)

         wz%x = comp%direc(7)
         wz%y = comp%direc(8)
         wz%z = comp%direc(9)

         a = epgetSProd(wx, wx)
         b = epgetSProd(wy, wy)
         c = epgetSProd(wz, wz)


         if(wx%x .eq. 1.d0 .and. wy%y  .eq.  1.d0) then
!                    this wuold be inclined cylinder etc
!              c should be 1.0
            if(abs(c-1.d0) .gt. MaxErrDirCos) then
               write(msg,*)' sum of direcos Z**2 for ',
     *         Det%nct, '-th comp. is wrong=',c
               call cerrorMsg(msg, 1)
               write(msg, *) ' Zx,Zy,Zz=',
     *         wz%x, wz%y, wz%z
               call cerrorMsg(msg, 0)
            else
!                make a, b
               if( wz%z  .eq. 1.d0) then
                  comp%rotation = .false.
               else
!                          force to  have sum**2=1
                  call epadjdir(comp%direc(7)) 
!                         make wx, wy
                  call epmakexy(comp%direc)
               endif
            endif
         else
            if(abs(a-1.d0) .gt. MaxErrDirCos .or.
     *         abs(b-1.d0) .gt. MaxErrDirCos ) then
               write(msg,*) ' sum dircos X**2=',a,
     *         ' or  sum dircos Y**2=',b, ' for ',
     *         Det%nct, '-th comp. wrong '
               call cerrorMsg(msg, 1)
               write(msg, *) ' Xx,Xy,Xz=', wx%x, wx%y, wx%z
               call cerrorMsg(msg, 1)
               write(msg,*)
     *          ' Yx,Yy,Yz=', wy%x, wy%y, wy%z
               call cerrorMsg(msg, 0)
            elseif(abs( epgetSProd(wx, wy)) .gt.  EpsFor90) then
               call cerrorMsg(
     *         ' direc cos X, Y is not perpendicular',1)
               write(msg,*) ' X*Y=',   epgetSProd(wx, wy)
               call cerrorMsg(msg, 1)
               write(msg, *) ' Xx,Xy,Xz=', wx%x, wx%y, wx%z
               call cerrorMsg(msg, 1)
               write(msg,*) ' Yx,Yy,Yzs=', wy%x, wy%y, wy%z
               call cerrorMsg(msg, 0)
            else
               call epadjdir(comp%direc(1))
               call epadjdir(comp%direc(4))
!                   form Z direction
               call epmakez(comp%direc(1))

            endif
         endif
      endif
      end
!     ********************
      real*8 function epgetSProd(a, b)
      implicit none
#include "ZepDirec.h"
       type(epDirec)::  a, b
      epgetSProd = a%x*b%x + a%y*b%y + a%z*b%z
      end

      subroutine epadjdir(m)
!         adjust direction cos so that its square sum
!         be 1.0
      implicit none
      real*8 m(3)
!         
      real*8 sum

      sum = m(1)**2 + m(2)**2 + m(3)**2
      if(sum .ne. 0.d0) then
         sum = sqrt(sum)
         m(1) = m(1)/sum
         m(2) = m(2)/sum
         m(3) = m(3)/sum
      else
         m(1) = 1.
         m(2) = 0.
         m(3) = 0.
      endif
      end

      subroutine epmakez(m)
      implicit none
      real*8 m(9)
!           m(1-3) are  dire.cos of X'
!           m(4-6) are  //          Y'
!           m(7-9) //               Z'
!        make Z'=X' x Y'
!       Z'x = X'y Y'z     X'z  Y'y
      m(7) = m(2)* m(6) - m(3)*m(5)
!       Z'y = X'z Y'x     X'x  Y'z
      m(8) = m(3)* m(4) - m(1)*m(6)
!       Z'z = X'x Y'y     X'y  Y'x
      m(9) = m(1)* m(5) - m(2)*m(4) 
      end
!     **********************
      subroutine epmakexy(m)
      implicit none
      real*8 m(9)
!         X' = Z' x z :  
!        X'x = Z'y 1 -Z'z 0
       m(1) = m(8)
!        X'y = Z'z 0 - Z'x 1
       m(2) = -m(7)
!        X'z = Z'x 0 - Z'y 0
       m(3) = 0.
!           normalize
       call epadjdir(m(1))
!        Y' = Z' x X'
!        Y'x = Z'y X'z - Z'z X'y   
       m(4) = m(8)*m(3) - m(9)*m(2)
!        Y'y = Z'z X'x - Z'x X'z   
       m(5) = m(9)*m(1) - m(7)*m(3)
!        Y'z = Z'x X'y - Z'y X'x   
       m(6) = m(7)*m(2) - m(8)*m(1)
       end
       subroutine eprbox(comp)
       use modNegativeC
       implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"

           integer  klena
           character*80 msg
           real(8):: a, b, c

           integer:: sidx
       type(Component)::  comp  ! input/output. where component data 
                                    !      is to read in.
!              read part of the structure
           call eprpst(comp, 3, 3, 1, 6 )
!!           check some values
           a = Volat( comp%vol + 1)
           b = Volat( comp%vol + 2)
           c = Volat( comp%vol + 3)
           if( a < 0. .or. b < 0. .or. c < 0. ) then
              negativeABC = negativeABC + 1
              sidx = Det%cmp(Det%nct)%subdidx 
              if( sidx > 0 ) then
                 write(0,*)
     *           '***** Info: negative box edge in=', SubdName(sidx)
              endif
           endif
           if(a < 0. ) then
              comp%orgx = comp%orgx + a
              a = -a
              Volat( comp%vol +1 ) =a
           endif
           if(b < 0.) then
              comp%orgy = comp%orgy + b
              b = -b
              Volat( comp%vol +2 ) = b
           endif
           if( c < 0.) then
              comp%orgz = comp%orgz + c     
              c = -c
              Volat( comp%vol +3 ) = c
           endif

           if( a == 0. .or. b == 0. .or. c==0. ) then
              if( trim(Field(2)) /= 'box_w'   ) then
                 write(0,*) ' box edge length a, b, c=',a,b,c
                 write(0,*) ' must be > 0'
                 stop
              endif
           endif
!           Det.nbox = Det.nbox+1
           return
!     ***************** cylinder
      entry eprcyl(comp)
!     ****************
!              read part of the structure

           call eprpst(comp, 2, 2, 7, 9)   ! 7,9 means only axis of cyl 
!!           check some values
           a = Volat( comp%vol + 1)  !  r
           b = Volat( comp%vol + 2)  !  h

           if( a== 0. .and. b == 0.) then
              if( index(Field(2), '_w' ) == 0 ) then
                 write(0,*) ' cyl r or heigth r, h=',a, b
                 write(0,*) 'invalid'
                 stop
              else
                 return !!!!!!!!!!!!
              endif
           endif
           if( a == 0. .or. b == 0. ) then
              write(0,*) ' cyl r or heigth r, h=',a, b
              write(0,*) 'invalid'
              stop
           elseif( a <  0. .or. b == 0. ) then
              write(0,*) ' cyl r or heigth r, h=',a, b
              write(0,*) ' invalid'
              stop
           endif

!           Det.ncyl=Det.ncyl+1
!              variants which is nothing but a rotation of 
!              the canonical form
!           call eprotation(comp)
           if( b < 0.d0 ) then
              sidx = Det%cmp(Det%nct)%subdidx 
              if( sidx > 0 ) then
                 write(0,*)
     *          '***** Info: negative Cyl height in=', 
     *          SubdName(sidx)
              endif
           endif
           if(index(confdata, 'cyl_z') .ne. 0 .or.
     *        Field(2)(1:klena(Field(2))) .eq. 'cyl'  .or.
     *        Field(2)(1:klena(Field(2))) .eq. 'cyl_w'   ) then
!               this is canonical form; nothing to do
              if( b < 0. ) then
                 comp%orgz = comp%orgz + b
                 b = -b
                 Volat( comp%vol +2 ) = b
                 negativeCylH = negativeCylH  + 1
              endif

           elseif(index(confdata, 'cyl_x') .ne. 0) then
!
!             canonical form is rotated by 90 deg around y

!              call eprotdirec(comp, 2,  0) 
              if( b < 0. ) then
                 comp%orgx = comp%orgx + b
                 b = -b
                 Volat( comp%vol +2 ) = b
                 negativeCylH = negativeCylH  + 1
              endif

           elseif(index(confdata,'cyl_y') .ne. 0) then
              if( b < 0. ) then
                 comp%orgy = comp%orgy + b
                 b = -b
                 Volat( comp%vol +2 ) = b
                 negativeCylH = negativeCylH  + 1
              endif
!               
!              call eprotdirec(comp, 3, 2)
           else
              write(msg,*) ' error structure=', confdata
              call cerrorMsg(msg, 0)
           endif
           return
!     ***************** pipe
      entry eprpip(comp)
!     ****************
!              read part of the structure
           call eprpst(comp, 3, 3, 7, 9)
           a = Volat( comp%vol + 1)
           b = Volat( comp%vol + 2)
           c = Volat( comp%vol + 3)
           if( a == 0. .and. b==0. .and. c==0.) then
              if( index(Field(2), '_w') == 0 ) then
                 write(0,*) ' pipe  attribute 0'
                 stop
              else
                 return  !!!!!!!!!!!
              endif
           endif
           if( c < 0.d0 ) then
              sidx = Det%cmp(Det%nct)%subdidx 
              if( sidx > 0 ) then
                 write(0,*)
     *          '***** Info: negative Pipe height in=',
     *           SubdName(sidx)
              endif
           endif

           if(index(confdata, 'pipe_z') .ne. 0 .or.
     *        Field(2)(1:klena(Field(2))) .eq. 'pipe'     ) then
!               this is canonical form; nothing to do
              if( c < 0. ) then
                 comp%orgz = comp%orgz + c     
                 c = -c
                 Volat( comp%vol +3 ) = c
                 negativePipeH = negativePipeH  + 1
              endif
           elseif(index(confdata, 'pipe_x') .ne. 0) then
              if( c < 0. ) then
                 comp%orgx = comp%orgx + c     
                 c = -c
                 Volat( comp%vol +3 ) = c
                 negativePipeH = negativePipeH  + 1
              endif
           elseif(index(confdata,'pipe_y') .ne. 0) then
              if( c < 0. ) then
                 comp%orgy = comp%orgy + c     
                 c = -c
                 Volat( comp%vol +3 ) = c
                 negativePipeH = negativePipeH  + 1
              endif
           else
              write(msg,*) ' error structure=', confdata
              call cerrorMsg(msg, 0)
           endif
              ! for safety. 
           if( Volat(comp%vol+pipeir) >= 
     *              Volat(comp%vol+pipeor) ) then
              write(0,*) 
     *       ' inner pipe radius must be < outer radius; they are:'
              write(0,'(1p,2g15.6)')  
     *        Volat(comp%vol+pipeir), Volat(comp%vol+pipeor)
              stop
           endif

           if(a <= 0. .or. b <=0. .or. c==0.) then
              write(0,*) ' error in pipe radius and/or height'
              write(0,*) ' r1, r2, h=', a, b, c
              stop
           endif


           return
!     ********************
      entry eprsph(comp)
!              read part of the structure
           call eprpst(comp, 1, 1, 1, 1)
!                this should be ... 1, 0, 0)  but it will 
!          generate array boundary error, so put 1,1 
!          they are neglected.
           a = Volat( comp%vol + 1)
           if( a < 0. ) then
              write(0,*) ' sphere radius <= 0', a
              stop
           elseif( a == 0. ) then
              if( trim(Field(2)) /= 'sphere_w'   ) then
                 write(0,*) ' sphere radius = 0', a
                 stop
              endif
           endif
!           Det.nsph = Det.nsph+1
           return

      end
      subroutine eprotdirec(comp, j2, j1)
       implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
       type(Component)::  comp  ! input/ouput
       integer j2, j1
       integer i

!        If the rotation matrix (direc) is  not a unit
!       matrix, this has priority and M1,M2 is neglected.
!        If it is the unit matrix, 
!       it is changed to
!       direc =  M2 x M1 where Mi is the rotation matrix
!       for the anti-clockwise rotation around ji-th
!       axis by 90 degrees.
!       if ji=0, no rotaion.  ji=1 ==>  x axis, ji=2 ==> y,
!                             ji=3 ==>  z axis.
!
 

       real*8 rmat(9, 4), M21(9)
       real*8 o, u, mu
       parameter(o=0.d0, u=1.d0, mu=-u)
!        Solaris refuse to use data x/-u/ so that mu is used.
       data  (rmat(i,2),i=1, 9)/   !  around x
     *   u, o,  o,
     *   o, o,  u,
     *   o, mu, o/
       data  (rmat(i,3), i=1, 9)/  !  around y
     *   o, o, mu,
     *   o, u, o,
     *   u, o, o/
       data  (rmat(i,4), i=1, 9)/   !  around z
     *    o, u, o,
     *   mu, o, o,
     *    o, o, u/
       data  (rmat(i,1), i=1, 9)/   !   unit matrix
     *    u, o, o,
     *    o, u, o,
     *    o, o, u/



       if(comp%direc(1) .eq. 1.d0 .and.
     *    comp%direc(5) .eq. 1.d0 .and.
     *    comp%direc(9) .eq. 1.d0 )  then
!            M21= M2 x M1
          call  epmatmul(rmat(1, j2+1),  rmat(1, j1+1),  M21)
          do i= 1, 9
             comp%direc(i) = M21(i)
          enddo
       endif
       end
      subroutine epirotdirec(comp, j2, j1)
       implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
       type(Component)::  comp  ! input/ouput
       integer j2, j1
       integer i
!          Inverse of eprotdirec
!        If the rotation matrix (direc) is  not a unit
!       matrix, this has priority and M1,M2 is neglected.
!        If it is the unit matrix, 
!       it is changed to
!       direc =  M2 x M1 where Mi is the rotation matrix
!       for the anti-clockwise rotation around ji-th
!       axis by 90 degrees.
!       if ji=0, no rotaion.  ji=1 ==>  x axis, ji=2 ==> y,
!                             ji=3 ==>  z axis.
!
 

       real*8 rmat(9, 4), M21(9)
       real*8 o, u, mu
       parameter(o=0.d0, u=-1.d0, mu=u)

       data  (rmat(i,2),i=1, 9)/   !  around x
     *   u, o,  o,
     *   o, o,  u,
     *   o, mu, o/
       data  (rmat(i,3), i=1, 9)/  !  around y
     *   o, o, mu,
     *   o, u, o,
     *   u, o, o/
       data  (rmat(i,4), i=1, 9)/   !  around z
     *    o, u, o,
     *   mu, o, o,
     *    o, o, u/
       data  (rmat(i,1), i=1, 9)/   !   unit matrix
     *    mu, o, o,
     *    o, mu, o,
     *    o, o, mu/



       if(comp%direc(1) .eq. 1.d0 .and.
     *    comp%direc(5) .eq. 1.d0 .and.
     *    comp%direc(9) .eq. 1.d0 )  then
!            M21= M2 x M1
          call  epmatmul(rmat(1, j2+1),  rmat(1, j1+1),  M21)
          do i= 1, 9
             comp%direc(i) = M21(i)
          enddo
       endif
       end
       subroutine epmatmul(m1, m2, m3)
       implicit none
       real*8 m1(9), m2(9), m3(9)
!           m3 = m2  x  m1
!        m3 must not be  m1 or m2
!
       integer i, j, k
       real*8 sum
       do i = 1, 3   !  row
          do j =1, 3  !  col.
             sum = 0.
             do k = 1, 3
                sum = sum + m2( (i-1)*3 + k ) *
     *                      m1( (j + (k-1)*3 ))
             enddo
             m3((i-1)*3 + j) = sum
          enddo
       enddo
       end


      subroutine eprpst(comp, nattrin, mattrin,  ndir1, ndir2)
       implicit none
!cc       #include  "ZepMaxdef.h"
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
!          process a part of config data; common to all 
!          structures.
!
       type(Component)::  comp  ! input/output.  
       integer nattrin  !  input.  number of volume attributes of this 
                      !          structure
       integer mattrin  !  input. nattr + number of some values  dereived
                      !        from the attributes and used later
       integer ndir1, ndir2 ! input. ndir1-th to ndir2-th component of
                            !  direction cosines which specify the
                            !  volume orientation are read

       integer ncn, l, klena, alen, l2, nattr, mattr
       character*80 msg

       character(MAX_STRUCCHR)::epparaphrase
       character(MAX_STRUCCHR)::tempph

       integer noldnames, i, icon, inttemp
       parameter (noldnames=17)
       character*8 translate(2, noldnames)
       logical kalpha
       real*8 realtemp
       integer:: uscl
!          names used in old version should be accepted as 
!          new names: translation table
       data translate/
     *  'pb',         'Pb',
     *  'air',        'Air',
     *  'c',          'Carbon',
     *  'h2o',        'H2O',
     *  'al',         'Al',
     *  'ar',         'Ar',
     *  'bgo',        'BGO',
     *  'csi',        'CsI',
     *  'cu',         'Cu',
     *  'fe',         'Fe',
     *  'g5',         'G5',
     *  'g10',        'Gten',
     *  'w',          'W',
     *  'scin',       'SCIN',
     *  'lixe',       'LiXe',
     *  'liar',       'LiAr',
     *  'si',         'Si'/
!           next 2 are to get fpolygon's height and npoly
       character*24 orgF(5)
       integer norgF
       integer npoly
       read(confdata, *) ncn
       if(comp%cn .ne. ncn .and. .not. Incused ) then
          write(msg, *) ' check comp. #=',ncn,
     *         ' which should be', comp%cn
          call cerrorMsg(msg, 1)
          write(0,*) ' input line: ',confdata
          stop
       endif

       l = index(confdata, '/')
       if(l .eq. 0) then
          write(msg, *) '/ is missing in',ncn,' line'
          call cerrorMsg(msg, 0)
       endif
       call kgetField(confdata, orgF, 2,  norgF)
       if(orgF(2)(1:8) == "fpolygon")  then
!               npoly are  not yet known; get it first
          call kgetField(confdata(l+1:), orgF, 4,  norgF)
          if( norgf /= 4  ) then
             write(0,*) ' strange: no npoly data for fpolygon'
             stop
          endif
          read(orgF(4),*) npoly
          nattr = npoly *2 + 2  !  npoly + height + verteces
          mattr = nattr + 4     ! xmin,xmax, ymin, ymax will be
                                ! added 
       else
          nattr= nattrin
          mattr = mattrin
       endif
       
       if(nattr .gt. 0)  then
          if(mattr .ge. nattr) then
             call epAttribWatch(mattr, comp%vol)
          else
             call cerrorMsg('mattr < nattr for eprpst',0) 
          endif
       endif
       comp%Nattributes = nattr
!              this check is moved to the top
!       l = index(confdata, '/')
!       if(l .eq. 0) then
!          write(msg, *) '/ is missing in',ncn,' line'
!          call cerrorMsg(msg, 0)
!       endif
!     keep Upper/lower case as they are.
   ! treat  Air*0.98 type;  0.98 specifies density calibration factor
!           =   (rho of this Air)/(rho of  Air without *)
       call epgetRhoc(FieldAsItis(3), comp%matter, comp%rhoc)
!!!!  followings are obso. replaced by epgetRhoc. (Jul. 2016)
!       l = index(FieldAsItis(3), '*')
!       if(l .gt. 0) then
!          comp%matter=  FieldAsItis(3)(1:l-1)
!         read(FieldAsItis(3)(l+1:klena(FieldAsItis(3))),*) comp%rhoc
!       else
!          comp%matter =
!     *      FieldAsItis(3)(1:klena(FieldAsItis(3)))
!          comp%rhoc = 1.0
!       endif

!      ^^^^^^^^^^^^^^
       if( comp%matter .eq. "0") then
!             should be subdetector citation. use the same matter  as the
!             definition
!             format is like this:    1 qring 0 0 0 0 / 0 0 0.
!             if format is like       1 qring Al 0 0 0 0 / 0 0 0
!             all matter in the subdetector will be replaced by Al
!             when it is expanded. However,
!             if format is like       1 qring sp 0 0 0 0 / 0 0 0
!             the original matter is kept. (same as the first case).

          if(nattr .eq. 0) then   ! subdetector is called
             comp%matter = " "
          else
             write(msg, *)
     *        ' matterial must be given for ',comp%cn, '-th comp.'
             write(*, '("data =", a)') confdata(1:klena(confdata))
             call cerrorMsg(msg, 0)
          endif
       endif
!      ^^^^^^^^^^      
       comp%MaxPathL = 0.
       comp%Modifier = 0
       comp%chno = NotGiven
!     4 box Pb c dE mp mo chno
!     1  2  3  4 5  6   7  8 / 
       if(Field(4) .eq. '=') then
          comp%CountIO = EqualFshort
       else
          if(kalpha(Field(4))) then
             call epexamEqI(Field(4), inttemp, icon)
!             call epexamEqI(Field(4), comp.CountIO, icon)
             comp%CountIO = inttemp
             if(icon .ne. 0) then
                call eperrEq(comp, Field(4))
             endif
          else
             read(Field(4), *) comp%CountIO
          endif
       endif
!     4 box Pb c dE mp mo chno
!     1  2  3  4 5  6   7  8 / 
       if(Field(5) .eq. '=') then
          comp%CountDE = EqualFshort
       else
          if(kalpha(Field(5))) then
!             call epexamEqI(Field(5), comp.CountDE, icon)
             call epexamEqI(Field(5), inttemp, icon)
             comp%CountDE = inttemp
             if(icon .ne. 0) then
                call eperrEq(comp, Field(5))
             endif
          else
             read(Field(5), *) comp%CountDE
          endif
       endif
!          up to countDE has been processed
!     4 box Pb c dE mp mo chno
!     1  2  3  4 5  6   7  8 / 
!          4  5 6  7  8       4 5  6    4 5     4 5  6  7   4 5  6 
!          c de mp mo ch /    c de /    c de/   c de mp /  c  de mp/
!          4 5 6 7  8       4  5  6  7   
!          c de mp mo /     c de mp mo/
       if( index(Field(5), "/" ) > 0 ) then
!               countDE/ type
!         no  more maxpath, modifier, chno
       elseif( Field(6) == '/' ) then
!               countDE / type
!             no more data
       elseif(kalpha(Field(6) )) then
!          call epexamEq(Field(6), comp.MaxPathL, icon)
          call epexamEq(Field(6), realtemp, icon)
          comp%MaxPathL = realtemp
          if(icon .ne. 0) then
             call eperrEq(comp, Field(6))
          endif
       elseif( Field(6) == "=" ) then
          comp%MaxPathL = EqualF
       else
!                   may be  num or num/
             read(Field(6), *) comp%MaxPathL
       endif
!     
!          up to maxpath has been processed
!     4 box Pb c dE mp mo chno
!     1  2  3  4 5  6   7  8 / 
       if( index( Field(6), '/')   > 0 ) then
!             MaxPathL/ ; so no more data
       elseif( Field(7) .eq.  '/' ) then
!             MaxPathL given. but modifier not given
       elseif( kalpha(Field(7)) ) then
          call epexamEqI(Field(7), inttemp, icon)
          comp%modifier = inttemp
          if(icon .ne. 0) then
             call eperrEq(comp, Field(7))
          endif
       elseif( Field(7) == '=' ) then
          comp%modifier = EqualFshort
       else
          read( Field(7), *) comp%modifier
       endif
!     
!          up to modifier has been processed
!     4 box Pb c dE mp mo chno
!     1  2  3  4 5  6   7  8 / 
       if( index( Field(7), '/')   > 0 ) then
!             modifier/ ; so no more data
       elseif( Field(8) .eq.  '/' ) then
!             modifiere given. but chno not given
       elseif( kalpha(Field(8)) ) then
          call epexamEqI(Field(8), inttemp, icon)
          comp%chno = inttemp
          if(icon .ne. 0) then
             call eperrEq(comp, Field(8))
          endif
       elseif( Field(8) == '=' ) then
          comp%chno = EqualFshort
       else
          read( Field(8), *) comp%chno
       endif


       comp%struc = Field(2)(1:klena(Field(2)))

!           set number for the structure
       do i = 1, NPreDefName
          if(PreDefName(i)(1:klena(PreDefName(i))) .eq.
     *         comp%struc(1:klena(PreDefName(i)))) then
             comp%strucNo = i
             goto 5
          endif
       enddo
       if(comp%struc(1:4) .eq. 'new-' ) then
          read(comp%struc(5:6), *) i
          comp%strucNo = i + NPreDefName

       else
          call epseeUnderScore(comp%struc, uscl)
           tempph=  epparaphrase(comp%struc(1:uscl))
           if(tempph(1:4) .eq. 'new-' ) then
              read(tempph(5:6),  *) i
              comp%strucNo = i + NPreDefName
           else
              comp%strucNo = 0 ! not yet determined
           endif
       endif
 5     continue

!

!       if(comp.struc .eq. 'box_w'  .or.
!     *      comp.struc .eq. 'sphere_w') then
                ! trim v9.14 
       if(index(trim(comp%struc), '_w') .ne. 0) then
          Det%nworld = Det%nworld + 1
       endif
!        for the compatibility of old config data which
!        uses only small letters for media name, old
!         names are translated into canonical form
      do i = 1,  noldnames
         if(translate(1,i) .eq. comp%matter) then
            comp%matter = translate(2,i)
            goto 10
         endif
      enddo
 10   continue

      alen = klena(confdata)

      l = index(confdata, '/')
!
      l2 = index(confdata(l+1:alen), '/')
      if(l2 .eq. 0) then
         confdata(alen+1:alen+2) = ' /'
         alen = alen + 2
      endif
!           process orign, attrib and direction cos.
      call epOrigDirec
     *  (comp, nattr, ndir1, ndir2, confdata(l+1:alen))
!     
      if(comp%strucNo .gt. NPreDefName) then
!          new-i structure; memorize the number of attrib
         nattrib(comp%strucNo) = nattr  

!             count such comp. with such a structure
!!!xxx         Det.nnew(comp.strucNo - NPreDefName) =
!!!xxx     *      Det.nnew(comp.strucNo - NPreDefName) + 1
      endif
      end
!     **************************
      subroutine epOrigDirec(comp, nattr, ndir1, ndir2, data)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!
!          process origin, attrib and direction cos part
!
!   data format:
!            n  (number; 3, 3.0,  -3, -3.0 -3.0d3 etc); as it is
!            a  (defined name)  replaced by it's value
!            +      previous value + previous thickness
!            ~      previous value - current thickness
!           +n      previous value + previous thickness + n
!           ~n      previous value - current thickness  + n
!           +a      previous value + previous thickness + a
!           ~a      previous value - current thickness  + a
!           =       previous value   (=*1)
!           =n      previous value + n   
!           =a      previous value + a
!           $i  (i integer > 0)    i-th comp.'s value
!           *i                    (current - i)-th value
!           $i:n            i-th value + n 
!           $i:a            i-th value + a
!           *i:n            (current -i)-th value + n
!           *i:a            (current -i)-th value + a
!
       type(Component)::  comp  ! input/output
       integer nattr  !  input.  number of volume attributes of this 
                      !          structure
       integer ndir1, ndir2 ! input. ndir1-th to ndir2-th component of
                            !  direction cosines which specify the
                            !  volume orientation are read
       character*(*) data   ! input orgx  origy orgz  {dirx  diry...} /

       integer  klena


       integer i, j, ix, alen, icon
       logical kalpha

       character*32 orgF(3)
       integer norgF
       integer lp
       real*8 realtemp

       lp = klena(data)

       orgF(1) = ' '
       orgF(2) = ' '
       orgF(3) = ' '
       norgF = 0
!
       call kgetField(data, orgF, 3,  norgF)

!          suppose original data has " 0 0 0 xx yy zz... /"
!
!        find xxx pos in   xxx yyy zzz /
!                                    0 0 0 xxx yyy zzz /
      call kgetCpos(data, i)   !     i
      call kgetBpos(data(i:lp),  j)  !   0  0  0 xxx yyy zzz /
!                                         j
      i = i + j -1
      call kgetCpos(data(i:lp), j)   !   0  0  0 xxx yyy zzz /   
!                                           j                                 
      i = i + j - 1
      call kgetBpos(data(i:lp), j) !    0  0  0 xxx yyy zzz /   
!                                           j
      i = i + j - 1
      call kgetCpos(data(i:lp), j)  !   0  0  0   xxx yyy zzz /   
!                                             j   
      i = i + j - 1
      call kgetBpos(data(i:lp), j) !    0  0  0   xxx yyy zzz /   
!                                              j
      i = i + j - 1
      call kgetCpos(data(i:lp), j)  !   0  0  0   xxx yyy zzz /   
!                                                 j
      
      ix = i + j - 1   ! ix is xxx pos in data;
!                      or   /'s pos + 1 as in   0 0 0/
      j = index(data(1:lp), '/')    ! since / exists, j !=0
      if(j .eq. 0) then
         call cerrorMsg(data, 1)
         call cerrorMsg('above data has no /', 0)
      endif

      if(j .lt. ix) ix = j

!        if origin part has =, +, treat them specially
!       ------    origx 
      comp%offsetx = 0.      
      if(orgF(1)(1:1) .eq.  '=' ) then
         if(klena( orgF(1)) .gt. 1 ) then
            if( kalpha(orgF(1)(2:24)) ) then
!               call epexamEq(orgF(1)(2:24), comp.offsetx, icon)
               call epexamEq(orgF(1)(2:24), realtemp, icon)
               comp%offsetx = realtemp
               if(icon .ne. 0) then
                  call eperrEq(comp, orgF(1)(2:24)) ! not come back
               endif
            else
!                 read, say, 23 in   "=23".  "=/"
               read( orgF(1)(2:24), *) comp%offsetx
            endif
         endif
         comp%orgx = EqualF
      elseif(orgF(1)(1:1) .eq.  '+' ) then
         if(klena( orgF(1)) .gt. 1 ) then
            if(kalpha( orgF(1)(2:24) )) then
!               call epexamEq(orgF(1)(2:24), comp.offsetx, icon)
               call epexamEq(orgF(1)(2:24), realtemp, icon)
               comp%offsetx = realtemp
               if(icon .ne. 0) then 
                  call eperrEq(comp, orgF(1)(2:24)) ! not come back
               endif
            else
!              read, say, 23 in   "+23".  "+/"
               read( orgF(1)(2:24), *) comp%offsetx
            endif
         endif
         comp%orgx = PlusF
      elseif(orgF(1)(1:1) .eq.  '~' ) then
         if(klena( orgF(1)) .gt. 1 ) then
            if(kalpha( orgF(1)(2:24) )) then
               call epexamEq(orgF(1)(2:24), realtemp, icon)
               comp%offsetx = realtemp
               if(icon .ne. 0) then 
                  call eperrEq(comp, orgF(1)(2:24)) ! not come back
               endif
            else
!              read, say, 23 in   "~23".  "~/" "~-34" 
               read( orgF(1)(2:24), *) comp%offsetx
            endif
         endif
         comp%orgx = MinusF
      elseif(kalpha(orgF(1))) then
         call epexamEq(orgF(1), comp%orgx, icon)
         if(icon .ne. 0) then 
            call eperrEq(comp, orgF(1)) ! not come back
         endif
      else
         read(orgF(1), *) comp%orgx
      endif
!      -------    origy 
      comp%offsety = 0.
      if(orgF(2)(1:1) .eq.  '=') then
         if(klena( orgF(2)) .gt. 1 ) then
            if( kalpha(orgF(2)(2:24)) ) then
!               call epexamEq(orgF(2)(2:24), comp.offsety, icon)
               call epexamEq(orgF(2)(2:24), realtemp, icon)
               comp%offsety = realtemp
               if(icon .ne. 0) then
                  call eperrEq(comp, orgF(2)(2:24)) ! not come back
               endif
            else
               read( orgF(2)(2:24), *) comp%offsety
            endif
         endif
         comp%orgy = EqualF
      elseif(orgF(2)(1:1) .eq.  '+') then
         if(klena( orgF(2)) .gt. 1 ) then
            if(kalpha( orgF(2)(2:24) )) then
!               call epexamEq(orgF(2)(2:24), comp.offsety, icon)
               call epexamEq(orgF(2)(2:24), realtemp, icon)
               comp%offsety = realtemp
               if(icon .ne. 0) then 
                  call eperrEq(comp, orgF(2)(2:24)) ! not come back
               endif
            else
               read( orgF(2)(2:24), *) comp%offsety
            endif
         endif
         comp%orgy = PlusF
      elseif(orgF(2)(1:1) .eq.  '~' ) then
         if(klena( orgF(2)) .gt. 1 ) then
            if(kalpha( orgF(2)(2:24) )) then
               call epexamEq(orgF(2)(2:24), realtemp, icon)
               comp%offsety = realtemp
               if(icon .ne. 0) then 
                  call eperrEq(comp, orgF(2)(2:24)) ! not come back
               endif
            else
!              read, say, 23 in   "~23".  "~/" "~-34" 
               read( orgF(2)(2:24), *) comp%offsety
            endif
         endif
         comp%orgy = MinusF
      elseif(kalpha(orgF(2))) then
         call epexamEq(orgF(2), comp%orgy, icon)
         if(icon .ne. 0) then
            call eperrEq(comp, orgF(2))  ! not come back
         endif
      else
         read(orgF(2), *) comp%orgy
      endif

!       ------- origz 
      comp%offsetz = 0.
      if(orgF(3)(1:1) .eq.  '=' ) then
         if(klena( orgF(3)) .gt. 1 ) then
            if( kalpha(orgF(3)(2:24)) ) then
!               call epexamEq(orgF(3)(2:24), comp.offsetz, icon)
               call epexamEq(orgF(3)(2:24), realtemp, icon)
               comp%offsetz = realtemp
               if(icon .ne. 0) then
                  call eperrEq(comp, orgF(3)(2:24)) ! not come back
               endif
            else
               read( orgF(3)(2:24), *) comp%offsetz
            endif
         endif
         comp%orgz = EqualF
      elseif(orgF(3)(1:1) .eq. '+') then
         if(klena( orgF(3)) .gt. 1 ) then
            if(kalpha( orgF(3)(2:24) )) then
!               call epexamEq(orgF(3)(2:24), comp.offsetz, icon)
               call epexamEq(orgF(3)(2:24), realtemp, icon)
               comp%offsetz = realtemp
               if(icon .ne. 0) then 
                  call eperrEq(comp, orgF(3)(2:24)) ! not come back
               endif
            else
               read( orgF(3)(2:24), *) comp%offsetz
            endif
         endif
         comp%orgz = PlusF
      elseif(orgF(3)(1:1) .eq.  '~' ) then
         if(klena( orgF(3)) .gt. 1 ) then
            if(kalpha( orgF(3)(2:24) )) then
               call epexamEq(orgF(3)(2:24), realtemp, icon)
               comp%offsetz = realtemp
               if(icon .ne. 0) then 
                  call eperrEq(comp, orgF(3)(2:24)) ! not come back
               endif
            else
!              read, say, 23 in   "~23".  "~" "~-34" 
               read( orgF(3)(2:24), *) comp%offsetz
            endif
         endif
         comp%orgz = MinusF
      elseif(kalpha(orgF(3))) then
         call epexamEq(orgF(3), comp%orgz, icon)
         if(icon .ne. 0) then
            call eperrEq(comp, orgF(3))  ! not come back
         endif
      else
         read(orgF(3), *) comp%orgz
      endif
!  -------------------------------
      alen = klena(data)

      if(nattr .le. 0) then
         call epgetDir(comp, data(ix:alen), ndir1, ndir2)
!               read(data(ix:alen), *)
!     *        (comp.direc(i), i = ndir1, ndir2)
      else
         call eprdVolattr(comp, ix, data, nattr, ndir1, ndir2)
      endif
      end
!     *********************
      subroutine eprdVolattr(comp, ix, data,  nattr, ndir1, ndir2)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!
!       process attribut after  origin part and direction cos part
!
       type(Component)::  comp  ! input/output
       character*(*) data   ! input vol-attr {dirx  diry...} /
                ! contains org part but see next
       integer ix   ! data(ix:*) is the data containing vol-attr.
       integer nattr  !  input.  number of volume attributes of this 
                      !          structure
       integer ndir1, ndir2  ! direction cos pos.

       character*100 msg
       integer  klena, mattrb
       logical kalpha

       integer i, j, k,  icon, lenattr
!       integer       lp
       real*8 realtemp

       character*24 attrb(maxiattr)

       
!      lp = klena(data)
       do i = 1, nattr
          attrb(i) = ' '
       enddo
       mattrb = 0
!
       call kgetField(data(ix:), attrb, nattr,  mattrb)
       if( mattrb .gt. 0  .and. attrb(mattrb) .eq. '/') then
          mattrb= mattrb -1
       elseif(nattr > 0 .and. attrb(1)=='/' .and.  ! v9.157
     *     index(trim(comp%struc), '_w') /= 0) then
          mattrb = 0
       endif

!       if(mattrb .lt. nattr  .and.
!     *      ( comp.struc .eq. 'box_w' .or.
!     *       comp.struc .eq. 'shpere_w' )) then
       if(mattrb .lt. nattr  .and.    ! trim 9.14
     *      index(trim(comp%struc), '_w') .ne. 0) then
!               box_w 0 0 0/  0 0 0 type. orn
!               box_w 0 0 0/  0 0 0 / contain list  type.
!                skip reading attribute; set 0. for the frist attrib.
          Volat( comp%vol + 1 ) = 0.
          goto 10
       elseif(mattrb .lt. nattr) then
          write(msg,*)
     *    comp%cn,'-th comp. has not enough attributs'
          call cerrorMsg(msg, 0)
       endif
       
!        see if vol attr has =
       do i = 1, nattr
!          if(index(attrb(i),  '+')  .gt. 0 ) then
          if(attrb(i) .eq.  '+' ) then
             write(msg,*)
     *       comp%cn,'-th comp. has  + for vol attribute'
             call cerrorMsg(msg, 1)
             call cerrorMsg(
     *          '"=" may be used but not "+" for vol. attr',0)
          endif


         if(attrb(i)(1:1) .eq.  '=') then
            lenattr= klena( attrb(i) ) 
            if( lenattr .gt. 1) then
               if( kalpha(attrb(i)(2:lenattr) ) ) then    
                  call epexamEq(attrb(i)(2:lenattr), realtemp, icon)
                  if(icon .ne. 0) then
                     call eperrEq(comp, attrb(i)(2:lenattr)) ! not come back
                  endif
                  VolatEq( comp%vol+i ) = realtemp
               else
!                   read, say, 23 in   "=23". 
                  read(attrb(i)(2:lenattr), *) VolatEq(comp%vol+i ) 
               endif
            else
               VolatEq(comp%vol+i ) = 0.
            endif
!((((((((
            Volat( comp%vol+i ) = EqualF
!))))))))
         elseif(kalpha(attrb(i))) then
!((((((((((
            call epexamEq(attrb(i), Volat(comp%vol+i), icon)
!)))))))))
            if(icon .ne. 0) then
               call eperrEq(comp, attrb(i))
            endif
         else
!(((((((((((((
            read(attrb(i), *) Volat(comp%vol+i)
!))))))))))
         endif
       enddo
!  -------------------------------
!         skip nattr attributes in data. 
!         get position of direction cos
!         
!        (if,  for nattr=3,  possible data is
!             a b c  d1 d2 d3 ../
!                or
!             a b c /       or
!             a b c/
!        position should be next to "c"

       j = ix
       do i = 1, min(nattr, mattrb)
          call kgetCpos(data(j:), k)   !   aaa  bbb  cccc
                                         !   k
          j = k + j -1                  !    j 
          call kgetBpos(data(j:), k) !  aaa  bbb ccc
                                       !     k
          j = k + j -1                 !     j
       enddo
!          aaa bbb ccc
!                     j   
       call epgetDir(comp, data(j:), ndir1, ndir2)
!       read(data(j:), *)
!     *        (comp.direc(i), i = ndir1, ndir2)
 10    continue
!         This is to replace "=" by it's corresponding
!          value. If this is delayed, newvol could not
!          be proccessed because usually, attribute values
!          are checked just after this routine is called;
!         The EqualF becomes the invalid value. So we must
!         replace "=" here. But, we don't know the 
!         previous value to be used as "=".  Therefore
!         common values are directly used  in the next
!         subroutine. This is not a recommended way.
!         Danger !  If the attribute values are not
!         checked inside the newvolume data processing
!         routine, we can delay the replacement at
!         epAttrb.

       call epAttrb1
       end
!     ****************************  
      subroutine epgetDir(comp, data, ndir1, ndir2)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
       type(Component)::  comp   ! input/output
      character*(*) data        ! input dirx  diry../
      integer ndir1, ndir2
      character*24 dirf(9)
      integer nf, need, i, icon
      character*128 msg
      logical kalpha


      call kgetField(data, dirf, 10, nf)
      
      need = ndir2 - ndir1 + 1  ! no. of field 

      if(nf .eq. 1) then
!          nothing to do; no direction cos.
      elseif( dirf(1) .eq. '/') then
!           nothing to do; no direction cos.
      elseif(nf .ge. need) then
         do i = 1, need
            if(kalpha(dirf(i))) then
               call epexamEq(dirf(i), comp%direc(ndir1+i-1), icon)
               if(icon .ne. 0) then
                  call eperrEq(comp, dirf(i))
               endif
            else
               read(dirf(i), *) comp%direc(ndir1+i-1)
            endif
         enddo
      elseif(nf .ge. need - 3) then
!           2nd group is assumed to be default (0 1 0)
!           so that read only first dir.
         do i = 1, need -3
            if(kalpha(dirf(i))) then
               call epexamEq(dirf(i), comp%direc(ndir1+i-1), icon)
               if(icon .ne. 0) then
                  call eperrEq(comp, dirf(i))
               endif
            else
               read(dirf(i), *) comp%direc(ndir1+i-1)
            endif
         enddo
      else
         call cerrorMsg(data, 1)
         call cerrorMsg(
     *   'above data has not enough direction cosines',1)
         write(msg,*) 'for ',comp%cn, '-th component'
         call cerrorMsg(msg, 1)
         call cerrorMsg('OR: some attributes are missing',0)
      endif
      end
!     ********************
      subroutine epcontain(fatal)
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"

           integer klena, ldat
!            see if 2nd /           
       integer l1, l2, i
       logical fatal
       character*2 bq/'\\'/  ! solaris cannot use '\'


!            get actual length of data
       ldat = klena(confdata)
!        find  first '/'  such as in 
!         1 box pb  0 0 0/
       l1=index(confdata, '/')
!        find 2nd /
       l2 = index(confdata(l1+1: ldat) , '/')
!          clear matreska counter for nct-th comp.
! (((((((((((((
!       do i = 1, MaxMatreska
!           Det.cmp(Det.nct).Contains(i) = 0
!           Det.cmp(Det.nct).ContainsR(i) = 0
!       enddo
           Det%cmp(Det%nct)%Contains =-1
           Det%cmp(Det%nct)%ContainsR =-1
! ))))))))))))))
       if(l2 .gt. 0) then
!           there is /
          if(confdata(ldat:ldat) .eq. bq(1:1)) then
!                   there is cont lines. not supported yet
             call cerrorMsg(
     *        'continuation by \ not supported yet', 0) 
          endif
!              put / in the last pos.
          confdata(ldat+1:ldat+1) = '/'
! ((((((((((
!          read(confdata(l1+l2+1:ldat +1), *)
!     *    (Det.cmp(Det.nct).ContainsR(i), i=1, MaxMatreska)
          read(confdata(l1+l2+1:ldat +1), *)
     *    ( CnArea(i+CnCounter), i=1, maxCnArea-CnCounter)
! ))))))))))))))))))
       endif

!   #        NMatreska ContainsR NContainer  PContained NPContainer
!   1              0      0         0           0          0                  
!   2  / 3 4       2      3 4       0           12         1
!   3  / 5 6 7     3     5 6 7      1           12         1
!   4              0      0         1           12         1
!   5              0      0         1           12         1
!   6              0      0         1           12         1 
!   7              0      0         1           12         1
!   8  / -9 -10    2    -9 -10      0           0          0
!   9              0      0         1           8          1
!  10              0      0         2          8 11        2
!  11  / -10       1     -10        0           0          0           
!  12  / -2        1      -2        0           0          0   
!                    Contains=|ContainsR|

       call epcontain1(Det%nct,  fatal)
       end
       subroutine epcontain1(n, fatal)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer n  ! input.  n-th component is  being examined
                 !    if  there is contain/pcontain  relation

      logical fatal
      integer i, j, k
      logical:: same
      character*100  msg

!          count matreska

       Det%cmp(n)%NMatreska = 0
!  (((((((((((
       do i = 1, maxCnArea - CnCounter
          if(CnArea( CnCounter+i ) .eq. 0) then
             goto 10
          endif
!!!!!!!!////  !  don't count the one already listed  V9.157<<<<<<<<<
          if( i == 1 ) then
             Det%cmp(n)%NMatreska = Det%cmp(n)%NMatreska  + 1
          else
             same = .false.
             do j = 1,  Det%cmp(n)%NMatreska 
                if( CnArea( CnCounter+ i) == 
     *                   CnArea( CnCounter+ j) ) then
                   same =.true.
                   exit
                endif
             enddo
             if(.not. same) then
                Det%cmp(n)%NMatreska = Det%cmp(n)%NMatreska  + 1
                if( i /=  Det%cmp(n)%NMatreska ) then
                   CnArea( CnCounter +  Det%cmp(n)%NMatreska ) =
     *                  CnArea( CnCounter+ i) 
                endif
             endif
          endif
!!!!!!!!!!!!! >>>>>>>>>  v9.157
       enddo
 10    continue
       if(Det%cmp(n)%NMatreska .gt. 0) then
          Det%cmp(n)%ContainsR = CnCounter
       endif
!!         next is needed to process partial containing
       CnCounter = CnCounter + Det%cmp(n)%NMatreska 
!  ))))))))))))))
       do i = 1, Det%cmp(n)%NMatreska
          if(i .eq. 1) then
             Det%cmp(n)%Contains = CnCounter
          endif
!          Det.cmp(n).Contains(i) =abs( Det.cmp(n).ContainsR(i))

          call epCnArea( abs(CnArea( i+ Det%cmp(n)%ContainsR)) )
          if(CnArea(CnCounter) .eq. n) then
!             cannot eat tail
             fatal  = .true.
             write(msg, *)
     *            ' ***** comp. ', n, ' is eating self'
             call cerrorMsg(msg, 1)
          endif
       enddo
!
!          next is to init cmp(j)'s PContained
!
       do i = 1, Det%cmp(n)%NMatreska
          if( CnArea(i+Det%cmp(n)%ContainsR) .lt. 0 ) then
!                  partially contained comp.
             j =  CnArea( i+Det%cmp(n)%Contains )
             Det%cmp(j)%PContained =  -1
          endif
       enddo


       do i = 1, Det%cmp(n)%NMatreska
          if( CnArea(i+Det%cmp(n)%ContainsR) .lt. 0 ) then
!                  partially contained comp.
!((((((((((((
!                j =  Det.cmp(n).Contains(i)
             j =  CnArea( i+Det%cmp(n)%Contains )
             if(Det%cmp(j)%PContained .eq. -1) then
                Det%cmp(j)%PContained = CnCounter
             endif
!))))))))))
!                  count no. of comp. which partially
!                  contains j
             Det%cmp(j)%NPContainer = Det%cmp(j)%NPContainer + 1
!                  j's  parent matreska is Det.nct
             k = Det%cmp(j)%NPContainer
             call epCnArea(n)
!              Det.cmp(j).PContained(k)
!     *           = n
!              To  identify inner comp. and outer comp.
!              Later,  NContainer(i)  = 0 meanas
!              that the  i-th comp. is not contained in another one.
           
!             j = Det.cmp(n).Contains(i)
             j = CnArea( Det%cmp(n)%Contains + i)
             Det%cmp(j)%NContainer = Det%cmp(j)%NContainer + 1
          endif
       enddo
       end

      subroutine epqcyl(rcyl)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!        returns radius of the cylinder 
!        when all comp. is cylinder of the
!        same radius.
      real*8  rcyl, rpipi, rpipo
!((((((((((
      rcyl = Volat( Det%cmp(1)%vol+cylr)
!)))))))
      return
!     **************
      entry  epqpip(rpipi, rpipo)
!         returns inner and outer    
!         radius of the pipe when
!         all comp. are the same pipe
!(((((((((((((
      rpipi = Volat( Det%cmp(1)%vol+pipeir)
      rpipo = Volat( Det%cmp(1)%vol+pipeor)
!))))))))))))
      end
!      character*200 x,  y
!       x = 
!     * ' 4 box sp 0 0 0 / .0000   .0000 .000000 / .00 .0 .0000000 '//
!     *   ' 0.000 0.0000  .10000 15.0000 1.330000 .0000 / 1 2 3'
!      call epcompZero(x,y)
!      write(*, *) y
!      end
      subroutine epcompZero(x, y)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
      character*(*) x, y
!          0.0000, .0000  1.0000 etc is
!        compressed as 0  0  1. etc.
!        but if it dose not have "."
!        such action is diabled
      integer maxp
      parameter (maxp = 16)
      character*20 pat(maxp)
      character*30 f(maxiattr)
      integer klena, i,j, nf, k

      data pat/
     * '0 ',
     * '00 ',
     * '000 ',
     * '0000 ',
     * '00000 ',
     * '000000 ',
     * '0000000 ',
     * '00000000 ',
     * '000000000 ',
     * '0000000000 ',
     * '00000000000 ',
     * '000000000000 ',
     * '0000000000000 ',
     * '00000000000000 ', 
     * '000000000000000 ',
     * '0000000000000000 '/

      call kgetField(x, f, maxiattr, nf)

      do i = 1, nf
         if(index( f(i), ".") .eq. 0) goto 10
         do j =  maxp, 1, -1
            k = index(f(i), pat(j)(1:(klena(pat(j))+1)) )
            if(k .gt. 0) then
               f(i)(k:klena(f(i))) = ' '
               if(f(i) .eq. ' ') then
                  f(i) = '0'
               elseif(f(i) .eq. '-.') then
                  f(i) = '0'
               else
                  k = klena(f(i))
                  if(f(i)(k:k) .eq. '.') then
                     f(i)(k:k) = ' '
                     if(f(i) .eq. ' ') f(i) = '0'
                  endif
                  goto 10
               endif
            endif
         enddo               
 10      continue
      enddo
      k = 0
      y = ' '
      do i = 1, nf
         j =klena(f(i))
         k = k + 1
         y(k:k+j) = ' '// f(i)(1:j)

         k = k + j 
      enddo
      end
!     ***************************
      subroutine epgetMedia
      implicit none
!        by looking at the configuration data, examine
!     a file for a given matter is ready or not. Then
!     if the file is ready, read it.  Count  number of
!     different media.
!
#include "ZepTrackv.h"
#include "Zcnfig.h"
#include "ZepManager.h"
      

      character*30 msg
      logical fatal
      integer i, j, icon

      character*8  missingMedia(Maxmedia)
      integer  missing, k

      character*24  tempc

      missing = 0
      MediaNo = 0

      fatal = .false.

      do i = 1,  Det%nct
!!!!!!     v9.14
         tempc = Det%cmp(i)%matter
         if( tempc == "world" ) exit
         if( tempc =="sp2" ) tempc="sp"  
!!!!
         call epSeeIfAlias(tempc)  ! if alias, replace it by true
         do j = 1,  MediaNo
!            if(Det.cmp(i).matter .eq. Media(j).name) then
            if(tempc .eq. Media(j)%name) then  ! v9.14
!               aleady read
!                i-th component has media index j.
               Det%Cn2media(i) = j
               goto 10
            endif
         enddo
!            this is new matter
         if(MediaNo .lt. Maxmedia) then
            MediaNo = MediaNo + 1
            j = MediaNo
         else
            write(msg, *) Maxmedia
            call cerrorMsg(
     *       'too many different media > '//msg, 0)
         endif
!            see if the matter file exists, and if exist, read it
!!!!!        tempc = Det.cmp(i).matter  v9.14
!            if we use Det.cmp(i).matter directly in below,
!           the string length cannot be transmitted correctly.
!           so that we replace it  by tempc which can  hold
!           upto 24 characters. 
         call eprdMFile(tempc, icon) 

         if(icon .ne. 0) then
!            see if already appeared
            do k = 1, missing
!!!               if(missingMedia(k) .eq. Det.cmp(i).matter) then
               if(missingMedia(k) .eq. tempc ) then  ! v9.14
                  goto 5
               endif
            enddo
            missing = missing + 1
            if(missing .le. Maxmedia) then
               write(msg, *) 'Media file:',
!!     *          Det.cmp(i).matter, ' not exist'
     *          tempc, ' not exist'   ! v9.14
               call cerrorMsg(msg, 1)
!!               missingMedia(missing) = Det.cmp(i).matter
               missingMedia(missing) = tempc   !  v9.14
            else
               call cerrorMsg(
     *              'Seems too many missing media names',1)
               missing = Maxmedia
            endif
 5          continue
            
            MediaNo = MediaNo - 1
            fatal = .true.
         else
!                i-th component has media index j.
            Det%Cn2media(i) = j
!            call epLen2Emin(Media(j))   ! make table for Length vs E
!                     E is minimum energy that can run Lenth
         endif
 10      continue
      enddo

      if(fatal) then
         call cerrorMsg(
     *    'Media files shown above are missing in any of the'//
     *    ' following directories: ',  1)
         do i = 1, MaxMediaDir
            if(MediaDir(i) .ne. ' ') then
               call cerrorMsg(MediaDir(i), 1)
            endif
         enddo
         call cerrorMsg('Fatal error, good bye', 0)
      endif

      NoOfMedia = MediaNo   ! memorize the # of media
      call epRangeAlloc( NoOfMedia )  ! allocate memory for Range calc.
      do i = 1, NoOfMedia
         call epRangeMkTbl(Media(i), i)   ! make Range tbl for i-th media
      enddo
 
      end
!     ***************************
      subroutine epwriteMedia
      use modAlias
      implicit none
!          print media info
#include "ZepTrackv.h"
#include "Zcnfig.h"
#include "ZepManager.h"

      integer:: i
      write(0,*) 'list of meida used'
      do i = 1, NoOfMedia
         write(0, *) Media(i)%name
      enddo

      if( NoOfAlias > 0 ) then
         write(0,*) 'alias  of  media' 
      endif
      do i = 1, NoOfAlias
         write(0,*) aliasName(i),"   ", trueName(i)
      enddo
      end       subroutine epwriteMedia
      
!     *******************************
      subroutine epqncp(numberOfComp)
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
        integer numberOfComp  ! output.  No of Component defined.

        numberOfComp = Det%nct
      end
!     ************************
      subroutine epqmat(n, mat)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer n
      character*(*) mat  ! output
      call epchkcmpn('epqmat', n)
      mat = Det%cmp(n)%matter
      end
!        Next is usable only within subdetector procssing time
!      after reading all config data, use the one in epquerry.f
!      subroutine epqSubdName(n, name)  
!      implicit none
!#include "Zep3Vec.h"
!#include  "Zcnfig.h"
!      integer,intent(in)::n ! component number
!      character(len=16),intent(out)::name
!
!      integer::epIsSubD
!      integer j
!      j = epIsSubD(n)
!      if(j > 0 ) then
!         name = SubDName(j)
!      else
!         name = ' '
!      endif
!      end

      subroutine epqmatrhoc(n, mat, lc, rhoc)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer n
!
!          This returns  Air*0.9345 etc as mat if
!          rhoc != 1.0  else  simply Air etc.
!
      character*(*) mat  ! output  length must be>15
      integer lc  !  output. length of 'mat' content
      real rhoc   !  rhoc. 1. or 0.9356 etc.

      character*20 matterAndrhoc
      integer klena
      call epchkcmpn('epqmatrhoc', n)       
      rhoc = Det%cmp(n)%rhoc
      if( rhoc .eq. 1.) then
         matterAndrhoc =Det%cmp(n)%matter
      else
         if(rhoc .ge. 0.01 .and. rhoc .lt. 10.0) then
!                                      9.99999 0.01001
            write(matterAndrhoc, "(a, a, f11.5)")
     *        Det%cmp(n)%matter(1:klena(Det%cmp(n)%matter)),
     *        "*",  Det%cmp(n)%rhoc 
         else
            write(matterAndrhoc, "(a, a, 1p, g12.4)")
     *        Det%cmp(n)%matter(1:klena(Det%cmp(n)%matter)),
     *        "*",  Det%cmp(n)%rhoc 
         endif
         call kseblk(matterAndrhoc, "&", lc)
      endif
!             compress if 2.030000E-5 etc to 2.03E-5
      call epsup00E(matterAndrhoc)
      call epsuplast0(matterAndrhoc)
      lc = klena(matterAndrhoc)
      mat = matterAndrhoc(1:lc)
      end
!     ************************
      subroutine epqstruc(n, struc)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer n
      character*(*) struc  ! output
      call epchkcmpn('epqstruc', n)
      struc = Det%cmp(n)%struc
      end

!     ************************
      subroutine epqCount(n, countn, counte)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer n
      integer countn, counte
      call epchkcmpn('epqCount', n)
      countn = Det%cmp(n)%CountIO
      counte = Det%cmp(n)%CountDE
      end
!     ***********************
      subroutine epqorg(n, origin)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

      integer n  ! input. component number
       type(epPos)::  origin  ! output. n-th comp.'s origin

!         inqurire the origin of n-th comp.

      call epchkcmpn('epqorg', n)
      origin%x= Det%cmp(n)%orgx
      origin%y= Det%cmp(n)%orgy
      origin%z= Det%cmp(n)%orgz
      end
!     ***************************
      subroutine epchkcmpn(id, n)
!        check n is valid comp. # or not
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

      character*(*) id  ! input. id to be printed
      integer n      !  input  comp. #

      character*50 msg

      if(n .le. 0 .or. n .gt. Det%nct) then
         write(msg, *)
     *   'comp.#=', n,' for ', id, ' is out of range'
         call cerrorMsg(msg, 0)
      endif
      end
!     *****************************
      subroutine epqvolatr(n, na, vol)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
!          inquire vol. attr. of n-th comp.
      integer n ! input. comp. number
      integer na ! output. # of attrib.
      real*8 vol(*)  ! output. at least na-dim.

      integer loc(maxiattr)

      integer i, j

      call epchkcmpn('epqvolatr', n)
      na = Det%cmp(n)%Nattributes

      call epatloc( Det%cmp(n), loc)
      do i = 1, na
!((((((((((((
         j = loc(i)
         vol(i) = Volat( Det%cmp(n)%vol+j )
!)))))))))
      enddo
      end
!     *****************************
      subroutine epqcmpdircos(n, dir)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
!          inquire vol. attr. of n-th comp.
      integer n ! input. comp. number
      real*8 dir(9)  ! output.  direction cos.

      integer i

      call epchkcmpn('epqcmpdircos', n)
      do i = 1, 9
         dir(i) = Det%cmp(n)%direc(i)
      enddo
      end

      subroutine epCnArea(cn)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
      integer cn  ! input.  component number to be stored.
!           CnCounter is never cleared except at the very first.
!       Each subdeterctor component is remembered with its  CnCounter
!       memorized 
      character*80 msg

      if(CnCounter .lt. maxCnArea) then
         CnCounter = CnCounter + 1
         CnArea(CnCounter) = cn
      else
         write(msg, *) 'MAX_CN_AREA in ZepMaxDef%h=',maxCnArea,
     *     ' exhausted '
         call cerrorMsg(msg, 0)
      endif
      end
      subroutine epAttribWatch(n, vol)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
      integer n  ! input.  number of attributes to be used for
                 !    a component to be treated
      integer vol ! output.   Volat(vol+i) can be used for the
                  !           i-th attribute of the component

      
      if(AttrCounter + n  .le.  maxattr ) then
         vol =  AttrCounter
         AttrCounter =  AttrCounter + n 
      else
         call cerrorMsg(
     *   'too many attributes; enlarge MAX_ATTRIB in ZepMaxdef%h', 0)
      endif
      end

      subroutine epnormvec(vec)
      implicit none
#include "Zep3Vec.h"
      real*8 norm

       type(ep3Vec)::  vec

      norm = sqrt(vec%x**2 + vec%y**2 + vec%z**2)
      if(norm .eq. 0.) then
         vec%x = 0.
         vec%y = 0.
         vec%z = 1.
      else
         vec%x = vec%x/norm
         vec%y = vec%y/norm
         vec%z = vec%z/norm
      endif
      end
      subroutine epGetThick(n, vec, origin )
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

      integer n    ! input.  comp. number
       type(ep3Vec)::  vec    !   output ; abc  of the
                             !   enveloping box in world coord.
      real(8),intent(out):: origin(3)


       type(epPos)::  posw, maxpos, org, orgx
       type(ep3Vec)::  abc

      integer i, imx
      real*8 dx(8), dy(8), dz(8)
      data dx/0., 1., 0., 1., 0., 1., 0., 1./
      data dy/0., 0., 1., 1., 0., 0., 1., 1./
      data dz/0., 0., 0., 0., 1., 1., 1., 1./

      call epqenvlper(n, org, abc)
      if(NVTX .eq. 0) then
         imx = 8
      else
         imx = NVTX
      endif

      do i = 1, imx
         if(NVTX .eq. 0) then
!             for all 8 corners of the box
            posw%x = org%x + dx(i) * abc%x
            posw%y = org%y + dy(i) * abc%y
            posw%z = org%z + dz(i) * abc%z
         else
            posw%x = org%x +  VTXx(i) 
            posw%y = org%y +  VTXy(i)
            posw%z = org%z +  VTXz(i)
         endif
!             now posw is one of imx corners in world coord.
         if(i .eq. 1) then
            orgx = org
            maxpos = posw
         else
            orgx%x = min(orgx%x, posw%x)
            orgx%y = min(orgx%y, posw%y)
            orgx%z = min(orgx%z, posw%z)
!             get max pos in world coord.
            maxpos%x = max(maxpos%x, posw%x)
            maxpos%y = max(maxpos%y, posw%y)
            maxpos%z = max(maxpos%z, posw%z)
         endif
      enddo
      vec%x = maxpos%x - orgx%x
      vec%y = maxpos%y - orgx%y
      vec%z = maxpos%z - orgx%z
      origin(1) =  orgx%x
      origin(2) =  orgx%y
      origin(3) =  orgx%z

      end
      subroutine epseeUnderScore(sname, loc)
      implicit none
      character(*),intent(in)::sname  ! structure name
      integer,intent(out):: loc ! if _ is used, loc+1  
                    ! is the _ location in sname
                    ! else length of sname 
      loc = index(sname,"_")
      if(loc > 0 ) then
         if(loc < 2 ) then
            write(0,*)
     *      ' component name error: wrong use of _ '
            write(0,*) ' it is ', sname
            stop
         endif
         if( len( trim( sname(loc+2:) ) ) > 2) then
            write(0,*)
     *      ' component name error: wrong use of _ '
            write(0,*) ' it is ', sname
            stop
         endif
         select case(trim( sname(loc+1:) ) )
         case ("x")
         case ("y")
         case ("z")
         case ("xy")
         case ("xz")
         case ("yz")
         case ("yx")
         case ("zx")
         case ("zy")
         case("w")   ! v9.164
         case default
            write(0,*)
     *      ' component shape error: wrong use of _ '
            write(0,*) ' it is ', sname
            call cerrorMsg('back trace',0)
         end select
         loc = loc -1
      else
         loc = len(sname)
      endif
      end
