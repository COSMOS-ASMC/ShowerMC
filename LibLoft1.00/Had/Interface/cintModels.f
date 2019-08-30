!      *******************
      subroutine cintModels(from)
!
!           establish the energy bound for the given interaction model.
!          and corresponding init.
      implicit none
#include  "Zmanagerp.h"
#include  "Zmanager.h"
#include  "Zevhnp.h"
#include  "Zptcl.h"
#include  "Zevhnv.h"
      
      character*(*)  from  ! input. cosmos or epics ; from which 

      integer i, j
      character*100 path

      CosOrEpi = from
      call cgetTopDir   ! get LibLoft top dir and set it in Zmanager.h

      do i = 1, MaxIntMdls
         ModelList(i) = ' '
         InteErg(i) = 2.d12  ! GeV 10^21 eV
         ModelList2(i) = ' '
         InteErg2(i) = 2.d12  ! GeV 10^21 eV
      enddo
!             ModelList       InteErg
!              "xxx"            100.
!              "yyy"            1.d12      ====>  xxx if E< 100  else yyy
!
!              "xxx"            1.d12      ====>  xxx at entier energy.
!
      IntModel(64:64) = '/'
      if( XsecModel == " ") then
         XsecModel = IntModel
      else
         XsecModel(64:64) = '/'
      endif
 
      read(IntModel, *, Err=100)
     *    (ModelList(i), InteErg(i), i=1, MaxIntMdls)
      read(XsecModel, *, Err=100)
     *    (ModelList2(i), InteErg2(i), i=1, MaxIntMdls)

      do i = 1, MaxIntMdls
         if(ModelList(i) .eq. ' ') then
            NoOfMdls = i-1
            goto 10
         endif
      enddo
      NoOfMdls = MaxIntMdls
 10   continue
      do i = 1, MaxIntMdls
         if(ModelList2(i) .eq. ' ') then
            NoOfMdls2 = i-1
            goto 20
         endif
      enddo
      NoOfMdls2 = MaxIntMdls
 20   continue

      do i = 1, NoOfMdls-1
!         if(InteErg(i) .ge. InteErg(i+1)) then
         if(InteErg(i) > InteErg(i+1)) then  ! v7.633
!                 energy region invalid
            call cerrorMsg(
     *           'IntModel is invalid; energy region not ascending', 1)
            call cerrorMsg(IntModel, 0)
         elseif(InteErg(i) == InteErg(i+1)) then ! v7.633
!             a model written between(i) and (i+1)
!             will not be used explicitly but initialization
!             is tried. It might be used as a rescuee when a
!             some model cannot cope with a specific particle.
            ! issue some message
            write(0,*) '================================='
            write(0,*)
     *        i,'-th and ',i+1,'-th energy in IntModel are same' 
            write(0,*) 'The model specified is ',ModelList(i+1), '.'
            write(0,*) 'It will not be used explicitly but might be'
            write(0,*) 'used when some model needs help for a '
            write(0,*) 'specific particle. dpmjet3 is a good  ' 
            write(0,*) 'candidate of such a model'
            write(0,*) '================================='
         endif

      enddo
      if( NoOfMdls2 /=  NoOfMdls .or.
     &    any(InteErg(1:NoOfMdls) /= InteErg2(1:NoOfMdls)) ) then
         do i = 1, NoOfMdls2-1
            if(InteErg2(i) > InteErg2(i+1)) then  
!                 energy region invalid
               call cerrorMsg(
     *           'XsecModel is invalid; energy region not ascending', 1)
               call cerrorMsg(XsecModel, 0)
            elseif(InteErg2(i) == InteErg2(i+1)) then
             write(0,*) '================================='
             write(0,*)
     *        i,'-th and ',i+1,'-th energy in IntModel are same' 
             write(0,*) 'The model specified is ',ModelList2(i+1), '.'
             write(0,*) 'It will not be used explicitly but might be'
             write(0,*) 'used when some model needs help for a '
             write(0,*) 'specific particle. dpmjet3 is a good  ' 
             write(0,*) 'candidate of such a model'
             write(0,*) '================================='
          endif
       enddo
      endif
!         exam models
      do i = 1, NoOfMdls
         do j = 1, nmdls
            if( ModelList(i) .eq.  RegMdls(j) ) goto 25
         enddo
         call cerrorMsg( ModelList(i),  1)
         call cerrorMsg('above model is not yet registered', 0)
 25      continue
      enddo
!         exam models
      do i = 1, NoOfMdls2
         do j = 1, nmdls
            if( ModelList2(i) .eq.  RegMdls(j) ) goto 28
         enddo
         call cerrorMsg( ModelList2(i),  1)
         call cerrorMsg('above model is not yet registered', 0)
 28      continue
      enddo

!
      if( index(IntModel, 'fritiof1.6') .gt. 0 .or.
     *    index(IntModel, 'nucrin') .gt. 0  .or.
     *    index(IntModel, 'dpmjet3') .gt. 0  .or.
     *    index(IntModel, 'incdpm3') .gt. 0  .or.
     *    index(IntModel, 'jam') .gt. 0  ) then

!           for Lund init.  dpm may use Lund at low energy or
!                                    for some particles
         call haddenC
         call chanwnC
      endif
      if( index(IntModel, 'phits') .gt. 0  .or.
     *    index(XsecModel, 'phits') .gt. 0 ) then
         call cprePhits
      endif
      if( index(IntModel, 'qgsjet2') .gt. 0 .or.
     *    index(XsecModel, 'qgsjet2') .gt. 0 ) then
         call ciniQGS
      endif
      if( index(IntModel, 'epos') .gt. 0 .or.
     *    index(XsecModel, 'epos') .gt. 0 ) then
         call ceposIniAll
      endif
      if( index(IntModel, 'sibyll') > 0 .or.
     *    index(XsecModel, 'sibyll') > 0 ) then
         call csibyllinit
      endif
      if( index(IntModel, 'incdpm3') .gt. 0 ) then
         call ciniincdpm3
      endif

      if( index(IntModel, 'gheisha' ) .gt. 0) then
!         Gheisha; init. Gheisha cannot use hadrin/nucrin in Cosmos
         call gpart
         call gheini
      endif
      if( index(IntModel, 'jam' ) .gt. 0) then
         ! nothing to do ?
      endif
      if( index(IntModel, 'dpmjet3') .gt. 0 .or.
     *    index(XsecModel, 'dpmjet3') .gt. 0  .or.
     *    index(IntModel, 'incdpm3') .gt. 0  ) then
         if(from .eq. 'cosmos') then
            ! init for dpmjet.  for Cosmos
            ! If DpmFile is ' ',
            ! $LIBLOFT/Data/DPM/atmos.inp is the control card for dpmjet3.
            If( DpmFile  .eq.  ' ' ) then
               call cformFullPath('atmos.inp', path)
               call cinidpmjet(path)
            else
               call cinidpmjet(DpmFile)
            endif
         elseif(from .eq. 'check') then
!           it is assumed that the Glauber file is
!           specified by prefix/...GLB

         elseif(from == 'epics' .or. from == 'gencol') then
!              dpmjet.inp is assumed in the same directory
!              as config file.
            call cformFullPath('dpmjet.inp', path)
            call cinidpmjet(path)
         else
            call cerrorMsg(from, 1)
            call cerrorMsg('above "from" in cintModels invalid',0)
         endif
      endif
      return

 100  continue
      call cerrorMsg(
     * 'IntModel syntax error; prob.missing " mark'//
     * ' for string; IntModel is', 1)
      call cerrorMsg(IntModel, 0)
      end

!           
!          next is for dpmjet file management.
!          dpmfiles for Epcis is assumed to be in the
!          same directory as config file resides so
!          extract the directory where the config file is.
      subroutine cfixPrefix(dsn)
      implicit none
#include "Zmanager.h"
      character*(*) dsn   ! input. config file name or path
!            output. PrefixConf= './' or '/.../.../' where
!                              config file exist.
!                    PrefixLeng is its length. 
!      
!
      integer i, j, klena, k

      if( index(dsn, '/') .eq. 0 ) then
         PrefixConf = './'
         PrefixLeng = 2
      else
!             find last '/'
         i = klena(dsn)
         do j = i, 1, -1
            k = index(dsn(j:j), '/')
            if( k .gt. 0 ) then
               PrefixConf = dsn(1:j)
               PrefixLeng = j
               goto 33
            endif
         enddo
      endif
 33   continue
      end
!     ***********  get Cosmos Top directory
!        ==>              LIBLOFT top directory 
!             and set it in TopDir in Zmanager.h
      subroutine cgetTopDir
      implicit none
#include "Zmanager.h"


      character*1 NULL
      integer  kgetenv

      NULL = char(0)
      TopDir = ' '
!     TopDirLeng = kgetenv("COSMOSTOP"//NULL, TopDir)
      TopDirLeng = kgetenv("LIBLOFT"//NULL, TopDir)      
      end

      subroutine csetCosOrEpi(from)
      implicit none
#include "Zmanager.h"
      character*(*) from   !  input. cosmos or epics
       CosOrEpi = from
       end

      subroutine cformFullPath(file, path)
      implicit none
#include "Zmanager.h"
      character*(*) file   !  input. file name
      character*(*) path   ! output. $COSMOSTOP/Data/DPM/file  
      integer klena

      path = ' '

      if(CosOrEpi .eq. 'cosmos') then
         path = TopDir(1:TopDirLeng)//'/Data/DPM/'
     *     //file(1:klena(file))
      elseif(CosOrEpi == 'epics' .or. CosOrEpi == 'check'
     * .or. CosOrEpi == 'gencol') then
         path = PrefixConf(1:PrefixLeng)//file(1:klena(file))
      else
         call cerrorMsg(
     *   "cintModels('cosmos') or cintModels('epics')"//
     *   " must have been called ", 0)
      endif

      end
      subroutine cputGencolCMS(cms)
!           this is introduced for EPOS
!     to save the single precision boost in EPOS
!     call afinal is avoided to convert cms to lab
!     conversion of which result becomes not accurate 
!     if we convet it to cms again. (Gencol needs to do so)
!
      implicit none
#include "Zptcl.h"
      type(ptcl):: cms ! input. Gencol must
                      ! inform the cms (/n)

      common /forGencol/ cmsSave
      type(ptcl):: cmsSave
      cmsSave = cms
      end
      subroutine cqGencolCMS(cms)
!      ask gencol cms
!
      implicit none
#include "Zptcl.h"
      type(ptcl):: cms  ! outut. 
      common /forGencol/ cmssave
      type(ptcl):: cmsSave
      cms = cmsSave 
      end
