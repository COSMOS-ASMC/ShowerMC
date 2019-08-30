      implicit none
      real(8)::E0n, xs
      real(8),external::SECTNU
      integer::iap, iat
      call qgs01init
      write(0,*) 'Enter Ap(2~64), At(1~64)'
      read(*,*) iap, iat
      
      E0n=10.
      do while(E0n<2.e11)
         xs = SECTNU(E0n, iap, iat)
         write(*,*)  E0n, xs, iap, iat
         E0n= E0n*10.**0.1
      enddo
      end
      subroutine cqQGS1File( file1, file2)
      character*(*) file1, file2
      file1='../QGSDAT01'
      file2='../SECTNU'
      end
      real*8 function  QSRAN(X)
      real*8  X                 !  not used                                     
      real*8 u
      call rndc(u)
      QSRAN = u
      end
      subroutine qgs01init
      implicit none
      integer MONIOU
      COMMON /AREA43/ MONIOU
            
      INTEGER          NSP
      COMMON /AREA12/  NSP
      INTEGER          ICH(95000)
      DOUBLE PRECISION ESP(4,95000)
      COMMON /AREA14/  ESP,ICH
      MONIOU = 0

      call PSASETC
      call XXASET
      call PSAINI
      end
