!     ****************************************************************
!     *                                                              *
!     * bremr:  used to solve equation for making sampling table     *
!     *         for brems                                            *
!     *                                                              *
!     * pairr: //  for pair                                          *
!     *                                                              *
!     ****************************************************************
!
!
      real*8 function bremr(v)
      implicit none
      real*8 v
!
      real*8 upsi, emain
      common /cbrem/upsi, emain

      real*8 ans

      call totcb(v, 1.d0, ans)
      bremr=ans/upsi-1.
      end
!
!     ***********
      real*8 function pairr(v)
      implicit none
      real*8 v
!
      real*8 upsi, emain
      common /cbrem/upsi, emain

      real*8 ans

      call totcp(v , 1.d0, ans)
      pairr=ans/upsi - 1.
      end
!     *****************
      subroutine totcb(vmin,vmax,ans)
!     *****************
      implicit none
      real*8 vmin, vmax, ans
!
!        integration of bremsung function from vmin to vmax.
!
      real*8 fbrem
      external fbrem
      real*8 v2, v1, ans1
!
      ans=0.
      v2=vmax

      do while(.true.)
         v1=v2/20.
         v1= max(v1, vmin)
         call k16pGaussLeg(fbrem, v1, v2, 16,  ans1)
         ans=ans+ans1
         if(v1 .eq. vmin) goto 5
         v2=v1
      enddo
    5 continue
      end
!
!
!
!     ***********
      subroutine totcp(vmin,vmax,ans)
!     ***********
      implicit none
      real*8 vmin, vmax, ans
!
!
!        integralation of pair-cre function
!
!
      real*8 d, vt, fpair, ans1, ans2
      external fpair  
!
      d=(vmax-vmin)/30.
      vt=vmax-d
      call k16pGaussLeg(fpair, vmin, vt, 10, ans1)
      call k16pGaussLeg(fpair, vt, vmax, 10, ans2)
      ans=ans1+ans2
      end
