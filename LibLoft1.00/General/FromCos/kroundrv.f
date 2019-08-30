!     ****************************************************************
!     *                                                              *
!     * kroundrv: get a round value of a given real value
!     *                                                              *
!     *********************** tested. 82.01.07 ***********************
!
!    /usage/
!           call kroundrv(u, am, j, un)
      subroutine kroundrv(u, am, j, un)
      implicit none
!
!       Let u be m*10**n in normalzed form.  this program  adjusts m so
!       that it is an integral of am as follows.  The  adjusted value
!       is put in un.
!
      real u !   input.  real value to be adjusted
      real am !  input.  reference constant  used for adjusting as follows
!       am may be 0.5, 1., 1.5, 2., 5., 10. or other values
!       if am=0.5, new m will be one of 1, 1.5, 2,...9.5, 10
!             1.    //                  1,2,3,.... 10
!             2.    //                  1,2,4,6,8,10
!             5.    //                  1,5,10
      integer j ! input.  one of -1,0 or 1 
!               to signify that 1) un be <= u, 2) un be neares
!               to u, or 3) un be >= u, respectively.
      real un !  adjusted value of u.
!
!       the new value of m will be adjusted by reference to j.
!       if u is negative, absolute is taken for adjusting and sign is added
!       after adjustment.
!       if u=0, un=0 results.  if am<=0, result is not guaranteed.
!
!
      logical small
      real ua, an, em, tmp, atmp, ux
      integer n
!
      if(u .eq. 0.) then 
         un = 0.
      else
         small=(j .lt. 0  .and.  u  .gt. 0.)  .or.
     *         (j .gt. 0  .and.  u  .lt. 0.)
         ua = abs(u)
!        decompose u
         an = log10(ua)
         if(an .lt. 0.) an = an-1.
         n = an
         em = ua/10.**n
!
         tmp = em/am
         atmp = aint(tmp)
         if(j .eq. 0)  then
            if(tmp .ne. atmp) em = aint(tmp+.5)*am
         elseif(small)  then
            em = atmp*am
         else
            if(tmp .ne. atmp) em = aint(tmp+1.)*am
         endif
         if(em .eq. 0.) em = 1.
         ux = em*10.**n
         if(u .lt. 0.) ux = -ux
         un=ux
      endif
      end
