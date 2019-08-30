!     ******************************************************************
!     *     psuedo-random number generators.   rndc, rndd, rnde        *
!     *                                                                *
!     *    In normal applications, the user may use only rndc as       *
!     *                                                                *
!     *       call rndc(u)                                             *
!     *                                                                *
!     *  where u is the output random number in (0., 1.) (excluding    *
!     *  the boundary values).  u is a real*8 variable.                *
!     *----------------------------------------------------------------*
!     *  Internally, there are 3 independent generators.               *
!     *  rndc and rndd use 2 generators: rnd1u or  rnd2u.  the default *
!     *  routine is rnd1u.  rnd1u/rnd2u can be specifed by calling the *
!     *  rndsw subroutine.  rnde uses always the 3rd generators.       *
!     *  the differece between rndc and rndd is the calling sequence.  *
!     *                                                                *
!     * if the user wants to have random numbers in an array,          *
!     *                                                                *
!     *       call rndd(ua, n)                                         *
!     *                                                                *
!     *  may be used, where ua  is an array  to get n random variables *
!     *  in (0., 1.).                                                  *
!     *                                                                *
!     *  for rnde, the calling sequence is                             *
!     *                                                                *
!     *      call rnde(ua, n)                                          *
!     *                                                                *
!     *  which is the same as that for rndd.                           *
!     *                                                                *
!     *    To switch the generators, the user may use                  *
!     *                                                                *
!     *      call rndsw(jold, jnew)                                    *
!     *                                                                *
!     *  where jold is the output value (1 or 2) indicating the        *
!     *  generator (rnd1u or rnd2u) currently used. jnew is the input  *
!     *  integer  (1 or 2) to specify the generator from the next call *
!     *  on.  This can be called any time.                             *
!     *    In some applications, the user may wish to initialize the   *
!     *  random number generators using explicit seeds, or the  user   *
!     *  may need to save the current status of the generators for     *
!     *  restarting the job in a later time, as if the job were        *
!     *  continued without halt.                                       *
!     *                                                                *
!     *   To give explicit seeds for generator 1,  use                 *
!     *                                                                *
!     *      call rnd1i(ir1)   for generator 1                         *
!     *      call rnd2i(ir2)   for generator 2                         *
!     *      call rnd3i(ir3)   for generator 3                         *
!     *                                                                *
!     *  where the argument is the input integer(s).                   *
!     *     ir1: ir1(1) and ir2(2) should contain integer seeds.       *
!     *     ir2: an integer seed.                                      *
!     *     ir3: an intgeger seed.                                     *
!     *  These can be called any time.                                 *
!     *                                                                *
!     *    To save the current status of the generator, use            *
!     *                                                                *
!     *      call rnd1s(ir1sv)     for generator 1                     *
!     *      call rnd2s(ir2sv)     for generator 2                     *
!     *      call rnd3s( r3sv)     for generator 3                     *
!     *                                                                *
!     *  where  ir1sv is an output integer array to get 2 integers,    *
!     *         ir2sv is an output intgger array to get 25 integers,   *
!     *         r3sv  is an output real*8 array to get 102 reals.      *
!     *  These can be called any time.                                 *
!     *                                                                *
!     *                                                                *
!     *     To restore the status when the rndis was called, (i=1,2,3),*
!     *  use,                                                          *
!     *                                                                *
!     *        call rnd1r(ir1sv)  for generator 1                      *
!     *        call rnd2r(ir2sv)  for generator 2                      *
!     *        call rnd3r(r3sv)   for generator 3                      *
!     *  where the arguments are the input which are obtained by       *
!     *  calling rndis (i=1,2,3).                                      *
!     *  These can be called any time.                                 *
!     *                                                                *
!     *                                                                *
!     *  ------------------------------------------------------------  *
!     *    General features of these generators.                       *
!     *                                                                *
!     *      They are all reliable generators with a long period and   *
!     *  free from undesirable correlations between successive random  *
!     *  numbers.  They are given in all standard fortran code so that *
!     *  they can be used in any macihne with integer/real expressed in*
!     *  32 or more bits.                                              *
!     *       Reference.                                               *
!     *                                                                *
!     ********************* tested 91.08.09 ******************k.k*******
      subroutine rndc(u)
         implicit none
         real*8 u
         integer n
         real*8 ua(n)
         real(8)::u1(1)
!
         integer iseed(25), ir1(2), ir2, jold, jnew
         logical first2/.true./
         integer jsw/1/
         save jsw, first2
!
         if(jsw .eq. 1) then
!!             call rnd1u(u, 1)  ! jaxa complains at execution time
                              ! since u is not array
             call rnd1u(u1, 1)
             u = u1(1)
         elseif(jsw .eq. 2) then
             if(first2) then
                first2=.false.
                call rnd2ix(31415926)
             endif
!!             call rnd2u(u, 1)
             call rnd2u(u1, 1)
             u = u1(1)
         else
             write(*,*) ' switch value error=',jsw, ' in rndc '
             stop 9999
         endif
         return
!     ***********
      entry rnd1r(ir1)
!     **********
         call rnd1i(ir1)
         return
!     ***********
      entry rndd(ua, n)
!     **********
         if(jsw .eq. 1) then
             call rnd1u(ua, n)
         elseif(jsw .eq. 2) then
             if(first2) then
                first2=.false.
                call rnd2ix(31415926)
             endif
             call rnd2u(ua, n)
         else
             write(*,*) ' switch value error=',jsw, ' in rndd '
             stop 9999
         endif
         return
!     ***********
      entry rnd2i(ir2)
!     ***********
         call rnd2ix(ir2)
         first2=.false.
         return
!     ***********
      entry rnd2r(iseed)
!     ***********
         call rnd2rx(iseed)
         first2=.false.
         return
!     ***********
      entry rndsw(jold, jnew)
!     ***********
         jold=jsw
         jsw=jnew
         return
      end
      subroutine rnd1u(ua, n)
!            random number generator given by L'ecuyer in
!            comm. acm vol 31, p.742, 1988
!            modified by f. James to return a vector of numbers
!            modified by K.K
!    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     calling sequences for rnd1u:                                   ++
!         call rnd1u (ua, n)        returns a vector ua of n         ++
!                      64-bit random floating point numbers between  ++
!                      zero and one.                                 ++
!         call rnd1i(irin)      initializes the generator from two   ++
!                      32-bit integer array irin                     ++
!         call rnd1s(irout)     outputs the current values of the    ++
!                      two integer seeds, to be used for restarting  ++
!    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          implicit none
          integer n
          real*8 ua(n)
! 
         integer irin(2), irout(2), iseed1, iseed2, i, k, iz
          save iseed1,iseed2
          data iseed1,iseed2 /12345,67890/
!
           do   i= 1, n
             k = iseed1/53668
             iseed1 = 40014*(iseed1 - k*53668) - k*12211
             if(iseed1 .lt. 0) iseed1=iseed1+2147483563
!
             k = iseed2/52774
             iseed2 = 40692*(iseed2 - k*52774) - k* 3791
             if(iseed2 .lt. 0) iseed2=iseed2+2147483399
!
             iz = iseed1 - iseed2
             if(iz .lt. 1) iz = iz + 2147483562
!
             ua(i) = iz * 4.656613d-10
           enddo
          return
!     ****************
      entry rnd1i(irin)
!     ****************
          iseed1 = irin(1)
          iseed2 = irin(2)
          return
!
!     ****************
      entry rnd1s(irout)
!     ****************
          irout(1)= iseed1
          irout(2)= iseed2
      end
      subroutine rnd2u(ua,n)
!         add-and-carry random number generator proposed by
!         Marsaglia and Zaman in siam j. scientific and statistical
!         computing, to appear probably 1990.
!         modified with enhanced initialization by F. James, 1990
!         modified by K.K
!    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     calling sequences for rnd2u:                                   ++
!         call rnd2u (ua, n)        returns a vector ua of n         ++
!                      64-bit random floating point numbers between  ++
!                      zero and one.                                 ++
!         call rnd2ix(int)     initializes the generator from one    ++
!                      32-bit integer int                            ++
!         call rnd2rx(ivec)    restarts the generator from vector    ++
!                      ivec of 25 32-bit integers (see rnd2s)        ++
!         call rnd2s(ivec)     outputs the current values of the 25  ++
!                    32-bit integer seeds, to be used for restarting ++
!    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
      integer n
      real*8 ua(n)
!
      integer ns, ns1, i, isd, k
      parameter (ns=24, ns1=ns+1)
      real*8 seeds(ns), uni
      integer iseeds(ns), isdin(ns1), isdout(ns1), icarry
      real*8 twop24, twom24
      integer itwo24, icons, i24, j24, ivec, jseed, inseed
      parameter (twop24=2.**24)
      parameter (twom24=2.**(-24), itwo24=2**24, icons=2147483563)
      real*8 carry
      save  i24, j24, carry, seeds
      data i24, j24, carry/ns, 10, 0./
!
!          the generator proper: "subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida state university, march, 1989
!
       do   ivec= 1, n
          uni = seeds(i24) - seeds(j24) - carry
          if(uni .lt. 0.d0) then
             uni = uni + 1.0d0
             carry = twom24
          else
             carry = 0.
          endif
!                 avoid exact zero
          if(uni .eq. 0.d0) then
              uni = twom24
          endif
          seeds(i24) = uni
          i24 = i24 - 1
          if(i24 .eq. 0) i24 = ns
          j24 = j24 - 1
          if(j24 .eq. 0) j24 = ns
          ua(ivec) = uni     
       enddo
       return
!           entry to restore the previous run
!     ******************
      entry rnd2rx(isdin)
!     ******************
       do   i= 1, ns
         seeds(i) =isdin(i)*twom24
       enddo
      carry = mod(isdin(ns1),10)*twom24
      isd = isdin(ns1)/10
      i24 = mod(isd,100)
      isd = isd/100
      j24 = isd
      return
!                    entry to get current status
!     ******************
      entry rnd2s(isdout)
!     ******************
       do   i= 1, ns
         isdout(i) = int(seeds(i)*twop24)
       enddo
      if(carry .gt. 0.) then
         icarry=1
      else
         icarry=0
      endif
      isdout(ns1) = 1000*j24 + 10*i24 + icarry
      return
!                    entry to initialize from one integer
!     ******************
      entry rnd2ix(inseed)
!     ******************
      jseed = inseed
      do   i= 1, ns
         k = jseed/53668
         jseed = 40014*(jseed-k*53668) -k*12211
         if(jseed .lt. 0) jseed = jseed+icons
         iseeds(i) = mod(jseed,itwo24)
      enddo
      do   i= 1,ns
         seeds(i) =iseeds(i)*twom24
      enddo
      i24 = ns
      j24 = 10
      carry = 0.
      if(seeds(ns) .lt. seeds(14)) carry = twom24
      end
      subroutine rnde(ua,n)
!           random number generator proposed by marsaglia and zaman
!           in report fsu-scri-87-50
!           modified by f. james, 1988 and 1989, to generate a vector
!           of pseudorandom numbers ua of length n.
!           modified by k.k
!    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     calling sequences for rnde:                                    ++
!         call rnde(ua, n)         returns a vector ua of n          ++
!                      32-bit random floating point numbers between  ++
!                      zero and one.                                 ++
!         call rnd3i(i1)          initializes the generator from one ++
!                      32-bit integer i1
!    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
!     integer  n
      integer,intent(in):: n
      real*8 ua(n), sina(102)
      logical first/.true./
      save  first
      integer ijkl, ijklin
!
      if(first) then
         first=.false.
!           default initialization. user has called rnde without rnd3i.
         ijkl = 54217137
         call rnd3ix(ijkl)
      endif
      call rnd3x(ua, n)
      return
!         initializing routine for rnde, may be called before
!         generating pseudorandom numbers with rnde. the input
!         values should be in the ranges:  0<=ijklin<=900 ooo ooo
!     **************
      entry rnd3i(ijklin)
!     *************
      first=.false.
      call rnd3ix(ijklin)
      return
!     ************
      entry rnd3r(sina)
!     ************
      first=.false.
      call rnd3rx(sina)
      end

      
      subroutine  rnd3ix(ijkl)
!             the standard values in marsaglia's paper, ijkl=54217137
          implicit none
!     integer n
          integer,intent(in):: n
!      If next is used Gfortran on Mac gets mad.
!     real*8 ua(n), u(97), uni, s, t, zuni
!     So we devide it next two.  ua(*) is essential
!         
          real(8):: ua(*)
          real(8):: u(97), uni, s, t, zuni
          real*8 sina(102), sout(102)
          integer  jj, m
          integer ns, ijkl, i97, j97, ij, kl, i, j, k, l, ii, ivec
          real*8 twom24, c, cd, cm
          parameter (ns=24, twom24=2.**(-24))
          save c, cd, cm, i97, j97
!
          ij = ijkl/30082
          kl = ijkl - 30082*ij
          i = mod(ij/177, 177) + 2
          j = mod(ij, 177) + 2
          k = mod(kl/169, 178) + 1
          l = mod(kl, 169)
           do   ii= 1, 97
              s = 0.
              t = .5
               do   jj= 1, ns
                  m = mod(mod(i*j,179)*k, 179)
                  i = j
                  j = k
                  k = m
                  l = mod(53*l+1, 169)
                  if(mod(l*m,64) .ge. 32) s = s+t
                  t = 0.5*t
               enddo
              u(ii) = s
           enddo
          c = 362436.*twom24
          cd = 7654321.*twom24
          cm = 16777213.*twom24
          i97 = 97
          j97 = 33
          return
!     ****************
      entry rnd3x(ua, n)
!     ****************
       do   ivec= 1, n
          uni = u(i97)-u(j97)
          if(uni .lt. 0.) uni=uni+1.
          u(i97) = uni
          i97 = i97-1
          if(i97 .eq. 0) i97=97
          j97 = j97-1
          if(j97 .eq. 0) j97=97
          c = c - cd
          if(c .lt. 0.) c=c+cm
          uni = uni-c
          if(uni .lt. 0.) uni=uni+1.
          ua(ivec) = uni
!                 replace exact zeros by uniform distr. *2**-24
          if(uni .eq. 0.) then
              zuni = twom24*u(2)
!               an exact zero here is very unlikely, but let's be safe.
              if(zuni .eq. 0.) zuni= twom24*twom24
              ua(ivec) = zuni
          endif
       enddo
      return
!     ****************** to get current status
      entry rnd3s(sout)
!     ***********
           do   i=1, 97
              sout(i)=u(i)
           enddo
          sout(98)=c
          sout(99)=cd
          sout(100)=cm
          sout(101)=i97
          sout(102)=j97
          return
!     ****************  to restore the old status
      entry rnd3rx(sina)
           do   i=1, 97
               u(i)=sina(i)
           enddo
          c=sina(98)
          cd=sina(99)
          cm=sina(100)
          i97=sina(101)
          j97=sina(102)
      end
