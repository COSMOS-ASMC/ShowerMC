*cmz :  3.14/16 01/02/89  17.15.14  by  nick van eijndhoven (cern)
*-- author :
      subroutine pcsdat(lunin,lunout,prflag)
c
c *** preparation of data stmts. for elastic and inelastic measured  ***
c *** cross sections of pion,kaon and proton/antiproton/neutron      ***
c *** on protons. calculate from this cross sections for strange     ***
c *** baryons on protons using additive quark-quark scattering model ***
c *** nve 22-feb-1988 cern geneva ***
c
c the program produces an output file which contains the data
c statements for all the needed cross-sections in sequence "/csdat"
c these data statements may be used for routine "ghesig"
c
c lunin  = unit number for cross-section data file
c lunout = unit number for data statements to be written
c prflag = printout flag (true/false)
c
c origin : h.fesefeldt 06-oct-87 (routine csdata)
c
c quark scattering amplitudes
c <pp |pp > = <nn |nn > = p
c <pn |pn > = <np |np > = p
c <pn |np > = <np |pn > = 0
c <pl |pl > = <nl |nl > = p - s
c <plb|plb> = <nlb|nlb> = p - s
c <pl |lp > = <nl |ln > = 0
c <pbn|pbn> = <nbp|nbp> = p
c <pbp|pbp> = <nbn|nbn> = p + a
c <pbp|nbn> = <nbn|pbp> = a
c <ll |ll > =           = (p - s)**2/p
c <lbl|lbl> =           = (p + a)*(p + s)**2/p**2
c <pbp|lbl> = <nbn|lbl> = a*(p - s)/p
c
c cross sections
c pi-  p  = 6p + 2a
c pi+  p  = 6p +  a
c k-   p  = 6p + 2a - 3s
c k0s  p  = 6p + a/2- 3s
c k0l  p  = 6p + a/2- 3s
c k0   p  = 6p      - 3s  ==> s(k0)=s(k+)
c k0b  p  = 6p + a  - 3s  ==> s(k0b)=s(k+)/4+s(k0l)/3+5s(k-)/12
c k+   p  = 6p      - 3s
c p    p  = 9p
c pb   p  = 9p + 5a
c n    p  = 9p
c
c for the following reactions we took the amplitudes
c
c p = s(p p)/9
c a = (s(pb p)-s(p p))/5
c s = 2s(p p)/9-s(k+ p)/3
c
c =====>
c
c nb   p  = 9p + 4a
c l    p  = 9p      - 3s
c lb   p  = 9p + 2a - 3s
c s-   p  = 9p      - 3s
c s+   p  = 9p      - 3s
c s-b  p  = 9p + 2a - 3s
c s+b  p  = 9p + 4a - 3s
c x0   p  = 9p      - 6s
c x-   p  = 9p      - 6s
c x0b  p  = 9p + 2a - 6s
c x-b  p  = 9p +  a - 6s
c
c --- dimension statements for cross section data ---
      dimension plab(41),csel(35,41),csin(35,41),cspiel(3,41),
     $          cspiin(3,41),cspnel(3,41),cspnin(3,41),
     $          elab(17),cnlwat(15),cnlwel(15,17),cnlwin(15,17),
     $          cscap(100),ekfiss(21),csfiss(4,21)
c
c
      logical prflag
c
      dimension a(6,41)
      character*72 ia
      dimension ax(9),lx(3),mx(3),bx(3),cx(9),ex(9)
      data ax/423.,78.,-54.,78.,34.25,-7.5,-54.,-7.5,27./
c
c --- initialize all arrays ---
      do 8000 j=1,41
      plab(j)=0.0
      do 8001 i=1,35
      csel(i,j)=0.0
      csin(i,j)=0.0
      if (i .gt. 3) go to 8001
      cspiel(i,j)=0.0
      cspiin(i,j)=0.0
      cspnel(i,j)=0.0
      cspnin(i,j)=0.0
 8001 continue
 8000 continue
c
      do 8002 j=1,17
      elab(j)=0.0
      do 8003 i=1,15
      cnlwel(i,j)=0.0
      cnlwin(i,j)=0.0
 8003 continue
      if (j .gt. 15) go to 8002
      cnlwat(j)=0.0
 8002 continue
c
      do 8004 j=1,21
      ekfiss(j)=0.0
      do 8005 i=1,4
      csfiss(i,j)=0.0
 8005 continue
 8004 continue
c
      do 8006 i=1,100
      cscap(i)=0.0
 8006 continue
c
      ip=0
      iat=-3
      cnlwat( 1)=   1.
      cnlwat( 2)=  16.
      cnlwat( 3)=  27.
      cnlwat( 4)=  56.
      cnlwat( 5)=  59.
      cnlwat( 6)=  64.
      cnlwat( 7)=  91.
      cnlwat( 8)= 112.
      cnlwat( 9)= 119.
      cnlwat(10)= 127.
      cnlwat(11)= 137.
      cnlwat(12)= 181.
      cnlwat(13)= 207.
      cnlwat(14)= 209.
      cnlwat(15)= 238.
      call ucopy(ax,cx,9)
      call minv(ax,3,det,lx,mx)
      ex(1)=ax(1)*cx(1)+ax(4)*cx(2)+ax(7)*cx(3)
      ex(2)=ax(1)*cx(4)+ax(4)*cx(5)+ax(7)*cx(6)
      ex(3)=ax(1)*cx(7)+ax(4)*cx(8)+ax(7)*cx(9)
      ex(4)=ax(2)*cx(1)+ax(5)*cx(2)+ax(8)*cx(3)
      ex(5)=ax(2)*cx(4)+ax(5)*cx(5)+ax(8)*cx(6)
      ex(6)=ax(2)*cx(7)+ax(5)*cx(8)+ax(8)*cx(9)
      ex(7)=ax(3)*cx(1)+ax(6)*cx(2)+ax(9)*cx(3)
      ex(8)=ax(3)*cx(4)+ax(6)*cx(5)+ax(9)*cx(6)
      ex(9)=ax(3)*cx(7)+ax(6)*cx(8)+ax(9)*cx(9)
      if(prflag) print 7001,ax,det,ex
 7001 format(1h ,'csdata matrix inversion ',9f10.6,2x,f10.1/
     *       1h ,'       unit matrix      ',9f10.6)
    1 read(lunin,1001) ia
      if (ia(1:3) .eq. 'end') go to 100
      if(prflag) print 1002,ia
      read(lunin,1001) ia
      if(prflag) print 1002,ia
      read(lunin,1001) ia
      if(prflag) print 1002,ia
 1001 format(a)
 1002 format(1h ,3x,a)
      ip=ip+1
      if(ip.gt.11) goto 95
      if(ip.gt.10) goto 90
      if(ip.gt. 5) goto 3
      do 2 i=1,41
      read(lunin,1003) plab(i),(a(j,i),j=1,6)
      if(prflag) print 1004,plab(i),(a(j,i),j=1,6)
 1003 format(7f10.2)
 1004 format(1h ,7f10.2)
    2 continue
      goto 5
    3 read(lunin,1001) ia
      if(prflag) print 1002,ia
      do 4 i=1,17
      read(lunin,1003) elab(i),(a(j,i),j=1,6)
      if(prflag) print 1010,elab(i),(a(j,i),j=1,6)
 1010 format(1h ,f10.4,6f10.2)
    4 continue
    5 if(ip.eq. 1) goto 20
      if(ip.eq. 2) goto 30
      if(ip.eq. 3) goto 40
      if(ip.eq. 4) goto 60
      if(ip.eq. 5) goto 70
      if(ip.ge. 6) goto 80
      goto 1
   20 do 21 i=1,41
      csel(9,i)=a(1,i)
      csin(9,i)=a(2,i)
      csel(7,i)=a(3,i)
   21 csin(7,i)=a(4,i)
      goto 1
   30 do 31 i=1,41
      csel(13,i)=a(1,i)
      csin(13,i)=a(2,i)
      csel(10,i)=a(3,i)
      csin(10,i)=a(4,i)
      csel(12,i)=a(5,i)
      csin(12,i)=a(6,i)
   31 continue
      goto 1
   40 if(prflag) print 1008
 1008 format(1h0,'quark parton scattering amplitudes'/
     *1h ,'    p(gev/c)     p         a         s          p          a
     *         s'/
     *1h ,'                 least square fit               calculated')
      do 41 i=1,41
      csel(14,i)=a(1,i)
      csin(14,i)=a(2,i)
      csel(15,i)=a(3,i)
      csin(15,i)=a(4,i)
      csel(16,i)=a(5,i)
      csin(16,i)=a(6,i)
      bx(1)=  6.*(csel( 7,i)+csin( 7,i))+ 6.*(csel( 9,i)+csin( 9,i))
     *      + 6.*(csel(10,i)+csin(10,i))+ 6.*(csel(12,i)+csin(12,i))
     *      + 6.*(csel(13,i)+csin(13,i))+ 9.*(csel(14,i)+csin(14,i))
     *      + 9.*(csel(15,i)+csin(15,i))+ 9.*(csel(16,i)+csin(16,i))
      bx(2)=     (csel( 7,i)+csin( 7,i))+    (csel( 9,i)+csin( 9,i))
     *      + .5*(csel(12,i)+csin(12,i))+ 5.*(csel(15,i)+csin(15,i))
      bx(3)=- 3.*(csel(10,i)+csin(10,i))- 3.*(csel(12,i)+csin(12,i))
     *      - 3.*(csel(13,i)+csin(13,i))
      pam =ax(1)*bx(1)+ax(4)*bx(2)+ax(7)*bx(3)
      aam =ax(2)*bx(1)+ax(5)*bx(2)+ax(8)*bx(3)
      sam =ax(3)*bx(1)+ax(6)*bx(2)+ax(9)*bx(3)
c     ratkp=csel(10,i)/(csel(10,i)+csin(10,i))
c     ratkm=csel(13,i)/(csel(13,i)+csin(13,i))
      ratp =csel(14,i)/(csel(14,i)+csin(14,i))
      ratpb=csel(15,i)/(csel(15,i)+csin(15,i))
c     csel(11,i)=(6.*pam           -3.*sam   )  * ratkp
c     csin(11,i)=(6.*pam           -3.*sam   )  *(1.-ratkp)
c     csel(12,i)=(6.*pam  +   aam  -3.*sam   )  * ratkm
c     csin(12,i)=(6.*pam  +   aam  -3.*sam   )  *(1.-ratkm)
      csel(11,i)= csel(10,i)
      csin(11,i)= csin(10,i)
      csel(12,i)= csel(10,i)/4.+csel(12,i)/3.+5.*csel(13,i)/12.
      csin(12,i)= csin(10,i)/4.+csin(12,i)/3.+5.*csin(13,i)/12.
      pamp= (csel(14,i)+csin(14,i))/9.
      aamp= (csel(15,i)+csin(15,i)-csel(14,i)-csin(14,i))/5.
      samp= 2.*(csel(14,i)+csin(14,i))/9.-(csel(10,i)+csin(10,i))/3.
      if(plab(i).lt.0.59) samp = 0.
      if(prflag) print 1009,plab(i),pam,aam,sam,pamp,aamp,samp
 1009 format(1h ,7f10.2)
      csel(17,i)=(9.*pamp +4.*aamp           )  *ratpb
      csin(17,i)=(9.*pamp +4.*aamp           )  *(1.-ratpb)
      csel(18,i)=(9.*pamp          -3.*samp  )  *ratp
      csin(18,i)=(9.*pamp          -3.*samp  )  *(1.-ratp)
      csel(19,i)=(9.*pamp +2.*aamp -3.*samp  )  *ratpb
      csin(19,i)=(9.*pamp +2.*aamp -3.*samp  )  *(1.-ratpb)
      csel(20,i)=(9.*pamp          -3.*samp  )  *ratp
      csin(20,i)=(9.*pamp          -3.*samp  )  *(1.-ratp)
      csel(22,i)=(9.*pamp          -3.*samp  )  *ratp
      csin(22,i)=(9.*pamp          -3.*samp  )  *(1.-ratp)
      csel(23,i)=(9.*pamp +4.*aamp -3.*samp  )  *ratpb
      csin(23,i)=(9.*pamp +4.*aamp -3.*samp  )  *(1.-ratpb)
      csel(25,i)=(9.*pamp +2.*aamp -3.*samp  )  *ratpb
      csin(25,i)=(9.*pamp +2.*aamp -3.*samp  )  *(1.-ratpb)
      csel(26,i)=(9.*pamp          -6.*samp  )  *ratp
      csin(26,i)=(9.*pamp          -6.*samp  )  *(1.-ratp)
      csel(27,i)=(9.*pamp          -6.*samp  )  *ratp
      csin(27,i)=(9.*pamp          -6.*samp  )  *(1.-ratp)
      csel(28,i)=(9.*pamp +2.*aamp -6.*samp  )  *ratpb
      csin(28,i)=(9.*pamp +2.*aamp -6.*samp  )  *(1.-ratpb)
      csel(29,i)=(9.*pamp +   aamp -6.*samp  )  *ratpb
      csin(29,i)=(9.*pamp +   aamp -6.*samp  )  *(1.-ratpb)
c --- set cross sections for the omega and anti-omega to the values ---
c --- as for xi- and anti xi- respectively ---
      csel(33,i)=csel(27,i)
      csin(33,i)=csin(27,i)
      csel(34,i)=csel(29,i)
      csin(34,i)=csin(29,i)
 41   continue
      do 42 i=1,41
      do 42 j=1,35
      if(csel(j,i).lt.0.) csel(j,i)=0.
      if(csin(j,i).lt.0.) csin(j,i)=0.
   42 continue
      goto 1
   60 do 61 i=1,41
      cspiel(1,i)=a(1,i)-a(2,i)
      cspiin(1,i)=a(2,i)
      cspiel(2,i)=a(3,i)-a(4,i)
      cspiin(2,i)=a(4,i)
      cspiel(3,i)=a(5,i)-a(6,i)
   61 cspiin(3,i)=a(6,i)
      goto 1
   70 do 71 i=1,41
      cspnel(1,i)=a(1,i)-a(2,i)
      cspnin(1,i)=a(2,i)
      cspnel(2,i)=a(3,i)-a(4,i)
      cspnin(2,i)=a(4,i)
      cspnel(3,i)=a(5,i)-a(6,i)
   71 cspnin(3,i)=a(6,i)
      goto 1
   80 iat=iat+3
      do 81 i=1,17
      cnlwel(iat+1,i)=a(1,i)-a(2,i)
      cnlwin(iat+1,i)=a(2,i)
      cnlwel(iat+2,i)=a(3,i)-a(4,i)
      cnlwin(iat+2,i)=a(4,i)
      cnlwel(iat+3,i)=a(5,i)-a(6,i)
   81 cnlwin(iat+3,i)=a(6,i)
      goto 1
   90 i2=0
      i1=0
      do 91 i=1,17
      i1=i2+1
      i2=i2+6
      if(i2.gt.100) i2=100
      read(lunin,1011) ia,(cscap(j),j=i1,i2)
      if(prflag)
     *print 1011,ia,(cscap(j),j=i1,i2)
 1011 format(a10,6f10.2)
   91 continue
      goto 1
   95 read(lunin,1001) ia
      if(prflag) print 1002,ia
      do 96 i=1,21
      read(lunin,1012) ekfiss(i),(csfiss(j,i),j=1,4)
      if(prflag) print 1012,ekfiss(i),(csfiss(j,i),j=1,4)
 1012 format(1x,f9.4,4f10.0)
   96 continue
c**
  100 if(.not.prflag) goto 200
      print 1005
 1005 format(1h0,'listings of all particle cross sections on protons mea
     *sured or calculated')
      do 105 k=1,7
      k1=(k-1)*5+1
      k2=k*5
      print 1006,k1,k2
 1006 format(1h0,'ipart1 - ipart2',2x,i3,'-',i3)
      do 104 i=1,41
      print 1007,plab(i),(csel(j,i),csin(j,i),j=k1,k2)
 1007 format(1h ,f10.2,2x,5(2x,2f10.2))
  104 continue
  105 continue
c
  200 continue
c
c --- check for un-physical values ---
      do 9000 j=1,41
      if (plab(j) .lt. 0.0) plab(j)=0.0
      do 9001 i=1,35
      if (csel(i,j) .lt. 0.0) csel(i,j)=0.0
      if (csin(i,j) .lt. 0.0) csin(i,j)=0.0
      if (i .gt. 3) go to 9001
      if (cspiel(i,j) .lt. 0.0) cspiel(i,j)=0.0
      if (cspiin(i,j) .lt. 0.0) cspiin(i,j)=0.0
      if (cspnel(i,j) .lt. 0.0) cspnel(i,j)=0.0
      if (cspnin(i,j) .lt. 0.0) cspnin(i,j)=0.0
 9001 continue
 9000 continue
c
      do 9002 j=1,17
      if (elab(j) .lt. 0.0) elab(j)=0.0
      do 9003 i=1,15
      if (cnlwel(i,j) .lt. 0.0) cnlwel(i,j)=0.0
      if (cnlwin(i,j) .lt. 0.0) cnlwin(i,j)=0.0
 9003 continue
      if (j .gt. 15) go to 9002
      if (cnlwat(j) .lt. 0.0) cnlwat(j)=0.0
 9002 continue
c
      do 9004 j=1,21
      if (ekfiss(j) .lt. 0.0) ekfiss(j)=0.0
      do 9005 i=1,4
      if (csfiss(i,j) .lt. 0.0) csfiss(i,j)=0.0
 9005 continue
 9004 continue
c
      do 9006 i=1,100
      if (cscap(i) .lt. 0.0) cscap(i)=0.0
 9006 continue
c
c --- write out the cross sections in the form of data stmts. ---
      write(lunout,2001)
 2001 format('+keep,/csdat.'/
     $ 'c --- cross-section data by "pcsdat" 01-feb-1989 ---'/'c')
c
      write(lunout,2010) (plab(i),i=1,41)
 2010 format(6x,'data plab /'/
     $ 8(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,'/')
c
      do 220 i=1,35
      write(lunout,2020) i,(csel(i,j),j=1,41)
 2020 format(6x,'data (csel(',i2,',j),j=1,41) /'/
     $ 8(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,'/')
 220  continue
c
      do 230 i=1,35
      write(lunout,2030) i,(csin(i,j),j=1,41)
 2030 format(6x,'data (csin(',i2,',j),j=1,41) /'/
     $ 8(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,'/')
 230  continue
c
      do 240 i=1,3
      write(lunout,2040) i,(cspiel(i,j),j=1,41)
 2040 format(6x,'data (cspiel(',i2,',j),j=1,41) /'/
     $ 8(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,'/')
 240  continue
c
      do 250 i=1,3
      write(lunout,2050) i,(cspiin(i,j),j=1,41)
 2050 format(6x,'data (cspiin(',i2,',j),j=1,41) /'/
     $ 8(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,'/')
 250  continue
c
      do 260 i=1,3
      write(lunout,2060) i,(cspnel(i,j),j=1,41)
 2060 format(6x,'data (cspnel(',i2,',j),j=1,41) /'/
     $ 8(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,'/')
 260  continue
c
      do 270 i=1,3
      write(lunout,2070) i,(cspnin(i,j),j=1,41)
 2070 format(6x,'data (cspnin(',i2,',j),j=1,41) /'/
     $ 8(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,'/')
 270  continue
c
      write(lunout,2080) (elab(i),i=1,17)
 2080 format(6x,'data elab /'/
     $ 3(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,',',g12.5,'/')
c
      write(lunout,2090) (cnlwat(i),i=1,15)
 2090 format(6x,'data cnlwat /'/
     $ 2(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',4(g12.5,','),g12.5,'/')
c
      do 280 i=1,15
      write(lunout,2100) i,(cnlwel(i,j),j=1,17)
 2100 format(6x,'data (cnlwel(',i2,',j),j=1,17) /'/
     $ 3(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,',',g12.5,'/')
 280  continue
c
      do 290 i=1,15
      write(lunout,2110) i,(cnlwin(i,j),j=1,17)
 2110 format(6x,'data (cnlwin(',i2,',j),j=1,17) /'/
     $ 3(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,',',g12.5,'/')
 290  continue
c
c --- write cscap in two parts because of continuation card limit ---
      write(lunout,2120) (cscap(i),i=1,100)
 2120 format(6x,'data (cscap(j),j=1,50) /'/
     $ 9(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',4(g12.5,','),g12.5,'/'/
     $       6x,'data (cscap(j),j=51,100) /'/
     $ 9(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',4(g12.5,','),g12.5,'/')
c
      write(lunout,2130) (ekfiss(i),i=1,21)
 2130 format(6x,'data ekfiss /'/
     $ 4(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,'/')
c
      do 300 i=1,4
      write(lunout,2140) i,(csfiss(i,j),j=1,21)
 2140 format(6x,'data (csfiss(',i2,',j),j=1,21) /'/
     $ 4(5x,'$ ',5(g12.5,',')/),
     $   5x,'$ ',g12.5,'/')
 300  continue
c
      write(lunout,2999)
 2999 format('c --- end of cross section data statements ---'/'c')
c
      return
      end
