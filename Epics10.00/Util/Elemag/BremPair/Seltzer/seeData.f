c       Seltzer & Berger's data brems01 ~ bremes10 are
c     expanded to data of each element ranging from Z=1 to 100
c     The upper energy may be specified. (GIVE Emax)
c
c      You can inspect each data file by gnuplot
c
c     files SZ1, SZ2, ... SZ100 will be created, each has data of
c     Z=1, Z=2,...Z=100.
c The format of each line is
c      
c        v   cross-section  Ek
c
c  where v = Eg/(Ee-me)
c  cross-section = (beta/Z)**2 v ds/dv (mb)
c  Ek = Electron kinetic energy in GeV
c
      integer ekmax, vmax
      parameter (ekmax = 57, vmax = 30)
      integer i, j, k, z1, z, ff
      character*4 name
      character*8 brem(10)
      real*4 xsec(vmax)
      real*4 ke(ekmax),  kv(vmax)
      real*4  Emax
      data Emax/10001./  !  give in MeV
c         kinetic energy of incident electron in MeV
       data ke/
     * 0.0010, 0.00150, 0.0020, 0.0030, 0.0040, 0.0050,
     * 0.0060, 0.008, 0.010, 0.01500, 0.020, 0.030,
     * 0.0400, 0.050, 0.060, 0.080, 0.100, 0.150,
     * 0.2000, 0.300, 0.400, 0.500, 0.600, 0.800,
     * 1.0000, 1.500, 2.000, 3.000, 4.000, 5.000,
     * 6.0000, 8.000, 10.00,15.000, 20.00, 30.00,
     * 40.000, 50.000, 60.000, 80.000,  100.00,   150.0,
     * 200.00, 300.00, 400.0,  500.0, 600.0, 800.0,
     * 1000.0, 1500.0,  2000.0, 3000.0,  4000.0,  5000.0,
     * 6000.0,  8000.0, 10000.0/
c         fractional energy of gamma Eg/Ek    
       data kv/
     * 0.0, 0.050, 0.10,     0.150, 0.20, 0.250,
     * 0.30, 0.350, 0.40, 0.450, 0.50, 0.550,
     * 0.60, 0.650, 0.70, 0.750, 0.80, 0.850,
     * 0.90, 0.9250, 0.950, 0.970, 0.990, 0.9950,
     * 0.9990, 0.99950, 0.99990, 0.99995, 0.99999, 1.00/

       data brem/
     * 'brems01', 'brems02', 'brems03', 'brems04', 'brems05',
     * 'brems06', 'brems07', 'brems08', 'brems09', 'brems10'
     * /      

c       brem01: containes data for Z=1 to 10
c       brem02:                    Z=11 to 20
c ..
c       brem10:                    Z=91 to 100
c
       do ff = 1, 10
          open(11,file=brem(ff))
          write(*, *) ' file ', brem(ff), ' opened'
          z1 = (ff-1)*10 + 1
          do z = z1, z1+9
             if(z .le. 9) then
                write(name, '("SZ",i1)') z
             elseif(z .le. 99) then
                write(name, '("SZ",i2)') z
             else
                write(name, '("SZ",i3)') z
             endif
             open(10, file=name)
             write(*,*) ' file=',name, ' opened'
             do i = 1,  ekmax
                read(11, '(6f12.5)') (xsec(k), k = 1,  vmax)
                if(ke(i) .le. Emax) then
                   do j =  1, vmax
                      write(10, 
     *                 '(2f12.5,e12.3)') kv(j), xsec(j), ke(i)/1000. ! to GeV
                   enddo
                   write(10,*)
                endif
             enddo
             close(10)
          enddo
          close(11)
       enddo
      end

      
