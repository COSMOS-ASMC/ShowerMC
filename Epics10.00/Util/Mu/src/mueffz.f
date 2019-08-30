       include  'musub.f'
!c         gives effecive z distribution  5000
         data dz/ 50000./
         dimension ztbl(100)
         open(13, file='c2s5001.#gd.data')
         call muroc0('c2g5100.rock.data(frejus)', a, z, zba, z2ba,rho,
     *   bh, bv, beta)
         call mueffz(dz, ztbl, 100, 1., zmin, zmax)
          do   i=1, 100
            write(13)( zmin+(i-1)*dz+dz/2)/100000, ztbl(i)
          enddo
       end
!
       subroutine mueffz(dz, ztbl, n, enf, zmin, zmax)
!
!    dz:  input.   bin of the ztbl (g/cm**2)
!  ztbl:  ouput.   ztbl(i), i=1, n.    ztbl(i) is the relative
!         weight of z=( (i-1)*dz+zmin, i*dz+zmin)
!     n:  input.   ztbl size
!  enf:   input.  sec(teta) enfancement factor (sec)**enf
!  zmin:  output.  min of the depth in ztbl
!  zmax:  ouptut.  max of the depth in ztbl
!
          parameter (pi=3.141592, Torad=pi/180.)
          dimension ztbl(n)
          call murocq(zmin, zmax)
          write(*,*) ' zmin=',zmin, ' zmax=',zmax
          call ksetrv(ztbl, 1, n, 0.)
           do   fai=0., 359.8, .5
              rfai=(fai+  .5/2)*Torad
               do   cs=1., .49999, -.001
                   csm=cs - .001/2
                   tet=acos(csm)
                   call murock(tet, rfai, z)
                   i=min((z-zmin)/dz +1, n)
                   ztbl(i)=ztbl(i)+ (1./csm)**enf
               enddo
           enddo
        end
