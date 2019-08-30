!           need Zcoord.h  Zobs.h  Zpos.h Zmagfield.h
         integer NoOfSites           ! No of particle observation sites
        integer NoOfASSites
!          
         real*8 CosLatitude          ! oos of Latitude of deepest obs. site
        real*8 SinLatitude          ! sin
        real*8 CosLongitude         ! cos of Longitude
        real*8 SinLongitude         ! sin of ..

        type(coord):: DetZaxis     ! detector's Z axis in 'xyz' system
        type(coord):: DetXaxis     !  //        X    // 
        type(coord):: DetYaxis     !  //        Y    // 

        type(coord):: Xprimary     ! primary system x axis in 'xyz'
        type(coord):: Yprimary     ! primary system y axis in 'xyz'
        type(coord):: Zprimary     ! primary system z axis in 'xyz'
                                    ! these are computed in cprimxyz in
                                    ! ciniTracking in ceventLoop
        real(8)::Txyz2prim(3,3)    ! matrix to transform vector in
                              ! E-xyz into primary system
                  ! vector must be given from the oriign of
                  ! the detecor
        real(8)::Tprim2xyz(3,3)  ! inverse of Txyz2prim
        real(8)::Txyz2det(3,3) ! xyz to detector system transform mat
        real(8)::Tdet2xyz(3,3) ! inverse of above
        type(coord):: PolarInjPos  ! polar angle of the injection point in xyz.

         type(magfield):: MagfieldNED     ! mag in 'ned' at deepest obs. site
        type(magfield):: MagfieldHVA     ! mag in 'hva' at //. both in T.
        type(magfield):: MagfieldXYZ     ! mag in 'xyz' at //. both in T.

         type site
           sequence
               type(position):: pos
               real*8  zpl           ! z value in 1ry system
               real*8  mu
               real*8  minitime
         end type site 
         type assite
           sequence
               type(position):: pos
               real*8  zpl 
               real*8  mu             ! Moliere Unit
               real*8  esize          ! electron size
               real*8  age            ! size weighted age
         end type assite

          type(site):: ObsSites(0:maxNoOfSites+1)
         type(assite):: ASObsSites(maxNoOfASSites)
!            to store Ne, age of a component shower for an electron
         real*8 CompASNe(maxNoOfASSites), CompASAge(maxNoOfASSites)

         common /Zobsvc/  ASObsSites,  ObsSites,
     *     MagfieldNED, MagfieldHVA, MagfieldXYZ, 
     *     CompASNe, CompASAge,
     *     DetZaxis, DetXaxis, DetYaxis,
     *     Xprimary, Yprimary, Zprimary,
     *     Txyz2det, Tdet2xyz, Txyz2prim, Tprim2xyz,
     *     PolarInjPos, 
     *     CosLatitude, SinLatitude, CosLongitude, SinLongitude,
     *     NoOfSites,    NoOfASSites
