*cmz :  3.14/16 28/09/90  10.36.58  by  nick van eijndhoven (cern)
*-- author :    nick van eijndhoven (cern)   02/02/89
      subroutine casax0(k,int,nfl)
c
c *** cascade of xi0 bar ***
c *** nve 17-jan-1989 cern geneva ***
c
c xi0 bar undergoes interaction with nucleon within nucleus.
c check if energetically possible to produce pions/kaons.
c if not, assume nuclear excitation occurs, degrade input particle
c in energy and no other particles are produced.
c if reaction is possible find correct number of pions/protons/
c neutrons produced using an interpolation to multiplicity data.
c replace some pions or protons/neutrons by kaons or strange baryons
c according to average multiplicity per inelastic reactions.
c
      common/prntfl/inbcd,newbcd,inbin,newbin,npevt,nevtp,lprt,nprt(10)
                    logical lprt,nprt
c
c --- initialization flags for various gheisha routines ---
      common /kginit/ kginit(50)
c
c
c *** not yet finished ==> take xi0 cascade instead ***
c
c --- initialization indicated by kginit(23) ---
      kginit(23)=1
c
c      if (nprt(4)) print 1000
 1000 format(' *casax0* not written yet ==> casx0 called instead')
c
      call casx0(k,int,nfl)
c
 9999 continue
      return
      end