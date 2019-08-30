*cmz :  3.14/16 28/09/90  10.38.04  by  nick van eijndhoven (cern)
*-- author :    nick van eijndhoven (cern)   02/02/89
      subroutine casaom(k,int,nfl)
c
c *** cascade of omega-  bar ***
c *** nve 17-jan-1989 cern geneva ***
c
c omega- bar undergoes interaction with nucleon within nucleus.
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
c *** not yet finished ==> take omega- cascade instead ***
c
c --- initialization indicated by kginit(24) ---
      kginit(24)=1
c
c      if (nprt(4)) print 1000
 1000 format(' *casaom* not written yet ==> casom called instead')
c
      call casom(k,int,nfl)
c
 9999 continue
      return
      end
