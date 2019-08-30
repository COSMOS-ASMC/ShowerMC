c...common block for buu 

c...Number of spatial grid.
      parameter(maxx=20,maxz=24)

c...Force acting on hadrons.
      common/jamevf/force(3,mxv)

c...Baryon density.
      common/jambuu1/rhob(-maxx:maxx,-maxx:maxx,-maxz:maxz)

c...Charge density.
      common/jambuu2/rhoz(-maxx:maxx,-maxx:maxx,-maxz:maxz)

c....Coulomb potential.
      common/jambuu3/cop(-maxx:maxx,-maxx:maxx,-maxz:maxz)

c...4-momentum of the grid.
      common/jambuu4/avp(4,-maxx:maxx,-maxx:maxx,-maxz:maxz)


c....end of jam3.inc

