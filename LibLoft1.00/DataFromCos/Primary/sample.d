# This is  
# an example table to show how you can specify primary particles.
# Any commnets may be placed bfore the line "#----..." below.
# The main body of the table consists of a number of blocks, each
# of which shows infomation on one type of primary.  
# Each block consists of a header line showing
# 1) the primary type, 
# 2) energy unit(MeV GeV etc),
# 3) energy type(particle energy, energy/neucleon, kinetic energy 
#    or total energy or moemtunm), 
# 4) differetial('d') or integral spectrum ('i')
#   ('d' is better)
# 5) power of the flux normalization factor.
#   (table data gives the flux value itself, or flux*E**2.0 etc)
# 6) if you want to put higher Emin than the one given in
#    the table, put a value for Emin here.
#    This enables the user to use the same table with diff. 
#    Emin easily. This is optional. If not listed, table min. is used.
#   This feature can be used only for differential spectrum.
#
# For example, the first line would look like
#
#  'CNO' 'MeV' 'KE/n'  'd' 1.5  /
#
# indicating the primary type is the Carbon/Nitrogen/Oxigen group,
# energy unit is in MeV,  enegy is kinetic energy per nucleon,
# spectrum is given in a differential form, and the flux value 
# is normalized as  dI/dE * E**1.5. 
# Or, If you give,
#
#    'CNO' 'MeV' 'KE/n'  'd' 1.5  100. /
#
#  The last 100. is to specify the minimum energy (in MeV, Kinetic 
#  energy / nucleon) to be sampled, though the table might contain
#  the energy lower than 100 MeV.
#  You can further give
#    'CNO' 'MeV' 'KE/n'  'd' 1.5  100. 1000. /
#
#  In this case the last 1000 is to specify the maximum energy
#  so that no event will be sampled over 1000. 
#
# The data is case insensitive so that you may write 'mev' instead
#             ~~~~~~~~~~~~~~~~
# of 'MeV'.  The available data types are listed below.
#
#  From the 2nd line of each block follow lines consisiting of
# 2 columns.  To be able to express a complex energy spectrum,
#  (If more than two, only first two are recognized).
# we  approximate a given primary spectrum by a number 
# of segment lines on a log-log scale graph.  The columns 
# express these segments. The maximum number of segments is 60.
#  (You can change it by chaging MAX_SEGMENTS in Zmaxdef.h).
# The first column is the left side  energy of a segment and
# the second column is the  flux there.
# The lines should continue in the ascending order of energy.
# To show the end of a block, you have to place two zero data.
#  **** next is obsolete ****
# The  spectrum is assumed to extend to infinity with a power 
# given in the last segment data.
# ******************

#  For example, to express a dI/dE = E**-3.0  specturm up to
#  100 PeV 
#   'P'  'PeV'  'E/P'  'd'  0.  /
#     10.0        1.
#     100.0       1.e-3
#      0           0
#  More complex data will be seen in the example data below.
#
# To express  mono energetic gamma and electrons, one may give
# as follows.
#
#     'gamma'  'TeV'  'E/P'  'd'  0  /
#        1.        2.
#        0          0
#     'e'      'GeV'  'E/P'   'd'  0  /
#        1000.     30.
#         0         0
#
# This case shows monoenergitic gamma rays with energy 1 TeV and
# electrons with 1000 GeV, with the 
# flux ratio (integral=differential), 2:30.
#
# Note: If two or more primaries are given and 
#      their minimum enegies are different, we assum that
#      the flux below each given minimum energy is absent.
#       
#
# The available symbols for primary type:
#         If multiple entries are shown, any one
#         can be used. As for the anti-particle,
#         one may put "~" after the particle 
#         symbol, say, 'e~' is the same as 'positron'.
#	  One may suspect why, say,  'rho' is listed
#	  as a primary.  This is to help debuggin the
#         program.
# From version 6, you can specify arbitrary isotopes
#( for  limited interaction models).
#
# 'gamma'  'photon'            gamma rays
# 'e'  'e-' 'electron'         electrons (not positrons)
# 'e+' 'positron'              positrons
# 'mu-'                        negative muons
# 'mu+'                        positive muons
# 'pi+'                        positive pions
# 'pi-'                        negative pions
# 'pi0'                        neutral pions
# 'k+'                         positive kaons
# 'k-'                         negative kaons
# 'k0l'                        k0 long
# 'k0s'                        k0 short
# 'p' 'proton'                 protons
# 'n' 'neutron'                neutrons
# 'rho'                        rho mesons
# 'omega'                      omega mesons
# 'bomega'                     Omega baryon
# 'phi'                        phi mesons
# 'D+'                         D +
# 'D-'                         D -
# 'D0'                         D 0
# 'DD~'                        D D~ pairs 
# 'nn~'                        n n~ pairs (pp~, nn~, pn~, p~n)
# 'neu_e'                      electron type neutrino (not anti)
# 'neu_mu'                     muon type neutrino (not anti)
# 'deuteron'                   deuteron
# 'alfa' 'He' 'alpha'          helium
# 'L' 'LiBeB'                  Light heavy. Li/Be/B group
# 'M' 'CNO'                    Medium heavy. C/N/O group
# 'H' 'NaMgSi'                 Heavy. Na/Mg/Si group
# 'VH' 'SClAr'                 Vergy heavy.  A/Cl/Ar group
# 'Fe' 'Iron'                  Iron group
# 
# 'iso 3 2'                    General isotope.  A=3 and Z=2 i.e.,
#                              3He. if 'iso 2 1', deuteron.   
#                              Of course, 'iso 4 2' can be used
#                              as 4He.
# Energy Unit. (Cosmos uses  GeV internally)
# (see also 'rig' for Energy type).
#
#   'eV'                  electron volt
#   'MeV'                10**6 eV
#   'GeV'                10**9 eV
#   'TeV'                10**12 eV
#   'PeV'                10**15 eV
#   'EeV'                10**18 eV
#
# Energy type.
# 
#    'E/P'  'E'            energy per particle
#   'KE/P'  'KE'           kinetic energy per particle
#    'E/n'                 energy per neucleon. If it is,
#                          say, for 'gamma', the same
#                          as E/P.
#    'KE/n'                kinetic erergy per neucleon. 
#    'p/P'   'p'           momentum per particle
#    'p/n'                 momentum per neucleon
#  
#    'rig'		   rigidity.  In this case, energy unit,
#                          say, 'GeV' is interpreted as 'GV'.
#
# The next is an example of a complex composition at low energy
#-----------------------------------------------------------
  'p'   'GeV'     'KE/n'  'd'     0   /
	0.1     1.2
        0.2     1.5
        .3      1.7
        .4      1.9
	.5	1.93
	.6	1.9
	.8	1.8
	1.5	1.5
	2.	1.25
	3.	.8
	4.	.55
	10.	.1
	20.	.02
	100.	2.8e-4
	0	0
  'He'  'GeV'    'KE/n'  'd'      0  /
	.1	.7
	.2	1.
	.4	1.2
	.6	1.25
	.8	1.2
	1.	1.15
	2.	.7
	5.	0.35
	10.	0.065
	30.	.008
	100.	2.e-4
         0        0
  'CNO' 'GeV'    'KE/n'  'd'      0 /
         .1	.013
	.2	.28
	.3	.4
	.5	.65
	.8	.8
	1.	.85
	1.3	.88
	2.0	.75
	4. 	.35
	6.	.2
	10.	.07
	20.	.012
         0        0
  'gamma'   'MeV'  'E'  'd'     1.  /
         100.      .1
         300.      .06
         500.      .04
         1000.     .01
         5000.     .003
         0         0
  'p~'    'Mev'   'E'    'i'    0  /
      1.      1.
      10.      .1
     100.     .005
	0 0
  'rho'  'TeV'   'E'     'd'     0 /
       100.,    1.
	0	0
