#  3   symbols.                  code, subcode, charge   
#----------------------------------------------------------
 'gamma'  'photon'  ' '         
 'e'  'e-' 'electron'         electrons (not positrons)
 'e+' 'positron'              positrons
 'mu-'                        negative muons
 'mu+'                        positive muons
 'pi+'                        positive pions
 'pi-'                        negative pions
 'pi0'                        neutral pions
 'k+'                         positive kaons
 'k-'                         negative kaons
 'k0l'                        k0 long
 'k0s'                        k0 short
 'p' 'proton'                 protons
 'n' 'neutron'                neutrons
 'rho'                        rho mesons
 'omega'                      omega mesons
 'phi'                        phi mesons
 'D+'                         D +
 'D-'                         D -
 'D0'                         D 0
 'DD~'                        D D~ pairs 
 'nn~'                        n n~ pairs (pp~, nn~, pn~, p~n)
 'neu_e'                      electron type neutrino (not anti)
 'neu_mu'                     muon type neutrino (not anti)
 'deuteron'                   deuteron
 'alfa' 'He' 'alpha'          helium
 'L' 'LiBeB'                  Light heavy. Li/Be/B group
 'M' 'CNO'                    Medium heavy. C/N/O group
 'H' 'NaMgSi'                 Heavy. Na/Mg/Si group
 'VH' 'SClAr'                 Vergy heavy.  A/Cl/Ar group
 'Fe' 'Iron'                  Iron group
 
Energy Unit. (Cosmos uses  GeV internally)

   'eV'                  electron volt
   'MeV'                10**6 eV
   'GeV'                10**9 eV
   'TeV'                10**12 eV
   'PeV'                10**15 eV
   'EeV'                10**18 eV

Energy type.
 
    'E/P'  'E'            energy per particle
   'KE/P'  'KE'           kinetic energy per particle
    'E/n'                 energy per neucleon. If it is
                          is, say, for 'gamma', the same
                          as E/P.
    'KE/n'                kinetic erergy per neucleon. 
    'p/P'   'p'           momentum per particle
    'p/n'                 momentum per neucleon

The next is an example of a complex composition at low energy
-----------------------------------------------------------
  'p'   'GeV'     'KE/n'  'd'     0  /
        0.5      1.
        1.0      .5
        2.0      .1
        5.0      .03
        0         0
  'He'  'GeV'    'KE/n'  'd'      0  /
        0.5      .2
        1.0      .08
        5.0      .01
         0        0
  'CNO' 'GeV'    'KE/n'  'd'      0   /
         .3      .1
         1.       .03
         10.       .005
         0        0
  'gamma'   'MeV'  'E'  'd'     1.5    /
         100.       1.
         300.      .8
         500.      1.0
         1000.     .9
         5000.     .8
         0         0
