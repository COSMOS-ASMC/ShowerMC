      subroutine  epZ2ChemSymb(Z, symbol)
      implicit none
      integer Z           ! input   1 to 100
      character*2  symbol  ! output  Ag H etc
!       
      character*2  sym(100)
      data sym/
     1 "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
     2 "Na","Mg", "Al", "Si", "P", "S", "Cl","Ar","K", "Ca",
     3 "Sc","Ti", "V",  "Cr", "Mn","Fe","Co","Ni","Cu","Zn",
     4 "Ga","Ge", "As", "Se", "Br","Kr","Rb","Sr","Y", "Zr",
     5 "Nb","Mo","Tc", "Ru","Rh","Pd", "Ag","Cd","In", "Sn",
     6 "Sb","Te","I", "Xe", "Cs","Ba", "La", "Ce","Pr","Nd",
     7 "Pm","Sm","Eu","Gd", "Tb","Dy", "Ho", "Er","Tm","Yb",
     8 "Lu","Hf","Ta","W", "Re", "Os", "Ir", "Pt","Au","Hg",
     9 "Tl","Pb","Bi","Po","At", "Rn", "Fr", "Ra","Ac","Th", 
     A "Pa","U", "Np","Pu","Am", "Cm", "Bk", "Cf","Es","Fm"/
      save sym
      symbol = sym(Z)
      end
