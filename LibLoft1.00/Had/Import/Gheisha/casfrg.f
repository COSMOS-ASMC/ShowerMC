*cmz :  3.14/16 13/03/89  14.48.45  by  nick van eijndhoven (cern)
*-- author :
      subroutine casfrg(nucflg,int,nfl)
c
c *** cascade of heavy fragments ***
c *** nve 11-may-1988 cern geneva ***
c
c origin : h.fesefeldt (02-dec-1986)
c
c --- nucflg is a flag to denote the nucrec action ---
c nucflg = 0 ==> no action taken by nucrec
c          1 ==> action taken by nucrec
      nucflg=1
      call nucrec(nopt,2)
      if (nopt .ne. 0) go to 9999
c
      nucflg=0
      call coscat
c
 9999 continue
      return
      end
