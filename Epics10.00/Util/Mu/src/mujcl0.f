//job job c2s5001
//*   to make library for mufolw test.  library is made in #load.load
//  exec fortcl,
//   parm.fort='nos,gostmt,num,name',
//*  parm.fort='nos,gostmt,num,name,debug(argchk,undef)',
//   parm.lked='ncal'
//fort.sysinc dd dsn=c2g5100.cosmos.gem,disp=shr
//fort.sysin dd *
       include  'mubrem.f'
       include  'mucset.f'
       include  'mudedx.f'
       include  'mufolw.f'
       include  'muhadr.f'
       include  'mumin.f'
       include  'mupair.f'
       include  'musub.f'
//*
//lked.syslmod dd dsn=c2s5001.#load.load,disp=shr
//
