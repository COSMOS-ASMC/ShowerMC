//job job c2s5091,class=a
//  exec fortclg,parm.fort='nos,num,gostmt'
//fort.sysinc dd dsn=c2g5100.cosmos.gem,disp=shr
//fort.sysin dd *
       include  'mufug.f'
/*
//lked.syslib dd dsn=c2s5091.#load.load,disp=shr
//  dd
//  dd
//  dd
//  dd dsn=c2g5600.lib.load,disp=shr
//go.ft05f001 dd *
    0/
   'c2g5100.rock.data(frejus)',55,23793,99999999/
   'c2s5091.#frj300.data', 'c2s5091.#gd2.data'/
//go.ft09f001 dd dsn[c2s5091.#cont.data,disp=shr
//
     above is for
     cont job information (0--> 1st, 1--> cont: if cont==>further
     lines not needed.
     rock profile dataset name, jobt,laste, nend /
     input and output dataset name /
