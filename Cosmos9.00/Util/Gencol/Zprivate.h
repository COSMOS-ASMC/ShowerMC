        integer maxn
         parameter ( maxn = 50000 )
        real*8  wzmin, wzmax, trackl, LabEquivE
         real*8  ke(maxn)
         integer indx(maxn)
        integer nevent
        integer outzero, outresul1
        integer  pjcode, pjsub, pjchg
         integer  tgcode, tgsub, tgchg
        integer  inpfileno, xyz
        real*8   pjpx, pjpy, pjpz
         real*8   tgpx, tgpy, tgpz
        real(8)::b
        character(8):: targetName
!        real  xpos, ypos, zpos  
        real(4)::xpos(3)     
        integer pjmnum, tgmnum  ! proj. target mass number
        type (ptcl):: pj, tg, plab
        common /Zotherc/ pj, tg, plab,  ke, wzmin, wzmax, 
     *    trackl, pjpx, pjpy, pjpz, tgpx, tgpy, tgpz,
!     *    xpos, ypos, zpos,  
     *    xpos, LabEquivE, b,
     *    indx, nevent, outzero, outresul1,
     *    pjmnum, tgmnum, inpfileno, xyz,
     *	  pjcode, pjsub, pjchg,
     *	  tgcode, tgsub, tgchg
        common /ZothercCha/ targetName



         



