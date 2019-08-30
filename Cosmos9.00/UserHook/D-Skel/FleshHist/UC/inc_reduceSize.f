!  This is the same as inc_readandput.f but this dose not read file
!  data is in buf.
!            this part is used only by rank 0
!      to read .dat created by each rank and reduce it
!      and to store data as a single file
      limit(1) = REALLIMITg
      limit(2) = REALLIMITe
      limit(3) = REALLIMITmu
      limit(4) = REALLIMITh

      do k = 1,  ansites
         do j = 1, 4
            do l = 1,nfai
               do i = 1, nrbin
                  nrfaiRec(i, l, j, k)=0
               enddo
            enddo
         enddo
      enddo
!         single file to store combined all .dat
!      msg=' '
!      msg =dir(1:lengdir)//"/"//execid(1:lengid)//".dat"
!///////////
!      write(0,*) ' asking main file opening'
!///////
!      call copenfw2(fnodat+1, msg, 1, icon)
!      if(icon .gt. 1) then
!         write(0,*) ' icon=', icon
!         call cerrorMsg(msg, 1)
!         call cerrorMsg('could not be opened', 0)
!      endif

      inbuf = 0
         do i = 1, bufc
            ldep = buf(i).ldep
            depidx = w2il(ldep)
            faiidx=  buf(i).faiidx
            ridx = buf(i).ridx
            codex = buf(i).code
            codex = min(codex, 4)
            wgt = buf(i).wgt
            if( nrfaiRecA(ridx, faiidx, codex, depidx) .gt.
     *           limit(codex) ) then
!                       accept with this prob.
               prob = limit(codex)/
     *              nrfaiRecA(ridx, faiidx, codex, depidx)

            else
               prob = 1.0
            endif
 
            if( prob .gt. 1.) then
               wwgt = wgt
               accept = .true.
            else
               prob = prob * wgt
               if(prob .gt. 1.) then
                  accept = .true.
                  wwgt = prob
               else
                  call rndc(u)
                  if(u .lt. prob) then
                     accept = .true.
                     wwgt= 1.
                  else
                     accept = .false.
                  endif
               endif
            endif
            if(accept) then
               nrfaiRec(ridx, faiidx, codex, depidx)=
     *              nrfaiRec(ridx, faiidx, codex, depidx) + wwgt
               inbuf = inbuf + 1
               buf(inbuf) = buf(i)
!               if(KeepWeight) then
!                  write(fnodat+1,
!     *  '(6i3, 1pE11.3, 0p,f6.1,1p2E11.3,0p, 2f8.4,f10.6,1pE11.3)')
!     *                 buf(i).ldep,  buf(i).code,  buf(i).subcode,
!     *                 buf(i).charge, buf(i).ridx, buf(i).faiidx,
!     *                 buf(i).rinmu, buf(i).fai, buf(i).Ek,
!     *                 buf(i).t, buf(i).wx, buf(i).wy, buf(i).wz,
!     *                 wwgt
!                  else
!                     write(fnodat+1,
!     *     '(6i3, 1pE11.3, 0p,f6.1,1p2E11.3,0p, 2f8.4,f10.6)')
!     *                 buf(i).ldep,  buf(i).code,  buf(i).subcode,
!     *                 buf(i).charge, buf(i).ridx, buf(i).faiidx,
!     *                 buf(i).rinmu, buf(i).fai, buf(i).Ek,
!     *                 buf(i).t, buf(i).wx, buf(i).wy, buf(i).wz
!                  endif
!               endif
            enddo
