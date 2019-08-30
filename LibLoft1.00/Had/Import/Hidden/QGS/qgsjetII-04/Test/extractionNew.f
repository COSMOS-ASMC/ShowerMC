c check diffractive cross sections          !so060413-beg
9     jdiff=0                             !non-diffractive
      nqst=0
      nint=0
      if(nbpom.ne.0)then
       do i=1,nbpom
        nqst=nqst+nqs(i)
        nint=nint+npomin(i)
       enddo
      endif
      if((nwp.ne.0.or.nwt.ne.0).and.nqst.eq.0)then   !not elastic nor ND
       lqat=0
       do ip=1,ia(1)
        lqat=lqat+lqa(ip)
       enddo
       lqbt=0
       do it=1,ia(2)
        lqbt=lqbt+lqb(it)
       enddo
       iwpt=0
       do ip=1,ia(1)
        if(iwp(ip).eq.1)then
         iwpt=1
         goto 10
        elseif(iwp(ip).ge.2)then
         iwpt=2
        endif
       enddo
10     continue
       iwtt=0
       do it=1,ia(2)
        if(iwt(it).eq.1)then
         iwtt=1
         goto 11
        elseif(iwt(it).ge.2)then
         iwtt=2
        endif
       enddo
11     continue
       if(lqat.eq.0.and.lqbt.eq.0)then
        if(nbpom.eq.0.or.nint.eq.0)then
         if(iwpt.eq.2.and.iwtt.ne.2)then
          jdiff=6                         !SD(LM)-proj
         elseif(iwpt.ne.2.and.iwtt.eq.2)then
          jdiff=7                         !SD(LM)-targ
         elseif(iwpt.eq.2.and.iwtt.eq.2)then
          jdiff=8                         !DD(LM)
         else
          goto 14
         endif
        else
         if(iwpt.ne.2.and.iwtt.ne.2)then
          jdiff=9                         !CD(DPE)
         else
          jdiff=10                        !CD+LMD
         endif
        endif
       elseif(lqat.gt.0.and.lqbt.eq.0.and.iwtt.ne.2)then
        jdiff=1                          !SD(HM)-proj
       elseif(lqat.eq.0.and.lqbt.gt.0.and.iwpt.ne.2)then
        jdiff=2                          !SD(HM)-targ
       elseif(lqat.gt.0.and.lqbt.eq.0.and.iwtt.eq.2)then
        jdiff=3                          !DD(LHM)-proj
       elseif(lqat.eq.0.and.lqbt.gt.0.and.iwpt.eq.2)then
        jdiff=4                          !DD(LHM)-targ

       elseif(lqat.gt.0.and.lqbt.gt.0)then
        if(nbpom.eq.0)stop'problem with nbpom!!!'
        xrapmax(1)=1.d0
        xrapmin(1)=1.d0/scm
        do ibpom=1,nbpom
         if(npompr(ibpom).gt.0)then
          do i=1,npompr(ibpom)
           ip=ilpr(i,ibpom)
           lpom=lnpr(i,ibpom)
           xrapmax(1)=min(xrapmax(1),1.d0/xpompr(lpom,ip)/scm)
          enddo
         endif
         if(npomtg(ibpom).gt.0)then
          do i=1,npomtg(ibpom)
           it=iltg(i,ibpom)
           lpom=lntg(i,ibpom)
           xrapmin(1)=max(xrapmin(1),xpomtg(lpom,it))
          enddo
         endif
        enddo
        if(xrapmin(1).gt..999d0*xrapmax(1))goto 14
        nraps=1
12      if(nraps.gt.90)stop'nraps>90'
        do ibpom=1,nbpom
         if(npomin(ibpom).gt.0)then
          do i=1,npomin(ibpom)
           if(nraps.eq.1)then
            if(1.d0/scm/xpomin(i,ibpom).lt..999d0*xrapmax(1)
     *      .and.xpopin(i,ibpom).gt.1.001d0*xrapmin(1))then !rap-gaps changed
             if(1.d0/scm/xpomin(i,ibpom).lt.1.001d0*xrapmin(1)
     *       .and.xpopin(i,ibpom).gt..999d0*xrapmax(1))then !no rap-gap (filled)
               goto 14
             elseif(xpopin(i,ibpom).gt..999d0*xrapmax(1))then
              xrapmax(1)=1.d0/scm/xpomin(i,ibpom)
             elseif(1.d0/scm/xpomin(i,ibpom).lt.1.001d0*xrapmin(1))then
              xrapmin(1)=xpopin(i,ibpom)
             else
              xrapmin(2)=xrapmin(1)
              xrapmin(1)=xpopin(i,ibpom)
              xrapmax(2)=1.d0/scm/xpomin(i,ibpom)
              nraps=2
              goto 12
             endif
            endif
           else
            if(1.d0/scm/xpomin(i,ibpom).lt..999d0*xrapmax(1)
     *      .and.xpopin(i,ibpom).gt.1.001d0*xrapmin(nraps))then !rap-gaps changed
             if(1.d0/scm/xpomin(i,ibpom).lt.1.001d0*xrapmin(nraps)
     *       .and.xpopin(i,ibpom).gt..999d0*xrapmax(1))then !no rap-gaps (filled)
              goto 14
             else
              do irap=1,nraps
               if(xpopin(i,ibpom).gt..999d0*xrapmax(irap).and.1.d0/scm
     *         /xpomin(i,ibpom).lt.1.001d0*xrapmin(irap))then !gap filled
                if(irap.lt.nraps)then
                 do j=irap,nraps-1
                  xrapmax(j)=xrapmax(j+1)
                  xrapmin(j)=xrapmin(j+1)
                 enddo
                endif
                nraps=nraps-1
                goto 12
               elseif(xpopin(i,ibpom).gt..999d0*xrapmax(irap))then
                xrapmax(irap)=min(1.d0/scm/xpomin(i,ibpom)
     *          ,xrapmax(irap))
               elseif(1.d0/scm/xpomin(i,ibpom)
     *         .lt.1.001d0*xrapmin(irap))then
                xrapmin(irap)=max(xpopin(i,ibpom),xrapmin(irap))
               elseif(1.d0/scm/xpomin(i,ibpom).gt.xrapmin(irap)
     *         .and.xpopin(i,ibpom).lt.xrapmax(irap))then
                xrapmin(irap)=max(xpopin(i,ibpom),xrapmin(irap))
                if(irap.lt.nraps)then
                 do j=1,nraps-irap
                  xrapmax(nraps-j+2)=xrapmax(nraps-j+1)
                  xrapmin(nraps-j+2)=xrapmin(nraps-j+1)
                 enddo
                endif
                xrapmin(irap+1)=xrapmin(irap)
                xrapmin(irap)=xpopin(i,ibpom)
                xrapmax(irap+1)=1.d0/scm/xpomin(i,ibpom)
                nraps=nraps+1
                goto 12
               endif
              enddo                       !end of irap-loop
             endif
            endif
           endif
          enddo                           !end of npin-loop
         endif
        enddo                             !end of ibpom-loop
        jdiff=5                          !DD(HM)
       endif
      endif                              !end of diffr. check
14    bdiff=b

ctp define collision type
      typevt=0                      !no interaction
      if(nwp.gt.0.or.nwt.gt.0)then         !so060413-end
       if(jdiff.eq.0)then                                  !ND (no rap-gaps)
        typevt=1
       elseif(jdiff.eq.8.or.jdiff.eq.10.or.
     *       (jdiff.gt.2.and.jdiff.lt.6))then !DD + (CD+LMD)
        typevt=2
       elseif(jdiff.eq.1.or.jdiff.eq.6)then                  !SD pro
        typevt=4
       elseif(jdiff.eq.2.or.jdiff.eq.7)then                  !SD tar
        typevt=-4
       elseif(jdiff.eq.9)then                                !CD
        typevt=3
       else
        stop'problem with typevt!'
       endif
      endif


c form projectile spectator part
      if(debug.ge.1)write (moniou,227)
      nspec=0
      do ip=1,ia(1)
       if(iwp(ip).eq.0)then
        if(debug.ge.2)write (moniou,228)ip
        nspec=nspec+1
        do l=1,3
         xas(nspec,l)=xa(ip,l)
        enddo
       endif
      enddo

      nspect=0
      do it=1,ia(2)
       if(iwt(it).eq.0)nspect=nspect+1
      enddo

c inelastic interaction: energy sharing and particle production
      if(nwp.ne.0.or.nwt.ne.0)then
       if(ia(1).eq.nspec.or.ia(2).eq.nspect)stop'ia(1)=nspec!!!'
       if(debug.ge.1)write (moniou,229)

       call qgsha(nbpom,ncola,ncolb,iret)
       if(iret.ne.0)goto 1
       if(nsp.le.nsp0+2)then
        if(debug.ge.1)write (moniou,230)
        goto 1
       endif
      else                                 !no interaction
       if(debug.ge.1)write (moniou,231)
       goto 1
      endif
      if(debug.ge.1)write (moniou,232)nsp

c fragmentation of the projectile spectator part
      if(debug.ge.1)write (moniou,233)
      call qgfrgm(nspec,xas)
      if(debug.ge.1)write (moniou,234)nsf
      if(debug.ge.1)write (moniou,235)

201   format(2x,'qgconf - configuration of the interaction')
202   format(2x,'qgconf: impact parameter b=',e10.3,' fm')
203   format(2x,'qgconf: nuclear configurations')
204   format(2x,'qgconf: check connections')
205   format(2x,'qgconf: ',i3,'-th proj. nucleon may interact with '
     *,i3,'-th target nucleon')
206   format(2x,'qgconf: ',i3,'-th projectile nucleon may interact')
207   format(2x,'qgconf: ',i3,'-th target nucleon may interact')
208   format(2x,'qgconf: ',i3,'-th rejection,'
     *,' redo Pomeron configuration')
209   format(2x,'qgconf: too many rejections,'
     *,' redo nuclear configuartions')
210   format(2x,'qgconf: Pomeron configuration')
211   format(2x,'qgconf: check ',i3,'-th projectile nucleon')
212   format(2x,'qgconf: interaction with ',i3,'-th target nucleon?')
213   format(2x,'qgconf: eikonals - total: ',e10.3,2x,'single: ',e10.3)
214   format(2x,'qgconf: ',i4,'-th Pomeron block connected to ',i3
     *,'-th proj. nucleon and'/4x,i3,'-th targ. nucleon;'
     *,' number of element. processes in the block: ',i3)
215   format(2x,'qgconf: ',i3
     *,'-th process in the block is single cut Pomeron')
219   format(2x,'qgconf: configuration of multi-Pomeron vertexes')
220   format(2x,'qgconf: ',i3,'-th process in the block contains '
     *,i3,' interm. Pomerons')
221   format(2x,'qgconf: ',i3,'-th process in the block contains '
     *,i3,' proj. legs')
222   format(2x,'qgconf: ',i3,'-th process in the block contains '
     *,i3,' targ. legs')
223   format(2x,'qgconf: check diffraction for ',i3,'-th proj. nucleon')
224   format(2x,'qgconf: diffr. of ',i3,'-th proj. nucleon,'
     *,' recoil of ',i3,'-th targ. nucleon')
225   format(2x,'qgconf: check diffraction for ',i3,'-th targ. nucleon')
226   format(2x,'qgconf: diffr. of ',i3,'-th targ. nucleon,'
     *,' recoil of ',i3,'-th proj. nucleon')
227   format(2x,'qgconf: projectile spectator part')
228   format(2x,'qgconf: ',i3,'-th proj. nucleon stays idle')
229   format(2x,'qgconf: inelastic interaction: energy sharing'
     *,' and particle production')
230   format(2x,'qgconf: no particle produced - rejection')
231   format(2x,'qgconf: no interaction - rejection')
232   format(2x,'qgconf: ',i5,' particles have been produced')
233   format(2x,'qgconf: fragmentation of the proj. spectator part')
234   format(2x,'qgconf: ',i3,' proj. fragments have been produced')
235   format(2x,'qgconf - end')
      return
      end
