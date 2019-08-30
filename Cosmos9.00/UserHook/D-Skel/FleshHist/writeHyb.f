
      if(ObserveAS) then
#if defined (DOMPI)
         msg = dir(1:lengdir)//"/"//execid(1:lengid)//".hyb"
#define NgXXX NgA
#define NeXXX NeA
#define NmuXXX NmuA
#define NhadXXX NhadA
#define SumElossXXX SumElossA
#else         
         msg = dir(1:lengdir)//"/"//execid(1:lengid)//
     *        "-@."//numb(1:lengn)//".hyb"
#define NgXXX Ng
#define NeXXX Ne
#define NmuXXX Nmu
#define NhadXXX Nhad
#define SumElossXXX SumEloss
#endif
         call copenfw2(fnoB, msg, 1, icon)
         if(icon .gt. 1) then
            write(0,*) ' icon=', icon
            call cerrorMsg(msg, 1)
            call cerrorMsg('could not be opened', 0)
         endif
         cog = 0.
         sumne = 0.

         do i = 1, NoOfASSites

            ASObsSites(i).esize = ASObsSites(i).esize* enhance

            if(i .gt. 1 .and. i  .lt. NoOfASSites ) then
               dd =(ASDepthList(i+1) - ASDepthList(i-1))/2.0
            elseif(i .eq. 1) then
               dd =(ASDepthList(2) - ASDepthList(1))
            else
               dd =(ASDepthList(NoOfASSites) -
     *              ASDepthList(NoOfASSites-1))
            endif
            cog = cog + ASObsSites(i).esize*dd*ASDepthList(i)
            sumne= sumne +ASObsSites(i).esize*dd
         enddo
!          0.1 is for g/cm2
         cog = cog*0.1/sumne

         cog2 = 0.
         sumne = 0.
         do i = 1, NoOfASSites
            if( ASObsSites(i).age .gt.
     *          (2.0-ASObsSites(NoOfASSites).age))  then
               if(i .gt. 1 .and. i  .lt. NoOfASSites ) then
                  dd =( ASDepthList(i+1) - ASDepthList(i-1))/2.0
               elseif(i .eq. 1) then
                  dd =(ASDepthList(2) - ASDepthList(1))
               else
                  dd =(ASDepthList(NoOfASSites) -
     *              ASDepthList(NoOfASSites-1))
               endif
               dd = dd
               cog2 = cog2 + ASObsSites(i).esize*ASDepthList(i)*dd
               sumne= sumne +ASObsSites(i).esize*dd
            endif
         enddo
         if(sumne .gt. 0.) then
            cog2 = cog2*0.1/sumne
         else
!              too deep penetration
            cog2 = ASDepthList(NoOfASSites)*0.1
         endif


         if(fnoB .ge. 0 )  then
            write(fnoB,
     *      '("h ", i4,  3i3, 1pE11.3, 0p 3f11.7, f7.2, 2f7.0)')
     *      EventNo,  inci.p.code,
     *      inci.p.subcode, inci.p.charge,
     *      inci.p.fm.e, -angle.r(1), -angle.r(2), -angle.r(3),
     *      firstz, cog, cog2
         else
            write(*,
     *      '("h ", i4,  3i3, 1pE11.3, 0p 3f11.7, f7.2, 2f7.0)')
     *      EventNo,  inci.p.code,
     *      inci.p.subcode, inci.p.charge,
     *      inci.p.fm.e, -angle.r(1), -angle.r(2), -angle.r(3),
     *      firstz, cog, cog2
         endif

         do i = 1, NoOfASSites 
            if(fnoB .ge. 0) then
               write(fnoB, '("t ", i3, 2f7.1,  2f6.3, 1p6E11.3)')
     *          i, 
     *          ASDepthList(i)*0.1,  ASObsSites(i).mu,
     *          ASObsSites(i).age,   ASDepthList(i)*0.1/cog2, 
     *          NgXXX(i), NeXXX(i),  NmuXXX(i), NhadXXX(i),
     *          ASObsSites(i).esize, SumElossXXX(i)  
            else
               write(*, '("t ", i3, 2f7.1,  2f6.3, 1p6E11.3)')
     *          i, 
     *          ASDepthList(i)*0.1,  ASObsSites(i).mu,
     *          ASObsSites(i).age,   ASDepthList(i)*0.1/cog2, 
     *          NgXXX(i), NeXXX(i), NmuXXX(i), NhadXXX(i),
     *          ASObsSites(i).esize, SumElossXXX(i)  
            endif
         enddo
         if(fnoB .gt. 0 ) then
            write(fnoB,*)
            close(fnoB)
         else
            write(*,*)
         endif
      endif
