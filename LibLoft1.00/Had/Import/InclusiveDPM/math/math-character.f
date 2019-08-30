      subroutine str2words(line, lmax0, ll, nword)
c Split words from a string given as line.
c lmax0 : possible maximum character number in line (input)
c ll    : words
c nword  : number of word
c
      character line*(*)
      character ll(100)*80

      lmax = lchend(line, lmax0)
      do l=1, lmax
         if(line(l:l).eq.',' .or. line(l:l).eq.';'
     &        .or. line(l:l).eq.':') then
            line(l:l) = ' '
         end if
      enddo

      i = 1
      j = 1
      nword = 0
      do 1 while (i.le.lmax .and. j.le.lmax)
         i=j
         do 10 while(line(i:i).eq.' ' .and. i.le.lmax)
            i=i+1
 10      continue

         j=i+1
         do 20 while(line(j:j).ne.' ' .and. j.le.lmax)
            j=j+1
 20      continue
         if(i.le.lmax) then
            nword = nword+1
            ll(nword) = line(i:j-1)
         end if
 1    continue
      end

      function numend(line, lmax)
      integer numend
      character line*(*)
      ll = lmax
      do 150 while (ichar(line(ll:ll)).lt.48
     &        .or. ichar(line(ll:ll)).gt.57)
         ll = ll-1
 150  continue
      numend = ll
      end
      
      function lchend(line, lmax)
      integer lchend
      character line*(*)
      ll=lmax
      do 150 while (ichar(line(ll:ll)).lt.33 
     &        .or. ichar(line(ll:ll)).eq.96
     &        .or. ichar(line(ll:ll)).gt.127)
         ll = ll-1
 150  continue
      lchend = ll
      end
      
      FUNCTION ISNUM(LINE)
      CHARACTER C(17)*1, LINE*(*)
      LOGICAL ISNUM
      DATA C/'0','1','2','3','4','5','6','7','8','9','.','+','-',' ',
     & 'E','D','e'/

      le = lchend(line,80)
      i = 1
      do 10 while(line(i:i).eq.' ' .and. i.lt.le) 
         i = i + 1
 10   continue

      do 18 j=1, 13
         IF(LINE(I:I).EQ.C(J)) THEN
            i = i + 1
            go to 19
         END IF
 18   continue
 19   continue

      do 20 while(line(i:i).ne.' ' .and. i.lt.le)
         DO 21 J=1, 17
            IF(LINE(I:I).EQ.C(J)) THEN
               i = i + 1
               go to 20
            END IF
 21      CONTINUE
         ISNUM=.FALSE.
         RETURN

 20   CONTINUE
c                At least 1st term is consist of above C(1..17).
      ISNUM=.true.
      RETURN
      END

