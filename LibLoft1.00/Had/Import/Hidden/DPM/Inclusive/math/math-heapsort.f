      SUBROUTINE HEAP_SORT(N,RA)
      implicit doubleprecision (a-h,o-z)
c Sorts anarray ra of length N into ascending numerical order using the
c heapsort algorithm. N is input: RA is replaced on output by its sorted 
c rearrangement
      DIMENSION RA(N)
      L=N/2+1
      IR=N
c The index L will be decremented from the inital value down to 1 during
c the "hiring" (heap creation) phase. Once it reaches 1, the index IR 
c will be decremented from its inital value down to 1 during the 
c "retirement-and-promotion" phase.
 10   CONTINUE
      IF(L.GT.1)THEN
         L=L-1
         RRA=RA(L)
      ELSE                      ! In retirement-and-promotion phase, 
         RRA=RA(IR)             ! Clear a space at end of array
         RA(IR)=RA(1)           ! Retire the top of the heap into it
         IR=IR-1                ! Decrease the size of the corporation.
         IF(IR.EQ.1)THEN        ! Done with the last promotion.
            RA(1)=RRA           ! The lease competent worker af all
            RETURN        
         ENDIF
      ENDIF
      I=L        ! Wheter we are in the hiring phase or promotion phase,
      J=L+L      ! we have set up to sift down element BRA to its proper
                 ! level.
 20   IF(J.LE.IR)THEN
         IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1 !compare to the better undrling
         ENDIF
         IF(RRA.LT.RA(J))THEN  ! demote BRA
            RA(I)=RA(J)
            I=J
            J=J+J
         ELSE    ! This is BRS's level, set J to termination of the
                 ! shift down.
            J=IR+1
         ENDIF
         GO TO 20
      ENDIF
      RA(I)=RRA                 ! Put BRA into its slot
      GO TO 10
      END
