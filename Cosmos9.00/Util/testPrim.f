      program testPrim
      implicit none
      character*80 fn
      call cerrorMsg(
     * 'Enter path to the primary data file(e.g; ../Data/prim.d)', 1)
      read(*, '(a)') fn
      call ciniSPrim(fn)
      call cprintPrim(5)
      end
