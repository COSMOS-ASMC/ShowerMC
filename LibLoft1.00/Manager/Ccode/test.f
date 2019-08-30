      integer leng, klena
      character*1  NULL
      character*128 path
      integer kgetenv
      NULL = char(0)
      leng = kgetenv("COSMOSTOP"//NULL, path)
      path = path(1:leng)
      write(*,*) klena(path),' ', path(1:klena(path))
      end
