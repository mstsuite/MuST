c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fstop(name)
c     =================================================================
c
c     this routine will stop the program when an error has occured
c     in a subroutine
c
      character*32 name,nam2
c
      data nam2/'null'/
c
      if(nam2.eq.name) then
         goto 10
      endif
c
      write(6,'(///)')
      write(6,'(3x,61(''*''))')
      write(6,'(3x,''* fstop called with index '',a32,t64,''*'')')name
      write(6,'(3x,61(''*''))')
c
 10   continue
cJL      call kill(0,9)
      stop
      end
