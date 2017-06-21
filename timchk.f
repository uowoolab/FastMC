      subroutine timchk(ktim,time)

c***********************************************************************
c     
c     timing routine for time elapsed in seconds
c
c***********************************************************************
      implicit none

      logical init
      character*12 dat,tim,zon
      integer idnode,mynode,ktim,day
      real(8) time,told,tsum,tnow
      integer info(8)

      save init,idnode,told,tsum,day

      data init/.true./

   10 format(/,' time elapsed since job start = ',f15.3,' seconds',/)

      call date_and_time(dat,tim,zon,info)
      
      if(init)then

         tsum=0.d0
         time=0.d0
         day=info(3)
         idnode=mynode()
         told=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))
         init=.false.

      else 

         tnow=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))
         if(day.ne.info(3))then
           told=told-86400.d0
           day=info(3)
         endif
         tsum=tsum+tnow-told
         told=tnow
         time=tsum

      endif

      if(ktim.gt.0.and.idnode.eq.0) write(nrite,10)time

      return
      end
