      module wang_landau_module

      use utility_pack
      implicit none

      integer, allocatable :: visit_hist(:,:)
      real(8), allocatable :: dos_hist(:,:)
      real(8), allocatable :: dlambda(:)
c     maxvar is the maximum number of variable 'bins' to sample
c     in the WL simulation. For now this just means the maximum
c     number of guests to attempt.
c     Ideally this should be a function of the liquid density of
c     a gas, and the available pore space in a material.
      integer, parameter :: maxvar=20
      save visit_hist
      save dos_hist,dlambda
      contains 
      subroutine alloc_wl_arrays
     &(idnode,nhist,ntpguest)
c************************************************************************
c
c     allocate the histogram and other arrays necessary
c     to carry out a Wang-Landau simulation.
c     
c************************************************************************
      implicit none
      integer idnode,i,nhist,ntpguest
      integer, parameter :: nwl=3
      integer, dimension(nwl) :: fail

      do i=1,nwl
        fail(i) = 0
      enddo
      allocate(visit_hist(maxvar+1,nhist), stat=fail(1))
      allocate(dos_hist(maxvar+1,nhist), stat=fail(2))
      allocate(dlambda(ntpguest), stat=fail(3))

      do i=1,nwl
        if(fail(i).gt.0)then
            if(idnode.eq.0)write(nrite,'(10i5)')fail(i)
            call error(idnode, 2317)
        endif
      enddo
      do i=1,nhist
        dos_hist(:,i)=1.d0
        visit_hist(:,i)=0
      enddo
      do i=1,ntpguest
        dlambda(i)=0.d0
      enddo
      
      end subroutine alloc_wl_arrays

      end module wang_landau_module
