      module wang_landau_module

      use utility_pack
      implicit none

      integer, allocatable :: weight_hist(:)

      save weight_hist

      contains 
      subroutine alloc_wl_arrays
     &(idnode, nhist)
c************************************************************************
c
c     allocate the histogram and other arrays necessary
c     to carry out a Wang-Landau simulation.
c     
c************************************************************************
      implicit none
      integer idnode,i,nhist
      integer, parameter :: nwl=1
      integer, dimension(nwl) :: fail

      do i=1,nwl
        fail(i) = 0
      enddo
      allocate(weight_hist(nhist), stat=fail(1))

      do i=1,nwl
        if(fail(i).gt.0)then
            if(idnode.eq.0)write(nrite,'(10i5)')fail(i)
            call error(idnode, 2317)
        endif
      enddo
      end subroutine alloc_wl_arrays

      end module wang_landau_module
