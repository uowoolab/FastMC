      module wang_landau_module

      use utility_pack
      implicit none

      integer, allocatable :: visit_hist(:,:)
      integer, allocatable :: tmat_vis(:,:)
      real(8), allocatable :: dos_hist(:,:)
      real(8), allocatable :: energy_hist(:,:,:)
      real(8), allocatable :: dlambda(:)
      real(8), allocatable :: tmat_c(:,:)
      !real(8), parameter :: mpdpw = 14.449439791871d0
      save visit_hist
      save dos_hist,dlambda,energy_hist
      save tmat_c,tmat_vis
      contains 
      subroutine alloc_wl_arrays
     &(idnode,nhist,ntpguest,maxn,ebins)
c************************************************************************
c
c     allocate the histogram and other arrays necessary
c     to carry out a Wang-Landau simulation.
c     
c************************************************************************
      implicit none
      integer idnode,i,j,k,nhist,ntpguest,maxn,ebins
      integer, parameter :: nwl=6
      integer, dimension(nwl) :: fail

      do i=1,nwl
        fail(i) = 0
      enddo
      allocate(visit_hist(maxn+1,nhist), stat=fail(1))
      allocate(dos_hist(maxn+1,nhist), stat=fail(2))
      allocate(dlambda(ntpguest), stat=fail(3))
      allocate(tmat_c(maxn+1,maxn+1), stat=fail(4))
      allocate(tmat_vis(maxn+1,maxn+1), stat=fail(5))
      allocate(energy_hist(ebins,maxn+1,nhist), stat=fail(6))

      do i=1,nwl
        if(fail(i).gt.0)then
            if(idnode.eq.0)write(nrite,'(10i5)')fail(i)
            call error(idnode, 2317)
        endif
      enddo
      do i=1,nhist
        do j=1,maxn+1
          dos_hist(j,i)=0.d0
          visit_hist(j,i)=0
          ! redundancy in outer i loop here, but
          ! currently nhist=1, so doesn't matter
          do k=1,j
            tmat_c(i,j)=0.d0
            tmat_c(j,i)=0.d0
            tmat_vis(i,j)=0
            tmat_vis(j,i)=0
          enddo
          do k=1,ebins
            energy_hist(k,j,i)=0.d0
          enddo
        enddo
      enddo

      do i=1,ntpguest
        dlambda(i)=0.d0
      enddo
      
      end subroutine alloc_wl_arrays

      end module wang_landau_module
