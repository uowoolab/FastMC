      module ewald_module

c***********************************************************************
c     
c     
c***********************************************************************
      use utility_pack

      implicit none
c     bunch of ewald arrays and parameters
      integer mxebuf
      real(8), allocatable :: ckc(:),cks(:),clm(:),slm(:)
      real(8), allocatable :: ckcsum(:,:),ckssum(:,:)
      real(8), allocatable :: ckcsorig(:,:),ckssorig(:,:) 
      real(8), allocatable :: ckcsnew(:,:),ckssnew(:,:)
      real(8), allocatable :: qfix_mol(:),qfix_molorig(:)
      real(8), allocatable :: elc(:,:),els(:,:)
      real(8), allocatable :: emc(:,:),ems(:,:)
      real(8), allocatable :: erc(:),fer(:)
      real(8), allocatable :: enc(:,:),ens(:,:)
      real(8), allocatable :: ewald1en(:)
      real(8), allocatable :: ewald2en(:)
      real(8), allocatable :: ewald3en(:)
      real(8), allocatable :: engsic(:)
      real(8), allocatable :: engsicprev(:)
   
      save mxebuf
      save ckc,cks,clm,slm
      save ckcsum,ckssum,ckcsorig,ckssorig
      save ckcsnew,ckssnew,qfix_mol,qfix_molorig
      save elc,els,emc,ems,erc,fer,enc,ens
      save ewald1en,ewald2en,ewald3en,engsic
      save engsicprev
      contains
      
      subroutine alloc_ewald_arrays
     &(idnode,maxmls,kmax1,kmax2,kmax3,rvdw,totatm,maxguest) 
      implicit none
      integer, parameter :: nv=26
      integer i,idnode,maxalloc,kmax1,kmax2,kmax3
      integer maxmls,totatm,maxguest
      integer, dimension(nv) :: fail
      real(8) rvdw
      
      maxalloc = totatm+maxguest
      mxegrd = max(1000, int(rvdw/0.01d0 + 0.5d0) + 4)
c     used to allocate ckcsum, ckssum etc.
      mxebuf = (2*kmax1+1)*(2*kmax2+1)*(2*kmax3+1)-1
c     N choose 2 + last entry is storage for total sum
      do i=1,nv
        fail(i) = 0
      enddo
      allocate(ckc(maxalloc),stat=fail(1))
      allocate(cks(maxalloc),stat=fail(2))
      allocate(ckcsum(maxmls+1,mxebuf),stat=fail(3))
      allocate(ckssum(maxmls+1,mxebuf),stat=fail(4))
      allocate(ckcsnew(maxmls+1,mxebuf),stat=fail(5))
      allocate(ckssnew(maxmls+1,mxebuf),stat=fail(6))
      allocate(ckcsorig(maxmls+1,mxebuf),stat=fail(7))
      allocate(ckssorig(maxmls+1,mxebuf),stat=fail(8))
      allocate(clm(maxalloc),stat=fail(9))
      allocate(slm(maxalloc),stat=fail(10))
      allocate(elc(maxalloc,0:1),stat=fail(11))
      allocate(els(maxalloc,0:1),stat=fail(12))
      allocate(emc(maxalloc,0:kmax2),stat=fail(13))
      allocate(ems(maxalloc,0:kmax2),stat=fail(14))
      allocate(enc(maxalloc,0:kmax3),stat=fail(15))
      allocate(ens(maxalloc,0:kmax3),stat=fail(16))
      allocate(erc(mxegrd),stat=fail(17))
      allocate(fer(mxegrd),stat=fail(18))
      allocate(ewald1en(maxmls+1),stat=fail(19))
      allocate(ewald2en(maxmls+1),stat=fail(20))
      allocate(ewald3en(maxmls+1),stat=fail(21))
      allocate(engsic(maxmls+9),stat=fail(23))
      allocate(engsicprev(maxmls+9),stat=fail(24))
      allocate(qfix_mol(maxmls+9),stat=fail(25))
      allocate(qfix_molorig(maxmls+9),stat=fail(26))

      do i=1,nv
        if(fail(i).gt.0)then
          if(idnode.eq.0)write(nrite,'(10i5)')fail(i)
          call error(idnode, 1002)
        endif
      enddo
      ewald1en(:) = 0.d0
      ewald2en(:) = 0.d0
      ewald3en(:) = 0.d0
      engsic(:) = 0.d0
      engsicprev(:) = 0.d0
      qfix_mol(:) = 0.d0
      qfix_molorig(:) = 0.d0
      end subroutine alloc_ewald_arrays

      end module ewald_module
