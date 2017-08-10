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
      real(8), allocatable :: chgsum_mol(:),chgsum_molorig(:)
      real(8), allocatable :: elc(:,:),els(:,:)
      real(8), allocatable :: emc(:,:),ems(:,:)
      real(8), allocatable :: erc(:),fer(:)
      real(8), allocatable :: enc(:,:),ens(:,:)
      real(8), allocatable :: ewald1en(:),ewald1entmp(:)
      real(8), allocatable :: ewald2en(:),ewald2entmp(:)
      real(8), allocatable :: ewald3en(:)
      real(8), allocatable :: engsic(:)
      real(8), allocatable :: engsicorig(:)
   
      save mxebuf
      save ckc,cks,clm,slm
      save ckcsum,ckssum,ckcsorig,ckssorig
      save ckcsnew,ckssnew,chgsum_mol,chgsum_molorig
      save elc,els,emc,ems,erc,fer,enc,ens
      save ewald1en,ewald2en,ewald3en,engsic
      save engsicorig,ewald1entmp,ewald2entmp
      contains
      
      subroutine alloc_ewald_arrays
     &(idnode,maxmls,kmax1,kmax2,kmax3,rvdw,totatm,maxguest) 
      implicit none
      integer, parameter :: nv=27
      integer i,idnode,maxalloc,kmax1,kmax2,kmax3
      integer maxmls,totatm,maxguest,mxcmls,ik,kk
      integer, dimension(nv) :: fail
      real(8) rvdw
      mxcmls=(maxmls)*(maxmls-1)/2+maxmls
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
      allocate(ewald1en(mxcmls),stat=fail(19))
      allocate(ewald1entmp(mxcmls),stat=fail(20))
      allocate(ewald2en(mxcmls),stat=fail(21))
      allocate(ewald2entmp(mxcmls),stat=fail(22))
      allocate(ewald3en(maxmls),stat=fail(23))
      allocate(engsic(mxcmls),stat=fail(24))
      allocate(engsicorig(mxcmls),stat=fail(25))
      allocate(chgsum_mol(maxmls),stat=fail(26))
      allocate(chgsum_molorig(maxmls),stat=fail(27))

      do i=1,nv
        if(fail(i).gt.0)then
          if(idnode.eq.0)write(nrite,'(10i5)')fail(i)
          call error(idnode, 1002)
        endif
      enddo
      do ik=1,mxcmls
        ewald1en(ik)=0.d0
        ewald1entmp(ik)=0.d0
        ewald2en(ik) = 0.d0
        ewald2entmp(ik)=0.d0
        engsic(ik) = 0.d0
        engsicorig(ik) = 0.d0
      enddo
      do ik=1,maxmls 
        ewald3en(ik) = 0.d0
        chgsum_mol(ik) = 0.d0
        chgsum_molorig(ik) = 0.d0
      enddo
      do kk=1,mxebuf
        do ik=1,maxmls
          ckcsum(ik,kk)=0.d0
          ckssum(ik,kk)=0.d0
          ckcsnew(ik,kk)=0.d0
          ckssnew(ik,kk)=0.d0
          ckcsorig(ik,kk)=0.d0
          ckssorig(ik,kk)=0.d0
        enddo
        ckcsum(maxmls+1,kk)=0.d0
        ckssum(maxmls+1,kk)=0.d0
        ckcsnew(maxmls+1,kk)=0.d0
        ckssnew(maxmls+1,kk)=0.d0
        ckcsorig(maxmls+1,kk)=0.d0
        ckssorig(maxmls+1,kk)=0.d0
      enddo
      do ik=1,maxalloc
        elc(ik,0:1)=0.d0
        els(ik,0:1)=0.d0
        ckc(ik)=0.d0
        cks(ik)=0.d0
        clm(ik)=0.d0
        slm(ik)=0.d0
        do kk=1,kmax2
          emc(ik,kk)=0.d0
          ems(ik,kk)=0.d0
        enddo
        do kk=1,kmax3
          enc(ik,kk)=0.d0
          ens(ik,kk)=0.d0
        enddo
      enddo
      do ik=1,mxegrd
        erc(ik)=0.d0
        fer(ik)=0.d0
      enddo
      end subroutine alloc_ewald_arrays

      end module ewald_module
