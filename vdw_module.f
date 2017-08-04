      module vdw_module

c***********************************************************************
c     
c     module for defining van der waald potential arrays
c     
c***********************************************************************
      use utility_pack

      implicit none

      integer, allocatable :: ltpvdw(:),lstvdw(:)
      integer, allocatable :: numfrz_mol(:,:),numtyp_mol(:,:) 
      integer, allocatable :: numfrz_gstmol(:,:),numtyp_gstmol(:,:) 
      real(8), allocatable :: vvv(:,:),prmvdw(:,:)
      real(8), allocatable :: steadd(:),vdwen(:)
      real(8), allocatable :: elrc_mol(:),origelrc_mol(:)
      real(8), allocatable :: delrc_mol(:),delrc_mol0(:)
      real(8) elrc,origelrc
c     maximum number of vdw parameters

      integer, parameter :: mxpvdw=5

      save ltpvdw,lstvdw,prmvdw,vvv,vdwen
      save elrc,origelrc,numfrz_mol,numtyp_mol
      save elrc_mol,origelrc_mol
      save delrc_mol,delrc_mol0
      save numtyp_gstmol,numfrz_gstmol
      contains
      
      subroutine alloc_vdw_arrays(idnode,maxvdw,maxmls,mxatyp)
      implicit none
      integer, parameter :: nv=14
      integer maxvdw,i,idnode,maxmls,mxatyp,mxcombo
      integer, dimension(nv) :: fail

      mxcombo=maxmls*(maxmls-1)/2+maxmls
      do i=1,nv
        fail(i) = 0
      enddo
      allocate (ltpvdw(maxvdw),stat=fail(1))
      allocate (lstvdw(maxvdw),stat=fail(2))
      allocate (steadd(maxvdw),stat=fail(3))
      allocate (prmvdw(maxvdw,mxpvdw),stat=fail(4))
      allocate (vvv(mxegrd,maxvdw),stat=fail(5))
      allocate (vdwen(mxcombo+1),stat=fail(6))
      allocate (elrc_mol(mxcombo),stat=fail(7))
      allocate (origelrc_mol(mxcombo),stat=fail(8))
      allocate (numfrz_mol(maxmls+1,mxatyp),stat=fail(9))
      allocate (numtyp_mol(maxmls+1,mxatyp),stat=fail(10))
      allocate (numfrz_gstmol(maxmls,mxatyp),stat=fail(11))
      allocate (numtyp_gstmol(maxmls,mxatyp),stat=fail(12))
      allocate (delrc_mol(mxcombo),stat=fail(13))
      allocate (delrc_mol0(mxcombo),stat=fail(14))

      do i=1,nv
        if(fail(i).gt.0)then
            if(idnode.eq.0)write(nrite,'(10i5)')fail(i)
            call error(idnode, 1002)
        endif
      enddo
      vdwen = 0.d0
      elrc=0.d0
      origelrc=0.d0
      elrc_mol=0.d0
      origelrc_mol=0.d0
      numfrz_mol=0
      numtyp_mol=0
      numfrz_gstmol=0
      numtyp_gstmol=0
      end subroutine alloc_vdw_arrays

      end module vdw_module
