      module utility_pack
      use parse_module
      implicit none

      integer nsite
      integer mxlist,mxegrd
      logical, allocatable :: lprobeng(:)
      real(8) volum, total_pressure
      real(8), allocatable :: xxx(:),molxxx(:,:),origmolxxx(:,:)
      real(8), allocatable :: yyy(:),molyyy(:,:),origmolyyy(:,:)
      real(8), allocatable :: zzz(:),molzzz(:,:),origmolzzz(:,:)
      real(8), allocatable :: newx(:),newy(:),newz(:)
      real(8), allocatable :: guestx(:,:),guesty(:,:),guestz(:,:)
      real(8), allocatable :: atmcharge(:),atmchg(:,:)
      real(8), allocatable :: atmweight(:),atmwght(:,:)
      real(8), allocatable :: frambuff(:,:)
      real(8), allocatable :: xdf(:),ydf(:),zdf(:),rsqdf(:) 
      character(8), allocatable :: unqatm(:),atomname(:),atmname(:,:)
      character*70, allocatable :: numgbuff(:,:)
      character*1, allocatable :: molnam(:,:)
      integer, allocatable :: lentry(:),list(:,:)
      integer, allocatable :: gstlentry(:),gstlist(:,:),locguest(:)
      integer, allocatable :: locfram(:),moldf(:),moltype(:)
      integer, allocatable :: lfreezesite(:),lfzsite(:,:)
      integer, allocatable :: nummols(:),numatoms(:)
      integer, allocatable :: guest_insert(:)  
      real(8), dimension(9) :: cell
      real(8), dimension(9) :: rcell
      real(8), dimension(9) :: ucell
      integer, allocatable :: lexatm(:,:),ltpsit(:,:)
      integer, allocatable :: ltype(:)
      integer, allocatable :: nexatm(:),noxatm(:),lexsit(:,:,:)
      integer, allocatable :: nexsit(:,:)
      integer, allocatable :: numfrz(:),numtyp(:),dens(:)
      integer, allocatable :: ins(:),del(:),dis(:),jmp(:),flx(:)
      integer, allocatable :: swp(:),swi(:)
      integer, allocatable :: switch_mol_count(:), switch_mols(:,:)
      integer, allocatable :: switch_chosen_guest(:)
      integer, allocatable :: ind(:),ilist(:),jlist(:),angdist(:)
      integer, allocatable :: nprob(:),nprobsites(:),lprobsites(:,:)

      real(8), allocatable :: grid(:,:)
      real(8), allocatable :: dbuff(:),delE(:)
      real(8), allocatable :: statbuff(:),chainstats(:)
      real(8), allocatable :: energy(:), node_avg(:,:),nodeweight(:)
      real(8), allocatable :: origenergy(:),molmass(:)
      real(8), allocatable :: surfacemols(:),origsurfmols(:)
      real(8), allocatable :: node_std(:,:)
      real(8), allocatable :: avgwindow(:), varwindow(:), sumwindowav(:)
c     Fugacity stuff
      real(8), allocatable :: gstpress(:)
      real(8), allocatable :: gstfuga(:), gstmolfract(:)
      real(8), allocatable :: Acc_factor(:), P_crit(:), T_crit(:)
      real(8), allocatable :: K_fug(:,:) 

c     FIXED PARAMETERS

      integer, parameter :: maxatm=20000

      integer, parameter :: mxexcl=50

c     maximum number of guests in the framework for gcmc sim
      integer, parameter :: maxguest=20000

c     max number of sites (atoms) for guest molecule
      integer, parameter :: mxguestsite=50
c     standard pi values

      real(8), parameter :: pi=3.141592653589793d0
      real(8), parameter :: sqrpi=1.7724538509055159d0

c     angstroms to bohr converter for .cube files

      real(8), parameter :: angs2bohr=1.889725989d0


c     min angle for gcmc rotation
      real(8), parameter :: minangle=pi/18.d0
c      real(8), parameter :: maxangle=pi/3.d0
c     min translation for gcmc displacement
      real(8), parameter :: delrmin=0.1d0
      real(8), parameter :: delrmax=3.d0

c     conversion factor for coulombic terms in internal units
c     i.e. (unit(charge)**2/(4 pi eps0 unit(length))/unit(energy))

      real(8), parameter :: r4pie0=138935.4835d0

c     boltzmann constant in internal units

      real(8), parameter :: iboltz=8.31451115d-1

c     boltzmann constant in kcal/mol/K units 

      real(8), parameter :: kboltz=1.9872041d-3

c     boltzmann constant in kg m^2 / (s^2 K)
      real(8), parameter :: boltz=1.3806503d-23 

c     Avogadro's number 
      real(8), parameter :: avo=6.022140857d23

c     Gas constant in terms of kcal / (mol K)
      real(8), parameter :: Rgas=8.314472d0/4184.d0

c     planck's constant in internal units

      real(8), parameter :: hbar=6.350780668d0

c     planck's constant in (m^2 kg) / s

      real(8), parameter :: hplanck = 6.62607004d-34

c     conversion factor for pressure from internal units to katm

      real(8), parameter :: prsunt=0.163882576d0

c     main input channel

      integer, parameter :: ncontrol=5

c     main output channel

      integer, parameter :: nrite=6

c     force field input channel

      integer, parameter :: nfield=9

c     configuration file input channel

      integer,parameter :: nconfig=10

c     history file input channel
 
      integer,parameter :: nhist=11

c     statistics file input channel

      integer,parameter :: nstats=12

      integer,parameter :: nrev=555


      save xxx,yyy,zzz
      save molxxx,molyyy,molzzz,frambuff
      save origmolxxx,origmolyyy,origmolzzz
      save atmcharge,atmchg,guestx,guesty,guestz
      save atmweight, atmwght
      save atomname,atmname
      save lentry,list
      save lfreezesite,lfzsite
      save nummols,numatoms,molnam
      save ucell,cell,rcell,volum,nsite,chainstats,dbuff
      save gstlentry,gstlist,locguest,locfram,statbuff
      save avgwindow, varwindow, sumwindowav,moldf,moltype
      save gstpress,angdist,node_avg,node_std,nodeweight
      save xdf,ydf,zdf,rsqdf
      save ins,del,dis,jmp,flx,swp,swi
      save switch_mol_count,switch_mols,switch_chosen_guest
      save lexatm,nexatm,noxatm,ilist,jlist,unqatm,ltpsit
      save mxlist,ltype,mxegrd,lexsit,nexsit
      save numtyp,numfrz,dens,ind,newx,newy,newz
      save nprob,nprobsites,lprobsites,grid
      save delE, energy, total_pressure, gstfuga, gstmolfract
      save origenergy,guest_insert,molmass
      save Acc_factor, P_crit, T_crit, K_fug
      save origsurfmols,surfacemols,lprobeng
      contains

      subroutine initscan
     &(idnode,imcon,volm,keyfce,rcut,eps,alpha,kmax1,kmax2,kmax3,lprob,
     &delr,rvdw,ntpguest,ntprob,ntpsite,ntpvdw,maxmls,mxatm,mxatyp,
     &griddim, gridfactor)
c**********************************************************************
c
c     scans input files for relevant maximums to allocate to arrays
c     communicate to all nodes??
c**********************************************************************
      
      implicit none
      integer, parameter :: mmk=1000

      logical loop,loop2,loop3,loop4,loop5,safe,lewald,lcut
      logical lrvdw,check,ldelr,kill,lprob
      character*8 name, chr(mmk)
      integer imcon,keyfce,idnode,idum,ntpvdw,maxmls,mxatm
      integer n,nummls,numsit,kmax1,kmax2,kmax3
      integer mxatyp,nrept,ifrz,nneu,ksite,isite
      integer j,ntpguest,ntprob,ntpsite,temp,gsite,iprob,qprob
      real(8) alpha,delr,rvdw,ppp,width
      real(8) fac,tol,tol1,rcut,eps,volm
      real(8), dimension(10) :: celprp
      real(8) griddim
      integer, dimension(3) :: gridfactor
      mxatm=0
      mxatyp=0
      keyfce=0
      ntprob=0
      ntpsite=0
      ntpguest=0
      data loop/.true./,loop2/.false./,loop4/.false./,lewald/.false./
      data lrvdw/.false./,ldelr/.false./,kill/.false./,loop5/.false./

      if(idnode.eq.0)open(nfield,file='FIELD',status='old')

      call getrec(safe,idnode,nfield)
      if(.not.safe)call abort_field_read(1,idnode,nfield)
      do while(loop)
        call getrec(safe,idnode,nfield)
        if(.not.safe)call abort_field_read(1,idnode,nfield)
        call lowcase(record, lenrec)
        call strip(record, lenrec)
        if(record(1).eq.'#'.or.record(1).eq.' ')then
        elseif (findstring('molecu',record,idum))then
           maxmls=intstr(record,lenrec,idum)

           do n=1,maxmls
             loop2=.true.
            
             do while(loop2)
               call getrec(safe,idnode,nfield)
               if(.not.safe)call abort_field_read(1,idnode,nfield)
               call lowcase(record,lenrec)
               call strip(record,lenrec)
               ksite=0

               if(record(1).eq.'#'.or.record(1).eq.' ')then
               elseif (findstring('nummol',record,idum))then
                 nummls=intstr(record,lenrec,idum)
               elseif (findstring('atoms',record,idum))then
                 numsit=intstr(record,lenrec,idum)
                 mxatm=mxatm+numsit*nummls
                 ksite=0

                 do isite=1,numsit
                   if (ksite.lt.numsit)then
                      call getrec(safe,idnode,nfield)
                      if(.not.safe)call abort_field_read
     &                    (1,idnode,nfield)
                      call getword(name,record,8,lenrec)
                      ppp=dblstr(record,lenrec,idum)
                      ppp=dblstr(record,lenrec,idum)
                      nrept=intstr(record,lenrec,idum)
                      ifrz=intstr(record,lenrec,idum)
                      nneu=intstr(record,lenrec,idum)
                      if(nrept.eq.0)nrept=1
                      ksite=ksite+nrept
                      if(mxatyp.eq.0)then
                         mxatyp=1
                         chr(1)=name
                      else
                         check=.true.
                         do j=1,mxatyp
                           if(name.eq.chr(j))check=.false.
                         enddo
                         if(check)then
                           mxatyp=mxatyp+1
                           if(mxatyp.lt.mmk)chr(mxatyp)=name
                         endif
                      endif
                   endif
                 enddo
               elseif (findstring('finish',record,idum))then
                 loop2=.false.
               endif
             enddo
           enddo
        elseif(findstring('vdw',record,idum))then
           ntpvdw=intstr(record,lenrec,idum)
        elseif(findstring('close',record,idum))then
           loop=.false.
        endif
      enddo
c     grab the cell vectors - need these to allocate
c     some of the arrays.
      if(idnode.eq.0)open(nconfig,file='CONFIG',status='old')
      call getrec(safe,idnode,nconfig)
      call getrec(safe,idnode,nconfig)
      imcon=intstr(record,lenrec,idum)
      imcon=intstr(record,lenrec,idum)
      call getrec(safe,idnode,nconfig)
      if(.not.safe)call abort_config_read(1,idnode,nconfig)
      cell(1)=dblstr(record,lenrec,idum)
      cell(2)=dblstr(record,lenrec,idum)
      cell(3)=dblstr(record,lenrec,idum)
      call getrec(safe,idnode,nconfig)
      if(.not.safe)call abort_config_read(1,idnode,nconfig)
      cell(4)=dblstr(record,lenrec,idum)
      cell(5)=dblstr(record,lenrec,idum)
      cell(6)=dblstr(record,lenrec,idum)
      call getrec(safe,idnode,nconfig)
      if(.not.safe)call abort_config_read(1,idnode,nconfig)
      cell(7)=dblstr(record,lenrec,idum)
      cell(8)=dblstr(record,lenrec,idum)
      cell(9)=dblstr(record,lenrec,idum)
      call dcell(cell,celprp)
      if(idnode.eq.0)open(ncontrol,file='CONTROL',status='old')
      loop3=.true.
      ntpguest=0
      do while(loop3)
        call getrec(safe,idnode,ncontrol)
        if(.not.safe)call abort_control_read(1,idnode,ncontrol)
        call lowcase(record,lenrec)
        call strip(record,lenrec)
        if (record(1).eq.'#'.or.record(1).eq.' ')then

        elseif(findstring('cut',record,idum))then
          rcut=dblstr(record,lenrec,idum)
          lcut=.true.
        elseif (findstring('&guest',record,idum))then
           ntpguest=ntpguest+1
           gsite=0 
           loop4=.true.
           
           do while(loop4)
         
             call getrec(safe,idnode,ncontrol)
             call lowcase(record,lenrec)
             call strip(record,lenrec)
             if(record(1).eq.'#'.or.record(1).eq.' ')then
c            record is commented out
             elseif(findstring('probability',record,idum))then
               iprob=0
               lprob=.true.
               qprob=intstr(record,lenrec,idum)
               ntprob=ntprob+qprob
               loop5=.true.
               do while(loop5) 
                 call getrec(safe,idnode,ncontrol)
                 if(.not.safe)call abort_control_read(1,idnode,ncontrol)
                 call strip(record,lenrec)
                 if(record(1).eq.'#'.or.record(1).eq.' ')then
                 elseif(record(1).eq.'e')then
                   ntprob=ntprob+1
                   iprob=iprob+1
                 else
                   iprob=iprob+1
                   temp=intstr(record,lenrec,idum)
                   gsite=gsite+temp
                 endif
                 if(iprob.eq.qprob)loop5=.false.
               enddo
               
             elseif(findstring('&end',record,idum))then
               loop4=.false.
               ntpsite=max(gsite,ntpsite)
             endif    
           enddo 

        elseif(findstring('delr',record,idum))then
          delr=dblstr(record,lenrec,idum)
          ldelr=.true.
        elseif(findstring('finish',record,idum))then
          loop3=.false.
        elseif(findstring('ewald',record,idum))then
          lewald=.true.
          keyfce=2
          if (findstring('precision',record,idum))then
             eps=dblstr(record,lenrec,idum)
             if (.not.lcut)then
                call error(idnode,-433)
                kill=.true. 
             else
               if(lewald)then
      
                  eps=min(abs(eps),0.5d0)
                  tol=sqrt(abs(log(eps*rcut)))
                  alpha=sqrt(abs(log(eps*rcut*tol)))/rcut
c                  print *, "sigma is, ", 1/sqrt(2.d0)/alpha
                  tol1=sqrt(-log(eps*rcut*(2.d0*tol*alpha)**2))
                  fac=1.d0
                  if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7)then
                     fac=2.d0**(1.d0/3.d0)
                  endif
                  kmax1=nint(0.25d0+fac*celprp(1)*alpha*tol1/pi)
                  kmax2=nint(0.25d0+fac*celprp(2)*alpha*tol1/pi)
                  kmax3=nint(0.25d0+fac*celprp(3)*alpha*tol1/pi)
               endif
             endif
          else
              alpha=dblstr(record,lenrec,idum)
              kmax1=intstr(record,lenrec,idum)
              kmax2=intstr(record,lenrec,idum)
              kmax3=intstr(record,lenrec,idum)
          endif
        elseif(findstring('rvdw',record,idum))then
          rvdw=dblstr(record,lenrec,idum)
          lrvdw=.true.
c Need to find the grid parameters in first scan, before allocating
        elseif (findstring('grid',record,idum))then
          if (findstring('spa',record,idum))then
            griddim = dblstr(record,lenrec,idum)
          elseif (findstring('fac',record,idum))then
            gridfactor(1) = intstr(record,lenrec,idum)
            gridfactor(2) = intstr(record,lenrec,idum)
            gridfactor(3) = intstr(record,lenrec,idum)
          endif

        endif
      enddo
      if(.not.ldelr)then
        call error(idnode,-433)
        kill=.true.
      endif

      width=min(celprp(7),celprp(8),celprp(9))/2.d0
      if(imcon.eq.4)width=sqrt(3.d0)*cell(1)/4.d0
      if(imcon.eq.5)width=cell(1)/2.d0
      if(imcon.eq.6)width=min(celprp(7),celprp(8))/2.d0

c     halt program if cutoff exceeds cell width

      if(rcut.gt.width)call error(idnode,95)
      
      call dcell(cell,celprp)
      if(imcon.eq.0) then
        volm=0.d0
      elseif(imcon.eq.4)then
        volm=0.5d0*celprp(10)
      elseif(imcon.eq.5)then
        volm=0.5d0*celprp(10)
      elseif(imcon.eq.7)then
        volm=0.5d0*celprp(10)

      else
        volm=celprp(10)
      endif

      if(.not.lrvdw)then
        if(rcut.gt.0d0)then
           rvdw=rcut
        endif
      endif
      if(.not.lewald)then
        call error(idnode,2311)
      endif
      if(idnode.eq.0)then
        close(ncontrol)  
        close(nconfig) 
        close(nfield)
      endif
      return
      end subroutine initscan 

      subroutine alloc_prob_arrays
     &(idnode,ntpguest,ntpsite,ntprob,gridsize)
c*********************************************************************
c
c      subroutine to allocate probability arrays
c
c*********************************************************************
      implicit none
      integer, parameter :: np=5
      integer ntpguest,ntprob,ntpsite,gridsize
      integer i,j, idnode
      integer, dimension(np) :: fail 
      do i=1,np
        fail(i) = 0
      enddo

      allocate(nprob(ntpguest),stat=fail(1))
      allocate(nprobsites(ntprob),stat=fail(2))
      allocate(lprobsites(ntprob,ntpsite),stat=fail(3))
      allocate(grid(ntprob,gridsize),stat=fail(4))
      allocate(lprobeng(ntpguest),stat=fail(5))
      do i=1,np
        if(fail(i).gt.0)then
            if(idnode.eq.0)write(nrite,'(10i5)')fail
            call error(idnode, 1001)
        endif
      enddo
      do i=1,ntpguest
        lprobeng(i)=.false.
      enddo
      do i=1,ntprob
        nprobsites(i)=0
        do j=1,gridsize
          grid(i,j)=0.d0
        enddo
      enddo

      return
      end subroutine alloc_prob_arrays

      subroutine storeenprob
     &(iguest,imol,rcell,ntpguest,ngrida,ngridb,ngridc,energy)
c*********************************************************************
c
c      subroutine to store energy in a grid
c
c*********************************************************************
      implicit none
      integer mol,imol,np,iguest,itprob
      integer jj,nmols,natms,itatm
      integer ngrida,ngridb,ngridc,ntpguest
      integer j,jguest
      real(8) comx,comy,comz,energy
      real(8), dimension(9) :: rcell 

      itprob=0

c     iterate over all molecules (outer loop) and atoms
c     (inner loop), store each guest in a temp array,
c     get fractional coordinates.
      np=nprob(iguest)
      mol=locguest(iguest)
      nmols=nummols(mol)
      natms=numatoms(mol)
c     commented out because the newx,y,and z arrays should be 
c     populated by the correct coordinates from the gcmc move.
      if(imol.gt.0)then
        jj = natms*(imol-1)
        do itatm=1,natms
          jj=jj+1
          newx(itatm)=molxxx(mol,jj)
          newy(itatm)=molyyy(mol,jj)
          newz(itatm)=molzzz(mol,jj)
        enddo
      endif

c     center of mass - calculated for the energy position.

      call com(natms,mol,newx,newy,newz,comx,comy,comz)
c     convert centre of mass to fractional coordinates

      call frac2(comx,comy,comz,rcell)
c      aidx=int(ngrida*(comx))
c      bidx=int(ngridb*(comy))
c      cidx=int(ngridc*(comz))

c      if(aidx.eq.0)aidx=aidx+1
c      if(bidx.eq.0)bidx=bidx+1
c      if(cidx.eq.0)cidx=cidx+1

      itprob=0
      do jguest=1,ntpguest+1
        np=nprob(jguest)
        itprob=itprob+np
        if (jguest.eq.iguest)exit
      enddo
c     Energy prob is the last two grids.
c     TODO(pboyd): check indices here...
      j=itprob-1
      call equitable_bin(comx,comy,comz,ngrida,ngridb,ngridc,
     &j,energy)
      call equitable_bin(comx,comy,comz,ngrida,ngridb,ngridc,
     &itprob,1.d0)
c      gidx=(aidx)*ngridb*ngridc+(bidx)*ngridc+(cidx)
c      grid(j,gidx)=grid(j,gidx)+energy
c      grid(j+1,gidx)=grid(j+1,gidx)+1.d0

      return
      end subroutine storeenprob

      subroutine storeprob
     &(ntpguest,rcell,ngrida,ngridb,ngridc)
c*********************************************************************
c
c      subroutine to store guest positions in a grid
c
c*********************************************************************
      implicit none
      integer mol,np,iprob,iguest,i,itprob
      integer jj,itmols,nmols,natms,itatm,jatm
      integer ngrida,ngridb,ngridc,ntpguest
      real(8) comx,comy,comz
      real(8), dimension(9) :: rcell 

      itprob=0

c     iterate over all molecules (outer loop) and atoms
c     (inner loop), store each guest in a temp array,
c     get fractional coordinates.
      do iguest=1,ntpguest
        np=nprob(iguest)
        if(lprobeng(iguest))np=np-2
        mol=locguest(iguest)
        nmols=nummols(mol)
        natms=numatoms(mol)
        jj=0
        do itmols=1,nmols
          do itatm=1,natms
            jj=jj+1
            newx(itatm)=molxxx(mol,jj)
            newy(itatm)=molyyy(mol,jj)
            newz(itatm)=molzzz(mol,jj)
          enddo

c         center of mass - calculated regardless of COM request for
c         now..

          call com(natms,mol,newx,newy,newz,comx,comy,comz)

c         convert centre of mass to fractional coordinates

          call frac2(comx,comy,comz,rcell)

c         convert guest atoms to fractional coordinates
          call frac(newx,newy,newz,natms,rcell)

c       iterate over number of probability plots
        
          do iprob=1,np

            do i=1,nprobsites(itprob+iprob)
              if(lprobsites(itprob+iprob,i).eq.0)then  
                call equitable_bin(comx,comy,comz,ngrida,ngridb,ngridc,
     &itprob+iprob,1.d0)
              else
                jatm=lprobsites(itprob+iprob,i)
                call equitable_bin(newx(jatm),newy(jatm),newz(jatm),
     &ngrida,ngridb,ngridc,itprob+iprob, 1.d0)
              endif
            enddo
          enddo
        enddo
        itprob=itprob+np
        if(lprobeng(iguest))itprob=itprob+2
      enddo

      return
      end subroutine storeprob
      
      subroutine writeenprob
     &(iguest,itprob,cell,ntpfram,gridsize,
     &ngrida,ngridb,ngridc,scell_factor)
c*********************************************************************
c
c     write gaussian .cube files for visualization
c
c*********************************************************************
      implicit none
      character*25 filename
      real(8) vol, val
      integer itprob,iguest,i,j
      integer igrid,ngrida,ngridb,ngridc,ntpfram
      integer nfram,framol,nmol,natms,iatm,m,n,o
      integer gridsize,ip,scell_factor
      real(8), dimension(9) :: cell

      vol=(cell(2)*cell(6)-cell(3)*cell(5))/(ngrida*ngridb)*
     &cell(7)/ngridc+
     &(cell(3)*cell(4)-cell(1)*cell(6))/(ngrida*ngridb)*
     &cell(8)/ngridc+
     &(cell(1)*cell(5)-cell(2)*cell(4))/(ngrida*ngridb)*
     &cell(9)/ngridc
c     get number of framework atoms
      nfram=0
      do i=1,ntpfram
        framol=locfram(i)
        nfram=nfram+nummols(framol)/dble(scell_factor)*numatoms(framol)
      enddo
 
      write(filename,"('Ehist_guest',i2.2,'.cube')")
     &iguest
      open(i,file=filename)

c     write header stuff for cube
      write(i,"('Energy histogram cube file',/,'outer loop a,
     & middle loop b, inner loop c')")
      write(i,'(i6, 3f12.6)')nfram,0.d0,0.d0,0.d0
      write(i,'(i6,3f12.6)')ngrida,cell(1)/ngrida,cell(2)/ngrida,
     &cell(3)/ngrida
      write(i,'(i6,3f12.6)')ngridb,cell(4)/ngridb,cell(5)/ngridb,
     &cell(6)/ngridb
      write(i,'(i6,3f12.6)')ngridc,cell(7)/ngridc,cell(8)/ngridc,
     &cell(9)/ngridc
c     write cartesians of framework
      do m=1,ntpfram
        framol=locfram(m)
        nmol=int(nummols(framol)/dble(scell_factor))
        natms=numatoms(framol)
        iatm=0
        do n=1,nmol
          do o=1,natms
            iatm=iatm+1
            write(i,'(i6,4f12.6)')atmnumber(atmwght(framol,o)),0.d0,
     &               molxxx(framol,iatm)*angs2bohr,
     &               molyyy(framol,iatm)*angs2bohr,
     &               molzzz(framol,iatm)*angs2bohr

          enddo
        enddo
      enddo
      igrid=0

c     inner loop is c vec, middle loop is b vec, outter loop is a vec.
c     this has been intrinsically stored in the array grid.
c     the cube file writes six entries per line *except* when it reaches
c     the end of the number of grid points in the c direction (fastest
c     loop).  In this case it will end the line before reaching the
c     sixth entry.  A nuisance.
      ip=0
      do j=1,gridsize
        ip=ip+1 
        if(grid(itprob+1,j).gt.0.d0)then
            val=grid(itprob,j)/grid(itprob+1,j)
        else
            val = 0.d0 
        endif
        write(i,'(e15.6)',advance='no')
     &(dble(val))

c       conditions for writing to a new line 
        if((mod(j,ngridc).eq.0).or.(ip.eq.6))then
           write(i,'(x)')
           ip=0
        endif
      enddo
      close(i)

      return
      end subroutine writeenprob

      subroutine writeprob
     &(iguest,itprob,iprob,cell,ntpfram,gridsize,
     &ngrida,ngridb,ngridc,steps,scell_factor)
c*********************************************************************
c
c     write gaussian .cube files for visualization
c
c*********************************************************************
      implicit none
      character*25 filename
      real(8) vol
      integer itprob,iguest,i,j,scell_factor
      integer igrid,ngrida,ngridb,ngridc,ntpfram
      integer nfram,framol,nmol,natms,iatm,m,n,o
      integer gridsize,steps,iprob,ip
      real(8), dimension(9) :: cell

      vol=(cell(2)*cell(6)-cell(3)*cell(5))/(ngrida*ngridb)*
     &cell(7)/ngridc+
     &(cell(3)*cell(4)-cell(1)*cell(6))/(ngrida*ngridb)*
     &cell(8)/ngridc+
     &(cell(1)*cell(5)-cell(2)*cell(4))/(ngrida*ngridb)*
     &cell(9)/ngridc
c     get number of framework atoms
      nfram=0
      do i=1,ntpfram
        framol=locfram(i)
        nfram=nfram+nummols(framol)/dble(scell_factor)*numatoms(framol)
      enddo
 
      write(filename,"('prob_guest',i2.2,'_prob_',i2.2,'.cube')")
     &iguest,iprob
      open(i,file=filename)

c     write header stuff for cube
      write(i,"('probability cube file',/,'outer loop a, middle loop
     & b, inner loop c')")
      write(i,'(i6, 3f12.6)')nfram,0.d0,0.d0,0.d0
      write(i,'(i6,3f12.6)')ngrida,cell(1)/ngrida,cell(2)/ngrida,
     &cell(3)/ngrida
      write(i,'(i6,3f12.6)')ngridb,cell(4)/ngridb,cell(5)/ngridb,
     &cell(6)/ngridb
      write(i,'(i6,3f12.6)')ngridc,cell(7)/ngridc,cell(8)/ngridc,
     &cell(9)/ngridc
c     write cartesians of framework
      do m=1,ntpfram
        framol=locfram(m)
        nmol=nummols(framol)/(dble(scell_factor))
        natms=numatoms(framol)
        iatm=0
        do n=1,nmol
          do o=1,natms
            iatm=iatm+1
            write(i,'(i6,4f12.6)')atmnumber(atmwght(framol,o)),0.d0,
     &               molxxx(framol,iatm)*angs2bohr,
     &               molyyy(framol,iatm)*angs2bohr,
     &               molzzz(framol,iatm)*angs2bohr

          enddo
        enddo
      enddo
      igrid=0

c     inner loop is c vec, middle loop is b vec, outter loop is a vec.
c     this has been intrinsically stored in the array grid.
c     the cube file writes six entries per line *except* when it reaches
c     the end of the number of grid points in the c direction (fastest
c     loop).  In this case it will end the line before reaching the
c     sixth entry.  A nuisance.
      ip=0
      do j=1,gridsize
        ip=ip+1 
        write(i,'(e15.6)',advance='no')
     &(dble(grid(itprob,j))/vol/dble(steps))

c       conditions for writing to a new line 
        if((mod(j,ngridc).eq.0).or.(ip.eq.6))then
           write(i,'(x)')
           ip=0
        endif
      enddo
      close(i)
  

      return
      end subroutine writeprob

      subroutine frac(x,y,z,natms,rcell)
c*********************************************************************
c
c     convert cartesian coordinates to fractional coordinates 
c     haven't checked if it works with other cell types   
c     than parallelpiped (imcon.eq.3) and rectangular (imcon.eq.1 or 2)
c     PB
c 
c*********************************************************************
      implicit none
      integer i,natms
      real(8) ssx,ssy,ssz
      real(8), dimension(natms) :: x,y,z
      real(8), dimension(9) :: rcell


      do i=1,natms

        ssx=x(i)*rcell(1)+y(i)*rcell(4)+z(i)*rcell(7)
        ssy=x(i)*rcell(2)+y(i)*rcell(5)+z(i)*rcell(8)
        ssz=x(i)*rcell(3)+y(i)*rcell(6)+z(i)*rcell(9)
        x(i)=modulo(ssx, 1.0)
        y(i)=modulo(ssy, 1.0)
        z(i)=modulo(ssz, 1.0)
      enddo

      return
      end subroutine frac

      subroutine frac2(x,y,z,rcell)
c*********************************************************************
c
c     convert cartesian coordinates to fractional coordinates
c     for one atom only (this is a cheap fix for the COM calcluation
c     in the subroutine storeprob)
c     haven't checked if it works with other cell types   
c     than parallelpiped (imcon.eq.3) and rectangular (imcon.eq.1 or 2)
c     PB
c 
c*********************************************************************
      implicit none
      real(8) ssx,ssy,ssz
      real(8) x,y,z
      real(8), dimension(9) :: rcell

      ssx=x*rcell(1)+y*rcell(4)+z*rcell(7)
      ssy=x*rcell(2)+y*rcell(5)+z*rcell(8)
      ssz=x*rcell(3)+y*rcell(6)+z*rcell(9)
      x=modulo(ssx, 1.0)
      y=modulo(ssy, 1.0)
      z=modulo(ssz, 1.0)
 
      return
      end subroutine frac2

      subroutine alloc_config_arrays
     &(idnode,mxnode,maxmls,mxatm,mxatyp,volm,
     &ntpguest,rcut,rvdw,delr)
c*********************************************************************
c
c     allocation of some arrays
c     
c     For reference some of the statistical arrays:
c      node_avg(1) = <N>
c      node_avg(2) = <E>
c      node_avg(3) = <EN>
c      node_avg(4) = <N2>
c      node_avg(5) = <E2>
c      node_avg(6) = <NF>  -- molecules close to the framework
c      node_avg(7) = Q_st
c      node_avg(8) = C_v
c      node_avg(9) = <exp(-E/kb/T)> 
c      
c      node_std(1) = N
c      node_std(2) = E
c      node_std(3) = EN
c      node_std(4) = N2
c      node_std(5) = E2
c      node_std(6) = NF
c      node_std(7) = Q_st
c      node_std(8) = C_v
c      node_std(9) = exp(-E/kb/T)
c      
c      chainstats(1) = prodcount
c      chainstats(2) = <N>
c      chainstats(3) = <E>
c      chainstats(4) = <EN>
c      chainstats(5) = <N2>
c      chainstats(6) = <E2>
c      chainstats(7) = <NF>
c      chainstats(8) = <exp(-E/kb/T)> 
c      chainstats(9) = stdN
c      chainstats(10) = stdE
c      chainstats(11) = stdEN
c      chainstats(12) = stdN2
c      chainstats(13) = stdE2
c      chainstats(14) = stdNF
c      chainstats(15) = stdQ_st
c      chainstats(16) = stdC_v
c      chainstats(17) = std[exp(-E/kb/T)]
c      
c      avgwindow(1) = <N>
c      avgwindow(2) = <E>
c      avgwindow(3) = <EN>
c      avgwindow(4) = <N2>
c      avgwindow(5) = <E2>
c      avgwindow(6) = <NF>
c      avgwindow(7) = <Qst>
c      avgwindow(8) = <Cv>
c      avgwindow(9) = <exp(-E/kb/T)>
c
c*********************************************************************
      implicit none
      integer, parameter :: na = 85 
      integer maxmls,mxatm,maxalloc,ntpguest,i,j
      integer mxatyp,idnode,mxnode,mxcmls
      real(8) density,ratio,cut,volm,rcut,rvdw,delr
      integer, dimension(na) :: fail

c     initialize fail check array
      do i=1,na
        fail(i) = 0
      enddo
      mxcmls=maxmls* (maxmls-1) / 2 + maxmls
      mxegrd=max(1000,int(rvdw/0.01d0+0.5d0)+4)
      maxalloc=mxatm+maxguest
c     the following line was added to test the size of maxalloc
c     which seemed to be giving insufficient virutal memory errors
c      if(idnode.eq.0)write(nrite,"('maxalloc: ', i9)")maxalloc
      density=dble(maxalloc)/volm
      cut=rcut+delr
      ratio=1.5d0*density*(4.d0*pi/3.d0)*cut**3
      mxlist=min(nint(ratio),(maxalloc+1)/2)
      allocate(locguest(maxmls),stat=fail(1))
      allocate(guestx(ntpguest,mxguestsite),stat=fail(2))
      allocate(guesty(ntpguest,mxguestsite),stat=fail(3))
      allocate(guestz(ntpguest,mxguestsite),stat=fail(4))
      allocate(locfram(maxmls),stat=fail(5))
      allocate(molmass(maxmls),stat=fail(6))
      allocate(list(maxalloc,mxlist),stat=fail(7))
      allocate(lentry(maxalloc),stat=fail(8))
      allocate(noxatm(maxalloc),stat=fail(9))
      allocate(nexatm(maxalloc),stat=fail(10))
      allocate(lexatm(maxalloc,mxexcl),stat=fail(11))
      allocate(nexsit(mxguestsite,mxexcl),stat=fail(12))
      allocate(lexsit(mxguestsite,mxguestsite,mxexcl),stat=fail(13))
      allocate(xxx(maxalloc),stat=fail(14))
      allocate(yyy(maxalloc),stat=fail(15))
      allocate(zzz(maxalloc),stat=fail(16))
      allocate(molnam(40,maxmls),stat=fail(17))
      allocate(nummols(maxmls),stat=fail(18))
      allocate(numatoms(maxmls),stat=fail(19))
      allocate(atmname(maxmls,maxalloc),stat=fail(20))
      allocate(atmchg(maxmls,maxalloc),stat=fail(21))
      allocate(atmwght(maxmls,maxalloc),stat=fail(22))
      allocate(lfzsite(maxmls,maxalloc),stat=fail(23))
      allocate(molxxx(maxmls,maxalloc),stat=fail(24))
      allocate(molyyy(maxmls,maxalloc),stat=fail(25))
      allocate(molzzz(maxmls,maxalloc),stat=fail(26))
      allocate(frambuff(maxmls,maxalloc),stat=fail(27))
      allocate(lfreezesite(maxalloc),stat=fail(28))
      allocate(atmcharge(maxalloc),stat=fail(29))
      allocate(atomname(maxalloc),stat=fail(30))
      allocate(atmweight(maxalloc),stat=fail(31))
      allocate(gstlentry(mxguestsite),stat=fail(32))
      allocate(gstpress(ntpguest),stat=fail(33))
      allocate(statbuff(1+ntpguest*16),stat=fail(34))
      allocate(chainstats(1+ntpguest*16),stat=fail(35))
      allocate(gstlist(mxguestsite,maxalloc),stat=fail(36))
      allocate(ilist(maxalloc),stat=fail(37))
      allocate(unqatm(maxalloc),stat=fail(38))
      allocate(jlist(maxalloc),stat=fail(39))
      allocate(ltpsit(maxmls,maxalloc),stat=fail(40))
      allocate(ltype(maxalloc),stat=fail(41))
      allocate(xdf(maxalloc),stat=fail(42))
      allocate(ydf(maxalloc),stat=fail(43))
      allocate(zdf(maxalloc),stat=fail(44))
      allocate(rsqdf(maxalloc),stat=fail(45))
      allocate(numtyp(mxatyp),stat=fail(46))
      allocate(numfrz(mxatyp),stat=fail(47))
      allocate(dens(mxatyp),stat=fail(48))
      allocate(newx(mxguestsite),stat=fail(49))
      allocate(newy(mxguestsite),stat=fail(50))
      allocate(newz(mxguestsite),stat=fail(51))
      allocate(ind(mxguestsite),stat=fail(52))
      allocate(energy(maxmls+1),stat=fail(53))
      allocate(origenergy(maxmls+1),stat=fail(54))
      allocate(delE(mxcmls),stat=fail(55))
      allocate(avgwindow(ntpguest*9),stat=fail(56))
      allocate(sumwindowav(ntpguest*9),stat=fail(57))
      allocate(varwindow(ntpguest*9),stat=fail(58))
      allocate(nodeweight(mxnode),stat=fail(59))
      allocate(node_avg(mxnode,ntpguest*9),stat=fail(60))
      allocate(node_std(mxnode,ntpguest*9),stat=fail(61))
      allocate(ins(ntpguest),stat=fail(62))
      allocate(del(ntpguest),stat=fail(63))
      allocate(dis(ntpguest),stat=fail(64))
      allocate(jmp(ntpguest),stat=fail(65))
      allocate(flx(ntpguest),stat=fail(66))
      allocate(swp(ntpguest),stat=fail(67))
      allocate(swi(ntpguest),stat=fail(68))
      allocate(switch_mol_count(ntpguest),stat=fail(69))
      allocate(switch_mols(ntpguest,maxguest),stat=fail(70))
      allocate(switch_chosen_guest(ntpguest),stat=fail(71))
      allocate(origmolxxx(maxmls,maxalloc),stat=fail(72))
      allocate(origmolyyy(maxmls,maxalloc),stat=fail(73))
      allocate(origmolzzz(maxmls,maxalloc),stat=fail(74))
      allocate(gstfuga(ntpguest),stat=fail(75))
      allocate(gstmolfract(ntpguest),stat=fail(76))
      allocate(Acc_factor(ntpguest),stat=fail(77))
      allocate(T_crit(ntpguest),stat=fail(78))
      allocate(P_crit(ntpguest),stat=fail(79))
      allocate(K_fug(ntpguest, ntpguest),stat=fail(80))
      allocate(surfacemols(maxmls), stat=fail(81))
      allocate(origsurfmols(maxmls), stat=fail(82))
      allocate(moltype(maxalloc), stat=fail(83))
      allocate(moldf(maxalloc), stat=fail(84))
      allocate(guest_insert(ntpguest), stat=fail(85))
      do i=1,na
        if(fail(i).gt.0)then
            if(idnode.eq.0)write(nrite,'(10i5)')fail(i)
            call error(idnode, 1000)
        endif
      enddo

c     initialize statistic arrays
c     probably put all this initialization stuff in a separate
c     module
      do i=1,mxnode
        nodeweight(i)=0.d0
        do j=1,ntpguest*9
          node_avg(i,j) = 0.d0
          node_std(i,j) = 0.d0
        enddo
      enddo
      do i = 1, ntpguest
        surfacemols(i) = 0.d0
        origsurfmols(i) = 0.d0
        do j = 1, ntpguest
            K_fug(i,j) = 0.d0
        enddo
      enddo
      do i=1,ntpguest*9
        avgwindow(i)=0.d0
        varwindow(i)=0.d0
        sumwindowav(i)=0.d0
      enddo
c     initialize the stat arrays
      chainstats(:)=0.d0
      do i=1,ntpguest
        guest_insert(i)=0
        ins(i) = 0.d0
        del(i) = 0.d0
        dis(i) = 0.d0
        jmp(i) = 0.d0
        flx(i) = 0.d0
        swp(i) = 0.d0
        swi(i) = 0.d0
      enddo
      do i=1,maxmls+1
        energy(i) = 0.d0
        origenergy(i) = 0.d0
        delE(i) = 0.d0
      enddo
      return
      end subroutine alloc_config_arrays


      subroutine fugacity(idnode, lfuga, temp, ntpguest)
c**********************************************************************
c
c     subroutine to compute fugacities from the supplied pressures 
c     of the gases 
c     
c*********************************************************************
      implicit none
      integer i, idnode, ntpguest
      logical lfuga,lfract,lpress,do_press
      real(8) temp,press_sum, fract_sum, no_fract_count
c     check if a mole fraction and total pressure is used, or just
c     guest partial pressures
      lfract=.false.
      lpress=.true.
      do_press=.true.
      fract_sum = 0.d0
      press_sum = 0.d0
      no_fract_count = 0.d0
      call fluid_properties(idnode, ntpguest, lfuga) 
      do i=1,ntpguest
        if(gstpress(i).lt.0.)then
            lpress=.false.
        elseif(gstpress(i).ge.0.)then
            press_sum = press_sum + gstpress(i)
        endif
        if(gstmolfract(i).gt.0.)then
            lfract=.true.
            fract_sum = fract_sum + gstmolfract(i)
        elseif(gstmolfract(i).le.0)then
            no_fract_count = no_fract_count + 1.d0
        endif
      enddo
      if ((lfract).and.(.not.lpress))then
        if(total_pressure.lt.0.)call error(idnode,2314)
        do_press=.false.
        if(ntpguest.eq.1)then
            gstmolfract(1) = 1.
            fract_sum = 1.
        endif
        if(fract_sum.lt.1.d0)then
            do i=1, ntpguest
              if(gstmolfract(i).lt.0)then
                  gstmolfract(i) = (1.d0 - fract_sum)/no_fract_count
              endif
            enddo
        endif
        fract_sum = 0.d0
        do i=1, ntpguest
            gstpress(i) = total_pressure*gstmolfract(i)
            fract_sum = fract_sum + gstmolfract(i)
        enddo        
        if((fract_sum.lt.0.999).or.(fract_sum.gt.1.001))then
            call error(idnode,2315)
        endif
      endif
      if((.not.lpress).and.(.not.lfract).and.(total_pressure.gt.0.))then
        do i=1, ntpguest
          gstpress(i)=total_pressure
        enddo
        total_pressure = ntpguest * total_pressure
        press_sum = total_pressure
      endif
      if((total_pressure.gt.0.).and.(.not.lfract))then
        if(press_sum.ne.total_pressure)then
          total_pressure = press_sum
          if(idnode.eq.0)write(nrite,"(/,a,f10.3,a,f10.3,/,a)")
     &'Warning - total pressure and sum of guess pressures do not 
     &match! total pressure: ',total_pressure, 'and sum of guest
     &pressures: ',press_sum, 'Using partial pressures supplied 
     &by the guests...'
        endif
        do i=1, ntpguest
            gstmolfract(i) = gstpress(i) / total_pressure
        enddo
        
      elseif((lpress).and.(.not.lfract))then
        total_pressure = press_sum 
        do i=1, ntpguest
            gstmolfract(i) = gstpress(i) / total_pressure
        enddo
      endif
c     NB: the way this is set up means that partial pressures
c     will override mole fractions
      if(lfuga)then
        if(idnode.eq.0)then
          write(nrite,"(/,'=============================================
     &==========')")
          write(nrite,"(/,' Calculating the fugacity with the Peng 
     &Robinsom EOS')")
          write(nrite,"(/,'=============================================
     &==========')")
        endif

c       main fugacity calculation
        do i=1,ntpguest
            call calc_fugacity(ntpguest, i, temp)
        enddo

      else
        do i=1,ntpguest
          gstfuga(i) = gstpress(i)
        enddo 
      endif
      if(idnode.eq.0)then
        write(nrite,'(a6,5x,a1,15x,a12,7x,a12)')'guest','y',
     & 'fugacity/bar','pressure/bar'
        do i=1, ntpguest
            write(nrite,'(i6,1x,f9.3,2(f19.6))') i, gstmolfract(i), 
     & gstfuga(i)/1.E5, gstpress(i)/1.E5
        enddo
      endif
      end subroutine fugacity

      subroutine guess_guest(iguess, ntpguest)
c**********************************************************************
c
c     Returns the guest type in order to extract fluid properties 
c     - currently just takes the total mass of the guest molecule..
c    
c     the guests are flagged by an integer value
c     0 - I dont know
c     1 - N2 (NB: Mass coincides with CO!)
c     2 - CO2
c     3 - CH4
c     4 - C2H6
c     5 - H2
c 
c*********************************************************************

      implicit none
      integer mol, i, j, natms,ntpguest
      integer iguess(ntpguest)
      real(8) mass

      do j=1, ntpguest
          mol = locguest(j)
          natms = numatoms(mol)
          mass=0.d0
          
          do i=1,natms
            mass=mass+atmwght(mol, i)
          enddo
c         N2
          if((mass.gt.27.0).and.(mass.lt.29.0))then
             iguess(j)=1
c         CO2
          elseif((mass.gt.43.0).and.(mass.lt.45.0))then
             iguess(j)=2
c         CH4
          elseif((mass.gt.15.5).and.(mass.lt.17.0))then
             iguess(j)=3
c         C2H6
          elseif((mass.gt.29.0).and.(mass.lt.31.0))then
             iguess(j)=4
c         H2
          elseif((mass.gt.1.0).and.(mass.lt.3.0))then
             iguess(j)=5
c         I don't know
          else
             iguess(j)=0
          endif
      enddo
      end subroutine guess_guest


      subroutine fluid_properties(idnode,ntpguest, lfuga)
c**********************************************************************
c
c     computes fluid properties based on the guest type
c     Original data found in:
c         'The properties of gases and liquids',
c         Poling, B.E.; Pausnitz, J.M.; O'Connell, J.P.
c         2001
c    
c*********************************************************************
      implicit none
      integer i, ntpguest, idnode
      logical lfuga
      integer, dimension(ntpguest) :: iguess

      call guess_guest(iguess, ntpguest) 
      if(lfuga)then
        do i = 1, ntpguest
          
          if(iguess(i).eq.1)then !case('N2')
              Acc_factor(i) = 0.037
              T_crit(i) = 126.2
              P_crit(i) = 3.398e6
          elseif(iguess(i).eq.2)then !case('CO2')
              Acc_factor(i) = 0.225
              T_crit(i) = 304.12
              P_crit(i) = 7.374e6
          elseif(iguess(i).eq.3)then!case('CH4')
              Acc_factor(i) = 0.011
              T_crit(i) = 190.56
              P_crit(i) = 4.599e6
          elseif(iguess(i).eq.4)then !case('C2H6')
              Acc_factor(i) = 0.099
              T_crit(i) = 305.32
              P_crit(i) = 4.872e6
          elseif(iguess(i).eq.5)then !case('H2')
              Acc_factor(i) = -0.217
              T_crit(i) = 32.98
              P_crit(i) = 1.293e6
          elseif(iguess(i).eq.0)then !case(default)
              if(idnode.eq.0)write(nrite, "(/,i4,a)")i,' guest was not 
     &found in the fluid properties, setting fugacities equal to 
     &pressures found in the CONTROL file'
              lfuga=.false.
          endif
        enddo
        if (ntpguest.eq.2)then
c           CH4 & C2H6
            if((iguess(1).eq.3).and.(iguess(2).eq.4).or.
     &          (iguess(1).eq.4).and.(iguess(2).eq.3))then
                K_fug(1,2) = -0.003
            else
                write(nrite, "(/,a,/,a)")'Binary mixture found, but 
     &data not available for guests. Assuming a mixing 
     &coefficient of 0.','If data is found, please edit the 
     &fluid_properties subroutine in utility_pack.f'
            endif
            K_fug(2,1) = K_fug(1,2)
        endif

      endif
      end subroutine fluid_properties
      subroutine calc_fugacity(ntpguest, iguest, temp)
c**********************************************************************
c
c     computes fugacity using Peng-Robinson EOS
c     
c    
c*********************************************************************
      implicit none
      real(8) fuga_coeff
      real(8) alpha, beta_f, gamma
      real(8), dimension(ntpguest) :: kappa, alphaT, agreek, bgreek
      real(8) agreek_sys, bgreek_sys
      real(8) A, B
      real(8) Z0, Z1, temp
      real(8) sum_za
      real(8) R, crit_conv
      integer i,j, ntpguest
      integer iguest 

      R = 8.314
      crit_conv = 1.d-5
      do i = 1, ntpguest
        kappa(i) =  0.37464 + 1.54226*acc_factor(i) -
     &              0.26992*acc_factor(i)*acc_factor(i)
        alphaT(i) = (1 + kappa(i) * (1 - sqrt(temp/T_crit(i))))**2
c       pure components
        agreek(i) = 0.45724 * (R * T_crit(i))**2 / P_crit(i)
     &* alphaT(i)
        bgreek(i) = 0.0778 * R * T_crit(i) / P_crit(i)
      enddo
      agreek_sys = 0.
      bgreek_sys = 0.
      if (ntpguest == 1) then
        agreek_sys = agreek(1)
        bgreek_sys = bgreek(1)
      else
        do i = 1, ntpguest 
          do j = 1, ntpguest
            agreek_sys = agreek_sys +
     &          gstmolfract(i)*gstmolfract(j)*sqrt(agreek(i)*agreek(j))
     &          *(1-K_fug(i,j))
          enddo
          bgreek_sys = bgreek_sys + gstmolfract(i)*bgreek(i)
        enddo 
      endif 

      A = agreek_sys * gstpress(iguest) / (R * temp)**2
      B = bgreek_sys * gstpress(iguest) / R / temp
      alpha = B - 1
      beta_f = A - 3*B*B - 2*B
      gamma = -A*B + B*B + B*B*B
      
c     Searching for zero using Newton-Rhapson
c     Z**3 + alpha*Z**2 + beta_f*Z + gamma = 0
c     Z: compressibility facor
      
      Z1 = 1.0
      do
         Z0 = Z1
         Z1 = Z0 - (Z0**3 + alpha*Z0**2 + beta_f*Z0 + gamma) / 
     &        (3*Z0**2 + 2*alpha*Z0 +beta_f)
         if (abs(Z1 - Z0) < crit_conv)  exit
      enddo
     
c     calculate the sum(aik*zi)
      sum_za = 0.0
      do i = 1, ntpguest
         sum_za = sum_za + gstmolfract(i)*sqrt(agreek(iguest)
     &           *agreek(i))*(1-K_fug(iguest,i))
      enddo 
      
c     calculate fugacity coefficient
      
      fuga_coeff = exp(bgreek(iguest)/bgreek_sys*(Z1-1) - log(Z1 - B) 
     &             - agreek_sys/(2.*SQRT(2.0)*bgreek_sys*R*temp)*
     &              (2*sum_za/agreek_sys-bgreek(iguest)/bgreek_sys) 
     &              * log( (Z1+(1+sqrt(2.0))*B)/(Z1+(1-sqrt(2.0))*B)))
      
      gstfuga(iguest) = fuga_coeff *gstmolfract(iguest) * 
     &                    gstpress(iguest)

      end subroutine calc_fugacity

      subroutine floorimages
     x  (imcon,natms,cell,xxx,yyy,zzz)
      
c***********************************************************************
c     this will adjust the coordinates to fit within the original
c     box.
c
c     for
c     imcon=0 no boundary conditions apply
c     imcon=1 standard cubic boundaries apply
c     imcon=2 orthorhombic boundaries apply
c     imcon=3 parallelepiped boundaries apply
c     imcon=4 truncated octahedron boundaries apply
c     imcon=5 rhombic dodecahedron boundaries apply
c     imcon=6 x-y parallelogram boundary conditions : no periodicity in z
c     imcon=7 hexagonal prism boundaries apply
c     
c     
c***********************************************************************
      
      implicit none

      integer imcon,natms,i
      real(8) aaa,bbb,ccc,det,rt2,rt3,ssx
      real(8) ssy,ssz,ddd,xss,yss,zss

      real(8), dimension(natms):: xxx,yyy,zzz
      real(8), dimension(9) :: cell,rcell

      data rt2/1.41421356623d0/,rt3/1.7320508075d0/
      if(imcon.eq.1)then

c     standard cubic boundary conditions
        
        
        aaa=1.d0/cell(1)

        do i=1,natms
          xxx(i)=xxx(i)-cell(1)*floor(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*floor(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*floor(aaa*zzz(i))
        enddo
        
      else if(imcon.eq.2)then

c     rectangular (slab) boundary conditions
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(5)
        ccc=1.d0/cell(9)
        
        do i=1,natms
          
          xxx(i)=xxx(i)-cell(1)*floor(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(5)*floor(bbb*yyy(i))
          zzz(i)=zzz(i)-cell(9)*floor(ccc*zzz(i))
          
        enddo
        
      else if(imcon.eq.3)then

c     parallelpiped boundary conditions

        call invert(cell,rcell,det)
        do i=1,natms
          ssx=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))
          ssy=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))
          ssz=(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i))
          
          xss=ssx-floor(ssx)
          yss=ssy-floor(ssy)
          zss=ssz-floor(ssz)
          
          xxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
          yyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
          zzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
          
        enddo
      else if(imcon.eq.4)then

c     truncated octahedral boundary conditions
        
c        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
c     x    abs(cell(5)-cell(9)).lt.1.d-6)) call error(idnode,130)
        
        aaa=1.d0/cell(1)
        
        do i=1,natms
          
          xxx(i)=xxx(i)-cell(1)*floor(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*floor(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*floor(aaa*zzz(i))
          
          if((abs(xxx(i))+abs(yyy(i))+abs(zzz(i))).ge.
     x      (0.75d0*cell(1)))then
            
            xxx(i)=xxx(i)-0.5d0*sign(cell(1),xxx(i))
            yyy(i)=yyy(i)-0.5d0*sign(cell(1),yyy(i))
            zzz(i)=zzz(i)-0.5d0*sign(cell(1),zzz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.5)then

c     rhombic dodecahedral boundary conditions
        
c        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
c     x    abs(cell(9)-cell(1)*rt2).lt.1.d-6)) 
c     x    call error(idnode,140)
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(9)
        
        do i=1,natms
          
          xxx(i)=xxx(i)-cell(1)*floor(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*floor(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(9)*floor(bbb*zzz(i))
          
          if((abs(xxx(i))+abs(yyy(i))+abs(rt2*zzz(i))).ge.
     x      cell(1))then
            
            xxx(i)=xxx(i)-0.5d0*sign(cell(1),xxx(i))
            yyy(i)=yyy(i)-0.5d0*sign(cell(1),yyy(i))
            zzz(i)=zzz(i)-0.5d0*sign(cell(9),zzz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.6) then

c     x-y boundary conditions 

        det = cell(1)*cell(5) - cell(2)*cell(4)

c        if(abs(det).lt.1.d-6)call error(idnode,120)
        
        det = 1.d0/det

        rcell(1) =  det*cell(5)
        rcell(2) = -det*cell(2)
        rcell(4) = -det*cell(4)
        rcell(5) =  det*cell(1)
        
        do i=1,natms

          ssx = rcell(1)*xxx(i) + rcell(4)*yyy(i)
          ssy = rcell(2)*xxx(i) + rcell(5)*yyy(i)

          xss = ssx - floor(ssx)
          yss = ssy - floor(ssy)

          xxx(i)=cell(1)*xss + cell(4)*yss
          yyy(i)=cell(2)*xss + cell(5)*yss

        enddo

      else if(imcon.eq.7) then

c     hexagonal prism boundary conditions
        
c        if(abs(cell(1)-rt3*cell(5)).ge.1.d-6)
c     x    call error(idnode,135)
        
        aaa=cell(1)/(rt3*2.d0)
        bbb=cell(1)/rt3
        ccc=rt3/cell(1)
        ddd=1.d0/cell(9)
        
        do i=1,natms
          
          yyy(i)=yyy(i)-bbb*floor(ccc*yyy(i))
          zzz(i)=zzz(i)-cell(9)*floor(ddd*zzz(i))
          
          if((abs(yyy(i))+abs(rt3*xxx(i))).ge.bbb)then
            
            xxx(i)=xxx(i)-rt3*sign(aaa,xxx(i))
            yyy(i)=yyy(i)-sign(aaa,yyy(i))
            
          endif
          
        enddo
        
      endif
      return
      end subroutine floorimages
      subroutine images
     x  (imcon,natms,cell,sxxx,syyy,szzz)
      
c***********************************************************************
c     
c     subroutine for calculating the minimum image
c     of atom pairs within a specified MD cell
c
c     note the following internal units apply everywhere
c     
c     unit of time      (to)    =          1 x 10**(-12) seconds
c     unit of length    (lo)    =          1 x 10**(-10) metres
c     unit of mass      (mo)    = 1.6605402  x 10**(-27) kilograms
c     unit of charge    (qo)    = 1.60217733 x 10**(-19) coulombs
c     unit of energy    (eo)    = 1.6605402  x 10**(-23) joules
c     unit of pressure  (po)    = 1.6605402  x 10**(  7) pascals
c     
c     for
c     imcon=0 no boundary conditions apply
c     imcon=1 standard cubic boundaries apply
c     imcon=2 orthorhombic boundaries apply
c     imcon=3 parallelepiped boundaries apply
c     imcon=4 truncated octahedron boundaries apply
c     imcon=5 rhombic dodecahedron boundaries apply
c     imcon=6 x-y parallelogram boundary conditions : no periodicity in z
c     imcon=7 hexagonal prism boundaries apply
c     
c     note: in all cases the centre of the cell is at (0,0,0)
c     warning - replicated data version: does not re-merge 
c     coordinate arrays
c     
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c     
c***********************************************************************
      
      implicit none

      integer imcon,natms,i
      real(8) cell,sxxx,syyy,szzz,aaa,bbb,ccc,det,rt2,rt3,ssx
      real(8) ssy,ssz,ddd,xss,yss,zss,rcell

      dimension sxxx(*),syyy(*),szzz(*)
      dimension cell(9),rcell(9)

      data rt2/1.41421356623d0/,rt3/1.7320508075d0/
      if(imcon.eq.1)then

c     standard cubic boundary conditions
        
        
        aaa=1.d0/cell(1)

        do i=1,natms
          sxxx(i)=sxxx(i)-cell(1)*nint(aaa*sxxx(i))
          syyy(i)=syyy(i)-cell(1)*nint(aaa*syyy(i))
          szzz(i)=szzz(i)-cell(1)*nint(aaa*szzz(i))
        enddo
        
      else if(imcon.eq.2)then

c     rectangular (slab) boundary conditions
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(5)
        ccc=1.d0/cell(9)
        
        do i=1,natms
          
          sxxx(i)=sxxx(i)-cell(1)*nint(aaa*sxxx(i))
          syyy(i)=syyy(i)-cell(5)*nint(bbb*syyy(i))
          szzz(i)=szzz(i)-cell(9)*nint(ccc*szzz(i))
          
        enddo
        
      else if(imcon.eq.3)then

c     parallelpiped boundary conditions
        call invert(cell,rcell,det)
        do i=1,natms
          ssx=(rcell(1)*sxxx(i)+rcell(4)*syyy(i)+rcell(7)*szzz(i))
          ssy=(rcell(2)*sxxx(i)+rcell(5)*syyy(i)+rcell(8)*szzz(i))
          ssz=(rcell(3)*sxxx(i)+rcell(6)*syyy(i)+rcell(9)*szzz(i))
          
          xss=ssx-nint(ssx)
          yss=ssy-nint(ssy)
          zss=ssz-nint(ssz)
          
          sxxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
          syyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
          szzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
          
        enddo
      else if(imcon.eq.4)then

c     truncated octahedral boundary conditions
        
c        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
c     x    abs(cell(5)-cell(9)).lt.1.d-6)) call error(idnode,130)
        
        aaa=1.d0/cell(1)
        
        do i=1,natms
          
          sxxx(i)=sxxx(i)-cell(1)*nint(aaa*sxxx(i))
          syyy(i)=syyy(i)-cell(1)*nint(aaa*syyy(i))
          szzz(i)=szzz(i)-cell(1)*nint(aaa*szzz(i))
          
          if((abs(sxxx(i))+abs(syyy(i))+abs(szzz(i))).ge.
     x      (0.75d0*cell(1)))then
            
            sxxx(i)=sxxx(i)-0.5d0*sign(cell(1),sxxx(i))
            syyy(i)=syyy(i)-0.5d0*sign(cell(1),syyy(i))
            szzz(i)=szzz(i)-0.5d0*sign(cell(1),szzz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.5)then

c     rhombic dodecahedral boundary conditions
        
c        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
c     x    abs(cell(9)-cell(1)*rt2).lt.1.d-6)) 
c     x    call error(idnode,140)
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(9)
        
        do i=1,natms
          
          sxxx(i)=sxxx(i)-cell(1)*nint(aaa*sxxx(i))
          syyy(i)=syyy(i)-cell(1)*nint(aaa*syyy(i))
          szzz(i)=szzz(i)-cell(9)*nint(bbb*szzz(i))
          
          if((abs(sxxx(i))+abs(syyy(i))+abs(rt2*szzz(i))).ge.
     x      cell(1))then
            
            sxxx(i)=sxxx(i)-0.5d0*sign(cell(1),sxxx(i))
            syyy(i)=syyy(i)-0.5d0*sign(cell(1),syyy(i))
            szzz(i)=szzz(i)-0.5d0*sign(cell(9),szzz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.6) then

c     x-y boundary conditions 

        det = cell(1)*cell(5) - cell(2)*cell(4)

c        if(abs(det).lt.1.d-6)call error(idnode,120)
        
        det = 1.d0/det

        rcell(1) =  det*cell(5)
        rcell(2) = -det*cell(2)
        rcell(4) = -det*cell(4)
        rcell(5) =  det*cell(1)
        
        do i=1,natms

          ssx = rcell(1)*sxxx(i) + rcell(4)*syyy(i)
          ssy = rcell(2)*sxxx(i) + rcell(5)*syyy(i)

          xss = ssx - nint(ssx)
          yss = ssy - nint(ssy)

          sxxx(i)=cell(1)*xss + cell(4)*yss
          syyy(i)=cell(2)*xss + cell(5)*yss

        enddo

      else if(imcon.eq.7) then

c     hexagonal prism boundary conditions
        
c        if(abs(cell(1)-rt3*cell(5)).ge.1.d-6)
c     x    call error(idnode,135)
        
        aaa=cell(1)/(rt3*2.d0)
        bbb=cell(1)/rt3
        ccc=rt3/cell(1)
        ddd=1.d0/cell(9)
        
        do i=1,natms
          
          syyy(i)=syyy(i)-bbb*nint(ccc*syyy(i))
          szzz(i)=szzz(i)-cell(9)*nint(ddd*szzz(i))
          
          if((abs(syyy(i))+abs(rt3*sxxx(i))).ge.bbb)then
            
            sxxx(i)=sxxx(i)-rt3*sign(aaa,sxxx(i))
            syyy(i)=syyy(i)-sign(aaa,syyy(i))
            
          endif
          
        enddo
        
      endif
      return
      end subroutine images

      subroutine dcell(aaa,bbb)

c***********************************************************************
c     
c     subroutine to calculate the dimensional properties of
c     a simulation cell specified by the input matrix aaa.
c     the results are returned in the array bbb, with :
c     
c     bbb(1 to 3) - lengths of cell vectors
c     bbb(4 to 6) - cosines of cell angles
c     bbb(7 to 9) - perpendicular cell widths
c     bbb(10)     - cell volume
c     
c***********************************************************************

      implicit none

      real(8) aaa,bbb,axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3

      dimension aaa(9),bbb(10)

c     calculate lengths of cell vectors

      bbb(1)=sqrt(aaa(1)*aaa(1)+aaa(2)*aaa(2)+aaa(3)*aaa(3))
      bbb(2)=sqrt(aaa(4)*aaa(4)+aaa(5)*aaa(5)+aaa(6)*aaa(6))
      bbb(3)=sqrt(aaa(7)*aaa(7)+aaa(8)*aaa(8)+aaa(9)*aaa(9))

c     calculate cosines of cell angles

      bbb(4)=(aaa(1)*aaa(4)+aaa(2)*aaa(5)+aaa(3)*aaa(6))/(bbb(1)*bbb(2))
      bbb(5)=(aaa(1)*aaa(7)+aaa(2)*aaa(8)+aaa(3)*aaa(9))/(bbb(1)*bbb(3))
      bbb(6)=(aaa(4)*aaa(7)+aaa(5)*aaa(8)+aaa(6)*aaa(9))/(bbb(2)*bbb(3))

c     calculate vector products of cell vectors

      axb1=aaa(2)*aaa(6)-aaa(3)*aaa(5)
      axb2=aaa(3)*aaa(4)-aaa(1)*aaa(6)
      axb3=aaa(1)*aaa(5)-aaa(2)*aaa(4)
      bxc1=aaa(5)*aaa(9)-aaa(6)*aaa(8)
      bxc2=aaa(6)*aaa(7)-aaa(4)*aaa(9)
      bxc3=aaa(4)*aaa(8)-aaa(5)*aaa(7)
      cxa1=aaa(8)*aaa(3)-aaa(2)*aaa(9)
      cxa2=aaa(1)*aaa(9)-aaa(3)*aaa(7)
      cxa3=aaa(2)*aaa(7)-aaa(1)*aaa(8)

c     calculate volume of cell

      bbb(10)=abs(aaa(1)*bxc1+aaa(2)*bxc2+aaa(3)*bxc3)

c     calculate cell perpendicular widths

      bbb(7)=bbb(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
      bbb(8)=bbb(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
      bbb(9)=bbb(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)

      return
      end subroutine dcell

      subroutine invert(a,b,d)

c***********************************************************************
c     
c     subroutine to invert a 3 * 3 matrix using cofactors
c     
c***********************************************************************

      implicit none

      real(8) a,b,d,r

      dimension a(9),b(9)

c     calculate adjoint matrix
      b(1)=a(5)*a(9)-a(6)*a(8)
      b(2)=a(3)*a(8)-a(2)*a(9)
      b(3)=a(2)*a(6)-a(3)*a(5)
      b(4)=a(6)*a(7)-a(4)*a(9)
      b(5)=a(1)*a(9)-a(3)*a(7)
      b(6)=a(3)*a(4)-a(1)*a(6)
      b(7)=a(4)*a(8)-a(5)*a(7)
      b(8)=a(2)*a(7)-a(1)*a(8)
      b(9)=a(1)*a(5)-a(2)*a(4)

c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d

c     complete inverse matrix
      b(1)=r*b(1)
      b(2)=r*b(2)
      b(3)=r*b(3)
      b(4)=r*b(4)
      b(5)=r*b(5)
      b(6)=r*b(6)
      b(7)=r*b(7)
      b(8)=r*b(8)
      b(9)=r*b(9)

      return
      end subroutine invert

      subroutine jacobi(a,v,n)

c***********************************************************************
c     
c     diagonalisation of real symmetric matices by jacobi method
c     
c     input parameters:
c     
c     a(n,n) is the matrix to be diagonalised
c     v(n,n) is the eigenvector matrix
c     n   is the dimension of the matrices
c     
c     jacobi processes lower triangle only (upper triangle unchanged)
c     
c     variable rho sets absolute tolerance on convergence
c     variable tes is a moving tolerance that diminishes
c     on each pass until at true convergence tes<rho
c     
c     author w.smith 1993
c     
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c     
c***********************************************************************

      implicit none

      logical pass
      integer n,i,j,k
      real(8) a,v,rho,tes,scl,v1,v2,v3,u,omg,s,c,tem

      dimension a(n,n),v(n,n)

      rho=1.0d-16
      tes=0.0d0
      scl=0.0d0

c     initialize eigenvectors

      do i=1,n
        do j=1,n
          v(i,j)=0.0d0
        enddo
        v(i,i)=1.0d0
      enddo

c     rescale matrix for optimal accuracy

      do i=1,n
        if(abs(a(i,i)).gt.scl)scl=abs(a(i,i))
      enddo
      do i=1,n
        do j=1,i
          a(i,j)=a(i,j)/scl
        enddo
      enddo

c     set initial value of moving tolerance

      do i=2,n
        do j=1,i-1
          tes=tes+2.0d0*a(i,j)*a(i,j)
        enddo
      enddo
      tes=sqrt(tes)

c     recycle until absolute tolerance satisfied

      do while(tes.gt.rho)

        tes=tes/dble(n)
        if(tes.lt.rho)tes=rho
        
c     jacobi diagonalisation
        
        pass=.true.
        
c     recycle until moving tolerance satisfied
        
        do while(pass)
          
          pass=.false.
          
          do i=2,n
            
            do j=1,i-1
              
              if(abs(a(i,j)).ge.tes)then
                pass=.true.
                v1=a(j,j)
                v2=a(i,j)
                v3=a(i,i)
                u=0.5d0*(v1-v3)
                if(abs(u).lt.rho)then
                  omg=-1.0d0
                else
                  omg=-v2/sqrt(v2*v2+u*u)
                  if(u.lt.0.0d0)omg=-omg
                endif
                s=omg/sqrt(2.0d0*(1.0d0+sqrt(1.0d0-omg*omg)))
                c=sqrt(1.0d0-s*s)
                do k=1,n
                  if(k.ge.i)then
                    tem=a(k,j)*c-a(k,i)*s
                    a(k,i)=a(k,j)*s+a(k,i)*c
                    a(k,j)=tem
                  else if(k.lt.j)then
                    tem=a(j,k)*c-a(i,k)*s
                    a(i,k)=a(j,k)*s+a(i,k)*c
                    a(j,k)=tem
                  else
                    tem=a(k,j)*c-a(i,k)*s
                    a(i,k)=a(k,j)*s+a(i,k)*c
                    a(k,j)=tem
                  endif
                  tem=v(k,j)*c-v(k,i)*s
                  v(k,i)=v(k,j)*s+v(k,i)*c
                  v(k,j)=tem
                enddo
                a(j,j)=v1*c*c+v3*s*s-2.0d0*v2*s*c
                a(i,i)=v1*s*s+v3*c*c+2.0d0*v2*s*c
                a(i,j)=(v1-v3)*s*c+v2*(c*c-s*s)
              endif
              
            enddo
            
          enddo
          
        enddo

      enddo

c     rescale matrix

      do i=1,n
        do j=1,i
          a(i,j)=scl*a(i,j)
        enddo
      enddo

      return
      end subroutine jacobi

      subroutine set_block(nnn,ccc,aaa)

c**********************************************************************
c
c     subroutine to initialise an array to a single value
c
c**********************************************************************

      implicit none

      integer i,nnn
      real(8) ccc,aaa(nnn)

      do i=1,nnn,2

        aaa(i)=ccc
        aaa(i+1)=ccc

      enddo
      
      return
      end subroutine set_block

      subroutine matmul(aaa,bbb,ccc)

c***********************************************************************
c     
c     utility to multiply 3x3 matrices
c
c**********************************************************************

      implicit none

      integer i
      real(8) aaa(9),bbb(9),ccc(9),tmp(9)

      tmp(1)=aaa(1)*bbb(1)+aaa(4)*bbb(2)+aaa(7)*bbb(3)
      tmp(2)=aaa(2)*bbb(1)+aaa(5)*bbb(2)+aaa(8)*bbb(3)
      tmp(3)=aaa(3)*bbb(1)+aaa(6)*bbb(2)+aaa(9)*bbb(3)

      tmp(4)=aaa(1)*bbb(4)+aaa(4)*bbb(5)+aaa(7)*bbb(6)
      tmp(5)=aaa(2)*bbb(4)+aaa(5)*bbb(5)+aaa(8)*bbb(6)
      tmp(6)=aaa(3)*bbb(4)+aaa(6)*bbb(5)+aaa(9)*bbb(6)

      tmp(7)=aaa(1)*bbb(7)+aaa(4)*bbb(8)+aaa(7)*bbb(9)
      tmp(8)=aaa(2)*bbb(7)+aaa(5)*bbb(8)+aaa(8)*bbb(9)
      tmp(9)=aaa(3)*bbb(7)+aaa(6)*bbb(8)+aaa(9)*bbb(9)
      
      do i=1,9
        ccc(i)=tmp(i)
      enddo
      
      return
      end subroutine matmul 

      subroutine debugging(idnode,natms,levcfg,imcon,cfgname,
     &eng,outdir,ins,del,dis,pass,delE,ewld1eng,ewld2sum,ewld3sum,
     &vdwsum,delrc,elrc,engunit,iguest)
c*********************************************************************
c     subroutine writes files necessary to get DL_poly energies
c     from insertions,deletions and displacements
c*********************************************************************
      implicit none
      logical ins,del,dis
      character*1 cfgname(80)
      character*8 outdir
      real(8) engunit,ewld1eng,ewld2sum,ewld3sum,vdwsum
      real(8) delrc,elrc
      real(8) eng,delE(*)
      integer pass,iguest,levcfg,imcon,idnode,natms

      idnode=idnode
      if(pass.eq.1)then
        call revive_debug
     &(natms,levcfg,imcon,cfgname,eng,outdir,ins,del,dis,pass)
      else

        call revive_debug
     &(natms,levcfg,imcon,cfgname,eng,outdir,ins,del,dis,pass)
        if(ins)then          
          open(39,file=outdir//'/debug_ins')
          write(39,'(3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a13,f20.10,/,
     &      3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a13,f20.10)') 
     &      "Total eng: ",delE(iguest),"ewald 1 sum:",ewld1eng,
     &      "ewald 2 sum: ",ewld2sum,"ewald 3 sum: ",ewld3sum,
     &      "sum of vdw: ",vdwsum,
     &      "lrng correct:",delrc/engunit
          close(39)
        elseif(del)then
          open(39,file=outdir//'/debug_del')
           
          write(39,'(3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a13,f20.10,/,
     &      3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a13,f20.10)') 
     &      "Total eng: ",delE(iguest),"ewald 1 sum:",ewld1eng,
     &      "ewald 2 sum: ",ewld2sum,"ewald 3 sum: ",ewld3sum,
     &      "sum of vdw: ",vdwsum,
     &      "lrng correct:",delrc/engunit
          close(39) 
        elseif(dis)then
          open(39,file=outdir//'/debug_dis')
           
          write(39,'(3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a13,f20.10,/,
     &      3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a19,f20.10)') 
     &      "Total eng: ",delE(iguest),"ewald 1 sum:",ewld1eng,
     &      "ewald 2 sum: ",ewld2sum,"ewald 3 sum: ",ewld3sum,
     &      "sum of vdw: ",vdwsum,
     &      "lrng correct total:",elrc/engunit
          close(39) 
        endif
      endif
      end subroutine debugging

      subroutine revive_debug
     &(natms,levcfg,imcon,cfgname,eng,outdir,ins,del,dis,pass)
c**********************************************************************
c 
c     subroutine to write the atom information to a file
c
c**********************************************************************
      implicit none
      integer, parameter :: nconfig=23
      logical ins,del,dis
      character*1 cfgname(80)
      character*8 outdir
      character*25 outfile
      integer i,natms,levcfg,imcon,pass
      real(8) eng


      if(ins)
     &write(outfile,'(a8,a11,i1)')outdir,'/REVCON_ins',pass
      if(del)
     &write(outfile,'(a8,a11,i1)')outdir,'/REVCON_del',pass
      if(dis)
     &write(outfile,'(a8,a11,i1)')outdir,'/REVCON_dis',pass
    
      open(nconfig,file=outfile,form='formatted')
      write(nconfig,'(80a1)')cfgname
      write(nconfig,'(2i10,e20.10)')levcfg,imcon,eng
      if(imcon.gt.0)write(nconfig,'(3f20.12)')cell

      do i=1,natms
        write(nconfig,'(a8,i10)')atomname(i),i
        write(nconfig,'(3g20.10)')xxx(i),yyy(i),zzz(i)
      enddo
      close(nconfig)
      return
      end subroutine revive_debug
      subroutine revive
     &(natms,levcfg,production,ntpguest,ntpmls,
     &imcon,cfgname,eng,outdir)
c**********************************************************************
c 
c     subroutine to write the atom information to a file
c
c**********************************************************************
      implicit none
      logical production
      character*1 cfgname(80)
      character*8 outdir
      character*25 outfile
      integer i,j,natms,levcfg,imcon,ntpguest,ntpmls,nsites
      real(8) eng

      write(outfile,'(a8,a7)')outdir,'/REVCON'
      open(nconfig,file=outfile,form='formatted')
      write(nconfig,'(80a1)')cfgname
      write(nconfig,'(2i10,e20.10)')levcfg,imcon,eng
      if(imcon.gt.0)write(nconfig,'(3f20.12)')cell

      do i=1,natms
        write(nconfig,'(a8,i10)')atomname(i),i
        write(nconfig,'(3g20.10)')xxx(i),yyy(i),zzz(i)
      enddo
      close(nconfig)

      write(outfile,'(a8,a7)')outdir,'/REVIVE'
      open(nrev,file=outfile,form='unformatted')
      write(nrev)production
      write(nrev)(nummols(i),i=1,ntpmls)
      write(nrev)(chainstats(i),i=1,ntpguest*16+1)
      write(nrev)natms
      do i=1,ntpmls
        nsites=nummols(i)*numatoms(i)
        write(nrev)(molxxx(i,j),j=1,nsites)
        write(nrev)(molyyy(i,j),j=1,nsites)
        write(nrev)(molzzz(i,j),j=1,nsites)
      enddo
      
   
      write(nrev)(energy(i),i=1,ntpguest)
      close(nrev)
      return
      end subroutine revive

      subroutine revscan(idnode,prevnodes)
c*********************************************************************
c
c     determine how many nodes were done in a previous run 
c
c*********************************************************************
      implicit none
      logical safe
      integer idnode,prevnodes,idum

      safe=.true.
      if(idnode.eq.0)then
         call system("ls | grep -c branch > 0110111")
         open(666,file="0110111")
      endif
      call getrec(safe,idnode,666)
      if(idnode.eq.0)close(666)
      prevnodes=intstr(record,lenrec,idum)
      if(idnode.eq.0)call system("rm 0110111")

      end subroutine revscan
      subroutine revread(localdir,production,ntpmls,totatm,ntpguest)
c*********************************************************************
c
c     read revive file and revcon file to get all the related data
c     to restart a gcmc simulation
c
c*********************************************************************
      implicit none
      logical production
      character*8 localdir
      character*15 outfile
      integer i,j,nsites,totatm,ntpguest,ios,ntpmls


      write(outfile,'(a8,a7)')localdir,'/REVIVE'
      open(nrev,file=outfile,form='unformatted',status='old',iostat=ios)

      if(ios.eq.29)then
c       value 29 means status=old is invalid ie. there is no existing
c       file called REVIVE.  Thus the restart will just start with the 
c       conditions in the CONFIG and FIELD file.

      else
        read(nrev)production
        read(nrev)(nummols(i),i=1,ntpmls)
        read(nrev)(chainstats(i),i=1,ntpguest*16+1)
        read(nrev)totatm
        do i=1,ntpmls
          nsites=nummols(i)*numatoms(i)
          read(nrev)(molxxx(i,j),j=1,nsites)
          read(nrev)(molyyy(i,j),j=1,nsites)
          read(nrev)(molzzz(i,j),j=1,nsites)
        enddo
   
        read(nrev)(energy(i),i=1,ntpguest)
        close(nrev)
      endif

      end subroutine revread
      
      subroutine equitable_bin(fposa,fposb,fposc,ngrida,ngridb,ngridc,
     &subgrid,val)
c***********************************************************************
c
c     update the global grid(subgrid,...) with the equitable binned 
c     atom position
c     
c**********************************************************************
      implicit none
c arguments
      integer ngrida,ngridb,ngridc,subgrid
      real(8) fposa,fposb,fposc
c internal variables
      integer bin,aidx,bidx,cidx
      integer alidx,aridx,blidx,bridx,clidx,cridx
      real(8) apart,bpart,cpart,gpta,gptb,gptc
      real(8) al,ar,bl,br,cl,cr,val

c NOTE: since fortran arrays start at 1 and probabilities are 1D
c a?idx and b?idx are 0 based index and c?idx is 1 based
      gpta = ngrida*(fposa)
      aidx = ceiling(gpta)
      apart = modulo(gpta, 1.d0)
      if(apart.gt.0.5)then
        if(aidx.ge.ngrida)then
          alidx = (ngrida-1)*ngridb*ngridc
          aridx = 0
        else
          alidx = (aidx-1)*ngridb*ngridc
          aridx = (aidx)*ngridb*ngridc
        endif
        al = 1.5 - apart
        ar = apart - 0.5
      else
        if(aidx.le.1)then
          alidx = (ngrida-1)*ngridb*ngridc
          aridx = 0
        else
          alidx = (aidx-2)*ngridb*ngridc
          aridx = (aidx-1)*ngridb*ngridc
        endif
        al = 0.5 - apart
        ar = 0.5 + apart
      endif
c equitable spread for b vector
      gptb = ngridb*(fposb)
      bidx = ceiling(gptb)
      bpart = modulo(gptb, 1.d0)
      if(bpart.gt.0.5)then
        if(bidx.ge.ngridb)then
          blidx = (ngridb-1)*ngridc
          bridx = 0
        else
          blidx = (bidx-1)*ngridc
          bridx = (bidx)*ngridc
        endif
        bl = 1.5 - bpart
        br = bpart - 0.5
      else
        if(bidx.le.1)then
          blidx = (ngridb-1)*ngridc
          bridx = 0
        else
          blidx = (bidx-2)*ngridc
          bridx = (bidx-1)*ngridc
        endif
        bl = 0.5 - bpart
        br = 0.5 + bpart
      endif
c equitable spread for c vector
      gptc = ngridc*(fposc)
      cidx = ceiling(gptc)
      cpart = modulo(gptc, 1.d0)
      if(cpart.gt.0.5)then
        if(cidx.ge.ngridc)then
          clidx = ngridc
          cridx = 1
        else
          clidx = cidx
          cridx = cidx + 1
        endif
        cl = 1.5 - cpart
        cr = cpart - 0.5
      else
        if(cidx.le.1)then
          clidx = ngridc
          cridx = 1
        else
          clidx = cidx-1
          cridx = cidx
        endif
        cl = 0.5 - cpart
        cr = 0.5 + cpart
      endif
      bin = alidx+blidx+clidx
      grid(subgrid,bin) = grid(subgrid,bin)+(al*bl*cl)*val
      bin = aridx+blidx+clidx
      grid(subgrid,bin) = grid(subgrid,bin)+(ar*bl*cl)*val
      bin = alidx+bridx+clidx
      grid(subgrid,bin) = grid(subgrid,bin)+(al*br*cl)*val
      bin = alidx+blidx+cridx
      grid(subgrid,bin) = grid(subgrid,bin)+(al*bl*cr)*val
      bin = aridx+bridx+clidx
      grid(subgrid,bin) = grid(subgrid,bin)+(ar*br*cl)*val
      bin = aridx+blidx+cridx
      grid(subgrid,bin) = grid(subgrid,bin)+(ar*bl*cr)*val
      bin = alidx+bridx+cridx
      grid(subgrid,bin) = grid(subgrid,bin)+(al*br*cr)*val
      bin = aridx+bridx+cridx
      grid(subgrid,bin) = grid(subgrid,bin)+(ar*br*cr)*val
      end subroutine equitable_bin

      integer function atmnumber(i) 
c*******************************************************************
c     generate an atomic number from an atomic mass 
c     EDIT (pb 09/01/13): this function reads the mass reported
c     on a standard periodic table and assigns an atomic number.
c     You will run into problems if you are using atomic masses
c     of isotopes in the FIELD file.
c*******************************************************************
      implicit none
      real(8) i
      if ((i.ge.0.0).and.(i.le.1.5))then
        atmnumber=1
      elseif((i.ge.3.9).and.(i.le.4.5))then
        atmnumber=2
      elseif((i.ge.6.5).and.(i.le.7.1))then
        atmnumber=3
      elseif((i.ge.8.9).and.(i.le.9.5))then
        atmnumber=4
      elseif((i.ge.10.5).and.(i.le.11.1))then
        atmnumber=5
      elseif((i.ge.11.9).and.(i.le.12.5))then
        atmnumber=6
      elseif((i.ge.13.9).and.(i.le.14.5))then
        atmnumber=7
      elseif((i.ge.15.5).and.(i.le.16.1))then
        atmnumber=8
      elseif((i.ge.18.5).and.(i.le.19.1))then
        atmnumber=9
      elseif((i.ge.19.9).and.(i.le.20.5))then
        atmnumber=10
      elseif((i.ge.22.5).and.(i.le.23.1))then
        atmnumber=11
      elseif((i.ge.23.9).and.(i.le.24.5))then
        atmnumber=12
      elseif((i.ge.26.5).and.(i.le.27.1))then
        atmnumber=13
      elseif((i.ge.27.9).and.(i.le.28.5))then
        atmnumber=14
      elseif((i.ge.30.5).and.(i.le.31.1))then
        atmnumber=15
      elseif((i.ge.31.9).and.(i.le.32.5))then
        atmnumber=16
      elseif((i.ge.34.9).and.(i.le.36.1))then
        atmnumber=17
c     Ar (18) has mass range that overlaps with Ca (20). Be careful 
c     with mass rounding here!
      elseif((i.ge.39.5).and.(i.le.39.9999))then
        atmnumber=18
      elseif((i.ge.38.9).and.(i.le.39.4))then
        atmnumber=19
      elseif((i.ge.40.0).and.(i.le.40.5))then
        atmnumber=20
      elseif((i.ge.44.5).and.(i.le.45.1))then
        atmnumber=21
      elseif((i.ge.47.5).and.(i.le.48.1))then
        atmnumber=22
      elseif((i.ge.50.5).and.(i.le.51.1))then
        atmnumber=23
      elseif((i.ge.51.5).and.(i.le.52.1))then
        atmnumber=24
      elseif((i.ge.54.5).and.(i.le.55.1))then
        atmnumber=25
      elseif((i.ge.55.5).and.(i.le.56.1))then
        atmnumber=26
c     Co (27) and Ni (28) have very close mass ranges
      elseif((i.ge.58.76).and.(i.le.59.1))then
        atmnumber=27
      elseif((i.ge.58.5).and.(i.le.59.75))then
        atmnumber=28
      elseif((i.ge.62.9).and.(i.le.64.1))then
        atmnumber=29
      elseif((i.ge.64.9).and.(i.le.66.1))then
        atmnumber=30
      elseif((i.ge.69.5).and.(i.le.70.1))then
        atmnumber=31
      elseif((i.ge.72.5).and.(i.le.73.1))then
        atmnumber=32
      elseif((i.ge.74.5).and.(i.le.75.1))then
        atmnumber=33
      elseif((i.ge.78.5).and.(i.le.79.1))then
        atmnumber=34
      elseif((i.ge.79.5).and.(i.le.80.1))then
        atmnumber=35
      elseif((i.ge.83.5).and.(i.le.84.1))then
        atmnumber=36
      elseif((i.ge.84.9).and.(i.le.86.1))then
        atmnumber=37
      elseif((i.ge.87.5).and.(i.le.88.1))then
        atmnumber=38
      elseif((i.ge.88.5).and.(i.le.89.1))then
        atmnumber=39
      elseif((i.ge.90.9).and.(i.le.91.5))then
        atmnumber=40
      elseif((i.ge.92.5).and.(i.le.93.1))then
        atmnumber=41
      elseif((i.ge.95.5).and.(i.le.96.1))then
        atmnumber=42
      elseif((i.ge.97.9).and.(i.le.98.1))then
        atmnumber=43
      elseif((i.ge.109.9).and.(i.le.101.5))then
        atmnumber=44
      elseif((i.ge.102.5).and.(i.le.103.1))then
        atmnumber=45
      elseif((i.ge.105.9).and.(i.le.106.5))then
        atmnumber=46
      elseif((i.ge.107.5).and.(i.le.108.1))then
        atmnumber=47
      elseif((i.ge.111.9).and.(i.le.112.5))then
        atmnumber=48
      elseif((i.ge.114.5).and.(i.le.115.1))then
        atmnumber=49
      elseif((i.ge.118.5).and.(i.le.119.1))then
        atmnumber=50
      elseif((i.ge.121.5).and.(i.le.122.1))then
        atmnumber=51
      elseif((i.ge.127.5).and.(i.le.128.1))then
        atmnumber=52
      elseif((i.ge.126.5).and.(i.le.127.1))then
        atmnumber=53
      elseif((i.ge.130.9).and.(i.le.131.5))then
        atmnumber=54
      elseif((i.ge.132.5).and.(i.le.133.1))then
        atmnumber=55
      elseif((i.ge.136.9).and.(i.le.137.5))then
        atmnumber=56
      elseif((i.ge.138.5).and.(i.le.139.1))then
        atmnumber=57
      elseif((i.ge.139.9).and.(i.le.140.5))then
        atmnumber=58
      elseif((i.ge.140.6).and.(i.le.141.1))then
        atmnumber=59
      elseif((i.ge.144.0).and.(i.le.144.5))then
        atmnumber=60
      elseif((i.ge.144.9).and.(i.le.145.1))then
        atmnumber=61
      elseif((i.ge.150.0).and.(i.le.150.6))then
        atmnumber=62
      elseif((i.ge.151.5).and.(i.le.152.1))then
        atmnumber=63
      elseif((i.ge.156.9).and.(i.le.157.5))then
        atmnumber=64
      elseif((i.ge.158.5).and.(i.le.159.1))then
        atmnumber=65
      elseif((i.ge.162.0).and.(i.le.163.1))then
        atmnumber=66
      elseif((i.ge.164.5).and.(i.le.165.1))then
        atmnumber=67
      elseif((i.ge.166.5).and.(i.le.167.9))then
        atmnumber=68
      elseif((i.ge.168.0).and.(i.le.169.1))then
        atmnumber=69
      elseif((i.ge.172.9).and.(i.le.173.5))then
        atmnumber=70
      elseif((i.ge.174.0).and.(i.le.175.1))then
        atmnumber=71
      elseif((i.ge.178.0).and.(i.le.179.1))then
        atmnumber=72
      elseif((i.ge.180.0).and.(i.le.181.1))then
        atmnumber=73
      elseif((i.ge.183.0).and.(i.le.184.1))then
        atmnumber=74
      elseif((i.ge.185.9).and.(i.le.186.5))then
        atmnumber=75
      elseif((i.ge.189.9).and.(i.le.190.5))then
        atmnumber=76
      elseif((i.ge.191.9).and.(i.le.192.5))then
        atmnumber=77
      elseif((i.ge.194.9).and.(i.le.195.5))then
        atmnumber=78
      elseif((i.ge.196.5).and.(i.le.197.1))then
        atmnumber=79
      elseif((i.ge.200.0).and.(i.le.201.1))then
        atmnumber=80
      elseif((i.ge.203.9).and.(i.le.204.6))then
        atmnumber=81
      elseif((i.ge.206.9).and.(i.le.207.6))then
        atmnumber=82
      elseif((i.ge.208.5).and.(i.le.209.1))then
        atmnumber=83
c     Po atomic number 84 has the same mass range as Bi (83)
      elseif((i.ge.209.9).and.(i.le.210.1))then
        atmnumber=85
      elseif((i.ge.221.9).and.(i.le.222.1))then
        atmnumber=86
      elseif((i.ge.222.9).and.(i.le.223.1))then
        atmnumber=87
      elseif((i.ge.225.9).and.(i.le.226.1))then
        atmnumber=88
      elseif((i.ge.226.9).and.(i.le.227.1))then
        atmnumber=89
      elseif((i.ge.231.9).and.(i.le.232.1))then
        atmnumber=90
      elseif((i.ge.230.9).and.(i.le.231.1))then
        atmnumber=91
      elseif((i.ge.237.9).and.(i.le.238.1))then
        atmnumber=92
c     Np atomic number 93 has the same mass range as U (92)
      elseif((i.ge.243.9).and.(i.le.244.1))then
        atmnumber=94
      elseif((i.ge.242.9).and.(i.le.243.1))then
        atmnumber=95
      elseif((i.ge.246.9).and.(i.le.247.1))then
        atmnumber=96
c     Bk atomic number 97 has the same mass range as Cm (96)
      elseif((i.ge.250.9).and.(i.le.251.1))then
        atmnumber=98
      elseif((i.ge.251.9).and.(i.le.252.1))then
        atmnumber=99
      elseif((i.ge.256.9).and.(i.le.257.1))then
        atmnumber=100
      elseif((i.ge.257.9).and.(i.le.258.1))then
        atmnumber=101
      elseif((i.ge.258.9).and.(i.le.259.1))then
        atmnumber=102
      elseif((i.ge.261.9).and.(i.le.262.1))then
        atmnumber=103
      elseif((i.ge.266.9).and.(i.le.267.1))then
        atmnumber=104
      elseif((i.ge.267.9).and.(i.le.268.1))then
        atmnumber=105
      elseif((i.ge.270.9).and.(i.le.271.1))then
        atmnumber=106
      elseif((i.ge.269.9).and.(i.le.270.1))then
        atmnumber=107
      elseif((i.ge.268.9).and.(i.le.269.1))then
        atmnumber=108
      elseif((i.ge.277.9).and.(i.le.278.1))then
        atmnumber=109
      elseif((i.ge.280.9).and.(i.le.281.1))then
        atmnumber=110
c     Rg atomic number 112 has the same mass range as Ds (110)
      elseif((i.ge.284.9).and.(i.le.285.1))then
        atmnumber=112
      elseif((i.ge.285.9).and.(i.le.286.1))then
        atmnumber=113
      elseif((i.ge.288.9).and.(i.le.289.1))then
        atmnumber=114
c     Uup atomic number 115 has the same mass range as Fl (114)
      elseif((i.ge.292.9).and.(i.le.293.1))then
        atmnumber=116
      elseif((i.ge.293.9).and.(i.le.294.1))then
        atmnumber=117
c     Uuo atomic number 118 has the same mass range as Uus (117)
      endif
      return
      end function atmnumber

      real(8) function duni(idnode)
 
c*********************************************************************
c     
c     random number generator based on the universal
c     random number generator of marsaglia, zaman and tsang
c     (stats and prob. lett. 8 (1990) 35-39.) it must be
c     called once to initialise parameters u,c,cd,cm
c     
c*********************************************************************

      implicit none

      logical new
      integer i,j,k,l,m,ii,jj,ir,jr,idnode
      integer,dimension(8) :: values
      real(4) s,t,c,cd,cm,uni
      real(4), dimension(97) :: u

      data new/.true./

      save u,c,cd,cm,uni,ir,jr,new
CVAM
CVAM      call VTBEGIN(134, ierr)
CVAM

      if(new)then

c       initial values of i,j,k must be in range 1 to 178 (not all 1)
c       initial value of l must be in range 0 to 168.

        call date_and_time(values=values)
c       note, these date_and_time values can be 0
        i=mod(values(8)*3, 177) + 1
        j=mod(values(7)*23, 177) + 1
        k=mod(values(8)*5, 177) + 1
        l=mod(values(7)*3, 168)
c       DEBUG
c        i=25  
c        j=89 
c        k=100
c        l=150
c       END DEBUG

c       This is in case we find the same problem arises
c        write(nrite,"(/,'values for date and time', 9i6,/)")values
        if(idnode.eq.0)write(nrite,"(/,'i,j,k,l values for duni()', 4i5,
     &/)")i,j,k,l
        ir=97
        jr=33
        new=.false.

        do ii=1,97
          s=0.0
          t=0.5
          do jj=1,24
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32)s=s+t
            t=0.5*t
          enddo
          u(ii)=s
        enddo
        c =  362436.0/16777216.0
        cd= 7654321.0/16777216.0
        cm=16777213.0/16777216.0
        duni=0.d0
      else

c       calculate random number
        uni=u(ir)-u(jr)
        if(uni.lt.0.0)uni=uni+1.0
        u(ir)=uni
        ir=ir-1
        if(ir.eq.0)ir=97
        jr=jr-1
        if(jr.eq.0)jr=97
        c=c-cd
        if(c.lt.0.0)c=c+cm
        uni=uni-c
        if(uni.lt.0.0)uni=uni+1.0
        duni=dble(uni)
      endif
CVAM
CVAM     call VTEND(134, ierr)
CVAM
      end function duni

      integer function loc2(i,j)

c*********************************************************************
c
c     calculates double index array minimum reference
c
c*********************************************************************

      integer i,j

      loc2=(max(i,j)*(max(i,j)-1))/2+min(i,j)

      return
      end function loc2


      function errorq(w,delw,x,delx,y,dely,z,delz)
c***********************************************************************
c
c     calculates the standard error based on the theory that
c     Q=f(w,x,y,z) has an error of
c
c    dQ^2 = (df/dw)^2 * delw^2 + (df/dx)^2 * delx^2 + (df/dy)^2 * dely^2
c            + (df/dz)^2 * delz^2
c
c
c***********************************************************************
      implicit none
      real(8) w,delw,x,delx,y,dely,z,delz
      real(8) numerator, denominator
      real(8) errorq

      numerator=w-x*y
      denominator=z-y*y

c      denomsq=denominator*denominator

c      dfdw=1.d0/denominator
      
c      dfdx=-y/denominator

c      dfdy=2.d0*y*numerator/(denomsq)-
c     & x/denominator

c      dfdz=-1.d0*numerator/denomsq

c      errorq=sqrt((dfdw*delw/w)**2+(dfdx*delx/x)**2+
c     & (dfdy*dely/y)**2+(dfdz*delz/z)**2)

      errorq=abs(numerator/denominator)*
     &sqrt((delw/w)**2+(delx/x)**2+(dely/y)**2+(delz/z)**2)
      

      return

      end function errorq

      function calc_Qst(E, N, N2, EN, temp)
c***********************************************************************
c
c     computes the isosteric heat of adsorption 
c
c***********************************************************************
      implicit none
      real(8) calc_Qst, E, N, N2, EN, temp

      calc_Qst = -1.d0*(EN - E*N)/(N2-N*N) + Rgas*temp
      end function calc_Qst

      function old_cv(E2, E, N, N2, EN, temp)
c***********************************************************************
c
c     computes the heat capacity from gcmc 75 paper 
c
c***********************************************************************
      implicit none
      real(8) old_cv, E2, E, N, N2, EN, temp

      old_cv =(E2 - (E*E) - (EN - E*N)**2
     &   / (N2 - N*N)) / (N * kboltz * temp*temp) 
      end function old_cv

      function calc_Cv(E2, E, N, N2, EN, temp)
c***********************************************************************
c
c     computes the heat capacity at constant volume from Tildesley
c
c***********************************************************************
      implicit none
      real(8) calc_Cv, E2, E, N, N2, EN, temp

      calc_Cv = (3.d0/2.d0 * Rgas) + (E2 - (E*E) - (EN - E*N)**2
     &   / (N2 - N*N)) / (N * kboltz * temp*temp) 
      end function calc_Cv

      function errorcv(w, delw, x, delx, y, dely, z, delz, q, delq, 
     & temp)
c***********************************************************************
c
c     calculates the standard error based on the theory that
c     Q=f(w,x) has an error of
c
c    dQ^2 = (df/dw)^2 * delw^2 + (df/dx)^2 * delx^2 
c
c    w = avgE2, x = avgE, y = avgN, z = avgN2, q = avgEN
c***********************************************************************
      implicit none
      real(8) w, delw, x, delx, temp, numerator, denominator
      real(8) y, dely, z, delz, q, delq
      real(8) errorcv
      numerator=w-x*x - (q-x*y)**2/(z - y*y)
      denominator=kboltz*temp*temp*y
      errorcv=abs(numerator/denominator) *
     &sqrt((delw/w)**2+(delx/x)**2 + (dely/y)**2 + (delz/z)**2 +
     & (delq/q)**2) 
      return
      end function errorcv

      function sdot0(n,aaa,bbb)

c***********************************************************************
c     
c     utility to calculate scalar product of two arrays
c
c**********************************************************************

      implicit none

      integer n,i
      real(8) sdot0,aaa,bbb

      dimension aaa(*),bbb(*)

      sdot0=0.d0

      do i=1,n
        sdot0=sdot0+aaa(i)*bbb(i)
      enddo

      return
      end function sdot0

      integer function loc8(i,j)

c*********************************************************************
c
c     calculates array reference for a triangular matrix
c
c*********************************************************************
      
      integer i,j
      
      loc8=(max(i,j)*(max(i,j)-1))/2+min(i,j)
      
      return
      end function loc8

      character*3 function intstr3(nnn)

c*********************************************************************
c
c     converts a 3 digit integer to a string "001" etc.
c
c*********************************************************************

      implicit none

      integer nnn

      write(intstr3,'(i3.3)')nnn

      return
      end function intstr3

      subroutine abort_field_read(kode,idnode,nfield)

c***********************************************************************
c     
c     subroutine for aborting FIELD file read
c     
c***********************************************************************
      
      implicit none

      integer kode,idnode,nfield

      if(idnode.eq.0) close (nfield)

      if(kode.eq.1)then

c     end of field file error exit

        call error(idnode,52)

      else if(kode.eq.2)then

c     unrecognised directive in field file

        call error(idnode,4)

      endif

      return
      end subroutine abort_field_read

      subroutine abort_control_read(kode,idnode,nread)

c***********************************************************************
c     
c     subroutine for aborting CONTROL file read
c     
c***********************************************************************
      
      implicit none

      integer kode,idnode,nread

      if(idnode.eq.0) close (nread)

      if(kode.eq.1)then

c     end of control file error exit

        call error(idnode,53)

      else if(kode.eq.2)then

c     general error exit from field file processing

        call error(idnode,0)

      endif

      return
      end subroutine abort_control_read

      subroutine abort_config_read(kode,idnode,nconf)

c***********************************************************************
c     
c     subroutine for aborting CONTROL file read
c     
c***********************************************************************
      
      implicit none

      integer kode,idnode,nconf

      if(idnode.eq.0) close (nconf)

      if(kode.eq.1)then

c     general error exit from field file processing

        call error(idnode,54)

      else if(kode.eq.2)then

c     end of config file error exit

        call error(idnode,55)

      endif

      return
      end subroutine abort_config_read

      subroutine translate
     &(natms,mol,newx,newy,newz,cell,rcell,comx,comy,comz,
     &a,b,c)
c***********************************************************************
c                                                                      *
c     routine to translate a guest in the cell to a specified position *
c     in fractional coordinates.                                       *
c                                                                      *
c***********************************************************************
      implicit none
      integer natms,mol,i
      real(8), dimension(natms) :: newx,newy,newz
c     these are the unit cell vectors, random moves
c     will be done in the dimensions of the cell.
c     only need to be calculated once at the beginning
c     of the simulation.
      real(8), dimension(9) :: cell,rcell 
      real(8) comx,comy,comz
      real(8) a,b,c,x,y,z
      real(8) ssx,ssy,ssz,xss,yss,zss,shiftx,shifty,shiftz
      do i=1,natms
        newx(i)=a+newx(i)
        newy(i)=b+newy(i)
        newz(i)=c+newz(i)
      enddo
c     check pbc for com (keep the molecule intact)
      call com(natms,mol,newx,newy,newz,comx,comy,comz)

c     all this to make sure the COM lies within the boundary of 
c     the simulation.
      x=comx
      y=comy
      z=comz
 
      ssx=(rcell(1)*x+rcell(4)*y+rcell(7)*z)
      ssy=(rcell(2)*x+rcell(5)*y+rcell(8)*z)
      ssz=(rcell(3)*x+rcell(6)*y+rcell(9)*z)
          
      xss=ssx-floor(ssx)
      yss=ssy-floor(ssy)
      zss=ssz-floor(ssz)
       
      x=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
      y=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
      z=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
 
      shiftx=x-comx
      shifty=y-comy
      shiftz=z-comz
 
      comx=x
      comy=y
      comz=z
      newx=newx+shiftx
      newy=newy+shifty
      newz=newz+shiftz

      end subroutine translate

      subroutine random_ins(idnode,natms,iguest,rcut,delr)
c*****************************************************************************
c
c     routine to randomly insert a guest into the cell 
c     ignore insertions placed too close to framework atoms?
c     the radial scan might be a bit costly
c
c*****************************************************************************
      implicit none
      logical good,insert,done
      integer i,iguest,idnode
      integer natms
      real(8) randa,randb,randc,rmin,rcut,rclim,delr
      real(8) xc,yc,zc,comx,comy,comz,theta
      real(8) u1,u2,u3,q1,q2,q3,q4
      real(8) sintheta,beta1,beta2,beta3,norm,rand1,rand2,rand3,rand4

      randa=duni(idnode)
      randb=duni(idnode)
      randc=duni(idnode)
      rand1=duni(idnode)
      rand2=duni(idnode)
      rand3=duni(idnode)
      rand4=duni(idnode)
 
      good=.false.
      insert=.true.
      done=.false.
      xc=0.d0
      yc=0.d0
      zc=0.d0
      rmin=0.1d0
      rclim=(rcut+delr)**2

c     convert rand fractions to cartesian coordinates in the cell 
      call cartesian(randa,randb,randc,xc,yc,zc)
c     xc,yc,zc are the coordinates of the com guestx,guesty,guestz 
c     are the positions of atom i relative to the com
      do i=1, natms
        newx(i)=guestx(iguest,i)
        newy(i)=guesty(iguest,i)
        newz(i)=guestz(iguest,i) 
      enddo
c     the guestx,guesty, and guestz positions are centered at the
c     origin, this is done when the values are read from the 
c     FIELD file
      comx=0.d0
      comy=0.d0
      comz=0.d0 
c     rotate

c    setup quaternion which will uniformly sample the sphere
      
c      u1=duni()
c      u2=duni()
c      u3=duni()

c      u1sqrt=sqrt(u1)
c      u1m1sqrt=sqrt(1.d0-u1)
c     not sure if this representation is correct.
c     ie. the Q matrix is assuming the scalar term comes first
c     q(w,x,y,z) -- taken from wikipedia entry for "Rotation matrix"
c     the values for w,x,y,z were taken from 
c     "planning.cs.uiuc.edu/node198.html".. no mention where the 
c     scalar "w" lies in this representation... assuming it's the first
c     entry.
c      q1=u1m1sqrt*sin(2.d0*pi*u2)
c      q2=u1m1sqrt*cos(2.d0*pi*u2)
c      q3=u1sqrt*sin(2.d0*pi*u3)
c      q4=u1sqrt*cos(2.d0*pi*u3)

      theta=rand1*2.d0*pi
      beta1=2.d0*rand2-1.d0
      beta2=2.d0*rand3-1.d0
      beta3=2.d0*rand4-1.d0

c      beta1=duni()
c      beta2=duni()
c      beta3=duni()
      sintheta=sin(theta/2.d0)

      norm=sqrt(beta1*beta1+beta2*beta2+beta3*beta3)
      u1=beta1/norm
      u2=beta2/norm
      u3=beta3/norm

c     quaternion below
      q1=cos(theta/2.d0)
      q2=sintheta*u1
      q3=sintheta*u2
      q4=sintheta*u3
      call rotation(newx,newy,newz,comx,comy,comz,natms,q1,q2,q3,q4)
c      call rotationeuler(newx,newy,newz,natms,2.d0*pi) 
c     add the com to each atom in the guest

      do i=1, natms
        newx(i)=newx(i)+xc
        newy(i)=newy(i)+yc
        newz(i)=newz(i)+zc
      enddo


c     test to see if radial overlap with all other atoms. 
c      call images(imcon,natms,cell,newx,newy,newz)
c      really messed right now, commented out for later
c      investigation.  Uncomment the do while loop at the
c      beginning once this radial evaluation is working
c      call radial_eval(imcon,totatm,insert,natms,newx,newy,newz,good)
      
c      enddo
      return
      end subroutine random_ins

      subroutine random_disp
     &(idnode,delr,randa,randb,randc)
c*****************************************************************************
c
c     routine to randomly displace newx,newy,newz atoms
c
c*****************************************************************************
      implicit none
      integer idnode
      real(8) delr
      real(8) randa,randb,randc,halfdelr

      randa=duni(idnode)
      randb=duni(idnode)
      randc=duni(idnode)
      halfdelr=delr*0.5d0
c      tol=-3.d0/2.d0*delr+1.d-5
c     generate random numbers in the a,b, and c directions
c      do while(.not.done)
      randa=(2.d0*randa-1.d0)*halfdelr
      randb=(2.d0*randb-1.d0)*halfdelr
      randc=(2.d0*randc-1.d0)*halfdelr

      return
      end subroutine random_disp

      subroutine random_jump
     &(idnode,randa,randb,randc,rotangle)
c*****************************************************************************
c
c     place molecule at a completely random position in the cell
c     should be similar to a simultaneous insertion and deletion
c
c*****************************************************************************
      implicit none
      integer idnode
      real(8) randa,randb,randc,rotangle

      randa=duni(idnode) 
      randb=duni(idnode) 
      randc=duni(idnode) 

c     rotate randomly up to 2 pies 
      rotangle = 8.D0*datan(1.D0)

      return
      end subroutine random_jump

      subroutine rotationeuler(xxr,yyr,zzr,natms,angx,angy,angz)
c*********************************************************************      
c     subroutine to rotate by the old euler angles.  singularities
c     are a problem with this, implemented to test robustness of 
c     quaternion rotation
c*********************************************************************
      implicit none
      integer i,natms
      real(8) alpha,beta,kappa,xtmp,ytmp,ztmp
      real(8) angx,angy,angz
      real(8), dimension(9) :: R1,R2,R3
      real(8), dimension(natms) :: xxr,yyr,zzr
c     alpha beta kappa generated from rotangle
c      alpha=(rand1*2.d0-1.d0)*rotangle
c      beta=(rand2*2.d0-1.d0)*rotangle
c      kappa=(rand3*2.d0-1.d0)*rotangle
      alpha=angx*pi/180.d0
      beta=angy*pi/180.d0
      kappa=angz*pi/180.d0

c      write(*,*)cos(alpha)
c      write(*,*)cos(beta)
c      write(*,*)cos(kappa)
c      write(*,*)sin(alpha)
c      write(*,*)sin(beta)
c      write(*,*)sin(kappa)
c     R1, rotation about the x axis by angle alpha

      R1(1)=1.d0
      R1(2)=0.d0
      R1(3)=0.d0
      R1(4)=0.d0
      R1(5)=cos(alpha)
      R1(6)=-1*sin(alpha)
      R1(7)=0.d0
      R1(8)=sin(alpha)
      R1(9)=cos(alpha)

c     R2, rotation about the y axis by angle beta

      R2(1)=cos(beta)
      R2(2)=0.d0
      R2(3)=sin(beta)
      R2(4)=0.d0
      R2(5)=1.d0
      R2(6)=0.d0
      R2(7)=-1*sin(beta)
      R2(8)=0.d0
      R2(9)=cos(beta)

c     R3, rotation about the z axis by angle kappa

      R3(1)=cos(kappa)
      R3(2)=-1*sin(kappa)
      R3(3)=0.d0
      R3(4)=sin(kappa)
      R3(5)=cos(kappa)
      R3(6)=0.d0
      R3(7)=0.d0
      R3(8)=0.d0
      R3(9)=1.d0

c     do some matrix multiplication.  This may have a significant
c     impact on the speed of each gcmc step

c      if(rand4.lt.0.5)then
c        call matmul(R2,R3,mult2)
c      else
c        call matmul(R3,R2,mult2)
c      endif

      do i=1,natms
        xtmp=xxr(i)*R1(1)+yyr(i)*R1(4)+zzr(i)*R1(7)
        ytmp=xxr(i)*R1(2)+yyr(i)*R1(5)+zzr(i)*R1(8)
        ztmp=xxr(i)*R1(3)+yyr(i)*R1(6)+zzr(i)*R1(9)
        xxr(i)=xtmp
        yyr(i)=ytmp
        zzr(i)=ztmp
      enddo

      do i=1,natms
        xtmp=xxr(i)*R2(1)+yyr(i)*R2(4)+zzr(i)*R2(7)
        ytmp=xxr(i)*R2(2)+yyr(i)*R2(5)+zzr(i)*R2(8)
        ztmp=xxr(i)*R2(3)+yyr(i)*R2(6)+zzr(i)*R2(9)
        xxr(i)=xtmp
        yyr(i)=ytmp
        zzr(i)=ztmp
      enddo
      
      do i=1,natms
        xtmp=xxr(i)*R3(1)+yyr(i)*R3(4)+zzr(i)*R3(7)
        ytmp=xxr(i)*R3(2)+yyr(i)*R3(5)+zzr(i)*R3(8)
        ztmp=xxr(i)*R3(3)+yyr(i)*R3(6)+zzr(i)*R3(9)
        xxr(i)=xtmp
        yyr(i)=ytmp
        zzr(i)=ztmp
      enddo
      end subroutine rotationeuler
      subroutine rotation(xxr,yyr,zzr,comx,comy,comz,natms,q1,q2,q3,q4)
c*****************************************************************************
c
c     quaternion rotation
c     the quaternion is chosen randomly 
c
c*****************************************************************************
      implicit none
      integer i,natms
      real(8) q1,q2,q3,q4,comx,comy,comz
      real(8), dimension(natms) :: xxr,yyr,zzr,tmpx,tmpy,tmpz
      real(8), dimension(9) :: Q

      do i=1,natms
        xxr(i)=xxr(i)-comx
        yyr(i)=yyr(i)-comy
        zzr(i)=zzr(i)-comz
      enddo 
      do i=1,9
        Q(i)=0.d0
      enddo
c     taken from wikipedia entry "Rotation matrix"
      Q(1)=2.d0*q2*q2+2.d0*q1*q1-1.d0
      Q(2)=2.d0*q2*q3-2.d0*q1*q4
      Q(3)=2.d0*q2*q4+2.d0*q1*q3
      Q(4)=2.d0*q2*q3+2.d0*q1*q4
      Q(5)=2.d0*q1*q1+2.d0*q3*q3-1.d0
      Q(6)=2.d0*q3*q4-2.d0*q1*q2
      Q(7)=2.d0*q2*q4-2.d0*q1*q3
      Q(8)=2.d0*q3*q4+2.d0*q1*q2
      Q(9)=2.d0*q4*q4+2.d0*q1*q1-1.d0

      do i=1,natms
        tmpx(i)=Q(1)*xxr(i)+Q(2)*yyr(i)+Q(3)*zzr(i)
        tmpy(i)=Q(4)*xxr(i)+Q(5)*yyr(i)+Q(6)*zzr(i)
        tmpz(i)=Q(7)*xxr(i)+Q(8)*yyr(i)+Q(9)*zzr(i)
      enddo
      do i=1,natms
        xxr(i)=tmpx(i)+comx
        yyr(i)=tmpy(i)+comy
        zzr(i)=tmpz(i)+comz
      enddo

      return
      end subroutine rotation

      subroutine random_rot
     &(idnode,rotangle,q1,q2,q3,q4)
c*****************************************************************************
c
c     routine to generate a random rotation quaternion 
c
c*****************************************************************************
      implicit none
      integer idnode
      real(8) rotangle,theta,beta1,beta2,beta3,norm
      real(8) sintheta,q1,q2,q3,q4,u1,u2,u3
      real(8) rand1,rand2,rand3,rand4
c    setup quaternion which will uniformly sample the sphere
      rand1=duni(idnode)
      rand2=duni(idnode)
      rand3=duni(idnode)
      rand4=duni(idnode) 
        
      theta=rand1*rotangle
      beta1=2.d0*rand2-1.d0
      beta2=2.d0*rand3-1.d0
      beta3=2.d0*rand4-1.d0

      sintheta=sin(theta/2.d0)

      norm=sqrt(beta1*beta1+beta2*beta2+beta3*beta3)
      u1=beta1/norm
      u2=beta2/norm
      u3=beta3/norm

c     quaternion below
      q1=cos(theta/2.d0)
      q2=sintheta*u1
      q3=sintheta*u2
      q4=sintheta*u3

      return
      end subroutine random_rot

      subroutine com(natms,mol,newx,newy,newz,comx,comy,comz)
c*****************************************************************************
c
c     routine to calculate the centre of mass 
c     given an index and number of atoms to iterate over
c
c*****************************************************************************
      implicit none
      real(8) weight,comx,comy,comz,newx,newy,newz,mass
      integer natms,i,mol
      dimension newx(*),newy(*),newz(*)


c     initialized values
      comx=0.d0
      comy=0.d0
      comz=0.d0 
      weight=0.d0
c     calculate the centre of mass of the molecule
      do i=1,natms
        mass=atmwght(mol,i)
        weight=weight+mass
        comx=comx+newx(i)*mass
        comy=comy+newy(i)*mass
        comz=comz+newz(i)*mass
      enddo
      
      comx=comx/weight
      comy=comy/weight
      comz=comz/weight
      return
      end subroutine com


      subroutine round(val,dum,base)
c*********************************************************************
c     rounds a quad precision number to the nearest base value
c*********************************************************************
      implicit none
      real(8) val,dum,mult
      integer base

      mult=10.d0**base
      dum=dble(anint(val*mult))/mult

      end subroutine round
      subroutine overlap_check(loverlap, natms,overlap)
c*****************************************************************************
c
c     subroutine evaluates whether to accept or reject 
c     a gcmc move based on atomic distances
c     DEPRECATED! - PB 07/08/2017
c*****************************************************************************
      implicit none

      logical loverlap
      real(8) overlap
      integer natms,i,j
      loverlap = .false.

      if (overlap.gt.0.d0)then
        do i = 1, natms
          do j = 1, gstlentry(i)
              if (rsqdf(j).lt.overlap)then
                  loverlap=.true.
                  return
              endif
          enddo
        enddo
      endif
      return
      end subroutine overlap_check
      subroutine energy_eval
     &(eng,rande,volm,iguest,jguest,temp,beta,
     &displace,insert,delete,swap,accepted)
c*****************************************************************************
c
c     subroutine evaluates whether to accept or reject 
c     a gcmc move
c
c*****************************************************************************
      implicit none
      logical accepted,displace,insert,delete,swap
      integer jngsts,ingsts,imol,jmol,iguest,jguest
      real(8) rande,test,edummy
      real(8) volm,ipress,jpress,temp,beta
      real(8) eng
      accepted=.false.
      test=0.d0

      ipress = gstfuga(iguest)
      imol = locguest(iguest) 
      ingsts = dble(nummols(imol))
c     round the energy to the nearest second decimal place.
c     *** WARNING - the higher the value of eng, the worse this
c        rounding is ***
c      call round(duni(),rande,2)
c     distance check for the neighbours of each guest atom.
      edummy=eng
c     first two if statements prevent exp overflow
      if (-1*beta*edummy.gt.700.d0)then
        test=exp(700.d0)
      elseif(beta*edummy.gt.700.d0)then
        test=exp(-700.d0)
      elseif(displace)then
        test=exp(-1.d0*beta*edummy)
      elseif(insert)then
        test=volm*ipress/
     &    ((ingsts+1d0)*boltz*temp)*exp(-1.d0*beta*edummy)
      elseif(delete)then
        if(abs(ipress).lt.1.d-7)then
            test = 1.d0
        else
            test=ingsts*boltz*temp/(volm*ipress)*exp(-1.d0*beta*edummy)
        endif
      elseif(swap)then
        jpress=gstfuga(jguest)
        jmol = locguest(jguest)
        jngsts = real(nummols(jmol))
c       ingsts and jngsts have been prematurely decremented and
c       incremented respectively in gcmc.f
        test=jpress*(ingsts+1.d0)
     &/((jngsts)*ipress)*exp(-1.d0*beta*edummy)
      endif
      if(rande.lt.test)then
         accepted=.true.
      endif
      return
      end subroutine energy_eval
      subroutine cartesian(fraca,fracb,fracc,cartx,carty,cartz)
c*****************************************************************************
c
c     subroutine takes 3 fractional coordinates and returns 
c     3 cartesian coordinates based on the cell vectors
c
c*****************************************************************************
      implicit none
      real(8) fraca,fracb,fracc,cartx,carty,cartz

      cartx=fraca*cell(1)+fracb*cell(4)+fracc*cell(7)
      carty=fraca*cell(2)+fracb*cell(5)+fracc*cell(8)
      cartz=fraca*cell(3)+fracb*cell(6)+fracc*cell(9)

      return
      end subroutine cartesian
      
      subroutine guestlistgen
     &(imcon,totatm,rcut,delr,
     &natms,newx,newy,newz)
c*****************************************************************************
c    
c     builds temporary neighbour list.. area for 
c     improvement.. especially finding the index
c     of the xnew,ynew,and znew components.    
c 
c*****************************************************************************
      implicit none
      logical chk
      integer imcon,natms,i,ii,j,totatm
      integer itatms
      real(8) rsq,rclim,delr,rcut
      real(8),dimension(natms) :: newx,newy,newz
  
      rclim=(rcut+delr)**2
c     initialize the counter array
c      do i=1,natms
c      enddo
      do i=1,natms
        gstlentry(i) = 0
        ii=0
       
        do j=1,totatm
          chk=.true.
          do itatms=1,natms
            if (j.eq.ind(itatms))chk=.false.
          enddo
          
          if(chk)then
            ii=ii+1
            xdf(ii)=newx(i)-xxx(j)
            ydf(ii)=newy(i)-yyy(j)
            zdf(ii)=newz(i)-zzz(j)
          endif
        enddo
        call images(imcon,ii,cell,xdf,ydf,zdf)
        ii=0
        do j=1,totatm
          chk=.true.
          do itatms=1,natms
            if (j.eq.ind(itatms))chk=.false.
          enddo
         
          if(chk)then
            ii=ii+1
            if(imcon.eq.6)then
               rsq=xdf(ii)*xdf(ii)+ydf(ii)*ydf(ii)
            else 
               rsq=xdf(ii)*xdf(ii)+ydf(ii)*ydf(ii)+zdf(ii)*zdf(ii)
            endif
            if(rsq.le.rclim)then
              gstlentry(i)=gstlentry(i)+1
              gstlist(i,gstlentry(i))=j
            endif
          endif
        enddo
    
      enddo

      return
      end subroutine guestlistgen

      subroutine get_guest(iguest, choice, mol, natms, nmols)
c***********************************************************************
c
c     Grabs a guest from a specified index and populates the 
c     newx newy and newz arrays.
c
c***********************************************************************
      implicit none
      integer atmadd,iatm,iguest,choice,at,natms,nmols
      integer mol,imol,i

      mol=locguest(iguest)
      natms=numatoms(mol)
      nmols=nummols(mol)
      atmadd=0
      if(iguest.gt.1)then
        do i=1,iguest-1
          imol=locguest(i)
          atmadd=atmadd+numatoms(imol)*nummols(imol)
        enddo
      endif
      at=(choice-1)*natms+1
      iatm=0
      do i=at,at-1+natms
        iatm=iatm+1 
        ind(iatm)=atmadd+i
        newx(iatm)=molxxx(mol,i)
        newy(iatm)=molyyy(mol,i)
        newz(iatm)=molzzz(mol,i)
      enddo

      end subroutine get_guest

      subroutine guest_exclude(ntpguest)
c************************************************************************
c                                                                       * 
c     subroutine populates the exclude lists for each molecule          *
c                                                                       *
c************************************************************************
      implicit none
      logical lchk
      integer ntpguest,mol,jz,ii,jj
      integer natms,m,i,k
 
c     make a general list of exclusion sites 
c     then we can build the neighbour lists based on these general lists
c     initialize arrays
            do i=1,ntpguest
        mol=locguest(i)
        natms=numatoms(mol)
        do ii=1,natms
          nexsit(i,ii)=0
          do jj=1,natms
            lexsit(i,ii,jj)=0
          enddo
        enddo
      enddo
      do i=1,ntpguest
        mol=locguest(i)
        natms=numatoms(mol)

        do m=1,natms-1
          ii=0

          do k=m+1,natms
            lchk=.true.
            do jz=1,min(nexsit(i,m),mxexcl)
              if(lexsit(i,m,jz).eq.k)lchk=.false.
            enddo
            if(lchk)then
              nexsit(i,m)=nexsit(i,m)+1
              nexsit(i,k)=nexsit(i,k)+1
              lexsit(i,m,nexsit(i,m))=k
              lexsit(i,k,nexsit(i,k))=m
            endif
            
          enddo
        enddo
      enddo
c     once the sites are built, then build the atom list of excluded
c     atoms, this can be found in subroutine exclude_atom
c     in exclude_terms.f

c     this is the old code which populates a list for immediate
c     calculation

c      last=natms
c      mpm2=natms/2
c      npm2=(natms-1)/2
c      ii=0
c 
c      do i=1,natms
c        ii=ii+1
c        nexatm(ii)=0
c      enddo
c 
c 
c      do m=1,mpm2
c        if(m.gt.npm2)last=mpm2
c 
c        ii=0
c 
c        do i=1,last
c          j=i+m
c          if(j.gt.natms)j=j-natms
c 
c          ii=ii+1
c          nexatm(ii)=nexatm(ii)+1
c          lexatm(ii,nexatm(ii))=j+indx
c 
c        enddo
c      enddo
      return
      end subroutine guest_exclude

      subroutine radial_eval
     &(imcon,totatm,insert,ngstatm,newx,newy,newz,good)
c*****************************************************************************
c
c     subroutine evaluates the radial distance between  
c     the cartesian coordinates in xnew,ynew,and znew
c     if they are above the minimum tolerance then returns
c     good=.true.
c     also builds temporary neighbour list.. area for 
c     improvement..
c
c*****************************************************************************
      implicit none
      logical good,done,insert
      integer imcon,natms,i,j,ngstatm,jatm,totatm
      real(8) rsq,rmin
      real(8),dimension(ngstatm) :: newx,newy,newz
      good=.true.
      done=.false.
      rmin=2.25d0
c     the do while loop is so that if an atom overlap is found
c     the loop termninates prematurely
c     so the gcmc move can be re-started

      do while(.not.done)
        do i=1,ngstatm
          if(.not.insert)then
            natms=gstlentry(i)
          else
            natms=totatm
          endif
          do j=1,natms
            if(.not.insert)then
              jatm=gstlist(i,j)
            else
              jatm=j
            endif
            xdf(j)=newx(i)-xxx(jatm)
            ydf(j)=newy(i)-yyy(jatm)
            zdf(j)=newz(i)-zzz(jatm)
          enddo
          call images(imcon,natms,cell,xdf,ydf,zdf)
          do j=1,natms
            if(imcon.eq.6)then
               rsq=xdf(j)*xdf(j)+ydf(j)*ydf(j)
            else 
               rsq=xdf(j)*xdf(j)+ydf(j)*ydf(j)+zdf(j)*zdf(j)
            endif
            if(rsq.lt.rmin)then
              good=.false.
              done=.true.
            endif
            if(i.eq.ngstatm.and.j.eq.natms)done=.true.
          enddo
        enddo
      enddo
      return

      end subroutine radial_eval

      subroutine condense(totatm,ntpfram,ntpguest)
c*****************************************************************************
c     
c     subroutine takes the current list of atoms and condenses it
c     to a 1-d array.  This makes energy calculations and neighbour list
c     generations easier and quicker to access. 
c
c*****************************************************************************
      implicit none
      integer i,itmls,ntpguest,itgst,itmol,ntpfram
      integer mol,indnam,nmol,ii,isite,newatm
      integer j,k,kk,lsite,jj,indexsite,natms,totatm,latom
c     guests are placed first
      jj=0
      indexsite=0
      lsite=0 
      do itgst=1,ntpguest
        mol=locguest(itgst)
        natms=numatoms(mol)
        nmol=nummols(mol)
        ii=0
        
        do itmol=1,nmol 
          do isite=1,natms
            ii=ii+1
            jj=jj+1
            xxx(jj)=molxxx(mol,ii)
            yyy(jj)=molyyy(mol,ii)
            zzz(jj)=molzzz(mol,ii)
            atmcharge(jj)=atmchg(mol,isite)
            lfreezesite(jj)=lfzsite(mol,isite)
            atomname(jj)=atmname(mol,isite)
            atmweight(jj)=atmwght(mol,isite)
            ltype(jj)=ltpsit(mol,isite)
            moltype(jj) = mol
            kk=0
c      set up exclusion sites for the atoms
c      this is based on the assumption that 
c      guest atoms are RIGID and pairwise
c      interactions WILL NOT be calculated
c      between them.
            
            do k=1,nexsit(itgst,isite)
              newatm=lexsit(itgst,isite,k)+lsite

              if(((newatm.gt.jj).and.
     &          (newatm-jj.le.totatm/2)).or.
     &          ((newatm.lt.jj)
     &          .and.(newatm+totatm-jj.le.(totatm-1)/2)
     &          ))then
                kk=kk+1
                

                lexatm(jj,kk)=newatm

                if(kk.gt.1)then
                  do j=kk,2,-1
                    if(lexatm(jj,j).lt.lexatm(jj,j-1))
     &                then
                      latom=lexatm(jj,j)
                      lexatm(jj,j)=lexatm(jj,j-1)
                      lexatm(jj,j-1)=latom
                    endif
                  enddo
                endif
              endif
            enddo

            nexatm(jj)=kk

          enddo
          lsite=lsite+natms
        enddo
      enddo    
c     framework atoms are added after
      indnam=0
      do itmls=1,ntpfram
        ii=0
        mol=locfram(itmls)
        natms=numatoms(mol)
        nmol=nummols(mol)
        do i=1,nmol
          do k=1,natms
            jj=jj+1
            ii=ii+1
            xxx(jj)=molxxx(mol,ii)
            yyy(jj)=molyyy(mol,ii)
            zzz(jj)=molzzz(mol,ii)
            atmcharge(jj)=atmchg(mol,k)
            lfreezesite(jj)=lfzsite(mol,k)
            atomname(jj)=atmname(mol,k)
            atmweight(jj)=atmwght(mol,k)
            ltype(jj)=ltpsit(mol,k)
            moltype(jj) = mol
      
          enddo
        enddo
      enddo 
      return
      end subroutine condense

      subroutine angular_dist(newx,newy,newz,natms,maxanglegrid)
c*****************************************************************************
c     
c     calculate the angular distribution of all guests at a particular 
c     step
c
c*****************************************************************************
      implicit none
      integer natms,iangle,maxanglegrid
      real(8) angle,comx,comy,comz
      real(8), dimension(3) :: zaxis,vector1
      real(8), dimension(natms) :: newx,newy,newz

      zaxis=(/0.d0,0.d0,1.d0/)
      angle=0.d0 
  
      comx=0.d0
      comy=0.d0
      comz=0.d0

c      call com(natms,mol,newx,newy,newz,comx,comy,comz)

c     the first vector with which the angle is calculated 
c     is between the second atom of the molecule and the com

      vector1=(/newx(2)-newx(1),newy(2)-newy(1),newz(2)-newz(1)/)

      call calc_angle(angle,vector1,zaxis)

      iangle=nint(angle/(pi/dble(maxanglegrid)))
   
      if(iangle.eq.0)iangle=1
      angdist(iangle)=angdist(iangle)+1

      end subroutine angular_dist
      
      subroutine jobcheck(idnode,jobsafe,ljob,production)
c*****************************************************************************
c     
c     checks a file called jobcontrol.in for job information 
c
c*****************************************************************************
      implicit none
      logical jobsafe, ljob, production, check, safe
      integer idum,idnode

      check=.true.
      safe=.true.

      if(idnode.eq.0)open(205,file='jobcontrol.in',status='old')

      call getrec(safe,idnode,205)
c     below is just a non mpi version of 'getrec' from parse_module.f
c     ...why??
c      read(205,'(a150)',end=100)line

c      do i=1,lenrec
c        lrec(i)=line(i:i)
c      enddo 

      call lowcase(record,lenrec)
      if(findstring('terminate',record,idum))then
        jobsafe=.false.
      elseif(findstring('start averaging',record,idum))then
        if((ljob).and.(.not.production))then 
           production=.true.
           if(idnode.eq.0)write(nrite,'(/,3x,a,/)')
     &"'START AVERAGING' found in jobcontrol.in, starting production
     & averaging"
        endif
      endif 
c     this is just nothing
c100   check=.false.
      if(idnode.eq.0)close(205)
      return
      end subroutine jobcheck


      subroutine calc_angle(angle,vect1,vect2) 
c*****************************************************************************
c     
c     calculate the angle between two vectors of dimension 3
c
c*****************************************************************************
      implicit none
      real(8) angle
      real(8), dimension(3) :: vect1,vect2

c     this is the numerator of the angle calc (dot product of vectors)
      angle=acos(dot_product(vect1,vect2)/
     &sqrt(dot_product(vect1,vect1))/sqrt(dot_product(vect2,vect2)))

      end subroutine calc_angle
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

   10 format(/,'time elapsed since job start = ',f15.3,' seconds',/)

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
      end subroutine timchk

      end module utility_pack
