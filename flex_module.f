      module flex_module
      
c***********************************************************************
c     
c     gcmc module for framework interchanges
c     author    - t. daff    nov 2010
c     
c***********************************************************************

      use utility_pack
c      use readinputs

      implicit none

c     global parameters defined in the module:
c     nfwks     - Number of frameworks
c     fwk_pos   - Atom postitions in each framework
c     fwk_chg  - Charges on each atom in each framework
c     fwk_cell - Unit cells of each framework
c     fwk_ener  - Relative energy of each framework

      real(8), allocatable :: fwk_posx(:,:)
      real(8), allocatable :: fwk_posy(:,:)
      real(8), allocatable :: fwk_posz(:,:)
      real(8), allocatable :: fwk_chg(:,:)
      real(8), allocatable :: fwk_cell(:,:)
      real(8), allocatable :: fwk_vol(:)
      real(8), allocatable :: fwk_ener(:)
      character*1, allocatable :: fwk_name(:,:)
      integer, allocatable :: fwk_counts(:)
      integer n_fwk, new_fwk, fwk_step_max
      integer curr_fwk, fwk_initial
      logical l_fwk_seq

c Temporary stores when flexing
      real(8), allocatable :: state_x(:,:)
      real(8), allocatable :: state_y(:,:)
      real(8), allocatable :: state_z(:,:)
      real(8), allocatable :: state_chg(:,:)
      real(8), allocatable :: state_cell(:)
      real(8) state_vol
      

      save fwk_posx, fwk_posy, fwk_posz, fwk_chg, fwk_cell
      save fwk_vol, fwk_ener, n_fwk, curr_fwk, fwk_counts
      save state_x, state_y, state_z, state_chg, state_cell
      save state_vol, l_fwk_seq, fwk_step_max, fwk_initial

      contains

      subroutine flex_init(idnode, natms, imcon, ntpmls, maxmls,
     &totatm, rcut, celprp, ntpguest, volm)

c***********************************************************************
c     
c     Read all the structures
c     
c***********************************************************************

      implicit none

c intent :: in
      integer idnode, natms, ntpguest
      integer imcon, ntpmls, maxmls, totatm
      real(8) rcut, volm
      real(8), dimension(10) :: celprp
      integer fwkd(8)

c local
      integer fwknum, statemaxmls, statemaxalloc
      integer k, l, m, indatm, indfwk, gstidx
      logical isguest

c these should be the same size as the incoming data
      statemaxmls = size(nummols)
      statemaxalloc = size(xxx)

      allocate (fwk_posx(n_fwk,natms))
      allocate (fwk_posy(n_fwk,natms))
      allocate (fwk_posz(n_fwk,natms))
      allocate (fwk_chg(n_fwk,natms))
      allocate (fwk_cell(n_fwk,9))
      allocate (fwk_vol(n_fwk))
      allocate (fwk_ener(n_fwk))
      allocate (fwk_counts(n_fwk))
      allocate (fwk_name(n_fwk,80))

      allocate(state_x(statemaxmls,statemaxalloc))
      allocate(state_y(statemaxmls,statemaxalloc))
      allocate(state_z(statemaxmls,statemaxalloc))
      allocate(state_chg(statemaxmls,statemaxalloc))
      allocate(state_cell(9))

c read in all the frameworks
      do fwknum=1,n_fwk
        call readfwk (idnode,fwknum,imcon,ntpmls,maxmls,
     &totatm,rcut,celprp)
        call readfwkfld (idnode,fwknum)
        fwk_counts(fwknum) = 0
      enddo

c random not initialized yet use clock instead
c will only sample first 1000 frameworks...
      if(fwk_initial.gt.1)then
        curr_fwk = fwk_initial
        if(idnode.eq.0)
     &write(nrite,"(/,/,' Selected starting framework: ',i5)")curr_fwk
      else
        call date_and_time(values=fwkd)
        curr_fwk = mod(fwkd(8),n_fwk)+1
        if(idnode.eq.0)
     &write(nrite,"(/,/,' Random starting framework: ',i5)")curr_fwk
      endif
      fwk_counts(curr_fwk) = 1
c Replace starting config with 
      cell = fwk_cell(curr_fwk,:)
c put in new framework positions
      do k = 1,maxmls
        isguest = .false.
        do gstidx=1, ntpguest
          if (k.eq.locguest(gstidx)) isguest=.true.
        enddo
        indatm = 0
        indfwk = 0
        do l = 1, nummols(k)
          do m = 1, numatoms(k)
            indatm=indatm+1
            if(.not.isguest) then
              indfwk = indfwk+1
              molxxx(k,indatm) = fwk_posx(curr_fwk,indfwk)
              molyyy(k,indatm) = fwk_posy(curr_fwk,indfwk)
              molzzz(k,indatm) = fwk_posz(curr_fwk,indfwk)
            endif
          enddo
        enddo
      enddo
c switch charges
      do l = 1,ntpmls
        isguest = .false.
        do gstidx=1, ntpguest
          if (l.eq.locguest(gstidx)) isguest=.true.
        enddo
        if(.not.isguest) then
          do m = 1,numatoms(l) 
            atmchg(l,m) = fwk_chg(curr_fwk, m)
          enddo
        endif
      enddo
c set the correct volume
      volm = fwk_vol(curr_fwk)
      return
      end subroutine flex_init


      subroutine readfwk
     &(idnode,fwk_id,imcon,ntpmls,maxmls,
     &totatm,rcut,celprp)
c*************************************************************************
c
c     Read framework files; slightly modified from readconfig to 
c     load up framework arrays instead
c
c*************************************************************************
      implicit none

      character*1 atname(8),cfgname(80)
      logical safe
      integer icfg,imcon,idnode,indatm,maxmls
      integer k,l,m,totatm,i,ntpmls
      integer idum,levcfg,mxnode
      real(8) xcoord,ycoord,zcoord,junk,volm,axx,test,rt3
      real(8) rcut,width
      real(8), dimension(10) :: celprp

      character*8 fwkname
      integer fwk_id

      safe=.true.
      write(fwkname,'(a,i5.5)'),'FWK',fwk_id
      if(idnode.eq.0)then
        open(nconfig,file=fwkname,status='old')
      endif
c     read header info
      call getrec(safe,idnode,nconfig)
c HEADER IS ENERGY
      fwk_ener(fwk_id) = dblstr(record,lenrec,idum)
      call copystring(record,cfgname,60)
      fwk_name(fwk_id,:) = cfgname(:)
      if(idnode.eq.0)then
        write(nrite,'(/,1x,a8," name    : ",60a1)')
     &fwkname,(cfgname(i),i=1,60)
        write(nrite,'(1x,a8," energy  : ",g22.12)')
     &fwkname,fwk_ener(fwk_id)
      endif

      call getrec(safe,idnode,nconfig)

      levcfg=intstr(record,lenrec,idum) 

      if(imcon.eq.0)then
        if(idnode.eq.0)write(nrite,'(3x, "WARNING - no periodic 
     &boundaries requested!")')
        do i=1,9
          fwk_cell(fwk_id,i)=0.d0
        enddo
        volum=0.d0
      else
        call getrec(safe,idnode,nconfig)
        if(.not.safe)call abort_config_read(1,idnode,nconfig)
        fwk_cell(fwk_id,1)=dblstr(record,lenrec,idum)
        fwk_cell(fwk_id,2)=dblstr(record,lenrec,idum)
        fwk_cell(fwk_id,3)=dblstr(record,lenrec,idum)
        call getrec(safe,idnode,nconfig)
        if(.not.safe)call abort_config_read(1,idnode,nconfig)
        fwk_cell(fwk_id,4)=dblstr(record,lenrec,idum)
        fwk_cell(fwk_id,5)=dblstr(record,lenrec,idum)
        fwk_cell(fwk_id,6)=dblstr(record,lenrec,idum)
        call getrec(safe,idnode,nconfig)
        if(.not.safe)call abort_config_read(1,idnode,nconfig)
        fwk_cell(fwk_id,7)=dblstr(record,lenrec,idum)
        fwk_cell(fwk_id,8)=dblstr(record,lenrec,idum)
        fwk_cell(fwk_id,9)=dblstr(record,lenrec,idum)
      endif
c     read atomic coordinates
      indatm=0
      safe=.true.

      do k=1,maxmls
        indatm=0
        do l=1,nummols(k)
          do m=1,numatoms(k)

            indatm=indatm+1

            fwk_posx(fwk_id,indatm)=0.d0
            fwk_posy(fwk_id,indatm)=0.d0
            fwk_posz(fwk_id,indatm)=0.d0

c try running on all nodes...
c            if(idnode.eq.0)then
            if(levcfg.eq.0)then
              call getrec(safe,idnode,nconfig)
              call copystring(record,atname,8)
              call getrec(safe,idnode,nconfig)
              xcoord=dblstr(record,lenrec,idum)
              ycoord=dblstr(record,lenrec,idum)
              zcoord=dblstr(record,lenrec,idum)
 
            elseif(levcfg.eq.1)then
              call getrec(safe,idnode,nconfig)
              call copystring(record,atname,8)
              call getrec(safe,idnode,nconfig)
              xcoord=dblstr(record,lenrec,idum)
              ycoord=dblstr(record,lenrec,idum)
              zcoord=dblstr(record,lenrec,idum)
c ignore velocities
              call getrec(safe,idnode,nconfig)
            
            elseif(levcfg.eq.2)then
              call getrec(safe,idnode,nconfig)
              call copystring(record,atname,8)
              call getrec(safe,idnode,nconfig)
              xcoord=dblstr(record,lenrec,idum)
              ycoord=dblstr(record,lenrec,idum)
              zcoord=dblstr(record,lenrec,idum)
c ignore velocities and forces
              call getrec(safe,idnode,nconfig)
              call getrec(safe,idnode,nconfig)
            endif


            call strip(atname,8)

            if (atmname(k,m).eq.mkwd8(atname))then
                fwk_posx(fwk_id,indatm)=xcoord
                fwk_posy(fwk_id,indatm)=ycoord
                fwk_posz(fwk_id,indatm)=zcoord

            else
                write(nrite,"(/,/,'unidentified atom label :',8a1,
     &             ': atom number ',i5)")atname,indatm
                safe=.false.
            endif

c            endif

            call gstate(safe)
            if(.not.safe)call error(idnode,25)
           
          enddo
        enddo
      enddo

c     sum coordinate arrays across all nodes
      if(mxnode.gt.1)then
        call gd2dsum(fwk_posx(fwk_id,:),ntpmls,nummols,numatoms,totatm)
        call gd2dsum(fwk_posy(fwk_id,:),ntpmls,nummols,numatoms,totatm)
        call gd2dsum(fwk_posz(fwk_id,:),ntpmls,nummols,numatoms,totatm)
      endif

c     check integrity of cell vectors
      if((imcon.eq.1).or.(imcon.eq.4).or.(imcon.eq.5))then
         axx=(abs(fwk_cell(fwk_id,1))+abs(fwk_cell(fwk_id,5)))/2.d0
         test=1.d-8*axx
         if(abs(fwk_cell(fwk_id,1)-axx).gt.test)call error(idnode,410)
         if(abs(fwk_cell(fwk_id,5)-axx).gt.test)call error(idnode,410)
         if(imcon.eq.5)then
           if(abs(fwk_cell(fwk_id,9)-axx*sqrt(2.d0)).gt.test)
     &       call error(idnode,410)
         else
           if(abs(fwk_cell(fwk_id,9)-axx).gt.test)call error(idnode,410)
         endif
      endif

      if(imcon.eq.7)then
        rt3=sqrt(3.d0)
        if(abs(fwk_cell(fwk_id,1)-rt3*fwk_cell(fwk_id,5)).ge.1.d-6)
     &    call error(idnode,410)
      endif

      if(imcon.eq.6)then
        if(abs(fwk_cell(fwk_id,3)).gt.1.d-10) call error(idnode,410)
        if(abs(fwk_cell(fwk_id,6)).gt.1.d-10) call error(idnode,410)
        if(abs(fwk_cell(fwk_id,7)).gt.1.d-10) call error(idnode,410)
        if(abs(fwk_cell(fwk_id,8)).gt.1.d-10) call error(idnode,410)
      endif

      if((imcon.eq.1).or.(imcon.eq.2).or.(imcon.eq.4).or.
     & (imcon.eq.5).or.(imcon.eq.7))then

        if(abs(fwk_cell(fwk_id,2)).gt.1.d-10) call error(idnode,410)
        if(abs(fwk_cell(fwk_id,3)).gt.1.d-10) call error(idnode,410)
        if(abs(fwk_cell(fwk_id,4)).gt.1.d-10) call error(idnode,410)
        if(abs(fwk_cell(fwk_id,6)).gt.1.d-10) call error(idnode,410)
        if(abs(fwk_cell(fwk_id,7)).gt.1.d-10) call error(idnode,410)
        if(abs(fwk_cell(fwk_id,8)).gt.1.d-10) call error(idnode,410)

      endif

      call dcell(fwk_cell(fwk_id,:),celprp)
      if(imcon.eq.0)then

        fwk_vol(fwk_id)=0.d0

      elseif(imcon.eq.4)then

        fwk_vol(fwk_id)=0.5d0*celprp(10)

      elseif(imcon.eq.5)then

        fwk_vol(fwk_id)=0.5d0*celprp(10)

      elseif(imcon.eq.7)then

        fwk_vol(fwk_id)=0.5d0*celprp(10)

      else

        fwk_vol(fwk_id)=celprp(10)

      endif

      if(idnode.eq.0)then
c        write(nrite,"(/,/,1x,'simulation cell vectors'/,/)")
c        write(nrite,"(21x,3f12.6)")fwk_cell(fwk_id,:)
        write(nrite,
     &     "(1x,a8,' volume  : ',g22.12)")fwkname,fwk_vol(fwk_id)
      endif

c     check value of cutoff and reset if necessary

      if(imcon.gt.0)then
        
        width=min(celprp(7),celprp(8),celprp(9))/2.d0
        if(imcon.eq.4)width=sqrt(3.d0)*cell(1)/4.d0
        if(imcon.eq.5)width=fwk_cell(fwk_id,1)/2.d0
        if(imcon.eq.6)width=min(celprp(7),celprp(8))/2.d0

c     halt program if potential cutoff exceeds cell width
   
        if(rcut.gt.width) call error(idnode,95)
      endif

      if(idnode.eq.0)close(nconfig)

      return

100   call abort_config_read(2,idnode,nconfig)

      end subroutine readfwk

      subroutine readfwkfld
     &(idnode,fwk_id) 
c     &(idnode,ntpvdw,maxvdw,ntpatm,ntpmls,ntpguest,
c     &ntpfram,totatm,rvdw,dlrpot,engunit)
c*************************************************************************
c
c     Read charges for framework
c
c*************************************************************************
      implicit none
      
c      logical loop,loop2,safe,atmchk,blank
c      logical lunits,lguest
      integer itmols,msite,jsite,ntpatm,totatm,maxvdw
      integer ifrz,irept,nrept,ksite,isite,ntpmls
      integer lsite,natms
      integer i,n,j,k,idum,idnode,ntpvdw,ntpguest,ntpfram
      character*8 junk
      character*8 atom1
      character*1 message(80)
      real(8) weight,charge,engunit,rvdw,dlrpot,sumchg
c      real(8) comx,comy,comz,mass,gpress

      character*8 fldname
      logical loop,loop2,safe
      logical lguest
      integer fwk_id, local_mols
      integer local_nsite


      local_mols=0


      natms=0
      ntpguest=0
      ntpfram=0
      lguest=.false.
      safe=.true.
c      blank=.true.
c      lunits=.false.
c      engunit=1.d0
c      totatm=0

      loop=.true.

      write(fldname,'(a,i5.5)'),'FLD',fwk_id
      if(idnode.eq.0)then
        open(nfield,file=fldname,status='old')
      endif
 
c     allocate guest pressures
c      do i=1,ntpguest
c        gstpress(i)=-1.d0
c      enddo

      call getrec(safe,idnode,nfield)
      if(.not.safe)call abort_field_read(1,idnode,nfield)

      do while(loop)
        call getrec(safe,idnode,nfield)
        if(.not.safe)call abort_field_read(1,idnode,nfield)

c        convert to lowercase and strip leading blanks

        call lowcase(record, lenrec)
        call strip(record, lenrec) 

        if(record(1).eq.'#'.or.record(1).eq.' ')then

c        molecular specification

        elseif (findstring('molecu',record,idum))then
          local_mols=intstr(record,lenrec,idum)
           
c          if(idnode.eq.0)
c     &      write(nrite,"(/,/,1x,'number of molecular types',6x,i10)")
c     &      local_mols           

          do n=1,local_mols

c            if(idnode.eq.0)
c     &        write(nrite,"(/,1x,'molecular species type',9x,i10)")
c     &        n

c        get name of molecular species
             
            call getrec(safe,idnode,nfield)
            if(.not.safe)call abort_field_read(1,idnode,nfield)

c                     
            call lowcase(record,lenrec)
            if(findstring('&guest',record,idum))then
              call getword(junk,record,6,lenrec)
              lguest=.true.
c               ntpguest=ntpguest+1
c               locguest(ntpguest)=n
c             else
c               ntpfram=ntpfram+1
c               locfram(ntpfram)=n
            endif
c             call strip(record,lenrec)
c             call copystring(record,molnam(1,n),40)
            
c             if(idnode.eq.0)
c     &         write(nrite,"(/,1x,'name of species:',13x,40a1)")
c     &         (molnam(i,n),i=1,40)

            loop2=.true.
            do while(loop2)
              call getrec(safe,idnode,nfield)
              if(.not.safe)call abort_field_read(1,idnode,nfield)

              call lowcase(record, lenrec)
              call strip(record, lenrec)

              ksite=0

              if (findstring('nummol',record,idum))then
c                  nummls=intstr(record,lenrec,idum)
c                  nummols(n)=nummls
c                  if(idnode.eq.0)
c     &              write(nrite,"(1x,'number of molecules ',
c     &              10x,i10)") nummls
                  
              elseif (findstring('atoms',record,idum))then
c        read in atom name, site number, mass, charge,  etc
                local_nsite=intstr(record,lenrec,idum)
c                  numatoms(n)=numsit

c                  if(idnode.eq.0)then
c                    write(nrite,"(/,1x,'number of atoms   ',10x,i10)")
c     &                   numsit
c                    write(nrite,"(/,/,1x,'atomic characteristics:',
c     &              /,21x,' site',5x,'name',10x,'mass',8x,'charge',
c     &              4x,'repeat',4x,'freeze'/)")
c                  endif

                do isite=1,local_nsite
c                    if(ksite.lt.numatoms(n))then

                  call getrec(safe,idnode,nfield)
                  if(.not.safe)call 
     &              abort_field_read(1,idnode,nfield)
                  call copystring(record,message,80)
                  call getword(atom1,record,8,lenrec)
                  weight=dblstr(record,lenrec,idum)
                  charge=dblstr(record,lenrec,idum)
                    
                  if(lguest)then
c                       nrept=1
c                       guestx(ntpguest,isite)=dblstr(record,lenrec,idum)
c                       guesty(ntpguest,isite)=dblstr(record,lenrec,idum)
c                       guestz(ntpguest,isite)=dblstr(record,lenrec,idum)
c                       ifrz=0
                  else
                    nrept=intstr(record,lenrec,idum)
                    ifrz=intstr(record,lenrec,idum)
                    if(nrept.eq.0)nrept=1

                    do irept=1,nrept
c                         nsite iterates over all atom sites 
c                      nsite=nsite+1
 
                      ksite=ksite+1
                      fwk_chg(fwk_id,ksite) = charge
c                          atmname(n,ksite)=atom1
c                          atmwght(n,ksite)=weight
c                          atmchg(n,ksite)=charge
c                          lfzsite(n,ksite)=ifrz

                    enddo
c                      if(idnode.eq.0)then
c                        write(nrite,"(21x,i5,5x,a8,2f12.5,3i10)")
c     &                     nsite+1,atom1,weight,charge,nrept,ifrz
                  endif
                    
c          establish a list of unique atom types (for vdw calc)                    
c                      atmchk=.true.
                      
c                      do jsite=1,ntpatm
c                        if (atom1.eq.unqatm(jsite))then
c                          atmchk=.false.
c                          do irept=ksite,ksite-nrept+1,-1
c                            ltpsit(n,irept)=jsite
c                          enddo

c                        endif
c                      enddo
c                      if (atmchk)then
c                        ntpatm=ntpatm+1
c                        unqatm(ntpatm)=atom1
c                        do irept=ksite,ksite-nrept+1,-1
c                          ltpsit(n,irept)=ntpatm
c                        enddo
c                      endif 
c                    endif
                enddo

              elseif (findstring('finish',record,idum))then
                loop2=.false.

c                  natms=natms+nummols(n)*numatoms(n)
c                  if(natms.gt.maxatm+maxguest)call error(idnode,75)

c                  if(lguest)then
c                   calculate com of the guest and subtract
c                   so the reference atoms are centred around
c                   the com.
c                    mass=0.d0
c                    comx=0.d0
c                    comy=0.d0
c                    comz=0.d0
                lguest=.false.
c                    do i=1,numsit
c                      comx=comx+guestx(ntpguest,i)*atmwght(n,i)
c                      comy=comy+guesty(ntpguest,i)*atmwght(n,i)
c                      comz=comz+guestz(ntpguest,i)*atmwght(n,i)
c                      mass=mass+atmwght(n,i)
c                    enddo
c                    comx=comx/mass
c                    comy=comy/mass
c                    comz=comz/mass
c                    do i=1,numsit
c                      guestx(ntpguest,i)=guestx(ntpguest,i)-comx
c                      guesty(ntpguest,i)=guesty(ntpguest,i)-comy
c                      guestz(ntpguest,i)=guestz(ntpguest,i)-comz
c                    enddo
c                    if(idnode.eq.0)then
c                      write(nrite,"(/,1x,'initial guest positions',/)")
c                      do i=1,numsit
c                        write(nrite,"(3x,'site:',2x,a4,5x,3f12.6)")
c     &                   atmname(n,i),guestx(ntpguest,i),
c     &                   guesty(ntpguest,i),guestz(ntpguest,i)
c                      enddo
c                    endif
c                  endif
              endif
            enddo 
c             totatm=totatm+numsit*nummls 
          enddo
         
         elseif(findstring('close',record,idum))then
           loop=.false.

c         else
c           if(idnode.eq.0)write(nrite,'(100a)')record
c           call abort_field_read(2,idnode,nfield) 

         endif
      enddo
c     check system charge
c      do itmols=1,ntpmls
c        do lsite=1,numatoms(itmols)  
c          sumchg=sumchg+dble(nummols(itmols)*atmchg(itmols,lsite)) 
c        enddo
c      enddo

c      if(abs(sumchg).gt.1.0d-6)then
c         write(nrite,'(/,1x,a,f12.6,a)')
c     &  ' *** warning - total system charge:',sumchg,' ***'
c      endif

      if(idnode.eq.0)close(nfield)

      return
      end subroutine readfwkfld

      subroutine flex_move()

c***********************************************************************
c     
c     Evaluate a flexing move
c
c***********************************************************************

      implicit none


      return
      end subroutine flex_move

      end module flex_module
