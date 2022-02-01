      module mc_moves
      use readinputs
      use ewald_module
      use vdw_module
      use utility_pack
      implicit none

      contains

      subroutine surface_check
     &(iatm,j,mol,surftol,overlap,loverlap,latmsurf)
c***********************************************************************
c
c     Checks if the given atom 'iatm' is near a surface atom 'jatm'.
c     This assumes that a 'surface' atom is frozen, and that
c     the array populated with the squared distance between neighbour
c     atoms (rsqdf) is correctly indexed for 'jatm'.
c
c***********************************************************************
      implicit none
      logical latmsurf,loverlap
      integer aa,ab,ivdw,iatm,j,jatm,mol
      real(8) ak,sig,req,surftol,surftolsq,overlap
      latmsurf=.false.

      jatm=ilist(j)

      aa = ltpsit(mol,iatm)
      ab = ltype(jatm)
      if(aa.gt.ab)then
        ak=(aa*(aa-1.d0)*0.5d0+ab+0.5d0)
      else
        ak=(ab*(ab-1.d0)*0.5d0+aa+0.5d0)
      endif
      ivdw=lstvdw(int(ak))
      sig = prmvdw(ivdw,2)
      req = sig
c      req = sig*(2.d0**(1.d0/6.d0))
      surftolsq = (surftol+req)**2
      if ((rsqdf(j).lt.surftolsq).and.
     &(lfreezesite(jatm).eq.1))latmsurf=.true.
      if (rsqdf(j).lt.overlap)loverlap=.true.
      return
      end subroutine surface_check
      
      subroutine poreblock_check
     &(iatm,imcon,numblocks,loverlap)
c***********************************************************************
c
c     Checks if the given atom 'iatm' is within a poreblocking radius.
c     Assumes that we are checking a 'guest' molecule and that
c     the coordinates of that guest have already been populated
c     in newx, newy, and newz
c
c***********************************************************************
      implicit none
      logical loverlap
      integer iatm,numblocks,iblk,ii,imcon,rsq
      real(8) overlap

      ii=0
c     assume newx, newy, and newz are properly populated with the 
c     guest in question
      do iblk=1,numblocks
        ii=ii+1
        bxdf(ii)=newx(iatm)-bxxx(iblk)
        bydf(ii)=newy(iatm)-byyy(iblk)
        bzdf(ii)=newz(iatm)-bzzz(iblk)
      enddo

c     shift by periodic image
      call images(imcon,ii,cell,bxdf,bydf,bzdf)

      do iblk=1,numblocks
        rsq=xdf(iblk)**2+ydf(iblk)**2+zdf(iblk)**2
        if(rsq-bradii(iblk).le.0.d0)then
            loverlap=.true.
            exit
        endif
      enddo
      return
      end subroutine poreblock_check

      subroutine insert_guests
     &(idnode,imcon,totatm,ntpguest,ntpfram,iguest,nguests,rcut,ddelr,
     &sumchg,surftol,overlap,keyfce,alpha,drewd,volm,newld,kmax1,kmax2,
     &kmax3,epsq,dlrpot,ntpatm,maxvdw,engunit,delrc,maxmls,iter,lblock,
     &numblocks)
c*****************************************************************************
c
c     routine that inserts guests to a desired number. 
c     Overlap checks are done to ensure that the particles are in 
c     a reasonable position. The energies of each guest insertion
c     are computed, but not used when keeping or discarding a particle.
c 
c     A maxiter counter is set to ensure that the loop doesn't continue
c     forever. If one reaches this counter, the program will exit with
c     an error. (currently set to 1,000 times the number of guests to
c     insert)
c
c     NB: none of this is coupled to an ensemble. 
c
c*****************************************************************************
      implicit none
      logical loverlap,lnewsurf,lblock
      integer idnode,imcon,ntpguest,iguest,nguests,keyfce
      integer maxfactor,maxiter,mol,nmol,iter,natms,totatm
      integer kmax1,kmax2,kmax3,ntpatm,maxvdw,maxmls,newld
      integer idum,ntpfram,ii,numblocks
      real(8) rcut,ddelr,surftol,estep,sumchg,overlap,alpha
      real(8) drewd,volm,epsq,dlrpot,engunit,delrc,chgtmp
      real(8) engsictmp
      maxfactor=1000
      maxiter=maxfactor*nguests
      mol=locguest(iguest)
      nmol=nummols(mol)
      natms=numatoms(mol)
      iter=0
      ! just return if the user doesn't know what they are doing.
      if(nmol.ge.nguests)return

      do while(nmol.lt.nguests)
        ! exit if the number of attempts exceeds the maximum number
        ! of allowed iterations.
        if(iter.ge.maxiter)call error(idnode, 2320)
        call random_ins(idnode,natms,iguest,rcut,ddelr)
        ! add guest
        estep=0.d0
        engsicorig = engsic
        call insertion
     &(imcon,iguest,keyfce,alpha,rcut,ddelr,drewd,totatm,
     &volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,delrc,estep,sumchg,chgtmp,engsictmp,maxmls,loverlap,
     &lnewsurf,surftol,overlap,newld,lblock,numblocks)
        ! if overlap, do nothing
        if(.not.loverlap)then
          call accept_move
     &(iguest,.true.,.false.,.false.,
     &lnewsurf,delrc,totatm,idum,ntpfram,ntpguest,
     &maxmls,sumchg,engsictmp,chgtmp,newld)
        else
          engsic = engsicorig
        endif
        delE(:) = 0.d0
        ! update the molecule count
        nmol=nummols(mol)
        iter=iter+1
      enddo 
      end subroutine insert_guests
      subroutine insertion
     &(imcon,iguest,keyfce,alpha,rcut,ddelr,drewd,totatm,
     &volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,delrc,estep,sumchg,chgtmp,engsictmp,maxmls,
     &loverlap,lnewsurf,surftol,overlap,newld,lblock,numblocks)
c***********************************************************************
c
c     inserts a particle in the framework and computes the 
c     energy of an additional particle.
c
c***********************************************************************
      implicit none
      logical loverlap,lnewsurf,latmsurf,lblock
      integer i,ik,j,kmax1,kmax2,kmax3,imcon,keyfce,numblocks
      integer totatm,sittyp,ntpatm,maxvdw,newld
      integer mol,natms,nmols,iguest
      integer jatm,l,mxcmls,maxmls
      real(8) drewd,dlrpot,volm,epsq,alpha,rcut,ddelr
      real(8) chg,engunit,sumchg,engsictmp
      real(8) ewld2sum,ewld3sum,vdwsum,chgtmp
      real(8) ewld1eng,ewld2eng,vdweng,delrc_tmp
      real(8) delrc,estep,surftol
      real(8) overlap
c     calculate the ewald sums ewald1 and ewald2(only for new mol)
c     store ewald1 and 2 separately for the next step.
c     calculate vdw interaction (only for new mol)
      mol=locguest(iguest)
      mxcmls=maxmls*(maxmls-1)/2 + maxmls
      lnewsurf=.false. 
      latmsurf=.false.
      loverlap = .false.
      natms=numatoms(mol)
      nmols=nummols(mol)
      ewld2sum=0.d0
      vdwsum=0.d0
      ind = 0 
      do ik=1,mxcmls
        delE(ik)=0.d0
        ewald1en(ik)=0.d0
        ewald2en(ik)=0.d0
        vdwen(ik)=0.d0
        delrc_mol(ik)=0.d0
      enddo
      call guestlistgen
     &(imcon,totatm,rcut,ddelr,
     &natms,newx,newy,newz)

c     calculate long range correction to vdw for the insertion
c     of an additional guest
      call gstlrcorrect(imcon,iguest,keyfce,ntpatm,maxvdw,
     &delrc,volm,maxmls,.false.)
c     do vdw and ewald2 energy calculations for the new atoms
      do i=1,natms
        ik=0
        do j=1,gstlentry(i)
          ik=ik+1
          jatm=gstlist(i,j)
          ilist(j)=jatm
          moldf(j)=moltype(jatm)
          xdf(ik)=newx(i)-xxx(jatm)
          ydf(ik)=newy(i)-yyy(jatm)
          zdf(ik)=newz(i)-zzz(jatm)
        enddo
        call images(imcon,ik,cell,xdf,ydf,zdf)

        do l=1,gstlentry(i)
          rsqdf(l)=xdf(l)**2+ydf(l)**2+zdf(l)**2
          call surface_check
     &(i,l,mol,surftol,overlap,loverlap,latmsurf)

          if(latmsurf)lnewsurf=.true.
          if(loverlap)return
        enddo
        if(lblock)call poreblock_check
     &(i,imcon,numblocks,loverlap)
        if(loverlap)return
c       figure out which index contains charges and ltype arrays
c       that match the guest...
        chg=atmchg(mol,i)
        sittyp=ltpsit(mol,i)
c       calc ewald2 interactions
        call ewald2
     & (chg,gstlentry(i),ewld2eng,mol,
     &  drewd,rcut,epsq)
        ewld2sum=ewld2sum+ewld2eng
c       calc vdw interactions 

        call srfrce
     & (sittyp,gstlentry(i),mol,vdweng,rcut,dlrpot)
        vdwsum=vdwsum+vdweng 
      enddo
c     the pairwise intramolecular coulombic correction
c     calculated for the guest at the begining of the 
c     program run.  Assumes constant bond distances.
      engsictmp=0.d0
      chgtmp=0.d0
      call ewald1_guest
     &(imcon,ewld1eng,natms,iguest,volm,alpha,sumchg,
     &chgtmp,engsictmp,kmax1,kmax2,kmax3,epsq,maxmls,newld,
     &.false.)
      ewld3sum=ewald3en(mol)
      estep= estep+ 
     &       (ewld1eng+ewld2sum-ewld3sum+vdwsum+delrc)/engunit
      
c      print *,"EWALD1 ",ewld1eng/engunit,-ewld3sum/engunit
c      print *,"EWALD2 ",ewld2sum/engunit
c      print *,"VDW    ",vdwsum/engunit
c      print *,"DELRC  ",delrc/engunit

      do i=1, maxmls
        ewld3sum=0.d0
        ik=loc2(mol,i)
        delrc_tmp=delrc_mol(ik)

        if(i.ne.mol)delE(i)=delE(i)+
     &(ewald1en(ik)+ewald2en(ik)+vdwen(ik)+delrc_tmp)
     &/engunit
c          write(*,*)nummols(i),
c     &ewald1en(i)/engunit,
c     &ewald2en(i)/engunit,
c     &vdwen(i)/engunit,
c     &delrc_tmp
c        endif
      enddo
      delE(mol)=delE(mol)+estep
      return
      end subroutine insertion
      subroutine guest_energy
     &(imcon,iguest,alpha,rcut,ddelr,drewd,
     &totatm,maxmls,volm,kmax1,kmax2,kmax3,epsq,dlrpot,
     &engunit,vdwsum,ewld2sum,ewld1eng,lsurf,newld,
     &surftol,sumchg,chgtmp,engsictmp,loverlap,overlap,estep,lexisting,
     &lblock,numblocks)
c***********************************************************************
c
c     compute the energy of a guest in it's current position.
c     Interaction energies between the guest and itself are avoided
c     by populating the ind() array with the indices of the guest.
c
c***********************************************************************
      implicit none
      logical lsurf,loverlap,latmsurf,lexisting,lblock
      integer i,ik,j,kmax1,kmax2,kmax3,imcon
      integer totatm,sittyp
      integer mol,natms,nmols,iguest,numblocks
      integer jatm,l,itatm
      integer maxmls,mxcmls,newld
      real(8) drewd,dlrpot,volm,epsq,alpha,rcut,ddelr
      real(8) chg,engunit,chgtmp,engsictmp
      real(8) ewld2sum,ewld3sum,vdwsum,overlap
      real(8) ewld1eng,ewld2eng,vdweng,sumchg
      real(8) estep,surftol
      mol=locguest(iguest)
      natms=numatoms(mol)
      nmols=nummols(mol)
      ewld2sum=0.d0
      vdwsum=0.d0
      lsurf=.false.
      loverlap=.false.
      mxcmls=maxmls*(maxmls-1)/2 + maxmls
      do ik=1,mxcmls
        delE(ik)=0.d0
        ewald1en(ik)=0.d0
        ewald2en(ik)=0.d0
        vdwen(ik)=0.d0
        delrc_mol(ik)=0.d0
      enddo
      call guestlistgen
     &(imcon,totatm,rcut,ddelr,
     &natms,newx,newy,newz)
      estep=0.d0
      chgtmp=0.d0
      engsictmp=0.d0

c     get the long range contributions of this guest
      call ewald1_guest
     &(imcon,ewld1eng,natms,iguest,volm,alpha,sumchg,
     &chgtmp,engsictmp,kmax1,kmax2,kmax3,epsq,maxmls,newld,
     &lexisting)
      
c     do vdw and ewald2 energy calculations for the atoms
      do i=1,natms
        itatm=ind(i)
        ik=0 
        do j=1,gstlentry(i)
          ik=ik+1
          jatm=gstlist(i,j)
          ilist(j)=jatm
          moldf(j)=moltype(jatm)
          xdf(ik)=newx(i)-xxx(jatm)
          ydf(ik)=newy(i)-yyy(jatm)
          zdf(ik)=newz(i)-zzz(jatm)
        enddo
        
        call images(imcon,ik,cell,xdf,ydf,zdf)
        do l=1,gstlentry(i)
          rsqdf(l)=xdf(l)**2+ydf(l)**2+zdf(l)**2
c         check if a surface atom
          call surface_check
     &(i,l,mol,surftol,overlap,loverlap,latmsurf)
          if(latmsurf)lsurf=.true.
          if(loverlap)return
        enddo
        if((.not.lexisting).and.lblock)call poreblock_check
     &(i,imcon,numblocks,loverlap)
        if(loverlap)return
        chg=atmcharge(itatm)
        sittyp=ltype(itatm)
        call ewald2
     & (chg,gstlentry(i),ewld2eng,mol,
     &  drewd,rcut,epsq)
        ewld2sum=ewld2sum+ewld2eng
c       calc vdw interactions
        call srfrce
     & (sittyp,gstlentry(i),mol,vdweng,rcut,dlrpot)
        vdwsum=vdwsum+vdweng
      enddo
      ewld3sum=ewald3en(mol)
      if(lexisting)then
        estep = (-ewld1eng + ewld2sum + vdwsum - ewld3sum)/engunit
      else
        estep = (ewld1eng + ewld2sum + vdwsum - ewld3sum)/engunit
      endif
      return
      end subroutine guest_energy
      subroutine deletion 
     &(imcon,keyfce,iguest,choice,alpha,rcut,ddelr,drewd,maxmls,
     &totatm,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,delrc,estep,linitsurf,surftol,sumchg,engsictmp,chgtmp,
     &overlap,newld)
c***********************************************************************
c
c     deletes a particle from the framework 
c
c***********************************************************************
      implicit none
      logical linitsurf,loverlap,latmsurf
      integer i,ik,j,kmax1,kmax2,kmax3,imcon,keyfce
      integer totatm,sittyp,ntpatm,maxvdw
      integer mol,natms,nmols,iguest
      integer jatm,l,itatm,choice
      integer mxcmls,maxmls,step,newld
      real(8) drewd,dlrpot,volm,epsq,alpha,rcut,ddelr
      real(8) chg,engunit,engsictmp,chgtmp
      real(8) ewld2sum,ewld3sum,vdwsum,sumchg,overlap
      real(8) ewld1eng,ewld2eng,vdweng,delrc_tmp
      real(8) delrc,estep,surftol
      
      ewld2sum=0.d0
      vdwsum=0.d0
      step=0
c     find which index the molecule "choice" is
      call get_guest(iguest,choice,mol,natms,nmols)
      linitsurf=.false.
      loverlap=.false.
      mxcmls=maxmls*(maxmls-1)/2 + maxmls
      do ik=1,mxcmls
        delE(ik)=0.d0
        ewald1en(ik)=0.d0
        ewald2en(ik)=0.d0
        vdwen(ik)=0.d0
        delrc_mol(ik)=0.d0
      enddo
      call gstlrcorrect(imcon,iguest,keyfce,ntpatm,maxvdw,
     &delrc,volm,maxmls,.true.)

      call guestlistgen
     &(imcon,totatm,rcut,ddelr,
     &natms,newx,newy,newz)
      estep=0.d0
      chgtmp=0.d0
      engsictmp=0.d0

c     get the long range contributions of this guest
      call ewald1_guest
     &(imcon,ewld1eng,natms,iguest,volm,alpha,sumchg,
     &chgtmp,engsictmp,kmax1,kmax2,kmax3,epsq,maxmls,newld,
     &.true.)
      
c     do vdw and ewald2 energy calculations for the atoms
      do i=1,natms
        itatm=ind(i)
        ik=0 
        do j=1,gstlentry(i)
          ik=ik+1
          jatm=gstlist(i,j)
          ilist(j)=jatm
          moldf(j)=moltype(jatm)
          xdf(ik)=newx(i)-xxx(jatm)
          ydf(ik)=newy(i)-yyy(jatm)
          zdf(ik)=newz(i)-zzz(jatm)
        enddo
        
        call images(imcon,ik,cell,xdf,ydf,zdf)
        do l=1,gstlentry(i)
          rsqdf(l)=xdf(l)**2+ydf(l)**2+zdf(l)**2
c         check if a surface atom
          call surface_check
     &(i,l,mol,surftol,overlap,loverlap,latmsurf)
          if(latmsurf)linitsurf=.true.
        enddo
        chg=atmcharge(itatm)
        sittyp=ltype(itatm)
        call ewald2
     & (chg,gstlentry(i),ewld2eng,mol,
     &  drewd,rcut,epsq)
        ewld2sum=ewld2sum+ewld2eng
c       calc vdw interactions
        call srfrce
     & (sittyp,gstlentry(i),mol,vdweng,rcut,dlrpot)
        vdwsum=vdwsum+vdweng
      enddo
      ewld3sum=ewald3en(mol)
      estep = (-ewld1eng + ewld2sum + vdwsum - ewld3sum)/engunit

c     calculate the pairwise intramolecular coulombic correction
c     (calculated at the begining - assumes constant bond distance)
c      estep= estep+ 
c     &     (ewld1eng-ewld2sum+ewld3sum-vdwsum-delrc)/engunit
      estep=estep-delrc/engunit
      do i=1, maxmls
        ik=loc2(mol,i)
        delrc_tmp=delrc_mol(ik)
        if(i.ne.mol)delE(i)=delE(i)-
     &(ewald1en(ik)+ewald2en(ik)+vdwen(ik)+delrc_tmp)
     &/engunit
      enddo
      delE(mol)=delE(mol)-
     &  estep
c      print *, delrc_mol(2), delrc_mol(3), delrc_mol(5),delrc_mol(8)
c      print *, ewald1en(2), ewald1en(3), ewald1en(5), ewald1en(8)
c      print *, ewald2en(2), ewald2en(3), ewald2en(5), ewald2en(8)
c      print *, vdwen(2), vdwen(3), vdwen(5), vdwen(8)
c      delE=-1.d0*delE
      end subroutine deletion
      subroutine accept_move
     &(iguest,insert,delete,displace,
     &lsurf,delrc,totatm,choice,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
c***********************************************************************
c
c     updates arrays for the rejection according to the type
c     of move.
c
c***********************************************************************
      implicit none
      integer iguest,natms,mol,i,mm,totatm,choice,at,k,newld
      integer ntpfram,ntpguest,maxmls,kk,ka,iatm
      logical insert,delete,displace,lsurf
      real(8) delrc,sumchg,chgtmp,engsictmp
c      ewaldaverage = ewaldaverage + abs(ewld1eng+ewld3sum)/engunit
      mol = locguest(iguest)
      kk = loc2(mol,mol)
      natms = numatoms(mol)
      at=(choice-1)*natms+1
c     update energy arrays
      do i=1,maxmls
         energy(i)=energy(i)+delE(i)
      enddo
c     update ewald1 sums
      if(insert)then
        do k=1,newld
          ckcsum(mol,k) = ckcsum(mol,k)+ckcsnew(mol,k)
          ckssum(mol,k) = ckssum(mol,k)+ckssnew(mol,k)
          ckcsum(maxmls+1,k) = ckcsum(maxmls+1,k)+ckcsnew(maxmls+1,k)
          ckssum(maxmls+1,k) = ckssum(maxmls+1,k)+ckssnew(maxmls+1,k)
c          ckcsnew(mol,k)=0.d0
c          ckssnew(mol,k)=0.d0
c          ckcsnew(maxmls+1,k)=0.d0
c          ckssnew(maxmls+1,k)=0.d0
        enddo
        mm=natms*nummols(mol)
        engsic(kk)=engsic(kk)+engsictmp
        chgsum_mol(mol)=chgsum_mol(mol)+chgtmp
        sumchg = sumchg+chgtmp
c       tally surface molecules
        if(lsurf)surfacemols(mol) = surfacemols(mol) + 1
c       update atomic coordinates
c       & increment type counters
        do i=1,natms
          molxxx(mol,mm+i)=newx(i)
          molyyy(mol,mm+i)=newy(i)
          molzzz(mol,mm+i)=newz(i)
          ka=ltpsit(mol,i)
          numtyp(ka)=numtyp(ka)+1
          numtyp_mol(mol,ka)=numtyp_mol(mol,ka)+1
          if(lfzsite(mol,i).ne.0)then
            numfrz(ka)=numfrz(ka)+1
            numfrz_mol(mol,ka)=numfrz_mol(mol,ka)+1
          endif
        enddo
c       update long range correction
        elrc=elrc+delrc
        elrc_mol(:)=elrc_mol(:)+delrc_mol(:)
c       update nummols,totatm, then condense everything to 1d arrays
        nummols(mol)=nummols(mol)+1
        totatm=totatm+natms
        call condense(totatm,ntpfram,ntpguest)
c       update the choice variable in case the user does something
c       with this after
        choice = nummols(mol)
      elseif(delete)then
        do k=1,newld
          ckcsum(mol,k) = ckcsum(mol,k)-ckcsnew(mol,k)
          ckssum(mol,k) = ckssum(mol,k)-ckssnew(mol,k)
          ckcsum(maxmls+1,k) = ckcsum(maxmls+1,k)-ckcsnew(maxmls+1,k)
          ckssum(maxmls+1,k) = ckssum(maxmls+1,k)-ckssnew(maxmls+1,k)
c          ckcsnew(mol,k)=0.d0
c          ckssnew(mol,k)=0.d0
c          ckcsnew(maxmls+1,k)=0.d0
c          ckssnew(maxmls+1,k)=0.d0
        enddo
        engsic(kk)=engsic(kk)-engsictmp
        chgsum_mol(mol)=chgsum_mol(mol)-chgtmp
        sumchg = sumchg-chgtmp
        mm=natms*nummols(mol)
c        mm=choice*natms
c       update surface molecules
        if(lsurf)surfacemols(mol) = surfacemols(mol) - 1
c       update atomic coordinates
c       & decrement counters
        iatm=0
c       This loops over all other molecules after the
c       deleted molecule to shift them back.
        do i=at,mm
          iatm=iatm+1
          molxxx(mol,i)=molxxx(mol,i+natms)
          molyyy(mol,i)=molyyy(mol,i+natms)
          molzzz(mol,i)=molzzz(mol,i+natms)
c         make sure we are only decrementing the 
c         site counts for the single molecule deleted.
          if(iatm.le.natms)then
            ka=ltpsit(mol,iatm)
            numtyp(ka)=numtyp(ka)-1
            numtyp_mol(mol,ka)=numtyp_mol(mol,ka)-1
            if(lfzsite(mol,iatm).ne.0)then
              numfrz(ka)=numfrz(ka)-1
              numfrz_mol(mol,ka)=numfrz_mol(mol,ka)-1
            endif
          endif
        enddo

c       update nummols,totatm, then condense everything to 1d arrays
        elrc=elrc-delrc
        elrc_mol(:)=elrc_mol(:)-delrc_mol(:)
        nummols(mol)=nummols(mol)-1
        totatm=totatm-natms
        call condense(totatm,ntpfram,ntpguest)
c        call images(imcon,totatm,cell,xxx,yyy,zzz)
      elseif(displace)then
c       this sums over all kpoints, should keep as one loop
        do k=1,newld
          ckcsum(mol,k) = ckcsum(mol,k)+ckcsnew(mol,k)
          ckssum(mol,k) = ckssum(mol,k)+ckssnew(mol,k)
          ckcsum(maxmls+1,k) = ckcsum(maxmls+1,k)+ckcsnew(maxmls+1,k)
          ckssum(maxmls+1,k) = ckssum(maxmls+1,k)+ckssnew(maxmls+1,k)
c          ckcsnew(mol,k)=0.d0
c          ckssnew(mol,k)=0.d0
c          ckcsnew(maxmls+1,k)=0.d0
c          ckssnew(maxmls+1,k)=0.d0
        enddo
        mm=0
        do i=at,at-1+natms
          mm=mm+1
          molxxx(mol,i)=newx(mm)
          molyyy(mol,i)=newy(mm)
          molzzz(mol,i)=newz(mm)
        enddo
        call condense(totatm,ntpfram,ntpguest)
      endif

      end subroutine accept_move
      
      subroutine reject_move
     &(iguest,jguest,insert,delete,displace,swap)
c***********************************************************************
c
c     updates arrays for the rejection according to the type
c     of move.
c
c***********************************************************************
      implicit none
      integer iguest,natms,mol,i,ka
      integer jguest,jnatms,jmol
      logical insert,delete,displace,swap

      mol = locguest(iguest)
      natms = numatoms(mol)
      delE(:)=0.d0
      if (insert) then
c       reset engsic
        engsic = engsicorig
      else if (delete)then 
c       reset engsic
        engsic = engsicorig
      elseif (displace)then
c       nothing yet..
      elseif(swap)then
        jmol=locguest(jguest)
        jnatms = numatoms(jmol)
        engsic=engsicorig
        elrc_mol=origelrc_mol

        do i=1,natms
          ka=ltpsit(mol,i)
          numtyp(ka)=numtyp(ka)+1
          numtyp_mol(mol,ka)=numtyp_mol(mol,ka)+1
          if(lfzsite(mol,i).ne.0)then
            numfrz(ka)=numfrz(ka)+1
            numfrz_mol(mol,ka)=numfrz_mol(mol,ka)+1
          endif
        enddo
        do i=1,jnatms
          ka=ltpsit(jmol,i)
          numtyp(ka)=numtyp(ka)-1
          numtyp_mol(jmol,ka)=numtyp_mol(jmol,ka)-1
          if(lfzsite(jmol,i).ne.0)then
            numfrz(ka)=numfrz(ka)-1
            numfrz_mol(jmol,ka)=numfrz_mol(jmol,ka)-1
          endif
        enddo
      endif 

      end subroutine reject_move

      end module mc_moves
