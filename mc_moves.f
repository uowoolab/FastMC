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
c     atoms is populated correctly for 'jatm'.
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

      subroutine insertion
     &(imcon,idnode,iguest,keyfce,alpha,rcut,delr,drewd,totatm,
     &ntpguest,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,delrc,estep,sumchg,chgtmp,engsictmp,maxmls,
     &loverlap,lnewsurf,surftol,overlap,newld)
c***********************************************************************
c
c     inserts a particle in the framework and computes the 
c     energy of an additional particle.
c
c***********************************************************************
      implicit none
      logical loverlap,lnewsurf,latmsurf
      integer i,ik,j,kmax1,kmax2,kmax3,imcon,keyfce
      integer totatm,sittyp,idnode,ntpatm,maxvdw,newld
      integer k,mol,ntpguest,natms,nmols,iguest,ivdw
      integer jatm,ka,aa,ak,ab,l,mxcmls,maxmls
      real(8) drewd,dlrpot,volm,epsq,alpha,rcut,delr
      real(8) chg,engsrp,engunit,sumchg,engsictmp
      real(8) ewld2sum,ewld3sum,vdwsum,chgtmp
      real(8) ewld1eng,ewld2eng,vdweng,delrc_tmp
      real(8) delrc,estep,sig,surftol,req
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
     &(imcon,iguest,totatm,rcut,delr,
     &natms,newx,newy,newz)

c     calculate long range correction to vdw for the insertion
c     of an additional guest
      call gstlrcorrect(idnode,imcon,iguest,keyfce,natms,ntpatm,maxvdw,
     &engunit,delrc,rcut,volm,maxmls,.false.)
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
c       figure out which index contains charges and ltype arrays
c       that match the guest...
        chg=atmchg(mol,i)
        sittyp=ltpsit(mol,i)
c       calc ewald2 interactions
        call ewald2
     & (chg,gstlentry(i),ewld2eng,mol,maxmls,
     &  drewd,rcut,epsq)
        ewld2sum=ewld2sum+ewld2eng
c       calc vdw interactions 

        call srfrce
     & (sittyp,gstlentry(i),mol,maxmls,vdweng,rcut,dlrpot)
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
     &(imcon,idnode,keyfce,iguest,choice,alpha,rcut,delr,drewd,
     &totatm,ntpguest,maxmls,volm,kmax1,kmax2,kmax3,epsq,dlrpot,
     &ntpatm,maxvdw,engunit,vdwsum,ewld2sum,ewld1eng,lsurf,newld,
     &surftol,sumchg,chgtmp,engsictmp,loverlap,overlap,estep,lexisting)
c***********************************************************************
c
c     compute the energy of a guest in it's current position
c
c***********************************************************************
      implicit none
      logical lsurf,loverlap,latmsurf,lexisting
      integer i,ik,j,kmax1,kmax2,kmax3,imcon,keyfce
      integer totatm,sittyp,idnode,ntpatm,maxvdw
      integer k,mol,ntpguest,natms,nmols,iguest,ivdw
      integer jatm,ka,aa,ak,ab,l,itatm,choice
      integer maxmls,mxcmls,newld
      real(8) drewd,dlrpot,volm,epsq,alpha,rcut,delr
      real(8) chg,engsrp,engunit,chgtmp,engsictmp
      real(8) ewld2sum,ewld3sum,vdwsum,overlap
      real(8) ewld1eng,ewld2eng,vdweng,sumchg
      real(8) estep,sig,surftol,req
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
     &(imcon,iguest,totatm,rcut,delr,
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
        chg=atmcharge(itatm)
        sittyp=ltype(itatm)
        call ewald2
     & (chg,gstlentry(i),ewld2eng,mol,maxmls,
     &  drewd,rcut,epsq)
        ewld2sum=ewld2sum+ewld2eng
c       calc vdw interactions
        call srfrce
     & (sittyp,gstlentry(i),mol,maxmls,vdweng,rcut,dlrpot)
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
     &(imcon,idnode,keyfce,iguest,choice,alpha,rcut,delr,drewd,maxmls,
     &totatm,ntpguest,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
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
      integer totatm,sittyp,idnode,ntpatm,maxvdw
      integer k,mol,ntpguest,natms,nmols,iguest,ivdw
      integer jatm,ka,aa,ak,ab,l,itatm,choice,atmadd
      integer iatm,at,imol,mxcmls,maxmls,step,newld
      real(8) drewd,dlrpot,volm,epsq,alpha,rcut,delr
      real(8) chg,engsrp,engunit,engsictmp,chgtmp
      real(8) ewld2sum,ewld3sum,vdwsum,sumchg,overlap
      real(8) ewld1eng,ewld2eng,vdweng,delrc_tmp
      real(8) delrc,estep,sig,surftol,req
      
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
      call gstlrcorrect(idnode,imcon,iguest,keyfce,natms,ntpatm,maxvdw,
     &engunit,delrc,rcut,volm,maxmls,.true.)

      call guestlistgen
     &(imcon,iguest,totatm,rcut,delr,
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
     & (chg,gstlentry(i),ewld2eng,mol,maxmls,
     &  drewd,rcut,epsq)
        ewld2sum=ewld2sum+ewld2eng
c       calc vdw interactions
        call srfrce
     & (sittyp,gstlentry(i),mol,maxmls,vdweng,rcut,dlrpot)
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

      end module mc_moves
