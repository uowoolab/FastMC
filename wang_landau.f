      subroutine wang_landau
     &(idnode,imcon,keyfce,alpha,rcut,delr,drewd,totatm,ntpguest,
     &ntpfram,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,engunit,
     &sumchg,ntpmls,maxmls,surftol,overlap,newld,outdir,levcfg,cfgname)
c*****************************************************************************
c     
c     Main routine for computing the weighted histogram of Wang and
c     Landau
c     PB - 21/08/2017
c
c*****************************************************************************
      use utility_pack
      use wang_landau_module
      use mc_moves

      implicit none
      logical lgchk,loverlap,lnewsurf,lprod
      character*8 outdir
      character*1 cfgname(80)      
      integer idnode,imcon,keyfce,ntpguest,kmax1,kmax2,kmax3
      integer ntpatm,maxvdw,newld,maxmls,totatm,steps,levcfg
      integer natms,iguest,mol,idum,ntpmls,ntpfram
      real(8) alpha,rcut,delr,drewd,volm,epsq,dlrpot,engunit
      real(8) sumchg,surftol,overlap,estep,chgtmp
      real(8) engsictmp,delrc,eng 
      lgchk=.true.
      do while(lgchk)
        ! randomly select guest
        iguest=1
        mol=locguest(iguest)
        natms=numatoms(mol)
        call random_ins(idnode,imcon,natms,totatm,iguest,rcut,delr)
        estep=0.d0
        call insertion
     & (imcon,idnode,iguest,keyfce,alpha,rcut,delr,drewd,totatm,
     & ntpguest,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     & engunit,delrc,estep,sumchg,chgtmp,engsictmp,maxmls,
     & loverlap,lnewsurf,surftol,overlap,newld)
        call accept_move
     & (imcon,idnode,iguest,.true.,.false.,.false.,estep,lnewsurf,
     & delrc,totatm,idum,ntpfram,ntpmls,ntpguest,maxmls,sumchg,
     & engsictmp,chgtmp,newld)
        lgchk=.false.
      enddo
      lprod=.true.
c      call revive
c     &(idnode,totatm,levcfg,lprod,ntpguest,maxmls,
c     &imcon,cfgname,eng,outdir)

      call error(idnode,2316)
      end subroutine wang_landau
