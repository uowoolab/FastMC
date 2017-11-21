      subroutine wang_landau
     &(idnode,imcon,keyfce,alpha,rcut,delr,drewd,totatm,ntpguest,
     &ntpfram,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,engunit,
     &sumchg,maxmls,surftol,overlap,newld,outdir,levcfg,cfgname,
     &wlprec)
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
      integer ntpatm,maxvdw,newld,maxmls,totatm,levcfg
      integer natms,iguest,mol,idum,ntpfram,i
      real(8) alpha,rcut,delr,drewd,volm,epsq,dlrpot,engunit
      real(8) sumchg,surftol,overlap,estep,chgtmp
      real(8) engsictmp,delrc,wlprec
      write(nrite, 
     &"(/'Entering main routine for Wang - Landau calculation',/)")
      write(nrite,
     &"('Initial coefficient set to',f9.5,/)")wlprec
      lgchk=.true.
c     Temporary: exit if more than one guest included in the
c     FIELD/CONTROL files. Make clear that this currently works
c     for estimating the partition function for a single guest
c     at a single temperature.
      if(ntpguest.gt.1)call error(idnode,2318)
      iguest=1
      mol=locguest(iguest)
      natms=numatoms(mol)

      do i=1,ntpguest
        ! compute reduced thermal debroglie wavelength for each guest

        if(guest_insert(i).gt.0)call insert_guests
     &(idnode,imcon,totatm,ntpguest,ntpfram,iguest,guest_insert(i),
     &rcut,delr,sumchg,surftol,overlap,keyfce,alpha,drewd,volm,newld,
     &kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,engunit,delrc,
     &maxmls)
      enddo

      do while(lgchk)
c       randomly select guest
c       chose an MC move to perform
c       accept/reject based on modified acceptance criteria
        call random_ins(idnode,natms,iguest,rcut,delr)
        estep=0.d0
        call insertion
     & (imcon,iguest,keyfce,alpha,rcut,delr,drewd,totatm,
     & volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     & engunit,delrc,estep,sumchg,chgtmp,engsictmp,maxmls,
     & loverlap,lnewsurf,surftol,overlap,newld)
        call accept_move
     & (iguest,.true.,.false.,.false.,lnewsurf,
     & delrc,totatm,idum,ntpfram,ntpguest,maxmls,sumchg,
     & engsictmp,chgtmp,newld)
        lgchk=.false.
      enddo
      lprod=.true.
      call revive
     &(totatm,levcfg,lprod,ntpguest,maxmls,
     &imcon,cfgname,0.d0,outdir)

      call error(idnode,2316)
      end subroutine wang_landau

      logical function convergence_check()
c**********************************************************************
c
c     check to see if the histogram is converged. 
c     PB - 15/11/2017
c
c**********************************************************************
      implicit none
       
      convergence_check=.false.
      return
      end function convergence_check

      real(8) function adjust_factor(f)
c**********************************************************************
c
c     Function to decrease the density of states scaling factor
c     PB - 15/11/2017
c
c**********************************************************************
      implicit none
      real(8) f
      adjust_factor = sqrt(f)
      return
      end function adjust_factor

      logical function accept_wl_move
     &(idnode, iguest, insert, delete, displace)
c**********************************************************************
c
c     Acceptance criteria for the Wang-Landau algorithm. Right now
c     this is for ajusting the number of molecules as a macrostate
c     variable.
c     PB - 15/11/2017
c
c**********************************************************************
      implicit none
      logical insert, delete, displace
      integer idnode, iguest

      accept_wl_move=.false.
      return
      end function accept_wl_move

