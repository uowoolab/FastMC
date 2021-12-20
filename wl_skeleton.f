      module wang_landau
      use utility_pack
      use mc_moves
      use ewald_module

      contains
      subroutine wang_landau_sim
     &(idnode,mxnode,imcon,keyfce,alpha,rcut,drewd,totatm,ntpguest,
     &ntpfram,volm,statvolm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,sumchg,maxmls,surftol,overlap,newld,outdir,levcfg,cfgname,
     &wlprec,prectol,
     &nnumg,temp,beta,mcsteps,eqsteps,flatcoeff,visittol,
     &maxn,minn,ebins)
c*****************************************************************************
c     
c     Skeleton so that at compile time one doesn't have to depend on
c     the mpfun library if gcmc is needed only
c     PB - 19/12/2021
c
c*****************************************************************************
      implicit none
      character*8 outdir
      character*1 cfgname(80)      
      integer idnode,mxnode,imcon,keyfce,totatm,ntpguest,ntpfram,nnumg
      integer kmax1,kmax2,kmax3,ntpatm,maxvdw,maxmls,newld,levcfg
      integer mcsteps,eqsteps,visittol,maxn,minn,ebins
      real(8) alpha,rcut,drewd,volm,statvolm,epsq,dlrpot,engunit,sumchg
      real(8) surftol,overlap,wlprec,prectol,temp,beta,flatcoeff
      end subroutine wang_landau_sim

      end module wang_landau
