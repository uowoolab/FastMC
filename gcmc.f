c*************************************************************
c      main gcmc program                                     *
c*************************************************************
c      personal reference                                    *
c      imcon                                                 *
c      0  -  no periodic boundaries                          *
c      1  -  cubic boundary conditions  (b.c.)               *
c      2  -  orthorhombic b.c.                               *
c      3  -  parallelpiped b.c.                              *
c      4  -  truncated octahedral b.c.                       *
c      5  -  rhombic dodecahedral b.c                        *
c      6  -  x-y parallelogram b.c. no periodicity in z      *
c      7  -  hexagonal prism b.c.                            *
c                                                            *
c                                                            *
c                                                            *
c                                                            *
c     avgwindow will reset to 0's once the average has been  *
c     carried over to varwindow.                             *
c     Varwindow will be a rolling variance calculation of    *
c     the windowed averages.                                 * 
c     currently the order is:                                *
c     chainstats(1) = production gcmc step count             *
c     chainstats(2) = rolling average number of guests <N>   *
c     chainstats(3) = rolling average energy <E>             *
c     chainstats(4) = rolling average for energy*guests <EN> *
c     chainstats(5) = rolling average for N^2 <N2>           *
c     chainstats(6) = rolling average for (energy)^2 <E2>    *
c     chainstats(7) = rolling average for <N surface>        *
c     chainstats(8) = rolling average for <exp(-E/kb/T)>     *
c     chainstats(9) = rolling stdev value for <N>            *
c     chainstats(10) = rolling stdev value for <E>           *
c     chainstats(11) = rolling stdev value for <EN>          *
c     chainstats(12) = rolling stdev value for <N2>          *
c     chainstats(13) = rolling stdev value for <E2>          *
c     chainstats(14) = rolling stdev value for <N surface>   *
c     chainstats(15) = rolling stdev value for <Q_st>        *
c     chainstats(16) = rolling stdev value for <C_v>         *
c     chainstats(17) = rolling stdev value for <exp(-E/kb/T)>*
c*************************************************************


      use utility_pack
      use readinputs
      use vdw_module
      use flex_module
      use ewald_module
      use mc_moves
      use wang_landau

      implicit none

      character*1 cfgname(80)      
      character*8 outdir,localdir
      character*2 mnth
      character(8):: date
      character(10) :: time
      character(5) :: zone
      character*25 outfile
      character*80 command
      character*25 outfile2
      character*18 outfile4
      logical lgchk,lspe,ljob,lprob,widchk
      logical lfuga, loverlap, lwidom, lwanglandau
      logical insert,delete,displace,lrestart,laccsample
      logical jump, flex, swap, switch
      logical tran, rota
      logical accepted,production,jobsafe,lnumg
      logical tick_tock_cycles(5)
      logical lnewsurf, linitsurf, lnewsurfj,linitsurfj
      integer nfo(8)
      integer accept_ins,accept_del,accept_disp,totaccept,nnumg,nhis
      integer accept_jump,accept_flex,accept_swap,scell_factor
      integer accept_switch,minchk,iter
      integer accept_tran,accept_rota
      integer ins_count,del_count,disp_count,buffer,levcfg,nstat
      integer jump_count,flex_count,swap_count,np,widcount
      integer tran_count,rota_count,switch_count
      integer totatm,jatm,imol,iatm,ntprob,ntpsite
      integer newld,gcmccount,prodcount,globalprod
      integer ibuff,iprob,cprob,prevnodes,nwind,rollstat,ci
      integer k,p,i,j,ik,jj,kk,l,ntpatm,gridsize,nwindsteps
      integer kmax1,kmax2,kmax3,ntpvdw,maxvdw,maxmls,mxnode
      integer ntpfram,randchoice,nmols,ngsts
      integer ichoice,jchoice,iswap,origtotatm,maxn,minn,ebins
      integer mxatm,imcon,keyfce,mxatyp,mcsteps,eqsteps,vstat
      integer iguest,ntpguest,ntpmls,natms,mol,rollcount
      integer jguest,jmol,jnatms,jnmols,ifram,nframmol
      integer ins_rollcount,nswapguest,num_swaps,swap_max
      integer swap_guest_max,idnode,visittol,nwidstep
      integer ngrida,ngridb,ngridc,istat,avcount,wngrida,wngridb,wngridc
      integer totalguests,globalnguests,globalsteps
      integer, allocatable :: fwksumbuff(:)
      real(8), allocatable :: gridbuff(:)
      real(8), dimension(10) :: celprp
      real(8), dimension(10) :: ucelprp
      real(8), dimension(9) :: rucell
      real(8) sumchg,chgtmp,engsictmp
      real(8) ewld1old
      real(8) det,engcpe,engacc,engac1,drewd,epsq,statvolm,volm
      real(8) engsrp,randmov,rande,delta
      real(8) stdQst,stdCv,rotangle,delrc
      real(8) tzero,timelp,engunit,rvdw,temp,beta
      real(8) dlrpot,rcut,eps,alpha,delr,delrdisp,gpress
      real(8) init,wlprec
      real(8) ewld1eng,flatcoeff
c     DEBUG
      real(8) pewld1,pewld2,pelrc,pvdw
c     DEBUG
      real(8) ecoul,evdw
      real(8) ecoulg,evdwg,overlap,surfmol
      real(8) ewld2sum,vdwsum,comx,comy,comz,ewaldaverage
      real(8) dmolecules,molecules2,energy2,Q_st,C_v
      real(8) weight,tw,spenergy,dlpeng,estep,thrd,twothrd
      real(8) estepi,estepj,sumweight,sweight
      real(8) E,aN,EN,E2,N2,H_const,avgH,stdH
      real(8) avgN,stdN,avgE,stdE,avgEN,stdEN,avgN2,stdN2,avgE2,stdE2
      real(8) avgNF,NF,avgCv,avgQst,stdNF,surftol,guest_toten
      real(8) griddim,griddima,griddimb,griddimc,grvol,prectol
      real(8) a,b,c,q1,q2,q3,q4 
      real(8) delE_fwk
      integer m, indatm, indfwk, gstidx
      logical isguest
c Cumulative move probabilities, remeber to zero and include in normalising
      real(8) mcinsf,mcdelf,mcdisf,mcswpf,mcflxf,mcjmpf,mcmvnorm
      real(8) mctraf,mcrotf,mcswif,mcmvsum
      real(8) disp_ratio,tran_ratio,rota_ratio,tran_delr
      real(8) rota_rotangle,jumpangle
      integer, dimension(3) :: gridfactor
      integer fwk_step_magnitude

      data lgchk/.true./,insert/.false./,delete/.false./,
     &lwidom/.false./,lwanglandau/.false./,
     &displace/.false./,accepted/.false./,production/.false./
      data jump/.false./,flex/.false./,swap/.false./,switch/.false./
      data tran/.false./,rota/.false./,laccsample/.false./
      data lspe/.false./,ljob/.false./,jobsafe/.true./,lrestart/.false./
     &,lnumg/.false./
      data lfuga/.false./,loverlap/.false./ ! change to false
      data linitsurf/.false./,lnewsurf/.false./
      data linitsurfj/.false./,lnewsurfj/.false./
      data accept_ins,accept_del,accept_disp,totaccept/0,0,0,0/
      data accept_jump, accept_flex, accept_swap, accept_switch/0,0,0,0/
      data ins_count,del_count,disp_count,gcmccount,prevnodes/0,0,0,0,0/
      data jump_count, flex_count, swap_count, switch_count/0,0,0,0/
      data tran_count, rota_count, accept_tran, accept_rota/0,0,0,0/
c     Default these to grand canonical.. can turn 'off' in CONTROL file
      data mcinsf/0.3333333/
      data mcdelf/0.3333333/
      data mcdisf/0.3333333/
      data mcjmpf, mcflxf, mcswpf, mctraf, mcrotf, mcswif/0,0,0,0,0,0/

      integer, parameter, dimension(3) :: revision = (/1, 4, 6 /)
      nwidstep=20
c     TODO(pboyd): include error checking for the number of guests
c     coinciding between the CONTROL file and the FIELD file.
      tw=0.d0
      ewaldaverage=0.d0
      ewld1old=0.d0
      guest_toten=0.d0
      ecoul=0.d0
      evdw=0.d0
c     default overlap to 2.0 Angstroms instead of 0.
      overlap=2.d0
c     surface tolerance default set to -1 angstroms (off)
      surftol=0.d0
      ntpatm=0
      mxatm=0
      mxegrd=0
      gridsize=0
      ibuff=0
      nnumg=1
      rollcount=0
      ins_rollcount=0
      avcount = 0
c     averaging window to calculate errors
      nwindsteps=100000
c     default flatness tolerance for Wang-Landau simulations.
      flatcoeff=0.7
c     default min number of visits for each histogram bin
c     in Wang-Landau simulations.
      visittol=10000
c     default Wang-Landau DOS filling coefficient is e. 
      wlprec = dexp(1.d0)
c     default Wang-Landau DOS precision tolerance
      prectol=1e-8
c     default maximum number of guest species to sample with WL sim
      maxn=200
      minn=0
c     ebins is the number of bins to store the energy in during
c     flat histogram sampling
      ebins=100
c     scoping issues
      delrc = 0

c     Default target acceptance ratios of 0.5      
      disp_ratio = 0.5d0
      tran_ratio = 0.5d0
      rota_ratio = 0.5d0

      thrd=1.d0/3.d0
      twothrd=2.d0/3.d0

c default length of a side of a grid box in angstroms (in CONTROL)
      griddim=0.1d0
c default supecell folding for grid points (in CONTROL)
      gridfactor = (/ 1, 1, 1 /)
      
c     global number of production steps over all nodes
      globalprod=0
c     local production steps (after equilibrium)
      prodcount=0
c     local insertion steps (after equilibrium) - to keep track 
c     newld is the number of ewald points in reciprocal space
      newld=0
c     mcsteps = number of production steps, energies and molecules
c counted towards final averages
      mcsteps=1
c     eqsteps = number of equilibrium steps, energies and molecules
c ignored
      eqsteps=0
c when running mc cycles keep track of if history to help averaging
      tick_tock_cycles = .false.
      
c     initialize communications
      call initcomms()
      call gsync()
      call timchk(0,tzero)

c     determine processor identities
      call machine(idnode,mxnode)

c     open main output file.
      if(idnode.eq.0)then

        open(nrite,file='OUTPUT')
        write(nrite,
     &"(/,20x,'FastMC version ',i1,'.',i1,'.',i1,/,/)")revision
        call date_and_time(date,time,zone,nfo)
        mnth=date(5:6)
        
        write(nrite,
     &"('Started : ',9x,a2,3x,a9,1x,a4,3x,a2,':',a2,a6,' GMT',/)")
     &date(7:8),month(mnth),date(1:4),time(1:2),time(3:4),zone
        write(nrite,"('Running on ',i4,' nodes',/,/)")mxnode
      endif
      call initscan
     &(idnode,imcon,volm,keyfce,rcut,eps,alpha,kmax1,kmax2,kmax3,lprob,
     & delr,rvdw,ntpguest,ntprob,ntpsite,ntpvdw,maxmls,mxatm,mxatyp,
     & griddim,gridfactor,nwind,nwindsteps)

      maxvdw=max(ntpvdw,(mxatyp*(mxatyp+1))/2)
      call alloc_config_arrays(idnode,mxnode,maxmls,mxatm,mxatyp,
     &volm,ntpguest,rcut,rvdw,delr,nwind)

      call alloc_vdw_arrays(idnode,maxvdw,maxmls,mxatyp)
      call readfield
     &(idnode,ntpvdw,maxvdw,ntpatm,ntpmls,ntpguest,
     &ntpfram,totatm,rvdw,dlrpot,engunit,sumchg)
     
      call readconfig(idnode,mxnode,imcon,cfgname,levcfg,
     &ntpmls,maxmls,totatm,volm,rcut,celprp)

      call alloc_ewald_arrays
     &(idnode,maxmls,kmax1,kmax2,kmax3,rvdw,totatm,maxguest)
c     volume reported in m^3 for statistical calculations
      statvolm=volm*1d-30
      call invert(cell,rcell,det)

      if(lprob)then
c      calculate number of grid points in a,b,and c directions
c      calculate the volume of a grid point (differs from griddim^3)
        griddima = celprp(1) 
     &/ceiling(celprp(1)/gridfactor(1)/griddim)
     &/gridfactor(1) 
        griddimb = celprp(2) 
     &/ceiling(celprp(2)/gridfactor(2)/griddim)
     &/gridfactor(2) 
        griddimc = celprp(3) 
     &/ceiling(celprp(3)/gridfactor(3)/griddim)
     &/gridfactor(3) 
c NB these are allocated before assigning to guests in readconfig
c        ngrida=gridfactor(1)*ceiling(celprp(1)/(griddim*gridfactor(1)))
c        ngridb=gridfactor(2)*ceiling(celprp(2)/(griddim*gridfactor(2)))
c        ngridc=gridfactor(3)*ceiling(celprp(3)/(griddim*gridfactor(3)))
        ngrida=ceiling(celprp(1)/(griddima*gridfactor(1)))
        ngridb=ceiling(celprp(2)/(griddimb*gridfactor(2)))
        ngridc=ceiling(celprp(3)/(griddimc*gridfactor(3)))
        grvol = celprp(10)/dble(ngrida*ngridb*ngridc)
        gridsize=ngrida*ngridb*ngridc
        if(idnode.eq.0)then
          write(nrite,"(/,' Probability grid size    :',i8,i8,i8)")
     &ngrida,ngridb,ngridc
          write(nrite,"(/,' Grid voxel dimensions (A):',
     &f9.4,f9.4,f9.4)")griddima,griddimb,griddimc
          write(nrite,"(/,' Grid voxel volume (A^3)  :  ',f12.5)")
     &grvol
        endif
        allocate(gridbuff(gridsize))
        call alloc_prob_arrays(idnode,ntpguest,ntpsite,ntprob,gridsize)
        grvol=celprp(1)/dble(ngrida)*celprp(2)/dble(ngridb)*
     &  celprp(3)/dble(ngridc)
      else
c       this is in case we run into allocation problems later on

        ntprob=0
        gridsize=1
        call alloc_prob_arrays(idnode,ntpguest,ntpsite,ntprob,gridsize)
      endif

c     produce unit cell for folding purposes
      do i=1, 9
        ci = ceiling(dble(i)/3.d0)
        ucell(i)=cell(i)/dble(gridfactor(ci))
      enddo
      call invert(ucell,rucell,det)
      call dcell(ucell,ucelprp)
      scell_factor = gridfactor(1)*gridfactor(2)*gridfactor(3)
      do i=1,ntpfram
        ifram=locfram(i)
        nframmol=nummols(ifram)
        if (scell_factor.gt.nframmol)then
          write(nrite,"(/a9,a44,/,a50,/,a45,/)")
     &      "WARNING: ","The grid factor setting in the CONTROL file ",
     &      "is larger than the number of framework molecules. ",
     &      "This may result in useless probability plots."
        endif
      enddo
      call readcontrol(idnode,lspe,temp,ljob,mcsteps,eqsteps,
     &ntpguest,lrestart,laccsample,lnumg,nnumg,nhis,
     &mcinsf,mcdelf,mcdisf,mcjmpf,mcflxf,mcswpf,swap_max,mcswif,
     &mctraf,mcrotf,disp_ratio,tran_ratio,rota_ratio,lfuga,overlap,
     &surftol,n_fwk,l_fwk_seq,fwk_step_max,fwk_initial,lwidom,nwidstep,
     &lwanglandau,wlprec,flatcoeff,visittol,prectol,maxn,minn,ebins)
c     square the overlap so that it can be compared to the rsqdf array
      overlap = overlap**2
c     square the surface tolerance so that it can be compared to the
c     rsqdf array
      call fugacity(idnode,lfuga,temp,ntpguest)
      if((.not.lwidom).and.(.not.lwanglandau))then
        if(idnode.eq.0)then
          write(nrite,'(a6,5x,a1,15x,a12,7x,a12)')'guest','y',
     &   'fugacity/bar','pressure/bar'
          do i=1, ntpguest
              write(nrite,'(i6,1x,f9.3,2(f19.6))') i, gstmolfract(i), 
     &   gstfuga(i)/1.E5, gstpress(i)/1.E5
          enddo
        endif
      endif

c     FLEX
      if(n_fwk.gt.0)then
        allocate(fwksumbuff(n_fwk))
        call flex_init(idnode, mxatm, imcon, ntpmls, maxmls,
     &totatm, rcut, celprp, ntpguest, volm, mxnode)
      endif
c Normalise the move frequencies
      mcmvnorm = mcinsf+mcdelf+mcdisf+mcjmpf+mcflxf+mcswpf+mctraf+mcrotf
     &+mcswif
      if(mcmvnorm.eq.0)then
        if(idnode.eq.0)
     &write(nrite,"(/,'No move frequencies specified, defaulting to 
     &Grand Canonical')")
        mcmvnorm = 1
        mcinsf = 1.d0/3.d0
        mcdelf = 1.d0/3.d0
        mcdisf = 1.d0/3.d0
        mcjmpf = 0.d0
        mcflxf = 0.d0
        mcswpf = 0.d0
        mctraf = 0.d0
        mcrotf = 0.d0
        mcswif = 0.d0
      endif
      if(idnode.eq.0)
     &write(nrite,"(/,'Normalised move frequencies:',/,
     &a18,f7.3,a18,f7.3,a18,f7.3,/,a18,f7.3,a18,f7.3,a18,f7.3,/,
     &a18,f7.3,a18,f7.3,a18,f7.3/)")
     &'insertion:',mcinsf/mcmvnorm, 'deletion:',mcdelf/mcmvnorm,
     &'displacement:',mcdisf/mcmvnorm,'jumping:',mcjmpf/mcmvnorm,
     &'flexing:',mcflxf/mcmvnorm,'swapping:',mcswpf/mcmvnorm,
     &'translation:',mctraf/mcmvnorm,'rotation:',mcrotf/mcmvnorm,
     &'switching:',mcswif/mcmvnorm
c Now we normalise
      mcmvsum=0.d0
      mcmvsum=mcmvsum+mcinsf
      if(mcinsf.gt.0.d0)mcinsf=mcmvsum/mcmvnorm
      mcmvsum=mcmvsum+mcdelf
      if(mcdelf.gt.0.d0)mcdelf=mcmvsum/mcmvnorm
      mcmvsum=mcmvsum+mcdisf
      if(mcdisf.gt.0.d0)mcdisf=mcmvsum/mcmvnorm
      mcmvsum=mcmvsum+mcjmpf
      if(mcjmpf.gt.0.d0)mcjmpf=mcmvsum/mcmvnorm
      mcmvsum=mcmvsum+mcflxf
      if(mcflxf.gt.0.d0)mcflxf=mcmvsum/mcmvnorm
      mcmvsum=mcmvsum+mcswpf
      if(mcswpf.gt.0.d0)mcswpf=mcmvsum/mcmvnorm
      mcmvsum=mcmvsum+mctraf
      if(mctraf.gt.0.d0)mctraf=mcmvsum/mcmvnorm
      mcmvsum=mcmvsum+mcrotf
      if(mcrotf.gt.0.d0)mcrotf=mcmvsum/mcmvnorm
      mcmvsum=mcmvsum+mcswif
      if(mcswif.gt.0.d0)mcswif=mcmvsum/mcmvnorm

      if(idnode.eq.0)
     &write(nrite,"(/,'Target acceptance ratios:',/,
     &a18,f7.3,a18,f7.3,a18,f7.3,/)")
     &'displacement:',disp_ratio,'translation:',tran_ratio,
     &'rotation:',rota_ratio

      if(.not.lspe)then
        if(idnode.gt.0)call sleep(idnode+1)
        init=duni(idnode)

c     initialize jobcontrol file
        if(idnode.eq.0)then
          open(205,file='jobcontrol.in')
          close(205)
        endif
      

c       initialize rotangle 
        rotangle=pi/3.d0
c        rotangle=2.d0*pi
        jumpangle=2.d0*pi
   
c==========================================================================        
c       if restart requested then descend into the branch
c       and read the REVIVE and REVCON for the appropriate
c       arrays
c==========================================================================

        if(lrestart)then
c         do a check to see if the number of branch directories
c         matches the number of nodes.  adjust accordingly depending
c         on the situation.  
c         This is a rather shitty way of scanning branch directories
c         if other files are called "branch" in the working directory
c         problems happen.

          call revscan(idnode,prevnodes)
          if(prevnodes.eq.0)then
            if(idnode.eq.0)write(nrite,"(/,a85,/)")
     &"No branches found in the working directory, starting from
     &the CONFIG file provided"
            write(outdir,"('branch',i2.2)")idnode+1
            write(command,"('mkdir 'a8)")outdir
            call system(command)
c         if the restart calculation has more nodes than the previous
c         calculation
          elseif(idnode+1.gt.prevnodes)then
c           create the new directory
            write(outdir,"('branch',i2.2)")idnode+1
            write(command,"('mkdir 'a8)")outdir
            call system(command)

c           apply the values from other branches to the new node.
            write(localdir,"('branch',i2.2)")(mod(idnode+1,prevnodes)+1)
            
            call revread(localdir,production,ntpmls,totatm,ntpguest)

          elseif(idnode+1.le.prevnodes)then 
c         read the values from the first mxnodes
            write(localdir,"('branch',i2.2)")idnode+1
            outdir=localdir
            call revread(localdir,production,ntpmls,totatm,ntpguest)
          endif
c     write warning if the number of nodes does not correspond with the 
c     previous calculation
          if(idnode.eq.0.and.prevnodes.ne.mxnode)then
           write(nrite,"(/,3x,a35,i2,a36,i2,/)")"WARNING - previous 
     &calculation had",prevnodes," branches while the current job has ",
     &mxnode
           if(mxnode.lt.prevnodes)write(nrite,"(3x,a9,i2,a4,i2,a17,/)")
     &"Branches ",mxnode+1,"  - ", prevnodes,"  will be ignored"
           if(mxnode.gt.prevnodes)
     & write(nrite,"(3x,a9,i2,a4,i2,a30,i2,a3,i2,/)")
     &"Branches ", prevnodes+1,"  - ", mxnode,"  will take data from
     & branches ",1," - ",prevnodes
          endif
          if(production)then
            if(eqsteps.gt.0)then
              production=.false.
              if(idnode.eq.0)write(nrite,"(3x,a41,i7,a6,/)")
     &"Production averaging will continue after ",eqsteps," steps"
            else
              if(idnode.eq.0)write(nrite,"(3x,a59,/)")
     &"Production averaging will start at the beginning of the run"
            endif
          endif

c     local output for each node
c     This may only work on system specific machines
        else 
          write(outdir,"('branch',i2.2)")idnode+1
          write(command,"('mkdir ',a8)")outdir
          call system(command)
        endif
c=========================================================================
c      open necessary archiving files
c=========================================================================

        if(ntpguest.gt.1)then
          do i=1,ntpguest
            if(lnumg)then
              write(outfile,"(a8,'/numguests',i2.2,'.out')")outdir,i
              open(400+i,file=outfile)
            endif
            if(abs(nhis).gt.0)then
              write(outfile4,"(a8,'/his',i2.2,'.xyz')")outdir,i
              open(500+i,file=outfile4)
            endif
          enddo
        else
          if(lnumg)then
            outfile=outdir // '/numguests.out'
            open(401,file=outfile)
          endif
          if(abs(nhis).gt.0)then
            outfile4=outdir // '/his.xyz'
            open(501,file=outfile4)
          endif
        endif
        outfile2=outdir // '/runningstats.out'
        open(202,file=outfile2)
c        outfile3=outdir // '/energies.out'
c        open(203,file=outfile3)
      endif

c     ins,del,dis store the last move made for each guest type
 
c     debugging.. need to see if all information gets to each
c     node
c      write(debug,"('debug',i2.2)")idnode
c      open(999,file=debug)
c      write(999,'(a20,i2,/)')'data for node ',idnode
c      write(999,'(a20,/,3f16.5,/,3f16.5,/,3f16.5)')'cell vectors',
c     & (cell(i),i=1,3),(cell(j),j=4,6),(cell(k),k=7,9)
      call guest_exclude(ntpguest)
      call condense(totatm,ntpfram,ntpguest)

c      write(999,"('atomic information',/)")
c      do i=1,totatm
c        write(999,'(3x,a4,3f16.6)')atomname(i),xxx(i),yyy(i),zzz(i)
c      enddo

c      write(999,"('some other info ',/)")
c      write(999,"('production?',3x,l)")production
c      write(999,"('number of guests',3x,i3)")ntpguest
c      write(999,"('number of equilibrium steps',3x,i9)")eqsteps
c      write(999,"('number of production steps ',3x,i9)")mcsteps
c      write(999,"('number of probability plots',3x,i9)")ntprob
c      close(999)
      engsrp=0.d0
      engcpe=0.d0
      engacc=0.d0
      engac1=0.d0

c     beta is a constant used in the acceptance criteria for the gcmc
c     moves
      beta=1.d0/(kboltz*temp)
c     this is the relative dielectric constant. default is 1
c     we don't need this....

      epsq=1.d0

c     create ewald interpolation arrays 
      call erfcgen(keyfce,alpha,rcut,drewd)
c     populate ewald3 arrays
      call single_point
     &(imcon,keyfce,alpha,drewd,rcut,delr,totatm,ntpfram,
     &ntpguest,ntpmls,volm,newld,kmax1,kmax2,kmax3,epsq,
     &dlrpot,ntpatm,maxvdw,spenergy,vdwsum,ecoul,dlpeng,
     &maxmls,surftol)
      if(idnode.eq.0)then
c        write(nrite,'(/,/,a35,f22.6)')'Configurational energy:
c     & ',spenergy
c        write(nrite,'(/,a35,f22.6)')'Initial framework energy :',
c     &(evdw+ecoul)/engunit
        do iguest=1,ntpguest
          mol=locguest(iguest)
          nmols=nummols(mol)
          spenergy=0.d0
          pewld1=0.d0; pewld2=0.d0; pelrc=0.d0;pvdw=0.d0
          do jmol=1,ntpmls
            ik=loc2(mol,jmol)
c            THIS EWALD1 AND 3 CALCULATION SHOULD BE THE SAME,
c            BUT THEY'RE CURRENTLY NOT!
c            print*, jmol,ewald1en(10)/engunit,ewald3en(jmol)/engunit
            pewld1=pewld1+ewald1en(ik)
            pewld2=pewld2+ewald2en(ik)
            pelrc=pelrc+elrc_mol(ik)
            pvdw=pvdw+vdwen(ik)
            spenergy=spenergy+(
     &ewald1en(ik) +
     &ewald2en(ik) +
     &elrc_mol(ik) +
     &vdwen(ik))/engunit
          enddo
c          print *,"EWALD1 ",pewld1/engunit
c          print *,"EWALD2 ",pewld2/engunit
c          print *,"VDW    ",pvdw/engunit
c          print *,"DELRC  ",pelrc/engunit
          if(nmols.gt.0)spenergy=spenergy-nmols*ewald3en(mol)/engunit
          energy(mol)=spenergy
          write(nrite,'(a23,i3,a9,f22.6)')'Initial guest ',
     &iguest,' energy :',spenergy
        enddo
        write(nrite,'(a35,f22.6)')'van der Waals energy :',
     &(vdwsum+elrc)/engunit
        write(nrite,'(a35,f22.6)')'Electrostatic energy :',
     &ecoul/engunit
        write(nrite,'(a35,f22.6)')'Energy reported by DL_POLY :',
     &dlpeng/engunit

      endif
     
      if(lspe)then
        if(lspe)lgchk=.false.
        call error(idnode,0)
      endif
   
      call timchk(1,tzero)
      if((eqsteps.eq.0).and.(.not.ljob))production=.true.

c******************************************************************
c
c       Widom insertion method begins
c
c******************************************************************
      if(lwidom)then
        lgchk=.false.
c        gcmccount=0 
        wngrida=ceiling(celprp(1)/(griddima))
        wngridb=ceiling(celprp(2)/(griddimb))
        wngridc=ceiling(celprp(3)/(griddimc))
        do iguest=1,ntpguest
          call widom_grid(idnode,iguest,nwidstep,wngrida,wngridb,
     &wngridc,rotangle,imcon,keyfce,alpha,rcut,delr,drewd,totatm,volm,
     &kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,engunit,delrc,maxmls,
     &overlap,newld,surftol,sumchg,rucell,prodcount,widcount,temp,
     &ngrida,ngridb,ngridc)
c          mol=locguest(i) 
c          istat=1+16*(iguest-1)
c          widchk=.true.
c          widcount=0
c          do while(widchk)
c            ins(iguest)=1
c            ins_count=ins_count+1
c            delE=0.d0
c            call random_ins(idnode,natms,iguest,rcut,delr)
c            estep = 0.d0
c            call insertion
c     & (imcon,iguest,keyfce,alpha,rcut,delr,drewd,totatm,
c     & volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
c     & engunit,delrc,estep,sumchg,chgtmp,engsictmp,maxmls,
c     & loverlap,lnewsurf,surftol,overlap,newld)
c            gcmccount = gcmccount + 1
c            call reject_move
c     &  (iguest,0,.true.,.false.,.false.,.false.)
c            if(.not.loverlap)then
c              if(lprobeng(iguest))then
c                call storeenprob(iguest,0,rucell,ntpguest,ngrida,
c     &            ngridb,ngridc,estep)
c              endif
c              H_const=dexp(-1.d0*estep/kboltz/temp)
cc             no rolling average here, just div by widcount at the end
c              chainstats(istat+7) = chainstats(istat+7)+H_const
c              widcount = widcount + 1
c            endif
c            if(widcount.ge.mcsteps)widchk=.false.
c          enddo
c          chainstats(istat+7) = chainstats(istat+7)/dble(widcount)
        enddo
c        prodcount=widcount
c        chainstats(1) = dble(widcount)
      endif
      if (lwanglandau)then
        lgchk=.false.
        call wang_landau_sim
     &(idnode,mxnode,imcon,keyfce,alpha,rcut,delr,drewd,totatm,ntpguest,
     &ntpfram,volm,statvolm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,sumchg,maxmls,surftol,overlap,newld,outdir,levcfg,cfgname,
     &wlprec,mcinsf,mcdelf,mcdisf,mcjmpf,mcflxf,mcswpf,mctraf,prectol,
     &mcrotf,mcswif,nnumg,temp,beta,mcsteps,eqsteps,flatcoeff,visittol,
     &maxn,minn,ebins)
      endif


c******************************************************************
c
c       gcmc begins
c
c******************************************************************

c     start delrdisp as delr (read from the CONTROL file)
      minchk=min(400,nnumg)
      delrdisp=delr
      tran_delr = delr
      rota_rotangle = rotangle
      call timchk(0,timelp)
c     DEBUG
c      call test
c     &(imcon,idnode,keyfce,alpha,rcut,delr,drewd,totatm,dlrpot,newld,
c     &ntpguest,volm,kmax1,kmax2,kmax3,epsq,ntpatm,maxvdw,surftol,
c     &engunit,ntpfram,maxmls,outdir,cfgname,levcfg,overlap)
c      call error(idnode,0)
c     END DEBUG
      if(lgchk)then 
        do i=1,ntpguest
          if(guest_insert(i).gt.0)then
c           probably should write to a different file than
c           'runningstats.out', but node specific so will do for now.
            imol=locguest(i)
            write(202,
     &"(1x,'Inserting ',i6,' guests of type ',i3)")
     &guest_insert(i),i
            call insert_guests
     &(idnode,imcon,totatm,ntpguest,ntpfram,i,guest_insert(i),
     &rcut,delr,sumchg,surftol,overlap,keyfce,alpha,drewd,volm,newld,
     &kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,engunit,delrc,
     &maxmls,iter)
            write(202,"(1x,'Successful Insertion of ',i6,
     &' guests of type ',i3)")
     &nummols(imol),i
            write(202,"(1x,i9,' trials. Success rate: ',f6.2,' %'/)")
     &iter,dble(nummols(imol))/dble(iter) * 100.d0
          endif
        enddo
      endif

      do while(lgchk)
        gcmccount=gcmccount+1
c     every so often, update the total number of prodcounts
c     across all nodes
        if(mod(gcmccount,minchk).eq.0)then

c         check jobcontrol
          call jobcheck(idnode,jobsafe,ljob,production)

          if(.not.jobsafe)then
            call gstate(jobsafe)
            lgchk=.false. 
            call error(idnode,-2312)
          endif
          if (production)then
            call gisum2(prodcount,1,globalprod)
            if(mcsteps.lt.0)then
c             Check for cycles here
              tick_tock_cycles = cshift(tick_tock_cycles, -1)
              tick_tock_cycles(1) = .true.
c             sum up for all guests
              totalguests = 0
              do i=1,ntpguest
                mol=locguest(i)
                totalguests = totalguests + nummols(mol)
              enddo
              call gisum2(totalguests,1,globalnguests)
c             If the total number of production steps over all nodes
c             is less than the total number of guests (plus one)
c             multiplied by desired cycles flag this check as not done
              if((globalprod).lt.(-mcsteps*(globalnguests+1)))then
                tick_tock_cycles(1) = .false.
              endif
c             safe to end if every check passes
              if(all(tick_tock_cycles))then
                if(idnode.eq.0)write(nrite, "('Completed at least ',
     &i10,' production cycles for each guest')")
     &-mcsteps
                lgchk=.false.
              endif
            else
              if(globalprod.ge.mcsteps)lgchk=.false.
            endif
          elseif(.not.ljob)then
c           not production yet; test if we are equilibratied
            if(eqsteps.gt.0)then
              if(gcmccount.ge.eqsteps)then
                if(idnode.eq.0)write(nrite, "('Completed at least ',
     &i10,' equlibration steps',/,'Starting production at ',i10)")
     &eqsteps, gcmccount
                production=.true.
              endif
            elseif(eqsteps.lt.0)then
              call gisum2(gcmccount,1,globalsteps)
c             Check for cycles here
              tick_tock_cycles = cshift(tick_tock_cycles, -1)
              tick_tock_cycles(1) = .true.
c             sum up for all guests
              totalguests = 0
              do i=1,ntpguest
                mol=locguest(i)
                totalguests = totalguests + nummols(mol)
              enddo
              call gisum2(totalguests,1,globalnguests)
c             If the total number of production steps over all nodes
c             is less than the total number of guests (plus one)
c             multiplied by desired cycles flag this check as not done
              if((globalsteps).lt.(-eqsteps*(globalnguests+1)))then
                tick_tock_cycles(1) = .false.
              endif
c             safe to end if every check passes
              if(all(tick_tock_cycles))then
                if(idnode.eq.0)write(nrite, "('Completed at least ',
     &i10,' equilibration cycles for each guest',/,
     &'Starting production at',i10,' over all nodes')")
     &-eqsteps, globalsteps
                production=.true.
c               reset the cycles for production
                tick_tock_cycles=.false.
              endif
            endif
          endif
        endif

c       randomly choose a guest type to move 
        if(ntpguest.gt.1)then
          iguest=floor(duni(idnode)*ntpguest)+1
        elseif(ntpguest.eq.1)then
          iguest=1
        endif
        mol=locguest(iguest)
        natms=numatoms(mol)
        nmols=nummols(mol)

        if(mod(gcmccount,1000).eq.0)then
          call revive
     &(totatm,0,production,ntpguest,ntpmls,imcon,cfgname,
     &   delE(mol),outdir)
        endif

c Randomly decide which MC move to do
        randmov=duni(idnode)

        if(randmov.lt.mcinsf)then
          insert = .true.
        elseif(randmov.lt.mcdelf)then
          delete = .true.
        elseif(randmov.lt.mcdisf)then
          displace = .true.
        elseif(randmov.lt.mcjmpf)then
          jump = .true.
        elseif(randmov.lt.mcflxf)then
          flex = .true.
        elseif(randmov.lt.mcswpf)then
          swap = .true.
        elseif(randmov.lt.mctraf)then
          displace = .true.
          tran = .true.
        elseif(randmov.lt.mcrotf)then
          displace = .true.
          rota = .true.
        elseif(randmov.lt.mcswif)then
          switch = .true.
        else
c Failover displace -- shouldn't reach here
          displace=.true.
        endif
        if((nhis.ne.0).and.(mod(gcmccount,abs(nhis)).eq.0))then
          write(202,'(a35,f20.15,a15,f15.10,a15,f15.10)')
     &'displacement acceptance ratio: ',
     &(dble(accept_disp)/dble(disp_count)),
     &'delr: ',delrdisp,'angle: ',rotangle
          if(nhis.lt.0)call hisarchive(ntpguest,gcmccount)
          if((nhis.gt.0).and.(prodcount.gt.0))call hisarchive
     &      (ntpguest,gcmccount) 
        endif
        if(nmols.ge.1.and.disp_count.ge.1)then
          if(mod(disp_count,100).eq.0)then
c update distance part of displacement on n00s
            if((accept_disp/dble(disp_count)).gt.disp_ratio)then 
              delrdisp=delrdisp*1.05d0
            else
              if(delrdisp.gt.delrmin)delrdisp=delrdisp*0.95d0
            endif
          elseif(mod(disp_count,50).eq.0)then
c update rotation part of displacement on n50s
            if((accept_disp/dble(disp_count)).gt.disp_ratio)then
              rotangle=rotangle*1.05d0
            else
              if(rotangle.gt.minangle)rotangle=rotangle*0.95d0
            endif
          endif
        endif
        if((nmols.ge.1).and.(tran_count.ge.1).and.
     &mod(tran_count,100).eq.0)then
c update distance moves on n00s
          if((accept_tran/dble(tran_count)).gt.tran_ratio)then 
            tran_delr=tran_delr*1.05d0
          else
            if(tran_delr.gt.delrmin)tran_delr=tran_delr*0.95d0
          endif
        endif
        if((nmols.ge.1).and.(rota_count.ge.1).and.
     &mod(rota_count,100).eq.0)then
c update rotation moves on n00s
          if((accept_rota/dble(rota_count)).gt.rota_ratio)then
            rota_rotangle=rota_rotangle*1.05d0
          else
            if(rota_rotangle.gt.minangle)
     &rota_rotangle=rota_rotangle*0.95d0
          endif
        endif
c       the following is added in case of initial conditions
        if (nmols.eq.0)then
          insert=.true.
          delete=.false.
          displace=.false.
          jump = .false.
          flex = .false.
          swap = .false.
          tran = .false.
          rota = .false.
          switch = .false.
        endif
c        if ((nmols.eq.0).or.(nmols.lt.mn_guest(iguest)))then
c          insert=.true.
c          delete=.false.
c          displace=.false.
c          jump = .false.
c          flex = .false.
c          swap = .false.
c          tran = .false.
c          rota = .false.
c          switch = .false.
c        endif
c        if(mx_guest(iguest).gt.0)then
c          if(nmols.gt.mx_guest(iguest))then
c            insert=.false.
c            delete=.true.
c            displace=.false.
c            jump = .false.
c            flex = .false.
c            swap = .false.
c            tran = .false.
c            rota = .false.
c            switch = .false.
c          endif
c        endif

c        ewald1en=0.d0
c        ewald2en=0.d0
c        vdwen=0.d0
c        chgsum_molorig=chgsum_mol
c        ckcsorig=ckcsum
c        ckssorig=ckssum
c        delE=0.d0
c***********************************************************************
c    
c       Insertion
c
c***********************************************************************
        if(insert)then
          call mc_insert
     &(idnode,imcon,keyfce,iguest,totatm,rcut,delr,ins_count,alpha,
     &drewd,ntpguest,ntpfram,ntpatm,volm,statvolm,kmax1,kmax2,kmax3,
     &epsq,dlrpot,maxvdw,newld,engunit,delrc,sumchg,maxmls,surftol,
     &overlap,accepted,temp,beta,accept_ins)
          insert=.false.
          

c***********************************************************************
c    
c             Deletion 
c
c***********************************************************************
        elseif(delete)then
          call mc_delete
     &(idnode,imcon,keyfce,iguest,totatm,rcut,delr,del_count,alpha,
     &drewd,ntpguest,ntpfram,ntpatm,volm,statvolm,kmax1,kmax2,kmax3,
     &epsq,dlrpot,maxvdw,newld,engunit,delrc,sumchg,maxmls,surftol,
     &overlap,accepted,temp,beta,accept_del)
          delete=.false.
          
            
c***********************************************************************
c    
c           Displacement 
c
c***********************************************************************

        elseif(displace)then
          call mc_displace
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,tran,
     &tran_count,tran_delr,rota,rota_count,rota_rotangle,disp_count,
     &delrdisp,rotangle,maxmls,kmax1,kmax2,kmax3,newld,alpha,rcut,delr,
     &drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,accepted,temp,
     &beta,delrc,ntpatm,maxvdw,accept_tran,accept_rota,accept_disp,
     &ntpguest,ntpfram)
          displace=.false.
          tran = .false.
          rota = .false.


c***********************************************************************
c    
c          Long jumps (eg insertion+deletion)
c
c***********************************************************************

        elseif(jump)then
          call mc_jump
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,
     &jump_count,jumpangle,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_jump,ntpfram,
     &ntpguest)
          jump=.false.


c***********************************************************************
c    
c         Framework flex 
c         PB - deprecated!
c
c***********************************************************************

        elseif(flex)then

          flex_count = flex_count + 1
          accepted = .false.
          delE_fwk = 0
          evdwg = 0
          ecoulg = 0

c FIXME(td): only works for a single guest
c Assume that energy of single guest is correct
          delE_fwk = energy(mol)+fwk_ener(curr_fwk)

c Select a new framework
c if sequential is requested move up to the maximum jump in a random
c direction
          new_fwk = curr_fwk
          if (l_fwk_seq)then
            do while ((new_fwk.eq.curr_fwk).or.(new_fwk.gt.n_fwk).or.
     &(new_fwk.lt.1))
              fwk_step_magnitude = 1+floor(duni(idnode)*
     &dble(fwk_step_max))
              if (duni(idnode).lt.0.5)then
                new_fwk = curr_fwk + fwk_step_magnitude
              else
                new_fwk = curr_fwk - fwk_step_magnitude
              endif
            enddo
c otherwise just pick any framework randomly
          else
            do while (new_fwk.eq.curr_fwk)
              new_fwk = 1+floor(duni(idnode)*dble(n_fwk))
            enddo
          endif
c store the state of the system
          state_cell = cell
          state_x = molxxx
          state_y = molyyy
          state_z = molzzz
          state_chg = atmchg
          state_vol = volm

c now start to switch out new configuration
          cell = fwk_cell(new_fwk,:)
c put in new positions for framework only
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
                  molxxx(k,indatm) = fwk_posx(new_fwk,indfwk)
                  molyyy(k,indatm) = fwk_posy(new_fwk,indfwk)
                  molzzz(k,indatm) = fwk_posz(new_fwk,indfwk)
                endif
              enddo
            enddo
          enddo
c switch charges for framework only
          do l = 1,ntpmls
            isguest = .false.
            do gstidx=1, ntpguest
              if (l.eq.locguest(gstidx)) isguest=.true.
            enddo
            if(.not.isguest) then
              do m = 1,numatoms(l) 
                atmchg(l,m) = fwk_chg(new_fwk, m)
              enddo
            endif
          enddo
c new volume
          volm = fwk_vol(new_fwk)
c calculate interations with new configuration
c without the condensing the energy breaks when wrapping on boundary
conditions
          ckcsnew = ckcsum
          ckssnew = ckssum

          call guest_exclude(ntpguest)
          call condense(totatm,ntpfram,ntpguest)
          call lrcorrect(imcon,keyfce,totatm,ntpatm,maxvdw,
     &rcut,volm,maxmls)
c Assume the single_point calculates the correct spenergy (evdwg+ecoulg)
          call single_point
     &(imcon,keyfce,alpha,drewd,rcut,delr,totatm,ntpfram,
     &ntpguest,ntpmls,volm,newld,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,
     &maxvdw,spenergy,vdwsum,ecoul,dlpeng,maxmls,surftol)

          delE_fwk = spenergy+fwk_ener(new_fwk)-delE_fwk

c energy acceptance criterion same as a displacement

          if(delE_fwk.lt.0.d0)then
            accepted=.true.
          else
            accepted=.false.
            rande=duni(idnode)
            call energy_eval
     &(estep,rande,statvolm,iguest,0,temp,beta,
     &.true.,.false.,.false.,.false.,accepted)
          endif

          if(accepted)then
c accept! Frameworks are fine as they are
            accept_flex = accept_flex + 1
            flx(iguest) = new_fwk
            energy(mol) = spenergy
            delE(mol) = delE_fwk
            cfgname(:) = fwk_name(new_fwk,:)
            curr_fwk = new_fwk
          else
c have to put everything back as it was
            flx(iguest) = -curr_fwk
            delE = 0.d0
            cell = state_cell
            molxxx = state_x
            molyyy = state_y
            molzzz = state_z
            atmchg = state_chg
            volm = state_vol
            ckcsum = ckcsnew
            ckssum = ckssnew
c rebuild all the atom lists
c            call guest_exclude(ntpguest)
            call condense(totatm,ntpfram,ntpguest)
c            call erfcgen(keyfce,alpha,rcut,drewd)
            call parlst(imcon,totatm,rcut,delr)
            call lrcorrect(imcon,keyfce,totatm,ntpatm,maxvdw,
     &rcut,volm,maxmls)
c             call single_point
c     &(imcon,keyfce,alpha,drewd,rcut,delr,totatm,ntpfram,
c     &ntpguest,ntpmls,volm,newld,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,
c     &maxvdw,spenergy,vdwsum,ecoul,dlpeng,maxmls,surftol)
          endif
          flex = .false.

        elseif(switch)then
c***********************************************************************
c    
c         Switch - only if multiple guests in sim
c         Takes one guest already in the framework (guesti) and
c         switches it's location with another guest already in the 
c         framework (guestj)
c
c***********************************************************************
c         first generate a maximum limit on the number of swaps
c         by determining the quantities of the different guest types
c         in the framework
          call mc_switch
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,
     &switch_count,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_switch,ntpfram,
     &ntpguest)

          switch = .false.
        elseif(swap)then
c***********************************************************************
c    
c         Swap - replace one guest in the simulation with another
c                of a different type
c
c***********************************************************************

          call mc_swap
     &(idnode,imcon,keyfce,iguest,jguest,totatm,volm,statvolm,
     &swap_count,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_swap,ntpfram,
     &ntpguest,rotangle)
          swap = .false.
        endif
c=========================================================================
c       once GCMC move is done, check if production is requested
c       if so, store averages, probability plots etc..
c=========================================================================
        if(production)then 
          if((laccsample).and.(accepted).or.(.not.laccsample))then
            if(n_fwk.gt.0)fwk_counts(curr_fwk)=fwk_counts(curr_fwk)+1
            prodcount=prodcount+1
            rollcount=rollcount+1
            chainstats(1) = dble(prodcount)
c           for the guest which was just perturbed and if the
c           move was accepted, update the list of framework 
c           adsorbed molecules and add to the average.
c           otherwise re-sum the old list.
            
            do i=1,ntpguest
              mol=locguest(i) 
              dmolecules=real(nummols(mol))
              molecules2=dmolecules*dmolecules
              energy2=energy(mol)*energy(mol)
              surfmol=real(surfacemols(mol))

              istat=1+16*(i-1)
c           avg <N>
              chainstats(istat+1) = chainstats(istat+1) + dmolecules
c           avg <E>
              chainstats(istat+2) = chainstats(istat+2) + energy(mol)
c           avg <EN>
              chainstats(istat+3) = chainstats(istat+3) +
     &energy(mol)*dmolecules 
c           avg <N^2>
              chainstats(istat+4) = chainstats(istat+4) + molecules2
c           avg <E^2>
              chainstats(istat+5) = chainstats(istat+5) + energy2
c           avg <NF>
              chainstats(istat+6) = chainstats(istat+6) + surfmol 
c           Sampling the Henry's coefficient (This requires an
c              energy calculation between the guest and framework
c              only.) Widom insertion is better for this purpos
c              as it more evenly samples configuration space.
              rollstat=9*(i-1)
            enddo
          endif
          if((mod(prodcount,nwindsteps).eq.0).or.(.not.lgchk))then
c         store averages for standard deviations
c         reset windows to zero
            avcount = avcount + 1
            avgwindow(1,avcount) = dble(rollcount)
            do i=1,ntpguest
              istat = 1+(i-1)*16
              rollstat = 1+(i-1)*9
c           get the averages for all the variables
              aN = chainstats(istat+1)/dble(rollcount) 
              E = chainstats(istat+2)/dble(rollcount)
              EN = chainstats(istat+3)/dble(rollcount)
              N2 = chainstats(istat+4)/dble(rollcount)
              E2 = chainstats(istat+5)/dble(rollcount)
              NF = chainstats(istat+6)/dble(rollcount)
c           store averages in separate array
              avgwindow(rollstat+1,avcount) = aN
              avgwindow(rollstat+2,avcount) = E
              avgwindow(rollstat+3,avcount) = EN
              avgwindow(rollstat+4,avcount) = N2
              avgwindow(rollstat+5,avcount) = E2
              avgwindow(rollstat+6,avcount) = NF
c           compute C_v and Q_st for the windowed averages
              avgwindow(rollstat+7,avcount) = 
     &calc_Qst(E, aN, N2, EN, temp)
              avgwindow(rollstat+8,avcount) = 
     &calc_Cv(E2, E, aN, N2, EN, temp)
c             reset chainstats 
              chainstats(istat+1)=0.d0
              chainstats(istat+2)=0.d0
              chainstats(istat+3)=0.d0
              chainstats(istat+4)=0.d0
              chainstats(istat+5)=0.d0
              chainstats(istat+6)=0.d0
            enddo 
            rollcount = 0
          endif
          if((lprob).and.(.not.lwidom))then
            call storeprob(ntpguest,rucell,ngrida,ngridb,ngridc)
            if((lprobeng(iguest)).and.(accepted))then
                call storeenprob(iguest,randchoice,rucell,ntpguest,
     &ngrida,ngridb,ngridc,guest_toten)
            endif
          endif
        endif

c       increment the numguests.out storage buffer

c        ibuff=ibuff+1
        if(lnumg.and.(mod(gcmccount,nnumg).eq.0))then
c          write(*,'(4e18.5)')ewld1eng/engunit,ewld2sum/engunit,
c     &ewld3sum/engunit,vdwsum/engunit
c          do i=1,natms
c            write(*,'(a3,3f15.5)')atmname(mol,i),newx(i),newy(i),newz(i)

c          enddo
          do i=1,ntpguest
            mol=locguest(i)
            write(400+i,"(i9,i7,2f20.6,7i4)")
     &        gcmccount,nummols(mol),energy(mol),
     &        delE(mol),ins(i),del(i),dis(i),jmp(i),flx(i),swp(i),swi(i)

          enddo
        endif
        do i=1,ntpguest 
          ins(i)=0
          del(i)=0
          dis(i)=0
          jmp(i)=0
          flx(i)=0
          swp(i)=0
          swi(i)=0
        enddo
      enddo
c*************************************************************************
c     END OF GCMC RUN
c*************************************************************************
      call timchk(0,timelp)

c     run statistics on uptake, energies locally
c     then add the sums globally for a global weighted
c     average
c      print *, "Ewald1 average", ewaldaverage/accept_ins
      if(prodcount.gt.0)then
        weight = chainstats(1)
        do i=1,ntpguest
          istat = 1+(i-1)*16
          vstat = 1+(i-1)*9
          avgn = 0.d0
          avge = 0.d0
          avgen = 0.d0
          avgn2 = 0.d0
          avge2 = 0.d0
          avgnf = 0.d0
          avgQst = 0.d0
          avgCv = 0.d0
          sumweight=0.d0
          do j = 1, avcount
            sweight = avgwindow(1,j)
            sumweight = sumweight+sweight
            avgn = avgn + sweight*avgwindow(vstat+1,j) 
            avge = avge + sweight*avgwindow(vstat+2,j)
            avgen = avgen + sweight*avgwindow(vstat+3,j) 
            avgn2 = avgn2 + sweight*avgwindow(vstat+4,j) 
            avge2 = avge2 + sweight*avgwindow(vstat+5,j) 
            avgnf = avgnf + sweight*avgwindow(vstat+6,j)
            avgQst = avgQst + sweight*avgwindow(vstat+7,j)
            avgCv = avgCv + sweight*avgwindow(vstat+8,j) 
          enddo
          avgn = avgn/sumweight
          avge = avge/sumweight
          avgen = avgen/sumweight
          avgn2 = avgn2/sumweight
          avge2 = avge2/sumweight
          avgnf = avgnf/sumweight
          avgQst = avgQst/sumweight
          avgCv = avgCv/sumweight
c         overwrite chainstats to communicate to other nodes
          chainstats(istat+1) = avgn
          chainstats(istat+2) = avge
          chainstats(istat+3) = avgen
          chainstats(istat+4) = avgn2
          chainstats(istat+5) = avge2
          chainstats(istat+6) = avgnf
          chainstats(istat+7) = avgQst
          chainstats(istat+8) = avgCv 

c         standard deviations.
          stdN = 0.d0
          stdE = 0.d0
          stdEN = 0.d0
          stdN2 = 0.d0
          stdE2 = 0.d0
          stdNF = 0.d0
          stdQst = 0.d0
          stdCv = 0.d0
          sumweight=0.d0
          do j = 1, avcount
            sweight = avgwindow(1,j)
            sumweight = sumweight+sweight
            stdN = stdN + sweight*(avgn - avgwindow(vstat+1,j))**2.d0
            stdE = stdE + sweight*(avge - avgwindow(vstat+2,j))**2.d0
            stdEN = stdEN + sweight*(avgen - avgwindow(vstat+3,j))**2.d0
            stdN2 = stdN2 + sweight*(avgn2 - avgwindow(vstat+4,j))**2.d0
            stdE2 = stdE2 + sweight*(avge2 - avgwindow(vstat+5,j))**2.d0
            stdNF = stdNF + sweight*(avgnf - avgwindow(vstat+6,j))**2.d0
            stdQst = stdQst + 
     & sweight*(avgQst - avgwindow(vstat+7,j))**2.d0
            stdCv = stdCv + 
     &sweight*(avgCv - avgwindow(vstat+8,j))**2.d0 
          enddo
          stdN = sqrt(stdN/sumweight)
          stdE = sqrt(stdE/sumweight)
          stdEN = sqrt(stdEN/sumweight)
          stdN2 = sqrt(stdN2/sumweight)
          stdE2 = sqrt(stdE2/sumweight)
          stdNF = sqrt(stdNF/sumweight)
          stdQst = sqrt(stdQst/sumweight)
          stdCv = sqrt(stdCv/sumweight)
          chainstats(istat+9) = stdN
          chainstats(istat+10) = stdE
          chainstats(istat+11) = stdEN
          chainstats(istat+12) = stdN2
          chainstats(istat+13) = stdE2
          chainstats(istat+14) = stdNF
          chainstats(istat+15) = stdQst
          chainstats(istat+16) = stdCv
        enddo
      endif

c     sum all variables for final presentation

c     first send all local statistics to the master node.
c     via csend([messagetag],[variable to send],[length],[destination],
c     [dummy var])
c     message tag for chainstats is even
      if(idnode.gt.0)then
        call csend(idnode*2+1,chainstats,1+ntpguest*16,0,1)
      endif


c     this is final file writing stuff.. probably should put this in a 
c     subroutine so it looks less messy in the main program. 
      if(idnode.eq.0)then

         nodeweight(1)=weight
         write(nrite,"(/,a100,/,3x,'Data reported from node ',
     &i3,/,a100,/)")repeat('*',100),0,repeat('*',100)
         do i=1,ntpguest
           istat=1+(i-1)*16
           mol=locguest(i)
           avgN = chainstats(istat+1)
           avgE = chainstats(istat+2)
           avgEN = chainstats(istat+3)
           avgN2 = chainstats(istat+4)
           avgE2 = chainstats(istat+5)
           avgNF = chainstats(istat+6)
           avgQst = chainstats(istat+7)
           avgCv = chainstats(istat+8)
           
           stdN = chainstats(istat+9) 
           stdE = chainstats(istat+10) 
           stdEN = chainstats(istat+11) 
           stdN2 = chainstats(istat+12) 
           stdE2 = chainstats(istat+13) 
           stdNF = chainstats(istat+14)
           stdQst = chainstats(istat+15)
           stdCv = chainstats(istat+16)

           write(nrite,"(5x,'guest ',i2,': ',40a,/)")i,
     &       (molnam(p,mol),p=1,40)
           if(.not.lwidom)then
             write(nrite,"(5x,a60,f20.9,/,5x,a60,f20.9,/,
     &         5x,a60,f20.9,/,5x,a60,f20.9,/,5x,a60,f20.9,/,
     &         5x,a60,i20,/,5x,a60,f20.9,/,5x,a60,f20.9,/,
     &         5x,a60,f20.9,/,5x,a60,f20.9,/,5x,a60,f20.9,/)")
     &         '<N>: ', avgN, 
     &         '<E>: ', avgE,
     &         '<E*N>: ', avgEN, 
     &         '<N*N>: ', avgN2,
     &         '<E*E>: ', avgE2,
     &         'Multiplier: ',prodcount,
     &         'Isosteric heat of adsorption (kcal/mol): ',avgQst,
     &         'Isosteric heat error: ', stdQst,
     &         'Heat capacity, Cv (kcal/mol/K): ', avgCv,
     &         'Heat capacity error: ', stdCv,
     &         '<surface adsorbed N>:', avgNF
           elseif(lwidom)then
             gcmccount=widcount
             avgH = chainstats(istat+7)
             write(nrite,"(5x,a60,E20.9,/)")
     &         "Henry's Constant (mol /kg /bar): ",
     &         avgH/avo/boltz/temp*1.d5
           endif
           nstat = (i-1)*9
           node_avg(1,(nstat+1):(nstat+8)) =
     &                    chainstats((istat+1):(istat+8))
           node_std(1,nstat+1:nstat+8) =
     &                   chainstats(istat+9:istat+16)
         enddo
         tw=tw+weight
         do i=1,mxnode-1
c        recieve all data from other nodes (via crecv)
c        prodcount used for weighting the mean and stdev
           statbuff(1:1+ntpguest*16) = 0.d0
           call crecv(i*2+1,statbuff,1+ntpguest*16,1)
           weight=statbuff(1)
           nodeweight(i+1)=weight
           tw=tw+weight
           prodcount=prodcount+int(weight)
           write(nrite,"(/,a100,/,3x,'Data reported from node ',i3,
     &/,a100,/)")repeat('*',100),i,repeat('*',100)
           do j=1,ntpguest
             mol=locguest(j)
             write(nrite,"(5x,'guest ',i2,': ',40a,/)")j,
     &       (molnam(p,mol),p=1,40)
             istat=1+(j-1)*16

             avgN = statbuff(istat+1)
             avgE = statbuff(istat+2)
             avgEN = statbuff(istat+3)
             avgN2 = statbuff(istat+4)
             avgE2 = statbuff(istat+5)
             avgNF = statbuff(istat+6)
             avgQst = statbuff(istat+7)
             avgCv = statbuff(istat+8)
            
             stdN = statbuff(istat+9)  
             stdE = statbuff(istat+10)
             stdEN = statbuff(istat+11)
             stdN2 = statbuff(istat+12)
             stdE2 = statbuff(istat+13)
             stdNF = statbuff(istat+14)
             stdQst = statbuff(istat+15)
             stdCv = statbuff(istat+16)

             !Q_st = calc_Qst(avgE, avgN, avgN2, avgEN, temp)
             !C_v = calc_cv(avgE2, avgE, avgN, avgN2, avgEN, temp)
             if((.not.lwidom).and.(.not.lwanglandau))then
               write(nrite,"(5x,a60,f20.9,/,5x,a60,f20.9,/,
     &           5x,a60,f20.9,/,5x,a60,f20.9,/,5x,a60,f20.9,/,
     &           5x,a60,i20,/,5x,a60,f20.9,/,5x,a60,f20.9,/,
     &           5x,a60,f20.9,/,5x,a60,f20.9,/,5x,a60,f20.9,/)")
     &           '<N>: ', avgN, 
     &           '<E>: ', avgE,
     &           '<E*N>: ', avgEN, 
     &           '<N*N>: ', avgN2,
     &           '<E*E>: ', avgE2,
     &           'Multiplier: ',prodcount,
     &           'Isosteric heat of adsorption (kcal/mol): ',avgQst,
     &           'Isosteric heat error: ', stdQst,
     &           'Heat capacity, Cv (kcal/mol/K): ', avgCv,
     &           'Heat capacity error: ', stdCv,
     &           '<surface adsorbed N>:', avgNF
             elseif(lwidom)then
                write(nrite,"(5x,a60,E20.9,/)")
     &            "Henry's Constant (mol /kg /bar): ", 
     &            avgH/avo/boltz/temp*1.d5
             endif
             nstat = (j-1)*9
             node_avg(i+1,nstat+1:nstat+8) = statbuff(istat+1:istat+8)
             node_std(i+1,nstat+1:nstat+8) = statbuff(istat+9:istat+16)
          enddo
        enddo
      endif

      call gisum(accept_ins,1,buffer)
      call gisum(ins_count,1,buffer)

      call gisum(accept_del,1,buffer)
      call gisum(del_count,1,buffer)

      call gisum(accept_disp,1,buffer)
      call gisum(disp_count,1,buffer)

      call gisum(accept_jump,1,buffer)
      call gisum(jump_count,1,buffer)

      call gisum(accept_flex,1,buffer)
      call gisum(flex_count,1,buffer)

      call gisum(accept_swap,1,buffer)
      call gisum(swap_count,1,buffer)

      call gisum(accept_tran,1,buffer)
      call gisum(tran_count,1,buffer)

      call gisum(accept_rota,1,buffer)
      call gisum(rota_count,1,buffer)
      
      call gisum(accept_switch,1,buffer)
      call gisum(switch_count,1,buffer)

      call gisum(gcmccount,1,buffer)

      if(n_fwk.gt.0)then 
        call gisum(fwk_counts,n_fwk,fwksumbuff)
      endif
c     write final probability cube files
      
      cprob=0
      cell=cell*angs2bohr
      ucell=ucell*angs2bohr
      if(lprob.and.prodcount.gt.0)then
        do i=1,ntpguest
          iprob=0
          if(lprobeng(i))then
              np=nprob(i)-2
          else
              np=nprob(i)
          endif
          do j=1,np
            cprob=cprob+1
            iprob=iprob+1
            if (.not.lwidom)then
              call gdsum3(grid,cprob,ntprob,gridsize,gridbuff)
              if(idnode.eq.0)call writeprob
     &  (i,cprob,iprob,ucell,ntpfram,gridsize,
     &  ngrida,ngridb,ngridc,prodcount,scell_factor)
            endif
          enddo
          if(lprobeng(i))then
              iprob=iprob+1
              cprob=cprob+1
              call gdsum3(grid,cprob,ntprob,gridsize,gridbuff)
c             add counters again for the tally grid.
              iprob=iprob+1
              cprob=cprob+1
              call gdsum3(grid,cprob,ntprob,gridsize,gridbuff)
              if(idnode.eq.0)call writeenprob
     &(i,cprob-1,ucell,ntpfram,gridsize,
     &ngrida,ngridb,ngridc,scell_factor)
          endif
        enddo
      endif
      
      if(idnode.eq.0)then
        write(nrite,"(/,a17,i9,a15,f13.3,a8)")
     &'time elapsed for ',gcmccount,' gcmc steps : ',timelp,' seconds'
        write(nrite,"(/,a30,i9)")
     &'total accepted steps : ',accept_ins+accept_del+accept_disp+
     &accept_jump+accept_flex+accept_swap+accept_tran+accept_rota
        if(ins_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'insertion ratio: ',dble(accept_ins)/dble(ins_count)
        if(del_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'deletion ratio: ',dble(accept_del)/dble(del_count)
        if(disp_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'displacement ratio: ',dble(accept_disp)/dble(disp_count)
        if(jump_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'jump ratio: ',dble(accept_jump)/dble(jump_count)
        if(flex_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'flex ratio: ',dble(accept_flex)/dble(flex_count)
        do jj=1,n_fwk
          write(nrite,"(/,6x,a26,i6,60a1,e13.6)")
     &'Framework population for: ',jj,(fwk_name(jj, kk),kk=1,60),
     & dble(fwk_counts(jj))/dble(prodcount)
        enddo
        if(swap_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'swap ratio: ',dble(accept_swap)/dble(swap_count)
        if(switch_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'switch ratio: ',dble(accept_switch)/dble(switch_count)
        if(tran_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'translation ratio: ',dble(accept_tran)/dble(tran_count)
        if(rota_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'rotation ratio: ',dble(accept_rota)/dble(rota_count)
        do i=1,ntpguest
          mol=locguest(i)

          avgN = 0.d0
          avgE = 0.d0
          avgEN = 0.d0
          avgN2 = 0.d0
          avgE2 = 0.d0
          avgNF = 0.d0
          avgQst = 0.d0
          avgCv = 0.d0
          stdN = 0.d0
          stdE = 0.d0
          stdEN = 0.d0
          stdN2 = 0.d0
          stdE2 = 0.d0
          stdNF = 0.d0
          stdQst = 0.d0
          stdCv = 0.d0
c         compute unions of averages and standard deviations
          call avunion(i,mxnode,avgN,avgE,avgEN,avgN2,
     &avgE2,avgNF,avgQst,avgCv)
c       isosteric heat of adsorption 
          !Q_st = calc_Qst(avgE, avgN, avgN2, avgEN,temp)
c       heat capacity
          !C_v = calc_Cv(avgE2, avgE, avgN, avgN2, avgEN,temp)
          call stdunion
     &(i,mxnode,stdN,stdE,stdEN,stdN2,stdE2,stdNF,stdQst,stdCv,
     &avgN,avgE,avgEN,avgN2,avgE2,avgNF,avgQst,avgCv)

          write(nrite,"(/,a100,/,'final stats for guest ',
     &          i2,3x,40a,/,a100,/)")
     &       repeat('=',100),
     &       i,(molnam(p,mol),p=1,40),
     &       repeat('=',100)

          if(.not.lwidom)then
            write(nrite,"(/,a36,f15.6,a5,f12.3,/,
     &        a36,f15.6,a5,f12.3,/,a36,f15.6,a5,f12.3,/,
     &        a36,f15.6,a5,f12.3,/,a36,f15.6,a5,f12.3,/)")
     &        '<N>: ',avgN, ' +/- ', stdN,
     &        '<E>: ',avgE, ' +/- ', stdE,
     &        '<E*N>: ',avgEN, ' +/- ', stdEN,
     &        '<N*N>: ',avgN2, ' +/- ', stdN2,
     &        '<E*E>: ',avgE2, ' +/- ', stdE2
            write(nrite,"(/,a60,f15.6,/,a60,f15.6)")
     &        'average number of guests: ',avgN,
     &        'standard error: ',stdN
            write(nrite,"(a60,f15.6,/,a60,f15.6)")
     &        'Isosteric heat of adsorption (kcal/mol): ',avgQst,
     &        'Isosteric heat error: ', stdQst
c           I added surface data here so that faps would read it in
c           as the [useless IMO] heat capacity of the guest. This
c           is a hack to force faps to report surface adsorbed data
c           in the C_v column of the faps results .csv file
            if(surftol.ge.0.d0)write(nrite,"(a60,f15.6,/,a60,f15.6)")
     &        'average surface adsorption: ',avgNF,
     &        'standard error: ',stdNF
            write(nrite,"(a60,f15.6,/,a60,f15.6,/,a60,i15,/)")
     &        'Heat capacity, Cv (kcal/mol/K): ', avgCv,
     &        'Heat capacity error: ', stdCv,
     &        'Total steps counted: ',int(tw)
          elseif(lwidom)then
            write(nrite,"(/a60,E15.6,/,a60,i15,/)")
     &        "Henry's Constant (mol /kg /bar): ", 
     &         avgH/avo/boltz/temp*1.d5,
     &        "Total steps counted: ",int(tw)
          endif
        enddo
      endif
      close(202)
      close(ncontrol)
      close(nconfig)
      close(nfield)
      do i=1,ntpguest
        if(lnumg)close(400+i)
        if(abs(nhis).gt.0)close(500+i)
      enddo
      if(idnode.eq.0)then
       close(nrite)
c       close(nang)
      endif
      call exitcomms()
      contains
      character*9 function month(date)
      implicit none
      character*2 date

      if(date.eq.'01')month='January'
      if(date.eq.'02')month='February'
      if(date.eq.'03')month='March'   
      if(date.eq.'04')month='April'   
      if(date.eq.'05')month='May'     
      if(date.eq.'06')month='June'
      if(date.eq.'07')month='July'
      if(date.eq.'08')month='August'
      if(date.eq.'09')month='September'
      if(date.eq.'10')month='October'
      if(date.eq.'11')month='November'
      if(date.eq.'12')month='December'

      return
      end function month
      subroutine avunion(iguest,mxnode,avgN,avgE,avgEN,avgN2,avgE2,
     &avgNF,avgQst,avgCv)

      implicit none
      real(8) avgE,avgN,avgEN,avgN2,avgE2,sumweight,weight
      real(8) avgNF,avgQst,avgCv
      integer iguest,node,mxnode,istat
      istat=(iguest-1)*9
      sumweight=0.d0
      do node=1,mxnode
        weight = nodeweight(node)
        sumweight = weight + sumweight
        avgN = avgN + weight*node_avg(node,istat+1)
        avgE = avgE + weight*node_avg(node,istat+2)
        avgEN = avgEN + weight*node_avg(node,istat+3)
        avgN2 = avgN2 + weight*node_avg(node,istat+4)
        avgE2 = avgE2 + weight*node_avg(node,istat+5)
        avgNF = avgNF + weight*node_avg(node,istat+6)
        avgQst = avgQst + weight*node_avg(node,istat+7)
        avgCv = avgCv + weight*node_avg(node,istat+8)
      enddo

      avgE = avgE/sumweight
      avgN = avgN/sumweight
      avgEN = avgEN/sumweight
      avgN2 = avgN2/sumweight
      avgE2 = avgE2/sumweight
      avgNF = avgNF/sumweight
      avgQst = avgQst/sumweight
      avgCv = avgCv/sumweight
      end subroutine avunion
       
      subroutine stdunion(iguest,mxnode,stdN,stdE,stdEN,stdN2,stdE2,
     &stdNF,stdQst,stdCv,avgN,avgE,avgEN,avgN2,avgE2,avgNF,
     &Q_st,C_v)

      implicit none
      real(8) stdN,stdE,stdEN,stdN2,stdE2,stdQst,stdCv
      real(8) weight,sumweight, stdNF
      real(8) avgN,avgE,avgEN,avgN2,avgE2,Q_st,C_v,avgNF
      integer mxnode,iguest,istat,node
      
      istat = (iguest-1)*9
      sumweight=0.d0
      do node = 1,mxnode 
        weight = nodeweight(node)
        sumweight = weight+sumweight

        stdN = stdN + weight*
     &(node_std(node,istat+1)**2+node_avg(node,istat+1)**2)
        stdE = stdE + weight*
     &(node_std(node,istat+2)**2+node_avg(node,istat+2)**2)
        stdEN = stdEN + weight*
     &(node_std(node,istat+3)**2+node_avg(node,istat+3)**2)
        stdN2 = stdN2 + weight*
     &(node_std(node,istat+4)**2+node_avg(node,istat+4)**2)
        stdE2 = stdE2 + weight*
     &(node_std(node,istat+5)**2+node_avg(node,istat+5)**2)
        stdNF = stdNF + weight*
     &(node_std(node,istat+6)**2+node_avg(node,istat+6)**2)
        stdQst = stdQst + weight*
     &(node_std(node,istat+7)**2+node_avg(node,istat+7)**2)
        stdCv = stdCv + weight*
     &(node_std(node,istat+8)**2+node_avg(node,istat+8)**2)
      enddo
      stdN = sqrt((stdN/sumweight) - avgN**2)
      stdE = sqrt((stdE/sumweight) - avgE**2)
      stdEN = sqrt((stdEN/sumweight) - avgEN**2)
      stdN2 = sqrt((stdN2/sumweight) - avgN2**2)
      stdE2 = sqrt((stdE2/sumweight) - avgE2**2)
      stdNF = sqrt((stdNF/sumweight) - avgNF**2)
      stdQst = sqrt((stdQst/sumweight) - Q_st**2)
      stdCv = sqrt((stdCv/sumweight) - C_v**2)
      end subroutine stdunion
      subroutine hisarchive(ntpguest,gcmccount)
c*****************************************************************************
c
c     writes an xyz file of the current geometries of the guest
c
c*****************************************************************************
      implicit none

      integer i,imol,nmols,natms,gcmccount,iatm,ntpguest

      do k=1,ntpguest
        imol=locguest(k)
        nmols=nummols(imol)
        natms=numatoms(imol)
        write(500+k,"(i6,/,'step ',i7)")nmols*natms,gcmccount

        iatm=0
        do i=1,nmols
           do j=1,natms
            iatm=iatm+1
            write(500+k,'(2x,a1,3f16.6)')atmname(imol,j),
     &  molxxx(imol,iatm),molyyy(imol,iatm),molzzz(imol,iatm)
          enddo
        enddo
      enddo

      end subroutine hisarchive

      subroutine mc_insert
     &(idnode,imcon,keyfce,iguest,totatm,rcut,delr,ins_count,alpha,
     &drewd,ntpguest,ntpfram,ntpatm,volm,statvolm,kmax1,kmax2,kmax3,
     &epsq,dlrpot,maxvdw,newld,engunit,delrc,sumchg,maxmls,surftol,
     &overlap,accepted,temp,beta,accept_ins)
c*******************************************************************************
c
c     keeps track of all the associated arrays and calls the 'insertion' 
c     subroutine. Acceptance/rejection criteria based on Metropolis
c     importance sampling coupled to a grand canonical ensemble.
c
c*******************************************************************************
      implicit none
      logical lnewsurf,loverlap,accepted
      integer iguest,idnode,imcon,natms,totatm,mol,nmol
      integer ins_count,keyfce,ntpguest,kmax1,kmax2,kmax3
      integer ntpatm,maxvdw,maxmls,newld,accept_ins
      integer ntpfram,randchoice
      real(8) rcut,delr,estep,alpha,drewd,volm
      real(8) epsq,dlrpot,engunit,delrc,sumchg,chgtmp,engsictmp
      real(8) surftol,overlap,rande,gpress,statvolm,temp,beta
      mol=locguest(iguest)
      nmol=nummols(mol)
      natms=numatoms(mol)
      lnewsurf = .false.
      engsicorig=engsic
      ins(iguest)=1
      ins_count=ins_count+1
      call random_ins(idnode,natms,iguest,rcut,delr)
      estep = 0.d0
      call insertion
     & (imcon,iguest,keyfce,alpha,rcut,delr,drewd,totatm,
     & volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     & engunit,delrc,estep,sumchg,chgtmp,engsictmp,maxmls,
     & loverlap,lnewsurf,surftol,overlap,newld)
      accepted=.false.

      if (.not.loverlap)then
        gpress=gstfuga(iguest)
        rande=duni(idnode)
        call energy_eval
     &(estep,rande,statvolm,iguest,0,temp,beta,
     &.false.,.true.,.false.,.false.,accepted)
      endif
c     DEBUG
c      accepted=.true.
c     END DEBUG
      if(accepted)then
        accept_ins=accept_ins+1
        randchoice=0
        call accept_move
     &(iguest,.true.,.false.,.false.,
     &lnewsurf,delrc,totatm,randchoice,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
      else
        call reject_move
     &(iguest,0,.true.,.false.,.false.,.false.)
      endif
      return
      end subroutine mc_insert

      subroutine mc_delete
     &(idnode,imcon,keyfce,iguest,totatm,rcut,delr,del_count,alpha,
     &drewd,ntpguest,ntpfram,ntpatm,volm,statvolm,kmax1,kmax2,kmax3,
     &epsq,dlrpot,maxvdw,newld,engunit,delrc,sumchg,maxmls,surftol,
     &overlap,accepted,temp,beta,accept_del)
c*******************************************************************************
c
c     keeps track of all the associated arrays and calls the 'deletion' 
c     subroutine. Acceptance/rejection criteria based on Metropolis
c     importance sampling coupled to a grand canonical ensemble.
c
c*******************************************************************************
      implicit none
      logical linitsurf,accepted
      integer iguest,idnode,imcon,totatm,mol,nmol
      integer del_count,keyfce,ntpguest,kmax1,kmax2,kmax3
      integer ntpatm,maxvdw,maxmls,newld,accept_del
      integer ntpfram,randchoice
      real(8) rcut,delr,estep,alpha,drewd,volm
      real(8) epsq,dlrpot,engunit,delrc,sumchg,chgtmp,engsictmp
      real(8) surftol,overlap,rande,gpress,statvolm,temp,beta
      mol=locguest(iguest)
      nmol=nummols(mol)

      engsicorig=engsic
      linitsurf = .false.
      del(iguest)=1
      del_count=del_count+1 
      randchoice=floor(duni(idnode)*nmol)+1
      estep = 0.d0
      call deletion 
     &(imcon,keyfce,iguest,randchoice,alpha,rcut,delr,drewd,
     &maxmls,totatm,volm,kmax1,kmax2,kmax3,epsq,dlrpot,
     &ntpatm,maxvdw,engunit,delrc,estep,linitsurf,surftol,sumchg,
     &engsictmp,chgtmp,overlap,newld)

      gpress=gstfuga(iguest)
      accepted=.false.

      rande=duni(idnode)
      call energy_eval
     &(-estep,rande,statvolm,iguest,0,temp,beta,
     &.false.,.false.,.true.,.false.,accepted)
         
c     the following occurs if the move is accepted.
      if(accepted)then
        accept_del=accept_del+1
        call accept_move
     &(iguest,.false.,.true.,.false.,
     &linitsurf,delrc,totatm,randchoice,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
      else
        call reject_move
     &(iguest,0,.false.,.true.,.false.,.false.)
      endif
      return
      end subroutine mc_delete

      subroutine mc_displace
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,tran,
     &tran_count,tran_delr,rota,rota_count,rota_rotangle,disp_count,
     &delrdisp,rotangle,maxmls,kmax1,kmax2,kmax3,newld,alpha,rcut,delr,
     &drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,accepted,temp,
     &beta,delrc,ntpatm,maxvdw,accept_tran,accept_rota,accept_disp,
     &ntpguest,ntpfram)
c*******************************************************************************
c
c     Keeps track of all the associated arrays and calls the 
c     'displace_guest' subroutine. 
c     The underlying mechanics of this code is a 'deletion/insertion'
c     move, where the molecule is only perturbed as far as delrdisp and
c     rotangle dictates.
c     Acceptance/rejection criteria based on Metropolis
c     importance sampling coupled to a canonical ensemble.
c
c*******************************************************************************
      implicit none
      logical tran,rota,loverlap,linitsurf,lnewsurf,accepted
      integer iguest,keyfce,imcon,idnode,tran_count,rota_count
      integer disp_count,mol,ik,newld,maxmls,totatm,maxvdw
      integer kmax1,kmax2,kmax3,randchoice,nmol,ntpatm
      integer accept_tran,accept_rota,accept_disp
      integer ntpguest,ntpfram
      real(8) a,b,c,q1,q2,q3,q4,tran_delr,rota_rotangle,epsq
      real(8) delrdisp,rotangle,volm,alpha,rcut,delr,drewd
      real(8) engunit,overlap,surftol,dlrpot,sumchg,estep,rande
      real(8) statvolm,temp,beta,delrc,engsictmp,chgtmp,guest_toten

      mol=locguest(iguest)
      nmol=nummols(mol)

      a=0.d0;b=0.d0;c=0.d0;q1=0.d0;q2=0.d0;q3=0.d0;q4=0.d0
      if(tran)then
        dis(iguest)=2
        tran_count = tran_count + 1
        call random_disp(idnode,tran_delr,a,b,c)
      elseif(rota)then
        dis(iguest)=3
        rota_count = rota_count + 1
c        a=0.5d0
c        b=0.5d0
c        c=0.5d0
        call random_rot(idnode,rota_rotangle,q1,q2,q3,q4)
      else
        dis(iguest)=1
        disp_count=disp_count+1
        call random_disp(idnode,delrdisp,a,b,c)
        call random_rot(idnode,rotangle,q1,q2,q3,q4)
      endif
      do ik=1,newld
        ckcsorig(mol,ik)=ckcsum(mol,ik) 
        ckssorig(mol,ik)=ckssum(mol,ik) 
        ckcsorig(maxmls+1,ik)=ckcsum(maxmls+1,ik) 
        ckssorig(maxmls+1,ik)=ckssum(maxmls+1,ik)
      enddo
c     choose a molecule from the list
      randchoice=floor(duni(idnode)*nmol)+1

      call displace_guest
     &(imcon,alpha,rcut,delr,drewd,totatm,newld,maxmls,volm,kmax1,
     &kmax2,kmax3,epsq,engunit,overlap,surftol,linitsurf,lnewsurf,
     &loverlap,iguest,randchoice,dlrpot,sumchg,a,b,c,q1,q2,q3,q4,estep)
      accepted=.false.
      if (.not.loverlap)then
        if(estep.lt.0.d0)then
          accepted=.true.
        else
          rande=duni(idnode)
          call energy_eval
     &(estep,rande,statvolm,iguest,0,temp,beta,
     &.true.,.false.,.false.,.false.,accepted)
        endif
      endif
c     DEBUG
c      accepted=.true.
c      accepted=.false.
c     END DEBUG
      if(accepted)then
        call accept_move
     &(iguest,.false.,.false.,.true.,
     &lnewsurf,delrc,totatm,randchoice,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
c       check if the energy probability grid is requested,
c       then one has to compute the standalone energy
c       at the new position (requires ewald1 calc)
        if(lprobeng(iguest))then
          call gstlrcorrect(imcon,iguest,keyfce,
     &                        ntpatm,maxvdw,delrc,
     &                        volm,maxmls,.true.) 
          guest_toten=estep-delrc/engunit
        endif            
c       tally surface molecules
        if((linitsurf).and.(.not.lnewsurf))then
          surfacemols(mol) = surfacemols(mol) - 1
          if(surfacemols(mol).lt.0)surfacemols(mol) = 0
        elseif((.not.linitsurf).and.(lnewsurf))then
          surfacemols(mol) = surfacemols(mol) + 1
        endif
        if(tran)then
          accept_tran = accept_tran + 1
        elseif(rota)then
          accept_rota = accept_rota + 1
        else
          accept_disp=accept_disp+1
        endif
      else
        do ik=1,newld
          ckcsum(mol,ik)=ckcsorig(mol,ik) 
          ckssum(mol,ik)=ckssorig(mol,ik) 
          ckcsum(maxmls+1,ik)=ckcsorig(maxmls+1,ik) 
          ckssum(maxmls+1,ik)=ckssorig(maxmls+1,ik)
c          ckcsnew(mol,ik)=0.d0
c          ckssnew(mol,ik)=0.d0
c          ckcsnew(maxmls+1,ik)=0.d0 
c          ckssnew(maxmls+1,ik)=0.d0
        enddo
        call reject_move
     &(iguest,0,.false.,.false.,.true.,.false.)
      endif
      return
      end subroutine mc_displace

      subroutine mc_jump
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,
     &jump_count,jumpangle,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_jump,ntpfram,
     &ntpguest)
c*******************************************************************************
c
c     Keeps track of all the associated arrays and calls the 
c     'displace_guest' subroutine. 
c     The underlying mechanics of this code is a 'deletion/insertion'
c     move, where the molecule is randomly perturbed to another point  
c     in the simulation cell.
c     Acceptance/rejection criteria based on Metropolis
c     importance sampling coupled to a canonical ensemble.
c
c*******************************************************************************
      implicit none
      logical accepted,linitsurf,lnewsurf,loverlap
      integer idnode,imcon,keyfce,iguest,totatm,jump_count,maxmls,kmax1
      integer kmax2,kmax3,newld,ntpatm,maxvdw,accept_jump,randchoice
      integer nmol,natms,mol,ik,ntpfram,ntpguest
      real(8) volm,statvolm,jumpangle,alpha,rcut,delr,drewd,epsq,engunit
      real(8) overlap,surftol,dlrpot,sumchg,temp,beta,delrc
      real(8) a,b,c,q1,q2,q3,q4,estep,gpress,rande,engsictmp,chgtmp
      real(8) xc,yc,zc,comx,comy,comz
      mol=locguest(iguest)
      nmol=nummols(mol)
      natms=numatoms(mol)

      jmp(iguest)=1
      jump_count=jump_count+1
c     choose a molecule from the list
      randchoice=floor(duni(idnode)*nmol)+1
c     find which index the molecule "randchoice" is
      do ik=1,newld
        ckcsorig(mol,ik)=ckcsum(mol,ik) 
        ckssorig(mol,ik)=ckssum(mol,ik) 
        ckcsorig(maxmls+1,ik)=ckcsum(maxmls+1,ik) 
        ckssorig(maxmls+1,ik)=ckssum(maxmls+1,ik)
      enddo
      call get_guest(iguest,randchoice,mol,natms,nmol)

      call com(natms,mol,newx,newy,newz,comx,comy,comz)
      call random_jump(idnode,xc,yc,zc,jumpangle)
      call random_rot(idnode,jumpangle,q1,q2,q3,q4)
c      call random_disp(idnode,10.d0,a,b,c)
c     have to shift the a,b,c values by the com of the current guest
c     so the displace_guest routine can be used here..
      a = xc - comx
      b = yc - comy
      c = zc - comz

c     displace is a translation + rotation of a guest based on
c     cartesian a,b,c and quaternion q1,q2,q3,q4.
c     Jump is technically a deletion + insertion so the above
c     manipulation of the randomly generated a,b,c values
c     converts it to a translation.
      call displace_guest
     &(imcon,alpha,rcut,delr,drewd,totatm,newld,maxmls,volm,kmax1,
     &kmax2,kmax3,epsq,engunit,overlap,surftol,linitsurf,lnewsurf,
     &loverlap,iguest,randchoice,dlrpot,sumchg,a,b,c,q1,q2,q3,q4,estep)

      accepted=.false.
      if (.not.loverlap)then
        gpress=gstfuga(iguest)
        if(estep.lt.0.d0)then
          accepted=.true.
        else
          accepted=.false.
          rande=duni(idnode)
          call energy_eval
     &(estep,rande,statvolm,iguest,0,temp,beta,
     &.true.,.false.,.false.,.false.,accepted)
        endif
      endif

      if(accepted)then
        accept_jump=accept_jump+1
        call accept_move
     &(iguest,.false.,.false.,.true.,
     &lnewsurf,delrc,totatm,randchoice,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
c       tally surface molecules
        if((linitsurf).and.(.not.lnewsurf))then
          surfacemols(mol) = surfacemols(mol) - 1
          if (surfacemols(mol).lt.0)surfacemols(mol) = 0
        elseif((.not.linitsurf).and.(lnewsurf))then
          surfacemols(mol) = surfacemols(mol) + 1
        endif
      else
        do ik=1,newld
          ckcsum(mol,ik)=ckcsorig(mol,ik) 
          ckssum(mol,ik)=ckssorig(mol,ik) 
          ckcsum(maxmls+1,ik)=ckcsorig(maxmls+1,ik) 
          ckssum(maxmls+1,ik)=ckssorig(maxmls+1,ik)
        enddo
        call reject_move
     &(iguest,0,.false.,.false.,.true.,.false.)
      endif
      return
      end subroutine mc_jump

      subroutine mc_switch
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,
     &switch_count,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_switch,ntpfram,
     &ntpguest)
c*******************************************************************************
c
c     Keeps track of all the associated arrays and calls the 
c     The underlying mechanics of this code is a 
c     'deletion / displacement / insertion'
c     move, where two molecules of different types in the simulation cell 
c     switch locations.
c     Acceptance/rejection criteria based on Metropolis
c     importance sampling coupled to a canonical ensemble.
c
c     *** Currently only one switch is performed, but one could possibly
c         include multi switches as long as (detailed) balance isn't
c         violated. ***
c*******************************************************************************
      implicit none
      logical accepted,loverlap,loverallap,linitsurf,linitsurfj
      logical lnewsurf,lnewsurfj
      integer switch_count,nswitchgst,ntpguest,i,j,mol,nmols,jguest
      integer ichoice,jchoice,imol,jmol,ik,kk,iatm,idnode,imcon,keyfce
      integer iguest,totatm,maxmls,kmax1,kmax2,kmax3,newld
      integer ntpatm,maxvdw,accept_switch,ntpfram,natms
      real(8) icomx,icomy,icomz,jcomx,jcomy,jcomz,estep
      real(8) estepi,estepj,volm,statvolm,alpha,rcut,delr,drewd,epsq
      real(8) engunit,overlap,surftol,dlrpot,sumchg,temp,beta,delrc
      real(8) engsictmp,chgtmp,rande 
      nswitchgst = 0
      do i=1,ntpguest
        mol=locguest(i)
        nmols=nummols(mol)
c       generate array of molecules to randomly select
        do j=1,nmols
          switch_mols(i,j) = j
        enddo
        switch_mol_count(i) = nmols
        if((nmols.gt.0).and.(i.ne.iguest))then
          nswitchgst=nswitchgst+1
          switch_chosen_guest(nswitchgst)=i
        endif
      enddo
c     return if one can't make a switch move
      if(nswitchgst.lt.2)return
      switch_count=switch_count+1
      engsicorig = engsic
c     loverallap is true if one or both of the guests 
c     overlap with other atoms in their new configurations.
      loverallap=.false.
c     random choice of the second guest to switch
      jguest = floor(duni(idnode)*nswitchgst)+1
      jguest = switch_chosen_guest(jguest)
c     molecule choice to switch on each
      swi(iguest)=swi(iguest)+1
      ichoice=floor(duni(idnode)*switch_mol_count(iguest))+1
      ichoice = switch_mols(iguest,ichoice)
      
      swi(jguest)=swi(jguest)+1
      jchoice=floor(duni(idnode)*switch_mol_count(jguest))+1
      jchoice= switch_mols(jguest,jchoice)
c     store original framework configuration if the move is rejected
      do i=1,maxmls
        origenergy(i) = energy(i)
        origsurfmols(i) = surfacemols(i)
      enddo
      jmol=locguest(jguest)
      imol=locguest(iguest)
c     back up original arrays in case move is rejected
      do kk=1,totatm
        origmolxxx(imol,kk) = molxxx(imol,kk)
        origmolyyy(imol,kk) = molyyy(imol,kk)
        origmolzzz(imol,kk) = molzzz(imol,kk)
        origmolxxx(jmol,kk) = molxxx(jmol,kk)
        origmolyyy(jmol,kk) = molyyy(jmol,kk)
        origmolzzz(jmol,kk) = molzzz(jmol,kk)
      enddo
      do ik=1,newld
        ckcsorig(imol,ik)=ckcsum(imol,ik) 
        ckssorig(imol,ik)=ckssum(imol,ik) 
        ckcsorig(jmol,ik)=ckcsum(jmol,ik) 
        ckssorig(jmol,ik)=ckssum(jmol,ik) 
        ckcsorig(maxmls+1,ik)=ckcsum(maxmls+1,ik) 
        ckssorig(maxmls+1,ik)=ckssum(maxmls+1,ik)
c        ckcsnew(imol,ik)=0.d0
c        ckssnew(imol,ik)=0.d0
c        ckcsnew(jmol,ik)=0.d0
c        ckssnew(jmol,ik)=0.d0
c        ckcsnew(maxmls+1,ik)=0.d0 
c        ckssnew(maxmls+1,ik)=0.d0
      enddo
c     keeping track of the guest orientations by
c     re-populating the 'template' configurations for
c     guestx,guesty,guestz
      call get_guest(jguest,jchoice,jmol,natms,nmols)
      call com(natms,jmol,newx,newy,newz,jcomx,jcomy,jcomz)
      do iatm=1,natms
        guestx(jguest,iatm)=newx(iatm) - jcomx
        guesty(jguest,iatm)=newy(iatm) - jcomy
        guestz(jguest,iatm)=newz(iatm) - jcomz
      enddo
      call get_guest(iguest,ichoice,imol,natms,nmols)
      call com(natms,imol,newx,newy,newz,icomx,icomy,icomz)
      do iatm=1,natms
        guestx(iguest,iatm)=newx(iatm) - icomx
        guesty(iguest,iatm)=newy(iatm) - icomy
        guestz(iguest,iatm)=newz(iatm) - icomz
      enddo
c************************************************************************
c       START SWITCH OF GUESTI AND GUESTJ
c       Delete it, since shifting it to the new position would
c       create infinite energies.
c************************************************************************
      estep = 0.d0
      call deletion 
     &(imcon,keyfce,iguest,ichoice,alpha,rcut,delr,drewd,maxmls,
     &totatm,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,delrc,estep,linitsurf,surftol,sumchg,engsictmp,chgtmp,
     &overlap,newld)
c     have to default accept move so the energy arrays are updated
c     and ichoice from iguest is actually deleted from the system
c     so that jchoice from jguest can be inserted there.
      call accept_move
     &(iguest,.false.,.true.,.false.,linitsurf,delrc,totatm,ichoice,
     &ntpfram,ntpguest,maxmls,sumchg,engsictmp,chgtmp,newld)
      estepi=-estep
      call get_guest(jguest,jchoice,jmol,natms,nmols)
      estep=0.d0
      call deletion 
     &(imcon,keyfce,jguest,jchoice,alpha,rcut,delr,drewd,maxmls,
     &totatm,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,delrc,estep,linitsurfj,surftol,sumchg,engsictmp,chgtmp,
     &overlap,newld)
      call accept_move
     &(jguest,.false.,.true.,.false.,linitsurfj,delrc,totatm,jchoice,
     &ntpfram,ntpguest,maxmls,sumchg,engsictmp,chgtmp,newld)
      estepj=-estep

c     now insert the guests in their new positions (jguest in iguests
c     position and vise versa)
      do iatm=1,natms
        newx(iatm) = guestx(jguest,iatm) + icomx
        newy(iatm) = guesty(jguest,iatm) + icomy
        newz(iatm) = guestz(jguest,iatm) + icomz
      enddo
      estep=0.d0
      call insertion
     & (imcon,jguest,keyfce,alpha,rcut,delr,drewd,totatm,
     & volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     & engunit,delrc,estep,sumchg,chgtmp,engsictmp,maxmls,
     & loverlap,lnewsurfj,surftol,overlap,newld)
      if(loverlap)loverallap=.true.
      estepj=estepj+estep
      call accept_move
     &(jguest,.true.,.false.,.false.,lnewsurfj,delrc,totatm,jchoice,
     &ntpfram,ntpguest,maxmls,sumchg,engsictmp,chgtmp,newld)

      call get_guest(iguest,ichoice,imol,natms,nmols)
      do iatm=1,natms
        newx(iatm) = guestx(iguest,iatm) + jcomx
        newy(iatm) = guesty(iguest,iatm) + jcomy
        newz(iatm) = guestz(iguest,iatm) + jcomz
      enddo
      estep=0.d0
      call insertion
     & (imcon,iguest,keyfce,alpha,rcut,delr,drewd,totatm,
     & volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     & engunit,delrc,estep,sumchg,chgtmp,engsictmp,maxmls,
     & loverlap,lnewsurf,surftol,overlap,newld)
      if(loverlap)loverallap=.true.
      estepi=estepi+estep
      call accept_move
     &(iguest,.true.,.false.,.false.,lnewsurf,delrc,totatm,ichoice,
     &ntpfram,ntpguest,maxmls,sumchg,engsictmp,chgtmp,newld)
c************************************************************************
c           END OF SWITCH
c************************************************************************
c     perform energy evaluation
      if((estepi+estepj).lt.0.d0)then
        accepted=.true.
      else
        accepted=.false.
        rande=duni(idnode)
c       assume this move is governed by the canonical ensemble
        call energy_eval
     &(estepi+estepj,rande,statvolm,iguest,jguest,temp,beta,
     &.true.,.false.,.false.,.false.,accepted)
      endif
      if(loverallap)accepted=.false.
c     DEBUG
c      accepted=.false.
c     END DEBUG
      if(accepted)then
        accept_switch=accept_switch+1
        call condense(totatm,ntpfram,ntpguest)
      else
        do ik=1,maxmls
          delE(ik)=0.d0
          energy(ik)=origenergy(ik)
          surfacemols(ik)=origsurfmols(ik)
        enddo
c       restore original framework if move is rejected
        do ik=1,totatm 
          molxxx(imol,ik)=origmolxxx(imol,ik)
          molyyy(imol,ik)=origmolyyy(imol,ik)
          molzzz(imol,ik)=origmolzzz(imol,ik)
          molxxx(jmol,ik)=origmolxxx(jmol,ik)
          molyyy(jmol,ik)=origmolyyy(jmol,ik)
          molzzz(jmol,ik)=origmolzzz(jmol,ik)
        enddo
c       restore original ewald1 sums if step is rejected
        do ik=1,newld
          ckcsum(imol,ik)=ckcsorig(imol,ik) 
          ckssum(imol,ik)=ckssorig(imol,ik) 
          ckcsum(jmol,ik)=ckcsorig(jmol,ik) 
          ckssum(jmol,ik)=ckssorig(jmol,ik) 
          ckcsum(maxmls+1,ik)=ckcsorig(maxmls+1,ik) 
          ckssum(maxmls+1,ik)=ckssorig(maxmls+1,ik)
c          ckcsnew(imol,ik)=0.d0
c          ckssnew(imol,ik)=0.d0
c          ckcsnew(jmol,ik)=0.d0
c          ckssnew(jmol,ik)=0.d0
c          ckcsnew(maxmls+1,ik)=0.d0 
c          ckssnew(maxmls+1,ik)=0.d0
        enddo
        elrc=origelrc
        elrc_mol=origelrc_mol
        engsic=engsicorig
c       restore original surfacemols if step is rejected
        call condense(totatm,ntpfram,ntpguest)
      endif
      return
      end subroutine mc_switch

      subroutine mc_swap
     &(idnode,imcon,keyfce,iguest,jguest,totatm,volm,statvolm,
     &swap_count,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_swap,ntpfram,
     &ntpguest,rotangle)
c*******************************************************************************
c
c     Keeps track of all the associated arrays and calls the 
c     'deletion' and 'insertion' subroutines. 
c     The underlying mechanics of this code is a 'deletion/insertion'
c     move, where a molecule in the simulation cell is swapped with a
c     molecule with a different identity.
c     Acceptance/rejection criteria based on Metropolis
c     importance sampling coupled to a grand canonical ensemble.
c
c*******************************************************************************
      implicit none
      logical accepted,linitsurf,lnewsurf,loverlap
      integer idnode,imcon,keyfce,iguest,jguest,totatm,swap_count
      integer maxmls,kmax1,kmax2,kmax3,origtotatm,ik,imol,jmol,iatm
      integer newld,ntpatm,maxvdw,accept_swap,randchoice
      integer nmols,natms,mol,ntpfram,ntpguest,ichoice,jchoice,j
      integer jnatms,jnmols
      real(8) volm,statvolm,alpha,rcut,delr,drewd,epsq,engunit
      real(8) overlap,surftol,dlrpot,sumchg,temp,beta,delrc
      real(8) a,b,c,q1,q2,q3,q4,estepi,estepj,gpress,rande,estep
      real(8) comx,comy,comz,engsictmp,chgtmp,rotangle

      swap_count = swap_count+1
      origtotatm = totatm 

      imol=locguest(iguest)
      natms=numatoms(imol)
      nmols=nummols(imol)
c     store original ewald1 sums in case the move is rejected
      do ik=1,maxmls
        chgsum_molorig(ik)=chgsum_mol(ik)
        origsurfmols(ik)=surfacemols(ik)
        origenergy(ik)=energy(ik)
      enddo
      engsicorig=engsic
      origelrc = elrc
      origelrc_mol=elrc_mol
      do j=1,iguest-1
        switch_chosen_guest(j) = j
      enddo 
      do j=iguest,ntpguest
        switch_chosen_guest(j)=j+1
      enddo
      ichoice=floor(duni(idnode)*nmols)+1
      call get_guest(iguest,ichoice,imol,natms,nmols)
      call com(natms,imol,newx,newy,newz,comx,comy,comz)
c     chose a second guest to swap with
      jguest = floor(duni(idnode) * (ntpguest-1)) + 1
      jguest = switch_chosen_guest(jguest)
      jmol=locguest(jguest)
c     store original framework configuration if the move is rejected
      do ik=1,totatm
        origmolxxx(imol,ik)=molxxx(imol,ik)
        origmolyyy(imol,ik)=molyyy(imol,ik)
        origmolzzz(imol,ik)=molzzz(imol,ik)
        origmolxxx(jmol,ik)=molxxx(jmol,ik)
        origmolyyy(jmol,ik)=molyyy(jmol,ik)
        origmolzzz(jmol,ik)=molzzz(jmol,ik)
      enddo
      do ik=1,newld
        ckcsorig(imol,ik)=ckcsum(imol,ik) 
        ckssorig(imol,ik)=ckssum(imol,ik) 
        ckcsorig(jmol,ik)=ckcsum(jmol,ik) 
        ckssorig(jmol,ik)=ckssum(jmol,ik) 
        ckcsorig(maxmls+1,ik)=ckcsum(maxmls+1,ik) 
        ckssorig(maxmls+1,ik)=ckssum(maxmls+1,ik)
c        ckcsnew(imol,ik)=0.d0
c        ckssnew(imol,ik)=0.d0
c        ckcsnew(jmol,ik)=0.d0
c        ckssnew(jmol,ik)=0.d0
c        ckcsnew(maxmls+1,ik)=0.d0
c        ckssnew(maxmls+1,ik)=0.d0 
      enddo
      swp(iguest) = 1
      swp(jguest) = 1     
c     delete first guest
      estepi=0.d0
      call deletion 
     &(imcon,keyfce,iguest,ichoice,alpha,rcut,delr,drewd,maxmls,
     &totatm,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,delrc,estepi,linitsurf,surftol,sumchg,engsictmp,chgtmp,
     &overlap,newld)
      call accept_move
     &(iguest,.false.,.true.,.false.,
     &linitsurf,delrc,totatm,ichoice,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
c     insert second guest
      jmol=locguest(jguest)
      jnatms=numatoms(jmol)
      jnmols=nummols(jmol)

      do iatm=1,jnatms
        newx(iatm) = guestx(jguest,iatm) + comx
        newy(iatm) = guesty(jguest,iatm) + comy
        newz(iatm) = guestz(jguest,iatm) + comz
      enddo
      call random_rot(idnode,rotangle,q1,q2,q3,q4)
      call rotation(newx,newy,newz,comx,comy,comz,
     & jnatms,q1,q2,q3,q4)

      estepj=0.d0
      call insertion
     & (imcon,jguest,keyfce,alpha,rcut,delr,drewd,totatm,
     & volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     & engunit,delrc,estepj,sumchg,chgtmp,engsictmp,maxmls,
     & loverlap,lnewsurf,surftol,overlap,newld)
      jchoice=0
      call accept_move
     &(jguest,.true.,.false.,.false.,
     &lnewsurf,delrc,totatm,jchoice,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)

c     acceptance criteria
      accepted=.false.
      if(.not.loverlap)then
        rande=duni(idnode)
        estep = estepj - estepi
        call energy_eval
     &(estep,rande,statvolm,iguest,jguest,temp,beta,
     &.false.,.false.,.false.,.true.,accepted)
      endif
c     DEBUG
c      accepted=.false.
c     END DEBUG
      if(accepted)then
        accept_swap=accept_swap+1
        call condense(totatm,ntpfram,ntpguest)
      else
c       restore original framework if move is rejected
        call reject_move
     &(iguest,jguest,.false.,.false.,.false.,.true.)
        do ik=1,totatm
          molxxx(imol,ik)=origmolxxx(imol,ik)
          molyyy(imol,ik)=origmolyyy(imol,ik)
          molzzz(imol,ik)=origmolzzz(imol,ik)
          molxxx(jmol,ik)=origmolxxx(jmol,ik)
          molyyy(jmol,ik)=origmolyyy(jmol,ik)
          molzzz(jmol,ik)=origmolzzz(jmol,ik)
        enddo
c       restore original ewald1 sums if step is rejected
        do ik=1,newld
          ckcsum(imol,ik)=ckcsorig(imol,ik) 
          ckssum(imol,ik)=ckssorig(imol,ik) 
          ckcsum(jmol,ik)=ckcsorig(jmol,ik) 
          ckssum(jmol,ik)=ckssorig(jmol,ik) 
          ckcsum(maxmls+1,ik)=ckcsorig(maxmls+1,ik) 
          ckssum(maxmls+1,ik)=ckssorig(maxmls+1,ik)
c          ckcsnew(imol,ik)=0.d0
c          ckssnew(imol,ik)=0.d0
c          ckcsnew(jmol,ik)=0.d0
c          ckssnew(jmol,ik)=0.d0
c          ckcsnew(maxmls+1,ik)=0.d0
c          ckssnew(maxmls+1,ik)=0.d0 
        enddo
        chgsum_mol(imol)=chgsum_molorig(imol)
        chgsum_mol(jmol)=chgsum_molorig(jmol)
c       restore original surfacemols if step is rejected
        surfacemols(imol)=origsurfmols(imol)
        surfacemols(jmol)=origsurfmols(jmol)
        elrc = origelrc
        totatm = origtotatm
        nummols(locguest(iguest)) = nummols(locguest(iguest)) + 1
        nummols(locguest(jguest)) = nummols(locguest(jguest)) - 1
        call condense(totatm,ntpfram,ntpguest)
        energy = origenergy
      endif
      return
      end subroutine mc_swap

      subroutine widom_grid
     &(idnode,iguest,nstep,wngrida,wngridb,wngridc,rotangle,imcon,
     &keyfce,alpha,rcut,delr,drewd,totatm,volm,kmax1,kmax2,kmax3,epsq,
     &dlrpot,ntpatm,maxvdw,engunit,delrc,maxmls,overlap,newld,surftol,
     &sumchg,rucell,prodcount,widomcount,temp,ngrida,ngridb,ngridc)
c*****************************************************************************
c
c     perform widom insertions on a grid of pre-defined resolution 
c
c*****************************************************************************
      implicit none
      logical loverlap,lnewsurf
      integer i,j,k,istep,nstep,ngrida,ngridb,ngridc,widomcount,idnode
      integer natms,iguest,mol,nmols,keyfce,imcon,totatm,ntpatm,maxvdw
      integer kmax1,kmax2,kmax3,maxmls,newld,prodcount,istat,aa
      integer wngrida,wngridb,wngridc,loopa,loopb,loopc
      real(8) afr,bfr,cfr,cartx,carty,cartz,q1,q2,q3,q4,rotangle,surftol
      real(8) alpha,rcut,delr,delrc,drewd,volm,epsq,dlrpot,engunit
      real(8) overlap,chgtmp,engsictmp,sumchg,H_const,temp,estep
      real(8) comx,comy,comz
      real(8), dimension(9) :: rucell
      widomcount=0
      mol=locguest(iguest)
      natms=numatoms(mol)
      istat=1+16*(iguest-1)

      rotangle=2.d0*pi
      comx=0.d0
      comy=0.d0
      comz=0.d0
      !TODO(pboyd): add option to iterate over super/or unit cell
      ! may improve sampling if supercell, but otherwise pointless
      loopa = ngrida
      loopb = ngridb
      loopc = ngridc
      do i=1,loopa
        do j=1,loopb
          do k=1,loopc
            ! I think this works!
            afr=dble(i)/dble(wngrida)
            bfr=dble(j)/dble(wngridb)
            cfr=dble(k)/dble(wngridc)
            call cartesian(afr,bfr,cfr,cartx,carty,cartz)
            do istep=1,nstep
c             populate newx,newy,newz with guestx,guesty,guestz
              do aa=1,natms
                newx(aa)=guestx(iguest,aa)
                newy(aa)=guesty(iguest,aa)
                newz(aa)=guestz(iguest,aa)
              enddo
c             random rotation
              call random_rot(idnode,rotangle,q1,q2,q3,q4)
              call rotation(newx,newy,newz,comx,comy,comz,natms,
     &q1,q2,q3,q4)
c             insert COM at gridpoint
              do aa=1,natms
                newx(aa)=newx(aa)+cartx
                newy(aa)=newy(aa)+carty
                newz(aa)=newz(aa)+cartz
              enddo
              estep = 0.d0
              loverlap=.false.
              lnewsurf=.false.
              widomcount = widomcount + 1
              call insertion
     & (imcon,iguest,keyfce,alpha,rcut,delr,drewd,totatm,
     & volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     & engunit,delrc,estep,sumchg,chgtmp,engsictmp,maxmls,
     & loverlap,lnewsurf,surftol,overlap,newld)
              call reject_move
     &  (iguest,0,.true.,.false.,.false.,.false.)

              if((.not.loverlap).and.(lprobeng(iguest)))then
                call storeenprob(iguest,0,rucell,ntpguest,ngrida,
     &            ngridb,ngridc,estep)
              endif
              if(.not.loverlap)then
                prodcount = prodcount + 1
                H_const=dexp(-1.d0*estep/kboltz/temp)
c               no rolling average here, just div by widcount at the end
                chainstats(istat+7) = chainstats(istat+7)+H_const
              endif
            enddo
          enddo
        enddo
      enddo
      chainstats(istat+7) = chainstats(istat+7)/dble(prodcount)
      chainstats(1) = dble(prodcount)

      end subroutine widom_grid
      subroutine single_point
     &(imcon,keyfce,alpha,drewd,rcut,delr,totatm,ntpfram,
     &ntpguest,ntpmls,volm,newld,kmax1,kmax2,kmax3,epsq,dlrpot,
     &ntpatm,maxvdw,spenergy,vdwsum,ecoul,dlpeng,maxmls,surftol)
c*****************************************************************************
c
c     does a single point energy calculation over all atoms in the
c      system
c
c*****************************************************************************
      implicit none
      logical latmsurf,lmolsurf
      integer i,ik,j,kmax1,kmax2,kmax3,imcon,keyfce
      integer totatm,newld,sittyp,ntpatm,maxvdw
      integer mol,ntpguest,ntpfram,maxmls
      integer jatm,iatm,ntpmls,step,iguest
      real(8) drewd,dlrpot,volm,epsq,alpha,rcut,delr,ecoul
      real(8) engsic,chg
      real(8) ewald1sum,ewald2sum,vdwsum
      real(8) ewald1eng,ewald2eng,ewald3eng,vdweng
      real(8) spenergy,dlpeng,surftol
      ewald2sum = 0.d0
      ewald3eng = 0.d0
      vdwsum = 0.d0
      ewald1eng = 0.d0
      ewald2eng = 0.d0
      vdweng = 0.d0 
      do ik=1,maxmls
        vdwen(ik)=0.d0
        ewald2en(ik)=0.d0
        ewald3en(ik)=0.d0
        ewald1en(ik)=0.d0
        delrc_mol(ik)=0.d0
        elrc_mol(ik)=0.d0
      enddo
      spenergy=0.d0
c     long range correction to short range forces.
c     this initializes arrays, a different subroutine
c     is called during the gcmc simulation
c     "gstlrcorrect"
      step=0
      iguest=1
      call lrcorrect(imcon,keyfce,totatm,ntpatm,maxvdw,
     &rcut,volm,maxmls,ntpguest)

c     generate neighbour list
      call condense(totatm,ntpfram,ntpguest)
      call parlst(imcon,totatm,rcut,delr)
c     reciprocal space ewald calculation
      call ewald1(imcon,ewald1sum,engsic,totatm,volm,alpha,sumchg,
     &kmax1,kmax2,kmax3,epsq,newld,maxmls)
c     compute for all molecules, the excluding contributions
c     to the reciprocal space sum.
      call excluding_ewald_charges
     &(imcon,ntpmls,ntpguest,ewald3eng,alpha,epsq)

c     compute short range interactions, molecule-by-molecule
c     integer 'i' counts from the atom beginning to the end.
      i=0

      do mol=1,ntpmls
        if(nummols(mol).gt.0)then
          do imol=1,nummols(mol)
            lmolsurf=.false.
            do iatm=1,numatoms(mol)
c             increment counter over all atoms
              i=i+1
c              aa = ltpsit(mol,i)
              ik=0
              if(lentry(i).gt.0)then
                do j=1,lentry(i)
                  ik=ik+1

                  jatm=list(i,j)
                  ilist(j)=jatm
                  moldf(j)=moltype(jatm)
                  xdf(ik)=xxx(i)-xxx(jatm)
                  ydf(ik)=yyy(i)-yyy(jatm)
                  zdf(ik)=zzz(i)-zzz(jatm)
                enddo
                call images(imcon,ik,cell,xdf,ydf,zdf)
                do j=1,lentry(i)
                  rsqdf(j)=xdf(j)**2+ydf(j)**2+zdf(j)**2
c                 check if a surface atom
                  call surface_check
     &(i,j,mol,surftol,overlap,loverlap,latmsurf)
                  if(latmsurf)lmolsurf=.true.
                enddo
                chg=atmcharge(i)
                sittyp=ltype(i)
                
                call ewald2
     &          (chg,lentry(i),ewald2eng,mol,
     &          drewd,rcut,epsq)
                ewald2sum=ewald2sum+ewald2eng
c               calc vdw interactions
                call srfrce
     &          (sittyp,lentry(i),mol,vdweng,rcut,dlrpot)
                vdwsum=vdwsum+vdweng
              endif
            enddo
            if(lmolsurf)surfacemols(mol)=surfacemols(mol)+1
          enddo
        endif
      enddo
c      print *,"EWALD1 ",ewald1sum/engunit,-ewald3eng/engunit
c      print *,"EWALD2 ",ewald2sum/engunit
c      print *,"VDW    ",vdwsum/engunit
c      print *,"DELRC  ",elrc/engunit
      dlpeng=(ewald1sum+ewald2sum+vdwsum-ewald3eng+elrc)
     
      ecoul = ewald2sum + ewald1sum - ewald3eng
      return
      end subroutine single_point
      subroutine excluding_ewald_charges
     &(imcon,ntpmls,ntpguest,ewald3sum,alpha,epsq)
c***********************************************************************
c                                                                      * 
c     Compute for each molecule type specified in the FIELD file       *
c     the intramolecular charge interactions that will be cancelled    *
c     from the long range ewald contribution.                          *
c     This includes 'frozen' atoms that will not interact with each    *
c     other in a framework, as well as rigid guest molecules which     *
c     move, but do not interact with each other.                       *
c                                                                      *
c***********************************************************************
      implicit none
      logical lguest,lfrzi,lfrzj
      integer nmols,mol,ntpmls,imol,iatm,itgst,gstmol
      integer ntpguest,imcon,katm
      real(8) ewald3mol,ewald3sum,engcpe,alpha,epsq
      real(8) chg
      ewald3sum=0.d0
      ewald3en(:)=0.d0
c     count total number of atoms
      do mol=1,ntpmls
        lguest=.false.
        ewald3mol=0.d0
        nmols=nummols(mol)
        natms=numatoms(mol)
c       check if a guest molecule, here we assume that
c       non-bonded interactions between atoms in a guest are 
c       excluded
        do itgst=1,ntpguest
          gstmol=locguest(itgst)
          if (gstmol.eq.mol)then
c           add the total number of atoms in the molecule
c           to the atom count. This is in case there are
c           guests already in the framework at the beginning
c           of the run.
            lguest=.true.
            do iatm=1,natms-1
              ik=0
              do jatm=iatm+1,natms
                ik=ik+1
                jlist(ik)=jatm
                xdf(ik)=guestx(itgst,iatm)-guestx(itgst,jatm)
                ydf(ik)=guesty(itgst,iatm)-guesty(itgst,jatm)
                zdf(ik)=guestz(itgst,iatm)-guestz(itgst,jatm)
              enddo
              call images(imcon,ik,cell,xdf,ydf,zdf)
              chg=atmchg(mol,iatm)
              call ewald3(chg,mol,ik,alpha,engcpe,epsq)
              ewald3mol=ewald3mol+engcpe
              ewald3en(mol)=ewald3en(mol)+engcpe
            enddo
            ewald3sum=ewald3sum+nmols*ewald3mol
          endif
        enddo
c       after this it will be framework atoms populating
c       the list.
c       NB This part doesn't work currently. It doesn't matter
c       since the framework atoms are fixed and the ewald1 energy
c       will not change. the SIC for this is constant. But this 
c       is driving me NUTS.
        iatm=0
        if(.not.lguest)then
          do imol=1,nmols
            do i=1,natms
              iatm=iatm+1
              jatm=iatm
              ik=0
              lfrzi=(lfzsite(mol,i).ne.0)
              do jmol=imol,nmols
c               do not count the last interaction 
c               it is iatm=jatm=natms
                katm=1
                if (imol.eq.jmol)katm=iatm+1
                do j=katm,natms
                  jatm=jatm+1 
                  lfrzj=(lfzsite(mol,j).ne.0)
c                 this is the same atom in the same image
                  if((lfrzi).and.(lfrzj))then
                    ik=ik+1
                    jlist(ik)=jatm
                    xdf(ik)=molxxx(mol,iatm)-molxxx(mol,jatm)
                    ydf(ik)=molyyy(mol,iatm)-molyyy(mol,jatm)
                    zdf(ik)=molzzz(mol,iatm)-molzzz(mol,jatm)
                  endif
                enddo
              enddo
              call images(imcon,ik,cell,xdf,ydf,zdf)
              chg=atmchg(mol,iatm)
              call ewald3(chg,mol,ik,alpha,engcpe,epsq)
              ewald3en(mol)=ewald3en(mol)+engcpe
c             don't add this contribution to the total.
c             It is a constant factor and will be cancelled
c             out when computing most interactions.
            enddo
          enddo
        endif
      enddo 
      return
      end subroutine excluding_ewald_charges
      subroutine insert_guest
     &(imcon,iguest,apos,bpos,cpos,angx,angy,angz,
     &keyfce,alpha,rcut,delr,drewd,totatm,volm,
     &kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,engunit,maxmls,
     &delrc,estep,sumchg,surftol,overlap,newld)
c***********************************************************************
c                                                                      *
c     insert a guest molecule at a given lattice coordinate (COM)      *
c     and orientation about the Eueler axes.                           *
c                                                                      *
c***********************************************************************
      implicit none
      logical loverlap,lnewsurf
      integer imcon,iguest,natms,mol,keyfce
      integer totatm,kmax1,kmax2,kmax3,idum
      integer ntpatm,maxvdw,i,maxmls,newld
      real(8) apos,bpos,cpos,angx,angy,angz,xc,yc,zc
      real(8) alpha,rcut,delr,drewd,volm,epsq,dlrpot
      real(8) engunit,delrc,estep,sumchg,surftol
      real(8) chgtmp,engsictmp,overlap
c     convert rand fractions to cartesian coordinates in the cell 
      call cartesian(apos,bpos,cpos,xc,yc,zc)
c     xc,yc,zc are the coordinates of the com guestx,guesty,guestz 
c     are the positions of atom i relative to the com
      mol=locguest(iguest)
      natms=numatoms(mol)
      estep=0.d0
      do i=1, natms
        newx(i)=guestx(iguest,i)
        newy(i)=guesty(iguest,i)
        newz(i)=guestz(iguest,i)
      enddo
c     rotate
      call rotationeuler(newx,newy,newz,natms,angx,angy,angz) 
c     add the com to each atom in the guest
      
      do i=1, natms
        newx(i)=newx(i)+xc
        newy(i)=newy(i)+yc
        newz(i)=newz(i)+zc
      enddo

      call insertion
     & (imcon,iguest,keyfce,alpha,rcut,delr,drewd,totatm,
     & volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     & engunit,delrc,estep,sumchg,chgtmp,engsictmp,maxmls,
     & loverlap,lnewsurf,surftol,overlap,newld)

      call accept_move
     &(iguest,.true.,.false.,.false.,
     &lnewsurf,delrc,totatm,idum,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)

      end subroutine insert_guest
      subroutine displace_guest
     &(imcon,alpha,rcut,delr,drewd,totatm,newld,
     &maxmls,volm,kmax1,kmax2,kmax3,
     &epsq,engunit,overlap,surftol,linitsurf,lnewsurf,loverlap,
     &iguest,imol,dlrpot,sumchg,a,b,c,q1,q2,q3,q4,estep)
c***********************************************************************
c                                                                      *
c     Displace a guest from its initial point to a pre-defined         *
c     coordinate (therefore useable for both random jumps, and         *
c     normal displacements)                                            *
c                                                                      *
c***********************************************************************
      implicit none
      logical loverlap,linitsurf,lnewsurf
      integer imcon,totatm,maxmls
      integer kmax1,kmax2,kmax3,iguest,imol,i,mol,newld
      integer nmols,natms,ik,mxcmls
      real(8) alpha,rcut,delr,drewd,volm,epsq,overlap,surftol
      real(8) dlrpot,sumchg,comx,comy,comz,q1,q2,q3,q4,a,b,c
      real(8) estep,enginit,engunit,ewld3sum,chgtmp,engsictmp
      real(8) ewld1eng,ewld2sum,vdwsum
      call get_guest(iguest,imol,mol,natms,nmols)
c     Calculate the energy for the chosen guest in
c     its orginal place
      mxcmls=maxmls*(maxmls-1)/2 + maxmls
      mol=locguest(iguest) 
      call guest_energy
     &(imcon,iguest,alpha,rcut,delr,drewd,
     &totatm,maxmls,volm,kmax1,kmax2,kmax3,epsq,dlrpot,
     &engunit,vdwsum,ewld2sum,ewld1eng,linitsurf,newld,
     &surftol,sumchg,chgtmp,engsictmp,loverlap,overlap,estep,.true.)
      ewld3sum=ewald3en(mol)
c      print *, "EWALD1:", (-ewld1eng-ewld3sum)/engunit
c      print *, "EWALD2:", ewld2sum/engunit
c      print *, "VDW   :", vdwsum/engunit
c      print *, "DELRC :", delrc/engunit
c      print *, "ESTEP :", estep - delrc/engunit

      enginit=estep
      do ik=1,mxcmls
        ewald1entmp(ik) = ewald1en(ik)
        ewald2entmp(ik) = ewald2en(ik)
        vdwentmp(ik) = vdwen(ik)
      enddo
      do ik=1,newld
        ckcsum(mol,ik)=ckcsum(mol,ik)-ckcsnew(mol,ik)
        ckssum(mol,ik)=ckssum(mol,ik)-ckssnew(mol,ik)
        ckcsum(maxmls+1,ik)=ckcsum(maxmls+1,ik)-ckcsnew(maxmls+1,ik)
        ckssum(maxmls+1,ik)=ckssum(maxmls+1,ik)-ckssnew(maxmls+1,ik)
c        ckcsnew(mol,ik)=0.d0
c        ckssnew(mol,ik)=0.d0
c        ckcsnew(maxmls+1,ik)=0.d0 
c        ckssnew(maxmls+1,ik)=0.d0
      enddo
      call translate
     &(natms,mol,newx,newy,newz,cell,rcell,comx,comy,comz,
     &a,b,c)
      call rotation(newx,newy,newz,comx,comy,comz,natms,q1,q2,q3,q4)
      call guest_energy
     &(imcon,iguest,alpha,rcut,delr,drewd,
     &totatm,maxmls,volm,kmax1,kmax2,kmax3,epsq,dlrpot,
     &engunit,vdwsum,ewld2sum,ewld1eng,lnewsurf,newld,
     &surftol,sumchg,chgtmp,engsictmp,loverlap,overlap,estep,.false.)
c      print *, "EWALD1:", (-ewld1eng-ewld3sum)/engunit
c      print *, "EWALD2:", ewld2sum/engunit
c      print *, "VDW   :", vdwsum/engunit
c      print *, "DELRC :", delrc/engunit
c      print *, "ESTEP :", estep + delrc/engunit 
c     total up the energy contributions.
      estep=estep-enginit
      do i=1, maxmls
        ik=loc2(mol,i)
        if(i.ne.mol)delE(i)=delE(i)+
     &(ewald1en(ik)-ewald1entmp(ik)
     &+ewald2en(ik)-ewald2entmp(ik)
     &+vdwen(ik)-vdwentmp(ik))
     &/engunit
      enddo
      delE(mol)=delE(mol)+
     &  estep

      return
      end subroutine displace_guest


      subroutine test
     &(imcon,idnode,keyfce,alpha,rcut,delr,drewd,totatm,dlrpot,newld,
     &ntpguest,volm,kmax1,kmax2,kmax3,epsq,ntpatm,maxvdw,surftol,
     &engunit,ntpfram,maxmls,outdir,cfgname,levcfg,overlap)
c***********************************************************************
c
c     testing for various bugs etc in fastmc. 
c
c***********************************************************************
      implicit none
      logical loverlap,lnewsurf,linitsurf
      integer imcon,idnode,iguest,keyfce,totatm,ntpguest
      integer kmax1,kmax2,kmax3,ntpatm,maxvdw,newld
      real(8) alpha,rcut,delr,drewd,volm,epsq,engunit
      integer i,ntpfram,randchoice,maxmls
      integer levcfg,imol
      real(8) apos,bpos,cpos
      real(8) angx,angy,angz,chgtmp,engsictmp
      real(8) delrc,estep,surftol
      real(8) sumchg,eng,overlap
      real(8) a,b,c,q1,q2,q3,q4
      real(8) rotangle,dlrpot
      character*8 outdir
      character*1 cfgname(80)      

      randchoice=0
      sumchg=0.d0
      write(*,*)"TESTING"
      iguest=1
      apos=0.5d0; bpos=0.5d0; cpos=0.5d0
      angx=0.d0; angy=90.d0; angz=0.d0
      call insert_guest
     &(imcon,iguest,apos,bpos,cpos,angx,angy,angz,
     &keyfce,alpha,rcut,delr,drewd,totatm,volm,
     &kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,engunit,maxmls,
     &delrc,estep,sumchg,surftol,overlap,newld)
      print *,"Insertion: ",iguest,estep
      do i=1,maxmls
        print *, "ENERGY ",i,energy(i)
      enddo
      iguest=1
      apos=0.5d0; bpos=0.5d0; cpos=0.7d0
      angx=0.d0; angy=90.d0; angz=0.d0
      call insert_guest
     &(imcon,iguest,apos,bpos,cpos,angx,angy,angz,
     &keyfce,alpha,rcut,delr,drewd,totatm,volm,
     &kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,engunit,maxmls,
     &delrc,estep,sumchg,surftol,overlap,newld)
      print *,"Insertion: ",iguest,estep
      do i=1,maxmls
        print *, "ENERGY ",i,energy(i)
      enddo
      iguest=2
      apos=0.5d0; bpos=0.5d0; cpos=0.3d0
      angx=0.d0; angy=90.d0; angz=0.d0
      call insert_guest
     &(imcon,iguest,apos,bpos,cpos,angx,angy,angz,
     &keyfce,alpha,rcut,delr,drewd,totatm,volm,
     &kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,engunit,maxmls,
     &delrc,estep,sumchg,surftol,overlap,newld)
      print *,"Insertion: ",iguest,estep
      do i=1,maxmls
        print *, "ENERGY ",i,energy(i)
      enddo

      iguest=1; imol=1
      call random_disp
     &(idnode,delr,a,b,c)
      rotangle=1.d0
      call random_rot
     &(idnode,rotangle,q1,q2,q3,q4)

c      print *, "a=",a,";b=",b,";c=",c
c      print *,"q1=",q1,";q2=",q2,";q3=",q3,";q4=",q4
      a=0.33263498544692993;b=0.40358364582061768
      c=-8.1818878650665283E-002
      q1=0.88821917901419711;q2=-0.22256139119037360
      q3=-0.36817090934697771;q4=0.16119335809322474    
      call displace_guest
     &(imcon,alpha,rcut,delr,drewd,totatm,newld,
     &maxmls,volm,kmax1,kmax2,kmax3,
     &epsq,engunit,overlap,surftol,linitsurf,lnewsurf,loverlap,
     &iguest,imol,dlrpot,sumchg,a,b,c,q1,q2,q3,q4,estep)
      print *, loverlap
      call accept_move
     &(iguest,.false.,.false.,.true.,
     &lnewsurf,delrc,totatm,imol,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
      print *, "Displacement: ",iguest,imol,estep
      do i=1,maxmls
        print *, "ENERGY ",i,energy(i)
      enddo

      iguest=1; imol=1
      call deletion 
     &(imcon,keyfce,iguest,imol,alpha,rcut,delr,drewd,maxmls,
     &totatm,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,delrc,estep,linitsurf,surftol,sumchg,engsictmp,chgtmp,
     &overlap,newld)
      call accept_move
     &(iguest,.false.,.true.,.false.,
     &lnewsurf,delrc,totatm,imol,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
      print *, "Deletion: ",iguest,imol,estep
      do i=1,maxmls
        print *, "ENERGY ",i,energy(i)
      enddo

      iguest=2; imol=1
      call deletion 
     &(imcon,keyfce,iguest,imol,alpha,rcut,delr,drewd,maxmls,
     &totatm,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,delrc,estep,linitsurf,surftol,sumchg,engsictmp,chgtmp,
     &overlap,newld)
      call accept_move
     &(iguest,.false.,.true.,.false.,
     &lnewsurf,delrc,totatm,imol,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
      print *, "Deletion: ",iguest,imol,estep
      do i=1,maxmls
        print *, "ENERGY ",i,energy(i)
      enddo
      iguest=1; imol=1
      call deletion 
     &(imcon,keyfce,iguest,imol,alpha,rcut,delr,drewd,maxmls,
     &totatm,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,delrc,estep,linitsurf,surftol,sumchg,engsictmp,chgtmp,
     &overlap,newld)
      call accept_move
     &(iguest,.false.,.true.,.false.,
     &lnewsurf,delrc,totatm,imol,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
      print *, "Deletion: ",iguest,imol,estep
      do i=1,maxmls
        print *, "ENERGY ",i,energy(i)
      enddo

      eng = 0.d0 
      call revive
     &(totatm,levcfg,production,ntpguest,ntpmls,
     &imcon,cfgname,eng,outdir)

      end subroutine test
      end
