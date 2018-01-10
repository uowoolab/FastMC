      module wang_landau
      use mpmodule
      use utility_pack
      use wang_landau_module
      use mc_moves
      use ewald_module

      contains
      subroutine wang_landau_sim
     &(idnode,mxnode,imcon,keyfce,alpha,rcut,delr,drewd,totatm,ntpguest,
     &ntpfram,volm,statvolm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,sumchg,maxmls,surftol,overlap,newld,outdir,levcfg,cfgname,
     &wlprec,mcinsf,mcdelf,mcdisf,mcjmpf,mcflxf,mcswpf,mctraf,prectol,
     &mcrotf,mcswif,nnumg,temp,beta,mcsteps,eqsteps,flatcoeff,visittol,
     &maxn,minn)
c*****************************************************************************
c     
c     Main routine for computing the weighted histogram of Wang and
c     Landau
c     PB - 21/08/2017
c
c*****************************************************************************
      implicit none
      logical wlchk,loverlap,lnewsurf,lprod,accepted
      logical insert,delete,displace,jump,swap,tran,rota,switch
      logical production,converge
      character*25 outfile,visfile
      character*31 dosfile
      character*8 outdir
      character*1 cfgname(80)      
      integer idnode,mxnode,imcon,keyfce,ntpguest,kmax1,kmax2,kmax3
      integer ntpatm,maxvdw,newld,maxmls,totatm,levcfg,nhist,ihist
      integer natms,iguest,jguest,mol,imol,idum,ntpfram,nnumg,minchk
      integer accept_ins,ins_count,accept_del,del_count,accept_tran
      integer tran_count,accept_disp,disp_count,accept_rota,rota_count
      integer totaccept,wlcount,accept_jump,jump_count,accept_switch
      integer switch_count,accept_swap,swap_count,mcsteps,eqsteps,tail
      integer wlstepcount,prod_count,minmol,maxmol,insmol,molidx,nmol
      integer varchunk,visittol,sweepcount,sweepsteps,i,j,k,n,maxn,minn
      integer ivarchunk,imaxmol,iminmol,iter
      real(8) alpha,rcut,delr,drewd,volm,epsq,dlrpot,engunit
      real(8) sumchg,surftol,overlap,estep,chgtmp,randmov
      real(8) engsictmp,delrc,logwlprec,wlprec,timelp,beta,statvolm
      real(8) mcinsf,mcdelf,mcdisf,mcjmpf,mcflxf,mcswpf,mctraf
      real(8) mcrotf,mcswif,lambda,temp,tran_delr,rota_rotangle
      real(8) delrdisp,rotangle,jumpangle,flatcoeff,prectol,mu
      if(idnode.eq.0)then
        write(nrite, 
     &"(/'Entering main routine for Wang - Landau calculation',/)")
        write(nrite,
     &"('Initial coefficient set to ',f9.5,/)")wlprec
        write(nrite,"('Convergence tolerance for the coefficient ',
     &E9.1,/)")prectol
      endif
      wlchk=.true.
      logwlprec=log(wlprec)
      data accept_ins,accept_del,accept_disp,totaccept/0,0,0,0/
      data accept_jump,accept_swap,accept_switch/0,0,0/
      data accept_tran,accept_rota/0,0/
      data ins_count,del_count,disp_count,wlcount/0,0,0,0/
      data jump_count,swap_count,switch_count/0,0,0/
      data tran_count,rota_count/0,0/
      data wlstepcount,sweepcount,sweepsteps,prod_count/0,0,0,0/ 

c     make a new file to write some Wang-Landau data to
      if(ntpguest.gt.1)then
        do i=1,ntpguest
          write(outfile,"(a8,'/wang_landau',i2.2,'.out')")outdir,i
          open(800+i,file=outfile)
          write(visfile,"(a8,'/visited_hist',i2.2,'.csv')")outdir,i
          open(700+i,file=visfile)
          write(dosfile,"(a8,'/denisty_of_states',i2.2,'.csv')")
     &outdir,i
          open(600+i,file=dosfile)
        enddo
      else
        outfile=outdir // '/wang_landau.out'
        open(801,file=outfile,status="replace")
        visfile=outdir // '/visited_hist.csv'
        open(701,file=visfile,status="replace")
        dosfile=outdir // '/density_of_states.csv'
        open(601,file=dosfile,status="replace")
      endif
      production=.false.
      rotangle=pi/3.d0
      delrdisp=delr
      tran_delr=delr
      rota_rotangle=rotangle
      jumpangle=rotangle
      ! set the number of concurrent histograms sampling.
      ! We will hard-code this to one and sample only the number of
      ! guests for now... 
      nhist=1
      ihist=nhist
      call alloc_wl_arrays(idnode,nhist,ntpguest,maxn)
c     Temporary: exit if more than one guest included in the
c     FIELD/CONTROL files. Make clear that this currently works
c     for estimating the partition function for a single guest
c     at a single temperature.
      if(ntpguest.gt.1)call error(idnode,2318)
      call timchk(0,timelp)
c     obtain the min/max number of molecules based on this node's identity
c     and maxn
      tail=0
      call divide_jobs
     &(idnode,mxnode,maxn,minn,maxmol,minmol,varchunk,tail)
      if(mxnode.eq.1)then
        write(nrite,
     &"('Sampling ',i5,' guest loadings on one node.',/)")
     &(varchunk-1)
      else
        if(idnode.eq.0)then
          write(nrite,"('Sampling split into ',i2,
     &' different jobs.')")mxnode
          do i=1,mxnode
            call divide_jobs
     &(i-1,mxnode,maxn,minn,imaxmol,iminmol,ivarchunk,tail)
            write(nrite, "(4x,'Rank ',i2,' running ',i4,' to ',i4)")
     &i-1,iminmol,imaxmol

          enddo
          call flush(nrite)
        endif
      endif
c     loop is kind of pointless for now but keep it for future
c     use.
      do i=1,ntpguest
        ! compute thermal debroglie wavelength for each guest
        ! assuming, 1) equilibrium ideal gas, 2) no internal degrees of freedom
        lambda = hplanck / sqrt(2.d0*pi*molmass(i)/1000.d0*boltz*temp/
     &avo)/1.d-10
        dlambda(i) = lambda**3
        if((guest_insert(i).gt.0).or.(minmol.gt.0))then
          imol=locguest(i)
          insmol=max(guest_insert(i), minmol)
          write(800+i,
     &"(1x,'Inserting ',i6,' guests of type ',i3)")
     &insmol,i
          call insert_guests
     &(idnode,imcon,totatm,ntpguest,ntpfram,i,insmol,
     &rcut,delr,sumchg,surftol,overlap,keyfce,alpha,drewd,volm,newld,
     &kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,engunit,delrc,
     &maxmls,iter)
          write(800+i,"(1x,'Successful Insertion of ',i6,
     &' guests of type ',i3)")
     &nummols(imol),i
          write(800+i,"(1x,i9,' trials. Success rate: ',f6.2,' %'/)")
     &iter,dble(nummols(mol))/dble(iter) * 100.d0
        endif
      enddo
      minchk=min(400,nnumg)
      do while(wlchk)
c       randomly select guest, should only be one for now
        iguest=floor(duni(idnode)*ntpguest)+1
        imol=locguest(iguest)
c       set all moves to false.
        insert=.false.
        delete=.false.
        displace=.false.
        jump=.false.
        swap=.false.
        tran=.false.
        rota=.false.
        switch=.false.
c       chose a move to perform
c       accept/reject based on modified acceptance criteria
        randmov=duni(idnode)
        if(randmov.lt.mcinsf)then
          insert = .true.
        elseif(randmov.lt.mcdelf)then
          delete = .true.
        elseif(randmov.lt.mcdisf)then
          displace = .true.
        elseif(randmov.lt.mcjmpf)then
          jump = .true.
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
c       Failover displace -- shouldn't reach here
          displace=.true.
        endif
c       Border conditions. Insertion when nmol is greater or equal to
c       maxmol
        nmol=nummols(imol)
        if(nmol.lt.minmol)then
          insert=.true.
          delete=.false.
          displace=.false.
          jump=.false.
          swap=.false.
          tran=.false.
          rota=.false.
          switch=.false.
c       Border condition. Deletion when nmol is less than or equal to
c       minmol
c       NB: by not changing the mcxyzf values, these border conditions
c       are favouring another move, which may not be intended by the
c       user. e.g. mcdisf has a large advantage now that mcdelf is
c       removed.
        elseif((nmol.le.minmol).and.(delete))then
c         update edge case
          visit_hist(nmol+1,ihist)=visit_hist(nmol+1,ihist)+1
          dos_hist(nmol+1,ihist)=dos_hist(nmol+1,ihist)+logwlprec
          tmat_c(nmol+1,nmol+1)=tmat_c(nmol+1,nmol+1)+1.d0
          delete=.false.
          if(randmov.lt.mcinsf)then
            insert = .true.
          elseif(randmov.lt.mcdisf)then
            displace = .true.
          elseif(randmov.lt.mcjmpf)then
            jump = .true.
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
c         Failover displace -- shouldn't reach here
            displace=.true.
          endif
        endif
        if(nmol.gt.maxmol)then
          insert=.false.
          delete=.true.
          displace=.false.
          jump=.false.
          swap=.false.
          tran=.false.
          rota=.false.
          switch=.false.
        elseif((nmol.ge.maxmol).and.(insert))then
c         update edge case
          visit_hist(nmol+1,ihist)=visit_hist(nmol+1,ihist)+1
          dos_hist(nmol+1,ihist)=dos_hist(nmol+1,ihist)+logwlprec
          tmat_c(nmol+1,nmol+1)=tmat_c(nmol+1,nmol+1)+1.d0
          insert=.false.
          if(randmov.lt.mcdelf)then
            delete = .true.
          elseif(randmov.lt.mcdisf)then
            displace = .true.
          elseif(randmov.lt.mcjmpf)then
            jump = .true.
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
c         Failover displace -- shouldn't reach here
            displace=.true.
          endif
        endif
c       fist two are controlled by new WL acceptance criteria
c       the rest are usual NVT.
c       NB: Swap and Switch not used yet in these calcs.
        if(insert)then
          call wl_insert
     &(idnode,imcon,keyfce,iguest,totatm,rcut,delr,ins_count,alpha,
     &drewd,ntpguest,ntpfram,ntpatm,volm,statvolm,kmax1,kmax2,kmax3,
     &epsq,dlrpot,maxvdw,newld,engunit,delrc,sumchg,maxmls,surftol,
     &overlap,accepted,temp,beta,accept_ins,minmol,ihist,logwlprec)
          insert=.false.
        elseif(delete)then
          call wl_delete
     &(idnode,imcon,keyfce,iguest,totatm,rcut,delr,del_count,alpha,
     &drewd,ntpguest,ntpfram,ntpatm,volm,statvolm,kmax1,kmax2,kmax3,
     &epsq,dlrpot,maxvdw,newld,engunit,delrc,sumchg,maxmls,surftol,
     &overlap,accepted,temp,beta,accept_del,minmol,ihist,logwlprec)
          delete=.false.
        elseif(displace)then
          call wl_displace
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,tran,
     &tran_count,tran_delr,rota,rota_count,rota_rotangle,disp_count,
     &delrdisp,rotangle,maxmls,kmax1,kmax2,kmax3,newld,alpha,rcut,delr,
     &drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,accepted,temp,
     &beta,delrc,ntpatm,maxvdw,accept_tran,accept_rota,accept_disp,
     &ntpguest,ntpfram)
          displace=.false.
          tran=.false.
          rota=.false.
        elseif(jump)then
          call wl_jump
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,
     &jump_count,jumpangle,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_jump,ntpfram,
     &ntpguest)
          jump=.false.
        elseif(switch)then
          call wl_switch
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,
     &switch_count,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_switch,ntpfram,
     &ntpguest)
          switch=.false.
        elseif(swap)then
          call wl_swap
     &(idnode,imcon,keyfce,iguest,jguest,totatm,volm,statvolm,
     &swap_count,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_swap,ntpfram,
     &ntpguest,rotangle)
          swap=.false.
        endif
c       offset molidx by one to count for the '0' bin.
        !molidx=(nummols(imol)-minmol)+1
        molidx=nummols(imol)+1
c       update the histograms (regardless of accept/reject)
c       DOS Histogram is updated by an increment of log(wlprec)
c       Visited state histogram incremented by 1.
        !write(*,*)visit_hist(:, ihist)
c        this creates a monotonically increasing DOS. Maybe just for
c        insert/deletes?
c        visit_hist(molidx,ihist)=visit_hist(molidx,ihist)+1
c        dos_hist(molidx,ihist)=dos_hist(molidx,ihist)+logwlprec
c       Convergence of a particular 'sweep' is determined when
c       the visited state histogram is 'flat'
        converge=convergence_check
     &(idnode,ihist,minmol,varchunk,flatcoeff,visittol)
        if(converge)then
          sweepcount=sweepcount+1
          ! update the precision factor
          wlprec=adjust_factor(wlprec)
          logwlprec=log(wlprec)
          write(800+iguest,"('Sweep ',i7)")sweepcount
          write(800+iguest,"('Histogram is flat, rescaling precision ',
     &'factor to ',f15.12)")wlprec
          totaccept=
     &accept_ins+accept_del+accept_disp+accept_jump+
     &accept_swap+accept_switch+accept_tran+accept_rota
          write(800+iguest,"('Number of accepted steps ',i9)")totaccept
          write(800+iguest,"('Total number of steps this sweep ',i9)")
     &sweepsteps
          sweepsteps=0
          call dump_wl_files
     &(idnode,ntpguest,nhist,minmol,maxmol,varchunk)
          call flush(800+iguest)
        endif
        if(accepted)then
          accepted=.false.
        endif
        wlstepcount=wlstepcount+1
        sweepsteps=sweepsteps+1
        if(production)then
c          if(prod_count.ge.mcsteps)wlchk=.false.
          prod_count=prod_count+1
        endif
c       done if the precision factor is less than the tolerance
        if((converge).and.((logwlprec).lt.prectol))then
          wlchk=.false.
        elseif(converge)then
          ! only reset histogram if we are not at the end
          ! of the run. This way the 'dump_wl_files' will
          ! not produce a '0' histogram.
          call reset_visit_hist(ihist,minmol,varchunk)
          
        endif
        if(wlstepcount.ge.eqsteps)production=.true.
      enddo
      call timchk(0,timelp)
      if(mxnode.gt.1)call report_timing
     &(idnode,mxnode,timelp,varchunk)
c     DUMP here, before everything gets merged to the 
c     mother node.
      call dump_wl_files
     &(idnode,ntpguest,nhist,minmol,maxmol,varchunk)

c     useful routine for summing DOS's from each node.
c     HOWEVER - the normalization of the DOS is in this routine
c     as well, so must keep and make sure it works for serial
c     executions as well.
      call stitch_branches
     &(idnode,mxnode,maxn,minn,minmol,maxmol,tail,ihist,varchunk)
      ! units = kcal/kg
      lprod=.true.
      !DOSFILE is 601, VISITED file is 701
      call revive
     &(totatm,levcfg,lprod,ntpguest,maxmls,
     &imcon,cfgname,0.d0,outdir)
      if(idnode.eq.0)then

        call compute_isotherm
     &(idnode,iguest,beta,ihist,maxn,minn,varchunk)
        write(nrite,"(/,a100,/,33x,'Wang-Landau calculation complete!',
     &/,a100,/)")repeat('*',100),repeat('*',100)
        write(nrite, "(a35,f15.12)")'tolerance met: ',(wlprec-1.d0)
        write(nrite, "(a35,i15)")
     &'total number of Monte Carlo steps: ',wlstepcount
        write(nrite,"(a35,i15)")
     &'number of W-L sweeps: ',sweepcount
        

        write(nrite, "(a35,f15.3,' seconds')")
     &'Time elapsed for W-L calculation: ',timelp
        ! write DOS histogram for guest i
        ! first two loops are currently pointless.
      endif 
      call terminate_wl
     &(idnode,ntpguest)
      end subroutine wang_landau_sim

      logical function convergence_check
     &(idnode,ihist,minmol,varchunk,flatcoeff,visittol)
c**********************************************************************
c                                                                     *
c     check to see if the histogram is converged.                     *
c     PB - 15/11/2017                                                 *
c                                                                     *
c**********************************************************************
      implicit none                                                   
      integer ihist,ik,ncount,varchunk,visittol,minmol,idnode 
      real(8) ave,var,std,vhist,delta,delta2,flatcoeff,mxvar,minvar 
      convergence_check=.false.                                    
      ave=0.d0                                                    
      var=0.d0
      ncount=0
      mxvar=0.d0
      minvar=1.d99 
      ! mean and std in one loop
      do ik=1,varchunk+1
        ncount=ncount+1
        
        vhist=dble(visit_hist(minmol+ik,ihist))
        if(vhist.lt.minvar)minvar=vhist
        if(vhist.gt.mxvar)mxvar=vhist
        ave=ave+vhist
        !delta=vhist-ave
        !ave=ave+delta/dble(ncount)
        !delta2=vhist-ave
        !var=var+delta*delta2
      enddo
      ave=ave/dble(ncount)
c     std could be another metric for convergence
c      std=sqrt(var)
c     return if there's at least one bin that hasn't
c     been visited yet.
      if(minvar.le.0.5d0)return
      convergence_check=
     &((minvar.ge.flatcoeff*ave).or.(int(minvar).ge.visittol))
      !write(*,*)"*****************************************"
      !write(*,*)convergence_check
      !write(*,*)"FLAT TOL      :",flatcoeff
      !write(*,*)"AVERAGE       :",ave
      !write(*,*)"MIN VAL       :",minvar
      !write(*,*)"AVE*TOL       :",flatcoeff*ave
      !write(*,*)"*****************************************"
      return
      end function convergence_check

      subroutine report_timing(idnode,mxnode,timelp,varchunk)
c**********************************************************************
c                                                                     * 
c     Report the parallel timings of the code. Useful data            *
c     for deciding on how to partition the N guest loadings, and      *
c     other convergence criteria options.                             *
c     PB - 14/12/2017                                                 *
c                                                                     *
c**********************************************************************
      implicit none
      integer idnode,mxnode,varchunk,varbuff,inode
      real(8) timelp,timebuff,maxtime,mintime
      real(8) :: timedata(mxnode)
      integer :: chunkdata(mxnode)

      maxtime=timelp
      mintime=timelp
      if(idnode.eq.0)then
        timedata(1)=timelp
        chunkdata(1)=varchunk
        do inode=1,mxnode-1
          call crecv(inode*mxnode*2+1,timebuff,1,2)  
          call crecv(inode*mxnode*3+1,varbuff,1,3)
          if(timebuff.gt.maxtime)maxtime=timebuff
          if(timebuff.lt.mintime)mintime=timebuff
          timedata(inode+1)=timebuff
          chunkdata(inode+1)=varbuff 
        enddo
        write(nrite,"(a100,/,40x,'Timing Data',/,a100)")
     &repeat('*',100),repeat('*',100)
        write(nrite,"('Load imbalance: ',f7.2,' %')")
     &(1.d0-mintime/maxtime)*100.d0
        write(nrite,"(/,5x,a5,10x,a10,10x,a8,10x,a17,/,a79)")
     &"node","data size","time(s)","load imbalance(%)",
     &repeat('-',74)
        do inode=1,mxnode
          write(nrite,"(5x,i3,10x,i7,10x,f14.5,12x,f7.3)")inode-1,
     &chunkdata(inode),timedata(inode),
     &(1.d0-timedata(inode)/maxtime)*100.d0
        enddo
        timelp=maxtime
      else
        call csend(idnode*mxnode*2+1,timelp,1,0,2)
        call csend(idnode*mxnode*3+1,varchunk,1,0,3)
      endif
      end subroutine report_timing

      subroutine reset_visit_hist(ihist,minmol,varchunk)
c**********************************************************************
c                                                                     *
c     Reset the visit histogram to '0'                                *
c     PB - 30/11/2017                                                 *
c                                                                     *
c**********************************************************************
      implicit none
      integer ik,ihist,varchunk,minmol
      do ik=1,varchunk+1
        visit_hist(minmol+ik,ihist)=0
      enddo
      end subroutine reset_visit_hist

      real(8) function adjust_factor(f)
c**********************************************************************
c                                                                     *
c     Function to decrease the density of states scaling factor       *
c     PB - 15/11/2017                                                 *
c                                                                     *
c**********************************************************************
      implicit none
      real(8) f
      adjust_factor = sqrt(f)
      return
      end function adjust_factor

      logical function eval_wl_move
     &(idnode,iguest,ihist,insert,delete,swap,estep,minmol,volm,
     &beta)
c**********************************************************************
c                                                                     *
c     Acceptance criteria for the Wang-Landau algorithm. Currently    *
c     hard-coded for ajusting the number of molecules as a macrostate *
c     variable.                                                       *
c     PB - 15/11/2017                                                 *
c                                                                     *
c**********************************************************************
      implicit none
      logical insert,delete,swap
      integer idnode,iguest,imol,nmol,molidx,minmol,ihist
      real(8) randn,contrib,estep,beta,volm,unbiased
      eval_wl_move=.false.
      imol=locguest(iguest)
      nmol=nummols(imol)
      !molidx=(nmol-minmol)+1
      molidx=nmol+1
      if(insert)then
        unbiased=exp(-1.d0*beta*estep)*volm/(1+nmol)
        contrib=unbiased/dlambda(iguest)*
     &exp(dos_hist(molidx,ihist)-dos_hist(molidx+1,ihist))
c       TODO(pboyd): Determine if biased contribution comes
c       from tmat_c or dos_hist.
        !write(*,*)"INSERT unbiased: ",unbiased
      elseif(delete)then
        unbiased=exp(-1.d0*beta*estep)*(nmol)/volm
        contrib=unbiased*dlambda(iguest)*
     &exp(dos_hist(molidx,ihist)-dos_hist(molidx-1,ihist))
        !write(*,*)"DELETE unbiased: ",unbiased
      else
        ! Nothing yet..
      endif
      randn=duni(idnode)
      !if(randn.lt.contrib)then
      !  write(*,*)"************************************"
      !  if(delete)write(*,*)"DELETION:  "
      !  if(insert)write(*,*)"INSERTION: "
      !  write(*,*)"ENERGY  : ",estep
      !  write(*,*)"UNBIASED: ",unbiased
      !  write(*,*)"BIASED  : ",contrib
      !  write(*,*)"LAMBDA^3: ",dlambda(iguest)
      !  write(*,*)"************************************"
      !endif
      eval_wl_move = (randn.lt.contrib)
      return
      end function eval_wl_move

      subroutine wl_insert
     &(idnode,imcon,keyfce,iguest,totatm,rcut,delr,ins_count,alpha,
     &drewd,ntpguest,ntpfram,ntpatm,volm,statvolm,kmax1,kmax2,kmax3,
     &epsq,dlrpot,maxvdw,newld,engunit,delrc,sumchg,maxmls,surftol,
     &overlap,accepted,temp,beta,accept_ins,minmol,ihist,logwlprec)
c*******************************************************************************
c                                                                              *
c     keeps track of all the associated arrays and calls the 'insertion'       *
c     subroutine. Acceptance/rejection criteria based on Metropolis            *
c     importance sampling coupled to a ?? ensemble.                            *
c                                                                              *
c*******************************************************************************
      implicit none                                                    
      logical lnewsurf,loverlap,accepted
      integer iguest,idnode,imcon,natms,totatm,mol,nmol
      integer ins_count,keyfce,ntpguest,kmax1,kmax2,kmax3
      integer ntpatm,maxvdw,maxmls,newld,accept_ins
      integer ntpfram,randchoice,minmol,ihist,molidx,old_molidx
      real(8) rcut,delr,estep,alpha,drewd,statvolm,logwlprec
      real(8) epsq,dlrpot,engunit,delrc,sumchg,chgtmp,engsictmp
      real(8) surftol,overlap,rande,volm,temp,beta,unbiased
      mol=locguest(iguest)
c     nmol is the original number of molecules. i.e. not N+1
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

      if (.not.loverlap)accepted=eval_wl_move
     &(idnode,iguest,ihist,.true.,.false.,.false.,estep,minmol,volm,
     &beta)
c     DEBUG
c      accepted=.true.
c     END DEBUG

c     update TM arrays regardless of acceptance/rejection
      !molidx=(nmol-minmol)+2
      !old_molidx=(nmol-minmol)+1
      molidx=nmol+2
      old_molidx=nmol+1
      !TODO(pboyd) look up proper unbiased probability..
      unbiased=exp(-1.d0*beta*estep)*volm/(1+nmol)
      if(unbiased.gt.1.d0)then
        tmat_c(old_molidx,molidx)=tmat_c(old_molidx,molidx)+1.d0
c        tmat_c(molidx,old_molidx)=tmat_c(molidx,old_molidx)+
c     &(1.d0-1.d0)
        tmat_vis(old_molidx,molidx)=tmat_vis(old_molidx,molidx)+1
        tmat_vis(old_molidx,old_molidx)=
     &tmat_vis(old_molidx,old_molidx)+1
c        tmat_vis(molidx,old_molidx)=
c     &tmat_vis(molidx,old_molidx)+1

      else
        tmat_c(old_molidx,molidx)=tmat_c(old_molidx,molidx)+unbiased
c        tmat_c(molidx,old_molidx)=tmat_c(molidx,old_molidx)+
c     &(1.d0-unbiased)
        tmat_c(old_molidx,old_molidx)=tmat_c(old_molidx,old_molidx)+
     &(1.d0-unbiased)
        tmat_vis(old_molidx,molidx)=tmat_vis(old_molidx,molidx)+1
        tmat_vis(old_molidx,old_molidx)=
     &tmat_vis(old_molidx,old_molidx)+1
c        tmat_vis(molidx,old_molidx)=
c     &tmat_vis(molidx,old_molidx)+1
      endif

      if(accepted)then
        accept_ins=accept_ins+1
        randchoice=0
        visit_hist(molidx,ihist)=visit_hist(molidx,ihist)+1
        dos_hist(molidx,ihist)=dos_hist(molidx,ihist)+logwlprec
        call accept_move
     &(iguest,.true.,.false.,.false.,
     &lnewsurf,delrc,totatm,randchoice,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
      else
        visit_hist(old_molidx,ihist)=visit_hist(old_molidx,ihist)+1
        dos_hist(old_molidx,ihist)=dos_hist(old_molidx,ihist)+logwlprec
        call reject_move
     &(iguest,0,.true.,.false.,.false.,.false.)
      endif
      return
      end subroutine wl_insert
      subroutine wl_delete
     &(idnode,imcon,keyfce,iguest,totatm,rcut,delr,del_count,alpha,
     &drewd,ntpguest,ntpfram,ntpatm,volm,statvolm,kmax1,kmax2,kmax3,
     &epsq,dlrpot,maxvdw,newld,engunit,delrc,sumchg,maxmls,surftol,
     &overlap,accepted,temp,beta,accept_del,minmol,ihist,logwlprec)
c*******************************************************************************
c                                                                              *
c     keeps track of all the associated arrays and calls the 'deletion'        *
c     subroutine. Acceptance/rejection criteria based on Metropolis            *
c     importance sampling coupled to                                           *
c                                                                              *
c*******************************************************************************
      implicit none
      logical linitsurf,accepted
      integer iguest,idnode,imcon,totatm,mol,nmol
      integer del_count,keyfce,ntpguest,kmax1,kmax2,kmax3
      integer ntpatm,maxvdw,maxmls,newld,accept_del
      integer ntpfram,randchoice,minmol,ihist,molidx,old_molidx
      real(8) rcut,delr,estep,alpha,drewd,statvolm,logwlprec
      real(8) epsq,dlrpot,engunit,delrc,sumchg,chgtmp,engsictmp
      real(8) surftol,overlap,rande,volm,temp,beta,unbiased
      mol=locguest(iguest)
c     nmol is the original number of molecules. i.e. not N-1
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

      accepted=eval_wl_move
     &(idnode,iguest,ihist,.false.,.true.,.false.,-estep,minmol,
     &volm,beta)

c     update TM arrays regardless of acceptance/rejection
      !molidx=(nmol-minmol)
      !old_molidx=(nmol-minmol)+1
      molidx=nmol
      old_molidx=nmol+1
      !TODO(pboyd) look up proper unbiased probability..
      ! NB: Witman has this as 1/exp(1.d0*beta*estep)
      unbiased=exp(-1.d0*beta*estep)*(nmol)/volm
      if(unbiased.gt.1.d0)then
        tmat_c(old_molidx,molidx)=tmat_c(old_molidx,molidx)+1.d0
c        tmat_c(molidx,old_molidx)=tmat_c(molidx,old_molidx)+
c     &(1.d0-1.d0)
        tmat_vis(old_molidx,molidx)=tmat_vis(old_molidx,molidx)+1
        tmat_vis(old_molidx,old_molidx)=
     &tmat_vis(old_molidx,old_molidx)+1
c        tmat_vis(molidx,old_molidx)=
c     &tmat_vis(molidx,old_molidx)+1

      else
        tmat_c(old_molidx,molidx)=tmat_c(old_molidx,molidx)+unbiased
c        tmat_c(molidx,old_molidx)=tmat_c(molidx,old_molidx)+
c     &(1.d0-unbiased)
        tmat_c(old_molidx,old_molidx)=tmat_c(old_molidx,old_molidx)+
     &(1.d0-unbiased)
        tmat_vis(old_molidx,molidx)=tmat_vis(old_molidx,molidx)+1
        tmat_vis(old_molidx,old_molidx)=
     &tmat_vis(old_molidx,old_molidx)+1
c        tmat_vis(molidx,old_molidx)=
c     &tmat_vis(molidx,old_molidx)+1
      endif
         
c     the following occurs if the move is accepted.
      if(accepted)then
        visit_hist(molidx,ihist)=visit_hist(molidx,ihist)+1
        dos_hist(molidx,ihist)=dos_hist(molidx,ihist)+logwlprec
        accept_del=accept_del+1
        call accept_move
     &(iguest,.false.,.true.,.false.,
     &linitsurf,delrc,totatm,randchoice,ntpfram,ntpguest,maxmls,
     &sumchg,engsictmp,chgtmp,newld)
      else
        visit_hist(old_molidx,ihist)=visit_hist(old_molidx,ihist)+1
        dos_hist(old_molidx,ihist)=dos_hist(old_molidx,ihist)+logwlprec
        call reject_move
     &(iguest,0,.false.,.true.,.false.,.false.)
      endif
      return
      end subroutine wl_delete

      subroutine wl_displace
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,tran,
     &tran_count,tran_delr,rota,rota_count,rota_rotangle,disp_count,
     &delrdisp,rotangle,maxmls,kmax1,kmax2,kmax3,newld,alpha,rcut,delr,
     &drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,accepted,temp,
     &beta,delrc,ntpatm,maxvdw,accept_tran,accept_rota,accept_disp,
     &ntpguest,ntpfram)
c*******************************************************************************
c                                                                              *
c     Keeps track of all the associated arrays and calls the                   *
c     'wl_displace_guest' subroutine.                                          *
c     The underlying mechanics of this code is a 'deletion/insertion'          *
c     move, where the molecule is only perturbed as far as delrdisp and        *
c     rotangle dictates.                                                       *
c     Acceptance/rejection criteria based on Metropolis                        *
c     importance sampling coupled to                                           *
c                                                                              *
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

      call wl_displace_guest
     &(imcon,alpha,rcut,delr,drewd,totatm,newld,
     &maxmls,volm,kmax1,kmax2,kmax3,
     &epsq,engunit,overlap,surftol,linitsurf,lnewsurf,loverlap,
     &iguest,randchoice,dlrpot,sumchg,a,b,c,q1,q2,q3,q4,estep)
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
      end subroutine wl_displace

      subroutine wl_jump
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,
     &jump_count,jumpangle,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_jump,ntpfram,
     &ntpguest)
c*******************************************************************************
c                                                                              *
c     Keeps track of all the associated arrays and calls the                   *
c     'wl_displace_guest' subroutine.                                          *
c     The underlying mechanics of this code is a 'deletion/insertion'          *
c     move, where the molecule is randomly perturbed to another point          *
c     in the simulation cell.                                                  *
c     Acceptance/rejection criteria based on Metropolis                        *
c     importance sampling coupled to                                           *
c                                                                              *
c*******************************************************************************
      implicit none
      logical accepted,linitsurf,lnewsurf,loverlap
      integer idnode,imcon,keyfce,iguest,totatm,jump_count,maxmls,kmax1
      integer kmax2,kmax3,newld,ntpatm,maxvdw,accept_jump,randchoice
      integer nmol,natms,mol,ik,ntpfram,ntpguest
      real(8) volm,statvolm,jumpangle,alpha,rcut,delr,drewd,epsq,engunit
      real(8) overlap,surftol,dlrpot,sumchg,temp,beta,delrc
      real(8) a,b,c,q1,q2,q3,q4,estep,gpress,rande,engsictmp,chgtmp

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

      call random_jump
     &(idnode,a,b,c,jumpangle)
      call random_rot(idnode,jumpangle,q1,q2,q3,q4)

      call wl_displace_guest
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
      end subroutine wl_jump

      subroutine wl_switch
     &(idnode,imcon,keyfce,iguest,totatm,volm,statvolm,
     &switch_count,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_switch,ntpfram,
     &ntpguest)
c********************************************************************************
c                                                                               *
c     Keeps track of all the associated arrays and calls the                    *
c     The underlying mechanics of this code is a                                *
c     'deletion / displacement / insertion'                                     *
c     move, where two molecules of different types in the simulation cell       *
c     switch locations.                                                         *
c     Acceptance/rejection criteria based on Metropolis                         *
c     importance sampling coupled to                                            *
c                                                                               *
c     *** Currently only one switch is performed, but one could possibly        *
c         include multi switches as long as (detailed) balance isn't            *
c         violated. ***                                                         *
c********************************************************************************
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
      end subroutine wl_switch

      subroutine wl_swap
     &(idnode,imcon,keyfce,iguest,jguest,totatm,volm,statvolm,
     &swap_count,maxmls,kmax1,kmax2,kmax3,newld,alpha,
     &rcut,delr,drewd,epsq,engunit,overlap,surftol,dlrpot,sumchg,
     &accepted,temp,beta,delrc,ntpatm,maxvdw,accept_swap,ntpfram,
     &ntpguest,rotangle)
c*******************************************************************************
c                                                                              *
c     Keeps track of all the associated arrays and calls the                   *
c     'deletion' and 'insertion' subroutines.                                  *
c     The underlying mechanics of this code is a 'deletion/insertion'          *
c     move, where a molecule in the simulation cell is swapped with a          *
c     molecule with a different identity.                                      *
c     Acceptance/rejection criteria based on Metropolis                        *
c     importance sampling coupled                                              *
c                                                                              *
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
      end subroutine wl_swap

      subroutine wl_displace_guest
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
      end subroutine wl_displace_guest

      subroutine dump_wl_files
     &(idnode,ntpguest,nhist,minmol,maxmol,varchunk)
c***********************************************************************
c                                                                      *
c     Write the DOS and visited histograms to their file channels      *
c     PB - 05/12/17                                                    *
c                                                                      *
c***********************************************************************
      implicit none
      integer idnode,ntpguest,i,j,k,n
      integer nhist,minmol,maxmol,varchunk
      do i=1,ntpguest
c       Overwrite old data...
        rewind(600+i)
        write(600+i,"('N,ln[Q(N)]')")
        rewind(700+i)
        write(700+i,"('N,Visit(N)')")
        do j=1,nhist
          do k=1,varchunk+1
c           subtract 1 to account for occupation of N=0
            n=(minmol+k)-1
            write(600+i,"(i6,',',f50.10)")n,dos_hist(n+1,j)
            write(700+i,"(i6,',',i20)")n,visit_hist(n+1,j)
          enddo
        enddo
        call flush(600+i)
        call flush(700+i)
      enddo
      end subroutine dump_wl_files

      real(8) function obtain_min_shift(obj,opt,nsample)
c***********************************************************************
c                                                                      *
c     Find the constant shift factor that will minimize the variance   *
c     between two sets of data, obj is the objective data to fit to,   *
c     opt is the data that will be shifted by the constant factor.     * 
c                                                                      *
c     NB: This should be a non-linear optimization..                   *
c     PB - 13/12/17                                                    *
c                                                                      *
c***********************************************************************
      implicit none
      integer nsample,i
      real(8) objav,optav
      real(8), dimension(nsample) :: obj, opt
      objav=0.d0
      optav=0.d0
      !do i=1,nsample
      !  objav=objav+obj(i)
      !  optav=optav+opt(i)
      !enddo
      !objav=objav/dble(nsample)
      !optav=optav/dble(nsample)

c      currently just a hack job.
      obtain_min_shift=obj(1)-opt(1)
c      obtain_min_shift=objav-optav
c     TODO(pboyd): report error in fit.
      end function obtain_min_shift

      subroutine stitch_branches
     &(idnode,mxnode,maxn,minn,minmol,maxmol,tail,ihist,varchunk)
c***********************************************************************
c                                                                      *
c     Routine to stitch together the density of states computed        *
c     on each cpu. This is a MPI blocking routine, as they must be     *
c     stitched incrementally, starting at node 0.                      * 
c                                                                      * 
c     PB - 13/12/17                                                    *
c                                                                      *
c***********************************************************************
      implicit none
      integer idnode,mxnode,inode,tail,ii,ik,ihist,idx,maxmol
      integer minmol,minmol2,maxmol2,ij,i,varchunk2,varchunk,minn,maxn
      real(8), allocatable :: dosbuff(:)
      ! objective and optimized
      real(8), allocatable :: obj(:), opt(:)
      real(8) const,shift
      allocate(dosbuff(maxn+1))
      allocate(obj(2*tail-2))
      allocate(opt(2*tail-2))

      if(idnode.eq.0)then
c       normalize the dos_hist on the main node
        shift=dos_hist(1,ihist)
        do i=1,varchunk
c         shift the DOS by the value recorded for N=0.
          dos_hist(i,ihist)=dos_hist(i,ihist)-shift
        enddo
        ! ranks start at 0, but FORTRAN starts at 1. Subtract
        ! mxnode by one to properly iterate over nodes.
        do inode=1,(mxnode-1)
c         crecv([messagetag],[variable to send],[length],[destination],
c     [dummy var])
c         WARNING: setting the last variable to 1 will send a
c         real(kind=16) type.
          call crecv(inode*2+1,dosbuff,maxn+1,2)
          ! ignore terminal entries
          ! merge the rest via min variance.
          ik=1
          ! obtain minmol and maxmol for the current node
          call divide_jobs
     &(inode,mxnode,maxn,minn,maxmol2,minmol2,varchunk2,tail)
          do ii = 1, 2*tail-2
            ik=ik+1
            opt(ii)=dosbuff(minmol2+ik)
            idx=minmol2+ik
            obj(ii)=dos_hist(idx,ihist)
          enddo
          const = obtain_min_shift(obj,opt,2*tail-2)
          ! report error in stitching with node inode
          ik=1
          do ij=minmol2+1,maxmol2+1
            ik=ik+1
            ! merge the tail portion by averaging.
            if((ik.ge.2).and.(ik.le.2*tail-1))then
              dos_hist(ij,ihist)=
     &(dos_hist(ij,ihist)+dosbuff(ij)+const)/2.d0
            else
              dos_hist(ij,ihist)=dosbuff(ij)+const
            endif
          enddo
        enddo
      else
c       csend([messagetag],[variable to send],[length],[destination],
c     [dummy var])
c       may have to reduce dos_hist to a 1D array..
        do ii=1,maxn+1
          dosbuff(ii)=dos_hist(ii,ihist)
        enddo
        call csend(idnode*2+1,dosbuff,maxn+1,0,2)
      endif

      end subroutine stitch_branches

      subroutine compute_isotherm
     &(idnode,iguest,beta,ihist,maxn,minn,varchunk)
c***********************************************************************
c                                                                      *
c     Determine the most probable adsorption value for the gas         *
c     for a range of pressures (currently treated as an ideal gas).    *
c     TODO(pboyd): make the pressure ranges user-definable             * 
c                                                                      * 
c     PB - 12/12/17                                                    *
c                                                                      *
c***********************************************************************
      use mpmodule
      implicit none
      character*12 outfile
      character*25 finaldos
      integer npress,idnode,iguest,ihist,ik,ip,maxn,minn,varchunk
      integer k,nid,nwds
      real(8) mu,pmin,pmax,pinterval,p,beta,niter
      !real(kind=16) z_uvt,n,dos_sum  
      type(mp_real) z_uvt,n,dos_sum,temp1,ktmp,dosmp,exp1 

      ! mpdpw is a parameter declared in the mpfuna.f90 file,
      ! but i've also declared it (commented out) in
      ! wang_landau_module.f in case the parameter in the mpfun package
      ! is not 'seen' by this program.
      nwds = int(1100/mpdpw + 2)
      ! varchunk is now the total range of N, not the amount
      ! allocated for each node
      varchunk = maxn-minn+1
      pinterval=1.d-1
      ! these are in bar.
      pmin=0.d0
      pmax=5.d0
      npress=int((pmax-pmin)/pinterval)

c     first write out the final dos in the main directory.
c     If run in parallel, this will be the 'stitched' version
c     of the DOS, performed on idnode=0
c     in subroutine 'stitch_branches'
      finaldos='final_dos.csv'
      open(16,file=finaldos,status="replace")
      write(16,"('N,ln[Q(N)]')")
      do k=1,maxn+1
c       subtract 1 to account for occupation of N=0
        nid=k-1
        write(16,"(i6,',',f50.10)")nid,dos_hist(k,ihist)
      enddo

c     compute isotherm data now.
      outfile='isotherm.csv'
      open(15,file=outfile,status="replace")

      write(15,"('p/bar,mmol/g')")
      do ip=1,npress+1
        ! convert from bar to Pa
        ! add 1 Pa to each pressure so there isn't a DIV0 error
        p=((dble(ip-1)+1d-5)*pinterval)*1.d5
        mu=log(p*dlambda(iguest)*beta)/beta 
        z_uvt = grand_canonical_partition
     &(idnode,ihist,beta,mu,minn,maxn,varchunk,nwds)
        dos_sum=mpreald(0.d0)
        do ik=1,maxn
          niter=dble(ik)-1.q0
          dosmp = mpreald(dos_hist(ik,ihist))
          ktmp = mpreald(beta*mu*niter)
          exp1 = exp(dosmp+ktmp)
          temp1 = mpreald(niter) * exp1
          dos_sum=dos_sum+temp1
        enddo
        n=dos_sum/z_uvt
        write(15,"(f15.6,',',f15.6)")p*1.d-5,n
       !write(15,*)p,n
      enddo 
      close(15)
      end subroutine compute_isotherm

      type(mp_real) function grand_canonical_partition
     &(idnode,ihist,beta,mu,minn,maxn,varchunk,nwds)
c***********************************************************************
c                                                                      *
c     Determine the grand canonical partition function from the        *
c     ln[Q(NVT)] histogram.                                            *
c                                                                      * 
c     PB - 12/12/17                                                    *
c                                                                      *
c***********************************************************************
      implicit none
      integer idnode,ihist,i,maxn,minn,varchunk,nwds
      real(8) mu,beta,niter
      type(mp_real) ktmp,dosmp

      grand_canonical_partition=mpreal(0.d0,nwds) 
      do i=1,varchunk
c       shift the DOS by the value recorded for N=0.
        niter=dble(i)-1.d0-dble(minn)
        dosmp = mpreald(dos_hist(i,ihist))
        ktmp = mpreald(beta*mu*niter)

        grand_canonical_partition = grand_canonical_partition 
     &+ exp(dosmp + ktmp)
      enddo

      end function grand_canonical_partition

      subroutine divide_jobs
     &(idnode,mxnode,maxn,minn,maxmol,minmol,varchunk,tail)
c***********************************************************************
c                                                                      *
c     Divide the jobs based on the node count. This is currently       *
c     an asymptotic function.                                          *
c     The function was designed to approach the max number of          *
c     molecules for the last compute node allocated.                   *
c     Thus the last node will sample large values for N, but very      *
c     few bins (1 or 2 + tail).                                        *
c     PB - 13/12/17                                                    *
c                                                                      *
c***********************************************************************
      implicit none
      integer idnode,mxnode,maxmol,minmol,varchunk,maxn,minn
      integer tail,nodex
c     5% addition to both sides of a N distribution so that one can 
c     stitch them together in post-processing.
      !tail=int(dble(maxn)*0.05d0)
      tail=3
      ! if run serially, just assume the user will add the tails
      ! if this calc will be stitched together in postprocessing.
      if(mxnode.eq.1)tail=0
      nodex=idnode+1
c     asymtotic function:
c     f(x) = (Nmax * x^2 + Nmax/max_cpu * x^2) / (x^2 + x + 1)
      maxmol=((maxn+1)*nodex**2+(maxn+1)/mxnode*nodex**2)/
     &(nodex**2+nodex+1)+tail
      minmol=((maxn+1)*idnode**2+(maxn+1)/mxnode*idnode**2)/
     &(idnode**2+idnode+1)-tail
c     TODO(pboyd): check if maxmol in N-1 is greater than N 
c     (due to tail)
      if(nodex.eq.1)minmol=minn
      if(nodex.eq.mxnode)maxmol=maxn

c      minmol=(maxn)/(mxnode)*(idnode)
c      maxmol=(maxn)/(mxnode)*(idnode+1)
      ! increment by 1 to account for the zeroth entry
      varchunk=maxmol-minmol
      end subroutine divide_jobs

      subroutine terminate_wl
     &(idnode,ntpguest)
c***********************************************************************
c                                                                      *
c     End the program within the Wang-Landau Module. This way all      *
c     the stuff reported for a typical GCMC sim is avoided.            *
c                                                                      * 
c     PB - 05/12/17                                                    *
c                                                                      *
c***********************************************************************
      implicit none
      integer idnode,ntpguest,i

      if(idnode.eq.0)write(nrite,"(/,/,'Successful termination of ',
     & 'Wang-Landau simulation!')")
c     close all i/o channels for each branch
      do i=1,ntpguest
        close (400+i)
        close (800+i)
        close (202)
c        close (203)
        close (500+i)
        close (205)
      enddo
      
      if(idnode.eq.0) then
        close (nrite)
        close (ncontrol)
        close (nconfig)
        close (nstats)
        close (nfield)
      endif
      call gsync()
      call exitcomms()
      return
      end subroutine terminate_wl
      end module wang_landau
