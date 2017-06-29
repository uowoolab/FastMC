      subroutine parlst(imcon,natms,fram,rcut,delr)

c***********************************************************************
c     
c     subroutine for constructing the verlet neighbour
c     list based on the brode-ahlrichs decomposition 
c     frozen atoms taken into account
c     
c     parallel replicated data version
c     
c***********************************************************************
      
      use utility_pack 
c      use exclude_module


      implicit none
      
      logical lchk,lfrzi,ldo,fram
      integer imcon,natms,ibig,last,mpm2
      integer npm2,idum,i,m,ii,j
      real(8) rclim,rsq,rcut,delr,rmin,rminsq

      
c     max size of verlet neighbour list for each atom
c      dens=dble(maxatm*ntpmls)/volm
c      cut=rcut+delr
c      ratio=1.5d0*dens*(4.d0*pi/3.d0)*cut**3
c      mxlist=min(nint(ratio),(maxatm*ntpmls+1)/2)
     
c     check size of work array

c      if(mxxdf.lt.(natms+1)/2)then
c        if(idnode.eq.0) write(nrite,*) 'mxxdf must be greater than ',
c     x    (natms+1)/2
c        call  error(idnode,474)
c      endif

c     set control variables
      last=natms
      lchk=.true.
      mpm2=natms/2
      npm2=(natms-1)/2
        
c     set cutoff radius
c     rmin to prevent framework atom interactions at too close a range
      if(fram)then
        rmin=3.d0
      else
        rmin=0.d0
      endif
      rminsq=rmin*rmin
      rclim=(rcut+delr)**2
c     construct pair force neighbour list
        
      do i=1,natms

        lentry(i)=0
        noxatm(i)=1
          
      enddo

c     outer loop over atoms
      do m=1,mpm2
          
        if(m.gt.npm2)last=mpm2
          
c     inner loop over atoms
          
        ii=0
          
        do i=1,last
            
c     calculate atom indices
            
          j=i+m
          if(j.gt.natms)j=j-natms
          ii=ii+1
          xdf(ii)=xxx(i)-xxx(j)
          ydf(ii)=yyy(i)-yyy(j)
          zdf(ii)=zzz(i)-zzz(j)
            
        enddo
c     apply minimum image convention
          
        call images(imcon,ii,cell,xdf,ydf,zdf)
          
c     allocate atoms to neighbour list
          
        ii=0
          
        do i=1,last

c******* commented out so that pairwise interactions *******
c******* will be calculated between frozen sites ********        
         lfrzi=(lfreezesite(i).ne.0)

c     calculate atom indices
            
          j=i+m
          if(j.gt.natms)j=j-natms
            
          ii=ii+1

c     reject atoms in excluded pair list
          if((nexatm(ii).gt.0).and.(lexatm(ii,noxatm(ii)).eq.j))
     x     then
            
            noxatm(ii)=min(noxatm(ii)+1,nexatm(ii))

c     reject frozen atom pairs

          else
          
            ldo=.true.

c******* commented out so that pairwise interactions *******
c******* will be calculated between frozen sites ********        
            if(lfrzi.and.(.not.fram))ldo=(lfreezesite(j).eq.0)
             
            if(ldo)then 

c     calculate interatomic distance
                  
              if(imcon.eq.6)then

                rsq=xdf(ii)*xdf(ii)+ydf(ii)*ydf(ii)

              else

                rsq=xdf(ii)*xdf(ii)+ydf(ii)*ydf(ii)+zdf(ii)*zdf(ii)

              endif

c     running check of neighbour list array capacity
                
              if((rsq.lt.rclim).and.(rsq.ge.rminsq))then
                lentry(ii)=lentry(ii)+1 
                    
                if(lentry(ii).gt.mxlist)then
                  lchk=.false.
c                  ibig=max(lentry(ii),ibig)
c
                endif
                  
c     compile neighbour list array
                  
                if(lchk)then
                  list(ii,lentry(ii))=j
                endif
                 
              endif                 
 
            endif

          endif

        enddo
          
      enddo
c     terminate job if neighbour list array exceeded
        
c      if(mxnode.gt.1) call gstate(lchk)
        
c      if(.not.lchk)then

c        call gimax(ibig,1,idum)
c        if(idnode.eq.0)then
c          write(nrite,*) ' mxlist must be at least  ',ibig
c          write(nrite,*) ' mxlist is currently ',mxlist
c        endif
c        call error(idnode,110)

c      endif

c     check all excluded atoms are accounted for
        
c      do i=1,ii
          
c        if(nexatm(i).gt.0.and.noxatm(i).ne.nexatm(i))lchk=.false.
          
c      enddo
      return
      end 
