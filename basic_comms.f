      subroutine initcomms()
      
c*********************************************************************
c     
c     communication harness initialisation
c     
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c     wl
c     2005/09/28 12:04:35
c     1.1
c     Exp
c     
c*********************************************************************
      
      implicit none
      
      include "comms.inc"
      
      integer ierr

CMPIU      define MPI_init MPI_init_

      call MPI_init(ierr)

      return
      end

      subroutine machine(idnode,mxnode)

c*********************************************************************
c     
c     subroutine for obtaining charcteristics of
c     the computer on which the program is being run
c     
c*********************************************************************

      implicit none

      integer idnode,mxnode,mynode,numnodes

      idnode=mynode()
      mxnode=numnodes()

      return
      end

      integer function mynode()

c*********************************************************************
c
c     routine to determine identity of processing node 
c
c     MPI version - t.forester may 1995
c
c     wl
c     2005/09/28 12:04:35
c     1.1
c     Exp
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer ierr

CMPIU define MPI_COMM_RANK MPI_COMM_RANK_

      call MPI_COMM_RANK(MPI_COMM_WORLD, mynode ,ierr)

      return
      end

      integer function nodedim()

c*********************************************************************
c
c     calculate dimension of hypercube
c
c     MPI version - t.forester may 1995
c
c     wl
c     2005/09/28 12:04:35
c     1.1
c     Exp
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer i,n,ierr,mxnode

CMPIU      define MPI_COMM_SIZE MPI_COMM_SIZE_

      call MPI_COMM_SIZE(MPI_COMM_WORLD, mxnode ,ierr)
      n=1
      nodedim = -1
      do i=0,16

         if(n.eq.mxnode)nodedim=i
         n=2*n

      enddo

      return
      end

      integer function numnodes()

c*********************************************************************
c
c     calculate number of nodes
c
c     MPI version - t.forester may 1995
c
c     wl
c     2005/09/28 12:04:35
c     1.1
c     Exp
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer ierr

CMPIU      define MPI_COMM_SIZE MPI_COMM_SIZE_

      call MPI_COMM_SIZE(MPI_COMM_WORLD, numnodes, ierr)

      return
      end

      subroutine csend(msgtag,buf,length,pe,idum)

c*********************************************************************
c
c     Intel-like  csend (double precision)
c
c     idum=1 send real(kind=16)
c     idum=2 send real(kind=8)
c     idum=3 send integer
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer msgtag,length,pe,idum

      integer len1,ierr
      real*8 buf(*)

CMPIU      define MPI_send MPI_send_

      len1 = length/Dlen
c      call MPI_send(buf,len1,MPI_DOUBLE_PRECISION,pe,msgtag,
c     x     MPI_COMM_WORLD,ierr)

c     set idum value to mean quadruple precision... dont think it 
c     was being used for anything
      if(idum.eq.1)then
        call MPI_send(buf,length,MPI_REAL16,pe,msgtag,
     x       MPI_COMM_WORLD,ierr)
      elseif(idum.eq.2)then
        call MPI_send(buf,length,MPI_DOUBLE_PRECISION,pe,msgtag,
     x       MPI_COMM_WORLD,ierr)
      else if(idum.eq.3)then
        call MPI_send(buf,length,MPI_INTEGER,pe,msgtag,
     x       MPI_COMM_WORLD,ierr)

      endif
      return
      end

      subroutine crecv(msgtag,buf,length,idum)

c*********************************************************************
c
c     Intel-like  crecv (double precision)
c
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c
c     wl
c     2005/09/28 12:04:35
c     1.1
c     Exp
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer msgtag,length

      integer len1,ierr,idum
      integer status(MPI_STATUS_SIZE)
      real*8 buf(*)

      len1 = length/Dlen

CMPIU      define MPI_RECV MPI_RECV_

c      call MPI_RECV(buf,len1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,
c     x     msgtag,MPI_COMM_WORLD,status,ierr)

      if(idum.eq.1)then
        call MPI_RECV(buf,length,MPI_REAL16,MPI_ANY_SOURCE,
     x       msgtag,MPI_COMM_WORLD,status,ierr)
      elseif(idum.eq.2)then
        call MPI_RECV(buf,length,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,
     x       msgtag,MPI_COMM_WORLD,status,ierr)
      elseif(idum.eq.3)then
        call MPI_RECV(buf,length,MPI_INTEGER,MPI_ANY_SOURCE,
     x       msgtag,MPI_COMM_WORLD,status,ierr)
      endif

      return 
      end

      subroutine gisum2(aaa,nnn,bbb)

c***********************************************************************
c     
c     Like gisum, but the local data isn't overwritten 
c
c***********************************************************************
      
      

      implicit none

      integer nnn,i,ierror,iii,kk,k,k0,k1,k2,msg1,msg2
      integer aaa(nnn),bbb(nnn)

      include "comms.inc"

      integer status(MPI_STATUS_SIZE)

CMPIU      define MPI_allreduce MPI_allreduce_

      call MPI_allreduce(aaa,bbb,nnn,MPI_INTEGER,
     x  MPI_SUM,MPI_COMM_WORLD,ierror)

      return
      end

      subroutine gisum(aaa,nnn,bbb)

c***********************************************************************
c     
c     global summation subroutine for hypercube - MPI version
c     integer version
c     
c***********************************************************************
      
      

      implicit none

      integer nnn,i,ierror,iii,kk,k,k0,k1,k2,msg1,msg2
      integer aaa(nnn),bbb(nnn)

      include "comms.inc"

      integer status(MPI_STATUS_SIZE)

CMPIU      define MPI_allreduce MPI_allreduce_

      call MPI_allreduce(aaa,bbb,nnn,MPI_INTEGER,
     x  MPI_SUM,MPI_COMM_WORLD,ierror)

      do i = 1,nnn
        aaa(i) = bbb(i)
      enddo

      return
      end

      subroutine gisum3(aaa,iii,mmm,nnn,ccc)

c***********************************************************************
c   
c     this subroutine sums an array with two dimensions   
c     to arrays with one dimension
c     - to be honest I don't know much about MPI, sums can probably
c       be done accross multiple dimensional arrays.  I'm just
c       too lazy to look this up, so 1d arrays it is.
c
c***********************************************************************
      
      

      implicit none

      integer nnn,i,mmm,ierror,iii,kk,k,k0,k1,k2,msg1,msg2
      integer aaa(mmm,nnn),ccc(nnn)
      integer, dimension(nnn) :: bbb

      include "comms.inc"

      integer status(MPI_STATUS_SIZE)

      do i=1,nnn
        bbb(i)=aaa(iii,i)
      enddo

CMPIU      define MPI_allreduce MPI_allreduce_

      call MPI_allreduce(bbb,ccc,nnn,MPI_INTEGER,
     x  MPI_SUM,MPI_COMM_WORLD,ierror)

      do i = 1,nnn
        aaa(iii,i) = ccc(i)
      enddo

      return
      end

      subroutine gdsum(aaa,nnn,bbb)

c***********************************************************************
c     
c     global summation subroutine for MPI - hypercube assumed
c     double precision version
c     
c***********************************************************************

      implicit none

      integer nnn,i,iii,kk,k1,k2,ierror
      real(8) aaa(nnn),bbb(nnn)

      include "comms.inc"

      integer status(MPI_STATUS_SIZE)

CMPIU      define MPI_allreduce MPI_allreduce_

      call MPI_allreduce(aaa,bbb,nnn,MPI_DOUBLE_PRECISION,
     x  MPI_SUM,MPI_COMM_WORLD,ierror)

        do i = 1,nnn
          aaa(i) = bbb(i)
        enddo

      return
      end

      subroutine gd2dsum(nnn,first,snd1,snd2,totatm)

c***********************************************************************
c     
c     global summation for a 2d array.. kind of messy
c
c***********************************************************************

      implicit none

      integer first,i,j,k,iatm,jatm,totatm
      integer kk,k1,k2,k0msg1,msg2,ierror
      integer snd1(first), snd2(first)
      real(8)  nnn(first,10000)
      real(8)  aaa(totatm),bbb(totatm)

     
      include "comms.inc"

      integer status(MPI_STATUS_SIZE)

CMPIU      define MPI_allreduce MPI_allreduce_

      iatm=0
      do i=1,first
        jatm=0
        do j=1,snd1(i)
          do k=1,snd2(i)
            iatm=iatm+1
            jatm=jatm+1
            aaa(iatm)=nnn(i,jatm)
          enddo
        enddo
      enddo
            
      call MPI_allreduce(aaa,bbb,totatm,MPI_DOUBLE_PRECISION,
     x  MPI_SUM,MPI_COMM_WORLD,ierror)

      iatm=0
      do i=1,first
        jatm=0
        do j=1,snd1(i)
          do k=1,snd2(i)
            iatm=iatm+1
            jatm=jatm+1
            nnn(i,jatm)=bbb(iatm)
          enddo
        enddo
      enddo 
    
      return
      end

      subroutine gdsum3(aaa,iii,mmm,nnn,ccc)

c***********************************************************************
c
c     this subroutine sums an array with two dimensions
c     to arrays with one dimension
c     double equivalent of gisum3
c
c***********************************************************************

      implicit none

      integer nnn,i,mmm,ierror,iii,kk,k,k0,k1,k2,msg1,msg2
      real(8) aaa(mmm,nnn),ccc(nnn)
      real(8), allocatable ::  bbb(:)
      integer fail
      include "comms.inc"
      integer status(MPI_STATUS_SIZE)
      allocate(bbb(nnn),stat=fail)
      do i=1,nnn
        bbb(i)=aaa(iii,i)
      enddo
CMPIU      define MPI_allreduce MPI_allreduce_
      call MPI_allreduce(bbb,ccc,nnn,MPI_DOUBLE_PRECISION,
     x  MPI_SUM,MPI_COMM_WORLD,ierror)
      do i = 1,nnn
        aaa(iii,i) = ccc(i)
      enddo
      return
      end

      subroutine gimax(aaa,nnn,bbb)

c***********************************************************************
c     
c     global maximum subroutine for hypercube - MPI version
c     integer version
c     
c***********************************************************************
      

      implicit none

      integer nnn,i,iii,kk,k1,k2,k,k0msg1,msg2,ierror
      integer aaa(nnn),bbb(nnn)

      include "comms.inc"

      integer status(MPI_STATUS_SIZE)
CMPIU      define MPI_allreduce MPI_allreduce_
      
      call MPI_allreduce(aaa,bbb,nnn,MPI_INTEGER,
     x   MPI_MAX,MPI_COMM_WORLD,ierror)
      
      do i = 1,nnn
        aaa(i) = bbb(i)
      enddo

      return
      end

      subroutine gstate(check)

c***********************************************************************
c     
c     global status subroutine : gisum version
c     
c***********************************************************************


      implicit none

      logical check
      integer i,idum

      i = 0
      if(.not.check) i = 1

      call gisum(i,1,idum)
      
      check = (i.eq.0)

      return
      end

      subroutine gsync()

c*********************************************************************
c     
c     barrier / synchronization routine
c
c     MPI version - t.forester may 1995
c     CPP version - w.smith
c
c     wl
c     2005/09/28 12:04:35
c     1.1
c     Exp
c
c*********************************************************************

      implicit none

      integer ierr

      include "comms.inc"

CMPIU      define MPI_BARRIER MPI_BARRIER_

      call  MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

      subroutine exitcomms()

c*********************************************************************
c
c     exitcomms: exit from communication harness
c
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c
c     wl
c     2005/09/28 12:04:35
c     1.1
c     Exp
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer ierr
CMPIU      define MPI_FINALIZE MPI_FINALIZE_

      call MPI_FINALIZE(ierr)
      call exit(0)

      return
      end
