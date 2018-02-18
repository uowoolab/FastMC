
c*************************************************************
c      DOS_TO_ISOTHERM                                       *
c*************************************************************
c      program to read in a .csv file containing the DOS     *
c      and output an isotherm.                               *
c*************************************************************
c     personal reference; channels for read/write
c     stderr 0
c     stdin  5
c     stdout 6
      use utility_pack
      use readinputs
      use wang_landau_module
      use parse_module
      use mpmodule

      implicit none
      logical loop,safe
      character*30 outval
      integer idnode,idum,ik,ihist,fail,nhist,n,i
      integer, allocatable :: n_hist(:)
      type(mp_real) z_uvt,dos_sum,temp1,ktmp,dosmp,exp1 

      idnode=0
      safe=.true.
      loop=.true.
      ! open the csv file
      open(33,file='final_dos.csv',status='old')
      ! read in the DOS
      call getrec(safe,idnode,33)
      if(.not.safe)call abort_csv_read(1,33)
      call getrec(safe,idnode,33)
      if(.not.safe)call abort_csv_read(1,33)
      nhist=0
      ihist=1
      do while(loop)
        call lowcase(record,lenrec)
        call strip(record, lenrec)
        call getrec(safe,idnode,33)
        if(record(1).eq.'#'.or.record(1).eq.' ')then
        else
            nhist=nhist+1
        endif
        if(.not.safe)loop=.false.
      enddo
      rewind(33)
      allocate(dos_hist(nhist,ihist), stat=fail)
      allocate(n_hist(nhist), stat=fail)
      call getrec(safe,idnode,33)
      if(.not.safe)call abort_csv_read(1,33)
      call getrec(safe,idnode,33)
      if(.not.safe)call abort_csv_read(1,33)
      loop=.true.
      ik=0
      do while(loop)
        call lowcase(record,lenrec)
        call strip(record, lenrec)
        if(record(1).eq.'#'.or.record(1).eq.' ')then
        else
            ik=ik+1
            n_hist(ik)=intstr(record,lenrec,idum)
            dos_hist(ik,ihist)=dblstr(record,lenrec,idum)
        endif
        call getrec(safe,idnode,33)
        if(.not.safe)loop=.false.
      enddo
      dos_sum=mpreald(0.d0)
      do i=1,nhist
        dos_sum = dos_sum + exp(mpreald(dos_hist(i,ihist)))
      enddo
      ! cannot use standard fortran formatters to write to output.
      ! must use 'mpwrite' from the mpfun package to output high
      ! precision numbers.
      ! mpwrite(output, number of digits, number of digits after
      ! decimal, mpreald)
      write(6, "(i7, 2x)", advance="no")i
      call mpwrite(6,20,7,dos_sum)
      close(33)

      ! COMPUTE THE isotherm with the functions in wang_landau
      contains
      subroutine abort_csv_read(kode,nfile)

c***********************************************************************
c     
c     subroutine for aborting FIELD file read
c     
c***********************************************************************
      
      implicit none

      integer kode,nfile,idnode

      close (nfile)
      if(kode.eq.1)then
c     end of field file error exit
        write(*,*)"End of file error"
        call exit(1)
      else if(kode.eq.2)then
c     unrecognised directive in field file
        write(*,*)"Unrecognized read error"
        call exit(1)
      endif

      return
      end subroutine abort_csv_read
      end
