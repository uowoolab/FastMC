
c*************************************************************
c      DOS_TO_ISOTHERM                                       *
c*************************************************************
c      program to read in a .csv file containing the DOS     *
c      and output an isotherm.                               *
c*************************************************************
      use utility_pack
      use readinputs
      use vdw_module
      use flex_module
      use ewald_module
      use mc_moves
      use wang_landau
      use parse_module
      use mpmodule

      implicit none
      logical loop,safe
      integer idnode,idum
      type(mp_real) z_uvt,n,dos_sum,temp1,ktmp,dosmp,exp1 

      idnode=0
      safe=.true.
      loop=.true.
      ! open the csv file
      open(33,file='final_dos.csv',status='old')
      ! read in the DOS
      call getrec(safe,idnode,33)
      if(.not.safe)call abort_csv_read(1,33)
      do while(loop)
        call lowcase(record,lenrec)
        call strip(record, lenrec)
        if(record(1).eq.'#'.or.record(1).eq.' ')then
        else
            N=intstr(record,lenrec,idum)
            dos=dblstr(record,lenrec,idum)
        endif
        call getrec(safe,idnode,33)
        if(.not.safe)loop=.false.
      enddo

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
      end subroutine abort_field_read
      end
