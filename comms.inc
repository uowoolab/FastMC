c*********************************************************************
c
c     Include file for MPI, PVM and SHMEM
c
c*********************************************************************

CMPIU define mpipriv mpipriv_
C$$$      integer Ilen,Dlen,Llen
C$$$      parameter (Ilen = 8, Dlen = 8, Llen = 8)

      include "mpif.h"

      integer Dlen
      parameter (Dlen = 8)

c     parameters for message tags

      integer Xgcopy_tag, Merge1_tag
      integer Dgsum_tag1, Dgsum_tag2
      integer Igmax_tag1, Igmax_tag2
      integer Igsum_tag1, Igsum_tag2
      integer Shmove_tag, Merge_tag
      integer Merge4_tag, Shell_tag
      integer Gstate_tag1, Gstate_tag2
      integer Passcon_tag

c     MPI tagsizes

      parameter(
     $     Xgcopy_tag = 105,
     $     Dgsum_tag1  = 6001,
     $     Dgsum_tag2  = 9001,
     $     Igsum_tag1  = 12001,
     $     Igsum_tag2  = 15001,
     $     Shmove_tag  = 18001,
     $     Merge_tag   = 6001,
     $     Merge1_tag  = 9001,
     $     Merge4_tag  = 12001,
     $     Gstate_tag1 = 15001,
     $     Gstate_tag2 = 18001,
     $     Shell_tag   = 21001,
     $     Passcon_tag =  3001,
     $     Igmax_tag1  = 24001,
     $     Igmax_tag2  = 27001)

