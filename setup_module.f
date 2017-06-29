      module setup_module

c***********************************************************************
c     
c     module for defining default array sizes
c     
c     note the following internal units apply everywhere
c     
c     unit of time      (to)    =          1 x 10**(-12) seconds
c     unit of length    (lo)    =          1 x 10**(-10) metres
c     unit of mass      (mo)    = 1.6605402  x 10**(-27) kilograms
c     unit of charge    (qo)    = 1.60217733 x 10**(-19) coulombs
c     unit of energy    (eo)    = 1.6605402  x 10**(-23) joules
c     unit of pressure  (po)    = 1.6605402  x 10**(  7) pascals
c     
c*********************************************************************

      implicit none

c     array allocation parameters (set by subroutine parset)

      integer kmaxa,kmaxb,kmaxc,minnode,msatms,msbad,msgrp
      integer mspmf,msteth,mxangl,mxatms,mxbond,mxbuff,mxcell
      integer mxcons,mxdihd,mxewld,mxexcl,mxfbp,mxfld,mxgatm,mxgrid
      integer mxgrp,mxinv,mxlist,mxlshp,mxneut,mxngp,mxnstk,mxpang
      integer mxpbnd,mxpdih,mxpfbp,mxpinv,mxpmf,mxproc,mxptbp,mxpvdw
      integer mxrdf,mxshl,mxsite,mxspmf,mxstak,mxtang,mxtbnd
      integer mxtbp,mxtcon,mxtdih,mxteth,mxtinv,mxtmls,mxtshl,mxungp
      integer mxvdw,mxxdf,mx2tbp,mx3fbp,mxebuf,mxquat,mxshak,mxspl
      integer kmaxd,kmaxe,kmaxf,mxspme,mxftab,mxhko,mxhke
      integer mxmet,mxsmet,mxpmet,mxter,mxpter,mxatyp,mxxtyp

      save kmaxa,kmaxb,kmaxc,minnode,msatms,msbad,msgrp
      save mspmf,msteth,mxangl,mxatms,mxbond,mxbuff,mxcell
      save mxcons,mxdihd,mxewld,mxexcl,mxfbp,mxfld,mxgatm,mxgrid
      save mxgrp,mxinv,mxlist,mxlshp,mxneut,mxngp,mxnstk,mxpang
      save mxpbnd,mxpdih,mxpfbp,mxpinv,mxpmf,mxproc,mxptbp,mxpvdw
      save mxrdf,mxshl,mxsite,mxspmf,mxstak,mxtang,mxtbnd
      save mxtbp,mxtcon,mxtdih,mxteth,mxtinv,mxtmls,mxtshl,mxungp
      save mxvdw,mxxdf,mx2tbp,mx3fbp,mxebuf,mxquat,mxshak,mxspl
      save kmaxd,kmaxe,kmaxf,mxspme,mxftab,mxhko,mxhke
      save mxmet,mxsmet,mxpmet,mxter,mxpter,mxatyp,mxxtyp

      end module setup_module
