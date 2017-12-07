      subroutine error(idnode,iode)

c***********************************************************************
c     
c     subroutine for printing error messages and bringing
c     about a controlled termination of the program
c     
c***********************************************************************
      use utility_pack 

      implicit none

      logical kill
      integer idnode,iode,kode

      kill=(iode.ge.0)
      kode = abs(iode)
      
      if(idnode.eq.0) then
        
        if(kill) then
          if(kode.eq.0)then
            write(nrite,'(/,/,1x,a)') 
     x      'normal termination of FASTMC'
          else
            write(nrite,'(/,/,1x,a,i5)') 
     x      'terminated due to ', kode
          endif
        endif

        if (kode.eq. 0) then

c     dummy entry

        elseif (kode.eq. 1) then
          write(nrite,'(/,/,1x,a)')
     x      'error - PVM_NODES unset'        
        elseif (kode.eq. 2) then
          write(nrite,'(/,/,1x,a)')
     x      'error - machine not a hypercube'        
        elseif (kode.eq. 3) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unknown directive found in CONTROL file'        
        elseif (kode.eq. 4) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unknown directive found in FIELD file'
        elseif (kode.eq. 5) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unknown energy unit requested'
        elseif (kode.eq. 6) then
          write(nrite,'(/,/,1x,a)')
     x      'error - energy unit not specified'
        elseif (kode.eq. 7) then
          write(nrite,'(/,/,1x,a)')
     x      'error - energy unit respecified'
        elseif (kode.eq. 8) then
          write(nrite,'(/,/,1x,a)')
     x      'error - time step not specified'
        elseif (kode.eq.10) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many molecule types specified'
        elseif (kode.eq.11) then
          write(nrite,'(/,/,1x,a)')
     x      'error - duplicate molecule directive in FIELD file'
        elseif (kode.eq.12) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unknown molecule directive in FIELD file'
        elseif (kode.eq.13) then
          write(nrite,'(/,/,1x,a)')
     x      'error - molecular species not yet specified'
        elseif (kode.eq.14) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many unique atom types specified'
        elseif (kode.eq.15) then
          write(nrite,'(/,/,1x,a)')
     x      'error - duplicate pair potential specified'
        elseif (kode.eq.16) then
          write(nrite,'(/,/,1x,a)')
     x      'error - strange exit from FIELD file processing'
        elseif (kode.eq.17) then
          write(nrite,'(/,/,1x,a)')
     x      'error - strange exit from CONTROL file processing'
        elseif (kode.eq.18) then
          write(nrite,'(/,/,1x,a)')
     x      'error - duplicate 3-body potential specified'
        elseif (kode.eq.19) then
          write(nrite,'(/,/,1x,a)')
     x      'error - duplicate 4-body potential specified'
        elseif (kode.eq.20) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many molecule sites specified'
        elseif (kode.eq.21) then
          write(nrite,'(/,/,1x,a)')
     x      'error - duplicate tersoff potential specified'
        elseif (kode.eq.22) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unsuitable radial increment in TABLE file'
        elseif (kode.eq.23) then
          write(nrite,'(/,/,1x,a)')
     x      'error - incompatible FIELD and TABLE file potentials'
        elseif (kode.eq.24) then
          write(nrite,'(/,/,1x,a)')
     x      'error - end of file encountered in TABLE file'
        elseif (kode.eq.25) then
          write(nrite,'(/,/,1x,a)')
     x      'error - wrong atom type found in CONFIG file'
        elseif (kode.eq.26) then
          write(nrite,'(/,/,1x,a)')
     x      'error - cutoff smaller than EAM potential range'
        elseif (kode.eq.27) then
          write(nrite,'(/,/,1x,a)')
     x      'error - incompatible FIELD and TABEAM file potentials'
        elseif (kode.eq.28) then
          write(nrite,'(/,/,1x,a)')
     x      'error - transfer buffer too small in mettab'
        elseif (kode.eq.29) then
          write(nrite,'(/,/,1x,a)')
     x      'error - end of file encountered in TABEAM file'
        elseif (kode.eq.30) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many chemical bonds specified'
        elseif (kode.eq.31) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many chemical bonds in system'
        elseif (kode.eq.32) then
          write(nrite,'(/,/,1x,a)')
     x      'error - integer array memory allocation failure'
        elseif (kode.eq.33) then
          write(nrite,'(/,/,1x,a)')
     x      'error - real array memory allocation failure'
        elseif (kode.eq.34) then
          write(nrite,'(/,/,1x,a)')
     x      'error - character array memory allocation failure'
        elseif (kode.eq.35) then
          write(nrite,'(/,/,1x,a)')
     x      'error - logical array  memory allocation failure'
        elseif (kode.eq.36) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed fmet array allocation in mettab'
        elseif (kode.eq.40) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many bond constraints specified'
        elseif (kode.eq.41) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many bond constraints in system'
        elseif (kode.eq.42) then
          write(nrite,'(/,/,1x,a)')
     x      'error - transfer buffer too small in merge1'
        elseif (kode.eq.45) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many atoms in CONFIG file'
        elseif (kode.eq.46) then
          write(nrite,'(/,/,1x,a)')
     x      'error - ewlbuf array too small in ewald1'
        elseif (kode.eq.47) then
          write(nrite,'(/,/,1x,a)')
     x      'error - transfer buffer too small in merge'
        elseif (kode.eq.48) then
          write(nrite,'(/,/,1x,a)')
     x      'error - transfer buffer too small in fortab'
        elseif (kode.eq.49) then
          write(nrite,'(/,/,1x,a)')
     x      'error - frozen core-shell unit specified'
        elseif (kode.eq.50) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many bond angles specified'
        elseif (kode.eq.51) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many bond angles in system'
        elseif (kode.eq.52) then
          write(nrite,'(/,/,1x,a)')
     x      'error - end of FIELD file encountered'
        elseif (kode.eq.53) then
          write(nrite,'(/,/,1x,a)')
     x      'error - end of CONTROL file encountered'
        elseif (kode.eq.54) then
          write(nrite,'(/,/,1x,a)')
     x      'error - problem reading CONFIG file'
        elseif (kode.eq.55) then
          write(nrite,'(/,/,1x,a)')
     x      'error - end of CONFIG file encountered'
        elseif (kode.eq.57) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many core-shell units specified'
        elseif (kode.eq.59) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many core-shell units in system'
        elseif (kode.eq.60) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many dihedral angles specified'
        elseif (kode.eq.61) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many dihedral angles in system'
        elseif (kode.eq.62) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many tethered atoms specified'
        elseif (kode.eq.63) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many tethered atoms in system'
        elseif (kode.eq.65) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many excluded pairs specified'
        elseif (kode.eq.66) then
          write(nrite,'(/,/,1x,a)')
     x      'error - incorrect boundary condition for HK ewald'
        elseif (kode.eq.67) then
          write(nrite,'(/,/,1x,a)')
     x      'error - incorrect boundary condition in thbfrc'
        elseif (kode.eq.69) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many link cells required in thbfrc'
        elseif (kode.eq.70) then
          write(nrite,'(/,/,1x,a)')
     x      'error - constraint bond quench failure'
        elseif (kode.eq.71) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many metal potentials specified'
        elseif (kode.eq.72) then
          write(nrite,'(/,/,1x,a)')
     x      'error - different metal potential types specified'
        elseif (kode.eq.73) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many inversion potentials specified'
        elseif (kode.eq.75) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many atoms in specified system'
        elseif (kode.eq.77) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many inversion potentials in system'
        elseif (kode.eq.79) then
          write(nrite,'(/,/,1x,a)')
     x      'error - incorrect boundary condition in fbpfrc'
        elseif (kode.eq.80) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many pair potentials specified'
        elseif (kode.eq.81) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unidentified atom in pair potential list'
        elseif (kode.eq.82) then
          write(nrite,'(/,/,1x,a)')
     x      'error - calculated pair potential index too large'
        elseif (kode.eq.83) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many three body potentials specified'
        elseif (kode.eq.84) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unidentified atom in 3-body potential list'
        elseif (kode.eq.85) then
          write(nrite,'(/,/,1x,a)')
     x      'error - required velocities not in CONFIG file'
        elseif (kode.eq.86) then
          write(nrite,'(/,/,1x,a)')
     x      'error - calculated 3-body potential index too large'
        elseif (kode.eq.87) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many link cells required in fbpfrc'
        elseif (kode.eq.88) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many tersoff potentials specified'
        elseif (kode.eq.89) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many four body potentials specified'
        elseif (kode.eq.90) then
          write(nrite,'(/,/,1x,a)')
     x      'error - system total electric charge nonzero'
        elseif (kode.eq.91) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unidentified atom in 4-body potential list'
        elseif (kode.eq.92) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unidentified atom in tersoff potential list'
        elseif (kode.eq.93) then
          write(nrite,'(/,/,1x,a)')
     x      'error - cannot use shell model with rigid molecules'
        elseif (kode.eq.95) then
          write(nrite,'(/,/,1x,a)')
     x      'error - potential cutoff exceeds half-cell width'
        elseif (kode.eq.97) then
          write(nrite,'(/,/,1x,a)')
     x      'error - cannot use shell model with neutral groups'
        elseif (kode.eq.99) then
          write(nrite,'(/,/,1x,a)')
     x      'error - cannot use shell model with constraints'
        elseif (kode.eq.100) then
          write(nrite,'(/,/,1x,a)')
     x      'error - forces working arrays too small'
        elseif (kode.eq.101) then
          write(nrite,'(/,/,1x,a)')
     x      'error - calculated 4-body potential index too large'
        elseif (kode.eq.102) then
          write(nrite,'(/,/,1x,a)')
     x      'error - parameter mxproc exceeded in shake arrays'
        elseif (kode.eq.103) then
          write(nrite,'(/,/,1x,a)')
     x      'error - parameter mxlshp exceeded in shake arrays'
        elseif (kode.eq.105) then
          write(nrite,'(/,/,1x,a)')
     x      'error - shake algorithm failed to converge'
        elseif (kode.eq.106) then
          write(nrite,'(/,/,1x,a,a)')
     x      'error - neighbour list array too small in parlink ',
     x      'subroutine'
        elseif (kode.eq.107) then
          write(nrite,'(/,/,1x,a,a)')
     x      'error - neighbour list array too small in parlinkneu ',
     x      'subroutine'
        elseif (kode.eq.108) then
          write(nrite,'(/,/,1x,a,a)')
     x      'error - neighbour list array too small in parneulst ',
     x      'subroutine'
        elseif (kode.eq.109) then
          write(nrite,'(/,/,1x,a,a)')
     x      'error - neighbour list array too small in parlst_nsq ',
     x      'subroutine'
        elseif (kode.eq.110) then
          write(nrite,'(/,/,1x,a,a)')
     x      'error - neighbour list array too small in parlst ',
     x      'subroutine'
        elseif (kode.eq.112) then
          write(nrite,'(/,/,1x,a)')
     x      'error - vertest array too small'
        elseif (kode.eq.120) then
          write(nrite,'(/,/,1x,a)')
     x      'error - invalid determinant in matrix inversion'
        elseif (kode.eq.130) then
          write(nrite,'(/,/,1x,a)')
     x      'error - incorrect octahedral boundary condition'
        elseif (kode.eq.135) then
          write(nrite,'(/,/,1x,a)')
     x      'error - incorrect hexagonal prism boundary condition'
        elseif (kode.eq.140) then
          write(nrite,'(/,/,1x,a)')
     x      'error - incorrect dodecahedral boundary condition'
        elseif (kode.eq.141) then
          write(nrite,'(/,/,1x,a)')
     x      'error - duplicate metal potential specified'
        elseif (kode.eq.145) then
          write(nrite,'(/,/,1x,a)')
     x      'error - no van der waals potentials defined'
        elseif (kode.eq.150) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unknown van der waals potential selected'
        elseif (kode.eq.151) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unknown metal potential selected'
        elseif (kode.eq.153) then
          write(nrite,'(/,/,1x,a)')
     x      'error - metals not permitted with multiple timestep'
        elseif (kode.eq.160) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unaccounted for atoms in exclude list '
        elseif (kode.eq.170) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many variables for statistic array '
        elseif (kode.eq.180) then
          write(nrite,'(/,/,1x,a)')
     x      'error - Ewald sum requested in non-periodic system'
        elseif (kode.eq.185) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many reciprocal space vectors'
        elseif (kode.eq.186) then
          write(nrite,'(/,/,1x,a)')
     x      'error - transfer buffer array too small in sysgen'
        elseif (kode.eq.190) then
          write(nrite,'(/,/,1x,a)')
     x      'error - buffer array too small in splice'
        elseif (kode.eq.200) then
          write(nrite,'(/,/,1x,a)')
     x      'error - rdf buffer array too small in revive'
        elseif (kode.eq.220) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many neutral groups in system'
        elseif (kode.eq.230) then
          write(nrite,'(/,/,1x,a)')
     x      'error - neutral groups improperly arranged'
        elseif (kode.eq.250) then
          write(nrite,'(/,/,1x,a)')
     x      'error - Ewald sum requested with neutral groups'
        elseif (kode.eq.260) then
          write(nrite,'(/,/,1x,a)')
     x      'error - parameter mxexcl exceeded in excludeneu routine'
        elseif (kode.eq.300) then
          write(nrite,'(/,/,1x,a)')
     x      'error - incorrect boundary condition in parlink'
        elseif (kode.eq.301) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many rigid body types '
        elseif (kode.eq.302) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many sites in rigid body '
        elseif (kode.eq.303) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many rigid bodies specified'
        elseif (kode.eq.304) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many rigid body sites in system '
        elseif (kode.eq.305) then
          write(nrite,'(/,/,1x,a)')
     x      'error - box size too small for link cells'
        elseif (kode.eq.306) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed to find principal axis system'
        elseif (kode.eq.310) then
          write(nrite,'(/,/,1x,a)')
     x      'error - quaternion setup failed '
        elseif (kode.eq.320) then
          write(nrite,'(/,/,1x,a)')
     x      'error - site in multiple rigid bodies'
        elseif (kode.eq.321) then
          write(nrite,'(/,/,1x,a)')
     x      'error - quaternion integrator failed'
        elseif (kode.eq.330) then
          write(nrite,'(/,/,1x,a)')
     x      'error - mxewld parameter incorrect'
        elseif (kode.eq.331) then
          write(nrite,'(/,/,1x,a)')
     x      'error - mxhke parameter incorrect'
        elseif (kode.eq.332) then
          write(nrite,'(/,/,1x,a)')
     x      'error - mxhko parameter too small'
        elseif (kode.eq.340) then
          write(nrite,'(/,/,1x,a)')
     x      'error - invalid integration option requested'
        elseif (kode.eq.350) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too few degrees of freedom'
        elseif (kode.eq.360) then
          write(nrite,'(/,/,1x,a)')
     x      'error - frozen atom found in rigid body'
        elseif (kode.eq.380) then
          write(nrite,'(/,/,1x,a)')
     x      'error - simulation temperature not specified'
        elseif (kode.eq.381) then
          write(nrite,'(/,/,1x,a)')
     x      'error - simulation timestep not specified'
        elseif (kode.eq.382) then
          write(nrite,'(/,/,1x,a)')
     x      'error - simulation cutoff not specified'
        elseif (kode.eq.383) then
          write(nrite,'(/,/,1x,a)')
     x      'error - simulation forces option not specified'
        elseif (kode.eq.384) then
          write(nrite,'(/,/,1x,a)')
     x      'error - verlet strip width not specified'
        elseif (kode.eq.385) then
          write(nrite,'(/,/,1x,a)')
     x      'error - primary cutoff not specified'
        elseif (kode.eq.386) then
          write(nrite,'(/,/,1x,a)')
     x      'error - primary cutoff larger than rcut'
        elseif (kode.eq.387) then
          write(nrite,'(/,/,1x,a)')
     x      'error - guest pressure not specified'
        elseif (kode.eq.388) then
          write(nrite,'(/,/,1x,a)')
     x      'error - npt incompatible with multiple timestep'
        elseif (kode.eq.389) then
          write(nrite,'(/,/,1x,a)')
     x      'number of pimd beads not specified in field file'
        elseif (kode.eq.390) then
          write(nrite,'(/,/,1x,a)')
     x      'error - npt ensemble requested in non-periodic system'
        elseif (kode.eq.391) then
          write(nrite,'(/,/,1x,a)')
     x      'error - incorrect number of pimd beads in config file'
        elseif (kode.eq.392) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many link cells requested'
        elseif (kode.eq.394) then
          write(nrite,'(/,/,1x,a)')
     x      'error - minimun image arrays exceeded'
        elseif (kode.eq.396) then
          write(nrite,'(/,/,1x,a)')
     x      'error - interpolation array exceeded'
        elseif (kode.eq.398) then
          write(nrite,'(/,/,1x,a)')
     x      'error - cutoff too small for rprim and delr'
        elseif (kode.eq.400) then
          write(nrite,'(/,/,1x,a)')
     x      'error - rvdw greater than cutoff'
        elseif (kode.eq.402) then
          write(nrite,'(/,/,1x,a)')
     x      'error - van der waals cutoff unset'
        elseif (kode.eq.410) then
          write(nrite,'(/,/,1x,a)')
     x      'error - cell not consistent with image convention'
        elseif (kode.eq.412) then
          write(nrite,'(/,/,1x,a)')
     x      'error - mxxdf parameter too small for shake routine'
        elseif (kode.eq.414) then
          write(nrite,'(/,/,1x,a)')
     x      'error - conflicting ensemble options in CONTROL file'
        elseif (kode.eq.416) then
          write(nrite,'(/,/,1x,a)')
     x      'error - conflicting force options in CONTROL file'
        elseif (kode.eq.418) then
          write(nrite,'(/,/,1x,a)')
     x      'error - bond vector work arrays too small in bndfrc'
        elseif (kode.eq.419) then
          write(nrite,'(/,/,1x,a)')
     x      'error - bond vector work arrays too small in angfrc'
        elseif (kode.eq.420) then
          write(nrite,'(/,/,1x,a)')
     x      'error - bond vector work arrays too small in tethfrc'
        elseif (kode.eq.421) then
          write(nrite,'(/,/,1x,a)')
     x      'error - bond vector work arrays too small in dihfrc'
        elseif (kode.eq.422) then
          write(nrite,'(/,/,1x,a)')
     x      'error - all-pairs must use multiple timestep'
        elseif (kode.eq.423) then
          write(nrite,'(/,/,1x,a)')
     x      'error - bond vector work arrays too small in shlfrc'
        elseif (kode.eq.424) then
          write(nrite,'(/,/,1x,a)')
     x      'error - electrostatics incorrect for all-pairs'
        elseif (kode.eq.425) then
          write(nrite,'(/,/,1x,a)')
     x      'error - transfer buffer array too small in shlmerge'
        elseif (kode.eq.426) then
          write(nrite,'(/,/,1x,a)')
     x      'error - neutral groups not permitted with all-pairs'
        elseif (kode.eq.427) then
          write(nrite,'(/,/,1x,a)')
     x      'error - bond vector work arrays too small in invfrc'
        elseif (kode.eq.430) then
          write(nrite,'(/,/,1x,a)')
     x      'error - integration routine not available'
        elseif (kode.eq.432) then
          write(nrite,'(/,/,1x,a)')
     x      'error - intlist failed to assign constraints '
        elseif (kode.eq.433) then
          write(nrite,'(/,/,1x,a)')
     x      'error - specify rcut before the Ewald sum precision'
        elseif (kode.eq.434) then
          write(nrite,'(/,/,1x,a)')
     x      'error - illegal entry into STRESS related routine'
        elseif (kode.eq.436) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unrecognised ensemble '
        elseif (kode.eq.438) then
          write(nrite,'(/,/,1x,a)')
     x      'error - PMF constraints failed to converge'
        elseif (kode.eq.440) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined angular potential'
        elseif (kode.eq.442) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined three body potential'
        elseif (kode.eq.443) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined four body potential'
        elseif (kode.eq.444) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined bond potential'
        elseif (kode.eq.445) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined many body potential'
        elseif (kode.eq.446) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined electrostatic key in dihfrc'
        elseif (kode.eq.447) then
          write(nrite,'(/,/,1x,a)')
     x      'error - 1-4 separation exceeds cutoff range'
        elseif (kode.eq.448) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined dihedral potential'
        elseif (kode.eq.449) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined inversion potential'
        elseif (kode.eq.450) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined tethering potential'
        elseif (kode.eq.451) then
          write(nrite,'(/,/,1x,a)')
     x      'error - three body potential cutoff undefined'
        elseif (kode.eq.452) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined pair potential'
        elseif (kode.eq.453) then
          write(nrite,'(/,/,1x,a)')
     x      'error - four body potential cutoff undefined'
        elseif (kode.eq.454) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined external field'
        elseif (kode.eq.456) then
          write(nrite,'(/,/,1x,a)')
     x      'error - core and shell in same rigid unit'
        elseif (kode.eq.458) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many PMF constraints - param. mspmf too small'
        elseif (kode.eq.460) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many PMF sites - parameter mxspmf too small'
        elseif (kode.eq.461) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined metal potential'
        elseif (kode.eq.462) then
          write(nrite,'(/,/,1x,a)')
     x      'error - PMF UNIT  record expected'
        elseif (kode.eq.463) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unidentified atom in metal potential list'
        elseif (kode.eq.464) then
          write(nrite,'(/,/,1x,a)')
     x      'error - thermostat time constant must be > 0.d0'
        elseif (kode.eq.465) then
          write(nrite,'(/,/,1x,a)')
     x      'error - calculated pair potential index too large'
        elseif (kode.eq.466) then
          write(nrite,'(/,/,1x,a)')
     x      'error - barostat time constant must be > 0.d0'
        elseif (kode.eq.468) then
          write(nrite,'(/,/,1x,a)')
     x      'error - r0 too large for snm potential with current cutoff'
        elseif (kode.eq.470) then
          write(nrite,'(/,/,1x,a)')
     x      'error - n<m in definition of n-m potential'
        elseif (kode.eq.474) then
          write(nrite,'(/,/,1x,a)')
     x      'error - mxxdf too small in parlst subroutine'
        elseif (kode.eq.475) then
          write(nrite,'(/,/,1x,a)')
     x      'error - mxxdf too small in parlst_nsq subroutine'
        elseif (kode.eq.476) then
          write(nrite,'(/,/,1x,a)')
     x      'error - mxxdf too small in parneulst subroutine'
        elseif (kode.eq.477) then
          write(nrite,'(/,/,1x,a)')
     x      'error - mxxdf too small in prneulst subroutine'
        elseif (kode.eq.478) then
          write(nrite,'(/,/,1x,a)')
     x      'error - mxxdf too small in forcesneu subroutine'
        elseif (kode.eq.479) then
          write(nrite,'(/,/,1x,a)')
     x      'error - mxxdf too small in multipleneu subroutine'
        elseif (kode.eq.484) then
          write(nrite,'(/,/,1x,a)')
     x      'error - only one potential of mean force permitted'
        elseif (kode.eq.486) then
          write(nrite,'(/,/,1x,a)')
     x      'error - HK real space screening function cutoff violation'
        elseif (kode.eq.487) then
          write(nrite,'(/,/,1x,a)')
     x      'error - HK recip space screening function cutoff violation'
        elseif (kode.eq.488) then
          write(nrite,'(/,/,1x,a)')
     x      'error - HK lattice control parameter set too large'
        elseif (kode.eq.490) then
          write(nrite,'(/,/,1x,a)')
     x      'error - PMF parameter mxpmf too small in passpmf'
        elseif (kode.eq.492) then
          write(nrite,'(/,/,1x,a)')
     x      'error - parameter mxcons < number of PMF constraints'
        elseif (kode.eq.494) then
          write(nrite,'(/,/,1x,a)')
     x      'error in csend: pvmfinitsend'
        elseif (kode.eq.496) then
          write(nrite,'(/,/,1x,a)')
     x      'error in csend: pvmfpack'
        elseif (kode.eq.498) then
          write(nrite,'(/,/,1x,a)')
     x      'error in csend: pvmfsend'
        elseif (kode.eq.500) then
          write(nrite,'(/,/,1x,a)')
     x      'error in crecv: pvmfrecv'
        elseif (kode.eq.502) then
          write(nrite,'(/,/,1x,a)')
     x      'error in crecv: pvmfunpack'
        elseif (kode.eq.504) then
          write(nrite,'(/,/,1x,a)')
     x      'error - cutoff too large for TABLE file'
        elseif (kode.eq.506) then
          write(nrite,'(/,/,1x,a)')
     x      'error - work arrays too small for quaternion integration'
        elseif (kode.eq.508) then
          write(nrite,'(/,/,1x,a)')
     x      'error - rigid bodies not permitted with RESPA algorithm'
        elseif (kode.eq.510) then
          write(nrite,'(/,/,1x,a)')
     x      'error - structure optimiser not permitted with RESPA'
        elseif (kode.eq.513) then
          write(nrite,'(/,/,1x,a)')
     x      'error - SPME not available for given boundary conditions'
        elseif (kode.eq.514) then
          write(nrite,'(/,/,1x,a)')
     x      'error - SPME routines have not been compiled in'
        elseif (kode.eq.1000) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of configuration arrays'
        elseif (kode.eq.1001) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of probability arrays'
        elseif (kode.eq.1002) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of van der Waals arrays'
        elseif (kode.eq.1010) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of angle arrays'
        elseif (kode.eq.1011) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of dihedral arrays'
        elseif (kode.eq.1012) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of exclude arrays'
        elseif (kode.eq.1013) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of rigid body arrays'
        elseif (kode.eq.1014) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of vdw arrays'
        elseif (kode.eq.1020) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of angle work arrays'
        elseif (kode.eq.1030) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of bond arrays'
        elseif (kode.eq.1040) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of bond work arrays'
        elseif (kode.eq.1050) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of dihedral arrays'
        elseif (kode.eq.1060) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of dihedral work arrays'
        elseif (kode.eq.1070) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of constraint arrays'
        elseif (kode.eq.1090) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of site arrays'
        elseif (kode.eq.1100) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of core_shell arrays'
        elseif (kode.eq.1110) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of core_shell work arrays'
        elseif (kode.eq.1120) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of inversion arrays'
        elseif (kode.eq.1130) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of inversion work arrays' 
        elseif (kode.eq.1140) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of four-body arrays'
        elseif (kode.eq.1150) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of four-body work arrays'
        elseif (kode.eq.1170) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of three-body arrays'
        elseif (kode.eq.1180) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of three-body work arrays'
        elseif (kode.eq.1200) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of external field arrays'
        elseif (kode.eq.1210) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of pmf arrays'
        elseif (kode.eq.1220) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of pmf_lf or pmf_vv work arrays'
        elseif (kode.eq.1230) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of pmf_shake work arrays'
        elseif (kode.eq.1240) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of ewald arrays'
        elseif (kode.eq.1250) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of excluded atom arrays'
        elseif (kode.eq.1260) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of tethering arrays'
        elseif (kode.eq.1270) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of tethering work arrays'
        elseif (kode.eq.1280) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of metal arrays'
        elseif (kode.eq.1290) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nvt_h0.f'
        elseif (kode.eq.1300) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of dens0 array in npt_b0.f'
        elseif (kode.eq.1310) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in npt_b0.f'
        elseif (kode.eq.1320) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of dens0 array in npt_h0.f'
        elseif (kode.eq.1330) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in npt_h0.f'
        elseif (kode.eq.1340) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of dens0 array in nst_b0.f'
        elseif (kode.eq.1350) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nst_b0.f'
        elseif (kode.eq.1360) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of dens0 array in nst_h0.f'
        elseif (kode.eq.1370) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nst_h0.f'
        elseif (kode.eq.1380) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nve_1.f'
        elseif (kode.eq.1390) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nvt_e1.f'
        elseif (kode.eq.1400) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nvt_b1.f'
        elseif (kode.eq.1410) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nvt_h1.f'
        elseif (kode.eq.1420) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in npt_b1.f'
        elseif (kode.eq.1430) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in npt_b1.f'
        elseif (kode.eq.1440) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in npt_h1.f'
        elseif (kode.eq.1450) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in npt_h1.f'
        elseif (kode.eq.1460) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nst_b1.f'
        elseif (kode.eq.1470) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in nst_b1.f'
        elseif (kode.eq.1480) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nst_h1.f'
        elseif (kode.eq.1490) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in nst_h1.f'
        elseif (kode.eq.1500) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nveq_1.f'
        elseif (kode.eq.1510) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nvtq_b1.f'
        elseif (kode.eq.1520) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nvtq_h1.f'
        elseif (kode.eq.1530) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nptq_b1.f'
        elseif (kode.eq.1540) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in nptq_b1.f'
        elseif (kode.eq.1550) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nptq_h1.f'
        elseif (kode.eq.1560) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in nptq_h1.f'
        elseif (kode.eq.1570) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nstq_b1.f'
        elseif (kode.eq.1580) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in nstq_b1.f'
        elseif (kode.eq.1590) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nstq_h1.f'
        elseif (kode.eq.1600) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in nstq_h1.f'
        elseif (kode.eq.1610) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in qshake.f'
        elseif (kode.eq.1615) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in qrattle_q.f'
        elseif (kode.eq.1620) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nveq_2.f'
        elseif (kode.eq.1625) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in qrattle_v.f'
        elseif (kode.eq.1630) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nvtq_b2.f'
        elseif (kode.eq.1640) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nvtq_h2.f'
        elseif (kode.eq.1650) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nptq_b2.f'
        elseif (kode.eq.1660) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in nptq_b2.f'
        elseif (kode.eq.1670) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nptq_h2.f'
        elseif (kode.eq.1680) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in nptq_h2.f'
        elseif (kode.eq.1690) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nstq_b2.f'
        elseif (kode.eq.1700) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in nstq_b2.f'
        elseif (kode.eq.1710) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of work arrays in nstq_h2.f'
        elseif (kode.eq.1720) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of density array in nstq_h2.f'
        elseif (kode.eq.1730) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of HK Ewald arrays'
        elseif (kode.eq.1740) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of property arrays'
        elseif (kode.eq.1750) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of spme arrays'
        elseif (kode.eq.1760) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of ewald_spme.f work arrays'
        elseif (kode.eq.1770) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of quench.f work arrays'
        elseif (kode.eq.1780) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of quatqnch.f work arrays'
        elseif (kode.eq.1790) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of quatbook.f work arrays'
        elseif (kode.eq.1800) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of intlist.f work arrays'
        elseif (kode.eq.1810) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of forces.f work arrays'
        elseif (kode.eq.1820) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of forcesneu.f work arrays'
        elseif (kode.eq.1830) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of neutlst.f work arrays'
        elseif (kode.eq.1840) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of multiple.f work arrays'
        elseif (kode.eq.1850) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of multipleneu.f work arrays'
        elseif (kode.eq.1860) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of multiple_nsq.f work arrays'
        elseif (kode.eq.1870) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of parlst_nsq.f work arrays'
        elseif (kode.eq.1880) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of parlst.f work arrays'
        elseif (kode.eq.1890) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of parlink.f work arrays'
        elseif (kode.eq.1900) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of parlinkneu.f work arrays'
        elseif (kode.eq.1910) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of parneulst.f work arrays'
        elseif (kode.eq.1920) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of strucopt.f work arrays'
        elseif (kode.eq.1930) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of vertest.f work arrays'
        elseif (kode.eq.1940) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of pair arrays'
        elseif (kode.eq.1945) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of tersoff arrays'
        elseif (kode.eq.1950) then
          write(nrite,'(/,/,1x,a)')
     x      'error - shell relaxation cycle limit exceeded'
        elseif (kode.eq.1953) then
          write(nrite,'(/,/,1x,a)')
     x      'error - tersoff radius of cutoff not defined'
        elseif (kode.eq.1955) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of tersoff work arrays'
        elseif (kode.eq.1960) then
          write(nrite,'(/,/,1x,a)')
     x      'error - conflicting shell option in FIELD file'
        elseif (kode.eq.1970) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of shell_relax work arrays'
        elseif (kode.eq.1972) then
          write(nrite,'(/,/,1x,a)')
     x      'error - unknown tersoff potential defined in FIELD file'
        elseif (kode.eq.1974) then
          write(nrite,'(/,/,1x,a)')
     x      'error - incorrect period boundary in tersoff.f'
        elseif (kode.eq.1976) then
          write(nrite,'(/,/,1x,a)')
     x      'error - too many link cells required in tersoff.f'
        elseif (kode.eq.1978) then
          write(nrite,'(/,/,1x,a)')
     x      'error - undefined potential in tersoff.f'
        elseif (kode.eq.1980) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nvevv_1.f work arrays'
        elseif (kode.eq.1990) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nvtvv_b1.f work arrays'
        elseif (kode.eq.2000) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nvtvv_e1.f work arrays'
        elseif (kode.eq.2010) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nvtvv_h1.f work arrays'
        elseif (kode.eq.2020) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptvv_b1.f dens0 array'
        elseif (kode.eq.2030) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptvv_b1.f work arrays'
        elseif (kode.eq.2040) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptvv_h1.f dens0 array'
        elseif (kode.eq.2050) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptvv_h1.f work arrays'
        elseif (kode.eq.2060) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstvv_b1.f dens0 array'
        elseif (kode.eq.2070) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstvv_b1.f work arrays'
        elseif (kode.eq.2080) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstvv_h1.f dens0 array'
        elseif (kode.eq.2090) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstvv_b1.f work arrays'
        elseif (kode.eq.2100) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nveqvv_1.f work arrays'
        elseif (kode.eq.2110) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nveqvv_2.f work arrays'
        elseif (kode.eq.2120) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nvtqvv_b1.f work arrays'
        elseif (kode.eq.2130) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nvtqvv_b2.f work arrays'
        elseif (kode.eq.2140) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nvtqvv_h1.f work arrays'
        elseif (kode.eq.2150) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nvtqvv_h2.f work arrays'
        elseif (kode.eq.2160) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptqvv_b1.f dens0 array'
        elseif (kode.eq.2170) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptqvv_b1.f work arrays'
        elseif (kode.eq.2180) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptqvv_b2.f dens0 array'
        elseif (kode.eq.2190) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptqvv_b2.f work arrays'
        elseif (kode.eq.2200) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptqvv_h1.f dens0 array'
        elseif (kode.eq.2210) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptqvv_h1.f work arrays'
        elseif (kode.eq.2220) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptqvv_h2.f dens0 array'
        elseif (kode.eq.2230) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nptqvv_h2.f work arrays'
        elseif (kode.eq.2240) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstqvv_b1.f dens0 array'
        elseif (kode.eq.2250) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstqvv_b1.f work arrays'
        elseif (kode.eq.2260) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstqvv_b2.f dens0 array'
        elseif (kode.eq.2270) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstqvv_b2.f work arrays'
        elseif (kode.eq.2280) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstqvv_h1.f dens0 array'
        elseif (kode.eq.2290) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstqvv_h1.f work arrays'
        elseif (kode.eq.2300) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstqvv_h2.f dens0 array'
        elseif (kode.eq.2310) then
          write(nrite,'(/,/,1x,a)')
     x      'error - failed allocation of nstqvv_h2.f work arrays'
        elseif (kode.eq.2311) then
          write(nrite,'(/,/,1x,a)')
     x      'error - GCMC sim needs ewald specified in CONTROL file.'
        elseif (kode.eq.2312) then
          write(nrite,'(/,/,a70,/,14x,a,/,a70,/,/)')
     &'*****************************************************************
     &********','premature termination of GCMC simulation.','**********
     &****************************************************************'

        elseif(kode.eq.2313)then
          write(nrite,'(/,/,1x,a)')
     &     'error - number of guests specified in the CONTROL file 
     &does not match number of guests in the FIELD file'
        elseif (kode.eq.2314) then
          write(nrite,'(/,/,1x,a,/,1x,a)')
     x      'error - No pressure information found in the CONTROL file',
     &'Either set individual guest partial pressures, or a total 
     &pressure and mole fractions for the guests'
        elseif (kode.eq.2315) then
          write(nrite,'(/,/,1x,a)')
     x      'error - The sum of mole fractions of your guests sums to 
     & greater than 1.'
        elseif (kode.eq.2316) then
          write(nrite, '(/,/,/,1x,a)')
     &      'error - General termination from Wang-Landau routine'
        elseif (kode.eq.2317) then
          write(nrite,'(/,/,1x,a)')
     &      'error - failed allocation of Wang-Landau arrays'
        elseif (kode.eq.2318) then
          write(nrite,'(/,/,1x,a,a)')
     &      'error - too many guests specified in input. Wang-Landau',
     &      ' currently only works with 1 guest'
        elseif (kode.eq.2319) then
          write(nrite,'(/,/,1x,a,a)')
     &      'error - the wl_precision coefficient set in the CONTROL',
     &      ' file must be greater than 1'
        elseif (kode.eq.2320) then
          write(nrite,'(/,/,1x,a,/,4x,a,/,4x,a)')
     &      'error - could not insert the number of guests requested.',
     &      'Try changing the maxfactor value in utility_pack.f to a',
     &      'larger value and recompiling.'
        else
          write(nrite,'(/,/,1x,a)')
     x      'error - unnamed error found'
        endif
        
      endif

      if (kill) then

c     close all i/o channels
        
        if(idnode.eq.0) then
          close (nrite)
          close (ncontrol)
          close (nconfig)
          close (nstats)
          close (nfield)
          close (401)
          close (202)
c          close (203)
          close (501)
          close (205)
        endif
        
        call gsync()
        call exitcomms()

      endif
      
      return
      end
