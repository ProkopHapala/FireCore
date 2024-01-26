program molguiapp
!use conf, only: natoms, nbonds, bonds, bond_orders, ufftypes, bond_orders_int
  implicit none
  character(20) :: tagcharge, rule_atom, conj
  character(200) :: input, output, output_qeq
!integer :: i
  
  read(*,*) input, output, output_qeq, tagcharge, rule_atom, conj
  if ( rule_atom /= 'simple' .and. rule_atom /= 'huckel' ) stop 'MAIN wrong rule for assigning atom types'
  if ( conj /= 'yes' .and. conj /= 'no' ) stop 'MAIN wrong switch for dealing with conjugation'
  
  ! read UFF parameters
  call read_uff_params

  ! read configuration
  call read_conf ( input )

  ! find neighbors
  call find_neigh

  ! set trivial cases
  call trivial_cases
  
  ! exception for nitro groups
  call special_case_nitro
  
  ! find a valid limit resonance structure
  call tree_walk

  if ( rule_atom == 'huckel' ) then
     ! find aromatic rings
     call find_rings
  else if ( rule_atom == 'simple' ) then
     ! assign resonant atoms
     call simple_rule
  end if
  
  ! set the rest of sp2 oxygens
  call assign_O2

  ! set the rest of sp3 oxygens
  call assign_O3

  ! set the rest of sp2 nitrogens
  call assign_N2

  ! set the rest of sp3 nitrogens
  call assign_N3

  ! set the rest of sp2 carbons
  call assign_C2

  ! set the rest of sp carbons
  call assign_C1

  ! check bond saturation
  call fix_saturation

  ! exception to avoid cumulenes
  call special_case_cumulene
  
  ! manually change sp3 nitrogen and oxygen to "resonant" when they are bonded to an sp2 atom (conjugation)
  if ( conj == 'yes' ) then
     call conjugation
  end if

  ! some check before assigning parameters
  call checks
  
  ! exception for amide groups
  call special_case_amide

  ! assign vdw parameters
  call assign_vdw
  
  ! assign bond parameters
  call assign_bonds
  
  ! assign angle parameters
  call assign_angles
  
  ! assign torsion parameters
  call assign_torsions
  
  ! assign inversion parameters
  call assign_inversions
  
  ! write data
  call write_data ( output, output_qeq, tagcharge )

!do i = 1, natoms
!   write(0,'(a,i2,a,a)') 'atom ', i, ' has type ', ufftypes(i)
!end do
!do i = 1, nbonds
!   write(0,'(a,i2,a,2i3,a,f7.1,i3)') 'bond ', i, ' between atoms ', bonds(i,1:2), ' has order ', bond_orders(i), bond_orders_int(i)
!end do
   
  stop
end program molguiapp
