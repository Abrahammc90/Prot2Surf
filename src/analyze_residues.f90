!> Analyze which residues are close to the surface across encounters.
!!
!! This program reads a PDB for solute2 and an associations/complexes
!! file, transforms residue centers for each encounter and counts how
!! many encounters bring residues within a distance threshold from the
!! surface. Results are written to `closest_residues.txt`.
program main
  USE read_input
  USE maths


  implicit none
  ! This program analyzes which residues in `pdb2` are brought
  ! within a given distance threshold of the surface across a set
  ! of encounter complexes stored in a complexes/assoc file.
  !
  ! Workflow:
  ! 1. Parse command-line arguments (pdb2, complexes, nb_encounters, threshold)
  ! 2. Read PDB and complexes files into memory
  ! 3. Compute center-of-geometry (COG) for each residue in pdb2
  ! 4. For each encounter: transform residue COGs and count residues
  !    whose transformed Z coordinate is below `dist_threshold`
  ! 5. Write results to `closest_residues.txt`
      

    character*128 :: pdb2_filename, complexes_filename

    ! File/structure sizes and containers
    integer :: tot_atoms2, tot_residues2, tot_chains2, tot_residue_atoms2
    integer :: tot_encounters, nb_encounters
    type ( type_pdb_file ) :: pdb2
    type ( type_assoc_file ) :: complexes

    ! Loop indices and temporary storage
    integer :: i, j, k
    real (kind=8), dimension(3) :: residue_cog
    real (kind=8), dimension(:, :), allocatable :: residues_cog, residue_crds, new_residues_cog
    real (kind=8) :: dist_threshold
    integer, dimension(:), allocatable :: encounters_count
    integer, dimension(:, :), allocatable :: all_closest_resid
    integer :: closest_resid
    character (len=12) :: residue_str
    character (len=4) :: str_resid
    character (len=128) :: argument
    integer :: ios, unit_number, count_arg, nb_argument
    logical :: pdb2_bool, complexes_bool, nb_encounters_bool

    ! Initialize flags/defaults
    pdb2_bool = .false.
    complexes_bool = .false.
    nb_encounters_bool = .false.
    dist_threshold = 6.0      ! default threshold in Angstroms
    nb_argument = 0
    count_arg = 1
    nb_argument = command_argument_count()

    ! Parse command-line arguments: expected flags are
    ! -pdb2 <file> : PDB file for solute2
    ! -complexes <file> : association/complexes file
    ! -nb_encounters <int> : optional number of encounters to use
    ! -threshold <real> : distance threshold in Angstroms
    do while ( count_arg <= nb_argument )
        call getarg( count_arg, argument )
        if ( trim(argument) == "-pdb2" ) then
            pdb2_bool = .true.
            call getarg( count_arg+1, argument )
            pdb2_filename = trim(argument)
            count_arg = count_arg + 1
        else if ( trim(argument) == "-complexes" ) then
          complexes_bool = .true.
          call getarg( count_arg+1, argument )
          complexes_filename = trim(argument)
          count_arg = count_arg + 1
        else if ( trim(argument) == "-nb_encounters" ) then
          nb_encounters_bool = .true.
          call getarg( count_arg+1, argument )
          read(argument, *, IOSTAT=ios) nb_encounters
          if (ios /= 0) then
            print *, "ERROR. Integer expected for the -nb_encounters ", &
            "argument."
          end if
          count_arg = count_arg + 1
        else if ( trim(argument) == "-threshold" ) then
            call getarg( count_arg+1, argument )
            read(argument, *, IOSTAT=ios) dist_threshold
            if (ios /= 0) then
              print *, "ERROR. Integer expected for the -threshold ", &
              "argument."
            end if
            count_arg = count_arg + 1
        end if

        count_arg = count_arg + 1
    end do


    ! Read input files into program structures
    ! - read_pdb allocates and fills `pdb2` and returns total atoms, residues and chains
    call read_pdb(pdb2, pdb2_filename, tot_atoms2, tot_residues2, tot_chains2)

    ! - read_assoc allocates and fills `complexes` and returns total encounters
    call read_assoc(complexes, complexes_filename, tot_encounters)

    ! Validate or set `nb_encounters` based on the complexes file
    if ( nb_encounters .gt. tot_encounters ) then
        write(*,*) ''
        write(*,*) 'WARNING: the number of encounters given is greater than'
        write(*,*) 'the number of encounters in ', complexes_filename
        write(*,*) 'Setting the number of encounters to the total'
        write(*,*) 'number of encounters available'
        nb_encounters = tot_encounters
        write(*,*) 'Total encounters available: ', nb_encounters
        write(*,*) ''
      else if (.not. nb_encounters_bool) then
        write(*,*) 'WARNING: number of encounters not given.'
        write(*,*) ''
        write(*,*) 'Setting the number of encounters to the total'
        write(*,*) 'number of encounters available.'
        nb_encounters = tot_encounters
        write(*,*) ''
        write(*,*) 'If you want a specific number of encounters,'
        write(*,*) 'please give it with the -nb_encounters option.'
        write(*,*) ''
      end if

    ! Compute center-of-geometry (COG) for every residue in pdb2
    ! `residues_cog(j,:)` will hold the COG of residue j
    allocate(residues_cog(tot_residues2, 3))

    do j = 1, tot_residues2
        tot_residue_atoms2 = pdb2 % residues(j) % natoms
        allocate(residue_crds(tot_residue_atoms2, 3))
        do k = 1, tot_residue_atoms2
          ! Copy atom coordinates for this residue into a temporary array
          residue_crds(k, :) = pdb2 % residues(j) % atoms(k) % coord
        end do
        ! calculate_cog computes centroid for the given residue coordinates
        call calculate_cog(residue_cog, residue_crds, tot_residue_atoms2)
        residues_cog(j, :) = residue_cog
        deallocate(residue_crds)
    end do
    
    ! Allocate arrays used to collect results per-residue across encounters
    allocate(new_residues_cog(tot_residues2, 3))
    allocate(encounters_count(tot_residues2))
    allocate(all_closest_resid(tot_residues2, tot_encounters))

    encounters_count = 0
    all_closest_resid = 0

    ! Transform residue centroids for each encounter and test threshold
    ! For each encounter i the transformations defined in `complexes`
    ! to the residue COGs is applied and transformed positions are stored in
    ! `new_residues_cog`. If the transformed Z coordinate is below the
    ! `dist_threshold`, the encounter index is recorded for that residue.
    do i = 1, tot_encounters

      call update_complex(complexes % xc1, complexes % xc2, &
        complexes % trans_vector (i, :), complexes % rot1 (i, :), complexes % rot2 (i, :), &
        tot_residues2, residues_cog, new_residues_cog)

      do j = 1, tot_residues2
        ! Check transformed Z coordinate against threshold
        if ( new_residues_cog(j, 3) .lt. dist_threshold ) then
          encounters_count(j) = encounters_count(j) + 1
          all_closest_resid(j, encounters_count(j)) = i
        end if

      end do

    end do


    ! Write a per-residue list of encounters that bring it within the threshold
    unit_number = 20
    open(unit=unit_number, file='closest_residues.txt', status='replace', action='write')
    do j = 1, tot_residues2

      closest_resid = pdb2 % residues(j) % resid
      write(str_resid, '(I4)') closest_resid
      residue_str = "Residue" // adjustr(str_resid) // ":"
      if (encounters_count(j) .gt. 0) then
        ! Write the list of encounter indexes for this residue on one line
        write(unit_number, "(1(A12)10000(I5,1X))") residue_str, (all_closest_resid(j, k), k = 1, encounters_count(j))
      end if

    end do

    contains
    !> Print usage/help text for `analyze_residues`.
    !!
    !! Shows required command-line options and examples for running the
    !! residue analysis tool.
    subroutine print_help()

      !STOP 1

      print *, ""
      print *, "This program receives as inputs the solute 2 pdb file, ", &
      "the complexes filename, the number of encounters to analyze and ", &
      "the distance threshold (A) to consider residues close to the surface"
      print *, ""
      print *, "Eg.: ./analyze_residues -pdb2 p2_noh.pdb -complexes assoc_complexes, ", &
      "-nb_encounters 5000 -threshold 6.0"
      print *, ""
      print *, "Default values: "
      print *, "* nb_encounters: Encounter complexes to process from the complexes file."
      print *, "* threshold: 6.0"
      print *, ""
      STOP
      
    end subroutine print_help
  

end program main