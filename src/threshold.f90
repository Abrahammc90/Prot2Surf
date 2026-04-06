!> \file threshold.f90
!! \brief Threshold-based encounter selection utility.
!!
!! threshold.f90 - Select encounters based on threshold criteria.
!!
!! Usage: ./threshold -help
!!
!! Computes thresholds (z_coord, atoms_dist, angles) and writes selected encounters to output. See print_help for usage.
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0
!! @date 2026-04-05
program main
  USE read_input
  USE mod_threshold


  implicit none

  ! High-level workflow:
  ! 1) Parse CLI options and validate threshold-related inputs.
  ! 2) Read PDB(s) and complexes structures.
  ! 3) Compute per-encounter scalar metric for selected threshold_type.
  ! 4) Filter encounters by cutoff and write reduced outputs.
      

    character*128 :: pdb1_filename, pdb2_filename, complexes_filename, output_complexes_filename, datadist_filename
    character*128 :: threshold_type
    character*128, dimension(:), allocatable :: arr_atoms1a, arr_atoms1b, arr_atoms2a, arr_atoms2b

    integer :: tot_atoms2, tot_encounters
    integer :: tot_chains2, tot_residues2
    integer :: tot_atoms1, tot_residues1, tot_chains1
    integer :: nb_encounters
    type ( type_pdb_file ) :: pdb1, pdb2
    type ( type_assoc_file ) :: complexes

    integer :: i
    real (kind=8) :: cutoff
    real (kind=8), dimension(:), allocatable :: distarray
    integer, dimension(:), allocatable :: encounter_indexes
    real (kind=8), dimension(:, :), allocatable :: atoms1a_coords, atoms1b_coords
    real (kind=8), dimension(:, :), allocatable :: atoms2a_coords, atoms2b_coords
    real (kind=8), dimension(3) :: cog1a, cog1b, cog2a, cog2b
    logical :: pdb1_bool, pdb2_bool, complexes_bool, array_bool
    logical :: arr_atoms1a_bool, arr_atoms1b_bool, arr_atoms2a_bool, arr_atoms2b_bool
    logical :: cutoff_bool, thresholdtype_bool, output_complexes_bool
    logical :: nb_encounters_bool, help_bool

    call parse_arguments( &
      pdb1_filename, pdb1_bool, &
      pdb2_filename, pdb2_bool, &
      complexes_filename, complexes_bool, &
      datadist_filename, array_bool, &
      output_complexes_filename, output_complexes_bool, &
      nb_encounters, nb_encounters_bool, &
      cutoff, cutoff_bool, &
      threshold_type, thresholdtype_bool, &
      arr_atoms1a, arr_atoms1a_bool, &
      arr_atoms1b, arr_atoms1b_bool, &
      arr_atoms2a, arr_atoms2a_bool, &
      arr_atoms2b, arr_atoms2b_bool, &
      help_bool)


    ! If help is requested, stop before reading files and filtering.
    if ( help_bool ) then
      if (thresholdtype_bool) then
        call print_help(threshold_type)
      else
        call print_help("main")
      end if
      STOP 0
    end if

    ! Check flags required in all threshold modes.
    ! Extra checks for each mode happen in each branch below.
    if (.not. pdb2_bool ) then
      print *, '[threshold] ERROR: PDB not provided. Use -pdb2 <file>'
      print *, '[threshold] For usage, run: ./threshold -help'
      STOP 1
    else if ( .not. arr_atoms2a_bool ) then 
      print *, '[threshold] ERROR: PDB2 provided but atoms2 not provided. Use -atoms2 <group>'
      print *, '[threshold] For usage, run: ./threshold -help'
      STOP 1
    else if (.not. array_bool) then
      print *, '[threshold] ERROR: Array output filename not given. Use -array <file>'
      print *, '[threshold] For usage, run: ./threshold -help'
      STOP 1
    else if ( .not. complexes_bool ) then
      print *, "ERROR. Encounter complexes file not provided."
      print *, "Please provide the complexes file with the '-complexes' option"
      print *, "For more information, please use the -help option:"
      print *, "./threshold -help"
      STOP 1
    else if ( .not. output_complexes_bool ) then
      print *, "ERROR. Output filename not provided."
      print *, "Please provide the complexes file with the '-complexes_output' option"
      print *, "For more information, please use the -help option:"
      print *, "./threshold -help"
      STOP 1
    else if ( .not. cutoff_bool ) then
      print *, "ERROR. Cut off not provided."
      print *, "Please provide the complexes file with the '-cutoff' option"
      print *, "For more information, please use the -help option:"
      print *, "./threshold -help"
      STOP 1
    else if ( .not. thresholdtype_bool ) then
      print *, "ERROR. threshold type not given."
      print *, "Please, provide the threshold type with the '-threshold_type' option"
      print *, "For more information, please use the -help option:"
      print *, "./threshold -help"
      STOP 1
    end if

    
    ! Read structural coordinates for solute 2 (always required).
    call read_pdb(pdb2, pdb2_filename, tot_atoms2, tot_residues2, tot_chains2)

    ! Read encounter transforms and metadata.
    call read_assoc(complexes, complexes_filename, tot_encounters)


    ! Keep number of encounters inside available range.
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

    allocate(distarray (nb_encounters))
    allocate(encounter_indexes (nb_encounters))

    if (trim(adjustl(threshold_type)) == "z_coord") then

      ! z_coord threshold mode:
      ! build one z-based scalar value per encounter and filter by cutoff.

      if (.not. output_complexes_bool) then
        print *, "ERROR. Complexes output filename not given."
        print *, "Please provide the complexes output file with the '-complexes_output option'"
        print *, "For more information, please use the -help option:"
        print *, "./threshold -threshold_type 3D_angle -help"
        STOP 1
      end if

      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call array_z_coord(distarray, nb_encounters, 1, complexes % xc1, complexes % xc2, &
      complexes % trans_vector, complexes % rot1, complexes % rot2, atoms2a_coords)

    else if (trim(adjustl(threshold_type)) == "3D_angle" .or. trim(adjustl(threshold_type)) == "2D_angle") then

      ! angle threshold modes:
      ! compute one 2D/3D angle value for each encounter.

      if (.not. pdb1_bool) then
        print *, "ERROR. Just one PDB provided. It is necessary to provide two pdbs ", &
        "to generate the angle threshold."
        print *, "Please provide the two pdb files with '-pdb1' and '-pdb2' options."
        print *, "For more information, please use the -help option:"
        print *, "./threshold -threshold_type 3D_angle -help"
        STOP 1
      else if (.not. arr_atoms1a_bool .or. .not. arr_atoms1b_bool .or. .not. arr_atoms2b_bool) then 
        print *, "ERROR. Not all necessary group of atoms given"
        print *, "Please provide all group of atoms with '-atoms1a', '-atoms1b', '-atoms2a' and '-atoms2b' options"
        print *, "For more information, please use the -help option:"
        print *, "./threshold -threshold_type 3D_angle -help"
        STOP 1
      else if (.not. output_complexes_bool) then
        print *, "ERROR. Complexes output filename not given."
        print *, "Please provide the complexes output file with the '-complexes_output option'"
        print *, "For more information, please use the -help option:"
        print *, "./threshold -threshold_type 3D_angle -help"
        STOP 1
      end if

      call read_pdb(pdb1, pdb1_filename, tot_atoms1, tot_residues1, tot_chains1)

      call read_atoms_coord(arr_atoms1a, atoms1a_coords, tot_atoms1, pdb1, pdb1_filename)
      call read_atoms_coord(arr_atoms1b, atoms1b_coords, tot_atoms1, pdb1, pdb1_filename)
      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)
      call read_atoms_coord(arr_atoms2b, atoms2b_coords, tot_atoms2, pdb2, pdb2_filename)

      call calculate_cog(cog1a, atoms1a_coords, size(atoms1a_coords(:, 3)))
      call calculate_cog(cog1b, atoms1b_coords, size(atoms1b_coords(:, 3)))
      call calculate_cog(cog2a, atoms2a_coords, size(atoms2a_coords(:, 3)))
      call calculate_cog(cog2b, atoms2b_coords, size(atoms2b_coords(:, 3)))
      
      if (trim(adjustl(threshold_type)) == "2D_angle") then
        ! 2D projected angle mode.
        call array_angle(distarray, nb_encounters, 2, &
        complexes % xc1, complexes % xc2, complexes % trans_vector, complexes % rot1, complexes % rot2, &
        cog1a, cog1b, cog2a, cog2b, 2)
      else if (trim(adjustl(threshold_type)) == "3D_angle") then
        ! Full 3D angle mode.
        call array_angle(distarray, nb_encounters, 2, &
        complexes % xc1, complexes % xc2, complexes % trans_vector, complexes % rot1, complexes % rot2, &
        cog1a, cog1b, cog2a, cog2b, 3)
      end if

    else if (trim(adjustl(threshold_type)) == "atoms_dist") then

      ! atoms_dist threshold mode:
      ! compute minimum inter-group distance per encounter and filter by cutoff.

      if (.not. pdb1_bool) then
        print *, "ERROR. Just one PDB provided. It is necessary to provide two pdbs ", &
        "to generate the angle threshold."
        print *, "Please provide the two pdb files with '-pdb1' and '-pdb2' options"
        print *, "For more information, please use the -help option:"
        print *, "./threshold -threshold atoms_dist -help"
        STOP 1
      else if (.not. arr_atoms1a_bool ) then 
        print *, "ERROR. Not all necessary group of atoms given"
        print *, "Please provide the group of atoms with '-atoms1' and '-atoms2' options"
        print *, "For more information, please use the -help option:"
        print *, "./threshold -threshold atoms_dist -help"
        STOP 1
      else if (.not. output_complexes_bool) then
        print *, "ERROR. Complexes output filename not given."
        print *, "Please provide the complexes output file with the '-complexes_output option'"
        print *, "For more information, please use the -help option:"
        print *, "./threshold -threshold_type atoms_dist -help"
        STOP 1
      end if

      call read_pdb(pdb1, pdb1_filename, tot_atoms1, tot_residues1, tot_chains1)
      call read_atoms_coord(arr_atoms1a, atoms1a_coords, tot_atoms1, pdb1, pdb1_filename)
      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call array_atoms_dist(distarray, nb_encounters, 1, complexes % xc1, complexes % xc2, &
      complexes % trans_vector, complexes % rot1, complexes % rot2, atoms1a_coords, atoms2a_coords)

    else
      print *, 'ERROR: ', trim(adjustl(threshold_type)), ' threshold not implemented'
      print *, "To know more about the type of matrices available, please run:"
      print *, "./threshold -help"
      STOP
    end if

    ! Build identity index array used to keep track of encounter ordering.
    encounter_indexes = 0
    do i = 1, nb_encounters
      encounter_indexes(i) = i
    end do

    ! Sort metric values and reorder encounter indexes in sync.
    call sort_array(distarray, encounter_indexes)

    call write_array(distarray, datadist_filename, cutoff)
    call write_cutoff_complexes(complexes, distarray, encounter_indexes, nb_encounters, cutoff, output_complexes_filename)

    print *, 'Threshold complete and encounters stored in ', output_complexes_filename
  

  contains
      !> Parse command-line arguments for `threshold`.
      !!
      !! Parses all CLI flags and fills program configuration values.
      !! Atom-group flags consume multiple tokens until the next flag.
      !! @param[out] p_pdb1_filename          Value parsed from `-pdb1`
      !! @param[out] p_pdb1_bool              Whether `-pdb1` was provided
      !! @param[out] p_pdb2_filename          Value parsed from `-pdb2`
      !! @param[out] p_pdb2_bool              Whether `-pdb2` was provided
      !! @param[out] p_complexes_filename     Value parsed from `-complexes`
      !! @param[out] p_complexes_bool         Whether `-complexes` was provided
      !! @param[out] p_datadist_filename      Value parsed from `-array`
      !! @param[out] p_array_bool             Whether `-array` was provided
      !! @param[out] p_output_complexes_filename Value parsed from `-complexes_output`
      !! @param[out] p_output_complexes_bool  Whether `-complexes_output` was provided
      !! @param[out] p_nb_encounters          Value parsed from `-nb_encounters`
      !! @param[out] p_nb_encounters_bool     Whether `-nb_encounters` was provided
      !! @param[out] p_cutoff                 Value parsed from `-cutoff`
      !! @param[out] p_cutoff_bool            Whether `-cutoff` was provided
      !! @param[out] p_threshold_type         Value parsed from `-threshold_type`
      !! @param[out] p_thresholdtype_bool     Whether `-threshold_type` was provided
      !! @param[out] p_arr_atoms1a            Atom list parsed from `-atoms1` or `-atoms1a`
      !! @param[out] p_arr_atoms1a_bool       Whether `-atoms1` or `-atoms1a` was provided
      !! @param[out] p_arr_atoms1b            Atom list parsed from `-atoms1b`
      !! @param[out] p_arr_atoms1b_bool       Whether `-atoms1b` was provided
      !! @param[out] p_arr_atoms2a            Atom list parsed from `-atoms2` or `-atoms2a`
      !! @param[out] p_arr_atoms2a_bool       Whether `-atoms2` or `-atoms2a` was provided
      !! @param[out] p_arr_atoms2b            Atom list parsed from `-atoms2b`
      !! @param[out] p_arr_atoms2b_bool       Whether `-atoms2b` was provided
      !! @param[out] p_help_bool              Whether `-help` was provided
      subroutine parse_arguments( &
        p_pdb1_filename, p_pdb1_bool, &
        p_pdb2_filename, p_pdb2_bool, &
        p_complexes_filename, p_complexes_bool, &
        p_datadist_filename, p_array_bool, &
        p_output_complexes_filename, p_output_complexes_bool, &
        p_nb_encounters, p_nb_encounters_bool, &
        p_cutoff, p_cutoff_bool, &
        p_threshold_type, p_thresholdtype_bool, &
        p_arr_atoms1a, p_arr_atoms1a_bool, &
        p_arr_atoms1b, p_arr_atoms1b_bool, &
        p_arr_atoms2a, p_arr_atoms2a_bool, &
        p_arr_atoms2b, p_arr_atoms2b_bool, &
        p_help_bool)
        implicit none

        character*128, intent(out) :: p_pdb1_filename
        logical, intent(out) :: p_pdb1_bool
        character*128, intent(out) :: p_pdb2_filename
        logical, intent(out) :: p_pdb2_bool
        character*128, intent(out) :: p_complexes_filename
        logical, intent(out) :: p_complexes_bool
        character*128, intent(out) :: p_datadist_filename
        logical, intent(out) :: p_array_bool
        character*128, intent(out) :: p_output_complexes_filename
        logical, intent(out) :: p_output_complexes_bool
        integer, intent(out) :: p_nb_encounters
        logical, intent(out) :: p_nb_encounters_bool
        real (kind=8), intent(out) :: p_cutoff
        logical, intent(out) :: p_cutoff_bool
        character*128, intent(out) :: p_threshold_type
        logical, intent(out) :: p_thresholdtype_bool
        character*128, dimension(:), allocatable, intent(out) :: p_arr_atoms1a
        logical, intent(out) :: p_arr_atoms1a_bool
        character*128, dimension(:), allocatable, intent(out) :: p_arr_atoms1b
        logical, intent(out) :: p_arr_atoms1b_bool
        character*128, dimension(:), allocatable, intent(out) :: p_arr_atoms2a
        logical, intent(out) :: p_arr_atoms2a_bool
        character*128, dimension(:), allocatable, intent(out) :: p_arr_atoms2b
        logical, intent(out) :: p_arr_atoms2b_bool
        logical, intent(out) :: p_help_bool

        integer :: nb_argument, count_arg, n_atoms, j, ios
        character*128 :: argument

        nb_argument = command_argument_count()
        if ( nb_argument == 0 ) then
          print *, '[threshold] ERROR: No options given.'
          print *, '[threshold] For usage, run: ./threshold -help'
          call print_help("main")
          STOP 1
        end if

        p_nb_encounters = 0
        p_cutoff = 0.d0
        p_pdb1_bool = .false.
        p_pdb2_bool = .false.
        p_complexes_bool = .false.
        p_cutoff_bool = .false.
        p_array_bool = .false.
        p_thresholdtype_bool = .false.
        p_output_complexes_bool = .false.
        p_nb_encounters_bool = .false.
        p_arr_atoms1a_bool = .false.
        p_arr_atoms1b_bool = .false.
        p_arr_atoms2a_bool = .false.
        p_arr_atoms2b_bool = .false.
        p_help_bool = .false.
        argument = ""
        count_arg = 1

        do while ( count_arg <= nb_argument )
          call getarg( count_arg, argument )

          if ( trim(argument) == "-pdb1" ) then
            p_pdb1_bool = .true.
            call getarg( count_arg+1, argument )
            p_pdb1_filename = trim(argument)
            count_arg = count_arg + 1
          else if ( trim(argument) == "-pdb2" ) then
            p_pdb2_bool = .true.
            call getarg( count_arg+1, argument )
            p_pdb2_filename = trim(argument)
            count_arg = count_arg + 1
          else if ( trim(argument) == "-complexes" ) then
            p_complexes_bool = .true.
            call getarg( count_arg+1, argument )
            p_complexes_filename = trim(argument)
            count_arg = count_arg + 1
          else if ( trim(argument) == "-array" ) then
            p_array_bool = .true.
            call getarg( count_arg+1, argument )
            p_datadist_filename = trim(argument)
            count_arg = count_arg + 1
          else if ( trim(argument) == "-complexes_output" ) then
            p_output_complexes_bool = .true.
            call getarg( count_arg+1, argument )
            p_output_complexes_filename = trim(argument)
            count_arg = count_arg + 1
          else if ( trim(argument) == "-nb_encounters" ) then
            p_nb_encounters_bool = .true.
            call getarg( count_arg+1, argument )
            read(argument, *, IOSTAT=ios) p_nb_encounters
            if (ios /= 0) then
              print *, "ERROR. Integer expected for the -nb_encounters argument."
            end if
            count_arg = count_arg + 1
          else if ( trim(argument) == "-cutoff" ) then
            p_cutoff_bool = .true.
            call getarg( count_arg+1, argument )
            read(argument, *, IOSTAT=ios) p_cutoff
            if (ios /= 0) then
              print *, "ERROR. Float expected for the -cutoff argument."
            end if
            count_arg = count_arg + 1
          else if ( trim(argument) == "-threshold_type" ) then
            p_thresholdtype_bool = .true.
            call getarg( count_arg+1, argument )
            p_threshold_type = trim(argument)
            count_arg = count_arg + 1
          else if ( trim(argument) == "-atoms1" .or. trim(argument) == "-atoms1a" ) then
            p_arr_atoms1a_bool = .true.

            call getarg( count_arg+1, argument )
            n_atoms = 0
            do while (argument(1:1) /= "-" .and. count_arg + n_atoms <= nb_argument)
              n_atoms = n_atoms + 1
              call getarg( count_arg+n_atoms, argument )
            end do
            n_atoms = n_atoms - 1
            allocate(p_arr_atoms1a(n_atoms))

            do j = 1, n_atoms
              call getarg( count_arg+j, argument )
              p_arr_atoms1a(j) = trim(argument)
            end do

            count_arg = count_arg + n_atoms
          else if ( trim(argument) == "-atoms1b" ) then
            p_arr_atoms1b_bool = .true.

            call getarg( count_arg+1, argument )
            n_atoms = 0
            do while (argument(1:1) /= "-" .and. count_arg + n_atoms <= nb_argument)
              n_atoms = n_atoms + 1
              call getarg( count_arg+n_atoms, argument )
            end do
            n_atoms = n_atoms - 1
            allocate(p_arr_atoms1b(n_atoms))

            do j = 1, n_atoms
              call getarg( count_arg+j, argument )
              p_arr_atoms1b(j) = trim(argument)
            end do

            count_arg = count_arg + n_atoms
          else if ( trim(argument) == "-atoms2" .or. trim(argument) == "-atoms2a" ) then
            p_arr_atoms2a_bool = .true.

            call getarg( count_arg+1, argument )
            n_atoms = 0
            do while (argument(1:1) /= "-" .and. count_arg + n_atoms <= nb_argument)
              n_atoms = n_atoms + 1
              call getarg( count_arg+n_atoms, argument )
            end do
            n_atoms = n_atoms - 1
            allocate(p_arr_atoms2a(n_atoms))

            do j = 1, n_atoms
              call getarg( count_arg+j, argument )
              p_arr_atoms2a(j) = trim(argument)
            end do

            count_arg = count_arg + n_atoms
          else if ( trim(argument) == "-atoms2b" ) then
            p_arr_atoms2b_bool = .true.

            call getarg( count_arg+1, argument )
            n_atoms = 0
            do while (argument(1:1) /= "-" .and. count_arg + n_atoms <= nb_argument)
              n_atoms = n_atoms + 1
              call getarg( count_arg+n_atoms, argument )
            end do
            n_atoms = n_atoms - 1
            allocate(p_arr_atoms2b(n_atoms))

            do j = 1, n_atoms
              call getarg( count_arg+j, argument )
              p_arr_atoms2b(j) = trim(argument)
            end do

            count_arg = count_arg + n_atoms
          else if ( trim(argument) == "-help" ) then
            p_help_bool = .true.
          else
            print *, '[threshold] ERROR: Argument ''', trim(argument), ''' not recognized.'
            print *, '[threshold] For usage, run: ./threshold -help'
            STOP 1
          end if

          count_arg = count_arg + 1
        end do

      end subroutine parse_arguments

    !> Print usage/help text for `threshold`.
    !!
    !! Explains command-line options and provides examples per
    !! threshold type.
    subroutine print_help(help_option)

      implicit none
      character(len=*), intent(in) :: help_option

      if (trim(help_option) == "main") then
        print *, ""
        print *, "This program receives as inputs the pdb files, atom groups (either indexes or atomnames)", &
        "encounter complexes file, cutoff to make the threshold distance, ", &
        "threshold type, array output file and complexes output file."
        print *, "If you want instructions on how to create a specific type of threshold, type:"
        print *, "./threshold -threshold_type <type of threshold> -help"
        print *, ""
        print *, "Eg.: ./threshold -threshold_type z_coord -help"
        print *, ""
        print *, "The threshold types currently available are:"
        print *, "* z_coord"
        print *, "* atoms_dist"
        print *, "* 2D_angle"
        print *, "* 3D_angle"
        print *, ""
        print *, "General required options: -pdb2 -atoms2 -complexes -array -complexes_output -cutoff -threshold_type"
        print *, ""
        STOP
      else if (trim(help_option) == "z_coord") then
        print *, ""
        print *, "To generate the z_coord threshold you will need to give the pdb of solute 2,", &
        "one group of atoms (either indexes or atomnames), the encounter complexes file, ", &
        "number of encounters, cutoff, threshold type, array output file and threshold output file."
        print *, ""
        print *, "Eg.: ./threshold -pdb2 p2_noh.pdb -atoms2 Cu -complexes assoc_complexes ", &
        "-threshold_type z_coord -cutoff 6.0 -array threshold_z_array.txt -complexes_output threshold_z.txt"
        print *, ""
        STOP
      else if (trim(help_option) == "atoms_dist") then
        print *, ""
        print *, "To generate the atoms_dist threshold you will need to give two pdbs, two groups of atoms ", &
        "(either indexes or atomnames), the encounter complexes file, number of encounters, cutoff, threshold type, ", &
        "array output file and threshold output file."
        print *, ""
        print *, "Eg.: ./threshold -pdb1 p1_noh.pdb -pdb2 p2_noh.pdb -complexes ", &
        "assoc_complexes -atoms1 C1 C4 -atoms2 Cu -threshold_type atoms_dist -cutoff 6.0 ", &
        "-array threshold_dist_array.txt -complexes_output threshold_dist.txt"
        print *, ""
        STOP
      else if (trim(help_option) == "2D_angle" .or. trim(help_option) == "3D_angle") then
        print *, ""
        print *, "To generate the angle threshold you will need to give two pdbs, four groups of atoms ", &
        "(either indexes or atomnames), the encounter complexes file, number of encounters, cutoff, threshold ", & 
        "type, array output file and threshold output file."
        print *, ""
        print *, "Eg.: ./threshold -pdb1 p1_noh.pdb -pdb2 p2_noh.pdb -complexes ", &
        "assoc_complexes -atoms1a 1 -atoms1b 411 -atoms2a Cu -atoms2b 4 5 6 7 8 -threshold_type 3D_angle ", &
        "-cutoff 30.0 -array threshold_angle_array.txt -complexes_output threshold_angle.txt"
        print *, ""
        print *, "For this specific threshold type you need to be aware that a center ", &
        "of geometry will be calculated for each group of atoms and indexed accordingly: "
        print *, "atoms1a -> cog1a, atoms1b -> cog1b, atoms2a -> cog2a, atoms2b -> cog2b"
        print *, "Then, the 1st vector will be calculated as (cog1b - cog1a) and the ", &
        "2nd vector will be calculated as (cog2b - cog2a)"
        print *, "Last, the angle between the two vectors will be computed."
        print *, ""
        STOP
      end if
      
    end subroutine print_help

    !> Read coordinates for a specified list of atoms (names or indices).
    !!
    !! Accepts either atom indices or atom names. If atom names are
    !! provided, all matching atoms in the PDB are returned (in order).
    !! @param[in]  atoms_array   Array of atom name strings or index strings
    !! @param[out] atoms_coord   Allocated array with atom coordinates (n,3)
    !! @param[in]  tot_atoms     Total number of atoms in `pdb`
    !! @param[in]  pdb           PDB structure holding atom coordinates
    !! @param[in]  pdb_filename  Reference filename for warnings
    subroutine read_atoms_coord(atoms_array, atoms_coord, tot_atoms, pdb, pdb_filename)
      implicit none

      character(len=*), dimension(:), intent(in) :: atoms_array
      real(kind=8), dimension(:, :), allocatable, intent(out) :: atoms_coord
      integer, intent(in) :: tot_atoms
      type ( type_pdb_file ), intent(in) :: pdb
      character(len=*), intent(in) :: pdb_filename

      character (len=4), dimension(:), allocatable :: array_atomnames
      character (len=4) :: atomname
      integer :: n, k, l, j, ios_local
      integer :: atom_index, atoms_count
      logical :: read_atomname, atom_found

      read_atomname = .false.

      n = size(atoms_array)


      read(atoms_array(1), *, IOSTAT=ios_local) atom_index

      if (ios_local /= 0) then
        read_atomname = .true.
      end if

      
      if (read_atomname) then

        allocate(array_atomnames(n))

        do j = 1, n
          read(atoms_array(j), *, IOSTAT=ios_local) atomname
          array_atomnames(j) = atomname
        end do
        
        atoms_count = 0
        do j = 1, n
          atom_found = .false.
          !print *, array_atomnames(j)
          do k = 1, tot_atoms
            if ( trim(adjustl(pdb % atoms(k) % name)) == &
            trim(adjustl(array_atomnames(j))) ) then
              atoms_count = atoms_count + 1
              atom_found = .true.
            end if
          end do
          if (.not. atom_found) then
            print *, "WARNING: atom ", trim(adjustl(array_atomnames(j))), &
            " not found in ", pdb_filename
          end if
          
        end do
        
        
        allocate(atoms_coord(atoms_count, 3))

        atoms_count = 0
        do j = 1, n
          do k = 1, tot_atoms
            if ( trim(adjustl(pdb % atoms(k) % name)) == &
            trim(adjustl(array_atomnames(j))) ) then
              atoms_count = atoms_count + 1
              atoms_coord(atoms_count, :) = pdb % atoms(k) % coord(:)
            end if
          end do
        end do
      
      else
      
        allocate(atoms_coord(n, 3))
        do j = 1, n
          read(atoms_array(j), *, IOSTAT=ios_local) l
          atoms_coord(j, :) = pdb % atoms(l) % coord(:)
        end do

      end if

    end subroutine read_atoms_coord

    !> Write selected encounter complexes to an output file.
    !!
    !! Writes the header lines from `complexes%head` and then appends
    !! encounter lines (from `complexes%lines`) in the order provided
    !! by `encounter_indexes` until a value exceeding `cutoff` is
    !! encountered.
    subroutine write_cutoff_complexes(complexes_arg, distarray_arg, encounter_indexes_arg, nb_encounters_arg, &
                                      cutoff_arg, output_name_arg)

      IMPLICIT NONE
      type(type_assoc_file), intent(in) :: complexes_arg
      real (kind=8), dimension(:), intent(in) :: distarray_arg
      integer, dimension(:), intent(in) :: encounter_indexes_arg
      integer, intent(in) :: nb_encounters_arg
      real(kind=8), intent(in) :: cutoff_arg
      character*128, intent(in) :: output_name_arg

      integer :: idx, unit_number = 15


      ! Write the head of the complexes file
      open(unit=unit_number, file=output_name_arg, status='replace', action='write')

      do idx = 1, 4
        write(unit_number, "(A)") complexes_arg % head(idx) 
      end do

      ! Write the encounters that are below the cutoff
      do idx = 1, nb_encounters_arg       
        write(unit_number, "(A)") complexes_arg % lines(encounter_indexes_arg(idx))
        if (distarray_arg(idx) > cutoff_arg) exit 
      end do

      ! Close the file
      close(unit_number)

    end subroutine write_cutoff_complexes

    !> Write array values up to a cutoff to a formatted file.
    !!
    !! Writes values from `array` (in ascending order) until `cutoff` is
    !! reached and stores them in `filename` using `F10.4` formatting.
    subroutine write_array(array_arg, filename_arg, cutoff_arg)

      IMPLICIT NONE
      real(kind=8), dimension(:), intent(in) :: array_arg
      character*128, intent(in) :: filename_arg
      real(kind=8), intent(in) :: cutoff_arg
      character*128 :: fmt
      
      integer :: idx, count
      integer :: unit = 15
      integer :: nelements
      
      nelements = size(array_arg(:))
      count = 0

      open(unit, file=trim(filename_arg), status='replace', action='write')

      write(fmt, '(I0, A)') nelements, 'F10.4'
      fmt = '(' // trim(fmt) // ')'

      do idx = 1, nelements
        if (array_arg(idx) > cutoff_arg) exit
        count = count + 1
      end do

      ! Write array
      write(unit, fmt) (array_arg(idx), idx = 1, count)
      
      close(unit)
      
    end subroutine write_array

END PROGRAM main