!> \file make_data.f90
!! \brief Build encounter-based matrix/array data from structural inputs.
!!
!! make_data.f90 - Build encounter-based matrices from PDB and complex files.
!!
!! Usage: ./make_data -help
!!
!! Reads PDB(s) and association/complex files, generates matrices (rmsd, z_coord, atoms_dist, angle) and arrays.
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0
!! @date 2026-04-05
program main
  USE read_input
  USE mod_matrix
  USE mod_array


  implicit none

  ! High-level workflow:
  ! 1) Parse CLI options and validate required inputs for selected data_type.
  ! 2) Read PDB/complexes data into in-memory derived types.
  ! 3) Build the requested metric (z_coord, rmsd, atoms_dist, 2D/3D angle).
  ! 4) Write resulting array/matrix and (optionally sorted) complexes.
      

    character*128 :: pdb1_filename, pdb2_filename, complexes_filename, sorted_complexes_filename
    character*128 :: data_type, datadist_filename
    character*128, dimension(:), allocatable :: arr_atoms1a, arr_atoms1b, arr_atoms2a, arr_atoms2b

    integer :: tot_atoms2, tot_encounters, tot_coords
    integer :: tot_chains2, tot_residues2
    integer :: tot_atoms1, tot_residues1, tot_chains1
    integer :: nb_encounters
    type ( type_pdb_file ) :: pdb1, pdb2
    type ( type_assoc_file ) :: complexes

    real (kind=8), dimension(:, :), allocatable :: distmatrix
    real (kind=8), dimension(:), allocatable :: distarray
    real (kind=8), dimension(:, :), allocatable :: atoms1a_coords, atoms1b_coords
    real (kind=8), dimension(:, :), allocatable :: atoms2a_coords, atoms2b_coords
    real (kind=8), dimension(3) :: cog1a, cog1b, cog2a, cog2b
    !real (kind=8), dimension(:, :), allocatable :: residues_cog, residue_crds
    logical :: pdb1_bool, pdb2_bool, complexes_bool
    logical :: matrix_bool, datatype_bool, nb_encounters_bool, help_bool
    logical :: arr_atoms1a_bool, arr_atoms1b_bool, arr_atoms2a_bool, arr_atoms2b_bool
    call parse_arguments( &
      pdb1_filename, pdb1_bool, &
      pdb2_filename, pdb2_bool, &
      complexes_filename, complexes_bool, &
      nb_encounters, nb_encounters_bool, &
      datadist_filename, matrix_bool, &
      data_type, datatype_bool, &
      arr_atoms1a, arr_atoms1a_bool, &
      arr_atoms1b, arr_atoms1b_bool, &
      arr_atoms2a, arr_atoms2a_bool, &
      arr_atoms2b, arr_atoms2b_bool, &
      help_bool)

    ! If help is requested, stop before running checks or calculations.
    if ( help_bool ) then
      if (datatype_bool) then
        call print_help(data_type)
      else
        call print_help("main")
      end if
      STOP 0
    end if

    ! Check flags required in all modes.
    ! Extra checks for each data type are done in each branch below.
    if (.not. pdb2_bool ) then
      print *, '[make_data] ERROR: PDB not provided. Use -pdb2 <file>'
      print *, '[make_data] For usage, run: ./make_data -help'
      STOP 1
    else if ( .not. arr_atoms2a_bool ) then 
      print *, '[make_data] ERROR: PDB2 provided but atoms2 not provided. Use -atoms2 <group>'
      print *, '[make_data] For usage, run: ./make_data -help'
      STOP 1
    else if ( .not. complexes_bool ) then
      print *, '[make_data] ERROR: Encounter complexes file not provided. Use -complexes <file>'
      print *, '[make_data] For usage, run: ./make_data -help'
      STOP 1
    else if ( .not. matrix_bool ) then
      print *, '[make_data] ERROR: Data output filename not provided. Use -input <file>'
      print *, '[make_data] For usage, run: ./make_data -help'
      STOP 1
    else if ( .not. datatype_bool ) then
      print *, '[make_data] ERROR: Data type not given. Use -data_type <type>'
      print *, '[make_data] For usage, run: ./make_data -help'
      STOP 1
    end if


    ! Read structural coordinates for solute 2 (always required).
    call read_pdb(pdb2, pdb2_filename, tot_atoms2, tot_residues2, tot_chains2)
    !call read_atoms_coord(arr_atoms2, atoms2_coords, tot_atoms2, pdb2, pdb2_filename)

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


    sorted_complexes_filename = "sorted_" // trim(complexes_filename)

    ! Branch by requested metric type and compute corresponding output data.
    print *, 'Building ', trim(data_type), ' matrix'

    if (trim(adjustl(data_type)) == "z_coord") then

      ! z_coord mode:
      ! - transform selected solute2 coordinates per encounter
      ! - compute one z value for each encounter

      allocate(distarray (nb_encounters))
      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call array_z_coord(complexes, distarray, nb_encounters, 1, complexes % xc1, complexes % xc2, &
      complexes % trans_vector, complexes % rot1, complexes % rot2, atoms2a_coords)

      ! z_coord creates an array (one value per encounter).
      call write_array(distarray, datadist_filename)
      call write_complexes(complexes, nb_encounters, complexes_filename)

    else if (trim(adjustl(data_type)) == "rmsd") then

      ! rmsd mode:
      ! - generate full pairwise matrix from transformed structures
      ! - matrix is consumed by clustering in a separate step/tool

      allocate(distmatrix (nb_encounters, nb_encounters))
      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call matrix_rmsd(distmatrix, nb_encounters, tot_coords, &
      complexes % xc1, complexes % xc2, complexes % trans_vector, &
      complexes % rot1, complexes % rot2, atoms2a_coords)

      ! rmsd creates a full pairwise matrix.
      call write_matrix(distmatrix, datadist_filename)

    else if (trim(adjustl(data_type)) == "3D_angle" .or. trim(adjustl(data_type)) == "2D_angle") then

      ! angle modes:
      ! - build vectors from COGs of two atom groups per solute
      ! - compute one 2D or 3D angle for each encounter

      if (.not. pdb1_bool) then
        print *, "ERROR. Just one PDB provided. It is necessary to provide two pdbs ", &
        "to generate the angle matrix."
        print *, "Please provide the two pdb files with '-pdb1' and '-pdb2' options."
        print *, "For more information, please use the -help option:"
        print *, "./make_data -data_type 3D_angle -help"
        STOP 1
      else if (.not. arr_atoms1a_bool .or. .not. arr_atoms1b_bool .or. .not. arr_atoms2b_bool) then 
        print *, "ERROR. Not all necessary group of atoms given"
        print *, "Please provide all group of atoms with '-atoms1a', '-atoms1b', '-atoms2a' and '-atoms2b' options"
        print *, "For more information, please use the -help option:"
        print *, "./make_data -data_type 3D_angle -help"
        STOP 1
      end if

      allocate(distarray (nb_encounters))
      call read_pdb(pdb1, pdb1_filename, tot_atoms1, tot_residues1, tot_chains1)

      call read_atoms_coord(arr_atoms1a, atoms1a_coords, tot_atoms1, pdb1, pdb1_filename)
      call read_atoms_coord(arr_atoms1b, atoms1b_coords, tot_atoms1, pdb1, pdb1_filename)
      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)
      call read_atoms_coord(arr_atoms2b, atoms2b_coords, tot_atoms2, pdb2, pdb2_filename)

      call calculate_cog(cog1a, atoms1a_coords, size(atoms1a_coords(:, 3)))
      call calculate_cog(cog1b, atoms1b_coords, size(atoms1b_coords(:, 3)))
      call calculate_cog(cog2a, atoms2a_coords, size(atoms2a_coords(:, 3)))
      call calculate_cog(cog2b, atoms2b_coords, size(atoms2b_coords(:, 3)))
      
      if (trim(adjustl(data_type)) == "2D_angle") then
        ! 2D projected angle mode (xy plane projection).
        call array_angle(complexes, distarray, nb_encounters, 2, &
        complexes % xc1, complexes % xc2, complexes % trans_vector, complexes % rot1, complexes % rot2, &
        cog1a, cog1b, cog2a, cog2b, 2)
      else if (trim(adjustl(data_type)) == "3D_angle") then
        ! Full 3D angle mode.
        call array_angle(complexes, distarray, nb_encounters, 2, &
        complexes % xc1, complexes % xc2, complexes % trans_vector, complexes % rot1, complexes % rot2, &
        cog1a, cog1b, cog2a, cog2b, 3)
      end if

      ! angle modes create an array (one angle per encounter).
      call write_array(distarray, datadist_filename)
      call write_complexes(complexes, nb_encounters, sorted_complexes_filename)

    else if (trim(adjustl(data_type)) == "atoms_dist") then

      ! atoms_dist mode:
      ! - compute minimum inter-group distance per encounter
      ! - write one value per encounter for later analysis/clustering

      if (.not. pdb1_bool) then
        print *, "ERROR. Just one PDB provided. It is necessary to provide two pdbs ", &
        "to generate the angle matrix."
        print *, "Please provide the two pdb files with '-pdb1' and '-pdb2' options"
        print *, "For more information, please use the -help option:"
        print *, "./make_data -data_type atoms_dist -help"
        STOP 1
      else if (.not. arr_atoms1a_bool ) then 
        print *, "ERROR. Not all necessary group of atoms given"
        print *, "Please provide the group of atoms with '-atoms1' and '-atoms2' options"
        print *, "For more information, please use the -help option:"
        print *, "./make_data -data_type atoms_dist -help"
        STOP 1
      end if

      allocate(distarray (nb_encounters))
      call read_pdb(pdb1, pdb1_filename, tot_atoms1, tot_residues1, tot_chains1)
      call read_atoms_coord(arr_atoms1a, atoms1a_coords, tot_atoms1, pdb1, pdb1_filename)
      !STOP 1
      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)
      !STOP 1
      
      !STOP 1

      call array_atoms_dist(complexes, distarray, nb_encounters, 1, complexes % xc1, complexes % xc2, &
      complexes % trans_vector, complexes % rot1, complexes % rot2, atoms1a_coords, atoms2a_coords)

      ! atoms_dist creates an array (minimum distance per encounter).
      call write_array(distarray, datadist_filename)
      call write_complexes(complexes, nb_encounters, sorted_complexes_filename)

    else
      print *, 'ERROR: ', trim(adjustl(data_type)), ' data not implemented'
      print *, "To know more about the type of matrices available, please run:"
      print *, "./make_data -help"
      STOP
    end if

    print *, 'Data processed and stored in ', datadist_filename
  

  contains
      !> Parse command-line arguments for `make_data`.
      !!
      !! Parses all CLI flags and fills program configuration values.
      !! Atom-group flags consume multiple tokens until the next flag.
      !! @param[out] p_pdb1_filename      Value parsed from `-pdb1`
      !! @param[out] p_pdb1_bool          Whether `-pdb1` was provided
      !! @param[out] p_pdb2_filename      Value parsed from `-pdb2`
      !! @param[out] p_pdb2_bool          Whether `-pdb2` was provided
      !! @param[out] p_complexes_filename Value parsed from `-complexes`
      !! @param[out] p_complexes_bool     Whether `-complexes` was provided
      !! @param[out] p_nb_encounters      Value parsed from `-nb_encounters`
      !! @param[out] p_nb_encounters_bool Whether `-nb_encounters` was provided
      !! @param[out] p_datadist_filename  Value parsed from `-input`
      !! @param[out] p_matrix_bool        Whether `-input` was provided
      !! @param[out] p_data_type          Value parsed from `-data_type`
      !! @param[out] p_datatype_bool      Whether `-data_type` was provided
      !! @param[out] p_arr_atoms1a        Atom list parsed from `-atoms1` or `-atoms1a`
      !! @param[out] p_arr_atoms1a_bool   Whether `-atoms1` or `-atoms1a` was provided
      !! @param[out] p_arr_atoms1b        Atom list parsed from `-atoms1b`
      !! @param[out] p_arr_atoms1b_bool   Whether `-atoms1b` was provided
      !! @param[out] p_arr_atoms2a        Atom list parsed from `-atoms2` or `-atoms2a`
      !! @param[out] p_arr_atoms2a_bool   Whether `-atoms2` or `-atoms2a` was provided
      !! @param[out] p_arr_atoms2b        Atom list parsed from `-atoms2b`
      !! @param[out] p_arr_atoms2b_bool   Whether `-atoms2b` was provided
      !! @param[out] p_help_bool          Whether `-help` was provided
      subroutine parse_arguments( &
        p_pdb1_filename, p_pdb1_bool, &
        p_pdb2_filename, p_pdb2_bool, &
        p_complexes_filename, p_complexes_bool, &
        p_nb_encounters, p_nb_encounters_bool, &
        p_datadist_filename, p_matrix_bool, &
        p_data_type, p_datatype_bool, &
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
        integer, intent(out) :: p_nb_encounters
        logical, intent(out) :: p_nb_encounters_bool
        character*128, intent(out) :: p_datadist_filename
        logical, intent(out) :: p_matrix_bool
        character*128, intent(out) :: p_data_type
        logical, intent(out) :: p_datatype_bool
        character*128, dimension(:), allocatable, intent(out) :: p_arr_atoms1a
        logical, intent(out) :: p_arr_atoms1a_bool
        character*128, dimension(:), allocatable, intent(out) :: p_arr_atoms1b
        logical, intent(out) :: p_arr_atoms1b_bool
        character*128, dimension(:), allocatable, intent(out) :: p_arr_atoms2a
        logical, intent(out) :: p_arr_atoms2a_bool
        character*128, dimension(:), allocatable, intent(out) :: p_arr_atoms2b
        logical, intent(out) :: p_arr_atoms2b_bool
        logical, intent(out) :: p_help_bool

        integer :: nb_argument, count_arg, n_atoms, i, ios
        character*128 :: argument

        nb_argument = command_argument_count()
        if ( nb_argument == 0 ) then
          print *, "[make_data] ERROR: No arguments parsed."
          print *, "[make_data] For usage, run: ./make_data -help"
          call print_help("main")
          STOP 1
        end if

        p_nb_encounters = 0
        p_pdb1_bool = .false.
        p_pdb2_bool = .false.
        p_complexes_bool = .false.
        p_matrix_bool = .false.
        p_datatype_bool = .false.
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
          else if ( trim(argument) == "-nb_encounters" ) then
            p_nb_encounters_bool = .true.
            call getarg( count_arg+1, argument )
            read(argument, *, IOSTAT=ios) p_nb_encounters
            if (ios /= 0) then
              print *, '[make_data] ERROR: Integer expected for the -nb_encounters argument.'
            end if
            count_arg = count_arg + 1
          else if ( trim(argument) == "-input" ) then
            p_matrix_bool = .true.
            call getarg( count_arg+1, argument )
            p_datadist_filename = trim(argument)
            count_arg = count_arg + 1
          else if ( trim(argument) == "-data_type" ) then
            p_datatype_bool = .true.
            call getarg( count_arg+1, argument )
            p_data_type = trim(argument)
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

            do i = 1, n_atoms
              call getarg( count_arg+i, argument )
              p_arr_atoms1a(i) = trim(argument)
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

            do i = 1, n_atoms
              call getarg( count_arg+i, argument )
              p_arr_atoms1b(i) = trim(argument)
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

            do i = 1, n_atoms
              call getarg( count_arg+i, argument )
              p_arr_atoms2a(i) = trim(argument)
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

            do i = 1, n_atoms
              call getarg( count_arg+i, argument )
              p_arr_atoms2b(i) = trim(argument)
            end do

            count_arg = count_arg + n_atoms
          else if ( trim(argument) == "-help" ) then
            p_help_bool = .true.
          else
            print *, "[make_data] ERROR: Argument '", trim(argument), "' not recognized."
            print *, "[make_data] For usage, run: ./make_data -help"
            STOP 1
          end if

          count_arg = count_arg + 1
        end do

      end subroutine parse_arguments

    !> Print usage/help text for `make_data`.
    !!
    !! Provides short usage examples for the main program and for each
    !! supported matrix type.
    subroutine print_help(help_option)

      implicit none
      character(len=*), intent(in) :: help_option

      !STOP 1

      if (trim(help_option) == "main") then
        print *, ""
        print *, "This program receives as inputs the pdb files, atom groups (either indexes or atomnames),", &
        " encounter complexes file, number of encounters, data type and output file."
        print *, "If you want instructions on how to create a specific type of data, type:"
        print *, "./make_data -data_type <type> -help"
        print *, ""
        print *, "Eg.: ./make_data -data_type rmsd -help"
        print *, ""
        print *, "The data types currently available are:"
        print *, "* rmsd"
        print *, "* z_coord"
        print *, "* atoms_dist"
        print *, "* 2D_angle"
        print *, "* 3D_angle"
        print *, ""
        print *, "General required options: -pdb2 -atoms2 -complexes -data_type -input"
        print *, ""
        STOP
      else if (trim(help_option) == "rmsd") then
        print *, ""
        print *, "To generate the rsmd matrix you will need to give pdb of solute 2, one group of ", &
        "atoms (either indexes or atomnames), the encounter complexes file, number of encounters, ", &
         "data type and output file."
        print *, ""
        print *, "Eg: ./make_data -pdb2 p2_noh.pdb -atoms2 C N O CA -complexes assoc_complexes ", &
        "-data_type rmsd -input matrix_rmsd.txt"
        print *, ""
        STOP
      else if (trim(help_option) == "z_coord") then
        print *, ""
        print *, "To generate z_coord data you will need to give the pdb of solute 2,", &
        "one group of atoms (either indexes or atomnames), the encounter complexes file, ", &
        "number of encounters, data type and output file."
        print *, ""
        print *, "Eg.: ./make_data -pdb2 p2_noh.pdb -atoms2 Cu -complexes assoc_complexes ", &
        "-data_type z_coord -input array_z.txt"
        print *, ""
        STOP
      else if (trim(help_option) == "atoms_dist") then
        print *, ""
        print *, "To generate atoms_dist data you will need to give two pdbs, two groups of atoms ", &
        "(either indexes or atomnames), the encounter complexes file, number of encounters, data type, ", &
        "and output file."
        print *, ""
        print *, "Eg.: ./make_data -pdb1 p1_noh.pdb -pdb2 p2_noh.pdb -complexes ", &
        "assoc_complexes -atoms1 C1 C4 -atoms2 Cu -data_type atoms_dist -input array_dist.txt"
        print *, ""
        STOP
      else if (trim(help_option) == "2D_angle" .or. trim(help_option) == "3D_angle") then
        print *, ""
        print *, "To generate angle data you will need to give two pdbs, four groups of atoms ", &
        "(either indexes or atomnames), the encounter complexes file, number of encounters, data ", & 
        "type and output file."
        print *, ""
        print *, "Eg.: ./make_data -pdb1 p1_noh.pdb -pdb2 p2_noh.pdb -complexes ", &
        "assoc_complexes -atoms1a 1 -atoms1b 411 -atoms2a Cu -atoms2b 4 5 6 7 8 -data_type 3D_angle ", &
        "-input array_angle.txt"
        print *, ""
        print *, "For this specific matrix type you need to be aware that a center ", &
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
      integer :: n, i, j, ios
      integer :: atom_index, atoms_count
      logical :: read_atomname, atom_found

      read_atomname = .false.

      n = size(atoms_array)


      read(atoms_array(1), *, IOSTAT=ios) atom_index

      if (ios /= 0) then
        read_atomname = .true.
      end if

      
      if (read_atomname) then

        allocate(array_atomnames(n))

        do j = 1, n
          read(atoms_array(j), *, IOSTAT=ios) atomname
          array_atomnames(j) = atomname
        end do
        
        atoms_count = 0
        do j = 1, n
          atom_found = .false.
          !print *, array_atomnames(j)
          do i = 1, tot_atoms
            if ( trim(adjustl(pdb % atoms(i) % name)) == &
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
          do i = 1, tot_atoms
            if ( trim(adjustl(pdb % atoms(i) % name)) == &
            trim(adjustl(array_atomnames(j))) ) then
              atoms_count = atoms_count + 1
              atoms_coord(atoms_count, :) = pdb % atoms(i) % coord(:)
            end if
          end do
        end do
      
      else
      
        allocate(atoms_coord(n, 3))
        do j = 1, n
          read(atoms_array(j), *, IOSTAT=ios) i
          atoms_coord(j, :) = pdb % atoms(i) % coord(:)
        end do

      end if

    end subroutine read_atoms_coord
      
    

END PROGRAM main