!> \file clust_all.f90
!! \brief Build matrices/arrays and perform clustering in memory (no disk I/O).
!! Usage: ./clust_all -help
!!
!! This driver combines two traditionally separate phases:
!! (1) metric generation from encounter transforms and atom selections,
!! and (2) hierarchical clustering on the generated metric.
!!
!! Compared with split workflows (`make_data` then `clust`), this program
!! keeps intermediate numeric data in memory to reduce temporary files and
!! simplify large-batch runs.
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0
!! @date 2026-04-05
program clust_all
  USE read_input
  USE mod_matrix
  USE mod_array
  USE mod_clust_algorithm
  USE mod_assoc

  implicit none

  ! High-level workflow:
  ! 1) Parse CLI options and validate inputs for selected data_type.
  ! 2) Read PDB/complexes information.
  ! 3) Build data representation in memory (array or matrix).
  ! 4) Run hierarchical clustering directly from in-memory data.
  ! 5) Emit clustering outputs using -output_name basename.
      
    character*128 :: pdb1_filename, pdb2_filename, complexes_filename, sorted_complexes_filename
    character*128 :: data_type, datadist_filename
    character*128 :: linkage_type, output_name
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
    logical :: pdb1_bool, pdb2_bool, complexes_bool, datadist_filename_bool
    logical :: datatype_bool, nb_encounters_bool, help_bool
    logical :: arr_atoms1a_bool, arr_atoms1b_bool, arr_atoms2a_bool, arr_atoms2b_bool

    call parse_arguments( &
      pdb1_filename, pdb1_bool, &
      pdb2_filename, pdb2_bool, &
      complexes_filename, complexes_bool, &
      nb_encounters, nb_encounters_bool, &
      data_type, datatype_bool, &
      datadist_filename, datadist_filename_bool, &
      linkage_type, output_name, &
      arr_atoms1a, arr_atoms1a_bool, &
      arr_atoms1b, arr_atoms1b_bool, &
      arr_atoms2a, arr_atoms2a_bool, &
      arr_atoms2b, arr_atoms2b_bool, &
      help_bool)

    ! Help path exits early and bypasses input validation/processing.
    if ( help_bool ) then
      if (datatype_bool) then
        call print_help(data_type)
      else
        call print_help("main")
      end if
    end if

    ! Global input validation shared by all metric modes.
    ! Metric-specific validation is performed inside each data_type branch.
    if (.not. pdb2_bool ) then
      print *, "ERROR. PDB not provided."
      print *, "Please provide a pdb file with '-pdb2' option"
      print *, "For more information, please use the -help option:"
      print *, "./clust_all -help"
      STOP 1
    else if ( .not. arr_atoms2a_bool ) then 
      print *, "ERROR. PDB2 provided but atoms2 not provided."
      print *, "Please provide a group of atoms with '-atoms2' option"
      print *, "For more information, please use the -help option:"
      print *, "./clust_all -help"
      STOP 1
    else if ( .not. complexes_bool ) then
      print *, "ERROR. Encounter complexes file not provided."
      print *, "Please provide the complexes file with the '-complexes' option"
      print *, "For more information, please use the -help option:"
      print *, "./clust_all -help"
      STOP 1
    else if ( .not. datatype_bool ) then
      print *, "ERROR. Data type not given."
      print *, "Please, provide the data type with the '-data_type' option"
      print *, "For more information, please use the -help option:"
      print *, "./clust_all -help"
      STOP 1
    else if ( .not. datadist_filename_bool ) then
      print *, "ERROR. Data output filename not given."
      print *, "Please, provide the output filename with the '-datadist' option"
      print *, "For more information, please use the -help option:"
      print *, "./clust_all -help"
      STOP 1
    end if

    ! Read structural coordinates for solute 2 (always required).
    call read_pdb(pdb2, pdb2_filename, tot_atoms2, tot_residues2, tot_chains2)

    ! Read encounter transforms and metadata.
    call read_assoc(complexes, complexes_filename, tot_encounters)

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

    ! Build data in memory according to selected metric.
    print *, 'Building ', trim(data_type), ' matrix in memory'

    if (trim(adjustl(data_type)) == "z_coord") then

      ! z_coord mode:
      ! - select atom group in solute 2
      ! - transform coordinates per encounter
      ! - extract one z-based scalar per encounter

      allocate(distarray (nb_encounters))

      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call array_z_coord(complexes, distarray, nb_encounters, 1, complexes % xc1, complexes % xc2, &
      complexes % trans_vector, complexes % rot1, complexes % rot2, atoms2a_coords)

      ! z_coord creates an array (one value per encounter).
      call write_array(distarray, datadist_filename)
      call write_complexes(complexes, nb_encounters, complexes_filename)

    else if (trim(adjustl(data_type)) == "rmsd") then

      ! rmsd mode:
      ! - build a full pairwise matrix from transformed structures
      ! - clustering later runs on this matrix representation

      allocate(distmatrix (nb_encounters, nb_encounters))

      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call matrix_rmsd(distmatrix, nb_encounters, tot_coords, &
      complexes % xc1, complexes % xc2, complexes % trans_vector, &
      complexes % rot1, complexes % rot2, atoms2a_coords)

    else if (trim(adjustl(data_type)) == "3D_angle" .or. trim(adjustl(data_type)) == "2D_angle") then

      ! angle modes (2D/3D):
      ! - require two atom groups per solute (a,b)
      ! - compute COG for each group to define two vectors
      ! - compute one angle per encounter between transformed vectors

      if (.not. pdb1_bool) then
        print *, "ERROR. Just one PDB provided. It is necessary to provide two pdbs ", &
        "to generate the angle matrix."
        print *, "Please provide the two pdb files with '-pdb1' and '-pdb2' options."
        print *, "For more information, please use the -help option:"
        print *, "./clust_all -data_type angle -help"
        STOP 1
      else if (.not. arr_atoms1a_bool .or. .not. arr_atoms1b_bool .or. .not. arr_atoms2b_bool) then 
        print *, "ERROR. Not all necessary group of atoms given"
        print *, "Please provide all group of atoms with '-atoms1a', '-atoms1b', '-atoms2a' and '-atoms2b' options"
        print *, "For more information, please use the -help option:"
        print *, "./clust_all -data_type angle -help"
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
        ! 2D projection angle (xy plane projection).
        call array_angle(complexes, distarray, nb_encounters, 2, &
        complexes % xc1, complexes % xc2, complexes % trans_vector, complexes % rot1, complexes % rot2, &
        cog1a, cog1b, cog2a, cog2b, 2)
      else if (trim(adjustl(data_type)) == "3D_angle") then
        ! Full 3D angle.
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
      ! - write scalar array used by array-based clustering path

      if (.not. pdb1_bool) then
        print *, "ERROR. Just one PDB provided. It is necessary to provide two pdbs ", &
        "to generate the atoms_dist matrix."
        print *, "Please provide the two pdb files with '-pdb1' and '-pdb2' options"
        print *, "For more information, please use the -help option:"
        print *, "./clust_all -data_type atoms_dist -help"
        STOP 1
      else if (.not. arr_atoms1a_bool ) then 
        print *, "ERROR. Not all necessary group of atoms given"
        print *, "Please provide the group of atoms with '-atoms1' and '-atoms2' options"
        print *, "For more information, please use the -help option:"
        print *, "./clust_all -data_type atoms_dist -help"
        STOP 1
      end if

      allocate(distarray (nb_encounters))

      call read_pdb(pdb1, pdb1_filename, tot_atoms1, tot_residues1, tot_chains1)
      call read_atoms_coord(arr_atoms1a, atoms1a_coords, tot_atoms1, pdb1, pdb1_filename)
      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call array_atoms_dist(complexes, distarray, nb_encounters, 1, &
                            complexes % xc1, complexes % xc2, complexes % trans_vector, &
                            complexes % rot1, complexes % rot2, atoms1a_coords, atoms2a_coords)

      ! atoms_dist creates an array (minimum distance per encounter).
      call write_array(distarray, datadist_filename)
      call write_complexes(complexes, nb_encounters, sorted_complexes_filename)

    else
      print *, 'ERROR: ', trim(adjustl(data_type)), ' matrix not implemented'
      print *, "To know more about the type of matrices available, please run:"
      print *, "./clust_all -help"
      STOP
    end if

    print *, 'Data generation complete'
    print *, 'Starting clustering with linkage type: ', trim(linkage_type)
    write(*,*)

    ! Perform clustering directly on in-memory data representation.
    ! - rmsd path uses matrix-based clustering
    ! - all other data types use array-based clustering
    if (trim(adjustl(data_type)) == "rmsd") then
      ! Matrix linkage clustering path.
      call linkage_clustering_from_matrix(distmatrix, nb_encounters, linkage_type, output_name, &
                              complexes)
    else
      ! Array linkage clustering path.
      call linkage_clustering_from_array(distarray, nb_encounters, linkage_type, output_name, &
                              complexes)
    end if

    print *, 'Clustering complete'
  

  contains
    !> Parse command-line arguments for `clust_all`.
    !!
    !! Parses all CLI flags and fills program configuration values.
    !! Atom-group flags consume multiple tokens until the next flag.
    !! @param[out] p_pdb1_filename         Value parsed from `-pdb1`
    !! @param[out] p_pdb1_bool             Whether `-pdb1` was provided
    !! @param[out] p_pdb2_filename         Value parsed from `-pdb2`
    !! @param[out] p_pdb2_bool             Whether `-pdb2` was provided
    !! @param[out] p_complexes_filename    Value parsed from `-complexes`
    !! @param[out] p_complexes_bool        Whether `-complexes` was provided
    !! @param[out] p_nb_encounters         Value parsed from `-nb_encounters`
    !! @param[out] p_nb_encounters_bool    Whether `-nb_encounters` was provided
    !! @param[out] p_data_type             Value parsed from `-data_type`
    !! @param[out] p_datatype_bool         Whether `-data_type` was provided
    !! @param[out] p_datadist_filename     Value parsed from `-datadist`
    !! @param[out] p_datadist_filename_bool Whether `-datadist` was provided
    !! @param[out] p_linkage_type          Value parsed from `-linkage` (default `min`)
    !! @param[out] p_output_name           Value parsed from `-output_name` (default `clust_all_output`)
    !! @param[out] p_arr_atoms1a           Atom list parsed from `-atoms1` or `-atoms1a`
    !! @param[out] p_arr_atoms1a_bool      Whether `-atoms1` or `-atoms1a` was provided
    !! @param[out] p_arr_atoms1b           Atom list parsed from `-atoms1b`
    !! @param[out] p_arr_atoms1b_bool      Whether `-atoms1b` was provided
    !! @param[out] p_arr_atoms2a           Atom list parsed from `-atoms2` or `-atoms2a`
    !! @param[out] p_arr_atoms2a_bool      Whether `-atoms2` or `-atoms2a` was provided
    !! @param[out] p_arr_atoms2b           Atom list parsed from `-atoms2b`
    !! @param[out] p_arr_atoms2b_bool      Whether `-atoms2b` was provided
    !! @param[out] p_help_bool             Whether `-help` was provided
    subroutine parse_arguments( &
      p_pdb1_filename, p_pdb1_bool, &
      p_pdb2_filename, p_pdb2_bool, &
      p_complexes_filename, p_complexes_bool, &
      p_nb_encounters, p_nb_encounters_bool, &
      p_data_type, p_datatype_bool, &
      p_datadist_filename, p_datadist_filename_bool, &
      p_linkage_type, p_output_name, &
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
      character*128, intent(out) :: p_data_type
      logical, intent(out) :: p_datatype_bool
      character*128, intent(out) :: p_datadist_filename
      logical, intent(out) :: p_datadist_filename_bool
      character*128, intent(out) :: p_linkage_type
      character*128, intent(out) :: p_output_name
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
        print *, ""
        print *, "ERROR. No arguments were parsed"
        call print_help("main")
        STOP 1
      end if

      p_nb_encounters = 0
      p_pdb1_bool = .false.
      p_pdb2_bool = .false.
      p_complexes_bool = .false.
      p_nb_encounters_bool = .false.
      p_datatype_bool = .false.
      p_datadist_filename_bool = .false.
      p_arr_atoms1a_bool = .false.
      p_arr_atoms1b_bool = .false.
      p_arr_atoms2a_bool = .false.
      p_arr_atoms2b_bool = .false.
      p_help_bool = .false.
      p_linkage_type = "min"
      p_output_name = "clust_all_output"
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
            print *, "ERROR. Integer expected for the -nb_encounters argument."
          end if
          count_arg = count_arg + 1
        else if ( trim(argument) == "-data_type" ) then
          p_datatype_bool = .true.
          call getarg( count_arg+1, argument )
          p_data_type = trim(argument)
          count_arg = count_arg + 1
        else if ( trim(argument) == "-datadist" ) then
          p_datadist_filename_bool = .true.
          call getarg( count_arg+1, argument )
          p_datadist_filename = trim(argument)
          count_arg = count_arg + 1
        else if ( trim(argument) == "-linkage" ) then
          call getarg( count_arg+1, argument )
          p_linkage_type = trim(argument)
          count_arg = count_arg + 1
        else if ( trim(argument) == "-output_name" ) then
          call getarg( count_arg+1, argument )
          p_output_name = trim(argument)
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
          print *, ""
          print *,  "ERROR. Argument ", trim(argument), " not recognized."
          print *, "For more information, please use the -help option:"
          print *, "./clust_all -help"
          STOP 1
        end if

        count_arg = count_arg + 1
      end do

    end subroutine parse_arguments

    !> Print usage/help text for `clust_all`.
    !!
    !! `help_option` controls whether general usage or mode-specific
    !! requirements are printed.
    subroutine print_help(help_option)

      implicit none
      character(len=*), intent(in) :: help_option

      if (trim(help_option) == "main") then
        print *, ""
        print *, "This program generates a matrix in memory and performs clustering"
        print *, "without writing the matrix to disk. Useful for large matrices."
        print *, ""
        print *, "Required arguments:"
        print *, "  -pdb2 <file>           PDB file for solute 2"
        print *, "  -atoms2 <atoms>        Group of atoms (names or indices)"
        print *, "  -complexes <file>      Encounter complexes file"
        print *, "  -data_type <type>      Type: z_coord, rmsd, atoms_dist, 2D_angle, 3D_angle"
        print *, "  -datadist <file>       Output data file (array for non-rmsd, matrix for rmsd)"
        print *, ""
        print *, "Optional arguments:"
        print *, "  -pdb1 <file>           PDB file for solute 1 (required for angle/atoms_dist)"
        print *, "  -atoms1a <atoms>       First group of atoms from pdb1"
        print *, "  -atoms1b <atoms>       Second group of atoms from pdb1"
        print *, "  -atoms2b <atoms>       Second group of atoms from pdb2"
        print *, "  -nb_encounters <n>     Number of encounters to process"
        print *, "  -linkage <type>        Linkage type: min, max, mean (default: min)"
        print *, "  -output_name <name>    Output basename for clustering outputs (default: clust_all_output)"
        print *, "  -help                  Show this help"
        print *, ""
        print *, "Example:"
        print *, "  ./clust_all -pdb2 p2.pdb -atoms2 CU -complexes assoc_complexes \\"
        print *, "    -data_type z_coord -datadist array_z.txt -linkage min \\" 
        print *, "    -output_name Cu_z_clusters"
        print *, ""
        STOP
      else if (trim(help_option) == "rmsd") then
        print *, ""
        print *, "RMSD matrix: measures RMSD between transformed solute positions"
        print *, ""
        print *, "Required: -pdb2, -atoms2, -complexes, -data_type rmsd, -datadist <file>"
        print *, "Optional: -nb_encounters, -linkage, -output_name"
        print *, ""
        STOP
      else if (trim(help_option) == "z_coord") then
        print *, ""
        print *, "Z-coordinate matrix: measures z-distance differences"
        print *, ""
        print *, "Required: -pdb2, -atoms2, -complexes, -data_type z_coord, -datadist <file>"
        print *, "Optional: -nb_encounters, -linkage, -output_name"
        print *, ""
        STOP
      else if (trim(help_option) == "atoms_dist") then
        print *, ""
        print *, "Atoms distance matrix: measures inter-atomic distances"
        print *, ""
        print *, "Required: -pdb1, -pdb2, -atoms1, -atoms2, -complexes,"
        print *, "          -data_type atoms_dist, -datadist <file>"
        print *, "Optional: -nb_encounters, -linkage, -output_name"
        print *, ""
        STOP
      else if (trim(help_option) == "2D_angle" .or. trim(help_option) == "3D_angle") then
        print *, ""
        print *, "Angle matrix: measures angles between vectors defined by atom groups"
        print *, ""
        print *, "Required: -pdb1, -pdb2, -atoms1a, -atoms1b, -atoms2a, -atoms2b,"
        print *, "          -complexes, -data_type 2D_angle or 3D_angle, -datadist <file>"
        print *, "Optional: -nb_encounters, -linkage, -output_name"
        print *, ""
        STOP
      end if
      
    end subroutine print_help

    !> Read coordinates for a specified list of atoms (names or indices).
    !!
    !! Accepts either atom indices or atom names. If atom names are
    !! provided, all matching atoms in the PDB are returned (in order).
    !! @param[in]  atoms_array  Array of atom name strings or index strings
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

      ! Detect whether user provided numeric atom indices or atom names.
      if (ios /= 0) then
        read_atomname = .true.
      end if

      if (read_atomname) then

        ! Name-based mode:
        ! 1) copy requested names
        ! 2) count all matching atoms in pdb (including repeated names)
        ! 3) allocate output to exact matched size
        ! 4) fill output coordinates preserving request-order grouping

        allocate(array_atomnames(n))

        do j = 1, n
          read(atoms_array(j), *, IOSTAT=ios) atomname
          array_atomnames(j) = atomname
        end do
        
        atoms_count = 0
        do j = 1, n
          atom_found = .false.
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

        ! Index-based mode:
        ! directly map each requested atom index to one coordinate row.
      
        allocate(atoms_coord(n, 3))
        do j = 1, n
          read(atoms_array(j), *, IOSTAT=ios) i
          atoms_coord(j, :) = pdb % atoms(i) % coord(:)
        end do

      end if

    end subroutine read_atoms_coord

END PROGRAM clust_all
