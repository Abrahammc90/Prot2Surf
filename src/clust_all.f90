!> Command-line utility to build matrices and perform clustering in memory.
!!
!! This program combines matrix generation and clustering without writing
!! the matrix to disk. Useful for large matrices that would exceed disk space.
!!
!! @author Abraham Muñiz-Chicharro
program clust_all
  USE read_input
  USE mod_matrix
  USE mod_array
  USE mod_clust_algorithm
  USE mod_assoc

  implicit none
      
    character*128 :: pdb1_filename, pdb2_filename, complexes_filename, sorted_complexes_filename
    character*128 :: data_type, datadist_filename
    character*128 :: linkage_type, output_name
    character*128 :: argument
    character*128, dimension(:), allocatable :: arr_atoms1a, arr_atoms1b, arr_atoms2a, arr_atoms2b

    integer :: tot_atoms2, tot_encounters, tot_coords
    integer :: tot_chains2, tot_residues2
    integer :: tot_atoms1, tot_residues1, tot_chains1
    integer :: nb_encounters
    type ( type_pdb_file ) :: pdb1, pdb2
    type ( type_assoc_file ) :: complexes

    integer :: i
    real (kind=8), dimension(:, :), allocatable :: distmatrix
    real (kind=8), dimension(:), allocatable :: distarray
    real (kind=8), dimension(:, :), allocatable :: atoms1a_coords, atoms1b_coords
    real (kind=8), dimension(:, :), allocatable :: atoms2a_coords, atoms2b_coords
    real (kind=8), dimension(3) :: cog1a, cog1b, cog2a, cog2b
    integer :: nb_argument, count_arg, n_atoms
    integer :: ios
    logical :: pdb1_bool, pdb2_bool, complexes_bool, atoms1_bool, atoms2_bool, datadist_filename_bool
    logical :: datatype_bool, nb_encounters_bool, help_bool
    logical :: arr_atoms1a_bool, arr_atoms1b_bool, arr_atoms2a_bool, arr_atoms2b_bool
    logical :: linkage_bool, output_bool, use_cuda_bool

    nb_argument = 0
    nb_encounters = 0
    count_arg = 1
    pdb1_bool = .false.
    pdb2_bool = .false.
    complexes_bool = .false.
    atoms1_bool = .false.
    atoms2_bool = .false.
    datatype_bool = .false.
    nb_encounters_bool = .false.
    arr_atoms1a_bool = .false.
    arr_atoms1b_bool = .false.
    arr_atoms2a_bool = .false.
    arr_atoms2b_bool = .false.
    help_bool = .false.
    linkage_bool = .false.
    output_bool = .false.
    use_cuda_bool = .false.
    argument = ""
    linkage_type = "min"
    output_name = "clust_all_output"

    nb_argument = command_argument_count()
    if ( nb_argument == 0 ) then
      print *, ""
      print *, "ERROR. No arguments were parsed"
      call print_help("main")
      STOP 1
    end if

    do while ( count_arg <= nb_argument )
      call getarg( count_arg, argument )
      
      if ( trim(argument) == "-pdb1" ) then
        pdb1_bool = .true.
        call getarg( count_arg+1, argument )
        pdb1_filename = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-pdb2" ) then
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
          print *, "ERROR. Integer expected for the -nb_encounters argument."
        end if
        count_arg = count_arg + 1
      else if ( trim(argument) == "-data_type" ) then
        datatype_bool = .true.
        call getarg( count_arg+1, argument )
        data_type = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-datadist" ) then
        datadist_filename_bool = .true.
        call getarg( count_arg+1, argument )
        datadist_filename = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-linkage" ) then
        linkage_bool = .true.
        call getarg( count_arg+1, argument )
        linkage_type = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-output_name" ) then
        output_bool = .true.
        call getarg( count_arg+1, argument )
        output_name = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-cuda" ) then
        use_cuda_bool = .true.
      else if ( trim(argument) == "-atoms1" .or. trim(argument) == "-atoms1a" ) then
        arr_atoms1a_bool = .true.
        
        call getarg( count_arg+1, argument )
        n_atoms = 0
        do while (argument(1:1) /= "-" .and. &
          count_arg + n_atoms <= nb_argument)
          n_atoms = n_atoms + 1
          call getarg( count_arg+n_atoms, argument )
        end do
        n_atoms = n_atoms-1
        allocate(arr_atoms1a(n_atoms))

        do i = 1, n_atoms
          call getarg( count_arg+i, argument )
          arr_atoms1a(i) = trim(argument)
        end do

        count_arg = count_arg + n_atoms

      else if ( trim(argument) == "-atoms1b" ) then
        arr_atoms1b_bool = .true.

        call getarg( count_arg+1, argument )
        n_atoms = 0
        do while (argument(1:1) /= "-" .and. &
          count_arg + n_atoms <= nb_argument)
          n_atoms = n_atoms + 1
          call getarg( count_arg+n_atoms, argument )
        end do
        n_atoms = n_atoms-1
        allocate(arr_atoms1b(n_atoms))

        do i = 1, n_atoms
          call getarg( count_arg+i, argument )
          arr_atoms1b(i) = trim(argument)
        end do

        count_arg = count_arg + n_atoms

      else if ( trim(argument) == "-atoms2" .or. trim(argument) == "-atoms2a" ) then
        arr_atoms2a_bool = .true.

        call getarg( count_arg+1, argument )
        n_atoms = 0
        do while (argument(1:1) /= "-" .and. &
          count_arg + n_atoms <= nb_argument)
          n_atoms = n_atoms + 1
          call getarg( count_arg+n_atoms, argument )
        end do
        n_atoms = n_atoms-1
        allocate(arr_atoms2a(n_atoms))

        do i = 1, n_atoms
          call getarg( count_arg+i, argument )
          arr_atoms2a(i) = trim(argument)
        end do

        count_arg = count_arg + n_atoms

      else if ( trim(argument) == "-atoms2b" ) then
        arr_atoms2b_bool = .true.

        call getarg( count_arg+1, argument )
        n_atoms = 0
        do while (argument(1:1) /= "-" .and. &
          count_arg + n_atoms <= nb_argument)
          n_atoms = n_atoms + 1
          call getarg( count_arg+n_atoms, argument )
        end do
        n_atoms = n_atoms-1
        allocate(arr_atoms2b(n_atoms))

        do i = 1, n_atoms
          call getarg( count_arg+i, argument )
          arr_atoms2b(i) = trim(argument)
        end do

        count_arg = count_arg + n_atoms

      else if ( trim(argument) == "-help" ) then
        help_bool = .true.
      else
        print *, ""
        print *,  "ERROR. Argument ", trim(argument), " not recognized."
        print *, "For more information, please use the -help option:"
        print *, "./clust_all -help"
        STOP 1
      end if

      count_arg = count_arg + 1

    end do

    if ( help_bool ) then
      if (datatype_bool) then
        call print_help(data_type)
      else
        call print_help("main")
      end if
    end if

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
      print *, "ERROR. Output array filename not given."
      print *, "Please, provide the output array filename with the '-input' option"
      print *, "For more information, please use the -help option:"
      print *, "./clust_all -help"
      STOP 1
    end if

    !Reads PDB
    call read_pdb(pdb2, pdb2_filename, tot_atoms2, tot_residues2, tot_chains2)

    !Reads Complexes
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

    allocate(distmatrix (nb_encounters, nb_encounters))
    allocate(distarray (nb_encounters))
    sorted_complexes_filename = "sorted_" // trim(complexes_filename)

    print *, 'Building ', trim(data_type), ' matrix in memory'

    if (trim(adjustl(data_type)) == "z_coord") then

      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call array_z_coord(complexes, distarray, nb_encounters, 1, complexes % xc1, complexes % xc2, &
      complexes % trans_vector, complexes % rot1, complexes % rot2, atoms2a_coords)

      call write_array(distarray, datadist_filename)
      call write_complexes(complexes, nb_encounters, complexes_filename)

    else if (trim(adjustl(data_type)) == "rmsd") then

      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call matrix_rmsd(distmatrix, nb_encounters, tot_coords, &
      complexes % xc1, complexes % xc2, complexes % trans_vector, &
      complexes % rot1, complexes % rot2, atoms2a_coords)

    else if (trim(adjustl(data_type)) == "3D_angle" .or. trim(adjustl(data_type)) == "2D_angle") then

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
        call array_angle(complexes, distarray, nb_encounters, 2, &
        complexes % xc1, complexes % xc2, complexes % trans_vector, complexes % rot1, complexes % rot2, &
        cog1a, cog1b, cog2a, cog2b, 2)
      else if (trim(adjustl(data_type)) == "3D_angle") then
        call array_angle(complexes, distarray, nb_encounters, 2, &
        complexes % xc1, complexes % xc2, complexes % trans_vector, complexes % rot1, complexes % rot2, &
        cog1a, cog1b, cog2a, cog2b, 3)
      end if

      call write_array(distarray, datadist_filename)
      call write_complexes(complexes, nb_encounters, sorted_complexes_filename)

    else if (trim(adjustl(data_type)) == "atoms_dist") then

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

      call read_pdb(pdb1, pdb1_filename, tot_atoms1, tot_residues1, tot_chains1)
      call read_atoms_coord(arr_atoms1a, atoms1a_coords, tot_atoms1, pdb1, pdb1_filename)
      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call array_atoms_dist(complexes, distarray, nb_encounters, 1, &
                            complexes % xc1, complexes % xc2, complexes % trans_vector, &
                            complexes % rot1, complexes % rot2, atoms1a_coords, atoms2a_coords)

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

    ! Perform clustering directly on the in-memory matrix
    if (trim(adjustl(data_type)) == "rmsd") then
      call linkage_clustering_from_matrix(distmatrix, nb_encounters, linkage_type, output_name, &
                              complexes, use_cuda=use_cuda_bool)
    else
      call linkage_clustering_from_array(distarray, nb_encounters, linkage_type, output_name, &
                              complexes)
    end if

    print *, 'Clustering complete'
  

  contains
    !> Print usage/help text for `clust_all`.
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
        print *, "  -matrix_type <type>    Type: z_coord, rmsd, atoms_dist, 2D_angle, 3D_angle"
        print *, ""
        print *, "Optional arguments:"
        print *, "  -pdb1 <file>           PDB file for solute 1 (required for angle/atoms_dist)"
        print *, "  -atoms1a <atoms>       First group of atoms from pdb1"
        print *, "  -atoms1b <atoms>       Second group of atoms from pdb1"
        print *, "  -atoms2b <atoms>       Second group of atoms from pdb2"
        print *, "  -array <file>          Output array filename"
        print *, "  -nb_encounters <n>     Number of encounters to process"
        print *, "  -linkage <type>        Linkage type: min, max, mean (default: min)"
        print *, "  -output_name <name>    Output basename (default: clust_all_output)"
        print *, "  -cuda                  Use CUDA acceleration"
        print *, "  -help                  Show this help"
        print *, ""
        print *, "Example:"
        print *, "  ./clust_all -pdb2 p2.pdb -atoms2 CU -complexes assoc_complexes \\"
        print *, "    -matrix_type z_coord -array array_z.txt -linkage min \\"
        print *, "    -output_name Cu_z_clusters"
        print *, ""
        STOP
      else if (trim(help_option) == "rmsd") then
        print *, ""
        print *, "RMSD matrix: measures RMSD between transformed solute positions"
        print *, ""
        print *, "Required: -pdb2, -atoms2, -complexes, -matrix_type rmsd"
        print *, "Optional: -nb_encounters, -linkage, -output_name, -cuda"
        print *, ""
        STOP
      else if (trim(help_option) == "z_coord") then
        print *, ""
        print *, "Z-coordinate matrix: measures z-distance differences"
        print *, ""
        print *, "Required: -pdb2, -atoms2, -complexes, -matrix_type z_coord, -array"
        print *, "Optional: -nb_encounters, -linkage, -output_name, -cuda"
        print *, ""
        STOP
      else if (trim(help_option) == "atoms_dist") then
        print *, ""
        print *, "Atoms distance matrix: measures inter-atomic distances"
        print *, ""
        print *, "Required: -pdb1, -pdb2, -atoms1, -atoms2, -complexes,"
        print *, "          -matrix_type atoms_dist, -array"
        print *, "Optional: -nb_encounters, -linkage, -output_name, -cuda"
        print *, ""
        STOP
      else if (trim(help_option) == "2D_angle" .or. trim(help_option) == "3D_angle") then
        print *, ""
        print *, "Angle matrix: measures angles between vectors defined by atom groups"
        print *, ""
        print *, "Required: -pdb1, -pdb2, -atoms1a, -atoms1b, -atoms2a, -atoms2b,"
        print *, "          -complexes, -matrix_type 2D_angle or 3D_angle, -array"
        print *, "Optional: -nb_encounters, -linkage, -output_name, -cuda"
        print *, ""
        STOP
      end if
      
    end subroutine print_help

    !> Read coordinates for a specified list of atoms (names or indices).
    subroutine read_atoms_coord(atoms_array, atoms_coord, tot_atoms, pdb, pdb_filename)
      implicit none

      character(len=*), dimension(:), intent(in) :: atoms_array
      character(len=*), intent(in) :: pdb_filename
      integer, intent(in) :: tot_atoms
      real(kind=8), dimension(:, :), allocatable, intent(out) :: atoms_coord
      type ( type_pdb_file ), intent(in) :: pdb

      character (len=4), dimension(:), allocatable :: array_atomnames
      character (len=4) :: atomname
      integer, dimension(:), allocatable :: array_atom_indexes
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

END PROGRAM clust_all
