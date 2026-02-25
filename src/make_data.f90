!> Command-line utility to build encounter-based matrices.
!!
!! Usage: see `print_help` below. This program reads PDB(s) and an
!! association/complex file and generates matrices (rmsd, z_coord,
!! atoms_dist, angle) and arrays derived from encounter transforms.
!!
!! @author Abraham Muñiz-Chicharro
program main
  USE read_input
  USE mod_matrix
  USE mod_array


  implicit none
      

    character*128 :: pdb1_filename, pdb2_filename, complexes_filename, sorted_complexes_filename
    character*128 :: data_type, datadist_filename
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
    !real (kind=8), dimension(:, :), allocatable :: residues_cog, residue_crds
    integer :: ios
    logical :: pdb1_bool, pdb2_bool, complexes_bool, atoms1_bool, atoms2_bool
    logical :: matrix_bool, datatype_bool, nb_encounters_bool, help_bool
    logical :: arr_atoms1a_bool, arr_atoms1b_bool, arr_atoms2a_bool, arr_atoms2b_bool

    nb_argument = 0
    nb_encounters = 0
    count_arg = 1
    pdb1_bool = .false.
    pdb2_bool = .false.
    complexes_bool = .false.
    atoms1_bool = .false.
    atoms2_bool = .false.
    matrix_bool = .false.
    datatype_bool = .false.
    nb_encounters_bool = .false.
    arr_atoms1a_bool = .false.
    arr_atoms1b_bool = .false.
    arr_atoms2a_bool = .false.
    arr_atoms2b_bool = .false.
    help_bool = .false.
    argument = ""
    !arr_atoms1 = ""


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
          print *, "ERROR. Integer expected for the -nb_encounters ", &
          "argument."
        end if
        count_arg = count_arg + 1
      else if ( trim(argument) == "-input" ) then
        matrix_bool = .true.
        call getarg( count_arg+1, argument )
        datadist_filename = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-data_type" ) then
        datatype_bool = .true.
        call getarg( count_arg+1, argument )
        data_type = trim(argument)
        count_arg = count_arg + 1

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
        print *, "./make_data -help"
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
      print *, "./make_data -help"
      STOP 1
    else if ( .not. arr_atoms2a_bool ) then 
      print *, "ERROR. PDB2 provided but atoms2 not provided."
      print *, "Please provide a group of atoms with '-atoms2' option"
      print *, "For more information, please use the -help option:"
      print *, "./make_data -help"
      STOP 1
    else if ( .not. complexes_bool ) then
      print *, "ERROR. Encounter complexes file not provided."
      print *, "Please provide the complexes file with the '-complexes' option"
      print *, "For more information, please use the -help option:"
      print *, "./make_data -help"
      STOP 1
    else if ( .not. matrix_bool ) then
      print *, "ERROR. Matrix output filename not provided."
      print *, "Please, provide the matrix outupt name with the '-matrix' option"
      print *, "For more information, please use the -help option:"
      print *, "./make_data -help"
      STOP 1
    else if ( .not. datatype_bool ) then
      print *, "ERROR. Matrix type not given."
      print *, "Please, provide the matrix type with the '-matrix_type' option"
      print *, "For more information, please use the -help option:"
      print *, "./make_data -help"
      STOP 1
    end if


    !Reads PDB
    call read_pdb(pdb2, pdb2_filename, tot_atoms2, tot_residues2, tot_chains2)
    !call read_atoms_coord(arr_atoms2, atoms2_coords, tot_atoms2, pdb2, pdb2_filename)

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


    sorted_complexes_filename = "sorted_" // trim(complexes_filename)

    print *, 'Building ', trim(data_type), ' matrix'

    if (trim(adjustl(data_type)) == "z_coord") then

      allocate(distarray (nb_encounters))
      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call array_z_coord(complexes, distarray, nb_encounters, 1, complexes % xc1, complexes % xc2, &
      complexes % trans_vector, complexes % rot1, complexes % rot2, atoms2a_coords)

      call write_array(distarray, datadist_filename)
      call write_complexes(complexes, nb_encounters, complexes_filename)

    else if (trim(adjustl(data_type)) == "rmsd") then

      allocate(distmatrix (nb_encounters, nb_encounters))
      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call matrix_rmsd(distmatrix, nb_encounters, tot_coords, &
      complexes % xc1, complexes % xc2, complexes % trans_vector, &
      complexes % rot1, complexes % rot2, atoms2a_coords)

      call write_matrix(distmatrix, datadist_filename)

    else if (trim(adjustl(data_type)) == "3D_angle" .or. trim(adjustl(data_type)) == "2D_angle") then

      if (.not. pdb1_bool) then
        print *, "ERROR. Just one PDB provided. It is necessary to provide two pdbs ", &
        "to generate the angle matrix."
        print *, "Please provide the two pdb files with '-pdb1' and '-pdb2' options."
        print *, "For more information, please use the -help option:"
        print *, "./make_data -matrix angle -help"
        STOP 1
      else if (.not. arr_atoms1a_bool .or. .not. arr_atoms1b_bool .or. .not. arr_atoms2b_bool) then 
        print *, "ERROR. Not all necessary group of atoms given"
        print *, "Please provide all group of atoms with '-atoms1a', '-atoms1b', '-atoms2a' and '-atoms2b' options"
        print *, "For more information, please use the -help option:"
        print *, "./make_data -matrix angle -help"
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
        "to generate the angle matrix."
        print *, "Please provide the two pdb files with '-pdb1' and '-pdb2' options"
        print *, "For more information, please use the -help option:"
        print *, "./make_data -matrix atoms_dist -help"
        STOP 1
      else if (.not. arr_atoms1a_bool ) then 
        print *, "ERROR. Not all necessary group of atoms given"
        print *, "Please provide the group of atoms with '-atoms1' and '-atoms2' options"
        print *, "For more information, please use the -help option:"
        print *, "./make_data -matrix atoms_dist -help"
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
        print *, "This program receives as inputs the pdb files, atom groups (either indexes or atomnames)", &
        "encounter complexes file, number of encounters to make the matrix distance, ", &
        "matrix type, array output file and matrix output file."
        print *, "If you want instructions on how to create a specific type of matrix, type:"
        print *, "./make_data -matrix_type <type of matrix> -help"
        print *, ""
        print *, "Eg.: ./make_data -matrix_type rmsd -help"
        print *, ""
        print *, "The type of matrices currently available are:"
        print *, "* rmsd"
        print *, "* z_coord"
        print *, "* atoms_dist"
        print *, "* angle"
        print *, ""
        STOP
      else if (trim(help_option) == "rmsd") then
        print *, ""
        print *, "To generate the rsmd matrix you will need to give pdb of solute 2, one group of ", &
        "atoms (either indexes or atomnames), the encounter complexes file, number of encounters, ", &
         "matrix type and matrix output file."
        print *, ""
        print *, "Eg: ./make_data -pdb2 p2_noh.pdb -atoms2 C N O CA -complexes assoc_complexes ", &
        "-matrix_type rmsd -matrix matrix_rmsd.txt"
        print *, ""
        STOP
      else if (trim(help_option) == "z_coord") then
        print *, ""
        print *, "To generate the z_coord matrix you will need to give the pdb of solute 2,", &
        "one group of atoms (either indexes or atomnames), the encounter complexes file, ", &
        "number of encounters, matrix type, array output file and matrix output file."
        print *, ""
        print *, "Eg.: ./make_data -pdb2 p2_noh.pdb -atoms2 Cu -complexes assoc_complexes ", &
        "-matrix_type z_coord -matrix matrix_z.txt -array array_z.txt"
        print *, ""
        STOP
      else if (trim(help_option) == "atoms_dist") then
        print *, ""
        print *, "To generate the atoms_dist matrix you will need to give two pdbs, two groups of atoms ", &
        "(either indexes or atomnames), the encounter complexes file, number of encounters, matrix type, ", &
        "array output file and matrix output file."
        print *, ""
        print *, "Eg.: ./make_data -pdb1 p1_noh.pdb -pdb2 p2_noh.pdb -complexes ", &
        "assoc_complexes -atoms1 C1 C4 -atoms2 Cu -matrix_type atoms_dist -matrix matrix_dist.txt ", &
        "-array array_dist.txt"
        print *, ""
        STOP
      else if (trim(help_option) == "2D_angle" .or. trim(help_option) == "3D_angle") then
        print *, ""
        print *, "To generate the angle matrix you will need to give two pdbs, four groups of atoms ", &
        "(either indexes or atomnames), the encounter complexes file, number of encounters, matrix ", & 
        "type, array output file and matrix output file."
        print *, ""
        print *, "Eg.: ./make_data -pdb1 p1_noh.pdb -pdb2 p2_noh.pdb -complexes ", &
        "assoc_complexes -atoms1a 1 -atoms1b 411 -atoms2a Cu -atoms2b 4 5 6 7 8 -matrix_type angle ", &
        "-matrix matrix_angle.txt -array array_angle.txt"
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
    !! @param[in]  atoms_array  Array of atom name strings or index strings
    !! @param[in]  pdb_filename  Reference filename for warnings
    !! @param[in]  tot_atoms     Total number of atoms in `pdb`
    !! @param[out] atoms_coord   Allocated array with atom coordinates (N,3)
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