!> Generate encounter selection based on threshold criteria.
!!
!! This executable computes thresholds (z_coord, atoms_dist, angles)
!! across encounter transforms and writes selected encounter records
!! to an output complexes file. See `print_help` for usage examples.
!!
!! @author Abraham Muñiz-Chicharro
program main
  USE read_input
  USE mod_threshold


  implicit none
      

    character*128 :: pdb1_filename, pdb2_filename, complexes_filename, output_complexes_filename, array_filename
    character*128 :: threshold_type
    character*128 :: argument
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
    integer :: nb_argument, count_arg, n_atoms
    integer :: ios
    logical :: pdb1_bool, pdb2_bool, complexes_bool, array_bool
    logical :: arr_atoms1a_bool, arr_atoms1b_bool, arr_atoms2a_bool, arr_atoms2b_bool
    logical :: cutoff_bool, thresholdtype_bool, output_complexes_bool
    logical :: nb_encounters_bool, help_bool

    nb_argument = 0
    nb_encounters = 0
    cutoff = 0.d0
    count_arg = 1
    pdb1_bool = .false.
    pdb2_bool = .false.
    complexes_bool = .false.
    cutoff_bool = .false.
    array_bool = .false.
    thresholdtype_bool = .false.
    output_complexes_bool = .false.
    nb_encounters_bool = .false.
    arr_atoms1a_bool = .false.
    arr_atoms1b_bool = .false.
    arr_atoms2a_bool = .false.
    arr_atoms2b_bool = .false.
    help_bool = .false.
    argument = ""

    
    nb_argument = command_argument_count()
    if ( nb_argument == 0 ) then
      print *, ""
      print *, "ERROR. Options were not given"
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
      else if ( trim(argument) == "-array" ) then
        array_bool = .true.
        call getarg( count_arg+1, argument )
        array_filename = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-complexes_output" ) then
        output_complexes_bool = .true.
        call getarg( count_arg+1, argument )
        output_complexes_filename = trim(argument)
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
      else if ( trim(argument) == "-cutoff" ) then
        cutoff_bool = .true.
        call getarg( count_arg+1, argument )
        read(argument, *, IOSTAT=ios) cutoff
        if (ios /= 0) then
          print *, "ERROR. Float expected for the -cutoff ", &
          "argument."
        end if
        count_arg = count_arg + 1
      else if ( trim(argument) == "-threshold_type" ) then
        thresholdtype_bool = .true.
        call getarg( count_arg+1, argument )
        threshold_type = trim(argument)
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
        print *, "./threshold -help"
        STOP 1
      end if

      count_arg = count_arg + 1

    end do


    if ( help_bool ) then
      if (thresholdtype_bool) then
        call print_help(threshold_type)
      else
        call print_help("main")
      end if
    end if

    if (.not. pdb2_bool ) then
      print *, "ERROR. PDB not provided."
      print *, "Please provide a pdb file with '-pdb2' option"
      print *, "For more information, please use the -help option:"
      print *, "./threshold -help"
      STOP 1
    else if ( .not. arr_atoms2a_bool ) then 
      print *, "ERROR. PDB2 provided but atoms2 not provided."
      print *, "Please provide a group of atoms with '-atoms2' option"
      print *, "For more information, please use the -help option:"
      print *, "./threshold -help"
      STOP 1
    else if (.not. array_bool) then
      print *, "ERROR. Array output filename not given."
      print *, "Please provide the array output file with the '-array option'"
      print *, "For more information, please use the -help option:"
      print *, "./threshold -help"
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

    allocate(distarray (nb_encounters))
    allocate(encounter_indexes (nb_encounters))

    if (trim(adjustl(threshold_type)) == "z_coord") then

      if (.not. output_complexes_bool) then
        print *, "ERROR. Complexes output filename not given."
        print *, "Please provide the complexes output file with the '-complexes_output option'"
        print *, "For more information, please use the -help option:"
        print *, "./threshold -threshold_type angle -help"
        STOP 1
      end if

      call read_atoms_coord(arr_atoms2a, atoms2a_coords, tot_atoms2, pdb2, pdb2_filename)

      call array_z_coord(distarray, nb_encounters, 1, complexes % xc1, complexes % xc2, &
      complexes % trans_vector, complexes % rot1, complexes % rot2, atoms2a_coords)

    else if (trim(adjustl(threshold_type)) == "3D_angle" .or. trim(adjustl(threshold_type)) == "2D_angle") then

      if (.not. pdb1_bool) then
        print *, "ERROR. Just one PDB provided. It is necessary to provide two pdbs ", &
        "to generate the angle threshold."
        print *, "Please provide the two pdb files with '-pdb1' and '-pdb2' options."
        print *, "For more information, please use the -help option:"
        print *, "./threshold -threshold_type angle -help"
        STOP 1
      else if (.not. arr_atoms1a_bool .or. .not. arr_atoms1b_bool .or. .not. arr_atoms2b_bool) then 
        print *, "ERROR. Not all necessary group of atoms given"
        print *, "Please provide all group of atoms with '-atoms1a', '-atoms1b', '-atoms2a' and '-atoms2b' options"
        print *, "For more information, please use the -help option:"
        print *, "./threshold -threshold_type angle -help"
        STOP 1
      else if (.not. output_complexes_bool) then
        print *, "ERROR. Complexes output filename not given."
        print *, "Please provide the complexes output file with the '-complexes_output option'"
        print *, "For more information, please use the -help option:"
        print *, "./threshold -threshold_type angle -help"
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
        call array_angle(distarray, nb_encounters, 2, &
        complexes % xc1, complexes % xc2, complexes % trans_vector, complexes % rot1, complexes % rot2, &
        cog1a, cog1b, cog2a, cog2b, 2)
      else if (trim(adjustl(threshold_type)) == "3D_angle") then
        call array_angle(distarray, nb_encounters, 2, &
        complexes % xc1, complexes % xc2, complexes % trans_vector, complexes % rot1, complexes % rot2, &
        cog1a, cog1b, cog2a, cog2b, 3)
      end if

    else if (trim(adjustl(threshold_type)) == "atoms_dist") then

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
        print *, "./threshold -threshold angle -help"
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

    encounter_indexes = 0
    do i = 1, nb_encounters
      encounter_indexes(i) = i
    end do

    call sort_array(distarray, encounter_indexes)

    call write_array(distarray, array_filename, cutoff)
    call write_complexes(complexes, distarray, encounter_indexes, nb_encounters, cutoff, output_complexes_filename)

    print *, 'Threshold complete and encounters stored in ', output_complexes_filename
  

  contains
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
        "threshold type and complexes output file."
        print *, "If you want instructions on how to create a specific type of threshold, type:"
        print *, "./threshold -threshold_type <type of threshold> -help"
        print *, ""
        print *, "Eg.: ./threshold -threshold_type rmsd -help"
        print *, ""
        print *, "The type of matrices currently available are:"
        print *, "* z_coord"
        print *, "* atoms_dist"
        print *, "* 2D_angle"
        print *, "* 3D_angle"
        print *, ""
        STOP
      else if (trim(help_option) == "z_coord") then
        print *, ""
        print *, "To generate the z_coord threshold you will need to give the pdb of solute 2,", &
        "one group of atoms (either indexes or atomnames), the encounter complexes file, ", &
        "number of encounters, threshold type and threshold output file."
        print *, ""
        print *, "Eg.: ./threshold -pdb2 p2_noh.pdb -atoms2 Cu -complexes assoc_complexes ", &
        "-threshold_type z_coord -complexes_output threshold_z.txt"
        print *, ""
        STOP
      else if (trim(help_option) == "atoms_dist") then
        print *, ""
        print *, "To generate the atoms_dist threshold you will need to give two pdbs, two groups of atoms ", &
        "(either indexes or atomnames), the encounter complexes file, number of encounters, threshold type ", &
        "and threshold output file."
        print *, ""
        print *, "Eg.: ./threshold -pdb1 p1_noh.pdb -pdb2 p2_noh.pdb -complexes ", &
        "assoc_complexes -atoms1 C1 C4 -atoms2 Cu -threshold_type atoms_dist -complexes_output threshold_dist.txt"
        print *, ""
        STOP
      else if (trim(help_option) == "2D_angle" .or. trim(help_option) == "3D_angle") then
        print *, ""
        print *, "To generate the angle threshold you will need to give two pdbs, four groups of atoms ", &
        "(either indexes or atomnames), the encounter complexes file, number of encounters, threshold ", & 
        "type and threshold output file."
        print *, ""
        print *, "Eg.: ./threshold -pdb1 p1_noh.pdb -pdb2 p2_noh.pdb -complexes ", &
        "assoc_complexes -atoms1a 1 -atoms1b 411 -atoms2a Cu -atoms2b 4 5 6 7 8 -threshold_type angle ", &
        "-complexes_output threshold_angle.txt"
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
    !! Similar to the helper in `make_matrix`. Returns an allocated
    !! `atoms_coord` array with matching coordinates for the atom names
    !! or indices provided in `atoms_array`.
    subroutine read_atoms_coord(atoms_array, atoms_coord, tot_atoms, pdb, pdb_filename)
      implicit none

      character(len=*), dimension(:), intent(in) :: atoms_array
      character(len=*), intent(in) :: pdb_filename
      integer, intent(in) :: tot_atoms
      real(kind=8), dimension(:, :), allocatable, intent(out) :: atoms_coord
      type ( type_pdb_file ), intent(in) :: pdb

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
    subroutine write_complexes(complexes_arg, distarray_arg, encounter_indexes_arg, nb_encounters_arg, cutoff_arg, output_name_arg)

      IMPLICIT NONE
      type(type_assoc_file), intent(in) :: complexes_arg
      real (kind=8), dimension(:), intent(in) :: distarray_arg
      integer, dimension(:), intent(in) :: encounter_indexes_arg
      integer, intent(in) :: nb_encounters_arg
      character*128, intent(in) :: output_name_arg
      real(kind=8), intent(in) :: cutoff_arg
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

    end subroutine write_complexes

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