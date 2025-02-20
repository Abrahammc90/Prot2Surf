program main
    USE read_input
    USE maths


    implicit none
      

    character*128 :: pdb2_filename, complexes_filename

    integer :: tot_atoms2, tot_residues2, tot_chains2, tot_residue_atoms2
    integer :: tot_encounters, nb_encounters
    type ( type_pdb_file ) :: pdb1, pdb2
    type ( type_assoc_file ) :: complexes

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

    pdb2_bool = .false.
    complexes_bool = .false.
    nb_encounters_bool = .false.
    dist_threshold = 6.0
    nb_argument = 0
    count_arg = 1
    nb_argument = command_argument_count()

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


    !Reads protein PDB
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

    allocate(residues_cog(tot_residues2, 3))

    do j = 1, tot_residues2
        tot_residue_atoms2 = pdb2 % residues(j) % natoms
        allocate(residue_crds(tot_residue_atoms2, 3))
        do k = 1, tot_residue_atoms2
          residue_crds(k, :) = pdb2 % residues(j) % atoms(k) % coord
        end do
        call calculate_cog(residue_cog, residue_crds, tot_residue_atoms2)
        residues_cog(j, :) = residue_cog
        deallocate(residue_crds)
    end do
    
    allocate(new_residues_cog(tot_residues2, 3))
    allocate(encounters_count(tot_residues2))
    allocate(all_closest_resid(tot_residues2, tot_encounters))

    encounters_count = 0
    all_closest_resid = 0

    

    do i = 1, tot_encounters

        call update_complex(complexes % xc1, complexes % xc2, &
            complexes % trans_vector (i, :), complexes % rot1 (i, :), complexes % rot2 (i, :), &
            tot_residues2, residues_cog, new_residues_cog)

        do j = 1, tot_residues2
            
            if ( new_residues_cog(j, 3) .lt. dist_threshold ) then
                encounters_count(j) = encounters_count(j) + 1
                all_closest_resid(j, encounters_count(j)) = i
            end if

        end do

    end do

    ! Open file to write the final cluster content
    unit_number = 20
    open(unit=unit_number, file='closest_residues.txt', status='replace', action='write')
    do j = 1, tot_residues2

        closest_resid = pdb2 % residues(j) % resid
        write(str_resid, '(I4)') closest_resid
        residue_str = "Residue" // adjustr(str_resid) // ":"
        if (encounters_count(j) .gt. 0) then
            write(unit_number, "(1(A12)10000(I5,1X))") residue_str, (all_closest_resid(j, k), k = 1, encounters_count(j))
        end if

    end do

    contains
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
      print *, "* nb_encounters: Maximum encounter complexes recorded in the complexes file."
      print *, "* thresold: 6.0"
      print *, ""
      STOP
      
    end subroutine print_help
  

end program main