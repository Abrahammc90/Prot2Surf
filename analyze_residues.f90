program main
    USE read_input
    USE maths


    implicit none
      

    character*128 :: pdb_filename, assoc_filename!, surface_filename
    !character*128 :: matrix_type, matrix_filename, array_filename
    character*128 :: str_n_enc

    integer :: tot_atoms, tot_encounters!, tot_surface_atoms
    integer :: tot_chains!, tot_surface_chains
    integer :: tot_residue_atoms
    integer :: tot_residues!, tot_surface_residues
    integer :: n_encounters
    type ( type_pdb_file ) :: pdb!, surface_pdb
    type ( type_assoc_file ) :: assoc

    integer :: i, j, k, min_j
    !real (kind=8), dimension(:, :), allocatable :: distmatrix
    !real (kind=8), dimension(:), allocatable :: distarray
    real (kind=8), dimension(3) :: residue_cog
    real (kind=8), dimension(:, :), allocatable :: residues_cog, residue_crds, new_residues_cog
    !real (kind=8), dimension(:), allocatable :: new_residues_cog_z
    real (kind=8) :: min_z
    !integer, dimension(:), allocatable :: encounters_closest_resid, encounters_count, counter
    integer, dimension(:), allocatable :: encounters_count
    integer, dimension(:, :), allocatable :: all_closest_resid
    integer :: closest_resid
    character (len=12) :: residue_str
    character (len=4) :: str_resid
    integer :: ios, unit_number

    pdb_filename = ''
    assoc_filename = ''
    !matrix_type = ''
    str_n_enc = ''
    !matrix_filename = ''

    call getarg ( 1, pdb_filename ) !Receives pdb input file
    call getarg ( 2, assoc_filename ) !Receives assoc input file
    !call getarg ( 3, matrix_type ) !Receives cluster mode
    call getarg ( 3, str_n_enc ) !Receives number of encounters
    !call getarg ( 5, matrix_filename ) !Receives name of the matrix output file
    !call getarg ( 6, array_filename ) !Receives name of the array output file
    !call getarg ( 7, surface_filename ) !Receives name of the array output file

    read(str_n_enc, *, iostat=ios) n_encounters


    !Reads protein PDB
    
    call read_pdb(pdb, pdb_filename, tot_atoms, tot_residues, tot_chains)
    
    !Reads surface PDB
    !call read_pdb(surface_pdb, surface_filename, tot_surface_atoms, tot_surface_residues, tot_surface_chains)

    !Reads Complexes
    call read_assoc(assoc, assoc_filename, tot_encounters)

    if ( n_encounters .gt. tot_encounters ) then
        write(*,*) 'WARNING: the number of encounters given is greater than'
        write(*,*) 'the number of encounters in ', assoc_filename
        write(*,*) 'Setting the number of encounters to the total'
        write(*,*) 'number of encounters available'
        n_encounters = tot_encounters
        write(*,*) 'Total encounters available: ', n_encounters
        write(*,*) ''
    end if

    !tot_encounters = 500

    !allocate(distmatrix (n_encounters, n_encounters))
    !allocate(distarray (n_encounters))


    !Calculate center of mass of all residues
    allocate(residues_cog(tot_residues, 3))

    !write(*,*) size(pdb % residues)
    !write(*,*) tot_residues
    !STOP 1
    !do j = 1, tot_residues
    !    write(*,*) pdb % residues(j) % natoms, j
    !end do
    !STOP 1

    do j = 1, tot_residues
        tot_residue_atoms = pdb % residues(j) % natoms
        allocate(residue_crds(tot_residue_atoms, 3))
        do k = 1, tot_residue_atoms
          residue_crds(k, :) = pdb % residues(j) % atoms(k) % coord
        end do
        call calculate_cog(residue_cog, residue_crds, tot_residue_atoms)
        residues_cog(j, :) = residue_cog
        deallocate(residue_crds)
    end do
    !STOP 1

    !Rotate and translate center of mass of all residues
    allocate(new_residues_cog(tot_residues, 3))
    !allocate(new_residues_cog_z(tot_residues))
    !allocate(encounters_closest_resid(tot_encounters))
    allocate(encounters_count(tot_residues))
    !allocate(counter(tot_residues))
    allocate(all_closest_resid(tot_residues, tot_encounters))

    encounters_count = 0
    all_closest_resid = 0

    

    do i = 1, tot_encounters

        min_j = 999999
        min_z = 999999.9

        call update_complex(assoc % xc1, assoc % xc2, &
            assoc % trans_vector (i, :), assoc % rot1 (i, :), assoc % rot2 (i, :), &
            tot_residues, residues_cog, new_residues_cog)

        do j = 1, tot_residues
            
            !if ( new_residues_cog(j, 3) .lt. min_z ) then
            !    min_z = new_residues_cog(j, 3)
            !    min_j = j
            !end if
            if ( new_residues_cog(j, 3) .lt. 6 ) then
                encounters_count(j) = encounters_count(j) + 1
                all_closest_resid(j, encounters_count(j)) = i
            end if
            !tot_residue_atoms = pdb % residues(j) % natoms
            !new_residues_cog_z(j) = new_residues_cog(j, 3)

        end do

        !encounters_count(min_j) = encounters_count(min_j) + 1
        !all_closest_resid(min_j, encounters_count(min_j)) = i

        
    end do

    !counter = 0
    !Write information

    ! Open file to write the final cluster content
    unit_number = 20
    open(unit=unit_number, file='closest_residues.txt', status='replace', action='write')
    do j = 1, tot_residues

        closest_resid = pdb % residues(j) % resid
        !write(*,*) closest_resid
        !allocate(closest_resid(tot_residues, encounters_count(i)))
        write(str_resid, '(I4)') closest_resid
        residue_str = "Residue" // adjustr(str_resid) // ":"
        !do j = 1, tot_encounters
        !    closest_resid = encounters_closest_resid(j)
        !    counter(closest_resid) = counter(closest_resid) + 1
        !    all_closest_resid(closest_resid, counter(closest_resid)) = j
        !end do
        if (encounters_count(j) .gt. 0) then
            write(unit_number, "(1(A12)10000(I5,1X))") residue_str, (all_closest_resid(j, k), k = 1, encounters_count(j))
        end if

    end do



end program main