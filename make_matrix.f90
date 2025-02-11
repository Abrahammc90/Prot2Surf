program main
    USE read_input
    USE mod_matrix


    implicit none
      

    character*128 :: pdb_filename, assoc_filename, surface_filename
    character*128 :: matrix_type, matrix_filename, array_filename
    character*128 :: str_n_enc

    integer :: tot_atoms, tot_encounters, tot_coords, tot_surface_atoms
    integer :: tot_chains, tot_surface_chains
    integer :: tot_chain_residues, tot_residue_atoms
    integer :: tot_residues, tot_surface_residues
    integer :: n_encounters
    type ( type_pdb_file ) :: pdb, surface_pdb
    type ( type_assoc_file ) :: assoc

    integer :: i, j, k
    real (kind=8), dimension(:, :), allocatable :: distmatrix
    real (kind=8), dimension(:), allocatable :: distarray
    real (kind=8), dimension(5, 3) :: HIX_ring_coords, HIY_ring_coords
    real (kind=8), dimension(:, :), allocatable :: protein_crds, surface_crds
    real (kind=8), dimension( 1, 3 ) :: Cu_coords
    real (kind=8), dimension( :, : ), allocatable :: bb_coords
    real (kind=8), dimension(3) :: HIS_ring_cog, residue_cog
    real (kind=8), dimension(:, :), allocatable :: residues_cog, residue_crds
    character(len=3), dimension(4) :: HIS_ring_atomnames
    character(len=1), dimension(3) :: bb_atomnames
    integer :: HIS_tot_coords, bb_index, position1, position2, position3
    integer :: crds_i
    integer :: ios
    integer :: last_index

    HIS_ring_atomnames = (/"CD2", "CE1", "ND1", "NE2"/)
    bb_atomnames = (/"C", "N", "O"/)

    pdb_filename = ''
    assoc_filename = ''
    matrix_type = ''
    str_n_enc = ''
    matrix_filename = ''

    call getarg ( 1, pdb_filename ) !Receives pdb input file
    call getarg ( 2, assoc_filename ) !Receives assoc input file
    call getarg ( 3, matrix_type ) !Receives cluster mode
    call getarg ( 4, str_n_enc ) !Receives number of encounters
    call getarg ( 5, matrix_filename ) !Receives name of the matrix output file
    call getarg ( 6, array_filename ) !Receives name of the array output file
    call getarg ( 7, surface_filename ) !Receives name of the array output file

    read(str_n_enc, *, iostat=ios) n_encounters


    !Reads protein PDB
    
    call read_pdb(pdb, pdb_filename, tot_atoms, tot_residues, tot_chains)
    
    !Reads surface PDB
    call read_pdb(surface_pdb, surface_filename, tot_surface_atoms, tot_surface_residues, tot_surface_chains)

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

    allocate(distmatrix (n_encounters, n_encounters))
    allocate(distarray (n_encounters))

    print *, 'Building ', trim(matrix_type), ' matrix'

    if (trim(adjustl(matrix_type)) == "Cu_z") then
      HIS_tot_coords = 0
      tot_coords = 0

      do i = 1, tot_atoms
        if ( trim(adjustl(pdb % atoms(i) % name)) == 'Cu' ) then
          tot_coords = tot_coords + 1
          Cu_coords(1, :) = pdb % atoms(i) % coord(:)
        end if
      end do

      call matrix_z(distmatrix, distarray, n_encounters, 1, assoc % xc1, assoc % xc2, &
      assoc % trans_vector, assoc % rot1, assoc % rot2, Cu_coords)

      call write_array(distarray, array_filename)

    else if (trim(adjustl(matrix_type)) == "rmsd") then
      tot_coords = 0

      do i = 1, tot_atoms
        !write(*,*) pdb % atomname(i)(1:3)
        !if ( ANY( trim(adjustl(pdb % atomname(i))) == &
        !  bb_atomnames ) .or. trim(adjustl(pdb % atomname(i))) == "CA" ) then
        position1 = INDEX(pdb % atoms(i) % name(:2), "C")
        position2 = INDEX(pdb % atoms(i) % name(:2), "N")
        position3 = INDEX(pdb % atoms(i) % name(:2), "O")
        if ( position1 > 0 .or. position2 > 0 .or. position3 > 0 ) then
          tot_coords = tot_coords + 1
        end if
      end do

      allocate(bb_coords(tot_coords, 3))
      bb_index = 0
      do i = 1, tot_atoms
        !if ( ANY( trim(adjustl(pdb % atomname(i))) == &
        !  bb_atomnames ) .or. trim(adjustl(pdb % atomname(i))) == "CA" ) then
        position1 = INDEX(pdb % atoms(i) % name(:2), "C")
        position2 = INDEX(pdb % atoms(i) % name(:2), "N")
        position3 = INDEX(pdb % atoms(i) % name(:2), "O")
        if ( position1 > 0 .or. position2 > 0 .or. position3 > 0 ) then
            bb_index = bb_index + 1
          bb_coords(bb_index, :) = pdb % atoms(i) % coord(:)
        end if
      end do

      !write(*,*) bb_coords(1, :)

      call matrix_rmsd(distmatrix, n_encounters, tot_coords, &
      assoc % xc1, assoc % xc2, assoc % trans_vector, &
      assoc % rot1, assoc % rot2, bb_coords)

      

    else if (trim(adjustl(matrix_type)) == "surface_chain_angle") then
      HIS_tot_coords = 0
      tot_coords = 0

      do i = 1, tot_atoms

        if ( trim(adjustl(pdb % atoms(i) % name)) == 'Cu' ) then
          tot_coords = tot_coords + 1
          Cu_coords(1, :) = pdb % atoms(i) % coord(:)
        else if ( pdb % atoms(i) % resname == "HIX" ) then

          if (ANY( trim(adjustl(pdb % atoms(i) % name)) == &
          HIS_ring_atomnames ) .or. &
          trim(adjustl(pdb % atoms(i) % name)) == "CG") then
            tot_coords = tot_coords + 1
            HIS_tot_coords = HIS_tot_coords + 1
            HIX_ring_coords(HIS_tot_coords, :) = pdb % atoms(i) % coord(:)
          end if
        end if
      end do

      call calculate_cog(HIS_ring_cog, HIX_ring_coords, HIS_tot_coords)

      allocate(protein_crds(2, 3))

      protein_crds(1, :) = Cu_coords(1, :)
      protein_crds(2, :) = HIS_ring_cog

      !v1(1) = Cu_coords(1, 1) - ring_cog(1)
      !v1(2) = Cu_coords(1, 2) - ring_cog(2)
      !v1(3) = Cu_coords(1, 3) - ring_cog(3)

      allocate(surface_crds(2, 3))

      last_index = size(surface_pdb % chains(1) % atoms)

      surface_crds(1, :) = surface_pdb % chains(1) % atoms(1) % coord
      surface_crds(2, :) = surface_pdb % chains(1) % atoms(last_index) % coord

      !v2(1) = surface_point1(1) - surface_point2(1)
      !v2(2) = surface_point1(2) - surface_point2(2)
      !v2(3) = surface_point1(3) - surface_point2(3)
      
      call matrix_chain_degree(distmatrix, distarray, n_encounters, 2, &
      assoc % xc1, assoc % xc2, assoc % trans_vector, assoc % rot1, assoc % rot2, &
      protein_crds, surface_crds)

      call write_array(distarray, array_filename)

    else if (trim(adjustl(matrix_type)) == "C1-C4_dist") then
      do i = 1, tot_atoms
        if ( trim(adjustl(pdb % atoms(i) % name)) == 'Cu' ) then
          Cu_coords(1, :) = pdb % atoms(i) % coord(:)
        end if
      end do

      tot_coords = 0
      do i = 1, tot_surface_atoms
        if ( trim(adjustl(surface_pdb % atoms(i) % name)) == 'C1' .or. &
        trim(adjustl(surface_pdb % atoms(i) % name)) == 'C4' ) then
          tot_coords = tot_coords + 1
        end if
      end do

      allocate(protein_crds(tot_coords, 3))

      crds_i = 0
      do i = 1, tot_surface_atoms
        if ( trim(adjustl(surface_pdb % atoms(i) % name)) == 'C1' .or. &
        trim(adjustl(surface_pdb % atoms(i) % name)) == 'C4' ) then
          crds_i = crds_i + 1
          protein_crds(crds_i, :) = surface_pdb % atoms(i) % coord(:)
        end if
      end do

      call matrix_dist(distmatrix, distarray, n_encounters, 1, assoc % xc1, assoc % xc2, &
      assoc % trans_vector, assoc % rot1, assoc % rot2, Cu_coords, protein_crds)

      call write_array(distarray, array_filename)

    else
      print *, 'ERROR: ', trim(adjustl(matrix_type)), ' matrix not implemented'
      STOP
    end if

    call write_matrix(distmatrix, matrix_filename)

    print *, 'Matrix complete and stored in ', matrix_filename

END PROGRAM main