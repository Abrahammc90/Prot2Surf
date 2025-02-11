program main

    USE mod_matrix
    USE mod_clust_algorithm
    USE read_input

    implicit none

    character*128 :: matrix_filename, str_n_enc, str_ncls, complexes_file, array_filename
    character*128 :: output_name
    integer :: ncls, n_encounters, tot_encounters
    type ( type_assoc_file ) :: assoc

    real (kind=8), dimension(:, :), allocatable :: matrix
    real (kind=8), dimension(:), allocatable :: array
    integer :: ios

    array_filename = ''

    call getarg ( 1, matrix_filename ) !Receives filename where matrix is store
    call getarg ( 2, str_n_enc ) !Receives number of encounters
    call getarg ( 3, str_ncls ) !Receives number of clusters to generate
    call getarg ( 4, complexes_file ) !Receives filename where encounters information is store
    call getarg ( 5, output_name ) !Receives filename of output file where cluster encounters will be recorded
    call getarg ( 6, array_filename ) !Receives filename of input file where values of enncounter complexes are recorded
    
    read(str_n_enc, *, iostat=ios) n_encounters
    read(str_ncls, *, iostat=ios) ncls

    !Reads Complexes
    call read_assoc(assoc, complexes_file, tot_encounters)

    if ( n_encounters .gt. tot_encounters ) then
      write(*,*) 'WARNING: the number of encounters given is greater than'
      write(*,*) 'the number of encounters in ', complexes_file
      write(*,*) 'Setting the number of encounters to the total'
      write(*,*) 'number of encounters available'
      n_encounters = tot_encounters
      write(*,*) 'Total encounters available: ', n_encounters
      write(*,*) ''
    end if

    if (ncls .ge. n_encounters) then
      write(*,*) "WARNING: number of clusters is equal or greater than number of encounter."
      write(*,*) "Setting the number of clusters to the half of number of encounters"
      ncls = n_encounters / 2
    end if

    call read_matrix(matrix, n_encounters, matrix_filename)

    if ( len(array_filename) .gt. 0 ) then
      call read_array(array, n_encounters, array_filename)
      call mean_linkage_clustering(matrix, n_encounters, ncls, output_name, assoc, array)
    else
      call mean_linkage_clustering(matrix, n_encounters, ncls, output_name, assoc)
    end if
  
  end program main