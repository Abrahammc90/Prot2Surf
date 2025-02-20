program main

    USE mod_matrix
    USE mod_clust_algorithm
    USE read_input

    implicit none

    character*128 :: matrix_filename, str_nb_enc, complexes_filename, array_filename
    character*128 :: output_name, argument
    integer :: nb_encounters, tot_encounters
    real(kind=8) :: k_factor
    type ( type_assoc_file ) :: complexes

    real (kind=8), dimension(:, :), allocatable :: matrix
    real (kind=8), dimension(:), allocatable :: array
    integer :: ios, nb_argument, count_arg
    logical :: complexes_bool, nb_encounters_bool, help_bool
    logical :: matrix_bool, array_bool, output_name_bool, k_factor_bool

    nb_argument = 0
    nb_encounters = 0
    count_arg = 1
    complexes_bool = .false.
    matrix_bool = .false.
    array_bool = .false.
    output_name_bool = .false.
    nb_encounters_bool = .false.
    k_factor_bool = .false.
    help_bool = .false.
    argument = ""
    array_filename = ''
    k_factor = 1.0

    nb_argument = command_argument_count()
    if ( nb_argument == 0 ) then
      print *, ""
      print *, "ERROR. No arguments were parsed"
      call print_help()
      STOP 1
    end if

    do while ( count_arg <= nb_argument )
      call getarg( count_arg, argument )
      
      if ( trim(argument) == "-complexes" ) then
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
      else if ( trim(argument) == "-matrix" ) then
        matrix_bool = .true.
        call getarg( count_arg+1, argument )
        matrix_filename = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-array" ) then
        array_bool = .true.
        call getarg( count_arg+1, argument )
        array_filename = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-output_name" ) then
        output_name_bool = .true.
        call getarg( count_arg+1, argument )
        output_name = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-k_factor" ) then
        k_factor_bool = .true.
        call getarg( count_arg+1, argument )
        read(argument, *, IOSTAT=ios) k_factor
        if (ios /= 0) then
          print *, ""
          print *, "ERROR. Float expected for the -k_factor ", &
          "argument."
        end if
        count_arg = count_arg + 1
      else if ( trim(argument) == "-help" ) then
        help_bool = .true.
      else
        print *,  "ERROR. Argument ", trim(argument), " not recognized."
        print *, "For more information, please use the -help option:"
        print *, "./make_matrix -help"
        STOP 1
      end if

      count_arg = count_arg + 1

    end do

    if ( help_bool ) then
      call print_help()
    end if

    if ( .not. complexes_bool ) then
      print *, "ERROR. Encounter complexes file not provided."
      print *, "Please provide the complexes file with the '-complexes' option"
      print *, "For more information, please use the -help option:"
      print *, "./make_matrix -help"
      STOP 1
    else if ( .not. array_bool ) then
      print *, "WARNING. Array input file not provided."
      print *, "Be aware, this file is necessary to correctly calculate the average", &
      "and standard deviation in case that the cluster analysis is not based on rmsd", &
      "between the different encounters."
      print *, "If you need to provide the array input file, please use the '-array' option"
      print *, "For more information, please use the -help option:"
      print *, "./make_matrix -help"
      STOP 1
    else if ( .not. matrix_bool ) then
      print *, "ERROR. Matrix input file not provided."
      print *, "Please, provide the matrix input file with the '-matrix' option"
      print *, "For more information, please use the -help option:"
      print *, "./make_matrix -help"
      STOP 1
    else if ( .not. output_name_bool ) then
      print *, "ERROR. Output name not provided."
      print *, "Please, provide the output_name with the '-output_name' option"
      print *, "For more information, please use the -help option:"
      print *, "./make_matrix -help"
      STOP 1
    end if

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

    !if (ncls .ge. n_encounters) then
    !  write(*,*) "WARNING: number of clusters is equal or greater than number of encounter."
    !  write(*,*) "Setting the number of clusters to the half of number of encounters"
    !  ncls = n_encounters / 2
    !end if

    call read_matrix(matrix, nb_encounters, matrix_filename)

    if ( len(array_filename) .gt. 0 ) then
      call read_array(array, nb_encounters, array_filename)
      call mean_linkage_clustering(matrix, nb_encounters, k_factor, output_name, complexes, array)
    else
      call mean_linkage_clustering(matrix, nb_encounters, k_factor, output_name, complexes)
    end if

    contains
    subroutine print_help()

      !STOP 1

      print *, ""
      print *, "This program receives as inputs the encounter complexes file, ", &
      "number of encounters to cluster, matrix input file, array input file and ", &
      "a k_factor number."
      print *, ""
      print *, "Eg.: ./clust -complexes assoc_complexes -matrix matrix_z.txt -array array_z.txt ", &
      "-nb_encounters 5000 -k_factor 1.0 -output_name Cu_z"
      print *, ""
      print *, "Default values: "
      print *, "* k_factor: 1.0. "
      print *, "* nb_encounters: Maximum encounter complexes recorded in the complexes file."
      print *, ""
      STOP
      
    end subroutine print_help
  
  end program main