!> Command-line program to perform hierarchical clustering on encounter matrices.
!!
!! Reads a pairwise matrix (and optional array) and runs hierarchical
!! clustering with user-specified linkage method (min, max, or mean).
!! Results are written to files using the `-output_name` base name.
!! See `print_help` below for usage.
!!
!! @author Abraham Muñiz-Chicharro
program main

  USE mod_matrix
  USE mod_array
  USE mod_clust_algorithm
  USE read_input

  implicit none

    character*128 :: datadist_filename, complexes_filename
    character*128 :: output_name, argument, linkage_type
    integer :: nb_encounters, tot_encounters
    type ( type_assoc_file ) :: complexes

    real (kind=8), dimension(:, :), allocatable :: matrix
    real (kind=8), dimension(:), allocatable :: array
    integer :: ios, nb_argument, count_arg
    logical :: complexes_bool, nb_encounters_bool, help_bool
    logical :: matrix_bool, array_bool, output_name_bool
    logical :: use_cuda_bool
    integer :: detect_unit, detect_stat, nlines
    character(len=1) :: dummy

    nb_argument = 0
    nb_encounters = 0
    count_arg = 1
    complexes_bool = .false.
    matrix_bool = .false.
    array_bool = .false.
    output_name_bool = .false.
    use_cuda_bool = .false.
    output_name_bool = .false.
    nb_encounters_bool = .false.
    help_bool = .false.
    argument = ""
    datadist_filename = ''
    linkage_type = 'mean'  ! Default to mean linkage

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
      else if ( trim(argument) == "-input" ) then
        call getarg( count_arg+1, argument )
        datadist_filename = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-output_name" ) then
        output_name_bool = .true.
        call getarg( count_arg+1, argument )
        output_name = trim(argument)
        count_arg = count_arg + 1
      else if ( trim(argument) == "-linkage" ) then
        call getarg( count_arg+1, argument )
        linkage_type = trim(argument)
        if (trim(linkage_type) /= 'min' .and. trim(linkage_type) /= 'max' .and. &
            trim(linkage_type) /= 'mean') then
          print *, "ERROR. Invalid linkage type: ", trim(linkage_type)
          print *, "Valid options are: 'min', 'max', 'mean'"
          STOP 1
        end if
        count_arg = count_arg + 1
      else if ( trim(argument) == "-help" ) then
        help_bool = .true.
      else if ( trim(argument) == "-cuda" ) then
        use_cuda_bool = .true.
      else
        print *,  "ERROR. Argument ", trim(argument), " not recognized."
        print *, "For more information, please use the -help option:"
        print *, "./clust -help"
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
      print *, "./clust -help"
      STOP 1
    end if

    if ( len_trim(datadist_filename) == 0 ) then
      print *, "ERROR. No input file provided."
      print *, "Please provide an input file (matrix or array) with '-input'"
      print *, "For more information, please use the -help option:"
      print *, "./clust -help"
      STOP 1
    end if

    if ( .not. output_name_bool ) then
      print *, "ERROR. Output name not provided."
      print *, "Please, provide the output_name with the '-output_name' option"
      print *, "For more information, please use the -help option:"
      print *, "./clust -help"
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

    ! Auto-detect input format: 1 line = array, >1 lines = matrix
    detect_unit = 99
    open(detect_unit, file=trim(datadist_filename), status='old', iostat=detect_stat)
    if (detect_stat /= 0) then
      print *, "ERROR. Cannot open input file: ", trim(datadist_filename)
      STOP 1
    end if
    nlines = 0
    do
      read(detect_unit, '(A)', iostat=detect_stat) dummy
      if (detect_stat /= 0) exit
      nlines = nlines + 1
    end do
    close(detect_unit)

    if (nlines == 0) then
      print *, "ERROR. Input file is empty: ", trim(datadist_filename)
      STOP 1
    else if (nlines == 1) then
      array_bool = .true.
      write(*,*) 'Detected input as 1D array (single line)'
    else
      matrix_bool = .true.
      write(*,*) 'Detected input as 2D matrix (', nlines, ' lines)'
    end if

    if ( array_bool ) then
      call read_array(array, nb_encounters, datadist_filename)
      call linkage_clustering_from_array(array, nb_encounters, linkage_type, output_name, complexes)
    else
      call read_matrix(matrix, nb_encounters, datadist_filename)
      call linkage_clustering_from_matrix(matrix, nb_encounters, linkage_type, output_name, complexes, use_cuda=use_cuda_bool)
    end if

    contains
    !> Print usage and help information for the `clust` program.
    !!
    !! Describes command-line options and provides an example invocation.
    subroutine print_help()

      print *, ""
      print *, "This program performs hierarchical clustering on encounter complexes."
      print *, "It receives as inputs the encounter complexes file, matrix input file,"
      print *, "optional array input file, and clustering parameters."
      print *, ""
      print *, "Usage:"
      print *, "  ./clust -complexes <file> -input <file> [OPTIONS]"
      print *, ""
      print *, "Required arguments:"
      print *, "  -complexes <file>       Encounter complexes file"
      print *, "  -input <file>           Input file (matrix or array, auto-detected)"
      print *, "  -output_name <name>     Base name for output files"
      print *, ""
      print *, "Optional arguments:"
      print *, "  -nb_encounters <N>      Number of encounters to cluster"
      print *, "                          (default: all encounters in complexes file)"
      print *, "  -linkage <type>         Linkage method: 'min', 'max', or 'mean'"
      print *, "                          (default: 'mean')"
      print *, "  -cuda                   Use hybrid CPU/GPU acceleration"
      print *, "  -help                   Display this help message"
      print *, ""
      print *, "Linkage types:"
      print *, "  min  - Minimum linkage (single linkage)"
      print *, "  max  - Maximum linkage (complete linkage)"
      print *, "  mean - Mean linkage (average linkage)"
      print *, ""
      print *, "The input format is auto-detected:"
      print *, "  - 1 line  -> treated as a 1D array"
      print *, "  - N lines -> treated as an NxN distance matrix"
      print *, ""
      print *, "Example:"
      print *, "  ./clust -complexes assoc_complexes -input matrix_z.txt \\"
      print *, "          -nb_encounters 5000 -output_name Cu_z -linkage mean"
      print *, ""
      STOP
      
    end subroutine print_help
  
  end program main