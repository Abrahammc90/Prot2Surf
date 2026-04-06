!> \file clust.f90
!! \brief Hierarchical clustering driver for encounter matrices/arrays.
!!
!! clust.f90 - Hierarchical clustering for encounter matrices.
!!
!! Usage: ./clust -help
!!
!! Reads a pairwise matrix (or array) and runs hierarchical clustering with user-specified linkage. Results written to output_name files.
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0
!! @date 2026-04-05
program main

  USE mod_matrix
  USE mod_array
  USE mod_clust_algorithm
  USE read_input

  implicit none

  ! High-level workflow:
  ! 1) Parse CLI options and validate required inputs.
  ! 2) Read encounter complexes metadata (transforms + headers).
  ! 3) Auto-detect whether -input is a 1D array or 2D matrix.
  ! 4) Run hierarchical clustering with requested linkage.
  ! 5) Write clustering outputs with basename given in -output_name.

    character*128 :: datadist_filename, complexes_filename
    character*128 :: output_name, linkage_type
    integer :: nb_encounters, tot_encounters
    type ( type_assoc_file ) :: complexes

    real (kind=8), dimension(:, :), allocatable :: matrix
    real (kind=8), dimension(:), allocatable :: array
    logical :: complexes_bool, nb_encounters_bool, help_bool
    logical :: matrix_bool, array_bool, output_name_bool
    integer :: detect_unit, detect_stat, nlines
    character(len=1) :: dummy

    call parse_arguments( &
      complexes_filename, complexes_bool, &
      nb_encounters, nb_encounters_bool, &
      datadist_filename, &
      output_name, output_name_bool, &
      linkage_type, help_bool)

    ! Help path exits before any heavy I/O or clustering work.
    if ( help_bool ) then
      call print_help()
      STOP 0
    end if

    ! Validate required options shared by both array/matrix input modes.
    if ( .not. complexes_bool ) then
      print *, '[clust] ERROR: Encounter complexes file not provided. Use -complexes <file>'
      print *, '[clust] For usage, run: ./clust -help'
      STOP 1
    end if

    if ( len_trim(datadist_filename) == 0 ) then
      print *, '[clust] ERROR: No input file provided. Use -input <file>'
      print *, '[clust] For usage, run: ./clust -help'
      STOP 1
    end if

    if ( .not. output_name_bool ) then
      print *, '[clust] ERROR: Output name not provided. Use -output_name <name>'
      print *, '[clust] For usage, run: ./clust -help'
      STOP 1
    end if

    ! Read encounter-complex records and transformation parameters.
    call read_assoc(complexes, complexes_filename, tot_encounters)

    ! Keep number of encounters inside available range.
    if ( nb_encounters .gt. tot_encounters ) then
      print *, '[clust] WARNING: -nb_encounters > available. Using all encounters in ', trim(complexes_filename)
      nb_encounters = tot_encounters
      print *, '[clust] Total encounters available: ', nb_encounters
    else if (.not. nb_encounters_bool) then
      print *, '[clust] WARNING: -nb_encounters not given. Using all available.'
      nb_encounters = tot_encounters
      print *, '[clust] If you want a specific number, use -nb_encounters <N>'
    end if

    ! Auto-detect input format:
    ! - single line -> 1D array path (array clustering)
    ! - multiple lines -> 2D matrix path (matrix clustering)
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

    ! Dispatch to the proper clustering routine based on detected input type.
    if ( array_bool ) then
      call read_array(array, nb_encounters, datadist_filename)
      call linkage_clustering_from_array(array, nb_encounters, linkage_type, output_name, complexes)
    else
      call read_matrix(matrix, nb_encounters, datadist_filename)
      call linkage_clustering_from_matrix(matrix, nb_encounters, linkage_type, output_name, complexes)
    end if

    contains
    !> Parse command-line arguments for `clust`.
    !!
    !! Parses all CLI flags and fills program configuration values.
    !! @param[out] p_complexes_filename Value parsed from `-complexes`
    !! @param[out] p_complexes_bool     Whether `-complexes` was provided
    !! @param[out] p_nb_encounters      Value parsed from `-nb_encounters`
    !! @param[out] p_nb_encounters_bool Whether `-nb_encounters` was provided
    !! @param[out] p_datadist_filename  Value parsed from `-input`
    !! @param[out] p_output_name        Value parsed from `-output_name`
    !! @param[out] p_output_name_bool   Whether `-output_name` was provided
    !! @param[out] p_linkage_type       Value parsed from `-linkage` (default `mean`)
    !! @param[out] p_help_bool          Whether `-help` was provided
    subroutine parse_arguments( &
      p_complexes_filename, p_complexes_bool, &
      p_nb_encounters, p_nb_encounters_bool, &
      p_datadist_filename, &
      p_output_name, p_output_name_bool, &
      p_linkage_type, p_help_bool)
      implicit none

      character*128, intent(out) :: p_complexes_filename
      logical, intent(out) :: p_complexes_bool
      integer, intent(out) :: p_nb_encounters
      logical, intent(out) :: p_nb_encounters_bool
      character*128, intent(out) :: p_datadist_filename
      character*128, intent(out) :: p_output_name
      logical, intent(out) :: p_output_name_bool
      character*128, intent(out) :: p_linkage_type
      logical, intent(out) :: p_help_bool

      integer :: nb_argument, count_arg, ios
      character*128 :: argument

      nb_argument = command_argument_count()
      if ( nb_argument == 0 ) then
        print *, '[clust] ERROR: No arguments parsed.'
        print *, '[clust] For usage, run: ./clust -help'
        call print_help()
        STOP 1
      end if

      p_nb_encounters = 0
      p_complexes_bool = .false.
      p_output_name_bool = .false.
      p_nb_encounters_bool = .false.
      p_help_bool = .false.
      p_datadist_filename = ''
      p_linkage_type = 'mean'
      argument = ""
      count_arg = 1

      do while ( count_arg <= nb_argument )
        call getarg( count_arg, argument )

        if ( trim(argument) == "-complexes" ) then
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
        else if ( trim(argument) == "-input" ) then
          call getarg( count_arg+1, argument )
          p_datadist_filename = trim(argument)
          count_arg = count_arg + 1
        else if ( trim(argument) == "-output_name" ) then
          p_output_name_bool = .true.
          call getarg( count_arg+1, argument )
          p_output_name = trim(argument)
          count_arg = count_arg + 1
        else if ( trim(argument) == "-linkage" ) then
          call getarg( count_arg+1, argument )
          p_linkage_type = trim(argument)
          if (trim(p_linkage_type) /= 'min' .and. trim(p_linkage_type) /= 'max' .and. &
              trim(p_linkage_type) /= 'mean') then
            print *, '[clust] ERROR: Invalid linkage type: ', trim(p_linkage_type)
            print *, '[clust] Valid options: min, max, mean'
            STOP 1
          end if
          count_arg = count_arg + 1
        else if ( trim(argument) == "-help" ) then
          p_help_bool = .true.
        else
          print *, '[clust] ERROR: Argument ''', trim(argument), ''' not recognized.'
          print *, '[clust] For usage, run: ./clust -help'
          STOP 1
        end if

        count_arg = count_arg + 1
      end do

    end subroutine parse_arguments

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