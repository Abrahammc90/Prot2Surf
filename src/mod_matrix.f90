MODULE mod_matrix

  !> Routines to build and read/write encounter-based matrices.
  !!
  !! This module contains helpers to construct pairwise matrices and
  !! arrays from transformed encounter coordinates (z-coordinate,
  !! inter-atomic distances, angles, RMSD), and to read/write matrix
  !! and array files. Several routines use OpenMP for parallelism.
  !!
  !! @author Abraham Muñiz-Chicharro
  !!
  USE maths
  USE OMP_LIB

    contains

    !> Build a matrix of absolute differences of Z-coordinates per encounter.
    !!
    !! Computes the Z coordinate for each encounter (as in `array_z_coord`)
    !! and fills `matrix(i,j) = |z_i - z_j|` for i<j (symmetrically).
    !!
    !! @param[out] matrix      Output (n,n) matrix of pairwise values
    !! @param[out] array       Output array with per-encounter Z values
    !! @param[in]  n           Number of encounters
    !! @param[in]  nb_atoms    Number of atoms in `solute_crds`
    !! @param[in]  xc1, xc2    Solute centers used for translation transformations
    !! @param[in]  trans_vector, rot1, rot2 Transform arrays (n,3)
    !! @param[in]  solute_crds Coordinates of solute atoms
    subroutine matrix_z_coord(matrix, array, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, solute_crds)

      IMPLICIT NONE
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(n, n), intent(out) :: matrix
      real(kind=8), dimension(n), intent(out) :: array
      real(kind=8), dimension(:, :), intent(in) :: solute_crds
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(nb_atoms, 3) :: new_coord
      integer :: i, j, progress_index
      real(kind=8) :: value_i, value_j

      progress_index = 0

      ! Parallelize the outer loop
      !$OMP PARALLEL DO PRIVATE(j, new_coord, value_i, value_j) SCHEDULE(DYNAMIC)
      do i = 1, n
        matrix(i, i) = 0
        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, solute_crds, new_coord)
        value_i = new_coord(1, 3)
        array(i) = value_i
        
        do j = i+1, n
          call update_complex(xc1, xc2, trans_vector(j, :), &
          rot1(j, :), rot2(j, :), nb_atoms, solute_crds, new_coord)
          value_j = new_coord(1, 3)
          
          matrix(i, j) = abs(value_i - value_j)
          matrix(j, i) = matrix(i, j)
        end do

        ! Ensure orderly printing using an atomic increment
        !$OMP CRITICAL
        progress_index = progress_index + 1
        ! Print the progress after the atomic increment
        write(*,*) 'Encounters processed: ', progress_index
        !$OMP END CRITICAL

      end do
      !$OMP END PARALLEL DO

    end subroutine matrix_z_coord

    !> Build matrix from minimum inter-atomic distances between two solutes.
    !!
    !! For each encounter, the second solute is transformed and the minimum
    !! distance to `solute1_crds` is computed; pairwise absolute
    !! differences of these minima populate `matrix`.
    !!
    !! @param[out] matrix         Output (n,n) matrix of pairwise values
    !! @param[out] array          Output (n) array with per-encounter minimum distances
    !! @param[in]  n              Number of encounters
    !! @param[in]  nb_atoms       Number of atoms in solute 2 to transform
    !! @param[in]  xc1, xc2       Solute centers used for translation transformations
    !! @param[in]  trans_vector   Transform translation vectors (n,3)
    !! @param[in]  rot1, rot2     Transform rotation arrays (n,3)
    !! @param[in]  solute1_crds   Coordinates of first solute atoms
    !! @param[in]  solute2_crds   Coordinates of second solute atoms
    subroutine matrix_atoms_dist(matrix, array, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, solute1_crds, solute2_crds)

      IMPLICIT NONE
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(n, n), intent(out) :: matrix
      real(kind=8), dimension(n), intent(out) :: array
      real(kind=8), dimension(:, :), intent(in) :: solute1_crds, solute2_crds
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(:), allocatable :: distances
      real(kind=8), dimension(nb_atoms, 3) :: new_coord
      integer :: i, j, k, tot_coords1, progress_index
      real(kind=8) :: min_i, min_j, dist

      tot_coords1 = size(solute1_crds(:, 1))
      allocate(distances(tot_coords1))
      progress_index = 0

      !$OMP PARALLEL DO PRIVATE(j, k, new_coord, distances, min_i, min_j, dist) SCHEDULE(DYNAMIC)
      do i = 1, n
        matrix(i, i) = 0
        distances = 999999.9

        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, solute2_crds, new_coord)

        do k = 1, tot_coords1
          call calculate_distance(new_coord, solute1_crds(k, :), dist)
          distances(k) = dist
        end do

        min_i = minval(distances)
        array(i) = min_i

        do j = i+1, n
          distances = 999999.9
          call update_complex(xc1, xc2, trans_vector(j, :), &
          rot1(j, :), rot2(j, :), nb_atoms, solute2_crds, new_coord)

          do k = 1, tot_coords1
            call calculate_distance(new_coord, solute1_crds(k, :), dist)
            distances(k) = dist
          end do

          min_j = minval(distances)
          matrix(i, j) = abs(min_i - min_j)
          matrix(j, i) = matrix(i, j)
        end do

        ! Ensure orderly printing using an atomic increment
        !$OMP CRITICAL
        progress_index = progress_index + 1
        ! Print the progress after the atomic increment
        write(*,*) 'Encounters processed: ', progress_index
        !$OMP END CRITICAL

      end do
      !$OMP END PARALLEL DO

      deallocate(distances)

    end subroutine matrix_atoms_dist

  
    !Subroutine matrix_plane_degree
! 
    !> Build matrix of angular differences between encounters.
    !!
    !! Computes the angle (2D or 3D) for each encounter and fills the
    !! pairwise absolute difference matrix.
    !!
    !! @param[out] matrix         Output (n,n) matrix of pairwise angular differences
    !! @param[out] array          Output (n) array with per-encounter angles
    !! @param[in]  n              Number of encounters
    !! @param[in]  nb_atoms       Number of atoms in solute 2 to transform
    !! @param[in]  xc1, xc2       Solute centers used for translation transformations
    !! @param[in]  trans_vector   Transform translation vectors (n,3)
    !! @param[in]  rot1, rot2     Transform rotation arrays (n,3)
    !! @param[in]  point1a        First point defining reference vector
    !! @param[in]  point1b        Second point defining reference vector
    !! @param[in]  point2a        First point of vector to be transformed
    !! @param[in]  point2b        Second point of vector to be transformed
    !! @param[in]  dimensions     Dimensionality of angle calculation (2 or 3)
    subroutine matrix_angle(matrix, array, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, &
      point1a, point1b, point2a, point2b, dimensions)

      IMPLICIT NONE
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(n, n), intent(out) :: matrix
      real(kind=8), dimension(n), intent(out) :: array
      real(kind=8), dimension(:), intent(in) :: point1a, point1b, point2a, point2b
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(nb_atoms, 3) :: new_coord_1, new_coord_2
      real(kind=8), dimension(3) :: v1, v2
      real(kind=8) :: theta1, theta2
      integer :: i, j, progress_index
      real(kind=8), dimension(2, 3) :: solute2_points
      integer, intent(in) :: dimensions

      progress_index = 0

      
      !print *, point1
      !print *, point2
      !print *, point3
      !print *, point4
      !STOP

      ! Initialize the reference vector
      v1(1) = point1b(1) - point1a(1)
      v1(2) = point1b(2) - point1a(2)
      v1(3) = point1b(3) - point1a(3)

      solute2_points(1, :) = point2a(:)
      solute2_points(2, :) = point2b(:)

      ! Parallelizing outer loop with OpenMP
      !$OMP PARALLEL DO PRIVATE(j, new_coord_1, new_coord_2, v2, theta1, theta2) SCHEDULE(DYNAMIC)
      do i = 1, n
        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, solute2_points, new_coord_1)
        
        !print *, point1
        !print *, point2
        !print *, new_coord_1(1, :)
        !print *, new_coord_1(2, :)
        !STOP

        matrix(i, i) = 0
        v2(1) = new_coord_1(1, 1) - new_coord_1(2, 1)
        v2(2) = new_coord_1(1, 2) - new_coord_1(2, 2)
        v2(3) = new_coord_1(1, 3) - new_coord_1(2, 3)

        if ( dimensions .eq. 2 ) then
          call vectors_angle_2D(v1, v2, theta1)
        else if (dimensions .eq. 3) then
          call vectors_angle_3D(v1, v2, theta1)
        else
          print *, "Dimensions indicated is ", dimensions
          print *, "Dimensions implemented are 2 or 3."
          STOP 1
        end if
        
        !print *, 'vector 1', v1
        !print *, 'vector 2', v2
        !print *, 'angle', theta1
        !STOP

        array(i) = theta1

        do j = i+1, n
          call update_complex(xc1, xc2, trans_vector(j, :), &
          rot1(j, :), rot2(j, :), nb_atoms, solute2_points, new_coord_2)
          
          v2(1) = new_coord_2(1, 1) - new_coord_2(2, 1)
          v2(2) = new_coord_2(1, 2) - new_coord_2(2, 2)
          v2(3) = new_coord_2(1, 3) - new_coord_2(2, 3)
          if ( dimensions .eq. 2 ) then
            call vectors_angle_2D(v1, v2, theta2)
          else if (dimensions .eq. 3) then
            call vectors_angle_3D(v1, v2, theta2)
          else
            print *, "Dimensions indicated is ", dimensions
            print *, "Dimensions implemented are 2 or 3."
            STOP 1
          end if
          
          matrix(i, j) = abs(theta1 - theta2)
          matrix(j, i) = matrix(i, j)
        end do

        ! Ensure orderly printing using an atomic increment
        !$OMP CRITICAL
        progress_index = progress_index + 1
        ! Print the progress after the atomic increment
        write(*,*) 'Encounters processed: ', progress_index
        !$OMP END CRITICAL

      end do
      !$OMP END PARALLEL DO

    end subroutine matrix_angle

    
    !> Build RMSD-based pairwise matrix between transformed coordinates.
    !!
    !! Applies transformation to the coordinate set for each encounter and
    !! computes pairwise RMSD values between transformed coordinate sets.
    !!
    !! @param[out] matrix         Output (n,n) matrix of pairwise RMSD values
    !! @param[in]  n              Number of encounters
    !! @param[in]  nb_atoms       Number of atoms in solute 2 to transform
    !! @param[in]  xc1, xc2       Solute centers used for translation transformations
    !! @param[in]  trans_vector   Transform translation vectors (n,3)
    !! @param[in]  rot1, rot2     Transform rotation arrays (n,3)
    !! @param[in]  coord          Coordinates to be transformed
    subroutine matrix_rmsd(matrix, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, coord)
      IMPLICIT NONE
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(n, n), intent(out) :: matrix
      real(kind=8), dimension(:, :), intent(in) :: coord
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(nb_atoms, 3) :: new_coord_1, new_coord_2
      real(kind=8) :: rmsd_value, percentage
      integer :: i, j

      ! Initialize the matrix based on rmsd
      do i = 1, n
        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, coord, new_coord_1)
        matrix(i, i) = 0
        do j = i+1, n
          call update_complex(xc1, xc2, trans_vector(j, :), &
        rot1(j, :), rot2(j, :), nb_atoms, coord, new_coord_2)
          call rmsd(rmsd_value, nb_atoms, new_coord_1, new_coord_2)
          matrix(i, j) = rmsd_value
          matrix(j, i) = matrix(i, j)
        end do
        if (mod(i, n/100) == 0) then
          percentage = (real(i)/real(n))*100 
        end if
        write(*,*) 'Encounters processed: ', i
      end do

    end subroutine matrix_rmsd

    !> Write a numeric matrix to a formatted text file.
    !!
    !! Writes `matrix` with each row as `F10.4` formatted floats. The
    !! output filename is provided by `filename`.
    !!
    !! @param[in]  matrix    Input matrix to write
    !! @param[in]  filename  Output filename
    subroutine write_matrix(matrix, filename)

      IMPLICIT NONE
      real(kind=8), dimension(:, :), intent(in) :: matrix
      character*128, intent(in) :: filename
      character*128 :: fmt
      
      integer :: i, j
      integer :: unit = 15
      integer :: nrows, ncols
      
      nrows = size(matrix(1, :))
      ncols = size(matrix(:, 1))

      open(unit, file=trim(filename), status='replace', action='write')

      write(fmt, '(I0, A)') nrows, 'F10.4'
      fmt = '(' // trim(fmt) // ')'

      ! Write matrix row by row
      do i = 1, nrows
        write(unit, fmt) (matrix(i, j), j = 1, ncols)
      end do
    
      close(unit)
      
    end subroutine write_matrix

    !> Write a numeric array to a formatted text file.
    !!
    !! Writes `array` as a single formatted line using `F10.4` floats.
    !!
    !! @param[in]  array     Input array to write
    !! @param[in]  filename  Output filename
    subroutine write_array(array, filename)

      IMPLICIT NONE
      real(kind=8), dimension(:), intent(in) :: array
      character*128, intent(in) :: filename
      character*128 :: fmt
      
      integer :: i
      integer :: unit = 15
      integer :: nelements
      
      nelements = size(array(:))

      open(unit, file=trim(filename), status='replace', action='write')

      write(fmt, '(I0, A)') nelements, 'F10.4'
      fmt = '(' // trim(fmt) // ')'

      ! Write matrix row by row
      write(unit, fmt) (array(i), i = 1, nelements)
    
      close(unit)
      
    end subroutine write_array

    !> Read a square numeric matrix from a formatted text file.
    !!
    !! Reads the file line-by-line to determine matrix size and parses
    !! rows with the expected `F10.4` formatting.
    !!
    !! @param[out]    matrix    Output matrix read from file (allocated)
    !! @param[inout]  n         Number of encounters (adjusted if file is smaller)
    !! @param[in]     filename  Input filename
    subroutine read_matrix(matrix, n, filename)
      implicit none
      real(kind=8), dimension(:, :), allocatable, intent(out) :: matrix
      integer, intent(inout) :: n
      character*128, intent(in) :: filename
      logical :: file_ex
      integer :: input_matrix, matrix_stat, matrix_size
      integer :: i, j, m
      !integer :: float_length, line_length, actual_length
      character(len=:), allocatable :: line_buffer
      real(kind=8), allocatable :: temp_row(:)
      character*128 :: fmt

      !float_length = 10
      !line_length = n*float_length

      inquire(FILE=filename,EXIST=file_ex)
      if (.NOT.file_ex) then
        write(*,*) "Matrix file not found"
        STOP 1
      end if

      open (input_matrix,FILE=filename,FORM='FORMATTED',STATUS='OLD',IOSTAT=matrix_stat)
      if (matrix_stat.NE.0) then
          write (*,*) "Error opening matrix file"
          STOP 1
      end if

      !Determine the full size of the matrix (m x m)
      m = 0
      do
          read(input_matrix, '(A)', iostat=matrix_stat) line_buffer
          if (matrix_stat /= 0) exit
          m = m + 1
      end do
      rewind(input_matrix)

      ! Check if n is valid
      if ( n .gt. m ) then
        write(*,*) "WARNING: number of encounters provided greater than"
        write(*,*) "the size of the matrix stored in ", filename
        write(*,*) "Setting the number of encounters to the size of the matrix"
        write(*,*) ""
        
        n = m
      end if

      ! Allocate matrix and temporary row
      allocate(matrix(n, n))
      allocate(temp_row(m))
      

      write(fmt, '(I0, A)') m, 'F10.4'
      fmt = '(' // trim(fmt) // ')'


      do i = 1, n
        read(input_matrix, fmt, IOSTAT=matrix_stat) temp_row
        if (matrix_stat /= 0) then
          write(*, *) "Error reading line ", i
          stop
        end if

        !if (i <= n) then
        matrix(i, :) = temp_row(1:n)
        !end if
      end do


      ! Clean up
      close(input_matrix)
      deallocate(temp_row)

    end subroutine read_matrix

    !> Read a numeric array from a formatted text file.
    !!
    !! Parses a single-row formatted array file produced by `write_array`.
    !!
    !! @param[out]    array     Output array read from file (allocated)
    !! @param[inout]  n         Number of encounters (adjusted if file is smaller)
    !! @param[in]     filename  Input filename
    subroutine read_array(array, n, filename)
      implicit none
      real(kind=8), dimension(:), allocatable, intent(out) :: array
      integer, intent(inout) :: n
      character*128, intent(in) :: filename
      logical :: file_ex
      integer :: input_array, array_stat, array_size
      integer :: i, j, m, line_length
      !integer :: float_length, line_length, actual_length
      character(len=:), allocatable :: line_buffer
      real(kind=8), allocatable :: temp_row(:)
      character*128 :: fmt

      !float_length = 10
      !line_length = n*float_length

      inquire(FILE=filename,EXIST=file_ex)
      if (.NOT.file_ex) then
        write(*,*) "Array file not found"
        STOP 1
      end if

      open (input_array,FILE=filename,FORM='FORMATTED',STATUS='OLD',IOSTAT=array_stat)
      if (array_stat.NE.0) then
          write (*,*) "Error opening array file"
          STOP 1
      end if
      
      !Determine the full size of the array (m)
      array_stat = 1
      m = n+1
      do while (array_stat.NE.0)
        m = m -1
        line_length = m*10
        allocate(character(len=line_length) :: line_buffer)
        write(fmt, '(I0, A)') m, 'F10.4'
        read(input_array, '(A)', IOSTAT=array_stat) line_buffer
        deallocate(line_buffer)
        !m = len(line_buffer) / 10 !Floats stored in string have 10.4F format
      end do

      rewind(input_array)

      ! Check if n is valid
      if ( n .gt. m ) then
        write(*,*) "WARNING: number of encounters provided greater than"
        write(*,*) "the size of the matrix stored in ", filename
        write(*,*) "Setting the number of encounters to the size of the array"
        write(*,*) ""
        
        n = m
      end if

      ! Allocate matrix and temporary row
      allocate(array(n))
      allocate(temp_row(m))
      

      write(fmt, '(I0, A)') m, 'F10.4'
      fmt = '(' // trim(fmt) // ')'


      !do i = 1, n
      read(input_array, fmt, IOSTAT=array_stat) temp_row
      if (array_stat /= 0) then
        write(*, *) "Error reading array "
        stop
      end if
      array(:) = temp_row(1:n)
      !end do


      ! Clean up
      close(input_array)
      deallocate(temp_row)

    end subroutine read_array
    
END MODULE mod_matrix