MODULE mod_array

  !> \file mod_array.f90
  !! \brief Build/read/write encounter-based arrays and per-encounter metrics.
  !!
  !! Provides routines to compute z-coordinate, distance, and angle-based arrays
  !! from encounter transformations, as well as array I/O and sorting helpers.
  !!
  !! @author Abraham Muñiz-Chicharro
  !! @version 1.0
  !! @date 2026-04-05
  USE maths
  USE mod_assoc
  USE OMP_LIB

  ! Workflow summary:
  ! - compute one scalar value per encounter from transforms and atom groups
  ! - keep encounter metadata aligned while sorting output arrays
  ! - provide array read/write and merge-sort helpers for stable ordering

    contains

    !> Compute the Z-coordinate (along transformed z-axis) for each encounter and sort.
    !!
    !! Computes the Z coordinate for each encounter and fills the output array.
    !! The resulting array is sorted from minimum to maximum value.
    !!
    !! @param[in,out] complexes   Association object reordered consistently with sorted outputs
    !! @param[out] array_sorted   Output array of length `n` with z-coordinates (sorted)
    !! @param[in]  n              Number of encounters
    !! @param[in]  nb_atoms       Number of atoms in `solute_crds`
    !! @param[in]  xc1, xc2       Solute centers used for translation transformations
    !! @param[in]  trans_vector   Translation vectors for each encounter (n, 3)
    !! @param[in]  rot1, rot2     Rotation vectors for each encounter (n, 3)
    !! @param[in]  solute_crds    Coordinates of solute atoms (nb_atoms, 3)
    subroutine array_z_coord(complexes, array_sorted, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, solute_crds)

      IMPLICIT NONE
      type(type_assoc_file), intent(inout) :: complexes
      real(kind=8), dimension(n), intent(out) :: array_sorted
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(:, :), intent(in) :: solute_crds

      real(kind=8), dimension(n) :: array
      real(kind=8), dimension(n) :: tmp_arr
      character(len=217), dimension(n) :: tmp_comp
      real(kind=8), dimension(nb_atoms, 3) :: new_coord
      integer :: i, ii, jj, t, num_threads, chunk_size
      integer :: start_idx, end_idx, left, mid, right, merge_size
      real(kind=8) :: value_i
      character(len=217) :: temp_complex
      real(kind=8), dimension(3) :: temp_trans, temp_rot1, temp_rot2

      ! Step 1: Compute z-coordinates in parallel
      !$OMP PARALLEL DO PRIVATE(new_coord, value_i) SCHEDULE(DYNAMIC)
      do i = 1, n
        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, solute_crds, new_coord)
        value_i = new_coord(1, 3)
        array(i) = value_i
      end do
      !$OMP END PARALLEL DO

      ! Step 2: Copy array values
      !$OMP PARALLEL DO SCHEDULE(STATIC)
      do i = 1, n
        array_sorted(i) = array(i)
      end do
      !$OMP END PARALLEL DO

      ! Step 3: Parallel merge sort
      ! 3a. Determine number of threads and chunk size
      num_threads = 1
      !$OMP PARALLEL
      !$OMP SINGLE
      num_threads = omp_get_num_threads()
      !$OMP END SINGLE
      !$OMP END PARALLEL

      chunk_size = (n + num_threads - 1) / num_threads

      ! 3b. Sort each chunk independently in parallel (insertion sort)
      !$OMP PARALLEL DO PRIVATE(start_idx, end_idx, ii, jj, value_i, &
      !$OMP& temp_complex, temp_trans, temp_rot1, temp_rot2) SCHEDULE(STATIC)
      do t = 0, num_threads - 1
        start_idx = t * chunk_size + 1
        end_idx = min((t + 1) * chunk_size, n)
        if (start_idx <= n) then
          do ii = start_idx + 1, end_idx
            value_i = array_sorted(ii)
            temp_complex = complexes%lines(ii)
            temp_trans = complexes%trans_vector(ii, :)
            temp_rot1 = complexes%rot1(ii, :)
            temp_rot2 = complexes%rot2(ii, :)
            jj = ii - 1
            do while (jj >= start_idx .and. array_sorted(jj) > value_i)
              array_sorted(jj + 1) = array_sorted(jj)
              complexes%lines(jj + 1) = complexes%lines(jj)
              complexes%trans_vector(jj + 1, :) = complexes%trans_vector(jj, :)
              complexes%rot1(jj + 1, :) = complexes%rot1(jj, :)
              complexes%rot2(jj + 1, :) = complexes%rot2(jj, :)
              jj = jj - 1
            end do
            array_sorted(jj + 1) = value_i
            complexes%lines(jj + 1) = temp_complex
            complexes%trans_vector(jj + 1, :) = temp_trans
            complexes%rot1(jj + 1, :) = temp_rot1
            complexes%rot2(jj + 1, :) = temp_rot2
          end do
        end if
      end do
      !$OMP END PARALLEL DO

      ! 3c. Iteratively merge sorted chunks in parallel
      merge_size = chunk_size
      do while (merge_size < n)
        !$OMP PARALLEL DO PRIVATE(left, mid, right) SCHEDULE(DYNAMIC)
        do t = 0, (n - 1) / (2 * merge_size)
          left = t * 2 * merge_size + 1
          mid = min(left + merge_size - 1, n)
          right = min(left + 2 * merge_size - 1, n)
          if (mid < right) then
            call merge_sorted_segments(array_sorted, complexes, &
              tmp_arr, tmp_comp, left, mid, right)
          end if
        end do
        !$OMP END PARALLEL DO
        merge_size = merge_size * 2
      end do

    end subroutine array_z_coord

    !> Compute minimum inter-atomic distance between two solutes for each encounter and sort.
    !!
    !! For each encounter, the second solute is transformed and the minimum
    !! distance to `solute1_crds` is computed and stored in the output array.
    !! The resulting array is sorted from minimum to maximum value.
    !!
    !! @param[in,out] complexes      Association object reordered consistently with sorted outputs
    !! @param[out] array_sorted   Output array of length `n` with minimum distances (sorted)
    !! @param[in]  n              Number of encounters
    !! @param[in]  nb_atoms       Number of atoms in `solute2_crds`
    !! @param[in]  xc1, xc2       Solute centers used for translation transformations
    !! @param[in]  trans_vector   Translation vectors for each encounter (n, 3)
    !! @param[in]  rot1, rot2     Rotation vectors for each encounter (n, 3)
    !! @param[in]  solute1_crds   Coordinates of solute1 atoms
    !! @param[in]  solute2_crds   Coordinates of solute2 atoms
    subroutine array_atoms_dist(complexes, array_sorted, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, solute1_crds, solute2_crds)

      IMPLICIT NONE
      type(type_assoc_file), intent(inout) :: complexes
      real(kind=8), dimension(n), intent(out) :: array_sorted
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(:, :), intent(in) :: solute1_crds, solute2_crds

      real(kind=8), dimension(:), allocatable :: distances
      real(kind=8), dimension(nb_atoms, 3) :: new_coord
      real(kind=8), dimension(n) :: array
      real(kind=8), dimension(n) :: tmp_arr
      character(len=217), dimension(n) :: tmp_comp
      integer :: i, j, k, tot_coords1, progress_index
      real(kind=8) :: min_i, min_j, dist
      integer :: ii, jj, t, num_threads, chunk_size
      integer :: start_idx, end_idx, left, mid, right, merge_size
      real(kind=8) :: value_i
      character(len=217) :: temp_complex
      real(kind=8), dimension(3) :: temp_trans, temp_rot1, temp_rot2

      tot_coords1 = size(solute1_crds(:, 1))
      allocate(distances(tot_coords1))
      progress_index = 0

      ! Compute minimum-distance for each encounter.
      !$OMP PARALLEL DO PRIVATE(j, k, new_coord, distances, min_i, min_j, dist) SCHEDULE(DYNAMIC)
      do i = 1, n
        distances = 999999.9

        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, solute2_crds, new_coord)

        do k = 1, tot_coords1
          call calculate_distance(new_coord, solute1_crds(k, :), dist)
          distances(k) = dist
        end do

        min_i = minval(distances)
        array(i) = min_i

      end do
      !$OMP END PARALLEL DO

      deallocate(distances)

      ! Step 2: Copy arrays in parallel
      !$OMP PARALLEL DO SCHEDULE(STATIC)
      do i = 1, n
        array_sorted(i) = array(i)
      end do
      !$OMP END PARALLEL DO

      ! Step 3: Parallel merge sort
      ! 3a. Determine number of threads and chunk size
      num_threads = 1
      !$OMP PARALLEL
      !$OMP SINGLE
      num_threads = omp_get_num_threads()
      !$OMP END SINGLE
      !$OMP END PARALLEL

      chunk_size = (n + num_threads - 1) / num_threads

      ! 3b. Sort each chunk independently in parallel (insertion sort)
      !$OMP PARALLEL DO PRIVATE(start_idx, end_idx, ii, jj, value_i, &
      !$OMP& temp_complex, temp_trans, temp_rot1, temp_rot2) SCHEDULE(STATIC)
      do t = 0, num_threads - 1
        start_idx = t * chunk_size + 1
        end_idx = min((t + 1) * chunk_size, n)
        if (start_idx <= n) then
          do ii = start_idx + 1, end_idx
            value_i = array_sorted(ii)
            temp_complex = complexes%lines(ii)
            temp_trans = complexes%trans_vector(ii, :)
            temp_rot1 = complexes%rot1(ii, :)
            temp_rot2 = complexes%rot2(ii, :)
            jj = ii - 1
            do while (jj >= start_idx .and. array_sorted(jj) > value_i)
              array_sorted(jj + 1) = array_sorted(jj)
              complexes%lines(jj + 1) = complexes%lines(jj)
              complexes%trans_vector(jj + 1, :) = complexes%trans_vector(jj, :)
              complexes%rot1(jj + 1, :) = complexes%rot1(jj, :)
              complexes%rot2(jj + 1, :) = complexes%rot2(jj, :)
              jj = jj - 1
            end do
            array_sorted(jj + 1) = value_i
            complexes%lines(jj + 1) = temp_complex
            complexes%trans_vector(jj + 1, :) = temp_trans
            complexes%rot1(jj + 1, :) = temp_rot1
            complexes%rot2(jj + 1, :) = temp_rot2
          end do
        end if
      end do
      !$OMP END PARALLEL DO

      ! 3c. Iteratively merge sorted chunks in parallel
      merge_size = chunk_size
      do while (merge_size < n)
        !$OMP PARALLEL DO PRIVATE(left, mid, right) SCHEDULE(DYNAMIC)
        do t = 0, (n - 1) / (2 * merge_size)
          left = t * 2 * merge_size + 1
          mid = min(left + merge_size - 1, n)
          right = min(left + 2 * merge_size - 1, n)
          if (mid < right) then
            call merge_sorted_segments(array_sorted, complexes, &
              tmp_arr, tmp_comp, left, mid, right)
          end if
        end do
        !$OMP END PARALLEL DO
        merge_size = merge_size * 2
      end do

    end subroutine array_atoms_dist

  
    !Subroutine matrix_plane_degree
! 
    !> Compute angle between two vectors (defined by atom groups) per encounter and sort.
    !!
    !! Computes the angle between two vectors (defined by atom groups) for each encounter.
    !! The resulting array is sorted from minimum to maximum value.
    !!
    !! @param[in,out] complexes      Association object reordered consistently with sorted outputs
    !! @param[out] array_sorted   Output array of length `n` with angles (sorted)
    !! @param[in]  n              Number of encounters
    !! @param[in]  nb_atoms       Number of atoms used for transformation
    !! @param[in]  xc1, xc2       Solute centers used for translation transformations
    !! @param[in]  trans_vector   Translation vectors for each encounter (n, 3)
    !! @param[in]  rot1, rot2     Rotation vectors for each encounter (n, 3)
    !! @param[in]  point1a        First point defining reference vector
    !! @param[in]  point1b        Second point defining reference vector
    !! @param[in]  point2a        First point of vector to be transformed
    !! @param[in]  point2b        Second point of vector to be transformed
    !! @param[in]  dimensions     Use 2 or 3 to select projection
    subroutine array_angle(complexes, array_sorted, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, &
      point1a, point1b, point2a, point2b, dimensions)

      IMPLICIT NONE
      type(type_assoc_file), intent(inout) :: complexes
      real(kind=8), dimension(n), intent(out) :: array_sorted
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(:), intent(in) :: point1a, point1b, point2a, point2b
      integer, intent(in) :: dimensions

      real(kind=8), dimension(nb_atoms, 3) :: new_coord_1, new_coord_2
      real(kind=8), dimension(3) :: v1, v2
      real(kind=8) :: theta1, theta2
      integer :: i, j
      real(kind=8), dimension(2, 3) :: solute2_points
      real(kind=8), dimension(n) :: array
      real(kind=8), dimension(n) :: tmp_arr
      character(len=217), dimension(n) :: tmp_comp
      integer :: ii, jj, t, num_threads, chunk_size
      integer :: start_idx, end_idx, left, mid, right, merge_size
      real(kind=8) :: value_i
      character(len=217) :: temp_complex
      real(kind=8), dimension(3) :: temp_trans, temp_rot1, temp_rot2

      
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
        v2(1) = new_coord_1(1, 1) - new_coord_1(2, 1)
        v2(2) = new_coord_1(1, 2) - new_coord_1(2, 2)
        v2(3) = new_coord_1(1, 3) - new_coord_1(2, 3)

        if ( dimensions .eq. 2 ) then
          call vectors_angle_2D(v1, v2, theta1)
        else if (dimensions .eq. 3) then
          call vectors_angle_3D(v1, v2, theta1)
        else
          print *, '[mod_array] ERROR: Dimensions indicated is', dimensions
          print *, '[mod_array] Only 2D or 3D implemented.'
          STOP 1
        end if
        
        !print *, 'vector 1', v1
        !print *, 'vector 2', v2
        !print *, 'angle', theta1
        !STOP

        array(i) = theta1

      end do
      !$OMP END PARALLEL DO

      ! Step 2: Copy arrays in parallel
      !$OMP PARALLEL DO SCHEDULE(STATIC)
      do i = 1, n
        array_sorted(i) = array(i)
      end do
      !$OMP END PARALLEL DO

      ! Step 3: Parallel merge sort
      ! 3a. Determine number of threads and chunk size
      num_threads = 1
      !$OMP PARALLEL
      !$OMP SINGLE
      num_threads = omp_get_num_threads()
      !$OMP END SINGLE
      !$OMP END PARALLEL

      chunk_size = (n + num_threads - 1) / num_threads

      ! 3b. Sort each chunk independently in parallel (insertion sort)
      !$OMP PARALLEL DO PRIVATE(start_idx, end_idx, ii, jj, value_i, &
      !$OMP& temp_complex, temp_trans, temp_rot1, temp_rot2) SCHEDULE(STATIC)
      do t = 0, num_threads - 1
        start_idx = t * chunk_size + 1
        end_idx = min((t + 1) * chunk_size, n)
        if (start_idx <= n) then
          do ii = start_idx + 1, end_idx
            value_i = array_sorted(ii)
            temp_complex = complexes%lines(ii)
            temp_trans = complexes%trans_vector(ii, :)
            temp_rot1 = complexes%rot1(ii, :)
            temp_rot2 = complexes%rot2(ii, :)
            jj = ii - 1
            do while (jj >= start_idx .and. array_sorted(jj) > value_i)
              array_sorted(jj + 1) = array_sorted(jj)
              complexes%lines(jj + 1) = complexes%lines(jj)
              complexes%trans_vector(jj + 1, :) = complexes%trans_vector(jj, :)
              complexes%rot1(jj + 1, :) = complexes%rot1(jj, :)
              complexes%rot2(jj + 1, :) = complexes%rot2(jj, :)
              jj = jj - 1
            end do
            array_sorted(jj + 1) = value_i
            complexes%lines(jj + 1) = temp_complex
            complexes%trans_vector(jj + 1, :) = temp_trans
            complexes%rot1(jj + 1, :) = temp_rot1
            complexes%rot2(jj + 1, :) = temp_rot2
          end do
        end if
      end do
      !$OMP END PARALLEL DO

      ! 3c. Iteratively merge sorted chunks in parallel
      merge_size = chunk_size
      do while (merge_size < n)
        !$OMP PARALLEL DO PRIVATE(left, mid, right) SCHEDULE(DYNAMIC)
        do t = 0, (n - 1) / (2 * merge_size)
          left = t * 2 * merge_size + 1
          mid = min(left + merge_size - 1, n)
          right = min(left + 2 * merge_size - 1, n)
          if (mid < right) then
            call merge_sorted_segments(array_sorted, complexes, &
              tmp_arr, tmp_comp, left, mid, right)
          end if
        end do
        !$OMP END PARALLEL DO
        merge_size = merge_size * 2
      end do

    end subroutine array_angle

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
      integer :: input_array, array_stat
      integer :: m, line_length
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

      input_array = 22
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
    !> Merge two adjacent sorted segments of arrays in place.
    !!
    !! Merges `arr(left:mid)` and `arr(mid+1:right)` into a single sorted
    !! segment `arr(left:right)`, reordering complexes and transform fields in the same way.
    !!
    !! @param[in,out] arr      Array of values being sorted
    !! @param[in,out] complexes Association object reordered alongside `arr`
    !! @param[out]    tmp_arr  Temporary buffer for values (at least size `right`)
    !! @param[out]    tmp_comp Temporary buffer for companion (at least size `right`)
    !! @param[in]     left     Start index of the first segment
    !! @param[in]     mid      End index of the first segment
    !! @param[in]     right    End index of the second segment
    subroutine merge_sorted_segments(arr, complexes, tmp_arr, tmp_comp, left, mid, right)

      IMPLICIT NONE
      real(kind=8), dimension(:), intent(inout) :: arr
      type(type_assoc_file), intent(inout) :: complexes
      real(kind=8), dimension(:), intent(inout) :: tmp_arr
      character(len=217), dimension(:), intent(inout) :: tmp_comp
      integer, intent(in) :: left, mid, right
      integer :: i, j, k
      real(kind=8), dimension(right, 3) :: tmp_trans, tmp_rot1, tmp_rot2

      ! Copy the segments to temporary arrays
      do i = left, right
        tmp_arr(i) = arr(i)
        tmp_comp(i) = complexes%lines(i)
        tmp_trans(i, :) = complexes%trans_vector(i, :)
        tmp_rot1(i, :) = complexes%rot1(i, :)
        tmp_rot2(i, :) = complexes%rot2(i, :)
      end do

      ! Merge back from temporary arrays into original
      i = left
      j = mid + 1
      k = left

      do while (i <= mid .and. j <= right)
        if (tmp_arr(i) <= tmp_arr(j)) then
          arr(k) = tmp_arr(i)
          complexes%lines(k) = tmp_comp(i)
          complexes%trans_vector(k, :) = tmp_trans(i, :)
          complexes%rot1(k, :) = tmp_rot1(i, :)
          complexes%rot2(k, :) = tmp_rot2(i, :)
          i = i + 1
        else
          arr(k) = tmp_arr(j)
          complexes%lines(k) = tmp_comp(j)
          complexes%trans_vector(k, :) = tmp_trans(j, :)
          complexes%rot1(k, :) = tmp_rot1(j, :)
          complexes%rot2(k, :) = tmp_rot2(j, :)
          j = j + 1
        end if
        k = k + 1
      end do

      ! Copy remaining elements from the left segment
      do while (i <= mid)
        arr(k) = tmp_arr(i)
        complexes%lines(k) = tmp_comp(i)
        complexes%trans_vector(k, :) = tmp_trans(i, :)
        complexes%rot1(k, :) = tmp_rot1(i, :)
        complexes%rot2(k, :) = tmp_rot2(i, :)
        i = i + 1
        k = k + 1
      end do

      ! Remaining elements from the right segment are already in place

    end subroutine merge_sorted_segments
    
END MODULE mod_array