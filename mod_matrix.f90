MODULE mod_matrix

    USE maths

    contains

    subroutine matrix_z(matrix, array, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, coord)
    
      IMPLICIT NONE
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(n, n), intent(out) :: matrix
      real(kind=8), dimension(n), intent(out) :: array
      real(kind=8), dimension(:, :), intent(in) :: coord
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(nb_atoms, 3) :: new_coord
      integer :: i, j
      real(kind=8) :: value_i, value_j

      !write(*,*) trans_vector(1, :)
      !STOP

      ! Initialize the matrix based on z-distance
      do i = 1, n
        matrix(i, i) = 0
        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, coord, new_coord)
        value_i = new_coord(1, 3)
        array(i) = value_i
        do j = i+1, n
          call update_complex(xc1, xc2, trans_vector(j, :), &
        rot1(j, :), rot2(j, :), nb_atoms, coord, new_coord)
          value_j = new_coord(1, 3)
          !write(*,*) 'new coord j', new_coord
          !write(*,*) value_i, value_j
          !write(*,*) 'trans i', trans_vector(i, :)
          !write(*,*) 'trans j', trans_vector(j, :)
          !write(*,*) 'rot1 i', rot1(i, :)
          !write(*,*) 'rot2 i', rot2(i, :)
          !write(*,*) 'rot1 j', rot1(j, :)
          !write(*,*) 'rot2 j', rot2(j, :)
          !if (j == 12) STOP
          matrix(i, j) = abs(value_i - value_j) !function
          matrix(j, i) = matrix(i, j)
        end do
      end do
    end subroutine matrix_z

    subroutine matrix_dist(matrix, array, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, protein_points, surface_points)
    
      IMPLICIT NONE
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(n, n), intent(out) :: matrix
      real(kind=8), dimension(n), intent(out) :: array
      real(kind=8), dimension(:, :), intent(in) :: protein_points, surface_points
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(:), allocatable :: C1_C4_distances
      real(kind=8), dimension(nb_atoms, 3) :: new_coord
      integer :: i, j, k, tot_coords
      real(kind=8) :: value_i, value_j, dist

      tot_coords = size(surface_points)
      allocate(C1_C4_distances(tot_coords))
      C1_C4_distances = 0

      ! Initialize the matrix based on the closest C1/C4 atom
      do i = 1, n
        matrix(i, i) = 0
        C1_C4_distances = 0
        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, protein_points, new_coord)
        do k = 1, tot_coords
          call calculate_distance(new_coord, surface_points(k, :), dist)
          C1_C4_distances(k) = dist
        end do
        value_i = minval(C1_C4_distances)
        array(i) = value_i
        do j = i+1, n
          C1_C4_distances = 0
          call update_complex(xc1, xc2, trans_vector(j, :), &
        rot1(j, :), rot2(j, :), nb_atoms, protein_points, new_coord)
          do k = 1, tot_coords
            call calculate_distance(new_coord, surface_points(k, :), dist)
            C1_C4_distances(k) = dist
          end do
          value_j = minval(C1_C4_distances)
          matrix(i, j) = abs(value_i - value_j) !function
          matrix(j, i) = matrix(i, j)
        end do
        write(*,*) 'Encounters processed: ', i
      end do
    end subroutine matrix_dist

    !Subroutine matrix_plane_degree
!
    subroutine matrix_chain_degree(matrix, array, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, protein_points, surface_points)
      
      IMPLICIT NONE
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(n, n), intent(out) :: matrix
      real(kind=8), dimension(n), intent(out) :: array
      real(kind=8), dimension(:, :), intent(in) :: protein_points, surface_points
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(nb_atoms, 3) :: new_coord_1, new_coord_2
      real(kind=8), dimension(3) :: v1, v2
      real(kind=8) :: theta1, theta2
      integer :: i, j

      ! Initialize the matrix based on rmsd
      v1(1) = surface_points(1, 1) - surface_points(2, 1)
      v1(2) = surface_points(1, 2) - surface_points(2, 2)
      v1(3) = surface_points(1, 3) - surface_points(2, 3)
      do i = 1, n
        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, protein_points, new_coord_1)
        
        matrix(i, i) = 0
        v2(1) = new_coord_1(1, 1) -  new_coord_1(2, 1)
        v2(2) = new_coord_1(1, 2) -  new_coord_1(2, 2)
        v2(3) = new_coord_1(1, 3) -  new_coord_1(2, 3)
        call vectors_angle(v1, v2, theta1)

        array(i) = theta1

        do j = i+1, n
          call update_complex(xc1, xc2, trans_vector(j, :), &
        rot1(j, :), rot2(j, :), nb_atoms, protein_points, new_coord_2)
          
          v2(1) = new_coord_2(1, 1) -  new_coord_2(2, 1)
          v2(2) = new_coord_2(1, 2) -  new_coord_2(2, 2)
          v2(3) = new_coord_2(1, 3) -  new_coord_2(2, 3)
          call vectors_angle(v1, v2, theta2)
          
          matrix(i, j) = abs(theta1 - theta2)
          matrix(j, i) = matrix(i, j)
        end do
      end do

    end subroutine matrix_chain_degree

    !subroutine matrix_residue_dist(matrix, array, n, nb_atoms, &
    !  xc1, xc2, trans_vector, rot1, rot2, coord)
    !
    !  IMPLICIT NONE
    !  real(kind=8), dimension(3), intent(in) :: xc1, xc2
    !  real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
    !  real(kind=8), dimension(n, n), intent(out) :: matrix
    !  real(kind=8), dimension(n), intent(out) :: array
    !  real(kind=8), dimension(:, :), intent(in) :: coord
    !  integer, intent(in) :: n, nb_atoms
    !  real(kind=8), dimension(nb_atoms, 3) :: new_coord
    !  integer :: i, j
    !  real(kind=8) :: value_i, value_j
!
    !  !write(*,*) trans_vector(1, :)
    !  !STOP
!
    !  ! Initialize the matrix based on z-distance
    !  do i = 1, n
    !    matrix(i, i) = 0
    !    call update_complex(xc1, xc2, trans_vector(i, :), &
    !    rot1(i, :), rot2(i, :), nb_atoms, coord, new_coord)
    !    value_i = new_coord(1, 3)
    !    array(i) = value_i
    !    do j = i+1, n
    !      call update_complex(xc1, xc2, trans_vector(j, :), &
    !    rot1(j, :), rot2(j, :), nb_atoms, coord, new_coord)
    !      value_j = new_coord(1, 3)
    !      !write(*,*) 'new coord j', new_coord
    !      !write(*,*) value_i, value_j
    !      !write(*,*) 'trans i', trans_vector(i, :)
    !      !write(*,*) 'trans j', trans_vector(j, :)
    !      !write(*,*) 'rot1 i', rot1(i, :)
    !      !write(*,*) 'rot2 i', rot2(i, :)
    !      !write(*,*) 'rot1 j', rot1(j, :)
    !      !write(*,*) 'rot2 j', rot2(j, :)
    !      !if (j == 12) STOP
    !      matrix(i, j) = abs(value_i - value_j) !function
    !      matrix(j, i) = matrix(i, j)
    !    end do
    !  end do
    !end subroutine matrix_residue_dist
    
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
        !print *, trans_vector(i, :)
        !print *, rot1(i, :)
        !print *, rot2(i, :)
        !print *, ""
        matrix(i, i) = 0
        do j = i+1, n
          call update_complex(xc1, xc2, trans_vector(j, :), &
        rot1(j, :), rot2(j, :), nb_atoms, coord, new_coord_2)
          !print *, trans_vector(j, :)
          !print *, rot1(j, :)
          !print *, rot2(j, :)
          !print *, ""
          !call rmsd(rmsd2, dist_max, & 
          !trans_vector(i, :), rot1(i, :), rot2(i, :), &
          !trans_vector(j, :), rot1(j, :), rot2(j, :))
          call rmsd(rmsd_value, nb_atoms, new_coord_1, new_coord_2)
          !print *, new_coord_1(1, :)
          !print *, new_coord_2(1, :)
          !print *, rmsd_value
          !stop
          matrix(i, j) = rmsd_value
          matrix(j, i) = matrix(i, j)
        end do
        if (mod(i, n/100) == 0) then
          percentage = (real(i)/real(n))*100 
          !write(*,'(F5.1, A)') percentage, '%'
        end if
        !if (mod(i, n/1000) == 0) write(*,*) i
      end do

    end subroutine matrix_rmsd

    subroutine write_matrix(matrix, filename)

      IMPLICIT NONE
      real(kind=8), dimension(:, :), intent(in) :: matrix
      character*128, intent(in) :: filename
      character*128 :: fmt
      
      integer :: i, j
      integer :: unit
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

    subroutine write_array(array, filename)

      IMPLICIT NONE
      real(kind=8), dimension(:), intent(in) :: array
      character*128, intent(in) :: filename
      character*128 :: fmt
      
      integer :: i
      integer :: unit
      integer :: nelements
      
      nelements = size(array(:))

      open(unit, file=trim(filename), status='replace', action='write')

      write(fmt, '(I0, A)') nelements, 'F10.4'
      fmt = '(' // trim(fmt) // ')'

      ! Write matrix row by row
      write(unit, fmt) (array(i), i = 1, nelements)
    
      close(unit)
      
    end subroutine write_array

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