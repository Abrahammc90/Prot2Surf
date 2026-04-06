MODULE mod_matrix

  !> \file mod_matrix.f90
  !! \brief Build/read/write encounter-based matrices and arrays.
  !!
  !! @author Abraham Muñiz-Chicharro
  !! @version 1.0
  !! @date 2026-04-05
  USE maths
  USE OMP_LIB

  ! Workflow summary:
  ! - compute one value per encounter (z, min-distance, angle, rmsd)
  ! - convert values into pairwise matrices when needed
  ! - provide matrix/array read/write helpers used by driver programs

    contains

    
    !> \brief Build RMSD-based pairwise matrix between transformed coordinates.
    !!
    !! Applies transformation to the coordinate set for each encounter and computes pairwise RMSD values between transformed coordinate sets.
    !!
    !! @param[out] matrix         Output (n,n) matrix of pairwise RMSD values
    !! @param[in]  n              Number of encounters
    !! @param[in]  nb_atoms       Number of atoms in solute2 to transform
    !! @param[in]  xc1, xc2       Solute centers for translation
    !! @param[in]  trans_vector   Transform translation vectors (n,3)
    !! @param[in]  rot1, rot2     Transform rotation arrays (n,3)
    !! @param[in]  coord          Coordinates to be transformed
    subroutine matrix_rmsd(matrix, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, coord)
      IMPLICIT NONE
      real(kind=8), dimension(n, n), intent(out) :: matrix
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(:, :), intent(in) :: coord

      real(kind=8), dimension(nb_atoms, 3) :: new_coord_1, new_coord_2
      real(kind=8) :: rmsd_value
      integer :: i, j, progress_step

      progress_step = max(1, n/100)

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
        if (mod(i, progress_step) == 0) then
          write(*,*) 'Encounters processed: ', i
        end if
      end do

    end subroutine matrix_rmsd

    !> \brief Write a numeric matrix to a formatted text file.
    !!
    !! Writes matrix with each row as F10.4 formatted floats. Output filename is provided by filename.
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
      integer :: input_matrix, matrix_stat
      integer :: i, m
      !integer :: float_length, line_length, actual_length
      character(len=1024) :: line_buffer
      real(kind=8), allocatable :: temp_row(:)
      character*128 :: fmt

      !float_length = 10
      !line_length = n*float_length

      inquire(FILE=filename,EXIST=file_ex)
      if (.NOT.file_ex) then
        write(*,*) "Matrix file not found"
        STOP 1
      end if

      input_matrix = 21
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
    
END MODULE mod_matrix