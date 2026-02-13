module mod_cuda

  !! CUDA GPU acceleration module for hierarchical clustering
  !! Provides Fortran interface to CUDA kernels for GPU-accelerated minimum pair search

  use iso_c_binding, only: c_double, c_int, c_loc, C_NULL_PTR
  implicit none

  private
  public :: cuda_find_min_pair

  interface
    subroutine cuda_find_min_pair_c(d_matrix, d_active_points, n, result_i, result_j, result_dist) &
        bind(C, name='cuda_find_min_pair_c')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: d_matrix
      type(c_ptr), value :: d_active_points
      integer(c_int), value :: n
      integer(c_int), intent(out) :: result_i, result_j
      real(c_double), intent(out) :: result_dist
    end subroutine cuda_find_min_pair_c
  end interface

contains

  !> Find minimum distance pair using CUDA GPU acceleration
  !!
  !! Offloads the minimum pair search to NVIDIA GPU using CUDA.
  !! Falls back to CPU if GPU is unavailable or on error.
  !!
  !! @param[in]     matrix          Distance matrix (n x n, column-major)
  !! @param[in]     n               Matrix dimension
  !! @param[out]    min_dist        Minimum distance value
  !! @param[out]    min_i, min_j    Indices of minimum pair (1-based)
  !! @param[in]     active_points   Logical array of active clusters
  
  subroutine cuda_find_min_pair(matrix, n, min_dist, min_i, min_j, active_points)

    implicit none

    real(kind=8), dimension(:, :), intent(in), target :: matrix
    integer, intent(in) :: n
    real(kind=8), intent(out) :: min_dist
    integer, intent(out) :: min_i, min_j
    logical, dimension(:), intent(in) :: active_points

    integer :: i, j
    integer(c_int) :: cuda_i, cuda_j, cuda_n
    real(c_double) :: cuda_dist
    integer(c_int), allocatable, target :: c_active(:)
    integer :: stat

    ! Convert active_points from logical to integer for C
    allocate(c_active(n), stat=stat)
    if (stat /= 0) then
      ! Fallback to CPU on allocation failure
      min_dist = huge(1.0d0)
      min_i = -1
      min_j = -1
      return
    end if

    do i = 1, n
      c_active(i) = merge(1, 0, active_points(i))
    end do

    ! Initialize outputs
    min_dist = huge(1.0d0)
    min_i = -1
    min_j = -1

    ! Convert dimension to C type
    cuda_n = int(n, c_int)

    ! Call CUDA kernel via C interop
    call cuda_find_min_pair_c( &
      c_loc(matrix(1, 1)), &
      c_loc(c_active(1)), &
      cuda_n, &
      cuda_i, cuda_j, &
      cuda_dist &
    )

    ! Convert results (C is 0-based, Fortran is 1-based)
    if (cuda_i >= 0 .and. cuda_j >= 0) then
      min_i = cuda_i + 1
      min_j = cuda_j + 1
      min_dist = cuda_dist
    end if

    deallocate(c_active)

  end subroutine cuda_find_min_pair

end module mod_cuda
