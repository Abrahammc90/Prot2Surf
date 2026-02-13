module mod_cuda

  !! CUDA GPU acceleration module for hierarchical clustering
  !! Uses persistent GPU memory to avoid repeated matrix transfers

  use iso_c_binding, only: c_double, c_int, c_loc, C_NULL_PTR
  implicit none

  private
  public :: cuda_init_clustering, cuda_find_and_merge, cuda_finalize_clustering
  public :: cuda_complete_clustering  ! New function for complete clustering outside loop

  interface
    function cuda_init_clustering_c(d_matrix, d_active_points, d_cluster_size, n) &
        bind(C, name='cuda_init_clustering_c')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: d_matrix
      type(c_ptr), value :: d_active_points
      type(c_ptr), value :: d_cluster_size
      integer(c_int), value :: n
      integer(c_int) :: cuda_init_clustering_c
    end function cuda_init_clustering_c

    function cuda_find_and_merge_c(result_i, result_j, result_dist, linkage_type) &
        bind(C, name='cuda_find_and_merge_c')
      use iso_c_binding
      implicit none
      integer(c_int), intent(out) :: result_i, result_j
      real(c_double), intent(out) :: result_dist
      integer(c_int), value :: linkage_type
      integer(c_int) :: cuda_find_and_merge_c
    end function cuda_find_and_merge_c

    function cuda_complete_clustering_c(h_matrix, h_active_points, h_cluster_size, &
        h_cluster_indexes, h_cluster_count, h_active_clusters, n, dist_threshold, linkage_type) &
        bind(C, name='cuda_complete_clustering_c')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: h_matrix
      type(c_ptr), value :: h_active_points
      type(c_ptr), value :: h_cluster_size
      type(c_ptr), value :: h_cluster_indexes
      type(c_ptr), value :: h_cluster_count
      type(c_ptr), value :: h_active_clusters
      integer(c_int), value :: n
      real(c_double), value :: dist_threshold
      integer(c_int), value :: linkage_type
      integer(c_int) :: cuda_complete_clustering_c
    end function cuda_complete_clustering_c

    subroutine cuda_finalize_clustering_c() &
        bind(C, name='cuda_finalize_clustering_c')
      use iso_c_binding
      implicit none
    end subroutine cuda_finalize_clustering_c
  end interface

contains

  !> Initialize GPU clustering - transfer matrix to GPU once
  !!
  !! @param[in]     matrix          Distance matrix (n x n, column-major)
  !! @param[in]     n               Matrix dimension
  !! @param[in]     active_points   Logical array of active clusters
  !! @param[in]     cluster_size    Integer array of cluster sizes
  !! @return        0 on success, -1 on error
  
  function cuda_init_clustering(matrix, n, active_points, cluster_size) result(status)

    implicit none

    real(kind=8), dimension(:, :), intent(in), target :: matrix
    integer, intent(in) :: n
    logical, dimension(:), intent(in) :: active_points
    integer, dimension(:), intent(in), target :: cluster_size
    integer :: status

    integer :: i
    integer(c_int) :: cuda_n, cuda_status
    integer(c_int), allocatable, target :: c_active(:)
    integer :: stat

    ! Convert active_points from logical to integer for C
    allocate(c_active(n), stat=stat)
    if (stat /= 0) then
      status = -1
      return
    end if

    do i = 1, n
      c_active(i) = merge(1, 0, active_points(i))
    end do

    ! Convert dimension to C type
    cuda_n = int(n, c_int)

    ! Call CUDA initialization
    cuda_status = cuda_init_clustering_c( &
      c_loc(matrix(1, 1)), &
      c_loc(c_active(1)), &
      c_loc(cluster_size(1)), &
      cuda_n &
    )

    deallocate(c_active)
    status = int(cuda_status)

  end function cuda_init_clustering

  !> Find minimum pair and merge on GPU
  !!
  !! @param[out]    min_dist        Minimum distance value
  !! @param[out]    min_i, min_j    Indices of minimum pair (1-based)
  !! @param[in]     linkage_str     Linkage type ('min', 'max', or 'mean')
  !! @return        0 on success, -1 on error
  
  function cuda_find_and_merge(min_dist, min_i, min_j, linkage_str) result(status)

    implicit none

    real(kind=8), intent(out) :: min_dist
    integer, intent(out) :: min_i, min_j
    character(len=*), intent(in) :: linkage_str
    integer :: status

    integer(c_int) :: cuda_i, cuda_j, cuda_status, linkage_type
    real(c_double) :: cuda_dist

    ! Convert linkage string to integer code
    select case (trim(linkage_str))
      case ('min')
        linkage_type = 0
      case ('max')
        linkage_type = 1
      case ('mean')
        linkage_type = 2
      case default
        linkage_type = 2  ! Default to mean
    end select

    ! Call CUDA find and merge
    cuda_status = cuda_find_and_merge_c( &
      cuda_i, cuda_j, &
      cuda_dist, &
      linkage_type &
    )

    if (cuda_status /= 0) then
      status = -1
      min_i = -1
      min_j = -1
      min_dist = huge(1.0d0)
      return
    end if

    ! Convert results (C is 0-based, Fortran is 1-based)
    if (cuda_i >= 0 .and. cuda_j >= 0) then
      min_i = cuda_i + 1
      min_j = cuda_j + 1
      min_dist = cuda_dist
      status = 0
    else
      min_i = -1
      min_j = -1
      min_dist = huge(1.0d0)
      status = -1
    end if

  end function cuda_find_and_merge

  !> Perform complete GPU clustering - all iterations on GPU
  !!
  !! @param[in]     matrix              Distance matrix (n x n, column-major)
  !! @param[in]     n                   Matrix dimension
  !! @param[in]     dist_threshold      Distance threshold for stopping
  !! @param[in]     linkage_str         Linkage type ('min', 'max', or 'mean')
  !! @param[in,out] cluster_indexes     Output cluster membership (n x n)
  !! @param[in,out] cluster_count       Output cluster sizes (n)
  !! @param[in,out] active_points       Output active clusters (n)
  !! @param[in,out] cluster_size        Cluster sizes for linkage computation
  !! @return        0 on success, -1 on error
  
  function cuda_complete_clustering(matrix, n, dist_threshold, linkage_str, &
                                     cluster_indexes, cluster_count, active_points, cluster_size) result(status)

    implicit none

    real(kind=8), dimension(:, :), intent(in), target :: matrix
    integer, intent(in) :: n
    real(kind=8), intent(in) :: dist_threshold
    character(len=*), intent(in) :: linkage_str
    integer, dimension(:, :), intent(inout), target :: cluster_indexes
    integer, dimension(:), intent(inout), target :: cluster_count
    logical, dimension(:), intent(inout) :: active_points
    integer, dimension(:), intent(in), target :: cluster_size
    integer :: status

    integer :: i
    integer(c_int) :: cuda_n, cuda_status, linkage_type
    integer(c_int), allocatable, target :: c_active_in(:), c_active_out(:)
    real(c_double) :: cuda_threshold
    integer :: stat

    ! Convert active_points from logical to integer for C
    allocate(c_active_in(n), c_active_out(n), stat=stat)
    if (stat /= 0) then
      status = -1
      return
    end if

    do i = 1, n
      c_active_in(i) = merge(1, 0, active_points(i))
    end do

    ! Convert linkage string to integer code
    select case (trim(linkage_str))
      case ('min')
        linkage_type = 0
      case ('max')
        linkage_type = 1
      case ('mean')
        linkage_type = 2
      case default
        linkage_type = 2  ! Default to mean
    end select

    ! Convert dimension and threshold to C types
    cuda_n = int(n, c_int)
    cuda_threshold = real(dist_threshold, c_double)

    ! Call CUDA complete clustering
    cuda_status = cuda_complete_clustering_c( &
      c_loc(matrix(1, 1)), &
      c_loc(c_active_in(1)), &
      c_loc(cluster_size(1)), &
      c_loc(cluster_indexes(1, 1)), &
      c_loc(cluster_count(1)), &
      c_loc(c_active_out(1)), &
      cuda_n, &
      cuda_threshold, &
      linkage_type &
    )

    ! Convert active clusters back to logical
    do i = 1, n
      active_points(i) = (c_active_out(i) /= 0)
    end do

    deallocate(c_active_in, c_active_out)
    status = int(cuda_status)

  end function cuda_complete_clustering

  !> Finalize GPU clustering - cleanup GPU memory
  
  subroutine cuda_finalize_clustering()

    implicit none

    call cuda_finalize_clustering_c()

  end subroutine cuda_finalize_clustering

end module mod_cuda
