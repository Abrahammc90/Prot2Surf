!> \file mod_clust_algorithm.f90
!! \brief Hierarchical clustering algorithms and output helpers.
!!
module mod_clust_algorithm

  !> Module: mod_clust_algorithm
  !!
  !! Clustering algorithms and helpers for processing encounter matrices.
  !!
  !! Provides hierarchical clustering (min, max, mean linkage) for distance/similarity matrices,
  !! sorting, and output routines. Supports OpenMP parallelism.
  !!
  !! @author Abraham Muñiz-Chicharro
  !! @version 1.0
  !! @date 2026-04-05

  USE read_input
  USE OMP_LIB

    ! Workflow summary:
    ! - cluster encounters from matrix or array inputs
    ! - track parent/count arrays during merges
    ! - support min/max/mean linkage
    ! - write cluster members and summary files to disk

  contains

    !> Perform hierarchical clustering on a distance matrix.
    !!
    !! This routine performs hierarchical clustering using minimum, maximum,
    !! or mean linkage on the input `matrix` (shape `(n,n)`) of
    !! pairwise distances or similarity scores. The routine mutates `matrix`
    !! during merging and terminates merging when clusters are farther apart
    !! than an internally computed threshold. Results are written to files
    !! named using `output_name` and complexes information is optionally
    !! recorded via the `complexes` structure.
    !!
    !! @param[in,out] matrix        Pairwise distance/similarity matrix (n, n)
    !!                              This matrix is updated during clustering
    !! @param[in]     n             Number of elements / encounters
    !! @param[in]     linkage_type  Type of linkage: 'min', 'max', or 'mean'
    !! @param[in]     output_name   Base name for output files
    !! @param[in]     complexes     Structure containing complex text lines/headers
    subroutine linkage_clustering_from_matrix(matrix, n, linkage_type, output_name, complexes)

      IMPLICIT NONE

      real (kind=8), dimension(:, :), allocatable, intent(inout) :: matrix
      integer, intent(in) :: n
      character(len=*), intent(in) :: linkage_type
      character*128, intent(in) :: output_name
      type(type_assoc_file) :: complexes

        integer :: i, j, k, min_i, min_j, merge_counter, remaining_clusters, num_dist
        integer, dimension(n) :: cluster_size ! To track the size of each cluster
        real (kind=8):: min_dist, dist, mean_dist, standard_deviation, dist_threshold
        real (kind=8):: sum_dist, sum_sq_dist
        logical, dimension(n) :: active_points

        integer, dimension(n) :: cluster_count  ! Number of elements in each cluster
        integer, dimension(n) :: cluster_parent ! Parent cluster for each point - O(n) tracking

        integer, dimension(n) :: representative_indexes  ! Store most representative value index for each cluster
        real (kind=8), dimension(n) :: representative_values  ! Store representative values for each cluster
        real (kind=8), dimension(n) :: cluster_average, cluster_sd

        real (kind=8) :: local_min_dist, min_val, max_val, dispersion_range_coeff
        integer :: local_min_i, local_min_j

        real (kind=8) :: alpha
        real (kind=8) :: start_time, end_time, elapsed_time

        ! Initialize cluster tracking (each point starts alone).
        active_points = .true.
        cluster_size = 1
        remaining_clusters = n
        merge_counter = 0

        ! Initialize cluster parent - each point starts in its own cluster
        do i = 1, n
          cluster_parent(i) = i
          cluster_count(i) = 1
        end do

        ! Compute matrix stats to build the stopping threshold.
        sum_dist = 0.0
        sum_sq_dist = 0.0
        num_dist = 0

        min_val = huge(1.0)
        max_val = 0.0d0

        !$OMP PARALLEL DO REDUCTION(+:sum_dist, sum_sq_dist, num_dist) PRIVATE(i, j)
        do i = 1, n
          do j = i + 1, n
                sum_dist = sum_dist + matrix(i, j)
                sum_sq_dist = sum_sq_dist + matrix(i, j) ** 2
                num_dist = num_dist + 1
                if (matrix(i, j) < min_val) then
                    min_val = matrix(i, j)
                end if
                if (matrix(i, j) > max_val) then
                    max_val = matrix(i, j)
                end if
          end do
        end do
        !$OMP END PARALLEL DO

        mean_dist = sum_dist / num_dist
        standard_deviation = sqrt((sum_sq_dist / num_dist) - mean_dist ** 2)

        dispersion_range_coeff = (max_val - min_val)/(max_val + min_val)  ! Relative range of values in the array
      
        alpha = standard_deviation / sqrt(standard_deviation**2 + mean_dist**2)
        dist_threshold = (alpha*mean_dist - (1-alpha)*standard_deviation) / dispersion_range_coeff  ! Scale threshold by relative range to adapt to different value distributions
      
        ! ====================================================================
        ! CPU PARALLEL CLUSTERING - OpenMP implementation
        ! ====================================================================
        ! Iterative nearest-cluster merges with selected linkage type and OpenMP loops.
        print *, 'mod_clust_algorithm v1.0 | Starting clustering... For help, see documentation or use linkage_type="help".'
        start_time = OMP_GET_WTIME()

        ! Main clustering loop
        do while (remaining_clusters > 1)

            
                min_dist = huge(1.0)  ! Set to a very large number
                min_i = -1
                min_j = -1

                !$OMP PARALLEL PRIVATE(i, j, dist, local_min_dist, local_min_i, local_min_j) SHARED(min_dist, min_i, min_j)
                    ! Initialize local private variables at the start of each thread
                    local_min_dist = huge(1.0)  ! Set to a very large number
                    local_min_i = -1
                    local_min_j = -1

                    !$OMP DO
                    do i = 1, n
                        if (.not. active_points(i)) cycle
                        do j = i + 1, n
                            if (.not. active_points(j)) cycle
                            dist = matrix(i, j)
                            if (dist < local_min_dist) then
                                local_min_dist = dist
                                local_min_i = i
                                local_min_j = j
                            end if
                        end do
                    end do
                    !$OMP END DO
                
                    ! After the parallel loop ends, update the global minimum in a thread-safe manner
                    !$OMP CRITICAL
                        if (local_min_dist < min_dist) then
                            min_dist = local_min_dist
                            min_i = local_min_i
                            min_j = local_min_j
                        end if
                    !$OMP END CRITICAL
                !$OMP END PARALLEL
          
                ! Stop merging if clusters are too far apart
                if (min_dist > dist_threshold) exit
          
                merge_counter = merge_counter + 1
          
                ! Merge clusters min_i and min_j
                if (merge_counter <= n) then
                    ! Update cluster parent only (not cluster_indexes)
                    !$OMP PARALLEL DO
                    do k = 1, n
                        if (cluster_parent(k) == min_j) then
                            cluster_parent(k) = min_i
                        end if
                    end do
                    !$OMP END PARALLEL DO
                    cluster_count(min_i) = cluster_count(min_i) + cluster_count(min_j)
                end if
          
                active_points(min_j) = .false.
                remaining_clusters = remaining_clusters - 1
          
                ! Update distances using selected linkage method
                !$OMP PARALLEL DO REDUCTION(+:sum_dist, sum_sq_dist) PRIVATE(i)
                do i = 1, n
                    if (i == min_i .or. i == min_j .or. .not. active_points(i)) cycle

                    sum_dist = sum_dist - matrix(min_i, i) - matrix(min_j, i)
                    sum_sq_dist = sum_sq_dist - matrix(min_i, i) ** 2 - matrix(min_j, i) ** 2

                    ! Apply the selected linkage method
                    if (trim(linkage_type) == 'min') then
                        ! Minimum linkage (single linkage)
                        matrix(min_i, i) = min(matrix(min_i, i), matrix(min_j, i))
                    else if (trim(linkage_type) == 'max') then
                        ! Maximum linkage (complete linkage)
                        matrix(min_i, i) = max(matrix(min_i, i), matrix(min_j, i))
                    else
                        ! Mean linkage (average linkage) - default
                        matrix(min_i, i) = (matrix(min_i, i) * cluster_size(min_i) + &
                                            matrix(min_j, i) * cluster_size(min_j)) / &
                                            (cluster_size(min_i) + cluster_size(min_j))
                    end if
                    matrix(i, min_i) = matrix(min_i, i)

                end do
                !$OMP END PARALLEL DO

                ! Update cluster_size AFTER distance update (must use pre-merge sizes above)
                cluster_size(min_i) = cluster_size(min_i) + cluster_size(min_j)
          
                if (mod(remaining_clusters, 100) == 0) then
                    print *, '[mod_clust_algorithm] Remaining clusters:', remaining_clusters
                end if
            end do
            
        end_time = OMP_GET_WTIME()
        elapsed_time = end_time - start_time
        print *, '========================================'
        print *, '[mod_clust_algorithm] Clustering complete.'
        print *, '  CPU Time:', elapsed_time, 'seconds'
        print *, '========================================'
        
        ! Find the most representative value using cluster_parent (O(n) memory)
        ! Representative: member with min avg distance to all other members
        !$OMP PARALLEL DO PRIVATE(i, j, k, dist, min_dist)
        do i = 1, n
            if (.not. active_points(i)) cycle

            min_dist = 1.0e30
            representative_indexes(i) = i

            ! For each member j of cluster i
            do j = 1, n
                if (cluster_parent(j) /= i) cycle
                dist = 0.0d0
                do k = 1, n
                    if (cluster_parent(k) /= i) cycle
                    dist = dist + matrix(j, k)
                end do
                dist = dist / cluster_count(i)

                if (dist < min_dist) then
                    min_dist = dist
                    representative_indexes(i) = j
                    representative_values(i) = min_dist
                end if
            end do
        end do
        !$OMP END PARALLEL DO

        ! Compute mean and standard deviation
        !$OMP PARALLEL DO PRIVATE(i, j, mean_dist, standard_deviation)
        do i = 1, n
            if (.not. active_points(i)) cycle
            mean_dist = 0.0d0
            do j = 1, n
                if (cluster_parent(j) /= i) cycle
                mean_dist = mean_dist + matrix(i, j)
            end do
            mean_dist = mean_dist / cluster_count(i)
            cluster_average(i) = mean_dist

            standard_deviation = 0.0d0
            do j = 1, n
                if (cluster_parent(j) /= i) cycle
                standard_deviation = standard_deviation + &
                    (cluster_average(i) - matrix(i, j))**2
            end do
            standard_deviation = sqrt(standard_deviation / cluster_count(i))
            cluster_sd(i) = standard_deviation
        end do
        !$OMP END PARALLEL DO

        call write_cluster_elements(n, cluster_parent, cluster_count, active_points, output_name)
        call write_clust_info(n, representative_indexes, active_points, cluster_count, cluster_average, cluster_sd, output_name)
        call write_cluster_complexes(n, cluster_parent, cluster_count, active_points, output_name, complexes)
          
          
        print *, '[mod_clust_algorithm] Clustering complete. See output files for results.'
    
    end subroutine linkage_clustering_from_matrix

    !> Perform hierarchical clustering on a distance array.
    !!
    !! This routine performs hierarchical clustering using minimum, maximum,
    !! or mean linkage on the input `array` (shape `(n)`) of
    !! pairwise distances or similarity scores. The routine mutates `array`
    !! during merging and terminates merging when clusters are farther apart
    !! than an internally computed threshold. Results are written to files
    !! named using `output_name` and complexes information is optionally
    !! recorded via the `complexes` structure.
    !!
    !! @param[in,out] array         Array of per-encounter values (n)
    !!                              This array is updated during clustering
    !! @param[in]     n             Number of elements / encounters
    !! @param[in]     linkage_type  Type of linkage: 'min', 'max', or 'mean'
    !! @param[in]     output_name   Base name for output files
    !! @param[in]     complexes     Structure containing complex text lines/headers
    subroutine linkage_clustering_from_array(array, n, linkage_type, output_name, complexes)

      IMPLICIT NONE

      real (kind=8), dimension(:), allocatable, intent(inout) :: array
      integer, intent(in) :: n
      character(len=*), intent(in) :: linkage_type
      character*128, intent(in) :: output_name
      type(type_assoc_file) :: complexes

      integer :: i, j, k, min_i, min_j, merge_counter, remaining_clusters, num_dist
      integer, dimension(n) :: cluster_size ! To track the size of each cluster
      real (kind=8):: min_dist, dist, mean_dist, standard_deviation, dist_threshold
      real (kind=8):: sum_dist, sum_sq_dist, old_val, dispersion_range_coeff, middle_range
      logical, dimension(n) :: active_points

      integer, dimension(n) :: cluster_count  ! Number of elements in each cluster
      integer, dimension(n) :: cluster_parent ! Parent cluster for each point - O(n) tracking

      integer, dimension(n) :: representative_indexes  ! Store most representative value index for each cluster
      real (kind=8), dimension(n) :: representative_values  ! Store representative values for each cluster
      real (kind=8), dimension(n) :: cluster_average, cluster_sd

      real (kind=8) :: local_min_dist
      integer :: local_min_i, local_min_j

      real (kind=8) :: alpha, factor
      real (kind=8) :: start_time, end_time, elapsed_time

      active_points = .true.
      cluster_size = 1
      remaining_clusters = n
      merge_counter = 0

      ! Initialize cluster parent - each point starts in its own cluster
      do i = 1, n
        cluster_parent(i) = i
        cluster_count(i) = 1
      end do

      ! Initialize statistics
      sum_dist = 0.0
      sum_sq_dist = 0.0
      num_dist = 0

      !$OMP PARALLEL DO REDUCTION(+:sum_dist, sum_sq_dist, num_dist) PRIVATE(i, j, dist)
      do i = 1, n
        do j = i + 1, n
            dist = array(j) - array(i)
            sum_dist = sum_dist + dist
            sum_sq_dist = sum_sq_dist + dist ** 2
            num_dist = num_dist + 1
        end do
      end do
      !$OMP END PARALLEL DO

      mean_dist = sum_dist / num_dist
      standard_deviation = sqrt((sum_sq_dist / num_dist) - mean_dist ** 2)
      middle_range = (array(n) + array(1)) / 2.0d0
      dispersion_range_coeff = (array(n) - array(1))/(array(n) + array(1))  ! Relative range of values in the array
      
      factor = 1 - abs(middle_range - (mean_dist + 2*standard_deviation)) / (middle_range + mean_dist + 2*standard_deviation)  ! Scale factor based on how close mean is to middle of range

      alpha = standard_deviation / sqrt(standard_deviation**2 + mean_dist**2)
    !   dist_threshold = (alpha*mean_dist - (1-alpha)*standard_deviation) / factor  ! Scale threshold by relative range to adapt to different value distributions
      dist_threshold = (alpha*mean_dist - (1-alpha)*standard_deviation) / dispersion_range_coeff  ! Scale threshold by relative range to adapt to different value distributions
    !   print *, 'mean_dist:', mean_dist
    !   print *, 'standard_deviation:', standard_deviation
    !   print *, 'middle_range:', middle_range
    !   print *, 'dispersion_range_coeff:', dispersion_range_coeff
    !   print *, 'factor:', factor
    !   print *, 'alpha:', alpha
    !   print *, 'dist_threshold:', dist_threshold
      
    !   STOP

      ! ====================================================================
      ! OPENMP PARALLEL CLUSTERING for sorted 1D array
      ! ====================================================================
    print *, 'mod_clust_algorithm v1.2 | Starting clustering ', &
             '(1D array)... For help, see documentation or use ', &
             'linkage_type="help".'
      start_time = OMP_GET_WTIME()

      ! Main clustering loop
      do while (remaining_clusters > 1)

          min_dist = huge(1.0d0)
          min_i = -1
          min_j = -1

          !$OMP PARALLEL PRIVATE(i, j, dist, local_min_dist, &
          !$OMP& local_min_i, local_min_j) SHARED(min_dist, min_i, min_j)
              local_min_dist = huge(1.0d0)
              local_min_i = -1
              local_min_j = -1

              !$OMP DO
              do i = 1, n
                  if (.not. active_points(i)) cycle
                  ! Find next active point after i
                  j = i + 1
                  do while (j <= n)
                      if (active_points(j)) exit
                      j = j + 1
                  end do
                  if (j > n) cycle
                  dist = array(j) - array(i)
                  if (dist < local_min_dist) then
                      local_min_dist = dist
                      local_min_i = i
                      local_min_j = j
                  end if
              end do
              !$OMP END DO

              !$OMP CRITICAL
                  if (local_min_dist < min_dist) then
                      min_dist = local_min_dist
                      min_i = local_min_i
                      min_j = local_min_j
                  end if
              !$OMP END CRITICAL
          !$OMP END PARALLEL

          ! Stop merging if clusters are too far apart
          if (min_dist > dist_threshold) exit

          merge_counter = merge_counter + 1

          ! Merge clusters min_i and min_j
          if (merge_counter <= n) then
              !$OMP PARALLEL DO
              do k = 1, n
                  if (cluster_parent(k) == min_j) then
                      cluster_parent(k) = min_i
                  end if
              end do
              !$OMP END PARALLEL DO
              cluster_count(min_i) = cluster_count(min_i) + cluster_count(min_j)
          end if

          active_points(min_j) = .false.
          remaining_clusters = remaining_clusters - 1

          ! Save old value before update
        !   old_val = array(min_i)

          cluster_average(min_i) = (array(min_i) * cluster_size(min_i) + array(min_j) * cluster_size(min_j)) / &
                          (cluster_size(min_i) + cluster_size(min_j))

          ! Update array value based on linkage type
          if (trim(linkage_type) == 'min') then
              ! Single linkage: keep minimum value (already the smaller)
              array(min_i) = array(min_i)
          else if (trim(linkage_type) == 'max') then
              ! Complete linkage: keep maximum value
              array(min_i) = array(min_j)
          else
              ! Average linkage: weighted mean
              array(min_i) = cluster_average(min_i)
          end if

          ! Update cluster size after array update (keeps weights correct)
          cluster_size(min_i) = cluster_size(min_i) + cluster_size(min_j)

          ! Re-sort: bubble the merged element rightward to maintain sorted order.
          ! The merged value can only stay the same or increase:
          !   min linkage -> keeps smaller value (no movement)
          !   max linkage -> takes larger value (moves right)
          !   avg linkage -> weighted mean between the two (moves right)
          j = min_i
          do while (.true.)
              ! Find the next active point after j
              i = j + 1
              do while (i <= n)
                  if (active_points(i)) exit
                  i = i + 1
              end do
              if (i > n) exit
              if (array(j) <= array(i)) exit
              ! Swap all data at positions j and i
              old_val = array(j);  array(j) = array(i);  array(i) = old_val
              k = cluster_size(j);  cluster_size(j) = cluster_size(i);  cluster_size(i) = k
              k = cluster_count(j);  cluster_count(j) = cluster_count(i);  cluster_count(i) = k
              ! Remap cluster_parent: all references j<->i
              !$OMP PARALLEL DO PRIVATE(k)
              do k = 1, n
                  if (cluster_parent(k) == j) then
                      cluster_parent(k) = i
                  else if (cluster_parent(k) == i) then
                      cluster_parent(k) = j
                  end if
              end do
              !$OMP END PARALLEL DO
              j = i
          end do

                    if (mod(remaining_clusters, 100) == 0) then
                        print *, '[mod_clust_algorithm] Remaining clusters:', remaining_clusters
                    end if
      end do

      end_time = OMP_GET_WTIME()
      elapsed_time = end_time - start_time
    print *, '========================================'
    print *, '[mod_clust_algorithm] Clustering complete (1D array).'
    print *, '  CPU Time:', elapsed_time, 'seconds'
    print *, '========================================'

      ! Find the most representative value: member closest to cluster mean
      ! Uses cluster_parent directly — O(n) memory
      !$OMP PARALLEL DO PRIVATE(i, j, dist, min_dist, mean_dist)
      do i = 1, n
          if (.not. active_points(i)) cycle

          ! Compute cluster mean from array
          mean_dist = 0.0d0
          do j = 1, n
              if (cluster_parent(j) /= i) cycle
              mean_dist = mean_dist + array(j)
          end do
          mean_dist = mean_dist / cluster_count(i)

          ! Find member closest to mean
          min_dist = huge(1.0d0)
          representative_indexes(i) = i
          do j = 1, n
              if (cluster_parent(j) /= i) cycle
              dist = abs(array(j) - mean_dist)
              if (dist < min_dist) then
                  min_dist = dist
                  representative_indexes(i) = j
                  representative_values(i) = array(j)
              end if
          end do
      end do
      !$OMP END PARALLEL DO

      ! Compute standard deviation per cluster
      !$OMP PARALLEL DO PRIVATE(i, j, mean_dist, standard_deviation)
      do i = 1, n
          standard_deviation = 0.0d0
          do j = 1, n
              if (cluster_parent(j) /= i) cycle
              standard_deviation = standard_deviation + &
                  (cluster_average(i) - array(j))**2
          end do
          standard_deviation = sqrt(standard_deviation / cluster_count(i))
          cluster_sd(i) = standard_deviation
      end do
      !$OMP END PARALLEL DO

      call write_cluster_elements(n, cluster_parent, cluster_count, active_points, output_name)
      call write_clust_info(n, representative_indexes, active_points, cluster_count, cluster_average, cluster_sd, output_name)
      call write_cluster_complexes(n, cluster_parent, cluster_count, active_points, output_name, complexes)

    print *, '[mod_clust_algorithm] Clustering complete. See output files for results.'

    end subroutine linkage_clustering_from_array

    !> Sort clusters and associated arrays by a given key array.
    !!
    !! Reorders `cluster_indexes` and all supplied per-cluster arrays so
    !! that they are sorted by ascending values in `array`. This is a
    !! stable reordering helper used prior to writing results.
    !!
    !! @param[in]     tot_encounters        Number of clusters/encounters
    !! @param[in,out] cluster_indexes       Matrix storing cluster member indexes
    !! @param[in,out] cluster_count         Number of members in each cluster
    !! @param[in,out] active_clusters       Logical flags marking active clusters
    !! @param[in,out] representative_indexes Representative index per cluster
    !! @param[in,out] cluster_average       Per-cluster average values
    !! @param[in,out] cluster_sd            Per-cluster standard deviation
    !! @param[in,out] array                 Key array used to sort clusters
    subroutine sort_complexes(tot_encounters, cluster_indexes, cluster_count, active_clusters, &
      representative_indexes, cluster_average, cluster_sd, array)

      integer, intent(in) :: tot_encounters
      integer, dimension(:, :), intent(inout) :: cluster_indexes ! Store indexes of each cluster
      integer, dimension(:), intent(inout) :: cluster_count  ! Number of elements in each cluster
      logical, dimension(:), intent(inout) :: active_clusters
      integer, dimension(:), intent(inout) :: representative_indexes
      real (kind=8), dimension(:), intent(inout) :: cluster_average, cluster_sd
      real (kind=8), dimension(:), intent(inout) :: array
      
      integer :: i, j, sorted_i, encounters_found, temp_j
      real (kind=8) :: temp
      integer, dimension(tot_encounters) :: sorted_indexes
      integer, dimension(tot_encounters, tot_encounters) :: temp_cluster_indexes
      integer, dimension(tot_encounters) :: temp_representative_indexes, temp_cluster_count
      real (kind=8), dimension(tot_encounters) :: temp_array
      real (kind=8), dimension(tot_encounters) :: temp_cluster_average, temp_cluster_sd
      logical, dimension(tot_encounters) :: temp_active_clusters

      temp_array = array
      temp_cluster_indexes = cluster_indexes
      temp_cluster_average = cluster_average
      temp_cluster_sd = cluster_sd
      temp_representative_indexes = representative_indexes

      ! Fill array indexes

      do i = 1, tot_encounters
        sorted_indexes(i) = i
      end do

      ! Bubble sort
      do i = 1, tot_encounters-1
        do j = 1, tot_encounters-i
            if (temp_array(j) > temp_array(j+1)) then
                temp = temp_array(j)
                temp_array(j) = temp_array(j+1)
                temp_array(j+1) = temp

                temp_j = sorted_indexes(j)
                sorted_indexes(j) = sorted_indexes(j+1)
                sorted_indexes(j+1) = temp_j
            end if
        end do
      end do

      do i = 1, tot_encounters
        sorted_i = sorted_indexes(i)
    
        ! Cache values for better memory access efficiency
        temp_cluster_count(i) = cluster_count(sorted_i)
        temp_active_clusters(i) = active_clusters(sorted_i)
        temp_cluster_average(i) = cluster_average(sorted_i)
        temp_cluster_sd(i) = cluster_sd(sorted_i)
        temp_representative_indexes(i) = representative_indexes(sorted_i)
    
        encounters_found = 0
        do j = 1, tot_encounters
            do temp_j = 1, cluster_count(sorted_i)
                if (cluster_indexes(sorted_i, temp_j) == sorted_indexes(j)) then
                    encounters_found = encounters_found + 1
                    temp_cluster_indexes(i, encounters_found) = sorted_indexes(j)
    
                    ! Exit early if all encounters are found
                    if (encounters_found == cluster_count(sorted_i)) exit
                end if
            end do
        end do
      end do

      
      array = temp_array
      cluster_indexes = temp_cluster_indexes
      cluster_average = temp_cluster_average
      cluster_sd = temp_cluster_sd
      representative_indexes = temp_representative_indexes
      active_clusters = temp_active_clusters
      cluster_count = temp_cluster_count
      
      
    end subroutine sort_complexes

    !> Write cluster membership lists to a text file.
    !!
    !! Outputs a file named `<output_name>_clusters.txt` containing the
    !! members of each active cluster in a matrix-like format.
    !!
    !! @param[in] tot_encounters   Number of possible clusters/encounters
    !! @param[in] cluster_parent   Parent cluster for each encounter
    !! @param[in] cluster_count    Number of members per cluster
    !! @param[in] active_clusters  Logical mask indicating active clusters
    !! @param[in] output_name      Base filename for output
    subroutine write_cluster_elements(tot_encounters, cluster_parent, cluster_count, active_clusters, output_name)

      IMPLICIT NONE
      integer, intent(in) :: tot_encounters
      integer, dimension(:), intent(in) :: cluster_parent
      integer, dimension(:), intent(in) :: cluster_count
      logical, dimension(:), intent(in) :: active_clusters
      character*128, intent(in) :: output_name
      integer :: i, j, unit_number, clust_number
      character(len=13):: clust_str
      character(len=4):: clust_n_str
      character*128 :: filename

            if (size(cluster_count) /= tot_encounters) then
                write(*,*) '[mod_clust_algorithm] WARNING: cluster_count size mismatch in write_cluster_elements.'
            end if

      ! Open file to write the final cluster content
      unit_number = 20
      filename = trim(adjustl(output_name)) // "_clusters.txt"
      open(unit=unit_number, file=filename, status='replace', action='write')

      ! Write the content of each active cluster
      clust_number = 1
      do i = 1, tot_encounters
        if (active_clusters(i)) then
          clust_str = ""
          write(clust_n_str, '(I4)') clust_number
          clust_str = "Cluster" // adjustr(clust_n_str) // ":"
          write(unit_number, "(A13)", advance='no') clust_str
          do j = 1, tot_encounters
            if (cluster_parent(j) == i) then
              write(unit_number, "(I5,1X)", advance='no') j
            end if
          end do
          write(unit_number, *)  ! newline
          clust_number = clust_number + 1
        end if
      end do

      ! Close the file
      close(unit_number)

    end subroutine write_cluster_elements

    !> Write summary information about clusters to a text file.
    !!
    !! Produces a summary file `<output_name>_clust_info.txt`
    !! listing cluster population, representative index, average and
    !! standard deviation for each active cluster.
    !!
    !! @param[in] tot_encounters        Total number of encounters
    !! @param[in] representative_indexes Representative index per cluster
    !! @param[in] active_clusters       Logical mask indicating active clusters
    !! @param[in] clust_count           Array with population of each cluster
    !! @param[in] cluster_average       Per-cluster average values
    !! @param[in] cluster_sd            Per-cluster standard deviations
    !! @param[in] output_name           Base filename for outputs
    subroutine write_clust_info(tot_encounters, representative_indexes, active_clusters, &
      clust_count, cluster_average, cluster_sd, output_name)

      IMPLICIT NONE
      integer, intent(in) :: tot_encounters
      integer, dimension(:), intent(in) :: representative_indexes
      logical, dimension(:), intent(in) :: active_clusters
      integer, dimension(tot_encounters), intent(in) :: clust_count
      real (kind=8), dimension(tot_encounters), intent(in) :: cluster_average, cluster_sd
      character*128, intent(in) :: output_name

      integer :: i, unit_number, clust_number
      character(len=13):: clust_str
      character(len=4):: clust_n_str
      character*128 :: filename

      ! Open file to write the final cluster content
      unit_number = 20
      filename = trim(adjustl(output_name)) // "_clust_info.txt"
      open(unit=unit_number, file=filename, status='replace', action='write')
      
      write(unit_number, "(1(A26,1X)1(A5,1X),2(A9,1X))") "Clust pop.", "Repr.", "Average", "SD"

      ! Write the content of each active cluster into the file in matrix format
      clust_number = 1
      do i = 1, tot_encounters
        if (active_clusters(i)) then
          clust_str = ""
          write(clust_n_str, '(I4)') clust_number
          clust_str = "Cluster" // adjustr(clust_n_str) // ":"
          write(unit_number, "(1(A12,3X)1(I10,1X)1(I5,1X)2(F10.3,1X))") clust_str, clust_count(i), &
          representative_indexes(i), &
          cluster_average(i), cluster_sd(i)
          clust_number = clust_number + 1
        end if
      end do

      ! Close the file
      close(unit_number)

    end subroutine write_clust_info

    !> Write complexes contained in each cluster to individual files.
    !!
    !! For each active cluster this routine creates a file
    !! `<output_name>_clusterXXXX_complexes` containing the header lines
    !! from `complexes%head` and the encounter lines for the cluster
    !! members taken from `complexes%lines`.
    !!
    !! @param[in] tot_encounters   Number of clusters/encounters
    !! @param[in] cluster_parent   Parent cluster for each encounter
    !! @param[in] cluster_count    Number of members per cluster
    !! @param[in] active_clusters  Logical mask indicating active clusters
    !! @param[in] output_name      Base filename for output files
    !! @param[in] complexes        Structure with header and lines to write
    subroutine write_cluster_complexes(tot_encounters, cluster_parent, cluster_count, active_clusters, &
                                output_name, complexes)

      IMPLICIT NONE
      integer, intent(in) :: tot_encounters
      integer, dimension(:), intent(in) :: cluster_parent
      integer, dimension(:), intent(in) :: cluster_count
      logical, dimension(:), intent(in) :: active_clusters
      character*128, intent(in) :: output_name
      type(type_assoc_file) :: complexes

      integer :: i, j, unit_number, clust_number
      character(len=4):: str_clust_number
      character(len=200) :: filename

            if (size(cluster_count) /= tot_encounters) then
                write(*,*) '[mod_clust_algorithm] WARNING: cluster_count size mismatch in write_cluster_complexes.'
            end if

      unit_number = 21
      clust_number = 1
      do i = 1, tot_encounters
        if (active_clusters(i)) then
          write(str_clust_number, '(I4.4)') clust_number
          filename = trim(adjustl(output_name)) // "_cluster" // str_clust_number // "_complexes"
          open(unit=unit_number, file=filename, status='replace', action='write')
          do j = 1, 4
            write(unit_number, "(A)") complexes % head(j)
          end do

          do j = 1, tot_encounters
            if (cluster_parent(j) == i) then
              write(unit_number, "(A)") complexes % lines(j)
            end if
          end do
          close(unit_number)
          clust_number = clust_number + 1
        end if
      end do

    end subroutine write_cluster_complexes
    !subroutine write_clust_info(min_i, min_j)

end module mod_clust_algorithm