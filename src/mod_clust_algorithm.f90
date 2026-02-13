module mod_clust_algorithm

  !> Clustering algorithms and helpers for processing encounter matrices.
  !!
  !! This module contains routines to perform hierarchical clustering
  !! (minimum, maximum, or mean linkage) on a distance/similarity matrix,
  !! helpers to sort cluster results, and routines to write cluster outputs
  !! and summaries to text files. Parallel regions use OpenMP for some
  !! heavy loops.
  !!
  !! Note: Several routines expect arrays sized to `n` or matrices sized to
  !! `(n,n)` where `n` is the number of encounters/complexes.
  !!
  !! See individual subroutines for parameter descriptions.
  !!
  !! @author Abraham Muñiz-Chicharro

  USE read_input
  USE OMP_LIB
  USE mod_cuda, ONLY: cuda_complete_clustering

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
    !! @param[in,out] matrix       Pairwise distance/similarity matrix (n,n).
    !!                              This matrix is updated during clustering.
    !! @param[in]     n            Number of elements / encounters (integer).
    !! @param[in]     linkage_type Type of linkage: 'min', 'max', or 'mean' (string).
    !! @param[in]     output_name  Base name for output files (string).
    !! @param[in]     complexes    Structure containing complex text lines/headers.
    !! @param[in]     opt_array    Optional array of per-encounter values used
    !!                              to compute cluster averages instead of
    !!                              reading values from `matrix`.
    subroutine linkage_clustering(matrix, n, linkage_type, output_name, complexes, opt_array, use_cuda)

      IMPLICIT NONE

      real (kind=8), dimension(:, :), intent(inout) :: matrix
      real (kind=8), dimension(:), intent(inout), optional :: opt_array
      character*128, intent(in) :: output_name
      character(len=*), intent(in) :: linkage_type
      integer, intent(in) :: n
      logical, intent(in), optional :: use_cuda
      type(type_assoc_file) :: complexes

        integer :: i, j, k, min_i, min_j, merge_counter, remaining_clusters, num_dist
        integer, dimension(n) :: cluster_size ! To track the size of each cluster
        real (kind=8):: min_dist, dist, mean_dist, standard_deviation, dist_threshold
        real (kind=8):: sum_dist, sum_sq_dist
        logical, dimension(n) :: active_points

        integer, dimension(n, n) :: cluster_indexes ! Store indexes of each cluster
        integer, dimension(n) :: cluster_count  ! Number of elements in each cluster

        integer, dimension(n) :: representative_indexes  ! Store most representative value index for each cluster
        real (kind=8), dimension(n) :: representative_values  ! Store representative values for each cluster
        real (kind=8), dimension(n) :: cluster_average, cluster_sd

        real (kind=8) :: local_min_dist
        integer :: local_min_i, local_min_j

        real (kind=8) :: alpha
        logical :: use_cuda_accel
        real (kind=8) :: start_time, end_time, elapsed_time

    
        ! Handle optional use_cuda parameter
        if (present(use_cuda)) then
            use_cuda_accel = use_cuda
        else
            use_cuda_accel = .false.
        end if

        active_points = .true.
        cluster_size = 1
        remaining_clusters = n
        merge_counter = 0

        ! Initialize cluster indexes with the initial clusters
        do i = 1, n
          cluster_indexes(i, 1) = i  ! Each initial cluster has itself as the only member
          cluster_count(i) = 1       ! Initially, each cluster has one member
        end do

        ! Initialize statistics
        sum_dist = 0.0
        sum_sq_dist = 0.0
        num_dist = 0

        !$OMP PARALLEL DO REDUCTION(+:sum_dist, sum_sq_dist, num_dist) PRIVATE(i, j)
        do i = 1, n
          do j = i + 1, n
              sum_dist = sum_dist + matrix(i, j)
              sum_sq_dist = sum_sq_dist + matrix(i, j) ** 2
              num_dist = num_dist + 1
          end do
        end do
        !$OMP END PARALLEL DO

        mean_dist = sum_dist / num_dist
        standard_deviation = sqrt((sum_sq_dist / num_dist) - mean_dist ** 2)

        alpha = standard_deviation / sqrt(standard_deviation**2 + mean_dist**2)
        dist_threshold = alpha*mean_dist - (1-alpha)*standard_deviation
      
        ! ====================================================================
        ! COMPLETE CUDA PARALLELIZATION - All clustering done on GPU
        ! ====================================================================
        if (use_cuda_accel) then
            print *, 'Starting GPU clustering (complete parallelization)...'
            start_time = OMP_GET_WTIME()
            
            ! Call complete GPU clustering - all iterations happen on GPU
            if (cuda_complete_clustering(matrix, n, dist_threshold, linkage_type, &
                                         cluster_indexes, cluster_count, active_points, cluster_size) /= 0) then
                print *, "ERROR: CUDA complete clustering failed"
                return
            end if
            
            ! GPU clustering is complete - active_points, cluster_count, and cluster_indexes are updated
            ! Note: cluster_size updates need to be reconstructed from merges
            ! For now, we skip updating clustering statistics as they're done on GPU
            remaining_clusters = count(active_points)
            
            end_time = OMP_GET_WTIME()
            elapsed_time = end_time - start_time
            print *, '========================================'
            print *, 'GPU clustering complete.'
            print *, 'Remaining clusters:', remaining_clusters
            print *, 'GPU Time:', elapsed_time, 'seconds'
            print *, '========================================'
            
        else
            ! ====================================================================
            ! CPU PARALLEL CLUSTERING - Original OpenMP implementation
            ! ====================================================================
            print *, 'Starting CPU clustering...'
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
                    !$OMP PARALLEL DO
                    do k = 1, cluster_count(min_j)
                        cluster_indexes(min_i, cluster_count(min_i) + k) = cluster_indexes(min_j, k)
                    end do
                    !$OMP END PARALLEL DO
              
                    cluster_count(min_i) = cluster_count(min_i) + cluster_count(min_j)
                end if
          
                active_points(min_j) = .false.
                cluster_size(min_i) = cluster_size(min_i) + cluster_size(min_j)
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

                    sum_dist = sum_dist + matrix(min_i, i)
                    sum_sq_dist = sum_sq_dist + matrix(min_i, i) ** 2

                end do
                !$OMP END PARALLEL DO
                sum_dist = sum_dist - min_dist
                sum_sq_dist = sum_sq_dist - min_dist ** 2
                num_dist = num_dist - 1
          
                if (mod(remaining_clusters, 100) == 0) then
                    print *, 'Remaining clusters:', remaining_clusters
                end if
            end do
            
            end_time = OMP_GET_WTIME()
            elapsed_time = end_time - start_time
            print *, '========================================'
            print *, 'CPU clustering complete.'
            print *, 'CPU Time:', elapsed_time, 'seconds'
            print *, '========================================'
        end if  ! End of if (use_cuda_accel) else block

        ! Find the most representative value
        !$OMP PARALLEL DO PRIVATE(i, j, k, dist, min_dist)
        do i = 1, n
            if (.not. active_points(i)) cycle
        
            min_dist = 1.0e30
            representative_indexes(i) = -1
        
            do j = 1, cluster_count(i)
                dist = 0.0
                do k = 1, cluster_count(i)
                    dist = dist + matrix(cluster_indexes(i, j), cluster_indexes(i, k))
                end do
                dist = dist / cluster_count(i)
              
                if (dist < min_dist) then
                    min_dist = dist
                    representative_indexes(i) = cluster_indexes(i, j)
                    representative_values(i) = min_dist
                end if
            end do
        end do
        !$OMP END PARALLEL DO
      
        ! Compute mean and standard deviation
        !$OMP PARALLEL DO PRIVATE(i, j, mean_dist, standard_deviation)
        do i = 1, n
            if (.not. active_points(i)) cycle
            mean_dist = 0.0
            do j = 1, cluster_count(i)
                if (present(opt_array)) then
                    mean_dist = mean_dist + opt_array(cluster_indexes(i, j))
                else
                    mean_dist = mean_dist + matrix(i, cluster_indexes(i, j))
                end if
            end do
            mean_dist = mean_dist / cluster_count(i)
            cluster_average(i) = mean_dist
          
            standard_deviation = 0.0
            do j = 1, cluster_count(i)
                if (present(opt_array)) then
                    standard_deviation = standard_deviation + (cluster_average(i) - &
                    opt_array(cluster_indexes(i, j)))**2
                else
                    standard_deviation = standard_deviation + (cluster_average(i) - &
                    matrix(i, cluster_indexes(i, j)))**2
                end if
            end do
            standard_deviation = sqrt(standard_deviation / cluster_count(i))
            cluster_sd(i) = standard_deviation
        end do
        !$OMP END PARALLEL DO
    
        !write(*,*) cluster_indexes(45, :)

        ! Optional sorting of clusters if opt_array is provided
        if ( present(opt_array) ) then
          call sort_complexes(n, cluster_indexes, cluster_count, active_points, &
          representative_indexes, cluster_average, cluster_sd, opt_array)
        end if
      
        call write_cluster_elements(n, cluster_indexes, cluster_count, active_points, output_name)
        call write_clust_info(n, representative_indexes, active_points, cluster_count, cluster_average, cluster_sd, output_name)
        call write_complexes(n, cluster_indexes, cluster_count, active_points, output_name, complexes)
          
          
        print *, 'Clustering complete.'
    
    end subroutine linkage_clustering

    !> Sort clusters and associated arrays by a given key array.
    !!
    !! Reorders `cluster_indexes` and all supplied per-cluster arrays so
    !! that they are sorted by ascending values in `array`. This is a
    !! stable reordering helper used prior to writing results.
    !!
    !! @param[in]     tot_encounters       Number of clusters/encounters
    !! @param[in,out] cluster_indexes      Matrix storing cluster member indexes
    !! @param[in,out] cluster_count        Number of members in each cluster
    !! @param[in,out] active_clusters      Logical flags marking active clusters
    !! @param[in,out] representative_indexes Representative index per cluster
    !! @param[in,out] cluster_average      Per-cluster average values
    !! @param[in,out] cluster_sd           Per-cluster standard deviation
    !! @param[in,out] array                Key array used to sort clusters
    subroutine sort_complexes(tot_encounters, cluster_indexes, cluster_count, active_clusters, &
      representative_indexes, cluster_average, cluster_sd, array)

      real (kind=8), dimension(:), intent(inout) :: array
      integer, intent(in) :: tot_encounters
      integer, dimension(:, :), intent(inout) :: cluster_indexes ! Store indexes of each cluster
      integer, dimension(:), intent(inout) :: cluster_count  ! Number of elements in each cluster
      integer, dimension(:), intent(inout) :: representative_indexes
      real (kind=8), dimension(:), intent(inout) :: cluster_average, cluster_sd
      logical, dimension(:), intent(inout) :: active_clusters
      
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
    !! @param[in] tot_encounters Number of possible clusters/encounters
    !! @param[in] cluster_indexes Matrix of cluster member indexes
    !! @param[in] cluster_count   Number of members per cluster
    !! @param[in] active_clusters Logical mask indicating active clusters
    !! @param[in] output_name     Base filename for output
    subroutine write_cluster_elements(tot_encounters, cluster_indexes, cluster_count, active_clusters, output_name)

      IMPLICIT NONE
      integer, intent(in) :: tot_encounters
      integer, dimension(:, :), intent(in) :: cluster_indexes
      logical, dimension(:), intent(in) :: active_clusters
      integer, dimension(:), intent(in) :: cluster_count
      character*128, intent(in) :: output_name
      integer :: i, j, unit_number, clust_number
      character(len=13):: clust_str
      character(len=4):: clust_n_str
      character*128 :: filename

      ! Open file to write the final cluster content
      unit_number = 20
      filename = trim(adjustl(output_name)) // "_clusters.txt"
      open(unit=unit_number, file=filename, status='replace', action='write')
      
      ! Write the content of each active cluster into the file in matrix format
      clust_number = 1
      do i = 1, tot_encounters
        if (active_clusters(i)) then
          clust_str = ""
          write(clust_n_str, '(I4)') clust_number
          clust_str = "Cluster" // adjustr(clust_n_str) // ":"
          write(unit_number, "(1(A13)10000(I5,1X))") clust_str, (cluster_indexes(i, j), j = 1, cluster_count(i))
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
    !! @param[in] tot_encounters     Total number of encounters
    !! @param[in] clust_count        Array with population of each cluster
    !! @param[in] active_clusters    Logical mask indicating active clusters
    !! @param[in] representative_indexes Representative index per cluster
    !! @param[in] output_name        Base filename for outputs
    !! @param[in] cluster_average    Per-cluster average values
    !! @param[in] cluster_sd         Per-cluster standard deviations
    subroutine write_clust_info(tot_encounters, representative_indexes, active_clusters, &
      clust_count, cluster_average, cluster_sd, output_name)

      IMPLICIT NONE
      integer, intent(in) :: tot_encounters
      integer, dimension(tot_encounters), intent(in) :: clust_count
      logical, dimension(:), intent(in) :: active_clusters
      integer, dimension(:), intent(in) :: representative_indexes
      character*128, intent(in) :: output_name
      integer :: i, unit_number, clust_number
      character(len=13):: clust_str
      character(len=4):: clust_n_str
      character*128 :: filename
      real (kind=8), dimension(tot_encounters) :: cluster_average, cluster_sd

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
    !! @param[in] tot_encounters Number of clusters/encounters
    !! @param[in] cluster_indexes Matrix of cluster member indexes
    !! @param[in] cluster_count   Number of members per cluster
    !! @param[in] active_clusters Logical mask indicating active clusters
    !! @param[in] output_name     Base filename for output files
    !! @param[in] complexes       Structure with header and lines to write
    subroutine write_complexes(tot_encounters, cluster_indexes, cluster_count, active_clusters, &
                                output_name, complexes)

      

      IMPLICIT NONE
      integer, intent(in) :: tot_encounters
      integer, dimension(:, :), intent(in) :: cluster_indexes
      logical, dimension(:), intent(in) :: active_clusters
      integer, dimension(:), intent(in) :: cluster_count
      character*128, intent(in) :: output_name
      integer :: i, j, unit_number, clust_number, encounter_index
      character(len=4):: str_clust_number
      character(len=200) :: filename
      type(type_assoc_file) :: complexes

      
      ! Write the content of each active cluster into the file in matrix format
      clust_number = 1
      do i = 1, tot_encounters
        if (active_clusters(i)) then
          write(str_clust_number, '(I4.4)') clust_number
          !write(*,*) str_clust_number
          filename = trim(adjustl(output_name)) // "_cluster" // str_clust_number // "_complexes"
          !write(*,*) filename
          open(unit=unit_number, file=filename, status='replace', action='write')
          do j = 1, 4
            write(unit_number, "(A)") complexes % head(j) 
          end do

          do j = 1, cluster_count(i)
            encounter_index = cluster_indexes(i, j)
            write(unit_number, "(A)") complexes % lines(encounter_index)
          end do
          close(unit_number)
          clust_number = clust_number + 1
        end if
      end do

      ! Close the file
      !close(unit_number)

    end subroutine write_complexes
    !subroutine write_clust_info(min_i, min_j)

    
end module mod_clust_algorithm