!> \file mod_clust_algorithm.f90
!! \brief Hierarchical clustering algorithms and output helpers.
!!
module mod_clust_algorithm

  !> Module: mod_clust_algorithm
  !!
  !! Clustering algorithms and helpers for processing encounter arrays.
  !!
  !! Provides hierarchical clustering (mean linkage) for sorted encounter arrays
  !! and output routines. Supports OpenMP parallelism.
  !!
  !! @author Abraham Muñiz-Chicharro
  !! @version 1.0
  !! @date 2026-04-05

  USE read_input
  USE OMP_LIB

    ! Workflow summary:
    ! - cluster encounters from sorted array inputs
    ! - track parent/count arrays during merges
    ! - support mean linkage
    ! - write cluster members and summary files to disk

  contains

    !> Perform hierarchical clustering on a distance array.
    !!
    !! This routine performs hierarchical clustering using mean linkage on
    !! the input `array` (shape `(n)`) of
    !! sorted per-encounter values. The routine mutates `array`
    !! during merging and terminates merging when clusters are farther apart
    !! than an internally computed threshold. Results are written to files
    !! named using `output_name` and complexes information is optionally
    !! recorded via the `complexes` structure.
    !!
    !! @param[in,out] array         Array of per-encounter values (n)
    !!                              This array is updated during clustering
    !! @param[in]     n             Number of elements / encounters
    !! @param[in]     output_name   Base name for output files
    !! @param[in]     complexes     Structure containing complex text lines/headers
    subroutine linkage_clustering(array, n, output_name, complexes)

      IMPLICIT NONE

      real (kind=8), dimension(:), allocatable, intent(inout) :: array
      integer, intent(in) :: n
      character*128, intent(in) :: output_name
      type(type_assoc_file) :: complexes

      integer :: i, j, k, min_i, min_j, merge_counter, remaining_clusters, num_dist
      integer, dimension(n) :: cluster_size ! To track the size of each cluster
      real (kind=8):: min_dist, dist, mean_dist, standard_deviation, dist_threshold
      real (kind=8):: sum_dist, sum_sq_dist, old_val, middle_range
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
    !   dispersion_range_coeff = (array(n) - array(1))/(array(n) + array(1))  ! Relative range of values in the array
      
      factor = 1 - abs(middle_range - (mean_dist + 2*standard_deviation)) / (middle_range + mean_dist + 2*standard_deviation)  ! Scale factor based on how close mean is to middle of range

      alpha = standard_deviation / sqrt(standard_deviation**2 + mean_dist**2)
      dist_threshold = (alpha*mean_dist - (1-alpha)*standard_deviation) / factor  ! Scale threshold by relative range to adapt to different value distributions
    !   dist_threshold = (alpha*mean_dist - (1-alpha)*standard_deviation)  ! Scale threshold by relative range to adapt to different value distributions
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
    print *, 'mod_clust_algorithm v1.0 | Starting mean-linkage clustering (1D array).'
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

          ! Mean linkage: weighted average of merged cluster values.
          array(min_i) = cluster_average(min_i)

          ! Update cluster size after array update (keeps weights correct)
          cluster_size(min_i) = cluster_size(min_i) + cluster_size(min_j)

          ! Re-sort: bubble the merged element rightward to maintain sorted order.
          ! The merged mean can only move right within a sorted ascending array.
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

    end subroutine linkage_clustering

    !> Write cluster membership lists to a text file.
    !!
    !! Outputs a file named `<output_name>_clusters.txt` containing the
    !! members of each active cluster in a tabular text format.
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

      ! Write the content of each active cluster into the file in tabular format
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
