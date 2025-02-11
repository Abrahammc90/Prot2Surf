module mod_clust_algorithm

  USE read_input

  contains

    subroutine mean_linkage_clustering(matrix, n, ncls, output_name, complexes, opt_array)

        IMPLICIT NONE

        real (kind=8), dimension(:, :), intent(inout) :: matrix
        real (kind=8), dimension(:), intent(inout), optional :: opt_array
        character*128, intent(in) :: output_name
        integer, intent(in) :: n, ncls
        type(type_assoc_file) :: complexes

        integer :: i, j, k, min_i, min_j, merge_counter, remaining_clusters
        integer, dimension(n) :: cluster_size ! To track the size of each cluster
        real (kind=8):: min_dist, dist, standard_deviation
        logical, dimension(n) :: active_points

        integer, dimension(n, n) :: cluster_indexes, cluster_indexes_sorted ! Store indexes of each cluster
        integer, dimension(n) :: cluster_count  ! Number of elements in each cluster

        integer, dimension(n) :: representative_indexes  ! Store most representative value index for each cluster
        real (kind=8), dimension(n) :: representative_values  ! Store representative values for each cluster
        real (kind=8), dimension(n) :: cluster_average, cluster_sd

    
        active_points = .true.
        cluster_size = 1
        remaining_clusters = n
        merge_counter = 0

        ! Initialize cluster indexes with the initial clusters
        do i = 1, n
          cluster_indexes(i, 1) = i  ! Each initial cluster has itself as the only member
          cluster_count(i) = 1       ! Initially, each cluster has one member
        end do

        ! Main loop for clustering
        do while (remaining_clusters > ncls)
          ! Find the closest pair of clusters
          min_dist = 1.0e30
          min_i = -1
          min_j = -1
        
          do i = 1, n
            if (.not. active_points(i)) cycle
            do j = i + 1, n
              if (.not. active_points(j)) cycle
              dist = matrix(i, j)
              if (dist < min_dist) then
                min_dist = dist
                min_i = i
                min_j = j
              end if
            end do
          end do
      
          ! Merge clusters min_i and min_j
          merge_counter = merge_counter + 1

          !if (remaining_clusters == 176) then
          !  call write_clusters(tot_encounters, cluster_indexes, active_clusters)
          !  STOP
          !end if

          ! Combine the indexes and values of the two clusters
          if (merge_counter <= n) then
            do k = 1, cluster_count(min_j)
              cluster_indexes(min_i, cluster_count(min_i) + k) = cluster_indexes(min_j, k)
            end do
            cluster_average(min_i) = cluster_average(min_i) + cluster_average(min_j)
            cluster_count(min_i) = cluster_count(min_i) + cluster_count(min_j)
          end if

          !if (remaining_clusters == 176) print *, 'patata'
      
          ! Mark cluster min_j as inactive and update the size of min_i
          active_points(min_j) = .false.
          cluster_size(min_i) = cluster_size(min_i) + cluster_size(min_j)
          remaining_clusters = remaining_clusters - 1

          

          ! Update distances for the new merged cluster using average linkage
          
          do i = 1, n
            if (i == min_i .or. i == min_j .or. .not. active_points(i)) cycle
            matrix(min_i, i) = (matrix(min_i, i) * cluster_size(min_i) + &
            matrix(min_j, i) * cluster_size(min_j)) / &
                                         (cluster_size(min_i) + cluster_size(min_j))
                                         matrix(i, min_i) = matrix(min_i, i)
          end do

          if ( mod(remaining_clusters, 100) == 0 ) then
            print *, 'remaining clusters', remaining_clusters
          end if

        end do


        ! Find the most representative value for each active cluster and calculate average and standard deviation
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
              !cluster_average(i) = dist
              !standard_deviation = 0
              !do k = 1, cluster_count(i)
              !  standard_deviation = standard_deviation + (cluster_average(i) - &
              !  matrix(cluster_indexes(i, j), cluster_indexes(i, k)))**2
              !end do
              !standard_deviation = sqrt(standard_deviation/cluster_count(i))
              !cluster_sd(i) = standard_deviation
              if (dist < min_dist) then
                  min_dist = dist
                  representative_indexes(i) = cluster_indexes(i, j)
                  representative_values(i) = min_dist
              end if
          end do
        end do

        !Calculate average and standard deviation
        do i = 1, n
          if (.not. active_points(i)) cycle
          dist = 0.0
          do j = 1, cluster_count(i)
              !do k = 1, cluster_count(i)
                  !dist = dist + matrix(cluster_indexes(i, j), cluster_indexes(i, k))
              !end do
              if ( present(opt_array) ) then
                dist = dist + opt_array(cluster_indexes(i, j))
              else
                dist = dist + matrix(i, cluster_indexes(i, j))
              end if
          end do
          dist = dist / cluster_count(i)
          cluster_average(i) = dist
          standard_deviation = 0.0
          do j = 1, cluster_count(i)
            if ( present(opt_array) ) then
              standard_deviation = standard_deviation + (cluster_average(i) - &
              opt_array(cluster_indexes(i, j)))**2
            else
              standard_deviation = standard_deviation + (cluster_average(i) - &
              matrix(i, cluster_indexes(i, j)))**2
            end if
          end do
          standard_deviation = sqrt(standard_deviation/cluster_count(i))
          cluster_sd(i) = standard_deviation
        end do

        !write(*,*) cluster_indexes(45, :)

        if ( present(opt_array) ) then
          call sort_complexes(n, cluster_indexes, cluster_count, active_points, &
          representative_indexes, cluster_average, cluster_sd, opt_array)
        !else
        !  call sort_complexes(n, cluster_indexes, cluster_count, cluster_indexes_sorted, matrix)
        end if

        call write_cluster_elements(n, cluster_indexes, cluster_count, active_points, output_name)
        call write_clust_info(n, representative_indexes, active_points, cluster_count, cluster_average, cluster_sd, output_name)
        call write_complexes(n, cluster_indexes, cluster_count, active_points, output_name, complexes)


        print *, 'Clustering complete.'
    
    end subroutine mean_linkage_clustering

    subroutine sort_complexes(tot_encounters, cluster_indexes, cluster_count, active_clusters, &
      representative_indexes, cluster_average, cluster_sd, array)

      !real (kind=8), dimension(:, :), intent(in) :: matrix
      real (kind=8), dimension(:), intent(inout) :: array
      integer, intent(in) :: tot_encounters
      integer, dimension(:, :), intent(inout) :: cluster_indexes ! Store indexes of each cluster
      integer, dimension(:), intent(inout) :: cluster_count  ! Number of elements in each cluster
      integer, dimension(:), intent(inout) :: representative_indexes
      real (kind=8), dimension(:), intent(inout) :: cluster_average, cluster_sd
      logical, dimension(:), intent(inout) :: active_clusters
      
      integer :: i, j, min_i, min_j, sorted_i, encounters_found, temp_j
      real (kind=8) :: temp
      integer, dimension(tot_encounters) :: sorted_indexes
      !logical, dimension(tot_encounters, tot_encounters) :: indexes_checked_matrix
      integer, dimension(tot_encounters, tot_encounters) :: temp_cluster_indexes
      integer, dimension(tot_encounters) :: temp_representative_indexes, temp_cluster_count
      real (kind=8), dimension(tot_encounters) :: temp_array
      real (kind=8), dimension(tot_encounters) :: temp_cluster_average, temp_cluster_sd
      logical, dimension(tot_encounters) :: temp_active_clusters
      !integer, dimension(:), allocatable :: temp_cluster
      !real (kind=8), dimension(:), allocatable :: temp_cluster_value

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
        temp_cluster_count(i) = cluster_count(sorted_i)
        temp_active_clusters(i) = active_clusters(sorted_i)
        temp_cluster_average(i) = cluster_average(sorted_i)
        temp_cluster_sd(i) = cluster_sd(sorted_i)
        temp_representative_indexes(i) = representative_indexes(sorted_i)

        encounters_found = 0
        do j = 1, tot_encounters

          if ( ANY( cluster_indexes(sorted_i,:) == sorted_indexes(j) ) ) then
            encounters_found = encounters_found + 1
            temp_cluster_indexes(i, encounters_found) = sorted_indexes(j)

            if ( encounters_found == cluster_count(sorted_i) ) EXIT
          end if
        
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

    subroutine write_clust_info(tot_encounters, representative_indexes, active_clusters, &
      clust_count, cluster_average, cluster_sd, output_name)

      IMPLICIT NONE
      integer, intent(in) :: tot_encounters
      integer, dimension(tot_encounters), intent(in) :: clust_count
      logical, dimension(:), intent(in) :: active_clusters
      integer, dimension(:), intent(in) :: representative_indexes
      character*128, intent(in) :: output_name
      integer :: i, j, unit_number, clust_number
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

    subroutine write_complexes(tot_encounters, cluster_indexes, cluster_count, active_clusters, &
                                output_name, complexes)

      

      IMPLICIT NONE
      integer, intent(in) :: tot_encounters
      integer, dimension(:, :), intent(in) :: cluster_indexes
      logical, dimension(:), intent(in) :: active_clusters
      integer, dimension(:), intent(in) :: cluster_count
      character*128, intent(in) :: output_name
      integer :: i, j, unit_number, clust_number, nfile, encounter_index
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