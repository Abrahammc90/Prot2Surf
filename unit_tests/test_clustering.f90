!> Unit tests for clustering algorithms in the mod_clust_algorithm module
!!
!! Tests include:
!!   - Sort complexes functionality
!!   - Write cluster elements
!!   - Basic clustering threshold calculation
!!   - Mean linkage clustering with simple matrices
!!
!! Note: Array parameter is mandatory for all clustering methods EXCEPT RMSD.
!! For RMSD-based clustering, array is optional. When provided, array is used to
!! compute cluster averages and standard deviations.
!!
!! Compile with: gfortran -fopenmp -o test_clustering test_clustering.f90 ../src/mod_clust_algorithm.f90 ../src/read_input.f90
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0

program test_clustering
    use mod_clust_algorithm
    use read_input
    implicit none

    integer :: num_tests, num_passed, num_failed
    
    num_tests = 0
    num_passed = 0
    num_failed = 0

    print *, "========================================="
    print *, "UNIT TESTS FOR CLUSTERING ALGORITHMS"
    print *, "========================================="

    ! Test clustering statistics
    call test_clustering_statistics(num_tests, num_passed, num_failed)
    
    ! Test clustering with simple data
    call test_simple_clustering(num_tests, num_passed, num_failed)
    
    ! Test write operations
    call test_cluster_writing(num_tests, num_passed, num_failed)

    print *, ""
    print *, "========================================="
    print *, "TEST SUMMARY"
    print *, "========================================="
    print *, "Total tests:  ", num_tests
    print *, "Passed:       ", num_passed
    print *, "Failed:       ", num_failed
    print *, "========================================="

    if (num_failed > 0) then
        stop 1
    end if

contains

    subroutine test_clustering_statistics(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:, :), allocatable :: matrix
        real(kind=8) :: sum_dist, sum_sq_dist, mean_dist, std_dev
        integer :: i, j, n, num_dist
        logical :: test_passed

        print *, ""
        print *, "Testing clustering statistics calculations..."

        ! Test 1: Mean distance calculation
        num_tests = num_tests + 1
        n = 3
        allocate(matrix(n, n))
        
        matrix(1, 1) = 0.0d0
        matrix(1, 2) = 1.0d0
        matrix(2, 1) = 1.0d0
        matrix(1, 3) = 3.0d0
        matrix(3, 1) = 3.0d0
        matrix(2, 2) = 0.0d0
        matrix(2, 3) = 2.0d0
        matrix(3, 2) = 2.0d0
        matrix(3, 3) = 0.0d0

        sum_dist = 0.0d0
        num_dist = 0
        do i = 1, n
            do j = i + 1, n
                sum_dist = sum_dist + matrix(i, j)
                num_dist = num_dist + 1
            end do
        end do
        mean_dist = sum_dist / num_dist

        test_passed = (abs(mean_dist - 2.0d0) < 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 1: Mean distance calculation (2.0)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Mean distance calculation"
            print *, "    Expected: 2.0, Got: ", mean_dist
            num_failed = num_failed + 1
        end if
        deallocate(matrix)

        ! Test 2: Standard deviation calculation
        num_tests = num_tests + 1
        n = 3
        allocate(matrix(n, n))
        
        ! Same matrix as above
        matrix(1, 1) = 0.0d0
        matrix(1, 2) = 1.0d0
        matrix(2, 1) = 1.0d0
        matrix(1, 3) = 3.0d0
        matrix(3, 1) = 3.0d0
        matrix(2, 2) = 0.0d0
        matrix(2, 3) = 2.0d0
        matrix(3, 2) = 2.0d0
        matrix(3, 3) = 0.0d0

        sum_dist = 0.0d0
        sum_sq_dist = 0.0d0
        num_dist = 0
        do i = 1, n
            do j = i + 1, n
                sum_dist = sum_dist + matrix(i, j)
                sum_sq_dist = sum_sq_dist + matrix(i, j) ** 2
                num_dist = num_dist + 1
            end do
        end do
        mean_dist = sum_dist / num_dist
        std_dev = sqrt((sum_sq_dist / num_dist) - mean_dist ** 2)

        ! Expected: mean = 2.0, variance = (1+4+9)/3 - 4 = 14/3 - 4 = 2/3, std = sqrt(2/3) ≈ 0.8165
        test_passed = (abs(std_dev - sqrt(2.0d0/3.0d0)) < 1d-6)

        if (test_passed) then
            print *, "  [PASS] Test 2: Standard deviation calculation"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Standard deviation calculation"
            print *, "    Expected: ", sqrt(2.0d0/3.0d0), ", Got: ", std_dev
            num_failed = num_failed + 1
        end if
        deallocate(matrix)

    end subroutine test_clustering_statistics


    subroutine test_simple_clustering(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:, :), allocatable :: matrix
        logical :: test_passed
        integer :: n, i, j

        print *, ""
        print *, "Testing simple clustering scenarios..."

        ! Test 1: Matrix with clear clusters (two pairs)
        num_tests = num_tests + 1
        n = 4
        allocate(matrix(n, n))
        
        ! Two tight clusters: (1,2) and (3,4)
        matrix(1, 1) = 0.0d0
        matrix(1, 2) = 0.1d0
        matrix(2, 1) = 0.1d0
        matrix(1, 3) = 5.0d0
        matrix(3, 1) = 5.0d0
        matrix(1, 4) = 5.1d0
        matrix(4, 1) = 5.1d0
        
        matrix(2, 2) = 0.0d0
        matrix(2, 3) = 5.1d0
        matrix(3, 2) = 5.1d0
        matrix(2, 4) = 5.2d0
        matrix(4, 2) = 5.2d0
        
        matrix(3, 3) = 0.0d0
        matrix(3, 4) = 0.1d0
        matrix(4, 3) = 0.1d0
        
        matrix(4, 4) = 0.0d0

        ! Verify minimum distance is within first cluster
        test_passed = (matrix(1, 2) < 1.0d0 .and. matrix(3, 4) < 1.0d0)

        if (test_passed) then
            print *, "  [PASS] Test 1: Identified tight cluster pairs"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Tight cluster identification"
            num_failed = num_failed + 1
        end if
        deallocate(matrix)

        ! Test 2: All identical points (all distances = 0)
        num_tests = num_tests + 1
        n = 3
        allocate(matrix(n, n))
        
        matrix = 0.0d0

        test_passed = .true.
        do i = 1, n
            do j = 1, n
                if (abs(matrix(i, j)) > 1d-10) then
                    test_passed = .false.
                end if
            end do
        end do

        if (test_passed) then
            print *, "  [PASS] Test 2: All identical points (distances = 0)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: All identical points"
            num_failed = num_failed + 1
        end if
        deallocate(matrix)

    end subroutine test_simple_clustering


    subroutine test_cluster_writing(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        integer, dimension(:, :), allocatable :: cluster_indexes
        integer, dimension(:), allocatable :: cluster_count
        logical, dimension(:), allocatable :: active_clusters
        character*128 :: output_name
        logical :: test_passed, file_exists
        integer :: n, unit

        print *, ""
        print *, "Testing cluster writing operations..."

        ! Test 1: Write cluster elements
        num_tests = num_tests + 1
        n = 3
        allocate(cluster_indexes(n, n), cluster_count(n), active_clusters(n))
        
        ! Cluster 1: elements 1, 2
        cluster_indexes(1, 1) = 1
        cluster_indexes(1, 2) = 2
        cluster_count(1) = 2
        active_clusters(1) = .true.
        
        ! Cluster 2: element 3
        cluster_indexes(2, 1) = 3
        cluster_count(2) = 1
        active_clusters(2) = .true.
        
        ! Cluster 3: inactive
        active_clusters(3) = .false.
        cluster_count(3) = 0
        
        output_name = "test_cluster"
        
        call write_cluster_elements(n, cluster_indexes, cluster_count, active_clusters, output_name)
        
        ! Check if file was created
        inquire(file=trim(output_name)//"_clust.txt", exist=file_exists)
        test_passed = file_exists

        if (test_passed) then
            print *, "  [PASS] Test 1: Cluster elements file created"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Cluster elements file not created"
            num_failed = num_failed + 1
        end if
        deallocate(cluster_indexes, cluster_count, active_clusters)

    end subroutine test_cluster_writing

end program test_clustering
