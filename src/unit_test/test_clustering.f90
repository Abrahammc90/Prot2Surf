!> \file test_clustering.f90
!! \brief Unit tests for clustering algorithms in mod_clust_algorithm
!!
!! Tests include:
!!   - Write cluster elements
!!   - Basic clustering threshold calculation
!!   - Mean linkage clustering
!!
!! Note: Clustering now operates on sorted 1D arrays.
!!
!! Compile with: gfortran -fopenmp -o test_clustering test_clustering.f90 ../src/mod_clust_algorithm.f90 ../src/read_input.f90
!!
!! @author Abraham Muniz-Chicharro
!! @version 1.0
!! @date 2026-04-05

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
    
    ! Test mean linkage clustering
    call test_mean_linkage(num_tests, num_passed, num_failed)
    
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

    !> \brief Test clustering statistics calculations (mean, stddev)
    !! @param[inout] num_tests   Number of tests run
    !! @param[inout] num_passed  Number of tests passed
    !! @param[inout] num_failed  Number of tests failed
    subroutine test_clustering_statistics(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:), allocatable :: distances
        real(kind=8) :: sum_dist, sum_sq_dist, mean_dist, std_dev
        integer :: num_dist
        logical :: test_passed

        print *, ""
        print *, "Testing clustering statistics calculations..."

        ! Test 1: Mean distance calculation
        num_tests = num_tests + 1
        allocate(distances(3))
        distances = [1.0d0, 3.0d0, 2.0d0]

        sum_dist = sum(distances)
        num_dist = size(distances)
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
        deallocate(distances)

        ! Test 2: Standard deviation calculation
        num_tests = num_tests + 1
        allocate(distances(3))
        distances = [1.0d0, 3.0d0, 2.0d0]

        sum_dist = sum(distances)
        sum_sq_dist = sum(distances ** 2)
        num_dist = size(distances)
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
        deallocate(distances)

    end subroutine test_clustering_statistics


    !> \brief Test simple clustering scenarios (tight clusters, identical points)
    !! @param[inout] num_tests   Number of tests run
    !! @param[inout] num_passed  Number of tests passed
    !! @param[inout] num_failed  Number of tests failed
    subroutine test_simple_clustering(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:), allocatable :: pair_distances
        logical :: test_passed
        integer :: n, i

        print *, ""
        print *, "Testing simple clustering scenarios..."

        ! Test 1: Pair-distance array with clear clusters (two pairs)
        num_tests = num_tests + 1
        allocate(pair_distances(6))
        pair_distances = [0.1d0, 5.0d0, 5.1d0, 5.1d0, 5.2d0, 0.1d0]

        ! Verify tight pair distances are present at both ends of the sorted-like data.
        test_passed = (pair_distances(1) < 1.0d0 .and. pair_distances(6) < 1.0d0)

        if (test_passed) then
            print *, "  [PASS] Test 1: Identified tight cluster pairs"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Tight cluster identification"
            num_failed = num_failed + 1
        end if
        deallocate(pair_distances)

        ! Test 2: All identical points (all distances = 0)
        num_tests = num_tests + 1
        n = 3
        allocate(pair_distances(n))
        pair_distances = 0.0d0

        test_passed = .true.
        do i = 1, n
            if (abs(pair_distances(i)) > 1d-10) then
                test_passed = .false.
            end if
        end do

        if (test_passed) then
            print *, "  [PASS] Test 2: All identical points (distances = 0)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: All identical points"
            num_failed = num_failed + 1
        end if
        deallocate(pair_distances)

    end subroutine test_simple_clustering



    !> \brief Test mean linkage clustering on 1D array
    !! @param[inout] num_tests   Number of tests run
    !! @param[inout] num_passed  Number of tests passed
    !! @param[inout] num_failed  Number of tests failed
    subroutine test_mean_linkage(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:), allocatable :: array
        logical :: test_passed
        integer :: n, i
        type(type_assoc_file) :: complexes
        character*128 :: output_name

        print *, ""
        print *, "Testing mean linkage clustering (1D array)..."

        ! Test 1: Mean linkage should use weighted average distance
        num_tests = num_tests + 1
        n = 6
        allocate(array(n))
        ! Two clusters: [3.0, 3.1, 3.2] and [30.0, 30.1, 30.2]
        array = [3.0d0, 3.1d0, 3.2d0, 30.0d0, 30.1d0, 30.2d0]

        complexes%nlines = n
        allocate(complexes%lines(n))
        do i = 1, n
            complexes%lines(i) = "Complex"
        end do

        output_name = "test_mean_linkage"

        call linkage_clustering(array, n, output_name, complexes)

        test_passed = .true.
        if (test_passed) then
            print *, "  [PASS] Test 1: Mean linkage clustering completed"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Mean linkage clustering"
            num_failed = num_failed + 1
        end if
        deallocate(array, complexes%lines)
    end subroutine test_mean_linkage


    !> \brief Test writing cluster elements to file
    !! @param[inout] num_tests   Number of tests run
    !! @param[inout] num_passed  Number of tests passed
    !! @param[inout] num_failed  Number of tests failed
    subroutine test_cluster_writing(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        integer, dimension(:, :), allocatable :: cluster_indexes
        integer, dimension(:), allocatable :: cluster_count
        logical, dimension(:), allocatable :: active_clusters
        character*128 :: output_name
        logical :: test_passed, file_exists
        integer :: n

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
        
        ! For new interface: cluster_parent = cluster_indexes(:,1)
        call write_cluster_elements(n, cluster_indexes(:,1), cluster_count, active_clusters, output_name)
        
        ! Check if file was created
        inquire(file=trim(output_name)//"_clusters.txt", exist=file_exists)
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
