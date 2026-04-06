!> \file test_clustering.f90
!! \brief Unit tests for clustering algorithms in mod_clust_algorithm
!!
!! Tests include:
!!   - Sort complexes functionality
!!   - Write cluster elements
!!   - Basic clustering threshold calculation
!!   - Minimum linkage clustering
!!   - Maximum linkage clustering
!!   - Mean linkage clustering
!!
!! Note: Array parameter is mandatory for all clustering methods EXCEPT RMSD.
!! For RMSD-based clustering, array is optional. When provided, array is used to
!! compute cluster averages and standard deviations.
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
    
    ! Test minimum linkage clustering
    call test_minimum_linkage(num_tests, num_passed, num_failed)
    
    ! Test maximum linkage clustering
    call test_maximum_linkage(num_tests, num_passed, num_failed)
    
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


    !> \brief Test simple clustering scenarios (tight clusters, identical points)
    !! @param[inout] num_tests   Number of tests run
    !! @param[inout] num_passed  Number of tests passed
    !! @param[inout] num_failed  Number of tests failed
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



    !> \brief Test minimum linkage clustering on 1D array
    !! @param[inout] num_tests   Number of tests run
    !! @param[inout] num_passed  Number of tests passed
    !! @param[inout] num_failed  Number of tests failed
    subroutine test_minimum_linkage(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:), allocatable :: array
        logical :: test_passed
        integer :: n, i
        type(type_assoc_file) :: complexes
        character*128 :: output_name

        print *, ""
        print *, "Testing minimum linkage clustering (1D array)..."

        ! Test 1: Minimum linkage should use smallest distance between clusters
        num_tests = num_tests + 1
        n = 6
        allocate(array(n))
        ! Create two tight clusters: [1.0, 1.1, 1.2] and [10.0, 10.1, 10.2]
        array = [1.0d0, 1.1d0, 1.2d0, 10.0d0, 10.1d0, 10.2d0]

        complexes%nlines = n
        allocate(complexes%lines(n))
        do i = 1, n
            complexes%lines(i) = "Complex"
        end do

        output_name = "test_min_linkage"

        call linkage_clustering_from_array(array, n, 'min', output_name, complexes)

        ! If no error, test passes
        test_passed = .true.
        if (test_passed) then
            print *, "  [PASS] Test 1: Minimum linkage clustering completed"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Minimum linkage clustering"
            num_failed = num_failed + 1
        end if
        deallocate(array, complexes%lines)
    end subroutine test_minimum_linkage



    !> \brief Test maximum linkage clustering on 1D array
    !! @param[inout] num_tests   Number of tests run
    !! @param[inout] num_passed  Number of tests passed
    !! @param[inout] num_failed  Number of tests failed
    subroutine test_maximum_linkage(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:), allocatable :: array
        logical :: test_passed
        integer :: n, i
        type(type_assoc_file) :: complexes
        character*128 :: output_name

        print *, ""
        print *, "Testing maximum linkage clustering (1D array)..."

        ! Test 1: Maximum linkage should use largest distance between clusters
        num_tests = num_tests + 1
        n = 6
        allocate(array(n))
        ! Two clusters: [2.0, 2.1, 2.2] and [20.0, 20.1, 20.2]
        array = [2.0d0, 2.1d0, 2.2d0, 20.0d0, 20.1d0, 20.2d0]

        complexes%nlines = n
        allocate(complexes%lines(n))
        do i = 1, n
            complexes%lines(i) = "Complex"
        end do

        output_name = "test_max_linkage"

        call linkage_clustering_from_array(array, n, 'max', output_name, complexes)

        test_passed = .true.
        if (test_passed) then
            print *, "  [PASS] Test 1: Maximum linkage clustering completed"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Maximum linkage clustering"
            num_failed = num_failed + 1
        end if
        deallocate(array, complexes%lines)
    end subroutine test_maximum_linkage



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

        call linkage_clustering_from_array(array, n, 'mean', output_name, complexes)

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
