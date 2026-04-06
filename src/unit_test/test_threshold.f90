!> Unit tests for threshold operations in the mod_threshold module
!!
!! Tests include:
!!   - Z-coordinate array calculation
!!   - Minimum distance calculation between atoms
!!   - Angle computation for different dimensions
!!   - Array sorting with encounter index tracking
!!
!! Compile with: gfortran -fopenmp -o test_threshold test_threshold.f90 ../src/mod_threshold.f90 ../src/maths.f90 ../src/mod_pdb.f90
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0

program test_threshold
    use mod_threshold
    use maths
    implicit none

    integer :: num_tests, num_passed, num_failed
    
    num_tests = 0
    num_passed = 0
    num_failed = 0

    print *, "========================================="
    print *, "UNIT TESTS FOR THRESHOLD MODULE"
    print *, "========================================="

    ! Test Z-coordinate calculation
    call test_array_z_coord(num_tests, num_passed, num_failed)
    
    ! Test minimum distance calculation
    call test_array_atoms_dist(num_tests, num_passed, num_failed)
    
    ! Test angle calculation
    call test_array_angle(num_tests, num_passed, num_failed)
    
    ! Test array sorting
    call test_sort_array(num_tests, num_passed, num_failed)

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

    !*******************************************************************************
    !> @brief Unit test for array_z_coord routine.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2024-06-09
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
    subroutine test_array_z_coord(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:), allocatable :: array
        real(kind=8), dimension(:, :), allocatable :: solute_crds, trans_vector, rot1, rot2
        real(kind=8), dimension(3) :: xc1, xc2
        integer :: n, nb_atoms, i
        logical :: test_passed

        print *, ""
        print *, "Testing array_z_coord..."

        ! Test 1: Identity transformation - Z values should match input
        num_tests = num_tests + 1
        n = 2
        nb_atoms = 1
        allocate(array(n), solute_crds(nb_atoms, 3), trans_vector(n, 3), rot1(n, 3), rot2(n, 3))
        
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        solute_crds(1, :) = [0.0d0, 0.0d0, 5.0d0]  ! Z = 5
        
        ! Identity transformations
        trans_vector(1, :) = [0.0d0, 0.0d0, 0.0d0]
        trans_vector(2, :) = [0.0d0, 0.0d0, 0.0d0]
        rot1(1, :) = [1.0d0, 0.0d0, 0.0d0]
        rot1(2, :) = [1.0d0, 0.0d0, 0.0d0]
        rot2(1, :) = [0.0d0, 1.0d0, 0.0d0]
        rot2(2, :) = [0.0d0, 1.0d0, 0.0d0]
        
        call array_z_coord(array, n, nb_atoms, xc1, xc2, trans_vector, rot1, rot2, solute_crds)
        
        test_passed = (abs(array(1) - 5.0d0) < 1d-8 .and. abs(array(2) - 5.0d0) < 1d-8)
        
        if (test_passed) then
            print *, "  [PASS] Test 1: Z-coordinate identity transformation"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Z-coordinate identity transformation"
            print *, "    Expected: [5.0, 5.0], Got: [", array(1), ", ", array(2), "]"
            num_failed = num_failed + 1
        end if
        deallocate(array, solute_crds, trans_vector, rot1, rot2)

        ! Test 2: Different Z values with translations
        num_tests = num_tests + 1
        n = 3
        nb_atoms = 1
        allocate(array(n), solute_crds(nb_atoms, 3), trans_vector(n, 3), rot1(n, 3), rot2(n, 3))
        
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        solute_crds(1, :) = [0.0d0, 0.0d0, 1.0d0]
        
        do i = 1, n
            trans_vector(i, :) = [0.0d0, 0.0d0, real(i-1, 8)]
            rot1(i, :) = [1.0d0, 0.0d0, 0.0d0]
            rot2(i, :) = [0.0d0, 1.0d0, 0.0d0]
        end do
        
        call array_z_coord(array, n, nb_atoms, xc1, xc2, trans_vector, rot1, rot2, solute_crds)
        
        test_passed = (abs(array(1) - 1.0d0) < 1d-8 .and. abs(array(2) - 2.0d0) < 1d-8 .and. abs(array(3) - 3.0d0) < 1d-8)
        
        if (test_passed) then
            print *, "  [PASS] Test 2: Z-coordinate with different translations"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Z-coordinate with translations"
            print *, "    Expected: [1.0, 2.0, 3.0], Got: [", array(1), ", ", array(2), ", ", array(3), "]"
            num_failed = num_failed + 1
        end if
        deallocate(array, solute_crds, trans_vector, rot1, rot2)

    end subroutine test_array_z_coord


    !*******************************************************************************
    !> @brief Unit test for array_atoms_dist routine.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2024-06-09
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
    subroutine test_array_atoms_dist(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:), allocatable :: array
        real(kind=8), dimension(:, :), allocatable :: solute1_crds, solute2_crds, trans_vector, rot1, rot2
        real(kind=8), dimension(3) :: xc1, xc2
        integer :: n, nb_atoms
        logical :: test_passed

        print *, ""
        print *, "Testing array_atoms_dist..."

        ! Test 1: Identical solutes - minimum distance should be 0
        num_tests = num_tests + 1
        n = 1
        nb_atoms = 1
        allocate(array(n), solute1_crds(nb_atoms, 3), solute2_crds(nb_atoms, 3))
        allocate(trans_vector(n, 3), rot1(n, 3), rot2(n, 3))
        
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        solute1_crds(1, :) = [1.0d0, 2.0d0, 3.0d0]
        solute2_crds(1, :) = [1.0d0, 2.0d0, 3.0d0]  ! Same point
        
        trans_vector(1, :) = [0.0d0, 0.0d0, 0.0d0]
        rot1(1, :) = [1.0d0, 0.0d0, 0.0d0]
        rot2(1, :) = [0.0d0, 1.0d0, 0.0d0]
        
        call array_atoms_dist(array, n, nb_atoms, xc1, xc2, trans_vector, rot1, rot2, &
                             solute1_crds, solute2_crds)
        
        test_passed = (abs(array(1)) < 1d-8)
        
        if (test_passed) then
            print *, "  [PASS] Test 1: Identical solutes (distance = 0)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Identical solutes"
            print *, "    Expected: 0.0, Got: ", array(1)
            num_failed = num_failed + 1
        end if
        deallocate(array, solute1_crds, solute2_crds, trans_vector, rot1, rot2)

        ! Test 2: Known distance with unit displacement
        num_tests = num_tests + 1
        n = 1
        nb_atoms = 1
        allocate(array(n), solute1_crds(nb_atoms, 3), solute2_crds(nb_atoms, 3))
        allocate(trans_vector(n, 3), rot1(n, 3), rot2(n, 3))
        
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        solute1_crds(1, :) = [0.0d0, 0.0d0, 0.0d0]
        solute2_crds(1, :) = [1.0d0, 0.0d0, 0.0d0]  ! 1 unit away
        
        trans_vector(1, :) = [0.0d0, 0.0d0, 0.0d0]
        rot1(1, :) = [1.0d0, 0.0d0, 0.0d0]
        rot2(1, :) = [0.0d0, 1.0d0, 0.0d0]
        
        call array_atoms_dist(array, n, nb_atoms, xc1, xc2, trans_vector, rot1, rot2, &
                             solute1_crds, solute2_crds)
        
        test_passed = (abs(array(1) - 1.0d0) < 1d-8)
        
        if (test_passed) then
            print *, "  [PASS] Test 2: Unit displacement (distance = 1.0)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Unit displacement"
            print *, "    Expected: 1.0, Got: ", array(1)
            num_failed = num_failed + 1
        end if
        deallocate(array, solute1_crds, solute2_crds, trans_vector, rot1, rot2)

    end subroutine test_array_atoms_dist


    !*******************************************************************************
    !> @brief Unit test for array_angle routine.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2024-06-09
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
    subroutine test_array_angle(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:), allocatable :: array
        real(kind=8), dimension(:, :), allocatable :: trans_vector, rot1, rot2
        real(kind=8), dimension(3) :: xc1, xc2, point1a, point1b, point2a, point2b
        integer :: n, nb_atoms
        logical :: test_passed

        print *, ""
        print *, "Testing array_angle..."

        ! Test 1: Anti-parallel vectors (180 degrees)
        num_tests = num_tests + 1
        n = 1
        nb_atoms = 2
        allocate(array(n), trans_vector(n, 3), rot1(n, 3), rot2(n, 3))
        
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        point1a = [0.0d0, 0.0d0, 0.0d0]
        point1b = [1.0d0, 0.0d0, 0.0d0]
        point2a = [0.0d0, 0.0d0, 0.0d0]
        point2b = [1.0d0, 0.0d0, 0.0d0]
        
        trans_vector(1, :) = [0.0d0, 0.0d0, 0.0d0]
        rot1(1, :) = [1.0d0, 0.0d0, 0.0d0]
        rot2(1, :) = [0.0d0, 1.0d0, 0.0d0]
        
        call array_angle(array, n, nb_atoms, xc1, xc2, trans_vector, rot1, rot2, &
                        point1a, point1b, point2a, point2b, 3)
        
        test_passed = (abs(array(1) - 180.0d0) < 1d-4)
        
        if (test_passed) then
            print *, "  [PASS] Test 1: Anti-parallel vectors (180 degrees)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Anti-parallel vectors"
            print *, "    Expected: 180.0, Got: ", array(1)
            num_failed = num_failed + 1
        end if
        deallocate(array, trans_vector, rot1, rot2)

        ! Test 2: Perpendicular vectors (90 degrees)
        num_tests = num_tests + 1
        n = 1
        nb_atoms = 2
        allocate(array(n), trans_vector(n, 3), rot1(n, 3), rot2(n, 3))
        
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        point1a = [0.0d0, 0.0d0, 0.0d0]
        point1b = [1.0d0, 0.0d0, 0.0d0]
        point2a = [0.0d0, 0.0d0, 0.0d0]
        point2b = [0.0d0, 1.0d0, 0.0d0]  ! Perpendicular
        
        trans_vector(1, :) = [0.0d0, 0.0d0, 0.0d0]
        rot1(1, :) = [1.0d0, 0.0d0, 0.0d0]
        rot2(1, :) = [0.0d0, 1.0d0, 0.0d0]
        
        call array_angle(array, n, nb_atoms, xc1, xc2, trans_vector, rot1, rot2, &
                        point1a, point1b, point2a, point2b, 3)
        
        test_passed = (abs(array(1) - 90.0d0) < 1d-4)
        
        if (test_passed) then
            print *, "  [PASS] Test 2: Perpendicular vectors (90 degrees)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Perpendicular vectors"
            print *, "    Expected: 90.0, Got: ", array(1)
            num_failed = num_failed + 1
        end if
        deallocate(array, trans_vector, rot1, rot2)

    end subroutine test_array_angle


    !*******************************************************************************
    !> @brief Unit test for sort_array routine.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2024-06-09
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
    subroutine test_sort_array(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:), allocatable :: arr
        integer, dimension(:), allocatable :: indexes
        integer :: n, i
        logical :: test_passed

        print *, ""
        print *, "Testing sort_array..."

        ! Test 1: Simple array sort
        num_tests = num_tests + 1
        n = 5
        allocate(arr(n), indexes(n))
        
        arr = [3.0d0, 1.0d0, 4.0d0, 1.5d0, 2.0d0]
        do i = 1, n
            indexes(i) = i
        end do
        
        call sort_array(arr, indexes)
        
        test_passed = (abs(arr(1) - 1.0d0) < 1d-10 .and. abs(arr(2) - 1.5d0) < 1d-10 .and. &
                      abs(arr(3) - 2.0d0) < 1d-10 .and. abs(arr(4) - 3.0d0) < 1d-10 .and. &
                      abs(arr(5) - 4.0d0) < 1d-10)
        
        if (test_passed) then
            print *, "  [PASS] Test 1: Array sorted correctly"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Array sort"
            print *, "    Got: [", arr(1), ", ", arr(2), ", ", arr(3), ", ", arr(4), ", ", arr(5), "]"
            num_failed = num_failed + 1
        end if
        deallocate(arr, indexes)

        ! Test 2: Indexes are reordered consistently
        num_tests = num_tests + 1
        n = 4
        allocate(arr(n), indexes(n))
        
        arr = [5.0d0, 2.0d0, 8.0d0, 1.0d0]
        do i = 1, n
            indexes(i) = i * 10  ! Original indexes: 10, 20, 30, 40
        end do
        
        call sort_array(arr, indexes)
        
        ! After sort: arr = [1.0, 2.0, 5.0, 8.0], indexes should match them
        test_passed = (indexes(1) == 40 .and. indexes(2) == 20 .and. &
                      indexes(3) == 10 .and. indexes(4) == 30)
        
        if (test_passed) then
            print *, "  [PASS] Test 2: Indexes reordered correctly"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Index reordering"
            print *, "    Expected: [40, 20, 10, 30]"
            print *, "    Got: [", indexes(1), ", ", indexes(2), ", ", indexes(3), ", ", indexes(4), "]"
            num_failed = num_failed + 1
        end if
        deallocate(arr, indexes)

        ! Test 3: Already sorted array
        num_tests = num_tests + 1
        n = 4
        allocate(arr(n), indexes(n))
        
        arr = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
        do i = 1, n
            indexes(i) = i
        end do
        
        call sort_array(arr, indexes)
        
        test_passed = (abs(arr(1) - 1.0d0) < 1d-10 .and. indexes(1) == 1)
        
        if (test_passed) then
            print *, "  [PASS] Test 3: Already sorted array unchanged"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 3: Already sorted array"
            num_failed = num_failed + 1
        end if
        deallocate(arr, indexes)

    end subroutine test_sort_array

end program test_threshold
