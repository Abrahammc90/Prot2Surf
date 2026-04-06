!> Unit tests for array operations in the mod_array module
!!
!! Tests include:
!!   - Write and read array operations
!!   - Z-coordinate array sorting with encounter metadata reordering
!!   - Minimum-distance array sorting with encounter metadata reordering
!!   - Angle array sorting with encounter metadata reordering
!!   - Low-level merge of sorted segments with metadata reordering
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0

program test_array
    use mod_array
    use mod_assoc
    use maths
    implicit none

    integer :: num_tests, num_passed, num_failed

    num_tests = 0
    num_passed = 0
    num_failed = 0

    print *, "========================================="
    print *, "UNIT TESTS FOR ARRAY MODULE"
    print *, "========================================="

    call test_write_read_array(num_tests, num_passed, num_failed)
    call test_array_z_coord_sorting(num_tests, num_passed, num_failed)
    call test_array_atoms_dist_sorting(num_tests, num_passed, num_failed)
    call test_array_angle_sorting(num_tests, num_passed, num_failed)
    call test_merge_sorted_segments(num_tests, num_passed, num_failed)

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
    !> @brief Unit test for write_array and read_array routines.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2026-04-06
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
    subroutine test_write_read_array(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:), allocatable :: array_orig, array_read
        character*128 :: filename
        integer :: n, i
        logical :: test_passed

        print *, ""
        print *, "Testing write_array and read_array..."

        num_tests = num_tests + 1
        n = 5
        allocate(array_orig(n))
        filename = "test_array_module.dat"
        array_orig = [1.5d0, 2.3d0, 4.7d0, 0.1d0, 9.9d0]

        call write_array(array_orig, filename)
        call read_array(array_read, n, filename)

        test_passed = .true.
        do i = 1, n
            if (abs(array_orig(i) - array_read(i)) > 1d-3) then
                test_passed = .false.
            end if
        end do

        if (test_passed) then
            print *, "  [PASS] Test 1: Write and read array"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Write and read array"
            num_failed = num_failed + 1
        end if

        deallocate(array_orig, array_read)

    end subroutine test_write_read_array


    !*******************************************************************************
    !> @brief Unit test for array_z_coord sorting and complex reordering.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2026-04-06
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
    subroutine test_array_z_coord_sorting(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        type(type_assoc_file) :: complexes
        real(kind=8), dimension(:), allocatable :: array_sorted
        real(kind=8), dimension(:, :), allocatable :: solute_crds
        real(kind=8), dimension(3) :: xc1, xc2
        integer :: n, nb_atoms
        logical :: test_passed

        print *, ""
        print *, "Testing array_z_coord sorting..."

        num_tests = num_tests + 1
        n = 3
        nb_atoms = 1
        allocate(array_sorted(n), solute_crds(nb_atoms, 3))
        call allocate_assoc_object(complexes, n)

        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        solute_crds(1, :) = [0.0d0, 0.0d0, 1.0d0]

        complexes%lines(1) = 'encounter_z3'
        complexes%lines(2) = 'encounter_z1'
        complexes%lines(3) = 'encounter_z2'

        complexes%trans_vector(1, :) = [0.0d0, 0.0d0, 2.0d0]
        complexes%trans_vector(2, :) = [0.0d0, 0.0d0, 0.0d0]
        complexes%trans_vector(3, :) = [0.0d0, 0.0d0, 1.0d0]
        complexes%rot1(1, :) = [1.0d0, 0.0d0, 0.0d0]
        complexes%rot1(2, :) = [1.0d0, 0.0d0, 0.0d0]
        complexes%rot1(3, :) = [1.0d0, 0.0d0, 0.0d0]
        complexes%rot2(1, :) = [0.0d0, 1.0d0, 0.0d0]
        complexes%rot2(2, :) = [0.0d0, 1.0d0, 0.0d0]
        complexes%rot2(3, :) = [0.0d0, 1.0d0, 0.0d0]

        call array_z_coord(complexes, array_sorted, n, nb_atoms, xc1, xc2, &
            complexes%trans_vector, complexes%rot1, complexes%rot2, solute_crds)

        test_passed = abs(array_sorted(1) - 1.0d0) < 1d-8 .and. &
                      abs(array_sorted(2) - 2.0d0) < 1d-8 .and. &
                      abs(array_sorted(3) - 3.0d0) < 1d-8 .and. &
                      trim(complexes%lines(1)) == 'encounter_z1' .and. &
                      trim(complexes%lines(2)) == 'encounter_z2' .and. &
                      trim(complexes%lines(3)) == 'encounter_z3'

        if (test_passed) then
            print *, "  [PASS] Test 1: Z-coordinate array sorted with metadata"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Z-coordinate array sorted with metadata"
            print *, "    array_sorted: ", array_sorted
            print *, "    lines: ", trim(complexes%lines(1)), trim(complexes%lines(2)), trim(complexes%lines(3))
            num_failed = num_failed + 1
        end if

        deallocate(array_sorted, solute_crds)
        deallocate(complexes%lines, complexes%trans_vector, complexes%rot1, complexes%rot2)

    end subroutine test_array_z_coord_sorting


    !*******************************************************************************
    !> @brief Unit test for array_atoms_dist sorting and complex reordering.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2026-04-06
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
    subroutine test_array_atoms_dist_sorting(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        type(type_assoc_file) :: complexes
        real(kind=8), dimension(:), allocatable :: array_sorted
        real(kind=8), dimension(:, :), allocatable :: solute1_crds, solute2_crds
        real(kind=8), dimension(3) :: xc1, xc2
        integer :: n, nb_atoms
        logical :: test_passed

        print *, ""
        print *, "Testing array_atoms_dist sorting..."

        num_tests = num_tests + 1
        n = 3
        nb_atoms = 1
        allocate(array_sorted(n), solute1_crds(1, 3), solute2_crds(nb_atoms, 3))
        call allocate_assoc_object(complexes, n)

        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        solute1_crds(1, :) = [0.0d0, 0.0d0, 0.0d0]
        solute2_crds(1, :) = [1.0d0, 0.0d0, 0.0d0]

        complexes%lines(1) = 'encounter_d3'
        complexes%lines(2) = 'encounter_d1'
        complexes%lines(3) = 'encounter_d2'

        complexes%trans_vector(1, :) = [2.0d0, 0.0d0, 0.0d0]
        complexes%trans_vector(2, :) = [0.0d0, 0.0d0, 0.0d0]
        complexes%trans_vector(3, :) = [1.0d0, 0.0d0, 0.0d0]
        complexes%rot1(1, :) = [1.0d0, 0.0d0, 0.0d0]
        complexes%rot1(2, :) = [1.0d0, 0.0d0, 0.0d0]
        complexes%rot1(3, :) = [1.0d0, 0.0d0, 0.0d0]
        complexes%rot2(1, :) = [0.0d0, 1.0d0, 0.0d0]
        complexes%rot2(2, :) = [0.0d0, 1.0d0, 0.0d0]
        complexes%rot2(3, :) = [0.0d0, 1.0d0, 0.0d0]

        call array_atoms_dist(complexes, array_sorted, n, nb_atoms, xc1, xc2, &
            complexes%trans_vector, complexes%rot1, complexes%rot2, solute1_crds, solute2_crds)

        test_passed = abs(array_sorted(1) - 1.0d0) < 1d-8 .and. &
                      abs(array_sorted(2) - 2.0d0) < 1d-8 .and. &
                      abs(array_sorted(3) - 3.0d0) < 1d-8 .and. &
                      trim(complexes%lines(1)) == 'encounter_d1' .and. &
                      trim(complexes%lines(2)) == 'encounter_d2' .and. &
                      trim(complexes%lines(3)) == 'encounter_d3'

        if (test_passed) then
            print *, "  [PASS] Test 1: Distance array sorted with metadata"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Distance array sorted with metadata"
            print *, "    array_sorted: ", array_sorted
            print *, "    lines: ", trim(complexes%lines(1)), trim(complexes%lines(2)), trim(complexes%lines(3))
            num_failed = num_failed + 1
        end if

        deallocate(array_sorted, solute1_crds, solute2_crds)
        deallocate(complexes%lines, complexes%trans_vector, complexes%rot1, complexes%rot2)

    end subroutine test_array_atoms_dist_sorting


    !*******************************************************************************
    !> @brief Unit test for array_angle sorting and complex reordering.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2026-04-06
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
    subroutine test_array_angle_sorting(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        type(type_assoc_file) :: complexes
        real(kind=8), dimension(:), allocatable :: array_sorted
        real(kind=8), dimension(3) :: xc1, xc2, point1a, point1b, point2a, point2b
        integer :: n, nb_atoms
        logical :: test_passed

        print *, ""
        print *, "Testing array_angle sorting..."

        num_tests = num_tests + 1
        n = 2
        nb_atoms = 2
        allocate(array_sorted(n))
        call allocate_assoc_object(complexes, n)

        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        point1a = [0.0d0, 0.0d0, 0.0d0]
        point1b = [1.0d0, 0.0d0, 0.0d0]
        point2a = [0.0d0, 0.0d0, 0.0d0]
        point2b = [1.0d0, 0.0d0, 0.0d0]

        complexes%lines(1) = 'encounter_a180'
        complexes%lines(2) = 'encounter_a90'

        complexes%trans_vector(1, :) = [0.0d0, 0.0d0, 0.0d0]
        complexes%trans_vector(2, :) = [0.0d0, 0.0d0, 0.0d0]
        complexes%rot1(1, :) = [1.0d0, 0.0d0, 0.0d0]
        complexes%rot2(1, :) = [0.0d0, 1.0d0, 0.0d0]
        complexes%rot1(2, :) = [0.0d0, 1.0d0, 0.0d0]
        complexes%rot2(2, :) = [-1.0d0, 0.0d0, 0.0d0]

        call array_angle(complexes, array_sorted, n, nb_atoms, xc1, xc2, &
            complexes%trans_vector, complexes%rot1, complexes%rot2, &
            point1a, point1b, point2a, point2b, 2)

        test_passed = abs(array_sorted(1) - 90.0d0) < 1d-4 .and. &
                      abs(array_sorted(2) - 180.0d0) < 1d-4 .and. &
                      trim(complexes%lines(1)) == 'encounter_a90' .and. &
                      trim(complexes%lines(2)) == 'encounter_a180'

        if (test_passed) then
            print *, "  [PASS] Test 1: Angle array sorted with metadata"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Angle array sorted with metadata"
            print *, "    array_sorted: ", array_sorted
            print *, "    lines: ", trim(complexes%lines(1)), trim(complexes%lines(2))
            num_failed = num_failed + 1
        end if

        deallocate(array_sorted)
        deallocate(complexes%lines, complexes%trans_vector, complexes%rot1, complexes%rot2)

    end subroutine test_array_angle_sorting


    !*******************************************************************************
    !> @brief Unit test for merge_sorted_segments routine.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2026-04-06
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
    subroutine test_merge_sorted_segments(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        type(type_assoc_file) :: complexes
        real(kind=8), dimension(:), allocatable :: arr, tmp_arr
        character(len=217), dimension(:), allocatable :: tmp_comp
        integer :: n
        logical :: test_passed

        print *, ""
        print *, "Testing merge_sorted_segments..."

        num_tests = num_tests + 1
        n = 6
        allocate(arr(n), tmp_arr(n), tmp_comp(n))
        call allocate_assoc_object(complexes, n)

        arr = [1.0d0, 3.0d0, 5.0d0, 2.0d0, 4.0d0, 6.0d0]

        complexes%lines(1) = 'enc1'
        complexes%lines(2) = 'enc3'
        complexes%lines(3) = 'enc5'
        complexes%lines(4) = 'enc2'
        complexes%lines(5) = 'enc4'
        complexes%lines(6) = 'enc6'

        complexes%trans_vector(1, :) = [1.0d0, 0.0d0, 0.0d0]
        complexes%trans_vector(2, :) = [3.0d0, 0.0d0, 0.0d0]
        complexes%trans_vector(3, :) = [5.0d0, 0.0d0, 0.0d0]
        complexes%trans_vector(4, :) = [2.0d0, 0.0d0, 0.0d0]
        complexes%trans_vector(5, :) = [4.0d0, 0.0d0, 0.0d0]
        complexes%trans_vector(6, :) = [6.0d0, 0.0d0, 0.0d0]
        complexes%rot1 = 0.0d0
        complexes%rot2 = 0.0d0
        complexes%rot1(:, 1) = [1.0d0, 3.0d0, 5.0d0, 2.0d0, 4.0d0, 6.0d0]
        complexes%rot2(:, 1) = [10.0d0, 30.0d0, 50.0d0, 20.0d0, 40.0d0, 60.0d0]

        call merge_sorted_segments(arr, complexes, tmp_arr, tmp_comp, 1, 3, 6)

        test_passed = abs(arr(1) - 1.0d0) < 1d-8 .and. &
                      abs(arr(2) - 2.0d0) < 1d-8 .and. &
                      abs(arr(3) - 3.0d0) < 1d-8 .and. &
                      abs(arr(4) - 4.0d0) < 1d-8 .and. &
                      abs(arr(5) - 5.0d0) < 1d-8 .and. &
                      abs(arr(6) - 6.0d0) < 1d-8 .and. &
                      trim(complexes%lines(1)) == 'enc1' .and. &
                      trim(complexes%lines(2)) == 'enc2' .and. &
                      trim(complexes%lines(3)) == 'enc3' .and. &
                      trim(complexes%lines(4)) == 'enc4' .and. &
                      trim(complexes%lines(5)) == 'enc5' .and. &
                      trim(complexes%lines(6)) == 'enc6' .and. &
                      abs(complexes%trans_vector(4, 1) - 4.0d0) < 1d-8 .and. &
                      abs(complexes%rot1(5, 1) - 5.0d0) < 1d-8 .and. &
                      abs(complexes%rot2(6, 1) - 60.0d0) < 1d-8

        if (test_passed) then
            print *, "  [PASS] Test 1: Merge keeps array and metadata aligned"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Merge keeps array and metadata aligned"
            print *, "    arr: ", arr
            print *, "    lines: ", trim(complexes%lines(1)), trim(complexes%lines(2)), &
                trim(complexes%lines(3)), trim(complexes%lines(4)), trim(complexes%lines(5)), trim(complexes%lines(6))
            num_failed = num_failed + 1
        end if

        deallocate(arr, tmp_arr, tmp_comp)
        deallocate(complexes%lines, complexes%trans_vector, complexes%rot1, complexes%rot2)

    end subroutine test_merge_sorted_segments

end program test_array
