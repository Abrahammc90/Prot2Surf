!> Unit tests for mathematical operations in the maths module
!!
!! Comprehensive test suite for all subroutines in the maths module including:
!!
!! - **Cross product**: 
!!    - Basis vectors (i x j = k)
!!    - Reverse order (anticommutative property)
!!    - Self cross product (v x v = 0)
!!    - Arbitrary vectors
!!
!! - **Coordinate transformations (update_complex)**:
!!    - Identity transformations
!!    - Pure translation
!!    - Only xc1 translation
!!    - Only xc2 translation
!!    - Only rotation (90 degrees)
!!    - Combined transformations
!!
!! - **3D angle calculations (vectors_angle_3D)**:
!!    - Parallel vectors (0 degrees)
!!    - Perpendicular vectors (90 degrees)
!!    - Opposite vectors (180 degrees)
!!    - 45 degree angles
!!
!! - **2D angle calculations (vectors_angle_2D)**:
!!    - Parallel vectors in XY plane (0 degrees)
!!    - Perpendicular vectors in XY plane (90 degrees)
!!    - Opposite vectors in XY plane (180 degrees)
!!    - Zero magnitude vectors handling
!!
!! - **RMSD calculation**:
!!    - Identical coordinates (RMSD = 0)
!!    - Simple translations
!!    - Known RMSD values
!!
!! - **Center of geometry (calculate_cog)**:
!!    - Single point
!!    - Points at origin
!!    - Symmetric point distributions
!!    - General cases
!!
!! - **Distance calculation (calculate_distance)**:
!!    - Distance to self (0)
!!    - Unit distances
!!    - 3D diagonal distances
!!    - Arbitrary point pairs
!!
!! Compile with: gfortran -o test_maths test_maths.f90 ../src/maths.f90 ../src/mod_pdb.f90
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0

program test_maths
    use maths
    implicit none

    integer :: num_tests, num_passed, num_failed
    
    num_tests = 0
    num_passed = 0
    num_failed = 0

    print *, "========================================="
    print *, "UNIT TESTS FOR MATHS MODULE"
    print *, "========================================="

    ! Test cross product
    call test_cross_product(num_tests, num_passed, num_failed)
    
    ! Test update_complex
    call test_update_complex(num_tests, num_passed, num_failed)
    
    ! Test angle calculations
    call test_vectors_angle_3D(num_tests, num_passed, num_failed)
    call test_vectors_angle_2D(num_tests, num_passed, num_failed)
    
    ! Test RMSD calculation
    call test_rmsd(num_tests, num_passed, num_failed)
    
    ! Test center of geometry
    call test_calculate_cog(num_tests, num_passed, num_failed)
    
    ! Test distance calculation
    call test_calculate_distance(num_tests, num_passed, num_failed)

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

    subroutine test_cross_product(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(3) :: v1, v2, v3, expected
        logical :: test_passed

        print *, ""
        print *, "Testing cross product..."

        ! Test 1: Standard basis vectors (i x j = k)
        num_tests = num_tests + 1
        v1 = [1.0d0, 0.0d0, 0.0d0]
        v2 = [0.0d0, 1.0d0, 0.0d0]
        expected = [0.0d0, 0.0d0, 1.0d0]
        call cross(v1, v2, v3)
        test_passed = assert_vector_equal(v3, expected, 1d-10)
        if (test_passed) then
            print *, "  [PASS] Test 1: i x j = k"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: i x j = k"
            print *, "    Expected: ", expected
            print *, "    Got:      ", v3
            num_failed = num_failed + 1
        end if

        ! Test 2: Reverse order gives negative result
        num_tests = num_tests + 1
        v1 = [1.0d0, 0.0d0, 0.0d0]
        v2 = [0.0d0, 1.0d0, 0.0d0]
        expected = [0.0d0, 0.0d0, 1.0d0]
        call cross(v2, v1, v3)
        test_passed = assert_vector_equal(v3, -expected, 1d-10)
        if (test_passed) then
            print *, "  [PASS] Test 2: j x i = -k"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: j x i = -k"
            print *, "    Expected: ", -expected
            print *, "    Got:      ", v3
            num_failed = num_failed + 1
        end if

        ! Test 3: Cross product with itself is zero
        num_tests = num_tests + 1
        v1 = [1.0d0, 2.0d0, 3.0d0]
        expected = [0.0d0, 0.0d0, 0.0d0]
        call cross(v1, v1, v3)
        test_passed = assert_vector_equal(v3, expected, 1d-10)
        if (test_passed) then
            print *, "  [PASS] Test 3: v x v = 0"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 3: v x v = 0"
            print *, "    Expected: ", expected
            print *, "    Got:      ", v3
            num_failed = num_failed + 1
        end if

        ! Test 4: Arbitrary vectors
        num_tests = num_tests + 1
        v1 = [1.0d0, 2.0d0, 3.0d0]
        v2 = [4.0d0, 5.0d0, 6.0d0]
        expected = [-3.0d0, 6.0d0, -3.0d0]
        call cross(v1, v2, v3)
        test_passed = assert_vector_equal(v3, expected, 1d-10)
        if (test_passed) then
            print *, "  [PASS] Test 4: Arbitrary vectors"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 4: Arbitrary vectors"
            print *, "    Expected: ", expected
            print *, "    Got:      ", v3
            num_failed = num_failed + 1
        end if

    end subroutine test_cross_product


    subroutine test_update_complex(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(3) :: xc1, xc2, trans, rot_vx, rot_vy
        real(kind=8), dimension(3, 3) :: coords, new_coords, expected
        integer :: i
        logical :: test_passed

        print *, ""
        print *, "Testing coordinate transformation (update_complex)..."

        ! Test 1: Identity transformation (no rotation, no translation)
        num_tests = num_tests + 1
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        trans = [0.0d0, 0.0d0, 0.0d0]
        rot_vx = [1.0d0, 0.0d0, 0.0d0]
        rot_vy = [0.0d0, 1.0d0, 0.0d0]
        coords(1, :) = [1.0d0, 2.0d0, 3.0d0]
        coords(2, :) = [4.0d0, 5.0d0, 6.0d0]
        coords(3, :) = [7.0d0, 8.0d0, 9.0d0]

        call update_complex(xc1, xc2, trans, rot_vx, rot_vy, 3, coords, new_coords)

        ! With identity rotation and no translation, output should be same as input
        test_passed = .true.
        do i = 1, 3
            if (.not. assert_vector_equal(new_coords(i, :), coords(i, :), 1d-8)) then
                test_passed = .false.
            end if
        end do

        if (test_passed) then
            print *, "  [PASS] Test 1: Identity transformation"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Identity transformation"
            num_failed = num_failed + 1
        end if

        ! Test 2: Pure translation
        num_tests = num_tests + 1
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        trans = [1.0d0, 2.0d0, 3.0d0]
        rot_vx = [1.0d0, 0.0d0, 0.0d0]
        rot_vy = [0.0d0, 1.0d0, 0.0d0]
        coords(1, :) = [0.0d0, 0.0d0, 0.0d0]

        call update_complex(xc1, xc2, trans, rot_vx, rot_vy, 1, coords, new_coords)
        expected(1, :) = [1.0d0, 2.0d0, 3.0d0]

        test_passed = assert_vector_equal(new_coords(1, :), expected(1, :), 1d-8)

        if (test_passed) then
            print *, "  [PASS] Test 2: Pure translation"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Pure translation"
            print *, "    Expected: ", expected(1, :)
            print *, "    Got:      ", new_coords(1, :)
            num_failed = num_failed + 1
        end if

        ! Test 3: Only xc1 translation
        num_tests = num_tests + 1
        xc1 = [1.0d0, 1.0d0, 1.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        trans = [0.0d0, 0.0d0, 0.0d0]
        rot_vx = [1.0d0, 0.0d0, 0.0d0]
        rot_vy = [0.0d0, 1.0d0, 0.0d0]
        coords(1, :) = [1.0d0, 1.0d0, 1.0d0]

        call update_complex(xc1, xc2, trans, rot_vx, rot_vy, 1, coords, new_coords)
        ! Expected: input shifted by xc1 center
        expected(1, :) = [0.0d0, 0.0d0, 0.0d0]

        test_passed = assert_vector_equal(new_coords(1, :), expected(1, :), 1d-8)

        if (test_passed) then
            print *, "  [PASS] Test 3: Only xc1 translation"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 3: Only xc1 translation"
            print *, "    Expected: ", expected(1, :)
            print *, "    Got:      ", new_coords(1, :)
            num_failed = num_failed + 1
        end if

        ! Test 4: Only xc2 translation
        num_tests = num_tests + 1
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [1.0d0, 1.0d0, 1.0d0]
        trans = [0.0d0, 0.0d0, 0.0d0]
        rot_vx = [1.0d0, 0.0d0, 0.0d0]
        rot_vy = [0.0d0, 1.0d0, 0.0d0]
        coords(1, :) = [0.0d0, 0.0d0, 0.0d0]

        call update_complex(xc1, xc2, trans, rot_vx, rot_vy, 1, coords, new_coords)
        ! Expected: shifted by xc2
        expected(1, :) = [1.0d0, 1.0d0, 1.0d0]

        test_passed = assert_vector_equal(new_coords(1, :), expected(1, :), 1d-8)

        if (test_passed) then
            print *, "  [PASS] Test 4: Only xc2 translation"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 4: Only xc2 translation"
            print *, "    Expected: ", expected(1, :)
            print *, "    Got:      ", new_coords(1, :)
            num_failed = num_failed + 1
        end if

        ! Test 5: Only rotation (90 degree rotation around z-axis)
        num_tests = num_tests + 1
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        trans = [0.0d0, 0.0d0, 0.0d0]
        ! Rotation matrix for 90 degrees around z-axis: x' = -y, y' = x, z' = z
        rot_vx = [0.0d0, -1.0d0, 0.0d0]  ! New x direction (was negative y)
        rot_vy = [1.0d0, 0.0d0, 0.0d0]   ! New y direction (was positive x)
        coords(1, :) = [1.0d0, 0.0d0, 0.0d0]

        call update_complex(xc1, xc2, trans, rot_vx, rot_vy, 1, coords, new_coords)
        ! Expected: [0, 1, 0] after 90 degree rotation
        expected(1, :) = [0.0d0, 1.0d0, 0.0d0]

        test_passed = assert_vector_equal(new_coords(1, :), expected(1, :), 1d-8)

        if (test_passed) then
            print *, "  [PASS] Test 5: Only rotation (90 degrees)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 5: Only rotation (90 degrees)"
            print *, "    Expected: ", expected(1, :)
            print *, "    Got:      ", new_coords(1, :)
            num_failed = num_failed + 1
        end if

        ! Test 6: Combined xc1 translation + rotation
        num_tests = num_tests + 1
        xc1 = [1.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        trans = [2.0d0, 0.0d0, 0.0d0]
        rot_vx = [0.0d0, -1.0d0, 0.0d0]
        rot_vy = [1.0d0, 0.0d0, 0.0d0]
        coords(1, :) = [1.0d0, 0.0d0, 0.0d0]

        call update_complex(xc1, xc2, trans, rot_vx, rot_vy, 1, coords, new_coords)
        ! Shift by xc1: [0, 0, 0], rotate 90 degrees: [0, 0, 0], translate by trans: [2, 0, 0]
        expected(1, :) = [2.0d0, 0.0d0, 0.0d0]

        test_passed = assert_vector_equal(new_coords(1, :), expected(1, :), 1d-8)

        if (test_passed) then
            print *, "  [PASS] Test 6: Combined xc1 translation + rotation"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 6: Combined xc1 translation + rotation"
            print *, "    Expected: ", expected(1, :)
            print *, "    Got:      ", new_coords(1, :)
            num_failed = num_failed + 1
        end if

        ! Test 7: Combined rotation + trans + xc2 translation
        num_tests = num_tests + 1
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [1.0d0, 0.0d0, 0.0d0]
        trans = [0.0d0, 1.0d0, 0.0d0]
        rot_vx = [0.0d0, -1.0d0, 0.0d0]
        rot_vy = [1.0d0, 0.0d0, 0.0d0]
        coords(1, :) = [1.0d0, 0.0d0, 0.0d0]

        call update_complex(xc1, xc2, trans, rot_vx, rot_vy, 1, coords, new_coords)
        ! Rotate 90 degrees: [0, 1, 0], translate by trans: [0, 2, 0], shift by xc2: [1, 2, 0]
        expected(1, :) = [1.0d0, 2.0d0, 0.0d0]

        test_passed = assert_vector_equal(new_coords(1, :), expected(1, :), 1d-8)

        if (test_passed) then
            print *, "  [PASS] Test 7: Combined rotation + trans + xc2 translation"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 7: Combined rotation + trans + xc2 translation"
            print *, "    Expected: ", expected(1, :)
            print *, "    Got:      ", new_coords(1, :)
            num_failed = num_failed + 1
        end if

    end subroutine test_update_complex


    subroutine test_vectors_angle_3D(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(3) :: v1, v2
        real(kind=8) :: angle
        logical :: test_passed

        print *, ""
        print *, "Testing 3D angle calculation (vectors_angle_3D)..."

        ! Test 1: Parallel vectors (0 degrees)
        num_tests = num_tests + 1
        v1 = [1.0d0, 0.0d0, 0.0d0]
        v2 = [2.0d0, 0.0d0, 0.0d0]
        call vectors_angle_3D(v1, v2, angle)
        test_passed = (abs(angle - 0.0d0) < 1d-6)

        if (test_passed) then
            print *, "  [PASS] Test 1: Parallel vectors (0 degrees)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Parallel vectors"
            print *, "    Expected: 0.0, Got: ", angle
            num_failed = num_failed + 1
        end if

        ! Test 2: Perpendicular vectors (90 degrees)
        num_tests = num_tests + 1
        v1 = [1.0d0, 0.0d0, 0.0d0]
        v2 = [0.0d0, 2.0d0, 0.0d0]
        call vectors_angle_3D(v1, v2, angle)
        test_passed = (abs(angle - 90.0d0) < 1d-6)

        if (test_passed) then
            print *, "  [PASS] Test 2: Perpendicular vectors (90 degrees)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Perpendicular vectors"
            print *, "    Expected: 90.0, Got: ", angle
            num_failed = num_failed + 1
        end if

        ! Test 3: Opposite vectors (180 degrees)
        num_tests = num_tests + 1
        v1 = [5.0d0, 0.0d0, 0.0d0]
        v2 = [-1.0d0, 0.0d0, 0.0d0]
        call vectors_angle_3D(v1, v2, angle)
        test_passed = (abs(angle - 180.0d0) < 1d-6)

        if (test_passed) then
            print *, "  [PASS] Test 3: Opposite vectors (180 degrees)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 3: Opposite vectors"
            print *, "    Expected: 180.0, Got: ", angle
            num_failed = num_failed + 1
        end if

        ! Test 4: 45 degree angle
        num_tests = num_tests + 1
        v1 = [1.0d0, 0.0d0, 0.0d0]
        v2 = [1.0d0, 1.0d0, 0.0d0]
        call vectors_angle_3D(v1, v2, angle)
        test_passed = (abs(angle - 45.0d0) < 1d-4)

        if (test_passed) then
            print *, "  [PASS] Test 4: 45 degree angle"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 4: 45 degree angle"
            print *, "    Expected: 45.0, Got: ", angle
            num_failed = num_failed + 1
        end if

    end subroutine test_vectors_angle_3D


    subroutine test_vectors_angle_2D(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(3) :: v1, v2
        real(kind=8) :: angle
        logical :: test_passed

        print *, ""
        print *, "Testing 2D angle calculation (vectors_angle_2D)..."

        ! Test 1: Parallel vectors in XY plane (0 degrees)
        num_tests = num_tests + 1
        v1 = [1.0d0, 0.0d0, 5.0d0]  ! Z component should be ignored
        v2 = [2.0d0, 0.0d0, -3.0d0]
        call vectors_angle_2D(v1, v2, angle)
        test_passed = (abs(angle - 0.0d0) < 1d-6)

        if (test_passed) then
            print *, "  [PASS] Test 1: Parallel vectors in XY plane (0 degrees)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Parallel vectors in XY plane"
            print *, "    Expected: 0.0, Got: ", angle
            num_failed = num_failed + 1
        end if

        ! Test 2: Perpendicular vectors in XY plane (90 degrees)
        num_tests = num_tests + 1
        v1 = [1.0d0, 0.0d0, 0.0d0]
        v2 = [0.0d0, 1.0d0, 0.0d0]
        call vectors_angle_2D(v1, v2, angle)
        test_passed = (abs(angle - 90.0d0) < 1d-6)

        if (test_passed) then
            print *, "  [PASS] Test 2: Perpendicular vectors in XY plane (90 degrees)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Perpendicular vectors in XY plane"
            print *, "    Expected: 90.0, Got: ", angle
            num_failed = num_failed + 1
        end if

        ! Test 3: Opposite vectors in XY plane (180 degrees)
        num_tests = num_tests + 1
        v1 = [1.0d0, 0.0d0, 10.0d0]
        v2 = [-1.0d0, 0.0d0, -5.0d0]
        call vectors_angle_2D(v1, v2, angle)
        test_passed = (abs(angle - 180.0d0) < 1d-6)

        if (test_passed) then
            print *, "  [PASS] Test 3: Opposite vectors in XY plane (180 degrees)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 3: Opposite vectors in XY plane"
            print *, "    Expected: 180.0, Got: ", angle
            num_failed = num_failed + 1
        end if

        ! Test 4: Zero magnitude in XY plane
        num_tests = num_tests + 1
        v1 = [0.0d0, 0.0d0, 5.0d0]  ! Zero in XY plane
        v2 = [1.0d0, 1.0d0, 0.0d0]
        call vectors_angle_2D(v1, v2, angle)
        test_passed = (abs(angle - 0.0d0) < 1d-6)

        if (test_passed) then
            print *, "  [PASS] Test 4: Zero magnitude vector in XY plane"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 4: Zero magnitude vector in XY plane"
            print *, "    Expected: 0.0, Got: ", angle
            num_failed = num_failed + 1
        end if

    end subroutine test_vectors_angle_2D


    subroutine test_rmsd(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:, :), allocatable :: coords1, coords2
        real(kind=8) :: rmsd_value
        logical :: test_passed
        integer :: n

        print *, ""
        print *, "Testing RMSD calculation..."

        ! Test 1: Identical coordinates (RMSD = 0)
        num_tests = num_tests + 1
        n = 3
        allocate(coords1(n, 3), coords2(n, 3))
        coords1(1, :) = [1.0d0, 2.0d0, 3.0d0]
        coords1(2, :) = [4.0d0, 5.0d0, 6.0d0]
        coords1(3, :) = [7.0d0, 8.0d0, 9.0d0]
        coords2 = coords1

        call rmsd(rmsd_value, n, coords1, coords2)
        test_passed = (abs(rmsd_value - 0.0d0) < 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 1: Identical coordinates (RMSD = 0)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Identical coordinates"
            print *, "    Expected: 0.0, Got: ", rmsd_value
            num_failed = num_failed + 1
        end if
        deallocate(coords1, coords2)

        ! Test 2: Simple translation (unit shift in x)
        num_tests = num_tests + 1
        n = 2
        allocate(coords1(n, 3), coords2(n, 3))
        coords1(1, :) = [0.0d0, 0.0d0, 0.0d0]
        coords1(2, :) = [1.0d0, 0.0d0, 0.0d0]
        coords2(1, :) = [1.0d0, 0.0d0, 0.0d0]
        coords2(2, :) = [2.0d0, 0.0d0, 0.0d0]

        call rmsd(rmsd_value, n, coords1, coords2)
        test_passed = (abs(rmsd_value - 1.0d0) < 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 2: Unit translation (RMSD = 1.0)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Unit translation"
            print *, "    Expected: 1.0, Got: ", rmsd_value
            num_failed = num_failed + 1
        end if
        deallocate(coords1, coords2)

        ! Test 3: Known RMSD value
        num_tests = num_tests + 1
        n = 1
        allocate(coords1(n, 3), coords2(n, 3))
        coords1(1, :) = [0.0d0, 0.0d0, 0.0d0]
        coords2(1, :) = [1.0d0, 1.0d0, 1.0d0]
        ! Expected RMSD = sqrt((1^2 + 1^2 + 1^2)/1) = sqrt(3)

        call rmsd(rmsd_value, n, coords1, coords2)
        test_passed = (abs(rmsd_value - sqrt(3.0d0)) < 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 3: Known RMSD (sqrt(3))"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 3: Known RMSD"
            print *, "    Expected: ", sqrt(3.0d0), ", Got: ", rmsd_value
            num_failed = num_failed + 1
        end if
        deallocate(coords1, coords2)

    end subroutine test_rmsd


    subroutine test_calculate_cog(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:, :), allocatable :: coords
        real(kind=8), dimension(3) :: cog, expected
        logical :: test_passed
        integer :: n

        print *, ""
        print *, "Testing center of geometry calculation..."

        ! Test 1: Single point
        num_tests = num_tests + 1
        n = 1
        allocate(coords(n, 3))
        coords(1, :) = [1.0d0, 2.0d0, 3.0d0]
        expected = [1.0d0, 2.0d0, 3.0d0]

        call calculate_cog(cog, coords, n)
        test_passed = assert_vector_equal(cog, expected, 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 1: Single point"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Single point"
            print *, "    Expected: ", expected
            print *, "    Got:      ", cog
            num_failed = num_failed + 1
        end if
        deallocate(coords)

        ! Test 2: Points at origin
        num_tests = num_tests + 1
        n = 3
        allocate(coords(n, 3))
        coords(1, :) = [1.0d0, 1.0d0, 1.0d0]
        coords(2, :) = [1.0d0, 1.0d0, 1.0d0]
        coords(3, :) = [1.0d0, 1.0d0, 1.0d0]
        expected = [1.0d0, 1.0d0, 1.0d0]

        call calculate_cog(cog, coords, n)
        test_passed = assert_vector_equal(cog, expected, 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 2: Points at (1,1,1)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Points at (1,1,1)"
            print *, "    Expected: ", expected
            print *, "    Got:      ", cog
            num_failed = num_failed + 1
        end if
        deallocate(coords)

        ! Test 3: Symmetric points around origin
        num_tests = num_tests + 1
        n = 2
        allocate(coords(n, 3))
        coords(1, :) = [1.0d0, 2.0d0, 3.0d0]
        coords(2, :) = [-1.0d0, -2.0d0, -3.0d0]
        expected = [0.0d0, 0.0d0, 0.0d0]

        call calculate_cog(cog, coords, n)
        test_passed = assert_vector_equal(cog, expected, 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 3: Symmetric points"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 3: Symmetric points"
            print *, "    Expected: ", expected
            print *, "    Got:      ", cog
            num_failed = num_failed + 1
        end if
        deallocate(coords)

        ! Test 4: General case
        num_tests = num_tests + 1
        n = 4
        allocate(coords(n, 3))
        coords(1, :) = [0.0d0, 0.0d0, 0.0d0]
        coords(2, :) = [2.0d0, 0.0d0, 0.0d0]
        coords(3, :) = [0.0d0, 2.0d0, 0.0d0]
        coords(4, :) = [0.0d0, 0.0d0, 2.0d0]
        expected = [0.5d0, 0.5d0, 0.5d0]

        call calculate_cog(cog, coords, n)
        test_passed = assert_vector_equal(cog, expected, 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 4: General case"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 4: General case"
            print *, "    Expected: ", expected
            print *, "    Got:      ", cog
            num_failed = num_failed + 1
        end if
        deallocate(coords)

    end subroutine test_calculate_cog


    subroutine test_calculate_distance(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(3) :: p1, p2
        real(kind=8) :: dist
        logical :: test_passed

        print *, ""
        print *, "Testing distance calculation..."

        ! Test 1: Distance to itself (0)
        num_tests = num_tests + 1
        p1 = [1.0d0, 2.0d0, 3.0d0]
        p2 = [1.0d0, 2.0d0, 3.0d0]

        call calculate_distance(p1, p2, dist)
        test_passed = (abs(dist - 0.0d0) < 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 1: Distance to itself (0)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Distance to itself"
            print *, "    Expected: 0.0, Got: ", dist
            num_failed = num_failed + 1
        end if

        ! Test 2: Unit distance along x-axis
        num_tests = num_tests + 1
        p1 = [0.0d0, 0.0d0, 0.0d0]
        p2 = [1.0d0, 0.0d0, 0.0d0]

        call calculate_distance(p1, p2, dist)
        test_passed = (abs(dist - 1.0d0) < 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 2: Unit distance along x-axis"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Unit distance along x-axis"
            print *, "    Expected: 1.0, Got: ", dist
            num_failed = num_failed + 1
        end if

        ! Test 3: 3D diagonal distance
        num_tests = num_tests + 1
        p1 = [0.0d0, 0.0d0, 0.0d0]
        p2 = [1.0d0, 1.0d0, 1.0d0]

        call calculate_distance(p1, p2, dist)
        test_passed = (abs(dist - sqrt(3.0d0)) < 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 3: 3D diagonal distance (sqrt(3))"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 3: 3D diagonal distance"
            print *, "    Expected: ", sqrt(3.0d0), ", Got: ", dist
            num_failed = num_failed + 1
        end if

        ! Test 4: Distance between arbitrary points
        num_tests = num_tests + 1
        p1 = [1.0d0, 2.0d0, 3.0d0]
        p2 = [4.0d0, 6.0d0, 8.0d0]
        ! Distance = sqrt((3^2 + 4^2 + 5^2)) = sqrt(50) = 5*sqrt(2)

        call calculate_distance(p1, p2, dist)
        test_passed = (abs(dist - sqrt(50.0d0)) < 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 4: Arbitrary points (sqrt(50))"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 4: Arbitrary points"
            print *, "    Expected: ", sqrt(50.0d0), ", Got: ", dist
            num_failed = num_failed + 1
        end if

    end subroutine test_calculate_distance


    !> Helper function: Check if two vectors are approximately equal
    logical function assert_vector_equal(v1, v2, tolerance)
        real(kind=8), dimension(3), intent(in) :: v1, v2
        real(kind=8), intent(in) :: tolerance
        integer :: i

        assert_vector_equal = .true.
        do i = 1, 3
            if (abs(v1(i) - v2(i)) > tolerance) then
                assert_vector_equal = .false.
                exit
            end if
        end do
    end function assert_vector_equal

end program test_maths
