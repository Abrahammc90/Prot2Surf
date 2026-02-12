!> Unit tests for matrix operations in the mod_matrix module
!!
!! Tests include:
!!   - Write and read matrix operations
!!   - Write and read array operations
!!   - Matrix Z-coordinate calculation
!!   - Matrix atoms distance calculation 
!!   - Matrix angle calculation
!!   - RMSD matrix calculation
!!
!! Compile with: gfortran -fopenmp -o test_matrix test_matrix.f90 ../src/mod_matrix.f90 ../src/maths.f90 ../src/mod_pdb.f90
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0

program test_matrix
    use mod_matrix
    use maths
    implicit none

    integer :: num_tests, num_passed, num_failed
    
    num_tests = 0
    num_passed = 0
    num_failed = 0

    print *, "========================================="
    print *, "UNIT TESTS FOR MATRIX MODULE"
    print *, "========================================="

    ! Test matrix write and read
    call test_write_read_matrix(num_tests, num_passed, num_failed)
    
    ! Test array write and read
    call test_write_read_array(num_tests, num_passed, num_failed)
    
    ! Test matrix_z_coord
    call test_matrix_z_coord(num_tests, num_passed, num_failed)
    
    ! Test matrix_rmsd
    call test_matrix_rmsd_calc(num_tests, num_passed, num_failed)

    ! Test matrix_atoms_dist
    call test_matrix_atoms_dist(num_tests, num_passed, num_failed)

    ! Test matrix_angle
    call test_matrix_angle_calc(num_tests, num_passed, num_failed)

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

    subroutine test_write_read_matrix(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:, :), allocatable :: matrix_orig, matrix_read
        character*128 :: filename
        integer :: n, i, j
        logical :: test_passed

        print *, ""
        print *, "Testing write_matrix and read_matrix..."

        ! Test 1: Write and read small matrix
        num_tests = num_tests + 1
        n = 3
        allocate(matrix_orig(n, n))
        filename = "test_matrix_small.dat"
        
        ! Create test matrix
        matrix_orig(1, :) = [1.0d0, 2.5d0, 3.7d0]
        matrix_orig(2, :) = [2.5d0, 0.0d0, 1.2d0]
        matrix_orig(3, :) = [3.7d0, 1.2d0, 0.0d0]
        
        call write_matrix(matrix_orig, filename)
        call read_matrix(matrix_read, n, filename)
        
        test_passed = .true.
        do i = 1, n
            do j = 1, n
                if (abs(matrix_orig(i, j) - matrix_read(i, j)) > 1d-3) then
                    test_passed = .false.
                end if
            end do
        end do
        
        if (test_passed) then
            print *, "  [PASS] Test 1: Write and read small matrix"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Write and read small matrix"
            num_failed = num_failed + 1
        end if
        deallocate(matrix_orig, matrix_read)

        ! Test 2: Verify matrix symmetry preserved
        num_tests = num_tests + 1
        n = 4
        allocate(matrix_orig(n, n))
        filename = "test_matrix_symmetric.dat"
        
        do i = 1, n
            matrix_orig(i, i) = 0.0d0
            do j = i+1, n
                matrix_orig(i, j) = abs(i - j) * 1.234d0
                matrix_orig(j, i) = matrix_orig(i, j)
            end do
        end do
        
        call write_matrix(matrix_orig, filename)
        call read_matrix(matrix_read, n, filename)
        
        test_passed = .true.
        do i = 1, n
            do j = i+1, n
                if (abs(matrix_read(i, j) - matrix_read(j, i)) > 1d-3) then
                    test_passed = .false.
                end if
            end do
        end do
        
        if (test_passed) then
            print *, "  [PASS] Test 2: Matrix symmetry preserved"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Matrix symmetry preserved"
            num_failed = num_failed + 1
        end if
        deallocate(matrix_orig, matrix_read)

    end subroutine test_write_read_matrix


    subroutine test_write_read_array(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:), allocatable :: array_orig, array_read
        character*128 :: filename
        integer :: n, i
        logical :: test_passed

        print *, ""
        print *, "Testing write_array and read_array..."

        ! Test 1: Write and read array
        num_tests = num_tests + 1
        n = 5
        allocate(array_orig(n))
        filename = "test_array.dat"
        
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

        ! Test 2: Array with negative and zero values
        num_tests = num_tests + 1
        n = 4
        allocate(array_orig(n))
        filename = "test_array_neg.dat"
        
        array_orig = [-1.5d0, 0.0d0, 3.3d0, -7.2d0]
        
        call write_array(array_orig, filename)
        call read_array(array_read, n, filename)
        
        test_passed = .true.
        do i = 1, n
            if (abs(array_orig(i) - array_read(i)) > 1d-3) then
                test_passed = .false.
            end if
        end do
        
        if (test_passed) then
            print *, "  [PASS] Test 2: Array with negative and zero values"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Array with negative and zero values"
            num_failed = num_failed + 1
        end if
        deallocate(array_orig, array_read)

    end subroutine test_write_read_array


    subroutine test_matrix_z_coord(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:, :), allocatable :: matrix, solute_crds, trans_vector, rot1, rot2
        real(kind=8), dimension(:), allocatable :: array
        real(kind=8), dimension(3) :: xc1, xc2
        integer :: n, nb_atoms, i, j
        logical :: test_passed

        print *, ""
        print *, "Testing matrix_z_coord..."

        ! Test 1: Identity transformation - Z values should remain same
        num_tests = num_tests + 1
        n = 2
        nb_atoms = 1
        allocate(matrix(n, n), array(n))
        allocate(solute_crds(nb_atoms, 3), trans_vector(n, 3), rot1(n, 3), rot2(n, 3))
        
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
        
        call matrix_z_coord(matrix, array, n, nb_atoms, xc1, xc2, &
                           trans_vector, rot1, rot2, solute_crds)
        
        test_passed = (abs(array(1) - 5.0d0) < 1d-8 .and. abs(array(2) - 5.0d0) < 1d-8 .and. &
                      abs(matrix(1, 2)) < 1d-8 .and. abs(matrix(1, 1)) < 1d-8)
        
        if (test_passed) then
            print *, "  [PASS] Test 1: Identity transformation"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Identity transformation"
            print *, "    Z values: ", array(1), array(2)
            print *, "    Matrix diff: ", matrix(1, 2)
            num_failed = num_failed + 1
        end if
        deallocate(matrix, array, solute_crds, trans_vector, rot1, rot2)

        ! Test 2: Matrix diagonal should be zero
        num_tests = num_tests + 1
        n = 3
        nb_atoms = 1
        allocate(matrix(n, n), array(n))
        allocate(solute_crds(nb_atoms, 3), trans_vector(n, 3), rot1(n, 3), rot2(n, 3))
        
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        solute_crds(1, :) = [0.0d0, 0.0d0, 1.0d0]
        
        do i = 1, n
            trans_vector(i, :) = [0.0d0, 0.0d0, real(i-1, 8)]
            rot1(i, :) = [1.0d0, 0.0d0, 0.0d0]
            rot2(i, :) = [0.0d0, 1.0d0, 0.0d0]
        end do
        
        call matrix_z_coord(matrix, array, n, nb_atoms, xc1, xc2, &
                           trans_vector, rot1, rot2, solute_crds)
        
        test_passed = .true.
        do i = 1, n
            if (abs(matrix(i, i)) > 1d-8) then
                test_passed = .false.
            end if
        end do
        
        if (test_passed) then
            print *, "  [PASS] Test 2: Matrix diagonal is zero"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Matrix diagonal is zero"
            num_failed = num_failed + 1
        end if
        deallocate(matrix, array, solute_crds, trans_vector, rot1, rot2)

    end subroutine test_matrix_z_coord


    subroutine test_matrix_rmsd_calc(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:, :), allocatable :: matrix, coords, trans_vector, rot1, rot2
        real(kind=8), dimension(3) :: xc1, xc2
        integer :: n, nb_atoms, i, j
        logical :: test_passed

        print *, ""
        print *, "Testing matrix_rmsd..."

        ! Test 1: Identical transformations should give RMSD = 0
        num_tests = num_tests + 1
        n = 2
        nb_atoms = 2
        allocate(matrix(n, n), coords(nb_atoms, 3), trans_vector(n, 3), rot1(n, 3), rot2(n, 3))
        
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        coords(1, :) = [1.0d0, 2.0d0, 3.0d0]
        coords(2, :) = [4.0d0, 5.0d0, 6.0d0]
        
        ! Same transformation for both
        do i = 1, n
            trans_vector(i, :) = [0.0d0, 0.0d0, 0.0d0]
            rot1(i, :) = [1.0d0, 0.0d0, 0.0d0]
            rot2(i, :) = [0.0d0, 1.0d0, 0.0d0]
        end do
        
        call matrix_rmsd(matrix, n, nb_atoms, xc1, xc2, trans_vector, rot1, rot2, coords)
        
        test_passed = (abs(matrix(1, 2)) < 1d-8)
        
        if (test_passed) then
            print *, "  [PASS] Test 1: Identical transformations (RMSD = 0)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Identical transformations"
            print *, "    RMSD: ", matrix(1, 2)
            num_failed = num_failed + 1
        end if
        deallocate(matrix, coords, trans_vector, rot1, rot2)

        ! Test 2: Matrix should be symmetric
        num_tests = num_tests + 1
        n = 3
        nb_atoms = 1
        allocate(matrix(n, n), coords(nb_atoms, 3), trans_vector(n, 3), rot1(n, 3), rot2(n, 3))
        
        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]
        coords(1, :) = [1.0d0, 0.0d0, 0.0d0]
        
        do i = 1, n
            trans_vector(i, :) = [real(i-1, 8), 0.0d0, 0.0d0]
            rot1(i, :) = [1.0d0, 0.0d0, 0.0d0]
            rot2(i, :) = [0.0d0, 1.0d0, 0.0d0]
        end do
        
        call matrix_rmsd(matrix, n, nb_atoms, xc1, xc2, trans_vector, rot1, rot2, coords)
        
        test_passed = .true.
        do i = 1, n
            do j = 1, n
                if (abs(matrix(i, j) - matrix(j, i)) > 1d-8) then
                    test_passed = .false.
                end if
            end do
        end do
        
        if (test_passed) then
            print *, "  [PASS] Test 2: RMSD matrix is symmetric"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: RMSD matrix is symmetric"
            num_failed = num_failed + 1
        end if
        deallocate(matrix, coords, trans_vector, rot1, rot2)

    end subroutine test_matrix_rmsd_calc


    subroutine test_matrix_atoms_dist(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:, :), allocatable :: matrix, solute1_crds, solute2_crds, trans_vector, rot1, rot2
        real(kind=8), dimension(:), allocatable :: array
        real(kind=8), dimension(3) :: xc1, xc2
        integer :: n, nb_atoms
        logical :: test_passed

        print *, ""
        print *, "Testing matrix_atoms_dist..."

        ! Test 1: Simple translations produce expected min distances
        num_tests = num_tests + 1
        n = 2
        nb_atoms = 1
        allocate(matrix(n, n), array(n))
        allocate(solute1_crds(1, 3), solute2_crds(nb_atoms, 3))
        allocate(trans_vector(n, 3), rot1(n, 3), rot2(n, 3))

        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]

        solute1_crds(1, :) = [0.0d0, 0.0d0, 0.0d0]
        solute2_crds(1, :) = [1.0d0, 0.0d0, 0.0d0]

        trans_vector(1, :) = [0.0d0, 0.0d0, 0.0d0]
        trans_vector(2, :) = [1.0d0, 0.0d0, 0.0d0]
        rot1(1, :) = [1.0d0, 0.0d0, 0.0d0]
        rot1(2, :) = [1.0d0, 0.0d0, 0.0d0]
        rot2(1, :) = [0.0d0, 1.0d0, 0.0d0]
        rot2(2, :) = [0.0d0, 1.0d0, 0.0d0]

        call matrix_atoms_dist(matrix, array, n, nb_atoms, xc1, xc2, trans_vector, rot1, rot2, solute1_crds, solute2_crds)

        test_passed = (abs(array(1) - 1.0d0) < 1d-8 .and. abs(array(2) - 2.0d0) < 1d-8 .and. &
                      abs(matrix(1, 2) - 1.0d0) < 1d-8 .and. abs(matrix(2, 1) - 1.0d0) < 1d-8)

        if (test_passed) then
            print *, "  [PASS] Test 1: Minimum distance per encounter"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Minimum distance per encounter"
            print *, "    array: ", array(1), array(2)
            print *, "    matrix(1,2): ", matrix(1, 2)
            num_failed = num_failed + 1
        end if

        deallocate(matrix, array, solute1_crds, solute2_crds, trans_vector, rot1, rot2)

    end subroutine test_matrix_atoms_dist


    subroutine test_matrix_angle_calc(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        real(kind=8), dimension(:, :), allocatable :: matrix, trans_vector, rot1, rot2
        real(kind=8), dimension(:), allocatable :: array
        real(kind=8), dimension(3) :: xc1, xc2
        real(kind=8), dimension(3) :: point1a, point1b, point2a, point2b
        integer :: n, nb_atoms, dimensions
        logical :: test_passed

        print *, ""
        print *, "Testing matrix_angle..."

        ! Test 1: 0 and 90 degree angles in 2D
        num_tests = num_tests + 1
        n = 2
        nb_atoms = 2
        dimensions = 2
        allocate(matrix(n, n), array(n))
        allocate(trans_vector(n, 3), rot1(n, 3), rot2(n, 3))

        xc1 = [0.0d0, 0.0d0, 0.0d0]
        xc2 = [0.0d0, 0.0d0, 0.0d0]

        point1a = [0.0d0, 0.0d0, 0.0d0]
        point1b = [1.0d0, 0.0d0, 0.0d0]
        point2a = [0.0d0, 0.0d0, 0.0d0]
        point2b = [1.0d0, 0.0d0, 0.0d0]

        trans_vector(1, :) = [0.0d0, 0.0d0, 0.0d0]
        rot1(1, :) = [1.0d0, 0.0d0, 0.0d0]
        rot2(1, :) = [0.0d0, 1.0d0, 0.0d0]

        trans_vector(2, :) = [0.0d0, 0.0d0, 0.0d0]
        rot1(2, :) = [0.0d0, 1.0d0, 0.0d0]
        rot2(2, :) = [-1.0d0, 0.0d0, 0.0d0]

        call matrix_angle(matrix, array, n, nb_atoms, xc1, xc2, trans_vector, rot1, rot2, &
                          point1a, point1b, point2a, point2b, dimensions)

        test_passed = (abs(array(1) - 0.0d0) < 1d-6 .and. abs(array(2) - 90.0d0) < 1d-6 .and. &
                      abs(matrix(1, 2) - 90.0d0) < 1d-6 .and. abs(matrix(2, 1) - 90.0d0) < 1d-6)

        if (test_passed) then
            print *, "  [PASS] Test 1: 2D angles (0 and 90 degrees)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: 2D angles (0 and 90 degrees)"
            print *, "    array: ", array(1), array(2)
            print *, "    matrix(1,2): ", matrix(1, 2)
            num_failed = num_failed + 1
        end if

        deallocate(matrix, array, trans_vector, rot1, rot2)

    end subroutine test_matrix_angle_calc

end program test_matrix
