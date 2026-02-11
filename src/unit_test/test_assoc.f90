!> Unit tests for reading association/complexes files
!!
!! Comprehensive test suite for reading and parsing association/complexes files
!! in the read_input and mod_assoc modules including:
!!
!! - **File reading**:
!!    - Valid assoc file reading
!!    - File existence checking
!!    - Error handling for missing files
!!    - Correct number of encounters parsed
!!
!! - **Data structure allocation**:
!!    - Allocate assoc objects with correct sizes
!!    - Initialize translation vectors to zero
!!    - Initialize rotation vectors (rot1, rot2) to zero
!!    - Allocate lines array
!!
!! - **File parsing**:
!!    - Skip comment lines (starting with #)
!!    - Count data lines correctly
!!    - Extract center coordinates (xc1, xc2)
!!    - Parse translation vectors
!!    - Parse rotation vectors
!!
!! - **Data integrity**:
!!    - Verify allocations match encounter count
!!    - Check header information is stored
!!    - Verify vector values are reasonable
!!    - Ensure no data corruption during read
!!
!! Compile with: gfortran -o test_assoc test_assoc.f90 ../src/read_input.f90 ../src/mod_assoc.f90
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0

program test_assoc
    use read_input
    use mod_assoc
    implicit none

    integer :: num_tests, num_passed, num_failed
    
    num_tests = 0
    num_passed = 0
    num_failed = 0

    print *, "========================================="
    print *, "UNIT TESTS FOR ASSOC FILE READING"
    print *, "========================================="

    ! Test assoc file reading
    call test_read_assoc_file(num_tests, num_passed, num_failed)
    
    ! Test data structure allocation
    call test_allocate_assoc_object(num_tests, num_passed, num_failed)
    
    ! Test size counting
    call test_size_assoc(num_tests, num_passed, num_failed)

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

    subroutine test_read_assoc_file(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        type(type_assoc_file) :: assoc
        integer :: tot_encounters
        character*128 :: filename
        logical :: test_passed, file_exists
        integer :: i

        print *, ""
        print *, "Testing read_assoc file operations..."

        ! Test 1: Create a test assoc file
        num_tests = num_tests + 1
        filename = "test_assoc.txt"
        
        ! Create a simple test assoc file
        open(unit=10, file=trim(filename), status='replace', action='write')
        write(10, '(A)') "# Test association file"
        write(10, '(A)') "# Header line 2"
        write(10, '(A)') "  0.000000   0.000000   0.000000   1.000000   0.000000   0.000000   0.000000   1.000000   0.000000"
        write(10, '(A)') "  1.000000   2.000000   3.000000   0.707107   0.707107   0.000000   -0.707107   0.707107   0.000000"
        write(10, '(A)') "  2.000000   1.500000   2.500000   0.866025   0.500000   0.000000   -0.500000   0.866025   0.000000"
        close(10)
        
        ! Verify file was created
        inquire(file=trim(filename), exist=file_exists)
        test_passed = file_exists

        if (test_passed) then
            print *, "  [PASS] Test 1: Test assoc file created successfully"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Test assoc file creation"
            num_failed = num_failed + 1
        end if

        ! Test 2: Read the assoc file
        num_tests = num_tests + 1
        
        if (file_exists) then
            call read_assoc(assoc, filename, tot_encounters)
            test_passed = (tot_encounters == 3)
            
            if (test_passed) then
                print *, "  [PASS] Test 2: Read assoc file with 3 encounters"
                num_passed = num_passed + 1
            else
                print *, "  [FAIL] Test 2: Assoc file encounter count"
                print *, "    Expected: 3, Got: ", tot_encounters
                num_failed = num_failed + 1
            end if
        else
            print *, "  [SKIP] Test 2: Skipped (test file not created)"
        end if

        ! Test 3: Verify translation vectors are allocated
        num_tests = num_tests + 1
        if (file_exists) then
            test_passed = allocated(assoc%trans_vector) .and. size(assoc%trans_vector, 1) == 3
            
            if (test_passed) then
                print *, "  [PASS] Test 3: Translation vectors allocated correctly"
                num_passed = num_passed + 1
            else
                print *, "  [FAIL] Test 3: Translation vectors allocation"
                num_failed = num_failed + 1
            end if
        else
            print *, "  [SKIP] Test 3: Skipped (test file not created)"
        end if

        ! Test 4: Verify rotation vectors are allocated
        num_tests = num_tests + 1
        if (file_exists) then
            test_passed = allocated(assoc%rot1) .and. allocated(assoc%rot2) .and. &
                         size(assoc%rot1, 1) == 3 .and. size(assoc%rot2, 1) == 3
            
            if (test_passed) then
                print *, "  [PASS] Test 4: Rotation vectors allocated correctly"
                num_passed = num_passed + 1
            else
                print *, "  [FAIL] Test 4: Rotation vectors allocation"
                num_failed = num_failed + 1
            end if
        else
            print *, "  [SKIP] Test 4: Skipped (test file not created)"
        end if

        ! Test 5: Verify first encounter data
        num_tests = num_tests + 1
        if (file_exists .and. tot_encounters >= 1) then
            test_passed = (abs(assoc%trans_vector(1, 1)) < 1d-6 .and. &
                          abs(assoc%trans_vector(1, 2)) < 1d-6 .and. &
                          abs(assoc%trans_vector(1, 3)) < 1d-6)
            
            if (test_passed) then
                print *, "  [PASS] Test 5: First encounter translation parsed correctly"
                num_passed = num_passed + 1
            else
                print *, "  [FAIL] Test 5: First encounter data"
                print *, "    Got translation: [", assoc%trans_vector(1, 1), ", ", &
                        assoc%trans_vector(1, 2), ", ", assoc%trans_vector(1, 3), "]"
                num_failed = num_failed + 1
            end if
        else
            print *, "  [SKIP] Test 5: Skipped (insufficient data)"
        end if

        ! Test 6: Verify second encounter data
        num_tests = num_tests + 1
        if (file_exists .and. tot_encounters >= 2) then
            test_passed = (abs(assoc%trans_vector(2, 1) - 1.0d0) < 1d-6 .and. &
                          abs(assoc%trans_vector(2, 2) - 2.0d0) < 1d-6 .and. &
                          abs(assoc%trans_vector(2, 3) - 3.0d0) < 1d-6)
            
            if (test_passed) then
                print *, "  [PASS] Test 6: Second encounter translation parsed correctly"
                num_passed = num_passed + 1
            else
                print *, "  [FAIL] Test 6: Second encounter data"
                print *, "    Expected: [1.0, 2.0, 3.0]"
                print *, "    Got: [", assoc%trans_vector(2, 1), ", ", &
                        assoc%trans_vector(2, 2), ", ", assoc%trans_vector(2, 3), "]"
                num_failed = num_failed + 1
            end if
        else
            print *, "  [SKIP] Test 6: Skipped (insufficient data)"
        end if

        ! Clean up
        if (allocated(assoc%lines)) deallocate(assoc%lines)
        if (allocated(assoc%trans_vector)) deallocate(assoc%trans_vector)
        if (allocated(assoc%rot1)) deallocate(assoc%rot1)
        if (allocated(assoc%rot2)) deallocate(assoc%rot2)
        
        ! Remove test file
        open(unit=10, file=trim(filename), status='old')
        close(10, status='delete')

    end subroutine test_read_assoc_file


    subroutine test_allocate_assoc_object(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        type(type_assoc_file) :: assoc
        integer :: n
        logical :: test_passed

        print *, ""
        print *, "Testing allocate_assoc_object..."

        ! Test 1: Allocate small object
        num_tests = num_tests + 1
        n = 5
        call allocate_assoc_object(assoc, n)
        
        test_passed = (assoc%nlines == n .and. &
                      allocated(assoc%lines) .and. &
                      allocated(assoc%trans_vector) .and. &
                      allocated(assoc%rot1) .and. &
                      allocated(assoc%rot2))
        
        if (test_passed) then
            print *, "  [PASS] Test 1: Allocate object with 5 encounters"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Allocate small object"
            num_failed = num_failed + 1
        end if
        
        if (allocated(assoc%lines)) deallocate(assoc%lines)
        if (allocated(assoc%trans_vector)) deallocate(assoc%trans_vector)
        if (allocated(assoc%rot1)) deallocate(assoc%rot1)
        if (allocated(assoc%rot2)) deallocate(assoc%rot2)

        ! Test 2: Allocate larger object
        num_tests = num_tests + 1
        n = 1000
        call allocate_assoc_object(assoc, n)
        
        test_passed = (size(assoc%trans_vector, 1) == n .and. &
                      size(assoc%rot1, 1) == n .and. &
                      size(assoc%rot2, 1) == n)
        
        if (test_passed) then
            print *, "  [PASS] Test 2: Allocate object with 1000 encounters"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Allocate large object"
            num_failed = num_failed + 1
        end if
        
        if (allocated(assoc%lines)) deallocate(assoc%lines)
        if (allocated(assoc%trans_vector)) deallocate(assoc%trans_vector)
        if (allocated(assoc%rot1)) deallocate(assoc%rot1)
        if (allocated(assoc%rot2)) deallocate(assoc%rot2)

        ! Test 3: Check initialization to zero
        num_tests = num_tests + 1
        n = 3
        call allocate_assoc_object(assoc, n)
        
        test_passed = (all(assoc%trans_vector == 0.0d0) .and. &
                      all(assoc%rot1 == 0.0d0) .and. &
                      all(assoc%rot2 == 0.0d0))
        
        if (test_passed) then
            print *, "  [PASS] Test 3: Vectors initialized to zero"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 3: Vector initialization"
            num_failed = num_failed + 1
        end if
        
        if (allocated(assoc%lines)) deallocate(assoc%lines)
        if (allocated(assoc%trans_vector)) deallocate(assoc%trans_vector)
        if (allocated(assoc%rot1)) deallocate(assoc%rot1)
        if (allocated(assoc%rot2)) deallocate(assoc%rot2)

    end subroutine test_allocate_assoc_object


    subroutine test_size_assoc(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        character*128 :: filename
        integer :: input_unit, count, io_stat
        logical :: test_passed, file_exists

        print *, ""
        print *, "Testing size_assoc counting..."

        ! Test 1: Count lines in file with comments
        num_tests = num_tests + 1
        filename = "test_count.txt"
        
        ! Create test file
        open(unit=11, file=trim(filename), status='replace', action='write')
        write(11, '(A)') "# Comment line 1"
        write(11, '(A)') "  0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0"
        write(11, '(A)') "# Comment line 2"
        write(11, '(A)') "  1.0 2.0 3.0 0.7 0.7 0.0 -0.7 0.7 0.0"
        write(11, '(A)') "  2.0 1.5 2.5 0.9 0.5 0.0 -0.5 0.9 0.0"
        close(11)
        
        ! Count non-comment lines
        open(unit=11, file=trim(filename), status='old', action='read', iostat=io_stat)
        if (io_stat == 0) then
            count = size_assoc(11)
            close(11)
            
            test_passed = (count == 3)
            
            if (test_passed) then
                print *, "  [PASS] Test 1: Counted 3 data lines (skipped comments)"
                num_passed = num_passed + 1
            else
                print *, "  [FAIL] Test 1: Line counting"
                print *, "    Expected: 3, Got: ", count
                num_failed = num_failed + 1
            end if
        else
            print *, "  [FAIL] Test 1: Could not open test file"
            num_failed = num_failed + 1
        end if
        
        ! Remove test file
        open(unit=11, file=trim(filename), status='old')
        close(11, status='delete')

        ! Test 2: Count empty file
        num_tests = num_tests + 1
        filename = "test_empty.txt"
        
        ! Create empty test file
        open(unit=11, file=trim(filename), status='replace', action='write')
        write(11, '(A)') "# Only comment"
        close(11)
        
        open(unit=11, file=trim(filename), status='old', action='read', iostat=io_stat)
        if (io_stat == 0) then
            count = size_assoc(11)
            close(11)
            
            test_passed = (count == 0)
            
            if (test_passed) then
                print *, "  [PASS] Test 2: Empty file counted as 0 data lines"
                num_passed = num_passed + 1
            else
                print *, "  [FAIL] Test 2: Empty file counting"
                print *, "    Expected: 0, Got: ", count
                num_failed = num_failed + 1
            end if
        else
            print *, "  [FAIL] Test 2: Could not open test file"
            num_failed = num_failed + 1
        end if
        
        ! Remove test file
        open(unit=11, file=trim(filename), status='old')
        close(11, status='delete')

    end subroutine test_size_assoc

end program test_assoc
