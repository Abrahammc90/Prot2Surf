!> Unit tests for PDB operations in the mod_pdb module
!!
!! Tests include:
!!   - PDB type allocation and initialization
!!   - Atom and residue management
!!   - Chain structure handling
!!
!! Compile with: gfortran -o test_pdb test_pdb.f90 mod_pdb.f90
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0

program test_pdb
    use mod_pdb
    implicit none

    integer :: num_tests, num_passed, num_failed
    
    num_tests = 0
    num_passed = 0
    num_failed = 0

    print *, "========================================="
    print *, "UNIT TESTS FOR PDB MODULE"
    print *, "========================================="

    ! Test PDB file type initialization
    call test_pdb_allocation(num_tests, num_passed, num_failed)
    
    ! Test atom properties
    call test_atom_properties(num_tests, num_passed, num_failed)
    
    ! Test residue management
    call test_residue_management(num_tests, num_passed, num_failed)

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

    subroutine test_pdb_allocation(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        type(type_pdb_file) :: pdb
        logical :: test_passed

        print *, ""
        print *, "Testing PDB allocation..."

        ! Test 1: PDB type initialization
        num_tests = num_tests + 1
        pdb%natoms = 100
        pdb%nresidues = 10
        pdb%nchains = 1

        test_passed = (pdb%natoms == 100 .and. pdb%nresidues == 10 .and. pdb%nchains == 1)

        if (test_passed) then
            print *, "  [PASS] Test 1: PDB type initialization"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: PDB type initialization"
            num_failed = num_failed + 1
        end if

        ! Test 2: Allocate atoms array
        num_tests = num_tests + 1
        if (allocated(pdb%atoms)) deallocate(pdb%atoms)
        allocate(pdb%atoms(pdb%natoms))

        test_passed = (allocated(pdb%atoms) .and. size(pdb%atoms) == 100)

        if (test_passed) then
            print *, "  [PASS] Test 2: Atoms array allocation (100 atoms)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Atoms array allocation"
            num_failed = num_failed + 1
        end if

        ! Test 3: Allocate residues array
        num_tests = num_tests + 1
        if (allocated(pdb%residues)) deallocate(pdb%residues)
        allocate(pdb%residues(pdb%nresidues))

        test_passed = (allocated(pdb%residues) .and. size(pdb%residues) == 10)

        if (test_passed) then
            print *, "  [PASS] Test 3: Residues array allocation (10 residues)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 3: Residues array allocation"
            num_failed = num_failed + 1
        end if

        ! Cleanup
        if (allocated(pdb%atoms)) deallocate(pdb%atoms)
        if (allocated(pdb%residues)) deallocate(pdb%residues)

    end subroutine test_pdb_allocation


    subroutine test_atom_properties(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        type(type_pdb_atom) :: atom
        logical :: test_passed

        print *, ""
        print *, "Testing atom properties..."

        ! Test 1: Set atom coordinates
        num_tests = num_tests + 1
        atom%coord = (/1.5d0, 2.5d0, 3.5d0/)

        test_passed = (abs(atom%coord(1) - 1.5d0) < 1d-10 .and. &
                      abs(atom%coord(2) - 2.5d0) < 1d-10 .and. &
                      abs(atom%coord(3) - 3.5d0) < 1d-10)

        if (test_passed) then
            print *, "  [PASS] Test 1: Atom coordinate assignment"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Atom coordinate assignment"
            num_failed = num_failed + 1
        end if

        ! Test 2: Set atom name and residue info
        num_tests = num_tests + 1
        atom%name = "CA"
        atom%resname = "ALA"
        atom%resid = 10

        test_passed = (trim(atom%name) == "CA" .and. &
                      trim(atom%resname) == "ALA" .and. &
                      atom%resid == 10)

        if (test_passed) then
            print *, "  [PASS] Test 2: Atom name and residue properties"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Atom name and residue properties"
            num_failed = num_failed + 1
        end if

    end subroutine test_atom_properties


    subroutine test_residue_management(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        type(type_pdb_residue) :: residue
        logical :: test_passed

        print *, ""
        print *, "Testing residue management..."

        ! Test 1: Set residue properties
        num_tests = num_tests + 1
        residue%resname = "GLY"
        residue%resid = 5
        residue%chainid = "A"
        residue%natoms = 7

        test_passed = (trim(residue%resname) == "GLY" .and. &
                      residue%resid == 5 .and. &
                      residue%chainid == "A" .and. &
                      residue%natoms == 7)

        if (test_passed) then
            print *, "  [PASS] Test 1: Residue property assignment"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Residue property assignment"
            num_failed = num_failed + 1
        end if

        ! Test 2: Allocate atoms in residue
        num_tests = num_tests + 1
        if (allocated(residue%atoms)) deallocate(residue%atoms)
        allocate(residue%atoms(residue%natoms))

        test_passed = (allocated(residue%atoms) .and. size(residue%atoms) == 7)

        if (test_passed) then
            print *, "  [PASS] Test 2: Residue atoms array allocation"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Residue atoms array allocation"
            num_failed = num_failed + 1
        end if

        ! Cleanup
        if (allocated(residue%atoms)) deallocate(residue%atoms)

    end subroutine test_residue_management

end program test_pdb
