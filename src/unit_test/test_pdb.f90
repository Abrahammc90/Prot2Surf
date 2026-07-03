!> Unit tests for PDB operations in the mod_pdb module
!!
!! Tests include:
!!   - PDB type allocation and initialization
!!   - Atom and residue management
!!   - Residue structure handling
!!   - Reading PDB files and filling the PDB data structure
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

    ! Test chain management
    call test_chain_management(num_tests, num_passed, num_failed)

    ! Test mod_pdb routines on a small PDB file
    call test_mod_pdb_routines(num_tests, num_passed, num_failed)

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
    !> @brief Unit test for PDB allocation and initialization routines.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2024-06-09
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
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
        allocate(pdb%residues(pdb%nresidues))
        allocate(pdb%chains(pdb%nchains))

        test_passed = (allocated(pdb%residues) .and. size(pdb%residues) == 10)

        if (test_passed) then
            print *, "  [PASS] Test 3: Residues array allocation (10 residues)"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 3: Residues array allocation"
            num_failed = num_failed + 1
        end if

        ! Cleanup
        deallocate(pdb%atoms)
        deallocate(pdb%residues)
        deallocate(pdb%chains)

    end subroutine test_pdb_allocation


    !*******************************************************************************
    !> @brief Unit test for atom property assignment in PDB structure.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2024-06-09
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
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


    !*******************************************************************************
    !> @brief Unit test for residue management in PDB structure.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2024-06-09
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
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


    !*******************************************************************************
    !> @brief Unit test for chain management in PDB structure.
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2024-06-09
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
    subroutine test_chain_management(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        type(type_pdb_chain) :: chain
        logical :: test_passed

        print *, ""
        print *, "Testing chain management..."

        ! Test 1: Set chain properties
        num_tests = num_tests + 1
        chain%chainid = "A"
        chain%natoms = 5
        chain%nresidues = 2

        test_passed = (chain%chainid == "A" .and. chain%natoms == 5 .and. chain%nresidues == 2)

        if (test_passed) then
            print *, "  [PASS] Test 1: Chain property assignment"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: Chain property assignment"
            num_failed = num_failed + 1
        end if

        ! Test 2: Allocate atoms and residues in chain
        num_tests = num_tests + 1
        if (allocated(chain%atoms)) deallocate(chain%atoms)
        if (allocated(chain%residues)) deallocate(chain%residues)
        allocate(chain%atoms(chain%natoms))
        allocate(chain%residues(chain%nresidues))

        test_passed = (allocated(chain%atoms) .and. size(chain%atoms) == 5 .and. &
                      allocated(chain%residues) .and. size(chain%residues) == 2)

        if (test_passed) then
            print *, "  [PASS] Test 2: Chain arrays allocation"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: Chain arrays allocation"
            num_failed = num_failed + 1
        end if

        ! Cleanup
        if (allocated(chain%atoms)) deallocate(chain%atoms)
        if (allocated(chain%residues)) deallocate(chain%residues)

    end subroutine test_chain_management


    !*******************************************************************************
    !> @brief Unit test for mod_pdb routines (count/allocate/fill).
    !>
    !> @author Abraham
    !> @version 1.0
    !> @date 2024-06-09
    !>
    !> @param[inout] num_tests   Number of tests executed (incremented)
    !> @param[inout] num_passed  Number of tests passed (incremented)
    !> @param[inout] num_failed  Number of tests failed (incremented)
    !*******************************************************************************
    subroutine test_mod_pdb_routines(num_tests, num_passed, num_failed)
        integer, intent(inout) :: num_tests, num_passed, num_failed
        type(type_pdb_file) :: pdb
        integer :: pdb_unit, natoms, nresidues, nchains
        character(len=128) :: filename
        logical :: test_passed

        print *, ""
        print *, "Testing mod_pdb routines (count/allocate/fill)..."

        filename = "test_pdb_temp.pdb"

        open(newunit=pdb_unit, file=filename, status="replace", action="write")
        write(pdb_unit,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)') "ATOM  ", 1, "CA  ", "ALA", "A", 1, 0.0d0, 0.0d0, 0.0d0
        write(pdb_unit,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)') "ATOM  ", 2, "CB  ", "ALA", "A", 1, 1.0d0, 0.0d0, 0.0d0
        write(pdb_unit,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)') "ATOM  ", 3, "CA  ", "GLY", "A", 2, 0.0d0, 1.0d0, 0.0d0
        write(pdb_unit,'(A3)') "TER"
        write(pdb_unit,'(A3)') "END"
        close(pdb_unit)

        open(newunit=pdb_unit, file=filename, status="old", action="read")
        natoms = count_atoms(pdb_unit)
        nresidues = count_residues(pdb_unit)
        nchains = count_chains(pdb_unit)

        num_tests = num_tests + 1
        test_passed = (natoms == 3 .and. nresidues == 2 .and. nchains == 1)
        if (test_passed) then
            print *, "  [PASS] Test 1: count_atoms/count_residues/count_chains"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 1: count_atoms/count_residues/count_chains"
            print *, "    natoms=", natoms, " nresidues=", nresidues, " nchains=", nchains
            num_failed = num_failed + 1
        end if

        call allocate_pdb_object(pdb, pdb_unit, natoms, nresidues, nchains)
        call fill_pdb_object(pdb, pdb_unit)

        num_tests = num_tests + 1
        test_passed = (pdb%natoms == 3 .and. pdb%nresidues == 2 .and. pdb%nchains == 1 .and. &
                      allocated(pdb%atoms) .and. allocated(pdb%residues) .and. allocated(pdb%chains))
        if (test_passed) then
            print *, "  [PASS] Test 2: allocate_pdb_object/fill_pdb_object"
            num_passed = num_passed + 1
        else
            print *, "  [FAIL] Test 2: allocate_pdb_object/fill_pdb_object"
            num_failed = num_failed + 1
        end if

        close(pdb_unit)
        if (allocated(pdb%atoms)) deallocate(pdb%atoms)
        if (allocated(pdb%residues)) deallocate(pdb%residues)
        if (allocated(pdb%chains)) deallocate(pdb%chains)
    end subroutine test_mod_pdb_routines

end program test_pdb
