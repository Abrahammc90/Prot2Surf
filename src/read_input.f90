!> \file read_input.f90
!! \brief Input readers for PDB and association/complexes files.
!!
!! @author Abraham Muñiz-Chicharro
!! @version 1.0
!! @date 2026-04-05
module read_input

    !> Routines to read input files (PDB and association/complex files).
    !!
    !! This module provides small wrappers to read PDB and association
    !! complex files into the project's Fortran derived types. The
    !! routines check file existence, allocate appropriate structures
    !! and fill them using helpers from `mod_pdb` and `mod_assoc`.
    !!
    !! @author Abraham Muñiz-Chicharro

    ! Workflow summary:
    ! - check that input files exist and can be opened
    ! - count data sizes from file contents
    ! - allocate output data structures
    ! - fill structures with parsed values

    USE mod_pdb
    USE mod_assoc

    contains

    !> Read a PDB file into a `type_pdb_file` object.
    !!
    !! Validates the existence of `filename`, opens it, determines the
    !! counts for atoms, residues and chains, allocates `pdb` and fills
    !! it using `allocate_pdb_object` and `fill_pdb_object`.
    !!
    !! @param[out] pdb           Allocated and filled PDB object
    !! @param[in]  filename      Path to PDB file
    !! @param[out] tot_atoms     Number of atoms found
    !! @param[out] tot_residues  Number of residues found
    !! @param[out] tot_chains    Number of chains found
    subroutine read_pdb(pdb, filename, tot_atoms, tot_residues, tot_chains)
        type ( type_pdb_file ) pdb
        integer :: input_pdb, pdb_stat
        logical :: pdb_ex
        character*128, intent(in) :: filename
        integer, intent(out) :: tot_atoms, tot_residues, tot_chains

        ! Validate file existence before attempting open/read.
        inquire(FILE=filename,EXIST=pdb_ex)
        if (.NOT.pdb_ex) then
            write(*,*) "PDB file not found"
            STOP 1
        end if

        ! Open the file and fail fast on I/O errors.
        open (input_pdb,FILE=filename,FORM='FORMATTED',IOSTAT=pdb_stat)
        if (pdb_stat.NE.0) then
            write (*,*) "Error opening PDB file"
            STOP 1
        end if
    
        ! Count atoms/residues/chains to size arrays.
        tot_atoms = count_atoms(input_pdb)
        tot_residues = count_residues(input_pdb)
        tot_chains = count_chains(input_pdb)

        ! Allocate and fill the PDB structure.
        call allocate_pdb_object(pdb,input_pdb,tot_atoms,tot_residues,tot_chains)
        call fill_pdb_object(pdb,input_pdb)

    end subroutine read_pdb

    !> Read an association/complexes file into a `type_assoc_file` object.
    !!
    !! Checks the file exists, opens it and uses `size_assoc`,
    !! `allocate_assoc_object` and `fill_assoc_object` to construct the
    !! `assoc` structure.
    !!
    !! @param[out] assoc          Allocated and filled assoc object
    !! @param[in]  filename       Path to complexes file
    !! @param[out] tot_encounters Number of encounter records read
    subroutine read_assoc(assoc, filename, tot_encounters)
        type ( type_assoc_file ) :: assoc
        integer :: input_assoc, assoc_stat
        logical :: assoc_ex
        character*128, intent(in) :: filename
        integer, intent(out) :: tot_encounters

        ! Validate file existence before opening.
        inquire(FILE=filename,EXIST=assoc_ex)
        if (.NOT.assoc_ex) then
            write(*,*) "Association complexes file not found"
            STOP 1
        end if

        input_assoc = 23
        ! Open file and stop on read/open errors.
        open (input_assoc,FILE=filename,FORM='FORMATTED',IOSTAT=assoc_stat)
        if (assoc_stat.NE.0) then
            write (*,*) "Error opening association file"
            STOP 1
        end if

        ! Count data rows (excluding 4 header lines).
        tot_encounters = max(0, size_assoc(input_assoc) - 4)

        ! Allocate arrays and read header/transform values.
        call allocate_assoc_object(assoc,tot_encounters)
        call fill_assoc_object(assoc,input_assoc)

    end subroutine read_assoc
end module read_input