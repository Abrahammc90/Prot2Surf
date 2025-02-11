module read_input

    USE mod_pdb
    USE mod_assoc

    contains

    subroutine read_pdb(pdb, filename, tot_atoms, tot_residues, tot_chains)
        type ( type_pdb_file ) pdb
        integer :: input_pdb, pdb_stat
        logical :: pdb_ex
        character*128, intent(in) :: filename
        integer, intent(out) :: tot_atoms, tot_residues, tot_chains

        inquire(FILE=filename,EXIST=pdb_ex)
        if (.NOT.pdb_ex) then
            write(*,*) "PDB file not found"
            STOP 1
        end if

        open (input_pdb,FILE=filename,FORM='FORMATTED',IOSTAT=pdb_stat)
        if (pdb_stat.NE.0) then
            write (*,*) "Error opening PDB file"
            STOP 1
        end if
    
        tot_atoms = count_atoms(input_pdb)
        tot_residues = count_residues(input_pdb)
        tot_chains = count_chains(input_pdb)
        call allocate_pdb_object(pdb,input_pdb,tot_atoms,tot_residues,tot_chains)
        call fill_pdb_object(pdb,input_pdb)

    end subroutine read_pdb

    subroutine read_assoc(assoc, filename, tot_encounters)
        type ( type_assoc_file ) :: assoc
        integer :: input_assoc, assoc_stat
        logical :: assoc_ex
        character*128, intent(in) :: filename
        integer, intent(out) :: tot_encounters

        inquire(FILE=filename,EXIST=assoc_ex)
        if (.NOT.assoc_ex) then
            write(*,*) "Association complexes file not found"
            STOP 1
        end if

        open (input_assoc,FILE=filename,FORM='FORMATTED',IOSTAT=assoc_stat)
        if (assoc_stat.NE.0) then
            write (*,*) "Error opening association file"
            STOP 1
        end if

        tot_encounters = size_assoc(input_assoc)
        call allocate_assoc_object(assoc,tot_encounters)
        call fill_assoc_object(assoc,input_assoc)

    end subroutine read_assoc
end module read_input