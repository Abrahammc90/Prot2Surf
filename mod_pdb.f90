!> Module to define a pdb type and pdb reading and writing routine for tools.
!! prepare_atom_protein routine has poor modularity making use in tools difficult.
!!  -- could be merged later
!!  -- need to add into protein object structure
!! @ingroup tools
!! @author Neil Bruce
!! 
module mod_pdb

    !> Data type for storing binned data (eg. Radial distribution functions)
    type type_pdb_file
       !> Number of atoms in pdb file
       integer :: natoms, nresidues, nchains
       !> Flag for record type, 1 = ATOM, 2 = HETATM
       type(type_pdb_atom), dimension (:), allocatable :: atoms
       type(type_pdb_residue), dimension (:), allocatable :: residues
       type(type_pdb_chain), dimension (:), allocatable :: chains
    end type type_pdb_file

    type type_pdb_chain
       !> Number of atoms in pdb file
       integer :: natoms, nresidues
       !> Chain id
       character :: chainid
       !> Flag for record type, 1 = ATOM, 2 = HETATM
       type(type_pdb_atom), dimension (:), allocatable :: atoms
       type(type_pdb_residue), dimension (:), allocatable :: residues
    end type type_pdb_chain

    type type_pdb_residue
       !> Number of atoms in pdb file
       integer :: natoms
       !> Chain id
       character :: chainid
       !> Residue name
       character (len=3) :: resname
       !> Residue number
       integer :: resid
       !> Flag for record type, 1 = ATOM, 2 = HETATM
       type(type_pdb_atom), dimension (:), allocatable :: atoms
    end type type_pdb_residue

    type type_pdb_atom
       !> Number of atoms in pdb file
       !> Flag for record type, 1 = ATOM, 2 = HETATM
       integer :: record
       !> Position, occupancy and temp factor - can be used to store other data eg in pqr
       real (kind=8), dimension (3)  :: coord
       real (kind=8) :: occ,beta
       !> Atom name, res name and chain name
       character (len=4) :: name
       character (len=3) :: resname
       character (len=1) :: chainid 
       !> Atom and residue numbers
       integer :: id, resid 
 
    end type type_pdb_atom
 
 contains
  
    !> Create a new histogram instance
    !! @param this : histogram instance to be created
    !! @param length : number of elements in histogram array
    !! @param unit_width : unit value of one bin width
    !! @param opt_initial_value : unit value of beginning of first bin (zero if not set)
    subroutine allocate_pdb_object(this, pdb_file, natoms, nresidues, nchains)
       IMPLICIT NONE
 
       ! Arguments
       type(type_pdb_file), intent(inout) :: this
       integer, intent(in) :: pdb_file
       integer, intent(in) :: natoms, nresidues, nchains
       integer :: stat
       integer :: tot_chain_atoms, tot_chain_residues, tot_residue_atoms
       integer :: residue_atoms_i, chain_residue_atoms_i
       integer :: residue_id, chain_residue_id, chain_id
       character(len=5) :: resid, prev_resid
       logical :: isopen, restarted_chain
       character(100) :: pdb_line
 
       ! Internal
       integer :: ierr
      
       stat = 0

       this%natoms = natoms
       this%nresidues = nresidues
       this%nchains = nchains
       

       allocate (this%atoms(natoms), STAT=ierr)
       if (ierr.NE.0) then
          write (*,*) "Error allocating PDB object"
          STOP 1
       end if
       allocate (this%residues(nresidues), STAT=ierr)
       if (ierr.NE.0) then
         write (*,*) "Error allocating PDB object"
         STOP 1
       end if
       allocate (this%chains(nchains), STAT=ierr)
       if (ierr.NE.0) then
          write (*,*) "Error allocating PDB object"
          STOP 1
       end if

       ! Check that file has been opened
       inquire(pdb_file,OPENED=isopen)
       if (.NOT.isopen) then
          write(*,*) "PDB file is not open for reading"
          STOP 1
       end if
       
       !Rewind pdb
       rewind(pdb_file)
     
       ! Read first line
       read(pdb_file, '(a)', IOSTAT=stat) pdb_line
       if (STAT.NE.0) then
          write (*,*) "Could not read first line of PDB file. Is file empty"
          STOP 1
       end if
     
       ! Allocate atoms and residues in chains

       tot_chain_residues = 1
       tot_chain_atoms = 0
       chain_id = 1
       
       residue_id = 0

       do while (stat.EQ.0)
       
         ! Check if line contains atomic data
         if ((pdb_line(1:4).EQ."ATOM").OR.(pdb_line(1:6).EQ."HETATM")) then
            
            tot_chain_atoms = tot_chain_atoms + 1
            resid = pdb_line(22:26)

            if ( residue_id == 0 ) then
               prev_resid = resid
               residue_id = residue_id + 1
            end if
            
            if ( resid .NE. prev_resid ) then
               tot_chain_residues = tot_chain_residues + 1
               prev_resid = resid
            end if

         else if (pdb_line(1:3).EQ."TER") then
           allocate(this%chains(chain_id)%atoms(tot_chain_atoms), STAT=ierr)
           if (ierr.NE.0) then
              write (*,*) "Error allocating PDB chain"
              STOP 1
           end if
           allocate(this%chains(chain_id)%residues(tot_chain_residues), STAT=ierr)
           if (ierr.NE.0) then
              write (*,*) "Error allocating PDB chain"
              STOP 1
           end if
           
           tot_chain_atoms = 0
           tot_chain_residues = 0
           chain_id = chain_id + 1
         end if
         ! Read next line
         read(pdb_file, '(a)', IOSTAT=stat) pdb_line
       end do

       rewind(pdb_file)
       stat = 0

       ! Read first atom
       do while (stat.EQ.0)
         read(pdb_file, '(a)', IOSTAT=stat) pdb_line
         if (STAT.NE.0) then
            write (*,*) "Could not read first line of PDB file. Is file empty"
            STOP 1
         end if
         if ((pdb_line(1:4).EQ."ATOM").OR.(pdb_line(1:6).EQ."HETATM")) then
            resid = pdb_line(22:26)
            prev_resid = resid
            EXIT
         end if
       end do

       ! Allocate atoms in residues
       tot_residue_atoms = 0

       residue_id = 0
       chain_residue_id = 0
       chain_id = 1
      
       restarted_chain = .false.
       !prev_resid = ''
       !resid = ''

       residue_atoms_i = 0 
       chain_residue_atoms_i = 0

       do while (stat.EQ.0)

          if ((pdb_line(1:4).EQ."ATOM").OR.(pdb_line(1:6).EQ."HETATM")) then
            resid = pdb_line(22:26)
          end if

          if ( ( resid .NE. prev_resid ) .or. ( restarted_chain ) ) then

            residue_id = residue_id + 1
            chain_residue_id = chain_residue_id + 1

            tot_residue_atoms = residue_atoms_i

            residue_atoms_i = 0
            chain_residue_atoms_i = 0

            !write(*,*) residue_id, tot_residue_atoms
            allocate(this%residues(residue_id)%atoms(tot_residue_atoms), STAT=ierr)
            if (ierr.NE.0) then
               write (*,*) "Error allocating PDB chain"
               STOP 1
            end if

            !write(*,*) chain_id, chain_residue_id, tot_residue_atoms
            !allocate(this%chains(chain_id)%residues(chain_residue_id)%atoms(tot_residue_atoms), STAT=ierr)
            !if (ierr.NE.0) then
            !   !STOP 1
            !   write (*,*) "Error allocating PDB chain"
            !   STOP 2
            !end if

            if ( restarted_chain ) then
               chain_residue_id = 0
               chain_id = chain_id + 1
            end if

            prev_resid = resid
            restarted_chain = .false.
          
          end if

          if ((pdb_line(1:4).EQ."ATOM").OR.(pdb_line(1:6).EQ."HETATM")) then

            residue_atoms_i = residue_atoms_i + 1
            chain_residue_atoms_i = chain_residue_atoms_i + 1

          else if ( pdb_line(1:3).EQ."TER" ) then

            restarted_chain = .true.
            !write(*,*) chain_residue_id, tot_chain_residue_atoms
            
          end if

          read(pdb_file, '(a)', IOSTAT=stat) pdb_line
       end do
       !STOP 1
       
    end subroutine allocate_pdb_object


 
    !> Function that runs through lines of a PDB file and calculates the number of 
    !! ATOM or HETATM records which is returned.
    !! @param pdb_file : record number for PDB file. Should already be opened in calling routine
    !! 
    function count_atoms(pdb_file) result(natoms)
  
    IMPLICIT NONE
 
    integer :: pdb_file
    logical :: isopen
    character(100) :: pdb_line
    integer :: stat
    integer :: natoms
 
    natoms = 0
 
    ! Check that file has been opened
    inquire(pdb_file,OPENED=isopen)
    if (.NOT.isopen) then
       write(*,*) "PDB file is not open for reading"
       STOP 1
    end if
    
    !Rewind pdb
    rewind(pdb_file)
 
    ! Read first line
    read(pdb_file, '(a)', IOSTAT=stat) pdb_line
    if (STAT.NE.0) then
       write (*,*) "Could not read first line of PDB file. Is file empty"
       STOP 1
    end if
 
    do while (STAT.EQ.0)
       ! Check if line contains atomic data
       if ((pdb_line(1:4).EQ."ATOM").OR.(pdb_line(1:6).EQ."HETATM")) then
          natoms = natoms + 1
       end if
       ! Read next line
       read(pdb_file, '(a)', IOSTAT=stat) pdb_line
    end do
    

    end function count_atoms

    !> Function that runs through lines of a PDB file and calculates the number of 
    !! ATOM or HETATM records which is returned.
    !! @param pdb_file : record number for PDB file. Should already be opened in calling routine
    !! 
    function count_residues(pdb_file) result(nresidues)
  
      IMPLICIT NONE
   
      integer :: pdb_file
      logical :: isopen
      character(100) :: pdb_line
      integer :: stat
      integer :: nresidues
      character(len=5) :: resid, prev_resid
   
      nresidues = 0
   
      ! Check that file has been opened
      inquire(pdb_file,OPENED=isopen)
      if (.NOT.isopen) then
         write(*,*) "PDB file is not open for reading"
         STOP 1
      end if
      
      !Rewind pdb
      rewind(pdb_file)
   
      ! Read first line
      read(pdb_file, '(a)', IOSTAT=stat) pdb_line
      if (STAT.NE.0) then
         write (*,*) "Could not read first line of PDB file. Is file empty"
         STOP 1
      end if
      
      do while (STAT.EQ.0)
         ! Check if line contains atomic data
         if ((pdb_line(1:4).EQ."ATOM").OR.(pdb_line(1:6).EQ."HETATM")) then
            resid = pdb_line(22:26)
            
            if (nresidues == 0) then
               prev_resid = resid
               nresidues = nresidues + 1
            end if

            if ( resid .NE. prev_resid ) then
               nresidues = nresidues + 1
               prev_resid = resid
            end if
         end if
         ! Read next line
         read(pdb_file, '(a)', IOSTAT=stat) pdb_line
      end do
  
      end function count_residues

    !> Function that runs through lines of a PDB file and calculates the number of 
    !! ATOM or HETATM records which is returned.
    !! @param pdb_file : record number for PDB file. Should already be opened in calling routine
    !! 
    function count_chains(pdb_file) result(nchains)
  
      IMPLICIT NONE
   
      integer :: pdb_file
      logical :: isopen
      character(100) :: pdb_line
      integer :: stat
      integer :: nchains
   
      nchains = 0
   
      ! Check that file has been opened
      inquire(pdb_file,OPENED=isopen)
      if (.NOT.isopen) then
         write(*,*) "PDB file is not open for reading"
         STOP 1
      end if
      
      !Rewind pdb
      rewind(pdb_file)
   
      ! Read first line
      read(pdb_file, '(a)', IOSTAT=stat) pdb_line
      if (STAT.NE.0) then
         write (*,*) "Could not read first line of PDB file. Is file empty"
         STOP 1
      end if
   
      do while (STAT.EQ.0)
         ! Check if line contains atomic data
         if (pdb_line(1:3).EQ."TER") then
            nchains = nchains + 1
         end if
         ! Read next line
         read(pdb_file, '(a)', IOSTAT=stat) pdb_line
      end do
      end function count_chains
 
    !>
    !!
    subroutine fill_pdb_object (this, pdb_file)
 
    IMPLICIT NONE
 
    ! Arguments
    type(type_pdb_file), intent(inout) :: this
    type(type_pdb_atom) :: atom
    !type(type_pdb_chain) :: chain
    integer, intent(in) :: pdb_file
 
    ! Variables
    integer :: chain_atom_i, chain_residue_i, chain_residue_atom_i, residue_atom_i
    integer :: residue_id, chain_residue_id, chain_id, atom_id
    integer :: prev_resid
    character (len=3) :: prev_resname
    integer :: stat
    logical :: isopen, restarted_chain
    character(100) :: pdb_line

    ! Check that file has been opened
    inquire(pdb_file,OPENED=isopen)
    if (.NOT.isopen) then
       write(*,*) "PDB file is not open for reading"
       STOP 1
    end if
    
    !Rewind pdb
    rewind(pdb_file)
    stat = 0
 
    ! Read first atom
    !do while (stat.EQ.0)
      read(pdb_file, '(a)', IOSTAT=stat) pdb_line
      if (STAT.NE.0) then
         write (*,*) "Could not read first line of PDB file. Is file empty"
         STOP 1
      end if
    !  if ((pdb_line(1:4).EQ."ATOM").OR.(pdb_line(1:6).EQ."HETATM")) then
    !     resid = pdb_line(22:26)
    !     prev_resid = resid
    !     EXIT
    !  end if
    !end do


 
    ! Initialise counter
    chain_atom_i = 0
    chain_residue_i = 1
    !chain_residue_atom_i = 0
    residue_atom_i = 0

    atom_id = 0
    residue_id = 0
    chain_residue_id = 1
    chain_id = 1

    restarted_chain = .false.


    do while (stat.EQ.0)
 
       ! Check if line contains atomic data
       if ((pdb_line(1:4).EQ."ATOM").OR.(pdb_line(1:6).EQ."HETATM")) then

          ! Scan line to fill data
          if (pdb_line(1:4).EQ."ATOM") then
             atom%record = 1
          else if (pdb_line(1:6).EQ."HETATM") then
             atom%record = 2
          end if
    
          read (pdb_line(7:26),'(I5,X,A4,X,A3,X,A,I4)') atom%id, atom%name, atom%resname, &
                                                         atom%chainid, atom%resid 
 
          !write (*,'(I5,"c  c",A4,"c  c",A3,"c  c",A,"c  c",I4)') atom%atomid(i), atom%atomname(i), atom%resname(i), &
          !                                               atom%chainid(i), atom%resid(i) 
         
          read (pdb_line(31:54),'(3F8.3)') atom%coord(1), atom%coord(2), atom%coord(3)

          if ( residue_id == 0 ) then
            prev_resid = atom%resid
            residue_id = residue_id + 1
            prev_resname = atom%resname
          end if
          
          if ( ( atom%resid .NE. prev_resid ) .and. ( .not. restarted_chain ) ) then

            !write(*,*) residue_atom_i, residue_id

            this % residues(residue_id) % natoms = residue_atom_i
            this % residues(residue_id) % resid = prev_resid
            this % residues(residue_id) % resname = prev_resname

            this % chains(chain_id) % residues(chain_residue_id) = this % residues(residue_id)

            residue_atom_i = 0
            !chain_residue_atom_i = 0
            
            chain_residue_i = chain_residue_i + 1
            chain_residue_id = chain_residue_id + 1
            residue_id = residue_id + 1

            prev_resid = atom%resid
            prev_resname = atom%resname

            !this % chains(chain_id) % residues(residue_id) % resid = atom%resid
            !this % chains(chain_id) % residues(residue_id) % resname = atom%resname

          end if

          !restarted_chain = .false.
          
          atom_id = atom_id + 1
          chain_atom_i = chain_atom_i + 1
          residue_atom_i = residue_atom_i + 1
          !chain_residue_atom_i = chain_residue_atom_i + 1
   
          !write(*,*) residue_atom_i
          if ( residue_id .gt. size(this % residues) ) STOP 1
          if ( chain_id .gt. size(this % chains) ) STOP 2
          if ( residue_atom_i .gt. size(this % residues(residue_id) % atoms) ) STOP 3
          if ( chain_residue_i .gt. size(this % chains(chain_id) % residues) ) STOP 4
          !if ( chain_residue_atom_i .gt. size(this % chains(chain_id) % residues(chain_residue_id) % atoms) ) STOP 5
          if ( chain_atom_i .gt. size( this % chains ( chain_id ) % atoms ) ) STOP 6


          this % atoms(atom_id) = atom
          this % residues(residue_id) % atoms(residue_atom_i) = atom
          this % chains(chain_id) % atoms(chain_atom_i) = atom
          !this % chains(chain_id) % residues(chain_residue_id) % atoms(chain_residue_atom_i) = atom

          
       else if (pdb_line(1:3).EQ."TER") then

         !write(*,*) residue_atom_i, residue_id
         this % residues(residue_id) % natoms = residue_atom_i
         this % residues(residue_id) % resid = prev_resid
         this % residues(residue_id) % resname = prev_resname
         this % chains(chain_id) % residues(chain_residue_id) = this % residues(residue_id)
         !residue_atom_i = 0

         this % chains(chain_id) % chainid = atom % chainid
         this % chains(chain_id) % nresidues = chain_residue_i
         this % chains(chain_id) % natoms = chain_atom_i

         !write(*,*) size( this%chains(chain_id)%residues(chain_residue_id)%atoms )

         chain_atom_i = 0
         chain_residue_i = 0
         chain_residue_id = 0
         chain_id = chain_id + 1

         !restarted_chain = .True.

       end if
 
       
       ! Read next line
       read(pdb_file, '(a)', IOSTAT=stat) pdb_line

    end do
    !STOP 1
 
    end subroutine fill_pdb_object
  
 end module mod_pdb