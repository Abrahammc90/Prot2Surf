!> Association/complexes file utilities.
!!
!! This module defines the `type_assoc_file` derived type and routines
!! to allocate and populate it from an association/complexes text
!! file. The `type_assoc_file` stores the raw file lines, header lines,
!! per-encounter translation vectors and rotation vectors, and the two
!! reference centers `xc1` and `xc2` extracted from the file header.
!!
!! NOTE: The routines in this file are adapted from the SDA7 project's
!! `mod_pdb.f90` (original author Neil Bruce). This file contains code
!! and ideas derived from that implementation; appropriate attribution
!! to the SDA7 sources is retained here. Consult the SDA7 project for
!! original authorship and licensing details.
!!
!! @author Abraham Muñiz-Chicharro
module mod_assoc

    !> Data type for storing binned data (eg. Radial distribution functions)
    type type_assoc_file
       !> Number of lines in assoc_file file
       integer :: nlines
       character(len=217), dimension(:), allocatable :: lines
       character(len=217), dimension(4) :: head
       !> Translational vector
       real (kind=8), dimension (:, :), allocatable :: trans_vector
       !> rotX, rotY
       real (kind=8), dimension (:,:), allocatable  :: rot1, rot2
       !> center of mass of p2
       real (kind=8), dimension (3) :: xc1, xc2
 
    end type type_assoc_file
 
 contains
  
   !> Allocate storage for an association/complexes object.
   !!
   !! Allocates the arrays inside `this` to hold `nlines` encounter
   !! records and initializes numeric arrays to zero.
   !!
   !! @param[in,out] this   Instance of `type_assoc_file` to initialize
   !! @param[in]     nlines Number of encounter lines to allocate
   subroutine allocate_assoc_object(this, nlines)
       IMPLICIT NONE
 
       ! Arguments
       type(type_assoc_file), intent(inout) :: this
       integer, intent(in) ::  nlines
 
       ! Internal
       integer :: ierr
 
       this%nlines = nlines

       allocate(this%lines(nlines))
       
       allocate (this%trans_vector(nlines,3), STAT=ierr)
       if (ierr.NE.0) then
          write (*,*) "Error allocating complexes object"
          STOP 1
       else
          this%trans_vector = 0.0
       end if
       allocate (this%rot1(nlines,3), STAT=ierr)
       if (ierr.NE.0) then
          write (*,*) "Error allocating complexes object"
          STOP 1
       else
          this%rot1 = 0.0
       end if
       allocate (this%rot2(nlines,3), STAT=ierr)
       if (ierr.NE.0) then
          write (*,*) "Error allocating complexes object"
          STOP 1
       else
          this%rot2 = 0.0
       end if
 
 
    end subroutine allocate_assoc_object
 
   !> Count encounter records in an association/complexes file.
   !!
   !! Scans `assoc_file` (a Fortran unit number) and counts the number
   !! of non-comment lines. Lines starting with `#` are treated as
   !! comments and ignored. Returns the number of encounter records.
   !!
   !! @param[in] assoc_file Fortran unit number of an already-open file
   !! @return nlines       Number of encounter records found
   function size_assoc(assoc_file) result(nlines)
  
    IMPLICIT NONE
 
    integer :: assoc_file
    logical :: isopen
    character(225) :: assoc_line
    integer :: stat
    integer :: nlines
 
    nlines = 0
 
    ! Check that file has been opened
    inquire(assoc_file,OPENED=isopen)
    if (.NOT.isopen) then
       write(*,*) "Complexes file is not open for reading"
       STOP 1
    end if
    
    !Rewind assoc
    rewind(assoc_file)
 
    ! Read first line
    read(assoc_file, '(a)', IOSTAT=stat) assoc_line
    if (STAT.NE.0) then
       write (*,*) "Could not read first line of complexes file. Is file empty"
       STOP 1
    end if
 


    do while (STAT.EQ.0)
       ! Check if line contains atomic data
       if ( assoc_line(1:1) .NE. "#" ) nlines = nlines + 1
       ! Read next line
       read(assoc_file, '(a)', IOSTAT=stat) assoc_line
    end do

    end function size_assoc
 
   !> Read an association/complexes file and populate `type_assoc_file`.
   !!
   !! Reads header lines from the open unit `assoc_file` and extracts
   !! the two reference centers (`xc1`, `xc2`) from header lines 3
   !! and 4. Subsequent non-header lines are stored into
   !! `this%lines` and the translation and rotation vectors are parsed
   !! into `this%trans_vector`, `this%rot1` and `this%rot2`.
   !!
   !! @param[in,out] this       Allocated `type_assoc_file` to populate
   !! @param[in]     assoc_file Fortran unit number of an already-open file
   subroutine fill_assoc_object (this, assoc_file)
 
    IMPLICIT NONE
 
    ! Arguments
    type(type_assoc_file), intent(inout) :: this
    integer, intent(in) :: assoc_file
 
    ! Variables
    integer :: i, assoc_i
    integer :: stat
    logical :: isopen
    character(len=217) :: assoc_line
 
    ! Check that file has been opened
    inquire(assoc_file,OPENED=isopen)
    if (.NOT.isopen) then
       write(*,*) "Complexes file is not open for reading"
       STOP 1
    end if
    
    !Rewind assoc
    rewind(assoc_file)
 
    ! Read first line
    read(assoc_file, '(a)', IOSTAT=stat) assoc_line
    if (STAT.NE.0) then
       write (*,*) "Could not read first line of complexes file. Is file empty"
       STOP 1
    end if
 
    ! Initialise counter
    i = 0
    assoc_i = 0
    do while (stat.EQ.0)
 
      i = i + 1
       ! Check if line contains atomic data
       if ( i .le. 4 ) then
         this % head(i) = assoc_line
       end if
       if ( i == 3 ) then
         read (assoc_line(2:26), '(3F8.3)') this%xc1(1), this%xc1(2), this%xc1(3)
       else if ( i == 4 ) then
          read (assoc_line(2:26), '(3F8.3)') this%xc2(1), this%xc2(2), this%xc2(3)
          !write(*,*) 'xc2', this%xc2
       else if ( i .gt. 4 ) then
          assoc_i = assoc_i + 1
          this % lines(assoc_i) = assoc_line
          read (assoc_line(17:98),'(9F9.3)') this%trans_vector(assoc_i,1), &
                                              this%trans_vector(assoc_i,2), &
                                              this%trans_vector(assoc_i,3), &
                                              this%rot1(assoc_i,1), &
                                              this%rot1(assoc_i,2), &
                                              this%rot1(assoc_i,3), &
                                              this%rot2(assoc_i,1), &
                                              this%rot2(assoc_i,2), &
                                              this%rot2(assoc_i,3)
       end if
       
       ! Read next line
       read(assoc_file, '(a)', IOSTAT=stat) assoc_line
    end do

    end subroutine fill_assoc_object

  
 end module mod_assoc
