MODULE mod_memory

  !> Utilities for memory management and batch processing.
  !!
  !! This module provides routines to detect available system memory,
  !! calculate batch sizes for large datasets, and manage memory-efficient
  !! processing of large matrices.
  !!
  !! @author Abraham Muñiz-Chicharro
  !!

  implicit none
  private
  public :: get_available_memory, calculate_max_batch_size

contains

  !> Get available system memory in bytes
  !!
  !! Detects available RAM on Linux systems by parsing /proc/meminfo.
  !! On other systems, returns -1.
  !!
  !! @return available memory in bytes, or -1 if cannot be determined
  function get_available_memory() result(memsize)
    implicit none
    integer(kind=8) :: memsize
    integer :: unit, ios
    character(len=256) :: line
    character(len=32) :: mem_str
    integer :: i, idx

    memsize = -1_8

    ! Try /proc/meminfo (Linux)
    open(newunit=unit, file='/proc/meminfo', status='old', action='read', &
         iostat=ios)
    
    if (ios == 0) then
      do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        
        ! Look for MemAvailable line
        if (index(line, 'MemAvailable') == 1) then
          ! Extract the numeric value
          idx = index(line, ':') + 1
          mem_str = adjustl(line(idx:))
          ! Remove 'kB' or other units
          do i = len_trim(mem_str), 1, -1
            if (mem_str(i:i) >= '0' .and. mem_str(i:i) <= '9') then
              mem_str = mem_str(1:i)
              exit
            end if
          end do
          read(mem_str, *, iostat=ios) memsize
          if (ios == 0) then
            memsize = memsize * 1024_8  ! Convert from kB to bytes
          end if
          exit
        end if
      end do
      close(unit)
    end if

  end function get_available_memory


  !> Calculate maximum batch size for matrix processing
  !!
  !! Given available memory and memory to leave free (in percentage),
  !! calculates the maximum number of encounters that can be processed.
  !!
  !! Formula:
  !!   usable_memory = available_memory * (1 - reserve_percent/100)
  !!   bytes_per_encounter = nb_total_encounters * 8 + other_overhead
  !!   max_batch = usable_memory / bytes_per_encounter
  !!
  !! @param[in]  available_memory    Available system RAM in bytes
  !! @param[in]  nb_total_encounters Total number of encounters to process
  !! @param[in]  reserve_percent     Percentage of memory to keep free (default 20)
  !!
  !! @return maximum batch size (number of encounters per batch)
  function calculate_max_batch_size(available_memory, nb_total_encounters, reserve_percent) &
           result(max_batch)
    implicit none
    integer(kind=8), intent(in) :: available_memory
    integer, intent(in) :: nb_total_encounters
    integer, intent(in), optional :: reserve_percent
    integer :: max_batch

    integer :: reserve
    real(kind=8) :: usable_memory, bytes_per_row, approx_batch
    real(kind=8) :: overhead_per_matrix

    ! Default reserve is 20%
    if (present(reserve_percent)) then
      reserve = reserve_percent
    else
      reserve = 20
    end if

    ! Calculate usable memory (leaving reserve% free)
    usable_memory = real(available_memory, kind=8) * (1.0d0 - real(reserve, kind=8) / 100.0d0)

    ! Each row of the full matrix requires:
    ! - 8 bytes per encounter (real(kind=8))
    ! - One full row for all encounters
    ! Each distance calculation also has some overhead
    
    ! For a safe estimate:
    ! bytes_per_row = nb_total_encounters * 8 (the full matrix row)
    ! But we also need space for temporary arrays
    ! Add 30% overhead for temp arrays and other data structures
    
    bytes_per_row = real(nb_total_encounters, kind=8) * 8.0d0
    overhead_per_matrix = bytes_per_row * 0.3d0  ! 30% overhead
    
    ! Calculate approximate batch size
    ! If we process batch_size encounters, we need:
    ! batch_size * nb_total_encounters * 8 bytes for the matrix rows
    ! Plus overhead for coordinates and other arrays
    
    approx_batch = usable_memory / (bytes_per_row + overhead_per_matrix)
    
    ! Ensure at least 1 batch
    max_batch = max(1, int(approx_batch))
    
    ! Cap at total encounters
    max_batch = min(max_batch, nb_total_encounters)

  end function calculate_max_batch_size

end MODULE mod_memory
