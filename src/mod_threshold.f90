MODULE mod_threshold

  !> \file mod_threshold.f90
  !! \brief Utilities to compute threshold arrays from encounter complexes.
  !!
  !! @author Abraham Muñiz-Chicharro
  !! @version 1.0
  !! @date 2026-04-05
  USE maths
  USE OMP_LIB

  ! Workflow summary:
  ! - compute one value per encounter for threshold filtering
  ! - support z-coordinate, minimum-distance, and angle values
  ! - provide sorting helpers for values and indexes

    contains

    !> Compute the Z-coordinate (along transformed z-axis) for each encounter.
    !!
    !! Recenters each encounter at `xc2`, applies rotation/translation using
    !! `rot1/rot2` and `trans_vector`, and extracts the Z coordinate of the
    !! first atom in `solute_crds` for every encounter into `array`.
    !!
    !! @param[out] array        Output array of length `n` with z-coordinates
    !! @param[in]  n            Number of encounters
    !! @param[in]  nb_atoms     Number of atoms in `solute_crds`
    !! @param[in]  xc1, xc2     Centers used for transform
    !! @param[in]  trans_vector Translation vectors for each encounter (n, 3)
    !! @param[in]  rot1, rot2   Rotation vectors for each encounter (n, 3)
    !! @param[in]  solute_crds  Coordinates of solute atoms (nb_atoms, 3)
    subroutine array_z_coord(array, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, solute_crds)

      IMPLICIT NONE
      real(kind=8), dimension(n), intent(out) :: array
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(:, :), intent(in) :: solute_crds

      real(kind=8), dimension(nb_atoms, 3) :: new_coord
      integer :: i, j, progress_index
      real(kind=8) :: value_i, value_j

      progress_index = 0

      ! Parallelize the outer loop
      !$OMP PARALLEL DO PRIVATE(j, new_coord, value_i, value_j) SCHEDULE(DYNAMIC)
      do i = 1, n
        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, solute_crds, new_coord)
        value_i = new_coord(1, 3)
        array(i) = value_i

        ! Ensure orderly printing using an atomic increment
        !$OMP CRITICAL
        progress_index = progress_index + 1
        ! Print the progress every 100 iterations
        if (mod(progress_index, 100) == 0) then
          write(*,*) 'Encounters processed: ', progress_index
        end if
        !$OMP END CRITICAL

      end do
      !$OMP END PARALLEL DO

    end subroutine array_z_coord

    !> Compute minimum inter-atomic distance between two solutes for each encounter.
    !!
    !! For each encounter, the second solute is transformed and the minimum
    !! distance between any atom of `solute1_crds` and the transformed
    !! `solute2_crds` is computed and stored in `array`.
    !!
    !! @param[out] array         Output array of length `n` with minimum distances
    !! @param[in]  n             Number of encounters
    !! @param[in]  nb_atoms      Number of atoms in `solute2_crds`
    !! @param[in]  xc1, xc2      Centers used for transform
    !! @param[in]  trans_vector  Translation vectors for each encounter (n, 3)
    !! @param[in]  rot1, rot2    Rotation vectors for each encounter (n, 3)
    !! @param[in]  solute1_crds  Coordinates of solute1 atoms
    !! @param[in]  solute2_crds  Coordinates of solute2 atoms
    subroutine array_atoms_dist(array, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, solute1_crds, solute2_crds)

      IMPLICIT NONE
      real(kind=8), dimension(n), intent(out) :: array
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(:, :), intent(in) :: solute1_crds, solute2_crds

      real(kind=8), dimension(:), allocatable :: distances
      real(kind=8), dimension(nb_atoms, 3) :: new_coord
      integer :: i, j, k, tot_coords1, progress_index
      real(kind=8) :: min_i, min_j, dist

      tot_coords1 = size(solute1_crds(:, 1))
      allocate(distances(tot_coords1))
      progress_index = 0

      ! Compute minimum-distance per encounter.
      !$OMP PARALLEL DO PRIVATE(j, k, new_coord, distances, min_i, min_j, dist) SCHEDULE(DYNAMIC)
      do i = 1, n
        distances = 999999.9

        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, solute2_crds, new_coord)

        do k = 1, tot_coords1
          call calculate_distance(new_coord, solute1_crds(k, :), dist)
          distances(k) = dist
        end do

        min_i = minval(distances)
        array(i) = min_i

        ! Ensure orderly printing using an atomic increment
        !$OMP CRITICAL
        progress_index = progress_index + 1
        ! Print the progress every 100 iterations
        if (mod(progress_index, 100) == 0) then
          write(*,*) 'Encounters processed: ', progress_index
        end if
        !$OMP END CRITICAL

      end do
      !$OMP END PARALLEL DO

      deallocate(distances)

    end subroutine array_atoms_dist

  
    !Subroutine array_angle
! 
    !> Compute angle between two vectors (defined by atom groups) per encounter.
    !!
    !! Transforms the pair of points from `solute2` for each encounter and
    !! computes the angle between the vector defined by `point1b-point1a`
    !! and the vector defined by the transformed `point2b-point2a`.
    !!
    !! @param[out] array         Output array of length `n` with angles in degrees
    !! @param[in]  n             Number of encounters
    !! @param[in]  nb_atoms      Number of atoms used for transformation
    !! @param[in]  xc1, xc2      Centers used for transform
    !! @param[in]  trans_vector  Translation vectors for each encounter (n, 3)
    !! @param[in]  rot1, rot2    Rotation vectors for each encounter (n, 3)
    !! @param[in]  point1a       First point defining reference vector
    !! @param[in]  point1b       Second point defining reference vector
    !! @param[in]  point2a       First point of vector to be transformed
    !! @param[in]  point2b       Second point of vector to be transformed
    !! @param[in]  dimensions    Use 2 or 3 to select projection
    subroutine array_angle(array, n, nb_atoms, &
      xc1, xc2, trans_vector, rot1, rot2, &
      point1a, point1b, point2a, point2b, dimensions)

      IMPLICIT NONE
      real(kind=8), dimension(n), intent(out) :: array
      integer, intent(in) :: n, nb_atoms
      real(kind=8), dimension(3), intent(in) :: xc1, xc2
      real(kind=8), dimension(:, :), intent(in) :: trans_vector, rot1, rot2
      real(kind=8), dimension(:), intent(in) :: point1a, point1b, point2a, point2b
      integer, intent(in) :: dimensions

      real(kind=8), dimension(nb_atoms, 3) :: new_coord_1, new_coord_2
      real(kind=8), dimension(3) :: v1, v2
      real(kind=8) :: theta
      integer :: i, j, progress_index
      real(kind=8), dimension(2, 3) :: solute2_points

      progress_index = 0

      
      !print *, point1
      !print *, point2
      !print *, point3
      !print *, point4
      !STOP

      ! Initialize the reference vector
      v1(1) = point1b(1) - point1a(1)
      v1(2) = point1b(2) - point1a(2)
      v1(3) = point1b(3) - point1a(3)

      solute2_points(1, :) = point2a(:)
      solute2_points(2, :) = point2b(:)

      ! Transform solute2 vector per encounter and compute
      ! the requested angular dimension (2D or 3D).
      ! Parallelizing outer loop with OpenMP
      !$OMP PARALLEL DO PRIVATE(j, new_coord_1, new_coord_2, v2, theta) SCHEDULE(DYNAMIC)
      do i = 1, n
        call update_complex(xc1, xc2, trans_vector(i, :), &
        rot1(i, :), rot2(i, :), nb_atoms, solute2_points, new_coord_1)
        
        !print *, point1
        !print *, point2
        !print *, new_coord_1(1, :)
        !print *, new_coord_1(2, :)
        !STOP

        v2(1) = new_coord_1(1, 1) - new_coord_1(2, 1)
        v2(2) = new_coord_1(1, 2) - new_coord_1(2, 2)
        v2(3) = new_coord_1(1, 3) - new_coord_1(2, 3)

        if ( dimensions .eq. 2 ) then
          call vectors_angle_2D(v1, v2, theta)
        else if (dimensions .eq. 3) then
          call vectors_angle_3D(v1, v2, theta)
        else
          print *, "Dimensions indicated is ", dimensions
          print *, "Dimensions implemented are 2 or 3."
          STOP 1
        end if
        
        !print *, 'vector 1', v1
        !print *, 'vector 2', v2
        !print *, 'angle', theta1
        !STOP

        array(i) = theta

        ! Ensure orderly printing using an atomic increment
        !$OMP CRITICAL
        progress_index = progress_index + 1
        ! Print the progress every 100 iterations
        if (mod(progress_index, 100) == 0) then
          write(*,*) 'Encounters processed: ', progress_index
        end if
        !$OMP END CRITICAL

      end do
      !$OMP END PARALLEL DO

    end subroutine array_angle

    !> Sort an array in ascending order and reorder encounter indexes accordingly.
    !!
    !! Simple selection sort used to reorder `arr` and the parallel
    !! `encounter_indexes` array. Intended for small to moderate sizes.
    !!
    !! @param[in,out] arr               Array to sort (modified in-place)
    !! @param[in,out] encounter_indexes Parallel index array reordered to match arr
    subroutine sort_array(arr, encounter_indexes)
      implicit none
      real (kind=8), intent(inout), dimension(:) :: arr
      integer, intent(inout), dimension(:) :: encounter_indexes
      integer :: i, j, min_index, tmp_index
      integer :: n
      real(kind=8) :: temp
      
      n = size(arr(:))

      do i = 1, n - 1
          min_index = i
          do j = i + 1, n
              if (arr(j) < arr(min_index)) then
                  min_index = j
              end if
          end do
          if (min_index /= i) then
              temp = arr(i)
              arr(i) = arr(min_index)
              arr(min_index) = temp

                tmp_index = encounter_indexes(i)
              encounter_indexes(i) = encounter_indexes(min_index)
                encounter_indexes(min_index) = tmp_index
          end if
      end do
    end subroutine sort_array
    
END MODULE mod_threshold
