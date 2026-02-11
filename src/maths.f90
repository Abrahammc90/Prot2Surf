MODULE maths

    USE mod_pdb

    contains

    !> Compute the cross product of two 3D vectors.
    !!
    !! Computes v3 = v1 x v2 using the standard right-hand rule.
    !!
    !! @param[in] v1  Left-hand operand vector (3 elements)
    !! @param[in] v2  Right-hand operand vector (3 elements)
    !! @param[out] v3 Resulting cross product (3 elements)
    subroutine cross(v1,v2,v3)
        IMPLICIT NONE
        
        real ( kind=8 ), dimension ( 3 ), intent ( in ) :: v1, v2
        real ( kind=8 ), dimension ( 3 ), intent ( out ) :: v3  
        v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
        v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
        v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
    end subroutine cross

    !> Apply translation and rotation to a set of atomic coordinates.
    !!
    !! Transforms the coordinates of a complex from a frame centered at
    !! `xc2` to a new frame centered at `xc1`. The routine first recenters
    !! coordinates by `xc2`, applies a rotation defined by `rot_vx` and
    !! `rot_vy` (their cross product defines the third axis), then applies
    !! an additional translation `trans`, and writes the results to
    !! `new_coords`.
    !!
    !! @param[in]  xc1        Target center position (3)
    !! @param[in]  xc2        Original center position (3)
    !! @param[in]  trans      Translation to apply after rotation (3)
    !! @param[in]  rot_vx     First rotation axis vector (3)
    !! @param[in]  rot_vy     Second rotation axis vector (3)
    !! @param[in]  nb_atoms   Number of atoms (integer)
    !! @param[in]  coords     Input coordinates, shape (nb_atoms,3)
    !! @param[out] new_coords Output coordinates, shape (nb_atoms,3)
    subroutine update_complex(xc1, xc2, trans, rot_vx, rot_vy, nb_atoms, coords, new_coords)
        IMPLICIT NONE

        real ( kind = 8 ), dimension (:, :), intent(in) :: coords
        real ( kind = 8 ), dimension (:, :), intent(out) :: new_coords
        real ( kind = 8 ), dimension ( 3 ), intent(in) :: trans, rot_vx, rot_vy
        real ( kind = 8 ), dimension ( 3 ), intent(in) :: xc1, xc2
        real ( kind = 8 ), dimension ( 3 ) :: xee, yee, zee
        real ( kind = 8 ), dimension ( 3 ) :: rot_vz, ozz, zoz, zzo
        real ( kind = 8 ), dimension ( 3 ) :: rot_vx_norm, rot_vy_norm, rot_vz_norm
        real ( kind = 8 ) :: magnitude_rot_vx, magnitude_rot_vy, magnitude_rot_vz
        integer :: nb_atoms, i
        real ( kind = 8 ), dimension ( nb_atoms, 3 ) :: trans_coords

        !write(*,*) 'trans', trans
        !write(*,*) 'rot X', rot_vx
        !write(*,*) 'rot Y', rot_vy
        !write(*,*) 'coords', coords(1, :)
        !stop

        !Normalizing vectors
        magnitude_rot_vx = sqrt ( dot_product ( rot_vx, rot_vx ) )
        rot_vx_norm = rot_vx / magnitude_rot_vx
        magnitude_rot_vy = sqrt ( dot_product ( rot_vy, rot_vy ) )
        rot_vy_norm = rot_vy / magnitude_rot_vy
        
        !Generate the z_ax rotation vector
        call cross(rot_vx_norm, rot_vy_norm, rot_vz) 
        magnitude_rot_vz = sqrt ( dot_product ( rot_vz, rot_vz ) )
        rot_vz_norm = rot_vz / magnitude_rot_vz

        !Creates rotation matrix
        xee = (/rot_vx_norm(1), rot_vy_norm(1), rot_vz_norm(1)/)
        yee = (/rot_vx_norm(2), rot_vy_norm(2), rot_vz_norm(2)/)
        zee = (/rot_vx_norm(3), rot_vy_norm(3), rot_vz_norm(3)/)

        !write(*,*) xee
        !write(*,*) yee
        !write(*,*) zee
        !stop

        !Recenter coordinates to origin
        do i = 1, nb_atoms
            new_coords(i, 1) = coords(i, 1) - xc2(1)
            new_coords(i, 2) = coords(i, 2) - xc2(2)
            new_coords(i, 3) = coords(i, 3) - xc2(3)
        end do

        !Apply rotation
        trans_coords = new_coords
        do i = 1, nb_atoms
            new_coords(i, 1) = dot_product(trans_coords(i,:), xee)
            new_coords(i, 2) = dot_product(trans_coords(i,:), yee)
            new_coords(i, 3) = dot_product(trans_coords(i,:), zee)
        end do

        !Recenter to xc1
        do i = 1, nb_atoms
            new_coords(i, 1) = new_coords(i, 1) + xc1(1) + trans(1)
            new_coords(i, 2) = new_coords(i, 2) + xc1(2) + trans(2)
            new_coords(i, 3) = new_coords(i, 3) + xc1(3) + trans(3)
        end do

    end subroutine update_complex

    !> Compute the angle in degrees between two 3D vectors.
    !!
    !! Calculates the angle (in degrees) between `vector1` and `vector2`
    !! using the dot product formula and `acos`. No modification of inputs.
    !!
    !! @param[in]  vector1       First 3D vector (3 elements)
    !! @param[in]  vector2       Second 3D vector (3 elements)
    !! @param[out] theta_degrees Angle between vectors in degrees
    subroutine vectors_angle_3D(vector1, vector2, theta_degrees)

        IMPLICIT NONE
        real (kind=8), intent(out) :: theta_degrees
        real( kind=8 ),dimension( 3 ), intent ( in ) :: vector1, vector2
        real (kind=8) :: magnitude_v1, magnitude_v2, magnitude_v12, cos_theta
        real (kind=8) :: theta_radians
        real (kind=8) :: dot_11, dot_22, dot_12
        integer :: i
        
        dot_11 = 0.0d0
        do i = 1, 3
            dot_11 = dot_11 + vector1(i)*vector1(i)
        end do
        magnitude_v1 = sqrt(dot_11)
        
        dot_22 = 0.0d0
        do i = 1, 3
            dot_22 = dot_22 + vector2(i)*vector2(i)
        end do
        magnitude_v2 = sqrt(dot_22)
        !print *, magnitude_v2
        
        dot_12 = 0.0d0
        do i = 1, 3
            dot_12 = dot_12 + vector1(i)*vector2(i)
        end do

        cos_theta = dot_12 / ( magnitude_v1 * magnitude_v2 )
        
        theta_radians = acos(cos_theta)
        theta_degrees = theta_radians * (180.0 / 3.141592653589793)

    end subroutine vectors_angle_3D

    !> Compute the angle in degrees between two vectors projected to 2D.
    !!
    !! Uses only the first two components of the input vectors to compute
    !! the angle in the XY plane. Returns 0 if either vector has zero
    !! length in the XY plane.
    !!
    !! @param[in]  vector1       First vector (use components 1 and 2)
    !! @param[in]  vector2       Second vector (use components 1 and 2)
    !! @param[out] theta_degrees Angle between projected vectors in degrees
    subroutine vectors_angle_2D(vector1, vector2, theta_degrees)

        IMPLICIT NONE
        real (kind=8), intent(out) :: theta_degrees
        real( kind=8 ),dimension( 3 ), intent ( in ) :: vector1, vector2
        real (kind=8) :: magnitude_v1, magnitude_v2, cos_theta
        real (kind=8) :: theta_radians
        real (kind=8) :: dot_12
        
        ! Compute the dot product of the vectors
        dot_12 = vector1(1)*vector2(1) + vector1(2)*vector2(2)

        ! Compute the norms of the vectors
        magnitude_v1 = sqrt(vector1(1)**2 + vector1(2)**2)
        magnitude_v2 = sqrt(vector2(1)**2 + vector2(2)**2)

        ! Guard against division by zero
        if (magnitude_v1 == 0.0d0 .or. magnitude_v2 == 0.0d0) then
            theta_radians = 0.0d0
            return
        end if

        ! Compute the cosine of the angle
        cos_theta = dot_12 / (magnitude_v1 * magnitude_v2)

        ! Clamp value to [-1, 1] to avoid numerical errors with acos
        cos_theta = max(-1.0d0, min(1.0d0, cos_theta))

        ! Compute the angle in radians
        theta_radians = acos(cos_theta)
        theta_degrees = theta_radians * (180.0 / 3.141592653589793)

        end subroutine vectors_angle_2D

    !> Compute the root-mean-square deviation (RMSD) between two sets of coordinates.
    !!
    !! Calculates RMSD = sqrt( (1/N) * sum_i ||coords1(i,:) - coords2(i,:)||^2 ).
    !!
    !! @param[out] rmsd_value Resulting RMSD (scalar)
    !! @param[in]  nb_atoms   Number of atoms / points (integer)
    !! @param[in]  coords1    First coordinate set, shape (nb_atoms,3)
    !! @param[in]  coords2    Second coordinate set, shape (nb_atoms,3)
    subroutine rmsd ( rmsd_value, nb_atoms, coords1, coords2)
        IMPLICIT NONE

        real(kind=8), intent(out) :: rmsd_value
        integer, intent(in) :: nb_atoms
        real(kind=8), dimension(:, :), intent(in) :: coords1, coords2
        real(kind=8) :: sum_sq_diff
        integer :: i, j

        !print *, coords1(1, :)
        !print *, coords2(1, :)
        !print*, nb_atoms
        !stop

        sum_sq_diff = 0
        do i = 1, nb_atoms
            sum_sq_diff = sum_sq_diff + (coords1(i, 1) - coords2(i, 1))**2 + &
                                        (coords1(i, 2) - coords2(i, 2))**2 + &
                                        (coords1(i, 3) - coords2(i, 3))**2
        end do
        rmsd_value = sqrt(sum_sq_diff/nb_atoms)


    end subroutine rmsd

    !> Compute the center of geometry (centroid) of a set of points.
    !!
    !! Sums the coordinates and divides by the number of points to produce
    !! the centroid.
    !!
    !! @param[out] cog      Centroid vector (3 elements)
    !! @param[in]  coords   Input coordinates, shape (nb_atoms,3)
    !! @param[in]  nb_atoms Number of atoms / points
    subroutine calculate_cog(cog, coords, nb_atoms)

        integer, intent(in) :: nb_atoms
        real (kind=8), dimension(3), intent(out) :: cog
        real (kind=8), dimension(nb_atoms, 3), intent(in) :: coords

        integer :: i

        cog = 0

        do i = 1, nb_atoms
            cog(1) = cog(1) + coords(i, 1)
            cog(2) = cog(2) + coords(i, 2)
            cog(3) = cog(3) + coords(i, 3)
        end do

        cog = cog / nb_atoms

    end subroutine calculate_cog

    !> Compute Euclidean distance between two 3D points.
    !!
    !! @param[in]  p1   First point (3 elements)
    !! @param[in]  p2   Second point (3 elements)
    !! @param[out] dist Euclidean distance ||p1 - p2||
    subroutine calculate_distance(p1, p2, dist)

        real (kind=8), dimension(3), intent(in) :: p1, p2
        real (kind=8), intent(out) :: dist

        dist = 0

        dist = dist + (p1(1) - p2(1))**2
        dist = dist + (p1(2) - p2(2))**2
        dist = dist + (p1(3) - p2(3))**2

        dist = sqrt(dist)

    end subroutine calculate_distance


END MODULE maths

