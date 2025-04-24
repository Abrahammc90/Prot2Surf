MODULE maths

    USE mod_pdb

    contains

    subroutine tr(vorig, x, y, z, vnew)

        IMPLICIT NONE
        real ( kind=8 ), dimension ( 3 ), intent ( inout ) :: vnew
        real( kind=8 ),dimension( 3 ), intent ( in ) :: vorig, x, y, z

        vnew(1) = dot_product(vorig, x)
        vnew(2) = dot_product(vorig, y)
        vnew(3) = dot_product(vorig, z)

    end subroutine tr

    subroutine norm(v)
        IMPLICIT NONE
        real ( kind = 8 ), dimension ( 3 ) :: v
        real ( kind = 8 ) :: magnitude_v
        magnitude_v = sqrt ( dot_product ( v, v ) )
        v = v / magnitude_v
    end subroutine norm

    subroutine cross(v1,v2,v3)
        IMPLICIT NONE
        
        real ( kind=8 ), dimension ( 3 ), intent ( in ) :: v1, v2
        real ( kind=8 ), dimension ( 3 ), intent ( out ) :: v3  
        v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
        v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
        v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
    end subroutine cross

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

        do i = 1, nb_atoms
            new_coords(i, 1) = coords(i, 1) - xc2(1)
            new_coords(i, 2) = coords(i, 2) - xc2(2)
            new_coords(i, 3) = coords(i, 3) - xc2(3)
        end do

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

        !ozz = (/1, 0, 0/)
        !zoz = (/0, 1, 0/)
        !zzo = (/0, 0, 1/)
!
!
!
        !call tr(ozz, rot_vx_norm, rot_vy_norm, rot_vz_norm, xee)
        !call tr(zoz, rot_vx_norm, rot_vy_norm, rot_vz_norm, yee)
        !call tr(zzo, rot_vx_norm, rot_vy_norm, rot_vz_norm, zee)

        !write(*,*) xee
        !write(*,*) yee
        !write(*,*) zee
        !stop

        trans_coords = new_coords
        do i = 1, nb_atoms

            new_coords(i, 1) = dot_product(trans_coords(i,:), xee)
            new_coords(i, 2) = dot_product(trans_coords(i,:), yee)
            new_coords(i, 3) = dot_product(trans_coords(i,:), zee)

            !call tr(trans_coords(i,:), xee, yee, zee, new_coords(i,:))
            !write(*,*) 'new coords 1', new_coords
            !call tr(yee, rot_vx_norm, rot_vy_norm, rot_vz_norm, new_coords(i,:))
            !write(*,*) 'new coords 2', new_coords
            !call tr(zee, rot_vx_norm, rot_vy_norm, rot_vz_norm, new_coords(i,:))
            !write(*,*) 'new coords 3', new_coords
            !STOP
        end do

        do i = 1, nb_atoms
            new_coords(i, 1) = new_coords(i, 1) + xc1(1) + trans(1)
            new_coords(i, 2) = new_coords(i, 2) + xc1(2) + trans(2)
            new_coords(i, 3) = new_coords(i, 3) + xc1(3) + trans(3)
        end do

    end subroutine update_complex

    subroutine vectors_angle(vector1, vector2, theta_degrees)

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

    end subroutine vectors_angle

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


    end subroutine

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

    end subroutine

    subroutine calculate_distance(p1, p2, dist)

        real (kind=8), dimension(3), intent(in) :: p1, p2
        real (kind=8), intent(out) :: dist

        dist = 0

        dist = dist + (p1(1) - p2(1))**2
        dist = dist + (p1(2) - p2(2))**2
        dist = dist + (p1(3) - p2(3))**2

        dist = sqrt(dist)

    end subroutine

END MODULE maths

