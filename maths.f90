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

    subroutine vector_perpendicular_to_plane(p1, p2, p3, n)
        IMPLICIT NONE
        real ( kind=8 ), dimension ( 3 ), intent ( out ) :: n
        real( kind=8 ),dimension( 3 ), intent ( in ) :: p1, p2, p3
        real( kind=8 ),dimension( 3 ) :: v1, v2
        real( kind=8 ) :: magnitude_n

        v1 = p3 - p1
        v2 = p2 - p1

        call cross(v1, v2, n)
        !Normalizing n
        magnitude_n = sqrt ( dot_product ( n, n ) )
        n = n / magnitude_n
    end subroutine vector_perpendicular_to_plane

    subroutine vectors_angle(vector1, vector2, theta_degrees)

        IMPLICIT NONE
        real (kind=8), intent(out) :: theta_degrees
        real( kind=8 ),dimension( 3 ), intent ( in ) :: vector1, vector2
        real (kind=8) :: magnitude_v1, magnitude_v2, cos_theta
        real (kind=8) :: theta_radians
        real(kind=8), dimension(3) :: v1, v2
        real (kind=8) :: test_v1, test_v2

        !v1 = (/1, 0, 0/)
        !v2 = (/0, 1, 0/)

        !test_v1 = sqrt ( dot_product ( v1,v1 ))
        !test_v2 = sqrt ( dot_product ( v2,v2 ))

        !cos_theta = dot_product(v1, v2) / ( test_v1 * test_v2 )
        !theta_radians = acos(cos_theta)
        !theta_degrees = theta_radians * (180.0 / 3.141592653589793)

        !write(*,*) theta_degrees
        !STOP 1

        magnitude_v1 = sqrt ( dot_product ( vector1, vector1 ) )
        magnitude_v2 = sqrt ( dot_product ( vector2, vector2 ) )
        cos_theta = dot_product(vector1, vector2) / ( magnitude_v1 * magnitude_v2 )

        !if (cos_theta > 1.0) cos_theta = 1.0
        !if (cos_theta < -1.0) cos_theta = -1.0
        theta_radians = acos(cos_theta)
        theta_degrees = theta_radians * (180.0 / 3.141592653589793)

    end subroutine vectors_angle

    subroutine rmsd_test ( rmsd2, dist_max, r1,y11,y12,r2,y21,y22 )
        IMPLICIT NONE
            real ( kind=8 ), intent ( out ) :: rmsd2
            real ( kind=8 ), intent ( in ) :: dist_max
            real ( kind=8 ), dimension ( 3 ) :: r1, y11, y12, y13
            real ( kind=8 ), dimension ( 3 ) :: r2, y21, y22, y23
            real ( kind=8 ), dimension ( 3 ) :: yy11, yy12, yy13
            real ( kind=8 ), dimension ( 3 ) :: yy21, yy22, yy23
            real ( kind=8 ), dimension ( 3 ) :: yy1m1, yy1m2, yy1m3
            real ( kind=8 ), dimension ( 3 ) :: yy2m1, yy2m2, yy2m3
    
            rmsd2 = 0.d0
            call cross( y13, y11, y12 )
            call cross( y23, y21, y22 )
    
            ! write(*,*) 'dismx ',dismx
            ! write(*,*) 'r1', (r1(m),m=1,3)
            ! write(*,*) 'y11', (y11(m),m=1,3)
            ! write(*,*) 'y12', (y12(m),m=1,3)
            ! write(*,*) 'y13', (y13(m),m=1,3)
            ! write(*,*) 'r2', (r2(m),m=1,3)
            ! write(*,*) 'y21', (y21(m),m=1,3)
            ! write(*,*) 'y13 ', (y13(m),m=1,3 )
            ! write(*,*) 'y23 ', (y23(m),m=1,3)
    
            yy11 = y11 * dist_max + r1
            yy12 = y12 * dist_max + r1
            yy13 = y13 * dist_max + r1
            yy21 = y21 * dist_max + r2
            yy22 = y22 * dist_max + r2
            yy23 = y23 * dist_max + r2
    
            yy1m1 = -y11 * dist_max + r1
            yy1m2 = -y12 * dist_max + r1
            yy1m3 = -y13 * dist_max + r1
            yy2m1 = -y21 * dist_max + r2
            yy2m2 = -y22 * dist_max + r2
            yy2m3 = -y23 * dist_max + r2
    
            rmsd2 = dot_product ( r1 - r2, r1 - r2 ) + &
                   2. * dot_product ( yy11 - yy21, yy11 - yy21 ) + &
                   2. * dot_product ( yy12 - yy22, yy12 - yy22 ) + &
                   2. * dot_product ( yy13 - yy23, yy13 - yy23 ) + &
                   dot_product ( yy1m1 - yy2m1, yy1m1 - yy2m1 ) + &
                   dot_product ( yy1m2 - yy2m2, yy1m2 - yy2m2 ) + &
                   dot_product ( yy1m3 - yy2m3, yy1m3 - yy2m3 )
    
            rmsd2 = rmsd2 / 7.
    end subroutine rmsd_test

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

