Module Forces_Module
    Implicit none
    
    Private
    Public :: find_force_LJ
    
    Contains
    
    Subroutine find_force_LJ(r, N, L, cutoff, F, pot)
        ! """"
        ! Calculates the forces applied to each particle of the system
        ! INPUTS: r, Nm L, cutoff, F
        ! OUTPUT: pot, F
        ! """"
        Implicit none
        real(8), dimension(N, 3), intent(in) :: r
        real(8), intent(in) :: L, cutoff
        real(8) :: d, f_ij
        real(8), dimension(3) :: d_r
        integer :: i, j, k
        integer, intent(in) :: N
        real(8), dimension(N, 3), intent(out) :: F
        real(8), intent(out) :: pot
    
        pot = 0.d0
    
        F = 0.d0
    
        do i = 1, N
            do j = i + 1, N
                d_r(:) = r(i, :) - r(j, :)
    
                do while (any(d_r(:) .gt. L/2.) .or. (any(d_r(:) .lt. (-L/2.))))
                    call pbc(d_r, L, size(d_r))
                end do
    
                d = (d_r(1)**2 + d_r(2)**2 + d_r(3)**2)**(1.d0/2.d0)
                if (d .le. cutoff) then
                    f_ij = 48.d0/d**14 - 24.d0/d**8
                    F(i, :) = F(i, :) + f_ij*d_r(:)
                    F(j, :) = F(j, :) - f_ij*d_r(:)
    
                    if (isnan(F(i, 1))) then
                        print *, i, j
                        stop
                    end if
    
                    pot = pot + 4.d0*(1.d0/d**12 - 1.d0/d**6) - 4.d0*(1/cutoff**12 - 1.d0/cutoff**6)
    
                end if
    
            end do
        end do
    
    End Subroutine find_force_LJ
    
End Module Forces_Module
