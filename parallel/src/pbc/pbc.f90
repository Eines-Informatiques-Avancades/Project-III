!> Module containing subroutines to apply periodic boundary conditions to a N-dimensional system.
module pbc_module

   implicit none

contains

!> Subroutine to apply periodic boundary conditions to a N-dimensional system, minimum image convention
!> @param vector Vector to apply PBC
!> @param L Box size
!> @param D Dimension of the system

   Subroutine pbc_mic(vector, L, D)
      ! """"
      ! Applies periodic boundary conditions to a N-dimensional system, minimum image convention
      ! INPUTS: vector, L, D (dimension of the system)
      ! OUTPUT: vector
      ! """"
      Implicit none
      integer :: i !> Loop variable
      integer, intent(in) :: D !> Dimension of the system
      real(8), dimension(D), intent(inout) :: vector !> Vector to apply PBC
      real(8), intent(in) :: L !> Box size

      do i = 1, D
         if (vector(i) .gt. L/2.) then
            vector(i) = vector(i) - L

         else if (vector(i) .lt. (-L/2.)) then
            vector(i) = vector(i) + L

         end if

         if (abs(vector(i)) .gt. 1000) then
            print *, "YOU HAVE INESTABILITIES!"
            stop
         end if

      end do

      Return
   End Subroutine

   !########################################################################################################

!> Subroutine to apply periodic boundary conditions to a N-dimensional system
!> @param vector Vector to apply PBC
!> @param L Box size
!> @param D Dimension of the system

   Subroutine pbc(vector, L, D)
      ! """"
      ! Applies periodic boundary conditions to a N-dimensional system, if a particle is outside the box, it is placed inside the box.
      ! (Use only if the box is not centered at the origin.)
      ! INPUTS: vector, L, D (dimension of the system)
      ! OUTPUT: vector
      ! """"
      Implicit none
      integer :: i !> Loop variable
      integer, intent(in) :: D !> Dimension of the system
      real(8), dimension(D), intent(inout) :: vector !> Vector to apply PBC
      real(8), intent(in) :: L !> Box size

      do i = 1, D
         if (vector(i) .gt. L) then
            vector(i) = vector(i) - L

         else if (vector(i) .lt. (-L)) then
            vector(i) = vector(i) + L

         end if

         if (abs(vector(i)) .gt. 1000) then
            print *, "YOU HAVE INESTABILITIES!"
            stop
         end if

      end do

      Return
   End Subroutine

   !########################################################################################################

end module pbc_module
