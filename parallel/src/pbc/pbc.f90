module pbc_module

   implicit none

contains

   Subroutine pbc_mic(vector, L, D)
      ! """"
      ! Applies periodic boundary conditions to a N-dimensional system, minimum image convention
      ! INPUTS: vector, L, D (dimension of the system)
      ! OUTPUT: vector
      ! """"
      Implicit none
      integer :: i
      integer, intent(in) :: D
      real(8), dimension(D), intent(inout) :: vector
      real(8), intent(in) :: L

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

   Subroutine pbc(vector, L, D)
      ! """"
      ! Applies periodic boundary conditions to a N-dimensional system, if a particle is outside the box, it is placed inside the box.
      ! (Use only if the box is not centered at the origin.)
      ! INPUTS: vector, L, D (dimension of the system)
      ! OUTPUT: vector
      ! """"
      Implicit none
      integer :: i
      integer, intent(in) :: D
      real(8), dimension(D), intent(inout) :: vector
      real(8), intent(in) :: L

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
