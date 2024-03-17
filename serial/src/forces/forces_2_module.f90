module forces_2_module
    implicit none
    
    private
    
    public :: get_forces
    public :: get_pot_ener
    
contains

    subroutine get_forces (forces,distances,positions,nbr_atomes,sigma,eps)
        implicit none
        
        real*8, allocatable, dimension (:,:),intent(inout) :: forces
        real*8, allocatable, dimension (:,:),intent(in)    :: distances
        real*8, allocatable, dimension (:,:),intent(in)    :: positions
        
        real*8,intent(in)  :: sigma,eps
        real*8             :: sigma_12,sigma_6
        real*8             :: dist2
        real*8             :: dr4,dr8,dr14
        real*8             :: part_6,part_12
        integer,intent(in) :: nbr_atomes
        integer            :: i,j,k
        
        do i=1,nbr_atomes
            do j=1,3
                forces(i,j)=0.0d0
            enddo
        enddo
        
        sigma_6=sigma**6
        sigma_12=sigma**12
        
        do i=1,nbr_atomes-1
            do k=i+1,nbr_atomes
                do j=1,3
                    dist2=(positions(i,j)-positions(k,j))
                    dr4=distances(i,k)*distances(i,k)
                    dr8=dr4*dr4
                    dr14=dr8*dr4*distances(i,k)
                    part_6  = sigma_6*6.0d0*(-(dist2)/dr8)
                    part_12 = sigma_12*6.d0*(-(2.0d0*dist2)/dr14) 
                    forces(i,j)=forces(i,j)-4.0d0*eps*(part_12-part_6)
                    forces(k,j)=forces(k,j)+4.0d0*eps*(part_12-part_6)
                enddo
            enddo
        enddo
        
    end subroutine get_forces
    
    subroutine get_pot_ener (distances,nbr_atomes,sigma,eps,energy)
        implicit none
        
        real*8, allocatable, dimension (:,:),intent(in)    :: distances
        
        real*8,intent(in)   :: sigma,eps
        real*8              :: sigma_6,sigma_12
        real*8              :: part_6,part_12
        real*8              :: dr6,dr12
        real*8,intent(inout):: energy
        integer,intent(in)  :: nbr_atomes
        integer             :: i,k
        
        energy=0.0d0
        
        sigma_6=sigma**6
        sigma_12=sigma**12
        
        do i=1,nbr_atomes-1
                do k=i+1,nbr_atomes
                    dr6=distances(i,k)*distances(i,k)*distances(i,k)
                    dr12=dr6*dr6
                    part_6  = sigma_6/dr6
                    part_12 = sigma_12/dr12
                    energy=energy+4*eps*(part_12-part_6)
                enddo
        enddo
        
    end subroutine get_pot_ener

end module forces_2_module
