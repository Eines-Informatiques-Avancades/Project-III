module file_writers

    contains


    subroutine read_param()

        

    end subroutine


    subroutine write_trajectory


    end subroutine


    subroutine write_thermodynamics(step, dt, pressure, Temp, pot, K_energy,p)
        
        write (96, *) step*dt, Pressure
        write (77, *) step*dt, Temp
        write (44, *) step*dt, pot, K_energy, pot + K_energy, p
    
    end subroutine 



end module 

