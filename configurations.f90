module configurations
    use globals
    contains

    subroutine generate_2h_configurations()
        integer :: i, j, counter

        counter = 1
        do i = h_f_min, h_f_max-1
            do j = i+1, h_f_max
                config_2h(counter, 1) = i
                config_2h(counter, 2) = j
                counter = counter + 1
            end do
        end do

    end subroutine generate_2h_configurations

    subroutine generate_1h1p_configurations()
        integer :: i, j, counter

        counter = 1
        do i = 1, n_hs_in
            do j = 1, n_ps_in
                config_1h1p(counter, 1) = holes_arr(i)
                config_1h1p(counter, 2) = particles_arr(j)
                counter = counter + 1
            end do
        end do

    end subroutine generate_1h1p_configurations
    
    subroutine generate_2h_config_str()
        integer :: i, j, counter

        counter = 1
        do i = 1, n_config
            write(config_2h_str(counter),'(A1,I3,A1,I3,A4)')"|", config_2h(i, 1), " ", config_2h(i, 2), ", I>"
            counter = counter + 1
        end do

        do i = h_f_min, h_f_max
            write(config_2h_str(counter),'(A1,I3,A1,I3,A4)')"|", i," ", i, ", I>"
            counter = counter + 1
        end do

        do i = 1, n_config
            write(config_2h_str(counter),'(A1,I3,A1,I3,A6)')"|", config_2h(i, 1), " ", config_2h(i, 2), ", III>"
            counter = counter + 1
        end do

    end subroutine generate_2h_config_str

end module configurations
